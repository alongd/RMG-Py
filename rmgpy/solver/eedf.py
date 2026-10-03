#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""Fingerprint-checked, branch-preserving tensor PCHIP EEDF tables.

The table returns one consistent row. It neither regenerates artifacts nor
selects a reactor operating point. Coordinates are u = ln(E/N / 1 Td).
"""
from dataclasses import dataclass
import itertools
import io
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.interpolate import PchipInterpolator

from rmgpy.tools.eedf.artifact import storage_items, validate_storage
from rmgpy.tools.eedf.schema import (
    EEDFError, FingerprintMismatch, check_model, content_hash,
    interpolant_identity,
)


class OutOfDomain(EEDFError):
    """An axis value is outside the qualified domain."""


class EnvelopeBreach(EEDFError):
    """A non-axis input left its screened envelope."""


class AmbiguousBranch(EEDFError):
    """An unambiguous branch was not explicitly selected."""


class IllConditionedCoordinate(EEDFError):
    """Mean energy folds or its u derivative is unusable."""


@dataclass(frozen=True)
class EEDFRow:
    """All chemistry, transport and energy quantities from one interpolation.

    Rates: m3/s; mean/attachment energy: eV; reduced transport: LoKI SI
    units; channel/group power: eV m3/s. iteration_count=-1 means the pinned
    solver does not expose the count. The manifest records this limitation.
    """
    u: float
    EN_Td: float
    branch_id: str
    composition: dict
    swarm: dict
    gas_fractions: dict
    state_populations: dict
    k_ine: np.ndarray
    k_sup: np.ndarray
    channel_power: np.ndarray
    power_groups: dict
    attachment_energy_eV: np.ndarray
    target_fractions: np.ndarray
    product_fractions: np.ndarray
    rate_floors: np.ndarray
    below_floor: np.ndarray
    f0: np.ndarray
    energy_eV: np.ndarray
    energy_edges_eV: np.ndarray
    converged: bool
    iteration_count: int
    denergy_du: float
    fingerprint: str

    def as_dict(self):
        """Mapping for validation against an independent direct solve."""
        return dict(vars(self))


class EEDFTable:
    """An immutable set of qualified branch tensors.

    ``load(path, current_model, artifact_sha256=expected_content_hash)`` requires the complete current row-input
    identity, usually constructed with ``model_inputs(spec)``. Database and
    engine SHAs are provenance only and must not appear in that identity.
    """

    @classmethod
    def load(cls, path, current_model, *, artifact_sha256, require_accepted=True):
        """Check provenance, physical inputs, implementation, and file bytes.

        ``require_accepted=False`` is reserved for offline qualification;
        ordinary loads refuse unsuccessful independent held-out verdicts."""
        path = Path(path)
        if path.is_file():
            path = path.parent
        try:
            with (path / 'manifest.json').open() as stream:
                manifest = json.load(stream)
        except FileNotFoundError as error:
            raise FingerprintMismatch('incomplete artifact: missing manifest') from error
        results = manifest.get('held_out_verdicts', [])
        accepted = bool(results) and all(result.get('passed') is True and bool(result.get('checks'))
            and all(check.get('passed') is True for check in result['checks']) for result in results)
        manifest['accepted'] = accepted
        if require_accepted and not accepted:
            raise FingerprintMismatch('held-out qualification')
        if require_accepted and any(channel.get('classification') not in ('A', 'B')
                                    for channel in manifest.get('channel_map', [])):
            raise FingerprintMismatch('channel classification is excluded from production')
        required = {'schema_version', 'row_inputs', 'fingerprint', 'axes', 'envelopes',
                    'tolerances', 'floors', 'channel_map', 'energy_eV',
                    'energy_edges_eV', 'branches', 'artifact_sha256', 'branch_certification'}
        if not required.issubset(manifest) or manifest['schema_version'] != 1:
            raise FingerprintMismatch('manifest schema')
        if 'state_properties' in manifest['row_inputs']:
            from rmgpy.tools.eedf.channels import validate_map
            validate_map(manifest['channel_map'])
        check_model(manifest, current_model, manifest['tolerances']['fingerprint_rtol'])
        if manifest['row_inputs']['interpolant'] != interpolant_identity():
            raise FingerprintMismatch('interpolant implementation')
        if manifest['axes'] != manifest['row_inputs']['axes']:
            raise FingerprintMismatch('axes')
        energy = np.asarray(manifest['energy_eV'], dtype=float)
        edges = np.asarray(manifest['energy_edges_eV'], dtype=float)
        if (energy.ndim != 1 or edges.ndim != 1 or len(edges) != len(energy) + 1
                or not np.all(np.isfinite(energy)) or not np.all(np.isfinite(edges))
                or np.any(energy < 0) or np.any(np.diff(energy) <= 0)
                or np.any(edges < 0) or np.any(np.diff(edges) <= 0)):
            raise FingerprintMismatch('invalid energy grid')
        h5_path = path / 'table.h5'
        with h5_path.open('rb') as stream:
            contents = stream.read()
        digest = hashlib.sha256(contents).hexdigest()
        if not artifact_sha256 or digest != artifact_sha256 or digest != manifest['artifact_sha256']:
            raise FingerprintMismatch('HDF5 artifact content SHA-256')
        table = cls()
        table.manifest, table.path = manifest, path
        table.axis_names = ['u'] + [name for name in manifest['axes'] if name != 'u']
        table.axes = [np.asarray(manifest['axes'][name], dtype=float) for name in table.axis_names]
        table._data, table._layout, table._symbolic = {}, {}, {}
        table.branch_certification = dict(manifest['branch_certification'])
        if set(table.branch_certification) != set(manifest['branches']):
            raise FingerprintMismatch('branch certification keys')
        with h5py.File(io.BytesIO(contents), 'r') as h5:
            validate_storage(h5)
            if h5.attrs['branch_certification_sha256'] != content_hash(table.branch_certification):
                raise FingerprintMismatch('branch certification')
            if h5.attrs['qualification_sha256'] != content_hash({key: manifest.get(key) for key in ('accepted', 'held_out_verdicts', 'screen_results', 'branch_detection')}):
                raise FingerprintMismatch('held-out qualification content')
            if h5.attrs['channel_map_content_sha256'] != content_hash(manifest['channel_map']):
                raise FingerprintMismatch('channel map content')
            for name in ('energy_eV', 'energy_edges_eV'):
                if not np.array_equal(h5[name][...], manifest[name]):
                    raise FingerprintMismatch(name)
            if h5.attrs['fingerprint'] != manifest['fingerprint']:
                raise FingerprintMismatch('HDF5 fingerprint')
            if h5.attrs['policy_sha256'] != content_hash({key: manifest[key] for key in ('tolerances', 'floors', 'envelopes')}):
                raise FingerprintMismatch('manifest policy')
            stored_axes = dict(storage_items(h5['axes']))
            stored_branches = dict(storage_items(h5['branches']))
            for name, axis in zip(table.axis_names, table.axes):
                if not np.array_equal(stored_axes[name][...], axis):
                    raise FingerprintMismatch('HDF5 axis ' + name)
            if set(stored_branches) != set(manifest['branches']):
                raise FingerprintMismatch('branches')
            for name in ('energy_eV', 'energy_edges_eV'):
                if name in h5['held_out']:
                    dataset = h5['held_out'][name]
                    if not np.array_equal(dataset[...], np.broadcast_to(manifest[name], dataset.shape)):
                        raise FingerprintMismatch('held-out energy grid: ' + name)
            for branch_id in manifest['branches']:
                table._read_branch(branch_id, stored_branches[branch_id])
            table._check_population_rows(h5)
        return table

    def _check_population_rows(self, h5):
        """Validate generated training and held-out densities against row inputs."""
        inputs = self.manifest['row_inputs']
        if 'state_properties' not in inputs:
            return  # Synthetic interpolation fixtures have no solver setup.
        from rmgpy.tools.eedf.integrity import resolved_properties, channel_fractions
        from rmgpy.tools.eedf.schema import SpecError
        population_spec = dict(inputs, input_files={}, solver_options={},
                               gas_properties={'fraction': inputs['gas_properties']['fraction']})
        channels = self.manifest['channel_map']
        for group_name, group in [('held_out', h5['held_out'])] + list(storage_items(h5['branches'])):
            shape = group['target_fractions'].shape[:-1]
            held_coordinates = json.loads(group.attrs['composition']) if group_name == 'held_out' else None
            for index in np.ndindex(shape):
                coordinates = (held_coordinates[index[0]] if held_coordinates is not None else
                               {name: float(axis[i]) for name, axis, i in zip(self.axis_names, self.axes, index)
                                if name != 'u'})
                state = {name: item['reference'] for name, item in inputs['envelopes'].items()}
                state.update(coordinates)
                try:
                    _, _, gas, populations = resolved_properties(population_spec, coordinates)
                    target, product = channel_fractions(channels, gas, populations, state)
                except (SpecError, KeyError) as error:
                    raise FingerprintMismatch('row populations cannot be resolved') from error
                for name, expected in (('target_fractions', target), ('product_fractions', product)):
                    actual = group[name][index]
                    if actual.shape != np.shape(expected) or actual.tobytes() != np.asarray(expected, dtype=actual.dtype).tobytes():
                        raise FingerprintMismatch('row ' + name + ' differs from LoKI populations')
                for name, expected in (('gas_fractions', gas), ('state_populations', populations)):
                    stored = dict(storage_items(group[name]))
                    if set(stored) != set(expected):
                        raise FingerprintMismatch('row population keys')
                    if any(float(stored[key][index]).hex() != value.hex() for key, value in expected.items()):
                        raise FingerprintMismatch('row ' + name + ' differs from setup populations')

    def _read_branch(self, branch_id, group):
        shape = tuple(len(axis) for axis in self.axes)
        rank = len(shape)
        arrays, layout, symbolic = [], {}, {}
        offset = 0

        def append(name, values):
            nonlocal offset
            values = np.asarray(values)
            if values.shape[:rank] != shape or not np.all(np.isfinite(values)):
                raise FingerprintMismatch('branch tensor ' + name)
            tail = values.shape[rank:]
            count = int(np.prod(tail)) if tail else 1
            layout[name] = (offset, offset + count, tail)
            arrays.append(values.reshape(shape + (count,)))
            offset += count

        for name, value in group.items():
            if isinstance(value, h5py.Group):
                for key, dataset in storage_items(value):
                    append(name + '/' + key, dataset[...])
                symbolic[name] = {}
                for key in set(value.attrs) - {'identifier_encoding'}:
                    raise FingerprintMismatch('symbolic populations need resolved output: ' + key)
            elif name in ('energy_eV', 'energy_edges_eV'):
                expected = np.broadcast_to(self.manifest[name], value.shape)
                if not np.array_equal(value[...], expected):
                    raise FingerprintMismatch('row energy grid: ' + name)
            elif name not in ('u', 'EN_Td', 'iteration_count'):
                append(name, value[...])
        for direction in ('k_ine', 'k_sup'):
            start, stop, _ = layout[direction]
            original = np.concatenate(arrays, axis=-1)[..., start:stop]
            if np.any(original < 0):
                raise FingerprintMismatch('negative rate')
            append(direction + '_positive', (original > 0).astype(float))
        packed = np.concatenate(arrays, axis=-1)
        floor_start, floor_stop, _ = layout['rate_floors']
        floors = packed[..., floor_start:floor_stop]
        if np.any(floors <= 0):
            raise FingerprintMismatch('nonpositive rate floor')
        for direction in ('k_ine', 'k_sup'):
            start, stop, _ = layout[direction]
            packed[..., start:stop] = np.log1p(packed[..., start:stop] / floors)
        self._data[branch_id], self._layout[branch_id], self._symbolic[branch_id] = packed, layout, symbolic

    def domain_check(self, u, c):
        """Refuse missing/unknown coordinates, extrapolation and envelope exits."""
        if not isinstance(c, dict):
            raise OutOfDomain('composition must be a named coordinate mapping')
        axis_names = set(self.axis_names[1:])
        envelopes = self.manifest['envelopes']
        unknown = set(c) - axis_names - set(envelopes)
        if unknown:
            raise EnvelopeBreach('unqualified coordinate ' + str(sorted(unknown)))
        if not axis_names.issubset(c):
            raise OutOfDomain('missing composition axis')
        point = [u] + [c[name] for name in self.axis_names[1:]]
        for name, axis, value in zip(self.axis_names, self.axes, point):
            if not np.isfinite(value) or not axis[0] <= value <= axis[-1]:
                raise OutOfDomain(name)
        for name, bounds in envelopes.items():
            value = c.get(name, bounds['reference'])
            if not np.isfinite(value) or not bounds['min'] <= value <= bounds['max']:
                raise EnvelopeBreach(name)
        return point

    def row(self, u, c, branch_id=None):
        """Tensor PCHIP for every quantity; logarithmic for positive rates.

        Composition dimensions are reduced first, in manifest axis order;
        u is reduced last. The returned derivative is exactly that of the
        returned mean energy, including the composition interpolation.
        """
        point = self.domain_check(u, c)
        if branch_id is None:
            if len(self._data) != 1:
                raise AmbiguousBranch('select one of ' + str(self.manifest['branches']))
            branch_id = next(iter(self._data))
        if branch_id not in self._data:
            raise AmbiguousBranch('unknown branch ' + str(branch_id))
        data = self._data[branch_id]
        for dimension in range(len(self.axes) - 1, 0, -1):
            data = PchipInterpolator(self.axes[dimension], data, axis=dimension, extrapolate=False)(point[dimension])
        layout = self._layout[branch_id]
        energy_start = layout['swarm/mean_energy_eV'][0]
        energy = data[:, energy_start]
        deltas = np.diff(energy)
        if np.any(deltas == 0) or (np.any(deltas > 0) and np.any(deltas < 0)):
            raise IllConditionedCoordinate('folded/flat mean-energy branch')
        interpolation = PchipInterpolator(self.axes[0], data, axis=0, extrapolate=False)
        values = interpolation(u)
        derivative = float(interpolation.derivative()(u)[energy_start])
        conditioning = self.manifest['tolerances']['u_condition']
        mean_energy = float(values[energy_start])
        if (not np.isfinite(derivative) or mean_energy <= 0 or
                abs(derivative) <= conditioning['min_abs_denergy_du'] or
                abs(derivative) / mean_energy <= conditioning['min_scaled_denergy_du']):
            raise IllConditionedCoordinate('d mean_energy / du')

        def take(name):
            start, stop, tail = layout[name]
            value = values[start:stop].reshape(tail)
            return float(value) if not tail else value.copy()

        mappings = {name: dict(self._symbolic[branch_id].get(name, {}))
                    for name in ('swarm', 'power_groups', 'gas_fractions', 'state_populations')}
        for name in layout:
            if '/' in name:
                parent, key = name.split('/', 1)
                mappings[parent][key] = take(name)
        result = {name: take(name) for name in ('k_ine', 'k_sup', 'channel_power',
                  'attachment_energy_eV', 'target_fractions', 'product_fractions', 'rate_floors', 'f0')}
        for direction in ('k_ine', 'k_sup'):
            result[direction] = np.maximum(0., np.expm1(result[direction])) * result['rate_floors']
            result[direction][take(direction + '_positive') <= 0] = 0
        if np.any(result['f0'] < 0):
            raise EEDFError('negative interpolated EEDF')
        normalization = np.sum(result['f0'] * np.sqrt(self.manifest['energy_eV']) * np.diff(self.manifest['energy_edges_eV']))
        if not np.isfinite(normalization) or normalization <= 0:
            raise EEDFError('zero interpolated EEDF normalization')
        result['f0'] /= normalization
        if not np.all(np.isfinite(result['f0'])):
            raise EEDFError('nonfinite normalized EEDF')
        for name in ('mobility_N', 'diffusion_N', 'energy_mobility_N', 'energy_diffusion_N'):
            if name in mappings['swarm'] and mappings['swarm'][name] <= 0:
                raise EEDFError('nonpositive transport ' + name)
        # A flag is conservative over the cell stencil and both rate directions.
        indices = []
        for axis, value in zip(self.axes, point):
            hi = int(np.searchsorted(axis, value, side='right'))
            hi = min(max(hi, 1), len(axis) - 1)
            indices.append([hi - 1, hi])
        start, stop, tail = layout['below_floor']
        flags = np.any([self._data[branch_id][index][start:stop].reshape(tail).astype(bool)
                        for index in itertools.product(*indices)], axis=0)
        for j, direction in enumerate(('k_ine', 'k_sup')):
            flags[:, j] |= result[direction] < result['rate_floors']
        return EEDFRow(u=float(u), EN_Td=float(np.exp(u)), branch_id=branch_id,
                       composition={name: c[name] for name in self.axis_names[1:]},
                       **mappings, **result, below_floor=flags,
                       energy_eV=np.asarray(self.manifest['energy_eV']),
                       energy_edges_eV=np.asarray(self.manifest['energy_edges_eV']),
                       converged=bool(take('converged') == 1), iteration_count=-1,
                       denergy_du=derivative, fingerprint=self.manifest['fingerprint'])
