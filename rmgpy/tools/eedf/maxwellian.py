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

"""Offline Maxwellian comparator rows from a byte-pinned native source table.

This module must never enter the runtime provider import graph. Historical
source rows are decoded here without invoking the runtime EEDFTable loader;
that loader's implementation fingerprint and production gates remain binding.
"""
import copy
import hashlib
import io
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.constants import elementary_charge, electron_mass

from rmgpy.tools.eedf.artifact import _storage_items, _validate_storage, NUMERIC, MAPPINGS
from rmgpy.tools.eedf.moments import _distribution_moments, _transport_moments
from rmgpy.tools.eedf.schema import (
    EEDFError, FingerprintMismatch, SpecError, _artifact_provider_kind, _content_hash,
)


COMPARATOR_BRANCH = 'maxwellian'


def _source_snapshot(path, *, artifact_sha256, branch_id):
    """Decode native tensors from the very bytes whose SHA was checked.

    The historical interpolant is provenance, not an executable dependency of
    conversion. Its fingerprint still binds all source row inputs to the HDF5.
    No historical qualification verdict is promoted to comparator validation.
    """
    path = Path(path)
    if path.is_file():
        path = path.parent
    manifest = json.loads((path / 'manifest.json').read_bytes())
    contents = (path / 'table.h5').read_bytes()
    digest = hashlib.sha256(contents).hexdigest()
    if digest != artifact_sha256 or digest != manifest.get('artifact_sha256'):
        raise FingerprintMismatch('source HDF5 artifact content SHA-256')
    if manifest.get('schema_version') != 1 or _artifact_provider_kind(manifest) != 'loki-table':
        raise FingerprintMismatch('source must be a schema-1 loki-table')
    if not isinstance(branch_id, str) or branch_id not in manifest['branches']:
        raise SpecError('select an existing source branch explicitly')
    if _content_hash(manifest['row_inputs']) != manifest['fingerprint']:
        raise FingerprintMismatch('source row-input fingerprint')
    axes = {'u': manifest['axes']['u'], **{
        name: values for name, values in manifest['axes'].items() if name != 'u'}}
    if axes != manifest['row_inputs']['axes']:
        raise FingerprintMismatch('source axes differ from row inputs')
    for name, values in axes.items():
        axis = np.asarray(values, dtype=float)
        if axis.ndim != 1 or len(axis) < 2 or not np.all(np.isfinite(axis)) or np.any(np.diff(axis) <= 0):
            raise SpecError('invalid source axis ' + name)
    shape = tuple(len(values) for values in axes.values())
    rows = np.empty(shape, dtype=object)
    with h5py.File(io.BytesIO(contents), 'r') as h5:
        _validate_storage(h5)
        if h5.attrs.get('schema_version') != manifest['schema_version']:
            raise FingerprintMismatch('source schema version')
        stored_kind = h5.attrs.get('provider_kind')
        if stored_kind != manifest.get('provider_kind'):
            raise FingerprintMismatch('source provider kind')
        bindings = {
            'fingerprint': manifest['fingerprint'],
            'channel_map_content_sha256': _content_hash(manifest['channel_map']),
            'branch_certification_sha256': _content_hash(manifest['branch_certification']),
            'policy_sha256': _content_hash({key: manifest[key] for key in ('tolerances', 'floors', 'envelopes')}),
        }
        for name, expected in bindings.items():
            if h5.attrs.get(name) != expected:
                raise FingerprintMismatch('source ' + name)
        stored_axes = dict(_storage_items(h5['axes']))
        if set(stored_axes) != set(axes):
            raise FingerprintMismatch('source axis names')
        for name, values in axes.items():
            if not np.array_equal(stored_axes[name][...], values):
                raise FingerprintMismatch('source axis ' + name)
        for name in ('energy_eV', 'energy_edges_eV'):
            if not np.array_equal(h5[name][...], manifest[name]):
                raise FingerprintMismatch('source ' + name)
        branches = dict(_storage_items(h5['branches']))
        if set(branches) != set(manifest['branches']):
            raise FingerprintMismatch('source branches')
        group = branches[branch_id]
        arrays = {name: group[name][...] for name in NUMERIC if name in group}
        mappings = {name: {key: value[...] for key, value in _storage_items(group[name])}
                    for name in MAPPINGS}
        for name, values in list(arrays.items()) + [
                (name + '/' + key, value) for name, mapping in mappings.items() for key, value in mapping.items()]:
            if values.shape[:len(shape)] != shape or not np.all(np.isfinite(values)):
                raise FingerprintMismatch('source tensor ' + name)
        for index in np.ndindex(shape):
            row = {name: value[index].copy() for name, value in arrays.items()}
            row.update({name: {key: float(value[index]) for key, value in mapping.items()}
                        for name, mapping in mappings.items()})
            point = {name: float(axes[name][i]) for name, i in zip(axes, index)}
            if row['u'] != point['u'] or row['EN_Td'] != np.exp(point['u']):
                raise FingerprintMismatch('source reference coordinate')
            for name in ('energy_eV', 'energy_edges_eV'):
                if not np.array_equal(row[name], manifest[name]):
                    raise FingerprintMismatch('source row ' + name)
            row['composition'] = {name: value for name, value in point.items() if name != 'u'}
            row['branch_id'] = branch_id
            rows[index] = row
    manifest['axes'] = axes
    return manifest, rows


def _maxwellian_distribution(energy, edges, temperature_eV, mean_rtol):
    """Normalize on the unchanged source grid, refusing unresolved means."""
    energy, edges = np.asarray(energy, dtype=float), np.asarray(edges, dtype=float)
    if (energy.ndim != 1 or len(energy) < 2 or edges.shape != (len(energy) + 1,)
            or not np.all(np.isfinite(energy)) or not np.all(np.isfinite(edges))
            or edges[0] != 0 or np.any(np.diff(edges) <= 0)
            or not np.allclose(np.diff(edges), edges[1], rtol=1e-10, atol=0)
            or not np.allclose(energy, (edges[:-1] + edges[1:]) / 2, rtol=1e-12, atol=0)):
        raise SpecError('Maxwellian requires a uniform source cell-centre energy grid')
    if (not np.isfinite(temperature_eV) or temperature_eV <= 0
            or not np.isfinite(mean_rtol) or mean_rtol <= 0):
        raise SpecError('Maxwellian temperature and mean tolerance must be positive')
    raw = 2 / np.sqrt(np.pi) * temperature_eV**-1.5 * np.exp(-energy / temperature_eV)
    norm = float(np.sum(raw * np.sqrt(energy) * np.diff(edges)))
    if not np.isfinite(norm) or norm <= 0:
        raise EEDFError('unresolved Maxwellian normalization')
    f0 = raw / norm
    mean = float(np.sum(f0 * energy**1.5 * np.diff(edges)))
    relative_error = mean / (1.5 * temperature_eV) - 1
    if abs(relative_error) > mean_rtol:
        raise EEDFError('unresolved Maxwellian mean: relative error ' + str(relative_error))
    return f0, {'temperature_eV': float(temperature_eV), 'raw_normalization': norm,
                'mean_energy_eV': mean, 'relative_mean_error': relative_error}


def _maxwellian_row(source, manifest):
    """Recompute collision and transport data, retaining only physical inputs."""
    source, manifest = copy.deepcopy((source, manifest))
    options = manifest['row_inputs']['solver_options']
    if (options.get('growthModelType') != 'temporal'
            or options.get('ionizationOperatorType') != 'usingSDCS'
            or options.get('includeEECollisions') is not False):
        raise SpecError('unsupported source collision/growth convention')
    channels = manifest['channel_map']
    energy, edges = source['energy_eV'], source['energy_edges_eV']
    temperature = 2 * source['swarm']['mean_energy_eV'] / 3
    f0, resolution = _maxwellian_distribution(
        energy, edges, temperature, manifest['tolerances']['H2']['rtol'])
    row = {name: copy.deepcopy(source[name]) for name in (
        'u', 'EN_Td', 'energy_eV', 'energy_edges_eV', 'composition',
        'gas_fractions', 'state_populations', 'target_fractions', 'product_fractions', 'rate_floors')}
    row.update(f0=f0, branch_id=COMPARATOR_BRANCH, converged=False, iteration_count=-1,
               setup='prescribed Maxwellian; no Boltzmann solve',
               gas_temperature_K=source['composition'].get('Tg_K', manifest['row_inputs']['Tg_K']))
    for name in ('target_fractions', 'product_fractions', 'rate_floors'):
        values = np.asarray(row[name])
        if values.shape != (len(channels),) or not np.all(np.isfinite(values)) or np.any(values < 0):
            raise SpecError('invalid source ' + name)
    if np.any(row['rate_floors'] <= 0):
        raise SpecError('nonpositive source rate floor')
    momentum = np.zeros_like(edges)
    for j, channel in enumerate(channels):
        section = channel.get('cross_section', {})
        x = np.asarray(section.get('energy_eV', []), dtype=float)
        sigma = np.asarray(section.get('sigma_m2', []), dtype=float)
        if (x.ndim != 1 or len(x) < 2 or sigma.shape != x.shape
                or not np.all(np.isfinite(x)) or not np.all(np.isfinite(sigma))
                or np.any(np.diff(x) <= 0) or np.any(sigma < 0) or x[0] < 0
                or x[-1] < edges[-1]):
            raise SpecError('unresolved source cross section: ' + channel['description'])
        if channel['kind'] == 'elastic' and (
                not np.isfinite(channel.get('mass_ratio', np.nan)) or channel['mass_ratio'] <= 0):
            raise SpecError('unresolved elastic mass ratio')
        forward = np.interp(edges, x, sigma, left=0, right=0)
        if channel['kind'] != 'elastic':
            forward[edges <= channel['threshold_eV']] = 0
        momentum += row['target_fractions'][j] * forward
        if '<->' in channel['description']:
            ratio = channel.get('statistical_weight_ratio', np.nan)
            if not np.isfinite(ratio) or ratio <= 0:
                raise SpecError('unresolved reverse statistical weight ratio')
            shifted = edges[1:] + channel['threshold_eV']
            reverse = ratio * shifted / edges[1:] * np.interp(shifted, x, sigma, left=0, right=0)
            momentum[1:] += row['product_fractions'][j] * reverse
    moments = _distribution_moments(row, channels)
    row.update({name: moments[name] for name in ('k_ine', 'k_sup', 'channel_power')})
    growth = sum(row['target_fractions'][j] * moments['k_ine'][j] * (
        1 if channel['kind'] == 'ionization' else -1 if channel['kind'] == 'attachment' else 0)
        for j, channel in enumerate(channels))
    gamma = np.sqrt(2 * elementary_charge / electron_mass)
    if np.any(momentum[1:] + growth / gamma / np.sqrt(edges[1:]) <= 0):
        raise EEDFError('nonpositive Maxwellian transport collision frequency')
    transport = _transport_moments(energy, edges, f0, momentum, growth)
    if any(not np.isfinite(value) or value <= 0 for value in transport.values()):
        raise EEDFError('nonpositive Maxwellian transport moment')
    row['swarm'] = dict(transport, mean_energy_eV=moments['mean_energy_eV'],
                        maxwellian_temperature_eV=temperature,
                        maxwellian_raw_normalization=resolution['raw_normalization'],
                        maxwellian_relative_mean_error=resolution['relative_mean_error'])
    row['attachment_energy_eV'] = np.zeros(len(channels))
    for j, channel in enumerate(channels):
        if channel['kind'] == 'attachment':
            # Evaluate the unweighted energy moment even for an empty target.
            unit = dict(row, target_fractions=np.ones(len(channels)))
            power = _distribution_moments(unit, channels)['channel_power'][j]
            if row['k_ine'][j] > 0:
                row['attachment_energy_eV'][j] = power / row['k_ine'][j]
    row['below_floor'] = np.stack([row['k_ine'] < row['rate_floors'],
                                   row['k_sup'] < row['rate_floors']], axis=-1)
    row['power_groups'] = {'collision_loss': float(np.sum(row['channel_power']))}
    return row
