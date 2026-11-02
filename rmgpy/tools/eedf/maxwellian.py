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


def _held_out_source_row(point, manifest, rows):
    """Interpolate only the reference mean; resolve populations from source inputs."""
    from scipy.interpolate import PchipInterpolator
    from rmgpy.tools.eedf.integrity import _resolved_properties, _channel_fractions

    names = list(manifest['axes'])
    means = np.asarray([r['swarm']['mean_energy_eV'] for r in rows.flat]).reshape(rows.shape)
    for dimension in range(len(names) - 1, -1, -1):
        name = names[dimension]
        means = PchipInterpolator(manifest['axes'][name], means, axis=dimension,
                                  extrapolate=False)(point[name])
    row = copy.deepcopy(rows.flat[0])
    row.update(u=point['u'], EN_Td=float(np.exp(point['u'])),
               composition={name: point[name] for name in names if name != 'u'},
               swarm={'mean_energy_eV': float(means)})
    inputs = manifest['row_inputs']
    population_spec = dict(inputs, input_files={}, solver_options={},
                           gas_properties={'fraction': inputs['gas_properties']['fraction']})
    _, _, gas, populations = _resolved_properties(population_spec, row['composition'])
    state = {name: bounds['reference'] for name, bounds in manifest['envelopes'].items()}
    state.update(row['composition'])
    target, product = _channel_fractions(manifest['channel_map'], gas, populations, state)
    row.update(gas_fractions=gas, state_populations=populations,
               target_fractions=np.asarray(target), product_fractions=np.asarray(product))
    # Rate floors are a collision-set policy, not an interpolated measurement.
    if any(not np.array_equal(r['rate_floors'], row['rate_floors']) for r in rows.flat):
        raise SpecError('source rate floors vary across nodes')
    return row


def _validate_comparator(manifest, rows, held_out):
    """Test the actual runtime interpolant against fresh Maxwellian quadratures.

    Assemble an offline in-memory candidate, never a loadable provisional artifact.
    A failed candidate is not published and cannot bypass runtime validation.
    """
    from rmgpy.solver.eedf import EEDFTable
    from rmgpy.tools.eedf.artifact import _write_rows

    table = EEDFTable()
    table._provider_kind = 'maxwellian-table'
    table._manifest = copy.deepcopy(manifest)
    table._axis_names = tuple(manifest['axes'])
    table._axes = tuple(np.asarray(axis) for axis in manifest['axes'].values())
    table._data, table._layout, table._symbolic = {}, {}, {}
    table._native_field_power, table._native_mobility, table._native_EN_Td = {}, {}, {}
    table._a6b_source = manifest['A6b-source']
    table._t10_interpolation = manifest['T10 interpolation qualification']
    with h5py.File(io.BytesIO(), 'w') as h5:
        group = h5.create_group(COMPARATOR_BRANCH)
        _write_rows(group, rows)
        table._read_branch(COMPARATOR_BRANCH, group)
    if not held_out:
        raise SpecError('comparator requires independent held-out quadratures')
    verdicts = []
    for direct in held_out:
        point = dict(direct['composition'], u=direct['u'])
        predicted = table.row(direct['u'], direct['composition'], COMPARATOR_BRANCH)
        checks = []

        def compare(name, actual, expected, tolerance, absolute=None):
            difference = np.abs(np.asarray(actual) - np.asarray(expected))
            allowed = tolerance['atol'] + tolerance['rtol'] * np.abs(expected)
            if absolute is not None:
                allowed = np.maximum(allowed, absolute)
            fractions = np.divide(difference, allowed, out=np.zeros_like(difference), where=np.asarray(allowed) > 0)
            fractions = np.where((allowed == 0) & (difference > 0), np.inf, fractions)
            maximum = float(np.max(fractions))
            checks.append({'quantity': name, 'max_tolerance_fraction': maximum,
                           'passed': bool(np.isfinite(maximum) and maximum <= 1)})

        tol = manifest['tolerances']
        compare('mean_energy_eV', predicted.swarm['mean_energy_eV'], direct['swarm']['mean_energy_eV'], tol['H1'])
        for name in ('mobility_N', 'diffusion_N', 'energy_mobility_N', 'energy_diffusion_N'):
            compare(name, predicted.swarm[name], direct['swarm'][name], tol['H2'])
        for name in ('k_ine', 'k_sup'):
            compare(name, getattr(predicted, name), direct[name], tol['H3'], direct['rate_floors'])
        compare('attachment_energy_eV', predicted.attachment_energy_eV, direct['attachment_energy_eV'], tol['H2'])
        compare('channel_power', predicted.channel_power, direct['channel_power'], tol['H4'])
        compare('collision_loss', sum(predicted.channel_power), sum(direct['channel_power']),
                {'rtol': tol['H4']['total_rtol'], 'atol': tol['H4']['atol']})
        distance = float(np.sum(np.abs(predicted.f0 - direct['f0']) * np.sqrt(direct['energy_eV'])
                                * np.diff(direct['energy_edges_eV'])))
        compare('f0_weighted_L1', distance, 0., {'rtol': 0., 'atol': tol['F0']['atol']})
        verdict = {'point': point, 'checks': checks, 'passed': all(c['passed'] for c in checks)}
        verdicts.append(verdict)
        if not verdict['passed']:
            failed = [c for c in checks if not c['passed']]
            raise EEDFError('comparator interpolation failed at ' + str(point) + ': ' + str(failed))
    resolution = max(abs(r['swarm']['maxwellian_relative_mean_error'])
                     for r in list(rows.flat) + held_out)
    return {'passed': True, 'training_count': rows.size, 'held_out_count': len(held_out),
            'max_relative_mean_error': float(resolution), 'verdicts': verdicts}


def generate_maxwellian(source_path, output_root, *, artifact_sha256, branch_id,
                        u_start=None, u_stop=None, command_line=None):
    """Publish a separately validated comparator from all or an explicit u slice.

    ``u_start`` and ``u_stop`` are inclusive native-node indices. No automatic
    trimming, grid extension, tolerance relaxation or source rewrite occurs.
    """
    from rmgpy.tools.eedf.artifact import _write_artifact
    from rmgpy.tools.eedf.schema import _interpolant_identity, COMPARATOR_NA
    from rmgpy.tools.eedf.validation import _held_out_points

    source, native = _source_snapshot(source_path, artifact_sha256=artifact_sha256, branch_id=branch_id)
    count = len(source['axes']['u'])
    if (u_start is None) != (u_stop is None):
        raise SpecError('select both inclusive u node indices')
    if u_start is None:
        u_start, u_stop = 0, count - 1
    if (type(u_start) is not int or type(u_stop) is not int
            or not 0 <= u_start < u_stop < count):
        raise SpecError('select at least two contiguous native u nodes')
    axes = copy.deepcopy(source['axes'])
    axes['u'] = axes['u'][u_start:u_stop + 1]
    provenance = {'artifact_sha256': source['artifact_sha256'], 'fingerprint': source['fingerprint'],
                  'branch_id': branch_id, 'source_u': source['axes']['u'],
                  'retained_indices': list(range(u_start, u_stop + 1)), 'retained_u': axes['u'],
                  'excluded_u': source['axes']['u'][:u_start] + source['axes']['u'][u_stop + 1:],
                  'energy_grid_sha256': _content_hash([source['energy_eV'], source['energy_edges_eV']]),
                  'channel_map_sha256': _content_hash(source['channel_map']),
                  'held_out_mean': 'tensor PCHIP of native source means; fresh Maxwellian quadrature'}
    inputs = copy.deepcopy(source['row_inputs'])
    inputs.update(provider_kind='maxwellian-table', comparator_source=provenance,
                  axes=axes, interpolant=_interpolant_identity())
    # Retain source physical input records; identify the new row source explicitly.
    inputs['row_source'] = 'prescribed Maxwellian quadrature; Te_eV = 2*source_mean_energy_eV/3'
    applicability = {name: {'status': 'NOT_APPLICABLE',
                           'reason': 'Prescribed Maxwellian comparator; no self-consistent field or Boltzmann solve.'}
                     for name in COMPARATOR_NA}
    manifest = {key: copy.deepcopy(source[key]) for key in (
        'schema_version', 'energy_eV', 'energy_edges_eV', 'envelopes', 'tolerances', 'floors', 'channel_map')}
    manifest.update(provider_kind='maxwellian-table', row_inputs=inputs, fingerprint=_content_hash(inputs),
                    axes=axes, accepted=False, branches=[COMPARATOR_BRANCH],
                    branch_certification={COMPARATOR_BRANCH: 'NOT_APPLICABLE: prescribed Maxwellian comparator'},
                    branch_detection=applicability['branch_scan'], held_out_verdicts=[],
                    generator_command=command_line,
                    units={'u': 'ln(source E/N / 1 Td); reference coordinate only',
                           'EN_Td': 'source reference Td; not a Maxwellian field',
                           'energy': 'eV', 'k': 'm^3/s', 'power': 'eV*m^3/s'})
    manifest.update({name: applicability[name] for name in (
        'A6b-source', 'A6b-runtime', 'T10 interpolation qualification')})
    manifest['comparator'] = {'source': provenance, 'applicability': applicability,
                              'label': 'MODEL-FORM COMPARATOR: Maxwellian-table energy-balance diagnostic'}
    selected = native[u_start:u_stop + 1]
    rows = np.empty(selected.shape, dtype=object)
    for index in np.ndindex(selected.shape):
        rows[index] = _maxwellian_row(selected[index], source)
    # Keep the source sampling policy, projected to the explicitly selected domain.
    policy = copy.deepcopy(inputs['held_out'])
    policy['off_centre'] = [p for p in policy.get('off_centre', [])
                            if axes['u'][0] <= p['u'] <= axes['u'][-1]]
    points = _held_out_points(axes, policy)
    held_out = [_maxwellian_row(_held_out_source_row(point, source, native), source) for point in points]
    manifest['held_out'] = points
    manifest['comparator']['numerical_validation'] = _validate_comparator(manifest, rows, held_out)
    return _write_artifact(output_root, manifest, {COMPARATOR_BRANCH: rows}, held_out)


def main(argv=None):
    """Offline CLI: source SHA and branch are always explicit."""
    import argparse
    import sys
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source')
    parser.add_argument('output_root')
    parser.add_argument('--artifact-sha256', required=True)
    parser.add_argument('--branch-id', required=True)
    parser.add_argument('--u-start', type=int)
    parser.add_argument('--u-stop', type=int)
    args = parser.parse_args(argv)
    path = generate_maxwellian(args.source, args.output_root, artifact_sha256=args.artifact_sha256,
                               branch_id=args.branch_id, u_start=args.u_start, u_stop=args.u_stop,
                               command_line=list(sys.argv if argv is None else argv))
    print('MODEL-FORM COMPARATOR: Maxwellian-table energy-balance diagnostic')
    print(path)


if __name__ == '__main__':
    main()
