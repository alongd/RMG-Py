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

"""Generate offline LoKI-B EEDF tables: python -m ...generate spec.yml.

Only explicit spec policy controls tolerances, floors, screening and refinement.
Generated tables are outside git. The LoKI binary/build are immutable inputs.
"""
import argparse
import copy
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import itertools
import json
import logging
from pathlib import Path
import sys
import time
import uuid

import numpy as np

from rmgpy.tools.eedf.artifact import write_artifact
from rmgpy.tools.eedf.channels import validate_map as _validate_map, validate_physical_map
from rmgpy.tools.eedf.integrity import check_solver, physical_properties, solver_environment
from rmgpy.tools.eedf.loki import LoKIDriver, enrich_row
from rmgpy.tools.eedf.schema import (
    EEDFError, SpecError, content_hash, file_hash, interpolant_identity,
    load_spec, row_inputs, validate_spec, FingerprintMismatch,
)
from rmgpy.tools.eedf.validation import (
    compare_row, freeze_accuracy_criteria, held_out_points, quantities,
    run_held_out_accuracy, split_branches,
)


class AxisScreenRefusal(EEDFError):
    """A significant coordinate has no axis or its screen is incomplete."""


def model_inputs(spec):
    """Compute current row-input identity after checking every physical file.

    The spec's axes are updated to the selected/refined grid by ``generate``.
    The saved generation_spec.json provides that exact spec for future loads.
    Paths and repository SHAs are excluded; file bytes and solver identity
    own the fingerprint. All property files, not only LXCat, shape rows.
    """
    shared = check_solver(spec)
    physical_properties(spec)
    for key in ('binary', 'cmake_cache', 'channel_map'):
        if file_hash(spec[key]['path']) != spec[key]['sha256']:
            raise FingerprintMismatch(key)
    for name, item in spec['input_files'].items():
        if file_hash(item['path']) != item['sha256']:
            raise FingerprintMismatch('input file ' + name)
    inputs = row_inputs(spec)
    # File closure is hashed above; store flattened constants/placeholders so
    # production loading can independently validate each row's populations.
    _, states = physical_properties(spec)
    inputs['state_properties'] = states
    inputs['solver_environment'] = solver_environment(spec)
    inputs['shared_objects'] = shared
    inputs['binary_sha256'] = file_hash(spec['binary']['path'])
    inputs['cross_sections'] = {name: file_hash(item['path']) for name, item in spec['input_files'].items()
                                if item['kind'] == 'cross_section'}
    inputs['auxiliary_inputs'] = {name: file_hash(item['path']) for name, item in spec['input_files'].items()
                                 if item['kind'] != 'cross_section'}
    inputs['channel_map_sha256'] = file_hash(spec['channel_map']['path'])
    channel_map = json.loads(Path(spec['channel_map']['path']).read_text())
    validate_physical_map(channel_map, spec, {name: values[0] for name, values in spec['axes'].items() if name != 'u'})
    inputs['reactions'] = [item['reaction'] for item in channel_map if item['reaction'] is not None]
    inputs['interpolant'] = interpolant_identity()
    return inputs


def _solve(driver, channel_map, spec, coordinate, fields, label):
    rows = driver.run(coordinate, fields, label)
    for row in rows:
        enrich_row(row, channel_map, spec)
        row['power_absolute'] = abs(row['power_groups']['field']) * spec['floors']['absolute_power_share']
    return rows


def qualify_interpolation(table, spec, *, criteria):
    """Run fresh pinned LoKI-B solves for an existing table's T10 points.

    This is an offline accuracy check, not a runtime identity check. The caller
    must supply frozen mobility and field-power criteria; omitting them refuses
    before any solver is started.
    """
    criteria = freeze_accuracy_criteria(criteria)
    validate_spec(spec)
    fingerprint = content_hash(model_inputs(spec))
    if fingerprint != table.manifest.get('fingerprint'):
        raise FingerprintMismatch('held-out direct-solver setup fingerprint')
    if not spec['held_out'].get('off_centre'):
        raise SpecError(
            'T10 qualification requires preselected held_out off_centre points')
    if len(table.manifest['branches']) != 1:
        raise EEDFError(
            'T10 multi-branch qualification requires a pinned direct-solver '
            'setup for each branch')
    points = held_out_points(table.manifest['axes'], spec['held_out'])
    channel_map = json.loads(Path(spec['channel_map']['path']).read_text())
    _validate_map(channel_map)
    driver = LoKIDriver(spec)
    counter = itertools.count()

    def direct_solver(point):
        coordinate = {
            name: point[name] for name in table.axis_names if name != 'u'}
        rows = _solve(
            driver, channel_map, spec, coordinate, [np.exp(point['u'])],
            't10_held_' + str(next(counter)))
        if len(rows) != 1:
            raise EEDFError('direct held-out solve returned multiple rows')
        return rows[0]

    return run_held_out_accuracy(
        table, direct_solver, points, criteria=criteria,
        setup_fingerprint=fingerprint)


def _screen(spec, driver, channel_map, prefix):
    """Record the frozen rule and direct perturbations over the candidate range.

    Screen every primary u node. Each proposed coordinate must be either a
    screened composition axis or a declared non-axis envelope. Frozen inputs
    with zero-width envelopes require no perturbation claim.
    """
    fields = np.exp(spec['axes']['u']).tolist()
    reference = {item['name']: item['reference'] for item in spec['screen']['candidates']}
    results = []
    baseline = _solve(driver, channel_map, spec, reference, fields, prefix + '_screen_ref')
    for i, candidate in enumerate(spec['screen']['candidates']):
        maximum, worst = 0., None
        for j, value in enumerate(candidate['values']):
            coordinates = dict(reference, **{candidate['name']: value})
            rows = _solve(driver, channel_map, spec, coordinates, fields, f'{prefix}_screen_{i}_{j}')
            for a, b in zip(baseline, rows):
                verdicts = compare_row(b, a, {'tolerances': spec['tolerances'], 'floors': spec['floors'], 'channel_map': channel_map, 'row_inputs': {'Tg_K': spec['Tg_K']}})
                for verdict in verdicts:
                    if verdict['rule'] in ('consistency', 'channel_power_sum'):
                        continue
                    ratio = (verdict['absolute_error'] / verdict['allowed_error'] if verdict['allowed_error']
                             else (float('inf') if verdict['absolute_error'] else 0.))
                    if ratio > maximum:
                        maximum, worst = ratio, {'u': a['u'], 'value': value, 'quantity': verdict['quantity']}
        significant = maximum > spec['screen']['threshold_fraction']
        name = candidate['name']
        if significant and name not in spec['axes']:
            raise AxisScreenRefusal(name + ' moves a quantity beyond the frozen screen threshold; supply its axis')
        if not significant and name in spec['axes']:
            if name not in spec['envelopes']:
                raise AxisScreenRefusal(name + ' needs a refusing envelope before dropping its axis')
            del spec['axes'][name]
        if significant and (min(candidate['values']) > spec['axes'][name][0] or max(candidate['values']) < spec['axes'][name][-1]):
            raise AxisScreenRefusal('screen does not cover axis ' + name)
        if not significant:
            bounds = spec['envelopes'][name]
            if min(candidate['values']) > bounds['min'] or max(candidate['values']) < bounds['max']:
                raise AxisScreenRefusal('screen does not cover envelope ' + name)
        results.append({'coordinate': name, 'reference': candidate['reference'], 'values': candidate['values'],
                        'max_fraction_of_held_out_tolerance': maximum, 'worst': worst,
                        'selected_axis': significant})
    screened = {result['coordinate'] for result in results}
    if any(bounds['min'] != bounds['max'] and name not in screened for name, bounds in spec['envelopes'].items()):
        raise AxisScreenRefusal('nonzero-width envelope has no sensitivity screen')
    if set(spec['axes']) - {'u'} - screened:
        raise AxisScreenRefusal('unscreened composition axes')
    return results


def _manifest(spec, driver, channel_map, screen_results, scans, row, command_line):
    inputs = model_inputs(spec)
    return {'schema_version': 1, 'row_inputs': inputs, 'fingerprint': content_hash(inputs),
            'axes': spec['axes'], 'energy_eV': row['energy_eV'].tolist(),
            'energy_edges_eV': row['energy_edges_eV'].tolist(), 'envelopes': spec['envelopes'],
            'screen_results': screen_results, 'tolerances': spec['tolerances'], 'floors': spec['floors'],
            'channel_map': channel_map, 'repositories': spec['repositories'],
            'solver': {'commit': spec['loki_commit'], 'binary_sha256': spec['binary']['sha256'],
                       'compiler': spec['compiler'], 'cmake_cache_sha256': spec['cmake_cache']['sha256'],
                       'options': spec['solver_options'], 'ee_setting': spec['solver_options']['includeEECollisions'],
                       'working_conditions': spec['working_conditions'], 'input_files': {
                           name: {'sha256': item['sha256'], 'kind': item['kind']} for name, item in spec['input_files'].items()}},
            'setups': dict(driver.setups), 'state_population_assumptions': spec['state_properties'],
            'density_convention': 'N = P/(k_B*Tg); P in Pa, Tg in K; electronDensity in m^-3',
            'units': {'u': 'ln(E/N / 1 Td)', 'EN_Td': 'Td', 'energy': 'eV', 'k': 'm^3/s',
                      'power': 'eV*m^3/s; signed LoKI groups, positive channel loss',
                      'mobility_N': '1/(m*s*V)', 'diffusion_N': '1/(m*s)',
                      'energy_mobility_N': 'eV/(m*s*V)', 'energy_diffusion_N': 'eV/(m*s)',
                      'townsend_N': 'm^2', 'attachment_N': 'm^2', 'f0': 'eV^(-3/2)'},
            'generator_command': command_line, 'generated_at': datetime.now(timezone.utc).isoformat(),
            'generation_seconds': 0., 'held_out': [], 'held_out_verdicts': [], 'accepted': False,
            'branch_detection': scans, 'refinement_history': [],
            'branch_certification': {branch: 'uncertified: unseeded scans' for branch in scans[0]['branch_sources']},
            'limitations': ['Pinned temporal-growth binary does not expose iteration count; -1 means unknown.',
                            'Pinned solver starts each solve with invertLinearMatrix; ordered scans do not prove seeded continuation.',
                            'This generator enumerates observed EEDF scan branches, not stability or P_abs continuation of reactor branches.']}


def generate(spec, command_line=None):
    """Generate, independently qualify and publish a content-addressed table.

    Failing cells get third-point nodes; old held-out points remain held out.
    A failed verdict publishes an explicitly unaccepted artifact and every
    failing quantity. No tolerance is relaxed. The caller's grid is updated
    to the exact screened/refined grid, for current-model fingerprint checks.
    """
    from rmgpy.solver.eedf import EEDFTable
    validate_spec(spec)
    started = time.monotonic()
    prefix = 'run_' + uuid.uuid4().hex
    channel_map = json.loads(Path(spec['channel_map']['path']).read_text())
    _validate_map(channel_map)
    driver = LoKIDriver(spec)
    spec['axes'] = {'u': spec['axes']['u'], **{name: spec['axes'][name] for name in sorted(spec['axes']) if name != 'u'}}
    screen_results = _screen(spec, driver, channel_map, prefix)
    held = {}
    history = []
    for round_index in range(spec['refinement']['max_rounds'] + 1):
        axes = spec['axes']
        names = list(axes)
        shape = tuple(len(axis) for axis in axes.values())
        fields = np.exp(axes['u']).tolist()
        composition_grid = list(itertools.product(*[axes[name] for name in names[1:]]))
        scans, per_node = [], []
        # Each job is isolated. This bound applies to all solves, including cold and held-out.
        with ThreadPoolExecutor(max_workers=spec['max_parallel']) as executor:
            pending = []
            for i, composition in enumerate(composition_grid):
                coordinate = dict(zip(names[1:], composition))
                base = f'{prefix}_r{round_index}_c{i}'
                pending.append((coordinate,
                    executor.submit(_solve, driver, channel_map, spec, coordinate, fields, base + '_up'),
                    executor.submit(_solve, driver, channel_map, spec, coordinate, fields[::-1], base + '_down'),
                    [executor.submit(_solve, driver, channel_map, spec, coordinate, [field], base + f'_cold{j}')
                     for j, field in enumerate(fields)]))
            for coordinate, up, down, cold in pending:
                branches, verdict = split_branches({'up': up.result(), 'down': down.result(),
                    'cold': [future.result()[0] for future in cold]}, spec['tolerances']['G1'])
                verdict['composition'] = coordinate
                scans.append(verdict)
                per_node.append(branches)
        branch_ids = set(per_node[0])
        if any(scan['branch_sources'] != scans[0]['branch_sources'] for scan in scans):
            raise EEDFError('ambiguous branch correspondence across composition; split the domain')
        if any(set(branches) != branch_ids for branches in per_node):
            raise EEDFError('branch count varies across composition; refine/split domain before tensor interpolation')
        tensors = {branch: np.empty(shape, dtype=object) for branch in sorted(branch_ids)}
        for composition_index, branches in enumerate(per_node):
            index = np.unravel_index(composition_index, shape[1:]) if len(shape) > 1 else ()
            for branch_id, rows in branches.items():
                for u_index, row in enumerate(rows):
                    row['branch_id'] = branch_id
                    tensors[branch_id][(u_index,) + index] = row
        points = held_out_points(axes, spec['held_out'])
        training = set(itertools.product(*axes.values()))
        if training.intersection(held):
            raise EEDFError('refinement would promote a held-out point')
        new_points = [point for point in points if tuple(point[name] for name in names) not in held]
        with ThreadPoolExecutor(max_workers=spec['max_parallel']) as executor:
            pending = [executor.submit(_solve, driver, channel_map, spec,
                {name: point[name] for name in names[1:]}, [np.exp(point['u'])],
                f'{prefix}_r{round_index}_held{i}') for i, point in enumerate(new_points)]
            for point, future in zip(new_points, pending):
                key = tuple(point[name] for name in names)
                held[key] = future.result()[0]
        sample = next(iter(tensors.values())).flat[0]
        manifest = _manifest(spec, driver, channel_map, screen_results, scans, sample, command_line or sys.argv)
        manifest['held_out'] = [dict(zip(names, key)) for key in held]
        provisional = write_artifact(Path(spec['scratch_root']) / f'{prefix}_r{round_index}_table', manifest, tensors, list(held.values()))
        table = EEDFTable.load(provisional, model_inputs(spec), artifact_sha256=file_hash(provisional / 'table.h5'), require_accepted=False)
        manifest['A6b-source'] = table.a6b_source
        manifest['A6b-runtime'] = {
            'check': 'A6b-runtime',
            'criterion': 'field power derived from interpolated mobility at query state',
            'passed_by_construction': True,
            'accuracy_claim': False,
        }
        failures, verdicts = [], []
        for key, direct in held.items():
            point = dict(zip(names, key))
            for branch_id in sorted(branch_ids):
                try:
                    predicted = table.row(point['u'], {name: point[name] for name in names[1:]}, branch_id).as_dict()
                    checks = compare_row(predicted, direct, manifest)
                    passed = all(check['passed'] for check in checks)
                    result = {'point': point, 'branch_id': branch_id, 'passed': passed, 'checks': checks}
                except EEDFError as error:
                    passed = False
                    result = {'point': point, 'branch_id': branch_id, 'passed': False, 'refusal': str(error), 'checks': []}
                verdicts.append(result)
                if not passed:
                    failures.append(point)
        history.append({'round': round_index, 'axes': copy.deepcopy(axes), 'held_out_count': len(held),
                        'failing_points': failures})
        manifest['held_out_verdicts'] = verdicts
        manifest['T10 interpolation qualification'] = {
            'check': 'T10 interpolation qualification',
            'criterion': 'separately frozen mobility and field-power accuracy criteria',
            'verdict_count': 0,
            'passed': False,
            'blocker': 'T10 accuracy criteria have not been supplied',
        }
        manifest['accepted'] = False
        manifest['refinement_history'] = history
        if not failures or round_index == spec['refinement']['max_rounds']:
            manifest['generation_seconds'] = time.monotonic() - started
            final = write_artifact(spec['output_root'], manifest, tensors, list(held.values()),
                                   extra_files={'generation_spec.json': json.dumps(spec, indent=2, allow_nan=False) + '\n',
                                                'validation.md': _validation_report(manifest)})
            logging.info('Generated %s; accepted=%s; %.3f s', final, manifest['accepted'], manifest['generation_seconds'])
            return final
        # Refine only failing cells, without promoting any existing independent point.
        for name in names:
            nodes = set(axes[name])
            for point in failures:
                axis = np.asarray(axes[name])
                i = min(max(np.searchsorted(axis, point[name], side='right') - 1, 0), len(axis) - 2)
                nodes.update([float(axis[i] + (axis[i+1] - axis[i]) / 3),
                              float(axis[i] + 2 * (axis[i+1] - axis[i]) / 3)])
            axes[name] = sorted(nodes)
    raise AssertionError('unreachable')


def _validation_report(manifest):
    lines = ['# EEDF held-out validation', '', 'Accepted: ' + str(manifest['accepted']),
             'Generation seconds: ' + str(manifest['generation_seconds']), '',
             'A6b-source: ' + ('PASS' if manifest['A6b-source']['passed'] else 'FAIL'),
             'A6b-runtime: PASS by construction (not an accuracy claim)',
             'T10 interpolation qualification: ' + (
                 'PASS' if manifest['T10 interpolation qualification']['passed']
                 else 'FAIL'), '',
             '| rule | checks | failures | maximum relative error |', '|---|---:|---:|---:|']
    for rule in ('H1', 'H2', 'H3', 'H4', 'F0', 'consistency', 'channel_power_sum'):
        checks = [check for result in manifest['held_out_verdicts'] for check in result['checks'] if check['rule'] == rule]
        relative = [check['relative_error'] for check in checks if check['relative_error'] is not None]
        lines.append(f"| {rule} | {len(checks)} | {sum(not check['passed'] for check in checks)} | {max(relative, default=0):.6g} |")
    lines.extend(['', '## Failing cells', ''])
    failures = [v for v in manifest['held_out_verdicts'] if not v['passed']]
    for result in failures:
        lines.append(canonical_failure(result))
    if not failures:
        lines.append('None.')
    lines.extend(['', '## Limits', ''] + manifest['limitations'])
    return '\n'.join(lines) + '\n'


def canonical_failure(result):
    return json.dumps({'point': result['point'], 'branch_id': result['branch_id'],
                       'refusal': result.get('refusal'), 'failing_checks': [v for v in result['checks'] if not v['passed']]}, sort_keys=True)


def main():
    """Run the generation CLI; unsuccessful qualification returns nonzero."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('spec')
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    spec = load_spec(args.spec)
    path = generate(spec)
    manifest = json.loads((path / 'manifest.json').read_text())
    return 0 if manifest['accepted'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
