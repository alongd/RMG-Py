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

"""Scan agreement, independent held-out points and validation verdicts."""
import itertools

import numpy as np
from scipy.stats import qmc

from rmgpy.tools.eedf.schema import EEDFError, FingerprintMismatch, SpecError


def quantities(row):
    """Numeric quantities shared by scan, sensitivity and held-out comparisons."""
    values = {'k_ine': np.asarray(row['k_ine']), 'k_sup': np.asarray(row['k_sup']),
              'channel_power': np.asarray(row['channel_power']),
              'attachment_energy_eV': np.asarray(row['attachment_energy_eV']),
              'f0': np.asarray(row['f0'])}
    values.update({'swarm.' + key: np.asarray(value) for key, value in row['swarm'].items()})
    values.update({'power.' + key: np.asarray(value) for key, value in row['power_groups'].items()})
    return values


def rows_agree(a, b, tolerance):
    """All row quantities, with rate/power absolute floors where available."""
    for name in ('energy_eV', 'energy_edges_eV'):
        if name in a or name in b:
            if name not in a or name not in b or not np.array_equal(a[name], b[name]):
                raise EEDFError('energy grid changes across solves')
    qa, qb = quantities(a), quantities(b)
    if set(qa) != set(qb):
        return False
    for key in qa:
        absolute = tolerance['atol']
        if key in ('k_ine', 'k_sup') and 'rate_floors' in a:
            absolute = np.maximum(absolute, np.maximum(a['rate_floors'], b['rate_floors']))
        if key == 'channel_power' and 'power_absolute' in a:
            absolute = max(absolute, a['power_absolute'], b['power_absolute'])
        if not np.all(abs(qa[key] - qb[key]) <= absolute + tolerance['rtol'] * np.maximum(abs(qa[key]), abs(qb[key]))):
            return False
    return True


def split_branches(passes, tolerance):
    """Keep distinct ordered scan trajectories; never average different EEDFs.

    Down-scan rows enter in descending order. A branch is a full observed scan
    trajectory. Stability of a coupled reactor is outside this scan's meaning.
    """
    ordered = {'up': passes['up'], 'down': passes['down'][::-1], 'cold': passes['cold']}
    if len({len(v) for v in ordered.values()}) != 1:
        raise EEDFError('incomplete scan path')
    branches, sources = {}, {}
    for name, rows in ordered.items():
        for branch_id, existing in branches.items():
            if all(rows_agree(a, b, tolerance) for a, b in zip(rows, existing)):
                sources[branch_id].append(name)
                break
        else:
            branch_id = 'branch_' + str(len(branches))
            branches[branch_id] = rows
            sources[branch_id] = [name]
    disagreements = [i for i in range(len(ordered['up']))
                     if not all(rows_agree(ordered['up'][i], ordered[name][i], tolerance)
                                for name in ('down', 'cold'))]
    return branches, {'agreement': not disagreements, 'disagreeing_nodes': disagreements,
                      'branch_sources': sources}


def held_out_points(axes, policy):
    """Every cell midpoint plus frozen off-centre and seeded points."""
    names = list(axes)
    grid = [np.asarray(v, dtype=float) for v in axes.values()]
    points = set()
    if policy['cell_midpoints']:
        for dimension, axis in enumerate(grid):
            options = list(grid)
            options[dimension] = (axis[:-1] + axis[1:]) / 2
            points.update(itertools.product(*options))
        points.update(itertools.product(*[(a[:-1] + a[1:]) / 2 for a in grid]))
    lhs = qmc.LatinHypercube(len(grid), seed=policy['seed']).random(policy['lhs_count'])
    if len(lhs):
        samples = qmc.scale(lhs, [v[0] for v in grid], [v[-1] for v in grid])
        points.update(map(tuple, samples))
    for point in policy.get('off_centre', []):
        if set(point) != set(names):
            raise SpecError('held_out off_centre coordinates')
        values = tuple(float(point[name]) for name in names)
        if any(not axis[0] <= value <= axis[-1]
               for axis, value in zip(grid, values)):
            raise SpecError('held_out off_centre domain')
        points.add(values)
    training = set(itertools.product(*grid))
    return [dict(zip(names, map(float, p))) for p in sorted(points - training)]


def freeze_accuracy_criteria(criteria):
    """Validate and copy the separately frozen T10 transport criteria."""
    if not isinstance(criteria, dict) or set(criteria) != {
            'mobility_N', 'field_power'}:
        raise SpecError('held-out accuracy criteria require mobility_N and field_power')
    frozen = {}
    for name, tolerance in criteria.items():
        if (not isinstance(tolerance, dict) or set(tolerance) != {'rtol', 'atol'}
                or any(not isinstance(value, (int, float)) or
                       not np.isfinite(value) or value < 0.
                       for value in tolerance.values())):
            raise SpecError('held-out accuracy criterion ' + name)
        frozen[name] = {key: float(value) for key, value in tolerance.items()}
    return frozen


def run_held_out_accuracy(table, direct_solver, points, *, criteria,
                          setup_fingerprint):
    """Compare table predictions with fresh direct-solver rows.

    ``criteria`` is deliberately mandatory and separate from A6b's identity
    tolerance. The direct-solver callback is invoked once per point; references
    are never obtained from the table under test.
    """
    frozen = freeze_accuracy_criteria(criteria)
    if not setup_fingerprint or setup_fingerprint != table.manifest.get('fingerprint'):
        raise FingerprintMismatch('held-out direct-solver setup fingerprint')

    names = tuple(table.axis_names)
    results = []

    def record(quantity, actual, reference, tolerance):
        difference = abs(float(actual) - float(reference))
        allowed = tolerance['atol'] + tolerance['rtol'] * abs(float(reference))
        return {
            'rule': 'T10',
            'quantity': quantity,
            'index': [],
            'absolute_error': difference,
            'relative_error': (difference / abs(float(reference))
                               if reference else None),
            'allowed_error': allowed,
            'passed': bool(difference <= allowed),
        }

    for point in points:
        if not isinstance(point, dict) or set(point) != set(names):
            raise SpecError('held-out point coordinates')
        composition = {name: point[name] for name in names[1:]}
        table.domain_check(point['u'], composition)
        direct = direct_solver(dict(point))
        if not isinstance(direct, dict):
            raise EEDFError('direct solver did not return a row mapping')
        predicted = table.row(point['u'], composition).as_dict()
        checks = compare_row(predicted, direct, table.manifest)
        checks.append(record(
            'mobility_N', predicted['swarm']['mobility_N'],
            direct['swarm']['mobility_N'], frozen['mobility_N']))
        checks.append(record(
            'field_power', predicted['power_groups']['field'],
            direct['power_groups']['field'], frozen['field_power']))
        results.append({
            'point': {name: float(point[name]) for name in names},
            'branch_id': predicted['branch_id'],
            'reference': 'fresh direct solver output',
            'criteria': frozen,
            'checks': checks,
            'passed': bool(checks) and all(check['passed'] for check in checks),
        })
    return results


def compare_row(predicted, direct, manifest):
    """Return named H1–H4 verdicts, using species flux and field power scales.

    Resolved target fractions supply conservative individual flux budgets. Rates
    below the dynamic-range floor use absolute tests. Declared flux groups
    cannot increase an allowance.
    No reactor flux or materiality budget is invented by this routine.
    """
    tol, floors = manifest['tolerances'], manifest['floors']
    for name in ('energy_eV', 'energy_edges_eV'):
        if not np.array_equal(predicted[name], direct[name]):
            raise EEDFError('held-out energy grid changes: ' + name)
    verdicts = []

    def record(rule, quantity, actual, reference, rtol, atol):
        actual, reference = np.asarray(actual), np.asarray(reference)
        difference = abs(actual - reference)
        allowed = np.asarray(atol) + rtol * abs(reference)
        passed = difference <= allowed
        for index in np.ndindex(difference.shape):
            ref, diff, allowance = float(reference[index]), float(difference[index]), float(allowed[index] if allowed.shape else allowed)
            verdicts.append({'rule': rule, 'quantity': quantity, 'index': list(index),
                             'absolute_error': diff, 'relative_error': diff / abs(ref) if ref else None,
                             'allowed_error': allowance, 'passed': bool(passed[index])})

    record('H1', 'mean_energy_eV', predicted['swarm']['mean_energy_eV'], direct['swarm']['mean_energy_eV'], **tol['H1'])
    for name, reference in direct['swarm'].items():
        if name != 'mean_energy_eV':
            record('H2', name, predicted['swarm'][name], reference, **tol['H2'])
    channels = manifest['channel_map']
    for direction in ('k_ine', 'k_sup'):
        targets = np.asarray(direct['target_fractions'] if direction == 'k_ine' else direct.get('product_fractions', np.zeros(len(channels))))
        fluxes = targets * direct[direction]
        for j, channel in enumerate(channels):
            # A declared flux-group label cannot increase an absolute budget.
            # Individual flux is a conservative bound on any species total.
            total = fluxes[j]
            # A rate coefficient has its own numerical scale even when its
            # target population (and therefore its flux) is zero.  Flux is
            # relevant to absolute budgets, but cannot suppress the declared
            # relative H3 tolerance for a coefficient above its own floor.
            significant = not direct['below_floor'][j, 0 if direction == 'k_ine' else 1]
            absolute = max(tol['H3']['atol'], direct['rate_floors'][j],
                           floors['absolute_flux_fraction'] * total / targets[j] if targets[j] else 0)
            record('H3', direction + ':' + str(j), predicted[direction][j], direct[direction][j],
                   tol['H3']['rtol'] if significant else 0, absolute)
    field = abs(direct['power_groups']['field'])
    for j, reference in enumerate(direct['channel_power']):
        significant = abs(reference) >= floors['relative_power_share'] * field
        record('H4', 'channel_power:' + str(j), predicted['channel_power'][j], reference,
               tol['H4']['rtol'] if significant else 0,
               max(tol['H4']['atol'], floors['absolute_power_share'] * field))
    total_keys = ['elastic_gain', 'elastic_loss', 'car_gain', 'car_loss', 'excitation_loss', 'excitation_gain',
                  'vibrational_loss', 'vibrational_gain', 'rotational_loss', 'rotational_gain', 'ionization', 'attachment']
    for name in total_keys + ['field', 'growth']:
        record('H4', 'power.' + name, predicted['power_groups'][name], direct['power_groups'][name],
               tol['H4']['rtol'], max(tol['H4']['atol'], floors['absolute_power_share'] * field))
    record('H4', 'total_loss', -sum(predicted['power_groups'][k] for k in total_keys),
           -sum(direct['power_groups'][k] for k in total_keys), tol['H4']['total_rtol'],
           max(tol['H4']['atol'], floors['absolute_power_share'] * field))
    # Distribution and attachment moments were previously absent from H1-H4.
    weights = np.sqrt(direct['energy_eV']) * np.diff(direct['energy_edges_eV'])
    distance = np.sum(abs(np.asarray(predicted['f0']) - direct['f0']) * weights)
    record('F0', 'f0.weighted_L1', distance, 0., 0., tol['F0']['atol'])
    record('F0', 'attachment_energy_eV', predicted['attachment_energy_eV'], direct['attachment_energy_eV'], **tol['H1'])
    from rmgpy.tools.eedf.moments import distribution_moments
    physical = dict(predicted, gas_temperature_K=predicted.get('composition', {}).get('Tg_K', manifest.get('row_inputs', {}).get('Tg_K', manifest.get('Tg_K'))))
    moments = distribution_moments(physical, manifest['channel_map'])
    record('consistency', 'f0.mean_energy_eV', predicted['swarm']['mean_energy_eV'], moments['mean_energy_eV'], **tol['H1'])
    for direction in ('k_ine', 'k_sup'):
        record('consistency', 'f0.' + direction, predicted[direction], moments[direction],
               tol['H3']['rtol'], np.maximum(tol['H3']['atol'], predicted['rate_floors']))
    field = abs(predicted['power_groups']['field'])
    record('consistency', 'f0.channel_power', predicted['channel_power'], moments['channel_power'],
           tol['H4']['rtol'], max(tol['H4']['atol'], floors['absolute_power_share'] * field))
    verdicts.extend(channel_power_checks(predicted, manifest))
    verdicts.extend(channel_power_checks(direct, manifest))
    return verdicts


def channel_power_checks(row, manifest):
    """Compare group sums with a material-channel and field-scaled error cap.

    A fixed absolute floor must not hide a missing channel at low fields. The
    relative scale is capped by the materiality threshold and smallest material
    channel in each group,
    so a large group cannot hide a large fractional error in that channel.
    All numerical budgets remain explicit spec/manifest fields.
    """
    groups = {'elastic': ('elastic_gain', 'elastic_loss'),
              'excitation': ('excitation_loss', 'excitation_gain'),
              'vibrational': ('vibrational_loss', 'vibrational_gain'),
              'rotational': ('rotational_loss', 'rotational_gain'),
              'ionization': ('ionization',), 'attachment': ('attachment',)}
    tolerance = manifest['tolerances']['channel_power_sum']
    verdicts = []
    for kind, names in groups.items():
        actual = sum(row['channel_power'][i] for i, channel in enumerate(manifest['channel_map']) if channel['kind'] == kind)
        reference = -sum(row['power_groups'].get(name, 0.) for name in names)
        powers = [abs(row['channel_power'][i]) for i, channel in enumerate(manifest['channel_map'])
                  if channel['kind'] == kind]
        field = abs(row['power_groups']['field'])
        floors = manifest['floors']
        material = [power for power in powers if power > 0 and power >= floors['relative_power_share'] * field]
        # Include the materiality threshold itself: an undercounted channel
        # might have fallen below it, and must not escape the error budget.
        relative_scale = min([abs(reference), floors['relative_power_share'] * field] + material)
        absolute = min(tolerance['atol'], floors['absolute_power_share'] * field)
        allowed = absolute + tolerance['rtol'] * relative_scale
        difference = abs(actual - reference)
        verdicts.append({'rule': 'channel_power_sum', 'quantity': 'power.sum.' + kind,
                         'index': [], 'absolute_error': float(difference),
                         'relative_error': float(difference / abs(reference)) if reference else None,
                         'allowed_error': float(allowed),
                         'passed': bool(np.isfinite(actual) and np.isfinite(reference) and difference <= allowed)})
    return verdicts
