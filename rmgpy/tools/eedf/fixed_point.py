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

"""Adaptive EEDF table/final-state fixed-point refinement."""

import copy
import json
from pathlib import Path

import numpy as np

from rmgpy.tools.eedf.generate import generate
from rmgpy.tools.eedf.schema import EEDFError, file_hash


class FixedPointRefusal(EEDFError):
    """The table/reactor iteration did not reach the declared fixed point."""


def _read(value):
    if isinstance(value, (str, Path)):
        path = Path(value)
        return json.loads(path.read_text()), path
    return copy.deepcopy(value), None


def _check(name, previous, current, tolerance, labels):
    previous = np.asarray(previous, dtype=float)
    current = np.asarray(current, dtype=float)
    if previous.shape != current.shape:
        raise FixedPointRefusal(name + ' operating-point shape changed')
    error = np.abs(current - previous)
    if set(tolerance) == {'rtol', 'atol'}:
        rtol = np.full(previous.shape, float(tolerance['rtol']))
        atol = np.full(previous.shape, float(tolerance['atol']))
    else:
        missing = set(labels) - set(tolerance)
        extra = set(tolerance) - set(labels)
        if missing or extra:
            raise FixedPointRefusal(
                name + ' quantity tolerances differ: missing {0}, extra {1}'.format(
                    sorted(missing), sorted(extra)))
        rtol = np.asarray([float(tolerance[label]['rtol']) for label in labels])
        atol = np.asarray([float(tolerance[label]['atol']) for label in labels])
    allowed = atol + rtol * np.abs(previous)
    if not error.size:
        raise FixedPointRefusal(name + ' has no declared quantities')
    index = int(np.argmax(error - allowed))
    return {
        'rule': name,
        'absolute_error': float(error.flat[index]),
        'allowed_error': float(allowed.flat[index]),
        'argmax': labels[index],
        'passed': bool(np.all(np.isfinite(error)) and np.all(error <= allowed)),
    }


def _operating_point(record):
    state = record['terminal_state']
    observables = state.get('observables', {})
    required = ('mean_energy_eV', 'electron_density')
    if any(name not in observables for name in required):
        raise FixedPointRefusal('terminal state lacks F2 mean energy/electron density')
    headline = state.get('headline_observables', {})
    if not headline:
        raise FixedPointRefusal('terminal state lacks F3 headline observables')
    return {
        'u': float(state['u']),
        'composition_names': sorted(state.get('composition', {})),
        'composition': [float(state['composition'][name])
                        for name in sorted(state.get('composition', {}))],
        'F2_names': list(required),
        'F2': [float(observables[name]) for name in required],
        'F3_names': sorted(headline),
        'F3': [float(headline[name]) for name in sorted(headline)],
    }


def _extend_axis(nodes, value):
    nodes = sorted(float(item) for item in nodes)
    if len(nodes) < 2 or not nodes[0] <= value <= nodes[-1]:
        raise FixedPointRefusal('final point lies outside table axis')
    if value == nodes[0]:
        spacing = nodes[1] - value
        return sorted(set(nodes + [value - spacing / 2., value + spacing / 2.]))
    if value == nodes[-1]:
        spacing = value - nodes[-2]
        return sorted(set(nodes + [value - spacing / 2., value + spacing / 2.]))
    upper = min(max(int(np.searchsorted(nodes, value, side='right')), 1), len(nodes) - 1)
    lower = upper - 1
    additions = [value, (nodes[lower] + value) / 2., (value + nodes[upper]) / 2.]
    return sorted(set(nodes + additions))


def _extend_composition_axis(nodes, value):
    """Refine a composition axis without leaving its physical domain."""
    nodes = sorted(float(item) for item in nodes)
    value = float(value)
    if (not nodes or any(node < 0. or node > 1. for node in nodes) or
            value < 0. or value > 1.):
        raise FixedPointRefusal('composition axis exceeds physical bounds [0, 1]')
    lower, upper = nodes[0], nodes[-1]
    return [node for node in _extend_axis(nodes, value)
            if 0. <= node <= 1. and lower <= node <= upper]


def refine_until_converged(failure_record, rerun, *, tolerances,
                           generator=generate, output_path=None):
    """Regenerate and re-run until direct PASS and F1--F3 all agree.

    ``rerun(artifact_path, previous_y)`` must return a qualification record.
    F1--F4 are caller-supplied; no fixed-point tolerance is embedded here.
    """
    initial, source = _read(failure_record)
    if initial.get('status') != 'FAIL':
        raise FixedPointRefusal('fixed-point refinement requires a FAIL record')
    missing = [name for name in ('F1', 'F2', 'F3', 'F4') if name not in tolerances]
    if missing:
        raise FixedPointRefusal('missing fixed-point tolerance ' + ', '.join(missing))
    spec_name = initial.get('generation_spec')
    if spec_name is None and source is not None:
        candidate = source.parent / 'generation_spec.json'
        if candidate.is_file():
            spec_name = str(candidate)
    if spec_name is None:
        raise FixedPointRefusal('FAIL record has no generation_spec')
    spec = json.loads(Path(spec_name).read_text())
    previous = initial
    previous_point = _operating_point(previous)
    previous_y = previous['terminal_state'].get('y')
    iterations = []
    limit = int(tolerances['F4']['max_iterations'])
    if limit < 1:
        raise FixedPointRefusal('F4 max_iterations must be positive')

    for iteration in range(1, limit + 1):
        point = _operating_point(previous)
        refined = copy.deepcopy(spec)
        refined['axes']['u'] = _extend_axis(refined['axes']['u'], point['u'])
        composition = dict(zip(point['composition_names'], point['composition']))
        missing_axes = set(refined['axes']) - {'u'} - set(composition)
        if missing_axes:
            raise FixedPointRefusal(
                'terminal state lacks composition axis: ' + ', '.join(sorted(missing_axes)))
        for name in set(refined['axes']) - {'u'}:
            value = composition[name]
            refined['axes'][name] = _extend_composition_axis(
                refined['axes'][name], value)
        refined['fixed_point_version'] = int(refined.get('fixed_point_version', 0)) + 1
        artifact = Path(generator(refined))
        current = rerun(artifact, previous_y)
        current_point = _operating_point(current)
        if current_point['composition_names'] != previous_point['composition_names']:
            raise FixedPointRefusal('composition identities changed across re-runs')
        if current_point['F2_names'] != previous_point['F2_names']:
            raise FixedPointRefusal('F2 identities changed across re-runs')
        if current_point['F3_names'] != previous_point['F3_names']:
            raise FixedPointRefusal('F3 identities changed across re-runs')
        checks = [
            _check('F1', [np.exp(previous_point['u'])],
                   [np.exp(current_point['u'])], tolerances['F1'], ['E/N_Td']),
            _check('F2', previous_point['F2'], current_point['F2'],
                   tolerances['F2'], previous_point['F2_names']),
            _check('F3', previous_point['F3'], current_point['F3'],
                   tolerances['F3'], previous_point['F3_names']),
        ]
        iterations.append({
            'iteration': iteration,
            'table_version': refined['fixed_point_version'],
            'artifact_path': str(artifact),
            'artifact_sha256': (file_hash(artifact / 'table.h5')
                                if (artifact / 'table.h5').is_file() else None),
            'qualification_status': current.get('status'),
            'drift_maxima': current.get('drift_maxima', {}),
            'agreement': checks,
            'operating_point': current_point,
        })
        result = {
            'status': 'PASS' if (current.get('status') == 'PASS' and
                                 all(check['passed'] for check in checks)) else 'REFINING',
            'iteration_count': iteration,
            'iterations': iterations,
            'final_record': current,
        }
        if result['status'] == 'PASS':
            if output_path is not None:
                Path(output_path).write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
            return result
        spec = refined
        previous = current
        previous_point = current_point
        previous_y = current['terminal_state'].get('y')

    result['status'] = 'REFUSED'
    result['refusal'] = 'F4 maximum iteration count exceeded'
    if output_path is not None:
        Path(output_path).write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    raise FixedPointRefusal(result['refusal'])
