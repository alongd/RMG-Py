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

import json
from pathlib import Path

import pytest

from rmgpy.tools.eedf.fixed_point import FixedPointRefusal, refine_until_converged


def failure_record(tmp_path):
    spec = tmp_path / 'generation_spec.json'
    spec.write_text(json.dumps({
        'axes': {'u': [0., 1., 2.], 'x': [0., .5, 1.]},
        'fixed_point_version': 0,
    }))
    return {
        'status': 'FAIL',
        'generation_spec': str(spec),
        'terminal_state': {
            'u': .8,
            'composition': {'x': .4},
            'observables': {'mean_energy_eV': 2., 'electron_density': 3.},
            'headline_observables': {'Ar+': .2},
            'y': [1., 2.],
        },
        'drift_maxima': {'A3': {'argmax': 'Ar excitation'}},
    }


def tolerances(limit=4):
    return {
        'F1': {'rtol': 0., 'atol': .03},
        'F2': {
            'mean_energy_eV': {'rtol': 0., 'atol': .02},
            'electron_density': {'rtol': .01, 'atol': 0.},
        },
        'F3': {'Ar+': {'rtol': 0., 'atol': .05}},
        'F4': {'max_iterations': limit},
    }


def synthetic_artifact(tmp_path, name):
    path = tmp_path / name
    path.mkdir()
    (path / 'table.h5').write_bytes(b'synthetic non-empty EEDF artifact')
    return path


def test_fixed_point_refines_cells_and_converges_on_smooth_arm(tmp_path):
    generated = []
    previous_states = []

    def generator(spec):
        generated.append(spec)
        return synthetic_artifact(
            tmp_path, 'artifact-' + str(len(generated)))

    def rerun(artifact, previous_y):
        iteration = len(generated)
        previous_states.append(previous_y)
        residual = .1 ** iteration
        return {
            'status': 'PASS',
            'terminal_state': {
                'u': .8 + residual,
                'composition': {'x': .4 + residual},
                'observables': {
                    'mean_energy_eV': 2. + residual,
                    'electron_density': 3. + residual,
                },
                'headline_observables': {'Ar+': .2 + residual},
                'y': [1. + iteration, 2.],
            },
            'drift_maxima': {'A3': {'argmax': 'Ar excitation',
                                    'absolute_error': residual}},
        }

    output = tmp_path / 'fixed_point.json'
    result = refine_until_converged(failure_record(tmp_path), rerun,
                                    tolerances=tolerances(), generator=generator,
                                    output_path=output)

    assert result['status'] == 'PASS'
    assert result['iteration_count'] == 3
    assert generated[0]['axes']['u'] == pytest.approx([0., .4, .8, .9, 1., 2.])
    assert generated[0]['axes']['x'] == pytest.approx([0., .2, .4, .45, .5, 1.])
    assert [item['fixed_point_version'] for item in generated] == [1, 2, 3]
    assert previous_states == [[1., 2.], [2., 2.], [3., 2.]]
    assert result['iterations'][-1]['agreement'][0]['argmax'] == 'E/N_Td'
    assert result['iterations'][-1]['agreement'][1]['argmax'] == 'mean_energy_eV'
    assert result['iterations'][-1]['agreement'][2]['argmax'] == 'Ar+'
    first_drift = json.loads(output.read_text())['iterations'][0]['drift_maxima']['A3']
    assert first_drift['argmax'] == 'Ar excitation'


def test_fixed_point_f4_limit_is_a_recorded_refusal(tmp_path):
    generated = []

    def generator(spec):
        generated.append(spec)
        return synthetic_artifact(
            tmp_path, 'artifact-' + str(len(generated)))

    def rerun(artifact, previous_y):
        return {
            'status': 'FAIL',
            'terminal_state': {
                'u': .9,
                'composition': {'x': .5},
                'observables': {'mean_energy_eV': 9., 'electron_density': 9.},
                'headline_observables': {'Ar+': 9.},
                'y': previous_y,
            },
            'drift_maxima': {'A3': {'argmax': 'never converges'}},
        }

    output = tmp_path / 'fixed_point.json'
    with pytest.raises(FixedPointRefusal, match='F4 maximum iteration count exceeded'):
        refine_until_converged(failure_record(tmp_path), rerun,
                               tolerances=tolerances(limit=2), generator=generator,
                               output_path=output)

    record = json.loads(output.read_text())
    assert record['status'] == 'REFUSED'
    assert record['iteration_count'] == 2
    assert len(record['iterations']) == 2


def test_composition_refinement_stays_inside_physical_and_axis_bounds():
    from rmgpy.tools.eedf.fixed_point import _extend_composition_axis
    assert _extend_composition_axis([0., .5, 1.], 0.) == pytest.approx(
        [0., .25, .5, 1.])
    assert _extend_composition_axis([0., .5, 1.], 1.) == pytest.approx(
        [0., .5, .75, 1.])
    with pytest.raises(FixedPointRefusal, match='physical bounds'):
        _extend_composition_axis([-.1, .5, 1.], .5)


@pytest.mark.parametrize('decisive_rule', ['F1', 'F2', 'F3'])
def test_each_fixed_point_agreement_rule_is_independently_decisive(
        tmp_path, decisive_rule):
    artifacts = []

    def generator(spec):
        artifact = synthetic_artifact(
            tmp_path, 'artifact-' + str(len(artifacts) + 1))
        artifacts.append(artifact)
        return artifact

    def rerun(artifact, previous_y):
        state = {
            'u': .8,
            'composition': {'x': .4},
            'observables': {'mean_energy_eV': 2., 'electron_density': 3.},
            'headline_observables': {'Ar+': .2},
            'y': previous_y,
        }
        if decisive_rule == 'F1':
            state['u'] = .9
        elif decisive_rule == 'F2':
            state['observables']['mean_energy_eV'] = 3.
        else:
            state['headline_observables']['Ar+'] = .3
        return {'status': 'PASS', 'terminal_state': state, 'drift_maxima': {}}

    output = tmp_path / 'fixed_point.json'
    with pytest.raises(FixedPointRefusal, match='F4 maximum iteration count exceeded'):
        refine_until_converged(
            failure_record(tmp_path), rerun, tolerances=tolerances(limit=1),
            generator=generator, output_path=output)

    record = json.loads(output.read_text())
    checks = {check['rule']: check['passed']
              for check in record['iterations'][0]['agreement']}
    assert checks == {
        'F1': decisive_rule != 'F1',
        'F2': decisive_rule != 'F2',
        'F3': decisive_rule != 'F3',
    }
    assert record['iterations'][0]['artifact_sha256'] is not None
    assert artifacts[0].joinpath('table.h5').stat().st_size > 0
