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

import copy
import json
import os
from pathlib import Path
import runpy
import numpy as np
import pytest

from rmgpy.tools.eedf.generate import generate, model_inputs, qualify_interpolation
from rmgpy.tools.eedf.schema import EEDFError, file_hash, FingerprintMismatch, SpecError
from rmgpy.tools.eedf.validation import held_out_points, run_held_out_accuracy
from rmgpy.solver.eedf import EEDFTable, OutOfDomain

fixture_table = runpy.run_path(
    str(Path(__file__).parents[2] / 'solver' / 'eedfTableTest.py'))['fixture_table']


def test_held_out_midpoints_and_seeded_lhs_never_enter_training_axes():
    axes = {'u': [0., 1., 2.], 'x': [0., 1.]}
    policy = {'cell_midpoints': True, 'lhs_count': 32, 'seed': 314159,
              'off_centre': [{'u': 0.2, 'x': 0.8}]}
    points = held_out_points(axes, policy)
    assert {'u': .5, 'x': .5} in points
    assert {'u': .5, 'x': 0.} in points
    assert {'u': 0., 'x': .5} in points
    assert {'u': .2, 'x': .8} in points
    assert points == held_out_points(axes, policy)
    assert all(not (p['u'] in axes['u'] and p['x'] in axes['x']) for p in points)


def test_mobility_negative_control_passes_runtime_identity_but_fails_held_out(
        tmp_path, monkeypatch):
    path, model = fixture_table(
        tmp_path, composition=True, transport_identity=True)
    table = EEDFTable.load(path, model, artifact_sha256=file_hash(path / 'table.h5'))
    point = {'u': 0.5, 'x': 0.4}
    direct = copy.deepcopy(table.row(point['u'], {'x': point['x']}).as_dict())
    calls = []

    def direct_solver(request):
        calls.append(request)
        return copy.deepcopy(direct)

    # Isolate this falsification to the new T10 mobility/field checks. The
    # existing H1-H4 comparison is exercised by the real integration run.
    monkeypatch.setattr(
        'rmgpy.tools.eedf.validation.compare_row',
        lambda predicted, reference, manifest: [{
            'rule': 'H1-H4-fixture', 'quantity': 'existing comparisons',
            'index': [], 'absolute_error': 0., 'relative_error': 0.,
            'allowed_error': 0., 'passed': True}])
    mobility = table._layout['branch_0']['swarm/mobility_N'][0]
    table._data['branch_0'][..., mobility] *= 1.1

    row = table.row(point['u'], {'x': point['x']})
    criteria = {
        'mobility_N': {'rtol': 2.e-3, 'atol': 0.},
        'field_power': {'rtol': 3.e-3, 'atol': 0.},
    }
    result = run_held_out_accuracy(
        table, direct_solver, [point], criteria=criteria,
        setup_fingerprint=table.manifest['fingerprint'])[0]

    assert row.a6b_runtime['passed']
    assert row.power_groups['field'] == pytest.approx(
        row.swarm['mobility_N'] * (row.EN_Td * 1.e-21) ** 2,
        rel=1.e-14, abs=0.)
    assert not result['passed']
    assert {check['quantity'] for check in result['checks'] if not check['passed']} == {
        'mobility_N', 'field_power'}
    assert calls == [point]


def test_held_out_accuracy_refuses_missing_criteria_setup_mismatch_and_domain(
        tmp_path, monkeypatch):
    path, model = fixture_table(tmp_path, transport_identity=True)
    table = EEDFTable.load(path, model, artifact_sha256=file_hash(path / 'table.h5'))
    solver = lambda point: table.row(point['u'], {}).as_dict()
    criteria = {
        'mobility_N': {'rtol': 2.e-3, 'atol': 0.},
        'field_power': {'rtol': 3.e-3, 'atol': 0.},
    }
    with pytest.raises(SpecError, match='criteria'):
        run_held_out_accuracy(
            table, solver, [{'u': .5}], criteria=None,
            setup_fingerprint=table.manifest['fingerprint'])
    with pytest.raises(FingerprintMismatch, match='setup fingerprint'):
        run_held_out_accuracy(
            table, solver, [{'u': .5}], criteria=criteria,
            setup_fingerprint='wrong')
    with pytest.raises(OutOfDomain):
        run_held_out_accuracy(
            table, solver, [{'u': -0.1}], criteria=criteria,
            setup_fingerprint=table.manifest['fingerprint'])


def test_qualify_interpolation_runs_fresh_solver_at_generated_midpoints(
        tmp_path, monkeypatch):
    path, model = fixture_table(tmp_path, transport_identity=True)
    table = EEDFTable.load(path, model, artifact_sha256=file_hash(path / 'table.h5'))
    channel_map = tmp_path / 'channels.json'
    channel_map.write_text('[]')
    spec = {
        'channel_map': {'path': str(channel_map)},
        'held_out': {
            'cell_midpoints': True, 'lhs_count': 0, 'seed': 7,
            'off_centre': [{'u': .25}],
        },
    }
    calls = []

    monkeypatch.setattr('rmgpy.tools.eedf.generate.validate_spec', lambda candidate: None)
    monkeypatch.setattr('rmgpy.tools.eedf.generate.model_inputs', lambda candidate: model)
    monkeypatch.setattr('rmgpy.tools.eedf.generate._validate_map', lambda candidate: None)
    monkeypatch.setattr('rmgpy.tools.eedf.generate.LoKIDriver', lambda candidate: object())
    monkeypatch.setattr(
        'rmgpy.tools.eedf.validation.compare_row',
        lambda predicted, reference, manifest: [])

    def solve(driver, mapping, candidate, coordinate, fields, label):
        calls.append((dict(coordinate), list(fields), label))
        return [table.row(np.log(fields[0]), coordinate).as_dict()]

    monkeypatch.setattr('rmgpy.tools.eedf.generate._solve', solve)
    criteria = {
        'mobility_N': {'rtol': 1.e-3, 'atol': 0.},
        'field_power': {'rtol': 1.e-3, 'atol': 0.},
    }
    results = qualify_interpolation(table, spec, criteria=criteria)

    assert len(results) == 4
    assert all(result['passed'] for result in results)
    assert [np.log(call[1][0]) for call in calls] == pytest.approx(
        [.25, .5, 1.5, 2.5])
    assert all(call[0] == {} and call[2].startswith('t10_held_') for call in calls)

    missing_off_centre = copy.deepcopy(spec)
    del missing_off_centre['held_out']['off_centre']
    with pytest.raises(SpecError, match='preselected.*off_centre'):
        qualify_interpolation(table, missing_off_centre, criteria=criteria)
    monkeypatch.setattr(
        'rmgpy.tools.eedf.generate.model_inputs', lambda candidate: {'wrong': True})
    with pytest.raises(FingerprintMismatch, match='setup fingerprint'):
        qualify_interpolation(table, spec, criteria=criteria)


def test_qualify_interpolation_refuses_unpinned_multibranch_reference(
        tmp_path, monkeypatch):
    path, model = fixture_table(tmp_path, branches=2)
    table = EEDFTable.load(
        path, model, artifact_sha256=file_hash(path / 'table.h5'))
    spec = {'held_out': {
        'cell_midpoints': True, 'lhs_count': 0, 'seed': 7,
        'off_centre': [{'u': .25}],
    }}
    criteria = {
        'mobility_N': {'rtol': 1.e-3, 'atol': 0.},
        'field_power': {'rtol': 1.e-3, 'atol': 0.},
    }
    monkeypatch.setattr('rmgpy.tools.eedf.generate.validate_spec', lambda candidate: None)
    monkeypatch.setattr('rmgpy.tools.eedf.generate.model_inputs', lambda candidate: model)

    with pytest.raises(EEDFError, match='multi-branch.*pinned direct-solver'):
        qualify_interpolation(table, spec, criteria=criteria)


def test_real_pure_ar_table_passes_legacy_checks_but_awaits_t10(tmp_path):
    spec_path = os.environ.get('EEDF_REAL_SPEC')
    if not spec_path:
        pytest.skip('EEDF_REAL_SPEC supplies the pinned real LoKI integration fixture')
    from rmgpy.tools.eedf.schema import load_spec
    spec = load_spec(spec_path)
    spec['axes']['u'] = np.linspace(np.log(16), np.log(18), 9).tolist()
    spec['held_out']['lhs_count'] = 4
    spec['scratch_root'] = str(tmp_path / 'loki')
    spec['output_root'] = str(tmp_path / 'tables')
    spec['refinement']['max_rounds'] = 0
    path = generate(spec, command_line=['pytest real LoKI fixture'])
    table = EEDFTable.load(
        path, model_inputs(spec), artifact_sha256=file_hash(path / 'table.h5'),
        require_accepted=False)
    assert not table.manifest['accepted']
    assert table.manifest['T10 interpolation qualification'] == {
        'check': 'T10 interpolation qualification',
        'criterion': 'separately frozen mobility and field-power accuracy criteria',
        'verdict_count': 0,
        'passed': False,
        'blocker': 'T10 accuracy criteria have not been supplied',
    }
    assert all(scan['agreement'] for scan in table.manifest['branch_detection'])
    assert len(table.manifest['held_out_verdicts']) >= 12
    assert table.row(np.log(17), {}).swarm['mean_energy_eV'] > 5
    assert table.manifest['row_inputs']['Tg_K'] == 298.15
    assert table.manifest['row_inputs']['P_Pa'] == pytest.approx(5 * 133.322368)
    expected = {'schema_version', 'row_inputs', 'fingerprint', 'axes', 'energy_eV', 'energy_edges_eV',
                'envelopes', 'screen_results', 'tolerances', 'floors', 'channel_map', 'repositories',
                'solver', 'setups', 'state_population_assumptions', 'density_convention', 'units',
                'generator_command', 'generated_at', 'generation_seconds', 'held_out',
                'held_out_verdicts', 'accepted', 'branch_detection', 'limitations',
                'refinement_history', 'branches', 'artifact_sha256', 'branch_certification',
                'A6b-source', 'A6b-runtime', 'T10 interpolation qualification'}
    assert set(table.manifest) == expected
    assert 'cold' in table.manifest['branch_detection'][0]['branch_sources']['branch_0']


@pytest.mark.parametrize('field', ['floors', 'tolerances', 'screen', 'refinement', 'held_out'])
def test_missing_generation_policy_refuses(field):
    from rmgpy.tools.eedf.schema import load_spec, validate_spec
    spec_path = os.environ.get('EEDF_REAL_SPEC')
    if not spec_path:
        pytest.skip('EEDF_REAL_SPEC supplies the complete fixture')
    spec = load_spec(spec_path)
    del spec[field]
    with pytest.raises(SpecError, match='spec fields'):
        validate_spec(spec)
