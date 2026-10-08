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

import rmgpy.tools.eedf.generate as generation
import rmgpy.tools.eedf.validation as validation
from rmgpy.tools.eedf.generate import (
    AxisScreenRefusal, generate, model_inputs, qualify_interpolation,
)
from rmgpy.tools.eedf.schema import EEDFError, file_hash, FingerprintMismatch, SpecError
from rmgpy.tools.eedf.validation import held_out_points
from rmgpy.solver.eedf import EEDFTable

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


def test_solve_preserves_requested_u_over_a_spread_of_valid_axis_values(monkeypatch):
    requested = [0.01, 0.1, 0.2, 0.7, 1.3, 4.0]
    received = []

    class Driver:
        def run(self, coordinates, fields, label):
            received.extend(fields)
            return [{
                'u': float(np.log(field)),
                'EN_Td': float(field),
                'power_groups': {'field': 1.},
            } for field in fields]

    monkeypatch.setattr(generation, 'enrich_row', lambda row, mapping, spec: None)
    rows = generation._solve(
        Driver(), [], {'floors': {'absolute_power_share': 0.}}, {},
        requested, 'canonical-u')

    assert received == pytest.approx(np.exp(requested), rel=0., abs=0.)
    assert [row['u'] for row in rows] == requested


def _screen_fixture(significant):
    spec = {
        'Tg_K': 300.,
        'axes': {'u': [0., 1.], 'x': [0., 1.]},
        'envelopes': {'x': {'reference': .5, 'min': 0., 'max': 1.}},
        'held_out': {
            'off_centre': [{'u': .25, 'x': .2}, {'u': .75, 'x': .8}],
        },
        'screen': {
            'threshold_fraction': .5,
            'candidates': [{'name': 'x', 'reference': .5, 'values': [0., 1.]}],
        },
        'tolerances': {},
        'floors': {},
    }

    def solve(driver, mapping, candidate, coordinates, values, label):
        return [{'u': float(i)} for i, unused in enumerate(values)]

    def compare(predicted, reference, manifest):
        return [{
            'rule': 'H1', 'quantity': 'fixture',
            'absolute_error': 1. if significant else 0.,
            'allowed_error': 1., 'passed': not significant,
        }]

    return spec, solve, compare


@pytest.mark.parametrize('significant', [False, True])
def test_screen_reconciles_off_centre_points_only_when_an_axis_is_removed(
        monkeypatch, significant):
    spec, solve, compare = _screen_fixture(significant)
    original = copy.deepcopy(spec['held_out']['off_centre'])
    monkeypatch.setattr(generation, '_solve', solve)
    monkeypatch.setattr(generation, 'compare_row', compare)

    generation._screen(spec, object(), [], 'screen')

    if significant:
        assert spec['axes']['x'] == [0., 1.]
        assert spec['held_out']['off_centre'] == original
    else:
        assert 'x' not in spec['axes']
        assert spec['held_out']['off_centre'] == [{'u': .25}, {'u': .75}]


def test_screen_axis_removal_is_a_noop_when_no_off_centre_points_are_declared(
        monkeypatch):
    spec, solve, compare = _screen_fixture(False)
    del spec['held_out']['off_centre']
    monkeypatch.setattr(generation, '_solve', solve)
    monkeypatch.setattr(generation, 'compare_row', compare)

    generation._screen(spec, object(), [], 'screen')

    assert 'x' not in spec['axes']
    assert 'off_centre' not in spec['held_out']


@pytest.mark.parametrize(
    'points, message', [
        ([{'u': .25, 'x': .2}, {'u': .25, 'x': .8}], 'collide'),
        ([{'u': 0., 'x': .2}], 'training grid'),
    ],
)
def test_screen_refuses_off_centre_points_that_lose_meaning_after_projection(
        monkeypatch, points, message):
    spec, solve, compare = _screen_fixture(False)
    spec['held_out']['off_centre'] = points
    original = copy.deepcopy(spec)
    monkeypatch.setattr(generation, '_solve', solve)
    monkeypatch.setattr(generation, 'compare_row', compare)

    with pytest.raises(AxisScreenRefusal, match='screening x.*' + message):
        generation._screen(spec, object(), [], 'screen')

    assert spec == original


def test_qualification_exposes_no_callback_comparator():
    assert not hasattr(validation, 'run_held_out_accuracy')


def test_mobility_negative_control_passes_runtime_identity_but_fails_held_out(
        tmp_path, monkeypatch):
    spec = {
        'channel_map': {'path': str(tmp_path / 'channels.json')},
        'held_out': {
            'cell_midpoints': True, 'lhs_count': 0, 'seed': 7,
            'off_centre': [{'u': .25}],
        },
    }
    retained = held_out_points({'u': [0., 1., 2., 3.]}, spec['held_out'])
    path, model = fixture_table(
        tmp_path / 'candidate', transport_identity=True,
        retained_held_out=retained)
    reference_path, _ = fixture_table(
        tmp_path / 'reference', transport_identity=True)
    table = EEDFTable.load(path, model, artifact_sha256=file_hash(path / 'table.h5'))
    reference_table = EEDFTable.load(
        reference_path, model, artifact_sha256=file_hash(reference_path / 'table.h5'))
    Path(spec['channel_map']['path']).write_text('[]')
    calls = []

    monkeypatch.setattr('rmgpy.tools.eedf.generate.validate_spec', lambda candidate: None)
    monkeypatch.setattr('rmgpy.tools.eedf.generate.model_inputs', lambda candidate: model)
    monkeypatch.setattr('rmgpy.tools.eedf.generate._validate_map', lambda candidate: None)
    monkeypatch.setattr('rmgpy.tools.eedf.generate.LoKIDriver', lambda candidate: object())

    def solve(driver, mapping, candidate, coordinate, u_values, label):
        calls.append((dict(coordinate), list(u_values), label))
        return [reference_table.row(u_values[0], coordinate).as_dict()]

    monkeypatch.setattr('rmgpy.tools.eedf.generate._solve', solve)

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

    criteria = {
        'mobility_N': {'rtol': 2.e-3, 'atol': 0.},
        'field_power': {'rtol': 3.e-3, 'atol': 0.},
    }
    results = qualify_interpolation(table, spec, criteria=criteria)
    row = table.row(.5, {})

    assert row.a6b_runtime['passed']
    assert row.power_groups['field'] == pytest.approx(
        row.swarm['mobility_N'] * (row.EN_Td * 1.e-21) ** 2,
        rel=1.e-14, abs=0.)
    assert results and all(not result['passed'] for result in results)
    assert all(
        {check['quantity'] for check in result['checks'] if not check['passed']} == {
            'mobility_N', 'field_power'}
        for result in results)
    assert all(result['reference'] == 'fresh direct solver output' for result in results)
    assert len(calls) == len(retained)


def test_qualify_interpolation_refuses_an_empty_pinned_point_set(
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
    criteria = {
        'mobility_N': {'rtol': 2.e-3, 'atol': 0.},
        'field_power': {'rtol': 3.e-3, 'atol': 0.},
    }
    table.manifest['held_out'] = []
    monkeypatch.setattr('rmgpy.tools.eedf.generate.validate_spec', lambda candidate: None)
    monkeypatch.setattr('rmgpy.tools.eedf.generate.model_inputs', lambda candidate: model)
    monkeypatch.setattr('rmgpy.tools.eedf.generate._validate_map', lambda candidate: None)
    monkeypatch.setattr('rmgpy.tools.eedf.generate.LoKIDriver', lambda candidate: object())

    with pytest.raises(SpecError, match='nonempty pinned held-out point set'):
        qualify_interpolation(table, spec, criteria=criteria)


def test_qualify_interpolation_refuses_declared_off_centre_training_node_before_solver(
        tmp_path, monkeypatch):
    path, model = fixture_table(tmp_path, transport_identity=True)
    table = EEDFTable.load(path, model, artifact_sha256=file_hash(path / 'table.h5'))
    channel_map = tmp_path / 'channels.json'
    channel_map.write_text('[]')
    spec = {
        'channel_map': {'path': str(channel_map)},
        'held_out': {
            'cell_midpoints': True, 'lhs_count': 0, 'seed': 7,
            'off_centre': [
                {'u': .25},
                {'u': table.manifest['axes']['u'][0]},
            ],
        },
    }
    criteria = {
        'mobility_N': {'rtol': 1.e-3, 'atol': 0.},
        'field_power': {'rtol': 1.e-3, 'atol': 0.},
    }
    monkeypatch.setattr('rmgpy.tools.eedf.generate.validate_spec', lambda candidate: None)
    monkeypatch.setattr('rmgpy.tools.eedf.generate.model_inputs', lambda candidate: model)
    monkeypatch.setattr('rmgpy.tools.eedf.generate._validate_map', lambda candidate: None)
    monkeypatch.setattr(
        'rmgpy.tools.eedf.generate.LoKIDriver',
        lambda candidate: pytest.fail('qualification solver was constructed'))

    with pytest.raises(
            SpecError, match='declared held_out off_centre point.*training grid'):
        qualify_interpolation(table, spec, criteria=criteria)


def test_qualify_interpolation_refuses_a_partial_evaluated_point_set(
        tmp_path, monkeypatch):
    path, model = fixture_table(
        tmp_path, transport_identity=True,
        retained_held_out=[{'u': .25}, {'u': .5}])
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
    monkeypatch.setattr('rmgpy.tools.eedf.generate.validate_spec', lambda candidate: None)
    monkeypatch.setattr('rmgpy.tools.eedf.generate.model_inputs', lambda candidate: model)
    monkeypatch.setattr('rmgpy.tools.eedf.generate._validate_map', lambda candidate: None)
    monkeypatch.setattr('rmgpy.tools.eedf.generate.LoKIDriver', lambda candidate: object())
    monkeypatch.setattr(
        'rmgpy.tools.eedf.generate._solve',
        lambda driver, mapping, candidate, coordinate, fields, label: [{}])
    monkeypatch.setattr(
        'rmgpy.tools.eedf.generate._held_out_accuracy',
        lambda table, direct, point, criteria, setup_fingerprint: {
            'point': {'u': -1.}, 'passed': True,
        })

    with pytest.raises(EEDFError, match='partial pinned held-out point set'):
        qualify_interpolation(table, spec, criteria={
            'mobility_N': {'rtol': 1.e-3, 'atol': 0.},
            'field_power': {'rtol': 1.e-3, 'atol': 0.},
        })


def test_qualify_interpolation_runs_fresh_solver_at_persisted_held_out_points(
        tmp_path, monkeypatch):
    persisted = [
        {'u': .25},
        {'u': .5},
        {'u': float(np.nextafter(.5, np.inf))},
        {'u': 1.5},
        {'u': 2.5},
    ]
    path, model = fixture_table(
        tmp_path, transport_identity=True, retained_held_out=persisted)
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

    def solve(driver, mapping, candidate, coordinate, u_values, label):
        calls.append((dict(coordinate), list(u_values), label))
        return [table.row(u_values[0], coordinate).as_dict()]

    monkeypatch.setattr('rmgpy.tools.eedf.generate._solve', solve)
    criteria = {
        'mobility_N': {'rtol': 1.e-3, 'atol': 0.},
        'field_power': {'rtol': 1.e-3, 'atol': 0.},
    }
    results = qualify_interpolation(table, spec, criteria=criteria)

    assert len(results) == len(persisted)
    assert all(result['passed'] for result in results)
    assert [result['point'] for result in results] == persisted
    assert [call[1][0] for call in calls] == [
        point['u'] for point in persisted]
    assert all(call[0] == {} and call[2].startswith('t10_') for call in calls)

    first_labels = {call[2] for call in calls}
    calls.clear()
    repeated = qualify_interpolation(table, spec, criteria=criteria)
    second_labels = {call[2] for call in calls}
    assert repeated == results
    assert first_labels.isdisjoint(second_labels)

    missing_off_centre = copy.deepcopy(spec)
    del missing_off_centre['held_out']['off_centre']
    with pytest.raises(SpecError, match='preselected.*off_centre'):
        qualify_interpolation(table, missing_off_centre, criteria=criteria)
    monkeypatch.setattr(
        'rmgpy.tools.eedf.generate.model_inputs', lambda candidate: {'wrong': True})
    with pytest.raises(FingerprintMismatch, match='setup fingerprint'):
        qualify_interpolation(table, spec, criteria=criteria)


def test_public_manifest_copy_cannot_reduce_loaded_qualification_points(
        tmp_path, monkeypatch):
    retained = [{'u': .25}, {'u': .5}]
    path, model = fixture_table(
        tmp_path, transport_identity=True, retained_held_out=retained)
    table = EEDFTable.load(
        path, model, artifact_sha256=file_hash(path / 'table.h5'))
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

    monkeypatch.setattr(
        'rmgpy.tools.eedf.generate.validate_spec', lambda candidate: None)
    monkeypatch.setattr(
        'rmgpy.tools.eedf.generate.model_inputs', lambda candidate: model)
    monkeypatch.setattr(
        'rmgpy.tools.eedf.generate._validate_map', lambda candidate: None)
    monkeypatch.setattr(
        'rmgpy.tools.eedf.generate.LoKIDriver', lambda candidate: object())
    monkeypatch.setattr(
        'rmgpy.tools.eedf.validation.compare_row',
        lambda predicted, reference, manifest: [])

    def solve(driver, mapping, candidate, coordinate, u_values, label):
        calls.append({'u': u_values[0], **coordinate})
        return [table.row(u_values[0], coordinate).as_dict()]

    monkeypatch.setattr('rmgpy.tools.eedf.generate._solve', solve)
    preview = table.manifest.copy()
    preview['held_out'].pop()

    results = qualify_interpolation(table, spec, criteria={
        'mobility_N': {'rtol': 1.e-3, 'atol': 0.},
        'field_power': {'rtol': 1.e-3, 'atol': 0.},
    })

    assert preview['held_out'] == retained[:-1]
    assert table.manifest['held_out'] == retained
    assert calls == retained
    assert [result['point'] for result in results] == retained
    assert all(result['passed'] for result in results)


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


def test_real_pure_ar_table_passes_legacy_checks_but_awaits_t10(
        tmp_path, monkeypatch):
    spec_path = os.environ.get('EEDF_REAL_SPEC')
    if not spec_path:
        pytest.skip('EEDF_REAL_SPEC supplies the pinned real LoKI integration fixture')
    from rmgpy.tools.eedf.schema import load_spec
    spec = load_spec(spec_path)
    spec['axes']['u'] = np.linspace(np.log(16), np.log(18), 9).tolist()
    spec['held_out']['lhs_count'] = 4
    spec['held_out']['off_centre'] = [{'u': spec['axes']['u'][0]}]
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
    monkeypatch.setattr(
        'rmgpy.tools.eedf.generate.LoKIDriver',
        lambda candidate: pytest.fail('qualification solver was constructed'))
    with pytest.raises(
            SpecError, match='declared held_out off_centre point.*training grid'):
        qualify_interpolation(table, spec, criteria={
            'mobility_N': {'rtol': 1.e-3, 'atol': 0.},
            'field_power': {'rtol': 1.e-3, 'atol': 0.},
        })
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
