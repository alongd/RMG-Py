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

import concurrent.futures
import builtins
import copy
import hashlib
import json
import os
from pathlib import Path
import re
import runpy
import threading
import time
from types import MappingProxyType
import numpy as np
import pytest

import rmgpy.tools.eedf.generate as generation
import rmgpy.tools.eedf.loki as loki
import rmgpy.tools.eedf.validation as validation
from rmgpy.tools.eedf.generate import (
    AxisScreenRefusal, generate, model_inputs, qualify_interpolation,
)
from rmgpy.tools.eedf.schema import EEDFError, _file_hash, FingerprintMismatch, SpecError
from rmgpy.tools.eedf.validation import _held_out_points
from rmgpy.solver.eedf import EEDFTable

fixture_table = runpy.run_path(
    str(Path(__file__).parents[2] / 'solver' / 'eedfTableTest.py'))['fixture_table']


def publication_fixture(tmp_path, monkeypatch, *, failure=None):
    """Exercise real storage, branch detection and T10 with analytic transport.

    Stub external setup and legacy comparisons, whose individual physical
    criteria have separate tests. Inject their failures at the publication
    boundary; keep the T10 comparator and production loader real.
    """
    from rmgpy.tools.eedf.schema import _content_hash

    path, model = fixture_table(tmp_path / 'seed')
    seed = EEDFTable.load(path, model, artifact_sha256=_file_hash(path / 'table.h5'))
    manifest = seed.manifest
    channels = manifest['channel_map']
    if failure == 'classification':
        channels[0]['classification'] = 'C'
    channel_path = tmp_path / 'channels.json'
    channel_path.write_text(json.dumps(channels))
    spec = {
        'axes': {'u': [0., 1., 2., 3.]},
        'held_out': {'cell_midpoints': True, 'lhs_count': 0, 'seed': 7,
                     'off_centre': [{'u': .25}]},
        'interpolation_accuracy': {name: {'rtol': .01, 'atol': 0.}
                                   for name in ('mobility_N', 'field_power')},
        'channel_map': {'path': str(channel_path), 'sha256': _file_hash(channel_path)},
        'tolerances': {'G1': {'rtol': 1.e-6, 'atol': 0.}},
        'refinement': {'max_rounds': 0}, 'max_parallel': 1,
        'scratch_root': str(tmp_path / 'scratch'),
        'output_root': str(tmp_path / 'tables'),
    }
    calls = []

    def inputs(candidate, mapping):
        return dict(model, axes=copy.deepcopy(candidate['axes']))

    def make_manifest(candidate, driver, mapping, screen, scans, row, command):
        result = copy.deepcopy(manifest)
        result.update(row_inputs=inputs(candidate, mapping), axes=candidate['axes'],
                      fingerprint=_content_hash(inputs(candidate, mapping)),
                      branch_detection=scans, limitations=[], generation_seconds=0.,
                      branch_certification={branch: 'uncertified: unseeded scans'
                                            for branch in scans[0]['branch_sources']})
        return result

    def solve(driver, mapping, candidate, coordinate, values, label):
        calls.append((label, list(values)))
        rows = []
        for u in values:
            row = seed.row(0., {})._as_dict()
            row.update(u=u, EN_Td=float(np.exp(u)), composition={}, setup=label,
                       converged=True, iteration_count=-1)
            row['swarm']['mean_energy_eV'] = 1. + u
            row['swarm']['mobility_N'] = (1. + .15 * u) * 1.e24
            if failure == 't10' and label.startswith('t10_'):
                row['swarm']['mobility_N'] *= 1.1
            if failure == 'branch' and label.endswith('_down'):
                row['swarm']['mobility_N'] *= 1.1
            row['power_groups']['field'] = row['swarm']['mobility_N'] * (row['EN_Td'] * 1.e-21) ** 2
            if failure == 'source' and '_held' not in label:
                row['power_groups']['field'] *= 1.01
            rows.append(row)
        return rows

    def compare(predicted, reference, candidate):
        passed = not (failure == 'legacy' and not reference['setup'].startswith('t10_'))
        return [{'rule': 'H1', 'quantity': 'fixture', 'index': [],
                 'absolute_error': 0. if passed else 1., 'relative_error': 0. if passed else 1.,
                 'allowed_error': 0., 'passed': passed}]

    monkeypatch.setattr(generation, '_validate_spec', lambda candidate: None)
    monkeypatch.setattr(generation, '_validate_map', lambda mapping: None)
    monkeypatch.setattr(generation, '_screen', lambda *args: [])
    monkeypatch.setattr(generation, 'LoKIDriver', lambda candidate: object())
    monkeypatch.setattr(generation, '_model_inputs', inputs)
    monkeypatch.setattr(generation, '_manifest', make_manifest)
    monkeypatch.setattr(generation, '_solve', solve)
    monkeypatch.setattr(generation, '_compare_row', compare)
    monkeypatch.setattr(validation, '_compare_row', compare)
    return spec, calls, inputs


@pytest.mark.parametrize('failure', [None, 'legacy', 't10', 'source', 'classification'])
def test_generate_publishes_only_when_all_qualification_gates_pass(
        tmp_path, monkeypatch, failure):
    spec, calls, inputs = publication_fixture(tmp_path, monkeypatch, failure=failure)
    path = generate(spec, command_line=['pytest publication'])
    manifest = json.loads((path / 'manifest.json').read_text())
    assert manifest['accepted'] is (failure is None)
    t10 = [v for v in manifest['held_out_verdicts'] if 'criteria' in v]
    legacy = [v for v in manifest['held_out_verdicts'] if 'criteria' not in v]
    assert len(t10) == len(legacy) == len(manifest['held_out']) == 4
    assert [v['point'] for v in t10] == manifest['held_out']
    assert all(v['criteria'] == spec['interpolation_accuracy'] for v in t10)
    assert all(v['reference'] == 'fresh direct solver output' for v in t10)
    assert len([label for label, _ in calls if label.startswith('t10_')]) == 4
    assert manifest['T10 interpolation qualification']['passed'] is (failure != 't10')
    assert '| T10 | 8 |' in (path / 'validation.md').read_text()
    assert json.loads((path / 'generation_spec.json').read_text()) == spec
    if failure is None:
        loaded = EEDFTable.load(path, inputs(spec, []), artifact_sha256=manifest['artifact_sha256'])
        assert loaded.manifest['accepted']
        # The fresh results are bound to the HDF5 qualification digest.
        manifest['held_out_verdicts'][-1]['criteria']['mobility_N']['rtol'] = .5
        (path / 'manifest.json').write_text(json.dumps(manifest))
        with pytest.raises(FingerprintMismatch, match='qualification'):
            EEDFTable.load(path, inputs(spec, []), artifact_sha256=manifest['artifact_sha256'])
    else:
        with pytest.raises(FingerprintMismatch, match='qualification'):
            EEDFTable.load(path, inputs(spec, []), artifact_sha256=manifest['artifact_sha256'])


def test_generate_refines_t10_failures_without_promoting_retained_points(tmp_path, monkeypatch):
    spec, calls, inputs = publication_fixture(tmp_path, monkeypatch, failure='t10')
    spec['refinement']['max_rounds'] = 1
    path = generate(spec)
    manifest = json.loads((path / 'manifest.json').read_text())
    history = manifest['refinement_history']
    assert len(history) == 2
    assert history[0]['failing_points'] == [{'u': .25}, {'u': .5}, {'u': 1.5}, {'u': 2.5}]
    assert len(history[1]['axes']['u']) > len(history[0]['axes']['u'])
    assert all(point['u'] not in spec['axes']['u'] for point in manifest['held_out'])
    assert len(manifest['held_out']) > 4
    assert not manifest['accepted']


def test_generate_without_frozen_criteria_stays_unaccepted(tmp_path, monkeypatch):
    spec, calls, inputs = publication_fixture(tmp_path, monkeypatch)
    del spec['interpolation_accuracy']
    path = generate(spec)
    manifest = json.loads((path / 'manifest.json').read_text())
    assert not manifest['accepted']
    assert not any(label.startswith('t10_') for label, _ in calls)
    assert manifest['T10 interpolation qualification']['verdict_count'] == 0


@pytest.mark.parametrize('failure', ['branch', 'off_centre', 'fingerprint'])
def test_generate_refuses_unqualified_branch_or_point_identity(tmp_path, monkeypatch, failure):
    spec, calls, inputs = publication_fixture(tmp_path, monkeypatch, failure=failure)
    message = 'multi-branch'
    if failure == 'off_centre':
        spec['held_out']['off_centre'] = [{'u': 0.}]
        message = 'training grid'
    elif failure == 'fingerprint':
        original = generation._manifest

        def changed_inputs(*args):
            result = original(*args)
            monkeypatch.setattr(generation, '_model_inputs', lambda *args: {'changed': True})
            return result

        monkeypatch.setattr(generation, '_manifest', changed_inputs)
        message = 'fingerprint'
    with pytest.raises(EEDFError, match=message):
        generate(spec)
    assert not Path(spec['output_root']).exists()


def test_channel_map_replacement_between_hash_and_parse(tmp_path, monkeypatch):
    mapping = tmp_path / 'channels.json'
    original = [{'reaction': None}]
    mapping.write_text(json.dumps(original))
    pin = _file_hash(mapping)
    other = tmp_path / 'other'
    other.write_text('other input')
    other_pin = {'path': str(other), 'sha256': _file_hash(other)}
    spec = {'channel_map': {'path': str(mapping), 'sha256': pin},
            'binary': other_pin, 'cmake_cache': other_pin,
            'input_files': {}, 'axes': {'u': [0., 1.]}}
    for name, replacement in {
            '_check_solver': lambda spec: {},
            '_physical_properties': lambda spec: ({}, {}),
            '_solver_environment': lambda spec: {},
            '_row_inputs': lambda spec: {},
            '_interpolant_identity': lambda: {},
            'qualification_setup_sha256': lambda spec: 'setup',
            '_validate_physical_map': lambda *args: None}.items():
        monkeypatch.setattr(generation, name, replacement)
    real_sha256 = hashlib.sha256
    replacements = []

    class ReplacingDigest:
        def __init__(self, *args, **kwargs):
            self.digest = real_sha256(*args, **kwargs)

        def update(self, value):
            self.digest.update(value)

        def hexdigest(self):
            digest = self.digest.hexdigest()
            if digest == pin and not replacements:
                mapping.write_text('[{"reaction": {"index": 999}}]')
                replacements.append(True)
            return digest

    monkeypatch.setattr(hashlib, 'sha256', ReplacingDigest)
    result = model_inputs(spec)
    assert replacements == [True], 'the replacement must actually occur'
    assert result['channel_map_sha256'] == pin
    assert result['reactions'] == []


@pytest.mark.parametrize('problem', ['missing', 'malformed', 'changed'])
def test_pinned_channel_map_refuses_invalid_input(tmp_path, problem):
    from rmgpy.tools.eedf.channels import _load_channel_map
    path = tmp_path / 'channels.json'
    path.write_text('broken JSON' if problem == 'malformed' else '[]')
    pin = {'path': str(path), 'sha256': _file_hash(path)}
    if problem == 'missing':
        path.unlink()
    elif problem == 'changed':
        path.write_text('[{}]')
    with pytest.raises((SpecError, FingerprintMismatch)):
        _load_channel_map(pin)


def test_generate_keeps_verified_map_when_driver_replaces_path(tmp_path, monkeypatch):
    path = tmp_path / 'channels.json'
    path.write_text('[{"reaction": null}]')
    spec = {'channel_map': {'path': str(path), 'sha256': _file_hash(path)},
            'axes': {'u': [0., 1.]}}
    seen = []

    class StopAfterSnapshot(Exception):
        pass

    def driver(candidate):
        path.write_text('[{"reaction": {"index": 999}}]')
        return object()

    def screen(candidate, driver, channels, prefix):
        seen.extend(channels)
        raise StopAfterSnapshot

    monkeypatch.setattr(generation, '_validate_spec', lambda spec: None)
    monkeypatch.setattr(generation, '_validate_map', lambda channels: None)
    monkeypatch.setattr(generation, 'LoKIDriver', driver)
    monkeypatch.setattr(generation, '_screen', screen)
    with pytest.raises(StopAfterSnapshot):
        generate(spec)
    assert seen == [{'reaction': None}]
    assert json.loads(path.read_text()) == [{'reaction': {'index': 999}}]


def test_generate_refuses_readonly_spec_before_import(monkeypatch):
    original_import = builtins.__import__

    def check_import(name, globals=None, locals=None, fromlist=(), level=0):
        assert name != 'rmgpy.solver.eedf', 'import preceded input refusal'
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, '__import__', check_import)
    with pytest.raises(TypeError, match='mutable'):
        generate(MappingProxyType({}), command_line=['pytest readonly spec'])


def test_generate_owns_inputs_before_waiting_for_import(monkeypatch):
    caller = {
        'Tg_K': 298.15,
        'envelopes': {
            'Tg_K': {'min': 298.15, 'max': 298.15, 'reference': 298.15},
        },
    }
    entered_import = threading.Event()
    release_import = threading.Event()
    observed = {}
    original_import = builtins.__import__

    class StopGeneration(Exception):
        pass

    def waiting_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == 'rmgpy.solver.eedf' and 'EEDFTable' in fromlist:
            entered_import.set()
            assert release_import.wait(5.)
        return original_import(name, globals, locals, fromlist, level)

    def record_validation(spec):
        observed['reference'] = spec['envelopes']['Tg_K']['reference']
        raise StopGeneration

    monkeypatch.setattr(builtins, '__import__', waiting_import)
    monkeypatch.setattr(generation, '_validate_spec', record_validation)
    errors = []

    def invoke():
        try:
            generate(caller, command_line=['ownership probe'])
        except Exception as error:
            errors.append(error)

    worker = threading.Thread(target=invoke)
    worker.start()
    assert entered_import.wait(5.)
    caller.copy()['envelopes']['Tg_K']['reference'] = 300.
    release_import.set()
    worker.join(5.)

    assert not worker.is_alive()
    assert len(errors) == 1 and isinstance(errors[0], StopGeneration)
    assert observed['reference'] == 298.15


def test_held_out_midpoints_and_seeded_lhs_never_enter_training_axes():
    axes = {'u': [0., 1., 2.], 'x': [0., 1.]}
    policy = {'cell_midpoints': True, 'lhs_count': 32, 'seed': 314159,
              'off_centre': [{'u': 0.2, 'x': 0.8}]}
    points = _held_out_points(axes, policy)
    assert {'u': .5, 'x': .5} in points
    assert {'u': .5, 'x': 0.} in points
    assert {'u': 0., 'x': .5} in points
    assert {'u': .2, 'x': .8} in points
    assert points == _held_out_points(axes, policy)
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

    monkeypatch.setattr(generation, '_enrich_row', lambda row, mapping, spec: None)
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
    monkeypatch.setattr(generation, '_compare_row', compare)

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
    monkeypatch.setattr(generation, '_compare_row', compare)

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
    monkeypatch.setattr(generation, '_compare_row', compare)

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
    retained = _held_out_points({'u': [0., 1., 2., 3.]}, spec['held_out'])
    path, model = fixture_table(
        tmp_path / 'candidate', transport_identity=True,
        retained_held_out=retained)
    reference_path, _ = fixture_table(
        tmp_path / 'reference', transport_identity=True)
    table = EEDFTable.load(path, model, artifact_sha256=_file_hash(path / 'table.h5'))
    reference_table = EEDFTable.load(
        reference_path, model, artifact_sha256=_file_hash(reference_path / 'table.h5'))
    Path(spec['channel_map']['path']).write_text('[]')
    spec['channel_map']['sha256'] = _file_hash(spec['channel_map']['path'])
    calls = []

    monkeypatch.setattr('rmgpy.tools.eedf.generate._validate_spec', lambda candidate: None)
    monkeypatch.setattr('rmgpy.tools.eedf.generate._model_inputs', lambda candidate, mapping: model)
    monkeypatch.setattr('rmgpy.tools.eedf.generate._validate_map', lambda candidate: None)
    monkeypatch.setattr('rmgpy.tools.eedf.generate.LoKIDriver', lambda candidate: object())

    def solve(driver, mapping, candidate, coordinate, u_values, label):
        calls.append((dict(coordinate), list(u_values), label))
        return [reference_table.row(u_values[0], coordinate)._as_dict()]

    monkeypatch.setattr('rmgpy.tools.eedf.generate._solve', solve)

    # Isolate this falsification to the new T10 mobility/field checks. The
    # existing H1-H4 comparison is exercised by the real integration run.
    monkeypatch.setattr(
        'rmgpy.tools.eedf.validation._compare_row',
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
    table = EEDFTable.load(path, model, artifact_sha256=_file_hash(path / 'table.h5'))
    channel_map = tmp_path / 'channels.json'
    channel_map.write_text('[]')
    spec = {
        'channel_map': {'path': str(channel_map), 'sha256': _file_hash(channel_map)},
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
    monkeypatch.setattr('rmgpy.tools.eedf.generate._validate_spec', lambda candidate: None)
    monkeypatch.setattr('rmgpy.tools.eedf.generate._model_inputs', lambda candidate, mapping: model)
    monkeypatch.setattr('rmgpy.tools.eedf.generate._validate_map', lambda candidate: None)
    monkeypatch.setattr('rmgpy.tools.eedf.generate.LoKIDriver', lambda candidate: object())

    with pytest.raises(SpecError, match='nonempty pinned held-out point set'):
        qualify_interpolation(table, spec, criteria=criteria)


def test_qualify_interpolation_refuses_declared_off_centre_training_node_before_solver(
        tmp_path, monkeypatch):
    path, model = fixture_table(tmp_path, transport_identity=True)
    table = EEDFTable.load(path, model, artifact_sha256=_file_hash(path / 'table.h5'))
    channel_map = tmp_path / 'channels.json'
    channel_map.write_text('[]')
    spec = {
        'channel_map': {'path': str(channel_map), 'sha256': _file_hash(channel_map)},
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
    monkeypatch.setattr('rmgpy.tools.eedf.generate._validate_spec', lambda candidate: None)
    monkeypatch.setattr('rmgpy.tools.eedf.generate._model_inputs', lambda candidate, mapping: model)
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
    table = EEDFTable.load(path, model, artifact_sha256=_file_hash(path / 'table.h5'))
    channel_map = tmp_path / 'channels.json'
    channel_map.write_text('[]')
    spec = {
        'channel_map': {'path': str(channel_map), 'sha256': _file_hash(channel_map)},
        'held_out': {
            'cell_midpoints': True, 'lhs_count': 0, 'seed': 7,
            'off_centre': [{'u': .25}],
        },
    }
    monkeypatch.setattr('rmgpy.tools.eedf.generate._validate_spec', lambda candidate: None)
    monkeypatch.setattr('rmgpy.tools.eedf.generate._model_inputs', lambda candidate, mapping: model)
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


@pytest.mark.parametrize('frozen_in_spec', [False, True])
def test_qualify_interpolation_runs_fresh_solver_at_persisted_held_out_points(
        tmp_path, monkeypatch, frozen_in_spec):
    persisted = [
        {'u': .25},
        {'u': .5},
        {'u': float(np.nextafter(.5, np.inf))},
        {'u': 1.5},
        {'u': 2.5},
    ]
    path, model = fixture_table(
        tmp_path, transport_identity=True, retained_held_out=persisted)
    table = EEDFTable.load(path, model, artifact_sha256=_file_hash(path / 'table.h5'))
    channel_map = tmp_path / 'channels.json'
    channel_map.write_text('[]')
    spec = {
        'channel_map': {'path': str(channel_map), 'sha256': _file_hash(channel_map)},
        'held_out': {
            'cell_midpoints': True, 'lhs_count': 0, 'seed': 7,
            'off_centre': [{'u': .25}],
        },
    }
    calls = []

    monkeypatch.setattr('rmgpy.tools.eedf.generate._validate_spec', lambda candidate: None)
    monkeypatch.setattr('rmgpy.tools.eedf.generate._model_inputs', lambda candidate, mapping: model)
    monkeypatch.setattr('rmgpy.tools.eedf.generate._validate_map', lambda candidate: None)
    monkeypatch.setattr('rmgpy.tools.eedf.generate.LoKIDriver', lambda candidate: object())
    monkeypatch.setattr(
        'rmgpy.tools.eedf.validation._compare_row',
        lambda predicted, reference, manifest: [])

    def solve(driver, mapping, candidate, coordinate, u_values, label):
        calls.append((dict(coordinate), list(u_values), label))
        return [table.row(u_values[0], coordinate)._as_dict()]

    monkeypatch.setattr('rmgpy.tools.eedf.generate._solve', solve)
    criteria = {
        'mobility_N': {'rtol': 1.e-3, 'atol': 0.},
        'field_power': {'rtol': 1.e-3, 'atol': 0.},
    }
    if frozen_in_spec:
        spec['interpolation_accuracy'] = criteria
        with pytest.raises(SpecError, match='frozen interpolation_accuracy'):
            qualify_interpolation(table, spec, criteria={
                key: {'rtol': .1, 'atol': 0.} for key in criteria})
    results = qualify_interpolation(table, spec, criteria=None if frozen_in_spec else criteria)

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
        'rmgpy.tools.eedf.generate._model_inputs', lambda candidate, mapping: {'wrong': True})
    with pytest.raises(FingerprintMismatch, match='setup fingerprint'):
        qualify_interpolation(table, spec, criteria=criteria)


def test_public_manifest_copy_cannot_reduce_loaded_qualification_points(
        tmp_path, monkeypatch):
    retained = [{'u': .25}, {'u': .5}]
    path, model = fixture_table(
        tmp_path, transport_identity=True, retained_held_out=retained)
    table = EEDFTable.load(
        path, model, artifact_sha256=_file_hash(path / 'table.h5'))
    channel_map = tmp_path / 'channels.json'
    channel_map.write_text('[]')
    spec = {
        'channel_map': {'path': str(channel_map), 'sha256': _file_hash(channel_map)},
        'held_out': {
            'cell_midpoints': True, 'lhs_count': 0, 'seed': 7,
            'off_centre': [{'u': .25}],
        },
    }
    calls = []

    monkeypatch.setattr(
        'rmgpy.tools.eedf.generate._validate_spec', lambda candidate: None)
    monkeypatch.setattr(
        'rmgpy.tools.eedf.generate._model_inputs', lambda candidate, mapping: model)
    monkeypatch.setattr(
        'rmgpy.tools.eedf.generate._validate_map', lambda candidate: None)
    monkeypatch.setattr(
        'rmgpy.tools.eedf.generate.LoKIDriver', lambda candidate: object())
    monkeypatch.setattr(
        'rmgpy.tools.eedf.validation._compare_row',
        lambda predicted, reference, manifest: [])

    def solve(driver, mapping, candidate, coordinate, u_values, label):
        calls.append({'u': u_values[0], **coordinate})
        return [table.row(u_values[0], coordinate)._as_dict()]

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


@pytest.mark.parametrize(
    'field, pinned, changed', [
        ('Tg_K', 298.15, 300.),
        ('P_Pa', 666.61184, 700.),
    ],
)
def test_qualification_owns_caller_spec_during_every_fresh_solve(
        tmp_path, monkeypatch, field, pinned, changed):
    retained = [{'u': .25}, {'u': .5}]
    path, model = fixture_table(
        tmp_path, transport_identity=True, retained_held_out=retained)
    table = EEDFTable.load(
        path, model, artifact_sha256=_file_hash(path / 'table.h5'))
    channel_map = tmp_path / 'channels.json'
    channel_map.write_text('[]')
    spec = {
        'Tg_K': 298.15,
        'P_Pa': 666.61184,
        'envelopes': {
            'Tg_K': {'min': 298.15, 'max': 298.15, 'reference': 298.15},
            'P_Pa': {
                'min': 666.61184, 'max': 666.61184,
                'reference': 666.61184,
            },
        },
        'channel_map': {'path': str(channel_map), 'sha256': _file_hash(channel_map)},
        'held_out': {
            'cell_midpoints': True, 'lhs_count': 0, 'seed': 7,
            'off_centre': [{'u': .25}],
        },
    }
    first_solve_started = threading.Event()
    caller_edit_finished = threading.Event()
    observed = []

    monkeypatch.setattr(
        'rmgpy.tools.eedf.generate._validate_spec', lambda candidate: None)
    monkeypatch.setattr(
        'rmgpy.tools.eedf.generate._model_inputs', lambda candidate, mapping: model)
    monkeypatch.setattr(
        'rmgpy.tools.eedf.generate._validate_map', lambda candidate: None)
    monkeypatch.setattr(
        'rmgpy.tools.eedf.generate.LoKIDriver', lambda candidate: object())
    monkeypatch.setattr(
        'rmgpy.tools.eedf.validation._compare_row',
        lambda predicted, reference, manifest: [])

    def solve(driver, mapping, candidate, coordinate, u_values, label):
        observed.append(candidate['envelopes'][field]['reference'])
        if len(observed) == 1:
            first_solve_started.set()
            assert caller_edit_finished.wait(5.)
        return [table.row(u_values[0], coordinate)._as_dict()]

    monkeypatch.setattr('rmgpy.tools.eedf.generate._solve', solve)
    criteria = {
        'mobility_N': {'rtol': 1.e-3, 'atol': 0.},
        'field_power': {'rtol': 1.e-3, 'atol': 0.},
    }

    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
        future = executor.submit(
            qualify_interpolation, table, spec, criteria=criteria)
        assert first_solve_started.wait(5.)
        next_spec = spec.copy()
        next_spec[field] = changed
        next_spec['envelopes'][field].update(
            min=changed, max=changed, reference=changed)
        caller_edit_finished.set()
        results = future.result()

    assert spec[field] == pinned
    assert observed == [pinned, pinned]
    assert [result['point'] for result in results] == retained
    assert all(result['passed'] for result in results)


@pytest.mark.parametrize(
    'field, pinned, changed', [
        ('Tg_K', 298.15, 300.),
        ('P_Pa', 666.61184, 700.),
    ],
)
def test_generate_owns_caller_spec_through_held_out_solves_and_publication(
        tmp_path, monkeypatch, field, pinned, changed):
    spec_path = os.environ.get('EEDF_REAL_SPEC')
    if not spec_path:
        pytest.skip('EEDF_REAL_SPEC supplies the complete pinned fixture')
    from rmgpy.tools.eedf.schema import load_spec

    spec = load_spec(spec_path)
    spec['axes'] = {'u': [0., 1.]}
    spec['held_out'].update(
        cell_midpoints=True, lhs_count=0, seed=7,
        off_centre=[{'u': .25}])
    spec['refinement']['max_rounds'] = 0
    spec['max_parallel'] = 1
    spec['scratch_root'] = str(tmp_path / 'scratch')
    spec['output_root'] = str(tmp_path / 'tables')

    fixture = Path(__file__).parent / 'fixtures' / 'ar'
    direct = loki.parse_output(fixture, 1.e-8)
    direct['composition'] = {}
    channel_map = json.loads(Path(spec['channel_map']['path']).read_text())
    loki._enrich_row(direct, channel_map, spec)
    direct['product_fractions'] = np.zeros(len(channel_map))
    direct['power_absolute'] = 0.
    first_held_started = threading.Event()
    caller_edit_finished = threading.Event()
    observed = []

    def solve(driver, mapping, local_spec, coordinates, u_values, label):
        if label.endswith('_held0'):
            observed.append(local_spec['envelopes'][field]['reference'])
            first_held_started.set()
            assert caller_edit_finished.wait(5.)
        elif '_held' in label:
            observed.append(local_spec['envelopes'][field]['reference'])
        rows = []
        for u in u_values:
            row = copy.deepcopy(direct)
            row['EN_Td'] = float(np.exp(u))
            row['u'] = float(u)
            row['composition'] = dict(coordinates)
            row['setup'] = label
            row['swarm']['mean_energy_eV'] = 1. + row['u']
            row['power_groups']['field'] = (
                row['swarm']['mobility_N'] * (row['EN_Td'] * 1.e-21) ** 2)
            rows.append(row)
        return rows

    monkeypatch.setattr(generation, '_solve', solve)
    monkeypatch.setattr(generation, '_screen', lambda *args: [])
    monkeypatch.setattr(generation, '_compare_row', lambda *args: [{
        'passed': True, 'rule': 'H1', 'relative_error': 0.,
        'quantity': 'fixture', 'absolute_error': 0., 'allowed_error': 1.,
    }])

    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
        future = executor.submit(
            generate, spec, command_line=['pytest generation ownership'])
        assert first_held_started.wait(5.)
        preview = spec.copy()
        preview[field] = changed
        preview['envelopes'][field].update(
            min=changed, max=changed, reference=changed)
        caller_edit_finished.set()
        path = future.result()

    saved = load_spec(path / 'generation_spec.json')
    manifest = json.loads((path / 'manifest.json').read_text())
    assert observed == [pinned, pinned]
    assert saved['envelopes'][field]['reference'] == pinned
    assert manifest['row_inputs']['envelopes'][field]['reference'] == pinned
    assert spec['envelopes'][field]['reference'] == pinned
    spec['envelopes'][field]['reference'] = changed
    assert saved['envelopes'][field]['reference'] == pinned


@pytest.mark.parametrize(
    'field, pinned, changed', [
        ('Tg_K', 298.15, 300.),
        ('P_Pa', 666.61184, 700.),
    ],
)
def test_loki_driver_owns_direct_caller_spec(field, pinned, changed):
    spec_path = os.environ.get('EEDF_REAL_SPEC')
    if not spec_path:
        pytest.skip('EEDF_REAL_SPEC supplies the complete pinned fixture')
    from rmgpy.tools.eedf.schema import load_spec

    spec = load_spec(spec_path)
    driver = loki.LoKIDriver(spec)
    preview = spec.copy()
    preview[field] = changed
    preview['envelopes'][field].update(
        min=changed, max=changed, reference=changed)

    assert driver.spec[field] == pinned
    assert driver.spec['envelopes'][field]['reference'] == pinned
    assert driver.spec is not spec
    assert driver.spec['envelopes'] is not spec['envelopes']


@pytest.mark.parametrize(
    'field, pinned, changed, setup_name', [
        ('Tg_K', 298.15, 300., 'gasTemperature'),
        ('P_Pa', 666.61184, 700., 'gasPressure'),
    ],
)
def test_loki_driver_owns_coordinates_through_real_solve_and_row_labelling(
        tmp_path, field, pinned, changed, setup_name):
    spec_path = os.environ.get('EEDF_REAL_SPEC')
    if not spec_path:
        pytest.skip('EEDF_REAL_SPEC supplies the complete pinned fixture')
    from rmgpy.tools.eedf.schema import load_spec

    spec = load_spec(spec_path)
    spec['scratch_root'] = str(tmp_path)
    driver = loki.LoKIDriver(spec)
    points = [{field: pinned}]
    fields = [16.99, 17.01]
    work = tmp_path / 'coordinate-ownership'

    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
        future = executor.submit(
            driver.run, points[0], fields, work.name)
        deadline = time.monotonic() + 60.
        while not (work / 'setup.in').exists():
            if future.done():
                future.result()
            assert time.monotonic() < deadline, 'setup timeout'
            time.sleep(.001)
        preview = points.copy()
        preview[0][field] = changed
        rows = future.result()

    setup = (work / 'setup.in').read_text()
    actual = float(re.search(
        rf'\b{setup_name}:\s*([^\n]+)', setup).group(1))
    assert actual == pinned
    assert [row['EN_Td'] for row in rows] == fields
    assert all(row['composition'][field] == pinned for row in rows)


def test_loki_driver_owns_field_sequence_through_real_solve(tmp_path):
    spec_path = os.environ.get('EEDF_REAL_SPEC')
    if not spec_path:
        pytest.skip('EEDF_REAL_SPEC supplies the complete pinned fixture')
    from rmgpy.tools.eedf.schema import load_spec

    spec = load_spec(spec_path)
    spec['scratch_root'] = str(tmp_path)
    driver = loki.LoKIDriver(spec)
    fields = [16.99, 17.01]
    requested = list(fields)
    work = tmp_path / 'field-ownership'

    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
        future = executor.submit(driver.run, {}, fields, work.name)
        deadline = time.monotonic() + 60.
        while not (work / 'setup.in').exists():
            if future.done():
                future.result()
            assert time.monotonic() < deadline, 'setup timeout'
            time.sleep(.001)
        fields[:] = [18.99, 19.01]
        rows = future.result()

    assert [row['EN_Td'] for row in rows] == requested


@pytest.mark.parametrize('shallow', [False, True], ids=['deep-control', 'shallow-preview'])
def test_loki_driver_owns_retained_spec_through_real_solve(tmp_path, shallow):
    spec_path = os.environ.get('EEDF_REAL_SPEC')
    if not spec_path:
        pytest.skip('EEDF_REAL_SPEC supplies the complete pinned fixture')
    from rmgpy.tools.eedf.schema import load_spec

    spec = load_spec(spec_path)
    spec['scratch_root'] = str(tmp_path)
    driver = loki.LoKIDriver(spec)
    fields = [16.99, 17.01]
    work = tmp_path / 'receiver-ownership'

    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
        future = executor.submit(driver.run, {}, fields, work.name)
        deadline = time.monotonic() + 60.
        while not (work / 'setup.in').exists():
            if future.done():
                future.result()
            assert time.monotonic() < deadline, 'setup timeout'
            time.sleep(.001)
        preview = driver.spec.copy() if shallow else copy.deepcopy(driver.spec)
        preview['solver_options']['numerics']['energyGrid']['cellNumber'] += 1
        rows = future.result()

    assert [row['EN_Td'] for row in rows] == fields
    assert len(list((work / 'output/solve').glob('reducedField_*/eedf.txt'))) == 2


def test_generate_owns_spec_before_first_validation(monkeypatch):
    spec_path = os.environ.get('EEDF_REAL_SPEC')
    if not spec_path:
        pytest.skip('EEDF_REAL_SPEC supplies the complete pinned fixture')
    from rmgpy.tools.eedf.schema import load_spec

    class StopValidation(Exception):
        pass

    spec = load_spec(spec_path)

    def validate(candidate):
        spec['envelopes']['Tg_K']['reference'] = 300.
        assert candidate is not spec
        assert candidate['envelopes']['Tg_K']['reference'] == 298.15
        raise StopValidation

    monkeypatch.setattr(generation, '_validate_spec', validate)
    with pytest.raises(StopValidation):
        generate(spec)


def test_loki_driver_owns_spec_before_first_check(monkeypatch):
    spec_path = os.environ.get('EEDF_REAL_SPEC')
    if not spec_path:
        pytest.skip('EEDF_REAL_SPEC supplies the complete pinned fixture')
    from rmgpy.tools.eedf.schema import load_spec

    class StopCheck(Exception):
        pass

    spec = load_spec(spec_path)

    def check(candidate):
        spec['envelopes']['Tg_K']['reference'] = 300.
        assert candidate is not spec
        assert candidate['envelopes']['Tg_K']['reference'] == 298.15
        raise StopCheck

    monkeypatch.setattr(loki, '_check_solver', check)
    with pytest.raises(StopCheck):
        loki.LoKIDriver(spec)


def test_qualify_interpolation_refuses_unpinned_multibranch_reference(
        tmp_path, monkeypatch):
    path, model = fixture_table(tmp_path, branches=2)
    table = EEDFTable.load(
        path, model, artifact_sha256=_file_hash(path / 'table.h5'))
    channel_map = tmp_path / 'channels.json'
    channel_map.write_text('[]')
    spec = {'channel_map': {'path': str(channel_map), 'sha256': _file_hash(channel_map)},
            'held_out': {
        'cell_midpoints': True, 'lhs_count': 0, 'seed': 7,
        'off_centre': [{'u': .25}],
    }}
    criteria = {
        'mobility_N': {'rtol': 1.e-3, 'atol': 0.},
        'field_power': {'rtol': 1.e-3, 'atol': 0.},
    }
    monkeypatch.setattr('rmgpy.tools.eedf.generate._validate_spec', lambda candidate: None)
    monkeypatch.setattr('rmgpy.tools.eedf.generate._model_inputs', lambda candidate, mapping: model)

    with pytest.raises(EEDFError, match='multi-branch.*pinned direct-solver'):
        qualify_interpolation(table, spec, criteria=criteria)


def test_real_pure_ar_table_passes_legacy_checks_but_awaits_t10(
        tmp_path, monkeypatch):
    spec_path = os.environ.get('EEDF_REAL_SPEC')
    if not spec_path:
        pytest.skip('EEDF_REAL_SPEC supplies the pinned real LoKI integration fixture')
    from rmgpy.tools.eedf.schema import load_spec
    spec = load_spec(spec_path)
    spec.pop('interpolation_accuracy', None)
    spec['axes']['u'] = np.linspace(np.log(16), np.log(18), 9).tolist()
    spec['held_out']['lhs_count'] = 4
    spec['held_out']['off_centre'] = [{'u': spec['axes']['u'][0]}]
    spec['scratch_root'] = str(tmp_path / 'loki')
    spec['output_root'] = str(tmp_path / 'tables')
    spec['refinement']['max_rounds'] = 0
    path = generate(spec, command_line=['pytest real LoKI fixture'])
    table = EEDFTable.load(
        path, model_inputs(spec), artifact_sha256=_file_hash(path / 'table.h5'),
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
    expected = {'schema_version', 'provider_kind', 'row_inputs', 'fingerprint', 'axes', 'energy_eV', 'energy_edges_eV',
                'envelopes', 'screen_results', 'tolerances', 'floors', 'channel_map', 'repositories',
                'solver', 'setups', 'state_population_assumptions', 'density_convention', 'units',
                'generator_command', 'generated_at', 'generation_seconds', 'held_out',
                'held_out_verdicts', 'accepted', 'branch_detection', 'limitations',
                'refinement_history', 'branches', 'artifact_sha256', 'branch_certification',
                'qualification_setup_sha256',
                'A6b-source', 'A6b-runtime', 'T10 interpolation qualification'}
    assert set(table.manifest) == expected
    assert table.manifest['provider_kind'] == 'loki-table'
    assert 'cold' in table.manifest['branch_detection'][0]['branch_sources']['branch_0']


@pytest.mark.parametrize('field', ['floors', 'tolerances', 'screen', 'refinement', 'held_out'])
def test_missing_generation_policy_refuses(field):
    from rmgpy.tools.eedf.schema import load_spec, _validate_spec
    spec_path = os.environ.get('EEDF_REAL_SPEC')
    if not spec_path:
        pytest.skip('EEDF_REAL_SPEC supplies the complete fixture')
    spec = load_spec(spec_path)
    del spec[field]
    with pytest.raises(SpecError, match='spec fields'):
        _validate_spec(spec)
