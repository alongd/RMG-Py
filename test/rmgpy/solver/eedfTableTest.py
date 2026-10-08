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
from pathlib import Path
import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator

from rmgpy.solver.eedf import (
    EEDFTable, OutOfDomain, EnvelopeBreach, AmbiguousBranch,
    IllConditionedCoordinate, field_power_from_mobility,
)
from rmgpy.tools.eedf.artifact import _write_artifact
from rmgpy.tools.eedf.schema import _file_hash, _content_hash, _interpolant_identity, FingerprintMismatch


def fixture_table(tmp_path, energies=(1., 2., 3., 4.), branches=1, composition=False,
                  transport_identity=True, corrupt_native_power=False,
                  include_field=True, composition_axis=(0., 1.),
                  nonlinear_mobility=False, retained_held_out=()):
    axes = {'u': [0., 1., 2., 3.]}
    if composition:
        axes['x'] = list(composition_axis)
    shape = tuple(len(v) for v in axes.values())
    def values(row):
        arr = np.empty(shape, dtype=object)
        for index in np.ndindex(shape):
            energy = energies[index[0]] + (index[1] if composition else 0)
            rate = (1., 2., 4., 8.)[index[0]] * 1e-15 * (1 + (index[1] if composition else 0))
            if nonlinear_mobility:
                mobility = (1. + 0.07 * index[0] ** 2 + 0.12 * index[1] ** 2 +
                            0.03 * index[0] * index[1]) * 1.e24
            else:
                mobility = (1.0 + 0.15 * index[0] +
                            (0.2 * index[1] if composition else 0.0)) * 1e24
            field = (mobility * np.exp(2. * axes['u'][index[0]]) * 1e-42
                     if transport_identity else rate * 2)
            if corrupt_native_power and index[0] == 1:
                field *= 1.01
            power_groups = {'field': field} if include_field else {}
            arr[index] = {'u': axes['u'][index[0]], 'EN_Td': np.exp(axes['u'][index[0]]),
                'swarm': {'mean_energy_eV': energy, 'mobility_N': mobility, 'diffusion_N': 2e24},
                'k_ine': np.array([rate]), 'k_sup': np.array([0.]),
                'channel_power': np.array([rate * 2]), 'attachment_energy_eV': np.array([0.]),
                'target_fractions': np.array([1.]), 'product_fractions': np.array([0.]),
                'rate_floors': np.array([1e-40]), 'below_floor': np.array([[False, True]]),
                'power_groups': power_groups, 'f0': np.array([1.]),
                'energy_eV': np.array([.5]), 'energy_edges_eV': np.array([0., 1.]),
                'gas_fractions': {'Ar': 1.}, 'state_populations': {'Ar(1S0)': 1.},
                'composition': {'x': axes['x'][index[1]]} if composition else {},
                'converged': True, 'iteration_count': -1, 'setup': 'fixture'}
        return arr
    inputs = {'axes': axes, 'Tg_K': 298.15, 'P_Pa': 666.61,
              'arm': {'gas': 'Ar'}, 'cross_sections': {'Ar': 'sha'},
              'channel_map_sha256': 'map', 'reactions': [], 'solver_options': {'growth': 'temporal'},
              'interpolant': _interpolant_identity(), 'state_populations': {'Ar(1S0)': 1.}}
    manifest = {'accepted': True,
        'T10 interpolation qualification': {
            'check': 'T10 interpolation qualification',
            'criterion': 'fixture frozen accuracy criteria',
            'verdict_count': 1, 'passed': True},
        'held_out_verdicts': [{'point': {'u': .5}, 'branch_id': 'branch_0', 'passed': True, 'checks': [{'passed': True}]}],
        'branch_certification': {f'branch_{i}': 'uncertified: unseeded scans' for i in range(branches)}, 'schema_version': 1, 'row_inputs': inputs, 'fingerprint': _content_hash(inputs),
        'axes': axes, 'envelopes': {'Tg_K': {'min': 298.15, 'max': 298.15, 'reference': 298.15},
                                'metastable': {'min': 0., 'max': 1e-5, 'reference': 0.}},
        'tolerances': {'fingerprint_rtol': 1e-9, 'u_condition': {'min_abs_denergy_du': 1e-6, 'min_scaled_denergy_du': 1e-6}},
        'floors': {'rate_absolute': 1e-40}, 'repositories': {'rmgpy': 'old', 'database': 'old'},
        'channel_map': [{'description': 'fixture', 'kind': 'excitation', 'classification': 'B'}],
        'energy_eV': [0.5], 'energy_edges_eV': [0., 1.]}
    held_out = []
    for coordinates in retained_held_out:
        row = copy.deepcopy(values(0).flat[0])
        row['u'] = coordinates['u']
        row['EN_Td'] = np.exp(coordinates['u'])
        row['composition'] = {name: coordinates[name] for name in axes if name != 'u'}
        held_out.append(row)
    manifest['held_out'] = copy.deepcopy(list(retained_held_out))
    path = _write_artifact(
        tmp_path, manifest,
        {f'branch_{i}': values(i) for i in range(branches)}, held_out)
    return path, inputs


def test_one_row_returns_monotone_positive_rates_and_consistent_derivative(tmp_path):
    path, model = fixture_table(tmp_path, composition=True)
    table = EEDFTable.load(path, model, artifact_sha256=_file_hash(path / 'table.h5'))
    rows = [table.row(u, {'x': .4}) for u in np.linspace(0, 3, 61)]
    assert all(row.k_ine[0] > 0 and row.k_sup[0] == 0 for row in rows)
    assert np.all(np.diff([r.k_ine[0] for r in rows]) > 0)
    assert rows[30].swarm['mean_energy_eV'] == pytest.approx(2.9)
    assert rows[30].denergy_du == pytest.approx(1.)
    assert rows[30].branch_id == 'branch_0'
    assert rows[30].below_floor[0, 1]


@pytest.mark.parametrize('shallow', [False, True], ids=['deep-control', 'shallow-preview'])
def test_public_axes_preview_cannot_mutate_table_state(tmp_path, shallow):
    path, model = fixture_table(tmp_path, composition=True)
    table = EEDFTable.load(
        path, model, artifact_sha256=_file_hash(path / 'table.h5'))
    expected = table.row(.5, {'x': .4})

    preview = table.axes.copy() if shallow else copy.deepcopy(table.axes)
    np.exp(preview[0], out=preview[0])
    actual = table.row(.5, {'x': .4})

    assert np.array_equal(table.axes[0], np.arange(4.))
    assert actual.swarm['mean_energy_eV'] == expected.swarm['mean_energy_eV']
    assert actual.k_ine == pytest.approx(expected.k_ine, rel=0., abs=0.)


def test_off_grid_field_power_is_derived_from_interpolated_mobility(tmp_path):
    path, model = fixture_table(tmp_path, composition=True, transport_identity=True)
    table = EEDFTable.load(path, model, artifact_sha256=_file_hash(path / 'table.h5'))

    row = table.row(0.5, {'x': 0.4})
    derived = row.swarm['mobility_N'] * (row.EN_Td * 1.e-21) ** 2
    u_nodes = np.arange(4., dtype=float)
    mobility_nodes = (1.08 + 0.15 * u_nodes) * 1.e24
    native_power_nodes = mobility_nodes * np.exp(2. * u_nodes) * 1.e-42
    independently_interpolated = PchipInterpolator(u_nodes, native_power_nodes)(0.5)
    independent_residual = abs(independently_interpolated - derived) / derived

    assert independent_residual > 1.e-6  # RED: independent PCHIP columns disagree off grid.
    assert row.EN_Td == pytest.approx(np.exp(0.5), rel=0., abs=0.)
    assert row.power_groups['field'] == pytest.approx(derived, rel=1.e-14, abs=0.)


def test_field_power_derivatives_match_finite_differences_in_u_and_composition(tmp_path):
    path, model = fixture_table(
        tmp_path, composition=True, transport_identity=True,
        composition_axis=(0., .5, 1.), nonlinear_mobility=True)
    table = EEDFTable.load(path, model, artifact_sha256=_file_hash(path / 'table.h5'))
    u, composition, step = 0.7, {'x': 0.4}, 1.e-6
    row = table.row(u, composition)

    u_difference = (table.row(u + step, composition).power_groups['field'] -
                    table.row(u - step, composition).power_groups['field']) / (2. * step)
    x_high, x_low = {'x': composition['x'] + step}, {'x': composition['x'] - step}
    x_difference = (table.row(u, x_high).power_groups['field'] -
                    table.row(u, x_low).power_groups['field']) / (2. * step)
    scale = (row.EN_Td * 1.e-21) ** 2
    explicit_formula = scale * (
        row.derivatives['u']['swarm']['mobility_N'] +
        2. * row.swarm['mobility_N'])

    assert row.derivatives['u']['power_groups']['field'] == pytest.approx(
        explicit_formula, rel=1.e-14, abs=0.)
    assert row.derivatives['u']['power_groups']['field'] == pytest.approx(
        u_difference, rel=1.e-9, abs=0.)
    assert row.derivatives['composition']['x']['power_groups']['field'] == pytest.approx(
        x_difference, rel=1.e-9, abs=0.)


def test_corrupted_native_power_fails_source_but_runtime_identity_still_passes(tmp_path):
    path, model = fixture_table(
        tmp_path, transport_identity=True, corrupt_native_power=True)
    with pytest.raises(FingerprintMismatch, match='A6b-source'):
        EEDFTable.load(
            path, model, artifact_sha256=_file_hash(path / 'table.h5'))
    table = EEDFTable.load(
        path, model, artifact_sha256=_file_hash(path / 'table.h5'),
        require_accepted=False)

    row = table.row(0.5, {})

    assert table.a6b_source['check'] == 'A6b-source'
    assert not table.a6b_source['passed']
    assert row.a6b_runtime == {
        'check': 'A6b-runtime',
        'criterion': 'field power derived from interpolated mobility at query state',
        'relative_error': 0.0,
        'passed': True,
    }
    assert row.power_groups['field'] == pytest.approx(
        row.swarm['mobility_N'] * (row.EN_Td * 1.e-21) ** 2,
        rel=1.e-14, abs=0.)


def test_missing_native_field_power_refuses(tmp_path):
    path, model = fixture_table(tmp_path, include_field=False)
    with pytest.raises(FingerprintMismatch, match='required native field power'):
        EEDFTable.load(path, model, artifact_sha256=_file_hash(path / 'table.h5'))


@pytest.mark.parametrize(
    'field, message', [
        ('EN_Td', 'native EN_Td state mapping'),
        ('u', 'native u state mapping'),
    ],
)
def test_native_field_state_mapping_refuses_mismatched_coordinate(
        tmp_path, field, message):
    import h5py
    path, model = fixture_table(tmp_path, transport_identity=True)
    with h5py.File(path / 'table.h5', 'r+') as table:
        table['branches/branch_0/' + field][1] += .01
    digest = _file_hash(path / 'table.h5')
    manifest_path = path / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    manifest['artifact_sha256'] = digest
    manifest_path.write_text(json.dumps(manifest))

    with pytest.raises(FingerprintMismatch, match=message):
        EEDFTable.load(path, model, artifact_sha256=digest)


def test_domain_envelope_branch_and_fold_refusals(tmp_path):
    path, model = fixture_table(tmp_path)
    table = EEDFTable.load(path, model, artifact_sha256=_file_hash(path / 'table.h5'))
    with pytest.raises(OutOfDomain): table.row(-.1, {})
    with pytest.raises(OutOfDomain): table.row(4., {})
    with pytest.raises(EnvelopeBreach): table.row(1., {'metastable': 2e-5})
    with pytest.raises(EnvelopeBreach): table.row(1., {'Tg_K': 300})
    with pytest.raises(EnvelopeBreach): table.row(1., {'unknown': 0.})
    path, model = fixture_table(tmp_path / 'two', branches=2)
    table = EEDFTable.load(path, model, artifact_sha256=_file_hash(path / 'table.h5'))
    with pytest.raises(AmbiguousBranch): table.row(1., {})
    assert table.row(1., {}, 'branch_1').branch_id == 'branch_1'
    path, model = fixture_table(tmp_path / 'fold', energies=(1., 2., 1., 3.))
    with pytest.raises(IllConditionedCoordinate): EEDFTable.load(path, model, artifact_sha256=_file_hash(path / 'table.h5')).row(1.5, {})
    path, model = fixture_table(tmp_path / 'flat', energies=(1., 1., 1., 1.))
    with pytest.raises(IllConditionedCoordinate): EEDFTable.load(path, model, artifact_sha256=_file_hash(path / 'table.h5')).row(1.5, {})


@pytest.mark.parametrize('field', ['axes', 'Tg_K', 'P_Pa', 'arm', 'cross_sections',
    'channel_map_sha256', 'reactions', 'solver_options', 'interpolant', 'state_populations'])
def test_each_model_fingerprint_mismatch_refuses(tmp_path, field):
    path, model = fixture_table(tmp_path)
    changed = copy.deepcopy(model)
    changed[field] = changed[field] + 1 if isinstance(changed[field], float) else {'changed': True}
    with pytest.raises(FingerprintMismatch, match=field): EEDFTable.load(path, changed, artifact_sha256=_file_hash(path / 'table.h5'))


def test_artifact_and_manifest_tampering_refuse(tmp_path):
    path, model = fixture_table(tmp_path)
    manifest_path = path / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    import h5py
    with h5py.File(path / 'table.h5', 'r+') as h5:
        h5['branches/branch_0/k_ine'][...] *= 1e6
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(FingerprintMismatch): EEDFTable.load(path, model, artifact_sha256=_file_hash(path / 'table.h5'))


def retained_points():
    """Historical held-out coordinates, including an exactly distinct neighbour."""
    midpoint = 0.5
    return [
        {'u': midpoint},
        {'u': np.nextafter(midpoint, np.inf)},
        {'u': 1.5},
    ]


def test_loader_preserves_historical_and_one_ulp_distinct_held_out_points(tmp_path):
    points = retained_points()
    path, model = fixture_table(tmp_path, retained_held_out=points)

    table = EEDFTable.load(path, model, artifact_sha256=_file_hash(path / 'table.h5'))

    assert table.manifest['held_out'] == points
    assert table.manifest['held_out'][0]['u'].hex() != table.manifest['held_out'][1]['u'].hex()


@pytest.mark.parametrize('corrupt', [
    pytest.param(lambda points: points.pop(0), id='dropped-midpoint'),
    pytest.param(lambda points: points[0].update(u=0.25), id='changed-coordinate'),
    pytest.param(lambda points: points[0].update(u=1.0), id='training-node-replacement'),
])
def test_loader_refuses_manifest_held_out_disagreement(tmp_path, corrupt):
    path, model = fixture_table(tmp_path, retained_held_out=retained_points())
    manifest_path = path / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    corrupt(manifest['held_out'])
    manifest_path.write_text(json.dumps(manifest))

    with pytest.raises(FingerprintMismatch, match='held-out coordinates'):
        EEDFTable.load(path, model, artifact_sha256=_file_hash(path / 'table.h5'))


@pytest.mark.parametrize('load_file', [False, True], ids=['directory', 'table-h5'])
def test_loader_refuses_duplicate_manifest_branch_entries(tmp_path, load_file):
    path, model = fixture_table(tmp_path)
    manifest_path = path / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    manifest['branches'].append(manifest['branches'][0])
    manifest_path.write_text(json.dumps(manifest))
    source = path / 'table.h5' if load_file else path

    with pytest.raises(FingerprintMismatch, match='duplicate branches'):
        EEDFTable.load(
            source, model, artifact_sha256=_file_hash(path / 'table.h5'),
            require_accepted=False)


@pytest.mark.parametrize('field', [
    'energy_eV', 'channel_map', 'tolerances', 'envelopes',
    'T10 interpolation qualification'])
def test_manifest_quantity_or_policy_tampering_refuses(tmp_path, field):
    path, model = fixture_table(tmp_path)
    manifest_path = path / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    if field == 'energy_eV': manifest[field] = [5.]
    elif field == 'channel_map': manifest[field][0]['description'] = 'wrong'
    elif field == 'tolerances': manifest[field]['fingerprint_rtol'] = 1
    elif field == 'envelopes': manifest[field]['metastable']['max'] = 1
    else: manifest[field]['passed'] = False
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(FingerprintMismatch):
        EEDFTable.load(
            path, model, artifact_sha256=_file_hash(path / 'table.h5'),
            require_accepted=field != 'T10 interpolation qualification')


def test_repository_shas_are_reported_without_invalidating_current_model(tmp_path):
    path, model = fixture_table(tmp_path)
    table = EEDFTable.load(path, model, artifact_sha256=_file_hash(path / 'table.h5'))
    assert table.manifest['repositories']['rmgpy'] == 'old'
    assert table.row(1., {}).k_ine[0] == pytest.approx(2e-15, rel=1e-12, abs=0)


def test_unqualified_held_out_verdict_refuses_production_load(tmp_path):
    path, model = fixture_table(tmp_path)
    manifest_path = path / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    manifest['accepted'] = False
    manifest['held_out_verdicts'][0]['passed'] = False
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(FingerprintMismatch, match='held-out qualification'):
        EEDFTable.load(path, model, artifact_sha256=_file_hash(path / 'table.h5'))


def test_loader_owns_acceptance_flag_before_manifest_io(tmp_path, monkeypatch):
    """A later caller edit cannot turn a production load into a provisional one."""
    path, model = fixture_table(tmp_path)
    manifest_path = path / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    manifest['accepted'] = False
    manifest['held_out_verdicts'][0]['passed'] = False
    manifest_path.write_text(json.dumps(manifest))
    import h5py
    qualification = {
        key: manifest.get(key) for key in (
            'accepted', 'held_out_verdicts', 'screen_results',
            'branch_detection', 'A6b-source', 'A6b-runtime',
            'T10 interpolation qualification')}
    with h5py.File(path / 'table.h5', 'r+') as table:
        table.attrs['qualification_sha256'] = _content_hash(qualification)
    digest = _file_hash(path / 'table.h5')
    manifest['artifact_sha256'] = digest
    manifest_path.write_text(json.dumps(manifest))
    storage = np.array([True])
    flag = storage[:1].reshape(())
    import rmgpy.solver.eedf as eedf
    original = eedf.json.load

    def load_and_edit(stream):
        result = original(stream)
        storage[0] = False
        return result

    monkeypatch.setattr(eedf.json, 'load', load_and_edit)
    with pytest.raises(FingerprintMismatch, match='held-out qualification'):
        EEDFTable.load(
            path, model, artifact_sha256=digest,
            require_accepted=flag)


def test_loader_refuses_nonscalar_acceptance_flag(tmp_path):
    path, model = fixture_table(tmp_path)
    with pytest.raises((TypeError, ValueError), match='require_accepted'):
        EEDFTable.load(
            path, model, artifact_sha256=_file_hash(path / 'table.h5'),
            require_accepted=np.array([True]))


def test_field_power_owns_each_input_before_later_conversion():
    mobility = np.array([2.])

    class EditingField:
        def __array__(self, dtype=None, copy=None):
            mobility[0] = 20.
            return np.array([3.], dtype=dtype)

    result = field_power_from_mobility(mobility, EditingField())

    assert result[0] == pytest.approx(2. * (3.e-21) ** 2)
    assert not result.flags.writeable
