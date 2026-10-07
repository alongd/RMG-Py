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
from collections.abc import Mapping
import gc
import json
from pathlib import Path
import threading

import numpy as np
import pytest

from rmgpy.solver.eedf import (
    AmbiguousBranch,
    EEDFError,
    EnvelopeBreach,
    IllConditionedCoordinate,
    OutOfDomain,
)
from rmgpy.solver.eedf_provider import (
    DEVELOPMENT_DIAGNOSTIC_STATUS,
    DEVELOPMENT_UNQUALIFIED_STATUS,
    EEDFProvider,
    development_unqualified_table_route,
)
from rmgpy.tools.eedf.artifact import write_artifact
from rmgpy.tools.eedf.loki import qualification_setup_sha256
from rmgpy.tools.eedf.schema import (
    FingerprintMismatch,
    content_hash,
    file_hash,
    interpolant_identity,
)


def fixture_provider_table(tmp_path, composition=False, branches=1, decreasing=False,
                           accepted=True, every_quantity_tolerance=True,
                           a2_rtol=1e-6, setup_fingerprint=True,
                           qualification_quantity_rules=True,
                           omitted_qualification_rule=None,
                           qualification_rule_overrides=None,
                           a6_atol=0., weight_spelling='1',
                           second_classification='A'):
    tmp_path.mkdir(parents=True, exist_ok=True)
    ar_input = tmp_path / 'Ar.fixture'
    mass_input = tmp_path / 'masses.fixture'
    ar_input.write_text('fixture cross section\n')
    mass_input.write_text('fixture masses\n')
    axes = {'u': [0., 1., 2., 3.]}
    if composition:
        axes['x'] = [0., 1.]
    shape = tuple(len(axis) for axis in axes.values())

    def branch_values():
        values = np.empty(shape, dtype=object)
        for index in np.ndindex(shape):
            u_index = index[0]
            x_index = index[1] if composition else 0
            scale = (1., 2., 4., 8.)[u_index] * (1. + x_index)
            composition_values = {'x': axes['x'][x_index]} if composition else {}
            values[index] = {
                'u': axes['u'][u_index],
                'EN_Td': np.exp(axes['u'][u_index]),
                'swarm': {'mean_energy_eV': ((4. - u_index) if decreasing else (1. + u_index)) + x_index,
                          'characteristic_energy_eV': .5 * scale,
                          'mobility_N': 1e24 * scale,
                          'diffusion_N': 2e24 * scale},
                'k_ine': np.array([scale, 10. * scale]) * 1e-15,
                'k_sup': np.array([1.5 * scale, 15. * scale]) * 1e-15,
                'channel_power': np.array([2. * scale, 20. * scale]) * 1e-15,
                'attachment_energy_eV': np.zeros(2),
                'target_fractions': np.ones(2),
                'product_fractions': np.zeros(2),
                'rate_floors': np.full(2, 1e-40),
                'below_floor': np.zeros((2, 2), dtype=bool),
                'power_groups': {'field': scale, 'balance': .01 * scale},
                'f0': np.array([1.]),
                'energy_eV': np.array([.5]),
                'energy_edges_eV': np.array([0., 1.]),
                'gas_fractions': {'Ar': 1.},
                'state_populations': {'Ar(1S0)': 1.},
                'composition': composition_values,
                'converged': True,
                'iteration_count': -1,
                'setup': 'provider fixture',
            }
        return values

    identities = [
        {'library': 'fixture', 'index': 1, 'repr': 'first'},
        {'library': 'fixture', 'index': 2, 'repr': 'second'},
    ]
    model = {
        'axes': axes,
        'Tg_K': 298.15,
        'P_Pa': 666.61,
        'arm': {'gas': 'Ar'},
        'loki_commit': 'fixture-loki-commit',
        'binary_sha256': 'fixture-binary-sha256',
        'cross_sections': {'Ar': file_hash(ar_input)},
        'auxiliary_inputs': {'masses': file_hash(mass_input)},
        'channel_map_sha256': 'map',
        'reactions': identities,
        'solver_options': {
            'growthModelType': 'temporal',
            'ionizationOperatorType': 'usingSDCS',
            'includeEECollisions': False,
        },
        'working_conditions': {'electronDensity': '1e16'},
        'interpolant': interpolant_identity(),
        'state_populations': {'Ar(1S0)': 1.},
    }
    rules = {
        'EN_Td': 'A2',
        'swarm.characteristic_energy_eV': 'A2',
        'power.field': 'A4',
        'power.balance': 'A4',
        'channel_power_fraction': 'A5',
        'attachment_energy_eV': 'A1',
        'target_fractions': 'A5',
        'product_fractions': 'A5',
        'f0.weighted_L1': 'A5',
    }
    rules.update(qualification_rule_overrides or {})
    if omitted_qualification_rule is not None:
        del rules[omitted_qualification_rule]
    if qualification_quantity_rules:
        model['qualification_quantity_rules'] = rules
    setup_spec = {
        'schema_version': 1,
        'arm': model['arm'],
        'Tg_K': model['Tg_K'],
        'P_Pa': model['P_Pa'],
        'loki_commit': model['loki_commit'],
        'binary': {'path': 'loki', 'sha256': model['binary_sha256']},
        'shared_objects': {'libfixture.so': 'fixture-shared-object-sha256'},
        'compiler': {'id': 'fixture-compiler'},
        'cmake_cache': {'path': 'CMakeCache.txt',
                        'sha256': 'fixture-cache-sha256'},
        'channel_map': {'path': 'channels.json', 'sha256': 'map'},
        'axes': axes,
        'envelopes': {
            'Tg_K': {'min': 298.15, 'max': 298.15, 'reference': 298.15},
            'metastable': {'min': 0., 'max': 1e-5, 'reference': 0.},
        },
        'gas_properties': {
            'fraction': ['Ar = 1'],
            'mass': 'masses',
        },
        'state_properties': {
            'population': ['Ar(1S0) = 1'],
            'energy': ['Ar(1S0) = 0'],
            'statisticalWeight': ['Ar(1S0) = ' + weight_spelling],
        },
        'solver_options': model['solver_options'],
        'working_conditions': model['working_conditions'],
        'input_files': {
            'Ar': {'path': str(ar_input), 'sha256': model['cross_sections']['Ar'],
                   'kind': 'cross_section'},
            'masses': {'path': str(mass_input),
                       'sha256': model['auxiliary_inputs']['masses'],
                       'kind': 'property'},
        },
        'timeout_s': 1.,
    }
    if qualification_quantity_rules:
        setup_spec['qualification_quantity_rules'] = rules
    tolerances = {
        'fingerprint_rtol': 1e-9,
        'u_condition': {'min_abs_denergy_du': 1e-6,
                        'min_scaled_denergy_du': 1e-6},
        'A1': {'rtol': 1e-6, 'atol': 1e-12},
        'A2': {'rtol': a2_rtol, 'atol': 0.},
        'A3': {'rtol': 1e-6, 'atol': 1e-40},
        'A4': {'rtol': 1e-6, 'atol': 1e-40},
        'A5': {'rtol': 1e-6, 'atol': 1e-12, 'share_min': 0.01},
        'A6': {'rtol': 1e-6, 'atol': a6_atol},
    }
    if every_quantity_tolerance:
        tolerances['G1'] = {'rtol': 1e-6, 'atol': 1e-12}
    manifest = {
        'accepted': accepted,
        'held_out_verdicts': [{'point': {'u': .5}, 'branch_id': 'branch_0',
                               'passed': accepted, 'checks': [{'passed': accepted}]}],
        'branch_certification': {
            f'branch_{index}': 'uncertified: unseeded scans'
            for index in range(branches)
        },
        'schema_version': 1,
        'row_inputs': model,
        'fingerprint': content_hash(model),
        'axes': axes,
        'envelopes': {
            'Tg_K': {'min': 298.15, 'max': 298.15, 'reference': 298.15},
            'metastable': {'min': 0., 'max': 1e-5, 'reference': 0.},
        },
        'tolerances': tolerances,
        'floors': {'rate_absolute': 1e-40, 'eedf_dynamic_range': 1.e-12,
                   'rate_flux_fraction': 1e-3,
                   'absolute_flux_fraction': 1e-4},
        'repositories': {'rmgpy': 'old', 'database': 'old'},
        'solver': {
            'commit': model['loki_commit'],
            'binary_sha256': model['binary_sha256'],
            'compiler': {'id': 'fixture-compiler'},
            'cmake_cache_sha256': 'fixture-cache-sha256',
            'options': model['solver_options'],
            'ee_setting': model['solver_options']['includeEECollisions'],
            'working_conditions': model['working_conditions'],
            'input_files': {
                'Ar': {'sha256': model['cross_sections']['Ar'],
                       'kind': 'cross_section'},
                'masses': {'sha256': model['auxiliary_inputs']['masses'],
                           'kind': 'property'},
            },
        },
        'channel_map': [
            {'description': 'first', 'kind': 'excitation', 'classification': 'A',
             'reaction': identities[0]},
            {'description': 'second', 'kind': 'excitation',
             'classification': second_classification,
             'reaction': identities[1]},
        ],
        'energy_eV': [.5],
        'energy_edges_eV': [0., 1.],
    }
    setup_spec['tolerances'] = tolerances
    setup_spec['floors'] = manifest['floors']
    if setup_fingerprint:
        model['qualification_setup_sha256'] = qualification_setup_sha256(setup_spec)
        manifest['fingerprint'] = content_hash(model)
        manifest['qualification_setup_sha256'] = model['qualification_setup_sha256']
    path = write_artifact(
        tmp_path,
        manifest,
        {f'branch_{index}': branch_values() for index in range(branches)},
        [],
    )
    setup_spec.update({
        'loki_commit': model['loki_commit'],
        'binary': {'path': str(path / 'loki'),
                   'sha256': model['binary_sha256']},
        'compiler': manifest['solver']['compiler'],
        'cmake_cache': {'sha256': manifest['solver']['cmake_cache_sha256']},
        'channel_map': {'path': str(path / 'channels.json'), 'sha256': 'map'},
    })
    (path / 'generation_spec.json').write_text(json.dumps(setup_spec))
    return path, model


def make_provider(path, model, **overrides):
    arguments = {
        'artifact_sha256': file_hash(path / 'table.h5'),
        'reaction_map': {4: (1, 'ine'), 7: (0, 'sup')},
        'branch': 'branch_0',
        'empirical_laws': {'Ar2+ recombination': {'evaluate_at': 'Te_eff'}},
    }
    arguments.update(overrides)
    return EEDFProvider(path, model, **arguments)


def test_pure_argon_and_composition_axis_use_the_same_provider_type(tmp_path):
    pure_path, pure_model = fixture_provider_table(tmp_path / 'pure')
    arm_path, arm_model = fixture_provider_table(tmp_path / 'arm', composition=True)

    pure = make_provider(pure_path, pure_model)
    arm = make_provider(arm_path, arm_model)

    assert type(pure) is EEDFProvider
    assert type(arm) is EEDFProvider
    assert pure.empirical_laws['Ar2+ recombination']['evaluate_at'] == 'Te_eff'
    assert arm.manifest['fingerprint'] == content_hash(arm_model)


def test_unqualified_table_is_development_context_only_and_carries_blockers(tmp_path):
    qualified_path, qualified_model = fixture_provider_table(
        tmp_path / 'qualified', accepted=True)
    with development_unqualified_table_route():
        with pytest.raises(EEDFError, match='already qualified'):
            make_provider(qualified_path, qualified_model)

    path, model = fixture_provider_table(tmp_path / 'unqualified', accepted=False)
    with pytest.raises(FingerprintMismatch, match='held-out qualification'):
        make_provider(path, model)

    with development_unqualified_table_route():
        provider = make_provider(path, model)

    notice = provider.development_notice
    assert notice['scientific_status'] == DEVELOPMENT_UNQUALIFIED_STATUS
    assert notice['artifact_blockers'] == [
        'artifact manifest is not accepted',
        '1 of 1 held-out verdicts failed',
        'branch_0: uncertified: unseeded scans',
    ]
    assert provider.manifest['scientific_status'] == DEVELOPMENT_UNQUALIFIED_STATUS
    with pytest.raises(EEDFError, match='refuses qualification'):
        provider.require_qualified('qualification')
    with pytest.raises(EEDFError, match='refuses export'):
        provider.require_qualified('export')

    # Leaving the context restores the production loader immediately.
    with pytest.raises(FingerprintMismatch, match='held-out qualification'):
        make_provider(path, model)


def test_complete_composition_identity_controls_the_one_row_cache(tmp_path):
    path, model = fixture_provider_table(tmp_path, composition=True)
    provider = make_provider(path, model)

    first = provider.row(1., {'x': 0.})
    first_identity = provider.cache_identity
    assert provider.row(1., {'x': 0.}) is first

    changed = provider.row(1., {'x': 1.})
    assert changed is not first
    assert provider.cache_identity != first_identity
    assert provider.reaction_rate(4, changed) == pytest.approx(2. * first.k_ine[1])
    assert provider.mean_energy_eV(changed) == pytest.approx(3.)
    assert provider.transport(changed)['mobility_N'] == pytest.approx(4e24)


def test_reaction_mapping_selects_ine_and_sup_and_refuses_bad_entries(tmp_path):
    path, model = fixture_provider_table(tmp_path)
    provider = make_provider(path, model)
    row = provider.row(1., {})

    assert provider.reaction_rate(4, row) == pytest.approx(row.k_ine[1])
    assert provider.reaction_rate(7, row) == pytest.approx(row.k_sup[0])
    with pytest.raises(KeyError, match='unmapped reaction index'):
        provider.reaction_rate(99, row)
    with pytest.raises(ValueError, match='channel column'):
        make_provider(path, model, reaction_map={0: (2, 'ine')})
    with pytest.raises(ValueError, match='channel side'):
        make_provider(path, model, reaction_map={0: (0, 'forward')})
    with pytest.raises(ValueError, match='reaction index'):
        make_provider(path, model, reaction_map={-1: (0, 'ine')})


def test_configuration_and_issued_rows_are_defensive(tmp_path):
    path, model = fixture_provider_table(tmp_path / 'one')
    other_path, other_model = fixture_provider_table(tmp_path / 'other')
    provider = make_provider(path, model)
    other = make_provider(other_path, other_model)
    row = provider.row(1., {})

    with pytest.raises(TypeError):
        provider.reaction_map[4] = (0, 'ine')
    manifest = provider.manifest
    manifest['fingerprint'] = 'mutated copy'
    assert provider.manifest['fingerprint'] != 'mutated copy'
    with pytest.raises(ValueError):
        row.k_ine[0] = 0.
    with pytest.raises(ValueError):
        row.k_ine.setflags(write=True)
    with pytest.raises(FingerprintMismatch, match='not issued'):
        other.reaction_rate(4, row)
    for u in np.linspace(0., 3., 50):
        provider.row(float(u), {})
    gc.collect()
    assert len(provider._issued_rows) <= 2


def test_cache_identity_and_row_are_published_atomically(tmp_path, monkeypatch):
    path, model = fixture_provider_table(tmp_path)
    provider = make_provider(path, model)
    first = provider.row(1., {})
    entered = threading.Event()
    release = threading.Event()
    original = provider._table.row

    def paused_row(u, composition, branch):
        if u == 2.:
            entered.set()
            assert release.wait(5.)
        return original(u, composition, branch)

    monkeypatch.setattr(provider._table, 'row', paused_row)
    writer = threading.Thread(target=lambda: provider.row(2., {}))
    writer.start()
    assert entered.wait(5.)

    reader_result = []
    reader = threading.Thread(target=lambda: reader_result.append(provider.row(1., {})))
    reader.start()
    assert reader.is_alive()
    release.set()
    writer.join(5.)
    reader.join(5.)

    assert len(reader_result) == 1
    assert reader_result[0].u == first.u
    assert provider.cache_identity == (1., ())
    assert provider.cached_row is reader_result[0]


def test_fingerprint_branch_and_named_domain_errors_are_preserved(tmp_path):
    path, model = fixture_provider_table(tmp_path, composition=True, branches=2)
    changed_model = copy.deepcopy(model)
    changed_model['Tg_K'] = 300.
    with pytest.raises(FingerprintMismatch, match='Tg_K'):
        make_provider(path, changed_model)
    with pytest.raises(AmbiguousBranch, match='declared branch missing'):
        make_provider(path, model, branch='missing')

    provider = make_provider(path, model)
    with pytest.raises(OutOfDomain, match=r'u=-0.1 outside \[0.0, 3.0\]'):
        provider.domain_check(-.1, {'x': .5}, y=np.zeros(1), context={'accepted': True})
    with pytest.raises(OutOfDomain, match=r'x=2.0 outside \[0.0, 1.0\]'):
        provider.domain_check(1., {'x': 2.})
    with pytest.raises(EnvelopeBreach, match='metastable'):
        provider.domain_check(1., {'x': .5, 'metastable': 2e-5})


def test_provider_refuses_a_monotone_decreasing_energy_branch(tmp_path):
    path, model = fixture_provider_table(tmp_path, decreasing=True)
    with pytest.raises(IllConditionedCoordinate, match='increase strictly'):
        make_provider(path, model)


def qualification_state(tmp_path, provider, *, composition=None):
    composition = composition or {}
    row = provider.row(1., composition)
    return {
        'u': 1.,
        'composition': composition,
        'record_path': str(tmp_path / 'eedf_qualification.json'),
        'gas_fractions': {'Ar': 1.},
        'state_populations': {'Ar(1S0)': 1.},
        'state_statistical_weights': {'Ar(1S0)': 1.},
        'energy_budget': {
            'P_abs': float(row.power_groups['field']),
            'Q_wall_electron': 0.,
            'Q_wall_ion': 0.,
            'Q_flow': 0.,
            'joule_power': float(row.power_groups['field']),
            'A6b_relative': 0.,
            'A6b_tolerance': 1e-6,
        },
    }


def direct_runner(provider, mutate=None):
    def run(spec, coordinates, fields, job_name):
        row = provider.row(float(np.log(fields[0])), coordinates).as_dict()
        direct = {}
        for key, value in row.items():
            if isinstance(value, np.ndarray):
                direct[key] = value.copy()
            elif isinstance(value, Mapping):
                direct[key] = dict(value)
            else:
                direct[key] = value
        if mutate is not None:
            mutate(direct, coordinates)
        return [direct]
    return run


def test_qualify_passes_only_when_direct_and_interpolated_rows_agree(tmp_path):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    record = provider.qualify(state, 2.5, runner=direct_runner(provider))

    assert record['status'] == 'PASS'
    assert record['terminal_state']['t'] == 2.5
    assert record['solver_row']['swarm']['mean_energy_eV'] == pytest.approx(2.)
    assert json.loads(Path(state['record_path']).read_text())['status'] == 'PASS'
    provider.require_qualified('export')


@pytest.mark.parametrize('quantity', [
    'EN_Td',
    'swarm.characteristic_energy_eV',
    'attachment_energy_eV',
    'target_fractions',
    'product_fractions',
    'power.balance',
    'f0',
])
def test_every_recorded_quantity_is_tolerance_checked(tmp_path, quantity):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    def mutate(direct, coordinates):
        if quantity.startswith('swarm.'):
            values = dict(direct['swarm'])
            values[quantity.split('.', 1)[1]] += 1.e6
            direct['swarm'] = values
        elif quantity.startswith('power.'):
            values = dict(direct['power_groups'])
            values[quantity.split('.', 1)[1]] += 1.e6
            direct['power_groups'] = values
        elif isinstance(direct[quantity], np.ndarray):
            direct[quantity][0] += 1.e6
        else:
            direct[quantity] += 1.e6

    from rmgpy.exceptions import PlasmaStateError
    rules = {
        'EN_Td': 'A2',
        'swarm.characteristic_energy_eV': 'A2',
        'attachment_energy_eV': 'A1',
        'target_fractions': 'A5',
        'product_fractions': 'A5',
        'power.balance': 'A4',
        'f0': 'A5',
    }
    check_quantities = {'f0': 'f0.weighted_L1'}
    with pytest.raises(
            PlasmaStateError,
            match='EEDF qualification failed: ' + rules[quantity]):
        provider.qualify(state, 2.5, runner=direct_runner(provider, mutate))

    record = json.loads(Path(state['record_path']).read_text())
    assert any(check['rule'] == rules[quantity] and
               check['quantity'] == check_quantities.get(quantity, quantity) and
               not check['passed'] for check in record['checks'])


def test_qualification_does_not_require_generation_tolerance(tmp_path):
    path, model = fixture_provider_table(
        tmp_path / 'table', every_quantity_tolerance=False)
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    record = provider.qualify(state, 2.5, runner=direct_runner(provider))

    assert record['status'] == 'PASS'
    assert all(check['rule'] != 'G1' for check in record['checks'])


@pytest.mark.parametrize('relative_drift, passes', [
    (1.e-4, True),
    (3.e-3, False),
])
def test_direct_transport_comparison_uses_A2_not_generation_G1(
        tmp_path, relative_drift, passes):
    path, model = fixture_provider_table(
        tmp_path / 'table', a2_rtol=2.e-3)
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    def mutate(direct, coordinates):
        swarm = dict(direct['swarm'])
        swarm['mobility_N'] *= 1. + relative_drift
        direct['swarm'] = swarm

    from rmgpy.exceptions import PlasmaStateError
    if passes:
        record = provider.qualify(
            state, 2.5, runner=direct_runner(provider, mutate))
        assert record['status'] == 'PASS'
        assert not any(check['rule'] == 'G1' for check in record['checks'])
    else:
        with pytest.raises(PlasmaStateError, match='EEDF qualification failed: A2'):
            provider.qualify(
                state, 2.5, runner=direct_runner(provider, mutate))
        record = json.loads(Path(state['record_path']).read_text())
        assert any(check['rule'] == 'A2' and not check['passed']
                   for check in record['checks'])


@pytest.mark.parametrize('identity, mutate', [
    ('operator', lambda spec: spec['solver_options'].__setitem__(
        'ionizationOperatorType', 'conservative')),
    ('binary', lambda spec: spec['binary'].__setitem__(
        'sha256', 'different-binary')),
    ('input', lambda spec: spec['input_files']['Ar'].__setitem__(
        'sha256', 'different-input')),
])
def test_generation_spec_must_match_the_loaded_manifest_before_running(
        tmp_path, identity, mutate):
    path, model = fixture_provider_table(tmp_path / 'table')
    spec_path = path / 'generation_spec.json'
    spec = json.loads(spec_path.read_text())
    mutate(spec)
    spec_path.write_text(json.dumps(spec))
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    calls = []

    def runner(*args):
        calls.append(args)
        return direct_runner(provider)(*args)

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='generation spec mismatch'):
        provider.qualify(state, 2.5, runner=runner)

    assert calls == []
    record = json.loads(Path(state['record_path']).read_text())
    assert record['status'] == 'FAIL'
    assert identity in record['error']


@pytest.mark.parametrize('mutate', [
    lambda spec: spec['state_properties']['energy'].__setitem__(
        0, 'Ar(1S0) = 123'),
    lambda spec: spec['state_properties']['statisticalWeight'].__setitem__(
        0, 'Ar(1S0) = 7'),
    lambda spec: spec['gas_properties'].__setitem__(
        'harmonicFrequency', 'masses'),
    lambda spec: spec['axes']['u'].__setitem__(1, .5),
], ids=['state-energy', 'statistical-weight', 'gas-property', 'grid'])
def test_qualification_refuses_any_unfingerprinted_physical_setup_change(
        tmp_path, mutate):
    path, model = fixture_provider_table(tmp_path / 'table')
    spec_path = path / 'generation_spec.json'
    spec = json.loads(spec_path.read_text())
    mutate(spec)
    spec_path.write_text(json.dumps(spec))
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    calls = []

    def runner(*args):
        calls.append(args)
        return direct_runner(provider)(*args)

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='physical setup fingerprint'):
        provider.qualify(state, 2.5, runner=runner)

    assert calls == []
    record = json.loads(Path(state['record_path']).read_text())
    assert record['status'] == 'FAIL'
    assert 'physical setup fingerprint' in record['error']


def test_qualification_refuses_an_artifact_without_a_setup_fingerprint(
        tmp_path):
    path, model = fixture_provider_table(
        tmp_path / 'table', setup_fingerprint=False)
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='no physical setup fingerprint'):
        provider.qualify(state, 2.5, runner=direct_runner(provider))

    record = json.loads(Path(state['record_path']).read_text())
    assert record['status'] == 'FAIL'


def test_qualification_refuses_an_artifact_without_quantity_rules(tmp_path):
    path, model = fixture_provider_table(
        tmp_path / 'table', qualification_quantity_rules=False)
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='qualification_quantity_rules'):
        provider.qualify(state, 2.5, runner=direct_runner(provider))

    record = json.loads(Path(state['record_path']).read_text())
    assert record['status'] == 'FAIL'
    assert 'qualification_quantity_rules' in record['error']


def test_qualification_refuses_an_unbudgeted_recorded_quantity(tmp_path):
    path, model = fixture_provider_table(
        tmp_path / 'table', omitted_qualification_rule='power.balance')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='no qualification rule for power.balance'):
        provider.qualify(state, 2.5, runner=direct_runner(provider))

    record = json.loads(Path(state['record_path']).read_text())
    assert record['status'] == 'FAIL'
    assert 'no qualification rule for power.balance' in record['error']


def test_qualification_refuses_a_terminal_statistical_weight_substitution(
        tmp_path):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    state['state_statistical_weights'] = {'Ar(1S0)': 7.}
    calls = []

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='physical setup fingerprint'):
        provider.qualify(
            state, 2.5, runner=lambda *args: calls.append(args))

    assert calls == []
    record = json.loads(Path(state['record_path']).read_text())
    assert record['status'] == 'FAIL'


def test_qualification_fingerprint_covers_renderer_owned_setup_fields(
        tmp_path, monkeypatch):
    import rmgpy.tools.eedf.loki as loki

    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    original = loki.setup_text

    def changed_renderer(*args, **kwargs):
        return original(*args, **kwargs).replace('isOn: true', 'isOn: false')

    monkeypatch.setattr(loki, 'setup_text', changed_renderer)
    calls = []
    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='physical setup fingerprint'):
        provider.qualify(
            state, 2.5, runner=lambda *args: calls.append(args))

    assert calls == []
    record = json.loads(Path(state['record_path']).read_text())
    assert record['status'] == 'FAIL'


def test_qualification_checks_the_actual_terminal_rendering(tmp_path, monkeypatch):
    import rmgpy.tools.eedf.loki as loki

    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    original = loki.setup_text

    def terminal_only_change(spec, coordinates, fields, folder):
        rendered = original(spec, coordinates, fields, folder)
        if fields != [1.]:
            return rendered.replace('ionizationOperatorType: usingSDCS',
                                    'ionizationOperatorType: conservative')
        return rendered

    monkeypatch.setattr(loki, 'setup_text', terminal_only_change)
    calls = []
    agreed = direct_runner(provider)

    def runner(*args):
        calls.append(args)
        return agreed(*args)

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='physical setup fingerprint'):
        provider.qualify(state, 2.5, runner=runner)

    assert calls == []
    assert json.loads(Path(state['record_path']).read_text())['status'] == 'FAIL'


@pytest.mark.parametrize('missing', ['gas_fractions', 'state_populations'])
def test_qualification_requires_every_terminal_setup_substitution(tmp_path, missing):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    del state[missing]
    calls = []
    agreed = direct_runner(provider)

    def runner(*args):
        calls.append(args)
        return agreed(*args)

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='terminal .* are incomplete'):
        provider.qualify(state, 2.5, runner=runner)

    assert calls == []
    assert json.loads(Path(state['record_path']).read_text())['status'] == 'FAIL'


def test_equivalent_terminal_weight_spellings_are_accepted(tmp_path):
    path, model = fixture_provider_table(
        tmp_path / 'table', weight_spelling='1.0')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    record = provider.qualify(state, 2.5, runner=direct_runner(provider))

    assert record['status'] == 'PASS'


def test_equivalent_terminal_only_weight_spelling_is_accepted(tmp_path, monkeypatch):
    import rmgpy.tools.eedf.loki as loki

    path, model = fixture_provider_table(tmp_path / 'table', weight_spelling='1')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    original = loki.setup_text

    def terminal_only_spelling(spec, coordinates, fields, folder):
        rendered = original(spec, coordinates, fields, folder)
        if fields != [1.]:
            return rendered.replace('Ar(1S0) = 1\n', 'Ar(1S0) = 1.0\n')
        return rendered

    monkeypatch.setattr(loki, 'setup_text', terminal_only_spelling)

    record = provider.qualify(state, 2.5, runner=direct_runner(provider))

    assert record['status'] == 'PASS'


def test_qualification_fingerprint_covers_acceptance_policy(tmp_path):
    path, model = fixture_provider_table(tmp_path / 'table')
    spec_path = path / 'generation_spec.json'
    spec = json.loads(spec_path.read_text())
    spec['floors']['rate_absolute'] = 1.
    spec_path.write_text(json.dumps(spec))
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    calls = []

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='physical setup fingerprint'):
        provider.qualify(state, 2.5, runner=lambda *args: calls.append(args))

    assert calls == []
    assert json.loads(Path(state['record_path']).read_text())['status'] == 'FAIL'


def test_direct_row_cannot_raise_its_own_A3_floor(tmp_path):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    def mutate(direct, coordinates):
        direct['k_ine'][0] = 1.e-3
        direct['rate_floors'][0] = 1.

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='EEDF qualification failed: A3'):
        provider.qualify(state, 2.5, runner=direct_runner(provider, mutate))

    assert json.loads(Path(state['record_path']).read_text())['status'] == 'FAIL'


def test_class_B_rate_material_on_interpolated_side_is_compared(tmp_path):
    path, model = fixture_provider_table(
        tmp_path / 'table', second_classification='B')
    provider = make_provider(path, model, reaction_map={7: (0, 'sup')})
    state = qualification_state(tmp_path, provider)

    def mutate(direct, coordinates):
        direct['k_ine'][1] = 0.

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='EEDF qualification failed: A3'):
        provider.qualify(state, 2.5, runner=direct_runner(provider, mutate))

    assert json.loads(Path(state['record_path']).read_text())['status'] == 'FAIL'


def test_fraction_quantity_cannot_use_a_power_budget(tmp_path):
    path, model = fixture_provider_table(
        tmp_path / 'table',
        qualification_rule_overrides={'target_fractions': 'A6'},
        a6_atol=10.)
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    def mutate(direct, coordinates):
        direct['target_fractions'][0] = 0.

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='invalid qualification rule'):
        provider.qualify(state, 2.5, runner=direct_runner(provider, mutate))

    assert json.loads(Path(state['record_path']).read_text())['status'] == 'FAIL'


@pytest.mark.parametrize('quantity, mutation', [
    ('gas_fractions', lambda row: row.pop('gas_fractions')),
    ('gas_fractions', lambda row: row['gas_fractions'].__setitem__('Ar', np.nan)),
    ('state_populations', lambda row: row.pop('state_populations')),
    ('state_populations', lambda row: row['state_populations'].__setitem__('Ar(1S0)', np.nan)),
    ('rate_floors', lambda row: row.pop('rate_floors')),
    ('rate_floors', lambda row: row['rate_floors'].__setitem__(0, np.nan)),
    ('converged', lambda row: row.pop('converged')),
    ('converged', lambda row: row.__setitem__('converged', False)),
])
def test_incomplete_or_nonfinite_solver_row_records_failure(
        tmp_path, quantity, mutation):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='invalid direct row'):
        provider.qualify(
            state, 2.5,
            runner=direct_runner(provider, lambda row, coordinates: mutation(row)))

    record = json.loads(Path(state['record_path']).read_text())
    assert record['status'] == 'FAIL'
    assert quantity in record['error']
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_time_conversion_refusal_is_recorded_and_latched(tmp_path):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    class BadTime:
        def __float__(self):
            raise EEDFError('bad terminal time')

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='bad terminal time'):
        provider.qualify(state, BadTime(), runner=direct_runner(provider))

    record = json.loads(Path(state['record_path']).read_text())
    assert record['status'] == 'FAIL'
    assert record['error'] == 'bad terminal time'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_qualify_composition_mismatch_refuses_and_invalidates_export(tmp_path):
    path, model = fixture_provider_table(tmp_path / 'table', composition=True)
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider, composition={'x': 0.25})

    def perturbed_composition_runner(spec, coordinates, fields, job_name):
        assert coordinates == {'x': 0.25}
        return [provider.row(float(np.log(fields[0])), {'x': 0.5}).as_dict()]

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='EEDF qualification failed: A2'):
        provider.qualify(state, 3., runner=perturbed_composition_runner)

    record = json.loads(Path(state['record_path']).read_text())
    assert record['status'] == 'FAIL'
    assert record['terminal_state']['composition'] == {'x': 0.25}
    assert record['solver_row']['composition'] == {'x': 0.5}
    assert any(not check['passed'] for check in record['checks'])
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_failed_qualification_is_terminal_and_cannot_requalify(tmp_path):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    def mismatch(direct, coordinates):
        direct['EN_Td'] += 1.e6

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='EEDF qualification failed: A2'):
        provider.qualify(state, 2.5, runner=direct_runner(provider, mismatch))

    calls = []

    def would_pass(*args):
        calls.append(args)
        return direct_runner(provider)(*args)

    with pytest.raises(PlasmaStateError, match='latched FAILED'):
        provider.qualify(state, 3., runner=would_pass)
    assert calls == []
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_state_resolution_refusal_is_recorded_and_latches_failure(
        tmp_path, monkeypatch):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    monkeypatch.chdir(tmp_path)
    provider.bind_qualification_state(
        lambda y, t: (_ for _ in ()).throw(EEDFError('resolver refusal')))
    calls = []

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='resolver refusal'):
        provider.qualify(np.zeros(1), 4., runner=lambda *args: calls.append(args))

    assert calls == []
    record = json.loads((tmp_path / 'eedf_qualification.json').read_text())
    assert record['status'] == 'FAIL'
    assert record['terminal_state'] == {'t': 4.}
    assert record['error'] == 'resolver refusal'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_qualify_missing_binary_is_a_recorded_failure(tmp_path):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='missing or unexecutable LoKI-B binary'):
        provider.qualify(state, 1.)

    record = json.loads(Path(state['record_path']).read_text())
    assert record['status'] == 'FAIL'
    assert 'missing or unexecutable' in record['error']


@pytest.mark.parametrize('nonfinite, serialized', [
    (float('nan'), 'NaN'),
    (float('inf'), 'Infinity'),
])
def test_nonfinite_solver_evidence_is_recorded_and_latches_failure(
        tmp_path, nonfinite, serialized):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    def mutate(direct, coordinates):
        direct['attachment_energy_eV'][0] = nonfinite

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='invalid direct row'):
        provider.qualify(state, 2.5, runner=direct_runner(provider, mutate))

    record = json.loads(Path(state['record_path']).read_text())
    assert record['status'] == 'FAIL'
    assert record['solver_row']['attachment_energy_eV'][0] == serialized
    assert record['error'] == 'invalid direct row: attachment_energy_eV'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_development_diagnostic_can_never_report_pass(tmp_path):
    path, model = fixture_provider_table(tmp_path / 'table', accepted=False)
    with development_unqualified_table_route():
        provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    record = provider.development_diagnostic(
        state, 4., runner=direct_runner(provider))

    assert record['status'] == DEVELOPMENT_DIAGNOSTIC_STATUS
    assert record['comparison_passed'] is True
    assert record['qualification_allowed'] is False
    with pytest.raises(EEDFError, match='refuses export'):
        provider.require_qualified('export')
