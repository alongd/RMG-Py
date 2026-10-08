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
from contextlib import contextmanager, nullcontext
import errno


from decimal import Decimal
import enum
from fractions import Fraction
import gc
import inspect
import json
import os
from pathlib import Path
import signal
import sys
import threading
from types import MappingProxyType, SimpleNamespace

import numpy as np
import pytest
import yaml

from rmgpy.solver.eedf import (
    AmbiguousBranch,
    EEDFError,
    EEDFTable,
    EnvelopeBreach,
    IllConditionedCoordinate,
    OutOfDomain,
)
from rmgpy.solver.eedf_provider import (
    DEVELOPMENT_DIAGNOSTIC_STATUS,
    DEVELOPMENT_UNQUALIFIED_STATUS,
    EEDFProvider,
    _development_unqualified_table_route,
)
from rmgpy.tools.eedf.artifact import _write_artifact
from rmgpy.tools.eedf.loki import qualification_setup_sha256
from rmgpy.tools.eedf.schema import (
    FingerprintMismatch,
    _content_hash,
    _file_hash,
    _interpolant_identity,
)


@contextmanager
def intended_signal_disposition(signum, handler):
    """Install one test-owned disposition and restore its Python/raw views."""
    from rmgpy.tools.eedf import signal_guard

    previous_view = signal.getsignal(signum)
    previous_raw = signal_guard._raw_sigaction(signum)
    signal.signal(signum, handler)
    try:
        yield
    finally:
        if previous_view is not None:
            signal.signal(signum, previous_view)
        if previous_raw is not None:
            signal_guard._restore_raw_sigaction(signum, previous_raw)


def fixture_provider_table(tmp_path, composition=False, branches=1, decreasing=False,
                           accepted=True, every_quantity_tolerance=True,
                           a2_rtol=1e-6, setup_fingerprint=True,
                           qualification_quantity_rules=True,
                           omitted_qualification_rule=None,
                           qualification_rule_overrides=None,
                           a6_atol=0., weight_spelling='1',
                           second_classification='A',
                           explicit_statistical_weights=True):
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
                'power_groups': {
                    'field': (1.e24 * scale) *
                             (np.exp(axes['u'][u_index]) * 1.e-21) ** 2,
                    'balance': .01 * scale},
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
        'cross_sections': {'Ar': _file_hash(ar_input)},
        'auxiliary_inputs': {'masses': _file_hash(mass_input)},
        'channel_map_sha256': 'map',
        'reactions': identities,
        'solver_options': {
            'growthModelType': 'temporal',
            'ionizationOperatorType': 'usingSDCS',
            'includeEECollisions': False,
        },
        'working_conditions': {'electronDensity': '1e16'},
        'interpolant': _interpolant_identity(),
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
        'scratch_root': str(tmp_path / 'qualification-runs'),
        'timeout_s': 1.,
    }
    if not explicit_statistical_weights:
        del setup_spec['state_properties']['statisticalWeight']
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
        'T10 interpolation qualification': {
            'check': 'T10 interpolation qualification',
            'criterion': 'fixture frozen accuracy criteria',
            'verdict_count': 1,
            'passed': accepted,
        },
        'held_out_verdicts': [{'point': {'u': .5}, 'branch_id': 'branch_0',
                               'passed': accepted, 'checks': [{'passed': accepted}]}],
        'branch_certification': {
            f'branch_{index}': 'uncertified: unseeded scans'
            for index in range(branches)
        },
        'schema_version': 1,
        'row_inputs': model,
        'fingerprint': _content_hash(model),
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
        manifest['fingerprint'] = _content_hash(model)
        manifest['qualification_setup_sha256'] = model['qualification_setup_sha256']
    path = _write_artifact(
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
        'artifact_sha256': _file_hash(path / 'table.h5'),
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
    assert pure._empirical_law_metadata['Ar2+ recombination']['evaluate_at'] == 'Te_eff'
    assert arm._manifest_snapshot['fingerprint'] == _content_hash(arm_model)


def test_unqualified_table_is_development_context_only_and_carries_blockers(tmp_path):
    qualified_path, qualified_model = fixture_provider_table(
        tmp_path / 'qualified', accepted=True)
    with _development_unqualified_table_route():
        with pytest.raises(EEDFError, match='already qualified'):
            make_provider(qualified_path, qualified_model)

    path, model = fixture_provider_table(tmp_path / 'unqualified', accepted=False)
    with pytest.raises(FingerprintMismatch, match='held-out qualification'):
        make_provider(path, model)

    with _development_unqualified_table_route():
        provider = make_provider(path, model)

    notice = provider.development_notice
    assert notice['scientific_status'] == DEVELOPMENT_UNQUALIFIED_STATUS
    assert notice['artifact_blockers'] == [
        'artifact manifest is not accepted',
        '1 of 1 held-out verdicts failed',
        'branch_0: uncertified: unseeded scans',
    ]
    assert provider._manifest_snapshot['scientific_status'] == DEVELOPMENT_UNQUALIFIED_STATUS
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
    first_identity = provider._cache_identity_snapshot
    assert provider.row(1., {'x': 0.}) is first

    changed = provider.row(1., {'x': 1.})
    assert changed is not first
    assert provider._cache_identity_snapshot != first_identity
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
        provider._reaction_mapping[4] = (0, 'ine')
    manifest = provider._manifest_snapshot
    manifest['fingerprint'] = 'mutated copy'
    assert provider._manifest_snapshot['fingerprint'] != 'mutated copy'
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
    assert provider._cache_identity_snapshot == (1., ())
    assert provider._cached_row_snapshot is reader_result[0]


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
        'scratch_root': str(tmp_path / 'qualification-runs'),
        'Tg_K': 298.15,
        'P_Pa': 666.61,
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
    def run(bundle):
        coordinates = bundle.coordinate_values
        fields = bundle.fields_Td
        row = provider.row(float(np.log(fields[0])), coordinates)._as_dict()
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
        direct['execution_identity'] = bundle.execution_identity
        return [direct]
    return run


def public_qualification_fixture(tmp_path, monkeypatch):
    """Exercise public qualification while replacing only the external solver."""
    import rmgpy.tools.eedf.loki as loki

    original_identity = qualification_setup_sha256

    def provision_driver_fields(spec):
        spec['solver_options']['numerics'] = {
            'energyGrid': {'maxEnergy': 1., 'cellNumber': 2},
            'maxPowerBalanceRelError': 1.e-3,
        }
        spec['tolerances']['normalization_atol'] = 1.e-9
        return original_identity(spec)

    monkeypatch.setattr(
        __import__(__name__, fromlist=['qualification_setup_sha256']),
        'qualification_setup_sha256', provision_driver_fields)
    path, model = fixture_provider_table(tmp_path / 'table')
    spec_path = path / 'generation_spec.json'
    spec = json.loads(spec_path.read_text())
    spec['binary']['path'] = '/usr/bin/true'
    spec['omp_threads'] = 1
    spec['scratch_root'] = str(tmp_path / 'jobs')
    spec_path.write_text(json.dumps(spec))
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    agreed = provider._json_value(provider.row(1., {})._as_dict())

    def initialize(driver, specification):
        driver.spec = specification
        driver.setups = {}

    def execute(args, **kwargs):
        output = Path(kwargs['cwd']) / 'output' / 'solve'
        output.mkdir(parents=True)
        (output / 'reducedField_1').write_text('fixture')
        return SimpleNamespace(returncode=0)

    def parsed(*args, **kwargs):
        row = copy.deepcopy(agreed)
        row['power_groups']['relative_balance_percent'] = 0.
        return row

    def enriched(row, channel_map, deck):
        result = copy.deepcopy(agreed)
        result['execution_identity'] = row['execution_identity']
        return result

    monkeypatch.setattr(loki.LoKIDriver, '__init__', initialize)
    monkeypatch.setattr(loki, '_check_solver', lambda specification: None)
    monkeypatch.setattr(loki, '_physical_properties', lambda specification: None)
    original_file_hash = loki._file_hash
    monkeypatch.setattr(
        loki, '_file_hash',
        lambda candidate: (spec['binary']['sha256']
                           if Path(candidate) == Path(spec['binary']['path'])
                           else original_file_hash(candidate)))
    monkeypatch.setattr(loki.subprocess, 'run', execute)
    monkeypatch.setattr(loki, 'parse_output', parsed)
    monkeypatch.setattr(loki, '_enrich_row', enriched)
    return provider, state


def qualify_provider(provider, y, t, *, runner=None):
    """Drive synthetic acceptance fixtures without exposing a production bypass."""
    if runner is None:
        return provider.qualify(y, t)
    from rmgpy.exceptions import PlasmaStateError
    from rmgpy.solver.eedf_provider import _qualification_attempt_locked

    @_qualification_attempt_locked
    def synthetic_attempt(provider, y, t):
        if provider._development_blockers:
            raise EEDFError(DEVELOPMENT_UNQUALIFIED_STATUS + ' refuses qualification')
        if provider._qualification_status == 'FAILED':
            raise PlasmaStateError(
                'EEDF qualification failed: provider is latched FAILED')
        return provider._run_qualification(y, t, runner, development=False)

    return synthetic_attempt(provider, y, t)


def supported_declaration_table(tmp_path, declaration):
    """Build a fixture with independently hand-audited supported declarations."""
    original, model = fixture_provider_table(tmp_path / 'original')
    original_provider = make_provider(original, model)
    spec = json.loads((original / 'generation_spec.json').read_text())
    manifest = original_provider._manifest_snapshot
    expected_gas = {'Ar': 1.}
    expected_populations = {'Ar(1S0)': 1.}
    expected_weights = {'Ar(1S0)': 1.}
    property_file = None
    if declaration == 'population-file':
        property_file = tmp_path / 'populations.txt'
        property_file.write_text('Ar(1S0) 1\n')
        entry = {
            'path': str(property_file),
            'sha256': _file_hash(property_file),
            'kind': 'property',
        }
        spec['input_files']['populations.txt'] = entry
        spec['state_properties']['population'] = ['populations.txt']
        model['auxiliary_inputs']['populations.txt'] = entry['sha256']
        manifest['solver']['input_files']['populations.txt'] = {
            'sha256': entry['sha256'], 'kind': 'property'}
    elif declaration == 'mixture':
        spec['envelopes'].update({
            'ar_fraction': {'reference': .9, 'min': .9, 'max': .9},
            'he_fraction': {'reference': .1, 'min': .1, 'max': .1},
        })
        spec['gas_properties']['fraction'] = [
            'Ar = {ar_fraction}', 'He = {he_fraction}']
        spec['state_properties']['population'].append('He(1S0) = 1')
        spec['state_properties']['energy'].append('He(1S0) = 0')
        spec['state_properties']['statisticalWeight'].append('He(1S0) = 1')
        manifest['envelopes'] = copy.deepcopy(spec['envelopes'])
        expected_gas = {'Ar': .9, 'He': .1}
        expected_populations = {'Ar(1S0)': 1., 'He(1S0)': 1.}
        expected_weights = {'Ar(1S0)': 1., 'He(1S0)': 1.}
        model['state_populations'] = copy.deepcopy(expected_populations)
    else:
        raise ValueError(declaration)
    identity = qualification_setup_sha256(spec)
    model['qualification_setup_sha256'] = identity
    manifest['row_inputs'] = model
    manifest['fingerprint'] = _content_hash(model)
    manifest['qualification_setup_sha256'] = identity
    values = np.empty((4,), dtype=object)
    for index, u in enumerate(model['axes']['u']):
        row = original_provider._json_value(original_provider.row(u, {})._as_dict())
        row['gas_fractions'] = copy.deepcopy(expected_gas)
        row['state_populations'] = copy.deepcopy(expected_populations)
        values[index] = row
    changed = _write_artifact(
        tmp_path / 'changed', manifest, {'branch_0': values}, [])
    (changed / 'generation_spec.json').write_text(json.dumps(spec))
    return (changed, model, expected_gas, expected_populations,
            expected_weights, property_file)


def test_qualify_passes_only_when_direct_and_interpolated_rows_agree(tmp_path):
    """Exercise acceptance machinery synthetically; this does not qualify a table."""
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    record = qualify_provider(provider, state, 2.5, runner=direct_runner(provider))

    assert record['status'] == 'PASS'
    assert record['terminal_state']['t'] == 2.5
    assert record['solver_row']['swarm']['mean_energy_eV'] == pytest.approx(2.)
    assert json.loads(Path(state['record_path']).read_text())['status'] == 'PASS'
    assert record['run_id']
    assert record['table_identity']
    assert record['execution_identity']
    assert record['terminal_state_identity']
    assert set(record['artifacts']) == {
        'setup.in', 'interpolated-row.json', 'solver-row.json',
        'qualification-evidence.json'}
    provider.require_qualified('export')


def test_authoritative_qualification_has_no_caller_supplied_runner(tmp_path):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    with pytest.raises(TypeError, match='runner'):
        provider.qualify(state, 2.5, runner=direct_runner(provider))


@pytest.mark.parametrize('reported_identity', [None, 'different-execution'])
def test_qualification_refuses_output_not_bound_to_prepared_execution(
        tmp_path, reported_identity):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    def unbound_runner(bundle):
        row = direct_runner(provider)(bundle)[0]
        if reported_identity is None:
            row.pop('execution_identity')
        else:
            row['execution_identity'] = reported_identity
        return [row]

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='output execution identity differs'):
        qualify_provider(provider, state, 2.5, runner=unbound_runner)

    assert json.loads(Path(state['record_path']).read_text())['status'] == 'FAIL'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


@pytest.mark.parametrize('declaration', ['population-file', 'mixture'])
def test_supported_property_declarations_reach_acceptance_and_invalid_neighbors_refuse(
        tmp_path, declaration):
    (path, model, gas, populations, weights,
     property_file) = supported_declaration_table(tmp_path, declaration)
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    state['gas_fractions'] = gas
    state['state_populations'] = populations
    state['state_statistical_weights'] = weights
    prepared = []

    def agreeing_runner(bundle):
        prepared.append(bundle)
        return direct_runner(provider)(bundle)

    record = qualify_provider(provider, state, 2.5, runner=agreeing_runner)

    assert record['status'] == 'PASS'
    provider.require_qualified('export')
    rendered = yaml.safe_load(prepared[0].setup_bytes)
    kinetics = rendered['electronKinetics']
    consumed_gas = {
        name.strip(): float(value)
        for name, value in (entry.split('=', 1)
                            for entry in kinetics['gasProperties']['fraction'])}
    consumed_populations = {
        name.strip(): float(value)
        for name, value in (entry.split('=', 1)
                            for entry in kinetics['stateProperties']['population'])}
    assert consumed_gas == gas
    assert consumed_populations == populations

    invalid_provider = make_provider(path, model)
    invalid_state = qualification_state(tmp_path / 'invalid', invalid_provider)
    invalid_state['gas_fractions'] = copy.deepcopy(gas)
    invalid_state['state_populations'] = copy.deepcopy(populations)
    invalid_state['state_statistical_weights'] = copy.deepcopy(weights)
    if declaration == 'population-file':
        property_file.write_text('Ar(1S0) 0.5\n')
        expected_error = 'input file populations.txt'
    else:
        invalid_state['gas_fractions'] = {'Ar': .8, 'He': .1}
        expected_error = 'normalize gas fractions to one'
    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match=expected_error):
        qualify_provider(invalid_provider,
            invalid_state, 2.5, runner=direct_runner(invalid_provider))
    with pytest.raises(EEDFError, match='FAILED.*export'):
        invalid_provider.require_qualified('export')


def test_omitted_optional_statistical_weight_declaration_reaches_acceptance(tmp_path):
    path, model = fixture_provider_table(
        tmp_path / 'table', explicit_statistical_weights=False)
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    record = qualify_provider(provider, state, 2.5, runner=direct_runner(provider))

    assert record['status'] == 'PASS'
    provider.require_qualified('export')


def test_qualification_renders_once_and_executes_the_tested_operator(
        tmp_path, monkeypatch):
    """The real qualification driver must execute its one prepared rendering."""
    import rmgpy.tools.eedf.loki as loki

    original_identity = qualification_setup_sha256

    def provision_driver_fields(spec):
        spec['solver_options']['numerics'] = {
            'energyGrid': {'maxEnergy': 1., 'cellNumber': 2},
            'maxPowerBalanceRelError': 1.e-3,
        }
        spec['tolerances']['normalization_atol'] = 1.e-9
        return original_identity(spec)

    monkeypatch.setattr(
        __import__(__name__, fromlist=['qualification_setup_sha256']),
        'qualification_setup_sha256', provision_driver_fields)
    path, model = fixture_provider_table(tmp_path / 'table')
    spec_path = path / 'generation_spec.json'
    spec = json.loads(spec_path.read_text())
    spec['binary']['path'] = '/usr/bin/true'
    spec['omp_threads'] = 1
    spec['scratch_root'] = str(tmp_path / 'jobs')
    spec_path.write_text(json.dumps(spec))
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    agreed = provider._json_value(provider.row(1., {})._as_dict())
    original_render = loki._setup_text
    render_folders = []
    executed_setups = []

    def renderer(specification, coordinates, fields, folder):
        rendered = original_render(specification, coordinates, fields, folder)
        render_folders.append(folder)
        if len(render_folders) > 2:
            rendered = rendered.replace(
                'ionizationOperatorType: usingSDCS',
                'ionizationOperatorType: conservative')
        return rendered

    def initialize(driver, specification):
        driver.spec = specification
        driver.setups = {}

    def execute(args, **kwargs):
        work = Path(kwargs['cwd'])
        executed_setups.append((work / 'setup.in').read_bytes())
        (work / 'output' / 'solve' / 'reducedField_1').mkdir(parents=True)
        return SimpleNamespace(returncode=0)

    def parsed(*args, **kwargs):
        row = copy.deepcopy(agreed)
        row['power_groups']['relative_balance_percent'] = 0.
        return row

    monkeypatch.setattr(loki, '_setup_text', renderer)
    monkeypatch.setattr(loki.LoKIDriver, '__init__', initialize)
    monkeypatch.setattr(loki, '_check_solver', lambda specification: None)
    monkeypatch.setattr(loki, '_physical_properties', lambda specification: None)
    original_file_hash = loki._file_hash
    monkeypatch.setattr(
        loki, '_file_hash',
        lambda candidate: (spec['binary']['sha256']
                           if Path(candidate) == Path(spec['binary']['path'])
                           else original_file_hash(candidate)))
    monkeypatch.setattr(loki.subprocess, 'run', execute)
    monkeypatch.setattr(loki, 'parse_output', parsed)
    def enrich_with_bound_identity(row, channel_map, deck):
        enriched = copy.deepcopy(agreed)
        enriched['execution_identity'] = row['execution_identity']
        return enriched

    monkeypatch.setattr(loki, '_enrich_row', enrich_with_bound_identity)

    record = qualify_provider(provider, state, 2.5)

    assert record['status'] == 'PASS'
    provider.require_qualified('export')
    assert render_folders == ['solve']
    assert len(executed_setups) == 1
    assert b'ionizationOperatorType: usingSDCS' in executed_setups[0]
    assert b'ionizationOperatorType: conservative' not in executed_setups[0]


def test_record_write_refusal_never_leaves_an_accepted_pass(
        tmp_path, monkeypatch, record_property):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    record_path = Path(state['record_path'])
    original_write = provider._stage_record
    writes = []

    def refuse_first_record_write(target, *args, **kwargs):
        is_publication = target.destination == record_path
        if is_publication:
            writes.append(target)
            if len(writes) == 1:
                raise OSError('qualification record publication refused')
        return original_write(target, *args, **kwargs)

    monkeypatch.setattr(provider, '_stage_record', refuse_first_record_write)
    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='publication refused'):
        qualify_provider(provider, state, 2.5, runner=direct_runner(provider))
    assert writes, 'record staging fault was not injected'

    persisted = [json.loads(candidate.read_text())
                 for candidate in tmp_path.rglob('*.json')]
    assert persisted
    assert all(record.get('status') != 'PASS' for record in persisted)
    assert provider._qualification_record.get('status') != 'PASS'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')

    record_property('hook_name', 'private-record-stage')
    record_property('hook_hits', int(bool(writes)))


def test_operator_mutation_after_preparation_cannot_change_execution(
        tmp_path):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    spec_path = path / 'generation_spec.json'
    executed = []

    def mutate_after_preparation(bundle):
        spec = json.loads(spec_path.read_text())
        spec['solver_options']['ionizationOperatorType'] = 'conservative'
        spec_path.write_text(json.dumps(spec))
        executed.append(bundle.setup_bytes)
        return direct_runner(provider)(bundle)

    record = qualify_provider(provider, state, 2.5, runner=mutate_after_preparation)

    assert record['status'] == 'PASS'
    assert b'ionizationOperatorType: usingSDCS' in executed[0]
    assert b'ionizationOperatorType: conservative' not in executed[0]
    provider.require_qualified('export')


def test_rendered_terminal_fraction_mismatch_refuses_before_acceptance(
        tmp_path, monkeypatch):
    import rmgpy.tools.eedf.loki as loki

    (path, model, gas, populations, weights,
     _) = supported_declaration_table(tmp_path / 'table', 'mixture')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    state['gas_fractions'] = gas
    state['state_populations'] = populations
    state['state_statistical_weights'] = weights
    original_render = loki._setup_text

    def wrong_fraction(spec, coordinates, fields, folder):
        rendered = original_render(spec, coordinates, fields, folder)
        return rendered.replace(
            'Ar = 0.9', 'Ar = 0.8', 1
        ).replace(
            'He = 0.1', 'He = 0.2', 1)

    monkeypatch.setattr(loki, '_setup_text', wrong_fraction)
    calls = []
    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='gas_fractions.Ar'):
        qualify_provider(provider, state, 2.5, runner=lambda bundle: calls.append(bundle))

    assert calls == []
    assert json.loads(provider.qualification_record_path.read_text())['status'] == 'FAIL'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_blank_rendered_terminal_value_names_the_refused_field(
        tmp_path, monkeypatch):
    import rmgpy.tools.eedf.loki as loki

    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    original_render = loki._setup_text

    def blank_fraction(spec, coordinates, fields, folder):
        return original_render(spec, coordinates, fields, folder).replace(
            '    fraction:\n      - Ar = 1\n', '    fraction: []\n')

    monkeypatch.setattr(loki, '_setup_text', blank_fraction)
    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='gasProperties.fraction'):
        qualify_provider(provider, state, 2.5, runner=direct_runner(provider))
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_output_write_failure_has_no_current_acceptance_manifest(
        tmp_path, monkeypatch):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    original_write = Path.write_text

    def refuse_solver_row(target, *args, **kwargs):
        if target.name == 'solver-row.json':
            raise OSError('solver-row output refused')
        return original_write(target, *args, **kwargs)

    monkeypatch.setattr(Path, 'write_text', refuse_solver_row)
    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='output refused'):
        qualify_provider(provider, state, 2.5, runner=direct_runner(provider))

    assert not Path(state['record_path']).is_file() or json.loads(
        Path(state['record_path']).read_text())['status'] != 'PASS'
    assert provider._acceptance_manifest_path is None
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_old_pass_cannot_satisfy_a_new_failed_attempt(tmp_path):
    old_path = tmp_path / 'old-acceptance.json'
    old_path.write_text(json.dumps({'status': 'PASS', 'run_id': 'old-run'}))
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path / 'new', provider)

    def mismatch(row, coordinates):
        row['EN_Td'] += 1.e6

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='EEDF qualification failed: A2'):
        qualify_provider(provider, state, 2.5, runner=direct_runner(provider, mismatch))

    assert json.loads(old_path.read_text())['run_id'] == 'old-run'
    assert json.loads(provider.qualification_record_path.read_text())['status'] == 'FAIL'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_acceptance_read_from_the_wrong_output_location_refuses(tmp_path):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    record = qualify_provider(provider, state, 2.5, runner=direct_runner(provider))
    wrong = tmp_path / 'copied-acceptance.json'
    wrong.write_text(Path(record['acceptance_manifest']).read_text())

    with pytest.raises(EEDFError, match='wrong output location'):
        provider._validate_acceptance_manifest(wrong, record['run_id'])
    provider.require_qualified('export')


@pytest.mark.parametrize('tamper', ['run-id', 'artifact'])
def test_export_revalidates_current_manifest_and_bound_artifacts(tmp_path, tamper):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    record = qualify_provider(provider, state, 2.5, runner=direct_runner(provider))

    if tamper == 'run-id':
        manifest_path = Path(record['acceptance_manifest'])
        manifest = json.loads(manifest_path.read_text())
        manifest['run_id'] = 'stale-run'
        manifest_path.write_text(json.dumps(manifest))
        expected = 'manifest bytes differ'
    else:
        artifact = Path(record['run_directory']) / 'solver-row.json'
        artifact.write_text('{}')
        expected = 'artifact bytes differ'

    with pytest.raises(EEDFError, match=expected):
        provider.require_qualified('export')
    assert provider._qualification_status == 'FAILED'


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
        qualify_provider(provider, state, 2.5, runner=direct_runner(provider, mutate))

    record = json.loads(Path(state['record_path']).read_text())
    assert any(check['rule'] == rules[quantity] and
               check['quantity'] == check_quantities.get(quantity, quantity) and
               not check['passed'] for check in record['checks'])


def test_qualification_does_not_require_generation_tolerance(tmp_path):
    path, model = fixture_provider_table(
        tmp_path / 'table', every_quantity_tolerance=False)
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    record = qualify_provider(provider, state, 2.5, runner=direct_runner(provider))

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
        record = qualify_provider(provider,
            state, 2.5, runner=direct_runner(provider, mutate))
        assert record['status'] == 'PASS'
        assert not any(check['rule'] == 'G1' for check in record['checks'])
    else:
        with pytest.raises(PlasmaStateError, match='EEDF qualification failed: A2'):
            qualify_provider(provider,
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
        qualify_provider(provider, state, 2.5, runner=runner)

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
        qualify_provider(provider, state, 2.5, runner=runner)

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
        qualify_provider(provider, state, 2.5, runner=direct_runner(provider))

    record = json.loads(Path(state['record_path']).read_text())
    assert record['status'] == 'FAIL'


def test_qualification_refuses_an_artifact_without_quantity_rules(tmp_path):
    path, model = fixture_provider_table(
        tmp_path / 'table', qualification_quantity_rules=False)
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='qualification_quantity_rules'):
        qualify_provider(provider, state, 2.5, runner=direct_runner(provider))

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
        qualify_provider(provider, state, 2.5, runner=direct_runner(provider))

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
        qualify_provider(provider,
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
    original = loki._setup_text

    def changed_renderer(*args, **kwargs):
        return original(*args, **kwargs).replace('isOn: true', 'isOn: false')

    monkeypatch.setattr(loki, '_setup_text', changed_renderer)
    calls = []
    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='physical setup fingerprint'):
        qualify_provider(provider,
            state, 2.5, runner=lambda *args: calls.append(args))

    assert calls == []
    record = json.loads(Path(state['record_path']).read_text())
    assert record['status'] == 'FAIL'


def test_qualification_checks_the_actual_terminal_rendering(tmp_path, monkeypatch):
    import rmgpy.tools.eedf.loki as loki

    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    original = loki._setup_text

    def terminal_only_change(spec, coordinates, fields, folder):
        rendered = original(spec, coordinates, fields, folder)
        if fields != [1.]:
            return rendered.replace('ionizationOperatorType: usingSDCS',
                                    'ionizationOperatorType: conservative')
        return rendered

    monkeypatch.setattr(loki, '_setup_text', terminal_only_change)
    calls = []
    agreed = direct_runner(provider)

    def runner(*args):
        calls.append(args)
        return agreed(*args)

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='physical setup fingerprint'):
        qualify_provider(provider, state, 2.5, runner=runner)

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
        qualify_provider(provider, state, 2.5, runner=runner)

    assert calls == []
    assert json.loads(Path(state['record_path']).read_text())['status'] == 'FAIL'


def test_equivalent_terminal_weight_spellings_are_accepted(tmp_path):
    path, model = fixture_provider_table(
        tmp_path / 'table', weight_spelling='1.0')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    record = qualify_provider(provider, state, 2.5, runner=direct_runner(provider))

    assert record['status'] == 'PASS'


def test_equivalent_terminal_only_weight_spelling_is_accepted(tmp_path, monkeypatch):
    import rmgpy.tools.eedf.loki as loki

    path, model = fixture_provider_table(tmp_path / 'table', weight_spelling='1')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    original = loki._setup_text

    def terminal_only_spelling(spec, coordinates, fields, folder):
        rendered = original(spec, coordinates, fields, folder)
        if fields != [1.]:
            return rendered.replace('Ar(1S0) = 1\n', 'Ar(1S0) = 1.0\n')
        return rendered

    monkeypatch.setattr(loki, '_setup_text', terminal_only_spelling)

    record = qualify_provider(provider, state, 2.5, runner=direct_runner(provider))

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
        qualify_provider(provider, state, 2.5, runner=lambda *args: calls.append(args))

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
        qualify_provider(provider, state, 2.5, runner=direct_runner(provider, mutate))

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
        qualify_provider(provider, state, 2.5, runner=direct_runner(provider, mutate))

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
        qualify_provider(provider, state, 2.5, runner=direct_runner(provider, mutate))

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
        qualify_provider(provider,
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
        qualify_provider(provider, state, BadTime(), runner=direct_runner(provider))

    record = json.loads(Path(state['record_path']).read_text())
    assert record['status'] == 'FAIL'
    assert record['error'] == 'bad terminal time'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_qualify_composition_mismatch_refuses_and_invalidates_export(tmp_path):
    path, model = fixture_provider_table(tmp_path / 'table', composition=True)
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider, composition={'x': 0.25})

    def perturbed_composition_runner(bundle):
        coordinates = bundle.coordinate_values
        fields = bundle.fields_Td
        assert coordinates == {'x': 0.25}
        row = provider.row(float(np.log(fields[0])), {'x': 0.5})._as_dict()
        row['execution_identity'] = bundle.execution_identity
        return [row]

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='EEDF qualification failed: A2'):
        qualify_provider(provider, state, 3., runner=perturbed_composition_runner)

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
        qualify_provider(provider, state, 2.5, runner=direct_runner(provider, mismatch))

    calls = []

    def would_pass(*args):
        calls.append(args)
        return direct_runner(provider)(*args)

    with pytest.raises(PlasmaStateError, match='latched FAILED'):
        qualify_provider(provider, state, 3., runner=would_pass)
    assert calls == []
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_state_resolution_refusal_is_recorded_and_latches_failure(
        tmp_path, monkeypatch):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    provider.bind_qualification_state(
        lambda y, t: (_ for _ in ()).throw(EEDFError('resolver refusal')))
    calls = []

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='resolver refusal'):
        qualify_provider(provider, np.zeros(1), 4., runner=lambda *args: calls.append(args))

    assert calls == []
    record = json.loads(provider.qualification_record_path.read_text())
    assert record['status'] == 'FAIL'
    assert record['terminal_state']['t'] == 4.
    assert record['terminal_state']['scratch_root'].endswith(
        '/table/qualification-runs')
    assert record['error'] == 'resolver refusal'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_qualify_missing_binary_is_a_recorded_failure(tmp_path):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='missing or unexecutable LoKI-B binary'):
        qualify_provider(provider, state, 1.)

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
        qualify_provider(provider, state, 2.5, runner=direct_runner(provider, mutate))

    record = json.loads(Path(state['record_path']).read_text())
    assert record['status'] == 'FAIL'
    assert record['solver_row']['attachment_energy_eV'][0] == serialized
    assert record['error'] == 'invalid direct row: attachment_energy_eV'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_development_diagnostic_can_never_report_pass(tmp_path):
    path, model = fixture_provider_table(tmp_path / 'table', accepted=False)
    with _development_unqualified_table_route():
        provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    record = provider.development_diagnostic(
        state, 4., runner=direct_runner(provider))

    assert record['status'] == DEVELOPMENT_DIAGNOSTIC_STATUS
    assert record['comparison_passed'] is True
    assert record['qualification_allowed'] is False
    with pytest.raises(EEDFError, match='refuses export'):
        provider.require_qualified('export')


def test_resolver_state_is_snapshotted_before_the_direct_run(tmp_path):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    resolved = qualification_state(tmp_path, provider)
    expected = copy.deepcopy(resolved)
    provider.bind_qualification_state(lambda y, t: resolved)

    def mutate_resolver_result(bundle):
        resolved['Tg_K'] = 900.
        resolved['state_populations']['Ar(1S0)'] = .5
        return direct_runner(provider)(bundle)

    record = qualify_provider(
        provider, np.zeros(1), 2.5, runner=mutate_resolver_result)

    assert record['status'] == 'PASS'
    assert record['terminal_state']['Tg_K'] == expected['Tg_K']
    assert record['terminal_state']['state_populations'] == expected['state_populations']
    provider.require_qualified('export')


@pytest.mark.parametrize('tamper', ['omit', 'retarget'])
def test_export_uses_committed_artifact_bytes_not_manifest_claims(tmp_path, tamper):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    record = qualify_provider(provider, state, 2.5, runner=direct_runner(provider))
    manifest_path = Path(record['acceptance_manifest'])
    manifest = json.loads(manifest_path.read_text())
    solver_row = Path(record['run_directory']) / 'solver-row.json'
    solver_row.write_text('{}\n')
    if tamper == 'omit':
        del manifest['artifacts']['solver-row.json']
    else:
        manifest['artifacts']['solver-row.json'] = _file_hash(solver_row)
    manifest_path.write_text(json.dumps(manifest))

    with pytest.raises(EEDFError, match='artifact bytes differ'):
        provider.require_qualified('export')


def test_publication_refuses_setup_bytes_changed_after_execution(
        tmp_path, monkeypatch):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)
    original_write = provider._write_json

    def tamper_before_first_row(path, value):
        path = Path(path)
        if path.name == 'interpolated-row.json':
            setup = path.parent / 'setup.in'
            setup.write_bytes(setup.read_bytes() + b'# changed after execution\n')
        return original_write(path, value)

    monkeypatch.setattr(provider, '_write_json', tamper_before_first_row)
    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='setup.in artifact bytes differ'):
        qualify_provider(provider, state, 2.5, runner=direct_runner(provider))

    assert not Path(state['record_path']).is_file() or json.loads(
        Path(state['record_path']).read_text())['status'] != 'PASS'


def test_publication_validates_before_the_atomic_pass_replace(
        tmp_path, monkeypatch):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    state = qualification_state(tmp_path, provider)

    def validation_failure(*args, **kwargs):
        raise OSError('acceptance validation failed')

    def cleanup_failure(*args, **kwargs):
        raise OSError('failure publication failed')

    monkeypatch.setattr(provider, '_validate_acceptance_manifest', validation_failure)
    monkeypatch.setattr(provider, '_write_failure', cleanup_failure)
    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='failure publication failed'):
        qualify_provider(provider, state, 2.5, runner=direct_runner(provider))

    assert not Path(state['record_path']).is_file() or json.loads(
        Path(state['record_path']).read_text())['status'] != 'PASS'


def test_pass_publication_has_no_fallible_step_after_atomic_replace(
        tmp_path, monkeypatch, record_property):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    record_directory = tmp_path / 'record-directory'
    record_directory.mkdir()
    record_path = record_directory / 'acceptance.json'
    state['record_path'] = str(record_path)
    original_replace = os.link
    hook_hits = []

    def replace_then_remove_access(source, destination, *args, **kwargs):
        result = original_replace(source, destination, *args, **kwargs)
        hook_hits.append('after-public-link')
        record_directory.chmod(0)
        return result

    monkeypatch.setattr(os, 'link', replace_then_remove_access)
    try:
        record = provider.qualify(state, 2.5)
    finally:
        record_directory.chmod(0o700)

    assert hook_hits
    assert record['status'] == 'PASS'
    assert json.loads(record_path.read_text())['status'] == 'PASS'
    provider.require_qualified('export')

    record_property('hook_name', 'after-public-link/permissions')
    record_property('hook_hits', len(hook_hits))


@pytest.mark.parametrize('failure_source', ['qualification', 'terminal'])
def test_failure_atomically_replaces_a_read_only_prior_pass(
        tmp_path, monkeypatch, failure_source):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    record = provider.qualify(state, 2.5)
    record_path = Path(record['acceptance_manifest'])
    record_path.chmod(0o400)

    try:
        if failure_source == 'qualification':
            invalid_state = copy.deepcopy(state)
            invalid_state['gas_fractions']['Ar'] = .5
            from rmgpy.exceptions import PlasmaStateError
            with pytest.raises(PlasmaStateError):
                provider.qualify(invalid_state, 2.5)
        else:
            provider.record_terminal_failure(
                state, 2.5, EEDFError('terminal qualification failed'))
    finally:
        if record_path.exists():
            record_path.chmod(0o600)

    assert json.loads(record_path.read_text())['status'] == 'FAIL'
    assert provider._qualification_status == 'FAILED'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_failure_atomically_replaces_a_read_only_pass_from_an_earlier_process(
        tmp_path, monkeypatch):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    record_path = Path(state['record_path'])
    record_path.write_text(json.dumps({'status': 'PASS', 'run_id': 'earlier-process'}))
    record_path.chmod(0o400)
    state['gas_fractions']['Ar'] = .5

    from rmgpy.exceptions import PlasmaStateError
    try:
        with pytest.raises(PlasmaStateError):
            provider.qualify(state, 2.5)
    finally:
        if record_path.exists():
            record_path.chmod(0o600)

    assert json.loads(record_path.read_text())['status'] == 'FAIL'
    assert provider._qualification_status == 'FAILED'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


@pytest.mark.parametrize('caller', ['qualification', 'terminal'])
@pytest.mark.parametrize('publication_failure', ['diagnostic-enospc', 'rename-eio'])
def test_fresh_provider_invalidates_prior_pass_before_failure_publication(
        tmp_path, monkeypatch, caller, publication_failure, record_property):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    provider.qualify(state, 2.5)
    record_path = Path(state['record_path'])
    fresh = make_provider(
        provider._table.path,
        provider._json_value(provider.runtime_metadata['row_inputs']))
    invalid_state = copy.deepcopy(state)
    invalid_state['gas_fractions']['Ar'] = .5
    original_write = fresh._stage_record
    original_replace = os.link
    hook_hits = []

    def fail_diagnostic(target, *args, **kwargs):
        if (publication_failure == 'diagnostic-enospc'
                and target.destination.name.endswith('.failure.json')):
            hook_hits.append('diagnostic-stage')
            raise OSError(errno.ENOSPC, 'failure diagnostic disk is full')
        return original_write(target, *args, **kwargs)

    def fail_failure_rename(source, destination, *args, **kwargs):
        if (publication_failure == 'rename-eio'
                and str(source).endswith('.failure.tmp')):
            hook_hits.append('failure-link')
            raise OSError(errno.EIO, 'failure record rename failed')
        return original_replace(source, destination, *args, **kwargs)

    monkeypatch.setattr(fresh, '_stage_record', fail_diagnostic)
    monkeypatch.setattr(os, 'link', fail_failure_rename)
    if caller == 'qualification':
        from rmgpy.exceptions import PlasmaStateError
        with pytest.raises(PlasmaStateError):
            fresh.qualify(invalid_state, 2.5)
    else:
        fresh.record_terminal_failure(
            invalid_state, 2.5, EEDFError('terminal qualification failed'))

    assert not record_path.exists() or json.loads(record_path.read_text())['status'] != 'PASS'
    assert fresh._qualification_status == 'FAILED'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        fresh.require_qualified('export')

    if publication_failure != 'incomplete-state':
        assert hook_hits, 'publication fault hook was not hit'

    record_property('hook_name', 'fresh-record/'+publication_failure)
    record_property('hook_hits', len(hook_hits))


@pytest.mark.parametrize('as_resolver', [False, True])
def test_terminal_failure_without_output_paths_replaces_known_pass(
        tmp_path, monkeypatch, as_resolver):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    provider.qualify(state, 2.5)
    record_path = Path(state['record_path'])
    default_state = copy.deepcopy(state)
    del default_state['record_path']
    del default_state['scratch_root']
    argument = default_state
    if as_resolver:
        provider.bind_qualification_state(lambda y, t: default_state)
        argument = np.array([1.])

    record = provider.record_terminal_failure(
        argument, 2.5, EEDFError('terminal power gate failed'))

    assert record['status'] == 'FAIL'
    assert json.loads(record_path.read_text())['status'] == 'FAIL'
    assert provider._qualification_status == 'FAILED'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


@pytest.mark.parametrize('caller', ['qualification', 'terminal'])
def test_wrong_typed_record_path_replaces_known_pass(
        tmp_path, monkeypatch, caller):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    provider.qualify(state, 2.5)
    record_path = Path(state['record_path'])
    invalid_state = copy.deepcopy(state)
    invalid_state['record_path'] = 17

    if caller == 'qualification':
        from rmgpy.exceptions import PlasmaStateError
        with pytest.raises(PlasmaStateError):
            provider.qualify(invalid_state, 2.5)
    else:
        provider.record_terminal_failure(
            invalid_state, 2.5, EEDFError('terminal power gate failed'))

    assert json.loads(record_path.read_text())['status'] == 'FAIL'
    assert provider._qualification_status == 'FAILED'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


@pytest.mark.parametrize('caller', ['qualification', 'terminal'])
@pytest.mark.parametrize(
    'failure', ['diagnostic-enospc', 'failure-rename-eio', 'incomplete-state'])
def test_fresh_resolver_path_is_invalidated_before_failure_finalization(
        tmp_path, monkeypatch, caller, failure, record_property):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    provider.qualify(state, 2.5)
    record_path = Path(state['record_path'])
    fresh = make_provider(
        provider._table.path,
        provider._json_value(provider.runtime_metadata['row_inputs']))
    invalid_state = copy.deepcopy(state)
    invalid_state['gas_fractions']['Ar'] = .5
    if failure == 'incomplete-state':
        del invalid_state['Tg_K']
    fresh.bind_qualification_state(lambda y, t: invalid_state)
    original_write = fresh._stage_record
    original_replace = os.link
    hook_hits = []

    def fail_diagnostic(target, *args, **kwargs):
        if (failure == 'diagnostic-enospc'
                and target.destination.name.endswith('.failure.json')):
            hook_hits.append('diagnostic-stage')
            raise OSError(errno.ENOSPC, 'failure diagnostic disk is full')
        return original_write(target, *args, **kwargs)

    def fail_failure_rename(source, destination, *args, **kwargs):
        if (failure == 'failure-rename-eio'
                and str(source).endswith('.failure.tmp')):
            hook_hits.append('failure-link')
            raise OSError(errno.EIO, 'failure record rename failed')
        return original_replace(source, destination, *args, **kwargs)

    monkeypatch.setattr(fresh, '_stage_record', fail_diagnostic)
    monkeypatch.setattr(os, 'link', fail_failure_rename)
    if caller == 'qualification':
        from rmgpy.exceptions import PlasmaStateError
        with pytest.raises(PlasmaStateError):
            fresh.qualify(np.array([1.]), 2.5)
    else:
        fresh.record_terminal_failure(
            np.array([1.]), 2.5, EEDFError('terminal power gate failed'))

    assert not record_path.exists() or json.loads(
        record_path.read_text())['status'] != 'PASS'
    assert fresh._qualification_status == 'FAILED'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        fresh.require_qualified('export')

    if failure != 'incomplete-state':
        assert hook_hits, 'publication fault hook was not hit'

    record_property('hook_name', 'resolver-record/'+failure)
    record_property('hook_hits', len(hook_hits))


@pytest.mark.parametrize('caller', ['qualification', 'terminal'])
def test_writable_prior_record_is_invalidated_when_parent_forbids_rename(
        tmp_path, monkeypatch, caller):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    provider.qualify(state, 2.5)
    record_path = Path(state['record_path'])
    record_directory = record_path.parent
    invalid_state = copy.deepcopy(state)
    invalid_state['gas_fractions']['Ar'] = .5
    record_directory.chmod(0o500)

    try:
        if caller == 'qualification':
            from rmgpy.exceptions import PlasmaStateError
            with pytest.raises(PlasmaStateError):
                provider.qualify(invalid_state, 2.5)
        else:
            provider.record_terminal_failure(
                invalid_state, 2.5, EEDFError('terminal qualification failed'))
    finally:
        record_directory.chmod(0o700)

    assert json.loads(record_path.read_text())['status'] != 'PASS'
    assert provider._qualification_status == 'FAILED'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


@pytest.mark.parametrize('failure_mode', ['read-only-directory', 'enospc'])
def test_failed_retry_withdraws_prior_pass_before_failure_publication(
        tmp_path, monkeypatch, failure_mode):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    record = provider.qualify(state, 2.5)
    record_path = Path(record['acceptance_manifest'])
    record_directory = record_path.parent
    original_write = Path.write_text

    def fail_failure_diagnostic(target, *args, **kwargs):
        if target.name.endswith('.failure.json'):
            if failure_mode == 'read-only-directory':
                record_directory.chmod(0o500)
                return original_write(target, *args, **kwargs)
            raise OSError(errno.ENOSPC, 'qualification failure disk is full')
        return original_write(target, *args, **kwargs)

    invalid_state = copy.deepcopy(state)
    invalid_state['gas_fractions']['Ar'] = .5
    monkeypatch.setattr(Path, 'write_text', fail_failure_diagnostic)
    from rmgpy.exceptions import PlasmaStateError
    try:
        with pytest.raises(PlasmaStateError):
            provider.qualify(invalid_state, 2.5)
    finally:
        record_directory.chmod(0o700)

    assert not record_path.exists() or json.loads(record_path.read_text())['status'] != 'PASS'
    assert provider._qualification_status == 'FAILED'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_already_unwritable_directory_invalidates_writable_prior_record(
        tmp_path, monkeypatch):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    record = provider.qualify(state, 2.5)
    record_path = Path(record['acceptance_manifest'])
    record_directory = record_path.parent
    invalid_state = copy.deepcopy(state)
    invalid_state['gas_fractions']['Ar'] = .5
    record_directory.chmod(0o500)

    from rmgpy.exceptions import PlasmaStateError
    try:
        with pytest.raises(PlasmaStateError):
            provider.qualify(invalid_state, 2.5)
    finally:
        record_directory.chmod(0o700)

    assert json.loads(record_path.read_text())['status'] != 'PASS'
    assert provider._qualification_status == 'FAILED'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


@pytest.mark.parametrize('rename_failure', ['record-path-directory', 'eio'])
def test_failed_pass_rename_removes_the_pass_labelled_temporary(
        tmp_path, monkeypatch, rename_failure, record_property):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    record_directory = tmp_path / 'acceptance-records'
    record_directory.mkdir()
    record_path = record_directory / 'acceptance.json'
    state['record_path'] = str(record_path)
    original_replace = os.link

    if rename_failure == 'record-path-directory':
        record_path.mkdir()
    else:
        failed = False

        def fail_first_rename(source, destination, *args, **kwargs):
            nonlocal failed
            if not failed:
                failed = True
                raise OSError(errno.EIO, 'acceptance rename failed')
            return original_replace(source, destination, *args, **kwargs)

        monkeypatch.setattr(os, 'link', fail_first_rename)

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError):
        provider.qualify(state, 2.5)

    if rename_failure == 'eio':
        assert failed, 'publication link hook was not hit'
    pass_temporaries = []
    for candidate in record_directory.iterdir():
        if not candidate.is_file():
            continue
        try:
            if json.loads(candidate.read_text()).get('status') == 'PASS':
                pass_temporaries.append(candidate)
        except (OSError, ValueError):
            pass
    assert pass_temporaries == []
    assert provider._qualification_status == 'FAILED'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')

    record_property('hook_name', 'PASS-link/'+rename_failure)
    record_property('hook_hits', int(failed) if rename_failure == 'eio' else 0)


@pytest.mark.parametrize('boundary', ['after-validation', 'final-rename'])
def test_keyboard_interrupt_before_pass_commit_cleans_and_latches_failure(
        tmp_path, monkeypatch, boundary, record_property):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    original_validate = provider._validate_acceptance_manifest
    original_replace = os.link
    hook_hits = []

    def interrupt_after_validation(*args, **kwargs):
        result = original_validate(*args, **kwargs)
        if boundary == 'after-validation' and kwargs.get('staged'):
            hook_hits.append(boundary)
            raise KeyboardInterrupt()
        return result

    def interrupt_final_rename(source, destination, *args, **kwargs):
        if (boundary == 'final-rename'
                and str(source).endswith('.tmp')
                and not str(source).endswith('.failure.tmp')):
            hook_hits.append(boundary)
            raise KeyboardInterrupt()
        return original_replace(source, destination, *args, **kwargs)

    monkeypatch.setattr(
        provider, '_validate_acceptance_manifest', interrupt_after_validation)
    monkeypatch.setattr(os, 'link', interrupt_final_rename)
    with pytest.raises(KeyboardInterrupt):
        provider.qualify(state, 2.5)

    assert hook_hits, 'interrupt hook was not hit'
    pass_files = []
    for candidate in tmp_path.rglob('*'):
        if not candidate.is_file() or candidate.suffix not in ('.json', '.tmp'):
            continue
        try:
            if json.loads(candidate.read_text()).get('status') == 'PASS':
                pass_files.append(candidate)
        except (OSError, ValueError):
            pass
    assert pass_files == []
    assert provider._qualification_status == 'FAILED'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')

    record_property('hook_name', 'interrupt/'+boundary)
    record_property('hook_hits', len(hook_hits))


def test_non_exception_abort_before_pass_commit_cleans_and_reraises(
        tmp_path, monkeypatch):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)

    class AttemptAbort(BaseException):
        pass

    def abort_validation(*args, **kwargs):
        raise AttemptAbort('stop before commit')

    monkeypatch.setattr(provider, '_validate_acceptance_manifest', abort_validation)
    with pytest.raises(AttemptAbort, match='stop before commit'):
        provider.qualify(state, 2.5)

    assert not any(
        json.loads(candidate.read_text()).get('status') == 'PASS'
        for candidate in tmp_path.rglob('*.tmp'))
    assert provider._qualification_status == 'FAILED'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_queued_attempt_after_interruption_sees_failed_latch(
        tmp_path, monkeypatch):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    original_validate = provider._validate_acceptance_manifest
    predecessor_staged = threading.Event()
    successor_requested = threading.Event()
    outcomes = {}

    def interrupt_predecessor(*args, **kwargs):
        result = original_validate(*args, **kwargs)
        if (kwargs.get('staged')
                and threading.current_thread().name == 'predecessor'):
            predecessor_staged.set()
            assert successor_requested.wait(5.)
            raise KeyboardInterrupt()
        return result

    def successor():
        assert predecessor_staged.wait(5.)
        successor_requested.set()
        try:
            outcomes['successor'] = provider.qualify(state, 2.5)
        except BaseException as error:
            outcomes['successor'] = error

    monkeypatch.setattr(
        provider, '_validate_acceptance_manifest', interrupt_predecessor)
    queued = threading.Thread(target=successor, name='successor')
    queued.start()
    original_name = threading.current_thread().name
    threading.current_thread().name = 'predecessor'
    try:
        with pytest.raises(KeyboardInterrupt):
            provider.qualify(state, 2.5)
    finally:
        threading.current_thread().name = original_name
    queued.join(5.)

    assert not queued.is_alive()
    assert 'latched FAILED' in str(outcomes['successor'])
    assert provider._qualification_status == 'FAILED'
    assert not any(
        json.loads(candidate.read_text()).get('status') == 'PASS'
        for candidate in tmp_path.rglob('*.tmp'))
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_terminal_failure_interruption_finalizes_before_queued_attempt(
        tmp_path, monkeypatch):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    entered = threading.Event()
    successor_requested = threading.Event()
    outcomes = {}
    original_name = threading.current_thread().name

    def interrupting_resolver(y, t):
        if threading.current_thread().name == 'terminal-predecessor':
            entered.set()
            assert successor_requested.wait(5.)
            raise KeyboardInterrupt()
        return state

    def successor():
        assert entered.wait(5.)
        successor_requested.set()
        try:
            outcomes['successor'] = provider.qualify(np.array([1.]), 2.5)
        except BaseException as error:
            outcomes['successor'] = error

    provider.bind_qualification_state(interrupting_resolver)
    queued = threading.Thread(target=successor, name='valid-successor')
    queued.start()
    threading.current_thread().name = 'terminal-predecessor'
    try:
        with pytest.raises(KeyboardInterrupt):
            provider.record_terminal_failure(
                np.array([1.]), 2.5,
                EEDFError('terminal power gate failed'))
    finally:
        threading.current_thread().name = original_name
    queued.join(5.)

    assert not queued.is_alive()
    assert 'latched FAILED' in str(outcomes['successor'])
    assert provider._qualification_status == 'FAILED'
    assert provider._qualification_record['status'] == 'FAIL'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


@pytest.mark.parametrize('caller,boundary', [
    ('qualification', 'clock'),
    ('qualification', 'run-id'),
    ('qualification', 'post-invalidation'),
    ('qualification', 'time-conversion'),
    ('qualification', 'state-resolution'),
    ('terminal', 'run-id'),
    ('terminal', 'post-invalidation'),
    ('terminal', 'time-conversion'),
    ('terminal', 'state-resolution'),
])
def test_arbitrary_pre_path_failure_replaces_known_pass(
        tmp_path, monkeypatch, caller, boundary):
    import rmgpy.solver.eedf_provider as provider_module

    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    provider.qualify(state, 2.5)
    record_path = Path(state['record_path'])
    argument = copy.deepcopy(state)
    del argument['record_path']
    del argument['scratch_root']
    attempt_time = 2.5
    injected = RuntimeError('injected pre-path failure at ' + boundary)

    if boundary == 'clock':
        original = provider_module.time.monotonic
        calls = [0]

        def fail_once():
            calls[0] += 1
            if calls[0] == 1:
                raise injected
            return original()

        monkeypatch.setattr(provider_module.time, 'monotonic', fail_once)
    elif boundary == 'run-id':
        original = provider_module.uuid.uuid4
        calls = [0]

        def fail_once():
            calls[0] += 1
            if calls[0] == 1:
                raise injected
            return original()

        monkeypatch.setattr(provider_module.uuid, 'uuid4', fail_once)
    elif boundary == 'post-invalidation':
        original = provider._invalidate_acceptance
        calls = [0]

        def fail_after_invalidation(*args, **kwargs):
            result = original(*args, **kwargs)
            calls[0] += 1
            if calls[0] == 1:
                raise injected
            return result

        monkeypatch.setattr(
            provider, '_invalidate_acceptance', fail_after_invalidation)
    elif boundary == 'time-conversion':
        class BadTime:
            def __float__(self):
                raise injected
        attempt_time = BadTime()
    else:
        monkeypatch.setattr(
            provider, '_qualification_state',
            lambda *args, **kwargs: (_ for _ in ()).throw(injected))

    from rmgpy.exceptions import PlasmaStateError
    if caller == 'qualification':
        with pytest.raises(PlasmaStateError):
            provider.qualify(argument, attempt_time)
    else:
        provider.record_terminal_failure(
            argument, attempt_time, EEDFError('terminal power gate failed'))

    assert not record_path.exists() or json.loads(
        record_path.read_text())['status'] == 'FAIL'
    assert provider._qualification_status == 'FAILED'
    assert provider._qualification_record['status'] == 'FAIL'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_export_read_failure_replaces_pass_before_terminal_refusal(
        tmp_path, monkeypatch):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    provider.qualify(state, 2.5)
    record_path = Path(state['record_path'])
    artifact = Path(provider._qualification_commit['run_directory']) / 'solver-row.json'
    original_read = Path.read_bytes

    def transient_read_error(path):
        if path == artifact:
            raise OSError(errno.EIO, 'temporary solver-evidence read error')
        return original_read(path)

    monkeypatch.setattr(Path, 'read_bytes', transient_read_error)
    with pytest.raises(EEDFError, match='manifest refuses export'):
        provider.require_qualified('export')
    monkeypatch.setattr(Path, 'read_bytes', original_read)

    assert json.loads(record_path.read_text())['status'] == 'FAIL'
    from rmgpy.exceptions import PlasmaStateError
    record = provider.record_terminal_failure(
        state, 2.5, PlasmaStateError('terminal power gate failed'))
    assert record['status'] == 'FAIL'
    assert json.loads(record_path.read_text())['status'] == 'FAIL'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def admission_test_consumer(admission, destination, payload):
    """Replay the old normal-return consumer too, so baseline failures do I/O."""
    if admission is not None:
        return admission.commit_bytes(destination, payload)
    source = destination.with_name('.legacy-export.tmp')
    source.write_bytes(payload)
    os.link(source, destination)
    source.unlink()
    return destination


def test_revocation_preserves_replacement_at_removal_boundary(tmp_path, monkeypatch):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    provider.qualify(state, 2.5)
    admission = provider.require_qualified('export')
    destination = tmp_path / 'export.json'
    foreign = tmp_path / 'foreign'
    foreign.write_bytes(b'other writer ordinary output')
    events = []
    original_link, original_unlink = os.link, Path.unlink
    original_rename, original_replace = os.rename, os.replace

    def install():
        if not events:
            before = destination.stat().st_ino
            original_replace(foreign, destination)
            events.append((before, destination.stat().st_ino))

    def unlink(path, *args, **kwargs):
        if path == destination:
            install()  # Baseline: after lstat, immediately before shared unlink.
        return original_unlink(path, *args, **kwargs)

    def rename(source, target, *args, **kwargs):
        if os.fspath(source) in (str(destination), destination.name):
            install()  # Replacement before atomic capture must also survive.
        return original_rename(source, target, *args, **kwargs)

    def link(source, target, *args, **kwargs):
        result = original_link(source, target, *args, **kwargs)
        if Path(target) == destination and not events:
            raise OSError(errno.EIO, 'ordinary failure after link')
        return result

    monkeypatch.setattr(Path, 'unlink', unlink)
    monkeypatch.setattr(os, 'rename', rename)
    monkeypatch.setattr(os, 'link', link)
    with pytest.raises((OSError, EEDFError)):
        admission.commit_bytes(destination, b'{"status":"PASS"}')
    assert len(events) == 1 and events[0][0] != events[0][1]
    assert destination.is_file(), 'cleanup deleted the raced-in foreign inode'
    assert destination.read_bytes() == b'other writer ordinary output'


def test_acceptance_gap_is_non_pass_without_unhandled_missing_file(tmp_path, monkeypatch):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    prior = provider.qualify(state, 2.5)
    destination = Path(state['record_path'])
    original_replace, original_link = os.replace, os.link
    observations = []

    def observe(source, target):
        if Path(target) == destination and str(source).endswith('.tmp'):
            absent = not destination.exists()
            try:
                provider.require_qualified('reader in publication gap')
            except EEDFError:
                observations.append((absent, provider._qualification_status))

    def replace(source, target, *args, **kwargs):
        observe(source, target)
        return original_replace(source, target, *args, **kwargs)

    def link(source, target, *args, **kwargs):
        observe(source, target)
        return original_link(source, target, *args, **kwargs)

    monkeypatch.setattr(os, 'replace', replace)
    monkeypatch.setattr(os, 'link', link)
    current = provider.qualify(state, 3.5)
    assert observations == [(True, 'CANDIDATE')]
    assert current['run_id'] != prior['run_id']
    assert current['status'] == 'PASS'


def test_successful_qualification_and_export_leave_no_private_directory(tmp_path, monkeypatch):
    previous = set(tmp_path.parent.rglob('.eedf-*.private'))
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    provider.qualify(state, 2.5)
    destination = tmp_path / 'export.json'
    provider.require_qualified('export').commit_bytes(destination, b'{"status":"PASS"}')
    assert destination.read_bytes() == b'{"status":"PASS"}'
    assert set(tmp_path.parent.rglob('.eedf-*.private')) == previous


@pytest.mark.parametrize('restoration', ['restore', 'occupied', 'no-hard-links', 'symlink'])
def test_foreign_capture_is_restored_or_retained_and_reported(tmp_path, monkeypatch, restoration):
    provider, state = public_qualification_fixture(tmp_path / 'fixture', monkeypatch)
    provider.qualify(state, 2.5)
    admission = provider.require_qualified('export')
    destination, foreign = tmp_path / 'output.json', tmp_path / 'foreign'
    if restoration == 'symlink':
        target = tmp_path / 'symlink-target'
        target.write_bytes(b'foreign data')
        foreign.symlink_to(target)
    else:
        foreign.write_bytes(b'foreign data')
    identity = foreign.lstat().st_ino
    before = set(tmp_path.parent.rglob('.eedf-*.private'))
    original_link, original_replace = os.link, os.replace
    events = []

    def handler(signum, frame):
        events.append(('signal', signum))
        provider.record_terminal_failure(state, 3., RuntimeError('shutdown'))

    def link(source, target, *args, **kwargs):
        if os.fspath(source) == 'captured':
            events.append(('restore',))
            if restoration == 'no-hard-links':
                raise OSError(errno.EOPNOTSUPP, 'hard links unsupported')
            if restoration == 'occupied' and not destination.exists():
                destination.write_bytes(b'new occupant')
        result = original_link(source, target, *args, **kwargs)
        if Path(target) == destination and not events:
            original_replace(foreign, destination)
            os.kill(os.getpid(), signal.SIGTERM)
        return result

    monkeypatch.setattr(os, 'link', link)
    with intended_signal_disposition(signal.SIGTERM, handler):
        with pytest.raises(EEDFError) as refused:
            admission.commit_bytes(destination, b'{"status":"PASS"}')
    assert events[0] == ('signal', signal.SIGTERM) and ('restore',) in events
    retained = set(tmp_path.parent.rglob('.eedf-*.private')) - before
    paths = [destination] + [path for directory in retained for path in directory.iterdir()]
    preserved = [path for path in paths if path.exists() and path.lstat().st_ino == identity]
    assert preserved and preserved[0].read_bytes() == b'foreign data'
    if restoration in ('restore', 'symlink'):
        assert preserved == [destination] and not retained
        if restoration == 'symlink':
            assert destination.is_symlink() and target.read_bytes() == b'foreign data'
    else:
        assert retained
        assert all(str(directory) in str(refused.value) for directory in retained)
        if restoration == 'occupied':
            assert destination.read_bytes() == b'new occupant'


@pytest.mark.parametrize('publication', ['acceptance', 'failure'])
def test_record_publication_preserves_raced_in_foreign_output(tmp_path, monkeypatch, publication):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    if publication == 'failure':
        provider.qualify(state, 2.5)
    destination = Path(state['record_path'])
    foreign = tmp_path / 'foreign-output'
    foreign.write_bytes(b'other writer record')
    original_link, original_replace = os.link, os.replace
    hits = []

    def link(source, target, *args, **kwargs):
        if Path(target) == destination and not hits:
            hits.append(source)
            original_replace(foreign, destination)
        return original_link(source, target, *args, **kwargs)

    monkeypatch.setattr(os, 'link', link)
    if publication == 'acceptance':
        from rmgpy.exceptions import PlasmaStateError
        with pytest.raises(PlasmaStateError):
            provider.qualify(state, 2.5)
    else:
        provider.record_terminal_failure(state, 3., RuntimeError('shutdown'))
    assert len(hits) == 1
    assert destination.read_bytes() == b'other writer record'
    assert provider._qualification_status == 'FAILED'


@pytest.mark.parametrize('lifecycle', ['fresh', 'reused'])
@pytest.mark.parametrize('signum', [signal.SIGINT, signal.SIGTERM])
def test_admission_refuses_failure_during_validation(tmp_path, monkeypatch, lifecycle, signum):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    provider.qualify(state, 2.5)
    if lifecycle == 'reused':
        provider.qualify(state, 3.5)
    events, failures = [], []

    def handler(received, frame):
        events.append(('received', received))
        failures.append(provider.record_terminal_failure(state, 4., RuntimeError('shutdown')))

    def trace(frame, event, arg):
        if (event == 'return' and frame.f_code is provider._validate_acceptance_manifest.__func__.__code__
                and not events):
            events.append(('sent', signum))
            os.kill(os.getpid(), signum)
        return trace

    raised = None
    with intended_signal_disposition(signum, handler):
        try:
            sys.settrace(trace)
            provider.require_qualified('export')
        except EEDFError as error:
            raised = error
        finally:
            sys.settrace(None)
    assert events == [('sent', signum), ('received', signum)]
    assert raised is not None, 'validation returned admission after terminal failure'
    disk = json.loads(Path(state['record_path']).read_text())
    assert disk['status'] == 'FAIL' and disk['run_id'] == failures[0]['run_id']


@pytest.mark.parametrize('boundary', ['gate-return', 'wrapper-return', 'caller', 'before-link', 'after-link'])
@pytest.mark.parametrize('signum', [signal.SIGINT, signal.SIGTERM])
def test_admission_revocation_prevents_consumer_commit(tmp_path, monkeypatch, boundary, signum):
    import rmgpy.solver.eedf_provider as module
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    snapshot = provider.qualify(state, 2.5)
    payload = json.dumps({'status': snapshot['status'], 'run_id': snapshot['run_id']}).encode()
    destination = tmp_path / 'export.json'
    events, failures = [], []

    def handler(received, frame):
        events.append(('received', received))
        failures.append(provider.record_terminal_failure(state, 4., RuntimeError('shutdown')))

    def deliver(site):
        if site == boundary and not events:
            events.append(('sent', signum, site))
            os.kill(os.getpid(), signum)

    original_link = os.link

    def link(source, target, *args, **kwargs):
        if Path(target) == destination:
            deliver('before-link')
        result = original_link(source, target, *args, **kwargs)
        if Path(target) == destination:
            deliver('after-link')
        return result

    def trace(frame, event, arg):
        if event == 'return':
            if frame.f_code is module.EEDFProvider.require_qualified.__wrapped__.__code__:
                deliver('gate-return')
            elif (frame.f_code is module.EEDFProvider.require_qualified.__code__
                  and frame.f_locals.get('provider') is provider):
                deliver('wrapper-return')
        return trace

    monkeypatch.setattr(os, 'link', link)
    refused = False
    with intended_signal_disposition(signum, handler):
        try:
            sys.settrace(trace)
            admission = provider.require_qualified('export')
            sys.settrace(None)
            deliver('caller')
            admission_test_consumer(admission, destination, payload)
        except (EEDFError, FileNotFoundError):
            refused = True
        finally:
            sys.settrace(None)
    assert len(events) == 2 and events[0][0] == 'sent' and events[1][0] == 'received'
    assert refused, 'normal-return admission committed output for a failed run'
    assert not destination.exists()
    assert not list(tmp_path.glob('.*export*'))
    disk = json.loads(Path(state['record_path']).read_text())
    assert disk['status'] == 'FAIL' and disk['run_id'] == failures[0]['run_id']


@pytest.mark.parametrize('snapshot_kind', ['record', 'status', 'manifest', 'path'])
def test_admission_snapshots_cannot_authorize_output(tmp_path, monkeypatch, snapshot_kind):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    record = provider.qualify(state, 2.5)
    admission = provider.require_qualified('export')
    snapshot = {'record': record, 'status': provider._qualification_status,
                'manifest': provider._manifest_snapshot, 'path': provider.qualification_record_path}[snapshot_kind]
    original = provider._json_value(snapshot)
    events = []

    def handler(signum, frame):
        events.append(signum)
        provider.record_terminal_failure(state, 4., RuntimeError('shutdown'))

    with intended_signal_disposition(signal.SIGTERM, handler):
        os.kill(os.getpid(), signal.SIGTERM)
    assert events == [signal.SIGTERM]
    assert provider._json_value(snapshot) == original  # C1: historical snapshots do not become live views.
    destination = tmp_path / 'snapshot-export.json'
    refused = False
    try:
        admission_test_consumer(admission, destination, json.dumps({'status': record['status'], 'run_id': record['run_id']}).encode())
    except EEDFError:
        refused = True
    assert refused and not destination.exists()
    with pytest.raises(EEDFError):
        provider.require_qualified('export from snapshot')


def test_admission_commits_bytes_once_without_boolean_or_overwrite(tmp_path, monkeypatch):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    snapshot = provider.qualify(state, 2.5)
    admission = provider.require_qualified('export')
    destination = tmp_path / 'export.json'
    payload = json.dumps({'status': snapshot['status'], 'run_id': snapshot['run_id']}).encode()
    admission_test_consumer(admission, destination, payload)
    assert destination.read_bytes() == payload
    assert admission is not None, 'normal return is not an admission capability'
    with pytest.raises(TypeError):
        bool(admission)
    with pytest.raises(TypeError):
        admission.commit_bytes(tmp_path / 'callback', lambda: payload)
    with pytest.raises((EEDFError, FileExistsError)):
        admission.commit_bytes(destination, b'unrelated replacement')
    assert destination.read_bytes() == payload


@pytest.mark.parametrize('second_policy', ['raising', 'returning'])
def test_admission_second_signal_during_commit_revocation(tmp_path, monkeypatch, second_policy, record_property):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    record = provider.qualify(state, 2.5)
    admission = provider.require_qualified('export')
    destination = tmp_path / 'export.json'
    events = []
    active = False
    original_read, original_capture, original_link = Path.read_bytes, os.rename, os.link

    def second(signum, frame):
        events.append(('second-received', signum))
        if second_policy == 'raising':
            raise KeyboardInterrupt('second signal during revocation')
        provider.record_terminal_failure(state, 4., RuntimeError('secondary shutdown'))

    def first(signum, frame):
        nonlocal active
        events.append(('first-received', signum))
        active = True
        try:
            provider.require_qualified('shutdown export')
        except BaseException as error:
            events.append(('caught', type(error).__name__))
        finally:
            active = False

    def read(path):
        if active and path == Path(state['record_path']):
            raise OSError(errno.EIO, 'shutdown manifest read')
        return original_read(path)

    def capture(path, target_name, *args, **kwargs):
        target = destination if admission is not None else Path(state['record_path'])
        if active and os.fspath(path) == target.name and not any(item[0] == 'second-sent' for item in events):
            events.append(('second-sent', signal.SIGINT))
            os.kill(os.getpid(), signal.SIGINT)
        return original_capture(path, target_name, *args, **kwargs)

    def link(source, target, *args, **kwargs):
        result = original_link(source, target, *args, **kwargs)
        if Path(target) == destination:
            events.append(('first-sent', signal.SIGTERM))
            os.kill(os.getpid(), signal.SIGTERM)
        return result

    monkeypatch.setattr(Path, 'read_bytes', read)
    monkeypatch.setattr(os, 'rename', capture)
    monkeypatch.setattr(os, 'link', link)
    refused = False
    with intended_signal_disposition(signal.SIGTERM, first):
        with intended_signal_disposition(signal.SIGINT, second):
            try:
                admission_test_consumer(admission, destination, json.dumps({'status': record['status'], 'run_id': record['run_id']}).encode())
            except (EEDFError, FileNotFoundError):
                refused = True
    assert sum(item[0] == 'first-received' for item in events) == 1
    assert sum(item[0] == 'second-received' for item in events) == 1
    assert refused and not destination.exists(), events
    disk = json.loads(Path(state['record_path']).read_text())
    assert disk['status'] == 'FAIL' and disk['run_id'] == record['run_id']

    record_property('hook_name', 'export-capture/second-signal')
    record_property('hook_hits', sum(event[0] == 'second-sent' for event in events))


def test_admission_superseded_generation_cannot_commit(tmp_path, monkeypatch):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    snapshot = provider.qualify(state, 2.5)
    admission = provider.require_qualified('export')
    nested = []
    events = []

    def handler(signum, frame):
        events.append(signum)
        nested.append(provider.qualify(state, 3.5))

    with intended_signal_disposition(signal.SIGTERM, handler):
        os.kill(os.getpid(), signal.SIGTERM)
    assert events == [signal.SIGTERM]
    assert snapshot['status'] == nested[0]['status'] == 'PASS'
    assert snapshot['run_id'] != nested[0]['run_id']
    destination = tmp_path / 'old-export.json'
    refused = False
    try:
        admission_test_consumer(admission, destination, b'{"status":"PASS"}')
    except EEDFError:
        refused = True
    assert refused and not destination.exists()
    provider.require_qualified('current export').commit_bytes(
        tmp_path / 'current-export.json', b'{"status":"PASS"}')


@pytest.mark.parametrize('changed_path', ['source', 'destination'])
def test_admission_ambiguous_inode_refuses_without_deleting_replacement(
        tmp_path, monkeypatch, changed_path):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    provider.qualify(state, 2.5)
    admission = provider.require_qualified('export')
    destination = tmp_path / 'export.json'
    replacement = []
    events = []
    original_link = os.link

    def handler(signum, frame):
        events.append(signum)
        provider.record_terminal_failure(state, 4., RuntimeError('shutdown'))

    def link(source, target, *args, **kwargs):
        result = original_link(source, target, *args, **kwargs)
        if Path(target) == destination:
            path = Path(source) if changed_path == 'source' else destination
            foreign = tmp_path / 'foreign-inode'
            foreign.write_bytes(b'belongs to another writer')
            os.replace(foreign, path)
            replacement.append(path)
            os.kill(os.getpid(), signal.SIGTERM)
        return result

    monkeypatch.setattr(os, 'link', link)
    refused = None
    with intended_signal_disposition(signal.SIGTERM, handler):
        try:
            admission_test_consumer(admission, destination, b'{"status":"PASS"}')
        except EEDFError as error:
            refused = error
    assert events == [signal.SIGTERM]
    assert refused is not None, 'ambiguous export ownership was treated as success'
    assert replacement[0].read_bytes() == b'belongs to another writer'
    if changed_path == 'source':
        assert not destination.exists()


@pytest.mark.parametrize('invalidation_forbidden', [False, True])
def test_admission_unlink_and_invalidation_limit_is_explicit(
        tmp_path, monkeypatch, invalidation_forbidden, record_property):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    provider.qualify(state, 2.5)
    admission = provider.require_qualified('export')
    destination = tmp_path / 'export.json'
    owned = set()
    events = []
    original_link, original_capture, original_truncate = os.link, os.rename, os.ftruncate
    capture_hits = []

    def handler(signum, frame):
        events.append(signum)
        provider.record_terminal_failure(state, 4., RuntimeError('shutdown'))

    def link(source, target, *args, **kwargs):
        result = original_link(source, target, *args, **kwargs)
        if Path(target) == destination:
            owned.update((Path(source), Path(target)))
            os.kill(os.getpid(), signal.SIGTERM)
        return result

    def capture(path, target_name, *args, **kwargs):
        if os.fspath(path) == destination.name:
            capture_hits.append(path)
            raise PermissionError('injected forbidden export capture')
        return original_capture(path, target_name, *args, **kwargs)

    def truncate(fd, length):
        if invalidation_forbidden:
            raise PermissionError('injected forbidden inode invalidation')
        return original_truncate(fd, length)

    monkeypatch.setattr(os, 'link', link)
    monkeypatch.setattr(os, 'rename', capture)
    monkeypatch.setattr(os, 'ftruncate', truncate)
    refused = None
    with intended_signal_disposition(signal.SIGTERM, handler):
        try:
            admission_test_consumer(admission, destination, b'{"status":"PASS"}')
        except (EEDFError, PermissionError) as error:
            refused = error
    assert events == [signal.SIGTERM]
    assert capture_hits, 'namespace-denial hook was not hit'
    assert isinstance(refused, EEDFError), 'filesystem limit did not produce an explicit admission refusal'
    assert provider._qualification_status == 'FAILED'
    if invalidation_forbidden:
        assert 'unlink and inode invalidation both failed' in str(refused)
        assert json.loads(destination.read_bytes())['status'] == 'PASS'  # Named filesystem limit.
    else:
        assert 'invalidated but unlink failed' in str(refused)
        assert json.loads(destination.read_bytes())['status'] == 'FAIL'
    assert 'publication_error' in provider._qualification_record

    record_property('hook_name', 'export-capture/denial')
    record_property('hook_hits', len(capture_hits))


@pytest.mark.parametrize('lifecycle', ['fresh', 'reused'])
@pytest.mark.parametrize('signum', [signal.SIGINT, signal.SIGTERM])
@pytest.mark.parametrize('action', ['terminal', 'export', 'qualify'])
@pytest.mark.parametrize('boundary', [
    'completion-entry', 'before-stage', 'after-stage', 'before-rename',
    'after-rename', 'completion-return', 'guard-return',
])
def test_returning_handler_cannot_revive_superseded_qualification(
        tmp_path, monkeypatch, lifecycle, signum, action, boundary, record_property):
    import rmgpy.solver.eedf_provider as provider_module
    from rmgpy.tools.eedf import signal_guard

    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    if lifecycle == 'reused':
        provider.qualify(state, 2.5)
    nested_state = dict(state, record_path=str(tmp_path / 'nested.json'))
    events = []
    nested = []
    fired = False
    outcome = {}

    def handler(received, frame):
        events.append(('received', received, provider._qualification_status))
        if action == 'terminal':
            nested.append(provider.record_terminal_failure(
                state, 2.5, RuntimeError('handler requested shutdown')))
        elif action == 'qualify':
            nested.append(provider.qualify(nested_state, 2.5))
        else:
            try:
                provider.require_qualified('handler export')
            except EEDFError:
                events.append(('export-refused', provider._qualification_status))

    def deliver(site):
        nonlocal fired
        if not fired and site == boundary:
            fired = True
            events.append(('sent', signum, site))
            os.kill(os.getpid(), signum)

    original_write = provider._stage_record
    original_replace = provider_module.os.link

    def stage(path, value):
        staged_pass = str(path.source).endswith('.tmp') and value.get('status') == 'PASS'
        if staged_pass:
            deliver('before-stage')
        result = original_write(path, value)
        if staged_pass:
            deliver('after-stage')
        return result

    def replace(source, target, *args, **kwargs):
        staged_pass = (str(source).endswith('.tmp')
                       and not str(source).endswith('.failure.tmp'))
        if staged_pass:
            # Runs inside the actual publication link, after any caller-side check.
            deliver('before-rename')
        result = original_replace(source, target, *args, **kwargs)
        if staged_pass:
            deliver('after-rename')
        return result

    def trace(frame, event, arg):
        if frame.f_code.co_filename == provider_module.__file__:
            if frame.f_code.co_name == 'complete':
                if event == 'line':
                    deliver('completion-entry')
                elif event == 'return':
                    deliver('completion-return')
        elif (frame.f_code.co_filename == signal_guard.__file__
              and frame.f_code.co_name == 'run_guarded_attempt'
              and event == 'return'):
            deliver('guard-return')
        return trace

    monkeypatch.setattr(provider, '_stage_record', stage)
    monkeypatch.setattr(provider_module.os, 'link', replace)
    with intended_signal_disposition(signum, handler):
        try:
            sys.settrace(trace)
            outcome['record'] = provider.qualify(state, 2.5)
        except BaseException as error:
            outcome['error'] = error
        finally:
            sys.settrace(None)
        assert signal.getsignal(signum) is handler

    assert fired, events
    assert len([event for event in events if event[0] == 'sent']) == 1
    assert len([event for event in events if event[0] == 'received']) == 1
    assert 'error' not in outcome, outcome.get('error')
    record = outcome['record']
    assert record == provider._qualification_record
    assert record['run_id'] == provider._qualification_run_id
    if action == 'terminal':
        assert provider._qualification_status == 'FAILED'
        assert record['status'] == nested[0]['status'] == 'FAIL'
        assert record['run_id'] == nested[0]['run_id']
        assert json.loads(Path(state['record_path']).read_text())['status'] == 'FAIL'
    else:
        assert provider._qualification_status == record['status'] == 'PASS'
        provider.require_qualified('export after handler')
        if action == 'qualify':
            assert record['run_id'] == nested[0]['run_id']
            assert provider.qualification_record_path == Path(nested_state['record_path'])
            old = Path(state['record_path'])
            assert not old.exists() or json.loads(old.read_text())['status'] != 'PASS'
        elif boundary in ('completion-entry', 'before-stage', 'after-stage',
                          'before-rename'):
            assert ('export-refused', 'CANDIDATE') in events
    assert not any(json.loads(path.read_text()).get('status') == 'PASS'
                   for path in tmp_path.rglob('*.tmp'))

    record_property('hook_name', 'reentry/'+boundary)
    record_property('hook_hits', sum(event[0] == 'sent' for event in events))


@pytest.mark.parametrize('serializable', [True, False])
def test_reentrant_failure_preserves_run_id_when_staged_unlink_is_forbidden(
        tmp_path, monkeypatch, serializable, record_property):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    record_path = Path(state['record_path'])
    original_replace = os.link
    sent = []
    failed_run = []

    def handler(signum, frame):
        record_path.parent.chmod(0o500)
        try:
            failed_state = dict(state)
            if not serializable:
                failed_state['diagnostic'] = object()
            record = provider.record_terminal_failure(
                failed_state, 2.5, RuntimeError('shutdown with read-only directory'))
            failed_run.append(record['run_id'])
        finally:
            record_path.parent.chmod(0o700)

    def replace(source, target, *args, **kwargs):
        if (not sent and str(source).endswith('.tmp')
                and not str(source).endswith('.failure.tmp')):
            sent.append(signal.SIGTERM)
            os.kill(os.getpid(), signal.SIGTERM)
        return original_replace(source, target, *args, **kwargs)

    monkeypatch.setattr(os, 'link', replace)
    with intended_signal_disposition(signal.SIGTERM, handler):
        record = provider.qualify(state, 2.5)
    assert sent == [signal.SIGTERM] and len(failed_run) == 1
    assert record['status'] == 'FAIL'
    assert record['run_id'] == provider._qualification_run_id == failed_run[0]
    disk = json.loads(record_path.read_text())
    assert disk['status'] == 'FAIL' and disk['run_id'] == failed_run[0]

    record_property('hook_name', 'before-public-link/read-only')
    record_property('hook_hits', len(sent))


@pytest.mark.parametrize('boundary', ['revocation', 'response'])
@pytest.mark.parametrize('second_policy', ['raising', 'returning'])
def test_second_signal_during_failure_bookkeeping_preserves_failure(
        tmp_path, monkeypatch, boundary, second_policy):
    import rmgpy.solver.eedf_provider as provider_module

    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    events = []
    primary = False
    second_sent = False
    expected_run = []
    original_validate = provider._validate_acceptance_manifest
    original_guard = provider_module.run_guarded_attempt

    def validate(*args, **kwargs):
        if primary:
            raise OSError(errno.EIO, 'handler export read failed')
        return original_validate(*args, **kwargs)

    if boundary == 'revocation':
        target = getattr(provider, '_revoke_publications', provider._enter_failed)
        marker = ('context.file.revoke(' if hasattr(provider, '_revoke_publications')
                  else 'for candidate in candidates:')
    else:
        target = getattr(provider, '_sync_qualification_responses', provider._enter_failed)
        marker = ('context.record.update(record)' if hasattr(provider, '_sync_qualification_responses')
                  else 'self._qualification_record = record')
    source, start = inspect.getsourcelines(target)
    target_line = next(start + index for index, line in enumerate(source)
                       if marker in line)

    def trace(frame, event, arg):
        nonlocal second_sent
        if (primary and not second_sent and event == 'line'
                and frame.f_code is target.__func__.__code__
                and frame.f_lineno == target_line):
            second_sent = True
            events.append(('second-sent', boundary, frame.f_lineno))
            os.kill(os.getpid(), signal.SIGINT)
        return trace

    def second_handler(signum, frame):
        events.append(('second-received', signum, frame.f_code.co_name))
        if second_policy == 'raising':
            raise KeyboardInterrupt('second signal inside bookkeeping')
        provider.record_terminal_failure(
            state, 2.5, RuntimeError('second handler requested shutdown'))

    def first_handler(signum, frame):
        nonlocal primary
        primary = True
        events.append(('first-received', signum))
        expected_run.append(provider._qualification_run_id)
        try:
            provider.require_qualified('shutdown export')
        except BaseException as error:
            # An ordinary returning handler may tolerate its export failure,
            # including the interrupt that arrived inside failure finalization.
            events.append(('first-caught', type(error).__name__))
        finally:
            primary = False

    def guard(lock, operation, finalize, complete):
        def after_complete(record):
            result = complete(record)
            if not events:
                events.append(('first-sent', signal.SIGTERM))
                os.kill(os.getpid(), signal.SIGTERM)
            return result
        return original_guard(lock, operation, finalize, after_complete)

    monkeypatch.setattr(provider, '_validate_acceptance_manifest', validate)
    monkeypatch.setattr(provider_module, 'run_guarded_attempt', guard)
    with intended_signal_disposition(signal.SIGTERM, first_handler):
        with intended_signal_disposition(signal.SIGINT, second_handler):
            try:
                sys.settrace(trace)
                record = provider.qualify(state, 2.5)
            finally:
                sys.settrace(None)
    assert second_sent, events
    assert len([event for event in events if event[0] == 'first-received']) == 1
    assert len([event for event in events if event[0] == 'second-received']) == 1
    assert provider._qualification_status == 'FAILED'
    assert record['status'] == provider._qualification_record['status'] == 'FAIL'
    assert record['run_id'] == provider._qualification_run_id == expected_run[0]
    disk = json.loads(Path(state['record_path']).read_text())
    assert disk['status'] == 'FAIL' and disk['run_id'] == expected_run[0]
    assert 'handler export read failed' in record['error']
    assert not any(json.loads(path.read_text()).get('status') == 'PASS'
                   for path in tmp_path.rglob('*.tmp'))


@pytest.mark.parametrize('lifecycle', ['fresh', 'reused'])
@pytest.mark.parametrize('schedule', [
    'TERM_TERM', 'INT_TERM', 'raw_INT_TERM', 'pending_INT_TERM',
])
def test_owned_cancellation_during_restore_cannot_leave_pass(
        tmp_path, monkeypatch, lifecycle, schedule):
    from rmgpy.tools.eedf import signal_guard

    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    provider.qualify(state, 2.5)
    if lifecycle == 'fresh':
        provider = make_provider(
            provider._table.path,
            provider._json_value(provider.runtime_metadata['row_inputs']))
    record_path = Path(state['record_path'])
    events_path = tmp_path / 'signal-events.json'
    original_signal = signal.signal
    original_raw_restore = signal_guard._restore_raw_sigaction
    events = []

    def send(signum):
        events.append(int(signum))
        events_path.write_text(json.dumps(events))
        os.kill(os.getpid(), signum)

    def restore(signum, handler):
        if signum != signal.SIGTERM or handler != signal.SIG_DFL:
            result = original_signal(signum, handler)
            if (schedule == 'pending_INT_TERM' and signum == signal.SIGINT
                    and getattr(handler, '__name__', '') == 'defer_signal'):
                send(signal.SIGINT)
            return result
        if schedule == 'TERM_TERM':
            send(signal.SIGTERM)
        result = original_signal(signum, handler)
        if schedule == 'INT_TERM':
            send(signal.SIGINT)
        if schedule != 'raw_INT_TERM':
            send(signal.SIGTERM)
        return result

    def restore_raw(signum, action):
        if signum == signal.SIGTERM and schedule == 'raw_INT_TERM':
            send(signal.SIGINT)
            send(signal.SIGTERM)
        return original_raw_restore(signum, action)

    child = os.fork()
    if child == 0:
        try:
            original_signal(signal.SIGTERM, signal.SIG_DFL)
            original_signal(signal.SIGINT, signal.default_int_handler)
            monkeypatch.setattr(signal, 'signal', restore)
            monkeypatch.setattr(signal_guard, '_restore_raw_sigaction', restore_raw)
            provider.qualify(state, 2.5)
        except BaseException:
            os._exit(2)
        os._exit(0)
    _, status = os.waitpid(child, 0)

    assert os.waitstatus_to_exitcode(status) == -signal.SIGTERM
    assert json.loads(events_path.read_text()) == (
        [15, 15] if schedule == 'TERM_TERM' else [2, 15])
    assert not record_path.exists() or json.loads(record_path.read_text())['status'] != 'PASS'
    assert not list(record_path.parent.glob('.*.tmp'))


def test_qualification_export_census_locator_matches_source():
    import ast
    import rmgpy.solver.eedf_provider as provider_module

    root = Path(__file__).resolve().parents[3]
    census = json.loads(
        (root / 'test/rmgpy/test_data/excited_export_census.json').read_text())
    row = census['search_exclusions'][
        'rmgpy/solver/eedf_provider.py:EEDFProvider._run_qualification']
    definitions = [node for node in ast.walk(ast.parse(
        Path(provider_module.__file__).read_text()))
        if isinstance(node, ast.FunctionDef) and node.name == '_run_qualification']
    assert len(definitions) == 1
    assert row['line'] == definitions[0].lineno


def test_returning_deferred_term_returns_finalized_failure(
        tmp_path, monkeypatch):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    original_prepare = provider._run_qualification
    handler_calls = []

    def prepare_then_interrupt(*args, **kwargs):
        result = original_prepare(*args, **kwargs)
        os.kill(os.getpid(), signal.SIGTERM)
        return result

    monkeypatch.setattr(provider, '_run_qualification', prepare_then_interrupt)
    with intended_signal_disposition(
            signal.SIGTERM,
            lambda signum, frame: handler_calls.append(signum)):
        record = provider.qualify(state, 2.5)

    assert handler_calls == [signal.SIGTERM]
    assert record['status'] == 'FAIL'
    assert provider._qualification_status == 'FAILED'
    assert json.loads(Path(state['record_path']).read_text())['status'] == 'FAIL'


def test_term_arriving_during_release_invalidates_before_redelivery(
        tmp_path, monkeypatch):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    original_prepare = provider._run_qualification
    original_signal = signal.signal
    handler_calls = []
    armed = [False]

    def returning_handler(signum, frame):
        handler_calls.append(signum)

    def prepare_then_arm(*args, **kwargs):
        result = original_prepare(*args, **kwargs)
        armed[0] = True
        return result

    def interrupt_during_restore(signum, handler):
        if (armed[0] and signum == signal.SIGTERM
                and handler is returning_handler):
            armed[0] = False
            os.kill(os.getpid(), signal.SIGTERM)
        return original_signal(signum, handler)

    monkeypatch.setattr(provider, '_run_qualification', prepare_then_arm)
    with intended_signal_disposition(signal.SIGTERM, returning_handler):
        monkeypatch.setattr(signal, 'signal', interrupt_during_restore)
        record = provider.qualify(state, 2.5)

    assert handler_calls == [signal.SIGTERM]
    assert not armed[0]
    assert record['status'] == 'FAIL'
    assert provider._qualification_status == 'FAILED'
    assert json.loads(Path(state['record_path']).read_text())['status'] == 'FAIL'


def test_default_term_during_release_invalidates_before_process_exit(
        tmp_path, monkeypatch):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    original_prepare = provider._run_qualification
    original_signal = signal.signal
    armed = [False]
    delivery = tmp_path / 'term-delivered'

    def prepare_then_arm(*args, **kwargs):
        result = original_prepare(*args, **kwargs)
        armed[0] = True
        return result

    def interrupt_during_restore(signum, handler):
        if (armed[0] and signum == signal.SIGTERM
                and handler == signal.SIG_DFL):
            armed[0] = False
            delivery.write_text('SIGTERM')
            os.kill(os.getpid(), signal.SIGTERM)
        return original_signal(signum, handler)

    monkeypatch.setattr(provider, '_run_qualification', prepare_then_arm)
    monkeypatch.setattr(signal, 'signal', interrupt_during_restore)
    child = os.fork()
    if child == 0:
        try:
            original_signal(signal.SIGTERM, signal.SIG_DFL)
            provider.qualify(state, 2.5)
        except BaseException:
            os._exit(2)
        os._exit(0)

    _, status = os.waitpid(child, 0)

    assert os.waitstatus_to_exitcode(status) == -signal.SIGTERM
    assert delivery.read_text() == 'SIGTERM'
    assert json.loads(Path(state['record_path']).read_text())['status'] == 'FAIL'


@pytest.mark.parametrize('caller', ['qualification', 'terminal'])
def test_signal_after_attempt_lock_acquire_refuses_queued_successor(
        tmp_path, monkeypatch, caller):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    provider.qualify(state, 2.5)
    underlying = provider._qualification_lock
    entered = threading.Event()
    successor_requested = threading.Event()
    fired = [False]
    outcomes = {}

    class InterruptingLock:
        def acquire(self):
            underlying.acquire()
            if (threading.current_thread() is threading.main_thread()
                    and not fired[0]):
                fired[0] = True
                entered.set()
                assert successor_requested.wait(5.)
                os.kill(os.getpid(), signal.SIGINT)
            return True

        def release(self):
            underlying.release()

        def __enter__(self):
            try:
                self.acquire()
            except BaseException:
                self.release()
                raise
            return self

        def __exit__(self, *args):
            self.release()

    provider._qualification_lock = InterruptingLock()

    def successor():
        assert entered.wait(5.)
        successor_requested.set()
        try:
            outcomes['successor'] = provider.qualify(state, 2.5)['status']
        except BaseException as error:
            outcomes['successor'] = error

    queued = threading.Thread(target=successor, name='queued-qualifier')
    queued.start()
    from rmgpy.exceptions import PlasmaStateError
    with intended_signal_disposition(signal.SIGINT, signal.default_int_handler):
        with pytest.raises(KeyboardInterrupt):
            if caller == 'qualification':
                provider.qualify(state, 2.5)
            else:
                provider.record_terminal_failure(
                    state, 2.5,
                    PlasmaStateError('terminal power gate failed'))
    queued.join(5.)

    assert fired[0]
    assert not queued.is_alive()
    assert 'latched FAILED' in str(outcomes['successor'])
    assert provider._qualification_status == 'FAILED'
    assert json.loads(Path(state['record_path']).read_text())['status'] == 'FAIL'


def test_interrupted_signal_install_restores_every_installed_handler(
        tmp_path, monkeypatch):
    from rmgpy.tools.eedf import signal_guard

    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    provider.qualify(state, 2.5)
    original_signal = signal.signal

    with intended_signal_disposition(signal.SIGTERM, signal.SIG_DFL):
        with intended_signal_disposition(signal.SIGINT, signal.default_int_handler):
            expected_term = signal_guard._raw_sigaction(signal.SIGTERM)
            interrupted = [False]

            def interrupt_int_install(signum, handler):
                if (signum == signal.SIGINT and callable(handler)
                        and handler is not signal.default_int_handler
                        and not interrupted[0]):
                    interrupted[0] = True
                    raise KeyboardInterrupt()
                return original_signal(signum, handler)

            monkeypatch.setattr(signal, 'signal', interrupt_int_install)
            with pytest.raises(KeyboardInterrupt):
                provider.qualify(state, 2.5)
            actual_term = signal_guard._raw_sigaction(signal.SIGTERM)

            assert interrupted[0]
            assert signal.getsignal(signal.SIGTERM) == signal.SIG_DFL
            assert actual_term.handler == expected_term.handler
            assert actual_term.flags == expected_term.flags
            assert provider._qualification_status == 'FAILED'
            assert json.loads(Path(state['record_path']).read_text())['status'] == 'FAIL'


@pytest.mark.parametrize('caller', ['qualification', 'terminal'])
@pytest.mark.parametrize('failure', ['sigint', 'oserror'])
def test_fresh_early_failure_invalidates_caller_named_record(
        tmp_path, monkeypatch, caller, failure):
    import rmgpy.solver.eedf_provider as provider_module

    owner, state = public_qualification_fixture(tmp_path, monkeypatch)
    owner.qualify(state, 2.5)
    record_path = Path(state['record_path'])
    fresh = make_provider(
        owner._table.path,
        owner._json_value(owner.runtime_metadata['row_inputs']))
    original_uuid = provider_module.uuid.uuid4
    fired = [False]

    def fail_run_id():
        if not fired[0]:
            fired[0] = True
            if failure == 'sigint':
                os.kill(os.getpid(), signal.SIGINT)
            raise OSError(errno.EIO, 'temporary run-id failure')
        return original_uuid()

    monkeypatch.setattr(provider_module.uuid, 'uuid4', fail_run_id)
    from rmgpy.exceptions import PlasmaStateError
    disposition = (intended_signal_disposition(
        signal.SIGINT, signal.default_int_handler)
        if failure == 'sigint' else nullcontext())
    with disposition:
        try:
            if caller == 'qualification':
                fresh.qualify(copy.deepcopy(state), 2.5)
            else:
                fresh.record_terminal_failure(
                    copy.deepcopy(state), 2.5,
                    PlasmaStateError('terminal power gate failed'))
        except BaseException:
            pass

    assert fired[0]
    assert json.loads(record_path.read_text())['status'] == 'FAIL'
    assert fresh.qualification_record_path == record_path
    assert fresh._qualification_status == 'FAILED'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        fresh.require_qualified('export')


@pytest.mark.parametrize('caller', ['qualification', 'terminal'])
def test_signal_at_first_attempt_line_finalizes_before_queued_attempt(
        tmp_path, monkeypatch, caller):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    provider.qualify(state, 2.5)
    record_path = Path(state['record_path'])
    entered = threading.Event()
    successor_requested = threading.Event()
    outcomes = {}
    target = (EEDFProvider._run_qualification if caller == 'qualification'
              else EEDFProvider.record_terminal_failure)
    target_code = inspect.unwrap(target).__code__
    fired = [False]

    def successor():
        assert entered.wait(5.)
        successor_requested.set()
        try:
            outcomes['successor'] = provider.qualify(state, 2.5)['status']
        except BaseException as error:
            outcomes['successor'] = error

    queued = threading.Thread(target=successor, name='queued-qualifier')
    queued.start()

    def trace(frame, event, argument):
        if event == 'line' and frame.f_code is target_code and not fired[0]:
            fired[0] = True
            entered.set()
            assert successor_requested.wait(5.)
            os.kill(os.getpid(), signal.SIGINT)
        return trace

    sys.settrace(trace)
    from rmgpy.exceptions import PlasmaStateError
    try:
        with intended_signal_disposition(
                signal.SIGINT, signal.default_int_handler):
            with pytest.raises(KeyboardInterrupt):
                if caller == 'qualification':
                    provider.qualify(state, 2.5)
                else:
                    provider.record_terminal_failure(
                        state, 2.5,
                        PlasmaStateError('terminal power gate failed'))
    finally:
        sys.settrace(None)
    queued.join(5.)

    assert fired[0]
    assert not queued.is_alive()
    assert 'latched FAILED' in str(outcomes['successor'])
    assert provider._qualification_status == 'FAILED'
    assert json.loads(record_path.read_text())['status'] == 'FAIL'


def test_run_prepared_refuses_an_executable_changed_by_the_run(
        tmp_path, monkeypatch):
    import rmgpy.tools.eedf.loki as loki

    path, _ = fixture_provider_table(tmp_path / 'table')
    spec = json.loads((path / 'generation_spec.json').read_text())
    binary = tmp_path / 'loki'
    binary.write_bytes(b'original executable')
    binary.chmod(0o755)
    spec['binary'] = {'path': str(binary), 'sha256': _file_hash(binary)}
    spec['scratch_root'] = str(tmp_path / 'jobs')
    spec['solver_options']['numerics'] = {
        'energyGrid': {'maxEnergy': 1., 'cellNumber': 2},
        'maxPowerBalanceRelError': 1.e-3,
    }
    spec['tolerances']['normalization_atol'] = 1.e-9
    driver = object.__new__(loki.LoKIDriver)
    driver.spec = spec
    driver.setups = {}
    bundle = driver.prepare({}, [1.], 'changed-binary', verify_dependencies=False)

    monkeypatch.setattr(loki, '_check_solver', lambda request: None)
    monkeypatch.setattr(loki, '_physical_properties', lambda request: None)
    monkeypatch.setattr(loki, '_solver_environment', lambda request: {})

    def execute(*args, **kwargs):
        binary.write_bytes(b'changed executable')
        output = Path(kwargs['cwd']) / 'output' / 'solve'
        output.mkdir(parents=True)
        (output / 'reducedField_1').write_text('fixture')
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(loki.subprocess, 'run', execute)
    monkeypatch.setattr(loki, 'parse_output', lambda *args: {
        'EN_Td': 1.,
        'power_groups': {'relative_balance_percent': 0.},
    })

    with pytest.raises(FingerprintMismatch, match='binary changed during execution'):
        driver.run_prepared(bundle)


def test_concurrent_qualification_attempts_do_not_share_a_latch(tmp_path):
    path, model = fixture_provider_table(tmp_path / 'table')
    provider = make_provider(path, model)
    first_state = qualification_state(tmp_path / 'first', provider)
    second_state = qualification_state(tmp_path / 'second', provider)
    first_entered = threading.Event()
    second_entered = threading.Event()
    release_first = threading.Event()
    calls = []
    errors = []
    calls_lock = threading.Lock()

    def serialized_runner(bundle):
        with calls_lock:
            calls.append(bundle.job_name)
            call_number = len(calls)
        if call_number == 1:
            first_entered.set()
            assert release_first.wait(5.)
        else:
            second_entered.set()
        return direct_runner(provider)(bundle)

    def attempt(state):
        try:
            qualify_provider(provider, state, 2.5, runner=serialized_runner)
        except Exception as error:
            errors.append(error)

    first = threading.Thread(target=attempt, args=(first_state,))
    second = threading.Thread(target=attempt, args=(second_state,))
    first.start()
    assert first_entered.wait(5.)
    second.start()
    overlapped = second_entered.wait(.25)
    release_first.set()
    first.join(5.)
    second.join(5.)

    assert not overlapped
    assert not first.is_alive() and not second.is_alive()
    assert errors == []
    assert len(calls) == 2


def test_queued_attempt_cannot_clear_predecessor_failure(
        tmp_path, monkeypatch):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    invalid_state = copy.deepcopy(state)
    invalid_state['gas_fractions']['Ar'] = .5
    first_entered = threading.Event()
    release_first = threading.Event()
    second_requested_transaction = threading.Event()
    original_state = provider._qualification_state
    original_run = provider._run_qualification
    outcomes = {}

    def held_snapshot(y, t):
        if threading.current_thread().name == 'invalid-first':
            first_entered.set()
            assert release_first.wait(5.)
        return original_state(y, t)

    def observe_transaction_request(*args, **kwargs):
        if threading.current_thread().name == 'valid-second':
            second_requested_transaction.set()
        return original_run(*args, **kwargs)

    def qualify(name, terminal_state):
        try:
            outcomes[name] = provider.qualify(terminal_state, 2.5)
        except Exception as error:
            outcomes[name] = error

    monkeypatch.setattr(provider, '_qualification_state', held_snapshot)
    monkeypatch.setattr(provider, '_run_qualification', observe_transaction_request)
    first = threading.Thread(
        target=qualify, args=('first', invalid_state), name='invalid-first')
    second = threading.Thread(
        target=qualify, args=('second', state), name='valid-second')
    first.start()
    assert first_entered.wait(5.)
    second.start()
    queued_before_release = second_requested_transaction.wait(1.)
    release_first.set()
    first.join(5.)
    second.join(5.)

    assert not first.is_alive() and not second.is_alive()
    assert 'normalize gas fractions to one' in str(outcomes['first'])
    assert 'latched FAILED' in str(outcomes['second'])
    assert not queued_before_release
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_null_optional_record_path_refuses_and_latches_inside_transaction(
        tmp_path, monkeypatch):
    provider, state = public_qualification_fixture(tmp_path, monkeypatch)
    state['record_path'] = None

    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError, match='invalid record_path'):
        provider.qualify(state, 2.5)

    assert provider._qualification_status == 'FAILED'
    assert provider.qualification_record_path is None
    assert not list((tmp_path / 'qualification').glob('*.json'))
    with pytest.raises(PlasmaStateError, match='latched FAILED'):
        provider.qualify(qualification_state(tmp_path, provider), 2.5)
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_working_condition_placeholders_are_terminal_fingerprint_fields(tmp_path):
    path, _ = fixture_provider_table(tmp_path / 'table')
    spec = json.loads((path / 'generation_spec.json').read_text())
    spec['envelopes']['electron_density'] = {
        'min': 1.e15, 'max': 1.e17, 'reference': 1.e16}
    spec['working_conditions']['electronDensity'] = '{electron_density}'

    reference = qualification_setup_sha256(
        spec, {'electron_density': 1.e16}, [1.])
    interior = qualification_setup_sha256(
        spec, {'electron_density': 2.e16}, [1.])

    assert interior == reference


def test_row_uses_one_normalized_query_for_cache_and_evaluation(tmp_path, monkeypatch):
    path, model = fixture_provider_table(tmp_path)
    provider = make_provider(path, model)
    storage = np.array([.75])
    query = storage[:1].reshape(())
    original = provider._table.row

    def edit_then_evaluate(u, composition, branch):
        storage[0] = 2.25
        return original(u, composition, branch)

    monkeypatch.setattr(provider._table, 'row', edit_then_evaluate)
    row = provider.row(query, {})

    assert row.u == .75
    assert provider._cache_identity_snapshot == (.75, ())


def test_row_owns_query_before_waiting_for_cache_lock(tmp_path, monkeypatch):
    path, model = fixture_provider_table(tmp_path)
    provider = make_provider(path, model)
    first_inside_table = threading.Event()
    release_first = threading.Event()
    second_waiting = threading.Event()
    original_row = provider._table.row
    original_lock = provider._cache_lock

    def paused_row(u, composition, branch):
        if u == 1.:
            first_inside_table.set()
            assert release_first.wait(5.)
        return original_row(u, composition, branch)

    class ObservedLock:
        def __enter__(self):
            if threading.current_thread().name == 'second-query':
                second_waiting.set()
            return original_lock.__enter__()

        def __exit__(self, *args):
            return original_lock.__exit__(*args)

    monkeypatch.setattr(provider._table, 'row', paused_row)
    monkeypatch.setattr(provider, '_cache_lock', ObservedLock())
    first = threading.Thread(target=lambda: provider.row(1., {}), name='first-query')
    first.start()
    assert first_inside_table.wait(5.)

    storage = np.array([.75])
    query = storage[:1].reshape(())
    returned = []
    second = threading.Thread(
        target=lambda: returned.append(provider.row(query, {})),
        name='second-query')
    second.start()
    assert second_waiting.wait(5.)
    storage[0] = 2.25
    release_first.set()
    first.join(5.)
    second.join(5.)

    assert not first.is_alive()
    assert not second.is_alive()
    assert returned[0].u == .75
    assert provider._cache_identity_snapshot == (.75, ())


def test_public_row_and_diagnostic_conversion_detach_read_only_values(tmp_path):
    path, model = fixture_provider_table(tmp_path)
    provider = make_provider(path, model)
    row = provider.row(1., {})

    row_values = row._as_dict()
    assert isinstance(row_values['swarm'], dict)
    assert row_values['k_ine'] is not row.k_ine
    row_values['swarm']['mean_energy_eV'] = -1.
    row_values['k_ine'][0] = -1.
    assert row.swarm['mean_energy_eV'] == 2.
    assert row.k_ine[0] > 0.

    read_only = np.array([1., 2.])
    read_only.setflags(write=False)
    diagnostic = {
        'transport': provider.transport(row),
        'metadata': provider.runtime_metadata,
        'nested': (MappingProxyType({'array': read_only}), 3.),
    }
    converted = provider.annotate_development_diagnostic(diagnostic)
    assert isinstance(converted['transport'], dict)
    assert isinstance(converted['metadata'], dict)
    assert isinstance(converted['nested'][0], dict)
    assert converted['nested'][0]['array'] is not read_only
    converted['nested'][0]['array'][0] = -1.
    assert read_only[0] == 1.


def test_public_diagnostic_conversion_refuses_unknown_values(tmp_path):
    path, model = fixture_provider_table(tmp_path)
    provider = make_provider(path, model)

    with pytest.raises(TypeError, match='diagnostic.*unsupported value type.*object'):
        provider.annotate_development_diagnostic({'unknown': object()})


def test_query_coordinates_normalize_extended_real_scalars(tmp_path):
    path, model = fixture_provider_table(tmp_path)
    provider = make_provider(path, model)
    expected = float(np.longdouble('1.25'))

    assert provider.row(np.longdouble('1.25'), {}).u == expected
    assert provider.row(np.array(np.longdouble('1.25')), {}).u == expected
    row = provider.row(1., {})
    assert provider.domain_check(row.u, row.composition) == [1.]


class _DiagnosticFlag(enum.Enum):
    READY = 'ready'


class _DiagnosticCode(enum.IntEnum):
    READY = 1


class _DiagnosticText(str, enum.Enum):
    READY = 'ready'


class _DiagnosticDict(dict):
    pass


class _DiagnosticTuple(tuple):
    pass


class _DiagnosticArray(np.ndarray):
    pass


class _DiagnosticFloat(float):
    pass


@pytest.mark.parametrize('value', [
    b'ok',
    np.longdouble('1.0000000000000000001'),
    np.clongdouble(1. + 2.j),
    np.str_('ok'),
    np.bytes_(b'ok'),
    np.datetime64('2026-10-25', 'D'),
    np.datetime64('2026-10-25', 'ns'),
    _DiagnosticFlag.READY,
    _DiagnosticCode.READY,
    _DiagnosticText.READY,
    Fraction(1, 3),
    Decimal('0.1'),
    _DiagnosticDict({'a': 1.}),
    _DiagnosticTuple((1.,)),
    np.array([1., 2.]).view(_DiagnosticArray),
    np.matrix([[1., 2.]]),
    np.ma.array([1., 99.], mask=[False, True]),
    np.array([np.longdouble('1.0000000000000000001')]),
    np.array([(np.array([1., 2.]), 2)], dtype=[('series', object), ('count', 'i4')]),
    np.array([((np.array([1., 2.]),), 2)],
             dtype=[('nested', [('series', object)]), ('count', 'i4')]),
    np.array([(1.25, 2)], dtype=[('energy', 'f8'), ('count', 'i4')]).view(np.recarray),
    np.array([(1.25, 2)], dtype=[('energy', 'f8'), ('count', 'i4')])[0],
    np.array([{'value': 1.}], dtype=object),
    _DiagnosticFloat(.1),
], ids=lambda value: type(value).__name__)
def test_public_diagnostic_conversion_refuses_every_non_allowlisted_type(
        tmp_path, value):
    path, model = fixture_provider_table(tmp_path)
    provider = make_provider(path, model)

    with pytest.raises(TypeError, match='diagnostic.*unsupported value type'):
        provider.annotate_development_diagnostic({'value': value})


def test_public_diagnostic_conversion_copies_numeric_structures(tmp_path):
    path, model = fixture_provider_table(tmp_path)
    provider = make_provider(path, model)
    record = np.array([(1.25, 2)], dtype=[('energy', 'f8'), ('count', 'i4')])
    source = np.array([1., 2.])
    source.setflags(write=False)

    converted = provider.annotate_development_diagnostic({
        'scalar': np.float64(1.25),
        'complex': np.complex128(1. + 2.j),
        'record': record,
        'nested': MappingProxyType({'array': source}),
    })

    assert converted['scalar'] == 1.25
    assert converted['complex'] == 1. + 2.j
    assert converted['record'].dtype == record.dtype
    assert not np.shares_memory(converted['record'], record)
    assert converted['nested']['array'].flags.writeable
    assert not np.shares_memory(converted['nested']['array'], source)


def test_provider_constructor_owns_configuration_before_table_io(
        tmp_path, monkeypatch):
    path, model = fixture_provider_table(tmp_path)
    reaction_map = {4: (1, 'ine')}
    original = EEDFTable.load

    def load_and_edit(*args, **kwargs):
        result = original(*args, **kwargs)
        reaction_map[4] = (0, 'sup')
        return result

    monkeypatch.setattr(EEDFTable, 'load', load_and_edit)
    provider = EEDFProvider(
        path, model, artifact_sha256=_file_hash(path / 'table.h5'),
        reaction_map=reaction_map, branch='branch_0')
    row = provider.row(1., {})

    assert provider.reaction_rate(4, row) == pytest.approx(row.k_ine[1])


def test_test_only_provider_accessors_are_private():
    for name in (
            'branch', 'reaction_map', 'empirical_laws', 'manifest',
            'cache_identity', 'cached_row'):
        assert not hasattr(EEDFProvider, name)
