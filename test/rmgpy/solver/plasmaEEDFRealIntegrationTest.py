"""Opt-in real LoKI-B artifact and pure-Ar EEDF harness contract."""

import json
import os
from types import MappingProxyType, SimpleNamespace

import pytest

import rmgpy.rmg.input as rmg_input
from rmgpy.rmg.main import RMG
from rmgpy.exceptions import PlasmaStateError
from rmgpy.tools.eedf.schema import file_hash
from real_eedf_harness import (
    DEVELOPMENT_BINDINGS,
    DEVELOPMENT_SWEEP_TOLERANCES,
    DEVELOPMENT_RETAINED_ENTRIES,
    compare_development_endpoints,
    development_argon_model,
    generate_real_artifact,
    load_real_spec,
    power_consistent_initialization,
    pure_argon_deck,
    reactor_readiness,
)
from rmgpy.solver.eedf_provider import DevelopmentWallBudgetExceeded
from rmgpy.solver.electronegative import manifest_values
import runRealEEDFArgon as real_eedf_runner
from runRealEEDFArgon import (
    development_diagnostic,
    development_failure,
    negative_fixture_result,
    select_restart_checkpoint,
)
from summarizeRealEEDFArgon import build_summary


def test_manifest_values_converts_read_only_nested_mappings():
    value = MappingProxyType({'endpoint': MappingProxyType({'value': 1.0})})
    assert json.dumps(manifest_values(value), allow_nan=False) == (
        '{"endpoint": {"value": 1.0}}')


def test_development_argon_channel_bindings_are_runtime_copies_only():
    sources = [
        SimpleNamespace(library='PlasmaArgon', entry=SimpleNamespace(index=index),
                        kinetics='original-{0}'.format(index))
        for index in sorted(DEVELOPMENT_RETAINED_ENTRIES)
    ]
    sources.append(SimpleNamespace(
        library='PlasmaRadiativeRecombination', entry=SimpleNamespace(index=1),
        kinetics='not in development deck'))

    reactions = development_argon_model(sources)

    assert len(reactions) == len(DEVELOPMENT_RETAINED_ENTRIES)
    assert all(reaction is not source for reaction, source in zip(reactions, sources))
    assert {reaction.entry.index for reaction in reactions} == DEVELOPMENT_RETAINED_ENTRIES
    for reaction in reactions:
        if reaction.entry.index in DEVELOPMENT_BINDINGS:
            assert reaction.kinetics.process == DEVELOPMENT_BINDINGS[reaction.entry.index]
            assert reaction.kinetics.side == 'ine'
        else:
            assert reaction.kinetics == 'original-{0}'.format(reaction.entry.index)
    assert all(source.kinetics == 'original-{0}'.format(source.entry.index)
               for source in sources[:-1])


def test_every_development_harness_diagnostic_carries_warning_and_blockers():
    blockers = ['artifact manifest is not accepted',
                'branch_0: uncertified: unseeded scans']
    assert development_diagnostic(blockers) == {
        'scientific_status': 'DEVELOPMENT ONLY \u2014 TABLE QUALIFICATION FAILED',
        'artifact_blockers': blockers,
        'export_allowed': False,
        'qualification_allowed': False,
    }


def test_wall_budget_failure_carries_stop_reason_and_last_progress():
    progress = {'t_s': 1.e-9, 'dt_s': 1.e-15, 'step_count': 12}
    reactor = SimpleNamespace(
        development_wall_budget_hit=True,
        development_last_progress=progress,
        energy_budget={'P_abs': 0.5},
        electron_energy_terms={'Q_inelastic': 0.125},
        y=None,
        y0=None,
    )

    failure = development_failure(
        reactor, ['artifact manifest is not accepted'],
        DevelopmentWallBudgetExceeded('budget exhausted'))

    assert failure['stop_reason'] == 'wall_budget'
    assert failure['last_progress'] == progress
    assert failure['exception_type'] == 'DevelopmentWallBudgetExceeded'
    assert failure['scientific_status'] == (
        'DEVELOPMENT ONLY \u2014 TABLE QUALIFICATION FAILED')
    assert failure['energy_budget']['P_abs'] == 0.5
    assert failure['latest_energy_terms']['Q_inelastic'] == 0.125


def test_negative_fixture_and_restart_checkpoint_keep_named_provenance():
    reactor = SimpleNamespace(
        y=[1.0, 0.0, 1.e-6], y0=[1.0, 0.0, 1.e-6],
        electron_kinetics={
            'provider': 'loki-table',
            'table': ('table.h5', 'a' * 64),
        },
    )
    message = (
        'accepted EEDF state is outside the table domain: u=1.0 outside [1.1, 4.1]')

    wrong_type = negative_fixture_result(reactor, RuntimeError(message))
    negative = negative_fixture_result(reactor, PlasmaStateError(message))
    checkpoint = select_restart_checkpoint([
        {'time_s': 0.0},
        {'time_s': 0.4, 'scientific_status': (
            'DEVELOPMENT ONLY \u2014 TABLE QUALIFICATION FAILED')},
        {'time_s': 1.0},
    ])

    assert negative['fixture'] == 'OVER-IONISED INITIAL-STATE DOMAIN-GUARD TEST'
    assert wrong_type['outcome'] == 'FAIL'
    assert wrong_type['checks']['plasma_state_error'] is False
    assert negative['outcome'] == 'PASS'
    assert all(negative['checks'].values())
    assert checkpoint['time_s'] == 0.4
    assert checkpoint['checkpoint_role'] == (
        'intermediate accepted state for restart test')
    assert checkpoint['scientific_status'] == (
        'DEVELOPMENT ONLY \u2014 TABLE QUALIFICATION FAILED')


@pytest.fixture(scope='session')
def real_pure_argon_artifact(tmp_path_factory):
    spec_path = os.environ.get('EEDF_REAL_SPEC')
    if not spec_path:
        pytest.skip('EEDF_REAL_SPEC is not set; it must name the pinned pure-Ar generation spec')
    try:
        load_real_spec(spec_path)
    except FileNotFoundError as error:
        pytest.skip(str(error))
    work = tmp_path_factory.mktemp('real-pure-ar-eedf')
    return generate_real_artifact(
        spec_path, work,
        command_line=['pytest', 'plasmaEEDFRealIntegrationTest.py'])


def test_full_real_artifact_is_unaccepted_and_fingerprinted(real_pure_argon_artifact):
    artifact, manifest = real_pure_argon_artifact
    assert manifest['accepted'] is False
    assert len([item for item in manifest['held_out_verdicts']
                if not item['passed']]) == 8
    assert file_hash(artifact / 'table.h5') == manifest['artifact_sha256']
    assert manifest['branches'] == ['branch_0']
    assert manifest['branch_detection'][0]['agreement'] is True
    assert manifest['row_inputs']['Tg_K'] == 298.15
    assert manifest['row_inputs']['P_Pa'] == pytest.approx(5.0 * 133.322368)
    assert reactor_readiness(manifest) == [
        'artifact manifest is not accepted',
        '8 of 188 held-out verdicts failed',
        'branch_0: uncertified: unseeded scans',
        "artifact has no reaction-owned channels; channel classifications are {'B': 39}",
    ]


def test_run_deck_is_reproducible_and_names_current_blockers(
        real_pure_argon_artifact, tmp_path):
    artifact, manifest = real_pure_argon_artifact
    initialization = power_consistent_initialization(artifact, manifest)
    deck = pure_argon_deck(artifact, manifest)
    compile(deck, 'input.py', 'exec')
    path = tmp_path / 'input.py'
    path.write_text(deck)

    rmg = RMG()
    rmg_input.read_input_file(str(path), rmg)
    reactor = rmg.reaction_systems[0]
    mole_fractions = {
        species.label: fraction
        for species, fraction in reactor.initial_mole_fractions.items()
    }

    assert "temperature=(298.15, 'K')" in deck
    assert "pressure=(666.61184, 'Pa')" in deck
    assert "'diffusionLength': (20.3, 'mm')" in deck
    assert "'absorbedPower': (0.5, 'W')" in deck
    assert str(artifact / 'table.h5') in deck
    assert manifest['artifact_sha256'] in deck
    assert reactor.electron_kinetics['table'][0] == str(artifact / 'table.h5')
    assert initialization['prescribed_absorbed_power_W'] == 0.5
    assert initialization['chamber_volume_m3'] == pytest.approx(2.356194490192345e-3)
    assert initialization['initial_reduced_field_Td'] == 17.0
    assert initialization['gas_density_m^-3'] == pytest.approx(1.6194013103625856e23)
    assert initialization['reference_electron_density_m^-3'] == pytest.approx(
        2.7715502016522555e13, rel=2.e-6)
    assert initialization['electron_mole_fraction'] == pytest.approx(
        1.711465992158105e-10, rel=2.e-6)
    assert initialization['reactor_reference_volume_m3'] == pytest.approx(
        3.7187455691997435, rel=2.e-6)
    assert initialization['reactor_absorbed_power_W'] == pytest.approx(
        789.1423192522975, rel=2.e-6)
    assert initialization['table_field_power_coefficient_eV_m3_s^-1'] == pytest.approx(
        2.9510139153579735e-16, rel=2.e-6)
    assert initialization['table_A6b_relative'] == pytest.approx(
        7.447369314637744e-6)
    assert initialization['table_A6b_tolerance'] == 1.e-6
    assert initialization['table_A6b_outcome'] == 'FAIL'
    assert initialization['over_ionised_field_power_W'] == pytest.approx(
        4.610238207728393e6, rel=2.e-4)
    assert mole_fractions == pytest.approx({
        'Ar': initialization['argon_mole_fraction'],
        'Arp': initialization['electron_mole_fraction'],
        'e-': initialization['electron_mole_fraction'],
        'Ars': 0.0,
    })
    assert sum(
        species.get_net_charge() * fraction
        for species, fraction in reactor.initial_mole_fractions.items()
    ) == pytest.approx(0.0)

    status = {
        'artifact': str(artifact),
        'artifact_sha256': manifest['artifact_sha256'],
        'accepted': manifest['accepted'],
        'branches': manifest['branches'],
        'branch_certification': manifest['branch_certification'],
        'blockers': reactor_readiness(manifest),
    }
    (tmp_path / 'harness-status.json').write_text(json.dumps(status, indent=2, sort_keys=True))
    assert status['blockers']


def test_sweep_comparison_declares_every_required_metric_before_execution():
    endpoint = {
        'converged_to_steady_state': True,
        'electron_density_m^-3': 2.4e14,
        'EN_Td': 5.77,
        'mean_energy_eV': 5.17,
        'composition': {'Ar': .9999999, 'Arp': 1.e-7},
        'wall_flux_mol_s': {'e-': 1.e-9, 'Arp': 1.e-9},
        'electron_power_partition_W': {
            'Q_inelastic': 310., 'Q_elastic': 474.,
            'Q_wall_electron': 1., 'Q_wall_ion': 3.,
        },
        'A6a_relative': 1.e-12,
        'A6b_relative': 3.12e-6,
        'convergence_time_s': 0.0165,
        'transient_extrema': {
            'electron_density_m^-3': {'min': 2.8e13, 'max': 2.4e14},
            'EN_Td': {'min': 5.77, 'max': 17.0},
        },
    }
    arms = [dict(endpoint, arm=name) for name in ('0.1x', '1x', '10x', '20x')]

    comparison = compare_development_endpoints(arms)

    assert comparison['tolerances'] == DEVELOPMENT_SWEEP_TOLERANCES
    assert comparison['all_endpoints_agree'] is True
    assert set(comparison['metrics']) == {
        'electron_density_m^-3', 'EN_Td', 'mean_energy_eV',
        'composition.Ar', 'composition.Arp',
        'wall_flux_mol_s.e-', 'wall_flux_mol_s.Arp',
        'electron_power_partition_W.Q_inelastic',
        'electron_power_partition_W.Q_elastic',
        'electron_power_partition_W.Q_wall_electron',
        'electron_power_partition_W.Q_wall_ion',
        'A6a_relative', 'A6b_relative',
    }
    assert comparison['convergence_time_comparison']['agreement_required'] is False
    assert comparison['convergence_time_comparison']['arms']['1x']['value_s'] == 0.0165
    assert comparison['transient_extrema_comparison']['agreement_required'] is False
    assert comparison['transient_extrema_comparison']['metrics']['EN_Td']['20x'] == {
        'min': 5.77, 'max': 17.0}

    restart = dict(endpoint, arm='restart')
    summary = build_summary(arms, restart)
    assert summary['scientific_status'] == (
        'DEVELOPMENT ONLY \u2014 TABLE QUALIFICATION FAILED')
    assert summary['export_allowed'] is False
    assert summary['qualification_allowed'] is False
    assert summary['different_branches_detected'] is False


def test_nonconverged_results_cannot_be_candidates_or_enter_summary():
    reactor = SimpleNamespace(steady_state_reached=False)
    with pytest.raises(RuntimeError, match='did not reach steady state'):
        real_eedf_runner.candidate_endpoint_classification(reactor)

    endpoint = {
        'converged_to_steady_state': True,
        'electron_density_m^-3': 2.4e14,
        'EN_Td': 5.77,
        'mean_energy_eV': 5.17,
        'composition': {'Ar': .9999999, 'Arp': 1.e-7},
        'wall_flux_mol_s': {'e-': 1.e-9, 'Arp': 1.e-9},
        'electron_power_partition_W': {
            'Q_inelastic': 310., 'Q_elastic': 474.,
            'Q_wall_electron': 1., 'Q_wall_ion': 3.,
        },
        'A6a_relative': 1.e-12,
        'A6b_relative': 3.12e-6,
        'convergence_time_s': 0.0165,
        'transient_extrema': {},
    }
    arms = [dict(endpoint, arm=name) for name in ('0.1x', '1x', '10x', '20x')]
    restart = dict(endpoint, arm='restart')
    arms[0]['converged_to_steady_state'] = False
    with pytest.raises(ValueError, match='0.1x.*did not converge'):
        build_summary(arms, restart)

    arms[0]['converged_to_steady_state'] = True
    restart['converged_to_steady_state'] = False
    with pytest.raises(ValueError, match='restart.*did not converge'):
        build_summary(arms, restart)
