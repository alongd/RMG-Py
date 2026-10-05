"""Opt-in real LoKI-B artifact and pure-Ar EEDF harness contract."""

import json
import os
from types import SimpleNamespace

import pytest

import rmgpy.rmg.input as rmg_input
from rmgpy.rmg.main import RMG
from rmgpy.tools.eedf.schema import file_hash
from real_eedf_harness import (
    DEVELOPMENT_BINDINGS,
    DEVELOPMENT_RETAINED_ENTRIES,
    development_argon_model,
    generate_real_artifact,
    load_real_spec,
    pure_argon_deck,
    reactor_readiness,
)
from rmgpy.solver.eedf_provider import DevelopmentWallBudgetExceeded
from runRealEEDFArgon import development_diagnostic, development_failure


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
        'scientific_status': 'DEVELOPMENT \u2014 UNQUALIFIED TABLE',
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
    assert failure['scientific_status'] == 'DEVELOPMENT \u2014 UNQUALIFIED TABLE'
    assert failure['energy_budget']['P_abs'] == 0.5
    assert failure['latest_energy_terms']['Q_inelastic'] == 0.125


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
    assert mole_fractions == pytest.approx({
        'Ar': 0.999998, 'Arp': 1.0e-6, 'e-': 1.0e-6, 'Ars': 0.0,
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
