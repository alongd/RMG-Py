"""Maxwellian comparator through the common reactor/provider energy path."""
import copy
import json
from pathlib import Path
import runpy

import numpy as np
import pytest

import rmgpy.constants as constants
import rmgpy.data.rmg as rmg_data_module
from rmgpy.exceptions import PlasmaStateError
from rmgpy.kinetics import EEDFChannel
from rmgpy.solver.eedf_provider import EEDFProvider
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.tools.eedf.artifact import _write_artifact
from rmgpy.tools.eedf.maxwellian import generate_maxwellian
from rmgpy.tools.eedf.schema import _content_hash, _interpolant_identity


@pytest.fixture(autouse=True)
def isolate_database():
    saved = rmg_data_module.database
    rmg_data_module.database = None
    try:
        yield
    finally:
        rmg_data_module.database = saved


def build_maxwellian_reactor(tmp_path, provider_kind='maxwellian-table'):
    helpers = runpy.run_path(str(Path(__file__).with_name('plasmaEEDFCoreTest.py')))
    source_data = runpy.run_path(str(Path(__file__).parents[1] / 'tools/eedf/maxwellianTest.py'))['source_data']
    manifest, rows = source_data((2., 2.02, 2.04))
    core, reaction = helpers['species_and_reaction']()
    reaction.kinetics = EEDFChannel(process=manifest['channel_map'][2]['description'],
                                    collision_set='fixture', side='ine')
    identity = {'library': 'Toy', 'index': 1, 'repr': repr(reaction.kinetics)}
    manifest['channel_map'][2].update(classification='A', reaction=identity)
    inputs = manifest['row_inputs']
    inputs.update(reactions=[identity], arm={'gases': ['Ar']},
                  P_Pa=helpers['PRESSURE'], interpolant=_interpolant_identity())
    manifest['floors']['absolute_power_share'] = 1e-4
    manifest['envelopes'] = {
        'Tg_K': {'min': inputs['Tg_K'], 'max': inputs['Tg_K'], 'reference': inputs['Tg_K']},
        'P_Pa': {'min': inputs['P_Pa'], 'max': inputs['P_Pa'], 'reference': inputs['P_Pa']},
        'electron_density': {'min': 0., 'max': 1e30, 'reference': 1e16},
    }
    inputs['envelopes'] = manifest['envelopes']
    manifest['fingerprint'] = _content_hash(inputs)
    # Synthetic LoKI-side rows and qualification metadata are test fixtures.
    # Both runs use this same source, map, grid, gas state, deck and closure.
    manifest['tolerances']['G1'] = {'rtol': 1e-6, 'atol': 1e-24}
    manifest.update(accepted=True, held_out=[{'u': .5}],
                    held_out_verdicts=[{'passed': True, 'checks': [{'passed': True}]}],
                    **{'T10 interpolation qualification': {'passed': True, 'verdict_count': 1}})
    for i, row in enumerate(rows):
        row['swarm'].update(mobility_N=(1 + .1*i)*1e24, diffusion_N=(2 + .2*i)*1e24)
        row['k_ine'] = np.array([8e-14, 4e-17, 1e-16]) * (1 + .1*i)
        row['channel_power'] = row['k_ine'] * np.array([.001, 11.548, 15.76])
        row['power_groups'] = {'field': row['swarm']['mobility_N'] * (row['EN_Td'] * 1e-21)**2}
    held = copy.deepcopy(rows[0])
    held.update(u=.5, EN_Td=float(np.exp(.5)))
    held['power_groups']['field'] = held['swarm']['mobility_N'] * (held['EN_Td'] * 1e-21)**2
    source = _write_artifact(tmp_path / 'source', manifest, {'branch_0': rows}, [held])
    path = (generate_maxwellian(source, tmp_path / 'comparator', artifact_sha256=source.name, branch_id='branch_0')
            if provider_kind == 'maxwellian-table' else source)
    config = {'provider': provider_kind, 'table': (str(path), path.name),
              'branch': 'maxwellian' if provider_kind == 'maxwellian-table' else 'branch_0', 'initial_reduced_field': (np.e, 'Td'), 'empirical_laws': {}}
    reactor = PlasmaReactor(
        (inputs['Tg_K'], 'K'), (inputs['P_Pa'], 'Pa'),
        {core[0]: 1e-6, core[1]: .999998, core[2]: 1e-6, core[3]: 0.},
        (inputs['Tg_K'], 'K'), termination=[], diffusion_length=(helpers['LAMBDA'], 'm'),
        ion_reduced_mobility=(1.535e-4, 'm^2/(V*s)'),
        wall_neutralization_products={'Ar+': 'Ar'}, thermo_source_assertions={'Ar+': 'ion'},
        electron_energy_balance={'absorbed_power': (.5, 'W'), 'chamber_volume': (1e-4, 'm^3'),
                                 'sheath': 'floating_wall'}, electron_kinetics=config)
    reactor.initialize_model(core, [reaction], [], [])
    return reactor, core, [reaction]


def test_maxwellian_rates_energy_equation_and_shared_provider(tmp_path):
    reactor, core, _ = build_maxwellian_reactor(tmp_path / 'M')
    loki, _, _ = build_maxwellian_reactor(tmp_path / 'E', provider_kind='loki-table')
    assert type(reactor) is type(loki) is PlasmaReactor
    assert type(reactor.eedf_provider) is type(loki.eedf_provider) is EEDFProvider
    assert reactor.eedf_provider.row.__func__ is loki.eedf_provider.row.__func__
    assert reactor.eedf_provider.reaction_rate.__func__ is loki.eedf_provider.reaction_rate.__func__
    assert reactor.maxwellian_comparator and not loki.maxwellian_comparator
    assert reactor.eedf_provider.runtime_metadata['channel_map'] == loki.eedf_provider.runtime_metadata['channel_map']
    assert reactor.eedf_provider.axis('u') == pytest.approx(loki.eedf_provider.axis('u'))
    assert not np.array_equal(reactor.eedf_row.k_ine, loki.eedf_row.k_ine)
    y = np.array(reactor.y0)
    delta = reactor.residual(0., y, np.zeros_like(y))[0]
    row = reactor.eedf_row
    terms = reactor.electron_energy_terms
    stored = constants.Na * constants.e * (
        row.swarm['mean_energy_eV'] * delta[reactor.electron_index]
        + y[reactor.electron_index] * row.denergy_du * delta[reactor.te_index])
    losses = sum(terms[key] for key in ('Q_inelastic', 'Q_elastic', 'Q_wall_electron', 'Q_wall_ion', 'Q_flow'))
    assert stored == pytest.approx(terms['P_abs'] - losses, rel=2e-8)
    volume = reactor.compute_volume(y)
    expected_rate = (row.k_ine[2] * constants.Na * y[reactor.electron_index] / volume
                     * y[reactor.species_index[core[1]]] / volume)
    assert reactor.core_reaction_rates[0] == pytest.approx(expected_rate, rel=1e-12)
    assert terms['joule_power'] is None and terms['power_mismatch_fraction'] is None
    assert terms['A6b-source']['status'] == terms['A6b-runtime']['status'] == 'NOT_APPLICABLE'
    assert terms['T10 interpolation qualification']['status'] == 'NOT_APPLICABLE'
    assert reactor.energy_budget['EN_Td'] is None
    assert reactor.energy_budget['reference_EN_Td'] == np.e
    assert reactor.energy_budget['A6a_passed']
    assert 'MODEL-FORM COMPARATOR: Maxwellian-table energy-balance diagnostic' in str(terms)
    assert 'EEDF energy-balance engineering intermediate' not in str(terms)


def test_maxwellian_terminal_keeps_physical_closure_without_qualification(tmp_path, monkeypatch):
    reactor, _, _ = build_maxwellian_reactor(tmp_path)
    monkeypatch.setattr(reactor.eedf_provider, 'qualify', lambda *a: pytest.fail('M production qualification'))
    reactor.energy_budget.update(A6a_relative=0., A6a_steady_relative=0.)
    reactor.validate_terminal_state()
    assert reactor.terminal_admission is None
    assert reactor.eedf_provider.qualification_record_path is None
    reactor.energy_budget['A6a_relative'] = .011
    with pytest.raises(PlasmaStateError, match='A6a'):
        reactor.validate_terminal_state()


def test_maxwellian_profile_and_manifests_have_literal_comparator_provenance(tmp_path):
    reactor, core, _ = build_maxwellian_reactor(tmp_path)
    helpers = runpy.run_path(str(Path(__file__).parents[1] / 'diagnosticModeTest.py'))
    reactor.snapshots = [[0., 1.] + [0.] * len(core)]
    profile = Path(helpers['_write_profile'](tmp_path, reactor, core))
    manifest_path = profile.with_suffix('.eedf.json')
    for path in (profile, manifest_path):
        text = path.read_text()
        assert 'MODEL-FORM COMPARATOR: Maxwellian-table energy-balance diagnostic' in text
        assert 'EEDF energy-balance engineering intermediate' not in text
    record = json.loads(manifest_path.read_text())
    assert record['provider'] == 'maxwellian-table'
    assert record['comparator']['numerical_validation']['passed'] is True
    wall = reactor.electronegative_wall_manifest()
    assert 'MODEL-FORM COMPARATOR: Maxwellian-table energy-balance diagnostic' in str(wall)
    assert 'EEDF energy-balance engineering intermediate' not in str(wall)
