"""Opt-in real LoKI binding of the full resolved argon electron library."""
import json
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import rmgpy.constants as constants
import rmgpy.data.rmg as data_module
from rmgpy.data.thermo import ThermoDatabase
from rmgpy.exceptions import PlasmaStateError
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.solver.eedf_provider import _development_unqualified_table_route
from rmgpy.tools.eedf.channel_map import _load_library
from rmgpy.tools.eedf.schema import _file_hash


@pytest.mark.parametrize('fraction,refused', [(1.e-7, False), (3.e-6, True)])
def test_real_argon_library_binds_and_advances(tmp_path, monkeypatch, fraction, refused):
    required = ('EEDF_ARGON_ARTIFACT', 'EEDF_DATABASE_ROOT')
    if any(not os.environ.get(key) for key in required):
        pytest.skip('requires a real, unqualified argon diagnostic table and database')
    artifact = Path(os.environ['EEDF_ARGON_ARTIFACT'])
    database = Path(os.environ['EEDF_DATABASE_ROOT'])
    manifest = json.loads((artifact / 'manifest.json').read_text())
    assert manifest['accepted'] is False
    assert manifest['excited_fraction_max'] == 1.e-5
    assert manifest['argon_elimination_electron_density_max_m3'] == 1.e18
    library, _ = _load_library(database, tmp_path)
    reactions = library.get_library_reactions()
    assert len(reactions) == 573
    by_label = {s.label: s for r in reactions for s in r.reactants + r.products}
    assert len(by_label) == 7
    species = list(by_label.values())
    thermo = ThermoDatabase()
    thermo.load_libraries(str(database / 'input/thermo/libraries'),
                          ['PlasmaExcitedNeutralThermo', 'PlasmaThermo', 'primaryThermoLibrary'])
    monkeypatch.setattr(data_module, 'database', SimpleNamespace(thermo=thermo, solvation=None))
    for s in species:
        if s.label != 'e-':
            s.thermo = thermo.get_thermo_data(s)
    for i, r in enumerate(reactions):
        r.index = i + 1
    initial = {s: fraction for s in species}
    initial[by_label['e-']] = initial[by_label['Arp']] = 1.e-7
    initial[by_label['Ar']] = 1. - 4. * fraction - 2.e-7
    tg, pressure = manifest['row_inputs']['Tg_K'], manifest['row_inputs']['P_Pa']
    reactor = PlasmaReactor(
        (tg, 'K'), (pressure, 'Pa'), initial, (tg, 'K'), termination=[],
        resonance_radiation={'model': 'complete', 'geometryModel': 'infiniteCylinder'},
        radiation_chamber_geometry={'shape': 'cylinder', 'radius': .05, 'length': .30},
        diffusion_length=(.02, 'm'), ion_reduced_mobility=(1.535e-4, 'm^2/(V*s)'),
        wall_single_bath_approximation=True,
        wall_neutralization_products={'Arp': 'Ar'},
        electron_energy_balance={'absorbed_power': (.5, 'W'), 'chamber_volume': (1.e-4, 'm^3'),
                                 'sheath': 'floating_wall'},
        electron_kinetics={'provider': 'loki-table',
                           'table': (str(artifact), _file_hash(artifact / 'table.h5')),
                           'branch': 'branch_0', 'initial_reduced_field': (17., 'Td'),
                           'empirical_laws': {}})
    with _development_unqualified_table_route():
        if refused:
            # Each individual 4s fraction is below 1e-5; only their sum exceeds it.
            with pytest.raises(PlasmaStateError, match='excited_fraction_max.*observed'):
                reactor.initialize_model(species, reactions, [], [])
            return
        reactor.initialize_model(species, reactions, [], [])
    reactor.residual(0., reactor.y0.copy(), np.zeros_like(reactor.y0))
    channels = {c['description']: i for i, c in enumerate(manifest['channel_map'])}
    counts = {'ordinary_ine': 0, 'sup': 0, 'branches': 0}
    for i, reaction in enumerate(reactions):
        marker = reaction.kinetics
        column = channels[marker.process]
        raw = getattr(reactor.eedf_row, 'k_' + marker.side)[column]
        weight = 1. if marker.branch_product is None else marker.branch_yield
        assert reactor.kf[i] / constants.Na == pytest.approx(raw * weight, rel=2.e-12, abs=1.e-40)
        counts['branches' if marker.branch_product is not None else
               'sup' if marker.side == 'sup' else 'ordinary_ine'] += 1
    assert counts == {'ordinary_ine': 15, 'sup': 10, 'branches': 548}
    assert len(reactor.eedf_reaction_map) == 573
    assert all(reactor.core_reaction_rates[i] > 0. for i, r in enumerate(reactions)
               if r.kinetics.branch_product is None)
    initial_y = reactor.y0.copy()
    reactor.advance(1.e-13)
    assert reactor.t >= 1.e-13
    assert np.all(np.isfinite(reactor.y))
    assert not np.array_equal(reactor.y, initial_y)
    charge = np.array([s.get_net_charge() for s in species])
    assert np.dot(charge, reactor.y[:len(species)]) == pytest.approx(0., abs=1.e-18)
    result = {
        'outcome': 'PASS', 'scientific_status': 'DEVELOPMENT ONLY - TABLE QUALIFICATION PENDING',
        'artifact': str(artifact), 'artifact_sha256': manifest['artifact_sha256'],
        'accepted': manifest['accepted'], 'species': len(species), 'reaction_bindings': len(reactions),
        'rate_comparisons': len(reactions), 'comparison_failures': 0, 'owner_counts': counts,
        'advanced_to_s': float(reactor.t), 'electron_density_m3': float(
            reactor.y[reactor.species_index[by_label['e-']]] * constants.Na / reactor.compute_volume(reactor.y)),
        'excited_fraction': float(sum(reactor.y[reactor.species_index[by_label[label]]]
                                      for label in ('Ars', 'Ar(1s4)', 'Ar(1s3)', 'Ar(1s2)')) /
                                  reactor.y[reactor.species_index[by_label['Ar']]]),
    }
    output = os.environ.get('EEDF_ARGON_REACTOR_REPORT')
    if output:
        Path(output).write_text(json.dumps(result, indent=2) + '\n')

