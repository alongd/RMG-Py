"""Engine-side acceptance tests for the PlasmaArgonDimer database library."""

import os
import re
import shutil
from collections import OrderedDict

import pytest

from rmgpy import constants, settings
from rmgpy.chemkin import (load_chemkin_file, save_chemkin_file,
                           save_species_dictionary)
from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.thermo import ThermoDatabase
from rmgpy.rmg.model import ReactionModel
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.transport import TransportData
from rmgpy.yaml_cantera2 import save_cantera_model


DATABASE_INPUT = settings['database.directory']
LIBRARY = 'PlasmaArgonDimer'

if not os.path.isdir(os.path.join(DATABASE_INPUT, 'kinetics', 'libraries', LIBRARY)):
    pytest.skip('kinetics library {0!r} is not present in the configured database ({1}); this '
                'suite exercises the real library and does not inject a substitute.'.format(
                    LIBRARY, DATABASE_INPUT), allow_module_level=True)


@pytest.fixture(scope='module')
def library():
    db = KineticsDatabase()
    db.load_libraries(os.path.join(DATABASE_INPUT, 'kinetics', 'libraries'), libraries=[LIBRARY])
    return db.libraries[LIBRARY]


@pytest.fixture(scope='module')
def rmg_database():
    db = RMGDatabase()
    db.load(
        DATABASE_INPUT,
        thermo_libraries=['primaryThermoLibrary', 'PlasmaThermo',
                          'PlasmaExcitedNeutralThermo', 'electrocatThermo'],
        reaction_libraries=[LIBRARY],
        kinetics_families=[],
        kinetics_depositories=[],
        solvation=False,
        surface=False,
        testing=True,
    )
    return db


def _mechanism(rmg_database):
    reactions = rmg_database.kinetics.libraries[LIBRARY].get_library_reactions()
    species = OrderedDict()
    for reaction in reactions:
        for spc in reaction.reactants + reaction.products:
            species.setdefault(spc.label, spc)
    for index, spc in enumerate(species.values(), start=1):
        spc.index = index
        spc.thermo = rmg_database.thermo.get_thermo_data(spc)
        if spc.thermo.__class__.__name__ != 'NASA':
            spc.thermo = spc.thermo.to_nasa(100, 5000, 1000)
        spc.transport_data = TransportData(
            shapeIndex=0 if len(spc.molecule[0].atoms) == 1 else 1,
            sigma=(3.0, 'angstrom'), epsilon=(100.0, 'K'),
            dipoleMoment=(0.0, 'De'), polarizability=(0.0, 'angstrom^3'),
            rotrelaxcollnum=0,
        )
    return list(species.values()), reactions


def _side_key(side):
    return sorted(spc.label for spc in side for _ in range(1))


def _find_reaction(reactions, reactants, products):
    for reaction in reactions:
        if _side_key(reaction.reactants) == sorted(reactants) and _side_key(reaction.products) == sorted(products):
            return reaction
    raise AssertionError((reactants, products, [str(r) for r in reactions]))


def test_real_library_loads_with_three_balanced_entries_and_explicit_electrons(library):
    assert len(library.entries) == 3
    def key(reaction):
        return (tuple(sorted(species.label for species in reaction.reactants)),
                tuple(sorted(species.label for species in reaction.products)))
    reactions = {key(reaction): reaction for reaction in library.get_library_reactions()}
    assert set(reactions) == {
        (('Ar', 'Ar', 'Arp'), ('Ar', 'Ar2p')),
        (('Ar2p', 'e-'), ('Ar', 'Ars')),
        (('Ar2p', 'e-'), ('Ar', 'Arp', 'e-')),
    }
    for reaction in reactions.values():
        assert reaction.is_balanced()
        assert reaction.electrons == 0
    dr = reactions[(('Ar2p', 'e-'), ('Ar', 'Ars'))]
    dei = reactions[(('Ar2p', 'e-'), ('Ar', 'Arp', 'e-'))]
    assert sum(species.is_electron() for species in dr.reactants) == 1
    assert sum(species.is_electron() for species in dr.products) == 0
    assert sum(species.is_electron() for species in dei.reactants) == 1
    assert sum(species.is_electron() for species in dei.products) == 1


def test_dr_unit_gate_is_preserved_by_loaded_engine(library):
    kinetics = library.entries[1].data
    expected = {300: 9.100000e-7, 8500: 1.183422e-7,
                10444: 1.043702e-7, 13462: 8.939815e-8}
    for te, source in expected.items():
        assert kinetics.get_rate_coefficient_two_temp(298.15, te) == pytest.approx(
            source * 1e-6 * constants.Na, rel=1e-6)


def test_dr_stores_the_shiu_biondi_te_law(library):
    kinetics = library.entries[1].data
    assert kinetics.n.value_si == pytest.approx(-0.61)
    assert kinetics.T0.value_si == pytest.approx(300.0)
    for te in (300.0, 8500.0, 13462.0):
        expected = 9.1e-7 * (300.0 / te) ** 0.61 * 1e-6 * constants.Na
        assert kinetics.get_rate_coefficient_two_temp(298.15, te) == pytest.approx(expected, rel=1e-12)


def test_loader_refuses_duplicate_explicit_and_metadata_electron(tmp_path):
    source = os.path.join(DATABASE_INPUT, 'kinetics', 'libraries', LIBRARY)
    target = tmp_path / LIBRARY
    shutil.copytree(source, target)
    reaction_file = target / 'reactions.py'
    content = reaction_file.read_text()
    content = content.replace("Tmax=(8500, 'K')),", "Tmax=(8500, 'K'), electrons=-1),", 1)
    reaction_file.write_text(content)
    db = KineticsDatabase()
    with pytest.raises(Exception, match='balanced|electron') as excinfo:
        db.load_libraries(str(tmp_path), libraries=[LIBRARY])
    message = str(excinfo.value)
    print('duplicate-electron refusal:', message)
    assert 'electron' in message.lower() or 'balanced' in message.lower()


def test_termolecular_reverse_fit_can_be_evaluated_at_1000_k(library):
    reaction = next(r for r in library.get_library_reactions() if r.reversible)
    thermo = ThermoDatabase()
    thermo.load_libraries(os.path.join(DATABASE_INPUT, 'thermo', 'libraries'),
                          libraries=['PlasmaThermo', 'primaryThermoLibrary'])
    thermo.load_groups(os.path.join(DATABASE_INPUT, 'thermo', 'groups'))
    for species in reaction.reactants + reaction.products:
        species.thermo = thermo.get_thermo_data(species)
    # The collision-limit diagnostic requests a reverse fit on its own grid.  The
    # source Tmax remains 300 K; no generic range guard may turn the 1000 K call into
    # an exception.
    reverse = reaction.generate_reverse_rate_coefficient()
    assert reverse.get_rate_coefficient(1000.0) > 0.0


def test_conversion_reverse_is_kf_over_keq_at_298_and_1000_k(rmg_database):
    species, reactions = _mechanism(rmg_database)
    conversion = _find_reaction(reactions, ['Arp', 'Ar', 'Ar'], ['Ar2p', 'Ar'])
    for temperature in (298.15, 1000.0):
        imf = {spc: 0.0 for spc in species}
        imf[next(s for s in species if s.label == 'Ar')] = 0.98
        imf[next(s for s in species if s.label == 'Arp')] = 0.01
        imf[next(s for s in species if s.label == 'Ar2p')] = 0.01
        imf[next(s for s in species if s.label == 'e-')] = 1e-12
        reactor = PlasmaReactor((temperature, 'K'), (5.0, 'torr'), imf,
                                (temperature, 'K'), n_sims=1, termination=[])
        reactor.initialize_model(species, reactions, [], [])
        index = reactor.reaction_index[conversion]
        assert reactor.kb[index] == pytest.approx(
            reactor.kf[index] / reactor.Keq[index], rel=1e-6)
        assert conversion.get_equilibrium_constant(temperature) == pytest.approx(
            reactor.Keq[index], rel=1e-12)


def test_rmg_database_wall_less_reactor_and_export_round_trips(rmg_database, tmp_path):
    species, reactions = _mechanism(rmg_database)
    ar = next(s for s in species if s.label == 'Ar')
    electron = next(s for s in species if s.label == 'e-')
    arp = next(s for s in species if s.label == 'Arp')
    ar2p = next(s for s in species if s.label == 'Ar2p')
    imf = {s: 0.0 for s in species}
    imf.update({ar: 0.9999997, electron: 1e-7, arp: 1e-7, ar2p: 1e-7})
    reactor = PlasmaReactor((298.15, 'K'), (5.0, 'torr'), imf,
                            (10444.0, 'K'), n_sims=1, termination=[])
    reactor.initialize_model(species, reactions, [], [])
    reactor.advance(1e-8)
    assert reactor.t == pytest.approx(1e-8)
    assert ar2p.thermo.label == '[Ar2p]'

    model = ReactionModel(species=species, reactions=reactions)
    chemkin = tmp_path / 'chem.inp'
    dictionary = tmp_path / 'dictionary.txt'
    cantera = tmp_path / 'chem.yaml'
    save_species_dictionary(str(dictionary), species)
    save_chemkin_file(str(chemkin), species, reactions)
    save_cantera_model(model, str(cantera))

    _, loaded_reactions = load_chemkin_file(str(chemkin), str(dictionary))
    loaded_conversion = _find_reaction(loaded_reactions, ['Arp', 'Ar', 'Ar'], ['Ar2p', 'Ar'])
    loaded_dr = _find_reaction(loaded_reactions, ['Ar2p', 'e-'], ['Ars', 'Ar'])
    assert loaded_conversion.reversible
    assert _side_key(loaded_conversion.reactants).count('Ar') == 2
    assert _side_key(loaded_conversion.products).count('Ar') == 1
    assert _side_key(loaded_dr.reactants).count('e-') == 1

    import cantera as ct
    gas = ct.Solution(str(cantera))
    cantera_conversion = next(r for r in gas.reactions()
                              if 'Arp' in r.equation and 'Ar2p' in r.equation)
    cantera_dr = next(r for r in gas.reactions()
                      if 'Ar2p' in r.equation and 'e-' in r.equation and 'Ars' in r.equation)
    assert cantera_conversion.reversible
    assert _equation_side_key(cantera_conversion.equation, left=True).count('Ar') == 2
    assert _equation_side_key(cantera_conversion.equation, left=False).count('Ar') == 1
    assert _equation_side_key(cantera_dr.equation, left=True).count('e-') == 1


def _equation_side_key(equation, left):
    sides = re.split(r'<=>|=>', equation)
    side = sides[0] if left else sides[1]
    return sorted(token.strip().split('(')[0] for token in side.split('+'))


def test_collision_limit_check_does_not_abort(rmg_database):
    species, reactions = _mechanism(rmg_database)
    conversion = _find_reaction(reactions, ['Arp', 'Ar', 'Ar'], ['Ar2p', 'Ar'])
    for temperature in (298.15, 1000.0):
        result = conversion.check_collision_limit_violation(
            temperature, temperature, 101325.0, 101325.0)
        assert isinstance(result, list)


def test_ar2p_thermo_wilhoit_and_nasa_conversion(rmg_database):
    data = rmg_database.thermo.libraries['PlasmaThermo'].entries['[Ar2p]'].data
    wilhoit = data.to_wilhoit()
    nasa = wilhoit.to_nasa(100.0, 5000.0, 1000.0)
    assert data.E0 is None
    for temperature in data.Tdata.value_si:
        assert nasa.get_heat_capacity(temperature) == pytest.approx(
            data.get_heat_capacity(temperature), rel=5e-3)
        assert nasa.get_enthalpy(temperature) == pytest.approx(
            data.get_enthalpy(temperature), abs=100.0)
        assert nasa.get_entropy(temperature) == pytest.approx(
            data.get_entropy(temperature), abs=0.1)
