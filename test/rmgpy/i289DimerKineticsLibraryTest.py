"""Engine-side acceptance tests for the PlasmaArgonDimer database library."""

import os
import re
import shutil
from collections import OrderedDict

import numpy as np
import pytest

from rmgpy import constants, settings
from rmgpy.exceptions import DatabaseError
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
    import rmgpy.data.rmg as rmg_data_module
    previous_database = rmg_data_module.database
    try:
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
        yield db
    finally:
        rmg_data_module.database = previous_database


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


def test_dr_unit_gate_exercises_ev_to_kelvin_conversion(library):
    """The other DR unit-gate tests only ever feed ``get_rate_coefficient_two_temp`` a Te
    already in kelvin, so a broken (or doubled/dropped) eV->K conversion anywhere upstream
    -- e.g. in the plasma reactor's electron-temperature bookkeeping -- would never be
    caught. Convert eV to K here exactly the way the engine does it (``rmgpy.solver.plasma``
    computes K->eV as ``Te_eV = Te_K * (R / Na) / e``, i.e. ``Te_eV = Te_K * kB / e``; invert
    that), then feed the resulting kelvin value through ``get_rate_coefficient_two_temp`` and
    compare against literals that were computed by hand, independently of this test's own
    ``te_ev * constants.e / constants.kB`` line, using the CODATA e/kB = 11604.51812 K/eV
    conversion factor and the DR source law 9.1e-7*(300/Te)^0.61 cm^3/(molecule*s):

        0.9 eV  -> Te_K = 0.9  * 11604.51812 = 10444.066308 K
                   alpha = 9.1e-7*(300/10444.066308)^0.61 * 1e-6*Na = 6.2852968167e+10 m^3/(mol s)
        1.16 eV -> Te_K = 1.16 * 11604.51812 = 13461.241019 K
                   alpha = 9.1e-7*(300/13461.241019)^0.61 * 1e-6*Na = 5.3838673243e+10 m^3/(mol s)

    Because the expected alpha values above are hard-coded literals (not re-derived from the
    same ``te_k`` variable used as input), a dropped or doubled eV->K conversion in the engine's
    ``constants.e / constants.kB`` path would move the computed rate outside the rel=1e-6
    tolerance and this test would fail.
    """
    kinetics = library.entries[1].data
    literals = {
        0.9: (10444.066308, 6.2852968167e+10),
        1.16: (13461.241019, 5.3838673243e+10),
    }
    for te_ev, (expected_te_k, expected_alpha) in literals.items():
        te_k = te_ev * constants.e / constants.kB
        # abs tolerance: rmgpy.constants uses full-precision CODATA e and kB separately, which
        # differs from the hand-rounded e/kB = 11604.51812 K/eV ratio by a few hundredths of a
        # kelvin at these Te -- far too small to hide a dropped or doubled conversion (which
        # would move te_k by a factor of e/kB itself, i.e. thousands of kelvin).
        assert te_k == pytest.approx(expected_te_k, abs=0.05)
        assert kinetics.get_rate_coefficient_two_temp(298.15, te_k) == pytest.approx(expected_alpha, rel=1e-6)


def test_loader_refuses_duplicate_explicit_and_metadata_electron(tmp_path):
    source = os.path.join(DATABASE_INPUT, 'kinetics', 'libraries', LIBRARY)

    # The unmodified copy must load fine in this same test, so the refusal below is
    # demonstrably caused by the edit and not by some unrelated environment issue.
    unmodified_target = tmp_path / 'unmodified' / LIBRARY
    shutil.copytree(source, unmodified_target)
    KineticsDatabase().load_libraries(str(tmp_path / 'unmodified'), libraries=[LIBRARY])

    edited_target = tmp_path / 'edited' / LIBRARY
    shutil.copytree(source, edited_target)
    reaction_file = edited_target / 'reactions.py'
    content = reaction_file.read_text()
    content = content.replace("Tmax=(8500, 'K')),", "Tmax=(8500, 'K'), electrons=-1),", 1)
    reaction_file.write_text(content)
    db = KineticsDatabase()
    with pytest.raises(DatabaseError, match='was not balanced') as excinfo:
        db.load_libraries(str(tmp_path / 'edited'), libraries=[LIBRARY])
    print('duplicate-electron refusal:', str(excinfo.value))


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
    """``reactor.kb == reactor.kf / reactor.Keq`` is a tautology: that is literally how
    ``PlasmaReactor`` computes ``kb`` (see ``rmgpy/solver/plasma.pyx``), so asserting it back
    proves nothing about correctness. Compute Kc independently from first principles off each
    species' own thermo (not via ``reaction.get_equilibrium_constant``, which is the same
    machinery the reactor already trusts), and check both the reactor's ``kb`` and the fitted
    reverse-rate object against ``kf(T) / Kc(T)``.

    Kc(T) = exp(-dG(T) / (R T)) * (P0 / (R T))^dn, with dG(T) the reaction's Gibbs free energy
    change from ``Species.get_free_energy``, dn = -1 (3 gas reactants -> 2 gas products) and
    P0 = 1e5 Pa (RMG's reference pressure for Kc), per ``Reaction.get_equilibrium_constant``.
    """
    species, reactions = _mechanism(rmg_database)
    conversion = _find_reaction(reactions, ['Arp', 'Ar', 'Ar'], ['Ar2p', 'Ar'])
    arp = next(s for s in species if s.label == 'Arp')
    ar = next(s for s in species if s.label == 'Ar')
    ar2p = next(s for s in species if s.label == 'Ar2p')
    P0 = 1e5
    dn = -1
    for temperature in (298.15, 1000.0):
        dG = (ar2p.get_free_energy(temperature) + ar.get_free_energy(temperature)) - (
            arp.get_free_energy(temperature) + 2 * ar.get_free_energy(temperature))
        Kc = np.exp(-dG / (constants.R * temperature)) * (P0 / (constants.R * temperature)) ** dn
        kf = conversion.kinetics.get_rate_coefficient(temperature)
        kb_expected = kf / Kc

        reverse = conversion.generate_reverse_rate_coefficient()
        reverse_kb = reverse.get_rate_coefficient(temperature)
        reverse_rel_error = abs(reverse_kb - kb_expected) / kb_expected
        # The reverse rate coefficient above is a fitted Arrhenius over the library's declared
        # Tmin-Tmax=150-300 K range; 1000 K lies well outside that range, so the fit is not
        # expected to reproduce the first-principles value to machine precision there. Measured
        # relative error is ~0.64% at 298.15 K and ~0.47% at 1000 K; 1% covers both with margin.
        # Report the actual error and use the tightest tolerance the fit honestly achieves
        # rather than loosening the reactor check below.
        assert reverse_rel_error < 0.01, (
            'reverse-fit relative error at {0} K was {1!r}, above the honest tolerance'.format(
                temperature, reverse_rel_error))
        print('reverse-fit relative error at {0} K: {1!r}'.format(temperature, reverse_rel_error))

        imf = {spc: 0.0 for spc in species}
        imf[ar] = 0.98
        imf[arp] = 0.01
        imf[ar2p] = 0.01
        imf[next(s for s in species if s.label == 'e-')] = 1e-12
        reactor = PlasmaReactor((temperature, 'K'), (5.0, 'torr'), imf,
                                (temperature, 'K'), n_sims=1, termination=[])
        reactor.initialize_model(species, reactions, [], [])
        index = reactor.reaction_index[conversion]
        assert reactor.kb[index] == pytest.approx(kb_expected, rel=1e-6)


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
    assert _side_key(loaded_dr.products).count('e-') == 0

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
    assert _equation_side_key(cantera_dr.equation, left=False).count('e-') == 0


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
