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

"""Admission, wrapped kinetics, and conversion regressions for resolved states."""

import pickle
from copy import deepcopy
from math import exp
from pathlib import Path
from types import SimpleNamespace

import pytest

from rmgpy.exceptions import ExcitedSpeciesThermoError, VibrationalManifoldError

from rmgpy.chemkin import load_chemkin_file, save_chemkin_file, save_species_dictionary
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import KineticsLibrary, LibraryReaction
from rmgpy.data.thermo import ThermoDatabase, ThermoLibrary
from rmgpy.exceptions import NetworkError, NonEquilibriumReverseRateError
from rmgpy.kinetics import (Arrhenius, Chebyshev, Lindemann, MultiArrhenius,
                            MultiPDepArrhenius, PDepArrhenius, ThirdBody, Troe,
                            TwoTemperaturePlasma)
from rmgpy.reaction import Reaction
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.rmg.pdep import PDepNetwork, PDepReaction
from rmgpy.solver.simple import SimpleReactor
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial

FIXTURE = Path(__file__).resolve().parents[2] / 'test_data' / 'superelastic' / 'thermal'


def exact_state_fixture_equilibrium_constant(temperature):
    """Return Kc for N2v1 -> N2 directly from the fixture NASA coefficients."""
    # N2: (a1, a6, a7) = (4.5, -1000, 3); N2v1: (3.5, 2400, 3).
    return temperature * exp(3400 / temperature - 1)


class ElectronArrhenius(Arrhenius):
    uses_electron_temperature = True


class ElectronChebyshev(Chebyshev):
    uses_electron_temperature = True


def reaction(resolved=True, reversible=True, kinetics=None, cls=Reaction):
    excited = Species(label='N2v1').from_adjacency_list(
        ('vibrationallevel 1\n' if resolved else '') +
        '1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}')
    ground = Species(label='N2').from_smiles('N#N')
    for index, sp in enumerate((excited, ground)):
        sp.thermo = NASA(polynomials=[NASAPolynomial(
            coeffs=[3.5, 0, 0, 0, 0, 1000 if index == 0 else 0, 1],
            Tmin=(200, 'K'), Tmax=(3000, 'K'))], Tmin=(200, 'K'), Tmax=(3000, 'K'),
            E0=(8314.462618 if index == 0 else 0, 'J/mol'))
    return cls(reactants=[excited], products=[ground], reversible=reversible,
               kinetics=kinetics if kinetics is not None else TwoTemperaturePlasma(A=(3e5, 's^-1'), n=1))


def install_exact_state_thermo(monkeypatch):
    import rmgpy.data.rmg as data_rmg
    database = ThermoDatabase()
    library = ThermoLibrary(label='ExactStateFixture')
    library.load(str(Path(__file__).resolve().parents[2] / 'test_data' / 'excited_states' / 'thermo.py'),
                 database.local_context, {})
    database.libraries = {library.label: library}
    database.library_order = [library.label]
    monkeypatch.setattr(data_rmg, 'database', SimpleNamespace(thermo=database, solvation=None))
    return database


def load_exact_state_thermo(reaction):
    for species in reaction.reactants + reaction.products:
        species.thermo = None
        species.get_thermo_data()


def wrapped(name, leaf):
    if name == 'multi':
        return MultiArrhenius(arrhenius=[leaf])
    if name == 'plog':
        return PDepArrhenius(pressures=([1, 100], 'bar'), arrhenius=[leaf, leaf])
    if name == 'multi-plog':
        return MultiPDepArrhenius(arrhenius=[wrapped('plog', leaf)])
    if name == 'nested':
        return wrapped('multi-plog', wrapped('multi', leaf))
    if name == 'third-body':
        return ThirdBody(arrheniusLow=leaf)
    if name == 'lindemann':
        return Lindemann(arrheniusLow=leaf, arrheniusHigh=Arrhenius(A=(3e5, 's^-1')))
    if name == 'troe':
        return Troe(arrheniusLow=leaf, arrheniusHigh=Arrhenius(A=(3e5, 's^-1')), alpha=0.5, T3=(100, 'K'), T1=(1000, 'K'))
    raise AssertionError(name)


WRAPPERS = ['multi', 'plog', 'multi-plog', 'nested', 'third-body', 'lindemann', 'troe']


@pytest.mark.parametrize('name', WRAPPERS)
def test_wrapped_rate_guard(name):
    leaf = ElectronArrhenius(A=(3e5, 's^-1')) if name in WRAPPERS[4:] else TwoTemperaturePlasma(A=(3e5, 's^-1'))
    rxn = reaction(kinetics=wrapped(name, leaf))
    with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
        rxn.check_reverse_from_equilibrium_supported()


@pytest.mark.parametrize('name', WRAPPERS)
def test_wrapped_rate_library_load(tmp_path, name):
    (tmp_path / 'dictionary.txt').write_text((FIXTURE / 'dictionary.txt').read_text())
    leaf = ElectronArrhenius(A=(3e5, 'm^3/(mol*s)')) if name in WRAPPERS[4:] else TwoTemperaturePlasma(A=(3e5, 'm^3/(mol*s)'))
    rate = wrapped(name, leaf)
    (tmp_path / 'reactions.py').write_text("entry(index=1, label='N2v1 + N2 <=> N2 + N2', kinetics=rate)\n")
    db = KineticsDatabase()
    context = dict(db.local_context, rate=rate)
    with pytest.raises(NonEquilibriumReverseRateError, match='Resolved species: N2v1'):
        KineticsLibrary(label=name).load(str(tmp_path / 'reactions.py'), context, db.global_context)


def test_duplicate_plog_merge_then_resolve(tmp_path):
    (tmp_path / 'dictionary.txt').write_text((FIXTURE / 'dictionary.txt').read_text().replace('vibrationallevel 1\n', ''))
    # Distinct molecular states are assigned after the initially unresolved load.
    (tmp_path / 'reactions.py').write_text("\n".join(
        "entry(index=%d, label='N2v1 + N2 <=> N2 + N2', duplicate=True, kinetics=PDepArrhenius(pressures=([1,100], 'bar'), arrhenius=[TwoTemperaturePlasma(A=(3e5,'m^3/(mol*s)'))]*2))" % i
        for i in (1, 2)))
    db = KineticsDatabase()
    lib = KineticsLibrary(label='merged')
    lib.load(str(tmp_path / 'reactions.py'), db.local_context, db.global_context)
    rxn = lib.get_library_reactions()[0]
    assert isinstance(rxn.kinetics, MultiPDepArrhenius)
    rxn.reactants[0].molecule[0].vibrational_level = 1
    with pytest.raises(NonEquilibriumReverseRateError, match='N2v1'):
        CoreEdgeReactionModel().add_reaction_to_core(rxn)


@pytest.mark.parametrize('thermal', [False, True], ids=['te', 'thermal'])
@pytest.mark.parametrize('reverse_first', [False, True])
def test_cross_library_explicit_reverse(monkeypatch, tmp_path, thermal, reverse_first):
    import rmgpy.data.rmg as data_rmg
    db = KineticsDatabase()
    libraries = []
    for name, equation, amplitude in [
            ('forward', 'N2v1 + N2 => N2 + N2', 1000),
            ('reverse', 'N2 + N2 => N2v1 + N2', 2000)]:
        directory = tmp_path / name
        directory.mkdir()
        (directory / 'dictionary.txt').write_text((FIXTURE / 'dictionary.txt').read_text())
        rate = ("Arrhenius(A=({0},'m^3/(mol*s)'))" if thermal else
                "TwoTemperaturePlasma(A=({0},'m^3/(mol*s)'))").format(amplitude)
        (directory / 'reactions.py').write_text('entry(index=1, label=%r, reversible=False, kinetics=%s)\n' % (equation, rate))
        lib = KineticsLibrary(label=name)
        lib.load(str(directory / 'reactions.py'), db.local_context, db.global_context)
        db.libraries[name] = lib
        libraries.append(lib)
    thermo = ThermoDatabase()
    thermo_library = ThermoLibrary(label='ExactStateFixture')
    thermo_library.load(str(Path(__file__).resolve().parents[2] / 'test_data' / 'excited_states' / 'thermo.py'),
                        thermo.local_context, {})
    thermo.libraries = {thermo_library.label: thermo_library}
    thermo.library_order = [thermo_library.label]
    monkeypatch.setattr(data_rmg, 'database', SimpleNamespace(
        kinetics=db, thermo=thermo, solvation=None))
    model = CoreEdgeReactionModel()
    model.declare_vibrational_manifold(Species(label='N2').from_smiles('N#N'))
    for lib in libraries[::(-1 if reverse_first else 1)]:
        rxn, is_new = model.make_new_reaction(
            lib.get_library_reactions()[0], generate_thermo=False, generate_kinetics=False)
        assert is_new
        model.add_reaction_to_core(rxn)
    assert len(model.core.reactions) == 2
    assert all(not rxn.reversible for rxn in model.core.reactions)
    for rxn in model.core.reactions:
        excited_reactant = any(molecule.has_resolved_state()
                               for species in rxn.reactants for molecule in species.molecule)
        excited_product = any(molecule.has_resolved_state()
                              for species in rxn.products for molecule in species.molecule)
        assert excited_reactant != excited_product
        expected = 1000 if excited_reactant else 2000
        assert rxn.kinetics.get_rate_coefficient(300, 12000) == pytest.approx(expected)


@pytest.mark.parametrize('thermal', [False, True], ids=['te', 'thermal'])
@pytest.mark.parametrize('reverse_first', [False, True])
def test_cross_library_explicit_reverse_refuses_undeclared_manifold(
        monkeypatch, tmp_path, thermal, reverse_first):
    import rmgpy.data.rmg as data_rmg
    db = KineticsDatabase()
    libraries = []
    for name, equation in [('forward', 'N2v1 + N2 => N2 + N2'), ('reverse', 'N2 + N2 => N2v1 + N2')]:
        directory = tmp_path / name
        directory.mkdir()
        (directory / 'dictionary.txt').write_text((FIXTURE / 'dictionary.txt').read_text())
        rate = "Arrhenius(A=(1000,'m^3/(mol*s)'))" if thermal else "TwoTemperaturePlasma(A=(1000,'m^3/(mol*s)'))"
        (directory / 'reactions.py').write_text(
            'entry(index=1, label=%r, reversible=False, kinetics=%s)\n' % (equation, rate))
        lib = KineticsLibrary(label=name)
        lib.load(str(directory / 'reactions.py'), db.local_context, db.global_context)
        db.libraries[name] = lib
        libraries.append(lib)
    monkeypatch.setattr(data_rmg, 'database', SimpleNamespace(kinetics=db))
    model = CoreEdgeReactionModel()
    for lib in libraries[::(-1 if reverse_first else 1)]:
        with pytest.raises(VibrationalManifoldError, match='N2v1'):
            model.make_new_reaction(
                lib.get_library_reactions()[0], generate_thermo=False, generate_kinetics=False)
    assert model.core.reactions == []


def imported_reaction(tmp_path, route):
    if route == 'constructed-template':
        return reaction(cls=TemplateReaction)
    if route == 'pickle':
        return pickle.loads(pickle.dumps(reaction()))
    if route == 'late-state':
        (tmp_path / 'dictionary.txt').write_text((FIXTURE / 'dictionary.txt').read_text().replace('vibrationallevel 1\n', ''))
        (tmp_path / 'reactions.py').write_text("entry(index=1, label='N2v1 + N2 <=> N2 + N2', kinetics=TwoTemperaturePlasma(A=(3e5,'m^3/(mol*s)')))\n")
        db = KineticsDatabase()
        lib = KineticsLibrary(label='late')
        lib.load(str(tmp_path / 'reactions.py'), db.local_context, db.global_context)
        rxn = lib.get_library_reactions()[0]
        rxn.reactants[0].molecule[0].vibrational_level = 1
        # Network admission needs a unimolecular channel; use the loaded state and rate.
        return Reaction(reactants=[rxn.reactants[0]], products=[rxn.products[0]], kinetics=rxn.kinetics)
    rxn = reaction(cls=LibraryReaction)
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    electron.thermo = deepcopy(rxn.products[0].thermo)
    rxn.kinetics = TwoTemperaturePlasma(A=(3e5, 's^-1'))
    if route == 'chemkin':
        save_chemkin_file(str(tmp_path / 'chem.inp'), [electron] + rxn.reactants + rxn.products, [rxn])
    else:
        from rmgpy.chemkin import get_species_identifier
        electron_label = get_species_identifier(electron)
        reactant_label = get_species_identifier(rxn.reactants[0])
        product_label = get_species_identifier(rxn.products[0])
        (tmp_path / 'chem.inp').write_text(
            'ELEMENTS N E END\nSPECIES {0} {1} {2} END\nREACTIONS\n'
            '{1}={2} 3e5 0 0\nTDEP/{0}/\nEND\n'.format(
                electron_label, reactant_label, product_label))
    save_species_dictionary(str(tmp_path / 'dictionary.txt'), [electron] + rxn.reactants + rxn.products)
    _, reactions = load_chemkin_file(str(tmp_path / 'chem.inp'), str(tmp_path / 'dictionary.txt'))
    assert len(reactions) == 1
    assert isinstance(reactions[0].kinetics, TwoTemperaturePlasma)
    assert reactions[0].reactants[0].molecule[0].has_resolved_state()
    return reactions[0]


@pytest.mark.parametrize('route', ['chemkin', 'tdep', 'pickle', 'constructed-template', 'late-state'])
@pytest.mark.parametrize('admission', ['core', 'edge', 'network', 'path'])
def test_admission_routes(tmp_path, route, admission):
    model = CoreEdgeReactionModel()
    expected = ExcitedSpeciesThermoError if route == 'chemkin' else NonEquilibriumReverseRateError
    with pytest.raises(expected, match='vibrationallevel 1'):
        # Chemkin now refuses at write/read; constructed and late-state routes
        # still exercise the individual admission guards.
        rxn = imported_reaction(tmp_path, route)
        if admission == 'network':
            model.add_reaction_to_unimolecular_networks(rxn, rxn.reactants[-1])
        elif admission == 'path':
            PDepNetwork().add_path_reaction(rxn)
        else:
            getattr(model, 'add_reaction_to_' + admission)(rxn)
    assert not model.core.reactions and not model.edge.reactions and not model.network_list


@pytest.mark.parametrize('direction', ['single', 'keep-forward', 'keep-reverse'])
def test_pdep_update_refuses_reversible_conversion(direction):
    model = CoreEdgeReactionModel()
    first = reaction(reversible=False, cls=PDepReaction)
    first.network = PDepNetwork(index=1)
    model.add_reaction_to_core(first)
    if direction != 'single':
        second = PDepReaction(reactants=first.products, products=first.reactants, reversible=False, kinetics=first.kinetics)
        second.network = PDepNetwork(index=2)
        model.add_reaction_to_core(second)
        first.products[0].thermo.polynomials[0].coeffs = [3.5, 0, 0, 0, 0, 2000 if direction == 'keep-reverse' else -2000, 1]
    original = list(model.core.reactions)
    with pytest.raises(NonEquilibriumReverseRateError, match='N2v1'):
        model.update_unimolecular_reaction_networks()
    assert model.core.reactions == original
    assert all(not rxn.reversible for rxn in original)


@pytest.mark.parametrize('reversible', [False, True])
def test_network_reverse_fitting_refuses_original_te_rate(reversible):
    rxn = reaction(reversible=reversible)
    rxn.network_kinetics = Arrhenius(A=(3e5, 's^-1'))
    with pytest.raises(NonEquilibriumReverseRateError, match='N2v1'):
        rxn.generate_reverse_rate_coefficient(network_kinetics=True)


@pytest.mark.parametrize('reversible', [False, True])
def test_chebyshev_conversion_refuses_before_fit(reversible):
    rate = ElectronChebyshev(coeffs=[[5, 0], [0, 0]], kunits='s^-1', Tmin=(300,'K'), Tmax=(2000,'K'), Pmin=(1,'bar'), Pmax=(100,'bar'))
    rxn = reaction(reversible=reversible, kinetics=rate, cls=LibraryReaction)
    rxn.elementary_high_p = True
    with pytest.raises(NonEquilibriumReverseRateError, match='N2v1'):
        rxn.generate_high_p_limit_kinetics()
    assert rxn.network_kinetics is None


def test_specific_collider_preserved_in_refusal(tmp_path):
    (tmp_path / 'dictionary.txt').write_text((FIXTURE / 'dictionary.txt').read_text())
    (tmp_path / 'reactions.py').write_text("entry(index=1, label='N2v1 (+N2) <=> N2 (+N2)', kinetics=TwoTemperaturePlasma(A=(3e5,'s^-1')))\n")
    db = KineticsDatabase()
    with pytest.raises(NonEquilibriumReverseRateError) as caught:
        KineticsLibrary(label='collider').load(str(tmp_path / 'reactions.py'), db.local_context, db.global_context)
    assert 'N2v1 (+N2) <=> N2 (+N2)' in str(caught.value)


@pytest.mark.parametrize('name', WRAPPERS)
def test_thermal_wrappers_remain_admissible(name):
    rxn = reaction(kinetics=wrapped(name, Arrhenius(A=(3e5, 's^-1'))))
    rxn.check_reverse_from_equilibrium_supported()
    model = CoreEdgeReactionModel()
    model.add_reaction_to_core(rxn)
    model.add_reaction_to_edge(rxn)
    PDepNetwork().add_path_reaction(rxn)
    if name in ('multi', 'plog', 'multi-plog', 'nested'):
        with pytest.raises(ExcitedSpeciesThermoError, match='N2v1'):
            rxn.generate_reverse_rate_coefficient()


@pytest.mark.parametrize('name', WRAPPERS[:4])
def test_thermal_wrappers_generate_quantitative_reverse_with_exact_state_thermo(monkeypatch, name):
    install_exact_state_thermo(monkeypatch)
    rxn = reaction(kinetics=wrapped(name, Arrhenius(A=(3e5, 's^-1'))))
    load_exact_state_thermo(rxn)
    reverse = rxn.generate_reverse_rate_coefficient()
    temperature, pressure = 1000, 1e5
    assert reverse.get_rate_coefficient(temperature, pressure) == pytest.approx(
        rxn.get_rate_coefficient(temperature, pressure) / rxn.get_equilibrium_constant(temperature),
        rel=1e-4)


def test_surface_thermal_control():
    from rmgpy.kinetics.surface import SurfaceArrhenius
    rxn = reaction(kinetics=SurfaceArrhenius(A=(3e5, 's^-1')))
    site = Species(label='X').from_adjacency_list('1 X u0 p0 c0')
    rxn.reactants.append(site)
    rxn.products.append(site)
    rxn.check_reverse_from_equilibrium_supported()
    CoreEdgeReactionModel().add_reaction_to_core(rxn)
    CoreEdgeReactionModel().add_reaction_to_edge(rxn)


def test_plog_simple_reactor_refuses_without_electron():
    rxn = reaction(kinetics=wrapped('plog', TwoTemperaturePlasma(A=(1000, 's^-1'), n=1)))
    reactor = SimpleReactor(T=(300,'K'), P=(1,'bar'), initial_mole_fractions={rxn.reactants[0]: 0.1, rxn.products[0]: 0.9}, termination=[])
    with pytest.raises(NonEquilibriumReverseRateError, match='N2v1'):
        reactor.initialize_model(rxn.reactants + rxn.products, [rxn], [], [])



def test_estimated_kinetics_checked_after_assignment(monkeypatch):
    model = CoreEdgeReactionModel()
    rxn = reaction(cls=TemplateReaction)
    rate = rxn.kinetics
    rxn.kinetics = None
    monkeypatch.setattr(model, 'generate_kinetics', lambda _: (rate, None, None, True))
    with pytest.raises(NonEquilibriumReverseRateError, match='N2v1'):
        model.apply_kinetics_to_reaction(rxn)


def test_participant_normalization_checked(monkeypatch):
    model = CoreEdgeReactionModel()
    rxn = reaction(resolved=False, cls=TemplateReaction)
    replacement = reaction().reactants[0]
    monkeypatch.setattr(model, 'make_new_species', lambda sp, **kwargs: (replacement if sp is rxn.reactants[0] else sp, True))
    with pytest.raises(NonEquilibriumReverseRateError, match='N2v1'):
        model.make_new_reaction(rxn, check_existing=False, generate_thermo=False, generate_kinetics=False)
    assert not model.new_reaction_list


@pytest.mark.parametrize('change', ['reversibility', 'kinetics', 'state'])
def test_repeated_admission_rechecks_changes(change):
    rxn = reaction(resolved=change != 'state', reversible=change != 'reversibility',
                   kinetics=Arrhenius(A=(3e5,'s^-1')) if change == 'kinetics' else None)
    model = CoreEdgeReactionModel()
    model.add_reaction_to_core(rxn)
    if change == 'reversibility':
        rxn.reversible = True
    elif change == 'kinetics':
        rxn.kinetics = TwoTemperaturePlasma(A=(3e5,'s^-1'))
    else:
        rxn.reactants[0].molecule[0].vibrational_level = 1
    with pytest.raises(NonEquilibriumReverseRateError, match='N2v1'):
        model.add_reaction_to_core(rxn)


def test_te_arrhenius_reverse_fitting_refused():
    rxn = reaction(kinetics=ElectronArrhenius(A=(3e5,'s^-1')))
    with pytest.raises(NonEquilibriumReverseRateError, match='N2v1'):
        rxn.generate_reverse_rate_coefficient()


def test_network_preserves_explicit_reverse():
    rxn = reaction(reversible=False)
    reverse = Reaction(reactants=rxn.products, products=rxn.reactants, reversible=False, kinetics=rxn.kinetics)
    network = PDepNetwork()
    network.add_path_reaction(rxn)
    network.add_path_reaction(reverse)
    assert network.path_reactions == [rxn, reverse]


@pytest.mark.parametrize('name', WRAPPERS)
def test_irreversible_wrappers_and_unresolved_controls(name):
    leaf = ElectronArrhenius(A=(3e5,'s^-1')) if name in WRAPPERS[4:] else TwoTemperaturePlasma(A=(3e5,'s^-1'))
    for resolved, reversible in [(True, False), (False, True)]:
        rxn = reaction(resolved=resolved, reversible=reversible, kinetics=wrapped(name, leaf))
        CoreEdgeReactionModel().add_reaction_to_core(rxn)
        CoreEdgeReactionModel().add_reaction_to_edge(rxn)
        network = PDepNetwork()
        if name in WRAPPERS[4:]:
            with pytest.raises(
                    NetworkError,
                    match=r'Electron reactions cannot enter pressure-dependent networks: N2v1.*N2'):
                network.add_path_reaction(rxn)
        else:
            network.add_path_reaction(rxn)
            assert network.path_reactions == [rxn]


def test_thermal_high_pressure_conversion_control(monkeypatch):
    install_exact_state_thermo(monkeypatch)
    rate = Chebyshev(coeffs=[[5,0],[0,0]], kunits='s^-1', Tmin=(300,'K'), Tmax=(2000,'K'), Pmin=(1,'bar'), Pmax=(100,'bar'))
    rxn = reaction(kinetics=rate, cls=LibraryReaction)
    load_exact_state_thermo(rxn)
    rxn.elementary_high_p = True
    assert rxn.generate_high_p_limit_kinetics()
    assert rxn.network_kinetics.get_rate_coefficient(1000) == pytest.approx(1e5)
    reverse = rxn.generate_reverse_rate_coefficient(network_kinetics=True)
    assert reverse.get_rate_coefficient(1000) == pytest.approx(
        1e5 / exact_state_fixture_equilibrium_constant(1000),
        rel=1e-4)


def test_thermal_high_pressure_conversion_refuses_unattested_thermo():
    rate = Chebyshev(coeffs=[[5,0],[0,0]], kunits='s^-1', Tmin=(300,'K'), Tmax=(2000,'K'), Pmin=(1,'bar'), Pmax=(100,'bar'))
    rxn = reaction(kinetics=rate, cls=LibraryReaction)
    rxn.elementary_high_p = True
    assert rxn.generate_high_p_limit_kinetics()
    with pytest.raises(ExcitedSpeciesThermoError, match='N2v1'):
        rxn.generate_reverse_rate_coefficient(network_kinetics=True)



def test_density_only_wrapped_leaf_refused():
    class DensityArrhenius(Arrhenius):
        uses_electron_density = True

    rxn = reaction(kinetics=wrapped('multi', DensityArrhenius(A=(3e5,'s^-1'))))
    with pytest.raises(NonEquilibriumReverseRateError, match='N2v1'):
        rxn.check_reverse_from_equilibrium_supported()


def test_barrier_conversion_refuses_before_erasing_te_declaration():
    from rmgpy.kinetics import ArrheniusEP

    class ElectronArrheniusEP(ArrheniusEP):
        uses_electron_temperature = True

    rate = ElectronArrheniusEP(A=(3e5,'s^-1'), n=0, alpha=0.5, E0=(10,'kJ/mol'))
    rxn = reaction(reversible=False, kinetics=rate)
    with pytest.raises(NonEquilibriumReverseRateError, match='N2v1'):
        rxn.fix_barrier_height()
    assert rxn.kinetics is rate


def test_model_tabulated_conversion_refuses_before_fit():
    from rmgpy.kinetics import KineticsData

    class ElectronKineticsData(KineticsData):
        uses_electron_temperature = True

        def to_arrhenius(self):
            return Arrhenius(A=(1e5, 's^-1'))

    rate = ElectronKineticsData(Tdata=([300,500,1000,1500],'K'), kdata=([1e5]*4,'s^-1'))
    rxn = reaction(reversible=False, kinetics=rate, cls=LibraryReaction)
    rxn.library = 'tabulated'
    with pytest.raises(VibrationalManifoldError, match='N2v1'):
        CoreEdgeReactionModel().make_new_reaction(rxn, check_existing=False, generate_thermo=False)
    assert rxn.kinetics is rate



def test_direct_duplicate_merge_rechecks_wrapped_rate():
    from rmgpy.data.base import Entry

    rate = wrapped('plog', TwoTemperaturePlasma(A=(3e5,'s^-1')))
    rxn = reaction()
    rxn.duplicate = True
    other = Reaction(reactants=rxn.reactants, products=rxn.products, duplicate=True)
    lib = KineticsLibrary(label='direct-merge')
    lib.entries = {1: Entry(index=1, item=rxn, data=rate),
                   2: Entry(index=2, item=other, data=rate)}
    with pytest.raises(NonEquilibriumReverseRateError, match='N2v1'):
        lib.convert_duplicates_to_multi()
    assert isinstance(lib.entries[1].data, MultiPDepArrhenius)



def thermal_reaction_with_te_network_surrogate(reversible=True):
    rate = Chebyshev(coeffs=[[5, 0], [0, 0]], kunits='s^-1',
                     Tmin=(300, 'K'), Tmax=(2000, 'K'),
                     Pmin=(1, 'bar'), Pmax=(100, 'bar'))
    rxn = reaction(reversible=reversible, kinetics=rate)
    rxn.network_kinetics = ElectronArrhenius(A=(3e5, 's^-1'))
    return rxn


@pytest.mark.parametrize('reversible', [False, True])
def test_selected_network_surrogate_reverse_fitting_refused(reversible):
    rxn = thermal_reaction_with_te_network_surrogate(reversible)
    with pytest.raises(NonEquilibriumReverseRateError, match='N2v1'):
        rxn.generate_reverse_rate_coefficient(network_kinetics=True)


@pytest.mark.parametrize('admission', ['model', 'path'])
def test_selected_network_surrogate_admission_refused(admission):
    rxn = thermal_reaction_with_te_network_surrogate()
    with pytest.raises(
            NetworkError,
            match=r'Electron reactions cannot enter pressure-dependent networks: N2v1.*N2'):
        if admission == 'model':
            CoreEdgeReactionModel().add_reaction_to_unimolecular_networks(rxn, rxn.reactants[0])
        else:
            PDepNetwork().add_path_reaction(rxn)


def test_selected_network_surrogate_pdep_update_refused():
    rxn = thermal_reaction_with_te_network_surrogate(reversible=False)
    network = PDepNetwork(source=rxn.reactants[:])
    network.path_reactions = [rxn]
    network.valid = True
    job = SimpleNamespace(
        output_file=None, Tmin=SimpleNamespace(value_si=300),
        Tmax=SimpleNamespace(value_si=2000), Pmin=SimpleNamespace(value_si=1e5),
        Pmax=SimpleNamespace(value_si=1e7), Tlist=SimpleNamespace(value_si=[300, 1000]),
        Plist=SimpleNamespace(value_si=[1e5, 1e7]), maximum_grain_size=None,
        minimum_grain_count=250, method='modified strong collision',
        interpolation_model=('Chebyshev', 2, 2), active_j_rotor=True,
        active_k_rotor=False, rmgmode=True)
    with pytest.raises(
            NetworkError,
            match=r'Electron reactions cannot enter pressure-dependent networks: N2v1.*N2'):
        network.update(CoreEdgeReactionModel(), job)


def test_unused_network_surrogate_keeps_core_thermal_rate(monkeypatch):
    install_exact_state_thermo(monkeypatch)
    rxn = thermal_reaction_with_te_network_surrogate()
    load_exact_state_thermo(rxn)
    rxn.check_reverse_from_equilibrium_supported()
    CoreEdgeReactionModel().add_reaction_to_core(rxn)
    reverse = rxn.generate_reverse_rate_coefficient(network_kinetics=False)
    assert reverse.get_rate_coefficient(1000, 1e5) == pytest.approx(
        rxn.kinetics.get_rate_coefficient(1000, 1e5) / rxn.get_equilibrium_constant(1000),
        rel=0.03)


def test_unused_network_surrogate_refuses_unattested_thermo():
    rxn = thermal_reaction_with_te_network_surrogate()
    rxn.check_reverse_from_equilibrium_supported()
    CoreEdgeReactionModel().add_reaction_to_core(rxn)
    with pytest.raises(ExcitedSpeciesThermoError, match='N2v1'):
        rxn.generate_reverse_rate_coefficient(network_kinetics=False)
