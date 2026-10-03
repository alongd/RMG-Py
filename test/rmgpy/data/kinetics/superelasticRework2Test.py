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

"""Resolved reversal and direction census, with baseline unresolved controls."""

from copy import deepcopy
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from rmgpy.chemkin import chemkin_duplicate_flags, mark_duplicate_reaction, mark_duplicate_reactions
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.common import find_degenerate_reactions
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import KineticsFamily, TemplateReaction
from rmgpy.data.kinetics.library import KineticsLibrary
from rmgpy.exceptions import NetworkError, NonEquilibriumReverseRateError, SpeciesIdentityError
from rmgpy.kinetics import (Arrhenius, MultiArrhenius, MultiPDepArrhenius,
                            PDepArrhenius, TwoTemperaturePlasma)
from rmgpy.kinetics.diffusionLimited import diffusion_limiter
from rmgpy.reaction import Reaction
from rmgpy.rmg.model import CoreEdgeReactionModel, ReactionModel, are_identical_species_references
from rmgpy.rmg.pdep import PDepNetwork
from rmgpy.rmg.pdep import PDepReaction
from rmgpy.solver.simple import SimpleReactor
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial
from rmgpy.tools.diffmodels import compare_model_reactions
from rmgpy.tools.isotopes import cluster, compare_isotopomers
from rmgpy.tools.mergemodels import combine_models


class ElectronArrhenius(Arrhenius):
    uses_electron_temperature = True


def make_reaction(resolved=True, reversible=True, kinetics=None, cls=Reaction):
    excited = Species(label='N2v1', index=1).from_adjacency_list(
        ('vibrationallevel 1\n' if resolved else '') +
        '1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}')
    ground = Species(label='N2', index=2).from_smiles('N#N')
    for index, sp in enumerate((excited, ground)):
        sp.thermo = NASA(polynomials=[NASAPolynomial(
            coeffs=[3.5, 0, 0, 0, 0, 1000 if index == 0 else 0, 1],
            Tmin=(200, 'K'), Tmax=(3000, 'K'))], Tmin=(200, 'K'), Tmax=(3000, 'K'),
            E0=(8314.462618 if index == 0 else 0, 'J/mol'))
    return cls(reactants=[excited], products=[ground], reversible=reversible,
               kinetics=kinetics if kinetics is not None else TwoTemperaturePlasma(A=(1e5, 's^-1'), n=1))


def explicit_pair(resolved=True, cls=Reaction):
    forward = make_reaction(resolved=resolved, reversible=False, cls=cls)
    reverse = cls(reactants=forward.products[:], products=forward.reactants[:],
                  reversible=False, kinetics=TwoTemperaturePlasma(A=(2e5, 's^-1'), n=1))
    return forward, reverse


def assert_merged_direction_and_rate(reactions, expected):
    matches = [reaction for reaction in reactions
               if are_identical_species_references(reaction, expected)
               and reaction.kinetics.is_identical_to(expected.kinetics)]
    assert len(matches) == 1
    actual = matches[0]
    assert actual.reversible is expected.reversible
    assert actual.kinetics.get_rate_coefficient(300, 12000) == expected.kinetics.get_rate_coefficient(300, 12000)


@pytest.mark.parametrize('reverse_first', [False, True])
def test_witness_model_merge_retains_both_directions(reverse_first):
    forward, reverse = explicit_pair()
    first, second = (reverse, forward) if reverse_first else (forward, reverse)
    species = forward.reactants + forward.products
    merged = ReactionModel(species=species, reactions=[first]).merge(
        ReactionModel(species=species, reactions=[second]))
    assert len(merged.reactions) == 2
    assert_merged_direction_and_rate(merged.reactions, first)
    assert_merged_direction_and_rate(merged.reactions, second)


def test_witness_direct_arrhenius_reverse_refuses_original():
    rxn = make_reaction(reversible=False)
    rxn.network_kinetics = Arrhenius(A=(3e5, 's^-1'))
    with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
        rxn.reverse_arrhenius_rate(rxn.network_kinetics, 's^-1')


def test_witness_unresolved_plog_initializes_like_base():
    rate = PDepArrhenius(pressures=([1, 100], 'bar'),
                         arrhenius=[TwoTemperaturePlasma(A=(1000, 's^-1'), n=1)] * 2)
    rxn = make_reaction(resolved=False, kinetics=rate)
    reactor = SimpleReactor(T=(300, 'K'), P=(1, 'bar'),
                            initial_mole_fractions={rxn.reactants[0]: 0.1, rxn.products[0]: 0.9}, termination=[])
    reactor.initialize_model(rxn.reactants + rxn.products, [rxn], [], [])
    assert reactor.kf[0] == pytest.approx(300000)
    assert reactor.kb[0] == pytest.approx(300000 / rxn.get_equilibrium_constant(300))
    assert rxn.get_reverse_from_equilibrium_refusal() is None


@pytest.mark.parametrize('resolved', [False, True], ids=['unresolved', 'resolved'])
def test_witness_molecule_participants_are_normalized(resolved):
    rxn = make_reaction(resolved=resolved, cls=TemplateReaction)
    rxn.reactants = [sp.molecule[0] for sp in rxn.reactants]
    rxn.products = [sp.molecule[0] for sp in rxn.products]
    model = CoreEdgeReactionModel()
    if resolved:
        with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
            model.make_new_reaction(rxn, check_existing=False, generate_thermo=False, generate_kinetics=False)
    else:
        added, is_new = model.make_new_reaction(rxn, check_existing=False, generate_thermo=False, generate_kinetics=False)
        assert is_new
        assert all(isinstance(sp, Species) for sp in added.reactants + added.products)
        assert all(not mol.has_resolved_state() for sp in added.reactants + added.products for mol in sp.molecule)


@pytest.mark.parametrize('options', [{}, {'check_identical': True}, {'check_only_label': True}])
def test_direction_census_isomorphism(options):
    forward, reverse = explicit_pair()
    assert forward.is_isomorphic(reverse, **options)
    assert not forward.is_same_reaction(reverse, **options)
    assert forward.is_isomorphic(deepcopy(forward), **options)


def test_direction_census_template_product_shortcut():
    forward, reverse = explicit_pair(cls=TemplateReaction)
    forward.is_forward = True
    reverse.is_forward = False
    assert forward.is_isomorphic(reverse, check_template_rxn_products=True)
    assert not forward.is_same_reaction(reverse, check_template_rxn_products=True)
    assert not forward.is_same_reaction(reverse, check_template_rxn_products=True, check_identical=True)


@pytest.mark.parametrize('site', ['has-template', 'matches-species', 'reference-identity'])
def test_direction_census_queries(site):
    forward, reverse = explicit_pair()
    if site == 'has-template':
        assert not forward.has_template(reverse.reactants, reverse.products)
        assert forward.has_template(forward.reactants, forward.products)
    elif site == 'matches-species':
        assert not forward.matches_species(reverse.reactants, reverse.products)
        assert forward.matches_species(forward.reactants, forward.products)
    else:
        assert not are_identical_species_references(forward, reverse)
        assert are_identical_species_references(forward, forward)


def test_direction_census_merge_models_tool():
    forward, reverse = explicit_pair()
    species = forward.reactants + forward.products
    merged = combine_models([ReactionModel(species=species, reactions=[forward]),
                              ReactionModel(species=species, reactions=[reverse])])
    assert len(merged.reactions) == 2
    assert_merged_direction_and_rate(merged.reactions, forward)
    assert_merged_direction_and_rate(merged.reactions, reverse)


def test_direction_census_diff_models():
    forward, reverse = explicit_pair()
    common, first, second = compare_model_reactions(
        ReactionModel(reactions=[forward]), ReactionModel(reactions=[reverse]))
    assert common == []
    assert first == [forward] and second == [reverse]


@pytest.mark.parametrize('site', ['compare-isotopomers', 'cluster'])
def test_direction_census_isotope_identity(site):
    forward, reverse = explicit_pair()
    if site == 'compare-isotopomers':
        assert not compare_isotopomers(forward, reverse)
    else:
        assert len(cluster([forward, reverse])) == 2


def test_direction_census_degeneracy():
    forward, reverse = explicit_pair(cls=TemplateReaction)
    forward.is_forward, reverse.is_forward = True, False
    forward.template = reverse.template = ['test-template']
    family = KineticsFamily(label='test')
    family.own_reverse = True
    found = find_degenerate_reactions([forward, reverse], same_reactants=0, kinetics_family=family)
    assert found == [forward, reverse]


def test_direction_census_depository():
    forward, reverse = explicit_pair()
    entry = Entry(index=1, label='forward', item=forward, data=forward.kinetics)
    depository = SimpleNamespace(entries={1: entry}, label='test/training')
    found = KineticsFamily(label='test').get_kinetics_from_depository(depository, reverse, [], 1)
    assert found == []


def test_direction_census_library_duplicates():
    forward, reverse = explicit_pair()
    library = KineticsLibrary(label='explicit-pair')
    library.entries = {index: Entry(index=index, item=rxn, data=rxn.kinetics)
                       for index, rxn in enumerate((forward, reverse), 1)}
    library.check_for_duplicates(mark_duplicates=True)
    assert not forward.duplicate and not reverse.duplicate
    library.mark_valid_duplicates([forward], [reverse])
    assert not forward.duplicate and not reverse.duplicate
    # Marked, supported PLOG contributions must enter the actual conversion
    # loop; unmarked or unsupported TwoTemperaturePlasma entries skip it.
    for entry in library.entries.values():
        entry.item.duplicate = True
        entry.data = PDepArrhenius(pressures=([1, 100], 'bar'),
                                  arrhenius=[entry.data] * 2)
        entry.item.kinetics = entry.data
    library.convert_duplicates_to_multi()
    assert len(library.entries) == 2


@pytest.mark.parametrize('site', ['pairwise', 'group', 'flags'])
def test_direction_census_chemkin_duplicate_marking(site):
    forward, reverse = explicit_pair()
    if site == 'pairwise':
        mark_duplicate_reaction(reverse, [forward])
    elif site == 'group':
        mark_duplicate_reactions([forward, reverse])
    else:
        assert chemkin_duplicate_flags([forward, reverse]) == [False, False]
    assert not forward.duplicate and not reverse.duplicate


@pytest.mark.parametrize('site', ['add-path', 'merge-path', 'merge-net'])
def test_direction_census_network_identity(site):
    forward, reverse = explicit_pair()
    network = PDepNetwork(source=forward.reactants)
    if site == 'add-path':
        network.add_path_reaction(forward)
        network.add_path_reaction(reverse)
        assert network.path_reactions == [forward, reverse]
    else:
        other = PDepNetwork(source=forward.reactants)
        attribute = 'path_reactions' if site == 'merge-path' else 'net_reactions'
        setattr(network, attribute, [forward])
        setattr(other, attribute, [reverse])
        network.merge(other)
        assert getattr(network, attribute) == [forward, reverse]


@pytest.mark.parametrize('resolved', [False, True])
def test_thermal_reverse_fitting_control(resolved):
    rate = Arrhenius(A=(3e5, 's^-1'))
    rxn = make_reaction(resolved=resolved, reversible=False, kinetics=rate)
    reverse = rxn.reverse_arrhenius_rate(rate, 's^-1')
    assert reverse.get_rate_coefficient(1000) == pytest.approx(rate.get_rate_coefficient(1000) / rxn.get_equilibrium_constant(1000))


def test_reverse_census_diffusion_limiter():
    rxn = make_reaction(reversible=False, kinetics=ElectronArrhenius(A=(3e5, 's^-1')))
    with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
        diffusion_limiter.get_effective_rate(rxn, 1000)


def test_unresolved_direct_flag_keeps_baseline_refusal():
    rxn = make_reaction(resolved=False)
    with pytest.raises(NonEquilibriumReverseRateError):
        rxn.check_reverse_from_equilibrium_supported()


def test_unresolved_model_merge_keeps_baseline_identity():
    forward, reverse = explicit_pair(resolved=False)
    species = forward.reactants + forward.products
    merged = ReactionModel(species=species, reactions=[forward]).merge(
        ReactionModel(species=species, reactions=[reverse]))
    assert merged.reactions == [forward]


def test_direction_census_labels_cannot_erase_resolved_direction():
    forward, reverse = explicit_pair()
    for sp in forward.reactants + forward.products:
        sp.label = 'same-label'
    assert forward.is_isomorphic(reverse, check_only_label=True)
    assert not forward.is_same_reaction(reverse, check_only_label=True)


def test_direction_census_thermal_explicit_pdep_pair():
    forward = make_reaction(reversible=False, kinetics=Arrhenius(A=(3e5, 's^-1')), cls=PDepReaction)
    reverse = PDepReaction(reactants=forward.products[:], products=forward.reactants[:],
                           reversible=False, kinetics=Arrhenius(A=(1e5, 's^-1')))
    forward.network, reverse.network = PDepNetwork(index=1), PDepNetwork(index=2)
    model = CoreEdgeReactionModel()
    model.add_reaction_to_core(forward)
    model.add_reaction_to_core(reverse)
    model.update_unimolecular_reaction_networks()
    assert model.core.reactions == [forward, reverse]
    assert not forward.reversible and not reverse.reversible


@pytest.mark.parametrize('name', ['multi', 'plog', 'multi-plog'])
def test_unresolved_executable_wrappers_keep_baseline_reactor(name):
    if name == 'multi':
        rate = MultiArrhenius(arrhenius=[ElectronArrhenius(A=(1000, 's^-1'))])
    else:
        rate = PDepArrhenius(pressures=([1, 100], 'bar'),
                             arrhenius=[TwoTemperaturePlasma(A=(1000, 's^-1'), n=1)] * 2)
        if name == 'multi-plog':
            rate = MultiPDepArrhenius(arrhenius=[rate])
    rxn = make_reaction(resolved=False, kinetics=rate)
    reactor = SimpleReactor(T=(300, 'K'), P=(1, 'bar'),
                            initial_mole_fractions={rxn.reactants[0]: 0.1, rxn.products[0]: 0.9}, termination=[])
    reactor.initialize_model(rxn.reactants + rxn.products, [rxn], [], [])
    expected = rate.get_rate_coefficient(300, 1e5)
    assert reactor.kf[0] == pytest.approx(expected)
    assert reactor.kb[0] == pytest.approx(expected / rxn.get_equilibrium_constant(300))


def surface_reaction(rate):
    gas = Species(label='CH3v1').from_smiles('[CH3]')
    gas.molecule[0].vibrational_level = 1
    vacant = Species(label='X').from_adjacency_list('1 X u0 p0 c0')
    adsorbed = Species(label='CH3X').from_adjacency_list(
        '1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}\n'
        '2 H u0 p0 c0 {1,S}\n3 H u0 p0 c0 {1,S}\n'
        '4 H u0 p0 c0 {1,S}\n5 X u0 p0 c0 {1,S}')
    for index, sp in enumerate((gas, vacant, adsorbed)):
        sp.thermo = NASA(polynomials=[NASAPolynomial(
            coeffs=[0 if index == 1 else 3.5, 0, 0, 0, 0, 1000 if index == 0 else 0, 1],
            Tmin=(200, 'K'), Tmax=(3000, 'K'))], Tmin=(200, 'K'), Tmax=(3000, 'K'))
    return Reaction(reactants=[gas, vacant], products=[adsorbed], reversible=False, kinetics=rate)


HELPERS = ['arrhenius', 'surface-arrhenius', 'sticking', 'surface-charge', 'gas-charge']


def helper_input(name, electron_dependent):
    from rmgpy.kinetics import ArrheniusChargeTransfer, SurfaceArrhenius, SurfaceChargeTransfer, StickingCoefficient
    bases = {'arrhenius': Arrhenius, 'surface-arrhenius': SurfaceArrhenius,
             'sticking': StickingCoefficient, 'surface-charge': SurfaceChargeTransfer,
             'gas-charge': ArrheniusChargeTransfer}
    base = bases[name]
    cls = type('DeclaredElectron' + base.__name__, (base,), {'uses_electron_temperature': True}) if electron_dependent else base
    if name == 'sticking':
        rate = cls(A=0.1)
    elif name == 'surface-arrhenius':
        rate = cls(A=(3e5, 'm^3/(mol*s)'))
    elif name in ('surface-charge', 'gas-charge'):
        rate = cls(A=(3e5, 's^-1'), electrons=-1, V0=(0, 'V'))
    else:
        rate = cls(A=(3e5, 's^-1'))
    rxn = surface_reaction(rate) if name in ('sticking', 'surface-arrhenius') else make_reaction(kinetics=rate, reversible=False)
    if name in ('surface-charge', 'gas-charge'):
        rxn.electrons = -1
    functions = {'arrhenius': rxn.reverse_arrhenius_rate,
                 'surface-arrhenius': rxn.reverse_surface_arrhenius_rate,
                 'sticking': rxn.reverse_sticking_coeff_rate,
                 'surface-charge': rxn.reverse_surface_charge_transfer_rate,
                 'gas-charge': rxn.reverse_arrhenius_charge_transfer_rate}
    args = (rate, 's^-1', 2.5e-5) if name == 'sticking' else (rate, 's^-1')
    return rxn, functions[name], args


@pytest.mark.parametrize('name', HELPERS)
def test_reverse_census_direct_helpers(name):
    rxn, helper, args = helper_input(name, True)
    with pytest.raises(NonEquilibriumReverseRateError, match='Resolved species:'):
        helper(*args)


@pytest.mark.parametrize('name', HELPERS)
def test_reverse_census_direct_helpers_selected_rate(name):
    rxn, helper, args = helper_input(name, True)
    rxn.kinetics = Arrhenius(A=(3e5, 's^-1'))
    with pytest.raises(NonEquilibriumReverseRateError, match='Resolved species:'):
        helper(*args)


@pytest.mark.parametrize('name', HELPERS)
def test_reverse_census_direct_helpers_thermal_control(name):
    rxn, helper, args = helper_input(name, False)
    result = helper(*args)
    if hasattr(result, 'V0'):
        value = result.get_rate_coefficient(350, result.V0.value_si)
    else:
        value = result.get_rate_coefficient(350)
    assert np.isfinite(value) and value > 0


@pytest.mark.parametrize('selected_surrogate', [False, True])
@pytest.mark.parametrize('through_method', [False, True])
def test_reverse_census_microcanonical(selected_surrogate, through_method):
    from rmgpy.pdep.reaction import calculate_microcanonical_rate_coefficient
    from rmgpy.species import TransitionState
    from rmgpy.statmech import Conformer
    rxn = make_reaction(reversible=False)
    rxn.transition_state = TransitionState(conformer=Conformer(E0=(20000, 'J/mol'), modes=[]))
    rxn.network_kinetics = Arrhenius(A=(3e5, 's^-1'), Ea=(20000, 'J/mol'))
    if selected_surrogate:
        rxn.kinetics = Arrhenius(A=(3e5, 's^-1'), Ea=(20000, 'J/mol'))
        rxn.network_kinetics = ElectronArrhenius(A=(3e5, 's^-1'), Ea=(20000, 'J/mol'))
    energies = np.linspace(0, 60000, 101)
    densities = np.ones((101, 1))
    with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
        if through_method:
            rxn.calculate_microcanonical_rate_coefficient(
                energies, np.array([0], dtype=np.int_), densities, densities, 1000)
        else:
            calculate_microcanonical_rate_coefficient(
                rxn, energies, np.array([0], dtype=np.int_), densities, densities, 1000)


@pytest.mark.parametrize('thermal', [False, True])
def test_reverse_census_interpolation_fit(thermal):
    from rmgpy.pdep.reaction import fit_interpolation_model
    rate = Arrhenius(A=(3e5, 's^-1')) if thermal else TwoTemperaturePlasma(A=(3e5, 's^-1'))
    rxn = make_reaction(reversible=False, kinetics=rate)
    temperatures = np.array([300, 500, 1000, 1500.])
    pressures = np.array([1e4, 1e5])
    rates = np.full((4, 2), 3e5)
    args = (rxn, temperatures, pressures, rates, ('PDepArrhenius',), 300, 1500, 1e4, 1e5)
    if thermal:
        result = fit_interpolation_model(*args)
        assert result.get_rate_coefficient(1000, 1e5) == pytest.approx(3e5)
    else:
        with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
            fit_interpolation_model(*args)


@pytest.fixture(scope='module')
def thermal_network():
    # Reuse the established statmech input, so entry tests start from a real
    # initialized network rather than from fabricated matrices or missing fields.
    path = Path(__file__).resolve().parents[2] / 'pdep' / 'networkTest.py'
    spec = importlib.util.spec_from_file_location('reversal_network_input', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    fixture = module.TestNetwork()
    fixture.setup_class()
    network = fixture.network
    for configuration in network.isomers + network.reactants + network.products:
        for species in configuration.species:
            species.thermo = None
    network.initialize(300, 1000, 1e4, 1e5, maximum_grain_size=5000, minimum_grain_count=20,
                       active_j_rotor=True, active_k_rotor=True)
    network.set_conditions(1000, 1e5)
    network.calculate_equilibrium_ratios()
    return network


NETWORK_SITES = ['initialize', 'conditions', 'microcanonical', 'rates', 'matrix',
                 'msc', 'rs', 'cse', 'cse-georgievskii', 'sls',
                 'direct-matrix', 'direct-msc', 'direct-rs', 'direct-cse',
                 'direct-cse-advanced', 'direct-cse-georgievskii', 'direct-sls', 'direct-sls-rates']


def run_network_site(network, name):
    from rmgpy.pdep import cse, me, msc, rs, sls
    calls = {
        'initialize': lambda: network.initialize(300, 1000, 1e4, 1e5, maximum_grain_size=5000, minimum_grain_count=20,
                                                 active_j_rotor=True, active_k_rotor=True),
        'conditions': lambda: network.set_conditions(1000, 1e5),
        'microcanonical': network.calculate_microcanonical_rates,
        'rates': lambda: network.calculate_rate_coefficients([1000], [1e5], 'modified strong collision'),
        'matrix': network.generate_full_me_matrix,
        'msc': network.apply_modified_strong_collision_method,
        'rs': network.apply_reservoir_state_method,
        'cse': network.apply_chemically_significant_eigenvalues_method,
        'cse-georgievskii': lambda: network.apply_chemically_significant_eigenvalues_method(method='georgievskii'),
        'sls': network.apply_simulation_least_squares_method,
        'direct-matrix': lambda: me.generate_full_me_matrix(network),
        'direct-msc': lambda: msc.apply_modified_strong_collision_method(network),
        'direct-rs': lambda: rs.apply_reservoir_state_method(network),
        'direct-cse': lambda: cse.apply_chemically_significant_eigenvalues_method(network),
        'direct-cse-advanced': lambda: cse.get_rate_coefficients_CSE_Advanced(network, 1000, 1e5),
        'direct-cse-georgievskii': lambda: cse.apply_chemically_significant_eigenvalues_method_georgievskii(network),
        'direct-sls': lambda: sls.apply_simulation_least_squares_method(network),
        'direct-sls-rates': lambda: sls.get_rate_coefficients_SLS(network, 1000, 1e5),
    }
    return calls[name]()


@pytest.mark.parametrize('name', NETWORK_SITES)
def test_reverse_census_network_after_late_declaration(thermal_network, name):
    from rmgpy.molecule import Molecule
    network = deepcopy(thermal_network)
    reaction = network.path_reactions[0]
    reaction.reactants[0].molecule = [Molecule().from_smiles('CCCCO')]
    reaction.reactants[0].molecule[0].vibrational_level = 1
    reaction.kinetics = TwoTemperaturePlasma(A=(3e5, 's^-1'))
    with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
        run_network_site(network, name)


def test_reverse_census_arkane_execute_reaches_rate_guard(thermal_network, monkeypatch):
    from arkane.pdep import PressureDependenceJob
    from rmgpy.molecule import Molecule

    network = deepcopy(thermal_network)
    reaction = network.path_reactions[0]
    reaction.reactants[0].molecule = [Molecule().from_smiles('CCCCO')]
    reaction.reactants[0].molecule[0].vibrational_level = 1
    reaction.kinetics = TwoTemperaturePlasma(A=(3e5, 's^-1'))
    job = PressureDependenceJob(
        network,
        Tlist=([300, 1000], 'K'),
        Plist=([1e4, 1e5], 'Pa'),
        method='modified strong collision',
        interpolationModel=('PDepArrhenius',),
    )
    monkeypatch.setattr(job, 'initialize', lambda: None)

    with pytest.raises(SpeciesIdentityError, match='PressureDependenceJob.execute'):
        job.execute(output_file=None, plot=False, print_summary=False)


@pytest.mark.parametrize('name', ['msc', 'rs', 'matrix'])
def test_reverse_census_network_thermal_control(thermal_network, name):
    network = deepcopy(thermal_network)
    output = run_network_site(network, name)
    assert output is not None


def test_reverse_census_rrkm_ignores_unused_network_surrogate(thermal_network):
    from rmgpy.molecule import Molecule
    network = deepcopy(thermal_network)
    reaction = network.path_reactions[0]
    reaction.reactants[0].molecule = [Molecule().from_smiles('CCCCO')]
    reaction.reactants[0].molecule[0].vibrational_level = 1
    reaction.kinetics = Arrhenius(A=(3e5, 's^-1'))
    reaction.network_kinetics = ElectronArrhenius(A=(3e5, 's^-1'))
    assert reaction.can_tst()
    with pytest.raises(
            NetworkError,
            match=r'Electron reactions cannot enter pressure-dependent networks: n-C4H10O.*n-C4H8.*H2O'):
        network.generate_full_me_matrix()


@pytest.mark.parametrize('site', ['leak', 'steady-state', 'rate-filter'])
def test_reverse_census_partial_network(site):
    from rmgpy.pdep.configuration import Configuration
    rxn = make_reaction(reversible=False, kinetics=ElectronArrhenius(A=(3e5, 's^-1')))
    network = PDepNetwork(source=rxn.products)
    if site == 'leak':
        network.path_reactions = [rxn]
        call = lambda: network.get_leak_coefficient(1000, 1e5)
    else:
        network.net_reactions = [rxn]
        network.isomers = [Configuration(*rxn.reactants), Configuration(*rxn.products)]
        call = (lambda: network.solve_ss_network(1000, 1e5)) if site == 'steady-state' else (
            lambda: network.get_rate_filtered_products(1000, 1e5, 0.1))
    with pytest.raises(
            NetworkError,
            match=r'Electron reactions cannot enter pressure-dependent networks: N2v1.*N2'):
        call()


def test_reverse_census_training_entry_fit(monkeypatch):
    rxn = make_reaction(reversible=False)
    data = rxn.kinetics
    entry_reaction = Reaction(reactants=[rxn.reactants[0].molecule[0]],
                              products=[rxn.products[0].molecule[0]], reversible=False)
    entry = Entry(index=1, item=entry_reaction, data=data, label='explicit-channel')
    reverse = TemplateReaction(reactants=rxn.products[:], products=rxn.reactants[:], template=['test'])
    database = KineticsDatabase()
    database.families['test'] = KineticsFamily(label='test')
    monkeypatch.setattr(database, 'generate_reactions_from_families', lambda *args, **kwargs: [reverse])
    thermo = SimpleNamespace(get_thermo_data=lambda sp: rxn.reactants[0].thermo)
    with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
        database.get_forward_reaction_for_family_entry(entry, 'test', thermo)


def test_reverse_census_isotope_equilibrium_flux():
    from rmgpy.tools.isotopes import ensure_correct_degeneracies
    rxn = make_reaction(kinetics=ElectronArrhenius(A=(3e5, 's^-1')), reversible=False)
    with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
        ensure_correct_degeneracies([rxn])


def test_reverse_census_model_comparison_before_removal(monkeypatch):
    from rmgpy.tools import diffmodels
    forward, reverse = explicit_pair()
    model1, model2 = ReactionModel(reactions=[forward]), ReactionModel(reactions=[reverse])
    monkeypatch.setattr(diffmodels.plt, 'show', lambda: None)
    diffmodels.compare_model_kinetics(model1, model2)
    assert model2.reactions == [reverse]
    diffmodels.plt.close('all')
    reverse.reversible = True
    with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
        diffmodels.compare_model_kinetics(model1, model2)
    assert model2.reactions == [reverse]


def test_reverse_census_arkane_output(tmp_path):
    from arkane.kinetics import KineticsJob
    from rmgpy.species import TransitionState
    from rmgpy.statmech import Conformer
    rxn = make_reaction(reversible=False, kinetics=Arrhenius(A=(3e5, 's^-1')))
    for sp in rxn.reactants + rxn.products:
        sp.conformer = Conformer(E0=sp.thermo.E0, modes=[])
    rxn.transition_state = TransitionState(conformer=Conformer(E0=(20000, 'J/mol'), modes=[]))
    job = KineticsJob(rxn)
    job.generate_kinetics()
    rxn.kinetics = ElectronArrhenius(A=(3e5, 's^-1'))
    job.usedTST = True
    with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
        job.write_output(str(tmp_path))


def test_direction_census_uncertainty_index():
    from rmgpy.tools.uncertainty import get_i_thing
    forward, reverse = explicit_pair()
    assert get_i_thing(reverse, [forward, reverse]) == 1
    assert get_i_thing(deepcopy(reverse), [forward, reverse]) == 1
    assert get_i_thing(deepcopy(forward.reactants[0]), forward.reactants) == 0


def test_direction_census_enlargement_summary():
    forward, reverse = explicit_pair()
    additions = [forward, reverse]
    CoreEdgeReactionModel().log_enlarge_summary([], additions, [], [], reactions_moved_from_edge=[forward])
    assert additions == [reverse]


def test_direction_census_family_reverse_attribute(monkeypatch):
    forward, reverse = explicit_pair(cls=TemplateReaction)
    forward.template = reverse.template = ['test']
    forward.is_forward = reverse.is_forward = True
    family = KineticsFamily(label='test')
    family.own_reverse = True
    monkeypatch.setattr(family, '_generate_reactions', lambda *args, **kwargs: [reverse])
    assert family.add_reverse_attribute(forward)
    assert forward.reverse is reverse
    assert forward.kinetics is not reverse.kinetics
    reverse.is_forward = False
    # Exercise both multiple-candidate comparisons, including the forbidden
    # retry, without letting degeneracy aggregation hide opposite candidates.
    import rmgpy.data.kinetics.family as family_module
    from rmgpy.exceptions import KineticsError
    monkeypatch.setattr(family_module, 'find_degenerate_reactions',
                        lambda reactions, *args, **kwargs: reactions)
    for retry in (False, True):
        candidates = iter([[], [forward, reverse]] if retry else [[forward, reverse]])
        monkeypatch.setattr(family, '_generate_reactions',
                            lambda *args, **kwargs: next(candidates))
        with pytest.raises(KineticsError):
            family.add_reverse_attribute(forward)


@pytest.mark.parametrize('resolved', [False, True])
def test_direction_census_arkane_network_export(tmp_path, resolved):
    from arkane.pdep import PressureDependenceJob
    from rmgpy.pdep.network import Network
    forward, reverse = explicit_pair(resolved=resolved)
    # Test exporting precomputed thermal network rates, without adding electron pdep support.
    forward.kinetics = Arrhenius(A=(3e5, 's^-1'), n=0, Ea=(0, 'J/mol'), T0=(1, 'K'))
    reverse.kinetics = Arrhenius(A=(1e5, 's^-1'), n=0, Ea=(0, 'J/mol'), T0=(1, 'K'))
    network = Network(label='explicit-pair', net_reactions=[forward, reverse])
    network.n_isom = 2
    network.n_reac = network.n_prod = 0
    job = PressureDependenceJob(network, Tlist=([300, 1000], 'K'), Plist=([1e4, 1e5], 'Pa'))
    job.K = np.full((2, 2, 2, 2), 1e5)
    path = tmp_path / 'output.py'
    if resolved:
        with pytest.raises(SpeciesIdentityError, match='N2v1'):
            job.save(str(path))
        assert not path.exists()
        return
    job.save(str(path))
    active = [line for line in path.read_text().splitlines() if line.startswith('pdepreaction(')]
    assert len(active) == 1


@pytest.mark.parametrize('declared_subclass', [False, True], ids=['native', 'declared-subclass'])
def test_arkane_pdep_output_refuses_resolved_electron_rates_before_writing(tmp_path, declared_subclass):
    from arkane.pdep import PressureDependenceJob
    from rmgpy.pdep.network import Network

    path_forward, path_reverse = explicit_pair()
    if declared_subclass:
        path_forward.kinetics = ElectronArrhenius(A=(1e5, 's^-1'), n=1, Ea=(0, 'J/mol'), T0=(1, 'K'))
        path_reverse.kinetics = ElectronArrhenius(A=(2e5, 's^-1'), n=1, Ea=(0, 'J/mol'), T0=(1, 'K'))
    forward = Reaction(reactants=path_forward.reactants[:], products=path_forward.products[:], reversible=False,
                       kinetics=Arrhenius(A=(1e5, 's^-1'), n=1, Ea=(0, 'J/mol'), T0=(1, 'K')))
    reverse = Reaction(reactants=path_reverse.reactants[:], products=path_reverse.products[:], reversible=False,
                       kinetics=Arrhenius(A=(2e5, 's^-1'), n=1, Ea=(0, 'J/mol'), T0=(1, 'K')))
    network = Network(label='explicit-pair', path_reactions=[path_forward, path_reverse],
                      net_reactions=[forward, reverse])
    network.n_isom = 2
    network.n_reac = network.n_prod = 0
    job = PressureDependenceJob(network, Tlist=([300, 1000], 'K'), Plist=([1e4, 1e5], 'Pa'))
    job.K = np.full((2, 2, 2, 2), 1e5)
    output = tmp_path / 'output.py'

    with pytest.raises(
            NonEquilibriumReverseRateError,
            match='reaction N2v1.*N2.*Resolved species: N2v1.*vibrationallevel 1'):
        job.save(str(output))

    assert not output.exists()
    assert not (tmp_path / 'chem.inp').exists()


def test_direction_census_library_lookup():
    forward, reverse = explicit_pair()
    library = KineticsLibrary(label='pair')
    library.entries = {i: Entry(index=i, item=r, data=r.kinetics)
                       for i, r in enumerate((forward, reverse), 1)}
    found = KineticsDatabase().generate_reactions_from_library(library, reverse.reactants[:], reverse.products[:])
    assert len(found) == 1
    assert found[0].reactants == reverse.reactants
    assert found[0].kinetics.get_rate_coefficient(300, 12000) == reverse.kinetics.get_rate_coefficient(300, 12000)


@pytest.mark.parametrize('site', ['core-overlap', 'edge-overlap', 'net-reuse'])
def test_direction_census_network_update(site, monkeypatch):
    from arkane.pdep import PressureDependenceJob
    from rmgpy.data.kinetics.library import LibraryReaction
    from rmgpy.statmech import Conformer, HarmonicOscillator
    forward = make_reaction(reversible=False, kinetics=Arrhenius(A=(3e5, 's^-1'), Ea=(10000, 'J/mol')))
    for species in forward.reactants + forward.products:
        species.conformer = Conformer(E0=species.thermo.E0, modes=[HarmonicOscillator(frequencies=([1000], 'cm^-1'))])
    bath = Species(label='Ar', reactive=False).from_smiles('[Ar]')
    model = CoreEdgeReactionModel()
    model.core.species = forward.reactants[:] + [bath]
    if site != 'edge-overlap':
        model.core.species += forward.products
    reverse_class = PDepReaction if site == 'net-reuse' else LibraryReaction
    reverse = reverse_class(reactants=forward.products[:], products=forward.reactants[:], reversible=False,
                            kinetics=Arrhenius(A=(1e5, 's^-1'), Ea=(0, 'J/mol')))
    network = PDepNetwork(index=1, source=forward.reactants[:])
    network.path_reactions = [forward]
    if site == 'net-reuse':
        network.net_reactions = [reverse]
    elif site == 'core-overlap':
        model.core.reactions = [reverse]
    else:
        model.edge.reactions = [reverse]
    job = PressureDependenceJob(network, Tlist=([300, 500, 1000, 1500], 'K'), Plist=([1e4, 1e5], 'Pa'),
                                method='modified strong collision', interpolationModel=('PDepArrhenius',))
    job.output_file = None
    # Isolate the real update's identity branches from its separately tested numerical solver.
    monkeypatch.setattr(network, 'initialize', lambda *args: None)
    monkeypatch.setattr(network, 'calculate_rate_coefficients', lambda *args: np.full((4, 2, 2, 2), 1e5))
    network.update(model, job)
    if site == 'net-reuse':
        assert len(network.net_reactions) == 2
    else:
        target = model.core if site == 'core-overlap' else model.edge
        assert len(target.reactions) == 2
        assert reverse in target.reactions
    added = network.net_reactions[-1]
    assert added.reactants == forward.reactants and added.products == forward.products
    assert not added.reversible and not reverse.reversible


@pytest.mark.parametrize('site', ['simple', 'liquid', 'mb-sampled', 'surface'])
def test_reverse_census_reactors(site):
    from rmgpy.solver.liquid import LiquidReactor
    from rmgpy.solver.mbSampled import MBSampledReactor
    from rmgpy.solver.surface import SurfaceReactor
    reaction = make_reaction(kinetics=ElectronArrhenius(A=(3e5, 's^-1')))
    fractions = {reaction.reactants[0]: 0.1, reaction.products[0]: 0.9}
    species = reaction.reactants + reaction.products
    if site == 'simple':
        reactor = SimpleReactor(T=(300, 'K'), P=(1, 'bar'), initial_mole_fractions=fractions, termination=[])
    elif site == 'liquid':
        reactor = LiquidReactor(T=(300, 'K'), initial_concentrations=fractions, termination=[])
    elif site == 'mb-sampled':
        reactor = MBSampledReactor(T=(300, 'K'), P=(1, 'bar'), initial_mole_fractions=fractions,
                                  k_sampling=(1, 's^-1'), constantSpeciesList=[], termination=[])
    else:
        vacancy = Species(label='X').from_adjacency_list('1 X u0 p0 c0')
        species.append(vacancy)
        reactor = SurfaceReactor(T=(300, 'K'), P_initial=(1, 'bar'), initial_gas_mole_fractions=fractions,
                                 initial_surface_coverages={vacancy: 1.0}, surface_volume_ratio=(100, 'm^-1'),
                                 surface_site_density=(2.5e-5, 'mol/m^2'), n_sims=1, termination=[])
    with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
        reactor.initialize_model(species, [reaction], [], [])


@pytest.mark.parametrize('reverse_first', [False, True])
@pytest.mark.parametrize('duplicate', [False, True])
def test_direction_census_same_family_admission(monkeypatch, reverse_first, duplicate):
    from rmgpy.rmg import model as model_module
    forward, reverse = explicit_pair(cls=TemplateReaction)
    family = KineticsFamily(label='test')
    family.own_reverse = True
    monkeypatch.setattr(model_module, 'get_family_library_object', lambda label: family)
    model = CoreEdgeReactionModel()
    for reaction in ((reverse, forward) if reverse_first else (forward, reverse)):
        reaction.family = 'test'
        reaction.duplicate = duplicate
        reaction.template = ['test']
        added, new = model.make_new_reaction(reaction, generate_thermo=False, generate_kinetics=False)
        assert new
        model.add_reaction_to_core(added)
    assert len(model.core.reactions) == 2
    assert all(not r.reversible for r in model.core.reactions)
