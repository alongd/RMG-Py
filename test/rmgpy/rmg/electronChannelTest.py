#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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

"""Independent irreversible electron channels retain both sourced rate laws."""

import importlib.util
from pathlib import Path

import pytest

from rmgpy import settings
from rmgpy.chemkin import mark_duplicate_reactions
from rmgpy.kinetics import Arrhenius, TwoTemperaturePlasma
from rmgpy.rmg.model import CoreEdgeReactionModel, ReactionModel
from rmgpy.thermo import ThermoData
# Load the sibling helper by path, as plasmaKineticsConsolidationTest loads its helper.
_helpers_spec = importlib.util.spec_from_file_location(
    'electronTestHelpers', Path(__file__).with_name('electronTestHelpers.py'),
)
_helpers = importlib.util.module_from_spec(_helpers_spec)
_helpers_spec.loader.exec_module(_helpers)
(database, ORDERS, LIBRARIES, controlled_thermo, admission_state,
 registry_contents, _COLLIDER_BOUNDARIES) = (
    _helpers.database, _helpers.ORDERS, _helpers.LIBRARIES, _helpers.controlled_thermo,
    _helpers.admission_state, _helpers.registry_contents, _helpers._COLLIDER_BOUNDARIES,
)


def library_pair(database, order):
    return [database.kinetics.libraries[name].get_library_reactions()[0] for name in order]


def admit_to_edge(database, order):
    model = CoreEdgeReactionModel()
    for reaction in library_pair(database, order):
        model.process_new_reactions(
            [reaction], None, generate_thermo=False, generate_kinetics=False,
        )
    return model


def assert_channels(reactions, database):
    assert len(reactions) == 2
    by_library = {reaction.library: reaction for reaction in reactions}
    assert set(by_library) == set(LIBRARIES)
    for library in LIBRARIES:
        reaction = by_library[library]
        entry = database.kinetics.libraries[library].entries[1]
        assert reaction.entry is entry
        assert reaction.kinetics is entry.data
        assert not reaction.reversible
        assert reaction.is_isomorphic(entry.item, either_direction=False)
    detachment = by_library['ElectronDetachment']
    attachment = by_library['ElectronAttachment']
    assert isinstance(detachment.kinetics, Arrhenius)
    assert isinstance(attachment.kinetics, TwoTemperaturePlasma)
    assert attachment.kinetics.Ea_e.value == pytest.approx(6.26, rel=0, abs=1e-12)
    assert detachment.kinetics.get_rate_coefficient(300) == pytest.approx(1.3850923748e8)
    assert attachment.kinetics.get_rate_coefficient_two_temp(300, 17406.78) > 0
    assert not detachment.duplicate and not attachment.duplicate


@pytest.mark.parametrize('order', ORDERS)
def test_independent_channels_survive_enlargement(database, order):
    model = admit_to_edge(database, order)
    for species in list(model.edge.species):
        model.enlarge(species)
    assert not model.edge.species and not model.edge.reactions
    assert_channels(model.core.reactions, database)


@pytest.mark.parametrize('order', ORDERS)
def test_independent_channels_survive_duplicate_marking(database, order):
    model = admit_to_edge(database, order)
    mark_duplicate_reactions(model.edge.reactions)
    assert_channels(model.edge.reactions, database)


@pytest.mark.parametrize('order', ORDERS)
def test_independent_channels_survive_edge_to_core_moves(database, order):
    model = admit_to_edge(database, order)
    moved = []
    for species in list(model.edge.species):
        moved.extend(model.add_species_to_core(species))
    assert not model.edge.reactions
    assert moved == model.core.reactions
    model.mark_chemkin_duplicates()
    assert_channels(model.core.reactions, database)


def separate_models(reactions):
    models = []
    for reaction in reactions:
        species = list(dict.fromkeys(reaction.reactants + reaction.products))
        for item in species:
            item.thermo = ThermoData(
                Tdata=([300, 1000], 'K'), Cpdata=([30, 30], 'J/(mol*K)'),
                H298=(0, 'kJ/mol'), S298=(100, 'J/(mol*K)'),
                Cp0=(30, 'J/(mol*K)'), CpInf=(30, 'J/(mol*K)'),
            )
        # Exercise merge's species mapping and participant order independently
        # of the shared, sorted references used by CoreEdgeReactionModel.
        reaction.reactants.reverse()
        reaction.products.reverse()
        models.append(ReactionModel(species=species, reactions=[reaction]))
    return models


@pytest.mark.parametrize('order', ORDERS)
def test_independent_channels_survive_model_merge(database, order):
    first, second = separate_models(library_pair(database, order))
    merged = first.merge(second)
    assert_channels(merged.reactions, database)
    assert len(merged.species) == 4
    assert all(item in merged.species for reaction in merged.reactions
               for item in reaction.reactants + reaction.products)
    assert sorted(merged.reactions[0].reactants) == sorted(merged.reactions[1].products)
    assert sorted(merged.reactions[0].products) == sorted(merged.reactions[1].reactants)


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('control', ('same_library', 'reversible', 'no_electron', 'same_direction'))
@pytest.mark.parametrize('operation', ('admission', 'merge'))
def test_existing_library_priority_controls(database, order, control, operation):
    reactions = library_pair(database, order)
    if control == 'same_library':
        reactions[1].library = reactions[1].family = reactions[0].library
    elif control == 'reversible':
        reactions[0].reversible = True
    elif control == 'no_electron':
        # Ordinary oxygen recombination/dissociation, with repeated participants.
        for reaction in reactions:
            oxygen = next(item for item in reaction.reactants + reaction.products if item.label == 'O')
            dioxygen = next(item for item in reaction.reactants + reaction.products if item.label == 'O2')
            if reaction.library == 'ElectronDetachment':
                reaction.reactants, reaction.products = [oxygen, oxygen], [dioxygen]
            else:
                reaction.reactants, reaction.products = [dioxygen], [oxygen, oxygen]
    else:
        reactions[1].reactants = reactions[0].reactants[:]
        reactions[1].products = reactions[0].products[:]

    if operation == 'merge':
        first, second = separate_models(reactions)
        admitted = first.merge(second).reactions
        expected = 1
    else:
        model = CoreEdgeReactionModel()
        for reaction in reactions:
            model.process_new_reactions([reaction], None, generate_thermo=False, generate_kinetics=False)
        admitted = model.edge.reactions
        # Same-library reverse pairs already coexist on the admission path.
        expected = 2 if control == 'same_library' else 1
    assert len(admitted) == expected
    assert admitted[0] is reactions[0]
    assert admitted[0].kinetics is reactions[0].kinetics


def model_snapshot(model):
    """Keep identities as well as values at the merge boundary."""
    return (tuple(model.species), tuple(model.reactions), [
        (reaction.reactants, tuple(reaction.reactants),
         reaction.products, tuple(reaction.products),
         reaction.pairs, None if reaction.pairs is None else tuple(reaction.pairs),
         reaction.specific_collider, reaction.kinetics, reaction.reversible,
         reaction.index, reaction.is_forward, reaction.rank, reaction.comment,
         reaction.label)
        for reaction in model.reactions
    ])


def assert_model_unchanged(model, snapshot):
    species, reactions, states = snapshot
    assert tuple(model.species) == species
    assert tuple(model.reactions) == reactions
    for reaction, state in zip(model.reactions, states):
        assert reaction.reactants is state[0]
        assert tuple(reaction.reactants) == state[1]
        assert reaction.products is state[2]
        assert tuple(reaction.products) == state[3]
        assert reaction.pairs is state[4]
        assert (None if reaction.pairs is None else tuple(reaction.pairs)) == state[5]
        assert (reaction.specific_collider, reaction.kinetics, reaction.reversible,
                reaction.index, reaction.is_forward, reaction.rank, reaction.comment,
                reaction.label) == state[6:]


@pytest.mark.parametrize('order', ORDERS)
def test_repeated_merges_preserve_channels_and_inputs(database, order):
    first, second = separate_models(library_pair(database, order))
    for model in (first, second):
        reaction = model.reactions[0]
        reaction.generate_pairs()
    snapshots = [model_snapshot(model) for model in (first, second)]
    outputs = [first.merge(second), second.merge(first)]
    for output in outputs:
        assert_channels(output.reactions, database)
        assert all(species in output.species for reaction in output.reactions
                   for species in reaction.reactants + reaction.products)
        for reaction in output.reactions:
            assert all(r in reaction.reactants and p in reaction.products
                       for r, p in reaction.pairs)
    for model, snapshot in zip((first, second), snapshots):
        assert_model_unchanged(model, snapshot)


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('kind', ('library', 'template', 'base'))
def test_merge_copies_retained_state_pairs_and_collider(database, order, kind):
    from rmgpy.data.kinetics.family import TemplateReaction
    from rmgpy.reaction import Reaction

    first, second = separate_models(library_pair(database, order))
    original = second.reactions[0]
    if kind == 'template':
        reaction = TemplateReaction(reactants=original.reactants[:], products=original.products[:],
                                    kinetics=original.kinetics, family='test')
    elif kind == 'base':
        reaction = Reaction(reactants=original.reactants[:], products=original.products[:],
                            kinetics=original.kinetics)
    else:
        reaction = original
    electron = next(s for s in second.species if s.is_electron())
    oxygen = next(s for s in second.species if s.label == 'O')
    dioxygen = next(s for s in second.species if s.label == 'O2')
    reaction.reactants = [oxygen, electron, electron]
    reaction.products = [dioxygen, electron, electron]
    reaction.pairs = [(oxygen, dioxygen), (electron, electron), (electron, electron)]
    reaction.specific_collider = oxygen
    reaction.is_forward = False
    reaction.rank, reaction.comment, reaction.label = 7, 'retained state', 'channel'
    reaction.allow_pdep_route = True
    reaction.allow_max_rate_violation = True
    second.reactions = [reaction]
    first.reactions = []
    snapshots = [model_snapshot(model) for model in (first, second)]
    merged = first.merge(second)
    for model, snapshot in zip((first, second), snapshots):
        assert_model_unchanged(model, snapshot)
    retained = merged.reactions[0]
    assert retained is not reaction
    assert type(retained) is type(reaction)
    assert retained.kinetics is reaction.kinetics
    assert retained.specific_collider in merged.species
    assert retained.specific_collider is next(s for s in first.species if s.label == 'O')
    assert retained.pairs == list(zip(retained.reactants, retained.products))
    assert retained.pairs is not reaction.pairs
    assert (retained.is_forward, retained.rank, retained.comment, retained.label) == (False, 7, 'retained state', 'channel')
    assert retained.allow_pdep_route and retained.allow_max_rate_violation
    if kind == 'library':
        assert retained.entry is reaction.entry
        assert retained.library == reaction.library


def unimolecular_pair(database, order):
    reactions = library_pair(database, order)
    for reaction in reactions:
        species = {s.label: s for s in reaction.reactants + reaction.products}
        if reaction.library == 'ElectronDetachment':
            reaction.reactants, reaction.products = [species['O-']], [species['O'], species['e-']]
            units = 's^-1'
        else:
            reaction.reactants, reaction.products = [species['O'], species['e-']], [species['O-']]
            units = 'm^3/(mol*s)'
        reaction.kinetics = Arrhenius(A=(1.0, units), n=0, Ea=(0, 'J/mol'))
        reaction.elementary_high_p = True
    return reactions


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('separate_sources', (False, True))
def test_pdep_admission_refuses_independent_channels(database, order, separate_sources):
    from types import SimpleNamespace
    from rmgpy.exceptions import NetworkError

    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    reactions = unimolecular_pair(database, order)
    sources = [next(s for s in r.reactants + r.products if s.label == 'O-') for r in reactions]
    if separate_sources:
        sources[1] = next(s for s in reactions[1].reactants + reactions[1].products if s.label == 'O')
    before = admission_state(model)
    for reaction, source in zip(reactions, sources):
        with pytest.raises(NetworkError, match='Electron reactions cannot enter pressure-dependent networks'):
            model.process_new_reactions([reaction], source, generate_thermo=False, generate_kinetics=False)
        assert admission_state(model) == before


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('operation', ('add_path', 'merge'))
def test_pdep_network_boundaries_refuse_channels(database, order, operation):
    from rmgpy.exceptions import NetworkError
    from rmgpy.pdep import Configuration
    from rmgpy.rmg.pdep import PDepNetwork

    reactions = unimolecular_pair(database, order)
    source = [next(s for s in reactions[0].reactants + reactions[0].products if s.label == 'O-')]
    first, second = PDepNetwork(source=source), PDepNetwork(source=source)
    second.path_reactions = [reactions[1]]  # Legacy network; new path insertion refuses immediately.
    second.products.append(Configuration(*reactions[1].products))
    before = (first.path_reactions[:], first.products[:], first.valid)
    with pytest.raises(NetworkError, match='Electron reactions cannot enter pressure-dependent networks'):
        if operation == 'add_path':
            first.add_path_reaction(reactions[0])
        else:
            first.merge(second)
    assert (first.path_reactions, first.products, first.valid) == before
    assert second.path_reactions == [reactions[1]]


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('participants', ('spectator', 'repeated', 'electron_only'))
def test_symmetric_same_direction_merge_deduplicates(database, order, participants):
    reactions = library_pair(database, order)
    for reaction in reactions:
        electron = next(s for s in reaction.reactants + reaction.products if s.is_electron())
        oxygen = next(s for s in reaction.reactants + reaction.products if s.label == 'O')
        side = {'spectator': [oxygen, electron], 'repeated': [oxygen, oxygen, electron],
                'electron_only': [electron]}[participants]
        reaction.reactants, reaction.products = side[:], side[::-1]
    first, second = separate_models(reactions)
    merged = first.merge(second)
    assert len(merged.reactions) == 1
    assert merged.reactions[0] is first.reactions[0]
    assert merged.reactions[0].kinetics is first.reactions[0].kinetics


@pytest.mark.parametrize('order', ORDERS)
def test_remerged_channels_initialize_plasma_and_load_cantera(database, order, tmp_path):
    import cantera as ct
    import numpy as np
    from rmgpy.solver.plasma import PlasmaReactor
    from rmgpy.species import Species
    from rmgpy.thermo import NASA, NASAPolynomial
    from rmgpy.transport import TransportData
    from collections import Counter
    from rmgpy.yaml_cantera2 import get_label, save_cantera_model

    first, second = separate_models(library_pair(database, order))
    first.merge(second)
    model = second.merge(first)
    assert_channels(model.reactions, database)
    cation = Species(label='Ar+').from_adjacency_list('1 Ar u1 p3 c+1')
    model.species.append(cation)
    for index, species in enumerate(model.species, 1):
        species.index = index
        species.transport_data = TransportData(
            shapeIndex=0 if len(species.molecule[0].atoms) == 1 else 1,
            epsilon=(100, 'K'), sigma=(3.0, 'angstrom'),
            dipoleMoment=(0, 'C*m'), polarizability=(0, 'angstrom^3'), rotrelaxcollnum=1)
        # Controlled constant-Cp NASA fixture; the electron's thermo is zero.
        coefficients = [0 if species.is_electron() else 3.5, 0, 0, 0, 0, 0, 0]
        species.thermo = NASA(polynomials=[NASAPolynomial(coeffs=coefficients, Tmin=(200, 'K'), Tmax=(3000, 'K'))],
                              Tmin=(200, 'K'), Tmax=(3000, 'K'))
    electron = next(s for s in model.species if s.is_electron())
    dioxygen = next(s for s in model.species if s.label == 'O2')
    reactor = PlasmaReactor(T=(300, 'K'), P=(5 * 133.322368, 'Pa'), Te=(17406.78, 'K'),
                            initial_mole_fractions={dioxygen: 1 - 2e-12, cation: 1e-12, electron: 1e-12},
                            termination=[], thermo_source_assertions={'O-': 'ion', 'Ar+': 'ion'})
    reactor.initialize_model(core_species=model.species, core_reactions=model.reactions,
                             edge_species=[], edge_reactions=[])
    # Controlled NASA fixtures are caller-asserted, not physical thermo evidence.
    assert reactor.thermo_provenance_diagnostics == {
        'O-': 'caller-asserted, not verified', 'Ar+': 'caller-asserted, not verified'}
    assert reactor.electron_species is electron
    assert len(reactor.kf) == 2
    assert np.all(np.isfinite(reactor.kf)) and np.all(reactor.kf > 0)
    assert np.all(reactor.kb == 0)
    for reaction, rate in zip(model.reactions, reactor.kf):
        kinetics = reaction.kinetics
        expected = (kinetics.get_rate_coefficient_two_temp(300, 17406.78)
                    if isinstance(kinetics, TwoTemperaturePlasma) else kinetics.get_rate_coefficient(300))
        assert rate == pytest.approx(expected)
    output = tmp_path / 'electron-channels.yaml'
    save_cantera_model(model, str(output))
    gas = ct.Solution(str(output))
    assert gas.n_reactions == 2
    assert all(not r.reversible and not r.duplicate for r in gas.reactions())
    gas.TP = 300, 5 * 133.322368
    for reaction in gas.reactions():
        source = next(r for r in model.reactions
                      if dict(Counter(get_label(s, model.species) for s in r.reactants)) == reaction.reactants)
        expected = (source.kinetics.get_rate_coefficient_two_temp(300, 17406.78)
                    if isinstance(source.kinetics, TwoTemperaturePlasma) else source.kinetics.get_rate_coefficient(300))
        rate = (reaction.rate(300, 17406.78) if isinstance(source.kinetics, TwoTemperaturePlasma)
                else reaction.rate(300))
        # Cantera uses kmol; its gas constant also differs slightly from RMG's.
        molar_rate = rate / 1000 ** (sum(reaction.reactants.values()) - 1)
        assert molar_rate == pytest.approx(expected, rel=1e-5)


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('operation', ('seed', 'library', 'seed_round_trip'))
def test_seed_library_restart_controls_keep_channels_and_sources(database, order, operation, monkeypatch, tmp_path):
    from rmgpy.data.base import Entry
    from rmgpy.data.kinetics.library import KineticsLibrary
    import rmgpy.rmg.model as model_module

    def fixture_thermo(species, solvent_name=None):
        species.thermo = ThermoData(Tdata=([300, 1000], 'K'), Cpdata=([30, 30], 'J/(mol*K)'),
                                   H298=(0, 'kJ/mol'), S298=(100, 'J/(mol*K)'),
                                   Cp0=(30, 'J/(mol*K)'), CpInf=(30, 'J/(mol*K)'))

    monkeypatch.setattr(model_module, 'submit', fixture_thermo)
    sources = library_pair(database, order)
    snapshots = [model_snapshot(ReactionModel(species=list(dict.fromkeys(r.reactants + r.products)), reactions=[r]))
                 for r in sources]
    entries = [database.kinetics.libraries[name].entries[1] for name in order]
    entry_sides = [(entry.item.reactants[:], entry.item.products[:], repr(entry.data)) for entry in entries]
    model = CoreEdgeReactionModel()
    if operation == 'seed_round_trip':
        # Restart loads a saved seed library, then calls add_seed_mechanism_to_core.
        # Exercise that boundary with production save/load and admission methods.
        seed = KineticsLibrary(name='restart', auto_generated=True)
        seed.entries = {i: Entry(index=i, label=r.to_labeled_str(), item=r, data=r.kinetics,
                                 long_desc='Originally from reaction library: ' + r.library + '\n')
                        for i, r in enumerate(sources, 1)}
        path = tmp_path / 'restart'
        path.mkdir()
        seed.save(str(path / 'reactions.py'))
        database.kinetics.load_libraries(str(tmp_path), libraries=['restart'], additive=True)
        model.add_seed_mechanism_to_core('restart')
    else:
        for name in order:
            if operation == 'seed':
                model.add_seed_mechanism_to_core(name)
            else:
                model.add_reaction_library_to_edge(name)
    admitted = model.core.reactions if operation != 'library' else model.edge.reactions
    if operation != 'seed_round_trip':
        assert_channels(admitted, database)
    else:
        assert len(admitted) == 2
        assert {r.library for r in admitted} == set(LIBRARIES)
        assert all(not r.reversible and not r.duplicate for r in admitted)
        for reaction in admitted:
            source = next(r for r in sources if r.library == reaction.library)
            assert reaction.is_isomorphic(source, either_direction=False)
            assert reaction.kinetics.is_identical_to(source.kinetics)
    for reaction, snapshot in zip(sources, snapshots):
        assert_model_unchanged(ReactionModel(species=list(snapshot[0]), reactions=[reaction]), snapshot)
    for entry, sides in zip(entries, entry_sides):
        assert (entry.item.reactants, entry.item.products, repr(entry.data)) == sides


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('control', ('same_library', 'reversible', 'no_electron', 'same_direction'))
def test_pdep_priority_controls(database, order, control):
    from rmgpy.rmg.pdep import PDepNetwork

    reactions = unimolecular_pair(database, order)
    if control == 'same_library':
        reactions[1].library = reactions[1].family = reactions[0].library
    elif control == 'reversible':
        reactions[0].reversible = True
    elif control == 'no_electron':
        for reaction in reactions:
            oxygen = next(s for s in reaction.reactants + reaction.products if s.label == 'O')
            dioxygen = next(s for s in database.kinetics.libraries[reaction.library].entries[1].item.reactants
                            + database.kinetics.libraries[reaction.library].entries[1].item.products if s.label == 'O2')
            reaction.reactants, reaction.products = ([oxygen, oxygen], [dioxygen]) if reaction.library == 'ElectronDetachment' else ([dioxygen], [oxygen, oxygen])
    else:
        reactions[1].reactants, reactions[1].products = reactions[0].reactants[:], reactions[0].products[:]
    model = CoreEdgeReactionModel()
    canonical = [model.make_new_reaction(r, check_existing=False, generate_thermo=False, generate_kinetics=False)[0]
                 for r in reactions]
    network = PDepNetwork(source=canonical[0].reactants[:])
    if control == 'no_electron':
        for reaction in canonical:
            network.add_path_reaction(reaction)
        assert network.path_reactions == [canonical[0]]
    else:
        from rmgpy.exceptions import NetworkError
        for reaction in canonical:
            with pytest.raises(NetworkError, match='Electron reactions cannot enter pressure-dependent networks'):
                network.add_path_reaction(reaction)
        assert not network.path_reactions


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('routed_first', (False, True))
def test_mixed_pdep_and_explicit_admission_refuses_channels(database, order, routed_first):
    from types import SimpleNamespace
    from rmgpy.exceptions import NetworkError

    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    reactions = unimolecular_pair(database, order)
    reactions[0].elementary_high_p = routed_first
    reactions[1].elementary_high_p = not routed_first
    if not routed_first:
        model.process_new_reactions([reactions[0]], None, generate_thermo=False, generate_kinetics=False)
    rejected = reactions[0] if routed_first else reactions[1]
    before = admission_state(model)
    with pytest.raises(NetworkError, match='Electron reactions cannot enter pressure-dependent networks'):
        model.process_new_reactions([rejected], None, generate_thermo=False, generate_kinetics=False)
    assert admission_state(model) == before


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('routed_first', (False, True))
@pytest.mark.parametrize('operation', ('seed', 'library'))
def test_seed_and_library_mixed_pdep_refusal(database, order, routed_first, operation, monkeypatch):
    from types import SimpleNamespace
    from rmgpy.exceptions import NetworkError
    import rmgpy.rmg.model as model_module

    monkeypatch.setattr(model_module, 'submit', controlled_thermo)
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    reactions = unimolecular_pair(database, order)
    reactions[0].elementary_high_p = routed_first
    reactions[1].elementary_high_p = not routed_first
    for reaction in reactions:
        entry = database.kinetics.libraries[reaction.library].entries[1]
        entry.item, entry.data = reaction, reaction.kinetics
    admission = model.add_seed_mechanism_to_core if operation == 'seed' else model.add_reaction_library_to_edge
    if not routed_first:
        admission(order[0])
    before = admission_state(model)
    with pytest.raises(NetworkError, match='Electron reactions cannot enter pressure-dependent networks'):
        admission(order[0] if routed_first else order[1])
    assert admission_state(model) == before


def source_state(reaction):
    species = list(dict.fromkeys(reaction.reactants + reaction.products))
    return (model_snapshot(ReactionModel(species=species, reactions=[reaction])),
            reaction.network_kinetics, repr(reaction.kinetics),
            tuple((s.index, s.label, tuple(atom.label for molecule in s.molecule for atom in molecule.atoms))
                  for s in species))


def assert_source_unchanged(reaction, state):
    species = list(dict.fromkeys(reaction.reactants + reaction.products))
    assert_model_unchanged(ReactionModel(species=species, reactions=[reaction]), state[0])
    assert reaction.network_kinetics is state[1]
    assert repr(reaction.kinetics) == state[2]
    assert tuple((s.index, s.label, tuple(atom.label for molecule in s.molecule for atom in molecule.atoms))
                 for s in species) == state[3]


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('operation', ('seed', 'library', 'saved_seed'))
def test_combined_electron_seed_preflight_is_atomic(database, order, operation, monkeypatch, tmp_path):
    from types import SimpleNamespace
    from rmgpy.data.base import Entry
    from rmgpy.data.kinetics.library import KineticsLibrary
    from rmgpy.exceptions import NetworkError
    import rmgpy.rmg.model as model_module

    monkeypatch.setattr(model_module, 'submit', controlled_thermo)
    reactions = unimolecular_pair(database, order)
    for reaction in reactions:
        reaction.elementary_high_p = reaction.library == 'ElectronAttachment'
    combined = KineticsLibrary(label='combined', name='combined', auto_generated=True)
    combined.entries = {i: Entry(index=i, label=r.to_labeled_str(), item=r, data=r.kinetics,
                                 long_desc='Originally from reaction library: ' + r.library + '\n')
                        for i, r in enumerate(reactions, 1)}
    if operation == 'saved_seed':
        path = tmp_path / 'combined'
        path.mkdir()
        combined.save(str(path / 'reactions.py'))
        database.kinetics.load_libraries(str(tmp_path), libraries=['combined'], additive=True)
    else:
        database.kinetics.libraries['combined'] = combined
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    # Even pending lists must keep their identities when the whole batch fails.
    model.new_reaction_list = [object()]
    model.new_species_list = [object()]
    before = admission_state(model)
    entries = list(database.kinetics.libraries['combined'].entries.values())
    for entry in entries:
        entry.data.comment = '  source comment must survive refusal  '
    source_states = [source_state(entry.item) for entry in entries]
    rate_states = [repr(entry.data) for entry in entries]
    method = model.add_reaction_library_to_edge if operation == 'library' else model.add_seed_mechanism_to_core
    with pytest.raises(NetworkError):
        method('combined')
    assert admission_state(model) == before
    for entry, state, rate_state in zip(entries, source_states, rate_states):
        assert_source_unchanged(entry.item, state)
        assert repr(entry.data) == rate_state


@pytest.mark.parametrize('order', ORDERS)
def test_corrected_electron_retry_admits_both_channels(database, order):
    from types import SimpleNamespace
    from rmgpy.exceptions import NetworkError

    first, second = unimolecular_pair(database, order)
    first.elementary_high_p = False
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    model.process_new_reactions([first], None, generate_thermo=False, generate_kinetics=False)
    with pytest.raises(NetworkError):
        model.process_new_reactions([second], None, generate_thermo=False, generate_kinetics=False)
    second.elementary_high_p = False
    model.process_new_reactions([second], None, generate_thermo=False, generate_kinetics=False)
    assert len(model.edge.reactions) == 2
    assert {r.library for r in model.edge.reactions} == set(order)
    assert model.reaction_counter == 2
    assert sorted(r.index for r in model.edge.reactions) == [1, 2]
    assert not model.network_list and not model.network_dict and model.network_count == 0


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('kind', ('library', 'template', 'molecule_template', 'cached_high_p', 'two_temperature', 'allow_route'))
def test_single_electron_routing_refuses_before_registration(database, order, kind):
    from types import SimpleNamespace
    from rmgpy.data.kinetics.family import KineticsFamily, TemplateReaction
    from rmgpy.exceptions import NetworkError

    reaction = unimolecular_pair(database, order)[0]
    source = next(s for s in reaction.reactants + reaction.products if s.label == 'O-')
    if kind in ('template', 'molecule_template'):
        database.kinetics.families['electron-fixture'] = KineticsFamily(label='electron-fixture')
        reaction = TemplateReaction(reactants=reaction.reactants[:], products=reaction.products[:],
                                    kinetics=reaction.kinetics, family='electron-fixture')
        if kind == 'molecule_template':
            reaction.reactants = [s.molecule[0] for s in reaction.reactants]
            reaction.products = [s.molecule[0] for s in reaction.products]
    elif kind == 'allow_route':
        reaction.elementary_high_p = False
        reaction.allow_pdep_route = True
    elif kind in ('cached_high_p', 'two_temperature'):
        reaction.kinetics = TwoTemperaturePlasma(A=(1, 'm^3/(mol*s)'), n=0,
                                                Ea_g=(0, 'J/mol'), Ea_e=(0, 'J/mol'))
        if kind == 'cached_high_p':
            reaction.elementary_high_p = False
            reaction.network_kinetics = Arrhenius(A=(1, 's^-1'), n=0, Ea=(0, 'J/mol'))
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    before = admission_state(model)
    participants = (reaction.reactants, tuple(reaction.reactants), reaction.products, tuple(reaction.products))
    index, network_kinetics = reaction.index, reaction.network_kinetics
    with pytest.raises(NetworkError, match='Electron reactions cannot enter pressure-dependent networks') as error:
        model.process_new_reactions([reaction], source, generate_thermo=False, generate_kinetics=False)
    assert admission_state(model) == before
    assert reaction.reactants is participants[0] and tuple(reaction.reactants) == participants[1]
    assert reaction.products is participants[2] and tuple(reaction.products) == participants[3]
    assert reaction.index == index and reaction.network_kinetics is network_kinetics
    assert str(reaction) in str(error.value)
    assert getattr(reaction, 'library', reaction.family) in str(error.value)


@pytest.mark.parametrize('order', ORDERS)
def test_process_electron_batch_refuses_before_any_registration(database, order):
    from types import SimpleNamespace
    from rmgpy.exceptions import NetworkError

    reactions = unimolecular_pair(database, order)
    reactions[0].elementary_high_p = False
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    before = admission_state(model)
    states = [source_state(r) for r in reactions]
    with pytest.raises(NetworkError):
        model.process_new_reactions(reactions, None, generate_thermo=False, generate_kinetics=False)
    assert admission_state(model) == before
    for reaction, state in zip(reactions, states):
        assert_source_unchanged(reaction, state)


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('operation', ('make_reaction', 'register', 'make_pdep', 'model_network', 'path', 'merge',
                                        'restore', 'configurations', 'update', 'explore'))
def test_electron_network_boundaries_refuse_without_mutation(database, order, operation, monkeypatch):
    from types import SimpleNamespace
    from rmgpy.exceptions import NetworkError
    from rmgpy.rmg.pdep import PDepNetwork
    from rmgpy.pdep import Configuration
    import rmgpy.rmg.pdep as pdep_module

    reaction = unimolecular_pair(database, order)[0]
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    source = next(s for s in reaction.reactants + reaction.products if s.label == 'O-')
    network = PDepNetwork(source=[source])
    other = PDepNetwork(source=[source])
    other.path_reactions = [reaction]  # Emulate a legacy serialized network.
    job = SimpleNamespace(network='unchanged')
    if operation == 'update':
        from arkane.pdep import PressureDependenceJob

        job = PressureDependenceJob(network=None, Tmin=(300, 'K'), Tmax=(1000, 'K'), Tcount=2,
                                    Pmin=(0.1, 'bar'), Pmax=(10, 'bar'), Pcount=2,
                                    method='modified strong collision', interpolationModel=('Chebyshev', 2, 2))
        job.network = 'unchanged'
        job.output_file = None
        # A complete, valid legacy network returns normally without the guard.
        network.valid = True
    isomer = next(s for s in reaction.reactants + reaction.products if s.label == 'O')
    if operation == 'explore':
        network.products = [Configuration(isomer)]
        monkeypatch.setattr(pdep_module, 'react_species', lambda reactants: [reaction])
    before, state = admission_state(model), source_state(reaction)
    network_state = {k: v[:] if isinstance(v, list) else v.copy() if isinstance(v, dict) else v
                     for k, v in network.__dict__.items()}
    if operation in ('configurations', 'update'):
        network.path_reactions = [reaction]
        network_state = {k: v[:] if isinstance(v, list) else v.copy() if isinstance(v, dict) else v
                     for k, v in network.__dict__.items()}
    with pytest.raises(NetworkError, match='Electron reactions cannot enter pressure-dependent networks'):
        if operation == 'make_reaction':
            model.make_new_reaction(reaction, generate_thermo=False, generate_kinetics=False)
        elif operation == 'register':
            model.register_reaction(reaction)
        elif operation == 'make_pdep':
            model.make_new_pdep_reaction(reaction)
        elif operation == 'model_network':
            model.add_reaction_to_unimolecular_networks(reaction, source)
        elif operation == 'path':
            network.add_path_reaction(reaction)
        elif operation == 'merge':
            network.merge(other)
        elif operation == 'restore':
            network.__setstate__({'index': 999, 'path_reactions': [reaction]})
        elif operation == 'configurations':
            network.update_configurations(model)
        elif operation == 'explore':
            network.explore_isomer(isomer)
        else:
            network.update(model, job)
    assert admission_state(model) == before
    assert_source_unchanged(reaction, state)
    assert network.__dict__ == network_state
    assert job.network == 'unchanged'


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('control', ('pdep_off', 'atom_limit'))
def test_explicit_electron_controls_remain_allowed(database, order, control):
    from types import SimpleNamespace

    reactions = unimolecular_pair(database, order)
    model = CoreEdgeReactionModel()
    if control == 'atom_limit':
        for reaction in reactions:
            reaction.elementary_high_p = False
        model.pressure_dependence = SimpleNamespace(maximum_atoms=0)
    model.process_new_reactions(reactions, None, generate_thermo=False, generate_kinetics=False)
    assert len(model.edge.reactions) == 2
    assert model.reaction_counter == 2
    assert not model.network_list


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('ineligible', ('atom_limit', 'bimolecular'))
def test_declared_electron_routes_refuse_even_when_ineligible(database, order, ineligible):
    from types import SimpleNamespace
    from rmgpy.exceptions import NetworkError

    reaction = (unimolecular_pair(database, order)[0] if ineligible == 'atom_limit'
                else library_pair(database, order)[0])
    reaction.elementary_high_p = True
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=0 if ineligible == 'atom_limit' else None)
    before, state = admission_state(model), source_state(reaction)
    with pytest.raises(NetworkError, match='Electron reactions cannot enter pressure-dependent networks'):
        model.process_new_reactions([reaction], None, generate_thermo=False, generate_kinetics=False)
    assert admission_state(model) == before
    assert_source_unchanged(reaction, state)


def generated_attachment(database, monkeypatch):
    """Return O + e- => O- with real family kinetics defined in reverse."""
    from rmgpy.data.kinetics.family import KineticsFamily, TemplateReaction

    species = database.kinetics.libraries['ElectronDetachment'].entries[1].item
    participants = {s.label: s for s in species.reactants + species.products}
    other = database.kinetics.libraries['ElectronAttachment'].entries[1].item
    participants.update({s.label: s for s in other.reactants + other.products})
    reaction = TemplateReaction(reactants=[participants['O'], participants['e-']],
                                products=[participants['O-']], family='electron-fixture')
    family = KineticsFamily(label=reaction.family)
    database.kinetics.families[reaction.family] = family
    monkeypatch.setattr(family, 'get_kinetics', lambda *args, **kwargs:
                        (Arrhenius(A=(1, 's^-1'), n=0, Ea=(0, 'J/mol')), 'rate rules', None, False))
    for participant in reaction.reactants + reaction.products:
        controlled_thermo(participant)
    return reaction


@pytest.mark.parametrize('order', ORDERS)
def test_generated_electron_direction_flip_refuses_before_registration(database, order, monkeypatch):
    from types import SimpleNamespace
    from rmgpy.exceptions import NetworkError

    reaction = generated_attachment(database, monkeypatch)
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=1)
    before, state = admission_state(model), source_state(reaction)
    with pytest.raises(NetworkError, match='Electron reactions cannot enter pressure-dependent networks') as error:
        model.process_new_reactions([reaction], None, generate_thermo=False)
    assert admission_state(model) == before
    assert_source_unchanged(reaction, state)
    assert str(reaction) in str(error.value)
    assert reaction.family in str(error.value)


@pytest.mark.parametrize('order', ORDERS)
def test_generated_electron_direction_flip_retry_matches_clean_run(database, order, monkeypatch):
    from types import SimpleNamespace
    from rmgpy.exceptions import NetworkError

    reaction = generated_attachment(database, monkeypatch)
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=1)
    with pytest.raises(NetworkError):
        model.process_new_reactions([reaction], None, generate_thermo=False)
    clean_reaction = generated_attachment(database, monkeypatch)
    clean = CoreEdgeReactionModel()
    # Lowering the atom limit cannot bypass the direction-independent refusal.
    for candidate_model, candidate in ((model, reaction), (clean, clean_reaction)):
        candidate_model.pressure_dependence = SimpleNamespace(maximum_atoms=0)
        try:
            candidate_model.process_new_reactions([candidate], None, generate_thermo=False)
        except NetworkError:
            pass
    assert registry_contents(model) == registry_contents(clean)
    assert model.species_counter == model.reaction_counter == 0
    # A corrected retry disables pressure dependence and must match clean admission.
    for candidate_model, candidate in ((model, reaction), (clean, clean_reaction)):
        candidate_model.pressure_dependence = None
        candidate_model.process_new_reactions([candidate], None, generate_thermo=False)
    assert registry_contents(model) == registry_contents(clean)
    assert model.species_counter == 3 and model.reaction_counter == 1
    assert len(model.edge.reactions) == 1
    assert len(model.edge.reactions[0].reactants) == 1  # Family kinetics actually flipped it.


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('filter_kind', ('atom_limit', 'shape', 'unreal_group'))
def test_generated_electron_filters_cannot_bypass_preflight(database, order, filter_kind, monkeypatch):
    from types import SimpleNamespace
    from rmgpy.exceptions import NetworkError
    from rmgpy.molecule.group import Group

    reaction = generated_attachment(database, monkeypatch)
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=0 if filter_kind == 'atom_limit' else None)
    if filter_kind == 'shape':
        reaction.products.append(reaction.products[0])  # Neither direction is unimolecular.
    elif filter_kind == 'unreal_group':
        model.unrealgroups = [Group().from_adjacency_list('1 O ux px cx')]
    before, state = admission_state(model), source_state(reaction)
    with pytest.raises(NetworkError, match='Electron reactions cannot enter pressure-dependent networks'):
        model.process_new_reactions([reaction], None, generate_thermo=False)
    assert admission_state(model) == before
    assert_source_unchanged(reaction, state)


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('boundary', ('make', 'register'))
@pytest.mark.parametrize('filter_kind', ('atom_limit', 'shape', 'unreal_group'))
def test_explicit_template_registration_is_not_automatic_routing(database, order, boundary, filter_kind, monkeypatch):
    from types import SimpleNamespace
    from rmgpy.molecule.group import Group

    reaction = generated_attachment(database, monkeypatch)
    reaction.kinetics = Arrhenius(A=(1, 'm^3/(mol*s)'), n=0, Ea=(0, 'J/mol'))
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=0 if filter_kind == 'atom_limit' else None)
    if filter_kind == 'shape':
        reaction.products.append(reaction.products[0])
    elif filter_kind == 'unreal_group':
        model.unrealgroups = [Group().from_adjacency_list('1 O ux px cx')]
    if boundary == 'make':
        admitted, is_new = model.make_new_reaction(reaction, generate_thermo=False, generate_kinetics=False)
        assert admitted is reaction and is_new
        assert model.species_counter == 3 and model.reaction_counter == 1
    else:
        model.register_reaction(reaction)
        assert model.reaction_dict
    assert not model.network_dict and not model.network_list


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('operation', ('seed', 'library'))
def test_stored_template_electron_rate_remains_explicit(database, order, operation, monkeypatch):
    from types import SimpleNamespace
    from rmgpy.data.base import Entry
    from rmgpy.data.kinetics.family import KineticsFamily, TemplateReaction
    from rmgpy.data.kinetics.library import KineticsLibrary
    from rmgpy.reaction import Reaction
    import rmgpy.rmg.model as model_module

    monkeypatch.setattr(model_module, 'submit', controlled_thermo)
    electron = library_pair(database, order)[0]
    participants = {s.label: s for s in electron.reactants + electron.products}
    ordinary = Reaction(reactants=[participants['O'], participants['O']], products=[participants['O2']],
                        reversible=False)
    family = 'electron-fixture'
    database.kinetics.families[family] = KineticsFamily(label=family)
    stored = KineticsLibrary(label='stored-template', name='stored-template', auto_generated=True)
    template_entry = Entry(index=2, item=electron, data=Arrhenius(A=(1, 'm^3/(mol*s)'), n=0, Ea=(0, 'J/mol')),
                           long_desc='Estimated using rate rule [fixture]\nfamily: ' + family)
    entries = [Entry(index=1, item=ordinary, data=Arrhenius(A=(1, 'm^3/(mol*s)'), n=0, Ea=(0, 'J/mol'))),
               template_entry]
    if order == ORDERS[1]:
        entries.reverse()
    stored.entries = {entry.index: entry for entry in entries}
    database.kinetics.libraries[stored.label] = stored
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=0)
    if operation == 'seed':
        model.add_seed_mechanism_to_core(stored.label)
    else:
        model.add_reaction_library_to_edge(stored.label)
    reactions = model.core.reactions if operation == 'seed' else model.edge.reactions
    assert len(reactions) == model.reaction_counter == 2
    recovered = next(r for r in reactions if isinstance(r, TemplateReaction))
    assert recovered.entry is template_entry and recovered.kinetics is template_entry.data
    assert recovered.is_isomorphic(electron, either_direction=False)
    assert not model.network_dict and not model.network_list


def real_metadata_attachment(database):
    """Generate canonical O2 => O2- using the shipped attachment family."""
    from rmgpy.data.kinetics.family import TemplateReaction
    from rmgpy.electron_balance import get_electron_placement_counts
    from rmgpy.species import Species

    name = 'Plasma_Electron_Attachment'
    families_path = str(Path(settings['database.directory']) / 'kinetics/families')
    database.kinetics.load_families(families_path, families=[name])
    family = database.kinetics.families[name]
    family.add_rules_from_training(thermo_database=None)
    family.fill_rules_by_averaging_up()
    oxygen = Species(label='O2').from_smiles('[O][O]')
    oxygen.generate_resonance_structures()
    reactions = database.kinetics.generate_reactions_from_families([oxygen], only_families=[name])
    assert len(reactions) == 1
    reaction = reactions[0]
    assert isinstance(reaction, TemplateReaction) and reaction.family == name
    assert reaction.electrons == -1
    assert get_electron_placement_counts(reaction) == (1, 0)
    assert len(reaction.reactants) == len(reaction.products) == 1
    assert not any(s.is_electron() for s in reaction.reactants + reaction.products)
    assert reaction.reactants[0].is_isomorphic(Species().from_smiles('[O][O]'))
    assert reaction.products[0].is_isomorphic(Species().from_smiles('[O-][O]'))
    reaction.kinetics = family.get_kinetics_for_template(
        family.retrieve_template(reaction.template), degeneracy=reaction.degeneracy)[0]
    for participant in reaction.reactants + reaction.products:
        controlled_thermo(participant)
    return reaction


@pytest.mark.database
@pytest.mark.parametrize('order', ORDERS)
def test_metadata_electron_real_family_refuses_before_registration(database, order):
    from types import SimpleNamespace
    from rmgpy.exceptions import NetworkError

    reaction = real_metadata_attachment(database)
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    before, state = admission_state(model), source_state(reaction)
    with pytest.raises(NetworkError, match='Electron reactions cannot enter pressure-dependent networks') as error:
        model.process_new_reactions([reaction], reaction.reactants[0], generate_thermo=False, generate_kinetics=False)
    assert admission_state(model) == before
    assert_source_unchanged(reaction, state)
    assert reaction.electrons == -1
    assert str(reaction) in str(error.value) and reaction.family in str(error.value)


@pytest.mark.database
@pytest.mark.parametrize('order', ORDERS)
def test_metadata_electron_real_family_retry_matches_clean_run(database, order, monkeypatch):
    from types import SimpleNamespace
    from rmgpy.exceptions import NetworkError

    import rmgpy.rmg.model as model_module

    def named_thermo(species, solvent_name=None):
        controlled_thermo(species)
        species.thermo.label = 'library_' + species.molecule[0].get_formula() + '_' + str(
            species.molecule[0].get_net_charge())

    monkeypatch.setattr(model_module, 'submit', named_thermo)
    reaction = real_metadata_attachment(database)
    for participant in reaction.reactants + reaction.products:
        assert not participant.label  # Preserve the real family's unnamed representation.
        participant.thermo = None
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    with pytest.raises(NetworkError):
        model.process_new_reactions([reaction], reaction.reactants[0], generate_kinetics=False)
    clean_reaction = real_metadata_attachment(database)
    for participant in clean_reaction.reactants + clean_reaction.products:
        assert not participant.label
        participant.thermo = None
    clean = CoreEdgeReactionModel()
    for candidate_model, candidate in ((model, reaction), (clean, clean_reaction)):
        candidate_model.pressure_dependence = None
        candidate_model.process_new_reactions([candidate], candidate.reactants[0], generate_kinetics=False)
    assert registry_contents(model) == registry_contents(clean)
    assert model.species_counter == 2 and model.reaction_counter == 1
    assert len(model.edge.reactions) == 1 and not model.network_list
    assert model.edge.reactions[0].electrons == -1
    assert {s.label for s in model.edge.species} == {'library_O2_0', 'library_O2_-1'}


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('representation', (
    'metadata_consumed', 'metadata_produced', 'declared_family_zero',
    'declared_library_zero', 'declared_conserved_zero',
))
@pytest.mark.parametrize('boundary', ('process', 'make', 'path', 'update'))
def test_metadata_electron_representations_refuse_without_mutation(
        database, order, representation, boundary, monkeypatch):
    from types import SimpleNamespace
    from arkane.pdep import PressureDependenceJob
    from rmgpy.data.kinetics.family import TemplateReaction
    from rmgpy.electron_placement import FAMILY_ELECTRON_PLACEMENT
    from rmgpy.exceptions import NetworkError
    from rmgpy.rmg.pdep import PDepNetwork
    from rmgpy.species import Species

    reaction = TemplateReaction(reactants=[Species(label='O2').from_smiles('[O][O]')],
                                products=[Species(label='O2-').from_smiles('[O-][O]')],
                                family='metadata-fixture',
                                kinetics=Arrhenius(A=(1, 's^-1'), n=0, Ea=(0, 'J/mol')))
    if representation == 'metadata_consumed':
        reaction.electrons = -1
    elif representation == 'metadata_produced':
        reaction.reactants, reaction.products = reaction.products, reaction.reactants
        reaction.electrons = 1
    elif representation == 'declared_family_zero':
        reaction.family = 'Plasma_Electron_Attachment'
    elif representation == 'declared_library_zero':
        # A generated seed's current family is silent; its provenance still declares electrons.
        from rmgpy.data.kinetics.library import LibraryReaction
        reaction = LibraryReaction(reactants=reaction.reactants, products=reaction.products,
                                   library='PlasmaRadiativeRecombination',
                                   elementary_high_p=True, kinetics=reaction.kinetics)
        reaction.family = 'seed'
    else:
        monkeypatch.setitem(FAMILY_ELECTRON_PLACEMENT, reaction.family, (1, 1))
    reaction.elementary_high_p = True
    from rmgpy.data.kinetics.family import KineticsFamily
    database.kinetics.families[reaction.family] = KineticsFamily(label=reaction.family)
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    network = PDepNetwork(source=reaction.reactants[:])
    job = PressureDependenceJob(network=None, Tmin=(300, 'K'), Tmax=(1000, 'K'), Tcount=2,
                               Pmin=(0.1, 'bar'), Pmax=(10, 'bar'), Pcount=2,
                               method='modified strong collision', interpolationModel=('Chebyshev', 2, 2))
    job.network = 'unchanged'
    job.output_file = None
    if boundary == 'update':
        network.path_reactions = [reaction]
        network.valid = True
    before, state = admission_state(model), source_state(reaction)
    net_before = {k: v[:] if isinstance(v, list) else v.copy() if isinstance(v, dict) else v
                  for k, v in network.__dict__.items()}
    electrons = reaction.electrons
    with pytest.raises(NetworkError, match='Electron reactions cannot enter pressure-dependent networks') as error:
        if boundary == 'process':
            model.process_new_reactions([reaction], reaction.reactants[0], generate_thermo=False, generate_kinetics=False)
        elif boundary == 'make':
            model.make_new_reaction(reaction, generate_thermo=False, generate_kinetics=False)
        elif boundary == 'path':
            network.add_path_reaction(reaction)
        else:
            network.update(model, job)
    assert str(reaction) in str(error.value)
    assert (getattr(reaction, 'library', None) or reaction.family) in str(error.value)
    assert admission_state(model) == before and network.__dict__ == net_before
    assert job.network == 'unchanged' and reaction.electrons == electrons
    assert_source_unchanged(reaction, state)


def parsed_specific_collider(database, order, electron=True):
    """Read the production Chemkin representation with a specific collider."""
    from rmgpy.chemkin import read_kinetics_entry
    from rmgpy.data.kinetics.library import LibraryReaction
    from rmgpy.species import Species

    oxygen = Species(label='O').from_smiles('[O]')
    dioxygen = Species(label='O2').from_smiles('[O][O]')
    if electron:
        collider = next(s for r in library_pair(database, order)
                        for s in r.reactants + r.products if s.is_electron()).copy(deep=True)
        collider.label = 'e-'
    else:
        collider = Species(label='N2').from_smiles('N#N')
    for species in (oxygen, dioxygen, collider):
        controlled_thermo(species)
    label = collider.label
    parsed = read_kinetics_entry(
        'O2(+{0})=>O+O(+{0}) 1.0e12 0.0 0.0'.format(label),
        {'O': oxygen, 'O2': dioxygen, label: collider},
        ['', 's^-1', 'cm^3/(mol*s)', 'cm^6/(mol^2*s)', 'cm^9/(mol^3*s)'], [], 'kcal/mol')
    assert parsed.specific_collider is collider and parsed.electrons == 0
    assert not any(s.is_electron() for s in parsed.reactants + parsed.products)
    return LibraryReaction(reactants=parsed.reactants, products=parsed.products,
                           specific_collider=parsed.specific_collider, kinetics=parsed.kinetics,
                           reversible=parsed.reversible, library=order[0], elementary_high_p=True)


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('operation', _COLLIDER_BOUNDARIES + ('seed', 'library', 'saved_seed'))
def test_collider_electron_boundaries_refuse_before_mutation(
        database, order, operation, monkeypatch, tmp_path):
    from types import SimpleNamespace
    from arkane.pdep import PressureDependenceJob
    from rmgpy.data.base import Entry
    from rmgpy.data.kinetics.library import KineticsLibrary
    from rmgpy.exceptions import NetworkError
    from rmgpy.pdep import Configuration
    from rmgpy.rmg.pdep import PDepNetwork
    import rmgpy.rmg.model as model_module
    import rmgpy.rmg.pdep as pdep_module

    monkeypatch.setattr(model_module, 'submit', controlled_thermo)
    reaction = parsed_specific_collider(database, order)
    library_name = order[0]
    if operation in ('seed', 'library', 'saved_seed'):
        library_name = 'electron-collider'
        stored = KineticsLibrary(label=library_name, name=library_name)
        stored.entries = {1: Entry(index=1, label=reaction.to_labeled_str(), item=reaction, data=reaction.kinetics)}
        if operation == 'saved_seed':
            path = tmp_path / library_name
            path.mkdir()
            stored.save(str(path / 'reactions.py'))
            database.kinetics.load_libraries(str(tmp_path), libraries=[library_name], additive=True)
            entries = database.kinetics.libraries[library_name].entries
            assert len(entries) == 1
            reaction = next(iter(entries.values())).item
            assert reaction.specific_collider.is_electron() and reaction.elementary_high_p
            assert not any(s.is_electron() for s in reaction.reactants + reaction.products)
        else:
            database.kinetics.libraries[library_name] = stored
    collider = reaction.specific_collider
    collider_state = (collider.index, collider.label, collider.molecule[0].to_adjacency_list())
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    source, isomer = reaction.reactants[0], reaction.products[0]
    network, other = PDepNetwork(source=[source]), PDepNetwork(source=[source])
    other.path_reactions = [reaction]
    job = PressureDependenceJob(network=None, Tmin=(300, 'K'), Tmax=(1000, 'K'), Tcount=2,
                               Pmin=(0.1, 'bar'), Pmax=(10, 'bar'), Pcount=2,
                               method='modified strong collision', interpolationModel=('Chebyshev', 2, 2))
    job.network, job.output_file = 'unchanged', None
    if operation in ('configurations', 'update'):
        network.path_reactions = [reaction]
        network.valid = True
    if operation == 'explore':
        network.products = [Configuration(isomer)]
        monkeypatch.setattr(pdep_module, 'react_species', lambda reactants: [reaction])
    before, state = admission_state(model), source_state(reaction)
    network_state = {k: v[:] if isinstance(v, list) else v.copy() if isinstance(v, dict) else v
                     for k, v in network.__dict__.items()}
    with pytest.raises(NetworkError, match='Electron reactions cannot enter pressure-dependent networks') as error:
        if operation == 'process':
            model.process_new_reactions([reaction], source, generate_thermo=False, generate_kinetics=False)
        elif operation == 'make_reaction':
            model.make_new_reaction(reaction, generate_thermo=False, generate_kinetics=False)
        elif operation == 'register':
            model.register_reaction(reaction)
        elif operation == 'make_pdep':
            model.make_new_pdep_reaction(reaction)
        elif operation == 'model_network':
            model.add_reaction_to_unimolecular_networks(reaction, source)
        elif operation == 'path':
            network.add_path_reaction(reaction)
        elif operation == 'merge':
            network.merge(other)
        elif operation == 'restore':
            network.__setstate__({'index': 999, 'path_reactions': [reaction]})
        elif operation == 'configurations':
            network.update_configurations(model)
        elif operation == 'update':
            network.update(model, job)
        elif operation == 'explore':
            network.explore_isomer(isomer)
        elif operation == 'library':
            model.add_reaction_library_to_edge(library_name)
        else:
            model.add_seed_mechanism_to_core(library_name)
    assert library_name in str(error.value) and 'e-' in str(error.value)
    assert admission_state(model) == before and network.__dict__ == network_state
    assert job.network == 'unchanged' and reaction.specific_collider is collider
    assert_source_unchanged(reaction, state)
    assert (collider.index, collider.label, collider.molecule[0].to_adjacency_list()) == collider_state


@pytest.mark.parametrize('order', ORDERS)
def test_ordinary_specific_collider_control(database, order):
    from types import SimpleNamespace

    reaction = parsed_specific_collider(database, order, electron=False)
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    model.process_new_reactions([reaction], reaction.reactants[0], generate_thermo=False, generate_kinetics=False)
    assert model.species_counter == 3 and model.reaction_counter == 1
    assert len(model.network_list) == 1 and model.network_list[0].path_reactions == [reaction]
    assert not reaction.specific_collider.is_electron()


def reconstructed_owner_library(database, provenance, ordinary_first):
    """Keep ownership only in long_desc, matching seed/restart reconstruction."""
    from rmgpy.data.base import Entry
    from rmgpy.data.kinetics.library import KineticsLibrary
    from rmgpy.electron_balance import get_placement_owner
    from rmgpy.reaction import Reaction
    from rmgpy.species import Species

    oxygen = Species(label='O').from_smiles('[O]')
    dioxygen = Species(label='O2').from_smiles('[O][O]')
    ordinary = Reaction(reactants=[oxygen, oxygen], products=[dioxygen], reversible=False)
    routed = Reaction(reactants=[dioxygen], products=[oxygen, oxygen], reversible=False)
    routed.elementary_high_p = True
    assert routed.electrons == 0 and get_placement_owner(routed) is None
    stored = KineticsLibrary(label='reconstructed-owner', name='reconstructed-owner', auto_generated=True)
    if provenance == 'library':
        owner = 'PlasmaRadiativeRecombination'
        description = 'Originally from reaction library: ' + owner
        database.kinetics.libraries[owner] = KineticsLibrary(label=owner, name=owner)
    else:
        owner = 'Plasma_Electron_Attachment'
        assert owner not in database.kinetics.families
        description = 'Estimated using rate rule [fixture]\nfamily: ' + owner
    entries = [Entry(index=1, item=ordinary, data=Arrhenius(A=(1, 'm^3/(mol*s)'), n=0, Ea=(0, 'J/mol'))),
               Entry(index=2, item=routed, data=Arrhenius(A=(2, 's^-1'), n=0, Ea=(0, 'J/mol')),
                     long_desc=description)]
    for entry in entries:
        entry.data.comment = '  source comment must survive refusal  '
    if not ordinary_first:
        entries.reverse()
    stored.entries = {entry.index: entry for entry in entries}
    database.kinetics.libraries[stored.label] = stored
    return stored


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('operation', ('seed', 'library'))
@pytest.mark.parametrize('provenance', ('library', 'unloaded_family'))
@pytest.mark.parametrize('ordinary_first', (True, False))
def test_reconstructed_electron_owner_preflight_is_atomic(
        database, order, operation, provenance, ordinary_first, monkeypatch):
    from types import SimpleNamespace
    from rmgpy.exceptions import NetworkError
    import rmgpy.rmg.model as model_module

    monkeypatch.setattr(model_module, 'submit', controlled_thermo)
    stored = reconstructed_owner_library(database, provenance, ordinary_first)
    getter = stored.get_library_reactions
    calls = []

    def reconstruct():
        calls.append(True)
        return getter()

    monkeypatch.setattr(stored, 'get_library_reactions', reconstruct)
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    model.new_species_list, model.new_reaction_list = [object()], [object()]
    before = admission_state(model)
    states = [(entry, source_state(entry.item), repr(entry.data)) for entry in stored.entries.values()]
    method = model.add_seed_mechanism_to_core if operation == 'seed' else model.add_reaction_library_to_edge
    with pytest.raises(NetworkError, match='Electron reactions cannot enter pressure-dependent networks') as error:
        method(stored.label)
    assert admission_state(model) == before
    assert stored.label in str(error.value) and calls == [True]
    for entry, state, rate in states:
        assert_source_unchanged(entry.item, state)
        assert repr(entry.data) == rate


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('operation', ('seed', 'library'))
@pytest.mark.parametrize('provenance', ('library', 'unloaded_family'))
@pytest.mark.parametrize('ordinary_first', (True, False))
def test_reconstructed_electron_owner_retry_matches_clean(
        database, order, operation, provenance, ordinary_first, monkeypatch):
    from types import SimpleNamespace
    from rmgpy.exceptions import NetworkError
    import rmgpy.rmg.model as model_module

    monkeypatch.setattr(model_module, 'submit', controlled_thermo)
    stored = reconstructed_owner_library(database, provenance, ordinary_first)
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    model.new_species_list, model.new_reaction_list = [object()], [object()]
    before = admission_state(model)
    method_name = 'add_seed_mechanism_to_core' if operation == 'seed' else 'add_reaction_library_to_edge'
    with pytest.raises(NetworkError):
        getattr(model, method_name)(stored.label)
    refused = admission_state(model)
    stored.entries[2].item.elementary_high_p = False
    getattr(model, method_name)(stored.label)
    clean_stored = reconstructed_owner_library(database, provenance, ordinary_first)
    clean_stored.entries[2].item.elementary_high_p = False
    clean = CoreEdgeReactionModel()
    clean.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    getattr(clean, method_name)(clean_stored.label)
    expected = (2, 2, 0, 0) if operation == 'seed' else (0, 0, 2, 2)
    assert model.get_model_size() == clean.get_model_size() == expected
    assert registry_contents(model) == registry_contents(clean)
    part, clean_part = ((model.core, clean.core) if operation == 'seed' else (model.edge, clean.edge))
    assert [repr(r.kinetics) for r in part.reactions] == [repr(r.kinetics) for r in clean_part.reactions]
    assert not model.network_list and before == refused
