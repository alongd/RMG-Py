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

"""Compiled cascade binding, single collision energy, and promotion refusals."""

import copy
import math
from types import SimpleNamespace

import numpy as np
import pytest

import rmgpy.constants as constants
import rmgpy.data.rmg as rmg_data_module
from rmgpy.data.base import Entry
from rmgpy.exceptions import PlasmaStateError, NonEquilibriumReverseRateError
from rmgpy.kinetics import EEDFChannel
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.species import Species
from rmgpy.tools.eedf.schema import _file_hash
from eedfBranchedTest import (
    ID6, branch_fixture, source_files, branch_reactions, family_from_reactions,
)
from plasmaEEDFCoreTest import artifact, thermo, TGAS, PRESSURE, VOLUME, LAMBDA, MU0, build_reactor


@pytest.fixture(autouse=True)
def isolate_database():
    saved = rmg_data_module.database
    rmg_data_module.database = None
    try:
        yield
    finally:
        rmg_data_module.database = saved


def cascade_model(tmp_path, branch_fixture, *, placement='core', mutate=None, zeros=False):
    snapshot, reactions, family = branch_fixture
    if zeros:
        snapshot['channels'][0]['product_yields'] = dict(zip(ID6, [1., 0., 0., 0.]))
        reactions = branch_reactions(snapshot)
        family = family_from_reactions(snapshot, reactions)
    electron, ground = reactions[0].reactants
    ion = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    ground.thermo, ion.thermo = thermo(), thermo(15.76)
    products = [reaction.products[1] for reaction in reactions]
    for species, energy in zip(products, [11.548, 11.624, 11.723, 11.828]):
        species.thermo = thermo(energy)
    core = [electron, ground, ion] + products
    path = artifact(tmp_path, reactions[0].kinetics, composition=False,
                    include_ar4s_envelope=False, family=family)
    if mutate is not None:
        mutate(core, reactions)
    # Real in-memory library provenance, including deliberately excessive
    # thermo in the radiation-refusal case. Never bypass the state guard.
    from rmgpy.data.thermo import ThermoDatabase, ThermoLibrary
    database = ThermoDatabase()
    library = ThermoLibrary(label='CascadeFixture')
    library.thermo_convention = 'ion'
    for i, species in enumerate(core[1:]):
        species.thermo.Cp0 = (20., 'J/(mol*K)')
        species.thermo.CpInf = (20., 'J/(mol*K)')
        library.entries[i] = Entry(index=i, label=species.label,
                                   item=species.molecule[0].copy(deep=True), data=copy.deepcopy(species.thermo))
    database.libraries = {library.label: library}
    database.library_order = [library.label]
    rmg_data_module.database = SimpleNamespace(thermo=database, solvation=None)
    initial = {species: 0. for species in core}
    initial.update({electron: 1e-6, ion: 1e-6, ground: .999998})
    reactor = PlasmaReactor(
        (TGAS, 'K'), (PRESSURE, 'Pa'), initial, (TGAS, 'K'), termination=[],
        diffusion_length=(LAMBDA, 'm'), ion_reduced_mobility=(MU0, 'm^2/(V*s)'),
        wall_single_bath_approximation=True,
        wall_neutralization_products={'Ar+': 'Ar'}, thermo_source_assertions={'Ar+': 'ion'},
        electron_energy_balance={'absorbed_power': (.5, 'W'), 'chamber_volume': (VOLUME, 'm^3'),
                                 'sheath': 'floating_wall'},
        electron_kinetics={'provider': 'loki-table', 'table': (str(path), _file_hash(path / 'table.h5')),
                           'branch': 'branch_0', 'initial_reduced_field': (1., 'Td'), 'empirical_laws': {}})
    if placement == 'core':
        core_reactions, edge_reactions = reactions, []
    elif placement == 'edge':
        core_reactions, edge_reactions = [], reactions
    else:
        core_reactions, edge_reactions = reactions[:2], reactions[2:]
    reactor.initialize_model(core, core_reactions, [], edge_reactions)
    return reactor, core, reactions


def test_id6_literal_rates_and_single_collision_power(tmp_path, branch_fixture):
    reactor, core, reactions = cascade_model(tmp_path, branch_fixture)
    reactor.residual(0., reactor.y0.copy(), np.zeros_like(reactor.y0))
    rates = np.asarray(reactor.kf)[:4] / constants.Na
    assert rates == pytest.approx([7.4204946996466445e-16, 2.120141342756184e-16,
                                  3.8476639183352967e-17, 7.45975657636435e-18], rel=2e-14, abs=0.)
    assert math.fsum(rates) == pytest.approx(1e-15, rel=2e-14, abs=0.)
    y = reactor.y0
    volume = reactor.compute_volume(y)
    n_e = y[reactor.species_index[core[0]]] * constants.Na / volume
    n_ar = y[reactor.species_index[core[1]]] * constants.Na / volume
    expected = 1e-15 * n_e * n_ar * 12.907 * constants.e * volume
    budget = reactor.energy_budget
    assert budget['Q_inelastic'] == pytest.approx(expected, rel=2e-14, abs=0.)
    assert budget['Q_inelastic_channels'][0] == pytest.approx(expected, rel=2e-14, abs=0.)
    assert budget['Q_inelastic_by_reaction_gross'] == pytest.approx(expected * np.array(list(ID6.values())))
    assert sum(budget['Q_inelastic_by_reaction']) == pytest.approx(expected, rel=2e-14, abs=0.)
    assert budget['Q_superelastic_channels'][0] == 0.
    radiation = []
    for reaction in reactions:
        dh = (reaction.products[1].thermo.get_enthalpy(TGAS)
              - reaction.reactants[1].thermo.get_enthalpy(TGAS)) / (constants.Na * constants.e)
        radiation.append(12.907 - dh)
    assert radiation == pytest.approx([1.359, 1.283, 1.184, 1.079], abs=1e-12)


def test_complete_edge_family_binds(tmp_path, branch_fixture):
    reactor, _, _ = cascade_model(tmp_path, branch_fixture, placement='edge')
    assert len(reactor.eedf_reaction_map) == 4
    assert [reactor.eedf_provider.reaction_yield(i) for i in range(4)] == list(ID6.values())


def test_split_core_edge_family_refuses(tmp_path, branch_fixture):
    with pytest.raises(PlasmaStateError, match='split cascade family'):
        cascade_model(tmp_path, branch_fixture, placement='split')


@pytest.mark.parametrize('zero_yield', [False, True])
def test_product_thermo_above_threshold_refuses(tmp_path, branch_fixture, zero_yield):
    def bad_thermo(core, reactions):
        reactions[-1].products[1].thermo = thermo(12.908)
    with pytest.raises(PlasmaStateError, match='cascade radiated energy is negative'):
        cascade_model(tmp_path, branch_fixture, mutate=bad_thermo, zeros=zero_yield)


def test_zero_branches_have_zero_rate_and_power(tmp_path, branch_fixture):
    reactor, _, _ = cascade_model(tmp_path, branch_fixture, zeros=True)
    reactor.residual(0., reactor.y0.copy(), np.zeros_like(reactor.y0))
    assert list(np.asarray(reactor.kf)[1:4]) == [0., 0., 0.]
    assert list(reactor.energy_budget['Q_inelastic_by_reaction'][1:4]) == [0., 0., 0.]


@pytest.mark.parametrize('problem', ['missing', 'extra', 'source', 'product', 'unweighted', 'reversible', 'yield'])
def test_reactor_refuses_bad_branch_deck(tmp_path, branch_fixture, problem):
    def mutate(core, reactions):
        if problem == 'missing':
            reactions.pop()
        elif problem == 'extra':
            reactions.append(copy.deepcopy(reactions[0]))
        elif problem == 'source':
            reactions[0].entry = Entry(index=9999)
        elif problem == 'product':
            reactions[0].products[1].molecule[0].electronic_state = '3P1'
        elif problem == 'reversible':
            reactions[0].reversible = True
        elif problem == 'yield':
            kin = reactions[0].kinetics
            reactions[0].kinetics = EEDFChannel(
                process=kin.process, collision_set=kin.collision_set,
                branch_product=kin.branch_product, branch_yield=kin.branch_yield + .001,
                yield_table_sha256=kin.yield_table_sha256)
        else:
            kin = reactions[0].kinetics
            reactions[0].kinetics = EEDFChannel(process=kin.process, collision_set=kin.collision_set)
    error = NonEquilibriumReverseRateError if problem == 'reversible' else PlasmaStateError
    with pytest.raises(error, match='cascade|EEDF|irreversible'):
        cascade_model(tmp_path, branch_fixture, mutate=mutate)


def test_unweighted_duplicate_preserves_error(tmp_path):
    def duplicate(core):
        from rmgpy.data.kinetics.library import LibraryReaction
        return LibraryReaction(reactants=[core[0], core[1]], products=[core[2], core[0], core[0]],
                               reversible=False, kinetics=EEDFChannel(
                                   process='e + Ar -> e + e + Ar+', collision_set='fixture'),
                               library='Toy', entry=Entry(index=1))
    with pytest.raises(PlasmaStateError, match='more than one EEDF reaction maps table channel 0 side.*double-count'):
        build_reactor(tmp_path, composition=False, extra_reaction=duplicate)


def test_model_enlarger_refuses_split_without_mutating_membership(branch_fixture, monkeypatch):
    from rmgpy.rmg.model import CoreEdgeReactionModel
    import rmgpy.rmg.main as main
    _, reactions, _ = branch_fixture
    model = CoreEdgeReactionModel()
    electron, ground = reactions[0].reactants
    products = [reaction.products[1] for reaction in reactions]
    model.core.species = [electron, ground]
    model.edge.species = products[:]
    model.edge.reactions = reactions[:]
    monkeypatch.setattr(main, 'determine_procnum_from_ram', lambda: 1)
    before = (model.core.species[:], model.edge.species[:], model.core.reactions[:], model.edge.reactions[:])
    with pytest.raises(PlasmaStateError, match='split cascade family promotion'):
        model.enlarge(products[0])
    assert (model.core.species, model.edge.species, model.core.reactions, model.edge.reactions) == before
    with pytest.raises(PlasmaStateError, match='split cascade family promotion'):
        model.add_reaction_to_core(reactions[0])
    assert (model.core.species, model.edge.species, model.core.reactions, model.edge.reactions) == before


def test_model_promotes_whole_family_and_allows_initial_seed(branch_fixture, monkeypatch):
    from rmgpy.rmg.model import CoreEdgeReactionModel
    import rmgpy.rmg.model as model_module
    _, reactions, _ = branch_fixture
    model = CoreEdgeReactionModel()
    electron, ground = reactions[0].reactants
    products = [reaction.products[1] for reaction in reactions]
    model.core.species = [electron] + products
    model.edge.species = [ground]
    model.edge.reactions = reactions[:]
    monkeypatch.setattr(model_module, 'get_db', lambda name: SimpleNamespace(is_molecule_forbidden=lambda molecule: False))
    assert model.add_species_to_core(ground) == reactions
    assert model.edge.reactions == []
    assert model.core.reactions == reactions
    seed = CoreEdgeReactionModel()
    for reaction in reactions:
        seed.add_reaction_to_core(reaction)
    assert seed.core.reactions == reactions
