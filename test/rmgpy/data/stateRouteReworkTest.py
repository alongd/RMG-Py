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

"""Resolved fallback templates, electron libraries and frequency estimation."""

import pytest

from rmgpy.data.base import Entry
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import KineticsFamily
from rmgpy.data.kinetics.groups import KineticsGroups
from rmgpy.data.kinetics.rules import KineticsRules
from rmgpy.data.statmech import StatmechDatabase, StatmechGroups, StatmechLibrary, GroupFrequencies
from rmgpy.data.thermo import ThermoDatabase, ThermoLibrary
from rmgpy.exceptions import ExcitedSpeciesThermoError, StateProvenanceError
from rmgpy.kinetics import Arrhenius, ArrheniusEP
from rmgpy.molecule import Group, Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.statmech import Conformer
from rmgpy.thermo import ThermoData

STATES = [('A', -1), ('', 1), ('A', 1)]
WILDCARD = 'electronicstate x\nvibrationallevel x\n'


def thermo():
    return ThermoData(Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
                      Cpdata=([75, 90, 100, 120, 140, 150, 170], 'J/(mol*K)'),
                      H298=(12345, 'J/mol'), S298=(12, 'J/(mol*K)'))


def fallback_family():
    root = Entry(label='Any', item=Group().from_adjacency_list(WILDCARD + '1 * R u0'))
    child = Entry(label='Others-O', item=Group().from_adjacency_list(WILDCARD + '1 * O u0'),
                  parent=root, nodal_distance=1)
    root.children = [child]
    family = KineticsFamily(label='Probe', allow_excited_reactants=True)
    family.groups = KineticsGroups(label='Probe/groups', entries={e.label: e for e in (root, child)}, top=[root])
    family.rules = KineticsRules(label='Probe/rules')
    family.rules.entries = {child.label: [Entry(label=child.label, item=[child], rank=1,
        data=ArrheniusEP(A=(123, 's^-1'), n=0, alpha=0, E0=(0, 'J/mol')))]}
    return family, root, child


@pytest.mark.parametrize('state', STATES)
@pytest.mark.parametrize('strict', [False, True])
def test_resolved_fallback_matches_original_molecule(state, strict):
    family, root, child = fallback_family()
    molecule = Molecule(smiles='C', electronic_state=state[0], vibrational_level=state[1])
    molecule.atoms[0].label = '*'
    atoms = {'*': molecule.atoms[0]}
    assert family.groups.match_node_to_structure(root, molecule, atoms, strict)
    assert not family.groups.match_node_to_structure(child, molecule, atoms, strict)
    with pytest.raises(StateProvenanceError, match='Others-O.*does not match'):
        family.groups.descend_tree(molecule, atoms, strict=strict)
    with pytest.raises(StateProvenanceError):
        template = family.get_reaction_template(Reaction(reactants=[molecule], products=[]))
        family.get_kinetics_for_template(template)


def test_unresolved_fallback_unchanged():
    family, root, child = fallback_family()
    molecule = Molecule(smiles='C')
    molecule.atoms[0].label = '*'
    template = family.get_reaction_template(Reaction(reactants=[molecule], products=[]))
    assert template == [child]
    assert family.get_kinetics_for_template(template)[0].A.value_si == 123


@pytest.mark.parametrize('state', STATES)
def test_resolved_matching_child_allowed(state):
    family, root, child = fallback_family()
    molecule = Molecule(smiles='O', electronic_state=state[0], vibrational_level=state[1])
    molecule.atoms[0].label = '*'
    template = family.get_reaction_template(Reaction(reactants=[molecule], products=[]))
    assert template == [child]
    assert family.get_kinetics_for_template(template)[0].A.value_si == 123


def electron_database(state=('', -1)):
    source = Molecule().from_adjacency_list('1 e u0 p0 c-1')
    source.electronic_state, source.vibrational_level = state
    library = ThermoLibrary(label='Electron')
    library.entries = {'Electron': Entry(label='Electron', item=source, data=thermo())}
    database = ThermoDatabase()
    database.libraries = {'Electron': library}
    database.library_order = ['Electron']
    return database, library


def electron(state, radicals):
    molecule = Molecule().from_adjacency_list('1 e u{} p0 c-1'.format(radicals))
    species = Species(molecule=[molecule])
    species.molecule[0].electronic_state, species.molecule[0].vibrational_level = state
    return species


@pytest.mark.parametrize('state', STATES)
@pytest.mark.parametrize('radicals', [0, 1])
@pytest.mark.parametrize('api', ['direct', 'full'])
def test_resolved_electron_refuses_unmatched_library(state, radicals, api):
    database, library = electron_database()
    species = electron(state, radicals)
    with pytest.raises(ExcitedSpeciesThermoError, match='electron'):
        if api == 'direct':
            database.get_thermo_data_from_library(species, library)
        else:
            database.get_thermo_data(species)


@pytest.mark.parametrize('state', [('', -1)] + STATES)
@pytest.mark.parametrize('radicals', [0, 1])
def test_electron_matching_state_keeps_multiplicity_convention(state, radicals):
    database, library = electron_database(state)
    species = electron(state, radicals)
    if state != ('', -1):
        with pytest.raises(ExcitedSpeciesThermoError, match='electron'):
            database.get_thermo_data_from_library(species, library)
        with pytest.raises(ExcitedSpeciesThermoError, match='electron'):
            database.get_thermo_data(species)
    else:
        assert database.get_thermo_data_from_library(species, library)[0].H298.value_si == 12345
        assert database.get_thermo_data(species).H298.value_si == 12345


@pytest.mark.parametrize('state', STATES)
@pytest.mark.parametrize('smiles', ['C', 'C1CC1'])
@pytest.mark.parametrize('api', ['node', 'frequency_groups', 'group_fit', 'database_groups', 'full'])
def test_resolved_statmech_estimation_refused(state, smiles, api):
    molecule = Molecule(smiles=smiles, electronic_state=state[0], vibrational_level=state[1])
    groups = StatmechGroups(label='groups')
    root = Entry(label='Any', item=Group().from_adjacency_list(WILDCARD + '1 * R u0'))
    child = Entry(label='Others-O', item=Group().from_adjacency_list(WILDCARD + '1 * O u0'),
                  parent=root, data=GroupFrequencies([(1234., 1234., 1)]))
    root.children = [child]
    groups.top = [root]
    groups.entries = {e.label: e for e in (root, child)}
    database = StatmechDatabase()
    database.groups = {'groups': groups}
    expected = ExcitedSpeciesThermoError if api in ('database_groups', 'full') else StateProvenanceError
    message = 'Library-only' if expected is ExcitedSpeciesThermoError else 'statmech.*resolved'
    with pytest.raises(expected, match=message):
        if api == 'node':
            groups._get_node(molecule, {'*': molecule.atoms[0]})
        elif api == 'frequency_groups':
            groups.get_frequency_groups(molecule)
        elif api == 'group_fit':
            groups.get_statmech_data(molecule, thermo())
        elif api == 'database_groups':
            database.get_statmech_data_from_groups(molecule, thermo())
        else:
            database.get_statmech_data(molecule, thermo())


def test_unresolved_ring_frequencies_unchanged():
    groups = StatmechGroups()
    selected = groups.get_frequency_groups(Molecule(smiles='C1CC1'))
    assert len(selected) == 1
    entry, count = next(iter(selected.items()))
    assert (entry.label, entry.item, count) == ('ringCH', None, 6)
    assert entry.data.generate_frequencies(count) == [
        2750., 2830., 2910., 2990., 3070., 3150., 900., 940., 980., 1020., 1060., 1100.]


@pytest.mark.parametrize('state', STATES)
def test_statmech_library_requires_matching_state(state):
    molecule = Molecule(smiles='C1CC1', electronic_state=state[0], vibrational_level=state[1])
    entry = Entry(label='ring', item=Molecule(smiles='C1CC1'), data=Conformer())
    library = StatmechLibrary()
    library.entries = {entry.label: entry}
    database = StatmechDatabase()
    database.libraries = {'ring': library}
    database.library_order = ['ring']
    with pytest.raises(ExcitedSpeciesThermoError):
        database.get_statmech_data_from_library(molecule, library)
    entry.item = molecule.copy(deep=True)
    with pytest.raises(ExcitedSpeciesThermoError):
        database.get_statmech_data_from_library(molecule, library)
    with pytest.raises(ExcitedSpeciesThermoError):
        database.get_statmech_data(molecule, thermo())
    ground = Molecule(smiles='C1CC1')
    entry.item = ground.copy(deep=True)
    assert database.get_statmech_data_from_library(ground, library)[0] is entry.data
    assert database.get_statmech_data(ground, thermo()) is entry.data


@pytest.mark.parametrize('state', STATES)
@pytest.mark.parametrize('source_kind', ['Training', 'Rate Rules', 'Library', 'PDep'])
def test_resolved_source_reconstruction_refused(state, source_kind):
    molecule = Molecule(smiles='C', electronic_state=state[0], vibrational_level=state[1])
    species = Species(molecule=[molecule], thermo=thermo())
    reaction = Reaction(reactants=[species], products=[],
                        kinetics=Arrhenius(A=(123, 's^-1'), n=0, Ea=(0, 'J/mol')))
    ground = Entry(label='Ground', item=Reaction(reactants=[Species(smiles='C')], products=[]),
                   data=reaction.kinetics)
    rule = Entry(label='Ground', item=None, data=ArrheniusEP(A=(123, 's^-1'), n=0, alpha=0, E0=(0, 'J/mol')))
    sources = {'Training': ('probe', ground, False),
               'Rate Rules': ('probe', {'rules': [(rule, 1)], 'training': [], 'degeneracy': 1}),
               'Library': 'probe', 'PDep': 1}
    with pytest.raises(ExcitedSpeciesThermoError, match='Library-only'):
        KineticsDatabase().reconstruct_kinetics_from_source(reaction, {source_kind: sources[source_kind]})
