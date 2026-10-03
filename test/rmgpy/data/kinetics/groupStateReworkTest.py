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

"""State constraints on every matching and recipe route."""
import pytest

from rmgpy.data.base import Entry, ForbiddenStructures, LogicOr, LogicAnd
from rmgpy.data.kinetics.family import KineticsFamily, ReactionRecipe
from rmgpy.data.kinetics.groups import KineticsGroups
from rmgpy.data.thermo import ThermoGroups
from rmgpy.data.transport import TransportGroups
from rmgpy.data.solvation import SoluteGroups
from rmgpy.exceptions import InvalidActionError, InvalidAdjacencyListError, UndeterminableKineticsError
from rmgpy.molecule import Molecule
from rmgpy.molecule.group import Group
from rmgpy.molecule.fragment import Fragment
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.data.kinetics.database import KineticsDatabase


N2 = '1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}'

GROUP = '1 *1 N u0 p1 c0 {2,T}\n2 *2 N u0 p1 c0 {1,T}'

LOWER_BOND = [['CHANGE_BOND', '*1', -1, '*2'], ['GAIN_RADICAL', '*1', 1], ['GAIN_RADICAL', '*2', 1]]

def write_family(path, label, declaration='', headers='', actions=None, group=GROUP, products=None):
    family_path = path / label
    family_path.mkdir()
    source = "name = 'State test'\nshortDesc = ''\nlongDesc = ''\n"
    source += "template(reactants=['Reactant'], products={}, ownReverse=False)\n".format(products or ['Product'])
    source += 'reversible = False\n' + declaration + '\n'
    source += 'recipe(actions={!r})\n'.format(LOWER_BOND if actions is None else actions)
    source += "entry(index=1, label='Reactant', kinetics=None, group={!r})\n".format(headers + '\n' + group)
    source += "tree('L1: Reactant')\n"
    (family_path / 'groups.py').write_text(source)
    (family_path / 'rules.py').write_text("name = 'State test rules'\nshortDesc = ''\nlongDesc = ''\n")
    return family_path

def load_family(path, label):
    database = KineticsDatabase()
    database.load_families(str(path), families=[label], depositories=[])
    return database.families[label]


@pytest.mark.parametrize('tree_type', [KineticsGroups, ThermoGroups, TransportGroups, SoluteGroups])
@pytest.mark.parametrize('state', [('A', -1), ('', 1)])
def test_others_fallback_respects_state(tree_type, state):
    tree = tree_type()
    root = Entry(label='Root', item=Group().from_adjacency_list(
        'electronicstate x\nvibrationallevel x\n1 *1 R ux'))
    fallback = Entry(label='Others-ground', item=Group().from_adjacency_list('1 *1 Ar u0'), parent=root)
    root.children = [fallback]
    tree.entries = {entry.label: entry for entry in [root, fallback]}
    tree.top = [root]
    molecule = Molecule().from_adjacency_list('1 *1 Ar u0 p4 c0')
    molecule.electronic_state, molecule.vibrational_level = state
    assert tree.descend_tree(molecule, molecule.get_all_labeled_atoms()) is root
    # A catch-all retains its structural fallback role for unresolved molecules.
    molecule.electronic_state, molecule.vibrational_level = '', -1
    assert tree.descend_tree(molecule, molecule.get_all_labeled_atoms()) is fallback
    molecule.electronic_state, molecule.vibrational_level = state
    fallback.item.electronic_state = [state[0]] if state[0] else []
    fallback.item.vibrational_level = [state[1]] if state[1] >= 0 else []
    assert tree.descend_tree(molecule, molecule.get_all_labeled_atoms()) is fallback


@pytest.mark.parametrize('state', [('A', -1), ('', 1)])
def test_reverse_preserves_unchanged_state_and_generates(tmp_path, state):
    headers = 'electronicstate [A]' if state[0] else 'vibrationallevel [1]'
    path = write_family(tmp_path, 'Reversible', 'allowExcitedReactants = True\nreversible = True', headers)
    family = load_family(tmp_path, 'Reversible')
    reactant = Molecule().from_adjacency_list(GROUP)
    reactant.electronic_state, reactant.vibrational_level = state
    products = family.apply_recipe([reactant])
    assert (products[0].electronic_state, products[0].vibrational_level) == state
    rebuilt = family.apply_recipe(products, forward=False)
    assert rebuilt[0].is_isomorphic(reactant)
    assert family.generate_reactions(products, prod_resonance=False)


@pytest.mark.parametrize('state', [('', -1), ('A', 1)])
def test_set_state_reversal_restores_declared_input(tmp_path, state):
    headers = 'electronicstate [A]\nvibrationallevel [1]' if state[0] else ''
    write_family(tmp_path, 'ExciteReverse', 'allowExcitedReactants = True\nreversible = True', headers,
                 actions=[['SET_STATE', '*1', 'B', 2]])
    family = load_family(tmp_path, 'ExciteReverse')
    molecule = Molecule().from_adjacency_list(GROUP)
    molecule.electronic_state, molecule.vibrational_level = state
    products = family.apply_recipe([molecule])
    assert products[0].electronic_state == 'B'
    assert family.apply_recipe(products, forward=False)[0].is_isomorphic(molecule)


def test_ambiguous_state_reversal_has_named_family_error(tmp_path):
    write_family(tmp_path, 'Ambiguous', 'allowExcitedReactants = True\nreversible = True', 'electronicstate x',
                 actions=[['SET_STATE', '*1', 'B', -1]])
    with pytest.raises(InvalidActionError) as error:
        load_family(tmp_path, 'Ambiguous')
    assert type(error.value).__name__ == 'StateReversalError'
    assert 'Ambiguous' in str(error.value)


@pytest.mark.parametrize('constraint', ['electronicstate x', 'electronicstate [A]\nvibrationallevel [1]'])
def test_split_template_carries_constraints_and_generates(tmp_path, constraint):
    group = '1 *1 C u1\n2 *2 C u1'
    write_family(tmp_path, 'Recombine', 'allowExcitedReactants = True', constraint,
                 actions=[['FORM_BOND', '*1', 1, '*2'], ['LOSE_RADICAL', '*1', 1], ['LOSE_RADICAL', '*2', 1]]
                 + ([['SET_STATE', '*1', 'A', 1]] if '[1]' in constraint else []),
                 group=group)
    family = load_family(tmp_path, 'Recombine')
    original = family.forward_template.reactants[0].item
    for part in original.split():
        assert part.electronic_state == original.electronic_state
        assert part.vibrational_level == original.vibrational_level
        assert part.electronic_state is not original.electronic_state
    first = Molecule(smiles='[CH3]', electronic_state='A', vibrational_level=1 if '[1]' in constraint else -1)
    assert family.generate_reactions([first, first.copy(deep=True)], prod_resonance=False)


@pytest.mark.parametrize('states', [('A', 'A'), ('A', '')])
def test_merged_template_tree_checks_each_reactant_state(states):
    groups = KineticsGroups(label='Merge/groups')
    root = Entry(label='Root', item=Group().from_adjacency_list('electronicstate x\n1 *1 Ar u0\n2 *2 Ar u0'))
    ground = Entry(label='Others-ground', item=Group().from_adjacency_list('1 *1 Ar u0\n2 *2 Ar u0'), parent=root)
    root.children = [ground]
    groups.top = [root]
    groups.entries = {entry.label: entry for entry in [root, ground]}
    molecules = [Molecule().from_adjacency_list('1 *{} Ar u0 p4 c0'.format(i+1)) for i in range(2)]
    for molecule, state in zip(molecules, states):
        molecule.electronic_state = state
    assert groups.get_reaction_template(Reaction(reactants=molecules)) == [root]


@pytest.mark.parametrize('states', [('A', 'A'), ('A', '')])
def test_reaction_matches_admits_and_checks_each_state(states):
    family = KineticsFamily(allow_excited_reactants=True)
    molecules = [Molecule().from_adjacency_list('1 *{} Ar u0 p4 c0'.format(i+1)) for i in range(2)]
    for molecule, state in zip(molecules, states):
        molecule.electronic_state = state
    reaction = Reaction(reactants=[Species(molecule=[molecule]) for molecule in molecules])
    group = Group().from_adjacency_list('electronicstate x\n1 *1 Ar u0\n2 *2 Ar u0')
    assert family.reaction_matches(reaction, group)
    group.electronic_state = ['A']
    assert family.reaction_matches(reaction, group) == (states == ('A', 'A'))
    group.electronic_state = []
    assert not family.reaction_matches(reaction, group)
    group.electronic_state = ['x']
    family.allow_excited_reactants = False
    assert not family.reaction_matches(reaction, group)


@pytest.mark.parametrize('multiplicity', [[], [1]])
def test_product_state_is_validated_without_reactant_gate(tmp_path, multiplicity):
    write_family(tmp_path, 'ExcitationProduct', '', actions=[['SET_STATE', '*1', 'A', -1]])
    family = load_family(tmp_path, 'ExcitationProduct')
    family.forward_template.products[0].item.multiplicity = multiplicity
    molecule = Molecule().from_adjacency_list(GROUP)
    assert family.generate_reactions([molecule], prod_resonance=False)
    family.forward_recipe = ReactionRecipe([])
    assert family.apply_recipe([molecule]) is None
    # Logical product templates obey their constituent Group state constraints too.
    product = family.forward_template.products[0]
    family.groups.entries['A-only'] = Entry(label='A-only', item=product.item)
    product.item = LogicOr(['A-only'], False)
    assert family.apply_recipe([molecule]) is None
    family.forward_recipe = ReactionRecipe([['SET_STATE', '*1', 'A', -1]])
    assert family.apply_recipe([Molecule().from_adjacency_list(GROUP)])


@pytest.mark.parametrize('keyword', ['electronicstate [A]', 'vibrationallevel [1]'])
def test_fragment_boolean_and_forbidden_respect_state(keyword):
    fragment = Fragment().from_smiles('C')
    group = Group().from_adjacency_list(keyword + '\n1 C u0')
    assert fragment.find_subgraph_isomorphisms(group) == []
    assert not fragment.is_subgraph_isomorphic(group)
    forbidden = ForbiddenStructures()
    forbidden.entries = {'Excited': Entry(label='Excited', item=group)}
    assert not forbidden.is_molecule_forbidden(fragment)


@pytest.mark.parametrize('route', ['constructor', 'setter', 'adjacency', 'recipe'])
def test_x_is_reserved_for_group_wildcards(route):
    with pytest.raises((ValueError, InvalidActionError, InvalidAdjacencyListError)) as error:
        if route == 'constructor':
            Molecule(electronic_state='x')
        elif route == 'setter':
            Molecule().electronic_state = 'x'
        elif route == 'adjacency':
            Molecule().from_adjacency_list('electronicstate x\n1 Ar u0 p4 c0')
        else:
            molecule = Molecule().from_adjacency_list('1 *1 Ar u0 p4 c0')
            ReactionRecipe([['SET_STATE', '*1', 'x', -1]]).apply_forward(molecule)
    assert 'reserved' in str(error.value)


@pytest.mark.parametrize('header', ['electronicstate [A]', 'vibrationallevel [1]', 'electronicstate x\nvibrationallevel x'])
def test_state_only_groups_round_trip(header):
    group = Group().from_adjacency_list(header)
    restored = Group().from_adjacency_list(group.to_adjacency_list())
    assert restored.electronic_state == group.electronic_state
    assert restored.vibrational_level == group.vibrational_level


@pytest.mark.parametrize('kwargs', [
    {'electronic_state': 'AB'}, {'electronic_state': [1]}, {'electronic_state': ['A/B']},
    {'vibrational_level': [True]}, {'vibrational_level': [1.0]}, {'vibrational_level': [-2]},
    {'vibrational_level': [2147483648]},
])
def test_group_constructor_rejects_invalid_constraints(kwargs):
    with pytest.raises(ValueError):
        Group(**kwargs)


def test_family_template_selection_enforces_opt_in():
    molecule = Molecule().from_adjacency_list('electronicstate A\n1 *1 Ar u0 p4 c0')
    root = Entry(label='Root', item=Group().from_adjacency_list('electronicstate x\n1 *1 Ar u0'))
    family = KineticsFamily(label='Unopted')
    family.groups = KineticsGroups(top=[root], label='Unopted/groups')
    with pytest.raises(UndeterminableKineticsError):
        family.get_reaction_template(Reaction(reactants=[molecule]))


def test_split_reactions_uses_state_aware_matching():
    molecule = Molecule().from_adjacency_list('electronicstate A\n1 *1 Ar u0 p4 c0')
    second = molecule.copy(deep=True)
    second.atoms[0].label = '*2'
    reaction = Reaction(reactants=[Species(molecule=[molecule]), Species(molecule=[second])])
    group = Group().from_adjacency_list('electronicstate [A]\n1 *1 Ar u0\n2 *2 Ar u0')
    family = KineticsFamily(allow_excited_reactants=True)
    assert family._split_reactions([reaction], group) == ([reaction], [], [0])
    group.electronic_state = []
    assert family._split_reactions([reaction], group) == ([], [reaction], [])


def test_ambiguous_species_boundary_reverse_is_refused_at_load(tmp_path):
    write_family(tmp_path, 'AmbiguousRecombine', 'allowExcitedReactants = True\nreversible = True',
                 'electronicstate x', group='1 *1 C u1\n2 *2 C u1',
                 actions=[['FORM_BOND', '*1', 1, '*2'], ['LOSE_RADICAL', '*1', 1], ['LOSE_RADICAL', '*2', 1]])
    with pytest.raises(InvalidActionError) as error:
        load_family(tmp_path, 'AmbiguousRecombine')
    assert type(error.value).__name__ == 'StateReversalError'
    assert 'AmbiguousRecombine' in str(error.value)


@pytest.mark.parametrize('tree_type', [KineticsGroups, ThermoGroups, TransportGroups, SoluteGroups])
@pytest.mark.parametrize('logic_type', [LogicOr, LogicAnd])
@pytest.mark.parametrize('state', [('A', -1), ('', 1)])
def test_negated_tree_nodes_keep_declared_state_domain(tree_type, logic_type, state):
    tree = tree_type()
    root = Entry(label='Root', item=Group().from_adjacency_list(
        'electronicstate x\nvibrationallevel x\n1 *1 R ux'))
    excluded = [Entry(label=element, item=Group().from_adjacency_list('1 *1 {} u0'.format(element)))
                for element in ['C', 'N']]
    negated = Entry(label='Negated-ground', item=logic_type([entry.label for entry in excluded], True), parent=root)
    root.children = [negated]
    tree.entries = {entry.label: entry for entry in [root, negated] + excluded}
    tree.top = [root]
    molecule = Molecule().from_adjacency_list('1 *1 Ar u0 p4 c0')
    molecule.electronic_state, molecule.vibrational_level = state
    assert tree.descend_tree(molecule, molecule.get_all_labeled_atoms()) is root
    molecule.electronic_state, molecule.vibrational_level = '', -1
    assert tree.descend_tree(molecule, molecule.get_all_labeled_atoms()) is negated
    for entry in excluded:
        entry.item.electronic_state = [state[0]] if state[0] else []
        entry.item.vibrational_level = [state[1]] if state[1] >= 0 else []
    molecule.electronic_state, molecule.vibrational_level = state
    assert tree.descend_tree(molecule, molecule.get_all_labeled_atoms()) is negated


@pytest.mark.parametrize('logic_type,invert', [(LogicOr, True), (LogicAnd, False)])
@pytest.mark.parametrize('state', [('A', -1), ('', 1)])
def test_empty_logical_patterns_keep_unresolved_default(logic_type, invert, state):
    tree = KineticsGroups()
    logic = logic_type([], invert)
    molecule = Molecule(smiles='C')
    molecule.electronic_state, molecule.vibrational_level = state
    assert not logic.match_to_structure(tree, molecule, {})
    molecule.electronic_state, molecule.vibrational_level = '', -1
    assert logic.match_to_structure(tree, molecule, {})


def write_boundary_family(path, label, right_header):
    family_path = path / label
    family_path.mkdir()
    actions = [['FORM_BOND', '*1', 1, '*2'], ['LOSE_RADICAL', '*1', 1],
               ['LOSE_RADICAL', '*2', 1], ['SET_STATE', '*1', 'B', -1]]
    source = "name = 'Boundary state test'\nshortDesc = ''\nlongDesc = ''\n"
    source += "template(reactants=['Left', 'Right'], products=['Product'], ownReverse=False)\n"
    source += 'reversible = True\nallowExcitedReactants = True\nrecipe(actions={!r})\n'.format(actions)
    for index, (name, header) in enumerate([
            ('Left', 'electronicstate [A]\nvibrationallevel [1]'), ('Right', right_header)], 1):
        adjacency = header + '\n1 *{} C u1 p0 c0'.format(index)
        source += "entry(index={}, label={!r}, kinetics=None, group={!r})\n".format(index, name, adjacency)
    source += "tree('L1: Left\\nL1: Right')\n"
    (family_path / 'groups.py').write_text(source)
    (family_path / 'rules.py').write_text("name = 'Boundary rules'\nshortDesc = ''\nlongDesc = ''\n")


@pytest.mark.parametrize('right_header,right_state', [
    ('electronicstate [A]', ('A', -1)), ('', ('', -1)), ('vibrationallevel [2]', ('', 2))])
def test_set_state_inverse_restores_every_original_species(tmp_path, right_header, right_state):
    write_boundary_family(tmp_path, 'BoundaryInverse', right_header)
    family = load_family(tmp_path, 'BoundaryInverse')
    states = [('A', 1), right_state]
    reactants = []
    for index, state in enumerate(states, 1):
        molecule = Molecule(smiles='[CH3]', electronic_state=state[0], vibrational_level=state[1])
        next(atom for atom in molecule.atoms if atom.is_carbon()).label = '*{}'.format(index)
        reactants.append(molecule)
    products = family.apply_recipe(reactants, relabel_atoms=False)
    assert [(m.electronic_state, m.vibrational_level) for m in products] == [('B', -1)]
    restored = family.apply_recipe(products, forward=False, relabel_atoms=False)
    assert [(m.electronic_state, m.vibrational_level) for m in restored] == states
    assert all(m.is_isomorphic(original) for m, original in zip(restored, reactants))


def test_set_state_boundary_inverse_refuses_unknown_other_species(tmp_path):
    write_boundary_family(tmp_path, 'BoundaryAmbiguous', 'electronicstate x')
    with pytest.raises(InvalidActionError) as error:
        load_family(tmp_path, 'BoundaryAmbiguous')
    assert type(error.value).__name__ == 'StateReversalError'
    assert 'BoundaryAmbiguous' in str(error.value)


def test_distinct_state_product_templates_need_distinct_matching_products(tmp_path):
    write_family(tmp_path, 'ProductAssignment', group='1 *1 C u1 p0 c0\n2 *2 C u1 p0 c0',
                 actions=[['SET_STATE', '*1', 'A', -1]], products=['First', 'Second'])
    family = load_family(tmp_path, 'ProductAssignment')
    for entry in family.forward_template.products:
        entry.item.electronic_state = ['A']
    reactants = []
    for index in [1, 2]:
        molecule = Molecule(smiles='[CH3]')
        next(atom for atom in molecule.atoms if atom.is_carbon()).label = '*{}'.format(index)
        reactants.append(molecule)
    assert family.apply_recipe(reactants) is None
