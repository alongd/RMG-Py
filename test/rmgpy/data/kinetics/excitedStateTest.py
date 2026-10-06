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

"""Family admission and product-state recipe behavior."""
import itertools
from pathlib import Path

import pytest

from rmgpy import settings
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.molecule import Molecule
from rmgpy.species import Species

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


@pytest.mark.parametrize('opt_in', [None, False, True])
@pytest.mark.parametrize('state', [('A3Su+', -1), ('', 1)])
@pytest.mark.parametrize('constraint', ['wildcard', 'explicit'])
def test_family_requires_opt_in_even_when_template_admits_state(tmp_path, opt_in, state, constraint):
    declaration = '' if opt_in is None else 'allowExcitedReactants = {}'.format(opt_in)
    if constraint == 'wildcard':
        headers = 'electronicstate x\nvibrationallevel x'
    else:
        electronic = '[A3Su+]' if state[0] else '[]'
        vibrational = '[1]' if state[1] >= 0 else '[]'
        headers = 'electronicstate {}\nvibrationallevel {}'.format(electronic, vibrational)
    write_family(tmp_path, 'StateTest', declaration, headers)
    family = load_family(tmp_path, 'StateTest')
    reactant = Molecule().from_adjacency_list(N2)
    reactant.electronic_state, reactant.vibrational_level = state
    reactions = family.generate_reactions([reactant], prod_resonance=False)
    assert bool(reactions) == (opt_in is True)
    assert (reactant.electronic_state, reactant.vibrational_level) == state


def test_family_opt_in_is_saved_and_does_not_leak_between_loads(tmp_path):
    headers = 'electronicstate x\nvibrationallevel x'
    write_family(tmp_path, 'AOptedIn', 'allowExcitedReactants = True', headers)
    write_family(tmp_path, 'BDefault', '', headers)
    database = KineticsDatabase()
    database.load_families(str(tmp_path), families=['AOptedIn', 'BDefault'], depositories=[])
    first, second = database.families['AOptedIn'], database.families['BDefault']
    assert first.allow_excited_reactants is True
    assert second.allow_excited_reactants is False
    saved = tmp_path / 'saved.py'
    first.save_groups(str(saved))
    assert 'allowExcitedReactants = True' in saved.read_text()


@pytest.mark.parametrize('state', [('A3Su+', -1), ('', 0), ('', 1), ('A3Su+', 1)])
def test_set_state_family_emits_resolved_product_and_template(tmp_path, state):
    write_family(tmp_path, 'Excitation', 'allowExcitedReactants = True',
                 actions=[['SET_STATE', '*1', state[0], state[1]]])
    family = load_family(tmp_path, 'Excitation')
    group = family.forward_template.products[0].item
    assert group.electronic_state == ([state[0]] if state[0] else [])
    assert group.vibrational_level == ([state[1]] if state[1] >= 0 else [])
    reactant = Molecule().from_adjacency_list(N2)
    reactions = family.generate_reactions([reactant], prod_resonance=False)
    assert reactions
    for reaction in reactions:
        assert (reaction.products[0].electronic_state, reaction.products[0].vibrational_level) == state
    assert not reactant.has_resolved_state()


def test_set_state_assigns_only_the_final_product_containing_its_atom(tmp_path):
    reactant = Molecule(smiles='NN')
    nitrogens = [atom for atom in reactant.atoms if atom.symbol == 'N']
    nitrogens[0].label, nitrogens[1].label = '*1', '*2'
    actions = [['SET_STATE', '*1', 'test', 2], ['BREAK_BOND', '*1', 1, '*2'],
               ['GAIN_RADICAL', '*1', 1], ['GAIN_RADICAL', '*2', 1]]
    write_family(tmp_path, 'SplitState', 'allowExcitedReactants = True', actions=actions,
                 group=reactant.to_group().to_adjacency_list(), products=['First', 'Second'])
    family = load_family(tmp_path, 'SplitState')
    products = family.apply_recipe([reactant])
    assert len(products) == 2
    first = next(product for product in products if product.contains_labeled_atom('*1'))
    second = next(product for product in products if product.contains_labeled_atom('*2'))
    assert (first.electronic_state, first.vibrational_level) == ('test', 2)
    assert not second.has_resolved_state()
    assert not reactant.has_resolved_state()
    assert reactant.has_bond(*nitrogens)


def test_set_state_direct_recipe_and_reverse_refusal():
    from rmgpy.data.kinetics.family import ReactionRecipe
    from rmgpy.exceptions import InvalidActionError
    reactant = Molecule().from_adjacency_list('1 *1 Ar u0 p4 c0')
    recipe = ReactionRecipe([['SET_STATE', '*1', '1s5', -1]])
    recipe.apply_forward(reactant)
    assert reactant.electronic_state == '1s5'
    with pytest.raises(InvalidActionError, match='SET_STATE'):
        recipe.get_reverse()


@pytest.mark.database
def test_all_database_families_exclude_resolved_nitrogen_and_argon():
    family_root = Path(settings['database.directory']) / 'kinetics/families'
    database = KineticsDatabase()
    database.load_families(str(family_root), families='all', depositories=[])
    expected_families = {path.parent.name for path in family_root.glob('*/groups.py')}
    assert set(database.families) == expected_families
    assert not any(family.allow_excited_reactants for family in database.families.values())
    adjacency = {'N2': N2, 'N2v1': 'vibrationallevel 1\n' + N2,
                 'N2A': 'electronicstate A3Su+\n' + N2, 'Ar': '1 Ar u0 p4 c0',
                 'Ar1s5': 'electronicstate 1s5\n1 Ar u0 p4 c0', 'e-': '1 e u1 p0 c-1'}
    molecules = {name: Molecule().from_adjacency_list(adj) for name, adj in adjacency.items()}
    cases = [(name,) for name in sorted(molecules)]
    cases += list(itertools.combinations_with_replacement(sorted(molecules), 2))
    resolved_names = {'N2v1', 'N2A', 'Ar1s5'}
    negative_queries, ground_reactions = 0, 0
    for label, family in database.families.items():
        for names in cases:
            reactants = [molecules[name].copy(deep=True) for name in names]
            reactions = family.generate_reactions(reactants, prod_resonance=False)
            if any(name in resolved_names for name in names):
                negative_queries += 1
                assert reactions == [], (label, names, reactions)
            else:
                ground_reactions += len(reactions)
    print('CENSUS: {} families; {} cases/family (6 unary, 21 pairwise); '
          '{} resolved queries; 0 resolved reactions; {} ground reactions'.format(
              len(database.families), len(cases), negative_queries, ground_reactions))


@pytest.mark.database
def test_pairing_argon_reaction_resolves_electron_placement():
    family_root = Path(settings['database.directory']) / 'kinetics/families'
    database = KineticsDatabase()
    database.load_families(str(family_root), families=['Plasma_Radiative_Recombination_Pairing'],
                           depositories=['training'])
    argon = Molecule().from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    family = database.families['Plasma_Radiative_Recombination_Pairing']
    reaction = family.generate_reactions([argon])[0]
    reaction.kinetics = family.get_kinetics(
        reaction, template_labels=reaction.template, degeneracy=reaction.degeneracy)[0][0]
    reaction.ensure_species()
    electron = Species(label='e').from_adjacency_list('1 e u1 p0 c-1')
    from rmgpy.electron_placement import resolve_electron_placement
    view = resolve_electron_placement(reaction, [electron] + reaction.reactants + reaction.products)
    assert view.electrons == 0
    assert sum(species.is_electron() for species in view.reactants) == 1


@pytest.mark.parametrize('state', [('A/B', -1), ('A', True), ('A', -2), ('A', 2147483648)])
def test_set_state_rejects_invalid_state_values(state):
    from rmgpy.data.kinetics.family import ReactionRecipe
    from rmgpy.exceptions import InvalidActionError
    molecule = Molecule().from_adjacency_list('1 *1 Ar u0 p4 c0')
    with pytest.raises(InvalidActionError, match='SET_STATE'):
        ReactionRecipe([['SET_STATE', '*1', state[0], state[1]]]).apply_forward(molecule)
    assert not molecule.has_resolved_state()


def test_set_state_can_clear_resolved_fields_and_later_assignment_wins(tmp_path):
    headers = 'electronicstate x\nvibrationallevel x'
    actions = [['SET_STATE', '*1', 'temporary', 2], ['SET_STATE', '*1', '', -1]]
    write_family(tmp_path, 'ClearState', 'allowExcitedReactants = True', headers, actions=actions)
    family = load_family(tmp_path, 'ClearState')
    reactant = Molecule().from_adjacency_list('electronicstate A3Su+\nvibrationallevel 1\n' + N2)
    reactions = family.generate_reactions([reactant], prod_resonance=False)
    assert reactions
    assert all(not reaction.products[0].has_resolved_state() for reaction in reactions)
    assert (reactant.electronic_state, reactant.vibrational_level) == ('A3Su+', 1)
