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

"""State constraints for reaction-family groups."""
import copy
import pickle

import pytest

from rmgpy.molecule import Molecule
from rmgpy.molecule.group import Group

N2 = '1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}'
GROUP = '1 N ux {2,T}\n2 N ux {1,T}'


@pytest.mark.parametrize('state,index', [(('', -1), 0), (('A3Su+', -1), 1),
                                        (('', 1), 2), (('A3Su+', 1), 3)])
@pytest.mark.parametrize('electronic,vibrational,matches', [
    ([], [], [True, False, False, False]),
    (['x'], ['x'], [True, True, True, True]),
    (['A3Su+'], [], [False, True, False, False]),
    ([], [1], [False, False, True, False]),
    (['A3Su+'], [1], [False, False, False, True]),
    (['B'], [2], [False, False, False, False]),
])
def test_molecule_matches_each_state_constraint_independently(state, index, electronic, vibrational, matches):
    molecule = Molecule().from_adjacency_list(N2)
    molecule.electronic_state, molecule.vibrational_level = state
    group = Group().from_adjacency_list(GROUP)
    group.electronic_state, group.vibrational_level = electronic, vibrational
    assert molecule.is_subgraph_isomorphic(group) == matches[index]
    assert bool(molecule.find_subgraph_isomorphisms(group)) == matches[index]
    assert molecule.is_subgraph_isomorphic(group, generate_initial_map=True) == matches[index]


@pytest.mark.parametrize('transport', ['copy', 'deep', 'python-copy', 'deepcopy', 'pickle'])
def test_group_state_lists_survive_transport_without_aliasing(transport):
    electronic, vibrational = ['A3Su+', 'B'], [0, 1]
    original = Group(atoms=Group().from_adjacency_list(GROUP).atoms,
                     electronic_state=electronic, vibrational_level=vibrational)
    electronic.append('C')
    vibrational.append(2)
    if transport == 'copy':
        restored = original.copy()
    elif transport == 'deep':
        restored = original.copy(deep=True)
    elif transport == 'python-copy':
        restored = copy.copy(original)
    elif transport == 'deepcopy':
        restored = copy.deepcopy(original)
    else:
        restored = pickle.loads(pickle.dumps(original))
    assert restored.electronic_state == ['A3Su+', 'B']
    assert restored.vibrational_level == [0, 1]
    restored.electronic_state.append('D')
    restored.vibrational_level.append(3)
    assert original.electronic_state == ['A3Su+', 'B']
    assert original.vibrational_level == [0, 1]


@pytest.mark.parametrize('headers,electronic,vibrational', [
    ('electronicstate [A3Su+,B3Pg]\nvibrationallevel [0,1]', ['A3Su+', 'B3Pg'], [0, 1]),
    ('electronicstate x\nvibrationallevel x', ['x'], ['x']),
    ('electronicstate []\nvibrationallevel []', [], []),
    ("electronicstate ['', 'a,b.c']\nvibrationallevel [-1,1]", ['', 'a,b.c'], [-1, 1]),
])
def test_group_adjacency_headers_round_trip_and_reset(headers, electronic, vibrational):
    original = Group().from_adjacency_list('label\nmultiplicity [1]\n' + headers + '\n' + GROUP)
    assert original.electronic_state == electronic
    assert original.vibrational_level == vibrational
    restored = Group().from_adjacency_list(original.to_adjacency_list(label='label'))
    assert restored.electronic_state == electronic
    assert restored.vibrational_level == vibrational
    restored.from_adjacency_list(GROUP)
    assert restored.electronic_state == []
    assert restored.vibrational_level == []


@pytest.mark.parametrize('headers', ['electronicstate A3Su+', 'vibrationallevel 1',
                                    'electronicstate [A/B]', 'electronicstate [A,]', 'vibrationallevel [-2]',
                                    'vibrationallevel [True]', 'vibrationallevel [2147483648]',
                                    'electronicstate [A]\nelectronicstate [B]'])
def test_malformed_group_state_headers_are_refused(headers):
    from rmgpy.exceptions import InvalidAdjacencyListError
    with pytest.raises(InvalidAdjacencyListError):
        Group().from_adjacency_list(headers + '\n' + GROUP)


@pytest.mark.parametrize('field,a,b,unresolved', [('electronic_state', 'A3Su+', 'B', ''),
                                                ('vibrational_level', 1, 2, -1)])
@pytest.mark.parametrize('left,right,equal,specific', [
    ('missing', 'missing', True, True), ('a', 'b', False, False),
    ('a', 'any', False, True), ('missing', 'any', False, True),
    ('any', 'missing', False, False), ('a', 'both', False, True),
    ('both', 'reordered', True, True), ('missing', 'unresolved', True, True),
])
def test_group_comparisons_respect_state_sets(field, a, b, unresolved, left, right, equal, specific):
    values = {'missing': [], 'a': [a], 'b': [b], 'any': ['x'],
              'both': [a, b], 'reordered': [b, a], 'unresolved': [unresolved]}
    labeled = GROUP.replace('1 N', '1 *1 N').replace('2 N', '2 *1 N')
    first = Group().from_adjacency_list(labeled)
    second = Group().from_adjacency_list(labeled)
    setattr(first, field, values[left])
    setattr(second, field, values[right])
    assert first.is_isomorphic(second) == equal
    assert bool(first.find_isomorphism(second)) == equal
    assert first.is_subgraph_isomorphic(second) == specific
    assert bool(first.find_subgraph_isomorphisms(second)) == specific
    assert first.is_subgraph_isomorphic(second, generate_initial_map=True) == specific


@pytest.mark.parametrize('state,electronic,vibrational', [
    (('', -1), [], []), (('A3Su+', -1), ['A3Su+'], []),
    (('', 0), [], [0]), (('', 1), [], [1]), (('A3Su+', 1), ['A3Su+'], [1]),
])
def test_to_group_carries_exact_state_constraints(state, electronic, vibrational):
    original = Molecule().from_adjacency_list(N2)
    original.electronic_state, original.vibrational_level = state
    group = original.to_group()
    assert group.electronic_state == electronic
    assert group.vibrational_level == vibrational
    assert original.is_subgraph_isomorphic(group)
    assert original.find_subgraph_isomorphisms(group)
