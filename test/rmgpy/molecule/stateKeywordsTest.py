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

"""Public identity and serialization contract for optional molecule states."""
import copy
import itertools
import pickle

import pytest

from rmgpy.exceptions import InvalidAdjacencyListError, SpeciesError
from rmgpy.molecule import Molecule
from rmgpy.molecule.adjlist import from_adjacency_list, to_adjacency_list
from rmgpy.molecule.filtration import filter_structures
from rmgpy.molecule.fragment import Fragment
from rmgpy.molecule.group import Group
from rmgpy.molecule.translator import to_inchi, to_inchi_key
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.species import Species

N2 = '1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}'
STATES = [('', -1), ('', 0), ('', 1), ('A3Su+', -1), ('A3Su+', 2)]


def molecule(electronic_state='', vibrational_level=-1):
    headers = ''
    if electronic_state:
        headers += f'electronicstate {electronic_state}\n'
    if vibrational_level >= 0:
        headers += f'vibrationallevel {vibrational_level}\n'
    return Molecule().from_adjacency_list(headers + N2)


@pytest.mark.parametrize('left_state,right_state', list(itertools.combinations(STATES, 2)))
@pytest.mark.parametrize('strict', [True, False])
def test_states_have_distinct_identity(left_state, right_state, strict):
    left, right = molecule(*left_state), molecule(*right_state)
    left.assign_atom_ids()
    for a, b in zip(left.atoms, right.atoms):
        b.id = a.id
    assert not left.is_isomorphic(right, strict=strict)
    assert not right.is_isomorphic(left, strict=strict)
    assert left.find_isomorphism(right, strict=strict) == []
    assert not left.is_identical(right, strict=strict)
    assert left != right
    assert left.fingerprint != right.fingerprint
    assert hash(left) != hash(right)


@pytest.mark.parametrize('state', STATES)
@pytest.mark.parametrize('transport', ['adjacency', 'shallow', 'deep', 'copy', 'deepcopy', 'pickle', 'repr'])
def test_state_round_trips(state, transport):
    original = molecule(*state)
    if transport == 'adjacency':
        restored = Molecule().from_adjacency_list(original.to_adjacency_list())
    elif transport == 'shallow':
        restored = original.copy()
        assert restored.atoms[0] is original.atoms[0]
    elif transport == 'deep':
        restored = original.copy(deep=True)
        assert restored.atoms[0] is not original.atoms[0]
    elif transport == 'copy':
        restored = copy.copy(original)
    elif transport == 'deepcopy':
        restored = copy.deepcopy(original)
    elif transport == 'pickle':
        restored = pickle.loads(pickle.dumps(original))
    else:
        restored = eval(repr(original), {'Molecule': Molecule})
    assert (restored.electronic_state, restored.vibrational_level) == state
    assert restored.is_isomorphic(original)
    assert restored.find_isomorphism(molecule(*state))
    assert restored.is_identical(original)
    assert restored.to_adjacency_list() == original.to_adjacency_list()


@pytest.mark.parametrize('state', STATES)
@pytest.mark.parametrize('source', ['smiles', 'inchi', 'atoms'])
def test_constructor_carries_state(state, source):
    ground = molecule()
    kwargs = {source: ground.to_smiles() if source == 'smiles' else ground.to_inchi()}
    if source == 'atoms':
        kwargs = {'atoms': ground.copy(deep=True).atoms, 'multiplicity': ground.multiplicity}
    original = Molecule(electronic_state=state[0], vibrational_level=state[1], **kwargs)
    assert (original.electronic_state, original.vibrational_level) == state
    assert original.is_isomorphic(molecule(*state))


@pytest.mark.parametrize('token', ['A3Su+', '2P1_2', 'Ar(4s)', 'A-B', 'a,b.c', '+-_.,()'])
def test_token_charset_is_accepted(token):
    assert molecule(token).electronic_state == token
    assert Molecule(electronic_state=token).electronic_state == token


@pytest.mark.parametrize('token', ['2P1/2', 'A|B', 'A:B', 'A B', 'A\\B', 'Σ', 'A\n', 'A\t', 'A@B'])
def test_invalid_electronic_token_is_refused(token):
    with pytest.raises(ValueError, match='electronic'):
        Molecule(electronic_state=token)
    # Trailing header whitespace is formatting, while constructor tokens are exact.
    if token not in ('A\n', 'A\t'):
        with pytest.raises(InvalidAdjacencyListError, match='electronicstate'):
            Molecule().from_adjacency_list('electronicstate ' + token + '\n' + N2)


@pytest.mark.parametrize('header', ['electronicstate', 'electronicstate A/B', 'electronicstate A|B',
                                  'electronicstate A\\B', 'electronicstate Σ', 'vibrationallevel',
                                  'vibrationallevel -1', 'vibrationallevel 1.5', 'vibrationallevel two'])
def test_invalid_state_header_is_refused(header):
    with pytest.raises(InvalidAdjacencyListError, match='electronicstate|vibrationallevel'):
        Molecule().from_adjacency_list(header + '\n' + N2)


@pytest.mark.parametrize('value', [-2, 1.5, '1', None, True])
def test_invalid_constructor_vibrational_level_is_refused(value):
    with pytest.raises(ValueError, match='vibrational'):
        Molecule(vibrational_level=value)


def test_electronic_tokens_are_case_sensitive():
    assert not molecule('A').is_isomorphic(molecule('a'))


@pytest.mark.parametrize('header', ['electronicstate A3Su+', 'vibrationallevel 1'])
def test_group_and_state_blind_parser_refuse_keywords(header):
    with pytest.raises(InvalidAdjacencyListError, match='group'):
        Group().from_adjacency_list(header + '\n1 N ux {2,T}\n2 N ux {1,T}')
    with pytest.raises(InvalidAdjacencyListError, match='state'):
        from_adjacency_list(header + '\n' + N2)


def test_low_level_parser_preserves_four_item_return_tuple():
    state = {}
    result = from_adjacency_list('electronicstate A3Su+\nvibrationallevel 2\n' + N2, state=state)
    assert len(result) == 4
    assert state == {'electronic_state': 'A3Su+', 'vibrational_level': 2}


@pytest.mark.parametrize('state', STATES[1:])
def test_resolved_molecules_only_match_their_state_constrained_group(state):
    group = Group().from_adjacency_list('1 N ux {2,T}\n2 N ux {1,T}')
    ground = molecule()
    assert ground.is_subgraph_isomorphic(group)
    assert ground.find_subgraph_isomorphisms(group)
    resolved = molecule(*state)
    assert not resolved.is_subgraph_isomorphic(group)
    assert resolved.find_subgraph_isomorphisms(group) == []
    resolved_group = resolved.to_group()
    assert resolved.is_subgraph_isomorphic(resolved_group)
    assert not ground.is_subgraph_isomorphic(resolved_group)


@pytest.mark.parametrize('state', STATES[1:])
def test_resolved_state_cannot_be_written_in_old_or_group_format(state):
    resolved = molecule(*state)
    with pytest.raises(InvalidAdjacencyListError, match='[Oo]ld'):
        resolved.to_adjacency_list(old_style=True)
    with pytest.raises(InvalidAdjacencyListError, match='[Oo]ld'):
        to_adjacency_list(resolved.atoms, resolved.multiplicity, old_style=True,
                          electronic_state=state[0], vibrational_level=state[1])
    with pytest.raises(InvalidAdjacencyListError, match='group'):
        to_adjacency_list(resolved.atoms, [1], group=True,
                          electronic_state=state[0], vibrational_level=state[1])


@pytest.mark.parametrize('header', ['electronicstate A3Su+', 'vibrationallevel 1'])
def test_fragment_parser_refuses_state(header):
    with pytest.raises(InvalidAdjacencyListError, match='state'):
        Fragment().from_adjacency_list(header + '\n' + N2)


@pytest.mark.parametrize('kwargs', [{'electronic_state': 'A3Su+'}, {'vibrational_level': 0}])
def test_fragment_constructor_refuses_state(kwargs):
    with pytest.raises(ValueError, match='Fragment.*state'):
        Fragment(**kwargs)


@pytest.mark.parametrize('state', STATES[1:])
@pytest.mark.parametrize('smiles', ['N#N', '[CH2]C=C', 'c1ccccc1', 'NC=O'])
def test_resonance_carries_state(state, smiles):
    original = Molecule(smiles=smiles, electronic_state=state[0], vibrational_level=state[1])
    species = Species(molecule=[original])
    species.generate_resonance_structures()
    assert species.molecule
    if smiles == '[CH2]C=C':
        assert len(species.molecule) >= 2
    assert all((m.electronic_state, m.vibrational_level) == state for m in species.molecule)


@pytest.mark.parametrize('state', STATES[1:])
def test_species_and_filter_reject_mixed_states(state):
    ground, resolved = molecule(), molecule(*state)
    with pytest.raises(SpeciesError, match='states'):
        Species(molecule=[ground, resolved])
    with pytest.raises(ValueError, match='states'):
        filter_structures([ground, resolved])


@pytest.mark.parametrize('state,suffix,key_suffix', [
    (('', -1), '', ''), (('', 0), '/v:0', '-v0'), (('', 1), '/v:1', '-v1'),
    (('A3Su+', -1), '/es:A3Su+', '-es413353752b'),
    (('A3Su+', 2), '/es:A3Su+/v:2', '-es413353752b-v2')])
def test_augmented_exports_include_state_and_standard_exports_do_not(state, suffix, key_suffix):
    ground, resolved = molecule(), molecule(*state)
    assert resolved.to_smiles() == ground.to_smiles()
    assert resolved.to_inchi() == ground.to_inchi()
    assert resolved.to_inchi_key() == ground.to_inchi_key()
    assert resolved.to_augmented_inchi() == ground.to_augmented_inchi() + suffix
    assert resolved.to_augmented_inchi_key() == ground.to_augmented_inchi_key() + key_suffix
    for level in (1, 2):
        assert to_inchi(resolved, aug_level=level) == to_inchi(ground, aug_level=level) + suffix
        assert to_inchi_key(resolved, aug_level=level) == to_inchi_key(ground, aug_level=level) + key_suffix


def test_model_deduplicates_a_species_list_by_state_in_cache_and_formula_bucket():
    model = CoreEdgeReactionModel()
    # Exercise identity indexing during input construction. Full-model manifold
    # partition validation is covered separately by vibrationalManifoldTest.
    model.defer_vibrational_validation = True
    specs = []
    for state in STATES:
        spec, is_new = model.make_new_species(molecule(*state), generate_thermo=False)
        assert is_new
        specs.append(spec)
    assert len(model.species_dict['N2']) == len(STATES)
    assert len(set(specs)) == len(STATES)
    for state, expected in zip(STATES, specs):
        spec, is_new = model.make_new_species(molecule(*state), generate_thermo=False)
        assert not is_new
        assert spec is expected
    model.species_cache = [None] * len(model.species_cache)
    for state, expected in zip(STATES, specs):
        assert model.check_for_existing_species(molecule(*state)) is expected
    restored = pickle.loads(pickle.dumps(specs))
    assert [(s.molecule[0].electronic_state, s.molecule[0].vibrational_level) for s in restored] == STATES


@pytest.mark.parametrize('header', ['electronicstate A3Su+', 'vibrationallevel 1'])
def test_unlabeled_state_headers_are_not_species_labels(header):
    assert Species().from_adjacency_list(header + '\n' + N2).label == ''
    assert Species().from_adjacency_list('N2_excited\n' + header + '\n' + N2).label == 'N2_excited'


def test_reusing_molecule_resets_state_and_fingerprint():
    item = molecule('A3Su+', 2)
    assert item.fingerprint.endswith('|es:A3Su+|v:2')
    item.from_adjacency_list(N2)
    assert (item.electronic_state, item.vibrational_level) == ('', -1)
    assert item.fingerprint == molecule().fingerprint


def test_public_state_changes_are_reflected_in_identity():
    item = molecule()
    assert item.fingerprint == molecule().fingerprint
    item.electronic_state = 'A3Su+'
    item.vibrational_level = 2
    assert item.fingerprint == molecule('A3Su+', 2).fingerprint
    assert hash(item) == hash(molecule('A3Su+', 2))
    assert item.is_isomorphic(molecule('A3Su+', 2))


@pytest.mark.parametrize('header', ['electronicstate A', 'vibrationallevel 1'])
def test_duplicate_state_headers_are_refused(header):
    with pytest.raises(InvalidAdjacencyListError, match='[Dd]uplicate'):
        Molecule().from_adjacency_list(header + '\n' + header + '\n' + N2)


@pytest.mark.parametrize('header', ['electronicstate A', 'vibrationallevel 1'])
def test_old_style_reader_refuses_state(header):
    old = '1 N 0 1 {2,T}\n2 N 0 1 {1,T}'
    with pytest.raises(InvalidAdjacencyListError, match='[Oo]ld'):
        Molecule().from_adjacency_list(header + '\n' + old)
