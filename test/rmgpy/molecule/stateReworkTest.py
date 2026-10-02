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

"""Regression contract for the six independent state-identity review findings."""
import hashlib
import json
import pickle
import subprocess
import sys
from pathlib import Path

import pytest

from rmgpy.data.kinetics.family import KineticsFamily, complete_round_trip
from rmgpy.exceptions import InvalidAdjacencyListError
from rmgpy.molecule import Molecule
from rmgpy.molecule.fragment import Fragment
from rmgpy.molecule.group import Group
from rmgpy.molecule import translator
from rmgpy.reaction import Reaction
from rmgpy.species import Species

N2 = '1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}'


def resolved():
    return Molecule().from_adjacency_list('electronicstate A\nvibrationallevel 2\n' + N2)


def assert_isolated_reader(script, payload):
    # Invalid InChI layers can abort the native backend in the rejected commit.
    # A subprocess records that failure without aborting all other regressions.
    checkout = Path(__file__).resolve().parents[3]
    guard = '''
import json, sys
from pathlib import Path
import rmgpy
import pytest
from rmgpy.molecule import Molecule
from rmgpy.molecule.inchi import AugmentedInChI
from rmgpy.molecule import translator
assert Path(rmgpy.__file__).resolve().parent == Path(sys.argv[2]) / 'rmgpy'
data = json.loads(sys.argv[1])
'''
    result = subprocess.run([sys.executable, '-c', guard + script, json.dumps(payload), str(checkout)],
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr


class TestFinding1:
    @pytest.mark.parametrize('token', ['A|v:2', 'A/B', 'Σ', None, 'A\n'])
    def test_electronic_assignment_validates_before_mutation(self, token):
        m = resolved()
        original = m.to_adjacency_list()
        with pytest.raises(ValueError, match='electronic'):
            m.electronic_state = token
        assert m.to_adjacency_list() == original

    @pytest.mark.parametrize('level', [-2, True, 0.5, '2', None, 2**40])
    def test_vibrational_assignment_validates_before_mutation(self, level):
        m = resolved()
        with pytest.raises(ValueError, match='vibrational'):
            m.vibrational_level = level
        assert m.vibrational_level == 2

    @pytest.mark.parametrize('field,value', [('electronic_state', 'A'), ('vibrational_level', 0),
                                           ('electronic_state', 'A|v:2'), ('vibrational_level', -2)])
    def test_fragment_assignment_cannot_introduce_a_state(self, field, value):
        f = Fragment()
        with pytest.raises(ValueError, match='state|vibrational'):
            setattr(f, field, value)
        assert (f.electronic_state, f.vibrational_level) == ('', -1)


class TestFinding2:
    def test_connected_split_preserves_the_same_resolved_identity(self):
        original = resolved()
        pieces = original.split()
        assert len(pieces) == 1
        assert pieces[0].is_isomorphic(original)
        assert pieces[0].to_adjacency_list() == original.to_adjacency_list()

    def test_update_preserves_resolved_spin_and_state(self):
        original = Molecule(smiles='[O][O]', multiplicity=1, electronic_state='a1Dg', vibrational_level=2)
        original.update()
        assert (original.electronic_state, original.vibrational_level, original.multiplicity) == ('a1Dg', 2, 1)

    def test_update_can_initialize_an_unset_multiplicity_without_dropping_state(self):
        original = Molecule(atoms=Molecule(smiles='N#N').atoms, electronic_state='A')
        original.update()
        assert (original.electronic_state, original.vibrational_level, original.multiplicity) == ('A', -1, 1)

    @pytest.mark.parametrize('excited_first', [True, False])
    def test_merge_with_a_resolved_molecule_is_refused(self, excited_first):
        first, second = resolved(), Molecule(smiles='[He]')
        if not excited_first:
            first, second = second, first
        with pytest.raises(NotImplementedError, match='resolved state'):
            first.merge(second)

    def test_fragment_merge_cannot_drop_a_molecule_state(self):
        with pytest.raises(NotImplementedError, match='resolved state'):
            Fragment(smiles='C').merge(resolved())

    def test_disconnected_split_is_refused_for_a_resolved_molecule(self):
        original = Molecule(smiles='N#N.[He]', electronic_state='A')
        with pytest.raises(NotImplementedError, match='resolved state'):
            original.split()

    def test_resonance_hybrid_carries_state(self):
        species = Species(molecule=[Molecule(smiles='[CH2]C=C', electronic_state='A', vibrational_level=2)])
        hybrid = species.get_resonance_hybrid()
        assert len(species.molecule) >= 2
        assert (hybrid.electronic_state, hybrid.vibrational_level, hybrid.multiplicity) == ('A', 2, 2)

    @pytest.mark.parametrize('excited_first', [True, False])
    def test_family_reaction_matches_cannot_match_merged_resolved_reactants(self, excited_first):
        family = KineticsFamily()
        group = Group().from_adjacency_list('1 N ux {2,T}\n2 N ux {1,T}')
        ground = Species(molecule=[Molecule(smiles='N#N')])
        helium = Species(molecule=[Molecule(smiles='[He]')])
        assert family.reaction_matches(Reaction(reactants=[ground, helium]), group)
        reactants = [Species(molecule=[resolved()]), helium]
        if not excited_first:
            reactants.reverse()
        assert family.reaction_matches(Reaction(reactants=reactants), group) is False


class TestFinding3:
    @pytest.mark.parametrize('smiles', ['[Ar]', '[He]', 'N#N'])
    @pytest.mark.parametrize('identifier_type', ['smiles', 'inchi'])
    def test_constructor_declared_state_survives_all_identifier_paths(self, smiles, identifier_type):
        ground = Molecule(smiles=smiles)
        value = smiles if identifier_type == 'smiles' else ground.to_inchi()
        m = Molecule(**{identifier_type: value}, electronic_state='Ar(4s)', vibrational_level=2)
        assert (m.electronic_state, m.vibrational_level) == ('Ar(4s)', 2)

    @pytest.mark.parametrize('smiles', ['N#N', '[Ar]'])
    @pytest.mark.parametrize('reader', ['from_smiles', 'from_inchi', 'translator_smiles', 'translator_inchi'])
    def test_state_blind_replacement_resets_resolved_state(self, smiles, reader):
        m = resolved()
        _ = m.fingerprint
        _ = m.inchi
        _ = m.smiles
        ground = Molecule(smiles=smiles)
        if reader == 'from_smiles':
            m.from_smiles(smiles)
        elif reader == 'from_inchi':
            m.from_inchi(ground.to_inchi())
        elif reader == 'translator_smiles':
            translator.from_smiles(m, smiles)
        else:
            translator.from_inchi(m, ground.to_inchi())
        assert (m.electronic_state, m.vibrational_level) == ('', -1)
        assert m.fingerprint == ground.fingerprint
        assert m.inchi == ground.inchi
        assert m.smiles == ground.smiles
        assert m.is_isomorphic(ground)


class TestFinding4:
    @pytest.mark.parametrize('smiles,multiplicity', [('N#N', 1), ('[CH2]C=C', 2), ('[CH2]', 1)])
    @pytest.mark.parametrize('state', [('A', -1), ('', 0), ('A(4s)+', 2)])
    @pytest.mark.parametrize('as_object', [False, True])
    def test_augmented_inchi_reads_back_state_layers(self, smiles, multiplicity, state, as_object):
        if smiles == '[CH2]' and multiplicity == 1:
            original = Molecule().from_adjacency_list('1 C u0 p1 c0 {2,S} {3,S}\n2 H u0 p0 c0 {1,S}\n3 H u0 p0 c0 {1,S}')
            original.electronic_state, original.vibrational_level = state
        else:
            original = Molecule(smiles=smiles, multiplicity=multiplicity,
                                electronic_state=state[0], vibrational_level=state[1])
        identifier = original.to_augmented_inchi()
        assert_isolated_reader('''
identifier = AugmentedInChI(data['identifier']) if data['as_object'] else data['identifier']
original = Molecule().from_adjacency_list(data['adjacency'])
restored = Molecule()
if data['as_object']:
    translator.from_augmented_inchi(restored, identifier)
else:
    restored.from_augmented_inchi(identifier)
assert [restored.electronic_state, restored.vibrational_level] == data['state']
assert restored.is_isomorphic(original)
''', {'identifier': identifier, 'as_object': as_object, 'adjacency': original.to_adjacency_list(), 'state': state})

    @pytest.mark.parametrize('layers', ['/es:A|v:2', '/es:A/B', '/es:', '/v:-2', '/v:1.5',
                                       '/es:A/es:B', '/v:1/v:2'])
    def test_malformed_augmented_state_layers_have_a_named_refusal(self, layers):
        assert_isolated_reader('''
with pytest.raises(ValueError, match='electronic|vibrational|state'):
    Molecule().from_augmented_inchi(data)
''', 'InChI=1S/N2/c1-2' + layers)

    def test_standard_inchi_reader_refuses_state_layers(self):
        assert_isolated_reader('''
with pytest.raises(ValueError, match='state'):
    Molecule().from_inchi(data)
''', 'InChI=1S/C3H5/c1-3-2/h3H,1-2H2/u1/es:A/v:2')


class TestFinding5:
    @pytest.mark.parametrize('initially_resolved', [False, True])
    def test_completed_species_copy_keeps_then_invalidates_identity_caches(self, initially_resolved):
        m = resolved() if initially_resolved else Molecule(smiles='N#N')
        original = Species(molecule=[m])
        fingerprint = original.fingerprint
        augmented = original.get_augmented_inchi()
        copied = complete_round_trip(original)
        assert copied.fingerprint == fingerprint
        assert copied.aug_inchi == augmented
        copied.molecule[0].electronic_state = 'B'
        assert copied.fingerprint == copied.molecule[0].fingerprint
        assert copied.aug_inchi is None
        assert copied.get_augmented_inchi() == copied.molecule[0].to_augmented_inchi()

    @pytest.mark.parametrize('field,value,suffix', [('electronic_state', 'A', '/es:A'),
                                                  ('vibrational_level', 2, '/v:2')])
    def test_species_fingerprint_cache_tracks_state_mutation(self, field, value, suffix):
        m = Molecule(smiles='N#N')
        species = Species(molecule=[m])
        old_fingerprint = species.fingerprint
        old_hash = hash(species)
        setattr(m, field, value)
        assert species.fingerprint == old_fingerprint + suffix.replace('/', '|')
        assert hash(species) != old_hash

    @pytest.mark.parametrize('field,value,suffix', [('electronic_state', 'A', '/es:A'),
                                                  ('vibrational_level', 2, '/v:2')])
    def test_species_augmented_inchi_cache_tracks_state_mutation(self, field, value, suffix):
        m = Molecule(smiles='N#N')
        species = Species(molecule=[m])
        previous = species.get_augmented_inchi()
        setattr(m, field, value)
        assert species.aug_inchi is None
        assert species.get_augmented_inchi() == previous + suffix

    def test_species_caches_follow_reset_to_unresolved(self):
        m = resolved()
        species = Species(molecule=[m])
        _ = species.fingerprint
        _ = species.get_augmented_inchi()
        m.electronic_state = ''
        m.vibrational_level = -1
        assert species.aug_inchi is None
        assert species.fingerprint == Molecule(smiles='N#N').fingerprint
        assert species.get_augmented_inchi() == Molecule(smiles='N#N').to_augmented_inchi()


class TestFinding6:
    @pytest.mark.parametrize('token', ['2P1/2', 'A|v:2', 'A:B', 'A B', 'A\\B', 'Σ', 'A@B',
                                     'A\vB', 'A\fB', 'A\u2028B'])
    def test_real_malformed_tokens_are_refused_without_an_extra_suffix(self, token):
        with pytest.raises(InvalidAdjacencyListError, match='electronicstate'):
            Molecule().from_adjacency_list('electronicstate ' + token + '\n' + N2)
        assert Molecule().from_adjacency_list('electronicstate A\n' + N2).electronic_state == 'A'

    def test_unresolved_reduce_shape_is_frozen_at_ccba0691c(self):
        m = Molecule(smiles='N#N')
        rebuild, args, state = m.__reduce__()
        assert (rebuild.__module__, rebuild.__name__) == ('rmgpy.molecule.molecule', '_rebuild_molecule')
        assert args == (Molecule,)
        assert len(state) == 3
        assert len(state[0]) == 8
        assert state[0][:2] == (m.vertices, m.ordered_vertices)
        assert state[0][2:] == (-1.0, 1, True, {}, '', '')
        assert state[1:] == ({}, [])

    def test_unresolved_pickle_bytes_are_frozen_at_ccba0691c(self):
        payload = pickle.dumps(Molecule(smiles='N#N'), protocol=4)
        assert len(payload) == 380
        assert hashlib.sha256(payload).hexdigest() == 'a4f97200e9e5d6e48f43adb863e3172b1c67c521e6b0f6bf39ae36d3988c2ae4'

    def test_unresolved_adjacency_text_is_frozen_at_ccba0691c(self):
        assert Molecule(smiles='N#N').to_adjacency_list() == '1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}\n'
