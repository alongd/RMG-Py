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

"""Saved library equation tokens must retain their declared full identities."""
import json
from pathlib import Path

import pytest

from rmgpy.data.base import Entry
from rmgpy.data.kinetics.depository import KineticsDepository
from rmgpy.data.kinetics.library import KineticsLibrary
from rmgpy.exceptions import SpeciesIdentityError
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.util import strip_generation_marker


def nitrogen(label, state=None, index=-1):
    molecule = Molecule(smiles='N#N')
    if state is not None:
        molecule.electronic_state = state
        molecule.vibrational_level = 1
    return Species(label=label, index=index, molecule=[molecule])


class TestExcitedExportRework7:
    @pytest.mark.parametrize('loader', [KineticsLibrary, KineticsDepository])
    @pytest.mark.parametrize('shadow_state', [None, 'C'])
    @pytest.mark.parametrize('side', ['reactant', 'product'])
    def test_saved_display_index_does_not_select_a_shadow_dictionary_species(self, tmp_path, loader, shadow_state, side):
        a, b = nitrogen('A', 'A', 1), nitrogen('B', 'B', 2)
        shadow_name = 'A(1)' if side == 'reactant' else 'B(2)'
        shadow, ground = nitrogen(shadow_name, shadow_state), nitrogen('C')
        library = loader()
        for index, label, reactant, product in [(1, 'A(1) <=> B(2)', a, b),
                                               (2, shadow_name + ' <=> C', shadow, ground)]:
            library.entries[index] = Entry(index=index, label=label,
                item=Reaction(reactants=[reactant], products=[product]), data=Arrhenius(A=(index, 's^-1')))
        path = tmp_path / 'reactions.py'
        library.save(str(path))
        restored = loader()
        restored.load(str(path), local_context={'Arrhenius': Arrhenius})
        first, second = sorted(restored.entries.values(), key=lambda entry: entry.index)
        assert first.item.reactants[0].is_isomorphic(a)
        assert first.item.products[0].is_isomorphic(b)
        assert second.item.reactants[0].is_isomorphic(shadow)
        selected = first.item.reactants[0] if side == 'reactant' else first.item.products[0]
        assert not selected.is_isomorphic(shadow)
        assert library.entries[1].label == 'A(1) <=> B(2)'
        print('SHADOW ROUNDTRIP', loader.__name__, side, shadow_state,
              selected.molecule[0].state_suffix(), first.label)

    def test_display_index_does_not_merge_with_exact_ground_collider_name(self, tmp_path):
        a, b, collider = nitrogen('A', 'A', 1), nitrogen('B', 'B', 2), nitrogen('A(1)', index=3)
        library = KineticsLibrary()
        label = 'A(1) (+A(1)) <=> B(2) (+A(1))'
        library.entries[1] = Entry(index=1, label=label,
            item=Reaction(reactants=[a], products=[b], specific_collider=collider),
            data=Arrhenius(A=(1, 's^-1')))
        path = tmp_path / 'reactions.py'
        library.save(str(path))
        restored = KineticsLibrary()
        restored.load(str(path), local_context={'Arrhenius': Arrhenius})
        reaction = next(iter(restored.entries.values())).item
        assert reaction.reactants[0].is_isomorphic(a)
        assert reaction.products[0].is_isomorphic(b)
        assert reaction.specific_collider.is_isomorphic(collider)
        assert reaction.reactants[0] is not reaction.specific_collider
        assert library.entries[1].label == label
        print('COLLIDER ROUNDTRIP', reaction.reactants[0].molecule[0].state_suffix(),
              'merged', reaction.reactants[0] is reaction.specific_collider)

    @pytest.mark.parametrize('loader', [KineticsLibrary, KineticsDepository])
    def test_safe_resolved_display_tokens_keep_their_equation_bytes(self, tmp_path, loader):
        a, b = nitrogen('A', 'A', 1), nitrogen('B', 'B', 2)
        library = loader()
        label = 'A(1)  <=>  B(2)'
        library.entries[1] = Entry(index=1, label=label,
            item=Reaction(reactants=[a], products=[b]), data=Arrhenius(A=(1, 's^-1')))
        path = tmp_path / 'reactions.py'
        library.save(str(path))
        restored = loader()
        restored.load(str(path), local_context={'Arrhenius': Arrhenius})
        entry = next(iter(restored.entries.values()))
        assert entry.label == label
        assert entry.item.reactants[0].is_isomorphic(a)
        assert entry.item.products[0].is_isomorphic(b)
        print('SAFE INDEX TOKENS RETAINED', loader.__name__, entry.label)

    def test_repeated_display_token_preserves_both_typed_participants(self, tmp_path):
        a, b, shadow = nitrogen('A', 'A', 1), nitrogen('B', 'B', 2), nitrogen('A(1)')
        library = KineticsLibrary()
        library.entries[1] = Entry(index=1, label='A(1) + A(1) <=> B(2) + A(1)',
            item=Reaction(reactants=[a, shadow], products=[b, shadow]),
            data=Arrhenius(A=(1, 'm^3/(mol*s)')))
        path = tmp_path / 'reactions.py'
        library.save(str(path))
        restored = KineticsLibrary()
        restored.load(str(path), local_context={'Arrhenius': Arrhenius})
        reaction = next(iter(restored.entries.values())).item
        assert sum(spc.is_isomorphic(a) for spc in reaction.reactants) == 1
        assert sum(spc.is_isomorphic(shadow) for spc in reaction.reactants) == 1
        assert sum(spc.is_isomorphic(b) for spc in reaction.products) == 1
        assert sum(spc.is_isomorphic(shadow) for spc in reaction.products) == 1
        print('REPEATED TOKEN PRESERVES BOTH IDENTITIES')

    @pytest.mark.parametrize('loader', [KineticsLibrary, KineticsDepository])
    def test_no_unique_dictionary_name_refuses_before_replacement(self, tmp_path, loader):
        a, shadow = nitrogen('A', 'A', 1), nitrogen('A', index=3)
        library = loader()
        library.entries[1] = Entry(index=1, label='A(1) <=> A(3)',
            item=Reaction(reactants=[a], products=[shadow]), data=Arrhenius(A=(1, 's^-1')))
        path = tmp_path / 'reactions.py'
        dictionary = tmp_path / 'dictionary.txt'
        path.write_text('previous reactions')
        dictionary.write_text('previous dictionary')
        with pytest.raises(SpeciesIdentityError) as error:
            library.save(str(path))
        assert 'identifier "A"' in str(error.value)
        assert 'index=1' in str(error.value) and 'index=3' in str(error.value)
        assert '|es:A|v:1' in str(error.value) and "state=['']" in str(error.value)
        assert path.read_text() == 'previous reactions'
        assert dictionary.read_text() == 'previous dictionary'
        print('NO UNIQUE NAME REFUSED', loader.__name__, str(error.value))

    @pytest.mark.parametrize('loader', [KineticsLibrary, KineticsDepository])
    def test_ground_shadow_keeps_exact_historical_files_and_reload(self, tmp_path, loader):
        baseline = json.loads((Path(__file__).resolve().parents[1] /
                               'test_data/excited_export_ground_shadow.json').read_text())[loader.__name__]
        a, b = nitrogen('A', index=1), nitrogen('B')
        shadow = Species(label='A(1)', molecule=[Molecule(smiles='[N]=[N]')])
        other = Species(label='C', molecule=[Molecule(smiles='[N]=[N]')])
        library = loader(name='ground shadow')
        for index, label, reactant, product in [(1, 'A(1)  <=>  B', a, b),
                                               (2, 'A(1) <=> C', shadow, other)]:
            library.entries[index] = Entry(index=index, label=label,
                item=Reaction(reactants=[reactant], products=[product]), data=Arrhenius(A=(index, 's^-1')))
        path = tmp_path / 'reactions.py'
        library.save(str(path))
        for name, expected in baseline['files'].items():
            assert strip_generation_marker((tmp_path / name).read_text()) == expected
        restored = loader()
        restored.load(str(path), local_context={'Arrhenius': Arrhenius})
        entries = sorted(restored.entries.values(), key=lambda entry: entry.index)
        assert [entry.label for entry in entries] == baseline['reload']['labels']
        assert [[spc.molecule[0].to_adjacency_list() for spc in entry.item.reactants + entry.item.products]
                for entry in entries] == baseline['reload']['identities']
        assert entries[0].item.reactants[0].is_isomorphic(shadow)
        assert not entries[0].item.reactants[0].is_isomorphic(a)
        print('GROUND SHADOW FILES AND RELOAD MATCH 0e69e741b', loader.__name__)
