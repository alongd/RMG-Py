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

"""Raw pattern refusal and ground index-alias compatibility."""
import pytest

from rmgpy.data.base import Database, Entry, ForbiddenStructures
from rmgpy.exceptions import DatabaseError, SpeciesIdentityError
from rmgpy.data.kinetics.library import KineticsLibrary
from rmgpy.data.kinetics.depository import KineticsDepository
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species


def nitrogen(resolved=False):
    mol = Molecule(smiles='N#N')
    if resolved:
        mol.electronic_state = 'A'
        mol.vibrational_level = 1
    return mol


def load_patterns(database, path):
    if isinstance(database, ForbiddenStructures):
        database.load_old(str(path))
    else:
        database.load_old_dictionary(str(path), pattern=True)


class TestExcitedExportRework6:
    @pytest.mark.parametrize('cls', [Database, ForbiddenStructures])
    @pytest.mark.parametrize('order', [(True, False), (False, True), (True,)])
    @pytest.mark.parametrize('final_blank', [False, True])
    def test_old_pattern_record_refuses_state_before_store(self, tmp_path, cls, order, final_blank):
        records = [nitrogen(state).to_adjacency_list(label='N2') for state in order]
        path = tmp_path / 'forbidden.txt'
        path.write_text('\n'.join(records) + ('\n' if final_blank else ''))
        database = cls()
        with pytest.raises(SpeciesIdentityError, match='N2') as error:
            load_patterns(database, path)
        assert 'Database.load_old_dictionary' in str(error.value)
        assert 'electronicstate A' in str(error.value)
        assert 'vibrationallevel 1' in str(error.value)
        if order[0]:
            assert not database.entries
        else:
            assert database.entries['N2'].item == records[0]
        print('RAW PATTERN REFUSED BEFORE STORE', cls.__name__, order, str(error.value))

    @pytest.mark.parametrize('cls', [Database, ForbiddenStructures])
    @pytest.mark.parametrize('scenario', ['different', 'repeat', 'discarded-invalid'])
    def test_ground_pattern_duplicates_keep_historical_last_record(self, tmp_path, cls, scenario):
        last = Molecule(smiles='O=O').to_adjacency_list(label='ordinary')
        first = nitrogen().to_adjacency_list(label='ordinary') if scenario == 'different' else last if scenario == 'repeat' else 'ordinary\nnot an adjacency list\n'
        path = tmp_path / 'ground.txt'
        path.write_text(first + '\n' + last + '\n')
        database = cls()
        load_patterns(database, path)
        assert set(database.entries) == {'ordinary'}
        assert database.entries['ordinary'].item.to_adjacency_list() == '1 O u0 p2 c0 {2,D}\n2 O u0 p2 c0 {1,D}\n'
        print('GROUND PATTERN LAST WINS', cls.__name__, scenario)

    def test_resolved_molecular_dictionary_still_loads(self, tmp_path):
        mol = nitrogen(True)
        path = tmp_path / 'molecules.txt'
        path.write_text(mol.to_adjacency_list(label='N2'))
        database = Database()
        database.load_old_dictionary(str(path), pattern=False)
        assert database.entries['N2'].item.is_isomorphic(mol)

    @pytest.mark.parametrize('resolved_first,missing', [(False, 'A(1)'), (True, 'B(2)')])
    def test_display_index_alias_refuses_ground_target(self, tmp_path, resolved_first, missing):
        first, last = nitrogen(resolved_first), nitrogen()
        path = tmp_path / 'reactions.py'
        path.write_text('entry(index=1, label="A(1) <=> B(2)", kinetics=Arrhenius(A=(1,"s^-1")))\n')
        (tmp_path / 'dictionary.txt').write_text(first.to_adjacency_list(label='A') + '\n' + last.to_adjacency_list(label='B') + '\n')
        library = KineticsLibrary()
        with pytest.raises(DatabaseError) as error:
            library.load(str(path), local_context={'Arrhenius': Arrhenius})
        assert str(error.value) == 'Species ' + missing + ' in kinetics library  is missing from its dictionary.'
        print('GROUND INDEX ALIAS REFUSED', str(error.value))

    @pytest.mark.parametrize('indexed', [False, True])
    def test_exact_ground_dictionary_keys_keep_loading(self, tmp_path, indexed):
        first, last = ('A(1)', 'B(2)') if indexed else ('A', 'B')
        path = tmp_path / 'reactions.py'
        path.write_text('entry(index=1, label=' + repr(first + ' <=> ' + last) + ', kinetics=Arrhenius(A=(1,"s^-1")))\n')
        (tmp_path / 'dictionary.txt').write_text(nitrogen().to_adjacency_list(label=first) + '\n' + nitrogen().to_adjacency_list(label=last) + '\n')
        library = KineticsLibrary()
        library.load(str(path), local_context={'Arrhenius': Arrhenius})
        assert len(library.entries) == 1
        assert library.entries[1].item.reactants[0].label == first
        assert library.entries[1].item.products[0].label == last
        print('GROUND EXACT KEYS ACCEPTED', indexed)

    def test_resolved_dictionary_index_aliases_keep_loading(self, tmp_path):
        first, last = nitrogen(True), nitrogen(True)
        last.electronic_state = 'B'
        path = tmp_path / 'reactions.py'
        path.write_text('entry(index=1, label="A(1) <=> B(2)", kinetics=Arrhenius(A=(1,"s^-1")))\n')
        (tmp_path / 'dictionary.txt').write_text(first.to_adjacency_list(label='A') + '\n' + last.to_adjacency_list(label='B') + '\n')
        library = KineticsLibrary()
        library.load(str(path), local_context={'Arrhenius': Arrhenius})
        reaction = library.entries[1].item
        assert reaction.reactants[0].molecule[0].is_isomorphic(first)
        assert reaction.products[0].molecule[0].is_isomorphic(last)

    def test_resolved_export_index_save_reload_roundtrip(self, tmp_path):
        first, last = nitrogen(True), nitrogen(True)
        last.electronic_state = 'B'
        reactant = Species(label='A', index=1, molecule=[first])
        product = Species(label='B', index=2, molecule=[last])
        reaction = Reaction(reactants=[reactant], products=[product])
        library = KineticsLibrary()
        library.entries[1] = Entry(index=1, label='A(1) <=> B(2)', item=reaction, data=Arrhenius(A=(1, 's^-1')))
        output = tmp_path / 'reactions.py'
        library.save(str(output))
        assert output.read_text().startswith('# RMG-PAIR-GENERATION ')
        restored = KineticsLibrary()
        restored.load(str(output), local_context={'Arrhenius': Arrhenius})
        assert len(restored.entries) == 1
        reaction = next(iter(restored.entries.values())).item
        assert reaction.reactants[0].molecule[0].is_isomorphic(first)
        assert reaction.products[0].molecule[0].is_isomorphic(last)
        print('RESOLVED INDEX SAVE RELOAD PRESERVES', reaction.reactants[0].molecule[0].state_suffix(), reaction.products[0].molecule[0].state_suffix())

    @pytest.mark.parametrize('loader', [KineticsLibrary, KineticsDepository])
    @pytest.mark.parametrize('resolved_first', [False, True])
    def test_mixed_indexed_export_uses_exact_ground_dictionary_names(self, tmp_path, loader, resolved_first):
        first, last = nitrogen(resolved_first), nitrogen(not resolved_first)
        reactant = Species(label='A', index=1, molecule=[first])
        product = Species(label='B', index=2, molecule=[last])
        library = loader()
        library.entries[1] = Entry(index=1, label='A(1)  <=>  B(2)',
                                  item=Reaction(reactants=[reactant], products=[product]),
                                  data=Arrhenius(A=(1, 's^-1')))
        output = tmp_path / 'reactions.py'
        library.save(str(output))
        restored = loader()
        restored.load(str(output), local_context={'Arrhenius': Arrhenius})
        entry = next(iter(restored.entries.values()))
        assert entry.label == 'A <=> B'
        assert entry.item.reactants[0].molecule[0].is_isomorphic(first)
        assert entry.item.products[0].molecule[0].is_isomorphic(last)
        assert library.entries[1].label == 'A(1)  <=>  B(2)'
        print('MIXED INDEX SAVE RELOAD PRESERVES', loader.__name__, resolved_first, entry.label)

    @pytest.mark.parametrize('loader', [KineticsLibrary, KineticsDepository])
    def test_ground_indexed_entry_in_resolved_library_uses_exact_names(self, tmp_path, loader):
        ground = Species(label='A', index=1, molecule=[nitrogen()])
        excited = Species(label='B', index=2, molecule=[nitrogen(True)])
        library = loader()
        for index, label, product in [(1, 'A(1) <=> B(2)', excited), (2, 'A(1) <=> A(1)', ground)]:
            library.entries[index] = Entry(index=index, label=label,
                item=Reaction(reactants=[ground], products=[product]), data=Arrhenius(A=(1, 's^-1')))
        output = tmp_path / 'reactions.py'
        library.save(str(output))
        restored = loader()
        restored.load(str(output), local_context={'Arrhenius': Arrhenius})
        assert len(restored.entries) == 2
        assert {entry.label for entry in restored.entries.values()} == {'A <=> B', 'A <=> A'}
        assert any(entry.item.products[0].is_isomorphic(excited) for entry in restored.entries.values())
        print('MIXED LIBRARY GROUND ENTRY RELOADS', loader.__name__)
