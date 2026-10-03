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

"""Entry equations must be checked against the accompanying full dictionary."""
import io
import json
from pathlib import Path

import pytest

from rmgpy.data.base import Entry
from rmgpy.data.kinetics.common import library_reaction_equation, save_entry
from rmgpy.data.kinetics.depository import KineticsDepository
from rmgpy.data.kinetics.family import KineticsFamily
from rmgpy.data.kinetics.library import KineticsLibrary
from rmgpy.data.kinetics.rules import KineticsRules
from rmgpy.exceptions import ResolvedStateTrainingError, SpeciesIdentityError
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


def entry(reactant, product, index=1, label=None):
    reaction = Reaction(reactants=[reactant], products=[product])
    return Entry(index=index, label=label or str(reaction), item=reaction, data=Arrhenius(A=(index, 's^-1')))


class TestExcitedExportRework8:
    @pytest.mark.parametrize('shadow_state', [None, 'C'])
    @pytest.mark.parametrize('side', ['reactant', 'product'])
    @pytest.mark.parametrize('loader', [KineticsLibrary, KineticsDepository])
    def test_public_entry_writer_uses_the_whole_dictionary(self, tmp_path, loader, shadow_state, side):
        a, b = nitrogen('A', 'A', 1), nitrogen('B', 'B', 2)
        shadow_name = 'A(1)' if side == 'reactant' else 'B(2)'
        shadow, c = nitrogen(shadow_name, shadow_state), nitrogen('C')
        library = loader()
        library.entries[1] = entry(a, b)
        library.entries[2] = entry(shadow, c, 2)
        stream = io.StringIO()
        for item in library.entries.values():
            library.save_entry(stream, item)
        path = tmp_path / 'reactions.py'
        path.write_text(stream.getvalue())
        library.save_dictionary(str(tmp_path / 'dictionary.txt'))
        restored = loader()
        restored.load(str(path), local_context={'Arrhenius': Arrhenius})
        first, second = sorted(restored.entries.values(), key=lambda item: item.index)
        assert first.item.reactants[0].is_isomorphic(a)
        assert first.item.products[0].is_isomorphic(b)
        assert second.item.reactants[0].is_isomorphic(shadow)
        assert library.entries[1].label == 'A(1) <=> B(2)'

    @pytest.mark.parametrize('shadow_state', [None, 'C'])
    @pytest.mark.parametrize('side', ['reactant', 'product'])
    def test_training_appender_checks_retained_dictionary_names(self, tmp_path, shadow_state, side):
        family = KineticsFamily(label=str(tmp_path / 'family'))
        old = KineticsDepository(label=family.label + '/training')
        shadow_name = 'N2-2(1)' if side == 'reactant' else 'N2-3(2)'
        shadow = nitrogen(shadow_name, shadow_state)
        old.entries[1] = entry(nitrogen('N2'), shadow)
        path = tmp_path / 'family/training/reactions.py'
        old.save(str(path))
        family.depositories.append(old)
        a, b = nitrogen('incoming', 'A', 1), nitrogen('product', 'B', 2)
        reaction = Reaction(reactants=[a], products=[b], kinetics=Arrhenius(A=(2, 's^-1')))
        before = {item.name: item.read_bytes() for item in path.parent.iterdir()}
        with pytest.raises(ResolvedStateTrainingError, match='do not support resolved species'):
            family.save_training_reactions([reaction])
        assert {item.name: item.read_bytes() for item in path.parent.iterdir()} == before
        restored = KineticsDepository()
        restored.load(str(path), local_context={'Arrhenius': Arrhenius})
        assert len(restored.entries) == 1
        first = next(iter(restored.entries.values()))
        assert first.item.products[0].is_isomorphic(shadow)

    @pytest.mark.parametrize('shadow_state', [None, 'C'])
    def test_training_dictionary_context_includes_unused_retained_species(self, tmp_path, shadow_state):
        family = KineticsFamily(label=str(tmp_path / 'family'))
        old = KineticsDepository(label=family.label + '/training')
        ground = nitrogen('N2')
        old.entries[1] = entry(ground, ground)
        path = tmp_path / 'family/training/reactions.py'
        old.save(str(path))
        # A retained dictionary record need not occur in the loaded depository.
        shadow = nitrogen('N2-2(1)', shadow_state)
        dictionary = path.parent / 'dictionary.txt'
        with dictionary.open('a') as stream:
            stream.write(shadow.molecule[0].to_adjacency_list(label=shadow.label) + '\n')
        family.depositories.append(old)
        a, b = nitrogen('incoming', 'A', 1), nitrogen('product', 'B', 2)
        reaction = Reaction(reactants=[a], products=[b], kinetics=Arrhenius(A=(2, 's^-1')))
        before = {item.name: item.read_bytes() for item in path.parent.iterdir()}
        with pytest.raises(ResolvedStateTrainingError, match='do not support resolved species'):
            family.save_training_reactions([reaction])
        assert {item.name: item.read_bytes() for item in path.parent.iterdir()} == before
        restored = KineticsDepository()
        restored.load(str(path), local_context={'Arrhenius': Arrhenius})
        assert len(restored.entries) == 1
        assert restored.get_species(str(dictionary))[shadow.label].is_isomorphic(shadow)

    @pytest.mark.parametrize('writer', ['common', 'family', 'rules', 'equation'])
    def test_resolved_entry_without_dictionary_context_refuses(self, writer):
        item = entry(nitrogen('A', 'A', 1), nitrogen('B', 'B', 2))
        stream = io.StringIO()
        calls = {
            'common': lambda: save_entry(stream, item),
            'family': lambda: KineticsFamily().save_entry(stream, item),
            'rules': lambda: KineticsRules().save_entry(stream, item),
            'equation': lambda: library_reaction_equation(item),
        }
        with pytest.raises(SpeciesIdentityError, match='emitted dictionary') as error:
            calls[writer]()
        assert 'A' in str(error.value) and '|es:A|v:1' in str(error.value)
        assert stream.getvalue() == ''

    @pytest.mark.parametrize('writer', ['common', 'family'])
    def test_explicit_dictionary_context_preserves_resolved_entry(self, tmp_path, writer):
        a, b, shadow = nitrogen('A', 'A', 1), nitrogen('B', 'B', 2), nitrogen('A(1)')
        library = KineticsLibrary()
        library.entries[1] = entry(a, b)
        library.entries[2] = entry(shadow, nitrogen('C'), 2)
        stream = io.StringIO()
        declarations = library.dictionary_references()
        for item in library.entries.values():
            if writer == 'common':
                save_entry(stream, item, declarations=declarations)
            else:
                KineticsFamily().save_entry(stream, item, declarations=declarations)
        path = tmp_path / 'reactions.py'
        path.write_text(stream.getvalue())
        library.save_dictionary(str(tmp_path / 'dictionary.txt'))
        restored = KineticsLibrary()
        restored.load(str(path), local_context={'Arrhenius': Arrhenius})
        assert restored.entries[1].item.reactants[0].is_isomorphic(a)
        assert restored.entries[1].item.products[0].is_isomorphic(b)

    @pytest.mark.parametrize('populated', [False, True])
    def test_standalone_ground_library_entry_retains_historical_bytes(self, populated):
        # Use literal historical bytes; the supplied entry need not be registered.
        fixture = json.loads((Path(__file__).resolve().parents[1] /
                              'test_data/excited_export_ground_conflict.json').read_text())
        payload = fixture['KineticsLibrary:False']['files']['reactions.py']
        expected = payload[payload.index('entry(\n'):]
        a = nitrogen('A')
        b = Species(label='A', molecule=[Molecule(smiles='O=O')])
        supplied = entry(a, b, index=0)
        supplied.data = Arrhenius(A=(1, 's^-1'))  # Match the historical fixture before save reindexed its entry.
        library = KineticsLibrary()
        if populated:
            library.entries[1] = entry(nitrogen('unrelated'), nitrogen('other'))
        stream = io.StringIO()
        library.save_entry(stream, supplied)
        assert stream.getvalue() == expected

    @pytest.mark.parametrize('loader', [KineticsLibrary, KineticsDepository])
    @pytest.mark.parametrize('reverse', [False, True])
    def test_ground_exact_name_conflict_matches_historical_save_and_reload(self, tmp_path, loader, reverse):
        fixture = json.loads((Path(__file__).resolve().parents[1] /
                              'test_data/excited_export_ground_conflict.json').read_text())
        expected = fixture[loader.__name__ + ':' + str(reverse)]
        graphs = ['N#N', 'O=O'][::(-1 if reverse else 1)]
        a, b = [Species(label='A', molecule=[Molecule(smiles=graph)]) for graph in graphs]
        library = loader()
        library.entries[1] = entry(a, b)
        path = tmp_path / 'reactions.py'
        library.save(str(path))
        for name, payload in expected['files'].items():
            assert strip_generation_marker((tmp_path / name).read_text()) == payload
        restored = loader()
        restored.load(str(path), local_context={'Arrhenius': Arrhenius})
        first = min(restored.entries.values(), key=lambda item: item.index)
        actual = {'count': len(restored.entries), 'label': first.label,
                  'identities': [spc.molecule[0].to_adjacency_list()
                                 for spc in first.item.reactants + first.item.products]}
        assert actual == expected['reload']
        assert first.item.reactants[0].is_isomorphic(a)
        assert first.item.products[0].is_isomorphic(a)

    @pytest.mark.parametrize('loader', [KineticsLibrary, KineticsDepository])
    @pytest.mark.parametrize('reverse', [False, True])
    def test_resolved_exact_name_conflict_refuses_before_replacing_pair(self, tmp_path, loader, reverse):
        a, b = nitrogen('A'), nitrogen('A', 'A')
        if reverse:
            a, b = b, a
        library = loader()
        library.entries[1] = entry(a, b)
        path, dictionary = tmp_path / 'reactions.py', tmp_path / 'dictionary.txt'
        path.write_text('previous reactions')
        dictionary.write_text('previous dictionary')
        with pytest.raises(SpeciesIdentityError, match='collides') as error:
            library.save(str(path))
        assert "state=['']" in str(error.value) and '|es:A|v:1' in str(error.value)
        assert path.read_text() == 'previous reactions'
        assert dictionary.read_text() == 'previous dictionary'
