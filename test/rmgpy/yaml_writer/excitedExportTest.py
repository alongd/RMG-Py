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

import itertools
import re

import pytest

from rmgpy.chemkin import get_species_identifier, load_chemkin_file, save_chemkin_file, save_species_dictionary
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.kinetics import Arrhenius
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial


@pytest.fixture
def state_species():
    species = []
    for index, (smiles, label, electronic, vibrational) in enumerate([
        ('N#N', 'N2', '', -1),
        ('N#N', 'N2(v=1)', '', 1),
        ('N#N', 'N2(A3Su+)', 'A3Su+', -1),
        ('[Ar]', 'Ar(1s5)', '1s5', -1),
    ], 1):
        mol = Molecule(smiles=smiles)
        mol.electronic_state = electronic
        mol.vibrational_level = vibrational
        thermo = NASA(polynomials=[
            NASAPolynomial(coeffs=[3.5, 0, 0, 0, 0, index * 1000, 2], Tmin=(200, 'K'), Tmax=(1000, 'K')),
            NASAPolynomial(coeffs=[3.5, 0, 0, 0, 0, index * 1000, 2], Tmin=(1000, 'K'), Tmax=(3000, 'K')),
        ], Tmin=(200, 'K'), Tmax=(3000, 'K'))
        species.append(Species(index=index, label=label, molecule=[mol], thermo=thermo))
    return species


class TestExcitedExport:
    def test_chemkin_state_labels_and_reimport(self, state_species, tmp_path):
        spcs = state_species
        assert [get_species_identifier(s) for s in spcs] == [
            'N2(1)', 'N2(v1)(2)', 'N2(eA3Su_p)(3)', 'Ar(e1s5)(4)',
        ]
        assert all(len(get_species_identifier(s)) <= 16 for s in spcs)
        assert all(re.fullmatch(r'[A-Za-z0-9\-_,()*#.:\[\]]+', get_species_identifier(s)) for s in spcs)
        rxn = Reaction(reactants=[spcs[0], spcs[2]], products=[spcs[1], spcs[1]],
                       kinetics=Arrhenius(A=(1e6, 'm^3/(mol*s)'), n=0, Ea=(0, 'J/mol')))
        chemkin = str(tmp_path / 'chem.inp')
        dictionary = str(tmp_path / 'species_dictionary.txt')
        save_chemkin_file(chemkin, spcs, [rxn])
        save_species_dictionary(dictionary, spcs)
        loaded, reactions = load_chemkin_file(chemkin, dictionary, use_chemkin_names=True)
        assert len(loaded) == 4
        assert [s.index for s in loaded] == [1, 2, 3, 4]
        assert [s.label for s in loaded] == ['N2', 'N2(v1)', 'N2(eA3Su_p)', 'Ar(e1s5)']
        assert [(s.molecule[0].electronic_state, s.molecule[0].vibrational_level) for s in loaded] == [
            ('', -1), ('', 1), ('A3Su+', -1), ('1s5', -1),
        ]
        assert all(a.is_isomorphic(b) for a, b in zip(spcs, loaded))
        assert all(not a.is_isomorphic(b) for a, b in itertools.combinations(loaded, 2))
        assert all(s.reactive for s in loaded[1:])
        assert len(reactions) == 1
        assert reactions[0].is_isomorphic(rxn)
        aliased, _ = load_chemkin_file(chemkin, dictionary)
        assert [s.label for s in aliased] == [s.label for s in spcs]
        assert all(a.is_isomorphic(b) for a, b in zip(spcs, aliased))

    @pytest.mark.parametrize('exporter', ['save_chemkin_file', 'save_chemkin_surface_file',
                                         'save_species_dictionary', 'save_transport_file'])
    def test_identifier_collision_refused_before_writing(self, state_species, tmp_path, exporter):
        import rmgpy.chemkin as chemkin
        from rmgpy.exceptions import ChemkinError
        first = Species(label='ethanol+', molecule=[Molecule(smiles='CCO')], thermo=state_species[0].thermo)
        second = Species(label='ether+', molecule=[Molecule(smiles='COC')], thermo=state_species[0].thermo)
        from rmgpy.transport import TransportData
        for species in (first, second):
            species.transport_data = TransportData(shapeIndex=2, sigma=(3.5, 'angstrom'),
                epsilon=(100, 'K'), dipoleMoment=(0, 'De'), polarizability=(0, 'angstrom^3'),
                rotrelaxcollnum=1)
        path = tmp_path / 'existing-file'
        path.write_text('sentinel')
        args = (str(path), [first, second])
        if exporter in ('save_chemkin_file', 'save_chemkin_surface_file'):
            args += ([],)
        with pytest.raises(ChemkinError, match='C2H6O.*ethanol.*ether') as error:
            getattr(chemkin, exporter)(*args)
        assert type(error.value).__name__ == 'ChemkinIdentifierCollisionError'
        assert path.read_text() == 'sentinel'

    def test_identifier_validator_accepts_state_convention(self, state_species):
        from rmgpy.chemkin import validate_species_identifiers
        validate_species_identifiers(state_species)

    @pytest.mark.parametrize('electronic,vibrational', [('X' * 32, -1), ('', 2147483647)])
    def test_identifier_refuses_overlong_state(self, electronic, vibrational):
        from rmgpy.exceptions import ChemkinError
        mol = Molecule(smiles='N#N')
        mol.electronic_state = electronic
        mol.vibrational_level = vibrational
        species = Species(label='N2', index=1, molecule=[mol])
        with pytest.raises(ChemkinError, match='16-character limit'):
            get_species_identifier(species)

    def test_dictionary_resolved_nitrogen_and_argon_are_reactive(self, state_species, tmp_path):
        from rmgpy.chemkin import load_species_dictionary
        path = str(tmp_path / 'species_dictionary.txt')
        save_species_dictionary(path, state_species)
        species = load_species_dictionary(path)
        assert not species[get_species_identifier(state_species[0])].reactive
        assert all(species[get_species_identifier(s)].reactive for s in state_species[1:])

    @pytest.mark.parametrize('labels', [['N2', 'N2', 'N2-2', 'Ar'], ['N2-2', 'N2', 'N2', 'Ar']])
    def test_rms_state_roundtrip_and_collider_names(self, state_species, tmp_path, labels):
        import yaml
        from rmgpy.kinetics import ThirdBody
        from rmgpy import yaml_rms
        # The pre-existing label N2-2 must not collide with an allocated alias.
        for spc, label in zip(state_species, labels):
            spc.label = label
        kinetics = ThirdBody(arrheniusLow=Arrhenius(A=(1e6, 'm^6/(mol^2*s)'), n=0, Ea=(0, 'J/mol')),
                             efficiencies={state_species[1].molecule[0]: 4})
        reaction = Reaction(reactants=[state_species[0], state_species[2]],
                            products=[state_species[1], state_species[1]], kinetics=kinetics)
        path = str(tmp_path / 'chem.rms')
        yaml_rms.write_rms(state_species, [reaction], path=path)
        data = yaml.safe_load(open(path))
        records = data['Phases'][0]['Species']
        assert len(records) == 4
        assert len({record['name'] for record in records}) == 4
        assert records[0]['smiles'] == 'N#N'
        assert 'adjlist' not in records[0]
        for record in records[1:]:
            assert 'adjlist' in record
            assert 'smiles' not in record
        assert 'vibrationallevel 1' in records[1]['adjlist']
        assert 'electronicstate A3Su+' in records[2]['adjlist']
        assert 'electronicstate 1s5' in records[3]['adjlist']
        assert data['Reactions'][0]['kinetics']['efficiencies'] == {records[1]['name']: 4.0}
        loaded = yaml_rms.load_rms_species(path)
        assert len(loaded) == 4
        assert [spc.label for spc in loaded] == [record['name'] for record in records]
        assert all(a.is_isomorphic(b) for a, b in zip(state_species, loaded))
        assert all(not a.is_isomorphic(b) for a, b in itertools.combinations(loaded, 2))
        assert kinetics.get_effective_collider_efficiencies(state_species).tolist() == [1, 4, 1, 1]

    def test_rms_reader_prefers_state_over_smiles_and_keeps_emitted_name(self, state_species, tmp_path):
        import yaml
        from rmgpy.yaml_rms import load_rms_species
        resolved = state_species[1]
        record = {'name': 'emitted-name', 'smiles': 'N#N',
                  'adjlist': resolved.molecule[0].to_adjacency_list(label='dictionary-name')}
        path = tmp_path / 'external.rms'
        path.write_text(yaml.safe_dump({'Phases': [{'Species': [record]}]}))
        loaded = load_rms_species(str(path))
        assert len(loaded) == 1
        assert loaded[0].label == 'emitted-name'
        assert loaded[0].is_isomorphic(resolved)
        assert not loaded[0].is_isomorphic(state_species[0])

    @pytest.mark.parametrize('writer', ['yaml_cantera1', 'yaml_cantera2'])
    def test_cantera_resolved_notes_and_unresolved_note_compatibility(self, state_species, writer):
        from importlib import import_module
        convert = import_module('rmgpy.' + writer).species_to_dict
        ground = convert(state_species[0], state_species)
        # At the baseline the attempted Species.to_smiles() is caught, so the
        # ground species has no note. Preserve that exact behavior.
        assert 'note' not in ground
        for spc, expected in zip(state_species[1:], [
            'vibrationallevel 1', 'electronicstate A3Su+', 'electronicstate 1s5',
        ]):
            note = convert(spc, state_species).get('note', '')
            assert expected in note
            parsed = Molecule().from_adjacency_list(note)
            assert parsed.is_isomorphic(spc.molecule[0])

    @pytest.mark.parametrize('old_style', [False, True])
    def test_smiles_keyed_library_refuses_resolved_colliders(self, state_species, tmp_path, old_style):
        import io
        from rmgpy.data.base import Entry
        from rmgpy.data.kinetics.common import save_entry
        from rmgpy.data.kinetics.library import KineticsLibrary
        from rmgpy.kinetics import ThirdBody
        ground, vibrational = state_species[:2]
        kinetics = ThirdBody(arrheniusLow=Arrhenius(A=(1e6, 'm^3/(mol*s)'), n=0, Ea=(0, 'J/mol')),
                             efficiencies={ground.molecule[0]: 2, vibrational.molecule[0]: 4})
        entry = Entry(index=1, label='N2 <=> N2', item=Reaction(reactants=[ground], products=[ground]), data=kinetics)
        stream = io.StringIO()
        if old_style:
            library = KineticsLibrary(label='excited_colliders')
            library.entries[1] = entry
            action = lambda: library.save_old(str(tmp_path / 'library'))
        else:
            action = lambda: save_entry(stream, entry)
        with pytest.raises(ValueError, match='resolved collider') as error:
            action()
        assert type(error.value).__name__ == 'SpeciesIdentityError'
        assert stream.getvalue() == ''
        assert not (tmp_path / 'library').exists()
        assert kinetics.get_effective_collider_efficiencies(state_species).tolist() == [2, 4, 1, 1]

    def test_condition_display_preserves_state_distinct_fractions(self, state_species):
        import ast
        from rmgpy.tools.canteramodel import CanteraCondition
        condition = CanteraCondition('IdealGasReactor', (1, 's'),
            dict(zip(state_species, [0.1, 0.2, 0.3, 0.4])), T0=(300, 'K'), P0=(1, 'bar'))
        text = str(condition).split('Initial Mole Fractions: ', 1)[1].splitlines()[0]
        fractions = ast.literal_eval(text)
        assert len(fractions) == 4
        assert sorted(fractions.values()) == [0.1, 0.2, 0.3, 0.4]
        assert fractions['N#N'] == 0.1
        assert any('/v:1' in key for key in fractions)
        assert any('/es:A3Su+' in key for key in fractions)

    def test_library_dictionary_refuses_shared_label_for_distinct_states(self, state_species, tmp_path):
        from rmgpy.data.base import Database, Entry
        ground, resolved = state_species[:2]
        ground.label = resolved.label = 'N#N'
        database = Database()
        database.entries[1] = Entry(index=1, item=Reaction(reactants=[ground], products=[resolved]))
        path = tmp_path / 'dictionary.txt'
        path.write_text('sentinel')
        with pytest.raises(ValueError, match='shared label.*N#N') as error:
            database.save_dictionary(str(path))
        assert type(error.value).__name__ == 'SpeciesIdentityError'
        assert path.read_text() == 'sentinel'

    @pytest.mark.parametrize('creation', ['species', 'model'])
    def test_unlabelled_species_names_retain_state(self, state_species, creation):
        molecules = [spc.molecule[0].copy(deep=True) for spc in state_species]
        if creation == 'species':
            species = [Species(molecule=[mol]) for mol in molecules]
            for spc in species:
                str(spc)
        else:
            from rmgpy.rmg.model import CoreEdgeReactionModel
            model = CoreEdgeReactionModel()
            species = []
            for mol in molecules:
                spc, is_new = model.make_new_species(mol, generate_thermo=False)
                assert is_new
                species.append(spc)
            for mol, spc in zip(molecules, species):
                assert model.check_for_existing_species(mol) is spc
        assert [spc.label for spc in species] == [
            'N#N', 'N#N|v:1', 'N#N|es:A3Su+', '[Ar]|es:1s5',
        ]
