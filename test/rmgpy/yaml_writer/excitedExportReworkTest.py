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

import copy
import os

import pytest
import yaml

from rmgpy.chemkin import (get_species_identifier, load_chemkin_file,
                          load_species_dictionary, save_chemkin_file,
                          save_species_dictionary)
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.library import KineticsLibrary
from rmgpy.data.kinetics.depository import KineticsDepository
from rmgpy.exceptions import ChemkinError, SpeciesIdentityError
from rmgpy.export import SpeciesReferences
from rmgpy.kinetics import Arrhenius, Lindemann, ThirdBody, Troe
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
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


def falloff(kind=Lindemann, efficiencies=None):
    low = Arrhenius(A=(1e6, 'm^3/(mol*s)'), n=0, Ea=(0, 'J/mol'))
    if kind is ThirdBody:
        return kind(arrheniusLow=low, efficiencies=efficiencies or {})
    high = Arrhenius(A=(1e3, 's^-1'), n=0, Ea=(0, 'J/mol'))
    kwargs = {'alpha': 0.5, 'T3': (100, 'K'), 'T1': (1000, 'K')} if kind is Troe else {}
    return kind(arrheniusHigh=high, arrheniusLow=low, efficiencies=efficiencies or {}, **kwargs)


def library_with(ground, product, collider=None, kinetics=None):
    library = KineticsLibrary(label='state_export')
    reaction = Reaction(reactants=[ground], products=[product], specific_collider=collider)
    collider_label = ' (+{0})'.format(collider.label) if collider else ''
    label = ground.label + collider_label + ' <=> ' + product.label + collider_label
    library.entries[1] = Entry(index=1, label=label, item=reaction,
                               data=kinetics or Arrhenius(A=(1e3, 's^-1'), Ea=(0, 'J/mol')))
    return library


class TestExcitedExportRework:
    @pytest.mark.parametrize('label', ['N2', 'unique_collider', 'N2(v1)(2)', 'N2(eA,v1)(2)'])
    def test_named_library_collider_round_trip(self, state_species, tmp_path, label):
        ground, resolved, _, argon = state_species
        resolved.label = label
        library = library_with(argon, argon, resolved, falloff())
        path = tmp_path / 'reactions.py'
        library.save(str(path))
        loaded = KineticsLibrary()
        loaded.load(str(path), local_context={'Arrhenius': Arrhenius, 'Lindemann': Lindemann})
        collider = next(iter(loaded.entries.values())).item.specific_collider
        assert collider.is_isomorphic(resolved)
        assert not collider.is_isomorphic(ground)

    @pytest.mark.parametrize('label', ['M', 'm'])
    def test_reserved_named_collider_save_refused(self, state_species, tmp_path, label):
        _, collider, _, argon = state_species
        collider.label = label
        library = library_with(argon, argon, collider, falloff())
        path = tmp_path / 'reactions.py'
        dictionary = tmp_path / 'dictionary.txt'
        path.write_text('existing reactions')
        dictionary.write_text('existing dictionary')
        with pytest.raises(SpeciesIdentityError, match='reserved'):
            library.save(str(path))
        assert path.read_text() == 'existing reactions'
        assert dictionary.read_text() == 'existing dictionary'

    @pytest.mark.parametrize('label', ['M', 'm'])
    @pytest.mark.parametrize('loader', [KineticsLibrary, KineticsDepository])
    def test_reserved_dictionary_collider_load_refused(self, state_species, tmp_path, label, loader):
        _, collider, _, argon = state_species
        collider.label = label
        library = library_with(argon, argon, collider, falloff())
        entry = library.entries[1]
        path = tmp_path / 'reactions.py'
        path.write_text('entry(index=1, label={0!r}, kinetics={1!r})\n'.format(entry.label, entry.data))
        (tmp_path / 'dictionary.txt').write_text(
            argon.molecule[0].to_adjacency_list(label=argon.label) + '\n' +
            collider.molecule[0].to_adjacency_list(label=collider.label) + '\n')
        with pytest.raises(SpeciesIdentityError, match='reserved'):
            loader().load(str(path), local_context={'Arrhenius': Arrhenius, 'Lindemann': Lindemann})

    def test_named_collider_dictionary_collision_refused(self, state_species, tmp_path):
        ground, resolved = state_species[:2]
        resolved.label = ground.label
        library = library_with(ground, ground, resolved, falloff())
        path = tmp_path / 'dictionary.txt'
        path.write_text('sentinel')
        with pytest.raises(SpeciesIdentityError, match='shared label'):
            library.save_dictionary(str(path))
        assert path.read_text() == 'sentinel'

    @pytest.mark.parametrize('later_rename', [False, True])
    def test_thermo_library_name_route_refuses_unproven_resolved_energy_transfer(
            self, state_species, monkeypatch, later_rename):
        from rmgpy.exceptions import StateProvenanceError
        from rmgpy.rmg import model as model_module
        model = model_module.CoreEdgeReactionModel()
        templates = {spc.molecule[0].state_suffix(): spc.thermo for spc in state_species[:3]}
        def submit(spc, solvent):
            spc.thermo = copy.deepcopy(templates[spc.molecule[0].state_suffix()])
            spc.thermo.label = 'N2'
        monkeypatch.setattr(model_module, 'submit', submit)
        source = state_species[1]
        if later_rename:
            spc, new = model.make_new_species(source.molecule[0].copy(deep=True), generate_thermo=False)
            assert new
            with pytest.raises(StateProvenanceError, match='energy transfer.*resolved-state provenance'):
                model.generate_thermo(spc, rename=True)
        else:
            with pytest.raises(StateProvenanceError, match='energy transfer.*resolved-state provenance'):
                model.make_new_species(source.molecule[0].copy(deep=True), generate_thermo=True)

    @pytest.mark.parametrize('failure', ['dictionary_collision', 'resolved_efficiency'])
    def test_library_refusal_preserves_existing_pair(self, state_species, tmp_path, failure):
        ground, resolved = state_species[:2]
        path = tmp_path / 'reactions.py'
        library_with(ground, ground).save(str(path))
        before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
        if failure == 'dictionary_collision':
            resolved.label = ground.label
            library = library_with(ground, resolved, kinetics=Arrhenius(A=(2e3, 's^-1'), Ea=(0, 'J/mol')))
        else:
            library = library_with(ground, resolved, kinetics=falloff(efficiencies={resolved.molecule[0]: 4}))
        with pytest.raises(SpeciesIdentityError):
            library.save(str(path))
        assert {p.name: p.read_bytes() for p in tmp_path.iterdir()} == before
        loaded = KineticsLibrary()
        loaded.load(str(path), local_context={'Arrhenius': Arrhenius, 'Lindemann': Lindemann})
        assert next(iter(loaded.entries.values())).item.reactants[0].is_isomorphic(ground)
        assert next(iter(loaded.entries.values())).item.products[0].is_isomorphic(ground)

    @pytest.mark.parametrize('existing', [True, False])
    @pytest.mark.parametrize('failure_stage', ['staging', 'second_rename'])
    def test_library_io_failure_preserves_both_files(self, state_species, tmp_path, monkeypatch,
                                                    existing, failure_stage):
        import tempfile
        ground, resolved = state_species[:2]
        path = tmp_path / 'reactions.py'
        if existing:
            library_with(ground, ground).save(str(path))
        before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
        library = library_with(ground, resolved)
        if failure_stage == 'staging':
            original = tempfile.mkstemp
            calls = []
            def fail(*args, **kwargs):
                calls.append(1)
                if len(calls) == 2:
                    raise OSError('injected staging failure')
                return original(*args, **kwargs)
            monkeypatch.setattr(tempfile, 'mkstemp', fail)
        else:
            original = os.replace
            calls = []
            def fail(src, dst):
                calls.append(1)
                if len(calls) == 2:
                    raise OSError('injected second rename failure')
                return original(src, dst)
            monkeypatch.setattr(os, 'replace', fail)
        with pytest.raises(OSError, match='injected'):
            library.save(str(path))
        assert {p.name: p.read_bytes() for p in tmp_path.iterdir()} == before

    @pytest.mark.parametrize('reference', ['reactant', 'product', 'collider', 'efficiency', 'coverage'])
    @pytest.mark.parametrize('collision', [False, True])
    @pytest.mark.parametrize('writer', ['gas', 'surface'])
    def test_chemkin_refuses_undeclared_references(self, state_species, tmp_path, reference, collision, writer):
        ground, resolved = state_species[:2]
        resolved.index = -1
        if collision:
            ground.index = -1
            ground.label = get_species_identifier(resolved)
        reaction = Reaction(reactants=[ground], products=[ground], kinetics=Arrhenius(A=(1, 's^-1'), Ea=(0, 'J/mol')))
        if reference == 'reactant':
            reaction.reactants = [resolved]
        elif reference == 'product':
            reaction.products = [resolved]
        elif reference == 'collider':
            reaction.specific_collider = resolved
            reaction.kinetics = falloff()
        elif reference == 'efficiency':
            reaction.kinetics = falloff(efficiencies={resolved.molecule[0]: 4})
        else:
            from rmgpy.kinetics import SurfaceArrhenius
            from rmgpy.quantity import Quantity
            reaction.kinetics = SurfaceArrhenius(A=(1, 's^-1'), Ea=(0, 'J/mol'), coverage_dependence={
                resolved: {'a': Quantity(0), 'm': Quantity(0), 'E': Quantity(0, 'J/mol')}})
        path = tmp_path / 'chem.inp'
        path.write_text('sentinel')
        with pytest.raises(ChemkinError, match='undeclared|collides'):
            from rmgpy.chemkin import save_chemkin_surface_file
            exporter = save_chemkin_file if writer == 'gas' else save_chemkin_surface_file
            exporter(str(path), [ground], [reaction])
        assert path.read_text() == 'sentinel'

    @pytest.mark.parametrize('index', [-1, 2])
    @pytest.mark.parametrize('electronic,vibrational', [('', -1), ('A', -1), ('', 1), ('A', 1)])
    @pytest.mark.parametrize('collider_position', [False, True])
    def test_full_identifier_space_round_trip(self, state_species, tmp_path, index, electronic,
                                              vibrational, collider_position):
        target = state_species[0]
        target.index = index
        target.molecule[0].electronic_state = electronic
        target.molecule[0].vibrational_level = vibrational
        argon = state_species[3]
        argon.molecule[0].electronic_state = ''
        argon.label = 'Ar'
        if collider_position:
            reaction = Reaction(reactants=[argon], products=[argon], specific_collider=target,
                                kinetics=falloff())
        else:
            reaction = Reaction(reactants=[target], products=[argon], kinetics=Arrhenius(A=(1, 's^-1'), Ea=(0, 'J/mol')))
        chemkin, dictionary = tmp_path / 'chem.inp', tmp_path / 'dictionary.txt'
        save_chemkin_file(str(chemkin), [target, argon], [reaction])
        save_species_dictionary(str(dictionary), [target, argon])
        loaded, reactions = load_chemkin_file(str(chemkin), str(dictionary), use_chemkin_names=True)
        restored = reactions[0].specific_collider if collider_position else reactions[0].reactants[0]
        assert restored.is_isomorphic(target)
        assert restored.index == index
        if target.molecule[0].has_resolved_state():
            assert get_species_identifier(restored) == get_species_identifier(target)
        assert any(spc is restored for spc in loaded)

    @pytest.mark.parametrize('token', ['InChI', 'AInChI'])
    def test_electronic_inchi_token_dictionary_round_trip(self, state_species, tmp_path, token):
        spc = state_species[0]
        spc.molecule[0].electronic_state = token
        path = tmp_path / 'dictionary.txt'
        save_species_dictionary(str(path), [spc])
        loaded = load_species_dictionary(str(path))
        restored = loaded[get_species_identifier(spc)]
        assert restored.is_isomorphic(spc)
        assert restored.molecule[0].electronic_state == token

    @pytest.mark.parametrize('terminal_newline', [True, False])
    def test_dictionary_duplicate_identity_refused(self, state_species, tmp_path, terminal_newline):
        ground, resolved = state_species[:2]
        path = tmp_path / 'dictionary.txt'
        text = ground.molecule[0].to_adjacency_list(label='N2') + '\n' + resolved.molecule[0].to_adjacency_list(label='N2')
        path.write_text(text + ('\n' if terminal_newline else ''))
        with pytest.raises(ChemkinError, match='N2.*conflicting|conflicting.*N2'):
            load_species_dictionary(str(path))

    @pytest.mark.parametrize('record_termination', ['\n\n', ''])
    def test_dictionary_ground_duplicate_keeps_historical_last_record_bytes(self, tmp_path, record_termination):
        first = Species(label='bath', molecule=[Molecule(smiles='CCO')])
        second = Species(label='bath', molecule=[Molecule(smiles='COC')])
        second.molecule[0].vertices.reverse()
        expected = second.to_adjacency_list()
        path = tmp_path / 'dictionary.txt'
        text = first.to_adjacency_list().strip() + '\n\n' + second.to_adjacency_list().strip()
        path.write_text(text + record_termination)

        loaded = load_species_dictionary(str(path), generate_resonance_structures=False)

        assert list(loaded) == ['bath']
        assert loaded['bath'].to_adjacency_list() == expected

    def test_species_references_ground_collision_does_not_reorder_inputs(self):
        first = Species(label='bath', molecule=[Molecule(smiles='CCO')])
        second = Species(label='bath', molecule=[Molecule(smiles='COC')])
        second.molecule[0].vertices.reverse()
        before = [species.to_adjacency_list() for species in (first, second)]

        SpeciesReferences([first, second], ['bath', 'bath'])

        assert [species.to_adjacency_list() for species in (first, second)] == before

    def test_native_cantera_yaml_preserves_state(self, state_species, tmp_path):
        import cantera as ct
        species = [spc.to_cantera(use_chemkin_identifier=True) for spc in state_species]
        gas = ct.Solution(thermo='ideal-gas', kinetics='gas', species=species, reactions=[])
        path = tmp_path / 'native.yaml'
        gas.write_yaml(str(path))
        data = yaml.safe_load(path.read_text())
        assert 'note' not in data['species'][0]
        for source, record in zip(state_species[1:], data['species'][1:]):
            assert Molecule().from_adjacency_list(record['note']).is_isomorphic(source.molecule[0])
        restored = ct.Solution(str(path))
        for source, converted in zip(state_species[1:], restored.species()[1:]):
            assert Molecule().from_adjacency_list(converted.input_data['note']).is_isomorphic(source.molecule[0])

    @pytest.mark.parametrize('kind', [ThirdBody, Lindemann, Troe])
    def test_ground_rms_efficiency_uses_emitted_name(self, state_species, kind):
        from rmgpy.yaml_rms import get_mech_dict
        nitrogen, _, _, argon = state_species
        argon.molecule[0].electronic_state = ''
        nitrogen.label = argon.label = 'bath'
        reaction = Reaction(reactants=[nitrogen], products=[argon],
                            kinetics=falloff(kind, {nitrogen.molecule[0]: 4}))
        data = get_mech_dict([nitrogen, argon], [reaction])
        records = data['Phases'][0]['Species']
        assert records[0]['name'] == 'bath-2'
        assert records[1]['name'] == 'bath'
        assert data['Reactions'][0]['kinetics']['efficiencies'] == {'bath-2': 4.0}
        assert records[0]['smiles'] == 'N#N'
        assert records[1]['smiles'] == '[Ar]'


    def test_depository_refuses_resolved_named_collider(self, state_species, tmp_path):
        from rmgpy.data.kinetics.depository import KineticsDepository
        ground, collider, _, argon = state_species
        collider.label = 'N2'
        library = library_with(argon, argon, collider, falloff())
        path = tmp_path / 'reactions.py'
        library.save(str(path))
        # A complete dictionary isolates the loader from the old omission in saving.
        dictionary = ''.join(spc.molecule[0].to_adjacency_list(label=spc.label) + '\n'
                             for spc in [argon, collider])
        (tmp_path / 'dictionary.txt').write_text(dictionary)
        # This deliberately hand-written dictionary belongs to a token-free fixture.
        from rmgpy.util import strip_generation_marker
        path.write_text(strip_generation_marker(path.read_text()))
        with pytest.raises(SpeciesIdentityError, match='resolved named collider'):
            KineticsDepository().load(str(path), local_context={'Arrhenius': Arrhenius, 'Lindemann': Lindemann})

    def test_rms_refuses_resolved_named_collider(self, state_species, tmp_path):
        from rmgpy.yaml_rms import write_rms
        ground, collider = state_species[:2]
        reaction = Reaction(reactants=[ground], products=[ground], specific_collider=collider,
                            kinetics=falloff())
        path = tmp_path / 'chem.rms'
        path.write_text('sentinel')
        with pytest.raises(SpeciesIdentityError, match='resolved named collider'):
            write_rms(state_species, [reaction], path=str(path))
        assert path.read_text() == 'sentinel'

    def test_rms_admission_refuses_shared_label_for_distinct_states(self, state_species, monkeypatch):
        from rmgpy.rmg import reactionmechanismsimulator_reactors as reactors
        monkeypatch.setattr(reactors, 'to_rms', lambda spc: spc.label)
        ground, resolved = state_species[:2]
        ground.label = resolved.label = 'N2'
        phase = reactors.Phase()
        phase.add_species(ground)
        with pytest.raises(SpeciesIdentityError, match='shared RMS label'):
            phase.add_species(resolved)
        assert phase.rmg_species == [ground]

    @pytest.mark.parametrize('side', ['reactants', 'products'])
    def test_ambiguous_resolved_participant_save_refused(self, state_species, tmp_path, side):
        ground, resolved = state_species[:2]
        resolved.label = 'A+B'
        library = library_with(ground, ground)
        setattr(library.entries[1].item, side, [resolved])
        library.entries[1].label = ('A+B <=> N2' if side == 'reactants' else 'N2 <=> A+B')
        path = tmp_path / 'reactions.py'
        dictionary = tmp_path / 'dictionary.txt'
        path.write_text('existing reactions')
        dictionary.write_text('existing dictionary')
        with pytest.raises(SpeciesIdentityError, match='ambiguous reaction separator'):
            library.save(str(path))
        assert path.read_text() == 'existing reactions'
        assert dictionary.read_text() == 'existing dictionary'

    @pytest.mark.parametrize('loader', [KineticsLibrary, KineticsDepository])
    def test_ambiguous_resolved_participant_load_refused(self, state_species, tmp_path, loader):
        ground, resolved, _, argon = state_species
        ground.label, resolved.label, argon.label = 'A', 'A+B', 'B'
        argon.molecule[0].electronic_state = ''
        path = tmp_path / 'reactions.py'
        path.write_text('entry(index=1, label="A+B <=> A+B", kinetics={0!r})\n'.format(
            Arrhenius(A=(1, 's^-1'), Ea=(0, 'J/mol'))))
        (tmp_path / 'dictionary.txt').write_text(''.join(
            spc.molecule[0].to_adjacency_list(label=spc.label) + '\n'
            for spc in [ground, resolved, argon]))
        with pytest.raises(SpeciesIdentityError, match='ambiguous reaction separator'):
            loader().load(str(path), local_context={'Arrhenius': Arrhenius})
