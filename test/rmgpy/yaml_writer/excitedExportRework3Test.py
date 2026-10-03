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
import io
import importlib.util
import types
from pathlib import Path
from unittest.mock import patch

import pytest

from rmgpy.chemkin import load_chemkin_file, save_chemkin
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.common import save_entry
from rmgpy.data.kinetics.depository import KineticsDepository
from rmgpy.data.kinetics.family import KineticsFamily
from rmgpy.data.kinetics.library import KineticsLibrary
from rmgpy.exceptions import GenerationMismatchError, ResolvedStateTrainingError, SpeciesIdentityError
from rmgpy.kinetics import Arrhenius, ThirdBody
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.species import Species
from rmgpy.statmech import Conformer


# This helper is outside the installed packages; load its sibling file directly.
_helpers_spec = importlib.util.spec_from_file_location(
    'excited_export_helpers', Path(__file__).with_name('excitedExportHelpers.py'),
)
_helpers = importlib.util.module_from_spec(_helpers_spec)
_helpers_spec.loader.exec_module(_helpers)
nitrogen = _helpers.nitrogen
database = _helpers.database


class TestExcitedExportRework3Review(_helpers.FluxDiagramWitness):
    def test_review_1_depository_rejects_mixed_generation(self, tmp_path):
        old = tmp_path / 'old'
        new = tmp_path / 'new'
        old.mkdir()
        new.mkdir()
        database(KineticsLibrary, -1, 1).save(str(old / 'reactions.py'))
        database(KineticsLibrary, 1, 2).save(str(new / 'reactions.py'))
        (new / 'dictionary.txt').write_bytes((old / 'dictionary.txt').read_bytes())
        with pytest.raises(GenerationMismatchError, match='dictionary.txt'):
            KineticsDepository().load(str(new / 'reactions.py'), local_context={'Arrhenius': Arrhenius})

    def test_review_2_failed_training_save_preserves_old_pair(self, tmp_path):
        family = KineticsFamily()
        family.save_depository(database(KineticsDepository, 1, 1), str(tmp_path))
        replacement = database(KineticsDepository, -1, 2)
        with patch.object(replacement, 'save', side_effect=OSError('injected before reactions')):
            with pytest.raises(OSError, match='injected'):
                family.save_depository(replacement, str(tmp_path))
        loaded = KineticsDepository()
        loaded.load(str(tmp_path / 'reactions.py'), local_context={'Arrhenius': Arrhenius})
        entry = next(iter(loaded.entries.values()))
        assert entry.data.A.value_si == 1
        assert entry.item.products[0].molecule[0].vibrational_level == 1

    def test_review_3_arkane_yaml_prefers_adjacency(self, tmp_path):
        from arkane.common import ArkaneSpecies
        from arkane.statmech import StatMechJob
        spc = nitrogen('excited_N2', 1)
        spc.molecule[0].electronic_state = 'A'
        ArkaneSpecies(species=spc).save_yaml(str(tmp_path))
        path = next((tmp_path / 'species').glob('*.yml'))
        loaded = ArkaneSpecies(conformer=Conformer())
        loaded.load_yaml(str(path))
        mol = Molecule().from_adjacency_list(loaded.adjacency_list)
        assert mol.state_suffix() == spc.molecule[0].state_suffix()
        consumer = nitrogen('excited_N2')
        StatMechJob(consumer, str(path)).load()
        assert consumer.molecule[0].state_suffix() == spc.molecule[0].state_suffix()

    def test_review_4_arkane_chemkin_refuses_resolved(self, tmp_path):
        from arkane.thermo import ThermoJob
        with pytest.raises(SpeciesIdentityError, match='Arkane.*Chemkin'):
            ThermoJob(nitrogen(level=1), 'NASA').write_chemkin(str(tmp_path))
        assert not (tmp_path / 'chem.inp').exists()

    @pytest.mark.parametrize('case', ['participants', 'efficiency', 'thermo-coverage', 'kinetic-coverage'])
    def test_review_5_native_cantera_preserves_or_refuses(self, case):
        ground, excited = nitrogen(), nitrogen(level=1)
        if case == 'participants':
            rxn = Reaction(reactants=[ground], products=[excited], kinetics=Arrhenius(A=(1, 's^-1')))
            with pytest.raises(SpeciesIdentityError, match='Reaction.to_cantera'):
                rxn.to_cantera(species_list=[ground, excited])
        elif case == 'efficiency':
            mol = Molecule(smiles='[Ar]'); mol.vibrational_level = 1
            rate = ThirdBody(arrheniusLow=Arrhenius(A=(1, 'm^3/(mol*s)')), efficiencies={mol: 4})
            rxn = Reaction(reactants=[ground], products=[ground], kinetics=rate)
            with pytest.raises(SpeciesIdentityError, match='undeclared'):
                rxn.to_cantera(species_list=[ground], use_chemkin_identifier=True)
        elif case == 'thermo-coverage':
            ground.molecule = [Molecule().from_adjacency_list('1 X u0 p0 c0')]
            mol = copy.deepcopy(ground.molecule[0]); mol.vibrational_level = 1
            ground.thermo.thermo_coverage_dependence = {mol.to_adjacency_list(): {
                'model': 'polynomial', 'enthalpy-coefficients': [(1, 'J/mol')],
                'entropy-coefficients': [(1, 'J/(mol*K)')]}}
            with pytest.raises(SpeciesIdentityError, match='undeclared'):
                ground.to_cantera(all_species=[ground], use_chemkin_identifier=True)
        else:
            from rmgpy.kinetics import SurfaceArrhenius
            rate = SurfaceArrhenius(A=(1, 's^-1'), coverage_dependence={excited: {
                'a': 1, 'm': 0, 'E': (1, 'J/mol')}})
            rxn = Reaction(reactants=[ground], products=[ground], kinetics=rate)
            with pytest.raises(SpeciesIdentityError, match='undeclared'):
                rxn.to_cantera(species_list=[ground], use_chemkin_identifier=True)

    def test_review_6_stamped_chemkin_requires_dictionary(self, tmp_path):
        model = CoreEdgeReactionModel()
        model.core.species = [nitrogen('A'), nitrogen('B', 1)]
        save_chemkin(model, str(tmp_path / 'chem.inp'), str(tmp_path / 'annotated.inp'), str(tmp_path / 'dictionary.txt'))
        with pytest.raises(GenerationMismatchError, match='missing.*dictionary'):
            load_chemkin_file(str(tmp_path / 'chem.inp'))

    def test_review_8_arkane_library_pair_reload(self, tmp_path):
        from arkane.output import save_kinetics_lib
        a, b = nitrogen('A'), Species(label='B', smiles='[N]=[N]')
        rxn = Reaction(reactants=[a], products=[b], kinetics=Arrhenius(A=(1, 's^-1')))
        path = tmp_path / 'new'
        save_kinetics_lib([rxn], str(path), 'audit', '')
        loaded = KineticsLibrary()
        loaded.load(str(path / 'reactions.py'), local_context={'Arrhenius': Arrhenius})
        assert next(iter(loaded.entries.values())).data.A.value_si == 1

    @pytest.mark.parametrize('indexed', [False, True])
    def test_review_9_valid_ground_label_bytes(self, indexed):
        import json
        fixture = json.loads((Path(__file__).resolve().parents[1] /
            'test_data/excited_export_ground_spacing.json').read_text())
        entry = database(KineticsLibrary, -1, 1).entries[1]
        entry.label = 'A  <=>  B'
        if indexed:
            entry.item.reactants[0].index = 1
            entry.item.products[0].index = 2
            entry.label = 'A(1)  <=>  B(2)'
        output = io.StringIO()
        save_entry(output, copy.deepcopy(entry))
        assert fixture['payloads'][str(indexed)] == output.getvalue()


class TestExcitedExportCensusRegressions:
    @pytest.mark.parametrize('reference', ['condition', 'sensitivity'])
    def test_cantera_rejects_undeclared_state_alias(self, tmp_path, reference):
        from rmgpy.tools.canteramodel import Cantera, CanteraCondition
        ground, excited = nitrogen('N2(v1)'), nitrogen(level=1)
        condition = CanteraCondition('IdealGasReactor', (1e-7, 's'),
            {excited if reference == 'condition' else ground: 1}, T0=(500, 'K'), P0=(1, 'bar'))
        simulation = Cantera(species_list=[ground], reaction_list=[], output_directory=str(tmp_path),
            conditions=[condition], sensitive_species=[excited] if reference == 'sensitivity' else [])
        simulation.load_model()
        with pytest.raises(SpeciesIdentityError, match='undeclared'):
            simulation.simulate()

    @pytest.mark.parametrize('legacy', [False, True])
    def test_arkane_archive_refuses_identity_collision(self, tmp_path, legacy):
        import yaml
        from arkane.common import ArkaneSpecies
        old, new = nitrogen('N2(v1)'), nitrogen(level=1)
        ArkaneSpecies(species=old).save_yaml(str(tmp_path))
        path = next((tmp_path / 'species').glob('*.yml'))
        if legacy:
            data = yaml.safe_load(path.read_text())
            del data['adjacency_list']
            path.write_text(yaml.safe_dump(data))
        original = path.read_bytes()
        with pytest.raises(SpeciesIdentityError, match='Arkane YAML archive'):
            ArkaneSpecies(species=new).save_yaml(str(tmp_path))
        assert path.read_bytes() == original

    def test_profile_csv_refuses_duplicate_identity_name(self, tmp_path):
        from rmgpy.rmg.listener import SimulationProfileWriter
        (tmp_path / 'solver').mkdir()
        writer = SimulationProfileWriter(str(tmp_path), 0, [nitrogen('N2(v1)'), nitrogen(level=1)])
        with pytest.raises(SpeciesIdentityError, match='simulation profile'):
            writer.update(types.SimpleNamespace(snapshots=[]))
        assert not list((tmp_path / 'solver').iterdir())

    @pytest.mark.parametrize('nested', [False, True])
    def test_rms_direct_rate_refuses_resolved_coverage(self, nested):
        from rmgpy.kinetics import SurfaceArrhenius, Lindemann
        from rmgpy.yaml_rms import obj_to_dict
        excited = nitrogen(level=1)
        rate = SurfaceArrhenius(A=(1, 's^-1'), Ea=(0, 'J/mol'),
            coverage_dependence={excited: {'a': 1, 'm': 0, 'E': (0, 'J/mol')}})
        if nested:
            rate = Lindemann(arrheniusLow=rate, arrheniusHigh=Arrhenius(A=(1, 's^-1')))
        with pytest.raises(SpeciesIdentityError, match='RMS kinetics coverage'):
            obj_to_dict(rate, [excited])

    def test_thermo_library_refuses_resolved_coverage(self):
        from rmgpy.data.thermo import save_entry as save_thermo_entry
        ground = nitrogen()
        resolved = Molecule().from_adjacency_list('1 X u0 p0 c0')
        resolved.vibrational_level = 1
        ground.thermo.thermo_coverage_dependence = {resolved.to_adjacency_list(): {'model': 'polynomial',
            'enthalpy-coefficients': [(1, 'J/mol')], 'entropy-coefficients': [(1, 'J/(mol*K)')]}}
        output = io.StringIO()
        with pytest.raises(SpeciesIdentityError, match='thermo.save_entry coverage'):
            save_thermo_entry(output, Entry(index=1, label='N2', item=ground.molecule[0], data=ground.thermo))
        assert not output.getvalue()

    @pytest.mark.parametrize('loader', [KineticsLibrary, KineticsDepository])
    def test_indexed_spacing_label_reloads_without_stripping_state(self, tmp_path, loader):
        db = database(loader, 1, 1)
        entry = db.entries[1]
        entry.item.reactants[0].index = 1
        entry.item.products[0].index = 2
        entry.label = 'A(1)  <=>  B(2)'
        db.save(str(tmp_path / 'reactions.py'))
        loaded = loader()
        loaded.load(str(tmp_path / 'reactions.py'), local_context={'Arrhenius': Arrhenius})
        restored = next(iter(loaded.entries.values()))
        # Ground display indices are not dictionary aliases. Mixed exports
        # emit exact dictionary names while retaining the resolved identity.
        assert restored.label == 'A <=> B'
        assert entry.label == 'A(1)  <=>  B(2)'
        assert restored.item.reactants[0].is_isomorphic(entry.item.reactants[0])
        assert restored.item.products[0].is_isomorphic(entry.item.products[0])
        assert restored.item.products[0].molecule[0].has_resolved_state()

    def test_training_append_refuses_resolved_input_before_mutating_stamped_pair(self, tmp_path):
        from rmgpy import settings
        family = KineticsFamily(label='probe')
        training = tmp_path / 'kinetics/families/probe/training'
        training.mkdir(parents=True)
        initial = database(KineticsDepository, -1, 1)
        initial.label = 'probe/training'
        family.depositories = [initial]
        family.save_depository(initial, str(training))
        before = {path.name: path.read_bytes() for path in training.iterdir()}
        reaction = Reaction(reactants=[nitrogen('A')], products=[nitrogen('B', 1)],
            kinetics=Arrhenius(A=(2, 's^-1'), Ea=(0, 'J/mol')))
        with patch.dict(settings, {'database.directory': str(tmp_path)}):
            with pytest.raises(ResolvedStateTrainingError, match='do not support resolved species'):
                family.save_training_reactions([reaction])
        assert {path.name: path.read_bytes() for path in training.iterdir()} == before
        loaded = KineticsDepository()
        loaded.load(str(training / 'reactions.py'), local_context={'Arrhenius': Arrhenius})
        assert len(loaded.entries) == 1


    def test_native_ground_participants_outside_phase_inventory(self):
        ground, product = nitrogen('A'), nitrogen('B')
        rxn = Reaction(reactants=[ground], products=[product], kinetics=Arrhenius(A=(1, 's^-1')))
        assert rxn.to_cantera(species_list=[], use_chemkin_identifier=True).equation == 'A <=> B'
        product.molecule[0].vibrational_level = 1
        with pytest.raises(SpeciesIdentityError, match='undeclared'):
            rxn.to_cantera(species_list=[], use_chemkin_identifier=True)

    def test_model_free_ground_translation(self, tmp_path):
        import cantera as ct
        from rmgpy.rmg.main import RMG
        model = CoreEdgeReactionModel()
        model.core.species = [nitrogen()]
        save_chemkin(model, str(tmp_path / 'chem.inp'), str(tmp_path / 'annotated.inp'),
                            str(tmp_path / 'species_dictionary.txt'))
        rmg = RMG(output_directory=str(tmp_path))
        assert rmg.reaction_model is None
        (tmp_path / 'tran.dat').write_text('N2 1 100.0 3.5 0.0 1.0 1.0\n')
        path = rmg.generate_cantera_files_from_chemkin(str(tmp_path / 'chem.inp'))
        assert ct.Solution(path).species_names == ['N2']

    @pytest.mark.parametrize('legacy', [False, True])
    def test_model_free_translation_refuses_resolved_pair(self, tmp_path, legacy):
        from rmgpy.rmg.main import RMG
        model = CoreEdgeReactionModel()
        model.core.species = [nitrogen(level=1)]
        paths = [tmp_path / 'chem.inp', tmp_path / 'species_dictionary.txt']
        save_chemkin(model, str(paths[0]), str(tmp_path / 'annotated.inp'), str(paths[1]))
        if legacy:
            from rmgpy.util import strip_generation_marker
            for path in paths:
                path.write_text(strip_generation_marker(path.read_text()))
        rmg = RMG(output_directory=str(tmp_path))
        with pytest.raises(SpeciesIdentityError, match='RMG.generate_cantera_files_from_chemkin'):
            rmg.generate_cantera_files_from_chemkin(str(paths[0]))
        assert not (tmp_path / 'cantera_from_ck').exists()

    def test_model_free_stamped_translation_requires_dictionary(self, tmp_path):
        from rmgpy.rmg.main import RMG
        model = CoreEdgeReactionModel()
        model.core.species = [nitrogen()]
        dictionary = tmp_path / 'species_dictionary.txt'
        save_chemkin(model, str(tmp_path / 'chem.inp'), str(tmp_path / 'annotated.inp'), str(dictionary))
        dictionary.unlink()
        rmg = RMG(output_directory=str(tmp_path))
        with pytest.raises(GenerationMismatchError, match='paired species dictionary'):
            rmg.generate_cantera_files_from_chemkin(str(tmp_path / 'chem.inp'))
        assert not (tmp_path / 'cantera_from_ck').exists()

    @pytest.mark.parametrize('missing_dictionary', [False, True])
    def test_translation_validates_surface_generation(self, tmp_path, missing_dictionary):
        from rmgpy.rmg.main import RMG
        from rmgpy.util import stamp_file_generation, strip_generation_marker
        model = CoreEdgeReactionModel()
        model.core.species = [nitrogen()]
        gas = tmp_path / 'chem.inp'
        surface = tmp_path / 'surface.inp'
        dictionary = tmp_path / 'species_dictionary.txt'
        save_chemkin(model, str(gas), str(surface), str(dictionary))
        if missing_dictionary:
            dictionary.unlink()
            gas.write_text(strip_generation_marker(gas.read_text()))
        else:
            content = strip_generation_marker(surface.read_text()) + '\n! separate generation\n'
            surface.write_text(stamp_file_generation([(str(surface), content, '!')])[0][1])
        rmg = RMG(output_directory=str(tmp_path))
        with pytest.raises(GenerationMismatchError, match='surface.inp'):
            rmg.generate_cantera_files_from_chemkin(str(gas), surface_file=str(surface))
        assert not (tmp_path / 'cantera_from_ck').exists()

    @pytest.mark.parametrize('surface', [False, True])
    def test_translation_consumes_validated_snapshots(self, tmp_path, monkeypatch, surface):
        from types import SimpleNamespace
        import rmgpy.chemkin as chemkin_module
        import rmgpy.rmg.main as main
        model = CoreEdgeReactionModel()
        model.core.species = [nitrogen()]
        gas = tmp_path / 'chem.inp'
        other = tmp_path / 'surface.inp'
        dictionary = tmp_path / 'species_dictionary.txt'
        save_chemkin(model, str(gas), str(other), str(dictionary))
        paths = [gas, other] if surface else [gas]
        original = {str(path): path.read_text() for path in paths}
        original_loader = chemkin_module.load_species_dictionary

        def swap_after_validation(*args, **kwargs):
            declarations = original_loader(*args, **kwargs)
            model.core.species = [nitrogen(level=1)]
            save_chemkin(model, str(gas), str(other), str(dictionary))
            return declarations

        seen = {}

        def translate(input_file, **kwargs):
            staged = [input_file] + ([kwargs['surface_file']] if surface else [])
            for path in staged:
                assert Path(path).name in [path.name for path in paths]
                seen[Path(path).name] = Path(path).read_text()

        monkeypatch.setattr(chemkin_module, 'load_species_dictionary', swap_after_validation)
        monkeypatch.setattr(main.ck2yaml, 'Parser', lambda: SimpleNamespace(convert_mech=translate))
        kwargs = {'surface_file': str(other)} if surface else {}
        main.RMG(output_directory=str(tmp_path)).generate_cantera_files_from_chemkin(str(gas), **kwargs)
        assert seen == {Path(path).name: content for path, content in original.items()}
        assert 'N2(v1)' in gas.read_text()

    def test_ground_arkane_drawing_uses_real_wells(self, tmp_path):
        from arkane.kinetics import KineticsDrawer, Well
        from rmgpy.species import TransitionState
        species = [nitrogen('A'), nitrogen('B')]
        for spc in species:
            spc.conformer.E0 = (0, 'kJ/mol')
        reaction = Reaction(reactants=species[:1], products=species[1:],
            transition_state=TransitionState(conformer=Conformer(E0=(10, 'kJ/mol'))))
        drawer = KineticsDrawer(options={'structures': False})
        assert drawer._get_label_size(Well(species[:1]), file_format='png')[2] > 0
        path = tmp_path / 'reaction.png'
        drawer.draw(reaction, file_format='png', path=str(path))
        assert path.read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
        assert len(drawer.wells) == 2
