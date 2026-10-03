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

from pathlib import Path
from types import SimpleNamespace

import pytest

import rmgpy.data.rmg as data_module
import rmgpy.rmg.input as inp
from rmgpy.data.thermo import ThermoDatabase, ThermoLibrary
from rmgpy.exceptions import SpeciesIdentityError, DatabaseError, StateProvenanceError
from rmgpy.rmg.main import RMG
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.species import Species

FIXTURES = Path(__file__).resolve().parents[1] / 'test_data' / 'excited_states'
N2 = '1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}'


@pytest.fixture
def database(monkeypatch):
    db = ThermoDatabase()
    lib = ThermoLibrary(label='ExcitedStateFixture')
    lib.load(str(FIXTURES / 'thermo.py'), db.local_context, {})
    db.libraries = {lib.label: lib}
    db.library_order = [lib.label]
    monkeypatch.setattr(data_module, 'database', SimpleNamespace(thermo=db, solvation=None))
    monkeypatch.setattr(inp, 'rmg', None)
    monkeypatch.setattr(inp, 'species_dict', {})
    return db


def load_deck(tmp_path, monkeypatch, declared=True, level=True):
    monkeypatch.setattr(data_module, 'database', None)
    deck = tmp_path / 'input.py'
    text = "database(thermoLibraries=[], kineticsFamilies='none')\n"
    text += "species(label='N2', structure=SMILES('N#N'))\n"
    if level:
        text += "species(label='N2v1', structure=adjacencyList(" + repr('vibrationallevel 1\n' + N2) + "))\n"
    if declared:
        text += "vibrationalManifold(species='N2')\n"
    text += 'simulator(atol=1e-16, rtol=1e-8)\nmodel(toleranceMoveToCore=0.1, toleranceInterruptSimulation=1.0)\n'
    deck.write_text(text)
    job = RMG()
    job.load_input(str(deck))
    return job


def install_database(monkeypatch, db):
    monkeypatch.setattr(data_module, 'database', SimpleNamespace(thermo=db, solvation=None))


def test_declared_unkeyed_species_gets_level_zero_thermo(database, tmp_path, monkeypatch, caplog):
    job = load_deck(tmp_path, monkeypatch)
    install_database(monkeypatch, database)
    ground = job.initial_species[0]
    with caplog.at_level('INFO'):
        job.reaction_model.generate_thermo(ground)
    assert ground.molecule[0].vibrational_level == -1
    assert ground.thermo.label == 'N2v0'
    assert ground.thermo.get_heat_capacity(1000) == pytest.approx(29.10161916, rel=0.0001)
    assert 'vibrationalManifold' in caplog.text
    assert 'N2' in caplog.text and 'v = 0' in caplog.text
    assert database.get_thermo_data(ground).label == 'N2v0'


def test_declared_species_refuses_another_thermo(database, tmp_path, monkeypatch):
    job = load_deck(tmp_path, monkeypatch)
    install_database(monkeypatch, database)
    ground = job.initial_species[0]
    ground.thermo = database.get_thermo_data(Species(label='thermal').from_smiles('N#N'))
    with pytest.raises(DatabaseError, match='N2') as caught:
        job.reaction_model.generate_thermo(ground)
    assert type(caught.value).__name__ == 'VibrationalManifoldError'


def test_missing_level_zero_library_entry_names_species(database, tmp_path, monkeypatch):
    job = load_deck(tmp_path, monkeypatch)
    install_database(monkeypatch, database)
    del database.libraries['ExcitedStateFixture'].entries['N2v0']
    with pytest.raises(DatabaseError, match='N2') as caught:
        job.reaction_model.generate_thermo(job.initial_species[0])
    assert type(caught.value).__name__ == 'VibrationalManifoldError'
    assert 'vibrationallevel 0' in str(caught.value)


def test_levels_without_declaration_are_refused(database, tmp_path, monkeypatch):
    with pytest.raises(DatabaseError, match='N2v1') as caught:
        load_deck(tmp_path, monkeypatch, declared=False)
    assert type(caught.value).__name__ == 'VibrationalManifoldError'


def test_model_without_levels_or_declaration_keeps_ensemble(database, tmp_path, monkeypatch):
    job = load_deck(tmp_path, monkeypatch, declared=False, level=False)
    install_database(monkeypatch, database)
    ground = job.initial_species[0]
    job.reaction_model.generate_thermo(ground)
    assert ground.thermo.label == 'N2'
    assert ground.thermo.get_heat_capacity(1000) == pytest.approx(37.41636749, rel=0.0001)


def test_generated_level_without_declaration_is_refused(database):
    model = CoreEdgeReactionModel()
    level = Species(label='N2v1').from_adjacency_list('vibrationallevel 1\n' + N2)
    with pytest.raises(DatabaseError, match='N2v1'):
        model.make_new_species(level, generate_thermo=False)


@pytest.mark.parametrize('label,state', [('Absent', ''), ('N2v1', 'vibrationallevel 1\n')])
def test_declaration_requires_existing_unresolved_input_species(database, label, state):
    inp.set_global_rmg(RMG())
    inp.rmg.reaction_model = CoreEdgeReactionModel()
    inp.rmg.initial_species = []
    if state:
        # Real input loading defers level validation until declarations are read.
        inp.rmg.reaction_model.defer_vibrational_validation = True
        inp.species(label, inp.adjacency_list(state + N2))
    with pytest.raises(DatabaseError, match=label) as caught:
        inp.vibrational_manifold(species=label)
    assert type(caught.value).__name__ == 'VibrationalManifoldError'


def test_declared_species_direct_estimation_is_refused(database, tmp_path, monkeypatch):
    job = load_deck(tmp_path, monkeypatch)
    ground = job.initial_species[0]
    for estimate in (lambda: database.get_thermo_data_from_groups(ground),
                     lambda: database.estimate_thermo_via_group_additivity(ground.molecule[0])):
        with pytest.raises(DatabaseError, match='N2') as caught:
            estimate()
        assert type(caught.value).__name__ == 'VibrationalManifoldError'


@pytest.mark.parametrize('state', ['vibrationallevel 1\n', 'electronicstate A3Su+\n'])
def test_unlabeled_resolved_species_label_is_distinct_from_ground(state):
    model = CoreEdgeReactionModel()
    ground, _ = model.make_new_species(Species().from_smiles('N#N'), generate_thermo=False)
    # The declaration is required only for vibrational levels.
    if state.startswith('vibrationallevel'):
        model.vibrational_manifolds = [ground]
    excited, _ = model.make_new_species(Species().from_adjacency_list(state + N2), generate_thermo=False)
    assert excited.label != ground.label
    assert excited.molecule[0].state_suffix() in excited.label


def test_declared_model_accepts_later_level_and_refuses_second_level_zero(database, tmp_path, monkeypatch):
    job = load_deck(tmp_path, monkeypatch, level=False)
    install_database(monkeypatch, database)
    level = Species(label='N2v1').from_adjacency_list('vibrationallevel 1\n' + N2)
    added, is_new = job.reaction_model.make_new_species(level, generate_thermo=False)
    with pytest.raises(StateProvenanceError, match="energy transfer"):
        job.reaction_model.generate_thermo(added)
    assert is_new and added.thermo.label == 'N2v1'
    zero = Species(label='ExplicitN2v0').from_adjacency_list('vibrationallevel 0\n' + N2)
    with pytest.raises(DatabaseError, match='ExplicitN2v0'):
        job.reaction_model.make_new_species(zero, generate_thermo=False)


def test_saved_input_preserves_declaration(database, tmp_path, monkeypatch):
    job = load_deck(tmp_path, monkeypatch, level=False)
    saved = tmp_path / 'saved.py'
    inp.save_input_file(str(saved), job)
    assert "vibrationalManifold(species='N2')" in saved.read_text()
    reloaded = RMG()
    reloaded.load_input(str(saved))
    install_database(monkeypatch, database)
    reloaded.reaction_model.generate_thermo(reloaded.initial_species[0])
    assert reloaded.initial_species[0].thermo.label == 'N2v0'


def test_declaration_refuses_duplicate_and_accepts_electronic_manifold(database, tmp_path, monkeypatch):
    job = load_deck(tmp_path, monkeypatch, level=False)
    with pytest.raises(DatabaseError, match='already declared'):
        inp.vibrational_manifold(species='N2')
    install_database(monkeypatch, database)
    job.reaction_model.defer_vibrational_validation = True
    try:
        inp.species('N2A', inp.adjacency_list('electronicstate A3Su+\n' + N2))
    finally:
        job.reaction_model.defer_vibrational_validation = False
    inp.vibrational_manifold(species='N2A')
    assert inp.species_dict['N2A'].props['vibrational_manifold'] == 'N2A'
    assert inp.species_dict['N2A'].molecule[0].electronic_state == 'A3Su+'


def test_declared_thermo_cannot_be_altered_after_assignment(database, tmp_path, monkeypatch):
    job = load_deck(tmp_path, monkeypatch)
    install_database(monkeypatch, database)
    ground = job.initial_species[0]
    job.reaction_model.generate_thermo(ground)
    job.reaction_model.generate_thermo(ground)  # legitimate processed library thermo remains valid
    ground.thermo.change_base_enthalpy(1000)
    with pytest.raises(DatabaseError, match='N2'):
        job.reaction_model.generate_thermo(ground)


def test_library_name_cannot_collide_with_ground_label(database):
    database.libraries['ExcitedStateFixture'].entries['N2A'].label = 'N2'
    model = CoreEdgeReactionModel()
    ground, _ = model.make_new_species(Species().from_smiles('N#N'))
    molecule = Species().from_adjacency_list('electronicstate A3Su+\n' + N2).molecule[0]
    excited, _ = model.make_new_species(molecule, generate_thermo=False)
    with pytest.raises(StateProvenanceError, match="energy transfer"):
        model.generate_thermo(excited, rename=True)
    assert excited.label != ground.label
    assert molecule.state_suffix() in excited.label


def test_model_qm_batch_leaves_resolved_species_to_libraries(database, monkeypatch):
    from rmgpy.qm.main import QMCalculator
    job = RMG()
    job.reaction_model = CoreEdgeReactionModel()
    job.quantum_mechanics = QMCalculator()
    monkeypatch.setattr(inp, 'rmg', job)
    species = Species(label='N2A').from_adjacency_list('electronicstate A3Su+\n' + N2)
    job.reaction_model.new_species_list = [species]
    with pytest.raises(StateProvenanceError, match="energy transfer"):
        job.reaction_model.apply_thermo_to_species(procnum=1)
    assert species.thermo.label == 'N2A'


def test_declaration_keeps_input_label_during_thermo_renaming(database, tmp_path, monkeypatch):
    job = load_deck(tmp_path, monkeypatch, level=False)
    install_database(monkeypatch, database)
    ground = job.initial_species[0]
    job.reaction_model.generate_thermo(ground, rename=True)
    assert ground.label == 'N2'
    assert ground.thermo.label == 'N2v0'


def test_resolved_input_label_survives_batch_thermo_assignment(database, monkeypatch):
    job = RMG()
    job.reaction_model = CoreEdgeReactionModel()
    monkeypatch.setattr(inp, 'rmg', job)
    molecule = Species().from_adjacency_list('electronicstate A3Su+\n' + N2).molecule[0]
    species, _ = job.reaction_model.make_new_species(molecule, label='DeclaredN2A', generate_thermo=False)
    with pytest.raises(StateProvenanceError, match="energy transfer"):
        job.reaction_model.apply_thermo_to_species(procnum=1)
    assert species.label == 'DeclaredN2A'
    assert species.thermo.label == 'N2A'


def test_saved_input_refuses_explicit_resolved_levels(database, tmp_path, monkeypatch):
    job = load_deck(tmp_path, monkeypatch)
    saved = tmp_path / 'refused.py'
    with pytest.raises(SpeciesIdentityError, match='save_input_file'):
        inp.save_input_file(str(saved), job)
    assert not saved.exists()
