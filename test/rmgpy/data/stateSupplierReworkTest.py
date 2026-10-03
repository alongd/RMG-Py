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

import rmgpy
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np
import pytest
from pysidt import MultiTargetSingleEvalSubgraphIsomorphicDecisionTree
from pysidt.sidt import Node, Datum, write_nodes

from rmgpy.data.thermo import ThermoDatabase
from rmgpy.exceptions import ExcitedSpeciesThermoError, StateProvenanceError
from rmgpy.ml.estimator import MLEstimator
from rmgpy.molecule import Molecule
from rmgpy.qm.main import QMCalculator
from rmgpy.qm.gaussian import GaussianMolPM3
from rmgpy.species import Species
from rmgpy.thermo import ThermoData

STATES = [('A', -1), ('', 1), ('A', 1)]

@pytest.mark.parametrize('state', STATES)
def test_loaded_fitted_sidt_root_refuses_before_tagging(state, tmp_path):
    ground = Molecule(smiles='[CH2]=*')
    values = np.array([12345., 12., 75., 90., 100., 120., 140., 150., 170.])
    root = Node(name='Root', group=None, items=[Datum(ground, values), Datum(ground.copy(deep=True), values)])
    tree = MultiTargetSingleEvalSubgraphIsomorphicDecisionTree(nodes={'Root': root}, target_weights=np.ones(9))
    tree.fit_node(root, skip_val=True)
    write_nodes(tree, str(tmp_path / 'Pt111_monodentate_adsorption_corrections.json'))
    db = ThermoDatabase()
    db.adsorption_groups = 'SIDT'
    db.load_sidts(str(tmp_path))
    mol = Molecule(smiles='[CH2]=*', electronic_state=state[0], vibrational_level=state[1])
    thermo = ThermoData(Tdata=([300,400,500,600,800,1000,1500], 'K'), Cpdata=([0]*7, 'J/(mol*K)'), H298=(0,'J/mol'), S298=(0,'J/(mol*K)'))
    before = [a.label for a in mol.atoms]
    with pytest.raises(StateProvenanceError, match='Adsorption'):
        db._add_adsorption_correction(thermo, None, mol, mol.get_surface_sites())
    assert thermo.H298.value_si == 0
    assert [a.label for a in mol.atoms] == before
    unresolved = ground.copy(deep=True)
    assert db._add_adsorption_correction(thermo, None, unresolved, unresolved.get_surface_sites())
    assert thermo.H298.value_si == 12345000

@pytest.mark.parametrize('state', STATES)
def test_qm_refuses_before_initialization(state):
    calculator = QMCalculator(software='gaussian', method='pm3')
    calculator.initialize = Mock(side_effect=AssertionError('QM initialization reached'))
    mol = Molecule(smiles='C=C', electronic_state=state[0], vibrational_level=state[1])
    with pytest.raises(StateProvenanceError, match='QM'):
        calculator.get_thermo_data(mol)
    calculator.initialize.assert_not_called()

@pytest.mark.parametrize('state', STATES)
def test_qm_job_refuses_before_cache_read(state, tmp_path):
    calculator = QMCalculator(software='gaussian', method='pm3', fileStore=str(tmp_path), scratchDirectory=str(tmp_path))
    job = GaussianMolPM3(Molecule(smiles='C=C'), calculator.settings)
    job.molecule.electronic_state, job.molecule.vibrational_level = state
    job.initialize = Mock()
    job.load_thermo_data = Mock(side_effect=AssertionError('QM cache reached'))
    with pytest.raises(StateProvenanceError, match='QM'):
        job.generate_thermo_data()
    job.initialize.assert_not_called()
    job.load_thermo_data.assert_not_called()

@pytest.mark.parametrize('state', STATES)
def test_qm_full_thermo_refuses_ground_result(state, tmp_path):
    from rmgpy.qm.symmetry import POINT_GROUP_DICTIONARY
    import rmgpy.rmg.input as input_module
    calculator = QMCalculator(software='gaussian', method='pm3', fileStore=str(tmp_path), scratchDirectory=str(tmp_path), onlyCyclics=False)
    calculator.initialize()
    mol = Molecule(smiles='C=C', electronic_state=state[0], vibrational_level=state[1])
    # Put the real ground Gaussian fixture at the would-be resolved cache key.
    ground_job = GaussianMolPM3(Molecule(smiles='C=C'), calculator.settings)
    ground_job.molecule = mol
    ground_job.unique_id = mol.to_augmented_inchi_key()
    ground_job.unique_id_long = mol.to_augmented_inchi()
    fixture = (Path(rmgpy.__file__).resolve().parents[1] / 'arkane/data/gaussian/ethylene.log')
    ground_job.check_file_names()
    Path(ground_job.output_file_path).write_text(fixture.read_text() + '\n' + ground_job.unique_id_long + '\n')
    with patch.object(input_module, 'rmg', SimpleNamespace(quantum_mechanics=calculator, ml_estimator=None, ml_settings=None)), patch.object(GaussianMolPM3, 'determine_point_group', lambda job: setattr(job, 'point_group', POINT_GROUP_DICTIONARY['D2h'])):
        with pytest.raises(ExcitedSpeciesThermoError, match='Library-only'):
            ThermoDatabase().get_thermo_data(Species(molecule=[mol]))

@pytest.mark.parametrize('state', STATES)
def test_ml_refuses_before_prediction(state):
    estimator = MLEstimator.__new__(MLEstimator)
    estimator.hf298_estimator = Mock(return_value=[[12.]])
    estimator.s298_cp_estimator = Mock(return_value=[[20.] + [5.]*7])
    mol = Molecule(smiles='CC', electronic_state=state[0], vibrational_level=state[1])
    with pytest.raises(StateProvenanceError, match='ML'):
        estimator.get_thermo_data(mol)
    with pytest.raises(StateProvenanceError, match='ML'):
        estimator.get_thermo_data_for_species(Species(molecule=[mol]))
    estimator.hf298_estimator.assert_not_called()
    estimator.s298_cp_estimator.assert_not_called()
    assert estimator.get_thermo_data('CC').H298.value_si == 12. * 4184.

@pytest.mark.parametrize('state', STATES)
def test_energy_transfer_estimation_refuses_without_overwriting(state):
    mol = Molecule(smiles='CC', electronic_state=state[0], vibrational_level=state[1])
    explicit = object()
    species = Species(molecule=[mol], energy_transfer_model=explicit)
    with pytest.raises(StateProvenanceError, match='energy transfer'):
        species.generate_energy_transfer_model()
    assert species.energy_transfer_model is explicit
    ground = Species(smiles='CC')
    ground.generate_energy_transfer_model()
    assert ground.energy_transfer_model.alpha0.value_si == 300 * 0.011962 * 1000


@pytest.mark.parametrize('state', STATES)
def test_arkane_statmech_refuses_before_file_read(state, tmp_path):
    from arkane.statmech import StatMechJob
    species = Species(molecule=[Molecule(smiles='CC', electronic_state=state[0], vibrational_level=state[1])])
    with pytest.raises(ExcitedSpeciesThermoError, match='Library-only'):
        StatMechJob(species, str(tmp_path / 'missing.py')).load()
