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

import pytest

from rmgpy.data.thermo import ThermoDatabase, ThermoLibrary
from rmgpy.exceptions import DatabaseError
from rmgpy.ml.estimator import MLEstimator
from rmgpy.qm.main import QMCalculator
from rmgpy.qm.molecule import QMMolecule
from rmgpy.species import Species

FIXTURES = Path(__file__).resolve().parents[1] / 'test_data' / 'excited_states'
N2 = '1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}'


def nitrogen(state='', label='N2'):
    return Species(label=label).from_adjacency_list(state + N2)


@pytest.fixture
def thermo_database():
    library = ThermoLibrary(label='ExcitedStateFixture')
    library.load(str(FIXTURES / 'thermo.py'), ThermoDatabase().local_context, {})
    database = ThermoDatabase()
    database.libraries = {library.label: library}
    database.library_order = [library.label]
    return database


@pytest.mark.parametrize('state,label,h298', [
    ('vibrationallevel 1\n', 'N2v1', 28637.360782),
    ('electronicstate A3Su+\n', 'N2A', 590691.029862),
])
def test_resolved_thermo_is_from_exact_library_entry(thermo_database, state, label, h298):
    thermo = thermo_database.get_thermo_data(nitrogen(state, label))
    assert thermo.label == label
    assert 'Thermo library: ExcitedStateFixture' in thermo.comment
    assert thermo.get_enthalpy(298.15) == pytest.approx(h298, rel=0.001)


@pytest.mark.parametrize('path', [
    'lookup', 'groups', 'group_additivity', 'compute_groups', 'hbi',
    'surface', 'database_ml', 'ml_molecule', 'ml_species', 'qm',
])
@pytest.mark.parametrize('state', ['vibrationallevel 2\n', 'electronicstate B3Pg\n'])
def test_resolved_thermo_estimators_refuse_with_named_error(path, state):
    database = ThermoDatabase()
    species = nitrogen(state, 'MissingN2State')
    molecule = species.molecule[0]
    ml = MLEstimator.__new__(MLEstimator)
    qm = QMCalculator()
    calls = {
        'lookup': lambda: database.get_thermo_data(species),
        'groups': lambda: database.get_thermo_data_from_groups(species),
        'group_additivity': lambda: database.estimate_thermo_via_group_additivity(molecule),
        'compute_groups': lambda: database.compute_group_additivity_thermo(molecule),
        'hbi': lambda: database.estimate_radical_thermo_via_hbi(molecule, None),
        'surface': lambda: database.get_thermo_data_for_surface_species(species),
        'database_ml': lambda: database.get_thermo_data_from_ml(species, ml, {}),
        'ml_molecule': lambda: ml.get_thermo_data(molecule),
        'ml_species': lambda: ml.get_thermo_data_for_species(species),
        'qm': lambda: qm.get_thermo_data(molecule),
    }
    with pytest.raises(DatabaseError) as caught:
        calls[path]()
    assert type(caught.value).__name__ == 'ExcitedSpeciesThermoError'
    assert 'library' in str(caught.value).lower()
    assert state.strip() in str(caught.value)
    if path in ('lookup', 'groups', 'surface', 'database_ml', 'ml_species'):
        assert 'MissingN2State' in str(caught.value)


@pytest.mark.parametrize('path', ['all', 'depository', 'qm_generate', 'qm_calculate', 'qm_cache'])
def test_other_thermo_entry_points_cannot_bypass_state_guard(path):
    database = ThermoDatabase()
    species = nitrogen('vibrationallevel 2\n', 'MissingN2State')
    qm = QMMolecule(species.molecule[0], None)
    calls = {
        'all': lambda: database.get_all_thermo_data(species),
        'depository': lambda: database.get_thermo_data_from_depository(species),
        'qm_generate': lambda: qm.generate_thermo_data(),
        'qm_calculate': lambda: qm.calculate_thermo_data(),
        'qm_cache': lambda: qm.load_thermo_data(),
    }
    with pytest.raises(DatabaseError) as caught:
        calls[path]()
    assert type(caught.value).__name__ == 'ExcitedSpeciesThermoError'


def test_all_thermo_query_returns_only_exact_state_libraries(thermo_database):
    matches = thermo_database.get_all_thermo_data(nitrogen('vibrationallevel 1\n', 'N2v1'))
    assert len(matches) == 1
    assert matches[0][2].label == 'N2v1'


def test_direct_qm_batch_refuses_resolved_species():
    with pytest.raises(DatabaseError) as caught:
        QMCalculator().run_jobs([nitrogen('electronicstate B3Pg\n', 'MissingN2State')])
    assert type(caught.value).__name__ == 'ExcitedSpeciesThermoError'
    assert 'MissingN2State' in str(caught.value)


def test_surface_binding_scaling_refuses_resolved_thermo(thermo_database):
    species = nitrogen('electronicstate A3Su+\n', 'N2A')
    thermo = thermo_database.get_thermo_data(species)
    with pytest.raises(DatabaseError) as caught:
        thermo_database.correct_binding_energy(thermo, species)
    assert type(caught.value).__name__ == 'ExcitedSpeciesThermoError'
