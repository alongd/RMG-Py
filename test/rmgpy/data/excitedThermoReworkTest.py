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

from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import pytest
from rmgpy.exceptions import VibrationalManifoldError

import rmgpy.data.rmg as data_module
import rmgpy.rmg.input as inp
from rmgpy import settings
from rmgpy.data.base import Entry
from rmgpy.data.solvation import SolvationDatabase, SoluteData
from rmgpy.data.thermo import ThermoDatabase, ThermoLibrary
from rmgpy.exceptions import DatabaseError, InputError, StateProvenanceError
from rmgpy.qm.main import QMCalculator
from rmgpy.rmg.main import RMG
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.species import Species
from rmgpy.thermo.thermoengine import generate_thermo_data, process_thermo_data

FIXTURES = Path(__file__).resolve().parents[1] / 'test_data' / 'excited_states'
N2 = '1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}'


def nitrogen(state='', label='N2'):
    return Species(label=label).from_adjacency_list(state + N2)


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


def assert_named_error(call, error_name, label):
    with pytest.raises(DatabaseError) as caught:
        call()
    assert type(caught.value).__name__ == error_name
    assert label in str(caught.value)
    assert 'library' in str(caught.value).lower()


@pytest.mark.parametrize('path', ['model', 'species'])
@pytest.mark.parametrize('state', ['electronicstate B3Pg\n', 'electronicstate A3Su+\n'])
def test_attached_ensemble_thermo_is_refused(database, path, state):
    spc = nitrogen(state, 'ExcitedN2')
    spc.thermo = database.get_thermo_data(nitrogen())
    call = spc.get_thermo_data if path == 'species' else lambda: CoreEdgeReactionModel().generate_thermo(spc)
    assert_named_error(call, 'ExcitedSpeciesThermoError', 'ExcitedN2')


@pytest.mark.parametrize('path', ['model', 'species'])
def test_altered_exact_state_thermo_is_refused(database, path):
    spc = nitrogen('electronicstate A3Su+\n', 'AlteredN2')
    spc.thermo = database.get_thermo_data(spc)
    spc.thermo.change_base_enthalpy(1000)
    call = spc.get_thermo_data if path == 'species' else lambda: CoreEdgeReactionModel().generate_thermo(spc)
    assert_named_error(call, 'ExcitedSpeciesThermoError', 'AlteredN2')


def test_state_mutation_reloads_exact_state_thermo(database):
    spc = nitrogen(label='MutableN2')
    spc.get_thermo_data()
    spc.molecule[0].electronic_state = 'A3Su+'
    thermo = spc.get_thermo_data()
    assert thermo.label == 'N2A'
    assert thermo.get_enthalpy(298.15) > 500000
    spc.molecule[0].electronic_state = ''
    assert spc.get_thermo_data().label == 'N2'


def test_state_mutation_to_missing_library_raises(database):
    spc = nitrogen(label='MissingMutableN2')
    spc.get_thermo_data()
    spc.molecule[0].electronic_state = 'B3Pg'
    assert_named_error(spc.get_thermo_data, 'ExcitedSpeciesThermoError', 'MissingMutableN2')


def test_declared_species_api_refuses_ensemble_cache(database):
    spc = nitrogen(label='DeclaredN2')
    model = CoreEdgeReactionModel()
    model.declare_vibrational_manifold(spc)
    spc.thermo = database.get_thermo_data(nitrogen())
    assert_named_error(spc.get_thermo_data, 'VibrationalManifoldError', 'DeclaredN2')


@pytest.mark.parametrize('level', [0, 1])
def test_fixed_level_library_requires_explicit_cpinf(level):
    db = ThermoDatabase()
    lib = ThermoLibrary(label='MissingCpInfFixture')
    lib.load(str(FIXTURES / 'missing_cpinf.py'), db.local_context, {})
    spc = nitrogen('vibrationallevel %d\n' % level, 'FixedN2')
    with pytest.raises(DatabaseError, match='CpInf') as caught:
        db.get_thermo_data_from_library(spc, lib)
    assert type(caught.value).__name__ == 'ExcitedSpeciesThermoError'
    assert 'FixedN2' in str(caught.value)
    assert lib.entries['N2v%d' % level].data.CpInf is None


@pytest.mark.parametrize('path', ['input', 'model'])
def test_duplicate_labels_name_both_states(database, path):
    model = CoreEdgeReactionModel()
    if path == 'input':
        inp.rmg = RMG()
        inp.rmg.reaction_model = model
        inp.rmg.initial_species = []
        inp.species('N2', nitrogen().molecule[0])
        call = lambda: inp.species('N2', nitrogen('electronicstate A3Su+\n').molecule[0])
    else:
        model.make_new_species(nitrogen(), generate_thermo=False)
        call = lambda: model.make_new_species(nitrogen('electronicstate A3Su+\n'), generate_thermo=False)
    with pytest.raises(InputError) as caught:
        call()
    assert type(caught.value).__name__ == 'DuplicateSpeciesLabelError'
    message = str(caught.value)
    assert 'N2' in message and 'electronicstate A3Su+' in message
    assert message.count('1 N u0 p1 c0') == 2
    if path == 'input':
        assert inp.species_dict['N2'].molecule[0].electronic_state == ''


@pytest.mark.parametrize('reverse', [False, True])
def test_automatic_state_label_is_unique_against_all_species(database, reverse):
    model = CoreEdgeReactionModel()
    state = 'electronicstate A3Su+\n'
    auto_label = nitrogen(state, '').smiles + nitrogen(state).molecule[0].state_suffix()
    if reverse:
        automatic, _ = model.make_new_species(nitrogen(state, ''), generate_thermo=False)
        explicit, _ = model.make_new_species(Species(label=auto_label).from_smiles('C'), generate_thermo=False)
    else:
        explicit, _ = model.make_new_species(Species(label=auto_label).from_smiles('C'), generate_thermo=False)
        automatic, _ = model.make_new_species(nitrogen(state, ''), generate_thermo=False)
    assert explicit.label == auto_label
    assert automatic.label != explicit.label
    assert automatic.molecule[0].state_suffix() in automatic.label


def test_ground_declaration_cannot_cover_excited_electronic_manifold(database):
    model = CoreEdgeReactionModel()
    model.declare_vibrational_manifold(nitrogen())
    model.make_new_species(nitrogen('electronicstate A3Su+\n', 'N2A'), generate_thermo=False)
    level = nitrogen('electronicstate A3Su+\nvibrationallevel 1\n', 'N2Av1')
    with pytest.raises(DatabaseError, match='N2Av1') as caught:
        model.make_new_species(level, generate_thermo=False)
    assert type(caught.value).__name__ == 'VibrationalManifoldError'


def test_electronic_manifold_declaration_keeps_electronic_state(database):
    model = CoreEdgeReactionModel()
    spc = nitrogen('electronicstate A3Su+\n', 'N2A')
    lib = database.libraries['ExcitedStateFixture']
    item = nitrogen('electronicstate A3Su+\nvibrationallevel 0\n', 'N2Av0').molecule[0]
    lib.entries['N2Av0'] = Entry(index=4, label='N2Av0', item=item,
                                data=deepcopy(lib.entries['N2v0'].data))
    try:
        model.declare_vibrational_manifold(spc)
    except DatabaseError as error:
        pytest.fail('Electronic manifold declaration refused: %s' % error)
    with pytest.raises(StateProvenanceError, match='energy transfer'):
        model.generate_thermo(spc)
    assert spc.thermo.label == 'N2Av0'
    assert spc.molecule[0].electronic_state == 'A3Su+'
    assert spc.molecule[0].vibrational_level == -1
    assert spc.get_thermo_data().label == 'N2Av0'
    model.make_new_species(nitrogen('electronicstate A3Su+\nvibrationallevel 1\n', 'N2Av1'),
                           generate_thermo=False)
    with pytest.raises(DatabaseError, match='duplicates v = 0'):
        model.make_new_species(nitrogen('electronicstate A3Su+\nvibrationallevel 0\n', 'N2Av0'),
                               generate_thermo=False)


@pytest.fixture
def water_database(database, monkeypatch):
    solvation = SolvationDatabase()
    solvation.load(str(Path(settings['database.directory']) / 'solvation'))
    monkeypatch.setattr(data_module, 'database', SimpleNamespace(thermo=database, solvation=solvation))
    return solvation


@pytest.mark.parametrize('declared', [False, True])
@pytest.mark.parametrize('path', ['engine', 'solute', 'groups'])
def test_solvation_requires_exact_state_library(database, water_database, declared, path, monkeypatch):
    spc = nitrogen('electronicstate A3Su+\n' if not declared else '', 'SolvatedN2')
    if declared:
        CoreEdgeReactionModel().declare_vibrational_manifold(spc)
    thermo = database.get_thermo_data(spc)
    calls = {'engine': lambda: process_thermo_data(spc, thermo, solvent_name='water'),
             'solute': lambda: water_database.get_solute_data(spc),
             'groups': lambda: water_database.get_solute_data_from_groups(spc)}
    assert_named_error(calls[path], 'VibrationalManifoldError' if declared else 'ExcitedSpeciesThermoError',
                       'SolvatedN2')
    if path == 'engine':
        monkeypatch.setattr(data_module, 'database', SimpleNamespace(thermo=database, solvation=None))
        assert_named_error(calls[path], 'VibrationalManifoldError' if declared else 'ExcitedSpeciesThermoError',
                           'SolvatedN2')


def test_declared_solute_lookup_uses_v0_library(database, water_database):
    spc = nitrogen(label='SolvatedDeclaredN2')
    CoreEdgeReactionModel().declare_vibrational_manifold(spc)
    expected = SoluteData(S=0.1, B=0.2, E=0.3, L=0.4, A=0.5, V=0.6)
    water_database.libraries['solute'].entries['N2v0'] = Entry(
        label='N2v0', item=nitrogen('vibrationallevel 0\n'), data=expected)
    actual = water_database.get_solute_data(spc)
    assert actual.S == expected.S
    assert actual.V == expected.V
    assert 'N2v0' in actual.comment
    gas = database.get_thermo_data(spc)
    # Rework 2: exact solute data does not authorize derived liquid thermo.
    with pytest.raises(VibrationalManifoldError, match='SolvatedDeclaredN2'):
        process_thermo_data(spc, gas, solvent_name='water')


@pytest.mark.parametrize('keyword', ['electronicstate Foo', 'vibrationallevel 0', 'vibrationallevel 1'])
def test_electron_state_keywords_are_refused(keyword):
    with pytest.raises(DatabaseError, match='electron') as caught:
        Species(label='electron').from_adjacency_list(keyword + '\n1 e u0 p0 c-1')
    assert type(caught.value).__name__ == 'ExcitedSpeciesThermoError'


@pytest.mark.parametrize('path', ['lookup', 'provenance'])
def test_mutated_electron_state_cannot_bypass_policy(database, path):
    spc = Species(label='electron').from_adjacency_list('1 e u0 p0 c-1')
    spc.molecule[0].electronic_state = 'Foo'
    if path == 'lookup':
        call = lambda: database.get_thermo_data_from_library(spc, database.libraries['ExcitedStateFixture'])
    else:
        from rmgpy.solver.plasma import PlasmaReactor
        call = lambda: PlasmaReactor(T=(300,'K'), P=(1,'bar'), initial_mole_fractions={}, Te=(10000,'K'))._check_charged_species_thermo_provenance([spc], [])
    assert_named_error(call, 'ExcitedSpeciesThermoError', 'electron')
    spc.molecule[0].electronic_state = ''
    spc.props['vibrational_manifold'] = 'electron'
    assert_named_error(call, 'ExcitedSpeciesThermoError', 'electron')


@pytest.mark.parametrize('declared', [False, True])
def test_missing_database_raises_named_error(database, monkeypatch, declared):
    spc = nitrogen('electronicstate B3Pg\n' if not declared else '', 'NoDatabaseN2')
    if declared:
        CoreEdgeReactionModel().declare_vibrational_manifold(spc)
    monkeypatch.setattr(data_module, 'database', None)
    assert_named_error(lambda: generate_thermo_data(spc),
                       'VibrationalManifoldError' if declared else 'ExcitedSpeciesThermoError', 'NoDatabaseN2')


def test_qm_batch_dispatches_unresolved_iterator_once(monkeypatch):
    qm = QMCalculator(software='mopac', method='pm3', onlyCyclics=True, maxRadicalNumber=0)
    jobs = []
    monkeypatch.setattr('rmgpy.qm.main._write_qm_files_star', lambda args: jobs.append(args[1]))
    spc = Species().from_smiles('C1CCCCC1')
    qm.run_jobs(iter([spc]), procnum=1)
    assert jobs == [spc.molecule[0]]


@pytest.mark.parametrize('property', ['get_heat_capacity', 'get_enthalpy', 'get_entropy', 'get_free_energy'])
def test_statmech_cannot_supply_resolved_thermo_without_library(database, property):
    from rmgpy.statmech import Conformer, IdealGasTranslation
    spc = nitrogen('electronicstate B3Pg\n', 'StatmechN2')
    spc.conformer = Conformer(E0=(0,'kJ/mol'), modes=[IdealGasTranslation(mass=(28,'amu'))])
    assert_named_error(lambda: getattr(spc, property)(300), 'ExcitedSpeciesThermoError', 'StatmechN2')
