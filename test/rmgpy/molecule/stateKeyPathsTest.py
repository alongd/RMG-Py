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

"""Resolved augmented keys must work as real QM filenames."""
from pathlib import Path
import json
import os
import pickle
import re
import string
import sys

import pytest
from rdkit import Chem

from rmgpy.molecule import Molecule
from rmgpy.exceptions import InvalidAdjacencyListError
from rmgpy.molecule.adjlist import from_adjacency_list
from rmgpy.molecule.translator import to_inchi_key
from rmgpy.qm.main import QMSettings
from rmgpy.qm.molecule import QMMolecule, Geometry
from rmgpy.qm.gaussian import GaussianMolPM3
from rmgpy.qm.mopac import MopacMolPM3
from rmgpy.qm.qmdata import QMData
from rmgpy.qm.symmetry import SymmetryJob, POINT_GROUP_DICTIONARY
from rmgpy.thermo import ThermoData


@pytest.fixture
def settings(tmp_path):
    return QMSettings(software='mopac', method='pm3',
                      fileStore=str(tmp_path / 'permanent'),
                      scratchDirectory=str(tmp_path / 'scratch'))


@pytest.mark.parametrize('state', [('', -1), ('', 1), ('A3Su+', -1)])
def test_qm_geometry_writes_real_state_specific_mol_files(settings, state):
    molecule = Molecule(smiles='N#N', electronic_state=state[0], vibrational_level=state[1])
    qm = QMMolecule(molecule, settings)
    geometry = qm.create_geometry()
    for name in (geometry.get_crude_mol_file_path(), geometry.get_refined_mol_file_path()):
        path = Path(name)
        assert path.parent == Path(settings.scratchDirectory)
        assert path.is_file()
        assert Chem.MolFromMolFile(str(path)).GetNumAtoms() == 2
    assert qm.unique_id_long == molecule.to_augmented_inchi()
    assert (molecule.electronic_state, molecule.vibrational_level) == state


RESOLVED_STATES = [('', 1), ('A3Su+', -1), ('A3Su+', 2)]


def nitrogen_qm_data():
    # Known input data tests file plumbing; no external QM calculation is run.
    return QMData(groundStateDegeneracy=1, numberOfAtoms=2, atomicNumbers=[7, 7],
                  molecularMass=(28, 'amu'), energy=(0, 'kJ/mol'),
                  atomCoords=([[0, 0, 0], [0, 0, 1.1]], 'angstrom'),
                  rotationalConstants=([60e9, 60e9], 'Hz'), frequencies=([2200], 'cm^-1'))


@pytest.mark.parametrize('state,key_suffix', [
    (('', -1), ''), (('', 0), '-v0'), (('', 1), '-v1'),
    (('A3Su+', -1), '-es413353752b'), (('A3Su+', 2), '-es413353752b-v2')])
@pytest.mark.parametrize('level', [1, 2])
def test_augmented_key_formats_are_safe_and_unresolved_key_is_frozen(state, key_suffix, level):
    molecule = Molecule(smiles='N#N', electronic_state=state[0], vibrational_level=state[1])
    base = 'IJGRMHOSHXDMSA-UHFFFAOYSA-N' + ('-mult1' if level == 1 else '')
    assert to_inchi_key(molecule, aug_level=level) == base + key_suffix
    if level == 2:
        assert molecule.to_augmented_inchi_key() == base + key_suffix


def test_all_token_characters_and_layer_boundaries_have_distinct_cache_filenames(tmp_path):
    tokens = ['', 'A', 'a', 'A-v1', 'A-es41']
    tokens += ['A' + character for character in string.ascii_letters + string.digits + '+-_.,()']
    tokens = list(dict.fromkeys(tokens))
    states = [(token, level) for token in tokens for level in (-1, 0, 1, 2, 10, 2147483647)]
    cache = {}
    filenames = []
    for state in states:
        key = Molecule(smiles='N#N', electronic_state=state[0], vibrational_level=state[1]).to_augmented_inchi_key()
        path = tmp_path / (key + '.json')
        path.write_text(json.dumps(state))
        assert path.parent == tmp_path
        assert re.fullmatch(r'[A-Za-z0-9-]+', key)
        cache[key] = state
        filenames.append(path.name.casefold())
    assert len(cache) == len(states)
    assert len(set(filenames)) == len(states)
    for key, state in cache.items():
        assert json.loads((tmp_path / (key + '.json')).read_text()) == list(state)


@pytest.mark.parametrize('smiles,multiplicity', [('[CH2]C=C', 2), ('[CH2]', 1)])
@pytest.mark.parametrize('state', RESOLVED_STATES)
def test_key_state_suffix_preserves_existing_unpaired_layers(smiles, multiplicity, state):
    ground = Molecule(smiles=smiles, multiplicity=multiplicity)
    molecule = Molecule(smiles=smiles, multiplicity=multiplicity,
                        electronic_state=state[0], vibrational_level=state[1])
    suffixes = {('', 1): '-v1', ('A3Su+', -1): '-es413353752b', ('A3Su+', 2): '-es413353752b-v2'}
    for level in (1, 2):
        assert to_inchi_key(molecule, aug_level=level) == to_inchi_key(ground, aug_level=level) + suffixes[state]


@pytest.mark.parametrize('state', RESOLVED_STATES)
def test_key_state_suffix_preserves_real_lone_pair_layer(state):
    ground = Molecule().from_adjacency_list(
        'multiplicity 1\n1 C u0 p1 c0 {2,S} {3,S}\n'
        '2 H u0 p0 c0 {1,S}\n3 H u0 p0 c0 {1,S}')
    assert ground.to_augmented_inchi_key() == 'HZVOZRGWRWCICA-UHFFFAOYSA-N-lp1'
    molecule = ground.copy(deep=True)
    molecule.electronic_state, molecule.vibrational_level = state
    suffixes = {('', 1): '-v1', ('A3Su+', -1): '-es413353752b', ('A3Su+', 2): '-es413353752b-v2'}
    for level in (1, 2):
        assert to_inchi_key(molecule, aug_level=level) == to_inchi_key(ground, aug_level=level) + suffixes[state]


def test_longest_electronic_state_round_trips():
    molecule = Molecule(smiles='N#N', electronic_state='A' * 32, vibrational_level=2147483647)
    restored = [Molecule().from_adjacency_list(molecule.to_adjacency_list()),
                Molecule().from_augmented_inchi(molecule.to_augmented_inchi()),
                pickle.loads(pickle.dumps(molecule)), eval(repr(molecule))]
    for other in restored:
        assert other.is_isomorphic(molecule)
        assert other.electronic_state == 'A' * 32
        assert other.vibrational_level == 2147483647


def test_longest_electronic_state_qm_files_fit_name_max(settings):
    molecule = Molecule(smiles='N#N', electronic_state='A' * 32, vibrational_level=2147483647)
    geometry = QMMolecule(molecule, settings).create_geometry()
    name_max = os.pathconf(settings.scratchDirectory, 'PC_NAME_MAX')
    for name in (geometry.get_crude_mol_file_path(), geometry.get_refined_mol_file_path()):
        path = Path(name)
        assert len(os.fsencode(path.name)) <= name_max
        assert path.parent == Path(settings.scratchDirectory)
        assert path.is_file()
        assert Chem.MolFromMolFile(str(path)).GetNumAtoms() == 2


@pytest.mark.parametrize('length', [33, 107])
@pytest.mark.parametrize('entry', ['constructor', 'assignment', 'adjacency', 'parser', 'aug_inchi', 'qm'])
def test_overlong_electronic_state_is_refused_before_file_writes(entry, length, settings, tmp_path):
    token = 'A' * length
    molecule = Molecule(smiles='N#N', electronic_state='B3Pg', vibrational_level=1)
    original = molecule.to_adjacency_list()
    adjacency = 'electronicstate ' + token + '\n1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}'
    error = InvalidAdjacencyListError if entry in ('adjacency', 'parser') else ValueError
    with pytest.raises(error, match='32'):
        if entry == 'constructor':
            Molecule(smiles='N#N', electronic_state=token)
        elif entry == 'assignment':
            molecule.electronic_state = token
        elif entry == 'adjacency':
            molecule.from_adjacency_list(adjacency)
        elif entry == 'parser':
            from_adjacency_list(adjacency, state={})
        elif entry == 'aug_inchi':
            molecule.from_augmented_inchi('InChI=1S/N2/c1-2/es:' + token)
        else:
            QMMolecule(Molecule(smiles='N#N', electronic_state=token), settings).create_geometry()
    assert molecule.to_adjacency_list() == original
    assert not list(tmp_path.rglob('*'))


@pytest.mark.parametrize('state', RESOLVED_STATES)
def test_qm_thermo_cache_saves_and_loads_resolved_identity(settings, state):
    molecule = Molecule(smiles='N#N', electronic_state=state[0], vibrational_level=state[1])
    qm = QMMolecule(molecule, settings)
    qm.check_paths()
    qm.thermo = ThermoData(Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
                          Cpdata=([29] * 7, 'J/(mol*K)'), H298=(1234, 'J/mol'), S298=(200, 'J/(mol*K)'))
    qm.point_group = POINT_GROUP_DICTIONARY['Dinfh']
    qm.qm_data = nitrogen_qm_data()
    qm.save_thermo_data()
    path = Path(qm.get_thermo_file_path())
    assert path.parent == Path(settings.fileStore)
    assert path.is_file()
    restored = QMMolecule(molecule.copy(deep=True), settings)
    assert restored.load_thermo_data().H298.value_si == pytest.approx(1234)
    assert restored.unique_id_long in path.read_text()
    assert restored.get_augmented_inchi_key() == qm.unique_id


@pytest.mark.parametrize('state', RESOLVED_STATES)
def test_symmetry_coordinate_file_uses_resolved_key(settings, state):
    molecule = Molecule(smiles='N#N', electronic_state=state[0], vibrational_level=state[1])
    qm = QMMolecule(molecule, settings)
    qm.check_paths()
    job = SymmetryJob(settings, qm.unique_id, nitrogen_qm_data())
    job.write_input_file()
    path = Path(job.input_file_path)
    assert path.parent == Path(settings.scratchDirectory)
    assert path.read_text().splitlines()[0] == '2'
    assert len(path.read_text().splitlines()) == 3


@pytest.mark.parametrize('state', RESOLVED_STATES)
@pytest.mark.parametrize('calculator', [GaussianMolPM3, MopacMolPM3])
def test_program_input_and_output_files_use_resolved_key(settings, state, calculator):
    molecule = Molecule(smiles='N#N', electronic_state=state[0], vibrational_level=state[1])
    qm = calculator(molecule, settings)
    qm.create_geometry()
    qm.write_input_file(1)
    input_path, output_path = Path(qm.input_file_path), Path(qm.output_file_path)
    assert input_path.parent == output_path.parent == Path(settings.scratchDirectory)
    assert qm.unique_id_long in input_path.read_text()
    output_path.write_text('file-path regression, no QM executable was run')
    assert output_path.read_text().startswith('file-path regression')


# Full filenames also include variable-length unpaired-electron layers.
def radical_chain(count, resolved=False, multiplicity=None):
    kwargs = {} if multiplicity is None else {'multiplicity': multiplicity}
    molecule = Molecule(smiles='[CH2]' + '[CH]' * (count - 2) + '[CH2]', **kwargs)
    if resolved:
        molecule.electronic_state = 'A' * 32
        molecule.vibrational_level = 2147483647
    return molecule


def assert_qm_filename_error(error, species, filename, byte_count, limit):
    assert type(error).__name__ == 'QMFileNameError'
    message = str(error)
    assert species in message
    assert filename in message
    assert str(byte_count) + ' bytes' in message
    assert 'NAME_MAX=' + str(limit) in message


@pytest.mark.parametrize('resolved,count,key_bytes', [(True, 49, 245), (False, 75, 244)])
@pytest.mark.parametrize('existing', [False, True])
@pytest.mark.parametrize('entry', ['qm', 'geometry', 'gaussian', 'mopac'])
def test_full_qm_filename_budget_refuses_before_any_file(settings, tmp_path, resolved, count, key_bytes,
                                                        existing, entry):
    molecule = radical_chain(count, resolved)
    key = molecule.to_augmented_inchi_key()
    assert len(os.fsencode(key)) == key_bytes
    if existing:
        Path(settings.scratchDirectory).mkdir(parents=True)
        Path(settings.fileStore).mkdir(parents=True)
    with pytest.raises(ValueError) as exc:
        if entry == 'geometry':
            geometry = Geometry(settings, key, molecule, molecule.to_augmented_inchi())
        else:
            calculator = {'qm': QMMolecule, 'gaussian': GaussianMolPM3, 'mopac': MopacMolPM3}[entry]
            qm = calculator(molecule, settings)
            geometry = Geometry(settings, qm.unique_id, molecule, qm.unique_id_long)
        rd_molecule, _ = geometry.rd_build()
        geometry.rd_embed(rd_molecule, 1)
    assert_qm_filename_error(exc.value, molecule.to_augmented_inchi(), key + '.refined.mol', key_bytes + 12, 255)
    for directory in (settings.scratchDirectory, settings.fileStore):
        assert not Path(directory).exists() or not list(Path(directory).iterdir())
    if not existing:
        assert not list(tmp_path.iterdir())


@pytest.mark.parametrize('resolved,count,longest_bytes', [(True, 48, 254), (False, 74, 253)])
@pytest.mark.parametrize('calculator', [GaussianMolPM3, MopacMolPM3])
def test_just_under_budget_key_writes_every_qm_file(settings, monkeypatch, resolved, count, longest_bytes, calculator):
    # Only file plumbing is under test; use the supported singlet MOPAC input keywords.
    molecule = radical_chain(count, resolved, multiplicity=1)
    qm = calculator(molecule, settings)
    qm.check_paths()
    geometry = Geometry(settings, qm.unique_id, molecule, qm.unique_id_long)
    rd_molecule, mapping = geometry.rd_build()
    rd_molecule, conformer = geometry.rd_embed(rd_molecule, 1)
    geometry.save_coordinates_from_rdmol(rd_molecule, conformer, mapping)
    qm.geometry = geometry
    qm.write_input_file(1)

    # The external executable boundary writes real output/log files; no QM job is run.
    class FileWritingProcess:
        def __init__(self, command, **kwargs):
            if calculator is MopacMolPM3:
                self.output = Path(command[1]).with_suffix('.out')
            else:
                self.output = Path(command[2])

        def communicate(self, **kwargs):
            self.output.write_text('QM filename regression fixture, no calculated results')
            return b'', b'ended normally'

    module = 'rmgpy.qm.mopac' if calculator is MopacMolPM3 else 'rmgpy.qm.gaussian'
    monkeypatch.setattr(module + '.Popen', FileWritingProcess)
    qm.executable_path = sys.executable
    assert qm.run() is False  # The fixture deliberately contains no successful QM result.
    coordinates = [[atom.coords[0], atom.coords[1], atom.coords[2]] for atom in molecule.atoms]
    qm.qm_data = QMData(groundStateDegeneracy=molecule.multiplicity, numberOfAtoms=len(molecule.atoms),
                        atomicNumbers=[atom.element.number for atom in molecule.atoms], molecularMass=(100, 'amu'),
                        energy=(0, 'kJ/mol'), atomCoords=(coordinates, 'angstrom'),
                        rotationalConstants=([60e9] * 3, 'Hz'), frequencies=([2200], 'cm^-1'))
    SymmetryJob(settings, qm.unique_id, qm.qm_data).write_input_file()
    qm.thermo = ThermoData(Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
                          Cpdata=([29] * 7, 'J/(mol*K)'), H298=(1234, 'J/mol'), S298=(200, 'J/(mol*K)'))
    qm.point_group = POINT_GROUP_DICTIONARY['C1']
    qm.save_thermo_data()
    scratch_suffixes = ('.crude.mol', '.refined.mol', '.symm', qm.input_file_extension, qm.output_file_extension)
    expected = {Path(settings.scratchDirectory) / (qm.unique_id + suffix) for suffix in scratch_suffixes}
    expected.add(Path(settings.fileStore) / (qm.unique_id + '.thermo'))
    actual = set(Path(settings.scratchDirectory).iterdir()) | set(Path(settings.fileStore).iterdir())
    assert actual == expected
    assert max(len(os.fsencode(path.name)) for path in actual) == longest_bytes
    for path in actual:
        assert path.read_text()
        assert len(os.fsencode(path.name)) <= os.pathconf(path.parent, 'PC_NAME_MAX')
    for path in (Path(geometry.get_crude_mol_file_path()), Path(geometry.get_refined_mol_file_path())):
        assert Chem.MolFromMolFile(str(path), removeHs=False).GetNumAtoms() == len(molecule.atoms)


@pytest.mark.parametrize('target,suffix,limit,byte_count', [
    ('scratchDirectory', '.refined.mol', 117, 118), ('fileStore', '.thermo', 112, 113)])
def test_qm_filename_budget_uses_each_existing_directory_limit(settings, monkeypatch, target, suffix, limit, byte_count):
    molecule = Molecule(smiles='N#N', electronic_state='A' * 32, vibrational_level=2147483647)
    for directory in (settings.scratchDirectory, settings.fileStore):
        Path(directory).mkdir(parents=True)
    actual_pathconf = os.pathconf

    def directory_limit(path, name):
        if Path(path) == Path(getattr(settings, target)):
            assert name == 'PC_NAME_MAX'
            return limit
        return actual_pathconf(path, name)

    monkeypatch.setattr(os, 'pathconf', directory_limit)
    with pytest.raises(ValueError) as exc:
        QMMolecule(molecule, settings)
    assert_qm_filename_error(exc.value, molecule.to_augmented_inchi(), molecule.to_augmented_inchi_key() + suffix,
                             byte_count, limit)
    assert getattr(settings, target) in str(exc.value)
    for directory in (settings.scratchDirectory, settings.fileStore):
        assert not list(Path(directory).iterdir())


@pytest.mark.parametrize('field', ['input_file_extension', 'output_file_extension'])
def test_qm_budget_includes_declared_program_input_and_log_before_geometry(settings, field):
    class LongFileGaussian(GaussianMolPM3):
        pass

    setattr(LongFileGaussian, field, '.log' + 'x' * 150)
    molecule = Molecule(smiles='N#N', electronic_state='A' * 32, vibrational_level=2147483647)
    with pytest.raises(ValueError) as exc:
        LongFileGaussian(molecule, settings)
    filename = molecule.to_augmented_inchi_key() + '.log' + 'x' * 150
    assert_qm_filename_error(exc.value, molecule.to_augmented_inchi(), filename, 260, 255)
    assert not Path(settings.scratchDirectory).exists()
    assert not Path(settings.fileStore).exists()


@pytest.mark.parametrize('resolved,count,byte_count', [(True, 51, 256), (False, 78, 258)])
def test_direct_symmetry_writer_refuses_overlong_filename(settings, resolved, count, byte_count):
    molecule = radical_chain(count, resolved)
    key = molecule.to_augmented_inchi_key()
    Path(settings.scratchDirectory).mkdir(parents=True)
    with pytest.raises(ValueError) as exc:
        SymmetryJob(settings, key, nitrogen_qm_data()).write_input_file()
    assert_qm_filename_error(exc.value, key, key + '.symm', byte_count, 255)
    assert not list(Path(settings.scratchDirectory).iterdir())


@pytest.mark.parametrize('entry', ['geometry', 'qm'])
def test_qm_budget_is_rechecked_after_target_directory_creation(settings, monkeypatch, entry):
    molecule = Molecule(smiles='N#N', electronic_state='A' * 32, vibrational_level=2147483647)
    if entry == 'geometry':
        job = Geometry(settings, molecule.to_augmented_inchi_key(), molecule, molecule.to_augmented_inchi())
        rd_molecule, _ = job.rd_build()
    else:
        job = QMMolecule(molecule, settings)
        for directory in (settings.scratchDirectory, settings.fileStore):
            Path(directory).mkdir(parents=True)
    actual_pathconf = os.pathconf
    monkeypatch.setattr(os, 'pathconf', lambda path, name: 117 if Path(path) == Path(settings.scratchDirectory)
                        else actual_pathconf(path, name))
    with pytest.raises(ValueError) as exc:
        if entry == 'geometry':
            job.rd_embed(rd_molecule, 1)
        else:
            job.check_paths()
    assert_qm_filename_error(exc.value, molecule.to_augmented_inchi(), molecule.to_augmented_inchi_key() + '.refined.mol',
                             118, 117)
    for directory in (settings.scratchDirectory, settings.fileStore):
        assert not list(Path(directory).iterdir())


def test_mopac_temporary_filename_budget_is_checked_before_geometry(settings, monkeypatch):
    molecule = radical_chain(78)
    for directory in (settings.scratchDirectory, settings.fileStore):
        Path(directory).mkdir(parents=True)
    # The configured directories may allow more than the not-yet-created MOPAC temp directory.
    monkeypatch.setattr(os, 'pathconf', lambda path, name: 400)
    with pytest.raises(ValueError) as exc:
        MopacMolPM3(molecule, settings)
    assert_qm_filename_error(exc.value, molecule.to_augmented_inchi(), molecule.to_augmented_inchi_key() + '.mop', 257, 255)
    for directory in (settings.scratchDirectory, settings.fileStore):
        assert not list(Path(directory).iterdir())
