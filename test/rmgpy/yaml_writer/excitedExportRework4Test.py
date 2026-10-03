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

"""Regressions for full-identity reloads, drawing refusal and ground coverage."""
import copy
import hashlib
import importlib.util
import io
from pathlib import Path
from unittest.mock import patch

import pytest

from rmgpy.data.base import Database, Entry
from rmgpy.data.thermo import ThermoDepository
from rmgpy.data.statmech import StatmechDepository, StatmechLibrary
from rmgpy.data.transport import TransportLibrary
from rmgpy.exceptions import SpeciesIdentityError
from rmgpy.kinetics import SurfaceArrhenius
from rmgpy.molecule import Molecule
from rmgpy.molecule.draw import MoleculeDrawer, ReactionDrawer
from rmgpy.reaction import Reaction
from rmgpy.species import Species

ROOT = Path(__file__).resolve().parents[3]
_spec = importlib.util.spec_from_file_location(
    'rework4_census', Path(__file__).with_name('excitedExportCensusTest.py'))
census = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(census)


def molecules():
    ground = Molecule(smiles='N#N')
    excited = copy.deepcopy(ground)
    excited.electronic_state = 'A'
    excited.vibrational_level = 1
    return ground, excited


LOADERS = [(ThermoDepository, 'thermo'), (StatmechDepository, 'statmech'),
           (StatmechLibrary, 'statmech'), (TransportLibrary, 'transport')]


class TestExcitedExportRework4:
    @pytest.mark.parametrize('reverse', [False, True])
    @pytest.mark.parametrize('trailing_blank', [False, True])
    def test_old_dictionary_collision(self, tmp_path, reverse, trailing_blank):
        pair = molecules()[::(-1 if reverse else 1)]
        # Separate adjacency records, including the final-record path without a blank line.
        content = '\n\n'.join(m.to_adjacency_list(label='N2').strip() for m in pair)
        path = tmp_path / 'dictionary.txt'
        path.write_text(content + ('\n\n' if trailing_blank else '\n'))
        db = Database()
        with pytest.raises(SpeciesIdentityError, match='N2') as error:
            db.load_old_dictionary(str(path), pattern=False)
        assert 'load_old_dictionary' in str(error.value)
        assert db.entries['N2'].item == pair[0].to_adjacency_list(label='N2')

    @pytest.mark.parametrize('cls,kind', LOADERS)
    @pytest.mark.parametrize('reverse', [False, True])
    @pytest.mark.parametrize('through_file', [False, True])
    def test_duplicate_reload_refuses(self, tmp_path, cls, kind, reverse, through_file):
        pair = molecules()[::(-1 if reverse else 1)]
        db = cls()
        if through_file:
            source = ''.join('entry(index={},label="N2",molecule={!r},{}=None)\n'.format(
                index, mol.to_adjacency_list(), kind) for index, mol in enumerate(pair, 1))
            path = tmp_path / 'records.py'
            path.write_text(source)
            load = lambda: db.load(str(path))
        else:
            db.load_entry(index=1, label='N2', molecule=pair[0].to_adjacency_list(), **{kind: None})
            load = lambda: db.load_entry(index=2, label='N2', molecule=pair[1].to_adjacency_list(), **{kind: None})
        with pytest.raises(SpeciesIdentityError, match='N2') as error:
            load()
        assert cls.__name__ + '.load_entry' in str(error.value)
        assert db.entries['N2'].item.is_isomorphic(pair[0])

    @pytest.mark.parametrize('cls,kind', LOADERS)
    def test_distinct_names_and_exact_repeats(self, cls, kind):
        ground, excited = molecules()
        db = cls()
        for index, label, molecule in [(1, 'N2', ground), (2, 'state', excited), (3, 'state', excited)]:
            db.load_entry(index=index, label=label, molecule=molecule.to_adjacency_list(), **{kind: None})
        assert len(db.entries) == 2 and db.entries['state'].index == 3
        other = copy.deepcopy(excited)
        other.electronic_state = 'B'
        with pytest.raises(SpeciesIdentityError, match='state'):
            db.load_entry(index=4, label='state', molecule=other.to_adjacency_list(), **{kind: None})

    def test_old_dictionary_exact_repeats_and_distinct_labels(self, tmp_path):
        ground, excited = molecules()
        content = '\n'.join(mol.to_adjacency_list(label=label) for label, mol in [
            ('N2', ground), ('state', excited), ('state', excited)])
        path = tmp_path / 'dictionary.txt'
        path.write_text(content)
        db = Database()
        db.load_old_dictionary(str(path), pattern=False)
        assert len(db.entries) == 2 and db.entries['state'].item.is_isomorphic(excited)

    @pytest.mark.parametrize('boundary', ['Molecule.draw', 'MoleculeDrawer.draw', 'ReactionDrawer.draw'])
    @pytest.mark.parametrize('backend_missing', [False, True])
    def test_drawing_refuses_before_output(self, tmp_path, boundary, backend_missing):
        ground, excited = molecules()
        target = tmp_path / 'image.png'
        target.write_bytes(b'existing image')
        before = target.read_bytes()
        if boundary == 'Molecule.draw':
            draw = lambda: excited.draw(str(target))
        elif boundary == 'MoleculeDrawer.draw':
            draw = lambda: MoleculeDrawer().draw(excited, 'png', str(target))
        else:
            reaction = Reaction(reactants=[Species(molecule=[ground])], products=[Species(molecule=[excited])])
            draw = lambda: ReactionDrawer().draw(reaction, 'png', str(target))
        with (patch('rmgpy.molecule.draw.cairo', None) if backend_missing else io.StringIO()):
            with pytest.raises(SpeciesIdentityError, match=boundary):
                draw()
        assert target.read_bytes() == before

    def test_ground_drawing_bytes(self, tmp_path):
        ground, _ = molecules()
        a, b = tmp_path / 'wrapper.png', tmp_path / 'drawer.png'
        ground.draw(str(a))
        MoleculeDrawer().draw(ground, 'png', str(b))
        assert a.read_bytes() == b.read_bytes()

    def test_ground_indexed_coverage_bytes(self):
        from rmgpy.data.kinetics.common import save_entry
        ground, _ = molecules()
        x = Species(label='X', index=3, molecule=[Molecule().from_adjacency_list('1 X u0 p0 c0')])
        a = Species(label='A', molecule=[ground])
        b = Species(label='B', molecule=[Molecule(smiles='[N]=[N]')])
        rate = SurfaceArrhenius(A=(1, 'm^3/(mol*s)'), Ea=(0, 'J/mol'),
                               coverage_dependence={x: {'a': 1, 'm': 0, 'E': (1, 'J/mol')}})
        entry = Entry(index=1, label='A + X <=> B + X',
                      item=Reaction(reactants=[a, x], products=[b, x]), data=rate)
        output = io.StringIO()
        save_entry(output, entry)
        assert output.getvalue().encode() == (ROOT / 'test/rmgpy/test_data/excited_export_ground_coverage.txt').read_bytes()

    def test_raw_bytes_and_snapshot_hashing(self, tmp_path):
        path = tmp_path / 'database.py'
        raw = b'name = "ground"\r\n'
        path.write_bytes(raw)
        db = ThermoDepository().load(str(path))
        assert db._loaded_file_sha256 == hashlib.sha256(raw).hexdigest()
        snapshot = 'name = "snapshot"\r\n'
        db.load(str(path), content=snapshot)
        assert db.name == 'snapshot'
        assert db._loaded_file_sha256 == hashlib.sha256(snapshot.encode()).hexdigest()

SINK_WRITERS = {
    'write': 'stream.write(payload)',
    'path-text': 'path.write_text(payload)',
    'path-bytes': 'path.write_bytes(payload)',
    'graph-dot': 'graph.write_dot(stream)',
    'graph-png': 'graph.write_png(stream)',
    'numpy-archive': 'np.savez_compressed(stream, payload)',
    'hdf5-dataset': 'stream.create_dataset("names", data=payload)',
    'hdf5-group': 'stream.create_group(payload)',
    'writelines': 'stream.writelines(payload)',
    'print-file': 'print(payload, file=stream)',
    'print-multiline': 'print(\n        payload,\n        file=stream,\n    )',
    'writerow': 'csv.writer(stream).writerow(payload)',
    'writerows': 'csv.writer(stream).writerows(payload)',
    'json-dump': 'json.dump(payload, stream)',
    'yaml-dump': 'yaml.dump(payload, stream)',
    'yaml-safe-dump': 'yaml.safe_dump(payload, stream)',
    'pickle-dump': 'pickle.dump(payload, stream)',
    'savefig': 'figure.savefig(stream)',
    'draw': 'drawer.draw(payload, "png", stream)',
    'cairo-png': 'surface.write_to_png(stream)',
    'cairo-pdf': 'cairo.PDFSurface(stream, 100, 100)',
    'cairo-svg': 'cairo.SVGSurface(stream, 100, 100)',
    'cairo-ps': 'cairo.PSSurface(stream, 100, 100)',
    'surface-finish': 'surface.finish()',
    'aliased-json': 'dump_json(payload, stream)',
    'module-output': 'stream.write(payload)',
    'class-output': 'stream.write(payload)',
    'cython-module': 'stream.write(payload)',
    'cython-nested-print': 'print(format_label(payload), file=stream)',
    'async-output': 'stream.write(payload)',
    'cython-output': 'stream.writelines(payload)',
}


class TestExcitedExportSinkDiscovery:
    @pytest.mark.parametrize('kind', SINK_WRITERS)
    def test_sink_without_naming_signal_requires_verdict(self, kind):
        statement = SINK_WRITERS[kind]
        imports = 'from json import dump as dump_json\n'
        if kind == 'module-output':
            source, suffix, target = imports + statement + '\n', '.py', '<module>'
        elif kind == 'class-output':
            source, suffix, target = 'class RawWriter:\n    ' + statement + '\n', '.py', 'RawWriter.<class>'
        elif kind == 'cython-module':
            source, suffix, target = statement + '\n', '.pyx', '<module>'
        elif kind in ('cython-output', 'cython-nested-print'):
            source, suffix, target = 'cpdef raw_save(object stream, object payload):\n    ' + statement + '\n', '.pyx', 'raw_save'
        else:
            prefix = 'async ' if kind == 'async-output' else ''
            source, suffix, target = imports + prefix + 'def raw_save(stream, payload):\n    ' + statement + '\n', '.py', 'raw_save'
        filename = 'rmgpy/new_sink' + suffix
        failures = census.violations(ROOT, census.POLICY, {filename: source})
        assert any(filename + ':' + target in failure for failure in failures), failures
        print('SINK', kind, 'REJECTED', [failure for failure in failures if filename in failure])

    def test_removed_verdict_is_rejected(self):
        policy = copy.deepcopy(census.POLICY)
        key = 'rmgpy/data/base.py:Database.load_old_dictionary'
        del policy['sites'][key]
        failures = census.violations(ROOT, policy)
        assert 'unrouted or changed candidate: ' + key in failures
        print('REMOVED VERDICT REJECTED:', key)
