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

"""Structural database routing and review regressions for resolved identities."""
import ast
import copy
import io
from pathlib import Path

import pytest

from rmgpy.data.base import Database, Entry, ForbiddenStructures
from rmgpy.data.solvation import SoluteLibrary, SolventLibrary, SoluteData, SolventData, save_entry
from rmgpy.data.statmech import StatmechDepository, StatmechLibrary
from rmgpy.data.thermo import ThermoDepository
from rmgpy.data.transport import TransportLibrary
from rmgpy.exceptions import SpeciesIdentityError
from rmgpy.molecule import Molecule
from rmgpy.species import Species

ROOT = Path(__file__).resolve().parents[3]
HELPER = '_store_entry'


def entry_mapping(node, aliases=()):
    return (isinstance(node, ast.Attribute) and node.attr == 'entries'
            or isinstance(node, ast.Name) and (node.id == 'entries' or node.id in aliases))


def database_loader_violations(root, overrides=None):
    """Discover all entry loaders and writers; allow indexed stores only in the helper.

    Resets of the whole mapping and edits to Entry metadata do not replace an
    identity. No explicit loader list or census signature is used here.
    """
    sources = {str(path.relative_to(root)): path.read_text()
               for path in (root / 'rmgpy/data').rglob('*.py')}
    sources.update(overrides or {})
    failures = []
    for filename, source in sorted(sources.items()):
        tree = ast.parse(source)
        methods = []

        class Methods(ast.NodeVisitor):
            def __init__(self):
                self.owners = []

            def visit_ClassDef(self, node):
                self.owners.append(node.name)
                self.generic_visit(node)
                self.owners.pop()

            def visit_FunctionDef(self, node):
                methods.append(('.'.join(self.owners + [node.name]), node))
                self.owners.append(node.name)
                self.generic_visit(node)
                self.owners.pop()

            visit_AsyncFunctionDef = visit_FunctionDef

        Methods().visit(tree)
        for name, method in methods:
            key = filename + ':' + name
            stores, calls = [], []
            nodes = list(ast.walk(method))
            aliases = set()
            # Follow local bindings, including aliases of aliases, conservatively
            # over the method. Copying a mapping does not alias the original.
            changed = True
            while changed:
                before = len(aliases)
                for node in nodes:
                    if isinstance(node, ast.Assign):
                        targets, value = node.targets, node.value
                    elif isinstance(node, (ast.AnnAssign, ast.NamedExpr)):
                        targets, value = [node.target], node.value
                    else:
                        continue
                    if entry_mapping(value, aliases):
                        aliases.update(target.id for target in targets if isinstance(target, ast.Name))
                changed = len(aliases) != before
            for node in nodes:
                if isinstance(node, ast.Subscript) and isinstance(node.ctx, ast.Store) and entry_mapping(node.value, aliases):
                    stores.append(node)
                if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == HELPER:
                    calls.append(node)
                if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and entry_mapping(node.func.value, aliases) and node.func.attr in ('update', 'setdefault', '__setitem__'):
                    stores.append(node)
            if key == 'rmgpy/data/base.py:Database.' + HELPER:
                continue
            if stores:
                failures.append('entry assignment bypasses collision helper: ' + key)
            if (method.name == 'load_entry' or stores) and not calls:
                failures.append('database loader missing collision helper: ' + key)
    return failures


def states():
    ground = Species(label='N2', molecule=[Molecule(smiles='N#N')])
    a = copy.deepcopy(ground)
    a.molecule[0].electronic_state = 'A'
    a.molecule[0].vibrational_level = 1
    b = copy.deepcopy(a)
    b.molecule[0].electronic_state = 'B'
    return ground, a, b


def saved_records(cls, pair):
    stream = io.StringIO()
    db = cls()
    for index, spc in enumerate(pair, 1):
        if cls is SoluteLibrary:
            entry = Entry(index=index, label='N2', item=spc, data=SoluteData())
        elif cls is SolventLibrary:
            entry = Entry(index=index, label='N2', item=[spc], data=SolventData())
        else:
            entry = Entry(index=index, label='N2', item=spc.molecule[0])
        if cls is ForbiddenStructures:
            db.save_entry(stream, entry)
        else:
            save_entry(stream, entry)
    return stream.getvalue()


class TestExcitedExportRework5:
    @pytest.mark.parametrize('cls', [SoluteLibrary, SolventLibrary, ForbiddenStructures])
    @pytest.mark.parametrize('pair', [(1, 0), (0, 1), (1, 2)])
    def test_review_reload_collision(self, tmp_path, cls, pair):
        species = states()
        path = tmp_path / 'records.py'
        path.write_text(saved_records(cls, [species[i] for i in pair]))
        db = cls()
        with pytest.raises(SpeciesIdentityError, match='N2') as error:
            db.load(str(path))
        assert cls.__name__ + '.load_entry' in str(error.value)
        item = db.entries['N2'].item
        if isinstance(item, list):
            item = item[0]
        assert item.is_isomorphic(species[pair[0]] if isinstance(item, Species) else species[pair[0]].molecule[0])
        print('RELOAD REFUSED', cls.__name__, pair, str(error.value))

    @pytest.mark.parametrize('cls', [SoluteLibrary, SolventLibrary, ForbiddenStructures])
    def test_exact_resolved_repeat_and_distinct_labels(self, tmp_path, cls):
        spc = states()[1]
        path = tmp_path / 'records.py'
        source = saved_records(cls, [spc, spc])
        path.write_text(source)
        db = cls()
        db.load(str(path))
        assert len(db.entries) == 1
        path.write_text(source.replace('label = "N2"', 'label = "state"', 1))
        restored = cls()
        restored.load(str(path))
        assert len(restored.entries) == 2

    @pytest.mark.parametrize('cls,kind', [(ThermoDepository, 'thermo'), (StatmechDepository, 'statmech'),
                                        (StatmechLibrary, 'statmech'), (TransportLibrary, 'transport'), (Database, None)])
    def test_review_ground_collision_keeps_last(self, tmp_path, cls, kind):
        pair = [Molecule(smiles='N#N'), Molecule(smiles='O=O')]
        db = cls()
        if kind is None:
            db.load_old_dictionary('unused', pattern=False, content='\n'.join(mol.to_adjacency_list(label='ordinary') for mol in pair))
        else:
            for index, mol in enumerate(pair, 1):
                db.load_entry(index=index, label='ordinary', molecule=mol.to_adjacency_list(), **{kind: None})
        assert len(db.entries) == 1
        assert db.entries['ordinary'].item.is_isomorphic(pair[-1])
        assert db.entries['ordinary'].index == (-1 if kind is None else 2)
        print('GROUND LAST WINS', cls.__name__)

    @pytest.mark.parametrize('cls', [SoluteLibrary, SolventLibrary, ForbiddenStructures])
    def test_new_loader_ground_collision_keeps_last(self, tmp_path, cls):
        pair = [states()[0], Species(label='N2', molecule=[Molecule(smiles='O=O')])]
        path = tmp_path / 'records.py'
        path.write_text(saved_records(cls, pair))
        db = cls()
        db.load(str(path))
        item = db.entries['N2'].item
        if isinstance(item, list):
            item = item[0]
        assert item.is_isomorphic(pair[-1] if isinstance(item, Species) else pair[-1].molecule[0])

    def test_review_forbidden_species_public_save(self, tmp_path):
        spc = states()[1]
        db = ForbiddenStructures()
        db.load_entry(label='N2', species=spc.to_adjacency_list())
        path = tmp_path / 'forbidden.py'
        db.save(str(path))
        content = path.read_text()
        assert 'electronicstate A' in content and 'vibrationallevel 1' in content
        restored = ForbiddenStructures().load(str(path)).entries['N2'].item
        assert isinstance(restored, Species) and restored.is_isomorphic(spc)
        print('PUBLIC FORBIDDEN SAVE PRESERVES', restored.molecule[0].state_suffix())

    def test_ground_forbidden_species_save_unchanged(self):
        db = ForbiddenStructures()
        db.load_entry(label='N2', species=states()[0].to_adjacency_list())
        stream = io.StringIO()
        db.save_entry(stream, db.entries['N2'])
        assert '    group = "N2",\n' in stream.getvalue()
        assert 'species =' not in stream.getvalue()

    def test_every_database_loader_routes_collision_helper(self):
        assert database_loader_violations(ROOT) == []

    def test_removed_helper_call_is_rejected(self):
        filename = 'rmgpy/data/solvation.py'
        source = (ROOT / filename).read_text()
        tree = ast.parse(source)
        cls = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == 'SolventLibrary')
        method = next(node for node in cls.body if isinstance(node, ast.FunctionDef) and node.name == 'load_entry')
        changed = False
        for node in ast.walk(method):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == HELPER:
                node.func.attr = 'unprotected_store'
                changed = True
        assert changed, 'positive loader has no helper call to remove'
        failures = database_loader_violations(ROOT, {filename: ast.unparse(tree)})
        assert any('SolventLibrary.load_entry' in failure for failure in failures)
        print('REMOVED HELPER REJECTED', failures)

    @pytest.mark.parametrize('method', ['load_entry', 'read_records'])
    def test_new_unprotected_loader_is_rejected(self, method):
        filename = 'rmgpy/data/new_loader.py'
        source = 'class NewLoader:\n    def {}(self, label, entry):\n        self.entries[label] = entry\n'.format(method)
        failures = database_loader_violations(ROOT, {filename: source})
        assert any('NewLoader.' + method in failure for failure in failures)
        print('NEW LOADER REJECTED', failures)

    def test_call_does_not_excuse_an_unprotected_assignment(self):
        source = 'class NewLoader:\n    def load_entry(self, label, entry):\n        self._store_entry(label, entry)\n        self.entries[label] = entry\n'
        failures = database_loader_violations(ROOT, {'rmgpy/data/new_loader.py': source})
        assert any('entry assignment bypasses' in failure for failure in failures)

    def test_top_level_loader_definition_is_discovered(self):
        source = 'def load_entry(db, label, entry):\n    db.entries[label] = entry\n'
        failures = database_loader_violations(ROOT, {'rmgpy/data/new_loader.py': source})
        assert any('new_loader.py:load_entry' in failure for failure in failures)

    def test_discarded_invalid_ground_dictionary_record_keeps_legacy_behavior(self):
        db = Database()
        content = 'ordinary\nnot an adjacency list\n\n' + Molecule(smiles='O=O').to_adjacency_list(label='ordinary')
        db.load_old_dictionary('unused', pattern=False, content=content)
        assert db.entries['ordinary'].item.is_isomorphic(Molecule(smiles='O=O'))

    def test_forbidden_species_and_molecule_repeat_share_full_identity(self):
        spc = states()[1]
        db = ForbiddenStructures()
        db.load_entry(label='N2', species=spc.to_adjacency_list())
        db.load_entry(label='N2', molecule=spc.to_adjacency_list())
        assert db.entries['N2'].item.is_isomorphic(spc.molecule[0])
        with pytest.raises(SpeciesIdentityError, match='N2'):
            db.load_entry(label='N2', species=states()[2].to_adjacency_list())

    @pytest.mark.parametrize('kind', ['subscript', 'update', 'other'])
    @pytest.mark.parametrize('call_helper', [False, True])
    def test_review_alias_write_bypasses_token_helper_is_rejected(self, kind, call_helper):
        receiver = 'other' if kind == 'other' else 'self'
        write = 'records.update({label: entry})' if kind == 'update' else 'records[label] = entry'
        helper = '        self._store_entry(label, entry, "probe")\n' if call_helper else ''
        source = ('class NewLoader:\n    def load_entry(self, label, entry, other=None):\n'
                  + helper + '        records = ' + receiver + '.entries\n        ' + write + '\n')
        failures = database_loader_violations(ROOT, {'rmgpy/data/new_loader.py': source})
        assert any('entry assignment bypasses collision helper' in failure for failure in failures)
        print('REVIEW ALIAS REJECTED', kind, call_helper, failures)

    @pytest.mark.parametrize('binding', ['cache = entries', 'records = self.entries; cache = records', 'cache: dict = self.entries'])
    def test_local_entries_binding_write_is_rejected(self, binding):
        source = 'class NewLoader:\n    def read_records(self, label, entry, entries):\n        self._store_entry(label, entry, "probe")\n        ' + binding + '\n        cache[label] = entry\n'
        failures = database_loader_violations(ROOT, {'rmgpy/data/new_loader.py': source})
        assert any('entry assignment bypasses collision helper' in failure for failure in failures)

    @pytest.mark.parametrize('use', ['entry = records[label]', 'records[label].index = 1', 'snapshot = dict(records); snapshot[label] = entry'])
    def test_read_metadata_and_copied_mapping_aliases_are_allowed(self, use):
        source = 'class NewLoader:\n    def load_entry(self, label, entry):\n        records = self.entries\n        ' + use + '\n        self._store_entry(label, entry, "probe")\n'
        assert database_loader_violations(ROOT, {'rmgpy/data/new_loader.py': source}) == []
        print('SAFE ALIAS CONTROL ACCEPTED', use)

    @pytest.mark.parametrize('method,body', [
        ('load_entry', 'self._store_entry(label, entry, "probe")'),
        ('read_records', 'records = self.entries; return records[label]'),
        ('load_entry', 'records = self.entries; self._store_entry(label, entry, "probe"); return records[label]'),
    ])
    def test_required_safe_alias_controls(self, method, body):
        source = 'class NewLoader:\n    def ' + method + '(self, label, entry):\n        ' + body + '\n'
        assert database_loader_violations(ROOT, {'rmgpy/data/new_loader.py': source}) == []
        print('REQUIRED SAFE CONTROL ACCEPTED', method, body)
