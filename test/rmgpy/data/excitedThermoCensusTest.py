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

import ast
import csv
import collections
import json
import re
from pathlib import Path

FIELDS = {'thermo', 'thermo_data', 'H298', 'conformer', 'E0', 'Cp0', 'CpInf'}
CALLS = {'get_thermo_data', 'get_all_thermo_data', 'process_thermo_data',
         'generate_thermo_data', 'generate_thermo', 'calculate_cp0', 'calculate_cpinf',
         'find_cp0_and_cpinf', 'change_base_enthalpy', 'change_base_entropy',
         'to_nasa', 'to_wilhoit', 'get_statmech_data', 'get_site_solute_data', 'get_solvation_correction', 'load_yaml', 'make_object', 'checked_thermo', 'checked_energy', 'require_state_thermo', 'require_atom_thermo_allowed', 'require_species_thermo_allowed', 'require_network_thermo_allowed', 'require_thermo_estimation_allowed', 'require_reference_thermo_allowed'}
MODELS = {'NASA', 'NASAPolynomial', 'ThermoData', 'Wilhoit', 'Conformer'}


def census(root):
    sites = collections.Counter()
    examples = {}
    for package in ('rmgpy', 'arkane', 'scripts'):
        for path in sorted((root / package).rglob('*')):
            if path.suffix not in ('.py', '.pyx'):
                continue
            rel = path.relative_to(root).as_posix()
            source = path.read_text()
            scope = []
            def record(kind, spelling, line):
                key = (rel, '.'.join(scope) or '<module>', kind, spelling)
                sites[key] += 1
                examples.setdefault(key, source.splitlines()[line - 1].strip())
            if path.suffix == '.py':
                class Visitor(ast.NodeVisitor):
                    def visit_ClassDef(self, node):
                        scope.append(node.name)
                        self.generic_visit(node)
                        scope.pop()
                    def visit_FunctionDef(self, node):
                        scope.append(node.name)
                        self.generic_visit(node)
                        scope.pop()
                    visit_AsyncFunctionDef = visit_FunctionDef
                    def visit_Attribute(self, node):
                        if node.attr in FIELDS:
                            record('write' if isinstance(node.ctx, ast.Store) else 'read',
                                   ast.unparse(node), node.lineno)
                        self.generic_visit(node)
                    def visit_Call(self, node):
                        name = node.func.id if isinstance(node.func, ast.Name) else node.func.attr if isinstance(node.func, ast.Attribute) else ''
                        if name in CALLS or name in MODELS or ('thermo' in name and name.startswith(('get_', 'estimate_', 'compute_', 'correct_'))):
                            record('call', ast.unparse(node.func), node.lineno)
                        if name in {'getattr', 'setattr'} and len(node.args) > 1 and isinstance(node.args[1], ast.Constant) and node.args[1].value in FIELDS:
                            record('dynamic', ast.unparse(node.func) + ':' + str(node.args[1].value), node.lineno)
                        for keyword in node.keywords:
                            if keyword.arg in {'thermo', 'conformer'}:
                                record('argument-write', ast.unparse(node.func) + ':' + keyword.arg, node.lineno)
                        self.generic_visit(node)
                    def visit_Constant(self, node):
                        if isinstance(node.value, str) and '{{' in node.value:
                            for match in re.finditer(r'\b(?:\w+(?:\[[^]\n]+\])?\.)+(?:thermo|conformer|E0|Cp0|CpInf)\b', node.value):
                                record('template', match.group(), node.lineno)
                Visitor().visit(ast.parse(source))
            else:
                frames = []
                triple = None
                for line_number, line in enumerate(source.splitlines(), 1):
                    stripped = line.lstrip()
                    indent = len(line) - len(stripped)
                    if triple:
                        if triple in stripped:
                            triple = None
                        continue
                    if stripped.startswith(('\"\"\"', "'" * 3)):
                        marker = stripped[:3]
                        if stripped.count(marker) == 1:
                            triple = marker
                        continue
                    if not stripped or stripped.startswith('#'):
                        continue
                    while frames and indent <= frames[-1][0]:
                        frames.pop()
                    declaration = re.match(r'(?:cdef\s+)?class\s+(\w+)|(?:def|cpdef|cdef)\s+(?:[\w.]+\s+)*(\w+)\s*\(|property\s+(\w+)\s*:', stripped)
                    if declaration:
                        frames.append((indent, next(group for group in declaration.groups() if group)))
                    scope[:] = [name for _, name in frames]
                    code = stripped.split('#', 1)[0]
                    for match in re.finditer(r'\b(?:\w+(?:\[[^]\n]+\])?\.)+(?:thermo|conformer|E0|Cp0|CpInf)\b', code):
                        record('cython-field', match.group(), line_number)
                    for match in re.finditer(r'\b([\w.]+)\s*\(', code):
                        name = match.group(1).rsplit('.', 1)[-1]
                        if name in CALLS or name in MODELS or ('thermo' in name and name.startswith(('get_', 'estimate_', 'compute_', 'correct_'))):
                            record('cython-call', match.group(1), line_number)
    return sites, examples




ROOT = Path(__file__).resolve().parents[3]
INVENTORY = Path(__file__).with_name('excitedThermoSites.csv')


def classified_sites():
    with INVENTORY.open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    return {(row['path'], row['scope'], row['kind'], row['expression']): int(row['count'])
            for row in rows}


def assert_classified(root, expected):
    actual, _ = census(root)
    assert actual == expected, ('Unclassified or changed thermo acquisition/consumer sites',
                                actual - expected, expected - actual)


def test_every_thermo_site_has_a_census_verdict():
    assert_classified(ROOT, collections.Counter(classified_sites()))


def test_new_python_consumer_and_acquisition_fail_census(tmp_path):
    import pytest
    package = tmp_path / 'rmgpy'
    package.mkdir()
    (package / 'new.py').write_text('def export(sp):\n    return sp.thermo\ndef fit(sp):\n    sp.thermo = NASA()\n')
    with pytest.raises(AssertionError, match='Unclassified'):
        assert_classified(tmp_path, collections.Counter())


def test_new_cython_indexed_energy_reader_fails_census(tmp_path):
    import pytest
    package = tmp_path / 'rmgpy'
    package.mkdir()
    (package / 'new.pyx').write_text('def energy(species_dict, label):\n    return species_dict[label].conformer.E0\n')
    actual, _ = census(tmp_path)
    assert any('conformer.E0' in key[3] for key in actual)
    with pytest.raises(AssertionError, match='Unclassified'):
        assert_classified(tmp_path, collections.Counter())


def test_new_site_in_already_classified_function_fails_census(tmp_path):
    import pytest
    package = tmp_path / 'rmgpy'
    package.mkdir()
    path = package / 'old.py'
    path.write_text('def export(sp):\n    return sp.thermo\n')
    expected, _ = census(tmp_path)
    path.write_text('def export(sp):\n    print(sp.thermo)\n    return sp.thermo\n')
    with pytest.raises(AssertionError, match='Unclassified'):
        assert_classified(tmp_path, expected)


def test_each_census_row_names_an_existing_behavioral_test():
    rows = json.loads(Path(__file__).with_name('excitedThermoRows.json').read_text())
    with INVENTORY.open(newline='') as stream:
        assigned = {row['row'] for row in csv.DictReader(stream)}
    assert assigned == set(rows)
    for name, row in rows.items():
        assert row['verdict'] and row['sites']
        path, test = row['test'].split(':')
        tree = ast.parse((ROOT / path).read_text())
        assert any(isinstance(node, ast.FunctionDef) and node.name == test for node in ast.walk(tree)), name
