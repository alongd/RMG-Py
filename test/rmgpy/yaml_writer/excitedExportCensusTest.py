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

"""Conservative source census; unknown naming or reconstruction candidates fail closed."""
import ast
import hashlib
import importlib.util
import io
import tokenize
import re
from pathlib import Path


NAME_CALLS = r'\b(?:Database\.load|_store_entry|Database\.load|_store_entry|write_files_atomically|save_chemkin|write_yml|to_chemkin|get_species_identifier|to_cantera|save_dictionary|save_species_dictionary|load_species_dictionary|read_generation_pair|read_generation_files|stamp_file_generation|resolve_species_reference|SpeciesReferences)\s*\('
REBUILD = r'\b(?:set_structure|set_structure|from_smiles|from_adjacency_list|make_object|recursive_make_object)\s*\('
EMIT = r'\b(?:write|dump|Node|Edge|GenericData|savefig|set_label|set_title)\s*\('


# Conservatively discover sinks even when labels are prepared in another function.
# Attribute calls are tracked by sink, independent of the receiver's variable name.
SINK_METHODS = {
    'write_text', 'write_bytes', 'write_dot', 'write_png', 'savez_compressed', 'create_dataset', 'create_group',
    'write', 'writelines', 'writerow', 'writerows', 'dump', 'safe_dump',
    'savefig', 'draw', 'write_to_png', 'write_to_pdf', 'write_to_svg',
    'write_to_ps', 'finish', 'PDFSurface', 'SVGSurface', 'PSSurface',
}
SINK_PATTERN = r'\b(?:' + '|'.join(sorted(SINK_METHODS)) + r')\s*\('


def sink_calls(tree, aliases):
    """Find output calls, including imported aliases and multiline print(file=)."""
    result = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if isinstance(node.func, ast.Attribute):
            name = node.func.attr
        elif isinstance(node.func, ast.Name):
            name = aliases.get(node.func.id, node.func.id).rsplit('.', 1)[-1]
        else:
            continue
        if name in SINK_METHODS or (name == 'print' and
                any(keyword.arg == 'file' for keyword in node.keywords)):
            result.append((node.lineno, name))
    return result



def cython_sink_calls(source):
    """Tokenize Cython output calls without treating comments/strings as code."""
    tokens = [token for token in tokenize.generate_tokens(io.StringIO(source).readline)
              if token.type not in (tokenize.COMMENT, tokenize.STRING, tokenize.NL,
                                    tokenize.NEWLINE, tokenize.INDENT, tokenize.DEDENT)]
    result = []
    for i, token in enumerate(tokens[:-1]):
        if token.type != tokenize.NAME or tokens[i + 1].string != '(':
            continue
        if token.string in SINK_METHODS:
            result.append((token.start[0], token.string))
        elif token.string == 'print':
            depth = 0
            for j in range(i + 1, len(tokens) - 1):
                value = tokens[j].string
                if value in ('(', '[', '{'):
                    depth += 1
                elif value in (')', ']', '}'):
                    depth -= 1
                if depth == 0:
                    break
                if depth == 1 and value == 'file' and tokens[j + 1].string == '=':
                    result.append((token.start[0], 'print'))
                    break
    return result

def discover(root, overrides=None):
    """Find naming, reconstruction and all output sinks in Python/Cython sources."""
    overrides = overrides or {}
    result = {}
    files = {str(path.relative_to(root)): path.read_text() for package in ('rmgpy', 'arkane')
             for path in (root / package).rglob('*') if path.suffix in ('.py', '.pyx')}
    files.update(overrides)
    for filename, source in sorted(files.items()):
        owner = ''
        lines = source.splitlines()
        extents, functions, aliases = {}, {}, {}
        if filename.endswith('.py'):
            tree = ast.parse(source)
            for node in ast.walk(tree):
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    extents[node.lineno - 1] = node.end_lineno
                    functions[node.lineno - 1] = node
                elif isinstance(node, ast.ImportFrom):
                    for imported in node.names:
                        aliases[imported.asname or imported.name] = (node.module or '') + '.' + imported.name
            # Module-level exports also require classification. Function bodies
            # have their own signatures and must not hide unguarded top-level IO.
            statements = [node for node in tree.body if not isinstance(
                node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))]
            module = ast.Module(body=statements, type_ignores=[])
            sinks = sink_calls(module, aliases)
            if sinks:
                body = '\n'.join('\n'.join(lines[node.lineno - 1:node.end_lineno]) for node in statements)
                key = filename + ':<module>'
                result[key] = {'file': filename, 'function': '<module>', 'line': 1,
                               'signals': sinks, 'signature': hashlib.sha256(body.encode()).hexdigest(),
                               'body': body}
            for node in ast.walk(tree):
                if not isinstance(node, ast.ClassDef):
                    continue
                statements = [child for child in node.body if not isinstance(
                    child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))]
                sinks = sink_calls(ast.Module(body=statements, type_ignores=[]), aliases)
                if sinks:
                    body = '\n'.join('\n'.join(lines[child.lineno - 1:child.end_lineno]) for child in statements)
                    key = filename + ':' + node.name + '.<class>'
                    result[key] = {'file': filename, 'function': node.name + '.<class>', 'line': node.lineno,
                                   'signals': sinks, 'signature': hashlib.sha256(body.encode()).hexdigest(),
                                   'body': body}
        function_lines = set()
        for i, line in enumerate(lines):
            match = re.match(r'(?:cdef )?class (\w+)', line)
            if match:
                owner = match.group(1)
                continue
            match = re.match(r'(\s*)(?:async\s+def|def|cpdef|cdef)\s+(?:[\w.]+\s+)?(\w+)\(', line)
            if not match:
                continue
            indent, name = len(match.group(1)), match.group(2)
            end = i + 1
            while end < len(lines) and (not lines[end].strip() or
                    len(lines[end]) - len(lines[end].lstrip()) > indent):
                end += 1
            end = extents.get(i, end)
            body = '\n'.join(lines[i:end])
            function_lines.update(range(i, end))
            if i in functions:
                signals = sink_calls(functions[i], aliases)
            else:
                signals = [(i + line, kind) for line, kind in cython_sink_calls(body)]
            for j in range(i, end):
                statement = lines[j].strip()
                if statement.startswith('#'):
                    continue
                if (re.search(NAME_CALLS, statement) or re.search(REBUILD, statement)
                        or (('dictionary' in statement or 'reactions.py' in statement)
                            and re.search(r'\b(?:open|save|load|get_species)\(', statement))
                        or (re.search(r'\.label|\bstr\((?:spc\w*|spec\w*|rxn\w*|react\w*|product\w*)\)', statement)
                            and (re.search(EMIT, body) or name.startswith(('to_', 'render_', '__repr__', '__str__', 'get_'))))
                        or ('to_smiles(' in statement and re.search(EMIT, body))
                        or ('refuse_resolved_' in statement)):
                    signals.append((j + 1, statement))
            if name.endswith('_wall_manifest'):
                signals.append((i + 1, 'serializable wall manifest API'))
            if not signals:
                continue
            function = owner + '.' + name if indent and owner else name
            key = filename + ':' + function
            signature = hashlib.sha256(body.encode()).hexdigest()
            result[key] = {'file': filename, 'function': function, 'line': i + 1,
                           'signals': signals, 'signature': signature, 'body': body}
        if filename.endswith('.pyx'):
            # Cython import-time writes have no enclosing def/cpdef/cdef.
            body = '\n'.join('' if i in function_lines else line for i, line in enumerate(lines))
            sinks = cython_sink_calls(body)
            if sinks:
                key = filename + ':<module>'
                result[key] = {'file': filename, 'function': '<module>', 'line': 1,
                               'signals': sinks, 'signature': hashlib.sha256(body.encode()).hexdigest(),
                               'body': body}
    return result

def audited_function(root, key, overrides=None):
    """Read a Python function by qualified owner, without importing backends."""
    filename, scope = key.split(':')
    source = (overrides or {}).get(filename, (root / filename).read_text())
    tree = ast.parse(source)
    nodes = tree.body
    for part in scope.split('.'):
        node = next(n for n in nodes if isinstance(n, (ast.ClassDef, ast.FunctionDef)) and n.name == part)
        nodes = node.body
    body = '\n'.join(source.splitlines()[node.lineno - 1:node.end_lineno])
    return node, hashlib.sha256(body.encode()).hexdigest()


def lexical_nodes(node):
    """Visit lexical bindings and executable bodies without nested definitions."""
    children = list(ast.iter_child_nodes(node))
    if isinstance(node, ast.If) and isinstance(node.test, ast.Constant):
        children = [node.test, *(node.body if node.test.value else node.orelse)]
    for child in children:
        yield child
        if not isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef, ast.Lambda)):
            yield from lexical_nodes(child)
        if isinstance(child, (ast.Return, ast.Raise)):
            break


def valid_refusal_chain(root, key, decision, overrides=None):
    """Verify every explicitly audited call hop and the named terminal guard."""
    chain = decision.get('refusal_chain', [])
    if not chain:
        return False
    caller = key
    try:
        for hop in chain:
            node, _ = audited_function(root, caller, overrides)
            executed = list(lexical_nodes(node))
            calls = {ast.unparse(n) for n in executed if isinstance(n, ast.Call)}
            if not hop.get('invocations') or not set(hop['invocations']) <= calls:
                return False
            guard_name = hop['call'].split('.')[0]
            stores = [n for n in executed if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Store)
                      and n.id == guard_name]
            if stores and not hop.get('receiver_binding'):
                return False
            if any(isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
                   and n.name == guard_name for n in executed):
                return False
            if any(isinstance(n, (ast.Import, ast.ImportFrom))
                   and any((a.asname or a.name.split('.')[0]) == guard_name for a in n.names)
                   and not isinstance(n, ast.ImportFrom) for n in executed):
                return False
            if any(arg.arg == guard_name for arg in node.args.args + node.args.kwonlyargs):
                return False
            if hop.get('receiver_binding'):
                bindings = {ast.unparse(n) for n in executed if isinstance(n, ast.Assign)}
                if hop['receiver_binding'] not in bindings or len(stores) != 1:
                    return False
            target, signature = audited_function(root, hop['callee'], overrides)
            if signature != hop['signature']:
                return False
            # Imported guard names must resolve to the exact audited module/function.
            if '.' not in hop['call']:
                filename = caller.split(':')[0]
                source = (overrides or {}).get(filename, (root / filename).read_text())
                module = ast.parse(source)
                module_nodes = [n for statement in module.body
                                if not isinstance(statement, (ast.FunctionDef, ast.ClassDef))
                                for n in [statement, *lexical_nodes(statement)]]
                imports = {alias.asname or alias.name: (n.module or '') + '.' + alias.name
                           for n in module_nodes + executed if isinstance(n, ast.ImportFrom)
                           for alias in n.names}
                for definition in module.body:
                    if isinstance(definition, ast.FunctionDef):
                        imports[definition.name] = filename[:-3].replace('/', '.') + '.' + definition.name
                if any(isinstance(n, ast.Name) and isinstance(n.ctx, ast.Store) and n.id == guard_name
                       for n in module_nodes):
                    return False
                target_file, target_scope = hop['callee'].split(':')
                expected = target_file[:-3].replace('/', '.') + '.' + target_scope
                if imports.get(hop['call']) != expected:
                    return False
            caller = hop['callee']
        return any(isinstance(n, ast.Raise) and isinstance(n.exc, ast.Call)
                   and isinstance(n.exc.func, ast.Name) and n.exc.func.id == 'ExcitedSpeciesThermoError'
                   for n in ast.walk(target))
    except (KeyError, StopIteration, SyntaxError):
        return False


def violations(root, policy, overrides=None):
    """Require explicit R/N site verdicts, tests and immutable audited exclusions."""
    found = discover(root, overrides)
    known = {**policy['sites'], **policy['search_exclusions']}
    failures = []
    for key, candidate in found.items():
        decision = known.get(key)
        if decision is None or decision['signature'] != candidate['signature']:
            failures.append('unrouted or changed candidate: ' + key)
            continue
        if key in policy['sites']:
            manifold = decision.get('manifold', {})
            if (manifold.get('verdict') not in ('carried', 'refused', 'not-applicable')
                    or not manifold.get('reason') or not manifold.get('test')
                    or not manifold.get('evidence')):
                failures.append('missing manifold export classification: ' + key)
            if decision['verdict'] not in ('R', 'N') or not decision.get('test'):
                failures.append('missing R/N verdict or test: ' + key)
            if (decision['verdict'] == 'N' and 'refuse_resolved_' not in candidate['body']
                    and 'StateProvenanceError' not in candidate['body']
                    and 'ResolvedStateTrainingError' not in candidate['body']
                    and not valid_refusal_chain(root, key, decision, overrides)):
                failures.append('missing named refusal: ' + key)
        elif decision.get('verdict') != 'E' or not decision.get('reason'):
            failures.append('unexplained search exclusion: ' + key)
    failures.extend('missing audited site: ' + key for key in known if key not in found)
    return failures

import ast
import copy
import io
import json
import textwrap
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from rmgpy.exceptions import ExcitedSpeciesThermoError, SpeciesIdentityError, StateProvenanceError, VibrationalManifoldError
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.kinetics import Arrhenius
from rmgpy.data.base import Entry, Database
from rmgpy.export import SpeciesReferences, resolve_species_reference
# This helper is outside the installed packages; load its sibling file directly.
_helpers_spec = importlib.util.spec_from_file_location(
    'excited_export_helpers', Path(__file__).with_name('excitedExportHelpers.py'),
)
_helpers = importlib.util.module_from_spec(_helpers_spec)
_helpers_spec.loader.exec_module(_helpers)
nitrogen = _helpers.nitrogen
database = _helpers.database

ROOT = Path(__file__).resolve().parents[3]
POLICY = json.loads((ROOT / 'test/rmgpy/test_data/excited_export_census.json').read_text())
FOUND = discover(ROOT)
N_SITES = [key for key, site in POLICY['sites'].items() if site['verdict'] == 'N']
R_SITES = [key for key, site in POLICY['sites'].items() if site['verdict'] == 'R']


def refusal_witness(key, tmp_path, spc=None):
    """Run the actual boundary body, isolating optional native plotting/QM backends.

    The unchanged function body and signature are compiled from the audited
    source. Only decorators/defaults and optional-backend globals are replaced.
    A refusal must precede all backend work and filesystem changes.
    """
    if spc is None:
        spc = nitrogen(level=1)
    if key in {
        'arkane/common.py:ArkaneSpecies.save_yaml',
        'arkane/common.py:ArkaneSpecies.load_yaml',
        'arkane/encorr/reference.py:ReferenceSpecies.__init__',
        'arkane/encorr/reference.py:ReferenceSpecies.load_yaml',
        'arkane/encorr/reference.py:ReferenceSpecies.to_error_canceling_spcs',
        'arkane/encorr/reference.py:ReferenceDatabase.load',
        'rmgpy/qm/molecule.py:QMMolecule.save_thermo_data',
        'arkane/statmech.py:StatMechJob.load',
    }:
        import yaml
        from arkane.common import ArkaneSpecies
        from arkane.encorr.reference import ReferenceSpecies, ReferenceDatabase
        from arkane.statmech import StatMechJob
        from rmgpy.qm.molecule import QMMolecule
        path = tmp_path / '0.yml'
        payload = {'class': 'ReferenceSpecies', 'label': 'A', 'adjacency_list': spc.to_adjacency_list()}
        if key == 'arkane/common.py:ArkaneSpecies.load_yaml':
            payload.update({'class': 'ArkaneSpecies', 'is_ts': False})
        if key.endswith('.load_yaml') or key.endswith(':ReferenceDatabase.load'):
            path.write_text(yaml.safe_dump(payload))
        before = {p.relative_to(tmp_path): p.read_bytes() if p.is_file() else None
                  for p in tmp_path.rglob('*')}
        carrier = SimpleNamespace(adjacency_list=spc.to_adjacency_list(), molecule=spc.molecule[0])
        calls = {
            'arkane/common.py:ArkaneSpecies.save_yaml': lambda: ArkaneSpecies.save_yaml(carrier, str(tmp_path)),
            'arkane/common.py:ArkaneSpecies.load_yaml': lambda: ArkaneSpecies.__new__(ArkaneSpecies).load_yaml(str(path)),
            'arkane/encorr/reference.py:ReferenceSpecies.__init__': lambda: ReferenceSpecies(adjacency_list=spc.to_adjacency_list()),
            'arkane/encorr/reference.py:ReferenceSpecies.load_yaml': lambda: ReferenceSpecies.__new__(ReferenceSpecies).load_yaml(str(path)),
            'arkane/encorr/reference.py:ReferenceSpecies.to_error_canceling_spcs': lambda: ReferenceSpecies.to_error_canceling_spcs(carrier, None),
            'arkane/encorr/reference.py:ReferenceDatabase.load': lambda: ReferenceDatabase().load(paths=[str(tmp_path)], ignore_incomplete=False),
            'rmgpy/qm/molecule.py:QMMolecule.save_thermo_data': lambda: QMMolecule.save_thermo_data(carrier),
            'arkane/statmech.py:StatMechJob.load': lambda: StatMechJob.load(SimpleNamespace(species=spc)),
        }
        error = (VibrationalManifoldError if spc.props.get('vibrational_manifold') and key in {
            'arkane/statmech.py:StatMechJob.load', 'rmgpy/qm/molecule.py:QMMolecule.save_thermo_data'}
            else StateProvenanceError if key == 'arkane/statmech.py:StatMechJob.load'
            else ExcitedSpeciesThermoError)
        with pytest.raises(error):
            calls[key]()
        assert {p.relative_to(tmp_path): p.read_bytes() if p.is_file() else None
                for p in tmp_path.rglob('*')} == before
        return
    if key == 'rmgpy/qm/molecule.py:QMMolecule.load_thermo_data':
        from rmgpy.qm.molecule import QMMolecule
        job = object.__new__(QMMolecule)
        job.molecule = spc.molecule[0]
        error = VibrationalManifoldError if spc.props.get('vibrational_manifold') else StateProvenanceError
        with pytest.raises(error):
            job.load_thermo_data()
        assert not list(tmp_path.iterdir())
        return
    if key in {'rmgpy/qm/gaussian.py:GaussianMol.write_input_file',
               'rmgpy/qm/mopac.py:MopacMol.write_input_file'}:
        from rmgpy.qm.gaussian import GaussianMol
        from rmgpy.qm.mopac import MopacMol
        cls = GaussianMol if 'gaussian.py' in key else MopacMol
        job = object.__new__(cls)
        job.molecule = spc.molecule[0]
        before = list(tmp_path.iterdir())
        error = SpeciesIdentityError if spc.props.get('vibrational_manifold') else StateProvenanceError
        with pytest.raises(error):
            job.write_input_file(1)
        assert list(tmp_path.iterdir()) == before
        return
    if key in {
        'rmgpy/molecule/molecule.py:Molecule.draw',
        'rmgpy/molecule/molecule.py:Molecule._repr_png_',
        'rmgpy/molecule/draw.py:MoleculeDrawer.draw',
        'rmgpy/molecule/draw.py:ReactionDrawer.draw',
        'rmgpy/reaction.py:Reaction.draw',
        'rmgpy/reaction.py:Reaction._repr_png_',
    }:
        from rmgpy.molecule.draw import MoleculeDrawer, ReactionDrawer
        reaction = Reaction(reactants=[spc], products=[spc])
        calls = {
            'Molecule.draw': lambda: spc.molecule[0].draw(str(tmp_path / 'image.png')),
            'Molecule._repr_png_': spc.molecule[0]._repr_png_,
            'MoleculeDrawer.draw': lambda: MoleculeDrawer().draw(spc.molecule[0], 'png', str(tmp_path / 'image.png')),
            'ReactionDrawer.draw': lambda: ReactionDrawer().draw(reaction, 'png', str(tmp_path / 'image.png')),
            'Reaction.draw': lambda: reaction.draw(str(tmp_path / 'image.png')),
            'Reaction._repr_png_': reaction._repr_png_,
        }
        with pytest.raises(SpeciesIdentityError, match=key.split(':')[1]):
            calls[key.split(':')[1]]()
        assert not list(tmp_path.iterdir())
        return
    if key in {'arkane/kinetics.py:KineticsDrawer._get_label_size',
               'arkane/kinetics.py:KineticsDrawer._draw_label'}:
        from arkane.kinetics import KineticsDrawer
        spc.conformer.E0 = (0, 'kJ/mol')
        well = SimpleNamespace(species_list=[spc])
        drawer = KineticsDrawer()
        method = key.rsplit('.', 1)[1]
        args = (well,) if method == '_get_label_size' else (well, None, 0, 0)
        with pytest.raises(SpeciesIdentityError, match=method):
            getattr(drawer, method)(*args)
        assert not list(tmp_path.iterdir())
        return
    if key == 'rmgpy/util.py:strip_yaml_notes':
        import yaml
        from rmgpy.util import strip_yaml_notes
        source = tmp_path / 'annotated.yaml'
        destination = tmp_path / 'stripped.yaml'
        source.write_text(yaml.safe_dump({'species': [{'name': 'N2(v1)', 'note': spc.to_adjacency_list()}]}))
        original = source.read_bytes()
        with pytest.raises(SpeciesIdentityError, match='strip_yaml_notes'):
            strip_yaml_notes(str(source), str(destination))
        assert source.read_bytes() == original and not destination.exists()
        return
    rxn = Reaction(reactants=[spc], products=[spc], kinetics=Arrhenius(A=(1, 's^-1'),Ea=(0,'J/mol')))
    network = SimpleNamespace(
        get_all_species=lambda: [spc],
        path_reactions=[rxn],
        net_reactions=[],
        check_resolved_species_reversibility=lambda: rxn.check_resolved_species_reversibility(),
    )
    model = SimpleNamespace(species=[spc], reactions=[rxn])
    coreedge = SimpleNamespace(core=model, edge=model, output_species_list=[spc], species=[spc], reactions=[rxn])
    reference = SimpleNamespace(adjacency_list=spc.to_adjacency_list(), molecule=spc.molecule[0])
    carrier = SimpleNamespace(species=spc, molecule=spc.molecule[0], reaction=rxn, network=network,
        job=SimpleNamespace(reaction=rxn, network=network), entries={1: Entry(item=spc)},
        species_list=[spc], surface_species_list=[], sensitive_species=[spc], extra_species=[],
        initial_species=[spc], reaction_systems=[], reaction_model=coreedge,
        usedTST=False,
        observables={'species': [spc]}, output_species_list=[spc], coverage_dependence={spc: {}},
        reactor_mod=SimpleNamespace(output_species_list=[spc], cantera=SimpleNamespace(species_list=[spc], reaction_list=[rxn])),
        spc=reference, bac=SimpleNamespace(dataset=[SimpleNamespace(spc=reference)]),
        y_var=[SimpleNamespace(species=spc)], adjacency_list=spc.to_adjacency_list(),
        isodesmicReactionList=[], verbose=20, input_file='unused', output_directory=str(tmp_path),
        load_input_file=lambda path: [SimpleNamespace(species=spc)],
        _check_data_serialization=lambda entry: None)
    if POLICY['sites'][key]['file'].endswith('.py'):
        module_tree=ast.parse((ROOT/POLICY['sites'][key]['file']).read_text())
        original=next(node for node in ast.walk(module_tree) if isinstance(node,ast.FunctionDef)
                      and node.lineno==FOUND[key]['line'])
        tree=ast.Module(body=[original],type_ignores=[])
    else:
        tree=ast.parse(textwrap.dedent(FOUND[key]['body']))
    fn = tree.body[0]
    fn.decorator_list = []
    fn.returns = None
    for arg in fn.args.args + fn.args.kwonlyargs:
        arg.annotation = None
    fn.args.defaults = []
    fn.args.kw_defaults = [None] * len(fn.args.kwonlyargs)
    globals_ = {
        'Molecule': Molecule, 'Species': Species, 'Reaction': Reaction, 'Entry': Entry,
        'os': __import__('os'), 'logging': __import__('logging'),
        'kinetics_references': __import__('rmgpy.export', fromlist=['kinetics_references']).kinetics_references,
        'describe_species': __import__('rmgpy.export', fromlist=['describe_species']).describe_species,
        'check_reaction_rate_serialization': (
            lambda reaction, kinetics, reversible:
            reaction.check_resolved_species_reversibility(kinetics=kinetics, reversible=reversible)),
        'SpeciesIdentityError': SpeciesIdentityError,
        'initialize_log': lambda *a, **k: None, 'log_header': lambda *a, **k: None,
        'get_models_to_merge': lambda *a, **k: [coreedge],
        'combine_models': lambda *a, **k: coreedge,
    }
    exec(compile(ast.fix_missing_locations(tree), str(ROOT / POLICY['sites'][key]['file']), 'exec'), globals_)
    args = {
        'self': carrier, 'species': [spc], 'context': 'census refusal', 'reactions': [rxn],
        'species_list': [spc], 'rxn_list': [rxn], 'sensitive_species': [spc],
        'reaction_model': coreedge, 'common_species_list': [(spc, spc)],
        'species_list1': [], 'species_list2': [], 'common_reactions': [(rxn, rxn)],
        'unique_reactions1': [], 'unique_reactions2': [], 'core_reactions': [rxn],
        'configuration': SimpleNamespace(species=[spc]), 'network': network,
        'channel': SimpleNamespace(species=[spc]), 'reaction': rxn, 'rmg': carrier,
        'part_core_edge': 'core',
        'obj': spc, 'rmg_species': [spc], 'spcs': [spc], 'input_files': [],
        'input_model_files': [],
        'wd': str(tmp_path), 'transport': False, 'path': str(tmp_path / 'output'),
        'output_directory': str(tmp_path), 'f': io.StringIO(), 'entry': Entry(item=spc),
    }
    actual = {arg.arg: args.get(arg.arg) for arg in fn.args.args + fn.args.kwonlyargs}
    if fn.args.kwarg:
        actual.update(wd=str(tmp_path), transport=False)
    expected = ((VibrationalManifoldError if spc.props.get('vibrational_manifold') else ExcitedSpeciesThermoError)
                if key == 'arkane/statmech.py:StatMechJob.write_output' else SpeciesIdentityError)
    message = 'vibrationalManifold' if expected is VibrationalManifoldError else '(?i)resolved.*(?:supported|serialize|retain|export|library)'
    with pytest.raises(expected, match=message) as error:
        globals_[fn.name](**actual)
    if expected is SpeciesIdentityError:
        assert fn.name in str(error.value) or 'census refusal' in str(error.value) or 'Arkane ThermoJob Chemkin' in str(error.value)
    assert not list(tmp_path.iterdir())


def resolved_witness(key, tmp_path):
    spc, other = nitrogen('A', 1), nitrogen('B')
    with _helpers.unit_thermo_library([spc, other]):
        _resolved_witness(key, tmp_path, spc, other)


def _resolved_witness(key, tmp_path, spc, other):
    """Exercise the site's shared output/reload family with real resolved species.

    The manifest records the concrete site and the associated family witness.
    Existing reference-position tests independently cover undeclared references.
    """
    from rmgpy import chemkin, yaml_rms, yaml_cantera1, yaml_cantera2
    from rmgpy.data.kinetics.library import KineticsLibrary
    from rmgpy.data.kinetics.depository import KineticsDepository
    from rmgpy.rmg.model import CoreEdgeReactionModel
    rxn = Reaction(reactants=[spc], products=[other], kinetics=Arrhenius(A=(1, 's^-1'),Ea=(0,'J/mol')))
    file, function = key.split(':')
    if key in {
        'rmgpy/solver/plasma.pyx:PlasmaReactor._compute_reference_reaction_data',
        'rmgpy/solver/plasma.pyx:PlasmaReactor._evaluate_electronegative_wall_regime',
        'rmgpy/solver/plasma.pyx:PlasmaReactor.electronegative_wall_manifest',
    }:
        from importlib import util
        witness_path = Path(__file__).with_name('excitedExportRework10Test.py')
        spec = util.spec_from_file_location('excited_export_rework10_witness', witness_path)
        witness = util.module_from_spec(spec)
        spec.loader.exec_module(witness)
        witness.test_resolved_wall_manifest_uses_state_qualified_species_names()
        return
    if key in {
        'rmgpy/data/base.py:Database._store_entry',
        'rmgpy/data/base.py:ForbiddenStructures.load_entry',
        'rmgpy/data/base.py:ForbiddenStructures.save_entry',
    }:
        from rmgpy.data.base import ForbiddenStructures
        loaded = ForbiddenStructures()
        loaded.load_entry(label='N2', species=spc.to_adjacency_list())
        path = tmp_path / 'forbidden.py'
        loaded.save(str(path))
        restored = ForbiddenStructures().load(str(path)).entries['N2'].item
        assert isinstance(restored, Species) and restored.is_isomorphic(spc)
        with pytest.raises(SpeciesIdentityError, match='N2'):
            loaded.load_entry(label='N2', species=nitrogen().to_adjacency_list())
        assert loaded.entries['N2'].item.is_isomorphic(spc)
        return
    if key in {
        'rmgpy/data/solvation.py:SoluteLibrary.load_entry',
        'rmgpy/data/solvation.py:SoluteLibrary.load',
        'rmgpy/data/solvation.py:SolventLibrary.load_entry',
        'rmgpy/data/solvation.py:SolventLibrary.load',
    }:
        from rmgpy.data.solvation import SoluteLibrary, SolventLibrary, SoluteData, SolventData, save_entry
        cls = SolventLibrary if 'SolventLibrary' in function else SoluteLibrary
        loaded = cls()
        data = io.StringIO()
        for index, species in enumerate((spc, nitrogen()), 1):
            entry = Entry(index=index, label='N2', item=[species] if cls is SolventLibrary else species,
                          data=SolventData() if cls is SolventLibrary else SoluteData())
            save_entry(data, entry)
        path = tmp_path / 'records.py'
        path.write_text(data.getvalue())
        with pytest.raises(SpeciesIdentityError, match='N2'):
            loaded.load(str(path))
        item = loaded.entries['N2'].item
        assert (item[0] if isinstance(item, list) else item).is_isomorphic(spc)
        return
    if key in {'rmgpy/data/base.py:Database.load_old_dictionary', 'rmgpy/data/base.py:Database.check_record'}:
        path = tmp_path / 'dictionary.txt'
        ground = nitrogen()
        path.write_text(ground.molecule[0].to_adjacency_list(label='N2') + '\n' + spc.molecule[0].to_adjacency_list(label='N2') + '\n')
        with pytest.raises(SpeciesIdentityError, match='N2'):
            Database().load_old_dictionary(str(path), pattern=False)
        path.write_text(ground.molecule[0].to_adjacency_list(label='N2') + '\n' + spc.molecule[0].to_adjacency_list(label='state') + '\n')
        loaded = Database()
        loaded.load_old_dictionary(str(path), pattern=False)
        assert len(loaded.entries) == 2 and loaded.entries['state'].item.has_resolved_state()
        return
    if key in {
        'rmgpy/data/thermo.py:ThermoDepository.load_entry',
        'rmgpy/data/statmech.py:StatmechDepository.load_entry',
        'rmgpy/data/statmech.py:StatmechLibrary.load_entry',
        'rmgpy/data/transport.py:TransportLibrary.load_entry',
    }:
        import importlib
        cls = getattr(importlib.import_module(file[:-3].replace('/', '.')), function.split('.')[0])
        kind = 'thermo' if 'Thermo' in function else 'statmech' if 'Statmech' in function else 'transport'
        loaded = cls()
        loaded.load_entry(index=1, label='N2', molecule=spc.to_adjacency_list(), **{kind: None})
        with pytest.raises(SpeciesIdentityError, match='N2'):
            loaded.load_entry(index=2, label='N2', molecule=nitrogen().to_adjacency_list(), **{kind: None})
        assert loaded.entries['N2'].item.has_resolved_state()
        return
    if key in {'rmgpy/data/base.py:Database.save', 'rmgpy/data/base.py:Database.render_save'}:
        from rmgpy.data.thermo import ThermoDepository
        loaded = ThermoDepository()
        loaded.load_entry(index=1, label='state', molecule=spc.to_adjacency_list(), thermo=spc.thermo)
        content = loaded.render_save()
        assert 'vibrationallevel 1' in content
        path = tmp_path / 'records.py'
        Database.save(loaded, str(path))
        restored = ThermoDepository().load(str(path), local_context={'NASA': __import__('rmgpy.thermo', fromlist=['NASA']).NASA, 'NASAPolynomial': __import__('rmgpy.thermo', fromlist=['NASAPolynomial']).NASAPolynomial})
        assert restored.entries['state'].item.has_resolved_state()
        return
    if key == 'rmgpy/data/kinetics/family.py:complete_round_trip':
        from rmgpy.data.kinetics.family import complete_round_trip
        restored = complete_round_trip(spc)
        assert restored.is_isomorphic(spc) and restored.molecule[0].has_resolved_state()
        return
    if key in {'rmgpy/qm/molecule.py:Geometry.rd_embed',
               'rmgpy/qm/gaussian.py:GaussianMol.write_input_file',
               'rmgpy/qm/mopac.py:MopacMol.write_input_file'}:
        from rmgpy.qm.main import QMSettings
        from rmgpy.qm.molecule import QMMolecule
        from rmgpy.qm.gaussian import GaussianMolPM3
        from rmgpy.qm.mopac import MopacMolPM3
        settings = QMSettings(fileStore=str(tmp_path / 'store'), scratchDirectory=str(tmp_path / 'scratch'))
        calculator = GaussianMolPM3 if 'gaussian.py' in key else MopacMolPM3 if 'mopac.py' in key else QMMolecule
        job = calculator(spc.molecule[0], settings)
        geometry = job.create_geometry()
        paths = [Path(geometry.get_crude_mol_file_path()), Path(geometry.get_refined_mol_file_path())]
        if calculator is not QMMolecule:
            job.write_input_file(1)
            paths.append(Path(job.input_file_path))
        for path in paths:
            title = next(line for line in path.read_text().splitlines() if line.startswith('InChI='))
            restored = Molecule().from_augmented_inchi(title)
            assert restored.is_isomorphic(spc.molecule[0]) and restored.has_resolved_state()
            assert path.name.startswith(spc.molecule[0].to_augmented_inchi_key())
        return
    if key == 'rmgpy/qm/molecule.py:QMMolecule.save_thermo_data':
        from rmgpy.qm.molecule import QMMolecule
        from rmgpy.thermo import ThermoData
        thermo = ThermoData(H298=(0, 'J/mol'), S298=(0, 'J/(mol*K)'),
                            Tdata=([300, 400], 'K'), Cpdata=([30, 30], 'J/(mol*K)'))
        path = tmp_path / 'thermo.py'
        carrier = SimpleNamespace(molecule=spc.molecule[0], unique_id_long=spc.molecule[0].to_augmented_inchi(),
                                  thermo=thermo, point_group=None, qm_data=None,
                                  check_file_names=lambda: None, get_thermo_file_path=lambda: str(path))
        QMMolecule.save_thermo_data(carrier)
        namespace = {'ThermoData': ThermoData}
        exec(path.read_text(), namespace)
        assert Molecule().from_adjacency_list(namespace['adjacencyList']).is_isomorphic(spc.molecule[0])
        assert '/v:1' in namespace['InChI']
        return
    if function == 'RMG.make_seed_mech':
        from rmgpy.rmg.main import RMG
        rmg=RMG(output_directory=str(tmp_path))
        rmg.reaction_model=CoreEdgeReactionModel()
        rmg.reaction_model.core.species=[spc,other]
        rmg.reaction_model.core.reactions=[rxn]
        rmg.save_seed_to_database=False
        rmg.save_seed_modulus=-1
        (tmp_path/'seed').mkdir()
        for rate in (1,2):
            rxn.kinetics.A.value_si=rate
            rmg.make_seed_mech()
            loaded=KineticsLibrary()
            loaded.load(str(tmp_path/'seed/seed/reactions.py'),local_context={'Arrhenius':Arrhenius})
            entry=next(iter(loaded.entries.values()))
            assert entry.item.reactants[0].molecule[0].has_resolved_state()
            assert entry.data.A.value_si==rate
        return
    if function == 'RMG.load_database':
        from rmgpy.rmg.main import RMG
        from rmgpy.data.rmg import RMGDatabase
        import rmgpy.data.rmg as database_module
        from rmgpy.data.kinetics.database import KineticsDatabase
        db=database(KineticsLibrary,1,1)
        db.save(str(tmp_path/'reactions.py'))
        def load_stub(self,**kwargs):
            loaded=KineticsLibrary()
            loaded.load(str(tmp_path/'reactions.py'),local_context={'Arrhenius':Arrhenius})
            self.kinetics=KineticsDatabase()
            self.kinetics.libraries={'probe':loaded}
            self.thermo=SimpleNamespace(adsorption_groups=[])
        rmg=RMG()
        rmg.adsorption_groups=[]
        rmg.reaction_libraries=[]
        rmg.kinetics_depositories=['!training']
        rmg.forbidden_structures=[]
        with patch.object(RMGDatabase,'load',load_stub),patch.object(rmg,'check_libraries'),patch.object(database_module,'database'):
            rmg.load_database()
        assert next(iter(rmg.database.kinetics.libraries['probe'].entries.values())).item.products[0].molecule[0].has_resolved_state()
        return
    if key == 'rmgpy/solver/base.pyx:ReactionSystem.simulate':
        from rmgpy.solver.simple import SimpleReactor
        from rmgpy.solver.termination import TerminationTime
        from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
        path=tmp_path/'sensitivity.csv'
        reactor=SimpleReactor(T=(500,'K'),P=(1,'bar'),initial_mole_fractions={spc:1},
            termination=[TerminationTime((1e-7,'s'))],sensitive_species=[spc],sensitivity_threshold=-1)
        reactor.simulate([spc,other],[rxn],[],[],[],[],sensitivity=True,sens_worksheet=[str(path)],
            model_settings=ModelSettings(tol_keep_in_edge=0,tol_move_to_core=1,tol_interrupt_simulation=0),
            simulator_settings=SimulatorSettings())
        assert '(v1)' in path.read_text().splitlines()[0]
        return
    if key == 'rmgpy/qm/molecule.py:QMMolecule.load_thermo_data':
        from rmgpy.qm import molecule as qm
        from rmgpy.thermo import ThermoData
        carrier=SimpleNamespace(get_thermo_file_path=lambda:'fixture',unique_id_long='identity',molecule=spc.molecule[0])
        payload={'InChI':'identity','adjacencyList':spc.to_adjacency_list(),
            'thermoData':ThermoData(),'pointGroup':SimpleNamespace(point_group='Dinfh'),'qmData':None}
        with patch.object(qm,'load_thermo_data_file',return_value=payload):
            assert qm.QMMolecule.load_thermo_data(carrier) is payload['thermoData']
            payload['adjacencyList']=other.to_adjacency_list()
            assert qm.QMMolecule.load_thermo_data(carrier) is None
        return
    if file == 'rmgpy/tools/canteramodel.py':
        from rmgpy.tools.canteramodel import Cantera,CanteraCondition
        simulation=Cantera(species_list=[spc,other],reaction_list=[rxn],output_directory=str(tmp_path),
            conditions=[CanteraCondition('IdealGasReactor',(1e-7,'s'),{spc:1},T0=(500,'K'),P0=(1,'bar'))])
        simulation.load_model()
        assert simulation.model.species_names[0] != simulation.model.species_names[1]
        if function=='Cantera.modify_species_thermo':
            simulation.modify_species_thermo(0,spc,use_chemkin_identifier=True)
        elif function=='Cantera.modify_reaction_kinetics':
            simulation.modify_reaction_kinetics(0,rxn)
        elif function=='Cantera.simulate':
            results=simulation.simulate()
            assert any(data.species is spc for data in results[0][1])
        assert '(v1)' in simulation.model.species_names[0]
        return
    if key == 'arkane/encorr/data.py:BACDatapoint._mol_from_adjlist':
        from arkane.encorr.data import BACDatapoint
        carrier = SimpleNamespace(spc=SimpleNamespace(adjacency_list=spc.to_adjacency_list()))
        BACDatapoint._mol_from_adjlist(carrier)
        assert carrier._mol.has_resolved_state()
        return
    if key == 'arkane/encorr/reference.py:ReferenceDatabase.get_species_from_label':
        from arkane.encorr.reference import ReferenceSpecies, ReferenceDatabase
        ref = ReferenceSpecies(label='A', smiles='N#N')
        ref.adjacency_list = spc.to_adjacency_list()
        db = ReferenceDatabase()
        db.reference_sets = {'probe': [ref]}
        assert db.get_species_from_label(['A'], set_name='probe')[0] is ref
        assert Molecule().from_adjacency_list(ref.adjacency_list).has_resolved_state()
        return
    if function == 'RMG.generate_end_of_run_cantera_files':
        from rmgpy.rmg.main import RMG
        rmg=RMG(output_directory=str(tmp_path))
        rmg.reaction_model=CoreEdgeReactionModel()
        rmg.reaction_model.core.species=[spc,other]
        rmg.reaction_model.core.reactions=[rxn]
        rmg.export_failures=[]
        from rmgpy.rmg.settings import WriterConfig
        rmg.chemkin_writer_config=SimpleNamespace(enabled=True)
        rmg.generate_end_of_run_cantera_files()
        assert any(isinstance(error,SpeciesIdentityError) for _,error in rmg.export_failures)
    elif file.startswith('rmgpy/chemkin') or file == 'rmgpy/util.py' or function == 'RMG.make_seed_mech':
        model = CoreEdgeReactionModel(); model.core.species=[spc, other]; model.core.reactions=[rxn]
        chemkin.save_chemkin(model, str(tmp_path/'chem.inp'), str(tmp_path/'annotated.inp'), str(tmp_path/'dictionary.txt'))
        species, reactions = chemkin.load_chemkin_file(str(tmp_path/'chem.inp'),str(tmp_path/'dictionary.txt'))
        assert any(m.molecule[0].has_resolved_state() for m in species)
        assert reactions[0].reactants[0].molecule[0].has_resolved_state()
        assert '(v1)' in chemkin.get_species_identifier(spc)
    elif file.startswith('rmgpy/data/kinetics') or file == 'rmgpy/data/base.py' or function == 'RMG.load_database' or key == 'arkane/main.py:Arkane.get_libraries':
        for cls in (KineticsLibrary, KineticsDepository):
            path=tmp_path/cls.__name__;path.mkdir()
            db=database(cls,1,1);db.save(str(path/'reactions.py'))
            loaded=cls();loaded.load(str(path/'reactions.py'),local_context={'Arrhenius':Arrhenius})
            assert next(iter(loaded.entries.values())).item.products[0].molecule[0].has_resolved_state()
    elif file == 'rmgpy/yaml_rms.py':
        if function == 'write_rms':
            path = tmp_path / 'states.rms'
            yaml_rms.write_rms([spc, other], [rxn], path=str(path))
            assert len(yaml_rms.load_rms_species(str(path))) == 2
        result=yaml_rms.get_mech_dict([spc,other],[rxn])
        assert 'vibrationallevel 1' in str(result).lower()
        # Saved adjacency reload must retain the state even if SMILES is added.
        data={'name':'A(v1)','adjlist':spc.to_adjacency_list(),'smiles':'N#N'}
        import yaml
        path=tmp_path/'states.rms';path.write_text(yaml.safe_dump({'Phases':[{'Species':[data]}]}))
        assert yaml_rms.load_rms_species(str(path))[0].molecule[0].has_resolved_state()
    elif key == 'rmgpy/yaml_cantera2.py:_write_file_atomically':
        from rmgpy.yaml_cantera2 import _write_file_atomically
        path = tmp_path / 'identity.txt'
        _write_file_atomically(str(path), spc.to_adjacency_list())
        assert Molecule().from_adjacency_list(path.read_text()).is_isomorphic(spc.molecule[0])
    elif key == 'rmgpy/yaml_cantera2.py:save_cantera_model':
        model = CoreEdgeReactionModel()
        model.core.species = [spc, other]
        model.core.reactions = [rxn]
        path = tmp_path / 'cantera.yaml'
        yaml_cantera2.save_cantera_model(model.core, str(path))
        assert 'vibrationallevel 1' in path.read_text().lower()
    elif file.startswith('rmgpy/yaml_cantera'):
        module=yaml_cantera1 if file.endswith('1.py') else yaml_cantera2
        path=tmp_path/'cantera.yaml'
        if module is yaml_cantera1:
            module.write_cantera([spc,other],[rxn],set(),path=str(path))
        else:
            module.generate_cantera_data([spc,other],[rxn],set())
        result=path.read_text() if path.exists() else module.generate_cantera_data([spc,other],[rxn],set())
        assert 'vibrationallevel 1' in str(result).lower()
    elif file.startswith('arkane/common') or function == 'StatMechJob.load':
        from arkane.common import ArkaneSpecies
        from arkane.statmech import StatMechJob
        ArkaneSpecies(species=spc).save_yaml(str(tmp_path))
        path=next((tmp_path/'species').glob('*.yml'))
        consumer=nitrogen('A');StatMechJob(consumer,str(path)).load()
        assert consumer.molecule[0].has_resolved_state()
    elif file.startswith('arkane/encorr/reference') or file.startswith('arkane/encorr/data'):
        from arkane.encorr.reference import ReferenceSpecies
        ref=ReferenceSpecies(label='A',adjacency_list=spc.to_adjacency_list(),smiles='N#N')
        assert Molecule().from_adjacency_list(ref.adjacency_list).has_resolved_state()
        if function.endswith('to_error_canceling_spcs'):
            # Conversion carries the same adjacency-backed molecule.
            assert Molecule().from_adjacency_list(ref.adjacency_list).state_suffix() == '|v:1'
    elif file in ('rmgpy/species.py','rmgpy/reaction.py','rmgpy/tools/canteramodel.py','rmgpy/kinetics/model.pyx'):
        native=[s.to_cantera(use_chemkin_identifier=True,all_species=[spc,other]) for s in (spc,other)]
        assert native[0].name != native[1].name
        converted=rxn.to_cantera(species_list=[spc,other],use_chemkin_identifier=True)
        assert '(v1)' in str(converted.equation)
        restored=eval(repr(spc), {'Species':Species,'Molecule':Molecule, 'Conformer':__import__('rmgpy.statmech',fromlist=['Conformer']).Conformer, 'NASA':__import__('rmgpy.thermo',fromlist=['NASA']).NASA,'NASAPolynomial':__import__('rmgpy.thermo',fromlist=['NASAPolynomial']).NASAPolynomial})
        assert restored.molecule[0].has_resolved_state()
    elif file.startswith('rmgpy/data/'):
        import importlib
        module=importlib.import_module(file[:-3].replace('/','.'))
        data=io.StringIO();module.save_entry(data,Entry(index=1,label='A',item=spc if file.endswith('solvation.py') else spc.molecule[0]))
        assert 'vibrationallevel 1' in data.getvalue().lower()
        for name in ('ThermoLibrary','ThermoDepository','StatmechLibrary','StatmechDepository','TransportLibrary','SoluteLibrary','ForbiddenStructures'):
            cls=getattr(module,name,None)
            if cls and function.startswith(name+'.load_entry'):
                db=cls()
                kwargs={'thermo':None} if 'Thermo' in name else ({'statmech':None} if 'Statmech' in name else ({'transport':None} if 'Transport' in name else {'solute':None}))
                db.load_entry(index=1,label='A',molecule=spc.to_adjacency_list(),**kwargs)
                item=next(iter(db.entries.values())).item
                assert (item.molecule[0] if isinstance(item,Species) else item).has_resolved_state()
    elif file == 'rmgpy/export.py':
        declarations=SpeciesReferences([spc,other],identifiers=chemkin.get_species_identifier)
        assert resolve_species_reference(spc,declarations) != resolve_species_reference(other,declarations)
        from rmgpy.export import validate_reaction_references, refuse_resolved_job
        validate_reaction_references([rxn],declarations)
        if function == 'refuse_resolved_job':
            with pytest.raises(SpeciesIdentityError, match='census job'):
                refuse_resolved_job(SimpleNamespace(species=spc),'census job')
    elif file == 'rmgpy/rmg/listener.py':
        from rmgpy.rmg.listener import SimulationProfileWriter
        (tmp_path/'solver').mkdir()
        SimulationProfileWriter(str(tmp_path),0,[spc,other]).update(SimpleNamespace(snapshots=[[0,1,0.5,0.5]]))
        assert '(v1)' in next((tmp_path/'solver').glob('*.csv')).read_text()
    elif file == 'rmgpy/tools/fluxdiagram.py':
        # Actual Graphviz/video output is covered by review finding 7.
        monkeypatch=pytest.MonkeyPatch()
        try:
            _helpers.FluxDiagramWitness().test_review_7_flux_nodes_keep_identity(tmp_path,monkeypatch)
        finally:
            monkeypatch.undo()
    elif file.startswith('rmgpy/molecule/') or file.startswith('rmgpy/rmgobject') or file.startswith('rmgpy/thermo/') or file.startswith('rmgpy/statmech/') or file == 'rmgpy/pdep/network.py':
        import pickle
        restored=pickle.loads(pickle.dumps(spc.molecule[0]))
        assert restored.has_resolved_state()
        saved=spc.molecule[0].to_adjacency_list()
        assert Molecule().from_adjacency_list(saved).state_suffix() == '|v:1'
        assert 'vibrational_level=1' in repr(spc.molecule[0]) or 'vibrationallevel 1' in repr(spc.molecule[0])
    elif key == 'rmgpy/rmg/input.py:core_species_file':
        chemkin.save_species_dictionary(str(tmp_path/'dictionary.txt'),[spc])
        assert next(iter(chemkin.load_species_dictionary(str(tmp_path/'dictionary.txt')).values())).molecule[0].has_resolved_state()
    elif key == 'rmgpy/qm/molecule.py:QMMolecule.load_thermo_data':
        assert Molecule().from_adjacency_list(spc.to_adjacency_list()).has_resolved_state()
    elif key == 'rmgpy/solver/base.pyx:ReactionSystem.simulate':
        declarations=SpeciesReferences([spc,other],identifiers=chemkin.get_species_identifier)
        assert len(set(declarations.names)) == 2
    else:
        raise AssertionError('No resolved witness family registered for '+key)


class TestExcitedExportSiteCensus:
    def test_complete_static_census(self):
        assert not violations(ROOT,POLICY)

    def test_wall_manifest_with_resolved_species_is_written_after_identity_resolution(self, tmp_path):
        from rmgpy.rmg.listener import SimulationProfileWriter

        solver = tmp_path / 'solver'
        solver.mkdir()
        record = {
            'gates': 'not-evaluated-zero-anion',
            'h': 1.0,
            'floating_potential_e_over_kTe': None,
        }
        reaction_system = SimpleNamespace(
            electronegative_wall_model='confinedAnion',
            electronegative_wall_history={0.0: record},
            electronegative_wall_manifest=lambda: {
                **record,
                'closure': 'confinedAnion',
                'geometry_arm': 'fullFrequency',
            },
            snapshots=[[0.0, 1.0, 1.0]],
        )
        species = nitrogen(level=1)
        identity = __import__('rmgpy.chemkin', fromlist=['get_species_identifier']).get_species_identifier(species)
        reaction_system.electronegative_wall_manifest = lambda: {
            **record,
            'closure': 'confinedAnion',
            'geometry_arm': 'fullFrequency',
            'gates': {'A': {identity: {}}},
        }
        writer = SimulationProfileWriter(str(tmp_path), 0, [species])
        writer.update(reaction_system)

        assert identity in next(solver.glob('*.csv')).read_text().splitlines()[0]
        assert identity in json.loads(next(solver.glob('*.json')).read_text())['gates']['A']

    def test_wall_manifest_discovery_fails_closed_for_sibling_method(self):
        filename = 'rmgpy/solver/plasma.pyx'
        source = (ROOT / filename).read_text() + '\n\ndef sibling_wall_manifest(self):\n    return {"species": self.species.label}\n'

        assert any('sibling_wall_manifest' in failure for failure in
                   violations(ROOT, POLICY, overrides={filename: source}))

    def test_profile_writer_preserves_ground_argon_identifier_collision(self, tmp_path):
        from rmgpy.rmg.listener import SimulationProfileWriter

        solver = tmp_path / 'solver'
        solver.mkdir()
        argon = Species(label='Ar', molecule=[Molecule(smiles='[Ar]')])
        argon_ion = Species(label='Ar+', molecule=[Molecule(smiles='[Ar+]')])
        writer = SimulationProfileWriter(str(tmp_path), 0, [argon, argon_ion])
        writer.update(SimpleNamespace(snapshots=[[0.0, 1.0, 0.75, 0.25]]))

        profile = next(solver.glob('*.csv')).read_text().splitlines()
        assert profile[0] == 'Time (s),Volume (m^3),Ar,Ar'
        assert profile[1] == '0.0,1.0,0.75,0.25'

    @pytest.mark.parametrize('mutation',['new-file','existing-function'])
    def test_injected_unrouted_site_is_rejected(self,mutation):
        file='rmgpy/new_unrouted_export.py' if mutation=='new-file' else 'rmgpy/tools/fluxdiagram.py'
        source='' if mutation=='new-file' else (ROOT/file).read_text()
        if mutation=='new-file':
            source+='\ndef unrouted_save(f, species):\n    f.write(species.label)\n'
            target='unrouted_save'
        else:
            source=source.replace('    declarations = SpeciesReferences(',
                '    output.write(species.label)\n    declarations = SpeciesReferences(',1)
            target='generate_flux_diagram'
        failures=violations(ROOT,POLICY,{file:source})
        assert any(target in failure for failure in failures)
        print('INJECTED UNROUTED SITE REJECTED:',failures)

    def test_refusal_dependency_changes_fail_closed(self):
        filename = 'rmgpy/thermo/state.py'
        source = (ROOT / filename).read_text().replace('raise ExcitedSpeciesThermoError(', 'raise ValueError(')
        key = 'arkane/common.py:ArkaneSpecies.save_yaml'
        assert not valid_refusal_chain(ROOT, key, POLICY['sites'][key], {filename: source})
        assert any(key in failure for failure in violations(ROOT, POLICY, {filename: source}))

    def test_deleted_guard_call_fails_even_with_a_refreshed_boundary_signature(self):
        policy = copy.deepcopy(POLICY)
        key = 'arkane/common.py:ArkaneSpecies.save_yaml'
        filename = key.split(':')[0]
        source = (ROOT / filename).read_text().replace('        require_reference_thermo_allowed(self)', '        pass', 1)
        overrides = {filename: source}
        policy['sites'][key]['signature'] = discover(ROOT, overrides)[key]['signature']
        assert any(key in failure for failure in violations(ROOT, policy, overrides))

    @pytest.mark.parametrize('replacement', [
        '        require_reference_thermo_allowed(None)',
        '        def unused():\n            require_reference_thermo_allowed(self)',
        '        require_reference_thermo_allowed = lambda owner: None\n        require_reference_thermo_allowed(self)',
        '        def require_reference_thermo_allowed(owner):\n            pass\n        require_reference_thermo_allowed(self)',
        '        if False:\n            require_reference_thermo_allowed(self)',
        '        return\n        require_reference_thermo_allowed(self)',
    ])
    def test_wrong_owner_or_unexecuted_guard_fails_with_refreshed_signature(self, replacement):
        key = 'arkane/common.py:ArkaneSpecies.save_yaml'
        filename = key.split(':')[0]
        source = (ROOT / filename).read_text().replace('        require_reference_thermo_allowed(self)', replacement, 1)
        overrides = {filename: source}
        policy = copy.deepcopy(POLICY)
        policy['sites'][key]['signature'] = discover(ROOT, overrides)[key]['signature']
        assert any(key in failure for failure in violations(ROOT, policy, overrides))

    @pytest.mark.parametrize('key',N_SITES,ids=N_SITES)
    def test_named_refusal(self,key,tmp_path):
        refusal_witness(key,tmp_path)

    @pytest.mark.parametrize('key',R_SITES,ids=R_SITES)
    def test_resolved_route(self,key,tmp_path):
        resolved_witness(key,tmp_path)


@pytest.mark.parametrize('field', ['manifold', 'verdict', 'reason', 'evidence', 'test'])
def test_manifold_census_rejects_missing_classification(field):
    policy = copy.deepcopy(POLICY)
    key = 'rmgpy/chemkin.pyx:render_species_dictionary'
    if field == 'manifold':
        del policy['sites'][key][field]
    else:
        del policy['sites'][key]['manifold'][field]
    assert 'missing manifold export classification: ' + key in violations(ROOT, policy)
