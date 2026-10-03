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

"""Enumerate scientific data suppliers without importing RMG.

Every Python/Cython function in rmgpy and Arkane is a candidate. This conservative
superset covers database wrappers, QM/ML, species/model assignment, energy
transfer, numeric models, ESS file loaders, and their transitive helpers without
relying on a potentially incomplete dynamic call graph. Reviewed classifications
are explicit records; discovery never assigns a verdict.
"""

import argparse
import ast
from collections import Counter
import hashlib
import json
from pathlib import Path

VERDICTS = {'carried', 'excluded', 'refused', 'not-reachable-for-resolved'}


def collect_sites(repository):
    """Return every Python/Cython method, callback and lambda in both packages."""
    repository = Path(repository)
    sites = {}
    for directory in (repository / 'rmgpy', repository / 'arkane'):
        for path in sorted(list(directory.rglob('*.py')) + list(directory.rglob('*.pyx'))):
            relative = path.relative_to(repository).as_posix()
            if path.suffix == '.pyx':
                sites.update(collect_cython_sites(path, relative))
                continue
            tree = ast.parse(path.read_text(), filename=relative)
            lambda_counts = Counter()
            occurrence_counts = Counter()

            def visit(node, scope=()):
                for child in ast.iter_child_nodes(node):
                    if isinstance(child, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)):
                        if isinstance(child, ast.Lambda):
                            lambda_counts[scope] += 1
                            name = '<lambda{}>'.format(lambda_counts[scope])
                        else:
                            name = child.name
                            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
                                for decorator in child.decorator_list:
                                    if isinstance(decorator, ast.Attribute) and decorator.attr in ('setter', 'deleter'):
                                        name += '[{}]'.format(decorator.attr)
                        qualified = scope + (name,)
                        if not isinstance(child, ast.ClassDef):
                            key = '{}:{}'.format(relative, '.'.join(qualified))
                            occurrence_counts[key] += 1
                            if occurrence_counts[key] > 1:
                                key += '[{}]'.format(occurrence_counts[key])
                            digest = hashlib.sha256(ast.dump(child, include_attributes=False).encode()).hexdigest()
                            sites[key] = {'file': relative, 'line': child.lineno, 'source_sha256': digest}
                        visit(child, qualified)
                    else:
                        visit(child, scope)

            visit(tree)
    return sites


def collect_cython_sites(path, relative):
    """Use Cython's parser for def/cdef/cpdef, properties and nested callbacks.

    For .pyx sources, hashes include the companion pxd and included pxi contents.
    Cythonized .py sources use the Python AST collector; their declarations and
    build parity require separate review. Parsing performs no compilation
    and imports no application modules. Failures propagate, rather than omit a
    supplier. All sources in this repository use local pxi includes.
    """
    from Cython.Compiler import Errors, Main
    from Cython.Compiler.Scanning import FileSourceDescriptor

    Errors.init_thread()
    context = Main.Context([], {}, options=Main.CompilationOptions(Main.default_options))
    module = relative.replace('/', '.').rsplit('.', 1)[0]
    source = FileSourceDescriptor(str(path.resolve()), relative)
    scope = context.find_module(module, need_pxd=False)
    tree = context.parse(source, scope, False, module)
    lines = path.read_text().splitlines(keepends=True)
    declarations = path.with_suffix('.pxd')
    dependencies = declarations.read_text() if declarations.exists() else ''
    # The parser resolves includes; hash them as dependencies too.
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('include '):
            included = path.parent / stripped.split('"')[1]
            dependencies += included.read_text()
    sites = {}
    counts = Counter()

    def visit(node, scope=()):
        kind = type(node).__name__
        if kind in ('CClassDefNode', 'PyClassDefNode', 'PropertyNode'):
            scope += (getattr(node, 'class_name', None) or node.name,)
        elif kind in ('DefNode', 'CFuncDefNode', 'LambdaNode'):
            if kind == 'LambdaNode':
                counts[scope + ('<lambda>',)] += 1
                name = '<lambda{}>'.format(counts[scope + ('<lambda>',)])
            elif kind == 'CFuncDefNode':
                declarator = node.declarator
                while not hasattr(declarator, 'name'):
                    declarator = declarator.base
                name = declarator.name
            else:
                name = node.name
            qualified = scope + (name,)
            key = '{}:{}'.format(relative, '.'.join(qualified))
            counts[key] += 1
            if counts[key] > 1:
                key += '[{}]'.format(counts[key])
            line = node.pos[1]
            indentation = len(lines[line - 1]) - len(lines[line - 1].lstrip())
            end = line
            # A whole function block, including nested definitions, is reviewed.
            while end < len(lines):
                text = lines[end]
                if text.strip() and len(text) - len(text.lstrip()) <= indentation:
                    break
                end += 1
            digest = hashlib.sha256((''.join(lines[line - 1:end]) + dependencies).encode()).hexdigest()
            sites[key] = {'file': relative, 'line': line, 'source_sha256': digest}
            scope = qualified
        for attr in getattr(node, 'child_attrs', ()):
            if attr == 'py_func_stat':
                # cpdef's synthetic Python wrapper is the same source site.
                continue
            child = getattr(node, attr, None)
            for item in child if isinstance(child, list) else [child]:
                if hasattr(item, 'child_attrs'):
                    visit(item, scope)

    visit(tree)
    return sites


def check_classifications(sites, classifications):
    """Fail on new, missing, stale or unclassified sites; no default verdict."""
    errors = []
    for site in sorted(sites.keys() - classifications.keys()):
        errors.append('unclassified new site: {}'.format(site))
    for site in sorted(classifications.keys() - sites.keys()):
        errors.append('removed site still classified: {}'.format(site))
    for site in sorted(sites.keys() & classifications.keys()):
        row = classifications[site]
        if row.get('verdict') not in VERDICTS:
            errors.append('missing or invalid verdict: {}'.format(site))
        if not isinstance(row.get('reason'), str) or not row['reason'].strip():
            errors.append('missing reason: {}'.format(site))
        if not isinstance(row.get('evidence'), list) or not row['evidence']:
            errors.append('missing evidence: {}'.format(site))
        if row.get('source_sha256') != sites[site]['source_sha256']:
            errors.append('changed source needs classification review: {}'.format(site))
    if errors:
        raise ValueError('\n'.join(errors))


def render_report(sites, classifications):
    """Generate a source-located table from the AST and checked verdict records."""
    check_classifications(sites, classifications)
    counts = Counter(row['verdict'] for row in classifications.values())
    lines = ['# Generated resolved-state database census', '',
             '{} function sites; Python AST and Cython parser traversal of every source under `rmgpy/` and `arkane/`, '
             'including methods, nested callbacks and lambdas. This is a conservative superset of data sites.'.format(len(sites)),
             '', 'Verdicts: {}.'.format(', '.join('{} {}'.format(counts[v], v) for v in sorted(counts))), '',
             'Counts by directory: {}.'.format(dict(sorted(Counter(str(Path(row['file']).parent) for row in sites.values()).items()))), '',
             '| Site | Source | Verdict | Reason | Evidence |', '| --- | --- | --- | --- | --- |']
    for key, site in sorted(sites.items()):
        row = classifications[key]
        lines.append('| `{}` | `{}:{}` | {} | {} | {} |'.format(
            key.split(':', 1)[1], site['file'], site['line'], row['verdict'],
            row['reason'].replace('|', '\\|'), ', '.join('`{}`'.format(e) for e in row['evidence'])))
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--manifest', type=Path)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    manifest = args.manifest or args.repository / 'test/rmgpy/data/stateDataCensus.json'
    classifications = json.loads(manifest.read_text())['sites']
    sites = collect_sites(args.repository)
    check_classifications(sites, classifications)
    if args.report:
        args.report.write_text(render_report(sites, classifications))
    print('Classified {} AST function sites: {}'.format(
        len(sites), dict(sorted(Counter(row['verdict'] for row in classifications.values()).items()))))


if __name__ == '__main__':
    main()
