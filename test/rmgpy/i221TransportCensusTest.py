#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
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

"""
I-221 round 112, HIGH 3 -- every process transport in the source tree carries the complete
reducers, checked mechanically rather than by a list someone keeps.

The reducers in `rmgpy/molecule/` lose state (`Atom.id`, `Atom.props`; a surface
`Molecule` or a `Fragment` raises), and fixing them there is out of this campaign's gates.
`family.install_complete_reducers()` hands `multiprocessing`'s pickler the completing
table instead. That is per-transport, so the question a reviewer has to be able to answer
is not "was it installed where the finding was" but "is there a transport it was not
installed at" -- rounds 95, 102, 107, 110 and 111 each closed the named site and left the
one beside it. This module answers it by parsing every source file:

* every module that starts a `multiprocessing` / `concurrent.futures` transport must call
  ``install_complete_reducers()`` at module level (so the parent has registered before it
  pickles the arguments, and a forked child inherits the registration);
* no module writes a pickle to storage -- counted, and zero is pinned as zero -- because a
  pickle written to a file is a transport with no dispatch table at all.

The census is also pinned by name, so a new transport turns this red and must be looked at
even if it does install the reducers.
"""

import ast
import os
import re
import subprocess
import sys
import textwrap

import pytest

import rmgpy

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(rmgpy.__file__)))
TREES = ('rmgpy', 'arkane', 'scripts')

#: Constructors that start another process or carry objects to one. `ForkingPickler`
#: itself is excluded: it is the pickler being completed, not a transport.
TRANSPORTS = {'Pool', 'Process', 'Pipe', 'Queue', 'SimpleQueue', 'JoinableQueue', 'Manager',
              'ProcessPoolExecutor', 'get_context'}
TRANSPORT_MODULES = ('multiprocessing', 'concurrent.futures')

#: The census at round 112 -- (module, spelling at the call). Extended from round 111's
#: "3 Pool sites", not re-derived: the same three Pools, plus the Pipe/Process pair in
#: family.py's `_spawn_tree_process` that a Pool-only count could not see. Both of those
#: pickle through ForkingPickler (Connection.send; Process args under spawn), so family.py's
#: module-level install covers them.
EXPECTED_TRANSPORTS = {
    ('rmgpy/data/kinetics/family.py', 'mp.Pool'),
    ('rmgpy/data/kinetics/family.py', 'mp.Pipe'),
    ('rmgpy/data/kinetics/family.py', 'mp.Process'),
    ('rmgpy/qm/main.py', 'Pool'),
    ('rmgpy/rmg/react.py', 'Pool'),
}

#: The only module that imports a pickler at all. This is what makes the plain-pickle path
#: -- `Reaction.__reduce__`, `LibraryReaction.__reduce__`, `ReactionModel.__reduce__`, all
#: of which nest the lossy reducers in `rmgpy/molecule/` -- accepted rather than a defect:
#: no code in the tree reaches it. family.py's own uses are an in-memory probe and
#: `_CompletePickler` on a BytesIO. A caller outside the tree who pickles a reaction with
#: plain `pickle` gets the lossy result; making that correct needs the molecule reducers,
#: which are outside this campaign's gates.
EXPECTED_PICKLE_IMPORTERS = {'rmgpy/data/kinetics/family.py'}

#: Pickle writes to storage: ``pickle.dump`` / ``Pickler(file)`` onto something that is
#: not an in-memory buffer. Round 111 counted zero; it is still zero.
EXPECTED_STORAGE = set()


def _sources():
    for tree in TREES:
        for dirpath, _, filenames in os.walk(os.path.join(ROOT, tree)):
            for name in sorted(filenames):
                if name.endswith(('.py', '.pyx')):
                    path = os.path.join(dirpath, name)
                    yield os.path.relpath(path, ROOT), path
    for name in sorted(os.listdir(ROOT)):
        if name.endswith('.py'):
            yield name, os.path.join(ROOT, name)


def _dotted(node):
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        inner = _dotted(node.value)
        return inner + '.' + node.attr if inner else None
    return None


def transports_in(source):
    """
    The transport calls in one module, as the spelling used at the call, and whether the
    module calls ``install_complete_reducers()`` at module level.

    Names are resolved through the module's own imports, so ``from multiprocessing import
    Pool``, ``import multiprocessing as mp`` and ``from concurrent import futures`` are all
    seen; a local class that happens to be called ``Pool`` is not.
    """
    tree = ast.parse(source)
    bound = {}  # local dotted prefix -> the transport module it names
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.startswith(TRANSPORT_MODULES):
                    bound[alias.asname or alias.name] = alias.name
        elif isinstance(node, ast.ImportFrom) and node.module:
            for alias in node.names:
                full = node.module + '.' + alias.name
                if node.module.startswith(TRANSPORT_MODULES) or full.startswith(TRANSPORT_MODULES):
                    bound[alias.asname or alias.name] = full
    calls = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        spelled = _dotted(node.func)
        if spelled is None:
            continue
        head, _, rest = spelled.partition('.')
        if head not in bound:
            continue
        resolved = bound[head] + ('.' + rest if rest else '')
        if resolved.rsplit('.', 1)[-1] in TRANSPORTS:
            calls.append(spelled)
    installs = any(isinstance(stmt, ast.Expr) and isinstance(stmt.value, ast.Call)
                   and _dotted(stmt.value.func) in ('install_complete_reducers',
                                                     'family.install_complete_reducers')
                   for stmt in tree.body)
    return calls, installs


def storage_in(source):
    """``pickle.dump(...)`` calls and pickler constructions given something other than a
    ``BytesIO`` -- a pickle that leaves the process through a file."""
    tree = ast.parse(source)
    found = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        spelled = _dotted(node.func) or ''
        if spelled.endswith(('pickle.dump', 'dill.dump', 'joblib.dump', 'shelve.open')):
            found.append(spelled)
        elif spelled.endswith('Pickler') and node.args:
            target = node.args[0]
            if not (isinstance(target, ast.Name) and 'buffer' in target.id.lower()) and \
                    not (isinstance(target, ast.Call) and (_dotted(target.func) or '').endswith('BytesIO')):
                found.append(spelled)
    return found


def pickle_importers():
    found = set()
    for rel, path in _sources():
        if rel.endswith('.pyx'):
            continue
        with open(path, encoding='utf-8') as f:
            tree = ast.parse(f.read())
        for node in ast.walk(tree):
            names = [a.name for a in node.names] if isinstance(node, ast.Import) else \
                [node.module] if isinstance(node, ast.ImportFrom) and node.module else []
            if any(n.split('.')[0] in ('pickle', 'cPickle', '_pickle', 'dill', 'joblib', 'shelve')
                   for n in names):
                found.add(rel)
    return found


def _census():
    transports, missing, storage, pyx = set(), [], set(), []
    for rel, path in _sources():
        with open(path, encoding='utf-8') as f:
            source = f.read()
        if rel.endswith('.pyx'):
            # Cython is not Python syntax; a mention is enough to make it someone's problem.
            if re.search(r'\b(multiprocessing|concurrent\.futures|pickle)\b', source):
                pyx.append(rel)
            continue
        calls, installs = transports_in(source)
        for spelled in calls:
            transports.add((rel, spelled))
        if calls and not installs:
            missing.append(rel)
        storage.update((rel, spelled) for spelled in storage_in(source))
    return transports, sorted(set(missing)), storage, pyx


class TestTransportCensus:

    def test_every_transport_installs_the_complete_reducers(self):
        _, missing, _, _ = _census()
        assert missing == [], (
            'these modules start a process transport without calling '
            'install_complete_reducers() at module level, so RMG objects crossing it go '
            'through the lossy reducers in rmgpy/molecule/: {0}'.format(missing))

    def test_the_census_is_the_one_on_record(self):
        transports, _, _, _ = _census()
        assert transports == EXPECTED_TRANSPORTS

    def test_no_pickle_is_written_to_storage(self):
        _, _, storage, _ = _census()
        assert storage == EXPECTED_STORAGE

    def test_plain_pickle_is_reached_from_one_module(self):
        assert pickle_importers() == EXPECTED_PICKLE_IMPORTERS

    def test_no_cython_module_touches_a_transport(self):
        _, _, _, pyx = _census()
        assert pyx == []


class TestTheCheckCanFail:
    """The census is only worth its green if each arm turns red on the shape it refuses."""

    @pytest.mark.parametrize('source', [
        'from multiprocessing import Pool\ndef f():\n    Pool(2).map(g, [])\n',
        'import multiprocessing as mp\ndef f():\n    mp.Pool(2)\n',
        'import multiprocessing\ndef f():\n    multiprocessing.get_context("spawn").Pool(2)\n',
        'from concurrent.futures import ProcessPoolExecutor as E\ndef f():\n    E(2)\n',
        # Installed inside a function is not installed before the parent pickles.
        'from multiprocessing import Pool\n'
        'def f():\n    install_complete_reducers()\n    Pool(2)\n',
    ])
    def test_uninstalled_transport_is_seen(self, source):
        calls, installs = transports_in(source)
        assert calls and not installs

    def test_installed_transport_passes(self):
        calls, installs = transports_in(
            'from multiprocessing import Pool\n'
            'from rmgpy.data.kinetics.family import install_complete_reducers\n'
            'install_complete_reducers()\n'
            'def f():\n    Pool(2)\n')
        assert calls == ['Pool'] and installs

    def test_a_local_pool_is_not_a_transport(self):
        assert transports_in('class Pool:\n    pass\nPool()\n') == ([], False)

    @pytest.mark.parametrize('source, found', [
        ('import pickle\npickle.dump(x, open("f", "wb"))\n', ['pickle.dump']),
        ('import pickle\npickle.Pickler(open("f", "wb"))\n', ['pickle.Pickler']),
        ('import pickle, io\npickle.Pickler(io.BytesIO())\n', []),
        ('import pickle\npickle.loads(pickle.dumps(x))\n', []),
    ])
    def test_storage_is_seen(self, source, found):
        assert storage_in(source) == found


class TestFreshImportRegisters:
    """Round 111's measurement, repeated in a fresh interpreter for each transport module:
    importing it alone must leave `Molecule` and `Species` in ForkingPickler's table."""

    @pytest.mark.parametrize('module', sorted({rel for rel, _ in EXPECTED_TRANSPORTS}))
    def test_fresh_import(self, module):
        name = module[:-3].replace('/', '.')
        probe = textwrap.dedent('''
            import importlib
            from multiprocessing.reduction import ForkingPickler
            importlib.import_module({0!r})
            from rmgpy.molecule.molecule import Molecule, Atom
            from rmgpy.species import Species
            table = ForkingPickler._extra_reducers
            print(all(cls in table for cls in (Atom, Molecule, Species)))
        ''').format(name)
        out = subprocess.run([sys.executable, '-c', probe], cwd=ROOT, capture_output=True,
                             text=True, check=True)
        assert out.stdout.strip().splitlines()[-1] == 'True', out.stderr
