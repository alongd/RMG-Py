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

"""Terminal acceptance effects through the compiled reactor and public provider."""

import ast
import json
import logging
import os
from pathlib import Path
import signal
import sys
from types import SimpleNamespace

import pytest

from eedfProviderTest import public_qualification_fixture, intended_signal_disposition
from plasmaEEDFCoreTest import build_reactor
from rmgpy.exceptions import MechanismWriterError
from rmgpy.rmg.main import RMG
from rmgpy.solver.eedf import EEDFError
from rmgpy.solver.eedf_provider import EEDFProvider


def completion_effect(rmg):
    """Run the actual execute completion block without a model-generation campaign."""
    source = Path(sys.modules[RMG.__module__].__file__)
    tree = ast.parse(source.read_text())
    execute = next(node for node in ast.walk(tree)
                   if isinstance(node, ast.FunctionDef) and node.name == 'execute')
    block = next(node for node in execute.body if isinstance(node, ast.If)
                 and any(isinstance(value, ast.Constant)
                         and value.value == 'MODEL GENERATION COMPLETED'
                         for value in ast.walk(node)))
    code = compile(ast.fix_missing_locations(ast.Module(body=[block], type_ignores=[])),
                   str(source), 'exec')
    exec(code, {'self': rmg, 'logging': logging, 'end_early': False})


@pytest.mark.parametrize('lifecycle', ['fresh', 'reused'])
@pytest.mark.parametrize('signum', [signal.SIGINT, signal.SIGTERM])
@pytest.mark.parametrize('boundary', ['gate-return', 'wrapper-return'])
@pytest.mark.parametrize('route', ['published', 'qualify'])
def test_revoked_terminal_has_no_acceptance_effect(
        tmp_path, monkeypatch, caplog, lifecycle, signum, boundary, route):
    reactor, _, _ = build_reactor(tmp_path / 'reactor')
    provider, state = public_qualification_fixture(tmp_path / 'qualification', monkeypatch)
    if route == 'published' or lifecycle == 'reused':
        provider.qualify(state, 2.5)
    if route == 'published' and lifecycle == 'reused':
        provider.qualify(state, 3.5)
    reactor.eedf_provider = provider
    provider.bind_qualification_state(lambda y, t: state)
    reactor.energy_budget = {'A6a_relative': 0., 'A6a_steady_relative': 0.,
                             'A6b_relative': 0., 'A6b_tolerance': 1.e-6}
    events, failures = [], []

    def handler(received, frame):
        events.append(('received', received))
        failures.append(provider.record_terminal_failure(state, 4., RuntimeError('shutdown')))

    target = (EEDFProvider.require_qualified.__wrapped__.__code__ if boundary == 'gate-return'
              else EEDFProvider.require_qualified.__code__)

    def trace(frame, event, arg):
        if (event == 'return' and frame.f_code is target and not events
                and (boundary == 'gate-return'
                     or frame.f_locals.get('method') is EEDFProvider.require_qualified.__wrapped__)):
            events.append(('sent', signum))
            os.kill(os.getpid(), signum)
        return trace

    with intended_signal_disposition(signum, handler):
        try:
            sys.settrace(trace)
            reactor.validate_terminal_state()
        finally:
            sys.settrace(None)
    assert events == [('sent', signum), ('received', signum)]
    context = SimpleNamespace(reaction_systems=[reactor],
                              reaction_model=SimpleNamespace(get_model_size=lambda: (3, 1, 0, 0)),
                              export_failures=[('example output', OSError('not written'))])
    with caplog.at_level(logging.INFO):
        completion_effect(context)
        with pytest.raises(MechanismWriterError):
            RMG.report_export_failures(context)
    destination = tmp_path / 'terminal-qualified.json'
    admission = getattr(reactor, 'terminal_admission', None)
    refused = False
    try:
        if admission is None:
            # Baseline protocol: normal reactor return authorized the consumer.
            destination.write_bytes(b'{"status":"PASS"}')
        else:
            admission.commit_bytes(destination, b'{"status":"PASS"}')
    except EEDFError:
        refused = True
    assert 'MODEL GENERATION COMPLETED' not in caplog.text
    assert 'Model generation itself completed' not in caplog.text
    assert 'valid; the output above' not in caplog.text
    assert refused and not destination.exists()
    disk = json.loads(Path(state['record_path']).read_text())
    assert disk['status'] == 'FAIL' and disk['run_id'] == failures[0]['run_id']


def test_non_eedf_completion_wording_is_unchanged(caplog):
    context = SimpleNamespace(reaction_systems=[SimpleNamespace(eedf_mode=False)],
                              reaction_model=SimpleNamespace(get_model_size=lambda: (3, 1, 0, 0)),
                              export_failures=[('example output', OSError('not written'))])
    with caplog.at_level(logging.INFO):
        completion_effect(context)
        with pytest.raises(MechanismWriterError):
            RMG.report_export_failures(context)
    assert 'MODEL GENERATION COMPLETED' in caplog.text
    assert 'Model generation itself completed and every file written before this point is' in caplog.text
    assert 'valid; the output above names what is missing.' in caplog.text


@pytest.mark.parametrize('kind', ['export', 'acceptance'])
@pytest.mark.parametrize('policy', ['returning', 'raising'])
@pytest.mark.parametrize('boundary', [
    'truncate-before', 'truncate-after', 'marker-after', 'capture-before',
    'capture-after', 'private-unlink', 'restore-before', 'restore-after', 'rmdir',
])
def test_second_signal_preserves_terminal_cleanup(
        tmp_path, monkeypatch, record_property, kind, policy, boundary):
    import rmgpy.solver.eedf_provider as module
    provider, state = public_qualification_fixture(tmp_path / 'fixture', monkeypatch)
    provider.qualify(state, 2.5)
    admission = provider.require_qualified('terminal output')
    destination = (tmp_path / 'output.json' if kind == 'export'
                   else Path(state['record_path']))
    foreign = tmp_path / 'foreign'
    foreign.write_bytes(b'other writer ordinary output')
    foreign_identity = foreign.stat().st_ino
    private_before = set(tmp_path.parent.rglob('.eedf-*.private'))
    native = {name: getattr(os, name) for name in
              ('link', 'rename', 'unlink', 'rmdir', 'ftruncate', 'pwrite', 'replace')}
    original_unlink = Path.unlink
    original_validate = provider._validate_acceptance_manifest
    events, traces = [], []
    active, second_sent, first_sent = [False], [False], [False]
    owned, expected = [None], [None]
    legacy = not hasattr(module, '_QualificationCapture')

    def send(site):
        if active[0] and not second_sent[0] and (site == boundary or site == 'legacy-unlink'):
            second_sent[0] = True
            if foreign.exists():
                native['replace'](foreign, destination)
            expected[0] = provider._qualification_failure['record']['run_id']
            events.append(('second-sent', site))
            os.kill(os.getpid(), signal.SIGINT)

    def second(signum, frame):
        import traceback
        events.append(('second-received', frame.f_code.co_name))
        traces.extend(traceback.format_stack(frame))
        if policy == 'raising':
            raise KeyboardInterrupt('second signal inside cleanup')
        try:
            provider.record_terminal_failure(state, 4., RuntimeError('second shutdown'))
        except EEDFError:
            # This really is a returning handler, including ownership refusals.
            events.append(('second-refusal-caught',))
        events.append(('second-returned',))

    def first(signum, frame):
        events.append(('first-received', signum))
        active[0] = True
        try:
            provider.require_qualified('shutdown output')
        except BaseException as error:
            events.append(('first-caught', type(error).__name__))
        finally:
            active[0] = False

    def validate(*args, **kwargs):
        if active[0]:
            raise OSError('first handler manifest read failed')
        return original_validate(*args, **kwargs)

    def start_first():
        if not first_sent[0]:
            first_sent[0] = True
            if boundary.startswith('restore-') and not legacy:
                native['replace'](foreign, destination)
            events.append(('first-sent', signal.SIGTERM))
            os.kill(os.getpid(), signal.SIGTERM)

    def is_owned(fd):
        info = os.fstat(fd)
        return (info.st_dev, info.st_ino) == owned[0]

    def truncate(fd, length):
        if is_owned(fd):
            send('truncate-before')
        result = native['ftruncate'](fd, length)
        if is_owned(fd):
            send('truncate-after')
        return result

    def pwrite(fd, payload, offset):
        result = native['pwrite'](fd, payload, offset)
        if is_owned(fd):
            send('marker-after')
        return result

    def capture(source, target, *args, **kwargs):
        match = os.fspath(source) == destination.name
        if match:
            send('capture-before')
        result = native['rename'](source, target, *args, **kwargs)
        if match:
            send('capture-after')
        return result

    def private_unlink(path, *args, **kwargs):
        if os.fspath(path) == 'captured':
            send('private-unlink')
        return native['unlink'](path, *args, **kwargs)

    def legacy_unlink(path, *args, **kwargs):
        if legacy and active[0] and path == destination and not second_sent[0]:
            # Baseline owns no private capture. Replace after its check, then
            # deliver the second signal at that actual destructive boundary.
            if foreign.exists():
                native['replace'](foreign, destination)
            result = original_unlink(path, *args, **kwargs)
            send('legacy-unlink')
            return result
        return original_unlink(path, *args, **kwargs)

    def rmdir(path, *args, **kwargs):
        if str(path).startswith('.eedf-'):
            send('rmdir')
        return native['rmdir'](path, *args, **kwargs)

    def link(source, target, *args, **kwargs):
        restoring = os.fspath(source) == 'captured'
        if restoring:
            send('restore-before')
        result = native['link'](source, target, *args, **kwargs)
        if restoring:
            send('restore-after')
        if Path(target) == destination and not first_sent[0]:
            info = destination.stat()
            owned[0] = (info.st_dev, info.st_ino)
            if kind == 'export':
                start_first()
        return result

    def trace(frame, event, arg):
        if kind == 'acceptance' and event == 'return' and frame.f_code is module.EEDFProvider._publish_acceptance.__code__:
            if owned[0] is None and destination.exists():
                info = destination.stat()
                owned[0] = (info.st_dev, info.st_ino)
            start_first()
        return trace

    monkeypatch.setattr(provider, '_validate_acceptance_manifest', validate)
    monkeypatch.setattr(Path, 'unlink', legacy_unlink)
    for name, function in [('link', link), ('rename', capture), ('unlink', private_unlink),
                           ('rmdir', rmdir), ('ftruncate', truncate), ('pwrite', pwrite)]:
        monkeypatch.setattr(os, name, function)
    outcome = None
    with intended_signal_disposition(signal.SIGTERM, first), intended_signal_disposition(signal.SIGINT, second):
        try:
            sys.settrace(trace)
            if kind == 'export':
                admission.commit_bytes(destination, b'{"status":"PASS"}')
                outcome = 'PASS'
            else:
                outcome = provider.qualify(state, 3.5)['status']
        except EEDFError:
            outcome = 'refused'
        finally:
            sys.settrace(None)
    record_property('events', repr(events))
    record_property('second_trace', ''.join(traces))
    assert first_sent[0] and second_sent[0], events
    assert sum(e[0] == 'second-received' for e in events) == 1
    assert next(e[1] for e in events if e[0] == 'second-received') == 'send'
    assert (('second-returned',) in events) == (policy == 'returning')
    assert outcome != 'PASS'
    assert provider._qualification_status == 'FAILED'
    assert provider._qualification_record['run_id'] == expected[0]
    private_after = set(tmp_path.parent.rglob('.eedf-*.private')) - private_before
    candidates = list(tmp_path.rglob('*'))
    for directory in private_after:
        assert any(directory.iterdir()), 'empty private directory leaked'
        candidates.extend(directory.rglob('*'))
    preserved = [path for path in candidates if path.is_file()
                 and path.stat().st_ino == foreign_identity
                 and path.read_bytes() == b'other writer ordinary output']
    assert preserved, 'foreign file lost during second-signal cleanup'
    for path in [destination] + [file for directory in private_after for file in directory.rglob('*')]:
        if path.is_file():
            try:
                assert json.loads(path.read_bytes()).get('status') != 'PASS'
            except (ValueError, UnicodeError):
                pass
