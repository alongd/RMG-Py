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

"""Public scan guard boundaries and repeated interruptions inside reporting."""
import ast
import dis
import inspect
import json
import logging
import os
from pathlib import Path
import signal
import sys
import tempfile
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from rmgpy.tools.eedf import branches as b
from rmgpy.tools.eedf import _native_context as c
n = c.native
MODE = sys.argv[1]
ROOT = os.getpid()
state = {}
held = []


class AppSignal(Exception):
    pass


class EntryFailure(Exception):
    pass


def signal_handler(*unused):
    state['raises'] += 1
    raise AppSignal('guard signal %d' % state['raises'])


signal.signal(signal.SIGUSR1, signal_handler)


def child_hook():
    if n.is_worker():
        if state['fixture'] == 'exec':
            os.execv('/usr/bin/sleep', ['sleep', '30'])
        if state['fixture'] == 'external':
            os._exit(0)
        time.sleep(30)


os.register_at_fork(after_in_child=child_hook)
original_process = c.Owner.process
original_fork = n.fork
original_cleanup = c.Owner.cleanup
original_report = b._native_cleanup_error
b.WORKER_KILL_TIMEOUT = .15


def process(owner, target, args):
    child = original_process(owner, target, args)
    state['record'] = child._native_record
    held.append(owner)
    return child


def fork(record):
    pid = original_fork(record)
    if pid:
        if state['fixture'] == 'exec':
            until = time.monotonic() + 2
            while os.readlink('/proc/%d/exe' % pid) != '/usr/bin/sleep':
                assert time.monotonic() < until
                time.sleep(.001)
        if state['fixture'] == 'external':
            assert os.waitpid(pid, 0) == (pid, 0)
        if MODE == 'early-fork':
            os.kill(ROOT, signal.SIGUSR1)
        raise EntryFailure('fork return before Popen assignment')
    return pid


def cleanup(owner, seconds=5.):
    if MODE == 'report':
        raise EntryFailure('force escape through retry bookkeeping')
    return original_cleanup(owner, seconds)


def trace(frame, event, arg):
    if os.getpid() != ROOT:
        return trace
    if MODE == 'report':
        if (frame.f_code is b._integrate_owned.__code__ and event == 'line'
                and lines[frame.f_lineno - 1].strip() == 'interrupted = exc'):
            state['escaped_cleanup'] = True
            raise EntryFailure('escaped Python retry')
        if (frame.f_code is original_report.__code__ and event == 'line'
                and frame.f_lineno == state['line'] and state['raises'] < state['count']):
            state['sites'].append((frame.f_code.co_name, frame.f_lineno))
            os.kill(ROOT, signal.SIGUSR1)
    elif frame.f_code is b._integrate_owned.__code__ and not state['sites']:
        selected = {'entry': 'call', 'exception': 'exception', 'return': 'return'}.get(MODE)
        if event == selected or (MODE == 'after-cleanup' and event == 'line'
                and lines[frame.f_lineno - 1].strip() == 'break'):
            state['sites'].append((event, frame.f_lineno))
            os.kill(ROOT, signal.SIGUSR1)
    return trace


def report(owner, wait_error=0):
    state['attempts'] += 1
    # Raising out of a Python trace disables tracing. Reinstall it at every
    # C-initiated reporting attempt, so each requested raise tests a fresh retry.
    sys.settrace(trace)
    return original_report(owner, wait_error)


class Log(logging.Handler):
    def emit(self, record):
        state['messages'].append(record.getMessage())


logger = logging.getLogger(b.__name__)
logger.addHandler(Log())
logger.propagate = False
c.Owner.process = process
c.Owner.cleanup = cleanup
n.fork = fork
if MODE == 'report':
    b._native_cleanup_error = report
lines = Path(b.__file__).read_text().splitlines()


class ReturnedReference:
    def __del__(self):
        state['release_facts'] = n.info(state['record'])


def boundary_body(owner, deferral, adapter, seeds, processes, timeout, power):
    child = owner.process(lambda: None, ())
    child.start()
    if MODE == 'copied-guard':
        app_pid = os.fork()
        if app_pid == 0:
            return []  # Return through the copied C guard, without its authority.
        state['app_pid'] = app_pid
        raise EntryFailure('parent returns directly to native guard')
    return ReturnedReference()


if MODE in ('copied-guard', 'reference-release'):
    n.fork = original_fork
    b._integrate_owned = boundary_body
if MODE == 'reference-release':
    original_guard = n.guard

    def release_result(*args):
        result = original_guard(*args)
        del result
        os.kill(ROOT, signal.SIGUSR1)

    n.guard = release_result


class Adapter:
    def integrate(self, seed, power=None):
        raise AssertionError('child hook must hold the worker')


def trial(directory, line=0, statement='', count=1, fixture='held'):
    state.clear()
    state.update(line=line, count=count, raises=0, fixture=fixture, record=None,
                 sites=[], messages=[], attempts=0, escaped_cleanup=False)
    path = directory/'branches.json'
    path.write_bytes(b'old artifact\n')
    before = path.read_bytes(), path.stat().st_mtime_ns
    error = None
    started = time.monotonic()
    sys.settrace(trace)
    try:
        b.scan(Adapter(), (0., 0.), path, u_steps=1, electron_density_decades=1,
               processes=1, timeout=2., lease_s=.5, reference_power=.5)
    except BaseException as exc:
        held.append(exc)
        error = exc
    finally:
        sys.settrace(None)
    if os.getpid() != ROOT:
        os._exit(0 if isinstance(error, RuntimeError) and 'another process' in str(error) else 91)
    if MODE == 'copied-guard':
        assert os.waitpid(state['app_pid'], 0) == (state['app_pid'], 0)
    facts = n.info(state['record']) if state['record'] is not None else None
    row = dict(mode=MODE, line=line, statement=statement, count=count,
               raises=state['raises'], sites=state['sites'], attempts=state['attempts'],
               messages=state['messages'], record=facts, error=repr(error),
               elapsed_s=time.monotonic()-started, escaped_cleanup=state['escaped_cleanup'])
    print('NATIVE_GUARD', json.dumps(row), flush=True)
    try:
        assert before == (path.read_bytes(), path.stat().st_mtime_ns)
        assert state['raises'] == (0 if MODE == 'copied-guard' else count), row
        if MODE == 'report':
            assert state['escaped_cleanup'] and state['attempts'] == count+1, row
            assert isinstance(error, c.NativeWorkerStatusError), row
            assert state['messages'] and str(facts['pid']) in state['messages'][-1], row
            assert facts['reaped'] == (2 if fixture == 'external' else 0), row
            assert row['elapsed_s'] < .75, row
        elif MODE == 'entry':
            assert facts is None and isinstance(error, AppSignal), row
        else:
            assert facts['reaped'] == 1 and not Path('/proc/%d' % facts['pid']).exists(), row
            if MODE == 'reference-release':
                assert state['release_facts']['reaped'] == 1, state
    finally:
        # Exact test exec rescue follows the recorded public verdict only.
        if facts and facts['pid'] and not facts['reaped']:
            os.kill(facts['pid'], signal.SIGKILL)
            until = time.monotonic()+2
            while not n.poll(state['record']):
                assert time.monotonic()<until
                time.sleep(.001)


def main():
    with tempfile.TemporaryDirectory() as directory:
        directory = Path(directory)
        if MODE != 'report':
            trial(directory)
            return
        source, first = inspect.getsourcelines(original_report)
        tree = ast.parse(''.join(source))
        starts = {line for _, line in dis.findlinestarts(original_report.__code__)}
        targets = {}
        for node in ast.walk(tree.body[0]):
            if isinstance(node, ast.stmt) and not isinstance(node, ast.FunctionDef):
                line = first + node.lineno - 1
                statement = source[node.lineno-1].strip()
                if line in starts and not statement.startswith('"""'):
                    targets[line] = statement
        for line, statement in sorted(targets.items()):
            if statement.startswith("message += '; waitpid errno"):
                # This branch requires an actual non-EINTR/non-ECHILD waitpid
                # failure; neither pending nor externally reaped fixtures fake it.
                continue
            fixture = 'external' if statement.startswith("message += '; externally") else 'exec'
            for count in (1, 2, 3, 5):
                trial(directory, line, statement, count, fixture)
    print('NATIVE_GUARD_PASS', MODE, flush=True)


if __name__ == '__main__':
    main()
