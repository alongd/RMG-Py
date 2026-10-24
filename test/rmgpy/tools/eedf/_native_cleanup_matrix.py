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

"""Public-scan cleanup statement coverage under real nested raising signals."""
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
ROOT_PID = os.getpid()
MODE = sys.argv[1]
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else None
state = {}
held = []


class AppSignal(Exception):
    pass


class EntryFailure(Exception):
    pass


def handler(signum, frame):
    state['entered'] += 1
    if state['entered'] < state['count']:
        try:
            os.kill(ROOT_PID, signal.SIGUSR1)
        except AppSignal:
            pass
    state['raised'] += 1
    raise AppSignal('real raise %d' % state['raised'])


signal.signal(signal.SIGUSR1, handler)


def hook():
    if n.is_worker():
        if state['fixture'] == 'exec':
            os.execv('/usr/bin/sleep', ['sleep', '30'])
        while True:
            time.sleep(30)


os.register_at_fork(after_in_child=hook)
original_fork = n.fork
original_process = c.Owner.process
original_cleanup = c.Owner.cleanup


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
            while True:
                if os.readlink('/proc/%d/exe' % pid) == '/usr/bin/sleep':
                    break
                assert time.monotonic() < until, 'exec fixture not ready'
                time.sleep(.001)
        raise EntryFailure('start exception enters cleanup')
    return pid


def cleanup(owner, seconds=5.):
    state['cleanups'] += 1
    if state['need_retry'] and state['cleanups'] == 1:
        raise EntryFailure('ordinary failure reaches retry bookkeeping')
    return original_cleanup(owner, seconds)


c.Owner.process = process
n.fork = fork
c.Owner.cleanup = cleanup
b.WORKER_KILL_TIMEOUT = .15


class Log(logging.Handler):
    def emit(self, record):
        state['messages'].append(record.getMessage())


logger = logging.getLogger(b.__name__)
logger.addHandler(Log())
logger.propagate = False


class Adapter:
    def integrate(self, seed, power=None):
        raise AssertionError('held child must not integrate')


def locations():
    function = getattr(b, '_integrate_owned', b._integrate_many)
    source, first = inspect.getsourcelines(function)
    module = ast.parse(''.join(source))
    blocks = [node for node in ast.walk(module) if isinstance(node, ast.Try)
              and any(isinstance(part, ast.Assign) and any(isinstance(t, ast.Name)
                      and t.id == 'interrupted' for t in part.targets) for part in node.finalbody)]
    assert len(blocks) == 1
    starts = {line for _, line in dis.findlinestarts(function.__code__)}
    rows = {}
    for part in blocks[0].finalbody:
        for node in ast.walk(part):
            if isinstance(node, (ast.stmt, ast.ExceptHandler)):
                line = first + node.lineno - 1
                if line in starts:
                    rows[line] = (function, line, source[node.lineno - 1].strip())
    return list(rows.values())


def trace(frame, event, arg):
    if (os.getpid() == ROOT_PID and event == 'line' and not state['injected']
            and frame.f_code is state['code'] and frame.f_lineno == state['line']):
        state['injected'] = True
        state['site_seen'] = (frame.f_code.co_name, frame.f_lineno)
        os.kill(ROOT_PID, signal.SIGUSR1)
    return trace


def trial(function, line, statement, count, directory):
    state.clear()
    # Retry setup reaches both except clauses and the final re-raise. A live
    # exec replacement reaches every pending-diagnostic and typed-error line.
    pending_words = ('facts =', 'pending =', 'message =', 'logging.getLogger',
                     'raise NativeWorkerStatusError', 'raise _native_cleanup_error',
                     'except (NativeWorkerStatusError')
    fixture = 'exec' if statement.startswith(pending_words) or statement == 'raise' else 'held'
    state.update(code=function.__code__, line=line, count=count, entered=0, raised=0,
                 fixture=fixture, messages=[], cleanups=0, need_retry=True,
                 injected=False, record=None)
    path = directory / 'branches.json'
    path.write_bytes(b'previous artifact\n')
    before = (path.read_bytes(), path.stat().st_mtime_ns)
    error = None
    started = time.monotonic()
    sys.settrace(trace)
    try:
        b.scan(Adapter(), (0., 0.), path, u_steps=1, electron_density_decades=1,
               processes=1, timeout=2., lease_s=.5, reference_power=.5)
    except BaseException as exc:
        error = exc
        held.append(exc)  # Retain traceback and Owner through every verdict.
    finally:
        sys.settrace(None)
    facts = n.info(state['record'])
    exact_report = isinstance(error, c.NativeWorkerStatusError) and any(
        'pending child PIDs [%d]' % facts['pid'] in message for message in state['messages'])
    passed = (facts['reaped'] == 1 or exact_report) and before == (path.read_bytes(), path.stat().st_mtime_ns)
    row = dict(statement=statement, line=line, function=function.__name__, count=count,
               entered=state['entered'], raised=state['raised'], injected=state['injected'],
               fixture=fixture, record=facts, messages=state['messages'], error=repr(error),
               elapsed_s=time.monotonic()-started, passed=passed)
    print('CLEANUP_SITE', json.dumps(row), flush=True)
    if OUT:
        with (OUT/('matrix-'+MODE+'.jsonl')).open('a') as out:
            out.write(json.dumps(row)+'\n')
    try:
        assert state['injected'] and state['raised'] == state['entered'] == count, row
        if MODE != 'red':
            assert passed, 'public cleanup escaped without exact reap/report: %r' % row
        return row
    finally:
        # Rescue only after recording/assessing the public verdict. Never call
        # Owner.cleanup; this does not contribute to the pass predicate.
        if facts['pid'] and not facts['reaped']:
            os.kill(facts['pid'], signal.SIGKILL)
            until = time.monotonic()+2
            while not n.poll(state['record']):
                assert time.monotonic()<until
                time.sleep(.001)


def main():
    rows = []
    with tempfile.TemporaryDirectory() as directory:
        for function, line, statement in locations():
            for count in (1, 2, 3, 5):
                rows.append(trial(function, line, statement, count, Path(directory)))
    if MODE == 'red':
        assert any(not row['passed'] for row in rows), 'baseline did not reproduce any lost cleanup'
    print('CLEANUP_MATRIX', json.dumps(dict(mode=MODE, rows=len(rows),
          passed=sum(row['passed'] for row in rows), failed=sum(not row['passed'] for row in rows))), flush=True)


if __name__ == '__main__':
    main()
