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

"""Standalone: reporter forks retain no authority over the parent ledger.

Run without a pytest supervisor: owner + one worker + one application child.
The child cache turns to ECHILD if even one foreign raw wait is attempted.
"""
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
MODE, OUTCOME = sys.argv[1].split('-')
ROOT = os.getpid()
state = dict(attempts=0, raises=0, fork_site=None, messages=[], signals=[])
held = []
read_fd, write_fd = os.pipe()


class AppSignal(Exception):
    pass


def signal_handler(*unused):
    state['raises'] += 1
    raise AppSignal('real reporting interruption')


def audit(event, args):
    if event == 'os.kill':
        state['signals'].append(list(args))


sys.addaudithook(audit)
signal.signal(signal.SIGUSR1, signal_handler)


def child_hook():
    if n.is_worker():
        os.execv('/usr/bin/sleep', ['sleep', '0.4' if OUTCOME == 'reap' else '30'])


os.register_at_fork(after_in_child=child_hook)


def app_fork(site):
    state['fork_site'] = site
    pid = os.fork()
    if pid:
        state['app_pid'] = pid
    else:
        os.close(read_fd)
        state['child_initial'] = n.info(state['record'])
        state['child_report_calls'] = 0
        state['signals'] = []


class ReporterError(Exception):
    def __del__(self):
        app_fork('exception-destruction')


class Log(logging.Handler):
    def emit(self, record):
        state['messages'].append(record.getMessage())


logger = logging.getLogger(b.__name__)
logger.addHandler(Log())
logger.propagate = False
original_report = b._native_cleanup_error


def child_observation(error):
    row = dict(error=repr(error), facts=n.info(state['record']),
               initial=state['child_initial'], report_calls=state['child_report_calls'],
               signals=state['signals'], messages=state['messages'])
    good = (isinstance(error, RuntimeError) and 'another process' in str(error)
            and row['facts']['reaped'] == row['initial']['reaped'] == 0
            and row['report_calls'] == 0 and not row['messages']
            and all(pid == os.getpid() for pid, sig in row['signals']))
    row['pass'] = good
    os.write(write_fd, (json.dumps(row)+'\n').encode())
    os._exit(0 if good else 91)


def report(owner, wait_error=0):
    if os.getpid() != ROOT:
        state['child_report_calls'] += 1
        # Bound the red case at the first illegal retry, with no test cleanup.
        child_observation(None)
    state['attempts'] += 1
    if state['attempts'] == 1:
        if MODE == 'callback':
            app_fork('callback')
            os.kill(os.getpid(), signal.SIGUSR1)
        raise ReporterError('fork when C clears this reporter exception')
    if state['attempts'] == 2:
        if OUTCOME == 'reap':
            # The exec fixture exits naturally. No wait/poll/kill assists C.
            time.sleep(.5)
        os.kill(ROOT, signal.SIGUSR1)
    return original_report(owner, wait_error)


def body(owner, deferral, adapter, seeds, processes, timeout, power):
    held.append(owner)
    child = owner.process(lambda: None, ())
    held.append(child)
    state['record'] = child._native_record
    child.start()
    pid = n.info(state['record'])['pid']
    until = time.monotonic()+2
    while os.readlink('/proc/%d/exe' % pid) != '/usr/bin/sleep':
        assert time.monotonic() < until
        time.sleep(.001)
    raise AppSignal('enter C guard without Python cleanup')


b._integrate_owned = body
b._native_cleanup_error = report
b.WORKER_KILL_TIMEOUT = .02


class Adapter:
    def integrate(self, seed, power=None):
        raise AssertionError('exec hook must replace worker')


def main():
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory)/'branches.json'
        path.write_bytes(b'old artifact\n')
        before = path.read_bytes(), path.stat().st_mtime_ns
        error = None
        try:
            b.scan(Adapter(), (0., 0.), path, u_steps=1, electron_density_decades=1,
                   processes=1, timeout=2., lease_s=.5, reference_power=.5)
        except BaseException as exc:
            error = exc
            held.append(exc)
        if os.getpid() != ROOT:
            child_observation(error)
        os.close(write_fd)
        child_row = json.loads(os.read(read_fd, 65536))
        os.close(read_fd)
        app_status = os.waitpid(state['app_pid'], 0)
        facts = n.info(state['record'])
        row = dict(mode=MODE, outcome=OUTCOME, child=child_row, app_status=app_status,
                   parent_facts=facts, error=repr(error), state=state.copy())
        row['state'].pop('record')
        print('REPORT_FORK', json.dumps(row), flush=True)
        try:
            assert before == (path.read_bytes(), path.stat().st_mtime_ns)
            assert state['fork_site'] == ('callback' if MODE == 'callback' else 'exception-destruction')
            assert state['raises'] == (2 if MODE == 'callback' else 1), row
            assert state['attempts'] == (2 if OUTCOME == 'reap' else 3), row
            if OUTCOME == 'reap':
                assert facts['reaped'] == 1 and isinstance(error, AppSignal), row
                assert not state['messages'] and not Path('/proc/%d' % facts['pid']).exists(), row
                # Raw identity is consumed exactly once; later reads retain it.
                assert n.info(state['record']) == facts
                try:
                    os.waitpid(facts['pid'], os.WNOHANG)
                except ChildProcessError:
                    pass
                else:
                    raise AssertionError('guard did not consume exact wait status')
            else:
                assert facts['reaped'] == 0 and isinstance(error, c.NativeWorkerStatusError), row
                assert len(state['messages']) == 1 and str(facts['pid']) in state['messages'][0], row
                assert Path('/proc/%d' % facts['pid']).exists(), row
            assert child_row['pass'] and app_status == (state['app_pid'], 0), row
        finally:
            # Only after recording and evaluating the public result: exact exec rescue.
            if not facts['reaped']:
                os.kill(facts['pid'], signal.SIGKILL)
                until = time.monotonic()+2
                while not n.poll(state['record']):
                    assert time.monotonic() < until
                    time.sleep(.001)
    print('REPORT_FORK_PASS', MODE, OUTCOME, flush=True)


if __name__ == '__main__':
    main()
