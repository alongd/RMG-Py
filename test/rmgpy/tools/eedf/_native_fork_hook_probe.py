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

"""Application callback forks cannot inherit scan ownership.

Run each mode standalone, without a pytest supervisor: at most owner, worker,
and one application copy. The alarm bounds a failing copied stack; no private
cleanup, polling or signaling assists its public scan verdict.
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

MODE = sys.argv[1]
assert MODE in ('before', 'parent', 'child', 'audit', 'control')
n = c.native
ROOT = os.getpid()
records, held, messages, signals = [], [], [], []
state = dict(forked=False, application_copy=False, app=None, status=None,
             initial=None, cleanup_calls=0, integration_calls=0, elapsed=None)
read_fd, write_fd = os.pipe()
original_process = c.Owner.process


def tracked_process(owner, target, args):
    process = original_process(owner, target, args)
    records.append(process._native_record)
    held.append(owner)
    return process


c.Owner.process = tracked_process


def packet(stage, **extra):
    row = dict(stage=stage, pid=os.getpid(), mode=MODE,
               facts=[n.info(record) for record in records], **extra)
    os.write(write_fd, (json.dumps(row) + '\n').encode())


def profile(frame, event, arg):
    if state['application_copy'] and event == 'call' and frame.f_code is c.Owner.cleanup.__code__:
        state['cleanup_calls'] += 1
        if state['cleanup_calls'] == 100:
            packet('100-foreign-cleanup-retries')


sys.setprofile(profile)


class Messages(logging.Handler):
    def emit(self, record):
        messages.append(record.getMessage())


logging.getLogger(b.__name__).addHandler(Messages())


def application_hook():
    if state['forked'] or not records:
        return
    if MODE == 'child':
        if not n.is_worker():
            return
    elif os.getpid() != ROOT:
        return
    facts = n.info(records[-1])
    if (MODE == 'parent') != bool(facts['pid']):
        return
    state['forked'] = True
    app = os.fork()
    if app == 0:
        os.close(read_fd)
        state.update(application_copy=True, initial=n.info(records[-1]), began=time.monotonic())
        signal.signal(signal.SIGALRM, signal.SIG_DFL)
        signal.alarm(1)
        return
    state['app'] = app
    if MODE in ('before', 'audit', 'child'):
        # A pre-fork copy can wrongly launch a grandchild on the red build.
        # Waiting here keeps even that failure within three live processes.
        state['status'] = os.waitpid(app, 0)
        packet('application-wait', app=app, status=state['status'])


def audit(event, args):
    if state['application_copy'] and event in ('os.kill', 'os.killpg'):
        signals.append([event, list(args)])
    if MODE == 'audit' and event == 'os.fork':
        application_hook()


sys.addaudithook(audit)
if MODE in ('before', 'parent', 'child'):
    os.register_at_fork(**{dict(before='before', parent='after_in_parent',
                               child='after_in_child')[MODE]: application_hook})


class Adapter:
    def integrate(self, seed, power=None):
        state['integration_calls'] += 1
        if state['application_copy']:
            packet('foreign-integration')
        raise ValueError('ordinary failed seed')


with tempfile.TemporaryDirectory(prefix='native-fork-hook-') as directory:
    output = Path(directory) / 'branches.json'
    error = None
    began = time.monotonic()
    try:
        b.scan(Adapter(), (0, 0), output, u_steps=1, electron_density_decades=1,
               processes=1, timeout=2, lease_s=2, reference_power=.5,
               recheck_timeouts=False)
    except BaseException as exc:
        error = exc
        held.append(exc)
    if state['application_copy']:
        packet('returned', error=repr(error), initial=state['initial'],
               cleanup_calls=state['cleanup_calls'], signals=signals,
               messages=messages, elapsed=time.monotonic()-state['began'],
               integration_calls=state['integration_calls'])
        signal.alarm(0)
        os._exit(0)
    sys.setprofile(None)
    if MODE == 'parent':
        state['status'] = os.waitpid(state['app'], 0)
        packet('application-wait', app=state['app'], status=state['status'])
    os.close(write_fd)
    data = b''
    while True:
        part = os.read(read_fd, 65536)
        if not part:
            break
        data += part
    os.close(read_fd)
    events = [json.loads(line) for line in data.decode().splitlines()]
    facts = [n.info(record) for record in records]
    row = dict(mode=MODE, root=ROOT, error=repr(error), facts=facts,
               events=events, artifact_exists=output.exists(),
               elapsed=time.monotonic()-began)
    print('FORK_HOOK_IDENTITY', json.dumps(row), flush=True)
    if MODE != 'control':
        assert not any(event['stage'] == 'foreign-integration' for event in events), 'application copy entered worker bootstrap'
        waits = [event for event in events if event['stage'] == 'application-wait']
        returned = [event for event in events if event['stage'] == 'returned']
        assert len(waits) == 1 and waits[0]['status'] == [waits[0]['app'], 0], 'copied scan did not escape before its alarm'
        assert len(returned) == 1, 'application copy never returned through public scan'
        child = returned[0]
        assert child['integration_calls'] == 0, 'application copy executed the adapter'
        assert child['error'] == "RuntimeError('native ownership belongs to another process or calling thread')", child
        assert child['elapsed'] < 1 and child['cleanup_calls'] <= 1, child
        assert not child['signals'] and not child['messages'], child
        assert len(child['facts']) == 1, child
        for key in ('pid', 'created_us', 'renewed_us', 'reaped', 'wait_status', 'retired'):
            assert child['facts'][0][key] == child['initial'][key], (key, child)
        if MODE in ('before', 'audit'):
            assert child['facts'][0] == child['initial'], 'pre-fork copy wrote the shared record'
    assert error is None and output.exists(), 'application fork damaged the original parent scan'
    assert len(facts) == 1 and facts[0]['phase'] == 'normal' and facts[0]['reaped'] == 1 and facts[0]['wait_status'] == 0, facts
    print('FORK_HOOK_PASS', MODE, flush=True)
