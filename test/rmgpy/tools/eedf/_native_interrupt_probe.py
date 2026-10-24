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

"""Fresh-owner regression for native ownership across real raising signals."""
import ctypes
import json
import logging
import os
from pathlib import Path
import signal
import sys
import tempfile
import threading
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from rmgpy.tools.eedf import _native_context as c
n = c.native
MODE = sys.argv[1]
record = None


def child_hook():
    if n.is_worker():
        if MODE == 'unreadable':
            assert ctypes.CDLL(None).prctl(4, 0, 0, 0, 0) == 0
        n._mark_hook(record)
        threading.Event().wait()
os.register_at_fork(after_in_child=child_hook)
from rmgpy.tools.eedf import branches as b

owners, records, sent = [], [], []
interrupts = post_interrupts = 0
post_loop = False
unwound = threading.Event()
stop = threading.Event()
OriginalOwner = c.Owner


class TrackedOwner(OriginalOwner):
    def __init__(self, lease_s):
        super().__init__(lease_s)
        owners.append(self)

    def process(self, target, args):
        global record
        child = super().process(target, args)
        record = child._native_record
        records.append(record)
        return child
c.Owner = TrackedOwner


class AppSignal(Exception):
    pass


def handler(signum, frame):
    global interrupts, post_interrupts
    while frame:
        if frame.f_code is b.scan.__code__:
            interrupts += 1
            raise AppSignal('during public scan')
        frame = frame.f_back
    if post_loop:
        post_interrupts += 1
        raise AppSignal('after scan unwind')
signal.signal(signal.SIGUSR1, handler)


def storm():
    begin = time.monotonic()
    after = None
    while not stop.is_set():
        now = time.monotonic()
        if MODE == 'storm':
            if now - begin >= .3:
                break
        elif unwound.is_set():
            if after is None:
                after = now
            if now - after >= .2:
                break
        if now - begin > 2:
            break
        sent.append(now)
        os.kill(os.getpid(), signal.SIGUSR1)
        time.sleep(.002)


original_fork = n.fork
sender = None


def fork(rec):
    global sender
    pid = original_fork(rec)
    if pid:
        if MODE.startswith('storm'):
            sender = threading.Thread(target=storm)
            sender.start()
            time.sleep(.05)
        os.kill(os.getpid(), signal.SIGUSR1)
    return pid
n.fork = fork
original_cleanup = OriginalOwner.cleanup
cleanup_calls = 0


def cleanup(owner, seconds=5.):
    global cleanup_calls
    cleanup_calls += 1
    if MODE in ('cleanup-before', 'bookkeeping') and cleanup_calls == 1:
        os.kill(os.getpid(), signal.SIGUSR1)
    result = original_cleanup(owner, seconds)
    if MODE == 'cleanup-after' and cleanup_calls == 1:
        os.kill(os.getpid(), signal.SIGUSR1)
    return result
OriginalOwner.cleanup = cleanup
injected = False


def trace(frame, event, arg):
    global injected
    if (MODE == 'bookkeeping' and not injected and event == 'line'
            and frame.f_code is b._integrate_owned.__code__):
        line = Path(b.__file__).read_text().splitlines()[frame.f_lineno - 1].strip()
        if line == 'interrupted = exc':
            injected = True
            os.kill(os.getpid(), signal.SIGUSR1)
    return trace


class Adapter:
    def integrate(self, seed, power=None):
        raise AssertionError('blocked child hook must not reach integration')


def main():
    global post_loop
    sys.settrace(trace)
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / 'branches.json'
        started = time.monotonic()
        kept = None
        try:
            b.scan(Adapter(), (0., 0.), path, u_steps=1, electron_density_decades=1,
                   processes=1, timeout=2., lease_s=.12, reference_power=.5)
        except AppSignal as exc:
            kept = exc  # Keep the owning frame/traceback alive through the verdict.
        else:
            raise AssertionError('application signal was lost')
        finally:
            sys.settrace(None)
        returned = time.monotonic()
        unwound.set()
        post_loop = True
        end = returned + (0 if MODE == 'exit-immediately' else .36)
        while True:
            try:
                if time.monotonic() >= end:
                    post_loop = False
                    break
                time.sleep(.002)
            except AppSignal:
                pass
        stop.set()
        if sender is not None:
            sender.join(1.)
            assert not sender.is_alive()
        before_reap = n.info(record)
        state = Path('/proc/%d/stat' % before_reap['pid'])
        # Observation must not finish product cleanup. In particular, a zombie
        # is a failure even if its retained Owner could be cleaned up privately.
        assert kept is not None and before_reap['reaped'] == 1 and not state.exists(), before_reap
        facts = n.info(record)
        assert facts['pid'] > 0 and facts['reaped'] == 1 and facts['reason'] in (4, 6), facts
        assert facts['ready_us'] > 0 and facts['ended_us'] > 0, facts
        assert returned - started < .75
        assert not path.exists()
        if MODE == 'bookkeeping':
            assert injected and interrupts == 3
            assert before_reap['reason'] == 6  # The guard cancels before lease expiry.
        if MODE.startswith('cleanup-'):
            assert cleanup_calls >= 2 and interrupts == 2
        if MODE.startswith('storm'):
            assert len(sent) >= 50 and sent[-1] - sent[0] >= .18
            if MODE == 'storm-unbounded':
                assert post_interrupts >= 30 and .18 <= sent[-1] - returned < .3
        row = dict(mode=MODE, interrupts=interrupts, post_interrupts=post_interrupts,
                   cleanup_calls=cleanup_calls, injected=injected, cleanup_s=returned-started,
                   sent=len(sent), storm_span=sent[-1]-sent[0] if sent else 0,
                   last_signal_after_unwind=sent[-1]-returned if sent else 0,
                   before_reap=before_reap, record=facts, artifact_absent=True)
        print('NATIVE_INTERRUPT', json.dumps(row), flush=True)


if __name__ == '__main__':
    main()
    if MODE == 'exit-immediately':
        os._exit(0)
