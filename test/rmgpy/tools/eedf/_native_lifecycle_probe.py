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

"""Standalone real owner interference: owner + native worker + ordinary app child.

Run directly, never as a pytest subprocess (that would exceed the process cap).
The app observer's /proc reads are test observations, never worker attribution.
"""
import json
import logging
import mmap
import multiprocessing as mp
import os
from pathlib import Path
import pickle
import signal
import struct
import sys
import threading
import time

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from rmgpy.tools.eedf import branches as b
from rmgpy.tools.eedf import _native_context as c
n = c.native
MODE = sys.argv[1]
OUT = Path(sys.argv[2])
OWNER = os.getpid()
FLAG = mmap.mmap(-1, 32)
REC = APP = None
signals = []
interrupts = 0
release = threading.Event()


class Interrupted(Exception):
    pass


def handler(sig, frame):
    global interrupts
    interrupts += 1
    while frame:
        if frame.f_code.co_name == 'wedge':
            raise Interrupted('real owner signal')
        frame = frame.f_back
signal.signal(signal.SIGUSR1, handler)
signal.signal(signal.SIGUSR2, lambda *unused: release.set())


def state(pid):
    try:
        return Path('/proc/%d/stat' % pid).read_text().rsplit(')', 1)[1].split()[0]
    except FileNotFoundError:
        return 'gone'


def app_job():
    assert not n.is_worker()
    if MODE == 'sigkill':
        os.setsid()
    begin = time.monotonic()
    while not FLAG[0] and time.monotonic() - begin < 5:
        time.sleep(.001)
    assert FLAG[0], 'worker was never admitted'
    worker = struct.unpack_from('q', FLAG, 8)[0]
    assert worker > 0 and worker not in (OWNER, os.getpid())
    initial = {str(pid): state(pid) for pid in (OWNER, os.getpid(), worker)}
    peak = sum(value != 'gone' for value in initial.values())
    if MODE == 'sigkill':
        while not FLAG[1] and time.monotonic() - begin < 5:
            time.sleep(.001)
        assert FLAG[1]
        os.kill(OWNER, signal.SIGKILL)
    attacked = time.monotonic()
    first = None
    heartbeats = 0
    while time.monotonic() - attacked < .55:
        heartbeats += 1
        if state(worker) in ('Z', 'gone') and first is None:
            first = time.monotonic()
        if MODE == 'raising-handler':
            os.kill(OWNER, signal.SIGUSR1)
        time.sleep(.002)
    if MODE == 'sigstop':
        os.kill(OWNER, signal.SIGCONT)
    elif MODE in ('audit-logger', 'finalizer-lock'):
        os.kill(OWNER, signal.SIGUSR2)
    info = n.info(REC)  # Shared state; the parent's private PID cache is not inherited.
    row = dict(case=MODE, owner=OWNER, app_pid=os.getpid(), worker=worker,
               initial_states=initial, measured_peak=peak, app_survived=True,
               first_worker_not_running_s=None if first is None else first-attacked,
               worker_state=state(worker), heartbeats=heartbeats, record=info)
    row['pass'] = (first is not None and peak == 3 and info['ready_us'] > 0
                   and info['reason'] == (5 if MODE == 'sigkill' else 4))
    if MODE != 'sigkill':
        row['pass'] &= (info['ended_us'] - info['deadline_us']) / 1e6 <= .25
    (OUT / ('interference-' + MODE + '-observer.json')).write_text(json.dumps(row, indent=2))
    print('APP_OBSERVER', json.dumps(row), flush=True)
    assert row['pass'], row


class Sink(logging.Handler):
    def emit(self, record):
        pass
sink = Sink()
log = logging.getLogger('native.app')
log.addHandler(sink)
log.setLevel(logging.INFO)


def slow_restore():
    time.sleep(.65)
    return 'restored'


class Slow:
    def __reduce__(self):
        return slow_restore, ()


class Finalizer:
    def __del__(self):
        log.warning('honest finalizer logger')


class Adapter:
    def integrate(self, seed, power=None):
        if MODE == 'sigkill':
            os.setsid()
        FLAG[1] = 1
        n._hold_gil(10.)
        raise AssertionError('watchdog did not terminate native integration')


def wedge():
    if MODE == 'raising-handler':
        until = time.monotonic() + .65
        while time.monotonic() < until:
            try:
                time.sleep(.02)
            except Interrupted:
                pass
    elif MODE == 'audit-logger':
        sys.audit('eedf.native_probe')
    elif MODE == 'finalizer-lock':
        value = Finalizer()
        del value
    elif MODE == 'slow-unpickle':
        assert pickle.loads(pickle.dumps(Slow())) == 'restored'
    elif MODE == 'sigstop':
        os.kill(OWNER, signal.SIGSTOP)
    elif MODE == 'sigkill':
        time.sleep(10)
    else:
        raise AssertionError(MODE)


def main():
    global APP, REC
    original_process = c.Owner.process
    def process(owner, target, args):
        global APP, REC
        child = original_process(owner, target, args)
        REC = child._native_record
        APP = mp.get_context('fork').Process(target=app_job)
        APP.start()
        return child
    c.Owner.process = process
    original_launch = c.NativePopen._launch
    def launch(popen, child):
        original_launch(popen, child)
        struct.pack_into('q', FLAG, 8, popen.pid)
        held = threading.Event()
        if MODE in ('audit-logger', 'finalizer-lock'):
            def hold():
                with sink.lock:
                    held.set()
                    release.wait(3)
            thread = threading.Thread(target=hold)
            thread.start()
            assert held.wait(1)
        FLAG[0] = 1
        wedge()
    c.NativePopen._launch = launch
    def audit(event, args):
        if os.getpid() == OWNER:
            if event == 'os.kill':
                signals.append(list(args))
            if event == 'eedf.native_probe' and MODE == 'audit-logger':
                log.warning('honest audit logger')
    sys.addaudithook(audit)
    path = OUT / ('interference-' + MODE + '-branches.json')
    try:
        b.scan(Adapter(), (0., 1.), path, u_steps=2, electron_density_decades=1,
               electron_density_limits=(1e16, 1e16), processes=1, reference_power=.5, timeout=2., lease_s=.25)
    except c.ScanLeaseExpired:
        pass
    else:
        raise AssertionError('stalled public scan did not cancel')
    APP.join(3)
    assert APP.exitcode == 0, APP.exitcode
    info = n.info(REC)
    assert info['reaped'] == 1 and os.waitstatus_to_exitcode(info['wait_status']) == 124, info
    assert all(item[0] == OWNER for item in signals), signals
    assert not path.exists()
    row = json.loads((OUT / ('interference-' + MODE + '-observer.json')).read_text())
    row.update(owner_record=info, owner_signal_audit=signals, interrupts=interrupts,
               app_exitcode=APP.exitcode)
    assert row['pass']
    (OUT / ('interference-' + MODE + '.json')).write_text(json.dumps(row, indent=2))
    print('INTERFERENCE', MODE, 'PASS', flush=True)


if __name__ == '__main__':
    main()
