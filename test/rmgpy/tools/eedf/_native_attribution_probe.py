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

"""A real handler forks an app child and raises across native PID publication."""
import json
import multiprocessing as mp
import os
from pathlib import Path
import signal
import sys
import time
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from rmgpy.tools.eedf import branches as b
from rmgpy.tools.eedf import _native_context as c
n = c.native
MODE, OUT = sys.argv[1], Path(sys.argv[2])
OWNER = os.getpid()
app = record = None
fired = False


class AppSignal(Exception):
    pass


def app_target():
    assert not n.is_worker()
    # Survive longer than the scan's explicit cancellation and native reap.
    time.sleep(.35)
    (OUT / ('attribution-' + MODE + '-app.json')).write_text(json.dumps(dict(pid=os.getpid(), survived=True)))


def handler(*unused):
    global app
    app = mp.get_context('fork').Process(target=app_target)
    app.start()
    raise AppSignal('honest application handler')
signal.signal(signal.SIGUSR1, handler)


def fire():
    global fired
    if not fired and os.getpid() == OWNER:
        fired = True
        os.kill(OWNER, signal.SIGUSR1)


class Adapter:
    def integrate(self, seed, power=None):
        n._hold_gil(10.)


def main():
    original_process = c.Owner.process
    def process(owner, target, args):
        global record
        child = original_process(owner, target, args)
        record = child._native_record
        return child
    c.Owner.process = process
    original_fork = n.fork
    if MODE == 'before':
        def audit(event, args):
            if event == 'os.fork':
                fire()
        sys.addaudithook(audit)
    elif MODE == 'after':
        def fork(rec):
            pid = original_fork(rec)
            if pid:
                fire()
            return pid
        n.fork = fork
    else:
        raise AssertionError(MODE)
    output = OUT / ('attribution-' + MODE + '-branches.json')
    try:
        b.scan(Adapter(), (0., 1.), output, u_steps=2, electron_density_decades=1,
               electron_density_limits=(1e16, 1e16), processes=1, reference_power=.5, timeout=2., lease_s=.25)
    except AppSignal:
        pass
    else:
        raise AssertionError('real handler exception was swallowed')
    app.join(2)
    assert app.exitcode == 0
    facts = n.info(record)
    assert facts['pid'] != app.pid
    if MODE == 'before':
        assert facts['pid'] == 0 and facts['reaped'] == 0, facts
    else:
        assert facts['pid'] > 0 and facts['reaped'] == 1, facts
        try:
            os.waitpid(facts['pid'], os.WNOHANG)
        except ChildProcessError:
            pass
        else:
            raise AssertionError('worker not reaped by the scan')
    assert not output.exists()
    row = dict(case=MODE, worker=facts, app=app.pid, app_survived=True, passed=True)
    (OUT / ('attribution-' + MODE + '.json')).write_text(json.dumps(row, indent=2))
    print('ATTRIBUTION', json.dumps(row), flush=True)


if __name__ == '__main__':
    main()
