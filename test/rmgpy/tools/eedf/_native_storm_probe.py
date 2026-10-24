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

"""Opt-in product qualification: 1,000 serial forks under real raising signals."""
import hashlib
import json
import os
from pathlib import Path
import signal
import sys
import threading
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from rmgpy.tools.eedf import _native_context as c
n = c.native
OWNER = os.getpid()
ACTIVE = None
BLOCK = False
interrupts = post_interrupts = 0
sent = []
stop = threading.Event()


class Storm(Exception):
    pass


def handler(signum, frame):
    global interrupts, post_interrupts
    while frame:
        if frame.f_code.co_name == 'one':
            interrupts += 1
            raise Storm()
        if frame.f_code.co_name == 'post_unwind':
            post_interrupts += 1
            raise Storm()
        frame = frame.f_back
signal.signal(signal.SIGUSR1, handler)


def child_hook():
    if n.is_worker():
        n._mark_hook(ACTIVE)
        if BLOCK:
            time.sleep(1.2)
os.register_at_fork(after_in_child=child_hook)


def sender():
    while not stop.is_set():
        sent.append(time.monotonic())
        os.kill(OWNER, signal.SIGUSR1)
        stop.wait(.002)


def churn():
    while not stop.is_set():
        data = [{str(i): i*i} for i in range(50)]
        with threading.RLock():
            assert len(data) == 50
        stop.wait(.001)


def interpreter(read_end, write_end):
    os.close(read_end)
    assert sum(range(1000)) == 499500
    assert json.loads(json.dumps({'native-child': [1, 2, 3]})) == {'native-child': [1, 2, 3]}
    with threading.RLock():
        assert hashlib.sha256(b'native-product').digest()
    os.write(write_end, b'interpreter-ok')
    os.close(write_end)
    # Return through standard multiprocessing bootstrap and native finish.


def one(owner, process):
    begin = time.monotonic()
    while True:
        try:
            facts = n.info(process._native_record)
            if facts['pid'] == 0:
                process.start()
            if not BLOCK:
                owner.checkpoint()
            if n.poll(process._native_record):
                return n.info(process._native_record)
            if time.monotonic() - begin > 2:
                raise TimeoutError('native fork/reap did not complete')
            time.sleep(.0005)
        except Storm:
            pass


def post_unwind():
    end = time.monotonic() + .22
    while True:
        try:
            if time.monotonic() >= end:
                return
            time.sleep(.002)
        except Storm:
            pass


def main():
    global ACTIVE, BLOCK
    output = Path(sys.argv[1])
    output.mkdir(exist_ok=True, parents=True)
    rows = []
    threads = [threading.Thread(target=sender), threading.Thread(target=churn)]
    for thread in threads:
        thread.start()
    try:
        for index in range(1000):
            BLOCK = index % 10 == 0
            owner = c.Owner(.5)
            rd, wr = os.pipe()
            process = owner.process(interpreter, (rd, wr))
            ACTIVE = process._native_record
            try:
                row = one(owner, process)
                os.close(wr)
                wr = None
                data = os.read(rd, 100)
                row.update(index=index, blocking_hook=BLOCK, child_payload=data.decode(),
                           exitcode=os.waitstatus_to_exitcode(row['wait_status']))
                row['pass'] = (row['thread_error'] == 0 and row['reaped'] == 1
                               and 0 < row['ready_us'] <= row['hook_us'])
                if BLOCK:
                    row['pass'] &= (row['exitcode'] == 124 and row['reason'] == 4
                                    and (row['ended_us']-row['deadline_us'])/1e6 <= .25)
                    try:
                        c.accepted(process, ('ok', True))
                    except c.ScanLeaseExpired:
                        pass
                    else:
                        raise AssertionError('expired storm result accepted')
                else:
                    row['pass'] &= (row['exitcode'] == 0 and row['phase'] == 'normal'
                                    and data == b'interpreter-ok')
                    assert c.accepted(process, ('ok', True)) == ('ok', True)
                rows.append(row)
                with (output/'trials.jsonl').open('a') as stream:
                    stream.write(json.dumps(row)+'\n')
                assert row['pass'], row
                if (index+1) % 100 == 0:
                    print('STORM_FORKS', index+1, 'raised', interrupts, flush=True)
            finally:
                # Never leave a failed qualification's exact child running.
                # This test-only rescue is a failed gate, not a product fallback.
                complete = owner.cleanup()
                if not complete:
                    facts = n.info(ACTIVE)
                    print('STORM_PENDING', json.dumps(facts), flush=True)
                    if facts['pid'] and not facts['reaped']:
                        os.kill(facts['pid'], signal.SIGKILL)
                        end = time.monotonic()+1
                        while not n.poll(ACTIVE) and time.monotonic() < end:
                            time.sleep(.001)
                if wr is not None:
                    os.close(wr)
                os.close(rd)
                assert complete, 'native cleanup failed during qualification'
        post_unwind()
    finally:
        stop.set()
        for thread in threads:
            thread.join(2)
    result = dict(trials=len(rows), unique_pids=len({r['pid'] for r in rows}),
                  normal_interpreter_passes=sum(not r['blocking_hook'] for r in rows),
                  expired_blocked_hooks=sum(r['blocking_hook'] for r in rows),
                  raised_interrupts=interrupts, post_unwind_interrupts=post_interrupts,
                  sent=len(sent), storm_span_s=sent[-1]-sent[0], all_pass=all(r['pass'] for r in rows),
                  interval_s=.002, hook_hold_s=1.2, lease_s=.5)
    assert result['trials'] == result['unique_pids'] == 1000 and interrupts > 0 and post_interrupts > 0
    (output/'summary.json').write_text(json.dumps(result, indent=2))
    print('NATIVE_STORM_RESULT', json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
