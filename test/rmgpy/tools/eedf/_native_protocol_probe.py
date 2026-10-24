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

"""Fresh-process native protocol assertions; called by nativeWorkerTest.py."""
import hashlib
import json
import logging
import os
from pathlib import Path
import signal
import sys
import tempfile
import threading
import time

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from rmgpy.tools.eedf import _native_context as c
n = c.native
MODE = sys.argv[1]
record = None
hook_read, hook_write = os.pipe() if MODE == 'blocking-hook' else (None, None)


def hook():
    if MODE.startswith('exec') and n.is_worker():
        os.execv('/usr/bin/sleep', ['sleep', '30'])
    if MODE == 'blocking-hook' and n.is_worker():
        # At-fork exceptions are ignored by CPython. Preserve every invocation
        # in order so the parent sees the first readiness value even if a later
        # hook overwrites the native timestamp after a swallowed assertion.
        os.write(hook_write, n.info(record)['ready_us'].to_bytes(8, 'little'))
        n._mark_hook(record)
        assert n.info(record)['ready_us'] > 0
        threading.Event().wait()
os.register_at_fork(after_in_child=hook)
# Register the application's blocking hook before the scanner's child hook.
from rmgpy.tools.eedf import branches as b


def target():
    if MODE in ('gil', 'healthy', 'lost-popen', 'expiry'):
        n._hold_gil(.35 if MODE == 'healthy' else 10.)
    elif MODE == 'bootstrap-finalizer':
        from multiprocessing.util import Finalize
        Finalize(None, threading.Event().wait, exitpriority=0)


def await_state(predicate, bound=1.):
    until = time.monotonic() + bound
    while not predicate():
        assert time.monotonic() < until, 'native state observation timed out'
        time.sleep(.001)


def check_reaped(process, expected):
    process.join(1.)
    facts = n.info(process._native_record)
    assert facts['reaped'] == 1, facts
    assert process.exitcode == expected, facts
    return facts


def protocol():
    global record, hook_write
    owner = c.Owner(.12, trace=True)
    process = owner.process(target, ())
    record = process._native_record
    if MODE in ('completion', 'race', 'expiry'):
        n._test_config(record, 'finish_gate', 1)
    if MODE in ('missing-handshake', 'create-error'):
        n._test_config(record, 'suppress_ready' if MODE == 'missing-handshake' else 'create_error',
                       1 if MODE == 'missing-handshake' else 11)
    original_fork = n.fork
    if MODE == 'lost-popen':
        def interrupted_fork(rec):
            pid = original_fork(rec)
            if pid:
                os.kill(os.getpid(), signal.SIGUSR1)
            return pid
        def interrupt(*unused):
            raise InterruptedError('return assignment interrupted')
        signal.signal(signal.SIGUSR1, interrupt)
        n.fork = interrupted_fork
    started = time.monotonic()
    try:
        try:
            process.start()
        except c.NativeWorkerStartupError:
            assert MODE in ('missing-handshake', 'create-error')
            assert owner.cleanup()
            facts = n.info(record)
            assert facts['reaped'] == 1 and facts['reason'] == 7, facts
        except InterruptedError:
            assert MODE == 'lost-popen' and process._popen is None
            assert n.info(record)['pid'] > 0
            assert owner.cleanup()
            assert n.info(record)['reaped'] == 1
        else:
            assert MODE not in ('missing-handshake', 'create-error', 'lost-popen')
            if MODE == 'healthy':
                while process.is_alive():
                    owner.checkpoint()
                    time.sleep(.01)
                facts = check_reaped(process, 0)
                assert facts['phase'] == 'normal' and time.monotonic() - started >= .35
                assert len(n.trace(record)) > 10
            elif MODE in ('blocking-hook', 'gil', 'expiry', 'bootstrap-finalizer'):
                facts = check_reaped(process, 124)
                assert facts['reason'] == 4 and time.monotonic() - started < .37, facts
                if MODE == 'blocking-hook':
                    os.close(hook_write)
                    hook_write = None
                    first_hook = os.read(hook_read, 8)
                    print('FIRST_HOOK_READY', first_hook.hex(), flush=True)
                    assert len(first_hook) == 8, 'first child hook witness missing'
                    assert int.from_bytes(first_hook, 'little') > 0, 'watchdog unarmed at first child hook'
                    assert 0 < facts['ready_us'] <= facts['hook_us'] < facts['ended_us']
                try:
                    owner.checkpoint()
                except c.ScanLeaseExpired:
                    pass
                else:
                    raise AssertionError('expired record revived')
                try:
                    c.accepted(process, ('ok', 1))
                except c.ScanLeaseExpired:
                    pass
                else:
                    raise AssertionError('expired result accepted')
            elif MODE in ('completion', 'race'):
                await_state(lambda: n.info(record)['finish_entered_us'])
                if MODE == 'race':
                    until = n.info(record)['deadline_us'] / 1e6
                    while time.monotonic() < until:
                        pass
                n._test_config(record, 'finish_gate', 0)
                process.join(1)
                facts = n.info(record)
                assert facts['reaped'] == 1
                if facts['phase'] == 'normal':
                    assert c.accepted(process, ('ok', 1)) == ('ok', 1)
                else:
                    assert MODE == 'race' and process.exitcode == 124
                    try:
                        c.accepted(process, ('ok', 1))
                    except c.ScanLeaseExpired:
                        pass
                    else:
                        raise AssertionError('expiry accepted')
            elif MODE == 'external-reap':
                pid, status = os.waitpid(process.pid, 0)
                assert pid == process.pid and status == 0
                process.join(1)
                assert n.info(record)['reaped'] == 2
                try:
                    c.accepted(process, ('ok', 1))
                except c.NativeWorkerStatusError:
                    pass
                else:
                    raise AssertionError('external reap accepted')
            elif MODE == 'reap-cache':
                await_state(lambda: n.info(record)['phase'] == 'normal')
                original = n.poll
                def interrupted_poll(rec):
                    result = original(rec)
                    if result:
                        raise InterruptedError('after exact wait, before Python store')
                    return result
                n.poll = interrupted_poll
                try:
                    process.join(1)
                except InterruptedError:
                    pass
                finally:
                    n.poll = original
                facts = check_reaped(process, 0)
                assert n.poll(record) == n.poll(record) == 1
                assert c.accepted(process, ('ok', 1)) == ('ok', 1)
            elif MODE == 'accept-before-reap':
                await_state(lambda: n.info(record)['phase'] == 'normal')
                assert n.info(record)['reaped'] == 0
                try:
                    c.accepted(process, ('ok', 1))
                except c.NativeWorkerStatusError:
                    pass
                else:
                    raise AssertionError('accepted without exact reap')
                check_reaped(process, 0)
            elif MODE == 'exited-owner-stall':
                await_state(lambda: n.info(record)['phase'] == 'normal')
                time.sleep(.16)
                check_reaped(process, 0)
                try:
                    owner.checkpoint()
                except c.ScanLeaseExpired:
                    pass
                else:
                    raise AssertionError('stalled owner accepted after child normal exit')
            elif MODE == 'foreign-thread':
                await_state(lambda: n.info(record)['phase'] == 'normal')
                observed = []
                thread = threading.Thread(target=lambda: observed.append(n.poll(record)))
                thread.start(); thread.join()
                assert observed == [0] and n.info(record)['reaped'] == 0
                # Dropping the last Python Owner reference on another thread
                # may cancel, but must not move exact reaping to that thread.
                holder = [owner]
                owner = None
                thread = threading.Thread(target=holder.clear)
                thread.start(); thread.join()
                assert n.info(record)['reaped'] == 0
                facts = check_reaped(process, 0)
            elif MODE == 'bad-envelope':
                facts = check_reaped(process, 0)
                for envelope in ({}, ('ok',), ('other', 1), ('lifecycle', ValueError())):
                    try:
                        c.accepted(process, envelope)
                    except c.NativeWorkerStatusError:
                        pass
                    else:
                        raise AssertionError('incomplete or invalid protocol envelope accepted')
            else:
                raise AssertionError(MODE)
        print('NATIVE_PROTOCOL', MODE, json.dumps(n.info(record)), flush=True)
    finally:
        n.fork = original_fork
        for fd in (hook_read, hook_write):
            if fd is not None:
                os.close(fd)
        if owner is not None:
            complete = owner.cleanup()
            if not complete:
                # Test-only rescue for intentionally broken watchdog mutants.
                # The C record is the exact unreaped child of this test owner.
                facts = n.info(record)
                if facts['pid'] and not facts['reaped']:
                    os.kill(facts['pid'], signal.SIGKILL)
                    await_state(lambda: n.poll(record), bound=1.)
            assert complete, 'native cleanup did not reap the worker'


def slow_restore():
    time.sleep(.2)
    return {'u': .5, 'n_e': 1e16}


class SlowResult:
    def __reduce__(self):
        return slow_restore, ()


class Adapter:
    def integrate(self, seed, power=None):
        return SlowResult()


def public_artifact():
    expected = {'artifact-lease': c.ScanLeaseExpired,
                'artifact-platform': c.NativeWorkerUnsupported,
                'artifact-startup': c.NativeWorkerStartupError}[MODE]
    if MODE == 'artifact-platform':
        original = n.check_platform
        n.check_platform = lambda: dict(original(), supported=False)
    if MODE == 'artifact-startup':
        original = c.Owner.process
        def broken_start(owner, target, args):
            process = original(owner, target, args)
            n._test_config(process._native_record, 'create_error', 11)
            return process
        c.Owner.process = broken_start
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / 'branches.json'
        for existing in (False, True):
            if existing:
                path.write_bytes(b'previous accepted artifact\n')
                before = (path.stat().st_mtime_ns, hashlib.sha256(path.read_bytes()).hexdigest())
            try:
                b.scan(Adapter(), (0., 1.), path, u_steps=2, electron_density_decades=1,
                       electron_density_limits=(1e16, 1e16), processes=1, reference_power=.5, timeout=2., lease_s=.08)
            except expected:
                pass
            else:
                raise AssertionError('public lifecycle refusal missing')
            if existing:
                assert (path.stat().st_mtime_ns, hashlib.sha256(path.read_bytes()).hexdigest()) == before
            else:
                assert not path.exists()
            assert list(Path(directory).iterdir()) == ([path] if existing else [])
    print('PUBLIC_ARTIFACT', MODE, 'PASS', flush=True)


def unsupported_context():
    import multiprocessing as mp
    method = MODE.split('-', 1)[1]
    if method == 'extension':
        c.native = None
    else:
        mp.set_start_method(method)
    before = mp.get_start_method()
    try:
        c.Owner(60.)
    except c.NativeWorkerUnsupported:
        pass
    else:
        raise AssertionError('unsupported native mode was admitted')
    assert mp.get_start_method() == before
    # The ordinary in-process route stays usable and selects no context.
    result = b._integrate_many(Adapter(), [b.Seed(0., 1e16)], 1, None)
    assert result[0][0] == 'ok' and isinstance(result[0][1], SlowResult)
    assert mp.get_start_method() == before
    print('NATIVE_PROTOCOL', MODE, 'PASS', flush=True)


def exec_replacement():
    """Exec ends the guarantee; refuse results and report the unreaped PID."""
    b.WORKER_KILL_TIMEOUT = .5
    records, messages = [], []
    original = c.Owner.process
    original_cleanup = c.Owner.cleanup
    cleanup_calls = 0
    warning_interrupted = False

    def interrupted(*unused):
        raise InterruptedError('real application signal during exec cleanup')
    signal.signal(signal.SIGUSR1, interrupted)

    def cleanup(owner, seconds=5.):
        nonlocal cleanup_calls
        cleanup_calls += 1
        result = original_cleanup(owner, seconds)
        if MODE == 'exec-retry' and cleanup_calls % 2:
            os.kill(os.getpid(), signal.SIGUSR1)
        return result
    c.Owner.cleanup = cleanup

    def tracked(owner, target, args):
        process = original(owner, target, args)
        records.append(process._native_record)
        return process

    class CaptureLog(logging.Handler):
        def emit(self, record):
            nonlocal warning_interrupted
            messages.append(record.getMessage())
            if MODE == 'exec-warning' and not warning_interrupted:
                warning_interrupted = True
                os.kill(os.getpid(), signal.SIGUSR1)

    c.Owner.process = tracked
    handler = CaptureLog()
    logger = logging.getLogger(b.__name__)
    logger.addHandler(handler)
    try:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'branches.json'
            for existing in (False, True):
                if existing:
                    path.write_bytes(b'previous accepted artifact\n')
                    before = (path.stat().st_mtime_ns, hashlib.sha256(path.read_bytes()).hexdigest())
                try:
                    b.scan(Adapter(), (0., 0.), path, u_steps=1, electron_density_decades=1,
                           processes=1, reference_power=.5, timeout=2., lease_s=.1)
                except c.NativeWorkerStatusError as exc:
                    message = str(exc)
                else:
                    raise AssertionError('exec replacement returned results')
                facts = n.info(records[-1])
                elapsed = time.monotonic() - facts['ended_us'] / 1e6
                assert facts['ended_us'] > 0 and facts['reaped'] == 0, facts
                assert facts['phase'] == 'cancelled', facts
                assert os.readlink('/proc/%d/exe' % facts['pid']) == '/usr/bin/sleep'
                assert 'native cleanup deadline reached' in message, message
                assert 'pending child PIDs [%d]' % facts['pid'] in message, message
                assert message in messages, messages
                assert .5 <= elapsed < .75, elapsed  # finally must not spend a second budget
                if existing:
                    assert (path.stat().st_mtime_ns, hashlib.sha256(path.read_bytes()).hexdigest()) == before
                else:
                    assert not path.exists()
                print('EXEC_PENDING', json.dumps(dict(facts=facts, cleanup_elapsed_s=elapsed,
                                                      artifact_unchanged=True)), flush=True)
                # Test-only teardown of the exact unreaped child; product never signals it.
                os.kill(facts['pid'], signal.SIGKILL)
                await_state(lambda: n.poll(records[-1]))
                assert n.info(records[-1])['reaped'] == 1
    finally:
        c.Owner.process = original
        c.Owner.cleanup = original_cleanup
        logger.removeHandler(handler)
        for rec in records:
            facts = n.info(rec)
            if facts['pid'] and not facts['reaped']:
                os.kill(facts['pid'], signal.SIGKILL)
                await_state(lambda: n.poll(rec))
    if MODE == 'exec-retry':
        assert cleanup_calls == 4, cleanup_calls
    if MODE == 'exec-warning':
        assert warning_interrupted
    print('NATIVE_PROTOCOL', MODE, 'PASS', flush=True)


if __name__ == '__main__':
    if MODE.startswith('unsupported-'):
        unsupported_context()
    elif MODE.startswith('exec'):
        exec_replacement()
    elif MODE.startswith('artifact-'):
        public_artifact()
    else:
        for _ in range(20 if MODE in ('race', 'completion') else 1):
            protocol()
