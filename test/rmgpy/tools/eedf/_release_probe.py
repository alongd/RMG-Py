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

"""Real application hooks/signals for the resumable-release regressions.

Standalone: owner + worker, with pytest as the third process. The lock case
runs in pytest itself so its application child does not become process four.
"""
import ctypes
import dis
import json
import logging
import os
from pathlib import Path
import signal
import sys
import threading
import time


def run(mode, output, in_pytest=False):
    owner = os.getpid()
    ready_r, ready_w = os.pipe()
    cleanup_event = threading.Event()
    app_entered = threading.Event()
    local = threading.local()
    worker = start = entered = returned = None
    interrupts = 0
    post_interrupts = 0
    injected = False
    attacks = []
    warnings, errors, reaps, sent = [], [], [], []
    application = []
    hook_time = 1.2 if mode.startswith(('storm', 'lease-')) else 5

    def before():
        if os.getpid() == owner and getattr(local, 'application', False):
            app_entered.set()
            time.sleep(2.5)

    def child_hook():
        if os.getpid() != owner:
            if getattr(local, 'application', False):
                os._exit(0)
            if mode == 'triple-unreadable':
                assert ctypes.CDLL(None).prctl(4, 0, 0, 0, 0) == 0
            os.write(ready_w, b'x')
            time.sleep(hook_time)

    if in_pytest:
        sys._eedf_application_before_fork = before
        sys._eedf_application_after_fork = child_hook
    else:
        os.register_at_fork(before=before, after_in_child=child_hook)
    from rmgpy.tools.eedf import branches as b
    original_timeout = b.WORKER_RECLAIM_TIMEOUT
    original_lease = getattr(b, 'RUN_LEASE', None)
    if mode.startswith('lease-'):
        b.RUN_LEASE = .25
    b.WORKER_RECLAIM_TIMEOUT = .5 if mode.startswith('storm') else 1

    class AppSignal(Exception):
        pass

    def handler(signum, frame):
        nonlocal interrupts, post_interrupts, returned
        interrupts += 1
        # Profiling itself can be interrupted and disabled by a raising signal.
        # Determine whether scan is still on the stack instead of depending on
        # a profile-return callback to mark the observation phase.
        while frame is not None:
            if frame.f_code is b.scan.__code__:
                raise AppSignal('signal %d' % interrupts)
            frame = frame.f_back
        if mode == 'storm-unbounded':
            if returned is None:
                returned = time.monotonic()
            post_interrupts += 1
            raise AppSignal('post-unwind signal %d' % post_interrupts)

    previous_handler = signal.signal(signal.SIGUSR1, handler)

    class Capture(logging.Handler):
        def filter(self, record):
            # Before the handler lock: logging.shutdown must not wait for this
            # daemon-only application filter and accidentally save its warning.
            if mode == 'exit-immediately' and threading.current_thread() is not threading.main_thread():
                time.sleep(.25)
            return True

        def emit(self, record):
            nonlocal injected
            message = record.getMessage()
            if message.startswith('releasing child pid '):
                if mode == 'warning-before' and not injected:
                    injected = True
                    os.kill(owner, signal.SIGUSR1)
                if mode == 'exit-immediately':
                    # A slow application sink makes consumer-only delivery
                    # deterministically lose its warning at interpreter exit.
                    sys.stderr.write(message + '\n')
                    sys.stderr.flush()
                warnings.append(message)
                if mode.startswith('lock'):
                    assert app_entered.is_set()
                    if mode == 'lock' or threading.current_thread() is threading.main_thread():
                        os.kill(owner, signal.SIGUSR1)  # original r273: every emit
            elif record.levelno >= logging.ERROR:
                errors.append(message)

    capture = Capture()
    logger = logging.getLogger(b.__name__)
    previous_level = logger.level
    logger.setLevel(logging.WARNING)
    logger.addHandler(capture)
    real_reap = b._kill_and_reap

    def reap(pid):
        reaps.append([pid, start.pid(), threading.current_thread().name])
        return real_reap(pid)

    b._kill_and_reap = reap
    logging_source = Path(logging.__file__).read_text().splitlines()
    source = Path(b.__file__).read_text().splitlines()
    opcode_target = None
    if mode.startswith('opcode-'):
        targets = {
            'publish': (b._reclaim, 'STORE_ATTR', 'releases'),
            'queued': (b._defer_reap, 'STORE_ATTR', 'queued'),
            'ready': (b._defer_reap, 'LOAD_ATTR', 'ready'),
            'handed': (b._release_candidate, 'STORE_ATTR', 'handed'),
            'warned': (b._warn_release, 'STORE_ATTR', 'warned'),
            'reclaimed': (b._reclaim, 'STORE_ATTR', 'reclaimed'),
            'delete': (b._integrate_many, 'DELETE_SUBSCR', None),
            'closed': (b._close_lifeline, 'POP_TOP', None),
            'lease': (b._close_reaper_lease, 'STORE_ATTR', 'closed'),
        }
        function, opname, argval = targets[mode[7:]]
        matches = [i.offset for i in dis.get_instructions(function)
                   if i.opname == opname and (argval is None or i.argval == argval)]
        assert matches, mode
        opcode_target = (function.__code__, matches[-1])

    def trace(frame, event, arg):
        nonlocal worker, start, injected, entered
        if os.getpid() != owner:
            return None
        if mode.startswith('lease-') and event == 'line':
            line = source[frame.f_lineno - 1].strip() if frame.f_code.co_filename == b.__file__ else ''
            attack = None
            if mode == 'lease-deadline-store':
                if (frame.f_code is b._reclaim.__code__ and not attacks
                        and line == 'start.deadline = time.monotonic() + WORKER_RECLAIM_TIMEOUT'):
                    attack = 'before first deadline store'
                elif line == 'interrupted = None' and len(attacks) == 1:
                    attack = 'outer cleanup entry'
            elif mode == 'lease-cleanup-entry':
                if (frame.f_code.co_name == 'close' and not attacks
                        and frame.f_code.co_filename.endswith('/multiprocessing/connection.py')):
                    attack = 'inner cleanup connection close'
                elif line == 'interrupted = None' and len(attacks) == 1:
                    attack = 'outer cleanup entry'
            elif mode == 'lease-retained-exception':
                if (frame.f_code is b._reclaim.__code__ and not attacks
                        and line == 'if start.deadline is None:'):
                    attack = 'before first deadline guard'
                elif line == 'interrupted = None' and len(attacks) == 1:
                    attack = 'outer finally before retry loop'
            if attack:
                if entered is None:
                    entered = time.monotonic()
                attacks.append(attack)
                os.kill(owner, signal.SIGUSR1)
        if (opcode_target and event == 'opcode' and not injected
                and frame.f_code is opcode_target[0] and frame.f_lasti == opcode_target[1]):
            injected = True
            attacks.append([frame.f_code.co_name, frame.f_lasti])
            os.kill(owner, signal.SIGUSR1)
        if event == 'line' and frame.f_code is b._WorkerPopen._launch.__code__ and worker is None:
            pid = getattr(frame.f_locals.get('self'), 'pid', None)
            if pid:
                worker, start = pid, b._STARTING.start
                assert os.read(ready_r, 1) == b'x'
                os.kill(owner, signal.SIGUSR1)
        if (event == 'line' and frame.f_code is b._defer_reap.__code__
                and mode.startswith('handoff') and not injected):
            line = source[frame.f_lineno - 1].strip()
            before_put = line == 'with _LIFELINES_LOCK:' or line == 'record.service.queue.put(record)'
            after_put = line.startswith('if _REAPER is not None') or line == 'if not record.service.ready:'
            if (mode == 'handoff-before' and before_put) or (mode == 'handoff-after' and after_put):
                injected = True
                os.kill(owner, signal.SIGUSR1)
        if event == 'line' and mode == 'retry-bookkeeping':
            line = source[frame.f_lineno - 1].strip() if frame.f_code.co_filename == b.__file__ else ''
            if frame.f_code is b._defer_reap.__code__ and line == 'record.service.queue.put(record)' and len(attacks) < 2:
                attacks.append('queue put interrupted')
                os.kill(owner, signal.SIGUSR1)
            if frame.f_code is b._integrate_many.__code__ and line == 'interrupted = exc':
                attacks.append('retry bookkeeping interrupted')
                os.kill(owner, signal.SIGUSR1)
        if (mode == 'warning-lock' and event == 'line' and not injected
                and frame.f_code is logging.Handler.handle.__code__
                and frame.f_locals.get('self') is capture
                and logging_source[frame.f_lineno - 1].strip() == 'try:'):
            injected = True
            attacks.append('Handler.handle after acquire before try')
            os.kill(owner, signal.SIGUSR1)
        return trace

    def profile(frame, event, arg):
        nonlocal entered
        if os.getpid() != owner:
            return
        if mode.startswith('lease-'):
            sys.settrace(trace)
            cursor = frame
            while cursor is not None:
                cursor.f_trace = trace
                cursor = cursor.f_back
        if opcode_target and event == 'call':
            sys.settrace(trace)
            cursor = frame
            while cursor is not None:
                if cursor.f_code is opcode_target[0]:
                    cursor.f_trace = trace
                    cursor.f_trace_opcodes = True
                    break
                cursor = cursor.f_back
        if mode == 'retry-bookkeeping' and frame.f_code is b._defer_reap.__code__:
            sys.settrace(trace)
            cursor = frame
            while cursor is not None:
                cursor.f_trace = trace
                cursor = cursor.f_back
        if (mode == 'warning-lock' and event == 'call'
                and frame.f_code is logging.Handler.handle.__code__):
            sys.settrace(trace)
            frame.f_trace = trace
        if event != 'call':
            return
        if frame.f_code is b._defer_reap.__code__ and mode.startswith('handoff'):
            # A raising trace callback disables itself; re-arm the second attack.
            sys.settrace(trace)
            frame.f_trace = trace
        if frame.f_code is b._reclaim.__code__ and entered is None:
            entered = time.monotonic()
            cleanup_event.set()

    def send():
        cleanup_event.wait()
        if mode.startswith('storm'):
            time.sleep(.1)
            until = time.monotonic() + .3
            while (returned is None or time.monotonic() < returned + .2
                   if mode == 'storm-unbounded' else time.monotonic() < until):
                sent.append(time.monotonic())
                os.kill(owner, signal.SIGUSR1)
                time.sleep(.002)
        else:
            for delay in (.2, .2):
                time.sleep(delay)
                os.kill(owner, signal.SIGUSR1)

    def fork_application():
        cleanup_event.wait()
        time.sleep(.2)
        local.application = True
        pid = os.fork()
        if pid == 0:
            os._exit(0)
        application.append(os.waitpid(pid, 0))

    sender = None
    if mode.startswith(('triple', 'storm')):
        sender = threading.Thread(target=send, daemon=True)
    elif mode.startswith('lock'):
        sender = threading.Thread(target=fork_application)
    if sender is not None:
        sender.start()

    class Slow:
        def integrate(self, seed, power=None):
            time.sleep(30)

    def state():
        try:
            return Path('/proc/%d/stat' % worker).read_text().rsplit(')', 1)[1].split()[0]
        except FileNotFoundError:
            return 'gone'

    begin = time.monotonic()
    previous_trace, previous_profile = sys.gettrace(), sys.getprofile()
    sys.setprofile(profile)
    sys.settrace(trace)
    try:
        try:
            b.scan(Slow(), (0, 0), output, u_steps=1, electron_density_decades=1,
                   processes=1, timeout=30, reference_power=.5, recheck_timeouts=False)
            returned = time.monotonic()
            outcome = 'return'
        except BaseException as exc:
            returned = time.monotonic()
            outcome = type(exc).__name__ + ': ' + str(exc)
            kept_exception = exc  # applications may retain errors and their tracebacks
            if mode == 'exit-immediately':
                print(json.dumps({'worker': worker, 'outcome': outcome}), flush=True)
                sys.exit(0)  # no wait or sleep between scan raising and interpreter exit
            if mode == 'storm-unbounded':
                # The application keeps catching real raising signals after
                # scan unwinds. Measurement/JSON waits until this storm ends.
                while True:
                    try:
                        if time.monotonic() >= returned + .22:
                            break
                        time.sleep(.01)
                    except AppSignal:
                        pass
        finally:
            sys.setprofile(previous_profile)
            sys.settrace(previous_trace)
        assert returned is not None
        row = dict(mode=mode, worker=worker, outcome=outcome, interrupts=interrupts,
                   injected=injected, cleanup_elapsed=returned-entered,
                   deadline=b.WORKER_RECLAIM_TIMEOUT, at_return=state(),
                   report_at_return=start.pid(), app_hook=app_entered.is_set())
        # Observe only; no waitpid on the worker until after recording evidence.
        until = begin + hook_time + 1
        while time.monotonic() < until and state() != 'gone':
            try:
                time.sleep(.005)
            except AppSignal:
                # The approved fallback can emit after scan raises. Its real
                # application signal must not discard the already-timed row.
                pass
        if sender is not None:
            sender.join(timeout=4)
        row.update(after_hook=state(), report=start.pid(), warnings=warnings,
                   interrupts_after_observation=interrupts,
                   errors=errors, reaps=reaps, application=application,
                   sent=len(sent), storm_span=(sent[-1]-sent[0] if sent else 0))
        row['post_interrupts'] = post_interrupts
        row['last_signal_after_unwind'] = sent[-1] - returned if sent else None
        row['attacks'] = attacks
        if mode.startswith('lease-'):
            row.update(absolute_deadline=start.deadline, start_done=start.done,
                       owner_token_alive=start.lease.owner() is not None,
                       consumer_alive=start.lease.service.thread.is_alive(),
                       registered=start in start.lease.service.starts,
                       run_lease=b.RUN_LEASE, retained_exception=type(kept_exception).__name__,
                       observed_after_entry=time.monotonic() - begin)
        if mode == 'warning-lock':
            row['main_owns_handler_lock'] = capture.lock._is_owned()
            row['consumer_alive'] = start.lease.service.thread.is_alive()
            frame = sys._current_frames().get(start.lease.service.thread.ident)
            stack = []
            while frame is not None:
                stack.append([frame.f_code.co_filename, frame.f_lineno, frame.f_code.co_name])
                frame = frame.f_back
            row['consumer_stack'] = stack
            emitter = getattr(start.lease.service, 'emitter', None)
            row['warning_pending'] = any(not r.warned and r.warning_pending
                                         for r in start.releases.values())
            row['emitter_outstanding'] = emitter is not None
            frame = sys._current_frames().get(emitter.ident) if emitter is not None else None
            stack = []
            while frame is not None:
                stack.append([frame.f_code.co_filename, frame.f_lineno, frame.f_code.co_name])
                frame = frame.f_back
            row['emitter_stack'] = stack
        return row
    finally:
        signal.signal(signal.SIGUSR1, previous_handler)
        b.WORKER_RECLAIM_TIMEOUT = original_timeout
        if original_lease is not None:
            b.RUN_LEASE = original_lease
        b._kill_and_reap = real_reap
        logger.removeHandler(capture)
        logger.setLevel(previous_level)
        if in_pytest:
            sys._eedf_application_before_fork = None
            sys._eedf_application_after_fork = None
        os.close(ready_r)
        os.close(ready_w)
        if mode != 'exit-immediately' and worker in (b._children() or set()):
            try:
                os.kill(worker, signal.SIGKILL)
                os.waitpid(worker, 0)
            except (ProcessLookupError, ChildProcessError):
                pass


if __name__ == '__main__':
    print(json.dumps(run(sys.argv[1], sys.argv[2])), flush=True)
