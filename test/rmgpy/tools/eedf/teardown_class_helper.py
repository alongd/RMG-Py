#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the "Software"),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""Real-process cleanup probes invoked by corroborateTest.py."""

import contextlib
import errno
import inspect
import json
import logging
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading
import time
from unittest import mock

from rmgpy.tools.eedf import corroborate as c


CASES = (
    "sigint-cgroup",
    "repeated-sigint-cgroup",
    "sigint-unit",
    "repeated-sigint-unit",
    "baseexception-cgroup",
    "baseexception-unit",
    "baseexception-verify",
    "sigint-live-launcher",
    "sigterm-cgroup",
    "sigterm-unit",
    "sigterm-verify",
    "unit-timeout",
    "unit-oserror",
    "verify-oserror",
    "sigint-inside-spawn",
    "sigterm-inside-spawn",
    "sigint-after-spawn",
    "idle-thread-sigint-cgroup",
    "idle-thread-sigterm-cgroup",
    "stat-permission",
    "stat-initial-eio",
    "stat-attempt-eio",
    "stat-recheck-eio",
    "mixed-replay",
    "mixed-custom-replay",
    "term-int-at-release-term",
    "term-int-at-release-deliver",
    "int-only-at-release-term",
    "int-only-at-release-int",
    "log-failure-control",
    "log-failure-removed-directory",
    "log-failure-permission-denied",
)
REAL_POPEN = subprocess.Popen
REAL_WRITE = Path.write_text
REAL_EXISTS = Path.exists


def worker(root, mode):
    case_root = root / mode
    case_root.mkdir(parents=True, exist_ok=True)
    ready = case_root / "owned-child.pid"
    marker = case_root / "survived"
    helper = case_root / "helper.py"
    helper.write_text(
        "import os,time,pathlib\n"
        "if os.fork()==0:\n"
        " os.setsid()\n"
        " if os.fork()==0:\n"
        "  f=os.open('/dev/null',os.O_RDWR);os.dup2(f,1);os.dup2(f,2)\n"
        f"  pathlib.Path({str(ready)!r}).write_text(str(os.getpid()))\n"
        "  time.sleep(.8)\n"
        f"  pathlib.Path({str(marker)!r}).write_text('survived')\n"
        "  os._exit(0)\n"
        " os._exit(0)\n"
        + ("time.sleep(5)\n" if mode == "sigint-live-launcher" else "")
    )
    group = case_root / "sealed" / "group"
    group.mkdir(parents=True, exist_ok=True)
    (group / "cgroup.events").write_text("populated 1\n")
    backend = c.SystemdScopeBackend()
    holder = {}
    counts = {
        "cgroup_calls": 0,
        "unit_calls": 0,
        "verify_calls": 0,
        "verified": False,
        "stat_calls": 0,
        "signals": [],
        "replayed": [],
    }
    original_empty = backend._empty

    def save():
        REAL_WRITE(case_root / "counts.json", json.dumps(counts))

    def send(sig):
        counts["signals"].append(signal.Signals(sig).name)
        save()
        if mode.startswith("idle-thread-"):
            # Thread-directed at an unmasked idle thread: the OS runs the C
            # handler there, and Python must still only record it.
            signal.pthread_kill(idle.ident, sig)
        else:
            os.kill(os.getpid(), sig)

    def terminate():
        try:
            os.kill(int(ready.read_text()), signal.SIGKILL)
        except ProcessLookupError:
            pass
        if holder["process"].poll() is None:
            holder["process"].kill()
        if mode != "stat-permission":
            REAL_WRITE(group / "cgroup.events", "populated 0\n")

    def launch(argv, **kwargs):
        holder["process"] = REAL_POPEN([sys.executable, str(helper)], **kwargs)
        if mode in ("sigint-inside-spawn", "sigterm-inside-spawn"):
            # The child exists but Popen has not returned to run() yet.
            deadline = time.monotonic() + 2.0
            while not ready.exists() and time.monotonic() < deadline:
                time.sleep(0.002)
            assert ready.exists()
            send(signal.SIGTERM if mode.startswith("sigterm-") else signal.SIGINT)
        return holder["process"]

    def control(unit):
        deadline = time.monotonic() + 2.0
        while not ready.exists() and time.monotonic() < deadline:
            time.sleep(0.002)
        assert ready.exists()
        if mode != "sigint-live-launcher":
            holder["process"].wait(timeout=2.0)
        if mode == "stat-permission":
            os.chmod(group.parent, 0)
        return group

    def write(path, *args, **kwargs):
        if path == group / "cgroup.kill":
            counts["cgroup_calls"] += 1
            save()
            if "unit" in mode or mode in ("stat-recheck-eio", "mixed-replay", "mixed-custom-replay"):
                raise PermissionError(errno.EACCES, "cgroup.kill unavailable")
            if mode == "baseexception-cgroup" and counts["cgroup_calls"] == 1:
                raise KeyboardInterrupt("one-shot kill interruption")
            if mode in ("sigint-cgroup", "sigint-live-launcher", "repeated-sigint-cgroup", "idle-thread-sigint-cgroup"):
                send(signal.SIGINT)
            if mode == "repeated-sigint-cgroup":
                send(signal.SIGINT)
            if mode in ("sigterm-cgroup", "idle-thread-sigterm-cgroup") or mode.startswith("term-int-"):
                send(signal.SIGTERM)
            if mode.startswith("log-failure-"):
                send(signal.SIGINT)
                send(signal.SIGTERM)
            terminate()
        return REAL_WRITE(path, *args, **kwargs)

    def unit_command(argv, **kwargs):
        assert "--kill-whom=all" in argv and "--signal=KILL" in argv
        counts["unit_calls"] += 1
        save()
        if mode in ("mixed-replay", "mixed-custom-replay"):
            if counts["unit_calls"] == 1:
                raise subprocess.TimeoutExpired(argv, 2.0)
            send(signal.SIGINT)
            send(signal.SIGTERM)
        if mode == "unit-timeout" and counts["unit_calls"] == 1:
            raise subprocess.TimeoutExpired(argv, 2.0)
        if mode == "unit-oserror" and counts["unit_calls"] == 1:
            raise OSError(errno.EIO, "one-shot user-bus I/O error")
        if mode == "baseexception-unit" and counts["unit_calls"] == 1:
            raise KeyboardInterrupt("one-shot unit interruption")
        if mode in ("sigint-unit", "repeated-sigint-unit"):
            send(signal.SIGINT)
        if mode == "repeated-sigint-unit":
            send(signal.SIGINT)
        if mode == "sigterm-unit":
            send(signal.SIGTERM)
        terminate()
        return subprocess.CompletedProcess(argv, 0, "", "")

    def empty(known_group):
        counts["verify_calls"] += 1
        save()
        if mode == "baseexception-verify" and counts["verify_calls"] == 1:
            raise SystemExit("one-shot verification interruption")
        if mode == "verify-oserror" and counts["verify_calls"] == 1:
            raise OSError(errno.EIO, "one-shot verification I/O error")
        if mode == "sigterm-verify":
            send(signal.SIGTERM)
        result = original_empty(known_group)
        counts["verified"] = result
        save()
        return result

    def exists(path):
        if path == group:
            counts["stat_calls"] += 1
            save()
            injected_at = {
                "stat-initial-eio": 1,
                "stat-attempt-eio": 2,
                "stat-recheck-eio": 3,
            }.get(mode)
            if counts["stat_calls"] == injected_at:
                raise OSError(errno.EIO, "one-shot " + mode)
        return REAL_EXISTS(path)

    def replay_handler(signum, frame):
        counts["replayed"].append(signal.Signals(signum).name)
        save()
        if signum == signal.SIGTERM:
            raise SystemExit("custom SIGTERM replay")
        raise KeyboardInterrupt("custom SIGINT replay")

    source, first = inspect.getsourcelines(c.SystemdScopeBackend.run)
    # The first statement after Popen has returned the launcher to run().
    target = next(first + i for i, line in enumerate(source) if "if process is not None:" in line)

    release_lines = {
        "term": first + next(
            i
            for i, line in enumerate(source)
            if "signal.signal(signal.SIGTERM, previous_signal_handlers[signal.SIGTERM])" in line
        ),
        "deliver": first + next(
            i
            for i, line in enumerate(source)
            if "if previous_actions[signal.SIGTERM] == signal.SIG_DFL:" in line
        ),
        "int": first + next(
            i
            for i, line in enumerate(source)
            if "signal.signal(signal.SIGINT, previous_signal_handlers[signal.SIGINT])" in line
        ),
    }
    release_fired = []

    def release_trace(frame, event, arg):
        # A real SIGINT lands as the release reaches the SIGTERM restoration,
        # the SIGTERM re-delivery, or the SIGINT restoration (deterministic hook).
        if frame.f_code is c.SystemdScopeBackend.run.__code__ and event == "line":
            where = mode.rsplit("-", 1)[1]
            if frame.f_lineno == release_lines[where] and not release_fired:
                release_fired.append(where)
                send(signal.SIGINT)
        return release_trace

    def trace(frame, event, arg):
        if (
            mode == "sigint-after-spawn"
            and frame.f_code is c.SystemdScopeBackend.run.__code__
            and event == "line"
            and frame.f_lineno == target
        ):
            sys.settrace(None)
            deadline = time.monotonic() + 2.0
            while not ready.exists() and time.monotonic() < deadline:
                time.sleep(0.002)
            assert ready.exists()
            send(signal.SIGINT)
        return trace

    exists_patch = (
        mock.patch.object(Path, "exists", exists)
        if mode.startswith("stat-") and mode != "stat-permission"
        else contextlib.nullcontext()
    )
    log_handler = None
    log_dir = case_root / "configured-logs"
    if mode.startswith("log-failure-"):
        # Stdlib FileHandler(delay=True) opens its stream inside emit(), outside
        # the try that routes errors to handleError, so a missing or unreadable
        # directory raises out of logger.info().
        log_dir.mkdir(exist_ok=True)
        log_handler = logging.FileHandler(log_dir / "solver.log", delay=True)
        c.logger.addHandler(log_handler)
        c.logger.setLevel(logging.INFO)
        if mode == "log-failure-removed-directory":
            log_dir.rmdir()
        if mode == "log-failure-permission-denied":
            log_dir.chmod(0)
    idle_stop = threading.Event()
    idle = threading.Thread(target=idle_stop.wait, daemon=True)
    if mode.startswith("idle-thread-"):
        idle.start()
    if mode == "mixed-custom-replay" or mode.startswith("log-failure-"):
        signal.signal(signal.SIGINT, replay_handler)
        signal.signal(signal.SIGTERM, replay_handler)
    with mock.patch.object(c.subprocess, "Popen", side_effect=launch), mock.patch.object(
        backend, "_control_group", side_effect=control
    ), mock.patch.object(c.subprocess, "run", side_effect=unit_command), mock.patch.object(
        Path, "write_text", write
    ), mock.patch.object(backend, "_empty", side_effect=empty), exists_patch:
        try:
            if mode == "sigint-after-spawn":
                sys.settrace(trace)
            if "-at-release-" in mode:
                sys.settrace(release_trace)
            result = backend.run(["finite-launcher"], cwd=case_root, timeout_s=0.1)
            outcome = "return cleanup_passed=" + str(result["cleanup_passed"])
        except BaseException as exc:
            outcome = type(exc).__name__ + ": " + str(exc)
        finally:
            sys.settrace(None)
            idle_stop.set()
            os.chmod(group.parent, 0o700)
            if log_handler is not None:
                c.logger.removeHandler(log_handler)
                if mode == "log-failure-permission-denied":
                    log_dir.chmod(0o700)
                counts["pending_after"] = sorted(
                    signal.Signals(x).name for x in signal.sigpending()
                )
                counts["blocked_after"] = sorted(
                    signal.Signals(x).name
                    for x in signal.pthread_sigmask(signal.SIG_BLOCK, set())
                )
                counts["sigint_handler_restored"] = (
                    signal.getsignal(signal.SIGINT) is replay_handler
                )
    holder["process"].wait(timeout=2.0)
    save()
    REAL_WRITE(case_root / "outcome.txt", outcome)


def main(root, cases):
    root.mkdir(parents=True, exist_ok=True)
    outside = REAL_POPEN(["/bin/sleep", "40"])
    try:
        for mode in cases:
            case_root = root / mode
            case_root.mkdir(exist_ok=True)
            result = subprocess.run(
                [sys.executable, __file__, str(root), "--case", mode],
                capture_output=True,
                text=True,
                timeout=7.0,
            )
            REAL_WRITE(case_root / "reviewer.stdout.log", result.stdout)
            REAL_WRITE(case_root / "reviewer.stderr.log", result.stderr)
            time.sleep(1.0)
            counts = json.loads((case_root / "counts.json").read_text())
            outcome_path = case_root / "outcome.txt"
            outcome = outcome_path.read_text() if outcome_path.exists() else "reviewer killed by signal"
            print(
                mode,
                "exit",
                result.returncode,
                "outcome",
                outcome,
                "finite_descendant_survived",
                (case_root / "survived").exists(),
                "unrelated_alive",
                outside.poll() is None,
                json.dumps(counts, sort_keys=True),
                flush=True,
            )
            assert outside.poll() is None
    finally:
        outside.terminate()
        outside.wait(timeout=2.0)


if __name__ == "__main__":
    probe_root = Path(sys.argv[1])
    if len(sys.argv) == 4 and sys.argv[2] == "--case":
        worker(probe_root, sys.argv[3])
    else:
        selected = tuple(sys.argv[2].split(",")) if len(sys.argv) == 3 else CASES
        main(probe_root, selected)
