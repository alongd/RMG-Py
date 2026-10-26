"""Regression tests for RDKit's process-level SIGINT disposition changes."""

import os
import platform
import signal
import subprocess
import sys

import pytest


_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))

_CHILD = r'''
import ctypes
import json
import os
import signal
import sys
import tempfile
import threading

from rdkit import Chem
from rdkit.Chem import AllChem
from rmgpy.qm.main import QMSettings
from rmgpy.qm.molecule import Geometry


class KernelSigaction(ctypes.Structure):
    _fields_ = [
        ("handler", ctypes.c_void_p),
        ("flags", ctypes.c_ulong),
        ("restorer", ctypes.c_void_p),
        ("mask", ctypes.c_ulong),
    ]


syscall_number = {"x86_64": 13, "aarch64": 134}[__import__("platform").machine()]
libc_syscall = ctypes.CDLL(None, use_errno=True).syscall
libc_syscall.restype = ctypes.c_long


def raw_sigaction():
    action = KernelSigaction()
    result = libc_syscall(
        ctypes.c_long(syscall_number),
        ctypes.c_long(signal.SIGINT),
        None,
        ctypes.byref(action),
        ctypes.c_long(ctypes.sizeof(ctypes.c_ulong)),
    )
    assert result == 0, ctypes.get_errno()
    return (action.handler, int(action.flags), action.restorer, int(action.mask))


def set_raw_sigaction_mask(mask):
    action = KernelSigaction()
    result = libc_syscall(
        ctypes.c_long(syscall_number),
        ctypes.c_long(signal.SIGINT),
        None,
        ctypes.byref(action),
        ctypes.c_long(ctypes.sizeof(ctypes.c_ulong)),
    )
    assert result == 0, ctypes.get_errno()
    action.mask = mask
    result = libc_syscall(
        ctypes.c_long(syscall_number),
        ctypes.c_long(signal.SIGINT),
        ctypes.byref(action),
        None,
        ctypes.c_long(ctypes.sizeof(ctypes.c_ulong)),
    )
    assert result == 0, ctypes.get_errno()


mode = sys.argv[1]
handled = []


def caller_handler(signum, frame):
    handled.append(signum)


if mode == "python_default" or mode in ("exception", "delayed"):
    signal.signal(signal.SIGINT, signal.default_int_handler)
elif mode == "nondefault":
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.siginterrupt(signal.SIGINT, False)
    set_raw_sigaction_mask(1 << (signal.SIGUSR1 - 1))
elif mode == "python_handler":
    signal.signal(signal.SIGINT, caller_handler)
elif mode == "ignore":
    signal.signal(signal.SIGINT, signal.SIG_IGN)
elif mode == "default":
    signal.signal(signal.SIGINT, signal.SIG_DFL)
else:
    raise AssertionError(mode)

with tempfile.TemporaryDirectory() as scratch:
    geometry = Geometry(QMSettings(scratchDirectory=scratch), "sigint-child", None)
    molecule = Chem.AddHs(Chem.MolFromSmiles("CCO"))

    original_embed = AllChem.EmbedMultipleConfs
    optimize_called = []
    original_optimize = AllChem.UFFOptimizeMolecule
    if mode == "exception":
        def embed_then_raise(*args, **kwargs):
            original_embed(*args, **kwargs)
            raise RuntimeError("injected after RDKit embedding")
        AllChem.EmbedMultipleConfs = embed_then_raise
    elif mode == "delayed":
        def embed_then_queue_sigint(*args, **kwargs):
            result = original_embed(*args, **kwargs)
            signal.pthread_kill(threading.get_ident(), signal.SIGINT)
            print("signal queued inside embed call", flush=True)
            return result
        def track_optimize(*args, **kwargs):
            optimize_called.append(True)
            return original_optimize(*args, **kwargs)
        AllChem.EmbedMultipleConfs = embed_then_queue_sigint
        AllChem.UFFOptimizeMolecule = track_optimize

    before = raw_sigaction()
    try:
        try:
            geometry.rd_embed(molecule, 1)
        except KeyboardInterrupt:
            if mode != "delayed":
                raise
            assert not optimize_called, "SIGINT was not delivered at embed-scope exit"
            after = raw_sigaction()
            assert after == before, (before, after)
            print("raw sigaction before/after:", json.dumps([before, after]), flush=True)
            print("delayed SIGINT delivered at embed-scope exit", flush=True)
        except RuntimeError as error:
            if mode != "exception" or str(error) != "injected after RDKit embedding":
                raise
            after = raw_sigaction()
            assert after == before, (before, after)
            print("raw sigaction before/after:", json.dumps([before, after]), flush=True)
            print("embedding exception restored raw sigaction", flush=True)
        else:
            if mode in ("exception", "delayed"):
                raise AssertionError("expected injected exception or delayed SIGINT")
            after = raw_sigaction()
            assert after == before, (before, after)
            print("raw sigaction before/after:", json.dumps([before, after]), flush=True)
            if mode == "nondefault":
                assert before[0] == 1, before
                assert before[1] & 0x10000000, before
                assert before[3] != 0, before
                print("non-default handler, SA_RESTART, and non-empty mask restored", flush=True)
    finally:
        AllChem.EmbedMultipleConfs = original_embed
        AllChem.UFFOptimizeMolecule = original_optimize

    if mode in ("python_default", "exception"):
        interrupts = 0
        for _ in range(2):
            try:
                os.kill(os.getpid(), signal.SIGINT)
            except KeyboardInterrupt:
                interrupts += 1
        assert interrupts == 2, interrupts
        print("two KeyboardInterrupts after embedding", flush=True)
    elif mode == "python_handler":
        for _ in range(2):
            os.kill(os.getpid(), signal.SIGINT)
        assert handled == [signal.SIGINT, signal.SIGINT], handled
        print("caller handler received both SIGINTs", flush=True)
    elif mode == "ignore":
        for _ in range(2):
            os.kill(os.getpid(), signal.SIGINT)
        print("ignored disposition survived both SIGINTs", flush=True)
    elif mode == "default":
        os.kill(os.getpid(), signal.SIGINT)
        raise AssertionError("SIG_DFL did not terminate the child")
'''


def _run_child(mode):
    environment = os.environ.copy()
    environment["PYTHONPATH"] = _REPO + os.pathsep + environment.get("PYTHONPATH", "")
    return subprocess.run(
        [sys.executable, "-c", _CHILD, mode],
        cwd=os.getcwd(),
        env=environment,
        capture_output=True,
        text=True,
        timeout=60,
    )


@pytest.fixture(autouse=True)
def _linux_sigaction_platform():
    if not sys.platform.startswith("linux") or platform.machine() not in ("x86_64", "aarch64"):
        pytest.skip("raw rt_sigaction evidence is only required on Linux x86_64/aarch64")


@pytest.mark.parametrize("mode", ["python_default", "python_handler", "ignore"])
def test_rd_embed_restores_sigint_disposition(mode):
    result = _run_child(mode)
    assert result.returncode == 0, result.stdout + result.stderr
    print(result.stdout, end="")
    assert "raw sigaction before/after:" in result.stdout
    if mode == "python_default":
        assert "two KeyboardInterrupts after embedding" in result.stdout
    elif mode == "python_handler":
        assert "caller handler received both SIGINTs" in result.stdout
    else:
        assert "ignored disposition survived both SIGINTs" in result.stdout


def test_rd_embed_restores_sig_dfl():
    result = _run_child("default")
    assert result.returncode == -signal.SIGINT, result.stdout + result.stderr
    print(result.stdout, end="")
    assert "raw sigaction before/after:" in result.stdout


def test_rd_embed_restores_nondefault_sigaction_flags_and_mask():
    result = _run_child("nondefault")
    assert result.returncode == 0, result.stdout + result.stderr
    print(result.stdout, end="")
    assert "non-default handler, SA_RESTART, and non-empty mask restored" in result.stdout


def test_rd_embed_restores_sigint_when_embedding_raises():
    result = _run_child("exception")
    assert result.returncode == 0, result.stdout + result.stderr
    print(result.stdout, end="")
    assert "embedding exception restored raw sigaction" in result.stdout
    assert "two KeyboardInterrupts after embedding" in result.stdout


def test_sigint_during_embed_is_delivered_at_scope_exit():
    result = _run_child("delayed")
    assert result.returncode == 0, result.stdout + result.stderr
    print(result.stdout, end="")
    assert "signal queued inside embed call" in result.stdout
    assert "delayed SIGINT delivered at embed-scope exit" in result.stdout
