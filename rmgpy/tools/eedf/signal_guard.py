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
# to deal in the Software without restriction, including without limitation  #
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

"""Defer owned terminal signals across one qualification transaction.

The ownership machinery is adapted from
``8b37dbaf06e5e651723dd2ccae3f75e008eb7e16:rmgpy/tools/eedf/corroborate.py``.
It reads the raw kernel action before any install and identifies CPython's
trampoline by the exact local ``signal_handler`` symbol in the mapped ELF
image.  Foreign C handlers are never touched.
"""

import ctypes
import os
import platform
import signal
import struct
import sys
import threading


class _KernelSigaction(ctypes.Structure):
    _fields_ = [
        ("handler", ctypes.c_void_p),
        ("flags", ctypes.c_ulong),
        ("restorer", ctypes.c_void_p),
        ("mask", ctypes.c_ulong),
    ]


_SA_SIGINFO = 4
_PYTHON_TRAMPOLINE = None
_SYS_RT_SIGACTION = {"x86_64": 13, "aarch64": 134}.get(platform.machine())
try:
    _LIBC_SYSCALL = (
        ctypes.CDLL(None, use_errno=True).syscall
        if sys.platform.startswith("linux") and _SYS_RT_SIGACTION is not None
        else None
    )
except (OSError, AttributeError):
    _LIBC_SYSCALL = None


def _rt_sigaction(signum, new, old):
    return _LIBC_SYSCALL(
        ctypes.c_long(_SYS_RT_SIGACTION), ctypes.c_long(signum),
        ctypes.byref(new) if new is not None else None,
        ctypes.byref(old) if old is not None else None,
        ctypes.c_long(ctypes.sizeof(ctypes.c_ulong)),
    )


def _raw_sigaction(signum):
    """Return the kernel action, or None when raw actions are unavailable."""
    if _LIBC_SYSCALL is None:
        return None
    action = _KernelSigaction()
    if _rt_sigaction(signum, None, action) != 0:
        return None
    return action


def _restore_raw_sigaction(signum, action):
    _rt_sigaction(signum, action, None)


class _DlInfo(ctypes.Structure):
    _fields_ = [
        ("dli_fname", ctypes.c_char_p),
        ("dli_fbase", ctypes.c_void_p),
        ("dli_sname", ctypes.c_char_p),
        ("dli_saddr", ctypes.c_void_p),
    ]


try:
    _DLADDR = ctypes.CDLL(None).dladdr
except (OSError, AttributeError):
    try:
        _DLADDR = ctypes.CDLL("libdl.so.2").dladdr
    except (OSError, AttributeError):
        _DLADDR = None


def _image_base(address):
    """Return the load address of the object containing ``address``."""
    if not address or _DLADDR is None:
        return None
    info = _DlInfo()
    if not _DLADDR(ctypes.c_void_p(address), ctypes.byref(info)):
        return None
    return info.dli_fbase


try:
    _PYTHON_IMAGE = _image_base(
        ctypes.cast(ctypes.pythonapi.PyOS_getsig, ctypes.c_void_p).value)
except (AttributeError, ValueError):
    _PYTHON_IMAGE = None


def _mapped_file(address):
    """Return the inode-checked mapped file containing ``address``."""
    try:
        with open("/proc/self/maps") as maps:
            for line in maps:
                fields = line.split(None, 5)
                if len(fields) < 6:
                    continue
                start, end = (
                    int(bound, 16) for bound in fields[0].split("-"))
                if not start <= address < end:
                    continue
                path, inode = fields[5].rstrip("\n"), int(fields[4])
                if (not path.startswith("/") or path.endswith(" (deleted)")
                        or not inode):
                    return None
                return path if os.stat(path).st_ino == inode else None
    except (OSError, ValueError):
        return None
    return None


def _elf_function_offset(path, name):
    """Read one ELF64 function offset from ``.symtab``, failing closed."""
    wanted = name.encode() + b"\0"
    try:
        with open(path, "rb") as image:
            header = image.read(64)
            if (len(header) < 64 or header[:4] != b"\x7fELF"
                    or header[4] != 2 or header[5] not in (1, 2)):
                return None
            order = "<" if header[5] == 1 else ">"
            phoff, shoff = struct.unpack_from(order + "QQ", header, 0x20)
            phentsize, phnum, shentsize, shnum = struct.unpack_from(
                order + "HHHH", header, 0x36)
            if phentsize != 56 or shentsize != 64 or not phnum or not shnum:
                return None
            image.seek(phoff)
            segments = image.read(56 * phnum)
            image.seek(shoff)
            sections = [
                struct.unpack_from(order + "IIQQQQII", image.read(64))
                for _ in range(shnum)
            ]
            loads = [
                vaddr
                for kind, _, _, vaddr in struct.iter_unpack(
                    order + "IIQQ32x", segments)
                if kind == 1
            ]
            symtabs = [section for section in sections if section[1] == 2]
            if not loads or len(symtabs) != 1 or symtabs[0][6] >= shnum:
                return None
            _, _, _, _, offset, size, link, _ = symtabs[0]
            image.seek(sections[link][4])
            strings = image.read(sections[link][5])
            image.seek(offset)
            symbols = image.read(size - size % 24)
        page = os.sysconf("SC_PAGE_SIZE")
    except (OSError, TypeError, ValueError, struct.error):
        return None
    values = {
        value
        for name_offset, info, _, section, value, _ in struct.iter_unpack(
            order + "IBBHQQ", symbols)
        if (info & 0xF == 2 and section
            and strings[name_offset:name_offset + len(wanted)] == wanted)
    }
    if len(values) != 1:
        return None
    return values.pop() - (min(loads) & ~(page - 1))


_SIGNAL_HANDLER_UNREAD = object()
_SIGNAL_HANDLER = _SIGNAL_HANDLER_UNREAD


def _cpython_signal_handler():
    """Return the exact mapped CPython ``signal_handler`` address."""
    global _SIGNAL_HANDLER
    if _SIGNAL_HANDLER is _SIGNAL_HANDLER_UNREAD:
        _SIGNAL_HANDLER = None
        if _PYTHON_IMAGE is not None:
            address = ctypes.cast(
                ctypes.pythonapi.PyOS_getsig, ctypes.c_void_p).value
            path = _mapped_file(address)
            offset = (_elf_function_offset(path, "signal_handler")
                      if path else None)
            if offset is not None:
                _SIGNAL_HANDLER = _PYTHON_IMAGE + offset
    return _SIGNAL_HANDLER


def _cpython_handler_by_image(python_views, actions):
    """Identify CPython's trampoline without changing a disposition."""
    trampoline = _cpython_signal_handler()
    if trampoline is None:
        return None
    found = set()
    for signum, action in actions.items():
        if action is None or not callable(python_views.get(signum)):
            continue
        handler = action.handler or 0
        if handler == trampoline and _image_base(handler) == _PYTHON_IMAGE:
            found.add(handler)
    return found.pop() if len(found) == 1 else None


def _effective_disposition(python_view, action):
    """Decide ownership from the current raw action, never the Python cache."""
    handler = action.handler or 0
    if handler == int(signal.SIG_DFL):
        return signal.SIG_DFL
    if handler == int(signal.SIG_IGN):
        return signal.SIG_IGN
    if action.flags & _SA_SIGINFO or not callable(python_view):
        return None
    if _PYTHON_TRAMPOLINE is None or handler != _PYTHON_TRAMPOLINE:
        return None
    return python_view


def run_guarded_attempt(lock, operation, finalize, complete):
    """Run one terminal attempt under I-333 signal ownership and lock order.

    ``operation`` must leave acceptance unpublished. ``finalize`` receives
    cancellation or failure and returns its record. Only a clean signal
    release calls ``complete`` to publish success. The lock covers all three
    phases, including returning signal handlers and final publication.
    """
    global _PYTHON_TRAMPOLINE
    deferred_signals = []
    previous_signal_handlers = {}
    previous_raw_actions = {}
    previous_actions = {}

    def defer_signal(signum, frame):
        if signum not in deferred_signals:
            deferred_signals.append(signum)

    result = None
    failure = None
    finalization_error = None
    restored_signals = set()

    def finalize_attempt(error):
        nonlocal result, finalization_error
        try:
            result = finalize(error)
        except BaseException as caught:
            if finalization_error is None:
                finalization_error = caught

    def signal_error(signum):
        return RuntimeError(
            '{0} interrupted qualification'.format(
                signal.Signals(signum).name))

    try:
        try:
            if threading.current_thread() is threading.main_thread():
                # This is the cleared I-333 ownership order.  Read both raw
                # actions before any install, put defaults first so an install
                # can teach the trampoline, then identify the exact local ELF
                # signal_handler only when no default install can do so.
                raw_actions = {}
                for guarded_signal in (signal.SIGTERM, signal.SIGINT):
                    raw_actions[guarded_signal] = _raw_sigaction(
                        guarded_signal)
                order = [
                    signum for signum, raw in raw_actions.items()
                    if raw is not None and not raw.handler
                ]
                order += [
                    signum for signum, raw in raw_actions.items()
                    if raw is not None and raw.handler
                ]
                if (_PYTHON_TRAMPOLINE is None
                        and not any(raw is not None and not raw.handler
                                    for raw in raw_actions.values())):
                    _PYTHON_TRAMPOLINE = _cpython_handler_by_image(
                        {signum: signal.getsignal(signum)
                         for signum in raw_actions},
                        raw_actions)
                for guarded_signal in order:
                    previous_raw = _raw_sigaction(guarded_signal)
                    if previous_raw is None:
                        continue
                    previous_view = signal.getsignal(guarded_signal)
                    previous_action = _effective_disposition(
                        previous_view, previous_raw)
                    if (previous_action is None
                            or previous_action == signal.SIG_IGN):
                        continue
                    previous_raw_actions[guarded_signal] = previous_raw
                    previous_actions[guarded_signal] = previous_action
                    previous_signal_handlers[guarded_signal] = previous_action
                    signal.signal(guarded_signal, defer_signal)
                    if _PYTHON_TRAMPOLINE is None:
                        installed = _raw_sigaction(guarded_signal)
                        if installed is not None and installed.handler:
                            _PYTHON_TRAMPOLINE = installed.handler
        except KeyboardInterrupt as error:
            failure = error
        except BaseException as error:
            failure = error

        with lock:
            if failure is None and not deferred_signals:
                try:
                    result = operation()
                except KeyboardInterrupt as error:
                    failure = error
                except BaseException as error:
                    failure = error
            if failure is not None:
                finalize_attempt(failure)
            elif deferred_signals:
                # Cancellation may have skipped operation altogether, leaving
                # an earlier acceptance to revoke before either signal restores.
                finalize_attempt(signal_error(
                    signal.SIGTERM if signal.SIGTERM in deferred_signals
                    else signal.SIGINT))

            # No PASS exists during release, including when TERM has its
            # original action while INT still defers. Late cancellation can
            # therefore terminate here without leaving acceptance behind.
            # Retain I-333's TERM priority and resample signals received inside
            # signal.signal before deciding whether completion is allowed.
            term_outcome = False
            int_outcome = False
            try:
                if signal.SIGTERM in previous_signal_handlers:
                    term_outcome = signal.SIGTERM in deferred_signals
                    if term_outcome:
                        finalize_attempt(signal_error(signal.SIGTERM))
                    signal.signal(
                        signal.SIGTERM,
                        previous_signal_handlers[signal.SIGTERM])
                    _restore_raw_sigaction(
                        signal.SIGTERM,
                        previous_raw_actions[signal.SIGTERM])
                    restored_signals.add(signal.SIGTERM)
                    late_term = signal.SIGTERM in deferred_signals
                    if late_term and not term_outcome:
                        finalize_attempt(signal_error(signal.SIGTERM))
                    term_outcome = late_term
                if term_outcome:
                    if previous_actions[signal.SIGTERM] == signal.SIG_DFL:
                        os.kill(os.getpid(), signal.SIGTERM)
                    else:
                        signal.raise_signal(signal.SIGTERM)
            finally:
                if signal.SIGINT in previous_signal_handlers:
                    int_outcome = (not term_outcome
                                   and signal.SIGINT in deferred_signals)
                    if int_outcome:
                        finalize_attempt(signal_error(signal.SIGINT))
                    signal.signal(
                        signal.SIGINT,
                        previous_signal_handlers[signal.SIGINT])
                    _restore_raw_sigaction(
                        signal.SIGINT,
                        previous_raw_actions[signal.SIGINT])
                    restored_signals.add(signal.SIGINT)
            if not term_outcome:
                late_int = signal.SIGINT in deferred_signals
                if late_int and not int_outcome:
                    finalize_attempt(signal_error(signal.SIGINT))
                int_outcome = late_int
            if not term_outcome and int_outcome:
                if previous_actions[signal.SIGINT] == signal.SIG_DFL:
                    os.kill(os.getpid(), signal.SIGINT)
                else:
                    signal.raise_signal(signal.SIGINT)
            if finalization_error is not None:
                raise finalization_error
            if failure is None and not deferred_signals:
                try:
                    result = complete(result)
                except BaseException as error:
                    failure = error
                    finalize_attempt(error)
                if finalization_error is not None:
                    raise finalization_error
    finally:
        # A failed partial install or lock entry still restores every action
        # installed so far.  Normal release above records each restored signal
        # so this fallback performs no duplicate delivery.
        for guarded_signal in (signal.SIGTERM, signal.SIGINT):
            if (guarded_signal in previous_signal_handlers
                    and guarded_signal not in restored_signals):
                signal.signal(
                    guarded_signal,
                    previous_signal_handlers[guarded_signal])
                _restore_raw_sigaction(
                    guarded_signal,
                    previous_raw_actions[guarded_signal])

    if failure is not None:
        raise failure
    return result
