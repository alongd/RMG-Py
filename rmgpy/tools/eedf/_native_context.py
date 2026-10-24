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

"""Private fork context for scanner workers; application contexts are unchanged.

Copyright (c) 2026 the RMG Team. Distributed under the MIT license in LICENSE.txt.
"""

import math
import multiprocessing
import os
import platform
import sys
import time
from multiprocessing import util
from rmgpy.data.kinetics.family import install_complete_reducers

install_complete_reducers()

try:
    from multiprocessing.context import ForkContext, ForkProcess
    from multiprocessing.popen_fork import Popen
except ImportError:  # Non-fork platforms retain the ordinary in-process path.
    ForkContext = ForkProcess = Popen = object

try:
    from . import _native_worker as native
except ImportError:
    native = None


class NativeWorkerUnsupported(RuntimeError):
    """The requested worker mode lacks the qualified native lifetime mechanism."""


class NativeWorkerStartupError(RuntimeError):
    """A native worker could not complete its watchdog admission handshake."""


class ScanLeaseExpired(RuntimeError):
    """The calling scan thread lost its ownership lease; no result is accepted."""

    def __init__(self, facts):
        self.facts = facts
        super().__init__('owner-unresponsive: scan ownership lease expired; native facts: {!r}'.format(facts))

    def __reduce__(self):
        return type(self), (self.facts,)


class NativeWorkerStatusError(RuntimeError):
    """An exact successful exit cannot be proved for a received worker result."""


def require_supported():
    """Refuse before allocating scan worker resources; never choose a fallback."""
    if native is None:
        raise NativeWorkerUnsupported('native worker extension is unavailable; run make build')
    facts = native.check_platform()
    if (sys.platform != 'linux' or platform.machine() != 'x86_64'
            or sys.implementation.name != 'cpython' or sys.version_info[:2] != (3, 9)
            or not facts['supported'] or facts['record_abi'] != 1
            or multiprocessing.get_context().get_start_method() != 'fork'):
        raise NativeWorkerUnsupported('requires Linux x86_64, glibc 2.39, CPython 3.9, '
                                      'lock-free ABI 1 and fork context; observed {!r}'.format(facts))
    return facts


class Owner:
    """Calling-thread facade; the C ledger outlives interrupted Python stores."""

    def __init__(self, lease_s, *, trace=False):
        if not math.isfinite(lease_s) or lease_s <= 0:
            raise ValueError('lease_s must be finite and positive')
        self.platform = require_supported()
        self.handle = native.owner(lease_s, trace)

    def checkpoint(self):
        if not native.renew(self.handle):
            facts = self.facts()
            if facts['reason'] == 7:
                raise NativeWorkerStartupError('watchdog admission failed: {!r}'.format(facts))
            if facts['reason'] == 4:
                raise ScanLeaseExpired(facts)
            raise NativeWorkerStatusError('native scan cancelled: {!r}'.format(facts))

    def process(self, target, args):
        process = NativeForkProcess(target=target, args=args)
        process._native_record = native.record(self.handle)
        return process

    def cleanup(self, seconds=5.0):
        # This one C operation publishes every cancellation before any Python
        # handler/logging callback can run, and caches every exact reap first.
        return native.cleanup(self.handle, seconds)

    def facts(self):
        return dict(native.owner_info(self.handle), platform=self.platform)


class NativePopen(Popen):
    """Standard fork pipe/bootstrap sequence with native identity and exit."""

    def _launch(self, process_obj):
        self.record = process_obj._native_record
        code = 1
        parent_r, child_w = os.pipe()
        child_r, parent_w = os.pipe()
        try:
            self.pid = native.fork(self.record)
        except BaseException:
            # The ledger, not self.pid, owns a child whose return was lost.
            for fd in (parent_r, child_w, child_r, parent_w):
                os.close(fd)
            facts = native.info(self.record)
            if facts['reason'] == 7:
                raise NativeWorkerStartupError('native watchdog admission failed: {!r}; platform: {!r}'.format(
                    facts, native.check_platform()))
            if facts['reason'] == 4:
                raise ScanLeaseExpired(facts)
            raise
        if self.pid == 0:
            try:
                os.close(parent_r)
                os.close(parent_w)
                code = process_obj._bootstrap(parent_sentinel=child_r)
            finally:
                native.finish(self.record, code)
        else:
            os.close(child_w)
            os.close(child_r)
            self.finalizer = util.Finalize(self, util.close_fds, (parent_r, parent_w))
            self.sentinel = parent_r

    def poll(self, flag=os.WNOHANG):
        if self.returncode is None:
            while True:
                reaped = native.poll(self.record)
                if reaped:
                    if reaped == 2:
                        # Unknown status must not become a successful outcome.
                        self.returncode = 255
                    else:
                        self.returncode = os.waitstatus_to_exitcode(native.info(self.record)['wait_status'])
                    break
                if flag == os.WNOHANG:
                    break
                time.sleep(0.001)
        return self.returncode

    def wait(self, timeout=None):
        # An adapter's descendants can hold a pipe sentinel after the exact
        # worker exits. Only waitpid on the recorded identity determines exit.
        until = None if timeout is None else time.monotonic() + timeout
        while self.poll() is None:
            if until is not None and time.monotonic() >= until:
                break
            time.sleep(0.001)
        return self.returncode

    def terminate(self):
        native.cancel(self.record)

    def kill(self):
        native.cancel(self.record)

    def _send_signal(self, sig):
        native.cancel(self.record)


class NativeForkProcess(ForkProcess):
    @staticmethod
    def _Popen(process_obj):
        return NativePopen(process_obj)


class NativeForkContext(ForkContext):
    Process = NativeForkProcess


def accepted(process, envelope):
    """Pipe readability never establishes successful completion by itself."""
    facts = native.info(process._native_record)
    if facts['reason'] == 4:
        raise ScanLeaseExpired(facts)
    if (facts['phase'] != 'normal' or facts['reaped'] != 1
            or not os.WIFEXITED(facts['wait_status']) or os.WEXITSTATUS(facts['wait_status']) != 0
            or not isinstance(envelope, tuple) or len(envelope) != 2
            or envelope[0] not in ('ok', 'error', 'lifecycle')):
        raise NativeWorkerStatusError('unverifiable worker result: {!r}'.format(facts))
    if envelope[0] == 'lifecycle':
        if isinstance(envelope[1], (ScanLeaseExpired, NativeWorkerUnsupported,
                                    NativeWorkerStartupError, NativeWorkerStatusError)):
            raise envelope[1]
        raise NativeWorkerStatusError('invalid lifecycle error envelope')
    return envelope
