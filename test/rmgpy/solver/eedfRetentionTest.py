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

"""Kernel resource regressions through the public qualification lifecycle."""

import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.parametrize('mode', ['default', 'no-gc', 'limit108', 'fixed', 'exports'])
def test_completed_ownership_has_constant_retention(tmp_path, mode, record_property):
    result = run_worker(tmp_path, mode)
    record_property('descriptor_counts', json.dumps(result['counts']))
    assert result['completed'] == 50, result
    assert len(set(result['counts'])) == 1, result
    assert result['registry'] == [1] * 50, result
    assert result['after_failure'] == result['before'], result
    assert result['retired_descriptors'] == 0, result


@pytest.mark.parametrize('route', ['supersede', 'failure', 'adoption', 'export', 'capture', 'settle'])
@pytest.mark.parametrize('signum', [2, 15])
@pytest.mark.parametrize('policy', ['returning', 'raising'])
def test_real_signal_at_retirement_releases_retained_objects(
        tmp_path, route, signum, policy, record_property):
    result = run_worker(tmp_path, 'signal', route, str(signum), policy)
    record_property('signal_events', json.dumps(result['events']))
    assert result['events'] == ['sent', 'received'], result
    assert result['after_failure'] == result['before'], result
    assert result['registry_after_failure'] == 0, result
    assert result['retired_descriptors'] == 0, result


@pytest.fixture(scope='module')
def close_interposer(tmp_path_factory):
    """Intercept the built helper's libc close; no product test hook is needed."""
    directory = tmp_path_factory.mktemp('native-close')
    source = directory / 'close.c'
    source.write_text(r'''
#define _GNU_SOURCE
#include <dlfcn.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <sys/syscall.h>
static int target = -1, position, after_close, notify_fd, ack_fd;
void retention_arm(int index, int after, int notify, int ack) {
    position = 0; target = index; after_close = after;
    notify_fd = notify; ack_fd = ack;
}
static void checkpoint(void) {
    char byte = 'x'; ssize_t result;
    target = -1;
    if (write(notify_fd, &byte, 1) != 1) _exit(91);
    do { result = read(ack_fd, &byte, 1); } while (result < 0 && errno == EINTR);
    if (result != 1) _exit(92);
}
int close(int fd) {
    Dl_info caller;
    int hit = target >= 0 && dladdr(__builtin_return_address(0), &caller)
        && strstr(caller.dli_fname, "_eedf_ownership") && position++ == target;
    if (hit && !after_close) checkpoint();
    int result = syscall(SYS_close, fd), saved_errno = errno;
    if (hit && after_close) checkpoint();
    errno = saved_errno;
    return result;
}
''')
    obj = directory / 'close.o'
    library = directory / 'close.so'
    # Separate compilation/linking keeps pytest plus compiler helpers <= 3.
    subprocess.run(['cc', '-fPIC', '-c', str(source), '-o', str(obj)], check=True)
    subprocess.run(['ld', '-shared', str(obj), '-o', str(library), '-ldl'], check=True)
    return library


@pytest.mark.parametrize('slot', range(5))
@pytest.mark.parametrize('after', [0, 1])
@pytest.mark.parametrize('signum', [2, 15])
@pytest.mark.parametrize('policy', ['returning', 'raising'])
def test_signal_inside_native_release_is_delivered_after_consistent_transition(
        tmp_path, close_interposer, slot, after, signum, policy, record_property):
    result = run_worker(tmp_path, 'native', str(slot), str(after), str(signum), policy,
                        preload=close_interposer)
    record_property('native_signal', json.dumps(result))
    assert result['sender_exit'] == 0 and result['deliveries'] == 1, result
    assert result['sentinel_survived'], result
    assert result['fd_delta'] == 0 and result['registry'] == 0, result


def run_worker(tmp_path, *arguments, preload=None):
    environment = dict(os.environ)
    environment.pop('PYTHONPATH', None)
    if preload is not None:
        environment['LD_PRELOAD'] = str(preload)
    completed = subprocess.run(
        [sys.executable, str(Path(__file__).resolve()), str(tmp_path), *arguments],
        cwd=Path(__file__).resolve().parents[3], env=environment,
        capture_output=True, text=True, timeout=180,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    return json.loads((tmp_path / 'retention.json').read_text())


def worker(root, mode, *arguments):
    import gc
    import resource
    import signal

    from eedfProviderTest import public_qualification_fixture, intended_signal_disposition
    from rmgpy.solver import eedf_provider as module

    if mode == 'native':
        return native_worker(root, *arguments)

    patches = pytest.MonkeyPatch()
    provider, state = public_qualification_fixture(root / 'fixture', patches)
    if mode != 'fixed':
        state.pop('record_path')
    owners, captures = [], []
    original_owner = module._QualificationExport.__init__
    original_capture = module._QualificationCapture.__init__

    def observe_owner(owner, *args, **kwargs):
        original_owner(owner, *args, **kwargs)
        owners.append(owner)

    def observe_capture(capture, *args, **kwargs):
        original_capture(capture, *args, **kwargs)
        captures.append(capture)

    patches.setattr(module._QualificationExport, '__init__', observe_owner)
    patches.setattr(module._QualificationCapture, '__init__', observe_capture)
    before = len(os.listdir('/proc/self/fd'))
    counts, registry, events = [], [], []
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    error = None
    completed = 0
    try:
        if mode == 'no-gc':
            gc.disable()
        if mode == 'limit108':
            resource.setrlimit(resource.RLIMIT_NOFILE, (108, hard))
        if mode != 'signal':
            if mode == 'exports':
                provider.qualify(state, 0.)
            for attempt in range(50):
                try:
                    if mode == 'exports':
                        provider.require_qualified('retention').commit_bytes(
                            root / ('export-' + str(attempt)), b'qualified output')
                    else:
                        assert provider.qualify(state, float(attempt))['status'] == 'PASS'
                except BaseException as caught:
                    error = repr(caught)
                    break
                completed += 1
                if mode != 'no-gc':
                    gc.collect()
                counts.append(len(os.listdir('/proc/self/fd')))
                registry.append(len(provider._qualification_files))
        else:
            route, number, policy = arguments
            signum = int(number)
            provider.qualify(state, 0.)
            if route == 'adoption':
                provider.record_terminal_failure(state, 1., RuntimeError('first failure'))
            targets = {
                'supersede': module._QualificationExport.revoke,
                'failure': module.EEDFProvider._write_failure,
                'adoption': module.EEDFProvider._invalidate_record,
                'export': module._QualificationAdmission.commit_bytes,
                'capture': module._QualificationCapture.capture,
                'settle': module._QualificationCapture.settle,
            }
            target = targets[route].__code__

            def handler(received, frame):
                events.append('received')
                # Every already released bundle must be completely inert.
                for owner in owners:
                    resources = getattr(owner, '_resources', None)
                    if resources is not None and resources.released:
                        assert owner.stream is None
                        assert owner.parent_fd is None
                        assert owner.private_parent_fd is None
                        assert owner.directory_fd is None
                if policy == 'raising':
                    raise KeyboardInterrupt('retirement boundary')

            def trace(frame, event, arg):
                if event == 'return' and frame.f_code is target and not events:
                    events.append('sent')
                    os.kill(os.getpid(), signum)
                return trace

            with intended_signal_disposition(signum, handler):
                try:
                    sys.settrace(trace)
                    if route == 'export':
                        provider.require_qualified('retention').commit_bytes(root / 'export', b'output')
                    elif route in ('failure', 'adoption'):
                        provider.record_terminal_failure(state, 2., RuntimeError('failure'))
                    else:
                        provider.qualify(state, 2.)
                except BaseException as caught:
                    error = repr(caught)
                finally:
                    sys.settrace(None)
    finally:
        resource.setrlimit(resource.RLIMIT_NOFILE, (soft, hard))
        gc.enable()
    # All observed objects stay alive: collecting them cannot satisfy this check.
    provider.record_terminal_failure(state, 99., RuntimeError('retention test complete'))
    after_failure = len(os.listdir('/proc/self/fd'))
    descriptors = sum(
        getattr(owner, field) is not None
        for owner in owners for field in ('parent_fd', 'private_parent_fd', 'directory_fd'))
    descriptors += sum(owner.stream is not None and not owner.stream.closed for owner in owners)
    descriptors += sum(capture.fd is not None for capture in captures)
    result = dict(before=before, counts=counts, registry=registry, completed=completed,
                  error=error, events=events, after_failure=after_failure,
                  registry_after_failure=len(provider._qualification_files),
                  retired_descriptors=descriptors)
    (root / 'retention.json').write_text(json.dumps(result, indent=2))
    patches.undo()


def native_worker(root, slot, after, number, policy):
    import ctypes
    import select
    import signal
    from eedfProviderTest import public_qualification_fixture, intended_signal_disposition

    patches = pytest.MonkeyPatch()
    provider, state = public_qualification_fixture(root / 'fixture', patches)
    state.pop('record_path')
    before = len(os.listdir('/proc/self/fd'))
    provider.qualify(state, 0.)
    owner = next(iter(provider._qualification_files.values()))
    notify_read, notify_write = os.pipe()
    ack_read, ack_write = os.pipe()
    deliveries, sentinels = [], []

    def handler(signum, frame):
        descriptors = [owner.parent_fd, owner.private_parent_fd, owner.directory_fd,
                       None if owner.stream is None else owner.stream.fileno()]
        released = owner._resources.released
        assert descriptors == [None] * 4 if released else all(fd is not None for fd in descriptors)
        if released:
            assert owner not in provider._qualification_files.values()
        deliveries.append(released)
        sentinels.append(os.open(root, os.O_RDONLY | os.O_DIRECTORY))
        owner.revoke()
        os.fstat(sentinels[-1])
        if policy == 'raising':
            raise KeyboardInterrupt('signal during native retirement')

    with intended_signal_disposition(int(number), handler):
        sender = os.fork()
        if sender == 0:
            os.close(notify_write)
            os.close(ack_read)
            ready, _, _ = select.select([notify_read], [], [], 15.)
            if not ready:
                os._exit(93)
            os.read(notify_read, 1)
            os.kill(os.getppid(), int(number))
            os.write(ack_write, b'x')
            os._exit(0)
        os.close(notify_read)
        os.close(ack_write)
        ctypes.CDLL(None).retention_arm(int(slot), int(after), notify_write, ack_read)
        try:
            try:
                provider.qualify(state, 1.)
            except BaseException:
                pass
        finally:
            ctypes.CDLL(None).retention_arm(-1, 0, -1, -1)
            _, status = os.waitpid(sender, 0)
            os.close(notify_write)
            os.close(ack_read)
    provider.record_terminal_failure(state, 2., RuntimeError('native probe complete'))
    survived = len(sentinels) == 1
    for descriptor in sentinels:
        os.fstat(descriptor)
        os.close(descriptor)
    result = dict(sender_exit=os.waitstatus_to_exitcode(status), deliveries=len(deliveries),
                  released_in_handler=deliveries, sentinel_survived=survived,
                  fd_delta=len(os.listdir('/proc/self/fd')) - before,
                  registry=len(provider._qualification_files))
    (root / 'retention.json').write_text(json.dumps(result, indent=2))
    patches.undo()


if __name__ == '__main__':
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    worker(Path(sys.argv[1]), *sys.argv[2:])
