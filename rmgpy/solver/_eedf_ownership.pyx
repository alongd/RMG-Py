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

"""Explicit fd ownership across actual Python call/return boundaries.

The release section holds the GIL, calls no Python code, and uses only owned
integer slots and exact builtin registries. Python signal delivery therefore
observes a live bundle or a fully released one, never close-then-clear state.
"""

import os
import weakref
from types import MethodType

from cpython.dict cimport PyDict_GetItem, PyDict_DelItem
from cpython.list cimport PyList_GET_SIZE, PyList_GET_ITEM, PyList_SetSlice
from cpython.ref cimport PyObject
from cpython.weakref cimport PyWeakref_GetObject
from libc.errno cimport errno
from posix.fcntl cimport openat, fcntl, F_DUPFD_CLOEXEC, O_CLOEXEC, AT_FDCWD
from posix.unistd cimport close, write


cdef class Handles


cdef class Registry:
    cdef public dict files
    cdef dict exports
    cdef list scopes
    cdef object __weakref__

    def __cinit__(self):
        self.files = {}
        self.exports = {}
        self.scopes = []

    @property
    def pending_exports(self):
        return self.exports.values()

    def track(self, owner, Handles parent=None):
        cdef Handles handles
        cdef list scope
        cdef object active
        if not self.scopes:
            raise RuntimeError('file ownership requires an active call scope')
        handles = Handles()
        handles.registry_ref = weakref.ref(self)
        handles.owner_ref = weakref.ref(owner)
        handles.parent = parent
        scope = self.scopes[-1][0]
        if parent is not None:
            # The outer owner operation may retry capture settlement after a
            # nested signal interrupts it. Keep captures through that actual
            # operation's exit, not merely through the inner capture call.
            for active in self.scopes:
                if active[1] is parent:
                    scope = active[0]
                    break
        scope.append((handles, owner))
        handles._pin()
        if parent is not None:
            parent.captures.append(owner)
        return handles

    def register_record(self, str key, owner, Handles handles, adopt=False):
        if adopt:
            winner = self.files.setdefault(key, owner)
            if winner is not owner:
                return winner
        else:
            self.files[key] = owner
        handles.key = key
        return owner

    def register_export(self, owner, Handles handles):
        handles.export_key = id(owner)
        self.exports[handles.export_key] = owner

    cdef object _invoke(self, function, subject, tuple args, dict kwargs,
                        Handles use, bint retiring):
        cdef list created = []
        cdef Handles handles
        cdef object entry, failure = None, result
        cdef int error = 0, current
        self.scopes.append((created, use))
        if use is not None:
            use._pin()
            if retiring:
                use.retired = True
                use.retained = False
        try:
            result = function(subject, *args, **kwargs)
            return result
        except BaseException as caught:
            failure = caught
            raise
        finally:
            # The Python body's return trace has already happened. All entries
            # keep their owners alive while registry deletion decrefs them.
            self.scopes.pop()
            for entry in reversed(created):
                handles = entry[0]
                current = handles._unpin()
                if current and not error:
                    error = current
            if use is not None:
                current = use._unpin()
                if current and not error:
                    error = current
            if error:
                close_error = OSError(error, os.strerror(error))
                if failure is None:
                    raise close_error
                failure.__context__ = close_error


cdef class Handles:
    cdef int descriptors[4]
    cdef unsigned int users
    cdef bint retained, retired
    cdef public bint released
    cdef object registry_ref, owner_ref, key, export_key
    cdef Handles parent
    cdef public list captures

    def __cinit__(self):
        cdef int slot
        for slot in range(4):
            self.descriptors[slot] = -1
        self.captures = []

    cdef void _pin(self) noexcept:
        self.users += 1
        if self.parent is not None:
            self.parent._pin()

    cdef int _unpin(self) noexcept:
        cdef int error = 0, parent_error
        self.users -= 1
        if not self.users and not self.retained:
            error = self._release()
        if self.parent is not None:
            parent_error = self.parent._unpin()
            if not error:
                error = parent_error
        return error

    cdef int _release(self) noexcept:
        cdef int slot, descriptor, error = 0
        cdef Py_ssize_t index
        cdef object owner
        cdef Registry registry
        if self.released:
            return 0
        # These weakrefs have no callbacks. The scope/borrower holds the owner
        # strongly; this local also guards against destructor re-entry on del.
        owner = <object>PyWeakref_GetObject(self.owner_ref)
        registry = <Registry><object>PyWeakref_GetObject(self.registry_ref)
        for slot in range(4):
            descriptor = self.descriptors[slot]
            self.descriptors[slot] = -1
            if descriptor >= 0 and close(descriptor) < 0 and not error:
                error = errno
        if registry is not None:
            if self.key is not None and PyDict_GetItem(registry.files, self.key) == <PyObject*>owner:
                PyDict_DelItem(registry.files, self.key)
            if self.export_key is not None and PyDict_GetItem(registry.exports, self.export_key) == <PyObject*>owner:
                PyDict_DelItem(registry.exports, self.export_key)
        if self.parent is not None:
            for index in range(PyList_GET_SIZE(self.parent.captures)):
                if PyList_GET_ITEM(self.parent.captures, index) == <PyObject*>owner:
                    PyList_SetSlice(self.parent.captures, index, index + 1, <object>NULL)
                    break
        self.retained = False
        self.retired = True
        self.released = True
        return error

    def retain(self):
        if not self.retired and not self.released:
            self.retained = True

    def descriptor(self, int slot):
        if not 0 <= slot < 4:
            raise IndexError(slot)
        if self.descriptors[slot] < 0:
            return None
        return self.descriptors[slot]

    def open(self, int slot, path, int flags, int mode=0o600, dir_fd=None):
        cdef bytes name = os.fsencode(path)
        cdef int parent_fd = AT_FDCWD if dir_fd is None else dir_fd
        cdef int descriptor
        if not 0 <= slot < 4:
            raise IndexError(slot)
        if self.released or not self.users:
            raise RuntimeError('released file ownership refuses acquisition')
        if self.descriptors[slot] >= 0:
            raise RuntimeError('file descriptor slot is already owned')
        descriptor = openat(parent_fd, name, flags | O_CLOEXEC, mode)
        if descriptor < 0:
            raise OSError(errno, os.strerror(errno), os.fsdecode(name))
        # No Python call/signal check between acquisition and ownership.
        self.descriptors[slot] = descriptor
        return descriptor

    def duplicate(self, int target, int source):
        cdef int descriptor, old, error = 0
        if not 0 <= target < 4 or not 0 <= source < 4:
            raise IndexError('descriptor slot')
        if self.released or self.descriptors[source] < 0:
            raise ValueError('closed descriptor')
        descriptor = fcntl(self.descriptors[source], F_DUPFD_CLOEXEC, 0)
        if descriptor < 0:
            raise OSError(errno, os.strerror(errno))
        old = self.descriptors[target]
        self.descriptors[target] = descriptor
        if old >= 0 and close(old) < 0:
            error = errno
        if error:
            raise OSError(error, os.strerror(error))

    def fileno(self):
        if self.descriptors[3] < 0:
            raise ValueError('closed ownership stream')
        return self.descriptors[3]

    @property
    def closed(self):
        return self.descriptors[3] < 0

    def write(self, bytes payload):
        cdef Py_ssize_t result
        if self.descriptors[3] < 0:
            raise ValueError('closed ownership stream')
        result = write(self.descriptors[3], <const char*>payload, len(payload))
        if result < 0:
            raise OSError(errno, os.strerror(errno))
        return result

    def read(self):
        chunks = []
        while True:
            chunk = os.read(self.fileno(), 65536)
            if not chunk:
                return b''.join(chunks)
            chunks.append(chunk)

    def __dealloc__(self):
        # Only abandonment of still-current authority can reach this fallback.
        # Every terminal transition closes explicitly, even with external refs.
        cdef int slot
        for slot in range(4):
            if self.descriptors[slot] >= 0:
                close(self.descriptors[slot])


cdef class _ScopedCall:
    cdef object function
    cdef str mode

    def __init__(self, function, mode):
        self.function = function
        self.mode = mode

    def __get__(self, subject, cls):
        if subject is None:
            return self
        return MethodType(self, subject)

    @property
    def __code__(self):
        return self.function.__code__

    @property
    def __wrapped__(self):
        return self.function

    @property
    def __name__(self):
        return self.function.__name__

    def __call__(self, subject, *args, **kwargs):
        cdef Registry registry
        cdef Handles use = None
        if self.mode == 'scope':
            provider = getattr(subject, '_provider', subject)
            registry = provider._file_registry
            # Serialize scope bookkeeping as well as the Python body. A second
            # thread must not push a scope while waiting for the provider lock.
            with provider._qualification_lock:
                return registry._invoke(self.function, subject, args, kwargs, None, False)
        else:
            use = subject._resources
            if use.released:
                if self.mode == 'retire' or self.function.__name__ in ('settle', 'finish', '_remove_owned'):
                    return None
                raise RuntimeError('released file ownership refuses use')
            registry = <Registry><object>PyWeakref_GetObject(use.registry_ref)
        if registry is None:
            raise RuntimeError('file ownership registry is unavailable')
        return registry._invoke(self.function, subject, args, kwargs, use, self.mode == 'retire')


def ownership_scope(function):
    return _ScopedCall(function, 'scope')


def ownership_use(function):
    return _ScopedCall(function, 'use')


def ownership_retire(function):
    return _ScopedCall(function, 'retire')
