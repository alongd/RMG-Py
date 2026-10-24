/*
 * Copyright (c) 2026 Prof. William H. Green, Prof. Richard H. West and the RMG Team.
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#define _GNU_SOURCE 1
#include <Python.h>
#include <errno.h>
#include <math.h>
#include <pthread.h>
#include <signal.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/syscall.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>
#ifdef __GLIBC__
#include <gnu/libc-version.h>
#endif

/* Shared state is one lock-free phase/deadline word. All ledger pointers and
 * reap caches are parent-only, touched with the GIL by the owning OS thread.
 * A ledger entry survives a lost Python return value or Popen assignment.
 * The child watchdog never touches that ledger or calls a Python API.
 */
enum { NEW, STARTING, RUNNING, NORMAL, EXPIRED, OWNER_DEAD, CANCELLED, START_FAILED };
#define RECORD_ABI 1
#define PHASE_MASK UINT64_C(15)
#define MAX_TIME (UINT64_MAX >> 4)
typedef struct {
    _Atomic uint64_t state, ready, ended, hook, finish_entered, finish_gate;
    _Atomic int tid, thread_error;
    pid_t owner;
    uint64_t generation;
    int suppress_ready, create_error;
} Shared;
typedef struct Owner Owner;
typedef struct Record {
    struct Record *next;
    Owner *owner;
    Shared *shared;
    pid_t pid;
    int refs, reaped, status, retired;
    uint64_t created, renewed, retired_at, max_gap;
    uint64_t *trace;
    size_t trace_used, trace_size;
} Record;
struct Owner {
    Owner *next;
    Record *records;
    pid_t pid, tid;
    uint64_t lease, deadline, cleanup_deadline, generation;
    int refs, closed, abandoned, cancelled, tracing;
};
static Owner *owners;
static pid_t native_child;
static uint64_t next_generation;
static const struct timespec tick = {0, 1000000};
static const char *phase_names[] = {
    "new", "starting", "running", "normal", "owner-unresponsive",
    "owner-dead", "cancelled", "startup-failed"
};

static uint64_t now_us(void) {
    struct timespec t;
    if (clock_gettime(CLOCK_MONOTONIC, &t)) _exit(127);
    return (uint64_t)t.tv_sec * UINT64_C(1000000) + t.tv_nsec / 1000;
}
static unsigned phase(uint64_t word) { return word & PHASE_MASK; }
static uint64_t deadline(uint64_t word) { return word >> 4; }
static uint64_t word(uint64_t until, unsigned state) { return (until << 4) | state; }
static int exit_code(unsigned state) {
    return state == EXPIRED ? 124 : state == OWNER_DEAD ? 125 :
           state == CANCELLED ? 126 : 127;
}
static int seconds_us(double seconds, uint64_t *out) {
    if (!isfinite(seconds) || seconds <= 0 || seconds * 1e6 > (double)(MAX_TIME - now_us())) {
        PyErr_SetString(PyExc_ValueError, "lease/deadline must be finite, positive and representable");
        return 0;
    }
    *out = (uint64_t)ceil(seconds * 1e6);
    return 1;
}
static unsigned stop(Shared *s, unsigned reason) {
    uint64_t old = atomic_load_explicit(&s->state, memory_order_acquire);
    while (phase(old) <= RUNNING) {
        uint64_t desired = word(deadline(old), reason);
        if (atomic_compare_exchange_weak_explicit(&s->state, &old, desired,
                memory_order_acq_rel, memory_order_acquire)) {
            atomic_store_explicit(&s->ended, now_us(), memory_order_release);
            return reason;
        }
    }
    return phase(old);
}
static void *watchdog(void *arg) {
    Shared *s = arg;
    atomic_store(&s->tid, (int)syscall(SYS_gettid));
    if (!s->suppress_ready)
        atomic_store_explicit(&s->ready, now_us(), memory_order_release);
    for (;;) {
        uint64_t old = atomic_load_explicit(&s->state, memory_order_acquire);
        unsigned p = phase(old);
        if (p == NORMAL) return NULL;
        if (p >= EXPIRED) _exit(exit_code(p));
        unsigned reason = 0;
        if (getppid() != s->owner) reason = OWNER_DEAD;
        else if (now_us() >= deadline(old)) reason = p == STARTING ? START_FAILED : EXPIRED;
        if (reason) {
            /* Renewal may have replaced the observed deadline. Never expire
             * that newer lease using an observation of the old word. */
            if (atomic_compare_exchange_strong_explicit(&s->state, &old,
                    word(deadline(old), reason), memory_order_acq_rel, memory_order_acquire)) {
                atomic_store_explicit(&s->ended, now_us(), memory_order_release);
                _exit(exit_code(reason));
            }
            continue;
        }
        nanosleep(&tick, NULL);
    }
}
static int owned(Owner *o) {
    if (o->pid != getpid() || o->tid != (pid_t)syscall(SYS_gettid)) {
        PyErr_SetString(PyExc_RuntimeError, "native ownership belongs to another process or calling thread");
        return 0;
    }
    return 1;
}
static Owner *get_owner(PyObject *obj) { return PyCapsule_GetPointer(obj, "eedf.owner"); }
static Record *get_record(PyObject *obj) { return PyCapsule_GetPointer(obj, "eedf.record"); }
static void trace(Record *r, uint64_t stamp) {
    if (r->renewed && stamp - r->renewed > r->max_gap) r->max_gap = stamp - r->renewed;
    r->renewed = stamp;
    if (!r->owner->tracing) return;
    if (r->trace_used == r->trace_size) {
        size_t size = r->trace_size ? r->trace_size * 2 : 64;
        uint64_t *data = realloc(r->trace, size * sizeof(*data));
        if (!data) { r->owner->tracing = -1; return; }
        r->trace = data;
        r->trace_size = size;
    }
    r->trace[r->trace_used++] = stamp;
}
static int poll_record(Record *r) {
    if (!r->pid || r->reaped) return r->reaped;
    int status;
    pid_t got = waitpid(r->pid, &status, WNOHANG);
    if (got == r->pid) { r->status = status; r->reaped = 1; }
    else if (got < 0 && errno == ECHILD) r->reaped = 2;
    else if (got < 0 && errno != EINTR) return -1;
    return r->reaped;
}
static void cancel_owner(Owner *o, unsigned reason) {
    if (!o->cancelled) o->cancelled = reason;
    for (Record *r = o->records; r; r = r->next)
        if (!r->retired) stop(r->shared, reason);
}
static int active(Owner *o) {
    for (Record *r = o->records; r; r = r->next)
        if (r->pid && !r->retired) return 1;
    return 0;
}
static int check_owner(Owner *o, uint64_t now) {
    if (!o->cancelled && o->deadline && active(o) && now >= o->deadline)
        cancel_owner(o, EXPIRED);
    for (Record *r = o->records; r && !o->cancelled; r = r->next) {
        unsigned p = phase(atomic_load_explicit(&r->shared->state, memory_order_acquire));
        if (!r->retired && (p == EXPIRED || p == OWNER_DEAD || p == START_FAILED))
            cancel_owner(o, p);
    }
    return !o->cancelled;
}
static void collect(Owner *o) {
    if (o->refs || !o->closed) return;
    for (Record *r = o->records; r; r = r->next)
        if (r->refs || (r->pid && !r->reaped)) return;
    Owner **link = &owners;
    while (*link && *link != o) link = &(*link)->next;
    if (*link) *link = o->next;
    Record *r = o->records;
    while (r) {
        Record *next = r->next;
        munmap(r->shared, sizeof(*r->shared));
        free(r->trace); free(r); r = next;
    }
    free(o);
}
static void drop_owner(PyObject *obj) {
    Owner *o = get_owner(obj);
    if (!o) { PyErr_Clear(); return; }
    o->refs--;
    if (o->pid == getpid() && !o->closed) {
        o->abandoned = o->closed = 1;
        cancel_owner(o, CANCELLED);
        if (o->tid == (pid_t)syscall(SYS_gettid))
            for (Record *r = o->records; r; r = r->next) poll_record(r);
    }
    collect(o);
}
static void drop_record(PyObject *obj) {
    Record *r = get_record(obj);
    if (!r) { PyErr_Clear(); return; }
    r->refs--;
    collect(r->owner);
}
static PyObject *platform_info(PyObject *self, PyObject *unused) {
    _Atomic uint64_t test64 = 0;
    _Atomic int testint = 0;
    int supported = 0;
    const char *libc = "not-glibc";
#ifdef __GLIBC__
    libc = gnu_get_libc_version();
#if defined(__linux__) && defined(__x86_64__) && PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION == 9
    supported = !strcmp(libc, "2.39") && atomic_is_lock_free(&test64) && atomic_is_lock_free(&testint);
#endif
#endif
    return Py_BuildValue("{s:O,s:s,s:s,s:s,s:i,s:O}", "supported", supported ? Py_True : Py_False,
        "glibc", libc, "python", Py_GetVersion(), "compiler", __VERSION__,
        "record_abi", RECORD_ABI, "lock_free", atomic_is_lock_free(&test64) && atomic_is_lock_free(&testint) ? Py_True : Py_False);
}
static PyObject *new_owner(PyObject *self, PyObject *args) {
    double seconds; int tracing = 0; uint64_t lease;
    if (!PyArg_ParseTuple(args, "d|p", &seconds, &tracing) || !seconds_us(seconds, &lease)) return NULL;
    /* Previously abandoned scans are reaped only by their actual parent thread. */
    for (Owner *o = owners, *next; o; o = next) {
        next = o->next;
        if (o->abandoned && o->pid == getpid() && o->tid == (pid_t)syscall(SYS_gettid)) {
            for (Record *r = o->records; r; r = r->next) poll_record(r);
            collect(o);
        }
    }
    Owner *o = calloc(1, sizeof(*o));
    if (!o) return PyErr_NoMemory();
    o->pid = getpid(); o->tid = syscall(SYS_gettid); o->lease = lease;
    o->refs = 1; o->tracing = tracing; o->generation = ++next_generation;
    o->next = owners; owners = o;
    PyObject *obj = PyCapsule_New(o, "eedf.owner", drop_owner);
    if (!obj) { o->refs = 0; o->closed = 1; collect(o); }
    return obj;
}
static PyObject *new_record(PyObject *self, PyObject *obj) {
    Owner *o = get_owner(obj);
    if (!o || !owned(o)) return NULL;
    if (o->closed || !check_owner(o, now_us())) {
        PyErr_SetString(PyExc_RuntimeError, "native owner is cancelled or closed"); return NULL;
    }
    Record *r = calloc(1, sizeof(*r));
    if (!r) return PyErr_NoMemory();
    r->shared = mmap(NULL, sizeof(*r->shared), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
    if (r->shared == MAP_FAILED) { free(r); return PyErr_SetFromErrno(PyExc_OSError); }
    r->owner = o; r->refs = 1; r->shared->owner = o->pid;
    r->shared->generation = ++next_generation;
    atomic_init(&r->shared->state, word(0, NEW));
    atomic_init(&r->shared->ready, 0); atomic_init(&r->shared->ended, 0);
    atomic_init(&r->shared->hook, 0); atomic_init(&r->shared->finish_entered, 0);
    atomic_init(&r->shared->finish_gate, 0); atomic_init(&r->shared->tid, 0);
    atomic_init(&r->shared->thread_error, 0);
    r->next = o->records; o->records = r;
    PyObject *record = PyCapsule_New(r, "eedf.record", drop_record);
    if (!record) r->refs = 0;
    return record;
}
static PyObject *native_fork(PyObject *self, PyObject *obj) {
    Record *r = get_record(obj);
    if (!r || !owned(r->owner)) return NULL;
    Owner *o = r->owner; Shared *s = r->shared;
    if (r->pid || phase(atomic_load(&s->state)) != NEW || o->closed || !check_owner(o, now_us())) {
        PyErr_SetString(PyExc_RuntimeError, "native record cannot start"); return NULL;
    }
    if (PySys_Audit("os.fork", NULL) < 0) return NULL;
    if (!owned(o)) return NULL;
    PyOS_BeforeFork();
    if (o->pid != getpid() || o->tid != (pid_t)syscall(SYS_gettid)) {
        /* Balance interpreter preparation before refusing an application copy.
         * Parent callbacks can fork too; validate again after they return. */
        PyOS_AfterFork_Parent();
        owned(o);
        return NULL;
    }
    uint64_t now = now_us();
    if (!check_owner(o, now)) {
        PyOS_AfterFork_Parent();
        PyErr_SetString(PyExc_RuntimeError, "owner lease expired in before-fork callbacks"); return NULL;
    }
    r->created = now;
    atomic_store_explicit(&s->state, word(now + 1000000, STARTING), memory_order_release);
    pid_t pid = fork(); int saved_errno = errno;
    if (pid == 0) {
        native_child = getpid();
        pthread_t thread;
        int rc = s->create_error ? s->create_error : pthread_create(&thread, NULL, watchdog, s);
        atomic_store(&s->thread_error, rc);
        if (rc) { stop(s, START_FAILED); _exit(127); }
        while (!atomic_load_explicit(&s->ready, memory_order_acquire)) nanosleep(&tick, NULL);
        while (phase(atomic_load_explicit(&s->state, memory_order_acquire)) == STARTING)
            nanosleep(&tick, NULL);
        const pid_t child_pid = getpid(), child_tid = (pid_t)syscall(SYS_gettid);
        PyOS_AfterFork_Child();
        if (getpid() != child_pid || (pid_t)syscall(SYS_gettid) != child_tid) {
            PyErr_SetString(PyExc_RuntimeError, "native ownership belongs to another process or calling thread");
            return NULL;
        }
        return PyLong_FromLong(0);
    }
    /* This store precedes every parent application hook and Python bytecode. */
    if (pid > 0) r->pid = pid;
    int admitted = 0;
    if (pid > 0) {
        while (!(admitted = atomic_load_explicit(&s->ready, memory_order_acquire) != 0)
                && phase(atomic_load(&s->state)) == STARTING && now_us() < now + 1000000)
            nanosleep(&tick, NULL);
        uint64_t old = atomic_load_explicit(&s->state, memory_order_acquire);
        uint64_t admission = now_us();
        admitted = admitted && phase(old) == STARTING && admission < deadline(old)
            && atomic_compare_exchange_strong_explicit(&s->state, &old,
                word(admission + o->lease, RUNNING), memory_order_acq_rel, memory_order_acquire);
        if (admitted) { o->deadline = admission + o->lease; trace(r, admission); }
        else { stop(s, START_FAILED); cancel_owner(o, START_FAILED); poll_record(r); }
    } else stop(s, START_FAILED);
    PyOS_AfterFork_Parent();
    if (pid < 0) { errno = saved_errno; return PyErr_SetFromErrno(PyExc_OSError); }
    if (PyErr_CheckSignals() < 0) return NULL;
    if (!admitted) { PyErr_SetString(PyExc_RuntimeError, "native watchdog startup handshake failed"); return NULL; }
    return PyLong_FromLong(pid);
}
static PyObject *renew_owner(PyObject *self, PyObject *obj) {
    Owner *o = get_owner(obj);
    if (!o || !owned(o)) return NULL;
    uint64_t now = now_us();
    if (o->closed || !check_owner(o, now)) Py_RETURN_FALSE;
    for (Record *r = o->records; r; r = r->next) {
        if (!r->pid || r->retired) continue;
        uint64_t old = atomic_load_explicit(&r->shared->state, memory_order_acquire);
        for (;;) {
            unsigned p = phase(old);
            if (p == CANCELLED) break;  /* Explicit per-seed timeout, not lease loss. */
            if ((p != RUNNING && p != NORMAL) || now >= deadline(old)) {
                cancel_owner(o, p >= EXPIRED ? p : EXPIRED); Py_RETURN_FALSE;
            }
            if (atomic_compare_exchange_weak_explicit(&r->shared->state, &old,
                    word(now + o->lease, p), memory_order_acq_rel, memory_order_acquire)) break;
        }
        trace(r, now);
    }
    o->deadline = active(o) ? now + o->lease : 0;
    Py_RETURN_TRUE;
}
static PyObject *finish(PyObject *self, PyObject *args) {
    PyObject *obj; int code;
    if (!PyArg_ParseTuple(args, "Oi", &obj, &code)) return NULL;
    Record *r = get_record(obj); if (!r) return NULL;
    if (native_child != getpid()) { PyErr_SetString(PyExc_RuntimeError, "finish is child-only"); return NULL; }
    if (r->shared->owner != getppid()) { stop(r->shared, OWNER_DEAD); _exit(125); }
    Shared *s = r->shared;
    atomic_store(&s->finish_entered, now_us());
    while (atomic_load_explicit(&s->finish_gate, memory_order_acquire)) nanosleep(&tick, NULL);
    uint64_t old = atomic_load_explicit(&s->state, memory_order_acquire);
    for (;;) {
        if (phase(old) != RUNNING) _exit(exit_code(phase(old)));
        if (now_us() >= deadline(old)) {
            if (atomic_compare_exchange_strong_explicit(&s->state, &old,
                    word(deadline(old), EXPIRED), memory_order_acq_rel, memory_order_acquire)) {
                atomic_store_explicit(&s->ended, now_us(), memory_order_release);
                _exit(124);
            }
            continue;
        }
        if (getppid() != s->owner) { _exit(exit_code(stop(s, OWNER_DEAD))); }
        if (atomic_compare_exchange_weak_explicit(&s->state, &old, word(deadline(old), NORMAL),
                memory_order_acq_rel, memory_order_acquire)) {
            /* No Python call, finalizer, signal dispatch or GIL reacquisition. */
            atomic_store_explicit(&s->ended, now_us(), memory_order_release);
            _exit(code);
        }
    }
}
static PyObject *info_record(Record *r) {
    Shared *s = r->shared; uint64_t value = atomic_load_explicit(&s->state, memory_order_acquire);
    return Py_BuildValue("{s:i,s:i,s:K,s:s,s:i,s:K,s:K,s:K,s:K,s:K,s:K,s:K,s:i,s:i,s:i,s:i,s:i,s:K,s:K}",
        "owner", s->owner, "pid", r->pid, "generation", s->generation,
        "phase", phase_names[phase(value)], "reason", phase(value), "deadline_us", deadline(value),
        "created_us", r->created, "renewed_us", r->renewed, "ready_us", atomic_load(&s->ready),
        "ended_us", atomic_load(&s->ended), "hook_us", atomic_load(&s->hook), "max_gap_us", r->max_gap,
        "reaped", r->reaped, "wait_status", r->status, "thread_error", atomic_load(&s->thread_error),
        "tid", atomic_load(&s->tid), "retired", r->retired, "retired_us", r->retired_at,
        "finish_entered_us", atomic_load(&s->finish_entered));
}
static PyObject *info(PyObject *self, PyObject *obj) {
    Record *r = get_record(obj); return r ? info_record(r) : NULL;
}
static PyObject *owner_info(PyObject *self, PyObject *obj) {
    Owner *o = get_owner(obj); if (!o || !owned(o)) return NULL;
    PyObject *records = PyList_New(0); if (!records) return NULL;
    for (Record *r = o->records; r; r = r->next) {
        PyObject *entry = info_record(r);
        if (!entry || PyList_Append(records, entry)) { Py_XDECREF(entry); Py_DECREF(records); return NULL; }
        Py_DECREF(entry);
    }
    return Py_BuildValue("{s:i,s:i,s:K,s:K,s:i,s:i,s:N}", "reason", o->cancelled, "closed", o->closed,
        "deadline_us", o->deadline, "cleanup_deadline_us", o->cleanup_deadline,
        "trace_complete", o->tracing >= 0, "pid", o->pid, "records", records);
}
static PyObject *poll(PyObject *self, PyObject *obj) {
    Record *r = get_record(obj); if (!r) return NULL;
    if (r->owner->pid != getpid()) {
        PyErr_SetString(PyExc_RuntimeError, "native record belongs to another process"); return NULL;
    }
    /* multiprocessing's global child census can run in an unrelated application
     * thread. It may observe our cached status, but only the scan caller reaps. */
    if (r->owner->tid != (pid_t)syscall(SYS_gettid)) return PyLong_FromLong(r->reaped);
    int status = poll_record(r); if (status < 0) return PyErr_SetFromErrno(PyExc_OSError);
    return PyLong_FromLong(status);
}
static PyObject *cancel(PyObject *self, PyObject *obj) {
    Record *r = get_record(obj); if (!r || !owned(r->owner)) return NULL;
    stop(r->shared, CANCELLED); Py_RETURN_NONE;
}
static PyObject *retire(PyObject *self, PyObject *obj) {
    Record *r = get_record(obj); if (!r || !owned(r->owner)) return NULL;
    if (!r->reaped) { PyErr_SetString(PyExc_RuntimeError, "cannot retire an unreaped native record"); return NULL; }
    uint64_t now = now_us();
    /* The last Python checkpoint can be interrupted by an honest callback.
     * Validate the lease in the same C operation that removes its final
     * active record; retirement must never erase an elapsed deadline. */
    if (r->owner->closed || !check_owner(r->owner, now)) Py_RETURN_FALSE;
    r->retired_at = now; trace(r, now); r->retired = 1;
    if (!active(r->owner)) r->owner->deadline = 0;
    Py_RETURN_TRUE;
}
typedef struct { int pending, external, error; } Sweep;

/* No Python API, allocation or user reference destruction in this region.
 * The guarded call enters it even when Python's cleanup retry itself raised. */
static Sweep cleanup_owner(Owner *o, uint64_t interval) {
    if (!o->cleanup_deadline) {
        o->cleanup_deadline = now_us() + interval;
        /* A local cancellation/reap attempt may already have spent this
         * allowance. Entering finally must not restart its deadline. */
        for (Record *r = o->records; r; r = r->next) {
            uint64_t ended = atomic_load_explicit(&r->shared->ended, memory_order_acquire);
            if (r->pid && !r->retired && !r->reaped && ended
                    && ended + interval < o->cleanup_deadline)
                o->cleanup_deadline = ended + interval;
        }
    }
    for (Record *r = o->records; r; r = r->next) {
        if (!r->retired) { cancel_owner(o, CANCELLED); break; }
    }
    Sweep sweep;
    do {
        sweep = (Sweep){0, 0, 0};
        for (Record *r = o->records; r; r = r->next) {
            if (poll_record(r) < 0 && !sweep.error) sweep.error = errno;
            if (r->pid && !r->reaped) sweep.pending++;
            if (r->reaped == 2) sweep.external++;
        }
        if (sweep.pending && now_us() < o->cleanup_deadline) nanosleep(&tick, NULL);
        else break;
    } while (sweep.pending);
    if (!sweep.pending) o->closed = 1;
    return sweep;
}
static PyObject *cleanup(PyObject *self, PyObject *args) {
    PyObject *obj; double seconds = 5.; uint64_t interval;
    if (!PyArg_ParseTuple(args, "O|d", &obj, &seconds) || !seconds_us(seconds, &interval)) return NULL;
    Owner *o = get_owner(obj); if (!o || !owned(o)) return NULL;
    Sweep sweep = cleanup_owner(o, interval);
    if (sweep.error) { errno = sweep.error; return PyErr_SetFromErrno(PyExc_OSError); }
    return PyBool_FromLong(!sweep.pending);
}
static PyObject *guard(PyObject *self, PyObject *args) {
    PyObject *obj, *body, *report; double seconds; uint64_t interval;
    if (!PyArg_ParseTuple(args, "OOOd", &obj, &body, &report, &seconds)
            || !seconds_us(seconds, &interval)) return NULL;
    Owner *o = get_owner(obj); if (!o || !owned(o)) return NULL;
    if (o->closed || !PyCallable_Check(body) || !PyCallable_Check(report)) {
        PyErr_SetString(PyExc_ValueError, "guard requires an open owner and callable body/reporter");
        return NULL;
    }
    /* args retains the capsule and callbacks until this C call returns. There
     * is no worker before the body is entered. Save its exception without
     * normalization or decrefs, then settle before Python can run again. */
    PyObject *result = PyObject_CallNoArgs(body);
    PyObject *error_type = NULL, *error_value = NULL, *error_tb = NULL;
    PyErr_Fetch(&error_type, &error_value, &error_tb);
    for (;;) {
        if (o->pid != getpid() || o->tid != (pid_t)syscall(SYS_gettid)) {
            /* A reporter or its exception destructor can fork too. Revalidate
             * before every raw sweep/retry: a copied stack has no authority. */
            Py_XDECREF(result); Py_XDECREF(error_type); Py_XDECREF(error_value); Py_XDECREF(error_tb);
            owned(o);
            return NULL;
        }
        Sweep sweep = cleanup_owner(o, interval);
        if (!sweep.pending && !sweep.external && !sweep.error) break;
        /* Reap precedes every reporting attempt. A raising Python logger or
         * signal cannot escape C retry bookkeeping or restart the deadline.
         * A permanently blocked/interrupted logger can delay return forever. */
        PyObject *diagnostic = PyObject_CallFunction(report, "i", sweep.error);
        if (diagnostic && PyExceptionInstance_Check(diagnostic)) {
            if (error_value && PyExceptionInstance_Check(error_value)) {
                if (error_tb) PyException_SetTraceback(error_value, error_tb);
                PyException_SetContext(diagnostic, error_value);
                error_value = NULL;  /* Context steals the saved value. */
            }
            Py_XDECREF(result); Py_XDECREF(error_type); Py_XDECREF(error_value); Py_XDECREF(error_tb);
            PyErr_SetObject((PyObject *)Py_TYPE(diagnostic), diagnostic);
            Py_DECREF(diagnostic);
            return NULL;
        }
        Py_XDECREF(diagnostic);
        PyErr_Clear();
    }
    /* All started records are accounted for before releasing any user object
     * or propagating an exception whose traceback may retain the whole scan. */
    if (error_type) {
        Py_XDECREF(result);
        PyErr_Restore(error_type, error_value, error_tb);
        return NULL;
    }
    return result;
}
static PyObject *record_trace(PyObject *self, PyObject *obj) {
    Record *r = get_record(obj); if (!r || !owned(r->owner)) return NULL;
    PyObject *list = PyList_New(r->trace_used); if (!list) return NULL;
    for (size_t i = 0; i < r->trace_used; i++) {
        PyObject *value = PyLong_FromUnsignedLongLong(r->trace[i]);
        if (!value) { Py_DECREF(list); return NULL; }
        PyList_SET_ITEM(list, i, value);
    }
    return list;
}
static PyObject *is_worker(PyObject *self, PyObject *unused) { return PyBool_FromLong(native_child == getpid()); }
/* Private protocol-test seams: no hooks are executed by the watchdog. */
static PyObject *test_config(PyObject *self, PyObject *args) {
    PyObject *obj; const char *key; int value;
    if (!PyArg_ParseTuple(args, "Osi", &obj, &key, &value)) return NULL;
    Record *r = get_record(obj); if (!r || !owned(r->owner)) return NULL;
    if (!strcmp(key, "finish_gate")) atomic_store(&r->shared->finish_gate, value);
    else if (r->pid) { PyErr_SetString(PyExc_RuntimeError, "startup seam must precede fork"); return NULL; }
    else if (!strcmp(key, "suppress_ready")) r->shared->suppress_ready = value;
    else if (!strcmp(key, "create_error")) r->shared->create_error = value;
    else { PyErr_SetString(PyExc_ValueError, "unknown native test seam"); return NULL; }
    Py_RETURN_NONE;
}
static PyObject *mark_hook(PyObject *self, PyObject *obj) {
    Record *r = get_record(obj); if (!r) return NULL;
    atomic_store(&r->shared->hook, now_us()); Py_RETURN_NONE;
}
static PyObject *hold_gil(PyObject *self, PyObject *args) {
    double seconds; uint64_t duration;
    if (!PyArg_ParseTuple(args, "d", &seconds) || !seconds_us(seconds, &duration)) return NULL;
    uint64_t until = now_us() + duration; volatile double y = 1.;
    while (now_us() < until) {
        for (int i = 0; i < 10000; i++) {
            double h = 1e-6, k1 = -y, k2 = -(y + h*k1/2), k3 = -(y + h*k2/2), k4 = -(y + h*k3);
            y += h*(k1 + 2*k2 + 2*k3 + k4)/6;
        }
    }
    return PyFloat_FromDouble(y);
}
static PyMethodDef methods[] = {
    {"check_platform", platform_info, METH_NOARGS, "Native runtime qualification facts."},
    {"owner", new_owner, METH_VARARGS, "Create a caller-thread ownership ledger."},
    {"record", new_record, METH_O, "Link a durable record before any fork."},
    {"fork", native_fork, METH_O, "Publish PID and admit watchdog before Python repair/hooks."},
    {"renew", renew_owner, METH_O, "Irreversibly reject expiry; renew this owner's active records."},
    {"finish", finish, METH_VARARGS, "Commit normal completion and _exit without returning to Python."},
    {"info", info, METH_O, "Snapshot a record, without changing ownership."},
    {"owner_info", owner_info, METH_O, "Snapshot the native ledger."},
    {"poll", poll, METH_O, "Exact-PID wait and native status cache."},
    {"cancel", cancel, METH_O, "Ask this worker to exit through its native watchdog."},
    {"retire", retire, METH_O, "Retire an exactly accounted-for record."},
    {"cleanup", cleanup, METH_VARARGS, "Publish cancellation then reap to a fixed deadline."},
    {"guard", guard, METH_VARARGS, "Run a worker-owning call and settle before any exception can escape."},
    {"trace", record_trace, METH_O, "Read opt-in native checkpoint timestamps."},
    {"is_worker", is_worker, METH_NOARGS, "Native child-role identity, never inherited by an app fork."},
    {"_test_config", test_config, METH_VARARGS, "Private startup/completion test seams."},
    {"_mark_hook", mark_hook, METH_O, "Private native hook timestamp."},
    {"_hold_gil", hold_gil, METH_VARARGS, "Private RK4 integration holding the GIL."},
    {NULL, NULL, 0, NULL}
};
static struct PyModuleDef module = {PyModuleDef_HEAD_INIT, "_native_worker", NULL, -1, methods};
PyMODINIT_FUNC PyInit__native_worker(void) { return PyModule_Create(&module); }
