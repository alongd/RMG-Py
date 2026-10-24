Operating branches
==================

The branch scanner records the steady states of a plasma reactor at one absorbed power, continues
the stable ones over power under a stated certificate, and matches a run's terminal state to one
recorded branch.  See :ref:`plasmaoperatingbranch` for the input option that uses it.

Native worker ownership
-----------------------

Worker scans require Linux x86_64, glibc 2.39, CPython 3.9.x, lock-free native
atomics and a successful per-worker watchdog handshake. Unsupported platforms,
multiprocessing modes or native startup failures raise ``NativeWorkerUnsupported``
or ``NativeWorkerStartupError``; there is no fallback lifetime mechanism. The
CPython patch and compiler versions are recorded in native startup diagnostics.
One process with no timeout runs in the caller, without worker isolation.

Each forked worker starts a pure-C watchdog before Python repairs its child state
or calls application child hooks. The calling scan thread renews ownership at
checkpoints; no background thread renews on its behalf. ``lease_s`` defaults to
60 seconds and must be finite and positive. Worker integration can run longer
while the caller keeps polling. A long honest callback or unpickle operation can
consume a lease: callers with known long operations can raise ``lease_s``.

Expiry raises ``ScanLeaseExpired`` with reason ``owner-unresponsive`` and native
record, deadline and reap facts. This is lifecycle cancellation, never a scientific
integration failure or a timeout recheck. The scan writes no new result artifact
and leaves any earlier file at the output path unchanged. Applications must not
interpret that older file as the failed call's result.

Results are accepted only after native normal completion, an exact-PID reap with
exit status zero, and a complete result envelope agree. Completion is committed
after the standard multiprocessing bootstrap and its finalizers finish. External
reaping makes the result unverifiable and raises ``NativeWorkerStatusError``.
Native cancellation requests worker self-exit; it does not signal a guessed PID.
Application child processes remain the application's responsibility.

An ``exec()`` inside a worker replaces its process image and ends the watchdog
guarantee. The replacement program can keep running past lease expiry or owner
death. RMG reports but does not kill that program. Results from it are rejected;
``NativeWorkerStatusError`` and the error log report the exact pending child PID
and bounded cleanup failure. A published terminal-state timestamp is not proof
of termination: only the exact reap records that the child exited. Earlier
artifacts remain unchanged.

A caller-thread native guard surrounds the complete worker-owning operation.
Before the scan returns or raises into caller code, every started worker is
exactly reaped or reported with its exact pending PID. Raising signals during
Python cleanup cannot bypass this final accounting, and retaining the exception
or its traceback does not postpone reaping. External reaping remains a typed
failure. No later scan, garbage collection or caller cleanup is required.

Native cleanup uses a single five-second allowance, never restarted by a retry.
Reaping precedes each reporting attempt. Application logging that blocks forever
or an endless raising-signal storm can delay scan return indefinitely; the guard
retains the reporting obligation throughout. Finite interruptions are retried.
A signal after a log append but before its return can produce duplicate lines.

A worker that expires while its owner is wedged stops running but can remain a
zombie until the owner resumes and reaps it, or exits and the operating system
reassigns it. Prompt reaping or logging by a wedged owner is not promised.
SIGKILL or power loss can leave no owner-side record or diagnostic; no durable
completion artifact is promised after those events.

Same-UID hostile tampering with records, the extension, workers or process state
is outside the guarantee. Ordinary application handlers, logging callbacks and
finalizers are covered by the tested native mechanism but can delay the caller
and cause lease cancellation. The startup evidence applies to the qualified
platform and successful handshake. A native startup hang stops further launches;
it is not an admitted worker with a lifetime guarantee. OS scheduling and
uninterruptible kernel waits can delay exit observation.

.. automodule:: rmgpy.tools.eedf.branches
    :members: scan, continue_branch, record_continuation, certificate_issue, tolerance_issue,
        reference_issue, matchable_issue, matching_branch, declared_branch, path_consistency,
        certified_coverage, aliasing_radius, PlasmaReactorAdapter, ReactorRun, Seed, ScanLeaseExpired,
        NativeWorkerUnsupported, NativeWorkerStartupError, NativeWorkerStatusError
