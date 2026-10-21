Operating branches
==================

The branch scanner records the steady states of a plasma reactor at one absorbed power, continues
the stable ones over power under a stated certificate, and matches a run's terminal state to one
recorded branch.  See :ref:`plasmaoperatingbranch` for the input option that uses it.

Interrupted worker starts
-------------------------

When an application signal interrupts a worker start, the scanner waits up to
``WORKER_RECLAIM_TIMEOUT`` from the first cleanup entry for the worker's PID report.
It releases undecided children at that deadline and warns; it never kills a child
based only on copied memory. The consumer registered before the fork adopts any
unfinished release after the deadline plus 0.1 s. A released child that reports
within ``WORKER_DEFERRED_REAP_TIMEOUT`` is killed and reaped by its exact PID if it
is still the scanner's child. Each registration also has a ``RUN_LEASE`` of
60 seconds, renewed by the poll loop and worker-start steps. Missing cleanup
entry or retaining an exception traceback cannot extend this lease. If the
calling thread stalls past it, the consumer adopts the worker and a resumed scan
raises ``RuntimeError`` after cleanup; it never returns results from an adopted
worker. A worker held in an application fork hook can
therefore remain alive when ``scan()`` raises, then be reaped after reporting.
A child that dies or execs before reporting is left alone.

Cleanup attempts release warnings synchronously and retries interrupted attempts
against the same deadline. A warning may appear twice if interrupted after the
write. The consumer reaps before handing fallback warnings to a single daemon
emitter, so a blocked logging handler cannot prevent reaping. If no consumer
thread can start, a released worker is left alone and logged and may remain a
zombie until the process exits. If repeated interruptions are followed by an
immediate process exit, a warning may be lost; the worker is re-parented on exit.
If an application signal leaks a logging handler's lock, fallback warnings may
be delayed indefinitely and one emitter may remain blocked until exit; the
released worker is still reaped.

.. automodule:: rmgpy.tools.eedf.branches
    :members: scan, continue_branch, record_continuation, certificate_issue, tolerance_issue,
        reference_issue, matchable_issue, matching_branch, declared_branch, path_consistency,
        certified_coverage, aliasing_radius, PlasmaReactorAdapter, ReactorRun, Seed
