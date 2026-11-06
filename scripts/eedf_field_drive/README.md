# Field-drive table job

Copy this recipe into the field-drive-job child of a fresh upload root. It uses
the existing frozen bundle: RMG/LoKI/CMake/Eigen/Cython source archives, scientific
inputs and spec-template, plus the runtime/bootstrap/build support scripts.
The PBS job enters that child beneath PBS_O_WORKDIR.

The original LoKI source archive is retained. apply_loki_patch.py verifies the
upstream source and approved diff, changes only temporal growth mixing, and checks
the resulting source hash before compilation. No native binaries are uploaded.

The approved arithmetic rewrite avoids a subnormal rounding stall without changing
alpha, convergence tolerances, the iteration predicate, or the adapter's refusal.
growth-mixing-provenance.json records the trace, validation bounds, per-observable
maxima and twenty screen fields that differed at roundoff level. It applies only
to field-drive. I-377 keeps its original pinned solver.

After LoKI builds, the job runs a single screen-reference solve at frozen index18
before native RMG compilation. Logs go to results/loki-preflight. Every solver
refusal retains the last64KiB of stdout, stderr and setup under results/loki-refusals
and preserves the original exception. Scientific inputs and table qualification
gates remain unchanged.

The normal runtime spec hashes the actual node binary, shared libraries and build
cache. results/loki-build-provenance.json binds those pins to the reviewed source.
The provenance archive includes the solver shared library and patch metadata.
Fetch results only. Recompute model_inputs and use normal table admission after
extracting the provenance and table at their recorded paths; never edit acceptance
flags or fingerprints.

Rebuild source-commit and bundle SHA manifests after committing. Keep the complete
upload below14,500,000 bytes. Submission authorization, queue and cadence are managed
by the dispatch; this recipe does not submit or retry jobs.

Execution uses three processes with five OpenMP threads each, in a fifteen-CPU
allocation. The template must declare max_parallel=3 and timeout_s=6200.
Training scans, cold starts, original held-out solves and fresh T10 share the
same bounded process-pool implementation. Results are consumed in input order;
worker setup text is returned to the parent for normal manifest publication.

After the pinned LoKI build, the first scientific step measures serial and
three-process execution at the same six low/mid/high-field points used locally.
results/node-scaling.json records those timings, bounded numeric equivalence,
and the recomputed round-0 projection. The job exits cleanly before native RMG
compilation/generation if the projection exceeds16.8h or equivalence refuses.
Elapsed setup/probe time and a1309s remaining-build reserve are included.
The projection assumes the sample scales like the full rows and does not bound
unknown future adaptive-refinement work.

The execution-only comparison uses rtol1e-6/atol1e-20 for every reported EEDF,
anisotropy, swarm, power and rate value. Each dimensional power-balance residual
must be <=1e-20 in absolute value; its derived relative percentage is not compared.
These bounds do not alter any table qualification criterion. Five-thread serial
repeats already exhibit two variants at about1e-12 from OpenMP reduction ordering;
byte identity is not an execution acceptance gate.

results/row-timings.jsonl records UTC start/end, coordinate/fields, wall seconds
and child user/system CPU for every call, including refusals. A process executes
one solver subprocess at a time, so CPU deltas are per row. Append locking protects
complete JSON records from simultaneous worker writes.
