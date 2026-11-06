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
