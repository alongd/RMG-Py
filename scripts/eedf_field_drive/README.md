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


## Reactor-impact tolerances

The optional generation-spec object `magnitude_tolerances` enables the
owner-approved strict rule. `magnitude-tolerances.json` records its parameters
and the selected 600 s horizon. Add that object to the measurement spec before
normal preparation/generation; it is copied into `row_inputs` and covered by
the normal table fingerprint and qualification setup identity. Old specs retain
their original comparisons. Existing artifacts are never relabelled accepted.

For nonelastic rates the absolute allowance is
`species_epsilon * n_floor / (ne_max * x * N * time_max_s)`, with
`n_floor = neutral_floor_fraction * N`,
`ne_max = electron_fraction_max * N`, and
`x = max(reference target, table target, target_fraction_floor)`.
Reverse rates use product fractions.
Here `N = P_Pa / (k_B * Tg_K)` uses query coordinates over envelope reference
values over nominal spec values, consistently for every impact scale.
Townsend/N uses `electron_epsilon / (abs(reference drift) * N * time_max_s)`.
Individual nonelastic powers use
`energy_epsilon * abs(reference field power)`.
The existing effective relative tolerance (including existing H3/H4
significance selection) is retained. Eligible checks use the larger of the original
`atol + rtol * abs(reference)` and the reactor-impact allowance plus the
relative term. Existing rate, flux and power floors therefore never tighten.
Zero reference drift or field power sets the impact allowance for that class
to zero. Underflow normalization still applies, so two values that normalize
to zero pass without a charge. Negative or nonfinite reference scales are rejected.

The original comparator always uses raw values, and every old pass costs zero
budget. Only the new allowance branch normalizes magnitudes below
`underflow_absolute` to zero, then compares against the larger of the old and
impact allowances, both evaluated with the normalized reference. A pass in
either branch admits the check.
Eligible channel kinds are excitation, vibrational,
rotational, ionization and attachment, including their reverse superelastic
rates. Elastic and effective coefficients/powers, aggregate power sums, T10,
H1, F0 and other transport checks keep their original arithmetic. Every check passing only
through the new absolute allowance contributes only its excess over the full
original allowance evaluated with the normalized reference to one of three
query budgets: species inventory divided by
`n_floor`, fractional electron change, or fractional field energy. All rate
constraints are conservatively summed together; duplicated H3 and consistency
constraints are counted independently. No cancellation is credited. A failed
cumulative check rejects the query through the existing all-checks aggregation.

The 18 prepared Ar+N2 reactor decks use `terminationTime=(600, 's')`.
The paired power controls in `probe_sweep.py` also end at 600 s.
**Horizon limitation:** future field/current waveform histories do not yet have
separate endpoint specifications. The largest available horizon is used;
recheck it before qualifying those histories, and regenerate identity if it
changes. The tolerance change alone does not establish a passing grid.
