# I048 reproduction

Run in the dispatched worktree. The commands use only the allowed scientific
fixtures and named scratch; never inspect an excluded dataset. `common.py`
contains the declarations. `legacy_cli.py` and `select_pool.py` adapt the
unchanged I043 scientific routines rather than maintaining another rotor model.
DFT uses native OpenMP threads within three assigned 3/3/2-core sets,
with one BLAS thread explicitly enforced through the already-installed pure
Python threadpool controller in `rmg_env`. Three DFT jobs can run concurrently.
The one- versus four-native-thread cost pilot reproduced the same energy and
showed only about 1.25 cores average CPU use in the four-core allocation.
The three-lane allocation uses the same declared energy level and selections.

```bash
cd /home/alon/Code/RMG-Py-kmc-i048-oligomer-series
I048_PYTHON=/home/alon/anaconda3/envs/rmg_env/bin/python
I048_SCRATCH=/home/alon/runs/i048-oligomer-series
I048_PROBE=$PWD/test/rmgpy/kmc/fixtures/i048_probe
I048_CPUS=0,2,4,6,8,10,12,14
export PYTHONPATH=$PWD
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR=$I048_SCRATCH/mpl
mkdir -p "$I048_SCRATCH/logs"
i048_run() {
  local I048_LOG=$1
  shift
  taskset -c "$I048_CPUS" "$@" \
    > >(taskset -c "$I048_CPUS" tee -a "$I048_SCRATCH/logs/$I048_LOG.stdout.log") \
    2> >(taskset -c "$I048_CPUS" tee -a "$I048_SCRATCH/logs/$I048_LOG.stderr.log" >&2)
}
i048_run build "$I048_PYTHON" setup.py build_ext --inplace
i048_run prepare "$I048_PYTHON" "$I048_PROBE/run_series.py" --prepare
i048_run baseline "$I048_PYTHON" "$I048_PROBE/baseline.py"
i048_run search "$I048_PYTHON" "$I048_PROBE/run_series.py" --stage search
i048_run minima "$I048_PYTHON" "$I048_PROBE/run_series.py" --stage minima
i048_run electronic "$I048_PYTHON" "$I048_PROBE/run_series.py" --stage electronic
i048_run rotors "$I048_PYTHON" "$I048_PROBE/run_series.py" --stage rotors
i048_run base-thermo "$I048_PYTHON" "$I048_PROBE/thermochemistry.py" --samples 1024 --no-roots
i048_run thermo "$I048_PYTHON" "$I048_PROBE/thermochemistry.py" --samples 2048
i048_run literature "$I048_PYTHON" "$I048_PROBE/literature.py"
i048_run report "$I048_PYTHON" "$I048_PROBE/report.py" --write-report
i048_run verify "$I048_PYTHON" "$I048_PROBE/verify_results.py" --replay
```

The initial run started `pipeline.py --after-search` while the searches were
running to queue the subsequent producer stages. When the first new search
completed, it was restarted with `--overlap-electronic`: one completed case at
a time receives minima checks and electronic energies on the same eight-core
set while other searches continue. The bulk stages follow synchronously,
so no two producers write the same case concurrently. Stages are restartable and
reuse their saved inputs; do not delete evidence. Preparation refuses changed
inputs or methods. Source-only and available-only options are explicitly
development checks and do not satisfy the dispatch Verifier.

The primary queue was later restarted while its preparation CLI was waiting
on a supplemental-owned case. The interrupted waiting job and original logs
are retained in scratch. The restarted queue writes `pipeline-routed` logs,
checks the per-case locks before selecting work, and uses `--if-idle` for its
preparation children. Status 75 retains a pending case without claiming it is
complete. Producing calculations and the original wall deadline continue.

After its electronic calculations completed, the first tetramer also received
its declared rotor integrations while searches continued. Its command was:

```bash
PYTHONPATH="$PWD:$I048_PROBE" i048_run early-rotor "$I048_PYTHON" \
  -c 'from run_series import rotor; rotor("ps4_0001",6)'
```

This uses core 6 within the same allocation. A per-case process lock prevents
overlapping rotor producers or replays from writing the same evidence.

The measured search CPU utilization supported adding a second preparation lane:

```bash
i048_run supplemental "$I048_PYTHON" "$I048_PROBE/supplemental.py"
```

It checks minima on core 10 and calculates electronic energies on cores 6, 8
and 10, within the same eight physical cores used by all production jobs.
Per-case locks serialize minimum and electronic producers. The supplemental
producer finishes its active case when all searches complete, then releases a
process lock that the bulk electronic stage must acquire before launching its
three lanes. The previously active `ps4_0000` is excluded from supplemental
scheduling because it started before these case locks were installed. No
scientific reduction, energy level or sequence weight changes.

Later primary preparation splits independent minimum checks across two
single-core lanes. Bulk minimum preparation divides the same eight cores
among cases whose pools are still missing. Candidate batches and their cores
are saved in each case's `checks/execution.json`; annotation of xTB execution
metadata is restricted to the candidates owned by that lane. The original
supplemental process retains its single minimum-check lane. Candidate selection,
optimization thresholds, Hessians and electronic levels are unchanged.

Prepared cases also received both declared rotor integrations during searches:

```bash
i048_run early-rotors "$I048_PYTHON" "$I048_PROBE/early_rotors.py"
```

This producer uses core 12 and the per-case rotor lock. It finishes its current
case when all searches complete; the bulk rotor stage waits on its process
lock before launching eight lanes. These scheduling changes spread the measured
rotor cost through the search period without changing either sampling count.

After searches and the early producer finish, the bulk rotor CLI distributes
individual missing basin/count pairs over those same eight single-core lanes.
`rotor_worker.py --seal` first closes symmetry serially and records the full
selection and input hashes. Workers call the unchanged I043 integration with
`only_indices=[index]`, retaining every basin center, fixed seed and proposal.
Process-local adapters make symmetry closure and the shared selection write
read-only during each worker. Different workers write different quadrature
files. `pipeline/bulk_rotors_execution.json` records their selection hashes,
counts and lanes; the Verifier checks these against the actual producing jobs.
Each fixed CPU lane takes another pending basin/count pair after finishing its
current one. Assignment is recorded before execution, under a supervisor lock;
the per-basin seed and scientific result do not depend on assignment order.
The canonical whole-case CLIs then verify the cached pools and record their
positive receipts. Early integrations retain their original scheduling.

During the final recovery period, the idle early producer was reloaded with
`early_rotors.py --basin-parallel` to use this same tested basin scheduler for
subsequently prepared cases, sharing the original eight physical cores with
the searches. The original execution record and streams are retained; the
reload receipt verifies that the old supervisor had no scientific children.
Both sample counts, the full partition and per-basin seeds remain fixed.
The bulk stage still waits on the early producer's global process lock.

The verifier replays producing CLIs from saved pools, independently reconstructs
every importance draw/density/basin mask, reconstructs projected modes from
the source Hessians, audits reflection and proper-permutation orbit closure
and energy logs, checks both gas and internal H/S identities, and reproduces
all report tables. It freshly calculates the pinned baseline, two representative
PBE energies, and representative pentamer rotor energies. It does not claim
independent fresh stochastic search convergence or a fresh DFT calculation for
every saved conformer. `monitor.py` records observed calculation thread affinity
and aggregate RSS; its snapshots are sampled resource evidence.

Production is bounded by the recorded preparation time plus 48 hours. The
author's replay shares that deadline. A later manager replay of completed
evidence gets a two-hour verification deadline so it remains runnable after
production has ended; missing production evidence still fails the full audit.

An original eight-hour per-search wrapper was shorter than the dispatch's total
budget. The progressing `ps5_00110` search received this guardian command:

```bash
i048_run guard-ps5_00110 "$I048_PYTHON" "$I048_PROBE/hold_timeout.py" ps5_00110 2090743
```

The guardian verifies the owned wrapper's command and holds only that wrapper;
the original CREST/xTB child continues unchanged. It resumes the wrapper only
after successful scientific-child evidence, or before the original production
deadline. GNU timeout then retains status 124 despite the already-completed
child's exit 0, as reproduced on a controlled six-second sleep process.
`timeout_hold.json` records both statuses, source hashes and stop/resume times.
The Verifier requires this evidence before accepting a guarded search. The two
idle preparation supervisors were reloaded with this distinction, retaining
their prior logs and `workflow_timeout_queue_restart.json`. Their new streams
are `pipeline-budgeted` and `supplemental-budgeted`. No scientific input or
original total deadline changed. Newly invoked searches use the remaining total
budget directly. Historical PIDs are documentation and must not be reused.
The full replay independently repeats the controlled protocol with
`hold_timeout.py --self-check`; its new evidence is retained under
`verification/timeout_protocol/`.

The bulk electronic CLI uses the original three 3/3/2-core lanes for individual
single points, so a final single case can use all three lanes. Its case locks
cover selection, direct points and the subsequent unchanged offset/image
materialization. `pipeline/bulk_electronic_execution.json` records the tasks
and their lanes; the source audit checks them against the frozen selections
and actual producing jobs. Cached energy files are preserved.

The original search supervisor assigns multiple cases to each fixed lane. A
held wrapper's preserved status can end that lane before its next case starts.
Once its current scientific child finishes, a remaining frozen case can use
the same lane through `run_series.py --stage search --species ps5_01010
--search-lane 0` (cores 0, 2, 4, 6) or the corresponding `--search-lane 1`
for cores 8, 10, 12, 14. Confirm that the case is not already running before
launching it. These selectors change allocation only; they retain the frozen
input, CREST settings and original total deadline.

If the original cap leaves cases unfinished, keep the frozen weights and
report only increments whose two lengths are fully covered. These commands
render and reproduce that explicitly incomplete evidence:

```bash
i048_run available-base "$I048_PYTHON" "$I048_PROBE/thermochemistry.py" --available-only --samples 1024 --no-roots
i048_run available-thermo "$I048_PYTHON" "$I048_PROBE/thermochemistry.py" --available-only --samples 2048
i048_run available-report "$I048_PYTHON" "$I048_PROBE/report.py" --available-only --refresh-cost --write-report
i048_run available-verify "$I048_PYTHON" "$I048_PROBE/verify_results.py" --available-only --replay
i048_run full-verify "$I048_PYTHON" "$I048_PROBE/verify_results.py" --replay
```

The available-only audit is labeled incomplete and does not satisfy the full
dispatch Verifier. The final command still requires all 25 species and fails
if any are missing. Use `I048_REPLAY_DEADLINE_UNIX=1791202600.173893` for the
author's replay, including a final failing coverage check after the cap;
this prevents accidentally using the later manager's verification allowance.
Fresh quantum checks belong inside the original cap. Report rendering and
cached numerical audits after it do not allocate more quantum time. New job
wrappers reserve their 30-second kill grace inside the same deadline.

Partial-source rows are also audited directly: completed CREST frame counts,
energy-window membership and stereochemistry, preserved wrapper receipts,
the frozen candidate reduction, minimum geometry hashes and Hessian frequency
units. `verification/partial_source_audit.json` records those reproduced counts.
This validates the report's partial search/minimum evidence without treating
it as completed thermochemistry or filling any missing atactic class.

The author also launched two operational entry points, each with pinned,
persisted streams: `deadline_guard.py` before the deadline, and
`finish_series.py` while producers ran. The guardian signals only this user's
probe calculation tree and reserves termination grace inside the cap. It must
not be launched for a later manager replay. The finalizer waits for completed
production or the original cap, runs the full Verifier even for incomplete
coverage, then requires a cached numerical audit and matching report before
committing the requested files locally. It never pushes. Its actual commands,
exit codes and final SHA are written to `verification/final_commands.json`
and `verification/final_summary.json` in scratch. These operational commands
belong to the explicitly authorized author dispatch; a later manager can run
the Verifier directly without making another commit.
