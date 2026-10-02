# I-038 incremental propensity maintenance

Code: `e182b22b531b371d503e67346820c04016123c39`,
`7ba88b6ba4b0a0e9a62a86d9b737ac8171a6300f`,
`b081dff5740fb5d5d1576b8de6fab3af26f69df3` on `i038-incremental-propensity`; no push.

## Design

`IsothermalSSA` defaults to incremental maintenance; `incremental=False` retains
uncached full propensity/MET reconstruction. Site-index revisions track exact
atom mappings, site positions, strand ownership/length, and component membership
and length. Only changed record propensities and MET site/pair contributions are
refreshed. Fixed-temperature rates and bounds are evaluated lazily, skipped for
ineligible records, and invalidated when temperature/configuration changes.
Canonical trees retain all zero leaves and the original binary64 addition and
draw order; changed leaves recompute ancestors. Sparse immutable reports retain
all canonical IDs. Active MET samplers preserve deterministic construction order.
Shared parsed graphs and graph-to-record dependencies avoid inactive candidate
construction. No chemistry, rates, RNG consumption, or Cython was changed.

## Reproduction

Python `/home/alon/anaconda3/envs/rmg_env/bin/python`; scratch/logs
`/home/alon/runs/i038/`. Build `PYTHONPATH=$PWD python setup.py build_ext --inplace`
passed. Database SHA `4a12d36fcdc193ede82c8d1ab5c1653495d445bc` remained read only.
The supplied 14,998-record artifact SHA256 is
`0883a1292e20708a17ee8dc7960cf5648b354e0c24ce3373dee83813cea85c9d`.
Its compiler fingerprint matches current sources. The scratch pytest plugin
`i038_inputs.py` supplies this authenticated artifact instead of recompiling
unchanged chemistry; assertions, matching, budgets, and oracles are unchanged.

Commands use this environment and append separate stdout/stderr logs with
`> >(tee -a RUN.stdout.log) 2> >(tee -a RUN.stderr.log >&2)`:

```bash
export I038_PYTHON=/home/alon/anaconda3/envs/rmg_env/bin/python
export PYTHONPATH=$PWD:/home/alon/runs/i038 PYTHONHASHSEED=0
export RMG_DATABASE_PATH=/home/alon/Code/RMG-database
export RMG_KMC_CACHE_ROOT=/home/alon/runs/i038/event-cache
export I038_PERIODIC_STACKS=0 TMPDIR=/home/alon/runs/i038/tmp
"$I038_PYTHON" -m pytest -p i038_inputs test/rmgpy/kmc -o addopts='' -q
# SSA partitions: rate oracle; other cases; full-catalogue ring and melt traces.
RMG_KMC_SLOW=1 "$I038_PYTHON" -m pytest -p i038_inputs test/rmgpy/kmc/ssaTest.py -o addopts='' -q -s -k real_k_act_is_pair_specific
RMG_KMC_SLOW=1 "$I038_PYTHON" -m pytest -p i038_inputs test/rmgpy/kmc/ssaTest.py -o addopts='' -q -s -k 'not real_k_act_is_pair_specific and not incremental_full_catalogue'
RMG_KMC_SLOW=1 "$I038_PYTHON" -m pytest -p i038_inputs test/rmgpy/kmc/ssaTest.py -o addopts='' -q -s -k 'incremental_full_catalogue and not 35000002'
RMG_KMC_SLOW=1 "$I038_PYTHON" -m pytest -p i038_inputs test/rmgpy/kmc/ssaTest.py -o addopts='' -q -s -k 'incremental_full_catalogue and 35000002'
RMG_KMC_SLOW=1 "$I038_PYTHON" -m pytest -p i038_inputs test/rmgpy/kmc/metTest.py -o addopts='' -q -s
RMG_KMC_SLOW=1 "$I038_PYTHON" -m pytest -p i038_inputs test/rmgpy/kmc/stateTest.py -o addopts='' -q -s
```

## Exact traces

Committed default tests compare events, hashes, binary64 times, RNG state, and
final counters; they cover three seeds, temperature changes, lazy rates, odd
reduction trees, report snapshots, mapping changes, and strand redistribution.
Final focused verification: **14 passed, 19 deselected, 1.17 s**; Ruff and
`git diff --check` passed.

Committed slow tests compare both engines step by step on all 14,998 records,
2,000 events per seed. Ring seeds 8013700/8013701 passed (**2 passed, 29 deselected,
1956.47 s**) on `e182b22b`; final code reproduced their complete trace digests.
Pristine comparison on `b081dff57`: **1 passed, 32 deselected, 3783.07 s**,
process wall 3797.38 s. All three seeds matched all 2,000 events and final counters.

| Seed | Initial state | SHA256 of event ID, state hash, and time.hex() per step |
|---|---|---|
| 8013700 | R1 ring, 700 K | `6d59b57955fba11bd350a859cca795d0a9be88488a38491090e228279f186381` |
| 8013701 | R1 ring, 700 K | `0f001b1ddd9bc4b599fad7322675de672edecd4669f22a716c95c250242b2e23` |
| 35000002 | Pristine melt, 800 K | `ddefbc3b9cb24a1344a0b0d01ba9331531fba72312d2630411b68e8878f1e7f1` |

The pristine case uses the named existing test's five chains (3/5/3/5/3 units),
1.05 g/cm³, volume `3.12948355531e-27 m³`. Logs: `traces-ring-final.*`,
`traces-melt-b081.*`, `replay-final.*`; per-event final ring traces: `final-traces/`.

## Benchmark

The supplied benchmark's repository pointer now names this worktree. Its scratch
copy changes implementation selection/output only, preserving inputs, seed,
20 warmups, 2,000 firings, 200 profile firings, and timing boundaries. Before-head
loads the original SSA snapshot; before-base uses archived code. After-base uses
current code on the unchanged 314-record input. Hardcoded JSON commit labels are
historical artifact labels. Before-head uses `0bbfe8ea3`; before-base uses
the supplied archive. After uses `b081dff57`.

```bash
I038_BEFORE=1 I038_COST_ROOT=/home/alon/runs/i038/cost-before "$I038_PYTHON" /home/alon/runs/i038/benchmark_common_ssa.py base --events 2000 --profile-events 200
I038_BEFORE=1 I038_COST_ROOT=/home/alon/runs/i038/cost-before "$I038_PYTHON" /home/alon/runs/i038/benchmark_common_ssa.py head --events 2000 --profile-events 200
I038_CURRENT_BASE=1 I038_COST_ROOT=/home/alon/runs/i038/cost-final "$I038_PYTHON" /home/alon/runs/i038/benchmark_common_ssa.py base --events 2000 --profile-events 200
I038_CURRENT_BASE=1 I038_COST_ROOT=/home/alon/runs/i038/cost-final "$I038_PYTHON" /home/alon/runs/i038/benchmark_common_ssa.py head --events 2000 --profile-events 200
```

| Artifact | Records | Before ms/event | After ms/event |
|---|---:|---:|---:|
| base | 314 | 10.2269 | 3.4768 |
| head | 14,998 | 519.7094 | 3.5108 |

Final head/base **1.010×**, meeting the ~2× target. Reproduced head speedup
**148.03×**; shared-host contention varied, so absolute timings are descriptive.
Both final warm profiles have **zero** record-rate/forward-rate/table interpolation
calls; before-base had 62,427 record-rate calls, before-head 1,672,827. Final head
profile: 1.466 s total, 1.324 s in indexed apply, 1.098 s in product apply including
0.733 s deepcopy. JSON/profiles: `cost-before/`, `cost-final/`. All four benchmark
event/state digests match `b0c033ad600a25e44b4480127b1e902d5bf34bb348bba94167cae135c3481c7f`.
The benchmark audits all records but selects six common ring records; mixed
production throughput requires separate measurement. Its digest excludes time;
the full-catalogue replay checks time explicitly.

## Suites and state oracle

Final complete default suite: **166 passed, 48 skipped, 458.45 s**
(`default-b081.*`), process wall 472.82 s, including the unchanged 10,000-step
state oracle.
Slow MET: **40 passed, 1 skipped, 649.42 s**, process wall 696.90 s. The skip is
the pre-existing deferred kernel reference. Slow SSA was partitioned into its
rate oracle, the three full-catalogue comparisons, and remaining cases. Rate
oracle retry: **1 passed, 31 deselected, 773.39 s**, process wall 787.22 s;
remaining cases: **28 passed, 4 deselected, 8029.84 s**, process wall
8042.87 s (`slow-ssa-rest.*`, loaded on `7ba88b6ba`). The final default/focused
runs cover the subsequent strand-length fix. Including the three trace cases, rate oracle,
and new strand-length regression, all **33** slow SSA cases passed across
partitions; this is not a one-invocation count. The 10,000-event pristine case
fired 306 non-inverse additions, retained ledger C=152/H=162, and passed.
Initial instrumented rate and pristine comparison runs exited 139; the latter
matched through 1,000 events.
Native crash cause is unconfirmed; retries disable periodic traceback dumps.

Slow state suite, including the unchanged million-step test: **23 passed in
11632.76 s (3 h 13 min 52.76 s)**; process wall **11642.66 s** (`slow-state.*`).
It completed exactly **1,000,000 steps**, across 24 chains and all five families,
with 13 cuts, 499950 component splits, and 499939 merges. The unchanged
trajectory directly calls independent graph discovery and `KMCState.apply`,
never SSA. A cProfile run of the actual first worker's first 2,000 steps (seed
8675309, Disproportionation) took **465.835 s**: candidate discovery 278.597 s,
including NetworkX matching 254.650 s, versus product apply 30.656 s. This
reproduces an independent oracle bottleneck outside propensity maintenance.
Profile: `state-oracle-long.profile`, logs `state-profile-long.*`. The oracle
and million-step budget were not weakened.
