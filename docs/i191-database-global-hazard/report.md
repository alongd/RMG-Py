# i191 — Is the process-global database a reachable hazard in production, or only in tests?

**Verdict: it is reachable in production.** All three questions answer *yes*, with one important
qualification on the third. The hazard is not a test artifact.

Branch `i191-database-global-hazard`, worktree
`/home/alon/Code/RMG-Py-i191-database-global-hazard`, base `546d1c727`.
Every measurement below was made against a `make build` of this worktree (exit 0; the loaded
`rmgpy.molecule.molecule` is
`/home/alon/Code/RMG-Py-i191-database-global-hazard/rmgpy/molecule/molecule.cpython-39-x86_64-linux-gnu.so`,
not a symlink to another tree).

---

## 0. The resolved database, and the trap firing once

```
database.directory = /home/alon/Code/RMG-database-plasma/input
resolved            = /home/alon/Code/RMG-database-plasma/input
cwd                 = /home/alon/Code/RMG-Py-i191-database-global-hazard
```

(`stdout.log`.) `./rmgrc` was created in this worktree from `rmgrc.template` (untracked by
design; `git check-ignore` confirms `.gitignore:97`). Every probe script prints this line first
and asserts on it.

**The cwd trap fired once and was caught by that assert.** The first run of
`q1b_dynamic_two_constructions.py` was launched from a different working directory and resolved
`database.directory = /home/alon/Code/RMG-database/input` — a *different database*. The script's
own assertion killed it (`q1b_dynamic_two_constructions.stderr.log:5-8`,
`AssertionError: wrong database.directory`), and it was re-run from the worktree root. Both runs
are in the `.stdout.log` because it is opened append; the first contributes only that one line.
Had the assert not been there, that run would have returned green about the wrong database.

---

## 1. Can a single process construct `RMGDatabase` more than once? — **Yes**

### The shape, confirmed at source

`rmgpy/data/rmg.py:48` holds `database = None` at module level. `RMGDatabase.__init__`
(`rmgpy/data/rmg.py:59-67`) ends with `global database; database = self`, preceded only by a
`logging.warning` when one already exists. Nothing anywhere resets it — `grep` finds no
assignment to `rmgpy.data.rmg.database` outside `__init__` and `test/conftest.py`.

### Production call paths

| Site | Constructs | Condition |
|---|---|---|
| `rmgpy/rmg/main.py:517` | every call to `RMG.load_database()` | unconditional |
| `rmgpy/tools/simulate.py:83` and `:88` | via `rmg.load_database()` | **inside** the per-reaction-system loop opened at `simulate.py:53` |
| `rmgpy/tools/fluxdiagram.py:604` | via `rmg.load_database()` | **inside** the comparable loop opened at `fluxdiagram.py:581` |
| `rmgpy/tools/uncertainty.py:478` | `Uncertainty.load_database()` | called from the five `ipython/*.ipynb` notebooks; `rmgpy/rmg/main.py:1727-1728` deliberately assigns `uncertainty.database = self.database` instead |
| `arkane/input.py:584` | `load_necessary_databases()` | first thing `load_input_file` does (`arkane/input.py:654`, before the `exec` at `:656`) |
| `arkane/input.py:111` | `get_db() or RMGDatabase()` | see "dead branch" below |

The decisive one is `rmgpy/tools/simulate.py`. `rmgpy/tools/loader.py`'s `load_rmg_job` does
**not** load a database, so the first construction happens *inside* the loop; with N reaction
systems and either a diffusion-limited `LiquidReactor` (line 83) or an `uncertainty()` block
(line 88), `rmg.load_database()` runs N times and constructs N databases. `scripts/simulate.py`
and `scripts/generateFluxDiagram.py` are both in `setup.py:140-155`'s `scripts=` list, i.e.
installed executables on `PATH`.

### Dynamic demonstration — **in-process** (necessarily)

`q1b_dynamic_two_constructions.py`, exit 0. One process, unmodified production code, calling
`RMG.load_database()` twice the way `simulate.py` does:

```
id(rmgpy.data.rmg.database) after call #1 = 125004333387248
id(rmgpy.data.rmg.database) after call #2 = 124999935625008
RESULT: identity changed across the two calls in ONE process: True
--- captured logging.warning output ---
An instance of RMGDatabase already exists. Re-initializing it.
```

The load was narrowed (one kinetics family, empty library lists) for speed; the narrowing bears
on load time only, not on whether the global is rebound.

### Two subsidiary findings

**`get_db('transport')` does not raise when the component is missing.** `q1c_transport_statmech_none.py`,
exit 0 — when the global `database` object exists but `.transport`/`.statmech` are `None`,
`get_db` returns `None` silently. `arkane/input.py:580-581` guards
`load_necessary_databases()` with `try: get_db('transport'); get_db('statmech') except
DatabaseError:`, so in that state the guard no-ops and transport/statmech are **never loaded**.
Reachability: in a standalone Arkane run this cannot happen, because `load_necessary_databases`
is the first thing `load_input_file` does. It is reachable in any process that loaded a database
some other way first — which is what the test suite does. **Latent, not currently live in
Arkane.**

**`arkane/input.py:111`'s `get_db() or RMGDatabase()` is a dead branch.** `get_db()` never
returns a falsy value — it raises `DatabaseError`. Measured in `q1a_arkane_input_line111_guard.py`
(exit 0): called on a fresh process, the `database(...)` directive raises rather than falling
back. It does not crash in real Arkane only because `load_necessary_databases()` has already
constructed one at line 584. The guard reads as defensive and is not.

---

## 2. Does the second construction change what an already-created species resolves to? — **Yes.** This is the teeth.

Four **separate processes**, so the global could not fool the instruments. Thermo-only loads;
methane; `primaryThermoLibrary` present (DB_A) vs absent (DB_B).

| Process | What it does | H298 (kJ/mol) | comment |
|---|---|---|---|
| P1 | DB_A only, create, resolve | −74.53377599999996 | `Thermo library: primaryThermoLibrary` |
| P2 | DB_B only, create, resolve | −74.89359999999998 | `Thermo group additivity estimation: group(Cs-HHHH)` |
| P4 (control) | DB_A only, create, resolve immediately | −74.53377599999996 | `Thermo library: primaryThermoLibrary` — matches P1, harness validated |
| **P3 (experiment)** | DB_A, create species, **do not resolve**, construct DB_B, *then* resolve | **−74.89359999999998** | **`Thermo group additivity estimation: group(Cs-HHHH)`** |

P3 also prints which object answered:

```
get_db('thermo') is db_a.thermo: False
get_db('thermo') is db_b.thermo: True
```

**A species created while one database was live resolves against whatever database is live at
resolution time, not the one that was live at creation time. Silently, with a different number.**

Kinetics gives the same answer with a louder failure mode. `q2_kinetics_teeth.py` (exit 0)
reproduces `rmgpy/rmg/model.py:1073`'s exact expression
`get_db('kinetics').families[reaction.family]` for a reaction whose `.family` was set under a
database that loaded `H_Abstraction`, after a second database that did not:

```
Lookup RAISED KeyError: 'H_Abstraction'
```

### Content divergence vs identity divergence

`simulate.py`'s two loads use the same `rmg` attributes both times, **except** that line 87 sets
`rmg.verbose_comments = True` before the line-88 reload. That flag reaches
`rmgpy/rmg/main.py:610` — `family.fill_rules_by_averaging_up(verbose=self.verbose_comments)` —
and from there `rmgpy/data/kinetics/family.py:1360-1371`'s rule averaging. So the second database
in that path differs from the first in *content*, not only identity. **This specific claim rests
on static source tracing, not on a running measurement** — see §6.

Identity divergence alone: a `get_db(` sweep across `rmgpy/` and `arkane/` — 27 call sites in 12
files, excluding the definition — found **no** comparison of database-derived objects by `is`,
`id()`, or identity-keyed container membership. `rmgpy/rmg/main.py:1734-1736` mutates the existing object
(`self.database.kinetics.load_families(...)`) rather than constructing a new one, so that
mutation is not lost. This negative is grep-based, not exhaustive.

### A third mechanism, found while tracing: `arkane/explorer.py` mutates the global in place

The ticket named this path and it turns out to carry its own teeth, with **no second construction
involved at all**. `arkane/explorer.py:111-118` does:

```python
thermo_database.libraries['thermojobs'] = thermo_library
thermo_database.library_order.insert(0, 'thermojobs')
kinetics_database.libraries['kineticsjobs'] = kinetics_library
kinetics_database.library_order.insert(0, ('kineticsjobs', 'Reaction Library'))
```

on the objects returned by `get_db('kinetics')` / `get_db('thermo')`. Nothing undoes it, and
`arkane/main.py:293-296` runs that block once per `ExplorerJob` in the input file, in one process.

`probe_explorer_global_mutation.py` (exit 0, **in-process** by nature — the hazard is in-process)
performs those two lines twice against a real loaded thermo database:

```
-- baseline, before any explorer job --
  H298 = -74.534 kJ/mol | comment = Thermo library: primaryThermoLibrary
-- after explorer job #1 --
  library_order = ['thermojobs', 'primaryThermoLibrary', 'GRI-Mech3.0']
  H298 = -74.605 kJ/mol | comment = Thermo library: thermojobs
-- after explorer job #2 --
  library_order = ['thermojobs', 'thermojobs', 'primaryThermoLibrary', 'GRI-Mech3.0']
  'thermojobs' appears 2 time(s) in library_order
```

Read this carefully, because half of it is *by design*: for one explorer job, preferring the
job's own computed thermo is the point. The defect is that the change is unscoped — it outlives
the job that made it. With two `explorer()` blocks in one input file, job #2's species resolve
through job #1's leftover library until job #2 overwrites the key, and `library_order` accumulates
a duplicate entry per job.

Honest limits on this one: the stand-in for the explorer's generated `thermojobs` library is a
real already-loaded library (GRI-Mech3.0) re-registered under that key — the *mechanism* is
verbatim, the *content* is a stand-in. Both shipped examples
(`examples/arkane/explorer/{methoxy,methyl+formaldehyde}/input.py`) contain exactly one
`explorer()` block, so the accumulation case is legal but unexercised in-repo. I did not run a
real two-explorer Arkane input end to end.

---

## 3. Does `fork` carry the database into workers, and can a worker's mutation be observed? — **Carried in: yes. Parent and siblings: no. The same worker's later tasks: yes.**

This is the question where reasoning would have got it half right. All of it is measured.

The production parallelism is `rmgpy/rmg/react.py`: `Pool(processes=procnum)` at line 69,
`p.map(_react_species_star, ...)` at line 70, and the worker reaches the database through
`get_db('kinetics')` at line 94 — it is never passed the database as an argument.

**Start method** (`probe_q3_1_start_method.py`, exit 0): a plain non-pytest interpreter here gets
`multiprocessing.get_start_method() == 'fork'` as the *platform default*, not because RMG sets it.
`test/conftest.py` forces `'fork'` for the suite, but on this Linux environment that forcing is
redundant with the default — production gets `fork` regardless.

**Carried in** (`probe_q3_2_carried_into_workers.py`, exit 0; parent **in-process**, workers
**separate-process**): a database loaded in the parent before `Pool()` is visible in every worker
with matching content (`library_keys`, entry counts) and the same `id()` value —
`137186340713712` in parent pid 7 and in worker pids 52/53/54. The matching `id()` is the same
virtual address under copy-on-write, not shared memory.

**Worker → parent: not observable** (`probe_q3_3_worker_mutation_to_parent.py`, exit 0). Two
distinct mutations tested, both negative at the parent after `close()`/`join()`:

```
(a) RESULT in_place_worker_mutation_visible_to_parent = False
(b) RESULT rebind_in_worker_visible_to_parent        = False
```

**Worker → sibling: not observable. Worker → its own later tasks: fully observable**
(`probe_q3_4_sibling_and_persistence.py`, exit 0). 40 tasks over 4 workers,
`maxtasksperchild=None` as production uses:

```
RESULT distinct_worker_pids   = [53, 54, 55, 56]
RESULT any_cross_pid_leak     = False
RESULT any_same_pid_persistence = True
```

with 95 concrete `(task_idx, pid, earlier_task_idx_same_pid)` examples — e.g. pid 54 handled tasks
`[1, 6, 7, 8, 10, 11, 13, 16, 19, 23, 26, 29, 32, 35, 39]` and by task 39 carried the sentinels
from all fourteen earlier tasks it had personally handled. The rebind experiment shows the same
thing directly: on pid 58, task 3's `old_id(database)` is the id task 0 rebound *to*.

This is the finding that matters. A `Pool` worker persists across tasks, so any database mutation
made while handling task 1 is still in that process when it handles task 40 — and **which tasks
land on which worker is not deterministic**. A mutating worker would therefore make results depend
on task-to-worker assignment, while being invisible to the parent that collects them.

**Does the real `react.py` path mutate anything?** Traced and measured, and the answer is no on
this evidence. `generate_reactions_from_families` (`rmgpy/data/kinetics/database.py:509-553`) →
`family.generate_reactions` (`rmgpy/data/kinetics/family.py:1845`), whose docstring states it does
not estimate kinetics. `KineticsRules.estimate_kinetics` (`rmgpy/data/kinetics/rules.py:361`) is
reachable only via `family.get_kinetics`, called solely from `rmgpy/rmg/model.py:1095`
(`CoreEdgeReactionModel.generate_kinetics`) — which runs in the **parent**, after `p.map` returns.
The one in-chain write, `self.forbidden = ForbiddenStructures()` in `add_reverse_attribute`
(`family.py:1940`), is restored by `try/finally` and executes only on a rare zero-reverse-reactions
branch. `probe_q3_5_production_path_mutation.py` (exit 0) drove the real call with real families
(`R_Recombination`, `H_Abstraction`) and real molecules over 24 tasks / 4 workers:

```
RESULT any_within_task_change              = False
RESULT any_same_pid_drift                  = False
RESULT parent_fingerprint_changed_after_pool = False
RESULT total_reactions_generated_across_all_tasks = 72   (probe_is_live = True)
```

So: the leak *channel* is real and wide open; today's `react.py` payload does not appear to put
anything into it. That is a property of the current call chain, not a property of the design.

---

## 4. The docstring correction

`get_db`'s docstring said *"First, the module level is queried. If this variable is empty, the
broadcasted variables are queried."* The body has no such fallback — it raises. The phrase is a
fossil of a `scoop`-based parallelisation (a `scoop` target still exists in the `Makefile`) whose
broadcast mechanism is gone.

`rmgpy/data/rmg.py:247-262` now describes what the function does: the module-level variable is the
only source; there is no fallback; an unset variable raises `DatabaseError`; a recognized name may
return `None` for a component that was never loaded, which is *not* an error here, so callers must
handle it; an unrecognized name raises `ValueError`.

**A project-wide sweep for the same fossil found one more, and it is corrected too.**
`rmgpy/rmg/input.py:2548-2555`, `get_input`, carried the identical sentence over an identically
fallback-free body. Fixing only the site the ticket named is how a partial fix gets read as total.
A case-insensitive `grep` for `broadcast` across the tree now returns only two hits, both in
`rmgpy/rmg/main.py` (lines 2300 and 2446) and both about numpy array broadcasting — unrelated.

`rmgpy/rmg/input.py` is outside the owner-gated `rmgpy/data/` tree. If the owner wants that second
correction reverted, it is a self-contained one-hunk revert.

---

## 5. Test counts

```
python -m pytest test/rmgpy/data/ -p no:cacheprovider --no-cov -q
```

| | result |
|---|---|
| before (base tree) | **388 passed, 7 skipped** in 17.99 s |
| after (both docstring corrections) | **388 passed, 7 skipped** in 16.03 s |

Logs: `tests-before.stdout.log`, `tests-after.stdout.log`. The only source changes in the worktree
are the two docstrings (`git diff --stat`: `rmgpy/data/rmg.py`, `rmgpy/rmg/input.py`).

---

## 6. What this investigation may and may not be claimed to show

**It may be claimed that:**

- The process-global database is reachable from production code, not only from tests. One process
  constructing `RMGDatabase` more than once happens on an installed console script's ordinary
  loop (`scripts/simulate.py`, `scripts/generateFluxDiagram.py`), not only under pytest.
- The hazard has teeth. An already-created species that has not yet resolved its thermo silently
  resolves against the *later* database and returns a different number with a different source
  comment (P1/P3/P4, three separate processes, harness validated by the control). The kinetics
  analogue fails loudly with `KeyError` instead.
- `fork` carries a loaded database into workers with identical content; worker mutations do **not**
  reach the parent or a sibling; they **do** persist into later tasks handled by the same worker
  process. All four of those are measured, with pids.
- `arkane/explorer.py`'s in-place mutation of the global thermo/kinetics databases changes what a
  subsequently-resolved species gets, is never undone, and accumulates duplicate `library_order`
  entries across explorer jobs.
- The `get_db` docstring promised a fallback mechanism that does not exist, and so did `get_input`.

**It may NOT be claimed that:**

- **That any of this has actually corrupted a real result.** A `grep -rIl` over `/home/alon/Code`
  for `"An instance of RMGDatabase already exists"`, restricted to `*.log`/`*.txt`/`*.out`, found
  six files: **five pytest runs and one ad-hoc probe script's stderr**
  (`plasma-pm2/archive/i027-i032-audit-scripts/i034_alpha_regen_probe_stderr.log`). **None is a
  production `rmg.py` or `Arkane.py` run.** That search covers only those three extensions and only
  this machine. The reachability is demonstrated; a field incident is not.
- **That the `verbose_comments` content-divergence claim was measured.** It is static tracing only
  (`main.py:610` → `family.py:1360-1371`). No script diffs a rule's comment across
  `verbose_comments=True`/`False` loads. Everything else in §2 has a running artifact.
- **That the identity-divergence negative is exhaustive.** "No production code compares
  database-derived objects by identity" is a `grep` sweep over `get_db(` call sites. A comparison
  behind an indirection would not have been caught.
- **That `react.py` is safe under parallelism in general.** The measurement covers
  `R_Recombination` and `H_Abstraction` with a minimal two-family kinetics load. A full-database
  run was not attempted, other families were not exercised, and the one candidate write site
  (`family.forbidden` in `add_reverse_attribute`) sits behind an error branch that was not forced.
  The claim is "no mutation observed on this path with these families", not "no mutation possible".
- **That the two-explorer accumulation happens in practice.** It is legal and unexercised in-repo;
  no real two-explorer Arkane input was run, and the `thermojobs` library in the probe is a
  real-library stand-in for a generated one.
- **That the `arkane/input.py` transport/statmech no-op (§1) is live.** In a standalone Arkane run
  it cannot be reached, because `load_necessary_databases()` runs first. It is latent.

**Not reached at all:**

- No end-to-end run of `scripts/simulate.py` on a real multi-reaction-system job with an
  `uncertainty()` block. The double construction there is proved by the loop structure plus an
  in-process reproduction of the two `load_database()` calls, not by driving the script.
- The IPython notebook path (`Uncertainty.load_database`, `rmgpy/tools/uncertainty.py:478`) is
  identified statically only. A Jupyter kernel is long-lived, so re-running one cell constructs a
  second database — plausible and unmeasured.
- No measurement of whether a *second* construction's cost (a full database reload per reaction
  system) has ever been noticed as a performance problem.

---

## 7. Recommendation, and where it stops

The evidence supports a design change: the resolution point (`get_db`) should not be able to
silently answer from a database other than the one a caller's objects were built against. That is a
dependency-injection question, and per the ticket's non-goals **it stops here as a conclusion, not
as work.** The blast radius — 27 `get_db(` call sites in 12 files across `rmgpy/` and `arkane/`, at
least the two in `rmgpy/rmg/react.py` executing inside `multiprocessing` workers — makes it the
owner's call.

Two much smaller, independent items fall out of the investigation and are *also* left unbuilt,
recorded here so they are not lost:

1. `rmgpy/tools/simulate.py:83,88` and `rmgpy/tools/fluxdiagram.py:604` reload the entire database
   once per reaction system. Hoisting the load out of the loop fixes both a correctness hazard and
   a large, silent performance cost, and is a self-contained change.
2. `arkane/explorer.py:111-118` should scope its `thermojobs`/`kineticsjobs` insertion to the job —
   or at least not push a duplicate `library_order` entry per job.

---

## Artifacts

All under `docs/i191-database-global-hazard/`. Every script prints the resolved
`database.directory` first and asserts on it; every run has a paired `.stdout.log` / `.stderr.log`.

| Script | Question | Measurement |
|---|---|---|
| `q1a_arkane_input_line111_guard.py` | 1 | separate-process |
| `q1b_dynamic_two_constructions.py` | 1 | in-process (necessarily) |
| `q1c_transport_statmech_none.py` | 1 | in-process |
| `q2_p1_thermo_dbA.py` | 2 | separate-process |
| `q2_p2_thermo_dbB.py` | 2 | separate-process |
| `q2_p3_thermo_second_construction.py` | 2 | separate-process (two constructions inside it) |
| `q2_p4_control.py` | 2 | separate-process |
| `q2_kinetics_teeth.py` | 2 | separate-process |
| `probe_explorer_global_mutation.py` | 2 | in-process (the hazard is in-process) |
| `probe_q3_1_start_method.py` | 3 | in-process |
| `probe_q3_2_carried_into_workers.py` | 3 | parent in-process, workers separate-process |
| `probe_q3_3_worker_mutation_to_parent.py` | 3 | mutation separate-process, check in-process |
| `probe_q3_4_sibling_and_persistence.py` | 3 | observations separate-process, aggregation in-process |
| `probe_q3_5_production_path_mutation.py` | 3 | worker fingerprints separate-process, diffing in-process |

`stdout.log` / `stderr.log` at the top of that directory carry the settings resolution and the
docstring verification; `build.stdout.log` / `build.stderr.log` carry the `make build`.
