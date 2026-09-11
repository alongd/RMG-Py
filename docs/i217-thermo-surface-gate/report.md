# i217 — RMGDatabase.load() gates thermo loading on `surface`

## Resolved database directory

`rmgrc` (git-ignored, copied from `rmgrc.template`) pins `database.directory = ../RMG-database-plasma/input`,
resolving to `/home/alon/Code/RMG-database-plasma/input` (confirmed to exist), printed via
`python -c "import rmgpy; print(rmgpy.settings['database.directory'])"` before any other measurement.

## The defect, reproduced

On the unmodified tree at `546d1c727`, in a fresh process:

```python
db = RMGDatabase()
db.load(path=..., surface=False, testing=True, depository=False, solvation=False, ...)
assert db.thermo is None   # true; nothing raised
```

`RMGDatabase.load()` called `self.load_thermo(...)` only inside `if surface:`
(`rmgpy/data/rmg.py:111-112`), and `RMGDatabase.__init__` sets `self.thermo = None` (line 60) — the
sentinel silently survives `load()` with `surface=False`.

**Downstream symptom and its distance from the cause.** Mirroring the real call site
`rmgpy/rmg/main.py:1980` (`list(self.database.thermo.libraries.keys())`) against a `db` built the
same way produces:

```
AttributeError: 'NoneType' object has no attribute 'libraries'
```

Distance: the raise happens in `rmgpy/rmg/main.py`, a different module from the one that set the
sentinel (`rmgpy/data/rmg.py`). The message names neither `thermo` nor `surface` — it only reports
`NoneType has no attribute 'libraries'`. Reconstructing the cause requires walking backward from
this `AttributeError` through `RMGDatabase.load()` to the `if surface:` gate — three files apart,
no shared vocabulary between the exception message and the trigger.

Logs: `docs/i217-thermo-surface-gate/repro_stdout.log` / `repro_stderr.log` (the `None`, no raise);
`repro2_stdout.log` / `repro2_stderr.log` (the downstream `AttributeError` and its traceback).

## Git archaeology

`git log -L112,113:rmgpy/data/rmg.py` shows three touching commits:

- `2a7998c7a9` introduced `if surface: self.load_surface(...)`. At that point `load_thermo(...)`
  was called unconditionally as the first statement of `load()`.
- `191253b486` ("Move the metal database to be attached to the thermo database", 2021-01-20)
  folded the metal database into `ThermoDatabase` and **replaced** the body of the `if surface:`
  branch with `self.load_thermo(...)`, deleting the unconditional call. This is the regression:
  thermo inherited a gate that existed for `load_surface`, not for itself. Confirmed via
  `git show 191253b486^:rmgpy/data/rmg.py`, which shows `self.load_thermo(...)` unconditional.
- `4b934ce39f` only threaded `adsorption_groups` through; it did not touch the gating.

The fix restores pre-2021 behaviour for thermo while leaving the `surface` argument's
adsorption-groups role, added by the same 2021 commit, intact.

## Caller search

Search for `surface=False` (keyword form) across the whole tree: exactly four sites, all tests
against the lightweight `testing_database`, matching the brief's preliminary list exactly:

- `test/rmgpy/data/kinetics/kineticsTest.py:69`
- `test/rmgpy/data/kinetics/familyTest.py:830` (`TestTreeGeneration`)
- `test/rmgpy/data/kinetics/familyTest.py:1004` (`TestGenerateReactions`)
- `test/rmgpy/tools/isotopesTest.py:70`

**Positional-argument search:** grepped every `*.load(`/`*.database.load(` call site in `rmgpy/`,
`arkane/`, `scripts/`, `test/`. No call site passes `surface` positionally — every other caller
either omits it (default `True`) or uses the keyword form. The four-site list is complete.

None of the four reads `database.thermo`, `database.thermo.groups`, or `database.thermo.surface`
after loading. `familyTest.py` additionally constructs its own separate `ThermoDatabase()` instance
(`cls.thermoDatabase`, line 833) rather than relying on `RMGDatabase.thermo`, so it is doubly
unaffected. No genuine correctness dependency on the skip was found.

**Correction to the brief's cost claim.** The brief states "two of those sites pass no
`thermo_libraries`." Measured: only **one** does.

- `kineticsTest.py:69` passes `thermo_libraries=["primaryThermoLibrary"]` — unaffected.
- `familyTest.py:830` and `:1004` pass `thermo_libraries=[]` explicitly — an empty list, not
  `None`. `ThermoDatabase.load_libraries` treats `None` as "load everything" and `[]` (or any
  non-`None` list) as "load exactly these," so an explicit `[]` loads zero libraries both before
  and after the fix. Unaffected.
- `isotopesTest.py:64-71` passes no `thermo_libraries` keyword at all -> default `None` -> after
  the fix this site now loads every thermo library under `testing_database/thermo/libraries`.

So the real cost site is `isotopesTest.py` alone, not two sites. Measured wall-clock,
`python -m pytest test/rmgpy/tools/isotopesTest.py -p no:cacheprovider --no-cov -q` under
`/usr/bin/time -v`, one run each (not averaged), each config in its own process:

| | before (baseline, `546d1c727`) | after (fixed) |
|---|---|---|
| wall clock | 4.08 s | 5.33 s |
| exit | 0 | 0 |

+~1.25 s (~30%) on this one file. Per the brief's non-goal, this slowdown is reported, not
optimised away.

## The fix

`rmgpy/data/rmg.py`, `RMGDatabase.load()`: moved the `self.load_thermo(...)` call out of the
`if surface:` branch so thermo loads unconditionally; `surface` is still forwarded as
`load_thermo`'s (and, through it, `ThermoDatabase.load`'s) adsorption-groups/metal-database switch.
Diff is exactly this call site, `if solvation:` and the `if not testing:` blocks untouched:

```diff
-        if surface:
-            self.load_thermo(os.path.join(path, 'thermo'), thermo_libraries, depository, surface, adsorption_groups)
+        self.load_thermo(os.path.join(path, 'thermo'), thermo_libraries, depository, surface, adsorption_groups)
```

Confirmed via `git diff --stat rmgpy/data/rmg.py`: `1 file changed, 1 insertion(+), 2 deletions(-)`.

## Regression test

`test/rmgpy/data/rmgTest.py` (new), class `RMGDatabaseLoadSurfaceTest`, two tests in one process
(deliberately ordered `surface=False` then `surface=True`, documented in the module docstring —
each `RMGDatabase()` instance owns its own `self.thermo`; the ordering is explicit, not
incidental). Both assert on values:

- `surface=False`: `db.thermo is not None`; `'group' in db.thermo.groups` with real entries;
  `'primaryThermoLibrary' in db.thermo.libraries`; `db.thermo.surface == {}`.
- `surface=True`: same thermo assertions, plus `'metal' in db.thermo.surface`, not `None`.

**Contradiction of the brief, discovered empirically and acted on per the brief's own
"measurement wins" instruction.** The brief's steps 3 and 5 predicted `db.thermo.groups` would
show the adsorption-groups entry present vs. absent depending on `surface`. `ThermoDatabase.
load_groups()` unconditionally loads the `adsorptionPt111`/`adsorptionLi` group categories into
`self.groups`, regardless of the `surface` argument — confirmed by reading `rmgpy/data/thermo.py`
and by two standalone `ThermoDatabase().load(..., surface=False)` / `load(..., surface=True)`
scripts (`docs/i217-thermo-surface-gate/groups_check_stdout.log`, `groups_check2_stdout.log`):
`db.thermo.groups` keys are identical in both cases. The value actually gated by `surface` is
`ThermoDatabase.load_surface()`, populating `self.surface` (`{}` vs. `{'metal': MetalDatabase(...)}`).
The test asserts on `db.thermo.surface`, not on presence/absence of a `db.thermo.groups` key.

RED on unmodified tree (`docs/i217-thermo-surface-gate/red_stdout.log`): `1 failed, 1 passed` — the
`surface=False` test fails with `AssertionError: RMGDatabase.load(surface=False) left db.thermo as
None ...`; the `surface=True` test passes (today's default path, unaffected).

GREEN after the fix (`docs/i217-thermo-surface-gate/green_stdout.log`): `2 passed`.

## Controls, each in a separate process

**Control A — `surface=False`, fixed tree** (`control_A_false_stdout.log`): `db.thermo` populated
(`group` category present, `primaryThermoLibrary` in libraries), `db.thermo.surface == {}` — the
metal/adsorption surface database is correctly not loaded.

**Control B — `surface=True` unchanged from `546d1c727`.** Captured comparable state (sorted key
sets of `db.thermo.groups`, `db.thermo.libraries`, `db.thermo.surface`, `MetalDatabase` type, and
the sorted list of adsorption-group entry labels) in two separate processes:

- `control_B_baseline_stdout.log` — captured on the reverted, pre-fix tree.
- `control_B_fixed_stdout.log` — captured on the fixed tree, same script.

```
$ diff docs/i217-thermo-surface-gate/control_B_baseline_stdout.log docs/i217-thermo-surface-gate/control_B_fixed_stdout.log
$ echo $?
0
```

No output, exit 0 — byte-identical. `surface=True` behaviour is unchanged by this fix.

## Test counts

`python -m pytest test/rmgpy/data/ -p no:cacheprovider --no-cov -q`:

| | baseline (`546d1c727`, no `rmgTest.py`) | after fix (with `rmgTest.py`) |
|---|---|---|
| passed | 387 | 390 |
| skipped | 8 | 7 |
| failed | 0 | 0 |
| collected | 395 | 397 |

No new failures (0 -> 0), which is the bar. Counts reconcile exactly: 390 = 387 (baseline) + 2 (new
`rmgTest.py` tests) + 1 (see below); skipped drops by exactly that 1.

The one skip->pass flip is `thermoTest.py::TestMolecularManipulationInvolvedInThermoEstimation::
test_deterministic_bicyclic_decomposition`. Its own docstring documents it as a known
non-deterministic test (RMG-Py issue #2562): it builds a `Molecule` directly and self-skips via
`pytest.skip(...)` inside an `except AssertionError` when a ring-decomposition ordering isn't
reproduced on that run. It does not touch `RMGDatabase`, `database.thermo`, or the `surface`
argument — its outcome depends on set/dict iteration order, not on this fix. Pre-existing
flakiness, unrelated to the change, not a new failure and not caused by it.

## Git-status / gitignore anomaly (flagged per the brief's own "report anything unexpected")

`git status` at commit time shows several stray untracked, non-gitignored files/directories at the
repo root: `.bash_profile`, `.bashrc`, `.claude/`, `.gitconfig`, `.gitmodules`, `.idea`, `.mcp.json`,
`.profile`, `.ripgreprc`, `.vscode`, `.zprofile`, `.zshrc`. None of these were created by this task;
they appear to be environment/worktree setup artifacts. `git check-ignore -v` returns nothing/exit 1
for all of them, confirmed against the actual `.gitignore` patterns (`.idea/*`, `*/.idea/*`,
`.vscode/*`, `*/.vscode/*` exist but do not match these literal root-level paths as currently
written). They are **not staged or committed** by this change.

Separately: `.gitignore` has no pattern covering `docs/i217-thermo-surface-gate/*` — its only log
patterns are narrowly scoped (`scripts/treegen.log`, `scripts/treegen_backup.log`,
`test/arkane/data/two_parameter_arrhenius_fit/arkane.log`, `test_log.txt`, `examples/**/arkane.log`).
So `stdout.log`/`stderr.log`/`report.md` under `docs/i217-thermo-surface-gate/` are **not**
auto-excluded by any ignore rule — the brief's non-goal ("do not commit ... if the repo's ignore
rules would otherwise pick them up") does not literally forbid committing them, since nothing picks
them up. The commit for this change is nonetheless kept scoped to the code fix and the regression
test only; `docs/i217-thermo-surface-gate/` (report + all scratch logs) is left untracked/uncommitted
on disk as evidence, not swept in by any `git add -A`/`.` (neither of which was used).

## What this may and may not be claimed to show

**May claim:** the specific silent-`None` defect in `RMGDatabase.load(surface=False)` is fixed by
an unconditional `load_thermo(...)` call; `surface=True` (the default, used by every non-test
caller found) is provably unchanged; the four existing test call sites that pass `surface=False`
do not depend on the previous skip; the only measurable cost is one test file's wall-clock, +30%
on a sub-second unit test.

**May not claim:** that every possible caller in downstream/user code (outside this tree) is safe
— only this repository was searched. That the process-global `RMGDatabase` singleton itself is
safe for concurrent or repeated loading in general (out of scope; flagged, not fixed, per the
non-goals). That the pytest-suite pass count of 390 is a "green suite" — 7 pre-existing skips
remain, one of which is a documented flaky test unrelated to this change. That timing was measured
with statistical rigor — one run each, not averaged; reported as a directional cost, not a
precise benchmark.
