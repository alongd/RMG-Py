# Contract — i229-r198-hardening

Opened: 2026-10-07 · Worktree: <worktree> · Base: i229-reland@54827dbe5

## Intent

Make the depository duplicate comparator conservative for units and runtime types:
values may collapse only when their exact supported classes, complete descriptor
inventories, dimensions, SI values, uncertainties, and enclosing kinetics fields all
match. Forward Arkane's `kineticsDepositories='all'` as inclusion of every depository,
and correct the family loader documentation to describe inclusion rather than precedence.

## Premise

`quantities.Quantity(...).simplified.dimensionality` retains dimensions while discarding
unit scale, so it distinguishes equal numeric SI values with incompatible dimensions but
continues to accept compatible representations such as cm³/(mol·s) and m³/(mol·s).
Exact runtime-type and descriptor-inventory checks make unmodelled subclass state fail
closed instead of being silently ignored. This was checked directly against the installed
RMG environment for scalar quantities, array quantities, and `RateUncertainty`.

## Verifier

Using `<python>`, the focused depository and Arkane
tests, `test/rmgpy/data`, `test/rmgpy/rmg/inputTest.py`, all Arkane input tests, and every
`*CensusTest.py` must complete with only the explicitly established baseline failures.
All commands run with this worktree first on `PYTHONPATH`; long runs are polled through
their persistent sessions and logged under `<run logs>`.

## Non-goals

- No new supported kinetics classes or persistence format.
- No scientific precedence among depositories; order remains deterministic inclusion only.
- No changes outside the named comparator, Arkane input mapping, family docstring, tests,
  census artifacts made stale by those edits, report, and this contract.
- Preserve the pre-existing untracked Chemkin YAML files.

## Gates

None. The user explicitly authorized amending the existing single commit and force-pushing
the named branch after verification.

## Evidence

- RED: all four incompatible-dimension shapes silently collapsed; both unknown
  subclasses silently collapsed; Arkane forwarded `'all'` as `None`.
- GREEN: focused comparator suite 52 passed; focused comparator plus Arkane 53 passed
  with 2 established skips; Arkane input files 6 passed with 4 established skips.
- Regression: RMG input 123 passed, 1 skipped; every census test 1049 passed; data
  suite 2747 passed, 8 skipped, with only the two established failures in
  `untypeableStructureTest` and `saturatedStructureTest`.
- Logs: `<run logs>`.

## Round r199 — NumPy storage hardening

### Intent

Make the reviewed `ArrayQuantity` identity branch reject unreviewed NumPy storage for
both values and uncertainties: ndarray subclasses, object arrays, and dtypes carrying
metadata. Ordinary exact float64 ndarrays must continue to prove independent records
identical.

### Premise

`np.array_equal` participates in NumPy's `__array_function__` dispatch, dtype equality
does not establish equal dtype metadata, and object arrays delegate element equality to
opaque Python objects. The supplied reproducer path had already been removed, but the
review records all three behaviors; RED tests reconstruct each through the public
depository conflict seam before production code changes.

### Verifier

With `<python>` and this worktree first on
`PYTHONPATH`, the focused conflict tests, `test/rmgpy/data`,
`test/rmgpy/rmg/inputTest.py`, and every `*CensusTest.py` must complete with only the
two established data-suite failures. Persistent long-running sessions are polled to
completion and logged under `<run logs>`.

### Non-goals

- Data migration / schema change: N/A; no persisted shape changes.
- External API or file-format compatibility: N/A; only duplicate proof becomes more
  conservative for unsafe in-memory array storage.
- Compute spend: N/A; local tests only.
- Shared or dirty checkouts: preserve the two pre-existing untracked Chemkin YAML files.
- Other people's files: N/A.
- Shared branches: amend and force-push only `origin/i229-reland`, explicitly authorized.
- Deletion: none.
- No changes to scalar comparison, dimension handling, supported kinetics classes,
  depository selection, Arkane, documentation, or census verdicts.

### Gates

The user explicitly authorized amending the pushed single commit and force-pushing with
lease after verification.

### Evidence

- The supplied `/dev/shm/r199-XFoeOQw5/reproduce_numpy.py` path no longer existed at
  execution time; the three recorded cases were reconstructed at the same public
  depository seam.
- RED: ndarray subclasses, unequal dtype metadata, and opaque object arrays silently
  collapsed for both `value_si` and `uncertainty_si` (6 failures); ordinary independent
  float64 arrays passed the identity control.
- GREEN: all 8 new storage cases passed; the full focused conflict suite passed 60.
- Regression: RMG input passed 123 with 1 skip; every census test passed 1049; data
  passed 2756 with 7 skips and only the two established failures in
  `untypeableStructureTest` and `saturatedStructureTest`.
- Logs: `<run logs>`.

## Round r200 — exact array dtype allowlist

### Intent

Replace the `ArrayQuantity` storage denylist with a closed allowlist for both stored
values and uncertainties. Exact base `np.ndarray` instances with exact native
`np.dtype('float64')` are the only arrays that can prove identity; structured, object,
metadata-bearing, non-native, narrow-float, complex, subclassed, or otherwise unknown
storage must refuse. Primitive leaves likewise use an explicit exact-type allowlist for
`None`, `bool`, `int`, `float`, and `str`.

### Premise

The only open choice was whether shipped depositories require integer arrays. An
inventory of all 141 shipped families, 183 loaded depositories, and 12,293 entries found
10 array fields, all exact native `float64`, and no integer storage. Therefore `int64`
is not included.

### Verifier

With `<python>` and this worktree first on
`PYTHONPATH`, reproduce each newly requested dtype through the public depository
loader-to-`get_kinetics` seam before and after the fix, then run the complete focused
conflict tests, `test/rmgpy/data`, `test/rmgpy/rmg/inputTest.py`, and every
`*CensusTest.py`. Long runs must be detached, polled to completion, and logged under
`<run logs>`; only the two established data-suite failures may
remain.

### Non-goals

- Data migration / schema change: N/A; no persisted data is rewritten.
- External API or file-format compatibility: unknown array storage becomes
  deliberately non-identical; accepted native float64 storage is unchanged.
- Compute spend: N/A; local database inventory and tests only.
- Shared or dirty checkouts: preserve all pre-existing untracked generated artifacts.
- Other people's files: N/A.
- Shared branches: amend and force-push only `origin/i229-reland`, explicitly authorized.
- Deletion: none.
- No depository ranking, new kinetics classes, database edits, or broader numeric
  coercion.

### Gates

The user explicitly authorized amending the existing pushed single commit and
force-pushing the named branch with lease protection after verification.

### Evidence

- Inventory: 141 shipped families, 183 depositories, and 12,293 entries contained 10
  array fields, all exact native float64; integer storage was absent.
- RED: the five requested unsupported dtypes silently collapsed across 16
  value/uncertainty and load-order cases; both native-float64 controls passed.
- GREEN: all 18 loader-path cases passed; the complete focused conflict suite passed 78.
- Regression: RMG input passed 123 with 1 skip; every census test passed 1049; data
  passed 2774 with 7 skips and only the two established failures in
  `untypeableStructureTest` and `saturatedStructureTest`.
- Census verdicts are unchanged; only the `_comparable` digest and mechanical reversal
  line numbers moved.
- Logs: `<run logs>`.
