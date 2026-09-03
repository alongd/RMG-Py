# i201-wrapper — report

Contract: `docs/contracts/i201-wrapper-v0.md`. Worktree `/home/alon/Code/RMG-Py-i201-wrapper`, branch
`i201-wrapper`, base `4910d8c52`.

## Completed

- Added module-level `_reverse_rate_at_reference_potential(kinetics, T, P)` in `rmgpy/reaction.py`,
  immediately above `class Reaction` (was line 86, now spans roughly lines 86-141 after the extended
  comment). Type-agnostic `getattr` chain, no isinstance/class tuple: `V0` present → charge-transfer
  leaf, evaluated at its own `v0.value_si`; `pressures` present → `PDepArrhenius`-shaped wrapper,
  bracketed and log-interpolated with each bracketing child recursed at its own potential; `arrhenius`
  present (and no `pressures`) → summing wrapper (`MultiArrhenius`/`MultiPDepArrhenius`), children
  summed via recursion; otherwise → unchanged `kinetics.get_rate_coefficient(T, P)`.
- Replaced `check_collision_limit_violation`'s old `reverse_v0 = getattr(...)` if/else (previously
  `reaction.py:2065-2069`) with a single call to the new helper, and extended the explanatory comment
  block above it to explain the wrapper-recursion rationale (kept the existing V0/Butler-Volmer
  explanation intact, per the brief).
- Added `TestReverseCollisionLimitPotentialWrappers` to `test/rmgpy/reactionTest.py`, reusing the i197
  infrastructure (`_ReverseSpyReaction`, `_RecordingReverseACT`, `_reverse_check_species`,
  `_REVERSE_V0 = 0.4`, `_REVERSE_PRESSURE = 1.0e5`): 5 tests, described under Verification.

## Deviation from the brief (disclosed, not a design change)

The brief's prescribed pseudocode calls `kinetics.get_adjacent_expressions(P)` on the `PDepArrhenius`
branch. `rmgpy/kinetics/arrhenius.pxd:177` declares that method `cdef`, not `cpdef` — it is a Cython
non-public method, unreachable via Python attribute access, and `rmgpy/reaction.pxd`'s cimport list
does not cimport `PDepArrhenius` from `arrhenius.pxd` either (confirmed by reading both files; only
`Arrhenius`, `SurfaceArrhenius`, `StickingCoefficient`, `SurfaceChargeTransfer` are cimported from
kinetics). Calling `kinetics.get_adjacent_expressions(P)` from `reaction.py` therefore raises
`AttributeError` — confirmed empirically in the prior session before this one.

The contract's own Non-goals section says "use the object's own `get_adjacent_expressions` for
pressure bracketing," so this is a premise failure in the contract/brief, not free license to redesign.
Per **"probe premises first"** and **"don't build to pass the check"**, I did not chase a `.pxd` change
(contract-forbidden) or invent a workaround that defeats the intent. Instead I reimplemented the same
`ilow`/`ihigh` bracket scan `get_adjacent_expressions` performs (`arrhenius.pyx:2628-2646`) directly
against `pressures.value_si` and `arrhenius`, both of which **are** Python-accessible (`arrhenius.pxd:174-175`,
`cdef public`). This reproduces the identical bracketing/interpolation logic — verified numerically
identical to the native wrapper for non-charge-transfer children (test 4, below) — without touching any
`.pxd` or `rmgpy/kinetics/` file. No `.pxd`/`kinetics` change was made; the fix lives entirely in
`rmgpy/reaction.py`.

## Verification

### 1. Build tail (`.so` rebuild confirmation)

Final build after the fix was reapplied (`build_final.log`), confirming `rmgpy/reaction.cpython-39-x86_64-linux-gnu.so`
was recompiled and copied:

```
gcc -pthread -B .../compiler_compat -shared ... build/temp.../rmgpy/reaction.o -o build/lib.../rmgpy/reaction.cpython-39-x86_64-linux-gnu.so
...
copying build/lib.linux-x86_64-cpython-39/rmgpy/reaction.cpython-39-x86_64-linux-gnu.so -> rmgpy
```

Three rebuilds were performed this session (RED, baseline-reconfirm, final/GREEN) — each preceded a
test run, per the brief's "rebuild before every test run" rule.

### 2. RED-before / GREEN-after (targeted wrapper tests)

Load at time of RED run: `/proc/loadavg` showed `3.38 11.39 14.54` (checked before timing).

**RED** — `rmgpy/reaction.py` reverted to HEAD (pre-fix, `git checkout -- rmgpy/reaction.py`), rebuilt,
then:
```
python -m pytest test/rmgpy/reactionTest.py -p no:cacheprovider -o addopts="" -q -k "TestReverseCollisionLimitPotentialWrappers"
```
Result: **4 failed, 1 passed, 92 deselected in 1.37s.**
- `test_multi_arrhenius_of_charge_transfer_evaluated_at_v0` (RED-before #1): `TypeError: Cannot convert
  _RecordingReverseACT to rmgpy.kinetics.arrhenius.Arrhenius` at `arrhenius.pyx:2790`, inside
  `MultiArrhenius.get_rate_coefficient`'s `for arrh in self.arrhenius:` loop — exactly the premise's
  predicted crash.
- `test_pdep_arrhenius_of_charge_transfer_evaluated_at_v0` (RED-before #2): `AssertionError: [0.0] ==
  [0.4]` — native PDepArrhenius handed the bracketing child only T, defaulting V to 0.0.
- `test_multi_pdep_arrhenius_of_pdep_of_charge_transfer_evaluated_at_v0` (RED-before #3): same
  `[0.0] == [0.4]` mismatch, two levels of wrapping.
- `test_helper_matches_native_for_non_charge_transfer_wrappers` (regression guard #4): failed with
  `ImportError: cannot import name '_reverse_rate_at_reference_potential' from 'rmgpy.reaction'`. This
  is expected, not a bug: the helper is a fix-introduced symbol, so this test cannot be green before the
  fix exists at all — it is a true regression guard *from the moment the fix lands onward*, not "green
  before and after" in the literal sense of "before this session wrote the fix into source." Flagged
  here per the dispatcher's instruction to disclose this nuance.
- `test_top_level_potential_tests_still_pass` (regression guard #5, sanity wrapper around the i197
  leaf-path tests): **passed** — confirms the untouched i197 top-level fix is unaffected.

**GREEN** — fix reapplied (`git apply` of the saved patch), rebuilt, then the same command:
```
python -m pytest test/rmgpy/reactionTest.py -p no:cacheprovider -o addopts="" -q -k "TestReverseCollisionLimitPotentialWrappers"
```
Result: **5 passed, 92 deselected in 0.68s.** All three RED-before tests now pass (`recorded == [0.4,
0.4]`, `[0.4]`, `[0.4]` respectively, with `0.0` and `1.0e5` asserted absent from each), and both
regression guards pass.

### 3. Full-file suite: baseline and after

**Baseline** (pre-fix `rmgpy/reaction.py`, i.e. HEAD, rebuilt), full file, no `-k` filter:
```
pytest test/rmgpy/reactionTest.py -p no:cacheprovider -o addopts="" -q
```
Result: **4 failed, 91 passed, 2 skipped in 0.97s.** The 4 failures are exactly the 3 RED-before tests
plus the import-based regression guard described above — no other test in the file is affected; the
pre-existing i197 tests and all other reactionTest.py tests pass unchanged.

**After** (fix reapplied, rebuilt), same command:
```
pytest test/rmgpy/reactionTest.py -p no:cacheprovider -o addopts="" -q
```
Result: **95 passed, 2 skipped in 0.78s.** (2 skips are pre-existing and unrelated to this change —
present in both runs, same count.)

## Contract questions

**Q2 (Premise — is "call the wrapper's own get_rate_coefficient" viable?).** No — confirmed FALSE for
`MultiArrhenius`: `rmgpy/kinetics/arrhenius.pyx:2787-2792` declares `cdef Arrhenius arrh` and iterates
`for arrh in self.arrhenius:`; `ArrheniusChargeTransfer` subclasses `KineticsModel`
(`arrhenius.pyx:2960`), not `Arrhenius` (`arrhenius.pyx:49`), so binding raises `TypeError` — reproduced
live in the RED run above. `PDepArrhenius`/`MultiPDepArrhenius` (`arrhenius.pyx:2647-2668`,
`2901-2919`) type children loosely enough to hold a charge-transfer child, but call
`alow.get_rate_coefficient(T)` with only `T`, so V silently defaults to 0.0 — also reproduced live.
Confirms the fix must recurse in `reaction.py` and bypass the native wrapper methods for
charge-transfer wrappers, which is what was implemented.

**Q3 (Reachability — can `generate_reverse_rate_coefficient` actually produce a wrapper-of-charge-
transfer today?).** Read in full (`rmgpy/reaction.py:1530` onward). Three branches recurse per child and
can structurally produce a wrapper whose children are charge-transfer:
- `PDepArrhenius` branch (`~reaction.py:1600-1607`): builds a fresh `PDepArrhenius`, and for each
  `kinetics` in `kf.arrhenius`, sets `rxn.kinetics = kinetics` and calls
  `rxn.generate_reverse_rate_coefficient(...)` recursively — if a child is `ArrheniusChargeTransfer`,
  the `elif isinstance(kf, ArrheniusChargeTransfer)` branch (~line 1560) returns
  `reverse_arrhenius_charge_transfer_rate`, itself an `ArrheniusChargeTransfer` (carries `V0`). So a
  `PDepArrhenius` wrapping charge-transfer children is structurally producible.
- `MultiArrhenius` branch (~1609-1615) and `MultiPDepArrhenius` branch (~1617-1623): same recursive
  shape, same conclusion.
No forward kinetics that is actually a `PDepArrhenius`/`MultiArrhenius`/`MultiPDepArrhenius` wrapping
`ArrheniusChargeTransfer`/`SurfaceChargeTransfer` children exists in the codebase or database today —
so the shape is *latent*: reachable by construction of `generate_reverse_rate_coefficient`, not by any
currently-generated forward kinetics. This matches the contract's Intent framing ("a guard against a
future widening of reverse generation, not a live-bug fix").

(The contract explicitly numbers only "contract question 2"; "q1" and "q3" are not separately numbered
in the contract text itself. Q3 above answers the brief's own labelled "contract q3" reachability ask.
No further distinct numbered question was found in the contract to answer as "q1" — the Intent and
Non-goals sections do not pose one.)

## Remaining work

None outstanding against this contract's scope. Not done (out of scope per "Do NOT" / Non-goals):
no `.pxd`/`rmgpy/kinetics/` change, no isinstance allowlist, no touch to the forward
`get_rate_coefficient` or the i197 leaf semantics.

## `.so` rebuild confirmation

Confirmed three times this session (RED build, baseline-reconfirm build, final/GREEN build) — the
final build tail is quoted in section 1 above; `rmgpy.cpython-39-x86_64-linux-gnu.so` for `reaction.py`
was recompiled and copied on each rebuild, and `from rmgpy.reaction import
_reverse_rate_at_reference_potential` succeeds against the final `.so` (used directly by the
regression-guard test, which imports it and passed in the GREEN/final runs).
