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

## Round 13 rework (adversarial review, base `efe072877`)

Rework brief: `docs/contracts/i201-wrapper-r13.md`. Four HIGH defects fixed in
`_reverse_rate_at_reference_potential` (module-level helper, `rmgpy/reaction.py`) and its call site in
`check_collision_limit_violation`; two pre-existing defects verified and written up report-only below
(not fixed, per brief).

### MUST-FIX 1 — restored the native `P == 0` raise

The reimplemented `PDepArrhenius` pressure branch had no `P == 0` guard: the bracket scan clamps
`ilow = ihigh = 0` when `P == 0` (`plist[i] >= 0` is true at `i=0`), silently returning the
lowest-pressure child instead of raising, unlike native `PDepArrhenius.get_rate_coefficient`
(`rmgpy/kinetics/arrhenius.pyx:2656-2657`). Added, as the first statement inside the pressure branch:

```python
if P == 0:
    raise ValueError(
        'No pressure specified to pressure-dependent PDepArrhenius.get_rate_coefficient().')
```

Not caught at the call site — it propagates, matching pre-`efe072877` behavior.

### MUST-FIX 2 — oracle interpolation tests

Every prior test used `P` equal to the grid midpoint, so `plow == phigh` and the log10-interpolation
line never ran. Added six oracle tests in `TestReverseRateAtReferencePotentialRound13` asserting the
helper agrees with native `PDepArrhenius.get_rate_coefficient` on ordinary (non-charge-transfer)
`Arrhenius` children: `P` strictly between two grid pressures, below the lowest, above the highest,
single-entry grid, duplicate-pressure grid, and (paired with MUST-FIX 1) `P == 0` raising `ValueError`.

### MUST-FIX 3 — malformed `PDepArrhenius(pressures=None, ...)` no longer summed

`quantity.Pressure(None)` returns `None`, so `PDepArrhenius(pressures=None, arrhenius=[...])` has
`.pressures is None`; the old `getattr(k, 'pressures', None) is not None` test was `False` for this
case and fell through to the summation branch, treating it like a `MultiArrhenius`. Fixed with a
module-level sentinel so *presence* of the attribute (not its truthiness) selects the branch:

```python
_MISSING = object()
...
pressures = getattr(kinetics, 'pressures', _MISSING)
if pressures is not _MISSING:
    ...
    plist = pressures.value_si   # None -> AttributeError
```

**Deviation from the brief, found and verified in this rework:** the brief's test plan asked to assert
that native's own `PDepArrhenius(pressures=None, ...).get_rate_coefficient(T, P)` also raises
`AttributeError`, as a same-failure-mode comparison. Verified interactively that it instead
**segfaults the process** — confirmed via a standalone `python -c` reproduction outside pytest (a
`pytest.raises` context cannot catch a segfault, so this had to be checked separately, and calling it
inside a test crashed the whole pytest run):

```
$ python -c "
from rmgpy.kinetics import Arrhenius, PDepArrhenius
a1 = Arrhenius(A=(1.0e10, 'cm^3/(mol*s)'), n=0.0, Ea=(20.0, 'kJ/mol'))
malformed = PDepArrhenius(pressures=None, arrhenius=[a1])
malformed.get_rate_coefficient(1000.0, 1.0e5)
"
Segmentation fault (core dumped)
```
whereas the helper itself (pure Python `getattr` access, not compiled Cython typed access) raises
cleanly:
```
>>> _reverse_rate_at_reference_potential(malformed, 1000.0, 1.0e5)
AttributeError: 'NoneType' object has no attribute 'value_si'
```
This is a separate, pre-existing, out-of-scope latent defect in native
`PDepArrhenius.get_rate_coefficient` (Cython typed access to `self.pressures.value_si` on a
None-valued attribute with no null check) — **new ticket, not fixed here**: a malformed
`PDepArrhenius(pressures=None, ...)` crashes the whole process instead of raising a catchable
exception. The round-13 test (`test_pdep_arrhenius_with_none_pressures_not_summed`) therefore only
exercises the helper, and documents the deviation in its own docstring.

### MUST-FIX 4 — narrowed the V0 leaf contract; BM/BEP refuse instead of fabricating

`getattr(k, 'V0', None) is not None` also matched the averaged charge-transfer forms
(`ArrheniusChargeTransferBM`, `arrhenius.pyx:3375`; `SurfaceChargeTransferBEP`, `surface.pyx:1094`)
whose plain `get_rate_coefficient`'s second argument is a reaction energy (dHrxn/dGrxn), not a
potential — evaluating them "at V0" fabricated a number. Distinguished by
`hasattr(kinetics, 'get_rate_coefficient_from_potential')` (present only on the averaged forms, absent
on the concrete evaluable `ArrheniusChargeTransfer`/`SurfaceChargeTransfer`); refuses via
`UnsupportedReverseRateError` instead. The refusal is folded into the same `try/except
UnsupportedReverseRateError` that already wraps `generate_reverse_rate_coefficient()` at the call
site, so the three parallel reverse lists (`collision_limit_r`, `kr_list`, `conditions_r`) stay in
lockstep: `kr` is computed inside the `try` and all three appends happen only after success.

### RED-before / GREEN-after

RED demonstrated against a rebuild of the *unmodified* `efe072877` source (swapped in via `git show
efe072877:rmgpy/reaction.py`, no `git stash` — restored afterward):

```
$ python -m pytest test/rmgpy/reactionTest.py -p no:cacheprovider -o addopts="" -q -k "TestReverseRateAtReferencePotentialRound13"
...
FAILED ...::test_pdep_arrhenius_p_zero_raises_matching_native            (MUST-FIX 1)
FAILED ...::test_pdep_arrhenius_with_none_pressures_not_summed           (MUST-FIX 3)
FAILED ...::test_arrhenius_charge_transfer_bm_leaf_refused               (MUST-FIX 4)
FAILED ...::test_surface_charge_transfer_bep_leaf_refused                (MUST-FIX 4)
FAILED ...::test_multi_arrhenius_wrapping_bm_child_refused               (MUST-FIX 4)
FAILED ...::test_check_collision_limit_violation_skips_reverse_for_wrapped_bm (MUST-FIX 4 integration)
6 failed, 8 passed, 97 deselected in 1.04s
```
(The 8 passes are the MUST-FIX 2 oracle-interpolation tests plus the MultiArrhenius/MultiPDepArrhenius
sentinel regression guards and the concrete-ACT/SCT-not-refused guard — that logic was already correct
in `efe072877`, just previously untested; expected to pass pre-rework too.)

GREEN after the rework (rebuilt):
```
$ python -m pytest test/rmgpy/reactionTest.py -p no:cacheprovider -o addopts="" -q -k "TestReverseRateAtReferencePotentialRound13"
14 passed, 97 deselected in 0.82s
```

Full-file suite, baseline (`efe072877`, unmodified, before this rework's test additions) reported
`95 passed, 2 skipped` (see contract `docs/contracts/closed/i201-wrapper-v0.md` Evidence section).
After this rework (14 new tests added, all four fixes applied):
```
$ python -m pytest test/rmgpy/reactionTest.py -p no:cacheprovider -o addopts="" -q
109 passed, 2 skipped in 0.84s
```
95 + 14 = 109 — no unrelated failures, all five pre-existing `TestReverseCollisionLimitPotentialWrappers`
tests and the i197 `TestReverseCollisionLimitPotential` leaf tests still pass unchanged.

### Report-only item 5 — dropped `electrons` in wrapper-recursion `Reaction(...)` construction

Verified against source. `Reaction.__init__`'s `electrons` parameter defaults to `electrons=0`
(`rmgpy/reaction.py:231`). The three wrapper-recursion branches in `generate_reverse_rate_coefficient`
each build a fresh helper `Reaction` without passing `electrons`:
```
rmgpy/reaction.py:1643   rxn = Reaction(reactants=self.reactants, products=self.products)   # PDepArrhenius branch
rmgpy/reaction.py:1652   rxn = Reaction(reactants=self.reactants, products=self.products)   # MultiArrhenius branch
rmgpy/reaction.py:1661   rxn = Reaction(reactants=self.reactants, products=self.products)   # MultiPDepArrhenius branch
```
so each `rxn.electrons` is `0`, regardless of `self.electrons`. Meanwhile the two per-child reverse-rate
builders that a recursive call into one of these branches can reach set electrons explicitly to the
negated original:
```
rmgpy/reaction.py:1536   kr = SurfaceChargeTransfer(alpha=kf.alpha.value, electrons=-1*self.electrons, V0=(V0,'V'))
rmgpy/reaction.py:1560   kr = ArrheniusChargeTransfer(alpha=kf.alpha.value, electrons=-1*self.electrons, V0=(V0,'V'))
```
but `self` inside those methods is the *helper* `rxn`, whose `electrons` is `0` (not the real
reaction's `electrons`), because the wrapper branch never propagated it. So a wrapped reverse
charge-transfer child generated through this real path gets `electrons = -1*0 = 0` instead of the
correct `-1*self.electrons` of the real reaction — the electron-transfer metadata on that generated
kinetics object is silently wrong (`0` instead of ∓1 etc.), even though the collision-limit check
itself (this round's four fixes) can look fine, since it only reads `V0` and the rate, not `electrons`.
**New ticket, not fixed here.**

### Report-only item 6 — positional-arg mismatch in `generate_reverse_rate_coefficient(kf.Tmin, kf.Tmax)`

Verified against source. `generate_reverse_rate_coefficient`'s signature:
```
rmgpy/reaction.py:1565   def generate_reverse_rate_coefficient(self, network_kinetics=False, Tmin=None, Tmax=None, surface_site_density=0):
```
The `PDepArrhenius` wrapper-recursion branch calls it positionally:
```
rmgpy/reaction.py:1646   kr.arrhenius.append(rxn.generate_reverse_rate_coefficient(kf.Tmin, kf.Tmax))
```
so `kf.Tmin` binds to the first positional parameter, `network_kinetics` (not `Tmin`), and `kf.Tmax`
binds to `Tmin` (not `Tmax`); `Tmax` is left at its default `None`. `kf.Tmin` is a `Quantity`, which is
truthy, so `network_kinetics` is effectively always "on" for this call, and the intended
temperature-range restriction (`Tmax`) is silently dropped. **New ticket, not fixed here.**

### Report-only item 7 — native `log10(0)` / `klow == 0` hazard is inherited, not introduced

Confirmed the same hazard class exists natively (`rmgpy/kinetics/arrhenius.pyx:2664-2667`):
```python
klow = alow.get_rate_coefficient(T)
khigh = ahigh.get_rate_coefficient(T)
if klow == khigh == 0.0: return 0.0
k = klow * 10 ** (log10(P / Plow) / log10(Phigh / Plow) * log10(khigh / klow))
```
Native only special-cases `klow == khigh == 0.0`; a `klow == 0.0` with `khigh != 0.0` (or vice versa)
still divides by/logs zero. The round-13 helper's own interpolation branch reproduces this same
formula for exact agreement with native (that is the point of the MUST-FIX 2 oracle tests). Left as-is
per the brief: "fixing" it here would break native-agreement semantics; it is native's own hazard,
not one this rework introduced.

### Build

`python setup.py build_ext --inplace -j 8` rebuilt `rmgpy/reaction.cpython-39-x86_64-linux-gnu.so`
four times this session (initial rework GREEN build; RED build against unmodified `efe072877`; RED
full-file suite confirmation on the same build; final restore-and-rebuild GREEN confirmation), each
completing with the standard `copying build/lib.../rmgpy/reaction.cpython-39-x86_64-linux-gnu.so ->
rmgpy` tail and exit 0. Full RED/GREEN and full-file-suite output for this round is in
`docs/contracts/i201-wrapper-r13.md`'s Evidence section.

### Note on the r13 brief's second contract-evidence target

The r13 brief also asked to append verifier output to `docs/contracts/i201-wrapper-v0.md`; the
actual path is `docs/contracts/closed/i201-wrapper-v0.md` — that is the prior, already-closed i201
contract (base `4910d8c52`), whose own Evidence section already documents its own RED/GREEN and
full-suite runs for the wrapper-recursion fix it covered. This round's evidence belongs to, and was
appended to, the *round-13* contract (`docs/contracts/i201-wrapper-r13.md`) instead — re-appending
duplicate round-13 evidence to a closed, unrelated contract would misattribute it. Flagging the path
discrepancy here rather than mutating a closed contract file.
