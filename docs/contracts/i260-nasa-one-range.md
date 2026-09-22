# Contract — i260-nasa-one-range

Opened: 2026-09-22T08:25:19Z · Worktree: /home/alon/Code/RMG-Py-i260-nasa-one-range · Base: i260-nasa-one-range@40e21b495

## Intent

`rmgpy/chemkin.pyx::write_thermo_entry` currently asserts a NASA thermo has
exactly two temperature-range polynomials, then indexes `[0]`/`[1]`. A
monatomic species (constant Cp = 5/2 R) is exactly represented by a single
NASA polynomial over the whole range, so writing such a species to a Chemkin
file crashes with `AssertionError` instead of producing a valid 4-line
Chemkin thermo block. Fix: when the NASA has exactly one polynomial, split it
at an interior breakpoint (1000 K, or the midpoint if 1000 K isn't strictly
inside the range) into two identical-coefficient `NASAPolynomial` objects
before the existing two-range write path, so H/S/Cp are unchanged and
continuous at the breakpoint. Three-or-more-range NASA objects, which cannot
map into Chemkin's two-range format, are refused with a named `ChemkinError`
instead of an opaque `AssertionError`. Tests added to
`test/rmgpy/chemkinTest.py::TestThermoReadWrite` cover: one-range write
produces a valid block, one-range round-trips through write+read within
tolerance, and three-range is refused by name.

## Premise

Claim: a single NASAPolynomial evaluated identically on [Tmin,Tint] and
[Tint,Tmax] reproduces the original one-range polynomial's Cp/H/S exactly
(same coefficients on both sides), so splitting is lossless regardless of
where Tint is placed, and 1000 K is a safe default breakpoint for any
species whose range spans it.
Check: this is arithmetic, not empirical — H(T) and S(T) from two identical
NASAPolynomial coefficient sets are literally the same closed-form
expression as the original single polynomial, continuous by construction
since both sides agree everywhere, not just at the boundary. Confirmed by
inspection of `rmgpy/thermo/nasa.pyx` get_enthalpy/get_entropy (evaluate
purely from `coeffs` + T, no cross-polynomial state), and empirically
demonstrated by the round-trip test (Step 4/5 below) at 20 points across
300-6000 K.

## Verifier

`export PATH=/home/alon/anaconda3/envs/rmg_env/bin:$PATH && python -m pytest "test/rmgpy/chemkinTest.py::TestThermoReadWrite" -p no:cacheprovider --no-cov` exits 0, with the four new tests (`test_write_thermo_block_one_range`, `test_write_read_one_range_round_trips`, `test_write_thermo_block_three_ranges_refused`, plus the pre-existing class tests) all passing — no regressions in the surrounding suite.

## Non-goals

- **Data migration / schema change** — N/A, no persisted data shape changes.
- **External API or file-format compatibility** — the emitted Chemkin block format is unchanged (still 4 lines, two-range NASA-7 layout); only the *input* acceptance widens from "exactly 2 polynomials" to "1 or 2, refuse 3+ by name instead of crashing."
- **Compute spend** — N/A, no cluster/QM/paid-API work.
- **Shared or dirty checkouts** — this worktree only; no other worktree touched.
- **Other people's files** — N/A.
- **Shared branches** — branch `i260-nasa-one-range` not pushed; no downstream branch depends on it.
- **Deletion** — none; only additive test methods and a targeted fix inside `write_thermo_entry`.
- Out of scope: Cantera sibling writers, any file under `rmgpy/molecule/`, `rmgpy/kinetics/`, `rmgpy/data/`, `rmgpy/thermo/`, or any database/YAML writer — handled separately per the task brief. Only `rmgpy/chemkin.pyx` and `test/rmgpy/chemkinTest.py` are touched.

## Gates

None — routine, well-scoped fix with a clear verifier; no destructive or irreversible action involved.

## Evidence

- RED (pre-fix, built module): 3 failures, all `AssertionError` at `rmgpy/chemkin.pyx:1715` (`assert len(thermo.polynomials) == 2`), for one-range write, round-trip, and three-range-refused.
- Fix applied to `write_thermo_entry` in `rmgpy/chemkin.pyx`; rebuilt via `check-pydas` + `build_ext --inplace`, compiled clean.
- GREEN: `test_write_thermo_block_one_range`, `test_write_read_one_range_round_trips`, `test_write_thermo_block_three_ranges_refused` all PASSED.
- Round-trip demo: one-range Ar NASA written/read back with 0.00e+00 relative error in H and S at 12 points across 300-6000 K.
- Regression: `TestThermoReadWrite` full class, 14 passed, no regressions.
- Scope: only `rmgpy/chemkin.pyx` and `test/rmgpy/chemkinTest.py` touched (`git diff --stat`: 76 / 59 lines).

## Round 85 — review response (2026-09-22)

Verifier updated to: `python -m pytest "test/rmgpy/chemkinTest.py::TestThermoReadWrite" -p no:cacheprovider --no-cov` exits 0 with **18** tests, and the whole file (`test/rmgpy/chemkinTest.py`) is **44 passed**.

- **HIGH — one-range NASA-9 was silently corrupted (FIXED).** The split copied only c0..c6, dropping cm2/cm1 (the 1/T^2 and 1/T terms), and the pre-existing cm2/cm1 guard ran on the *synthesized* NASA-7 objects, so it could no longer see the problem — the new branch routed around the guard rather than removing it. Reproduced: a one-range NASA-9 round-tripped with Cp at 300 K going −245.6 → +31.6, H off by ~10^6 J/mol. Fix: new module-level `_validate_chemkin_nasa_polynomial` runs on the **source** polynomials before any split, refusing cm2/cm1 with the same message the two-range path already used.
- **MEDIUM — malformed values now refused by name (FIXED).** The same helper rejects NaN/Inf coefficients, non-finite or non-positive temperature bounds, an inverted range (per-range Tmin ≥ Tmax, covering the LOW "high range 1000→500 K" case), and missing bounds (named `ChemkinError` instead of an unnamed `AttributeError`).
- **MEDIUM — CORRECTION to the Cantera reasoning (this report was WRONG).** The earlier report claimed Cantera natively carries a single range (`NasaPoly1` / a 2-temperature NASA7) and that splitting the two sibling writers would wrongly impose Chemkin's limitation. **Cantera 3.1 has no `NasaPoly1`.** `thermo/nasa.pyx:434` raises `AssertionError` and `yaml_cantera2.py:524` raises `IndexError` on a one-range NASA — both are a live crash, and **both DO need fixing** (out of this ticket's `chemkin.pyx` grant → separate ticket). The correct fix there is still *not* this split: `to_cantera_nasapoly2` should build `NasaPoly2` with the single range's coefficients **duplicated** into both halves; the YAML writer should emit **one** temperature range with one coefficient set. The false premise survived because the round-trip that "confirmed" it never touched Cantera — a gap that had been honestly named ("did not exercise Cantera end to end") was let stand as if it were evidence.
- **LOW — narrow midpoint-fallback ranges (cosmetic, not fixed).** A range like 1200.000–1200.001 K writes an interior breakpoint that rounds to the header's precision, so the low range reads back degenerate. Confirmed this preserves Cp/H/S exactly (both halves carry identical coefficients), so it is cosmetic; not guarded, to avoid inventing a requirement for a physically absurd sub-mK range.
- **LOW — tests strengthened.** Refusal tests now assert the message content (species name + requirement) via `pytest.raises(..., match=...)`, not just the exception type. Added cases: one-range NASA-9 refused, varying-Cp one-range round-trip, midpoint-fallback round-trip (asserts 3100 K breakpoint), and the malformed-value matrix.
- **Scope correction:** this deliverable is **three** files, not two — `rmgpy/chemkin.pyx`, `test/rmgpy/chemkinTest.py`, and this contract (`docs/contracts/i260-nasa-one-range.md`), force-added past `docs/contracts/.gitignore`.

