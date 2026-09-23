# i264 — The Cantera writer assumes every NASA polynomial has two temperature ranges

## Summary

`rmgpy/yaml_cantera2.py::species_to_dict` crashed with a bare `IndexError` when a species carried a
**one-range** NASA polynomial (a single exact fit, as a monatomic species like metastable argon
gets). The Cantera YAML writer now emits a genuine one-range thermo block; the object-API sibling
`rmgpy/thermo/nasa.pyx::NASA.to_cantera`, which had the same defect, duplicates the single
coefficient set across two ranges because Cantera's Python object API has no single-range class.

Two source commits plus an evidence commit on `i264-cantera-one-range` (base `plasma@98d465d3b`):

- `1476e50ba` — Handle one-range NASA thermo in the Cantera YAML writer
- `d5d7d0940` — Handle a one-range NASA polynomial in `NASA.to_cantera()`
- `adda44907` — Add i264 verifier evidence (probes / RED-GREEN drivers / round-trip / tests)

## The load-bearing premise, probed against the installed Cantera

**Claim:** Cantera's YAML loader accepts a genuine one-range NASA7 block (`temperature-ranges` of
length 2, one coefficient set), even though its Python object API has no single-range NASA class.
A prior attempt on the Chemkin side left both Cantera sites alone on the unverified reasoning that
"Cantera natively carries a single range, so splitting would impose a limitation the format lacks."

**Command:** `PYTHONPATH=<worktree> /home/alon/anaconda3/envs/rmg_env/bin/python docs/i264-cantera-one-range/probe_cantera.py`

**Result — CONFIRMED (not inverted):**

```
=== cantera 3.1.0 ===
has NasaPoly1 attr: False
has NasaPoly2 attr: True
[OK]   one-range (2 breakpoints, 1 data): loaded. thermo type = NasaPoly2
[OK]   two-range duplicated (3 breakpoints, 2 data): loaded. thermo type = NasaPoly2
```

`NasaPoly1` does not exist in Cantera 3.1.0; `NasaPoly2` is the only NASA7 species-thermo object the
Python API exposes. But the YAML *schema* accepts a one-range block, loads it (as `NasaPoly2`
internally), and returns Cp/H/S byte-identical to a duplicated two-range block.

Consequently the two sites need **different** fixes:
- **YAML writer** → emit a genuine one-range block (option 1). Cleanest; preserves the fact that
  the fit is exact rather than pretending to a meaningless breakpoint.
- **Object API (`to_cantera`)** → option 1 impossible (no `NasaPoly1`), so duplicate the single
  coefficient set across two identical-content ranges, mirroring `rmgpy/chemkin.pyx`.

## Commit 1 — `rmgpy/yaml_cantera2.py::species_to_dict`

- New `CanteraThermoWriteError(OutputError)` in `rmgpy/exceptions.py`, fatal like
  `MechanismWriterError`: a silently skipped species would leave a hole in the exported file while
  the export reported success.
- Construction generalised over `sorted_polys` of any length ≥ 1. Empty list → `CanteraThermoWriteError`
  by name instead of a fall-through `IndexError`.
- `temperature-ranges` = `[Tmin of first poly] + [Tmax of each poly]`; `data` = each polynomial's 7
  coefficients in order. Reduces to the previous 2-range output byte-for-byte when `len == 2`.

**RED (base `/home/alon/Code/RMG-Py-plasma`):** one-range species → `IndexError`.
**GREEN (this worktree):** one-range species writes
`{'model': 'NASA7', 'temperature-ranges': [200.0, 6000.0], 'data': [[2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.37967]]}`;
two-range species byte-identical to old output; zero-polynomial species →
`CanteraThermoWriteError: Cannot write Cantera thermo for species 'Ar-empty': its NASA thermo carries no polynomials.`

## Round-trip proof (verifier step 2)

`driver_roundtrip.py` writes a full Cantera YAML phase via `save_cantera_model`, loads it back with
`ct.Solution`, and compares Cp/H/S from the loaded object against the **original RMG thermo object**
(never a synthesized intermediate), for a one-range and a two-range species at
T = 250/500/1000/2000/4000/5900 K.

```
--- one-range (Ar1R) ---
   T (K)       Cp_RMG        Cp_CT          H_RMG           H_CT        S_RMG         S_CT
   250.0      20.7862      20.7862      -1000.855      -1000.853     151.1847     151.1846
  1000.0      20.7862      20.7862      14588.780      14588.764     180.0005     180.0003
  5900.0      20.7862      20.7862     116441.062     116440.931     216.8950     216.8947
max relative error across Cp/H/S at all spot temps: 1.128e-06
--- two-range (Ar2R, non-regression) --- (identical numbers, same 1.128e-06)
```

Max relative error **1.128e-06**, identical for one-range and two-range (i.e. not fix-specific). Root
cause: `cantera.gas_constant/1000 = 8.31446261815324` J/mol/K vs `rmgpy.constants.R = 8.314472`
J/mol/K — a pre-existing cross-library constant-vintage mismatch, ~1.13e-6 relative, not a defect in
the writer. The emitted one-range YAML carries `temperature-ranges: [200.0, 6000.0]` (length 2), so
the round-trip proves a real one-range block, not a duplicated one. Independently re-run by the
overseer session: exit 0, same numbers.

## Commit 2 — `rmgpy/thermo/nasa.pyx::NASA.to_cantera`

When `len(polys) == 1`, duplicate the single polynomial's 7 coefficients across both halves of a
`NasaPoly2`, splitting at 1000 K when that is strictly inside `(Tmin, Tmax)` else at the midpoint —
the same rule as `chemkin.pyx`'s existing one-range branch (read only, not modified). Both halves are
the identical polynomial, so Cp/H/S are continuous across the breakpoint and Tint is
thermodynamically irrelevant. The 7-coefficient assertion is kept.

`nasa.pyx` is cythonized: rebuilt via `python utilities.py check-pydas` then
`python setup.py build_ext --inplace`; the loaded `.so` confirmed by `__file__` to be the freshly
built one in this worktree. `rmgpy/solver/settings.pxi` not committed.

**RED (base):** `AssertionError: Cantera NasaPoly2 objects only accept 2 polynomials`.
**GREEN:** `to_cantera()` succeeds; Cp/H/S agree to 1.128e-06 (same R-constant cause); two-range
unaffected.

## Positional-index census

`grep -rn -E "polynomials\[|sorted_polys\[|polys\[0\]|polys\[1\]" --include=*.py --include=*.pyx rmgpy arkane scripts`
→ exactly **3 sites**, all accounted for:

1. `rmgpy/yaml_cantera2.py` — target, fixed (commit 1); the remaining `sorted_polys[0].Tmin` seeds
   the breakpoint list and is correct for any N ≥ 1 (empty list raised above it).
2. `rmgpy/thermo/nasa.pyx::to_cantera` — object API, fixed (commit 2); new code is length-branched.
3. `rmgpy/chemkin.pyx` (~1787/1803/1804) — pre-existing one-range/two-range branch, already
   length-guarded before indexing; examined only, unchanged.

No other call site indexes a NASA polynomial list positionally. The grep spans `.py` and `.pyx`
across `rmgpy`, `arkane`, and `scripts`.

## Tests

Per-directory pytest with the rmg_env interpreter (the whole suite cannot collect in one invocation
because `test/rmgpy/data/rmgTest.py` and `test/rmgpy/rmg/rmgTest.py` share a basename):

- `yaml_cantera1Test.py yaml_cantera2Test.py i167CanteraExportPathTest.py tools/canteramodelTest.py
  yaml_writer/ thermo/nasaTest.py` → **106 passed, 10 skipped** (skips are the opt-in
  `TestRecentlyGeneratedCanteraYaml2GasOnly` suite).
- `chemkinTest.py` → **44 passed** (includes the existing one-range Chemkin-writer tests — confirms
  `chemkin.pyx` unaffected).

## What could not be reached

- The whole-tree pytest run (duplicate-basename collection failure — pre-existing, environmental).
  Per-directory runs cover every test file touching the Cantera writer, object-API NASA writer, and
  Chemkin writer.
- The ~1.13e-6 RMG/Cantera gas-constant discrepancy is documented but not fixed — pre-existing,
  unrelated, out of scope.
- No end-to-end argon+metastable RMG run was executed here; the fix was proven at the writer/object
  boundary with constructed one-range thermo, which is where the crash lived.

## Gate compliance

No edits under `rmgpy/molecule/`, `rmgpy/kinetics/`, `rmgpy/data/`, RMG-database, or thermo
fitting/storage. No push/merge/rebase. Source diff limited to `rmgpy/exceptions.py`,
`rmgpy/thermo/nasa.pyx`, `rmgpy/yaml_cantera2.py`.

## Round-98 hardening

The initial fix (above) made `species_to_dict` and `NASA.to_cantera` handle a one-range NASA
object, but left several adjacent failure modes silent or crash-by-`assert`/`IndexError` rather
than refused-by-name, and had only driver-script coverage rather than real pytest coverage. This
round closes both gaps without touching the fix's actual behavior for valid input.

**Construction-time validation.** Both `species_to_dict` (`rmgpy/yaml_cantera2.py`) and
`NASA.to_cantera` (`rmgpy/thermo/nasa.pyx`) now raise the named `CanteraThermoWriteError`
(`rmgpy/exceptions.py`) — not a bare `assert` (which vanishes under `python -O`) and not a
positional `IndexError`/`TypeError` — for every one of: more than 2 polynomials, a polynomial that
does not have exactly 7 coefficients (e.g. NASA9 data), a non-finite coefficient, an
inverted/degenerate range (`Tmin >= Tmax`), and (in `species_to_dict` only, since only there are
multiple polynomials stitched together) a non-contiguous range gap between consecutive
polynomials. `nasa.pyx` was rebuilt (`make build`) after editing; the `.so` was confirmed to
resolve from this worktree. `rmgpy/solver/settings.pxi` (auto-written by `check-pydas`) was not
committed.

**Real pytest coverage.** `test/rmgpy/thermo/nasaTest.py::TestNASA` gained
`test_to_cantera_one_range` (headline — this exact input crashed with a bare `AssertionError`
before the original fix) plus one refusal test per failure mode:
`test_to_cantera_refuses_more_than_two_polynomials`,
`test_to_cantera_refuses_non_seven_coefficient_polynomial`,
`test_to_cantera_refuses_non_finite_coefficients`. `test/rmgpy/yaml_cantera2Test.py::TestCanteraWriter2`
gained the mirror set for the dict-writer path: `test_species_to_dict_one_range` (headline) plus
`test_species_to_dict_refuses_more_than_two_polynomials`,
`..._refuses_non_seven_coefficient_polynomial`, `..._refuses_gapped_ranges`,
`..._refuses_inverted_range`, `..._refuses_non_finite_values`. A round-trip test,
`test_one_range_species_roundtrips_through_cantera`, writes a full Cantera YAML phase for a
one-range argon species via `save_cantera_model`, loads it back with `cantera.Solution(...)`, and
compares Cp/H/S against the original RMG thermo object at six spot temperatures (250-5900 K) to a
5e-6 relative tolerance — the same ~1.13e-6 RMG/Cantera gas-constant gap documented above, not a
new defect. Every comparison value is checked with `math.isfinite` before use, so a NaN comparison
value cannot be silently swallowed by `max(0.0, nan) == 0.0` the way a naive
`max(max_rel_err, rel_err)` accumulator would.

Verified (see `docs/i264-cantera-one-range/round98/`):
- `full_nasaTest.log`: `pytest test/rmgpy/thermo/nasaTest.py --no-cov` → **24 passed**.
- `full_yaml_cantera2Test.log`: `pytest test/rmgpy/yaml_cantera2Test.py --no-cov` → **36 passed, 5
  skipped** (skips are the pre-existing opt-in `TestRecentlyGeneratedCanteraYaml2GasOnly` suite,
  unaffected by this change).
- `nasa_headline_red_on_base.log` / `yaml_cantera2_headline_red_on_base.log`: run against the BASE
  checkout (`/home/alon/Code/RMG-Py-plasma`, commit `98d465d3b`, pre-fix). `NASA.to_cantera()` on a
  one-range object raises a bare `AssertionError: Cantera NasaPoly2 objects only accept 2
  polynomials`; `species_to_dict` on the same input raises a bare `IndexError: list index out of
  range`; and `from rmgpy.exceptions import CanteraThermoWriteError` itself fails with
  `ImportError` — i.e. the new refusal tests cannot even be collected against base, and the
  headline tests fail for the documented pre-fix reason. Confirmed empirically rather than assumed
  (an earlier characterization of the `to_cantera` path had guessed `IndexError`; direct testing
  showed the two paths fail differently — `to_cantera` via `assert`, `species_to_dict` via raw
  indexing — both are real, and both are named `CanteraThermoWriteError` now).

**Evidence drivers exit non-zero on failure.** All three drivers under
`docs/i264-cantera-one-range/` (`driver_yaml_writer.py`, `driver_nasa_to_cantera.py`,
`driver_roundtrip.py`) track every arm's pass/fail and call `sys.exit(1)` if any arm's outcome
does not match its expectation — printed `[OK]`/`[FAIL]` text is not the only signal. Each
driver's exit-code behavior was proven live: a one-line tamper (weakening an assertion's
tolerance, or replacing an expected `CanteraThermoWriteError` with `"ok"`) was introduced, the
driver was re-run and observed to exit 1 with the tampered arm named as `[FAIL]`
(`*_TAMPERED_stdout.log`, `*_TAMPERED_exitcode.txt`), the tamper was reverted from a saved copy,
`git diff --stat`/`grep` confirmed a clean revert back to the real fix, and the driver was re-run
to confirm exit 0 (`*_GREEN_stdout.log`, `*_GREEN_stderr.log`, `*_GREEN_exitcode.txt`). This
demonstrates the exit code is a real signal tied to the arm outcomes, not printed theater. The
tamper itself was never left in the committed driver files.

**`python -O` demonstration (Task D).** `optimize_flag_demo.log` runs the same script twice, once
under `python` and once under `python -O`, feeding `NASA.to_cantera()` a 3-range object and a
9-coefficient (NASA9-shaped) object. Both refusals are raised identically as
`CanteraThermoWriteError` under both invocations — proving the refusal is a real `raise`, not an
`assert` that `-O`'s `__debug__ == False` would strip.

**What did not change.** No behavior for valid one- or two-range NASA input changed. No template,
atom type, or database content was touched. `save_cantera_model` was not given reload-on-export
behavior. The ~1.13e-6 RMG/Cantera gas-constant mismatch was not chased — it is bounded by the 5e-6
round-trip tolerance and documented, not fixed, consistent with the original hard constraint.
