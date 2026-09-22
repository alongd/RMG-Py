"""
Driver for the i264 commit 2 verifier (RED/GREEN): rmgpy.thermo.nasa.NASA.to_cantera().

Builds a one-range NASA thermo object (single polynomial, 200-6000 K) and a
two-range one (non-regression), calls .to_cantera() on each, and reports the
outcome. Against the BASE checkout, the one-range case raises an
AssertionError ("only accept 2 polynomials"); against this worktree it
returns a cantera.NasaPoly2 whose Cp/H/S match the original RMG polynomial
exactly (both ranges carry duplicated coefficients of the same polynomial).

Round-98 hardening additionally exercises every named refusal that
``to_cantera`` now raises (as ``CanteraThermoWriteError``, not a bare
``assert`` that vanishes under ``python -O``): more than 2 polynomials,
non-7-coefficient (NASA9-shaped) data, and non-finite coefficients. Every
arm's outcome is tracked and the driver exits non-zero if any arm's
outcome does not match its expectation -- eyeballing printed [OK]/[FAIL]
text is no longer how this is graded.
"""
import argparse
import math
import sys
import traceback

from rmgpy.exceptions import CanteraThermoWriteError
from rmgpy.thermo.nasa import NASA, NASAPolynomial

AR_COEFFS = [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.37967]


def make_nasa(two_range: bool) -> NASA:
    if two_range:
        polys = [
            NASAPolynomial(coeffs=AR_COEFFS, Tmin=(200, "K"), Tmax=(1000, "K")),
            NASAPolynomial(coeffs=AR_COEFFS, Tmin=(1000, "K"), Tmax=(6000, "K")),
        ]
    else:
        polys = [NASAPolynomial(coeffs=AR_COEFFS, Tmin=(200, "K"), Tmax=(6000, "K"))]
    return NASA(polynomials=polys, Tmin=(200, "K"), Tmax=(6000, "K"))


def check_roundtrip(nasa_obj, ct_poly):
    """Returns (ok, max_rel_err). Rejects non-finite comparison values
    explicitly via math.isfinite instead of letting max(0.0, rel_err)
    silently swallow a NaN."""
    max_rel_err = 0.0
    for T in (250.0, 500.0, 1000.0, 2000.0, 4000.0, 5900.0):
        cp_rmg = nasa_obj.get_heat_capacity(T)
        h_rmg = nasa_obj.get_enthalpy(T)
        s_rmg = nasa_obj.get_entropy(T)
        # cantera's cp/h/s are per kmol; convert to per mol for comparison.
        cp_ct, h_ct, s_ct = ct_poly.cp(T) / 1000.0, ct_poly.h(T) / 1000.0, ct_poly.s(T) / 1000.0
        for name, rmg_val, ct_val in (("Cp", cp_rmg, cp_ct), ("H", h_rmg, h_ct), ("S", s_rmg, s_ct)):
            if not (math.isfinite(rmg_val) and math.isfinite(ct_val)):
                print(f"  [FAIL] non-finite comparison value at T={T}: {name}_RMG={rmg_val}, {name}_CT={ct_val}")
                return False, float("nan")
            rel_err = abs(rmg_val - ct_val) / max(abs(rmg_val), 1e-8)
            max_rel_err = max(max_rel_err, rel_err)
        print(f"  T={T:7.1f} Cp_RMG={cp_rmg:10.4f} Cp_CT={cp_ct:10.4f} "
              f"H_RMG={h_rmg:12.3f} H_CT={h_ct:12.3f} "
              f"S_RMG={s_rmg:10.4f} S_CT={s_ct:10.4f}")
    print(f"  max relative error: {max_rel_err:.3e}")
    return max_rel_err < 5e-6, max_rel_err


def run_arm(name, fn, expect):
    print(f"=== {name} ===")
    try:
        result = fn()
    except Exception as e:
        if expect == "ok":
            print(f"[FAIL] unexpected {type(e).__name__}: {e}")
            traceback.print_exc()
            ok = False
        elif isinstance(e, expect):
            print(f"[OK] raised {type(e).__name__} as expected: {e}")
            ok = True
        else:
            print(f"[FAIL] expected {expect.__name__}, got {type(e).__name__}: {e}")
            traceback.print_exc()
            ok = False
    else:
        if expect == "ok":
            print(f"[OK] succeeded: {result!r}")
            ok = True
        else:
            print(f"[FAIL] expected {expect.__name__} to be raised, but succeeded: {result!r}")
            ok = False
    print()
    return ok


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", action="store_true",
                         help="Running against the pre-fix BASE checkout: the "
                              "one-range arm is EXPECTED to raise AssertionError, "
                              "excluded from the exit-code decision there.")
    args = parser.parse_args()

    results = {}

    nasa_one = make_nasa(two_range=False)
    if args.base:
        print("=== one-range NASA.to_cantera() (base: AssertionError expected) ===")
        try:
            nasa_one.to_cantera()
            print("[FAIL-FOR-BASE] one-range case did NOT raise on base")
        except AssertionError as e:
            print(f"[EXPECTED-RED] AssertionError raised (base/pre-fix behavior): {e}")
        except Exception as e:
            print(f"[UNEXPECTED] {type(e).__name__}: {e}")
            traceback.print_exc()
        print()
    else:
        def _one():
            ct_poly = nasa_one.to_cantera()
            ok, _ = check_roundtrip(nasa_one, ct_poly)
            if not ok:
                raise AssertionError("round-trip comparison failed or exceeded tolerance")
            return ct_poly
        results["one-range NASA.to_cantera()"] = run_arm(
            "one-range NASA.to_cantera()", _one, "ok")

    nasa_two = make_nasa(two_range=True)

    def _two():
        ct_poly2 = nasa_two.to_cantera()
        ok, _ = check_roundtrip(nasa_two, ct_poly2)
        if not ok:
            raise AssertionError("round-trip comparison failed or exceeded tolerance")
        return ct_poly2
    results["two-range NASA.to_cantera() (non-regression)"] = run_arm(
        "two-range NASA.to_cantera() (non-regression)", _two, "ok")

    nasa_three = NASA(polynomials=[
        NASAPolynomial(coeffs=AR_COEFFS, Tmin=(200, "K"), Tmax=(1000, "K")),
        NASAPolynomial(coeffs=AR_COEFFS, Tmin=(1000, "K"), Tmax=(3000, "K")),
        NASAPolynomial(coeffs=AR_COEFFS, Tmin=(3000, "K"), Tmax=(6000, "K")),
    ], Tmin=(200, "K"), Tmax=(6000, "K"))
    results[">2 polynomials (named refusal)"] = run_arm(
        ">2 polynomials (named refusal)", nasa_three.to_cantera, CanteraThermoWriteError)

    nasa_nasa9 = NASA(polynomials=[
        NASAPolynomial(coeffs=list(AR_COEFFS) + [0.0, 0.0], Tmin=(200, "K"), Tmax=(6000, "K")),
    ], Tmin=(200, "K"), Tmax=(6000, "K"))
    results["non-7-coefficient / NASA9 (named refusal)"] = run_arm(
        "non-7-coefficient / NASA9 (named refusal)", nasa_nasa9.to_cantera, CanteraThermoWriteError)

    nasa_nan = NASA(polynomials=[
        NASAPolynomial(coeffs=[float("nan")] * 7, Tmin=(200, "K"), Tmax=(6000, "K")),
    ], Tmin=(200, "K"), Tmax=(6000, "K"))
    results["NaN coefficient (named refusal)"] = run_arm(
        "NaN coefficient (named refusal)", nasa_nan.to_cantera, CanteraThermoWriteError)

    n_fail = sum(1 for ok in results.values() if not ok)
    print(f"=== summary: {len(results) - n_fail}/{len(results)} arms passed ===")
    if n_fail:
        print(f"[FAIL] {n_fail} arm(s) failed")
        sys.exit(1)
    print("[OK] all arms passed")
    sys.exit(0)


if __name__ == "__main__":
    main()
