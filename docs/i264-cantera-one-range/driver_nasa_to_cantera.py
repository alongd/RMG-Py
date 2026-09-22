"""
Driver for the i264 commit 2 verifier (RED/GREEN): rmgpy.thermo.nasa.NASA.to_cantera().

Builds a one-range NASA thermo object (single polynomial, 200-6000 K) and a
two-range one (non-regression), calls .to_cantera() on each, and reports the
outcome. Against the BASE checkout, the one-range case must raise an
AssertionError ("only accept 2 polynomials"); against this worktree it must
return a cantera.NasaPoly2 whose Cp/H/S match the original RMG polynomial
exactly (both ranges carry duplicated coefficients of the same polynomial).
"""
import sys
import traceback

from rmgpy.thermo.nasa import NASA, NASAPolynomial

AR_COEFFS = [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.37967]


def make_nasa(two_range: bool) -> NASA:
    if two_range:
        polys = [
            NASAPolynomial(coeffs=AR_COEFFS, Tmin=(200, "K"), Tmax=(1000, "K")),
            NASAPolynomial(coeffs=AR_COEFFS, Tmin=(1000, "K"), Tmax=(6000, "K")),
        ]
    else:
        polys = [
            NASAPolynomial(coeffs=AR_COEFFS, Tmin=(200, "K"), Tmax=(6000, "K")),
        ]
    return NASA(polynomials=polys, Tmin=(200, "K"), Tmax=(6000, "K"))


def main():
    print("=== one-range NASA.to_cantera() ===")
    nasa_one = make_nasa(two_range=False)
    try:
        ct_poly = nasa_one.to_cantera()
        print(f"[OK] to_cantera() succeeded: {ct_poly!r}")
        max_rel_err = 0.0
        for T in (250.0, 500.0, 1000.0, 2000.0, 4000.0, 5900.0):
            cp_rmg = nasa_one.get_heat_capacity(T)
            h_rmg = nasa_one.get_enthalpy(T)
            s_rmg = nasa_one.get_entropy(T)
            # cantera's cp/h/s are per kmol; convert to per mol for comparison.
            cp_ct, h_ct, s_ct = ct_poly.cp(T) / 1000.0, ct_poly.h(T) / 1000.0, ct_poly.s(T) / 1000.0
            for rmg_val, ct_val in ((cp_rmg, cp_ct), (h_rmg, h_ct), (s_rmg, s_ct)):
                rel_err = abs(rmg_val - ct_val) / max(abs(rmg_val), 1e-8)
                max_rel_err = max(max_rel_err, rel_err)
            print(f"  T={T:7.1f} Cp_RMG={cp_rmg:10.4f} Cp_CT={cp_ct:10.4f} "
                  f"H_RMG={h_rmg:12.3f} H_CT={h_ct:12.3f} "
                  f"S_RMG={s_rmg:10.4f} S_CT={s_ct:10.4f}")
        print(f"  max relative error: {max_rel_err:.3e}")
    except AssertionError as e:
        print(f"[EXPECTED-RED] AssertionError raised (base/pre-fix behavior): {e}")
        traceback.print_exc()
    except Exception as e:
        print(f"[UNEXPECTED] {type(e).__name__}: {e}")
        traceback.print_exc()

    print()
    print("=== two-range NASA.to_cantera() (non-regression) ===")
    nasa_two = make_nasa(two_range=True)
    try:
        ct_poly2 = nasa_two.to_cantera()
        print(f"[OK] to_cantera() succeeded: {ct_poly2!r}")
    except Exception as e:
        print(f"[UNEXPECTED] {type(e).__name__}: {e}")
        traceback.print_exc()


if __name__ == "__main__":
    main()
