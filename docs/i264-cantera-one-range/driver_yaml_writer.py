"""
Driver for the i264 verifier, step 1 (RED/GREEN), extended for round-98
hardening: exercises every named refusal that ``species_to_dict`` now raises,
in addition to the original one-range/two-range/zero-polynomial arms.

Run against the BASE plasma checkout (RMG-Py-plasma) to reproduce the
IndexError on the one-range arm (that arm is expected to FAIL there -- pass
``--base`` so this script does not also fail on that account), and against
this worktree, where every arm must pass and the driver exits 0.

Unlike the pre-round-98 version of this file, EVERY arm's outcome is tracked;
the driver calls ``sys.exit(1)`` if any arm did not match its expectation, so
a regression here cannot be missed by eyeballing printed text.
"""
import argparse
import sys
import traceback

from rmgpy.exceptions import CanteraThermoWriteError
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.thermo.nasa import NASA, NASAPolynomial
from rmgpy.yaml_cantera2 import species_to_dict

# A real single NASA7 coefficient set for monatomic Ar (Cp = 5/2 R constant,
# from NIST-style Ar thermo; a0=2.5, a1..a4=0, a5/a6 chosen for continuity).
AR_COEFFS = [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.37967]


def _species(label, polys):
    mol = Molecule(smiles="[Ar]")
    nasa = NASA(polynomials=polys, Tmin=(200, "K"), Tmax=(6000, "K"))
    return Species(label=label, molecule=[mol], thermo=nasa)


def make_species(two_range: bool) -> Species:
    if two_range:
        polys = [
            NASAPolynomial(coeffs=AR_COEFFS, Tmin=(200, "K"), Tmax=(1000, "K")),
            NASAPolynomial(coeffs=AR_COEFFS, Tmin=(1000, "K"), Tmax=(6000, "K")),
        ]
        return _species("Ar-two-range", polys)
    polys = [NASAPolynomial(coeffs=AR_COEFFS, Tmin=(200, "K"), Tmax=(6000, "K"))]
    return _species("Ar-one-range", polys)


def run_arm(name, fn, expect):
    """expect: "ok" -- fn() must succeed; or an exception class -- fn() must
    raise exactly (an instance of) that class. Returns True on match."""
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
            print(f"[OK] succeeded: thermo={result['thermo']}")
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
                              "one-range arm is EXPECTED to raise IndexError, "
                              "so that arm's outcome is excluded from the "
                              "exit-code decision.")
    args = parser.parse_args()

    results = {}

    sp_one = make_species(two_range=False)
    if args.base:
        # On base this is a known, expected IndexError -- report it but do
        # not let it flip the exit code.
        print("=== one-range species (base: IndexError expected) ===")
        try:
            species_to_dict(sp_one, [sp_one])
            print("[FAIL-FOR-BASE] one-range species did NOT raise IndexError on base")
        except IndexError:
            print("[EXPECTED-RED] IndexError raised (base/pre-fix behavior), as expected on base")
        except Exception as e:
            print(f"[UNEXPECTED] {type(e).__name__}: {e}")
            traceback.print_exc()
        print()
    else:
        results["one-range species"] = run_arm(
            "one-range species", lambda: species_to_dict(sp_one, [sp_one]), "ok")

    sp_two = make_species(two_range=True)
    results["two-range species (non-regression)"] = run_arm(
        "two-range species (non-regression)", lambda: species_to_dict(sp_two, [sp_two]), "ok")

    sp_empty = _species("Ar-empty", [])
    results["zero-polynomial species (named refusal)"] = run_arm(
        "zero-polynomial species (named refusal)",
        lambda: species_to_dict(sp_empty, [sp_empty]), CanteraThermoWriteError)

    sp_three = _species("Ar-three-range", [
        NASAPolynomial(coeffs=AR_COEFFS, Tmin=(200, "K"), Tmax=(1000, "K")),
        NASAPolynomial(coeffs=AR_COEFFS, Tmin=(1000, "K"), Tmax=(3000, "K")),
        NASAPolynomial(coeffs=AR_COEFFS, Tmin=(3000, "K"), Tmax=(6000, "K")),
    ])
    results[">2 polynomials (named refusal)"] = run_arm(
        ">2 polynomials (named refusal)",
        lambda: species_to_dict(sp_three, [sp_three]), CanteraThermoWriteError)

    sp_nasa9 = _species("Ar-nasa9", [
        NASAPolynomial(coeffs=list(AR_COEFFS) + [0.0, 0.0], Tmin=(200, "K"), Tmax=(6000, "K")),
    ])
    results["non-7-coefficient / NASA9 (named refusal)"] = run_arm(
        "non-7-coefficient / NASA9 (named refusal)",
        lambda: species_to_dict(sp_nasa9, [sp_nasa9]), CanteraThermoWriteError)

    sp_gapped = _species("Ar-gapped", [
        NASAPolynomial(coeffs=AR_COEFFS, Tmin=(200, "K"), Tmax=(1000, "K")),
        NASAPolynomial(coeffs=AR_COEFFS, Tmin=(1500, "K"), Tmax=(6000, "K")),
    ])
    results["gapped ranges (named refusal)"] = run_arm(
        "gapped ranges (named refusal)",
        lambda: species_to_dict(sp_gapped, [sp_gapped]), CanteraThermoWriteError)

    sp_inverted = _species("Ar-inverted", [
        NASAPolynomial(coeffs=AR_COEFFS, Tmin=(1000, "K"), Tmax=(200, "K")),
    ])
    results["inverted range (named refusal)"] = run_arm(
        "inverted range (named refusal)",
        lambda: species_to_dict(sp_inverted, [sp_inverted]), CanteraThermoWriteError)

    sp_nan = _species("Ar-nan", [
        NASAPolynomial(coeffs=[float("nan")] * 7, Tmin=(200, "K"), Tmax=(6000, "K")),
    ])
    results["NaN coefficient (named refusal)"] = run_arm(
        "NaN coefficient (named refusal)",
        lambda: species_to_dict(sp_nan, [sp_nan]), CanteraThermoWriteError)

    n_fail = sum(1 for ok in results.values() if not ok)
    print(f"=== summary: {len(results) - n_fail}/{len(results)} arms passed ===")
    if n_fail:
        print(f"[FAIL] {n_fail} arm(s) failed")
        sys.exit(1)
    print("[OK] all arms passed")
    sys.exit(0)


if __name__ == "__main__":
    main()
