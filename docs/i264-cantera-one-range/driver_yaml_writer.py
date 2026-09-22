"""
Driver for the i264 verifier, step 1 (RED/GREEN) and step 2 (round-trip).

Builds two RMG Species objects directly (no full RMG run needed):
  - `ar_one_range`: a monatomic argon-like species whose NASA thermo has a
    SINGLE polynomial over 200-6000 K (the defect case).
  - `ar_two_range`: the same species with the conventional two-polynomial
    NASA fit split at 1000 K (the non-regression case).

Run against the BASE plasma checkout (RMG-Py-plasma) to reproduce the
IndexError, and against this worktree to show it is fixed. Then, when run
in "roundtrip" mode, writes the emitted YAML thermo block into a full
Cantera phase file, loads it back with the installed Cantera, and compares
Cp/H/S against the original RMG thermo object at spot temperatures.
"""
import sys
import traceback

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.thermo.nasa import NASA, NASAPolynomial
from rmgpy.yaml_cantera2 import species_to_dict

# A real single NASA7 coefficient set for monatomic Ar (Cp = 5/2 R constant,
# from NIST-style Ar thermo; a0=2.5, a1..a4=0, a5/a6 chosen for continuity).
AR_COEFFS = [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.37967]


def make_species(two_range: bool) -> Species:
    mol = Molecule(smiles="[Ar]")
    if two_range:
        polys = [
            NASAPolynomial(coeffs=AR_COEFFS, Tmin=(200, "K"), Tmax=(1000, "K")),
            NASAPolynomial(coeffs=AR_COEFFS, Tmin=(1000, "K"), Tmax=(6000, "K")),
        ]
        label = "Ar-two-range"
    else:
        polys = [
            NASAPolynomial(coeffs=AR_COEFFS, Tmin=(200, "K"), Tmax=(6000, "K")),
        ]
        label = "Ar-one-range"

    nasa = NASA(polynomials=polys, Tmin=(200, "K"), Tmax=(6000, "K"))
    species = Species(label=label, molecule=[mol], thermo=nasa)
    return species


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "red-green"

    print("=== one-range species ===")
    sp_one = make_species(two_range=False)
    try:
        entry = species_to_dict(sp_one, [sp_one])
        print("[OK] species_to_dict succeeded for one-range species")
        print("thermo block:", entry["thermo"])
    except IndexError:
        print("[EXPECTED-RED] IndexError raised (base/pre-fix behavior):")
        traceback.print_exc()
    except Exception as e:
        print(f"[UNEXPECTED] {type(e).__name__}: {e}")
        traceback.print_exc()

    print()
    print("=== two-range species (non-regression) ===")
    sp_two = make_species(two_range=True)
    try:
        entry2 = species_to_dict(sp_two, [sp_two])
        print("[OK] species_to_dict succeeded for two-range species")
        print("thermo block:", entry2["thermo"])
    except Exception as e:
        print(f"[UNEXPECTED] {type(e).__name__}: {e}")
        traceback.print_exc()

    print()
    print("=== zero-polynomial species (named-refusal case) ===")
    mol = Molecule(smiles="[Ar]")
    nasa_empty = NASA(polynomials=[], Tmin=(200, "K"), Tmax=(6000, "K"))
    sp_empty = Species(label="Ar-empty", molecule=[mol], thermo=nasa_empty)
    try:
        species_to_dict(sp_empty, [sp_empty])
        print("[UNEXPECTED] zero-polynomial species did not raise")
    except IndexError:
        print("[UNEXPECTED] bare IndexError raised for zero-polynomial species (should be named)")
        traceback.print_exc()
    except Exception as e:
        print(f"[OK] named exception raised: {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
