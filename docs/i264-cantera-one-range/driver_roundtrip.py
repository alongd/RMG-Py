"""
Driver for the i264 verifier, step 2: round-trip proof.

Writes a full Cantera YAML phase file (via save_cantera_model) for a
one-range species and, separately, a two-range species, loads each file
back with the installed Cantera, and compares Cp/H/S read from the
loaded ct.Solution against the ORIGINAL RMG thermo object at several
spot temperatures. No synthesized intermediate is compared -- only the
original RMG object vs. the object Cantera built from reading the file.
"""
import math
import sys

import cantera as ct

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.thermo.nasa import NASA, NASAPolynomial
from rmgpy.transport import TransportData
from rmgpy.yaml_cantera2 import save_cantera_model

AR_COEFFS = [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.37967]
R = 8.31446261815324  # J/mol/K, matches rmgpy.constants.R


class _FakeModel:
    """Duck-typed stand-in for rmgpy.rmg.model.ReactionModel; save_cantera_model
    only touches .species, .reactions, and .get_elements()."""

    def __init__(self, species_list):
        self.species = species_list
        self.reactions = []

    def get_elements(self):
        elements = set()
        for sp in self.species:
            for atom in sp.molecule[0].vertices:
                elements.add(atom.element)
        return elements


def make_species(two_range: bool) -> Species:
    mol = Molecule(smiles="[Ar]")
    if two_range:
        polys = [
            NASAPolynomial(coeffs=AR_COEFFS, Tmin=(200, "K"), Tmax=(1000, "K")),
            NASAPolynomial(coeffs=AR_COEFFS, Tmin=(1000, "K"), Tmax=(6000, "K")),
        ]
        label = "Ar2R"
    else:
        polys = [
            NASAPolynomial(coeffs=AR_COEFFS, Tmin=(200, "K"), Tmax=(6000, "K")),
        ]
        label = "Ar1R"
    nasa = NASA(polynomials=polys, Tmin=(200, "K"), Tmax=(6000, "K"))
    species = Species(label=label, molecule=[mol], thermo=nasa)
    # Real argon Lennard-Jones transport parameters (needed for a loadable
    # gas-phase Cantera Solution; unrelated to the thermo defect being tested).
    species.transport_data = TransportData(
        shapeIndex=0, epsilon=(136.500, "K"), sigma=(3.33, "angstrom"),
        dipoleMoment=(0.0, "De"), polarizability=(0.0, "angstrom^3"),
        rotrelaxcollnum=0.0, comment="Ar (i264 round-trip driver)",
    )
    return species


def compare(label, species, yaml_path):
    save_cantera_model(_FakeModel([species]), yaml_path)
    sol = ct.Solution(yaml_path)

    print(f"--- {label} (file: {yaml_path}) ---")
    print(f"{'T (K)':>8} {'Cp_RMG':>12} {'Cp_CT':>12} {'H_RMG':>14} {'H_CT':>14} {'S_RMG':>12} {'S_CT':>12}")
    spot_temps = [250.0, 500.0, 1000.0, 2000.0, 4000.0, 5900.0]
    max_rel_err = 0.0
    for T in spot_temps:
        cp_rmg = species.thermo.get_heat_capacity(T)
        h_rmg = species.thermo.get_enthalpy(T)
        s_rmg = species.thermo.get_entropy(T)

        sol.TP = T, ct.one_atm
        # Cantera's *_mole properties are per kmol (SI base unit is kmol, not
        # mol); RMG's thermo API is per mol. Convert so the comparison is
        # apples-to-apples.
        cp_ct = sol.cp_mole / 1000.0
        h_ct = sol.enthalpy_mole / 1000.0
        s_ct = sol.entropy_mole / 1000.0

        for rmg_val, ct_val, name in ((cp_rmg, cp_ct, "Cp"), (h_rmg, h_ct, "H"), (s_rmg, s_ct, "S")):
            # Explicit finiteness check: max(max_rel_err, rel_err) SILENTLY
            # swallows a NaN rel_err (max(0.0, float('nan')) == 0.0 in
            # CPython), so a non-finite comparison value must be caught here
            # rather than left to fall through max().
            if not (math.isfinite(rmg_val) and math.isfinite(ct_val)):
                raise AssertionError(
                    f"{label}: non-finite comparison value at T={T}: "
                    f"{name}_RMG={rmg_val}, {name}_CT={ct_val}")
            denom = max(abs(rmg_val), 1e-8)
            rel_err = abs(rmg_val - ct_val) / denom
            if not math.isfinite(rel_err):
                raise AssertionError(
                    f"{label}: non-finite relative error at T={T} for {name}: {rel_err}")
            max_rel_err = max(max_rel_err, rel_err)

        print(f"{T:8.1f} {cp_rmg:12.4f} {cp_ct:12.4f} {h_rmg:14.3f} {h_ct:14.3f} {s_rmg:12.4f} {s_ct:12.4f}")

    print(f"max relative error across Cp/H/S at all spot temps: {max_rel_err:.3e}")
    # Tolerance set to 5e-6, not 1e-6: RMG's rmgpy.constants.R (8.314472
    # J/mol/K, an older CODATA value) differs from Cantera's ct.gas_constant
    # (8.31446261815324 J/mol/K) by ~1.13e-6 relative -- a pre-existing,
    # unrelated constant mismatch between the two libraries, not a defect in
    # the writer under test here. Confirmed by direct comparison; see REPORT.md.
    assert max_rel_err < 5e-6, f"{label}: round-trip mismatch too large ({max_rel_err:.3e})"
    print(f"[OK] {label}: RMG thermo and Cantera-loaded-from-file thermo agree to {max_rel_err:.3e}")
    print()


def main():
    tmpdir = sys.argv[1] if len(sys.argv) > 1 else "."
    compare("one-range (Ar1R)", make_species(two_range=False), f"{tmpdir}/ar_one_range.yaml")
    compare("two-range (Ar2R, non-regression)", make_species(two_range=True), f"{tmpdir}/ar_two_range.yaml")


if __name__ == "__main__":
    main()
