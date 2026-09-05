#!/usr/bin/env python
"""i099: what, if anything, could settle the unzip barrier without a fit or a QM job.

`k_unzip = 100.0` is declared in every polymer deck since poly_101 with no cited
basis anywhere in the project.  Its temperature-dependent successor in the
solver is `k_depropagation` (an {A, n, Ea} triplet evaluated at runtime T,
rmgpy/solver/polymer.pyx:7764-7768), so "what should k_unzip be" is really
"what is the chain-end depropagation barrier".

This probe does NOT answer that.  It measures the one quantity a Blowers-Masel /
Evans-Polanyi barrier would be built from -- the reaction enthalpy of the
elementary depropagation step, computed from RMG's own thermochemistry -- and
then asks RMG's kinetics families for their own estimate of that step.  Both are
no-fit and CPU-only.  The point is to establish whether the settling route
exists, not to set a deck value from it.

Model reaction (the elementary beta-scission a PS chain-end radical undergoes):

    Ph-CH(.)-CH2-CH(Ph)-CH3  ->  CH2=CH-Ph  +  CH3-CH(.)-Ph
    1,3-diphenylbut-1-yl          styrene      1-phenylethyl

The radical is on C1; beta-scission cuts the C2-C3 bond, releasing one styrene
unit and leaving a shorter chain-end radical.  That is exactly one unzip event.

Run:  python test/rmgpy/solver/i099_unzip_basis_probe.py
"""
import os
import sys

import rmgpy

assert os.path.realpath(rmgpy.__file__).startswith(
    "/home/alon/Code/RMG-Py-i099-arrhenius-empty"), rmgpy.__file__

import rmgpy.solver.base as _b  # noqa: E402
print("SOLVER BACKEND (runtime MRO):",
      [c.__module__ + "." + c.__name__ for c in _b.ReactionSystem.__mro__])

from rmgpy import settings  # noqa: E402
from rmgpy.data.rmg import RMGDatabase  # noqa: E402
from rmgpy.species import Species  # noqa: E402

CHAIN_END_RADICAL = "[CH](c1ccccc1)CC(C)c1ccccc1"   # 1,3-diphenylbut-1-yl
STYRENE = "C=Cc1ccccc1"
SHORT_END_RADICAL = "C[CH]c1ccccc1"                  # 1-phenylethyl

# Temperatures the polymer decks actually run at.
TEMPERATURES = (298.15, 700.0, 900.0, 1100.0)


def _smiles(obj):
    """SMILES of a Molecule or of a Species' first resonance structure."""
    mol = obj if hasattr(obj, "to_smiles") else obj.molecule[0]
    return mol.to_smiles()


def main():
    db = RMGDatabase()
    db.load_thermo(os.path.join(settings["database.directory"], "thermo"))

    reactant = Species().from_smiles(CHAIN_END_RADICAL)
    styrene = Species().from_smiles(STYRENE)
    short = Species().from_smiles(SHORT_END_RADICAL)

    for spc in (reactant, styrene, short):
        spc.thermo = db.thermo.get_thermo_data(spc)

    print("\n--- thermochemistry (RMG group additivity, no fit, no QM) ---")
    for name, spc in (("chain-end radical", reactant),
                      ("styrene", styrene),
                      ("shorter end radical", short)):
        print(f"  {name:22s} {spc.molecule[0].to_smiles():34s} "
              f"H298 = {spc.thermo.get_enthalpy(298.15) / 1000.0:9.2f} kJ/mol   "
              f"source = {spc.thermo.comment[:60]}")

    print("\n--- depropagation (one unzip event) dHrxn(T) ---")
    for T in TEMPERATURES:
        dH = (styrene.thermo.get_enthalpy(T) + short.thermo.get_enthalpy(T)
              - reactant.thermo.get_enthalpy(T))
        dS = (styrene.thermo.get_entropy(T) + short.thermo.get_entropy(T)
              - reactant.thermo.get_entropy(T))
        print(f"  T = {T:7.2f} K   dHrxn = {dH / 1000.0:8.2f} kJ/mol   "
              f"dSrxn = {dS:7.2f} J/mol/K   dGrxn = {(dH - T * dS) / 1000.0:8.2f} kJ/mol")

    print("\nNOTE: dHrxn is a LOWER BOUND on the depropagation barrier for an "
          "endothermic step;\n      it is NOT the barrier. Converting it to one "
          "needs a family's own\n      Blowers-Masel/Evans-Polanyi parameters -- "
          "which is what the kinetics\n      leg below asks RMG for.")

    if "--kinetics" not in sys.argv:
        print("\n(skipping the kinetics-family leg; pass --kinetics to run it, "
              "it loads the\n full kinetics database and takes several minutes)")
        return

    print("\n--- RMG's own family estimate for the same elementary step ---")
    db.load_kinetics(os.path.join(settings["database.directory"], "kinetics"),
                     kinetics_families="default")
    from rmgpy.data.kinetics.common import find_degenerate_reactions  # noqa: E402
    rxns = []
    for family in db.kinetics.families.values():
        try:
            rxns.extend(family.generate_reactions([reactant.molecule[0]]))
        except Exception as exc:  # noqa: BLE001
            print(f"  (family {family.label} raised {type(exc).__name__}: {exc})")
    rxns = find_degenerate_reactions(rxns, set(), kinetics_database=db.kinetics)

    # POSITIVE CONTROL. "No family produced the step" is only a finding if
    # generation worked at all; if this count is 0 the negative result is a
    # probe artifact and must be reported as such, not as a fact about RMG.
    print(f"  positive control: {len(rxns)} unimolecular reaction(s) generated "
          f"from the model reactant across all families")
    by_family = {}
    for rxn in rxns:
        by_family.setdefault(rxn.family, []).append(rxn)
    for fam in sorted(by_family):
        print(f"    {fam}: {len(by_family[fam])}")
    print("  product sets generated:")
    for rxn in rxns:
        print(f"    {rxn.family:32s} -> {' + '.join(sorted(_smiles(p) for p in rxn.products))}")

    target = {_smiles(styrene), _smiles(short)}
    hits = 0
    for rxn in rxns:
        prods = {_smiles(p) for p in rxn.products}
        if prods != target:
            continue
        hits += 1
        try:
            rxn.kinetics = rxn.family.get_kinetics_for_template(
                rxn.template, degeneracy=rxn.degeneracy)[0]
        except Exception as exc:  # noqa: BLE001
            print(f"  family {rxn.family.label}: kinetics lookup failed "
                  f"({type(exc).__name__}: {exc})")
            continue
        print(f"  family {rxn.family.label}: {rxn.kinetics}")
        for T in TEMPERATURES:
            try:
                print(f"      k({T:7.2f} K) = "
                      f"{rxn.kinetics.get_rate_coefficient(T):.4e}")
            except Exception as exc:  # noqa: BLE001
                print(f"      k({T:7.2f} K) unavailable: {exc}")
    if hits == 0:
        print("  NO family produced this elementary step from the model "
              "reactant. That is itself\n  the finding: the settling route is "
              "not available off the shelf here.")


if __name__ == "__main__":
    main()
