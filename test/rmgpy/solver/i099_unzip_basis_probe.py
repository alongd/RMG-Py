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
import math
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


def _family_obj(db, fam):
    """`rxn.family` is a label string after find_degenerate_reactions."""
    return db.kinetics.families[fam] if isinstance(fam, str) else fam


def _family_label(fam):
    return fam if isinstance(fam, str) else fam.label


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
    by_family = {_family_label(k): v for k, v in by_family.items()}
    for fam in sorted(by_family):
        print(f"    {fam}: {len(by_family[fam])}")
    print("  product sets generated:")
    for rxn in rxns:
        print(f"    {_family_label(rxn.family):32s} -> {' + '.join(sorted(_smiles(p) for p in rxn.products))}")

    target = {_smiles(styrene), _smiles(short)}
    hits = 0
    for rxn in rxns:
        prods = {_smiles(p) for p in rxn.products}
        if prods != target:
            continue
        hits += 1
        try:
            _fam = _family_obj(db, rxn.family)
            _tmpl = (rxn.template if not isinstance(rxn.template[0], str)
                     else _fam.retrieve_template(rxn.template))
            rxn.kinetics = _fam.get_kinetics_for_template(
                _tmpl, degeneracy=rxn.degeneracy)[0]
        except Exception as exc:  # noqa: BLE001
            print(f"  family {_family_label(rxn.family)}: kinetics lookup failed "
                  f"({type(exc).__name__}: {exc})")
            continue
        print(f"  family {_family_label(rxn.family)}: {rxn.kinetics}")
        for T in TEMPERATURES:
            try:
                print(f"      k({T:7.2f} K) = "
                      f"{rxn.kinetics.get_rate_coefficient(T):.4e}")
            except Exception as exc:  # noqa: BLE001
                print(f"      k({T:7.2f} K) unavailable: {exc}")
    if hits == 0:
        print("  NO family produced this step in the SCISSION direction from "
              "the model reactant.\n  RMG's families are written in the "
              "addition direction, so ask the reverse question.")

    # Depropagation is the exact reverse of propagation, and propagation IS a
    # template of R_Addition_MultipleBond (radical + monomer -> chain radical).
    # Asking in that direction is asking RMG for its own rate rule; the
    # depropagation rate then follows from that rate and the equilibrium
    # constant built out of the thermochemistry printed above. No fit, no QM.
    print("\n--- the reverse (propagation) direction: radical + styrene ---")
    prop = []
    for family in db.kinetics.families.values():
        try:
            prop.extend(family.generate_reactions(
                [short.molecule[0], styrene.molecule[0]]))
        except Exception as exc:  # noqa: BLE001
            print(f"  (family {family.label} raised {type(exc).__name__}: {exc})")
    prop = find_degenerate_reactions(prop, set(), kinetics_database=db.kinetics)
    print(f"  positive control: {len(prop)} bimolecular reaction(s) generated")
    want = {_smiles(reactant)}
    for rxn in prop:
        if {_smiles(p) for p in rxn.products} != want:
            continue
        try:
            _fam = _family_obj(db, rxn.family)
            _tmpl = (rxn.template if not isinstance(rxn.template[0], str)
                     else _fam.retrieve_template(rxn.template))
            rxn.kinetics = _fam.get_kinetics_for_template(
                _tmpl, degeneracy=rxn.degeneracy)[0]
        except Exception as exc:  # noqa: BLE001
            print(f"  family {_family_label(rxn.family)}: kinetics lookup failed "
                  f"({type(exc).__name__}: {exc})")
            continue
        print(f"  PROPAGATION via {_family_label(rxn.family)}: {rxn.kinetics}")
        for T in TEMPERATURES:
            try:
                kf = rxn.kinetics.get_rate_coefficient(T)
                print(f"      k_propagation({T:7.2f} K) = {kf:.4e} "
                      f"(family units)")
            except Exception as exc:  # noqa: BLE001
                print(f"      k({T:7.2f} K) unavailable: {exc}")
        # Close the loop: k_depropagation = k_propagation / Kc, with Kc built
        # from the SAME thermochemistry printed above. This is the arithmetic
        # the settling route consists of -- reported to scope the route, NOT
        # adopted as a deck value (see the caveat printed after it).
        R_SI = 8.314462618
        P_STD = 1.0e5
        print("      -- closing the loop through Kc (per-chain-end s^-1) --")
        for T in TEMPERATURES:
            dG_dep = ((styrene.thermo.get_enthalpy(T) + short.thermo.get_enthalpy(T)
                       - reactant.thermo.get_enthalpy(T))
                      - T * (styrene.thermo.get_entropy(T) + short.thermo.get_entropy(T)
                             - reactant.thermo.get_entropy(T)))
            # propagation is the reverse of depropagation; dn = -1 (2 -> 1)
            Kp_prop = math.exp(dG_dep / (R_SI * T))     # = exp(-dG_prop/RT)
            Kc_prop = Kp_prop * (R_SI * T / P_STD)      # m^3/mol
            kf_si = rxn.kinetics.get_rate_coefficient(T) * 1.0e-6   # cm^3 -> m^3
            print(f"      k_depropagation({T:7.2f} K) = {kf_si / Kc_prop:.4e} s^-1"
                  f"   [Kc_prop = {Kc_prop:.4e} m^3/mol]")
        print("      CAVEAT: the rate rule above is a generic [R_R;YJ] fallback "
              "at Euclidean\n      distance 7.81 with E0 = 0.5 kcal/mol and "
              "A = 1e13 cm^3/(mol s). Literature\n      styrene propagation is "
              "nearer A ~ 1e10 cm^3/(mol s), Ea ~ 2 kcal/mol, so this\n      "
              "number scopes the ROUTE and must not be adopted as a deck value.")


if __name__ == "__main__":
    main()
