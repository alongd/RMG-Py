#!/usr/bin/env python
"""
Kill-check for the i067 brief's central premise:

    "generatePolymerConstraints is not set, so polymer proxy chemistry is
     unbounded, and the polymer-family reaction-generation leg has been
     running unbounded."

The premise is load-bearing: if proxy species are ALREADY bounded by the gas
tier, then adding a generatePolymerConstraints block does not *bound an
unbounded leg* -- it REROUTES the >=polymerSizeThreshold band out of the gas
tier and into a new one, which can loosen as easily as tighten.

Run under an activated rmg_env, from the i067 worktree.
"""
import sys, os
sys.path.insert(0, "/home/alon/Code/RMG-Py-i067-enlarge-bound")

from rmgpy.molecule import Molecule
from rmgpy.species import Species
import rmgpy.rmg.input as rmg_input
from rmgpy.rmg.main import RMG
from rmgpy.constraints import (
    fails_species_constraints,
    is_polymer_constraint_member,
    DEFAULT_POLYMER_SIZE_THRESHOLD,
)

print("rmgpy imported from:", os.path.dirname(sys.modules["rmgpy"].__file__))
print("DEFAULT_POLYMER_SIZE_THRESHOLD =", DEFAULT_POLYMER_SIZE_THRESHOLD)
print()

# --- reproduce the poly_105 deck's gas tier verbatim -------------------------
rmg = RMG()
rmg_input.rmg = rmg
rmg.species_constraints = dict(
    allowed=["input species", "seed mechanisms", "reaction libraries"],
    maximumCarbonAtoms=30,
    maximumOxygenAtoms=5,
    maximumNitrogenAtoms=0,
    maximumSiliconAtoms=0,
    maximumSulfurAtoms=0,
    maximumHeavyAtoms=35,
    maximumRadicalElectrons=2,
    maximumSingletCarbenes=0,
    maximumCarbeneRadicals=0,
    allowSingletO2=True,
)
rmg.polymer_constraints = None          # <-- exactly the poly_105 deck

MONOMER = "Cc1ccc(C)c(C)c1O"            # a novolac-ish repeat surrogate
# real proxy species lifted verbatim from poly_105/RMG.log (edge, species 80)
TRIMER_PROXY = "Cc1ccc(CCc2c(C)ccc(CCc3c(O)[c]ccc3C)c2O)c(O)c1C"
# a synthetic OVERSIZE proxy: five linked units, well past the gas tier
PENTAMER = "Cc1ccc(CCc2c(C)ccc(CCc3c(C)ccc(CCc4c(C)ccc(CCc5c(C)ccc(C)c5O)c4O)c3O)c2O)c(O)c1C"

def describe(smiles, tag, proxy_flag):
    m = Molecule(smiles=smiles)
    m.is_polymer_proxy = proxy_flag
    heavy = m.get_num_atoms() - m.get_num_atoms("H")
    c = m.get_num_atoms("C")
    member = is_polymer_constraint_member(m, rmg.polymer_constraints)
    reason = fails_species_constraints(m)
    print(f"{tag:24s} C={c:3d} heavy={heavy:3d} is_polymer_proxy={proxy_flag!s:5s} "
          f"routes_to_polymer_tier={member!s:5s} rejected={bool(reason)!s:5s} {reason or ''}")
    return heavy, member, bool(reason)

print("--- A. polymer_constraints = None (the poly_105 deck as it ships) ---")
describe(MONOMER,      "monomer surrogate",   False)
h_t, m_t, r_t = describe(TRIMER_PROXY,  "trimer proxy (poly_105)", True)
h_p, m_p, r_p = describe(PENTAMER,      "pentamer proxy",          True)
print()

fail = 0
if m_p:
    print("PREMISE HOLDS: an oversize proxy routes to the (absent) polymer tier and bypasses.")
else:
    print("PREMISE INVERTED: the oversize proxy does NOT route to the polymer tier; "
          "with polymer_constraints=None it falls through to the GAS tier.")
if r_p:
    print(f"  -> and the gas tier REJECTS it ({h_p} heavy > 35). The leg is already bounded.")
else:
    print(f"  -> and it is NOT rejected. The leg really is unbounded.")
    fail = 1
print()

print("--- B. add a generatePolymerConstraints block IDENTICAL to the gas tier ---")
rmg.polymer_constraints = dict(rmg.species_constraints)
rmg.polymer_constraints.pop("allowed", None)
h2, m2, r2 = describe(PENTAMER, "pentamer proxy", True)
h3, m3, r3 = describe(TRIMER_PROXY, "trimer proxy", True)
print("  size routing is now ON: every species with heavy >= "
      f"{DEFAULT_POLYMER_SIZE_THRESHOLD} left the gas tier for the polymer tier.")
print()

print("--- C. what a species in the 15..35-heavy band loses/gains on the reroute ---")
# a species the GAS tier rejects on maximumCarbonAtoms but a loose polymer tier would not
C32 = "C" * 32
rmg.polymer_constraints = None
print("   polymer_constraints=None:")
describe(C32, "C32 n-alkane", False)
rmg.polymer_constraints = {"maximumHeavyAtoms": 60}   # a 'bound' that is looser than the gas tier
print("   polymer_constraints={'maximumHeavyAtoms': 60}  (a legal, validation-passing block):")
h4, m4, r4 = describe(C32, "C32 n-alkane", False)
if not r4:
    print("  -> ADDING a polymer block ADMITTED a species the gas tier refused. "
         "A 'bound' block can enlarge the generated space.")
else:
    print("  -> still refused.")

sys.exit(fail)
