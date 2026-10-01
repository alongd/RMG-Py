"""One arm of the i075 before/after harness for the i061 moment error-weight
floor. Run twice, in separate processes, and diff the two JSON dumps:

    python i075_ab_arm.py 100.0 on.json     # floor ACTIVE (shipped behaviour)
    python i075_ab_arm.py 0.0   off.json    # floor INERT  (pre-fix behaviour)

The "off" arm is not a second build. MOMENT_EWT_FLOOR_K is a module-level
constant in rmgpy/solver/polymer.pyx, so it is a plain module global that the
compiled code looks up at call time; setting it to 0.0 makes the floor
comparison `atol_array[i] < 0.0` false for every slot, i.e. NOTHING is floored
and atol_array reaches DASPK exactly as it did before the fix. The arm asserts
that (floored == [] with K=0) rather than assuming it, so a Cython change that
made the constant unpatchable would fail here instead of silently turning the
A/B into a comparison of the same thing with itself.

Everything numeric is dumped as float.hex() so the consumer can assert BITWISE
equality where the original i061 closure claimed bit-identity, rather than
approximate agreement.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TREE = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, TREE)
sys.path.insert(0, HERE)

import numpy as np  # noqa: E402
from rmgpy.reaction import Reaction  # noqa: E402
import rmgpy.solver.polymer as polymer_mod  # noqa: E402
from rmgpy.solver.polymer import (HybridPolymerSystem,  # noqa: E402
                                  PolymerPoolConfig)

FLOOR_K = float(sys.argv[1])
OUTFILE = sys.argv[2]
polymer_mod.MOMENT_EWT_FLOOR_K = FLOOR_K

from solverPolymerTest import (_two_pool_species, _KIN,  # noqa: E402
                               TestRtolNearFloorConviction)

ATOL = 1.0e-14
RTOL = 1.0e-4
result = {"floor_k": FLOOR_K,
          "polymer_module": polymer_mod.__file__,
          "tree": TREE}


def _hex(a):
    return [float(v).hex() for v in np.asarray(a, dtype=float).ravel()]


def _pool(lbl, mi):
    return PolymerPoolConfig(label=lbl, xs=2, explicit_dp_to_species_index={},
                             mu_indices=mi, monomer_poly_index=None,
                             k_scission=0.0, k_unzip=0.0, tail_kinetics=None)


def _healthy_system():
    """Two pools whose moments stay far above the floor for the whole run, so
    the floor is INERT on the physics and the trajectory must not move."""
    sp, core, mask = _two_pool_species()
    feed = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
    feed.polymer_flux_archetype = 2
    foldback = Reaction(reactants=[sp["B"]], products=[sp["B"], sp["G"]],
                        **_KIN)
    rs = HybridPolymerSystem(
        T=800.0, P=1.0e5, initial_mole_fractions={core[8]: 0.0}, V_poly=1.0,
        polymer_pools=[_pool("A", (1, 2, 3)), _pool("B", (5, 6, 7))],
        mass_transfer=[], gas_species_mask=mask.copy(),
        constant_gas_volume=False,
        initial_polymer_moments={"A": (1.0, 5.0, 30.0),
                                 "B": (0.8, 4.0, 24.0)},
        termination=[], allow_default_prospective_edge=True,
        allow_unstamped_proxy_rows=True)
    rs.initialize_model(core, [feed, foldback], [], [], atol=ATOL, rtol=RTOL)
    return rs, core


# ---- scope census + every consumer that must stay bitwise unchanged --------
rs, core = _healthy_system()
atol_arr = np.asarray(rs.atol_array, dtype=float)
chem_arr = np.asarray(rs._chem_atol_array, dtype=float)
inc = np.asarray(rs._char_rate_include_mask, dtype=bool)
mu_idx = sorted({int(i) for p in rs.polymer_pools for i in p.mu_indices})

floored = [i for i in range(rs.neq) if atol_arr[i] != chem_arr[i]]
result["neq"] = int(rs.neq)
result["num_core_species"] = int(rs.num_core_species)
result["mu_indices"] = mu_idx
result["include_mask"] = [bool(v) for v in inc]
result["floored_slots"] = floored
result["floored_values"] = _hex([atol_arr[i] for i in floored])
result["atol_array"] = _hex(atol_arr)
result["chem_atol_array"] = _hex(chem_arr)
result["pool_mu_floors"] = _hex(getattr(rs, "_pool_mu_floors", []))
result["softclamp_lam"] = float(getattr(rs, "_softclamp_lam", float("nan"))).hex()
jac = getattr(rs, "_jac_wt_atol", None)
result["jac_wt_atol"] = None if jac is None else _hex(jac)

# ---- healthy trajectory ----------------------------------------------------
traj = []
for t in (1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0):
    rs.advance(t)
    traj.append(_hex(np.asarray(rs.y)))
result["healthy_traj"] = traj

# ---- the scalar-init near-floor canary ------------------------------------
# The same fixture the i061 closure called bit-identical, driven directly
# instead of scraped out of a pytest failure message, so the comparison is on
# the numbers themselves.
canary = TestRtolNearFloorConviction()
crs = canary._fixture(rtol=1.0e-4)
crs.advance(100.0)
cy = np.asarray(crs.y)
cdn = np.asarray(crs.residual(100.0, cy.copy(), np.zeros_like(cy))[0])
result["canary_y"] = _hex(cy)
result["canary_dn"] = _hex(cdn)
result["canary_y2"] = float(cy[2]).hex()
result["canary_dn2"] = float(cdn[2]).hex()

with open(OUTFILE, "w") as fh:
    json.dump(result, fh, indent=1)

print("arm K=%r -> %s" % (FLOOR_K, OUTFILE))
print("  polymer module : %s" % polymer_mod.__file__)
print("  floored slots  : %s" % floored)
print("  canary y[2]    : %r" % float(cy[2]))
