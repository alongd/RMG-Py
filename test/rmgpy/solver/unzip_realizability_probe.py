#!/usr/bin/env python3
"""
Standalone probe for the legacy k_unzip realizability hole.

Constructs a single polymer pool directly -- no database, no reaction
generation -- and integrates the REAL solver residual (rmgpy.solver.polymer
HybridPolymerSystem.residual) as an ODE right-hand side. The full RMG
generation that first exposed this takes ~11 minutes to reach the failure;
this reaches it in milliseconds.

The instrument is deliberately the RHS, not DASSL: the defect is that the
vector field points OUT of the realizable cone (mu1 >= mu0 >= 0 for any
chain-length distribution with k >= 1). No integrator can repair a vector
field that leaves the feasible set.

Usage:
    python test/rmgpy/solver/unzip_realizability_probe.py
"""

import contextlib
import io
import sys

import numpy as np
from scipy.integrate import solve_ivp

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig


MU0, MU1, MU2, MONOMER = 1, 2, 3, 4
INITIAL_MOMENTS = (1.0, 5.0, 30.0)


def _spc(smiles, label):
    s = Species(molecule=[Molecule().from_smiles(smiles)])
    s.label = label
    return s


def build_system(k_unzip, moments=INITIAL_MOMENTS, k_scission=0.0, V_poly=1.0):
    """One pool, one gas monomer target, one inert. V_poly = 1 so moles ==
    concentration and every moment quantity below is directly comparable."""
    inert = _spc("N#N", "N2")
    mu0_s = _spc("CO", "poly_mu0")
    mu1_s = _spc("C=O", "poly_mu1")
    mu2_s = _spc("C#N", "poly_mu2")
    monomer = _spc("C", "M")
    core = [inert, mu0_s, mu1_s, mu2_s, monomer]
    mask = np.array([True, False, False, False, True], dtype=bool)

    pool = PolymerPoolConfig(
        label="poly", xs=2, explicit_dp_to_species_index={},
        mu_indices=(MU0, MU1, MU2), monomer_poly_index=MONOMER,
        k_scission=k_scission, k_unzip=k_unzip, tail_kinetics=None,
    )
    rs = HybridPolymerSystem(
        T=800.0, P=1.0e5, initial_mole_fractions={inert: 1.0}, V_poly=V_poly,
        polymer_pools=[pool], mass_transfer=[], gas_species_mask=mask,
        constant_gas_volume=False,
        initial_polymer_moments={"poly": moments}, termination=[],
    )
    # initialize_model prints the multi-page POLYMER SOLVER DIAGNOSTIC banner
    # to stdout; swallow it so this probe's own numbers stay readable.
    with contextlib.redirect_stdout(io.StringIO()):
        rs.initialize_model(core, [], [], [])
    return rs


def integrate(k_unzip, t_end, moments=INITIAL_MOMENTS, rtol=1e-10, atol=1e-14):
    """Integrate the real residual. Returns (solution, first cone-exit time)."""
    rs = build_system(k_unzip, moments=moments)
    y0 = rs.y.copy()
    y0[MU0], y0[MU1], y0[MU2] = moments
    zeros = np.zeros_like(y0)

    def rhs(t, y):
        return rs.residual(t, y, zeros)[0]

    def cone_event(t, y):
        # mu1 - mu0: positive inside the cone, 0 at the all-monomer boundary,
        # negative outside it. Scaled by mu0(0) so the event tolerance is not
        # chased down to zero as the pool drains.
        return y[MU1] - y[MU0]

    cone_event.terminal = False
    cone_event.direction = -1

    sol = solve_ivp(rhs, (0.0, t_end), y0, method="LSODA",
                    rtol=rtol, atol=atol, events=cone_event,
                    max_step=t_end / 500.0)
    crossing = float(sol.t_events[0][0]) if len(sol.t_events[0]) else None
    return sol, crossing


def summarize(k_unzip, t_end, moments=INITIAL_MOMENTS):
    sol, crossing = integrate(k_unzip, t_end, moments=moments)
    mu0 = sol.y[MU0]
    mu1 = sol.y[MU1]
    gap = mu1 - mu0
    # Mass ledger: repeat units drained out of mu1 must equal gas monomer gained.
    drained = moments[1] - mu1[-1]
    released = sol.y[MONOMER][-1] - sol.y[MONOMER][0]
    return {
        "k_unzip": k_unzip,
        "t_end": t_end,
        "status": sol.status,
        "crossing": crossing,
        "min_gap": float(gap.min()),
        "min_mu1": float(mu1.min()),
        "min_mu0": float(mu0.min()),
        "final_mu0": float(mu0[-1]),
        "final_mu1": float(mu1[-1]),
        "drained_mu1": float(drained),
        "released_monomer": float(released),
        "mass_defect": float(drained - released),
    }


SWEEP = ((0.05, 200.0), (1.0, 20.0), (100.0, 0.2))

# The poly_104-like arm: k_scission runs alongside k_unzip, so mu0 is no
# longer constant and the (mu0, mu1) subsystem is the 2x2 linear ODE
#   d/dt [mu0, mu1] = [[-k_s, k_s], [-k_u, 0]] [mu0, mu1].
# Neither row reads mu2 or mu3, so the crossing times below are independent
# of the mu2 initial condition. Cone exit (mu1 < mu0) is structural in k_u;
# mu1 < 0 is NOT -- it needs k_u above k_s/4, where the eigenvalues go
# complex and the trajectory spirals through zero.
SCISSION_ARM_KS = 1.0
SCISSION_ARM_IC = (0.099, 0.7378, 8.0)
SCISSION_ARM_SWEEP = ((0.05, 20.0), (1.0, 20.0), (100.0, 0.5))


def scission_arm():
    """Reproduce the coupled scission+unzip case at k_scission = 1.0."""
    print()
    print("=" * 78)
    print(f"coupled arm: k_scission={SCISSION_ARM_KS}, "
          f"(mu0, mu1, mu2) = {SCISSION_ARM_IC}")
    print("=" * 78)
    out = []
    for k, t_end in SCISSION_ARM_SWEEP:
        rs = build_system(k, moments=SCISSION_ARM_IC, k_scission=SCISSION_ARM_KS)
        y0 = rs.y.copy()
        y0[MU0], y0[MU1], y0[MU2] = SCISSION_ARM_IC
        zeros = np.zeros_like(y0)

        def rhs(t, y):
            return rs.residual(t, y, zeros)[0]

        def cone(t, y):
            return y[MU1] - y[MU0]

        def zero_mu1(t, y):
            return y[MU1]

        cone.terminal = False
        cone.direction = -1
        zero_mu1.terminal = False
        zero_mu1.direction = -1
        sol = solve_ivp(rhs, (0.0, t_end), y0, method="LSODA",
                        rtol=1e-10, atol=1e-14, events=(cone, zero_mu1),
                        max_step=t_end / 2000.0)
        t_cone = float(sol.t_events[0][0]) if len(sol.t_events[0]) else None
        t_neg = float(sol.t_events[1][0]) if len(sol.t_events[1]) else None
        out.append((k, t_cone, t_neg, float(sol.y[MU1].min())))
        print(f"k_unzip={k:<8g} k_unzip/k_scission={k / SCISSION_ARM_KS:<8g} "
              f"cone t={('%.4g' % t_cone) if t_cone else 'none':<10} "
              f"mu1<0 t={('%.4g' % t_neg) if t_neg else 'none':<10} "
              f"min(mu1)={sol.y[MU1].min():+.4e}")
    print(f"  (k_scission/4 = {SCISSION_ARM_KS / 4.0:g} is the mu1 < 0 "
          f"threshold; cone exit has no threshold.)")
    return out


def main():
    print("=" * 78)
    print("probe: legacy k_unzip realizability, require mu1 >= mu0 >= 0")
    print(f"initial moments (mu0, mu1, mu2) = {INITIAL_MOMENTS}  [mean DP = 5]")
    print("=" * 78)
    rows = []
    for k, t_end in SWEEP:
        r = summarize(k, t_end)
        rows.append(r)
        cross = f"{r['crossing']:.6g}" if r["crossing"] is not None else "none"
        print(f"k_unzip={k:<8g} t_end={t_end:<7g} "
              f"mu1->mu0 crossing t={cross:<12} "
              f"min(mu1-mu0)={r['min_gap']:+.6e}  min(mu1)={r['min_mu1']:+.6e}")
        print(f"{'':>28}mass: drained mu1={r['drained_mu1']:.9e}  "
              f"released={r['released_monomer']:.9e}  "
              f"defect={r['mass_defect']:+.3e}")

    scission_arm()

    print("-" * 78)
    crossed = [r for r in rows if r["crossing"] is not None]
    if crossed:
        print(f"VERDICT: cone exit at {len(crossed)}/{len(rows)} sampled k_unzip values.")
        for r in crossed:
            k, t = r["k_unzip"], r["crossing"]
            print(f"  k={k:<8g} t_cross={t:.6g}   k*t_cross={k * t:.6g}")
        print("  A constant k*t_cross means the crossing time scales as "
              "1/k_unzip: STRUCTURAL, no threshold.")
        return 1
    print("VERDICT: no cone exit at any sampled k_unzip.")
    for r in rows:
        print(f"  k={r['k_unzip']:<8g} min(mu1-mu0)={r['min_gap']:+.6e}  "
              f"min(mu0)={r['min_mu0']:+.6e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
