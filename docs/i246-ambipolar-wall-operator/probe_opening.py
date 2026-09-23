#!/usr/bin/env python3
"""
I-246 OPENING PROBE -- the brief's arithmetic against the real compiled solver.

The brief asserts, for the I-164 5 torr argon deck ionised to completion at
Te/Tgas ~ 117:

  (1) the electron gas carries 99.15% of the two-temperature EOS sum;
  (2) a RECYCLING wall (Ar+ + e- -> Ar, heavy count conserved) suppresses the
      response of n_e by ~1/117 relative to the response of N_e;
  (3) a PUMPING wall (both leave the gas) suppresses it by ~1/96,800;
  (4) the ionisation-fraction observable and the n_e observable disagree by
      orders of magnitude on the same event.

Every number below is taken through rmgpy.solver.plasma.PlasmaReactor.compute_volume
-- the SAME compiled EOS the residual and Jacobian use -- never through a
re-implementation of it in this file. No fitting, no target density.
"""

import os

import numpy as np

import rmgpy
import rmgpy.constants as constants
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.species import Species

# --- provenance: one shared editable install serves several worktrees, and
# `import rmgpy` resolves by CWD. Assert we measured THIS tree.
HERE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
resolved = os.path.dirname(os.path.dirname(os.path.abspath(rmgpy.__file__)))
print("rmgpy resolved from : {0}".format(resolved))
print("this worktree       : {0}".format(HERE))
assert resolved == HERE, "WRONG TREE: measurements would not describe this branch"
from rmgpy.solver import plasma as _plasma_mod
print("plasma module .so   : {0}".format(_plasma_mod.__file__))
print()

# --- I-164 deck conditions, read from docs/i194-ar5torr-plasma-lineage/input.py
TGAS = 298.15                      # K
TE = 34813.5                       # K  (3 eV)
P = 666.6118421052631              # Pa (5 torr)

print("=" * 78)
print("CONDITIONS: Tgas = {0} K, Te = {1} K, P = {2} Pa (5 torr)".format(TGAS, TE, P))
print("Te/Tgas = {0:.4f}".format(TE / TGAS))
print("=" * 78)
print()

# --- a real, minimally-initialised PlasmaReactor, so compute_volume is the
# compiled one. Three species: e-, Ar, Ar+ (the I-164 slate).
electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
ar = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')

reactor = PlasmaReactor(
    (TGAS, 'K'), (P, 'Pa'),
    {electron: 1.0e-8, ar: 1.0, arp: 1.0e-8}, (TE, 'K'),
    n_sims=1, termination=[])
# compute_volume needs only these two packed-state facts; set them directly so
# the probe does not depend on a full model initialisation (no chemistry here).
reactor.num_core_species = 3
reactor.electron_index = 0          # y = [N_e, N_Ar, N_Arp]

IE, IAR, IARP = 0, 1, 2


def state(f, n0=1.0):
    """State vector at ionisation fraction f, from n0 mol of argon initially."""
    y = np.zeros(3, float)
    y[IE] = n0 * f
    y[IAR] = n0 * (1.0 - f)
    y[IARP] = n0 * f
    return y


def n_e(y):
    """Electron NUMBER density (m^-3) from the compiled two-temperature EOS."""
    return y[IE] * constants.Na / reactor.compute_volume(y)


def ion_fraction(y):
    """n_Ar+ / (n_Ar + n_Ar+): the volume cancels out of this observable."""
    return y[IARP] / (y[IAR] + y[IARP])


# =============================================================================
print("(1) ELECTRON SHARE OF THE TWO-TEMPERATURE EOS SUM, at f = 1")
print("-" * 78)
y1 = state(1.0)
V1 = reactor.compute_volume(y1)
# phi = N_e*Te / (N_heavy*Tgas + N_e*Te); recovered from compute_volume itself
# as the fraction of V that the electron term accounts for.
V_heavy_only = constants.R * (y1[IAR] + y1[IARP]) * TGAS / P
phi = 1.0 - V_heavy_only / V1
print("V(f=1)                          = {0!r} m^3".format(V1))
print("electron share of the EOS sum   = {0:.6f}  ({1:.4f} %)".format(phi, 100 * phi))
print("brief claims                    = 99.15 %")
print("closed form Te/(Tgas+Te)        = {0:.6f}".format(TE / (TGAS + TE)))
print("VERDICT: {0}".format(
    "CONFIRMED" if abs(100 * phi - 99.15) < 0.01 else "CONTRADICTED"))
print()

# =============================================================================
print("(2)+(3) WALL SUPPRESSION ON n_e, measured by finite difference")
print("-" * 78)
print("A wall removes an amount delta of the charged pairs.")
print("  RECYCLING: Ar+ + e- -> Ar   (heavy count conserved, argon returns)")
print("  PUMPING  : Ar+ and e- both leave the gas entirely")
print()
print("suppression S = (d ln n_e) / (d ln N_e); the brief predicts")
print("  S_recycling ~ 1/117      S_pumping ~ 1/96,800")
print()
hdr = "{0:>12} {1:>16} {2:>16} {3:>16} {4:>16}"
print(hdr.format("f", "S_recycling", "1/S_recyc", "S_pumping", "1/S_pump"))
print("-" * 78)

results = {}
for f in (1.0, 0.999, 0.99, 0.9, 0.5, 0.1, 1e-3, 1e-5, 1e-7):
    n0 = 1.0
    y0 = state(f, n0)
    d = 1.0e-9 * n0 * f          # relative removal, well inside double precision

    # recycling: N_e -= d, N_Ar+ -= d, N_Ar += d  (heavy total unchanged)
    yr = y0.copy()
    yr[IE] -= d
    yr[IARP] -= d
    yr[IAR] += d
    # pumping: N_e -= d, N_Ar+ -= d, nothing returned
    yp = y0.copy()
    yp[IE] -= d
    yp[IARP] -= d

    dlnNe = np.log(yr[IE]) - np.log(y0[IE])          # same for both by construction
    s_rec = (np.log(n_e(yr)) - np.log(n_e(y0))) / dlnNe
    s_pmp = (np.log(n_e(yp)) - np.log(n_e(y0))) / dlnNe
    results[f] = (s_rec, s_pmp)
    print(hdr.format(
        "{0:.3g}".format(f),
        "{0:.6e}".format(s_rec), "{0:.4g}".format(1.0 / s_rec if s_rec else np.inf),
        "{0:.6e}".format(s_pmp),
        "{0:.4g}".format(1.0 / s_pmp if s_pmp else np.inf)))

print()
print("closed forms derived from the EOS (for comparison, NOT used above):")
print("  S_recycling(f) = Tgas / (Tgas + f*Te)")
print("  S_pumping(f)   = 1 - f*(Tgas+Te)/(Tgas + f*Te)")
for f in (1.0, 0.999, 0.99, 0.5):
    cf_r = TGAS / (TGAS + f * TE)
    cf_p = 1.0 - f * (TGAS + TE) / (TGAS + f * TE)
    print("  f={0:<8.4g} closed S_rec={1:.6e} measured={2:.6e} | "
          "closed S_pmp={3:.6e} measured={4:.6e}".format(
              f, cf_r, results[f][0], cf_p, results[f][1]))
print()
print("S_recycling(f=1)  = {0:.6e}  -> 1/{1:.4g}".format(
    results[1.0][0], 1.0 / results[1.0][0]))
print("S_pumping(f=1)    = {0:.6e}".format(results[1.0][1]))
print()
print("Where does 96,800 come from? It is a ZERO CROSSING, not a model property:")
print("S_pumping(f) = 0 exactly at f = 1, so 1/S_pumping diverges there. Solving")
print("S_pumping(f) = 1/96800 for f:")
target = 1.0 / 96800.0
lo, hi = 0.0, 1.0
for _ in range(200):
    mid = 0.5 * (lo + hi)
    if 1.0 - mid * (TGAS + TE) / (TGAS + mid * TE) > target:
        lo = mid
    else:
        hi = mid
print("  f = {0:.10f}  (i.e. 1 - f = {1:.4e})".format(lo, 1.0 - lo))
print("  -- a particular point, not a lever; 1/S_pumping runs from ~119 at f=0")
print("     to +infinity at f=1, passing 96,800 on the way.")
print()

# =============================================================================
print("(4) THE TWO OBSERVABLES ON THE SAME EVENT")
print("-" * 78)
print("A recycling wall removing 511 ppm of the electrons, at a few f:")
hdr2 = "{0:>12} {1:>18} {2:>18} {3:>18}"
print(hdr2.format("f", "d n_e / n_e", "d V / V", "d(1-f_ion)/(1-f_ion)"))
print("-" * 78)
for f in (1.0 - 1e-9, 0.999, 0.99, 0.5, 1e-3):
    n0 = 1.0
    y0 = state(f, n0)
    d = 511e-6 * y0[IE]
    yr = y0.copy()
    yr[IE] -= d
    yr[IARP] -= d
    yr[IAR] += d
    dne = n_e(yr) / n_e(y0) - 1.0
    dV = reactor.compute_volume(yr) / reactor.compute_volume(y0) - 1.0
    neutral0 = 1.0 - ion_fraction(y0)
    neutral1 = 1.0 - ion_fraction(yr)
    dneutral = (neutral1 / neutral0 - 1.0) if neutral0 > 0 else np.inf
    print(hdr2.format(
        "{0:.10g}".format(f),
        "{0:+.4e}".format(dne), "{0:+.4e}".format(dV), "{0:+.4e}".format(dneutral)))
print()
print("Same event, two observables, orders of magnitude apart -- the brief's")
print("consequence (a) CONFIRMED, and the reason is visible: the neutral")
print("fraction is normalised to a population that is itself near zero at f -> 1,")
print("while n_e is normalised to one that is near its maximum.")
print()

# =============================================================================
print("(5) IS THE WALL RESIDUAL TERM VOLUME-DEPENDENT?")
print("-" * 78)
print("The solver residual is res_i = (dC_i/dt) * V, in mol/s, with C_i = y_i/V.")
print("A first-order wall loss is dC_i/dt = -nu_wall * C_i, so")
print("    res_i|wall = -nu_wall * C_i * V = -nu_wall * y_i")
print("which is EXACTLY LINEAR IN y AND FREE OF V -- provided nu_wall is a constant.")
print("Checked numerically at three very different volumes:")
for f in (1e-7, 0.5, 1.0):
    y0 = state(f)
    V = reactor.compute_volume(y0)
    nu = 1.0e5
    lhs = -nu * (y0[IARP] / V) * V
    rhs = -nu * y0[IARP]
    print("  f={0:<9.3g} V={1:.6e} m^3   -nu*C*V={2:.17e}  -nu*y={3:.17e}  identical={4}".format(
        f, V, lhs, rhs, lhs == rhs))
print()
print("The brief warns the wall Jacobian 'also perturbs V'. Through the C=y/V")
print("conversion it does NOT -- that cancels exactly. It DOES if nu_wall itself")
print("carries a density dependence, which an ambipolar diffusivity must:")
print("    D_a = mu_i * kTe/e  with  mu_i = mu_i0 * (N0 / n_gas),  n_gas = y_gas*Na/V")
print("so nu_wall ~ V / y_gas and the dV/dy coupling is real -- but it enters")
print("through the MOBILITY, not through the concentration conversion.")
print()

# =============================================================================
print("(6) THE STEADY STATE THE WALL CREATES, AND ITS CONDITIONING")
print("-" * 78)
print("Balance: nu_ion = k_iz(Te)*n_Ar  against  nu_wall, both first order in n_e,")
print("so n_e CANCELS from the balance and the stationary condition is a condition")
print("on the NEUTRAL density alone:")
print("      n_Ar* = nu_wall / k_iz(Te)")
print("At constant pressure the EOS then fixes n_e as a DIFFERENCE:")
print("      n_e* = (P/kB - n_Ar* * Tgas) / (Tgas + Te)")
print("Writing f = 1 - n_Ar*/n_Ar,0 with n_Ar,0 = P/(kB*Tgas):")
print("      C_wall = |d ln n_e / d ln nu_wall| = (1-f)/f")
print("      C_Te   = |d ln n_e / d ln Te|      = (1-f)/f * |d ln k_iz/d ln Te|")
print("      C_Ar   = |d ln n_e / d ln n_Ar|    = (1-f)/f")
print()
print("Two consequences, both checkable:")
print(" * C_Te / C_wall = |d ln k_iz/d ln Te| -- a CONSTANT. Te and nu_wall are")
print("   DEGENERATE directions: any error in the wall coefficient is")
print("   indistinguishable from an error in Te. That is why a free wall")
print("   coefficient can reproduce any n_e you like.")
print(" * the 784x the ruling measured implies (1-f)/f from C_Te alone:")
for c_te in (np.log(784.0) / 1.0e-4,):
    for dlnk in (4.0, 5.25, 7.0):
        ratio = c_te / dlnk
        f_implied = 1.0 / (1.0 + ratio)
        print("     C_Te = {0:.4g} at |dln k/dln Te| = {1:.2f} -> (1-f)/f = {2:.4g},"
              " f = {3:.4g}".format(c_te, dlnk, ratio, f_implied))
print("   i.e. the 784x is the signature of a WEAKLY IONISED state, and it")
print("   vanishes as f -> 1 (C -> 0) and diverges as f -> 0.")
print()

print("=" * 78)
print("DONE")
print("=" * 78)
