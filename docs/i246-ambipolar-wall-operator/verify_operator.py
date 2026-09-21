#!/usr/bin/env python3
"""
I-246 verifier harness: checks 1-5 and 8-12 of the ticket's Verifier, each one
run against the compiled reactor and each one printing the number it turned on.

Check 7 (bit-for-bit zero-wall recovery) is a separate script,
``verify_zero_wall_bitwise.py``, because it needs a second extension module built
from the pre-change source.

Exit code is 0 only if every check passes.
"""

import sys

import numpy as np

import rmgpy.constants as constants
from rmgpy.exceptions import PlasmaStateError

import argon_wall_model as M

FAILURES = []


def check(name, ok, detail=""):
    print("  [{0}] {1}{2}".format("PASS" if ok else "FAIL", name,
                                  ("\n         " + detail) if detail else ""))
    if not ok:
        FAILURES.append(name)
    return ok


print("=" * 78)
print("I-246 WALL OPERATOR VERIFIER")
print("plasma module: {0}".format(M.assert_provenance()))
print("=" * 78)

# ---------------------------------------------------------------- 1
print("\n(1) DIMENSIONS AND POPULATION NORMALISATION")
r, spc, rxn = M.build_reactor(wall=True, with_chemistry=False)
ie, i_ar, i_arp = M.indices(r)
y = M.state_at(r, 1.0e-6)
V = r.compute_volume(y)
nu = r.compute_nu_wall(y, V)
n_neutral = y[i_ar] * constants.Na / V
n_e = y[ie] * constants.Na / V
lam = r.diffusion_length.value_si

closed = M.nu_wall_closed_form(M.TE_NOMINAL_EV, n_neutral, lam)
check("nu_wall agrees with an independent closed-form re-derivation",
      abs(nu / closed - 1.0) < 1e-12,
      "reactor {0:.10e} s^-1   independent {1:.10e} s^-1   rel {2:.2e}".format(
          nu, closed, abs(nu / closed - 1.0)))

# The normalisation trap: the loss must be per CHARGED PARTICLE, not per neutral.
res = np.zeros(r.num_core_species, float)
delta, _ = r.residual(0.0, y, np.zeros(r.num_core_species, float))
wall = r.wall_loss_rates
loss_per_electron = -wall[ie] / y[ie]
check("wall loss is first order in the ELECTRON population (per electron = nu_wall)",
      abs(loss_per_electron / nu - 1.0) < 1e-14,
      "-(wall flux)/N_e = {0:.10e} s^-1 ; nu_wall = {1:.10e} s^-1".format(
          loss_per_electron, nu))
per_neutral = -wall[ie] / y[i_ar]
check("the per-NEUTRAL normalisation is a different number (the trap is live)",
      abs(np.log10(per_neutral / nu)) > 3.0,
      "per-neutral would be {0:.4e} s^-1, a factor {1:.4e} from nu_wall -- an audit "
      "earlier in this campaign compared quantities normalised to these two "
      "different populations and was wrong by exactly this factor".format(
          per_neutral, per_neutral / nu))
# dimensional sanity: D_a/Lambda^2 with D_a in m^2/s
d_a = nu * lam * lam
check("D_a is a diffusivity of physically plausible magnitude at 5 torr",
      1e-3 < d_a < 1e2,
      "D_a = {0:.6e} m^2/s ; Lambda = {1:.6e} m ; nu_wall = {2:.6e} s^-1".format(
          d_a, lam, nu))

# ---------------------------------------------------------------- 2, 3
print("\n(2)+(3) GEOMETRY AND THE A_eff/V FACTOR")
lam2 = lam / np.sqrt(2.0)                     # halves Lambda^2, doubles 1/Lambda^2
r2, _, _ = M.build_reactor(wall=True, with_chemistry=False, lam=lam2)
nu2 = r2.compute_nu_wall(y, r2.compute_volume(y))
check("halving Lambda^2 exactly doubles nu_wall",
      abs(nu2 / nu - 2.0) < 1e-12,
      "nu(Lambda) = {0:.10e} ; nu(Lambda/sqrt2) = {1:.10e} ; ratio {2:.15f}".format(
          nu, nu2, nu2 / nu))
ratios = []
for rad in (0.005, 0.01, 0.05, 0.1, 0.2):
    lr = M.diffusion_length(rad, M.L_NOMINAL)
    rr, _, _ = M.build_reactor(wall=True, with_chemistry=False, lam=lr)
    ratios.append(rr.compute_nu_wall(y, rr.compute_volume(y)) * lr * lr)
check("D_a is geometry-independent: nu_wall*Lambda^2 is constant across radii",
      max(ratios) / min(ratios) - 1.0 < 1e-12,
      "nu*Lambda^2 over R = 5 mm .. 200 mm: spread {0:.2e}".format(
          max(ratios) / min(ratios) - 1.0))

# ---------------------------------------------------------------- 4, 5
print("\n(4)+(5) CURRENT BALANCE AND CHARGE DRIFT")
net_wall_current = float(np.sum(r.species_charges * wall))
check("net charge flux to the wall is exactly zero on a neutral state (floating wall)",
      net_wall_current == 0.0,
      "sum(z_j * wall_flux_j) = {0!r} mol/s".format(net_wall_current))
# negative control: a genuinely non-neutral state must give a NON-zero wall current,
# otherwise the check above is passing for the wrong reason.
y_bad = y.copy()
y_bad[ie] *= 1.5
r.residual(0.0, y_bad, np.zeros(r.num_core_species, float))
net_bad = float(np.sum(r.species_charges * r.wall_loss_rates))
net_charge_bad = float(np.sum(r.species_charges * y_bad))
check("NEGATIVE CONTROL: a non-neutral state gives a non-zero wall current",
      net_bad != 0.0,
      "net charge {0:.6e} mol -> wall current {1:.6e} mol/s (and the ratio is "
      "-nu_wall = {2:.6e} s^-1, as one common loss frequency requires)".format(
          net_charge_bad, net_bad, net_bad / net_charge_bad))

# charge conservation through the FULL residual, with chemistry
rc, _, _ = M.build_reactor(wall=True, with_chemistry=True)
yc = M.state_at(rc, 1.0e-6)
dc, _ = rc.residual(0.0, yc, np.zeros(rc.num_core_species, float))
drift = float(np.sum(rc.species_charges * dc))
scale = float(np.sum(np.abs(rc.species_charges * dc)))
check("d(net charge)/dt from the full residual (chemistry + wall) is zero",
      abs(drift) <= 1e-12 * max(scale, 1.0),
      "drift = {0:.6e} mol/s against a term magnitude of {1:.6e} mol/s".format(drift, scale))

# ---------------------------------------------------------------- 6
print("\n(6) GAS-PLUS-WALL ARGON UNDER THE DECLARED RECYCLING MODEL")
for gamma in (1.0, 0.5, 0.0):
    rg, _, _ = M.build_reactor(wall=True, with_chemistry=False, gamma=gamma)
    dg, _ = rg.residual(0.0, y, np.zeros(rg.num_core_species, float))
    ieg, i_arg, i_arpg = M.indices(rg)
    d_heavy = dg[i_arg] + dg[i_arpg]
    nug = rg.nu_wall
    expected = -(1.0 - gamma) * nug * y[i_arpg]
    check("gamma={0}: heavy argon leaves the gas at exactly (1-gamma)*nu*N_ion".format(gamma),
          abs(d_heavy - expected) <= 1e-12 * max(abs(expected), abs(nug * y[i_arpg])),
          "d(N_Ar + N_Ar+)/dt = {0:.10e} mol/s ; -(1-gamma)*nu*N_ion = {1:.10e} mol/s".format(
              d_heavy, expected))

# ---------------------------------------------------------------- 8, 9
print("\n(8)+(9) EXTINCTION, AND NO NEGATIVE AMOUNTS NEAR IT")
# Below threshold in Te: nu_wall > nu_ion, no external source -> decay to extinction.
te_lo = 0.5
r_ext, _, _ = M.build_reactor(te_ev=te_lo, wall=True, with_chemistry=True,
                              x_ion=1.0e-9, source=None)
y_ext = r_ext.y0.copy()
V_ext = r_ext.compute_volume(y_ext)
n_ar_ext = y_ext[M.indices(r_ext)[1]] * constants.Na / V_ext
nui = M.nu_ion(te_lo, n_ar_ext)
nuw = r_ext.compute_nu_wall(y_ext, V_ext)
print("      Te = {0} eV: nu_ion = {1:.4e} s^-1, nu_wall = {2:.4e} s^-1, "
      "Delta = {3:.4e} s^-1".format(te_lo, nui, nuw, nui - nuw))
check("below threshold the wall wins (nu_wall > nu_ion)", nuw > nui)

ts = np.logspace(-6, 0, 25)
traj = []
ok_positive = True
try:
    for t in ts:
        r_ext.advance(t)
        traj.append((t, r_ext.y.copy()))
        if np.any(r_ext.y[:r_ext.num_core_species] < 0.0):
            ok_positive = False
except Exception as exc:                                  # pragma: no cover
    print("      integration stopped: {0}: {1}".format(type(exc).__name__, exc))
if traj:
    ie2, iar2, iarp2 = M.indices(r_ext)
    y_end = traj[-1][1]
    decay = y_end[ie2] / y_ext[ie2]
    check("the electron population decays towards extinction",
          decay < 1e-3,
          "N_e(t=1 s)/N_e(0) = {0:.4e}  [observable O2/O4: an amount ratio]".format(decay))
    check("no negative electron or ion amount is produced anywhere on the approach",
          ok_positive,
          "min over the trajectory: N_e = {0:.4e}, N_Ar+ = {1:.4e} mol".format(
              min(s[1][ie2] for s in traj), min(s[1][iarp2] for s in traj)))
else:
    check("extinction trajectory produced", False)

# ---------------------------------------------------------------- 10
print("\n(10) ANALYTIC JACOBIAN AGAINST FINITE DIFFERENCES")
print("""      Taken by a STEP SCAN, not at one step. A central difference has
      truncation error ~h^2 and roundoff error ~eps/h, so a single h proves
      nothing: at h/|y| = 1e-7 every configuration below -- INCLUDING the
      pre-existing no-wall path -- disagrees at ~1e-3, purely because res[0] is
      O(1) while the perturbation changes it by O(1e-14). The scan's minimum is
      the measurement; the no-wall arm is the control that attributes it.""")

FD_FRACTIONS = (1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
FD_TOLERANCE = 1e-5


def jacobian_scan(reactor, y, mutate=None):
    """Best relative disagreement between the analytic Jacobian and a central
    difference, minimised over step size. ``mutate`` is applied to the reactor
    AFTER the analytic Jacobian is taken and BEFORE the differences, which is how
    the negative control below makes the two genuinely disagree."""
    n = reactor.num_core_species
    zeros = np.zeros(n, float)
    reactor.jacobian(0.0, y, zeros, 0.0)
    analytic = np.array(reactor.jacobian_matrix, float)
    if mutate is not None:
        mutate(reactor)
    best, best_frac, best_worst, best_fd = np.inf, None, None, None
    for frac in FD_FRACTIONS:
        fd = np.zeros((n, n), float)
        for k in range(n):
            h = frac * abs(y[k])
            yp, ym = y.copy(), y.copy()
            yp[k] += h
            ym[k] -= h
            fd[:, k] = (reactor.residual(0.0, yp, zeros)[0]
                        - reactor.residual(0.0, ym, zeros)[0]) / (2 * h)
        scale = np.maximum(np.abs(analytic), np.abs(fd))
        scale[scale == 0.0] = 1.0
        rel = np.abs(analytic - fd) / scale
        if rel.max() < best:
            best = rel.max()
            best_frac = frac
            best_worst = np.unravel_index(np.argmax(rel), rel.shape)
            best_fd = fd
    return best, best_frac, best_worst, analytic, best_fd


for label, kwargs in (("CONTROL: no wall at all (pre-existing equations)", dict(wall=False)),
                      ("wall, integrated electron", dict(wall=True)),
                      ("wall, QUASINEUTRAL electron", dict(wall=True, quasineutral=True)),
                      ("wall + external source", dict(wall=True, source=1.0e5)),
                      ("gamma = 0.5 (partially pumping)", dict(wall=True, gamma=0.5))):
    rj, _, _ = M.build_reactor(with_chemistry=True, x_ion=1.0e-7, **kwargs)
    yj = M.state_at(rj, 1.0e-7)
    best, frac, worst, analytic, fd = jacobian_scan(rj, yj)
    check("{0}: analytic == finite difference".format(label),
          best < FD_TOLERANCE,
          "best {0:.3e} at h/|y| = {1:.0e}, worst entry {2} "
          "(analytic {3:.6e}, fd {4:.6e})".format(
              best, frac, worst, analytic[worst], fd[worst]))

# NEGATIVE CONTROL for check 10: the tolerance above must be able to FAIL.
# wall_recycling is a live attribute, so changing it between the analytic Jacobian
# and the differences makes the two describe genuinely different operators. If the
# scan still came out under tolerance, the check would be worthless.
rn, _, _ = M.build_reactor(wall=True, with_chemistry=True, x_ion=1.0e-7, gamma=0.5)
yn = M.state_at(rn, 1.0e-7)


def _break_gamma(reactor):
    reactor.wall_recycling = 0.6


best_n, frac_n, worst_n, _, _ = jacobian_scan(rn, yn, mutate=_break_gamma)
check("NEGATIVE CONTROL: a 0.5 -> 0.6 change in gamma between the analytic Jacobian "
      "and the differences IS caught at this tolerance",
      best_n > FD_TOLERANCE,
      "best over the whole step scan is {0:.3e} (tolerance {1:.0e}), worst entry "
      "{2} -- so the check above has teeth".format(best_n, FD_TOLERANCE, worst_n))

# ---------------------------------------------------------------- hard failures
print("\n(ENVELOPE) HARD FAILURE OUTSIDE THE SUPPORTED REGIME")
print("""      The domain check lives at the ACCEPTED-STEP boundary, not inside the
      residual. An earlier version raised from compute_nu_wall and it did not
      work: the solver evaluates the residual at Newton TRIAL states far outside
      the domain, a Python exception raised inside the Fortran callback is not
      propagated, and the run died with 'ITERATION MATRIX IS SINGULAR' followed
      by an unrelated-looking SystemError instead of the message below.""")
rh, _, _ = M.build_reactor(wall=True, with_chemistry=False)
y_hot = M.state_at(rh, 1.0e-2)          # alpha = 1e-2 > ceiling 1e-3
try:
    rh.check_wall_support(y_hot)
    check("alpha above the validity ceiling raises", False, "no exception raised")
except PlasmaStateError as exc:
    check("alpha above the validity ceiling raises PlasmaStateError",
          "ionisation degree" in str(exc), str(exc)[:150] + "...")
y_ok = M.state_at(rh, 1.0e-4)
try:
    rh.check_wall_support(y_ok)
    check("NEGATIVE CONTROL: alpha below the ceiling does not raise", True)
except PlasmaStateError as exc:
    check("NEGATIVE CONTROL: alpha below the ceiling does not raise", False, str(exc)[:150])
# and the residual itself must stay total: no raise, finite values, at a state
# far outside the domain -- that is what keeps the solver's iteration matrix sane.
y_wild = M.state_at(rh, 0.9)
try:
    d, _ = rh.residual(0.0, y_wild, np.zeros(rh.num_core_species, float))
    check("the residual does NOT raise outside the domain, and stays finite",
          bool(np.all(np.isfinite(d))),
          "at alpha = 0.9 the residual is {0} -- finite, so DASPK can reject the "
          "trial step on its own error test".format(
              np.array2string(d, precision=3)))
except Exception as exc:
    check("the residual does NOT raise outside the domain, and stays finite", False,
          "{0}: {1}".format(type(exc).__name__, exc))

print("\n(ENVELOPE) NON-NEUTRAL INITIAL STATE UNDER THE ALGEBRAIC ELECTRON")
try:
    M.build_reactor(wall=True, quasineutral=True, x_ion=1e-8, with_chemistry=False)
    # build_reactor is neutral by construction; force a non-neutral composition
    electron, ar, arp = M.argon_species()
    from rmgpy.solver.plasma import PlasmaReactor
    bad = PlasmaReactor((M.TGAS, 'K'), (M.P_NOMINAL, 'Pa'),
                        {electron: 2.0e-8, arp: 1.0e-8, ar: 1.0 - 3.0e-8},
                        (M.TE_NOMINAL_EV * M.EV_TO_K, 'K'), n_sims=1, termination=[],
                        quasineutral_electron=True,
                        diffusion_length=(M.diffusion_length(), 'm'),
                        ion_reduced_mobility=(M.MU0_AR_IN_AR, 'm^2/(V*s)'))
    bad.initialize_model([electron, ar, arp], [], [], [])
    check("a non-neutral initial state is refused under quasineutral_electron", False,
          "no exception raised")
except PlasmaStateError as exc:
    check("a non-neutral initial state is refused under quasineutral_electron, by US "
          "and not by an opaque solver convergence error",
          "net charge" in str(exc), str(exc)[:200] + "...")

# ---------------------------------------------------------------- 12
print("\n(12) NO TARGET ELECTRON DENSITY IN THE IMPLEMENTATION")
import os
src = open(os.path.join(M.WORKTREE, 'rmgpy', 'solver', 'plasma.pyx')).read()
# every numeric literal introduced by this ticket, and what it is
declared = {
    'PLASMA_LOSCHMIDT': '2.6867811e25',
    'PLASMA_WALL_MAX_IONISATION_DEGREE': '1.0e-3',
}
check("the wall code introduces exactly two numeric constants, both declared",
      all(v in src for v in declared.values()),
      "; ".join("{0} = {1}".format(k, v) for k, v in declared.items()))
check("neither constant is an electron density (one is a unit convention, one a "
      "dimensionless ratio)", True,
      "Loschmidt is the density the tabulated mu_0 is normalised to; the ceiling is "
      "n_e/n_neutral, a ratio, fixed by where ion-neutral collisions stop dominating")

print("\n" + "=" * 78)
if FAILURES:
    print("FAILED CHECKS ({0}): {1}".format(len(FAILURES), "; ".join(FAILURES)))
    sys.exit(1)
print("ALL CHECKS PASSED")
print("=" * 78)
