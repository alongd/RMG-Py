#!/usr/bin/env python3
"""
I-246: what does the algebraic (quasineutral) electron do to the conditioning,
and what happened to the 784x?

The governing ruling measured that near the ionisation/wall-loss crossing a 0.01%
perturbation in the electron temperature moved the computed electron density by
~784x. Two separate things are measured here, because they are separate:

  A  the CONDITION NUMBER of the Jacobian, integrated electron vs algebraic one,
     on matched states. This is what "removes one stiff direction from the
     Jacobian" has to mean if it means anything.

  B  the SENSITIVITY d ln n_e / d ln Te, taken by bracketing ALONG each branch
     and never across the extinction boundary, as envelope.md predeclared.

A and B answer different questions. A is about the linear algebra the solver
does; B is about whether the answer is determined. Improving A does not improve
B, and saying otherwise would be the substitution taking credit for something it
does not do.
"""

import numpy as np

import rmgpy.constants as constants

import argon_wall_model as M
from sweep import grid_point

print("plasma module: {0}".format(M.assert_provenance()))

# ====================================================================== A
print("\n" + "=" * 92)
print("A  JACOBIAN CONDITION NUMBER: INTEGRATED vs ALGEBRAIC ELECTRON")
print("=" * 92)
print("Same chemistry, same wall, same state, same cj. cond_2 of d(delta)/dy.")
print()
print("""The RAW condition number is not the quantity to look at and is reported only
so that saying so is checkable: the state vector spans twelve decades (y_Ar ~ 1
against y_e ~ 1e-12), so cond_2 of the bare matrix is dominated by that spread
and comes out at 1e23 or infinite for BOTH configurations. What the solver
actually solves is the equilibrated system, columns scaled by the error weight
w_j = atol + rtol*|y_j| and rows by the same. That is the number below.""")
print()
print("""Three measures, because the choice of row scaling for an ALGEBRAIC row is a real
judgement call and the answer must not depend on my making it:
  raw   -- no scaling. This is literally the matrix DASPK factors (LU with
           partial pivoting, no equilibration), so it is not merely academic.
  ewt   -- rows and columns scaled by the error weight w = atol + rtol*|y|.
           Arguable for the algebraic row, whose residual has units of y rather
           than of dy/dt.
  equil -- full equilibration: columns by w, then each row by its own largest
           entry. Chooses nothing; it is the scaling-neutral comparison.""")
print()
ATOL, RTOL = 1e-16, 1e-8


def _cond_equilibrated(J, w):
    A = J * w[None, :]
    rownorm = np.max(np.abs(A), axis=1)
    rownorm[rownorm == 0.0] = 1.0
    return np.linalg.cond(A / rownorm[:, None])


print("{0:>9} {1:>11} {2:>11} {3:>11} {4:>11} {5:>11} {6:>11}".format(
    "alpha", "raw integ", "raw quasi", "ewt integ", "ewt quasi", "eq integ", "eq quasi"))
print("-" * 92)

ratios = []
equil_ratios = []
for alpha in (1e-12, 1e-10, 1e-8, 1e-6, 1e-5, 1e-4, 5e-4):
    raw, scaled, equil = {}, {}, {}
    for qn in (False, True):
        r, _, _ = M.build_reactor(wall=True, with_chemistry=True, x_ion=1e-10,
                                  quasineutral=qn)
        y = M.state_at(r, alpha)
        n = r.num_core_species
        # cj as DASPK forms it for a step h: cj ~ 1/h. Use a step comparable to
        # the fastest time scale in the problem, 1/nu_wall.
        cj = r.compute_nu_wall(y, r.compute_volume(y)) or 1.0
        r.jacobian(0.0, y, np.zeros(n), cj)
        J = np.array(r.jacobian_matrix, float) - cj * np.identity(n)
        if r.quasineutral_electron:          # the algebraic row carries no -cj
            J[r.electron_index, :] = np.array(r.jacobian_matrix, float)[r.electron_index, :]
        raw[qn] = np.linalg.cond(J)
        w = ATOL + RTOL * np.abs(y)
        scaled[qn] = np.linalg.cond(J * (w[None, :] / w[:, None]))
        equil[qn] = _cond_equilibrated(J, w)
    ratios.append(scaled[False] / scaled[True])
    equil_ratios.append(equil[False] / equil[True])
    print("{0:>9.1e} {1:>11.3e} {2:>11.3e} {3:>11.3e} {4:>11.3e} {5:>11.3e} {6:>11.3e}".format(
        alpha, raw[False], raw[True], scaled[False], scaled[True],
        equil[False], equil[True]))

print()
print("ratio (integrated / quasineutral); BELOW 1 means the algebraic electron is")
print("WORSE conditioned, above 1 means better:")
for name, rs in (("error-weight scaled", ratios), ("equilibrated", equil_ratios)):
    f = [x for x in rs if np.isfinite(x)]
    print("  {0:<22} min {1:.4g}   median {2:.4g}   max {3:.4g}".format(
        name, min(f), float(np.median(f)), max(f)))
print()
print("""The three scalings DISAGREE, by a lot, and the disagreement is the finding.

Equilibrated, the algebraic system is about 1e5 times better conditioned
(cond ~ 2.6 against ~2.5e5): the charge row is exactly linear and perfectly
conditioned in itself, so intrinsically the substitution does what it was
supposed to do. Raw and error-weight-scaled, it is about the same as the
integrated system, and at first it was 175x WORSE.

That 175x was a pure row-scaling artefact, and it is worth stating because it
was a live defect rather than a curiosity. The charge row's entries are the
charges themselves, exactly +1 and -1, while the differential rows carry entries
of order of the rates, 1e2-1e5. DASPK factors the iteration matrix with partial
pivoting and NO equilibration, so leaving an order-1 row among them costs real
accuracy. Scaling that row by a constant -- an exact operation on a DAE, which
changes the factored matrix and nothing about the solution -- removes all of it.
PlasmaReactor now does that automatically (charge_row_scale, set once at
initialization from the differential rows), which is why the raw column above
reads ~6.6e5 rather than ~9.4e7.

So the honest summary of A: with its row scaled, the algebraic electron is
conditioning-NEUTRAL in the matrix the solver actually factors, and a large
improvement in the matrix an equilibrating solver would factor. It is also,
independently, CORRECT: the two configurations agree on the ionisation degree at
the stopping point to 13 significant figures (Q3 of confirm_branches). And the
electron stops being a differential unknown, so its error test and its step-size
vote are gone. None of that is the 784x.""")

# ====================================================================== B
print("\n" + "=" * 92)
print("B  d ln n_e / d ln Te, BY BRACKETING ALONG EACH BRANCH")
print("=" * 92)
print("A 0.01% perturbation in Te, exactly as the governing ruling applied it,")
print("with the resulting change in n_e stated as a FACTOR so it is comparable")
print("to the 784x directly.")
print()
h = 1.0e-4                         # the ruling's 0.01%

print("--- B1  the unbounded growth branch: no wall, fixed integration time ---")
print("This is the regime the 784x was measured in: nothing stops the growth, so")
print("n_e at a fixed time is an exponential of a rate that depends on Te.")
print("{0:>10} {1:>16} {2:>16} {3:>14} {4:>12}".format(
    "Te (eV)", "N_e(t) at Te", "N_e(t) at Te*1.0001", "factor", "dln/dln"))
print("-" * 92)
for te in (1.2, 1.5):
    vals = []
    for scale in (1.0, 1.0 + h):
        r, _, _ = M.build_reactor(te_ev=te * scale, wall=False, with_chemistry=True,
                                  x_ion=1.0e-12)
        ie = r.electron_index
        r.advance(3.0e-5)
        vals.append(r.y[ie])
    factor = vals[1] / vals[0]
    print("{0:>10.2f} {1:>16.6e} {2:>16.6e} {3:>14.4e} {4:>12.4e}".format(
        te, vals[0], vals[1], factor, np.log(factor) / h))

print("\n--- B2  the sub-threshold source-sustained branch, WITH the wall ---")
print("n_e = S_ext / (nu_wall - nu_ion): a bounded, unique steady state. The")
print("sensitivity is taken by bracketing along the branch, never across the")
print("boundary at the top of it.")
print("{0:>10} {1:>15} {2:>15} {3:>14} {4:>12}  {5}".format(
    "Te (eV)", "n_e (m^-3)", "n_e at Te*1.0001", "factor", "dln/dln", "distance to boundary"))
print("-" * 110)
te_thr = None
lo, hi = 0.3, 6.0
for _ in range(80):
    mid = 0.5 * (lo + hi)
    if grid_point(mid, M.P_NOMINAL, M.R_NOMINAL)['dnu'] < 0:
        lo = mid
    else:
        hi = mid
te_thr = hi
for te in (0.5, 0.7, 0.8, 0.84, 0.85, 0.8520):
    if te >= te_thr:
        continue
    g0 = grid_point(te, M.P_NOMINAL, M.R_NOMINAL)
    g1 = grid_point(te * (1 + h), M.P_NOMINAL, M.R_NOMINAL)
    if not (np.isfinite(g0['ne_sub_hi']) and np.isfinite(g1['ne_sub_hi'])):
        continue
    factor = g1['ne_sub_hi'] / g0['ne_sub_hi']
    print("{0:>10.4f} {1:>15.6e} {2:>15.6e} {3:>14.6f} {4:>12.4e}  {5:.3e} eV".format(
        te, g0['ne_sub_hi'], g1['ne_sub_hi'], factor, np.log(factor) / h, te_thr - te))

print("\nthreshold Te = {0:.6f} eV at 5 torr, R = 0.05 m".format(te_thr))

print("\n--- B3  where the 784x actually lives ---")
print("Sweeping the distance to the boundary, on the branch, to find where a")
print("0.01% Te perturbation moves n_e by a factor of 784:")
print("{0:>14} {1:>16} {2:>16}".format("Te_thr - Te (eV)", "factor on n_e", "C_Te"))
print("-" * 92)
for gap in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
    te = te_thr - gap
    if te <= 0:
        continue
    g0 = grid_point(te, M.P_NOMINAL, M.R_NOMINAL)
    g1 = grid_point(te * (1 + h), M.P_NOMINAL, M.R_NOMINAL)
    if not (np.isfinite(g0['ne_sub_hi']) and np.isfinite(g1['ne_sub_hi'])):
        print("{0:>14.1e} {1:>16} {2:>16}".format(gap, "branch ends", "-"))
        continue
    factor = g1['ne_sub_hi'] / g0['ne_sub_hi']
    print("{0:>14.1e} {1:>16.6e} {2:>16.4e}".format(gap, factor, np.log(factor) / h))

print("\n" + "=" * 92)
print("READING")
print("=" * 92)
print("""
A 0.01% move in Te produces a factor on n_e that is ~1 far from the boundary and
diverges as the boundary is approached. The 784x is therefore not a property of
the model, of the solver, or of the Jacobian: it is a reading taken at a
particular distance from the sustainment boundary, and it can be made any number
at all by moving along the branch. What IS a property of the model is that the
sensitivity diverges at the boundary and is O(1) away from it -- so the useful
statement is a distance, not a factor.

Whatever the algebraic electron does to A, it does nothing for B, and that is
correct rather than disappointing: B is set by the physics of the balance, not
by how the state vector is packed. No repacking of the unknowns can make a
quantity determined that the equations do not determine.

So: nothing here beats the 784x, and a wall operator was never going to. What
would is closing the discharge power balance -- which would determine Te instead
of prescribing it, and so fix the distance to the boundary, which is the only
thing the sensitivity actually depends on.
""")
