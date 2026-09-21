#!/usr/bin/env python3
"""
Verifier check 7: with no wall, the reactor reproduces the PRE-CHANGE volume-only
equations BIT-FOR-BIT -- not approximately, and not against a profile this branch
generated for itself.

The only honest way to show that is to run the pre-change code. This script builds
a second extension module from `git show <base>:rmgpy/solver/plasma.pyx`, imports
both, drives them through identical states, and compares the raw doubles with
`==`. If one bit of the residual or of the Jacobian moved on a wall-less reactor,
this fails.

Run it through `build_and_verify_zero_wall.sh`, which does the build and the
cleanup around it.
"""

import sys

import numpy as np

import argon_wall_model as M

print("=" * 78)
print("CHECK 7: ZERO-WALL BIT-FOR-BIT RECOVERY OF THE PRE-CHANGE EQUATIONS")
print("plasma module (this branch): {0}".format(M.assert_provenance()))

from rmgpy.solver.plasma import PlasmaReactor as NewReactor
from rmgpy.solver.plasma_base_i246 import PlasmaReactor as BaseReactor
import rmgpy.solver.plasma_base_i246 as base_mod

print("plasma module (base 311818121): {0}".format(base_mod.__file__))
print("=" * 78)

failures = []


def build(cls, **kwargs):
    electron, ar, arp = M.argon_species()
    x = 1.0e-7
    imf = {electron: x, arp: x, ar: 1.0 - 2.0 * x}
    r = cls((M.TGAS, 'K'), (M.P_NOMINAL, 'Pa'), imf,
            (M.TE_NOMINAL_EV * M.EV_TO_K, 'K'), n_sims=1, termination=[], **kwargs)
    r.initialize_model([electron, ar, arp],
                       [M.ionisation_reaction(electron, ar, arp)], [], [])
    return r


new = build(NewReactor)            # no wall parameters at all -> has_wall is False
base = build(BaseReactor)
print("\nnew reactor has_wall = {0}, quasineutral_electron = {1}".format(
    new.has_wall, new.quasineutral_electron))

n = new.num_core_species
assert n == base.num_core_species

print("\n--- identical initial packing -----------------------------------------")
same_y0 = np.array_equal(new.y0, base.y0)
print("y0 identical                 : {0}".format(same_y0))
print("  new  y0 = {0}".format(list(new.y0)))
print("  base y0 = {0}".format(list(base.y0)))
if not same_y0:
    failures.append("y0")
same_dydt0 = np.array_equal(new.dydt0, base.dydt0)
print("dydt0 identical              : {0}".format(same_dydt0))
if not same_dydt0:
    failures.append("dydt0")
same_V = new.compute_volume(new.y0) == base.compute_volume(base.y0)
print("compute_volume identical     : {0}  ({1!r})".format(same_V, new.compute_volume(new.y0)))
if not same_V:
    failures.append("compute_volume")

print("\n--- residual and Jacobian over a range of states ----------------------")
print("{0:>12} {1:>14} {2:>14} {3:>10} {4:>10}".format(
    "alpha", "max|d residual|", "max|d jacobian|", "res bitwise", "jac bitwise"))
print("-" * 78)
rng = np.random.default_rng(20260921)
all_res_exact = True
all_jac_exact = True
for alpha in (0.0, 1e-12, 1e-9, 1e-7, 1e-5, 1e-4, 5e-4):
    y = M.state_at(new, alpha)
    dydt = rng.normal(size=n) * 1e-3     # a non-zero dydt exercises the -dydt term too
    for cj in (0.0, 1.0, 1.0e6):
        dn, in_ = new.residual(0.0, y.copy(), dydt.copy())
        db, ib = base.residual(0.0, y.copy(), dydt.copy())
        pn = np.array(new.jacobian(0.0, y.copy(), dydt.copy(), cj), float)
        pb = np.array(base.jacobian(0.0, y.copy(), dydt.copy(), cj), float)
        res_exact = np.array_equal(dn, db) and in_ == ib
        jac_exact = np.array_equal(pn, pb)
        all_res_exact = all_res_exact and res_exact
        all_jac_exact = all_jac_exact and jac_exact
        if cj == 1.0:
            print("{0:>12.1e} {1:>14.3e} {2:>14.3e} {3:>10} {4:>10}".format(
                alpha, float(np.max(np.abs(dn - db))), float(np.max(np.abs(pn - pb))),
                str(res_exact), str(jac_exact)))
        if not res_exact:
            failures.append("residual at alpha={0:g}, cj={1:g}".format(alpha, cj))
        if not jac_exact:
            failures.append("jacobian at alpha={0:g}, cj={1:g}".format(alpha, cj))

# jacobian_matrix (the derived d res/d y the rest of RMG reads) must match too
new.jacobian(0.0, M.state_at(new, 1e-7), np.zeros(n), 1.0)
base.jacobian(0.0, M.state_at(new, 1e-7), np.zeros(n), 1.0)
jm_exact = np.array_equal(np.array(new.jacobian_matrix, float),
                          np.array(base.jacobian_matrix, float))
print("\njacobian_matrix identical    : {0}".format(jm_exact))
if not jm_exact:
    failures.append("jacobian_matrix")

print("\n--- NEGATIVE CONTROL --------------------------------------------------")
print("The comparison above must be able to FAIL. Declare a wall on the new")
print("reactor, with the SAME chemistry and the SAME state, and re-run it:")
walled = build(NewReactor,
               diffusion_length=(M.diffusion_length(), 'm'),
               ion_reduced_mobility=(M.MU0_AR_IN_AR, 'm^2/(V*s)'))
y = M.state_at(new, 1e-7)
dw = walled.residual(0.0, y.copy(), np.zeros(n))[0]
db = base.residual(0.0, y.copy(), np.zeros(n))[0]
differs = not np.array_equal(dw, db)
print("  a walled reactor DOES differ from base: {0}   max|delta| = {1:.6e} mol/s".format(
    differs, float(np.max(np.abs(dw - db)))))
if not differs:
    failures.append("NEGATIVE CONTROL: a walled reactor did not differ from base")

print("\n" + "=" * 78)
if failures:
    print("CHECK 7 FAILED: {0}".format("; ".join(failures)))
    sys.exit(1)
print("CHECK 7 PASSED -- residual, Jacobian, jacobian_matrix, y0, dydt0 and the EOS")
print("are bit-for-bit identical to the pre-change build on a wall-less reactor,")
print("across 7 ionisation degrees x 3 cj values, with a non-zero dydt, and the")
print("negative control confirms the comparison can fail.")
print("=" * 78)
