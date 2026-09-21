#!/usr/bin/env python3
"""
I-246: where does a real argon glow discharge sit on this model's branch, and
what is the conditioning there?

READ THIS FIRST. The electron density used below is a REPORTING COMPARISON and
nothing else. It is not a threshold, not a fixture, not an acceptance criterion,
and no parameter anywhere in this ticket was chosen by reference to it. The
ticket's brief states that the model is five to seven orders from a real argon
glow discharge and asks whether that gap closes; answering that requires naming
the number the gap is measured against, once, here, in a script that computes
nothing the implementation depends on.

The question it answers is the ticket's headline one, stated in the only way
that is actually decidable: not "is n_e right?" but "at the place where n_e
would be right, how well determined is it?"
"""

import numpy as np

import rmgpy.constants as constants

import argon_wall_model as M
from sweep import grid_point, wall_constant, CONDITIONING_LIMIT, ALPHA_CEILING

print("plasma module: {0}".format(M.assert_provenance()))
KB = constants.kB

# Typical positive-column electron density for a low-pressure argon glow
# discharge, order of magnitude only: 1e16 - 1e17 m^-3. Comparison point.
N_E_GLOW_LOW, N_E_GLOW_HIGH = 1.0e16, 1.0e17

print("\n" + "=" * 94)
print("WHERE A REAL GLOW DISCHARGE SITS ON THE ABOVE-THRESHOLD BRANCH")
print("=" * 94)
print("Comparison point (NOT a criterion): n_e = {0:.0e} - {1:.0e} m^-3.".format(
    N_E_GLOW_LOW, N_E_GLOW_HIGH))
print()

for p_torr, radius in ((5.0, 0.05), (1.0, 0.05), (5.0, 0.01)):
    pa = p_torr * M.TORR_TO_PA
    lam = M.diffusion_length(radius, M.L_NOMINAL)
    n_ar0 = pa / (KB * M.TGAS)

    # threshold
    lo, hi = 0.3, 6.0
    for _ in range(90):
        mid = 0.5 * (lo + hi)
        if grid_point(mid, pa, radius)['dnu'] < 0:
            lo = mid
        else:
            hi = mid
    te_thr = hi
    te_k = te_thr * M.EV_TO_K

    print("-" * 94)
    print("p = {0} torr, R = {1} m, Lambda = {2:.4f} cm, n_Ar,0 = {3:.4e} m^-3, "
          "Te_thr = {4:.6f} eV".format(p_torr, radius, lam * 100, n_ar0, te_thr))
    print()
    print("{0:>14} {1:>12} {2:>12} {3:>13} {4:>13} {5:>11}".format(
        "n_e (m^-3)", "f*", "alpha*", "Te needed(eV)", "Te-Te_thr(eV)", "C_wall"))
    print("-" * 94)

    for n_e in (N_E_GLOW_LOW, N_E_GLOW_HIGH, 1e18, 1e20, ALPHA_CEILING * n_ar0):
        # invert n_e* = n_Ar,0*Tgas*f/(Tgas+Te) for f, then n_Ar* = (1-f)*n_Ar,0,
        # then k_iz = K/n_Ar*^2 and finally the Te that gives that k_iz.
        f = n_e * (M.TGAS + te_k) / (n_ar0 * M.TGAS)
        if not (0.0 < f < 1.0):
            continue
        n_ar_star = (1.0 - f) * n_ar0
        alpha = n_e / n_ar_star
        # solve for Te: nu_ion(Te) = nu_wall(Te) at THIS neutral density
        a, b = 0.3, 6.0
        for _ in range(90):
            mid = 0.5 * (a + b)
            K = wall_constant(mid, lam)
            if M.k_ionisation(mid) * n_ar_star < K / n_ar_star:
                a = mid
            else:
                b = mid
        te_needed = b
        c_wall = (1.0 - f) / (2.0 * f)
        print("{0:>14.2e} {1:>12.4e} {2:>12.4e} {3:>13.6f} {4:>13.3e} {5:>11.4e}".format(
            n_e, f, alpha, te_needed, te_needed - te_thr, c_wall))

print("-" * 94)
print("""
READING

A real glow-discharge electron density IS on this model's branch. It is inside
the transport model's validity ceiling, so the model is not "five to seven orders
away" in the sense of being unable to reach the number -- it reaches it.

What it cannot do is DETERMINE it. At n_e = 1e16 m^-3 the conditioning number
C_wall = |d ln n_e / d ln nu_wall| is of order 1e5, thousands of times past the
limit predeclared in envelope.md. Concretely: the ion mobility this model is
built on is known to +/-3 per cent (Ellis et al. 1976), and a 3 per cent error in
nu_wall there moves n_e by a factor of exp(0.03 * C_wall) -- which is not a
correction, it is a different universe.

Equivalently, and this is the more useful statement: the electron temperature
needed to sit at a glow-discharge density differs from the threshold temperature
by a few parts in ten thousand of an eV. Te is PRESCRIBED in this model, not
solved for. Nobody knows it that well, and no wall operator can supply it.

So the gap does not close by choosing a wall coefficient -- it closes, if at all,
by closing the discharge power balance, which is what would determine Te instead
of prescribing it. That is the next milestone and it is out of scope here, by
design. This is the result.
""")
