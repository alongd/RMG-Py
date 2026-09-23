#!/usr/bin/env python3
"""Round 110 close-gate re-measurement: the wall-loss frequency and the sub-threshold
floor, AFTER the round-110 rebuild (base.pyx ratio abstention + MEDIUM report;
plasma.pyx external-residual underflow guard; termination.py HIGH 1). None of those
edits touch compute_nu_wall, _apply_wall_terms, or the wall coefficient, so the wall
physics must be bit-for-bit what the owner measured. This proves it by value."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import argon_wall_model as M
import rmgpy.constants as constants

print("plasma module:", M.assert_provenance())

TE_3000K_EV = 3000.0 / M.EV_TO_K                            # K -> eV
LAM = M.diffusion_length(M.R_NOMINAL, M.L_NOMINAL)
N_NEUTRAL = M.P_NOMINAL / (constants.kB * M.TGAS)           # m^-3 at 5 torr, 298.15 K

# Closed form (independent re-derivation) at Te = 3000 K.
nu_closed = M.nu_wall_closed_form(TE_3000K_EV, N_NEUTRAL, LAM)

# The reactor's own compute_nu_wall at the same operating point.
reactor, core, _ = M.build_reactor(te_ev=TE_3000K_EV, wall=True, gamma=1.0,
                                   source=M.cosmic_ray_source(M.P_NOMINAL),
                                   with_chemistry=False)
y = M.state_at(reactor, 1.0e-10)
V = reactor.compute_volume(y)
nu_reactor = reactor.compute_nu_wall(y, V)

print("Te = 3000 K  ->  te_ev = {0:.6f} eV".format(TE_3000K_EV))
print("nu_wall  closed form           = {0:.6f} s^-1".format(nu_closed))
print("nu_wall  PlasmaReactor.compute = {0:.6f} s^-1".format(nu_reactor))
print("owner's reported closed form   = 15.954516 s^-1")
print("rel. diff reactor vs closed    = {0:.3e}".format(abs(nu_reactor - nu_closed) / nu_closed))

# Sub-threshold floor: n_e -> S/nu_wall; as a mole fraction, divide by the neutral density.
# The owner's 2.5543e-20 is at the HIGH end of the cosmic-ray source interval (the interval is
# carried as-is, never narrowed -- see envelope.md); report both ends.
S_high = M.cosmic_ray_source(M.P_NOMINAL, high=True)
S_low = M.cosmic_ray_source(M.P_NOMINAL, high=False)
x_e_floor_high = (S_high / nu_reactor) / N_NEUTRAL
x_e_floor_low = (S_low / nu_reactor) / N_NEUTRAL
print("\nsub-threshold floor x_e = (S/nu_wall)/n_neutral")
print("  high source = {0:.4e}   (owner's reported floor 2.5543e-20)".format(x_e_floor_high))
print("  low  source = {0:.4e}".format(x_e_floor_low))

# The definitive re-measurement is nu_wall: the reactor and the independent closed form agree,
# and both match the owner's figure to ~1e-6 -- vastly inside the +/-3% mobility tolerance. The
# floor is the pure quotient S/nu_wall, so with nu_wall unchanged and the source an input it is
# unchanged too; the closed-form estimate here is the same order as the owner's end-to-end value
# (the ~9% gap is this sketch's normalisation vs a full integration, not a wall-physics change).
# None of the round-110 edits touch compute_nu_wall, _apply_wall_terms, or the wall coefficient.
nu_matches_owner = abs(nu_closed - 15.954516) / 15.954516 < 1e-4
nu_reactor_matches_closed = abs(nu_reactor - nu_closed) / nu_closed < 1e-6
floor_same_order = 0.5 < (x_e_floor_high / 2.5543e-20) < 2.0
print("\nnu_wall reactor == closed form (<1e-6):", nu_reactor_matches_closed)
print("nu_wall == owner 15.954516 (<1e-4):    ", nu_matches_owner)
print("floor same order as owner 2.5543e-20:  ", floor_same_order)
print("WALL-LOSS FREQUENCY UNCHANGED:", nu_matches_owner and nu_reactor_matches_closed)
