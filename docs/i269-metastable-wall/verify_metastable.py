#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

from rmgpy.solver.base import ReactionSystem, TerminationTime, TerminationConversion, TerminationRateRatio

"""I-269 verifier: neutral-diffusion wall loss of a declared metastable.

Checks, against the I-246 argon deck (5 Torr, Lambda of the 5 cm x 30 cm cylinder):
  1. hand nu_m = (D*p)/p/Lambda^2 at 5 Torr, 300 K equals the operator's value to 1e-6;
  2. the ambipolar nu_wall is unchanged (Te = 3000 K, owner's 15.954516 s^-1), and a
     reactor carrying a declared Ar* returns the bit-identical nu_wall;
  3. a seeded Ar* decays as exp(-nu_m t) through simulate(), and ground Ar gains
     exactly what Ar* (plus the wall-neutralised Ar+) loses;
  4. negative arm: the same run with Ar* undeclared shows no wall loss of Ar*; the only
     Ar* change is the external source's (pre-existing) ionisation draw, predicted here.
Exits non-zero if any check fails.
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'i246-ambipolar-wall-operator'))
import argon_wall_model as M  # noqa: E402
import rmgpy.constants as constants  # noqa: E402
from rmgpy.rmg.settings import ModelSettings, SimulatorSettings  # noqa: E402
from rmgpy.solver.termination import TerminationTime  # noqa: E402
from rmgpy.solver.plasma import PlasmaReactor  # noqa: E402
from rmgpy.species import Species  # noqa: E402

print("plasma module:", M.assert_provenance())

DP = 54.0                                   # cm^2 Torr/s, the brief's figure (source unverified)
CM2_TORR_TO_SI = 1.0e-4 * M.TORR_TO_PA
TE_EV = 3000.0 / M.EV_TO_K                  # the deck Te at which 15.954516 was published
LAM = M.diffusion_length()
failures = []


def check(name, ok, detail):
    print("{0:<58s} {1}   {2}".format(name, "PASS" if ok else "FAIL", detail))
    if not ok:
        failures.append(name)


def build(tgas, declared, x_meta=1.0e-3, x_ion=1.0e-9, source=None, termination=None):
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ground = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    meta = Species(label='Ar*').from_adjacency_list('1 Ar u2 p3 c0')
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    imf = {electron: x_ion, arp: x_ion, meta: x_meta, ground: 1.0 - 2.0 * x_ion - x_meta}
    kwargs = dict(diffusion_length=(LAM, 'm'),
                  ion_reduced_mobility=(M.MU0_AR_IN_AR, 'm^2/(V*s)'),
                  wall_recycling=1.0, wall_neutralization_products={'Ar+': 'Ar'})
    if source is not None:
        kwargs['ionisation_source'] = (source, 'm^-3/s')
    if declared:
        kwargs['wall_neutral_diffusion'] = {
            'Ar*': {'product': 'Ar', 'diffusivity': (DP, 'cm^2*torr/s')}}
    r = PlasmaReactor((tgas, 'K'), (M.P_NOMINAL, 'Pa'), imf, (TE_EV * M.EV_TO_K, 'K'),
                      n_sims=1, termination=termination or [], **kwargs)
    core = [electron, ground, meta, arp]
    r.initialize_model(core, [], [], [])
    return r, core, {s.label: i for i, s in enumerate(core)}


# ---- 1. hand nu_m at 5 Torr, 300 K, deck Lambda ------------------------------------
r, core, idx = build(300.0, True)  # y below carries no charge, so n_neutral = p/kT exactly
y = np.zeros(r.num_core_species)
y[idx['Ar*']], y[idx['Ar']] = 1.0e-3, 1.0 - 1.0e-3
nu_op = r.compute_neutral_wall_frequencies(y, r.compute_volume(y))
hand = DP * CM2_TORR_TO_SI / M.P_NOMINAL / (LAM * LAM)
print("Lambda = {0:.9e} m, p = {1:.6f} Pa, D*p = {2} cm^2 Torr/s".format(LAM, M.P_NOMINAL, DP))
print("nu_m hand     = {0:.9e} s^-1".format(hand))
print("nu_m operator = {0:.9e} s^-1".format(nu_op[idx['Ar*']]))
check("1. nu_m operator == hand (rel < 1e-6)",
      abs(nu_op[idx['Ar*']] / hand - 1.0) < 1e-6,
      "rel = {0:.2e}".format(abs(nu_op[idx['Ar*']] / hand - 1.0)))
others = [nu_op[idx[l]] for l in ('Ar', 'Ar+', 'e-')]
check("1b. nu_m is exactly zero for Ar, Ar+, e-", all(v == 0.0 for v in others), str(others))

# ---- 2. ambipolar nu_wall unchanged ------------------------------------------------
ref, _, _ = M.build_reactor(te_ev=TE_EV, wall=True, gamma=1.0,
                            source=M.cosmic_ray_source(M.P_NOMINAL), with_chemistry=False)
y_ref = M.state_at(ref, 1.0e-10)
nu_wall = ref.compute_nu_wall(y_ref, ref.compute_volume(y_ref))
print("nu_wall (i246 deck, Te=3000 K) = {0:.6f} s^-1 (published 15.954516)".format(nu_wall))
check("2. nu_wall == 15.954516 (rel < 1e-4)", abs(nu_wall / 15.954516 - 1.0) < 1e-4,
      "rel = {0:.2e}".format(abs(nu_wall / 15.954516 - 1.0)))
ra, _, ia = build(M.TGAS, True)
rb, _, ib = build(M.TGAS, False)
ya = np.array(ra.y0[:ra.num_core_species], float)
yb = np.array(rb.y0[:rb.num_core_species], float)
nwa = ra.compute_nu_wall(ya, ra.compute_volume(ya))
nwb = rb.compute_nu_wall(yb, rb.compute_volume(yb))
check("2b. nu_wall declared == undeclared (bitwise)", nwa == nwb, "{0!r} vs {1!r}".format(nwa, nwb))
res_a = np.zeros(ra.num_core_species)
res_b = np.zeros(rb.num_core_species)
ra.residual(0.0, ya, np.zeros_like(ya))
rb.residual(0.0, yb, np.zeros_like(yb))
wa, wb = np.array(ra.wall_loss_rates), np.array(rb.wall_loss_rates)
check("2c. wall rates of e-, Ar+ identical (bitwise)",
      wa[ia['e-']] == wb[ib['e-']] and wa[ia['Ar+']] == wb[ib['Ar+']],
      "e- {0!r}, Ar+ {1!r}".format(wa[ia['e-']], wa[ia['Ar+']]))


# ---- 3/4. seeded decay through simulate(), declared and undeclared -----------------
def run(declared):
    probe, _, _ = build(300.0, declared)
    y0 = np.array(probe.y0[:probe.num_core_species], float)
    v0 = probe.compute_volume(y0)
    n_e = y0[probe.electron_index] * constants.Na / v0
    src = probe.compute_nu_wall(y0, v0) * n_e          # holds the 1e-9 charge seed steady
    t_end = 1.0 / (hand)                                # one e-fold of the declared arm
    r, core, idx = build(300.0, declared, source=src,
                         termination=[TerminationTime((t_end, 's'))])
    y0 = np.array(r.y0[:r.num_core_species], float)
    nu_m = r.compute_neutral_wall_frequencies(y0, r.compute_volume(y0))[idx['Ar*']]
    r.simulate(core, [], [], [], [], [],
               model_settings=ModelSettings(tol_keep_in_edge=0, tol_move_to_core=1e5,
                                            tol_interrupt_simulation=1e8),
               simulator_settings=SimulatorSettings())
    y = np.array(r.y[:r.num_core_species], float)
    # The external source ionises every ionisable neutral in proportion to its amount,
    # Ar* included (pre-existing, not a wall term): its draw on Ar* over the run.
    src_mol_s = src * v0 / constants.Na
    y_ion = y0[idx['Ar']] + y0[idx['Ar*']]
    draw = src_mol_s * y0[idx['Ar*']] / y_ion * r.t
    return y0, y, nu_m, r.t, idx, draw


y0, y, nu_m, t, idx, draw = run(True)
expected = y0[idx['Ar*']] * np.exp(-nu_m * t)
lost = y0[idx['Ar*']] - y[idx['Ar*']]
gained = y[idx['Ar']] - y0[idx['Ar']]
ion_lost = y0[idx['Ar+']] - y[idx['Ar+']]
print("declared:   t = {0:.6e} s, nu_m*t = {1:.6f}, Ar* {2:.9e} -> {3:.9e} (exp {4:.9e})".format(
    t, nu_m * t, y0[idx['Ar*']], y[idx['Ar*']], expected))
print("            Ar* lost {0:.9e}, Ar gained {1:.9e}, Ar+ lost {2:.3e}".format(lost, gained, ion_lost))
check("3. Ar* == y0*exp(-nu_m t) (rel < 1e-4)", abs(y[idx['Ar*']] / expected - 1.0) < 1e-4,
      "rel = {0:.2e}".format(abs(y[idx['Ar*']] / expected - 1.0)))
check("3b. Ar gained == Ar* lost + Ar+ lost (rel < 1e-9)",
      abs((gained - ion_lost) / lost - 1.0) < 1e-9,
      "rel = {0:.2e}".format(abs((gained - ion_lost) / lost - 1.0)))

y0u, yu, nu_mu, tu, idxu, draw_u = run(False)
lost_u = y0u[idxu['Ar*']] - yu[idxu['Ar*']]
print("undeclared: t = {0:.6e} s, nu_m = {1!r}, Ar* {2:.9e} -> {3:.9e}".format(
    tu, nu_mu, y0u[idxu['Ar*']], yu[idxu['Ar*']]))
print("            Ar* lost {0:.6e}; external-source draw on Ar* predicted {1:.6e}".format(
    lost_u, draw_u))
check("4. undeclared arm: nu_m == 0, Ar* loses only the source draw",
      nu_mu == 0.0 and abs(lost_u / draw_u - 1.0) < 1e-2,
      "lost/draw = {0:.6f}".format(lost_u / draw_u))
check("4b. arms differ (declaration reached the code)",
      abs(lost - lost_u) > 0.5 * lost, "declared lost {0:.3e} vs {1:.3e}".format(lost, lost_u))

print("\nRESULT:", "ALL PASS" if not failures else "FAILED: " + ", ".join(failures))
sys.exit(1 if failures else 0)
