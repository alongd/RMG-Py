#!/usr/bin/env python3
"""Charge and particle conservation across the wall operator, metastable present.

Builds the actual collision deck -- ground Ar, metastable Ar*, Ar+, e- -- with
wall_recycling = 1.0, resolves Ar+ -> ground Ar by the energy rule, then reads the
per-species wall flux the reactor records in wall_loss_rates and checks:

  * net charge rate  = sum_j z_j * wall_flux_j  == 0  (a common nu removes charge at
    rate nu*(net charge), which is zero in a neutral gas: no charge made at the wall);
  * net heavy rate   = sum over heavy species of wall_flux_j == 0 at gamma = 1 (every
    neutralised Ar+ returns as ground Ar; the metastable is untouched);
  * at gamma = 0.5 the pumping wall removes exactly half the heavy flux while charge
    is still conserved.

Also prints the mole-fraction bookkeeping: where the recycled Ar+ population lands.
"""
import numpy as np

import rmgpy.constants as constants
from red_states import electron, ground_ar, meta_ar, ar_cation, build, _idx


def deck(gamma):
    e, g, m, c = electron(), ground_ar(0.0), meta_ar('Ar*', 11.5), ar_cation()
    imf = {e: 1e-6, c: 1e-6, m: 1e-4, g: 1.0 - 1e-6 - 1e-6 - 1e-4}
    r = build([e, g, m, c], imf, wall_recycling=gamma)
    return r, [e, g, m, c]


def run(gamma):
    r, core = deck(gamma)
    z = r.species_charges
    ie = r.electron_index
    ig = _idx(r, 0)          # ground Ar (lowest-enthalpy neutral)
    ic = _idx(r, 1)          # Ar+
    # metastable index: the neutral that is not ground
    i_meta = [j for j in range(len(z)) if z[j] == 0 and j != ig][0]
    tgt = int(r.wall_recycle_target[ic])

    # a charge-neutral state with the metastable populated
    y = np.zeros(r.num_core_species, float)
    y[ig] = 1.0
    y[i_meta] = 1.0e-2
    y[ic] = 1.0e-4
    y[ie] = 1.0e-4          # net charge zero
    dydt = np.zeros(r.num_core_species, float)
    r.residual(0.0, y.copy(), dydt.copy())
    wall = np.asarray(r.wall_loss_rates, float)

    net_charge_rate = float(sum(z[j] * wall[j] for j in range(len(z))))
    heavy = [j for j in range(len(z)) if j != ie]
    net_heavy_rate = float(sum(wall[j] for j in heavy))
    nu = r.nu_wall

    print('=' * 72)
    print('gamma = {0}'.format(gamma))
    print('  Ar+ recycle target      : {0!r}  (ground state)'.format(core[tgt].label))
    print('  nu_wall                 : {0:.6g} s^-1'.format(nu))
    print('  wall flux, ground Ar    : {0:+.6e} mol/s'.format(wall[ig]))
    print('  wall flux, metastable   : {0:+.6e} mol/s   (must be 0: neutral)'.format(wall[i_meta]))
    print('  wall flux, Ar+          : {0:+.6e} mol/s'.format(wall[ic]))
    print('  wall flux, e-           : {0:+.6e} mol/s'.format(wall[ie]))
    print('  net CHARGE rate         : {0:+.3e} mol/s   (must be ~0)'.format(net_charge_rate))
    print('  net HEAVY rate          : {0:+.3e} mol/s'.format(net_heavy_rate))
    expect_heavy = -(1.0 - gamma) * nu * y[ic]
    print('  expected heavy rate     : {0:+.3e} mol/s   -(1-gamma)*nu*y_Ar+'.format(expect_heavy))
    # ground Ar should receive gamma*nu*y_Ar+ ; Ar+ should lose nu*y_Ar+
    print('  Ar+ loss  == -nu*y_Ar+  : {0}'.format(
        np.isclose(wall[ic], -nu * y[ic], rtol=1e-12, atol=0.0)))
    print('  ground gain== +g*nu*yAr+: {0}'.format(
        np.isclose(wall[ig], gamma * nu * y[ic], rtol=1e-12, atol=0.0)))
    print('  metastable untouched    : {0}'.format(wall[i_meta] == 0.0))
    print('  net charge conserved    : {0}'.format(abs(net_charge_rate) <= 1e-18 * max(1.0, abs(nu))))
    print('  net heavy == expected   : {0}'.format(
        np.isclose(net_heavy_rate, expect_heavy, rtol=1e-12, atol=1e-30)))


if __name__ == '__main__':
    run(1.0)
    run(0.5)
    print('=' * 72)
