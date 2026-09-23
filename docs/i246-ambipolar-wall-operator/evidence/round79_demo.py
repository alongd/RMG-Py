#!/usr/bin/env python
"""Round-79 demonstration: the skeleton key, the full source delivery, the latched
wall interface and its availability map. Prints numbers a reviewer can read against
the reproduction tests. Run with rmg_env on PATH, PYTHONPATH=worktree, MPLCONFIGDIR set.
"""
import numpy as np

import rmgpy.constants as constants
from rmgpy.molecule import Molecule
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.species import Species
from rmgpy.thermo import ThermoData

EV_TO_K = 1.0 / 8.617333262e-5
TORR = 101325.0 / 760.0
TGAS, P, TE_EV = 298.15, 5.0 * TORR, 3.0
MU0 = 1.535e-4
LAM = 1.0 / np.sqrt((2.405 / 0.05) ** 2 + (np.pi / 0.30) ** 2)


def th(h_kj):
    return ThermoData(Tdata=([298, 400, 600, 800, 1000, 1500, 2000], 'K'),
                      Cpdata=([4.0 * constants.R] * 7, 'J/(mol*K)'),
                      H298=(h_kj, 'kJ/mol'), S298=(200.0, 'J/(mol*K)'))


def reactor(core, imf, **kw):
    r = PlasmaReactor((TGAS, 'K'), (P, 'Pa'), imf, (TE_EV * EV_TO_K, 'K'),
                      n_sims=1, termination=[],
                      diffusion_length=(LAM, 'm'), ion_reduced_mobility=(MU0, 'm^2/(V*s)'), **kw)
    r.initialize_model(core, [], [], [])
    return r


print("== 1. the skeleton key sits between formula and is_isomorphic ==")
Ar = Molecule().from_adjacency_list('1 Ar u0 p4 c0')
Ars = Molecule().from_adjacency_list('multiplicity 3\n1 Ar u2 p3 c0')
dme = Molecule().from_smiles('COC')
eth = Molecule().from_smiles('CCO')


def skel(m):
    s = m.to_inchi()
    for sep in ('/q', '/p'):
        i = s.find(sep)
        if i != -1:
            s = s[:i]
    return s


print("  Ar   skeleton:", skel(Ar))
print("  Ar*  skeleton:", skel(Ars), " -> same as Ar (electronic states unified)")
print("  DME  skeleton:", skel(dme))
print("  EtOH skeleton:", skel(eth), " -> differs from DME (isomers separated)")

print("\n== 2. HIGH 1: DME+ neutralises to DME, not the lower-enthalpy isomer ethanol ==")
e = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
d = Species(label='DME').from_smiles('COC'); d.thermo = th(-184.0)
et = Species(label='EtOH').from_smiles('CCO'); et.thermo = th(-235.0)
dp = Species(label='DME+').from_smiles('[CH3][O+][CH3]')
r = reactor([e, d, et, dp], {e: 1e-6, dp: 1e-6, et: 1e-6, d: 1 - 3e-6}, wall_recycling=1.0)
z = r.species_charges
ic = [j for j in range(len(z)) if z[j] == 1 and j != r.electron_index][0]
print("  ethanol H298 < DME H298, yet DME+ recycles to:",
      [e, d, et, dp][int(r.wall_recycle_target[ic])].label)

print("\n== 3. HIGH 3: declared source delivered in FULL with a He bath gas ==")
he = Species(label='He').from_adjacency_list('1 He u0 p1 c0')
ar = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
S = 1.0e18
r = reactor([e, ar, he, arp], {e: 1e-6, arp: 1e-6, he: 0.5, ar: 0.5 - 2e-6},
            wall_recycling=0.0, ionisation_source=(S, 'm^-3/s'))
z = r.species_charges
ie = r.electron_index
i_ar = [j for j in range(len(z)) if z[j] == 0 and [e, ar, he, arp][j].label == 'Ar'][0]
i_he = [j for j in range(len(z)) if z[j] == 0 and [e, ar, he, arp][j].label == 'He'][0]
y = np.zeros(r.num_core_species, float); y[i_ar] = 0.5; y[i_he] = 0.5
V = r.compute_volume(y)
delta, _ = r.residual(0.0, y.copy(), np.zeros_like(y))
print("  declared S*V/Na = {0:.6e} mol/s;  delivered to e- = {1:.6e} mol/s;  ratio = {2:.4f}"
      .format(S * V / constants.Na, delta[ie], delta[ie] / (S * V / constants.Na)))

print("\n== 4. HIGH 5: latched wall interface + availability map (Ar/Ar*/Ar+/e-, gamma=1) ==")
ar2 = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0'); ar2.thermo = th(0.0)
ars = Species(label='Ar*').from_adjacency_list('1 Ar u2 p3 c0'); ars.thermo = th(1110.0)
arp2 = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
# Round 83 note: this Ar/Ar* deck no longer auto-resolves -- two electronic states share
# Ar+'s skeleton and nothing certifies the ground, so a declaration is now REQUIRED (see
# round83). Round 79 auto-picked ground Ar here; that pick is exactly what round 83 refuses.
r = reactor([e, ar2, ars, arp2], {e: 1e-6, arp2: 1e-6, ars: 1e-6, ar2: 1 - 3e-6},
            wall_recycling=1.0, wall_neutralization_products={'Ar+': 'Ar'})
print("  nu_wall_latched      = {0:.4f} s^-1".format(r.nu_wall_latched))
print("  wall_flux (mol/s)    =", np.array2string(np.array(r.wall_flux), precision=3))
print("  electron energy flux = {0:.4e} W   [{1}]".format(
    r.wall_electron_energy_flux, r.wall_energy_availability['wall_electron_energy_flux']))
print("  neutralisation flux  = {0!r}   [{1}]".format(
    r.wall_neutralization_energy_flux, r.wall_energy_availability['wall_neutralization_energy_flux']))
print("  ion directed flux    = {0!r}   [{1}]  <- power balance owns this".format(
    r.wall_ion_energy_flux, r.wall_energy_availability['wall_ion_energy_flux']))

print("\n== 5. a rejected Newton trial cannot leak into the latch ==")
latched = np.array(r.wall_flux, float)
r.residual(0.0, np.array(r.y0, float) * 1e6, np.zeros(r.num_core_species))
print("  wall_flux unchanged after a wild residual:", np.array_equal(np.array(r.wall_flux, float), latched))
print("  residual scratch (wall_loss_rates) moved: ",
      not np.array_equal(np.array(r.wall_loss_rates, float), latched))
