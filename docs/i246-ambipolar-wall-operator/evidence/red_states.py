#!/usr/bin/env python3
"""Red/green evidence harness for the I-246 rework.

Run against the *currently built* rmgpy.solver.plasma .so. It exercises every
repair's state and prints, per check, what the loaded module does RIGHT NOW.
Run once before the .pyx changes (banks the red states) and once after the
rebuild (shows them absent). It asserts nothing about which build it is looking
at; it reports observed behaviour and the caller reads the transition.

Usage:
    PYTHONPATH=<worktree> python docs/i246-ambipolar-wall-operator/evidence/red_states.py
"""
import traceback

import numpy as np

import rmgpy.constants as constants
from rmgpy.exceptions import PlasmaStateError
from rmgpy.species import Species
from rmgpy.thermo import ThermoData
from rmgpy.solver.plasma import PlasmaReactor

EV_J_PER_MOL = 96485.33212
TGAS = 298.15
TE_K = 3.0 / 8.617333262e-5
P = 5.0 * 101325.0 / 760.0
LAM = 0.02
MU0 = 1.535e-4


def argon_thermo(excitation_eV=0.0):
    Cp = 2.5 * constants.R
    return ThermoData(
        Tdata=([298, 400, 600, 800, 1000, 1500, 2000], 'K'),
        Cpdata=([Cp] * 7, 'J/(mol*K)'),
        H298=(excitation_eV * EV_J_PER_MOL / 1000.0, 'kJ/mol'),
        S298=(154.8, 'J/(mol*K)'))


def sp(label, adj, thermo_eV=None):
    s = Species(label=label).from_adjacency_list(adj)
    if thermo_eV is not None:
        s.thermo = argon_thermo(thermo_eV)
    return s


def electron():
    return sp('e-', '1 e u1 p0 c-1')


def ground_ar(thermo_eV=None):
    return sp('Ar', '1 Ar u0 p4 c0', thermo_eV)


def meta_ar(label='Ar*', thermo_eV=None):
    return sp(label, '1 Ar u2 p3 c0', thermo_eV)


def ar_cation():
    return sp('Ar+', 'multiplicity 2\n1 Ar u1 p3 c+1')


def build(core, imf, **wall):
    kwargs = dict(diffusion_length=(LAM, 'm'), ion_reduced_mobility=(MU0, 'm^2/(V*s)'))
    kwargs.update(wall)
    r = PlasmaReactor((TGAS, 'K'), (P, 'Pa'), imf, (TE_K, 'K'),
                      n_sims=1, termination=[], **kwargs)
    r.initialize_model(core, [], [], [])
    return r


def _idx(r, want_charge):
    z = r.species_charges
    for j in range(len(z)):
        if j != r.electron_index and z[j] == want_charge:
            return j
    return -1


def report(tag, fn):
    print('=' * 72)
    print(tag)
    try:
        outcome = fn()
        print('  OUTCOME: NO RAISE ->', outcome)
    except PlasmaStateError as exc:
        msg = str(exc)
        # Distinct phrases the NEW messages carry, chosen so the species repr (which
        # embeds "ThermoData(...)") cannot pollute the flag the way a bare "thermo" would.
        flags = {
            'ambiguous(old)': 'is ambiguous' in msg,
            'thermochem': 'no usable thermochemistry' in msg,
            'degeneracy': 'degeneracy threshold' in msg or 'thermal energy k_B' in msg,
            'anion': 'does not support anions' in msg,
            'floor': 'numerical floor' in msg,
            'electron-pop': 'electron population of' in msg,
            'declaration-syntax': 'wallNeutralizationProducts={' in msg,
        }
        print('  OUTCOME: PlasmaStateError |',
              ' '.join('{0}={1}'.format(k, v) for k, v in flags.items() if v) or '(no target phrase)')
        print('  MESSAGE:', ' '.join(msg[:420].split()))
    except Exception as exc:  # noqa
        print('  OUTCOME: OTHER', type(exc).__name__, exc)
        traceback.print_exc()


def h1_headline_bookkeeping():
    e, g, m, c = electron(), ground_ar(0.0), meta_ar('Ar*', 11.5), ar_cation()
    imf = {e: 1e-6, c: 1e-6, m: 1e-6, g: 1.0 - 3e-6}
    r = build([e, g, m, c], imf, wall_recycling=1.0)
    core = [e, g, m, c]
    tgt = int(r.wall_recycle_target[_idx(r, 1)])
    return 'Ar+ recycles to {0!r}'.format(core[tgt].label)


def h2_thermo_absent():
    e, g, m, c = electron(), ground_ar(None), meta_ar('Ar*', None), ar_cation()
    imf = {e: 1e-6, c: 1e-6, m: 1e-6, g: 1.0 - 3e-6}
    build([e, g, m, c], imf, wall_recycling=1.0)
    return 'initialised (picked one silently)'


def h3_degeneracy_outside():
    rt = constants.R * TGAS
    e, g, m, c = electron(), ground_ar(0.0), meta_ar('Ar_hi', 2.0 * rt / EV_J_PER_MOL), ar_cation()
    imf = {e: 1e-6, c: 1e-6, m: 1e-6, g: 1.0 - 3e-6}
    r = build([e, g, m, c], imf, wall_recycling=1.0)
    core = [e, g, m, c]
    tgt = int(r.wall_recycle_target[_idx(r, 1)])
    return 'resolved; Ar+ recycles to {0!r}'.format(core[tgt].label)


def h4_degeneracy_inside():
    rt = constants.R * TGAS
    e, g, m, c = electron(), ground_ar(0.0), meta_ar('Ar_near', 0.5 * rt / EV_J_PER_MOL), ar_cation()
    imf = {e: 1e-6, c: 1e-6, m: 1e-6, g: 1.0 - 3e-6}
    build([e, g, m, c], imf, wall_recycling=1.0)
    return 'initialised (should have refused near-degeneracy)'


def h5_declaration():
    e, g, m, c = electron(), ground_ar(0.0), meta_ar('Ar*', 11.5), ar_cation()
    imf = {e: 1e-6, c: 1e-6, m: 1e-6, g: 1.0 - 3e-6}
    r = build([e, g, m, c], imf, wall_recycling=1.0,
              wall_neutralization_products={'Ar+': 'Ar'})
    core = [e, g, m, c]
    tgt = int(r.wall_recycle_target[_idx(r, 1)])
    return 'declared; Ar+ recycles to {0!r}'.format(core[tgt].label)


def floor_accepts():
    e, g, c = electron(), ground_ar(None), ar_cation()
    imf = {e: 1e-6, c: 1e-6, g: 1.0 - 2e-6}
    r = build([e, g, c], imf, wall_recycling=0.0)
    ie, ig, ic = r.electron_index, _idx(r, 0), _idx(r, 1)
    y = np.zeros(r.num_core_species, float)
    y[ig] = 1.0e-9
    y[ic] = 1.0e-15
    y[ie] = 1.0e-15
    r.check_wall_support(y)
    return 'check_wall_support ACCEPTED a sub-floor state (floor={0:.3g})'.format(
        r.wall_neutral_floor)


def anion_accepts():
    e, g = electron(), ground_ar(None)
    cl = sp('Cl', '1 Cl u1 p3 c0')
    cln = sp('Cl-', '1 Cl u0 p4 c-1')
    c = ar_cation()
    imf = {e: 1e-6, c: 1e-6, cln: 1e-6, cl: 1e-6, g: 1.0 - 4e-6}
    build([e, g, cl, cln, c], imf, wall_recycling=1.0)
    return 'anion Cl- admitted; it receives the Ar+ mobility'


def electron_nan_accepts():
    e, g, c = electron(), ground_ar(None), ar_cation()
    imf = {e: 1e-6, c: 1e-6, g: 1.0 - 2e-6}
    r = build([e, g, c], imf, wall_recycling=0.0)
    ie, ig = r.electron_index, _idx(r, 0)
    y = np.zeros(r.num_core_species, float)
    y[ig] = 1.0
    y[ie] = float('nan')
    r.check_wall_support(y)
    return 'check_wall_support ACCEPTED a NaN electron population'


if __name__ == '__main__':
    report('H1  deliverable deck (Ar, Ar*, Ar+, e-) wall_recycling=1.0', h1_headline_bookkeeping)
    report('H2  two neutral Ar, NO thermo  (must name thermo)', h2_thermo_absent)
    report('H3  gap = 2*R*T  (must RESOLVE to ground Ar)', h3_degeneracy_outside)
    report('H4  gap = 0.5*R*T (must REFUSE near-degeneracy)', h4_degeneracy_inside)
    report('H5  wallNeutralizationProducts={"Ar+":"Ar"} declaration', h5_declaration)
    report('F   sub-floor accepted state (must be refused)', floor_accepts)
    report('A   anion in the core with a wall (must be refused)', anion_accepts)
    report('E   NaN electron population (must be refused)', electron_nan_accepts)
    print('=' * 72)
