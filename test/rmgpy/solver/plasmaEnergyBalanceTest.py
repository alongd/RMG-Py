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

"""
Tests for the electron energy balance on :class:`PlasmaReactor` (I-274, M8-B).

With ``electron_energy_balance`` declared, Te is a solved state variable:

    d(3/2 N_e R Te)/dt = P_abs - Q_inelastic - Q_elastic - Q_wall - Q_flow

Every loss term is checked against a value computed here by hand from its own
formula, never against the reactor's output for another term. The integration
tests use a toy argon whose particle-balance root Te* can be solved in closed
form, so "Te is set by particle balance, n_e by power" is checked exactly.

No assertion fits an absorbed power to an electron density.
"""

import copy
import pickle

import numpy as np
import pytest
from scipy.optimize import brentq

import rmgpy.constants as constants
from rmgpy.exceptions import PlasmaStateError
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.kinetics import TwoTemperaturePlasma
from rmgpy.reaction import Reaction
from rmgpy.solver.plasma import PLASMA_LOSCHMIDT, PlasmaReactor
from rmgpy.species import Species
from rmgpy.thermo import ThermoData

EV_TO_K = 1.0 / 8.617333262e-5           # K per eV
EV_J_PER_MOL = 96485.33212               # J/mol per eV
# The hand values use CODATA-2018 eV/K conversions; the engine converts through rmgpy.constants
# and the quantities package, an older vintage. Hand-vs-engine agreement is therefore held to
# CONST_RTOL, the size of that gap propagated through exp(-E/kTe), not to machine precision.
CONST_RTOL = 5e-5
TORR_TO_PA = 101325.0 / 760.0
TGAS = 298.15
P_NOMINAL = 5.0 * TORR_TO_PA
MU0_AR_IN_AR = 1.535e-4
RADIUS, LENGTH = 0.05, 0.30
LAMBDA = 1.0 / np.sqrt((2.405 / RADIUS) ** 2 + (np.pi / LENGTH) ** 2)
CHAMBER_VOLUME = np.pi * RADIUS ** 2 * LENGTH

E_IZ_EV, E_EX_EV = 15.76, 11.55
# Toy rate laws k = A exp(-E/kTe) (TwoTemperaturePlasma with Ea_g = Ea_e = E, n = 0).
A_IZ = 2.34e-14                           # m^3/s
A_EX = 5.0e-15                            # m^3/s
A_SE = 4.3e-16                            # m^3/s, superelastic Ar* + e -> Ar + e
A_RC = 1.0e-18                            # m^3/s, Ar+ + e -> Ar (radiative, toy)

# Lieberman & Lichtenberg (2005) Table 3.3, argon elastic scattering:
#     K_el = 2.336e-14 Te^1.609 exp(0.0618 (ln Te)^2 - 0.1171 (ln Te)^3)  m^3/s, Te in eV
LL_AR_ELASTIC = {'A': (2.336e-14, 'm^3/s'), 'n': 1.609, 'b': 0.0618, 'c': -0.1171}


def _thermo(h_ev):
    Cp = 2.5 * constants.R
    return ThermoData(Tdata=([298, 400, 600, 800, 1000, 1500, 2000], 'K'),
                      Cpdata=([Cp] * 7, 'J/(mol*K)'),
                      H298=(h_ev * EV_J_PER_MOL / 1000.0, 'kJ/mol'),
                      S298=(154.8, 'J/(mol*K)'))


def _two_temp(a_m3s, e_ev):
    e = (e_ev, 'eV/molecule')
    return TwoTemperaturePlasma(A=(a_m3s, 'm^3/(molecule*s)'), n=0.0, Ea_g=e, Ea_e=e,
                                T0=(1.0, 'K'))


def _k(a_m3s, e_ev, te_k):
    """Hand rate coefficient, m^3/(mol s)."""
    return a_m3s * constants.Na * np.exp(-e_ev * EV_TO_K / te_k)


# The toy library's declared electron energies, eV per event. 'Toy:5' is a lumped
# metastable-to-resonance mixing proxy written Ars + e- => Ar + e-: its thermo dH is
# -11.55 eV but the electron pays the m->r gap, +0.076 eV -- only a declaration says so.
TOY_ENERGIES = {'Toy:1': (E_IZ_EV, 'eV'), 'Toy:2': (E_EX_EV, 'eV'),
                'Toy:3': (-E_EX_EV, 'eV'), 'Toy:4': (0.0, 'eV'), 'Toy:5': (0.076, 'eV')}


class ToyLibraryReactor(PlasmaReactor):
    """Resolves 'Toy:<index>' against the toy reactions built below instead of a loaded
    database; everything else is the production reactor."""
    toy_entries = {}

    def _declared_entry_kinetics(self, library, index):
        if library != 'Toy' or index not in self.toy_entries:
            raise PlasmaStateError('no entry {0}:{1}'.format(library, index))
        return self.toy_entries[index]


def _energy_decl(power_w=0.5, elastic=True, volume=CHAMBER_VOLUME, energies=None):
    return {'absorbed_power': (power_w, 'W'),
            'chamber_volume': (volume, 'm^3'),
            'sheath': 'floating_wall',
            'electron_energies': dict(TOY_ENERGIES if energies is None else energies),
            'elastic_collisions': {'Ar': dict(LL_AR_ELASTIC)} if elastic else {}}


def _toy(index, reactants, products, kinetics, entries):
    entries[index] = kinetics
    return LibraryReaction(reactants=reactants, products=products, reversible=False,
                           kinetics=kinetics, library='Toy')


def _build(te_ev=1.0, x_ion=1.0e-6, power_w=0.5, energy=True, metastable=False,
           recombination=False, elastic=True, source=None, radius=RADIUS, mixing=False,
           energies=None):
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    ar.thermo, arp.thermo = _thermo(0.0), _thermo(E_IZ_EV)
    core = [electron, ar, arp]
    entries = {}
    rxns = [_toy(1, [electron, ar], [arp, electron, electron], _two_temp(A_IZ, E_IZ_EV), entries)]
    imf = {electron: x_ion, arp: x_ion, ar: 1.0 - 2.0 * x_ion}
    if metastable:
        ars = Species(label='Ars').from_adjacency_list('1 Ar u2 p3 c0')
        ars.thermo = _thermo(E_EX_EV)
        core.append(ars)
        imf[ars] = 1.0e-7
        imf[ar] -= 1.0e-7
        rxns.append(_toy(2, [electron, ar], [ars, electron], _two_temp(A_EX, E_EX_EV), entries))
        rxns.append(_toy(3, [electron, ars], [ar, electron], _two_temp(A_SE, 0.0), entries))
        if mixing:
            rxns.append(_toy(5, [electron, ars], [ar, electron], _two_temp(100 * A_SE, 0.0), entries))
    if recombination:
        rxns.append(_toy(4, [arp, electron], [ar], _two_temp(A_RC, 0.0), entries))
    lam = 1.0 / np.sqrt((2.405 / radius) ** 2 + (np.pi / LENGTH) ** 2)
    kwargs = dict(diffusion_length=(lam, 'm'), ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'))
    if metastable:
        kwargs['wall_neutralization_products'] = {'Ar+': 'Ar'}
    if source is not None:
        kwargs['ionisation_source'] = (source, 'm^-3/s')
    if energy:
        kwargs['electron_energy_balance'] = _energy_decl(
            power_w, elastic, volume=np.pi * radius ** 2 * LENGTH,
            energies={k: v for k, v in (TOY_ENERGIES if energies is None else energies).items()
                      if int(k[4:]) in entries})
    ToyLibraryReactor.toy_entries = entries
    reactor = ToyLibraryReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf, (te_ev * EV_TO_K, 'K'),
                            n_sims=1, termination=[], **kwargs)
    reactor.initialize_model(core, rxns, [], [])
    return reactor, core, rxns


def _idx(reactor, core, label):
    return reactor.species_index[[s for s in core if s.label == label][0]]


def _nu_wall_hand(te_k, n_neutral, lam=LAMBDA):
    mu_i = MU0_AR_IN_AR * PLASMA_LOSCHMIDT / n_neutral
    return mu_i * (constants.R / constants.Na) * te_k / constants.e / lam ** 2


def _sheath_factor(mass_kg):
    """Ion energy at a floating wall per lost ion, in units of k_B Te:
    Te/2 (Bohm presheath) + (Te/2) ln(M / (2 pi m_e)) (floating sheath)."""
    return 0.5 + 0.5 * np.log(mass_kg / (2.0 * np.pi * constants.m_e))


# ------------------------------------------------------------------ state layout


def test_energy_balance_off_keeps_the_species_only_state():
    """Default path: no declaration, no Te slot; neq is the species count."""
    r, core, _ = _build(energy=False)
    assert not r.energy_balance
    assert r.neq == r.num_core_species
    assert len(r.y0) == r.num_core_species


def test_te_is_a_state_variable_when_declared():
    """Clause 1: Te is carried in the integrated state, after the species."""
    r, core, _ = _build(te_ev=1.3)
    assert r.energy_balance
    assert r.neq == r.num_core_species + 1
    assert r.te_index == r.num_core_species
    assert r.y0[r.te_index] == pytest.approx(1.3 * EV_TO_K, rel=1e-14)
    assert r.absorbed_power_density == pytest.approx(0.5 / CHAMBER_VOLUME, rel=1e-14)


# ------------------------------------------------------------------ loss terms


def _rxn_index(reactor, rxns, reactants, products):
    for rxn in rxns:
        if (sorted(s.label for s in rxn.reactants) == sorted(reactants)
                and sorted(s.label for s in rxn.products) == sorted(products)):
            return reactor.reaction_index[rxn]
    raise KeyError((reactants, products))


def test_electron_energies_are_declared_and_thermo_is_only_a_cross_check():
    """Q_inelastic energy per reaction is the DECLARED value (ionisation +15.76, excitation
    +11.55, superelastic -11.55 eV); the thermo dH is carried beside it as a cross-check and
    agrees here. A reaction that consumes an electron (recombination, declared 0) is charged
    the consumed electron's mean energy 3/2 R Te on top."""
    r, core, rxns = _build(metastable=True, recombination=True)
    cross = np.asarray(r.energy_threshold_thermo, float) / EV_J_PER_MOL
    thr = np.asarray(r.energy_threshold, float) / EV_J_PER_MOL
    consumed = np.asarray(r.energy_electrons_consumed, int)
    iz = _rxn_index(r, rxns, ['e-', 'Ar'], ['Ar+', 'e-', 'e-'])
    ex = _rxn_index(r, rxns, ['e-', 'Ar'], ['Ars', 'e-'])
    se = _rxn_index(r, rxns, ['e-', 'Ars'], ['Ar', 'e-'])
    rc = _rxn_index(r, rxns, ['Ar+', 'e-'], ['Ar'])
    assert thr[iz] == pytest.approx(E_IZ_EV, rel=CONST_RTOL)
    assert thr[ex] == pytest.approx(E_EX_EV, rel=CONST_RTOL)
    assert thr[se] == pytest.approx(-E_EX_EV, rel=CONST_RTOL)
    assert thr[rc] == 0.0 and consumed[rc] == 1
    assert consumed[iz] == consumed[ex] == consumed[se] == 0
    for j, expect in ((iz, E_IZ_EV), (ex, E_EX_EV), (se, -E_EX_EV), (rc, -E_IZ_EV)):
        assert cross[j] == pytest.approx(expect, rel=CONST_RTOL)


def test_lumped_mixing_proxy_takes_its_declared_energy_not_its_enthalpy():
    """PM correction 2026-09-24: 'Ars + e- => Ar + e-' written for m->r mixing followed by
    radiation costs the electron the m->r gap (+0.076 eV). Thermo says -11.55 eV (a gain);
    the balance must use the declaration, and the thermo value stays a cross-check."""
    r, core, rxns = _build(metastable=True, mixing=True)
    j = r.energy_reaction_keys.index('Toy:5')
    assert r.energy_threshold[j] / EV_J_PER_MOL == pytest.approx(0.076, rel=CONST_RTOL)
    assert r.energy_threshold_thermo[j] / EV_J_PER_MOL == pytest.approx(-E_EX_EV, rel=CONST_RTOL)
    assert r.energy_reaction_keys[j] == 'Toy:5'


def test_an_undeclared_electron_reaction_is_refused_naming_its_thermo_value():
    energies = dict(TOY_ENERGIES)
    del energies['Toy:2']
    with pytest.raises(PlasmaStateError, match=r'undeclared: .*Ars.*thermo dH = 11\.55'):
        _build(metastable=True, energies=energies)


def test_inelastic_loss_hand_value():
    """Q_inelastic = sum_j r_j eps_j, hand-computed from each toy rate law."""
    te_k = 1.0 * EV_TO_K
    r, core, rxns = _build(te_ev=1.0, metastable=True, recombination=True)
    y = np.array(r.y0, float)
    r.residual(0.0, y, np.zeros_like(y))
    V = r.compute_volume(y)
    C = y[:r.num_core_species] / V
    ie, iar, iarp, iars = (_idx(r, core, l) for l in ('e-', 'Ar', 'Ar+', 'Ars'))
    R = constants.R
    hand = V * C[ie] * (_k(A_IZ, E_IZ_EV, te_k) * C[iar] * E_IZ_EV * EV_J_PER_MOL
                        + _k(A_EX, E_EX_EV, te_k) * C[iar] * E_EX_EV * EV_J_PER_MOL
                        - _k(A_SE, 0.0, te_k) * C[iars] * E_EX_EV * EV_J_PER_MOL
                        + _k(A_RC, 0.0, te_k) * C[iarp] * 1.5 * R * te_k)
    assert r.electron_energy_terms['Q_inelastic'] == pytest.approx(hand, rel=CONST_RTOL)


def test_elastic_loss_hand_value_at_one_ev():
    """Q_elastic = 3 (m_e/M) K_el(Te) n_Ar N_e k_B (Te - Tg). At Te = 1 eV the L&L fit
    reduces to K_el = A = 2.336e-14 m^3/s exactly (ln 1 = 0): a literal hand value."""
    te_k = 1.0 * EV_TO_K
    r, core, _ = _build(te_ev=1.0)
    y = np.array(r.y0, float)
    r.residual(0.0, y, np.zeros_like(y))
    V = r.compute_volume(y)
    ie, iar = _idx(r, core, 'e-'), _idx(r, core, 'Ar')
    m_ar = core[[s.label for s in core].index('Ar')].molecular_weight.value_si
    n_ar = y[iar] * constants.Na / V
    hand = 3.0 * constants.m_e / m_ar * 2.336e-14 * n_ar * y[ie] * constants.R * (te_k - TGAS)
    assert r.electron_energy_terms['Q_elastic'] == pytest.approx(hand, rel=CONST_RTOL)


def test_elastic_fit_hand_value_at_0p9_ev():
    """The L&L fit away from its trivial point: K_el(0.9 eV) by hand is 1.9734e-14 m^3/s
    (2.336e-14 * 0.9^1.609 * exp(0.0618*ln(0.9)^2 - 0.1171*ln(0.9)^3))."""
    te_k = 0.9 * EV_TO_K
    r, core, _ = _build(te_ev=0.9)
    y = np.array(r.y0, float)
    r.residual(0.0, y, np.zeros_like(y))
    V = r.compute_volume(y)
    ie, iar = _idx(r, core, 'e-'), _idx(r, core, 'Ar')
    m_ar = core[[s.label for s in core].index('Ar')].molecular_weight.value_si
    n_ar = y[iar] * constants.Na / V
    per = r.electron_energy_terms['Q_elastic'] / (
        3.0 * constants.m_e / m_ar * n_ar * y[ie] * constants.R * (te_k - TGAS))
    assert per == pytest.approx(1.9734e-14, rel=1e-4)


def test_wall_loss_hand_value_uses_the_particle_balance_flux():
    """Q_wall = nu N_e 2 R Te + nu N_Ar+ R Te (1/2 + 1/2 ln(M/2 pi m_e)), with nu the SAME
    number the particle balance applies (clause 8): the electron wall-loss scratch of the
    same residual evaluation is -nu N_e."""
    te_k = 1.0 * EV_TO_K
    r, core, _ = _build(te_ev=1.0)
    y = np.array(r.y0, float)
    r.residual(0.0, y, np.zeros_like(y))
    V = r.compute_volume(y)
    ie, iar, iarp = (_idx(r, core, l) for l in ('e-', 'Ar', 'Ar+'))
    n_ar = y[iar] * constants.Na / V
    nu = _nu_wall_hand(te_k, n_ar)
    assert r.nu_wall == pytest.approx(nu, rel=1e-9)
    assert -r.wall_loss_rates[ie] == pytest.approx(nu * y[ie], rel=1e-12)
    m_arp = core[iarp].molecular_weight.value_si
    R = constants.R
    assert r.electron_energy_terms['Q_wall_electron'] == pytest.approx(
        2.0 * R * te_k * nu * y[ie], rel=1e-9)
    assert r.electron_energy_terms['Q_wall_ion'] == pytest.approx(
        _sheath_factor(m_arp) * R * te_k * nu * y[iarp], rel=1e-9)
    # the energy term is built from the identical flux the species row uses
    assert r.electron_energy_terms['Q_wall_electron'] == pytest.approx(
        2.0 * R * te_k * (-r.wall_loss_rates[ie]), rel=1e-15)


def test_flow_loss_is_zero_in_the_closed_batch():
    r, core, _ = _build()
    y = np.array(r.y0, float)
    r.residual(0.0, y, np.zeros_like(y))
    assert r.electron_energy_terms['Q_flow'] == 0.0


def test_energy_row_conserves_electron_energy():
    """The Te row is d(3/2 N_e R Te)/dt = P_abs - sum Q, written for Te: with the
    electron's own net rate dN_e/dt from the same residual,
    3/2 R (N_e dTe/dt + Te dN_e/dt) equals P_abs - sum Q."""
    r, core, _ = _build(te_ev=1.2, metastable=True, recombination=True)
    y = np.array(r.y0, float)
    dydt = np.zeros_like(y)
    delta, _ = r.residual(0.0, y, dydt)
    t = r.electron_energy_terms
    ie = r.electron_index
    te = y[r.te_index]
    lhs = 1.5 * constants.R * (y[ie] * delta[r.te_index] + te * delta[ie])
    losses = t['Q_inelastic'] + t['Q_elastic'] + t['Q_wall_electron'] + t['Q_wall_ion'] + t['Q_flow']
    assert t['P_abs'] == pytest.approx(r.absorbed_power_density * r.compute_volume(y), rel=1e-14)
    assert lhs == pytest.approx(t['P_abs'] - losses, rel=1e-9, abs=1e-12 * t['P_abs'])


def test_energy_mode_jacobian_matches_residual_differences():
    """The energy-mode Jacobian (species + Te) agrees with a central difference of the
    residual, including the Te column that the fixed-Te analytic Jacobian lacks."""
    r, core, _ = _build(te_ev=1.0, metastable=True, recombination=True)
    y = np.array(r.y0, float)
    dydt = np.zeros_like(y)
    pd = np.array(r.jacobian(0.0, y, dydt, 0.0), float)
    col = r.te_index
    h = 1e-6 * y[col]
    yp, ym = y.copy(), y.copy()
    yp[col] += h
    ym[col] -= h
    fd = (r.residual(0.0, yp, dydt)[0] - r.residual(0.0, ym, dydt)[0]) / (2 * h)
    assert np.allclose(pd[:, col], fd, rtol=1e-4, atol=1e-12 * np.max(np.abs(fd)))
    assert pd[_idx(r, core, 'Ar+'), col] > 0.0          # hotter electrons ionise faster


def test_reduce_round_trips_the_energy_balance():
    r, core, _ = _build()
    for clone in (copy.deepcopy(r), pickle.loads(pickle.dumps(r))):
        assert clone.electron_energy_balance == r.electron_energy_balance


def test_energy_balance_requires_a_wall():
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    with pytest.raises(PlasmaStateError, match='wall'):
        PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), {electron: 1e-6, ar: 1 - 1e-6},
                      (EV_TO_K, 'K'), electron_energy_balance=_energy_decl())


@pytest.mark.parametrize('bad, match', [
    ({'absorbed_power': (-1.0, 'W')}, 'absorbed_power'),
    ({'absorbed_power': (1.0, 'm')}, 'absorbed_power'),
    ({'sheath': 'collisional'}, 'sheath'),
    ({'discharge_current': (1.0, 'A')}, 'discharge_current'),
    ({'elastic_collisions': {'Ar': {'A': (1e-14, 'm^3/s'), 'n': 1.0}}}, 'elastic'),
])
def test_energy_balance_declaration_is_validated(bad, match):
    decl = _energy_decl()
    decl.update(bad)
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    with pytest.raises(PlasmaStateError, match=match):
        PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), {electron: 1e-6, ar: 1 - 1e-6},
                      (EV_TO_K, 'K'), diffusion_length=(LAMBDA, 'm'),
                      ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                      electron_energy_balance=decl)


def test_undeclared_elastic_partner_is_refused_at_initialisation():
    decl = _energy_decl(energies={})
    decl['elastic_collisions'] = {'He': dict(LL_AR_ELASTIC)}
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    ar.thermo, arp.thermo = _thermo(0.0), _thermo(E_IZ_EV)
    r = PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), {electron: 1e-6, arp: 1e-6, ar: 1 - 2e-6},
                      (EV_TO_K, 'K'), diffusion_length=(LAMBDA, 'm'),
                      ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                      electron_energy_balance=decl)
    with pytest.raises(PlasmaStateError, match="'He'"):
        r.initialize_model([electron, ar, arp], [], [], [])


# ------------------------------------------------------------------ integration


def _te_star_hand(radius=RADIUS):
    """Particle-balance root of the toy: k_iz(Te) n_Ar = nu_wall(Te)."""
    n_ar = P_NOMINAL / (constants.R / constants.Na * TGAS)
    lam = 1.0 / np.sqrt((2.405 / radius) ** 2 + (np.pi / LENGTH) ** 2)
    f = lambda te: (A_IZ * np.exp(-E_IZ_EV * EV_TO_K / te) * n_ar
                    - _nu_wall_hand(te, n_ar, lam))
    return brentq(f, 0.3 * EV_TO_K, 3.0 * EV_TO_K, xtol=1e-10)


def _run(r, t_end=40.0):
    for t in np.logspace(-9, np.log10(t_end), 120):
        r.advance(t)
    return r


def _n_e(r, y=None):
    y = r.y if y is None else y
    return y[r.electron_index] * constants.Na / r.compute_volume(y)


def test_steady_te_is_the_particle_balance_root_from_two_initial_temperatures():
    """Clauses 1, 9 and verifier item 2: two initial Te values reach the same steady Te,
    and that Te is the closed-form particle-balance root (the energy equation does not
    overturn it)."""
    te_star = _te_star_hand()
    finals = []
    for te0 in (0.6, 2.0):
        r, core, _ = _build(te_ev=te0, x_ion=1e-9, power_w=0.5)
        _run(r)
        finals.append((r.y[r.te_index], _n_e(r)))
    for te, _ in finals:
        assert te == pytest.approx(te_star, rel=1e-5)
    assert finals[0][1] == pytest.approx(finals[1][1], rel=1e-4)


def test_steady_electron_density_is_linear_in_power_and_matches_the_global_formula():
    """Clauses 2 and 3: n_e scales with P_abs at fixed wall, and one point matches
    P_abs/V_ch = n_e nu_wall (E_c + E_e + E_i) with every term by hand."""
    dens = {}
    for p in (0.05, 0.5):
        r, core, _ = _build(te_ev=1.0, x_ion=1e-9, power_w=p)
        _run(r)
        dens[p] = _n_e(r)
    assert dens[0.5] / dens[0.05] == pytest.approx(10.0, rel=1e-4)

    te = _te_star_hand()
    n_ar = P_NOMINAL / (constants.R / constants.Na * TGAS)
    nu_w = _nu_wall_hand(te, n_ar)
    kB = constants.R / constants.Na
    m_ar = Species().from_adjacency_list('1 Ar u0 p4 c0').molecular_weight.value_si
    te_ev = te / EV_TO_K
    k_el = 2.336e-14 * te_ev ** 1.609 * np.exp(0.0618 * np.log(te_ev) ** 2 - 0.1171 * np.log(te_ev) ** 3)
    # per ionisation (= per wall loss at steady state, nu_iz = nu_w):
    e_c = E_IZ_EV * constants.e + 3 * constants.m_e / m_ar * k_el * n_ar * kB * (te - TGAS) / nu_w
    e_e = 2 * kB * te
    e_i = _sheath_factor(m_ar) * kB * te
    n_e_hand = 0.5 / CHAMBER_VOLUME / (nu_w * (e_c + e_e + e_i))
    assert dens[0.5] == pytest.approx(n_e_hand, rel=1e-3)


def test_geometry_moves_te_star_through_the_wall_term():
    """Clause 4: a larger chamber loses less to the wall, so Te* falls (and rises for a
    smaller one) -- the steady Te tracks the hand root at each radius."""
    tes = {}
    for radius in (0.5 * RADIUS, 1.5 * RADIUS):
        r, core, _ = _build(te_ev=1.0, x_ion=1e-9, radius=radius)
        _run(r)
        tes[radius] = r.y[r.te_index]
        assert tes[radius] == pytest.approx(_te_star_hand(radius), rel=1e-5)
    assert tes[0.5 * RADIUS] > tes[1.5 * RADIUS]


def test_zero_power_is_reported_as_extinction():
    """Clause 5: with no absorbed power the electrons cool to the gas and the discharge
    dies to the source-held floor; the reactor labels that state extinct rather than
    presenting it as a sustained steady state."""
    r, core, _ = _build(te_ev=1.0, x_ion=1e-9, power_w=0.0, source=6.6e4)
    _run(r, t_end=40.0)
    assert r.y[r.te_index] < 0.2 * EV_TO_K
    assert r.discharge_state() == 'extinct'
    r2, _, _ = _build(te_ev=1.0, x_ion=1e-9, power_w=0.5, source=6.6e4)
    _run(r2)
    assert r2.discharge_state() == 'sustained'


def test_latched_budget_closes_at_steady_state():
    """Clause 7: at steady state P_abs = sum of losses + d(U_e)/dt (solver derivative),
    and the latched wall electron energy flux equals the budget's wall electron term
    (clause 8, the same flux)."""
    r, core, _ = _build(te_ev=1.0, x_ion=1e-9, power_w=0.5, metastable=True, recombination=True)
    _run(r)
    b = r.energy_budget
    losses = b['Q_inelastic'] + b['Q_elastic'] + b['Q_wall_electron'] + b['Q_wall_ion'] + b['Q_flow']
    assert abs(b['P_abs'] - losses - b['dU_dt']) <= 1e-6 * b['P_abs']
    assert abs(b['dU_dt']) <= 1e-4 * b['P_abs']
    assert b['Q_wall_electron'] == pytest.approx(r.wall_electron_energy_flux, rel=1e-12)
