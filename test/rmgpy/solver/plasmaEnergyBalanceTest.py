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
from rmgpy.kinetics import Arrhenius, TwoTemperaturePlasma
from rmgpy.reaction import Reaction
from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
import rmgpy.solver.plasma as plasma_module
from rmgpy.solver.plasma import PLASMA_LOSCHMIDT, PlasmaReactor
from rmgpy.solver.termination import TerminationTime
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
A_POOL = 6.2e-16                          # m^3/s, Ars + Ars -> Ar + Ar+ + e (pooling, toy)
E_POOL_EV = 2.0 * E_EX_EV - E_IZ_EV       # 7.34 eV, credited to the ejected electron

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
                'Toy:3': (-E_EX_EV, 'eV'), 'Toy:4': (0.0, 'eV'), 'Toy:5': (0.076, 'eV'),
                'Toy:6': (-E_POOL_EV, 'eV')}
ARS_IGNORE = {'ignore': 'toy metastable: its elastic loss is not part of these checks'}


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
            'elastic_collisions': ({'Ar': dict(LL_AR_ELASTIC)} if elastic else
                                   {'Ar': {'ignore': 'elastic loss switched off for this check'}})}


def _toy(index, reactants, products, kinetics, entries):
    entries[index] = kinetics
    return LibraryReaction(reactants=reactants, products=products, reversible=False,
                           kinetics=kinetics, library='Toy')


def _build(te_ev=1.0, x_ion=1.0e-6, power_w=0.5, energy=True, metastable=False,
           recombination=False, elastic=True, source=None, radius=RADIUS, mixing=False,
           energies=None, termination=None, pooling=False, elastic_decl=None,
           cls=ToyLibraryReactor):
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
        if pooling:
            # heavy-heavy: gas-temperature kinetics (no electron reacts)
            rxns.append(_toy(6, [ars, ars], [ar, arp, electron],
                             Arrhenius(A=(A_POOL * constants.Na, 'm^3/(mol*s)'), n=0.0, Ea=(0.0, 'J/mol'),
                                       T0=(1.0, 'K')), entries))
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
        if metastable:
            kwargs['electron_energy_balance']['elastic_collisions']['Ars'] = dict(ARS_IGNORE)
        if elastic_decl is not None:
            kwargs['electron_energy_balance']['elastic_collisions'] = elastic_decl
    ToyLibraryReactor.toy_entries = entries
    reactor = cls((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf, (te_ev * EV_TO_K, 'K'),
                            n_sims=1, termination=termination or [], **kwargs)
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
    assert t['P_abs'] == pytest.approx(r.absorbed_power_reactor, rel=1e-14)
    assert lhs == pytest.approx(t['P_abs'] - losses, rel=1e-9, abs=1e-12 * t['P_abs'])


def _fd_jacobian(r, y, dydt, cj):
    """dres/dy + cj dres/d(dydt), both by central differences of the residual: the
    iteration matrix DASPK asks for, built without the engine's Jacobian code."""
    n = len(y)
    out = np.zeros((n, n))
    for k in range(n):
        h = 1e-5 * abs(y[k]) if y[k] != 0.0 else 1e-30
        yp, ym = y.copy(), y.copy()
        yp[k] += h
        ym[k] -= h
        dres_dy = (r.residual(0.0, yp, dydt)[0] - r.residual(0.0, ym, dydt)[0]) / (2 * h)
        dp, dm = dydt.copy(), dydt.copy()
        dp[k] += 1.0
        dm[k] -= 1.0
        dres_dyp = (r.residual(0.0, y, dp)[0] - r.residual(0.0, y, dm)[0]) / 2.0
        out[:, k] = dres_dy + cj * dres_dyp
    return out


def test_energy_mode_jacobian_matches_residual_differences():
    """Every column -- electron, ion, neutral, metastable, Te -- and the mass matrix at
    cj != 0 agree with a finite difference of the residual itself."""
    r, core, _ = _build(te_ev=1.0, metastable=True, recombination=True, pooling=True)
    y = np.array(r.y0, float)
    dydt = np.zeros_like(y)
    for cj in (0.0, 3.7e3):
        pd = np.array(r.jacobian(0.0, y, dydt, cj), float)
        fd = _fd_jacobian(r, y, dydt, cj)
        for k in range(len(y)):
            scale = np.max(np.abs(fd[:, k])) or 1.0
            assert np.allclose(pd[:, k], fd[:, k], rtol=1e-4, atol=1e-7 * scale), (cj, k)
    assert pd[_idx(r, core, 'Ar+'), r.te_index] > 0.0          # hotter electrons ionise faster


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
    decays -- 'extinct' -- until the external source holds it at the floor n_e =
    S / nu_wall(Tg), where it is 'source-supported'. Neither is presented as sustained."""
    r, core, _ = _build(te_ev=1.0, x_ion=1e-9, power_w=0.0, source=6.6e4)
    _run(r, t_end=0.05)
    assert r.discharge_state() == 'extinct'
    for t in np.logspace(np.log10(0.06), np.log10(400.0), 60):
        r.advance(t)
    assert r.y[r.te_index] < 0.2 * EV_TO_K
    assert r.discharge_state() == 'source-supported'
    n_ar = P_NOMINAL / (constants.R / constants.Na * TGAS)
    assert _n_e(r) == pytest.approx(6.6e4 / _nu_wall_hand(r.y[r.te_index], n_ar), rel=1e-2)
    r2, _, _ = _build(te_ev=1.0, x_ion=1e-9, power_w=0.5, source=6.6e4)
    _run(r2)
    assert r2.discharge_state() == 'self-sustained'


def _simulate(r, core, rxns):
    """The production entry, as the model builder drives it."""
    return r.simulate(core, rxns, [], [], [], [],
                      model_settings=ModelSettings(tol_keep_in_edge=0, tol_move_to_core=1e5,
                                                   tol_interrupt_simulation=1e8),
                      simulator_settings=SimulatorSettings())


def test_zero_power_run_terminates_as_extinct_on_simulate():
    """Clause 5 (PM ruling): extinction is a TERMINAL state, judged over physical time.
    Once the discharge is extinct and stays so -- Te within 1 K of Tg and nu_iz/nu_loss
    < 1e-6 -- for the declared multiple of the slower of 1/nu_loss and the electron
    energy-relaxation time, simulate() stops and records it."""
    r, core, rxns = _build(te_ev=1.0, x_ion=1e-9, power_w=0.0, source=6.6e4,
                           termination=[TerminationTime((40.0, 's'))])
    terminated, _res, _inv, _ss, _sr, t_final, _conv = _simulate(r, core, rxns)
    assert terminated
    assert r.terminal_state() == 'extinct'
    rec = r.energy_terminal
    assert rec['termination'] == 'extinct'
    assert rec['t'] == t_final < 40.0
    assert abs(rec['Te'] - TGAS) <= 1.0
    assert rec['nu_ionisation'] < 1e-6 * rec['nu_loss']
    assert rec['duration'] >= rec['persistence_time'] > 0.0
    assert rec['n_e'] == pytest.approx(_n_e(r), rel=1e-12)


def test_sustained_run_never_triggers_extinction():
    r, core, rxns = _build(te_ev=1.0, x_ion=1e-9, power_w=0.5, source=6.6e4,
                           termination=[TerminationTime((40.0, 's'))])
    terminated, _res, _inv, _ss, _sr, t_final, _conv = _simulate(r, core, rxns)
    assert terminated and t_final >= 40.0
    assert r.terminal_state() is None and r.energy_terminal is None
    assert r.discharge_state() == 'self-sustained'


def test_cold_start_powered_run_never_terminates_during_its_heating_transient():
    """Round 3 fix 2: a powered run started with Te at the gas temperature looks extinct
    until it heats. It has never been self-sustained and P_abs > 0, so it may not stop."""
    r, core, rxns = _build(te_ev=(TGAS + 0.5) / EV_TO_K, x_ion=1e-9, power_w=0.5, source=6.6e4,
                           termination=[TerminationTime((40.0, 's'))])
    terminated, _res, _inv, _ss, _sr, t_final, _conv = _simulate(r, core, rxns)
    assert r.energy_terminal is None and t_final >= 40.0
    assert r.discharge_state() == 'self-sustained'


def _extinct_toy():
    r, core, _ = _build(te_ev=1.0, x_ion=1e-9, power_w=0.0, source=6.6e4)
    _run(r, t_end=0.05)
    assert r.discharge_state() == 'extinct' and abs(r.y[r.te_index] - TGAS) <= 1.0
    return r, core, np.array(r.y, float)


def test_extinction_persistence_time_is_the_slower_of_loss_and_energy_relaxation():
    """The persistence time is the declared multiple of max(1/nu_loss, tau_E), with
    tau_E = 1 / (3 sum_p (m_e/M_p) K_p(Te) n_p) the elastic energy-relaxation time."""
    r, core, y = _extinct_toy()
    b = r.energy_budget
    V = r.compute_volume(y)
    te_ev = y[r.te_index] / EV_TO_K
    lt = np.log(te_ev)
    k = 2.336e-14 * te_ev ** 1.609 * np.exp(0.0618 * lt * lt - 0.1171 * lt ** 3)
    m_ar = core[[s.label for s in core].index('Ar')].molecular_weight.value_si
    nu_e = 3.0 * constants.m_e / m_ar * k * y[_idx(r, core, 'Ar')] * constants.Na / V
    hand = plasma_module.PLASMA_EXTINCT_PERSIST_MULTIPLE * max(1.0 / b['nu_loss'], 1.0 / nu_e)
    assert r.extinction_persistence_time(y) == pytest.approx(hand, rel=CONST_RTOL)


def test_extinction_counts_physical_time_not_calls():
    """Repeated calls at one time never count; the condition must hold over increasing
    time for the persistence time, and one miss restarts the clock."""
    r, core, y = _extinct_toy()
    tau = r.extinction_persistence_time(y)
    r.energy_terminal = None
    r.energy_extinct_since = float('nan')
    for _ in range(1000):
        r._update_terminal_state(y, 1.0)
    assert r.terminal_state() is None
    r._update_terminal_state(y, 1.0 + 0.5 * tau)
    r._update_terminal_state(y, 1.0 + 0.4 * tau)          # time going backwards: ignored
    assert r.terminal_state() is None
    miss = y.copy()
    miss[r.te_index] = TGAS + 2.0
    r._update_terminal_state(miss, 1.0 + 0.6 * tau)        # restarts the clock
    r._update_terminal_state(y, 1.0 + 0.7 * tau)
    r._update_terminal_state(y, 1.0 + 1.69 * tau)
    assert r.terminal_state() is None
    r._update_terminal_state(y, 1.0 + 1.71 * tau)
    assert r.terminal_state() == 'extinct'
    assert r.energy_terminal['duration'] == pytest.approx(1.01 * tau, rel=1e-9)


def test_powered_run_may_terminate_only_after_it_has_been_self_sustained():
    r0, core, y = _extinct_toy()
    r, _, _ = _build(te_ev=1.0, x_ion=1e-9, power_w=0.5, source=6.6e4)
    r.energy_budget = dict(r0.energy_budget)
    r.energy_extinct_since = float('nan')
    tau = r0.extinction_persistence_time(y)
    for i in range(50):
        r._update_terminal_state(y, 1.0 + i * tau)
    assert r.terminal_state() is None                       # never sustained, P_abs > 0
    r.energy_was_self_sustained = True
    for i in range(50, 60):
        r._update_terminal_state(y, 1.0 + i * tau)
    assert r.terminal_state() == 'extinct'


def test_energy_off_reactor_has_no_terminal_state():
    """The criterion belongs to the energy balance: with Te prescribed the reactor never
    evaluates it, so the energy-off path is unchanged."""
    r, core, rxns = _build(te_ev=1.0, x_ion=1e-9, energy=False, source=6.6e4,
                           termination=[TerminationTime((1e-3, 's'))])
    _simulate(r, core, rxns)
    assert r.terminal_state() is None and r.energy_terminal is None


# ------------------------------------------------------------------ discharge states

TOL = getattr(plasma_module, 'PLASMA_SELF_SUSTAINED_RTOL', 1e-3)


@pytest.mark.parametrize('nu_iz, nu_src, expect', [
    (1.0 - TOL, 0.0, 'self-sustained'),                          # at the cutoff
    ((1.0 - TOL) * (1.0 - 1e-9), 0.0, 'extinct'),                # just below, no source
    (0.6, 0.4, 'source-supported'),                              # the source makes up the deficit
    (0.0, 1.0, 'source-supported'),                              # a source-held floor
    (0.6, (1.0 - TOL) * (1.0 - 1e-9) - 0.6, 'extinct'),          # the source falls just short
    (0.6, (1.0 - TOL) - 0.6, 'source-supported'),                # ... or just makes it
    (0.5, 0.0, 'extinct'),                                       # the old >= 1/2 cutoff
    (1.3, 0.0, 'self-sustained'),                                # growing
])
def test_three_discharge_states_at_their_boundaries(nu_iz, nu_src, expect):
    assert plasma_module.classify_discharge(nu_iz, nu_src, 1.0) == expect


def test_latched_state_reads_the_external_source():
    """At P_abs = 0 the discharge decays (extinct); the latched budget carries the source
    frequency the classification used."""
    r, core, y = _extinct_toy()
    b = r.energy_budget
    ie = r.electron_index
    V = r.compute_volume(y)
    assert b['nu_source'] == pytest.approx(6.6e4 * V / constants.Na / y[ie], rel=1e-12)
    assert b['discharge_state'] == plasma_module.classify_discharge(b['nu_ionisation'], b['nu_source'], b['nu_loss'])


# ------------------------------------------------------------------ round 3 boundary fixes

def _metastable_noise():
    r, core, _ = _build(te_ev=1.0, metastable=True)
    y = np.array(r.y0, float)
    return r, core, y, _idx(r, core, 'Ars')


def test_energy_mode_clamps_a_neutral_negative_within_its_own_atol():
    """At an ACCEPTED step in energy mode a decayed neutral's noise within its own atol is
    clamped to exactly zero."""
    r, core, y, i = _metastable_noise()
    y[i] = -0.5 * r.atol_array[i]
    r.check_wall_support(y)
    assert y[i] == 0.0


def test_neutral_noise_tolerance_is_the_species_own_atol():
    r, core, y, i = _metastable_noise()
    r.atol_array[i] = 1.0e-10
    y[i] = -1.0e-12
    r.check_wall_support(y)
    assert y[i] == 0.0
    r, core, y, i = _metastable_noise()
    r.atol_array[i] = 1.0e-20
    y[i] = -1.0e-18
    with pytest.raises(PlasmaStateError, match='negative'):
        r.check_wall_support(y)
    assert y[i] == -1.0e-18


def test_an_initial_negative_is_refused_not_repaired():
    """The clamp repairs integrator noise at accepted steps only; the initial composition is
    the modeller's input and is never repaired."""
    r, core, y, i = _metastable_noise()
    y[i] = -0.5 * r.atol_array[i]
    with pytest.raises(PlasmaStateError, match='negative'):
        r.check_wall_support(y, accepted=False)
    assert y[i] == -0.5 * r.atol_array[i]


def test_energy_mode_refuses_zero_initial_electrons():
    """Round 3 fix 4: the Te row is an energy per electron. With no electrons it is
    undefined, so energy mode refuses N_e <= 0 at initialisation even when a source would
    ignite the gas."""
    with pytest.raises(PlasmaStateError, match='electron temperature is undefined'):
        _build(x_ion=0.0, source=6.6e4)


def test_nonpositive_electrons_in_the_residual_are_an_unreachable_state():
    r, core, _ = _build()
    y = np.array(r.y0, float)
    y[r.electron_index] = 0.0
    with pytest.raises(PlasmaStateError, match='unreachable'):
        r.residual(0.0, y, np.zeros_like(y))


class _SpyReactor(ToyLibraryReactor):
    """Records the smallest mole amount of every residual evaluation."""
    seen = []

    def residual(self, t, y, dydt, senpar=np.zeros(1, float)):
        _SpyReactor.seen.append(float(np.min(y[:self.num_core_species])))
        return ToyLibraryReactor.residual(self, t, y, dydt, senpar)


def test_jacobian_never_steps_a_mole_amount_below_zero():
    """A central difference at a near-zero amount would evaluate the residual at a negative
    one (a negative electron count for the Te row). One-sided there."""
    r, core, _ = _build(te_ev=1.0, cls=_SpyReactor)
    y = np.array(r.y0, float)
    y[r.electron_index] = 1.0e-40
    _SpyReactor.seen = []
    r.jacobian(0.0, y, np.zeros_like(y), 0.0)
    assert _SpyReactor.seen and min(_SpyReactor.seen) >= 0.0


def test_absorbed_power_is_a_constant_total_on_the_reactor_inventory():
    """Round 3 fix 5: P_abs is a total power. The reactor is a fixed-inventory image of the
    chamber, scaled once at t0 by V_ref = N_heavy(t0) R Tg / P, so its heating is the
    constant P_abs V_ref / V_chamber -- not (P_abs/V_chamber) V(state), which moves with
    the electron contribution to the constant-pressure volume."""
    r, core, _ = _build()
    y = np.array(r.y0, float)
    n_heavy = sum(y[j] for j in range(r.num_core_species) if j != r.electron_index)
    v_ref = n_heavy * constants.R * TGAS / P_NOMINAL
    expect = 0.5 * v_ref / CHAMBER_VOLUME
    assert r.absorbed_power_reactor == pytest.approx(expect, rel=1e-12)
    r.residual(0.0, y, np.zeros_like(y))
    p1 = r.electron_energy_terms['P_abs']
    y2 = y.copy()
    y2[r.electron_index] *= 1.0e3                           # V(state) moves
    assert r.compute_volume(y2) != r.compute_volume(y)
    r.residual(0.0, y2, np.zeros_like(y2))
    assert r.electron_energy_terms['P_abs'] == p1 == pytest.approx(expect, rel=1e-12)


def test_every_core_neutral_needs_an_elastic_rate_or_an_ignore():
    """Round 3 fix 6: an undeclared neutral is refused, not logged."""
    with pytest.raises(PlasmaStateError, match="'Ars'"):
        _build(metastable=True, elastic_decl={'Ar': dict(LL_AR_ELASTIC)})


@pytest.mark.parametrize('bad', [{'ignore': ''}, {'ignore': '   '}, {'ignore': 3},
                                 {'ignore': 'x', 'A': (1e-14, 'm^3/s')}])
def test_an_elastic_ignore_needs_a_reason_and_nothing_else(bad):
    with pytest.raises(PlasmaStateError, match='ignore'):
        _build(metastable=True, elastic_decl={'Ar': dict(LL_AR_ELASTIC), 'Ars': bad})


def test_an_elastic_ignore_is_recorded_and_contributes_nothing():
    r, core, _ = _build(metastable=True)
    assert r.energy_elastic_ignored == {'Ars': ARS_IGNORE['ignore']}
    assert r.energy_elastic_labels == ['Ar']


@pytest.mark.parametrize('te_frac, ok', [(0.49, False), (0.51, True)])
def test_initial_te_below_half_the_gas_temperature_is_refused(te_frac, ok):
    """Round 3 fix 7: an initial Te below the evaluation floor (Tg/2) would silently be
    evaluated at the floor; refuse it instead."""
    te_ev = te_frac * TGAS / EV_TO_K
    if ok:
        _build(te_ev=te_ev)
    else:
        with pytest.raises(PlasmaStateError, match='initial electron temperature'):
            _build(te_ev=te_ev)


def test_pooling_credits_its_declared_energy_to_the_electron():
    """Entry 90's shape: Ars + Ars -> Ar + Ar+ + e- ejects an electron with the excess
    2 E_m - E_iz, declared as a NEGATIVE loss (a gain). No electron is consumed, so no
    3/2 kTe is charged."""
    te_k = 1.0 * EV_TO_K
    r, core, rxns = _build(te_ev=1.0, metastable=True, pooling=True)
    y = np.array(r.y0, float)
    r.residual(0.0, y, np.zeros_like(y))
    j = r.energy_reaction_keys.index('Toy:6')
    assert r.energy_threshold[j] / EV_J_PER_MOL == pytest.approx(-E_POOL_EV, rel=CONST_RTOL)
    assert r.energy_electrons_consumed[j] == 0 and r.energy_electron_net[j] == 1
    V = r.compute_volume(y)
    c_ars = y[_idx(r, core, 'Ars')] / V
    hand = V * A_POOL * constants.Na * c_ars * c_ars * (-E_POOL_EV * EV_J_PER_MOL)
    assert r.electron_energy_terms['Q_inelastic_by_reaction'][j] == pytest.approx(hand, rel=CONST_RTOL)


def _independent_closure(sabotage=None):
    """Run the toy to steady state and return (closure, p_abs, terms), every loss term
    recomputed HERE from the solver's rates, the declared toy energies, the hand elastic fit
    and the wall-loss scratch -- never by re-reading the engine's energy terms.
    ``sabotage(reactor)`` may damage the engine before the run."""
    r, core, rxns = _build(te_ev=1.0, x_ion=1e-9, power_w=0.5, metastable=True,
                           recombination=True, pooling=True)
    if sabotage is not None:
        sabotage(r)
    _run(r)
    y = np.array(r.y, float)
    r.residual(r.t, y, np.zeros_like(y))
    V = r.compute_volume(y)
    R = constants.R
    te = y[r.te_index]
    ie = r.electron_index
    q_inel = 0.0
    for rxn in rxns:
        j = r.reaction_index[rxn]
        key = 'Toy:{0}'.format([k for k, v in ToyLibraryReactor.toy_entries.items()
                                if v is rxn.kinetics][0])
        eps = TOY_ENERGIES[key][0] * EV_J_PER_MOL
        n_e_in = sum(1 for sp in rxn.reactants if sp.label == 'e-')
        n_e_out = sum(1 for sp in rxn.products if sp.label == 'e-')
        consumed = max(n_e_in - n_e_out, 0)
        q_inel += r.core_reaction_rates[j] * V * (eps + consumed * 1.5 * R * te)
    te_ev = te / EV_TO_K
    lt = np.log(te_ev)
    k_el = 2.336e-14 * te_ev ** 1.609 * np.exp(0.0618 * lt * lt - 0.1171 * lt ** 3)
    m_ar = core[[s.label for s in core].index('Ar')].molecular_weight.value_si
    q_el = 3 * constants.m_e / m_ar * k_el * y[_idx(r, core, 'Ar')] * constants.Na / V * y[ie] * R * (te - TGAS)
    iarp = _idx(r, core, 'Ar+')
    q_we = 2 * R * te * (-r.wall_loss_rates[ie])
    q_wi = _sheath_factor(core[iarp].molecular_weight.value_si) * R * te * (-r.wall_loss_rates[iarp])
    n_heavy0 = sum(r.y0[j] for j in range(r.num_core_species) if j != ie)
    p_abs = 0.5 * (n_heavy0 * R * TGAS / P_NOMINAL) / CHAMBER_VOLUME
    du_dt = 1.5 * R * (te * r.dydt[ie] + y[ie] * r.dydt[r.te_index])
    closure = p_abs - q_inel - q_el - q_we - q_wi - du_dt
    return closure, p_abs, dict(r=r, q_inel=q_inel, q_we=q_we, q_wi=q_wi, du_dt=du_dt)


# The hand Te in eV uses CODATA-2018, the engine an older vintage; they differ by 1.0e-6 in
# Te, and the elastic term (~99 % of P here) goes as ~Te^1.6, so the independent sum carries
# ~1.6e-6 of P from the constants alone. 5e-6 sits above that and far below every term it
# guards (the smallest, the wall electrons, is ~1e-3 of P).
CLOSURE_RTOL = 5e-6


def test_latched_budget_closes_at_steady_state():
    """Clause 7, able to fail: P - sum Q - dU/dt from INDEPENDENTLY recomputed terms must
    vanish at the solver's steady state. An engine that dropped or mis-signed a term would
    steady-state where ITS sum vanishes, and this one would not (next test)."""
    closure, p_abs, t = _independent_closure()
    assert min(t['q_inel'], t['q_we'], t['q_wi']) > 100 * CLOSURE_RTOL * p_abs
    assert abs(closure) <= CLOSURE_RTOL * p_abs
    assert abs(t['du_dt']) <= 1e-4 * p_abs
    r = t['r']
    assert r.energy_budget['Q_wall_electron'] == pytest.approx(r.wall_electron_energy_flux, rel=1e-12)


def test_the_budget_check_catches_a_dropped_loss_term():
    """The closure above is not an identity: drop the engine's ion sheath energy (0.25 % of
    P) and the independent closure misses by about that much."""
    def drop_ion_sheath(r):
        r.energy_ion_sheath_factor = np.zeros_like(r.energy_ion_sheath_factor)
    closure, p_abs, t = _independent_closure(drop_ion_sheath)
    assert abs(closure) > 100 * CLOSURE_RTOL * p_abs
