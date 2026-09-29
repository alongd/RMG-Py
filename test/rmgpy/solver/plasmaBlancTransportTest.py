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
Composition-weighted (Blanc's-law) multi-bath wall transport on :class:`PlasmaReactor`
(I-294).

With a reduced mobility declared per (ion, bath gas) pair, the ion's mobility in the
mixture is

    1/mu_i = sum_b x_b / mu_{i,b},     mu_{i,b} = K0_{i,b} * N0 / N,

x_b the mole fraction of bath b among the declared baths (read from the reactor STATE),
N the total neutral density. A declared metastable's diffusivity combines the same way,
1/D_m = sum_b x_b / D_{m,b}. Bath identity is the declared neutral species
label; electronic states such as Ar and Ar* remain distinct unless explicitly
lumped.

Every expected value below is computed here, from the declared K0 / D*p numbers and the
state, never read back from another reactor output.

PLACEHOLDER DATA: the He-bath values (Ar+ in He, Ar2+ in He, Ar* in He) and the Ar2+ in
Ar value are NOT SOURCED. They are round test-only stand-ins chosen to differ strongly
from the Ar-bath values so the composition weighting is visible; they must never be
copied into a deck or a database.
"""

import copy
import json
import pickle
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import brentq

import rmgpy.constants as constants
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.exceptions import PlasmaStateError
from rmgpy.kinetics import Arrhenius, TwoTemperaturePlasma
from rmgpy.reaction import Reaction
from rmgpy.rmg.input import plasma_reactor
from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
from rmgpy.solver.base import ReactionSystem
from rmgpy.solver.plasma import PLASMA_LOSCHMIDT, PlasmaReactor
from rmgpy.solver.termination import TerminationSteadyState, TerminationTime
from rmgpy.species import Species
from rmgpy.thermo import ThermoData

EV_TO_K = 1.0 / 8.617333262e-5
EV_J_PER_MOL = 96485.33212
TORR_TO_PA = 101325.0 / 760.0
TGAS = 298.15
P_AR_TORR, P_HE_TORR = 5.0, 0.5               # 0.5 torr He in a 5 torr argon glow
RADIUS, LENGTH = 0.05, 0.30
LAMBDA = 1.0 / np.sqrt((2.405 / RADIUS) ** 2 + (np.pi / LENGTH) ** 2)
TE_EV = 1.0

MU = 'm^2/(V*s)'
# Ar+ in Ar: Ellis et al., ADNDT 17, 177 (1976) -- the value the argon decks carry.
K0_ARP_AR = 1.535e-4
# ---- PLACEHOLDER, NOT SOURCED (test-only) --------------------------------------------
K0_ARP_HE_PLACEHOLDER = 2.0e-3
K0_AR2P_AR_PLACEHOLDER = 1.9e-4
K0_AR2P_HE_PLACEHOLDER = 2.6e-3
DP_ARS_HE_PLACEHOLDER = 150.0                  # cm^2*torr/s
# ---------------------------------------------------------------------------------------
DP_ARS_AR = 47.0                               # cm^2*torr/s, the argon decks' value

FD_FRACTIONS = (1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
FD_TOLERANCE = 1e-5


def _thermo(h_ev):
    cp = 2.5 * constants.R
    return ThermoData(Tdata=([298, 400, 600, 800, 1000, 1500, 2000], 'K'),
                      Cpdata=([cp] * 7, 'J/(mol*K)'),
                      H298=(h_ev * EV_J_PER_MOL / 1000.0, 'kJ/mol'),
                      S298=(154.8, 'J/(mol*K)'))


def _species():
    """Fresh species (the solver indexes by identity, so never shared)."""
    e = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    he = Species(label='He').from_adjacency_list('1 He u0 p1 c0')
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    ar2p = Species(label='Ar2+').from_adjacency_list(
        'multiplicity 2\n1 Ar u0 p3 c+1 {2,S}\n2 Ar u1 p3 c0 {1,S}')
    ars = Species(label='Ars').from_adjacency_list('1 Ar u2 p3 c0')
    ar.thermo, he.thermo, ars.thermo = _thermo(0.0), _thermo(0.0), _thermo(11.55)
    arp.thermo, ar2p.thermo = _thermo(15.76), _thermo(14.5)
    return {'e-': e, 'Ar': ar, 'He': he, 'Ar+': arp, 'Ar2+': ar2p, 'Ars': ars}


def _arp_two_bath():
    return {'Ar+': {'perBath': {'Ar': (K0_ARP_AR, MU), 'He': (K0_ARP_HE_PLACEHOLDER, MU)}}}


def _two_ion_two_bath():
    return {'Ar+': {'perBath': {'Ar': (K0_ARP_AR, MU), 'He': (K0_ARP_HE_PLACEHOLDER, MU)}},
            'Ar2+': {'perBath': {'Ar': (K0_AR2P_AR_PLACEHOLDER, MU), 'He': (K0_AR2P_HE_PLACEHOLDER, MU)}}}


def _ars_two_bath():
    return {'Ars': {'product': 'Ar',
                    'diffusivity': {'Ar': (DP_ARS_AR, 'cm^2*torr/s'),
                                    'He': (DP_ARS_HE_PLACEHOLDER, 'cm^2*torr/s')}}}


def _reactor(labels=('e-', 'Ar', 'He', 'Ar+'), x=None, mobilities=None, pressure_torr=None,
             te_ev=TE_EV, cls=PlasmaReactor, rxns=None, spc=None, **kwargs):
    """A wall reactor over the named species, initialised. ``x`` maps label -> mole
    fraction; the default is 0.5 torr He in 5 torr Ar with a 1e-6 charge-neutral seed."""
    spc = spc or _species()
    if x is None:
        x = {'Ar': P_AR_TORR / (P_AR_TORR + P_HE_TORR), 'He': P_HE_TORR / (P_AR_TORR + P_HE_TORR)}
        x = {k: v for k, v in x.items() if k in labels}
        if 'Ar2+' in labels:
            x.update({'e-': 2e-6, 'Ar+': 1e-6, 'Ar2+': 1e-6})
        else:
            x.update({'e-': 1e-6, 'Ar+': 1e-6})
        if 'Ars' in labels:
            x['Ars'] = 1e-6
    imf = {spc[k]: v for k, v in x.items()}
    kwargs.setdefault('diffusion_length', (LAMBDA, 'm'))
    kwargs.setdefault('wall_neutralization_products',
                      {k: v for k, v in (('Ar+', 'Ar'), ('Ar2+', ('Ar', 2))) if k in labels})
    if 'Ars' in labels:
        # Electronic states are distinct collision partners unless this deck explicitly
        # adopts the named approximation that Ar* carries Ar's transport weight.
        kwargs.setdefault('wall_bath_lumping', {'Ars': 'Ar'})
    if mobilities is None and 'ion_reduced_mobility' not in kwargs:
        mobilities = _arp_two_bath() if 'Ar2+' not in labels else _two_ion_two_bath()
    if mobilities is not None:
        kwargs['ion_reduced_mobilities'] = mobilities
    p = (P_AR_TORR + P_HE_TORR) if pressure_torr is None else pressure_torr
    reactor = cls((TGAS, 'K'), (p * TORR_TO_PA, 'Pa'), imf, (te_ev * EV_TO_K, 'K'),
                  n_sims=1, termination=kwargs.pop('termination', []), **kwargs)
    core = [spc[k] for k in labels]
    reactor.initialize_model(core, rxns or [], [], [])
    return reactor, core, spc


def _i(reactor, spc, label):
    return reactor.species_index[spc[label]]


def _state(reactor, spc, amounts):
    y = np.zeros(reactor.num_core_species, float)
    for label, value in amounts.items():
        y[_i(reactor, spc, label)] = value
    return y


def _nu_hand(k0_eff, te_k, n_total):
    """nu = D_a/Lambda^2 with D_a = mu*kTe/e, mu = K0_eff*N0/N (independent re-derivation)."""
    mu = k0_eff * PLASMA_LOSCHMIDT / n_total
    return mu * (constants.R / constants.Na) * te_k / constants.e / (LAMBDA * LAMBDA)


def _blanc(pairs):
    """Hand Blanc combination: pairs = [(x_b, value_b)], returns 1/sum(x_b/value_b)."""
    return 1.0 / sum(xb / vb for xb, vb in pairs)


def _jacobian_scan(reactor, y):
    n = reactor.num_core_species
    zeros = np.zeros(n, float)
    reactor.jacobian(0.0, y, zeros, 0.0)
    analytic = np.array(reactor.jacobian_matrix, float)
    best = np.inf
    for frac in FD_FRACTIONS:
        fd = np.zeros((n, n), float)
        for k in range(n):
            h = frac * abs(y[k])
            yp, ym = y.copy(), y.copy()
            yp[k] += h
            ym[k] -= h
            fd[:, k] = (reactor.residual(0.0, yp, zeros)[0]
                        - reactor.residual(0.0, ym, zeros)[0]) / (2 * h)
        scale = np.maximum(np.abs(analytic), np.abs(fd))
        scale[scale == 0.0] = 1.0
        best = min(best, (np.abs(analytic - fd) / scale).max())
    return best


# ------------------------------------------------------------------ arithmetic


def test_blanc_mixture_mobility_matches_hand_values():
    """1/K_eff = 0.9/1.535e-4 + 0.1/2.0e-3 = 5913.192182... -> K_eff = 1.691133941e-4
    m^2/(V s) at x_He = 0.1, and nu_wall = K_eff*N0/N * kTe/e / Lambda^2 at that state --
    in the public helper, in the residual's wall flux, and in the latched per-ion flux."""
    r, core, spc = _reactor()
    y = _state(r, spc, {'Ar': 0.9, 'He': 0.1, 'Ar+': 1e-6, 'e-': 1e-6})
    i = _i(r, spc, 'Ar+')
    k_eff = r.compute_mixture_reduced_mobilities(y)[i]
    assert k_eff == pytest.approx(1.6911339411132838e-4, rel=1e-14)
    assert k_eff == pytest.approx(_blanc([(0.9, K0_ARP_AR), (0.1, K0_ARP_HE_PLACEHOLDER)]), rel=1e-14)

    V = r.compute_volume(y)
    n_total = 1.0 * constants.Na / V
    expected = _nu_hand(1.6911339411132838e-4, r.Te.value_si, n_total)
    nu = r.compute_ion_wall_frequencies(y, V)
    assert nu[i] == pytest.approx(expected, rel=1e-12)
    for j in range(r.num_core_species):
        if j != i:
            assert nu[j] == 0.0
    r.residual(0.0, y, np.zeros(r.num_core_species))
    assert -r.wall_loss_rates[i] / y[i] == pytest.approx(expected, rel=1e-12)
    # the pure-Ar and pure-He limits are the declared pair values themselves
    assert r.compute_mixture_reduced_mobilities(
        _state(r, spc, {'Ar': 1.0, 'Ar+': 1e-6, 'e-': 1e-6}))[i] == K0_ARP_AR
    assert r.compute_mixture_reduced_mobilities(
        _state(r, spc, {'He': 1.0, 'Ar+': 1e-6, 'e-': 1e-6}))[i] == K0_ARP_HE_PLACEHOLDER


def test_blanc_two_ions_each_combine_their_own_pairs():
    """Each cation combines ITS OWN (ion, bath) values over the same composition."""
    r, core, spc = _reactor(labels=('e-', 'Ar', 'He', 'Ar+', 'Ar2+'))
    y = _state(r, spc, {'Ar': 0.7, 'He': 0.3, 'Ar+': 1e-6, 'Ar2+': 2e-6, 'e-': 3e-6})
    k = r.compute_mixture_reduced_mobilities(y)
    assert k[_i(r, spc, 'Ar+')] == pytest.approx(
        _blanc([(0.7, K0_ARP_AR), (0.3, K0_ARP_HE_PLACEHOLDER)]), rel=1e-14)
    assert k[_i(r, spc, 'Ar2+')] == pytest.approx(
        _blanc([(0.7, K0_AR2P_AR_PLACEHOLDER), (0.3, K0_AR2P_HE_PLACEHOLDER)]), rel=1e-14)
    V = r.compute_volume(y)
    n_total = 1.0 * constants.Na / V
    r.residual(0.0, y, np.zeros(r.num_core_species))
    for label in ('Ar+', 'Ar2+'):
        j = _i(r, spc, label)
        assert -r.wall_loss_rates[j] / y[j] == pytest.approx(
            _nu_hand(k[j], r.Te.value_si, n_total), rel=1e-12)
    ie = r.electron_index
    assert r.wall_loss_rates[ie] == pytest.approx(
        r.wall_loss_rates[_i(r, spc, 'Ar+')] + r.wall_loss_rates[_i(r, spc, 'Ar2+')], rel=1e-14)


def test_blanc_metastable_diffusivity_matches_hand_values():
    """1/(D N)_eff = x_Ar/(D N)_Ar + x_He/(D N)_He with D*N = D*p/(k_B T_gas); Ars counts
    in the declared Ar bath. nu_m = (D N)_eff / (N Lambda^2)."""
    labels = ('e-', 'Ar', 'He', 'Ar+', 'Ars')
    r, core, spc = _reactor(labels=labels, wall_neutral_diffusion=_ars_two_bath())
    y = _state(r, spc, {'Ar': 0.8, 'Ars': 0.05, 'He': 0.15, 'Ar+': 1e-6, 'e-': 1e-6})
    V = r.compute_volume(y)
    kt = (constants.R / constants.Na) * TGAS
    cm2torr = 1e-4 * TORR_TO_PA
    dn_ar, dn_he = DP_ARS_AR * cm2torr / kt, DP_ARS_HE_PLACEHOLDER * cm2torr / kt
    x_ar, x_he = 0.85, 0.15
    dn_eff = _blanc([(x_ar, dn_ar), (x_he, dn_he)])
    n_total = 1.0 * constants.Na / V
    expected = dn_eff / (n_total * LAMBDA * LAMBDA)
    nu_m = r.compute_neutral_wall_frequencies(y, V)
    assert nu_m[_i(r, spc, 'Ars')] == pytest.approx(expected, rel=1e-12)
    r.residual(0.0, y, np.zeros(r.num_core_species))
    assert -r.wall_loss_rates[_i(r, spc, 'Ars')] / y[_i(r, spc, 'Ars')] == pytest.approx(expected, rel=1e-12)


def test_blanc_weights_track_the_reactor_state():
    """The composition is read from the state passed in, not frozen at t0: consuming He
    moves the mobility toward the argon value."""
    r, core, spc = _reactor()
    i = _i(r, spc, 'Ar+')
    for x_he in (0.5, 0.2, 0.01):
        y = _state(r, spc, {'Ar': 1.0 - x_he, 'He': x_he, 'Ar+': 1e-6, 'e-': 1e-6})
        assert r.compute_mixture_reduced_mobilities(y)[i] == pytest.approx(
            _blanc([(1.0 - x_he, K0_ARP_AR), (x_he, K0_ARP_HE_PLACEHOLDER)]), rel=1e-14)


# ------------------------------------------------------------ bit identity


def _single_bath_pair(per_bath):
    """Ar/Ars/Ar+/e- with a metastable wall loss and a pair source, in the flat map form
    or the per-bath form naming only the Ar bath."""
    labels = ('e-', 'Ar', 'Ars', 'Ar+')
    x = {'Ar': 1.0 - 3e-6, 'Ars': 1e-6, 'Ar+': 1e-6, 'e-': 1e-6}
    if per_bath:
        mob = {'Ar+': {'perBath': {'Ar': (K0_ARP_AR, MU)}}}
        wnd = {'Ars': {'product': 'Ar', 'diffusivity': {'Ar': (DP_ARS_AR, 'cm^2*torr/s')}}}
    else:
        mob = {'Ar+': (K0_ARP_AR, MU)}
        wnd = {'Ars': {'product': 'Ar', 'diffusivity': (DP_ARS_AR, 'cm^2*torr/s')}}
    return _reactor(labels=labels, x=x, mobilities=mob, pressure_torr=P_AR_TORR,
                    wall_neutral_diffusion=wnd,
                    wall_bath_lumping=({'Ars': 'Ar'} if per_bath else None),
                    ionisation_source=(6.6e4, 'm^-3/s'))


def test_single_bath_per_bath_form_is_bitwise_identical_to_the_flat_map():
    """A per-bath declaration naming one bath is the flat declaration, bit for bit: the
    residual, the Jacobian, and an integrated trajectory are identical, not merely close."""
    flat, _, _ = _single_bath_pair(False)
    blanc, _, _ = _single_bath_pair(True)
    n = flat.num_core_species
    y = flat.y0.copy()
    y[1] *= 0.999
    res_f = flat.residual(0.0, y, np.zeros(n))[0].copy()
    res_b = blanc.residual(0.0, y, np.zeros(n))[0].copy()
    assert np.array_equal(res_f, res_b)
    flat.jacobian(0.0, y, np.zeros(n), 0.0)
    blanc.jacobian(0.0, y, np.zeros(n), 0.0)
    assert np.array_equal(np.array(flat.jacobian_matrix), np.array(blanc.jacobian_matrix))
    for t in np.logspace(-8, -2, 13):
        flat.advance(t)
        blanc.advance(t)
        assert np.array_equal(flat.y, blanc.y), t
    assert np.array_equal(flat.wall_flux, blanc.wall_flux)


# ------------------------------------------------------------------ refusals


def test_undeclared_pair_is_refused_naming_the_pair_without_a_threshold():
    """No threshold declared: a neutral bath in the core with no (ion, bath) value is
    refused whatever its amount, naming the missing pair."""
    with pytest.raises(PlasmaStateError, match=r"\('Ar\+', 'He'\)"):
        _reactor(mobilities={'Ar+': {'perBath': {'Ar': (K0_ARP_AR, MU)}}})
    with pytest.raises(PlasmaStateError, match=r"\('Ar\+', 'He'\)"):
        _reactor(mobilities={'Ar+': {'perBath': {'Ar': (K0_ARP_AR, MU)}}},
                 x={'Ar': 1.0 - 2e-6, 'He': 0.0, 'Ar+': 1e-6, 'e-': 1e-6})


def test_undeclared_pair_above_the_threshold_is_refused_at_the_initial_state():
    with pytest.raises(PlasmaStateError) as exc:
        _reactor(mobilities={'Ar+': {'perBath': {'Ar': (K0_ARP_AR, MU)}}}, wall_bath_threshold=1e-3)
    assert "('Ar+', 'He')" in str(exc.value)
    assert 'threshold' in str(exc.value)


def test_undeclared_bath_below_an_explicit_threshold_is_ignored():
    """He at 1e-4 under a declared 1e-3 threshold is ignored in the weights (not in N):
    the mobility is the Ar value applied to the TOTAL neutral density."""
    x = {'Ar': 1.0 - 1e-4 - 2e-6, 'He': 1e-4, 'Ar+': 1e-6, 'e-': 1e-6}
    r, core, spc = _reactor(mobilities={'Ar+': {'perBath': {'Ar': (K0_ARP_AR, MU)}}}, x=x,
                            wall_bath_threshold=1e-3)
    y = r.y0.copy()
    i = _i(r, spc, 'Ar+')
    assert r.compute_mixture_reduced_mobilities(y)[i] == K0_ARP_AR
    V = r.compute_volume(y)
    n_total = (y[_i(r, spc, 'Ar')] + y[_i(r, spc, 'He')]) * constants.Na / V
    assert r.compute_ion_wall_frequencies(y, V)[i] == pytest.approx(
        _nu_hand(K0_ARP_AR, r.Te.value_si, n_total), rel=1e-12)


def test_threshold_is_documented_as_a_heuristic_composition_cutoff():
    docs = (Path(__file__).resolve().parents[3]
            / 'documentation/source/users/rmg/input.rst').read_text()
    for public_text in (plasma_reactor.__doc__, docs):
        normalized = ' '.join(public_text.split())
        assert 'heuristic composition cutoff' in normalized
        assert 'bounds the error' not in normalized
        assert 'limits the error' not in normalized

    normalized_docs = ' '.join(docs.split())
    assert 'Ar and Ar* are distinct baths' in normalized_docs
    assert 'scatter the ion identically' not in normalized_docs
    assert 'exact for an electronic ground state and its metastables' not in (
        normalized_docs)
    assert 'Ar with its metastables -- is exact' not in normalized_docs
    assert 'with its metastables) is exact' not in normalized_docs

    with pytest.raises(PlasmaStateError,
                       match='heuristic composition cutoff'):
        _reactor(wall_bath_threshold=0.02)


def test_threshold_label_requires_live_pair():
    """Only a live species in a live omitted bath changes its transport."""
    mobilities = {
        'Ar+': {'perBath': {'Ar': (K0_ARP_AR, MU)}},
        'Ar2+': {'perBath': {'Ar': (K0_AR2P_AR_PLACEHOLDER, MU),
                             'He': (K0_AR2P_HE_PLACEHOLDER, MU)}},
    }
    labels = ('e-', 'Ar', 'He', 'Ar+', 'Ar2+')
    for case, he_amount, arp_amount, expected in (
            ('neither', 0.0, 0.0, 'available'),
            ('bath', 1e-4, 0.0, 'available'),
            ('species', 0.0, 1e-6, 'available'),
            ('both', 1e-4, 1e-6,
             'available-blanc-threshold-approximation')):
        x = {'Ar': 1.0 - he_amount - 2e-6 - arp_amount,
             'He': he_amount, 'Ar+': arp_amount, 'Ar2+': 1e-6,
             'e-': 1e-6 + arp_amount}
        r, _, _ = _reactor(labels=labels, x=x, mobilities=mobilities,
                            wall_bath_threshold=1e-3)
        assert r.wall_energy_availability['wall_flux'] == expected


def _public_wall_output(reactor):
    """The public wall attributes consumed by the deck-side JSON writer."""
    return {
        'wall_flux': (None if reactor.wall_flux is None
                      else reactor.wall_flux.tolist()),
        'wall_electron_energy_flux': reactor.wall_electron_energy_flux,
        'wall_neutralization_energy_flux': (
            reactor.wall_neutralization_energy_flux),
        'wall_ion_energy_flux': reactor.wall_ion_energy_flux,
        'wall_energy_availability': reactor.wall_energy_availability,
    }


def _energy_balance(*bath_labels):
    """Minimal floating-wall energy balance for public-output regressions."""
    return {
        'absorbed_power': (0.0, 'W'),
        'chamber_volume': (np.pi * RADIUS ** 2 * LENGTH, 'm^3'),
        'sheath': 'floating_wall',
        'electron_energies': {},
        'elastic_collisions': {
            label: {'ignore': 'wall-output regression'}
            for label in bath_labels
        },
    }


def test_energy_balance_preserves_missing_ion_pair_provenance():
    """The sheath latch retains an omitted live mobility-pair qualifier."""
    x = {'Ar': 1.0 - 1e-4 - 2e-6, 'He': 1e-4,
         'Ar+': 1e-6, 'e-': 1e-6}
    mobilities = {
        'Ar+': {'perBath': {'Ar': (K0_ARP_AR, MU)}},
    }
    r, _, _ = _reactor(
        x=x, mobilities=mobilities, wall_bath_threshold=1e-3,
        electron_energy_balance=_energy_balance('Ar', 'He'))
    r.advance(1e-8)

    output = _public_wall_output(r)
    approximation = 'available-blanc-threshold-approximation'
    assert output['wall_energy_availability'] == {
        'wall_flux': approximation,
        'wall_electron_energy_flux': approximation,
        'wall_neutralization_energy_flux': approximation,
        'wall_ion_energy_flux': (
            'available-floating-wall-sheath-model-'
            'blanc-threshold-approximation'),
    }
    assert np.isfinite(output['wall_ion_energy_flux'])
    json.dumps(output, allow_nan=False)


def test_energy_balance_preserves_single_bath_provenance():
    """The sheath latch retains the explicit scalar-transport qualifier."""
    r, _, _ = _reactor(
        mobilities=None, ion_reduced_mobility=(K0_ARP_AR, MU),
        wall_single_bath_approximation=True,
        electron_energy_balance=_energy_balance('Ar', 'He'))
    r.advance(1e-8)

    output = _public_wall_output(r)
    approximation = 'available-single-bath-approximation'
    assert output['wall_energy_availability'] == {
        'wall_flux': approximation,
        'wall_electron_energy_flux': approximation,
        'wall_neutralization_energy_flux': approximation,
        'wall_ion_energy_flux': (
            'available-floating-wall-sheath-model-'
            'single-bath-approximation'),
    }
    assert np.isfinite(output['wall_ion_energy_flux'])
    json.dumps(output, allow_nan=False)


def test_energy_balance_public_output_rejects_nonfinite_ion_term():
    """A broken sheath term cannot be advertised or serialized as available."""
    r, _, spc = _reactor(
        labels=('e-', 'Ar', 'Ar+'), mobilities=None,
        ion_reduced_mobility=(K0_ARP_AR, MU),
        electron_energy_balance=_energy_balance('Ar'))
    r.energy_ion_sheath_factor[_i(r, spc, 'Ar+')] = np.inf
    r._latch_energy_budget(r.y0, 0.0)

    output = _public_wall_output(r)
    assert output['wall_ion_energy_flux'] is None
    assert output['wall_energy_availability'][
        'wall_ion_energy_flux'] == 'unavailable'
    json.dumps(output, allow_nan=False)


def test_real_solve_public_output_marks_missing_ion_pair_fields():
    """A deck writer serializes affected ion transport after a real solve."""
    x = {'Ar': 1.0 - 1e-4 - 2e-6, 'He': 1e-4,
         'Ar+': 1e-6, 'e-': 1e-6}
    mobilities = {
        'Ar+': {'perBath': {'Ar': (K0_ARP_AR, MU)}},
    }
    r, _, _ = _reactor(x=x, mobilities=mobilities,
                        wall_bath_threshold=1e-3)
    r.advance(1e-8)

    output = _public_wall_output(r)
    availability = output['wall_energy_availability']
    expected = 'available-blanc-threshold-approximation'
    assert availability == {
        'wall_flux': expected,
        'wall_electron_energy_flux': expected,
        'wall_neutralization_energy_flux': expected,
        'wall_ion_energy_flux': 'declared-absent',
    }
    assert 'wall_neutralization_energy_flux' in output
    assert output['wall_neutralization_energy_flux'] is not None
    assert output['wall_ion_energy_flux'] is None
    json.dumps(output, allow_nan=False)


def test_real_solve_public_output_marks_missing_neutral_pair_fields():
    """An omitted neutral pair affects only neutral flux after a real solve."""
    x = {'Ar': 1.0 - 1e-4 - 3e-6, 'He': 1e-4, 'Ars': 1e-6,
         'Ar+': 1e-6, 'e-': 1e-6}
    diffusion = {
        'Ars': {'product': 'Ar',
                'diffusivity': {'Ar': (DP_ARS_AR, 'cm^2*torr/s')}},
    }
    r, _, _ = _reactor(
        labels=('e-', 'Ar', 'He', 'Ar+', 'Ars'), x=x,
        mobilities=_arp_two_bath(), wall_neutral_diffusion=diffusion,
        wall_bath_threshold=1e-3)
    r.advance(1e-8)

    output = _public_wall_output(r)
    availability = output['wall_energy_availability']
    assert availability == {
        'wall_flux': 'available-blanc-threshold-approximation',
        'wall_electron_energy_flux': 'available',
        'wall_neutralization_energy_flux': 'available',
        'wall_ion_energy_flux': 'declared-absent',
    }
    assert 'wall_neutralization_energy_flux' in output
    assert output['wall_neutralization_energy_flux'] is not None
    assert output['wall_ion_energy_flux'] is None
    json.dumps(output, allow_nan=False)


def test_undeclared_bath_crossing_the_threshold_at_an_accepted_state_is_refused():
    x = {'Ar': 1.0 - 1e-4 - 2e-6, 'He': 1e-4, 'Ar+': 1e-6, 'e-': 1e-6}
    r, core, spc = _reactor(mobilities={'Ar+': {'perBath': {'Ar': (K0_ARP_AR, MU)}}}, x=x,
                            wall_bath_threshold=1e-3)
    y = r.y0.copy()
    y[_i(r, spc, 'He')] = 0.01 * y[_i(r, spc, 'Ar')]
    with pytest.raises(PlasmaStateError, match=r"\('Ar\+', 'He'\)"):
        r.check_wall_support(y)


def _transient_bath_reactor():
    """Make a real-solver Ar -> Ars -> He transient undeclared bath."""
    spc = _species()
    k1 = Arrhenius(A=(10.0, 's^-1'), n=0.0, Ea=(0.0, 'kJ/mol'))
    k2 = Arrhenius(A=(100.0, 's^-1'), n=0.0, Ea=(0.0, 'kJ/mol'))
    reactions = [Reaction(reactants=[spc['Ar']], products=[spc['Ars']],
                          kinetics=k1, reversible=False),
                 Reaction(reactants=[spc['Ars']], products=[spc['He']],
                          kinetics=k2, reversible=False)]
    labels = ('e-', 'Ar', 'Ars', 'He', 'Ar+')
    x = {'Ar': 1.0 - 2e-6 - 2e-8, 'Ars': 1e-8, 'He': 1e-8,
         'Ar+': 1e-6, 'e-': 1e-6}
    mobilities = {'Ar+': {'perBath': {'Ar': (K0_ARP_AR, MU),
                                      'He': (K0_ARP_HE_PLACEHOLDER, MU)}}}
    return _reactor(labels=labels, x=x, spc=spc, rxns=reactions,
                    mobilities=mobilities, wall_bath_lumping={},
                    wall_bath_threshold=0.01, diffusion_length=(1000.0, 'm'))


def test_advance_refuses_a_threshold_crossing_that_recedes_before_tout():
    """The real solver must expose transient Ars even though a single legacy
    ``advance`` lands after the bath recedes below the threshold."""
    control, _, spc = _transient_bath_reactor()
    ReactionSystem.advance(control, 1e-8)
    ReactionSystem.advance(control, 2.0)
    neutral_total = sum(control.y[_i(control, spc, label)]
                        for label in ('Ar', 'Ars', 'He'))
    assert control.y[_i(control, spc, 'Ars')] / neutral_total < 0.01

    reactor, _, _ = _transient_bath_reactor()
    reactor.advance(1e-8)
    with pytest.raises(PlasmaStateError, match=r"\('Ar\+', 'Ars'\)"):
        reactor.advance(2.0)


def test_thresholded_advance_lands_on_legacy_endpoint_when_never_approached():
    """An inactive threshold changes observation, not the solver endpoint."""
    legacy, _, _ = _reactor()
    thresholded, _, _ = _reactor(wall_bath_threshold=1e-3)
    legacy.advance(1e-4)
    thresholded.advance(1e-4)
    assert legacy.t == thresholded.t == 1e-4
    assert np.array_equal(legacy.y, thresholded.y)


class _NoProgressReactor(PlasmaReactor):
    def step(self, tout):
        return 0


def test_thresholded_advance_refuses_a_solver_step_that_makes_no_progress():
    reactor, _, _ = _reactor(cls=_NoProgressReactor, wall_bath_threshold=1e-3)
    with pytest.raises(PlasmaStateError, match='no progress'):
        reactor.advance(1e-4)


def test_missing_metastable_pair_is_refused_naming_the_pair():
    wnd = {'Ars': {'product': 'Ar', 'diffusivity': {'Ar': (DP_ARS_AR, 'cm^2*torr/s')}}}
    with pytest.raises(PlasmaStateError, match=r"\('Ars', 'He'\)"):
        _reactor(labels=('e-', 'Ar', 'He', 'Ar+', 'Ars'), wall_neutral_diffusion=wnd)


@pytest.mark.parametrize('kwargs,message', [
    (dict(wall_single_bath_approximation=True), 'wall_single_bath_approximation'),
    (dict(mobilities={'Ar+': {'perBath': {'Ar': (K0_ARP_AR, MU)}}, 'Ar2+': (K0_AR2P_AR_PLACEHOLDER, MU)},
          labels=('e-', 'Ar', 'Ar+', 'Ar2+'), x={'Ar': 1.0, 'Ar+': 1e-6, 'Ar2+': 1e-6, 'e-': 2e-6}),
     'per-bath'),
    (dict(mobilities={'Ar+': (K0_ARP_AR, MU)}, labels=('e-', 'Ar', 'Ar+'),
          x={'Ar': 1.0, 'Ar+': 1e-6, 'e-': 1e-6}, wall_bath_threshold=1e-3), 'wall_bath_threshold'),
    (dict(wall_bath_threshold=1.5), 'wall_bath_threshold'),
    (dict(wall_bath_threshold=-0.1), 'wall_bath_threshold'),
    (dict(wall_bath_threshold=float('nan')), 'wall_bath_threshold'),
    (dict(mobilities={'Ar+': {'perBath': {'Ar': (K0_ARP_AR, 'm^2/s'), 'He': (K0_ARP_HE_PLACEHOLDER, MU)}}}),
     "'Ar'"),
    (dict(mobilities={'Ar+': {'perBath': {'Ar': (K0_ARP_AR, MU), 'He': (-1.0, MU)}}}), "'He'"),
    (dict(mobilities={'Ar+': {'perBath': {'Ar': (K0_ARP_AR, MU), 'Ar+': (K0_ARP_AR, MU),
                                         'He': (K0_ARP_HE_PLACEHOLDER, MU)}}}), 'not a neutral'),
    (dict(ion_reduced_mobility=(K0_ARP_AR, MU), mobilities=None, labels=('e-', 'Ar', 'Ars', 'Ar+'),
          x={'Ar': 1.0, 'Ars': 1e-6, 'Ar+': 1e-6, 'e-': 1e-6},
          wall_neutral_diffusion={'Ars': {'product': 'Ar',
                                          'diffusivity': {'Ar': (DP_ARS_AR, 'cm^2*torr/s')}}}),
     'per-bath'),
    (dict(labels=('e-', 'Ar', 'He', 'Ars', 'Ar+'),
          wall_neutral_diffusion={'Ars': {'product': 'Ar',
                                          'diffusivity': (DP_ARS_AR, 'cm^2*torr/s')}}),
     'per-bath'),
])
def test_blanc_configuration_errors_are_refused_by_name(kwargs, message):
    kwargs = dict(kwargs)
    with pytest.raises(PlasmaStateError, match=message):
        _reactor(**kwargs)


def test_direct_empty_per_bath_map_has_a_named_refusal():
    with pytest.raises(PlasmaStateError, match='perBath.*non-empty'):
        _reactor(mobilities={'Ar+': {'perBath': {}}})


def test_direct_unconvertible_threshold_has_a_named_refusal():
    with pytest.raises(PlasmaStateError, match='wall_bath_threshold'):
        _reactor(wall_bath_threshold=10**400)


@pytest.mark.parametrize('mapping,message', [
    ({'Xe': 'Ar'}, 'source.*core neutral'),
    ({'Ar+': 'Ar'}, 'source.*core neutral'),
    ({'He': 'Xe'}, 'target.*declared bath'),
    ({'He': 'Ar+'}, 'target.*declared bath'),
    ({'He': 'He'}, 'self'),
    ({'He': 'Ar', 'Ar': 'He'}, 'cycle'),
])
def test_wall_bath_lumping_validates_complete_label_graph(mapping, message):
    """Lumping is a one-hop relation between real neutral core species."""
    with pytest.raises(PlasmaStateError, match=message):
        _reactor(wall_bath_lumping=mapping)


def test_wall_bath_lumping_refuses_chains():
    labels = ('e-', 'Ar', 'He', 'Ars', 'Ar+')
    x = {'Ar': 0.8, 'He': 0.1, 'Ars': 0.1 - 2e-6, 'Ar+': 1e-6, 'e-': 1e-6}
    mobilities = {'Ar+': {'perBath': {'Ar': (K0_ARP_AR, MU),
                                      'He': (K0_ARP_HE_PLACEHOLDER, MU)}}}
    with pytest.raises(PlasmaStateError, match='chain'):
        _reactor(labels=labels, x=x, mobilities=mobilities,
                 wall_bath_lumping={'Ars': 'He', 'He': 'Ar'})


def test_wall_bath_lumping_is_refused_outside_blanc_mode():
    with pytest.raises(PlasmaStateError, match='only.*per-bath'):
        _reactor(mobilities={'Ar+': (K0_ARP_AR, MU)},
                 wall_single_bath_approximation=True,
                 wall_bath_lumping={'He': 'Ar'})


def test_per_bath_map_refuses_two_labels_that_lump_to_one_group():
    labels = ('e-', 'Ar', 'Ars', 'Ar+')
    x = {'Ar': 1.0 - 3e-6, 'Ars': 1e-6, 'Ar+': 1e-6, 'e-': 1e-6}
    mobilities = {'Ar+': {'perBath': {'Ar': (K0_ARP_AR, MU),
                                      'Ars': (1.7e-4, MU)}}}
    with pytest.raises(PlasmaStateError, match='both.*lump.*same'):
        _reactor(labels=labels, x=x, mobilities=mobilities,
                 wall_bath_lumping={'Ars': 'Ar'})


def test_per_bath_diffusivity_refuses_two_labels_that_lump_to_one_group():
    labels = ('e-', 'Ar', 'Ars', 'Ar+')
    x = {'Ar': 1.0 - 3e-6, 'Ars': 1e-6, 'Ar+': 1e-6, 'e-': 1e-6}
    mobilities = {'Ar+': {'perBath': {'Ar': (K0_ARP_AR, MU)}}}
    diffusion = {'Ars': {'product': 'Ar',
                         'diffusivity': {'Ar': (DP_ARS_AR, 'cm^2*torr/s'),
                                         'Ars': (50.0, 'cm^2*torr/s')}}}
    with pytest.raises(PlasmaStateError, match='both.*lump.*same'):
        _reactor(labels=labels, x=x, mobilities=mobilities,
                 wall_neutral_diffusion=diffusion,
                 wall_bath_lumping={'Ars': 'Ar'})


def test_lumped_group_reports_its_declared_target_label():
    mobilities = {'Ar+': {'perBath': {'He': (K0_ARP_HE_PLACEHOLDER, MU)}}}
    r, _, _ = _reactor(mobilities=mobilities, wall_bath_lumping={'Ar': 'He'})
    assert r.wall_bath_labels == ['He']


def test_blanc_source_comments_define_baths_by_species_label():
    source = (Path(__file__).resolve().parents[3]
              / 'rmgpy/solver/plasma.pyx').read_text()
    normalized = ' '.join(source.split())
    assert 'Each neutral species label defines a distinct bath' in normalized
    assert 'Blanc bath identity is label-based' in normalized
    for stale in ('A bath is a heavy skeleton',
                  'bath (skeleton) index',
                  'Moles per bath (skeleton) group',
                  'are one bath exactly',
                  'summing them is summing one bath',
                  'Keyed on the heavy skeleton, the identity'):
        assert stale not in normalized


@pytest.mark.parametrize('ambiguous_label', ['', 'X'])
def test_unlabelled_or_duplicate_label_neutrals_do_not_merge(ambiguous_label):
    spc = _species()
    first = Species(label=ambiguous_label).from_adjacency_list('1 He u0 p1 c0')
    second = Species(label=ambiguous_label).from_adjacency_list(
        '1 Ne u0 p4 c0')
    first.label = second.label = ambiguous_label
    imf = {spc['Ar']: 1.0 - 2e-6, first: 0.0, second: 0.0,
           spc['Ar+']: 1e-6, spc['e-']: 1e-6}
    reactor = PlasmaReactor(
        (TGAS, 'K'), ((P_AR_TORR + P_HE_TORR) * TORR_TO_PA, 'Pa'), imf,
        (TE_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
        diffusion_length=(LAMBDA, 'm'),
        ion_reduced_mobilities={
            'Ar+': {'perBath': {'Ar': (K0_ARP_AR, MU)}}},
        wall_neutralization_products={'Ar+': 'Ar'}, wall_bath_threshold=0.01)
    core = [spc['e-'], spc['Ar'], first, second, spc['Ar+']]
    reactor.initialize_model(core, [], [], [])
    assert reactor.wall_bath_group[2] != reactor.wall_bath_group[3]
    assert reactor.wall_bath_labels[reactor.wall_bath_group[2]] != \
        reactor.wall_bath_labels[reactor.wall_bath_group[3]]


def test_ar_and_metastable_are_distinct_baths_without_explicit_lumping():
    """Ar* is not silently folded into Ar: only a deck-declared lumping does that."""
    mob = {'Ar+': {'perBath': {'Ar': (K0_ARP_AR, MU), 'Ars': (1.7e-4, MU),
                               'He': (K0_ARP_HE_PLACEHOLDER, MU)}}}
    r, _, _ = _reactor(labels=('e-', 'Ar', 'Ars', 'He', 'Ar+'), mobilities=mob,
                       wall_bath_lumping={})
    assert set(r.wall_bath_labels) == {'Ar', 'Ars', 'He'}


def test_flat_decks_keep_the_single_bath_refusal():
    """The legacy flat form on a two-gas bath still refuses unless the approximation is
    opted into -- the per-bath form is the alternative, not a silent default."""
    with pytest.raises(PlasmaStateError, match='spans more than one gas'):
        _reactor(mobilities={'Ar+': (K0_ARP_AR, MU)})
    r, _, _ = _reactor(mobilities={'Ar+': (K0_ARP_AR, MU)}, wall_single_bath_approximation=True)
    assert r.wall_bath_is_mixture


# ------------------------------------------------------------------ Jacobian / copies


def test_blanc_wall_jacobian_matches_finite_difference():
    """Two ions, two baths, a two-bath metastable and a pair source at prescribed Te: the
    analytic Jacobian carries the composition derivative of every Blanc weight."""
    labels = ('e-', 'Ar', 'He', 'Ars', 'Ar+', 'Ar2+')
    r, core, spc = _reactor(labels=labels, wall_neutral_diffusion=_ars_two_bath(),
                            ionisation_source=(6.6e4, 'm^-3/s'), wall_recycling=0.5)
    y = _state(r, spc, {'Ar': 0.7, 'He': 0.25, 'Ars': 0.05, 'Ar+': 2e-6, 'Ar2+': 3e-6, 'e-': 5e-6})
    best = _jacobian_scan(r, y)
    assert best < FD_TOLERANCE, best


def test_blanc_jacobian_has_one_sided_zero_bath_derivative_and_single_bath_limit():
    """A declared He bath produced from zero has the same one-sided derivative that
    connects the exact Ar-only limit to the two-bath Blanc expression."""
    r, _, spc = _reactor()
    y = _state(r, spc, {'Ar': 1.0, 'He': 0.0, 'Ar+': 1e-6, 'e-': 1e-6})
    n = r.num_core_species
    zeros = np.zeros(n)
    r.jacobian(0.0, y, zeros, 0.0)
    analytic = np.array(r.jacobian_matrix, float)
    he = _i(r, spc, 'He')
    h = 1e-8
    forward = (r.residual(0.0, y + np.eye(1, n, he)[0] * h, zeros)[0]
               - r.residual(0.0, y, zeros)[0]) / h
    scale = np.maximum(np.maximum(np.abs(analytic[:, he]), np.abs(forward)), 1.0)
    assert np.max(np.abs(analytic[:, he] - forward) / scale) < 2e-5
    ar = _i(r, spc, 'Ar')
    assert r.compute_mixture_reduced_mobilities(y)[_i(r, spc, 'Ar+')] == K0_ARP_AR
    assert forward[_i(r, spc, 'Ar+')] != 0.0


def test_blanc_jacobian_matches_the_residual_at_a_negative_bath_trial_amount():
    """The derivative follows Blanc clipping for a negative Newton trial."""
    r, _, spc = _reactor()
    y = _state(r, spc, {'Ar': 1.1, 'He': -0.1, 'Ar+': 1e-6, 'e-': 1e-6})
    best = _jacobian_scan(r, y)
    assert best < FD_TOLERANCE, best


@pytest.mark.parametrize('copier', [copy.deepcopy, lambda r: pickle.loads(pickle.dumps(r))])
def test_blanc_declarations_survive_a_copy(copier):
    labels = ('e-', 'Ar', 'He', 'Ars', 'Ar+')
    r, core, spc = _reactor(labels=labels, wall_neutral_diffusion=_ars_two_bath(),
                            wall_bath_threshold=1e-6)
    clone = copier(r)
    assert clone.ion_reduced_mobilities == r.ion_reduced_mobilities
    assert clone.wall_neutral_diffusion == r.wall_neutral_diffusion
    assert clone.wall_bath_threshold == 1e-6
    # a copy carries copies of the species: initialise it on its own, in the same order
    by_label = {s.label: s for s in clone.initial_mole_fractions}
    clone.initialize_model([by_label[s.label] for s in core], [], [], [])
    y = _state(r, spc, {'Ar': 0.7, 'He': 0.25, 'Ars': 0.05, 'Ar+': 2e-6, 'e-': 2e-6})
    n = r.num_core_species
    assert np.array_equal(r.residual(0.0, y, np.zeros(n))[0], clone.residual(0.0, y, np.zeros(n))[0])


# ------------------------------------------------------ two-bath steady state


A_IZ, E_IZ_EV = 2.34e-14, 15.76      # toy ionisation law k = A exp(-E/kTe), m^3/s
LL_AR_ELASTIC = {'A': (2.336e-14, 'm^3/s'), 'n': 1.609, 'b': 0.0618, 'c': -0.1171}


class _ToyReactor(PlasmaReactor):
    toy_entries = {}

    def _declared_entry_kinetics(self, library, index):
        if library != 'Toy' or index not in self.toy_entries:
            raise PlasmaStateError('no entry {0}:{1}'.format(library, index))
        return self.toy_entries[index]


def test_two_bath_energy_balance_deck_reaches_steady_state_on_blanc_transport():
    """0.5 torr He in 5 torr Ar, Ar+ with per-bath mobilities, Te solved by the energy
    balance. The run reaches a steady state; there the latched per-ion wall frequency
    matches the hand value from the declared K0s and the state's composition to 1e-6, and
    Te satisfies the particle balance k_iz(Te) n_Ar = nu_wall(Te) under Blanc transport."""
    spc = _species()
    e, ar, arp = spc['e-'], spc['Ar'], spc['Ar+']
    kin = TwoTemperaturePlasma(A=(A_IZ, 'm^3/(molecule*s)'), n=0.0,
                               Ea_g=(E_IZ_EV, 'eV/molecule'), Ea_e=(E_IZ_EV, 'eV/molecule'),
                               T0=(1.0, 'K'))
    _ToyReactor.toy_entries = {1: kin}
    rxn = LibraryReaction(reactants=[e, ar], products=[arp, e, e], reversible=False,
                          kinetics=kin, library='Toy')
    energy = {'absorbed_power': (0.5, 'W'), 'chamber_volume': (np.pi * RADIUS ** 2 * LENGTH, 'm^3'),
              'sheath': 'floating_wall', 'electron_energies': {'Toy:1': (E_IZ_EV, 'eV')},
              'elastic_collisions': {
                  'Ar': dict(LL_AR_ELASTIC),
                  'He': {'ignore': 'PLACEHOLDER test: He elastic loss is not part of this '
                                   'transport check'}}}
    r, core, spc = _reactor(cls=_ToyReactor, rxns=[rxn], spc=spc,
                            electron_energy_balance=energy,
                            termination=[TerminationSteadyState(tolerance=1e-6),
                                         TerminationTime((10.0, 's'))])
    terminated = r.simulate(core, [rxn], [], [], [], [],
                            model_settings=ModelSettings(tol_keep_in_edge=0, tol_move_to_core=1e5,
                                                         tol_interrupt_simulation=1e8),
                            simulator_settings=SimulatorSettings())[0]
    assert terminated
    assert r.steady_state_reached
    assert r.wall_energy_availability['wall_flux'] == 'available'

    y = np.array(r.y[:r.num_core_species], float)
    te = float(r.y[r.te_index])
    V = r.compute_volume(np.array(r.y, float))
    i_ar, i_he, i_arp = _i(r, spc, 'Ar'), _i(r, spc, 'He'), _i(r, spc, 'Ar+')
    n_total = (y[i_ar] + y[i_he]) * constants.Na / V
    x_ar = y[i_ar] / (y[i_ar] + y[i_he])
    k_eff = _blanc([(x_ar, K0_ARP_AR), (1.0 - x_ar, K0_ARP_HE_PLACEHOLDER)])
    nu_hand = _nu_hand(k_eff, te, n_total)
    nu_latched = -r.wall_flux[i_arp] / y[i_arp]
    assert nu_latched == pytest.approx(nu_hand, rel=1e-6)
    assert x_ar == pytest.approx(P_AR_TORR / (P_AR_TORR + P_HE_TORR), rel=1e-6)

    # Te from the particle balance alone, with the reactor's own eV <-> K conversion.
    n_ar = y[i_ar] * constants.Na / V
    kte = lambda t: t * (constants.R / constants.Na) / constants.e       # eV

    def balance(t):
        return A_IZ * np.exp(-E_IZ_EV / kte(t)) * n_ar - _nu_hand(k_eff, t, n_total)
    te_star = brentq(balance, 0.3 * EV_TO_K, 5.0 * EV_TO_K, xtol=1e-9)
    assert te == pytest.approx(te_star, rel=1e-4)
