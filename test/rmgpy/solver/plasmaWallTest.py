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
Tests for the charged-particle wall boundary operator on :class:`PlasmaReactor`.

This is the pytest port of the standalone harness in
``docs/i246-ambipolar-wall-operator/verify_operator.py`` and
``verify_zero_wall_bitwise.py``, using the same argon model as
``docs/i246-ambipolar-wall-operator/argon_wall_model.py``. It exists so the
wall operator's behaviour is checked by CI, not only by a script someone has
to remember to run.

No assertion below states a target electron density. Everything is expressed
as a loss frequency, a charge, a ratio, or a conservation statement.
"""

import copy
import logging
import os
import pickle
import re

import numpy as np
import pytest

import rmgpy.constants as constants
import rmgpy.solver.plasma
from rmgpy import settings
from rmgpy.exceptions import PlasmaStateError
from rmgpy.kinetics import Arrhenius, VoronovEIArrhenius
from rmgpy.reaction import Reaction
from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
from rmgpy.solver.plasma import PLASMA_LOSCHMIDT, PlasmaReactor
from rmgpy.solver.termination import TerminationSteadyState, TerminationTime
from rmgpy.species import Species
from rmgpy.thermo import ThermoData

EV_TO_K = 1.0 / 8.617333262e-5           # K per eV
TORR_TO_PA = 101325.0 / 760.0

TGAS = 298.15                             # K
P_NOMINAL = 5.0 * TORR_TO_PA              # Pa
TE_NOMINAL_EV = 3.0

R_NOMINAL = 0.05                          # m, cylinder radius
L_NOMINAL = 0.30                          # m, cylinder length

# Ellis, Pai, McDaniel, Mason & Viehland, At. Data Nucl. Data Tables 17, 177 (1976),
# doi:10.1016/0092-640X(76)90001-2. The exact table entry behind 1.535 is unverified.
MU0_AR_IN_AR = 1.535e-4                   # m^2/(V s) at the Loschmidt density

# Resolve the Voronov coefficients from the configured database, the way the rest of
# the suite locates database files, rather than hard-coding one machine's checkout.
# ``database.directory`` comes from the worktree's rmgrc (see docs/rmgrc.md).
VORONOV_YAML = os.path.join(settings['database.directory'], 'kinetics', 'voronov.yaml')


def _diffusion_length(radius=R_NOMINAL, length=L_NOMINAL):
    """Lowest diffusion eigenmode of a finite cylinder: 1/Lambda^2 = (2.405/R)^2 + (pi/L)^2."""
    return 1.0 / np.sqrt((2.405 / radius) ** 2 + (np.pi / length) ** 2)


def _nu_wall_closed_form(te_ev, n_neutral, lam, mu0=MU0_AR_IN_AR):
    """Independent re-derivation of nu_wall = D_a/Lambda^2, D_a = mu_i*kTe/e,
    mu_i = mu0*N0/n_neutral. Kept separate from compute_nu_wall so the two can
    be compared rather than one checking itself."""
    mu_i = mu0 * PLASMA_LOSCHMIDT / n_neutral
    d_a = mu_i * (constants.R / constants.Na) * (te_ev * EV_TO_K) / constants.e
    return d_a / (lam * lam)


def _argon_species():
    """Fresh species objects (solver indexing is by identity, so never shared)."""
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    return electron, ar, arp


def _ionisation_reaction(electron, ar, arp):
    """e + Ar -> Ar+ + e + e, from the repository's own Voronov coefficients."""
    kinetics = VoronovEIArrhenius(Z=18, N=18, yaml_path_or_obj=VORONOV_YAML)
    return Reaction(reactants=[electron, ar], products=[arp, electron, electron],
                    reversible=False, kinetics=kinetics)


def _build_reactor(te_ev=TE_NOMINAL_EV, pressure=P_NOMINAL, tgas=TGAS,
                   x_ion=1.0e-6, wall=True, gamma=1.0, source=None,
                   quasineutral=False, with_chemistry=True, max_alpha=None,
                   lam=None, mu0=MU0_AR_IN_AR, termination=None, reactor_cls=None):
    """A fully initialised argon PlasmaReactor, with or without a wall.

    ``x_ion`` is the charge-neutral seed mole fraction of Ar+ (and of e-); it
    is an initial condition only."""
    electron, ar, arp = _argon_species()
    imf = {electron: x_ion, arp: x_ion, ar: 1.0 - 2.0 * x_ion}
    kwargs = {}
    if wall:
        kwargs['diffusion_length'] = (lam if lam is not None
                                      else _diffusion_length(), 'm')
        kwargs['ion_reduced_mobility'] = (mu0, 'm^2/(V*s)')
        kwargs['wall_recycling'] = gamma
        if source is not None:
            kwargs['ionisation_source'] = (source, 'm^-3/s')
        if max_alpha is not None:
            kwargs['max_ionisation_degree'] = max_alpha
    reactor = (reactor_cls or PlasmaReactor)(
        (tgas, 'K'), (pressure, 'Pa'), imf, (te_ev * EV_TO_K, 'K'),
        n_sims=1, termination=termination or [],
        quasineutral_electron=quasineutral, **kwargs)
    core_species = [electron, ar, arp]
    core_reactions = [_ionisation_reaction(electron, ar, arp)] if with_chemistry else []
    reactor.initialize_model(core_species, core_reactions, [], [])
    return reactor, core_species, core_reactions


def _state_at(reactor, alpha, n_total_mol=1.0):
    """Packed state at ionisation degree alpha = n_e/n_neutral, charge neutral."""
    y = np.zeros(reactor.num_core_species, float)
    ie, i_ar, i_arp = _indices(reactor)
    y[i_ar] = n_total_mol
    y[i_arp] = alpha * n_total_mol
    y[ie] = alpha * n_total_mol
    return y


def _indices(reactor):
    z = reactor.species_charges
    ie = reactor.electron_index
    i_ar = [j for j in range(len(z)) if z[j] == 0][0]
    i_arp = [j for j in range(len(z)) if z[j] == 1 and j != ie][0]
    return ie, i_ar, i_arp


# ---------------------------------------------------------------- item 1


def test_nu_wall_matches_closed_form_rederivation():
    """nu_wall = D_a/Lambda**2 matches an independent closed-form re-derivation
    to machine precision."""
    r, _, _ = _build_reactor(wall=True, with_chemistry=False)
    ie, i_ar, i_arp = _indices(r)
    y = _state_at(r, 1.0e-6)
    V = r.compute_volume(y)
    nu = r.compute_nu_wall(y, V)
    n_neutral = y[i_ar] * constants.Na / V
    lam = r.diffusion_length.value_si
    closed = _nu_wall_closed_form(TE_NOMINAL_EV, n_neutral, lam)
    assert abs(nu / closed - 1.0) < 1e-12


# ---------------------------------------------------------------- item 2


def test_nu_wall_geometry_scaling():
    """Halving Lambda**2 exactly doubles nu_wall; nu_wall*Lambda**2 is constant
    across a range of chamber radii (the diffusivity itself carries no geometry)."""
    r, _, _ = _build_reactor(wall=True, with_chemistry=False)
    y = _state_at(r, 1.0e-6)
    lam = r.diffusion_length.value_si
    nu = r.compute_nu_wall(y, r.compute_volume(y))

    lam_half_sq = lam / np.sqrt(2.0)      # halves Lambda^2
    r2, _, _ = _build_reactor(wall=True, with_chemistry=False, lam=lam_half_sq)
    nu2 = r2.compute_nu_wall(y, r2.compute_volume(y))
    assert abs(nu2 / nu - 2.0) < 1e-12

    ratios = []
    for radius in (0.005, 0.01, 0.05, 0.1, 0.2):
        lr = _diffusion_length(radius, L_NOMINAL)
        rr, _, _ = _build_reactor(wall=True, with_chemistry=False, lam=lr)
        ratios.append(rr.compute_nu_wall(y, rr.compute_volume(y)) * lr * lr)
    assert max(ratios) / min(ratios) - 1.0 < 1e-12


# ---------------------------------------------------------------- item 3


def test_wall_loss_is_first_order_in_electron_population():
    """The wall loss is first order in the ELECTRON population:
    -(wall flux)/N_e == nu_wall exactly. Normalising to the neutral population
    instead gives a number orders of magnitude away -- a mis-normalisation
    error made earlier in this project."""
    r, _, _ = _build_reactor(wall=True, with_chemistry=False)
    ie, i_ar, i_arp = _indices(r)
    y = _state_at(r, 1.0e-6)
    V = r.compute_volume(y)
    nu = r.compute_nu_wall(y, V)
    r.residual(0.0, y, np.zeros(r.num_core_species, float))
    wall = r.wall_loss_rates

    loss_per_electron = -wall[ie] / y[ie]
    assert abs(loss_per_electron / nu - 1.0) < 1e-14

    per_neutral = -wall[ie] / y[i_ar]
    assert abs(np.log10(per_neutral / nu)) > 3.0, (
        "normalising the electron wall flux to the neutral population instead "
        "of the electron population should be orders of magnitude off nu_wall, "
        "but was only {0!r}".format(per_neutral / nu))


# ---------------------------------------------------------------- item 4


def test_net_charge_flux_zero_on_neutral_nonzero_on_nonneutral():
    """Net charge flux to the wall is exactly 0.0 on a charge-neutral state,
    and non-zero (both directions) on a deliberately non-neutral one."""
    r, _, _ = _build_reactor(wall=True, with_chemistry=False)
    ie, i_ar, i_arp = _indices(r)
    y = _state_at(r, 1.0e-6)
    r.residual(0.0, y, np.zeros(r.num_core_species, float))
    net_wall_current = float(np.sum(r.species_charges * r.wall_loss_rates))
    assert net_wall_current == 0.0

    # wall_loss_rates is signed as a residual contribution (negative = loss), so
    # net_wall_current = sum(charge * wall) = -nu_wall * net_charge(y). Excess
    # electrons make net_charge(y) more negative, hence net_wall_current > 0.
    y_extra_e = y.copy()
    y_extra_e[ie] *= 1.5
    r.residual(0.0, y_extra_e, np.zeros(r.num_core_species, float))
    net_extra_e = float(np.sum(r.species_charges * r.wall_loss_rates))
    assert net_extra_e != 0.0
    assert net_extra_e > 0.0

    # deficient electrons (excess positive charge): net_charge(y) goes positive,
    # so net_wall_current < 0.
    y_short_e = y.copy()
    y_short_e[ie] *= 0.5
    r.residual(0.0, y_short_e, np.zeros(r.num_core_species, float))
    net_short_e = float(np.sum(r.species_charges * r.wall_loss_rates))
    assert net_short_e != 0.0
    assert net_short_e < 0.0


# ---------------------------------------------------------------- item 5


def test_net_charge_drift_zero_under_full_residual():
    """d(net charge)/dt from the full residual (chemistry + wall) is zero."""
    r, _, _ = _build_reactor(wall=True, with_chemistry=True)
    y = _state_at(r, 1.0e-6)
    dc, _ = r.residual(0.0, y, np.zeros(r.num_core_species, float))
    drift = float(np.sum(r.species_charges * dc))
    scale = float(np.sum(np.abs(r.species_charges * dc)))
    assert abs(drift) <= 1e-12 * max(scale, 1.0)


# ---------------------------------------------------------------- item 6


def test_heavy_atom_recycling_conservation():
    """Heavy atoms leave the gas at exactly (1-gamma)*nu_wall*N_ion, checked
    at gamma = 1.0, 0.5, 0.0."""
    r0, _, _ = _build_reactor(wall=True, with_chemistry=False)
    y = _state_at(r0, 1.0e-6)
    for gamma in (1.0, 0.5, 0.0):
        rg, _, _ = _build_reactor(wall=True, with_chemistry=False, gamma=gamma)
        ieg, i_arg, i_arpg = _indices(rg)
        dg, _ = rg.residual(0.0, y, np.zeros(rg.num_core_species, float))
        d_heavy = dg[i_arg] + dg[i_arpg]
        nug = rg.nu_wall
        expected = -(1.0 - gamma) * nug * y[i_arpg]
        assert abs(d_heavy - expected) <= 1e-12 * max(abs(expected), abs(nug * y[i_arpg]))


# ---------------------------------------------------------------- item 7


def test_wall_operator_residual_inert_and_jacobian_is_the_wall_linearization():
    """Two invariants on the wall operator, at an all-neutral state, each of which a
    deterministic corruption of the wall path breaks -- unlike the earlier version of
    this test, which compared two WALL-LESS reactors and so never evaluated the wall
    path at all (see docs/i246-ambipolar-wall-operator/evidence/bitwise_mutation.log).

    RESIDUAL is bit-for-bit inert: with no charged particles present every wall loss
    ``nu*y[j]`` is identically zero, so a wall reactor's residual EQUALS a wall-less
    one's (``==``, not approx). A stray constant or wrong index in ``_apply_wall_terms``
    breaks this.

    JACOBIAN is NOT inert, and must not be claimed to be: the wall loss is ``nu*y[j]``,
    whose derivative in ``y[j]`` is ``nu`` regardless of ``y[j]``, so the linearization
    carries the wall slope even where the value is zero. The invariant here is exact
    instead: the wall reactor's Jacobian MINUS the wall-less one's equals the analytic
    wall linearization -- ``-nu`` on every charged diagonal and ``+gamma*nu`` from each
    cation into its recycle target. A wrong sign, a wrong nu, or a mis-placed recycle
    term breaks this."""
    r_wall, _, _ = _build_reactor(wall=True, with_chemistry=True)
    r_none, _, _ = _build_reactor(wall=False, with_chemistry=True)
    assert r_wall.has_wall is True
    assert r_none.has_wall is False

    ie, i_ar, i_arp = _indices(r_wall)
    # All-neutral: no electrons, no ions. The ionisation reaction is first order in the
    # electron, so it too is inert; the only wall-attributable difference is the wall term.
    y = np.zeros(r_wall.num_core_species, float)
    y[i_ar] = 1.0
    dydt = np.zeros(r_wall.num_core_species, float)

    dw, iw = r_wall.residual(0.0, y.copy(), dydt.copy())
    dn, in_ = r_none.residual(0.0, y.copy(), dydt.copy())
    assert np.array_equal(dw, dn), (
        "wall residual not inert on a chargeless state; differs by {0!r}".format(
            np.asarray(dw) - np.asarray(dn)))
    assert iw == in_

    pw = np.array(r_wall.jacobian(0.0, y.copy(), dydt.copy(), 1.0), float)
    pn = np.array(r_none.jacobian(0.0, y.copy(), dydt.copy(), 1.0), float)

    nu = r_wall.compute_nu_wall(y, r_wall.compute_volume(y))
    assert nu > 0.0
    expected = np.zeros_like(pw)
    gamma = r_wall.wall_recycling
    for j in range(r_wall.num_core_species):
        if r_wall.species_charges[j] == 0:
            continue
        expected[j, j] -= nu
        if j == ie:
            continue
        target = int(r_wall.wall_recycle_target[j])
        if target >= 0 and gamma > 0.0:
            expected[target, j] += gamma * nu
    diff = pw - pn
    assert np.allclose(diff, expected, rtol=1e-10, atol=1e-9), (
        "wall Jacobian contribution does not match the analytic linearization; "
        "max abs residual {0:.3e}".format(np.max(np.abs(diff - expected))))
    # And the difference is non-trivial: the wall genuinely changed the Jacobian.
    assert np.max(np.abs(diff)) > 1.0


# ---------------------------------------------------------------- item 8


def test_below_threshold_decays_to_extinction_without_negative_amounts():
    """Below threshold (Te = 0.5 eV, no external source) the electron
    population decays toward extinction and no negative species amount is
    produced anywhere on the way."""
    te_lo = 0.5
    r, _, _ = _build_reactor(te_ev=te_lo, wall=True, with_chemistry=True,
                             x_ion=1.0e-9, source=None)
    ie, i_ar, i_arp = _indices(r)
    y0 = r.y0.copy()
    V0 = r.compute_volume(y0)
    nu_wall = r.compute_nu_wall(y0, V0)
    n_ar = y0[i_ar] * constants.Na / V0
    n_e = y0[ie] * constants.Na / V0
    kinetics = VoronovEIArrhenius(Z=18, N=18, yaml_path_or_obj=VORONOV_YAML)
    k_ion = kinetics.get_rate_coefficient(te_lo * EV_TO_K) / constants.Na
    nu_ion = k_ion * n_ar
    assert nu_wall > nu_ion, (
        "this test's premise is that the wall wins below threshold; "
        "nu_wall={0!r}, nu_ion={1!r}".format(nu_wall, nu_ion))

    ts = np.logspace(-6, 0, 25)
    trajectory = []
    for t in ts:
        r.advance(t)
        trajectory.append(r.y.copy())
        assert np.all(r.y[:r.num_core_species] >= 0.0), (
            "negative species amount at t={0!r}: {1!r}".format(t, r.y[:r.num_core_species]))

    decay = trajectory[-1][ie] / y0[ie]
    assert decay < 1e-3


# ---------------------------------------------------------------- item 9

FD_FRACTIONS = (1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
FD_TOLERANCE = 1e-5


def _jacobian_scan(reactor, y, mutate=None):
    """Best relative disagreement between the analytic Jacobian and a central
    difference, minimised over step size. ``mutate`` is applied AFTER the
    analytic Jacobian is taken and BEFORE the differences."""
    n = reactor.num_core_species
    zeros = np.zeros(n, float)
    reactor.jacobian(0.0, y, zeros, 0.0)
    analytic = np.array(reactor.jacobian_matrix, float)
    if mutate is not None:
        mutate(reactor)
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
        rel = np.abs(analytic - fd) / scale
        best = min(best, rel.max())
    return best


@pytest.mark.parametrize("label,kwargs", [
    ("no wall at all (pre-existing equations)", dict(wall=False)),
    ("wall, integrated electron", dict(wall=True)),
    ("wall, quasineutral electron", dict(wall=True, quasineutral=True)),
    ("wall + external source", dict(wall=True, source=1.0e5)),
    ("gamma = 0.5 (partially pumping)", dict(wall=True, gamma=0.5)),
])
def test_analytic_jacobian_matches_finite_difference_scan(label, kwargs):
    """The analytic Jacobian agrees with a central finite difference, scanned
    over step size (a single step size is not diagnostic: truncation ~h^2,
    roundoff ~eps/h). Includes a no-wall control arm."""
    r, _, _ = _build_reactor(with_chemistry=True, x_ion=1.0e-7, **kwargs)
    y = _state_at(r, 1.0e-7)
    best = _jacobian_scan(r, y)
    assert best < FD_TOLERANCE, "{0}: best relative disagreement {1!r}".format(label, best)


def test_jacobian_scan_negative_control_can_fail():
    """The scan above must be able to fail: changing wall_recycling between
    the analytic Jacobian and the finite differences makes the two describe
    genuinely different operators."""
    r, _, _ = _build_reactor(wall=True, with_chemistry=True, x_ion=1.0e-7, gamma=0.5)
    y = _state_at(r, 1.0e-7)

    def _break_gamma(reactor):
        reactor.wall_recycling = 0.6

    best = _jacobian_scan(r, y, mutate=_break_gamma)
    assert best > FD_TOLERANCE


# ---------------------------------------------------------------- item 10


def test_ionisation_degree_ceiling_raises_naming_ionisation_degree():
    """An ionisation degree above max_ionisation_degree raises PlasmaStateError
    naming the ionisation degree; one below it does not raise.

    The ceiling is enforced by :meth:`check_wall_support`, called on ACCEPTED
    states after every step/advance -- deliberately not by
    :meth:`compute_nu_wall`, which the solver also evaluates at unphysical
    Newton trial states where a Python exception cannot propagate through the
    Fortran callback. See the docstring of ``compute_nu_wall`` in plasma.pyx.
    """
    r, _, _ = _build_reactor(wall=True, with_chemistry=False)
    y_hot = _state_at(r, 1.0e-2)   # alpha = 1e-2 > default ceiling 1e-3
    with pytest.raises(PlasmaStateError, match="ionisation degree"):
        r.check_wall_support(y_hot)

    y_ok = _state_at(r, 1.0e-4)    # alpha = 1e-4 < default ceiling 1e-3
    r.check_wall_support(y_ok)  # must not raise


# ---------------------------------------------------------------- item 11


def test_nonneutral_initial_state_raises_under_quasineutral_electron():
    """A non-charge-neutral initial composition with quasineutral_electron=True
    raises PlasmaStateError naming the net charge."""
    electron, ar, arp = _argon_species()
    bad = PlasmaReactor(
        (TGAS, 'K'), (P_NOMINAL, 'Pa'),
        {electron: 2.0e-8, arp: 1.0e-8, ar: 1.0 - 3.0e-8},
        (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
        quasineutral_electron=True,
        diffusion_length=(_diffusion_length(), 'm'),
        ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'))
    with pytest.raises(PlasmaStateError, match="net charge"):
        bad.initialize_model([electron, ar, arp], [], [], [])


# ---------------------------------------------------------------- item 12


def test_diffusion_length_and_mobility_must_be_declared_together():
    """Declaring a diffusion length without an ion mobility, or vice versa, raises."""
    electron, ar, arp = _argon_species()
    imf = {electron: 1.0e-6, arp: 1.0e-6, ar: 1.0 - 2.0e-6}

    with pytest.raises(PlasmaStateError):
        PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf,
                      (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                      diffusion_length=(_diffusion_length(), 'm'))

    electron2, ar2, arp2 = _argon_species()
    imf2 = {electron2: 1.0e-6, arp2: 1.0e-6, ar2: 1.0 - 2.0e-6}
    with pytest.raises(PlasmaStateError):
        PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf2,
                      (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                      ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'))


# ---------------------------------------------------------------- item 13


def test_wall_parameter_bounds_each_raise():
    """wall_recycling outside [0,1], a non-positive diffusion length, a
    non-positive mobility, and a negative ionisation source each raise."""

    def _make(**overrides):
        electron, ar, arp = _argon_species()
        imf = {electron: 1.0e-6, arp: 1.0e-6, ar: 1.0 - 2.0e-6}
        kwargs = dict(diffusion_length=(_diffusion_length(), 'm'),
                      ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'))
        kwargs.update(overrides)
        return PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf,
                             (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1,
                             termination=[], **kwargs)

    with pytest.raises(PlasmaStateError):
        _make(wall_recycling=1.5)

    with pytest.raises(PlasmaStateError):
        _make(diffusion_length=(-_diffusion_length(), 'm'))

    with pytest.raises(PlasmaStateError):
        _make(ion_reduced_mobility=(-MU0_AR_IN_AR, 'm^2/(V*s)'))

    with pytest.raises(PlasmaStateError):
        _make(ionisation_source=(-1.0e5, 'm^-3/s'))


# ---------------------------------------------------------------- item 14


def test_quasineutral_electron_row_is_net_charge_with_no_dydt_dependence():
    """With quasineutral_electron=True, the electron row of the residual is
    the net charge and carries no dydt dependence; and the same reactor
    integrates successfully."""
    r, _, _ = _build_reactor(wall=True, with_chemistry=True, x_ion=1.0e-7,
                             quasineutral=True)
    ie, i_ar, i_arp = _indices(r)
    y = _state_at(r, 1.0e-7)

    dydt_a = np.zeros(r.num_core_species, float)
    delta_a, _ = r.residual(0.0, y.copy(), dydt_a)
    net_charge = float(np.sum(r.species_charges * y))
    assert delta_a[ie] == net_charge

    dydt_b = dydt_a.copy()
    dydt_b[ie] = 1.0e7   # a large, arbitrary perturbation of dydt at the electron row
    delta_b, _ = r.residual(0.0, y.copy(), dydt_b)
    assert delta_b[ie] == delta_a[ie]

    # the same reactor integrates successfully
    r.advance(1.0e-8)


# ---------------------------------------------------------------- item 15


def test_charge_row_scale_is_positive_finite_and_leaves_the_solution_unchanged():
    """charge_row_scale is strictly positive and finite once initialize_model
    has run with quasineutral_electron=True, and is exactly 1.0 when it is
    False. It is an exact row scaling of the DAE residual and Jacobian, so it
    must not change the solution: two otherwise-identical reactors, given two
    different positive charge_row_scale values right after initialize_model,
    integrate to the same trajectory within solver tolerance."""
    r_true, _, _ = _build_reactor(wall=True, with_chemistry=True, x_ion=1.0e-7,
                                  quasineutral=True)
    assert r_true.charge_row_scale > 0.0
    assert np.isfinite(r_true.charge_row_scale)

    r_false, _, _ = _build_reactor(wall=True, with_chemistry=True, x_ion=1.0e-7,
                                   quasineutral=False)
    assert r_false.charge_row_scale == 1.0

    ra, _, _ = _build_reactor(wall=True, with_chemistry=True, x_ion=1.0e-7,
                              quasineutral=True)
    rb, _, _ = _build_reactor(wall=True, with_chemistry=True, x_ion=1.0e-7,
                              quasineutral=True)
    ra.charge_row_scale = 1.0
    rb.charge_row_scale = 1.0e6
    assert ra.charge_row_scale != rb.charge_row_scale

    for t in np.logspace(-9, -7, 6):
        ra.advance(t)
        rb.advance(t)
        ya = ra.y[:ra.num_core_species]
        yb = rb.y[:rb.num_core_species]
        scale = np.maximum(np.abs(ya), np.abs(yb))
        scale[scale == 0.0] = 1.0
        assert np.all(np.abs(ya - yb) / scale < 1e-6), (
            "t={0!r}: ya={1!r} yb={2!r}".format(t, ya, yb))


# ---------------------------------------------------------------- item 17


def _wall_neutral_moles(reactor, y):
    """Total moles of neutral heavy species in ``y`` -- the quantity nu_wall
    divides by, before any floor is applied."""
    return float(sum(y[j] for j in range(reactor.num_core_species)
                     if reactor.neutral_heavy_mask[j]))


def _neutral_density(reactor, y, V):
    """The neutral-heavy number density (m^-3) -- the intensive quantity nu_wall uses,
    and the quantity the floor is defined in after round 90."""
    return _wall_neutral_moles(reactor, y) * constants.Na / V


def test_jacobian_matches_residual_below_the_neutral_floor():
    """The Jacobian must differentiate the residual the solver ACTUALLY evaluates,
    including on the clamped branch below the neutral-DENSITY floor.

    ``compute_nu_wall`` floors the neutral number density, so below the floor nu_wall is
    a CONSTANT -- pinned to the floor density, with no dependence on the neutral amount
    AND none on V -- and its derivative with respect to every species is exactly zero. A
    Jacobian that keeps either the -1/y_neutral term or the dV/V term there reports a
    dependence the residual does not have. (The old moles clamp pinned the neutral amount
    but left n_neutral = floor*Na/V, so nu kept a V dependence and the Jacobian kept the
    dV/V term; pinning the density removes both.)

    The floor is intensive, so a below-floor state is reached by driving the neutral
    DENSITY down -- a charge-dominated composition that inflates the volume -- not by
    shrinking the neutral moles at fixed density, which an isobaric reactor holds pinned
    at ~P/kT.
    """
    r, _, _ = _build_reactor(wall=True, with_chemistry=False)
    ie, i_ar, i_arp = _indices(r)
    floor = r.wall_neutral_density_floor
    assert floor > 0.0, "the floor must exist for this test to mean anything"

    # A charge-dominated state whose neutral number density is an order of magnitude
    # BELOW the floor, so compute_nu_wall is on its clamped branch.
    y = np.zeros(r.num_core_species, float)
    y[ie] = 1.0
    y[i_arp] = 1.0
    y[i_ar] = 1.0                                  # provisional, to fix the volume scale
    V0 = r.compute_volume(y)
    y[i_ar] = 0.1 * floor * V0 / constants.Na
    V = r.compute_volume(y)
    assert _neutral_density(r, y, V) < floor, "state is not below the density floor"

    nu_here = r.compute_nu_wall(y, V)

    # Precondition: nu is CLAMPED -- constant in y. Two below-floor states give the
    # identical nu (not merely the identical nu/V the old moles clamp gave). That
    # constancy is exactly what the Jacobian must reproduce as a zero derivative.
    y_half = y.copy()
    y_half[i_ar] = 0.5 * y[i_ar]
    V_half = r.compute_volume(y_half)
    assert _neutral_density(r, y_half, V_half) < floor
    assert np.isclose(r.compute_nu_wall(y_half, V_half), nu_here, rtol=1e-12, atol=0.0), (
        "precondition failed: nu_wall is not clamped constant below the floor")

    zeros = np.zeros(r.num_core_species, float)
    jac = np.asarray(r.jacobian(0.0, y, zeros, 0.0), float)

    # Full-matrix central finite difference of the residual, every column, at the
    # below-floor state. Small relative steps keep every perturbed state clamped.
    worst = 0.0
    for col in range(r.num_core_species):
        h = 1.0e-6 * abs(y[col]) if y[col] != 0.0 else 1.0e-6 * y[i_ar]
        yp = y.copy(); yp[col] += h
        ym = y.copy(); ym[col] -= h
        assert _neutral_density(r, yp, r.compute_volume(yp)) < floor, "step left the clamp"
        assert _neutral_density(r, ym, r.compute_volume(ym)) < floor, "step left the clamp"
        rp = np.asarray(r.residual(0.0, yp, zeros)[0], float)
        rm = np.asarray(r.residual(0.0, ym, zeros)[0], float)
        fd = (rp - rm) / (2.0 * h)
        err = np.max(np.abs(jac[:, col] - fd)) / max(1.0, np.max(np.abs(fd)))
        worst = max(worst, err)

    assert worst < 1.0e-6, (
        "analytic Jacobian disagrees with the residual it claims to differentiate "
        "below the neutral density floor: relative error {0:.3e}".format(worst))


def test_jacobian_finite_when_neutrals_are_exhausted():
    """A trial state with NO neutrals left must still produce a finite Jacobian.

    ``compute_nu_wall`` never raises there by design -- the floor keeps it finite
    so DASPK can reject the state on its own terms. The Jacobian has to honour the
    same contract; dividing by an unfloored y_neutral makes it raise instead, and
    an exception inside the Fortran callback is exactly the failure mode the
    accepted-step domain check was moved out of the residual to avoid.
    """
    r, _, _ = _build_reactor(wall=True, with_chemistry=False)
    ie, i_ar, i_arp = _indices(r)

    y = _state_at(r, 1.0e-6)
    y[i_ar] = 0.0
    assert _wall_neutral_moles(r, y) == 0.0

    V = r.compute_volume(y)
    assert np.isfinite(r.compute_nu_wall(y, V)), "residual side already broken"

    zeros = np.zeros(r.num_core_species, float)
    res = np.asarray(r.residual(0.0, y, zeros)[0], float)
    assert np.all(np.isfinite(res))

    jac = np.asarray(r.jacobian(0.0, y, zeros, 0.0), float)
    assert np.all(np.isfinite(jac)), "Jacobian is not finite with no neutrals left"


# ---------------------------------------------------------------- item 19


_WALL_ATTRS_BY_VALUE = (
    # (attribute, how to read a comparable value out of it)
    ('diffusion_length', lambda r: None if r.diffusion_length is None
     else r.diffusion_length.value_si),
    ('ion_reduced_mobility', lambda r: None if r.ion_reduced_mobility is None
     else r.ion_reduced_mobility.value_si),
    ('ionisation_source', lambda r: None if r.ionisation_source is None
     else r.ionisation_source.value_si),
    ('mobility_reference_density', lambda r: r.mobility_reference_density),
    ('wall_recycling', lambda r: r.wall_recycling),
    ('max_ionisation_degree', lambda r: r.max_ionisation_degree),
    ('quasineutral_electron', lambda r: r.quasineutral_electron),
    ('has_wall', lambda r: r.has_wall),
)


def _wall_fingerprint(reactor):
    return {name: read(reactor) for name, read in _WALL_ATTRS_BY_VALUE}


@pytest.mark.parametrize("copy_name", ["pickle", "deepcopy"])
def test_round_trip_preserves_the_wall_by_value(copy_name):
    """A copied wall-enabled reactor must still HAVE its wall, checked by value.

    ``__reduce__`` enumerates the constructor arguments by hand. A hand-written
    enumeration that predates a parameter silently drops it: the copy is a valid
    PlasmaReactor, it initialises, it integrates, and it has no wall at all. So
    asserting that the object survives the trip proves nothing -- the assertion
    has to read the wall parameters back out and compare them.

    deepcopy is covered alongside pickle because it routes through
    ``__reduce_ex__`` as well, so one omission breaks both.
    """
    source_rate = 1.0e5
    gamma = 0.5
    r, _, _ = _build_reactor(wall=True, with_chemistry=False, source=source_rate,
                             gamma=gamma, quasineutral=True, max_alpha=7.5e-4)

    before = _wall_fingerprint(r)
    # Guard the test against itself: the fixture must really carry a wall, or a
    # copy that drops it would trivially "match".
    assert before['has_wall'] is True
    assert before['diffusion_length'] is not None and before['diffusion_length'] > 0.0
    assert before['ion_reduced_mobility'] is not None and before['ion_reduced_mobility'] > 0.0
    assert before['ionisation_source'] == source_rate
    assert before['wall_recycling'] == gamma
    assert before['max_ionisation_degree'] == 7.5e-4
    assert before['quasineutral_electron'] is True

    if copy_name == "pickle":
        clone = pickle.loads(pickle.dumps(r))
    else:
        clone = copy.deepcopy(r)

    after = _wall_fingerprint(clone)

    mismatched = {k: (before[k], after[k]) for k in before if before[k] != after[k]}
    assert not mismatched, (
        "{0} round trip did not preserve the wall. before -> after: {1}".format(
            copy_name, mismatched))


@pytest.mark.parametrize("copy_name", ["pickle", "deepcopy"])
def test_round_trip_of_a_wall_less_reactor_stays_wall_less(copy_name):
    """The negative arm: a reactor with no wall must not acquire one.

    Without this, a ``__reduce__`` that hard-coded a wall would pass the test
    above. The two together pin the round trip in both directions.
    """
    r, _, _ = _build_reactor(wall=False, with_chemistry=False)
    assert r.has_wall is False

    clone = (pickle.loads(pickle.dumps(r)) if copy_name == "pickle"
             else copy.deepcopy(r))

    assert clone.has_wall is False
    assert clone.diffusion_length is None
    assert clone.ion_reduced_mobility is None


def _plasma_pyx_source():
    """The .pyx this module is compiled from.

    Read from source deliberately: ``PlasmaReactor`` is a Cython ``cdef class``,
    and ``inspect.signature`` reports its ``__init__`` as ``(*args, **kwargs)``
    -- the real parameter names are not introspectable at runtime at all. A
    first draft of the test below compared argument COUNTS via ``inspect`` and
    was red for that reason rather than for the defect, which is exactly the
    false-positive this file is supposed to avoid.
    """
    path = os.path.join(os.path.dirname(rmgpy.solver.plasma.__file__), 'plasma.pyx')
    with open(path) as f:
        return f.read()


def test_reduce_enumerates_every_constructor_parameter():
    """``__reduce__`` must mention every parameter ``__init__`` accepts.

    The value checks above catch today's omission. This one is aimed at the NEXT
    one: a parameter added to ``__init__`` without being added to ``__reduce__``
    fails here, naming the parameter, instead of silently producing copies that
    have quietly lost it.

    This is a source-level check, and says so: see ``_plasma_pyx_source``.
    """
    src = _plasma_pyx_source()

    init = re.search(r'\n    def __init__\(self,(.*?)\):\n', src, re.S)
    assert init, "could not locate PlasmaReactor.__init__ in plasma.pyx"
    params = re.findall(r'([A-Za-z_][A-Za-z_0-9]*)\s*(?:=[^,]*)?(?:,|$)',
                        re.sub(r'#.*', '', init.group(1)))
    params = [p for p in dict.fromkeys(params) if p not in ('self',)]
    assert 'diffusion_length' in params, (
        "parameter scrape failed; got {0!r}".format(params))

    reduce_body = re.search(r'\n    def __reduce__\(self\):(.*?)\n    (?:cpdef|def|cdef) ',
                            src, re.S)
    assert reduce_body, "could not locate PlasmaReactor.__reduce__ in plasma.pyx"
    body = reduce_body.group(1)
    # Search only the RETURNED reconstruction tuple, with the method docstring and all
    # comments stripped first. A parameter named merely in a comment or in the docstring of
    # __reduce__ (both are present -- the docstring lists parameters by name, and the return
    # tuple carries an explanatory comment) must NOT satisfy the check: the property is that
    # the parameter is actually CARRIED in the reconstruction arguments, not that its name
    # appears somewhere in the method. This is what makes the test more than a lexical
    # surrogate for the property it names.
    body = re.sub(r'(?s)"""(?:.*?)"""', '', body)     # docstring
    body = re.sub(r'#.*', '', body)                   # comments
    ret = re.search(r'return\s*\((.*)\)\s*\Z', body.strip(), re.S)
    assert ret, "could not locate the return tuple of __reduce__ after stripping docstring/comments"
    returned = ret.group(1)

    missing = [p for p in params if not re.search(r'\bself\.%s\b' % re.escape(p), returned)]
    assert not missing, (
        "__reduce__ does not carry {0} of __init__'s parameters in its returned "
        "reconstruction tuple, so a pickle or deepcopy silently drops them: {1}".format(
            len(missing), missing))


# ---------------------------------------------------------------- I-246 rework
#
# The ground/metastable resolution, the thermo precondition, the near-degeneracy
# threshold, the wall_neutralization_products declaration, the neutral-floor and
# electron-population validators, the single-cation enforcement, and the direct
# diffusionLength dimension check.
#
# Two kinds of test live here, and which kind a test is can be read off its DOCSTRING --
# the brief is scoped by an observable marker, not asserted in bulk over every function.
#
#  - A DEFECT REPRODUCTION names, in its docstring, the finding and the round it closes
#    ("Round 96 HIGH 1", "MEDIUM: ...", "HIGH 3: ..."). THAT is the red-state claim, and it
#    is banked: the test was run red on the built module BEFORE its fix, in that round's
#    before.log (round90/93/96/100_before.log). If a test names a round/finding, its red
#    state exists in that round's log; if it does not, it makes no such claim.
#  - An INVARIANT / PROPERTY / ACCEPTANCE check -- a closed-form re-derivation of nu_wall, a
#    geometry scaling law, a charge/heavy-atom conservation statement, a parameterisation
#    invariance, a "the four verdicts hold together" acceptance -- asserts a positive
#    property directly. It has no "red state" in the defect sense and NONE is claimed for
#    it. Such a test either states the property it holds or carries no round/finding tag.
#
# So the standard is met by construction: the only tests that claim a banked red are the
# ones whose docstring names a finding, and those are exactly the ones reproduced red in a
# before.log. A reader auditing the claim checks the docstring tag against the round log,
# not a separate per-test manifest that could fall out of step with either.

EV_J_PER_MOL = 96485.33212


def _argon_thermo(excitation_eV):
    """Monatomic-argon thermo with an electronic offset. The offset is all that the
    ground-state rule reads; Cp and S are the same for both states, so the enthalpy
    difference is purely the excitation energy."""
    Cp = 2.5 * constants.R
    return ThermoData(
        Tdata=([298, 400, 600, 800, 1000, 1500, 2000], 'K'),
        Cpdata=([Cp] * 7, 'J/(mol*K)'),
        H298=(excitation_eV * EV_J_PER_MOL / 1000.0, 'kJ/mol'),
        S298=(154.8, 'J/(mol*K)'))


def _metastable_species(label='Ar*', excitation_eV=None):
    """A second neutral argon state, constructed locally (no database). It differs
    from the ground state by two unpaired electrons (multiplicity 3), so it is a
    distinct species that nonetheless shares argon's heavy composition."""
    s = Species(label=label).from_adjacency_list('1 Ar u2 p3 c0')
    if excitation_eV is not None:
        s.thermo = _argon_thermo(excitation_eV)
    return s


def _ground_species(label='Ar', excitation_eV=None):
    s = Species(label=label).from_adjacency_list('1 Ar u0 p4 c0')
    if excitation_eV is not None:
        s.thermo = _argon_thermo(excitation_eV)
    return s


def _metastable_reactor(ground_eV=0.0, meta_eV=11.5, meta_label='Ar*', gamma=1.0,
                        neutralization=None, with_meta_thermo=True,
                        with_ground_thermo=True):
    """A wall reactor whose core carries ground Ar, a second neutral Ar state, Ar+
    and e-. Thermo is attached locally so the ground-state energy rule can run."""
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ground = _ground_species('Ar', ground_eV if with_ground_thermo else None)
    meta = _metastable_species(meta_label, meta_eV if with_meta_thermo else None)
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    imf = {electron: 1.0e-6, arp: 1.0e-6, meta: 1.0e-6, ground: 1.0 - 3.0e-6}
    kwargs = dict(diffusion_length=(_diffusion_length(), 'm'),
                  ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                  wall_recycling=gamma)
    if neutralization is not None:
        kwargs['wall_neutralization_products'] = neutralization
    reactor = PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf, (TE_NOMINAL_EV * EV_TO_K, 'K'),
                            n_sims=1, termination=[], **kwargs)
    core = [electron, ground, meta, arp]
    reactor.initialize_model(core, [], [], [])
    return reactor, core


def _cation_recycle_label(reactor, core):
    z = reactor.species_charges
    i_arp = [j for j in range(len(z)) if z[j] == 1 and j != reactor.electron_index][0]
    return core[int(reactor.wall_recycle_target[i_arp])].label


def test_wall_requires_a_declaration_for_two_electronic_states():
    """Ground Ar, metastable Ar*, Ar+, e-, wall_recycling=1.0, and NO declaration.
    Two neutral states share Ar+'s heavy skeleton, and nothing in the deck certifies
    which is the ground-state product the wall returns: formation enthalpy orders them
    but an excited-only deck (true ground absent) is indistinguishable from
    ground+metastable. So the wall REFUSES and names the declaration, rather than
    picking the lower enthalpy. Energy is no longer a vote about identity -- it survives
    only in the wall-energy interface. (Round 79 auto-resolved this deck to ground; that
    was the narrowing round 83 closes.)"""
    with pytest.raises(PlasmaStateError) as exc:
        _metastable_reactor(ground_eV=0.0, meta_eV=11.5, gamma=1.0)
    msg = str(exc.value)
    assert 'ambiguous' in msg and 'heavy skeleton' in msg
    assert 'wallNeutralizationProducts' in msg


def test_wall_ground_state_charge_and_particle_conservation_with_metastable():
    """With the metastable present and gamma=1 and the product DECLARED, the wall
    conserves both net charge and heavy atoms; the metastable itself takes no wall
    flux. (The declaration is required now that two electronic states are present.)"""
    reactor, core = _metastable_reactor(gamma=1.0, neutralization={'Ar+': 'Ar'})
    z = reactor.species_charges
    ie = reactor.electron_index
    i_ground = [j for j in range(len(z)) if z[j] == 0
                and core[j].label == 'Ar'][0]
    i_meta = [j for j in range(len(z)) if z[j] == 0 and core[j].label == 'Ar*'][0]
    i_arp = [j for j in range(len(z)) if z[j] == 1 and j != ie][0]
    y = np.zeros(reactor.num_core_species, float)
    y[i_ground] = 1.0
    y[i_meta] = 1.0e-2
    y[i_arp] = 1.0e-4
    y[ie] = 1.0e-4
    reactor.residual(0.0, y.copy(), np.zeros_like(y))
    wall = np.asarray(reactor.wall_loss_rates, float)
    nu = reactor.nu_wall
    assert wall[i_meta] == 0.0
    assert np.isclose(wall[i_arp], -nu * y[i_arp], rtol=1e-12, atol=0.0)
    assert np.isclose(wall[i_ground], nu * y[i_arp], rtol=1e-12, atol=0.0)
    net_charge_rate = sum(z[j] * wall[j] for j in range(len(z)))
    assert abs(net_charge_rate) <= 1e-18 * max(1.0, abs(nu))
    net_heavy_rate = sum(wall[j] for j in range(len(z)) if j != ie)
    assert np.isclose(net_heavy_rate, 0.0, atol=1e-18)


def test_wall_refuses_two_neutral_states_by_multiplicity_not_thermo():
    """Two neutral Ar states sharing Ar+'s skeleton are refused because there are TWO,
    not because their thermo is missing: the product is a declaration, not an energy
    inference, so the refusal fires whether or not thermo is attached (here it is
    absent) and BEFORE any enthalpy is read. It names the declaration to supply."""
    with pytest.raises(PlasmaStateError) as exc:
        _metastable_reactor(with_meta_thermo=False, with_ground_thermo=False)
    msg = str(exc.value)
    assert 'ambiguous' in msg and 'heavy skeleton' in msg
    assert 'wallNeutralizationProducts' in msg


def test_wall_refuses_two_states_even_across_a_clear_energy_gap():
    """A clear energy gap does NOT license a pick. Two neutral states well separated in
    enthalpy are still refused, because a clear gap between two PRESENT states is
    exactly what an excited-only deck (5 and 11.5 eV, true ground absent) also shows --
    the ordering is real but the anchor to the absolute ground is not. Only a
    declaration resolves it. (Round 79 resolved this to the lower state; round 83 does
    not, because the gap cannot tell ground+excited from excited+excited.)"""
    rt_ev = constants.R * TGAS / EV_J_PER_MOL
    with pytest.raises(PlasmaStateError) as exc:
        _metastable_reactor(ground_eV=0.0, meta_eV=2.0 * rt_ev, meta_label='Ar_hi')
    assert 'ambiguous' in str(exc.value) and 'heavy skeleton' in str(exc.value)


def test_wall_refuses_two_states_regardless_of_their_separation():
    """Whatever the energy separation -- here a near-degenerate pair -- two neutral
    states sharing the ion's skeleton are refused with the declaration named. Energy no
    longer picks at any separation; multiplicity alone requires the declaration."""
    rt_ev = constants.R * TGAS / EV_J_PER_MOL
    with pytest.raises(PlasmaStateError) as exc:
        _metastable_reactor(ground_eV=0.0, meta_eV=0.5 * rt_ev, meta_label='Ar_near')
    msg = str(exc.value)
    assert 'ambiguous' in msg and 'heavy skeleton' in msg
    assert 'wallNeutralizationProducts' in msg


def test_wall_neutralization_products_declaration_resolves_ambiguity():
    """The declaration names the ground-state product explicitly and wins."""
    reactor, core = _metastable_reactor(neutralization={'Ar+': 'Ar'})
    assert _cation_recycle_label(reactor, core) == 'Ar'


def test_wall_neutralization_products_missing_neutral_is_refused():
    """A declaration naming a neutral that is not in the core is refused by name."""
    with pytest.raises(PlasmaStateError) as exc:
        _metastable_reactor(neutralization={'Ar+': 'Ar_ground_not_here'})
    assert 'Ar_ground_not_here' in str(exc.value)


def test_check_wall_support_refuses_subfloor_inventory():
    """An accepted state whose neutral number DENSITY is at or below the intensive floor
    is refused. After round 90 the floor is a fraction of the reference density, so the
    refusable state is a genuinely depleted (over-ionised) gas -- reached by inflating the
    volume with a charge-dominated composition -- and NOT a small-but-normal-density
    inventory, which the old extensive moles floor wrongly refused on its history."""
    r, _, _ = _build_reactor(wall=True, gamma=0.0, with_chemistry=False)
    ie, i_ar, i_arp = _indices(r)
    floor = r.wall_neutral_density_floor
    y = np.zeros(r.num_core_species, float)
    y[ie] = 1.0
    y[i_arp] = 1.0
    y[i_ar] = 1.0                                  # provisional, to fix the volume scale
    V0 = r.compute_volume(y)
    y[i_ar] = 0.1 * floor * V0 / constants.Na      # drive n_neutral an order below floor
    assert _neutral_density(r, y, r.compute_volume(y)) < floor
    with pytest.raises(PlasmaStateError) as exc:
        r.check_wall_support(y)
    assert 'numerical floor' in str(exc.value)


def test_check_wall_support_refuses_nonfinite_electron():
    """A NaN or negative electron population is refused, not carried into the
    ionisation-degree ratio (nan > ceiling is False, so it used to pass)."""
    r, _, _ = _build_reactor(wall=True, gamma=0.0, with_chemistry=False)
    ie, i_ar, i_arp = _indices(r)
    y = np.zeros(r.num_core_species, float)
    y[i_ar] = 1.0
    y[ie] = float('nan')
    with pytest.raises(PlasmaStateError) as exc:
        r.check_wall_support(y)
    assert 'electron population of' in str(exc.value)
    y[ie] = -1.0e-6
    with pytest.raises(PlasmaStateError):
        r.check_wall_support(y)


def test_wall_refuses_anion_in_core():
    """A single mobility applied to an anion is the wrong sign of physics (the sheath
    confines anions), so an anion in the core with a wall is refused."""
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ground = _ground_species('Ar')
    cl = Species(label='Cl').from_adjacency_list('1 Cl u1 p3 c0')
    cln = Species(label='Cl-').from_adjacency_list('1 Cl u0 p4 c-1')
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    imf = {electron: 1e-6, arp: 1e-6, cln: 1e-6, cl: 1e-6, ground: 1.0 - 4e-6}
    # Ar and Cl are distinct heavy skeletons, so the neutral bath is a mixture; opt into
    # the single-bath approximation so this deck reaches the anion check it is about.
    reactor = PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf,
                            (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                            diffusion_length=(_diffusion_length(), 'm'),
                            ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                            wall_single_bath_approximation=True)
    with pytest.raises(PlasmaStateError) as exc:
        reactor.initialize_model([electron, ground, cl, cln, arp], [], [], [])
    assert 'does not support anions' in str(exc.value)


def test_direct_diffusion_length_dimension_is_checked():
    """A directly-stated diffusionLength must be a length: (2, 's') is refused rather
    than silently reinterpreted as (2, 'm')."""
    from rmgpy.rmg.input import _plasma_wall_kwargs
    from rmgpy.exceptions import InputError
    with pytest.raises(InputError) as exc:
        _plasma_wall_kwargs({'diffusionLength': (2.0, 's')},
                            (MU0_AR_IN_AR, 'm^2/(V*s)'), None, 1.0, None, None)
    assert 'diffusionLength' in str(exc.value) and 'length' in str(exc.value)
    # and the valid spelling still resolves
    kwargs = _plasma_wall_kwargs({'diffusionLength': (2.03, 'cm')},
                                 (MU0_AR_IN_AR, 'm^2/(V*s)'), None, 1.0, None, None)
    assert abs(kwargs['diffusion_length'][0] - 0.0203) < 1e-12


# ============================================================ round-79 repairs
# Each test below reproduces a defect the earlier rework exposed (see
# docs/i246-ambipolar-wall-operator/rework.md). The recurring root is a species
# identified by a key coarser than the physics distinguishes: the wall recycle and
# the ionisation source both matched on element COUNT, which collides constitutional
# isomers and cannot separate a molecular ion from an isomeric neutral. The correct
# key is the heavy-atom skeleton -- the standard InChI truncated before its charge
# layer -- which unifies electronic states (Ar and Ar* share it) yet separates
# isomers (DME and ethanol do not).


def _thermo_with_h298(h298_kj):
    """Minimal thermo carrying a chosen formation enthalpy at 298 K; Cp and S are
    placeholders, since only the enthalpy ordering is read by the ground-state rule."""
    return ThermoData(
        Tdata=([298, 400, 600, 800, 1000, 1500, 2000], 'K'),
        Cpdata=([4.0 * constants.R] * 7, 'J/(mol*K)'),
        H298=(h298_kj, 'kJ/mol'),
        S298=(200.0, 'J/(mol*K)'))


def _isomer_reactor(gamma=1.0, dme_kj=-184.0, eth_kj=-235.0, neutralization=None):
    """Core carrying a molecular cation DME+ alongside TWO same-formula (C2H6O)
    neutrals: its own neutral DME and the constitutional isomer ethanol, ethanol the
    lower in formation enthalpy so a lowest-enthalpy rule would 'crown' it."""
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    dme = Species(label='DME').from_smiles('COC')
    dme.thermo = _thermo_with_h298(dme_kj)
    eth = Species(label='EtOH').from_smiles('CCO')
    eth.thermo = _thermo_with_h298(eth_kj)
    dmep = Species(label='DME+').from_smiles('[CH3][O+][CH3]')
    imf = {electron: 1.0e-6, dmep: 1.0e-6, eth: 1.0e-6, dme: 1.0 - 3.0e-6}
    kwargs = dict(diffusion_length=(_diffusion_length(), 'm'),
                  ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                  wall_recycling=gamma,
                  # DME and ethanol are distinct skeletons: opt into the single-bath
                  # approximation so the recycle-target behaviour under test is reached.
                  wall_single_bath_approximation=True)
    if neutralization is not None:
        kwargs['wall_neutralization_products'] = neutralization
    reactor = PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf,
                            (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[], **kwargs)
    core = [electron, dme, eth, dmep]
    reactor.initialize_model(core, [], [], [])
    return reactor, core


def _cation_index(reactor):
    z = reactor.species_charges
    return [j for j in range(len(z)) if z[j] == 1 and j != reactor.electron_index][0]


def test_wall_does_not_transmute_a_molecular_ion_across_isomers():
    """HIGH 1: DME+ neutralises to DME, never to the constitutional isomer ethanol
    just because ethanol's enthalpy is lower. Element-count matching admitted the
    isomer and 'lowest enthalpy' then crowned it -- a transmutation. The skeleton
    (InChI) excludes ethanol from the candidate set entirely."""
    reactor, core = _isomer_reactor()
    assert core[int(reactor.wall_recycle_target[_cation_index(reactor)])].label == 'DME'


def test_multi_state_deck_refuses_before_reading_enthalpy():
    """Round-79 HIGH 1 was 'NaN enthalpy initialises on a coin-toss'. Under round 83
    that failure mode is gone by construction: identity never reads enthalpy, so a NaN
    formation enthalpy cannot influence the pick. This deck -- two neutral Ar states,
    the ground one carrying a NaN enthalpy -- refuses on MULTIPLICITY, before any
    enthalpy is examined, and the message is about the ambiguity, not thermochemistry.
    (Whether NaN thermo is handled correctly where it DOES still matter -- the
    wall-energy interface -- is test_nonfinite_recycle_thermo_leaves_energy_unavailable
    below.)"""
    with pytest.raises(PlasmaStateError) as exc:
        _metastable_reactor(ground_eV=float('nan'), meta_eV=11.5)
    msg = str(exc.value)
    assert 'ambiguous' in msg and 'heavy skeleton' in msg
    assert 'thermochemistry' not in msg


def test_wall_refuses_electron_without_a_cation():
    """HIGH 2: a wall with an electron but NO positive ion would remove electrons
    with no charge-conserving partner, driving the net charge -- not floating-wall
    behaviour. 'Exactly one cation' rejected two and accepted zero; zero is refused."""
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = _ground_species('Ar')
    imf = {electron: 1.0e-6, ar: 1.0 - 1.0e-6}
    reactor = PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf,
                            (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                            diffusion_length=(_diffusion_length(), 'm'),
                            ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'))
    with pytest.raises(PlasmaStateError) as exc:
        reactor.initialize_model([electron, ar], [], [], [])
    assert 'cation' in str(exc.value) or 'positive ion' in str(exc.value)


def test_ceiling_uses_ion_inventory_not_electron_alone():
    """HIGH 2: the ionisation-degree ceiling gates on the charged inventory, not n_e
    alone. 1% Ar+ with 1 ppm electrons is a 1e-2 ion fraction and must be refused
    under a 1e-3 ceiling, though n_e/n_neutral = 1e-6 slips under it."""
    r, _, _ = _build_reactor(wall=True, gamma=0.0, with_chemistry=False, max_alpha=1.0e-3)
    ie, i_ar, i_arp = _indices(r)
    y = np.zeros(r.num_core_species, float)
    y[i_ar] = 1.0
    y[i_arp] = 1.0e-2
    y[ie] = 1.0e-6
    with pytest.raises(PlasmaStateError) as exc:
        r.check_wall_support(y)
    assert 'ionisation degree' in str(exc.value)


def test_declared_source_is_delivered_in_full_in_a_mixture():
    """HIGH 3: the external pair source is apportioned over IONISABLE neutrals. A
    non-ionisable bath gas (He, with no He+ in the core) must not sit in the
    denominator and silently swallow half the declared source."""
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = _ground_species('Ar')
    he = Species(label='He').from_adjacency_list('1 He u0 p1 c0')
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    source = 1.0e18
    imf = {electron: 1.0e-6, arp: 1.0e-6, he: 0.5, ar: 0.5 - 2.0e-6}
    reactor = PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf,
                            (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                            diffusion_length=(_diffusion_length(), 'm'),
                            ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                            wall_recycling=0.0, ionisation_source=(source, 'm^-3/s'),
                            wall_single_bath_approximation=True)
    core = [electron, ar, he, arp]
    reactor.initialize_model(core, [], [], [])
    z = reactor.species_charges
    ie = reactor.electron_index
    i_ar = [j for j in range(len(z)) if z[j] == 0 and core[j].label == 'Ar'][0]
    i_he = [j for j in range(len(z)) if z[j] == 0 and core[j].label == 'He'][0]
    i_arp = [j for j in range(len(z)) if z[j] == 1][0]
    y = np.zeros(reactor.num_core_species, float)
    y[i_ar] = 0.5
    y[i_he] = 0.5
    V = reactor.compute_volume(y)
    delta, _ = reactor.residual(0.0, y.copy(), np.zeros_like(y))
    source_total = source * V / constants.Na       # mol of pairs per second
    # The pair source is one whole reaction Ar -> Ar+ + e-: it produces the electron AND the
    # cation and consumes the neutral, all at the same rate. Asserting only the electron leg
    # would pass a source that created electrons from nothing without an ion or a consumed
    # neutral. At this all-neutral state no charged species is present, so the wall loss (which
    # scales with n_e / n_ion) is zero and each leg equals the full source exactly.
    assert np.isclose(delta[ie], source_total, rtol=1e-9, atol=0.0)      # electron produced
    assert np.isclose(delta[i_arp], source_total, rtol=1e-9, atol=0.0)   # cation produced
    assert np.isclose(delta[i_ar], -source_total, rtol=1e-9, atol=0.0)   # neutral consumed
    # He is not ionisable: it neither produces nor consumes, so it stays out of the balance.
    assert np.isclose(delta[i_he], 0.0, atol=1e-30)


def test_direct_construction_checks_diffusion_length_dimension():
    """MEDIUM: direct PlasmaReactor(...) construction must reject a diffusion length
    that is not a length -- (2, 's') became a 2 m length silently, the input-file path
    already guards this."""
    electron, ar, arp = _argon_species()
    imf = {electron: 1.0e-6, arp: 1.0e-6, ar: 1.0 - 2.0e-6}
    with pytest.raises(PlasmaStateError) as exc:
        PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf,
                      (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                      diffusion_length=(2.0, 's'),
                      ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'))
    assert 'length' in str(exc.value)


def test_direct_construction_checks_mobility_dimension():
    """MEDIUM: a diffusivity (m^2/s) is not a mobility; direct construction rejects it."""
    electron, ar, arp = _argon_species()
    imf = {electron: 1.0e-6, arp: 1.0e-6, ar: 1.0 - 2.0e-6}
    with pytest.raises(PlasmaStateError) as exc:
        PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf,
                      (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                      diffusion_length=(_diffusion_length(), 'm'),
                      ion_reduced_mobility=(1.0e-4, 'm^2/s'))
    assert 'mobility' in str(exc.value)


def test_wall_diagnostics_latched_only_at_accepted_states():
    """HIGH 5: the wall-flux interface is latched at accepted states (init, and after
    each accepted step) and NEVER from inside the residual, so no rejected Newton
    trial can leak in. A residual evaluation at a wild state moves the internal
    scratch (wall_loss_rates) but must leave the latched wall_flux untouched."""
    r, core = _metastable_reactor(gamma=1.0, neutralization={'Ar+': 'Ar'})
    latched = np.array(r.wall_flux, float)
    assert latched.shape[0] == r.num_core_species
    assert np.all(np.isfinite(latched))
    y = np.array(r.y0, float) * 1.0e6                 # a state the solver would reject
    r.residual(0.0, y, np.zeros_like(y))
    assert np.array_equal(np.array(r.wall_flux, float), latched)      # latch untouched
    assert not np.array_equal(np.array(r.wall_loss_rates, float), latched)  # scratch moved


def test_wall_energy_interface_declares_ion_term_absent_not_broken():
    """HIGH 5 + addition: a consumer must distinguish DECLARED-ABSENT (by design)
    from UNAVAILABLE (could not compute) from AVAILABLE, in code, without reading a
    docstring. The ion directed/sheath term is declared-absent (a sheath model is a
    contract non-goal); the electron thermal term is available from T_e; the
    neutralisation term is unavailable here because Ar+ carries no thermo."""
    r, core = _metastable_reactor(gamma=1.0, neutralization={'Ar+': 'Ar'})
    avail = r.wall_energy_availability
    assert avail['wall_ion_energy_flux'] == 'declared-absent'
    assert np.isnan(r.wall_ion_energy_flux)
    assert avail['wall_electron_energy_flux'] == 'available'
    assert np.isfinite(r.wall_electron_energy_flux) and r.wall_electron_energy_flux > 0.0
    assert avail['wall_neutralization_energy_flux'] == 'unavailable'
    assert np.isnan(r.wall_neutralization_energy_flux)


def test_saved_input_preserves_wall_neutralization_products():
    """HIGH 4: the input writer emits every other wall keyword but dropped
    wallNeutralizationProducts, so the declaration fallback did not survive a
    save/reload -- the reloaded deck could refuse, infer a different product, or undo
    a deliberate override. It must be written back."""
    from rmgpy.rmg.input import _format_plasma_wall
    reactor, core = _metastable_reactor(neutralization={'Ar+': 'Ar'})
    text = _format_plasma_wall(reactor)
    assert 'wallNeutralizationProducts' in text
    assert "'Ar+': 'Ar'" in text


def test_source_apportionment_jacobian_matches_fd_in_a_mixture():
    """HIGH 3 changed both the residual denominator and its Jacobian to the ionisable
    total; they must still agree to finite-difference precision. A He/Ar/Ar+ mixture
    with a source makes y_ionisable (Ar only) differ from the full neutral total
    (Ar+He), the path the argon-only analytic-Jacobian arm never exercises."""
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = _ground_species('Ar')
    he = Species(label='He').from_adjacency_list('1 He u0 p1 c0')
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    imf = {electron: 1.0e-7, arp: 1.0e-7, he: 0.4, ar: 0.6 - 2.0e-7}
    reactor = PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf,
                            (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                            diffusion_length=(_diffusion_length(), 'm'),
                            ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                            wall_recycling=1.0, ionisation_source=(1.0e5, 'm^-3/s'),
                            wall_single_bath_approximation=True)
    core = [electron, ar, he, arp]
    reactor.initialize_model(core, [], [], [])
    z = reactor.species_charges
    ie = reactor.electron_index
    i_ar = [j for j in range(len(z)) if z[j] == 0 and core[j].label == 'Ar'][0]
    i_he = [j for j in range(len(z)) if z[j] == 0 and core[j].label == 'He'][0]
    i_arp = [j for j in range(len(z)) if z[j] == 1 and j != ie][0]
    y = np.zeros(reactor.num_core_species, float)
    y[i_ar] = 0.6
    y[i_he] = 0.4
    y[i_arp] = 1.0e-7
    y[ie] = 1.0e-7
    best = _jacobian_scan(reactor, y)
    assert best < FD_TOLERANCE, best


# ============================================================ round-83 repairs
# A guard that reads a NEARBY quantity instead of the GOVERNING one. Presence where
# inventory governs (HIGH 1, HIGH 4); aggregate where per-species governs (MED); a
# projection where constitution governs (HIGH 2); relative where absolute governs
# (HIGH 3); residual-scratch on the advance() path where the live simulate()/step()
# path governs (the latch). Each test below reproduces one such site, red on the
# pre-round-83 module.


def test_step_publishes_diagnostics_on_the_production_path():
    """HIGH 1: ReactionSystem.simulate drives the reactor through step(), never
    advance(); the latch was wired only into advance(), so a production run published
    wall diagnostics frozen at initialisation -- a guarantee that never ran. After a
    real step() the published diagnostic time must track the solver time, not fall
    behind it, and the latched flux must reflect the advanced state."""
    r, _, _ = _build_reactor(wall=True, gamma=1.0, with_chemistry=False, x_ion=1.0e-6)
    t_init = r.wall_diagnostics_time             # latched at t = 0 by initialize_model
    assert t_init == 0.0
    flux0 = np.array(r.wall_flux, float)
    r.step(1.0e-6)                                # the entry simulate() actually uses
    assert r.t > 0.0
    # the check that would have caught it: published diagnostic time must not lag t
    assert r.wall_diagnostics_time == r.t
    assert r.wall_diagnostics_time > t_init
    assert not np.array_equal(np.array(r.wall_flux, float), flux0)


def test_wall_refuses_two_excited_states_with_no_ground():
    """HIGH 3: two neutral Ar states at 5 and 11.5 eV, the true ground (0 eV) absent.
    'Lowest present' (5 eV) is NOT the ground state, but nothing in the deck certifies
    that -- there is no absolute energy reference. Round 79 crowned the 5 eV state;
    round 83 refuses and names the declaration, because a clear ordering among PRESENT
    states cannot anchor the absolute ground."""
    with pytest.raises(PlasmaStateError) as exc:
        _metastable_reactor(ground_eV=5.0, meta_eV=11.5, meta_label='Ar_hi')
    msg = str(exc.value)
    assert 'ambiguous' in msg and 'heavy skeleton' in msg
    assert 'wallNeutralizationProducts' in msg


def test_skeleton_key_merges_tautomers_which_forces_the_declaration():
    """HIGH 2 root cause. Standard InChI deliberately merges tautomers, so the skeleton
    key yields the SAME value for 2-pyridone and 2-hydroxypyridine. A cation sharing
    that key matches BOTH neutrals, and an energy pick transmutes one tautomer into the
    other (round 79 did exactly this on pyridone+). The key cannot be refined to split
    them without breaking charge-independence: the FixedH /f layer that distinguishes
    tautomers sits BEHIND the /q,/p the key must truncate to unite an ion with its
    neutral (verified in evidence/round83_key_probe.log). So the fix is not a finer
    key -- that would be the fourth projection in a row -- but the multiplicity refusal:
    two matches, whatever they are, require a declaration. This pins the merge that
    makes that refusal necessary, and the isomer split that must survive it."""
    r, _, _ = _build_reactor(wall=True, gamma=0.0, with_chemistry=False)
    pyridone = Species(label='pyridone').from_smiles('O=c1cccc[nH]1')
    hydroxypyridine = Species(label='hydroxypyridine').from_smiles('Oc1ccccn1')
    assert r._skeleton_key(pyridone) == r._skeleton_key(hydroxypyridine)   # merged
    dme = Species(label='DME').from_smiles('COC')
    ethanol = Species(label='EtOH').from_smiles('CCO')
    assert r._skeleton_key(dme) != r._skeleton_key(ethanol)                # still split


def test_wall_refuses_electrons_without_commensurate_ion_inventory():
    """HIGH 4: the zero-cation guard checks a cation SPECIES exists, not that it is
    PRESENT in inventory. A declared Ar+ at zero moles with n_e > 0 passed, and the
    common wall loss frequency then removed electrons with no ion partner, driving the
    net charge. Quasineutrality is a precondition of this wall, knowable only at a
    state, so check_wall_support refuses a grossly non-neutral accepted state -- whether
    or not the electron is carried algebraically."""
    r, _, _ = _build_reactor(wall=True, gamma=0.0, with_chemistry=False)
    ie, i_ar, i_arp = _indices(r)
    y = np.zeros(r.num_core_species, float)
    y[i_ar] = 1.0
    y[i_arp] = 0.0                       # cation species present, zero inventory
    y[ie] = 1.0e-6                       # electrons with no ion partner
    with pytest.raises(PlasmaStateError) as exc:
        r.check_wall_support(y)
    assert 'net charge' in str(exc.value)


def test_declared_source_refused_when_no_ionisable_inventory_remains():
    """HIGH 5: the external pair source is apportioned over IONISABLE neutrals. When
    those are consumed but an unsupported bath gas (He) remains, y_ionisable -> 0 and
    the residual silently deposits ZERO -- the declared source vanishes with no error. A
    user who declares a source and receives none cannot detect it from the output, so an
    accepted state there is refused."""
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = _ground_species('Ar')
    he = Species(label='He').from_adjacency_list('1 He u0 p1 c0')
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    imf = {electron: 1.0e-7, arp: 1.0e-7, he: 0.5, ar: 0.5 - 2.0e-7}
    reactor = PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf,
                            (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                            diffusion_length=(_diffusion_length(), 'm'),
                            ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                            wall_recycling=1.0, ionisation_source=(1.0e18, 'm^-3/s'),
                            wall_single_bath_approximation=True)
    core = [electron, ar, he, arp]
    reactor.initialize_model(core, [], [], [])
    z = reactor.species_charges
    ie = reactor.electron_index
    i_ar = [j for j in range(len(z)) if z[j] == 0 and core[j].label == 'Ar'][0]
    i_he = [j for j in range(len(z)) if z[j] == 0 and core[j].label == 'He'][0]
    i_arp = [j for j in range(len(z)) if z[j] == 1 and j != ie][0]
    y = np.zeros(reactor.num_core_species, float)
    y[i_he] = 0.5                        # only unsupported bath gas remains
    y[i_ar] = 0.0                        # ionisable neutral consumed
    y[i_arp] = 1.0e-7
    y[ie] = 1.0e-7
    with pytest.raises(PlasmaStateError) as exc:
        reactor.check_wall_support(y)
    assert 'ionisable' in str(exc.value)


def test_check_wall_support_refuses_an_individual_negative_neutral():
    """MED: support validation summed the neutral inventory, so an individual negative
    population offset by a positive one passed the aggregate check -- and the residual
    then apportioned the wall source over per-species populations, turning the negative
    one into a negative (injecting) source allocation. The governing quantity is the
    per-species inventory, not the sum."""
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = _ground_species('Ar')
    he = Species(label='He').from_adjacency_list('1 He u0 p1 c0')
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    imf = {electron: 1.0e-7, arp: 1.0e-7, he: 0.5, ar: 0.5 - 2.0e-7}
    reactor = PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf,
                            (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                            diffusion_length=(_diffusion_length(), 'm'),
                            ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                            wall_single_bath_approximation=True)
    core = [electron, ar, he, arp]
    reactor.initialize_model(core, [], [], [])
    z = reactor.species_charges
    ie = reactor.electron_index
    i_ar = [j for j in range(len(z)) if z[j] == 0 and core[j].label == 'Ar'][0]
    i_he = [j for j in range(len(z)) if z[j] == 0 and core[j].label == 'He'][0]
    i_arp = [j for j in range(len(z)) if z[j] == 1 and j != ie][0]
    y = np.zeros(reactor.num_core_species, float)
    y[i_he] = 1.0
    y[i_ar] = -1.0e-3                    # individual negative; aggregate stays ~0.999
    y[i_arp] = 1.0e-7
    y[ie] = 1.0e-7
    with pytest.raises(PlasmaStateError) as exc:
        reactor.check_wall_support(y)
    assert 'negative' in str(exc.value)


def _noise_reactor():
    """Ar/He/Ar+/e- with a wall, and a state vector the caller perturbs one entry of."""
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = _ground_species('Ar')
    he = Species(label='He').from_adjacency_list('1 He u0 p1 c0')
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    imf = {electron: 1.0e-7, arp: 1.0e-7, he: 0.5, ar: 0.5 - 2.0e-7}
    reactor = PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf,
                            (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                            diffusion_length=(_diffusion_length(), 'm'),
                            ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                            wall_single_bath_approximation=True)
    core = [electron, ar, he, arp]
    reactor.initialize_model(core, [], [], [])
    z = reactor.species_charges
    ie = reactor.electron_index
    idx = {'e-': ie,
           'Ar': [j for j in range(len(z)) if z[j] == 0 and core[j].label == 'Ar'][0],
           'He': [j for j in range(len(z)) if z[j] == 0 and core[j].label == 'He'][0],
           'Ar+': [j for j in range(len(z)) if z[j] == 1 and j != ie][0]}
    y = np.zeros(reactor.num_core_species, float)
    y[idx['He']] = 0.5
    y[idx['Ar']] = 0.5
    y[idx['Ar+']] = 1.0e-7
    y[idx['e-']] = 1.0e-7
    return reactor, idx, y


def test_energy_off_wall_check_refuses_a_neutral_negative_within_atol():
    """I-274 round 3: the neutral-noise clamp belongs to the energy balance. With Te
    prescribed the wall check is exactly the pre-I-274 code: an individual negative neutral,
    however small, is refused and left as it was."""
    reactor, idx, y = _noise_reactor()
    assert not reactor.energy_balance
    atol = reactor.atol_array[idx['Ar']]
    y[idx['Ar']] = -0.5 * atol
    with pytest.raises(PlasmaStateError) as exc:
        reactor.check_wall_support(y)
    assert 'negative' in str(exc.value)
    assert y[idx['Ar']] == -0.5 * atol


def test_neutral_negative_beyond_its_atol_is_still_refused_and_left_unclamped():
    reactor, idx, y = _noise_reactor()
    atol = reactor.atol_array[idx['Ar']]
    y[idx['Ar']] = -2.0 * atol
    with pytest.raises(PlasmaStateError) as exc:
        reactor.check_wall_support(y)
    assert 'negative' in str(exc.value)
    assert y[idx['Ar']] == -2.0 * atol


@pytest.mark.parametrize('label', ['e-', 'Ar+'])
def test_charged_negative_within_atol_is_never_accepted(label):
    """The noise tolerance is for NEUTRALS only. A charged population below zero, however
    small, corrupts the charge bookkeeping the wall operator rests on and stays refused."""
    reactor, idx, y = _noise_reactor()
    atol = reactor.atol_array[idx[label]]
    y[idx['Ar+']] = 0.0
    y[idx['e-']] = 0.0
    y[idx[label]] = -0.5 * atol
    with pytest.raises(PlasmaStateError):
        reactor.check_wall_support(y)
    assert y[idx[label]] == -0.5 * atol


def test_non_finite_neutral_is_still_refused():
    reactor, idx, y = _noise_reactor()
    y[idx['Ar']] = float('nan')
    with pytest.raises(PlasmaStateError) as exc:
        reactor.check_wall_support(y)
    assert 'non-finite' in str(exc.value)


def test_direct_construction_checks_mobility_reference_density_dimension():
    """MED: a directly-constructed reactor must reject mobility_reference_density given
    in the wrong dimension. (3, 'kg') is a mass, not a number density; taking its SI
    value as m^-3 would scale the ion mobility from the wrong quantity."""
    with pytest.raises(PlasmaStateError) as exc:
        PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'),
                      {Species(label='e-').from_adjacency_list('1 e u1 p0 c-1'): 1.0},
                      (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                      diffusion_length=(_diffusion_length(), 'm'),
                      ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                      mobility_reference_density=(3.0, 'kg'))
    msg = str(exc.value)
    assert 'mobility_reference_density' in msg and 'number density' in msg


def test_direct_construction_checks_ionisation_source_dimension():
    """MED: ionisation_source given as (7, 'kg') is a mass, not a rate density; the
    constructor must reject it rather than take 7 as 7 m^-3 s^-1."""
    with pytest.raises(PlasmaStateError) as exc:
        PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'),
                      {Species(label='e-').from_adjacency_list('1 e u1 p0 c-1'): 1.0},
                      (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                      diffusion_length=(_diffusion_length(), 'm'),
                      ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                      ionisation_source=(7.0, 'kg'))
    assert 'ionisation_source' in str(exc.value)


def test_wall_only_options_without_a_wall_are_refused():
    """MED: wall-only options (here an ionisation source) supplied without the
    diffusion_length/mobility that DECLARE a wall were stored and silently ignored --
    the deck would run as a plain volume reactor with the source doing nothing. Refuse,
    naming the option."""
    with pytest.raises(PlasmaStateError) as exc:
        PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'),
                      {Species(label='e-').from_adjacency_list('1 e u1 p0 c-1'): 1.0},
                      (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                      ionisation_source=(1.0e18, 'm^-3/s'))
    msg = str(exc.value)
    assert 'no wall' in msg and 'ionisation_source' in msg


def test_wall_neutralization_products_unknown_ion_key_is_refused():
    """MED: a wallNeutralizationProducts key that names no cation in the core (a typo,
    'Ar+2' for 'Ar+') was silently ignored, so the declaration the user wrote to be safe
    did nothing and the ion fell back to inference. The unknown key is refused. (Single
    neutral Ar here, so Ar+ itself resolves cleanly -- the only fault is the typo.)"""
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = _ground_species('Ar')
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    imf = {electron: 1.0e-6, arp: 1.0e-6, ar: 1.0 - 2.0e-6}
    reactor = PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf,
                            (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                            diffusion_length=(_diffusion_length(), 'm'),
                            ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                            wall_neutralization_products={'Ar+2': 'Ar'})
    with pytest.raises(PlasmaStateError) as exc:
        reactor.initialize_model([electron, ar, arp], [], [], [])
    msg = str(exc.value)
    assert 'Ar+2' in msg and 'not a positive ion' in msg


def test_nonfinite_recycle_thermo_leaves_energy_unavailable():
    """Enthalpy left identity but still matters in the wall-energy interface: a NaN
    formation enthalpy no longer blocks initialisation (the declaration sets identity),
    and the neutralisation energy term is reported 'unavailable', not a silent finite
    number and not a crash."""
    reactor, core = _metastable_reactor(ground_eV=float('nan'), meta_eV=11.5,
                                        neutralization={'Ar+': 'Ar'})
    avail = reactor.wall_energy_availability['wall_neutralization_energy_flux']
    assert avail == 'unavailable'
    assert np.isnan(reactor.wall_neutralization_energy_flux)


# ================================ round 88 ================================
# Four HIGH + two MEDIUM. Three of the HIGH are the round-83 class -- a guard, a
# key, or a diagnostic reading a quantity adjacent to the one the physics governs.
# HIGH 2 is the wall made invisible to the termination path.


def _simulate(reactor, core, rxns, edge=None, edge_rxns=None):
    """Drive a reactor through the production simulate() entry, as the model builder
    does. Settings are the inert-friendly ones the steady-state suite uses."""
    return reactor.simulate(
        core, rxns, edge or [], edge_rxns or [], [], [],
        model_settings=ModelSettings(tol_keep_in_edge=0, tol_move_to_core=1e5,
                                     tol_interrupt_simulation=1e8),
        simulator_settings=SimulatorSettings())


# ---- HIGH 2: a wall-only system must not declare itself inert on simulate() ----

def test_wall_only_deck_integrates_past_t0_on_the_production_path():
    """HIGH 2: a wall-only deck (no gas-phase chemistry, gamma=0) carries no
    char_rate, but the wall is removing ~1.9e-2 mol/s of ion-electron pairs. The
    inert test read char_rate -- the gas-phase CHEMISTRY diagnostic -- as the
    system's total flux and terminated at t=0 claiming 'the composition cannot
    change'. Demonstrated on simulate(), the production entry, not a residual probe.
    The inert test must consult the reactor's TOTAL flux; the wall term moves it."""
    term = [TerminationSteadyState(tolerance=1e-8), TerminationTime((1.0, 's'))]
    r, core, rxns = _build_reactor(wall=True, gamma=0.0, with_chemistry=False,
                                   x_ion=1.0e-4, termination=term)
    y0 = np.array(r.y0[:r.num_core_species], float)
    terminated, _res, _inv, _ss, _sr, t_final, _conv = _simulate(r, core, rxns)
    # The wall is depleting the plasma, so the run must integrate, not stop at t=0.
    assert t_final > 0.0
    assert not np.array_equal(np.array(r.y[:r.num_core_species], float), y0)


def test_genuinely_inert_deck_still_terminates_at_t0_unchanged():
    """The other half of HIGH 2's ruling: a reactor with no non-chemical flux must
    behave EXACTLY as before. A wall-less plasma deck with no chemistry has zero
    total flux and must still be recognised as inert and terminate at t=0."""
    term = [TerminationSteadyState(tolerance=1e-8), TerminationTime((1.0, 's'))]
    r, core, rxns = _build_reactor(wall=False, with_chemistry=False, termination=term)
    terminated, _res, _inv, _ss, _sr, t_final, _conv = _simulate(r, core, rxns)
    assert terminated and t_final == 0.0


def test_adapter_does_not_launder_a_source_channels_nonfinite_reading_to_absence():
    """Round 113 BLOCKING 1, THROUGH the production path (the base.pyx adapter), not by calling
    update() directly. Round 112 closed the reproductions inside update(), but base.pyx mapped
    every UNARMED non-finite hook value back to None -- 'no channel' -- so a source-driven
    electron channel that reported nan/+inf/-inf while unarmed was laundered to absence and a flat
    generic channel could fire on it. A reactor with an active ionisation source has a channel IN
    PLAY every step; a non-finite reading it produces must poison, armed or not. Driven through the
    reactor's own simulate(): the adapter must hand the poison to the criterion (never None), and
    the run must not reach steady state on it. (A wall deck with NO source keeps reporting nan as
    its documented 'no such channel' sentinel, and is still judged on the generic channel -- pinned
    by ``test_wall_only_deck_integrates_past_t0_on_the_production_path``.)"""
    orig = TerminationSteadyState.update
    for poison in (float('nan'), float('inf'), float('-inf')):
        class _PoisonChannel(PlasmaReactor):
            def steady_state_external_residual(self, t_now, y_now, t_prev, y_prev):
                return poison
            def steady_state_external_armed(self, t_now, y_now):
                return False

        term = [TerminationSteadyState(tolerance=1e-8), TerminationTime((50.0, 's'))]
        r, core, rxns = _build_reactor(wall=True, gamma=0.5, with_chemistry=False,
                                       x_ion=1.0e-4, source=1.0e22, termination=term,
                                       reactor_cls=_PoisonChannel)
        seen = []

        def spy(self, *a, **k):
            ret = orig(self, *a, **k)
            seen.append((k.get('external_residual'), ret))
            return ret

        TerminationSteadyState.update = spy
        try:
            _simulate(r, core, rxns)
        finally:
            TerminationSteadyState.update = orig

        supplied = [v for v, ret in seen if v is not None]
        assert supplied, (poison, "the adapter laundered every source-channel reading to None")
        assert all(not np.isfinite(v) for v in supplied), (poison, supplied[:3])
        assert not r.steady_state_reached, (
            poison, "a non-finite source-channel reading was laundered into a steady state")



# ================================ round 114 ================================
# BLOCKING: the adapter decided whether an external channel is present by reading the
# PlasmaReactor's own ``has_wall and ionisation_source > 0``, so ANY other reactor that
# overrides the documented hooks had its reading silently discarded. Presence is now
# declared by the hook contract (``steady_state_external_channel``: None = no channel,
# a float = a supplied reading), one rule at the adapter for every reactor.


def _wallless_hook_run(value, tolerance, hook='legacy'):
    """A WALL-LESS PlasmaReactor subclass whose steady-state hook returns ``value``, driven
    through simulate(). ``hook='legacy'`` overrides the documented double hook
    ``steady_state_external_residual``; ``'channel'`` overrides the object hook."""
    if hook == 'legacy':
        class _Hooked(PlasmaReactor):
            def steady_state_external_residual(self, t_now, y_now, t_prev, y_prev):
                return value
    else:
        class _Hooked(PlasmaReactor):
            def steady_state_external_channel(self, t_now, y_now, t_prev, y_prev):
                return value
    term = [TerminationSteadyState(tolerance=tolerance), TerminationTime((50.0, 's'))]
    r, core, rxns = _build_reactor(wall=False, termination=term, reactor_cls=_Hooked)
    terminated, _res, _inv, _ss, _sr, t_final, _conv = _simulate(r, core, rxns)
    return r, terminated, t_final


def test_wallless_control_reaches_steady_state_without_a_hook():
    """The control the round-114 tests stand on: the same wall-less deck, no hook override,
    DOES terminate as steady before the 50 s backstop. Without it, 'not steady' below could
    mean the deck simply never settles rather than that the hook's reading was honoured."""
    term = [TerminationSteadyState(tolerance=1e-6), TerminationTime((50.0, 's'))]
    r, core, rxns = _build_reactor(wall=False, termination=term)
    _t, _res, _inv, _ss, _sr, t_final, _conv = _simulate(r, core, rxns)
    assert r.steady_state_reached and t_final < 50.0


@pytest.mark.parametrize('poison', [float('inf'), float('nan'), float('-inf')])
def test_wallless_subclass_nonfinite_legacy_hook_poisons_through_simulate(poison):
    """Round 114 BLOCKING, the reported reproduction: a wall-less subclass whose documented
    hook returns a non-finite value terminated as steady (1.4049 s, residual 3.97e-11) because
    the adapter consulted ``has_wall`` and discarded the reading. A reactor that overrides the
    hook supplies a channel; a non-finite supplied reading must poison the decision."""
    r, _terminated, t_final = _wallless_hook_run(poison, 1e-6)
    assert not r.steady_state_reached, (poison, t_final, r.steady_state_residual)


def test_wallless_subclass_finite_legacy_hook_is_honoured_through_simulate():
    """Round 114 BLOCKING, the finite half (a regression against a0de2256d): a supplied finite
    residual 2.0 is far above tolerance 1e-6, so it must hold the run open -- not be discarded."""
    r, _terminated, t_final = _wallless_hook_run(2.0, 1e-6)
    assert not r.steady_state_reached, (t_final, r.steady_state_residual)


@pytest.mark.parametrize('poison', [float('inf'), float('nan'), float('-inf')])
def test_channel_hook_nonfinite_poisons_through_simulate(poison):
    """The object hook itself: a supplied non-finite reading poisons, for any reactor."""
    r, _terminated, t_final = _wallless_hook_run(poison, 1e-6, hook='channel')
    assert not r.steady_state_reached, (poison, t_final)


def test_channel_hook_none_is_absence_through_simulate():
    """None is the contract's 'no channel': the generic criterion decides alone, exactly as
    the no-hook control does."""
    r, _terminated, t_final = _wallless_hook_run(None, 1e-6, hook='channel')
    assert r.steady_state_reached and t_final < 50.0


def test_channel_hook_default_declares_absence_for_an_ordinary_and_a_sourceless_plasma():
    """The default contract: a reactor that overrides nothing has no channel (None), and the
    frozen PlasmaReactor implementation keeps its documented absence -- no source-driven wall
    -- so the plasma deck path is unchanged. With a source it supplies its reading every step
    (nan included), which is what round 113 pinned."""
    r, _core, _rxns = _build_reactor(wall=False)
    y = np.array(r.y0[:r.num_core_species], float)
    assert r.steady_state_external_channel(2.0, y, 1.0, y) is None
    r, _core, _rxns = _build_reactor(wall=True, gamma=0.5, source=None)
    assert r.steady_state_external_channel(2.0, y, 1.0, y) is None
    r, _core, _rxns = _build_reactor(wall=True, gamma=0.5, source=1.0e22)
    reading = r.steady_state_external_channel(2.0, y, 1.0, y)
    assert reading is not None
    assert reading == r.steady_state_external_residual(2.0, y, 1.0, y) or (
        np.isnan(reading) and np.isnan(r.steady_state_external_residual(2.0, y, 1.0, y)))


# ---- HIGH 1: the skeleton key must key charged and neutral species alike ----

def _dme_isotope_species():
    e = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    dme12 = Species(label='DME').from_smiles('COC')
    dme12.thermo = _thermo_with_h298(-184.0)
    dme13 = Species(label='DME-13C').from_smiles('[13CH3]O[CH3]')
    dme13.thermo = _thermo_with_h298(-184.0)
    dmep13 = Species(label='DME+-13C').from_smiles('[13CH3][O+][CH3]')
    return e, dme12, dme13, dmep13


def test_skeleton_key_keeps_isotopes_so_the_wall_cannot_transmute_nuclei():
    """HIGH 1: the key truncated the standard InChI at the first /q or /p, but the
    isotope layer /i sits AFTER /q, so a charged species lost it and a neutral did
    not -- the key was a different rule depending on charge. A 13C cation then keyed
    identically to an ordinary-carbon neutral, and the wall recycled 13C into 12C: a
    transmuted nucleus. Removing the /q and /p layers (not truncating at them) keeps
    /i, so isotopes discriminate on both charge states while Ar and Ar* still
    coincide."""
    r, _, _ = _build_reactor(wall=True, with_chemistry=False)
    _e, dme12, dme13, dmep13 = _dme_isotope_species()
    # charge independence preserved: the 13C neutral and 13C cation still coincide
    assert r._skeleton_key(dme13) == r._skeleton_key(dmep13)
    # the nucleus discriminates: the 13C cation must NOT key as the 12C neutral
    assert r._skeleton_key(dmep13) != r._skeleton_key(dme12)


def test_declared_isotopic_neutral_is_not_hidden_by_key_truncation():
    """HIGH 1, the escape hatch: a correct declaration naming the isotopic neutral
    must resolve. Truncation stripped the neutral's /i layer from the ion's key, so
    the isotopic neutral never entered the candidate 'matches' list and the
    declaration was refused as 'not sharing the heavy composition'. With the layer
    kept, ion and its isotopic neutral share a key and the declaration resolves."""
    e, _dme12, dme13, dmep13 = _dme_isotope_species()
    imf = {e: 1.0e-6, dmep13: 1.0e-6, dme13: 1.0 - 2.0e-6}
    reactor = PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf,
                            (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                            diffusion_length=(_diffusion_length(), 'm'),
                            ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                            wall_recycling=1.0,
                            wall_neutralization_products={'DME+-13C': 'DME-13C'})
    core = [e, dme13, dmep13]
    reactor.initialize_model(core, [], [], [])
    i_cat = _cation_index(reactor)
    assert core[int(reactor.wall_recycle_target[i_cat])].label == 'DME-13C'


# ---- HIGH 3: quasineutrality is a ratio, not an absolute mole floor ----

def test_quasineutrality_bound_is_relative_not_an_absolute_mole_floor():
    """HIGH 3: quasineutrality is a RATIO, but the guard compared the net charge to a
    fixed 1e-12 mol. An electron inventory of 1e-13 mol with no ion partner at all is
    100% non-neutral, yet slipped under the absolute floor. The bound must be a
    fraction of the charged inventory, independent of the system's size and units."""
    e = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = _ground_species('Ar')
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    imf = {e: 1.0e-13, arp: 0.0, ar: 1.0 - 1.0e-13}
    reactor = PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf,
                            (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                            diffusion_length=(_diffusion_length(), 'm'),
                            ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                            wall_recycling=0.0)
    with pytest.raises(PlasmaStateError) as exc:
        reactor.initialize_model([e, ar, arp], [], [], [])
    assert 'quasineutrality' in str(exc.value) or 'net charge' in str(exc.value)


def test_initial_quasineutrality_bound_is_relative_under_algebraic_electron():
    """HIGH 3 at the second site: with the electron carried algebraically, the packed
    initial state must SATISFY quasineutrality, and 'satisfy' is again relative. Ar+
    at 5e-13 mol against e- at 1e-13 mol is a net +4e-13 -- a two-thirds charge
    imbalance -- that the absolute 1e-12 mol floor admitted."""
    e = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = _ground_species('Ar')
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    imf = {e: 1.0e-13, arp: 5.0e-13, ar: 1.0 - 6.0e-13}
    reactor = PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf,
                            (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                            quasineutral_electron=True)
    with pytest.raises(PlasmaStateError) as exc:
        reactor.initialize_model([e, ar, arp], [], [], [])
    assert 'net charge' in str(exc.value)


# ---- HIGH 4: neutralisation energy is owed per ion lost, not scaled by gamma ----

def _energy_reactor(gamma):
    e = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = _ground_species('Ar', excitation_eV=0.0)
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    arp.thermo = _argon_thermo(15.76)          # ionisation energy as the enthalpy offset
    imf = {e: 1.0e-6, arp: 1.0e-6, ar: 1.0 - 2.0e-6}
    reactor = PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf,
                            (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                            diffusion_length=(_diffusion_length(), 'm'),
                            ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                            wall_recycling=gamma)
    reactor.initialize_model([e, ar, arp], [], [], [])
    return reactor


def test_neutralisation_energy_is_owed_per_ion_lost_not_scaled_by_recycling():
    """HIGH 4: wall_recycling is a MASS-return fraction. The neutralisation enthalpy
    is deposited when the ion recombines with an electron at the wall, which happens
    whether or not the neutral returns to the gas. The energy term was gated on and
    multiplied by gamma, so a fully-pumping wall (gamma=0) reported 0 W -- and, worse,
    labelled it 'available'. The energy must be independent of gamma."""
    r0 = _energy_reactor(gamma=0.0)
    r1 = _energy_reactor(gamma=1.0)
    e0 = r0.wall_neutralization_energy_flux
    e1 = r1.wall_neutralization_energy_flux
    assert r0.wall_energy_availability['wall_neutralization_energy_flux'] == 'available'
    assert e0 > 0.0
    # gamma governs where the neutral goes, not whether the enthalpy is deposited
    assert np.isclose(e0, e1, rtol=1e-12, atol=0.0)


# ---- MEDIUM: a duplicate declaration label, and a truthy-string flag ----

def test_duplicate_neutral_label_in_declaration_is_refused():
    """MEDIUM: a declaration names ONE product, but if two core species carry that
    label the name identifies two things and core ordering silently picked the first.
    Round 83's whole point is that the modeller names the product; an ambiguous label
    defeats it. Refuse."""
    with pytest.raises(PlasmaStateError) as exc:
        _metastable_reactor(meta_label='Ar', neutralization={'Ar+': 'Ar'})
    msg = str(exc.value)
    assert "'Ar'" in msg and ('more than one' in msg or 'ambiguous' in msg)


def test_quasineutral_electron_flag_parses_boolean_strings_strictly():
    """MEDIUM: quasineutral_electron went through bool(value), so the STRING 'False'
    -- any non-empty string -- enabled quasineutral mode. Parse boolean-like strings
    by value, and refuse a string that is not boolean-like rather than reading it as
    True."""
    e, ar, arp = _argon_species()
    imf = {e: 1.0e-6, arp: 1.0e-6, ar: 1.0 - 2.0e-6}
    r = PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf,
                      (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                      quasineutral_electron='False')
    assert r.quasineutral_electron is False
    with pytest.raises(PlasmaStateError):
        PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf,
                      (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                      quasineutral_electron='maybe')


# ====================================================================== round 90
#
# Two HIGH (external-source ignition; the neutral floor's dimension), two MEDIUM
# (availability honesty; strict boolean coercion) and a LOW (a char_rate=0 divide),
# each reproduced red-first on the built module in evidence/round90_before.log.


def test_external_source_ignites_a_zero_electron_deck_through_simulate():
    """HIGH 1: a discharge driven by an external ionisation_source must start from
    neutral gas -- exactly zero electrons. The source is a zeroth-order electron-ion pair
    source (its rate does not depend on n_e), so it seeds the first electrons and the
    integrated electron population rises from 0. Demonstrated through the production
    simulate() entry, not a constructor probe: before round 90 the reactor refused to
    initialise at all."""
    r, core, rxns = _build_reactor(wall=True, with_chemistry=False, x_ion=0.0, source=1.0e20)
    ie = r.electron_index
    assert r.y0[ie] == 0.0, "the deck must start at exactly zero electrons"
    r.termination = [TerminationTime((1.0e-4, 's'))]
    _simulate(r, core, rxns)
    assert r.y[ie] > 0.0, "the ionisation source did not seed any electrons"


def test_zero_electron_deck_without_a_source_is_still_refused():
    """HIGH 1 boundary: with no ionisation_source the only electron production is the
    n_e-proportional gas-phase chemistry, so a zero-electron state is a fixed point that
    can never ignite. The guard must still refuse it, by name."""
    with pytest.raises(PlasmaStateError) as exc:
        _build_reactor(wall=True, with_chemistry=False, x_ion=0.0, source=None)
    msg = str(exc.value)
    assert 'ignite' in msg and 'ionisation_source' in msg


def test_neutral_floor_is_a_density_not_a_history_dependent_mole_count():
    """HIGH 2: nu_wall depends on the neutral number DENSITY (intensive), so acceptance
    must not depend on the neutral MOLES (the deck's absolute inventory). A state with
    small neutral moles but normal density -- which an isobaric reactor holds at ~P/kT no
    matter how much mass the wall has pumped out -- must be accepted and its nu_wall read
    from the density, unclamped. The old extensive 1e-6-mol floor refused it and clamped
    nu_wall on the reactor's history (nu 18.5 vs 185 for the same intensive state)."""
    r, _, _ = _build_reactor(wall=True, with_chemistry=False, x_ion=1.0e-9)   # ~1 mol init
    ie, i_ar, i_arp = _indices(r)
    y = np.zeros(r.num_core_species, float)
    y[i_ar] = 1.0e-7          # neutral MOLES far below the old 1e-6-mol floor
    y[i_arp] = 1.0e-13
    y[ie] = 1.0e-13           # alpha ~ 1e-6, far under the ceiling
    V = r.compute_volume(y)
    n_neutral = _neutral_density(r, y, V)
    assert n_neutral > r.wall_neutral_density_floor, "state should be at normal density"
    r.check_wall_support(y)                    # accepted, not refused
    nu = r.compute_nu_wall(y, V)
    closed = _nu_wall_closed_form(TE_NOMINAL_EV, n_neutral, r.diffusion_length.value_si)
    assert abs(nu / closed - 1.0) < 1.0e-12, "nu_wall was clamped on an in-domain density"


def test_non_finite_wall_flux_is_not_reported_available():
    """MEDIUM: an availability flag must be checked at least as hard as the number it
    vouches for -- a non-finite wall_flux must read 'unavailable', not assert 'available'
    over an inf. (Round 93's MEDIUM 1 now refuses the extreme mobility this once used at
    CONSTRUCTION, so a non-finite wall_flux is instead produced by latching a hand-built
    degenerate state -- a non-finite electron population -- which the accepted-step domain
    check would reject but the latch must still describe honestly.)"""
    r, _, _ = _build_reactor(wall=True, with_chemistry=False)
    y = np.array(r.y0[:r.num_core_species], float)
    V = r.compute_volume(y)
    y[r.electron_index] = float('inf')
    r._latch_wall_diagnostics(y, V, 0.0)
    assert not np.isfinite(np.asarray(r.wall_flux)).all()
    assert r.wall_energy_availability['wall_flux'] == 'unavailable'


def test_sub_underflow_diffusion_length_is_refused_by_name():
    """MEDIUM: a finite, positive diffusion length whose SQUARE underflows to zero would
    divide by zero in nu_wall = D_a/Lambda^2 and surface as a raw ZeroDivisionError from
    the Fortran callback. Refuse it by name at construction instead."""
    with pytest.raises(PlasmaStateError) as exc:
        _build_reactor(wall=True, with_chemistry=False, lam=1.0e-200)
    msg = str(exc.value)
    assert 'diffusion_length' in msg and 'square' in msg


def test_coerce_bool_flag_refuses_uninterpretable_values():
    """MEDIUM: _coerce_bool_flag used to end in bool(value), so 2, 0.5, NaN and object()
    enabled the flag while [] and {} disabled it -- the truthiness of the type, not the
    meaning. Refuse anything that is not a bool, None, or a recognised boolean string;
    the interpretable cases still pass through."""
    from rmgpy.solver.plasma import _coerce_bool_flag
    for v in (2, 0.5, float('nan'), object(), [], {}):
        with pytest.raises(PlasmaStateError):
            _coerce_bool_flag(v, 'quasineutral_electron', 'id')
    assert _coerce_bool_flag(True, 'q', 'id') is True
    assert _coerce_bool_flag(None, 'q', 'id') is False
    assert _coerce_bool_flag('on', 'q', 'id') is True


def test_wall_only_run_does_not_divide_by_zero_chemistry_rate():
    """LOW: base.pyx builds core/edge/network rate RATIOS by dividing by the CHEMISTRY
    char_rate. In a supported wall-only run every core rate is exactly 0, so char_rate is
    0 and the ratios were 0/0 NaN (an 'invalid value encountered in divide' warning). The
    ratios are chemistry-relative enlargement signals and keep reading char_rate -- only
    the inert / termination gates read the total -- but the 0/0 must not launder a NaN.
    Under strict float-error handling the wall-only run must not raise."""
    r, core, rxns = _build_reactor(wall=True, with_chemistry=False, x_ion=1.0e-6)
    r.termination = [TerminationTime((1.0e-4, 's'))]
    # Premise: with no gas chemistry the characteristic chemistry rate is exactly zero, so
    # every core/edge/network/surface rate ratio divides by zero.
    y0 = np.array(r.y0, float)
    r.residual(0.0, y0.copy(), np.zeros_like(y0))
    assert np.all(np.asarray(r.core_species_rates, float) == 0.0), "the deck must carry no gas chemistry"
    # Under strict float-error handling the wall-only run must complete without a 0/0 divide
    # (the ratios abstain rather than laundering a NaN through argmax and the branching numbers).
    with np.errstate(invalid='raise', divide='raise'):
        result = _simulate(r, core, rxns)
    assert result[0] is True, "the wall-only run did not reach its termination-time backstop"


def test_neutral_density_floor_is_independent_of_initial_inventory():
    """HIGH 2, the two-histories statement made structural: the floor no longer reads the
    initial inventory at all, so two reactors whose initial neutral moles differ get the
    IDENTICAL density floor and identical nu_wall at the same intensive state. The old
    1e-6*(initial neutral moles) floor differed with the inventory, which is exactly how
    two histories reaching the same intensive state were judged differently."""
    r1, _, _ = _build_reactor(wall=True, with_chemistry=False, x_ion=1.0e-9)
    r2, _, _ = _build_reactor(wall=True, with_chemistry=False, x_ion=1.0e-4)
    ie1, i_ar1, i_arp1 = _indices(r1)
    ie2, i_ar2, i_arp2 = _indices(r2)
    assert r1.y0[i_ar1] != r2.y0[i_ar2], "the two decks must differ in initial inventory"
    # ...yet the intensive floor is identical
    assert r1.wall_neutral_density_floor == r2.wall_neutral_density_floor
    # ...and the same intensive state yields identical nu_wall from both
    y = np.zeros(r1.num_core_species, float)
    y[i_ar1] = 1.0e-7
    y[i_arp1] = 1.0e-13
    y[ie1] = 1.0e-13
    assert r1.compute_nu_wall(y, r1.compute_volume(y)) == r2.compute_nu_wall(y, r2.compute_volume(y))


# ====================================================================== round 93
#
# Ignition works (round 90) -- but it is a new dynamical regime, and two mechanisms
# correct for a SEEDED run are wrong at the zero boundary. Two HIGH (algebraic-mode
# ignition; a reached steady state the criterion cannot recognise), three MEDIUM
# (finite inputs -> infinite nu_wall; a pure-parent mobility against a gas mixture;
# a floor that moves under a physics-preserving reparameterisation) and a LOW (a
# nonzero source flux lost to a squaring underflow). Each reproduced red-first on the
# built module in evidence/round93_before.log.


def test_external_source_ignites_a_zero_electron_deck_in_quasineutral_mode():
    """HIGH 1: round 90 proved ignition from exactly zero electrons in the INTEGRATED
    mode; the algebraic quasineutralElectron=True mode still refused. The algebraic
    charge row is driven to the solver's ABSOLUTE accuracy, but check_wall_support tests
    quasineutrality RELATIVELY. At the first microstep the entire charged inventory sits
    ~2e-33 mol -- seventeen orders below the integrator's atol -- so the row's absolute
    convergence noise reads as an ~11% relative imbalance and the guard refuses a state
    the solver has merely not yet resolved. The relative guard must stand down while the
    whole charged inventory is below the integrator's resolution; above it the row holds
    net/magnitude at machine epsilon on its own. Demonstrated across the source range."""
    for src in (1.0e5, 1.0e12, 1.0e20):
        r, core, rxns = _build_reactor(wall=True, with_chemistry=False, x_ion=0.0,
                                       source=src, quasineutral=True)
        assert r.y0[r.electron_index] == 0.0, "the deck must start at exactly zero electrons"
        r.termination = [TerminationTime((1.0e-4, 's'))]
        _simulate(r, core, rxns)   # must not raise the quasineutrality refusal
        assert r.y[r.electron_index] > 0.0, (
            "source={0:g} did not ignite the algebraic-electron deck".format(src))


def test_saturating_discharge_reports_a_steady_state_and_a_dead_one_does_not():
    """HIGH 2: a source-driven discharge reaches n_e = S/nu_wall and holds it, yet the
    criterion never armed. A saturating-from-zero trajectory has log-log slope bounded by
    1 -- ``nu*t/(exp(nu*t)-1)`` -- so R>=1 never fires; and the electron saturates far
    below the mole floor, so the generic residual reads only the inert neutrals (which a
    weak discharge never perturbs, residual exactly 0 from t=0). The reactor supplies the
    electron's own slope so firing waits for it to saturate, and arms on ``t*nu_wall>=1``
    (the R=t/tau>=1 standard evaluated from the known relaxation time). Same criterion, a
    model that never started -- no source -- must still report NOT a steady state."""
    term = [TerminationSteadyState(tolerance=1.0e-8), TerminationTime((200.0, 's'))]
    r, core, rxns = _build_reactor(wall=True, with_chemistry=False, x_ion=0.0,
                                   source=1.0e5, termination=term)
    _simulate(r, core, rxns)
    yv = np.asarray(r.y[:r.num_core_species], float)
    V = r.compute_volume(yv)
    nu = r.compute_nu_wall(yv, V)
    n_e_density = r.y[r.electron_index] * constants.Na / V
    assert r.steady_state_reached, "the saturated discharge was not recognised as steady"
    assert abs(n_e_density / (1.0e5 / nu) - 1.0) < 1.0e-3, (
        "settled electron density {0!r} m^-3 is not S/nu_wall {1!r}".format(
            n_e_density, 1.0e5 / nu))

    term2 = [TerminationSteadyState(tolerance=1.0e-8), TerminationTime((200.0, 's'))]
    r2, core2, rxns2 = _build_reactor(wall=False, with_chemistry=False, termination=term2)
    _simulate(r2, core2, rxns2)
    assert not r2.steady_state_reached, (
        "a model that never started must not report a steady-state result")


def test_finite_wall_inputs_that_make_nu_wall_infinite_are_refused_at_construction():
    """MEDIUM 1: a finite reduced mobility (1e308) and a finite diffusion length whose
    SQUARE is subnormal (Lambda=1e-160 -> Lambda^2=1e-320) both give nu_wall=inf, yet were
    admitted -- the run then carries a non-finite residual. Marking the wall term
    'unavailable' is not the same as refusing a state that cannot be integrated. Refuse at
    construction: nu_wall evaluated at the reference density must be finite, and Lambda^2
    must be a normal (not subnormal) divisor."""
    with pytest.raises(PlasmaStateError):
        _build_reactor(wall=True, with_chemistry=False, mu0=1.0e308)
    with pytest.raises(PlasmaStateError):
        _build_reactor(wall=True, with_chemistry=False, lam=1.0e-160)


def test_neutral_mixture_is_refused_unless_the_single_bath_approximation_is_opted_into(caplog):
    """MEDIUM 2 (round 100): the wall carries ONE ion reduced mobility (Ar+ in Ar) but
    n_neutral sums every neutral heavy species, so a 50/50 Ar/He mixture applies the
    Ar+-in-Ar mobility to the He fraction too -- transport that does not describe the gas. A
    warning is not enough (a consumer of the latched fluxes cannot read a log line) and a
    passive availability label is not enough (it documents the wrong transport rather than
    gating on it). So the mixture is REFUSED at construction unless the user opts in via
    wall_single_bath_approximation=True, consciously accepting the approximation. Refusing
    outright would forbid every multi-species plasma, and running silently is wrong; the
    opt-in is the honest middle. With the opt-in the run proceeds and WARNS, naming the
    gases; input.rst says the same at the mobility keyword. Ar and its metastable Ar* share
    one skeleton (one bath, exact) and neither refuse nor warn (separate test)."""
    e, ar, arp = _argon_species()
    he = Species(label='He').from_adjacency_list('1 He u0 p1 c0')
    he.thermo = _thermo_with_h298(0.0)
    imf = {e: 1.0e-6, arp: 1.0e-6, ar: 0.5 - 1.0e-6, he: 0.5 - 1.0e-6}

    # Without the opt-in: refused at construction, naming the gases and the keyword.
    with pytest.raises(PlasmaStateError) as exc:
        r = PlasmaReactor(
            (TGAS, 'K'), (P_NOMINAL, 'Pa'), imf, (TE_NOMINAL_EV * EV_TO_K, 'K'),
            n_sims=1, termination=[],
            diffusion_length=(_diffusion_length(), 'm'),
            ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'), wall_recycling=1.0)
        r.initialize_model([e, ar, he, arp], [], [], [])
    assert 'spans more than one gas' in str(exc.value)
    assert 'wall_single_bath_approximation' in str(exc.value)

    # With the opt-in: constructs and warns.
    with caplog.at_level(logging.WARNING):
        reactor = PlasmaReactor(
            (TGAS, 'K'), (P_NOMINAL, 'Pa'), imf, (TE_NOMINAL_EV * EV_TO_K, 'K'),
            n_sims=1, termination=[],
            diffusion_length=(_diffusion_length(), 'm'),
            ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'), wall_recycling=1.0,
            wall_single_bath_approximation=True)
        reactor.initialize_model([e, ar, he, arp], [], [], [])   # accepted on opt-in
    assert 'spans more than one gas' in caplog.text, \
        "an opted-in bath mixture must still warn that the single mobility is an approximation"

    # the single-bath Ar/Ar* deliverable shares one skeleton and must NOT warn
    caplog.clear()
    with caplog.at_level(logging.WARNING):
        _metastable_reactor(gamma=1.0, neutralization={'Ar+': 'Ar'})
    assert 'spans more than one gas' not in caplog.text, \
        "Ar and Ar* are one bath and must not warn"


def test_density_floor_is_invariant_under_mobility_reparameterisation():
    """MEDIUM 3: the transport law reads mu0*Nref as a product, so scaling
    (Nref -> Nref*c, mu0 -> mu0/c) leaves nu_wall bit-identical -- but the acceptance floor
    was FRACTION*Nref, which scaled by c. Physically identical inputs were then accepted or
    refused differently. This is round 90's HIGH 2 in a new coordinate: the floor stopped
    depending on history and now depended on parameterisation. It must be built from a
    parameterisation-invariant density."""
    def floor_of(c):
        e, ar, arp = _argon_species()
        imf = {e: 1.0e-6, arp: 1.0e-6, ar: 1.0 - 2.0e-6}
        r = PlasmaReactor(
            (TGAS, 'K'), (P_NOMINAL, 'Pa'), imf, (TE_NOMINAL_EV * EV_TO_K, 'K'),
            n_sims=1, termination=[],
            diffusion_length=(_diffusion_length(), 'm'),
            ion_reduced_mobility=(MU0_AR_IN_AR / c, 'm^2/(V*s)'),
            mobility_reference_density=(PLASMA_LOSCHMIDT * c, 'm^-3'), wall_recycling=1.0)
        r.initialize_model([e, ar, arp], [], [], [])
        return r.wall_neutral_density_floor
    base = floor_of(1.0)
    for c in (1.0e3, 1.0e-3):
        assert floor_of(c) == base, (
            "the density floor moved under a physics-preserving reparameterisation (c={0:g})".format(c))


def test_small_but_nonzero_source_flux_is_not_lost_to_a_squaring_underflow():
    """LOW: get_non_chemical_char_rate squared the delivered per-species rates to form an
    L2 norm, so a small source whose flux is representable but whose SQUARE underflows
    (below ~sqrt(DBL_MIN)) reported exactly zero -- reading as an inert reactor while a
    source was declared and admitted. The norm must be scale-robust so a nonzero flux
    stays nonzero.

    The scalar norm alone is a weak witness -- it is positive if EITHER pair member is
    injected. Round 100 strengthens it: a pair source seeds an ion AND an electron, so
    check the per-species residual delivers a positive rate into BOTH the electron and the
    cation, not merely that some aggregate is nonzero."""
    r, core, rxns = _build_reactor(wall=True, with_chemistry=False, x_ion=0.0, source=1.0e-140)
    assert r.get_non_chemical_char_rate() > 0.0
    ie, _i_ar, i_arp = _indices(r)
    y = np.array(r.y0[:r.num_core_species], float)
    delta, _ = r.residual(0.0, y, np.zeros(r.num_core_species, float))
    assert delta[ie] > 0.0, "the source did not inject the electron of the pair"
    assert delta[i_arp] > 0.0, "the source did not inject the cation of the pair"


# ---------------------------------------------------------------- round 96
#
# The steady-state REPORT was unsound in both directions, the round-93 MEDIUMs were
# sharpened, and a subnormal source slipped the ignition guard. Each red state below was
# reproduced first on the built module and is banked in
# docs/i246-ambipolar-wall-operator/evidence/round96_before.log.


def _slow_neutral_drift_reactor(slow_k, xseed, source=1.0e5, gamma=1.0,
                                te_ev=TE_NOMINAL_EV, termination=None):
    """A source-driven discharge (electron ignites from exactly zero and saturates to
    S/nu_wall) plus a slow, charge-decoupled neutral drain ``Ar -> X`` at rate constant
    ``slow_k`` on top of a live ``xseed`` fraction of X. X stands in for a slow neutral
    channel (a metastable excitation, an isomerisation); it drains on a ~1/slow_k timescale,
    so with ``slow_k`` tiny the electron reaches steady density while X is still filling --
    a saturated electron sitting over a neutral that has NOT reached steady state. X carries
    a distinct heavy skeleton so the Ar+ wall recycle stays unambiguous (Ar+ -> Ar)."""
    e = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    x = Species(label='X').from_adjacency_list('1 He u0 p1 c0')
    x.thermo = _thermo_with_h298(0.0)
    imf = {e: 0.0, arp: 0.0, ar: 1.0 - xseed, x: xseed}
    r = PlasmaReactor(
        (TGAS, 'K'), (P_NOMINAL, 'Pa'), imf, (te_ev * EV_TO_K, 'K'), n_sims=1,
        termination=termination or [], diffusion_length=(_diffusion_length(), 'm'),
        ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'), wall_recycling=gamma,
        ionisation_source=(source, 'm^-3/s'),
        # Ar and X carry distinct skeletons: opt into the single-bath approximation, the
        # transport is not what this drifting-neutral criterion test is about.
        wall_single_bath_approximation=True)
    core = [e, ar, arp, x]
    slow = Reaction(reactants=[ar], products=[x], reversible=False,
                    kinetics=Arrhenius(A=(slow_k, 's^-1'), n=0, Ea=(0, 'J/mol')))
    r.initialize_model(core, [slow], [], [])
    return r, core, [slow]


def test_a_saturated_electron_does_not_vouch_for_a_still_drifting_neutral():
    """Round 96 HIGH 1: steady_state_external_armed armed the WHOLE criterion once the
    electron passed t*nu_wall >= 1, so a slow neutral reaction that had not run through its
    own timescale was declared stationary the moment the electron saturated. Arming is
    per-quantity: the electron's arm vouches only for the electron and may license the
    system only while the generic (neutral) residual is not still RISING toward its own arm.
    A neutral mid-transit -- X filling on a ~1e12 s timescale, its residual below tolerance
    only because it has not yet reached that scale -- must keep the system NOT steady, while
    a deck with no such neutral still reports steady once the electron saturates. Both, in
    one build, so the fix cannot buy one by breaking the other."""
    term = [TerminationSteadyState(tolerance=1.0e-8), TerminationTime((300.0, 's'))]
    r, core, rxns = _slow_neutral_drift_reactor(slow_k=1.0e-12, xseed=1.0e-2, termination=term)
    _simulate(r, core, rxns)
    assert not r.steady_state_reached, (
        "a saturated electron wrongly vouched for a neutral still draining Ar -> X")

    term2 = [TerminationSteadyState(tolerance=1.0e-8), TerminationTime((300.0, 's'))]
    r2, core2, rxns2 = _build_reactor(wall=True, with_chemistry=False, x_ion=0.0,
                                      source=1.0e5, termination=term2)
    _simulate(r2, core2, rxns2)
    assert r2.steady_state_reached, (
        "the pure saturating discharge with no drifting neutral was not recognised steady")


def test_a_stationary_composition_over_a_shrinking_inventory_terminates():
    """Round 96 HIGH 2: the generic criterion measures mole FRACTIONS while the external
    electron residual measured electron MOLES. On a pumped discharge whose fractions go
    stationary but whose absolute inventory shrinks, the moles residual reads large motion
    where the fraction is flat, poisoning the MAX fold so the run never terminated on a
    genuinely stationary composition. Measured on the same intensive quantity (the electron
    mole fraction), the state is recognised as steady and the run stops before the backstop."""
    term = [TerminationSteadyState(tolerance=1.0e-8), TerminationTime((50.0, 's'))]
    r, core, rxns = _build_reactor(wall=True, gamma=0.0, with_chemistry=False,
                                   x_ion=1.0e-4, source=1.0e22, termination=term)
    _simulate(r, core, rxns)
    assert r.steady_state_reached, (
        "a stationary composition over a shrinking inventory was not recognised as steady")
    assert r.t < 50.0, "the run reached the time backstop instead of terminating on steady state"


def test_the_wall_guard_evaluates_the_runtime_expression_not_a_reference_proxy():
    """Round 96 MEDIUM 1: the construction guard evaluated nu_wall at the reference density,
    where mu_i is exactly mu0, but run time forms mu_i = mu0*Nref/n_neutral. A finite
    mu0=1e20 with a finite Nref=1e308 leaves mu0 finite yet mu0*Nref = inf, so nu_wall is
    infinite for every real state while the reference-density proxy read finite and admitted
    it. The guard must evaluate the run-time expression at its worst case (the neutral
    floor). The physical single-bath case must still construct."""
    e, ar, arp = _argon_species()
    imf = {e: 1.0e-6, arp: 1.0e-6, ar: 1.0 - 2.0e-6}

    def _make(mu0, nref):
        return PlasmaReactor(
            (TGAS, 'K'), (P_NOMINAL, 'Pa'), imf, (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1,
            termination=[], diffusion_length=(_diffusion_length(), 'm'),
            ion_reduced_mobility=(mu0, 'm^2/(V*s)'),
            mobility_reference_density=(nref, 'm^-3'), wall_recycling=1.0)

    with pytest.raises(PlasmaStateError):
        _make(1.0e20, 1.0e308)          # each finite; the product mu0*Nref overflows
    _make(MU0_AR_IN_AR, PLASMA_LOSCHMIDT)   # the physical case constructs


def test_multi_gas_bath_records_an_availability_state_not_just_a_warning():
    """Round 96 MEDIUM 2, kept under round 100's opt-in: once the single-bath approximation
    is opted into (wall_single_bath_approximation=True), a consumer reading the latched wall
    fluxes must still see the approximation as a queryable STATE, not merely a log line. When
    the neutral bath spans more than one heavy skeleton the single ion mobility is applied to
    the summed density, so every nu_wall-derived flux is downgraded from 'available' to
    'available-single-bath-approximation' in the availability dict. A single-skeleton bath
    (Ar, or Ar and its metastable) reports plain 'available'. The downgrade never launders a
    genuinely unavailable (NaN) field into a usable one."""
    e, ar, arp = _argon_species()
    he = Species(label='He').from_adjacency_list('1 He u0 p1 c0')
    he.thermo = _thermo_with_h298(0.0)
    imf = {e: 1.0e-6, arp: 1.0e-6, ar: 0.5 - 1.0e-6, he: 0.5 - 1.0e-6}
    reactor = PlasmaReactor(
        (TGAS, 'K'), (P_NOMINAL, 'Pa'), imf, (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1,
        termination=[], diffusion_length=(_diffusion_length(), 'm'),
        ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'), wall_recycling=1.0,
        wall_single_bath_approximation=True)
    reactor.initialize_model([e, ar, he, arp], [], [], [])
    assert reactor.wall_bath_is_mixture is True
    y = np.array(reactor.y0[:reactor.num_core_species], float)
    reactor.residual(0.0, y, np.zeros(reactor.num_core_species, float))
    avail = reactor.wall_energy_availability
    assert avail['wall_flux'] == 'available-single-bath-approximation'
    assert avail['wall_electron_energy_flux'] == 'available-single-bath-approximation'

    # single-gas control: plain 'available', no downgrade
    r1, _, _ = _build_reactor(wall=True, with_chemistry=False)
    y1 = _state_at(r1, 1.0e-6)
    r1.residual(0.0, y1, np.zeros(r1.num_core_species, float))
    assert r1.wall_bath_is_mixture is False
    assert r1.wall_energy_availability['wall_flux'] == 'available'


def test_a_subnormal_source_that_injects_nothing_is_refused():
    """Round 96 LOW: a positive but subnormal ionisation_source (5e-324, 1e-320) -- or any
    value whose volumetric molar rate source/Na underflows -- reads as source > 0.0 and so
    switches off the zero-electron ignition guard, yet source*V/Na injects exactly zero: the
    deck declares ignition-from-zero it can never achieve. Refuse it at construction. A
    normal source injects and is admitted."""
    for bad in (5.0e-324, 1.0e-320, 1.0e-290):
        with pytest.raises(PlasmaStateError):
            _build_reactor(wall=True, with_chemistry=False, x_ion=0.0, source=bad)
    _build_reactor(wall=True, with_chemistry=False, x_ion=0.0, source=1.0e5)


# ---------------------------------------------------------------- round 100
#
# The steady-state MECHANISM was unsound even though its behaviour was right: the
# external arm was a two-point comparison with a permanent latch, persistence still
# depended on accepted solver steps, and the physical window was anchored to absolute
# log time instead of the system relaxation time. The source guard evaluated a
# reference-density proxy, and the multi-gas mobility was documented rather than gated.
# Each red state below was reproduced first on the built module and is banked in
# docs/i246-ambipolar-wall-operator/evidence/round100_before.log.


def test_the_source_guard_evaluates_the_runtime_expression_at_the_actual_volume():
    """Round 100 MEDIUM: the __init__ source guard checks ionisation_source/Na, but the
    residual injects ionisation_source*V/Na -- a different expression, the same shape as
    round 96's nu_wall overflow. At an extreme-but-finite volume the runtime product
    underflows to exactly zero (a declared source that injects nothing, switching off the
    zero-electron ignition guard while n_e can never leave zero) though source/Na is a
    normal double that passes __init__. The guard must evaluate the run-time expression at
    the actual initial volume. The physical case at a normal volume still constructs."""
    e, ar, arp = _argon_species()
    imf = {e: 0.0, arp: 0.0, ar: 1.0}

    def make(P, source):
        r = PlasmaReactor((TGAS, 'K'), (P, 'Pa'), imf, (TE_NOMINAL_EV * EV_TO_K, 'K'),
                          n_sims=1, termination=[], diffusion_length=(_diffusion_length(), 'm'),
                          ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                          wall_recycling=0.0, ionisation_source=(source, 'm^-3/s'))
        r.initialize_model([e, ar, arp], [], [], [])
        return r

    # source/Na = 1.66e-284 is a normal double (passes __init__), but at V ~ 2.5e-37 m^3
    # the runtime source*V/Na ~ 4e-321 underflows below the smallest normal double.
    with pytest.raises(PlasmaStateError):
        make(1.0e40, 1.0e-260)
    make(P_NOMINAL, 1.0e5)          # the physical case constructs


def test_wall_relaxation_time_is_one_over_nu_wall_and_absent_without_a_wall():
    """Round 100 HIGH 3: the criterion anchors persistence to the SYSTEM relaxation time,
    which the reactor supplies. A wall deck reports 1/nu_wall (finite, positive, evaluated
    at the state); a wall-less deck reports nan, sending the criterion to its absolute-time
    fallback so ordinary reactors are unchanged."""
    r, core, _ = _build_reactor(wall=True, with_chemistry=False)
    y = _state_at(r, 1.0e-6)
    V = r.compute_volume(y)
    nu = r.compute_nu_wall(y, V)
    tau = r.steady_state_relaxation_time(1.0, y)
    assert np.isfinite(tau) and tau > 0.0
    assert np.isclose(tau, 1.0 / nu, rtol=1e-12)

    r2, core2, _ = _build_reactor(wall=False, with_chemistry=False)
    y2 = _state_at(r2, 1.0e-6)
    assert not np.isfinite(r2.steady_state_relaxation_time(1.0, y2))


def test_all_four_steady_state_verdicts_hold_together_in_one_build():
    """Round 100 verifier: the four steady-state verdicts must coexist under ONE build -- a
    fix for one must not buy itself by breaking another. (1) a slow-drifting neutral is NOT
    steady; (2) a stationary composition over a shrinking inventory IS steady; (3) a
    zero-seed ignition saturates to a steady state; (4) a genuinely inert deck is NOT
    steady."""
    t1 = [TerminationSteadyState(tolerance=1.0e-8), TerminationTime((300.0, 's'))]
    r1, c1, x1 = _slow_neutral_drift_reactor(slow_k=1.0e-12, xseed=1.0e-2, termination=t1)
    _simulate(r1, c1, x1)
    assert not r1.steady_state_reached, "(1) a still-drifting neutral was wrongly called steady"

    t2 = [TerminationSteadyState(tolerance=1.0e-8), TerminationTime((50.0, 's'))]
    r2, c2, x2 = _build_reactor(wall=True, gamma=0.0, with_chemistry=False,
                                x_ion=1.0e-4, source=1.0e22, termination=t2)
    _simulate(r2, c2, x2)
    assert r2.steady_state_reached, "(2) a stationary composition over a shrinking inventory was missed"

    t3 = [TerminationSteadyState(tolerance=1.0e-8), TerminationTime((300.0, 's'))]
    r3, c3, x3 = _build_reactor(wall=True, with_chemistry=False, x_ion=0.0,
                                source=1.0e5, termination=t3)
    _simulate(r3, c3, x3)
    assert r3.steady_state_reached, "(3) a zero-seed ignition did not reach steady state"

    t4 = [TerminationSteadyState(tolerance=1.0e-8), TerminationTime((1.0, 's'))]
    r4, c4, x4 = _build_reactor(wall=False, with_chemistry=False, termination=t4)
    _simulate(r4, c4, x4)
    assert not r4.steady_state_reached, "(4) a genuinely inert deck was wrongly called steady"


# ---------------------------------------------------------------- round 106 MEDIUM


def _wall_reactor_with_flag(flag):
    """A wall PlasmaReactor built with an explicit wallSingleBathApproximation value,
    bypassing _build_reactor (which does not expose the flag)."""
    electron, ar, arp = _argon_species()
    imf = {electron: 1.0e-6, arp: 1.0e-6, ar: 1.0 - 2.0e-6}
    return PlasmaReactor(
        (TGAS, 'K'), (P_NOMINAL, 'Pa'), imf, (TE_NOMINAL_EV * EV_TO_K, 'K'),
        n_sims=1, termination=[],
        diffusion_length=(_diffusion_length(), 'm'),
        ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
        wall_recycling=1.0, wall_single_bath_approximation=flag)


def test_wall_single_bath_approximation_is_coerced_by_value_not_truthiness():
    """Round 106 MEDIUM: wall_single_bath_approximation was set with ``bool(value)``, which
    reads the STRING "False" as True (bool of any non-empty string is True), reads 2 and NaN
    as True, and so opts a deck into the single-bath approximation it wrote "False" to
    decline -- silently. It must be coerced by VALUE, the same discipline quasineutralElectron
    already gets: a genuine bool/None passes through, a boolean-like string parses by its
    meaning, and anything else is refused by name rather than read as truthy."""
    assert _wall_reactor_with_flag(True).wall_single_bath_approximation is True
    assert _wall_reactor_with_flag(False).wall_single_bath_approximation is False
    # the exact trap: the string "False" must NOT enable the flag
    assert _wall_reactor_with_flag('False').wall_single_bath_approximation is False
    assert _wall_reactor_with_flag('True').wall_single_bath_approximation is True
    # a non-boolean-like value is refused, not silently read as True
    for bad in (2, float('nan'), 0.5):
        with pytest.raises(PlasmaStateError):
            _wall_reactor_with_flag(bad)


def test_coerce_bool_flag_accept_reject_set_is_enumerated_and_the_refusal_reaches_the_user():
    """Round 109 MEDIUM 3: enumerate exactly what _coerce_bool_flag accepts and rejects, and
    confirm a rejection is a REFUSAL the user sees (a raised PlasmaStateError that propagates out
    of the constructor), not a log line the run continues past. Includes the motivating case
    ("False"), plus 2, NaN, None and the empty string."""
    from rmgpy.solver.plasma import _coerce_bool_flag

    # ACCEPTED, coerced by value (never by truthiness). None and "" mean opt-out, not refusal;
    # "False" -- the case that motivated the helper -- must be False, not the True that bool()
    # would give.
    accepted = {
        True: True, False: False, None: False,
        'True': True, 'true': True, 'TRUE': True, '1': True, 'yes': True, 'on': True,
        'False': False, 'false': False, '0': False, 'no': False, 'off': False,
        '': False, '  Off  ': False,
    }
    for value, expected in accepted.items():
        result = _coerce_bool_flag(value, 'flag', 'id')
        assert result is expected, (value, result, expected)

    # REJECTED: refused by name with a PlasmaStateError. A number opts in on its truthiness
    # (bool(2), bool(0.5), bool(NaN) are all True), an unrecognised string is a likely typo, and
    # a container/object has no boolean meaning here -- all must fail loudly, never coerce.
    for value in (2, 0, 0.5, float('nan'), float('inf'), 'maybe', 'Flase', [], {}, object()):
        with pytest.raises(PlasmaStateError):
            _coerce_bool_flag(value, 'flag', 'id')

    # ...and the refusal actually reaches the user: it propagates out of the reactor
    # constructor rather than being swallowed into a log line the run ignores.
    for bad in (2, float('nan'), 'maybe'):
        with pytest.raises(PlasmaStateError):
            _wall_reactor_with_flag(bad)


def test_the_ionisation_source_is_validated_at_the_evolved_volume_not_only_the_initial_one():
    """Round 106 MEDIUM: the __init__/set_initial guard validated source*V/Na only at the
    INITIAL volume, while the residual and Jacobian recompute source*V/Na at the evolved and
    Newton-trial volume -- guard and computation on different expressions at different points,
    the campaign's recurring shape. A source finite and positive at V0 can overflow (or
    underflow to zero) once scaled by an extreme trial volume, and was then applied as an
    infinite (or vanishing) rate silently. Now every site forms the product through one
    validated helper, so a volume that makes it non-finite is caught wherever the solve
    reaches it."""
    r, _, _ = _build_reactor(wall=True, with_chemistry=False, source=1.0e18, x_ion=1.0e-6)
    ie, i_ar, i_arp = _indices(r)
    # Valid at the initial volume: initialisation already succeeded above.
    y0 = _state_at(r, 1.0e-6)
    source_mol_0 = r.ionisation_source.value_si * r.compute_volume(y0) / constants.Na
    assert np.isfinite(source_mol_0) and source_mol_0 > 0.0

    # A wild Newton-trial neutral amount inflates the volume until source*V/Na overflows.
    # Both the residual and the Jacobian must refuse it, not propagate an infinite source.
    y = _state_at(r, 1.0e-6)
    y[i_ar] = 1.0e300
    V = r.compute_volume(y)
    assert not np.isfinite(r.ionisation_source.value_si * V / constants.Na), \
        "the trial volume must actually overflow source*V/Na for this to test the guard"
    dydt = np.zeros_like(y)
    with pytest.raises(PlasmaStateError):
        r.residual(0.0, y, dydt)
    with pytest.raises(PlasmaStateError):
        r.jacobian(0.0, y, dydt, 0.0)


def _solver_test_functions(source):
    """(name, FunctionDef) for every ``test_*`` in a test-file source."""
    import ast
    return [(n.name, n) for n in ast.walk(ast.parse(source))
            if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and n.name.startswith('test_')]


def _tests_without_a_real_assertion(source):
    """Names of ``test_*`` functions in `source` that carry NO non-vacuous, REACHABLE assertion
    (round 111 LOW 2, tightened round 112 LOW). A docstring -- even one naming a round/finding tag
    -- with a ``pass`` body asserts nothing; a banked red or a property check must contain a real
    assertion, so a tag can no longer stand in for one (the round-106 census accepted 'tag OR
    assert', which a tagged ``pass`` satisfied).

    ``ast.walk`` over the whole function is too permissive -- it counted four things that back no
    claim, closed here (round 112 LOW):
      * a vacuous ``assert`` whose test is decided at PARSE TIME, referencing no runtime value:
        ``assert True`` but also ``assert 1 == 1`` (the round-111 census only caught a bare
        ``ast.Constant``, so a constant COMPARE slipped through);
      * an ``assert`` in a statically DEAD branch (``if False:`` / ``if 0:``) that never runs;
      * an ``assert`` inside an UNCALLED nested function, which ``ast.walk`` reaches but the test
        body never executes;
      * a context manager whose AST merely CONTAINS the substring ``raises``/``warns`` (e.g. a
        helper named ``a_thing_that_raises``), rather than an actual ``pytest.raises`` / ``.warns``.
    So the detector walks reachable statements only -- never into a nested def/class or a dead
    branch -- and recognises ``raises``/``warns`` by the call TARGET, not a substring."""
    import ast

    def is_vacuous(test):
        # An assert whose test references no runtime value (only literals and operators on them)
        # is decided at parse time: ``assert True``, ``assert 1 == 1``, ``assert 1 < 2 and 3``.
        for sub in ast.walk(test):
            if isinstance(sub, (ast.Name, ast.Call, ast.Attribute, ast.Subscript, ast.Starred)):
                return False
        return True

    def is_pytest_raises_or_warns(item):
        # Only a genuine ``pytest.raises`` / ``pytest.warns`` context manager backs a claim
        # (round 113 LOW). A bare ``raises(...)`` need not be pytest's, and an unrelated
        # ``fake.raises(...)`` merely contains the word -- match the call TARGET as
        # ``pytest.<raises|warns>``, not any attribute or name that ends in raises/warns.
        expr = item.context_expr
        if not isinstance(expr, ast.Call):
            return False
        func = expr.func
        return (isinstance(func, ast.Attribute) and func.attr in ('raises', 'warns')
                and isinstance(func.value, ast.Name) and func.value.id == 'pytest')

    def const_truth(test):
        # True/False if `test` is a compile-time constant truth value, else None. A test that
        # references no runtime value folds at parse time, so ``if 1 == 0:`` is as dead as
        # ``if False:`` (round 114 LOW); it holds only literals and operators, so eval is safe.
        if not is_vacuous(test):
            return None
        try:
            return bool(eval(compile(ast.Expression(test), '<census>', 'eval'), {'__builtins__': {}}))
        except Exception:
            return None

    def empty_iter(node):
        # A statically-empty iterable, so a ``for`` body over it never runs: [], (), {}, an
        # empty string/bytes literal, or ``range(0)`` (round 113 LOW).
        if isinstance(node, (ast.List, ast.Tuple, ast.Set)):
            return len(node.elts) == 0
        if isinstance(node, ast.Dict):
            return len(node.keys) == 0
        if isinstance(node, ast.Constant):
            return isinstance(node.value, (str, bytes)) and len(node.value) == 0
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id == 'range'):
            return (len(node.args) == 1 and isinstance(node.args[0], ast.Constant)
                    and node.args[0].value == 0)
        return False

    def walk(stmts):
        for stmt in stmts:
            # Nothing after a return/raise/break/continue in the same block runs (round 114 LOW).
            if isinstance(stmt, (ast.Return, ast.Raise, ast.Break, ast.Continue)):
                return False
            # An uncalled nested definition never runs in the test body: do not credit its asserts.
            if isinstance(stmt, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                continue
            if isinstance(stmt, ast.Assert):
                if not is_vacuous(stmt.test):
                    return True
                continue
            if isinstance(stmt, ast.With):
                if any(is_pytest_raises_or_warns(it) for it in stmt.items):
                    return True
                if walk(stmt.body):
                    return True
                continue
            if isinstance(stmt, ast.If):
                truth = const_truth(stmt.test)
                if truth is True:
                    if walk(stmt.body):     # the else-branch is statically dead
                        return True
                elif truth is False:
                    if walk(stmt.orelse):   # the if-branch is statically dead
                        return True
                elif walk(stmt.body) or walk(stmt.orelse):
                    return True
                continue
            if isinstance(stmt, ast.While):
                truth = const_truth(stmt.test)
                if truth is False:          # ``while False:`` -- body never runs, else does
                    if walk(stmt.orelse):
                        return True
                elif walk(stmt.body) or walk(stmt.orelse):
                    return True
                continue
            if isinstance(stmt, ast.For):
                if empty_iter(stmt.iter):   # ``for _ in []:`` -- body never runs, else does
                    if walk(stmt.orelse):
                        return True
                elif walk(stmt.body) or walk(stmt.orelse):
                    return True
                continue
            # Any other compound statement (try/with-less block): descend into its
            # reachable child bodies, but never into a nested def/class (handled above).
            for field in ('body', 'orelse', 'finalbody'):
                child = getattr(stmt, field, None)
                if child and walk(child):
                    return True
            for handler in getattr(stmt, 'handlers', []) or []:
                if walk(handler.body):
                    return True
        return False

    # An ``async def`` test is never awaited without an async plugin, which this suite does not
    # load: pytest skips it with a warning, so none of its asserts run (round 114 LOW).
    return [name for name, node in _solver_test_functions(source)
            if isinstance(node, ast.AsyncFunctionDef) or not walk(node.body)]


def test_every_solver_test_asserts_a_property():
    """Round 106 census, TIGHTENED for round 111 (LOW 2). The census made the evidence standard
    executable -- a test is a banked defect reproduction or a property check -- but accepted a
    docstring naming a round/finding tag OR an assertion, so a tagged docstring with a ``pass``
    body satisfied it while asserting nothing. A tag is prose about intent; only an assertion
    backs a claim. Every ``test_*`` function in the two solver test files must now contain a
    real (non-vacuous) assertion, tag or no tag."""
    here = os.path.dirname(os.path.abspath(__file__))
    files = [os.path.join(here, 'plasmaWallTest.py'),
             os.path.join(here, 'steadyStateTest.py')]
    total = 0
    offenders = []
    for path in files:
        source = open(path).read()
        total += len(_solver_test_functions(source))
        offenders += ['{0}::{1}'.format(os.path.basename(path), name)
                      for name in _tests_without_a_real_assertion(source)]
    assert total > 0, "no test functions were discovered -- the census enforcement is vacuous"
    assert not offenders, (
        "{0} of {1} test functions across the two solver test files contain no non-vacuous "
        "assertion, so they back no claim -- neither a banked red nor a property: {2}".format(
            len(offenders), total, offenders))


def test_the_census_rejects_a_tagged_docstring_with_no_assertion():
    """Round 111 LOW 2, the census's own tripwire. Before the tightening a ``test_*`` whose body
    was a round/finding-tagged docstring and ``pass`` satisfied the census, because a tag counted
    as evidence -- the ``is_tagged or is_asserting`` rule. Feed exactly that and confirm the
    detector now reports it, alongside a test that only asserts a bare constant; and confirm a
    real assertion and a ``pytest.raises`` context manager are still accepted. This is what the
    round-106 census could not do to itself."""
    rejected = (
        'def test_tagged_but_empty():\n'
        '    """Round 111 HIGH 1 -- names a finding but asserts nothing."""\n'
        '    pass\n'
        'def test_vacuous_only():\n'
        '    """LOW 2 tag."""\n'
        '    assert True\n'
        # round 112 LOW: four holes the round-111 census still had.
        'def test_vacuous_compare():\n'                       # `1 == 1` is decided at parse time
        '    assert 1 == 1\n'
        'def test_dead_branch_assert():\n'                    # under `if False:` -- never runs
        '    if False:\n'
        '        assert compute() == 1\n'
        'def test_uncalled_nested_assert():\n'                # the assert is in an uncalled inner def
        '    def inner():\n'
        '        assert compute() == 1\n'
        'def test_unrelated_raises_context_manager():\n'      # AST holds "raises" but is not pytest.raises
        '    with a_helper_that_raises():\n'
        '        do_work()\n'
        # round 113 LOW: four more the round-112 census still credited.
        'def test_bare_raises_cm():\n'                        # bare raises(...) -- need not be pytest's
        '    with raises(ValueError):\n'
        '        do_work()\n'
        'def test_fake_raises_cm():\n'                        # fake.raises(...) -- not pytest.raises
        '    with fake.raises(ValueError):\n'
        '        do_work()\n'
        'def test_while_false_assert():\n'                    # while False -- body never runs
        '    while False:\n'
        '        assert compute() == 1\n'
        'def test_empty_loop_assert():\n'                     # for _ in [] -- body never runs
        '    for _ in []:\n'
        '        assert compute() == 1\n'
        # round 114 LOW: three more the round-113 census still credited.
        'def test_folded_dead_branch_assert():\n'            # `1 == 0` folds to False -- never runs
        '    if 1 == 0:\n'
        '        assert compute() == 1\n'
        'def test_assert_after_return():\n'                  # unreachable after the return
        '    return\n'
        '    assert compute() == 1\n'
        'async def test_async_assert():\n'                   # no async plugin: pytest never awaits it
        '    assert compute() == 1\n'
    )
    accepted = (
        'def test_real_runtime_value():\n'                    # references a runtime value -> real
        '    assert len([1, 2]) == 2\n'
        'def test_context_manager():\n'
        '    with pytest.raises(ValueError):\n'
        '        raise ValueError()\n'
        'def test_warns_context_manager():\n'
        '    with pytest.warns(UserWarning):\n'
        '        do_work()\n'
        'def test_real_assert_inside_a_live_loop():\n'        # a reachable assert in a loop body
        '    for i in range(3):\n'
        '        assert step(i) is True\n'
    )
    assert set(_tests_without_a_real_assertion(rejected)) == {
        'test_tagged_but_empty', 'test_vacuous_only', 'test_vacuous_compare',
        'test_dead_branch_assert', 'test_uncalled_nested_assert',
        'test_unrelated_raises_context_manager', 'test_bare_raises_cm', 'test_fake_raises_cm',
        'test_while_false_assert', 'test_empty_loop_assert', 'test_folded_dead_branch_assert',
        'test_assert_after_return', 'test_async_assert'}
    assert _tests_without_a_real_assertion(accepted) == []


def test_a_subnormal_electron_fraction_does_not_make_the_external_residual_infinite():
    """Round 110 HIGH 1 addendum. ``steady_state_external_residual`` guards ``ne_now > 0.0``
    on the electron MOLES, then divides by the total and takes ``log()`` of the FRACTION. A
    positive subnormal electron moles passes the moles guard, yet ``xe = ne/tot`` underflows to
    exactly 0.0, so ``log(0) = -inf`` and the compiled hook returns ``+inf``. rework.md claimed
    this channel can never be infinite; the compiled hook proves otherwise. An unresolvable
    fraction is no usable number -- like a sub-floor electron -- so the hook must return ``nan``
    (its documented 'no information' sentinel), keeping the channel vocabulary to {nan, finite}
    so a generic ``inf`` can never arrive on the external side."""
    reactor = _build_reactor(source=1.0e18, wall=True, gamma=1.0)[0]
    ie = reactor.electron_index
    n = reactor.num_core_species

    y_prev = np.full(n, 1.0)
    y_prev[ie] = 1.0e-4                        # a normal, finite previous electron fraction
    y_now = np.full(n, 10.0)                   # order-10 total inventory
    y_now[ie] = 5.0e-324                        # smallest positive subnormal: passes ne_now > 0

    # Premise, empirically: the moles are positive but the fraction underflows to exact zero.
    assert y_now[ie] > 0.0
    assert y_now[ie] / y_now.sum() == 0.0

    val = reactor.steady_state_external_residual(2.0, y_now, 1.0, y_prev)
    assert not np.isinf(val), "a subnormal electron fraction still yields an infinite residual"
    assert np.isnan(val), "an unresolvable electron fraction must read as no information (nan)"


def _char_rate_division_lines(src):
    """Every line of a Cython source on which it divides by the bare name ``char_rate`` in
    CODE -- comments and string/docstring literals excluded via the tokenizer, so a mention of
    ``char_rate`` in prose or dead code cannot register and a behaviourally-identical rewrite
    cannot hide. base.pyx is not valid Python (``cdef`` etc.), so this tokenizes rather than
    ``ast.parse``-s; the token stream is enough to find a ``/`` immediately followed by the
    ``char_rate`` NAME."""
    import io
    import tokenize

    toks = tokenize.generate_tokens(io.StringIO(src).readline)
    code = [t for t in toks if t.type not in (
        tokenize.COMMENT, tokenize.STRING, tokenize.NL, tokenize.NEWLINE,
        tokenize.INDENT, tokenize.DEDENT)]
    return [b.start[0] for a, b in zip(code, code[1:])
            if a.type == tokenize.OP and a.string == '/'
            and b.type == tokenize.NAME and b.string == 'char_rate']


def _method_line_span(src, header_regex):
    """The [start, end) line span of the def whose header matches ``header_regex``, ended by the
    next def/cpdef/cdef at the same or lower indent. Used to whitelist a guarded, display-only
    region structurally rather than by hard-coded line numbers."""
    lines = src.splitlines()
    start = next(i for i, l in enumerate(lines, 1) if re.match(header_regex, l))
    indent = len(lines[start - 1]) - len(lines[start - 1].lstrip())
    for i in range(start, len(lines)):
        l = lines[i]
        if l.strip() and (len(l) - len(l.lstrip())) <= indent and re.match(r'\s*(cpdef|cdef|def)\b', l):
            return start, i + 1
    return start, len(lines) + 1


def test_no_enlargement_ratio_divides_by_char_rate_outside_the_abstaining_helper():
    """Round 110 census, REPAIRED for round 111 (LOW 1). The round-110 cross-channel census
    asserted source SUBSTRINGS for the folds it already knew about -- so a comment carrying the
    string satisfied it, a behaviourally-identical rewrite broke it, and, decisively, it could
    not see a division site nobody had listed. That is why it MISSED the surface-species ratio
    (round 111 HIGH 2), which divided ``max(production, consumption) / char_rate`` directly:
    the census's own round finding it could not detect.

    The repair ENUMERATES rather than lists. Every code division by ``char_rate`` in base.pyx is
    found structurally (via the tokenizer, comments and strings excluded) and must fall inside
    ``log_rates`` -- the only legitimate site, a display line guarded by ``if char_rate == 0.0``.
    Any division outside it is an enlargement/pruning/promotion decision dividing by a rate that
    can be zero, which must instead route through ``_rate_ratios_or_zero`` (finite, positive
    denominator or abstain). The surface-ratio bug lands OUTSIDE ``log_rates`` and fails here;
    after the fix the count outside is zero."""
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(os.path.dirname(os.path.dirname(here)))
    with open(os.path.join(root, 'rmgpy', 'solver', 'base.pyx')) as fh:
        base_src = fh.read()
    with open(os.path.join(root, 'rmgpy', 'solver', 'plasma.pyx')) as fh:
        plasma_src = fh.read()

    division_lines = _char_rate_division_lines(base_src)
    log_start, log_end = _method_line_span(base_src, r'\s*cpdef\s+log_rates\b')
    outside = [ln for ln in division_lines if not (log_start <= ln < log_end)]
    assert not outside, (
        "base.pyx divides by char_rate outside the guarded display method log_rates, at "
        "line(s) {0} -- an enlargement/pruning/promotion ratio dividing by a rate that can be "
        "zero. Route it through _rate_ratios_or_zero (round 111 HIGH 2); every such division "
        "found: {1}".format(outside, division_lines))
    # The sanctioned route exists, and the FORBIDDEN dimensional denominators (owner's ruling:
    # 1.0, an epsilon, or any floor) stay gone -- a floored denominator compares a dimensional
    # rate against a dimensionless tolerance (round 110 HIGH 2).
    assert '_rate_ratios_or_zero' in base_src, "the abstaining ratio helper is missing"
    assert 'char_rate if char_rate > 0.0 else 1.0' not in base_src, \
        "the dimensional edge-ratio denominator (round 110 HIGH 2) was reintroduced"
    assert 'ratio_denom' not in base_src, \
        "a floored ratio denominator was reintroduced; the ratio criterion must abstain, not floor"
    # The chemistry / non-chemical flux folds this census also catalogues remain present.
    assert 'char_rate * char_rate' in base_src
    assert 'non_chemical_char_rate * non_chemical_char_rate' in base_src
    assert '_apply_wall_terms(y, V, res)' in plasma_src


def test_zero_core_flux_with_wall_loss_abstains_on_ratio_but_keeps_an_absolute_criterion():
    """Round 110 HIGH 2 (owner's ruling; matrix case 'non-zero wall loss with zero core flux').
    A wall-driven plasma with no gas-phase reactions has char_rate == 0, so the enlargement
    RATIO criterion abstains (see the denominator matrix in steadyStateTest). But the wall is
    moving the composition, so get_non_chemical_char_rate() > 0 and total_char_rate > 0: the
    dimensioned ABSOLUTE criterion is present and governs, and the run is NOT dismissed as one
    that never started on the strength of an undefined ratio. The wall physics itself is
    unchanged -- this only asserts which criterion is in force when the ratio is undefined."""
    r, core, rxns = _build_reactor(wall=True, with_chemistry=False, source=1.0e18)
    y = _state_at(r, 1.0e-6)                       # a live discharge state (electrons present)
    r.residual(0.0, y.copy(), np.zeros_like(y))
    char_rate = float(np.sqrt(np.sum(np.asarray(r.core_species_rates, float) ** 2)))

    # No gas chemistry: the ratio denominator is zero, so the ratio criterion abstains...
    assert char_rate == 0.0, "a chemistry-free deck must carry no gas-phase characteristic rate"
    # ...but the wall carries real flux, so the absolute criterion is live (total_char_rate > 0).
    assert r.get_non_chemical_char_rate() > 0.0, "the wall loss must supply the absolute criterion"

    # Drive the PRODUCTION path (simulate + total_char_rate), not just static probes: with the
    # ratio abstaining on the zero denominator, the run must still terminate deterministically on
    # the dimensioned absolute criterion rather than divide by zero or hang. Under strict float
    # handling a reverted ratio (bare ``rates / char_rate`` = 0/0) would raise here instead.
    r.termination = [TerminationSteadyState(tolerance=1e-6, window=2), TerminationTime((1.0e-4, 's'))]
    with np.errstate(invalid='raise', divide='raise'):
        result = _simulate(r, core, rxns)
    assert result[0] is True, "the wall-only discharge did not terminate through the absolute criterion"


# ---- I-269: neutral diffusion wall loss for DECLARED excited neutrals ----
#
# A declared excited neutral (Ar* below) is lost to the wall at nu_m = D_m/Lambda^2,
# with D_m = (D*N)_ref/n_gas scaled exactly like the ion mobility, and returns as the
# declared ground-state product. The reference value is supplied per species, as D*p
# or as D*N; nothing is inferred from the electronic state.

# cm^2 Torr/s, the deck's reference D*p: Wieme & Lenaerts, D(1 Torr) = 3.20e-3*T^1.68,
# 46.4 at 300 K and 47.0 at 302.2 K; full citation in documentation/source/users/rmg/
# input.rst (wallNeutralDiffusion). (The 54 once quoted corresponds to ~330 K.)
DP_AR_META = 47.0
AR_META_EV = 11.55                        # Ar(1s5) excitation energy, eV
CM2_TORR_TO_SI = 1.0e-4 * TORR_TO_PA      # (cm^2 Torr/s) -> (m^2 Pa/s)


def _meta_declaration(diffusivity=(DP_AR_META, 'cm^2*torr/s'), product='Ar', label='Ar*'):
    return {label: {'product': product, 'diffusivity': diffusivity}}


def _neutral_diffusion_reactor(declaration='default', tgas=TGAS, pressure=P_NOMINAL,
                               gamma=1.0, x_meta=1.0e-3, x_ion=1.0e-6, source=None,
                               ground_eV=0.0, meta_eV=AR_META_EV):
    """Ground Ar, metastable Ar*, Ar+, e- on a wall, with the wall's ion product pinned
    to ground Ar, and (by default) Ar* declared to diffuse to the wall as Ar. The
    states carry thermo H298 = ground_eV / meta_eV (None for no thermo)."""
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ground = _ground_species('Ar', ground_eV)
    meta = _metastable_species('Ar*', meta_eV)
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    imf = {electron: x_ion, arp: x_ion, meta: x_meta, ground: 1.0 - 2.0 * x_ion - x_meta}
    kwargs = dict(diffusion_length=(_diffusion_length(), 'm'),
                  ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                  wall_recycling=gamma,
                  wall_neutralization_products={'Ar+': 'Ar'})
    if source is not None:
        kwargs['ionisation_source'] = (source, 'm^-3/s')
    if declaration == 'default':
        declaration = _meta_declaration()
    if declaration is not None:
        kwargs['wall_neutral_diffusion'] = declaration
    reactor = PlasmaReactor((tgas, 'K'), (pressure, 'Pa'), imf, (TE_NOMINAL_EV * EV_TO_K, 'K'),
                            n_sims=1, termination=[], **kwargs)
    core = [electron, ground, meta, arp]
    reactor.initialize_model(core, [], [], [])
    return reactor, core


def _meta_state(reactor, core, n_meta=1.0e-3, n_ion=0.0, n_total=1.0):
    y = np.zeros(reactor.num_core_species, float)
    idx = {s.label: i for i, s in enumerate(core)}
    y[idx['Ar*']] = n_meta
    y[idx['Ar+']] = n_ion
    y[idx['e-']] = n_ion
    y[idx['Ar']] = n_total - n_meta - 2.0 * n_ion
    return y, idx


def test_neutral_wall_frequency_matches_hand_arithmetic():
    """nu_m = (D*p)/p/Lambda^2 at a state with no charged population, where the
    operator's neutral density is exactly p/(k_B T). Hand value from the declared
    47 cm^2 Torr/s, 5 Torr and the deck geometry; agreement to 1e-12 relative."""
    for tgas in (300.0, TGAS, 600.0):
        r, core = _neutral_diffusion_reactor(tgas=tgas)
        y, idx = _meta_state(r, core, n_ion=0.0)
        nu = r.compute_neutral_wall_frequencies(y, r.compute_volume(y))
        lam = _diffusion_length()
        hand = DP_AR_META * CM2_TORR_TO_SI / P_NOMINAL / (lam * lam)
        assert abs(nu[idx['Ar*']] / hand - 1.0) < 1e-12, (tgas, nu[idx['Ar*']], hand)
        for label in ('Ar', 'Ar+', 'e-'):
            assert nu[idx[label]] == 0.0


def test_neutral_wall_frequency_scales_inversely_with_pressure():
    lam = _diffusion_length()
    for p_torr in (1.0, 5.0, 20.0):
        r, core = _neutral_diffusion_reactor(pressure=p_torr * TORR_TO_PA)
        y, idx = _meta_state(r, core)
        nu = r.compute_neutral_wall_frequencies(y, r.compute_volume(y))[idx['Ar*']]
        assert abs(nu * p_torr * lam * lam / (DP_AR_META * 1.0e-4) - 1.0) < 1e-12


def test_neutral_wall_dn_form_equals_dp_form():
    """(D*N) = (D*p)/(k_B T) at the reactor gas temperature: the two spellings of the
    same reference coefficient give the same frequency."""
    # R/Na, not constants.kB: the engine converts with the same Boltzmann constant its
    # EOS uses, and the two differ at the 6e-8 level (constant vintages).
    dn_si = DP_AR_META * CM2_TORR_TO_SI / ((constants.R / constants.Na) * TGAS)
    r1, core1 = _neutral_diffusion_reactor()
    r2, core2 = _neutral_diffusion_reactor(
        declaration=_meta_declaration(diffusivity=(dn_si, '1/(m*s)')))
    y, idx = _meta_state(r1, core1, n_ion=1.0e-6)
    nu1 = r1.compute_neutral_wall_frequencies(y, r1.compute_volume(y))[idx['Ar*']]
    nu2 = r2.compute_neutral_wall_frequencies(y, r2.compute_volume(y))[idx['Ar*']]
    assert abs(nu1 / nu2 - 1.0) < 1e-12


def test_neutral_wall_residual_moves_metastable_to_ground_only():
    """The residual loses nu_m*y[Ar*] from Ar* and adds exactly that to Ar. Every
    charged row, and nu_wall, are bit-identical to the undeclared reactor."""
    rd, core_d = _neutral_diffusion_reactor()
    ru, core_u = _neutral_diffusion_reactor(declaration=None)
    y, idx = _meta_state(rd, core_d, n_ion=1.0e-6)
    zeros = np.zeros(rd.num_core_species, float)
    res_d = np.array(rd.residual(0.0, y, zeros)[0], float)
    wall_d = np.array(rd.wall_loss_rates, float)
    res_u = np.array(ru.residual(0.0, y, zeros)[0], float)
    wall_u = np.array(ru.wall_loss_rates, float)
    V = rd.compute_volume(y)
    nu_m = rd.compute_neutral_wall_frequencies(y, V)[idx['Ar*']]
    loss = nu_m * y[idx['Ar*']]
    assert loss > 0.0

    assert wall_u[idx['Ar*']] == 0.0
    assert wall_d[idx['Ar*']] == -loss
    assert abs((wall_d[idx['Ar']] - wall_u[idx['Ar']]) / loss - 1.0) < 1e-12
    for label in ('e-', 'Ar+'):
        assert wall_d[idx[label]] == wall_u[idx[label]]
        assert res_d[idx[label]] == res_u[idx[label]]
    assert rd.compute_nu_wall(y, V) == ru.compute_nu_wall(y, V)
    assert abs((res_u[idx['Ar*']] - res_d[idx['Ar*']]) / loss - 1.0) < 1e-12
    assert abs((res_d[idx['Ar']] - res_u[idx['Ar']]) / loss - 1.0) < 1e-12


def test_neutral_wall_loss_conserves_heavy_particles():
    rd, core = _neutral_diffusion_reactor(gamma=1.0)
    y, idx = _meta_state(rd, core, n_ion=1.0e-6)
    rd.residual(0.0, y, np.zeros(rd.num_core_species, float))
    wall = np.array(rd.wall_loss_rates, float)
    heavy = wall[idx['Ar']] + wall[idx['Ar*']] + wall[idx['Ar+']]
    assert abs(heavy) <= 1e-15 * abs(wall[idx['Ar*']])


def test_undeclared_metastable_has_no_neutral_wall_loss():
    r, core = _neutral_diffusion_reactor(declaration=None)
    y, idx = _meta_state(r, core, n_ion=1.0e-6)
    assert not np.any(r.compute_neutral_wall_frequencies(y, r.compute_volume(y)))
    r.residual(0.0, y, np.zeros(r.num_core_species, float))
    assert r.wall_loss_rates[idx['Ar*']] == 0.0


def test_neutral_wall_jacobian_matches_finite_difference():
    r, core = _neutral_diffusion_reactor(gamma=0.5)
    y, _ = _meta_state(r, core, n_ion=1.0e-6)
    best = _jacobian_scan(r, y)
    assert best < FD_TOLERANCE, best


def test_neutral_wall_jacobian_scan_can_fail():
    """Negative control: dropping the declared coefficient between the analytic
    Jacobian and the finite differences must break agreement, so the scan above
    actually sees the neutral wall term."""
    r, core = _neutral_diffusion_reactor(gamma=0.5)
    y, _ = _meta_state(r, core, n_ion=1.0e-6)

    def drop(reactor):
        reactor.wall_neutral_dn = np.zeros(reactor.num_core_species, float)

    assert _jacobian_scan(r, y, mutate=drop) > FD_TOLERANCE


def test_neutral_wall_declaration_survives_pickle_and_deepcopy():
    r, core = _neutral_diffusion_reactor()
    for clone in (pickle.loads(pickle.dumps(r)), copy.deepcopy(r)):
        assert set(clone.wall_neutral_diffusion) == {'Ar*'}
        assert clone.wall_neutral_diffusion['Ar*']['product'] == 'Ar'


@pytest.mark.parametrize("value", [float('nan'), float('inf'), float('-inf'), 0.0, -1.0,
                                   5.0e-324, 1.0e308])
def test_neutral_wall_refuses_nonfinite_or_nonpositive_diffusivity(value):
    """One mechanism, one message: NaN, +-inf, zero, negative, a subnormal, and a value
    whose worst-case frequency overflows are all refused by the same check."""
    with pytest.raises(PlasmaStateError) as exc:
        _neutral_diffusion_reactor(declaration=_meta_declaration(diffusivity=(value, 'cm^2*torr/s')))
    assert 'finite, positive' in str(exc.value)
    assert 'Ar*' in str(exc.value)


@pytest.mark.parametrize("label,declaration,fragment", [
    ("wrong dimension", _meta_declaration(diffusivity=(54.0, 'cm^2/s')), 'dimension'),
    ("not a dict", [('Ar*', 'Ar')], 'dict'),
    ("entry missing diffusivity", {'Ar*': {'product': 'Ar'}}, 'diffusivity'),
    ("entry missing product", {'Ar*': {'diffusivity': (54.0, 'cm^2*torr/s')}}, 'product'),
    ("unknown entry key", {'Ar*': {'product': 'Ar', 'diffusivity': (54.0, 'cm^2*torr/s'),
                                    'gamma': 1.0}}, 'gamma'),
    ("charged species declared", _meta_declaration(label='Ar+'), 'neutral'),
    ("product is the source", _meta_declaration(product='Ar*'), 'itself'),
    ("product not in core", _meta_declaration(product='Xe'), 'Xe'),
    ("product charged", _meta_declaration(product='Ar+'), 'neutral'),
])
def test_neutral_wall_refuses_malformed_declaration(label, declaration, fragment):
    with pytest.raises(PlasmaStateError) as exc:
        _neutral_diffusion_reactor(declaration=declaration)
    assert fragment in str(exc.value), (label, str(exc.value))


def test_neutral_wall_refuses_element_changing_product():
    """The wall returns the SAME atoms: a declared product with a different element
    composition is refused."""
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ground = _ground_species('Ar')
    meta = _metastable_species('Ar*')
    he = Species(label='He').from_adjacency_list('1 He u0 p1 c0')
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    imf = {electron: 1e-6, arp: 1e-6, meta: 1e-3, he: 0.1, ground: 0.9 - 1e-3 - 2e-6}
    reactor = PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf, (TE_NOMINAL_EV * EV_TO_K, 'K'),
                            n_sims=1, termination=[],
                            diffusion_length=(_diffusion_length(), 'm'),
                            ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                            wall_recycling=1.0, wall_neutralization_products={'Ar+': 'Ar'},
                            wall_single_bath_approximation=True,
                            wall_neutral_diffusion=_meta_declaration(product='He'))
    with pytest.raises(PlasmaStateError) as exc:
        reactor.initialize_model([electron, ground, meta, he, arp], [], [], [])
    assert 'element' in str(exc.value)


def test_neutral_wall_refuses_declaration_without_a_wall():
    electron, ar, arp = _argon_species()
    with pytest.raises(PlasmaStateError) as exc:
        PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), {ar: 1.0}, (TE_NOMINAL_EV * EV_TO_K, 'K'),
                      n_sims=1, termination=[], wall_neutral_diffusion=_meta_declaration())
    assert 'wall_neutral_diffusion' in str(exc.value)


def test_seeded_metastable_decays_exponentially_on_the_wall():
    """A seeded Ar* in an argon deck with no chemistry decays as exp(-nu_m t) through
    the production integrator, and argon atoms are conserved: ground Ar gains exactly
    what Ar* and Ar+ lose (gamma=1). The charged population is a 1e-9 seed held at its
    steady state by the external source (S = nu_wall * n_e); without it the seed decays
    at nu_wall >> nu_m through the integrator's tolerance and the accepted-state
    electron guard stops the run. Ar* -> Ar leaves the neutral total unchanged, so nu_m
    is constant along the trajectory and the exponential is exact to O(1e-7)."""
    probe, core = _neutral_diffusion_reactor(x_ion=1.0e-9)
    y0 = np.array(probe.y0[:probe.num_core_species], float)
    v0 = probe.compute_volume(y0)
    n_e = y0[0] * constants.Na / v0
    r, core = _neutral_diffusion_reactor(x_ion=1.0e-9,
                                         source=probe.compute_nu_wall(y0, v0) * n_e)
    y0 = np.array(r.y0[:r.num_core_species], float)
    idx = {s.label: i for i, s in enumerate(core)}
    nu_m = r.compute_neutral_wall_frequencies(y0, r.compute_volume(y0))[idx['Ar*']]
    r.termination = [TerminationTime((1.0 / nu_m, 's'))]
    _simulate(r, core, [])
    y = np.array(r.y[:r.num_core_species], float)
    expected = y0[idx['Ar*']] * np.exp(-nu_m * r.t)
    assert abs(y[idx['Ar*']] / expected - 1.0) < 1e-4, (y[idx['Ar*']], expected, r.t)
    gained = y[idx['Ar']] - y0[idx['Ar']]
    lost = y0[idx['Ar*']] - y[idx['Ar*']]
    ion_lost = y0[idx['Ar+']] - y[idx['Ar+']]
    assert abs((gained - ion_lost) / lost - 1.0) < 1e-9


# ---- I-269 round 2: the declaration must be a same-nuclei, strictly downhill map ----

def test_neutral_wall_refuses_isotope_transmuting_product():
    """HIGH 1: element counts cannot tell 13C-DME from DME, so the wall would turn a
    13C nucleus into a 12C one. Source and product must share the skeleton key
    (InChI minus charge layers), which keeps the isotope layer."""
    e, dme12, dme13, dmep13 = _dme_isotope_species()
    dme13.thermo = _thermo_with_h298(-100.0)       # above DME, so only the nuclei differ
    imf = {e: 1.0e-6, dmep13: 1.0e-6, dme13: 0.1, dme12: 0.9 - 2.0e-6}
    reactor = PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf,
                            (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                            diffusion_length=(_diffusion_length(), 'm'),
                            ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                            wall_recycling=1.0, wall_single_bath_approximation=True,
                            wall_neutralization_products={'DME+-13C': 'DME-13C'},
                            wall_neutral_diffusion=_meta_declaration(label='DME-13C',
                                                                     product='DME'))
    with pytest.raises(PlasmaStateError) as exc:
        reactor.initialize_model([e, dme12, dme13, dmep13], [], [], [])
    assert 'nuclei' in str(exc.value), str(exc.value)


@pytest.mark.parametrize('case, kwargs, fragment', [
    ("ground declared to pump up to the metastable",
     dict(declaration=_meta_declaration(label='Ar', product='Ar*')), 'not strictly higher'),
    ("degenerate source and product",
     dict(meta_eV=0.0), 'not strictly higher'),
    ("source thermo absent", dict(meta_eV=None), 'thermo'),
    ("product thermo absent", dict(ground_eV=None), 'thermo'),
])
def test_neutral_wall_refuses_a_declaration_that_is_not_strictly_downhill(case, kwargs, fragment):
    """HIGH 2: the wall only de-excites. A declaration whose source is not strictly
    above its product in H298 is refused, and so is one whose ordering cannot be
    established because a thermo is missing -- no default, no inference."""
    with pytest.raises(PlasmaStateError) as exc:
        _neutral_diffusion_reactor(**kwargs)
    assert fragment in str(exc.value), (case, str(exc.value))


def test_neutral_wall_refuses_a_cyclic_declaration_at_construction():
    """HIGH 2: Ar* -> Ar together with Ar -> Ar* is a cycle; it is refused by label at
    construction, before any species or thermo exists."""
    declaration = {'Ar*': {'product': 'Ar', 'diffusivity': (DP_AR_META, 'cm^2*torr/s')},
                   'Ar': {'product': 'Ar*', 'diffusivity': (DP_AR_META, 'cm^2*torr/s')}}
    electron, ar, arp = _argon_species()
    with pytest.raises(PlasmaStateError) as exc:
        PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), {ar: 1.0}, (TE_NOMINAL_EV * EV_TO_K, 'K'),
                      n_sims=1, termination=[],
                      diffusion_length=(_diffusion_length(), 'm'),
                      ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                      wall_neutral_diffusion=declaration)
    assert 'cycle' in str(exc.value), str(exc.value)


def test_neutral_wall_declared_label_absent_from_core_warns_once(caplog):
    """MEDIUM: a declared label not (yet) in the core stays tolerated, for core growth,
    but is reported at WARNING, once per label, across repeated model initialisation."""
    declaration = dict(_meta_declaration())
    declaration['Kr*'] = {'product': 'Kr', 'diffusivity': (DP_AR_META, 'cm^2*torr/s')}
    with caplog.at_level(logging.WARNING):
        r, core = _neutral_diffusion_reactor(declaration=declaration)
        r.initialize_model(core, [], [], [])
    hits = [rec for rec in caplog.records
            if rec.levelno >= logging.WARNING and "'Kr*'" in rec.getMessage()]
    assert len(hits) == 1, [rec.getMessage() for rec in caplog.records]


def test_neutral_wall_frequency_at_the_deck_value_47_torr_cm2_per_s():
    """LOW: D*p = 47 cm^2 Torr/s (Wieme & Lenaerts at 300 K), 5 Torr, 300 K, R = 5 cm,
    L = 30 cm. By hand: D = 47e-4/5 m^2/s = 9.4e-4 m^2/s; 1/Lambda^2 = (2.405/0.05)^2
    + (pi/0.30)^2 = 2313.6100 + 109.6623 = 2423.2723 m^-2; nu_m = 2.2778759 s^-1."""
    hand = (47.0e-4 / 5.0) * ((2.405 / 0.05) ** 2 + (np.pi / 0.30) ** 2)
    assert abs(hand - 2.2778759348558) < 1e-12
    r, core = _neutral_diffusion_reactor(
        declaration=_meta_declaration(diffusivity=(47.0, 'cm^2*torr/s')), tgas=300.0,
        pressure=5.0 * 101325.0 / 760.0)
    y, idx = _meta_state(r, core, n_ion=0.0)       # no charge, so n_neutral = p/(k_B T)
    nu = r.compute_neutral_wall_frequencies(y, r.compute_volume(y))[idx['Ar*']]
    assert abs(nu / hand - 1.0) < 1e-12, (nu, hand)


# ---- I-279: gas-temperature laws for the two wall-transport coefficients ----
#
# Both laws are opt-in declarations. Absent, the ion reduced mobility mu0*N_L and the
# metastable D*p are held at their declared values at every Tg, exactly as before; the
# tests below pin that the keys-absent arithmetic is today's, operation for operation.
#   ion:        mu0*N_L -> mu0*N_L * (Tg/T_ref)^m       (mobilityReferenceTemperature,
#                                                       mobilityTemperatureExponent)
#   ion temp.:  D_a     -> D_a * (1 + Tg/Te)           (ambipolarIonTemperature='gas')
#   metastable: D*p     -> D*p * (Tg/T_ref)^m          ('referenceTemperature',
#                                                       'temperatureExponent')

TGAS_HOT = 1000.0                         # K
MU0_T_REF = 300.0                         # K, the mobility law's reference temperature
MU0_T_EXP = -0.35
DP_T_REF = 302.2                          # K, where Wieme & Lenaerts give D*p = 47
DP_T_EXP = 1.68
KB_ENGINE = constants.R / constants.Na    # the Boltzmann constant the engine's EOS uses


def _mobility_law(t_ref=MU0_T_REF, m=MU0_T_EXP):
    return dict(mobility_reference_temperature=(t_ref, 'K'), mobility_temperature_exponent=m)


def _meta_law_declaration(t_ref=DP_T_REF, m=DP_T_EXP, label='Ar*'):
    decl = _meta_declaration(label=label)
    decl[label]['referenceTemperature'] = (t_ref, 'K')
    decl[label]['temperatureExponent'] = m
    return decl


def _ion_law_reactor(tgas=TGAS, te_ev=TE_NOMINAL_EV, declaration=None, **ion_kwargs):
    """Ground Ar, Ar*, Ar+, e- on the nominal wall, with optional ion-law keywords."""
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ground = _ground_species('Ar', 0.0)
    meta = _metastable_species('Ar*', AR_META_EV)
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    imf = {electron: 1.0e-6, arp: 1.0e-6, meta: 1.0e-3, ground: 1.0 - 2.0e-6 - 1.0e-3}
    kwargs = dict(diffusion_length=(_diffusion_length(), 'm'),
                  ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'),
                  wall_recycling=1.0, wall_neutralization_products={'Ar+': 'Ar'})
    if declaration is not None:
        kwargs['wall_neutral_diffusion'] = declaration
    kwargs.update(ion_kwargs)
    reactor = PlasmaReactor((tgas, 'K'), (P_NOMINAL, 'Pa'), imf, (te_ev * EV_TO_K, 'K'),
                            n_sims=1, termination=[], **kwargs)
    core = [electron, ground, meta, arp]
    reactor.initialize_model(core, [], [], [])
    return reactor, core


def _engine_n_neutral(reactor, y, V):
    """n_neutral exactly as compute_nu_wall forms it: the same summation order, the same
    Avogadro product, the same floor."""
    y_neutral = 0.0
    for j in range(reactor.num_core_species):
        if reactor.neutral_heavy_mask[j]:
            y_neutral += y[j]
    n = y_neutral * constants.Na / V
    return n if n > reactor.wall_neutral_density_floor else reactor.wall_neutral_density_floor


def test_tgas_laws_declared_at_the_gas_temperature_are_bit_identical_to_no_law():
    """I-279 test 1a: both laws declared with T_ref equal to Tg give (Tg/T_ref)^m =
    pow(1.0, m) = 1.0 exactly, so nu_wall and nu_m are == to the reactor with no law."""
    decl = _meta_law_declaration(t_ref=TGAS)
    r_law, core = _ion_law_reactor(declaration=decl, **_mobility_law(t_ref=TGAS))
    r_off, _ = _ion_law_reactor(declaration=_meta_declaration())
    y, _ = _meta_state(r_law, core, n_ion=1.0e-6)
    V = r_law.compute_volume(y)
    assert r_law.compute_nu_wall(y, V) == r_off.compute_nu_wall(y, V)
    assert np.array_equal(r_law.compute_neutral_wall_frequencies(y, V),
                          r_off.compute_neutral_wall_frequencies(y, V))


def test_tgas_laws_absent_keep_todays_arithmetic_at_1000_k():
    """I-279 test 1b: at Tg = 1000 K with no law declared, nu_wall and nu_m are == to
    today's expressions, re-evaluated here in the engine's own operation order."""
    r, core = _ion_law_reactor(tgas=TGAS_HOT, declaration=_meta_declaration())
    y, idx = _meta_state(r, core, n_ion=1.0e-6)
    V = r.compute_volume(y)
    n = _engine_n_neutral(r, y, V)
    lam = r.diffusion_length.value_si
    mu_i = MU0_AR_IN_AR * PLASMA_LOSCHMIDT / n
    d_a = mu_i * (constants.R / constants.Na) * r.Te.value_si / constants.e
    assert r.compute_nu_wall(y, V) == d_a / (lam * lam)
    dn = (DP_AR_META * CM2_TORR_TO_SI) / ((constants.R / constants.Na) * TGAS_HOT)
    assert r.wall_neutral_dn_by_label['Ar*'] == dn
    assert r.compute_neutral_wall_frequencies(y, V)[idx['Ar*']] == dn * (1.0 / (n * lam * lam))


def test_ion_mobility_law_at_1000_k_and_5_torr():
    """I-279 test 2: with mobilityReferenceTemperature = 300 K and exponent -0.35, the
    ion mobility at 1000 K is mu0*N_L/n * (1000/300)^-0.35, at a state with no charged
    population where n = p/(k_B T) exactly. Red before I-279: the keyword is unknown."""
    r, core = _ion_law_reactor(tgas=TGAS_HOT, **_mobility_law())
    y, _ = _meta_state(r, core, n_ion=0.0)
    n = P_NOMINAL / (KB_ENGINE * TGAS_HOT)
    lam = _diffusion_length()
    hand = (MU0_AR_IN_AR * PLASMA_LOSCHMIDT / n * (TGAS_HOT / MU0_T_REF) ** MU0_T_EXP
            * KB_ENGINE * (TE_NOMINAL_EV * EV_TO_K) / constants.e / (lam * lam))
    nu = r.compute_nu_wall(y, r.compute_volume(y))
    assert abs(nu / hand - 1.0) < 1e-12, (nu, hand)
    # and the law actually moved it: (1000/300)^-0.35 = 0.65613
    assert abs(nu / _nu_wall_closed_form(TE_NOMINAL_EV, n, lam) - 0.65613) < 1e-5


def test_metastable_diffusion_law_at_1000_k_and_5_torr():
    """I-279 test 3: with referenceTemperature = 302.2 K and exponent 1.68, D*p = 47 at
    302.2 K becomes 47*(1000/302.2)^1.68 at 1000 K, so nu_m = 47e-4*(1000/302.2)^1.68/5
    /Lambda^2 = 17.0073 s^-1. Red before I-279: the entry key is refused as unknown."""
    lam = _diffusion_length()
    hand = 47.0e-4 * (TGAS_HOT / DP_T_REF) ** DP_T_EXP / 5.0 / (lam * lam)
    assert abs(hand - 17.0073) < 1e-4, hand
    r, core = _neutral_diffusion_reactor(declaration=_meta_law_declaration(), tgas=TGAS_HOT)
    y, idx = _meta_state(r, core, n_ion=0.0)
    nu = r.compute_neutral_wall_frequencies(y, r.compute_volume(y))[idx['Ar*']]
    assert abs(nu / hand - 1.0) < 1e-12, (nu, hand)


def test_metastable_diffusion_law_in_dn_form_scales_dn_itself():
    """For a D*N declaration the exponent is that of D*N itself (m_p - 1): D*N quoted at
    T_ref with exponent 0.68 is the same law as D*p at T_ref with exponent 1.68."""
    dn_ref = DP_AR_META * CM2_TORR_TO_SI / (KB_ENGINE * DP_T_REF)
    decl = {'Ar*': {'product': 'Ar', 'diffusivity': (dn_ref, '1/(m*s)'),
                    'referenceTemperature': (DP_T_REF, 'K'), 'temperatureExponent': DP_T_EXP - 1.0}}
    r_dn, core = _neutral_diffusion_reactor(declaration=decl, tgas=TGAS_HOT)
    r_dp, _ = _neutral_diffusion_reactor(declaration=_meta_law_declaration(), tgas=TGAS_HOT)
    assert abs(r_dn.wall_neutral_dn_by_label['Ar*'] / r_dp.wall_neutral_dn_by_label['Ar*']
               - 1.0) < 1e-12


@pytest.mark.parametrize('te_ev', [1.0, TE_NOMINAL_EV])
def test_ambipolar_ion_temperature_gas_multiplies_nu_wall_by_one_plus_tg_over_te(te_ev):
    """I-279 test 4: ambipolarIonTemperature='gas' restores the (1 + Ti/Te) factor with
    Ti = Tg, so nu_wall/nu_wall(off) = 1 + Tg/Te, at two electron temperatures (to
    rounding: the factor multiplies D_a before the division by Lambda^2)."""
    r_on, core = _ion_law_reactor(tgas=TGAS_HOT, te_ev=te_ev, ambipolar_ion_temperature='gas')
    r_off, _ = _ion_law_reactor(tgas=TGAS_HOT, te_ev=te_ev)
    y, _ = _meta_state(r_on, core, n_ion=1.0e-6)
    V = r_on.compute_volume(y)
    ratio = r_on.compute_nu_wall(y, V) / r_off.compute_nu_wall(y, V)
    expected = 1.0 + TGAS_HOT / (te_ev * EV_TO_K)
    assert abs(ratio / expected - 1.0) < 1e-15, (ratio, expected)


def test_ambipolar_ion_temperature_rereads_te_on_every_call():
    """The (1 + Tg/Te) factor follows the CURRENT Te, which the energy balance solves;
    it is not frozen at the constructor's value."""
    r, core = _ion_law_reactor(tgas=TGAS_HOT, ambipolar_ion_temperature='gas')
    y, _ = _meta_state(r, core, n_ion=1.0e-6)
    V = r.compute_volume(y)
    te0 = r.Te.value_si
    nu0 = r.compute_nu_wall(y, V)
    r.Te.value_si = 2.0 * te0
    nu1 = r.compute_nu_wall(y, V)
    expected = 2.0 * (1.0 + TGAS_HOT / (2.0 * te0)) / (1.0 + TGAS_HOT / te0)
    assert abs(nu1 / nu0 / expected - 1.0) < 1e-14


@pytest.mark.parametrize('case, ion_kwargs, fragment', [
    ("half ion law, T_ref only", dict(mobility_reference_temperature=(300.0, 'K')),
     'mobility_temperature_exponent'),
    ("half ion law, exponent only", dict(mobility_temperature_exponent=-0.35),
     'mobility_reference_temperature'),
    ("dimensioned ion exponent", _mobility_law(m=(-0.35, 'K')), 'mobility_temperature_exponent'),
    ("bool ion exponent", _mobility_law(m=True), 'mobility_temperature_exponent'),
    ("non-finite ion exponent", _mobility_law(m=float('nan')), 'mobility_temperature_exponent'),
    ("ion exponent too large for a float", _mobility_law(m=10 ** 400), 'mobility_temperature_exponent'),
    ("zero ion T_ref", _mobility_law(t_ref=0.0), 'mobility_reference_temperature'),
    ("negative ion T_ref", _mobility_law(t_ref=-300.0), 'mobility_reference_temperature'),
    ("ion T_ref not a temperature", dict(mobility_reference_temperature=(300.0, 'm'),
                                         mobility_temperature_exponent=-0.35),
     'mobility_reference_temperature'),
    ("unknown ambipolar value", dict(ambipolar_ion_temperature='ion'), 'ambipolar_ion_temperature'),
    ("boolean ambipolar value", dict(ambipolar_ion_temperature=True), 'ambipolar_ion_temperature'),
])
def test_ion_tgas_law_refusals(case, ion_kwargs, fragment):
    """I-279 test 5 (ion, reactor level): half a law, a dimensioned or non-numeric
    exponent, T_ref <= 0 or not a temperature, and an unknown ambipolarIonTemperature
    value are each refused, naming the parameter."""
    with pytest.raises(PlasmaStateError) as exc:
        _ion_law_reactor(**ion_kwargs)
    assert fragment in str(exc.value), (case, str(exc.value))


@pytest.mark.parametrize('case, entry_update, fragment', [
    ("half law, referenceTemperature only", {'referenceTemperature': (DP_T_REF, 'K')},
     'temperatureExponent'),
    ("half law, temperatureExponent only", {'temperatureExponent': DP_T_EXP},
     'referenceTemperature'),
    ("dimensioned exponent", {'referenceTemperature': (DP_T_REF, 'K'),
                              'temperatureExponent': (1.68, 'K')}, 'temperatureExponent'),
    ("zero T_ref", {'referenceTemperature': (0.0, 'K'), 'temperatureExponent': DP_T_EXP},
     'referenceTemperature'),
    ("negative T_ref", {'referenceTemperature': (-302.2, 'K'), 'temperatureExponent': DP_T_EXP},
     'referenceTemperature'),
    ("T_ref not a temperature", {'referenceTemperature': (302.2, 's'),
                                 'temperatureExponent': DP_T_EXP}, 'referenceTemperature'),
    ("bare float T_ref", {'referenceTemperature': DP_T_REF, 'temperatureExponent': DP_T_EXP},
     'referenceTemperature'),
    ("bare int T_ref", {'referenceTemperature': 302, 'temperatureExponent': DP_T_EXP},
     'referenceTemperature'),
    ("exponent too large for a float", {'referenceTemperature': (DP_T_REF, 'K'),
                                        'temperatureExponent': 10 ** 400}, 'temperatureExponent'),
])
def test_metastable_tgas_law_refusals(case, entry_update, fragment):
    """I-279 test 5 (metastable): half a law, a dimensioned exponent and T_ref <= 0 or
    not a temperature are each refused, naming the entry key -- by the law's own
    validation, not by the unknown-key refusal that named the key before I-279."""
    decl = _meta_declaration()
    decl['Ar*'].update(entry_update)
    with pytest.raises(PlasmaStateError) as exc:
        _neutral_diffusion_reactor(declaration=decl)
    msg = str(exc.value)
    assert fragment in msg and 'Ar*' in msg and 'unknown' not in msg, (case, msg)


def test_ion_overflow_guard_refuses_a_subnormal_worst_case_from_the_law_factor():
    """I-279: the mobility law's factor (10/1)^-10 = 1e-10 is itself a normal number, but
    with mu0 = 1e-280 and Lambda = 1.7e15 m it takes the worst-case nu_wall (at the
    neutral-density floor) from a normal ~1e-300 to a positive SUBNORMAL ~1e-310, which
    the construction guard must refuse by the same finite-normal-positive predicate the
    other wall inputs use. The same deck without the law is admitted (control). The
    smallness comes from Lambda^2, not mu0: in the guard's operation order a tiny mu_i
    underflows to exactly 0 at the k_B product, which was always refused."""
    mu0, tg, lam = 1.0e-280, TGAS, 1.7e15
    law = dict(mobility_reference_temperature=(tg / 10.0, 'K'), mobility_temperature_exponent=-10.0)
    geometry = dict(ion_reduced_mobility=(mu0, 'm^2/(V*s)'), diffusion_length=(lam, 'm'))
    # The guard's expression, in its operation order, as a precondition on both arms.
    worst_n = rmgpy.solver.plasma.PLASMA_NEUTRAL_DENSITY_FLOOR_FRACTION * PLASMA_LOSCHMIDT
    tiny = np.finfo(np.float64).tiny
    for factor, subnormal in ((1.0, False), ((tg / (tg / 10.0)) ** -10.0, True)):
        worst_mu = mu0 * PLASMA_LOSCHMIDT / worst_n * factor
        worst_nu = (worst_mu * KB_ENGINE * (TE_NOMINAL_EV * EV_TO_K) / constants.e) / (lam * lam)
        assert (0.0 < worst_nu < tiny) if subnormal else (worst_nu >= tiny), (factor, worst_nu)
    _ion_law_reactor(tgas=tg, **geometry)
    with pytest.raises(PlasmaStateError) as exc:
        _ion_law_reactor(tgas=tg, **geometry, **law)
    assert 'wall loss frequency' in str(exc.value), str(exc.value)


@pytest.mark.parametrize('name, value', [
    ('mobility_reference_temperature', (300.0, 'K')),
    ('mobility_temperature_exponent', -0.35),
    ('ambipolar_ion_temperature', 'gas'),
])
def test_ion_tgas_keys_without_a_wall_are_refused(name, value):
    """I-279 test 5: each ion key given to a wall-less reactor is refused by name, at the
    reactor and at the input directive, rather than stored and silently ignored."""
    from rmgpy.exceptions import InputError
    from rmgpy.rmg.input import _plasma_wall_kwargs
    electron, ar, arp = _argon_species()
    with pytest.raises(PlasmaStateError) as exc:
        PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), {ar: 1.0}, (TE_NOMINAL_EV * EV_TO_K, 'K'),
                      n_sims=1, termination=[], **{name: value})
    assert name in str(exc.value)
    deck_name = {'mobility_reference_temperature': 'mobilityReferenceTemperature',
                 'mobility_temperature_exponent': 'mobilityTemperatureExponent',
                 'ambipolar_ion_temperature': 'ambipolarIonTemperature'}[name]
    with pytest.raises(InputError) as exc:
        _plasma_wall_kwargs(None, None, None, 1.0, None, None, **{deck_name: value})
    assert deck_name in str(exc.value)


@pytest.mark.parametrize('case, deck_kwargs, fragment', [
    ("half law", dict(mobilityReferenceTemperature=(300.0, 'K')), 'mobilityTemperatureExponent'),
    ("dimensioned exponent", dict(mobilityReferenceTemperature=(300.0, 'K'),
                                  mobilityTemperatureExponent=(-0.35, 'K')),
     'mobilityTemperatureExponent'),
    ("zero T_ref", dict(mobilityReferenceTemperature=(0.0, 'K'), mobilityTemperatureExponent=-0.35),
     'mobilityReferenceTemperature'),
    ("T_ref in eV", dict(mobilityReferenceTemperature=(0.025, 'eV'),
                         mobilityTemperatureExponent=-0.35), 'mobilityReferenceTemperature'),
    ("exponent too large for a float", dict(mobilityReferenceTemperature=(300.0, 'K'),
                                            mobilityTemperatureExponent=10 ** 400),
     'mobilityTemperatureExponent'),
    ("unknown ambipolar value", dict(ambipolarIonTemperature='ion'), 'ambipolarIonTemperature'),
])
def test_ion_tgas_law_refusals_at_the_input_directive(case, deck_kwargs, fragment):
    """I-279 test 5 (ion, input level): the deck keywords are validated while the input
    file is read, naming the keyword."""
    from rmgpy.exceptions import InputError
    from rmgpy.rmg.input import _plasma_wall_kwargs
    with pytest.raises(InputError) as exc:
        _plasma_wall_kwargs({'diffusionLength': (2.03, 'cm')}, (MU0_AR_IN_AR, 'm^2/(V*s)'),
                            None, 1.0, None, None, **deck_kwargs)
    assert fragment in str(exc.value), (case, str(exc.value))


def test_ion_tgas_keys_resolve_at_the_input_directive():
    """The input directive converts the deck keywords to the reactor's arguments, the
    reference temperature in K, and emits none of them when they are absent."""
    from rmgpy.rmg.input import _plasma_wall_kwargs
    kwargs = _plasma_wall_kwargs({'diffusionLength': (2.03, 'cm')}, (MU0_AR_IN_AR, 'm^2/(V*s)'),
                                 None, 1.0, None, None,
                                 mobilityReferenceTemperature=(300.0, 'K'),
                                 mobilityTemperatureExponent=-0.35,
                                 ambipolarIonTemperature='gas')
    assert kwargs['mobility_reference_temperature'] == 300.0
    assert kwargs['mobility_temperature_exponent'] == -0.35
    assert kwargs['ambipolar_ion_temperature'] == 'gas'
    bare = _plasma_wall_kwargs({'diffusionLength': (2.03, 'cm')}, (MU0_AR_IN_AR, 'm^2/(V*s)'),
                               None, 1.0, None, None)
    for key in ('mobility_reference_temperature', 'mobility_temperature_exponent',
                'ambipolar_ion_temperature'):
        assert key not in bare


def test_tgas_laws_survive_deepcopy_pickle_and_the_input_writer():
    """I-279 test 6: __reduce__ (so deepcopy and pickle) and _format_plasma_wall carry
    all three ion keys and both new neutral-diffusion keys. The neutral half fails if
    the stored declaration is rebuilt from 'product' and 'diffusivity' alone."""
    from rmgpy.rmg.input import _format_plasma_wall
    r, core = _ion_law_reactor(tgas=TGAS_HOT, declaration=_meta_law_declaration(),
                               ambipolar_ion_temperature='gas', **_mobility_law())
    assert r.mobility_T_factor != 1.0
    for clone in (copy.deepcopy(r), pickle.loads(pickle.dumps(r))):
        assert clone.mobility_reference_temperature == MU0_T_REF
        assert clone.mobility_temperature_exponent == MU0_T_EXP
        assert clone.ambipolar_ion_temperature == 'gas'
        assert clone.mobility_T_factor == r.mobility_T_factor
        entry = clone.wall_neutral_diffusion['Ar*']
        assert entry['referenceTemperature'] == (DP_T_REF, 'K')
        assert entry['temperatureExponent'] == DP_T_EXP
        assert clone.wall_neutral_dn_by_label == r.wall_neutral_dn_by_label
    text = _format_plasma_wall(r)
    assert 'mobilityReferenceTemperature = (300.0,"K")' in text, text
    assert 'mobilityTemperatureExponent = -0.35' in text, text
    assert "ambipolarIonTemperature = 'gas'" in text, text
    assert "'referenceTemperature': (302.2, 'K')" in text, text
    assert "'temperatureExponent': 1.68" in text, text


def test_undeclared_tgas_laws_reconstruct_as_none():
    """I-279 test 6: a wall reactor with no mobility law, and a wall-less reactor, both
    reconstruct with None for the three ion keys -- no internal default or sentinel
    leaks into __reduce__ (a wall-less copy would otherwise trip the no-wall refusal) --
    and the writer emits none of the new keywords."""
    from rmgpy.rmg.input import _format_plasma_wall
    r_wall, _ = _ion_law_reactor(declaration=_meta_declaration())
    electron, ar, arp = _argon_species()
    r_bare = PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), {ar: 1.0},
                           (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[])
    for r in (r_wall, r_bare):
        args = r.__reduce__()[1]
        assert args[-3:] == (None, None, None), args[-3:]
        clone = copy.deepcopy(r)
        assert clone.mobility_reference_temperature is None
        assert clone.mobility_temperature_exponent is None
        assert clone.ambipolar_ion_temperature is None
        text = _format_plasma_wall(r)
        for keyword in ('mobilityReferenceTemperature', 'mobilityTemperatureExponent',
                        'ambipolarIonTemperature', 'referenceTemperature', 'temperatureExponent'):
            assert keyword not in text
