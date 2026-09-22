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
import os
import pickle
import re

import numpy as np
import pytest

import rmgpy.constants as constants
import rmgpy.solver.plasma
from rmgpy import settings
from rmgpy.exceptions import PlasmaStateError
from rmgpy.kinetics import VoronovEIArrhenius
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

# Ellis, McDaniel & Albritton, At. Data Nucl. Data Tables 17 (1976) 177.
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
                   lam=None, mu0=MU0_AR_IN_AR, termination=None):
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
    reactor = PlasmaReactor(
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


def test_jacobian_matches_residual_below_the_neutral_floor():
    """The Jacobian must differentiate the residual the solver ACTUALLY evaluates,
    including on the clamped branch below ``wall_neutral_floor``.

    ``compute_nu_wall`` floors the neutral moles, so below the floor nu_wall is
    CONSTANT in the neutral amount and its derivative with respect to a neutral
    species is exactly zero. A Jacobian that keeps differentiating the unclamped
    1/y_neutral there reports a derivative for a dependence the residual does not
    have, and the discrepancy diverges as y_neutral -> 0.

    This samples strictly BELOW the floor on purpose. The finite-difference check
    in verify_operator.py sampled a normal state and is structurally blind to this:
    above the floor the two branches agree exactly, so any state the earlier check
    could reach reproduces the blind spot rather than closing it.
    """
    r, _, _ = _build_reactor(wall=True, with_chemistry=False)
    ie, i_ar, i_arp = _indices(r)
    floor = r.wall_neutral_floor
    assert floor > 0.0, "the floor must exist for this test to mean anything"

    # A state whose neutral content is an order of magnitude BELOW the floor,
    # so compute_nu_wall is on its clamped branch.
    y = _state_at(r, 1.0e-6)
    y[i_ar] = 0.1 * floor
    assert _wall_neutral_moles(r, y) < floor

    V = r.compute_volume(y)
    nu_here = r.compute_nu_wall(y, V)

    # Precondition: the EXPLICIT 1/y_neutral is clamped here. nu_wall is NOT flat
    # in the neutral amount below the floor -- nu ~ V/y_neutral_eff, and V keeps its
    # EOS dependence on that species -- but with y_neutral_eff pinned to the floor,
    # nu/V is identical at any two below-floor states. That ratio is the precise
    # statement of "clamped", and the clamped term is the one the Jacobian omits.
    y_half = y.copy()
    y_half[i_ar] = 0.05 * floor
    V_half = r.compute_volume(y_half)
    assert np.isclose(r.compute_nu_wall(y_half, V_half) / V_half, nu_here / V,
                      rtol=1e-12, atol=0.0), (
        "precondition failed: nu_wall is not clamped below the floor")

    zeros = np.zeros(r.num_core_species, float)
    jac = np.asarray(r.jacobian(0.0, y, zeros, 0.0), float)

    # Central finite difference of the residual with respect to the neutral,
    # over steps small enough to stay under the floor.
    col = i_ar
    best = None
    for h in (1.0e-3 * floor, 1.0e-4 * floor, 1.0e-5 * floor):
        yp = y.copy(); yp[col] += h
        ym = y.copy(); ym[col] -= h
        assert _wall_neutral_moles(r, yp) < floor, "step left the clamped branch"
        rp = np.asarray(r.residual(0.0, yp, zeros)[0], float)
        rm = np.asarray(r.residual(0.0, ym, zeros)[0], float)
        fd = (rp - rm) / (2.0 * h)
        err = np.max(np.abs(jac[:, col] - fd)) / max(1.0, np.max(np.abs(fd)))
        best = err if best is None else min(best, err)

    assert best < 1.0e-6, (
        "analytic Jacobian disagrees with the residual it claims to differentiate "
        "below the neutral floor: relative error {0:.3e}".format(best))


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

    missing = [p for p in params if not re.search(r'\bself\.%s\b' % re.escape(p), body)]
    assert not missing, (
        "__reduce__ does not carry {0} of __init__'s parameters, so a pickle or "
        "deepcopy silently drops them: {1}".format(len(missing), missing))


# ---------------------------------------------------------------- I-246 rework
#
# The ground/metastable resolution, the thermo precondition, the near-degeneracy
# threshold, the wall_neutralization_products declaration, the neutral-floor and
# electron-population validators, the single-cation enforcement, and the direct
# diffusionLength dimension check. Every red state here is banked, on the built
# module, in docs/i246-ambipolar-wall-operator/evidence/before.log.

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
    """An accepted state whose neutral inventory is at or below the numerical floor
    is refused: the ceiling bounds n_e/n_neutral, so a collapsed-inventory state
    (low alpha, tiny n_neutral) passed before this rework while nu_wall silently came
    off the clamp."""
    r, _, _ = _build_reactor(wall=True, gamma=0.0, with_chemistry=False)
    ie, i_ar, i_arp = _indices(r)
    y = np.zeros(r.num_core_species, float)
    y[i_ar] = 1.0e-9                      # below floor (~1e-6 of the initial neutral moles)
    y[i_arp] = 1.0e-15
    y[ie] = 1.0e-15                       # alpha = 1e-6, far under the ceiling
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
    reactor = PlasmaReactor((TGAS, 'K'), (P_NOMINAL, 'Pa'), imf,
                            (TE_NOMINAL_EV * EV_TO_K, 'K'), n_sims=1, termination=[],
                            diffusion_length=(_diffusion_length(), 'm'),
                            ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'))
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
                  wall_recycling=gamma)
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
                            wall_recycling=0.0, ionisation_source=(source, 'm^-3/s'))
    core = [electron, ar, he, arp]
    reactor.initialize_model(core, [], [], [])
    z = reactor.species_charges
    ie = reactor.electron_index
    i_ar = [j for j in range(len(z)) if z[j] == 0 and core[j].label == 'Ar'][0]
    i_he = [j for j in range(len(z)) if z[j] == 0 and core[j].label == 'He'][0]
    y = np.zeros(reactor.num_core_species, float)
    y[i_ar] = 0.5
    y[i_he] = 0.5
    V = reactor.compute_volume(y)
    delta, _ = reactor.residual(0.0, y.copy(), np.zeros_like(y))
    source_total = source * V / constants.Na       # mol of pairs per second
    assert np.isclose(delta[ie], source_total, rtol=1e-9, atol=0.0)


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
                            wall_recycling=1.0, ionisation_source=(1.0e5, 'm^-3/s'))
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
                            wall_recycling=1.0, ionisation_source=(1.0e18, 'm^-3/s'))
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
                            ion_reduced_mobility=(MU0_AR_IN_AR, 'm^2/(V*s)'))
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
