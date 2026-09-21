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

import numpy as np
import pytest

import rmgpy.constants as constants
from rmgpy.exceptions import PlasmaStateError
from rmgpy.kinetics import VoronovEIArrhenius
from rmgpy.reaction import Reaction
from rmgpy.solver.plasma import PLASMA_LOSCHMIDT, PlasmaReactor
from rmgpy.species import Species

EV_TO_K = 1.0 / 8.617333262e-5           # K per eV
TORR_TO_PA = 101325.0 / 760.0

TGAS = 298.15                             # K
P_NOMINAL = 5.0 * TORR_TO_PA              # Pa
TE_NOMINAL_EV = 3.0

R_NOMINAL = 0.05                          # m, cylinder radius
L_NOMINAL = 0.30                          # m, cylinder length

# Ellis, McDaniel & Albritton, At. Data Nucl. Data Tables 17 (1976) 177.
MU0_AR_IN_AR = 1.535e-4                   # m^2/(V s) at the Loschmidt density

VORONOV_YAML = '/home/alon/Code/RMG-database-plasma/input/kinetics/voronov.yaml'


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


def test_no_wall_residual_and_jacobian_bitwise_reproducible():
    """A reactor with no wall declared produces residual and Jacobian values
    bit-for-bit identical to one whose wall parameters are absent -- built
    twice, independently, from the same equations. Uses ==, not approx."""
    ra, _, _ = _build_reactor(wall=False, with_chemistry=True)
    rb, _, _ = _build_reactor(wall=False, with_chemistry=True)
    assert ra.has_wall is False
    assert rb.has_wall is False

    y = _state_at(ra, 1.0e-6)
    dydt = np.zeros(ra.num_core_species, float)
    da, ia = ra.residual(0.0, y.copy(), dydt.copy())
    db, ib = rb.residual(0.0, y.copy(), dydt.copy())
    assert np.array_equal(da, db)
    assert ia == ib

    pa = np.array(ra.jacobian(0.0, y.copy(), dydt.copy(), 1.0), float)
    pb = np.array(rb.jacobian(0.0, y.copy(), dydt.copy(), 1.0), float)
    assert np.array_equal(pa, pb)


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
