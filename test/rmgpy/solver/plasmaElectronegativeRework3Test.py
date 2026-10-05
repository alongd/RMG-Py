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

"""High-precision EN derivatives and role-scaled trace-component regressions."""
from decimal import Decimal, localcontext
from pathlib import Path
import copy

import numpy as np
import pytest
import rmgpy

from plasmaElectronegativeWallTest import model, KR, KZ
from rmgpy import constants
from rmgpy.exceptions import ElectronegativeWallRegimeError
from rmgpy.solver.electronegative import closure_factor_gradient


@pytest.fixture(autouse=True)
def local_import():
    assert Path.cwd() in Path(rmgpy.__file__).resolve().parents
    import rmgpy.data.rmg as data_module
    saved = data_module.database
    data_module.database = None
    yield
    data_module.database = saved


def decimal_wall_energy_derivatives(reactor, state, columns):
    """150-digit constitutive reference, independent of the production Jacobian.

    Calibrate constant transport prefactors at the sample, then reevaluate the
    EOS, V/N transport, geometry factor and current balance in Decimal. The
    derivative uses central differences at 1e-20 relative steps in Decimal;
    doubling the step checks convergence without double-precision cancellation.
    """
    def dec(value):
        return Decimal.from_float(float(value))

    with localcontext() as context:
        context.prec = 150
        y = [dec(v) for v in state]
        electron = reactor.electron_index
        neutral = [j for j in range(reactor.num_core_species) if reactor.neutral_heavy_mask[j]]
        cations = [j for j in range(reactor.num_core_species)
                   if j != electron and reactor.species_charges[j] > 0]
        anions = [j for j in range(reactor.num_core_species)
                  if j != electron and reactor.species_charges[j] < 0]
        gas = dec(reactor.T.value_si)
        pressure = dec(reactor.P.value_si)
        constant = dec(constants.R)
        ion_temperature = gas if reactor.ambipolar_ion_temperature is not None else Decimal(0)
        weight = (Decimal(1) if reactor.electronegative_wall_geometry == 'fullFrequency'
                  else dec(KR)/(dec(KR)+dec(KZ)))

        def volume(values):
            heavy = sum(values[j] for j in range(reactor.num_core_species) if j != electron)
            return constant*(gas*heavy + values[-1]*values[electron])/pressure

        v0 = volume(y)
        n0 = sum(y[j] for j in neutral)
        ep = reactor.compute_ion_wall_components(state, reactor.compute_volume(state))['total_ep']
        prefactor = {j: dec(ep[j])*n0/((y[-1]+ion_temperature)*v0) for j in cations}
        mass = {j: dec(np.exp(reactor.energy_ion_sheath_factor[j]-0.5)) for j in cations}

        def energy(values):
            ne = values[electron]
            te = values[-1]
            n = sum(values[j] for j in neutral)
            v = volume(values)
            currents = {j: prefactor[j]*(te+ion_temperature)*v/n*values[j] for j in cations}
            gamma = sum(currents.values())
            mean_mass = sum(currents[j]*mass[j] for j in cations)/gamma
            plus = sum(values[j]*dec(reactor.species_charges[j]) for j in cations)
            minus = sum(values[j] for j in anions)
            factor = (Decimal(1) if reactor.electronegative_wall_model == 'o2ReferenceQualifiedUnity'
                      else weight*ne/(ne+minus) + (Decimal(1)-weight))
            phi = mean_mass.ln() + (ne/plus).ln() - factor.ln()
            return -te*(Decimal(1)+phi)*(factor*gamma/ne)/Decimal('1.5')

        result = []
        for j in columns:
            estimates = []
            for multiplier in (Decimal(1), Decimal(2)):
                step = abs(y[j])*Decimal('1e-20')*multiplier
                a, b = y.copy(), y.copy()
                a[j] += step
                b[j] -= step
                estimates.append((energy(a)-energy(b))/(2*step))
            assert abs(estimates[0]-estimates[1]) <= abs(estimates[0])*Decimal('1e-30')
            result.append(float(estimates[0]))
        return np.array(result)


@pytest.mark.parametrize('alpha', [0.5, 1.e3, 1.e12, 1.e23, 5.e23, 1.e30])
@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
@pytest.mark.parametrize('closure', ['confinedAnion', 'o2ReferenceQualifiedUnity'])
@pytest.mark.parametrize('ion_temperature', [None, 'gas'])
def test_extreme_alpha_wall_energy_columns_match_high_precision(alpha, arm, closure, ion_temperature):
    r, _, _ = model(arm=arm, energy=True, closure=closure,
                    ambipolar_ion_temperature=ion_temperature)
    y = r.y.copy()
    y[4] = 5.e-7
    y[0] = y[4]/alpha
    y[2] = y[0]+y[4]
    if y[0]/r.compute_volume(y) < 1.e-30:
        with pytest.raises(ElectronegativeWallRegimeError,match='domain.*e-'):
            r.monitor_electronegative_wall(y,0.)
        return
    r.monitor_electronegative_wall(y, 0.)
    r.kf[:] = 0.
    columns = list(range(len(y)))
    expected = decimal_wall_energy_derivatives(r, y, columns)
    wall = np.asarray(r.compute_electronegative_wall_jacobian(y))
    full = np.asarray(r.jacobian(0., y, np.zeros_like(y), 0.))
    assert np.isfinite(wall).all() and np.isfinite(full).all()
    np.testing.assert_allclose(wall[-1, columns], expected, rtol=5.e-8, atol=0.)
    np.testing.assert_allclose(full[-1, columns], expected, rtol=5.e-8, atol=0.)


@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
def test_subnormal_trace_anion_is_refused_outside_domain(arm):
    r, _, _ = model(arm=arm, energy=True)
    y = r.y.copy()
    y[4] = 1.e-320
    y[2] = y[0]+y[4]
    with pytest.raises(ElectronegativeWallRegimeError,match='domain.*Cl-'):
        r.monitor_electronegative_wall(y,7.)


@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
def test_trace_neutral_resolves_bulk_chemistry_derivative(arm):
    r, _, _ = model(arm=arm, energy=True)
    y = r.y.copy()
    y[1], y[3] = 1.e-29, 1.
    r.monitor_electronegative_wall(y, 7.)
    matrix = np.asarray(r.jacobian(7., y, np.zeros_like(y), 0.))
    estimates = []
    for step in (1.e-5, 2.e-5):
        a, b = y.copy(), y.copy()
        a[1] += step
        b[1] += 2.*step
        f0 = r.residual(7., y, np.zeros_like(y))[0][0]
        f1 = r.residual(7., a, np.zeros_like(y))[0][0]
        f2 = r.residual(7., b, np.zeros_like(y))[0][0]
        estimates.append((-3.*f0+4.*f1-f2)/(2.*step))
    assert estimates[0] == pytest.approx(estimates[1], rel=5.e-7)
    assert estimates[0] > 250.
    assert matrix[0, 1] == pytest.approx(estimates[0], rel=2.e-5)


@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
def test_small_alpha_factor_gradient_does_not_subtract_h_from_one(arm):
    r, _, _ = model(arm=arm)
    y = r.y.copy()
    y[4] = 1.e-30
    _, _, gradient = closure_factor_gradient(y, 0, [4], 'confinedAnion', arm, (KR, KZ))
    weight = 1. if arm == 'fullFrequency' else KR/(KR+KZ)
    expected = weight*y[4]/(y[0]+y[4])**2
    assert gradient[0] == pytest.approx(expected, rel=1.e-14, abs=0.)


@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
def test_neutral_wall_derivative_keeps_tiny_charged_pressure(arm):
    r, _, _ = model(arm=arm, energy=True)
    y = r.y.copy()
    y[0], y[4], y[2] = 1.e-29, 5.e-30, 1.5e-29
    expected = decimal_wall_energy_derivatives(r, y, [1, 3])
    actual = r.compute_electronegative_wall_jacobian(y)[-1, [1, 3]]
    assert np.all(expected > 0.)
    np.testing.assert_allclose(actual, expected, rtol=5.e-8, atol=0.)


def test_fixed_te_radial_only_unrepresentable_electron_frequency_refuses():
    r, _, _ = model(arm='radialOnly')
    y = r.y.copy()
    y[0], y[4], y[2] = 1.e-320, 5.e-7, 5.e-7
    old = copy.deepcopy(r.electronegative_wall_last_valid_state)
    with pytest.raises(ElectronegativeWallRegimeError, match='domain.*e-|electron wall.*frequency|wall frequency'):
        r.monitor_electronegative_wall(y, 7.)
    assert r.electronegative_wall_last_valid_state[0] == old[0]
    np.testing.assert_array_equal(r.electronegative_wall_last_valid_state[1], old[1])


@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
def test_monitor_checks_every_operator_before_publication(monkeypatch, arm):
    from plasmaElectronegativeWallTest import FixtureEnergyReactor
    calls = []
    original = FixtureEnergyReactor._energy_jacobian

    def counted(self, time, state, derivative, cj):
        calls.append(time)
        return original(self, time, state, derivative, cj)

    monkeypatch.setattr(FixtureEnergyReactor, '_energy_jacobian', counted)
    r, _, _ = model(arm=arm, energy=True)
    for time in range(1, 66):
        r.monitor_electronegative_wall(r.y.copy(), float(time))
    assert calls == list(map(float,range(66)))
    old = copy.deepcopy(r.electronegative_wall_last_valid_state)
    y = r.y.copy()
    y[0], y[4], y[2] = 1.e-300, 5.e-301, 1.5e-300
    with pytest.raises(ElectronegativeWallRegimeError, match='domain.*e-|energy Jacobian.*range'):
        r.monitor_electronegative_wall(y, 66.)
    assert calls[-1] == 65.  # The domain refusal precedes evaluating the operator.
    assert r.electronegative_wall_last_valid_state[0] == old[0]
    np.testing.assert_array_equal(r.electronegative_wall_last_valid_state[1], old[1])


@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
def test_fixed_te_neutral_wall_column_keeps_tiny_charged_pressure(arm):
    r, _, _ = model(arm=arm)
    y = r.y.copy()
    y[0], y[4], y[2] = 1.e-29, 5.e-30, 1.5e-29
    r.kf[:] = 0.
    matrix = np.asarray(r.jacobian(0., y, np.zeros_like(y), 0.))
    wall = r.compute_electronegative_wall_jacobian(y)
    charged_pressure = constants.R/r.P.value_si*(
        r.T.value_si*(y[2]+y[4])+r.Te.value_si*y[0])
    neutral = y[1]+y[3]
    loss = r.compute_ion_wall_frequencies(y, r.compute_volume(y))[2]*y[2]
    expected = loss*(charged_pressure/r.compute_volume(y))/neutral
    assert expected > 0.
    np.testing.assert_allclose(matrix[[0, 2], 1], expected, rtol=1.e-14, atol=0.)
    np.testing.assert_allclose(wall[[0, 2], 1], expected, rtol=1.e-14, atol=0.)
