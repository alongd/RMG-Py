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

"""Published Holstein limits and Walsh closed-form domain."""
import math
import pytest
from rmgpy.solver.radiation import walsh_escape_factor


@pytest.mark.parametrize('geometry, coefficient', [('slab', 1.875), ('infiniteCylinder', 1.60)])
def test_holstein_1951_doppler_normalized_value(geometry, coefficient):
    # Holstein (1951), Tables I/III, m=1, and Eqs.6.2/6.3.
    tau = 1.e6
    result = walsh_escape_factor(tau, 0., geometry)
    assert result * tau * math.sqrt(math.pi * math.log(tau)) == pytest.approx(coefficient)


@pytest.mark.parametrize('geometry, coefficient', [('slab', 1.150), ('infiniteCylinder', 1.115)])
def test_holstein_1951_pressure_normalized_value(geometry, coefficient):
    # Tables I/III, m=1/2. tau_L=tau_D/(sqrt(pi)*a).
    tau, a = 1.e8, .1
    result = walsh_escape_factor(tau, a, geometry)
    assert result * math.sqrt(math.sqrt(math.pi) * tau / a) == pytest.approx(coefficient)


@pytest.mark.parametrize('tau, a', [(0., .1), (1., .1), (9999., .1), (float('nan'), .1), (float('inf'), .1), (1e5, -.1), (1e5, .5001), (1e5, float('nan'))])
def test_refuses_outside_thick_profile_domain(tau, a):
    with pytest.raises(ValueError, match='radiation'):
        walsh_escape_factor(tau, a, 'slab')


def test_refuses_unknown_geometry():
    with pytest.raises(ValueError, match='geometry'):
        walsh_escape_factor(1e5, .1, 'sphere')


@pytest.mark.parametrize('geometry', ['slab', 'infiniteCylinder'])
@pytest.mark.parametrize('tau,a', [(1e4, 1e-9), (2e4, 1e-3), (2e5, .03), (1e6, .1)])
def test_derivatives_match_independent_finite_differences(geometry, tau, a):
    from rmgpy.solver.radiation import walsh_escape_derivatives
    # Inward tau difference at the declared boundary; central elsewhere.
    step = tau * 1e-5
    if tau == 1e4:
        d_tau = (-3*walsh_escape_factor(tau,a,geometry)
                 +4*walsh_escape_factor(tau+step,a,geometry)
                 -walsh_escape_factor(tau+2*step,a,geometry))/(2*step)
    else:
        d_tau = (walsh_escape_factor(tau+step,a,geometry)
                 -walsh_escape_factor(tau-step,a,geometry))/(2*step)
    step_a = a * 1e-4
    d_a = (walsh_escape_factor(tau,a+step_a,geometry)
           -walsh_escape_factor(tau,a-step_a,geometry))/(2*step_a)
    got_tau, got_a = walsh_escape_derivatives(tau, a, geometry)
    assert got_tau == pytest.approx(d_tau, rel=1e-5, abs=1e-18)
    assert got_a == pytest.approx(d_a, rel=1e-4, abs=1e-9)


@pytest.fixture
def nist_lines():
    from pathlib import Path
    from rmgpy import settings
    from rmgpy.solver.radiation import load_argon_lines
    return load_argon_lines(Path(settings.require_database_directory())/'eedf/argon_resonance_radiation.json')


def test_sourced_line_identity_and_vacuum_lifetime(nist_lines):
    # NIST ASD archived rows, not values calculated by the rate implementation.
    expected = {'3P1':('018001.000003','1s4',106.6659865e-9,1.320e8),
                '1P1':('018001.000005','1s2',104.8219924e-9,5.32e8)}
    assert nist_lines['query']['sha256'] == '015756d4605ff3bd24eed07c12826815d8b1ea3fcacdd1190eb6cfaa22ea36af'
    for state,(identity,paschen,wavelength,rate) in expected.items():
        row = nist_lines['lines'][state]
        assert row['upper_nist_id'] == identity
        assert row['paschen'] == paschen
        assert row['wavelength_m'] == pytest.approx(wavelength,rel=1e-12)
        assert row['A_s_1'] == rate
        assert row['g_ratio'] == 3


@pytest.mark.parametrize('temperature', [249.,501.,float('nan')])
def test_temperature_refusal(nist_lines, temperature):
    from rmgpy.solver.radiation import argon_escape
    with pytest.raises(ValueError,match='temperature'):
        argon_escape(nist_lines['lines']['3P1'],1.6e23,temperature,.05,'infiniteCylinder')


@pytest.mark.parametrize('state', ['3P1','1P1'])
def test_ground_density_derivative_includes_self_broadening(nist_lines, state):
    from rmgpy.solver.radiation import argon_escape
    line = nist_lines['lines'][state]
    density = 1.6e23
    value, slope, tau, a = argon_escape(line,density,298.15,.05,'infiniteCylinder')
    h = density*1e-4
    difference = (argon_escape(line,density+h,298.15,.05,'infiniteCylinder')[0]-
                  argon_escape(line,density-h,298.15,.05,'infiniteCylinder')[0])/(2*h)
    assert slope == pytest.approx(difference,rel=1e-5,abs=1e-35)
    assert 0 < value < .01 and tau > 1e4 and 0 < a < .5


def test_nitrogen_impact_width_omission_envelope(nist_lines):
    # Yatom et al. (2023), Eq.15, DOI 10.1088/1361-6595/ad0741:
    # FWHM(cm)=8.18e-12 lambda(cm)^2 (alpha(cm^3) R^2(a0^2))^.4
    #          * (T(K)/mu(u))^.3 N(cm^-3).
    # NIST Coulomb radius estimate gives R^2=25.17 and 28.04; use 30
    # as an envelope, alpha(N2)=1.710 A^3 (NIST CCCBDB), T<=500 K.
    # This bounds that impact estimate, not unknown molecular corrections.
    from rmgpy.solver.radiation import argon_escape
    mu = 39.948*28.014/(39.948+28.014)
    frequency_per_density = 2.99792458e10*8.18e-12*(1.710e-24*30)**.4*(500/mu)**.3
    assert frequency_per_density < 8.3e-10  # Hz cm^3
    for state in ('3P1','1P1'):
        line = nist_lines['lines'][state]
        # 20% N2 at 5 torr, deliberately at the maximum temperature.
        density_total = (5.*101325./760.)/(1.380649e-23*500.)
        g,_,tau,a = argon_escape(line,.8*density_total,500.,.05,'infiniteCylinder')
        width_d = math.sqrt(2*1.380649e-23*500/(39.948*1.66053906660e-27))/line['wavelength_m']
        a_n2 = frequency_per_density*.2*density_total/1e6/(2*width_d)
        broadened = walsh_escape_factor(tau,a+a_n2,'infiniteCylinder')
        assert 0 < broadened/g-1 < .10
