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

"""Source mapping and owner-ruling limit tests for the declared axial extension.

The finite cylinder is not required to reproduce a radial-only source in
absolute units. Source validation is radial; the no-radial synthetic limit
below explicitly tests the model extension, not the paper.
"""

import pytest
from scipy.integrate import quad
from scipy.special import j0

from verify_mapping import (BETA, CHI, FROZEN_ALPHA, RADIUS, assert_tested_tree,
                            complete_source, engine_reference)


@pytest.fixture(autouse=True)
def tested_tree():
    assert_tested_tree()


def test_center_to_volume_normalization():
    # Independent spatial integral of the source's equation (4) profile.
    average, _ = quad(lambda radius: 2.0 * radius * j0(CHI * radius / RADIUS)
                      / RADIUS ** 2, 0.0, RADIUS)
    assert average == pytest.approx(0.4317548070196805, rel=1e-13)
    assert average == pytest.approx(BETA, rel=1e-13)


@pytest.mark.parametrize('alpha,h', [(0.0, 1.0), (0.5, 2.0 / 3.0), (0.6, 0.625)])
def test_complete_radial_source_ratio(alpha, h):
    # This proves a ratio WITHIN the source; it does not prove engine mapping.
    positive = complete_source(0.05, 0.0)
    negative = complete_source(0.05, alpha)
    assert negative['nu_center_converted'] / positive['nu_center_converted'] == \
        pytest.approx(h, rel=1e-14)


@pytest.mark.parametrize('p_o2', FROZEN_ALPHA)
def test_radial_only_diffusion_limit(p_o2):
    alpha = FROZEN_ALPHA[p_o2]
    reference = complete_source(p_o2, alpha)
    # The genuinely matching geometry recovers the sourced asymptote.
    candidate = engine_reference(p_o2, radial_only=True) / (1.0 + alpha)
    assert reference['nu_center_converted'] == pytest.approx(candidate, rel=1e-6)


def active_frequency(alpha, radial, axial, arm):
    import numpy as np
    from rmgpy.solver.electronegative import closure_factor_gradient
    h, factor, derivative = closure_factor_gradient(
        np.array([1.,alpha]),0,[1],'confinedAnion',arm,(radial,axial))
    return factor*(radial+axial)


@pytest.mark.parametrize('alpha',[0.5,0.6])
@pytest.mark.parametrize('arm',['fullFrequency','radialOnly'])
def test_infinite_column_limit_recovers_radial_source(alpha,arm):
    # z goes to zero as L goes to infinity. Compare the actual production
    # factor helper to the independently derived complete sourced expression.
    source=complete_source(0.05,alpha)
    radial_ep=complete_source(0.05,0.)['nu_radial_asymptote']
    assert active_frequency(alpha,radial_ep,0.,arm) == pytest.approx(source['nu_center_converted'],rel=1.e-6)


@pytest.mark.parametrize('arm',['fullFrequency','radialOnly'])
def test_alpha_zero_preserves_complete_frequency_bits(arm):
    import struct
    radial,axial=2345.67,123.456
    assert struct.pack('d',active_frequency(0.,radial,axial,arm)) == struct.pack('d',radial+axial)


@pytest.mark.parametrize('alpha',[0.5,0.6])
def test_no_axial_loss_arms_identical(alpha):
    assert active_frequency(alpha,35.,0.,'fullFrequency') == active_frequency(alpha,35.,0.,'radialOnly')


@pytest.mark.parametrize('alpha,h',[(0.5,2./3.),(0.6,0.625)])
def test_no_radial_loss_synthetic_model_extension(alpha,h):
    # MODEL-EXTENSION TEST ONLY: the radial source has no axial-loss channel.
    assert active_frequency(alpha,0.,10.,'fullFrequency') == h*10.
    assert active_frequency(alpha,0.,10.,'radialOnly') == 10.
