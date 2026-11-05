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

"""Declared total argon 4s population and its interpolation coordinate."""
import itertools
import numpy as np

AXIS = 'f_Ar4s'
NODES = [0., 1.e-7, 1.e-6, 1.e-5, 1.e-4]
STATES = ('Ar(3P2)', 'Ar(3P1)', 'Ar(3P0)', 'Ar(1P1)')
CONVENTION = {'axis': AXIS, 'ratio': 'total-4s/ground',
              'weights': [5, 3, 1, 3],
              'coordinate': 'linear-to-first-positive-then-log',
              'first_positive': 1.e-7}


def validate_population_axis(document):
    """Require the complete frozen convention whenever the axis is present."""
    from rmgpy.tools.eedf.schema import SpecError
    axes = document.get('axes', {})
    declaration = document.get('argon_4s_population')
    if AXIS not in axes and declaration is None:
        return
    if declaration != CONVENTION or AXIS not in axes:
        raise SpecError('argon 4s population axis requires its declared convention')
    nodes = np.asarray(axes[AXIS], dtype=float)
    if (nodes.ndim != 1 or len(nodes) < len(NODES) or not np.all(np.isfinite(nodes))
            or nodes[0] != 0. or nodes[-1] != 1.e-4
            or not set(NODES).issubset(nodes) or np.any(np.diff(nodes) <= 0)):
        raise SpecError('f_Ar4s requires the five base nodes from zero to top node 1e-4')
    # Refinement may add nodes below the original first positive node; the
    # declared coordinate transition remains fixed at 1e-7.
    if 'excited_fraction_max' in document:
        raise SpecError('f_Ar4s replaces the zero-population excited_fraction_max premise')


def populations(fraction):
    """Return normalized state populations for the declared ratio to ground."""
    from rmgpy.tools.eedf.schema import SpecError
    if not np.isfinite(fraction) or not 0. <= fraction <= NODES[-1]:
        raise SpecError('f_Ar4s outside zero to top node 1e-4')
    ground = 1. / (1. + fraction)
    return {'Ar(1S0)': ground, **{
        state: fraction * ground * weight / 12.
        for state, weight in zip(STATES, CONVENTION['weights'])}}


def coordinate(values):
    """PCHIP coordinate: f/1e-7 below 1e-7, 1+ln(f/1e-7) above."""
    values = np.asarray(values, dtype=float)
    return np.where(values <= 1.e-7, values / 1.e-7,
                    1. + np.log(np.maximum(values, 1.e-7) / 1.e-7))


def coordinate_derivative(value):
    return 1. / max(float(value), 1.e-7)


def check_zero_target_nodes(spec, channels):
    """Refuse the pinned LoKI inverse-assembly defect at any production node.

    Operators.cpp:578 nests inverse assembly inside targetDensity != 0.
    A zero target with a populated inverse source is therefore unsupported.
    """
    from rmgpy.tools.eedf.integrity import _resolved_properties, _channel_fractions
    from rmgpy.tools.eedf.schema import SpecError
    # Field does not affect prescribed populations; check every other axis.
    names = [name for name in spec['axes'] if name != 'u']
    checked = 0
    for values in itertools.product(*(spec['axes'][name] for name in names)):
        point = dict(zip(names, values))
        _, _, gas, states = _resolved_properties(spec, point)
        target, product = _channel_fractions(channels, gas, states, point)
        for channel, source, reverse_source in zip(channels, target, product):
            if source == 0. and reverse_source > 0.:
                raise SpecError('LoKI zero-population target with populated product: '
                                + channel['description'] + ' at ' + repr(point))
        checked += 1
    return checked
