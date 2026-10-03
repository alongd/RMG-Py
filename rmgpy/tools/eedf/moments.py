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

"""Independent moments from f0 using the pinned solver's cell quadrature."""
import numpy as np
from scipy.constants import elementary_charge, electron_mass, Boltzmann
from rmgpy.tools.eedf.schema import EEDFError


def distribution_moments(row, channels):
    """Integrate rates and collision powers using stored collision sections.

    The table format uses a fixed uniform grid. Reverse channels require
    explicitly resolved target/product statistical-weight ratios.
    Source: CrossSection.cpp:80–116, EedfCollisions.cpp:445–482,
    Operators.cpp:178–185,238–267. SI gamma=sqrt(2 e/m_e).
    """
    energy = np.asarray(row['energy_eV'])
    edges = np.asarray(row['energy_edges_eV'])
    widths = np.diff(edges)
    f0 = np.asarray(row['f0'])
    gamma = np.sqrt(2 * elementary_charge / electron_mass)
    result = {'mean_energy_eV': float(np.sum(f0 * energy**1.5 * widths)),
              'k_ine': np.zeros(len(channels)), 'k_sup': np.zeros(len(channels)),
              'channel_power': np.zeros(len(channels))}
    for j, channel in enumerate(channels):
        if 'cross_section' not in channel:
            raise EEDFError('cross section missing for EEDF consistency: ' + channel['description'])
        section = channel['cross_section']
        sigma = np.interp(edges, section['energy_eV'], section['sigma_m2'], left=0, right=0)
        threshold = channel['threshold_eV']
        if channel['kind'] != 'elastic':
            sigma[edges <= threshold] = 0
        cells = (sigma[:-1]+sigma[1:])/2
        lmin = min(int(threshold / widths[0]), len(energy))
        cells[:lmin] = 0
        kernel = gamma * energy * cells * widths
        result['k_ine'][j] = kernel.dot(f0)
        if '<->' in channel['description']:
            if 'statistical_weight_ratio' not in channel:
                raise EEDFError('reverse cross-section statistical weights not resolved')
            ratio = channel['statistical_weight_ratio']
            result['k_sup'][j] = (gamma * ratio * widths[0]
                * np.dot(cells[lmin:] * energy[lmin:], f0[:len(energy)-lmin]))
        target = row['target_fractions'][j]
        if channel['kind'] == 'elastic':
            if 'mass_ratio' not in channel:
                raise EEDFError('elastic mass ratio not resolved')
            g = 2 * edges**2 * sigma * channel['mass_ratio'] * target
            g[[0, -1]] = 0
            Tg = row['gas_temperature_K']
            # Tg is supplied by compare_row from the explicit physical state.
            thermal = Boltzmann / elementary_charge * Tg
            power_kernel = -gamma * (g[1:]*(thermal-widths/2)-g[:-1]*(thermal+widths/2))
            result['channel_power'][j] = power_kernel.dot(f0)
        elif channel['kind'] == 'attachment':
            result['channel_power'][j] = target * (kernel * energy).dot(f0)
        elif channel['kind'] == 'ionization':
            opb = channel.get('opb_eV', threshold)
            if opb <= 0:
                opb = threshold
            auxiliary = np.cumsum(1 / (1 + (energy / opb)**2))
            tics = np.zeros(len(energy))
            for k in range(1, len(energy)):
                kmax = int((k + 1 - lmin) / 2)
                if kmax > 0:
                    tics[k] = cells[k] * auxiliary[kmax - 1] / (opb * np.arctan((energy[k] - threshold) / (2 * opb)))
            result['channel_power'][j] = (gamma * target * energy[lmin] * widths[0]**2
                                         * np.dot(f0, energy * tics)) if lmin < len(energy) else 0.
        else:
            # LoKI represents excitation energy losses at the lower threshold face.
            effective_threshold = edges[lmin] if lmin < len(energy) else threshold
            result['channel_power'][j] = (target * result['k_ine'][j] - row['product_fractions'][j] * result['k_sup'][j]) * effective_threshold
    return result
