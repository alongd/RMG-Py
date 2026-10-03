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

import numpy as np
from rmgpy.tools.eedf.validation import split_branches


def test_two_different_scan_paths_are_not_averaged():
    def row(energy):
        return {'swarm': {'mean_energy_eV': energy}, 'k_ine': np.array([energy * 1e-15]),
                'k_sup': np.array([0.]), 'power_groups': {'field': energy * 1e-16},
                'f0': np.array([energy]), 'channel_power': np.array([energy * 1e-16]),
                'attachment_energy_eV': np.array([0.])}
    low = [row(1.), row(2.)]
    high = [row(3.), row(4.)]
    branches, verdict = split_branches({'up': low, 'down': high[::-1], 'cold': low},
                                      {'rtol': 1e-6, 'atol': 0.})
    assert len(branches) == 2
    assert not verdict['agreement']
    assert branches['branch_0'][1]['swarm']['mean_energy_eV'] == 2.
    assert branches['branch_1'][1]['swarm']['mean_energy_eV'] == 4.
