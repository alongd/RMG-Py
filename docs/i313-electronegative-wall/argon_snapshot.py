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

"""Capture the unchanged map-mode argon wall API for base/tip byte comparison.

This is a synthetic no-chemistry deck with solved electron energy. It does not
exercise either new closure (neither was implemented after qualification failed).
"""

import argparse
import json
import math
import struct
from pathlib import Path

import numpy as np

import rmgpy.constants as constants
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.species import Species
from rmgpy.thermo import ThermoData

from verify_mapping import assert_tested_tree


def encode(value):
    """Represent each floating value exactly, including diagnostic NaNs."""
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [encode(item) for item in value]
    if isinstance(value, (float, np.floating)):
        return 'f64:' + struct.pack('>d', float(value)).hex()
    if isinstance(value, np.integer):
        return int(value)
    return value


def snapshot(directory):
    assert_tested_tree()
    directory.mkdir(parents=True, exist_ok=True)
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    argon = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    cation = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    cation.thermo = ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([2.5 * constants.R] * 7, 'J/(mol*K)'),
        H298=(1520.0, 'kJ/mol'), S298=(154.8, 'J/(mol*K)'))
    radius, length = 0.05, 0.30
    volume = math.pi * radius ** 2 * length
    lam = 1.0 / math.sqrt((2.405 / radius) ** 2 + (math.pi / length) ** 2)
    te_kelvin = 1.02 * constants.e * constants.Na / constants.R
    reactor = PlasmaReactor(
        (298.15, 'K'), (5.0 * 133.322368, 'Pa'),
        {electron: 1e-7, argon: 1.0 - 2e-7, cation: 1e-7},
        (te_kelvin, 'K'), diffusion_length=(lam, 'm'),
        ion_reduced_mobilities={'Ar+': (1.535e-4, 'm^2/(V*s)')},
        thermo_source_assertions={'Ar+': 'ion'},
        electron_energy_balance={
            'absorbed_power': (1e-5, 'W'), 'chamber_volume': (volume, 'm^3'),
            'sheath': 'floating_wall', 'electron_energies': {},
            'elastic_collisions': {'Ar': {'ignore': 'isolated wall identity fixture'}},
        })
    reactor.initialize_model([electron, argon, cation], [], [], [])
    channels = {key: [] for key in ('wall_flux', 'nu_wall_latched',
                                   'energy_budget', 'trajectory')}
    for time in (1e-7, 1e-6, 1e-5, 1e-4):
        reactor.advance(time)
        channels['wall_flux'].append(encode(reactor.wall_flux))
        channels['nu_wall_latched'].append(encode(reactor.nu_wall_latched))
        channels['energy_budget'].append(encode(reactor.energy_budget))
        channels['trajectory'].append(encode({'t': reactor.t, 'y': reactor.y}))
    for name, values in channels.items():
        path = directory / (name + '.json')
        path.write_text(json.dumps(values, sort_keys=True, indent=2) + '\n')
        print(path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    snapshot(parser.parse_args().directory)
