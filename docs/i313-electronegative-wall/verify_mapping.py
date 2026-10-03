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

"""Reproduce the source-to-engine qualification gate, without changing transport.

Sources: Kemaneci et al., arXiv:1612.07268, equations (2), (4)-(6),
(14), (25), and Appendix A. Frozen conditions and mobility estimate:
/home/alon/runs/i311-electronegative-wall/alpha_0d.py.

This is an offline qualification tool, not the independent Gate C reference.
The transport reference state is electropositive. Source comparisons remain
radial; the owner-authorised axial multiplier is a declared geometric extension.
This driver measures the exact deck axial share and actual closure difference;
it does not qualify either independent reference interface.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import j1, jn_zeros

import rmgpy
import rmgpy.constants as constants
from rmgpy.solver.plasma import PLASMA_LOSCHMIDT, PlasmaReactor
from rmgpy.species import Species
from rmgpy.thermo import ThermoData

RADIUS = 0.05
LENGTH = 0.30
TGAS = 298.15
TE_EV = 1.02
NE = 5.2e15
TORR = 133.322368
K0 = 0.00025697516730204976  # I-311 Langevin estimate, NOT production mobility data.
CHI = float(jn_zeros(0, 1)[0])
BETA = float(2.0 * j1(CHI) / CHI)
MASS_O2 = 31.998 * 1.66053906660e-27
FROZEN_ALPHA = {
    0.05: 0.4965113464303668,
    0.10: 0.5321488783122146,
    0.25: 0.5639292236211563,
    0.50: 0.5856769812319058,
    1.00: 0.6123329596064239,
}


def assert_tested_tree():
    """Refuse evidence from an editable install pointing at another worktree."""
    root = Path(__file__).resolve().parents[2]
    resolved = Path(rmgpy.__file__).resolve()
    if root not in resolved.parents:
        raise RuntimeError('rmgpy resolved outside tested tree: {}'.format(resolved))
    return root


def complete_source(p_o2, alpha, include_ion_temperature=True):
    """Full equation (6), with the qualified equation (14) Bohm ratio equal to one.

    Return both the published global equation (25), which multiplies its h by
    a volume-average density, and equation (5) converted consistently using
    the equation (4) Bessel profile. Do not silently identify those conventions.
    """
    kb = constants.R / constants.Na
    neutral_density = (5.0 + p_o2) * TORR / (kb * TGAS)
    mobility = K0 * PLASMA_LOSCHMIDT / neutral_density
    da = mobility * (TE_EV + (kb * TGAS / constants.e
                             if include_ion_temperature else 0.0))
    bohm = math.sqrt(constants.e * TE_EV / MASS_O2)
    correction = 1.0 / (1.0 + alpha)
    edge_to_center = correction / math.sqrt(
        1.0 + (RADIUS * bohm / (CHI * da * j1(CHI))) ** 2)
    nu_published = (2.0 / RADIUS) * bohm * edge_to_center
    return {
        'diffusivity': da,
        'bohm_speed': bohm,
        'edge_to_center': edge_to_center,
        'nu_published_eq25': nu_published,
        'nu_center_converted': nu_published / BETA,
        'nu_radial_asymptote': correction * da * CHI ** 2 / RADIUS ** 2,
    }


def engine_reference(p_o2, include_ion_temperature=True, radial_only=False):
    """Measure the public map-mode API at the frozen neutral density and Te.

    Synthetic thermo supports this no-chemistry transport fixture only. Charged
    pressure is included in the EOS so it does not perturb the frozen neutral
    density. No chemistry or mobility data is added to the database.
    """
    assert_tested_tree()
    kb = constants.R / constants.Na
    te_kelvin = constants.e * TE_EV / kb
    neutral_density = (5.0 + p_o2) * TORR / (kb * TGAS)
    ne = NE
    pressure = kb * ((neutral_density + ne) * TGAS + ne * te_kelvin)
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    oxygen = Species(label='O2').from_adjacency_list(
        'multiplicity 3\n1 O u1 p2 c0 {2,S}\n2 O u1 p2 c0 {1,S}')
    cation = Species(label='O2+').from_adjacency_list(
        'multiplicity 2\n1 O u1 p1 c+1 {2,D}\n2 O u0 p2 c0 {1,D}')
    cation.thermo = ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([30.0] * 7, 'J/(mol*K)'), H298=(0.0, 'J/mol'),
        S298=(150.0, 'J/(mol*K)'))
    core = [electron, ar, oxygen, cation]
    fractions = {electron: ne, ar: neutral_density * 5.0 / (5.0 + p_o2),
                 oxygen: neutral_density * p_o2 / (5.0 + p_o2), cation: ne}
    total = sum(fractions.values())
    fractions = {species: amount / total for species, amount in fractions.items()}
    lam = (RADIUS / CHI if radial_only else
           1.0 / math.sqrt((2.405 / RADIUS) ** 2 + (math.pi / LENGTH) ** 2))
    kwargs = {'ambipolar_ion_temperature': 'gas'} if include_ion_temperature else {}
    reactor = PlasmaReactor(
        (TGAS, 'K'), (pressure, 'Pa'), fractions, (te_kelvin, 'K'),
        diffusion_length=(lam, 'm'),
        ion_reduced_mobilities={'O2+': (K0, 'm^2/(V*s)')},
        wall_recycling=0.0, wall_single_bath_approximation=True,
        thermo_source_assertions=['O2+'], **kwargs)
    reactor.initialize_model(core, [], [], [])
    y = np.asarray(reactor.y0)
    volume = reactor.compute_volume(y)
    measured_neutral = (y[1] + y[2]) * constants.Na / volume
    if not math.isclose(measured_neutral, neutral_density, rel_tol=1e-12):
        raise RuntimeError('fixture changed the frozen neutral density')
    return float(reactor.compute_ion_wall_frequencies(y, volume)[3])


def measurements():
    """Record every nominal frozen study point under both temperature choices."""
    assert_tested_tree()
    rows = []
    for include_ion_temperature in (False, True):
        for pressure, alpha in FROZEN_ALPHA.items():
            source = complete_source(pressure, alpha, include_ion_temperature)
            nu_engine = engine_reference(pressure, include_ion_temperature)
            candidate = nu_engine / (1.0 + alpha)
            source_ep = complete_source(pressure, 0.0, include_ion_temperature)
            rows.append(dict(
                p_o2_torr=pressure, alpha=alpha, h=1.0 / (1.0 + alpha),
                ion_temperature='gas' if include_ion_temperature else 'neglected',
                nu_engine_ep=nu_engine, nu_engine_times_h=candidate,
                f_z=(math.pi / LENGTH)**2 / ((2.405/RADIUS)**2+(math.pi/LENGTH)**2),
                nu_radial_ep=nu_engine*((2.405/RADIUS)**2 / ((2.405/RADIUS)**2+(math.pi/LENGTH)**2)),
                nu_axial_ep=nu_engine*((math.pi/LENGTH)**2 / ((2.405/RADIUS)**2+(math.pi/LENGTH)**2)),
                closure_difference_relative_ep=(alpha/(1.+alpha)) * ((math.pi/LENGTH)**2 / ((2.405/RADIUS)**2+(math.pi/LENGTH)**2)),
                closure_difference_relative_full=alpha * ((math.pi/LENGTH)**2 / ((2.405/RADIUS)**2+(math.pi/LENGTH)**2)),
                source_ratio=source['nu_center_converted'] /
                source_ep['nu_center_converted'],
                ratio_to_candidate=source['nu_center_converted'] / candidate,
                **source))
    return {'rmgpy_file': str(Path(rmgpy.__file__).resolve()),
            'chi': CHI, 'volume_to_center': BETA, 'rows': rows}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    arguments = parser.parse_args()
    result = measurements()
    payload = json.dumps(result, indent=2, allow_nan=False) + '\n'
    if arguments.output:
        arguments.output.write_text(payload)
    print(payload, end='')
