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

from pathlib import Path
import numpy as np
import pytest
from rmgpy.tools.eedf.loki import parse_output


def test_stored_ar_output_has_independent_known_rates_and_transport():
    row = parse_output(Path(__file__).parent / 'fixtures' / 'ar', normalization_atol=1e-8)
    assert row['swarm']['mean_energy_eV'] == pytest.approx(5.63869451471794)
    assert row['swarm']['mobility_N'] == pytest.approx(1.01724123397156e24)
    assert row['k_ine'][-1] == pytest.approx(7.44149321789454e-21, rel=1e-12, abs=0)
    assert row['power_groups']['field'] == pytest.approx(3.21679922919037e-16, rel=1e-12, abs=0)
    assert len(row['channels']) == 39
    assert np.sum(row['f0'] * np.sqrt(row['energy_eV']) * np.diff(row['energy_edges_eV'])) == pytest.approx(1)


def test_attachment_power_uses_same_eedf_weighted_energy(tmp_path):
    from rmgpy.tools.eedf.loki import enrich_row
    from scipy.constants import elementary_charge, electron_mass
    gamma = np.sqrt(2 * elementary_charge / electron_mass)
    rate = gamma * 1e-20 * (.5 * .5 + 1.5)
    row = {'channels': ['e + X(gnd) -> X(-,gnd), Attachment'], 'composition': {}, 'energy_eV': np.array([.5, 1.5]),
           'energy_edges_eV': np.array([0., 1., 2.]), 'f0': np.ones(2),
           'k_ine': np.array([rate]), 'k_sup': np.array([0.])}
    channel = {'description': 'e + X(gnd) -> X(-,gnd), Attachment', 'kind': 'attachment', 'threshold_eV': 0.,
               'target_fraction': .25, 'product_fraction': 0., 'sigma_max_m2': 1e-20,
               'cross_section': {'energy_eV': [0., 2.], 'sigma_m2': [1e-20, 1e-20]}}
    spec = {'Tg_K': 298.15, 'envelopes': {}, 'input_files': {}, 'solver_options': {}, 'gas_properties': {'fraction': ['X = 1.0']},
            'state_properties': {'population': ['X(gnd) = 0.25']},
            'floors': {'rate_absolute': 1e-40, 'eedf_dynamic_range': 1e-30}}
    from rmgpy.tools.eedf.schema import file_hash
    section = tmp_path / 'attachment.txt'
    section.write_text('PARAM.: E = 0 eV\nCOMMENT: [X(gnd) + e -> X(-,gnd), Attachment]\n--\n0 1e-20\n2 1e-20\n--\n')
    spec['input_files']['attachment.txt'] = {'path': str(section), 'sha256': file_hash(section), 'kind': 'cross_section'}
    channel.update(classification='B', flux_group=channel['description'], reaction=None)
    enriched = enrich_row(row, [channel], spec)
    # The zero-energy face is zero for attachment. Nodal averaging gives
    # cell sections [0.5, 1]*sigma, so the energy moment is 19/14 eV.
    assert enriched['attachment_energy_eV'][0] == pytest.approx(19 / 14)
    assert enriched['channel_power'][0] == pytest.approx(.25 * rate * 19 / 14, rel=1e-12, abs=0)
