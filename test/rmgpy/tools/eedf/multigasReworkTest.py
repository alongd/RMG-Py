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

"""Regressions for map order, Effective realization and diagnostic policies."""
import copy
import json
import os
from pathlib import Path

import numpy as np
import pytest

from rmgpy.tools.eedf.channel_map import _diagnostic_d6
from rmgpy.tools.eedf.channels import _effective_elastic_section, _validate_physical_map
from rmgpy.tools.eedf.generate import _solve
from rmgpy.tools.eedf.integrity import _channel_fractions, _resolved_properties
from rmgpy.tools.eedf.loki import parse_output
from rmgpy.tools.eedf.n2_vibrational import LadderError, _policy
from rmgpy.tools.eedf.schema import SpecError
from rmgpy.tools.eedf.validation import _compare_row
from rmgpy.tools.eedf.vibrational import qualify_ladder


@pytest.mark.parametrize('reverse', [False, True])
def test_populated_electronic_excitation(reverse):
    arrow, product = ('<->', 'Ar(3P2)') if reverse else ('->', 'Ar(4p,group)')
    channel = {'description': f'e + Ar(1S0) {arrow} e + {product}, Excitation',
               'kind': 'excitation', 'target_fraction': .8,
               'product_fraction': .2 if reverse else 0.}
    populations = {'Ar(1S0)': .8, 'Ar(3P2)': .2}
    assert _channel_fractions([channel], {'Ar': 1.}, populations, {}) == ([.8], [.2 if reverse else 0.])
    channel['product_fraction'] = 0. if reverse else .2
    with pytest.raises(SpecError, match='product_fraction differs'):
        _channel_fractions([channel], {'Ar': 1.}, populations, {})


@pytest.mark.parametrize('spoof_owner', [False, True])
def test_diagnostic_d6_cannot_qualify(spoof_owner):
    policy = _diagnostic_d6()
    assert policy['set_by'] == 'diagnostic'
    if spoof_owner:
        policy['set_by'] = 'owner'
    _policy(policy, allow_diagnostic=True)
    with pytest.raises(LadderError, match='diagnostic D6'):
        _policy(policy)
    with pytest.raises(LadderError, match='diagnostic D6'):
        qualify_ladder({'K': 10, 'axis': 'Tv_K', 'D6': policy},
                       [1.] + [0.] * 10, {'forward_power': 1., 'top_leak_power': 0.})


def test_diagnostic_d6_refuses_supplied_owner_limits():
    owner = dict(_diagnostic_d6(), set_by='owner', status='production')
    before = copy.deepcopy(owner)
    with pytest.raises(SpecError, match='diagnostic D6 conflicts'):
        _diagnostic_d6(owner)
    assert owner == before
    _policy(owner)
    assert qualify_ladder({'K': 10, 'axis': 'Tv_K', 'D6': owner},
                          [1.] + [0.] * 10, {'forward_power': 1., 'top_leak_power': 0.})['action'] == 'accept'


@pytest.mark.parametrize('permutation', [False, True])
def test_saved_real_effective_row_through_solve_and_compare(permutation):
    """Use pinned external inputs/output; never vendor collision payloads.

    EEDF_MULTIGAS_SPEC is the converter's spec and EEDF_MULTIGAS_OUTPUT is
    the saved reducedField_60 directory at x_N2=.1, Tv=3000 K.
    """
    spec_path, output = os.environ.get('EEDF_MULTIGAS_SPEC'), os.environ.get('EEDF_MULTIGAS_OUTPUT')
    if not spec_path or not output:
        pytest.skip('requires saved real Effective+inelastic mixture output')
    spec = json.loads(Path(spec_path).read_text())
    channels = json.loads(Path(spec['channel_map']['path']).read_text())
    coordinate = {'x_N2': .1, 'Tv_K': 3000.,
                  'electron_density': spec['envelopes']['electron_density']['reference']}
    _, states, _, _ = _resolved_properties(spec, coordinate)
    # Saved pre-fix maps retain the Effective section on their derived row.
    for channel in channels:
        if channel['kind'] == 'effective':
            derived = next(c for c in channels if c['description'] ==
                           channel['description'].replace(', Effective', ', Elastic'))
            section = _effective_elastic_section(channel, channels, states)
            derived.update(cross_section=section, sigma_max_m2=max(section['sigma_m2']))
    grid = spec['solver_options']['numerics']['energyGrid']
    row = parse_output(Path(output), spec['tolerances']['normalization_atol'],
                       np.linspace(0, grid['maxEnergy'], grid['cellNumber'] + 1))
    row.update(composition=coordinate, branch_id='only', converged=True)
    original_rates = dict(zip(row['channels'], row['k_ine']))
    assert len(set(value for value in original_rates.values() if value > 0)) > 10
    if permutation:
        channels = channels[::-1]
    original_map = copy.deepcopy(channels)

    class SavedDriver:
        def run(self, coordinates, fields, label):
            return [copy.deepcopy(row)]

    result = _solve(SavedDriver(), channels, spec, coordinate, [np.log(60.)], 'saved')[0]
    assert channels == original_map
    assert result['channels'] == [c['description'] for c in channels]
    for description, rate in zip(result['channels'], result['k_ine']):
        assert rate == original_rates[description]
    manifest = {'channel_map': channels, 'tolerances': spec['tolerances'],
                'floors': spec['floors'], 'Tg_K': spec['Tg_K']}
    failed = [check for check in _compare_row(result, result, manifest) if not check['passed']]
    assert not failed, failed
    source = next(c for c in channels if c['kind'] == 'effective')
    derived = next(c for c in channels if c['description'] ==
                   source['description'].replace(', Effective', ', Elastic'))
    assert derived['cross_section'] != source['cross_section']
    assert _effective_elastic_section(source, channels[::-1], states) == derived['cross_section']
    # A forged source copy must not pass provenance validation on the derived row.
    derived.update(cross_section=copy.deepcopy(source['cross_section']),
                   sigma_max_m2=source['sigma_max_m2'])
    with pytest.raises(SpecError, match='derived Elastic differs'):
        _validate_physical_map(channels, spec, coordinate)
