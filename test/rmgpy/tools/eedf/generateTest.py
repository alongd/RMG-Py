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

import copy
import json
import os
from pathlib import Path
import numpy as np
import pytest

from rmgpy.tools.eedf.generate import generate, model_inputs
from rmgpy.tools.eedf.schema import file_hash, SpecError
from rmgpy.tools.eedf.validation import held_out_points
from rmgpy.solver.eedf import EEDFTable


def test_held_out_midpoints_and_seeded_lhs_never_enter_training_axes():
    axes = {'u': [0., 1., 2.], 'x': [0., 1.]}
    policy = {'cell_midpoints': True, 'lhs_count': 32, 'seed': 314159}
    points = held_out_points(axes, policy)
    assert {'u': .5, 'x': .5} in points
    assert {'u': .5, 'x': 0.} in points
    assert {'u': 0., 'x': .5} in points
    assert points == held_out_points(axes, policy)
    assert all(not (p['u'] in axes['u'] and p['x'] in axes['x']) for p in points)


def test_real_pure_ar_table_passes_held_out_and_three_scans(tmp_path):
    spec_path = os.environ.get('EEDF_REAL_SPEC')
    if not spec_path:
        pytest.skip('EEDF_REAL_SPEC supplies the pinned real LoKI integration fixture')
    from rmgpy.tools.eedf.schema import load_spec
    spec = load_spec(spec_path)
    spec['axes']['u'] = np.linspace(np.log(16), np.log(18), 9).tolist()
    spec['held_out']['lhs_count'] = 4
    spec['scratch_root'] = str(tmp_path / 'loki')
    spec['output_root'] = str(tmp_path / 'tables')
    spec['refinement']['max_rounds'] = 0
    path = generate(spec, command_line=['pytest real LoKI fixture'])
    table = EEDFTable.load(path, model_inputs(spec), artifact_sha256=file_hash(path / 'table.h5'))
    assert table.manifest['accepted']
    assert all(scan['agreement'] for scan in table.manifest['branch_detection'])
    assert len(table.manifest['held_out_verdicts']) >= 12
    assert table.row(np.log(17), {}).swarm['mean_energy_eV'] > 5
    assert table.manifest['row_inputs']['Tg_K'] == 298.15
    assert table.manifest['row_inputs']['P_Pa'] == pytest.approx(5 * 133.322368)
    expected = {'schema_version', 'row_inputs', 'fingerprint', 'axes', 'energy_eV', 'energy_edges_eV',
                'envelopes', 'screen_results', 'tolerances', 'floors', 'channel_map', 'repositories',
                'solver', 'setups', 'state_population_assumptions', 'density_convention', 'units',
                'generator_command', 'generated_at', 'generation_seconds', 'held_out',
                'held_out_verdicts', 'accepted', 'branch_detection', 'limitations',
                'refinement_history', 'branches', 'artifact_sha256', 'branch_certification',
                'qualification_setup_sha256'}
    assert set(table.manifest) == expected
    assert 'cold' in table.manifest['branch_detection'][0]['branch_sources']['branch_0']


@pytest.mark.parametrize('field', ['floors', 'tolerances', 'screen', 'refinement', 'held_out'])
def test_missing_generation_policy_refuses(field):
    from rmgpy.tools.eedf.schema import load_spec, validate_spec
    spec_path = os.environ.get('EEDF_REAL_SPEC')
    if not spec_path:
        pytest.skip('EEDF_REAL_SPEC supplies the complete fixture')
    spec = load_spec(spec_path)
    del spec[field]
    with pytest.raises(SpecError, match='spec fields'):
        validate_spec(spec)
