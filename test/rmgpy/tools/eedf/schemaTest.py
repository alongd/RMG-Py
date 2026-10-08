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
import pytest

from rmgpy.tools.eedf.schema import _content_hash, _check_model, FingerprintMismatch


def test_fingerprint_ignores_repository_commits_but_checks_row_inputs():
    _row_inputs = {'cross_sections': {'ar': 'abc'}, 'Tg_K': 298.15, 'P_Pa': 666.61}
    manifest = {'row_inputs': _row_inputs, 'fingerprint': _content_hash(_row_inputs),
                'repositories': {'rmgpy': 'old', 'database': 'old'}}
    _check_model(manifest, dict(_row_inputs), relative_tolerance=1e-9)
    changed = dict(_row_inputs, cross_sections={'ar': 'def'})
    with pytest.raises(FingerprintMismatch, match='cross_sections'):
        _check_model(manifest, changed, relative_tolerance=1e-9)


def test_json_spec_preserves_exponential_numeric_tolerances(tmp_path):
    from rmgpy.tools.eedf.schema import load_spec, SpecError
    import json
    import os
    path = tmp_path / 'spec.json'
    spec = json.loads(Path(os.environ['EEDF_REAL_SPEC']).read_text())
    spec['tolerances']['H1']['rtol'] = 1e-6
    path.write_text(json.dumps(spec))
    loaded = load_spec(path)
    assert loaded['tolerances']['H1']['rtol'] == 1e-6
    spec['tolerances']['H1']['rtol'] = 'not numeric'
    path.write_text(json.dumps(spec))
    with pytest.raises(SpecError, match='numeric'):
        load_spec(path)


def test_spec_accepts_complete_runtime_qualification_policy(tmp_path):
    from rmgpy.tools.eedf.schema import load_spec, SpecError
    import json
    import os
    path = tmp_path / 'spec.json'
    spec = json.loads(Path(os.environ['EEDF_REAL_SPEC']).read_text())
    spec['engine_state_map'] = [{
        'formula': 'Ar', 'electronic_state': '', 'vibrational_level': None,
        'multiplicity': 1, 'loki_state': 'Ar(1S0)', 'statistical_weight': 1.,
    }]
    for key in ('A1', 'A2', 'A3', 'A4', 'A6'):
        spec['tolerances'][key] = {'rtol': 1.e-3, 'atol': 0.}
    spec['tolerances']['A5'] = {'rtol': 1.e-2, 'atol': 1.e-3,
                                'share_min': 1.e-2}
    spec['qualification_quantity_rules'] = {'EN_Td': 'A2'}
    path.write_text(json.dumps(spec))
    assert load_spec(path)['engine_state_map'][0]['loki_state'] == 'Ar(1S0)'

    del spec['tolerances']['A6']
    path.write_text(json.dumps(spec))
    with pytest.raises(SpecError, match='tolerances fields'):
        load_spec(path)


def test_generator_only_helpers_are_not_public_entries():
    import rmgpy.tools.eedf.artifact as artifact
    import rmgpy.tools.eedf.channels as channels
    import rmgpy.tools.eedf.integrity as integrity
    import rmgpy.tools.eedf.loki as loki
    import rmgpy.tools.eedf.moments as moments
    import rmgpy.tools.eedf.schema as schema
    import rmgpy.tools.eedf.validation as validation

    retired = {
        artifact: ('storage_key', 'storage_items', 'write_rows',
                   'validate_storage', 'write_artifact'),
        channels: ('validate_map', 'validate_physical_map'),
        integrity: ('solver_environment', 'check_solver',
                    'physical_properties', 'legacy_string',
                    'neutral_state_key', 'resolved_properties',
                    'channel_fractions'),
        loki: ('setup_text', 'enrich_row'),
        moments: ('distribution_moments',),
        schema: ('canonical_json', 'content_hash', 'file_hash', 'check_model',
                 'validate_spec', 'interpolant_identity', 'row_inputs'),
        validation: ('quantities', 'rows_agree', 'split_branches',
                     'held_out_points', 'freeze_accuracy_criteria',
                     'compare_row', 'channel_power_checks'),
    }
    for module, names in retired.items():
        assert all(not hasattr(module, name) for name in names)
