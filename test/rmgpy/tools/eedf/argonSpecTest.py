"""Production-spec refresh preserves the grid and nitrogen inputs."""
import json
import os
from pathlib import Path

import pytest


def test_real_spec_refresh_preserves_grid_and_reports_all_argon_owners(tmp_path):
    from rmgpy.tools.eedf.argon import refresh_spec
    required = ('EEDF_REAL_SPEC', 'EEDF_DATABASE_ROOT', 'EEDF_ARGON_REALIZATION')
    if any(not os.environ.get(key) for key in required):
        pytest.skip('pinned local EEDF inputs are supplied for this integration check')
    base_path = Path(os.environ['EEDF_REAL_SPEC'])
    original_bytes = base_path.read_bytes()
    original = json.loads(original_bytes)
    destination = tmp_path / 'refreshed'
    report = refresh_spec(base_path, os.environ['EEDF_DATABASE_ROOT'],
                          os.environ['EEDF_ARGON_REALIZATION'], destination)
    assert report['outcome'] == 'converted'
    assert len(report['paired_owners']) == 10
    exceptions = report['missing_owners']
    assert len(exceptions) == 137
    assert all(row['ruled_exception'] for row in exceptions)
    assert sum(row['status'] == 'eliminated-upper-level' for row in exceptions) == 33
    assert len(report['branch_families']) == 137
    spec_path = destination / ('ar-n2-spec.json' if 'N2' in original['arm']['gases'] else 'pure-ar-spec.json')
    result = json.loads(spec_path.read_text())
    from rmgpy.tools.eedf.argon_population import AXIS, NODES, CONVENTION
    assert {k: v for k, v in result['axes'].items() if k != AXIS} == original['axes']
    assert result['axes'][AXIS] == NODES
    assert result['argon_4s_population'] == CONVENTION
    assert 'excited_fraction_max' not in result
    assert report['zero_target_node_check']['passed']
    assert result['interpolation_accuracy'] == {k: {'rtol': .01, 'atol': 0.} for k in ('mobility_N', 'field_power')}
    assert any(0. < p[AXIS] < 1.e-7 for p in result['held_out']['off_centre'])
    assert result['argon_elimination_electron_density_max_m3'] == 1.e18
    assert all(f'Ar({state}) = 0' in result['state_properties']['population']
               for state in ('3P2', '3P1', '3P0', '1P1'))
    for key, value in original['input_files'].items():
        if 'Ar_LXCat' not in key and 'Ar4s_assembled' not in key:
            assert result['input_files'][key] == value
    old_channels = json.loads((base_path.parent / original['channel_map']['path']).read_text())
    new_channels = json.loads(Path(result['channel_map']['path']).read_text())
    old_nitrogen = {c['description']: c for c in old_channels if c['description'].startswith('e + N2(')}
    new_nitrogen = {c['description']: c for c in new_channels if c['description'].startswith('e + N2(')}
    assert old_nitrogen == new_nitrogen
    assert base_path.read_bytes() == original_bytes
