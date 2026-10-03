#!/usr/bin/env python3

"""Resolved wall manifests retain identity without changing ground records."""

import ast
import importlib.util
import inspect
import json
from pathlib import Path

from rmgpy.chemkin import get_species_identifier
from rmgpy.rmg.listener import SimulationProfileWriter


ROOT = Path(__file__).resolve().parents[3]
_wall_spec = importlib.util.spec_from_file_location(
    'rework10_wall_fixture', ROOT / 'test/rmgpy/solver/plasmaElectronegativeWallTest.py')
_wall = importlib.util.module_from_spec(_wall_spec)
_wall_spec.loader.exec_module(_wall)


def resolved_wall_model():
    """Construct the review's resolved Cl- before reactor initialization."""
    namespace = dict(_wall.model.__globals__)
    tree = ast.parse(inspect.getsource(_wall.model))
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and node.value == '1 Cl u0 p4 c-1':
            node.value = 'electronicstate 1S0\n1 Cl u0 p4 c-1'
    exec(compile(ast.fix_missing_locations(tree), 'resolved_wall_model', 'exec'), namespace)
    return namespace['model']()


def test_resolved_wall_manifest_uses_state_qualified_species_names():
    reactor, species, _ = resolved_wall_model()
    reactor.monitor_electronegative_wall(reactor.y.copy(), 7.0)
    manifest = reactor.electronegative_wall_manifest()
    anion = next(spc for spc in species if spc.label == 'Cl-')
    identity = get_species_identifier(anion)

    assert list(manifest['gates']['A']) == [identity]
    assert manifest['extrema']['min_conf']['anion'] == identity
    assert all(identity in equation for equation, _ in
               manifest['gates']['A'][identity]['destruction_channels'])
    json.dumps(manifest, allow_nan=False)


def test_resolved_wall_profile_writes_state_in_json_and_csv(tmp_path):
    reactor, species, _ = resolved_wall_model()
    reactor.monitor_electronegative_wall(reactor.y.copy(), 7.0)
    manifest = reactor.electronegative_wall_manifest()
    reactor.snapshots = [[7.0, reactor.compute_volume(reactor.y)] + list(reactor.y[:len(species)])]
    reactor.electronegative_wall_history = {7.0: manifest}
    solver = tmp_path / 'solver'
    solver.mkdir()

    SimulationProfileWriter(str(tmp_path), 0, species).update(reactor)

    anion = next(spc for spc in species if spc.label == 'Cl-')
    identity = get_species_identifier(anion)
    profile = next(solver.glob('*.csv')).read_text()
    written = json.loads(next(solver.glob('*.json')).read_text())
    assert identity in profile.splitlines()[0]
    assert identity in written['gates']['A']
    assert written['extrema']['min_conf']['anion'] == identity


def test_ground_wall_manifest_keeps_historical_bytes(tmp_path):
    reactor, species, _ = _wall.model()
    reactor.monitor_electronegative_wall(reactor.y.copy(), 7.0)
    manifest = reactor.electronegative_wall_manifest()
    reactor.snapshots = [[7.0, reactor.compute_volume(reactor.y)] + list(reactor.y[:len(species)])]
    reactor.electronegative_wall_history = {7.0: manifest}
    solver = tmp_path / 'solver'
    solver.mkdir()
    expected = ROOT / 'test/rmgpy/test_data/excited_export_ground_wall_manifest.json'

    SimulationProfileWriter(str(tmp_path), 0, species).update(reactor)

    assert next(solver.glob('*.json')).read_bytes() == expected.read_bytes()
