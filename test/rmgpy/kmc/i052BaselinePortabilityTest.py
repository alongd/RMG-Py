"""The frozen I-052 baseline tests skip when their inputs are absent and stay
strict when the inputs are present (fresh clones at another prefix must not
fail them)."""

import json

import pytest

import i052_baseline_freeze as freeze
import i052BaselineFreezeTest as frozen_tests


def test_missing_artifact_skips_with_reason(monkeypatch, tmp_path):
    monkeypatch.setattr(freeze, "ARTIFACT", tmp_path / "absent.json")
    with pytest.raises(pytest.skip.Exception, match="artifact"):
        frozen_tests.test_manifest_hashes_are_immutable()
    with pytest.raises(pytest.skip.Exception, match="artifact"):
        frozen_tests.test_replay_frozen_observations()


def test_missing_database_or_design_skips_with_reason(monkeypatch, tmp_path):
    monkeypatch.setattr(freeze, "DATABASE", tmp_path / "no_database")
    with pytest.raises(pytest.skip.Exception, match="database"):
        frozen_tests.test_manifest_hashes_are_immutable()
    monkeypatch.undo()
    manifest = json.loads(freeze.MANIFEST.read_text())
    manifest["design"]["path"] = str(tmp_path / "no_design.md")
    assert "design" in freeze.baseline_skip_reason(manifest)


def test_absent_historical_commit_skips_with_reason(monkeypatch):
    monkeypatch.setattr(freeze, "RMG_PY_PIN", "0" * 40)
    with pytest.raises(pytest.skip.Exception, match="commit"):
        frozen_tests.test_manifest_hashes_are_immutable()


def test_environment_difference_skips_with_reason(monkeypatch):
    monkeypatch.setattr(
        freeze, "environment_pins", lambda: {"python_version": "0.0.0"})
    with pytest.raises(pytest.skip.Exception, match="environment"):
        frozen_tests.test_manifest_hashes_are_immutable()


def test_present_inputs_stay_strict():
    manifest = json.loads(freeze.MANIFEST.read_text())
    if freeze.baseline_skip_reason(manifest) is not None:
        pytest.skip("frozen baseline inputs absent on this machine")
    corrupted = json.loads(json.dumps(manifest))
    corrupted["artifact"]["sha256"] = "0" * 64
    with pytest.raises(AssertionError):
        freeze.verify_manifest(corrupted)
