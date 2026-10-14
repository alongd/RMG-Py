"""The frozen I-052 baseline tests skip when their inputs are absent and stay
strict when the inputs are present (fresh clones at another prefix must not
fail them).

Every test builds a complete synthetic set of inputs in a temp directory and
removes exactly one, so each exercises only its own condition on any machine.
"""

import json

import pytest

import i052_baseline_freeze as freeze
import i052BaselineFreezeTest as frozen_tests


@pytest.fixture
def inputs(monkeypatch, tmp_path):
    """Point every probed input at a present synthetic stand-in."""
    artifact = tmp_path / "artifact.json"
    design = tmp_path / "design.md"
    database = tmp_path / "database"
    artifact.write_text("{}")
    design.write_text("design")
    database.mkdir()
    environment = {"python_version": "synthetic"}
    manifest = {
        "artifact": {"path": str(artifact), "sha256": "0" * 64},
        "design": {"path": str(design), "sha256": "0" * 64},
        "environment": environment,
    }
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest))
    monkeypatch.setattr(freeze, "ARTIFACT", artifact)
    monkeypatch.setattr(freeze, "DATABASE", database)
    monkeypatch.setattr(freeze, "_has_commit", lambda root, commit: True)
    monkeypatch.setattr(freeze, "environment_pins", lambda: environment)
    monkeypatch.setattr(frozen_tests, "MANIFEST", manifest_path)
    return manifest, manifest_path, tmp_path


def test_complete_inputs_have_no_skip_reason(inputs):
    manifest, _path, _tmp = inputs
    assert freeze.baseline_skip_reason(manifest) is None


def test_missing_artifact_skips_with_reason(inputs, monkeypatch):
    _manifest, _path, tmp = inputs
    monkeypatch.setattr(freeze, "ARTIFACT", tmp / "absent.json")
    with pytest.raises(pytest.skip.Exception, match="artifact"):
        frozen_tests.test_manifest_hashes_are_immutable()
    with pytest.raises(pytest.skip.Exception, match="artifact"):
        frozen_tests.test_replay_frozen_observations()


def test_missing_database_skips_with_reason(inputs, monkeypatch):
    _manifest, _path, tmp = inputs
    monkeypatch.setattr(freeze, "DATABASE", tmp / "no_database")
    with pytest.raises(pytest.skip.Exception, match="database"):
        frozen_tests.test_manifest_hashes_are_immutable()


def test_missing_design_reports_reason(inputs):
    manifest, _path, tmp = inputs
    manifest["design"]["path"] = str(tmp / "no_design.md")
    assert "design" in freeze.baseline_skip_reason(manifest)


def test_absent_historical_commit_skips_with_reason(inputs, monkeypatch):
    monkeypatch.setattr(freeze, "_has_commit", lambda root, commit: False)
    with pytest.raises(pytest.skip.Exception, match="commit"):
        frozen_tests.test_manifest_hashes_are_immutable()


def test_environment_difference_skips_with_reason(inputs, monkeypatch):
    monkeypatch.setattr(
        freeze, "environment_pins", lambda: {"python_version": "other"})
    with pytest.raises(pytest.skip.Exception, match="environment"):
        frozen_tests.test_manifest_hashes_are_immutable()


def test_present_inputs_stay_strict(inputs, monkeypatch):
    manifest, _path, _tmp = inputs
    monkeypatch.setattr(freeze, "_require_git_commit", lambda root, commit: None)
    monkeypatch.setattr(freeze, "_verify_frozen_json_bytes", lambda: None)
    manifest["pins"] = {}  # disagrees with the harness pins
    assert freeze.baseline_skip_reason(manifest) is None
    with pytest.raises(AssertionError):
        freeze.verify_manifest(manifest)
