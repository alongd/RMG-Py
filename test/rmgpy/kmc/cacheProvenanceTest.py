"""Mutation controls for expensive-generation cache provenance."""

import hashlib
import json

import pytest

from cache_provenance import generator_code_unchanged, supplied_artifact


def test_explicit_stale_artifact_is_rejected_without_cache_fallback(tmp_path, monkeypatch):
    payload = json.dumps({"provenance": {"compiler_sources_sha256": "stale"}}).encode()
    path = tmp_path / (hashlib.sha256(payload).hexdigest() + ".json")
    path.write_bytes(payload)
    monkeypatch.setenv("RMG_KMC_ARTIFACT", str(path))
    with pytest.raises(AssertionError):
        supplied_artifact()


def test_explicit_missing_artifact_is_rejected_without_cache_fallback(tmp_path, monkeypatch):
    monkeypatch.setenv("RMG_KMC_ARTIFACT", str(tmp_path / "missing.json"))
    with pytest.raises(FileNotFoundError):
        supplied_artifact()


def test_cache_reuse_accepts_only_kmc_changes(monkeypatch):
    generator_code_unchanged.cache_clear()
    monkeypatch.setattr(
        "cache_provenance.subprocess.check_output",
        lambda *args, **kwargs: "rmgpy/kmc/compiler.py\nrmgpy/kmc/state.py\n",
    )
    assert generator_code_unchanged("repository", "origin", "current")


def test_cache_reuse_rejects_generator_change(monkeypatch):
    generator_code_unchanged.cache_clear()
    monkeypatch.setattr(
        "cache_provenance.subprocess.check_output",
        lambda *args, **kwargs: "rmgpy/data/kinetics/family.py\n",
    )
    assert not generator_code_unchanged("repository", "origin", "current")


def test_cache_reuse_checks_uncommitted_generator_changes(monkeypatch):
    generator_code_unchanged.cache_clear()
    replies = iter(("", "rmgpy/molecule/resonance.py\n"))
    monkeypatch.setattr(
        "cache_provenance.subprocess.check_output",
        lambda *args, **kwargs: next(replies),
    )
    assert not generator_code_unchanged("repository", "origin", "current")
