"""Mutation controls for expensive-generation cache provenance."""

import hashlib
import json
from pathlib import Path

import pytest

from cache_provenance import generator_code_unchanged, supplied_artifact


def test_explicit_stale_artifact_is_rejected_without_cache_fallback(tmp_path, monkeypatch):
    monkeypatch.delenv("RMG_KMC_ALLOW_STALE_ARTIFACT", raising=False)
    payload = json.dumps({"provenance": {"compiler_sources_sha256": "stale"}}).encode()
    path = tmp_path / (hashlib.sha256(payload).hexdigest() + ".json")
    path.write_bytes(payload)
    monkeypatch.setenv("RMG_KMC_ARTIFACT", str(path))
    with pytest.raises(AssertionError):
        supplied_artifact()


def test_stale_diagnostic_mode_still_validates_structure(tmp_path, monkeypatch):
    payload = json.dumps({"provenance": {"compiler_sources_sha256": "stale"}}).encode()
    path = tmp_path / (hashlib.sha256(payload).hexdigest() + ".json")
    path.write_bytes(payload)
    monkeypatch.setenv("RMG_KMC_ARTIFACT", str(path))
    monkeypatch.setenv("RMG_KMC_ALLOW_STALE_ARTIFACT", "1")
    with pytest.warns(UserWarning, match="phase-2b"):
        with pytest.raises(ValueError, match="schema"):
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


def test_oracle_cache_supports_the_explicit_materialized_database_pin(monkeypatch, tmp_path):
    import completeness_oracle as oracle
    monkeypatch.setenv("RMG_DATABASE_SHA", "materialized-pin")
    calls = []
    def repository_head(command, **kwargs):
        calls.append(command)
        return "repository-head\n"
    monkeypatch.setattr(oracle.subprocess, "check_output", repository_head)
    repository = Path(__file__).resolve().parents[3]
    database = tmp_path / "database"
    (database / "input/kinetics").mkdir(parents=True)
    (database / "input/thermo").mkdir(parents=True)
    (database / "input/kinetics/rules.py").write_text("rules = 1\n")
    (database / "input/thermo/library.py").write_text("library = 1\n")
    key = oracle.oracle_cache_key(repository, database, 1)
    assert key.startswith("repository-head-materialized-pin-")
    assert len(calls) == 1 and calls[0][2] == str(repository)
