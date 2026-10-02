"""Mutation controls for expensive-generation cache provenance."""

from cache_provenance import generator_code_unchanged


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
