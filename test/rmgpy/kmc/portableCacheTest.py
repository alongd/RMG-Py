import json
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile

import pytest

from portable_cache import (
    artifact_cache_key,
    compile_environment_options,
    database_identity,
    export_cache,
    identity,
    import_cache,
    load_or_generate,
    migrate,
    validate_artifact_cache_entry,
)


def _rewrite_archive(archive, mutate):
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary) / "kmc-cache"
        with tarfile.open(archive, "r:gz") as tar:
            tar.extractall(temporary, filter="data")
        root = Path(temporary) / "kmc-cache"
        mutate(root)
        replacement = Path(temporary) / "replacement.tgz"
        with tarfile.open(replacement, "w:gz") as tar:
            tar.add(root, arcname="kmc-cache")
        shutil.copy2(replacement, archive)


def _make_database(path):
    (path / "input/kinetics").mkdir(parents=True)
    (path / "input/thermo").mkdir(parents=True)
    (path / "input/kinetics/marker.txt").write_text("kinetics")
    (path / "input/thermo/marker.txt").write_text("thermo")
    return path


def test_identity_excludes_kmc_and_changes_for_other_rmgpy(tmp_path, monkeypatch):
    repo = tmp_path / "repo"
    (repo / "rmgpy/kmc").mkdir(parents=True)
    (repo / "rmgpy/data").mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / "rmgpy/kmc/x.py").write_text("one")
    (repo / "rmgpy/data/x.py").write_text("one")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    db = _make_database(tmp_path / "db")
    before = identity(repo, db)
    (repo / "rmgpy/kmc/x.py").write_text("two")
    assert identity(repo, db) == before
    (repo / "rmgpy/data/x.py").write_text("two")
    assert identity(repo, db) != before


def test_identity_hashes_database_once_per_fixture_call(tmp_path, monkeypatch):
    repo = tmp_path / "repo"
    (repo / "rmgpy").mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / "rmgpy/x.py").write_text("x")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    db = _make_database(tmp_path / "db")
    calls = []
    monkeypatch.setattr("portable_cache.database_content_digest", lambda path: calls.append(path) or ("digest", ()))

    identity(repo, db)

    assert calls == [db]


def test_export_import_and_manifest_mismatch(tmp_path):
    repo = tmp_path / "repo"
    (repo / "rmgpy").mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / "rmgpy/x.py").write_text("x")
    for relative in ("test/rmgpy/kmc/compile_event_set_fixture.py", "test/rmgpy/kmc/portable_cache.py"):
        path = repo / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("driver")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    db = _make_database(tmp_path / "db")
    source = tmp_path / "cache"
    (source / "generated-reactions").mkdir(parents=True)
    (source / "generated-reactions/entry.pickle").write_bytes(b"entry")
    (source / "manifest.json").write_text(json.dumps(identity(repo, db)))
    archive = tmp_path / "cache.tgz"
    export_cache(source, archive, repo, db)
    destination = tmp_path / "destination"
    import_cache(archive, destination, repo, db)
    assert (destination / "generated-reactions/entry.pickle").read_bytes() == b"entry"
    (repo / "rmgpy/x.py").write_text("changed")
    with pytest.raises(ValueError, match="manifest"):
        import_cache(archive, tmp_path / "refused", repo, db)


def test_import_rejects_tampered_file(tmp_path):
    repo = tmp_path / "repo"
    (repo / "rmgpy").mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / "rmgpy/x.py").write_text("x")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    db = _make_database(tmp_path / "db")
    source = tmp_path / "cache/generated-reactions"
    source.mkdir(parents=True)
    (source / "entry.pickle").write_bytes(b"entry")
    (source.parent / "manifest.json").write_text(json.dumps(identity(repo, db)))
    archive = tmp_path / "cache.tgz"
    export_cache(source.parent, archive, repo, db)
    _rewrite_archive(archive, lambda root: (root / "generated-reactions/entry.pickle").write_bytes(b"tampered"))
    with pytest.raises(ValueError, match="hash"):
        import_cache(archive, tmp_path / "refused", repo, db)


def test_import_rejects_stray_file(tmp_path):
    repo = tmp_path / "repo"
    (repo / "rmgpy").mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / "rmgpy/x.py").write_text("x")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    db = _make_database(tmp_path / "db")
    source = tmp_path / "cache/generated-reactions"
    source.mkdir(parents=True)
    (source / "entry.pickle").write_bytes(b"entry")
    (source.parent / "manifest.json").write_text(json.dumps(identity(repo, db)))
    archive = tmp_path / "cache.tgz"
    export_cache(source.parent, archive, repo, db)
    _rewrite_archive(archive, lambda root: (root / "stray.txt").write_bytes(b"stray"))
    with pytest.raises(ValueError, match="unknown"):
        import_cache(archive, tmp_path / "refused", repo, db)


def test_import_rejects_symlink_before_copying(tmp_path):
    repo = tmp_path / "repo"
    (repo / "rmgpy").mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / "rmgpy/x.py").write_text("x")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    db = _make_database(tmp_path / "db")
    source = tmp_path / "cache/generated-reactions"
    source.mkdir(parents=True)
    (source / "entry.pickle").write_bytes(b"entry")
    (source.parent / "manifest.json").write_text(json.dumps(identity(repo, db)))
    archive = tmp_path / "cache.tgz"
    export_cache(source.parent, archive, repo, db)
    _rewrite_archive(archive, lambda root: (root / "generated-reactions/link.pickle").symlink_to("entry.pickle"))
    destination = tmp_path / "refused"
    with pytest.raises(ValueError, match="symlink"):
        import_cache(archive, destination, repo, db)
    assert not destination.exists()


def test_migration_copies_old_named_entries(tmp_path, monkeypatch):
    monkeypatch.setenv("PYTHONHASHSEED", "0")
    repo = tmp_path / "repo"
    (repo / "rmgpy").mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / "rmgpy/x.py").write_text("x")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    db = _make_database(tmp_path / "db")
    origin = subprocess.check_output(
        ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True
    ).strip()
    seed = os.environ["PYTHONHASHSEED"]
    old = tmp_path / f"generated-reactions/{origin}-declared-{database_identity(db)}/{seed}"
    old.mkdir(parents=True)
    (old / "x.pickle").write_bytes(b"x")
    target = migrate(tmp_path, repo, db, "declared")
    assert (target / "generated-reactions/x.pickle").read_bytes() == b"x"
    assert (old / "x.pickle").exists()


def test_migration_rejects_stale_database_digest(tmp_path, monkeypatch):
    monkeypatch.setenv("PYTHONHASHSEED", "0")
    repo = tmp_path / "repo"
    (repo / "rmgpy").mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / "rmgpy/x.py").write_text("x")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    db = _make_database(tmp_path / "db")
    old_digest = database_identity(db)
    origin = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    old = tmp_path / f"generated-reactions/{origin}-declared-{old_digest}/0"
    old.mkdir(parents=True)
    (old / "x.pickle").write_bytes(b"x")
    (db / "input/kinetics/marker.txt").write_text("changed")
    target = migrate(tmp_path, repo, db, "declared")
    assert not (target / "generated-reactions/x.pickle").exists()


def test_export_rejects_changed_database_against_bucket_manifest(tmp_path):
    repo = tmp_path / "repo"
    (repo / "rmgpy").mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / "rmgpy/x.py").write_text("x")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    db = _make_database(tmp_path / "db")
    source = tmp_path / "cache"
    (source / "generated-reactions").mkdir(parents=True)
    (source / "manifest.json").write_text(json.dumps(identity(repo, db)))
    (db / "input/kinetics/marker.txt").write_text("changed")
    with pytest.raises(ValueError, match="identity"):
        export_cache(source, tmp_path / "cache.tgz", repo, db)


def test_migration_skips_mismatched_database(tmp_path):
    repo = tmp_path / "repo"
    (repo / "rmgpy").mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / "rmgpy/x.py").write_text("x")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    db = _make_database(tmp_path / "db")
    origin = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    wrong = tmp_path / f"generated-reactions/{origin}-wrong-db/0"
    wrong.mkdir(parents=True)
    (wrong / "x.pickle").write_bytes(b"x")
    target = migrate(tmp_path, repo, db)
    assert not (target / "generated-reactions/x.pickle").exists()


def test_git_and_archive_copies_have_same_database_identity(tmp_path):
    git_db = tmp_path / "git-db"
    (git_db / "input/kinetics").mkdir(parents=True)
    (git_db / "input/thermo").mkdir(parents=True)
    (git_db / "input/kinetics/data.txt").write_text("same")
    (git_db / "input/thermo/data.txt").write_text("same")
    subprocess.run(["git", "init", "-q", str(git_db)], check=True)
    subprocess.run(["git", "-C", str(git_db), "add", "."], check=True)
    subprocess.run(["git", "-C", str(git_db), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    archive_db = tmp_path / "archive-db"
    shutil.copytree(git_db / "input", archive_db / "input")
    assert database_identity(git_db) == database_identity(archive_db)


def test_database_content_change_misses_portable_cache_with_same_declaration(tmp_path, monkeypatch):
    repo = tmp_path / "repo"
    (repo / "rmgpy").mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / "rmgpy/x.py").write_text("x")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    db = _make_database(tmp_path / "db")
    monkeypatch.setenv("RMG_DATABASE_SHA", "declared-revision")
    before = identity(repo, db)
    (db / "input/kinetics/marker.txt").write_text("changed")
    assert identity(repo, db) != before


def test_generation_cache_hit_does_not_call_stubbed_generator(tmp_path):
    path = tmp_path / "entry.pickle"
    calls = []
    generate = lambda: calls.append("called") or {"entry": 1}
    dump = lambda value: json.dumps(value).encode()
    load = lambda payload: json.loads(payload)
    assert load_or_generate(path, generate, dump, load) == {"entry": 1}
    assert load_or_generate(path, generate, dump, load) == {"entry": 1}
    assert calls == ["called"]


@pytest.mark.parametrize("seed", [None, "random"])
def test_uncacheable_generation_never_reads_persistent_entries(tmp_path, monkeypatch, seed):
    if seed is None:
        monkeypatch.delenv("PYTHONHASHSEED", raising=False)
    else:
        monkeypatch.setenv("PYTHONHASHSEED", seed)
    path = tmp_path / "uncacheable" / "same-pid" / "entry.pickle"
    path.parent.mkdir(parents=True)
    path.write_bytes(b"old")
    calls = []
    value = load_or_generate(
        path,
        lambda: calls.append("called") or {"new": True},
        lambda value: json.dumps(value).encode(),
        json.loads,
        read_cache=False,
    )
    assert value == {"new": True}
    assert calls == ["called"]


def test_artifact_key_includes_compile_options(tmp_path):
    repo = tmp_path / "repo"
    (repo / "rmgpy").mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / "rmgpy/x.py").write_text("x")
    for relative in ("test/rmgpy/kmc/compile_event_set_fixture.py", "test/rmgpy/kmc/portable_cache.py"):
        path = repo / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("driver")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    db = _make_database(tmp_path / "db")
    base = {"temperature_grid": [300, 325], "use_plpsec_library": "1"}
    disabled = {**base, "use_plpsec_library": "0"}
    assert artifact_cache_key(repo, db, base) != artifact_cache_key(repo, db, disabled)


def test_artifact_key_changes_when_compile_driver_changes(tmp_path):
    repo = tmp_path / "repo"
    (repo / "rmgpy").mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / "rmgpy/x.py").write_text("x")
    driver = repo / "test/rmgpy/kmc/compile_event_set_fixture.py"
    helper = repo / "test/rmgpy/kmc/portable_cache.py"
    driver.parent.mkdir(parents=True, exist_ok=True)
    driver.write_text("driver")
    helper.write_text("helper")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    db = _make_database(tmp_path / "db")
    before = artifact_cache_key(repo, db, {})
    driver.write_text("changed driver")
    assert artifact_cache_key(repo, db, {}) != before


def test_artifact_hit_rejects_misplaced_and_corrupt_entries(tmp_path):
    repo = tmp_path / "repo"
    (repo / "rmgpy").mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / "rmgpy/x.py").write_text("x")
    for relative in ("test/rmgpy/kmc/compile_event_set_fixture.py", "test/rmgpy/kmc/portable_cache.py"):
        path = repo / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("driver")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    db = _make_database(tmp_path / "db")
    bucket = tmp_path / "artifact"
    bucket.mkdir()
    artifact = bucket / "artifact.json"
    artifact.write_text("{}")
    manifest = bucket / "manifest.json"
    manifest.write_text(json.dumps({"identity": "wrong", "artifact": artifact.name, "artifact_sha256": "wrong"}))
    with pytest.raises(ValueError, match="identity"):
        validate_artifact_cache_entry(manifest, "expected", db, None)
    manifest.write_text(json.dumps({"identity": "expected", "artifact": artifact.name, "artifact_sha256": __import__("hashlib").sha256(artifact.read_bytes()).hexdigest()}))
    with pytest.raises(ValueError):
        validate_artifact_cache_entry(manifest, "expected", db, None)


def test_artifact_hit_checks_hash_before_json_validation(tmp_path, monkeypatch):
    repo = tmp_path / "repo"
    (repo / "rmgpy").mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / "rmgpy/x.py").write_text("x")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    db = _make_database(tmp_path / "db")
    bucket = tmp_path / "artifact"
    bucket.mkdir()
    artifact = bucket / "artifact.json"
    original = b'{"provenance": {"compiler_sources_sha256": "sources"}}\n'
    artifact.write_bytes(original)
    manifest = bucket / "manifest.json"
    manifest.write_text(json.dumps({
        "identity": "expected",
        "artifact": artifact.name,
        "artifact_sha256": hashlib.sha256(original).hexdigest(),
    }))
    artifact.write_bytes(b'{ "provenance": {"compiler_sources_sha256": "sources"} }\n')
    monkeypatch.setattr("rmgpy.kmc.compiler.compiler_source_hash", lambda: "sources")
    monkeypatch.setattr("rmgpy.kmc.compiler.validate_artifact", lambda value: None)
    monkeypatch.setattr("rmgpy.kmc.database_provenance.provenance_matches_database", lambda *args: True)
    with pytest.raises(ValueError, match="artifact hash mismatch"):
        validate_artifact_cache_entry(manifest, "expected", db, None)


def test_driver_main_cold_hit_and_generation_only_modes(tmp_path, monkeypatch, capsys):
    import compile_event_set_fixture as driver

    database_path = _make_database(tmp_path / "db")
    (database_path / "input/kinetics/families/stub").mkdir(parents=True)
    (database_path / "input/kinetics/families/stub/groups.py").write_text("")
    cache_root = tmp_path / "cache"
    calls = []

    class FakeKinetics:
        def generate_reactions_from_families(self, *args):
            calls.append("generate")
            return []

    class FakeDatabase:
        def __init__(self):
            self.kinetics = FakeKinetics()
            self.thermo = object()

        def load_kinetics(self, *args, **kwargs):
            pass

        def load_thermo(self, *args, **kwargs):
            pass

    class FakeCompiler:
        @staticmethod
        def discover_family_reactions(kinetics, proxies, families, family_universe=None):
            return {}, set(), {"stub": kinetics.generate_reactions_from_families([], [])}

        def __init__(self, *args, **kwargs):
            pass

        def write_artifact(self, directory):
            path = Path(directory) / "stub-artifact.json"
            path.write_text(json.dumps({"provenance": {"compiler_sources_sha256": "sources"}}))
            return path, None

    monkeypatch.setenv("PYTHONHASHSEED", "0")
    monkeypatch.setenv("RMG_KMC_CACHE_ROOT", str(cache_root))
    monkeypatch.setattr(driver, "RMGDatabase", FakeDatabase)
    monkeypatch.setattr(driver, "EventSetCompiler", FakeCompiler)
    monkeypatch.setattr(driver, "prepare_rate_rules", lambda *args, **kwargs: None)
    monkeypatch.setattr(driver, "ps_proxy_set", lambda units: ())
    monkeypatch.setattr(driver, "resolve_database_declaration", lambda value: value)
    monkeypatch.setattr(driver, "database_content_digest", lambda path: ("database", {}))
    monkeypatch.setattr(driver, "identity", lambda repo, db: {"identity": "generation"})
    monkeypatch.setattr(driver, "identity_name", lambda value: "generation")
    monkeypatch.setattr(driver, "artifact_cache_key", lambda repo, db, options: "artifact")
    monkeypatch.setattr(driver, "migrate", lambda *args: cache_root / "portable/generation")
    monkeypatch.setattr(
        driver,
        "validate_artifact_cache_entry",
        lambda manifest, key, db, declaration: manifest.parent / json.loads(manifest.read_text())["artifact"],
    )
    (tmp_path / "repo").mkdir()
    monkeypatch.chdir(tmp_path / "repo")
    output = tmp_path / "output"

    def run(disable=False):
        if disable:
            monkeypatch.setenv("RMG_KMC_DISABLE_ARTIFACT_CACHE", "1")
        else:
            monkeypatch.delenv("RMG_KMC_DISABLE_ARTIFACT_CACHE", raising=False)
        old_argv = sys.argv
        sys.argv = ["driver", str(database_path), str(output)]
        try:
            driver.main()
        finally:
            sys.argv = old_argv

    run()
    first = next(output.glob("*.json")).read_bytes()
    assert calls == ["generate"]
    run()
    assert calls == ["generate"]
    run(disable=True)
    assert calls == ["generate"]
    assert next(output.glob("*.json")).read_bytes() == first

    capsys.readouterr()
    monkeypatch.setenv("PYTHONHASHSEED", "random")
    monkeypatch.setattr(driver, "migrate", lambda *args: pytest.fail("random seed migrated cache"))
    monkeypatch.setattr(
        driver,
        "validate_artifact_cache_entry",
        lambda *args: pytest.fail("random seed read artifact cache"),
    )
    read_cache = []
    original_load_or_generate = driver.load_or_generate

    def load_or_generate(*args, **kwargs):
        read_cache.append(kwargs["read_cache"])
        return original_load_or_generate(*args, **kwargs)

    monkeypatch.setattr(driver, "load_or_generate", load_or_generate)
    calls.clear()
    run()
    assert read_cache == [False]
    assert calls == ["generate"]
    assert "cached public RMG generation" not in capsys.readouterr().out
    calls.clear()
    read_cache.clear()
    run()
    assert read_cache == [False]
    assert calls == ["generate"]
    assert "cached public RMG generation" not in capsys.readouterr().out


def test_driver_main_publishes_database_policy_and_provider_provenance(
    tmp_path, monkeypatch
):
    import compile_event_set_fixture as driver
    from rmgpy.kmc.barrier_e0 import FixedBBarrierE0Provider

    database_path = _make_database(tmp_path / "db")
    (database_path / "input/kinetics/families/stub").mkdir(parents=True)
    (database_path / "input/kinetics/families/stub/groups.py").write_text("")
    observed = []

    class FakeDatabase:
        def __init__(self):
            self.kinetics = type(
                "Kinetics", (),
                {"generate_reactions_from_families": lambda *args: []},
            )()
            self.thermo = object()

        def load_kinetics(self, *args, **kwargs):
            pass

        def load_thermo(self, *args, **kwargs):
            pass

    class FakeCompiler:
        @staticmethod
        def discover_family_reactions(*args, **kwargs):
            return {}, set(), {"stub": []}

        def __init__(self, *args, **kwargs):
            observed.append(kwargs)

        def write_artifact(self, directory):
            path = Path(directory) / "provenance.json"
            provider = observed[-1]["barrier_e0_provider"]
            path.write_text(json.dumps({
                "provenance": {
                    "database_sha": observed[-1]["rmg_database_sha"],
                    "applicability_policy_version": "persistent-carbene-applicability/1",
                    **({"barrier_e0_provider": provider.provenance} if provider else {}),
                }
            }))
            return path, None

    monkeypatch.setattr(driver, "RMGDatabase", FakeDatabase)
    monkeypatch.setattr(driver, "EventSetCompiler", FakeCompiler)
    monkeypatch.setattr(driver, "prepare_rate_rules", lambda *args, **kwargs: None)
    monkeypatch.setattr(driver, "ps_proxy_set", lambda units: ())
    monkeypatch.setattr(driver, "resolve_database_declaration", lambda value: value)
    monkeypatch.setattr(driver, "database_content_digest", lambda path: ("content", {}))
    monkeypatch.setattr(driver, "identity", lambda repo, db: {"identity": "driver"})
    monkeypatch.setattr(driver, "identity_name", lambda value: "driver")
    monkeypatch.setattr(driver, "artifact_cache_key", lambda *args: "artifact")
    monkeypatch.setattr(
        driver, "migrate", lambda *args: tmp_path / "cache" / "portable" / "driver"
    )
    monkeypatch.setenv("PYTHONHASHSEED", "0")
    monkeypatch.setenv("RMG_KMC_CACHE_ROOT", str(tmp_path / "cache"))
    monkeypatch.chdir(tmp_path)
    output = tmp_path / "output"
    old_argv = sys.argv
    sys.argv = [
        "driver", str(database_path), str(output),
        "--database-sha", "declared-db",
        "--barrier-e0-fixed-b", "900",
    ]
    try:
        driver.main()
    finally:
        sys.argv = old_argv
    artifact = json.loads(next(output.glob("*.json")).read_text())
    provenance = artifact["provenance"]
    assert provenance["database_sha"] == "declared-db"
    assert provenance["applicability_policy_version"] == (
        "persistent-carbene-applicability/1"
    )
    assert provenance["barrier_e0_provider"] == FixedBBarrierE0Provider(900).provenance

    observed.clear()
    output_off = tmp_path / "output-off"
    sys.argv = ["driver", str(database_path), str(output_off), "--database-sha", "declared-db"]
    try:
        driver.main()
    finally:
        sys.argv = old_argv
    off = json.loads(next(output_off.glob("*.json")).read_text())["provenance"]
    assert "barrier_e0_provider" not in off


def test_portable_cache_replays_in_fresh_depth_one_clone(tmp_path):
    repo = Path(__file__).resolve().parents[3]
    database = _make_database(tmp_path / "db")
    source = tmp_path / "source-cache"
    source.mkdir()
    identity_value = identity(repo, database)
    (source / "manifest.json").write_text(json.dumps(identity_value))
    artifact = {"artifact": "bounded-fixture", "source": "generation"}
    (source / "artifact.json").write_text(json.dumps(artifact, sort_keys=True))
    archive = tmp_path / "cache.tgz"
    export_cache(source, archive, repo, database)

    clone = tmp_path / "clone"
    subprocess.run(
        ["git", "clone", "--depth", "1", f"file://{repo}", str(clone)],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    imported = tmp_path / "imported"
    import_cache(archive, imported, clone, database)
    replay_script = (
        "import json,sys; from pathlib import Path; "
        "from portable_cache import load_or_generate; "
        "p=Path(sys.argv[1]); calls=[]; "
        "v=load_or_generate(p, lambda: calls.append(1) or {'artifact':'bounded-fixture'}, "
        "lambda x: json.dumps(x, sort_keys=True).encode(), json.loads); "
        "print(len(calls), json.dumps(v, sort_keys=True))"
    )
    env = {
        **os.environ,
        "PYTHONPATH": str(clone / "test/rmgpy/kmc"),
        "PYTHONDONTWRITEBYTECODE": "1",
    }
    replay = subprocess.check_output(
        [sys.executable, "-c", replay_script, str(imported / "artifact.json")],
        cwd=clone, env=env, text=True,
    ).strip()
    assert replay.startswith("0 ")
    assert hashlib.sha256((imported / "artifact.json").read_bytes()).hexdigest() == hashlib.sha256(
        (source / "artifact.json").read_bytes()
    ).hexdigest()

    (clone / "rmgpy/data/replay_marker.py").write_text("changed")
    miss = subprocess.check_output(
        [sys.executable, "-c", replay_script, str(tmp_path / "miss.json")],
        cwd=clone, env=env, text=True,
    ).strip()
    assert miss.startswith("1 ")


def test_new_rmg_kmc_environment_setting_changes_artifact_key(tmp_path, monkeypatch):
    repo = tmp_path / "repo"
    (repo / "rmgpy").mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / "rmgpy/x.py").write_text("x")
    for relative in ("test/rmgpy/kmc/compile_event_set_fixture.py", "test/rmgpy/kmc/portable_cache.py"):
        path = repo / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("driver")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    db = _make_database(tmp_path / "db")
    options = {"environment": compile_environment_options()}
    before = artifact_cache_key(repo, db, options)
    monkeypatch.setenv("RMG_KMC_NEW_COMPILER_SWITCH", "enabled")
    after = artifact_cache_key(repo, db, {"environment": compile_environment_options()})
    assert after != before
