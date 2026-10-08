import json
import os
from pathlib import Path
import shutil
import subprocess
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
