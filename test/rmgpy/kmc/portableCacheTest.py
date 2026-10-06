import json
import os
import shutil
import subprocess

import pytest

from portable_cache import (
    artifact_cache_key,
    database_identity,
    export_cache,
    identity,
    import_cache,
    migrate,
)


def test_identity_excludes_kmc_and_changes_for_other_rmgpy(tmp_path, monkeypatch):
    repo = tmp_path / "repo"
    (repo / "rmgpy/kmc").mkdir(parents=True)
    (repo / "rmgpy/data").mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / "rmgpy/kmc/x.py").write_text("one")
    (repo / "rmgpy/data/x.py").write_text("one")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    db = tmp_path / "db"
    (db / "input").mkdir(parents=True)
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
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    db = tmp_path / "db"
    (db / "input").mkdir(parents=True)
    source = tmp_path / "cache"
    (source / "generated-reactions").mkdir(parents=True)
    (source / "generated-reactions/entry.pickle").write_bytes(b"entry")
    archive = tmp_path / "cache.tgz"
    export_cache(source, archive, repo, db)
    destination = tmp_path / "destination"
    import_cache(archive, destination, repo, db)
    assert (destination / "generated-reactions/entry.pickle").read_bytes() == b"entry"
    (repo / "rmgpy/x.py").write_text("changed")
    with pytest.raises(ValueError, match="manifest"):
        import_cache(archive, tmp_path / "refused", repo, db)


def test_migration_copies_old_named_entries(tmp_path):
    repo = tmp_path / "repo"
    (repo / "rmgpy").mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / "rmgpy/x.py").write_text("x")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    db = tmp_path / "db"
    (db / "input").mkdir(parents=True)
    origin = subprocess.check_output(
        ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True
    ).strip()
    seed = os.environ.get("PYTHONHASHSEED", "default")
    old = tmp_path / f"generated-reactions/{origin}-{database_identity(db)}/{seed}"
    old.mkdir(parents=True)
    (old / "x.pickle").write_bytes(b"x")
    target = migrate(tmp_path, repo, db)
    assert (target / "generated-reactions/x.pickle").read_bytes() == b"x"
    assert (old / "x.pickle").exists()


def test_migration_skips_mismatched_database(tmp_path):
    repo = tmp_path / "repo"
    (repo / "rmgpy").mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / "rmgpy/x.py").write_text("x")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    db = tmp_path / "db"
    (db / "input").mkdir(parents=True)
    origin = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    wrong = tmp_path / f"generated-reactions/{origin}-wrong-db/0"
    wrong.mkdir(parents=True)
    (wrong / "x.pickle").write_bytes(b"x")
    target = migrate(tmp_path, repo, db)
    assert not (target / "generated-reactions/x.pickle").exists()


def test_git_and_archive_copies_have_same_database_identity(tmp_path):
    git_db = tmp_path / "git-db"
    (git_db / "input/sub").mkdir(parents=True)
    (git_db / "input/sub/data.txt").write_text("same")
    subprocess.run(["git", "init", "-q", str(git_db)], check=True)
    subprocess.run(["git", "-C", str(git_db), "add", "."], check=True)
    subprocess.run(["git", "-C", str(git_db), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    archive_db = tmp_path / "archive-db"
    shutil.copytree(git_db / "input", archive_db / "input")
    assert database_identity(git_db) == database_identity(archive_db)


def test_artifact_key_includes_compile_options(tmp_path):
    repo = tmp_path / "repo"
    (repo / "rmgpy").mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / "rmgpy/x.py").write_text("x")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b", "-c", "user.name=a", "commit", "-qm", "x"], check=True)
    db = tmp_path / "db"
    (db / "input").mkdir(parents=True)
    base = {"temperature_grid": [300, 325], "use_plpsec_library": "1"}
    disabled = {**base, "use_plpsec_library": "0"}
    assert artifact_cache_key(repo, db, base) != artifact_cache_key(repo, db, disabled)
