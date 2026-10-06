import json
import subprocess

import pytest

from portable_cache import export_cache, identity, identity_name, import_cache, migrate


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
    old = tmp_path / "generated-reactions/old/0"
    old.mkdir(parents=True)
    (old / "x.pickle").write_bytes(b"x")
    target = migrate(tmp_path, repo, db)
    assert (target / "generated-reactions/x.pickle").read_bytes() == b"x"
    assert (old / "x.pickle").exists()
