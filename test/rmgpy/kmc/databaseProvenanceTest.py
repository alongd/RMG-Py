"""Database revision and content identity tests for kMC artifacts."""

import hashlib
import shutil
import subprocess
from pathlib import Path

import pytest

from rmgpy.kmc.compiler import EventSetCompiler
from rmgpy.kmc.database_provenance import (
    DATABASE_IDENTITY_KEYS,
    database_content_digest,
    database_provenance,
    provenance_matches_database,
)


def _database(root):
    kinetics = root / "input" / "kinetics" / "families"
    thermo = root / "input" / "thermo" / "libraries"
    kinetics.mkdir(parents=True)
    thermo.mkdir(parents=True)
    (kinetics / "rules.py").write_text("rule = 1\n")
    (thermo / "primary.py").write_text("thermo = 1\n")


def _git_database(root):
    _database(root)
    subprocess.run(["git", "init", "-q", str(root)], check=True)
    subprocess.run(["git", "-C", str(root), "add", "."], check=True)
    subprocess.run(
        ["git", "-C", str(root), "-c", "user.name=Test", "-c", "user.email=test@example.invalid", "commit", "-qm", "initial"],
        check=True,
    )
    return subprocess.check_output(
        ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
    ).strip()


def test_git_database_records_clean_head_and_declared_match(tmp_path):
    root = tmp_path / "database"
    head = _git_database(root)

    identity = database_provenance(root, head)

    assert identity["rmg_database_sha_declared"] == head
    assert identity["rmg_database_sha_actual"] == head
    assert identity["rmg_database_tree_clean"] is True
    assert identity["rmg_database_sha_matches_declared"] is True
    assert identity["rmg_database_sha"] == head


def test_git_database_records_dirty_tree_and_declared_mismatch(tmp_path):
    root = tmp_path / "database"
    head = _git_database(root)
    (root / "input" / "kinetics" / "families" / "rules.py").write_text("rule = 2\n")

    identity = database_provenance(root, "declared-but-not-head")

    assert identity["rmg_database_sha_actual"] == head
    assert identity["rmg_database_tree_clean"] is False
    assert identity["rmg_database_sha_matches_declared"] is False
    assert identity["rmg_database_sha"] is None


def test_database_copy_without_git_has_stable_content_digest(tmp_path):
    source = tmp_path / "source"
    copy = tmp_path / "copy"
    _database(source)
    shutil.copytree(source, copy)

    first, files = database_content_digest(copy)
    second, copied_files = database_content_digest(copy)
    identity = database_provenance(copy, "declared-copy")

    assert first == second == identity["rmg_database_content_sha256"]
    assert files == copied_files == (
        "input/kinetics/families/rules.py",
        "input/thermo/libraries/primary.py",
    )
    assert identity["rmg_database_sha_actual"] is None
    assert identity["rmg_database_tree_clean"] is None
    assert identity["rmg_database_has_git_metadata"] is False
    assert identity["rmg_database_sha"] is None


def test_database_digest_changes_when_kinetics_input_changes(tmp_path):
    root = tmp_path / "database"
    _database(root)
    before, _ = database_content_digest(root)
    (root / "input" / "kinetics" / "families" / "rules.py").write_text("rule = 2\n")

    after, _ = database_content_digest(root)

    assert before != after
    assert after == hashlib.sha256(
        b"input/kinetics/families/rules.py\0rule = 2\n\0"
        b"input/thermo/libraries/primary.py\0thermo = 1\n\0"
    ).hexdigest()


def test_database_nested_in_unrelated_git_repo_has_no_git_identity(tmp_path):
    outer = tmp_path / "outer"
    database = outer / "database"
    _database(database)
    subprocess.run(["git", "init", "-q", str(outer)], check=True)
    subprocess.run(["git", "-C", str(outer), "add", "."], check=True)
    subprocess.run(
        ["git", "-C", str(outer), "-c", "user.name=Test", "-c", "user.email=test@example.invalid", "commit", "-qm", "outer"],
        check=True,
    )

    identity = database_provenance(database, None)

    assert identity["rmg_database_sha_actual"] is None
    assert identity["rmg_database_has_git_metadata"] is False
    assert identity["rmg_database_sha"] is None


@pytest.mark.parametrize("database_path", [None])
def test_digest_is_unknown_without_a_database_path(database_path):
    assert database_content_digest(database_path) == (None, ())


@pytest.mark.parametrize("kind", ["missing", "unsupported", "empty"])
def test_digest_rejects_paths_without_input_files(tmp_path, kind):
    path = tmp_path / kind
    if kind == "unsupported":
        (path / "other").mkdir(parents=True)
    elif kind == "empty":
        (path / "input" / "kinetics").mkdir(parents=True)
        (path / "input" / "thermo").mkdir()
    with pytest.raises(ValueError):
        database_content_digest(path)


def test_digest_rejects_input_symlink(tmp_path):
    root = tmp_path / "database"
    _database(root)
    (root / "input" / "kinetics" / "external.py").symlink_to(
        root / "input" / "thermo" / "libraries" / "primary.py"
    )

    with pytest.raises(ValueError, match="symlink"):
        database_content_digest(root)


def test_compiler_rejects_external_library_path(tmp_path):
    root = tmp_path / "database"
    _database(root)
    external = tmp_path / "external.py"
    external.write_text("external = True\n")
    database = type(
        "Database",
        (),
        {"families": {}, "library_order": [str(external)]},
    )()

    with pytest.raises(ValueError, match="outside hashed input roots"):
        EventSetCompiler(
            kinetics_database=database,
            proxies=[],
            families=[],
            database_path=root,
            thermo_database=None,
        )


def test_compiler_and_shared_thermo_provenance_are_identical(tmp_path):
    root = tmp_path / "database"
    head = _git_database(root)
    compiler = EventSetCompiler(
        kinetics_database=type("Database", (), {"families": {}})(),
        proxies=[],
        families=[],
        database_path=root,
        rmg_database_sha=head,
        thermo_database=None,
    )

    artifact = compiler.compile()
    expected = {
        key: compiler.database_provenance[key] for key in DATABASE_IDENTITY_KEYS
    }
    artifact_identity = {
        key: artifact["provenance"][key] for key in DATABASE_IDENTITY_KEYS
    }
    thermo_identity = {
        key: compiler.reference_thermo_provider.provenance[key]
        for key in DATABASE_IDENTITY_KEYS
    }

    assert artifact_identity == expected
    assert thermo_identity == expected


def test_injected_assignment_with_different_identity_is_rejected(tmp_path):
    root = tmp_path / "database"
    head = _git_database(root)
    from rmgpy.kmc.reference_thermo import SharedThermoAssignment

    assignment = SharedThermoAssignment(
        None,
        "wrong",
    )
    provider = type("Provider", (), {"assignment": assignment, "provenance": assignment.provenance})()
    with pytest.raises(ValueError, match="identity differs"):
        EventSetCompiler(
            kinetics_database=type("Database", (), {"families": {}})(),
            proxies=[],
            families=[],
            database_path=root,
            rmg_database_sha=head,
            thermo_database=None,
            reference_thermo_provider=provider,
        )


@pytest.mark.parametrize("provider_kind", ["null-sha", "empty"])
def test_injected_identity_is_complete_for_gitless_database(tmp_path, provider_kind):
    root = tmp_path / "database"
    _database(root)
    from rmgpy.kmc.reference_thermo import SharedThermoAssignment

    if provider_kind == "null-sha":
        assignment = SharedThermoAssignment(None, None)
        provider = type(
            "Provider", (), {"assignment": assignment, "provenance": assignment.provenance}
        )()
    else:
        provider = type("Provider", (), {"provenance": {}})()
    with pytest.raises(ValueError, match="identity differs"):
        EventSetCompiler(
            kinetics_database=type("Database", (), {"families": {}})(),
            proxies=[],
            families=[],
            database_path=root,
            thermo_database=None,
            reference_thermo_provider=provider,
        )


def test_external_loaded_library_objects_are_rejected(tmp_path):
    root = tmp_path / "database"
    _database(root)
    from rmgpy.data.kinetics.database import KineticsDatabase
    from rmgpy.data.thermo import ThermoDatabase

    external_thermo = tmp_path / "external_thermo.py"
    shutil.copy(
        "test/rmgpy/test_data/testing_database/thermo/libraries/primaryThermoLibrary.py",
        external_thermo,
    )
    thermo = ThermoDatabase()
    thermo.load_libraries(
        str(root / "input/thermo/libraries"), [str(external_thermo)]
    )
    external_kinetics = tmp_path / "external_kinetics"
    shutil.copytree(
        "test/rmgpy/test_data/testing_database/kinetics/libraries/lib_net",
        external_kinetics,
    )
    kinetics = KineticsDatabase()
    kinetics.load_libraries(
        str(root / "input/kinetics/libraries"), [str(external_kinetics)]
    )
    kinetics.library_order = [(external_kinetics.name, "Reaction Library")]

    with pytest.raises(ValueError, match="outside hashed input roots"):
        EventSetCompiler(
            kinetics_database=kinetics,
            proxies=[],
            families=[],
            database_path=root,
            thermo_database=thermo,
        )

    external_thermo.unlink()
    shutil.rmtree(external_kinetics)
    with pytest.raises(ValueError, match="recorded loaded library source is missing"):
        EventSetCompiler(
            kinetics_database=kinetics,
            proxies=[],
            families=[],
            database_path=root,
            thermo_database=thermo,
        )


def test_consumer_identity_accepts_fresh_digest_and_rejects_mismatch(tmp_path):
    root = tmp_path / "database"
    _database(root)
    from rmgpy.kmc.database_provenance import provenance_matches_database

    digest, _ = database_content_digest(root)
    fresh = {
        "rmg_database_sha": None,
        "rmg_database_sha_declared": "declared-copy",
        "rmg_database_content_sha256": digest,
    }
    mismatch = {**fresh, "rmg_database_content_sha256": "wrong"}

    assert provenance_matches_database(fresh, root, "declared-copy")
    assert not provenance_matches_database(mismatch, root, "declared-copy")


def test_fixture_declaration_path_accepts_gitless_stub_artifact(tmp_path, monkeypatch):
    root = tmp_path / "database"
    _database(root)
    digest, _ = database_content_digest(root)
    from rmgpy.kmc.database_provenance import resolve_database_declaration

    for declared in (None, "declared-copy"):
        if declared is None:
            monkeypatch.delenv("RMG_DATABASE_SHA", raising=False)
        else:
            monkeypatch.setenv("RMG_DATABASE_SHA", declared)
        resolved = resolve_database_declaration()
        artifact = {
            "provenance": {
                "rmg_database_sha": None,
                "rmg_database_sha_declared": resolved,
                "rmg_database_content_sha256": digest,
            }
        }
        assert provenance_matches_database(
            artifact["provenance"], root, resolved
        )


def test_generation_cache_key_changes_with_database_contents(tmp_path, monkeypatch):
    root = tmp_path / "database"
    _database(root)
    monkeypatch.delenv("RMG_DATABASE_SHA", raising=False)
    from completeness_oracle import oracle_cache_key

    repository = Path(__file__).resolve().parents[3]
    before = oracle_cache_key(repository, root, 1)
    (root / "input/kinetics/families/rules.py").write_text("rule = changed\n")

    assert oracle_cache_key(repository, root, 1) != before
