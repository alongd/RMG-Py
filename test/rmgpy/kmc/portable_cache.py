"""Portable identity, migration, and archive operations for kMC caches."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile

from rmgpy.kmc.database_provenance import database_content_digest

SCHEMA = 1
ARTIFACT_ENV_EXCLUSIONS = frozenset({
    "RMG_KMC_CACHE_ROOT",
    "RMG_KMC_ARTIFACT",
})


def _run(*args: str, cwd: Path | None = None) -> str:
    return subprocess.check_output(args, cwd=cwd, text=True, stderr=subprocess.DEVNULL).strip()


def tree_identity(repository: Path) -> str:
    """Hash tracked rmgpy content, excluding the kMC implementation itself."""
    entries = []
    tracked = _run("git", "-C", str(repository), "ls-files", "rmgpy").splitlines()
    for relative in tracked:
        if relative.startswith("rmgpy/kmc/"):
            continue
        path = repository / relative
        entries.append((path.relative_to(repository).as_posix(), hashlib.sha256(path.read_bytes()).hexdigest()))
    return hashlib.sha256(json.dumps(entries, sort_keys=True).encode()).hexdigest()


def database_identity(database: Path) -> str:
    digest, _ = database_content_digest(database)
    if digest is None:
        raise ValueError(f"database content digest unavailable: {database}")
    return digest


def identity(repository: Path, database: Path) -> dict:
    import rdkit
    hash_seed = os.environ.get("PYTHONHASHSEED")
    return {
        "schema": SCHEMA,
        "rmgpy_tree_sha256": tree_identity(repository),
        "database": database_identity(database),
        "rmg_database_content_sha256": database_identity(database),
        "python": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        "rdkit": getattr(rdkit, "__version__", "unknown"),
        "pythonhashseed": hash_seed,
    }


def identity_name(value: dict) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def atomic_copy(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.{os.getpid()}.tmp")
    try:
        shutil.copy2(source, temporary)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        temporary.write_bytes(payload)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def load_or_generate(path: Path, generate, dump, load):
    """Load a migrated generation entry or call the expensive producer once."""
    if path.is_file():
        return load(path.read_bytes())
    value = generate()
    atomic_write(path, dump(value))
    return value


def _file_hashes(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(root.rglob("*"))
        if path.is_file() and not path.is_symlink() and path != root / "manifest.json"
    }


def _source_hash(repository: Path, relative_paths: tuple[str, ...]) -> str:
    entries = []
    for relative in relative_paths:
        path = repository / relative
        entries.append((relative, hashlib.sha256(path.read_bytes()).hexdigest()))
    return hashlib.sha256(json.dumps(entries, sort_keys=True).encode()).hexdigest()


def artifact_identity(repository: Path, database: Path) -> dict:
    """Artifact identity additionally includes the kMC compiler implementation."""
    value = identity(repository, database)
    tracked = _run("git", "-C", str(repository), "ls-files", "rmgpy/kmc").splitlines()
    value["kmc_tree_sha256"] = hashlib.sha256(json.dumps(
        [(p, hashlib.sha256((repository / p).read_bytes()).hexdigest()) for p in tracked],
        sort_keys=True,
    ).encode()).hexdigest()
    value["compile_driver_sha256"] = _source_hash(repository, (
        "test/rmgpy/kmc/compile_event_set_fixture.py",
        "test/rmgpy/kmc/portable_cache.py",
    ))
    return value


def artifact_cache_key(repository: Path, database: Path, compile_options: dict) -> str:
    return identity_name({
        **artifact_identity(repository, database),
        "compile_options": compile_options,
    })


def validate_artifact_cache_entry(
    manifest_path: Path,
    expected_key: str,
    database: Path,
    database_declaration: str | None,
) -> Path:
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("identity") != expected_key:
        raise ValueError("identity mismatch")
    artifact_path = manifest_path.parent / manifest["artifact"]
    payload = artifact_path.read_bytes()
    if manifest.get("artifact_sha256") != hashlib.sha256(payload).hexdigest():
        raise ValueError("artifact hash mismatch")
    artifact = json.loads(payload)
    from rmgpy.kmc.compiler import compiler_source_hash, validate_artifact
    from rmgpy.kmc.database_provenance import provenance_matches_database
    validate_artifact(artifact)
    if artifact["provenance"].get("compiler_sources_sha256") != compiler_source_hash():
        raise ValueError("compiler provenance mismatch")
    if not provenance_matches_database(
        artifact["provenance"], database, database_declaration
    ):
        raise ValueError("database provenance mismatch")
    return artifact_path


def compile_environment_options() -> dict[str, str]:
    """Return every content-affecting RMG_KMC_ setting in stable order."""
    return {
        name: value
        for name, value in sorted(os.environ.items())
        if name.startswith("RMG_KMC_") and name not in ARTIFACT_ENV_EXCLUSIONS
    }


def tree_identity_at_commit(repository: Path, commit: str) -> str:
    """Hash the committed rmgpy blobs, excluding rmgpy/kmc, via Git objects."""
    entries = []
    listing = _run("git", "-C", str(repository), "ls-tree", "-r", "-z", commit, "--", "rmgpy")
    for item in listing.split("\0"):
        if not item:
            continue
        metadata, relative = item.split("\t", 1)
        if relative.startswith("rmgpy/kmc/"):
            continue
        blob = metadata.split()[2]
        content = subprocess.check_output(
            ["git", "-C", str(repository), "cat-file", "blob", blob]
        )
        entries.append((relative, hashlib.sha256(content).hexdigest()))
    return hashlib.sha256(json.dumps(entries, sort_keys=True).encode()).hexdigest()


def migrate(
    cache_root: Path,
    repository: Path,
    database: Path,
    database_sha: str | None = None,
) -> Path:
    generation_identity = identity(repository, database)
    target = cache_root / "portable" / identity_name(generation_identity)
    target.mkdir(parents=True, exist_ok=True)
    atomic_write(target / "manifest.json", (json.dumps(generation_identity, sort_keys=True, indent=2) + "\n").encode())
    old = cache_root / "generated-reactions"
    out = target / "generated-reactions"
    out.mkdir(exist_ok=True)
    current_tree = generation_identity["rmgpy_tree_sha256"]
    current_database = generation_identity["database"]
    current_seed = os.environ.get("PYTHONHASHSEED")
    if current_seed is None:
        print("skipping legacy migration: PYTHONHASHSEED is unset")
        return target
    for directory in old.iterdir() if old.is_dir() else ():
        if not directory.is_dir():
            continue
        try:
            origin_commit, declaration, old_database = directory.name.rsplit("-", 2)
        except ValueError:
            print(f"skipping legacy cache without database content digest: {directory}")
            continue
        if old_database != current_database:
            print(f"skipping cache with mismatched database content: {directory}")
            continue
        try:
            if tree_identity_at_commit(repository, origin_commit) != current_tree:
                print(f"skipping cache with mismatched rmgpy tree: {directory}")
                continue
        except subprocess.CalledProcessError:
            print(f"skipping cache with missing origin commit: {directory}")
            continue
        for seed in directory.iterdir():
            if seed.name != current_seed:
                print(f"skipping cache with mismatched hash seed: {directory}/{seed.name}")
                continue
            if seed.is_dir():
                for entry in seed.iterdir():
                    if entry.is_file() and entry.suffix == ".pickle":
                        destination = out / entry.name
                        if not destination.exists():
                            atomic_copy(entry, destination)
    return target


def export_cache(cache: Path, archive: Path, repository: Path, database: Path) -> None:
    expected = identity(repository, database)
    source_manifest_path = cache / "manifest.json"
    if not source_manifest_path.is_file():
        raise ValueError("cache bucket has no identity manifest")
    source_manifest = json.loads(source_manifest_path.read_text())
    if source_manifest != expected:
        raise ValueError("cache bucket identity does not match this checkout")
    manifest = dict(source_manifest)
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary) / "kmc-cache"
        shutil.copytree(cache, root)
        manifest["files"] = _file_hashes(root)
        (root / "manifest.json").write_text(json.dumps(manifest, sort_keys=True, indent=2) + "\n")
        archive.parent.mkdir(parents=True, exist_ok=True)
        fd, temporary_archive = tempfile.mkstemp(prefix=f".{archive.name}.", dir=archive.parent)
        os.close(fd)
        try:
            with tarfile.open(temporary_archive, "w:gz") as tar:
                tar.add(root, arcname="kmc-cache")
            os.replace(temporary_archive, archive)
        finally:
            Path(temporary_archive).unlink(missing_ok=True)


def import_cache(archive: Path, destination: Path, repository: Path, database: Path) -> None:
    expected = identity(repository, database)
    with tempfile.TemporaryDirectory() as temporary:
        with tarfile.open(archive, "r:gz") as tar:
            tar.extractall(temporary, filter="data")
        root = Path(temporary) / "kmc-cache"
        actual = json.loads((root / "manifest.json").read_text())
        if {key: value for key, value in actual.items() if key != "files"} != expected:
            raise ValueError("cache manifest identity does not match this checkout")
        files = actual.get("files")
        observed = _file_hashes(root)
        if not isinstance(files, dict) or set(files) != set(observed):
            raise ValueError("cache archive contains unknown or missing files")
        if any(observed[path] != digest for path, digest in files.items()):
            raise ValueError("cache archive file hash mismatch")
        destination.mkdir(parents=True, exist_ok=True)
        for item in root.iterdir():
            if item.name != "manifest.json":
                target = destination / item.name
                if item.is_dir():
                    for source in item.rglob("*"):
                        if source.is_file():
                            atomic_copy(source, destination / source.relative_to(root))
                else:
                    atomic_copy(item, target)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("identity", "migrate", "export", "import"))
    parser.add_argument("cache")
    parser.add_argument("database")
    parser.add_argument("extra", nargs="?")
    args = parser.parse_args()
    repository = Path.cwd()
    database = Path(args.database)
    cache = Path(args.cache)
    if args.command == "identity":
        print(json.dumps(identity(repository, database), sort_keys=True, indent=2))
    elif args.command == "migrate":
        print(migrate(cache, repository, database))
    elif args.command == "export":
        export_cache(cache, Path(args.extra), repository, database)
    else:
        import_cache(Path(args.extra), cache, repository, database)


if __name__ == "__main__":
    main()
