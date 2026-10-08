"""Identity and content provenance for the RMG database inputs used by kMC.

The legacy ``rmg_database_sha`` field retains its old meaning only when it is
a clean, verified git HEAD (and, when supplied, matches the declaration).
``None`` means that the revision was not verified; consumers should use the
new declared/actual/tree/digest fields to explain why.
"""

from __future__ import annotations

import hashlib
import os
import subprocess
from pathlib import Path


DATABASE_IDENTITY_KEYS = (
    "rmg_database_sha_declared",
    "rmg_database_sha_actual",
    "rmg_database_tree_clean",
    "rmg_database_sha_matches_declared",
    "rmg_database_content_sha256",
    "rmg_database_content_files",
    "rmg_database_has_git_metadata",
    "rmg_database_sha",
)


def _input_roots(database_path: str | Path | None) -> tuple[Path, ...]:
    if database_path is None:
        return ()
    root = Path(database_path)
    candidates = (
        (root / "input" / "kinetics", root / "input" / "thermo"),
        (root / "kinetics", root / "thermo"),
    )
    for roots in candidates:
        if all(path.is_dir() for path in roots):
            return tuple(path for path in roots if path.is_dir())
    return ()


def database_content_digest(
    database_path: str | Path | None,
) -> tuple[str | None, tuple[str, ...]]:
    """Hash sorted relative paths and bytes of the kinetics/thermo input files.

    The input set is every regular file below ``input/kinetics`` and
    ``input/thermo`` (or, for a small standalone database, ``kinetics`` and
    ``thermo``), recursively.  ``__pycache__`` files are excluded.  The digest
    stream is ``relative POSIX path + NUL + file bytes + NUL`` per file in
    sorted path order; directory names and metadata never affect it.  Missing
    paths, unsupported layouts, empty trees, and symlinks are errors rather
    than silently receiving an identity for a different input set.
    """
    if database_path is None:
        return None, ()
    root = Path(database_path).resolve()
    if not root.exists():
        raise ValueError(f"database path does not exist: {database_path}")
    input_roots = _input_roots(root)
    if not input_roots:
        raise ValueError(
            f"database has no supported kinetics/thermo input layout: {database_path}"
        )
    files: list[tuple[str, Path]] = []

    def raise_walk_error(error):
        raise error

    for input_root in input_roots:
        if input_root.is_symlink():
            raise ValueError(f"symlink under database input roots: {input_root}")
        for current, directories, filenames in os.walk(
            input_root,
            followlinks=False,
            onerror=raise_walk_error,
        ):
            current_path = Path(current)
            for directory in directories:
                path = current_path / directory
                if path.is_symlink():
                    raise ValueError(f"symlink under database input roots: {path}")
            for filename in filenames:
                path = current_path / filename
                if path.is_symlink():
                    raise ValueError(f"symlink under database input roots: {path}")
                if "__pycache__" in path.parts:
                    continue
                relative = path.relative_to(root).as_posix()
                files.append((relative, path))
    files.sort(key=lambda item: item[0])
    if not files:
        raise ValueError(f"database input trees are empty: {database_path}")
    digest = hashlib.sha256()
    for relative, path in files:
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest(), tuple(relative for relative, _ in files)


def _git_identity(database_path: str | Path | None) -> tuple[str | None, bool | None, bool]:
    if database_path is None:
        return None, None, False
    path = Path(database_path).resolve()
    try:
        top_level = Path(
            subprocess.check_output(
                ["git", "-C", str(path), "rev-parse", "--show-toplevel"],
                text=True,
                stderr=subprocess.DEVNULL,
            ).strip()
        ).resolve()
        if top_level != path:
            return None, None, False
        subprocess.check_output(
            ["git", "-C", str(path), "rev-parse", "--is-inside-work-tree"],
            text=True,
            stderr=subprocess.DEVNULL,
        )
        actual = subprocess.check_output(
            ["git", "-C", str(path), "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
        status = subprocess.check_output(
            ["git", "-C", str(path), "status", "--porcelain", "--untracked-files=all", "--"],
            text=True,
            stderr=subprocess.DEVNULL,
        )
        return actual, not status, True
    except (OSError, subprocess.CalledProcessError):
        return None, None, False


def database_provenance(database_path: str | Path | None, declared: str | None) -> dict:
    """Return one authoritative database identity for compiler and thermo."""
    actual, clean, has_git = _git_identity(database_path)
    digest, files = database_content_digest(database_path)
    verified = (
        actual is not None
        and clean is True
        and (declared is None or actual == declared)
    )
    return {
        "rmg_database_sha_declared": declared,
        "rmg_database_sha_actual": actual,
        "rmg_database_tree_clean": clean,
        "rmg_database_sha_matches_declared": (
            actual == declared if actual is not None and declared is not None else None
        ),
        "rmg_database_content_sha256": digest,
        "rmg_database_content_files": list(files),
        "rmg_database_has_git_metadata": has_git,
        # Kept for old artifact validators only when the revision was checked.
        "rmg_database_sha": actual if verified else None,
    }


def check_database_identity(candidate: dict, expected: dict) -> None:
    """Reject an injected assignment/provider carrying a different identity."""
    for key in DATABASE_IDENTITY_KEYS:
        if key not in candidate or candidate[key] != expected.get(key):
            raise ValueError(f"injected thermo database identity differs: {key}")


def reject_external_library_paths(
    database_path: str | Path | None, *databases: object
) -> None:
    """Reject configured library files outside the hashed database tree."""
    if database_path is None:
        return
    root = Path(database_path).resolve()
    candidates: list[object] = []
    recorded_sources: list[object] = []
    for database in databases:
        if database is None:
            continue
        candidates.extend(getattr(database, "library_order", ()) or ())
        candidates.extend((getattr(database, "external_library_labels", {}) or {}).keys())
        for library in (getattr(database, "libraries", {}) or {}).values():
            recorded_source = getattr(library, "source_path", None)
            if recorded_source is not None:
                recorded_sources.append(recorded_source)
            candidates.extend(
                getattr(library, attribute, None)
                for attribute in ("path", "file", "filename")
            )
    allowed_roots = tuple(path.resolve() for path in _input_roots(root))

    def is_under(path: Path, parent: Path) -> bool:
        try:
            path.relative_to(parent)
            return True
        except ValueError:
            return False

    pending = [*candidates, *recorded_sources]
    while pending:
        candidate = pending.pop(0)
        if isinstance(candidate, tuple):
            pending.extend(candidate[:1])
            continue
        if not isinstance(candidate, (str, Path)):
            continue
        path = Path(candidate).expanduser()
        if not path.exists():
            if candidate in recorded_sources:
                raise ValueError(f"recorded loaded library source is missing: {candidate}")
            continue
        resolved = path.resolve()
        if not any(is_under(resolved, allowed) for allowed in allowed_roots):
            raise ValueError(
                f"kinetics/thermo library is outside hashed input roots: {candidate}"
            ) from None


def provenance_matches_database(
    provenance: dict,
    database_path: str | Path | None,
    declared: str | None,
) -> bool:
    """Match fresh declared+digest identity, retaining old-artifact fallback."""
    current = database_provenance(database_path, declared)
    if "rmg_database_content_sha256" in provenance:
        return (
            provenance.get("rmg_database_sha_declared")
            == current["rmg_database_sha_declared"]
            and provenance.get("rmg_database_content_sha256")
            == current["rmg_database_content_sha256"]
        )
    expected_legacy = current["rmg_database_sha"] or declared
    return provenance.get("rmg_database_sha") == expected_legacy


def resolve_database_declaration(explicit: str | None = None) -> str | None:
    """Resolve one database declaration for producers, consumers, and caches."""
    return explicit if explicit is not None else os.environ.get("RMG_DATABASE_SHA") or None


def legacy_identity_provenance(database_commit: str | None) -> dict:
    """Build a complete identity shape for standalone thermo test providers."""
    identity = {key: None for key in DATABASE_IDENTITY_KEYS}
    identity["rmg_database_sha"] = database_commit
    return identity
