"""Identities for the generated-reaction and compilation caches.

Every identity is a full SHA-256 over canonical JSON (sorted keys, no
whitespace, floats as ``float.hex``) that carries a domain tag and a schema
version.  Where a choice exists the hash covers more rather than less:
over-invalidation costs a recompute, a collision serves wrong chemistry.

Identities provided here
------------------------
* ``generation_request_identity``: one public generator call, over the
  observable input state (not a canonical form that discards it).
* ``compilation_request_identity``: one ``CompileRequest`` with ``output_dir``
  and the database location excluded and every default resolved.
* ``published_content_sha``: the existing artifact content hash.
* ``generation_dependency_manifest`` / ``compilation_dependency_manifest``: two
  fixed, explicit manifests of interpreter, libraries and RMG-Py bytes.
* ``database_snapshot_manifest`` / ``verify_database_snapshot``: a content
  manifest of the database input trees and its re-verification.
* ``build_envelope`` / ``verify_envelope``: the request-binding envelope.

Nothing here reads, writes or locks a cache.  Storage and reuse are built on
top of these functions.

Items that could not be pinned robustly are listed in ``UNPINNED_ITEMS``.
"""

from __future__ import annotations

import functools
import hashlib
import importlib.util
import json
import math
import os
import platform
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence

SCHEMA_VERSION = 1

UNPINNED_ITEMS = """\
Items the manifests do not (or cannot) pin robustly:

* Native libraries below the Python wheels: libstdc++/libc, BLAS/LAPACK (only
  the NumPy/SciPy build configuration strings are recorded), Boost, the RDKit
  shared objects (only RDKit's version, Boost version and build string),
  DASSL/DASPK/Fortran runtimes.  A conda environment rebuilt with different
  binaries but identical version strings is not detected.
* Loaded-versus-disk drift: bytes are read from disk when a manifest is first
  built in a process (then cached for that process).  A module imported before
  an on-disk edit, or edited after the first manifest call, is hashed as the
  disk bytes at that moment, not as the loaded code.  The compiler's own
  ``compiler_sources_sha256`` has the opposite property (hash at import) and
  is hashed into the compilation identity as well.
* Modules outside the two fixed lists: anything the lists do not name is not
  hashed even when it is imported (the lists are explicit by design).
* The git sha and dirty flag cover tracked files only; untracked files in the
  RMG-Py checkout are not part of the dirty flag.
* CPU features and thread counts: results are assumed independent of them.
* ``PYTHONHASHSEED`` is checked to be a fixed integer; set iteration order
  that depends on object addresses rather than string hashes is not detected.
* The generator reads the process-global database registry
  (``rmgpy.data.rmg.get_db('kinetics')``, ``rmgpy/data/kinetics/common.py``)
  in one path; the family fingerprint covers the database object passed to
  the call, so a different global database is not detected.
"""


class IdentityError(Exception):
    """Base class of identity failures."""


class UnstableHashSeedError(IdentityError):
    """``PYTHONHASHSEED`` is unset or random, so no reusable identity exists."""


class SnapshotMismatchError(IdentityError):
    """A database directory no longer matches its manifest."""


class EnvelopeError(IdentityError):
    """An envelope is malformed or bound to a different request."""


# ---------------------------------------------------------------------------
# Canonical encoding
# ---------------------------------------------------------------------------


def _canon(value: Any) -> Any:
    """Convert ``value`` to a JSON tree with exact floats and tagged non-JSON types."""
    if value is None or type(value) is bool or type(value) is str:
        return value
    if type(value) is int:
        return value
    if type(value) is float:
        if not math.isfinite(value):
            raise IdentityError(f"non-finite float {value!r} cannot be part of an identity")
        return {"$float": value.hex()}
    if isinstance(value, (list, tuple)):
        return [_canon(item) for item in value]
    if isinstance(value, dict):
        result = {}
        for key, item in value.items():
            if type(key) is not str:
                raise IdentityError(f"non-string key {key!r} cannot be part of an identity")
            result[key] = _canon(item)
        return result
    if hasattr(value, "item") and type(value).__module__ == "numpy":
        return _canon(value.item())
    raise IdentityError(f"unsupported value of type {type(value).__name__} in an identity")


def canonical_bytes(value: Any) -> bytes:
    """Canonical JSON bytes of ``value`` (floats as ``float.hex``)."""
    return json.dumps(
        _canon(value), sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False
    ).encode("ascii")


def sha256_canonical(kind: str, value: Any) -> str:
    """Full SHA-256 hex of ``value`` under domain tag ``kind`` and the schema version."""
    envelope = {"kind": kind, "schema_version": SCHEMA_VERSION, "value": value}
    return hashlib.sha256(canonical_bytes(envelope)).hexdigest()


def reset_caches() -> None:
    """Forget per-process manifest caches (tests, long-lived workers)."""
    _cached_manifest.cache_clear()


# ---------------------------------------------------------------------------
# Hash seed
# ---------------------------------------------------------------------------


def require_stable_hash_seed() -> str:
    """Return ``PYTHONHASHSEED`` or refuse: unset/random seeds give no reusable identity."""
    value = os.environ.get("PYTHONHASHSEED")
    if value is None or not re.fullmatch(r"[0-9]{1,10}", value) or int(value) > 4294967295:
        raise UnstableHashSeedError(
            "PYTHONHASHSEED must be set to a fixed integer in [0, 4294967295] "
            f"to produce a reusable identity (got {value!r})"
        )
    return value


# ---------------------------------------------------------------------------
# Dependency manifests
# ---------------------------------------------------------------------------

GENERATION_MODULES: tuple[str, ...] = (
    "rmgpy",
    "rmgpy.constants",
    "rmgpy.constraints",
    "rmgpy.exceptions",
    "rmgpy.quantity",
    "rmgpy.reaction",
    "rmgpy.species",
    "rmgpy.data.base",
    "rmgpy.data.rmg",
    "rmgpy.data.kinetics.common",
    "rmgpy.data.kinetics.database",
    "rmgpy.data.kinetics.depository",
    "rmgpy.data.kinetics.family",
    "rmgpy.data.kinetics.groups",
    "rmgpy.data.kinetics.library",
    "rmgpy.data.kinetics.rules",
    "rmgpy.kinetics.arrhenius",
    "rmgpy.kinetics.model",
    "rmgpy.kinetics.uncertainties",
    "rmgpy.molecule.adjlist",
    "rmgpy.molecule.atomtype",
    "rmgpy.molecule.converter",
    "rmgpy.molecule.element",
    "rmgpy.molecule.filtration",
    "rmgpy.molecule.fragment",
    "rmgpy.molecule.graph",
    "rmgpy.molecule.group",
    "rmgpy.molecule.inchi",
    "rmgpy.molecule.kekulize",
    "rmgpy.molecule.molecule",
    "rmgpy.molecule.pathfinder",
    "rmgpy.molecule.resonance",
    "rmgpy.molecule.symmetry",
    "rmgpy.molecule.translator",
    "rmgpy.molecule.util",
    "rmgpy.molecule.vf2",
    "rmgpy.rmg.input",
    "rmgpy.thermo.model",
    "rmgpy.thermo.nasa",
    "rmgpy.thermo.thermodata",
    "rmgpy.thermo.wilhoit",
)

COMPILATION_EXTRA_MODULES: tuple[str, ...] = (
    "rmgpy.kmc.atom_map",
    "rmgpy.kmc.barrier_e0",
    "rmgpy.kmc.compile",
    "rmgpy.kmc.compiler",
    "rmgpy.kmc.database_provenance",
    "rmgpy.kmc.event_record",
    "rmgpy.kmc.kinetics_library",
    "rmgpy.kmc.proxy_padding",
    "rmgpy.kmc.reference_thermo",
    "rmgpy.data.thermo",
    "rmgpy.data.solvation",
    "rmgpy.thermo.thermoengine",
)

_SIBLING_SUFFIXES = (".py", ".pyx", ".pxd")


def module_files(name: str) -> dict[str, str]:
    """``{file name: sha256}`` for a module's loaded file and its source siblings.

    Located from the import system (``find_spec``), never from the working
    directory.  Extension modules contribute their ``.so`` and any ``.py``,
    ``.pyx`` and ``.pxd`` with the same stem beside it.
    """
    try:
        spec = importlib.util.find_spec(name)
    except (ImportError, ValueError, AttributeError) as error:
        raise IdentityError(f"cannot locate module {name!r}: {error}") from error
    if spec is None or not spec.origin or spec.origin in ("built-in", "frozen"):
        raise IdentityError(f"module {name!r} has no file to hash")
    origin = Path(spec.origin).resolve()
    paths = {origin}
    stem = origin.name.split(".")[0]
    for suffix in _SIBLING_SUFFIXES:
        sibling = origin.with_name(stem + suffix)
        if sibling.is_file():
            paths.add(sibling.resolve())
    result = {}
    for path in sorted(paths):
        result[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def module_manifest(names: Iterable[str]) -> dict[str, dict[str, str]]:
    return {name: module_files(name) for name in sorted(set(names))}


def git_state(root: str | Path) -> dict[str, Any]:
    """``{"sha", "dirty"}`` of the git checkout rooted exactly at ``root``, else both None."""
    none = {"sha": None, "dirty": None}
    try:
        root = Path(root).resolve()
        top = Path(
            subprocess.check_output(
                ["git", "-C", str(root), "rev-parse", "--show-toplevel"],
                text=True, stderr=subprocess.DEVNULL,
            ).strip()
        ).resolve()
        if top != root:
            return none
        sha = subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
        status = subprocess.check_output(
            ["git", "-C", str(root), "status", "--porcelain", "--untracked-files=no"],
            text=True, stderr=subprocess.DEVNULL,
        )
        return {"sha": sha, "dirty": bool(status.strip())}
    except (OSError, subprocess.CalledProcessError):
        return none


def _package_root() -> Path:
    spec = importlib.util.find_spec("rmgpy")
    if spec is None or not spec.origin:
        raise IdentityError("cannot locate the rmgpy package")
    return Path(spec.origin).resolve().parent.parent


def _library_versions() -> dict[str, Any]:
    import numpy
    import scipy

    numpy_build: Any = None
    try:
        numpy_build = numpy.__config__.show(mode="dicts")  # numpy >= 1.25
    except Exception:  # pragma: no cover - older NumPy
        numpy_build = None
    try:
        from rdkit import rdBase, __version__ as rdkit_version

        rdkit = {
            "version": rdkit_version,
            "build": str(getattr(rdBase, "rdkitBuild", "")),
            "boost": str(getattr(rdBase, "boostVersion", "")),
        }
    except ImportError as error:
        raise IdentityError(f"RDKit is required for an identity: {error}") from error
    return {
        "numpy": {
            "version": numpy.__version__,
            "git_revision": getattr(numpy.version, "git_revision", None),
            "build": _jsonable(numpy_build),
        },
        "scipy": {
            "version": scipy.__version__,
            "git_revision": getattr(scipy.version, "git_revision", None),
        },
        "rdkit": rdkit,
    }


def _jsonable(value: Any) -> Any:
    try:
        return json.loads(json.dumps(value, default=str, sort_keys=True))
    except (TypeError, ValueError):
        return str(value)


@functools.lru_cache(maxsize=None)
def _cached_manifest(kind: str, seed: str) -> str:
    modules = GENERATION_MODULES + (COMPILATION_EXTRA_MODULES if kind == "compilation" else ())
    manifest = {
        "git": git_state(_package_root()),
        "hash_seed": seed,
        "kind": kind,
        "libraries": _library_versions(),
        "modules": module_manifest(modules),
        "platform": {"machine": platform.machine(), "system": platform.system()},
        "python": {
            "cache_tag": sys.implementation.cache_tag,
            "hexversion": sys.hexversion,
            "implementation": platform.python_implementation(),
            "version": sys.version,
        },
    }
    return json.dumps(manifest, sort_keys=True)


def generation_dependency_manifest() -> dict[str, Any]:
    """The fixed manifest a generated-reaction payload depends on."""
    return json.loads(_cached_manifest("generation", require_stable_hash_seed()))


def compilation_dependency_manifest() -> dict[str, Any]:
    """The fixed manifest an artifact depends on (generation modules plus the compiler)."""
    return json.loads(_cached_manifest("compilation", require_stable_hash_seed()))


def manifest_sha(manifest: dict[str, Any]) -> str:
    return sha256_canonical("dependency-manifest", manifest)


def generation_producer_identity() -> str:
    return manifest_sha(generation_dependency_manifest())


def compilation_producer_identity() -> str:
    return manifest_sha(compilation_dependency_manifest())


# ---------------------------------------------------------------------------
# Generic object canonicalisation (family content, constraints)
# ---------------------------------------------------------------------------

_ADDRESS = re.compile(r" at 0x[0-9a-fA-F]+")


def _object_state(value: Any) -> Any:
    """Deterministic JSON tree for a database/constraint object, or refuse."""
    if value is None or type(value) in (bool, int, float, str):
        return _canon(value)
    if isinstance(value, (list, tuple)):
        return [_object_state(item) for item in value]
    if isinstance(value, dict):
        return {"$dict": [[_object_state(k), _object_state(v)] for k, v in value.items()]}
    if isinstance(value, (set, frozenset)):
        return {"$set": sorted(json.dumps(_object_state(item), sort_keys=True) for item in value)}
    if hasattr(value, "item") and type(value).__module__ == "numpy":
        return _object_state(value.item())
    if hasattr(value, "reactants") and hasattr(value, "products"):
        return {
            "$reaction": {
                "family": getattr(value, "family", None),
                "kinetics": _object_state(getattr(value, "kinetics", None)),
                "label": getattr(value, "label", None),
                "products": [_object_state(item) for item in value.products],
                "reactants": [_object_state(item) for item in value.reactants],
                "reversible": getattr(value, "reversible", None),
                "template": _object_state(getattr(value, "template", None)),
            },
            "$type": type(value).__name__,
        }
    if hasattr(value, "molecule") and isinstance(getattr(value, "molecule"), list):
        return {
            "$species": [m.to_adjacency_list() for m in value.molecule],
            "label": getattr(value, "label", None),
            "$type": type(value).__name__,
        }
    if hasattr(value, "to_adjacency_list"):
        return {"$adjacency_list": value.to_adjacency_list(), "$type": type(value).__name__}
    text = str(value) if type(value).__name__.startswith("Logic") else repr(value)
    if _ADDRESS.search(text):
        raise IdentityError(f"{type(value).__name__} has no deterministic representation")
    return {"$repr": text, "$type": type(value).__name__}


def _entry_state(entry: Any) -> Any:
    state = {}
    for name, value in sorted(vars(entry).items()):
        if name == "parent":
            state[name] = getattr(value, "label", None)
        elif name == "children":
            state[name] = [getattr(child, "label", None) for child in (value or [])]
        else:
            state[name] = _object_state(value)
    return state


def _database_state(database: Any) -> Any:
    """Ordered entries of a Database-like object (order is generation-relevant)."""
    entries = []
    for label, entry in database.entries.items():
        if isinstance(entry, list):
            entries.append([label, [_entry_state(item) for item in entry]])
        else:
            entries.append([label, _entry_state(entry)])
    top = [getattr(entry, "label", None) for entry in (getattr(database, "top", None) or [])]
    return {"entries": entries, "top": top}


_FAMILY_SCALARS = (
    "label", "name", "solvent", "short_desc", "long_desc", "metal", "site", "facet", "reverse", "reversible",
    "own_reverse", "boundary_atoms", "tree_distances", "save_order", "reverse_map", "reactant_num",
    "product_num", "auto_generated", "allow_charged_species", "electrons",
)


def _template_state(template: Any) -> Any:
    if template is None:
        return None
    return {
        "reactants": [getattr(item, "label", None) for item in template.reactants],
        "products": [getattr(item, "label", None) for item in template.products],
    }


def _recipe_state(recipe: Any) -> Any:
    return None if recipe is None else _object_state(list(recipe.actions))


def family_content_fingerprint(family: Any) -> str:
    """SHA-256 over a family's recipe, groups/tree, rules, forbidden structures and flags."""
    state = {name: _object_state(getattr(family, name, None)) for name in _FAMILY_SCALARS}
    state["forward_template"] = _template_state(family.forward_template)
    state["reverse_template"] = _template_state(family.reverse_template)
    state["forward_recipe"] = _recipe_state(family.forward_recipe)
    state["reverse_recipe"] = _recipe_state(family.reverse_recipe)
    state["groups"] = _database_state(family.groups)
    state["rules"] = _database_state(family.rules)
    forbidden = getattr(family, "forbidden", None)
    state["forbidden"] = None if forbidden is None else _database_state(forbidden)
    return sha256_canonical("family-content", state)


# ---------------------------------------------------------------------------
# Concrete generation request
# ---------------------------------------------------------------------------


class _Numbering:
    """First-visit numbering so object sharing (aliasing) is part of the hash."""

    def __init__(self) -> None:
        self.seen: dict[int, int] = {}
        self.keep: list[Any] = []

    def number(self, obj: Any) -> int:
        key = id(obj)
        if key not in self.seen:
            self.seen[key] = len(self.seen)
            self.keep.append(obj)
        return self.seen[key]


def _atom_state(atom: Any, index: dict[int, int], numbering: _Numbering) -> list:
    return [
        numbering.number(atom),
        atom.element.symbol,
        atom.element.isotope,
        atom.radical_electrons,
        atom.charge,
        atom.label,
        atom.lone_pairs,
        atom.id,
        None if atom.atomtype is None else atom.atomtype.label,
        atom.site,
        atom.morphology,
        _object_state(dict(atom.props)),
        [
            [index[id(neighbour)], _object_state(bond.order), index[id(bond.vertex1)]]
            for neighbour, bond in atom.edges.items()
        ],
    ]


def _molecule_state(molecule: Any, numbering: _Numbering) -> dict[str, Any]:
    atoms = list(molecule.atoms)
    index = {id(atom): position for position, atom in enumerate(atoms)}
    for atom in atoms:
        for neighbour in atom.edges:
            if id(neighbour) not in index:
                raise IdentityError("molecule has a bond to an atom outside the molecule")
    return {
        "atoms": [_atom_state(atom, index, numbering) for atom in atoms],
        "facet": molecule.facet,
        "metal": molecule.metal,
        "multiplicity": molecule.multiplicity,
        "number": numbering.number(molecule),
        "polymer_proxy": molecule.is_polymer_proxy,
        "props": _object_state(dict(molecule.props)),
        "reactive": molecule.reactive,
        "symmetry_number": _object_state(float(molecule.symmetry_number)),
        "type": type(molecule).__name__,
    }


def participant_state(participant: Any, numbering: Optional[_Numbering] = None) -> dict[str, Any]:
    """Observable state of one generator participant (a Species or a Molecule)."""
    numbering = numbering or _Numbering()
    if hasattr(participant, "molecule") and not hasattr(participant, "atoms"):
        thermo = participant.thermo
        return {
            "aug_inchi": participant.aug_inchi,
            "creation_iteration": participant.creation_iteration,
            "explicitly_allowed": participant.explicitly_allowed,
            "has_conformer": participant.conformer is not None,
            "has_transport_data": participant.transport_data is not None,
            "index": participant.index,
            "is_moment_dummy": participant.is_moment_dummy,
            "is_solvent": participant.is_solvent,
            "kind": "Species",
            "label": participant.label,
            "molecules": [_molecule_state(item, numbering) for item in participant.molecule],
            "number": numbering.number(participant),
            "props": _object_state(dict(participant.props)),
            "reactive": participant.reactive,
            "symmetry_number": _object_state(float(participant.symmetry_number)),
            "thermo": None if thermo is None else _object_state(thermo),
            "type": type(participant).__name__,
        }
    return {"kind": "Molecule", "molecule": _molecule_state(participant, numbering)}


def generation_policy() -> dict[str, Any]:
    """Effective global state that can change generator output without changing arguments.

    Found by reading ``rmgpy/data/kinetics/family.py``, ``common.py`` and
    ``rmgpy/constraints.py``: the missing-reverse tolerance and its audit
    path/running counts (family.py), the process-wide species and polymer size
    constraints read through ``rmgpy.rmg.input.get_input`` (constraints.py),
    the polymer-size routing default, and the profiling switch (output
    neutral, hashed anyway).
    """
    from rmgpy import constraints
    from rmgpy.data.kinetics import family as family_module

    tolerance = family_module.TOLERATE_MISSING_REVERSE
    policy: dict[str, Any] = {
        "default_polymer_size_threshold": constraints.DEFAULT_POLYMER_SIZE_THRESHOLD,
        "dropped_reverse_log_path": family_module.DROPPED_REVERSE_LOG_PATH,
        "profile_wasted_builds": bool(family_module.PROFILE_WASTED_BUILDS),
        "tolerate_missing_reverse": {"type": type(tolerance).__name__, "value": _object_state(tolerance)},
    }
    policy["dropped_reverse_counts"] = (
        None
        if tolerance is False
        else _object_state(dict(sorted(family_module._dropped_reverse_counts.items())))
    )
    try:
        from rmgpy.rmg.input import get_input

        for name in ("species_constraints", "polymer_constraints"):
            try:
                policy[name] = {"configured": True, "value": _object_state(get_input(name))}
            except Exception:
                policy[name] = {"configured": False, "value": None}
    except ImportError as error:  # pragma: no cover
        raise IdentityError(f"cannot read the generation constraints: {error}") from error
    return policy


def _resolve_families(kinetics_database: Any, only_families: Optional[Sequence[str]]) -> tuple[list[str], list[str]]:
    available = list(kinetics_database.families)
    if only_families is None:
        return available, []
    requested = list(only_families)
    selected = [label for label in available if label in requested]
    unknown = sorted(label for label in set(requested) if label not in kinetics_database.families)
    return selected, unknown


def generation_request_description(
    kinetics_database: Any,
    reactants: Sequence[Any],
    products: Optional[Sequence[Any]] = None,
    only_families: Optional[Sequence[str]] = None,
    resonance: bool = True,
) -> dict[str, Any]:
    """The exact structure hashed by ``generation_request_identity``."""
    numbering = _Numbering()
    selected, unknown = _resolve_families(kinetics_database, only_families)
    return {
        "dependency_manifest_sha": manifest_sha(generation_dependency_manifest()),
        "families": {
            "fingerprints": {label: family_content_fingerprint(kinetics_database.families[label]) for label in selected},
            "selected": selected,
            "unknown": unknown,
        },
        "policy": generation_policy(),
        "products": (
            None if products is None else [participant_state(item, numbering) for item in products]
        ),
        "reactants": [participant_state(item, numbering) for item in reactants],
        "resonance": bool(resonance) if type(resonance) is bool else _object_state(resonance),
    }


def generation_request_identity(
    kinetics_database: Any,
    reactants: Sequence[Any],
    products: Optional[Sequence[Any]] = None,
    only_families: Optional[Sequence[str]] = None,
    resonance: bool = True,
) -> str:
    """SHA-256 of one concrete call to ``generate_reactions_from_families``.

    Call this BEFORE the generator runs: the generator mutates its inputs
    (atom ids, labels), and the identity is of the state the call observes.
    """
    require_stable_hash_seed()
    description = generation_request_description(
        kinetics_database, reactants, products, only_families, resonance
    )
    return sha256_canonical("generation-request", description)


# ---------------------------------------------------------------------------
# Compilation request
# ---------------------------------------------------------------------------


def compile_environment_switches() -> dict[str, str]:
    """Effective environment switches that change a compile (resolved, validated).

    ``RMG_KMC_PLPSEC_LIBRARY`` (compiler.py ``EventSetCompiler.__init__``) is the
    only one the compile path reads; ``RMG_DATABASE_SHA`` is resolved into the
    database declaration instead.  ``RMG_I067_PROFILE`` is output neutral and
    is hashed through ``generation_policy``.
    """
    selection = os.environ.get("RMG_KMC_PLPSEC_LIBRARY", "1")
    if selection not in {"0", "1"}:
        raise ValueError("RMG_KMC_PLPSEC_LIBRARY must be 0 or 1")
    return {"RMG_KMC_PLPSEC_LIBRARY": selection}


def compilation_request_description(request: Any) -> dict[str, Any]:
    """The exact structure hashed by ``compilation_request_identity``."""
    from rmgpy.kmc import compiler
    from rmgpy.kmc.compile import enumerate_family_universe
    from rmgpy.kmc.database_provenance import database_provenance, resolve_database_declaration

    database_path = Path(request.database_path)
    declared = resolve_database_declaration(request.database_sha)
    family_universe = (
        list(request.family_universe)
        if request.family_universe is not None
        else enumerate_family_universe(database_path)
    )
    rmgpy_path = Path(compiler.__file__).resolve().parents[2]
    return {
        "database": database_provenance(database_path, declared),
        "dependency_manifest_sha": manifest_sha(compilation_dependency_manifest()),
        "environment": compile_environment_switches(),
        "producer": {
            "compiler_sha256": compiler._LOADED_COMPILER_HASH,
            "compiler_sources_sha256": compiler.compiler_source_hash(),
            "rmgpy_sha": compiler._git_sha(rmgpy_path),
        },
        "request": {
            "barrier_e0_fixed_b": request.barrier_e0_fixed_b,
            "database_sha": declared,
            "families": list(request.families),
            "family_universe": family_universe,
            "proxy_units": (
                compiler.PS_PROXY_UNITS if request.proxy_units is None else request.proxy_units
            ),
            "temperature_grid": [float(item) for item in request.temperature_grid],
        },
    }


def compilation_request_identity(request: Any) -> str:
    """SHA-256 of a normalised ``CompileRequest`` (no ``output_dir``, no database path)."""
    require_stable_hash_seed()
    return sha256_canonical("compilation-request", compilation_request_description(request))


def published_content_sha(artifact: dict[str, Any]) -> str:
    """The artifact content hash: SHA-256 of its canonical JSON, i.e. the artifact id."""
    from rmgpy.kmc.compiler import sha256_json

    return sha256_json(artifact)


# ---------------------------------------------------------------------------
# Snapshot binding
# ---------------------------------------------------------------------------


def _snapshot_files(database_path: str | Path) -> dict[str, str]:
    from rmgpy.kmc.database_provenance import _input_roots

    root = Path(database_path)
    if not root.is_dir():
        raise IdentityError(f"database path is not a directory: {database_path}")
    root = root.resolve()
    roots = _input_roots(root)
    if not roots:
        raise IdentityError(f"database has no supported kinetics/thermo input layout: {database_path}")
    files: dict[str, str] = {}

    def raise_walk_error(error):
        raise IdentityError(f"cannot walk database: {error}") from error

    for input_root in roots:
        if input_root.is_symlink():
            raise IdentityError(f"symlink under database input roots: {input_root}")
        for current, directories, filenames in os.walk(input_root, followlinks=False, onerror=raise_walk_error):
            current_path = Path(current)
            for name in [*directories, *filenames]:
                if (current_path / name).is_symlink():
                    raise IdentityError(f"symlink under database input roots: {current_path / name}")
            for name in filenames:
                path = current_path / name
                if "__pycache__" in path.parts:
                    continue
                files[path.relative_to(root).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    if not files:
        raise IdentityError(f"database input trees are empty: {database_path}")
    return dict(sorted(files.items()))


def database_snapshot_manifest(database_path: str | Path) -> dict[str, Any]:
    """Relative path to sha256 for every kinetics/thermo input file, plus a manifest hash.

    These are exactly the files ``database_content_digest`` reads, i.e. what
    generation and compilation consume.  Files outside the two input trees
    are neither read by the compiler nor listed.
    """
    files = _snapshot_files(database_path)
    return {
        "files": files,
        "manifest_sha256": sha256_canonical("database-snapshot", files),
        "schema_version": SCHEMA_VERSION,
    }


def verify_database_snapshot(manifest: dict[str, Any], database_path: str | Path) -> None:
    """Re-hash ``database_path`` and raise ``SnapshotMismatchError`` on any difference.

    Detects an edited, replaced (including by a symlink), added or removed
    file between the manifest and this call.
    """
    if (
        not isinstance(manifest, dict)
        or set(manifest) != {"files", "manifest_sha256", "schema_version"}
        or manifest["schema_version"] != SCHEMA_VERSION
        or not isinstance(manifest["files"], dict)
    ):
        raise SnapshotMismatchError("malformed snapshot manifest")
    if manifest["manifest_sha256"] != sha256_canonical("database-snapshot", manifest["files"]):
        raise SnapshotMismatchError("snapshot manifest hash does not match its file list")
    try:
        current = _snapshot_files(database_path)
    except IdentityError as error:
        raise SnapshotMismatchError(f"database no longer matches its snapshot: {error}") from error
    expected = manifest["files"]
    added = sorted(set(current) - set(expected))
    removed = sorted(set(expected) - set(current))
    changed = sorted(name for name in set(current) & set(expected) if current[name] != expected[name])
    if added or removed or changed:
        raise SnapshotMismatchError(
            f"database differs from its snapshot: added={added[:5]} removed={removed[:5]} changed={changed[:5]}"
        )


# ---------------------------------------------------------------------------
# Request-binding envelope
# ---------------------------------------------------------------------------

_ENVELOPE_KEYS = {"schema_version", "request_sha", "producer_identity", "payload_sha256", "payload"}
_HEX64 = re.compile(r"[0-9a-f]{64}")
DEFAULT_MAX_ENVELOPE_BYTES = 64 << 20
_ENVELOPE_MAX_DEPTH = 2
_ENVELOPE_MAX_KEYS = 16
_ENVELOPE_TOKEN = re.compile(rb'(?s)"[^"\\]*(?:\\.[^"\\]*)*"|[\[\]{}:,]|"')


def _require_hex(value: Any, name: str) -> str:
    if type(value) is not str or not _HEX64.fullmatch(value):
        raise EnvelopeError(f"{name} must be a 64-character lowercase hex SHA-256")
    return value


def build_envelope(request_sha: str, producer_identity: str, payload: bytes) -> dict[str, Any]:
    """Bind ``payload`` to ``request_sha`` and the producing manifest."""
    _require_hex(request_sha, "request_sha")
    _require_hex(producer_identity, "producer_identity")
    if not isinstance(payload, (bytes, bytearray)):
        raise EnvelopeError("payload must be bytes")
    try:
        text = bytes(payload).decode("utf-8")
    except UnicodeDecodeError as error:
        raise EnvelopeError("payload must be UTF-8 text") from error
    return {
        "payload": text,
        "payload_sha256": hashlib.sha256(bytes(payload)).hexdigest(),
        "producer_identity": producer_identity,
        "request_sha": request_sha,
        "schema_version": SCHEMA_VERSION,
    }


def verify_envelope(
    envelope: Any, *, request_sha: str, producer_identity: Optional[str] = None
) -> bytes:
    """Return the payload bytes iff the envelope is intact and bound to ``request_sha``."""
    if type(envelope) is not dict or set(envelope) != _ENVELOPE_KEYS:
        raise EnvelopeError("envelope must be an object with exactly the envelope keys")
    if envelope["schema_version"] != SCHEMA_VERSION or type(envelope["schema_version"]) is not int:
        raise EnvelopeError("unsupported envelope schema version")
    _require_hex(request_sha, "expected request_sha")
    _require_hex(envelope["request_sha"], "request_sha")
    _require_hex(envelope["producer_identity"], "producer_identity")
    _require_hex(envelope["payload_sha256"], "payload_sha256")
    if type(envelope["payload"]) is not str:
        raise EnvelopeError("payload must be a string")
    if envelope["request_sha"] != request_sha:
        raise EnvelopeError("envelope is bound to a different request")
    if producer_identity is not None:
        _require_hex(producer_identity, "expected producer_identity")
        if envelope["producer_identity"] != producer_identity:
            raise EnvelopeError("envelope was produced by a different producer")
    try:
        payload = envelope["payload"].encode("utf-8")
    except UnicodeEncodeError as error:
        raise EnvelopeError("payload is not valid UTF-8 text") from error
    if hashlib.sha256(payload).hexdigest() != envelope["payload_sha256"]:
        raise EnvelopeError("payload does not match its digest")
    return payload


def envelope_to_bytes(envelope: dict[str, Any]) -> bytes:
    return json.dumps(envelope, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode("ascii")


def _no_duplicates(pairs: list) -> dict:
    result: dict = {}
    for key, value in pairs:
        if key in result:
            raise EnvelopeError(f"duplicate key {key!r}")
        result[key] = value
    return result


def _preflight_envelope(data: bytes) -> None:
    """Bound nesting depth and key count by a token scan, before ``json.loads`` allocates anything."""
    depth = 0
    keys = 0
    for match in _ENVELOPE_TOKEN.finditer(data):
        token = match.group()
        first = token[:1]
        if first == b'"':
            if len(token) < 2:
                raise EnvelopeError("unterminated string")
        elif first == b"[" or first == b"{":
            depth += 1
            if depth > _ENVELOPE_MAX_DEPTH:
                raise EnvelopeError("envelope nesting depth exceeds the limit")
        elif first == b"]" or first == b"}":
            depth -= 1
            if depth < 0:
                raise EnvelopeError("unbalanced brackets")
        elif first == b":":
            keys += 1
            if keys > _ENVELOPE_MAX_KEYS:
                raise EnvelopeError("envelope key count exceeds the limit")
    if depth != 0:
        raise EnvelopeError("unbalanced brackets")


def envelope_from_bytes(data: bytes, *, max_bytes: int = DEFAULT_MAX_ENVELOPE_BYTES) -> dict[str, Any]:
    """Parse envelope bytes strictly (duplicate keys, non-finite numbers and non-objects refused).

    The size, nesting depth and key count are bounded before ``json.loads`` runs.
    """
    if not isinstance(data, (bytes, bytearray)) or len(data) > max_bytes:
        raise EnvelopeError("envelope bytes missing or too large")
    data = bytes(data)
    _preflight_envelope(data)

    def reject_constant(token):
        raise EnvelopeError(f"non-finite constant {token}")

    try:
        value = json.loads(data.decode("utf-8"), object_pairs_hook=_no_duplicates, parse_constant=reject_constant)
    except EnvelopeError:
        raise
    except (ValueError, RecursionError) as error:
        raise EnvelopeError(f"cannot parse envelope: {error}") from error
    if type(value) is not dict:
        raise EnvelopeError("envelope must be a JSON object")
    return value
