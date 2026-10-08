"""Compile the real PS event set in an isolated hash-seeded interpreter."""

import argparse
import hashlib
import inspect
import json
import logging
import os
import pickle
import subprocess
import tempfile
from pathlib import Path

from portable_cache import (
    atomic_copy,
    atomic_write,
    artifact_cache_key,
    compile_environment_options,
    identity,
    identity_name,
    migrate,
)

from rmgpy.data.rmg import RMGDatabase
from rmgpy.kmc.compiler import (
    EventSetCompiler,
    PS_FAMILY_CANDIDATES,
    PS_PROXY_UNITS,
    prepare_rate_rules,
    ps_proxy_set,
)
from rmgpy.kmc.database_provenance import (
    database_content_digest,
    resolve_database_declaration,
)


def dump_generated_reactions(reactions):
    pending = list(reactions)
    visited = set()
    atoms = {}
    while pending:
        reaction = pending.pop()
        if id(reaction) in visited:
            continue
        visited.add(id(reaction))
        for species in reaction.reactants + reaction.products:
            for molecule in species.molecule:
                for atom in molecule.atoms:
                    atoms[id(atom)] = atom, atom.id
        reverse = getattr(reaction, "reverse", None)
        if reverse is not None:
            pending.append(reverse)
    return pickle.dumps((reactions, list(atoms.values())), protocol=4)


def load_generated_reactions(payload):
    reactions, atom_ids = pickle.loads(payload)
    for atom, identifier in atom_ids:
        atom.id = identifier
    return reactions


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("database")
    parser.add_argument("output")
    parser.add_argument("--database-sha", help="commit of a materialized pinned snapshot")
    parser.add_argument("--family-universe", type=Path, help="JSON list from the pinned git tree")
    args = parser.parse_args()
    logger = logging.getLogger("rmgpy.kmc.compiler")
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)
    logger.propagate = False
    database_path = Path(args.database)
    family_root = database_path / "input/kinetics/families"
    family_universe = sorted(
        path.name for path in family_root.iterdir() if (path / "groups.py").is_file()
    )
    if args.family_universe:
        family_universe = json.loads(args.family_universe.read_text())

    database = RMGDatabase()
    print("loading pinned RMG families", flush=True)
    database.load_kinetics(
        str(database_path / "input/kinetics"),
        reaction_libraries=[],
        seed_mechanisms=None,
        kinetics_families=list(PS_FAMILY_CANDIDATES),
        kinetics_depositories=["training"],
    )
    database.load_thermo(
        str(database_path / "input/thermo"),
        thermo_libraries=["primaryThermoLibrary"],
        depository=True,
    )
    print("preparing non-auto-generated rate rules from training", flush=True)
    prepare_rate_rules(database.kinetics, database.thermo, verbose=True)
    database_commit = resolve_database_declaration(args.database_sha)
    database_content_sha, _ = database_content_digest(database_path)
    cache_root = Path(os.environ.get("RMG_KMC_CACHE_ROOT", str(Path.cwd() / ".kmc-cache")))
    cache_identity = identity(Path.cwd(), database_path)
    portable_root = cache_root / "portable" / identity_name(cache_identity)
    if not (portable_root / "manifest.json").exists():
        migrate(cache_root, Path.cwd(), database_path, database_commit)
    generated_cache = portable_root / "generated-reactions"
    generated_cache.mkdir(parents=True, exist_ok=True)
    compile_options = {
        "database_sha": database_commit,
        "database_content_sha256": database_content_sha,
        "family_universe": family_universe,
        "temperature_grid": list(__import__("rmgpy.kmc.compiler", fromlist=["DEFAULT_T_GRID"]).DEFAULT_T_GRID),
        "environment": compile_environment_options(),
    }
    artifact_key = artifact_cache_key(Path.cwd(), database_path, compile_options)
    artifact_cache = cache_root / "portable-artifacts" / artifact_key
    artifact_manifest = artifact_cache / "manifest.json"
    if artifact_manifest.is_file():
        artifact_name = json.loads(artifact_manifest.read_text())["artifact"]
        cached_artifact = artifact_cache / artifact_name
        if not cached_artifact.is_file():
            raise RuntimeError(f"artifact manifest names missing file: {cached_artifact}")
        output = Path(args.output)
        output.mkdir(parents=True, exist_ok=True)
        destination = output / cached_artifact.name
        atomic_copy(cached_artifact, destination)
        print(destination)
        return
    generate = database.kinetics.generate_reactions_from_families
    generation_source = hashlib.sha256(inspect.getsource(generate).encode()).hexdigest()

    def cached_generate(reactants, products=None, only_families=None, resonance=True):
        parameters = {
            "cache_schema": 1,
            "generator_source": generation_source,
            "reactants": [
                [
                    molecule.to_adjacency_list(remove_h=False)
                    for molecule in species.molecule
                ]
                for species in reactants
            ],
            "products": (
                None if products is None else [str(species) for species in products]
            ),
            "families": only_families,
            "resonance": resonance,
        }
        key = hashlib.sha256(
            json.dumps(parameters, sort_keys=True).encode()
        ).hexdigest()
        path = generated_cache / (key + ".pickle")
        if path.is_file():
            print(f"cached public RMG generation: {only_families}", flush=True)
            return load_generated_reactions(path.read_bytes())
        reactions = generate(reactants, products, only_families, resonance)
        temporary = path.with_suffix(f".{os.getpid()}.tmp")
        temporary.write_bytes(dump_generated_reactions(reactions))
        temporary.replace(path)
        return reactions

    database.kinetics.generate_reactions_from_families = cached_generate
    print("discovering PS chemistry", flush=True)
    proxies = ps_proxy_set(PS_PROXY_UNITS)
    active, excluded, reactions = EventSetCompiler.discover_family_reactions(
        database.kinetics,
        proxies,
        PS_FAMILY_CANDIDATES,
        family_universe=family_universe,
    )
    compiler = EventSetCompiler(
        database.kinetics,
        proxies,
        active,
        excluded_families=excluded,
        database_path=database_path,
        thermo_database=database.thermo,
        reaction_cache=reactions,
        rmg_database_sha=database_commit,
    )
    print(
        f"compiling {sum(len(value) for value in reactions.values())} generated reactions",
        flush=True,
    )
    artifact_cache.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=artifact_cache.parent) as temporary:
        path, _ = compiler.write_artifact(temporary)
        final_path = artifact_cache / path.name
        os.replace(path, final_path)
    path = final_path
    atomic_write(artifact_cache / "manifest.json", json.dumps({
        "identity": artifact_cache_key(Path.cwd(), database_path, compile_options),
        "artifact": path.name,
    }, sort_keys=True, indent=2).encode() + b"\n")
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    destination = output / path.name
    if destination != path:
        atomic_copy(path, destination)
    path = destination
    print(path)


if __name__ == "__main__":
    main()
