"""The one cold-compile procedure for a kMC event-set artifact.

``compile_artifact`` loads the database, prepares rate rules, discovers family
reactions, compiles, and writes the content-addressed artifact.  Callers (this
module's CLI, the test fixture driver) share it, so there is exactly one
definition of what a cold compile does.

Nothing here is specific to one polymer: the family candidate list, the family
universe and the proxy size are inputs.  With a deck (``--deck``) the proxies
are the declared polymer's own (``declared_proxies.declared_proxy_set``) and
every requested family is a candidate on them; without one the built-in proxy
set is compiled as before.  A schema-v2 (copolymer) deck's proxies
(``declared_proxies.sequence_proxy_set``) pair its structures over the
family-tree classes of the loaded families, so they are built once the
families load, and ``proxy_source`` records the pairs.  The artifact's provenance names the source
(``proxy_source``: the deck's declaration digests, or the built-in set).  Caching is not part of this module;
the only seam is ``generation_provider``, a callable with the signature of
``RMGDatabase.kinetics.generate_reactions_from_families`` that every reaction
generation of the compile goes through (default: the database's own method).

CLI::

    python -m rmgpy.kmc.compile DATABASE OUTPUT (--families FILE.json | --family F)
        [--family-universe FILE.json] [--database-sha SHA]
        [--barrier-e0-fixed-b K] [--proxy-units N] [--deck DECK.py]

``--family F`` (exactly one family) compiles that family's artifact alone: the
one-element form of ``--families``.  ``F`` must be one family directory of the
database (``input/kinetics/families/F/groups.py``); a selector that RMG would
expand (``all``, ``none``, a recommended set such as ``default``, ``!F``) or
any other name is refused before anything loads.

The CLI runs the compile in a spawned child process that owns its own
``RMGDatabase``; the child exits promptly if the parent dies.
"""

from __future__ import annotations

import argparse
import json
import logging
import multiprocessing
import os
import sys
import threading
import time
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional, Sequence

from rmgpy.data.rmg import RMGDatabase
from rmgpy.kmc.barrier_e0 import FixedBBarrierE0Provider
from rmgpy.kmc.compiler import (
    DEFAULT_T_GRID,
    REQUESTED_FAMILY_FILTER_REASON,
    EventSetCompiler,
    prepare_rate_rules,
    ps_proxy_set,
)
from rmgpy.kmc.database_provenance import resolve_database_declaration
from rmgpy.molecule.molecule import Atom, Molecule

# Seconds between parent-liveness polls.
PARENT_POLL_SECONDS = 0.5
# RMG's atom-id counter (rmgpy.molecule.molecule): a fresh process's first id,
# and the number of ids before it wraps
FRESH_ATOM_ID = -(2 ** 15)
ATOM_ID_SPAN = 2 ** 16


@dataclass(frozen=True)
class CompileRequest:
    """Plain-value description of one cold compile."""

    database_path: str
    output_dir: str
    families: tuple[str, ...]
    # None enumerates every family directory of the database that has groups.py.
    family_universe: Optional[tuple[str, ...]] = None
    database_sha: Optional[str] = None
    barrier_e0_fixed_b: Optional[float] = None
    # None uses the compiler's own proxy-set default.
    proxy_units: Optional[int] = None
    temperature_grid: tuple[float, ...] = tuple(DEFAULT_T_GRID)
    # a deck file declaring the polymer; None compiles the built-in proxy set
    deck: Optional[str] = None


def enumerate_family_universe(database_path: str | Path) -> list[str]:
    """Every family of a database tree, from the directory layout alone."""
    family_root = Path(database_path) / "input/kinetics/families"
    return sorted(
        path.name for path in family_root.iterdir() if (path / "groups.py").is_file()
    )


def compile_artifact(
    request: CompileRequest,
    *,
    generation_provider_factory: Optional[Callable[[Callable], Callable]] = None,
) -> tuple[Path, str]:
    """Cold-compile ``request`` and return ``(artifact path, artifact id)``.

    The database is created inside this call, so a caller cannot hold its
    generator beforehand.  ``generation_provider_factory(default_generate)``
    is called once with the database's own ``generate_reactions_from_families``
    and returns the provider the compiler uses instead: a callable with that
    method's exact signature ``(reactants, products=None, only_families=None,
    resonance=True)``.  Every reaction generation of the compile (discovery,
    per-proxy generation, padded-root regeneration) goes through it.  Without
    a factory the database's own method is used and nothing is cached.
    """
    inputs = load_compile_inputs(request)
    database = inputs.database
    generation_provider = (
        generation_provider_factory(database.kinetics.generate_reactions_from_families)
        if generation_provider_factory is not None
        else None
    )
    print("discovering chemistry", flush=True)
    active, excluded, reactions = EventSetCompiler.discover_family_reactions(
        database.kinetics,
        inputs.proxies,
        request.families,
        family_universe=inputs.family_universe,
        generation_provider=generation_provider,
        **inputs.discovery_options,
    )
    compiler = make_compiler(
        request, inputs, active, excluded, reactions, generation_provider
    )
    print(
        f"compiling {sum(len(value) for value in reactions.values())} generated reactions",
        flush=True,
    )
    path, _ = compiler.write_artifact(request.output_dir)
    return Path(path), Path(path).stem


@dataclass
class CompileInputs:
    """What a compile reads before it discovers anything: the loaded database,
    the proxies and their source, and the options the request implies."""

    database: RMGDatabase
    proxies: Sequence
    proxy_source: dict
    compiler_options: dict
    discovery_options: dict
    family_universe: list
    barrier_e0_provider: Optional[FixedBBarrierE0Provider]
    database_commit: Optional[str]


def next_atom_id() -> int:
    """The atom id RMG assigns next (reading it assigns it)."""
    probe = Molecule(atoms=[Atom(element="C")])
    probe.assign_atom_ids()
    return probe.atoms[0].id


def advance_atom_ids(count: int) -> None:
    """Assign ``count`` atom ids (modulo the counter's span) to nothing."""
    count %= ATOM_ID_SPAN
    if count:
        Molecule(atoms=[Atom(element="C") for _ in range(count)]).assign_atom_ids()


def seek_atom_ids(value: int) -> None:
    """Make ``value`` the atom id RMG assigns next."""
    advance_atom_ids(value - next_atom_id() - 1)


def peek_atom_id() -> int:
    """The atom id RMG assigns next, leaving it unassigned."""
    value = next_atom_id()
    seek_atom_ids(value)
    return value


def replay_atom_ids(database, step: str, start: int) -> None:
    """Leave the atom-id counter where ``step``, begun at ``start``, left it the
    first time it ran on ``database``, whatever it did this time (a reused
    database skips its load; a second rate-rule preparation may assign fewer)."""
    recorded = getattr(database, "kmc_atom_ids_assigned", None)
    if recorded is None:
        recorded = {}
        setattr(database, "kmc_atom_ids_assigned", recorded)
    if step not in recorded:
        recorded[step] = (peek_atom_id() - start) % ATOM_ID_SPAN
    seek_atom_ids(start)
    advance_atom_ids(recorded[step])


def load_compile_inputs(request: CompileRequest, database=None) -> CompileInputs:
    """Read the declaration, load the pinned database (or reuse ``database``,
    one this function loaded for the same request) and prepare its rate rules.

    RMG numbers the atoms of every generated reaction from one process-wide
    counter (16 bits, wrapping), and the compiled records order their bond
    rewrites by those ids; the counter is first set to where a fresh process
    starts it, so a compile does not depend on what ran before it."""
    seek_atom_ids(FRESH_ATOM_ID)
    database_path = Path(request.database_path)
    # the declaration is read before anything loads
    proxies, proxy_source, compiler_options, discovery_options = _proxy_inputs(request)
    family_universe = (
        list(request.family_universe)
        if request.family_universe is not None
        else enumerate_family_universe(database_path)
    )
    barrier_e0_provider = (
        FixedBBarrierE0Provider(request.barrier_e0_fixed_b)
        if request.barrier_e0_fixed_b is not None
        else None
    )
    start = peek_atom_id()
    if database is None:
        database = RMGDatabase()
        print("loading pinned RMG families", flush=True)
        database.load_kinetics(
            str(database_path / "input/kinetics"),
            reaction_libraries=[],
            seed_mechanisms=None,
            kinetics_families=list(request.families),
            kinetics_depositories=["training"],
        )
        database.load_thermo(
            str(database_path / "input/thermo"),
            thermo_libraries=["primaryThermoLibrary"],
            depository=True,
        )
    replay_atom_ids(database, "load", start)
    if callable(proxies):
        # a schema-v2 deck: the pairs are read from the loaded family trees
        print("pairing structures over family-tree classes", flush=True)
        proxies = proxies(database.kinetics.families)
    print("preparing non-auto-generated rate rules from training", flush=True)
    start = peek_atom_id()
    prepare_rate_rules(database.kinetics, database.thermo)
    replay_atom_ids(database, "prepare", start)
    database_commit = resolve_database_declaration(request.database_sha)
    return CompileInputs(
        database=database,
        proxies=proxies,
        proxy_source=proxy_source,
        compiler_options=compiler_options,
        discovery_options=discovery_options,
        family_universe=family_universe,
        barrier_e0_provider=barrier_e0_provider,
        database_commit=database_commit,
    )


def make_compiler(
    request: CompileRequest,
    inputs: CompileInputs,
    active,
    excluded,
    reactions,
    generation_provider=None,
) -> EventSetCompiler:
    """The compiler of ``request`` over the discovered ``active`` families."""
    return EventSetCompiler(
        inputs.database.kinetics,
        inputs.proxies,
        active,
        excluded_families=excluded,
        database_path=Path(request.database_path),
        thermo_database=inputs.database.thermo,
        reaction_cache=reactions,
        rmg_database_sha=inputs.database_commit,
        barrier_e0_provider=inputs.barrier_e0_provider,
        temperature_grid=request.temperature_grid,
        generation_provider=generation_provider,
        proxy_source=inputs.proxy_source,
        **inputs.compiler_options,
    )


def _proxy_inputs(request: CompileRequest):
    """(proxies, proxy source, compiler options, discovery options) of a
    request: the declared polymer's proxies and its requested families as the
    candidates, or the built-in set and its own candidates."""
    if request.deck is None:
        proxies = (
            ps_proxy_set(request.proxy_units)
            if request.proxy_units is not None
            else ps_proxy_set()
        )
        source = {
            "kind": "built-in",
            "proxy_set": f"{ps_proxy_set.__module__}.{ps_proxy_set.__qualname__}",
            "proxy_units": request.proxy_units,
            "reason": "no deck was given",
        }
        return proxies, source, {}, {}
    from rmgpy.kmc.declared_proxies import (
        DEFAULT_PROXY_UNITS,
        declared_proxy_set,
        sequence_proxy_set,
        sequence_structures,
    )
    from rmgpy.kmc.run import read_deck
    from rmgpy.polymer_spec import PolymerSpec

    deck = read_deck(request.deck)
    units = request.proxy_units if request.proxy_units is not None else DEFAULT_PROXY_UNITS
    build = sequence_proxy_set if isinstance(deck.spec, PolymerSpec) else declared_proxy_set
    source = {
        "kind": "declaration",
        "label": deck.label,
        "deck_sha256": deck.sha256,
        **deck.digests(),
        "proxy_units": units,
        "proxy_set": f"{build.__module__}.{build.__qualname__}",
    }
    if build is sequence_proxy_set:
        # the structures are built (and refused) before anything loads
        sequence_structures(deck.spec, units)

        def proxies(families):
            built, record = sequence_proxy_set(deck.spec, units, families)
            source.update(pairs=record["pairs"], pair_rule=record["rule"],
                          pair_families=record["families"], pair_classes=record["classes"])
            return built
    else:
        proxies = declared_proxy_set(deck.spec, units)
    compiler_options = {
        "family_candidates": tuple(request.families),
        "family_filter_reason": REQUESTED_FAMILY_FILTER_REASON,
        # the built-in set's archived ring-junction oracle is not this polymer's
        "archived_junctions": False,
    }
    return proxies, source, compiler_options, {"proxy_set_name": "declared proxy set"}


def _watch_parent(parent_pid: int) -> None:
    while os.getppid() == parent_pid:
        time.sleep(PARENT_POLL_SECONDS)
    os._exit(1)


def _child_main(request: CompileRequest, parent_pid: int, connection) -> None:
    """Compile in the child; report ``("ok", ...)`` or ``("error", traceback)``."""
    threading.Thread(
        target=_watch_parent, args=(parent_pid,), daemon=True, name="parent-watch"
    ).start()
    try:
        _configure_logging()
        path, artifact_id = compile_artifact(request)
        connection.send(("ok", (str(path), artifact_id)))
    except BaseException:  # the parent must always learn why the child failed
        connection.send(("error", traceback.format_exc()))
        sys.stdout.flush()
        sys.stderr.flush()
        os._exit(1)
    finally:
        connection.close()


def _configure_logging() -> None:
    logger = logging.getLogger("rmgpy.kmc.compiler")
    if not logger.handlers:
        logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)
    logger.propagate = False


def compile_in_child(request: CompileRequest) -> tuple[Path, str]:
    """Run ``compile_artifact`` in a spawned child with its own database.

    The child's stdout and stderr are the parent's.  A child exception is
    re-raised here as ``RuntimeError`` carrying the child's traceback.
    """
    context = multiprocessing.get_context("spawn")
    receiver, sender = context.Pipe(duplex=False)
    child = context.Process(
        target=_child_main, args=(request, os.getpid(), sender), name="kmc-compile"
    )
    child.start()
    sender.close()
    message = None
    try:
        while message is None:
            if receiver.poll(PARENT_POLL_SECONDS):
                try:
                    message = receiver.recv()
                except EOFError:
                    break
            elif not child.is_alive():
                break
    finally:
        child.join()
        receiver.close()
    if message is None:
        raise RuntimeError(
            f"compile child exited with code {child.exitcode} without a result"
        )
    status, payload = message
    if status == "error":
        raise RuntimeError("compile child failed:\n" + payload)
    path, artifact_id = payload
    return Path(path), artifact_id


def one_family_problem(database_path, name: str) -> Optional[str]:
    """Why ``name`` is not exactly one family of the database (None if it is):
    it must be the name of a family directory with ``groups.py``, never a
    selector RMG expands into several families or none."""
    if name.startswith("!") or name.lower() in ("all", "none"):
        return (f"--family {name!r} is a family selector, not one family: name exactly one family "
                "directory of the database")
    try:
        universe = enumerate_family_universe(database_path)
    except OSError as error:
        return f"--family {name!r}: the database has no readable input/kinetics/families ({error})"
    if name not in universe:
        return (f"--family {name!r} is not a family of {database_path} (no "
                f"input/kinetics/families/{name}/groups.py; recommended sets such as 'default' "
                "expand into several families and are refused)")
    return None


def _read_list(path: str) -> tuple[str, ...]:
    value = json.loads(Path(path).read_text())
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ValueError(f"{path} must hold a JSON list of strings")
    return tuple(value)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m rmgpy.kmc.compile")
    parser.add_argument("database")
    parser.add_argument("output")
    add_request_arguments(parser)
    args = parser.parse_args(argv)
    try:
        request = request_from_arguments(parser, args, args.output)
        path, artifact_id = compile_in_child(request)
    except Exception:
        traceback.print_exc()
        return 1
    print(path)
    print(artifact_id)
    return 0


def add_request_arguments(parser: argparse.ArgumentParser) -> None:
    """The arguments naming a compile request (all but the database and output)."""
    which = parser.add_mutually_exclusive_group(required=True)
    which.add_argument("--families", help="JSON list of family candidates")
    which.add_argument("--family", action="append",
                       help="exactly one family: the artifact of that family alone")
    parser.add_argument("--family-universe", help="JSON list from the pinned git tree")
    parser.add_argument("--database-sha", help="commit of a materialized pinned snapshot")
    parser.add_argument(
        "--barrier-e0-fixed-b",
        type=float,
        help="opt in to fixed-B Wilhoit E0 values for compiler barrier floors only",
    )
    parser.add_argument("--proxy-units", type=int, help="proxy chain length")
    parser.add_argument("--deck", help="RMG deck declaring the polymer: compile its own "
                        "proxies with every requested family as a candidate")


def request_from_arguments(
    parser: argparse.ArgumentParser, args: argparse.Namespace, output: str
) -> CompileRequest:
    """The request ``args`` name (``parser`` reports a bad ``--family``)."""
    if args.family is not None:
        if len(args.family) != 1 or not args.family[0]:
            parser.error(f"--family takes exactly one family name, got {args.family}")
        problem = one_family_problem(args.database, args.family[0])
        if problem:
            parser.error(problem)
    return CompileRequest(
        database_path=str(Path(args.database).resolve()),
        output_dir=str(Path(output).resolve()),
        families=(
            tuple(args.family) if args.family is not None else _read_list(args.families)
        ),
        family_universe=(
            _read_list(args.family_universe) if args.family_universe else None
        ),
        database_sha=args.database_sha,
        barrier_e0_fixed_b=args.barrier_e0_fixed_b,
        proxy_units=args.proxy_units,
        deck=str(Path(args.deck).resolve()) if args.deck else None,
    )


if __name__ == "__main__":
    sys.exit(main())
