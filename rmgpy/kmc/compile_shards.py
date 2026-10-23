"""The cold compile (``rmgpy.kmc.compile``) split into deterministic shards.

A single-process compile discovers every (proxy, family) *unit* with RMG's own
``generate_reactions_from_families`` and then compiles every candidate pair in
one loop.  On large proxies most of the time goes to one step inside RMG's
generation: ``find_degenerate_reactions`` recounts the degeneracy of every
reaction a family generated through its reverse template
(``family.calculate_degeneracy`` regenerates the forward reactions of the new
reactant pair), and one unit can hold thousands of such recounts.  So the
shard unit is finer than a proxy pair: a block of recounts, and a block of
candidates.

Stages (each shard a separate process; ``K`` of each set independently):

``census``   shard k generates the units ``k, k+K, ...`` (in compile order:
             proxies by site type, families sorted) without recounting, and
             records each unit's reaction count and, per recount call in
             order, a content key of the reaction recounted.
``recount``  the recount calls of all units, in order, are split into K
             contiguous blocks; shard k regenerates the units its block
             touches and computes only its own recounts.
``build``    the candidates of all units, in compile order, are split into K
             contiguous blocks; shard k regenerates the units its block
             touches with every recount served from the recount table, and
             compiles the pairs its candidates initiate.
``merge``    replays the single-process candidate loop over the build blocks
             in order: the first build of a pair key is the published one, a
             later candidate of the same pair is resolved against it; then the
             single-process tail (archived junctions, compaction, canonical
             order, certification, ceiling pairs, validation) writes the
             artifact.

The single-process artifact depends on more than the recount values, and the
build reproduces each of those dependencies rather than assuming them away:

- RMG numbers a generated reaction's atoms from one process-wide counter, and
  records order their bond rewrites by those ids.  The census and recount
  stages measure how many ids each unit's generation and each recount assign;
  the build sets the counter to where the single-process discovery had it
  before each unit (``_generation_starts``) and checks where it ends.
- A recount changes the reaction it recounts: ``calculate_degeneracy`` makes
  its reactants independent (new atom ids, resonance structures generated) and
  its regeneration leaves those structures' atoms reordered and labelled, and
  the compiled graphs show it.  The recount stage runs the recount in the same
  two steps (``_Recount.prefix`` and ``rest``) and records what the second did
  to each structure (``recount_effect``); the build runs the first and
  re-applies the record, refusing any structure whose result differs.
- The build generates all its units before it builds its compiler and
  compiles them in order, as the single-process compile does.

Every regeneration is checked against the census: a recount at an ordinal
whose reaction key differs, an extra or missing recount, or a different
reaction count is refused.  Every stage reads the complete output of the
stage before it and refuses a missing, unreadable, incomplete or foreign
(another request, code or shard count) shard by name; outputs are written
atomically, so a failed shard leaves no file.  Nothing is chosen differently:
the units, their order, the recount values and the pair compiles are the
single-process ones, so the merged artifact is the single-process artifact.

CLI::

    python -m rmgpy.kmc.compile_shards {census,recount,build} DATABASE WORKDIR
        --shard k --census-shards K0 --recount-shards K1 --build-shards K2
        (--families FILE.json | --family F) [--deck DECK.py] [--proxy-units N] ...
    python -m rmgpy.kmc.compile_shards merge DATABASE WORKDIR OUTPUT --census-shards K0
        --recount-shards K1 --build-shards K2 (same request arguments)
    python -m rmgpy.kmc.compile_shards plan DATABASE WORKDIR --census-shards K0
        --recount-shards K1 --build-shards K2 (same request arguments)

``plan`` reads the census and prints, per recount and build shard, its units,
recounts or candidates, and an estimated time from the census's measured
generation seconds and the per-recount and per-candidate costs given, and with
``--max-hours H`` exits 1 if any shard is estimated over H hours.
"""

from __future__ import annotations

import argparse
import gzip
import itertools
import json
import os
import pickle
import resource
import sys
import tempfile
import time
import traceback
from collections import defaultdict
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any, Callable, Optional, Sequence

from rmgpy.data.kinetics.common import (
    check_for_same_reactants,
    ensure_independent_atom_ids,
    find_degenerate_reactions,
    generate_molecule_combos,
)
from rmgpy.exceptions import KineticsError
from rmgpy.kmc.compile import (
    ATOM_ID_SPAN,
    FRESH_ATOM_ID,
    advance_atom_ids,
    CompileRequest,
    add_request_arguments,
    load_compile_inputs,
    make_compiler,
    next_atom_id,
    request_from_arguments,
    seek_atom_ids,
)
from rmgpy.kmc.compiler import (
    EventSetCompiler,
    _graph_adjacencies,
    _template,
    compile_candidates,
    compiler_source_hash,
    discovery_entry,
    initiating_record,
    sha256_json,
)

STAGES = ("census", "recount", "build")
# what a recount that is not this shard's returns; the reactions of that
# generation are only counted, never compiled
PLACEHOLDER_DEGENERACY = 1.0


class ShardInputError(RuntimeError):
    """A stage's input is missing, failed, incomplete, foreign or not reproduced."""


@dataclass(frozen=True)
class ShardCounts:
    census: int
    recount: int
    build: int

    def __post_init__(self):
        for stage in STAGES:
            if getattr(self, stage) < 1:
                raise ValueError(f"{stage} shard count must be at least 1")


# --- partitions -------------------------------------------------------------------------


def blocks(n: int, count: int) -> list[tuple[int, int]]:
    """``count`` contiguous [start, stop) blocks covering range(n) in order, sizes
    differing by at most one."""
    return [((k * n) // count, ((k + 1) * n) // count) for k in range(count)]


def round_robin(n: int, shard: int, count: int) -> list[int]:
    """The indices of range(n) shard ``shard`` of ``count`` owns: shard, shard+count, ..."""
    return list(range(shard, n, count))


# --- stage files ------------------------------------------------------------------------


def stage_path(workdir, stage: str, shard: int, count: int) -> Path:
    return Path(workdir) / stage / f"{stage}.{shard}of{count}.pkl.gz"


def write_stage(workdir, stage: str, shard: int, count: int, header, payload) -> Path:
    """Write one shard's output atomically (nothing is left if writing fails)."""
    target = stage_path(workdir, stage, shard, count)
    target.parent.mkdir(parents=True, exist_ok=True)
    document = {"stage": stage, "shard": shard, "of": count, "header": header,
                "payload": payload, "complete": True}
    descriptor, temporary = tempfile.mkstemp(prefix=f".{target.name}.", dir=target.parent)
    try:
        with os.fdopen(descriptor, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb") as stream:
            pickle.dump(document, stream, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(temporary, target)
    finally:
        Path(temporary).unlink(missing_ok=True)
    return target


def read_stage(workdir, stage: str, count: int, header) -> list:
    """The payloads of all ``count`` shards of ``stage``, in shard order; any
    missing, unreadable, incomplete or foreign shard is refused by name.
    ``header`` None skips the request/code comparison."""
    payloads = []
    for shard in range(count):
        path = stage_path(workdir, stage, shard, count)
        name = f"{stage} shard {shard} of {count} ({path})"
        if not path.is_file():
            raise ShardInputError(f"{name} is missing: the shard failed or never ran")
        try:
            with gzip.open(path, "rb") as stream:
                document = pickle.load(stream)
        except Exception as error:
            raise ShardInputError(f"{name} is unreadable ({type(error).__name__}: {error})")
        if not isinstance(document, dict) or document.get("complete") is not True:
            raise ShardInputError(f"{name} is not complete")
        if document.get("of") != count:
            raise ShardInputError(f"{name} was written for {document.get('of')} shards")
        if document.get("stage") != stage or document.get("shard") != shard:
            raise ShardInputError(
                f"{name} holds {document.get('stage')} shard {document.get('shard')}")
        if header is not None:
            for field, value in header.items():
                if document["header"].get(field) != value:
                    raise ShardInputError(
                        f"{name} has another {field} than this stage: it belongs to "
                        "another compile")
        payloads.append(document["payload"])
    return payloads


# --- the recount interception -----------------------------------------------------------


@contextmanager
def intercept_recounts(family, handler: Callable):
    """Route ``family.calculate_degeneracy`` through ``handler(ordinal, reaction,
    compute)`` while the context is open: ``ordinal`` counts this context's
    top-level recounts, ``compute()`` runs the family's own recount (a recount it
    makes itself is not intercepted) and ``compute.prefix()`` only what that
    recount does to ``reaction`` itself."""
    own = vars(family)
    had, previous = "calculate_degeneracy" in own, own.get("calculate_degeneracy")
    original = family.calculate_degeneracy
    ordinals = itertools.count()
    depth = [0]

    def recount(reaction, resonance=True):
        if depth[0]:
            return original(reaction, resonance=resonance)
        return handler(next(ordinals), reaction, _Recount(original, reaction, resonance, depth))

    family.calculate_degeneracy = recount
    try:
        yield
    finally:
        if had:
            family.calculate_degeneracy = previous
        else:
            del family.calculate_degeneracy


class _Recount:
    """One intercepted recount: call it to run the family's own recount, or run
    it as ``prefix()`` then ``rest(same)``: the same steps of
    ``KineticsFamily.calculate_degeneracy`` in the same order, split where the
    recount stops changing what the reaction is and starts regenerating it."""

    def __init__(self, original, reaction, resonance, depth):
        self.original, self.reaction, self.resonance, self.depth = (
            original, reaction, resonance, depth)
        self.family = original.__self__

    def __call__(self):
        self.depth[0] += 1
        try:
            return self.original(self.reaction, resonance=self.resonance)
        finally:
            self.depth[0] -= 1

    def prefix(self):
        """Its reactant list made distinct, their atom ids made independent
        (renumbered from RMG's counter, resonance structures generated), the
        polymer-proxy flag carried; returns the identical-reactant count."""
        reaction = self.reaction
        if reaction.is_charge_transfer_reaction():
            return None
        reactants = reaction.reactants
        poly = any(item.is_polymer_proxy for item in reaction.reactants + reaction.products)
        reactants, same = check_for_same_reactants(reactants)
        ensure_independent_atom_ids(reactants, resonance=self.resonance)
        if poly:
            for item in reactants:
                item.is_polymer_proxy = True
        return same

    def rest(self, same):
        """The forward regeneration from the reactants and its degeneracy."""
        reaction = self.reaction
        if reaction.is_charge_transfer_reaction():
            return 1
        self.depth[0] += 1
        try:
            reactions = []
            for combo in generate_molecule_combos(reaction.reactants):
                reactions.extend(self.family._generate_reactions(
                    combo, products=reaction.products, forward=True,
                    prod_resonance=self.resonance, react_non_reactive=True,
                    apply_species_constraints=False))
            reactions = find_degenerate_reactions(
                reactions, same, template=reaction.template, kinetics_family=self.family,
                resonance=self.resonance)
        finally:
            self.depth[0] -= 1
        if len(reactions) != 1:
            raise KineticsError(
                f"Unable to calculate degeneracy for reaction {reaction} in reaction family "
                f"{self.family.label}. Expected 1 reaction but generated {len(reactions)}")
        return reactions[0].degeneracy


def _participants(reaction):
    return [("reactants", position, item) for position, item in enumerate(reaction.reactants)] + [
        ("products", position, item) for position, item in enumerate(reaction.products)]


def _snapshot(reaction):
    """The participants' structures, their atoms in order, sorting labels and labels."""
    return [(item, [(molecule, list(molecule.atoms), [atom.sorting_label for atom in molecule.atoms],
                     [atom.label for atom in molecule.atoms])
                    for molecule in getattr(item, "molecule", [])])
            for _, _, item in _participants(reaction)]


def _adjacency_digest(molecule) -> str:
    return sha256(molecule.to_adjacency_list().encode("utf-8")).hexdigest()


def recount_effect(reaction, snapshot) -> list:
    """What the rest of a recount did to its reaction's structures since
    ``snapshot``: per changed structure its atom order (as positions in the
    snapshot), atom labels and sorting labels, whether it was re-sorted (its
    connectivity values recomputed), and the digest of the result.  A recount
    that replaced a structure is refused."""
    effect = []
    for (side, position, item), (before, structures) in zip(_participants(reaction), snapshot):
        molecules = getattr(item, "molecule", [])
        if item is not before or len(molecules) != len(structures) or any(
                molecule is not old for molecule, (old, *_) in zip(molecules, structures)):
            raise ShardInputError(
                f"the recount replaced the structures of {side} {position}; the build "
                "cannot reproduce it")
        for index, (molecule, (_, atoms, sorted_before, labels_before)) in enumerate(
                zip(molecules, structures)):
            where = {id(atom): place for place, atom in enumerate(atoms)}
            order = [where[id(atom)] for atom in molecule.atoms]
            labels = [atom.label for atom in molecule.atoms]
            sorting = [atom.sorting_label for atom in molecule.atoms]
            if (order == list(range(len(atoms))) and sorting == sorted_before
                    and labels == labels_before):
                continue
            resorted = any(label < 0 for label in sorted_before) and all(
                label >= 0 for label in sorting)
            effect.append((side, position, index, order, labels, sorting, resorted,
                           _adjacency_digest(molecule)))
    return effect


def apply_recount_effect(label, ordinal, reaction, effect) -> None:
    """Leave the reaction's structures as the full recount did (``recount_effect``)."""
    sides = {"reactants": reaction.reactants, "products": reaction.products}
    for side, position, index, order, labels, sorting, resorted, digest in effect:
        molecule = sides[side][position].molecule[index]
        atoms = list(molecule.atoms)
        molecule.atoms = [atoms[place] for place in order]
        for atom, atom_label, sorting_label in zip(molecule.atoms, labels, sorting):
            atom.label = atom_label
            atom.sorting_label = sorting_label
        if resorted:
            molecule.update_connectivity_values()
        if _adjacency_digest(molecule) != digest:
            raise ShardInputError(
                f"{label}: recount {ordinal} leaves {side} {position} structure {index} "
                "otherwise than the full recount did: the recount has an effect the build "
                "does not reproduce")


def recount_key(reaction) -> str:
    """Content key of a reaction whose degeneracy is recounted: family, template
    and every resonance structure of every reactant and product."""
    sides = [
        [[molecule.to_adjacency_list() for molecule in getattr(item, "molecule", [item])]
         for item in side]
        for side in (reaction.reactants, reaction.products)
    ]
    text = json.dumps([reaction.family, sorted(map(str, reaction.template)), sides])
    return sha256(text.encode("utf-8")).hexdigest()


def _checked(label, keys, ordinal, key):
    if ordinal >= len(keys):
        raise ShardInputError(
            f"{label}: recount {ordinal} was not in the census (census recorded "
            f"{len(keys)}): the generation is not reproducible")
    if keys[ordinal] != key:
        raise ShardInputError(
            f"{label}: recount {ordinal} is not the census's (key {key[:12]} != "
            f"{keys[ordinal][:12]}): the generation is not reproducible")


def recording(keys: list, prefix_ids: Optional[list] = None, key=recount_key):
    """Census: record each recount's key (and the atom ids its prefix assigns);
    recount nothing."""
    def handler(ordinal, reaction, compute):
        keys.append(key(reaction))
        if prefix_ids is not None:
            before = next_atom_id()
            compute.prefix()
            prefix_ids.append(_ids_between(before, next_atom_id()))
        return PLACEHOLDER_DEGENERACY
    return handler


def computing(label, keys, owned, results: dict, consumed: Optional[dict] = None,
              effects: Optional[dict] = None, key=recount_key):
    """Recount: check every recount against the census; compute the ``owned``
    ordinals into ``results`` (with the atom ids each assigned into ``consumed``
    and what its regeneration did to its reaction into ``effects``); every other
    recount only changes its reaction as its prefix does."""
    def handler(ordinal, reaction, compute):
        _checked(label, keys, ordinal, key(reaction))
        if effects is None:
            if ordinal in owned:
                results[ordinal] = compute()
                return results[ordinal]
            return PLACEHOLDER_DEGENERACY
        if ordinal not in owned:
            compute.prefix()
            return PLACEHOLDER_DEGENERACY
        before = next_atom_id()
        same = compute.prefix()
        snapshot = _snapshot(reaction)
        results[ordinal] = compute.rest(same)
        consumed[ordinal] = _ids_between(before, next_atom_id())
        effects[ordinal] = recount_effect(reaction, snapshot)
        return results[ordinal]
    return handler


def serving(label, keys, table, replay: Optional[dict] = None, key=recount_key):
    """Build: check every recount against the census; serve its value from
    ``table``.  With ``replay`` ({ordinal: (prefix ids, recount ids, effect)}),
    leave the reaction and RMG's atom-id counter as the full recount did."""
    def handler(ordinal, reaction, compute):
        _checked(label, keys, ordinal, key(reaction))
        if replay is not None:
            prefix_ids, recount_ids, effect = replay[ordinal]
            compute.prefix()
            apply_recount_effect(label, ordinal, reaction, effect)
            if recount_ids < prefix_ids:
                raise ShardInputError(
                    f"{label}: recount {ordinal} assigned fewer atom ids ({recount_ids}) "
                    f"than its prefix ({prefix_ids})")
            advance_atom_ids(recount_ids - prefix_ids)
        return table[ordinal]
    return handler


def _rss() -> str:
    """This process's peak resident memory (for the per-shard memory request)."""
    return f"peak rss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.0f} MB"


def _ids_between(before: int, after: int) -> int:
    """Atom ids assigned between two ``next_atom_id`` readings (not the readings)."""
    return (after - before - 1) % ATOM_ID_SPAN


def _atom_id(value: int) -> int:
    return (value - FRESH_ATOM_ID) % ATOM_ID_SPAN + FRESH_ATOM_ID


# --- the plan: units, census, recount table -------------------------------------------------


@dataclass
class _Plan:
    inputs: Any
    units: list            # [(site_type, family)] in compile order
    discovery_order: list  # unit indices in the single-process discovery order
    discovery_start: int   # the atom id the single-process discovery starts from
    proxies: dict          # site_type -> the discovery proxy
    universe: set
    candidates: list
    header: dict


def _request_identity(request: CompileRequest) -> dict:
    data = asdict(request)
    data.pop("output_dir")
    return json.loads(json.dumps(data))


def _plan(request: CompileRequest, database=None) -> _Plan:
    inputs = load_compile_inputs(request, database)
    # where the single-process compile is when its discovery starts
    discovery_start = next_atom_id()
    universe, candidates = EventSetCompiler.discovery_candidates(
        inputs.database.kinetics, request.families, inputs.family_universe)
    proxies = {}
    for proxy in inputs.proxies:
        if proxy.site_type in proxies:
            raise ValueError(f"two proxies share the site type {proxy.site_type!r}")
        proxies[proxy.site_type] = proxy
    units = [
        (site, family)
        for site in sorted(proxies)
        for family in EventSetCompiler.proxy_family_candidates(proxies[site], candidates)
    ]
    index = {unit: position for position, unit in enumerate(units)}
    discovery_order = [
        index[(proxy.site_type, family)]
        for proxy in inputs.proxies
        for family in EventSetCompiler.proxy_family_candidates(proxy, candidates)
    ]
    header = {"request": _request_identity(request), "code": compiler_source_hash(),
              "units": [list(unit) for unit in units], "discovery_start": discovery_start}
    return _Plan(inputs, units, discovery_order, discovery_start, proxies, universe,
                 candidates, header)


def _label(plan, index) -> str:
    site, family = plan.units[index]
    return f"unit {index} ({site}, {family})"


def _generate(plan, index, handler):
    site, family = plan.units[index]
    kinetics = plan.inputs.database.kinetics
    with intercept_recounts(kinetics.families[family], handler):
        return EventSetCompiler.generate_on_proxy(kinetics, plan.proxies[site], family)


def _census(workdir, counts: ShardCounts, plan) -> list[dict]:
    """The census entry of every unit, by unit index (each exactly once)."""
    entries = [None] * len(plan.units)
    for shard, payload in enumerate(read_stage(workdir, "census", counts.census, plan.header)):
        owned = round_robin(len(plan.units), shard, counts.census)
        if [unit["index"] for unit in payload["units"]] != owned:
            raise ShardInputError(
                f"census shard {shard} of {counts.census} does not hold its units {owned}")
        for unit in payload["units"]:
            entries[unit["index"]] = unit
    return entries


def _recount_jobs(census) -> list[tuple[int, int]]:
    return [(index, ordinal) for index, unit in enumerate(census)
            for ordinal in range(len(unit["recount_keys"]))]


def _candidates(census) -> list[tuple[int, int]]:
    return [(index, position) for index, unit in enumerate(census)
            for position in range(unit["reactions"])]


def _owned(items, block) -> dict[int, list[int]]:
    owned = defaultdict(list)
    for index, item in items[block[0]:block[1]]:
        owned[index].append(item)
    return owned


def _recount_table(workdir, counts: ShardCounts, plan, census):
    """(every census recount's value, and its replay (prefix ids, recount ids,
    effect), by unit and ordinal), each recount held by exactly its own shard."""
    jobs = _recount_jobs(census)
    table, replay = defaultdict(dict), defaultdict(dict)
    parts = read_stage(workdir, "recount", counts.recount, plan.header)
    for shard, (payload, block) in enumerate(zip(parts, blocks(len(jobs), counts.recount))):
        expected = {index: sorted(ordinals) for index, ordinals in _owned(jobs, block).items()}
        held = [{index: sorted(values) for index, values in payload[field].items()}
                for field in ("recounts", "consumed", "effects")]
        if payload["block"] != list(block) or any(found != expected for found in held):
            raise ShardInputError(
                f"recount shard {shard} of {counts.recount} does not hold its recounts "
                f"{block[0]}..{block[1]}")
        for index, values in payload["recounts"].items():
            table[index].update(values)
            prefix_ids = census[index]["prefix_ids"]
            for ordinal in values:
                replay[index][ordinal] = (prefix_ids[ordinal], payload["consumed"][index][ordinal],
                                          payload["effects"][index][ordinal])
    return table, replay


def _generation_starts(plan, census, replay) -> dict[int, tuple[int, int]]:
    """(start, end): the atom ids each unit's generation starts from and leaves
    the counter at in the single-process discovery, units in its order; a
    unit assigns its own ids and its recounts' in full."""
    starts, position = {}, plan.discovery_start
    for index in plan.discovery_order:
        extra = sum(recount - prefix for prefix, recount, _ in replay[index].values())
        end = _atom_id(position + census[index]["atom_ids"] + extra)
        starts[index] = (position, end)
        position = end
    return starts


def _check_generation(plan, index, census, reactions, seen):
    unit = census[index]
    if len(reactions) != unit["reactions"] or seen != len(unit["recount_keys"]):
        raise ShardInputError(
            f"{_label(plan, index)} generated {len(reactions)} reactions with {seen} "
            f"recounts; the census recorded {unit['reactions']} with "
            f"{len(unit['recount_keys'])}: the generation is not reproducible")


def _counted(handler, counter):
    def counting(ordinal, reaction, compute):
        counter[0] = ordinal + 1
        return handler(ordinal, reaction, compute)
    return counting


def _compiler(request, plan, census):
    """The compiler a single-process compile builds after discovering ``census``."""
    active, scheduled, generated = set(), set(), defaultdict(set)
    for (site, family), unit in zip(plan.units, census):
        scheduled.add(family)
        if unit["reactions"]:
            active.add(family)
            generated[site].add(family)
    for proxy in plan.inputs.proxies:
        proxy.metadata["generated_families"] = sorted(generated[proxy.site_type])
    active, excluded = EventSetCompiler.family_accounting(
        plan.universe, plan.candidates, scheduled, active,
        plan.inputs.discovery_options.get("proxy_set_name"))
    compiler = make_compiler(request, plan.inputs, active, excluded, {})
    compiler._reset_applicability()
    return compiler


# --- stages -------------------------------------------------------------------------------


def run_census(request, workdir, shard, counts: ShardCounts, *, database=None):
    plan = _plan(request, database)
    units = []
    for index in round_robin(len(plan.units), shard, counts.census):
        keys, prefix_ids = [], []
        started, before = time.time(), next_atom_id()
        reactions = _generate(plan, index, recording(keys, prefix_ids))
        units.append({"index": index, "reactions": len(reactions), "recount_keys": keys,
                      "prefix_ids": prefix_ids,
                      # less the two readings around each prefix
                      "atom_ids": _ids_between(before, next_atom_id()) - 2 * len(prefix_ids),
                      "seconds": round(time.time() - started, 3)})
        print(f"census {_label(plan, index)}: {len(reactions)} reactions, {len(keys)} "
              f"recounts, {units[-1]['seconds']:.1f}s, {_rss()}", flush=True)
    write_stage(workdir, "census", shard, counts.census, plan.header, {"units": units})
    return plan.inputs.database


def run_recount(request, workdir, shard, counts: ShardCounts, *, database=None):
    plan = _plan(request, database)
    census = _census(workdir, counts, plan)
    jobs = _recount_jobs(census)
    block = blocks(len(jobs), counts.recount)[shard]
    recounts, consumed, effects = {}, {}, {}
    for index, ordinals in sorted(_owned(jobs, block).items()):
        values, ids, done, seen = {}, {}, {}, [0]
        started = time.time()
        reactions = _generate(plan, index, _counted(
            computing(_label(plan, index), census[index]["recount_keys"], set(ordinals),
                      values, ids, done), seen))
        _check_generation(plan, index, census, reactions, seen[0])
        recounts[index] = {ordinal: values[ordinal] for ordinal in ordinals}
        consumed[index] = {ordinal: ids[ordinal] for ordinal in ordinals}
        effects[index] = {ordinal: done[ordinal] for ordinal in ordinals}
        print(f"recount {_label(plan, index)}: {len(ordinals)} of "
              f"{len(census[index]['recount_keys'])} in {time.time() - started:.1f}s, {_rss()}",
              flush=True)
    write_stage(workdir, "recount", shard, counts.recount, plan.header,
                {"block": list(block), "recounts": recounts, "consumed": consumed,
                 "effects": effects})
    return plan.inputs.database


def run_build(request, workdir, shard, counts: ShardCounts, *, database=None):
    plan = _plan(request, database)
    census = _census(workdir, counts, plan)
    table, replay = _recount_table(workdir, counts, plan, census)
    starts = _generation_starts(plan, census, replay)
    candidates = _candidates(census)
    block = blocks(len(candidates), counts.build)[shard]
    owned = sorted(_owned(candidates, block).items())
    # as in a single-process compile: every generation, then the compiler, then
    # the candidate loop
    generated = {}
    for index, _ in owned:
        seen = [0]
        started = time.time()
        # the reaction atoms get the ids the single-process discovery gave them,
        # and each recount leaves its reaction as the full recount did
        start, end = starts[index]
        seek_atom_ids(start)
        reactions = _generate(plan, index, _counted(
            serving(_label(plan, index), census[index]["recount_keys"], table[index],
                    replay[index]), seen))
        _check_generation(plan, index, census, reactions, seen[0])
        reached = next_atom_id()
        if reached != end:
            raise ShardInputError(
                f"{_label(plan, index)}: the generation left RMG's atom ids at {reached}, "
                f"the single-process discovery at {end}")
        generated[index] = compile_candidates(reactions)
        print(f"build {_label(plan, index)}: generated in {time.time() - started:.1f}s, {_rss()}",
              flush=True)
    compiler = _compiler(request, plan, census)
    proxy_inputs, provenance = compiler._compile_provenance()
    proxies = {proxy.site_type: proxy for proxy in compiler.proxies}
    entries, built = [], set()
    for index, positions in owned:
        started = time.time()
        proxy = proxies[plan.units[index][0]]
        for position in positions:
            reaction, sides = generated[index][position]
            entries.append(_build_entry(compiler, proxy, reaction, sides, provenance, built))
        print(f"build {_label(plan, index)}: {len(positions)} of {len(generated[index])} "
              f"candidates in {time.time() - started:.1f}s, {_rss()}", flush=True)
    write_stage(workdir, "build", shard, counts.build, plan.header,
                {"block": list(block), "provenance_sha256": sha256_json(provenance),
                 "entries": entries})
    return plan.inputs.database


def _build_entry(compiler, proxy, reaction, sides, provenance, built):
    """One candidate of the single-process loop, its steps in that loop's order;
    the pair is compiled when this block has not compiled it yet (the merge
    keeps the first compile in order)."""
    oriented, coverage_site, excluded = compiler._candidate_channel(proxy, reaction)
    if excluded is not None:
        return {"excluded_channel": excluded}
    pair_key = (reaction.family, _template(reaction), proxy.frontier, sides)
    entry = {"pair_key": pair_key}
    if pair_key not in built:
        built.add(pair_key)
        start = len(compiler._applicability_refusals)
        pair_records, initiating = compiler._linked_family_pair(proxy, reaction, provenance)
        entry["build"] = {
            "records": pair_records,
            "initiating": initiating.event_id if initiating else None,
            "refusals": compiler._applicability_refusals[start:],
        }
    entry["reactant_graphs"] = _graph_adjacencies(oriented.reactants)
    entry["discovery"] = discovery_entry(proxy, reaction, coverage_site, None)
    return entry


def run_merge(request, workdir, counts: ShardCounts, *, database=None):
    """Write the artifact of the build shards; returns (path, artifact id, database)."""
    plan = _plan(request, database)
    census = _census(workdir, counts, plan)
    parts = read_stage(workdir, "build", counts.build, plan.header)
    candidates = _candidates(census)
    compiler = _compiler(request, plan, census)
    proxy_inputs, provenance = compiler._compile_provenance()
    for shard, (payload, block) in enumerate(zip(parts, blocks(len(candidates), counts.build))):
        if payload["block"] != list(block) or len(payload["entries"]) != block[1] - block[0]:
            raise ShardInputError(
                f"build shard {shard} of {counts.build} does not hold its candidates "
                f"{block[0]}..{block[1]}")
        if payload["provenance_sha256"] != sha256_json(provenance):
            raise ShardInputError(
                f"build shard {shard} of {counts.build} compiled under another provenance")
    records, discovery, excluded_channels, refusals, pairs = [], [], [], [], {}
    for shard, payload in enumerate(parts):
        for entry in payload["entries"]:
            if "excluded_channel" in entry:
                excluded_channels.append(entry["excluded_channel"])
                continue
            key = entry["pair_key"]
            if key not in pairs:
                build = entry.get("build")
                if build is None:
                    raise ShardInputError(
                        f"build shard {shard} of {counts.build} met pair {key[:2]} first "
                        "without compiling it")
                pairs[key] = build["records"]
                records.extend(build["records"])
                refusals.extend(build["refusals"])
                event_id = build["initiating"]
            else:
                initiating = initiating_record(pairs[key], entry["reactant_graphs"])
                event_id = initiating.event_id if initiating else None
            discovery.append({**entry["discovery"], "event_id": event_id})
    compiler._applicability_refusals = refusals
    compiler._finish_artifact(proxy_inputs, provenance, records, discovery, excluded_channels)
    print(f"merging {len(candidates)} candidates into {len(records)} pair records, {_rss()}", flush=True)
    path, artifact_id = compiler._write_artifact_file(request.output_dir)
    return Path(path), artifact_id, plan.inputs.database


RUNNERS = {"census": run_census, "recount": run_recount, "build": run_build}


def run_sharded_in_process(request, workdir, counts: ShardCounts, *, database=None):
    """Every shard of every stage in turn, then the merge, in this process
    (``database``: reuse one loaded for the same request)."""
    for stage in STAGES:
        for shard in range(getattr(counts, stage)):
            database = RUNNERS[stage](request, workdir, shard, counts, database=database)
    return run_merge(request, workdir, counts, database=database)


# --- plan: per-shard estimate ----------------------------------------------------------------


def estimate(workdir, counts: ShardCounts, plan, recount_seconds: float,
             candidate_seconds: float) -> dict[str, list[dict]]:
    """Per recount and build shard: units touched, items, estimated seconds
    (census generation seconds of every unit touched + items x unit cost)."""
    census = _census(workdir, counts, plan)
    out = {}
    for stage, items, count, cost in (
            ("recount", _recount_jobs(census), counts.recount, recount_seconds),
            ("build", _candidates(census), counts.build, candidate_seconds)):
        rows = []
        for block in blocks(len(items), count):
            owned = _owned(items, block)
            generation = sum(census[index]["seconds"] for index in owned)
            rows.append({"units": len(owned), "items": block[1] - block[0],
                         "seconds": generation + (block[1] - block[0]) * cost})
        out[stage] = rows
    return out


# --- CLI ----------------------------------------------------------------------------------


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m rmgpy.kmc.compile_shards")
    parser.add_argument("stage", choices=STAGES + ("merge", "plan"))
    parser.add_argument("database")
    parser.add_argument("workdir")
    parser.add_argument("output", nargs="?", help="artifact directory (merge)")
    parser.add_argument("--shard", type=int)
    parser.add_argument("--census-shards", type=int, required=True)
    parser.add_argument("--recount-shards", type=int, required=True)
    parser.add_argument("--build-shards", type=int, required=True)
    parser.add_argument("--recount-seconds", type=float, default=10.5,
                        help="plan: estimated seconds per recount")
    parser.add_argument("--candidate-seconds", type=float, default=2.0,
                        help="plan: estimated seconds per candidate compile")
    parser.add_argument("--max-hours", type=float,
                        help="plan: exit 1 if any shard is estimated over this")
    add_request_arguments(parser)
    args = parser.parse_args(argv)
    counts = ShardCounts(args.census_shards, args.recount_shards, args.build_shards)
    if args.stage in STAGES:
        if args.shard is None or not 0 <= args.shard < getattr(counts, args.stage):
            parser.error(f"{args.stage} needs --shard in 0..{getattr(counts, args.stage) - 1}")
    if args.stage == "merge" and args.output is None:
        parser.error("merge needs OUTPUT")
    try:
        request = request_from_arguments(parser, args, args.output or args.workdir)
        if args.stage in STAGES:
            RUNNERS[args.stage](request, args.workdir, args.shard, counts)
        elif args.stage == "merge":
            path, artifact_id, _ = run_merge(request, args.workdir, counts)
            print(path)
            print(artifact_id)
        else:
            rows = estimate(args.workdir, counts, _plan(request), args.recount_seconds,
                            args.candidate_seconds)
            over = []
            for stage, shards in rows.items():
                worst = max(shards, key=lambda row: row["seconds"])
                print(f"{stage}: {len(shards)} shards, at most {worst['items']} items over "
                      f"{worst['units']} units, estimated at most {worst['seconds'] / 3600:.2f} h")
                for shard, row in enumerate(shards):
                    print(f"  {stage} {shard}: {row['units']} units, {row['items']} items, "
                          f"{row['seconds'] / 3600:.2f} h")
                    if args.max_hours is not None and row["seconds"] > 3600 * args.max_hours:
                        over.append(f"{stage} {shard}")
            if over:
                print(f"over {args.max_hours} h: {', '.join(over)}; raise the shard counts",
                      file=sys.stderr)
                return 1
    except Exception:
        traceback.print_exc()
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
