"""Certify an artifact once, streaming, and write every emission shard's inputs (I-295).

Every emission shard used to parse the whole artifact and certify every
record, then emit its ``1/N`` slice: its memory and time grew with the whole
catalog (two parses of the text, every record certificate held), and more
shards did not help.  A 5.3 GB artifact failed every shard inside the JSON
parse.

Certification cannot simply be split: a stored certificate records each
family-tree node as its atoms in their in-memory order, and RMG's matcher
re-sorts a group's atoms in place when it matches it (I-251).  So a record's
certificate bytes depend on every record certified before it, in artifact
order, and a shard certifying only its own records publishes other bytes.
The certification is therefore run once, here, in artifact order (the history
every whole-artifact shard had), but streaming:

1. ``artifact_stream.build_index`` (or a verified ``emit index`` file): the
   head, and per record its byte span and links, in one pass.
2. Certify each record in artifact order (``certificates.certify_one``), read
   by span with its rate owner; the certificate is spilled to disk and only
   ``certified`` and the first gap are kept.  The certifier's caches are
   bounded (PP: 1.9 GB after 1800 of 2410 records, flat after).
3. ``certificates.refusal_reasons`` over the records' status, links and rate
   sources; each refused record's new id (``republish_refused``).
4. The shard assignment (``assign_shards``): each record's reference-case key
   (``emit.record_key``: family, template-side species, template, other-side
   species) is read once more, and the records of one key go to one shard, the
   keys placed largest first on the least loaded shard by an estimated cost
   (records, and the roster subsets whose participants the key admits:
   ``emit.CaseFilter``).  A shard then builds only its own keys' reference
   cases, where ``sorted(ids)[k::n]`` scattered every key over every shard.  A
   key costing more than one shard's share is split over several shards (each
   rebuilds that key's cases), named in the manifest; a record whose key
   cannot be read is its own unit.
5. Each record again, in order: republished, relinked and its certificate's
   owner bound (``relink``, ``publish_certificate``), written to the input
   of the shard owning it with its certificate, and to the inputs of the
   shards owning a record that links to it (emission reads a record's rate
   owner in its published form).

These are the steps of ``certificates.certify_records``, which runs them on a
list in memory.  A shard (``emit --prepared DIR``) reads the manifest and its
own input only: its records, their certificates and the records they link,
then emits; its rows are the whole-artifact shard's rows.

The manifest binds the inputs to the artifact's bytes, the roster, the
database content, the code and the shard count; a shard refuses any other.  It
carries the assignment (the shard of every artifact id), each shard input's
sha256 (checked on the bytes read, before anything is unpickled) and a digest
of each shard's (source id, published id) pairs, which the writer computes
from what it wrote and the reader from what it read.  A ``kmc-emit-prepared/1``
manifest (``sorted(ids)[k::n]``, no digests) is refused by name.

**In blocks (I-304b).**  The same steps, the certification split into J
processes: ``write_plan`` (``emit prepare-plan``) fixes J contiguous ranges of
artifact order; ``prepare_block`` (``emit prepare-block``, one fresh process
per block) certifies its range (rate owners read from the artifact wherever
they are) and computes its records' keys (steps 2 and 4 without the
placement), then seals its files exclusively; ``prepare_finish`` (``emit
prepare-finish``) proves every block's exact range from its bytes, runs steps
3, 4 (placement) and 5 over all records in artifact order, and writes a new,
never rewritten generation directory.  A block's certificate differs from the
serial one at most in its tree nodes' atom order (no presort); every other
byte of the shard inputs is the serial prepare's (``compare_prepared``).  The
design is ``pk-i304b/DESIGN.md``.
"""

from __future__ import annotations

import gzip
import hashlib
import heapq
import io
import json
import math
import os
import pickle
import shutil
import tempfile
import time
import zlib
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from rmgpy.kmc import artifact_stream as astream
from rmgpy.kmc import certificates as cm
from rmgpy.kmc import emit
from rmgpy.kmc.emit import PreparedRefusal

MANIFEST_SCHEMA = "kmc-emit-prepared/2"
OLD_MANIFEST_SCHEMAS = ("kmc-emit-prepared/1",)
MANIFEST = "manifest.json"
# the estimated cost of a shard's work: one per record emitted, and per roster
# subset whose member cases a key makes the shard build (the balance only; no
# row depends on it)
COST_WEIGHTS = {"record": 1.0, "subset": 0.05}
# records between resets of the key pass's species cache (its size only)
SPECIES_RESET = 512


def shard_file(k: int, n: int) -> str:
    return f"shard.{k}of{n}.pkl.gz"


def _light(d: dict) -> dict:
    """What ``refusal_reasons`` reads of a record: id, status, link, rate source."""
    out = {key: d[key] for key in ("event_id", "status", "reverse_of") if key in d}
    rs = d.get("rate_source")
    if rs is not None:
        out["rate_source"] = {key: rs[key] for key in ("kind", "reference") if key in rs} \
            if isinstance(rs, dict) else rs
    return out


def identity(database, roster) -> dict:
    """What a shard's rows depend on besides the artifact: the roster, the
    database content and the code."""
    return {"roster_sha256": roster.sha256,
            "database_content_sha256": emit.database_content_sha256(Path(database)),
            "code_sources_sha256": emit.code_sources_sha256(),
            "code_extensions_sha256": emit.code_extensions_sha256()}


def _rss() -> int:
    return emit.rss_mb()


def _one(artifact, entry):
    from rmgpy.kmc.event_record import EventRecord
    return EventRecord.from_dict(astream.read_records(artifact, [entry])[0])


@dataclass
class Unit:
    """Records placed together: those of one reference-case key (``name``,
    its canonical text) or one record whose key could not be read."""
    name: str
    ids: list
    subsets: int = 0
    pieces: list = field(default_factory=list)


def unit_cost(records: int, subsets: int, weights=None) -> float:
    w = weights or COST_WEIGHTS
    return records * w["record"] + subsets * w["subset"]


def assign_shards(units, shards: int, weights=None):
    """(shard of each record id, plan): every unit on one shard, largest
    estimated cost first, each on the least loaded shard (ties: the lowest
    shard; order: cost, then name).  A unit costing more than a shard's share
    (the total over ``shards``) with more than one record is split into
    ``min(records, shards, ceil(cost / share))`` contiguous runs of its sorted
    ids, each run carrying the unit's whole case cost (each shard holding a run
    builds the key's cases) and each on a distinct shard; the plan names every
    split unit (with its shards) and every unit
    whose case cost alone exceeds the share.  Deterministic in the units."""
    if shards < 1:
        raise ValueError("the shard count must be at least 1")
    units = sorted(units, key=lambda u: u.name)
    names = [u.name for u in units]
    if len(set(names)) != len(names):
        raise ValueError("two assignment units share a name")
    total = sum(unit_cost(len(u.ids), u.subsets, weights) for u in units)
    share = total / shards if total else 0.0
    pieces, split, case_bound = [], [], []
    for u in units:
        ids = sorted(u.ids)
        cost = unit_cost(len(ids), u.subsets, weights)
        parts = 1
        if share and cost > share and len(ids) > 1:
            parts = min(len(ids), shards, math.ceil(cost / share))
        if unit_cost(0, u.subsets, weights) > share and share:
            case_bound.append(u.name)
        size = math.ceil(len(ids) / parts) if ids else 0
        runs = [ids[i:i + size] for i in range(0, len(ids), size)] if size else [ids]
        if len(runs) > 1:
            split.append({"unit": u.name, "records": len(ids), "parts": len(runs)})
        for part, run in enumerate(runs):
            pieces.append((-unit_cost(len(run), u.subsets, weights), u.name, part, run))
    pieces.sort(key=lambda p: (p[0], p[1], p[2]))
    loads = [(0.0, k) for k in range(shards)]
    heapq.heapify(loads)
    shard_of, load_of, placed = {}, [0.0] * shards, defaultdict(list)
    for negcost, name, part, run in pieces:
        # the least loaded shard not already holding a run of this unit (a
        # split unit has at most ``shards`` runs, so one exists)
        held = []
        load, k = heapq.heappop(loads)
        while k in placed[name]:
            held.append((load, k))
            load, k = heapq.heappop(loads)
        for item in held:
            heapq.heappush(loads, item)
        for i in run:
            shard_of[i] = k
        load_of[k] = load - negcost
        placed[name].append(k)
        heapq.heappush(loads, (load - negcost, k))
    for entry in split:
        entry["shards"] = sorted(placed[entry["unit"]])
    plan = {"units": len(units), "split": split, "case_bound": case_bound,
            "loads": [round(x, 3) for x in load_of], "share": round(share, 3),
            "weights": dict(weights or COST_WEIGHTS)}
    return shard_of, plan


def subset_shapes(fam, systems) -> Counter:
    """(family, direction, component formulas) -> the number of roster subsets
    (``emit.roster_subsets``) with those participants' formulas."""
    formulas, out = {}, Counter()
    for name, _family, direction, chosen, _subset in emit.roster_subsets(fam, systems):
        out[(name, direction, emit.subset_formulas(chosen, formulas))] += 1
    return out


def window_start(position: int) -> int:
    """The first artifact position of the key pass's species-cache window
    holding ``position``: the pass starts a fresh cache at 0 and before every
    record whose 1-based count is a multiple of ``SPECIES_RESET``."""
    if position < SPECIES_RESET - 1:
        return 0
    return (position + 1) // SPECIES_RESET * SPECIES_RESET - 1


def record_keys(artifact, universe, certifier, lo: int = 0, hi=None):
    """(entry, key name or None, template-side molecules, other-side
    molecules) of ``universe[lo:hi]``, in order: ``emit.record_parts`` with the
    species cache reset where the whole-universe pass resets it.  The cache
    decides whether a key collision raises (the record is then unkeyed), so a
    pass starting at ``lo`` first re-keys, silently, from the reset point at or
    before it: its cache then holds what the whole pass's held there."""
    from rmgpy.kmc import channels as chm
    hi = len(universe) if hi is None else hi
    species = chm.SpeciesKeys()
    for position in range(window_start(lo), hi):
        if (position + 1) % SPECIES_RESET == 0:
            species = chm.SpeciesKeys()
        e = universe[position]
        d = _one(artifact, e).to_dict()
        try:
            key, mols, others = emit.record_parts(certifier, d, species)
        except Exception as error:   # its own unit; the shard's row reports it
            emit.reraise_resource(error)
            key = mols = others = None
        if position >= lo:
            name = None if key is None else json.dumps(key, separators=(",", ":"))
            yield e, name, mols, others


def key_subsets(name, mols, others, shapes) -> int:
    """The roster subsets a key admits (``subset_shapes``), from one of its
    records' molecules (a key fixes their species, so any of its records)."""
    family = json.loads(name)[0]
    return shapes[(family, "forward", emit.component_formulas(mols))] + \
        shapes[(family, "reverse", emit.component_formulas(others))]


def units_of(keyed, subsets) -> list:
    """``Unit`` list of (event id, key name or None) in artifact order: one
    per key (``subsets``: key name -> its subset count), one per unkeyed
    record."""
    by_key, units = defaultdict(list), []
    for event_id, name in keyed:
        if name is None:
            units.append(Unit(f"unkeyed:{event_id}", [event_id]))
        else:
            by_key[name].append(event_id)
    return units + [Unit(name, ids, subsets[name]) for name, ids in by_key.items()]


def assignment_units(artifact, universe, certifier, fam, systems) -> list:
    """The units of ``universe`` (``Unit``): records grouped by
    ``emit.record_key``, each key's subsets counted from the roster."""
    shapes = subset_shapes(fam, systems) if systems else Counter()
    keyed, subsets = [], {}
    for e, name, mols, others in record_keys(artifact, universe, certifier):
        keyed.append((e.event_id, name))
        if name is not None and name not in subsets:
            subsets[name] = key_subsets(name, mols, others, shapes)
    return units_of(keyed, subsets)


def own_digest(pairs) -> str:
    """sha256 of a shard's sorted (source id, published id) pairs."""
    return hashlib.sha256(json.dumps(sorted(pairs), separators=(",", ":")).encode()).hexdigest()


def file_sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


@dataclass
class Context:
    """What every prepare step reads first: the roster, the verified index,
    the declaration binding and the records selected (artifact order)."""
    roster: object
    index: object
    bound_by: tuple
    universe: list


def open_context(artifact, roster_path, families=None, index_path=None,
                 legacy_artifact: bool = False) -> Context:
    roster = emit.load_roster(roster_path)
    index = astream.index_for(artifact, index_path)
    bound_by = emit.bind_artifact(index.head, roster.digests(), f"the roster {roster_path}",
                                  legacy_artifact)
    emit.check_roster_for_artifact(roster, index.head)
    universe = [e for e in index.entries if not families or e.family in families]
    return Context(roster, index, bound_by, universe)


def certify_range(artifact, universe, entry_of, certifier, lo, hi, spill, light, summaries,
                  label="certified"):
    """Certify ``universe[lo:hi]`` in order, each through its rate owner read
    from the artifact (a partner anywhere in ``entry_of``, never its
    certificate); pickle each certificate to ``spill`` (with ``item``'s
    wrapping, if given) and keep its record's ``_light`` and summary."""
    t = time.time()
    for count, e in enumerate(universe[lo:hi], 1):
        record = _one(artifact, e)
        d = record.to_dict()
        light[e.event_id] = _light(d)
        owner = None
        if (d.get("rate_source") or {}).get("kind") == "reference-thermo reverse":
            partner = entry_of.get(d.get("reverse_of"))
            owner = None if partner is None else _one(artifact, partner)
        rc = cm.certify_one(certifier, record, d, owner)
        spill(e, rc)
        summaries.append({"certified": rc.get("certified"),
                          "gaps": list(rc.get("gaps") or ())[:1]})
        if count % 1000 == 0:
            print(f"PREPARE {label} {count}/{hi - lo} "
                  f"{time.time() - t:.0f}s rss={_rss()}MB", flush=True)


def gate(artifact, universe, entry_of, light, summaries):
    """(reasons, replacements, final ids): ``refusal_reasons`` over every
    record in artifact order, each refused record's republished id."""
    reasons = cm.refusal_reasons(
        ((light[e.event_id], s) for e, s in zip(universe, summaries)), light)
    replacements = {i: cm.republish_refused(_one(artifact, entry_of[i]), reason).event_id
                    for i, reason in reasons.items()}
    final_ids = {replacements.get(e.event_id, e.event_id) for e in universe}
    return reasons, replacements, final_ids


def publish(artifact, universe, entry_of, certificates, reasons, replacements, final_ids,
            shard_of, shards, work):
    """Write every shard input into ``work``, in artifact order: each record
    republished, relinked and its certificate bound (``certificates``: one per
    record of ``universe``, in order), to the shard owning it and to the
    shards owning a record linking to it.  Returns (files, items written,
    owned (source, published) pairs per shard)."""
    readers = defaultdict(set)
    for e in universe:
        for target in e.links:
            if target in entry_of:
                readers[target].add(shard_of[e.event_id])
    files = [work / shard_file(k, shards) for k in range(shards)]
    written = [0] * shards
    owned = [[] for _ in range(shards)]
    writers = [gzip.open(f, "wb", compresslevel=1) for f in files]
    try:
        for e, rc in zip(universe, certificates):
            record = _one(artifact, e)
            if e.event_id in reasons:
                record = cm.republish_refused(record, reasons[e.event_id])
            if reasons:
                record = cm.relink(record, replacements)
            published, rc, _ref = cm.publish_certificate(record, rc, replacements, final_ids)
            d = published.to_dict()
            own = shard_of[e.event_id]
            owned[own].append([e.event_id, d["event_id"]])
            for k in sorted({own} | readers[e.event_id]):
                pickle.dump({"source": e.event_id, "own": k == own, "record": d,
                             "certificate": rc if k == own else None},
                            writers[k], protocol=pickle.HIGHEST_PROTOCOL)
                written[k] += 1
        for k, writer in enumerate(writers):
            pickle.dump({"end": written[k]}, writer, protocol=pickle.HIGHEST_PROTOCOL)
    finally:
        for writer in writers:
            writer.close()
    ids = [e.event_id for e in universe]
    if sorted(i for pairs in owned for i, _ in pairs) != sorted(ids):
        raise RuntimeError("the shard inputs do not own every record once")
    return files, written, owned


def place(units, shards, entry_of):
    t = time.time()
    shard_of, plan = assign_shards(units, shards)
    if set(shard_of) != set(entry_of):
        raise RuntimeError("the assignment does not place every record once")
    plan["seconds"] = round(time.time() - t, 1)
    return shard_of, plan


def manifest_of(artifact, ctx, families, fams, shards, files, written, owned, shard_of, plan,
                database, replacements, certify_s, t0) -> dict:
    index, ids = ctx.index, [e.event_id for e in ctx.universe]
    return {
        "schema": MANIFEST_SCHEMA,
        "artifact": {"id": Path(artifact).stem, "sha256": index.sha256, "bytes": index.size},
        "artifact_records": len(index.entries), "families": fams,
        "only": sorted(families) if families else None,
        "artifact_ids": ids, "shards": shards,
        "files": [f.name for f in files], "items": written,
        "files_sha256": [file_sha256(f) for f in files],
        "own": [len(pairs) for pairs in owned],
        "own_sha256": [own_digest(pairs) for pairs in owned],
        "assignment": {"basis": "record_key", "shard_of": [shard_of[i] for i in ids], **plan},
        "head": index.head, "bound_by": list(ctx.bound_by),
        "identity": identity(database, ctx.roster),
        "republished": len(replacements), "certify_s": certify_s,
        "seconds": round(time.time() - t0, 1), "peak_rss_mb": _rss()}


def prepare(artifact, roster_path, database, shards: int, out, *, families=None,
            index_path=None, legacy_artifact: bool = False) -> dict:
    """Certify ``artifact`` once and write the inputs of ``shards`` emission
    shards and their manifest into ``out``; returns the manifest."""
    if shards < 1:
        raise ValueError("the shard count must be at least 1")
    t0 = time.time()
    ctx = open_context(artifact, roster_path, families, index_path, legacy_artifact)
    roster, index, universe = ctx.roster, ctx.index, ctx.universe
    print("PREPARE index", len(index.entries), "records", index.size, "bytes",
          f"{time.time() - t0:.0f}s rss={_rss()}MB", flush=True)
    fams = sorted({e.family for e in universe})
    db = emit.load_database(fams, database)
    fam = db.kinetics.families
    certifier = emit.make_certifier(fam, emit.certification_rules(index.head, roster.polymer),
                                    emit.reference_thermo(db))
    entry_of = {e.event_id: e for e in universe}
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    # the manifest goes last: from here until it is rewritten, shards refuse the directory
    (out / MANIFEST).unlink(missing_ok=True)
    work = Path(tempfile.mkdtemp(dir=out, prefix=".prepare."))
    try:
        # 2. certify in artifact order; spill the certificates
        light, summaries = {}, []
        t = time.time()
        spill_path = work / "certificates.pkl.gz"
        with gzip.open(spill_path, "wb", compresslevel=1) as spill:
            certify_range(artifact, universe, entry_of, certifier, 0, len(universe),
                          lambda e, rc: pickle.dump(rc, spill, protocol=pickle.HIGHEST_PROTOCOL),
                          light, summaries)
        certify_s = round(time.time() - t, 1)
        # 3. the certification gate and the republished ids
        reasons, replacements, final_ids = gate(artifact, universe, entry_of, light, summaries)
        del summaries
        print("PREPARE gate", len(reasons), "republished", f"{certify_s}s certification",
              f"rss={_rss()}MB", flush=True)
        # 4. the assignment: one key's records on one shard
        t = time.time()
        systems, _refusals = emit.build_reference_systems(roster)
        units = assignment_units(artifact, universe, certifier, fam, systems)
        shard_of, plan = place(units, shards, entry_of)
        plan["seconds"] = round(time.time() - t, 1)
        print("PREPARE assignment", plan["units"], "units", len(plan["split"]), "split",
              len(plan["case_bound"]), "case-bound", f"{plan['seconds']}s rss={_rss()}MB",
              flush=True)
        # 5. publish each record into the inputs of the shards that read it
        with gzip.open(spill_path, "rb") as spill:
            files, written, owned = publish(
                artifact, universe, entry_of, (pickle.load(spill) for _ in universe), reasons,
                replacements, final_ids, shard_of, shards, work)
        spill_path.unlink()
        manifest = manifest_of(artifact, ctx, families, fams, shards, files, written, owned,
                               shard_of, plan, database, replacements, certify_s, t0)
        for f in files:
            os.replace(f, out / f.name)
        tmp = work / MANIFEST
        tmp.write_text(json.dumps(manifest, sort_keys=True, separators=(",", ":")) + "\n")
        os.replace(tmp, out / MANIFEST)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print("PREPARE done", shards, "shards", written, "items",
          f"{manifest['seconds']}s peak rss={manifest['peak_rss_mb']}MB", flush=True)
    return manifest


# ---------------------------------------------------------------------------
# prepare in blocks (I-304b): plan, J independent blocks, one finishing step
# ---------------------------------------------------------------------------
PLAN = "plan.json"
PLAN_SCHEMA = "kmc-emit-prepare-plan/1"
BLOCK_SCHEMA = "kmc-emit-prepare-block/1"
SEAL_SCHEMA = "kmc-emit-prepare-seal/1"
RESOURCE_EXIT = 3


class BlockRefusal(ValueError):
    """A block run, its plan or its blocks that cannot be published."""

    def __init__(self, reason):
        super().__init__(f"refused:blocks: {reason}")
        self.reason = reason


def canonical_bytes(value) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def ids_sha256(ids) -> str:
    return hashlib.sha256(json.dumps(list(ids), separators=(",", ":")).encode()).hexdigest()


def block_name(block: int, blocks: int) -> str:
    return f"b{block}of{blocks}"


def block_ranges(records: int, blocks: int) -> list:
    """(lo, hi) of each of ``blocks`` contiguous ranges of ``records``
    positions, as equal as integers allow."""
    return [(records * j // blocks, records * (j + 1) // blocks) for j in range(blocks)]


def _resolved(path):
    return None if path is None else str(Path(path).resolve())


def make_plan(artifact, roster_path, database, blocks: int, *, families=None, index_path=None,
              legacy_artifact: bool = False):
    """(the plan of a block prepare, its ``Context``): the inputs (paths and
    digests), the identity and each block's range of artifact order with the
    digest of its ids.  A pure function of the inputs' bytes and the code, so
    every stage recomputes it and compares it with the plan file."""
    ctx = open_context(artifact, roster_path, families, index_path, legacy_artifact)
    universe = ctx.universe
    if not universe:
        raise BlockRefusal("no records to prepare")
    if not 1 <= blocks <= len(universe):
        raise BlockRefusal(f"{blocks} blocks for {len(universe)} records: at least 1, at most "
                           "one per record")
    ids = [e.event_id for e in universe]
    ranges = [{"block": j, "lo": lo, "hi": hi, "first": ids[lo], "last": ids[hi - 1],
               "ids_sha256": ids_sha256(ids[lo:hi])}
              for j, (lo, hi) in enumerate(block_ranges(len(ids), blocks))]
    plan = {"schema": PLAN_SCHEMA,
            "artifact": {"path": _resolved(artifact), "id": Path(artifact).stem,
                         "sha256": ctx.index.sha256, "bytes": ctx.index.size},
            "index": _resolved(index_path),
            "roster": {"path": _resolved(roster_path), "sha256": ctx.roster.sha256},
            "database": _resolved(database), "only": sorted(families) if families else None,
            "legacy_artifact": bool(legacy_artifact), "bound_by": list(ctx.bound_by),
            "identity": identity(database, ctx.roster),
            "artifact_records": len(ctx.index.entries), "records": len(ids),
            "families": sorted({e.family for e in universe}), "ids_sha256": ids_sha256(ids),
            "species_reset": SPECIES_RESET, "blocks": blocks, "ranges": ranges}
    return plan, ctx


def _exclusive_write(path: Path, data: bytes) -> bool:
    """Write ``data`` at ``path`` complete or not at all, only if nothing is
    there: a temporary file linked into place (``link`` fails on an existing
    name, atomically).  False when ``path`` already exists."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(tmp, path)
        except FileExistsError:
            return False
        return True
    finally:
        os.unlink(tmp)


def write_plan(run, artifact, roster_path, database, blocks: int, *, families=None,
               index_path=None, legacy_artifact: bool = False) -> dict:
    """Write ``run``/plan.json once (a rerun with the same plan is a no-op;
    with another plan it is refused); returns the plan."""
    plan, _ctx = make_plan(artifact, roster_path, database, blocks, families=families,
                           index_path=index_path, legacy_artifact=legacy_artifact)
    path = Path(run) / PLAN
    data = canonical_bytes(plan)
    if not _exclusive_write(path, data) and path.read_bytes() != data:
        raise BlockRefusal(f"{path} holds another plan (other inputs, code or blocks): use "
                           "another run directory")
    print("PLAN", path, plan["records"], "records", blocks, "blocks",
          hashlib.sha256(data).hexdigest(), flush=True)
    return plan


def load_plan(run):
    """(plan, its sha256, ``Context``): ``run``'s plan, refused unless the
    plan recomputed from its inputs' current bytes, the database and this
    code is the file's."""
    path = Path(run) / PLAN
    try:
        data = path.read_bytes()
        plan = json.loads(data)
    except (OSError, ValueError) as error:
        raise BlockRefusal(f"{path}: {error}")
    if not isinstance(plan, dict) or plan.get("schema") != PLAN_SCHEMA:
        raise BlockRefusal(f"{path} is not a {PLAN_SCHEMA} plan")
    now, ctx = make_plan(plan["artifact"]["path"], plan["roster"]["path"], plan["database"],
                         plan["blocks"], families=plan["only"], index_path=plan["index"],
                         legacy_artifact=plan["legacy_artifact"])
    if canonical_bytes(now) != data:
        differ = sorted(k for k in set(now) | set(plan) if now.get(k) != plan.get(k))
        raise BlockRefusal(f"{path} is not the plan of its inputs as they are now (differs in "
                           f"{differ})")
    return plan, hashlib.sha256(data).hexdigest(), ctx


def _write_json(path: Path, value) -> None:
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_bytes(canonical_bytes(value))
    os.replace(tmp, path)


def _attempt_id() -> str:
    import socket
    return f"{socket.gethostname()}.{os.getpid()}.{time.time_ns()}"


class _Framed:
    """A gzip pickle stream written at ``path`` as ``.part``, synced and
    renamed on a clean exit (``info``: its file name, sha256 and size); an
    exception leaves only the ``.part``."""

    def __init__(self, path: Path):
        self.path, self.part, self.info = path, path.with_name(path.name + ".part"), None

    def __enter__(self):
        self.raw = open(self.part, "wb")
        self.handle = gzip.GzipFile(fileobj=self.raw, mode="wb", compresslevel=1)
        return self

    def dump(self, item) -> None:
        pickle.dump(item, self.handle, protocol=pickle.HIGHEST_PROTOCOL)

    def __exit__(self, kind, value, trace):
        self.handle.close()
        self.raw.flush()
        os.fsync(self.raw.fileno())
        self.raw.close()
        if kind is None:
            os.replace(self.part, self.path)
            self.info = {"file": self.path.name, "sha256": file_sha256(self.path),
                         "bytes": self.path.stat().st_size}
        return False


def prepare_block(run, block: int) -> dict:
    """Certify block ``block`` of ``run``'s plan in a fresh certifier and
    compute its records' reference-case keys; write the attempt's files and
    seal them, unless the block is sealed already.  Returns the attempt's
    status (``state``: sealed, superseded or already-sealed).  A resource
    error marks the attempt failed and propagates: it never seals."""
    t0 = time.time()
    run = Path(run)
    plan, plan_sha, ctx = load_plan(run)
    blocks = plan["blocks"]
    if not 0 <= block < blocks:
        raise BlockRefusal(f"block {block} of a {blocks}-block plan")
    name = block_name(block, blocks)
    seal = run / "seals" / f"{name}.json"
    if seal.exists():
        print("BLOCK", name, "already sealed:", json.loads(seal.read_text())["attempt"],
              flush=True)
        return {"state": "already-sealed", "block": block}
    attempt = _attempt_id()
    attempts = run / "attempts"
    attempts.mkdir(parents=True, exist_ok=True)
    status = attempts / f"{name}.{attempt}.status.json"
    base = {"block": block, "blocks": blocks, "attempt": attempt, "plan_sha256": plan_sha,
            "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    _write_json(status, {**base, "state": "running"})
    try:
        state = _run_block(run, plan, plan_sha, ctx, block, name, attempt, seal, t0)
    except BaseException as error:
        from rmgpy.kmc.resource_errors import RESOURCE_ERRORS
        _write_json(status, {**base, "state": "failed",
                             "resource": isinstance(error, RESOURCE_ERRORS),
                             "error": f"{type(error).__name__}: {error}"[:2000]})
        raise
    _write_json(status, {**base, **state})
    print("BLOCK", name, state["state"], attempt, f"{state['records']} records",
          f"{round(time.time() - t0, 1)}s peak rss={_rss()}MB", flush=True)
    return state


def _run_block(run, plan, plan_sha, ctx, block, name, attempt, seal, t0) -> dict:
    rng = plan["ranges"][block]
    lo, hi = rng["lo"], rng["hi"]
    universe, roster, head = ctx.universe, ctx.roster, ctx.index.head
    entry_of = {e.event_id: e for e in universe}
    db = emit.load_database(plan["families"], plan["database"])
    fam = db.kinetics.families
    certifier = emit.make_certifier(fam, emit.certification_rules(head, roster.polymer),
                                    emit.reference_thermo(db))
    stem = run / "attempts" / f"{name}.{attempt}"
    # 1. certify universe[lo:hi] in artifact order (prepare's step 2)
    light, summaries, sources = {}, [], []
    t = time.time()

    with _Framed(Path(f"{stem}.certs.pkl.gz")) as certs:
        def spill(e, rc):
            sources.append(e.event_id)
            certs.dump({"source": e.event_id, "certificate": rc})
        certify_range(plan["artifact"]["path"], universe, entry_of, certifier, lo, hi, spill,
                      light, summaries, label=f"{name} certified")
        certs.dump({"end": len(sources)})
    certify_s = round(time.time() - t, 1)
    # 2. the keys of its records (prepare's step 4, without the placement)
    t = time.time()
    systems, _refusals = emit.build_reference_systems(roster)
    shapes = subset_shapes(fam, systems) if systems else Counter()
    keyed, subsets = [], {}
    for e, key, mols, others in record_keys(plan["artifact"]["path"], universe, certifier,
                                            lo, hi):
        keyed.append(key)
        if key is not None and key not in subsets:
            subsets[key] = key_subsets(key, mols, others, shapes)
    keys_s = round(time.time() - t, 1)
    if sources != [e.event_id for e in universe[lo:hi]] or len(keyed) != hi - lo:
        raise RuntimeError(f"block {name} did not certify and key exactly its range")
    header = {"schema": BLOCK_SCHEMA, "plan_sha256": plan_sha, "block": block,
              "blocks": plan["blocks"], "lo": lo, "hi": hi, "ids_sha256": rng["ids_sha256"],
              "attempt": attempt}
    items = [header] + [{"source": i, "light": light[i], "summary": s, "key": k}
                        for i, s, k in zip(sources, summaries, keyed)] + \
        [{"subsets": subsets}, {"end": len(sources), "ids_sha256": ids_sha256(sources)}]
    with _Framed(Path(f"{stem}.meta.pkl.gz")) as meta:
        for item in items:
            meta.dump(item)
    for f in (meta.path, certs.path):
        os.chmod(f, 0o444)
    sealed = {"schema": SEAL_SCHEMA, "plan_sha256": plan_sha, "block": block,
              "blocks": plan["blocks"], "attempt": attempt, "records": len(sources),
              "ids_sha256": ids_sha256(sources), "meta": meta.info, "certs": certs.info,
              "certify_s": certify_s, "keys_s": keys_s,
              "seconds": round(time.time() - t0, 1), "peak_rss_mb": _rss()}
    state = "sealed" if _exclusive_write(seal, canonical_bytes(sealed)) else "superseded"
    return {"state": state, "records": len(sources), "certify_s": certify_s, "keys_s": keys_s,
            "meta": meta.info, "certs": certs.info}


def _stream(source):
    """The items of a gzip pickle stream (a path, or an open binary file read
    from its start), then None if it ends unframed."""
    if not isinstance(source, (str, Path)):
        source.seek(0)
    with gzip.open(source, "rb") as handle:
        while True:
            try:
                yield pickle.load(handle)
            except (EOFError, pickle.UnpicklingError, OSError, zlib.error):
                yield None
                return


def _attempt_states(run: Path, name: str) -> list:
    out = []
    for path in sorted((run / "attempts").glob(f"{name}.*.status.json")):
        try:
            s = json.loads(path.read_text())
        except (OSError, ValueError):
            out.append(f"{path.name}: unreadable status")
            continue
        state = s.get("state")
        if state == "failed":
            state = ("failed: resource " if s.get("resource") else "failed: ") + \
                str(s.get("error"))[:200]
        elif state == "running":
            state = "running, or killed (no seal)"
        elif state == "sealed":
            state = "sealed, but its seal is gone"
        out.append(f"attempt {s.get('attempt')}: {state}")
    return out or ["no attempt"]


@dataclass
class Blocks:
    """What the finishing step reads of the sealed blocks: per record (in
    artifact order) its light row, summary and key; each key's subsets; each
    block's certificate stream."""
    light: dict
    summaries: list
    keyed: list
    subsets: dict
    certs: list
    seals: list


def read_blocks(run, plan, plan_sha, universe) -> Blocks:
    """Every block of ``plan``, refused by name unless each has exactly one
    seal of this plan whose files are the sealed bytes and hold exactly its
    range of artifact order (section 2 of the design)."""
    run = Path(run)
    blocks = plan["blocks"]
    names = [block_name(j, blocks) for j in range(blocks)]
    seals_dir = run / "seals"
    present = {p.name[:-len(".json")] for p in seals_dir.glob("*.json")} \
        if seals_dir.is_dir() else set()
    problems = [f"foreign seal {n}.json (not a block of this {blocks}-block plan)"
                for n in sorted(present - set(names))]
    problems += [f"block {n} is missing: no sealed attempt ({'; '.join(_attempt_states(run, n))})"
                 for n in names if n not in present]
    if problems:
        raise BlockRefusal("; ".join(problems))
    ids = [e.event_id for e in universe]
    position = {i: k for k, i in enumerate(ids)}
    ranges = plan["ranges"]

    def owner_of(i):
        k = position[i]
        return next(r["block"] for r in ranges if r["lo"] <= k < r["hi"])

    out = Blocks({}, [], [], {}, [], [])
    try:
        _read_each(run, plan, plan_sha, names, ids, position, owner_of, out)
    except BaseException:
        for _n, handle, _expected in out.certs:
            handle.close()
        raise
    return out


def _read_each(run, plan, plan_sha, names, ids, position, owner_of, out) -> None:
    blocks, ranges, seals_dir = plan["blocks"], plan["ranges"], run / "seals"
    for j, n in enumerate(names):
        rng = ranges[j]
        try:
            seal = json.loads((seals_dir / f"{n}.json").read_text())
        except (OSError, ValueError) as error:
            raise BlockRefusal(f"block {n}: its seal is unreadable ({error})")
        if (seal.get("schema"), seal.get("plan_sha256"), seal.get("block"),
                seal.get("blocks")) != (SEAL_SCHEMA, plan_sha, j, blocks):
            raise BlockRefusal(f"block {n}: foreign seal (another plan, block or schema)")
        # each file is opened once: what is hashed is what is unpickled (a
        # file replaced after the check is not read)
        opened = {}
        for part in ("meta", "certs"):
            f = run / "attempts" / str(seal[part]["file"])
            handle = open(f, "rb") if f.parent == run / "attempts" and f.is_file() else None
            if handle is not None:
                h = hashlib.sha256()
                for chunk in iter(lambda: handle.read(1 << 20), b""):
                    h.update(chunk)
            if handle is None or h.hexdigest() != seal[part]["sha256"]:
                if handle is not None:
                    handle.close()
                for other in opened.values():
                    other.close()
                raise BlockRefusal(f"block {n}: its {part} file {f.name} is missing or is not "
                                   "the sealed bytes (altered, truncated or replaced)")
            opened[part] = handle
        expected = ids[rng["lo"]:rng["hi"]]
        out.certs.append((n, opened["certs"], expected))   # closed by the caller on refusal
        with opened["meta"] as handle:
            handle.seek(0)
            items = _stream(io.BytesIO(handle.read()))
        header = next(items)
        want = {"schema": BLOCK_SCHEMA, "plan_sha256": plan_sha, "block": j, "blocks": blocks,
                "lo": rng["lo"], "hi": rng["hi"], "ids_sha256": rng["ids_sha256"],
                "attempt": seal["attempt"]}
        if not isinstance(header, dict) or {k: header.get(k) for k in want} != want:
            raise BlockRefusal(f"block {n}: foreign block (its header is not this plan's block "
                               f"{j} [{rng['lo']}, {rng['hi']}) of attempt {seal['attempt']})")
        got, trailer, subsets = [], None, None
        for item in items:
            if item is None:
                break
            if "end" in item:
                trailer = item
                break
            if "subsets" in item:
                subsets = item["subsets"]
                continue
            got.append(item["source"])
            out.light[item["source"]] = item["light"]
            out.summaries.append(item["summary"])
            out.keyed.append((item["source"], item["key"]))
        if got != expected:
            seen, why = set(), []
            for i in got:
                if i in seen:
                    why.append(f"duplicate: holds {i} twice")
                elif i not in position:
                    why.append(f"foreign: holds {i}, not a record of the plan")
                elif not rng["lo"] <= position[i] < rng["hi"]:
                    why.append(f"duplicate: holds {i}, a record of block "
                               f"{block_name(owner_of(i), blocks)}")
                seen.add(i)
            lost = [i for i in expected if i not in seen]
            if lost:
                why.append(f"short: {len(lost)} of its {len(expected)} records missing "
                           f"(first {lost[0]})")
            raise BlockRefusal(f"block {n}: " + ("; ".join(why) or "its records are out of "
                                                  "artifact order"))
        if trailer is None or subsets is None or trailer.get("end") != len(got) or \
                trailer.get("ids_sha256") != rng["ids_sha256"]:
            raise BlockRefusal(f"block {n}: truncated (no trailer, or one that does not count "
                               "its records)")
        for key, count in subsets.items():
            if out.subsets.setdefault(key, count) != count:
                raise BlockRefusal(f"block {n}: counts {count} roster subsets for the key "
                                   f"{key[:120]}, another block {out.subsets[key]}")
        out.seals.append({"block": j, "attempt": seal["attempt"],
                          "meta_sha256": seal["meta"]["sha256"],
                          "certs_sha256": seal["certs"]["sha256"],
                          "certify_s": seal.get("certify_s"), "keys_s": seal.get("keys_s")})
    missing_keys = {k for _i, k in out.keyed if k is not None} - set(out.subsets)
    if missing_keys:
        raise BlockRefusal(f"{len(missing_keys)} keys without a roster subset count")


def _block_certificates(certs):
    """Each block's certificates in artifact order, refused (before anything
    is published) unless each stream holds its block's records in order."""
    for n, handle, expected in certs:
        items = _stream(handle)
        for i in expected:
            item = next(items)
            if item is None or item.get("source") != i:
                raise BlockRefusal(f"block {n}: its certificates are not its records' in order "
                                   f"(at {i})")
            yield item["certificate"]
        end = next(items)
        if end is None or end.get("end") != len(expected):
            raise BlockRefusal(f"block {n}: its certificate stream is truncated or too long")
        handle.close()


def prepare_finish(run, shards: int):
    """Publish ``run``'s sealed blocks as the inputs of ``shards`` emission
    shards: the certification gate over every record, the assignment of the
    blocks' keys, every record republished, relinked and written to the shards
    that read it, and the manifest, all in a new generation directory
    ``run/prepared/<generation>`` (never written again), then ``run/current``
    pointed at it.  Returns (manifest, generation directory)."""
    if shards < 1:
        raise ValueError("the shard count must be at least 1")
    t0 = time.time()
    run = Path(run).resolve()
    plan, plan_sha, ctx = load_plan(run)
    blocks = read_blocks(run, plan, plan_sha, ctx.universe)
    try:
        return _finish(run, plan, plan_sha, ctx, blocks, shards, t0)
    finally:
        for _n, handle, _expected in blocks.certs:
            handle.close()


def _finish(run, plan, plan_sha, ctx, blocks, shards, t0):
    universe = ctx.universe
    artifact = plan["artifact"]["path"]
    entry_of = {e.event_id: e for e in universe}
    print("FINISH blocks", plan["blocks"], "records", len(universe),
          f"{time.time() - t0:.0f}s rss={_rss()}MB", flush=True)
    reasons, replacements, final_ids = gate(artifact, universe, entry_of, blocks.light,
                                            blocks.summaries)
    blocks.summaries = None
    print("FINISH gate", len(reasons), "republished", f"rss={_rss()}MB", flush=True)
    shard_of, assignment = place(units_of(blocks.keyed, blocks.subsets), shards, entry_of)
    blocks.keyed = None
    print("FINISH assignment", assignment["units"], "units", len(assignment["split"]), "split",
          len(assignment["case_bound"]), "case-bound", f"rss={_rss()}MB", flush=True)
    generations = run / "prepared"
    generations.mkdir(exist_ok=True)
    generation = f"{shards}shards.{time.strftime('%Y%m%dT%H%M%SZ', time.gmtime())}." \
                 f"{os.getpid()}.{time.time_ns() % 10**9}"
    work = generations / f".tmp.{generation}"
    work.mkdir()
    try:
        files, written, owned = publish(artifact, universe, entry_of,
                                        _block_certificates(blocks.certs), reasons,
                                        replacements, final_ids, shard_of, shards, work)
        certify_s = round(sum(s["certify_s"] or 0 for s in blocks.seals), 1)
        manifest = manifest_of(artifact, ctx, plan["only"], plan["families"], shards, files,
                               written, owned, shard_of, assignment, plan["database"],
                               replacements, certify_s, t0)
        manifest["blocks"] = {"plan_sha256": plan_sha, "blocks": plan["blocks"],
                              "seals": blocks.seals}
        (work / MANIFEST).write_bytes(canonical_bytes(manifest))
        for f in list(files) + [work / MANIFEST]:
            os.chmod(f, 0o444)
        final = generations / generation
        os.rename(work, final)
    except BaseException:
        shutil.rmtree(work, ignore_errors=True)
        raise
    link = run / f".current.{os.getpid()}"
    link.unlink(missing_ok=True)
    os.symlink(Path("prepared") / generation, link)
    os.replace(link, run / "current")
    print("FINISH done", shards, "shards", written, "items", final,
          f"{manifest['seconds']}s peak rss={manifest['peak_rss_mb']}MB", flush=True)
    print("PREPARED", final, flush=True)
    return manifest, final


# ---------------------------------------------------------------------------
# two prepared directories compared (serial against blocks)
# ---------------------------------------------------------------------------
# manifest fields that are a run's own (timing, gzip bytes, how it was made;
# and ``assignment.seconds``)
RUN_OWN = ("seconds", "certify_s", "peak_rss_mb", "files_sha256", "blocks")


def _paths(a, b, path=""):
    if type(a) is not type(b):
        yield path
    elif isinstance(a, dict):
        for k in sorted(set(a) | set(b), key=str):
            if a.get(k, _paths) is _paths or b.get(k, _paths) is _paths:
                yield f"{path}.{k}"
            else:
                yield from _paths(a[k], b[k], f"{path}.{k}")
    elif isinstance(a, (list, tuple)):
        if len(a) != len(b):
            yield f"{path}[len]"
        for k, (x, y) in enumerate(zip(a, b)):
            yield from _paths(x, y, f"{path}[{k}]")
    elif a != b:
        yield path


def _items(path: Path) -> list:
    out = []
    for item in _stream(path):
        if item is None or "end" in item:
            break
        out.append(item)
    return out


def compare_prepared(a, b) -> dict:
    """Two prepared directories (manifest and shard inputs) compared under
    the block identity criterion: every manifest field but a run's own
    (``RUN_OWN``) equal; every shard input equal item by item, except that a
    certificate may differ in its tree nodes' atom order
    (``certificates._same_certificate``), and the record's reference digest
    with it (each side's must be its own certificate's).  Returns
    {"differences": [path, ...], "node_order": [source id, ...],
    "items": count}."""
    a, b = Path(a), Path(b)
    ma, mb = (json.loads((d / MANIFEST).read_text()) for d in (a, b))

    def comparable(m):
        out = {k: v for k, v in m.items() if k not in RUN_OWN}
        out["assignment"] = {k: v for k, v in m["assignment"].items() if k != "seconds"}
        return out

    diffs = [f"manifest{p}" for p in _paths(comparable(ma), comparable(mb))]
    node_order, count = [], 0
    for k, (fa, fb) in enumerate(zip(ma["files"], mb["files"])):
        ia, ib = _items(a / fa), _items(b / fb)
        if len(ia) != len(ib):
            diffs.append(f"shard{k}[len] {len(ia)} != {len(ib)}")
        for n, (x, y) in enumerate(zip(ia, ib)):
            count += 1
            where = f"shard{k}[{n}]({x.get('source', '')[:16]})"
            ca, cb = x.get("certificate"), y.get("certificate")
            ra, rb = dict(x.get("record") or {}), dict(y.get("record") or {})
            refa, refb = dict(ra.pop("rate_certificate", None) or {}), \
                dict(rb.pop("rate_certificate", None) or {})
            sha_a, sha_b = refa.pop("sha256", None), refb.pop("sha256", None)
            rest_a = {key: v for key, v in x.items() if key not in ("record", "certificate")}
            rest_b = {key: v for key, v in y.items() if key not in ("record", "certificate")}
            diffs += [f"{where}{p}" for p in _paths(rest_a, rest_b)]
            diffs += [f"{where}.record{p}" for p in _paths(ra, rb)]
            diffs += [f"{where}.record.rate_certificate{p}" for p in _paths(refa, refb)]
            if (ca is None) != (cb is None):
                diffs.append(f"{where}.certificate (present on one side only)")
                continue
            if ca is None:
                if sha_a != sha_b:
                    diffs.append(f"{where}.record.rate_certificate.sha256 (no certificate)")
                continue
            for side, c, sha in (("a", ca, sha_a), ("b", cb, sha_b)):
                if cm.certificate_sha256(c) != sha:
                    diffs.append(f"{where}.record.rate_certificate.sha256 is not side {side}'s "
                                 "certificate's")
            if cm.canonical_json(ca) == cm.canonical_json(cb):
                continue
            outer_a = {key: v for key, v in ca.items() if key != "certificate"}
            outer_b = {key: v for key, v in cb.items() if key != "certificate"}
            if cm.canonical_json(outer_a) == cm.canonical_json(outer_b) and \
                    cm._same_certificate(ca.get("certificate"), cb.get("certificate")):
                node_order.append(x.get("source"))
            else:
                diffs += [f"{where}.certificate{p}" for p in _paths(ca, cb)]
    return {"differences": diffs, "node_order": sorted(set(node_order)), "items": count}


@dataclass
class PreparedShard:
    """One shard's prepared input (``Emitter.init``'s selection)."""
    head: dict
    sha256: str
    artifact_records: int
    families: list
    artifact_ids: list
    records: list
    store: dict
    replacements: dict
    picked: list


def open_manifest(prepared, artifact, roster, database, chunk, only=None) -> dict:
    """The manifest of ``prepared``, refused unless it was prepared from these
    artifact bytes, roster, database content and code, for this shard count
    and family selection."""
    path = Path(prepared) / MANIFEST
    if not path.is_file():
        raise PreparedRefusal("refused:prepared", f"{prepared} has no {MANIFEST} (an "
                                                  "incomplete or failed prepare)")
    manifest = json.loads(path.read_text())
    if manifest.get("schema") in OLD_MANIFEST_SCHEMAS:
        raise PreparedRefusal("refused:prepared", f"{path} is a {manifest['schema']} manifest "
                                                  "(the sorted(ids)[k::n] assignment, no shard "
                                                  f"input digests); this code reads "
                                                  f"{MANIFEST_SCHEMA}: prepare again")
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise PreparedRefusal("refused:prepared", f"{path} is not a {MANIFEST_SCHEMA} manifest")
    k, n = emit.parse_chunk(chunk)
    if n != manifest["shards"]:
        raise PreparedRefusal("refused:prepared", f"prepared for {manifest['shards']} shards, "
                                                  f"asked for shard {k}/{n}")
    if (sorted(only) if only else None) != manifest["only"]:
        raise PreparedRefusal("refused:prepared", f"prepared for families {manifest['only']}, "
                                                  f"asked for {only}")
    sha, size = astream.artifact_digest(artifact)
    if (sha, size) != (manifest["artifact"]["sha256"], manifest["artifact"]["bytes"]):
        raise PreparedRefusal("refused:prepared", f"prepared from another artifact "
                                                  f"({manifest['artifact']['sha256'][:16]}...) "
                                                  f"than {artifact} ({sha[:16]}...)")
    now = identity(database, roster)
    differ = [key for key, value in now.items() if manifest["identity"].get(key) != value]
    if differ:
        raise PreparedRefusal("refused:prepared", f"prepared under another {', '.join(differ)} "
                                                  "than this shard's")
    return manifest


def assigned_ids(manifest, k) -> list:
    """The artifact ids the manifest assigns to shard ``k``, sorted."""
    ids, shard_of = manifest["artifact_ids"], manifest["assignment"]["shard_of"]
    if len(shard_of) != len(ids) or any(not 0 <= s < manifest["shards"] for s in shard_of):
        raise PreparedRefusal("refused:prepared", "the manifest's assignment does not place "
                                                  "every artifact id on one of its shards")
    return sorted(i for i, s in zip(ids, shard_of) if s == k)


def read_shard(prepared, manifest, chunk) -> PreparedShard:
    """Shard ``chunk``'s input: its published records (with the records they
    link), their certificate store, their republished ids and its records in
    source-id order.  The input's bytes are read once and checked against the
    manifest's sha256 before anything is unpickled."""
    k, n = emit.parse_chunk(chunk)
    records, store, replacements, own = [], {}, {}, {}
    count = None
    name = manifest["files"][k]
    try:
        data = (Path(prepared) / name).read_bytes()
    except OSError as error:
        raise PreparedRefusal("refused:prepared", f"shard input {name}: {error}")
    if hashlib.sha256(data).hexdigest() != manifest["files_sha256"][k]:
        raise PreparedRefusal("refused:prepared", f"shard input {name}: its sha256 is not the "
                                                  "manifest's (another prepare's, altered or "
                                                  "truncated)")
    with gzip.GzipFile(fileobj=io.BytesIO(data), mode="rb") as handle:
        while True:
            try:
                item = pickle.load(handle)
            except (EOFError, pickle.UnpicklingError, OSError, zlib.error):
                break
            if "end" in item:
                count = item["end"]
                break
            records.append(item["record"])
            if item["own"]:
                record = item["record"]
                if item["source"] in own:
                    raise PreparedRefusal("refused:prepared", f"shard input {name} owns "
                                                              f"{item['source']} twice")
                own[item["source"]] = record
                store[record["rate_certificate"]["sha256"]] = item["certificate"]
                if record["event_id"] != item["source"]:
                    replacements[item["source"]] = record["event_id"]
    del data
    if count is None or count != len(records) or count != manifest["items"][k]:
        raise PreparedRefusal("refused:prepared", f"shard input {name} is truncated")
    wanted = assigned_ids(manifest, k)
    if sorted(own) != wanted:
        raise PreparedRefusal("refused:prepared", f"shard input {name} does not hold the "
                                                  f"records the manifest assigns to shard {k}")
    if own_digest([s, r["event_id"]] for s, r in own.items()) != manifest["own_sha256"][k]:
        raise PreparedRefusal("refused:prepared", f"shard input {name}: its records' source and "
                                                  "published ids are not the manifest's")
    return PreparedShard(head=manifest["head"], sha256=manifest["artifact"]["sha256"],
                         artifact_records=manifest["artifact_records"],
                         families=manifest["families"], artifact_ids=manifest["artifact_ids"],
                         records=records, store=store, replacements=replacements,
                         picked=[own[i] for i in wanted])


def main(argv=None) -> int:
    import argparse
    parser = argparse.ArgumentParser(prog="python -m rmgpy.kmc.emit prepare")
    parser.add_argument("artifact")
    parser.add_argument("roster")
    parser.add_argument("--database", required=True)
    parser.add_argument("--shards", type=int, required=True)
    parser.add_argument("--out", required=True, help="the prepared directory")
    parser.add_argument("--family", action="append", default=None)
    parser.add_argument("--index", default=None, help="the artifact's index (emit index)")
    parser.add_argument(emit.LEGACY_ARTIFACT_FLAG, action="store_true")
    args = parser.parse_args(argv)
    try:
        prepare(args.artifact, args.roster, args.database, args.shards, args.out,
                families=args.family, index_path=args.index,
                legacy_artifact=args.legacy_artifact)
    except (emit.ArtifactBindingRefusal, emit.RosterRefusal, astream.IndexMismatch) as error:
        print("REFUSED", error, flush=True)
        return 1
    return 0


def main_plan(argv=None) -> int:
    import argparse
    parser = argparse.ArgumentParser(prog="python -m rmgpy.kmc.emit prepare-plan")
    parser.add_argument("artifact")
    parser.add_argument("roster")
    parser.add_argument("--database", required=True)
    parser.add_argument("--blocks", type=int, required=True)
    parser.add_argument("--run", required=True, help="the block prepare's run directory")
    parser.add_argument("--family", action="append", default=None)
    parser.add_argument("--index", default=None, help="the artifact's index (emit index)")
    parser.add_argument(emit.LEGACY_ARTIFACT_FLAG, action="store_true")
    args = parser.parse_args(argv)
    try:
        write_plan(args.run, args.artifact, args.roster, args.database, args.blocks,
                   families=args.family, index_path=args.index,
                   legacy_artifact=args.legacy_artifact)
    except (BlockRefusal, emit.ArtifactBindingRefusal, emit.RosterRefusal,
            astream.IndexMismatch) as error:
        print("REFUSED", error, flush=True)
        return 1
    return 0


def main_block(argv=None) -> int:
    """Exit 0: sealed, superseded or already sealed; 1: refused; 3: a
    resource error (the attempt is marked failed and never seals)."""
    import argparse
    from rmgpy.kmc.resource_errors import RESOURCE_ERRORS
    parser = argparse.ArgumentParser(prog="python -m rmgpy.kmc.emit prepare-block")
    parser.add_argument("run")
    parser.add_argument("--block", type=int, required=True)
    args = parser.parse_args(argv)
    try:
        prepare_block(args.run, args.block)
    except (BlockRefusal, emit.ArtifactBindingRefusal, emit.RosterRefusal,
            astream.IndexMismatch) as error:
        print("REFUSED", error, flush=True)
        return 1
    except RESOURCE_ERRORS as error:
        print("BLOCK FAILED resource", type(error).__name__, flush=True)
        return RESOURCE_EXIT
    return 0


def main_finish(argv=None) -> int:
    import argparse
    parser = argparse.ArgumentParser(prog="python -m rmgpy.kmc.emit prepare-finish")
    parser.add_argument("run")
    parser.add_argument("--shards", type=int, required=True)
    args = parser.parse_args(argv)
    try:
        prepare_finish(args.run, args.shards)
    except (BlockRefusal, emit.ArtifactBindingRefusal, emit.RosterRefusal,
            astream.IndexMismatch) as error:
        print("REFUSED", error, flush=True)
        return 1
    return 0
