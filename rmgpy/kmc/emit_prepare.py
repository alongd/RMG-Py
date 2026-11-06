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


def assignment_units(artifact, universe, certifier, fam, systems) -> list:
    """The units of ``universe`` (``Unit``): records grouped by
    ``emit.record_key``, each key's subsets counted from the roster."""
    from rmgpy.kmc import channels as chm
    shapes = subset_shapes(fam, systems) if systems else Counter()
    by_key, subsets, units = defaultdict(list), {}, []
    species = chm.SpeciesKeys()
    for count, e in enumerate(universe, 1):
        if count % SPECIES_RESET == 0:
            species = chm.SpeciesKeys()
        d = _one(artifact, e).to_dict()
        try:
            key, mols, others = emit.record_parts(certifier, d, species)
        except Exception as error:   # its own unit; the shard's row reports it
            emit.reraise_resource(error)
            units.append(Unit(f"unkeyed:{e.event_id}", [e.event_id]))
            continue
        name = json.dumps(key, separators=(",", ":"))
        by_key[name].append(e.event_id)
        if name not in subsets:
            subsets[name] = shapes[(key[0], "forward", emit.component_formulas(mols))] + \
                shapes[(key[0], "reverse", emit.component_formulas(others))]
    units += [Unit(name, ids, subsets[name]) for name, ids in by_key.items()]
    return units


def own_digest(pairs) -> str:
    """sha256 of a shard's sorted (source id, published id) pairs."""
    return hashlib.sha256(json.dumps(sorted(pairs), separators=(",", ":")).encode()).hexdigest()


def file_sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def prepare(artifact, roster_path, database, shards: int, out, *, families=None,
            index_path=None, legacy_artifact: bool = False) -> dict:
    """Certify ``artifact`` once and write the inputs of ``shards`` emission
    shards and their manifest into ``out``; returns the manifest."""
    if shards < 1:
        raise ValueError("the shard count must be at least 1")
    t0 = time.time()
    roster = emit.load_roster(roster_path)
    index = astream.index_for(artifact, index_path)
    head = index.head
    bound_by = emit.bind_artifact(head, roster.digests(), f"the roster {roster_path}",
                                  legacy_artifact)
    emit.check_roster_for_artifact(roster, head)
    print("PREPARE index", len(index.entries), "records", index.size, "bytes",
          f"{time.time() - t0:.0f}s rss={_rss()}MB", flush=True)
    universe = [e for e in index.entries if not families or e.family in families]
    fams = sorted({e.family for e in universe})
    db = emit.load_database(fams, database)
    fam = db.kinetics.families
    certifier = emit.make_certifier(fam, emit.certification_rules(head, roster.polymer),
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
            for count, e in enumerate(universe, 1):
                record = _one(artifact, e)
                d = record.to_dict()
                light[e.event_id] = _light(d)
                owner = None
                if (d.get("rate_source") or {}).get("kind") == "reference-thermo reverse":
                    partner = entry_of.get(d.get("reverse_of"))
                    owner = None if partner is None else _one(artifact, partner)
                rc = cm.certify_one(certifier, record, d, owner)
                pickle.dump(rc, spill, protocol=pickle.HIGHEST_PROTOCOL)
                summaries.append({"certified": rc.get("certified"),
                                  "gaps": list(rc.get("gaps") or ())[:1]})
                if count % 1000 == 0:
                    print(f"PREPARE certified {count}/{len(universe)} "
                          f"{time.time() - t:.0f}s rss={_rss()}MB", flush=True)
        certify_s = round(time.time() - t, 1)
        # 3. the certification gate and the republished ids
        reasons = cm.refusal_reasons(
            ((light[e.event_id], s) for e, s in zip(universe, summaries)), light)
        del summaries
        replacements = {i: cm.republish_refused(_one(artifact, entry_of[i]), reason).event_id
                        for i, reason in reasons.items()}
        final_ids = {replacements.get(e.event_id, e.event_id) for e in universe}
        print("PREPARE gate", len(reasons), "republished", f"{certify_s}s certification",
              f"rss={_rss()}MB", flush=True)
        # 4. the assignment: one key's records on one shard
        t = time.time()
        systems, _refusals = emit.build_reference_systems(roster)
        units = assignment_units(artifact, universe, certifier, fam, systems)
        shard_of, plan = assign_shards(units, shards)
        if set(shard_of) != set(entry_of):
            raise RuntimeError("the assignment does not place every record once")
        plan["seconds"] = round(time.time() - t, 1)
        print("PREPARE assignment", plan["units"], "units", len(plan["split"]), "split",
              len(plan["case_bound"]), "case-bound", f"{plan['seconds']}s rss={_rss()}MB",
              flush=True)
        # 5. publish each record into the inputs of the shards that read it
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
            with gzip.open(spill_path, "rb") as spill:
                for e in universe:
                    rc = pickle.load(spill)
                    record = _one(artifact, e)
                    if e.event_id in reasons:
                        record = cm.republish_refused(record, reasons[e.event_id])
                    if reasons:
                        record = cm.relink(record, replacements)
                    published, rc, _ref = cm.publish_certificate(record, rc, replacements,
                                                                 final_ids)
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
        spill_path.unlink()
        ids = [e.event_id for e in universe]
        if sorted(i for pairs in owned for i, _ in pairs) != sorted(ids):
            raise RuntimeError("the shard inputs do not own every record once")
        manifest = {
            "schema": MANIFEST_SCHEMA,
            "artifact": {"id": Path(artifact).stem, "sha256": index.sha256, "bytes": index.size},
            "artifact_records": len(index.entries), "families": fams,
            "only": sorted(families) if families else None,
            "artifact_ids": ids, "shards": shards,
            "files": [f.name for f in files], "items": written,
            "files_sha256": [file_sha256(f) for f in files],
            "own": [len(pairs) for pairs in owned],
            "own_sha256": [own_digest(pairs) for pairs in owned],
            "assignment": {"basis": "record_key", "shard_of": [shard_of[i] for i in ids],
                           **plan},
            "head": head, "bound_by": list(bound_by), "identity": identity(database, roster),
            "republished": len(replacements), "certify_s": certify_s,
            "seconds": round(time.time() - t0, 1), "peak_rss_mb": _rss()}
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
