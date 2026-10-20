"""A catalog set (S-g): per-family template catalogs that run as one.

One family is the unit of work: a family is compiled, emitted and assembled
alone into a ``kmc-template-catalog/1`` file.  A ``kmc-catalog-set/1`` manifest
names the member catalogs and what ties them together::

    {"schema": "kmc-catalog-set/1",
     "shared": {database_path, database_content_sha256, code_sources_sha256,
                code_extensions_sha256, roster_sha256, database_head, code_sha},
     "catalogs": [{"path", "sha256", "families", "complete"}, ...],
     "dependencies": [{"kind": "reverse_of", "record_event_id", "family",
                       "requires_record_event_id", "requires_family"}, ...]}

``shared`` is what every member was emitted against; members that differ in
one of its fields are refused by field name, and so is a member whose
database path, database content digest, rmgpy source or extension digest or roster sha256 is
null or missing (unknown is not shared).  ``database_head`` and ``code_sha``
are null outside git; they must be known for every member or for none (then
the content digests identify the inputs).

A set gives exactly what one joint assemble of its members' records gives, or
it is refused when it is built (and again when it is loaded).  The joint
assemble differs from separate ones only through collision resolution, which
groups equal executable keys across every family it holds: reverse links,
derived primaries and the balance check stay inside one family.  Each catalog
lists every key that entered its collision resolution
(``gates.collisions.entry_keys``, refused copies included); a key in two
members is refused by name.  (Resolving it across members instead would need
what a catalog does not keep: base flags, counts, calibration thermo and the
refused copies, i.e. a joint assemble of the rows.)  ``dependencies`` lists every link
from a payload of one member to a record of another.  Only ``reverse_of`` (a
reciprocal pair, admitted or refused together by the emitter) can span
catalogs: a ``merged_into`` or ``derived`` template reference is loaded inside
its own catalog (``load_catalogs`` loads each member alone), so such a link
across members is refused.  A set that lacks one side of a listed dependency is
refused naming both records.

Built by ``python -m rmgpy.kmc.emit catalog-set CATALOG... --out SET.json``;
read by ``run.load_catalogs``, which replaces a set by its verified members.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

CATALOG_SET_SCHEMA = "kmc-catalog-set/1"
CATALOG_SCHEMA = "kmc-template-catalog/1"     # rmgpy.kmc.emit.SCHEMA
KNOWN_INPUTS = ("database_path", "database_content_sha256", "code_sources_sha256",
                "code_extensions_sha256", "roster_sha256")
GIT_INPUTS = ("database_head", "code_sha")
SHARED_INPUTS = KNOWN_INPUTS + GIT_INPUTS
SET_KEYS = ("schema", "shared", "catalogs", "dependencies")
TEMPLATE_LINKS = ("merged_into", "covered_by", "derived_primary")


class CatalogSetError(ValueError):
    """A catalog set that cannot be built or trusted."""


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def payload_links(payloads) -> List[tuple]:
    """Every reference of a payload to another record: (kind, record id,
    required record id, required template key or None)."""
    out = []
    for p in payloads:
        rid = p.get("record_event_id")
        if p.get("reverse_of") is not None:
            out.append(("reverse_of", rid, p["reverse_of"], None))
        for t in p.get("templates", []):
            if t.get("status") == "merged" and t.get("merged_into"):
                out.append(("merged_into", rid, t["merged_into"]["record"], t["merged_into"].get("key")))
            elif t.get("status") == "covered":
                out += [("covered_by", rid, c["record"], c.get("key")) for c in t.get("covered_by", [])]
            elif t.get("status") == "emitted" and t.get("derived"):
                primary = t["derived"]["primary"]
                out.append(("derived_primary", rid, primary["record"], primary.get("key")))
    return out


def read_bytes(path) -> bytes:
    """The one read of a member file: its sha256 is checked and its content
    parsed from these bytes, and the runner loads these bytes."""
    return Path(path).read_bytes()


def _read_member(path) -> dict:
    path = Path(path)
    try:
        raw = read_bytes(path)
    except OSError as error:
        raise CatalogSetError(f"catalog {path}: unreadable ({error})") from error
    return member_from_bytes(path, raw)


def member_from_bytes(path, raw: bytes) -> dict:
    """The member view of catalog ``path`` whose bytes are ``raw`` (never
    read again): its sha256, provenance, payloads, families, completeness and
    collision-entry keys."""
    path = Path(path)
    try:
        doc = json.loads(raw)
    except ValueError as error:
        raise CatalogSetError(f"catalog {path}: unreadable ({error})") from error
    if not isinstance(doc, dict) or doc.get("schema") != CATALOG_SCHEMA:
        raise CatalogSetError(f"catalog {path}: not a {CATALOG_SCHEMA} file")
    provenance = doc.get("provenance")
    payloads = doc.get("payloads")
    if not isinstance(provenance, dict) or not isinstance(payloads, list):
        raise CatalogSetError(f"catalog {path}: no provenance or payloads")
    families = provenance.get("families")
    if not isinstance(families, list) or not families or not all(isinstance(f, str) and f for f in families):
        raise CatalogSetError(f"catalog {path}: provenance.families must be a non-empty list of family names")
    entry_keys = ((doc.get("gates") or {}).get("collisions") or {}).get("entry_keys") \
        if isinstance(doc.get("gates"), dict) else None
    if not isinstance(entry_keys, list) or not all(isinstance(k, str) and k for k in entry_keys):
        raise CatalogSetError(f"catalog {path}: no gates.collisions.entry_keys (the executable keys that "
                              "entered its collision resolution); a set compares them across members, so "
                              "re-assemble the catalog with this code")
    emitted = {t.get("key") for p in payloads for t in p.get("templates", []) if t.get("status") == "emitted"}
    unlisted = sorted(str(k) for k in emitted - set(entry_keys))
    if unlisted:
        raise CatalogSetError(f"catalog {path}: loads key {unlisted[0]} that its gates.collisions.entry_keys "
                              "does not list")
    return {"path": str(path), "raw": raw, "sha256": sha256_bytes(raw), "provenance": provenance,
            "payloads": payloads, "families": sorted(families), "complete": doc.get("complete") is True,
            "entry_keys": sorted(set(entry_keys))}


def _check_shared(members: Sequence[dict]) -> Dict[str, object]:
    for field in KNOWN_INPUTS:
        for m in members:
            value = m["provenance"].get(field)
            if value is None or value == "":
                how = "missing" if field not in m["provenance"] else "null"
                raise CatalogSetError(f"catalog {m['path']}: the shared input {field} is {how}; unknown "
                                      "provenance cannot be shared, so the catalog cannot join a set")
    # GIT_INPUTS are compared like the rest: known for one member and null
    # for another differ, so only all-known-and-equal or all-null pass
    shared = {}
    for field in SHARED_INPUTS:
        values = {json.dumps(m["provenance"].get(field), sort_keys=True) for m in members}
        if len(values) > 1:
            where = ", ".join(f"{m['path']}={m['provenance'].get(field)!r}" for m in members)
            raise CatalogSetError(f"catalogs differ in the shared input {field}: {where}")
        shared[field] = members[0]["provenance"].get(field)
    return shared


def _check_collisions(members: Sequence[dict]) -> None:
    """Refuse a key that entered the collision resolution of two members: one
    joint assemble resolves equal keys together (merged or refused), which
    separate catalogs cannot reproduce."""
    owner: Dict[str, int] = {}
    for i, m in enumerate(members):
        for key in m["entry_keys"]:
            if key not in owner:
                owner[key] = i
                continue
            first = members[owner[key]]

            def holders(member):
                rids = sorted({str(p.get("record_event_id")) for p in member["payloads"]
                               for t in p.get("templates", []) if t.get("key") == key})
                return f"records {rids}" if rids else "no loaded record (refused inside its catalog)"
            raise CatalogSetError(
                f"executable key {key} entered the collision resolution of two catalogs: {first['path']} "
                f"(families {first['families']}; {holders(first)}) and {m['path']} (families {m['families']}; "
                f"{holders(m)}); one joint assemble resolves equal keys together (merged or refused), which "
                "separate catalogs cannot reproduce: assemble these families into one catalog")


def _check_disjoint(members: Sequence[dict]) -> Dict[str, tuple]:
    """record id -> (member index, family); families and records are owned once."""
    owner_family: Dict[str, str] = {}
    for m in members:
        for family in m["families"]:
            if family in owner_family:
                raise CatalogSetError(f"catalogs {owner_family[family]} and {m['path']} share the family "
                                      f"{family!r}; a set holds one catalog per family")
            owner_family[family] = m["path"]
    records: Dict[str, tuple] = {}
    for i, m in enumerate(members):
        for p in m["payloads"]:
            rid = p.get("record_event_id")
            if rid in records:
                raise CatalogSetError(f"record {rid} is in two catalogs of the set "
                                      f"({members[records[rid][0]]['path']} and {m['path']})")
            records[rid] = (i, p.get("family"))
    return records


def _cross_links(members: Sequence[dict], records: Dict[str, tuple]) -> List[dict]:
    """The reverse_of links between members; a template link between members,
    or a link to a record in no member of a complete member, is an error."""
    found = {}
    for i, m in enumerate(members):
        for kind, rid, req, key in payload_links(m["payloads"]):
            target = records.get(req)
            if target is None:
                if kind in TEMPLATE_LINKS or m["complete"]:
                    raise CatalogSetError(
                        f"{kind}: record {rid} (family {records[rid][1]}) of {m['path']} requires record "
                        f"{req}, which is in no catalog of the set")
                continue
            if target[0] == i:
                continue
            if kind in TEMPLATE_LINKS:
                raise CatalogSetError(
                    f"{kind}: record {rid} (family {records[rid][1]}) of {m['path']} requires record {req} "
                    f"(family {target[1]}) of {members[target[0]]['path']}; a merged, covered or derived template loads "
                    "only inside its own catalog, so its owner must be in the same catalog")
            found[(kind, rid, req)] = {"kind": kind, "record_event_id": rid, "family": records[rid][1],
                                       "requires_record_event_id": req, "requires_family": target[1]}
    return [found[k] for k in sorted(found)]


def check_union(members: Sequence[dict]) -> List[dict]:
    """Everything that must hold for catalogs used together, whether one set,
    several sets or plain catalogs given to the runner: shared inputs known
    and equal, families and records owned once, no collision-entry key in two
    members, and no template link across members or dangling link of a
    complete member.  Returns the reverse_of links between members.
    Completeness (``--allow-incomplete``) relaxes none of it."""
    _check_shared(members)
    records = _check_disjoint(members)
    _check_collisions(members)
    return _cross_links(members, records)


def build_catalog_set(paths: Sequence, out) -> dict:
    """Write the ``kmc-catalog-set/1`` manifest of the catalogs ``paths`` to
    ``out`` and return it.  Refused (``CatalogSetError``): a catalog that is
    not readable, a shared input unknown for a member or different between
    members, members that share a family, a record or a collision-entry key,
    or a link to a record outside the set."""
    if not paths:
        raise CatalogSetError("a catalog set needs at least one catalog")
    members = [_read_member(p) for p in paths]
    shared = _check_shared(members)
    dependencies = check_union(members)
    base = Path(out).resolve().parent
    doc = {"schema": CATALOG_SET_SCHEMA, "shared": shared,
           "catalogs": [{"path": os.path.relpath(Path(m["path"]).resolve(), base), "sha256": m["sha256"],
                         "families": m["families"], "complete": m["complete"]} for m in members],
           "dependencies": dependencies}
    Path(out).write_text(json.dumps(doc, sort_keys=True, indent=1, allow_nan=False) + "\n")
    return doc


def is_catalog_set(doc) -> bool:
    return isinstance(doc, dict) and doc.get("schema") == CATALOG_SET_SCHEMA


def read_set_members(set_path, doc) -> List[Tuple[Path, bytes]]:
    """(path, bytes) of each member catalog of the set ``doc`` read from
    ``set_path``, after verifying the members' sha256, the shared inputs, the
    collision-entry keys and the dependencies (``CatalogSetError``).  Each
    member is read once: the bytes returned are the bytes verified, and the
    runner loads them (never the file again)."""
    set_path = Path(set_path)
    missing = [k for k in SET_KEYS if k not in doc]
    if missing:
        raise CatalogSetError(f"catalog set {set_path}: missing fields {missing}")
    entries = doc["catalogs"]
    if not isinstance(entries, list) or not entries or not all(
            isinstance(e, dict) and isinstance(e.get("path"), str) for e in entries):
        raise CatalogSetError(f"catalog set {set_path}: catalogs must be a non-empty list of {{path, sha256, ...}}")
    paths = [set_path.parent / e["path"] for e in entries]
    members = []
    for entry, path in zip(entries, paths):
        member = _read_member(path)
        if member["sha256"] != entry.get("sha256"):
            raise CatalogSetError(f"catalog set {set_path}: member {path} has sha256 {member['sha256'][:16]}..., "
                                  f"the set records {str(entry.get('sha256'))[:16]}...")
        members.append(member)
    shared = _check_shared(members)
    for field in SHARED_INPUTS:
        if doc["shared"].get(field) != shared[field]:
            raise CatalogSetError(f"catalog set {set_path}: shared input {field} of the members "
                                  f"({shared[field]!r}) differs from the set's ({doc['shared'].get(field)!r})")
    records = _check_disjoint(members)
    _check_collisions(members)
    listed = set()
    for dep in doc["dependencies"]:
        rid, req = dep.get("record_event_id"), dep.get("requires_record_event_id")
        absent = [(side, r, f) for side, r, f in (("the dependent", rid, dep.get("family")),
                                                  ("the required", req, dep.get("requires_family")))
                  if r not in records]
        if absent:
            raise CatalogSetError(
                f"catalog set {set_path} lacks {absent[0][0]} side of the dependency {dep.get('kind')}: record "
                f"{rid} (family {dep.get('family')}) requires record {req} (family {dep.get('requires_family')}); "
                f"record {absent[0][1]} (family {absent[0][2]}) is in no member")
        listed.add((dep.get("kind"), rid, req))
    for dep in _cross_links(members, records):
        if (dep["kind"], dep["record_event_id"], dep["requires_record_event_id"]) not in listed:
            raise CatalogSetError(
                f"catalog set {set_path} does not list the dependency {dep['kind']}: record "
                f"{dep['record_event_id']} (family {dep['family']}) requires record "
                f"{dep['requires_record_event_id']} (family {dep['requires_family']})")
    return [(path, m["raw"]) for path, m in zip(paths, members)]
