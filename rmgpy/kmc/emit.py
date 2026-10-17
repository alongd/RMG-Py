"""Template catalog emission (I-155): a compiled event-set artifact plus a roster
of declared reference systems becomes one versioned, provenance-stamped,
gate-checked ``kmc-template-catalog/1`` file the melt engine can load.

Two steps, both from the command line::

    python -m rmgpy.kmc.emit ARTIFACT ROSTER --database DB --out rows.jsonl
        [--chunk i/n] [--family F]... [--limit N] [--no-contexts]
    python -m rmgpy.kmc.emit assemble ROWS... --out catalog.json [--database DB]

**Emit.** For each record (or the shard ``--chunk i/n``: records ``i::n`` of the
event-id-sorted list): certify the records, match the reference cases of the
roster's systems, ``template_emission.emit_record`` with the production rates
and the reference Kc.  One JSON line per record: the payload plus status,
refusal and timing.  An error is a row (``status: error``), never silent.  The
first line of a rows file is its header (artifact, database, code, roster and
shard identity).

**Roster.** A JSON file ``{"schema": "kmc-reference-roster/1", "polymer":
<LinearHomopolymerSpec mapping>, "systems": {name: [participant, ...]}}``.  A
participant is ``{"chain": {"dp", "head", "tail", "radical_unit"}}`` (a chain of
the declared repeat unit; ``head``/``tail`` name an end variant of the matched
provider's ``END_CHEMISTRIES`` vocabulary, ``null`` the artificial H-capped end),
``{"smiles": ...}`` or ``{"adjacency_list": ...}``.  The declared repeat unit is
matched to a registered provider (``boundaries.PROVIDERS``, ``REPEAT_UNIT`` and
``spec_from_declaration``) by graph isomorphism; a declaration no provider
matches is refused (``refused:no_provider``), never dropped.

**Assemble.** Merge the shards, then ``resolve_collisions``, ``derive_weights``
and ``check_cycles`` over the whole set.  The catalog is written even when a gate
fails; ``complete`` is true only when every gate passes and no record errored.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import traceback
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any, Callable, Optional

from rmgpy.kmc import boundaries as b
from rmgpy.kmc import template_emission as te
from rmgpy.kmc import windows as win

SCHEMA = "kmc-template-catalog/1"
ROWS_SCHEMA = "kmc-emit-rows/1"
ROSTER_SCHEMA = "kmc-reference-roster/1"
# Kept equal to the acceptance harness' labels so that rows of the two agree
# (neither string selects data: the Kc label tags the reference thermo, the
# rate label tags the production-rate provenance)
REFERENCE_THERMO_LABEL = "i147-reference"
REFERENCE_RATES_LABEL = "i138-reference"
PAYLOAD_FIELDS = ("key", "name", "family", "status", "direction", "executable", "weight", "derived")
HEADER_MATCH = ("schema", "artifact_id", "artifact_sha256", "artifact_records", "families",
                "database_path", "database_head", "code_sha", "roster_sha256", "temperatures",
                "selection")
CYCLE_REFUSED_FRACTION = 0.05


# ---------------------------------------------------------------------------
# reference systems from a declaration
# ---------------------------------------------------------------------------
class RosterRefusal(ValueError):
    """A roster participant that cannot be built (``status`` names why)."""

    def __init__(self, status, reason):
        super().__init__(f"{status}: {reason}")
        self.status = status
        self.reason = reason


@dataclass(frozen=True)
class Roster:
    polymer: Any
    systems: dict
    sha256: str
    path: str


@dataclass
class Participant:
    """One reference participant: the spec the windows machinery reads and the
    builder of its fresh proxy build (a build is made per use, as the
    acceptance harness does)."""

    label: str
    spec: Any
    make: Callable

    def build(self, name):
        return self.make(name)


def load_roster(path) -> Roster:
    from rmgpy.polymer_spec import normalized_linear_homopolymer_spec
    raw = Path(path).read_bytes()
    data = json.loads(raw)
    if not isinstance(data, dict) or data.get("schema") != ROSTER_SCHEMA:
        raise RosterRefusal("refused:roster", f"roster schema must be {ROSTER_SCHEMA!r}")
    if set(data) != {"schema", "polymer", "systems"}:
        raise RosterRefusal("refused:roster", "a roster has exactly schema, polymer and systems")
    polymer = normalized_linear_homopolymer_spec(data["polymer"])
    systems = data["systems"]
    if not isinstance(systems, dict) or not all(isinstance(v, list) for v in systems.values()):
        raise RosterRefusal("refused:roster", "systems must map a name to a participant list")
    return Roster(polymer, systems, hashlib.sha256(raw).hexdigest(), str(path))


def _normal_graph(graph):
    """(atom label by key, bond order by key pair, port atom by name) of a graph
    declaration with explicit single-bonded terminal H atoms folded into their
    neighbour's hydrogen count."""
    ports = {name: port.atom for name, port in graph.ports}
    edges = defaultdict(dict)
    for bond in graph.bonds:
        left, right = bond.atoms
        edges[left][right] = bond.order
        edges[right][left] = bond.order
    hydrogens = {atom.key: atom.hydrogens for atom in graph.atoms}
    folded = set()
    for atom in graph.atoms:
        if (atom.element == "H" and atom.isotope == -1 and atom.key not in ports.values()
                and len(edges[atom.key]) == 1 and atom.radical_electrons == 0):
            (other, order), = edges[atom.key].items()
            if order == 1.0:
                hydrogens[other] += 1
                folded.add(atom.key)
    labels = {a.key: (a.element, a.isotope, a.charge, a.radical_electrons, hydrogens[a.key])
              for a in graph.atoms if a.key not in folded}
    bonds = {k: {o: v for o, v in edges[k].items() if o not in folded} for k in labels}
    return labels, bonds, ports


def isomorphic(first, second) -> bool:
    """Whether two graph declarations are the same repeat unit: a bijection of
    atoms preserving element, isotope, charge, radicals, hydrogen count, bond
    orders, and the head and tail ports."""
    (la, ba, pa), (lb, bb, pb) = _normal_graph(first), _normal_graph(second)
    if set(pa) != set(pb) or len(la) != len(lb) or sorted(la.values()) != sorted(lb.values()):
        return False
    if any(first_port.bond_order != second_port.bond_order for (_n, first_port), (_m, second_port)
           in zip(first.ports, second.ports)):
        return False
    order, seen = [], set()
    stack = [pa["head"]]
    while stack:
        key = stack.pop()
        if key not in seen:
            seen.add(key)
            order.append(key)
            stack.extend(sorted(set(ba[key]) - seen, reverse=True))
    if len(order) != len(la):
        return False

    def extend(index, mapping, used):
        if index == len(order):
            return True
        key = order[index]
        for candidate in sorted(lb):
            if candidate in used or lb[candidate] != la[key]:
                continue
            if index == 0 and candidate != pb["head"]:
                continue
            if len(ba[key]) != len(bb[candidate]):
                continue
            if any(mapping[other] not in bb[candidate] or bb[candidate][mapping[other]] != order_
                   for other, order_ in ba[key].items() if other in mapping):
                continue
            if sum(1 for other in ba[key] if other in mapping) != \
                    sum(1 for other in bb[candidate] if other in used):
                continue
            mapping[key] = candidate
            used.add(candidate)
            if extend(index + 1, mapping, used):
                return True
            del mapping[key]
            used.discard(candidate)
        return False
    mapping = {}
    if not extend(0, mapping, set()):
        return False
    return mapping[pa["tail"]] == pb["tail"]


def match_provider(graph):
    """The registered provider whose ``REPEAT_UNIT`` is the declared repeat unit,
    or None."""
    from rmgpy.polymer_spec import GraphDeclaration
    for name in sorted(b.PROVIDERS):
        unit = getattr(b.PROVIDERS[name], "REPEAT_UNIT", None)
        if unit is not None and isomorphic(
                graph, GraphDeclaration.from_mapping(unit, required_ports=("head", "tail"))):
            return b.PROVIDERS[name]
    return None


def explicit_build(name, smiles=None, adjacency_list=None):
    """A provider-free build of an explicit small molecule: no ports, one
    repeat coordinate per heavy atom."""
    from rmgpy.molecule.molecule import Molecule
    if smiles is not None:
        molecule = Molecule(smiles=smiles)
    else:
        molecule = Molecule().from_adjacency_list(adjacency_list)
    objects = list(molecule.atoms)
    molecule.update()
    molecule.assign_atom_ids()
    rows = dict(b.rows_of_molecule(molecule))
    heavy = sorted(a for a, row in rows.items() if row[0] != "H")
    metadata = b.BoundaryMetadata("small", (), tuple((a, k, 0) for k, a in enumerate(heavy)))
    metadata.check(rows, name)
    return b.ProxyBuild(name, molecule, metadata, tuple(atom.id for atom in objects))


CHAIN_KEYS = {"dp", "head", "tail", "radical_unit"}


def _participant(provider, graph, entry, label):
    if not isinstance(entry, dict) or len(entry) != 1:
        raise RosterRefusal("refused:roster", f"{label}: exactly one of chain, smiles, adjacency_list")
    (kind, value), = entry.items()
    if kind == "chain":
        if not isinstance(value, dict) or "dp" not in value or set(value) - CHAIN_KEYS:
            raise RosterRefusal("refused:roster", f"{label}: a chain declares dp and optionally "
                                                  f"head, tail, radical_unit")
        if provider is None:
            raise RosterRefusal("refused:no_provider",
                                f"{label}: no registered provider declares this repeat unit "
                                f"(registered: {sorted(b.PROVIDERS)})")
        try:
            spec = provider.spec_from_declaration(value["dp"], value.get("head"), value.get("tail"),
                                                  value.get("radical_unit"))
            participant = Participant(label, spec, lambda name, s=spec: win.build_spec(s, name))
            participant.build("roster_check")
        except RosterRefusal:
            raise
        except Exception as error:   # reported as a refusal, never dropped
            raise RosterRefusal("refused:variant", f"{label}: {type(error).__name__}: {error}")
        return participant
    if kind in ("smiles", "adjacency_list"):
        try:
            def make(name, kind=kind, value=value):
                return explicit_build(name, **{kind: value})
            make("roster_check")
        except Exception as error:
            raise RosterRefusal("refused:variant", f"{label}: {type(error).__name__}: {error}")
        return Participant(label, ("small", f"{kind}:{value}"), make)
    raise RosterRefusal("refused:roster", f"{label}: unknown participant kind {kind!r}")


def build_reference_systems(roster: Roster):
    """(systems, refusals): ``systems`` maps a name to its participants; a
    system with a refused participant is left out and each refusal is a row
    (``kind: roster_refusal``, system, participant, status, reason)."""
    from rmgpy.polymer_spec import GraphDeclaration  # noqa: F401  (documents the declaration)
    provider = match_provider(roster.polymer.repeat_unit)
    systems, refusals = {}, []
    for name, entries in roster.systems.items():
        built, bad = [], []
        for index, entry in enumerate(entries):
            try:
                built.append(_participant(provider, roster.polymer.repeat_unit, entry,
                                          f"{name}[{index}]"))
            except RosterRefusal as error:
                bad.append({"kind": "roster_refusal", "system": name, "participant": index,
                            "status": error.status, "reason": error.reason})
        if bad:
            refusals.extend(bad)
        elif built:
            systems[name] = built
        else:
            refusals.append({"kind": "roster_refusal", "system": name, "participant": None,
                             "status": "refused:roster", "reason": "an empty system"})
    return systems, refusals


# ---------------------------------------------------------------------------
# emission (the acceptance harness' path, on declared systems)
# ---------------------------------------------------------------------------
def load_database(families, path):
    """The compile's own load (``rmgpy.kmc.compile.compile_artifact``)."""
    from rmgpy.data.rmg import RMGDatabase
    from rmgpy.kmc.compiler import prepare_rate_rules
    path = Path(path)
    db = RMGDatabase()
    db.load_kinetics(str(path / "input/kinetics"), reaction_libraries=[], seed_mechanisms=None,
                     kinetics_families=list(families), kinetics_depositories=["training"])
    db.load_thermo(str(path / "input/thermo"), thermo_libraries=["primaryThermoLibrary"],
                   depository=True)
    prepare_rate_rules(db.kinetics, db.thermo, verbose=True)
    return db


def rss_mb():
    try:
        with open("/proc/self/status") as handle:
            for line in handle:
                if line.startswith("VmHWM:"):
                    return int(line.split()[1]) // 1024
    except OSError:
        pass
    return -1


def molecularity(family, direction):
    from rmgpy.kmc import template_emission as te
    return te._molecularity(family, direction)


def reference_cases(certifier, fam, systems, species):
    """Every member case of every family direction on every subset of every
    system's participants of the direction's molecularity, keyed by (family,
    template-side species, template, other-side species)."""
    out = defaultdict(list)
    stats = Counter()
    for name, family in sorted(fam.items()):
        directions = ["forward"] + (["reverse"] if family.reverse_template is not None else [])
        for direction in directions:
            n = molecularity(family, direction)
            for participants in systems.values():
                for subset in sorted(set(combinations(range(len(participants)), n))):
                    chosen = [participants[i] for i in subset]
                    builds = [p.build(f"p{k}") for k, p in enumerate(chosen)]
                    specs = [p.spec for p in chosen]
                    cases = win.member_cases(certifier, family, builds, direction,
                                             label=f"{name}|{direction}|{subset}")
                    stats[f"{name}|{direction}"] += len(cases)
                    for case in cases:
                        mols = win.template_molecules(case.witness, case.direction)
                        other = win.template_molecules(
                            case.witness, "reverse" if case.direction == "forward" else "forward")
                        key = (name, tuple(sorted(species(m) for m in mols)),
                               ";".join(case.pathway), tuple(sorted(species(m) for m in other)))
                        out[key].append((case, specs))
    return out, dict(stats)


def record_key(certifier, record, species):
    """(family, template-side species, template, other-side species): the
    record's channel (members of other channels on the same species and
    template are other sites)."""
    _side, mols, others = certifier.template_side(record)
    return (record["family"], tuple(sorted(species(m) for m in mols)), record["template"],
            tuple(sorted(species(m) for m in others)))


def bind(certifier, record, rc, case, owner):
    """``recertify`` the record on the member's template side through an RMG
    isomorphism of the template sides (``windows_acceptance.bind_record``,
    any direction)."""
    from rmgpy.kmc import certificates as cm
    from rmgpy.kmc import windows as win
    _side, mols, _others = certifier.template_side(record)
    proxy = win.template_molecules(case.witness, case.direction)
    roots_r = [cm.aromatic_root(m) for m in mols]
    roots_p = [cm.aromatic_root(m) for m in proxy]
    amap, order = {}, []
    for mr in roots_r:
        for k, mp in enumerate(roots_p):
            if k in order or len(mp.atoms) != len(mr.atoms):
                continue
            maps = mr.find_isomorphism(mp, strict=False)
            if maps:
                order.append(k)
                amap.update({a.id: c.id for a, c in maps[0].items()})
                break
        else:
            return {"same_selection": False, "reason": "no isomorphic member fragment"}
    # the member's molecules in the record's reactant order: a certificate
    # mapping's identity carries the reactant index, so an unequal order
    # (the member's components go by least atom id) swaps the roles
    proxy = [proxy[k] for k in order]
    res = certifier.recertify(record, proxy, amap, certificate=rc, rate_owner=owner)
    return {"same_selection": res["same_selection"], "same_support": res["same_support"],
            "refusals": res["refusals"][:2]}


def select_records(records, chunk=None, limit=None):
    """The event-id-sorted records, then shard ``chunk`` ("i/n": ``i::n``), then
    at most ``limit``."""
    picked = sorted(records, key=lambda r: r["event_id"])
    if chunk:
        k, n = parse_chunk(chunk)
        picked = picked[k::n]
    if limit:
        picked = picked[:limit]
    return picked


def parse_chunk(chunk):
    if chunk is None:
        return 0, 1
    k, n = (int(x) for x in str(chunk).split("/"))
    if not 0 <= k < n:
        raise ValueError(f"bad chunk {chunk!r}: need 0 <= i < n")
    return k, n


def git_head(path) -> Optional[str]:
    try:
        out = subprocess.run(["git", "-C", str(path), "rev-parse", "HEAD"], capture_output=True,
                             text=True, timeout=60)
        return out.stdout.strip() or None if out.returncode == 0 else None
    except (OSError, subprocess.SubprocessError):
        return None


def code_sha() -> Optional[str]:
    """HEAD of the code repository, ``-dirty`` when tracked files differ."""
    here = Path(__file__).resolve().parent
    head = git_head(here)
    if head is None:
        return None
    try:
        dirty = subprocess.run(["git", "-C", str(here), "status", "--porcelain", "-uno"],
                               capture_output=True, text=True, timeout=60).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        dirty = ""
    return head + ("-dirty" if dirty else "")


class Emitter:
    """The emission session of one artifact against a roster's systems."""

    def __init__(self, artifact, database, systems, families=None, contexts=True):
        self.artifact = Path(artifact)
        self.database = Path(database)
        self.systems = systems
        self.only = list(families) if families else None
        self.contexts = contexts
        self.state = {}

    def init(self):
        from rmgpy.kmc import certificates as cm
        from rmgpy.kmc import channels as chm
        from rmgpy.kmc.event_record import EventRecord
        from rmgpy.kmc.reference_thermo import GasPhaseRMGReferenceThermo
        S = self.state
        art = json.loads(self.artifact.read_text())
        S["artifact_records"] = len(art["records"])
        if self.only:
            # one family's records (certified among themselves: a reversed
            # record's rate owner is in its own family)
            art["records"] = [r for r in art["records"] if r["family"] in self.only]
        families = sorted({r["family"] for r in art["records"]})
        t = time.time()
        db = load_database(families, self.database)
        fam = db.kinetics.families
        S.update(db=db, fam=fam, families=families, load_s=round(time.time() - t, 1))
        certifier = cm.Certifier(fam, max_channel_sets=8, max_entries=1024)
        t = time.time()
        records = [EventRecord.from_dict(r) for r in art["records"]]
        certified, store, _replacements = cm.certify_records(records, fam, certifier, required=True)
        certified = [r.to_dict() if hasattr(r, "to_dict") else r for r in certified]
        S.update(certify_s=round(time.time() - t, 1), records=certified, store=store,
                 by_id={r["event_id"]: r for r in certified}, certifier=certifier)
        species = chm.SpeciesKeys()
        t = time.time()
        cases, stats = reference_cases(certifier, fam, self.systems, species)
        S.update(cases=cases, case_stats=stats, cases_s=round(time.time() - t, 1), species=species,
                 bounds={n: win.span_bound(win.family_constants(f)) for n, f in fam.items()},
                 counter=te.RmgCounter(fam),
                 kc=te.ThermoKc(GasPhaseRMGReferenceThermo(db.thermo, REFERENCE_THERMO_LABEL)),
                 rates=te.ProductionRates(db, families, database_path=self.database,
                                          database_sha=REFERENCE_RATES_LABEL))

    @property
    def records(self):
        return self.state["records"]

    def run_record(self, record):
        from rmgpy.kmc import certificates as cm
        S = self.state
        contexts = self.contexts
        certifier, fam = S["certifier"], S["fam"]
        start, cpu = time.time(), time.process_time()
        rc = S["store"].get((record.get("rate_certificate") or {}).get("sha256"))
        row = {"kind": "row", "event_id": record["event_id"][:24], "family": record["family"],
               "orientation": record["orientation"], "template": record["template"][:100],
               "record_status": record["status"], "rate_kind": (record.get("rate_source") or {})
               .get("kind"), "certified": (rc or {}).get("certified")}
        family = fam[record["family"]]
        why = te.upstream_refusal(record, rc, S["by_id"])
        members, member_rows = [], []
        try:
            if why is None:
                key = record_key(certifier, record, S["species"])
                found = S["cases"].get(key, [])
                row["members"] = len(found)
                if found:
                    owner = cm._rate_owner(record, S["by_id"])
                    row["bind"] = bind(certifier, record, rc, found[0][0], owner)
                    if not row["bind"].get("same_selection"):
                        why = f"refused:recertify: {row['bind']}"[:300]
                if why is None:
                    for case, specs in found:
                        out = []
                        keep = []
                        res = win.run_member(certifier, family, case,
                                             bound=S["bounds"][record["family"]],
                                             contexts=contexts, specs=specs, keep=keep,
                                             contexts_out=out)
                        if out:
                            ctx = out[0]
                        elif keep:
                            fp = keep[0]
                            ctx = ((win.ContextResult(tuple((p, None) for p in
                                                            win.b.vector_ports(case.witness)),
                                                      fp.status, fp.reason, fp.witness, fp.window,
                                                      fp.observation, fp),), ())
                        else:
                            ctx = ((), ())
                        members.append((case, res, ctx))
                        member_rows.append({k: res.get(k) for k in (
                            "label", "status", "reason", "window_atoms", "context_classes",
                            "covered", "seconds")})
            if why is not None:
                payload = te.finish({"schema": te.SCHEMA, "status": "refused",
                                     "record_event_id": record["event_id"],
                                     "refusal": {"by": "upstream", "status": te.status_of(why),
                                                 "reason": why}, "templates": [],
                                     "contexts": []})
            else:
                payload = te.emit_record(record, family=family, members=members,
                                         counter=S["counter"], rates=S["rates"],
                                         record_certificate=rc, by_id=S["by_id"],
                                         check_upstream=False, kc=S["kc"])
            row["payload"] = payload
            row["status"] = payload["status"]
            row["refusal"] = (payload.get("refusal") or {}).get("status")
        except te.EmissionError as error:
            row.update(status="error", refusal="error:emission", reason=str(error)[:400])
        except Exception:   # reported, never hidden
            row.update(status="error", refusal="error:exception",
                       reason=traceback.format_exc()[-1500:])
        row["member_rows"] = member_rows
        row["seconds"] = round(time.time() - start, 1)
        row["cpu"] = round(time.process_time() - cpu, 1)
        row["rss_mb"] = rss_mb()
        return row


def shard_header(emitter, roster, refusals, picked, chunk, partial):
    from rmgpy.kmc.compiler import DEFAULT_T_GRID
    S = emitter.state
    artifact = emitter.artifact
    return {"kind": "header", "schema": ROWS_SCHEMA, "artifact_id": artifact.stem,
            "artifact_sha256": hashlib.sha256(artifact.read_bytes()).hexdigest(),
            "artifact_records": S["artifact_records"], "families": S["families"],
            "database_path": str(emitter.database), "database_head": git_head(emitter.database),
            "code_sha": code_sha(), "roster_sha256": roster.sha256,
            "temperatures": [float(t) for t in DEFAULT_T_GRID], "chunk": chunk,
            "selection": "partial" if partial else "all",
            "assigned": [r["event_id"] for r in picked], "roster_refusals": refusals,
            "systems": {n: [p.label for p in ps] for n, ps in emitter.systems.items()}}


def run_emit(args):
    roster = load_roster(args.roster)
    systems, refusals = build_reference_systems(roster)
    for refusal in refusals:
        print("ROSTER", json.dumps(refusal), flush=True)
    emitter = Emitter(args.artifact, args.database, systems, families=args.family,
                      contexts=not args.no_contexts)
    t0 = time.time()
    emitter.init()
    picked = select_records(emitter.records, args.chunk, args.limit)
    header = shard_header(emitter, roster, refusals, picked, args.chunk,
                          partial=bool(args.family or args.limit))
    print("INIT", {k: emitter.state.get(k) for k in ("load_s", "certify_s", "cases_s", "case_stats")},
          len(emitter.records), "records, selected", len(picked), flush=True)
    errors = 0
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as log:
        log.write(json.dumps(header, default=str) + "\n")
        for refusal in refusals:
            log.write(json.dumps(refusal) + "\n")
        for k, record in enumerate(picked):
            row = emitter.run_record(record)
            errors += row["status"] == "error"
            log.write(json.dumps(row, default=str) + "\n")
            log.flush()
            print(f"[{time.time() - t0:6.0f}s] {k + 1}/{len(picked)} {row['event_id'][:16]} "
                  f"{row['family']} {row['orientation']} -> {row['status']} "
                  f"{row.get('refusal') or ''} t={row['seconds']}s rss={row['rss_mb']}MB",
                  flush=True)
    return 1 if errors else 0


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------
@dataclass
class Services:
    """What the balance functions need from the database."""

    counter: Any
    kc: Any


def load_services(database, families):
    from rmgpy.kmc.reference_thermo import GasPhaseRMGReferenceThermo
    db = load_database(families, database)
    return Services(te.RmgCounter(db.kinetics.families),
                    te.ThermoKc(GasPhaseRMGReferenceThermo(db.thermo, REFERENCE_THERMO_LABEL)))


def dumps(catalog) -> str:
    return json.dumps(catalog, sort_keys=True, indent=1, allow_nan=False) + "\n"


def _canon(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"))


def read_rows(paths):
    """([(path, header or None, rows)], failures): rows files in a stable order."""
    files, failures = [], []
    for path in paths:
        header, rows, refusals = None, [], []
        for number, line in enumerate(Path(path).read_text().splitlines(), 1):
            if not line.strip():
                continue
            try:
                item = json.loads(line)
            except ValueError as error:
                failures.append(f"{path}:{number}: unreadable row ({error})")
                continue
            kind = item.get("kind")
            if kind == "header":
                if header is not None:
                    failures.append(f"{path}: more than one header")
                header = item
            elif kind == "roster_refusal":
                refusals.append(item)
            else:
                rows.append(item)
        files.append({"path": str(path), "header": header, "rows": rows})
        if header is None:
            failures.append(f"{path}: no shard header (artifact, database, roster and shard "
                            f"identity unknown)")

    def order(f):
        h = f["header"] or {}
        return (parse_chunk(h.get("chunk")) if h else (10 ** 9, 1), _canon(h), f["path"])
    files.sort(key=order)
    return files, failures


def stripped_payloads(payloads):
    """The fields ``loaded_templates`` reads, emitted templates only."""
    out = []
    for p in payloads:
        templates = [{k: t[k] for k in PAYLOAD_FIELDS if k in t}
                     for t in p.get("templates", []) if t.get("status") == "emitted"]
        if templates:
            out.append({"record_event_id": p["record_event_id"], "family": p.get("family"),
                        "reverse_of": p.get("reverse_of"), "templates": templates})
    return sorted(out, key=lambda p: p["record_event_id"])


def _refusal_reason(row):
    payload = row.get("payload") or {}
    if row.get("reason"):
        return row["reason"]
    if (payload.get("refusal") or {}).get("reason"):
        return payload["refusal"]["reason"]
    refused = Counter((t.get("refusal") or {}).get("status") for t in payload.get("templates", [])
                      if t.get("status") == "refused")
    return ", ".join(f"{k} x{v}" for k, v in sorted(refused.items(), key=str)) or None


def assemble(paths, *, database=None, services=None) -> dict:
    """The catalog of the rows files ``paths`` (see the module docstring)."""
    failures = []
    files, read_failures = read_rows(paths)
    failures += read_failures
    headers = [f["header"] for f in files if f["header"]]
    first = headers[0] if headers else {}
    for field in HEADER_MATCH:
        values = {_canon(h.get(field)) for h in headers}
        if len(values) > 1:
            failures.append(f"shard headers disagree on {field}: {sorted(values)}")
    # coverage: every shard once, every assigned record a row, the whole artifact assigned
    ns = {parse_chunk(h.get("chunk"))[1] for h in headers}
    complete_coverage = bool(headers)
    if len(ns) > 1:
        failures.append(f"shard counts disagree: {sorted(ns)}")
        complete_coverage = False
    elif ns:
        (n,) = ns
        seen = Counter(parse_chunk(h.get("chunk"))[0] for h in headers)
        missing = [k for k in range(n) if not seen[k]]
        if len(missing) <= 8:
            failures += [f"missing shard {k}/{n}" for k in missing]
        elif missing:
            failures.append(f"missing shards {missing[0]}/{n} .. {missing[-1]}/{n} "
                            f"({len(missing)} of {n})")
        complete_coverage = complete_coverage and not missing
        for k in range(n):
            if seen[k] > 1:
                failures.append(f"duplicate shard {k}/{n} ({seen[k]} files)")
                complete_coverage = False
    if first.get("selection", "all") != "all":
        failures.append("the rows are a partial selection (--family or --limit)")
        complete_coverage = False
    by_id = {}
    for f in files:
        for row in f["rows"]:
            by_id.setdefault(row["event_id"], row)
    assigned = {i[:24] for h in headers for i in h.get("assigned", [])}
    if headers and first.get("selection", "all") == "all" and \
            len(assigned) != first.get("artifact_records"):
        failures.append(f"the shards assign {len(assigned)} of the artifact's "
                        f"{first.get('artifact_records')} records")
        complete_coverage = False
    if headers:
        absent = sorted(assigned - set(by_id))
        if absent:
            failures.append(f"{len(absent)} assigned record(s) have no row, e.g. {absent[0]}")
            complete_coverage = False
        extra = sorted(set(by_id) - assigned)
        if extra:
            failures.append(f"{len(extra)} row(s) were not assigned to any shard, e.g. {extra[0]}")
    roster_refusals = {}
    for f in files:
        for h in ([f["header"]] if f["header"] else []):
            for item in h.get("roster_refusals", []):
                roster_refusals[_canon(item)] = item
    for item in sorted(roster_refusals.values(), key=_canon):
        failures.append(f"roster {item.get('status')}: system {item.get('system')} "
                        f"participant {item.get('participant')}: {item.get('reason')}")
    rows = [by_id[k] for k in sorted(by_id)]
    for row in rows:
        if row.get("status") == "error":
            failures.append(f"error {row['event_id']}: {str(row.get('reason', ''))[:300]}")
    if rows and not any(r.get("status") in ("emitted", "partial") for r in rows):
        failures.append("every record was refused (a vacuous catalog)")

    have = [k for k, r in enumerate(rows) if r.get("payload")]
    payloads = [rows[k]["payload"] for k in have]
    collisions = []
    try:
        merged, collisions = te.resolve_collisions(payloads)
    except Exception as error:    # a split reciprocal pair, a bad weight
        failures.append(f"catalog consistency: {type(error).__name__}: {error}"[:600])
        merged = payloads
    refused_collisions = [c for c in collisions if c["disposition"] != "merged"]
    for c in refused_collisions:
        failures.append(f"{c['disposition']} {c['key'][:16]}: {c['copies']} copies, "
                        f"log10 diameter {c['log10_diameter']} (records {c['records']})")
    payloads = merged
    cycle_table = []
    if payloads:
        try:
            services = services or load_services(database or first.get("database_path"),
                                                 first.get("families", []))
            payloads = te.derive_weights(payloads, counter=services.counter)
            payloads, cycle_table = te.check_cycles(payloads, counter=services.counter,
                                                    kc=services.kc, complete=complete_coverage)
        except Exception as error:   # reported, never hidden
            failures.append(f"balance check could not run: {type(error).__name__}: "
                            f"{str(error)[:500]}")
    for k, payload in zip(have, payloads):
        rows[k]["payload"] = payload
        rows[k]["status"] = payload["status"]
        rows[k]["refusal"] = (payload.get("refusal") or {}).get("status")
    refused_pairs = [r for r in cycle_table if r.get("status") == "refused:cycle"]
    failures += te.cycle_gates(cycle_table, payloads, complete_coverage)
    if cycle_table and len(refused_pairs) > CYCLE_REFUSED_FRACTION * len(cycle_table):
        failures.append(f"more than {CYCLE_REFUSED_FRACTION:.0%} of the loaded reciprocal pairs are "
                        f"refused by the cycle check: {len(refused_pairs)} of {len(cycle_table)}")
    calibration = []
    for row in rows:
        p = row.get("payload") or {}
        for where, cal in [("record", p.get("calibration"))] + \
                [(f"template {t.get('key', '')[:12]}", t.get("calibration"))
                 for t in p.get("templates", []) if t.get("status") == "emitted"]:
            for c in cal or []:
                calibration.append(c)
                if not c.get("ok"):
                    failures.append(f"calibration {row['event_id']} {where} at {c.get('T')} K: {c}")
    checked = [r for r in cycle_table if r.get("status") != "not_evaluable" and r.get("rows")]
    balance = {
        "calibration_rows": len(calibration),
        "calibration_max_relative_error": max([c["relative_error"] for c in calibration
                                               if c.get("relative_error") is not None], default=None),
        "calibration_ok": all(c.get("ok") for c in calibration),
        "worst_ratio_error": max([r["max_ratio_error"] for r in checked], default=None),
        "worst_grid_error": max([r.get("max_grid_error", 0.0) for r in cycle_table], default=0.0),
        "worst_offgrid_error": max([r.get("max_offgrid_error", 0.0) for r in cycle_table],
                                   default=0.0),
        "self_inverse_templates": sum(1 for p in payloads for t in p.get("templates", [])
                                      if t.get("status") == "emitted" and t.get("self_inverse"))}
    out_payloads = stripped_payloads(payloads)
    expected = sum(1 for p in payloads for t in p.get("templates", [])
                   if t.get("status") == "emitted")
    temperatures = [float(t) for t in first.get("temperatures", [])]
    round_trip = {"temperatures": len(temperatures), "templates": expected, "ok": True}
    for T in temperatures:
        try:
            loaded = len(te.loaded_templates(out_payloads, T))
        except Exception as error:
            round_trip["ok"] = False
            failures.append(f"round trip at {T} K: {type(error).__name__}: {str(error)[:300]}")
            break
        if loaded != expected:
            round_trip["ok"] = False
            failures.append(f"round trip at {T} K loads {loaded} templates, rows emitted {expected}")
            break
    if not temperatures:
        round_trip["ok"] = False
        failures.append("round trip: no provenance temperatures")
    completeness = defaultdict(Counter)
    refusals = []
    for row in rows:
        completeness[f"{row.get('family')}|{row.get('orientation')}"][row.get("status")] += 1
        if row.get("status") != "emitted":
            full = (row.get("payload") or {}).get("record_event_id") or row["event_id"]
            refusals.append({"record_event_id": full, "family": row.get("family"),
                             "orientation": row.get("orientation"), "status": row.get("status"),
                             "refusal": row.get("refusal"), "reason": _refusal_reason(row)})
    completeness = {k: {s: v.get(s, 0) for s in ("emitted", "partial", "refused", "error")}
                    for k, v in sorted(completeness.items())}
    failures = list(dict.fromkeys(failures))
    return {
        "schema": SCHEMA,
        "provenance": {k: first.get(k) for k in (
            "artifact_id", "artifact_sha256", "artifact_records", "database_path", "database_head",
            "code_sha", "roster_sha256", "families", "temperatures", "systems")},
        "payloads": out_payloads,
        "refusals": refusals,
        "gates": {
            "collisions": {"keys": len(collisions),
                           "merged": len(collisions) - len(refused_collisions),
                           "refused": len(refused_collisions),
                           "copies": sum(c["copies"] for c in collisions),
                           "tolerance_log10": te.LOCALITY_TOLERANCE_LOG10},
            "cycles": {"pairs": len(cycle_table), "refused": len(refused_pairs),
                       "not_evaluable": len(te.cycle_not_evaluable(cycle_table)),
                       "by_status": dict(Counter(r.get("status") for r in cycle_table))},
            "balance": balance,
            "round_trip": round_trip,
            "coverage": {"complete": complete_coverage, "rows": len(rows),
                         "assigned": len(assigned)},
            "midpoint": {"evaluated": False,
                         "reason": "the primary-weight midpoint check regenerates rates from the "
                                   "artifact records and is not part of assembly"},
            "failures": failures},
        "completeness": completeness,
        "complete": not failures,
    }


def run_assemble(args):
    services = None
    catalog = assemble(args.rows, database=args.database, services=services)
    Path(args.out).write_text(dumps(catalog))
    print("CATALOG", args.out, "complete" if catalog["complete"] else "INCOMPLETE",
          len(catalog["payloads"]), "payloads", flush=True)
    for failure in catalog["gates"]["failures"]:
        print("GATE", failure, flush=True)
    return 0 if catalog["complete"] else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "assemble":
        parser = argparse.ArgumentParser(prog="python -m rmgpy.kmc.emit assemble")
        parser.add_argument("rows", nargs="+", help="rows files (shards) of one artifact")
        parser.add_argument("--out", required=True)
        parser.add_argument("--database", default=None,
                            help="RMG-database path for the balance check (default: the "
                                 "shard headers')")
        return run_assemble(parser.parse_args(argv[1:]))
    parser = argparse.ArgumentParser(prog="python -m rmgpy.kmc.emit")
    parser.add_argument("artifact")
    parser.add_argument("roster")
    parser.add_argument("--database", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--chunk", default=None, help="i/n: records i::n of the sorted list")
    parser.add_argument("--family", action="append", default=None,
                        help="only these families' records (repeatable; a partial catalog)")
    parser.add_argument("--limit", type=int, default=None, help="at most N records (a partial catalog)")
    parser.add_argument("--no-contexts", action="store_true")
    return run_emit(parser.parse_args(argv[1:] if argv and argv[0] == "run" else argv))


if __name__ == "__main__":
    sys.exit(main())
