"""I-138 acceptance for rmgpy/kmc/template_emission.py (not collected).

One reference per run (its compile artifact from the I-127 probe recipe,
``/home/alon/runs/i138/compile/reference_compile.py``):

* ``dp8_closed`` (homolysis; the unequal-fragment recombination records),
  ``dp8_headrad`` (unimolecular families), ``dp4_closed+pe`` (H abstraction),
  ``dp4rad+pe`` (disproportionation, recombination of unequal radicals),
  ``pe2_identical`` (two 1-phenylethyl-ended chains: identical radical
  recombination), ``c16`` (the alkane provider on C16).

Per record: the upstream status (record status, record certificate from
``certificates.certify_records(required=True)``), the provider members of
its channel (every member case of the reference specs, any family direction
and subset, whose template-side species, template and other-side species
equal the record's),
the binding (``recertify`` on the member witness), ``windows.run_member``
with contexts, and ``template_emission.emit_record`` with the production
rate path (``ProductionRates``) and the independent ``rmg_moves`` count.
Per family and direction whose reactants are the reference species: the
matcher union over the emitted templates against ``rmg_moves`` on the whole
molecules (``template_emission.classify_reference``), RMG-only moves traced
to their record's status, the per-key overlap multiplicity.  ``--dp12`` (the
dp8 reference): the homolysis set comparison on PS DP12.

Every record row goes to ``rows.<ref>.jsonl``; ``summary.<ref>.json`` holds the
tables, timings, the projection of a sampled run and the gate failures
(non-empty: exit 1).

``ladder:<family>`` (item 4): the family's ``windows.ladder_plan`` (the pair
roster of ``windows.pair_roster``; the generated class roster is the other
branch's: ``--roster classes.<family>.pkl`` when given, else none), every
member of every item at DP <= L as a whole-species template (production
rate on the species itself, both counts, no locality test), and every DP
L + 1 control member checked to be covered by a local template.

Usage: python template_emission_acceptance.py OUTDIR REF [--artifact PATH]
       [--sample K] [--chunk i/n] [--contexts 0|1] [--dp12] [--limit N]
       python template_emission_acceptance.py OUTDIR ladder:<family> [--sample K]
       [--limit N] [--roster PATH]
"""

import argparse
import json
import os
import sys
import time
import traceback
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

DATABASE = Path(os.environ.get("RMG_KMC_I138_DATABASE", "/home/alon/runs/phase2b/database"))
COMPILES = Path(os.environ.get("RMG_KMC_I138_COMPILES", "/home/alon/runs/i138/compile"))
_S = {}


def references():
    from rmgpy.kmc import boundaries as b
    S = b.PsChainSpec
    pe = S(1, tail_end="benzylic_radical")
    import interior_radical_specs
    interior_radical_specs.install()
    return {
        "dp8_closed": [S(8)],
        "dp8_headrad": [S(8, head_end="primary_radical")],
        "dp4_closed+pe": [S(4), pe],
        "dp4rad+pe": [S(4, tail_end="benzylic_radical"), pe],
        "pe2_identical": [S(2, tail_end="benzylic_radical"), S(2, tail_end="benzylic_radical")],
        "c16": [("alkane", 16)],
        # I-150: linear C25 alkyl radical, radical at C11 (0-based 10); compiled by
        # /home/alon/runs/i150/ref_compile.py c25r11 (intra_H_migration)
        "c25r11": [("alkane", 25, None, None, 10)],
    }


def load_database(families, path=DATABASE):
    """The compile's own load (``rmgpy.kmc.compile.compile_artifact``)."""
    from rmgpy.data.rmg import RMGDatabase
    from rmgpy.kmc.compiler import prepare_rate_rules
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


def reference_cases(certifier, fam, specs, species):
    """Every member case of every family direction on every subset of the
    reference specs of the direction's molecularity, keyed by (family,
    template-side species, template, other-side species)."""
    from itertools import combinations
    from rmgpy.kmc import windows as win
    out = defaultdict(list)
    stats = Counter()
    for name, family in sorted(fam.items()):
        directions = ["forward"] + (["reverse"] if family.reverse_template is not None else [])
        for direction in directions:
            n = molecularity(family, direction)
            for subset in sorted(set(combinations(range(len(specs)), n))):
                chosen = [specs[i] for i in subset]
                cases = win.member_cases(certifier, family, win.builds_of(chosen), direction,
                                         label=f"{name}|{direction}|{subset}")
                stats[f"{name}|{direction}"] += len(cases)
                for case in cases:
                    mols = win.template_molecules(case.witness, case.direction)
                    other = win.template_molecules(
                        case.witness, "reverse" if case.direction == "forward" else "forward")
                    key = (name, tuple(sorted(species(m) for m in mols)),
                           ";".join(case.pathway), tuple(sorted(species(m) for m in other)))
                    out[key].append((case, chosen))
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


def init(ref, artifact, only=None):
    from rmgpy.kmc import certificates as cm
    from rmgpy.kmc import channels as chm
    from rmgpy.kmc import template_emission as te
    from rmgpy.kmc import windows as win
    from rmgpy.kmc.event_record import EventRecord
    from rmgpy.kmc.reference_thermo import GasPhaseRMGReferenceThermo
    art = json.loads(Path(artifact).read_text())
    if only:
        # one family's records (certified among themselves: a reversed
        # record's rate owner is in its own family)
        art["records"] = [r for r in art["records"] if r["family"] in only]
    families = sorted({r["family"] for r in art["records"]})
    t = time.time()
    db = load_database(families)
    fam = db.kinetics.families
    _S.update(db=db, fam=fam, load_s=round(time.time() - t, 1))
    certifier = cm.Certifier(fam, max_channel_sets=8, max_entries=1024)
    t = time.time()
    records = [EventRecord.from_dict(r) for r in art["records"]]
    certified, store, replacements = cm.certify_records(records, fam, certifier, required=True)
    certified = [r.to_dict() if hasattr(r, "to_dict") else r for r in certified]
    _S.update(certify_s=round(time.time() - t, 1), records=certified, store=store,
              by_id={r["event_id"]: r for r in certified}, certifier=certifier)
    species = chm.SpeciesKeys()
    t = time.time()
    cases, stats = reference_cases(certifier, fam, references()[ref], species)
    _S.update(cases=cases, case_stats=stats, cases_s=round(time.time() - t, 1), species=species,
              bounds={n: win.span_bound(win.family_constants(f)) for n, f in fam.items()},
              counter=te.RmgCounter(fam),
              kc=te.ThermoKc(GasPhaseRMGReferenceThermo(db.thermo, "i147-reference")),
              rates=te.ProductionRates(db, families, database_path=DATABASE,
                                       database_sha="i138-reference"))


def upstream_of(record):
    """``open`` or the upstream refusal status of a record (cheap: the
    certificate store and the record status)."""
    from rmgpy.kmc import template_emission as te
    rc = _S["store"].get((record.get("rate_certificate") or {}).get("sha256"))
    why = te.upstream_refusal(record, rc, _S["by_id"])
    return "open" if why is None else te.status_of(why)


def upstream_counts(records):
    """Per family and orientation, every record's upstream status."""
    out = defaultdict(Counter)
    for r in records:
        out[f"{r['family']}|{r['orientation']}"][upstream_of(r)] += 1
    return {k: dict(v) for k, v in sorted(out.items())}


def run_record(record, contexts=True):
    from rmgpy.kmc import certificates as cm
    from rmgpy.kmc import template_emission as te
    from rmgpy.kmc import windows as win
    certifier, fam = _S["certifier"], _S["fam"]
    start, cpu = time.time(), time.process_time()
    rc = _S["store"].get((record.get("rate_certificate") or {}).get("sha256"))
    row = {"event_id": record["event_id"][:24], "family": record["family"],
           "orientation": record["orientation"], "template": record["template"][:100],
           "record_status": record["status"], "rate_kind": (record.get("rate_source") or {})
           .get("kind"), "certified": (rc or {}).get("certified")}
    family = fam[record["family"]]
    why = te.upstream_refusal(record, rc, _S["by_id"])
    members, member_rows = [], []
    try:
        if why is None:
            key = record_key(certifier, record, _S["species"])
            found = _S["cases"].get(key, [])
            row["members"] = len(found)
            if found:
                owner = cm._rate_owner(record, _S["by_id"])
                row["bind"] = bind(certifier, record, rc, found[0][0], owner)
                if not row["bind"].get("same_selection"):
                    why = f"refused:recertify: {row['bind']}"[:300]
            if why is None:
                for case, specs in found:
                    out = []
                    keep = []
                    res = win.run_member(certifier, family, case,
                                         bound=_S["bounds"][record["family"]],
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
                                     counter=_S["counter"], rates=_S["rates"],
                                     record_certificate=rc, by_id=_S["by_id"],
                                     check_upstream=False, kc=_S["kc"])
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


# --------------------------------------------------------------------------
# reference tables
# --------------------------------------------------------------------------
def reference_molecules(ref):
    from rmgpy.kmc import boundaries as b
    from rmgpy.kmc import windows as win
    mols, _parts = b.witness_reactants(win.builds_of(references()[ref]))
    return list(mols)


def comparison(ref, rows, molecules=None, label=None, owners_default=None):
    """Per family and direction: matcher union over the emitted templates of
    the records whose reactants are the reference species, against rmg_moves
    on the whole molecules."""
    from rmgpy.kmc import template_emission as te
    fam, counter, species = _S["fam"], _S["counter"], _S["species"]
    molecules = molecules or reference_molecules(ref)
    ref_species = tuple(sorted(species(m) for m in molecules))
    index = defaultdict(list)
    templates = defaultdict(list)
    loaded = set()
    # a merged template is carried by the record that names it: loading a
    # subset of records must not lose a template whose owner is outside it
    owned = {(row["payload"].get("record_event_id"), t["key"]): t
             for row in rows for t in (row.get("payload") or {}).get("templates", [])
             if t.get("status") == "emitted"}
    for row in rows:
        payload = row.get("payload") or {}
        rec = next((r for r in _S["records"] if r["event_id"].startswith(row["event_id"])), None)
        if rec is None:
            continue
        family = fam[rec["family"]]
        direction = te.record_direction(rec, family)
        reac = counter.species_of_rows(te.record_graph_rows(rec, "reactant_graphs"))
        prod = counter.species_of_rows(te.record_graph_rows(rec, "product_graphs"))
        status = row.get("status")
        refusal = row.get("refusal") or ""
        index[(rec["family"], direction, reac, prod)].append(
            (rec["template"], status, refusal, (payload.get("refusal") or {}).get("reason", ""),
             payload.get("weight_source")))
        if not set(reac) <= set(ref_species):
            continue
        for t in payload.get("templates", []):
            if t.get("status") == "merged":
                into = t["merged_into"]
                t = owned.get((into["record"], into["key"]))
                if t is None:
                    raise te.EmissionError(f"orphaned merged template of {rec['event_id'][:16]}")
                owner_id = into["record"]
            elif t.get("status") == "emitted":
                owner_id = rec["event_id"]
            else:
                continue
            name = f"{owner_id[4:16]}:{t['key'][:12]}"
            if name not in loaded:
                loaded.add(name)
                d = t.get("derived")
                primary = owned.get((d["primary"]["record"], d["primary"]["key"])) if d else None
                templates[(rec["family"], direction)].append(
                    te.matcher_from_payload(t, 700.0, name=name, primary=primary))
    out = {}
    for (name, direction), tlist in sorted(templates.items()):
        family = fam[name]

        def owners(mv, name=name, direction=direction):
            ids = set(mv.product_rows())
            reac = tuple(sorted(species(m) for m in molecules
                                if {x.id for x in m.atoms} <= ids))
            prod = counter.species_of_rows(mv.product_rows())
            hits = [x for xs in [index.get((name, direction, reac, prod), [])] for x in xs]
            tmpl = ";".join(mv.template or ())
            same = [x for x in hits if x[0] == tmpl] or hits
            if not same:
                return owners_default or ("unowned", "no record of the channel")
            status, refusal, reason, source = same[0][1], same[0][2], same[0][3], same[0][4]
            if status in ("emitted", "partial") and source == "emitted_witness":
                return (te.BEYOND_PROXY, "the record's templates need the chain beyond the "
                                         "reference's artificial ends")
            if status in ("emitted", "partial"):
                return ("emitted", "the record's templates should produce it")
            return (refusal or f"refused:{status}", reason or "")
        table = te.classify_reference(tlist, molecules, family=family, directions=(direction,),
                                      counter=counter, owners=owners)
        table["templates"] = len(tlist)
        try:
            te.check_class_overlap(tlist, te.melt_state(
                {a: r for m in molecules for a, r in te.rm.molecule_rows(m).items()}))
            table["class_overlaps"] = 0
        except te.ClassOverlapError as error:
            found = te.class_overlaps(tlist, te.melt_state(
                {a: r for m in molecules for a, r in te.rm.molecule_rows(m).items()}))
            table["class_overlaps"] = len(found)
            table.setdefault("failures", []).append(f"ClassOverlapError: {error}")
        out[f"{name}|{direction}"] = table
    # families with records on the reference species but no emitted template
    for key, xs in index.items():
        name, direction, reac, _prod = key
        if set(reac) <= set(ref_species) and f"{name}|{direction}" not in out:
            out.setdefault(f"{name}|{direction}", {"templates": 0, "failures": [
                f"{name} {direction}: every record of the family was refused "
                f"({dict(Counter(x[2] for k2, v in index.items() if k2[:2] == (name, direction) for x in v))})"]})
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("outdir")
    parser.add_argument("ref")
    parser.add_argument("--artifact", default=None)
    parser.add_argument("--sample", type=int, default=1)
    parser.add_argument("--chunk", default=None, help="i/n: records i::n of the sampled list")
    parser.add_argument("--contexts", type=int, default=1)
    parser.add_argument("--dp12", action="store_true")
    parser.add_argument("--limit", type=int, default=None, help="at most N records (timing)")
    parser.add_argument("--roster", default=None)
    parser.add_argument("--event", action="append", default=None,
                        help="only the records whose event id starts with this (repeatable)")
    parser.add_argument("--family", action="append", default=None,
                        help="only these families' records (repeatable)")
    parser.add_argument("--controls-only", action="store_true",
                        help="ladder: only the L+1 control items")
    parser.add_argument("--open-only", action="store_true",
                        help="sample only the records no upstream status refuses")
    args = parser.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    if args.ref.startswith("ladder:"):
        return ladder(out, args)
    if args.ref not in references():
        raise SystemExit(f"unknown reference {args.ref}")
    artifact = args.artifact
    if artifact is None:
        base = COMPILES / args.ref
        artifact = str(base / (base / "artifact.path").read_text().strip())
    t0 = time.time()
    init(args.ref, artifact, only=args.family)
    records = sorted(_S["records"], key=lambda r: r["event_id"])
    _S["upstream"] = upstream_counts(records)
    print("UPSTREAM", json.dumps(_S["upstream"]), flush=True)
    pool = [r for r in records if upstream_of(r) == "open"] if args.open_only else records
    if args.event:
        pool = [r for r in pool if any(r["event_id"].startswith(e) for e in args.event)]
    picked = pool[::args.sample]
    if args.chunk:
        k, n = (int(x) for x in args.chunk.split("/"))
        picked = picked[k::n]
    if args.limit:
        picked = picked[:args.limit]
    print("INIT", {k: _S.get(k) for k in ("load_s", "certify_s", "cases_s", "case_stats")},
          len(_S["records"]), "records", "selected", len(picked),
          [r["event_id"][:16] for r in picked], flush=True)
    for problem in selection_problems(args, records, picked):
        print("SELECTION", problem, flush=True)
    rows = []
    tag = args.ref + (f".{args.chunk.replace('/', 'of')}" if args.chunk else "") + \
        ("." + "+".join(args.family) if args.family else "")
    with open(out / f"rows.{tag}.jsonl", "w") as log:
        for k, record in enumerate(picked):
            row = run_record(record, contexts=bool(args.contexts))
            rows.append(row)
            log.write(json.dumps(row, default=str) + "\n")
            log.flush()
            print(f"[{time.time() - t0:6.0f}s] {k + 1}/{len(picked)} {row['event_id'][:16]} "
                  f"{row['family']} {row['orientation']} -> {row['status']} "
                  f"{row.get('refusal') or ''} t={row['seconds']}s rss={row['rss_mb']}MB",
                  flush=True)
    summary = summarize(args, rows, records, picked, time.time() - t0)
    (out / f"summary.{tag}.json").write_text(json.dumps(summary, indent=1, default=str))
    with open(out / f"rows.{tag}.resolved.jsonl", "w") as log:   # after the collision gate
        for row in rows:
            log.write(json.dumps(row, default=str) + "\n")
    print("SUMMARY", json.dumps({k: v for k, v in summary.items() if k != "comparison"},
                                default=str)[:6000], flush=True)
    print("GATES", len(summary["gates"]), json.dumps(summary["gates"][:40], default=str),
          flush=True)
    return 1 if summary["gates"] else 0


def collide(rows):
    """The catalog collision gate (``template_emission.resolve_collisions``)
    over every row's payload: rows carry the resolved payloads (status and
    refusal follow); returns the collision table."""
    from rmgpy.kmc import template_emission as te
    have = [k for k, r in enumerate(rows) if r.get("payload")]
    merged, table = te.resolve_collisions([rows[k]["payload"] for k in have])
    for k, payload in zip(have, merged):
        rows[k]["payload"] = payload
        rows[k]["status"] = payload["status"]
        rows[k]["refusal"] = (payload.get("refusal") or {}).get("status")
    return table


def interpolation_table(payloads):
    """Worst leave-one-out interpolation error of the stored weight tables
    (``template_emission.interpolation_error``) over the emitted templates.
    REPORTED, NOT GATED: each interior point is predicted from its two grid
    neighbours (twice the spacing), and with three or fewer interior points
    that measures the coarseness of the grid, not the runtime error; the gate
    is ``midpoint_table``."""
    from rmgpy.kmc import template_emission as te
    worst, at, who, n = 0.0, None, None, 0
    for p in payloads:
        for t in p.get("templates", []):
            if t.get("status") == "emitted" and t.get("weight") and not t.get("derived"):
                err, T = te.interpolation_error(t["weight"]["T"], t["weight"]["w"])
                n += 1
                if err > worst:
                    worst, at, who = err, T, f"{p.get('record_event_id', '')[:16]}:{t['key'][:12]}"
    return {"templates": n, "worst": worst, "worst_T": at, "worst_template": who,
            "method": "each interior grid point predicted from its two grid neighbours "
                      "(2x the spacing; records carry only the grid k table); reported only"}


def midpoint_table(payloads):
    """Worst error of the runtime interpolation at the grid midpoints against
    the record's OWN kinetics there (gated at 1%): the weight table
    interpolated linearly in 1/T at (T_i + T_i+1) / 2 against
    ``ssa * k / T_ch`` with ``k`` of the record's own rate regenerated at the
    midpoints (``ProductionRates`` on the midpoint grid, the emitted witness's
    channel).  Primary templates with a stored emitted witness only; the rest
    (derived, merged, no witness, no regenerated rate) are counted as
    ``unavailable`` with the reason."""
    from rmgpy.kmc import template_emission as te
    worst, at, who, n = 0.0, None, None, 0
    unavailable = Counter()
    unverified = []      # loaded primaries whose midpoint check could not be evaluated: gated
    cache = {}
    for p in payloads:
        rid = p.get("record_event_id")
        record = _S["by_id"].get(rid)
        for t in p.get("templates", []):
            if t.get("status") != "emitted" or not t.get("weight"):
                continue
            if t.get("derived"):
                unavailable["derived"] += 1
                continue
            wit = ((t.get("maps") or {}).get("witnesses") or {}).get("emitted")
            lrow = next((r for r in (t.get("locality") or {}).get("rows", [])
                         if r.get("label") == "emitted" and r.get("status") == "ok"), None)
            name = f"{(rid or '')[:16]}:{t['key'][:12]}"
            if record is None or wit is None or lrow is None:
                unavailable["no emitted witness"] += 1
                unverified.append((name, "no emitted witness or record to regenerate on"))
                continue
            temps = t["weight"]["T"]
            mids = tuple((a + b) / 2.0 for a, b in zip(temps, temps[1:]))
            if mids not in cache:
                cache[mids] = te.ProductionRates(
                    _S["db"], sorted({r["family"] for r in _S["records"]}), database_path=DATABASE,
                    database_sha="i138-reference", temperature_grid=mids)
            r_rows, p_rows = te.rows_from_json(wit["r"]), te.rows_from_json(wit["p"])
            gen, prods = (p_rows, r_rows) if t.get("member_reverse") else (r_rows, p_rows)
            # the regenerated k is the EMITTED witness's, so its own count (the one the
            # weight was calibrated with) is the denominator, not the record's
            tch = (lrow.get("counts") or {}).get("matcher")
            try:
                result = cache[mids](te.RateRequest(record["family"], record["template"],
                                                    record["orientation"], gen, prods))
                if result.problems or tuple(result.temperatures) != mids:
                    raise te.EmissionError("; ".join(result.problems) or "grid mismatch")
                truth = te.weight_table(mids, result.k, result.ssa_multiplier, tch)
            except Exception as error:   # reported, never hidden
                reason = f"{type(error).__name__}: {str(error)[:120]}"
                unavailable[reason] += 1
                unverified.append((name, reason))
                continue
            n += 1
            for T, w in zip(mids, truth):
                err = abs(te._interpolate(temps, t["weight"]["w"], T) / w - 1.0)
                if err > worst:
                    worst, at, who = err, T, f"{rid[:16]}:{t['key'][:12]}"
    return {"templates": n, "worst": worst, "worst_T": at, "worst_template": who,
            "unavailable": dict(unavailable), "unverified": unverified,
            "method": "weight table interpolated at the grid midpoints against ssa*k/T_ch of the "
                      "record's own kinetics regenerated there"}


def midpoint_gates(mid):
    """The gate failures of ``midpoint_table``: an error above 1%, and every
    loaded primary template whose midpoint check could not be evaluated (a
    gate that evaluates nothing must not pass)."""
    out = []
    if mid["worst"] > 0.01:
        out.append(f"primary-weight midpoint error {mid['worst']:.3g} > 1% at "
                   f"{mid['worst_T']} K ({mid['worst_template']})")
    for name, reason in mid.get("unverified", ()):
        out.append(f"primary-weight midpoint check not evaluable for {name}: {reason}")
    return out


def is_complete(args):
    """Whether the run covers the whole catalog (the cycle check then refuses an
    absent partner; a subset reports it ``not_evaluable``)."""
    return args.sample == 1 and not args.chunk and not args.limit and not args.open_only \
        and not args.event


def cycle_gate(rows, complete):
    """The equilibrium-cycle check of the resolved, loaded catalog against
    independent thermo (``template_emission.check_cycles``); refusals follow
    into the rows, and every failure is returned as a gate."""
    from rmgpy.kmc import template_emission as te
    have = [k for k, r in enumerate(rows) if r.get("payload")]
    payloads = [rows[k]["payload"] for k in have]
    payloads = te.derive_weights(payloads, counter=_S["counter"])
    out, table = te.check_cycles(payloads, counter=_S["counter"], kc=_S["kc"],
                                 complete=complete)
    for k, payload in zip(have, out):
        rows[k]["payload"] = payload
        rows[k]["status"] = payload["status"]
        rows[k]["refusal"] = (payload.get("refusal") or {}).get("status")
    checked = [r for r in table if r["status"] != "not_evaluable"]
    worst_row = max([r for r in checked if r.get("rows")], key=lambda r: r["max_ratio_error"],
                    default=None)
    loaded_pairs = len(table)
    refused = [r for r in table if r["status"] == "refused:cycle"]
    interp = interpolation_table(payloads)
    mid = midpoint_table(payloads)
    self_inverse = sum(1 for p in out for t in p.get("templates", [])
                       if t.get("status") == "emitted" and t.get("self_inverse"))
    return {"pairs": loaded_pairs, "refused": len(refused), "interpolation": interp, "midpoint": mid,
            "self_inverse_templates": self_inverse,
            "not_evaluable": te.cycle_not_evaluable(table),
            "by_status": dict(Counter(r["status"] for r in table)),
            "worst_ratio_error": worst_row["max_ratio_error"] if worst_row else None,
            "worst_pair": {"forward": worst_row["forward"], "reverse": worst_row["reverse"],
                           **(worst_row.get("worst") or {})} if worst_row else None,
            "worst_offgrid_error": max([r.get("max_offgrid_error", 0.0) for r in table] or [0.0]),
            "worst_grid_error": max([r.get("max_grid_error", 0.0) for r in table] or [0.0]),
            "gates": te.cycle_gates(table, out, complete)
            + midpoint_gates(mid),
            "table": table}


def orbit_table(rows):
    """Templates refused ``refused:non_transferable`` and emitted, per row."""
    refused = Counter()
    emitted = 0
    for row in rows:
        for t in (row.get("payload") or {}).get("templates", []):
            if t.get("status") == "emitted":
                emitted += 1
            elif (t.get("refusal") or {}).get("status") == "refused:non_transferable":
                refused["refused:non_transferable"] += 1
    return {"emitted_templates": emitted, "refused": dict(refused)}


def selection_problems(args, records, picked):
    """Gate failures of the record selection: an ``--event`` prefix that matches
    no record (event ids differ between compiles), or a selection that processes
    nothing (a run that checks nothing must not pass)."""
    out = [f"--event prefix {e} matches no record" for e in (getattr(args, "event", None) or ())
           if not any(r["event_id"].startswith(e) for r in records)]
    if not picked:
        out.append("the selection processed 0 records (a vacuous run)")
    return out


def summarize(args, rows, records, picked, wall):
    from rmgpy.kmc import template_emission as te
    gates = selection_problems(args, records, picked)
    by_family = defaultdict(Counter)
    calibration_rows = []
    try:
        collisions = collide(rows)
    except te.EmissionError as error:      # a split reciprocal pair, a bad weight
        gates.append(f"catalog consistency: {type(error).__name__}: {error}"[:600])
        collisions = []
    refused_collisions = [c for c in collisions if c["disposition"] != "merged"]
    for c in refused_collisions:
        gates.append(f"{c['disposition']} {c['key'][:16]}: {c['copies']} copies, "
                     f"log10 diameter {c['log10_diameter']} (records {c['records']})")
    complete = is_complete(args)
    cycles = cycle_gate(rows, complete)
    gates += cycles["gates"]
    if cycles["pairs"] and cycles["refused"] > 0.05 * cycles["pairs"]:
        gates.append(f"more than 5% of the loaded reciprocal pairs are refused by the cycle "
                     f"check: {cycles['refused']} of {cycles['pairs']}")
    for row in rows:
        fam_key = f"{row['family']}|{row['orientation']}"
        by_family[fam_key][row.get("refusal") or row["status"]] += 1
        p = row.get("payload") or {}
        for t in p.get("templates", []):
            by_family[fam_key]["templates:" + t["status"]] += 1
            if t["status"] == "refused":
                by_family[fam_key]["template " + t["refusal"]["status"]] += 1
        cals = [("record", p.get("calibration"))] + \
            [(f"template {t.get('key', '')[:12]}", t.get("calibration"))
             for t in p.get("templates", []) if t.get("status") == "emitted"]
        for where, cal in cals:
            for c in cal or []:
                calibration_rows.append((row["event_id"], c["T"], c.get("relative_error"),
                                         c["ok"], where))
                if not c["ok"]:
                    gates.append(f"calibration {row['event_id']} {where} at {c['T']} K: {c}")
        by_family[fam_key]["weight_source:" + str(p.get("weight_source"))] += 1
        if row["status"] == "error":
            gates.append(f"error {row['event_id']}: {row.get('reason', '')[:300]}")
    out = {"collisions": {"keys": len(collisions),
                          "merged": len(collisions) - len(refused_collisions),
                          "refused": len(refused_collisions),
                          "copies": sum(c["copies"] for c in collisions),
                          "max_log10_diameter": max([c["log10_diameter"] or 0.0
                                                     for c in collisions] or [0.0]),
                          "tolerance_log10": te.LOCALITY_TOLERANCE_LOG10, "table": collisions},
           "cycles": cycles,
           "orbit": orbit_table(rows),
           "ref": args.ref, "families": args.family or "all", "records": len(records),
           "upstream_all": _S.get("upstream"), "open_only": args.open_only,
           "processed": len(picked),
           "sample": args.sample, "chunk": args.chunk, "contexts": args.contexts,
           "load_s": _S.get("load_s"), "certify_s": _S.get("certify_s"),
           "cases_s": _S.get("cases_s"), "case_stats": _S.get("case_stats"),
           "per_family": {k: dict(v) for k, v in sorted(by_family.items())},
           "calibration": {"rows": len(calibration_rows),
                           "max_relative_error": max([c[2] for c in calibration_rows
                                                      if c[2] is not None] or [None])
                           if calibration_rows else None,
                           "all_ok": all(c[3] for c in calibration_rows)},
           "production": {"calls": _S["rates"].calls, "generated": _S["rates"].generated},
           "wall_s": round(wall), "record_seconds": round(sum(r["seconds"] for r in rows), 1),
           "record_cpu": round(sum(r["cpu"] for r in rows), 1)}
    if rows:
        per = out["record_seconds"] / len(rows)
        n_open = sum(v.get("open", 0) for v in (_S.get("upstream") or {}).values())
        opened = [r for r in rows if upstream_of(next(
            x for x in records if x["event_id"].startswith(r["event_id"]))) == "open"]
        per_open = sum(r["seconds"] for r in opened) / len(opened) if opened else None
        out["projection"] = {"per_record_s": round(per, 1),
                             "full_record_s": round(per * len(records)),
                             "open_records": n_open,
                             "per_open_record_s": round(per_open, 1) if per_open else None,
                             "full_open_s": round(per_open * n_open) if per_open else None,
                             "note": "sampled records are a deterministic stride of the "
                                     "event-id-sorted list (of the open ones with "
                                     "--open-only); upstream-refused records cost ~0 s"}
    emitted = sum(1 for r in rows if r["status"] in ("emitted", "partial"))
    if rows and not emitted:
        gates.append("every processed record was refused (a vacuous run)")
    if complete:
        out["comparison"] = comparison(args.ref, rows)
        for key, table in out["comparison"].items():
            for f in table.get("failures", []):
                gates.append(f"{args.ref} {key}: {f}")
            if table.get("templates", 0) and not table.get("matched"):
                gates.append(f"{args.ref} {key}: no compared move matched (vacuous)")
        if args.dp12 and args.ref == "dp8_closed":
            out["dp12"] = dp12(rows)
    else:
        out["comparison"] = "skipped: a sampled or chunked run (needs every record)"
    out["gates"] = gates
    return out


def dp12(rows):
    """The DP12 homolysis set comparison with the DP8 templates (RMG moves on
    PS DP12, no DP12 compile): matched / RMG-only / matcher-only."""
    from rmgpy.kmc import boundaries as b
    from rmgpy.kmc import windows as win
    t = time.time()
    mols, _parts = b.witness_reactants(win.builds_of([b.PsChainSpec(12)]))
    table = comparison("dp8_closed", rows, molecules=list(mols),
                       owners_default=("refused:no_dp12_record", "set comparison only"))
    for key, v in table.items():
        v["note"] = "RMG-only moves have no DP12 record: reported, not traced"
    return {"seconds": round(time.time() - t, 1), "tables": table}


# --------------------------------------------------------------------------
# item 4: the ladder
# --------------------------------------------------------------------------
def ladder(out, args):
    from dataclasses import replace
    from rmgpy.kmc import certificates as cm
    from rmgpy.kmc import template_emission as te
    from rmgpy.kmc import windows as win
    name = args.ref.split(":", 1)[1]
    t0 = time.time()
    db = load_database([name])
    fam = db.kinetics.families
    family = fam[name]
    certifier = cm.Certifier(fam, max_channel_sets=8, max_entries=1024)
    counter = te.RmgCounter(fam)
    rates = te.ProductionRates(db, [name], database_path=DATABASE,
                               database_sha="i138-reference")
    constants = win.family_constants(family)
    classes = ()
    if args.roster:
        import pickle
        classes = tuple(pickle.loads(Path(args.roster).read_bytes()))
    roster = win.pair_roster(certifier, family)
    constants = replace(constants, pair_kinds=roster, classes=classes, merge=None)
    bound = win.span_bound(constants)
    plan = win.ladder_plan(constants, bound)
    molecularity = dict(constants.directions)
    print("LADDER", name, "L", bound.L, "items", len(plan.items), "roster", roster,
          "classes", len(classes), round(time.time() - t0, 1), "s", flush=True)
    items = [x for x in plan.items if x.control or not args.controls_only][::args.sample]
    if args.limit:
        items = items[:args.limit]
    rows, gates = [], []
    with open(out / f"rows.ladder_{name}.jsonl", "w") as log:
        for k, item in enumerate(items):
            directions = [d for d, n in sorted(molecularity.items()) if n == 1] \
                if len(item.specs) == 1 else list(item.directions)
            for direction in directions:
                start = time.time()
                try:
                    cases = win.member_cases(certifier, family, win.builds_of(item.specs),
                                             direction, label=item.label)
                except Exception as error:   # reported, never hidden
                    rows.append({"item": item.label, "direction": direction, "status": "error",
                                 "reason": f"{type(error).__name__}: {error}"[:300]})
                    continue
                for case in cases:
                    start = time.time()
                    keep = []
                    res = win.run_member(certifier, family, case, bound=bound, contexts=False,
                                         specs=list(item.specs), keep=keep)
                    row = {"item": item.label, "dps": list(item.dps or (item.dp,)),
                           "control": item.control, "direction": direction,
                           "member": case.label[-120:], "window_status": res["status"]}
                    try:
                        if item.control:
                            row["control_check"] = te.ladder_control(
                                case, res, keep[0] if keep else None, family=family,
                                counter=counter)
                            v = row["control_check"]["verdict"]
                            row["status"] = v
                            if res["status"] == "fixed_point" and \
                                    v != "local template covers it":
                                gates.append(f"L+1 control {case.label[-80:]}: {v}")
                        else:
                            entry = te.ladder_member(
                                case, res, family=family, counter=counter, rates=rates,
                                orientation="as_generated" if direction == "forward"
                                else "reversed")
                            row["template"] = entry
                            row["status"] = entry["status"]
                            row["refusal"] = (entry.get("refusal") or {}).get("status")
                            for c in entry.get("calibration", []) or []:
                                if not c["ok"]:
                                    gates.append(f"ladder calibration {case.label[-60:]}: {c}")
                    except te.EmissionError as error:
                        row.update(status="error", reason=str(error)[:300])
                        gates.append(f"ladder error {case.label[-60:]}: {error}")
                    row["seconds"] = round(time.time() - start, 1)
                    rows.append(row)
                    log.write(json.dumps(row, default=str) + "\n")
                    log.flush()
            print(f"[{time.time() - t0:6.0f}s] {k + 1}/{len(items)} {item.label[:60]}", flush=True)
    statuses = Counter((r.get("control"), r.get("status"), r.get("refusal")) for r in rows)
    emitted = sum(1 for r in rows if r.get("status") == "emitted")
    if rows and not emitted and not all(r.get("control") for r in rows):
        gates.append("every ladder member was refused (a vacuous run)")
    summary = {"family": name, "L": bound.L, "items": len(plan.items), "processed": len(items),
               "sample": args.sample, "limit": args.limit, "classes": len(classes),
               "statuses": [[list(map(str, k)), v] for k, v in statuses.items()],
               "emitted": emitted, "wall_s": round(time.time() - t0),
               "seconds": round(sum(r.get("seconds", 0) for r in rows), 1),
               "projection_full_s": round(sum(r.get("seconds", 0) for r in rows) *
                                          len(plan.items) / max(1, len(items))),
               "gates": gates}
    (out / f"summary.ladder_{name}.json").write_text(json.dumps(summary, indent=1, default=str))
    print("SUMMARY", json.dumps(summary, default=str)[:4000], flush=True)
    return 1 if gates else 0


if __name__ == "__main__":
    sys.exit(main())
