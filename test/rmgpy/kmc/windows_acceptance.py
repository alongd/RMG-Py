"""I-137 acceptance for rmgpy/kmc/windows.py (not collected).

Parts (each its own run, both streams persisted by the caller):

* ``required``: the required supported cases (H abstraction by a benzylic
  radical from a closed chain; chain-end recombination of two chains,
  identical and unequal; disproportionation; beta-scission of a mid-chain
  radical; an intramolecular H shift; a synthetic ``unchanged_species``
  chain-to-chain transfer) and the alkane provider on C16 homolysis; full
  context enumeration.
* ``dp8``: the closed PS DP8 reference compile (the I-127 probe artifact,
  ``--artifact``), certified with ``certify_records(required=True)``; every
  record by status; every certified record mapped onto the PS provider's DP8
  witness (template species, template, RMG isomorphism of the template side),
  bound by ``recertify`` (rate owner supplied), and every member of its
  channel run through items 1-5 and 7; full context enumeration.
* ``refs``: the DP8 head-radical and DP4 + 1-phenylethyl witnesses (the I-135
  reference families and directions, single molecules and pairs); contexts
  near the window (radical ends everywhere).
* ``ladder:<family>``: the family's ladder (``windows.ladder_plan``): the
  unimolecular directions on every single-chain item, the bimolecular ones on
  every pair; contexts near the window.

Every member row goes to ``rows.jsonl``; ``summary.json`` holds the tables and
the gate failures (non-empty: exit 1).  A span-bound violation fails the run.

Usage: python windows_acceptance.py OUTDIR PART [--workers N] [--artifact PATH]
"""

import argparse
import hashlib
import json
import os
import sys
import time
import traceback
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

_S = {}
REFERENCE_SETS = {"H_Abstraction": ("forward",), "Disproportionation": ("forward",),
                  "R_Recombination": ("forward", "reverse"),
                  "R_Addition_MultipleBond": ("forward", "reverse"),
                  "intra_H_migration": ("forward",)}


def _init():
    """Load the database once (the parent, before the pool forks) and reset
    the per-worker state (the pool initializer: bounded caches, recycled
    workers keep memory bounded)."""
    import channelsTest as ct
    from rmgpy.kmc import certificates as cm
    from rmgpy.kmc import windows as win
    if "fam" not in _S:
        db = ct.database()
        _S["fam"] = db.kinetics.families
        _S["bounds"] = {name: win.span_bound(win.family_constants(f))
                        for name, f in _S["fam"].items()}
    _S["certifier"] = cm.Certifier(_S["fam"], max_channel_sets=8, max_entries=1024)
    _S["cases"] = {}


def _rss_mb():
    try:
        with open("/proc/self/status") as handle:
            for line in handle:
                if line.startswith("VmHWM:"):
                    return int(line.split()[1]) // 1024
    except OSError:
        pass
    return -1


# --------------------------------------------------------------------------
# witness specs
# --------------------------------------------------------------------------
def spec_dict(spec):
    from dataclasses import asdict
    if isinstance(spec, tuple) and spec[0] == "alkane":
        return {"alkane": spec[1]}
    return asdict(spec)


def spec_from(d):
    from rmgpy.kmc import boundaries as b
    if "alkane" in d:
        return ("alkane", d["alkane"])
    d = dict(d)
    for key in ("radical_units", "secondary_radical_units", "pendant_loss_units"):
        d[key] = tuple(d.get(key) or ())
    return b.PsChainSpec(**d)


def builds_of(specs):
    from rmgpy.kmc import boundaries as b
    out = []
    for k, spec in enumerate(specs):
        if isinstance(spec, tuple) and spec[0] == "alkane":
            out.append(b.PROVIDERS["linear_alkane"].build(spec[1]))
        else:
            out.append(b.PROVIDERS["polystyrene"].build(f"p{k}", spec))
    return out


def cases_for(specs, family_name, direction):
    from rmgpy.kmc import windows as win
    key = (json.dumps([spec_dict(s) for s in specs], sort_keys=True), family_name, direction)
    if key not in _S["cases"]:
        if len(_S["cases"]) > 2:
            _S["cases"].clear()
        _S["cases"][key] = win.member_cases(_S["certifier"], _S["fam"][family_name],
                                            builds_of(specs), direction,
                                            label="+".join(_label(s) for s in specs))
    return _S["cases"][key]


def _label(spec):
    from rmgpy.kmc import windows as win
    if isinstance(spec, tuple):
        return f"C{spec[1]}"
    return win._ends_label(spec)


def _coords(case):
    return {a: (p.index, u, k) for p in case.witness.participants
            for a, u, k in p.metadata.coordinates}


def _ring(case):
    out = set()
    for m in case.reactants:
        out |= {a.id for a in m.atoms if m.is_atom_in_cycle(a)}
    return out


def _labelled(case):
    return [a for _s, _l, a in case.moves[0].labels]


# --------------------------------------------------------------------------
# required cases
# --------------------------------------------------------------------------
def required_sets():
    from rmgpy.kmc import boundaries as b
    S = b.PsChainSpec
    return [
        ("benzylic H abstraction from a closed chain", "H_Abstraction", "forward",
         [S(1, tail_end="benzylic_radical"), S(4)], "habs_interior"),
        ("chain-end recombination (identical chains)", "R_Recombination", "forward",
         [S(3, tail_end="benzylic_radical"), S(3, tail_end="benzylic_radical")], "rec_ends"),
        ("chain-end recombination (unequal chains)", "R_Recombination", "forward",
         [S(3, tail_end="benzylic_radical"), S(3, head_end="primary_radical")], "rec_ends"),
        ("disproportionation", "Disproportionation", "forward",
         [S(3, tail_end="benzylic_radical"), S(3, head_end="primary_radical")], "backbone"),
        ("beta-scission of a mid-chain radical", "R_Addition_MultipleBond", "reverse",
         [S(5, radical_units=(2,))], "cc_scission"),
        ("intramolecular H shift", "intra_H_migration", "forward",
         [S(4, head_end="primary_radical")], "backbone"),
        ("synthetic unchanged_species chain-to-chain transfer", "H_Abstraction", "forward",
         [S(3, radical_units=(1,)), S(3)], "unchanged"),
        ("C16 alkane homolysis", "R_Recombination", "reverse", [("alkane", 16)], "all"),
    ]


def select(name, case):
    from rmgpy.kmc import windows as win
    ring = _ring(case)
    labelled = _labelled(case)
    coords = _coords(case)
    if name == "all":
        return True
    if any(a in ring for a in labelled):
        return False
    if win.species_domain_refusal(case.witness):
        return False
    if name == "backbone":
        return True
    if name == "habs_interior":
        donor = case.moves[0].atom("*1")
        radical = case.moves[0].atom("*3")
        return coords.get(donor, (None,))[0] == 1 and coords[donor][1] in (1, 2) and \
            coords.get(radical) == (0, 0, 1)
    if name == "rec_ends":
        return all(coords.get(a) in ((0, 2, 1), (1, 2, 1), (1, 0, 0)) for a in labelled)
    if name == "cc_scission":
        rows = case.witness.reactant_rows()
        return all(rows[a][0] == "C" for a in labelled)
    if name == "unchanged":
        return case.kind == "synthetic:unchanged_species"
    raise ValueError(name)


# --------------------------------------------------------------------------
# tasks
# --------------------------------------------------------------------------
def run_cases_task(task):
    """(part, label, specs, family, direction, selection, contexts, near, chunk, of)."""
    from rmgpy.kmc import windows as win
    part, label, specs, family_name, direction, selection, contexts, near, chunk, of = task
    specs = [spec_from(s) for s in specs]
    start = time.time()
    rows, failures = [], []
    try:
        cases = cases_for(specs, family_name, direction)
        picked = [c for c in cases if selection is None or select(selection, c)]
        for k, case in enumerate(picked):
            if k % of != chunk:
                continue
            try:
                res = win.run_member(_S["certifier"], _S["fam"][family_name], case,
                                     bound=_S["bounds"][family_name], contexts=contexts,
                                     near_only=near)
            except win.SpanBoundError as error:
                failures.append(f"SPAN BOUND: {case.label}: {error}")
                res = {"label": case.label, "status": "span_bound_violation",
                       "reason": str(error)}
            res.update(part=part, case=label, kind=case.kind,
                       ladder_dp=None)
            rows.append(res)
    except Exception:   # reported, never hidden
        failures.append(f"{label} {family_name} {direction}: {traceback.format_exc()[-1500:]}")
    return {"task": [part, label, family_name, direction, chunk, of], "rows": rows,
            "failures": failures, "seconds": round(time.time() - start, 1),
            "rss_mb": _rss_mb()}


def run_dp8_task(task):
    """(record index chunk, of, artifact path, store path)."""
    from rmgpy.kmc import boundaries as b
    from rmgpy.kmc import certificates as cm
    from rmgpy.kmc import channels as chm
    from rmgpy.kmc import windows as win
    chunk, of, records_path, store_path = task
    start = time.time()
    rows, failures = [], []
    try:
        records = json.loads(Path(records_path).read_text())
        store = json.loads(Path(store_path).read_text())
        by_id = {r["event_id"]: r for r in records}
        certifier, fam = _S["certifier"], _S["fam"]
        family = fam["R_Recombination"]
        proxy = cases_for([b.PsChainSpec(8)], "R_Recombination", "reverse")
        species = chm.SpeciesKeys()
        proxy_keys = []
        for case in proxy:
            mols = win.template_molecules(case.witness, "reverse")
            proxy_keys.append((tuple(sorted(species(m) for m in mols)), ";".join(case.pathway)))
        done = {}
        for record in records:
            group = min(record["event_id"], record.get("reverse_of") or record["event_id"])
            if int(hashlib.sha256(group.encode()).hexdigest(), 16) % of != chunk:
                continue
            rc = store[record["rate_certificate"]["sha256"]]
            row = {"part": "dp8", "event_id": record["event_id"][:20],
                   "orientation": record["orientation"], "record_status": record["status"],
                   "rate_source": record["rate_source"]["kind"],
                   "certified": rc.get("certified"), "template": record["template"][:90]}
            if record["status"] not in cm.FIRING_STATUSES:
                # published refused by status (certification failed: e.g. a
                # domain refusal, an aryl-radical channel)
                row.update(status="refused:record_status", reason=record.get("status_reason",
                                                                             "")[:240])
                rows.append(row)
                continue
            if not cm.requires_certification(record, by_id):
                row.update(status="outside the required population",
                           reason=(rc.get("gaps") or [""])[0][:200])
                rows.append(row)
                continue
            if not rc.get("certified"):
                row.update(status="uncertified firing record", reason=(rc.get("gaps") or
                                                                      [""])[0][:200])
                failures.append(f"dp8 record {record['event_id'][:20]} fires uncertified")
                rows.append(row)
                continue
            _side, mols, _others = certifier.template_side(record)
            key = (tuple(sorted(species(m) for m in mols)), record["template"])
            members = [c for c, pk in zip(proxy, proxy_keys) if pk == key]
            row["proxy_members"] = len(members)
            if not members:
                row.update(status="unmapped", reason="no DP8 proxy member has the record's "
                                                     "template species and template")
                failures.append(f"dp8 record {record['event_id'][:20]} unmapped")
                rows.append(row)
                continue
            owner = cm._rate_owner(record, by_id)
            bind = bind_record(certifier, record, rc, members[0], owner)
            row["recertify"] = bind
            if not bind.get("same_selection"):
                failures.append(f"dp8 record {record['event_id'][:20]}: recertify on the proxy "
                                f"witness: {bind}")
            results = []
            for case in members:
                ck = (case.label,)
                if ck not in done:
                    done[ck] = win.run_member(certifier, family, case,
                                              bound=_S["bounds"]["R_Recombination"],
                                              contexts=True, near_only=False)
                results.append(done[ck])
            row["members"] = results
            row["status"] = ",".join(sorted({r["status"] for r in results}))
            rows.append(row)
    except win.SpanBoundError as error:
        failures.append(f"SPAN BOUND: {error}")
    except Exception:
        failures.append(f"dp8 chunk {chunk}: {traceback.format_exc()[-1500:]}")
    return {"task": ["dp8", chunk, of], "rows": rows, "failures": failures,
            "seconds": round(time.time() - start, 1), "rss_mb": _rss_mb()}


def bind_record(certifier, record, rc, case, owner):
    """``recertify`` the record on the proxy member's template side through an
    RMG isomorphism of the template sides."""
    from rmgpy.kmc import certificates as cm
    from rmgpy.kmc import windows as win
    _side, mols, _others = certifier.template_side(record)
    proxy = win.template_molecules(case.witness, "reverse")
    roots_r = [cm.aromatic_root(m) for m in mols]
    roots_p = [cm.aromatic_root(m) for m in proxy]
    amap, used = {}, set()
    for mr in roots_r:
        for k, mp in enumerate(roots_p):
            if k in used or len(mp.atoms) != len(mr.atoms):
                continue
            maps = mr.find_isomorphism(mp, strict=False)
            if maps:
                used.add(k)
                amap.update({a.id: c.id for a, c in maps[0].items()})
                break
        else:
            return {"same_selection": False, "reason": "no isomorphic proxy fragment"}
    res = certifier.recertify(record, proxy, amap, certificate=rc, rate_owner=owner)
    return {"same_selection": res["same_selection"], "same_support": res["same_support"],
            "fixed_point": res["fixed_point"], "refusals": res["refusals"][:2],
            "flips": len(res["flips"])}


# --------------------------------------------------------------------------
# task lists
# --------------------------------------------------------------------------
def required_tasks(workers):
    tasks = []
    for label, family, direction, specs, selection in required_sets():
        of = 2 if family in ("H_Abstraction", "R_Recombination") else 1
        for chunk in range(of):
            tasks.append(("required", label, [spec_dict(s) for s in specs], family, direction,
                          selection, True, False, chunk, of))
    return tasks


def ref_tasks():
    from rmgpy.kmc import boundaries as b
    S = b.PsChainSpec
    refs = {"DP8 head radical": [S(8, head_end="primary_radical")],
            "DP4 + 1-phenylethyl": [S(4), S(1, tail_end="benzylic_radical")]}
    tasks = []
    for label, specs in refs.items():
        groups = [[s] for s in specs] + ([specs] if len(specs) == 2 else [])
        for family, directions in sorted(REFERENCE_SETS.items()):
            for direction in directions:
                for group in groups:
                    of = 4 if len(group) == 2 else 2
                    for chunk in range(of):
                        tasks.append(("refs", label, [spec_dict(s) for s in group], family,
                                      direction, None, True, True, chunk, of))
    return tasks


def ladder_tasks(family_name):
    from rmgpy.kmc import windows as win
    fam = _S["fam"]
    constants = win.family_constants(fam[family_name])
    bound = win.span_bound(constants)
    plan = win.ladder_plan(constants, bound)
    molecularity = dict(constants.directions)
    tasks = []
    for item in plan.items:
        for direction, n in sorted(molecularity.items()):
            if n != len(item.specs):
                continue
            tasks.append(("ladder", f"{item.label}|{'/'.join(item.kinds)}|dp{item.dp}"
                          f"{'|control' if item.control else ''}",
                          [spec_dict(s) for s in item.specs], family_name, direction, None,
                          True, True, 0, 1))
    return tasks, plan, bound


# --------------------------------------------------------------------------
# summaries and gates
# --------------------------------------------------------------------------
def member_rows(rows):
    for row in rows:
        if "members" in row:
            for m in row["members"]:
                yield dict(m, part=row["part"], record=row["event_id"])
        elif row.get("part") != "dp8" or "label" in row:
            yield row


def summarize(part, rows, failures, extra=None):
    gates = list(failures)
    members = list(member_rows(rows))
    statuses = Counter(r.get("status") for r in members)
    if members and not statuses.get("fixed_point"):
        gates.append("every member was refused (a vacuous run)")
    for r in members:
        if r.get("status") == "fixed_point":
            if r.get("mapping_problems"):
                gates.append(f"{r['label']}: mapping edits outside the window "
                             f"{r['mapping_problems'][:2]}")
            if r.get("boundary_bond_problems"):
                gates.append(f"{r['label']}: boundary bonds {r['boundary_bond_problems'][:2]}")
            if not r.get("span_bound_held"):
                gates.append(f"{r['label']}: span bound not checked")
        elif not r.get("status") or r["status"] == "span_bound_violation":
            gates.append(f"{r.get('label')}: unaccounted status {r.get('status')}")
    contexts = Counter()
    for r in members:
        for c in r.get("contexts", []):
            contexts["merged" if c["merged_into"] is not None else c["status"]] += 1
    table = {"part": part, "members": len(members), "statuses": dict(statuses),
             "contexts": dict(contexts), "gates": gates}
    if extra:
        table.update(extra)
    return table


def required_gates(rows):
    out, table = [], {}
    by_case = defaultdict(list)
    for r in rows:
        by_case[r["case"]].append(r)
    for label, _f, _d, _s, _sel in required_sets():
        got = by_case.get(label, [])
        reached = [r for r in got if r["status"] == "fixed_point"]
        table[label] = {"members": len(got), "fixed_point": len(reached),
                        "statuses": dict(Counter(r["status"] for r in got)),
                        "mappings": sum(r.get("mappings", 0) for r in got),
                        "slot_orders": sorted({tuple(map(tuple, r.get("slot_orders", [])))
                                               for r in got})}
        if not reached:
            out.append(f"required case '{label}': no member reached a fixed point "
                       f"({dict(Counter(r['status'] for r in got))})")
        if label.startswith("C16"):
            continue
        for r in got:
            # a reacting atom on an artificial port is a cut-point event (its
            # physical-end analogue is a context, its long-chain analogue an
            # interior member): accounted, not required
            if r["status"] not in ("fixed_point", "refused:reacting_port"):
                out.append(f"required case '{label}': member {r['label'][-60:]} "
                           f"{r['status']}: {r.get('reason', '')[:160]}")
    return out, table


def ladder_gates(rows, plan, bound):
    out = []
    per_dp = defaultdict(Counter)
    for r in rows:
        dp = int(r["case"].split("|dp")[1].split("|")[0])
        per_dp[dp]["members"] += 1
        per_dp[dp][r.get("status")] += 1
        if r.get("status") == "fixed_point":
            both = any(v.get("both_ends") for v in r.get("span", {}).values())
            per_dp[dp]["both_ends"] += int(both)
            if both and dp > bound.L:
                out.append(f"control DP{dp} > L={bound.L}: a window holds both ends: "
                           f"{r['label'][-80:]}")
    for dp in (bound.L, bound.L + 1):
        # the L+1 control (and L itself) is vacuous without a fixed point
        if not per_dp[dp].get("fixed_point"):
            out.append(f"DP{dp} (L={bound.L}): no member reached a fixed point")
    return out, {str(k): dict(v) for k, v in sorted(per_dp.items())}


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def run_pool(tasks, workers, fn, log):
    results = []
    started = time.time()
    if workers > 1:
        import multiprocessing as mp
        context = mp.get_context("fork")
        _init()
        with context.Pool(workers, initializer=_init, maxtasksperchild=3) as pool:
            for res in pool.imap_unordered(fn, tasks):
                results.append(res)
                _log(res, log, started)
    else:
        _init()
        for task in tasks:
            res = fn(task)
            results.append(res)
            _log(res, log, started)
    return results


def _log(res, log, started):
    for row in res["rows"]:
        log.write(json.dumps(row, default=str) + "\n")
    log.flush()
    print(f"[{time.time() - started:6.0f}s] {res['task']} rows={len(res['rows'])} "
          f"fail={len(res['failures'])} t={res['seconds']}s rss={res['rss_mb']}MB", flush=True)
    for f in res["failures"]:
        print("FAILURE", f[:600], flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("outdir")
    parser.add_argument("part")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--artifact", default=None)
    args = parser.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(out / "mpl"))
    t0 = time.time()
    extra = {}
    with open(out / f"rows.{args.part.replace(':', '_')}.jsonl", "w") as log:
        if args.part == "required":
            results = run_pool(required_tasks(args.workers), args.workers, run_cases_task, log)
        elif args.part == "refs":
            results = run_pool(ref_tasks(), args.workers, run_cases_task, log)
        elif args.part.startswith("ladder:"):
            _init()
            tasks, plan, bound = ladder_tasks(args.part.split(":", 1)[1])
            extra = {"L": bound.L, "L_plus_1": bound.L + 1, "span_bound": bound.__dict__,
                     "ladder_items": len(plan.items),
                     "ladder_kinds": dict(Counter(k for i in plan.items for k in i.kinds)),
                     "ladder_excluded": len(plan.excluded),
                     "excluded_reasons": dict(Counter(r.split(":")[0]
                                                      for _k, r in plan.excluded)),
                     "tasks": len(tasks)}
            print("LADDER", json.dumps(extra, default=str), flush=True)
            results = run_pool(tasks, args.workers, run_cases_task, log)
        elif args.part == "dp8":
            results = dp8(out, args, log)
            extra = _S.get("dp8_extra", {})
        else:
            raise SystemExit(f"unknown part {args.part}")
    rows = [r for res in results for r in res["rows"]]
    failures = [f for res in results for f in res["failures"]]
    if args.part == "required":
        gates, table = required_gates(rows)
        failures += gates
        extra["required"] = table
    if args.part.startswith("ladder:"):
        gates, per_dp = ladder_gates(rows, None, _S["bounds"][args.part.split(":", 1)[1]])
        failures += gates
        extra["per_dp"] = per_dp
        Ls = [int(k) for k, v in per_dp.items() if v.get("both_ends")]
        extra["largest_dp_with_both_ends_in_a_window"] = max(Ls) if Ls else None
    summary = summarize(args.part, rows, failures, extra)
    summary["wall_s"] = round(time.time() - t0)
    (out / f"summary.{args.part.replace(':', '_')}.json").write_text(
        json.dumps(summary, indent=1, default=str))
    print("SUMMARY", json.dumps(summary, default=str)[:6000], flush=True)
    return 1 if summary["gates"] else 0


def dp8(out, args, log):
    """Certify the probe artifact's records (required) and run the records."""
    from rmgpy.kmc import certificates as cm
    _init()
    from rmgpy.kmc.event_record import EventRecord
    artifact = json.loads(Path(args.artifact).read_text())
    records = [EventRecord.from_dict(r) for r in artifact["records"]]
    t = time.time()
    certified, store, replacements = cm.certify_records(records, _S["fam"], _S["certifier"],
                                                         required=True)
    certified = [r.to_dict() if hasattr(r, "to_dict") else r for r in certified]
    statuses = Counter((r["status"], store[r["rate_certificate"]["sha256"]]["certified"])
                       for r in certified)
    _S["dp8_extra"] = {"records": len(records), "certify_s": round(time.time() - t),
                       "status_certified": [[list(k), v] for k, v in sorted(statuses.items())],
                       "replacements": len(replacements)}
    print("DP8 records", _S["dp8_extra"], flush=True)
    records_path, store_path = out / "dp8_records.json", out / "dp8_store.json"
    records_path.write_text(json.dumps(certified, default=str))
    store_path.write_text(json.dumps(store, default=str))
    of = max(1, args.workers * 3)
    tasks = [(k, of, str(records_path), str(store_path)) for k in range(of)]
    return run_pool(tasks, args.workers, run_dp8_task, log)


if __name__ == "__main__":
    sys.exit(main())
