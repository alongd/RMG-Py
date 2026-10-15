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
* ``roster:<family>``: the family's generated class roster alone
  (``windows.probe_pattern`` per declared pattern through the pool): the
  class table per direction and slot, class and generated item counts.
* ``ladder:<family>``: the family's ladder (``windows.ladder_plan`` with the
  kind roster of ``windows.pair_roster`` and the class roster, probed first
  through the pool): the class closings first (``windows.close_classes``,
  one task per direction, slot and spec, against ``windows.closing_partners``;
  ``--closures`` reuses saved ones, ``--close-chunk i/n`` runs one chunk of
  them and stops), the merge (``windows.merge_classes``), then the unimolecular
  directions on every single-chain item, each pair in the directions its
  kinds or classes are rostered for; coverage checked per direction
  (``windows.ladder_direction_gaps``), per class (``windows.ladder_class_gaps``,
  ``windows.closure_gaps``) and per planned leaf tuple
  (``windows.leaf_tuple_gaps``), the L + 1 control per participant.

Every reacting-port refusal carries one cover per event it stands for (the
witness specs are passed to ``run_member``); a member with an uncovered event
fails the run.  ``--sample K``
runs every K-th task only (a timing projection: ``projection`` in the
summary).

Every member row goes to ``rows.jsonl``; ``summary.json`` holds the tables and
the gate failures (non-empty: exit 1).  A span-bound violation fails the run.

Usage: python windows_acceptance.py OUTDIR PART [--workers N] [--artifact PATH]
       [--sample K] [--roster classes.<family>.pkl] [--closures closures.<family>*.pkl ...]
       [--close-chunk i/n]

Task times: ``seconds`` (wall) and ``cpu`` (the worker's process time; the
projection under a loaded machine).
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


def _cpu():
    return time.process_time()


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
        return {"alkane": spec[1], "ends": list(spec[2:4])}
    if isinstance(spec, tuple) and spec[0] == "small":
        return {"small": spec[1]}
    return asdict(spec)


def spec_from(d):
    from rmgpy.kmc import boundaries as b
    if "alkane" in d:
        return ("alkane", d["alkane"]) + tuple(d.get("ends") or ())
    if "small" in d:
        return ("small", d["small"])
    d = dict(d)
    for key in ("radical_units", "secondary_radical_units", "pendant_loss_units"):
        d[key] = tuple(d.get(key) or ())
    return b.PsChainSpec(**d)


def builds_of(specs):
    from rmgpy.kmc import windows as win
    return win.builds_of(specs)


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
    if isinstance(spec, tuple) and spec[0] == "alkane":
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
    """(part, label, specs, family, direction, selection, contexts, meta, chunk, of);
    ``meta``: ladder fields copied to every row (dps, pair_kinds, long_index)."""
    from rmgpy.kmc import windows as win
    part, label, specs, family_name, direction, selection, contexts, meta, chunk, of = task
    specs = [spec_from(s) for s in specs]
    start, cpu = time.time(), _cpu()
    rows, failures = [], []
    try:
        cases = cases_for(specs, family_name, direction)
        picked = [c for c in cases if selection is None or select(selection, c)]
        if meta and meta.get("sites"):
            picked = [c for c in picked if win.site_select(c, specs, meta["sites"])]
        for k, case in enumerate(picked):
            if k % of != chunk:
                continue
            try:
                res = win.run_member(_S["certifier"], _S["fam"][family_name], case,
                                     bound=_S["bounds"][family_name], contexts=contexts,
                                     specs=specs)
            except win.SpanBoundError as error:
                failures.append(f"SPAN BOUND: {case.label}: {error}")
                res = {"label": case.label, "status": "span_bound_violation",
                       "reason": str(error)}
            res.update(part=part, case=label, kind=case.kind,
                       **{k: v for k, v in (meta or {}).items() if k != "sites"})
            rows.append(res)
    except Exception:   # reported, never hidden
        failures.append(f"{label} {family_name} {direction}: {traceback.format_exc()[-1500:]}")
    return {"task": [part, label, family_name, direction, chunk, of], "rows": rows,
            "failures": failures, "seconds": round(time.time() - start, 1),
            "cpu": round(_cpu() - cpu, 1), "rss_mb": _rss_mb()}


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
                                              contexts=True, specs=[b.PsChainSpec(8)])
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
                          selection, True, None, chunk, of))
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
                                      direction, None, True, None, chunk, of))
    return tasks


def run_probe_task(task):
    """(family, pattern): the site classes the family's templates accept for
    one declared pattern (``windows.probe_pattern``)."""
    from rmgpy.kmc import windows as win
    family_name, pattern, partners = task
    start = time.time()
    failures, classes = [], ()
    try:
        family = _S["fam"][family_name]
        classes = win.probe_pattern(_S["certifier"], family, pattern,
                                    win.family_constants(family), _S["bounds"][family_name],
                                    partners=partners)
    except Exception:   # reported, never hidden
        failures.append(f"probe {family_name} {pattern.label}: {traceback.format_exc()[-1500:]}")
    return {"task": ["probe", family_name, pattern.label], "rows": [], "classes": classes,
            "failures": failures, "seconds": round(time.time() - start, 1),
            "rss_mb": _rss_mb()}


def class_roster(family_name, workers, log):
    """The family's generated class roster, one probe task per declared
    pattern through the pool."""
    from rmgpy.kmc import windows as win
    t = time.time()
    _init()
    family = _S["fam"][family_name]
    partners = win.probe_partners(_S["certifier"], family, win.family_constants(family))
    calibrate = round(time.time() - t, 1)
    print("PROBE PARTNERS", {d: {s: n for s, (n, _x) in v.items()} for d, v in partners.items()},
          calibrate, "s", flush=True)
    for pattern in win.participant_patterns():    # cached before the pool forks
        for dp in range(1, win.PATTERN_DP_LIMIT + 1):
            win.pattern_specs(pattern, dp)
    tasks = [(family_name, p, partners) for p in win.participant_patterns()]
    results = run_pool(tasks, workers, run_probe_task, log)
    classes = tuple(sorted((c for res in results for c in res["classes"]),
                           key=lambda c: (c.direction, c.slot, repr(c.site), c.pattern)))
    failures = [f for res in results for f in res["failures"]]
    return classes, failures, {"probe_wall_s": round(time.time() - t),
                               "probe_task_s": round(sum(r["seconds"] for r in results)),
                               "probe_calibration_s": calibrate,
                               "probe_partners": {d: {str(s): n for s, (n, _x) in v.items()}
                                                  for d, v in partners.items()}}


def run_close_task(task):
    """(family, classes of one direction, slot and spec, partner (name, spec)
    or None): ``windows.close_classes``."""
    from rmgpy.kmc import windows as win
    family_name, classes, partner = task
    start, cpu = time.time(), _cpu()
    rows, failures, closures = [], [], ()
    first = classes[0]
    try:
        closures, rows = win.close_classes(_S["certifier"], _S["fam"][family_name], list(classes),
                                           partner, bound=_S["bounds"][family_name])
        for r in rows:
            r.update(part="ladder", case=f"close|{first.direction}|slot{first.slot}|"
                                         f"{win._spec_name(first.spec)}",
                     kind=r.get("kind", "kept"), dps=[first.dp], pair_kinds=None,
                     long_index=None, ladder_kinds=["class:close"],
                     closing_partner=partner[0] if partner else None)
    except Exception:   # reported, never hidden
        failures.append(f"close {family_name} {win.class_label(first)}: "
                        f"{traceback.format_exc()[-1500:]}")
    return {"task": ["close", family_name, first.direction, first.slot,
                     win._spec_name(first.spec), partner[0] if partner else None],
            "rows": rows, "closures": closures, "failures": failures,
            "seconds": round(time.time() - start, 1), "cpu": round(_cpu() - cpu, 1),
            "rss_mb": _rss_mb()}


def class_closings(family_name, classes, workers, log, chunk=None):
    """Every class closed once (one task per direction, slot and spec,
    against the closing partner of the other slot)."""
    from rmgpy.kmc import windows as win
    _init()
    family = _S["fam"][family_name]
    constants = win.family_constants(family)
    t = time.time()
    partners = win.closing_partners(_S["certifier"], family, constants)
    partner_s = round(time.time() - t, 1)
    print("CLOSING PARTNERS", {d: {s: n for s, (n, _x) in v.items()} for d, v in
                               partners.items()}, partner_s, "s", flush=True)
    molecularity = dict(constants.directions)
    fallback = dict(win.CLOSING_CANDIDATES)
    groups = defaultdict(list)
    for c in classes:
        if molecularity.get(c.direction) == 2:     # unimolecular classes run once, alone
            groups[(c.direction, c.slot, repr(c.spec))].append(c)
    tasks, failures = [], []
    for key in sorted(groups):
        direction, slot, _spec = key
        partner = None
        if molecularity.get(direction) == 2:
            partner = partners.get(direction, {}).get(1 - slot)
            if partner is None:
                name = groups[key][0].partner
                partner = (name, fallback[name]()) if name in fallback else None
                failures.append(f"no closing partner for {direction} slot {1 - slot}; "
                                f"probe partner {name}")
        tasks.append((family_name, tuple(groups[key]), partner))
    if chunk:
        k, n = chunk
        tasks = tasks[k::n]
    results = run_pool(tasks, workers, run_close_task, log)
    closures = [c for res in results for c in res.get("closures", ())]
    for res in results:
        res.pop("closures", None)
    timing = {"closing_partners": {d: {str(s): n for s, (n, _x) in v.items()}
                                   for d, v in partners.items()},
              "closing_partner_s": partner_s, "closing_tasks": len(tasks),
              "closing_task_s": round(sum(r["seconds"] for r in results)),
              "closing_cpu_s": round(sum(r.get("cpu", 0) for r in results))}
    results.append({"rows": [], "failures": failures, "seconds": 0, "cpu": 0})
    return closures, results, timing


def ladder_tasks(family_name, classes=(), merge=None):
    from dataclasses import replace
    from rmgpy.kmc import windows as win
    fam = _S["fam"]
    family = fam[family_name]
    constants = win.family_constants(family)
    t = time.time()
    roster = win.pair_roster(_S["certifier"], family)
    constants = replace(constants, pair_kinds=roster, classes=tuple(classes), merge=merge)
    bound = win.span_bound(constants)
    plan = win.ladder_plan(constants, bound)
    _S["roster_s"] = round(time.time() - t, 1)
    molecularity = dict(constants.directions)
    tasks = []
    for item in plan.items:
        if len(item.specs) == 1:
            directions = [d for d, n in sorted(molecularity.items()) if n == 1]
        else:
            directions = list(item.directions)
        dps = list(item.dps) + ([None] if item.long_index is not None else [])
        meta = {"dps": dps, "pair_kinds": list(item.pair_kinds) or None,
                "long_index": item.long_index, "ladder_kinds": list(item.kinds)}
        if item.sites:
            meta["sites"] = [[k, s, [list(x) for x in site]] for k, s, site in item.sites]
        of = 2 if item.long_index is not None else 1
        for direction in directions:
            for chunk in range(of):
                tasks.append(("ladder", f"{item.label}|{'/'.join(item.kinds)}|dp{item.dp}"
                              f"{'|control' if item.control else ''}",
                              [spec_dict(s) for s in item.specs], family_name, direction, None,
                              True, meta, chunk, of))
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
    """Per direction coverage (``windows.ladder_direction_gaps``) and the
    L + 1 control per participant: no window holds both ends of a short
    participant above L."""
    from rmgpy.kmc import windows as win
    out = list(win.ladder_direction_gaps(plan, rows))
    out += win.ladder_class_gaps(plan, _S.get("classes", ()))
    out += win.closure_gaps(_S.get("classes", ()), _S.get("merge"),
                            {d for d, n in plan.pair_kinds} | set(_S.get("closed_directions", ())))
    out += win.leaf_tuple_gaps(plan, rows)
    per_dp = defaultdict(Counter)
    per_pair = defaultdict(Counter)
    for r in rows:
        dps = r.get("dps") or []
        top = max([d for d in dps if d is not None] or [0])
        per_dp[top]["members"] += 1
        per_dp[top][r.get("status")] += 1
        kinds = "+".join(r["pair_kinds"]) if r.get("pair_kinds") else "single"
        group = f"{r.get('direction')}|{kinds}|{'long' if r.get('long_index') else 'short'}"
        per_pair[group]["members"] += 1
        per_pair[group][r.get("status")] += 1
        if r.get("status") == "fixed_point":
            span = r.get("span", {})
            both = any(v.get("both_ends") for v in span.values())
            per_dp[top]["both_ends"] += int(both)
            for index, dp in enumerate(dps):
                if dp is not None and dp > bound.L and span.get(str(index), {}).get("both_ends"):
                    out.append(f"control DP{dp} > L={bound.L}: a window holds both ends of "
                               f"participant {index}: {r['label'][-80:]}")
    return out, {str(k): dict(v) for k, v in sorted(per_dp.items())}, \
        {k: dict(v) for k, v in sorted(per_pair.items())}


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
          f"fail={len(res['failures'])} t={res['seconds']}s cpu={res.get('cpu')}s "
          f"rss={res['rss_mb']}MB", flush=True)
    for f in res["failures"]:
        print("FAILURE", f[:600], flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("outdir")
    parser.add_argument("part")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--artifact", default=None)
    parser.add_argument("--sample", type=int, default=1)
    parser.add_argument("--roster", default=None,
                        help="classes.<family>.pkl of a roster:<family> run (else probed)")
    parser.add_argument("--closures", nargs="*", default=None,
                        help="closures.<family>*.pkl of earlier closing runs (else closed)")
    parser.add_argument("--task-chunk", default=None,
                        help="i/n: run chunk i of n of the (sampled) plan tasks")
    parser.add_argument("--close-chunk", default=None,
                        help="i/n: run closing chunk i of n, save it and stop")
    args = parser.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(out / "mpl"))
    t0 = time.time()
    extra = {}
    close_results = []
    with open(out / f"rows.{args.part.replace(':', '_')}.jsonl", "w") as log:
        if args.part == "required":
            results = run_pool(required_tasks(args.workers)[::args.sample], args.workers,
                               run_cases_task, log)
        elif args.part == "refs":
            results = run_pool(ref_tasks()[::args.sample], args.workers, run_cases_task, log)
        elif args.part.startswith("roster:"):
            name = args.part.split(":", 1)[1]
            classes, failures0, timing = class_roster(name, args.workers, log)
            import pickle
            (out / f"classes.{name}.pkl").write_bytes(pickle.dumps(classes))
            from rmgpy.kmc import windows as win
            extra = dict(timing, **win.class_summary(classes, win.family_constants(
                _S["fam"][name])), failures=failures0)
            results = [{"rows": [], "failures": failures0, "seconds": 0}]
        elif args.part.startswith("ladder:"):
            from rmgpy.kmc import windows as win
            _init()
            name = args.part.split(":", 1)[1]
            if args.roster:
                import pickle
                classes = pickle.loads(Path(args.roster).read_bytes())
                failures0, timing = [], {"roster_from": args.roster}
            else:
                classes, failures0, timing = class_roster(name, args.workers, log)
                import pickle
                (out / f"classes.{name}.pkl").write_bytes(pickle.dumps(classes))
            import pickle
            family = _S["fam"][name]
            constants = win.family_constants(family)
            chunk = tuple(int(x) for x in args.close_chunk.split("/")) if args.close_chunk \
                else None
            close_results = []
            _S["plan"] = None
            if args.closures:
                closures = [c for path in args.closures
                            for c in pickle.loads(Path(path).read_bytes())]
                timing["closures_from"] = args.closures
            else:
                closures, close_results, ctiming = class_closings(name, classes, args.workers,
                                                                  log, chunk)
                timing.update(ctiming)
                suffix = f".{chunk[0]}of{chunk[1]}" if chunk else ""
                (out / f"closures.{name}{suffix}.pkl").write_bytes(pickle.dumps(closures))
            if chunk:
                results = close_results
                extra.update(timing)
                _S["plan"] = None
            else:
                t = time.time()
                merged_root = tuple(d for d, n in constants.directions
                                    if n == 2 and d not in constants.per_reactant)
                # merged roots: every predicate of any descent reads within r bonds
                # of the labelled atoms (the tree's largest predicate radius)
                merge = win.merge_classes(classes, closures, merged_root=merged_root,
                                          reach=constants.radius)
                timing["merge_s"] = round(time.time() - t, 1)
                timing["merge"] = win.merge_summary(classes, merge)
                timing["per_reactant"] = list(constants.per_reactant)
                timing["merged_root"] = list(merged_root)
                timing["joined"] = [list(x) for x in constants.joined]
                _S["merge"] = merge
                _S["closed_directions"] = tuple(d for d, n in constants.directions if n == 2)
            if not chunk:
                tasks, plan, bound = ladder_tasks(name, classes, _S["merge"])
                _S["classes"] = classes
                pairs = Counter((d, "+".join(i.pair_kinds), "long" if i.long_index else "short")
                                for i in plan.items for d in i.directions)
                extra = {"L": bound.L, "L_plus_1": bound.L + 1, "span_bound": bound.__dict__,
                         "roster": [[d, list(k)] for d, k in plan.pair_kinds],
                         "roster_s": _S.get("roster_s"),
                         "single_directions": list(plan.single_directions),
                         "ladder_items": len(plan.items),
                         "ladder_kinds": dict(Counter(k for i in plan.items for k in i.kinds)),
                         "pair_items": {"|".join(k): v for k, v in sorted(pairs.items())},
                         "ladder_excluded": len(plan.excluded),
                         "excluded_reasons": dict(Counter(r.split(":")[0]
                                                          for _k, r in plan.excluded)),
                         "tasks": len(tasks), "sample": args.sample}
                extra.update(timing)
                extra.update(win.class_summary(classes, win.family_constants(_S["fam"][name])))
                _S["plan"] = plan
                print("LADDER", json.dumps(extra, default=str), flush=True)
                picked = tasks[::args.sample]
                if args.task_chunk:
                    k, n = (int(x) for x in args.task_chunk.split("/"))
                    picked = picked[k::n]
                    extra["task_chunk"] = args.task_chunk
                results = run_pool(picked, args.workers, run_cases_task, log)
                _S["plan_task_s"] = sum(r["seconds"] for r in results)
                _S["plan_cpu_s"] = sum(r.get("cpu", 0) for r in results)
                results += close_results
                results.append({"rows": [], "failures": failures0, "seconds": 0})
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
    if args.part.startswith("ladder:") and _S.get("plan") is not None:
        gates, per_dp, per_pair = ladder_gates(rows, _S["plan"],
                                               _S["bounds"][args.part.split(":", 1)[1]])
        if args.sample == 1:
            failures += gates
        else:
            extra["sampled_gates"] = gates
        extra["per_dp"] = per_dp
        extra["per_direction_kind"] = per_pair
        Ls = [dp for r in rows if r.get("status") == "fixed_point"
              for index, dp in enumerate(r.get("dps") or [])
              if dp is not None and r.get("span", {}).get(str(index), {}).get("both_ends")]
        extra["largest_dp_with_both_ends_in_a_window"] = max(Ls) if Ls else None
    from rmgpy.kmc import windows as win
    failures += win.reacting_port_gaps(list(member_rows(rows)))
    covers = [r for r in member_rows(rows) if r.get("status") == "refused:reacting_port"]
    extra["reacting_port_covers"] = [
        {"label": r["label"][-110:], "covered": win.cover_verdict(r.get("covers") or []),
         "events": [[c.get("event"), c.get("route"), c.get("verdict"), c.get("status"),
                     c.get("window_atoms"), (c.get("reason") or "")[:120]]
                    for c in r.get("covers") or []]} for r in covers]
    extra["reacting_port_events"] = dict(Counter(
        (c.get("event") or "").split("@")[0] + ":" + str(c.get("verdict"))
        for r in covers for c in r.get("covers") or []))
    summary = summarize(args.part, rows, failures, extra)
    summary["wall_s"] = round(time.time() - t0)
    task_s = sum(res["seconds"] for res in results)
    summary["task_seconds"] = round(task_s)
    summary["projection"] = {"sample": args.sample, "task_seconds_full": round(
        task_s * args.sample), "wall_s_full_at_workers": round(
        task_s * args.sample / max(1, args.workers)), "workers": args.workers}
    if args.part.startswith("ladder:"):
        # closings run whole (never sampled); the plan's tasks are sampled
        close_s = sum(r["seconds"] for r in close_results)
        close_cpu = sum(r.get("cpu", 0) for r in close_results)
        plan_s, plan_cpu = _S.get("plan_task_s", 0), _S.get("plan_cpu_s", 0)
        summary["projection"].update(
            closing_task_s=round(close_s), closing_cpu_s=round(close_cpu),
            plan_task_s_full=round(plan_s * args.sample),
            plan_cpu_s_full=round(plan_cpu * args.sample),
            task_seconds_full=round(close_s + plan_s * args.sample),
            cpu_seconds_full=round(close_cpu + plan_cpu * args.sample))
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
    return run_pool(tasks[::args.sample], args.workers, run_dp8_task, log)


if __name__ == "__main__":
    sys.exit(main())
