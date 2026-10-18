"""Legacy/declared member parity (I-165 S2b): the window machinery run on a
legacy witness and on its declared twin must give the same members, paired
one to one, with equal outcomes.  Test-side, not collected.

``compare`` is pure: row dicts in, a ``Report`` out (no database, no RMG
objects).  The driver builds the legacy specs and their declared twins
(``declared_parity_support.legacy_to_chain`` / ``legacy_alkane_to_chain``),
runs ``member_cases`` and ``run_member(contexts=True)`` on both (the twins
under ``policy=LEGACY_V1``: parity tests the plumbing, not the policy) and
turns each member into a row whose fields are written in the declared
coordinate space (legacy positions through the position map ``pi``):

``pair``      the member's mapped labelled atoms (``windows.event_signatures``)
              and pathway, by coordinate;
``status``    ``run_member``'s status;
``window``    the closed window's marked graph (``window_labels``) by
              coordinate;
``view``      the member's ``certificate_view`` on the window, every atom
              (reacting ones included) by coordinate;
``span``      ``window_span`` extents per participant;
``contexts``  status per context, keyed by port and end class
              (``{methyl, benzyl} -> declared``,
              ``{primary_radical, benzylic_radical} -> radical``);
``deferred``  legacy contexts using ``vinylidene``/``styryl`` (no declared twin
              yet): listed, never compared.

Failures (exit 1): an unpaired member; a pair disagreeing on status, window,
view, span or the contexts both sides have; a member guarded
(``refused:whole_species_search``) on one side only; a reference/direction
with zero members on either side or a count off ``expected``; a guarded set
off ``expected`` (on the members run).  Listed: members guarded on both
sides, deferred contexts, numbered corrections C-1 (the depository growth
unit follows the provider; before frozen in ``data/c1_before.json``) and C-2
(alkane 5-tuples refuse, from S2a).

CLI (one family per process)::

    python declared_windows_parity.py --reference c16 --family R_Recombination \\
        --direction reverse [--limit N] [--offset K] [--sides both|legacy|declared] \\
        [--rows-out rows.json] [--out report.json]
    python declared_windows_parity.py compare legacy.json... --declared declared.json... \\
        [--out report.json]
    python declared_windows_parity.py expected rows.json... [--write]
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
DATA = HERE / "data"
EXPECTED = DATA / "declared_parity_expected.json"
C1_BEFORE = DATA / "c1_before.json"

GUARD = "refused:whole_species_search"
GUARD_RE = re.compile(r"refused:whole_species_search: (?P<node>\S+) label (?P<label>\S+) "
                      r"candidates (?P<candidates>.*)$")
END_CLASS = {"methyl": "declared", "benzyl": "declared", "declared": "declared",
             "primary_radical": "radical", "benzylic_radical": "radical",
             "radical": "radical"}
DEFERRED_ENDS = ("vinylidene", "styryl")
COMPARED = ("status", "window", "view", "span")


# --------------------------------------------------------------------------
# the pure comparator
# --------------------------------------------------------------------------
@dataclass
class Report:
    failures: list = field(default_factory=list)
    listed: dict = field(default_factory=lambda: {"guarded": [], "deferred": []})
    compared: int = 0
    notes: list = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.failures

    @property
    def exit_code(self) -> int:
        return 0 if self.ok else 1

    def as_dict(self) -> dict:
        return {"ok": self.ok, "failures": list(self.failures), "listed": self.listed,
                "compared": self.compared, "notes": list(self.notes)}


def is_guarded(row) -> bool:
    return str(row.get("status") or "").startswith(GUARD)


def guard_entry(row) -> dict:
    m = GUARD_RE.match(str(row.get("status") or ""))
    parts = m.groupdict() if m else {"node": None, "label": None, "candidates": None}
    return {"key": row["key"], "pair": row["pair"], "family": row.get("family"), **parts}


def _by_pair(rows, side, key, failures):
    out = {}
    for row in rows:
        if row["pair"] in out:
            failures.append(f"{key}: {side} has two members with pair {row['pair']}")
        out[row["pair"]] = row
    return out


def compare(legacy_rows, declared_rows, expected) -> Report:
    """Compare legacy and declared member rows against ``expected``
    (key -> {"members", "guarded", "deferred"})."""
    report = Report()
    groups = defaultdict(lambda: ([], []))
    for row in legacy_rows:
        groups[row["key"]][0].append(row)
    for row in declared_rows:
        groups[row["key"]][1].append(row)
    keys = sorted(set(groups) | set(expected))
    for key in keys:
        legacy, declared = groups[key] if key in groups else ([], [])
        want = expected.get(key)
        if want is None:
            report.failures.append(f"{key}: not in expected")
            continue
        for side, rows in (("legacy", legacy), ("declared", declared)):
            if not rows:
                report.failures.append(f"{key}: zero members on the {side} side")
            elif len(rows) != want["members"]:
                report.failures.append(f"{key}: {side} member count {len(rows)} != expected "
                                       f"{want['members']}")
        if not legacy or not declared:
            continue
        lmap = _by_pair(legacy, "legacy", key, report.failures)
        dmap = _by_pair(declared, "declared", key, report.failures)
        for pair in sorted(set(lmap) ^ set(dmap)):
            side = "legacy" if pair in lmap else "declared"
            report.failures.append(f"{key}: {side} member {pair} without a partner")
        run_pairs, guarded = set(), set()
        for pair in sorted(set(lmap) & set(dmap)):
            lrow, drow = lmap[pair], dmap[pair]
            for ctx in lrow.get("deferred") or ():
                report.listed["deferred"].append({"key": key, "pair": pair, "context": ctx})
            if not lrow.get("run", True) and not drow.get("run", True):
                continue
            if lrow.get("run", True) != drow.get("run", True):
                report.failures.append(f"{key}: pair {pair} run on one side only")
                continue
            run_pairs.add(pair)
            lg, dg = is_guarded(lrow), is_guarded(drow)
            if lg != dg:
                report.failures.append(f"{key}: pair {pair} guarded on one side only "
                                       f"({'legacy' if lg else 'declared'})")
                continue
            if lg:
                guarded.add(pair)
                report.listed["guarded"].append(guard_entry(lrow))
                continue
            report.compared += 1
            for name in COMPARED:
                if lrow.get(name) != drow.get(name):
                    report.failures.append(f"{key}: pair {pair} disagrees on {name}: "
                                           f"{str(lrow.get(name))[:120]} != "
                                           f"{str(drow.get(name))[:120]}")
            lctx, dctx = lrow.get("contexts") or {}, drow.get("contexts") or {}
            for ctx in sorted(set(lctx) & set(dctx)):
                if lctx[ctx] != dctx[ctx]:
                    report.failures.append(f"{key}: pair {pair} disagrees on contexts {ctx}: "
                                           f"{lctx[ctx]} != {dctx[ctx]}")
        want_guarded = {g for g in want.get("guarded", ()) if g in run_pairs}
        if guarded != want_guarded:
            report.failures.append(f"{key}: guarded set {sorted(guarded)} != expected "
                                   f"{sorted(want_guarded)} (on the members run)")
        deferred = sum(len(r.get("deferred") or ()) for r in legacy if r.get("run", True))
        if all(r.get("run", True) for r in legacy) and deferred != want.get("deferred", 0):
            report.notes.append(f"{key}: deferred contexts {deferred}, expected "
                                f"{want.get('deferred', 0)} (listed, not a failure)")
    return report


# --------------------------------------------------------------------------
# references
# --------------------------------------------------------------------------
def references() -> dict:
    from rmgpy.kmc import boundaries as b
    S = b.PsChainSpec
    return {
        "dp8_closed": ([S(8)], [("R_Recombination", "reverse")]),
        "dp8_headrad": ([S(8, head_end="primary_radical")],
                        [("intra_H_migration", "forward"),
                         ("R_Addition_MultipleBond", "reverse")]),
        "dp4_closed+benzyl": ([S(4), S(1, tail_end="benzylic_radical")],
                              [("H_Abstraction", "forward")]),
        "dp4_headrad+benzyl": ([S(4, head_end="primary_radical"),
                                S(1, tail_end="benzylic_radical")],
                               [("Disproportionation", "forward")]),
        "c16": ([("alkane", 16, None, None)], [("R_Recombination", "reverse")]),
    }


def twins(specs):
    """(declared twin specs, per participant index the legacy -> declared
    coordinate map)."""
    import declared_parity_support as sup
    from rmgpy.kmc import boundaries as b
    from rmgpy.kmc import declared_chain as dc
    from rmgpy.kmc import windows as win
    out, maps = [], []
    for spec in specs:
        if isinstance(spec, b.PsChainSpec):
            ctx = sup.ps_context()
            out.append(win.DeclaredSpec(ctx.provider.name, sup.legacy_to_chain(spec, ctx)))
            maps.append(ctx.coord_map)
        elif isinstance(spec, tuple) and spec and spec[0] == "alkane":
            pe = dc.DeclaredChainProvider(sup.pe_declaration())
            out.append(win.DeclaredSpec(pe.name, sup.legacy_alkane_to_chain(
                spec, pe.geometry.keys[0])))
            maps.append(lambda c: c)
        else:
            raise ValueError(f"no declared twin for {spec!r}")
    return out, maps


# --------------------------------------------------------------------------
# rows
# --------------------------------------------------------------------------
def _canon(value) -> str:
    return json.dumps(value, sort_keys=True, default=repr, separators=(",", ":"))


def coordinate_names(witness, maps=None) -> dict:
    """atom id -> its coordinate (participant, unit, declared position); a
    hydrogen -> ("H", its heavy neighbour's coordinate); else ("?", id)."""
    coords = {}
    for p in witness.participants:
        cmap = maps[p.index] if maps else (lambda c: c)
        for a, u, k in p.metadata.coordinates:
            uu, kk = cmap((u, k))
            coords[a] = (p.index, uu, kk)
    rows = witness.reactant_rows()
    out = {}
    for a, row in rows.items():
        if a in coords:
            out[a] = coords[a]
        elif row[0] == "H":
            heavy = [n for n, _o in row[5] if n in coords]
            out[a] = ("H", coords[heavy[0]]) if heavy else ("?", a)
        else:
            out[a] = ("?", a)
    return out


def pair_id(case, maps=None) -> str:
    from rmgpy.kmc import windows as win
    sigs = win.event_signatures(case)
    cmaps = maps or {}

    def mapped(x):
        if x[0] == "H":
            return ("H", mapped(x[1]))
        if x[0] == "?" or not cmaps:
            return x
        p, u, k = x
        uu, kk = cmaps[p]((u, k))
        return (p, uu, kk)
    canon = sorted(sorted(repr((lab, mapped(x))) for lab, x in sig) for sig in sigs)
    return _canon([canon, list(case.pathway)])


def number_pairs(rows) -> list:
    """Members whose mapped labelled atoms coincide by coordinate (hydrogens
    of one carbon) get ``#k`` in member order: paired by occurrence."""
    seen = Counter()
    out = []
    for row in rows:
        base = row["pair"].split("#")[0]
        out.append({**row, "pair": f"{base}#{seen[(row['key'], base)]}"})
        seen[(row["key"], base)] += 1
    return out


def _translate_labels(labels, names):
    nodes, edges = labels
    return {"nodes": sorted((repr(names.get(a, a)), repr(v)) for a, v in nodes.items()),
            "edges": sorted((repr(sorted((repr(names.get(a, a)), repr(names.get(b, b))))),
                             repr(v)) for (a, b), v in edges.items())}


def context_key(assignment, end_class=END_CLASS) -> tuple:
    """(key, deferred?) of a context row's assignment ([[port, end], ...])."""
    parts, deferred = [], False
    for port, end in assignment:
        if end in DEFERRED_ENDS:
            deferred = True
        parts.append(f"({port[0]},{port[1]!r}):{end_class.get(end, end)}")
    return ";".join(sorted(parts)), deferred


def member_row(key, family, case, result, fp, maps, run=True) -> dict:
    from rmgpy.kmc import windows as win
    row = {"key": key, "pair": pair_id(case, maps), "run": run, "family": family,
           "label": case.label[-160:], "status": result.get("status") if run else "not_run",
           "window": None, "view": None, "span": None, "contexts": {}, "deferred": []}
    if not run:
        return row
    if fp is not None and result.get("status") == "fixed_point":
        names = coordinate_names(fp.witness, maps)
        atoms = set(fp.window.atoms)
        row["window"] = _canon(_translate_labels(win.window_labels(fp.witness, atoms), names))
        cert = fp.observation.certificate
        named = {a: repr(n) for a, n in names.items()}      # one sortable type
        row["view"] = _canon(win.certificate_view(cert, case.wanted, atoms, named,
                                                  transport_reacting=True))
        row["span"] = {str(i): [s["backbone_atoms"], s["units"],
                                list(s["physical_ends_inside"])]
                       for i, s in (result.get("span") or {}).items()}
    for ctx in result.get("contexts") or ():
        ckey, deferred = context_key(ctx["assignment"])
        if deferred:
            row["deferred"].append(ckey)
        else:
            row["contexts"][ckey] = ctx["status"]
    return row


def run_side(family_name, specs, direction, *, key, side, policy=None, maps=None,
             limit=None, offset=0, stats=None) -> list:
    """Rows of one side: every member (members outside [offset, offset+limit)
    are rows that were not run)."""
    import one_family
    from rmgpy.kmc import certificates as cm
    from rmgpy.kmc import windows as win
    db = one_family.family_db(family_name)
    fams = db.kinetics.families
    family = fams[family_name]
    certifier = cm.Certifier(fams, max_channel_sets=16)
    kwargs = {} if policy is None else {"policy": policy}
    stats = {} if stats is None else stats
    builds = win.builds_of(specs)
    cases = win.member_cases(certifier, family, builds, direction, label=f"{key}|{side}",
                             stats=stats, **kwargs)
    provider = win.provider_object(builds[0].metadata.provider)
    bound = win.span_bound(win.family_constants(family), provider)
    rows = []
    for k, case in enumerate(cases):
        run = k >= offset and (limit is None or k < offset + limit)
        if not run:
            rows.append({**member_row(key, family_name, case, {}, None, maps, run=False),
                         "index": k})
            continue
        keep = []
        t0 = time.perf_counter()
        result = win.run_member(certifier, family, case, bound=bound, contexts=True,
                                specs=specs, keep=keep, **kwargs)
        row = member_row(key, family_name, case, result, keep[0] if keep else None, maps)
        row["seconds"] = round(time.perf_counter() - t0, 1)
        row["index"] = k
        rows.append(row)
    return number_pairs(rows)


def c1_table(indices) -> dict:
    """C-1 rows (member, before frozen, after live) for the c16 legacy
    members run (``indices``: their positions in ``member_cases`` order)."""
    import capture_c1_before as cap
    import one_family
    from rmgpy.kmc import certificates as cm
    from rmgpy.kmc import windows as win
    frozen = json.loads(C1_BEFORE.read_text())
    db = one_family.family_db(frozen["family"])
    fams = db.kinetics.families
    family = fams[frozen["family"]]
    certifier = cm.Certifier(fams, max_channel_sets=16)
    cases = win.member_cases(certifier, family, win.builds_of(cap.C16), frozen["direction"],
                             label="c16")
    cases = [c for k, c in enumerate(cases) if k in indices]
    after = cap.c1_rows(certifier, family, cases)
    out = []
    for label in sorted(after):
        out.append({"member": label, "before": frozen["members"][label], "after": after[label]})

    def verdicts(side):
        return dict(Counter(ch["verdict"] for r in out for ch in r[side]["channels"]))

    def units(side):
        return sorted({repr(ch["repeat_unit"]) for r in out for ch in r[side]["channels"]})
    changed = [r["member"] for r in out
               if [c["verdict"] for c in r["before"]["channels"]] !=
               [c["verdict"] for c in r["after"]["channels"]]]
    if changed:
        summary = (f"C-1: {len(changed)} of {len(out)} members change their "
                   f"depository_exclusion verdict; before {verdicts('before')} (repeat "
                   f"{units('before')}), after {verdicts('after')} (repeat {units('after')})")
    else:
        summary = (f"C-1: no depository_exclusion verdict changes on {len(out)} members; "
                   f"before {verdicts('before')} (repeat {units('before')}), after "
                   f"{verdicts('after')} (repeat {units('after')})")
    return {"rows": out, "changed": changed, "summary": summary}


def c2_table() -> dict:
    """C-2 (S2a): an alkane 5-tuple (interior radical) spec is refused."""
    from rmgpy.kmc import windows as win
    spec = ("alkane", 16, None, None, 7)
    try:
        win.builds_of([spec])
        after = "built"
    except Exception as error:   # the refusal is the after side
        after = str(error)[:160]
    return {"spec": repr(spec), "before": "built an interior radical (59246f18b)",
            "after": after}


def run_reference(reference, family, direction, *, limit=None, offset=0,
                  sides=("legacy", "declared")) -> dict:
    """Both sides of one reference/direction, compared against ``expected``."""
    from rmgpy.kmc import windows as win
    specs, pairs = references()[reference]
    if (family, direction) not in pairs:
        raise ValueError(f"{reference} has no {family} {direction} run")
    key = f"{reference}|{family}|{direction}"
    twin_specs, maps = twins(specs)
    stats = {"legacy": {}, "declared": {}}
    rows = {"legacy": [], "declared": []}
    seconds = {}
    for side in sides:
        t0 = time.perf_counter()
        if side == "legacy":
            rows[side] = run_side(family, specs, direction, key=key, side=side, maps=maps,
                                  limit=limit, offset=offset, stats=stats[side])
        else:
            rows[side] = run_side(family, twin_specs, direction, key=key, side=side,
                                  policy=win.LEGACY_V1, limit=limit, offset=offset,
                                  stats=stats[side])
        seconds[side] = round(time.perf_counter() - t0, 1)
    report = None
    if set(sides) == {"legacy", "declared"}:
        expected = json.loads(EXPECTED.read_text())
        report = compare(rows["legacy"], rows["declared"],
                         {key: expected[key]} if key in expected else {})
    corrections = {"C-2": c2_table()}
    if reference == "c16":
        corrections["C-1"] = c1_table({r["index"] for r in rows["legacy"] if r["run"]})
    return {"key": key, "rows": rows, "members": {s: len(rows[s]) for s in rows},
            "stats": stats, "seconds": seconds,
            "report": None if report is None else report.as_dict(),
            "corrections": corrections}


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def _write(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, sort_keys=True, indent=1, default=repr) + "\n")


def expected_from(rows) -> dict:
    """The expected entry of each key from full legacy rows."""
    out = {}
    by = defaultdict(list)
    for r in rows:
        by[r["key"]].append(r)
    for key, rs in sorted(by.items()):
        if not all(r["run"] for r in rs):
            raise ValueError(f"{key}: expected needs a full legacy run "
                             f"({sum(r['run'] for r in rs)} of {len(rs)} run)")
        out[key] = {"members": len(rs), "guarded": sorted(r["pair"] for r in rs if is_guarded(r)),
                    "deferred": sum(len(r["deferred"]) for r in rs)}
    return out


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "compare":
        parser = argparse.ArgumentParser(prog="declared_windows_parity.py compare")
        parser.add_argument("legacy", nargs="+")
        parser.add_argument("--declared", nargs="+", required=True)
        parser.add_argument("--out")
        args = parser.parse_args(argv[1:])

        def load(paths, side):
            merged = {}
            for p in paths:
                for r in number_pairs(json.loads(Path(p).read_text())["rows"][side]):
                    if r["run"] or (r["key"], r["pair"]) not in merged:
                        merged[(r["key"], r["pair"])] = r
            return list(merged.values())
        report = compare(load(args.legacy, "legacy"), load(args.declared, "declared"),
                         json.loads(EXPECTED.read_text()))
        if args.out:
            _write(args.out, report.as_dict())
        print(json.dumps({"ok": report.ok, "failures": report.failures[:20],
                          "compared": report.compared,
                          "guarded": len(report.listed["guarded"]),
                          "deferred": len(report.listed["deferred"])}, indent=1))
        return report.exit_code
    if argv and argv[0] == "expected":
        parser = argparse.ArgumentParser(prog="declared_windows_parity.py expected")
        parser.add_argument("rows", nargs="+")
        parser.add_argument("--write", action="store_true")
        args = parser.parse_args(argv[1:])
        merged = {}
        for p in args.rows:
            for r in number_pairs(json.loads(Path(p).read_text())["rows"]["legacy"]):
                if r["run"] or (r["key"], r["pair"]) not in merged:
                    merged[(r["key"], r["pair"])] = r
        entries = expected_from(list(merged.values()))
        print(json.dumps(entries, indent=1))
        if args.write:
            current = json.loads(EXPECTED.read_text()) if EXPECTED.is_file() else {}
            current.update(entries)
            _write(EXPECTED, current)
        return 0
    parser = argparse.ArgumentParser(prog="declared_windows_parity.py")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--family", required=True)
    parser.add_argument("--direction", required=True)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--offset", type=int, default=0)
    parser.add_argument("--sides", default="both", choices=("both", "legacy", "declared"))
    parser.add_argument("--out", default=None)
    args = parser.parse_args(argv)
    sides = ("legacy", "declared") if args.sides == "both" else (args.sides,)
    result = run_reference(args.reference, args.family, args.direction, limit=args.limit,
                           offset=args.offset, sides=sides)
    if args.out:
        _write(args.out, result)
    report = result["report"]
    print(json.dumps({"key": result["key"], "members": result["members"],
                      "seconds": result["seconds"], "stats": result["stats"],
                      "report": None if report is None else
                      {k: report[k] for k in ("ok", "compared", "notes")} |
                      {"failures": report["failures"][:20],
                       "guarded": report["listed"]["guarded"],
                       "deferred": len(report["listed"]["deferred"])},
                      "C-1": (result["corrections"].get("C-1") or {}).get("summary"),
                      "C-2": result["corrections"]["C-2"]}, indent=1, default=repr))
    return 0 if report is None else (0 if report["ok"] else 1)


if __name__ == "__main__":
    sys.exit(main())
