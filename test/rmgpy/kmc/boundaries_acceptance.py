"""I-136 acceptance sweep for rmgpy/kmc/boundaries.py (not collected).

For every PS proxy of ``compiler.ps_proxy_set(3)`` (the 3-unit set plus its
``@5`` copies), every family in the proxy's ``family_candidates`` and every
direction the family has, RMG generates the reaction witnesses with
``rmg_moves.generate_moves`` (no kinetics).  Generation is streamed one
reactant resonance-form combination at a time (``resonance="given"`` on each
form ``reactant_forms`` returns); that is exactly the loop of
``resonance="all"`` (checked on the small tasks with ``--self-check``) but
keeps memory per task bounded.

* Every kept move becomes a witness; metadata survival is asserted on every
  one (``witness_from_move``).
* Extension checks run on one representative per mapped reaction: moves that
  share (family, direction, labelled atom ids) differ only in the resonance
  form of the reactants; the representative is the move with the least
  ``sort_key``.  Every member's connectivity change (bonds made or broken,
  H included) must equal the representative's, else the class is reported
  as split.  Representatives whose witness graphs are identical (e.g. the
  two template slot orders of a symmetric recombination) are counted as
  exact duplicates and checked once.
* Vectors per representative: the full box {0,1,2}^k when k <= 2 artificial
  ports, otherwise +1 and +2 at each port alone plus all ports at +1 and at
  +2 jointly; the first representative of every task gets the full box.
  Every vector is verified (composite map, sequential equality, reversed
  port order); marked duplicates are searched among each witness's vectors.

Alkanes: C8 and C16 homolysis witnesses (R_Recombination, reverse), full box
{0,1,2}^2, plus marked duplicates across all homolysis witnesses.

Usage: python boundaries_acceptance.py OUTDIR [--workers N] [--self-check]
"""

import argparse
import json
import sys
import time
import traceback
from collections import Counter
from dataclasses import replace
from itertools import product
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

FULL_BOX_MAX_PORTS = 2
_STATE = {}


def _init():
    import rmg_moves_harness as harness
    from rmgpy.kmc.compiler import ps_proxy_set

    _STATE["families"] = harness.families()
    _STATE["proxies"] = ps_proxy_set(3)


def _rss_mb():
    try:
        with open("/proc/self/status") as handle:
            for line in handle:
                if line.startswith("VmHWM:"):
                    return int(line.split()[1]) // 1024
    except OSError:
        pass
    return -1


def _connectivity(witness):
    before, after = witness.reactant_bonds(), witness.product_bonds()
    return tuple(sorted(pair for pair in set(before) | set(after)
                        if (pair in before) != (pair in after)))


def vectors_for(witness, full_box):
    from rmgpy.kmc import boundaries as b

    k = len(b.vector_ports(witness))
    if full_box or k <= FULL_BOX_MAX_PORTS:
        return list(product(range(3), repeat=k)), True
    vectors = [tuple([0] * k)]
    for position in range(k):
        for depth in (1, 2):
            vector = [0] * k
            vector[position] = depth
            vectors.append(tuple(vector))
    vectors += [tuple([1] * k), tuple([2] * k)]
    return vectors, False


def check_witness(witness, tally, full_box=False):
    """Evaluate the vector set of one witness into ``tally``."""
    from rmgpy.kmc import boundaries as b

    vectors, boxed = vectors_for(witness, full_box)
    tally["full_box_witnesses"] += int(boxed)
    plan = b._vector_plan(witness)
    tally["refused_ports"] += len(plan[0])
    results = []
    for vector in vectors:
        try:
            result = b.evaluate_vector(witness, vector, plan=plan)
        except b.BoundaryError as error:
            tally["failures"] += 1
            tally["failure_reasons"].append(f"{witness.label} {vector}: {error}")
            continue
        results.append(result)
        tally["vectors"] += 1
        tally[f"status:{result.status}"] += 1
        if result.status == "refused":
            tally["refused_reasons"][_reason_kind(result.reason)] += 1
    duplicates = b.vector_duplicates(results)
    tally["duplicates"] += len(duplicates)
    for first, second in duplicates:
        tally["duplicate_examples"].append(f"{witness.label}: {first} ~ {second}")


def _reason_kind(reason):
    if "reacting atom" in reason:
        return "port is a reacting atom"
    if "no extension path" in reason:
        return "no extension path"
    return reason


def new_tally():
    return {"kept_moves": 0, "witnesses_metadata_ok": 0, "metadata_failures": [],
            "classes": 0, "class_splits": [], "exact_duplicate_witnesses": 0,
            "checked_witnesses": 0, "full_box_witnesses": 0, "vectors": 0,
            "status:ok": 0, "status:refused": 0, "status:ceiling_exceeded": 0,
            "refused_ports": 0, "refused_reasons": Counter(), "failures": 0,
            "failure_reasons": [], "duplicates": 0, "duplicate_examples": [],
            "statuses": Counter()}


def run_task(task):
    from rmgpy.kmc import boundaries as b
    from rmgpy.kmc import rmg_moves as rm

    proxy_index, family_name, direction, self_check = task
    proxy = _STATE["proxies"][proxy_index]
    family = _STATE["families"][family_name]
    start = time.time()
    tally = new_tally()
    try:
        molecules, participants = b.witness_reactants(proxy.boundary_builds)
        forms = [rm.reactant_forms(m) for m in molecules]
        classes, seen_graphs = {}, set()
        first = True
        streamed = []
        for combo in product(*(range(len(f)) for f in forms)):
            given = [forms[k][combo[k]] for k in range(len(forms))]
            moves = rm.generate_moves(family, given, direction, resonance="given")
            for move in moves:
                tally["statuses"][move.status] += 1
                if self_check:
                    streamed.append(replace(move, forms=combo).sort_key())
                if move.status != "kept":
                    continue
                move = replace(move, forms=combo)
                tally["kept_moves"] += 1
                try:
                    witness = b.witness_from_move(move, molecules, participants, forms=forms)
                except b.BoundaryError as error:
                    tally["metadata_failures"].append(f"{family_name}/{direction}: {error}")
                    continue
                tally["witnesses_metadata_ok"] += 1
                key = move.labels
                connectivity = _connectivity(witness)
                if key in classes:
                    if classes[key] != connectivity:
                        tally["class_splits"].append(f"{witness.label}")
                    continue
                classes[key] = connectivity
                tally["classes"] += 1
                graph_key = (witness.reactants, witness.products)
                if graph_key in seen_graphs:
                    tally["exact_duplicate_witnesses"] += 1
                    continue
                seen_graphs.add(graph_key)
                tally["checked_witnesses"] += 1
                check_witness(witness, tally, full_box=first)
                first = False
        if self_check:
            reference = sorted(m.sort_key() for m in rm.generate_moves(family, molecules,
                                                                        direction))
            if sorted(streamed) != reference:
                tally["failure_reasons"].append("streamed generation != resonance='all'")
                tally["failures"] += 1
            tally["self_checked"] = True
    except Exception:  # reported, never hidden
        tally["failures"] += 1
        tally["failure_reasons"].append(traceback.format_exc()[-1500:])
    tally["seconds"] = round(time.time() - start, 1)
    tally["max_rss_mb"] = _rss_mb()
    return proxy.site_type, family_name, direction, tally


def alkane_part():
    from rmgpy.kmc import boundaries as b
    from rmgpy.kmc import rmg_moves as rm

    families = _STATE["families"]
    rows = []
    for n in (8, 16):
        start = time.time()
        tally = new_tally()
        build = b.PROVIDERS["linear_alkane"].build(n)
        molecules, participants = b.witness_reactants([build])
        forms = [rm.reactant_forms(m) for m in molecules]
        moves = rm.generate_moves(families["R_Recombination"], molecules, "reverse")
        witnesses = []
        for move in moves:
            tally["statuses"][move.status] += 1
            if move.status != "kept":
                continue
            tally["kept_moves"] += 1
            witness = b.witness_from_move(move, molecules, participants, forms=forms)
            tally["witnesses_metadata_ok"] += 1
            witnesses.append(witness)
            tally["checked_witnesses"] += 1
            check_witness(witness, tally, full_box=True)
        across = b.marked_duplicates([(w.label, w) for w in witnesses])
        tally["cross_witness_duplicates"] = len(across)
        tally["seconds"] = round(time.time() - start, 1)
        rows.append((f"C{n} alkane", "R_Recombination", "reverse", tally))
    return rows


def _jsonable(tally):
    out = {}
    for key, value in tally.items():
        out[key] = dict(value) if isinstance(value, Counter) else value
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("outdir")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--self-check", action="store_true")
    parser.add_argument("--only", default=None, help="comma-separated site types")
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    started = time.time()
    _init()
    tasks = []
    only = set(args.only.split(",")) if args.only else None
    for index, proxy in enumerate(_STATE["proxies"]):
        if only and proxy.site_type not in only:
            continue
        for family_name in proxy.metadata["family_candidates"]:
            family = _STATE["families"][family_name]
            for direction in ("forward", "reverse"):
                if direction == "reverse" and family.reverse_template is None:
                    continue
                small = len(proxy.reactants) == 1 and "@5" not in proxy.site_type
                tasks.append((index, family_name, direction, args.self_check and small))
    # heaviest first: pairs, then 5-unit contexts
    tasks.sort(key=lambda t: (-len(_STATE["proxies"][t[0]].reactants),
                              "@5" not in _STATE["proxies"][t[0]].site_type, t))
    results = []
    with open(outdir / "tasks.jsonl", "w") as log:
        if args.workers > 1:
            import multiprocessing as mp

            context = mp.get_context("fork")
            with context.Pool(args.workers, initializer=_init) as pool:
                for row in pool.imap_unordered(run_task, tasks):
                    results.append(row)
                    log.write(json.dumps([row[0], row[1], row[2], _jsonable(row[3])]) + "\n")
                    log.flush()
                    print(f"[{time.time() - started:7.0f}s] {row[0]} {row[1]} {row[2]} "
                          f"kept={row[3]['kept_moves']} checked={row[3]['checked_witnesses']} "
                          f"vectors={row[3]['vectors']} fail={row[3]['failures']} "
                          f"t={row[3]['seconds']}s rss={row[3]['max_rss_mb']}MB", flush=True)
        else:
            for task in tasks:
                row = run_task(task)
                results.append(row)
                log.write(json.dumps([row[0], row[1], row[2], _jsonable(row[3])]) + "\n")
                log.flush()
                print(row[0], row[1], row[2], row[3]["kept_moves"], row[3]["seconds"],
                      flush=True)
        for row in alkane_part():
            results.append(row)
            log.write(json.dumps([row[0], row[1], row[2], _jsonable(row[3])]) + "\n")
    write_table(results, outdir, time.time() - started)


def write_table(results, outdir, seconds):
    by_proxy = {}
    for site_type, family, direction, tally in results:
        by_proxy.setdefault(site_type, []).append((family, direction, tally))
    lines = ["| proxy | kept witnesses (metadata ok) | mapped reactions | exact dup | "
             "checked | vectors (ok/refused/ceiling) | duplicates | failures |",
             "|---|---|---|---|---|---|---|---|"]
    totals = Counter()
    reasons = Counter()
    failure_lines = []
    for site_type in sorted(by_proxy, key=lambda s: ("alkane" in s, "@5" in s, s)):
        total = Counter()
        for family, direction, tally in by_proxy[site_type]:
            for key in ("kept_moves", "witnesses_metadata_ok", "classes",
                        "exact_duplicate_witnesses", "checked_witnesses", "vectors",
                        "status:ok", "status:refused", "status:ceiling_exceeded",
                        "duplicates", "failures"):
                total[key] += tally.get(key, 0)
            total["metadata_failures"] += len(tally["metadata_failures"])
            total["class_splits"] += len(tally["class_splits"])
            total["cross"] += tally.get("cross_witness_duplicates", 0)
            reasons.update(tally["refused_reasons"])
            for reason in tally["failure_reasons"][:3] + tally["metadata_failures"][:3] + \
                    tally["class_splits"][:3]:
                failure_lines.append(f"{site_type} {family}/{direction}: {reason}")
        totals.update(total)
        classes = total["classes"] if "alkane" not in site_type else total["checked_witnesses"]
        lines.append(
            f"| {site_type} | {total['kept_moves']} ({total['witnesses_metadata_ok']}) | "
            f"{classes} | {total['exact_duplicate_witnesses']} | {total['checked_witnesses']} | "
            f"{total['vectors']} ({total['status:ok']}/{total['status:refused']}/"
            f"{total['status:ceiling_exceeded']}) | {total['duplicates']}"
            + (f" (+{total['cross']} cross-witness)" if total["cross"] else "")
            + f" | {total['failures'] + total['metadata_failures'] + total['class_splits']} |")
    lines.append(
        f"| **total** | {totals['kept_moves']} ({totals['witnesses_metadata_ok']}) | "
        f"{totals['classes']} | {totals['exact_duplicate_witnesses']} | "
        f"{totals['checked_witnesses']} | {totals['vectors']} ({totals['status:ok']}/"
        f"{totals['status:refused']}/{totals['status:ceiling_exceeded']}) | "
        f"{totals['duplicates']} | {totals['failures'] + totals['metadata_failures'] + totals['class_splits']} |")
    text = "\n".join(lines) + "\n\nRefusal reasons (vectors): " + json.dumps(dict(reasons)) + \
        "\n\nFailures:\n" + ("\n".join(failure_lines) or "none") + \
        f"\n\nWall time: {seconds:.0f} s\n"
    (outdir / "table.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
