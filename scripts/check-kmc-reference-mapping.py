#!/usr/bin/env python
"""Check compiled-record keys against declaration-generated reference cases.

This deliberately stops before binding or rate emission: it builds the same
reference-case index as emit, computes each record's member key, and reports
whether that key has a provider.
"""

import argparse
import json
from collections import Counter
from pathlib import Path

from rmgpy.kmc import channels as chm
from rmgpy.kmc import certificates as cm
from rmgpy.kmc import emit
from rmgpy.kmc.run import read_deck


def prior_unmapped(path):
    if path is None:
        return set()
    out = set()
    with open(path) as handle:
        for line in handle:
            row = json.loads(line)
            if row.get("kind") == "row" and row.get("refusal") == "refused:unmapped":
                out.add((row.get("payload") or {}).get("record_event_id") or row["event_id"])
    return out


def compact_systems(systems):
    """Drop systems that differ only by pair order.

    ``member_cases`` tries both template-slot orders once both species are
    present, so the checker needs one unordered pair.  The written roster
    retains both declared role orders; only this read-only index build is
    compacted.
    """
    out, seen = {}, set()
    for name, entries in systems.items():
        encoded = [json.dumps(entry, sort_keys=True, separators=(",", ":"))
                   for entry in entries]
        key = tuple(sorted(encoded)) if len(entries) == 2 else tuple(encoded)
        if key not in seen:
            seen.add(key)
            out[name] = entries
    return out


def named_refusal(record):
    status = record.get("status")
    if status and status != "enabled":
        return f"refused:record_status:{status}"
    reason = record.get("status_reason")
    return None if not reason else f"refused:record_status:{reason}"


def check(args):
    artifact = json.loads(Path(args.artifact).read_text())
    records = artifact["records"]
    if args.limit is not None:
        records = records[:args.limit]
    deck = read_deck(args.deck)
    provider = emit.match_provider(deck.spec.repeat_unit)
    roster = emit.reference_systems(deck.spec, provider, dps=args.dps)
    if args.roster_out:
        roster = emit.write_roster(roster, args.roster_out)
    compact = emit.Roster(roster.polymer, compact_systems(roster.systems),
                          roster.sha256, roster.path)
    systems, refusals = emit.build_reference_systems(compact)
    if refusals:
        raise RuntimeError(f"generated roster has refusals: {refusals}")

    families = sorted({record["family"] for record in records})
    database = emit.load_database(families, args.database)
    fam = database.kinetics.families
    certifier = cm.Certifier(fam, max_channel_sets=8, max_entries=1024)
    species = chm.SpeciesKeys()
    cases, case_stats = emit.reference_cases(certifier, fam, systems, species)

    old = prior_unmapped(args.old_rows)
    selected_ids = {record["event_id"] for record in records}
    required = old & selected_ids
    counts = Counter()
    required_missing = []
    for record in records:
        refusal = named_refusal(record)
        try:
            key = emit.record_key(certifier, record, species)
            found = bool(cases.get(key))
        except Exception as error:
            found = False
            refusal = f"refused:key:{type(error).__name__}:{error}"
        group = (record.get("site_type") or "<none>", record.get("orientation") or "<none>")
        counts[group + (("found" if found else "not_found"),)] += 1
        if not found and refusal:
            counts[group + ("named_refusal",)] += 1
        if record["event_id"] in required:
            counts[group + (("prior_unmapped_found" if found else "prior_unmapped_missing"),)] += 1
            if not found and refusal == "refused:unmapped":
                required_missing.append(record["event_id"])
            elif not found and refusal is None:
                required_missing.append(record["event_id"])

    print("site_type\torientation\tfound\tnot_found\tnamed_refusal\tprior_unmapped_found\tprior_unmapped_missing")
    groups = sorted({key[:2] for key in counts})
    for site_type, orientation in groups:
        values = [counts[(site_type, orientation, field)] for field in
                  ("found", "not_found", "named_refusal", "prior_unmapped_found",
                   "prior_unmapped_missing")]
        print("\t".join(map(str, (site_type, orientation, *values))))
    print(json.dumps({
        "records": len(records),
        "found": sum(v for key, v in counts.items() if key[2] == "found"),
        "not_found": sum(v for key, v in counts.items() if key[2] == "not_found"),
        "prior_unmapped_in_rows": len(old),
        "prior_unmapped_selected": len(required),
        "prior_unmapped_missing": required_missing,
        "case_stats": case_stats,
        "reference_systems": len(systems),
        "written_systems": len(roster.systems),
        "roster_sha256": roster.sha256,
    }, sort_keys=True))
    return 1 if required_missing else 0


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact")
    parser.add_argument("deck")
    parser.add_argument("--database", required=True)
    parser.add_argument("--old-rows")
    parser.add_argument("--roster-out")
    parser.add_argument("--dps", type=lambda value: tuple(int(x) for x in value.split(",")),
                        default=(1, 2, 3, 8))
    parser.add_argument("--limit", type=int)
    return parser.parse_args()


if __name__ == "__main__":
    raise SystemExit(check(parse_args()))
