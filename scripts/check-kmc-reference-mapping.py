#!/usr/bin/env python
"""Check compiled-record keys against declaration-generated reference cases.

This deliberately stops before binding or rate emission: it certifies the
artifact and refuses records upstream exactly as emit does, builds the same
reference-case index, computes each firing record's member key, and reports
whether that key has a provider.  Every selected record is one of: a named
refusal (production's own upstream refusal, never emitted), found, or
uncovered; the exit status is non-zero when any record is uncovered.

The deck must be the declaration the artifact was compiled from: a carried
declaration digest that differs is refused by name, and an artifact carrying
none (a compile without a deck) is refused unless ``--legacy-artifact``
(``emit.bind_artifact``, as emit binds its roster), before anything is built.
"""

import argparse
import json
import re
from collections import Counter
from pathlib import Path

from rmgpy.kmc import channels as chm
from rmgpy.kmc import certificates as cm
from rmgpy.kmc import emit
from rmgpy.kmc import template_emission as te
from rmgpy.kmc.run import read_deck


def prior_unmapped(path):
    out = set()
    with open(path) as handle:
        for line in handle:
            row = json.loads(line)
            if row.get("kind") == "row" and row.get("refusal") == "refused:unmapped":
                out.add((row.get("payload") or {}).get("record_event_id") or row["event_id"])
    return out


# what the check reads of every record (certification, upstream refusal, key)
RECORD_FIELDS = ("event_id", "family", "status", "orientation", "template", "reactant_graphs",
                 "product_graphs")


def artifact_problems(records):
    """Why the artifact cannot be checked: no records, a record lacking a
    field the check reads, or an event id carried twice (empty: none)."""
    if not records:
        return ["refused:artifact: the artifact has no records"]
    out, seen = [], set()
    for index, record in enumerate(records):
        lacking = [field for field in RECORD_FIELDS if field not in record]
        if lacking:
            out.append(f"refused:artifact: record {record.get('event_id', f'#{index}')} lacks "
                       f"{', '.join(lacking)}")
        event_id = record.get("event_id")
        if event_id in seen:
            out.append(f"refused:artifact: duplicate event id {event_id}")
        seen.add(event_id)
    return out


def named_refusal(record, store, by_id):
    """Production's upstream refusal of a certified record (None: it fires)."""
    return te.upstream_refusal(record, emit.record_certificate(record, store), by_id)


def reference_dps(records, requested):
    """Requested DPs plus numeric proxy depths recorded as ``site_type@N``."""
    depths = set(requested)
    for record in records:
        match = re.search(r"@(\d+)$", record.get("site_type") or "")
        if match:
            depths.add(int(match.group(1)))
    return tuple(sorted(depths))


def check(args):
    artifact = json.loads(Path(args.artifact).read_text())
    all_records = artifact["records"]
    problems = artifact_problems(all_records)
    if problems:
        raise RuntimeError("; ".join(problems[:5]))
    deck = read_deck(args.deck)
    bound_by = emit.bind_artifact(artifact, deck.digests(), f"the deck {args.deck}",
                                  args.legacy_artifact)
    provider = emit.match_provider(deck.spec.repeat_unit)
    roster = emit.reference_systems(deck.spec, provider, dps=reference_dps(all_records, args.dps))
    if args.roster_out:
        roster = emit.write_roster(roster, args.roster_out)
    systems, refusals = emit.build_reference_systems(roster)
    if refusals:
        raise RuntimeError(f"generated roster has refusals: {refusals}")

    # certification decides which records fire, and a reverse record's
    # requirement reads its rate owner: certify the whole artifact, as emit
    # does, whatever --limit selects
    families = sorted({record["family"] for record in all_records})
    database = emit.load_database(families, args.database)
    fam = database.kinetics.families
    certifier = cm.Certifier(fam, max_channel_sets=8, max_entries=1024)
    certified, store, replacements = emit.certified_records(all_records, fam, certifier)
    by_id = {record["event_id"]: record for record in certified}
    records = certified if args.limit is None else certified[:args.limit]
    species = chm.SpeciesKeys()
    cases, case_stats = emit.reference_cases(certifier, fam, systems, species)

    # prior rows are optional: without them there is no prior-unmapped report
    old = None
    if args.old_rows is not None:
        old = {replacements.get(event_id, event_id) for event_id in prior_unmapped(args.old_rows)}
        unknown_old = old - set(by_id)
        if not old or unknown_old:
            raise RuntimeError(f"old rows identify {len(old)} prior-unmapped records; "
                               f"{len(unknown_old)} are absent from the artifact")
    selected_ids = {record["event_id"] for record in records}
    required = old & selected_ids if old is not None else set()
    counts = Counter()
    by_status, by_reason = Counter(), Counter()
    key_errors = {}
    required_missing = []
    missing = []
    for record in records:
        group = (record.get("site_type") or "<none>", record.get("orientation") or "<none>")
        why = named_refusal(record, store, by_id)
        if why is not None:
            outcome = "named_refusal"
            by_status[te.status_of(why)] += 1
            by_reason[why] += 1
        else:
            try:
                key = emit.record_key(certifier, record, species)
                outcome = "found" if cases.get(key) else "uncovered"
            except Exception as error:   # production's row is error:exception
                outcome = "uncovered"
                key_errors[record["event_id"]] = f"{type(error).__name__}: {error}"[:300]
        counts[group + (outcome,)] += 1
        if record["event_id"] in required:
            counts[group + ("prior_unmapped_" + outcome,)] += 1
            if outcome == "uncovered":
                required_missing.append(record["event_id"])
        if outcome == "uncovered":
            missing.append(record["event_id"])

    fields = ("found", "uncovered", "named_refusal")
    if old is not None:
        fields += ("prior_unmapped_found", "prior_unmapped_uncovered",
                   "prior_unmapped_named_refusal")
    print("\t".join(("site_type", "orientation") + fields))
    for site_type, orientation in sorted({key[:2] for key in counts}):
        values = [counts[(site_type, orientation, field)] for field in fields]
        print("\t".join(map(str, (site_type, orientation, *values))))
    summary = {
        "records": len(records),
        "artifact_records": len(all_records),
        "found": sum(v for key, v in counts.items() if key[2] == "found"),
        "uncovered": sum(v for key, v in counts.items() if key[2] == "uncovered"),
        "named_refused": sum(by_status.values()),
        "named_refused_by_status": dict(by_status),
        "named_refused_by_reason": dict(by_reason),
        "key_errors": key_errors,
        "missing": missing,
        "case_stats": case_stats,
        "reference_systems": len(systems),
        "written_systems": len(roster.systems),
        "roster_sha256": roster.sha256,
        "artifact_bound_by": list(bound_by),
        "legacy_artifact": not bound_by,
    }
    if old is not None:
        summary.update(prior_unmapped_in_rows=len(old), prior_unmapped_selected=len(required),
                       prior_unmapped_missing=required_missing)
    print(json.dumps(summary, sort_keys=True))
    return 1 if missing else 0


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact")
    parser.add_argument("deck")
    parser.add_argument("--database", required=True)
    parser.add_argument("--old-rows", help="earlier emission rows; their refused:unmapped "
                        "records are reported as prior-unmapped")
    parser.add_argument("--roster-out")
    parser.add_argument("--dps", type=lambda value: tuple(int(x) for x in value.split(",")),
                        default=(1, 2, 3, 8))
    parser.add_argument("--limit", type=int)
    parser.add_argument(emit.LEGACY_ARTIFACT_FLAG, action="store_true",
                        help="check an artifact that carries no declaration digest (a compile "
                             "without --deck) against the deck, unbound")
    return parser.parse_args()


if __name__ == "__main__":
    raise SystemExit(check(parse_args()))
