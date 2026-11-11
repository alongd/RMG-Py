#!/usr/bin/env python
"""The member-case table of an emit catalog (I-329): ``win.member_cases`` of
every subset of ``emit.roster_subsets``, computed once, split over K jobs, for
``scripts/emit_resume.py --resume --member-cases TABLE`` to read.

Nearly all of a shard's reference-case build is ``member_cases``, whose result
depends on the family, the direction and the chosen participants, not on the
shard; every shard of a catalog builds nearly the same subsets.  Here they are
computed once:

    python scripts/emit_member_cases.py build ARTIFACT ROSTER --database DB \\
        --subsets k/K --out part-k.table [--family F]... [--legacy-artifact] \\
        [--save-seconds S]
    python scripts/emit_member_cases.py merge PART... --out TABLE [--partial]

``build`` computes subsets i % K == k of the enumeration (round-robin: the cost
clusters by family), each as the shard's loop calls it, and prints one
``MEMBER i/N part k/K <seconds>s cases=c <system> <label>`` line per subset
after a ``MEMBER part k/K: m of N subsets`` line.  It saves the part every
``--save-seconds`` (600) and at the end; a rerun with the same inputs keeps a
saved part's subsets and computes the rest.  ``merge`` joins parts of one
identity into a table, and refuses parts of other identities, of another K,
the same part twice, or (without ``--partial``) a missing or incomplete part.

The table's identity (``emit_resume.table_identity``) is the code, database,
artifact, roster, systems and the subset enumeration; a shard uses a table only
under its own, and computes whatever the table lacks.  Nothing under ``rmgpy/``
is changed, so the rows header's code identity is that of the checkout.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))

import emit_resume as R  # noqa: E402

from rmgpy.kmc import emit  # noqa: E402
from rmgpy.kmc import windows as win  # noqa: E402

SAVE_SECONDS = 600.0
clock = time.monotonic


def _part(text):
    try:
        k, K = (int(n) for n in text.split("/"))
    except ValueError:
        raise argparse.ArgumentTypeError(f"{text!r} is not k/K")
    if not (K >= 1 and 0 <= k < K):
        raise argparse.ArgumentTypeError(f"{text!r} is not k/K with 0 <= k < K")
    return k, K


def build_fields(artifact, roster_sha256, database) -> dict:
    """The case identity's ``TABLE_FIELDS``, by the functions
    ``emit_resume.case_identity`` reads them with."""
    artifact, database = Path(artifact), Path(database)
    return {"artifact_id": artifact.stem,
            "artifact_sha256": emit.astream.artifact_digest(artifact)[0],
            "roster_sha256": roster_sha256, "database_head": emit.git_head(database),
            "database_content_sha256": emit.database_content_sha256(database),
            "code_sha": emit.code_sha(), "code_sources_sha256": emit.code_sources_sha256(),
            "code_extensions_sha256": emit.code_extensions_sha256()}


def build_part(args):
    """Compute part ``args.subsets`` of the table: the certifier, families
    and systems ``Emitter.init`` builds, without certifying any record."""
    part, K = args.subsets
    out = Path(args.out)
    roster = emit.load_roster(args.roster)
    art = json.loads(Path(args.artifact).read_text())
    try:
        emit.bind_artifact(art, roster.digests(), f"the roster {args.roster}",
                           args.legacy_artifact)
        emit.check_roster_for_artifact(roster, art)
    except (emit.ArtifactBindingRefusal, emit.RosterRefusal) as error:
        print("REFUSED", error, flush=True)
        return 1
    if args.family:
        art["records"] = [r for r in art["records"] if r["family"] in args.family]
    families = sorted({r["family"] for r in art["records"]})
    systems, _refusals = emit.build_reference_systems(roster)
    db = emit.load_database(families, args.database)
    fam = db.kinetics.families
    certifier = emit.make_certifier(fam, emit.certification_rules(art, roster.polymer),
                                    emit.reference_thermo(db))
    del art
    identity = R.table_identity(build_fields(args.artifact, roster.sha256, args.database),
                                fam, systems)
    unknown = [k for k in R.TABLE_REQUIRED if not identity.get(k)]
    if unknown:
        print(f"REFUSED unknown {', '.join(unknown)}: a shard never uses such a table",
              flush=True)
        return 1
    names = R.system_names(systems)
    subsets = list(emit.roster_subsets(fam, systems))
    total = len(subsets)
    mine = [i for i in range(total) if i % K == part]
    print(f"MEMBER part {part}/{K}: {len(mine)} of {total} subsets", flush=True)
    done = {}
    if out.exists():
        table, why = R.open_table(out, identity)
        if table is not None and table.meta.get("part") != [part, K]:
            why = f"another part: {table.meta.get('part')!r}"
            table.close()
            table = None
        if table is None:
            print(f"MEMBER part not resumed: {why}", flush=True)
        else:
            try:
                for key in table.entries:
                    done[key] = table.blob(key)
            except Exception as error:
                emit.reraise_resource(error)
                done = {}
                print(f"MEMBER part not resumed: {type(error).__name__}: {error}", flush=True)
            else:
                print(f"MEMBER part resumed: {len(done)} done", flush=True)
            table.close()
    out.parent.mkdir(parents=True, exist_ok=True)
    began = last = clock()
    first = True

    def save(complete):
        nonlocal first, last
        size = R.write_table(out, identity, {"part": [part, K], "complete": complete},
                             done.items(), stale=first)
        first, last = False, clock()
        return size

    for i in mine:
        name, family, direction, chosen, subset = subsets[i]
        label = f"{name}|{direction}|{subset}"
        key = R.subset_key(names, label, chosen)
        if key in done:
            continue
        t = clock()
        builds = [p.build(f"p{k}") for k, p in enumerate(chosen)]
        cases = win.member_cases(certifier, family, builds, direction,
                                 label=f"{name}|{direction}|{subset}")
        done[key] = R.dump_cases(cases)
        print(f"MEMBER {i}/{total} part {part}/{K} {clock() - t:.1f}s cases={len(cases)} "
              f"{key[0]} {label}", flush=True)
        if clock() - last >= args.save_seconds:
            save(False)
    size = save(True)
    print(f"MEMBER part {part}/{K} saved {out}: {len(done)} of {len(mine)} subsets, "
          f"{size} bytes, {clock() - began:.0f}s", flush=True)
    return 0


def merge(args):
    """Join the parts into one table (blobs streamed part to table)."""
    out = Path(args.out).resolve()
    if any(Path(p).resolve() == out for p in args.parts):
        print(f"REFUSED {args.out} is one of the parts", flush=True)
        return 1
    tables = []
    try:
        for path in args.parts:
            identity = R.peek_identity(path)
            table, why = (None, "not a member-case table") if identity is None else \
                R.open_table(path, identity)
            if table is None:
                print(f"REFUSED {path}: {why}", flush=True)
                return 1
            tables.append((path, identity, table))
        identity = tables[0][1]
        for path, other, _table in tables[1:]:
            keys = sorted(k for k in set(identity) | set(other)
                          if identity.get(k, ...) != other.get(k, ...))
            if keys:
                print(f"REFUSED {path}: another identity than {tables[0][0]}; the keys that "
                      f"differ: {', '.join(keys)}", flush=True)
                return 1
        parts = []
        for path, _identity, table in tables:
            meta = table.meta.get("part")
            if not (isinstance(meta, list) and len(meta) == 2 and
                    all(type(n) is int for n in meta)):
                print(f"REFUSED {path}: not a part (part {meta!r})", flush=True)
                return 1
            parts.append(tuple(meta))
        Ks = sorted({K for _k, K in parts})
        if len(Ks) != 1:
            print(f"REFUSED parts of different K: {', '.join(map(str, Ks))}", flush=True)
            return 1
        K = Ks[0]
        seen = set()
        for (k, _K), (path, _i, _t) in zip(parts, tables):
            if k in seen:
                print(f"REFUSED duplicate part {k}/{K}: {path}", flush=True)
                return 1
            seen.add(k)
        missing = sorted(set(range(K)) - seen)
        incomplete = [str(path) for path, _i, table in tables
                      if table.meta.get("complete") is not True]
        if (missing or incomplete) and not args.partial:
            why = []
            if missing:
                why.append(f"missing parts: {', '.join(map(str, missing))} of {K}")
            if incomplete:
                why.append(f"incomplete: {', '.join(incomplete)}")
            print(f"REFUSED {'; '.join(why)} (--partial merges them)", flush=True)
            return 1
        owner = {}
        for path, _i, table in tables:
            for key in table.entries:
                if key in owner:
                    print(f"REFUSED the subset {key!r} is in {owner[key]} and {path}",
                          flush=True)
                    return 1
                owner[key] = path

        def items():
            for _path, _i, table in tables:
                for key in table.entries:
                    yield key, table.blob(key)

        try:
            size = R.write_table(out, identity, {"part": None, "complete": not (
                missing or incomplete)}, items())
        except ValueError as error:   # an entry that does not match its digest
            print(f"REFUSED {error}", flush=True)
            return 1
    finally:
        for _path, _i, table in tables:
            table.close()
    total = identity["subsets"]
    print(f"TABLE {total} subsets, {len(owner)} present, {size} bytes: {out}", flush=True)
    if len(owner) < total:
        print(f"TABLE missing {total - len(owner)} of {total} subsets", flush=True)
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(prog="scripts/emit_member_cases.py")
    sub = parser.add_subparsers(dest="command", required=True)
    b = sub.add_parser("build", help="compute one part of the member-case table")
    b.add_argument("artifact")
    b.add_argument("roster")
    b.add_argument("--database", required=True)
    b.add_argument("--subsets", type=_part, required=True,
                   help="k/K: subsets i of the enumeration with i %% K == k")
    b.add_argument("--out", required=True)
    b.add_argument("--family", action="append", default=None,
                   help="only these families' records, as the shards' --family")
    b.add_argument(emit.LEGACY_ARTIFACT_FLAG, action="store_true")
    b.add_argument("--save-seconds", type=R._seconds, default=SAVE_SECONDS,
                   help="save the part every this many seconds (0: after every subset)")
    m = sub.add_parser("merge", help="join parts into one member-case table")
    m.add_argument("parts", nargs="+")
    m.add_argument("--out", required=True)
    m.add_argument("--partial", action="store_true",
                   help="merge even when parts are missing or incomplete (a shard computes "
                        "the subsets the table lacks)")
    args = parser.parse_args(argv)
    return build_part(args) if args.command == "build" else merge(args)


if __name__ == "__main__":
    sys.exit(main())
