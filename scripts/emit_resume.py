#!/usr/bin/env python
"""``python -m rmgpy.kmc.emit`` (one shard's rows) with ``--resume`` (I-321),
for a checkout of polymer-kmc 7a9979d44, outside the ``rmgpy`` package.

A shard killed at any point (a cluster job's wall-time limit) and run again
with ``--resume`` and the same inputs keeps the complete rows of its ``--out``,
when that file's header is this run's own, and emits only the records after
them: the result is the uninterrupted run's rows, the wall-clock fields aside.
Lines from the first incomplete or out-of-order one on are dropped (appended
to ``<out>.resume-dropped``) and the kept prefix is rewritten atomically.  A
file with another header, or no rows header, is refused by the keys that
differ and left untouched.  Without ``--resume`` it is ``emit``'s run.

It lives outside ``rmgpy/`` on purpose: the rows header and the prepared
manifest carry ``code_sources_sha256``, a digest of every source file under
``rmgpy/``, so a change there would refuse every prepared directory and every
rows file written at 7a9979d44.  This file changes no byte under ``rmgpy/``;
its ``run_emit_resume`` is that commit's ``emit.run_emit`` with the resume
branch of the polymer-kmc tip, and the helpers are the tip's, verbatim.

    python scripts/emit_resume.py ARTIFACT ROSTER --database DB --out rows.jsonl
        [--chunk i/n] [--family F]... [--limit N] [--no-contexts] [--legacy-artifact]
        [--prepared DIR] [--resume]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

from rmgpy.kmc import emit


class ResumeRefusal(ValueError):
    """A rows file ``--resume`` must not continue: written with another
    header, not a rows file, or changed while it was being resumed."""


@dataclass
class ResumePrefix:
    """What ``resume_prefix`` keeps of an existing rows file: the first
    ``kept`` picked records' rows (``errors`` of them error rows), the bytes to
    keep (``prefix``) and to drop (``dropped``), and the file's ``stat`` when
    it was read (None: no file)."""

    kept: int
    errors: int
    prefix: bytes
    dropped: bytes
    stat: Any


def _is_row_of(line: bytes, record) -> Optional[dict]:
    """The row ``line`` when it is a complete row of ``record``, else None."""
    if not line.endswith(b"\n"):
        return None
    try:
        row = json.loads(line)
    except ValueError:
        return None
    if not isinstance(row, dict) or row.get("kind") != "row":
        return None
    if (row.get("event_id"), row.get("orientation")) != (record["event_id"][:24],
                                                         record["orientation"]):
        return None
    return row


def resume_prefix(out, fixed_lines, picked) -> ResumePrefix:
    """The part of the rows file ``out`` a resumed run keeps: its header and
    roster refusal lines (``fixed_lines``, as this run writes them) and then
    the complete rows of ``picked`` in order, up to the first incomplete or
    out-of-order line.  Refuses (``ResumeRefusal``) a file whose first line is
    not this run's header, naming the header keys that differ.  An absent or
    empty file, or one holding only the start of this run's header, keeps
    nothing."""
    path = Path(out)
    fixed = [line.encode() for line in fixed_lines]
    try:
        stat = path.stat()
    except FileNotFoundError:
        return ResumePrefix(0, 0, b"".join(fixed), b"", None)
    data = path.read_bytes()
    if not data:
        return ResumePrefix(0, 0, b"".join(fixed), b"", stat)
    if not data.endswith(b"\n") and b"\n" not in data and fixed[0].startswith(data):
        # this run's own header cut short: a long header is on disk before the
        # first row flushes the rest
        return ResumePrefix(0, 0, b"".join(fixed), data, stat)
    lines = data.splitlines(keepends=True)
    if lines[0] != fixed[0]:
        try:
            old = json.loads(lines[0]) if lines[0].endswith(b"\n") else None
        except ValueError:
            old = None
        if not isinstance(old, dict) or old.get("kind") != "header":
            raise ResumeRefusal(f"{out}: its first line is not a complete rows header")
        new = json.loads(fixed[0])
        keys = sorted(k for k in set(old) | set(new)
                      if k not in old or k not in new or old[k] != new[k])
        raise ResumeRefusal(f"{out} was written with another header; the keys that differ: "
                            f"{', '.join(keys) or '(none: the same values, other bytes)'}")
    kept, errors, at = 0, 0, len(fixed)
    if lines[:len(fixed)] != fixed:
        at = 1     # a roster refusal line cut short: no row is kept
    else:
        for line in lines[at:]:
            row = kept < len(picked) and _is_row_of(line, picked[kept])
            if not row:
                break
            kept += 1
            errors += row.get("status") == "error"
        at += kept
    prefix = b"".join(fixed) + b"".join(lines[len(fixed):len(fixed) + kept])
    return ResumePrefix(kept, errors, prefix, b"".join(lines[at:]), stat)


def rewrite_prefix(out, resume: ResumePrefix) -> None:
    """Make ``out`` hold exactly ``resume.prefix``: unchanged when it already
    does; else, after checking the file is as ``resume_prefix`` read it, the
    dropped bytes appended to ``<out>.resume-dropped`` and the prefix written
    to a temporary file moved onto ``out`` (atomic)."""
    import tempfile
    path = Path(out)
    if resume.stat is not None and not resume.dropped and \
            resume.stat.st_size == len(resume.prefix):
        return
    try:
        now = path.stat()
    except FileNotFoundError:
        now = None
    same = (now is None and resume.stat is None) or (
        now is not None and resume.stat is not None and
        (now.st_ino, now.st_size, now.st_mtime_ns) ==
        (resume.stat.st_ino, resume.stat.st_size, resume.stat.st_mtime_ns))
    if not same:
        raise ResumeRefusal(f"{out} changed while it was being resumed")
    if resume.dropped:
        with open(path.with_name(path.name + ".resume-dropped"), "ab") as dropped:
            dropped.write(resume.dropped)
            dropped.flush()
            os.fsync(dropped.fileno())
    fd, temp = tempfile.mkstemp(dir=path.parent, prefix=path.name + ".", suffix=".resume-tmp")
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(resume.prefix)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp, path)
    except BaseException:
        if os.path.exists(temp):
            os.unlink(temp)
        raise


def run_emit_resume(args):
    roster = emit.load_roster(args.roster)
    prepared_dir = getattr(args, "prepared", None)
    prepared = None
    # the roster's declaration must be the artifact's, before anything is built
    try:
        if prepared_dir is None:
            artifact = json.loads(Path(args.artifact).read_text())
        else:
            from rmgpy.kmc import emit_prepare
            manifest = emit_prepare.open_manifest(prepared_dir, args.artifact, roster,
                                                  args.database, args.chunk, args.family)
            artifact = manifest["head"]
        bound_by = emit.bind_artifact(artifact, roster.digests(), f"the roster {args.roster}",
                                      args.legacy_artifact)
        emit.check_roster_for_artifact(roster, artifact)
        if prepared_dir is not None:
            prepared = emit_prepare.read_shard(prepared_dir, manifest, args.chunk)
    except (emit.ArtifactBindingRefusal, emit.RosterRefusal, emit.PreparedRefusal) as error:
        print("REFUSED", error, flush=True)
        return 1
    print("BINDING", {"bound_by": list(bound_by), "legacy_artifact": not bound_by}, flush=True)
    # the partition: the artifact's own record ids (of the selected families),
    # read before anything this shard computes
    if prepared is None:
        artifact_ids = [r["event_id"] for r in artifact.get("records", [])
                        if not args.family or r.get("family") in args.family]
    else:
        artifact_ids = prepared.artifact_ids
    del artifact
    systems, refusals = emit.build_reference_systems(roster)
    for refusal in refusals:
        print("ROSTER", json.dumps(refusal), flush=True)
    emitter = emit.Emitter(args.artifact, args.database, systems, families=args.family,
                           contexts=not args.no_contexts, polymer=roster.polymer)
    t0 = time.time()
    if prepared is None:
        emitter.init()
        picked = emit.select_records(emitter.records, args.chunk, args.limit,
                                     artifact_ids=artifact_ids,
                                     replacements=emitter.state.get("replacements"))
        selected = len(emitter.records)
    else:
        emitter.init(prepared)
        picked = prepared.picked[:args.limit] if args.limit else prepared.picked
        selected = len(artifact_ids)
        del prepared
    header = emit.shard_header(emitter, roster, refusals, picked, args.chunk,
                               partial=emit.is_partial_selection(
                                   args.family, args.limit, emitter.state["artifact_records"],
                                   selected),
                               bound_by=bound_by, artifact_ids=artifact_ids)
    print("INIT", {k: emitter.state.get(k) for k in ("load_s", "certify_s", "cases_s", "case_stats")},
          len(emitter.records), "records, selected", len(picked), flush=True)
    errors = 0
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    start, mode = 0, "w"
    if getattr(args, "resume", False):
        fixed = [json.dumps(header, default=str) + "\n"] + \
            [json.dumps(refusal) + "\n" for refusal in refusals]
        try:
            resume = resume_prefix(out, fixed, picked)
            rewrite_prefix(out, resume)
        except ResumeRefusal as error:
            print("REFUSED resume:", error, flush=True)
            return 1
        print(f"RESUME kept={resume.kept} of {len(picked)}", flush=True)
        start, mode, errors = resume.kept, "a", resume.errors
    with open(out, mode) as log:
        if mode == "w":
            log.write(json.dumps(header, default=str) + "\n")
            for refusal in refusals:
                log.write(json.dumps(refusal) + "\n")
        for k, record in enumerate(picked[start:], start):
            row = emitter.run_record(record)
            errors += row["status"] == "error"
            log.write(json.dumps(row, default=str) + "\n")
            log.flush()
            print(f"[{time.time() - t0:6.0f}s] {k + 1}/{len(picked)} {row['event_id'][:16]} "
                  f"{row['family']} {row['orientation']} -> {row['status']} "
                  f"{row.get('refusal') or ''} t={row['seconds']}s rss={row['rss_mb']}MB",
                  flush=True)
    return 1 if errors else 0


def main(argv=None):
    parser = argparse.ArgumentParser(prog="scripts/emit_resume.py")
    parser.add_argument("artifact")
    parser.add_argument("roster")
    parser.add_argument("--database", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--chunk", default=None, help="i/n: records i::n of the sorted list")
    parser.add_argument("--family", action="append", default=None,
                        help="only these families' records (repeatable; a partial catalog unless "
                             "they are all the artifact's records)")
    parser.add_argument("--limit", type=int, default=None, help="at most N records (a partial catalog)")
    parser.add_argument("--no-contexts", action="store_true")
    parser.add_argument(emit.LEGACY_ARTIFACT_FLAG, action="store_true",
                        help="emit an artifact that carries no declaration digest (a compile "
                             "without --deck) against the roster's declaration, unbound")
    parser.add_argument("--prepared", default=None,
                        help="the directory emit prepare wrote: read this shard's certified "
                             "records from it instead of parsing and certifying the artifact")
    parser.add_argument("--resume", action="store_true",
                        help="keep the complete rows of an existing --out written with this "
                             "run's header and emit only the records after them; refuse a "
                             "file with another header")
    return run_emit_resume(parser.parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())
