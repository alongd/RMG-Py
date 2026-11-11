#!/usr/bin/env python
"""``python -m rmgpy.kmc.emit`` (one shard's rows) with ``--resume`` (I-321,
I-324), for a checkout of polymer-kmc 119e4e477, outside the ``rmgpy`` package.

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
rows file written at 119e4e477.  This file changes no byte under ``rmgpy/``;
its ``run_emit_resume`` is that commit's ``emit.run_emit`` with the resume
branch, statement for statement (``emitResumeTest`` pins it), and reads a
serial prepared directory or a block-prepared one (``--prepared RUN/current``)
through that commit's ``emit_prepare``.

    python scripts/emit_resume.py ARTIFACT ROSTER --database DB --out rows.jsonl
        [--chunk i/n] [--family F]... [--limit N] [--no-contexts] [--legacy-artifact]
        [--prepared DIR] [--resume]

With ``--resume`` the reference cases (``Emitter.init``'s ``reference_cases``,
``cases_s``: 13 h for one cluster shard) are saved, with their stats and the
species registry they filled, to ``<out>.cases.pickle`` before the first row
(I-326), and a later ``--resume`` of the same shard loads them instead of
building them again.  The sidecar's first line is its identity: the rows
header fields that bind the run (artifact, roster, database, code, chunk), the
prepared manifest's digest and the shard's wanted case keys, then the payload's
length and sha256.  A sidecar whose identity is not this run's, or whose
payload does not match its digest, is moved to ``<out>.cases.pickle.stale``
and the cases are built again; a run with an unknown identity field neither
saves nor loads one.  Without ``--resume`` no sidecar is read or written.  The
cache wraps ``Emitter.init`` from ``main`` and leaves ``run_emit_resume``
statement for statement as it was.

A shard whose case build is cut before it ends (all of one day's job) keeps
its progress (I-327): with ``--resume`` the cases are built by
``checkpointed_reference_cases``, ``emit.reference_cases`` statement for
statement (``emitResumeTest`` pins it), which saves the subsets done so far,
their cases and stats and the species registry to ``<out>.cases.checkpoint``
every ``--case-checkpoint-seconds`` (600) of build time, under the sidecar's
identity, and prints one ``CASES checkpoint i/n subsets`` line per save.  A
later ``--resume`` with a checkpoint of its identity continues from its
subset and ends with the uninterrupted build's cases; any other checkpoint is
moved to ``<out>.cases.checkpoint.stale`` and the build starts from its first
subset.  Once the sidecar is saved the checkpoint is removed.
"""
from __future__ import annotations

import argparse
import contextlib
import copyreg
import dataclasses
import hashlib
import io
import json
import math
import os
import pickle
import sys
import tempfile
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

from rmgpy.kmc import emit
from rmgpy.kmc import windows as win


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


CASES_FORMAT = "emit-resume-cases/1"
CASES_SUFFIX = ".cases.pickle"
# identity fields without which a sidecar cannot be told from another run's
CASES_REQUIRED = ("artifact_sha256", "roster_sha256", "database_content_sha256",
                  "code_sources_sha256", "code_extensions_sha256")
# the payload fields of the identity line, filled in after the payload
PAYLOAD_KEYS = ("payload_bytes", "payload_sha256")


def cases_path(out) -> Path:
    out = Path(out)
    return out.with_name(out.name + CASES_SUFFIX)


def _digest(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     default=str).encode()).hexdigest()


def case_identity(emitter, args, roster_sha256, wanted) -> dict:
    """What the reference cases of this run depend on, by the rows header's
    own functions and values, read when the cases are built: the header
    fields that bind the run, the prepared manifest's bytes and the shard's
    wanted keys (``case_filter``; None on the whole-artifact path)."""
    S = emitter.state
    artifact = emitter.artifact
    prepared = getattr(args, "prepared", None)
    manifest = None
    if prepared is not None:
        from rmgpy.kmc import emit_prepare
        manifest = hashlib.sha256((Path(prepared) / emit_prepare.MANIFEST)
                                  .read_bytes()).hexdigest()
    return {"format": CASES_FORMAT, "artifact_id": artifact.stem,
            "artifact_sha256": S.get("artifact_sha256") or
            emit.astream.artifact_digest(artifact)[0],
            "artifact_records": S["artifact_records"], "families": S["families"],
            "family": args.family,
            "database_path": str(emitter.database),
            "database_head": emit.git_head(emitter.database),
            "database_content_sha256": emit.database_content_sha256(emitter.database),
            "code_sha": emit.code_sha(), "code_sources_sha256": emit.code_sources_sha256(),
            "code_extensions_sha256": emit.code_extensions_sha256(),
            "roster_sha256": roster_sha256,
            "systems": {n: [p.label for p in ps] for n, ps in emitter.systems.items()},
            "chunk": args.chunk, "contexts": not args.no_contexts,
            "prepared_manifest_sha256": manifest,
            "wanted_sha256": None if wanted is None else _digest(
                [sorted(wanted.keys, key=repr),
                 sorted(([k, sorted(v, key=repr)] for k, v in wanted.sides.items()), key=repr)])}


class _CasesPickler(pickle.Pickler):
    """Writes a ``boundaries.Witness`` as its dataclass fields only: the views
    it memoises on itself (``reactant_rows`` ..., mappingproxies that cannot be
    pickled) are derived, and rebuilt on first use after loading."""

    def reducer_override(self, obj):
        from rmgpy.kmc.boundaries import Witness
        if type(obj) is Witness:
            return (copyreg.__newobj__, (Witness,),
                    {f.name: obj.__dict__[f.name] for f in dataclasses.fields(Witness)})
        return NotImplemented


class _Hashing(io.RawIOBase):
    def __init__(self, handle):
        self.handle, self.sha, self.size = handle, hashlib.sha256(), 0

    def writable(self):
        return True

    def write(self, data):
        self.sha.update(data)
        self.size += len(data)
        return self.handle.write(data)


def _identity_line(identity, size, sha) -> bytes:
    return (json.dumps({**identity, "payload_bytes": f"{size:020d}", "payload_sha256": sha},
                       sort_keys=True, default=str) + "\n").encode()


def write_cases(path, identity, payload, stale=True) -> None:
    """Write the sidecar ``path`` atomically: its identity line and the
    pickled ``payload``, through a temporary file moved onto ``path``; a
    sidecar already there is moved to ``<path>.stale`` first (``stale``; else
    it is replaced: a checkpoint by its own successor)."""
    path = Path(path)
    fd, temp = tempfile.mkstemp(dir=path.parent, prefix=path.name + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "w+b") as handle:
            head = _identity_line(identity, 0, "0" * 64)
            handle.write(head)
            hashing = _Hashing(handle)
            _CasesPickler(hashing, protocol=pickle.HIGHEST_PROTOCOL).dump(payload)
            line = _identity_line(identity, hashing.size, hashing.sha.hexdigest())
            assert len(line) == len(head)
            handle.seek(0)
            handle.write(line)
            handle.flush()
            os.fsync(handle.fileno())
        if stale and path.exists():
            os.replace(path, path.with_name(path.name + ".stale"))
        os.replace(temp, path)
    except BaseException:
        if os.path.exists(temp):
            os.unlink(temp)
        raise


def read_cases(path, identity):
    """(payload, None) when ``path`` is a sidecar of ``identity`` whose payload
    matches its digest, else (None, why); never raises on a bad sidecar.  The
    digest checked is of the bytes read from the same open file."""
    try:
        with open(path, "rb") as handle:
            line = handle.readline()
            try:
                head = json.loads(line) if line.endswith(b"\n") else None
            except ValueError:
                head = None
            if not isinstance(head, dict):
                return None, "its first line is not an identity"
            stored = {k: v for k, v in head.items() if k not in PAYLOAD_KEYS}
            expected = json.loads(json.dumps(identity, default=str))
            if stored != expected:
                keys = sorted(k for k in set(stored) | set(expected)
                              if stored.get(k, ...) != expected.get(k, ...))
                return None, f"another identity; the keys that differ: {', '.join(keys)}"
            try:
                size, sha = int(head["payload_bytes"]), head["payload_sha256"]
            except (KeyError, TypeError, ValueError):
                return None, "no payload digest"
            start, digest, read = handle.tell(), hashlib.sha256(), 0
            for block in iter(lambda: handle.read(1 << 20), b""):
                digest.update(block)
                read += len(block)
            if (read, digest.hexdigest()) != (size, sha):
                return None, f"its payload ({read} bytes) does not match its digest"
            handle.seek(start)
            return pickle.load(handle), None
    except FileNotFoundError:
        return None, "no sidecar"
    except Exception as error:   # a bad sidecar is rebuilt, never fatal
        emit.reraise_resource(error)
        return None, f"unreadable: {type(error).__name__}: {error}"


CHECKPOINT_SUFFIX = ".cases.checkpoint"
CHECKPOINT_SECONDS = 600.0


def checkpoint_path(out) -> Path:
    out = Path(out)
    return out.with_name(out.name + CHECKPOINT_SUFFIX)


class CaseCheckpoint:
    """The reference-case build's progress, saved to ``path`` under the
    sidecar's ``identity`` (``write_cases``: atomic, digested) after each
    subset once ``interval`` seconds of build time have passed since the last
    save: the index of the next subset, the cases and stats so far and the
    species registry.  The ``formulas`` memo is not saved: it is keyed by
    ``id`` of this process's participants (another process's ids could name
    other participants) and is recomputed on demand."""

    clock = staticmethod(time.monotonic)

    def __init__(self, path, identity, interval):
        self.path, self.identity, self.interval = Path(path), identity, interval
        self.start = self.total = 0
        self.began = self.last = None

    def restore(self, fam, systems, out, stats, species) -> int:
        """The index of the first subset to build: a checkpoint of this
        identity, whose digest and subset count match, fills ``out``,
        ``stats`` and ``species``; any other file is moved to ``<path>.stale``
        and the build starts from 0."""
        self.total = sum(1 for _ in emit.roster_subsets(fam, systems))
        start = 0
        if self.path.exists():
            payload, why = read_cases(self.path, self.identity)
            if payload is not None:
                why = self._invalid(payload)
            if why is None:
                out.update(payload["out"])
                stats.update(payload["stats"])
                species.__dict__.update(payload["species"])
                start = payload["next"]
                print(f"CASES checkpoint resumed at {start}/{self.total} subsets", flush=True)
            else:
                try:
                    os.replace(self.path, self.path.with_name(self.path.name + ".stale"))
                except OSError as error:
                    why += f" (not moved: {error})"
                print(f"CASES checkpoint ignored: {why}", flush=True)
        self.start = start
        print(f"CASES build start {start}/{self.total} subsets", flush=True)
        self.began = self.last = self.clock()
        return start

    def _invalid(self, payload) -> Optional[str]:
        if not isinstance(payload, dict) or set(payload) != {"next", "subsets", "out", "stats",
                                                               "species"}:
            return "not a checkpoint payload"
        if payload["subsets"] != self.total:
            return f"subsets: {payload['subsets']!r}, this build has {self.total}"
        index = payload["next"]
        if type(index) is not int or not 0 <= index <= self.total:
            return f"next: {index!r} is not a subset index of {self.total}"
        if not (isinstance(payload["out"], dict) and isinstance(payload["stats"], dict) and
                isinstance(payload["species"], dict)):
            return "not a checkpoint payload"
        return None

    def save(self, index, out, stats, species) -> None:
        """Save before subset ``index`` (``index`` subsets done) when the
        interval has passed and something was built since the start."""
        now = self.clock()
        if index <= self.start or now - self.last < self.interval:
            return
        payload = {"next": index, "subsets": self.total, "out": dict(out), "stats": dict(stats),
                   "species": species.__dict__}
        try:
            write_cases(self.path, self.identity, payload, stale=False)
        except Exception as error:   # the build does not depend on the save
            emit.reraise_resource(error)
            print(f"CASES checkpoint not saved at {index}/{self.total}: "
                  f"{type(error).__name__}: {error}", flush=True)
        else:
            done = self.clock()
            print(f"CASES checkpoint {index}/{self.total} subsets, build {now - self.began:.0f}s, "
                  f"{self.path.stat().st_size} bytes, write {done - now:.2f}s", flush=True)
        self.last = self.clock()

    def finish(self) -> None:
        """Remove the checkpoint: the saved sidecar holds the whole build (a
        checkpoint left behind is never read while the sidecar loads)."""
        try:
            self.path.unlink()
        except FileNotFoundError:
            pass
        except OSError as error:
            print(f"CASES checkpoint not removed: {error}", flush=True)


def checkpointed_reference_cases(certifier, fam, systems, species, wanted=None, *, checkpoint):
    """``emit.reference_cases`` statement for statement (``emitResumeTest``
    pins it), with a ``CaseCheckpoint``: the build starts at the subset the
    checkpoint restored and saves before each subset."""
    out = defaultdict(list)
    stats = Counter()
    formulas = {}
    if wanted is not None:
        names = {key[0] for key in wanted.keys}
        fam = {name: family for name, family in fam.items() if name in names}
    start = checkpoint.restore(fam, systems, out, stats, species)
    for index, (name, family, direction, chosen, subset) in enumerate(
            emit.roster_subsets(fam, systems)):
        if index < start:
            continue
        checkpoint.save(index, out, stats, species)
        if wanted is not None:
            if not wanted.admits(name, direction, emit.subset_formulas(chosen, formulas)):
                stats["subsets_skipped"] += 1
                continue
            stats["subsets_enumerated"] += 1
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
            if wanted is not None and key not in wanted.keys:
                stats["cases_dropped"] += 1
                continue
            out[key].append((case, specs))
    return out, dict(stats)


# this tree's own function: the copy above stands in for it, and only for it
_REFERENCE_CASES = emit.reference_cases


@contextlib.contextmanager
def cached_cases(args):
    """For one ``run_emit_resume`` call: ``Emitter.init`` builds its reference
    cases through the sidecar of ``args.out`` (load when it is this run's,
    else build and save).  Both patched names are restored on exit."""
    init, build = emit.Emitter.init, emit.reference_cases
    path = cases_path(args.out)
    roster_sha256 = emit.load_roster(args.roster).sha256
    interval = getattr(args, "case_checkpoint_seconds", CHECKPOINT_SECONDS)

    def cached_init(self, prepared=None):

        def reference_cases(certifier, fam, systems, species, wanted=None):
            identity = case_identity(self, args, roster_sha256, wanted)
            unknown = [k for k in CASES_REQUIRED if not identity.get(k)]
            if unknown:
                print(f"CASES not cached: unknown {', '.join(unknown)}", flush=True)
                return build(certifier, fam, systems, species, wanted=wanted)
            payload, why = read_cases(path, identity)
            if payload is not None:
                species.__dict__.update(payload["species"])
                print(f"CASES loaded {path}", flush=True)
                return payload["cases"], payload["stats"]
            print(f"CASES rebuilt: {why}", flush=True)
            checkpoint = None
            if build is _REFERENCE_CASES:
                path.parent.mkdir(parents=True, exist_ok=True)
                checkpoint = CaseCheckpoint(checkpoint_path(args.out), identity, interval)
                cases, stats = checkpointed_reference_cases(certifier, fam, systems, species,
                                                            wanted=wanted, checkpoint=checkpoint)
            else:   # a replaced reference_cases (a test's): no checkpoint
                cases, stats = build(certifier, fam, systems, species, wanted=wanted)
            try:
                path.parent.mkdir(parents=True, exist_ok=True)
                write_cases(path, identity, {"cases": cases, "stats": stats,
                                             "species": species.__dict__})
            except Exception as error:   # the rows do not depend on the save
                emit.reraise_resource(error)
                print(f"CASES not saved: {type(error).__name__}: {error}", flush=True)
            else:
                print(f"CASES saved {path}", flush=True)
                if checkpoint is not None:
                    checkpoint.finish()
            return cases, stats

        emit.reference_cases = reference_cases
        try:
            return init(self, prepared)
        finally:
            emit.reference_cases = build

    emit.Emitter.init = cached_init
    try:
        yield
    finally:
        emit.Emitter.init, emit.reference_cases = init, build


def _seconds(text) -> float:
    try:
        value = float(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{text!r} is not a number of seconds")
    if not math.isfinite(value) or value < 0:
        raise argparse.ArgumentTypeError(f"{text!r} is not a finite number of seconds >= 0")
    return value


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
    parser.add_argument("--case-checkpoint-seconds", type=_seconds, default=CHECKPOINT_SECONDS,
                        help="with --resume: save the reference-case build's progress to "
                             "<out>.cases.checkpoint every this many seconds of build time "
                             "(0: after every subset); a later --resume continues from it")
    args = parser.parse_args(argv)
    if not args.resume:
        return run_emit_resume(args)
    with cached_cases(args):
        return run_emit_resume(args)


if __name__ == "__main__":
    sys.exit(main())
