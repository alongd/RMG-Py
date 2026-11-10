"""A killed emit shard resumes from its first record without a row (I-321, I-324).

An emit shard of a large catalog needs days, and a cluster job is capped at
24 h, so one shard runs as a chain of jobs.  ``--resume`` keeps the complete
rows of an existing rows file whose header is this run's own, drops anything
from the first incomplete or out-of-order line on (saved beside the file), and
emits the remaining records:

- stopped after any record, with or without a cut last line or a stray row,
  the resumed file is the uninterrupted one apart from the wall-clock fields;
- a file written with another header (chunk, roster, artifact) or not a rows
  file at all is refused by the keys that differ, and left untouched; so is a
  file that changes while it is being resumed;
- without ``--resume`` the rows file is overwritten, as before;
- the script's ``run_emit_resume`` is this tree's ``emit.run_emit`` plus the
  resume branch, line for line.

At 119e4e477 the resuming run is ``scripts/emit_resume.py``; every
uninterrupted reference run is ``emit.main``, that commit's own code.  Each
shard test runs on the whole artifact, on a serial ``prepare`` directory and
on a directory prepared in two blocks and finished (``prepare-plan``,
``prepare-block``, ``prepare-finish``; ``--prepared RUN/current``).

Emission is the digest stub of ``emitShardMemoryTest`` on the PE toy subset,
with wall-clock fields that differ per call, so the comparison strips them for
a reason; the last test emits for real against the database.
"""
import difflib
import importlib.util
import inspect
import json
import os
import sys
from pathlib import Path

import pytest

from rmgpy.kmc import emit

from emitShardMemoryTest import (DATABASE, _prepare, artifact, digest_row,  # noqa: F401
                                 needs_database, roster_and_doc, stubbed)
from emitPrepareBlocksTest import _blocks

# the resuming entry point and the module that holds its helpers: at
# 7a9979d44 the script outside ``rmgpy/`` (whose sources the header digests)
SCRIPT = Path(__file__).resolve().parents[3] / "scripts" / "emit_resume.py"
_spec = importlib.util.spec_from_file_location("emit_resume", SCRIPT)
resumer = importlib.util.module_from_spec(_spec)
sys.modules["emit_resume"] = resumer     # its dataclass resolves annotations there
_spec.loader.exec_module(resumer)
resume_main = resumer.main
WALL_CLOCK = ("seconds", "cpu", "rss_mb")


class Killed(BaseException):
    """A job killed in the middle of a shard (not an error a row reports)."""


@pytest.fixture
def clocked(monkeypatch, stubbed):
    """The stub emission, with wall-clock fields that change on every call
    and a stop after ``stop_after`` records."""
    state = {"calls": 0, "stop_after": None, "ran": []}

    def run(self, record):
        if state["stop_after"] is not None and len(state["ran"]) >= state["stop_after"]:
            raise Killed()
        state["calls"] += 1
        state["ran"].append(record["event_id"][:24])
        row = digest_row(self, record)
        row.update(seconds=float(state["calls"]), cpu=0.5 * state["calls"],
                   rss_mb=1000 + state["calls"])
        return row

    monkeypatch.setattr(emit.Emitter, "run_record", run)
    return state


def _args(artifact, roster, out, flags):
    return [str(artifact), str(roster), "--database", str(DATABASE), "--out", str(out), *flags]


def _reference(artifact, roster, out, flags):
    assert emit.main(_args(artifact, roster, out, flags)) == 0
    return out.read_bytes()


def _stopped(clocked, artifact, roster, out, flags, k):
    clocked["stop_after"], clocked["ran"] = k, []
    with pytest.raises(Killed):
        emit.main(_args(artifact, roster, out, flags))
    clocked["stop_after"], clocked["ran"] = None, []


def _resume(artifact, roster, out, flags):
    return resume_main(_args(artifact, roster, out, [*flags, "--resume"]))


def _stripped(data: bytes):
    """The lines of a rows file without the fields that measure the run: a
    row's ``seconds``, ``cpu`` and ``rss_mb`` and each member row's
    ``seconds``."""
    rows = []
    for line in data.decode().splitlines():
        row = json.loads(line)
        for key in WALL_CLOCK:
            row.pop(key, None)
        for member in row.get("member_rows") or ():
            member.pop("seconds", None)
        rows.append(row)
    return rows


def _rows(data: bytes):
    return [json.loads(line) for line in data.decode().splitlines()][1:]


@pytest.fixture(params=["whole", "prepared", "blocks"])
def shard(request, tmp_path, artifact, roster_and_doc, clocked):
    """(artifact, roster, flags) of shard 0/2, on the whole-artifact path, a
    serial prepared directory and one prepared in two blocks and finished
    (``RUN/current``), and its uninterrupted rows."""
    roster = roster_and_doc[0]
    flags = ["--chunk", "0/2"]
    if request.param == "prepared":
        flags += ["--prepared", str(_prepare(artifact, roster, tmp_path, 2))]
    elif request.param == "blocks":
        flags += ["--prepared", str(_blocks(artifact, roster, tmp_path, 2, shards=2) / "current")]
    whole = _reference(artifact, roster, tmp_path / "uninterrupted.jsonl", flags)
    clocked["ran"] = []
    return artifact, roster, flags, whole


def _picked(whole):
    return len(_rows(whole))


@pytest.mark.parametrize("k", [0, 1, 3])
def test_resume_after_a_stop_is_the_uninterrupted_rows(tmp_path, shard, clocked, capsys, k):
    artifact, roster, flags, whole = shard
    out = tmp_path / "rows.jsonl"
    _stopped(clocked, artifact, roster, out, flags, k)
    assert len(_rows(out.read_bytes())) == k
    kept = out.read_bytes()
    capsys.readouterr()
    assert _resume(artifact, roster, out, flags) == 0
    resumed = out.read_bytes()
    assert f"RESUME kept={k} of {_picked(whole)}" in capsys.readouterr().out
    # the kept rows are the bytes the first job wrote, the rest are new
    assert resumed.startswith(kept)
    assert resumed != whole, "the wall-clock fields must differ, or the strip proves nothing"
    assert _stripped(resumed) == _stripped(whole)
    assert len(clocked["ran"]) == _picked(whole) - k
    assert not (tmp_path / "rows.jsonl.resume-dropped").exists()


def test_resume_after_a_cut_last_line(tmp_path, shard, clocked, capsys):
    artifact, roster, flags, whole = shard
    out = tmp_path / "rows.jsonl"
    _stopped(clocked, artifact, roster, out, flags, 4)
    data = out.read_bytes()
    last = data.rstrip(b"\n").rfind(b"\n") + 1
    fragment = data[last:last + (len(data) - last) // 2]
    out.write_bytes(data[:last] + fragment)
    capsys.readouterr()
    assert _resume(artifact, roster, out, flags) == 0
    assert f"RESUME kept=3 of {_picked(whole)}" in capsys.readouterr().out
    assert _stripped(out.read_bytes()) == _stripped(whole)
    assert (tmp_path / "rows.jsonl.resume-dropped").read_bytes() == fragment


def test_a_row_out_of_order_ends_the_kept_prefix(tmp_path, shard, clocked, capsys):
    artifact, roster, flags, whole = shard
    out = tmp_path / "rows.jsonl"
    _stopped(clocked, artifact, roster, out, flags, 3)
    stray = whole.decode().splitlines(keepends=True)[-1].encode()
    after = b'{"kind": "row", "event_id": "x"}\n'
    out.write_bytes(out.read_bytes() + stray + after)
    capsys.readouterr()
    assert _resume(artifact, roster, out, flags) == 0
    assert f"RESUME kept=3 of {_picked(whole)}" in capsys.readouterr().out
    assert _stripped(out.read_bytes()) == _stripped(whole)
    assert (tmp_path / "rows.jsonl.resume-dropped").read_bytes() == stray + after


def test_a_finished_file_resumes_to_itself(tmp_path, shard, clocked, capsys):
    artifact, roster, flags, whole = shard
    out = tmp_path / "rows.jsonl"
    out.write_bytes(whole)
    before = os.stat(out)
    capsys.readouterr()
    assert _resume(artifact, roster, out, flags) == 0
    n = _picked(whole)
    assert f"RESUME kept={n} of {n}" in capsys.readouterr().out
    assert out.read_bytes() == whole and clocked["ran"] == []
    assert os.stat(out).st_ino == before.st_ino


@pytest.mark.parametrize("empty", [False, True])
def test_resume_without_rows_is_a_fresh_run(tmp_path, shard, clocked, capsys, empty):
    artifact, roster, flags, whole = shard
    out = tmp_path / "rows.jsonl"
    if empty:
        out.write_bytes(b"")
    capsys.readouterr()
    assert _resume(artifact, roster, out, flags) == 0
    assert f"RESUME kept=0 of {_picked(whole)}" in capsys.readouterr().out
    assert _stripped(out.read_bytes()) == _stripped(whole)


def test_a_cut_header_of_this_run_keeps_no_row(tmp_path, shard, clocked, capsys):
    """A long header (the cluster's holds every artifact id) reaches the disk
    before the first row flushes the rest: a job killed then leaves the start
    of this run's own header, which is resumed from nothing, not refused."""
    artifact, roster, flags, whole = shard
    out = tmp_path / "rows.jsonl"
    head = whole.split(b"\n")[0]
    out.write_bytes(head[:len(head) // 2])
    capsys.readouterr()
    assert _resume(artifact, roster, out, flags) == 0
    assert f"RESUME kept=0 of {_picked(whole)}" in capsys.readouterr().out
    assert _stripped(out.read_bytes()) == _stripped(whole)
    assert (tmp_path / "rows.jsonl.resume-dropped").read_bytes() == head[:len(head) // 2]


def test_a_cut_refusal_line_keeps_no_row(tmp_path, artifact, roster_and_doc, clocked, capsys,
                                         monkeypatch):
    """Roster refusal lines sit between the header and the rows; one cut short
    keeps nothing after the header, and the refusals are written again."""
    roster = roster_and_doc[0]
    refusal = {"kind": "roster_refusal", "system": "s", "status": "refused:test"}
    monkeypatch.setattr(emit, "build_reference_systems", lambda roster: ({}, [refusal]))
    flags = ["--chunk", "0/2"]
    whole = _reference(artifact, roster, tmp_path / "uninterrupted.jsonl", flags)
    out = tmp_path / "rows.jsonl"
    head, line = whole.decode().splitlines(keepends=True)[:2]
    assert json.loads(line) == refusal
    out.write_bytes((head + line[:10]).encode())
    capsys.readouterr()
    assert _resume(artifact, roster, out, flags) == 0
    assert f"RESUME kept=0 of {len(whole.splitlines()) - 2}" in capsys.readouterr().out
    assert _stripped(out.read_bytes()) == _stripped(whole)


def _untouched(out, data, stat):
    after = os.stat(out)
    assert out.read_bytes() == data
    assert (after.st_ino, after.st_mtime_ns) == (stat.st_ino, stat.st_mtime_ns)


@pytest.mark.parametrize("key", ["roster_sha256", "artifact_sha256", "artifact_id", "missing"])
def test_a_foreign_header_is_refused_by_its_keys(tmp_path, shard, clocked, capsys, key):
    artifact, roster, flags, whole = shard
    out = tmp_path / "rows.jsonl"
    lines = whole.decode().splitlines(keepends=True)
    head = json.loads(lines[0])
    if key == "missing":
        del head["families"]
        key = "families"
    else:
        head[key] = "another"
    data = (json.dumps(head) + "\n" + "".join(lines[1:3])).encode()
    out.write_bytes(data)
    stat = os.stat(out)
    capsys.readouterr()
    assert _resume(artifact, roster, out, flags) == 1
    text = capsys.readouterr().out
    refused = [line for line in text.splitlines() if line.startswith("REFUSED resume:")]
    assert len(refused) == 1 and key in refused[0], text
    assert "RESUME kept" not in text
    _untouched(out, data, stat)
    assert clocked["ran"] == []


def test_another_chunks_file_is_refused_by_its_keys(tmp_path, shard, clocked, capsys):
    artifact, roster, flags, whole = shard
    out = tmp_path / "rows.jsonl"
    other = [f if f != "0/2" else "1/2" for f in flags]
    _stopped(clocked, artifact, roster, out, other, 2)
    data, stat = out.read_bytes(), os.stat(out)
    capsys.readouterr()
    assert _resume(artifact, roster, out, flags) == 1
    refused = [line for line in capsys.readouterr().out.splitlines()
               if line.startswith("REFUSED resume:")]
    assert len(refused) == 1
    for key in ("chunk", "assigned", "partition"):
        assert key in refused[0]
    _untouched(out, data, stat)


@pytest.mark.parametrize("first", [b"", b"not json\n", b'{"kind": "row"}\n', b'{"kind": "hea'])
def test_a_file_that_is_not_rows_is_refused(tmp_path, shard, clocked, capsys, first):
    artifact, roster, flags, whole = shard
    out = tmp_path / "rows.jsonl"
    data = first + b"tail\n" if first else b"\n" + whole
    out.write_bytes(data)
    stat = os.stat(out)
    capsys.readouterr()
    assert _resume(artifact, roster, out, flags) == 1
    assert "REFUSED resume:" in capsys.readouterr().out
    _untouched(out, data, stat)


def test_a_file_changed_while_resuming_is_refused(tmp_path, shard, clocked, capsys, monkeypatch):
    artifact, roster, flags, whole = shard
    out = tmp_path / "rows.jsonl"
    _stopped(clocked, artifact, roster, out, flags, 3)
    out.write_bytes(out.read_bytes() + b'{"cut')
    scan = resumer.resume_prefix

    def raced(*args, **kwargs):
        found = scan(*args, **kwargs)
        with open(out, "ab") as other:   # another job writes after the scan
            other.write(b"more")
        return found

    monkeypatch.setattr(resumer, "resume_prefix", raced)
    data = out.read_bytes() + b"more"
    capsys.readouterr()
    assert _resume(artifact, roster, out, flags) == 1
    assert "REFUSED resume:" in capsys.readouterr().out
    assert out.read_bytes() == data
    assert not (tmp_path / "rows.jsonl.resume-dropped").exists()
    assert sorted(p.name for p in tmp_path.iterdir() if p.name.startswith("rows.jsonl")) == \
        ["rows.jsonl"]


def test_without_resume_the_file_is_overwritten(tmp_path, shard, clocked, capsys):
    artifact, roster, flags, whole = shard
    out = tmp_path / "rows.jsonl"
    _stopped(clocked, artifact, roster, out, flags, 3)
    capsys.readouterr()
    assert emit.main(_args(artifact, roster, out, flags)) == 0
    assert "RESUME" not in capsys.readouterr().out
    assert _stripped(out.read_bytes()) == _stripped(whole)
    assert len(clocked["ran"]) == _picked(whole)


def test_an_error_row_kept_sets_the_exit_code(tmp_path, artifact, roster_and_doc, clocked,
                                              monkeypatch):
    """The exit code is the uninterrupted run's: an error row among the kept
    ones still fails the shard."""
    roster = roster_and_doc[0]
    flags = ["--chunk", "0/2"]
    run = emit.Emitter.run_record

    def first_fails(self, record):
        row = run(self, record)
        if clocked["calls"] == 1:
            row["status"] = "error"
        return row

    monkeypatch.setattr(emit.Emitter, "run_record", first_fails)
    out = tmp_path / "rows.jsonl"
    _stopped(clocked, artifact, roster, out, flags, 2)
    assert _rows(out.read_bytes())[0]["status"] == "error"
    assert _resume(artifact, roster, out, flags) == 1
    assert sum(r["status"] == "error" for r in _rows(out.read_bytes())) == 1


@needs_database
def test_real_emission_resumes_to_the_uninterrupted_rows(tmp_path, artifact, roster_and_doc,
                                                         monkeypatch, capsys):
    """Real reference systems, cases and emission: two records, stopped after one."""
    roster = roster_and_doc[0]
    flags = ["--chunk", "0/2", "--limit", "2"]
    whole = _reference(artifact, roster, tmp_path / "uninterrupted.jsonl", flags)
    assert len(_rows(whole)) == 2
    run, ran = emit.Emitter.run_record, []

    def stop_after_one(self, record):
        if ran:
            raise Killed()
        ran.append(record["event_id"])
        return run(self, record)

    out = tmp_path / "rows.jsonl"
    monkeypatch.setattr(emit.Emitter, "run_record", stop_after_one)
    with pytest.raises(Killed):
        emit.main(_args(artifact, roster, out, flags))
    monkeypatch.setattr(emit.Emitter, "run_record", run)
    capsys.readouterr()
    assert _resume(artifact, roster, out, flags) == 0
    assert "RESUME kept=1 of 2" in capsys.readouterr().out
    assert _stripped(out.read_bytes()) == _stripped(whole)


def _statements(source: str, name: str):
    """The statements of function ``name`` in ``source``, one per line as
    ``ast.unparse`` writes them, with ``emit.X`` read as ``X`` and the
    function named ``run_emit``: line breaks and the module qualifier aside,
    the code itself."""
    import ast

    class Unqualify(ast.NodeTransformer):
        def visit_Attribute(self, node):
            self.generic_visit(node)
            if isinstance(node.value, ast.Name) and node.value.id == "emit":
                return ast.copy_location(ast.Name(id=node.attr, ctx=node.ctx), node)
            return node

    tree = ast.parse(source)
    found = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name]
    assert len(found) == 1, f"{name} is defined {len(found)} times"
    function = Unqualify().visit(found[0])
    function.name = "run_emit"
    return ast.unparse(function).splitlines()


# what ``run_emit_resume`` replaces of ``run_emit``: the fresh write ...
REPLACED = [
    "    with open(out, 'w') as log:",
    "        log.write(json.dumps(header, default=str) + '\\n')",
    "        for refusal in refusals:",
    "            log.write(json.dumps(refusal) + '\\n')",
    "        for (k, record) in enumerate(picked):",
]
# ... by the resume branch and a write that starts at the first record without a row
RESUME_BRANCH = [
    "    (start, mode) = (0, 'w')",
    "    if getattr(args, 'resume', False):",
    "        fixed = [json.dumps(header, default=str) + '\\n'] + "
    "[json.dumps(refusal) + '\\n' for refusal in refusals]",
    "        try:",
    "            resume = resume_prefix(out, fixed, picked)",
    "            rewrite_prefix(out, resume)",
    "        except ResumeRefusal as error:",
    "            print('REFUSED resume:', error, flush=True)",
    "            return 1",
    "        print(f'RESUME kept={resume.kept} of {len(picked)}', flush=True)",
    "        (start, mode, errors) = (resume.kept, 'a', resume.errors)",
    "    with open(out, mode) as log:",
    "        if mode == 'w':",
    "            log.write(json.dumps(header, default=str) + '\\n')",
    "            for refusal in refusals:",
    "                log.write(json.dumps(refusal) + '\\n')",
    "        for (k, record) in enumerate(picked[start:], start):",
]


def test_run_emit_resume_is_this_trees_run_emit_plus_the_resume_branch():
    """The script copies ``emit.run_emit``: a change to it under ``rmgpy/``
    that the script does not take fails here."""
    ours = _statements(inspect.getsource(emit), "run_emit")
    theirs = _statements(SCRIPT.read_text(), "run_emit_resume")
    removed, added = [], []
    for op, a0, a1, b0, b1 in difflib.SequenceMatcher(None, ours, theirs,
                                                     autojunk=False).get_opcodes():
        if op != "equal":
            removed += ours[a0:a1]
            added += theirs[b0:b1]
    assert (removed, added) == (REPLACED, RESUME_BRANCH)
    # and the replaced lines are one contiguous block of run_emit
    start = ours.index(REPLACED[0])
    assert ours[start:start + len(REPLACED)] == REPLACED
