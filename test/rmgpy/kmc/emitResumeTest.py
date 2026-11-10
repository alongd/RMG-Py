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
    # nothing dropped or temporary is left; the reference-case sidecar (I-326)
    # is saved before the rows file is read, is this run's own, and stays
    assert sorted(p.name for p in tmp_path.iterdir() if p.name.startswith("rows.jsonl")) == \
        ["rows.jsonl", "rows.jsonl.cases.pickle"]


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


# --- the reference-case sidecar (I-326) ------------------------------------
#
# Before its first row a shard builds its reference cases (``Emitter.init``:
# ``reference_cases``, ``cases_s``), 13 h for one cluster shard.  With
# ``--resume`` the script saves the cases, their stats and the species
# registry that step filled to ``<out>.cases.pickle`` and a later ``--resume``
# of the same shard loads them instead, only when the sidecar's identity is
# this run's; anything else is rebuilt.

INIT = emit.Emitter.init
REFERENCE_CASES = emit.reference_cases


def _witness(n):
    """A real ``Witness`` with its memoised (mappingproxy) views populated."""
    from rmgpy.kmc.boundaries import Witness
    row = ("C", 0, 0, 0, 0, ())
    witness = Witness(reactants=((1, row),), products=((1, row),), participants=(),
                      id_base=n, label=f"w{n}")
    witness.reactant_rows(), witness.product_rows(), witness.reactant_bonds()
    assert "_reactant_rows" in vars(witness)
    return witness


def _species_state(species):
    return dict(species._by_sig), {k: m.to_smiles() for k, m in species._first.items()}


@pytest.fixture
def built(monkeypatch, clocked):
    """A counting ``reference_cases`` that returns cases holding a real
    ``Witness`` and fills the species registry; and, per emitted record, the
    cases and species state the emitter read."""
    from rmgpy.molecule import Molecule
    state = {"builds": 0, "made": [], "seen": [], "unpicklable": False}

    def build(certifier, fam, systems, species, wanted=None):
        state["builds"] += 1
        species._by_sig["sig-CC"] = "CC"
        species._first["CC"] = Molecule(smiles="CC")
        cases = {("fam", ("CC",), "tmpl", ("C=C",)): [(_witness(1), {"n": 1})],
                 ("fam", ("CC",), "tmpl2", ("C",)): [(_witness(2), {"n": 2}),
                                                     (_witness(3), {"n": 3})]}
        if state["unpicklable"]:
            cases["bad"] = [(lambda: None, {})]
        state["made"].append((cases, _species_state(species)))
        return cases, {"fam|forward": 3}

    monkeypatch.setattr(emit, "reference_cases", build)
    run = emit.Emitter.run_record

    def spy(self, record):
        state["seen"].append((self.state["cases"], _species_state(self.state["species"]),
                              dict(self.state["case_stats"])))
        return run(self, record)

    monkeypatch.setattr(emit.Emitter, "run_record", spy)
    return state


def _sidecar(out):
    return out.with_name(out.name + ".cases.pickle")


def _stopped_resume(clocked, artifact, roster, out, flags, k):
    clocked["stop_after"], clocked["ran"] = k, []
    with pytest.raises(Killed):
        _resume(artifact, roster, out, flags)
    clocked["stop_after"], clocked["ran"] = None, []
    assert emit.Emitter.init is INIT


def _cases_lines(text):
    return [line for line in text.splitlines() if line.startswith("CASES")]


def test_a_valid_sidecar_skips_the_case_build(tmp_path, shard, clocked, built, capsys):
    artifact, roster, flags, whole = shard
    out = tmp_path / "rows.jsonl"
    _stopped_resume(clocked, artifact, roster, out, flags, 2)
    assert built["builds"] == 1 and _sidecar(out).is_file()
    text = capsys.readouterr().out
    assert any(line.startswith("CASES saved") for line in _cases_lines(text)), text
    made_cases, made_species = built["made"][0]
    built["seen"].clear()
    assert _resume(artifact, roster, out, flags) == 0
    text = capsys.readouterr().out
    assert built["builds"] == 1, "the resumed run rebuilt the reference cases"
    assert _cases_lines(text) == [f"CASES loaded {_sidecar(out)}"], text
    assert built["seen"], "no record was emitted"
    for cases, species, stats in built["seen"]:
        assert cases == made_cases and species == made_species
        assert stats["fam|forward"] == 3
    witness = built["seen"][0][0][("fam", ("CC",), "tmpl", ("C=C",))][0][0]
    assert dict(witness.reactant_rows()) == dict(_witness(1).reactant_rows())
    assert _stripped(out.read_bytes()) == _stripped(whole)
    assert emit.Emitter.init is INIT


def _edit_identity(path, key, value):
    data = path.read_bytes()
    end = data.index(b"\n") + 1
    head = json.loads(data[:end])
    assert key in head
    head[key] = value
    path.write_bytes((json.dumps(head) + "\n").encode() + data[end:])


@pytest.mark.parametrize("foreign", ["chunk", "artifact_sha256", "code_sources_sha256", "format"])
def test_a_foreign_sidecar_is_rebuilt(tmp_path, shard, clocked, built, capsys, foreign):
    artifact, roster, flags, whole = shard
    out = tmp_path / "rows.jsonl"
    if foreign == "chunk":
        other = tmp_path / "other.jsonl"
        assert _resume(artifact, roster, other,
                       [f if f != "0/2" else "1/2" for f in flags]) == 0
        _sidecar(other).replace(_sidecar(out))
    else:
        _stopped_resume(clocked, artifact, roster, out, flags, 1)
        _edit_identity(_sidecar(out), foreign, "another")
    old = _sidecar(out).read_bytes()
    builds = built["builds"]
    capsys.readouterr()
    assert _resume(artifact, roster, out, flags) == 0
    lines = _cases_lines(capsys.readouterr().out)
    assert built["builds"] == builds + 1
    assert lines[0].startswith("CASES rebuilt:") and foreign in lines[0], lines
    assert lines[1] == f"CASES saved {_sidecar(out)}"
    assert _stripped(out.read_bytes()) == _stripped(whole)
    assert out.with_name(out.name + ".cases.pickle.stale").read_bytes() == old
    # the rebuilt sidecar is this run's: the next resume loads it
    out.write_bytes(out.read_bytes()[:out.read_bytes().index(b"\n") + 1])
    assert _resume(artifact, roster, out, flags) == 0
    assert built["builds"] == builds + 1
    assert _cases_lines(capsys.readouterr().out) == [f"CASES loaded {_sidecar(out)}"]


@pytest.mark.parametrize("damage", ["empty", "cut_identity", "truncated", "flipped", "garbage"])
def test_a_corrupt_sidecar_is_rebuilt(tmp_path, shard, clocked, built, capsys, damage):
    artifact, roster, flags, whole = shard
    out = tmp_path / "rows.jsonl"
    _stopped_resume(clocked, artifact, roster, out, flags, 1)
    path = _sidecar(out)
    data = path.read_bytes()
    end = data.index(b"\n") + 1
    path.write_bytes({"empty": b"",
                      "cut_identity": data[:end // 2],
                      "truncated": data[:end + (len(data) - end) // 2],
                      "flipped": data[:-5] + bytes([data[-5] ^ 0x01]) + data[-4:],
                      "garbage": b"not json\n" + data[end:]}[damage])
    builds = built["builds"]
    capsys.readouterr()
    assert _resume(artifact, roster, out, flags) == 0
    lines = _cases_lines(capsys.readouterr().out)
    assert built["builds"] == builds + 1
    assert lines[0].startswith("CASES rebuilt:") and lines[1].startswith("CASES saved"), lines
    assert _stripped(out.read_bytes()) == _stripped(whole)


def test_without_resume_no_sidecar_is_read(tmp_path, shard, clocked, built, capsys):
    artifact, roster, flags, whole = shard
    out = tmp_path / "rows.jsonl"
    _stopped_resume(clocked, artifact, roster, out, flags, 1)
    path = _sidecar(out)
    data, before = path.read_bytes(), os.stat(path)
    builds = built["builds"]
    capsys.readouterr()
    assert resume_main(_args(artifact, roster, out, flags)) == 0
    assert built["builds"] == builds + 1
    assert _cases_lines(capsys.readouterr().out) == []
    _untouched(path, data, before)
    assert _stripped(out.read_bytes()) == _stripped(whole)
    assert emit.Emitter.init is INIT


def test_an_unsaveable_build_still_writes_rows(tmp_path, shard, clocked, built, capsys):
    artifact, roster, flags, whole = shard
    out = tmp_path / "rows.jsonl"
    built["unpicklable"] = True
    capsys.readouterr()
    assert _resume(artifact, roster, out, flags) == 0
    lines = _cases_lines(capsys.readouterr().out)
    assert lines[0].startswith("CASES rebuilt:") and lines[1].startswith("CASES not saved:"), lines
    assert _stripped(out.read_bytes()) == _stripped(whole)
    assert sorted(p.name for p in tmp_path.iterdir() if p.name.startswith("rows.jsonl")) == \
        ["rows.jsonl"]


def test_an_unknown_identity_is_never_trusted(tmp_path, artifact, roster_and_doc, clocked, built,
                                              capsys, monkeypatch):
    """A run whose database digest cannot be read cannot tell its sidecar from
    another run's: it neither saves nor loads one.  (The whole-artifact path:
    a prepared directory refuses an unknown database digest before any case.)"""
    roster = roster_and_doc[0]
    flags = ["--chunk", "0/2"]
    monkeypatch.setattr(emit, "database_content_sha256", lambda path: None)
    whole = _reference(artifact, roster, tmp_path / "uninterrupted.jsonl", flags)
    built["builds"] = 0
    out = tmp_path / "rows.jsonl"
    _stopped_resume(clocked, artifact, roster, out, flags, 1)
    assert not _sidecar(out).exists()
    capsys.readouterr()
    assert _resume(artifact, roster, out, flags) == 0
    lines = _cases_lines(capsys.readouterr().out)
    assert len(lines) == 1 and lines[0].startswith("CASES not cached:") and \
        "database_content_sha256" in lines[0], lines
    assert built["builds"] == 2
    assert not _sidecar(out).exists()
    assert _stripped(out.read_bytes()) == _stripped(whole)


def test_the_sidecar_identity_is_the_rows_header(tmp_path, shard, clocked, built):
    """The identity reads the header's own functions: every key they share
    holds the header's value (a copy that drifted would fail here)."""
    artifact, roster, flags, whole = shard
    out = tmp_path / "rows.jsonl"
    assert _resume(artifact, roster, out, flags) == 0
    data = _sidecar(out).read_bytes()
    identity = json.loads(data[:data.index(b"\n")])
    header = json.loads(out.read_bytes().split(b"\n")[0])
    shared = set(identity) & set(header)
    assert shared >= {"artifact_id", "artifact_sha256", "artifact_records", "families",
                      "database_path", "database_head", "database_content_sha256", "code_sha",
                      "code_sources_sha256", "code_extensions_sha256", "roster_sha256",
                      "systems", "chunk"}
    assert {k: identity[k] for k in shared} == {k: header[k] for k in shared}
    prepared = "--prepared" in flags
    assert (identity["prepared_manifest_sha256"] is not None) == prepared
    assert (identity["wanted_sha256"] is not None) == prepared


def test_a_witness_round_trips_without_its_views(tmp_path):
    import dataclasses
    from rmgpy.kmc.boundaries import Witness
    witness = _witness(7)
    path = tmp_path / "x.cases.pickle"
    identity = {"format": "test"}
    resumer.write_cases(path, identity, {"cases": {"k": [(witness, {})]}, "stats": {},
                                         "species": {}})
    payload, why = resumer.read_cases(path, identity)
    assert why is None
    loaded = payload["cases"]["k"][0][0]
    assert type(loaded) is Witness and loaded == witness and loaded is not witness
    assert set(vars(loaded)) == {f.name for f in dataclasses.fields(Witness)}
    assert dict(loaded.reactant_rows()) == dict(witness.reactant_rows())
    assert dict(loaded.reactant_bonds()) == dict(witness.reactant_bonds())
    assert sorted(p.name for p in tmp_path.iterdir()) == ["x.cases.pickle"]


@needs_database
def test_real_emission_resumes_from_the_saved_cases(tmp_path, artifact, roster_and_doc,
                                                    monkeypatch, capsys):
    """Real reference systems, cases and emission: a resumed run stopped after
    one record saves its cases; the next resume loads them, builds nothing,
    and writes the uninterrupted rows."""
    roster = roster_and_doc[0]
    flags = ["--chunk", "0/2", "--limit", "2"]
    whole = _reference(artifact, roster, tmp_path / "uninterrupted.jsonl", flags)
    builds, run, ran = [], emit.Emitter.run_record, []

    def counting(*args, **kwargs):
        builds.append(1)
        return REFERENCE_CASES(*args, **kwargs)

    def stop_after_one(self, record):
        if ran:
            raise Killed()
        ran.append(record["event_id"])
        return run(self, record)

    monkeypatch.setattr(emit, "reference_cases", counting)
    out = tmp_path / "rows.jsonl"
    monkeypatch.setattr(emit.Emitter, "run_record", stop_after_one)
    with pytest.raises(Killed):
        _resume(artifact, roster, out, flags)
    assert builds == [1] and _sidecar(out).is_file()
    monkeypatch.setattr(emit.Emitter, "run_record", run)
    capsys.readouterr()
    assert _resume(artifact, roster, out, flags) == 0
    text = capsys.readouterr().out
    assert builds == [1], "the resumed run rebuilt the reference cases"
    assert f"CASES loaded {_sidecar(out)}" in text and "RESUME kept=1 of 2" in text
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
