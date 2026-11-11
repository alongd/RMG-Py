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


def _statements(source: str, name: str, rename: str = "run_emit", docstring: bool = True):
    """The statements of function ``name`` in ``source``, one per line as
    ``ast.unparse`` writes them, with ``emit.X`` read as ``X`` and the
    function named ``rename`` (and without its docstring unless
    ``docstring``): line breaks and the module qualifier aside, the code
    itself."""
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
    function.name = rename
    if not docstring and ast.get_docstring(function) is not None:
        function.body = function.body[1:]
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


# ---------------------------------------------------------------------------
# the case build's checkpoint (I-327)
# ---------------------------------------------------------------------------
# A shard's case build can outlast a job: with ``--resume`` the script builds
# the cases through ``checkpointed_reference_cases`` (``emit.reference_cases``
# plus the checkpoint lines), saving the build's progress to
# ``<out>.cases.checkpoint``; a later ``--resume`` continues from it and ends
# with the uninterrupted build's cases, stats and species registry.  Here
# ``reference_cases`` is this tree's own, on eight synthetic subsets (cases
# with real ``Witness``es and real species keys), killed inside the build.

import hashlib  # noqa: E402
import re  # noqa: E402
import subprocess  # noqa: E402
from types import SimpleNamespace  # noqa: E402

from rmgpy.kmc import windows as win  # noqa: E402

SUBSET_SMILES = ("CC", "C=C", "[CH3]", "C")
SUBSETS = [("famA" if i < 5 else "famB", "forward" if i % 2 == 0 else "reverse", i)
           for i in range(8)]
SKIPPED = 3     # the subset a prepared shard's case filter does not admit


def _checkpoint(out):
    return out.with_name(out.name + ".cases.checkpoint")


def _case_view(cases, stats, species):
    """What a build returns, as plain comparable values: per key in order, its
    members in order (every ``MemberCase`` field, the witness, the wanted
    labels, the number of moves, the specs); the stats in order; the species
    registry (``species``: its ``__dict__``)."""
    members = [(key, [((c.label, c.family, c.direction, c.pathway, c.support, c.species,
                        c.kind), c.witness, sorted(c.wanted, key=repr), len(c.moves), specs)
                      for c, specs in value]) for key, value in cases.items()]
    first = [(k, [f.to_smiles() for f in (v if isinstance(v, (list, tuple)) else [v])])
             for k, v in species["_first"].items()]
    return members, list(stats.items()), list(species["_by_sig"].items()), first


@pytest.fixture
def subsetted(monkeypatch, clocked):
    """This tree's ``reference_cases`` on ``SUBSETS``: subset i has i % 3
    cases; on a prepared shard the case filter skips subset ``SKIPPED`` and
    drops the keys of famB's pathway t1.  ``kill_at``: the subset whose
    ``member_cases`` kills the job; ``tick``: build seconds per subset on the
    checkpoint's clock.  Each row carries the digest of the cases, stats and
    species emission read (``seen``: their views)."""
    from rmgpy.kmc.channels import SpeciesKeys
    from rmgpy.molecule import Molecule
    state = {"calls": 0, "kill_at": None, "tick": 0.0, "now": 0.0, "unpicklable": False,
             "seen": []}

    def participant(i):
        return SimpleNamespace(spec={"participant": i}, build=lambda label: ("build", i, label))

    def roster_subsets(fam, systems):
        for name, direction, i in SUBSETS:
            yield name, name, direction, [participant(i)], (i,)

    def side(n, direction):
        return [Molecule(smiles=SUBSET_SMILES[(n + (direction == "reverse")) % 4])]

    def member_cases(certifier, family, builds, direction, label=""):
        i = builds[0][1]
        if state["kill_at"] is not None and i >= state["kill_at"]:
            raise Killed()
        state["calls"] += 1
        state["now"] += state["tick"]
        return [win.MemberCase(label=f"{label}#{j}", family=family, direction=direction,
                               pathway=(f"t{(10 * i + j) % 2}",), witness=_witness(10 * i + j),
                               moves=(lambda: None,) if state["unpicklable"] else
                               (("move", i, j),), wanted=frozenset({f"w{j}"}))
                for j in range(i % 3)]

    def case_filter(certifier, records, species):
        keys = SpeciesKeys()
        wanted = set()
        for name, direction, i in SUBSETS:
            for j in range(i % 3):
                n = 10 * i + j
                other = "reverse" if direction == "forward" else "forward"
                key = (name, tuple(keys(m) for m in side(n, direction)), f"t{n % 2}",
                       tuple(keys(m) for m in side(n, other)))
                if not (name == "famB" and n % 2):
                    wanted.add(key)
        sides = {(name, d): frozenset(f"F{i}" for _n, _d, i in SUBSETS if i != SKIPPED)
                 for name in ("famA", "famB") for d in ("forward", "reverse")}
        return emit.CaseFilter(frozenset(wanted), sides)

    monkeypatch.setattr(emit, "reference_cases", REFERENCE_CASES)
    monkeypatch.setattr(emit, "roster_subsets", roster_subsets)
    monkeypatch.setattr(emit, "subset_formulas",
                        lambda chosen, memo: f"F{chosen[0].spec['participant']}")
    monkeypatch.setattr(emit, "case_filter", case_filter)
    monkeypatch.setattr(win, "member_cases", member_cases)
    monkeypatch.setattr(win, "template_molecules", lambda witness, direction:
                        side(witness.id_base, direction))
    if hasattr(resumer, "CaseCheckpoint"):   # absent: the build has no checkpoint at all
        monkeypatch.setattr(resumer.CaseCheckpoint, "clock", staticmethod(lambda: state["now"]))
    run = emit.Emitter.run_record

    def spy(self, record):
        S = self.state
        view = _case_view(S["cases"], S["case_stats"], S["species"].__dict__)
        state["seen"].append(view)
        row = run(self, record)
        row["cases_sha256"] = hashlib.sha256(repr(view).encode()).hexdigest()
        return row

    monkeypatch.setattr(emit.Emitter, "run_record", spy)
    return state


def _enumerated(flags):
    """The subsets that reach ``member_cases`` on this shard."""
    prepared = "--prepared" in flags
    return [i for _n, _d, i in SUBSETS if not (prepared and i == SKIPPED)]


def _fresh(tmp_path, shard, subsetted, capsys):
    """The uninterrupted resumed run's rows, the view emission read and its
    ``member_cases`` calls."""
    artifact, roster, flags, whole = shard
    out = tmp_path / "fresh.jsonl"
    subsetted["calls"], subsetted["seen"] = 0, []
    assert _resume(artifact, roster, out, flags) == 0
    capsys.readouterr()
    assert subsetted["seen"], "no record was emitted"
    seen = subsetted["seen"][0]
    assert all(view == seen for view in subsetted["seen"])
    assert subsetted["calls"] == len(_enumerated(flags))
    assert seen[0], "the build kept no case"
    subsetted["calls"], subsetted["seen"] = 0, []
    return out.read_bytes(), seen


def _killed_resume(subsetted, artifact, roster, out, flags, at):
    subsetted["kill_at"] = at
    with pytest.raises(Killed):
        _resume(artifact, roster, out, flags)
    subsetted["kill_at"] = None
    assert emit.Emitter.init is INIT


def _finished(subsetted, shard, out, fresh, seen):
    """The resumed run finishes with the uninterrupted cases and rows, and
    leaves the sidecar and no checkpoint."""
    artifact, roster, flags, whole = shard
    subsetted["seen"] = []
    assert _resume(artifact, roster, out, flags) == 0
    assert subsetted["seen"] and all(view == seen for view in subsetted["seen"])
    assert _stripped(out.read_bytes()) == _stripped(fresh)
    assert [{k: v for k, v in row.items() if k != "cases_sha256"}
            for row in _stripped(out.read_bytes())] == _stripped(whole)
    assert _sidecar(out).is_file() and not _checkpoint(out).exists()
    assert emit.Emitter.init is INIT


def _lines(text, prefix="CASES"):
    return [line for line in text.splitlines() if line.startswith(prefix)]


@pytest.mark.parametrize("k", [0, 1, 4])
def test_a_stopped_case_build_resumes_to_the_uninterrupted_cases(tmp_path, shard, subsetted,
                                                                  capsys, k):
    artifact, roster, flags, whole = shard
    flags = [*flags, "--case-checkpoint-seconds", "0"]
    shard = artifact, roster, flags, whole
    fresh, seen = _fresh(tmp_path, shard, subsetted, capsys)
    out = tmp_path / "rows.jsonl"
    _killed_resume(subsetted, artifact, roster, out, flags, k)
    done = [i for i in _enumerated(flags) if i < k]
    assert subsetted["calls"] == len(done)
    assert not _sidecar(out).exists()
    assert _checkpoint(out).exists() == (k > 0)
    text = capsys.readouterr().out
    assert _lines(text)[:2] == ["CASES rebuilt: no sidecar", "CASES build start 0/8 subsets"]
    assert [line.split(",")[0] for line in _lines(text, "CASES checkpoint ")] == \
        [f"CASES checkpoint {i}/8 subsets" for i in range(1, k + 1)]
    subsetted["calls"] = 0
    _finished(subsetted, shard, out, fresh, seen)
    assert subsetted["calls"] == len(_enumerated(flags)) - len(done), "a subset was built twice"
    lines = _lines(capsys.readouterr().out)
    if k:
        assert lines[1:3] == [f"CASES checkpoint resumed at {k}/8 subsets",
                              f"CASES build start {k}/8 subsets"], lines
    else:
        assert lines[:2] == ["CASES rebuilt: no sidecar", "CASES build start 0/8 subsets"], lines


def test_a_twice_stopped_build_resumes(tmp_path, shard, subsetted, capsys):
    artifact, roster, flags, whole = shard
    flags = [*flags, "--case-checkpoint-seconds", "0"]
    shard = artifact, roster, flags, whole
    fresh, seen = _fresh(tmp_path, shard, subsetted, capsys)
    out = tmp_path / "rows.jsonl"
    _killed_resume(subsetted, artifact, roster, out, flags, 2)
    _killed_resume(subsetted, artifact, roster, out, flags, 5)
    assert "CASES checkpoint resumed at 2/8 subsets" in _lines(capsys.readouterr().out)
    _finished(subsetted, shard, out, fresh, seen)
    assert subsetted["calls"] == len(_enumerated(flags))
    assert "CASES checkpoint resumed at 5/8 subsets" in _lines(capsys.readouterr().out)


def _rewrite_payload(path, change):
    data = path.read_bytes()
    head = json.loads(data[:data.index(b"\n")])
    identity = {k: v for k, v in head.items() if k not in resumer.PAYLOAD_KEYS}
    payload, why = resumer.read_cases(path, identity)
    assert why is None, why
    change(payload)
    resumer.write_cases(path, identity, payload, stale=False)


def _restarts(tmp_path, shard, subsetted, capsys, damage):
    """A checkpoint left at subset 4 and damaged by ``damage(path)`` is
    ignored (moved to ``.stale``): the build starts from subset 0 and the rows
    are the uninterrupted ones.  Returns the ignored line."""
    artifact, roster, flags, whole = shard
    flags = [*flags, "--case-checkpoint-seconds", "0"]
    shard = artifact, roster, flags, whole
    fresh, seen = _fresh(tmp_path, shard, subsetted, capsys)
    out = tmp_path / "rows.jsonl"
    _killed_resume(subsetted, artifact, roster, out, flags, 4)
    path = _checkpoint(out)
    damage(path)
    old = path.read_bytes()
    capsys.readouterr()
    subsetted["calls"] = 0
    _finished(subsetted, shard, out, fresh, seen)
    assert subsetted["calls"] == len(_enumerated(flags)), "the build did not restart"
    lines = _lines(capsys.readouterr().out)
    assert lines[1].startswith("CASES checkpoint ignored: ") and \
        lines[2] == "CASES build start 0/8 subsets", lines
    assert path.with_name(path.name + ".stale").read_bytes() == old
    return lines[1]


@pytest.mark.parametrize("foreign", ["chunk", "artifact_sha256", "code_sources_sha256",
                                     "subsets"])
def test_a_foreign_checkpoint_restarts_the_build(tmp_path, shard, subsetted, capsys, foreign):
    def damage(path):
        if foreign == "subsets":
            _rewrite_payload(path, lambda payload: payload.update(subsets=9))
        else:
            _edit_identity(path, foreign, "another")

    assert foreign in _restarts(tmp_path, shard, subsetted, capsys, damage)


@pytest.mark.parametrize("damage", ["empty", "cut_identity", "truncated", "flipped", "garbage",
                                    "no_next", "next_past_the_end"])
def test_a_corrupt_checkpoint_restarts_the_build(tmp_path, shard, subsetted, capsys, damage):
    def corrupt(path):
        if damage == "no_next":
            return _rewrite_payload(path, lambda payload: payload.pop("next"))
        if damage == "next_past_the_end":
            return _rewrite_payload(path, lambda payload: payload.update(next=9))
        data = path.read_bytes()
        end = data.index(b"\n") + 1
        path.write_bytes({"empty": b"",
                          "cut_identity": data[:end // 2],
                          "truncated": data[:end + (len(data) - end) // 2],
                          "flipped": data[:-5] + bytes([data[-5] ^ 0x01]) + data[-4:],
                          "garbage": b"not json\n" + data[end:]}[damage])

    _restarts(tmp_path, shard, subsetted, capsys, corrupt)


def test_without_resume_no_checkpoint_is_read_or_written(tmp_path, shard, subsetted, capsys):
    artifact, roster, flags, whole = shard
    out = tmp_path / "rows.jsonl"
    _killed_resume(subsetted, artifact, roster, out, [*flags, "--case-checkpoint-seconds", "0"], 4)
    path = _checkpoint(out)
    data, before = path.read_bytes(), os.stat(path)
    names = sorted(p.name for p in tmp_path.iterdir() if p.name.startswith("rows.jsonl."))
    capsys.readouterr()
    subsetted["calls"] = 0
    assert resume_main(_args(artifact, roster, out, [*flags, "--case-checkpoint-seconds", "0"])) \
        == 0
    assert subsetted["calls"] == len(_enumerated(flags))
    assert _lines(capsys.readouterr().out) == []
    _untouched(path, data, before)
    assert sorted(p.name for p in tmp_path.iterdir() if p.name.startswith("rows.jsonl.")) == names
    assert _stripped(out.read_bytes()) == _stripped(_fresh(tmp_path, shard, subsetted, capsys)[0])


def test_checkpoints_follow_the_interval(tmp_path, artifact, roster_and_doc, clocked, subsetted,
                                         capsys):
    """At the default 600 s and 250 s per subset, the build saves before
    subsets 3 and 6, one progress line each; a job killed at subset 7 leaves
    the checkpoint of subset 6."""
    roster = roster_and_doc[0]
    flags = ["--chunk", "0/2"]
    out = tmp_path / "rows.jsonl"
    subsetted["tick"] = 250.0
    _killed_resume(subsetted, artifact, roster, out, flags, 7)
    lines = _lines(capsys.readouterr().out, "CASES checkpoint ")
    assert [re.sub(r"\d+ bytes", "N bytes", line) for line in lines] == [
        "CASES checkpoint 3/8 subsets, build 750s, N bytes, write 0.00s",
        "CASES checkpoint 6/8 subsets, build 1500s, N bytes, write 0.00s"], lines
    assert int(re.search(r"(\d+) bytes", lines[-1]).group(1)) == _checkpoint(out).stat().st_size
    data = _checkpoint(out).read_bytes()
    head = json.loads(data[:data.index(b"\n")])
    payload, why = resumer.read_cases(_checkpoint(out), {
        k: v for k, v in head.items() if k not in resumer.PAYLOAD_KEYS})
    assert why is None and (payload["next"], payload["subsets"]) == (6, 8)
    assert list(payload["stats"].items()) == [("famA|forward", 0 + 2 + 1), ("famA|reverse", 1 + 0),
                                              ("famB|reverse", 2)]


@pytest.mark.parametrize("value", ["-1", "nan", "inf", "x"])
def test_the_case_checkpoint_option_refuses_bad_values(tmp_path, shard, subsetted, capsys, value):
    artifact, roster, flags, whole = shard
    out = tmp_path / "rows.jsonl"
    with pytest.raises(SystemExit) as stop:
        _resume(artifact, roster, out, [*flags, "--case-checkpoint-seconds", value])
    assert stop.value.code == 2
    # the option's own refusal, not argparse's "unrecognized arguments"
    assert f"argument --case-checkpoint-seconds: '{value}' is not" in capsys.readouterr().err
    assert not out.exists() and not _checkpoint(out).exists()


def test_an_unsaveable_checkpoint_does_not_stop_the_build(tmp_path, shard, subsetted, capsys):
    """Subset 0 has no case, so the save before subset 1 succeeds; every later
    save holds an unpicklable case and fails, and so does the sidecar, which
    leaves the checkpoint of subset 1 in place."""
    artifact, roster, flags, whole = shard
    flags = [*flags, "--case-checkpoint-seconds", "0"]
    out = tmp_path / "rows.jsonl"
    subsetted["unpicklable"] = True
    capsys.readouterr()
    assert _resume(artifact, roster, out, flags) == 0
    lines = _lines(capsys.readouterr().out)
    assert [re.sub(r"\d+ bytes", "N bytes", line.split(":")[0]) for line in lines
            if line.startswith("CASES checkpoint")] == \
        ["CASES checkpoint 1/8 subsets, build 0s, N bytes, write 0.00s"] + \
        [f"CASES checkpoint not saved at {i}/8" for i in range(2, 8)], lines
    assert lines[-1].startswith("CASES not saved:"), lines
    assert subsetted["calls"] == len(_enumerated(flags))
    assert [{k: v for k, v in row.items() if k != "cases_sha256"}
            for row in _stripped(out.read_bytes())] == _stripped(whole)
    assert sorted(p.name for p in tmp_path.iterdir() if p.name.startswith("rows.jsonl")) == \
        ["rows.jsonl", "rows.jsonl.cases.checkpoint"]


# what ``checkpointed_reference_cases`` changes of ``reference_cases``: its
# signature and the subset loop's head ...
CASES_REPLACED = [
    "def reference_cases(certifier, fam, systems, species, wanted: Optional[CaseFilter]=None):",
    "    for (name, family, direction, chosen, subset) in roster_subsets(fam, systems):",
]
# ... by a checkpoint argument, its restore, and a loop that starts at the
# restored subset and saves before each subset
CASES_CHECKPOINT = [
    "def reference_cases(certifier, fam, systems, species, wanted=None, *, checkpoint):",
    "    start = checkpoint.restore(fam, systems, out, stats, species)",
    "    for (index, (name, family, direction, chosen, subset)) in "
    "enumerate(roster_subsets(fam, systems)):",
    "        if index < start:",
    "            continue",
    "        checkpoint.save(index, out, stats, species)",
]
# ... and (I-329) the one ``member_cases`` call, which reads the subset's
# result from a member-case table (``--member-cases``) when it has it
CASES_MEMBER = (
    "        cases = win.member_cases(certifier, family, builds, direction, "
    "label=f'{name}|{direction}|{subset}')",
    "        cases = checkpoint.member_cases(certifier, family, builds, direction, "
    "label=f'{name}|{direction}|{subset}', chosen=chosen)",
)


def test_the_copy_is_this_trees_reference_cases_plus_the_checkpoint_lines():
    """The script copies ``emit.reference_cases``: a change to it under
    ``rmgpy/`` that the copy does not take fails here."""
    ours = _statements(inspect.getsource(emit), "reference_cases", "reference_cases",
                       docstring=False)
    theirs = _statements(SCRIPT.read_text(), "checkpointed_reference_cases", "reference_cases",
                         docstring=False)
    removed, added = [], []
    for op, a0, a1, b0, b1 in difflib.SequenceMatcher(None, ours, theirs,
                                                     autojunk=False).get_opcodes():
        if op != "equal":
            removed += ours[a0:a1]
            added += theirs[b0:b1]
    assert (removed, added) == ([*CASES_REPLACED, CASES_MEMBER[0]],
                                [*CASES_CHECKPOINT, CASES_MEMBER[1]])
    # the loop body after the checkpoint lines is the original's, whole, but
    # for the one member_cases statement
    body = ours[ours.index(CASES_REPLACED[1]) + 1:]
    assert body.count(CASES_MEMBER[0]) == 1
    body[body.index(CASES_MEMBER[0])] = CASES_MEMBER[1]
    assert theirs[theirs.index(CASES_CHECKPOINT[-1]) + 1:] == body


# a job killed by the driver below exits with this code
KILLED = 77
DRIVER = """
import json, os, sys
sys.path.insert(0, sys.argv[1])
cfg = json.loads(sys.argv[2])
from rmgpy.kmc import windows as win
import emitResumeTest as T
if cfg["kill_at"] is not None:
    member_cases, done = win.member_cases, []

    def killing(*args, **kwargs):
        if len(done) >= cfg["kill_at"]:
            sys.stdout.flush()
            os._exit(T.KILLED)
        done.append(1)
        return member_cases(*args, **kwargs)

    win.member_cases = killing
sys.exit(T.resume_main(cfg["argv"]))
"""


def _job(tmp_path, argv, kill_at=None):
    """One cluster job: the script in a new process (``kill_at``: killed in
    the case build when that many subsets are done)."""
    driver = tmp_path / "driver.py"
    driver.write_text(DRIVER)
    root = str(SCRIPT.parents[1])
    env = {**os.environ, "PYTHONPATH": os.pathsep.join(
        [root, *filter(None, [os.environ.get("PYTHONPATH")])])}
    done = subprocess.run([sys.executable, str(driver), str(Path(__file__).resolve().parent),
                           json.dumps({"argv": argv, "kill_at": kill_at})],
                          cwd=tmp_path, env=env, capture_output=True, text=True, timeout=3600)
    return done.returncode, done.stdout, done.stderr


def _sidecar_view(path):
    data = path.read_bytes()
    head = json.loads(data[:data.index(b"\n")])
    payload, why = resumer.read_cases(path, {k: v for k, v in head.items()
                                             if k not in resumer.PAYLOAD_KEYS})
    assert why is None, why
    return _case_view(payload["cases"], payload["stats"], payload["species"])


@needs_database
def test_real_case_build_resumes_in_a_fresh_process(tmp_path, artifact, roster_and_doc, capsys):
    """Real reference systems, cases and emission, each job its own process: a
    build killed after 40 subsets (a checkpoint after every one) is continued
    by the next job from subset 40, which ends with the uninterrupted job's
    cases, stats, species registry and rows."""
    roster = roster_and_doc[0]
    flags = ["--chunk", "0/2", "--limit", "1"]
    whole = tmp_path / "uninterrupted.jsonl"
    code, text, err = _job(tmp_path, _args(artifact, roster, whole, [*flags, "--resume"]))
    assert code == 0, err[-4000:]
    total = int(re.search(r"CASES build start 0/(\d+) subsets", text).group(1))
    assert total > 40
    out = tmp_path / "rows.jsonl"
    code, text, err = _job(tmp_path, _args(artifact, roster, out, [
        *flags, "--resume", "--case-checkpoint-seconds", "0"]), kill_at=40)
    assert code == KILLED, err[-4000:]
    saves = re.findall(r"CASES checkpoint (\d+)/\d+ subsets, build \d+s, (\d+) bytes, "
                       r"write ([\d.]+)s", text)
    assert [int(i) for i, _b, _w in saves] == list(range(1, 41)), text
    assert _checkpoint(out).is_file() and not _sidecar(out).exists()
    with capsys.disabled():
        print(f"\nI-327 toy checkpoint at 40/{total} subsets: "
              f"{max(int(b) for _i, b, _w in saves)} bytes, "
              f"write peak {max(float(w) for _i, _b, w in saves):.2f}s")
    code, text, err = _job(tmp_path, _args(artifact, roster, out, [*flags, "--resume"]))
    assert code == 0, err[-4000:]
    assert f"CASES checkpoint resumed at 40/{total} subsets" in text, text
    assert f"CASES build start 40/{total} subsets" in text, text
    assert _sidecar_view(_sidecar(out)) == _sidecar_view(_sidecar(whole))
    assert not _checkpoint(out).exists()
    assert _stripped(out.read_bytes()) == _stripped(whole.read_bytes())


def test_a_checkpoint_that_cannot_be_removed_does_not_fail_the_job(tmp_path, monkeypatch,
                                                                   capsys):
    """Once the sidecar is saved, a checkpoint the job cannot remove is
    reported, not raised: the shard's cases are safe in the sidecar."""
    path = tmp_path / "rows.jsonl.cases.checkpoint"
    path.write_bytes(b"x")

    def refuse(self, missing_ok=False):
        raise PermissionError(13, "Permission denied", str(self))

    monkeypatch.setattr(Path, "unlink", refuse)
    resumer.CaseCheckpoint(path, {}, 0).finish()
    assert capsys.readouterr().out.startswith("CASES checkpoint not removed: ")
    monkeypatch.undo()
    resumer.CaseCheckpoint(path, {}, 0).finish()
    assert not path.exists()
    resumer.CaseCheckpoint(path, {}, 0).finish()     # already gone: nothing to say
    assert capsys.readouterr().out == ""


# ---------------------------------------------------------------------------
# the member-case table (I-329)
# ---------------------------------------------------------------------------
# Nearly all of a shard's case build is ``win.member_cases``, whose result
# depends on the family, the direction and the chosen participants, not on the
# shard.  ``scripts/emit_member_cases.py build --subsets k/K`` computes it for
# subsets i % K == k of ``emit.roster_subsets`` (K independent jobs), ``merge``
# joins the parts into one table, and ``emit_resume.py --resume --member-cases
# TABLE`` has its loop copy read each subset's result from the table instead of
# calling it.  The rest of the loop is the shard's own, so its cases, stats,
# species registry and rows are those of the in-shard build; a subset the
# table lacks, or a table of another identity, is computed in the shard.

import pickle  # noqa: E402
import time  # noqa: E402

BUILDER = SCRIPT.with_name("emit_member_cases.py")
TABLE_LINE = re.compile(r"CASES member-case table: (\d+) read, (\d+) computed of (\d+) enumerated")
# every field of the table identity: a table with another value of any is not used
TABLE_KEYS = ["format", "artifact_id", "artifact_sha256", "roster_sha256", "database_head",
              "database_content_sha256", "code_sha", "code_sources_sha256",
              "code_extensions_sha256", "systems", "subsets", "subsets_sha256"]


def _builder():
    """The build and merge script (it imports ``emit_resume``: the module above)."""
    spec = importlib.util.spec_from_file_location("emit_member_cases", BUILDER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _system_of(i):
    return "sysA" if i < 4 else "sysB"


@pytest.fixture
def tabled(monkeypatch, subsetted):
    """``subsetted`` with real systems: subset i is participant i % 4 of
    ``sysA`` (i < 4) or ``sysB``, so subsets 0 and 4 share the label
    ``famA|forward|(0,)`` and only the system tells them apart."""
    parts = [SimpleNamespace(spec={"participant": i}, label=f"part{i}",
                             build=lambda label, i=i: ("build", i, label)) for i in range(8)]
    systems = {"sysA": parts[:4], "sysB": parts[4:]}

    def roster_subsets(fam, systems):
        for name, direction, i in SUBSETS:
            yield name, name, direction, [systems[_system_of(i)][i % 4]], (i % 4,)

    monkeypatch.setattr(emit, "build_reference_systems", lambda roster: (systems, []))
    monkeypatch.setattr(emit, "roster_subsets", roster_subsets)
    subsetted["builder"] = _builder()
    return subsetted


def _build(tabled, artifact, roster, out, k, K, *flags, code=0):
    assert tabled["builder"].main(["build", str(artifact), str(roster), "--database",
                                   str(DATABASE), "--subsets", f"{k}/{K}", "--out", str(out),
                                   *flags]) == code
    return out


def _merge(tabled, parts, out, *flags):
    return tabled["builder"].main(["merge", *map(str, parts), "--out", str(out), *flags])


def _table(tabled, shard, tmp_path, K, parts=None, *flags):
    """The table of the ``parts`` (default all) of K, each built in turn."""
    artifact, roster, _flags, _whole = shard
    built = [_build(tabled, artifact, roster, tmp_path / f"part{k}of{K}.table", k, K)
             for k in (range(K) if parts is None else parts)]
    out = tmp_path / f"table{K}.table"
    assert _merge(tabled, built, out, *flags) == 0
    tabled["calls"] = 0
    return out


def _tabled_resume(tabled, shard, out, table, *flags):
    """A resumed run with ``table``: its rows, printed text and table line."""
    artifact, roster, shard_flags, _whole = shard
    tabled["calls"], tabled["seen"] = 0, []
    assert _resume(artifact, roster, out, [*shard_flags, *flags, "--member-cases",
                                           str(table)]) == 0
    return out.read_bytes()


def _same(tabled, rows, fresh, seen):
    assert tabled["seen"] and all(view == seen for view in tabled["seen"])
    assert _stripped(rows) == _stripped(fresh)


def _table_line(text):
    found = TABLE_LINE.findall(text)
    assert len(found) == 1, text
    return tuple(int(n) for n in found[0])


@pytest.mark.parametrize("K", [1, 2, 3])
def test_a_merged_table_is_the_in_shard_build(tmp_path, shard, tabled, capsys, K):
    fresh, seen = _fresh(tmp_path, shard, tabled, capsys)
    table = _table(tabled, shard, tmp_path, K)
    capsys.readouterr()
    rows = _tabled_resume(tabled, shard, tmp_path / "rows.jsonl", table)
    text = capsys.readouterr().out
    _same(tabled, rows, fresh, seen)
    n = len(_enumerated(shard[2]))
    assert _table_line(text) == (n, 0, n) and tabled["calls"] == 0, text


def test_without_member_cases_nothing_changes(tmp_path, shard, tabled, capsys):
    artifact, roster, flags, _whole = shard
    tabled["calls"], tabled["seen"] = 0, []
    assert _resume(artifact, roster, tmp_path / "rows.jsonl", flags) == 0
    text = capsys.readouterr().out
    assert "member-case table" not in text
    assert tabled["calls"] == len(_enumerated(flags))


def _edit_table_identity(path, key, value):
    data = path.read_bytes()
    cut = data.index(b"\n") + 1
    head = json.loads(data[:cut])
    head[key] = value
    path.write_bytes(json.dumps(head, sort_keys=True).encode() + b"\n" + data[cut:])


def _computes_all(tabled, shard, tmp_path, capsys, table, fresh, seen, why):
    rows = _tabled_resume(tabled, shard, tmp_path / "rows.jsonl", table)
    text = capsys.readouterr().out
    _same(tabled, rows, fresh, seen)
    n = len(_enumerated(shard[2]))
    assert tabled["calls"] == n
    assert any(line.startswith("CASES member-case table not used: ") and why in line
               for line in text.splitlines()), text
    assert _table_line(text) == (0, n, n), text


@pytest.mark.parametrize("key", TABLE_KEYS)
def test_a_foreign_table_is_not_used(tmp_path, shard, tabled, capsys, key):
    fresh, seen = _fresh(tmp_path, shard, tabled, capsys)
    table = _table(tabled, shard, tmp_path, 2)
    _edit_table_identity(table, key, "another")
    capsys.readouterr()
    _computes_all(tabled, shard, tmp_path, capsys, table, fresh, seen, key)


@pytest.mark.parametrize("damage", ["truncated", "flipped", "not a table", "empty", "missing"])
def test_a_damaged_table_is_not_used(tmp_path, shard, tabled, capsys, damage):
    fresh, seen = _fresh(tmp_path, shard, tabled, capsys)
    table = _table(tabled, shard, tmp_path, 2)
    data = table.read_bytes()
    if damage == "truncated":
        table.write_bytes(data[:-10])
    elif damage == "flipped":
        middle = (data.index(b"\n") + len(data)) // 2
        table.write_bytes(data[:middle] + bytes([data[middle] ^ 0xFF]) + data[middle + 1:])
    elif damage == "not a table":
        table.write_bytes(b'{"kind": "row"}\n')
    elif damage == "empty":
        table.write_bytes(b"")
    else:
        table.unlink()
    capsys.readouterr()
    _computes_all(tabled, shard, tmp_path, capsys, table, fresh, seen, "")


def test_an_unknown_identity_never_uses_a_table(tmp_path, shard, tabled, capsys, monkeypatch):
    fresh, seen = _fresh(tmp_path, shard, tabled, capsys)
    table = _table(tabled, shard, tmp_path, 2)
    identity = resumer.case_identity
    monkeypatch.setattr(resumer, "case_identity", lambda *a: {
        **identity(*a), "code_extensions_sha256": None})
    capsys.readouterr()
    rows = _tabled_resume(tabled, shard, tmp_path / "rows.jsonl", table)
    text = capsys.readouterr().out
    _same(tabled, rows, fresh, seen)
    assert tabled["calls"] == len(_enumerated(shard[2]))
    assert "CASES member-case table not used: unknown code_extensions_sha256" in text, text


def test_a_change_under_rmgpy_refuses_the_table(tmp_path, shard, tabled, capsys, monkeypatch):
    """The table identity carries ``code_sources_sha256``: a source change under
    ``rmgpy/`` (which could change ``roster_subsets`` or ``member_cases``)
    makes every table built before it unused."""
    fresh, seen = _fresh(tmp_path, shard, tabled, capsys)
    with monkeypatch.context() as before:
        before.setattr(emit, "code_sources_sha256", lambda: "c" * 64)
        table = _table(tabled, shard, tmp_path, 2)
    capsys.readouterr()
    _computes_all(tabled, shard, tmp_path, capsys, table, fresh, seen, "code_sources_sha256")


def _entries(path):
    data = path.read_bytes()
    cut = data.index(b"\n") + 1
    head = json.loads(data[:cut])
    body = data[cut:]
    index = pickle.loads(body[int(head["index_offset"]):])
    return head, cut, index


def test_a_bad_entry_is_computed(tmp_path, shard, tabled, capsys):
    """An entry whose bytes do not match its own digest (the file's digest
    made to match) is computed; every other entry is read."""
    fresh, seen = _fresh(tmp_path, shard, tabled, capsys)
    table = _table(tabled, shard, tmp_path, 2)
    head, cut, index = _entries(table)
    key = ("sysA", "famA|reverse|(1,)")      # subset 1: one case, enumerated on every path
    off, length, _sha = index["entries"][key]
    body = bytearray(table.read_bytes()[cut:])
    body[off + length // 2] ^= 0xFF
    head["payload_sha256"] = hashlib.sha256(bytes(body)).hexdigest()
    table.write_bytes(json.dumps(head, sort_keys=True).encode() + b"\n" + bytes(body))
    capsys.readouterr()
    rows = _tabled_resume(tabled, shard, tmp_path / "rows.jsonl", table)
    text = capsys.readouterr().out
    _same(tabled, rows, fresh, seen)
    n = len(_enumerated(shard[2]))
    assert tabled["calls"] == 1 and _table_line(text) == (n - 1, 1, n), text
    assert "CASES member-case table entry sysA famA|reverse|(1,) not used" in text, text


def test_a_partial_table_computes_the_missing_subsets(tmp_path, shard, tabled, capsys):
    fresh, seen = _fresh(tmp_path, shard, tabled, capsys)
    table = _table(tabled, shard, tmp_path, 3, [0, 1], "--partial")
    assert "TABLE missing 2 of 8 subsets" in capsys.readouterr().out
    rows = _tabled_resume(tabled, shard, tmp_path / "rows.jsonl", table)
    text = capsys.readouterr().out
    _same(tabled, rows, fresh, seen)
    enumerated = _enumerated(shard[2])
    missing = [i for i in enumerated if i % 3 == 2]
    assert tabled["calls"] == len(missing)
    assert _table_line(text) == (len(enumerated) - len(missing), len(missing), len(enumerated))


def test_a_table_build_stopped_in_the_shard_resumes(tmp_path, shard, tabled, capsys):
    """A shard killed in its case build with a partial table (it computes the
    subsets the table lacks) resumes from its checkpoint to the in-shard cases."""
    artifact, roster, flags, whole = shard
    fresh, seen = _fresh(tmp_path, shard, tabled, capsys)
    table = _table(tabled, shard, tmp_path, 3, [0, 1], "--partial")
    out = tmp_path / "rows.jsonl"
    extra = ["--case-checkpoint-seconds", "0", "--member-cases", str(table)]
    _killed_resume(tabled, artifact, roster, out, [*flags, *extra], 5)
    assert tabled["calls"] == 1 and _checkpoint(out).is_file()     # subset 2, then killed at 5
    capsys.readouterr()
    rows = _tabled_resume(tabled, shard, out, table, "--case-checkpoint-seconds", "0")
    text = capsys.readouterr().out
    _same(tabled, rows, fresh, seen)
    assert tabled["calls"] == 1, "a subset was computed twice"
    assert "CASES checkpoint resumed at 5/8 subsets" in text, text
    assert _sidecar(out).is_file() and not _checkpoint(out).exists()


def test_member_cases_needs_resume(tmp_path, shard, tabled, capsys):
    artifact, roster, flags, _whole = shard
    with pytest.raises(SystemExit):
        resume_main(_args(artifact, roster, tmp_path / "rows.jsonl",
                          [*flags, "--member-cases", str(tmp_path / "t.table")]))
    assert not (tmp_path / "rows.jsonl").exists()


@pytest.mark.parametrize("value", ["1/1", "2", "a/b", "-1/2", "0/0", "3/2"])
def test_the_subsets_option_refuses_bad_values(tmp_path, shard, tabled, value):
    artifact, roster, _flags, _whole = shard
    with pytest.raises(SystemExit):
        tabled["builder"].main(["build", str(artifact), str(roster), "--database",
                                str(DATABASE), "--subsets", value, "--out",
                                str(tmp_path / "p.table")])
    assert not (tmp_path / "p.table").exists()


def _member_lines(text):
    return [line for line in text.splitlines() if line.startswith("MEMBER ")]


@pytest.mark.parametrize("K", [1, 3])
def test_the_build_enumerates_roster_subsets(tmp_path, shard, tabled, capsys, K):
    """Part k of K builds subsets i % K == k of ``emit.roster_subsets`` (here
    the fixture's), in its order, one progress line each after its share."""
    artifact, roster, _flags, _whole = shard
    expected = [(i, _system_of(i), f"{name}|{direction}|({i % 4},)")
                for name, direction, i in SUBSETS]
    for k in range(K):
        capsys.readouterr()
        _build(tabled, artifact, roster, tmp_path / f"p{k}.table", k, K)
        lines = _member_lines(capsys.readouterr().out)
        mine = [e for e in expected if e[0] % K == k]
        assert lines[0] == f"MEMBER part {k}/{K}: {len(mine)} of 8 subsets", lines
        done = [re.match(rf"MEMBER (\d+)/8 part {k}/{K} [\d.]+s cases=\d+ (\S+) (.+)$", line)
                for line in lines[1:1 + len(mine)]]
        assert all(done), lines
        assert [(int(m.group(1)), m.group(2), m.group(3)) for m in done] == mine


def test_the_build_calls_member_cases_as_the_loop_does():
    """The build step's builds and ``member_cases`` call are the statements of
    ``emit.reference_cases`` itself: a change to them under ``rmgpy/`` fails
    here (and the label of a table entry is the shard's)."""
    ours = [line.strip() for line in _statements(inspect.getsource(emit), "reference_cases",
                                                 "reference_cases", docstring=False)]
    theirs = [line.strip() for line in _statements(BUILDER.read_text(), "build_part",
                                                   "build_part", docstring=False)]
    calls = [CASES_MEMBER[0].strip(),
             "builds = [p.build(f'p{k}') for (k, p) in enumerate(chosen)]"]
    for line in calls:
        assert ours.count(line) == 1 and theirs.count(line) == 1, line
    # and the build enumerates through the module attribute, never a copy
    assert "def roster_subsets" not in BUILDER.read_text()
    assert any("roster_subsets(fam, systems)" in line for line in theirs)


def _table_blobs(path):
    head, cut, index = _entries(path)
    data = path.read_bytes()[cut:]
    return {key: data[off:off + n] for key, (off, n, _sha) in index["entries"].items()}, \
        index["meta"]


def test_a_part_rerun_keeps_its_done_subsets(tmp_path, shard, tabled, capsys):
    artifact, roster, _flags, _whole = shard
    whole = _build(tabled, artifact, roster, tmp_path / "whole.table", 0, 1)
    tabled["calls"], tabled["kill_at"] = 0, 4
    out = tmp_path / "part.table"
    with pytest.raises(Killed):
        _build(tabled, artifact, roster, out, 0, 1, "--save-seconds", "0")
    tabled["kill_at"] = None
    blobs, meta = _table_blobs(out)
    assert len(blobs) == 4 and meta == {"part": [0, 1], "complete": False}
    tabled["calls"] = 0
    capsys.readouterr()
    _build(tabled, artifact, roster, out, 0, 1)
    assert tabled["calls"] == 4, "a done subset was built again"
    assert "MEMBER part resumed: 4 done" in capsys.readouterr().out
    assert _table_blobs(out) == _table_blobs(whole)


@pytest.mark.parametrize("case", ["identity", "duplicate", "K", "missing", "incomplete", "out"])
def test_merge_refuses(tmp_path, shard, tabled, capsys, case):
    artifact, roster, _flags, _whole = shard
    parts = [_build(tabled, artifact, roster, tmp_path / f"p{k}.table", k, 2) for k in (0, 1)]
    out = tmp_path / "table.table"
    partial_ok = False
    if case == "identity":
        _edit_table_identity(parts[1], "code_sha", "another")
        why = "code_sha"
    elif case == "duplicate":
        parts = [parts[0], parts[0]]
        why = "part 0/2"
    elif case == "K":
        parts[1] = _build(tabled, artifact, roster, tmp_path / "q.table", 1, 3)
        why = "K"
    elif case == "missing":
        parts, why, partial_ok = parts[:1], "missing parts: 1", True
    elif case == "incomplete":
        parts[1].unlink()
        tabled["kill_at"] = 3
        with pytest.raises(Killed):
            _build(tabled, artifact, roster, parts[1], 1, 2, "--save-seconds", "0")
        tabled["kill_at"] = None
        why, partial_ok = "incomplete", True
    else:
        out, why = parts[0], "is one of the parts"
    before = {p: p.read_bytes() for p in parts}
    capsys.readouterr()
    assert _merge(tabled, parts, out) == 1
    text = capsys.readouterr().out
    assert text.startswith("REFUSED") and why in text, text
    assert {p: p.read_bytes() for p in parts} == before
    if out not in parts:
        assert not out.exists()
    if partial_ok:
        assert _merge(tabled, parts, out, "--partial") == 0 and out.is_file()


def _script(cwd, argv):
    root = str(SCRIPT.parents[1])
    env = {**os.environ, "PYTHONPATH": os.pathsep.join(
        [root, *filter(None, [os.environ.get("PYTHONPATH")])])}
    done = subprocess.run([sys.executable, *map(str, argv)], cwd=cwd, env=env,
                          capture_output=True, text=True, timeout=3600)
    return done.returncode, done.stdout, done.stderr


@pytest.fixture(scope="module")
def real_tables(tmp_path_factory, roster_and_doc):
    """The PE toy's member-case table built in 1, 2 and 3 parts, each part a
    job in its own process, and merged; with each K's build seconds (summed
    over its parts) and table bytes."""
    if not (DATABASE / "input/kinetics/families").is_dir():
        pytest.skip("RMG_DATABASE_PATH is not a database checkout")
    roster, doc = roster_and_doc
    work = tmp_path_factory.mktemp("tables")
    raw = json.dumps(doc, sort_keys=True, separators=(",", ":")).encode()
    art = work / "art" / f"{hashlib.sha256(raw).hexdigest()}.json"
    art.parent.mkdir()
    art.write_bytes(raw)
    builder = _builder()
    tables, info = {}, {}
    for K in (1, 2, 3):
        parts, seconds = [], 0.0
        for k in range(K):
            part = work / f"part{k}of{K}.table"
            t = time.time()
            code, text, err = _script(work, [BUILDER, "build", art, roster, "--database",
                                             DATABASE, "--subsets", f"{k}/{K}", "--out", part])
            seconds += time.time() - t
            assert code == 0, err[-4000:]
            assert re.search(rf"^MEMBER part {k}/{K}: \d+ of \d+ subsets$", text, re.M), text
            parts.append(part)
        tables[K] = work / f"table{K}.table"
        assert builder.main(["merge", *map(str, parts), "--out", str(tables[K])]) == 0
        info[K] = (round(seconds, 1), tables[K].stat().st_size)
    return tables, info


@needs_database
@pytest.mark.parametrize("path", ["whole", "prepared", "blocks"])
def test_a_table_built_in_parts_gives_every_shard_the_in_shard_cases(
        tmp_path, artifact, roster_and_doc, real_tables, monkeypatch, capsys, path):
    """Real reference systems, database and member cases (emission is the
    digest of the record; the cases, stats and species it reads are compared): a
    table built in 1, 2 or 3 jobs, each in its own process, and merged gives
    the shard the in-shard build's cases, stats, species registry and rows,
    with no subset computed in the shard."""
    roster = roster_and_doc[0]
    flags = ["--chunk", "0/2", "--limit", "2"]
    if path == "prepared":
        flags += ["--prepared", str(_prepare(artifact, roster, tmp_path, 2))]
    elif path == "blocks":
        flags += ["--prepared", str(_blocks(artifact, roster, tmp_path, 2, shards=2) / "current")]
    seen = []

    def spy(self, record):
        S = self.state
        view = _case_view(S["cases"], S["case_stats"], S["species"].__dict__)
        seen.append(view)
        # (no digest of the view in the row: real objects' reprs name their ids)
        return digest_row(self, record)

    monkeypatch.setattr(emit.Emitter, "run_record", spy)
    capsys.readouterr()
    base = tmp_path / "in-shard.jsonl"
    assert _resume(artifact, roster, base, flags) == 0
    assert "member-case table" not in capsys.readouterr().out
    assert seen and seen[0][0], "the in-shard build kept no case"
    expected, rows = seen[0], _stripped(base.read_bytes())
    tables, info = real_tables
    for K, table in tables.items():
        seen.clear()
        out = tmp_path / f"with-table{K}.jsonl"
        assert _resume(artifact, roster, out, [*flags, "--member-cases", str(table)]) == 0
        text = capsys.readouterr().out
        read, computed, enumerated = _table_line(text)
        assert computed == 0 and read == enumerated > 0, text
        assert seen and all(view == expected for view in seen), K
        assert _stripped(out.read_bytes()) == rows, K
    with capsys.disabled():
        print(f"\nI-329 toy table ({path}): " + ", ".join(
            f"K={K}: {s}s build, {b} bytes" for K, (s, b) in info.items()))
