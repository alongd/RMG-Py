"""Policy tests for the declaration-generated reference mapping checker."""

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace


SCRIPT = Path(__file__).parents[3] / "scripts" / "check-kmc-reference-mapping.py"
SPEC = importlib.util.spec_from_file_location("check_kmc_reference_mapping", SCRIPT)
CHECKER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECKER)


def test_firing_statuses_are_not_named_refusals_even_with_a_reason():
    assert CHECKER.named_refusal({"status": "enabled", "status_reason": "diagnostic"}) is None
    assert CHECKER.named_refusal({"status": "irreversible",
                                  "status_reason": "RMG reaction is irreversible"}) is None
    assert CHECKER.named_refusal({"status": "refused", "status_reason": "policy"}) == \
        "refused:record_status:refused"


def test_numeric_proxy_depths_are_added_from_artifact_site_types():
    records = [{"site_type": "pristine"}, {"site_type": "end_radical@5"},
               {"site_type": "other@12"}, {"site_type": "not-a-depth@x"}]
    assert CHECKER.reference_dps(records, (1, 2, 3, 8)) == (1, 2, 3, 5, 8, 12)


def test_proxy_depth_generates_the_corresponding_declared_chain():
    from rmgpy.kmc import emit
    from melt_engine_support import pe_mapping
    from rmgpy.polymer_spec import LinearHomopolymerSpec

    spec = LinearHomopolymerSpec.from_mapping(pe_mapping({8: 1}))
    dps = CHECKER.reference_dps([{"site_type": "end_radical@5"}], (1, 2, 3, 8))
    roster = emit.reference_systems(spec, emit.match_provider(spec.repeat_unit), dps=dps,
                                    small=False, monomer=False, molecularity={1})
    assert "single:chain:dp5:capped" in roster.systems


def test_zero_matching_cases_for_required_records_exits_nonzero(tmp_path, monkeypatch):
    records = [{"event_id": f"event-{index}", "status": "enabled",
                "site_type": "pristine", "orientation": "as_generated",
                "family": "family"}
               for index in range(762)]
    artifact = tmp_path / "artifact.json"
    artifact.write_text(json.dumps({"records": records}))
    spec = SimpleNamespace(repeat_unit=object())
    monkeypatch.setattr(CHECKER, "read_deck", lambda _path: SimpleNamespace(spec=spec))
    monkeypatch.setattr(CHECKER.emit, "match_provider", lambda _graph: object())
    monkeypatch.setattr(CHECKER.emit, "reference_systems", lambda *_args, **_kwargs:
                        SimpleNamespace(systems={}, sha256="hash"))
    monkeypatch.setattr(CHECKER.emit, "build_reference_systems", lambda _roster: ({}, []))
    monkeypatch.setattr(CHECKER.emit, "load_database", lambda *_args:
                        SimpleNamespace(kinetics=SimpleNamespace(families={})))
    monkeypatch.setattr(CHECKER.cm, "Certifier", lambda *_args, **_kwargs: object())
    monkeypatch.setattr(CHECKER.chm, "SpeciesKeys", lambda: object())
    monkeypatch.setattr(CHECKER.emit, "reference_cases", lambda *_args: ({}, {}))
    monkeypatch.setattr(CHECKER.emit, "record_key", lambda *_args: ("missing",))
    monkeypatch.setattr(CHECKER, "prior_unmapped",
                        lambda _path: {record["event_id"] for record in records})
    args = SimpleNamespace(artifact=artifact, deck="deck.py", database="database",
                           old_rows="rows.jsonl", roster_out=None, dps=(1, 2, 3, 8),
                           limit=None)
    assert CHECKER.check(args) == 1


def test_one_prior_match_does_not_hide_other_required_records(tmp_path, monkeypatch, capsys):
    records = [{"event_id": f"event-{index}", "status": "enabled",
                "site_type": "pristine", "orientation": "as_generated",
                "family": "family"}
               for index in range(762)]
    artifact = tmp_path / "artifact.json"
    artifact.write_text(json.dumps({"records": records}))
    spec = SimpleNamespace(repeat_unit=object())
    monkeypatch.setattr(CHECKER, "read_deck", lambda _path: SimpleNamespace(spec=spec))
    monkeypatch.setattr(CHECKER.emit, "match_provider", lambda _graph: object())
    monkeypatch.setattr(CHECKER.emit, "reference_systems", lambda *_args, **_kwargs:
                        SimpleNamespace(systems={}, sha256="hash"))
    monkeypatch.setattr(CHECKER.emit, "build_reference_systems", lambda _roster: ({}, []))
    monkeypatch.setattr(CHECKER.emit, "load_database", lambda *_args:
                        SimpleNamespace(kinetics=SimpleNamespace(families={})))
    monkeypatch.setattr(CHECKER.cm, "Certifier", lambda *_args, **_kwargs: object())
    monkeypatch.setattr(CHECKER.chm, "SpeciesKeys", lambda: object())
    monkeypatch.setattr(CHECKER.emit, "reference_cases",
                        lambda *_args: ({("found",): [object()]}, {}))
    monkeypatch.setattr(CHECKER.emit, "record_key",
                        lambda _certifier, record, _species:
                        ("found",) if record["event_id"] == "event-0" else ("missing",))
    monkeypatch.setattr(CHECKER, "prior_unmapped", lambda _path: {"event-0"})
    args = SimpleNamespace(artifact=artifact, deck="deck.py", database="database",
                           old_rows="rows.jsonl", roster_out=None, dps=(1, 2, 3, 8),
                           limit=None)

    assert CHECKER.check(args) == 1
    summary = json.loads(capsys.readouterr().out.splitlines()[-1])
    assert (summary["found"], summary["not_found"], summary["named_refused"]) == (1, 761, 0)


def test_status_refusal_stays_named_when_its_record_key_raises(tmp_path, monkeypatch, capsys):
    records = [{"event_id": f"event-{index}",
                "status": "refused" if index == 0 else "enabled",
                "site_type": "pristine", "orientation": "as_generated",
                "family": "family"}
               for index in range(762)]
    artifact = tmp_path / "artifact.json"
    artifact.write_text(json.dumps({"records": records}))
    spec = SimpleNamespace(repeat_unit=object())
    monkeypatch.setattr(CHECKER, "read_deck", lambda _path: SimpleNamespace(spec=spec))
    monkeypatch.setattr(CHECKER.emit, "match_provider", lambda _graph: object())
    monkeypatch.setattr(CHECKER.emit, "reference_systems", lambda *_args, **_kwargs:
                        SimpleNamespace(systems={}, sha256="hash"))
    monkeypatch.setattr(CHECKER.emit, "build_reference_systems", lambda _roster: ({}, []))
    monkeypatch.setattr(CHECKER.emit, "load_database", lambda *_args:
                        SimpleNamespace(kinetics=SimpleNamespace(families={})))
    monkeypatch.setattr(CHECKER.cm, "Certifier", lambda *_args, **_kwargs: object())
    monkeypatch.setattr(CHECKER.chm, "SpeciesKeys", lambda: object())
    monkeypatch.setattr(CHECKER.emit, "reference_cases",
                        lambda *_args: ({("found",): [object()]}, {}))

    def record_key(_certifier, record, _species):
        if record["event_id"] == "event-0":
            raise ValueError("a refused record needs no key")
        return ("found",)

    monkeypatch.setattr(CHECKER.emit, "record_key", record_key)
    monkeypatch.setattr(CHECKER, "prior_unmapped", lambda _path: {"event-0"})
    args = SimpleNamespace(artifact=artifact, deck="deck.py", database="database",
                           old_rows="rows.jsonl", roster_out=None, dps=(1, 2, 3, 8),
                           limit=None)

    assert CHECKER.check(args) == 0
    summary = json.loads(capsys.readouterr().out.splitlines()[-1])
    assert (summary["found"], summary["not_found"], summary["named_refused"]) == (761, 1, 1)
    assert summary["missing"] == []


# -- I-191: a firing record is one production certifies and does not refuse upstream


def _records(refused_status=()):
    return [{"event_id": f"event-{index}",
             "status": "refused" if index in refused_status else "enabled",
             "site_type": "pristine", "orientation": "as_generated", "family": "family",
             "rate_source": {"kind": "RMG family estimate"}}
            for index in range(762)]


def _certify(refuse=(), seen=None):
    """``certify_records`` with ``required``: a refused record of the required
    population is republished REFUSED by status under a new event id."""
    def certify(records, _families, _certifier, required=True):
        assert required
        if seen is not None:
            seen.extend(record["event_id"] for record in records)
        out, store, replacements = [], {}, {}
        for record in records:
            record = dict(record)
            certified = record["event_id"] not in refuse
            if not certified:
                replacements[record["event_id"]] = record["event_id"] + "-refused"
                record.update(status="refused", event_id=record["event_id"] + "-refused",
                              status_reason="rate certification refused: no tree descent")
            sha = "sha-" + record["event_id"]
            store[sha] = {"certified": certified, "gaps": [] if certified else ["no tree descent"]}
            out.append(dict(record, rate_certificate={"sha256": sha}))
        return out, store, replacements
    return certify


def _check(tmp_path, monkeypatch, records, *, certify, record_key, prior, limit=None):
    from rmgpy.kmc.event_record import EventRecord

    artifact = tmp_path / "artifact.json"
    artifact.write_text(json.dumps({"records": records}))
    spec = SimpleNamespace(repeat_unit=object())
    monkeypatch.setattr(CHECKER, "read_deck", lambda _path: SimpleNamespace(spec=spec))
    monkeypatch.setattr(CHECKER.emit, "match_provider", lambda _graph: object())
    monkeypatch.setattr(CHECKER.emit, "reference_systems", lambda *_args, **_kwargs:
                        SimpleNamespace(systems={}, sha256="hash"))
    monkeypatch.setattr(CHECKER.emit, "build_reference_systems", lambda _roster: ({}, []))
    monkeypatch.setattr(CHECKER.emit, "load_database", lambda *_args:
                        SimpleNamespace(kinetics=SimpleNamespace(families={})))
    monkeypatch.setattr(CHECKER.cm, "Certifier", lambda *_args, **_kwargs: object())
    monkeypatch.setattr(CHECKER.cm, "certify_records", certify)
    monkeypatch.setattr(EventRecord, "from_dict", staticmethod(dict))
    monkeypatch.setattr(CHECKER.chm, "SpeciesKeys", lambda: object())
    monkeypatch.setattr(CHECKER.emit, "reference_cases",
                        lambda *_args: ({("found",): [object()]}, {}))
    monkeypatch.setattr(CHECKER.emit, "record_key",
                        lambda _certifier, record, _species: record_key(record))
    monkeypatch.setattr(CHECKER, "prior_unmapped", lambda _path: set(prior))
    args = SimpleNamespace(artifact=artifact, deck="deck.py", database="database",
                           old_rows="rows.jsonl", roster_out=None, dps=(1, 2, 3, 8),
                           limit=limit)
    return CHECKER.check(args)


def test_an_enabled_record_certification_refuses_is_a_named_refusal(tmp_path, monkeypatch,
                                                                     capsys):
    records = _records()
    status = _check(tmp_path, monkeypatch, records, certify=_certify(refuse={"event-0"}),
                    record_key=lambda record: ("missing",) if record["event_id"].startswith(
                        "event-0") else ("found",),
                    prior={"event-0"})
    summary = json.loads(capsys.readouterr().out.splitlines()[-1])
    assert status == 0
    assert (summary["found"], summary["uncovered"], summary["named_refused"]) == (761, 0, 1)
    assert summary["named_refused_by_status"] == {"refused:record_status": 1}
    assert summary["named_refused_by_reason"] == {
        "refused:record_status: rate certification refused: no tree descent": 1}
    assert summary["missing"] == []


def test_a_certified_record_without_a_reference_system_still_exits_nonzero(tmp_path,
                                                                           monkeypatch, capsys):
    records = _records()
    status = _check(tmp_path, monkeypatch, records, certify=_certify(refuse={"event-0"}),
                    record_key=lambda record: ("missing",) if record["event_id"] == "event-1"
                    else ("found",),
                    prior={"event-1"})
    summary = json.loads(capsys.readouterr().out.splitlines()[-1])
    assert status == 1
    assert (summary["found"], summary["uncovered"], summary["named_refused"]) == (760, 1, 1)
    assert summary["missing"] == ["event-1"]
    assert summary["prior_unmapped_missing"] == ["event-1"]


def test_a_firing_record_outside_the_required_population_is_refused_as_production_does(
        tmp_path, monkeypatch, capsys):
    from rmgpy.kmc import template_emission as te

    records = _records()
    records[2]["rate_source"] = {"kind": "k_library"}
    status = _check(tmp_path, monkeypatch, records, certify=_certify(),
                    record_key=lambda record: ("missing",) if record["event_id"] == "event-2"
                    else ("found",),
                    prior={"event-2"})
    summary = json.loads(capsys.readouterr().out.splitlines()[-1])
    assert status == 0
    why = te.upstream_refusal(dict(records[2], rate_certificate={"sha256": "sha-event-2"}),
                              {"certified": True, "gaps": []},
                              {record["event_id"]: record for record in records})
    assert summary["named_refused_by_reason"] == {why: 1}
    assert summary["named_refused_by_status"] == {"refused:certificate": 1}


def test_limit_certifies_every_record_and_reports_only_the_selection(tmp_path, monkeypatch,
                                                                     capsys):
    records = _records()
    seen = []
    status = _check(tmp_path, monkeypatch, records, certify=_certify(refuse={"event-1"}, seen=seen),
                    record_key=lambda record: ("found",), prior={"event-1"}, limit=3)
    summary = json.loads(capsys.readouterr().out.splitlines()[-1])
    assert status == 0
    assert seen == [record["event_id"] for record in records]
    assert summary["records"] == 3
    assert (summary["found"], summary["uncovered"], summary["named_refused"]) == (2, 0, 1)
    assert summary["prior_unmapped_selected"] == 1
