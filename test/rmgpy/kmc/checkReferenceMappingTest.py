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
