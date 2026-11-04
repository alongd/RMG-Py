"""The streaming artifact reader (I-295).

Emission's prepare stage reads the artifact through ``artifact_stream``: one
pass indexes every record's byte span and links, and records are then read
one at a time by span.  The reader must give exactly what ``json.loads``
gives, at any block size, and refuse what is not an artifact.
"""
import hashlib
import json
import random

import pytest

from rmgpy.kmc import artifact_stream as astream


def _id(k):
    return "evt_" + hashlib.sha256(str(k).encode()).hexdigest()


def _artifact(n=12, seed=0):
    rng = random.Random(seed)
    records = []
    for k in range(n):
        partner = k + 1 if k % 2 == 0 else k - 1
        records.append({"event_id": _id(k), "family": f"F{k % 3}", "reverse_of": _id(partner),
                        "k_table": {"T": [300.0, 1e3],
                                    "k": [rng.random() * 10 ** rng.randint(-9, 9), -1.5e-300]},
                        "note": "é中 \"quoted\" \\ evt_" + "0" * 64,
                        "nested": [{"deep": [None, True, False, 12345678901234567890]}],
                        "status": "enabled"})
    return {"provenance": {"proxy_source": {"kind": "declaration", "pairs": {"1": []}}},
            "records": records, "schema_version": 19, "zz_after": [1.25, {"x": "y"}]}


def _write(path, doc, **dump):
    text = json.dumps(doc, **dump)
    path.write_bytes(text.encode("utf-8"))
    return path


DUMPS = [dict(sort_keys=True, separators=(",", ":"), ensure_ascii=True),
         dict(indent=2, ensure_ascii=False),
         dict(separators=(" , ", " : "), ensure_ascii=False)]


@pytest.mark.parametrize("dump", DUMPS)
@pytest.mark.parametrize("block", [1, 7, 64, 1 << 20])
def test_the_index_and_spans_reproduce_json_loads(tmp_path, dump, block):
    doc = _artifact()
    path = _write(tmp_path / "a.json", doc, **dump)
    index = astream.build_index(path, block=block)
    assert index.sha256 == hashlib.sha256(path.read_bytes()).hexdigest()
    assert index.size == path.stat().st_size
    assert index.head == {k: v for k, v in doc.items() if k != "records"}
    assert [e.event_id for e in index.entries] == [r["event_id"] for r in doc["records"]]
    assert astream.read_records(path, index.entries) == doc["records"]
    assert astream.read_records(path, index.entries[::-3]) == doc["records"][::-3]


def test_links_are_the_event_id_string_values_of_other_artifact_records(tmp_path):
    doc = _artifact(4)
    doc["records"][0]["deep"] = {"x": [_id(3)], _id(2): "a member name is not a link",
                                 "absent": _id(99)}
    index = astream.build_index(_write(tmp_path / "a.json", doc))
    assert index.entries[0].links == tuple(sorted({_id(1), _id(3)}))
    assert index.entries[2].links == (_id(3),)


@pytest.mark.parametrize("text, match", [
    ('{"provenance": {}}', "no records member"),
    ('[1, 2]', "expected one of"),
    ('{"records": [{"event_id": "a"}]} x', "data after"),
    ('{"records": [{"event_id": "a"}, {"event_id": "a"}]}', "two records of event id"),
    ('{"records": [{"family": "F"}]}', "has no event id"),
    ('{"records": [{"event_id": "a"}, ]}', "unreadable value"),
    ('{"records": [{"event_id": "a"}', "expected one of"),
])
def test_what_is_not_an_artifact_is_refused(tmp_path, text, match):
    path = tmp_path / "a.json"
    path.write_text(text)
    with pytest.raises(astream.ArtifactFormatError, match=match):
        astream.build_index(path, block=3)


def test_an_index_of_other_bytes_is_refused(tmp_path):
    path = _write(tmp_path / "a.json", _artifact(), **DUMPS[0])
    index_path = astream.write_index(astream.build_index(path), tmp_path / "a.index")
    assert astream.load_index(index_path, path).entries == astream.build_index(path).entries
    raw = bytearray(path.read_bytes())
    raw[-3] = ord("7") if raw[-3] != ord("7") else ord("8")
    path.write_bytes(bytes(raw))
    with pytest.raises(astream.IndexMismatch, match="refused:index"):
        astream.load_index(index_path, path)


def test_a_span_that_is_not_its_record_is_refused(tmp_path):
    path = _write(tmp_path / "a.json", _artifact(), **DUMPS[0])
    entry, other = astream.build_index(path).entries[1:3]
    with pytest.raises(astream.IndexMismatch, match="are not record"):
        astream.read_records(path, [astream.Entry(entry.event_id, None, other.offset,
                                                  other.length, ())])
