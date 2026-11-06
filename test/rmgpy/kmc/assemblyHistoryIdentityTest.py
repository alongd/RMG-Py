#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""I-308: the identity that orders every catalog-assembly choice (the
representative of a collision, the source pair of a reciprocal unit, the
cycle-check temperatures) does not depend on the history in which the emit
process certified the member contexts.

A template's ``certificate`` is the sha256 of a member context's observation
certificate, certified inside the emit process; its node atom order follows
RMG's in-place group sort at a group's first match in that process, so its
bytes (and the payload ``sha256`` over them) move with the shard count and
the record order (pk-i304: novolac fs-ihm whole N=2 vs N=3, 3 rows; PE toy,
4).  Two payload sets that differ only there must assemble alike."""

import hashlib
import itertools
import json

from rmgpy.kmc import emit
from rmgpy.kmc import template_emission as te

import templateEmissionTest as tet

SALTS = [f"history-{k}" for k in range(12)]


def _rehistory(payloads, salt):
    """``payloads`` as another emit history would store them: every
    template's ``certificate`` replaced by a salted digest, resealed (so
    ``payload.sha256`` moves too); nothing else changed."""
    out = []
    for p in payloads:
        p = json.loads(json.dumps(p))
        for j, t in enumerate(p.get("templates", [])):
            t["certificate"] = hashlib.sha256(
                f"{salt}|{p['record_event_id']}|{j}".encode()).hexdigest()
        out.append(te.finish(p))
    return out


def _history_free(value):
    """``value`` (payloads, or a table) as canonical JSON without the two
    history-dependent digests, for comparing outputs."""
    value = json.loads(json.dumps(value))

    def walk(v):
        if isinstance(v, dict):
            v.pop("sha256", None)
            for t in v.get("templates", []) if isinstance(v.get("templates"), list) else []:
                if isinstance(t, dict):
                    t.pop("certificate", None)
            for x in v.values():
                walk(x)
        elif isinstance(v, list):
            for x in v:
                walk(x)
    walk(value)
    return te.canonical_json(value)


def _diff_paths(a, b, path=""):
    if type(a) != type(b):
        return {path}
    if isinstance(a, dict):
        return set().union(*[_diff_paths(a.get(k), b.get(k), f"{path}.{k}")
                             for k in set(a) | set(b)])
    if isinstance(a, list):
        if len(a) != len(b):
            return {path + "[len]"}
        return set().union(*[_diff_paths(x, y, path + "[]") for x, y in zip(a, b)])
    return set() if a == b else {path}


def test_the_fixtures_differ_in_exactly_the_history_digests():
    """Non-vacuity: the rehistoried payloads differ from the originals at
    ``sha256`` and ``templates[].certificate`` and nowhere else, as the
    pk-i304 rows do."""
    payloads, _counter = tet._two_pairs(0.004)
    for salt in SALTS[:3]:
        moved = _rehistory(payloads, salt)
        paths = set().union(*[_diff_paths(a, b) for a, b in zip(payloads, moved)])
        assert paths == {".sha256", ".templates[].certificate"}


def test_the_identity_does_not_follow_the_template_certificates():
    payloads, _counter = tet._toy_pair()
    p = payloads[0]
    for salt in SALTS:
        q = _rehistory([p], salt)[0]
        assert q["sha256"] != p["sha256"]
        assert te.partition_free_identity(q) == te.partition_free_identity(p)
        assert te.partition_free_identity(q, tet.base_template(q)) == \
            te.partition_free_identity(p, tet.base_template(p))


def test_the_identity_still_follows_every_other_template_field():
    """Partition- and history-free, but not content-free: changing any stored
    template field other than the certificate changes both identities."""
    payloads, _counter = tet._toy_pair()
    p = payloads[0]
    t = tet.base_template(p)
    fields = sorted(set(t) - set(te.TEMPLATE_HISTORY_FIELDS))
    assert len(fields) > 10
    for field_ in fields:
        q = json.loads(json.dumps(p))
        tet.base_template(q)[field_] = {"changed": field_}
        assert te.partition_free_identity(q) != te.partition_free_identity(p), field_
        assert te.partition_free_identity(q, tet.base_template(q)) != \
            te.partition_free_identity(p, t), field_


def test_the_collision_representative_does_not_follow_the_certificate_history():
    """Three copies 0.004 decades apart merge into one copy's weight; under
    every certificate history the merged payloads and the collision table are
    the same, in any input order."""
    members = tet._group([0.0, 0.004, 0.008])
    weights = {te.canonical_json(tet.base_template(p)["weight"]) for p in members}
    assert len(weights) == 3      # a different representative would show
    seen_out, seen_table = set(), set()
    for salt in SALTS:
        moved = _rehistory(members, salt)
        for order in itertools.permutations(moved):
            out, table = te.resolve_collisions(list(order))
            assert [r["disposition"] for r in table] == ["merged"]
            seen_out.add(_history_free(sorted(out, key=lambda p: p["record_event_id"])))
            seen_table.add(_history_free(table))
    assert len(seen_out) == 1 and len(seen_table) == 1


def test_the_representative_pair_of_a_reciprocal_unit_does_not_follow_the_certificate_history():
    kf, kr = "a" * 64, "b" * 64
    payloads = tet._reciprocal([kf, kf], [kr, kr], (0.003, 0.0))
    seen = set()
    for salt in SALTS:
        out, table = te.resolve_collisions(_rehistory(payloads, salt))
        assert sorted(r["disposition"] for r in table) == ["merged", "merged"]
        seen.add((_history_free(sorted(out, key=lambda p: p["record_event_id"])),
                  _history_free(table)))
    assert len(seen) == 1


def test_the_cycle_check_temperatures_do_not_follow_the_certificate_history(monkeypatch):
    payloads, counter = tet._toy_pair()
    original = te.check_temperatures
    used = []

    def spy(grid, pair, extra=()):
        used[-1].add(tuple(pair))
        return original(grid, pair, extra)
    monkeypatch.setattr(te, "check_temperatures", spy)
    for salt in SALTS:
        used.append(set())
        _out, table = tet._check(tet._res(_rehistory(payloads, salt), counter), counter)
        assert table[0]["status"] == "ok"
    assert used[0] and all(u == used[0] for u in used)


def test_the_assembled_catalog_does_not_follow_the_certificate_history():
    """The assembly chain (collisions, reciprocal weights, cycle check) on two
    pairs that merge as one unit: the loaded catalog, the collision table and
    the cycle table are the same under every certificate history."""
    payloads, counter = tet._two_pairs(0.004)
    seen = set()
    for salt in SALTS:
        out, collisions = te.resolve_collisions(_rehistory(payloads, salt))
        out = te.derive_weights(out, counter=counter)
        out, cycles = tet._check(out, counter)
        assert sorted(r["disposition"] for r in collisions) == ["merged", "merged"]
        assert cycles
        seen.add((te.canonical_json(emit.stripped_payloads(out)),
                  _history_free(collisions), _history_free(cycles)))
    assert len(seen) == 1


def test_copies_with_different_chemistry_still_order_deterministically():
    """Copies that differ in stored content (here their weights) order by
    that content alone: the same representative under every certificate
    history, every input order and every assignment of record ids."""
    members = tet._group([0.0, 0.004, 0.008])
    ids = [p["record_event_id"] for p in members]
    loaded = set()
    for salt, perm in zip(itertools.cycle(SALTS), itertools.permutations(ids)):
        moved = _rehistory(tet._relabel_records(members, dict(zip(ids, perm))), salt)
        for order in (moved, moved[::-1]):
            out, _table = te.resolve_collisions(list(order))
            live = [tet.base_template(p) for p in out
                    if tet.base_template(p)["status"] == "emitted"]
            assert len(live) == 1
            loaded.add(te.canonical_json(live[0]["weight"]))
    assert len(loaded) == 1
