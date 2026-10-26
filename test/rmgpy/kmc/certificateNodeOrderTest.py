"""A stored record certificate re-derives whatever atom order RMG's matcher left
in the family tree's groups (I-251).

The certificate records each tree node it evaluates (``nodes``), a group as its
atoms in their in-memory order.  RMG's matcher re-sorts a group's atoms in place
(by connectivity) when it matches it, so a record certified before some other
record's match and re-derived after it carries the same node with its atoms
permuted.  Emission's recertification then refused it ``refused:recertify``
("the production binding rerun ... differs in ['certificate']"): on the PP
catalog, 12 intra_H_migration records, the first few each shard certified, and
their reciprocal pairs became one-way channels at the cycle check.  The
stored node is the same node: ``same_node_semantics`` already accepts it for the
tree check, and the whole-certificate comparison must too, while a node whose
content changed is still refused.
"""

import json
import os
from pathlib import Path

import pytest

from rmgpy.kmc import compile as compile_module

REPOSITORY = Path(__file__).resolve().parents[3]
DECKS = REPOSITORY / "test/rmgpy/kmc/fixtures/decks"
DATABASE = Path(os.environ.get("RMG_DATABASE_PATH", "/nonexistent"))
needs_database = pytest.mark.skipif(
    not (DATABASE / "input/kinetics/families").is_dir(),
    reason="RMG_DATABASE_PATH is not a database checkout",
)
FAMILIES = ("intra_H_migration",)
SITES = ("interior_radical:0", "head_end_radical")


@pytest.fixture(scope="module")
def certified(tmp_path_factory):
    """The PE toy deck's intra_H_migration records on two sites, compiled, then
    certified as emission certifies them (``Emitter.init``)."""
    from rmgpy.kmc import emit
    from rmgpy.kmc.run import read_deck
    deck = DECKS / "pe_toy_deck.py"
    original = compile_module._proxy_inputs

    def subset(request):
        proxies, source, options, discovery = original(request)
        return [p for p in proxies if p.site_type in SITES], source, options, discovery

    compile_module._proxy_inputs = subset
    try:
        path, _ = compile_module.compile_artifact(compile_module.CompileRequest(
            database_path=str(DATABASE), output_dir=str(tmp_path_factory.mktemp("art")),
            families=FAMILIES, deck=str(deck), database_sha="certificate-node-order"))
        emitter = emit.Emitter(path, DATABASE, {}, families=FAMILIES,
                               polymer=read_deck(deck).spec)
        emitter.init()
    finally:
        compile_module._proxy_inputs = original
    return emitter.state


def _groups(state, record, body):
    """(label, group) of every group node the record's stored certificate evaluates."""
    entries = state["fam"][record["family"]].groups.entries
    return [(label, entries[label].item) for label, node in
            sorted(body["certificate"].get("nodes", {}).items())
            if node.get("kind") == "group" and label in entries]


def _fresh(certifier):
    """A certifier like ``certifier`` with empty caches: what a long certification
    run gives a record whose cached descents it evicted (the caches are bounded,
    least recently used)."""
    return type(certifier)(certifier.families, resonance=certifier.resonance,
                           max_channel_sets=certifier.max_channel_sets, domain=certifier.domain,
                           repeat=certifier.repeat, growth=certifier.growth,
                           thermo=certifier.thermo)


def _resort(group, i, j):
    """Swap atoms ``i`` and ``j`` of ``group`` as a re-sort by RMG's matcher leaves
    them (``Graph.sort_vertices``: the list order, numbered by sorting label)."""
    group.atoms[i], group.atoms[j] = group.atoms[j], group.atoms[i]
    for index, atom in enumerate(group.atoms):
        atom.sorting_label = index


def _reordered(state):
    """An enabled record, its stored certificate, and one of its group nodes whose
    atoms are re-sorted in place (two swapped) so that the certificate re-derived
    without caches lists that node's atoms in another order.  Returns (record,
    body, label, restore)."""
    from rmgpy.kmc import certificates as cm
    from rmgpy.kmc.emit import record_certificate
    for record in state["records"]:
        body = record_certificate(record, state["store"])
        if record["status"] != "enabled" or not body or not body.get("certificate"):
            continue
        for label, group in _groups(state, record, body):
            if len(group.atoms) < 2:
                continue
            order = list(group.atoms)
            labels = [a.sorting_label for a in order]

            def restore(g=group, o=order, s=labels):
                g.atoms[:] = o
                for atom, lab in zip(o, s):
                    atom.sorting_label = lab
            _resort(group, 0, 1)
            fresh = cm.rederived_certificate(record, state["by_id"], _fresh(state["certifier"]))
            if (fresh.get("certificate") or {}).get("nodes", {}).get(label) != \
                    body["certificate"]["nodes"][label]:
                return record, body, label, restore
            restore()
    raise AssertionError("no group node whose re-derived atom order a re-sort changes")


def _problems(state, record, body):
    from rmgpy.kmc import certificates as cm
    return cm.stored_certificate_problems(record, body, state["by_id"],
                                          _fresh(state["certifier"]))


@needs_database
def test_a_node_whose_atoms_the_matcher_reordered_still_re_derives(certified):
    state = certified
    record, body, label, restore = _reordered(state)
    try:
        assert _problems(state, record, body) == []
    finally:
        restore()
    assert _problems(state, record, body) == []


@needs_database
def test_a_node_whose_content_changed_is_still_refused(certified):
    state = certified
    record, body, label, restore = _reordered(state)
    restore()
    group = state["fam"][record["family"]].groups.entries[label].item
    atom = next(a for a in group.atoms if a.label)
    before = list(atom.radical_electrons)
    atom.radical_electrons = before + [max(before or [0]) + 7]
    try:
        problems = _problems(state, record, body)
    finally:
        atom.radical_electrons = before
    assert problems and any(label in p or "certificate" in p for p in problems), problems
