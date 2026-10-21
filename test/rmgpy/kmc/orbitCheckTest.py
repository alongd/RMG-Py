"""Tests for the orbit check of rmgpy/kmc/template_emission.py
(``orbit_problem``): whether a channel's transitions are one orbit of the
witness reactants' automorphism group.

A symmetric hydrogen-rich witness (reduced from an emitted alkane
H-abstraction witness) is decided within a time bound; on a corpus of orbit
decisions (every call of local emissions of polyethylene and polypropylene
records, and synthetic witnesses: phenyl rings, heteroatoms, charges and lone
pairs, radical sites, hydrogens breaking a heavy-atom symmetry, an A + A pair,
H2 and asymmetric controls) the check gives the decision of the full graph,
recorded in the fixture and recomputed here with a full-graph matcher where
that is fast.
"""

import gzip
import json
import signal
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path

import pytest

from rmgpy.kmc import melt_state as ms
from rmgpy.kmc import template_emission as te

FIXTURES = Path(__file__).parent / "fixtures" / "orbit_check"


class Late(Exception):
    pass


@contextmanager
def time_bound(seconds):
    def alarm(*_):
        raise Late(f"over {seconds} s")
    previous = signal.signal(signal.SIGALRM, alarm)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)


def _case(doc):
    return te.rows_from_json(doc["rows"]), [tuple(getattr(ms, e[0])(*e[1:]) for e in es) for es in doc["edits"]]


def _first_unrelated(problem):
    """The transition the check names as not related to transition 0 (None: one orbit)."""
    return None if problem is None else int(problem.split("transition ")[-1].split(" ")[0])


def _full_first_unrelated(rows, edit_sets):
    """The full-graph decision: every atom (hydrogens included) a node,
    labelled with its row and the edits written on it, matched by networkx
    VF2 (the reference; exponential on symmetric hydrogen-rich witnesses)."""
    import networkx as nx
    from networkx.algorithms import isomorphism as iso

    def annotated(edits):
        g = nx.Graph()
        for a, row in rows.items():
            g.add_node(a, lab=repr(tuple(row[:5])))
        for a, row in rows.items():
            for n, o in row[5]:
                if a < n:
                    g.add_edge(a, n, order=float(o), elab=repr((float(o), None)))
        marks = defaultdict(list)
        for e in edits:
            if isinstance(e, ms.SetBond):
                before = g.edges[e.a, e.b]["order"] if g.has_edge(e.a, e.b) else 0.0
                g.add_edge(e.a, e.b, order=before, elab=repr((before, float(e.new))))
            else:
                marks[e.atom].append((type(e).__name__, e.new))
        for a, m in marks.items():
            g.nodes[a]["lab"] = g.nodes[a]["lab"] + repr(tuple(sorted(m)))
        return g

    first = annotated(edit_sets[0])
    for k, edits in enumerate(edit_sets[1:], 1):
        if not iso.GraphMatcher(first, annotated(edits), node_match=lambda x, y: x["lab"] == y["lab"],
                                edge_match=lambda x, y: x["elab"] == y["elab"]).is_isomorphic():
            return k
    return None


def test_the_orbit_check_of_a_symmetric_hydrogen_rich_witness_is_fast():
    """n-Pentadecane + 1-butyl, H abstraction from the mirror-image carbons 7
    and 9 (the atom ids of the emitted witness): one orbit, decided well
    within the bound (the full-graph matcher takes over ten seconds here and
    did not finish in minutes on the emitted 78- and 84-atom witnesses)."""
    rows, edit_sets = _case(json.loads((FIXTURES / "slow_witness.json").read_text()))
    with time_bound(5):
        assert te.orbit_problem(rows, edit_sets) is None


def _corpus():
    with gzip.open(FIXTURES / "corpus.json.gz", "rt") as fh:
        return json.load(fh)


def test_the_corpus_holds_both_decisions_and_every_kind_of_case():
    corpus = _corpus()
    assert len(corpus) >= 300
    assert {c["full"] is None for c in corpus} == {True, False}
    sources = {c["source"].split()[0] for c in corpus}
    assert {"pe", "pp_disp", "pp_recomb", "pp_intra", "pp_habs", "synthetic"} <= sources
    assert {c["method"] for c in corpus} == {"legacy", "vf2pp", "lift"}


def test_the_orbit_check_gives_the_full_graph_decision_on_the_corpus():
    """Per case, the check's decision (the first transition not related to
    transition 0) is the full graph's, as recorded (networkx VF2, else VF2++,
    else a verified lift of the folded isomorphism to the full graph); where
    the full-graph matcher took under 50 ms it is recomputed here."""
    wrong, live = [], 0
    for i, doc in enumerate(_corpus()):
        rows, edit_sets = _case(doc)
        with time_bound(5):
            got = _first_unrelated(te.orbit_problem(rows, edit_sets))
        want = doc["full"]
        if doc["method"] == "legacy" and doc["legacy_seconds"] < 0.05:
            live += 1
            assert _full_first_unrelated(rows, edit_sets) == want, (i, doc["source"])
        if got != want:
            wrong.append((i, doc["source"], got, want))
    assert wrong == []
    assert live >= 250
