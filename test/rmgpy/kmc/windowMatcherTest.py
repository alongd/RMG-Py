"""Acceptance for the window matcher (rmgpy/kmc/window_matcher.py).

Every expected transition set comes from ``transition_oracle.py``, which
enumerates elementary moves from chemistry rules over the raw MeltState graph
and decides resonance equivalence by brute-force pi placement.  It shares no
code with the matcher.

Template sets.  Each family has a hand-written edit rule (``FAMILY_RULES``:
which atoms are edited and their final values).  Its windows are cut by hand
rule (edited atoms closed, their neighbours added, heavy neighbours as ports
with their exterior) around every move the oracle finds on small WITNESS
states (DP 6, 7, 9).  Isomorphic windows are one class.  The acceptance states
(DP 8, 12, 21, 30) are never used to build templates, so they test transfer
from short witnesses to long chains as well as the matcher itself.
"""

import itertools
import os
import pickle
import subprocess
import sys
import time
from collections import Counter
from dataclasses import replace
from fractions import Fraction
from pathlib import Path

import networkx as nx
import pytest

import transition_oracle
from transition_oracle import Oracle, canonical

from rmgpy.kmc import window_matcher as wm
from rmgpy.kmc.melt_state import (
    MeltState,
    Node,
    Origin,
    SetBond,
    SetRadical,
    materialize,
)
from rmgpy.kmc.window_matcher import (
    RESONANCE_FORM_LIMIT,
    AtomEdit,
    BondEdit,
    ConjugatedRingError,
    PairSampler,
    PathPredicate,
    ResonanceLimitError,
    StaleTransitionError,
    Template,
    TemplateError,
    WindowAtom,
    apply_transition,
    component_forms,
    form_rows,
    match,
    product_rows,
)
from rmgpy.polymer_spec import LinearHomopolymerSpec

TEST_DIR = Path(__file__).resolve().parent
ROOT = TEST_DIR.parents[2]
SLOW = os.environ.get("RMG_KMC_SLOW") == "1"

HA = "h_abstraction_benzylic"
BS = "beta_scission_midchain"
RC = "recombination_chain_end"
DS = "disproportionation"
HS = "h_shift"
HH = "h_atom_recombination"


# --------------------------------------------------------------------------
# states
# --------------------------------------------------------------------------
def _atom(key, element="C", hydrogens=0):
    return {"key": key, "element": element, "hydrogens": hydrogens,
            "electronic_state": {"radical_electrons": 0}}


def ps_spec(counts, tail="H"):
    """PS chains: head cap H; tail cap H or CH3 (asymmetric chain ends)."""
    ring = [_atom(f"p{i}", hydrogens=0 if i == 1 else 1) for i in range(1, 7)]
    ring_bonds = [{"atoms": [f"p{i}", f"p{i + 1}"], "order": 1.5} for i in range(1, 6)]
    ring_bonds.append({"atoms": ["p1", "p6"], "order": 1.5})
    repeat = {
        "atoms": [_atom("c1", hydrogens=2), _atom("c2", hydrogens=1)] + ring,
        "bonds": [{"atoms": ["c1", "c2"], "order": 1},
                  {"atoms": ["c2", "p1"], "order": 1}] + ring_bonds,
        "ports": {"head": {"atom": "c1", "bond_order": 1},
                  "tail": {"atom": "c2", "bond_order": 1}},
        "electronic_state": {"multiplicity": 1},
    }

    def cap(element, hydrogens):
        return {"atoms": [_atom("h", element, hydrogens)], "bonds": [],
                "ports": {"cap": {"atom": "h", "bond_order": 1}},
                "electronic_state": {"multiplicity": 1}}

    styrene_ring = [_atom(f"p{i}", hydrogens=0 if i == 1 else 1) for i in range(1, 7)]
    styrene = {
        "atoms": [_atom("a", hydrogens=2), _atom("b", hydrogens=1)] + styrene_ring,
        "bonds": [{"atoms": ["a", "b"], "order": 2},
                  {"atoms": ["b", "p1"], "order": 1}] + ring_bonds,
        "electronic_state": {"multiplicity": 1},
    }
    mapping = {
        "schema_version": 1, "label": "PS", "repeat_unit": repeat,
        "end_groups": {"head": cap("H", 0),
                       "tail": cap("H", 0) if tail == "H" else cap("C", 3)},
        "small_species": {"styrene": {"graph": styrene, "disposition": "mobile"}},
        "initial_ensemble": {"dp_counts": dict(counts)},
        "initial_mobile_counts": {},
        "density_basis": "whole_initial_mixture",
        "initial_mass_kg": 1.25, "density_kg_m3": 1050.0,
        "locality_policy": {"radius": 2, "bounds": {
            "elements": ["C", "H"], "isotopes": [-1], "formal_charges": [0],
            "radical_electrons": [0, 1], "lone_pairs": [0],
            "bond_orders": [1, 1.5, 2, 3], "spin_multiplicities": [1, 2, 3],
            "maximum_heavy_atom_degree": 4,
            "maximum_explicit_hydrogens_per_heavy_atom": 4,
            "maximum_total_hydrogens_on_free_small_species": 12,
            "maximum_family_arity": 2, "maximum_root_tuple_size": 4,
            "maximum_free_small_species_heavy_atoms": 8,
            "attachment_classes": [], "topology_metadata_schema": {}}},
        "family_manifest": {}, "library_manifest": {},
        "transport_model": {}, "diagnostics": {},
    }
    return LinearHomopolymerSpec.from_mapping(mapping)


def merge(*states):
    """One MeltState holding the molecules of several (ids renumbered in order,
    chain ordinals made unique)."""
    atoms, origin = {}, {}
    offset = chains = 0
    for state in states:
        ids = sorted(state.atoms)
        shift = {a: offset + k for k, a in enumerate(ids)}
        top = -1
        for a in ids:
            node = state.atoms[a]
            atoms[shift[a]] = Node(node.element, node.isotope, node.charge, node.radical,
                                   node.lone_pairs, {shift[b]: o for b, o in node.edges.items()})
            o = state.origin[a]
            ordinal = o.ordinal + chains if o.kind == "chain" else o.ordinal
            top = max(top, o.ordinal if o.kind == "chain" else -1)
            origin[shift[a]] = Origin(o.kind, ordinal, o.segment, o.key, o.is_backbone)
        offset += len(ids)
        chains += top + 1
    return MeltState(atoms, origin)


def find(state, chain, segment, key):
    hits = [a for a, o in state.origin.items()
            if (o.kind, o.ordinal, o.segment, o.key) == ("chain", chain, segment, key)]
    assert len(hits) == 1, (chain, segment, key, hits)
    return hits[0]


def remove_h(state, atom):
    """Homolyse the lowest-id C-H bond of ``atom``: a carbon radical plus a free H."""
    h = min(b for b in state.atoms[atom].edges if state.atoms[b].element == "H")
    state.apply([SetBond(atom, h, 1.0, 0.0), SetRadical(atom, 0, 1), SetRadical(h, 0, 1)])


def cleave(state, a, b):
    state.apply([SetBond(a, b, 1.0, 0.0), SetRadical(a, 0, 1), SetRadical(b, 0, 1)])


def rich_state(d):
    """The acceptance recipe at DP ``d``.

    chain 0 (DP d, H head / CH3 tail, asymmetric): two mid-chain benzylic
    radicals (units d//3 and 2d//3); chain 1 (DP d-1): head primary radical and
    tail benzylic chain-end radical; chain 2 (DP d): one secondary mid-chain
    radical (unit d//2); chain 3 (DP d+2): cut between units k and k+1, giving
    a benzylic chain end and a primary chain end.  Each H removal leaves a free
    H atom (5 in total).
    """
    state = merge(materialize(ps_spec({d: 1}, tail="CH3")),
                  materialize(ps_spec({d - 1: 1, d: 1, d + 2: 1})))
    remove_h(state, find(state, 0, f"unit {d // 3}", "c2"))
    remove_h(state, find(state, 0, f"unit {2 * d // 3}", "c2"))
    remove_h(state, find(state, 1, "unit 1", "c1"))
    remove_h(state, find(state, 1, f"unit {d - 1}", "c2"))
    remove_h(state, find(state, 2, f"unit {d // 2}", "c1"))
    k = (d + 2) // 2
    cleave(state, find(state, 3, f"unit {k}", "c2"), find(state, 3, f"unit {k + 1}", "c1"))
    return state


_STATES = {}


def state_at(d):
    if d not in _STATES:
        _STATES[d] = rich_state(d)
    return pickle.loads(pickle.dumps(_STATES[d]))


# --------------------------------------------------------------------------
# hand-built small molecules
# --------------------------------------------------------------------------
def molecules(*specs):
    """Build a state from (atoms, bonds) specs: atoms are (element, radical) or
    (element, radical, charge, lone_pairs); bonds are (i, j, order) local."""
    atoms, origin, offset = {}, {}, 0
    for mol, (atom_list, bonds) in enumerate(specs):
        edges = {offset + i: {} for i in range(len(atom_list))}
        for i, j, order in bonds:
            edges[offset + i][offset + j] = float(order)
            edges[offset + j][offset + i] = float(order)
        for i, spec in enumerate(atom_list):
            element, radical = spec[0], spec[1]
            charge, lone_pairs = (spec[2], spec[3]) if len(spec) > 2 else (0, 0)
            atoms[offset + i] = Node(element, None, charge, radical, lone_pairs, edges[offset + i])
            origin[offset + i] = Origin("mobile", mol, f"m{mol}", f"a{i}", False)
        offset += len(atom_list)
    return MeltState(atoms, origin)


def with_h(heavy, bonds, hydrogens):
    """Append explicit H atoms: ``hydrogens[i]`` H on heavy atom i."""
    atom_list = list(heavy)
    bonds = list(bonds)
    for i, n in enumerate(hydrogens):
        for _ in range(n):
            atom_list.append(("H", 0))
            bonds.append((i, len(atom_list) - 1, 1))
    return atom_list, bonds


def ethane():
    return with_h([("C", 0), ("C", 0)], [(0, 1, 1)], [3, 3])


def propane():
    return with_h([("C", 0)] * 3, [(0, 1, 1), (1, 2, 1)], [3, 2, 3])


def methyl():
    return with_h([("C", 1)], [], [3])


def benzyl(form="aromatic"):
    """Benzyl radical: 0 = CH2 carbon, 1 = ipso, 2/6 ortho, 3/5 meta, 4 para."""
    ring = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]
    if form == "aromatic":
        heavy = [("C", 1)] + [("C", 0)] * 6
        bonds = [(0, 1, 1)] + [(a, b, 1.5) for a, b in ring]
    else:
        radical_at = {"ortho": 2, "para": 4, "ortho2": 6}[form]
        heavy = [("C", 1 if i == radical_at else 0) for i in range(7)]
        doubles = {"ortho": {(3, 4), (5, 6)}, "para": {(2, 3), (5, 6)},
                   "ortho2": {(2, 3), (4, 5)}}[form]
        bonds = [(0, 1, 2)] + [(a, b, 2 if (a, b) in doubles else 1) for a, b in ring]
    return with_h(heavy, bonds, [2, 0, 1, 1, 1, 1, 1])


# --------------------------------------------------------------------------
# hand-written templates (unit tests)
# --------------------------------------------------------------------------
def methyl_acceptor_ethane_donor():
    """CH3* + H-CH2CH3: acceptor window = CH3* (closed); donor window = the
    H-bearing carbon (closed), its three H (closed) and the other carbon as a
    port whose exterior is three H."""
    atoms = [WindowAtom("C", radical=1), WindowAtom("H"), WindowAtom("H"), WindowAtom("H"),
             WindowAtom("C", participant=1), WindowAtom("H", participant=1),
             WindowAtom("H", participant=1), WindowAtom("H", participant=1),
             WindowAtom("C", participant=1, kind="port",
                        exterior=((1.0, "H", False),) * 3)]
    bonds = [(0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0), (4, 5, 1.0), (4, 6, 1.0), (4, 7, 1.0),
             (4, 8, 1.0)]
    edits = [AtomEdit(0, radical=0), BondEdit(0, 5, 1.0), BondEdit(4, 5, 0.0),
             AtomEdit(4, radical=1)]
    return Template("ch3_abstraction", HA, atoms, bonds, edits, weight=1.0)


def benzyl_acceptor_ethane_donor():
    """Benzyl radical as a fully closed window (ring included) + the ethane donor."""
    base = methyl_acceptor_ethane_donor()
    donor = [a for a in base.atoms if a.participant == 1]
    atoms = ([WindowAtom("C", radical=1)] + [WindowAtom("C")] * 6 + [WindowAtom("H")] * 7 + donor)
    bonds = [(0, 1, 1.0)] + [(a, b, 1.5) for a, b in
                             [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]]
    h = 7
    for heavy, n in [(0, 2), (2, 1), (3, 1), (4, 1), (5, 1), (6, 1)]:
        for _ in range(n):
            bonds.append((heavy, h, 1.0))
            h += 1
    d = 14  # donor offset
    bonds += [(d, d + 1, 1.0), (d, d + 2, 1.0), (d, d + 3, 1.0), (d, d + 4, 1.0)]
    edits = [AtomEdit(0, radical=0), BondEdit(0, d + 1, 1.0), BondEdit(d, d + 1, 0.0),
             AtomEdit(d, radical=1)]
    return Template("benzyl_abstraction", HA, atoms, bonds, edits, weight=1.0)


def recombination_methyl():
    atoms = [WindowAtom("C", radical=1)] + [WindowAtom("H")] * 3
    atoms += [WindowAtom("C", radical=1, participant=1)] + [WindowAtom("H", participant=1)] * 3
    bonds = [(0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0), (4, 5, 1.0), (4, 6, 1.0), (4, 7, 1.0)]
    edits = [AtomEdit(0, radical=0), AtomEdit(4, radical=0), BondEdit(0, 4, 1.0)]
    return Template("methyl_recombination", RC, atoms, bonds, edits, weight=1.0)


def h_atoms():
    atoms = [WindowAtom("H", radical=1), WindowAtom("H", radical=1, participant=1)]
    edits = [AtomEdit(0, radical=0), AtomEdit(1, radical=0), BondEdit(0, 1, 1.0)]
    return Template("h_plus_h", HH, atoms, [], edits, weight=1.0)


def heterolysis():
    """CH3-CH3 -> CH3+ + CH3-: two outcomes that differ only by charge allocation."""
    atoms = [WindowAtom("C")] + [WindowAtom("H")] * 3 + [WindowAtom("C")] + [WindowAtom("H")] * 3
    bonds = [(0, 4, 1.0), (0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0), (4, 5, 1.0), (4, 6, 1.0),
             (4, 7, 1.0)]
    edits = [BondEdit(0, 4, 0.0), AtomEdit(0, charge=1), AtomEdit(4, charge=-1, lone_pairs=1)]
    return Template("heterolysis", "heterolysis", atoms, bonds, edits, weight=1.0)


# --------------------------------------------------------------------------
# family rules and window cutting (test-side; no matcher code)
# --------------------------------------------------------------------------
def family_rule(desc):
    """(participant of each edited atom, [edits], path predicate or None) for
    an oracle move description.  Edits: ("rad", a, value) / ("bond", a, b, order)."""
    tag = desc[0]
    if tag == "HA":
        _, x, y, h = desc
        return ({x: 0, y: 1, h: 1},
                [("rad", x, 0), ("bond", x, h, 1.0), ("bond", y, h, 0.0), ("rad", y, 1)], None)
    if tag == "BS":
        _, x, a, b = desc
        return ({x: 0, a: 0, b: 0},
                [("bond", x, a, 2.0), ("bond", a, b, 0.0), ("rad", x, 0), ("rad", b, 1)], None)
    if tag == "RC":
        _, x, y = desc
        return {x: 0, y: 1}, [("rad", x, 0), ("rad", y, 0), ("bond", x, y, 1.0)], None
    if tag == "DS":
        _, x, z, y, h = desc
        return ({x: 0, z: 1, y: 1, h: 1},
                [("rad", x, 0), ("bond", x, h, 1.0), ("bond", y, h, 0.0), ("bond", z, y, 2.0),
                 ("rad", z, 0)], None)
    if tag == "HS":
        _, n, x, y, h = desc
        return ({x: 0, y: 0, h: 0},
                [("rad", x, 0), ("bond", x, h, 1.0), ("bond", y, h, 0.0), ("rad", y, 1)],
                (x, y, n - 1, n - 1))
    if tag == "HH":
        _, a, b = desc
        return {a: 0, b: 1}, [("rad", a, 0), ("rad", b, 0), ("bond", a, b, 1.0)], None
    raise AssertionError(desc)


def cut_window(state, oracle, desc):
    """Template data (atoms, bonds, edits, paths) for the window of one move."""
    parts, edits, path = family_rule(desc)
    edited = list(dict.fromkeys(a for a in desc[1:] if isinstance(a, int) and a in parts))
    window = dict(parts)
    for a in edited:
        for b in sorted(state.atoms[a].edges):
            if b in window:
                assert window[b] == parts[a]
            else:
                window[b] = parts[a]
    order = edited + sorted(a for a in window if a not in parts)
    index = {a: k for k, a in enumerate(order)}
    atoms = []
    for a in order:
        node = state.atoms[a]
        outside = sorted(b for b in node.edges if b not in window)
        closed = a in parts or not outside
        exterior = tuple((node.edges[b], state.atoms[b].element, b in oracle.in_ring)
                         for b in outside)
        atoms.append(WindowAtom(node.element, participant=window[a],
                                kind="closed" if closed else "port", exterior=exterior,
                                isotope=node.isotope, charge=node.charge, radical=node.radical,
                                lone_pairs=node.lone_pairs))
    bonds = [(index[a], index[b], state.atoms[a].edges[b])
             for a in order for b in sorted(state.atoms[a].edges) if b in index and a < b]
    out = []
    for edit in edits:
        if edit[0] == "rad":
            out.append(AtomEdit(index[edit[1]], radical=edit[2]))
        else:
            out.append(BondEdit(index[edit[1]], index[edit[2]], edit[3]))
    paths = () if path is None else (PathPredicate(index[path[0]], index[path[1]], path[2], path[3]),)
    return atoms, bonds, out, paths


def _window_graph(atoms, bonds, edits, paths, swap=False):
    g = nx.Graph()
    final = {}
    for e in edits:
        if isinstance(e, AtomEdit):
            final[e.atom] = (e.radical, e.charge, e.lone_pairs)
    for i, a in enumerate(atoms):
        part = 1 - a.participant if swap else a.participant
        g.add_node(i, label=repr((a.element, a.isotope, a.charge, a.radical, a.lone_pairs, a.kind,
                                  a.exterior, part, final.get(i))))
    for i, j, o in bonds:
        g.add_edge(i, j, label=("bond", float(o), None))
    for e in edits:
        if isinstance(e, BondEdit):
            old = g.edges[e.a, e.b]["label"][1] if g.has_edge(e.a, e.b) else 0.0
            g.add_edge(e.a, e.b, label=("bond", old, float(e.order)))
    for p in paths:
        g.add_edge(p.a, p.b, label=("path", p.minimum, p.maximum))
    return g


def _same(g1, g2):
    return nx.is_isomorphic(g1, g2, node_match=lambda x, y: x["label"] == y["label"],
                            edge_match=lambda x, y: x["label"] == y["label"])


def cut_templates(family, witnesses, n=None):
    """One class per isomorphism class of windows cut around every oracle move
    on the witness states (participants may be exchanged)."""
    classes = []  # (graph, hash, template data)
    for state in witnesses:
        oracle = Oracle(state)
        for desc in oracle.moves(family, n).values():
            data = cut_window(state, oracle, desc)
            g = _window_graph(*data)
            h = nx.weisfeiler_lehman_graph_hash(g, node_attr="label", edge_attr="label")
            two = any(a.participant == 1 for a in data[0])
            swapped = _window_graph(*data, swap=True) if two else None
            hs = (nx.weisfeiler_lehman_graph_hash(swapped, node_attr="label", edge_attr="label")
                  if two else None)
            if any((h == ch and _same(g, cg)) or (two and hs == ch and _same(swapped, cg))
                   for cg, ch, _ in classes):
                continue
            classes.append((g, h, data))
    tag = family if n is None else f"{family}{n}"
    return [Template(f"{tag}#{k}", family, *data[:3], weight=1.0, paths=data[3])
            for k, (_g, _h, data) in enumerate(classes)]


WITNESS_DPS = (6, 7, 9)
CASES = [(HA, None), (BS, None), (RC, None), (DS, None), (HS, 4), (HS, 5), (HS, 6), (HH, None)]
_TEMPLATES = {}


def templates_for(family, n=None):
    if (family, n) not in _TEMPLATES:
        witnesses = [rich_state(d) for d in WITNESS_DPS]
        _TEMPLATES[(family, n)] = cut_templates(family, witnesses, n)
    return _TEMPLATES[(family, n)]


# --------------------------------------------------------------------------
# comparison with the oracle
# --------------------------------------------------------------------------
def oracle_record(rows):
    return {a: (r[0], r[1] or 0, r[2], r[3], r[4], r[5]) for a, r in rows.items()}


def matcher_keys(state, result):
    """(family, oracle canonical product) per matcher transition (a Counter)."""
    keys = Counter()
    for tr in result.transitions:
        keys[(tr.family, canonical(oracle_record(product_rows(state, tr))))] += 1
    return keys


REPORT = []


def compare(family, n, state, label):
    templates = templates_for(family, n)
    t0 = time.perf_counter()
    result = match(templates, state)
    t1 = time.perf_counter()
    oracle = Oracle(state)
    expected = oracle.moves(family, n)
    t2 = time.perf_counter()
    got = matcher_keys(state, result)
    t3 = time.perf_counter()
    stats = result.classes
    embeddings = sum(s.embeddings for s in stats.values())
    evaluated = sum(s.evaluated for s in stats.values())
    row = (label, family if n is None else f"{family} n={n}", len(templates), embeddings,
           evaluated, len(result.transitions), oracle.last_raw, len(expected),
           round(t1 - t0, 2), round(t2 - t1, 2), round(t3 - t2, 2))
    REPORT.append(row)
    print("WMROW", *row)
    duplicates = [k for k, v in got.items() if v != 1]
    assert not duplicates, f"{len(duplicates)} oracle moves found more than once"
    assert set(got) == set(expected), (
        f"matcher-only {len(set(got) - set(expected))}, oracle-only {len(set(expected) - set(got))}")
    for name, s in stats.items():
        if s.pair_count is not None:
            assert s.pair_count == s.transitions, (name, s.pair_count, s.transitions)
    return result


# --------------------------------------------------------------------------
# template contract
# --------------------------------------------------------------------------
def test_template_refusals():
    base = methyl_acceptor_ethane_donor()
    atoms, bonds, edits = list(base.atoms), list(base.bonds), list(base.edits)
    with pytest.raises(TemplateError, match="port"):
        Template("x", HA, atoms, bonds, edits + [AtomEdit(8, radical=1)], weight=1.0)
    with pytest.raises(TemplateError, match="connectivity"):
        Template("x", HA, atoms, bonds, [AtomEdit(0, radical=0), AtomEdit(4, radical=1)],
                 weight=1.0)
    with pytest.raises(TemplateError, match="participant"):
        Template("x", HA, atoms, bonds + [(0, 8, 1.0)], edits, weight=1.0)
    with pytest.raises(TemplateError, match="weight"):
        Template("x", HA, atoms, bonds, edits, weight=0.0)
    with pytest.raises(TemplateError, match="exterior"):
        Template("x", HA, atoms[:8] + [WindowAtom("C", participant=1, kind="port")],
                 bonds, edits, weight=1.0)
    with pytest.raises(TemplateError, match="no-op"):
        Template("x", HA, atoms, bonds, edits + [AtomEdit(1, radical=0)], weight=1.0)
    with pytest.raises(TemplateError, match="valence"):
        Template("x", HA, atoms, bonds, [AtomEdit(0, radical=0), BondEdit(0, 5, 2.0),
                                         BondEdit(4, 5, 0.0), AtomEdit(4, radical=1)], weight=1.0)
    with pytest.raises(TemplateError, match="unique"):
        match([base, base], molecules(methyl(), ethane()))


def test_role_swap_is_computed_from_the_template():
    assert recombination_methyl().role_swap is True
    assert h_atoms().role_swap is True
    assert methyl_acceptor_ethane_donor().role_swap is False
    assert {t.role_swap for t in templates_for(DS)} == {False}
    assert {t.role_swap for t in templates_for(HA)} == {False}
    swaps = Counter(t.role_swap for t in templates_for(RC))
    assert swaps[True] >= 2 and swaps[False] >= 1, swaps  # (P,P), (B,B) vs (P,B)


# --------------------------------------------------------------------------
# embeddings, dedup and termini
# --------------------------------------------------------------------------
def test_ch3_abstraction_gives_three_transitions_per_methyl():
    state = molecules(methyl(), ethane())
    result = match([methyl_acceptor_ethane_donor()], state)
    stats = result.classes["ch3_abstraction"]
    # acceptor: 3! H permutations; donor: 2 carbons x 3 moved H x 2 spectator orders
    assert stats.embeddings == 6 * 12
    assert stats.transitions == 6 == len(result.transitions)
    moved = Counter(tr.mapping[4] for tr in result.transitions)
    assert sorted(moved.values()) == [3, 3]
    oracle = Oracle(state)
    assert len({canonical(oracle_record(product_rows(state, t))) for t in result.transitions}) == 6
    assert stats.sites[0] == {0: ((0,),)}
    assert stats.pair_count == 6


def test_phenyl_flips_and_spectator_h_collapse():
    state = molecules(benzyl(), ethane())
    result = match([benzyl_acceptor_ethane_donor()], state)
    stats = result.classes["benzyl_abstraction"]
    # acceptor: 2 ring flips x 2 CH2 H orders, matched in 1 of the 4 forms
    assert stats.embeddings == 4 * 12
    assert stats.transitions == 6
    assert len(stats.sites[0][0]) == 1


def test_closed_atoms_and_ports_respect_termini():
    ch3 = methyl_acceptor_ethane_donor()
    assert len(match([ch3], molecules(methyl(), ethane())).transitions) == 6
    # in propane the port's exterior is (C, H, H), not (H, H, H): no site
    assert len(match([ch3], molecules(methyl(), propane())).transitions) == 0
    # a closed carbon whose window degree is 3 never matches a 4-coordinate carbon
    atoms = [a for k, a in enumerate(ch3.atoms) if k != 7]
    bonds = [(0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0), (4, 5, 1.0), (4, 6, 1.0), (4, 7, 1.0)]
    short = Template("short", HA, atoms, bonds, ch3.edits, weight=1.0)
    assert len(match([short], molecules(methyl(), ethane())).transitions) == 0
    ethyl = with_h([("C", 1), ("C", 0)], [(0, 1, 1)], [2, 3])
    assert len(match([short], molecules(methyl(), ethyl)).transitions) == 0  # radical differs


# --------------------------------------------------------------------------
# resonance
# --------------------------------------------------------------------------
def _oracle_placements(rows):
    """Integer placements (brute force, oracle side) of every region of a component."""
    record = oracle_record(rows)
    regions, free = transition_oracle._regions(record)
    out = []
    for region in regions:
        bonds, placements = transition_oracle._placements(region, record, free)
        out.append((region, bonds, placements))
    return out


def assert_forms_equal_oracle(state, cid):
    forms = component_forms(state, cid)
    rows = [form_rows(state, cid, f) for f in forms]
    assert len(set(map(lambda r: tuple(sorted(r.items())), rows))) == len(forms)
    keys = {canonical(oracle_record(r)) for r in rows}
    assert len(keys) == 1  # every form is a valid rendering of the same species
    regions = _oracle_placements(rows[0])
    total = 1
    for region, bonds, placements in regions:
        # each integer placement is rendered by exactly one distinct region form
        rendered = {}
        for r in rows:
            rec = oracle_record(r)
            sig = (tuple(rec[a][3] for a in region),
                   tuple(dict(rec[a][5])[b] for a, b in bonds))
            rendered[sig] = rec
        for placement in placements:
            hits = [s for s, rec in rendered.items()
                    if transition_oracle._consistent(region, bonds, placement, rec)]
            assert len(hits) == 1, (region, placement, hits)
        for sig, rec in rendered.items():
            assert any(transition_oracle._consistent(region, bonds, p, rec) for p in placements)
        total *= len(rendered)
    assert total == len(forms)
    return forms


def test_benzyl_forms_equal_brute_force_and_do_not_depend_on_the_start_form():
    reference = None
    for form in ("aromatic", "ortho", "para", "ortho2"):
        state = molecules(benzyl(form))
        forms = assert_forms_equal_oracle(state, 0)
        assert len(forms) == 4
        listing = [form_rows(state, 0, f) for f in forms]
        if reference is None:
            reference = listing
        assert listing == reference
        first = listing[0]  # aromatic form first: radical on CH2, ring bonds 1.5
        assert first[0][3] == 1
        assert dict(first[1][5]) == {0: 1.0, 2: 1.5, 6: 1.5}


def test_cross_form_matching_gives_one_transition_per_outcome():
    template = benzyl_acceptor_ethane_donor()
    results = {}
    for form in ("aromatic", "ortho", "para"):
        state = molecules(benzyl(form), ethane())
        result = match([template], state)
        assert len(result.transitions) == 6
        results[form] = sorted((t.edit_key, t.identity) for t in result.transitions)
    assert results["aromatic"] == results["ortho"] == results["para"]


def _polyenyl(n_carbons):
    heavy = [("C", 1 if i == n_carbons - 1 else 0) for i in range(n_carbons)]
    bonds = [(i, i + 1, 2 if i % 2 == 0 and i + 1 < n_carbons - 1 else 1)
             for i in range(n_carbons - 1)]
    hydrogens = [2] + [1] * (n_carbons - 2) + [2]
    return with_h(heavy, bonds, hydrogens)


def test_resonance_limit_raises_at_form_65():
    state = molecules(_polyenyl(127))
    assert len(component_forms(state, 0)) == RESONANCE_FORM_LIMIT == 64
    state = molecules(_polyenyl(129))
    with pytest.raises(ResonanceLimitError) as info:
        component_forms(state, 0)
    assert info.value.count == RESONANCE_FORM_LIMIT + 1
    with pytest.raises(ResonanceLimitError):
        match([h_atoms()], state)


def test_conjugated_ring_larger_than_eight_is_refused():
    cot = with_h([("C", 0)] * 8, [(i, (i + 1) % 8, 2 if i % 2 == 0 else 1) for i in range(8)],
                 [1] * 8)
    state = molecules(cot)
    assert len(assert_forms_equal_oracle(state, 0)) == 2
    ten = with_h([("C", 0)] * 10, [(i, (i + 1) % 10, 2 if i % 2 == 0 else 1) for i in range(10)],
                 [1] * 10)
    with pytest.raises(ConjugatedRingError):
        component_forms(molecules(ten), 0)


def test_dp30_radical_component_form_counts():
    """Premise P2: DP30 PS chains with one and two benzylic radicals stay under 64."""
    state = state_at(30)
    counts = {}
    for cid in state.components.ids():
        if state.components.radical_count(cid) and len(state.components.atoms(cid)) > 1:
            counts[cid] = len(component_forms(state, cid))
    print("WMFORMS", sorted(counts.values()))
    assert sorted(counts.values()) == [1, 1, 4, 4, 16]
    two = max(counts, key=counts.get)
    assert_forms_equal_oracle(state, two)


def test_identity_distinguishes_charge_allocation():
    state = molecules(ethane())
    result = match([heterolysis()], state)
    assert result.classes["heterolysis"].embeddings == 2 * 36
    assert len(result.transitions) == 2
    charges = sorted(tuple(sorted((r[2] for r in product_rows(state, t).values() if r[2])))
                     for t in result.transitions)
    assert charges == [(-1, 1), (-1, 1)]
    assert len({t.identity for t in result.transitions}) == 2


# --------------------------------------------------------------------------
# acceptance against the oracle
# --------------------------------------------------------------------------
@pytest.mark.parametrize("family,n", CASES, ids=[f"{f}{n or ''}" for f, n in CASES])
@pytest.mark.parametrize("dp", [8, 12, 21, pytest.param(30, marks=pytest.mark.slow)])
def test_transition_sets_equal_the_oracle(family, n, dp):
    if dp == 30 and not SLOW and os.environ.get("RMG_KMC_DP30") != "1":
        pytest.skip("DP30: set RMG_KMC_SLOW=1 (or RMG_KMC_DP30=1)")
    compare(family, n, state_at(dp), f"DP{dp}")


def test_positive_control_dp3_and_the_legacy_one_per_anchor_matcher():
    """DP3: the matcher equals the oracle; the legacy SSA matcher keeps one
    match per anchor, so on a C-H window it sees one site per H-bearing carbon
    where the oracle (and the new matcher) see one per C-H bond."""
    from rmgpy.kmc.ssa import _subgraph_mappings

    state = merge(materialize(ps_spec({3: 2})))
    remove_h(state, find(state, 0, "unit 3", "c2"))
    for family in (HA, HH):
        compare(family, None, state, "DP3")
    component = {str(a): {"element": state.atoms[a].element, "radical": state.atoms[a].radical,
                          "edges": {str(b): o for b, o in state.atoms[a].edges.items()}}
                 for a in state.components.atoms(state.components.component_of(
                     find(state, 1, "unit 1", "c1")))}
    expected = {0: {"element": "C", "edges": {1: 1.0}}, 1: {"element": "H", "edges": {0: 1.0}}}
    legacy = _subgraph_mappings(expected, component)
    carbons_with_h = sum(1 for a in component.values() if a["element"] == "C"
                         and any(component[b]["element"] == "H" for b in a["edges"]))
    ch_bonds = sum(1 for a in component.values() if a["element"] == "C"
                   for b in a["edges"] if component[b]["element"] == "H")
    assert len(legacy) == carbons_with_h < ch_bonds
    # the new matcher: one HA transition per (benzylic acceptor, sp3 C-H bond)
    result = match(templates_for(HA), state)
    chain1 = set(state.components.atoms(state.components.component_of(find(state, 1, "unit 1", "c1"))))
    sp3_ch = sum(1 for a in chain1 if state.atoms[a].element == "C" and len(state.atoms[a].edges) == 4
                 for b in state.atoms[a].edges if state.atoms[b].element == "H")
    by_name = {t.name: t for t in templates_for(HA)}
    donors = [t for t in result.transitions
              if t.mapping[by_name[t.template].edited[1][0]] in chain1]
    assert len(donors) == sp3_ch


def test_identical_windows_on_dp10_and_dp20_give_equal_site_counts():
    def two_chains(a, b):
        state = materialize(ps_spec({a: 1, b: 1} if a != b else {a: 2}))
        for chain, dp in ((0, a), (1, b)):
            remove_h(state, find(state, chain, f"unit {dp}", "c2"))
        return state

    long = two_chains(10, 20)
    short = two_chains(10, 10)
    results = {}
    for label, state in (("10+20", long), ("10+10", short)):
        res = {}
        for family in (RC, HA):
            result = compare(family, None, state, label)
            for t in templates_for(family):
                s = result.classes[t.name]
                if s.transitions:
                    res[(family, t.name)] = (s.pair_count, tuple(sorted(
                        len(v) for v in s.sites[0].values())))
        results[label] = res
    rc = {k: v for k, v in results["10+20"].items() if k[0] == RC}
    assert rc == {k: v for k, v in results["10+10"].items() if k[0] == RC}
    assert list(rc.values()) == [(1, (1, 1))]
    ha_long = {k: v[1] for k, v in results["10+20"].items() if k[0] == HA}
    ha_short = {k: v[1] for k, v in results["10+10"].items() if k[0] == HA}
    assert ha_long == ha_short  # one acceptor site per chain in both states


# --------------------------------------------------------------------------
# participant quotient and the exact pair sampler
# --------------------------------------------------------------------------
class _NeedMore(Exception):
    def __init__(self, n):
        super().__init__(n)
        self.n = n


class _ScriptedRng:
    def __init__(self, prefix):
        self.prefix = prefix
        self.sizes = []

    def randrange(self, n):
        k = len(self.sizes)
        if k >= len(self.prefix):
            raise _NeedMore(n)
        self.sizes.append(n)
        return self.prefix[k]


def exact_distribution(draw):
    """Exact output distribution of ``draw(rng)`` by enumerating every
    sequence of ``rng.randrange`` outcomes."""
    dist = Counter()
    stack = [()]
    while stack:
        prefix = stack.pop()
        rng = _ScriptedRng(prefix)
        try:
            out = draw(rng)
        except _NeedMore as more:
            stack.extend(prefix + (k,) for k in range(more.n))
            continue
        p = Fraction(1)
        for n in rng.sizes:
            p /= n
        dist[out] += p
    return dist


def _sites(spec):
    return {c: tuple((c, k) for k in range(n)) for c, n in spec.items() if n}


def test_pair_sampler_exact_probabilities_ordered():
    n0 = {0: 1, 1: 0, 2: 2, 3: 1, 4: 3, 5: 0}
    n1 = {0: 2, 1: 3, 2: 0, 3: 1, 4: 2, 5: 1}
    sampler = PairSampler(_sites(n0), _sites(n1), unordered=False)
    pairs = [(s0, s1) for c0, v0 in _sites(n0).items() for s0 in v0
             for c1, v1 in _sites(n1).items() if c1 != c0 for s1 in v1]
    total0, total1 = sum(n0.values()), sum(n1.values())
    assert sampler.total == len(pairs) == total0 * total1 - sum(n0[c] * n1[c] for c in n0)
    dist = exact_distribution(sampler.draw)
    assert dist == {p: Fraction(1, len(pairs)) for p in pairs}
    sampler.update(5, 2, 1)
    n0[5] = 2
    pairs = [(s0, s1) for c0, v0 in _sites(n0).items() for s0 in v0
             for c1, v1 in _sites(n1).items() if c1 != c0 for s1 in v1]
    assert exact_distribution(sampler.draw) == {p: Fraction(1, len(pairs)) for p in pairs}


def test_pair_sampler_exact_probabilities_unordered():
    n = {0: 1, 1: 2, 2: 0, 3: 3, 4: 1, 5: 2}
    sites = _sites(n)
    sampler = PairSampler(sites, sites, unordered=True)
    pairs = {tuple(sorted((a, b))) for c0, v0 in sites.items() for a in v0
             for c1, v1 in sites.items() if c1 != c0 for b in v1}
    assert sampler.total == len(pairs)
    assert exact_distribution(sampler.draw) == {p: Fraction(1, len(pairs)) for p in pairs}


def test_pair_sampler_on_matched_sites():
    state = molecules(methyl(), methyl(), ethane(), methyl(), ethane(), propane())
    result = match([methyl_acceptor_ethane_donor(), recombination_methyl()], state)
    for name in ("ch3_abstraction", "methyl_recombination"):
        stats = result.classes[name]
        sampler = result.sampler(name)
        joint = {t.edit_key for t in result.transitions if t.template == name}
        assert sampler.total == stats.pair_count == len(joint)
        dist = exact_distribution(sampler.draw)
        assert len(dist) == stats.pair_count
        assert set(dist.values()) == {Fraction(1, stats.pair_count)}


# --------------------------------------------------------------------------
# apply
# --------------------------------------------------------------------------
def _component_rows(state, atoms):
    cids = sorted({state.components.component_of(a) for a in atoms})
    return {a: (state.atoms[a].element, state.atoms[a].isotope, state.atoms[a].charge,
                state.atoms[a].radical, state.atoms[a].lone_pairs,
                tuple(sorted((b, float(o)) for b, o in state.atoms[a].edges.items())))
            for c in cids for a in state.components.atoms(c)}


@pytest.mark.parametrize("family,n", CASES, ids=[f"{f}{n or ''}" for f, n in CASES])
def test_apply_reaches_the_predicted_and_oracle_product(family, n):
    base = state_at(8)
    result = match(templates_for(family, n), base)
    expected = Oracle(base).moves(family, n)
    assert result.transitions
    picks = result.transitions[:: max(1, len(result.transitions) // 4)][:5]
    for tr in picks:
        state = pickle.loads(pickle.dumps(base))
        predicted = product_rows(state, tr)
        revision = state.revision
        apply_transition(state, tr)
        assert state.revision == revision + 1
        rows = _component_rows(state, list(predicted))
        assert rows == predicted
        assert (family, canonical(oracle_record(rows))) in expected
        for cid in sorted({state.components.component_of(a) for a in predicted}):
            assert form_rows(state, cid, component_forms(state, cid)[0]) == {
                a: rows[a] for a in state.components.atoms(cid)}


def test_apply_from_a_quinoid_start_form():
    state = molecules(benzyl("para"), ethane())
    tr = match([benzyl_acceptor_ethane_donor()], state).transitions[0]
    predicted = product_rows(state, tr)
    apply_transition(state, tr)
    assert _component_rows(state, list(predicted)) == predicted


def test_apply_rolls_back_and_refuses_stale_and_pure_pi():
    state = state_at(8)
    result = match(templates_for(BS), state)
    tr = result.transitions[0]
    before = (state.full_hash(), state.revision)
    for stage in ("after_form", "after_edits", "after_normalize"):
        def fault(name, stage=stage):
            if name == stage:
                raise RuntimeError(stage)
        with pytest.raises(RuntimeError):
            apply_transition(state, tr, fault=fault)
        assert (state.full_hash(), state.revision) == before
    pure_pi = replace(tr, edits=tuple(e for e in tr.edits if not (
        isinstance(e, SetBond) and (e.old == 0 or e.new == 0))))
    with pytest.raises(wm.WindowMatcherError, match="connectivity"):
        apply_transition(state, pure_pi)
    assert (state.full_hash(), state.revision) == before
    apply_transition(state, tr)
    with pytest.raises(StaleTransitionError):
        apply_transition(state, result.transitions[1])


# --------------------------------------------------------------------------
# determinism
# --------------------------------------------------------------------------
def fingerprint():
    import hashlib

    state = rich_state(8)
    digest = hashlib.sha256()
    for family, n in CASES:
        result = match(templates_for(family, n), state)
        for tr in result.transitions:
            digest.update(repr((tr.template, tr.edit_key, tr.identity)).encode())
        for name, s in result.classes.items():
            digest.update(repr((name, s.embeddings, s.evaluated, s.transitions, s.pair_count,
                                s.sites)).encode())
    return digest.hexdigest()


def test_output_is_identical_across_processes_and_hash_seeds():
    code = ("import sys; sys.path.insert(0, %r); import windowMatcherTest as t; "
            "print(t.fingerprint())" % str(TEST_DIR))
    outs = set()
    for seed in ("0", "4242"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        proc = subprocess.run([sys.executable, "-c", code], cwd=ROOT, env=env,
                              capture_output=True, text=True, timeout=600)
        assert proc.returncode == 0, proc.stderr[-2000:]
        outs.add(proc.stdout.strip().splitlines()[-1])
    assert len(outs) == 1


def test_transitions_are_ordered_by_integer_edit_keys():
    state = state_at(8)
    for family, n in CASES:
        result = match(templates_for(family, n), state)
        keys = [t.edit_key for t in result.transitions]
        assert all(isinstance(k, int) for key in keys for k in key)
        assert keys == sorted(keys) and len(set(keys)) == len(keys)
