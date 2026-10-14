"""Acceptance checks for the atom-native melt state (rmgpy/kmc/melt_state.py).

The state keeps atoms as the only authority; molecules are connected
components.  Every expected count below is derived independently from the PS
repeat-unit / end-group formulas, never from the module under test.
"""

import os
import pickle
import subprocess
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import networkx as nx
import pytest

from rmgpy.kmc import melt_state
from rmgpy.kmc.melt_state import (
    STAGES,
    EditError,
    MeltState,
    MeltStateError,
    Node,
    Origin,
    SetBond,
    SetCharge,
    SetLonePairs,
    SetRadical,
    materialize,
)
from rmgpy.molecule.element import get_element
from rmgpy.polymer_spec import LinearHomopolymerSpec

ROOT = Path(__file__).resolve().parents[3]
TEST_DIR = Path(__file__).resolve().parent

# PS: repeat unit C8H8, two H end groups -> chain DP n is C(8n) H(8n+2).
STYRENE = {"C": 8, "H": 8}


def _atom(key, element="C", hydrogens=0, radical=0):
    return {"key": key, "element": element, "hydrogens": hydrogens,
            "electronic_state": {"radical_electrons": radical}}


def _ring(prefix_hydrogens_first=0):
    atoms = [_atom(f"p{i}", hydrogens=prefix_hydrogens_first if i == 1 else 1)
             for i in range(1, 7)]
    bonds = [{"atoms": [f"p{i}", f"p{i + 1}"], "order": 1.5} for i in range(1, 6)]
    bonds.append({"atoms": ["p1", "p6"], "order": 1.5})
    return atoms, bonds


def ps_mapping(counts=None, mobile=None, end_element="H"):
    ring_atoms, ring_bonds = _ring()
    repeat = {
        "atoms": [_atom("c1", hydrogens=2), _atom("c2", hydrogens=1)] + ring_atoms,
        "bonds": [{"atoms": ["c1", "c2"], "order": 1},
                  {"atoms": ["c2", "p1"], "order": 1}] + ring_bonds,
        "ports": {"head": {"atom": "c1", "bond_order": 1},
                  "tail": {"atom": "c2", "bond_order": 1}},
        "electronic_state": {"multiplicity": 1},
    }

    def cap():
        return {"atoms": [_atom("h", end_element)], "bonds": [],
                "ports": {"cap": {"atom": "h", "bond_order": 1}},
                "electronic_state": {"multiplicity": 1}}

    sty_ring, sty_ring_bonds = _ring()
    styrene = {
        "atoms": [_atom("a", hydrogens=2), _atom("b", hydrogens=1)] + sty_ring,
        "bonds": [{"atoms": ["a", "b"], "order": 2},
                  {"atoms": ["b", "p1"], "order": 1}] + sty_ring_bonds,
        "electronic_state": {"multiplicity": 1},
    }
    mapping = {
        "schema_version": 1, "label": "PS", "repeat_unit": repeat,
        "end_groups": {"head": cap(), "tail": cap()},
        "small_species": {"styrene": {"graph": styrene, "disposition": "mobile"}},
        "initial_ensemble": {"dp_counts": dict(counts or {10: 1, 20: 1, 30: 1})},
        "initial_mobile_counts": dict(mobile if mobile is not None else {"styrene": 2}),
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
    return mapping


def ps_spec(counts=None, mobile=None, end_element="H"):
    return LinearHomopolymerSpec.from_mapping(ps_mapping(counts, mobile, end_element))


@pytest.fixture(scope="module")
def spec():
    return ps_spec()


@pytest.fixture()
def state(spec):
    return materialize(spec)


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def fingerprint():
    """Process-independent digest used by the cross-process test."""
    s = materialize(ps_spec())
    rows = [(i, s.origin[i].kind, s.origin[i].ordinal, s.origin[i].segment,
             s.origin[i].key, s.origin[i].is_backbone) for i in sorted(s.origin)]
    return s.state_hash() + "|" + repr(rows)


def find(state, kind, ordinal, segment, key):
    hits = [i for i, o in state.origin.items()
            if (o.kind, o.ordinal, o.segment, o.key) == (kind, ordinal, segment, key)]
    assert len(hits) == 1, (kind, ordinal, segment, key, hits)
    return hits[0]


def chain_atom(state, chain, unit, key):
    return find(state, "chain", chain, f"unit {unit}", key)


def hydrogens(state, atom):
    return sorted(j for j in state.atoms[atom].edges if state.atoms[j].element == "H")


def comp_counter(state, cid):
    return Counter({k: v for k, v in state.components.composition(cid).items()})


def element_mass(element, isotope=None):
    return Fraction(get_element(element, isotope if isotope else -1).mass)


def independent_mass(counter):
    return sum((n * element_mass(el, iso) for (el, iso), n in counter.items()),
               Fraction(0))


def assert_canonical_ids(state):
    seen = set()
    for cid in state.components.ids():
        members = state.components.atoms(cid)
        assert list(members) == sorted(members)
        assert cid == min(members)
        assert not seen & set(members)
        seen |= set(members)
        for a in members:
            assert state.components.component_of(a) == cid
    assert seen == set(state.atoms)


def snapshot(state):
    return state.state_hash(), pickle.dumps(state)


def homolysis_edits(state, chain=1, unit=10):
    a = chain_atom(state, chain, unit, "c2")
    b = chain_atom(state, chain, unit + 1, "c1")
    return a, b, [SetBond(a, b, 1.0, 0.0), SetRadical(a, 0, 1), SetRadical(b, 0, 1)]


# --------------------------------------------------------------------------
# materialization
# --------------------------------------------------------------------------
def test_materialization_counts_match_independent_formulas(state):
    expected_atoms = sum(8 * dp + 8 * dp + 2 for dp in (10, 20, 30)) + 2 * 16
    assert len(state.atoms) == expected_atoms
    assert sorted(state.atoms) == list(range(expected_atoms))
    assert all(n.element in ("C", "H") for n in state.atoms.values())
    cids = state.components.ids()
    assert len(cids) == 3 + 2
    by_size = sorted(len(state.components.atoms(c)) for c in cids)
    assert by_size == sorted([16, 16] + [16 * dp + 2 for dp in (10, 20, 30)])
    total = Counter()
    for cid in cids:
        total += comp_counter(state, cid)
    assert total == Counter({("C", None): 8 * 60 + 16, ("H", None): 8 * 60 + 6 + 16})
    # each chain composition from the formula C(8n)H(8n+2)
    chain_comps = sorted(
        (comp_counter(state, c) for c in cids if state.origin[c].kind == "chain"),
        key=lambda c: c[("C", None)])
    for dp, comp in zip((10, 20, 30), chain_comps):
        assert comp == Counter({("C", None): 8 * dp, ("H", None): 8 * dp + 2})
    # mobile molecules are styrene
    mobiles = [comp_counter(state, c) for c in cids if state.origin[c].kind == "mobile"]
    assert mobiles == [Counter({("C", None): 8, ("H", None): 8})] * 2
    # every hydrogen is an explicit node; no implicit counters survive
    assert not hasattr(next(iter(state.atoms.values())), "implicit_hydrogens")
    for i, n in state.atoms.items():
        if n.element == "H":
            assert len(n.edges) == 1
            assert state.atoms[next(iter(n.edges))].element in ("C", "H")


def test_order_and_origin_fields(state):
    # chains by (DP, ordinal), then mobiles; ids contiguous per molecule
    first = {}
    for i in sorted(state.origin):
        o = state.origin[i]
        first.setdefault((o.kind, o.ordinal), i)
    assert [k for k, _ in sorted(first.items(), key=lambda kv: kv[1])] == [
        ("chain", 0), ("chain", 1), ("chain", 2), ("mobile", 0), ("mobile", 1)]
    sizes = [len(state.components.atoms(state.components.component_of(first[("chain", n)])))
             for n in range(3)]
    assert sizes == [16 * 10 + 2, 16 * 20 + 2, 16 * 30 + 2]
    # canonical order inside chain 0: head_end, unit 1..10, tail_end
    ids = [i for i in sorted(state.origin) if state.origin[i].ordinal == 0
           and state.origin[i].kind == "chain"]
    segments = [state.origin[i].segment for i in ids]
    collapsed = [s for k, s in enumerate(segments) if k == 0 or s != segments[k - 1]]
    assert collapsed == ["head_end"] + [f"unit {k}" for k in range(1, 11)] + ["tail_end"]
    # within a unit: declaration key order, each heavy atom's H right after it
    unit1 = [i for i in ids if state.origin[i].segment == "unit 1"]
    heavy_then_h = [(state.atoms[i].element, state.origin[i].key) for i in unit1]
    assert heavy_then_h[:4] == [("C", "c1"), ("H", "c1.H0"), ("H", "c1.H1"), ("C", "c2")]
    # backbone: c1/c2 of each unit only; H and ring atoms and end groups are not
    backbone = [i for i in ids if state.origin[i].is_backbone]
    assert len(backbone) == 2 * 10
    assert {state.origin[i].key for i in backbone} == {"c1", "c2"}
    mobile_keys = {state.origin[i].segment for i in state.origin
                   if state.origin[i].kind == "mobile"}
    assert mobile_keys == {"styrene"}
    assert not any(state.origin[i].is_backbone for i in state.origin
                   if state.origin[i].kind == "mobile")
    # origin is immutable
    with pytest.raises(Exception):
        state.origin[0].ordinal = 5
    with pytest.raises(Exception):
        state.atoms[0].radical = 1


def test_ids_identical_across_processes_and_hash_seeds():
    code = ("import sys; sys.path.insert(0, %r); import meltStateTest as m; "
            "print(m.fingerprint())" % str(TEST_DIR))
    outs = []
    for seed in ("1", "12345"):
        env = dict(os.environ, PYTHONHASHSEED=seed, PYTHONPATH=str(ROOT),
                   MPLCONFIGDIR=os.environ.get("TMPDIR", "/tmp"))
        done = subprocess.run([sys.executable, "-c", code], cwd=ROOT, env=env,
                              capture_output=True, text=True, timeout=170)
        assert done.returncode == 0, done.stderr[-2000:]
        outs.append(done.stdout.strip().splitlines()[-1])
    assert outs[0] == outs[1]
    assert outs[0] == fingerprint()


def test_materialization_graph_equals_spec_graph():
    """The explicit-H chain is the spec's own authoritative capped graph."""
    from rmgpy.polymer_shared_support import _authoritative_capped_graph

    spec = ps_spec({3: 1}, {})
    reference = _authoritative_capped_graph(spec, 3)
    s = materialize(spec)
    assert len(s.components.ids()) == 1
    graph = nx.Graph()
    for i, n in s.atoms.items():
        graph.add_node(i, element=n.element)
    for i, n in s.atoms.items():
        for j, order in n.edges.items():
            graph.add_edge(i, j, order=float(order))
    ref = nx.Graph()
    index = {a: k for k, a in enumerate(reference.atoms)}
    for a in reference.atoms:
        ref.add_node(index[a], element=a.element.symbol)
    for a in reference.atoms:
        for b, bond in a.edges.items():
            ref.add_edge(index[a], index[b], order=float(bond.order))
    assert nx.is_isomorphic(
        graph, ref,
        node_match=lambda x, y: x["element"] == y["element"],
        edge_match=lambda x, y: x["order"] == y["order"])


def test_refuses_unsupported_declarations():
    with pytest.raises(MeltStateError):
        materialize(object())
    with pytest.raises(MeltStateError):
        materialize({"not": "a spec"})
    mapping = ps_mapping()
    # a hypervalent / unsupported element is refused with a clear error
    mapping["end_groups"]["head"]["atoms"][0]["element"] = "Si"
    mapping["end_groups"]["head"]["atoms"][0]["hydrogens"] = 3
    mapping["locality_policy"]["bounds"]["elements"] = ["C", "H", "Si"]
    spec = LinearHomopolymerSpec.from_mapping(mapping)
    with pytest.raises(MeltStateError, match="Si"):
        materialize(spec)


# --------------------------------------------------------------------------
# components and observables
# --------------------------------------------------------------------------
def test_component_index_invariants_and_exact_mass(state):
    assert_canonical_ids(state)
    for cid in state.components.ids():
        comp = comp_counter(state, cid)
        assert state.components.mass(cid) == independent_mass(comp)
        assert isinstance(state.components.mass(cid), Fraction)
        assert state.components.radical_count(cid) == 0
        assert len(state.components.atoms(cid)) > 0
    with pytest.raises(KeyError):
        state.components.atoms(10 ** 9)


def test_isotope_aware_mass():
    nodes = {
        0: Node("H", 2, 0, 0, 0, {1: 1.0}),
        1: Node("H", None, 0, 0, 0, {0: 1.0}),
    }
    origin = {i: Origin("mobile", 0, "HD", f"h{i}", False) for i in nodes}
    s = MeltState(nodes, origin)
    (cid,) = s.components.ids()
    assert s.components.composition(cid) == {("H", 2): 1, ("H", None): 1}
    assert s.components.mass(cid) == element_mass("H", 2) + element_mass("H")


def test_observables(state):
    cids = state.components.ids()
    chains = lambda cid: state.origin[cid].kind == "chain"
    n, m1, m2, m3 = state.mass_moments(chains)
    masses = [state.components.mass(c) for c in cids if chains(c)]
    assert n == 3 == len(masses)
    assert m1 == sum(masses)
    assert m2 == sum(m * m for m in masses)
    assert m3 == sum(m ** 3 for m in masses)
    assert state.mass_moments(lambda cid: False) == (0, 0, 0, 0)
    everything = state.mass_moments(lambda cid: True)
    assert everything[0] == 5
    # backbone-origin counts are plain integers, one per atom flagged is_backbone
    counts = {state.origin[c].ordinal: state.backbone_count(c)
              for c in cids if chains(c)}
    assert counts == {0: 20, 1: 40, 2: 60}
    assert all(isinstance(v, int) for v in counts.values())
    assert [state.backbone_count(c) for c in cids if not chains(c)] == [0, 0]


# --------------------------------------------------------------------------
# edits
# --------------------------------------------------------------------------
def test_midchain_homolysis_splits_and_conserves(state):
    a, b, edits = homolysis_edits(state, chain=1, unit=10)
    cid = state.components.component_of(a)
    assert state.components.component_of(b) == cid
    before = comp_counter(state, cid)
    mass_before = state.components.mass(cid)
    backbone_before = state.backbone_count(cid)
    journal = state.apply(edits)
    left, right = state.components.component_of(a), state.components.component_of(b)
    assert left != right
    assert left == cid  # the head part keeps the minimum id
    assert right == min(state.components.atoms(right))
    assert comp_counter(state, left) + comp_counter(state, right) == before
    assert state.components.mass(left) + state.components.mass(right) == mass_before
    assert state.backbone_count(left) + state.backbone_count(right) == backbone_before == 40
    assert state.components.radical_count(left) == 1
    assert state.components.radical_count(right) == 1
    assert_canonical_ids(state)
    # the journal records both endpoints of the changed edge and old membership
    assert set(journal.old_nodes) == {a, b}
    assert journal.old_nodes[a].edges[b] == 1.0
    assert cid in dict(journal.old_components)
    assert dict(journal.old_components)[cid] == tuple(sorted(
        state.components.atoms(left) + state.components.atoms(right)))


def test_reforming_the_bond_restores_the_hash(state):
    h0 = state.state_hash()
    a, b, edits = homolysis_edits(state)
    cid = state.components.component_of(a)
    state.apply(edits)
    assert state.state_hash() != h0
    state.apply([SetBond(b, a, 0.0, 1.0), SetRadical(a, 1, 0), SetRadical(b, 1, 0)])
    assert state.state_hash() == h0
    assert state.components.component_of(a) == state.components.component_of(b) == cid
    assert_canonical_ids(state)


def test_join_gets_minimum_id_regardless_of_edit_order(state):
    # join the radical right fragment of chain 2 to the head of chain 0:
    # the merged component must carry the minimum id of all its atoms.
    a, b, edits = homolysis_edits(state, chain=2, unit=5)
    state.apply(edits)
    head = chain_atom(state, 0, 1, "c1")  # CH2 head carbon of the DP 10 chain (non-radical)
    hy = hydrogens(state, head)[-1]  # a CH2 hydrogen, not the end cap
    # make a radical on chain 0's unit-1 c1 by removing one H, then bond to b
    state.apply([SetBond(head, hy, 1.0, 0.0), SetRadical(head, 0, 1), SetRadical(hy, 0, 1)])
    state.apply([SetBond(b, head, 0.0, 1.0), SetRadical(b, 1, 0), SetRadical(head, 1, 0)])
    assert_canonical_ids(state)
    merged = state.components.component_of(b)
    assert merged == state.components.component_of(head) == min(state.components.atoms(merged))
    assert merged < head  # head_end atoms of chain 0 precede


def test_h_abstraction_between_chains(state):
    # radical on chain 2 via homolysis, abstract the tertiary H of chain 0 unit 5
    a, b, edits = homolysis_edits(state, chain=2, unit=5)
    state.apply(edits)
    donor = chain_atom(state, 0, 5, "c2")
    hy = hydrogens(state, donor)[0]
    cid_a = state.components.component_of(donor)
    cid_b = state.components.component_of(b)
    assert cid_a != cid_b
    ca, cb = comp_counter(state, cid_a), comp_counter(state, cid_b)
    state.apply([
        SetBond(donor, hy, 1.0, 0.0),
        SetBond(b, hy, 0.0, 1.0),
        SetRadical(donor, 0, 1),
        SetRadical(b, 1, 0),
    ])
    assert state.components.component_of(hy) == state.components.component_of(b)
    na = comp_counter(state, state.components.component_of(donor))
    nb = comp_counter(state, state.components.component_of(b))
    assert ca - na == Counter({("H", None): 1})
    assert nb - cb == Counter({("H", None): 1})
    assert state.components.radical_count(state.components.component_of(donor)) == 1
    assert state.components.radical_count(state.components.component_of(b)) == 0
    assert_canonical_ids(state)


def _release_h(state, chain, unit):
    carbon = chain_atom(state, chain, unit, "c2")
    h = hydrogens(state, carbon)[0]
    state.apply([SetBond(carbon, h, 1.0, 0.0), SetRadical(carbon, 0, 1), SetRadical(h, 0, 1)])
    return carbon, h


def test_free_hydrogen_is_its_own_component_and_h2_forms(state):
    n0 = len(state.components.ids())
    _, h1 = _release_h(state, 0, 3)
    cid = state.components.component_of(h1)
    assert cid == h1 and state.components.atoms(cid) == (h1,)
    assert comp_counter(state, cid) == Counter({("H", None): 1})
    assert state.components.radical_count(cid) == 1
    assert state.components.mass(cid) == element_mass("H")
    assert len(state.components.ids()) == n0 + 1
    _, h2 = _release_h(state, 1, 7)
    state.apply([SetBond(h1, h2, 0.0, 1.0), SetRadical(h1, 1, 0), SetRadical(h2, 1, 0)])
    cid = state.components.component_of(h1)
    assert cid == min(h1, h2) == state.components.component_of(h2)
    assert comp_counter(state, cid) == Counter({("H", None): 2})
    assert state.components.radical_count(cid) == 0
    assert len(state.components.ids()) == n0 + 1
    assert_canonical_ids(state)


def _ring_atoms(state, chain=0, unit=4):
    """Two adjacent ring CH carbons (ortho pair) and the rest of the ring in order."""
    ipso = chain_atom(state, chain, unit, "p1")
    ring = [ipso]
    while len(ring) < 6:
        nxt = [j for j in state.atoms[ring[-1]].edges
               if j not in ring and state.atoms[j].element == "C"
               and state.origin[j].segment == f"unit {unit}"
               and state.origin[j].key.startswith("p")]
        ring.append(min(nxt))
    return ring


def test_non_bridge_ring_cut_no_split_valence_gate(state):
    ring = _ring_atoms(state)
    r1, r2 = ring[1], ring[2]
    cid = state.components.component_of(r1)
    snap = snapshot(state)
    # refused: breaking an aromatic bond without fixing radicals / orders
    with pytest.raises(MeltStateError):
        state.apply([SetBond(r1, r2, 1.5, 0.0)])
    assert snapshot(state) == snap
    with pytest.raises(MeltStateError):
        state.apply([SetBond(r1, r2, 1.5, 0.0), SetRadical(r1, 0, 1), SetRadical(r2, 0, 1)])
    assert snapshot(state) == snap
    # accepted: open the ring to the Kekule triene diradical, radicals set
    path = [r2, ring[3], ring[4], ring[5], ring[0], r1]
    orders = [(path[0], path[1], 2.0), (path[1], path[2], 1.0), (path[2], path[3], 2.0),
              (path[3], path[4], 1.0)]
    # r2=ring3 ; ring3-ring4 single ; ring4=ring5 ; ring5-ipso single ; ipso=r1
    edits = [SetBond(r1, r2, 1.5, 0.0), SetRadical(r1, 0, 1), SetRadical(r2, 0, 1),
             SetBond(path[0], path[1], 1.5, 2.0), SetBond(path[1], path[2], 1.5, 1.0),
             SetBond(path[2], path[3], 1.5, 2.0), SetBond(path[3], path[4], 1.5, 1.0),
             SetBond(path[4], path[5], 1.5, 2.0)]
    before = comp_counter(state, cid)
    journal = state.apply(edits)
    assert state.components.component_of(r1) == state.components.component_of(r2) == cid
    assert comp_counter(state, cid) == before
    assert len(state.components.ids()) == 5
    # the union of the touched component was recomputed even though nothing split
    assert cid in dict(journal.old_components)
    assert dict(journal.new_components) == dict(journal.old_components)
    assert state.components.radical_count(cid) == 2
    assert_canonical_ids(state)


def test_ring_cut_combined_with_backbone_cut_still_splits(state):
    ring = _ring_atoms(state, chain=1, unit=3)
    a, b, edits = homolysis_edits(state, chain=1, unit=10)
    # open the ring on unit 3 in the same edit set as the backbone cut
    r1, r2 = ring[1], ring[2]
    path = [r2, ring[3], ring[4], ring[5], ring[0], r1]
    ring_edits = [SetBond(r1, r2, 1.5, 0.0), SetRadical(r1, 0, 1), SetRadical(r2, 0, 1)]
    for k, order in enumerate((2.0, 1.0, 2.0, 1.0, 2.0)):
        ring_edits.append(SetBond(path[k], path[k + 1], 1.5, order))
    cid = state.components.component_of(a)
    state.apply(edits + ring_edits)
    assert state.components.component_of(a) != state.components.component_of(b)
    assert state.components.component_of(r1) == state.components.component_of(r2) == cid
    assert len(state.components.ids()) == 6
    assert_canonical_ids(state)


def test_charge_and_lone_pair_edits_obey_valence(state):
    carbon = chain_atom(state, 0, 2, "c2")
    h = hydrogens(state, carbon)[0]
    # heterolysis to a carbocation + hydride: C+ (u0 p0 c+1), H- (u0 p1 c-1)
    snap = snapshot(state)
    with pytest.raises(MeltStateError):
        state.apply([SetBond(carbon, h, 1.0, 0.0), SetCharge(carbon, 0, 1)])
    assert snapshot(state) == snap
    state.apply([SetBond(carbon, h, 1.0, 0.0), SetCharge(carbon, 0, 1),
                 SetCharge(h, 0, -1), SetLonePairs(h, 0, 1)])
    assert state.components.atoms(h) == (h,)
    assert state.atoms[h].charge == -1 and state.atoms[h].lone_pairs == 1
    assert state.atoms[carbon].charge == 1
    assert_canonical_ids(state)


# --------------------------------------------------------------------------
# rollback
# --------------------------------------------------------------------------
@pytest.mark.parametrize("stage", STAGES)
def test_injected_exception_rolls_back_byte_identically(state, stage):
    assert set(STAGES) == {"after_apply", "after_components", "after_postcheck"}
    a, b, edits = homolysis_edits(state)
    # make the set richer than a single edge: also an H move on another chain
    donor = chain_atom(state, 0, 5, "c2")
    hy = hydrogens(state, donor)[0]
    target = chain_atom(state, 2, 6, "c1")
    h2 = hydrogens(state, target)[0]
    edits = edits + [SetBond(donor, hy, 1.0, 0.0), SetRadical(donor, 0, 1),
                     SetBond(target, h2, 1.0, 0.0), SetRadical(target, 0, 1),
                     SetBond(hy, h2, 0.0, 1.0)]
    snap = snapshot(state)
    seen = []

    def fault(name):
        seen.append(name)
        if name == stage:
            raise RuntimeError("injected")

    with pytest.raises(RuntimeError, match="injected"):
        state.apply(edits, fault=fault)
    assert seen[-1] == stage
    assert snapshot(state) == snap
    assert_canonical_ids(state)
    # the state is still usable and the same edits then succeed
    state.apply(edits)
    assert state.state_hash() != snap[0]


def test_failed_precondition_rolls_back(state):
    a, b, edits = homolysis_edits(state)
    snap = snapshot(state)
    bad = edits + [SetRadical(chain_atom(state, 0, 1, "c1"), 1, 0)]  # expects an old radical of 1
    with pytest.raises(EditError):
        state.apply(bad)
    assert snapshot(state) == snap
    with pytest.raises(EditError):
        state.apply([SetBond(a, b, 2.0, 0.0), SetRadical(a, 0, 1), SetRadical(b, 0, 1)])
    assert snapshot(state) == snap


def test_postcondition_failure_rolls_back(state):
    a, b, _ = homolysis_edits(state)
    snap = snapshot(state)
    with pytest.raises(MeltStateError):  # valence violated: bond broken, no radicals
        state.apply([SetBond(a, b, 1.0, 0.0)])
    assert snapshot(state) == snap
    with pytest.raises(MeltStateError):  # valence violated: radical without losing a bond
        state.apply([SetRadical(a, 0, 1)])
    assert snapshot(state) == snap


def test_journal_mutant_missing_an_endpoint_is_detected(state, monkeypatch):
    """Skipping one edge endpoint in the journal must break rollback visibly."""
    original = MeltState._journal_nodes

    def mutant(self, journal, ids):
        ids = sorted(ids)
        return original(self, journal, ids[:-1])

    a, b, edits = homolysis_edits(state)
    snap = snapshot(state)
    monkeypatch.setattr(MeltState, "_journal_nodes", mutant)

    def fault(name):
        if name == "after_apply":
            raise RuntimeError("injected")

    with pytest.raises(Exception):
        state.apply(edits, fault=fault)
    # the real implementation restores byte-identically; the mutant must not
    assert snapshot(state) != snap


# --------------------------------------------------------------------------
# refusals
# --------------------------------------------------------------------------
def test_edits_on_missing_atoms_and_create_attempts_are_refused(state):
    snap = snapshot(state)
    n = len(state.atoms)
    c = chain_atom(state, 0, 1, "c1")
    for bad in (
        [SetBond(c, n, 0.0, 1.0)],
        [SetBond(n + 5, c, 0.0, 1.0)],
        [SetRadical(n, 0, 1)],
        [SetRadical(-1, 0, 1)],
        [{"op": "create_atom", "element": "H"}],
        [("create", "H")],
        [SetBond(c, c, 0.0, 1.0)],
    ):
        with pytest.raises(EditError):
            state.apply(bad)
        assert snapshot(state) == snap
    assert len(state.atoms) == n


def test_duplicate_targets_and_bad_orders_are_refused(state):
    a, b, _ = homolysis_edits(state)
    snap = snapshot(state)
    with pytest.raises(EditError):
        state.apply([SetBond(a, b, 1.0, 0.0), SetBond(b, a, 1.0, 0.0)])
    with pytest.raises(EditError):
        state.apply([SetRadical(a, 0, 1), SetRadical(a, 0, 1)])
    with pytest.raises(EditError):
        state.apply([SetBond(a, b, 1.0, 2.5)])
    with pytest.raises(EditError):
        state.apply([SetRadical(a, 0, -1)])
    assert snapshot(state) == snap


def test_asymmetric_adjacency_is_refused_on_construction():
    ok = {0: Node("H", None, 0, 0, 0, {1: 1.0}), 1: Node("H", None, 0, 0, 0, {0: 1.0})}
    origin = {i: Origin("mobile", 0, "H2", f"h{i}", False) for i in ok}
    MeltState(ok, origin)
    one_sided = {0: Node("H", None, 0, 0, 0, {1: 1.0}), 1: Node("H", None, 0, 0, 0, {})}
    with pytest.raises(MeltStateError):
        MeltState(one_sided, origin)
    mismatched = {0: Node("H", None, 0, 0, 0, {1: 1.0}), 1: Node("H", None, 0, 0, 0, {0: 2.0})}
    with pytest.raises(MeltStateError):
        MeltState(mismatched, origin)
    dangling = {0: Node("H", None, 0, 0, 0, {7: 1.0})}
    with pytest.raises(MeltStateError):
        MeltState(dangling, {0: origin[0]})
    with pytest.raises(MeltStateError):  # valence: H with two bonds
        MeltState({0: Node("H", None, 0, 0, 0, {1: 1.0, 2: 1.0}),
                   1: Node("H", None, 0, 0, 0, {0: 1.0}),
                   2: Node("H", None, 0, 0, 0, {0: 1.0})},
                  {i: Origin("mobile", 0, "x", f"h{i}", False) for i in range(3)})
    with pytest.raises(MeltStateError):  # no empty state, atoms and origin must agree
        MeltState({}, {})
    with pytest.raises(MeltStateError):
        MeltState(ok, {0: origin[0]})


def test_state_hash_is_deterministic_and_content_sensitive(spec):
    s1, s2 = materialize(spec), materialize(spec)
    assert s1.state_hash() == s2.state_hash()
    assert pickle.dumps(s1) == pickle.dumps(s2)
    a, b, edits = homolysis_edits(s2)
    s2.apply(edits)
    assert s1.state_hash() != s2.state_hash()


# --------------------------------------------------------------------------
# state_hash sensitivity, one field at a time
# --------------------------------------------------------------------------
def _h2(isotope):
    nodes = {0: Node("H", isotope, 0, 0, 0, {1: 1.0}), 1: Node("H", None, 0, 0, 0, {0: 1.0})}
    origin = {i: Origin("mobile", 0, "H2", f"h{i}", False) for i in nodes}
    return MeltState(nodes, origin)


def test_isotope_changes_hash_composition_and_mass():
    plain, deuterated = _h2(None), _h2(2)
    assert plain.state_hash() != deuterated.state_hash()
    assert plain.components.composition(0) != deuterated.components.composition(0)
    assert plain.components.mass(0) != deuterated.components.mass(0)
    assert plain.components.mass(0) == 2 * element_mass("H")
    assert deuterated.components.mass(0) == element_mass("H", 2) + element_mass("H")
    carbon13 = materialize(ps_spec({2: 1}, {}))
    c = next(i for i, n in carbon13.atoms.items() if n.element == "C")
    base = carbon13.state_hash()
    carbon13._atoms[c] = melt_state.replace(carbon13.atoms[c], isotope=13)
    assert carbon13.state_hash() != base


@pytest.mark.parametrize("field", ["isotope", "charge", "radical", "lone_pairs", "bond_order"])
def test_state_hash_distinguishes_each_field(field):
    s = materialize(ps_spec({2: 1}, {}))
    base = s.state_hash()
    c = next(i for i, n in s.atoms.items() if n.element == "C")
    if field == "bond_order":
        j = next(iter(s.atoms[c].edges))
        for x, y in ((c, j), (j, c)):
            edges = dict(s.atoms[x].edges)
            edges[y] = 2.0
            s._atoms[x] = melt_state.replace(s.atoms[x], edges=edges)
    else:
        value = {"isotope": 13, "charge": 1, "radical": 1, "lone_pairs": 1}[field]
        s._atoms[c] = melt_state.replace(s.atoms[c], **{field: value})
    assert s.state_hash() != base


# --------------------------------------------------------------------------
# hardening: public transactions, conservation, install atomicity, identity
# --------------------------------------------------------------------------
def _h_move_edits(state):
    donor = chain_atom(state, 0, 5, "c2")
    hy = hydrogens(state, donor)[0]
    target = chain_atom(state, 2, 6, "c1")
    h2 = hydrogens(state, target)[0]
    return [SetBond(donor, hy, 1.0, 0.0), SetRadical(donor, 0, 1),
            SetBond(target, h2, 1.0, 0.0), SetRadical(target, 0, 1),
            SetBond(hy, h2, 0.0, 1.0)]


def test_multi_apply_transaction_rolls_back_on_exception(state):
    _a, _b, backbone_cut = homolysis_edits(state)
    snap = snapshot(state)
    with pytest.raises(RuntimeError, match="later stage"):
        with state.transaction() as tx:
            tx.apply(backbone_cut)
            assert state.state_hash() != snap[0]  # visible inside the transaction
            state.apply(_h_move_edits(state))      # state.apply accumulates too
            assert len(tx.journal.edits) == len(backbone_cut) + 5
            raise RuntimeError("later stage failed")
    assert snapshot(state) == snap
    assert state.revision == 0
    assert_canonical_ids(state)
    # explicit rollback and explicit commit
    tx = state.transaction()
    tx.apply(backbone_cut)
    tx.rollback()
    assert snapshot(state) == snap
    with state.transaction() as tx:
        tx.apply(backbone_cut)
    assert state.state_hash() != snap[0]


def test_failed_apply_inside_transaction_undoes_only_itself(state):
    _a, _b, backbone_cut = homolysis_edits(state)
    snap = snapshot(state)
    with state.transaction() as tx:
        tx.apply(backbone_cut)
        mid = snapshot(state)
        with pytest.raises(EditError):
            tx.apply([SetBond(0, 1, 2.0, 0.0)])
        assert snapshot(state) == mid
        tx.rollback()
    assert snapshot(state) == snap


def test_transaction_guards_and_revision(state, spec):
    a, b, edits = homolysis_edits(state)
    assert state.revision == 0
    journal = state.apply(edits)
    assert state.revision == 1
    with pytest.raises(MeltStateError):
        state.rollback(journal)  # already committed
    with pytest.raises(MeltStateError):
        journal.rollback_into(state) if hasattr(journal, "rollback_into") else state.rollback(journal)
    # one revision per commit, however many applies it holds
    with state.transaction() as tx:
        tx.apply(_h_move_edits(state))
        tx.apply([SetRadical(chain_atom(state, 0, 1, "c1"), 0, 0)])
    assert state.revision == 2
    # nested transactions are refused
    with state.transaction() as tx:
        with pytest.raises(MeltStateError):
            state.transaction()
        with pytest.raises(MeltStateError):
            state.begin()
    assert state.revision == 3
    # rollback never bumps, failed one-shot apply never bumps
    with pytest.raises(EditError):
        state.apply([SetBond(a, b, 1.0, 0.0)])
    with pytest.raises(RuntimeError):
        with state.transaction():
            raise RuntimeError("x")
    assert state.revision == 3
    # a journal belongs to the state that made it
    other = materialize(spec)
    foreign = other.transaction()
    with pytest.raises(MeltStateError):
        state.rollback(foreign.journal)
    foreign.rollback()
    # an open transaction's journal can be rolled back through the state
    snap = snapshot(state)
    tx = state.transaction()
    tx.apply(edits)
    state.rollback(tx.journal)
    assert snapshot(state) == snap
    with pytest.raises(MeltStateError):
        state.rollback(tx.journal)  # already rolled back
    # an explicit commit finalizes; the transaction is then closed
    tx = state.transaction()
    tx.apply(edits)
    tx.commit()
    assert state.revision == 4
    with pytest.raises(MeltStateError):
        tx.rollback()
    with pytest.raises(MeltStateError):
        tx.apply(edits)


def test_charge_creation_is_refused(state):
    carbon = chain_atom(state, 0, 2, "c2")
    h = hydrogens(state, carbon)[0]
    snap = snapshot(state)
    # both local valence equations hold (C+: 3 bonds, charge +1; H.: radical 1)
    with pytest.raises(PostconditionError, match="charge"):
        state.apply([SetBond(carbon, h, 1.0, 0.0), SetCharge(carbon, 0, 1), SetRadical(h, 0, 1)])
    assert snapshot(state) == snap
    # a conserving heterolysis is still accepted
    state.apply([SetBond(carbon, h, 1.0, 0.0), SetCharge(carbon, 0, 1),
                 SetCharge(h, 0, -1), SetLonePairs(h, 0, 1)])


def test_expected_radical_delta(state):
    a, b, edits = homolysis_edits(state)
    snap = snapshot(state)
    with pytest.raises(PostconditionError, match="radical"):
        state.apply(edits, expected_radical_delta=0)
    assert snapshot(state) == snap
    state.apply(edits, expected_radical_delta=2)
    snap = snapshot(state)
    with pytest.raises(PostconditionError, match="radical"):
        with state.transaction(expected_radical_delta=0) as tx:
            tx.apply(_h_move_edits(state))  # +2 radicals, expected 0
    assert snapshot(state) == snap
    with state.transaction(expected_radical_delta=2) as tx:
        tx.apply(_h_move_edits(state))


def test_install_is_atomic_with_respect_to_rollback(state):
    a, b, edits = homolysis_edits(state)
    snap = snapshot(state)

    def fault(name):
        if name == melt_state.INSTALL_STAGE:
            raise RuntimeError("inside install")

    with pytest.raises(RuntimeError, match="inside install"):
        state.apply(edits, fault=fault)
    assert snapshot(state) == snap
    assert_canonical_ids(state)
    with pytest.raises(RuntimeError, match="inside install"):
        with state.transaction() as tx:
            tx.apply(_h_move_edits(state))
            tx.apply(edits, fault=fault)
    assert snapshot(state) == snap
    assert state.revision == 0


def test_full_hash_covers_origin_and_revision_is_separate(spec):
    nodes = {0: Node("H", None, 0, 0, 0, {1: 1.0}), 1: Node("H", None, 0, 0, 0, {0: 1.0})}
    mobile = {i: Origin("mobile", 0, "H2", f"h{i}", False) for i in nodes}
    chain = {0: Origin("chain", 0, "unit 1", "c1", True), 1: Origin("chain", 0, "unit 1", "c2", False)}
    s1, s2 = MeltState(nodes, mobile), MeltState(nodes, chain)
    assert s1.state_hash() == s2.state_hash()  # graph identity only
    assert s1.full_hash() != s2.full_hash()    # cache identity includes ancestry
    flipped = {0: chain[0], 1: Origin("chain", 0, "unit 1", "c2", True)}
    assert MeltState(nodes, flipped).full_hash() != s2.full_hash()
    s = materialize(spec)
    f0 = s.full_hash()
    assert f0 == materialize(spec).full_hash()
    a, b, edits = homolysis_edits(s)
    s.apply(edits)
    assert s.full_hash() != f0
    s.apply([SetBond(a, b, 0.0, 1.0), SetRadical(a, 1, 0), SetRadical(b, 1, 0)])
    assert s.state_hash() == materialize(spec).state_hash()
    assert s.revision == 2  # the hash is back, the revision is not


def test_mass_units_are_explicit(state):
    for cid in state.components.ids():
        kg = state.components.mass_kg_per_mol(cid)
        g = state.components.molar_mass_g_mol(cid)
        assert isinstance(kg, Fraction) and isinstance(g, Fraction)
        assert g == kg * 1000
        assert kg == independent_mass(comp_counter(state, cid))
    styrene = next(c for c in state.components.ids() if state.origin[c].kind == "mobile")
    assert abs(float(state.components.molar_mass_g_mol(styrene)) - 104.15) < 0.01
    n, m1, _m2, _m3 = state.mass_moments(lambda c: True, unit="g_per_mol")
    assert m1 == sum(state.components.molar_mass_g_mol(c) for c in state.components.ids())
    with pytest.raises(MeltStateError):
        state.mass_moments(lambda c: True, unit="lb")


def test_mobile_origin_query_and_zero_ancestry(state):
    mobile = [c for c in state.components.ids() if state.origin[c].is_mobile]
    assert len(mobile) == 2
    for c in mobile:
        assert state.origin_kinds(c) == frozenset({"mobile"})
        assert state.backbone_count(c) == 0  # undefined in v1, not "not polymer-derived"
    chains = [c for c in state.components.ids() if not state.origin[c].is_mobile]
    assert all(state.origin_kinds(c) == frozenset({"chain"}) for c in chains)
    assert "undefined" in MeltState.backbone_count.__doc__
