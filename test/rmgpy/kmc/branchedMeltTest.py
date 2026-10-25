"""The kMC initial state and the declared linearity rule of a branched
declaration (S6).

``materialize`` builds the melt with ``random_bond_v1`` (below the gel
point, which the declaration checks), realizes the declared link counts
exactly and reports the branch points against their expectation; the
generator refuses a short cycle by drawing again.  At runtime a formed bond
between two residues is allowed only at two vacant sites of a declared link
type, within each residue's site count, closing no cycle shorter than
``minimum_cycle_bonds``; every other inter-residue bond is leak."""

from collections import deque
from fractions import Fraction

import numpy as np
import pytest

import melt_engine_support as S
from branched_support import branched_spec, counts
from meltEngineTest import _pair_counts, _radical_window, _recombination

from rmgpy.kmc import melt_engine as me


def _residue_of(state, a):
    from rmgpy.kmc.melt_state import segment_residue
    o = state.origin[a]
    return (o.ordinal, o.segment) if segment_residue(o.segment) else None


def _inter_residue_bonds(state):
    out = []
    for a, node in state.atoms.items():
        for b in node.edges:
            if a < b and node.element != "H" and state.atoms[b].element != "H" \
                    and _residue_of(state, a) and _residue_of(state, b) \
                    and _residue_of(state, a) != _residue_of(state, b):
                out.append((a, b))
    return out


def _cycle_through(state, a, b):
    """Bonds of the shortest cycle through the bond a-b (None: none)."""
    dist, queue = {a: 0}, deque([a])
    while queue:
        x = queue.popleft()
        for y in state.atoms[x].edges:
            if (x, y) in ((a, b), (b, a)) or y in dist or state.atoms[y].element == "H":
                continue
            dist[y] = dist[x] + 1
            queue.append(y)
    return dist[b] + 1 if b in dist else None


# --------------------------------------------------------------------------
# random_bond_v1
# --------------------------------------------------------------------------
def test_materialize_realizes_the_declared_links_exactly():
    from rmgpy.kmc.melt_state import materialize, segment_residue
    spec = branched_spec()
    state = materialize(spec)
    site_class = {(r, site.atom): site.site_class for r, residue in spec.residues for _n, site in residue.sites}
    links = {}
    for a, b in _inter_residue_bonds(state):
        pair = tuple(sorted(site_class[(segment_residue(state.origin[x].segment), state.origin[x].key)]
                            for x in (a, b)))
        links[pair] = links.get(pair, 0) + 1
    assert links == {("methylene", "ortho"): 30, ("methylene", "para"): 10}


def test_materialize_is_deterministic():
    from rmgpy.kmc.melt_state import materialize
    one, two = materialize(branched_spec()), materialize(branched_spec())
    assert one.atoms == two.atoms and one.origin == two.origin
    other = materialize(branched_spec(counts(seed=8)))
    assert other.atoms != one.atoms


def test_the_realized_melt_is_reported_against_the_declaration():
    """Branch points (a ring with all three sites linked) against their exact
    expectation under uniform site occupancy: 40 rings x P(both ortho of 30
    in 80) x P(para of 10 in 40)."""
    from rmgpy.kmc import melt_branched as mb
    from rmgpy.kmc.melt_state import materialize
    spec = branched_spec()
    report = mb.realized(spec)
    assert report["link_counts"] == report["declared_link_counts"] == {"m_o": 30, "m_p": 10}
    assert report["expected_branch_points"] == 40 * Fraction(30 * 29, 80 * 79) * Fraction(10, 40)
    assert isinstance(report["branch_points"], int)
    assert sum(size * n for size, n in report["molecule_sizes"].items()) == 60
    assert sum(report["molecule_sizes"].values()) == len(materialize(spec).components)
    assert report["gel_branching_ratio"] < 1


def test_no_cycle_in_the_melt_is_shorter_than_declared():
    from rmgpy.kmc.melt_state import materialize
    for seed in range(4):
        state = materialize(branched_spec(counts(seed=seed)))
        for a, b in _inter_residue_bonds(state):
            cycle = _cycle_through(state, a, b)
            assert cycle is None or cycle >= 8


def test_the_generator_refuses_a_short_cycle_and_draws_again():
    """One bridge, four rings: the bridge's ortho and para draws land on one
    ring now and then (a four-bond cycle through it, and two links between
    one residue pair); each is drawn again, so the melt never holds one."""
    from rmgpy.kmc import melt_branched as mb
    refused = 0
    for seed in range(12):
        spec = branched_spec(counts(residues={"Ph": 4, "M": 1}, links={"m_o": 1, "m_p": 1}, seed=seed))
        result = mb.random_bond(spec)
        refused += result.refused
        assert len({(i, j) for i, _s, j, _t in result.links}) == 2
    assert refused > 0


def test_a_declared_longer_minimum_cycle_shapes_the_melt():
    """Three rings, two bridges at ortho sites: both bridges on one ring pair
    close an 8-bond cycle, refused at a declared minimum of 10, so the melt
    is the open three-ring chain on every seed."""
    from rmgpy.kmc import melt_branched as mb

    def longer(mapping):
        mapping["minimum_cycle_bonds"] = 10
        counts(residues={"Ph": 3, "M": 2}, links={"m_o": 4, "m_p": 0})(mapping)
    for seed in range(8):
        def mutate(mapping, seed=seed):
            longer(mapping)
            mapping["initial_ensemble"]["seed"] = seed
        result = mb.random_bond(branched_spec(mutate))
        assert sorted(len(m) for m in result.molecules) == [5]


# --------------------------------------------------------------------------
# the declared linearity rule
# --------------------------------------------------------------------------
def _met(mapping):
    mapping["transport_model"] = {"met": S.pe_met()}


def _engine(spec, state, templates):
    return me.MeltEngine(spec, S.ToyCatalog(builder=lambda _t: list(templates)), 700.0,
                         np.random.default_rng(3), state=state)


def _melt(mutate=None):
    """Ph1-M1-Ph2 (o2-o2), Ph3-M2-Ph4 (o2-o2) and Ph5-M3-Ph6-M4-Ph7 (o2-o2,
    o6-o2); then M2-Ph4 and M4-Ph7 broken: two benzylic radicals with a vacant
    site and two aryl radicals at a vacant ortho site."""
    from rmgpy.kmc import melt_branched as mb
    from rmgpy.kmc.melt_state import SetBond, SetRadical

    def both(mapping):
        _met(mapping)
        if mutate is not None:
            mutate(mapping)
    spec = branched_spec(both)
    state = mb.build_molecules(spec, [
        (["Ph", "M", "Ph"], [(1, "a", 0, "o2"), (1, "b", 2, "o2")]),
        (["Ph", "M", "Ph"], [(1, "a", 0, "o2"), (1, "b", 2, "o2")]),
        (["Ph", "M", "Ph", "M", "Ph"], [(1, "a", 0, "o2"), (1, "b", 2, "o2"),
                                        (3, "a", 2, "o6"), (3, "b", 4, "o2")])])
    atom = mb.residue_atoms(state)

    def cut(molecule, bridge, ring, key):
        m, r = atom[(molecule, bridge)]["m"], atom[(molecule, ring)][key]
        state.apply([SetBond(m, r, 1.0, 0.0), SetRadical(m, 0, 1), SetRadical(r, 0, 1)])
    cut(1, 1, 2, "r2")
    cut(2, 3, 4, "r2")
    return spec, state, atom


def test_a_bond_at_two_vacant_sites_of_a_declared_link_is_allowed():
    """Bridge radical + a vacant para site on another molecule: methylene-para."""
    spec, state, atom = _melt()
    ring = atom[(0, 0)]["r4"]
    S.strip_h(state, ring)
    engine = _engine(spec, state, [_recombination(state, ring, atom[(1, 1)]["m"], "m_p")])
    allowed, leak = _pair_counts(engine, ring, atom[(1, 1)]["m"], "m_p")
    assert allowed > 0 and leak == 0


def test_a_bond_at_an_undeclared_site_is_leak():
    """Bridge radical + the ring carbon r3, which declares no site."""
    spec, state, atom = _melt()
    ring = atom[(0, 2)]["r3"]
    S.strip_h(state, ring)
    engine = _engine(spec, state, [_recombination(state, ring, atom[(1, 1)]["m"], "off_site")])
    assert _pair_counts(engine, ring, atom[(1, 1)]["m"], "off_site") == (0, 1)


def test_a_bond_of_an_undeclared_link_type_is_leak():
    """Two vacant ring sites (para + ortho) joined directly: no ring-ring link."""
    spec, state, atom = _melt()
    ring = atom[(0, 0)]["r4"]
    S.strip_h(state, ring)
    other = atom[(1, 2)]["r2"]
    engine = _engine(spec, state, [_recombination(state, ring, other, "ring_ring")])
    assert _pair_counts(engine, ring, other, "ring_ring") == (0, 1)


def test_a_bond_beyond_a_residues_site_count_is_leak():
    """A bridge carbon whose two sites are linked, made radical by an H loss,
    joined to a vacant para site: a third link on a two-site residue."""
    spec, state, atom = _melt()
    bridge, ring = atom[(0, 1)]["m"], atom[(1, 0)]["r4"]
    for a in (bridge, ring):
        S.strip_h(state, a)
    engine = _engine(spec, state, [_recombination(state, ring, bridge, "third_link")])
    assert _pair_counts(engine, ring, bridge, "third_link") == (0, 1)


def _closure(state, a, b, name):
    """A ring closure joining the radicals ``a`` and ``b`` of one molecule:
    one window (both radicals and their hydrogens closed, every other heavy
    neighbour a port), induced over all its atoms, since the two radicals'
    neighbourhoods may touch (the matcher embeds induced windows and reads a
    port's exterior outside the whole image)."""
    from rmgpy.kmc import window_matcher as wm
    from rmgpy.kmc.window_matcher import AtomEdit, BondEdit, Template, WindowAtom
    from meltEngineTest import _ring_members
    cid = state.components.component_of(a)
    rows = next(r for r in (wm.form_rows(state, cid, f) for f in wm.component_forms(state, cid))
                if r[a][3] == 1 and r[b][3] == 1)
    ring = _ring_members(state)
    order = [a, b]
    for centre in (a, b):
        order += [x for x, _o in rows[centre][5] if x not in order]
    index = {x: i for i, x in enumerate(order)}
    atoms = []
    for x in order:
        element, _iso, charge, radical, lone_pairs, edges = rows[x]
        if x in (a, b) or element == "H":
            atoms.append(WindowAtom(element, 0, radical=radical, charge=charge, lone_pairs=lone_pairs))
        else:
            exterior = tuple((o, rows[y][0], y in ring) for y, o in edges if y not in index)
            atoms.append(WindowAtom(element, 0, "port", exterior, charge=charge, radical=radical,
                                    lone_pairs=lone_pairs))
    bonds = sorted({(min(index[x], index[y]), max(index[x], index[y]), o)
                    for x in order for y, o in rows[x][5] if y in index})
    return Template(name, "toy_recombination", atoms, bonds,
                    [BondEdit(0, 1, 1), AtomEdit(0, radical=0), AtomEdit(1, radical=0)], weight=1.0e-2)


def _uni(engine, name):
    index = [t.name for t in engine.templates].index(name)
    assert engine.uni_allowed[index] + engine.uni_leak[index] > 0, f"{name} did not match"
    return engine.uni_allowed[index], engine.uni_leak[index]


def test_a_ring_closure_shorter_than_declared_is_leak():
    """The bridge radical of Ph3-M2 closing on Ph3's other vacant ortho site:
    a 4-bond cycle at two vacant sites of a declared link, within both site
    counts: leak by the cycle rule alone."""
    spec, state, atom = _melt()
    ring = atom[(1, 0)]["r6"]
    S.strip_h(state, ring)
    engine = _engine(spec, state, [_closure(state, ring, atom[(1, 1)]["m"], "short")])
    allowed, leak = _uni(engine, "short")
    assert allowed == 0 and leak > 0


def test_a_ring_closure_of_the_declared_length_is_allowed_and_a_longer_rule_refuses_it():
    """M4 radical of Ph5-M3-Ph6-M4 closing on Ph5's vacant para site: an
    8-bond cycle, allowed at the declared minimum 8, leak at 10."""
    for minimum, expect_leak in ((8, False), (10, True)):
        def longer(mapping, minimum=minimum):
            mapping["minimum_cycle_bonds"] = minimum
        spec, state, atom = _melt(longer)
        ring = atom[(2, 0)]["r4"]
        S.strip_h(state, ring)
        engine = _engine(spec, state, [_closure(state, ring, atom[(2, 3)]["m"], "closure")])
        allowed, leak = _uni(engine, "closure")
        assert (allowed == 0 and leak > 0) if expect_leak else (allowed > 0 and leak == 0)
