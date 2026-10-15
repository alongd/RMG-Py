"""Tests for rmgpy/kmc/windows.py: recertified closed windows (items 1-3),
the span bound (4), chain-end contexts (5), the short-chain ladder (6) and the
spin/resonance diagnostics (7).

Fast part: stub observers on toy witnesses (no RMG database) for the fixed
point; the provider's end reconstruction; context enumeration and merging;
the span bound, L and the ladder plan; support and closure on real
certificates (the RMG database is loaded once).
Slow part (``RMG_KMC_SLOW=1``): the acceptance tables (``windows_acceptance``).
"""

import math
import os
from pathlib import Path

import pytest

from rmgpy.kmc import boundaries as b
from rmgpy.kmc import domain as dom
from rmgpy.kmc import rmg_moves as rm
from rmgpy.kmc import windows as win

SLOW = os.environ.get("RMG_KMC_SLOW") == "1"
HEAD, TAIL = (0, "head"), (0, "tail")


# --------------------------------------------------------------------------
# toy witnesses (no database): an alkane with one edited atom
# --------------------------------------------------------------------------
def toy(carbons=8, centre=3):
    """C<carbons> with both ends artificial; the product side puts a radical
    on ``centre`` (an edit; rows are not checked chemically)."""
    build = b.PROVIDERS["linear_alkane"].build(carbons)
    _mols, parts = b.witness_reactants([build])
    rows = rm.molecule_rows(_mols[0])
    product = dict(rows)
    row = product[centre]
    product[centre] = row[:3] + (1,) + row[4:]
    return b.make_witness(tuple(sorted(rows.items())), tuple(sorted(product.items())), parts,
                          label="toy")


def ball(witness, centre, radius):
    rows = witness.reactant_rows()
    seen = {centre}
    frontier = [centre]
    for _ in range(radius):
        nxt = []
        for a in frontier:
            for n, _o in rows[a][5]:
                if n not in seen:
                    seen.add(n)
                    nxt.append(n)
        frontier = nxt
    return frozenset(seen)


def stub(view_of, radius=lambda w: 1, centre=3, log=None):
    """An observer whose support is the ball of ``radius(w)`` bonds around
    ``centre`` and whose certificate view is ``view_of(w)``."""
    def observe(witness):
        if log is not None:
            log.append(dict(witness.extension_counts))
        support = win.Support(atoms=ball(witness, centre, radius(witness)),
                              reacting=frozenset({centre}))
        value = view_of(witness)
        return win.Observation(support=support, viewer=lambda inside, v=value: v)
    return observe


# --------------------------------------------------------------------------
# item 3: recertified fixed point
# --------------------------------------------------------------------------
def test_fixed_point_pins():
    assert win.LOOP_STATUSES == ("fixed_point", "no_fixed_point", "ceiling_exceeded",
                                 "certificate_unstable")
    assert isinstance(win.PORT_DEPTH_LIMIT, int) and win.PORT_DEPTH_LIMIT >= 2
    assert isinstance(win.FIXED_POINT_RETRIES, int) and win.FIXED_POINT_RETRIES >= 1
    res = win.fixed_point(toy(), stub(lambda w: "same"))
    assert res.status == "fixed_point", res.reason
    # ports in vector order (participant, port id); no extension was needed
    assert [k for k, _d in res.depths] == list(b.vector_ports(toy()))
    assert dict(res.depths) == {HEAD: 0, TAIL: 0}
    # one independent +1 check per artificial port, both passed
    assert sorted(c["port"] for c in res.checks[-1]) == sorted([HEAD, TAIL])
    assert all(c["same_window"] and c["same_certificate"] for c in res.checks[-1])
    assert res.window.atoms >= ball(toy(), 3, 1)


def test_fixed_point_checks_every_port_independently():
    """The certificate changes when the TAIL is extended to depth 1 only: a
    single +1 check (or a joint extension of every port) accepts too early;
    the per-port independent check extends the tail to depth 2."""
    res = win.fixed_point(toy(), stub(lambda w: "s1" if w.depth(TAIL) == 1 else "s0"))
    assert res.status == "fixed_point", res.reason
    assert dict(res.depths) == {HEAD: 0, TAIL: 2}
    failed = [c for round_ in res.checks for c in round_ if not c["same_certificate"]]
    assert {c["port"] for c in failed} == {TAIL}


def test_changing_certificate_is_not_a_fixed_point_until_it_stabilizes():
    seen = []
    res = win.fixed_point(toy(), stub(lambda w: min(w.depth(HEAD) + w.depth(TAIL), 4), log=seen))
    assert res.status == "fixed_point", res.reason
    assert dict(res.depths) == {HEAD: 2, TAIL: 2}
    # every extension was re-observed (recertified), the check witnesses too
    assert {HEAD: 2, TAIL: 2} in seen and {HEAD: 3, TAIL: 2} in seen
    assert res.retries == 2


def test_never_stabilizing_certificate_is_certificate_unstable():
    res = win.fixed_point(toy(), stub(lambda w: w.depth(HEAD) + w.depth(TAIL)))
    assert res.status == "certificate_unstable", (res.status, res.reason)
    assert res.retries == win.FIXED_POINT_RETRIES + 1 or res.retries == win.FIXED_POINT_RETRIES


def test_ceiling_refuses_instead_of_truncating():
    """The support reaches the head port: it must be extended; the ceiling
    forbids it, so the witness is refused, never truncated."""
    w = toy(centre=1)
    res = win.fixed_point(w, stub(lambda w_: "s", radius=lambda w_: 2, centre=1),
                          ceiling=w.heavy_count())
    assert res.status == "ceiling_exceeded", (res.status, res.reason)
    assert res.window is None
    # with room it reaches a fixed point after extending the head
    ok = win.fixed_point(w, stub(lambda w_: "s", radius=lambda w_: 2, centre=1))
    assert ok.status == "fixed_point" and dict(ok.depths)[HEAD] >= 1


def test_growing_window_is_no_fixed_point():
    res = win.fixed_point(toy(centre=1),
                          stub(lambda w: "s", radius=lambda w: 2 + w.depth(HEAD), centre=1))
    assert res.status == "no_fixed_point", (res.status, res.reason)
    assert "limit" in res.reason


# --------------------------------------------------------------------------
# item 5: physical contexts (provider reconstruction), enumeration, merging
# --------------------------------------------------------------------------
def ps_witness(units=3, **kw):
    build = b.PROVIDERS["polystyrene"].build("ps", b.PsChainSpec(units, **kw))
    _mols, parts = b.witness_reactants([build])
    rows = tuple(sorted(rm.molecule_rows(_mols[0]).items()))
    return b.make_witness(rows, rows, parts, label="ps")


def _partner(witness, key):
    """The backbone neighbour of a port atom in its own unit."""
    meta = witness.participant(key[0]).metadata
    port = witness.port(key)
    unit, position = meta.coordinate(port.atom_id)
    other = 1 - position
    return next(a for a, u, k in meta.coordinates if (u, k) == (unit, other))


def _h_count(rows, atom):
    return sum(1 for n, _o in rows[atom][5] if rows[n][0] == "H")


@pytest.mark.parametrize("key,end,expect", [
    (HEAD, "methyl", "methyl"), (HEAD, "primary_radical", "primary_radical"),
    (HEAD, "vinylidene", "vinylidene"), (TAIL, "benzyl", "benzyl"),
    (TAIL, "benzylic_radical", "benzylic_radical"), (TAIL, "styryl", "styryl")])
def test_physical_end_is_reconstructed_not_relabelled(key, end, expect):
    w = ps_witness()
    port = w.port(key)
    (cap,) = port.cap_atom_ids
    new, record = b.reconstruct_end(w, key, end)
    for rows in (new.reactant_rows(), new.product_rows()):
        assert cap not in rows
    built = new.port(key)
    assert built.physicality == "physical" and built.cap_atom_ids == ()
    assert built.atom_id == port.atom_id
    rows = new.reactant_rows()
    partner = _partner(w, key)
    order = dict(rows[port.atom_id][5]).get(partner)
    if end in ("vinylidene", "styryl"):
        assert order == 2.0
        assert _h_count(rows, partner) == _h_count(w.reactant_rows(), partner) - 1
        assert record.deleted and set(record.deleted) - {cap}
    else:
        assert order == 1.0
    if end in ("methyl", "benzyl"):
        fresh = [a for a in record.fresh_ids]
        assert fresh and all(a not in w.reactant_rows() for a in fresh)
        assert any(n in fresh for n, _o in rows[port.atom_id][5])
    if end.endswith("radical"):
        assert rows[port.atom_id][3] == 1
    verdict = dom.check_rows(dict(rows))
    assert verdict.accepted, verdict.reason
    assert expect in verdict.details["ends"]


def test_end_reconstruction_never_deletes_a_reacting_atom():
    w = toy(centre=0)  # the head port atom itself is edited
    with pytest.raises(b.ExtensionRefused):
        b.reconstruct_end(w, HEAD, "methyl")


def test_context_assignments_policy():
    w = ps_witness()
    alone = win.context_assignments(w, near_ports=())
    assert alone[0] == ((HEAD, None), (TAIL, None))
    singles = [a for a in alone if sum(e is not None for _k, e in a) == 1]
    assert len(singles) == 6
    # radical end chemistries interact through the radical inventory: their
    # combinations are always enumerated
    assert ((HEAD, "primary_radical"), (TAIL, "benzylic_radical")) in alone
    assert len(alone) == 8
    near = win.context_assignments(w, near_ports=(HEAD, TAIL))
    assert len(near) == 16 and len(set(near)) == 16
    physical_head = ps_witness(radical_end=True)
    only_tail = win.context_assignments(physical_head, near_ports=(TAIL,))
    assert {k for a in only_tail for k, _e in a} == {TAIL}


def _ctx(witness, atoms, view):
    window = win.Window(atoms=frozenset(atoms), support=win.Support(atoms=frozenset(atoms)))
    obs = win.Observation(support=window.support, viewer=lambda inside, v=view: v)
    return win.ContextResult(assignment=(), status="fixed_point", witness=witness, window=window,
                             observation=obs)


def test_contexts_with_equal_windows_but_unequal_certificates_do_not_merge():
    w = toy()
    atoms = ball(w, 3, 1)
    same = win.merge_contexts([_ctx(w, atoms, "a"), _ctx(w, atoms, "a")])
    assert same == ((0, 1),)
    different = win.merge_contexts([_ctx(w, atoms, "a"), _ctx(w, atoms, "b")])
    assert different == ()
    # unequal windows never merge, even with equal certificates
    assert win.merge_contexts([_ctx(w, atoms, "a"), _ctx(w, ball(w, 3, 2), "a")]) == ()
    # refusals never merge
    refused = win.ContextResult(assignment=(), status="refused:domain")
    assert win.merge_contexts([refused, refused]) == ()


# --------------------------------------------------------------------------
# items 4 and 6: span bound, L, ladder
# --------------------------------------------------------------------------
def test_span_bound_formula():
    c = win.FamilyConstants("X", radius=4, reacting_distance=1, free_label_search=True)
    bound = win.span_bound(c)
    assert bound.reach == 5
    assert bound.backbone_atoms == 1 + 2 * 5 + 1 + (2 * 5 + 1)
    assert bound.repeats == math.ceil((bound.backbone_atoms + 1) / 2)
    assert bound.L == (bound.backbone_atoms + 1) // 2
    anchored = win.span_bound(win.FamilyConstants("Y", radius=2, reacting_distance=0,
                                                  free_label_search=False))
    assert anchored.backbone_atoms == 0 + 2 * 3 + 1 and anchored.L == 4
    assert "reach" in bound.formula and "L" in bound.formula


@pytest.mark.parametrize("radius,distance,free", [(2, 1, False), (4, 0, False), (4, 1, True),
                                                  (5, 2, False)])
def test_L_and_L_plus_one(radius, distance, free):
    bound = win.span_bound(win.FamilyConstants("Z", radius, distance, free))
    L, A = bound.L, bound.backbone_atoms
    phases = [dict(head_shift=h, tail_shift=t) for h in (False, True) for t in (False, True)]
    fits = [win.chain_backbone_atoms(b.PsChainSpec(L, **p)) for p in phases]
    assert min(fits) <= A  # some phase of DP L fits in one window component
    over = [win.chain_backbone_atoms(b.PsChainSpec(L + 1, **p)) for p in phases]
    assert min(over) > A   # no phase of DP L+1 does


def test_ladder_covers_cut_phases_and_controls():
    c = win.FamilyConstants("R_Recombination", radius=2, reacting_distance=1,
                            free_label_search=False, bimolecular=True)
    bound = win.span_bound(c)
    plan = win.ladder_plan(c, bound)
    items = plan.items
    assert max(i.dp for i in items) == bound.L + 1
    assert all(i.control for i in items if i.dp == bound.L + 1)
    assert not any(i.control for i in items if i.dp <= bound.L)
    kinds = {i.kind for i in items}
    for kind in ("intact", "scission:phase_A", "scission:phase_B", "pendant_loss",
                 "ch_cleavage:benzylic", "ch_cleavage:secondary", "pair:short_short",
                 "pair:short_long"):
        assert kind in kinds, kind
    specs = [s for i in items for s in i.specs]
    assert any(s.head_shift for s in specs) and any(s.tail_shift for s in specs)
    # every declared end chemistry at every DP 1..L+1 of the intact chains
    intact = {(i.dp, i.specs[0].head_end, i.specs[0].tail_end) for i in items if i.kind == "intact"}
    for dp in range(1, bound.L + 2):
        for head in ("methyl", "primary_radical", "vinylidene"):
            for tail in ("benzyl", "benzylic_radical", "styryl"):
                assert (dp, head, tail) in intact or (dp, head, tail) in plan.excluded_keys
    # every planned species is in the declared domain (the rest is excluded
    # with its domain reason, never silently dropped)
    for item in items:
        for spec in item.specs:
            if item.kind == "pair:short_long" and spec is item.specs[-1]:
                continue
            build = b.PROVIDERS["polystyrene"].build("x", spec)
            assert dom.check(build.molecule).accepted, (item.label, dom.check(build.molecule).reason)
    assert all(reason for _key, reason in plan.excluded)
