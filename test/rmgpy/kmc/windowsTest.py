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
from collections import defaultdict
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
    # every round: +1 and +2 for every artificial port alone and for all
    # ports jointly, all passed
    final = res.checks[-1]
    assert sorted((c["ports"], c["step"]) for c in final if not c["joint"]) == sorted(
        ((p,), s) for p in (HEAD, TAIL) for s in (1, 2))
    assert sorted((c["ports"], c["step"]) for c in final if c["joint"]) == [
        ((HEAD, TAIL), 1), ((HEAD, TAIL), 2)]
    assert all(c["same_window"] and c["same_certificate"] for c in final)
    assert res.window.atoms >= ball(toy(), 3, 1)


def test_fixed_point_checks_every_port_independently():
    """The certificate changes when the TAIL is extended to depth 1 only: a
    single +1 check (or a joint extension of every port) accepts too early;
    the per-port independent check extends the tail to depth 2."""
    res = win.fixed_point(toy(), stub(lambda w: "s1" if w.depth(TAIL) == 1 else "s0"))
    assert res.status == "fixed_point", res.reason
    assert dict(res.depths) == {HEAD: 0, TAIL: 2}
    failed = [c for round_ in res.checks for c in round_ if not c["same_certificate"]]
    assert {c["ports"] for c in failed if not c["joint"]} == {(TAIL,)}


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


def test_fixed_point_checks_plus_two_for_every_port():
    """The certificate changes once the tail is two repeats longer: a +1-only
    check accepts (0, 0); the +2 check extends the tail until +1 and +2 agree."""
    res = win.fixed_point(toy(), stub(lambda w: "s1" if w.depth(TAIL) >= 2 else "s0"))
    assert res.status == "fixed_point", res.reason
    assert dict(res.depths) == {HEAD: 0, TAIL: 2}
    first = res.checks[0]
    assert any(c["ports"] == (TAIL,) and c["step"] == 2 and not c["same_certificate"]
               for c in first)


def test_fixed_point_checks_all_ports_jointly_at_plus_one_and_plus_two():
    """Certificates that change only when every port is extended: the
    independent checks pass, the joint +1 (resp. +2) check does not."""
    both = win.fixed_point(toy(), stub(lambda w: "b" if w.depth(HEAD) and w.depth(TAIL)
                                       else "a"))
    assert both.status == "fixed_point", both.reason
    assert dict(both.depths) == {HEAD: 1, TAIL: 1}
    far = win.fixed_point(toy(), stub(lambda w: "b" if w.depth(HEAD) >= 2 and
                                      w.depth(TAIL) >= 2 else "a"))
    assert far.status == "fixed_point", far.reason
    assert dict(far.depths) == {HEAD: 2, TAIL: 2}
    assert any(c["joint"] and c["step"] == 2 and not c["same_certificate"]
               for c in far.checks[0])


def test_fixed_point_keeps_the_independent_checks_beside_the_joint_one():
    """A change seen only when the head alone is one repeat longer: the joint
    checks agree, only the independent head +1 check sees it."""
    res = win.fixed_point(toy(), stub(lambda w: "x" if (w.depth(HEAD), w.depth(TAIL)) == (1, 0)
                                      else "s"))
    assert res.status == "fixed_point", res.reason
    # (1, 0) differs from both of its +1 neighbours: the head and the tail
    # grow once more
    assert dict(res.depths) == {HEAD: 2, TAIL: 1}
    assert not res.checks[0][0]["same_certificate"] and res.checks[0][0]["ports"] == (HEAD,)


def test_fixed_point_has_no_switch_for_its_checks():
    import inspect
    params = set(inspect.signature(win.fixed_point).parameters)
    assert not params & {"check_ports", "joint_check", "steps", "checks"}, params
    assert win.CHECK_STEPS == (1, 2)


def test_reacting_atom_on_an_artificial_port_is_named():
    """The window reaches the head port, whose atom is edited: the port cannot
    be extended; the refusal names the reacting port (a cut-point event)."""
    res = win.fixed_point(toy(centre=0), stub(lambda w: "s", radius=lambda w: 1, centre=0))
    assert res.status == "refused:reacting_port", (res.status, res.reason)
    assert res.window is None


# --------------------------------------------------------------------------
# item 5: physical contexts (provider reconstruction), enumeration, merging
# --------------------------------------------------------------------------
def ps_witness(units=3, **kw):
    build = b.PROVIDERS["polystyrene"].build("ps", b.PsChainSpec(units, **kw))
    _mols, parts = b.witness_reactants([build])
    rows = tuple(sorted(rm.molecule_rows(_mols[0]).items()))
    return b.make_witness(rows, rows, parts, label="ps")


def ps_pair():
    builds = [b.PROVIDERS["polystyrene"].build("pe", b.PsChainSpec(1, tail_end="benzylic_radical")),
              b.PROVIDERS["polystyrene"].build("dp3", b.PsChainSpec(3))]
    _mols, parts = b.witness_reactants(builds)
    rows = tuple(sorted(rm.molecule_rows(_mols[0]).items()) + sorted(
        rm.molecule_rows(_mols[1]).items()))
    return b.make_witness(rows, rows, parts, label="pair")


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


def test_end_reconstruction_never_touches_a_reacting_atom():
    w = toy(centre=0)  # the head port atom itself is edited
    with pytest.raises(b.EndConflict):
        b.reconstruct_end(w, HEAD, "primary_radical")
    # a methyl end only replaces the (unreacting) cap by a fresh end H
    built, record = b.reconstruct_end(w, HEAD, "methyl")
    assert b.edits(built) == b.edits(w) and record.fresh_ids
    # a reacting cap is never deleted
    capped = toy(centre=2)
    (cap,) = capped.port(HEAD).cap_atom_ids
    product = dict(capped.product_rows())
    product[cap] = product[cap][:3] + (1,) + product[cap][4:]
    edited = b.make_witness(capped.reactants, tuple(sorted(product.items())),
                            capped.participants)
    with pytest.raises(b.EndConflict):
        b.reconstruct_end(edited, HEAD, "methyl")


def test_context_assignments_policy():
    """Near ports: every declared end and the artificial cut; far ports: the
    artificial cut and their radical ends (read whole-species); the full
    product over all ports, across participants, whatever the certificate
    reads; then one single per far non-radical end (the reduction check)."""
    import inspect
    assert "cross" not in inspect.signature(win.context_assignments).parameters
    w = ps_witness()
    alone = win.context_assignments(w, near_ports=())
    assert alone[0] == ((HEAD, None), (TAIL, None))
    # 2 x 2 product of (cut, radical end) + 4 far singles
    assert ((HEAD, "primary_radical"), (TAIL, "benzylic_radical")) in alone
    singles = [a for a in alone if sum(e is not None for _k, e in a) == 1]
    assert sorted(e for a in singles for _k, e in a if e) == sorted(
        ["primary_radical", "benzylic_radical", "methyl", "vinylidene", "benzyl", "styryl"])
    # + the 2 x 2 joint combinations of the far non-radical ends
    assert len(alone) == 4 + 4 + 4 and len(set(alone)) == len(alone)
    one = win.context_assignments(w, near_ports=(HEAD,))
    assert len(one) == 4 * 2 + 2
    assert ((HEAD, "methyl"), (TAIL, "benzylic_radical")) in one
    near = win.context_assignments(w, near_ports=(HEAD, TAIL))
    assert len(near) == 16 and len(set(near)) == 16
    assert win.context_assignments(w) == near
    # two participants: the product runs across participants
    pair = ps_pair()
    ports = b.vector_ports(pair)
    assert ports == ((0, "head"), (1, "head"), (1, "tail"))
    every = win.context_assignments(pair)
    assert len(every) == 4 * 4 * 4
    assert (((0, "head"), "methyl"), ((1, "head"), "methyl"), ((1, "tail"), None)) in every
    heads = win.context_assignments(pair, near_ports=((0, "head"), (1, "head")))
    assert len(heads) == 4 * 4 * 2 + 2
    assert (((0, "head"), "methyl"), ((1, "head"), "methyl"), ((1, "tail"), None)) in heads
    # a near end with a far radical end
    assert (((0, "head"), None), ((1, "head"), "methyl"),
            ((1, "tail"), "benzylic_radical")) in heads
    physical_head = ps_witness(radical_end=True)
    only_tail = win.context_assignments(physical_head, near_ports=(TAIL,))
    assert {k for a in only_tail for k, _e in a} == {TAIL}


def test_far_port_reduction_is_checked_and_promotes_on_a_change():
    """A far port's non-radical ends are reduced to singles; a single whose
    window or certificate differs from the all-artificial context promotes
    the port to near (the reduction is a checked claim, never assumed)."""
    w = toy()
    atoms = ball(w, 3, 1)

    def at(assignment, view, status="fixed_point"):
        r = _ctx(w, atoms, view)
        return win.ContextResult(assignment, status, "", r.witness, r.window, r.observation)
    base = at(((HEAD, None), (TAIL, None)), "a")
    same = at(((HEAD, None), (TAIL, "methyl")), "a")
    other = at(((HEAD, "methyl"), (TAIL, None)), "b")
    refused = at(((HEAD, "vinylidene"), (TAIL, None)), None, "refused:domain")
    verdicts = win.far_port_verdicts([base, same, other, refused], far=(HEAD, TAIL))
    assert verdicts == {(TAIL, "methyl"): True, (HEAD, "methyl"): False,
                        (HEAD, "vinylidene"): False}
    assert win.far_port_promotions(verdicts) == (HEAD,)


def test_contexts_with_equal_atoms_but_different_marks_do_not_merge():
    """Equal atom sets and equal certificate views, but the windows differ as
    marked graphs (a reacting mark): no merge."""
    first, second = toy(centre=3), toy(centre=4)
    atoms = ball(first, 3, 1) | ball(first, 4, 1)
    assert win.merge_contexts([_ctx(first, atoms, "a"), _ctx(second, atoms, "a")]) == ()
    assert win.merge_contexts([_ctx(first, atoms, "a"), _ctx(first, atoms, "a")]) == ((0, 1),)


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
    assert win.CLOSURE_EXTENSION == 1
    assert bound.reach == 5   # predicate radius + one closure bond
    # a free-label search adds no extent: its candidates are reacting atoms
    # (otherwise the member is refused)
    assert bound.backbone_atoms == 1 + 2 * 5 + 1 == 12
    assert bound.repeats == math.ceil((bound.backbone_atoms + 1) / 2) == 7
    assert bound.L == (bound.backbone_atoms + 1) // 2 == 6
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


ROSTER = (("forward", ("radical", "radical")), ("reverse", ("alkene", "closed")))


def _species(spec):
    from rmgpy.kmc import channels as chm
    return chm.SpeciesKeys()(b.PROVIDERS["polystyrene"].build("x", spec).molecule)


def test_ladder_covers_cut_phases_and_controls():
    c = win.FamilyConstants("Disproportionation", radius=2, reacting_distance=1,
                            free_label_search=False, bimolecular=True,
                            directions=(("forward", 2), ("reverse", 2)), pair_kinds=ROSTER,
                            classes=())
    bound = win.span_bound(c)
    plan = win.ladder_plan(c, bound)
    items = plan.items
    assert max(i.dp for i in items) == bound.L + 1
    assert all(i.control for i in items if i.dp == bound.L + 1)
    assert not any(i.control for i in items if max(i.dps) <= bound.L)
    assert all(i.control == (max(i.dps) > bound.L) for i in items)
    kinds = {k for i in items for k in i.kinds}
    for kind in ("intact", "scission:phase_A", "scission:phase_B", "pendant_loss",
                 "ch_cleavage:benzylic", "ch_cleavage:secondary", "pair:short_short",
                 "pair:short_long"):
        assert kind in kinds, kind
    specs = [s for i in items for s in i.specs]
    assert any(s.head_shift for s in specs) and any(s.tail_shift for s in specs)
    # every declared end chemistry at every DP 1..L+1 of the intact chains
    intact = {(i.dp, i.specs[0].head_end, i.specs[0].tail_end) for i in items
              if "intact" in i.kinds}
    for dp in range(1, bound.L + 2):
        for head in ("methyl", "primary_radical", "vinylidene"):
            for tail in ("benzyl", "benzylic_radical", "styryl"):
                assert (dp, head, tail) in intact or (dp, head, tail) in plan.excluded_keys
    # every planned species is in the declared domain (the rest is excluded
    # with its domain reason, never silently dropped)
    for item in items:
        for spec in item.specs:
            build = b.PROVIDERS["polystyrene"].build("x", spec)
            assert dom.check(build.molecule).accepted, (item.label, dom.check(build.molecule).reason)
    assert all(reason for _key, reason in plan.excluded)


def test_ladder_pairs_cover_every_dp_pair_kind_and_direction():
    """Short/short pairs: every DP pair (a, b), 1 <= a, b <= L + 1 (any DP
    above L is a control), for every kind combination of the family's
    roster in its own direction; reverse disproportionation pairs a saturated
    chain with an alkene (ethylbenzene + styrene at DP 1)."""
    c = win.FamilyConstants("Disproportionation", radius=2, reacting_distance=1,
                            free_label_search=False, bimolecular=True,
                            directions=(("forward", 2), ("reverse", 2)), pair_kinds=ROSTER,
                            classes=())
    bound = win.span_bound(c)
    L = bound.L
    pairs = [i for i in win.ladder_plan(c, bound).items if "pair:short_short" in i.kinds]
    seen = defaultdict(set)
    for item in pairs:
        assert len(item.specs) == 2 and len(item.dps) == 2 and item.directions
        for direction in item.directions:
            seen[(direction, item.pair_kinds)].add(tuple(sorted(item.dps)))
    grid = {tuple(sorted((a, c_))) for a in range(1, L + 2) for c_ in range(1, L + 2)}
    assert seen[("forward", ("radical", "radical"))] == grid
    assert seen[("reverse", ("alkene", "closed"))] == grid
    assert set(seen) == {("forward", ("radical", "radical")), ("reverse", ("alkene", "closed"))}
    # unequal DPs, e.g. DP2 + DP3 (L = 4), and the controls (L + 1, b)
    assert any(sorted(i.dps) == [2, 3] for i in pairs)
    assert all(i.control for i in pairs if max(i.dps) == L + 1)
    S = b.PsChainSpec
    eb, st = _species(S(1, head_end="methyl", tail_end="benzyl")), \
        _species(S(1, head_end="methyl", tail_end="styryl"))
    rev = [i for i in pairs if "reverse" in i.directions]
    assert any(sorted(_species(x) for x in i.specs) == sorted([eb, st]) for i in rev)
    assert not any("reverse" in i.directions and
                   any(dom.check(b.PROVIDERS["polystyrene"].build("x", x).molecule).details
                       .get("radical_sites") for x in i.specs) for i in pairs)
    fwd = [i for i in pairs if "forward" in i.directions]
    benzylic2 = _species(S(2, head_end="methyl", tail_end="benzylic_radical"))
    primary3 = _species(S(3, head_end="primary_radical", tail_end="benzyl"))
    assert any(sorted(_species(x) for x in i.specs) == sorted([benzylic2, primary3])
               for i in fwd)
    # no pair is planned twice
    keys = [(tuple(sorted(_species(x) for x in i.specs)), i.directions) for i in pairs]
    assert len(keys) == len(set(keys))


def test_ladder_coverage_is_checked_per_direction():
    c = win.FamilyConstants("Disproportionation", radius=2, reacting_distance=1,
                            free_label_search=False, bimolecular=True,
                            directions=(("forward", 2), ("reverse", 2)), pair_kinds=ROSTER,
                            classes=())
    bound = win.span_bound(c)
    plan = win.ladder_plan(c, bound)
    fwd = [{"direction": "forward", "status": "fixed_point", "pair_kinds": ["radical", "radical"],
            "dps": [a, a]} for a in (bound.L, bound.L + 1)]
    gaps = win.ladder_direction_gaps(plan, fwd)
    assert gaps and all("reverse" in g for g in gaps)
    rev = [dict(r, direction="reverse", pair_kinds=["alkene", "closed"]) for r in fwd]
    assert win.ladder_direction_gaps(plan, fwd + rev) == []


@pytest.mark.parametrize("spec,smiles", [
    (dict(units=1, head_shift=True, head_end="benzyl", tail_end="benzyl"), "Cc1ccccc1"),
    (dict(units=1, head_shift=True, head_end="benzyl", tail_end="benzylic_radical"),
     "[CH2]c1ccccc1"),
    (dict(units=0, tail_shift=True, head_end="methyl", tail_end="methyl"), "C"),
    (dict(units=0, tail_shift=True, head_end="primary_radical", tail_end="methyl"), "[CH3]")])
def test_one_carbon_chains_are_built(spec, smiles):
    from rmgpy.molecule import Molecule
    build = b.PROVIDERS["polystyrene"].build("one", b.PsChainSpec(**spec))
    assert build.molecule.is_isomorphic(Molecule(smiles=smiles))
    assert dom.check(build.molecule).accepted
    assert all(not p.artificial for p in build.metadata.ports)
    assert win.chain_backbone_atoms(b.PsChainSpec(**spec)) == 1


def test_one_carbon_species_are_in_the_ladder():
    c = win.FamilyConstants("X", radius=2, reacting_distance=1, free_label_search=False,
                            classes=())
    plan = win.ladder_plan(c, win.span_bound(c))
    keys = {_species(s) for i in plan.items for s in i.specs if i.dp == 1}
    from rmgpy.kmc import channels as chm
    from rmgpy.molecule import Molecule
    for smiles in ("Cc1ccccc1", "[CH2]c1ccccc1", "C", "[CH3]"):
        assert chm.SpeciesKeys()(Molecule(smiles=smiles)) in keys, smiles
    # only the alkene ends of a one-carbon chain are not constructible
    one = [(k, reason) for k, reason in plan.excluded if k[0] in ("C0", "C1")]
    assert one and all("alkene end needs two backbone carbons" in reason or
                       "outside the declared domain" in reason or "duplicate" in reason
                       for _k, reason in one), one
    with pytest.raises(b.MetadataError):
        b.PROVIDERS["polystyrene"].build("x", b.PsChainSpec(1, head_shift=True,
                                                          head_end="styryl",
                                                          tail_end="benzyl"))


def test_span_bound_raises_on_an_oversized_window():
    w = ps_witness(units=8)
    bound = win.span_bound(win.FamilyConstants("X", radius=4, reacting_distance=0,
                                               free_label_search=False))
    with pytest.raises(win.SpanBoundError):
        win.check_span(w, set(w.reactant_rows()), bound)
    small = {a for a, u, k in w.participants[0].metadata.coordinates if u in (3, 4)}
    assert win.check_span(w, small, bound)[0]["units"] == 2


def test_span_bound_is_in_each_providers_own_units():
    """The bound is A backbone atoms; a PS repeat holds two of them, an alkane
    repeat one: seven alkane carbons are seven units, inside A = 11."""
    w = toy(carbons=16, centre=8)
    bound = win.span_bound(win.FamilyConstants("X", radius=4, reacting_distance=0,
                                               free_label_search=False))
    assert bound.backbone_atoms == 11
    seven = set(range(5, 12))
    span = win.check_span(w, seven, bound)
    assert span[0]["backbone_atoms"] == 7 and span[0]["units"] == 7
    assert win.repeats_bound(bound, "linear_alkane") == 11
    assert win.repeats_bound(bound, "polystyrene") == 6
    with pytest.raises(win.SpanBoundError):
        win.check_span(w, set(range(2, 14)), bound)


def test_uncovered_reacting_port_member_fails_the_gate():
    """Every event a refused member stands for needs its own cover: a fixed
    point or a domain refusal.  Certificate failures, ceiling hits, incomplete
    descents and port or end refusals are gaps, never cover, and a covered
    event never stands in for another."""
    fixed = {"kind": "interior", "status": "fixed_point"}
    port = {"kind": "physical", "status": "refused:reacting_port"}
    end = {"kind": "physical", "port": [1, "head"], "end": "vinylidene", "status": "refused:end"}
    assert win.cover_verdict([fixed])
    assert not win.cover_verdict([fixed, end])
    assert not win.cover_verdict([end, port])
    assert not win.cover_verdict([])
    assert not win.cover_verdict([{"kind": "interior", "status": "refused:incomplete"}])
    for status in ("refused:incomplete", "ceiling_exceeded", "certificate_unstable",
                   "no_fixed_point", "refused:recertification: x", "refused:extension",
                   "refused:rematch", "refused:unsupported:spin", "unmatched", "no_witness",
                   "refused:guard: element_count C", None):
        assert win.cover_status_verdict(status) == "gap", status
        assert not win.cover_verdict([{"kind": "interior", "status": status}]), status
    assert win.cover_status_verdict("refused:domain") == "out_of_domain"
    # the family generates no mapping on the event's atoms: no such reaction
    assert win.cover_status_verdict("not_generated") == "no_event"
    assert win.cover_verdict([fixed, {"kind": "physical", "port": [1, "head"],
                                      "end": "vinylidene", "status": "not_generated"}])
    # a physical cover never hides a failed interior event, nor the reverse
    methyl = {"kind": "physical", "port": [1, "head"], "end": "methyl", "status": "fixed_point"}
    assert not win.cover_verdict([methyl, {"kind": "interior", "status": "ceiling_exceeded"}])
    assert not win.cover_verdict([fixed, dict(methyl, status="refused:incomplete")])
    # an event outside the declared domain needs no rate
    outside = dict(end, status="refused:domain")
    assert win.cover_verdict([fixed, methyl, outside])
    # event-keyed rows: every matching member of an event must be covered
    rows = [{"event": "interior@1:head", "status": "fixed_point"},
            {"event": "methyl@1:head", "status": "fixed_point,ceiling_exceeded",
             "members": [{"status": "fixed_point"}, {"status": "ceiling_exceeded"}]}]
    assert not win.cover_verdict(rows)
    assert win.cover_gaps(rows) == ["methyl@1:head: ceiling_exceeded,fixed_point"]
    rows[1]["members"][1]["status"] = "refused:domain"
    assert win.cover_verdict(rows)
    rows = [{"label": "a", "status": "refused:reacting_port", "covers": [fixed, end],
             "covered": win.cover_verdict([fixed, end])},
            {"label": "b", "status": "refused:reacting_port", "covers": [fixed],
             "covered": win.cover_verdict([fixed])},
            {"label": "c", "status": "fixed_point"}]
    gaps = win.reacting_port_gaps(rows)
    assert len(gaps) == 1 and "a" in gaps[0] and "refused:end" in gaps[0]
    # a row without a cover search is a gap too
    assert win.reacting_port_gaps([{"label": "d", "status": "refused:reacting_port"}])


def test_far_port_joint_reduction_is_checked():
    """Single far ends merging with the all-artificial context does not imply
    that their combinations merge: every combination of far non-radical ends
    per participant and across participants is a context, checked, and a
    failed joint check promotes those ports to near."""
    w = ps_witness()
    alone = win.context_assignments(w, near_ports=())
    joint = [a for a in alone if sum(e is not None for _k, e in a) == 2 and
             not any(e in b.RADICAL_ENDS for _k, e in a)]
    assert sorted(joint) == sorted(((HEAD, h), (TAIL, t)) for h in ("methyl", "vinylidene")
                                   for t in ("benzyl", "styryl"))
    pair = ps_pair()
    ports = b.vector_ports(pair)
    assert win.far_port_groups(ports, ()) == (((1, "head"), (1, "tail")), ports)
    across = win.context_assignments(pair, near_ports=())
    assert (((0, "head"), "methyl"), ((1, "head"), "vinylidene"), ((1, "tail"), "styryl")) \
        in across
    toyw = toy()
    atoms = ball(toyw, 3, 1)

    def at(assignment, view):
        r = _ctx(toyw, atoms, view)
        return win.ContextResult(assignment, "fixed_point", "", r.witness, r.window,
                                 r.observation)
    base = at(((HEAD, None), (TAIL, None)), "a")
    first = at(((HEAD, "methyl"), (TAIL, None)), "a")
    second = at(((HEAD, None), (TAIL, "methyl")), "a")
    both = at(((HEAD, "methyl"), (TAIL, "methyl")), "b")
    assert win.contexts_merge(base, first) and win.contexts_merge(base, second)
    assert not win.contexts_merge(base, both)
    verdicts = win.far_port_verdicts([base, first, second, both], far=(HEAD, TAIL))
    assert verdicts[(HEAD, "methyl")] and verdicts[(TAIL, "methyl")]
    assert verdicts[((HEAD, "methyl"), (TAIL, "methyl"))] is False
    assert win.far_port_promotions(verdicts) == (HEAD, TAIL)
    merged = win.far_port_verdicts([base, first, second, at(both.assignment, "a")],
                                   far=(HEAD, TAIL))
    assert all(merged.values()) and win.far_port_promotions(merged) == ()


def test_port_at_reach_plus_one_is_near_and_takes_every_end():
    """A port whose atom is exactly reach + 1 heavy bonds from a reacting atom
    is near (it can lie in a context's window) and takes every declared end;
    one bond farther it is far."""
    w = ps_witness(units=8)
    bound = win.span_bound(win.FamilyConstants("X", radius=2, reacting_distance=0,
                                               free_label_search=False))
    assert bound.reach == 3
    coords = {(u, k): a for a, u, k in w.participants[0].metadata.coordinates}
    edge = coords[(2, 0)]       # backbone index 4 = reach + 1 bonds from the head port
    beyond = coords[(2, 1)]     # index 5
    assert win.near_ports(w, {edge}, bound) == (HEAD,)
    assert win.near_ports(w, {beyond}, bound) == ()
    near = win.near_ports(w, {edge}, bound)
    assignments = win.context_assignments(w, near)
    assert {e for a in assignments for k, e in a if k == HEAD and e} == set(b.HEAD_ENDS)
    # every declared head end combined with the far tail's radical ends
    for end in b.HEAD_ENDS:
        assert ((HEAD, end), (TAIL, "benzylic_radical")) in assignments, end


def test_dp1_styrene_is_built_from_both_alkene_ends():
    """Styrene is vinylidene at one end and styryl at the other: the DP1
    declaration with both alkene ends builds it (and the ladder merges it
    with the other declaration that builds styrene)."""
    from rmgpy.molecule import Molecule
    build = b.PROVIDERS["polystyrene"].build("st", b.PsChainSpec(1, head_end="vinylidene",
                                                               tail_end="styryl"))
    assert build.molecule.is_isomorphic(Molecule(smiles="C=Cc1ccccc1"))
    assert dom.check(build.molecule).details["ends"] == ["styryl", "vinylidene"]
    c = win.FamilyConstants("X", radius=2, reacting_distance=1, free_label_search=False,
                            classes=())
    plan = win.ladder_plan(c, win.span_bound(c))
    assert not any("not constructible" in reason for key, reason in plan.excluded
                   if key[:1] == (1,) and "vinylidene" in key and "styryl" in key)
    with pytest.raises(b.MetadataError):     # a shared carbon of two double bonds
        b.PROVIDERS["polystyrene"].build("x", b.PsChainSpec(2, head_shift=True,
                                                          head_end="styryl",
                                                          tail_end="styryl"))


def test_alkane_proxy_builds_declared_physical_ends():
    """An alkane event at a reacting port needs its physical-end witnesses
    too: the alkane proxy builds a declared end (methyl, primary radical) as
    a physical port, and the event specs carry it."""
    from rmgpy.molecule import Molecule
    build = win.build_spec(("alkane", 5, "primary_radical", "methyl"))
    assert build.molecule.is_isomorphic(Molecule(smiles="[CH2]CCCC"))
    assert all(not p.artificial and not p.cap_atom_ids for p in build.metadata.ports)
    half = win.build_spec(("alkane", 5, "methyl", None))
    assert [p.artificial for p in half.metadata.ports] == [False, True]
    assert win.event_specs([("alkane", 16)], (((0, "head"), "primary_radical"),)) == \
        ([("alkane", 16, "primary_radical", None, None)], {})
    assert win.event_specs([("alkane", 16, None, "methyl")], (((0, "head"), "interior"),)) == \
        ([("alkane", 17, None, "methyl", None)], {0: 1})
    with pytest.raises(b.MetadataError):
        win.build_spec(("alkane", 5, "vinylidene", None))


# --------------------------------------------------------------------------
# the generated roster (reacting-site environment classes)
# --------------------------------------------------------------------------
def P(*args):
    return win.Pattern(*args)


def test_site_classes_derive_from_the_domain_declarations(monkeypatch):
    """The class set is derived from domain.py's declarations: every end group,
    radical kind and small species maps to reacting-site kinds, every site kind
    is carried by some class, the patterns pair every end group with every
    other and add interior radicals up to MAX_RADICAL_SITES; a new declared
    feature without a mapping fails the derivation."""
    kinds = win.declared_site_kinds()
    for end in dom.END_GROUPS:
        assert set(win.END_SITE_KINDS[end]) <= set(kinds), end
    for kind in dom.RADICAL_KINDS:
        assert set(win.RADICAL_SITE_KINDS[kind]) <= set(kinds), kind
    for name in dom.SMALL_SPECIES:
        assert set(win.SMALL_SITE_KINDS[name]) <= set(kinds), name
    for site in ("terminal_primary", "terminal_benzylic", "interior_benzylic",
                 "interior_secondary", "methyl_radical", "benzyl_radical", "h_atom",
                 "backbone_ch2", "backbone_chph", "end_methyl", "end_benzyl",
                 "vinylidene_terminal", "styryl_terminal", "h2"):
        assert site in kinds, site
    patterns = win.participant_patterns()
    chains = [p for p in patterns if p.kind == "chain"]
    assert {e for p in chains for e in p.ends} == set(dom.END_GROUPS)
    assert {tuple(sorted(p.ends)) for p in chains} == {
        tuple(sorted((x, y))) for x in dom.END_GROUPS for y in dom.END_GROUPS}
    assert {k for p in chains for k in p.interior} == set(win.INTERIOR_RADICALS)
    assert all(p.radical_sites <= dom.MAX_RADICAL_SITES for p in patterns)
    assert max(p.radical_sites for p in chains) == dom.MAX_RADICAL_SITES
    assert P("chain", ("benzyl", "methyl"), ("secondary_backbone",)) in patterns
    # mixed species: an unsaturated end and a radical on one chain
    assert P("chain", ("benzylic_radical", "vinylidene"), ()) in patterns
    assert P("chain", ("methyl", "vinylidene"), ("benzylic",)) in patterns
    assert {p.ends[0] for p in patterns if p.kind == "small"} == set(dom.SMALL_SPECIES)
    classes = win.site_classes()
    realized = {k for c in classes for k in c.site}
    assert set(kinds) <= realized, sorted(set(kinds) - realized)
    for c in classes:
        assert c.dp == win.pattern_min_dp(c.pattern)
    secondary = [c for c in classes if c.site == ("interior_secondary",)]
    assert {c.pattern for c in secondary} >= {
        P("chain", ("benzyl", "methyl"), ("secondary_backbone",)),
        P("chain", ("methyl", "vinylidene"), ("secondary_backbone",)),
        P("chain", ("benzylic_radical", "methyl"), ("secondary_backbone",))}
    mixed = [c for c in classes if c.pattern.radical_sites and
             any(e in b.ALKENE_ENDS for e in c.pattern.ends)]
    assert {c.site[0] for c in mixed} >= {"terminal_benzylic", "interior_benzylic",
                                          "interior_secondary", "backbone_chph",
                                          "vinylidene_terminal"}
    monkeypatch.setitem(dom.END_GROUPS, "allyl", "CH2=CH-CH2- (undeclared here)")
    with pytest.raises(win.WindowError):
        win.declared_site_kinds()
    with pytest.raises(win.WindowError):
        win.participant_patterns()
    monkeypatch.delitem(dom.END_GROUPS, "allyl")
    monkeypatch.setattr(dom, "RADICAL_KINDS", dom.RADICAL_KINDS + ("vinylic",))
    with pytest.raises(win.WindowError):
        win.declared_site_kinds()
    monkeypatch.undo()
    monkeypatch.setattr(dom, "SMALL_SPECIES", dom.SMALL_SPECIES + ("CH4",))
    with pytest.raises(win.WindowError):
        win.declared_site_kinds()


def test_pattern_min_dp_is_the_smallest_in_domain_dp():
    S = b.PsChainSpec
    sec = P("chain", ("benzyl", "methyl"), ("secondary_backbone",))
    assert win.pattern_specs(sec, 1) == ()
    assert win.pattern_min_dp(sec) == 2
    for spec in win.pattern_specs(sec, 2):
        assert "interior_secondary" in win.pattern_site_kinds(sec, spec)
    mixed = P("chain", ("benzylic_radical", "vinylidene"), ())
    assert win.pattern_min_dp(mixed) == 2
    assert win.pattern_specs(mixed, 5) == (S(5, head_end="benzylic_radical", head_shift=True,
                                             tail_end="vinylidene", tail_shift=True),)
    two = P("chain", ("methyl", "methyl"), ("benzylic", "benzylic"))
    for spec in win.pattern_specs(two, win.pattern_min_dp(two)):
        verdict = dom.check(b.PROVIDERS["polystyrene"].build("x", spec).molecule)
        assert verdict.accepted and len(verdict.details["radical_sites"]) == 2
    assert win.pattern_specs(P("small", ("H",)), 1) == (("small", "H"),)
    assert dom.check(win.small_build("H").molecule).kind == "small:H"
    assert dom.check(win.small_build("H2").molecule).kind == "small:H2"


def _class(direction, slot, site, pattern, dp, spec):
    return win.SiteClass(direction, slot, site, pattern, dp, spec)


def test_ladder_takes_the_product_of_the_slot_classes():
    """Generated pairs: the product of the slot-0 and slot-1 class
    representatives of a bimolecular direction, each at its smallest DP;
    unimolecular classes alone; every class reaches the plan.  A plan without
    its rosters raises."""
    S = b.PsChainSpec
    sec = _class("forward", 1, (("*3", "interior_secondary"),),
                 P("chain", ("benzyl", "methyl"), ("secondary_backbone",)), 2,
                 S(2, head_end="methyl", tail_end="benzyl", secondary_radical_units=(1,)))
    pe = _class("forward", 1, (("*3", "terminal_benzylic"),),
                P("chain", ("benzylic_radical", "methyl"), ()), 1,
                S(1, head_end="methyl", tail_end="benzylic_radical"))
    eb = _class("forward", 0, (("*1", "end_benzyl"),), P("chain", ("benzyl", "methyl"), ()), 1,
                S(1, head_end="methyl", tail_end="benzyl"))
    mixed = _class("forward", 0, (("*1", "backbone_chph"),),
                   P("chain", ("benzylic_radical", "vinylidene"), ()), 5,
                   S(5, head_end="vinylidene", tail_end="benzylic_radical"))
    c = win.FamilyConstants("H_Abstraction", radius=4, reacting_distance=0,
                            free_label_search=False, bimolecular=True,
                            directions=(("forward", 2),),
                            pair_kinds=(("forward", ("closed", "radical")),),
                            classes=(sec, pe, eb, mixed))
    plan = win.ladder_plan(c)
    generated = [i for i in plan.items if "pair:class" in i.kinds]
    got = {tuple(sorted(_species(x) for x in i.specs)) for i in generated}
    want = {tuple(sorted((_species(x.spec), _species(y.spec)))) for x in (eb, mixed)
            for y in (sec, pe)}
    assert want <= got and len(generated) == len(got)
    item = next(i for i in generated if {repr(x) for x in i.specs} ==
                {repr(mixed.spec), repr(sec.spec)})
    assert sorted(item.dps) == [2, 5] and item.directions == ("forward",)
    assert win.ladder_class_gaps(plan, c.classes) == []
    # the long artificial chains carry every interior site kind with no
    # feature within the span bound: an interior secondary radical too
    longs = [i.specs[1] for i in plan.items if "pair:short_long" in i.kinds]
    assert any(x.secondary_radical_units and not x.head_end and not x.tail_end
               for x in longs)
    assert any(x.radical_units and not x.head_end and not x.tail_end for x in longs)
    assert win.ladder_class_gaps(plan, c.classes + (
        _class("forward", 0, (("*1", "end_methyl"),), P("chain", ("methyl", "methyl"), ()), 1,
               S(1, head_end="methyl", tail_end="methyl", tail_shift=True)),))
    uni = _class("reverse", 0, (("*1", "interior_secondary"),), sec.pattern, 2, sec.spec)
    plan2 = win.ladder_plan(win.FamilyConstants("R", 2, 1, False, False, (("reverse", 1),),
                                                classes=(uni,)))
    assert any(win.class_label(uni) in i.classes and len(i.specs) == 1 for i in plan2.items)
    with pytest.raises(win.WindowError):
        win.ladder_plan(win.FamilyConstants("H_Abstraction", 4, 0, False, True,
                                            (("forward", 2),), classes=(sec,)))
    with pytest.raises(win.WindowError):
        win.ladder_plan(win.FamilyConstants("X", 2, 1, False))


# --------------------------------------------------------------------------
# rework 3: class equivalence, per-reactant reduction, pre-filter, covers
# --------------------------------------------------------------------------
def _cctx(witness, atoms, view, status="fixed_point", centre=3):
    """A closed class member on a toy witness (a stub certificate view that
    ignores the transport)."""
    support = win.Support(atoms=frozenset(atoms), reacting=frozenset({centre}))
    window = win.Window(atoms=frozenset(atoms), support=support)
    obs = win.Observation(support=support, viewer=lambda inside, transport=None, v=view: v)
    return win.ContextResult((), status, "", witness, window, obs)


def _closure(label, results, pathways, *, slot=0, partner="p", direction="forward",
             unpaired=0, elements=(), training=False):
    return win.ClassClosure(label, direction, slot, partner, tuple(results),
                            tuple(sorted(set(pathways))),
                            tuple(sorted({win.own_leaf(p, slot) for p in pathways})),
                            tuple(sorted({r.status for r in results})), unpaired, elements,
                            training)


def test_class_equivalence_merges_beyond_the_window_never_across_leaves():
    """Two classes whose closed windows and certificate views are equal (here
    two chains differing only beyond the window) merge; different leaves
    never merge, even with equal windows and views; unequal views, windows,
    member statuses or partners, a member without a partner in the other
    class, or no fixed point at all: no merge.  A merged-root tree compares
    the windows widened to the span bound (``reach``) and the whole-species
    reads too."""
    near, far = toy(8, 3), toy(12, 3)
    a = _closure("A", [_cctx(near, ball(near, 3, 1), "v")], [("L1", "Y")])
    same = _closure("B", [_cctx(far, ball(far, 3, 1), "v")], [("L1", "Y")])
    assert win.classes_merge(a, same) and win.classes_merge(same, a)
    leaf = _closure("C", [_cctx(far, ball(far, 3, 1), "v")], [("L2", "Y")])
    assert not win.classes_merge(a, leaf)
    assert not win.classes_merge(a, _closure("D", [_cctx(far, ball(far, 3, 1), "w")],
                                             [("L1", "Y")]))
    assert not win.classes_merge(a, _closure("E", [_cctx(far, ball(far, 3, 2), "v")],
                                             [("L1", "Y")]))
    refused = win.ContextResult((), "refused:domain")
    assert not win.classes_merge(a, _closure("F", [_cctx(far, ball(far, 3, 1), "v"), refused],
                                             [("L1", "Y")]))
    assert not win.classes_merge(a, _closure("G", [_cctx(far, ball(far, 3, 1), "v")],
                                             [("L1", "Y")], partner="q"))
    assert not win.classes_merge(a, _closure("H", [_cctx(far, ball(far, 3, 1), "v")],
                                             [("L1", "Y")], slot=1))
    nothing = _closure("R", [refused], [])
    assert not win.classes_merge(nothing, nothing)
    two = _closure("T", [_cctx(near, ball(near, 3, 1), "v"), _cctx(near, ball(near, 3, 2), "v")],
                   [("L1", "Y")])
    assert not win.classes_merge(a, two) and not win.classes_merge(two, a)
    # a merged root: the windows widened to the span bound decide (``reach``)
    assert win.classes_merge(a, same, reach=2)
    assert not win.classes_merge(a, same, reach=4)       # toy(8)'s end lies within 4 bonds
    assert not win.classes_merge(a, _closure("U", [_cctx(far, ball(far, 3, 1), "v")],
                                             [("L1", "Y")], unpaired=1), reach=2)


def test_merge_classes_maps_each_class_to_its_cheapest_equivalent():
    """Classes merge per direction and slot into the cheapest equivalent class
    (smallest DP first); a class with another leaf, in another slot or
    without a closing stays its own representative; the summary names the
    map."""
    S = b.PsChainSpec
    near, far = toy(8, 3), toy(12, 3)
    eb3 = _class("forward", 0, (("*1", "end_methyl"),), P("chain", ("benzyl", "methyl"), ()), 3,
                 S(3, head_end="methyl", tail_end="benzyl"))
    st2 = _class("forward", 0, (("*1", "end_methyl"),), P("chain", ("methyl", "styryl"), ()), 2,
                 S(2, head_end="methyl", tail_end="styryl"))
    ch2 = _class("forward", 0, (("*1", "backbone_ch2"),), P("chain", ("benzyl", "methyl"), ()),
                 3, S(3, head_end="methyl", tail_end="benzyl"))
    other = _class("forward", 1, (("*3", "terminal_benzylic"),),
                   P("chain", ("benzylic_radical", "methyl"), ()), 1,
                   S(1, head_end="methyl", tail_end="benzylic_radical"))
    lone = _class("forward", 0, (("*1", "end_benzyl"),), P("chain", ("benzyl", "methyl"), ()), 3,
                  S(3, head_end="methyl", tail_end="benzyl"))
    L = win.class_label
    closures = [_closure(L(eb3), [_cctx(far, ball(far, 3, 1), "v")], [("C/H3/Cs", "Y")]),
                _closure(L(st2), [_cctx(near, ball(near, 3, 1), "v")], [("C/H3/Cs", "Y")]),
                _closure(L(ch2), [_cctx(far, ball(far, 3, 1), "v")], [("C/H2/NonDeC", "Y")]),
                _closure(L(other), [_cctx(near, ball(near, 3, 1), "v")], [("X", "C_rad")],
                         slot=1)]
    merge = win.merge_classes([eb3, st2, ch2, other, lone], closures)
    assert merge.rep(L(eb3)) == L(st2) and merge.rep(L(st2)) == L(st2)
    assert merge.rep(L(ch2)) == L(ch2) and merge.rep(L(other)) == L(other)
    assert merge.rep(L(lone)) == L(lone) and L(lone) not in merge.closed
    assert win.closure_gaps([eb3, st2, ch2, other, lone], merge) == \
        [f"class {L(lone)}: no closing"]
    table = win.merge_summary([eb3, st2, ch2, other, lone], merge)
    assert table["forward|slot0"]["classes"] == 4
    assert table["forward|slot0"]["representatives"] == 3
    assert sorted(table["forward|slot0"]["map"][L(st2)]) == sorted([L(eb3), L(st2)])


def _reduced_constants(per_reactant):
    """Toy H_Abstraction roster: slot 0 spec s1 carries two representative
    classes (leaves A, B) and one merged class; spec s2 leaf B; slot 1 spec
    t1 leaves P, Q; spec t2 leaf Q."""
    S = b.PsChainSpec
    s1, s2 = S(2, head_end="methyl", tail_end="benzyl"), S(1, head_end="methyl", tail_end="benzyl")
    s3 = S(3, head_end="methyl", tail_end="benzyl")
    t1 = S(1, head_end="methyl", tail_end="benzylic_radical")
    t2 = S(1, head_end="primary_radical", tail_end="benzyl")
    eb = P("chain", ("benzyl", "methyl"), ())
    x1 = _class("forward", 0, (("*1", "end_methyl"),), eb, 2, s1)
    x2 = _class("forward", 0, (("*1", "backbone_ch2"),), eb, 2, s1)
    x3 = _class("forward", 0, (("*1", "end_benzyl"),), eb, 1, s2)
    xm = _class("forward", 0, (("*1", "backbone_chph"),), eb, 3, s3)
    y1 = _class("forward", 1, (("*3", "terminal_benzylic"),),
                P("chain", ("benzylic_radical", "methyl"), ()), 1, t1)
    y2 = _class("forward", 1, (("*3", "interior_benzylic"),),
                P("chain", ("benzylic_radical", "methyl"), ()), 1, t1)
    y3 = _class("forward", 1, (("*3", "terminal_primary"),),
                P("chain", ("benzyl", "primary_radical"), ()), 1, t2)
    L = win.class_label
    classes = (x1, x2, x3, xm, y1, y2, y3)
    leaves = {L(x1): ("A",), L(x2): ("B",), L(x3): ("B",), L(xm): ("A",), L(y1): ("P",),
              L(y2): ("Q",), L(y3): ("Q",)}
    merge = win.ClassMerge(tuple(sorted((L(c), L(x1) if c is xm else L(c)) for c in classes)),
                           tuple(sorted(leaves.items())), frozenset(leaves), frozenset(),
                           frozenset(leaves))
    c = win.FamilyConstants("H_Abstraction", radius=4, reacting_distance=0,
                            free_label_search=False, bimolecular=True,
                            directions=(("forward", 2),),
                            pair_kinds=(("forward", ("closed", "radical")),), classes=classes,
                            per_reactant=per_reactant, merge=merge)
    return c, classes, (s1, s2, s3, t1, t2)


def test_per_reactant_plan_covers_leaf_tuples_merged_root_keeps_the_product():
    """Per-reactant tree: no class product; the classes are covered by their
    closings, and pairs of representative specs cover every leaf tuple that
    occurs (A, B x P, Q: one pair here, s1 + t1).  The same roster on a
    merged-root direction (not in ``per_reactant``) keeps the product of the
    representative specs (s1, s2 x t1, t2) and never expands a merged-away
    class's spec (s3)."""
    from dataclasses import replace
    c, classes, (s1, s2, s3, t1, t2) = _reduced_constants(("forward",))
    plan = win.ladder_plan(c)
    leaf = [i for i in plan.items if "pair:leaf" in i.kinds]
    assert {t for i in leaf for t in i.leaf_tuples} == {("A", "P"), ("A", "Q"), ("B", "P"),
                                                        ("B", "Q")}
    assert len(leaf) == 1 and {_species(x) for x in leaf[0].specs} == {_species(s1),
                                                                       _species(t1)}
    assert not [i for i in plan.items if "pair:class" in i.kinds]
    assert win.ladder_class_gaps(plan, classes) == []
    assert {label for label, _how in plan.class_cover} == {win.class_label(x) for x in classes}
    merged_root = win.ladder_plan(replace(c, per_reactant=()))
    product = [i for i in merged_root.items if "pair:class" in i.kinds]
    got = {tuple(sorted(_species(x) for x in i.specs)) for i in product}
    want = {tuple(sorted((_species(x), _species(y)))) for x in (s1, s2) for y in (t1, t2)}
    assert got == want
    assert not [i for i in merged_root.items if "pair:leaf" in i.kinds]
    assert not any(_species(s3) in {_species(x) for x in i.specs} for i in product)
    assert win.ladder_class_gaps(merged_root, classes) == []
    # a class pair expands the members carrying a representative site in each
    # slot (the merged-away class's site, backbone_chph, is not one of them);
    # merged with a grid pair of the same species it runs every member
    entries, _excluded = win._class_plan(replace(c, per_reactant=()))
    s1_t1 = next(e for e in entries if e[3] == "pair:class" and
                 [repr(x) for x, _dp in e[1]] == [repr(s1), repr(t1)])
    assert all(not i.sites for i in product if "pair:short_short" in i.kinds)
    assert {(slot, site) for _k, slot, site in s1_t1[5]} == {
        (0, (("*1", "end_methyl"),)), (0, (("*1", "backbone_ch2"),)),
        (1, (("*3", "terminal_benzylic"),)), (1, (("*3", "interior_benzylic"),))}
    # without closings the per-reactant direction keeps the unreduced product
    unreduced = win.ladder_plan(replace(c, merge=None))
    assert len([i for i in unreduced.items if "pair:class" in i.kinds]) == 3 * 2


def test_leaf_tuple_coverage_check_needs_a_fixed_point_per_tuple():
    """Every planned leaf tuple needs a fixed-point member with that pathway
    in the item's direction; a refused member, another direction or a
    missing tuple is a gap."""
    item = win.LadderItem("p", ("pair:leaf",), 1, (), directions=("forward",),
                          leaf_tuples=(("A", "P"), ("B", "P")))
    plan = win.LadderPlan("F", 4, (item,), ())
    ok = {"status": "fixed_point", "direction": "forward", "pathway": ["A", "P"]}
    refused = {"status": "refused:domain", "direction": "forward", "pathway": ["B", "P"]}
    gaps = win.leaf_tuple_gaps(plan, [ok, refused])
    assert len(gaps) == 1 and "B;P" in gaps[0]
    assert win.leaf_tuple_gaps(plan, [ok, dict(refused, status="fixed_point")]) == []
    assert len(win.leaf_tuple_gaps(plan, [ok, dict(refused, status="fixed_point",
                                                   direction="reverse")])) == 1
    assert len(win.leaf_tuple_gaps(plan, [])) == 2
    chosen, left = win.leaf_cover(({"s1": {"A", "B"}, "s2": {"B"}}, {"t1": {"P"}, "t2": {"Q"}}))
    assert {t for _x, _y, ts in chosen for t in ts} == {("A", "P"), ("A", "Q"), ("B", "P"),
                                                        ("B", "Q")} and left == ()
    assert len(chosen) == 2
    _chosen, left = win.leaf_cover(({"s1": {"A"}}, {"t1": {"P"}}), lambda k0, k1: "dropped")
    assert left == (("A", "P"),)


def test_pair_prefilter_drops_only_proved_refusals():
    """One-product templates: the product carries both reactants' unpaired
    electrons plus the recipe's delta; more than MAX_RADICAL_SITES drops the
    pair, exactly MAX keeps it; two-product templates never drop."""
    S = b.PsChainSpec
    c = win.FamilyConstants("R_Addition_MultipleBond", 5, 2, False, True,
                            (("forward", 2), ("reverse", 1)), joined=(("forward", 0),))
    alkene_radical, = win.pattern_specs(P("chain", ("benzylic_radical", "vinylidene"), ()), 2)
    styrene = S(1, head_end="vinylidene", tail_end="styryl")
    biradical = S(3, head_end="primary_radical", tail_end="benzylic_radical")
    methyl = win.ONE_CARBON["methyl_radical"]()
    assert win.unpaired_electrons(alkene_radical) == 1 and win.unpaired_electrons(biradical) == 2
    assert win.pair_prefilter(c, "forward", alkene_radical, biradical)
    assert win.pair_prefilter(c, "forward", styrene, biradical) is None
    assert win.pair_prefilter(c, "forward", alkene_radical, methyl) is None
    assert win.pair_prefilter(c, "reverse", alkene_radical, biradical) is None
    two = win.FamilyConstants("H_Abstraction", 4, 0, False, True, (("forward", 2),))
    assert win.pair_prefilter(two, "forward", biradical, biradical) is None


def test_no_event_only_when_no_mapping_touches_the_event_atoms():
    """Manager ruling (round 2): ``no_event`` is accepted when RMG generates no
    mapping on the event's atoms under any labels; the atoms mapped under
    other labels, or a heavy atom missing from the witness, stay a gap.  A
    product-side domain refusal is out of domain (covered)."""
    event = {frozenset({("*1", (1, 2, 0)), ("*2", ("H", (1, 2, 0))), ("*3", (0, 0, 1))})}
    present = {(1, 2, 0), (0, 0, 1), (1, 1, 1)}
    elsewhere = [("m0", {frozenset({("*1", (1, 1, 1)), ("*2", ("H", (1, 1, 1))),
                                    ("*3", (0, 0, 1))})})]
    status, _why = win.unmatched_event(event, elsewhere, present)
    assert status == "not_generated" and win.cover_status_verdict(status) == "no_event"
    relabelled = [("m1", {frozenset({("*3", (1, 2, 0)), ("*2", ("H", (1, 2, 0))),
                                     ("*1", (0, 0, 1))})})]
    status, why = win.unmatched_event(event, elsewhere + relabelled, present)
    assert status == "unmatched" and "1 map its atoms under other labels" in why
    assert win.cover_status_verdict(status) == "gap"
    status, _why = win.unmatched_event(event, elsewhere, {(0, 0, 1)})
    assert status == "unmatched"
    assert win.cover_verdict([{"event": "e", "status": "not_generated"}])
    assert not win.cover_verdict([{"event": "e", "status": "unmatched"}])
    products = {"status": "refused:domain",
                "reason": "refused:domain: products species at atom 3: a u2 site"}
    assert win.cover_verdict([{"event": "e", "members": [products]}])
    assert win._event_verdict([products["status"]]) == "out_of_domain"


# --------------------------------------------------------------------------
# boundaries: declared sites, end verification
# --------------------------------------------------------------------------
@pytest.mark.parametrize("flags", [dict(secondary_radical_units=(3,)),
                                   dict(pendant_loss_units=(3,)), dict(radical_units=(3,)),
                                   dict(secondary_radical_units=(0,))])
def test_requested_sites_that_do_not_exist_are_refused(flags):
    with pytest.raises(b.MetadataError):
        b.PROVIDERS["polystyrene"].build("x", b.PsChainSpec(3, head_end="methyl",
                                                          tail_end="benzyl", **flags))


def _tampered(built, atom, *, radical=None, bond=None):
    rows = dict(built.reactant_rows())
    prows = dict(built.product_rows())
    for side in (rows, prows):
        if radical is not None:
            row = side[atom]
            side[atom] = row[:3] + (radical,) + row[4:]
        if bond is not None:
            x, y, order = bond
            for a, n in ((x, y), (y, x)):
                row = side[a]
                side[a] = row[:5] + (tuple(sorted((m, order if m == n else o)
                                                  for m, o in row[5])),)
    return b.Witness(tuple(sorted(rows.items())), tuple(sorted(prows.items())),
                     built.participants, built.id_base, built.extension_counts, built.label)


def test_verify_end_checks_survivor_electrons_and_bonds():
    w = ps_witness()
    built, record = b.reconstruct_end(w, HEAD, "methyl")
    b.verify_end(w, record, built)
    partner = _partner(w, HEAD)
    with pytest.raises(b.MapCheckError):
        b.verify_end(w, record, _tampered(built, partner, radical=1))
    with pytest.raises(b.MapCheckError):
        b.verify_end(w, record, _tampered(built, None, bond=(w.port(HEAD).atom_id, partner, 2.0)))


def test_end_never_sets_a_bond_the_reaction_edits():
    """A provider whose end edit sets a bond the reaction changes (without the
    partner guard of the standard edit): reconstruct_end refuses it."""
    w = toy(centre=2)
    product = dict(w.product_rows())
    for a, n in ((0, 1), (1, 0)):
        row = product[a]
        product[a] = row[:5] + (tuple(sorted((m, 2.0 if m == n else o) for m, o in row[5])),)
    edited = b.make_witness(w.reactants, tuple(sorted(product.items())), w.participants)

    class Loose(b.LinearAlkaneProvider):
        def end_edit(self, metadata, port, end, rows, protected):
            return b.EndEdit(end, tuple(port.cap_atom_ids), bonds=((0, 1, 2.0),))
    providers = dict(b.PROVIDERS, linear_alkane=Loose())
    with pytest.raises(b.EndConflict):
        b.reconstruct_end(edited, HEAD, "vinylidene", providers)


# --------------------------------------------------------------------------
# item 2: closure (no database)
# --------------------------------------------------------------------------
def test_closure_runs_on_both_reaction_sides():
    """The product side makes a 3=4 double bond: the pi system {3, 4} exists
    on the product side only, so atom 4 (and its H) is in the window only if
    the products are closed too."""
    w = toy()
    product = dict(w.product_rows())
    for a, n in ((3, 4), (4, 3)):
        row = product[a]
        edges = tuple(sorted(((x, 2.0 if x == n else o) for x, o in row[5])))
        product[a] = row[:3] + (0,) + row[4:5] + (edges,)
    w2 = b.make_witness(w.reactants, tuple(sorted(product.items())), w.participants)
    window = win.close(w2, win.Support(atoms=frozenset({3}), reacting=frozenset({3, 4})))
    assert window.refusal is None
    assert 4 in window.atoms
    h4 = {n for n, _o in w2.reactant_rows()[4][5] if w2.reactant_rows()[n][0] == "H"}
    assert h4 <= window.atoms
    reactant_only = win.close(w2, win.Support(atoms=frozenset({3}), reacting=frozenset({3})),
                              sides=("reactants",))
    assert 4 not in reactant_only.atoms


def ring_witness(size, centre=0):
    from rmgpy.molecule import Molecule
    m = Molecule(smiles="C1" + "C" * (size - 1) + "1")
    for k, atom in enumerate(m.atoms):
        atom.id = k
    rows = rm.molecule_rows(m)
    product = dict(rows)
    product[centre] = rows[centre][:3] + (1,) + rows[centre][4:]
    return b.make_witness(tuple(sorted(rows.items())), tuple(sorted(product.items())), ())


def test_ring_closure_and_refusals():
    eight = win.close(ring_witness(8), win.Support(atoms=frozenset({0}),
                                                   reacting=frozenset({0})))
    assert eight.refusal is None and set(range(8)) <= eight.atoms
    nine = win.close(ring_witness(9), win.Support(atoms=frozenset({0}),
                                                  reacting=frozenset({0})))
    assert nine.refusal and nine.refusal.startswith("refused:ring")
    assert nine.atoms == frozenset({0})       # refused, never truncated to a partial ring
    touching = win.close(ring_witness(9), win.Support(atoms=frozenset({0, 1}),
                                                      reacting=frozenset({20})))
    assert touching.refusal and touching.refusal.startswith("refused:large_ring")


def test_fixed_point_reports_observation_refusals():
    """A domain (or any certificate) refusal stops the loop by name, without
    a window: never a truncated one."""
    def observe(witness):
        return win.Observation(win.Support(frozenset()), lambda inside: None, None,
                               "refused:domain: an aryl radical product")
    res = win.fixed_point(toy(), observe)
    assert res.status == "refused:domain: an aryl radical product"
    assert res.window is None


# --------------------------------------------------------------------------
# items 1, 3, 5 on real certificates (RMG database)
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def rmg():
    import rmg_moves_harness as harness
    if not harness.DATABASE.is_dir():
        pytest.skip(f"needs the RMG database at {harness.DATABASE}")
    import channelsTest as ct
    from rmgpy.kmc import certificates as cm
    fam = ct.database().kinetics.families
    return fam, cm.Certifier(fam, max_channel_sets=16)


@pytest.fixture(scope="module")
def habs_cases(rmg):
    fam, certifier = rmg
    builds = [b.PROVIDERS["polystyrene"].build("pe", b.PsChainSpec(1, tail_end="benzylic_radical")),
              b.PROVIDERS["polystyrene"].build("dp3", b.PsChainSpec(3))]
    return win.member_cases(certifier, fam["H_Abstraction"], builds, "forward", label="PE+DP3")


def _interior(case, unit=1, position=0):
    coords = {a: (p.index, u, k) for p in case.witness.participants
              for a, u, k in p.metadata.coordinates}
    return coords.get(case.moves[0].atom("*1")) == (1, unit, position) and \
        coords.get(case.moves[0].atom("*3")) == (0, 0, 1)


def test_failed_predicate_outside_the_ball_pulls_its_atoms_into_the_window(rmg, habs_cases):
    fam, certifier = rmg
    case = next(c for c in habs_cases if _interior(c))
    cert = certifier.certify(fam["H_Abstraction"], win.template_molecules(case.witness,
                                                                         "forward"),
                             member_labels=case.wanted)
    full = win.certificate_support(cert, case.wanted, classifier=win.SearchClassifier(
        fam["H_Abstraction"]))
    matched = win.certificate_support(cert, case.wanted, include_failed=False,
                                      classifier=win.SearchClassifier(fam["H_Abstraction"]))
    assert full.refusal is None
    extra = set(full.atoms) - set(matched.atoms)
    assert extra, "no failed predicate reaches beyond the matched balls"
    rows = case.witness.reactant_rows()
    beyond = {a for a in extra if any(n in matched.atoms for n, _o in rows[a][5])}
    assert beyond      # one bond outside the matched predicates' support
    failed_only = {a for node, result, atoms in full.sources if result == "failed"
                   for a in atoms} - set(matched.atoms)
    assert beyond <= failed_only
    window = win.close(case.witness, full)
    assert beyond <= window.atoms
    observed = win.MemberObserver(certifier, fam["H_Abstraction"], "forward", case.wanted)
    assert beyond <= set(observed(case.witness).support.atoms)


def test_real_member_reaches_a_recertified_fixed_point(rmg, habs_cases):
    fam, certifier = rmg
    case = next(c for c in habs_cases if _interior(c))
    bound = win.span_bound(win.family_constants(fam["H_Abstraction"]))
    res = win.run_member(certifier, fam["H_Abstraction"], case, bound=bound, contexts=False)
    assert res["status"] == "fixed_point", res
    assert res["span_bound_held"] and not res["mapping_problems"]
    assert not res["boundary_bond_problems"]
    assert res["observations"] >= 1 + len(b.vector_ports(case.witness))
    for index, span in res["span"].items():
        assert span["backbone_atoms"] <= bound.backbone_atoms


def test_domain_refused_member_is_refused_not_truncated(rmg, habs_cases):
    fam, certifier = rmg
    builds = [b.PROVIDERS["polystyrene"].build("pe", b.PsChainSpec(1, tail_end="benzylic_radical")),
              b.PROVIDERS["polystyrene"].build("dp3", b.PsChainSpec(3))]
    stats = {}
    every = win.member_cases(certifier, fam["H_Abstraction"], builds, "forward",
                             skip_ring_changes=False)
    win.member_cases(certifier, fam["H_Abstraction"], builds, "forward", stats=stats)
    ring = win.ring_atoms(every[0].reactants)
    family = fam["H_Abstraction"]

    def touches_ring(case):
        change = win.connectivity_change(family, next(iter(case.wanted)))
        return any(a in ring for _k, pair in change for a in pair)
    aryls = [c for c in every if touches_ring(c)]
    assert all(win.species_domain_refusal(c.witness) for c in aryls)
    # the ring-change mappings the enumeration skips are exactly these members'
    assert stats["ring_change_mappings"] == sum(len(c.moves) for c in aryls) > 0
    assert len(every) - len(aryls) == len(habs_cases)
    aryl = next(c for c in aryls if c.moves[0].atom("*1") in ring)
    res = win.run_member(certifier, fam["H_Abstraction"], aryl, contexts=False)
    assert res["status"] == "refused:domain", res
    assert "window_atoms" not in res
    # the certificate path refuses it too (eligibility)
    observer = win.MemberObserver(certifier, fam["H_Abstraction"], "forward", aryl.wanted)
    obs = observer(aryl.witness)
    assert obs.refusal and obs.refusal.startswith("refused:domain")
    assert win.fixed_point(aryl.witness, observer).window is None


def test_scoped_certificate_equals_the_full_one_on_the_member(rmg, habs_cases):
    """The scoped certifier keeps every mapping of the member and of any other
    pathway with its outcome; its view of the member equals the full
    certificate's."""
    fam, certifier = rmg
    family = fam["H_Abstraction"]
    for case in [c for c in habs_cases if _interior(c)][:2]:
        mols = win.template_molecules(case.witness, "forward")
        full = certifier.certify(family, mols, member_labels=case.wanted)
        scoped = win.scoped_certifier(certifier, family, case.wanted).certify(
            family, mols, member_labels=case.wanted)
        inside = frozenset(case.witness.reactant_rows())
        assert scoped["complete"] and scoped["eligible"]
        assert win.certificate_view(full, case.wanted, inside) == \
            win.certificate_view(scoped, case.wanted, inside)
        assert win.member_entries(full, case.wanted)[0][1]["mappings"] == \
            win.member_entries(scoped, case.wanted)[0][1]["mappings"]
        def ran(cs):
            return sum(1 for m in cs.moves if m.status != "excluded")
        assert ran(certifier.channel_set(family, mols)) > ran(
            win.scoped_certifier(certifier, family, case.wanted).channel_set(family, mols))


def _coords(case):
    return {a: (p.index, u, k) for p in case.witness.participants
            for a, u, k in p.metadata.coordinates}


def _habs(rmg, donor, at):
    fam, certifier = rmg
    S = b.PsChainSpec
    specs = [S(1, tail_end="benzylic_radical"), donor]
    builds = [b.PROVIDERS["polystyrene"].build(f"p{k}", x) for k, x in enumerate(specs)]
    cases = win.member_cases(certifier, fam["H_Abstraction"], builds, "forward")
    picked = [c for c in cases if _coords(c).get(c.moves[0].atom("*1")) == at and
              _coords(c).get(c.moves[0].atom("*3")) == (0, 0, 1)]
    assert picked, at
    return specs, picked


def test_contexts_run_across_participants_and_with_far_radical_ends(rmg):
    """H abstraction by 1-phenylethyl from DP2 (donor unit 0 position 1): both
    participants' heads physical methyl ends is a context (whatever the
    descent reads), it reaches a fixed point and merges with neither
    single-head context.  With a DP8 donor, the near head methyl end with the
    far tail benzylic radical is a context too."""
    fam, certifier = rmg
    family = fam["H_Abstraction"]
    bound = win.span_bound(win.family_constants(family))
    _specs, (case, *_rest) = _habs(rmg, b.PsChainSpec(2), (1, 0, 1))
    near = win.member_near_ports(case, bound)
    assert set(near) == set(b.vector_ports(case.witness))
    both = (((0, "head"), "methyl"), ((1, "head"), "methyl"), ((1, "tail"), None))
    assert both in win.context_assignments(case.witness, near)
    base = win.run_member(certifier, family, case, bound=bound, contexts=False)
    assert base["status"] == "fixed_point", base
    depths = tuple((tuple(p), d) for p, d in base["depths"])
    runs = [win.run_context(certifier, family, case, a, bound=bound, start_depths=depths)
            for a in (both, (((0, "head"), "methyl"), ((1, "head"), None), ((1, "tail"), None)),
                      (((0, "head"), None), ((1, "head"), "methyl"), ((1, "tail"), None)))]
    assert [r.status for r in runs] == ["fixed_point"] * 3, [(r.status, r.reason) for r in runs]
    assert not win.contexts_merge(runs[1], runs[0]) and not win.contexts_merge(runs[2], runs[0])
    _specs, (long_case, *_rest) = _habs(rmg, b.PsChainSpec(8), (1, 0, 1))
    near = win.member_near_ports(long_case, bound)
    assert (1, "head") in near and (1, "tail") not in near
    mixed = (((0, "head"), None), ((1, "head"), "methyl"), ((1, "tail"), "benzylic_radical"))
    assert mixed in win.context_assignments(long_case.witness, near)
    run = win.run_context(certifier, family, long_case, mixed, bound=bound)
    assert run.status == "fixed_point", (run.status, run.reason)


def test_alkane_window_is_bounded_in_alkane_units(rmg):
    """1-phenylethyl abstracting H from C8 of C16: the window spans seven
    alkane carbons (units), inside the family's bound in alkane units."""
    fam, certifier = rmg
    family = fam["H_Abstraction"]
    bound = win.span_bound(win.family_constants(family))
    builds = [b.PROVIDERS["polystyrene"].build("pe", b.PsChainSpec(1, tail_end="benzylic_radical")),
              b.PROVIDERS["linear_alkane"].build(16)]
    cases = win.member_cases(certifier, family, builds, "forward")
    case = next(c for c in cases if _coords(c).get(c.moves[0].atom("*1")) == (1, 8, 0))
    res = win.run_member(certifier, family, case, bound=bound, contexts=False)
    assert res["status"] == "fixed_point", res
    assert res["span_bound_held"]
    alkane = res["span"]["1"]
    assert alkane["units"] == alkane["backbone_atoms"] <= bound.backbone_atoms


def test_pair_roster_is_derived_from_the_templates_per_direction(rmg):
    fam, certifier = rmg
    roster = set(win.pair_roster(certifier, fam["Disproportionation"]))
    assert ("reverse", ("alkene", "closed")) in roster
    assert ("forward", ("radical", "radical")) in roster
    habs = set(win.pair_roster(certifier, fam["H_Abstraction"]))
    assert {("forward", ("closed", "radical")), ("forward", ("alkene", "radical")),
            ("forward", ("radical", "radical"))} <= habs
    assert all(d == "forward" for d, _k in habs)
    add = set(win.pair_roster(certifier, fam["R_Addition_MultipleBond"]))
    assert ("forward", ("alkene", "radical")) in add and ("forward", ("closed", "radical")) \
        not in add
    # unimolecular directions carry no pair kinds
    assert not any(d == "reverse" for d, _k in add)


def test_reverse_disproportionation_of_ethylbenzene_and_styrene(rmg):
    """Saturated + alkene (reverse disproportionation): RMG generates in-domain
    members; the first reaches a fixed point holding both ends of both
    participants."""
    fam, certifier = rmg
    family = fam["Disproportionation"]
    bound = win.span_bound(win.family_constants(family))
    S = b.PsChainSpec
    builds = [b.PROVIDERS["polystyrene"].build("eb", S(1, head_end="methyl", tail_end="benzyl")),
              b.PROVIDERS["polystyrene"].build("st", S(1, head_end="methyl", tail_end="styryl"))]
    cases = [c for c in win.member_cases(certifier, family, builds, "reverse")
             if not win.species_domain_refusal(c.witness)]
    assert cases
    res = win.run_member(certifier, family, cases[0], bound=bound, contexts=False)
    assert res["status"] == "fixed_point", res
    assert all(s["both_ends"] for s in res["span"].values())


def test_reacting_port_member_has_a_covering_witness(rmg):
    """Abstraction from the DP2 head carbon (an artificial port atom): refused
    as a reacting port; it stands for the interior event (the chain goes on:
    a methylene H, ``C/H2/NonDeC``) and one event per declared head end, each
    covered by its own witness: interior and methyl reach fixed points, the
    primary-radical and vinylidene ends make species the domain refuses."""
    fam, certifier = rmg
    family = fam["H_Abstraction"]
    bound = win.span_bound(win.family_constants(family))
    specs, cases = _habs(rmg, b.PsChainSpec(2), (1, 0, 0))
    for case in cases:
        res = win.run_member(certifier, family, case, bound=bound, contexts=False, specs=specs)
        assert res["status"] == "refused:reacting_port", res
        assert res["covered"], res["covers"]
        events = {c["event"]: c for c in res["covers"]}
        assert set(events) == {f"{e}@1:head" for e in ("interior",) + b.HEAD_ENDS}
        assert events["interior@1:head"]["verdict"] == "covered"
        assert {tuple(m["pathway"]) for m in events["interior@1:head"]["members"]} == {
            ("C/H2/NonDeC", "C_rad/H/CbCs")}
        assert events["methyl@1:head"]["verdict"] == "covered"
        assert events["primary_radical@1:head"]["verdict"] == "out_of_domain"
        assert "u2" in events["primary_radical@1:head"]["reason"]
        # a product-side domain refusal is out of domain (manager ruling, round 2)
        assert "products species" in events["primary_radical@1:head"]["reason"]
        assert events["vinylidene@1:head"]["verdict"] == "out_of_domain"
    # a ceiling hit on the interior event is a gap the physical covers never hide
    capped = win.run_member(certifier, family, cases[0], bound=bound, contexts=False,
                            specs=specs, ceiling=80)
    events = {c["event"]: c for c in capped["covers"]}
    assert events["interior@1:head"]["status"] == "ceiling_exceeded"
    assert events["methyl@1:head"]["verdict"] == "covered"
    assert capped["covered"] is False
    assert win.reacting_port_gaps([capped])


def test_generated_roster_carries_the_reproduced_classes(rmg):
    """The two environments no round-1 ladder pair produced: an interior
    secondary radical abstracting from ethylbenzene (RMG leaf
    ``C/H2/CbCs;C_rad/H/NonDeC``), and 1-phenylethyl abstracting a benzylic H
    from a chain with a vinylidene end and a benzylic end radical (in domain
    from DP5: at DP4 the product radicals sit four bonds apart).  Both classes
    come out of the family's probe with those smallest DPs, and the ladder's
    product pairs carry them to fixed points."""
    from dataclasses import replace
    fam, certifier = rmg
    family = fam["H_Abstraction"]
    constants = win.family_constants(family)
    S = b.PsChainSpec
    sec = P("chain", ("benzyl", "methyl"), ("secondary_backbone",))
    mixed = P("chain", ("benzylic_radical", "vinylidene"), ())
    eb, pe = P("chain", ("benzyl", "methyl"), ()), P("chain", ("benzylic_radical", "methyl"), ())
    classes = win.class_roster(certifier, family, patterns=[sec, mixed, eb, pe],
                               constants=constants)
    radical = [c for c in classes if c.pattern == sec and
               ("*3", "interior_secondary") in c.site]
    assert radical and min(c.dp for c in radical) == 2, radical
    donor = [c for c in classes if c.pattern == mixed and ("*1", "backbone_chph") in c.site]
    assert donor and min(c.dp for c in donor) == 5, donor
    plan = win.ladder_plan(replace(constants, pair_kinds=(("forward", ("closed", "radical")),),
                                   classes=classes))
    assert win.ladder_class_gaps(plan, classes) == []
    generated = [i for i in plan.items if "pair:class" in i.kinds]

    def item_of(*specs):
        want = sorted(_species(x) for x in specs)
        return next(i for i in generated if sorted(_species(x) for x in i.specs) == want)
    bound = win.span_bound(constants)
    first = item_of(min(radical, key=lambda c: c.dp).spec, S(1, head_end="methyl",
                                                           tail_end="benzyl"))
    rows = [win.run_member(certifier, family, c, bound=bound, contexts=False)
            for c in win.member_cases(certifier, family, win.builds_of(first.specs), "forward")
            if c.pathway == ("C/H2/CbCs", "C_rad/H/NonDeC")]
    assert rows and all(r["status"] == "fixed_point" for r in rows), rows
    second = item_of(min(donor, key=lambda c: c.dp).spec, S(1, head_end="methyl",
                                                           tail_end="benzylic_radical"))
    cases = [c for c in win.member_cases(certifier, family, win.builds_of(second.specs),
                                         "forward")
             if not win.species_domain_refusal(c.witness) and
             any(("*1", "backbone_chph") in sites
                 for k, v in win._slot_sites(c, _mixed_index(second, mixed)).items()
                 for sites in v)]
    assert cases
    res = win.run_member(certifier, family, cases[0], bound=bound, contexts=False)
    assert res["status"] == "fixed_point", res


def _mixed_index(item, pattern):
    for k, spec in enumerate(item.specs):
        if isinstance(spec, b.PsChainSpec) and \
                sorted(e for e in (spec.head_end, spec.tail_end)) == sorted(pattern.ends):
            return k
    raise AssertionError(item)


# --------------------------------------------------------------------------
# rework 3 on the real families: per-reactant trees, closings, the
# full-vs-reduced comparison, the pre-filter
# --------------------------------------------------------------------------
def test_per_reactant_tree_and_joined_templates_are_read_from_the_loaded_family(rmg):
    """H_Abstraction and R_Addition_MultipleBond have one tree per reactant
    (two roots in ``groups.top``, not auto-generated): their forward
    direction is per-reactant.  Disproportionation and R_Recombination are
    auto-generated with one merged ``Root``: never.  One-product templates:
    R_Recombination forward (two radicals lost) and R_Addition_MultipleBond
    forward (net zero)."""
    fam, _certifier = rmg
    got = {name: win.family_constants(fam[name]) for name in
           ("H_Abstraction", "R_Addition_MultipleBond", "Disproportionation",
            "R_Recombination", "intra_H_migration")}
    assert got["H_Abstraction"].per_reactant == ("forward",)
    assert got["R_Addition_MultipleBond"].per_reactant == ("forward",)
    for name in ("Disproportionation", "R_Recombination", "intra_H_migration"):
        assert got[name].per_reactant == (), name
    for name in ("Disproportionation", "R_Recombination"):
        assert [e.label for e in fam[name].groups.top] == ["Root"]
        assert not win.per_reactant_tree(fam[name], "forward")
        assert not win.per_reactant_tree(fam[name], "reverse")
    assert len(fam["H_Abstraction"].groups.top) == 2
    assert not win.per_reactant_tree(fam["H_Abstraction"], "reverse")
    assert dict(got["R_Recombination"].joined) == {"forward": -2}
    assert dict(got["R_Addition_MultipleBond"].joined) == {"forward": 0}
    assert got["H_Abstraction"].joined == () and got["Disproportionation"].joined == ()


@pytest.fixture(scope="module")
def habs_closing(rmg):
    fam, certifier = rmg
    family = fam["H_Abstraction"]
    constants = win.family_constants(family)
    return family, constants, win.span_bound(constants), \
        win.closing_partners(certifier, family, constants)


def test_class_closings_merge_beyond_the_span_bound_and_split_leaves(rmg, habs_closing):
    """Real H_Abstraction closings against the closing partner: abstraction
    from the methyl end of two DP4 chains whose other ends differ (benzyl /
    styryl, seven backbone bonds away, beyond the span bound) merge; the
    methylene class of the same chain (leaf ``C/H2/NonDeC``, not
    ``C/H3/Cs``) does not; the windows widened to the span bound (``reach``)
    still merge, widened past the far end they do not."""
    fam, certifier = rmg
    family, constants, bound, partners = habs_closing
    assert partners["forward"][1][0] not in ("H", "methyl_radical")
    eb = P("chain", ("benzyl", "methyl"), ())
    st = P("chain", ("methyl", "styryl"), ())
    (spec_eb,), (spec_st,) = win.pattern_specs(eb, 4), win.pattern_specs(st, 4)
    methyl_eb = _class("forward", 0, (("*1", "end_methyl"),), eb, 4, spec_eb)
    ch2_eb = _class("forward", 0, (("*1", "backbone_ch2"),), eb, 4, spec_eb)
    methyl_st = _class("forward", 0, (("*1", "end_methyl"),), st, 4, spec_st)
    first, rows = win.close_classes(certifier, family, [methyl_eb, ch2_eb], partners["forward"][1],
                                    bound=bound)
    second, _rows = win.close_classes(certifier, family, [methyl_st], partners["forward"][1],
                                      bound=bound)
    by = {c.label: c for c in first + second}
    a, ch2, b_ = (by[win.class_label(x)] for x in (methyl_eb, ch2_eb, methyl_st))
    assert a.leaves == ("C/H3/Cs",) and b_.leaves == ("C/H3/Cs",)
    assert ch2.leaves == ("C/H2/NonDeC",)
    assert a.fixed and b_.fixed and ch2.fixed and all(r["status"] == "fixed_point" for r in rows)
    assert win.classes_merge(a, b_) and win.classes_merge(b_, a)
    assert not win.classes_merge(a, ch2)
    assert win.classes_merge(a, b_, reach=bound.reach)
    assert not win.classes_merge(a, b_, reach=8)
    merge = win.merge_classes([methyl_eb, ch2_eb, methyl_st], first + second)
    assert merge.rep(win.class_label(methyl_st)) == merge.rep(win.class_label(methyl_eb))
    assert merge.rep(win.class_label(ch2_eb)) == win.class_label(ch2_eb)


def _grid_classes(certifier, family, specs, slot, partner):
    """The classes ``specs`` carry in ``slot`` against ``partner`` (a probe on
    the given specs only; one pattern per spec, so labels stay unique)."""
    out = []
    for spec in specs:
        pattern = P("chain", (win._spec_name(spec), "grid"), ())
        for site in win.probe_slots(certifier, family, [spec, partner[1]],
                                    "forward").get(slot, ()):
            out.append(win.SiteClass("forward", slot, site, pattern, spec.units, spec,
                                     partner[0]))
    return out


def test_reduced_per_reactant_plan_covers_every_leaf_tuple_of_the_full_plan(rmg, habs_closing):
    """H_Abstraction (per-reactant) on a small grid (closed donors DP1-3,
    radicals DP1-2): the full plan is the product of every class spec per
    slot; its leaf tuples are the pathways of its supported in-domain
    members.  The reduced plan (closings, merge, one pair per leaf tuple)
    has fewer pairs, its closings carry exactly the full plan's leaves per
    slot, and its pairs reach a fixed point on every leaf tuple of the full
    plan."""
    from dataclasses import replace
    fam, certifier = rmg
    family, constants, bound, partners = habs_closing
    S = b.PsChainSpec
    donors = [S(1, head_end="methyl", tail_end="benzyl"), S(2, head_end="methyl", tail_end="benzyl"),
              S(3, head_end="methyl", tail_end="benzyl")]
    radicals = [S(1, head_end="methyl", tail_end="benzylic_radical"),
                S(2, head_end="primary_radical", tail_end="benzyl")]
    classes = _grid_classes(certifier, family, donors, 0, partners["forward"][1]) + \
        _grid_classes(certifier, family, radicals, 1, partners["forward"][0])
    full = set()
    for x in donors:
        for y in radicals:
            for case in win.member_cases(certifier, family, win.builds_of([x, y]), "forward"):
                if case.support == "supported" and not win.species_domain_refusal(case.witness):
                    full.add(case.pathway)
    assert len({p[0] for p in full}) >= 3 and len({p[1] for p in full}) >= 2
    groups = defaultdict(list)
    for c in classes:
        groups[(c.slot, repr(c.spec))].append(c)
    closures = []
    for (slot, _key), cs in sorted(groups.items()):
        got, _rows = win.close_classes(certifier, family, cs, partners["forward"][1 - slot],
                                       bound=bound)
        closures += got
    merge = win.merge_classes(classes, closures)
    for slot in (0, 1):
        closed = {leaf for c in closures if c.slot == slot for leaf in c.leaves}
        assert closed == {p[slot] for p in full}, (slot, closed)
    plan = win.ladder_plan(replace(constants, pair_kinds=(("forward", ("closed", "radical")),),
                                   classes=tuple(classes), merge=merge))
    leaf = [i for i in plan.items if "pair:leaf" in i.kinds]
    assert leaf and len(leaf) < len(donors) * len(radicals)
    assert not [i for i in plan.items if "pair:class" in i.kinds]
    assert {tuple(t) for i in leaf for t in i.leaf_tuples} >= full
    reached = set()
    for item in leaf:
        todo = {tuple(t) for t in item.leaf_tuples} - reached
        for case in win.member_cases(certifier, family, win.builds_of(item.specs), "forward"):
            if case.pathway in todo and case.support == "supported" and \
                    not win.species_domain_refusal(case.witness):
                res = win.run_member(certifier, family, case, bound=bound, contexts=False)
                if res["status"] == "fixed_point":
                    reached.add(case.pathway)
                    todo.discard(case.pathway)
    assert full <= reached, sorted(full - reached)
    print("GRID full pairs", len(donors) * len(radicals), "reduced pairs", len(leaf),
          "classes", len(classes), "representatives",
          len({merge.rep(win.class_label(c)) for c in classes}), "leaf tuples", sorted(full))
    assert win.leaf_tuple_gaps(plan, [{"status": "fixed_point", "direction": "forward",
                                       "pathway": list(p)} for p in reached]) == []


def test_prefilter_never_drops_a_pair_with_a_fixed_point_member(rmg):
    """R_Addition_MultipleBond forward (one product): on a small grid of
    alkenes (styrene; a DP2 chain with a vinylidene end and a benzylic end
    radical) and radicals (methyl; a DP3 biradical), every pair the
    pre-filter drops has no supported in-domain member (so none can reach a
    fixed point), and the pairs it keeps at exactly MAX_RADICAL_SITES
    product electrons have a fixed-point member."""
    fam, certifier = rmg
    family = fam["R_Addition_MultipleBond"]
    constants = win.family_constants(family)
    bound = win.span_bound(constants)
    S = b.PsChainSpec
    styrene = S(1, head_end="vinylidene", tail_end="styryl")
    alkene_radical, = win.pattern_specs(P("chain", ("benzylic_radical", "vinylidene"), ()), 2)
    methyl = win.ONE_CARBON["methyl_radical"]()
    biradical = S(3, head_end="primary_radical", tail_end="benzylic_radical")
    dropped, at_max = 0, []
    for x in (styrene, alkene_radical):
        for y in (methyl, biradical):
            why = win.pair_prefilter(constants, "forward", x, y)
            cases = [c for c in win.member_cases(certifier, family, win.builds_of([x, y]),
                                                 "forward")
                     if c.support == "supported" and not win.species_domain_refusal(c.witness)]
            if why:
                dropped += 1
                assert not cases, (why, [c.label for c in cases])
                continue
            if win.unpaired_electrons(x) + win.unpaired_electrons(y) == dom.MAX_RADICAL_SITES:
                at_max.append(any(win.run_member(certifier, family, c, bound=bound,
                                                 contexts=False)["status"] == "fixed_point"
                                  for c in cases[:4]))
    assert dropped == 1 and len(at_max) == 2 and any(at_max), at_max


# --------------------------------------------------------------------------
# slow: the required supported cases and the R_Recombination ladder control
# --------------------------------------------------------------------------
@pytest.mark.slow
@pytest.mark.skipif(not SLOW, reason="set RMG_KMC_SLOW=1")
def test_required_supported_cases_reach_recertified_fixed_points(rmg):
    """Every selected member of every required case reaches a fixed point (a
    reacting atom on an artificial port is accounted by its status), every
    mapping's edits lie in the window on both sides, the span bound holds,
    the A+A recombination member covers both template slot orders."""
    import windows_acceptance as wa
    fam, certifier = rmg
    reached = {}
    for label, family_name, direction, specs, selection in wa.required_sets():
        family = fam[family_name]
        bound = win.span_bound(win.family_constants(family))
        cases = [c for c in win.member_cases(certifier, family, wa.builds_of(specs), direction)
                 if wa.select(selection, c)]
        assert cases, label
        rows = [win.run_member(certifier, family, c, bound=bound, contexts=False)
                for c in cases]
        statuses = {r["status"] for r in rows}
        assert statuses <= {"fixed_point", "refused:reacting_port"}, (label, statuses)
        for r in rows:
            if r["status"] == "fixed_point":
                assert not r["mapping_problems"] and r["span_bound_held"], (label, r)
        reached[label] = sum(r["status"] == "fixed_point" for r in rows)
        if "identical" in label:
            assert any(len(r["slot_orders"]) == 2 for r in rows)
        print("REQUIRED", label, len(rows), reached[label])
    assert all(reached.values()), reached
    assert len(reached) == 8


@pytest.mark.slow
@pytest.mark.skipif(not SLOW, reason="set RMG_KMC_SLOW=1")
def test_recombination_ladder_control(rmg):
    """Homolysis of the intact closed chains of DP L and L + 1 (R_Recombination):
    no window of DP L + 1 holds both physical ends (non-overlap control)."""
    fam, certifier = rmg
    family = fam["R_Recombination"]
    bound = win.span_bound(win.family_constants(family))
    both = {}
    for dp in (bound.L, bound.L + 1):
        spec = b.PsChainSpec(dp, head_end="methyl", tail_end="benzyl")
        build = b.PROVIDERS["polystyrene"].build("ladder", spec)
        rows = [win.run_member(certifier, family, c, bound=bound, contexts=False)
                for c in win.member_cases(certifier, family, [build], "reverse")]
        reached = [r for r in rows if r["status"] == "fixed_point"]
        assert reached
        both[dp] = sum(any(s["both_ends"] for s in r["span"].values()) for r in reached)
    print("LADDER CONTROL", bound.L, both)
    assert both[bound.L + 1] == 0
