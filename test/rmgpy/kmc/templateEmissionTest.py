"""Tests for rmgpy/kmc/template_emission.py (I-138): executable local window
templates cut from recertified closed windows.

Fast part, no database (toy alkane homolysis on the linear-alkane provider,
a graph oracle for the channel count, stub production rates): the cut
(closed atoms, ports, edits never on a port, construction-carried port
labels under permuted atom ids), end and interior templates, the count and
weight (count mismatch, T_ch = 0, whole-channel normalization and the
calibration identity), locality over extension vectors (stub rates moving at
one port step only), upstream refusals, the schema round trip and
determinism (two subprocesses, two hash seeds, permuted atoms).
Fast part with the RMG database (``channelsTest.database``): the
forbidden-move exclusion on the H-atom abstraction reference, the DP3
chain-to-chain transfer produced, and the bimolecular pair counts on real
windows (identical and unequal radicals).
Slow part (``RMG_KMC_SLOW=1``): ``template_emission_acceptance.py`` on the
reference compiles.
"""

import json
import math
import os
import subprocess
import sys
from collections import defaultdict, namedtuple
from pathlib import Path

import pytest

from rmgpy.kmc import boundaries as b
from rmgpy.kmc import rmg_moves as rm
from rmgpy.kmc import template_emission as te
from rmgpy.kmc import window_matcher as wm
from rmgpy.kmc import windows as win

SLOW = os.environ.get("RMG_KMC_SLOW") == "1"
HERE = Path(__file__).resolve().parent
GRID = tuple(float(t) for t in range(300, 1001, 25))
K_RULE = tuple(1.0e9 * math.exp(-15000.0 / t) for t in GRID)
HEAD, TAIL = (0, "head"), (0, "tail")


# --------------------------------------------------------------------------
# toy chemistry: C-C homolysis / recombination of linear alkanes
# --------------------------------------------------------------------------
class ToyFamily:
    label = "toy_cc"
    reverse_template = object()

    class forward_template:
        reactants = (0, 1)

    class reverse_template_:
        reactants = (0,)


def _break(rows, a, c):
    rows = dict(rows)
    for x, y in ((a, c), (c, a)):
        row = rows[x]
        rows[x] = row[:3] + (row[3] + 1,) + row[4:5] + (tuple(e for e in row[5] if e[0] != y),)
    return rows


def _join(rows, a, c):
    rows = dict(rows)
    for x, y in ((a, c), (c, a)):
        row = rows[x]
        rows[x] = row[:3] + (row[3] - 1,) + row[4:5] + (tuple(sorted(row[5] + ((y, 1.0),))),)
    return rows


class ToyCounter(te.RmgCounter):
    """The graph oracle of the toy family: homolysis of every C-C bond
    (reverse; both orientations are directed mappings) and recombination of
    every radical pair across two molecules (forward; both slot orders).
    ``extra`` adds phantom outcomes (the count-mismatch fixture)."""

    def __init__(self, extra=0):
        from rmgpy.kmc import channels as chm
        self.species = chm.SpeciesKeys()
        self._product_species = {}
        self.extra = extra
        self.calls = 0

    def family(self, name):
        return ToyFamily

    def channel(self, family, molecules, direction, pathway, product_species):
        self.calls += 1
        rows = {}
        for m in molecules:
            rows.update(rm.molecule_rows(m))
        outcomes = []
        if direction == "reverse":
            for a, row in rows.items():
                for c, _o in row[5]:
                    if a < c and row[0] == "C" and rows[c][0] == "C":
                        outcomes.append((_break(rows, a, c), ((0,), (0,))))
        else:
            parts = [rm.molecule_rows(m) for m in molecules]
            rad = [[a for a, r in p.items() if r[3]] for p in parts]
            if len(parts) == 2:
                for a in rad[0]:
                    for c in rad[1]:
                        outcomes.append((_join(rows, a, c), ((0, 1), (1, 0))))
        keys, orders, directed = [], [], 0
        for prod, slot in outcomes:
            if self.species_of_rows(prod) != tuple(product_species):
                continue
            directed += 2
            keys.append(rm.product_key(prod))
            orders.append(slot)
        keys = sorted(set(keys)) + [("phantom", k) for k in range(self.extra)]
        return te.RmgChannel(tuple(keys), directed, tuple(orders), (("kept", directed),))


def alkane(n, head=None, tail=None):
    build = b.PROVIDERS["linear_alkane"].build(n, head_end=head, tail_end=tail)
    mols, parts = b.witness_reactants([build])
    return rm.molecule_rows(mols[0]), parts


def homolysis(n, bond, head=None, tail=None):
    rows, parts = alkane(n, head, tail)
    return b.make_witness(tuple(sorted(rows.items())), tuple(sorted(_break(rows, *bond).items())),
                          parts, label=f"C{n}|{bond}")


def recombination(n1, n2, head1=None):
    """Two primary-radical-ended alkanes (tail physical radical, head
    artificial unless ``head1``) recombining at their tails."""
    b1 = b.PROVIDERS["linear_alkane"].build(n1, head_end=head1, tail_end="primary_radical")
    b2 = b.PROVIDERS["linear_alkane"].build(n2, tail_end="primary_radical")
    mols, parts = b.witness_reactants([b1, b2])
    rows = {}
    for m in mols:
        rows.update(rm.molecule_rows(m))
    a, c = n1 - 1, len(b1.construction_order) + n2 - 1
    return b.make_witness(tuple(sorted(rows.items())), tuple(sorted(_join(rows, a, c).items())),
                          parts, label=f"C{n1}+C{n2}")


def heavy_ball(witness, radius):
    rows = witness.reactant_rows()
    support = b.edit_support(witness)
    seen = {a for a in support if rows[a][0] != "H"}
    frontier = list(seen)
    for _ in range(radius):
        nxt = []
        for a in frontier:
            for n, _o in rows[a][5]:
                if rows[n][0] != "H" and n not in seen:
                    seen.add(n)
                    nxt.append(n)
        frontier = nxt
    return seen


def window(witness, radius=2):
    support = win.Support(atoms=frozenset(heavy_ball(witness, radius)),
                          reacting=frozenset(b.edit_support(witness)))
    return win.close(witness, support)


def cut(witness, radius=2, direction="reverse", reverse=False):
    return te.cut_template(witness, window(witness, radius).atoms, family="toy_cc",
                           pathway=("toy",), direction=direction, reverse=reverse)


Move = namedtuple("Move", "labels product_labels")


def member(witness, radius=2, direction="reverse"):
    w = window(witness, radius)
    pair = sorted(a for a in b.edit_support(witness) if witness.reactant_rows()[a][0] == "C")
    move = Move(((0, "*1", pair[0]), (0, "*2", pair[1])), (("*1", pair[0]), ("*2", pair[1])))
    case = win.MemberCase(f"toy|{witness.label}", "toy_cc", direction, ("toy",), witness,
                          (move,), frozenset())
    ctx = win.ContextResult(tuple((p, None) for p in b.vector_ports(witness)), "fixed_point", "",
                            witness, w, None)
    return case, {"status": "fixed_point"}, ((ctx,), ())


def graph(rows):
    return [m.to_adjacency_list() for m in te.molecules_of(rows)]


def toy_record(witness, transitions, *, status="enabled", orientation="reversed"):
    """The record of the witness's homolysis channel: k = k_rule * (channel
    members), ssa 1."""
    r, p = dict(witness.reactant_rows()), dict(witness.product_rows())
    return {"event_id": "evt_" + "0" * 64, "family": "toy_cc", "template": "toy",
            "orientation": orientation, "status": status,
            "status_reason": "" if status == "enabled" else "domain",
            "reactant_graphs": graph(r), "product_graphs": graph(p),
            "reactant_multiplicities": [1], "ssa_multiplier": 1.0,
            "k_table": {"T": list(GRID), "k": [k * transitions for k in K_RULE],
                        "interpolation": "linear-ln-k", "extrapolation": "refuse"},
            "rate_source": {"kind": "RMG family estimate", "source": "rate rules",
                            "entry": "toy", "rank": 11}}


def stub_rates(counter, shift=lambda req: 0.0):
    """Production k on a witness: k_rule times the RMG (oracle) degeneracy of
    the witness's channel, moved by ``shift(request)`` decades."""
    def rates(req):
        species = counter.species_of_rows(req.products)
        deg = counter.channel(ToyFamily, te.molecules_of(req.generation), "reverse", None,
                              species).transitions
        d = shift(req)
        return te.RateResult(GRID, tuple(k * deg * 10 ** d for k in K_RULE), 1.0,
                             {"selection": {"kind": "stub"}}, ())
    return rates


def emit(witness, *, radius=2, counter=None, shift=lambda req: 0.0, record=None, **kw):
    counter = counter or ToyCounter()
    case, row, ctx = member(witness, radius)
    if record is None:
        oracle = ToyCounter().channel(ToyFamily, te.molecules_of(dict(witness.reactant_rows())),
                                      "reverse", None,
                                      ToyCounter().species_of_rows(witness.product_rows()))
        record = toy_record(witness, oracle.transitions)
    return te.emit_record(record, family=ToyFamily, members=[(case, row, ctx)], counter=counter,
                          rates=kw.pop("rates", None) or stub_rates(ToyCounter(), shift),
                          check_upstream=kw.pop("check_upstream", False), **kw)


def base_template(payload):
    return next(t for t in payload["templates"] if t["base"])


# --------------------------------------------------------------------------
# item 1: the cut
# --------------------------------------------------------------------------
def test_cut_closed_atoms_ports_and_edits():
    w = homolysis(8, (3, 4))
    c = cut(w)
    rows = w.reactant_rows()
    heavy = sorted(c.witness_atoms[i] for i, a in enumerate(c.atoms) if a[1] == "C")
    assert heavy == [1, 2, 3, 4, 5, 6]
    for i, a in enumerate(c.atoms):
        x = c.witness_atoms[i]
        closed = all(n in set(c.witness_atoms) for n, _o in rows[x][5])
        assert (a[6] == "closed") == closed
        if a[6] == "port":
            assert a[7] == ((1.0, "C", False),)
    edited = {c.witness_atoms[i] for i in c.edited()}
    assert edited == {3, 4}
    assert all(c.atoms[i][6] == "closed" for i in c.edited())
    t = c.matcher(1.0)
    assert len(t.edits) == 3      # the broken bond, two radicals
    assert c.participants == 1 and c.participant_connected()
    assert te.template_refusal(c) is None
    assert set(c.boundary_labels().values()) == {"artificial"}


def test_edit_on_a_port_is_an_error():
    w = homolysis(8, (3, 4))
    atoms = {3, 4} | {h for a in (3, 4) for h, _o in w.reactant_rows()[a][5]
                      if w.reactant_rows()[h][0] == "H"}
    with pytest.raises(te.EmissionError, match="edit touches the window port"):
        te.cut_template(w, atoms, family="toy_cc", pathway=("toy",), direction="reverse")


def _relabel(witness, perm):
    rows_r = {perm[a]: r[:5] + (tuple(sorted((perm[n], o) for n, o in r[5])),)
              for a, r in witness.reactant_rows().items()}
    rows_p = {perm[a]: r[:5] + (tuple(sorted((perm[n], o) for n, o in r[5])),)
              for a, r in witness.product_rows().items()}
    parts = [b.Participant(p.index, p.metadata.relabel(perm), tuple(perm[a] for a in p.atom_ids))
             for p in witness.participants]
    return b.make_witness(tuple(sorted(rows_r.items())), tuple(sorted(rows_p.items())), parts,
                          label=witness.label)


def _perm(witness, seed=7):
    import random
    ids = sorted(witness.reactant_rows())
    shuffled = list(ids)
    random.Random(seed).shuffle(shuffled)
    return {a: 1000 + x for a, x in zip(ids, shuffled)}


def _anchors(witness, ids):
    """Per mapped atom: itself (heavy) or its heavy neighbour (a hydrogen:
    interchangeable H leaves map up to their automorphism)."""
    rows = witness.reactant_rows()
    return [a if rows[a][0] != "H" else ("H", rows[a][5][0][0]) for a in ids]


def test_port_labels_come_from_construction_under_permuted_atoms():
    """C10 with a PHYSICAL methyl head and an artificial tail: the head-side
    port reaches only the physical end (``physical``: refused as a
    whole-species guard), the tail-side one the artificial tail.  Relabelling
    every atom id changes neither the labels nor the executable key, and the
    witness map follows the permutation."""
    w = homolysis(10, (4, 5), head="methyl")
    c = cut(w)
    labels = {c.witness_atoms[i]: lab for i, lab in c.boundary_labels().items()}
    assert labels == {2: "physical", 7: "artificial"}
    assert te.template_refusal(c).startswith("refused:whole_species_guard")
    perm = _perm(w)
    pw = _relabel(w, perm)
    pc = te.cut_template(pw, {perm[a] for a in window(w).atoms}, family="toy_cc",
                         pathway=("toy",), direction="reverse")
    assert pc.key == c.key and pc.meta == c.meta and pc.atoms == c.atoms
    assert _anchors(pw, pc.witness_atoms) == _anchors(pw, [perm[a] for a in c.witness_atoms])
    # the same window with an artificial head: both ports artificial
    a = cut(homolysis(10, (4, 5)))
    assert set(a.boundary_labels().values()) == {"artificial"}
    assert te.template_refusal(a) is None


def test_end_template_does_not_match_mid_chain():
    """The window of the bond next to a physical methyl end holds the end
    (closed CH3): on a long chain it fires only beside a terminus."""
    w = homolysis(10, (1, 2), head="methyl")
    c = cut(w)
    end = c.witness_atoms.index(0)
    assert c.atoms[end][6] == "closed" and c.meta[end][2] == ("head", "backbone_head", "physical")
    rows, _parts = alkane(16, head="methyl", tail="methyl")
    state = te.melt_state(rows)
    result = wm.match([c.matcher(1.0)], state)
    broken = sorted(tuple(sorted((e.a, e.b))) for tr in result.transitions for e in tr.edits
                    if isinstance(e, wm.SetBond))
    assert broken == [(1, 2), (13, 14)]


def test_interior_template_does_not_match_at_a_physical_terminus():
    c = cut(homolysis(12, (5, 6)))
    rows, _parts = alkane(10, head="methyl", tail="methyl")
    state = te.melt_state(rows)
    result = wm.match([c.matcher(1.0)], state)
    broken = sorted(tuple(sorted((e.a, e.b))) for tr in result.transitions for e in tr.edits
                    if isinstance(e, wm.SetBond))
    assert broken == [(3, 4), (4, 5), (5, 6)]   # never a bond within three of an end


def test_disconnected_participant_window_is_refused_by_name():
    w = homolysis(12, (5, 6))
    atoms = set(window(w).atoms) | {0} | {h for h, _o in w.reactant_rows()[0][5]
                                          if w.reactant_rows()[h][0] == "H"}
    c = te.cut_template(w, atoms, family="toy_cc", pathway=("toy",), direction="reverse")
    assert te.template_refusal(c).startswith("refused:path_predicate")


# --------------------------------------------------------------------------
# item 2: count and weight
# --------------------------------------------------------------------------
def test_whole_channel_weight_and_calibration():
    """C10 bond 3-4: its mirror bond 5-6 gives the same species, so the channel
    has two members, T_R = 2, w = k / 2, and the calibration identity holds
    at 600/700/800 K."""
    w = homolysis(10, (3, 4))
    out = emit(w)
    assert out["status"] == "emitted", out.get("refusal")
    assert out["counts"]["record"]["matcher"] == out["counts"]["record"]["rmg"] == 2
    k = out["rate"]["k"]
    assert out["weight"]["w"] == pytest.approx([x / 2 for x in k], rel=1e-15)
    assert [c["T"] for c in out["calibration"]] == [600.0, 700.0, 800.0]
    assert all(c["ok"] and c["relative_error"] <= 1e-9 for c in out["calibration"])
    t = base_template(out)
    assert t["status"] == "emitted" and t["locality"]["accepted"]
    assert t["counts"]["emitted"]["matcher"] == t["counts"]["emitted"]["rmg"] == 2


def test_count_mismatch_is_refused_with_both_numbers():
    w = homolysis(10, (3, 4))
    out = emit(w, counter=ToyCounter(extra=1))
    assert out["status"] == "refused"
    assert out["refusal"]["status"] == "refused:count_mismatch"
    assert (out["refusal"]["matcher"], out["refusal"]["rmg"]) == (2, 3)
    assert "matcher 2" in out["refusal"]["reason"] and "rmg_moves 3" in out["refusal"]["reason"]


def test_zero_transitions_is_a_compile_error():
    """The interior template of C10 bond 3-4 under a record whose channel is
    the methyl loss of C10: RMG has two members, the template none."""
    w = homolysis(10, (3, 4))
    record = toy_record(homolysis(10, (0, 1)), 2)
    with pytest.raises(te.EmissionError, match="T_ch = 0"):
        emit(w, record=record)
    with pytest.raises(te.EmissionError, match="T_ch = 0"):
        te.weight_table(GRID, K_RULE, 1.0, 0)
    with pytest.raises(te.EmissionError, match="non-finite or non-positive"):
        te.weight_table(GRID, [0.0] * len(GRID), 1.0, 1)


# --------------------------------------------------------------------------
# item 3: locality over extension vectors
# --------------------------------------------------------------------------
def _shift_at(vector, decades):
    return lambda req: decades if tuple(req.vector) == vector and req.label != "context" else 0.0


@pytest.mark.parametrize("vector,decades,accepted", [
    ((0, 2), 0.02, False),     # only at tail +2
    ((0, 1), 0.02, False),     # only at tail +1
    ((0, 2), 0.005, True),
    ((2, 0), 0.02, False),     # only at head +2
])
def test_locality_spread_per_port_step(vector, decades, accepted):
    out = emit(homolysis(10, (3, 4)), shift=_shift_at(vector, decades))
    t = base_template(out)
    assert t["locality"]["accepted"] is accepted
    if not accepted:
        assert t["refusal"]["status"] == "refused:non_local_kinetics"
        assert t["locality"]["max_spread_log10"] == pytest.approx(decades, rel=1e-9)
    labels = [r["label"] for r in t["locality"]["rows"]]
    for want in ("record", "context", "emitted", "port 0:head +1", "port 0:head +2",
                 "port 0:tail +1", "port 0:tail +2", "joint +1", "joint +2"):
        assert want in labels, (want, labels)


def test_locality_needs_the_emitted_witness():
    shift = lambda req: 0.02 if req.label == "emitted" else 0.0
    t = base_template(emit(homolysis(10, (3, 4)), shift=shift))
    assert t["refusal"]["status"] == "refused:non_local_kinetics"


def test_locality_recounts_the_channel_on_every_witness():
    """C8 central bond: one member; one carbon more makes its mirror bond a
    second member of the same channel (k doubles with RMG's degeneracy).  The
    weight is local only when T_ch is recounted on every witness."""
    w = homolysis(8, (3, 4))
    out = emit(w)
    t = base_template(out)
    assert t["status"] == "emitted", t.get("refusal")
    rows = {r["label"]: r for r in t["locality"]["rows"]}
    assert rows["emitted"]["counts"]["matcher"] == 1
    assert rows["port 0:tail +1"]["counts"]["matcher"] == 2
    assert rows["port 0:tail +1"]["counts"]["rmg"] == 2


def test_locality_witnesses_are_distinct_marked_graphs():
    w = homolysis(10, (3, 4))
    ws = te.locality_witnesses(w)
    labels = [x.label for x in ws]
    assert labels[0] == "emitted"
    assert "participant 0 +1" not in labels       # equal to the joint vector: folded
    assert not b.marked_duplicates([(x.label, x.witness) for x in ws if x.status == "ok"])


# --------------------------------------------------------------------------
# upstream refusals, schema, determinism
# --------------------------------------------------------------------------
def test_template_reaching_the_padding_weighs_on_its_emitted_witness():
    """A context whose window sits on the padded witness (C8 (2,3), both
    artificial ends extended by one unit) with a port exterior on the
    padding: the template cannot apply to the proxy record, so the record
    and unextended-context rows are not applicable (never a zero count) and
    the weight and calibration come from the emitted witness."""
    w0 = homolysis(8, (2, 3))
    n = len(list(b.vector_ports(w0)))
    padded = b.evaluate_vector(w0, (1,) * n).witness
    win_ = window(padded, 2)
    case, row, _ctx = member(w0, 2)
    ctx = win.ContextResult(tuple((p, None) for p in b.vector_ports(w0)), "fixed_point", "",
                            padded, win_, None)
    oracle = ToyCounter().channel(ToyFamily, te.molecules_of(dict(w0.reactant_rows())),
                                  "reverse", None, ToyCounter().species_of_rows(w0.product_rows()))
    record = toy_record(w0, oracle.transitions)
    cut_ = te.cut_template(padded, win_.atoms, family="toy_cc", pathway=("toy",),
                           direction="reverse")
    assert te.reaches_beyond(cut_, padded, w0)
    assert not te.reaches_beyond(cut(w0), w0, w0)
    out = te.emit_record(record, family=ToyFamily, members=[(case, row, ((ctx,), ()))],
                         counter=ToyCounter(), rates=stub_rates(ToyCounter()),
                         check_upstream=False)
    assert out["status"] == "emitted", out.get("refusal")
    assert out["weight_source"] == "emitted_witness" and out["counts"]["record"] is None
    t = base_template(out)
    rows = {r["label"]: r for r in t["locality"]["rows"]}
    assert rows["record"]["status"] == rows["context"]["status"] == "not_applicable"
    assert t["weight_source"] == "emitted_witness"
    assert t["weight"]["w"] == rows["emitted"]["w"]
    assert [c["ok"] for c in t["calibration"]] == [True] * 3
    te.validate_payload(out)


def test_padding_reach_through_a_port_exterior_alone():
    """The padded C8 (2,3) window holds only original atoms: the template
    reaches the padding through a port's exterior neighbour alone, so the
    port-neighbour reach is load-bearing.  A port is never edited, so its
    neighbours are the same on both witness sides."""
    w0 = homolysis(8, (2, 3))
    padded = b.evaluate_vector(w0, (1,) * len(list(b.vector_ports(w0)))).witness
    cut_ = te.cut_template(padded, window(padded, 2).atoms, family="toy_cc",
                           pathway=("toy",), direction="reverse")
    assert set(cut_.witness_atoms) <= set(dict(w0.reactant_rows()))
    assert te.reaches_beyond(cut_, padded, w0)


def test_a_refused_emitted_witness_leaves_the_template_without_a_weight():
    """Weight from the emitted witness (the template reaches the padding) but
    the emitted witness's own row is refused (its rate is not rebound): the
    template is refused and carries no weight and no calibration."""
    w0 = homolysis(8, (2, 3))
    padded = b.evaluate_vector(w0, (1,) * len(list(b.vector_ports(w0)))).witness
    win_ = window(padded, 2)
    case, row, _ctx = member(w0, 2)
    ctx = win.ContextResult(tuple((p, None) for p in b.vector_ports(w0)), "fixed_point", "",
                            padded, win_, None)
    oracle = ToyCounter().channel(ToyFamily, te.molecules_of(dict(w0.reactant_rows())),
                                  "reverse", None, ToyCounter().species_of_rows(w0.product_rows()))
    base = stub_rates(ToyCounter())

    def rates(req):
        res = base(req)
        if req.label == "emitted":
            return te.RateResult(res.temperatures, res.k, res.ssa_multiplier, res.provenance,
                                 ("stub: the emitted witness's selection is not rebound",))
        return res

    out = te.emit_record(toy_record(w0, oracle.transitions), family=ToyFamily,
                         members=[(case, row, ((ctx,), ()))], counter=ToyCounter(),
                         rates=rates, check_upstream=False)
    t = base_template(out)
    rows = {r["label"]: r for r in t["locality"]["rows"]}
    assert rows["emitted"]["status"] == "refused:rebind"
    assert t["status"] == "refused" and t["weight"] is None and "calibration" not in t
    assert out["status"] == "refused"
    te.validate_payload(out)


# --------------------------------------------------------------------------
# catalog-level collisions
# --------------------------------------------------------------------------
def _toy_catalog():
    """C12 homolysis at (3,4) and at (4,5): two records (two channels), one
    local template (radius 2, artificial ends far)."""
    out = []
    for k, bond in enumerate([(3, 4), (4, 5)]):
        w = homolysis(12, bond)
        oracle = ToyCounter().channel(ToyFamily, te.molecules_of(dict(w.reactant_rows())),
                                      "reverse", None,
                                      ToyCounter().species_of_rows(w.product_rows()))
        rec = toy_record(w, oracle.transitions)
        rec["event_id"] = "evt_" + str(k + 1) * 64
        out.append(emit(w, record=rec))
    return out


def _per_outcome(templates, rows):
    state = te.melt_state(rows)
    counts = defaultdict(int)
    for tr in wm.match(templates, state).transitions:
        counts[rm.product_key(wm.product_rows(state, tr))] += 1
    return counts


def test_equal_weight_collision_is_one_template():
    """Two records emit the same template with equal weights: loaded together
    every matched cleavage fires twice; after the catalog collision gate it
    fires once, and the surviving template lists both source records."""
    payloads = _toy_catalog()
    keys = [base_template(p)["key"] for p in payloads]
    assert keys[0] == keys[1]
    rows = alkane(12)[0]
    before = _per_outcome(te.loaded_templates(payloads, 700.0), rows)
    assert before and set(before.values()) == {2}
    merged, table = te.resolve_collisions(payloads)
    after = _per_outcome(te.loaded_templates(merged, 700.0), rows)
    assert set(after) == set(before) and set(after.values()) == {1}
    assert [r["disposition"] for r in table] == ["merged"] and table[0]["copies"] == 2
    # the representative is the copy of the least partition-free identity,
    # never the least record id, and whatever the input order
    by_identity = [p["record_event_id"] for p in sorted(payloads, key=lambda p: te._template_identity(
        te.partition_free_identity(p), p, base_template(p)))]
    assert sorted(by_identity) == ["evt_" + "1" * 64, "evt_" + "2" * 64]
    assert te.resolve_collisions(payloads[::-1])[1] == table
    by_record = {p["record_event_id"]: base_template(p) for p in merged}
    owner, other = by_record[by_identity[0]], by_record[by_identity[1]]
    assert owner["status"] == "emitted"
    assert owner["collision"]["records"] == by_identity
    assert other["status"] == "merged" and other["merged_into"]["record"] == by_identity[0]
    assert [p["status"] for p in merged] == ["emitted", "emitted"]
    for p in merged:
        te.validate_payload(p)


def test_unequal_weight_collision_is_refused():
    """The same template with weights 10 % (0.041 decades) apart: no merge,
    every copy and its record refused ``refused:collision`` with keys,
    records and weights; nothing of the key is loaded."""
    payloads = _toy_catalog()
    t = base_template(payloads[1])
    t["weight"]["w"] = [w * 1.1 for w in t["weight"]["w"]]
    payloads[1] = te.finish(payloads[1])
    out, table = te.resolve_collisions(payloads)
    assert [r["disposition"] for r in table] == ["refused:collision"]
    assert table[0]["log10_diameter"] == pytest.approx(math.log10(1.1), rel=1e-9)
    # copies are listed by partition-free identity, never by record id (the
    # ids carry the artifact's family list), and whatever the input order
    by_identity = [p["record_event_id"] for p in sorted(payloads, key=lambda p: te._template_identity(
        te.partition_free_identity(p), p, base_template(p)))]
    assert sorted(by_identity) == ["evt_" + "1" * 64, "evt_" + "2" * 64]
    assert te.resolve_collisions(payloads[::-1])[1] == table
    for p in out:
        t = base_template(p)
        assert t["status"] == "refused" and t["refusal"]["status"] == "refused:collision"
        assert t["refusal"]["keys"] == [t["key"]] and len(t["refusal"]["weights"]) == 2
        assert t["refusal"]["records"] == by_identity
        assert p["status"] == "refused" and p["refusal"]["status"] == "refused:collision"
        te.validate_payload(p)
    assert te.loaded_templates(out, 700.0) == []


def test_a_payload_attaches_only_to_its_own_record():
    """A payload emitted for one record cannot ride on another (another event
    id, family or orientation): attach and EventRecord.validate refuse it."""
    from dataclasses import replace
    w = homolysis(10, (3, 4))
    payload = emit(w)
    rec = _event_record(w)
    payload = te.finish(dict(payload, record_event_id=rec.event_id))
    te.attach(rec, payload).validate()
    other = te.finish(dict(payload, record_event_id="evt_" + "9" * 64))
    with pytest.raises(te.EmissionError, match="belongs to"):
        te.attach(rec, other)
    flipped = te.finish(dict(payload, orientation="as_generated"))
    with pytest.raises(te.EmissionError, match="belongs to"):
        te.attach(rec, flipped)
    smuggled = replace(rec, schema_version=te.RECORD_SCHEMA, executable_template=other)
    with pytest.raises(ValueError, match="belongs to"):
        smuggled.validate()


def test_upstream_refused_records_stay_refused_with_their_status():
    w = homolysis(10, (3, 4))
    rec = toy_record(w, 2, status="refused")
    out = emit(w, record=rec, check_upstream=True, record_certificate={"certified": True})
    assert out["status"] == "refused" and out["templates"] == []
    assert out["refusal"]["status"] == "refused:record_status"
    rec = toy_record(w, 2)
    out = emit(w, record=rec, check_upstream=True,
               record_certificate={"certified": False, "gaps": ["binding: no member"]})
    assert out["refusal"]["status"] == "refused:certificate" and out["templates"] == []
    case, _row, ctx = member(w)
    out = te.emit_record(rec, family=ToyFamily,
                         members=[(case, {"status": "refused:reacting_port", "covered": True},
                                   ctx)],
                         counter=ToyCounter(), rates=stub_rates(ToyCounter()),
                         check_upstream=False)
    assert out["refusal"]["status"] == "refused:window" and out["templates"] == []
    assert out["refusal"]["window_status"] == "refused:reacting_port"


def _event_record(witness):
    from rmgpy.kmc.event_record import EventRecord
    rec = toy_record(witness, 2)
    return EventRecord(family="toy_cc", template="toy", arity=1,
                       participant_site_types=["ref0"], reactant_multiplicities=[1],
                       k_table=rec["k_table"], reactant_graphs=rec["reactant_graphs"],
                       product_graphs=rec["product_graphs"], rate_source=rec["rate_source"],
                       orientation="reversed")


def test_schema_round_trip_and_hashing():
    from rmgpy.kmc.event_record import EventRecord
    w = homolysis(10, (3, 4))
    old = _event_record(w)
    payload = te.finish(dict(emit(w), record_event_id=old.event_id))
    te.validate_payload(payload)
    old.validate()
    assert "executable_template" not in old.to_dict()
    assert old.to_dict()["schema_version"] == "rewrite_record/0.1"
    new = te.attach(old, payload)
    assert new.schema_version == "rewrite_record/0.2"
    assert new.event_id == old.event_id          # the template is not part of the id
    new.validate()
    again = EventRecord.from_json(new.to_json())
    again.validate()
    assert again.executable_template == payload and again.event_id == old.event_id
    assert EventRecord.from_json(old.to_json()).to_dict() == old.to_dict()
    t = base_template(payload)
    m = te.matcher_from_payload(dict(t, family="toy_cc"), 650.0)
    assert isinstance(m, wm.Template) and m.weight > 0
    lo, hi = t["weight"]["w"][GRID.index(625.0)], t["weight"]["w"][GRID.index(650.0)]
    assert m.weight == pytest.approx(hi)
    with pytest.raises(te.EmissionError, match="extrapolation"):
        te.matcher_from_payload(t, 1200.0)
    tampered = json.loads(json.dumps(payload))
    base_template(tampered)["weight"]["w"][0] *= 2
    with pytest.raises(te.EmissionError, match="sha256"):
        te.validate_payload(tampered)
    bad = json.loads(json.dumps(new.to_dict()))
    bad["executable_template"] = tampered
    with pytest.raises(ValueError, match="executable_template"):
        EventRecord.from_dict(bad).validate()
    lone = json.loads(json.dumps(old.to_dict()))
    lone["schema_version"] = "rewrite_record/0.2"
    with pytest.raises(ValueError, match="0.2"):
        EventRecord.from_dict(lone).validate()


_DETERMINISM = r"""
import sys
sys.path.insert(0, {here!r})
import templateEmissionTest as t
from rmgpy.kmc import template_emission as te
payload = t.emit(t.homolysis(10, (3, 4)))
print(te.canonical_json(payload))
"""


def test_determinism_across_processes_hash_seeds_and_permuted_atoms():
    outs = []
    for seed in ("0", "12345"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        res = subprocess.run([sys.executable, "-c", _DETERMINISM.format(here=str(HERE))],
                             capture_output=True, text=True, env=env, timeout=300)
        assert res.returncode == 0, res.stderr[-2000:]
        outs.append(res.stdout.strip().splitlines()[-1])
    assert outs[0] == outs[1]
    w = homolysis(10, (3, 4))
    one = emit(w)
    perm = _perm(w, seed=3)
    two = emit(_relabel(w, perm), record=toy_record(w, 2))
    assert te.canonical_json(te._strip(one)) == te.canonical_json(te._strip(two))
    assert one["sha256"] == two["sha256"]
    m1, m2 = base_template(one)["maps"], base_template(two)["maps"]
    pw = _relabel(w, perm)
    assert _anchors(pw, m2["witness"]) == _anchors(pw, [perm[a] for a in m1["witness"]])


# --------------------------------------------------------------------------
# bimolecular counts (toy): matcher pair count = joint transition count
# --------------------------------------------------------------------------
@pytest.mark.parametrize("n1,n2,head1", [(5, 5, None), (3, 6, "methyl")])
def test_pair_count_equals_joint_transitions(n1, n2, head1):
    """Identical radicals (five C5 primary radicals: 10 unordered pairs) and
    unequal ones (a propyl radical whose methyl end is in its window, three
    copies, against two hexyl radicals: 6 ordered pairs)."""
    w = recombination(n1, n2, head1)
    c = cut(w, direction="forward")
    assert c.participants == 2
    t = c.matcher(1.0, "rec")
    assert t.role_swap is (head1 is None)
    copies = [alkane(n1, head=head1, tail="primary_radical")[0] for _ in range(3)] + \
        [alkane(n2, tail="primary_radical")[0] for _ in range(2)]
    rows, offset = {}, 0
    for r in copies:
        shift = {a: a + offset for a in r}
        rows.update({shift[a]: x[:5] + (tuple((shift[n], o) for n, o in x[5]),)
                     for a, x in r.items()})
        offset += 1000
    state = te.melt_state(rows)
    joint = wm.match([t], state)
    quotient = wm.match([t], state, joint=False)
    assert quotient.classes["rec"].pair_count == len(joint.transitions) > 0
    expected = 10 if head1 is None else 3 * 2
    assert len(joint.transitions) == expected


# --------------------------------------------------------------------------
# with the RMG database: forbidden moves excluded, the DP3 transfer produced,
# real bimolecular pair counts
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def rmg():
    import rmg_moves_harness as harness
    if not harness.DATABASE.is_dir():
        pytest.skip(f"needs the RMG database at {harness.DATABASE}")
    import channelsTest as ct
    from rmgpy.kmc import certificates as cm
    fam = ct.database().kinetics.families
    return fam, cm.Certifier(fam, max_channel_sets=16), te.RmgCounter(fam)


def _real_cut(rmg, family, specs, direction, pick, contexts=False):
    fam, certifier, _counter = rmg
    cases = win.member_cases(certifier, fam[family], win.builds_of(specs), direction)
    for case in cases:
        if not pick(case) or case.support != "supported" or \
                win.species_domain_refusal(case.witness):
            continue
        keep = []
        row = win.run_member(certifier, fam[family], case, bound=win.span_bound(
            win.family_constants(fam[family])), contexts=contexts, specs=specs, keep=keep)
        if row["status"] != "fixed_point":
            continue
        fp = keep[0]
        cut_ = te.cut_template(fp.witness, fp.window.atoms, family=family,
                               pathway=case.pathway, direction=direction,
                               labels=te.labels_of(case.moves[0]))
        return case, fp, cut_
    raise AssertionError(f"no supported fixed-point {family} member on {specs}")


def _copies(rows, n, offset=0):
    """``n`` copies of a component's rows, ids shifted by 10**7 per copy."""
    out = {}
    for c in range(n):
        shift = offset + (c + 1) * 10 ** 7
        out.update({a + shift: x[:5] + (tuple((m + shift, o) for m, o in x[5]),)
                    for a, x in rows.items()})
    return out


@pytest.fixture(scope="module")
def dp3(rmg):
    """The DP3 benzylic radical chain + DP3 pair: the synthetic
    unchanged_species transfer member and one kept member, each cut from its
    fixed-point window (whose padding the templates may need)."""
    S = b.PsChainSpec
    specs = [S(3, radical_units=(1,)), S(3)]
    transfer = _real_cut(rmg, "H_Abstraction", specs, "forward",
                         lambda x: x.kind == "synthetic:unchanged_species")
    # a kept member on a CH2 donor (the methine donor at the radical's own
    # position is the transfer again on the symmetric padded pair; the
    # benzylic CH2 members are refused upstream at a reacting port)
    kept = _real_cut(rmg, "H_Abstraction", specs, "forward",
                     lambda x: x.kind == "kept" and x.pathway[0] == "C/H2/NonDeC")
    return transfer, kept


def test_dp3_transfer_is_produced_and_required(rmg, dp3):
    """H abstraction between a DP3 benzylic radical chain and DP3: the
    synthetic unchanged_species transfer's template produces the transfer on
    its witness pair (required, matcher-only by policy); without it the
    required move is reported missing.  The window reaches the padding, so
    the template cannot apply to the bare DP3 proxies."""
    _fam, _certifier, counter = rmg
    (case, fp, c), _kept = dp3
    assert te.reaches_beyond(c, fp.witness, case.witness)
    t = c.matcher(1.0, "transfer")
    mols = te.molecules_of(te.template_sides(fp.witness)[0])
    table = te.classify_reference([t], mols, family="H_Abstraction", directions=("forward",),
                                  counter=counter,
                                  owners=lambda mv: ("refused:outside_this_test", "one template"))
    assert table["rmg_unchanged_species"] >= 1 and table["unchanged_produced"] >= 1
    assert not [f for f in table["failures"] if "forbidden" in f or "unchanged" in f], \
        table["failures"]
    empty = te.classify_reference([], mols, family="H_Abstraction", directions=("forward",),
                                  counter=counter,
                                  owners=lambda mv: ("emitted", "its record emitted"))
    assert any("required unchanged_species" in f for f in empty["failures"])


def test_forbidden_outcomes_are_never_produced(rmg, dp3):
    """The kept and transfer templates on the kept member's padded pair plus a
    second copy of its radical chain: RMG files the beta-H abstractions
    between the two radicals forbidden (disprop1_base_case); the templates
    match kept moves and never produce a forbidden-only outcome."""
    _fam, _certifier, counter = rmg
    (_c1, _fp1, transfer), (_c2, fp, kept) = dp3
    r = te.template_sides(fp.witness)[0]
    parts = rm.split_rows(dict(r))
    radical = next(x for x in parts if any(row[3] for row in x.values()))
    rows = dict(r)
    rows.update(_copies(radical, 1))
    mols = te.molecules_of(rows)
    assert len(mols) == 3
    templates = [kept.matcher(1.0, "kept"), transfer.matcher(1.0, "transfer")]
    table = te.classify_reference(templates, mols, family="H_Abstraction",
                                  directions=("forward",), counter=counter,
                                  owners=lambda mv: ("refused:outside_this_test", ""))
    assert table["rmg_forbidden"] >= 1 and table["matched"] >= 1
    assert not [f for f in table["failures"] if "forbidden" in f], table["failures"]


@pytest.fixture(scope="module")
def compile_db():
    """H abstraction and disproportionation loaded the compile's way (the
    production rate path checks the rate-rule preparation policy).  Not the
    process's cached load (``channelsTest.database``): that one is prepared
    without the verbose kinetics comments the compile keeps, and the
    compiler refuses it; the tests using this load run in their own process
    (``in_own_process``)."""
    import template_emission_acceptance as acc
    if not acc.DATABASE.is_dir():
        pytest.skip(f"needs the RMG database at {acc.DATABASE}")
    return acc.load_database(["H_Abstraction", "Disproportionation"])


# A second database beside the cached one does not fit the 6 GB test limit:
# a test needing both runs alone in a child process (the child sets the flag)
OWN_PROCESS = os.environ.get("RMG_KMC_OWN_PROCESS") == "1"
in_own_process = pytest.mark.skipif(
    not OWN_PROCESS, reason="runs in its own process (test_compile_way_database_test_alone)")


def _rates(compile_db, family):
    import template_emission_acceptance as acc
    return te.ProductionRates(compile_db, [family], database_path=acc.DATABASE,
                              database_sha="i138-test")


def _request(family, case, witness, orientation="as_generated"):
    r, p = te.template_sides(witness)
    return te.RateRequest(family, ";".join(case.pathway), orientation, r, p, witness,
                          "emitted", ())


@in_own_process
def test_identity_transfer_gets_an_executable_rate(rmg, dp3, compile_db):
    """The DP3 benzylic H transfer (central benzylic H, C/H/Cs2Cb;C_rad/CbCs2)
    whose products are its reactants: RMG's generator drops it, the explicit
    labelled reaction carries the family's estimate, and its emitted witness
    gets an ok row with a weight (as does the whole-species ladder member)."""
    fam, _certifier, counter = rmg
    (case, fp, c), _kept = dp3
    rates = _rates(compile_db, "H_Abstraction")
    req = _request("H_Abstraction", case, fp.witness)
    assert rates.identity(req)
    res = rates(req)
    assert not res.problems and res.provenance["path"] == "explicit_identity"
    assert all(k > 0 for k in res.k)
    row = te.witness_row(c.matcher(1.0, "t"), fp.witness, label="emitted", vector=(),
                         family=fam["H_Abstraction"], direction="forward",
                         pathway=te.count_pathway(fam["H_Abstraction"], case.pathway, False),
                         reverse=False, counter=counter, rates=rates,
                         record={"template": ";".join(case.pathway),
                                 "orientation": "as_generated"})
    assert row["status"] == "ok", row
    assert row["counts"]["matcher"] == row["counts"]["rmg"] == 1
    assert row["w"] == pytest.approx([k * res.ssa_multiplier for k in res.k], rel=1e-12)
    entry = te.ladder_member(case, {"status": "fixed_point"}, family=fam["H_Abstraction"],
                             counter=counter, rates=rates, orientation="as_generated")
    assert entry["status"] == "emitted", entry.get("refusal")
    assert all(w > 0 for w in entry["weight"]["w"])
    assert [x["ok"] for x in entry["calibration"]] == [True] * 3


@in_own_process
def test_explicit_rate_equals_production_rate_off_identity(rmg, compile_db):
    """The explicit labelled-reaction path against RMG's generator on
    non-identity transfers (CH2 and CH3 donors: degeneracy 2 and 3) and on a
    disproportionation of two identical radicals (the identical-reactant
    halving): equal k to 1e-9 at 600/700/800 K, equal ssa."""
    _fam, certifier, _counter = rmg
    S = b.PsChainSpec
    sets = [("H_Abstraction", [S(3, radical_units=(1,)), S(3)], 3),
            ("Disproportionation", [S(2, tail_end="benzylic_radical")] * 2, 1)]
    checked = 0
    for family, specs, n in sets:
        rates = _rates(compile_db, family)
        seen = []
        for case in win.member_cases(certifier, compile_db.kinetics.families[family],
                                     win.builds_of(specs), "forward"):
            if case.kind != "kept" or case.support != "supported" or case.pathway in seen:
                continue
            req = _request(family, case, case.witness)
            if rates.identity(req):
                continue
            seen.append(case.pathway)
            produced = rates(req)
            assert not produced.problems and produced.provenance["path"] == "generator"
            records, why = rates.explicit(req)
            hit = [x for x in records or () if x.orientation == "as_generated" and x.k_table
                   and x.template == req.template]
            assert len(hit) == 1, why
            k = dict(zip(hit[0].k_table["T"], hit[0].k_table["k"]))
            for t, value in zip(produced.temperatures, produced.k):
                if t in te.CALIBRATION_TEMPERATURES:
                    assert k[t] == pytest.approx(value, rel=1e-9)
            assert hit[0].ssa_multiplier == produced.ssa_multiplier
            checked += 1
            if len(seen) == n:
                break
    assert checked == 4


@pytest.mark.parametrize("name", ["test_identity_transfer_gets_an_executable_rate",
                                  "test_explicit_rate_equals_production_rate_off_identity"])
def test_compile_way_database_test_alone(name, tmp_path):
    """The compile-way database tests, each run alone in a child process
    (whatever other database this process holds), and required to pass."""
    if OWN_PROCESS:
        pytest.skip("the child runs the test itself")
    import rmg_moves_harness as harness
    if not harness.DATABASE.is_dir():
        pytest.skip(f"needs the RMG database at {harness.DATABASE}")
    out = subprocess.run(
        [sys.executable, "-m", "pytest", "-p", "no:cacheprovider", "--no-cov", "-q",
         f"--basetemp={tmp_path / 'child'}", f"{__file__}::{name}"],
        env=dict(os.environ, RMG_KMC_OWN_PROCESS="1"), capture_output=True, text=True,
        timeout=580)
    tail = out.stdout.strip().splitlines()[-1] if out.stdout.strip() else ""
    # the summary line names exactly one passed test and nothing else
    assert out.returncode == 0 and tail.strip("= ").split(" in ")[0] == "1 passed", \
        out.stdout[-3000:] + out.stderr[-2000:]


def test_c16_two_record_collision_fires_once_after_the_gate(rmg):
    """Records evt_1417a43aae8f and evt_f00a90e83259 of the C16 acceptance emit
    one C-C homolysis template (key f5338ad0...) with equal weights: loaded
    together nine cleavages of C16 fire twice and the reference comparison
    fails; after the catalog gate every cleavage fires once, at the single
    template's weight, and the comparison has no multiplicity failure."""
    fam, _certifier, counter = rmg
    payloads = json.loads((HERE / "fixtures" / "i138_c16_collision_payloads.json").read_text())
    mols, _parts = b.witness_reactants(win.builds_of([("alkane", 16)]))
    rows = {}
    for m in mols:
        rows.update(rm.molecule_rows(m))

    def fired(loaded):
        state = te.melt_state(rows)
        out = defaultdict(list)
        for tr in wm.match(loaded, state).transitions:
            out[rm.product_key(wm.product_rows(state, tr))].append(tr)
        return out

    before = fired(te.loaded_templates(payloads, 700.0))
    assert sum(1 for v in before.values() if len(v) == 2) == 9
    table = te.classify_reference(te.loaded_templates(payloads, 700.0), list(mols),
                                  family="R_Recombination", directions=("reverse",),
                                  counter=counter, owners=lambda mv: ("refused:other", ""))
    assert any("produced 2 times" in f for f in table["failures"])
    merged, collisions = te.resolve_collisions(payloads)
    assert [(r["copies"], r["disposition"]) for r in collisions] == [(2, "merged")]
    assert collisions[0]["key"].startswith("f5338ad0")
    after = fired(te.loaded_templates(merged, 700.0))
    assert set(after) == set(before) and {len(v) for v in after.values()} == {1}
    w700 = 5.513978431788572e-10
    assert all(v[0].weight == pytest.approx(w700, rel=1e-12) for v in after.values())
    table = te.classify_reference(te.loaded_templates(merged, 700.0), list(mols),
                                  family="R_Recombination", directions=("reverse",),
                                  counter=counter, owners=lambda mv: ("refused:other", ""))
    assert not [f for f in table["failures"] if "produced" in f]


@pytest.mark.parametrize("second", ["identical", "unequal"])
def test_real_recombination_pair_count_equals_joint_transitions(rmg, second):
    """A real recombination template, cut from its padded fixed-point window,
    on copies of its witness's reactant species: the quotient pair count
    equals the joint transition count (3 identical copies: 3 unordered pairs;
    2 + 2 unequal: 4 ordered pairs), the role swap set by the matcher."""
    S = b.PsChainSpec
    first = S(1, tail_end="benzylic_radical")
    other = first if second == "identical" else S(2, head_end="primary_radical")
    case, fp, c = _real_cut(rmg, "R_Recombination", [first, other], "forward", lambda x: True)
    t = c.matcher(1.0, "rec")
    assert t.role_swap is (second == "identical")
    parts = rm.split_rows(dict(te.template_sides(fp.witness)[0]))
    assert len(parts) == 2
    if second == "identical":
        assert rmg[2].species_of_rows(parts[0]) == rmg[2].species_of_rows(parts[1])
        rows = _copies(parts[0], 3)
    else:
        rows = _copies(parts[0], 2)
        rows.update(_copies(parts[1], 2, offset=5 * 10 ** 6))
    state = te.melt_state(rows)
    joint = wm.match([t], state)
    quotient = wm.match([t], state, joint=False)
    assert quotient.classes["rec"].pair_count == len(joint.transitions) > 0
    assert len(joint.transitions) == (3 if second == "identical" else 4)
    # the reference comparison runs RMG on the pair (R_Recombination's
    # template is one tree root: its forward direction is still bimolecular)
    assert te._molecularity(rmg[0]["R_Recombination"], "forward") == 2
    pair = te.molecules_of(dict(te.template_sides(fp.witness)[0]))
    table = te.classify_reference([t], pair, family="R_Recombination", directions=("forward",),
                                  counter=rmg[2], owners=lambda mv: ("refused:other_site", ""))
    assert table["rmg_kept"] >= 1 and table["matched"] == 1 and not table["failures"], table


# --------------------------------------------------------------------------
# slow: the acceptance on the reference compiles
# --------------------------------------------------------------------------
@pytest.mark.slow
@pytest.mark.skipif(not SLOW, reason="set RMG_KMC_SLOW=1")
def test_reference_acceptance():
    out = Path(os.environ.get("RMG_KMC_I138_OUT", "/tmp/i138_acceptance"))
    res = subprocess.run([sys.executable, str(HERE / "template_emission_acceptance.py"),
                          str(out), "c16", "--sample", "1"], capture_output=True, text=True)
    assert res.returncode == 0, res.stdout[-3000:] + res.stderr[-3000:]


# --------------------------------------------------------------------------
# I-147: weighted collision merge, reciprocal pairs, single orbit, equilibrium
# cycle, class overlap, subset loading, the identity-transfer guard
# --------------------------------------------------------------------------
import itertools


def _rid(k):
    return "evt_" + format(k, "x") * 64


def _variant(payload, rid, shift=0.0, *, key=None, reverse_of=None, grid=None):
    """A copy of ``payload`` as record ``rid``: every emitted weight moved by
    ``shift`` decades, optionally another key, reverse link or T grid."""
    p = json.loads(json.dumps(payload))
    p["record_event_id"] = rid
    if reverse_of is not None:
        p["reverse_of"] = reverse_of
    for t in p["templates"]:
        if t["status"] != "emitted":
            continue
        t["weight"]["w"] = [w * 10 ** shift for w in t["weight"]["w"]]
        if grid is not None:
            t["weight"]["T"] = list(grid)
            t["weight"]["w"] = t["weight"]["w"][:len(grid)]
        if key is not None:
            t["key"] = key
    return te.finish(p)


def _group(shifts, **kw):
    base = _toy_catalog()[0]
    return [_variant(base, _rid(k + 1), s, **kw) for k, s in enumerate(shifts)]


def _canon(payloads):
    return te.canonical_json(sorted(payloads, key=lambda p: p["record_event_id"]))


def test_collision_diameter_is_not_pairwise_or_chained():
    """Three copies 0.006 decades apart in a row: every consecutive pair is
    within 0.01 but the diameter is 0.012: refused, in any order."""
    members = _group([0.0, 0.006, 0.012])
    for order in itertools.permutations(members):
        out, table = te.resolve_collisions(list(order))
        assert [r["disposition"] for r in table] == ["refused:collision"]
        assert table[0]["log10_diameter"] == pytest.approx(0.012, abs=1e-9)
        for p in out:
            t = base_template(p)
            assert t["status"] == "refused" and t["refusal"]["status"] == "refused:collision"


def test_weighted_collision_merges_within_the_diameter_in_any_order():
    """AMENDED (I-162 S-g rework 1: the representative no longer follows the
    record ids): the same three copies 0.004 apart merge into the template of
    the copy of least ``partition_free_identity``, at its weight; every input
    order gives one merged catalog."""
    members = _group([0.0, 0.004, 0.008])
    owner = min(members, key=lambda p: te.partition_free_identity(p, base_template(p)))
    rid = owner["record_event_id"]
    others = sorted({_rid(1), _rid(2), _rid(3)} - {rid})
    canon, owner_w = set(), None
    for order in itertools.permutations(members):
        out, table = te.resolve_collisions(list(order))
        canon.add(_canon(out) + te.canonical_json(table))
        assert [r["disposition"] for r in table] == ["merged"] and table[0]["copies"] == 3
        assert table[0]["log10_diameter"] == pytest.approx(0.008, abs=1e-9)
        by = {p["record_event_id"]: base_template(p) for p in out}
        assert by[rid]["status"] == "emitted"
        assert sorted(by[rid]["collision"]["records"]) == [_rid(1), _rid(2), _rid(3)]
        assert by[rid]["collision"]["records"][0] == rid
        assert by[rid]["weight"] == base_template(owner)["weight"]
        assert by[others[0]]["status"] == by[others[1]]["status"] == "merged"
        assert by[others[0]]["merged_into"]["record"] == rid
        for p in out:
            te.validate_payload(p)
    assert len(canon) == 1


def test_collision_on_different_temperature_grids_is_a_named_refusal():
    """Copies on a 29-point and a 3-point grid cannot be compared: refused by
    name, no exception."""
    base = _toy_catalog()[0]
    a = _variant(base, _rid(1))
    b = _variant(base, _rid(2), grid=[300.0, 650.0, 1000.0])
    out, table = te.resolve_collisions([a, b])
    assert [r["disposition"] for r in table] == ["refused:collision_grid"]
    for p in out:
        t = base_template(p)
        assert t["status"] == "refused" and t["refusal"]["status"] == "refused:collision_grid"
        assert "grid" in t["refusal"]["reason"]
        te.validate_payload(p)


def test_catalog_refuses_a_zero_or_unusable_weight_or_channel():
    """T_ch = 0, a non-finite and a non-positive weight at any grid
    temperature fail the catalog gate (compile errors)."""
    base = _toy_catalog()[0]
    for bad in (0.0, -1.0, float("inf"), float("nan")):
        p = json.loads(json.dumps(base))
        base_template(p)["weight"]["w"][7] = bad
        with pytest.raises(te.EmissionError, match="weight"):
            te.resolve_collisions([p])
    p = json.loads(json.dumps(base))
    base_template(p)["counts"]["context"]["matcher"] = 0
    with pytest.raises(te.EmissionError, match="T_ch"):
        te.resolve_collisions([p])


def _reciprocal(forward_keys, reverse_keys, reverse_shift=(0.0, 0.0, 0.0)):
    """Forward records 1.. and their reverse records 5..: ``*_keys`` give
    each record's base template key."""
    base = _toy_catalog()[0]
    out = []
    for k, key in enumerate(forward_keys):
        out.append(_variant(base, _rid(k + 1), key=key, reverse_of=_rid(k + 5)))
    for k, key in enumerate(reverse_keys):
        out.append(_variant(base, _rid(k + 5), reverse_shift[k], key=key,
                            reverse_of=_rid(k + 1)))
    return out


def test_a_merge_that_splits_a_reciprocal_pair_fails():
    """Two forward records merge but their reverses carry different keys."""
    kf, k1, k2 = "a" * 64, "b" * 64, "c" * 64
    with pytest.raises(te.ReciprocalSplitError, match="reciprocal"):
        te.resolve_collisions(_reciprocal([kf, kf], [k1, k2]))
    # reverses merged, forward not
    with pytest.raises(te.ReciprocalSplitError, match="reciprocal"):
        te.resolve_collisions(_reciprocal([k1, k2], [kf, kf]))
    # a reverse copy of a third forward record lands in the merged reverse
    # template while its forward partner stays alone
    payloads = _reciprocal([kf, kf, k1], [k2, k2, k2])
    with pytest.raises(te.ReciprocalSplitError, match="reciprocal"):
        te.resolve_collisions(payloads)


def test_a_reciprocal_pair_merges_together_with_a_bijective_link():
    kf, kr = "a" * 64, "b" * 64
    out, table = te.resolve_collisions(_reciprocal([kf, kf], [kr, kr], (0.003, 0.0)))
    assert sorted(r["disposition"] for r in table) == ["merged", "merged"]
    links = {}
    for p in out:
        for t in p["templates"]:
            if t["status"] == "emitted":
                links[t["key"]] = t["reverse"]["key"]
                assert t["reverse"]["records"] == sorted(
                    {_rid(1), _rid(2)} if t["key"] == kr else {_rid(5), _rid(6)})
    assert links == {kf: kr, kr: kf}


def test_orbit_problem_on_a_witness():
    """Transitions at the mirror bonds of a chain are one orbit; two inequivalent
    bonds are not."""
    w = homolysis(12, (5, 6))
    rows = dict(w.reactant_rows())
    t = cut(w).matcher(1.0, "t")
    state = te.melt_state(rows)
    edits = [x.edits for x in wm.match([t], state).transitions]
    assert len(edits) >= 4
    assert te.orbit_problem(rows, [edits[0]]) is None
    assert te.orbit_problem(rows, [edits[0], edits[0]]) is None
    why = te.orbit_problem(rows, [edits[0], edits[1]])
    assert why and "orbit" in why
    # the chain's reflection maps the first cleavage onto the last
    assert te.orbit_problem(rows, [edits[0], edits[-1]]) is None


class ForcedCounter(ToyCounter):
    """The channel is every transition of ``template`` (the species check of
    the toy oracle dropped): a class whose transitions are not symmetric."""

    def __init__(self, template):
        super().__init__()
        self.template = template

    def channel(self, family, molecules, direction, pathway, product_species):
        rows = {}
        for m in molecules:
            rows.update(rm.molecule_rows(m))
        state = te.melt_state(rows)
        keys = sorted({rm.product_key(wm.product_rows(state, x))
                       for x in wm.match([self.template], state).transitions})
        return te.RmgChannel(tuple(keys), 2 * len(keys), tuple(((0,), (0,)) for _ in keys),
                             (("kept", 2 * len(keys)),))


def test_a_class_that_is_not_one_orbit_is_refused_non_transferable():
    w = homolysis(12, (5, 6))
    template = cut(w).matcher(1.0, "t")
    counter = ForcedCounter(template)
    case, row, ctx = member(w)
    record = toy_record(w, 1)
    out = te.emit_record(record, family=ToyFamily, members=[(case, row, ctx)],
                         counter=counter, rates=stub_rates(counter), check_upstream=False)
    t = base_template(out)
    assert t["status"] == "refused" and t["refusal"]["status"] == "refused:non_transferable"
    assert out["status"] == "refused"
    # the ordinary mirror class of the same bond stays emitted
    assert base_template(emit(w))["status"] == "emitted"


# -- reciprocal pair on a witness ------------------------------------------
def _td(cp, h298, s298, cp0, cpinf):
    return {"Tdata": [300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0], "Cpdata": cp,
            "H298": h298, "S298": s298, "Cp0": cp0, "CpInf": cpinf}


# a toy thermo with dCp != 0, so that ln Kc is NOT linear in 1/T: A <-> B + C
TOY_A = _td([200.0, 220.0, 240.0, 260.0, 290.0, 310.0, 340.0], -100000.0, 500.0, 100.0, 700.0)
TOY_B = _td([90.0, 100.0, 108.0, 115.0, 130.0, 140.0, 155.0], -30000.0, 250.0, 60.0, 400.0)
TOY_C = _td([125.0, 135.0, 147.0, 160.0, 175.0, 185.0, 200.0], -50000.0, 330.0, 60.0, 400.0)


def toy_form(homolysis_direction):
    return {"kind": "ThermoData-Kc/1", "P0": 1e5,
            "reactants": [TOY_A] if homolysis_direction else [TOY_B, TOY_C],
            "products": [TOY_B, TOY_C] if homolysis_direction else [TOY_A]}


class ToyKc:
    """The toy thermo through the same exact evaluation the payload stores."""

    def form(self, reactant_rows, product_rows, temperatures):
        return toy_form(len(rm.split_rows(dict(reactant_rows))) == 1)

    def __call__(self, reactant_rows, product_rows, temperatures):
        f = self.form(reactant_rows, product_rows, temperatures)
        return [te.kc_at(f, t) for t in temperatures]


toy_kc = ToyKc()
def _kc_grid():
    return tuple(te.kc_at(toy_form(False), t) for t in GRID)


KC_CACHE = []


def _kc():
    if not KC_CACHE:
        KC_CACHE.append(_kc_grid())
    return KC_CACHE[0]


def _pair_rates(by_orientation):
    def rates(req):
        return te.RateResult(GRID, by_orientation[req.orientation], 1.0, {}, ())
    return rates


def _toy_pair(hid=None, wid=None, scale=1.0, h_scale=1.0, counter=None, certificate=None):
    """Homolysis H of C12 at (3,4) (reversed record) and the recombination W
    of its two products (as_generated), linked by ``reverse_of``;
    ``k_W = Kc * k_H`` times ``scale`` both ways, ``h_scale`` on H alone."""
    h = homolysis(12, (3, 4))
    r = recombination(4, 8)
    k_h = tuple(k * scale * h_scale for k in K_RULE)
    k_w = tuple(k * scale * c for k, c in zip(K_RULE, _kc()))
    hid, wid = hid or _rid(1), wid or _rid(2)
    rec_h = toy_record(h, 1)
    rec_h.update(event_id=hid, reverse_of=wid)
    rec_h["k_table"]["k"] = list(k_h)
    rec_w = toy_record(r, 1, orientation="as_generated")
    rec_w.update(event_id=wid, reverse_of=hid, reactant_multiplicities=[1, 1])
    rec_w["k_table"]["k"] = list(k_w)
    rates = _pair_rates({"reversed": k_h, "as_generated": k_w})
    counter = counter or ToyCounter()
    case_h, row_h, ctx_h = member(h)
    p_h = te.emit_record(rec_h, family=ToyFamily, members=[(case_h, row_h, ctx_h)],
                         counter=counter, rates=rates, check_upstream=False, kc=toy_kc,
                         record_certificate=certificate)
    case_w, row_w, ctx_w = member(r, direction="forward")
    p_w = te.emit_record(rec_w, family=ToyFamily, members=[(case_w, row_w, ctx_w)],
                         counter=counter, rates=rates, check_upstream=False, kc=toy_kc,
                         record_certificate=certificate)
    return [p_h, p_w], counter


def _check(payloads, counter, **kw):
    return te.check_cycles(payloads, counter=counter, kc=toy_kc, **kw)


def _res(payloads, counter):
    """The catalog gate: collisions resolved, reciprocal weights derived."""
    resolved = te.resolve_collisions([te.finish(p) for p in payloads])[0]
    return te.derive_weights(resolved, counter=counter)


def test_toy_reciprocal_pair_satisfies_the_cycle_against_thermo():
    """AMENDED (it encoded the circular ratio): the pair passes against the
    independent thermo Kc, which the table reports."""
    payloads, counter = _toy_pair()
    assert [p["status"] for p in payloads] == ["emitted", "emitted"]
    out, table = _check(_res(payloads, counter), counter)
    assert len(table) == 1 and table[0]["status"] == "ok", table
    assert table[0]["max_ratio_error"] <= 1e-12 and table[0]["rows"] >= 4
    assert table[0]["max_offgrid_error"] <= 1e-12 and table[0]["temperatures"] > len(GRID)
    T = table[0]["worst"]["T"]
    assert table[0]["worst"]["Kc"] in (pytest.approx(te.kc_at(toy_form(True), T), rel=1e-12),
                                       pytest.approx(te.kc_at(toy_form(False), T), rel=1e-12))
    assert [p["status"] for p in out] == ["emitted", "emitted"]
    assert te.one_way_problems(out) == [] and te.cycle_gates(table, out) == []


def test_homolysis_rate_and_weight_doubled_is_detected_against_thermo():
    """k_H and its weights doubled together keep the records' own k_F / k_R
    self-consistent (and the calibration exact); the independent Kc is
    unchanged, so the loaded pair fails."""
    payloads, counter = _toy_pair(h_scale=2.0)
    assert all(c["ok"] for c in payloads[0]["calibration"])
    out, table = _check(_res(payloads, counter), counter)
    assert table[0]["status"] == "refused:cycle"
    assert "reverse-by-Kc" in table[0]["reason"] and te.cycle_gates(table, out)


def test_perturbing_one_weight_by_1e_6_fails_the_cycle_and_is_gated():
    """AMENDED: perturb the loaded recombination weight at one temperature.
    Both records are refused, the failure is an acceptance gate, and nothing
    of the pair stays loaded (no one-way channel)."""
    payloads, counter = _toy_pair()
    base_template(payloads[1])["weight"]["w"][10] *= 1 + 1e-6
    out, table = _check(_res(payloads, counter), counter)
    assert table[0]["status"] == "refused:cycle"
    assert "reverse-by-Kc" in table[0]["reason"]
    for p in out:
        t = base_template(p)
        assert t["status"] == "refused" and t["refusal"]["status"] == "refused:cycle"
        assert p["status"] == "refused"
    gates = te.cycle_gates(table, out)
    assert gates and "refused:cycle" in gates[0]
    assert te.one_way_problems(out) == [] and te.loaded_templates(out, 700.0) == []


def test_a_cycle_whose_inverse_edit_is_not_unique_is_refused():
    """AMENDED (signature): the reverse record carries a second base template
    making the same cleavage: each forward transition is undone twice."""
    payloads, counter = _toy_pair()
    twin = json.loads(json.dumps(base_template(payloads[1])))
    twin["key"] = "d" * 64
    payloads[1]["templates"].append(twin)
    out, table = _check(_res(payloads, counter), counter)
    assert table[0]["status"] == "refused:cycle" and table[0]["bijection_failures"] >= 1
    assert "inverse edit" in table[0]["reason"]


def test_a_cycle_without_a_channel_transition_is_refused():
    """AMENDED (signature): the forward record's template cannot make the
    reverse record's cleavage on its witness: named refusal of both."""
    payloads, counter = _toy_pair()
    other = emit(homolysis(12, (5, 6)), radius=4)
    payloads[0] = te.finish(dict(other, record_event_id=payloads[0]["record_event_id"],
                                 reverse_of=payloads[0]["reverse_of"]))
    out, table = _check(_res(payloads, counter), counter)
    assert table[0]["status"] == "refused:cycle"
    for p in out:
        assert base_template(p)["refusal"]["status"] == "refused:cycle"


def test_a_pair_with_an_emitted_witness_weight_is_checked():
    """Templates reaching beyond the proxy graph weigh on their emitted
    witness (``counts.record`` absent): checked there, not skipped."""
    payloads, counter = _toy_pair()
    for p in payloads:
        p["weight_source"], p["counts"] = "emitted_witness", {"record": None}
    out, table = _check(_res(payloads, counter), counter)
    assert table[0]["status"] == "ok" and table[0]["rows"] >= 4
    base_template(payloads[0])["weight"]["w"][3] *= 1.001
    out, table = _check(_res(payloads, counter), counter)
    assert table[0]["status"] == "refused:cycle"


def test_a_pair_of_context_templates_is_checked():
    """No base template on either side (only context templates load): the pair
    is checked on its rows (grid, midpoints, random temperatures), its second
    direction derived like any other, and a bad context weight is caught."""
    payloads, counter = _toy_pair()
    for p in payloads:
        p["templates"][0]["base"] = False
    out, table = _check(_res(payloads, counter), counter)
    assert table[0]["status"] == "ok" and table[0]["rows"] >= 4
    assert table[0]["max_offgrid_error"] <= 1e-12
    payloads[0]["templates"][0]["weight"]["w"][5] *= 1.01
    out, table = _check(_res(payloads, counter), counter)
    assert table[0]["status"] == "refused:cycle" and te.cycle_gates(table, out)


def test_a_live_pair_without_its_reverse_record_is_refused_when_complete():
    payloads, counter = _toy_pair()
    out, table = _check(payloads[:1], counter)
    assert table[0]["status"] == "refused:cycle" and "not in the catalog" in table[0]["reason"]
    assert te.loaded_templates(out, 700.0) == []
    out, table = _check(payloads[:1], counter, complete=False)
    assert table[0]["status"] == "not_evaluable" and base_template(out[0])["status"] == "emitted"
    assert te.cycle_gates(table, out, complete=False) == []
    assert len(te.cycle_not_evaluable(table)) == 1


def _two_pairs(shift_b=0.0, b_first=False):
    """Pair A (records 1/6) and pair B (2/5), B's rates scaled by
    10**shift_b both ways (Kc unchanged); the least forward id is A's, the
    least reverse id is B's: independent least-id representatives come from
    different pairs."""
    a, counter = _toy_pair(_rid(1), _rid(6), counter=None)
    b, _ = _toy_pair(_rid(2), _rid(5), scale=10 ** shift_b, counter=counter)
    return a + b, counter


def test_unit_representatives_come_from_one_source_pair_and_the_loaded_catalog_is_exact():
    """The C16 repro: scaling a pair both ways by 10**0.004 (Kc, calibration
    preserved) must not leave a loaded forward weight from one pair and a
    loaded reverse weight from the other."""
    payloads, counter = _two_pairs(0.004)
    out, table = te.resolve_collisions(payloads)
    out = te.derive_weights(out, counter=counter)
    assert sorted(r["disposition"] for r in table) == ["merged", "merged"]
    live = {p["record_event_id"]: base_template(p)["status"] for p in out}
    # AMENDED (I-162 S-g rework 2): the source pair is chosen by content, not
    # by the least ids; whichever it is, both representatives are one pair's
    assert live in ({_rid(1): "emitted", _rid(6): "emitted", _rid(2): "merged", _rid(5): "merged"},
                    {_rid(2): "emitted", _rid(5): "emitted", _rid(1): "merged", _rid(6): "merged"})
    assert all(len(r["pair_deviations"]) == 2 for r in table)
    out2, cycles = _check(out, counter)
    assert all(c["status"] == "ok" for c in cycles), cycles
    assert all(c["max_offgrid_error"] <= 1e-12 for c in cycles
               if c["worst"]["status"] == "emitted")
    assert max(c["max_ratio_error"] for c in cycles) <= 1e-12 or \
        all(c["worst"]["status"] == "merged" for c in cycles if c["max_ratio_error"] > 1e-12)
    assert te.one_way_problems(out2) == []


def test_reciprocal_groups_are_admitted_or_refused_together():
    """Forward copies agree, reverse copies are 10 % apart: the forward group
    is not merged alone (and no half-merged unit reaches the catalog)."""
    payloads, counter = _two_pairs(0.0)
    base_template(payloads[3])["weight"]["w"] = [w * 1.1 for w in base_template(payloads[3])["weight"]["w"]]
    payloads[3] = te.finish(payloads[3])
    out, table = te.resolve_collisions(payloads)
    assert {r["disposition"] for r in table} == {"refused:collision"}
    assert len(table) == 2
    assert [base_template(p)["status"] for p in out] == ["refused"] * 4
    assert te.loaded_templates(out, 700.0) == []


def test_merge_beyond_the_pair_tolerance_is_refused_with_every_deviation_reported():
    """Each group's diameter is 0.008 (inside), but the second pair's
    loaded ratio is 0.016 decades off: refused, deviations in the table."""
    payloads, counter = _two_pairs(0.0)
    for p in payloads:        # legacy payloads: the expectation is each pair's own ratio
        for t in p["templates"]:
            t.pop("thermo")
    base_template(payloads[2])["weight"]["w"] = [w * 10 ** 0.008 for w in base_template(payloads[2])["weight"]["w"]]
    base_template(payloads[3])["weight"]["w"] = [w * 10 ** -0.008 for w in base_template(payloads[3])["weight"]["w"]]
    payloads[2], payloads[3] = te.finish(payloads[2]), te.finish(payloads[3])
    out, table = te.resolve_collisions(payloads)
    assert {r["disposition"] for r in table} == {"refused:collision"}
    assert table[0]["max_pair_deviation_log10"] == pytest.approx(0.016, abs=1e-9)
    assert sorted(d["deviation_log10"] for d in table[0]["pair_deviations"]) == \
        pytest.approx([0.0, 0.016], abs=1e-9)


def test_a_cycle_failure_of_the_representative_refuses_the_merged_partner_too():
    """A 1e-6 error in the loaded representative of a merged reciprocal unit
    refuses the representatives AND their merged copies on both sides, is
    gated, and leaves no loaded reverse link to a non-loaded template."""
    payloads, counter = _two_pairs(0.0)
    out, _t = te.resolve_collisions(payloads)
    out = te.derive_weights(out, counter=counter)
    rep = next(p for p in out if p["record_event_id"] == _rid(1))
    base_template(rep)["weight"]["w"][10] *= 1 + 1e-6
    out = [te.finish(p) for p in out]
    out2, table = _check(out, counter)
    assert {p["record_event_id"]: base_template(p)["status"] for p in out2} == \
        {_rid(k): "refused" for k in (1, 2, 5, 6)}
    assert te.loaded_templates(out2, 700.0) == [] and te.one_way_problems(out2) == []
    gates = te.cycle_gates(table, out2)
    assert gates and all("refused:cycle" in g for g in gates)


def test_a_loaded_template_with_a_non_loaded_reverse_is_a_one_way_gate():
    payloads, counter = _toy_pair()
    base_template(payloads[1])["status"] = "refused"
    base_template(payloads[1])["refusal"] = {"status": "refused:window", "reason": "x"}
    payloads[1] = te.finish(payloads[1])
    assert [r for r, _w in te.one_way_problems(payloads)] == [payloads[0]["record_event_id"]]
    out, table = _check(payloads, counter)
    assert te.loaded_templates(out, 700.0) == []
    assert te.cycle_gates(table, out)


def test_thermo_kc_is_whole_species_rmg_thermo_in_concentration_units(rmg):
    """ThermoKc on the real thermo path: ethane <-> 2 CH3 against the
    reaction's own Kc, and the A + A convention ssa ratio 2."""
    import channelsTest as ct
    from rmgpy.kmc.reference_thermo import GasPhaseRMGReferenceThermo
    from rmgpy.molecule import Molecule
    from rmgpy.reaction import Reaction
    from rmgpy.species import Species
    db = ct.database()
    prov = GasPhaseRMGReferenceThermo(db.thermo, "test-only-database-commit")
    kc = te.ThermoKc(prov)
    def rows_of(smiles, offset):
        m = Molecule().from_smiles(smiles)
        m.assign_atom_ids()
        return {a + offset: (x[0], x[1], x[2], x[3], x[4], tuple((n + offset, o) for n, o in x[5]))
                for a, x in rm.molecule_rows(m).items()}
    ethane = rows_of("CC", 0)
    two = dict(rows_of("[CH3]", 100))
    two.update(rows_of("[CH3]", 200))
    got = kc(ethane, two, [700.0, 800.0])
    sp = lambda s: Species(molecule=[Molecule().from_smiles(s)])
    rxn = Reaction(reactants=[sp("CC")], products=[sp("[CH3]"), sp("[CH3]")])
    rxn_kc = list(prov.evaluate(rxn, [700.0, 800.0]).equilibrium_constants)
    assert got == pytest.approx(rxn_kc, rel=1e-9)
    assert te.ssa_of_rows(two, ToyCounterReal(rmg)) == 2.0
    assert te.ssa_of_rows(ethane, ToyCounterReal(rmg)) == 1.0


def ToyCounterReal(rmg):
    return rmg[2]


# --------------------------------------------------------------------------
# I-147 rework 2: detailed balance at every runtime temperature
# --------------------------------------------------------------------------
OFF_GRID = (312.5, 437.3, 612.5, 888.8)


def test_interpolation_is_linear_in_inverse_temperature():
    """A table of w = exp(a - b / T) is reproduced exactly between its grid
    points (linear-in-T interpolation of ln w is not)."""
    ts = [300.0, 400.0, 600.0]
    ws = [math.exp(30.0 - 20000.0 / t) for t in ts]
    for t in (333.3, 350.0, 512.5):
        assert te._interpolate(ts, ws, t) == pytest.approx(math.exp(30.0 - 20000.0 / t), rel=1e-12)
    with pytest.raises(te.EmissionError, match="extrapolation"):
        te._interpolate(ts, ws, 299.0)


def test_interpolation_error_of_a_stored_table():
    ts = [float(t) for t in range(300, 1001, 25)]
    smooth = [math.exp(30.0 - 20000.0 / t) for t in ts]
    err, _ = te.interpolation_error(ts, smooth)
    assert err <= 1e-12
    kinked = list(smooth)
    kinked[10] *= 1.05
    err, at = te.interpolation_error(ts, kinked)
    assert err == pytest.approx(0.05, rel=0.05) and at == ts[10]


def test_stored_thermo_form_reproduces_kc_at_any_temperature():
    """``kc_at`` of the toy form against the closed expression of
    ``Reaction.get_equilibrium_constant`` (Ka = exp(-dG/RT), C0 = P0/RT)."""
    from rmgpy import constants
    from rmgpy.thermo import ThermoData
    td = lambda d: ThermoData(Tdata=(d["Tdata"], "K"), Cpdata=(d["Cpdata"], "J/(mol*K)"),
                              H298=(d["H298"], "J/mol"), S298=(d["S298"], "J/(mol*K)"),
                              Cp0=(d["Cp0"], "J/(mol*K)"), CpInf=(d["CpInf"], "J/(mol*K)"))
    for T in (301.7, 512.5, 999.0):
        dg = td(TOY_B).get_free_energy(T) + td(TOY_C).get_free_energy(T) - td(TOY_A).get_free_energy(T)
        want = math.exp(-dg / constants.R / T) * (1e5 / constants.R / T)
        assert te.kc_at(toy_form(True), T) == pytest.approx(want, rel=1e-13)
        assert te.kc_at(toy_form(True), T) * te.kc_at(toy_form(False), T) == pytest.approx(1.0, rel=1e-12)


def _primary_derived(out):
    """(primary index, derived index) of a resolved toy pair."""
    (d,) = [k for k, p in enumerate(out) if base_template(p).get("derived")]
    return 1 - d, d


def test_the_second_direction_of_a_unit_is_derived_from_the_primary_at_the_requested_T():
    """AMENDED (I-162 S-g rework 1): primary = the emitted template of the
    as_generated record (here the recombination), then of the lesser template
    key; never by record id (the ids carry the artifact's family list).  The
    other direction's weight at an off-grid T is w_P(T) * (N_P / N_D) /
    (ssa_ratio * Kc_P(T)) with Kc evaluated exactly at T, not interpolated,
    not taken at the nearest grid T."""
    payloads, counter = _toy_pair()
    out = _res(payloads, counter)
    ip, idd = _primary_derived(out)
    h, w = base_template(out[ip]), base_template(out[idd])
    assert out[ip]["orientation"] == "as_generated" and out[idd]["orientation"] == "reversed"
    assert "derived" not in h and w["derived"]["primary"] == {"record": out[ip]["record_event_id"],
                                                              "key": h["key"]}
    assert "as_generated record, then of the lesser template key" in w["derived"]["rule"]
    form = toy_form(out[ip]["orientation"] == "reversed")        # the homolysis is unimolecular
    for T in OFF_GRID:
        wh = te.weight_at(h, T)
        got = te.weight_at(w, T, h)
        want = wh * w["derived"]["counts_ratio"] / (w["derived"]["ssa_ratio"] * te.kc_at(form, T))
        assert got == pytest.approx(want, rel=1e-12)
        # and the loaded sums satisfy S_P / S_D = ssa_ratio * Kc_P exactly
        nh = te._count_of(out[ip], h)
        nw = te._count_of(out[idd], w)
        assert (nh * wh) / (nw * got) == pytest.approx(
            h["thermo"]["ssa_ratio"] * te.kc_at(form, T), rel=1e-12)
    with pytest.raises(te.EmissionError, match="primary"):
        te.weight_at(w, 500.0)
    with pytest.raises(te.EmissionError, match="extrapolation"):
        te.weight_at(w, 1200.0, h)
    # the catalog loads it through the primary
    assert len(te.loaded_templates(out, 437.3)) == 2


def test_toy_pair_balances_at_off_grid_temperatures():
    payloads, counter = _toy_pair()
    out, table = _check(_res(payloads, counter), counter, extra_temperatures=OFF_GRID)
    assert table[0]["status"] == "ok"
    assert set(OFF_GRID) <= set(table[0]["extra"])
    assert max(table[0]["extra"].values()) <= 1e-12
    assert table[0]["max_offgrid_error"] <= 1e-12 and table[0]["temperatures"] >= 2 * len(GRID) - 1


def test_independently_interpolated_directions_are_gated_off_the_grid():
    """Both weight tables are exact at every grid temperature but the second
    direction is not derived (independent interpolation): the midpoint and
    random temperatures expose it."""
    payloads, counter = _toy_pair()
    out = _res(payloads, counter)
    base_template(out[_primary_derived(out)[1]]).pop("derived")
    out = [te.finish(p) for p in out]
    out2, table = _check(out, counter)
    assert table[0]["status"] == "refused:cycle" and table[0]["max_grid_error"] <= 1e-12
    assert table[0]["max_offgrid_error"] > 1e-10 and "off the grid" in table[0]["reason"]
    gates = te.cycle_gates(table, out2)
    assert gates and te.loaded_templates(out2, 500.0) == []


def test_a_merged_only_failure_refuses_the_owner_and_every_copy():
    """A merged copy's record that cannot be evaluated (no stored witness)
    refuses its representative and the other copies too."""
    payloads, counter = _two_pairs(0.0)
    out, _t = te.resolve_collisions(payloads)
    out = te.derive_weights(out, counter=counter)
    copy = next(p for p in out if p["record_event_id"] == _rid(2))
    assert base_template(copy)["status"] == "merged"
    base_template(copy)["maps"].pop("witnesses")
    out = [te.finish(p) for p in out]
    out2, table = _check(out, counter)
    assert {p["record_event_id"]: base_template(p)["status"] for p in out2} == \
        {_rid(k): "refused" for k in (1, 2, 5, 6)}
    assert te.loaded_templates(out2, 700.0) == [] and te.one_way_problems(out2) == []
    assert te.cycle_gates(table, out2)


def test_an_absent_reverse_on_a_partial_catalog_is_reported_not_gated():
    payloads, counter = _toy_pair()
    out, table = _check(_res(payloads, counter)[:1], counter, complete=False)
    assert table[0]["status"] == "not_evaluable"
    assert te.cycle_gates(table, out, complete=False) == []
    assert te.cycle_not_evaluable(table) and base_template(out[0])["status"] == "emitted"


def test_a_subset_carries_the_primary_of_a_derived_template():
    payloads, counter = _toy_pair()
    out = _res(payloads, counter)
    alone = [out[_primary_derived(out)[1]]]
    with pytest.raises(te.EmissionError, match="primary"):
        te.loaded_templates(alone, 500.0)
    carried = te.carry_owners(alone, out)
    assert len(te.loaded_templates(carried, 500.0)) == 2


def test_thermo_form_from_rmg_matches_the_provider_at_any_temperature(rmg):
    """The stored form of real RMG ThermoData reproduces
    ``GasPhaseRMGReferenceThermo`` Kc at off-grid temperatures to 1e-12."""
    import channelsTest as ct
    from rmgpy.kmc.reference_thermo import GasPhaseRMGReferenceThermo
    from rmgpy.molecule import Molecule
    prov = GasPhaseRMGReferenceThermo(ct.database().thermo, "test-only-database-commit")
    kc = te.ThermoKc(prov)

    def rows_of(smiles, offset):
        m = Molecule().from_smiles(smiles)
        m.assign_atom_ids()
        return {a + offset: (x[0], x[1], x[2], x[3], x[4], tuple((n + offset, o) for n, o in x[5]))
                for a, x in rm.molecule_rows(m).items()}
    ethane = rows_of("CC", 0)
    two = dict(rows_of("[CH3]", 100))
    two.update(rows_of("[CH3]", 200))
    form = kc.form(ethane, two, [700.0])
    assert form["kind"] == te.KC_FORM
    for T in (312.5, 437.3, 612.5, 888.8):
        assert te.kc_at(form, T) == pytest.approx(kc(ethane, two, [T])[0], rel=1e-12)


def test_c16_n_octyl_pair_balances_at_312_5_and_612_5_kelvin(rmg):
    """The reviewer's repro: the resolved C16 catalog, 2 n-octyl <-> n-hexadecane
    (record evt_7b9bec086f..., merged into evt_19427fbeec...; reverse
    evt_1417a43aae... merged into evt_f00a90e832...), loaded S_recomb / S_hom
    against 2 Kc(T) from independent RMG thermo at off-grid temperatures."""
    import channelsTest as ct
    from rmgpy.kmc.reference_thermo import GasPhaseRMGReferenceThermo
    payloads = json.loads((HERE / "fixtures" / "i147_c16_unit_payloads.json").read_text())
    kc = te.ThermoKc(GasPhaseRMGReferenceThermo(ct.database().thermo, "test-only-database-commit"))
    out, table = te.check_cycles(payloads, counter=rmg[2], kc=kc, complete=True,
                                 extra_temperatures=(312.5, 612.5))
    ids = {p["record_event_id"][:16] for p in payloads}
    assert {"evt_7b9bec086f4e", "evt_1417a43aae8f"} <= ids
    pair = next(r for r in table if r["forward"].startswith("evt_1417a43aae8f")
                or r["reverse"].startswith("evt_1417a43aae8f"))
    assert pair["status"] == "ok", pair
    for T in (312.5, 612.5):
        assert pair["extra"][T] <= 1e-12, (T, pair["extra"])
    assert max(r["max_offgrid_error"] for r in table) <= 1e-12 or all(
        r["status"] == "ok" for r in table)
    assert te.cycle_gates(table, out) == []


# --------------------------------------------------------------------------
# I-147 rework 3: self-inverse classes, merged-copy kinetics, interpolation report
# --------------------------------------------------------------------------
def test_a_self_inverse_template_is_primary_and_never_derived(monkeypatch):
    """A class that is its own reverse after merging (R3H on a linear chain,
    C11 <-> C13) finds ITSELF as the partner of its transitions: it is a
    primary with its production weight, never derived from itself.  (The real
    C25 case needs the C25 compile; the partner discovery is forced here.)"""
    payloads, counter = _toy_pair()
    resolved = te.resolve_collisions([te.finish(p) for p in payloads])[0]
    h_id, h_key = resolved[0]["record_event_id"], base_template(resolved[0])["key"]
    monkeypatch.setattr(te, "_inverse_partners", lambda *a, **k: {(h_id, h_key)})
    out = te.derive_weights(resolved, counter=counter)
    h = base_template(out[0])
    assert h.get("self_inverse") is True and "derived" not in h
    assert te.weight_at(h, 437.3) == pytest.approx(
        te._interpolate(h["weight"]["T"], h["weight"]["w"], 437.3), rel=1e-15)
    assert len(te.loaded_templates(out, 437.3)) == 2


def test_a_merged_copy_must_match_the_loaded_weight_on_its_own_kinetics():
    """After the merge the copy's own table drifts: 1e-6 passes, 3 % is gated
    and refuses the whole unit."""
    for factor, ok in ((1 + 1e-6, True), (1.03, False)):
        payloads, counter = _two_pairs(0.0)
        out, _t = te.resolve_collisions(payloads)
        out = te.derive_weights(out, counter=counter)
        copy = next(p for p in out if p["record_event_id"] == _rid(2))
        assert base_template(copy)["status"] == "merged"
        base_template(copy)["weight"]["w"] = [w * factor for w in base_template(copy)["weight"]["w"]]
        out = [te.finish(p) for p in out]
        out2, table = _check(out, counter)
        if ok:
            assert all(r["status"] == "ok" for r in table) and not te.cycle_gates(table, out2)
        else:
            assert any("own kinetics" in r["reason"] for r in table)
            assert te.cycle_gates(table, out2) and te.loaded_templates(out2, 700.0) == []


def test_the_interpolation_report_ignores_derived_tables():
    import template_emission_acceptance as acc
    payloads, counter = _toy_pair()
    out = _res(payloads, counter)
    clean = acc.interpolation_table(out)
    base_template(out[_primary_derived(out)[1]])["weight"]["w"][7] *= 1.5   # never interpolated
    rep = acc.interpolation_table(out)
    # AMENDED (I-162 S-g rework 1): the primary is now the recombination, whose
    # own table is not exact under interpolation; the derived table is ignored
    assert rep["templates"] == 1 and rep == clean


# -- class overlap -----------------------------------------------------------
def test_two_templates_of_one_family_on_one_site_are_a_class_overlap():
    w = homolysis(12, (5, 6))
    t2, t3 = cut(w, 2).matcher(1.0, "r2"), cut(w, 3).matcher(1.0, "r3")
    assert t2.family == t3.family == "toy_cc"
    rows = alkane(12)[0]
    found = te.class_overlaps([t2, t3], rows)
    assert found and all(set(x["templates"]) == {"r2", "r3"} for x in found)
    with pytest.raises(te.ClassOverlapError, match="r2"):
        te.check_class_overlap([t2, t3], rows)
    # different families producing one net edit are distinct pathways: they sum
    other = wm.Template("r3", "other_family", t3.atoms, t3.bonds, t3.edits, weight=1.0)
    assert te.class_overlaps([t2, other], rows) == []
    te.check_class_overlap([t2, other], rows)
    result = wm.match([t2, other], te.melt_state(rows))
    assert len(result.transitions) > len({tr.identity[2] for tr in result.transitions})


# -- leftovers ---------------------------------------------------------------
def test_a_subset_without_the_owner_of_a_merged_template_is_refused_or_carried():
    """C16 n-octyl recombination evt_7b9bec086f4e... is merged into the
    template of evt_19427fbeeceb...: loading it alone must not silently load
    nothing."""
    payloads = json.loads((HERE / "fixtures" / "i147_c16_orphan_payloads.json").read_text())
    by = {p["record_event_id"][:16]: p for p in payloads}
    merged_record = by["evt_7b9bec086f4e"]
    assert [t["status"] for t in merged_record["templates"]] == ["merged"]
    assert merged_record["status"] == "emitted"
    with pytest.raises(te.EmissionError, match="orphan"):
        te.loaded_templates([merged_record], 700.0)
    carried = te.carry_owners([merged_record], payloads)
    loaded = te.loaded_templates(carried, 700.0)
    assert len(loaded) == 1
    both = te.carry_owners(payloads, payloads)
    assert len(te.loaded_templates(both, 700.0)) == 1


def test_a_non_identity_request_without_a_generated_reaction_stays_refused(monkeypatch):
    """The generator returns nothing for a non-identity request: no explicit
    identity transfer is tried and the result is a refusal."""
    from rmgpy.kmc import compiler
    rates = te.ProductionRates.__new__(te.ProductionRates)
    rates.compiler = type("C", (), {"kinetics_database": None, "generation_provider": None})()
    rates._cache, rates._max, rates.calls, rates.generated = {}, 8, 0, 0
    rates.explicit_calls = 0
    rates._counter = ToyCounter()
    monkeypatch.setattr(compiler, "_resolve_generator",
                        lambda *a, **k: (lambda *x, **y: []))

    def boom(self, request):
        raise AssertionError("the explicit identity path must not run")
    monkeypatch.setattr(te.ProductionRates, "explicit", boom)
    w = homolysis(8, (3, 4))
    r, p = dict(w.reactant_rows()), dict(w.product_rows())
    req = te.RateRequest("toy_cc", "toy", "reversed", r, p, w)
    assert not rates.identity(req)
    res = rates(req)
    assert res.problems and "0 reactions" in res.problems[0]
    assert res.k == ()


# --------------------------------------------------------------------------
# I-150: a self-reverse family's channel counted on its reverse record
# --------------------------------------------------------------------------
@pytest.mark.parametrize("decades,ok", [(0.00994, True), (-0.00994, True),
                                        (0.01006, False), (-0.01006, False)])
def test_merged_copy_guard_and_collision_diameter_share_one_boundary(decades, ok):
    """Both sides of the 0.01-decade boundary, in the band the old linear
    guard (|ratio - 1| <= 10**0.01 - 1) got wrong on the low side: a copy 0.01006
    decades BELOW the loaded weight (ratio 0.97710, |ratio - 1| = 0.0229)
    passed the guard and failed the diameter.  The guard refuses exactly what
    the diameter refuses."""
    # the same two weights, as a collision diameter
    payloads, _counter = _two_pairs(0.0)
    copy = [p for p in payloads if p["record_event_id"] == _rid(2)][0]
    base_template(copy)["weight"]["w"] = [w * 10 ** decades for w in
                                         base_template(copy)["weight"]["w"]]
    payloads = [te.finish(p) if p is copy else p for p in payloads]
    for p in payloads:
        for t in p["templates"]:
            t.pop("thermo", None)
    _out, table = te.resolve_collisions(payloads, tolerance_log10=te.LOCALITY_TOLERANCE_LOG10)
    diameter_ok = all(r["disposition"] == "merged" for r in table)
    # and as the merged-copy guard
    payloads, counter = _two_pairs(0.0)
    out, _t = te.resolve_collisions(payloads)
    out = te.derive_weights(out, counter=counter)
    merged = next(p for p in out if p["record_event_id"] == _rid(2))
    base_template(merged)["weight"]["w"] = [w * 10 ** decades for w in
                                           base_template(merged)["weight"]["w"]]
    out = [te.finish(p) for p in out]
    _out2, cycles = _check(out, counter)
    guard_ok = not any("own kinetics" in r["reason"] for r in cycles)
    assert guard_ok is ok
    assert diameter_ok is ok
    assert te.beyond_tolerance(abs(decades)) is (not ok)


class _Fam:
    label = "toy_h_shift"
    reverse_template = None


class _HShiftCounter(te.RmgCounter):
    """RMG's moves replaced by a toy own-reverse family on linear alkyl
    radicals: the radical carbon takes an H from a carbon 1 or 2 bonds away
    (ring label R2H / R3H); the template is (ring, class of the radical
    carbon, class of the donor), the class being primary or secondary.  The
    reading depends on the side it is made from: primary -> secondary is
    (ring, p_rad, s_H), the same event read backwards is (ring, s_rad, p_H)."""

    def __init__(self):
        super().__init__({})

    def moves(self, family, molecules, direction):
        (mol,) = molecules
        rows = rm.molecule_rows(mol)
        carbons = [a for a, r in rows.items() if r[0] == "C"]
        nbrs = {a: [n for n, _o in rows[a][5] if rows[n][0] == "C"] for a in carbons}
        cls = {a: "p" if len(nbrs[a]) == 1 else "s" for a in carbons}
        (rad,) = [a for a in carbons if rows[a][3]]
        out = []
        for donor in carbons:
            d = abs(donor - rad)           # carbons are numbered along the chain
            if d not in (1, 2):
                continue
            for h in [n for n, _o in rows[donor][5] if rows[n][0] == "H"]:
                new = dict(rows)
                r_row, d_row, h_row = rows[rad], rows[donor], rows[h]
                new[rad] = r_row[:3] + (0,) + r_row[4:5] + (tuple(sorted(r_row[5] + ((h, 1.0),))),)
                new[donor] = d_row[:3] + (1,) + d_row[4:5] + (
                    tuple(e for e in d_row[5] if e[0] != h),)
                new[h] = h_row[:5] + (((rad, 1.0),),)
                out.append(rm.Move(
                    "toy_h_shift", "forward", (0,), (0,), (True,), ((0, "*1", rad), (0, "*3", donor)),
                    (), "kept", "", (f"R{d + 1}H", cls[rad] + "_rad", cls[donor] + "_H"),
                    tuple(sorted(new.items())), ()))
        return out


def _butyl(radical):
    """n-butyl radical rows, carbons 0..3 along the chain, the radical on ``radical``."""
    rows = dict(alkane(4)[0])
    h = next(n for n, _o in rows[radical][5] if rows[n][0] == "H")
    r_row = rows[radical]
    rows[radical] = r_row[:3] + (1,) + r_row[4:5] + (tuple(e for e in r_row[5] if e[0] != h),)
    del rows[h]
    return {a: (r[:5] + (tuple(e for e in r[5] if e[0] != h),)) for a, r in rows.items()}


def _h_shift_case():
    """Record R3H, C3 radical -> C1 radical, read from its reverse side: the
    reactant is the C1 radical, the product the primary radical (C3, and the
    mirror image C0 reached by R2H: the same species)."""
    counter = _HShiftCounter()
    reactant = _butyl(1)
    product = _butyl(3)
    species = counter.species_of_rows(product)
    pathway = ("R3H", "p_rad", "s_H")
    return counter, te.molecules_of(reactant), species, pathway


def test_the_reverse_reading_of_a_self_reverse_family_is_filtered_to_the_records_recipe():
    """R3H and R2H both reach the primary-radical species from the C1 radical
    (3 + 3 outcomes); the record's channel is the three R3H outcomes only.
    No filter takes all six."""
    counter, mols, species, pathway = _h_shift_case()
    unfiltered = counter.channel(_Fam, mols, "forward", None, species)
    assert unfiltered.transitions == 6
    flt = te.count_pathway(_Fam, pathway, True)
    assert counter.channel(_Fam, mols, "forward", flt, species).transitions == 3


def test_a_reverse_reading_on_another_node_of_the_same_recipe_still_counts():
    """The moves of the reverse reading carry (R3H, s_rad, p_H), not the
    record's (R3H, p_rad, s_H): a plain template filter finds none, which is
    why the filter is the inverse move's template."""
    counter, mols, species, pathway = _h_shift_case()
    assert counter.channel(_Fam, mols, "forward", tuple(pathway), species).transitions == 0
    flt = te.count_pathway(_Fam, pathway, True)
    assert counter.channel(_Fam, mols, "forward", flt, species).transitions == 3


def test_a_template_side_reading_and_other_families_keep_the_plain_filter():
    assert te.count_pathway(_Fam, ("a", "b"), False) == ("a", "b")
    assert te.count_pathway(_Fam, (), True) is None

    class Other:
        reverse_template = object()
    assert te.count_pathway(Other, ("a", "b"), True) == ("a", "b")


# --------------------------------------------------------------------------
# I-150 rework 1: the midpoint gate, --event subsets, the boundary's float slack
# --------------------------------------------------------------------------
class _ArrheniusRates:
    """ProductionRates stand-in: k = A exp(-E/T) of the emitted witness's channel
    (one transition) at the requested grid, or a named failure."""
    fail = None

    def __init__(self, db, families, *, database_path=None, database_sha=None,
                 temperature_grid=None, **_kw):
        self.grid = tuple(temperature_grid)

    def __call__(self, request):
        if self.fail:
            return te.RateResult((), (), 0.0, {}, (self.fail,))
        return te.RateResult(self.grid, tuple(1.0e9 * math.exp(-15000.0 / t) for t in self.grid),
                             1.0, {}, ())


def _midpoint_payload(record_count, witness_count, *, derived=False):
    """A payload whose emitted witness has ``witness_count`` transitions and
    whose record has ``record_count`` (a homolysis padded once at the head:
    the record's witness sees 2 equivalent cleavages, the emitted one 1); the
    weight table is the emitted witness's rate over its count, an exact
    Arrhenius interpolant."""
    rows = te.rows_json(alkane(4)[0])
    grid = [300.0, 400.0, 500.0, 600.0]
    w = [1.0e9 * math.exp(-15000.0 / t) / witness_count for t in grid]
    t = {"status": "emitted", "key": "k" * 64, "base": True, "weight": {"T": grid, "w": w},
         "maps": {"witnesses": {"emitted": {"r": rows, "p": rows}}},
         "locality": {"rows": [{"label": "emitted", "status": "ok",
                                "counts": {"matcher": witness_count}}]}}
    if derived:
        t["derived"] = {"primary": {"record": "evt_x", "key": "k"}}
    record = {"event_id": "evt_" + "a" * 64, "family": "toy_cc", "template": "toy",
              "orientation": "reversed"}
    payload = {"record_event_id": record["event_id"], "member_reverse": False,
               "weight_source": "record", "counts": {"record": {"matcher": record_count}},
               "templates": [t]}
    return payload, record


@pytest.fixture
def midpoint_env(monkeypatch):
    import template_emission_acceptance as acc
    monkeypatch.setattr(te, "ProductionRates", _ArrheniusRates)
    monkeypatch.setattr(_ArrheniusRates, "fail", None)
    monkeypatch.setitem(acc._S, "db", None)
    monkeypatch.setitem(acc._S, "records", [{"family": "toy_cc"}])
    return acc


def test_the_midpoint_check_divides_by_the_emitted_witness_count(midpoint_env):
    """C13 homolysis padded once at the head: record count 2, emitted-witness
    count 1.  The regenerated k is the emitted witness's, so its count (the one
    the weight was calibrated with) is the denominator: an exact interpolant
    has no midpoint error."""
    acc = midpoint_env
    payload, record = _midpoint_payload(2, 1)
    acc._S["by_id"] = {record["event_id"]: record}
    mid = acc.midpoint_table([payload])
    assert mid["templates"] == 1 and not mid.get("unverified")
    assert mid["worst"] < 1e-6
    assert acc.midpoint_gates(mid) == []


def test_an_unevaluable_loaded_primary_gates_and_a_derived_one_does_not(midpoint_env,
                                                                         monkeypatch):
    acc = midpoint_env
    payload, record = _midpoint_payload(1, 1)
    acc._S["by_id"] = {record["event_id"]: record}
    monkeypatch.setattr(_ArrheniusRates, "fail", "no generated reaction")
    mid = acc.midpoint_table([payload])
    assert mid["templates"] == 0 and len(mid["unverified"]) == 1
    gates = acc.midpoint_gates(mid)
    assert len(gates) == 1 and "no generated reaction" in gates[0] and "kkkkkkkkkkkk" in gates[0]
    derived, record = _midpoint_payload(1, 1, derived=True)
    mid = acc.midpoint_table([derived])
    assert not mid.get("unverified") and acc.midpoint_gates(mid) == []
    assert mid["unavailable"] == {"derived": 1}


def test_an_event_subset_is_not_a_complete_catalog():
    import argparse
    import template_emission_acceptance as acc
    base = dict(sample=1, chunk=None, limit=None, open_only=False)
    assert acc.is_complete(argparse.Namespace(**base, event=None)) is True
    assert acc.is_complete(argparse.Namespace(**base, event=["evt_9c441910"])) is False


def test_a_subset_of_a_reciprocal_pair_refuses_nothing_for_the_unselected_partner():
    payloads, counter = _toy_pair()
    _out, table = _check(payloads[:1], counter, complete=False)
    assert [r["status"] for r in table] == ["not_evaluable"]


def test_the_boundary_accepts_an_exact_scaling_and_refuses_beyond_it():
    """A copy scaled by exactly 10**0.01 measures 0.010000000000000675 decades
    in floating point: inside; 10**0.01001 is outside (float slack 1e-12)."""
    w = 1.2345e-9
    for sign in (1, -1):
        exact = abs(math.log10(w / (w * 10 ** (sign * 0.01))))
        assert not te.beyond_tolerance(exact)
        assert te.beyond_tolerance(abs(math.log10(w / (w * 10 ** (sign * 0.01001)))))
    assert te.beyond_tolerance(0.01 + 1e-9)


def test_an_acceptance_selection_that_matches_nothing_gates():
    """A cluster compile numbers its records differently: --event prefixes that
    match no record, or a selection that processes none, must not pass."""
    import argparse
    import template_emission_acceptance as acc
    records = [{"event_id": "evt_08c9cfdc6831aa"}, {"event_id": "evt_1561c939b769bb"}]
    args = argparse.Namespace(event=["evt_9c441910", "evt_08c9cf"])
    assert acc.selection_problems(args, records, [records[0]]) == [
        "--event prefix evt_9c441910 matches no record"]
    none = acc.selection_problems(argparse.Namespace(event=["evt_9c441910"]), records, [])
    assert any("evt_9c441910" in g for g in none) and any("processed 0 records" in g for g in none)
    assert acc.selection_problems(argparse.Namespace(event=None), records, records) == []
    assert acc.selection_problems(argparse.Namespace(event=None), records, []) == [
        "the selection processed 0 records (a vacuous run)"]


# --------------------------------------------------------------------------
# I-162 S-g rework 1: no choice follows the record ids.  A record id hashes
# the artifact's provenance (its family list among it), so the same chemistry
# compiled alone or with other families carries other ids; the primary of a
# pair, the representative of a collision and the cycle-check temperatures
# must not change with them.
# --------------------------------------------------------------------------
def _relabel_records(payloads, mapping):
    """``payloads`` with every record id (own, reverse link, rate provenance)
    renamed by ``mapping``: the same chemistry under other ids."""
    out = []
    for p in payloads:
        p = json.loads(json.dumps(p))
        p["record_event_id"] = mapping.get(p["record_event_id"], p["record_event_id"])
        if p.get("reverse_of") is not None:
            p["reverse_of"] = mapping.get(p["reverse_of"], p["reverse_of"])
        prov = (p.get("rate") or {}).get("provenance") or {}
        if prov.get("event_id") in mapping:
            prov["event_id"] = mapping[prov["event_id"]]
        out.append(te.finish(p))
    return out


def _loaded_weights(out, temperatures):
    """key -> (derived?, weights at ``temperatures``) of every loaded template
    (a derived one through its primary)."""
    index = {(p["record_event_id"], t["key"]): t for p in out for t in p["templates"]
             if t.get("status") == "emitted"}
    got = {}
    for (_rid_, key), t in index.items():
        d = t.get("derived")
        primary = index[(d["primary"]["record"], d["primary"]["key"])] if d else None
        got[key] = (bool(d), [te.weight_at(t, T, primary) for T in temperatures])
    return got


def test_the_primary_of_a_reciprocal_pair_does_not_follow_the_record_ids():
    """The pair under swapped record ids: the same template is derived from
    the same primary, and every loaded weight is the same (the reviewer's
    repro moved a loaded weight at 650 K from 1.002965894 to 1)."""
    payloads, counter = _toy_pair()
    swapped = _relabel_records(payloads, {_rid(1): _rid(2), _rid(2): _rid(1)})
    temps = list(OFF_GRID) + [650.0]
    a = _loaded_weights(_res(payloads, counter), temps)
    b = _loaded_weights(_res(swapped, counter), temps)
    assert sorted(derived for derived, _w in a.values()) == [False, True]
    assert a == b


def test_the_collision_representative_does_not_follow_the_record_ids():
    """Three copies 0.004 decades apart merge; whichever record carries which
    id, the loaded weight is the same copy's."""
    members = _group([0.0, 0.004, 0.008])
    ids = [p["record_event_id"] for p in members]
    seen = set()
    for perm in itertools.permutations(ids):
        out, table = te.resolve_collisions(_relabel_records(members, dict(zip(ids, perm))))
        assert [r["disposition"] for r in table] == ["merged"]
        loaded = [base_template(p) for p in out if base_template(p)["status"] == "emitted"]
        assert len(loaded) == 1
        seen.add(te.canonical_json(loaded[0]["weight"]))
    assert len(seen) == 1


def test_the_representative_pair_of_a_reciprocal_unit_does_not_follow_the_record_ids():
    """Two pairs merge as one unit (the reverse copies 0.003 decades apart):
    the loaded forward and reverse weights are one source pair's, the same
    pair under any ids."""
    kf, kr = "a" * 64, "b" * 64
    payloads = _reciprocal([kf, kf], [kr, kr], (0.003, 0.0))
    swap = {_rid(1): _rid(2), _rid(2): _rid(1), _rid(5): _rid(6), _rid(6): _rid(5)}
    seen = []
    for ps in (payloads, _relabel_records(payloads, swap)):
        out, table = te.resolve_collisions(ps)
        assert sorted(r["disposition"] for r in table) == ["merged", "merged"]
        seen.append(sorted((t["key"], te.canonical_json(t["weight"])) for p in out
                           for t in p["templates"] if t.get("status") == "emitted"))
    assert seen[0] == seen[1]


def test_the_cycle_check_temperatures_do_not_follow_the_record_ids(monkeypatch):
    payloads, counter = _toy_pair()
    renamed = _relabel_records(payloads, {_rid(1): _rid(10), _rid(2): _rid(11)})
    original = te.check_temperatures
    used = []

    def spy(grid, pair, extra=()):
        temps = original(grid, pair, extra)
        used[-1].update(temps)
        return temps
    monkeypatch.setattr(te, "check_temperatures", spy)
    for ps in (payloads, renamed):
        used.append(set())
        _out, table = _check(_res(ps, counter), counter)
        assert table[0]["status"] == "ok"
    assert used[0] and used[0] == used[1]


# --------------------------------------------------------------------------
# I-162 S-g rework 2: the identity is the whole stored content less the ids
# --------------------------------------------------------------------------
def test_the_identity_is_every_stored_field_but_the_record_ids():
    """Renaming the record ids leaves the identity unchanged; changing any
    other stored field of the payload or of the template changes it."""
    payloads, _counter = _toy_pair()
    p = payloads[0]
    t = base_template(p)
    renamed = _relabel_records(payloads, {_rid(1): _rid(12), _rid(2): _rid(13)})[0]
    assert te.partition_free_identity(renamed) == te.partition_free_identity(p)
    assert te.partition_free_identity(renamed, base_template(renamed)) == \
        te.partition_free_identity(p, t)
    stripped = {"record_event_id", "reverse_of", "sha256", "certificate"}
    for field_ in sorted(set(p) - stripped):
        q = json.loads(json.dumps(p))
        q[field_] = {"changed": field_}
        assert te.partition_free_identity(q) != te.partition_free_identity(p), field_
    q = json.loads(json.dumps(p))
    q["rate"]["k"] = [2 * k for k in q["rate"]["k"]]
    assert te.partition_free_identity(q) != te.partition_free_identity(p)
    for field_ in sorted(t):
        q = json.loads(json.dumps(p))
        base_template(q)[field_] = {"changed": field_}
        assert te.partition_free_identity(q, base_template(q)) != \
            te.partition_free_identity(p, t), field_


def _thermo_variant_pairs():
    """Pair A (records 1/6) and pair B (2/5), equal in everything but the
    stored thermo of B's recombination (ssa_ratio 0.1 % larger)."""
    payloads, counter = _two_pairs(0.0)
    w_b = next(p for p in payloads if p["record_event_id"] == _rid(5))
    th = base_template(w_b)["thermo"]
    th["ssa_ratio"] *= 1.001
    payloads = [te.finish(p) for p in payloads]
    return payloads, counter


def test_copies_that_differ_only_in_thermo_choose_by_thermo_never_by_id():
    """The spar's repro: two source pairs tie on everything the earlier
    identity read, so the record ids chose the representative recombination,
    and the derived homolysis took that copy's thermo.  Under any ids the
    loaded weights (and the representative's thermo) are the same."""
    payloads, counter = _thermo_variant_pairs()
    swap = {_rid(1): _rid(2), _rid(2): _rid(1), _rid(5): _rid(6), _rid(6): _rid(5)}
    temps = list(OFF_GRID) + [650.0]
    seen = []
    for ps in (payloads, _relabel_records(payloads, swap)):
        out = _res(ps, counter)
        thermo = sorted(te.canonical_json(t.get("thermo")) for p in out for t in p["templates"]
                        if t.get("status") == "emitted")
        seen.append((_loaded_weights(out, temps), thermo))
    assert seen[0] == seen[1]


def test_the_as_generated_side_is_primary_even_with_the_greater_key():
    """The recombination's key is the lesser in the toy pair; with the
    orientations exchanged the as_generated homolysis (the greater key) must
    be the primary."""
    payloads, counter = _toy_pair()
    for p in payloads:
        p["orientation"] = "as_generated" if p["orientation"] == "reversed" else "reversed"
    payloads = [te.finish(p) for p in payloads]
    out = _res(payloads, counter)
    ip, idd = _primary_derived(out)
    assert out[ip]["orientation"] == "as_generated"
    assert base_template(out[ip])["key"] > base_template(out[idd])["key"]


# --------------------------------------------------------------------------
# I-162 S-g rework 3: only the id slots and the id forms the code writes are
# masked; every other string is kept as stored
# --------------------------------------------------------------------------
def _with_selection_entry(payload, label):
    p = json.loads(json.dumps(payload))
    p["rate"]["provenance"].setdefault("selection", {})["entry"] = label
    return te.finish(p)


def test_hex_like_labels_are_not_record_ids():
    payloads, _counter = _toy_pair()
    p = payloads[0]
    a, b = (json.loads(json.dumps(p)) for _ in range(2))
    a["pathway"], b["pathway"] = ["evt_dead"], ["evt_beef"]
    assert te.partition_free_identity(a) != te.partition_free_identity(b)
    a, b = _with_selection_entry(p, "evt_dead"), _with_selection_entry(p, "evt_beef")
    assert te.partition_free_identity(a) != te.partition_free_identity(b)
    assert te.partition_free_identity(a, base_template(a)) != \
        te.partition_free_identity(b, base_template(b))


def test_free_text_masks_exactly_the_id_forms_the_code_writes():
    """A refusal reason quoting the record's own id (whole, or its 16-character
    ``str(id)[:16]`` form of check_cycles / one_way_problems) is masked; a
    reason quoting a hex-like label is kept."""
    payloads, _counter = _toy_pair()
    p = payloads[0]

    def with_reason(payload, text):
        q = json.loads(json.dumps(payload))
        q["refusal"] = {"status": "refused:cycle", "reason": text}
        return q
    rid = p["record_event_id"]
    renamed = _relabel_records([p], {rid: _rid(14)})[0]
    for form in (lambda r: r, lambda r: r[:16]):
        assert te.partition_free_identity(with_reason(p, f"reverse record {form(rid)} x")) == \
            te.partition_free_identity(with_reason(renamed, f"reverse record {form(_rid(14))} x"))
    assert te.partition_free_identity(with_reason(p, "rule evt_dead")) != \
        te.partition_free_identity(with_reason(p, "rule evt_beef"))


def test_the_retained_rule_does_not_follow_the_ids_when_labels_look_like_ids():
    """The spar's repro: two copies equal but for a rate-entry label that
    looks like an id (evt_dead, evt_beef).  Masking every hex run made them
    tie, so the record ids chose which rule the catalog kept."""
    members = _group([0.0, 0.0])
    members = [_with_selection_entry(members[0], "evt_dead"),
               _with_selection_entry(members[1], "evt_beef")]
    ids = [p["record_event_id"] for p in members]
    kept = set()
    for mapping in ({}, dict(zip(ids, reversed(ids)))):
        out, table = te.resolve_collisions(_relabel_records(members, mapping))
        assert [r["disposition"] for r in table] == ["merged"]
        (rep,) = [p for p in out if base_template(p)["status"] == "emitted"]
        kept.add(rep["rate"]["provenance"]["selection"]["entry"])
    assert len(kept) == 1


# --------------------------------------------------------------------------
# I-188: partial payloads, duplicate reciprocal copies, T_ch = 0 off the record
# --------------------------------------------------------------------------
I188_PAYLOADS = HERE / "fixtures" / "i188_reciprocal_partial_payloads.json"


def _i188_payloads():
    """Four records of a full emit: two reciprocal record pairs whose base
    templates are the same two executables (each forward copy on both
    forward records, each reverse copy on both reverse records), two of the
    records ``partial`` with a base template refused at emission (no key)."""
    return json.loads(I188_PAYLOADS.read_text())


def _keyless(payloads):
    return [t for p in payloads for t in p["templates"]
            if t["status"] == "refused" and "key" not in t]


def _emitted_keys(payloads):
    from collections import Counter
    return Counter(t["key"] for p in payloads for t in p["templates"] if t["status"] == "emitted")


def test_a_keyless_refused_template_of_a_partial_payload_resolves_by_status():
    """A template refused at emission was never cut, so it carries no key and
    takes no part in a collision unit; resolving a catalog that holds one
    reads its status, not its missing key."""
    payloads = _i188_payloads()
    refused = _keyless(payloads)
    assert refused and all(t["base"] for t in refused)
    assert {p["status"] for p in payloads} == {"partial", "emitted"}
    out, table = te.resolve_collisions(payloads)
    assert _keyless(out) == refused
    assert table and {c["disposition"] for c in table} == {"merged"}


def test_duplicate_copies_of_a_reciprocal_template_merge_to_one_loaded_template():
    payloads = _i188_payloads()
    copies = _emitted_keys(payloads)
    duplicated = {k for k, n in copies.items() if n > 1}
    assert len(duplicated) == 2
    out, table = te.resolve_collisions(payloads)
    assert {c["key"] for c in table} == duplicated
    loaded = _emitted_keys(out)
    assert set(loaded) == set(copies) and set(loaded.values()) == {1}
    merged = [t for p in out for t in p["templates"] if t["status"] == "merged"]
    assert sorted(t["key"] for t in merged) == sorted(duplicated)


@pytest.fixture(scope="module")
def i188_balance(rmg):
    """The fixture family's RMG counter and the independent thermo Kc."""
    import channelsTest as ct
    from rmgpy.kmc.reference_thermo import GasPhaseRMGReferenceThermo
    return rmg[2], te.ThermoKc(GasPhaseRMGReferenceThermo(ct.database().thermo,
                                                          "test-only-database-commit"))


def _derived_pairs(payloads):
    """{(record, key) of a derived template: (record, key) of its primary}."""
    return {(p["record_event_id"], t["key"]): (t["derived"]["primary"]["record"],
                                               t["derived"]["primary"]["key"])
            for p in payloads for t in p["templates"]
            if t["status"] == "emitted" and t.get("derived")}


def test_merged_reciprocal_copies_derive_one_reverse_from_one_primary(i188_balance):
    """After merging, every loaded template of a reciprocal record pair has
    exactly one inverse partner: one side is primary, the other derived
    from it, as a bijection between the two records' loaded templates."""
    counter, _kc = i188_balance
    out = te.derive_weights(te.resolve_collisions(_i188_payloads())[0], counter=counter)
    by_record = {p["record_event_id"]: p for p in out}
    pairs = _derived_pairs(out)
    assert pairs
    loaded = {(p["record_event_id"], t["key"]) for p in out for t in p["templates"]
              if t["status"] == "emitted"}
    primaries = list(pairs.values())
    assert len(set(primaries)) == len(primaries)
    assert set(pairs) | set(primaries) == loaded and not set(pairs) & set(primaries)
    for (rid, _key), (prid, _pkey) in pairs.items():
        assert by_record[rid]["reverse_of"] == prid


def test_a_template_with_several_inverse_partners_stays_independent(i188_balance):
    """Unmerged, each forward copy is undone by two identical reverse
    copies: no single primary exists, so neither direction is derived."""
    counter, _kc = i188_balance
    payloads = _i188_payloads()
    duplicated = {k for k, n in _emitted_keys(payloads).items() if n > 1}
    out = te.derive_weights(payloads, counter=counter)
    assert not [key for _rid, key in _derived_pairs(out) if key in duplicated]


def test_a_derived_reverse_balances_off_the_grid_within_the_cycle_tolerance(i188_balance):
    """The merged catalog passes its own cycle check at the grid, the
    midpoints and the random off-grid temperatures at ``CYCLE_RTOL``: the
    derived weight evaluates Kc exactly at T."""
    counter, kc = i188_balance
    out = te.derive_weights(te.resolve_collisions(_i188_payloads())[0], counter=counter)
    out, table = te.check_cycles(out, counter=counter, kc=kc, complete=True,
                                 extra_temperatures=(312.5, 412.5, 612.5))
    assert len(table) == 2 and {r["status"] for r in table} == {"ok"}, table
    assert max(r["max_offgrid_error"] for r in table) <= te.CYCLE_RTOL
    assert max(r["max_ratio_error"] for r in table) <= te.CYCLE_RTOL
    assert te.cycle_gates(table, out) == []


class _ChannelUpTo(ToyCounter):
    """The toy oracle whose channel is empty on a witness of more than ``n``
    carbons: RMG reads no move of the record's channel there."""

    def __init__(self, n):
        super().__init__()
        self.n = n

    def channel(self, family, molecules, direction, pathway, product_species):
        carbons = sum(1 for m in molecules for x in rm.molecule_rows(m).values() if x[0] == "C")
        if carbons > self.n:
            self.calls += 1
            return te.RmgChannel((), 0, (), ())
        return super().channel(family, molecules, direction, pathway, product_species)


def test_a_witness_outside_the_record_channel_is_a_compile_error():
    """RETRACTED (I-188 read T_ch = 0 with rmg_moves 0 as a witness outside the
    channel and refused the template): the catalog rows that showed it were
    candidates counted in another member's orientation.  T_ch = 0 stays a
    compile error, whatever RMG counts."""
    with pytest.raises(te.EmissionError, match="T_ch = 0 on witness 'port 0:head"):
        emit(homolysis(10, (3, 4)), counter=_ChannelUpTo(10))


def _selecting(witness, pathway):
    """The toy member whose fixed-point certificate selects ``pathway`` (the
    member channel of its window's own recertification)."""
    from dataclasses import replace
    from types import SimpleNamespace
    case, row, ((ctx,), merges) = member(witness)
    cert = {"channels": [{"pathway": list(pathway), "members": []}]}
    return case, row, ((replace(ctx, observation=SimpleNamespace(certificate=cert)),), merges)


def test_a_fixed_point_selected_through_another_node_is_refused_by_name():
    """A proxy's artificial cap next to a labelled atom can make RMG read a
    node that the real chain never reaches: grown back (the fixed point keeps
    a recertification flip), the window selects a deeper node and the
    record's template makes no transition there (the PE 1,3-H shift of
    oct-4-yl, Cs_H_out_H/NonDeC on the proxy, Cs_H_out_H/(NonDeC/Cs) on
    undecan-5-yl).  The candidate is refused by name before any count, so
    the record is refused, not an emission error; a window that still
    selects the member's node is a candidate, and T_ch = 0 on it stays a
    compile error."""
    w = homolysis(10, (4, 5))
    case, _row, (results, merges) = _selecting(w, ("deeper",))
    (cand,), _rows = te.context_candidates(case, results, merges, family="toy_cc",
                                           pathway=("toy",), direction="reverse", reverse=False)
    assert cand.refusal and cand.refusal.startswith("refused:selection_flip"), cand.refusal
    assert "deeper" in cand.refusal and "toy" in cand.refusal
    out = te.emit_record(toy_record(w, 1), family=ToyFamily,
                         members=[_selecting(w, ("deeper",))], counter=_ChannelUpTo(0),
                         rates=stub_rates(ToyCounter()), check_upstream=False)
    assert out["status"] == "refused"
    assert "refused:selection_flip" in json.dumps(out["refusal"]), out["refusal"]
    case, _row, (results, merges) = _selecting(w, ("toy",))
    (cand,), _rows = te.context_candidates(case, results, merges, family="toy_cc",
                                           pathway=("toy",), direction="reverse", reverse=False)
    assert cand.refusal is None
    with pytest.raises(te.EmissionError, match="T_ch = 0"):
        te.emit_record(toy_record(w, 1), family=ToyFamily, members=[_selecting(w, ("toy",))],
                       counter=_ChannelUpTo(0), rates=stub_rates(ToyCounter()),
                       check_upstream=False)


def test_each_candidate_is_counted_in_its_own_members_orientation():
    """One record, two members of opposite orientation (its own homolysis
    witness and the recombination that forms the same bond): each candidate
    is counted on its member's witness read in that member's orientation, in
    either member order."""
    h, r = homolysis(12, (3, 4)), recombination(4, 8)
    for members in ([member(h), member(r, direction="forward")],
                    [member(r, direction="forward"), member(h)]):
        out = te.emit_record(toy_record(h, 1), family=ToyFamily, members=members,
                             counter=ToyCounter(), rates=_pair_rates({"reversed": K_RULE}),
                             check_upstream=False)
        assert out["status"] == "emitted", out["templates"]
        t = base_template(out)
        assert len(t["sources"]) == 2
        assert {r["status"] for r in t["locality"]["rows"]} <= {"ok", "not_applicable"}


# --------------------------------------------------------------------------
# I-190: one executable is one key; a refused copy never stands for it
# --------------------------------------------------------------------------
def test_the_canonical_executable_is_blind_to_the_metadata():
    """A five-carbon path: marking one atom next to the centre in the
    metadata (a boundary label, an RMG label) must not change which executable
    serialization is canonical, so equal executables get one key."""
    edges = {(0, 1): "1", (1, 2): "1", (0, 3): "1", (3, 4): "1"}
    plain = ["m"] * 5
    for v in range(5):
        marked = list(plain)
        marked[v] = "physical"
        assert te._canonical_order(5, ["C"] * 5, marked, edges)[1][0] == \
            te._canonical_order(5, ["C"] * 5, plain, edges)[1][0], v


def _record_candidates(witness):
    case, _row, (results, merges) = member(witness)
    return te.context_candidates(case, results, merges, family="toy_cc", pathway=("toy",),
                                 direction="reverse", reverse=False)[0]


def test_one_executable_has_one_key_whatever_its_port_reaches():
    """C10 bond 4-5 cut with the record's RMG labels on two witnesses whose
    windows are equal and whose head port reaches a physical end on one, the
    artificial end on the other: one executable, one key."""
    (phys,) = _record_candidates(homolysis(10, (4, 5), head="methyl"))
    (art,) = _record_candidates(homolysis(10, (4, 5)))
    ex_p, ex_a = phys.cut.executable(), art.cut.executable()
    assert sorted(map(json.dumps, ex_p["atoms"])) == sorted(map(json.dumps, ex_a["atoms"]))
    assert len(ex_p["bonds"]) == len(ex_a["bonds"])
    assert phys.refusal.startswith("refused:whole_species_guard") and art.refusal is None
    assert phys.cut.key == art.cut.key


def _two_member_emit(witnesses):
    plain = homolysis(10, (4, 5))
    oracle = ToyCounter().channel(ToyFamily, te.molecules_of(dict(plain.reactant_rows())),
                                  "reverse", None, ToyCounter().species_of_rows(plain.product_rows()))
    return te.emit_record(toy_record(plain, oracle.transitions), family=ToyFamily,
                          members=[member(w) for w in witnesses], counter=ToyCounter(),
                          rates=stub_rates(ToyCounter()), check_upstream=False)


def _partition(payload):
    return sorted((t["status"], t.get("key"), t["base"], sorted(map(tuple, t["sources"])))
                  for t in payload["templates"])


def test_a_refused_copy_does_not_stand_for_an_accepted_executable():
    """The whole-species guard reads the candidate's own witness (its port
    reaches a physical end); an equal executable accepted on another witness
    is the template, and the refused copy is one more of its sources."""
    out = _two_member_emit([homolysis(10, (4, 5), head="methyl"), homolysis(10, (4, 5))])
    assert out["status"] == "emitted", out["refusal"]
    (t,) = out["templates"]
    assert t["status"] == "emitted" and t["base"]
    assert len(t["sources"]) == 2


def test_the_member_order_does_not_change_the_emitted_partition():
    """The two records of a reciprocal pair draw their candidates from the
    same members; the templates, their statuses and sources must not follow
    the order the members come in."""
    witnesses = [homolysis(10, (4, 5), head="methyl"), homolysis(10, (4, 5))]
    one, two = _two_member_emit(witnesses), _two_member_emit(witnesses[::-1])
    assert _partition(one) == _partition(two)
    assert (one["status"], one["refusal"]) == (two["status"], two["refusal"])


# --------------------------------------------------------------------------
# I-192: an own-reverse record with members of both orientations; no
# non-finite spread; an order-free representative
# --------------------------------------------------------------------------
class _ShiftCounter(_HShiftCounter):
    def __init__(self):
        super().__init__()
        self.families = {_Fam.label: _Fam}


SHIFT = ("R3H", "p_rad", "s_H")


def _shifted(rows, rad, donor):
    """``rows`` with one H moved from carbon ``donor`` to the radical ``rad``."""
    rows = dict(rows)
    h = next(n for n, _o in rows[donor][5] if rows[n][0] == "H")
    r_row, d_row, h_row = rows[rad], rows[donor], rows[h]
    rows[rad] = r_row[:3] + (r_row[3] - 1,) + r_row[4:5] + (tuple(sorted(r_row[5] + ((h, 1.0),))),)
    rows[donor] = d_row[:3] + (d_row[3] + 1,) + d_row[4:5] + (
        tuple(e for e in d_row[5] if e[0] != h),)
    rows[h] = h_row[:5] + (((rad, 1.0),),)
    return rows


def _shift_witness(primary_side, n=8):
    """The C``n`` 1,3 H shift between the primary radical at the tail and the
    secondary radical two carbons in, generated on the primary side (the
    record's recipe ``SHIFT`` is read there) or on the secondary side (RMG
    reads ``(R3H, s_rad, p_H)`` there)."""
    alkane_ = b.PROVIDERS["linear_alkane"]
    if primary_side:
        build, rad, donor = alkane_.build(n, tail_end="primary_radical"), n - 1, n - 3
    else:
        build = alkane_.build(n, tail_end="methyl", interior_radical=n - 3)
        rad, donor = n - 3, n - 1
    mols, parts = b.witness_reactants([build])
    rows = rm.molecule_rows(mols[0])
    return b.make_witness(tuple(sorted(rows.items())),
                          tuple(sorted(_shifted(rows, rad, donor).items())), parts,
                          label=f"C{n}|{'primary' if primary_side else 'secondary'}")


def _shift_member(witness, label):
    w = window(witness, 2)
    pair = sorted(a for a in b.edit_support(witness) if witness.reactant_rows()[a][0] == "C")
    move = Move(((0, "*1", pair[0]), (0, "*3", pair[1])), (("*1", pair[0]), ("*3", pair[1])))
    case = win.MemberCase(label, _Fam.label, "forward", SHIFT, witness, (move,), frozenset())
    ctx = win.ContextResult(tuple((p, None) for p in b.vector_ports(witness)), "fixed_point", "",
                            witness, w, None)
    return case, {"status": "fixed_point"}, ((ctx,), ())


class _ShiftKc:
    """The toy thermo of the shift: the primary radical TOY_B, the secondary TOY_C."""

    def form(self, reactant_rows, product_rows, temperatures):
        rows = dict(reactant_rows)
        primary = any(r[0] == "C" and r[3] and sum(rows[n][0] == "C" for n, _o in r[5]) == 1
                      for r in rows.values())
        return {"kind": "ThermoData-Kc/1", "P0": 1e5,
                "reactants": [TOY_B if primary else TOY_C],
                "products": [TOY_C if primary else TOY_B]}

    def __call__(self, reactant_rows, product_rows, temperatures):
        f = self.form(reactant_rows, product_rows, temperatures)
        return [te.kc_at(f, t) for t in temperatures]


def _shift_pair(members):
    """The as_generated shift P (primary -> secondary, its recipe read on its
    reactants) and its reversed record Q (read on Q's products), both drawn
    from ``members``; ``k_Q = k_P / Kc`` so the pair balances."""
    pw, sw = _shift_witness(True), _shift_witness(False)
    kc = _ShiftKc()
    k_p = tuple(K_RULE)
    k_q = tuple(k / c for k, c in zip(k_p, kc(pw.reactant_rows(), pw.product_rows(), GRID)))
    out = []
    for witness, orientation, rid, rev, k in ((pw, "as_generated", _rid(1), _rid(2), k_p),
                                              (sw, "reversed", _rid(2), _rid(1), k_q)):
        record = toy_record(witness, 1, orientation=orientation)
        record.update(family=_Fam.label, template=";".join(SHIFT), event_id=rid, reverse_of=rev)
        record["k_table"]["k"] = list(k)
        out.append(te.emit_record(record, family=_Fam, members=members, counter=_ShiftCounter(),
                                  rates=_pair_rates({"as_generated": k_p, "reversed": k_q}),
                                  check_upstream=False, kc=kc))
    return out, kc


def _mixed_shift_members():
    """One member generated on each side; the secondary-side one sorts first."""
    return [_shift_member(_shift_witness(True), "b|primary"),
            _shift_member(_shift_witness(False), "a|secondary")]


def test_an_own_reverse_record_counts_its_recipe_on_its_template_side():
    """An own-reverse family whose readings differ by side, one record pair,
    members generated on both sides.  Every count of a record's channel is on
    its reactants in its direction, whichever member the witness came from:
    the recipe filter follows the record's template side (plain for P,
    through the inverse for Q), never a member's.  Both records emit in either
    member order, and the loaded pair passes the cycle gate on and off the
    grid."""
    for members in (_mixed_shift_members(), _mixed_shift_members()[::-1]):
        payloads, kc = _shift_pair(members)
        assert [p["status"] for p in payloads] == ["emitted", "emitted"], \
            [p.get("refusal") for p in payloads]
        out, table = te.check_cycles(_res(payloads, _ShiftCounter()), counter=_ShiftCounter(),
                                     kc=kc)
        assert len(table) == 1 and table[0]["status"] == "ok", table
        assert table[0]["max_offgrid_error"] <= te.CYCLE_RTOL
        assert te.cycle_gates(table, out) == []


def test_each_template_carries_its_own_members_orientation():
    """The template keeps the orientation of the member it was cut from (its
    witnesses' generation side), not the record's last member's."""
    payloads, _kc = _shift_pair(_mixed_shift_members())
    for p in payloads:
        for t in p["templates"]:
            secondary = t["member"].endswith("secondary")
            want = secondary if p["orientation"] == "as_generated" else not secondary
            assert t["member_reverse"] is want, (p["orientation"], t["member"])
        assert "member_reverse" not in p


def test_a_candidate_with_no_ok_locality_row_is_refused_by_name_and_serializes():
    """No row ok and no record weight (every count mismatched): the spread is
    unmeasured, written null, the candidate refused by name, and the payload
    seals (a non-finite spread used to make ``finish`` raise)."""
    w = homolysis(10, (4, 5))
    (cand,) = _record_candidates(w)
    entry = te.emit_candidate(cand, record=toy_record(w, 1), family=ToyFamily,
                              direction="reverse", pathway=("toy",), cpath=("toy",),
                              reverse=False, counter=ToyCounter(extra=1),
                              rates=stub_rates(ToyCounter()), weights=None,
                              temperatures=list(GRID))
    assert entry["status"] == "refused" and entry["refusal"]["status"]
    loc = entry["locality"]
    assert loc["max_spread_log10"] is None and set(loc["spread_log10"]) == {None}
    assert loc["accepted"] is False
    sealed = te.finish({"schema": te.SCHEMA, "status": "refused",
                        "refusal": {"status": "refused:base_template"}, "templates": [entry]})
    te.validate_payload(sealed)


def test_an_unmeasured_or_non_finite_spread_is_null():
    """No ok row measures nothing, and an infinite weight measures nothing at
    its temperature: null there and as the maximum, never inf or NaN."""
    per, worst = te.spread([], list(GRID))
    assert worst is None and set(per) == {None}
    assert te.spread([{"w": [1.0, float("inf")]}, {"w": [1.0, 10.0]}], [1.0, 2.0]) == \
        ([0.0, None], None)
    assert te.spread([{"w": [1.0, 1.0]}, {"w": [10.0, 1.0]}], [1.0, 2.0]) == ([1.0, 0.0], 1.0)


def test_reversing_the_member_order_gives_a_byte_identical_payload():
    """Two accepted copies of one executable on different witnesses: which
    copy represents the template (its witnesses, weights and maps) is chosen
    by a rule on the candidates themselves, and the payload's lists are
    ordered by content, so the member order leaves no trace."""
    witnesses = [homolysis(10, (3, 4)), homolysis(10, (5, 6))]
    oracle = ToyCounter().channel(ToyFamily, te.molecules_of(dict(witnesses[0].reactant_rows())),
                                  "reverse", None,
                                  ToyCounter().species_of_rows(witnesses[0].product_rows()))

    def emit_(ws):
        return te.emit_record(toy_record(witnesses[0], oracle.transitions), family=ToyFamily,
                              members=[member(w) for w in ws], counter=ToyCounter(),
                              rates=stub_rates(ToyCounter()), check_upstream=False)
    one, two = emit_(witnesses), emit_(witnesses[::-1])
    assert one["status"] == "emitted"
    (t,) = one["templates"]
    assert len(t["sources"]) == 2
    assert te.canonical_json(one) == te.canonical_json(two)
    mixed = [_shift_pair(m)[0] for m in (_mixed_shift_members(), _mixed_shift_members()[::-1])]
    assert te.canonical_json(mixed[0]) == te.canonical_json(mixed[1])


# --------------------------------------------------------------------------
# I-193: a degenerate identity exchange keeps the record's own orientation
# --------------------------------------------------------------------------
DEGENERATE = ("R3H", "s_rad", "s_H")


def _degenerate_shift_witness(n=15):
    """The C``n`` (odd ``n``, both ends artificial) 1,3 H shift between the
    two secondary carbons one off the middle (n = 15: carbon 8 takes an H from
    carbon 6, the window clear of both ends): reactant and product are one
    species, so only the atom map tells the two sides apart."""
    alkane_ = b.PROVIDERS["linear_alkane"]
    rad, donor = (n + 1) // 2, (n - 3) // 2
    mols, parts = b.witness_reactants([alkane_.build(n, interior_radical=rad)])
    rows = rm.molecule_rows(mols[0])
    return b.make_witness(tuple(sorted(rows.items())),
                          tuple(sorted(_shifted(rows, rad, donor).items())), parts,
                          label=f"C{n}|degenerate")


def _channel_rates(counter, pathway):
    """Toy production rates as RMG's path degeneracy makes them: ``K_RULE`` per
    transition of the witness's channel read in the record's own direction (a
    reversed record's k_f / Kc carries the symmetry numbers, so its degeneracy
    is the one on its own reactants, the generation's products); a symmetric
    extension doubles it."""
    def rates(req):
        source, target = (req.products, req.generation) if req.orientation == "reversed" else \
            (req.generation, req.products)
        target = counter.species_of_rows(target)
        n = sum(tuple(m.template) == pathway and counter.species_of_rows(m.product_rows()) == target
                for m in counter.moves(_Fam, te.molecules_of(source), "forward"))
        return te.RateResult(GRID, tuple(k * n for k in K_RULE), 1.0, {}, ())
    return rates


def test_a_degenerate_exchange_keeps_the_records_own_orientation():
    """An own-reverse family's member is generated on the record's template
    side: the as_generated record's reactants, the reversed record's
    products.  When both sides of the record are one species the species
    cannot say which side the member's witness stands on; the record's own
    orientation does, and the template carries it."""
    counter = _ShiftCounter()
    w = _degenerate_shift_witness()
    assert counter.species_of_rows(dict(w.reactant_rows())) == \
        counter.species_of_rows(dict(w.product_rows()))
    for orientation in ("as_generated", "reversed"):
        record = toy_record(w, 2, orientation=orientation)    # one donor carbon, two H
        record.update(family=_Fam.label, template=";".join(DEGENERATE), event_id=_rid(1))
        if orientation == "reversed":
            # the template side (the member's witness reactants) is the record's products
            record["reactant_graphs"], record["product_graphs"] = \
                record["product_graphs"], record["reactant_graphs"]
        out = te.emit_record(record, family=_Fam, members=[_shift_member(w, "a|degenerate")],
                             counter=counter, rates=_channel_rates(counter, DEGENERATE),
                             check_upstream=False)
        assert out["status"] == "emitted", out.get("refusal")
        assert {t["member_reverse"] for t in out["templates"]} == {orientation == "reversed"}


# --------------------------------------------------------------------------
# depository overrides (I-211): a certified channel whose analogues are
# depository reactions at finite sizes fires there at RMG's own production
# rate (the depository entry; its reverse through Kc), the rule elsewhere
# --------------------------------------------------------------------------
class ToyDepository:
    label = "toy_cc/training"

    def __init__(self, entries):
        self.entries = {e.index: e for e in entries}

    def __str__(self):
        return f"Kinetics Depository {self.label}"


def toy_entry(index, reactants, products, rank=3):
    from types import SimpleNamespace
    from rmgpy.reaction import Reaction
    from rmgpy.species import Species
    return SimpleNamespace(index=index, rank=rank, data=None,
                           label=" + ".join(reactants) + " <=> " + " + ".join(products),
                           item=Reaction(reactants=[Species().from_smiles(x) for x in reactants],
                                         products=[Species().from_smiles(x) for x in products]))


class DepositoryCounter(ToyCounter):
    """The toy oracle on a toy family that has a depository."""

    def __init__(self, depositories, **kw):
        super().__init__(**kw)
        self._family = type("ToyDepositoryFamily", (ToyFamily,),
                            {"depositories": list(depositories)})

    def family(self, name):
        return self._family


# n-pentyl + n-octyl -> n-tridecane: one size the C12 pair's templates reach
C13_R, C13_P = ["[CH2]CCCC", "[CH2]CCCCCCC"], ["CCCCCCCCCCCCC"]


def _c13_match(index=7):
    from rmgpy.molecule import Molecule
    entry = toy_entry(index, C13_R, C13_P)
    return entry, {"depository": "toy_cc/training", "index": index, "label": entry.label,
                   "orientation": "forward", "rank": 3, "kinetics": "None",
                   "reaction": {"reactants": [Molecule(smiles=x).to_adjacency_list() for x in C13_R],
                                "products": [Molecule(smiles=x).to_adjacency_list() for x in C13_P]}}


def _certificate(*matches):
    return {"certified": True, "certificate": {"channels": [{"depository_exclusion": {
        "verdict": "excluded:except_matched", "matched": [f"toy_cc/training#{m['index']}"
                                                          for m in matches],
        "matches": list(matches)}}]}}


K_DEP = tuple(30.0 * k for k in K_RULE)      # the toy depository rate (not the rule's)


def _depository_rates(counter, base, k_dep=K_DEP, reactants=None, products=None):
    """Production rates: the C13 reaction (generated on its recombination
    side; or ``reactants`` -> ``products``) from the toy depository, its
    reverse k_dep / Kc; anything else ``base``."""
    reactants, products = reactants or C13_R, products or C13_P
    a = counter.species_of_rows(_rows(reactants))
    fwd = {"kind": "RMG family estimate", "source": "Kinetics Depository toy_cc/training",
           "entry": " + ".join(reactants) + " <=> " + " + ".join(products), "rank": 3,
           "family_template_direction": "forward"}

    def rates(req):
        if counter.species_of_rows(req.generation) != a:
            return base(req)
        if req.orientation == "as_generated":
            return te.RateResult(GRID, k_dep, 1.0, {"selection": fwd}, ())
        return te.RateResult(GRID, tuple(k / c for k, c in zip(k_dep, _kc())), 1.0,
                             {"selection": {"kind": "reference-thermo reverse",
                                            "reference": "k_family / Kc", "forward": fwd}}, ())
    return rates


def _rows(smiles):
    from rmgpy.molecule import Molecule
    rows, k = {}, 0
    for x in smiles:
        m = Molecule(smiles=x)
        for a in m.atoms:
            a.id, k = k, k + 1
        rows.update(rm.molecule_rows(m))
    return rows


def test_a_locality_witness_rated_from_a_depository_is_a_named_row():
    """An extension witness of the C12 homolysis (n-tridecane -> n-pentyl +
    n-octyl) IS a depository reaction: RMG rates it from the depository, so
    it is no locality row of the rule weight.  It is named (``depository``,
    the entry) and leaves the template emitted on the other rows; nothing
    else changes."""
    entry, _match = _c13_match()
    counter = DepositoryCounter([ToyDepository([entry])])
    payload = emit(homolysis(12, (3, 4)), counter=counter)
    t = base_template(payload)
    dep = [r for r in t["locality"]["rows"] if r["status"] == "depository"]
    assert len(dep) == 1, [(r["label"], r["status"]) for r in t["locality"]["rows"]]
    assert dep[0]["depository"]["index"] == 7 and dep[0]["counts"]["matcher"] == 2
    assert t["status"] == "emitted" and t["locality"]["accepted"]
    plain = base_template(emit(homolysis(12, (3, 4))))
    assert t["weight"] == plain["weight"] and t["key"] == plain["key"]
    assert [r for r in plain["locality"]["rows"] if r["status"] == "depository"] == []


def _override_pair():
    entry, match = _c13_match()
    counter = DepositoryCounter([ToyDepository([entry])])
    payloads, counter = _toy_pair(counter=counter, certificate=_certificate(match))
    resolved = _res(payloads, counter)
    rates = _depository_rates(counter, _pair_rates({"reversed": K_RULE, "as_generated": K_RULE}))
    out, table, failures = te.depository_overrides(resolved, counter=counter, rates=rates,
                                                   kc=toy_kc)
    by = {p["orientation"]: base_template(p) for p in out}
    return out, table, failures, by["reversed"], by["as_generated"]


def test_a_matched_size_takes_the_depository_rate_and_its_reverse_through_kc():
    """At catalog assembly the certificate's matched analogue (n-pentyl +
    n-octyl -> n-tridecane) is rated by RMG's production path: the
    recombination template, whose one transition makes it, carries
    w = ssa k_depository / T; the homolysis template, two transitions on
    n-tridecane, carries the reverse derived through the stored Kc, so the
    pair balances at any temperature (S_f / S_r = ssa_ratio Kc) and equals
    the production reverse (k_depository / Kc) on the grid."""
    out, table, failures, h, w = _override_pair()
    assert failures == [], failures
    (ow,), (oh,) = w["overrides"], h["overrides"]
    assert ow["counts"] == {"matcher": 1, "rmg": 1} and oh["counts"] == {"matcher": 2, "rmg": 2}
    assert ow["entry"]["index"] == 7 and oh["entry"]["index"] == 7
    assert ow["weight"]["w"] == pytest.approx(K_DEP, rel=1e-12)
    assert "derived" in oh and "derived" not in ow
    for T in (612.5, 777.7, 999.0):
        wf, wr = te.override_weight(ow, T), te.override_weight(oh, T)
        assert wf * 1 / (wr * 2) == pytest.approx(te.kc_at(toy_form(False), T), rel=1e-12)
    rev = [k / c for k, c in zip(K_DEP, _kc())]
    assert [x * 2 for x in oh["weight"]["w"]] == pytest.approx(rev, rel=1e-9)
    assert len(table) == 2 and all(r["status"] == "ok" for r in table), table
    # the loaded matcher carries it at the requested temperature
    loaded = te.loaded_templates(out, 700.0)
    assert sorted(len(t.overrides) for t in loaded) == [1, 1]
    got = sorted(o.weight for t in loaded for o in t.overrides)
    want = sorted([te.override_weight(ow, 700.0), te.override_weight(oh, 700.0)])
    assert got == pytest.approx(want, rel=1e-12)


def test_overrides_name_only_the_matched_sizes():
    """Every override is the matched reaction itself (one unit longer, or any
    other size the templates reach, keeps the rule weight: no override)."""
    out, _table, _failures, h, w = _override_pair()
    for t in (h, w):
        for o in t["overrides"]:
            sides = [sorted(m.to_smiles() for m in te.molecules_of(te.rows_from_json(o["reaction"][k])))
                     for k in ("reactants", "products")]
            assert sorted(map(tuple, sides)) == sorted([tuple(sorted(C13_R)), tuple(C13_P)])


def test_an_override_needs_the_depository_selected_by_production():
    """The production path must select the depository entry the match names;
    a rule selection there (RMG would not use the depository) is a named
    gate failure, never an override at the rule rate."""
    entry, match = _c13_match()
    counter = DepositoryCounter([ToyDepository([entry])])
    payloads, counter = _toy_pair(counter=counter, certificate=_certificate(match))
    resolved = _res(payloads, counter)
    rule = _pair_rates({"reversed": K_RULE, "as_generated": K_RULE})
    out, _table, failures = te.depository_overrides(resolved, counter=counter, rates=rule,
                                                    kc=toy_kc)
    assert failures and all("depository" in f for f in failures), failures
    assert all(not t.get("overrides") for p in out for t in p["templates"])


def test_a_catalog_without_matches_is_untouched():
    """No certified match and no depository row: the payloads come back
    byte-identical (the PS catalogs)."""
    payloads, counter = _toy_pair()
    resolved = _res(payloads, counter)
    out, table, failures = te.depository_overrides(resolved, counter=counter, rates=None,
                                                   kc=toy_kc)
    assert out == resolved and table == [] and failures == []


# n-heptyl + n-octyl -> n-pentadecane: no certificate match and no locality
# row of the pair's templates (their rows reach C13, C14 and C16 only)
C14_R, C14_P = ["[CH2]CCCCCC", "[CH2]CCCCCCC"], ["CCCCCCCCCCCCCCC"]


def test_every_depository_reaction_the_loaded_templates_make_is_overridden():
    """A family with a certified match is rated over its whole depository:
    n-heptyl + n-octyl -> n-pentadecane is an entry no certificate or
    locality row named, yet the loaded templates make it, so RMG would rate
    it from the depository there too; it gets its override (the reverse
    derived)."""
    entry, match = _c13_match()
    other = toy_entry(8, C14_R, C14_P)
    counter = DepositoryCounter([ToyDepository([entry, other])])
    payloads, counter = _toy_pair(counter=counter, certificate=_certificate(match))
    resolved = _res(payloads, counter)
    c13, c14 = _depository_rates(counter, None), _depository_rates(
        counter, None, reactants=C14_R, products=C14_P)
    base = _pair_rates({"reversed": K_RULE, "as_generated": K_RULE})

    def rates(req):
        for f, smiles in ((c13, C13_R), (c14, C14_R)):
            if counter.species_of_rows(req.generation) == counter.species_of_rows(_rows(smiles)):
                return f(req)
        return base(req)
    out, _table, failures = te.depository_overrides(resolved, counter=counter, rates=rates,
                                                    kc=toy_kc)
    assert failures == [], failures
    for p in out:
        t = base_template(p)
        named = sorted(o["entry"]["index"] for o in t["overrides"])
        assert named == [7, 8], (p["orientation"], named)


def test_the_whole_depository_is_read_as_a_melt_holds_it():
    """The family's depository also holds a fused aromatic entry (a ring
    carbon with three benzene bonds, valence 4.5 as written) and a sulfur
    entry the loaded templates would make.  The melt holds rings in a Kekule
    form and never holds sulfur: the aromatic entry is read Kekule, the
    sulfur one can never occur, and neither stops the overrides."""
    entry, match = _c13_match()
    aromatic = toy_entry(9, ["c1ccc2ccccc2c1", "[CH3]"], ["Cc1cccc2ccccc12"])
    sulfur = toy_entry(10, ["CS[CH2]", "[CH3]"], ["CCSC"])
    counter = DepositoryCounter([ToyDepository([entry, aromatic, sulfur])])
    payloads, counter = _toy_pair(counter=counter, certificate=_certificate(match))
    resolved = _res(payloads, counter)
    out, _table, failures = te.depository_overrides(
        resolved, counter=counter, rates=_depository_rates(counter, None), kc=toy_kc)
    assert failures == [], failures
    for p in out:
        assert sorted(o["entry"]["index"] for o in base_template(p)["overrides"]) == [7]


def _production_pair(estimates):
    """ProductionRates on a stubbed pair: the forward estimate (its source the
    toy depository) and the compiler's reference-thermo reverse, which names
    no forward source (the compiler attaches ``forward_rate_source`` only to
    a source with a verbose comment)."""
    from types import SimpleNamespace
    fwd = {"kind": "RMG family estimate", "reference": "forward RMG family estimate",
           "source": "Kinetics Depository toy_cc/training", "entry": "C5H11 + C8H17 <=> C13H28",
           "rank": 6, "family_template_direction": "forward"}

    def record(orientation, rate_source, k):
        return SimpleNamespace(orientation=orientation, template="T", rate_source=rate_source,
                               k_table={"T": list(GRID), "k": list(k)}, ssa_multiplier=1.0,
                               event_id=orientation, status="emitted")
    records = [record("as_generated", dict(fwd), K_DEP) for _ in range(estimates)] + [
        record("reversed", {"kind": "reference-thermo reverse", "reference": "k_family / Kc"}, K_RULE)]
    rates = te.ProductionRates.__new__(te.ProductionRates)
    rates.calls = 0
    rates._pair = lambda request: (tuple(records), "", "generator")
    return rates, fwd


def test_the_production_reverse_names_the_forward_selection_it_inverts():
    """The reverse of a depository-rated pair is k_depository / Kc: the
    probe names the pair's forward estimate as its forward selection, so the
    gate sees the depository RMG rates it from."""
    rates, fwd = _production_pair(1)
    req = te.RateRequest("toy_cc", "T", "reversed", {}, {})
    got = rates(req)
    assert got.problems == ()
    assert got.provenance["selection"]["kind"] == "reference-thermo reverse"
    assert got.provenance["selection"]["forward"] == te._selection(fwd)
    ahead = rates(te.RateRequest("toy_cc", "T", "as_generated", {}, {}))
    assert ahead.provenance["selection"] == te._selection(fwd) and "forward" not in ahead.provenance["selection"]


def test_a_production_reverse_without_one_estimate_names_no_forward_selection():
    """No (or no single) forward estimate in the pair: no forward selection is
    invented, so the gate refuses the override by name."""
    for estimates in (0, 2):
        rates, _fwd = _production_pair(estimates)
        got = rates(te.RateRequest("toy_cc", "T", "reversed", {}, {}))
        assert "forward" not in got.provenance["selection"]
