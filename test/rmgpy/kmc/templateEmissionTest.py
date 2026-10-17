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
    owner, other = (base_template(p) for p in merged)
    assert owner["status"] == "emitted"
    assert owner["collision"]["records"] == ["evt_" + "1" * 64, "evt_" + "2" * 64]
    assert other["status"] == "merged" and other["merged_into"]["record"] == "evt_" + "1" * 64
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
    for p in out:
        t = base_template(p)
        assert t["status"] == "refused" and t["refusal"]["status"] == "refused:collision"
        assert t["refusal"]["keys"] == [t["key"]] and len(t["refusal"]["weights"]) == 2
        assert t["refusal"]["records"] == ["evt_" + "1" * 64, "evt_" + "2" * 64]
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
    production rate path checks the rate-rule preparation policy)."""
    import template_emission_acceptance as acc
    if not acc.DATABASE.is_dir():
        pytest.skip(f"needs the RMG database at {acc.DATABASE}")
    return acc.load_database(["H_Abstraction", "Disproportionation"])


def _rates(compile_db, family):
    import template_emission_acceptance as acc
    return te.ProductionRates(compile_db, [family], database_path=acc.DATABASE,
                              database_sha="i138-test")


def _request(family, case, witness, orientation="as_generated"):
    r, p = te.template_sides(witness)
    return te.RateRequest(family, ";".join(case.pathway), orientation, r, p, witness,
                          "emitted", ())


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
    """The same three copies 0.004 apart merge into the least record id's
    template at its weight; every input order gives one merged catalog."""
    members = _group([0.0, 0.004, 0.008])
    canon, owner_w = set(), None
    for order in itertools.permutations(members):
        out, table = te.resolve_collisions(list(order))
        canon.add(_canon(out) + te.canonical_json(table))
        assert [r["disposition"] for r in table] == ["merged"] and table[0]["copies"] == 3
        assert table[0]["log10_diameter"] == pytest.approx(0.008, abs=1e-9)
        by = {p["record_event_id"]: base_template(p) for p in out}
        assert by[_rid(1)]["status"] == "emitted"
        assert by[_rid(1)]["collision"]["records"] == [_rid(1), _rid(2), _rid(3)]
        assert by[_rid(1)]["weight"] == base_template(members[0])["weight"]
        assert by[_rid(2)]["status"] == by[_rid(3)]["status"] == "merged"
        assert by[_rid(2)]["merged_into"]["record"] == _rid(1)
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
KC = tuple(2.0e-3 * math.exp(t / 400.0) for t in GRID)


def toy_kc(reactant_rows, product_rows, temperatures):
    """The toy thermo: Kc of a recombination of two species is ``KC``, of the
    homolysis it undoes ``1 / KC`` (concentration units)."""
    assert list(temperatures) == list(GRID)
    one = len(rm.split_rows(dict(reactant_rows))) == 1
    return [1.0 / k if one else k for k in KC]


def _pair_rates(by_orientation):
    def rates(req):
        return te.RateResult(GRID, by_orientation[req.orientation], 1.0, {}, ())
    return rates


def _toy_pair(hid=None, wid=None, scale=1.0, h_scale=1.0, counter=None):
    """Homolysis H of C12 at (3,4) (reversed record) and the recombination W
    of its two products (as_generated), linked by ``reverse_of``;
    ``k_W = Kc * k_H`` times ``scale`` both ways, ``h_scale`` on H alone."""
    h = homolysis(12, (3, 4))
    r = recombination(4, 8)
    k_h = tuple(k * scale * h_scale for k in K_RULE)
    k_w = tuple(k * scale * c for k, c in zip(K_RULE, KC))
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
                         counter=counter, rates=rates, check_upstream=False, kc=toy_kc)
    case_w, row_w, ctx_w = member(r, direction="forward")
    p_w = te.emit_record(rec_w, family=ToyFamily, members=[(case_w, row_w, ctx_w)],
                         counter=counter, rates=rates, check_upstream=False, kc=toy_kc)
    return [p_h, p_w], counter


def _check(payloads, counter, **kw):
    return te.check_cycles(payloads, counter=counter, kc=toy_kc, **kw)


def test_toy_reciprocal_pair_satisfies_the_cycle_against_thermo():
    """AMENDED (it encoded the circular ratio): the pair passes against the
    independent thermo Kc, which the table reports."""
    payloads, counter = _toy_pair()
    assert [p["status"] for p in payloads] == ["emitted", "emitted"]
    out, table = _check(payloads, counter)
    assert len(table) == 1 and table[0]["status"] == "ok", table
    assert table[0]["max_ratio_error"] <= 1e-12 and table[0]["rows"] >= 4
    assert table[0]["worst"]["Kc"] in (pytest.approx(1 / KC[GRID.index(table[0]["worst"]["T"])]),
                                       pytest.approx(KC[GRID.index(table[0]["worst"]["T"])]))
    assert [p["status"] for p in out] == ["emitted", "emitted"]
    assert te.one_way_problems(out) == [] and te.cycle_gates(table, out) == []


def test_homolysis_rate_and_weight_doubled_is_detected_against_thermo():
    """k_H and its weights doubled together keep the records' own k_F / k_R
    self-consistent (and the calibration exact); the independent Kc is
    unchanged, so the loaded pair fails."""
    payloads, counter = _toy_pair(h_scale=2.0)
    assert all(c["ok"] for c in payloads[0]["calibration"])
    out, table = _check(payloads, counter)
    assert table[0]["status"] == "refused:cycle"
    assert table[0]["max_ratio_error"] == pytest.approx(1.0, rel=1e-9)
    assert "Kc" in table[0]["reason"] and te.cycle_gates(table, out)


def test_perturbing_one_weight_by_1e_6_fails_the_cycle_and_is_gated():
    """AMENDED: perturb the loaded recombination weight at one temperature.
    Both records are refused, the failure is an acceptance gate, and nothing
    of the pair stays loaded (no one-way channel)."""
    payloads, counter = _toy_pair()
    base_template(payloads[1])["weight"]["w"][10] *= 1 + 1e-6
    payloads[1] = te.finish(payloads[1])
    out, table = _check(payloads, counter)
    assert table[0]["status"] == "refused:cycle"
    assert table[0]["max_ratio_error"] == pytest.approx(1e-6, rel=1e-3)
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
    payloads[1] = te.finish(payloads[1])
    out, table = _check(payloads, counter)
    assert table[0]["status"] == "refused:cycle" and table[0]["bijection_failures"] >= 1
    assert "inverse edit" in table[0]["reason"]


def test_a_cycle_without_a_channel_transition_is_refused():
    """AMENDED (signature): the forward record's template cannot make the
    reverse record's cleavage on its witness: named refusal of both."""
    payloads, counter = _toy_pair()
    other = emit(homolysis(12, (5, 6)), radius=4)
    payloads[0] = te.finish(dict(other, record_event_id=payloads[0]["record_event_id"],
                                 reverse_of=payloads[0]["reverse_of"]))
    out, table = _check(payloads, counter)
    assert table[0]["status"] == "refused:cycle"
    for p in out:
        assert base_template(p)["refusal"]["status"] == "refused:cycle"


def test_a_pair_with_an_emitted_witness_weight_is_checked():
    """Templates reaching beyond the proxy graph weigh on their emitted
    witness (``counts.record`` absent): checked there, not skipped."""
    payloads, counter = _toy_pair()
    for p in payloads:
        p["weight_source"], p["counts"] = "emitted_witness", {"record": None}
    out, table = _check([te.finish(p) for p in payloads], counter)
    assert table[0]["status"] == "ok" and table[0]["rows"] >= 4
    base_template(payloads[0])["weight"]["w"][3] *= 1.001
    out, table = _check([te.finish(p) for p in payloads], counter)
    assert table[0]["status"] == "refused:cycle"


def test_a_pair_of_context_templates_is_checked():
    """No base template on either side (only context templates load): the
    pair is still checked, and a bad context weight is caught."""
    payloads, counter = _toy_pair()
    for p in payloads:
        base_template(p)["base"] = False
    out, table = _check([te.finish(p) for p in payloads], counter)
    assert table[0]["status"] == "ok" and table[0]["rows"] >= 4
    base_template(payloads[0])["weight"]["w"][5] *= 1.01
    out, table = _check([te.finish(p) for p in payloads], counter)
    assert table[0]["status"] == "refused:cycle" and te.cycle_gates(table, out)


def test_a_live_pair_without_its_reverse_record_is_refused_when_complete():
    payloads, counter = _toy_pair()
    out, table = _check(payloads[:1], counter)
    assert table[0]["status"] == "refused:cycle" and "not in the catalog" in table[0]["reason"]
    assert te.loaded_templates(out, 700.0) == []
    out, table = _check(payloads[:1], counter, complete=False)
    assert table[0]["status"] == "not_evaluable" and base_template(out[0])["status"] == "emitted"


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
    assert sorted(r["disposition"] for r in table) == ["merged", "merged"]
    live = {p["record_event_id"]: base_template(p)["status"] for p in out}
    assert live == {_rid(1): "emitted", _rid(6): "emitted", _rid(2): "merged", _rid(5): "merged"}
    assert all(len(r["pair_deviations"]) == 2 for r in table)
    out2, cycles = _check(out, counter)
    assert all(c["status"] == "ok" for c in cycles), cycles
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
    for k, sh in ((1, 0.008), (2, -0.0)):
        pass
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
    kc = te.ThermoKc(GasPhaseRMGReferenceThermo(db.thermo, "test-only-database-commit"))
    ethane = rm.molecule_rows(Molecule().from_smiles("CC"))
    # the methyl pair: two components of 4 atoms, ids disjoint
    one = rm.molecule_rows(Molecule().from_smiles("[CH3]"))
    two = dict(one)
    two.update({a + 100: (x[0], x[1], x[2], x[3], x[4], tuple((n + 100, o) for n, o in x[5]))
                for a, x in one.items()})
    got = kc(ethane, two, [700.0, 800.0])
    sp = lambda s: Species(molecule=[Molecule().from_smiles(s)])
    rxn = Reaction(reactants=[sp("CC")], products=[sp("[CH3]"), sp("[CH3]")])
    rxn_kc = [rxn.get_equilibrium_constant(t, type="Kc") for t in (700.0, 800.0)]
    assert got == pytest.approx(rxn_kc, rel=1e-9)
    assert te.ssa_of_rows(two, ToyCounterReal(rmg)) == 2.0
    assert te.ssa_of_rows(ethane, ToyCounterReal(rmg)) == 1.0


def ToyCounterReal(rmg):
    return rmg[2]


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
