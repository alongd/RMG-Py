"""Red tests (I-252): the probe states a product law reads, on a declaration
whose period is not symmetric (``rmgpy.kmc.kernel_rates``).

The fixture is a reduced PP catalog (the payloads of three templates of the
PP catalog assembled at 29df6852c, with their reciprocal twins) and the PP deck
it is bound to.  Cutting a PP backbone bond leaves a CH end on the head-side
fragment and a CH2 end on the tail-side one, and interior backbone carbons
alternate between CH2 and CH; the interior positions of one probe carry two
classes.  (Until I-279 the law probes were head-side fragments only, so the
CH2 end was a site class no probe state carried; it is now probed at the
start of a piece, ``kmc-kernel-rates/4``.)
"""

import gzip
import json
import re
from pathlib import Path

import pytest

from rmgpy.kmc import kernel_rates as kr_
from rmgpy.kmc import run as kr

FIXTURES = Path(__file__).parent / "fixtures" / "i252"
T = 700.0
CH2_END_DONOR = "02007c0487b8:f029c5581706"        # H_Abstraction: radical_end:CH2 + alkene_terminal
INTERIOR_ACCEPTOR = "1190f86ca177:480e026ea139"    # Disproportionation: a radical_interior:CH participant


def pp_catalog():
    return json.loads(gzip.decompress((FIXTURES / "pp_catalog.json.gz").read_bytes()))


def pp_deck(tmp_path=None, chains=None):
    path = FIXTURES / "pp_deck.py"
    if chains is not None:
        text = path.read_text().replace("'dp_counts': {30: 100}", "'dp_counts': {30: %d}" % chains)
        assert text != path.read_text()
        path = tmp_path / f"pp_deck_{chains}.py"
        path.write_text(text)
    return kr.read_deck(path)


@pytest.fixture(scope="module")
def pp():
    states = []
    original = kr_._Laws._states

    def spy(self, label):
        out = original(self, label)
        states.append((self.probes, label, out))
        return out
    kr_._Laws._states = spy
    try:
        table = kr_.kernel_rates(pp_catalog(), pp_deck().spec, T)
    finally:
        kr_._Laws._states = original
    return {"table": table, "states": states}


def terms(table, template):
    return [t for k in table["kernels"].values() for t in k["templates"] if t["template"] == template]


def carried(probes, label, n, x):
    family, tau = kr_.family_of(label), kr_.type_of(label)
    state, paths = probes.melt()
    _seg, feature = probes.put(state, paths[0], family, n, x, tau=tau)
    return (family if family == "chain" else f"{family}:{kr_._label(state, feature)}") + (f"/{tau}" if tau else "")


def test_the_tail_side_end_class_is_probed_at_the_start_of_a_piece(pp):
    """The CH2 end has site counts and (since I-279) probes: the start of a piece of every length.  Its terms
    have a law or say why, and the reason is never a missing probe."""
    mine = [t for t in terms(pp["table"], CH2_END_DONOR)]
    assert mine
    for t in mine:
        assert t["law"] is not None or t.get("law_refused") or t.get("law_inert"), t
        assert not re.search(r"no probe of class", t.get("law_refused") or ""), t["law_refused"]
    states = next(s for _p, label, s in pp["states"] if label == "radical_end:CH2")
    assert states and all(x == "start" for _n, x in states)


def test_every_state_a_law_reads_carries_the_laws_class(pp):
    """Interior probes alternate CH2 and CH positions on PP: the law of a class
    reads only the positions that carry it."""
    checked = 0
    for probes, label, states in pp["states"]:
        for n, x in states:
            assert carried(probes, label, n, x) == label, (label, n, x)
            checked += 1
    assert checked and any(kr_.family_of(label) == "radical_interior" for _p, label, _s in pp["states"])


def test_the_interior_participant_has_a_law(pp):
    """One term per kernel (since I-279 the partner's classes include the other chain types' ends)."""
    mine = terms(pp["table"], INTERIOR_ACCEPTOR)
    laws = [t["law"][s] for t in mine if t["law"] for s in ("S0", "S1")]
    assert any(law["class"] == "radical_interior:CH" and law["far"]["outcomes"] for law in laws), mine


def test_every_term_has_a_law_or_says_why(pp):
    for k in pp["table"]["kernels"].values():
        for t in k["templates"]:
            assert t["law"] is not None or t.get("law_refused") or t.get("law_inert"), t


def test_the_table_does_not_read_the_decks_ensemble(tmp_path, pp):
    """The probes are built from the declaration alone: the deck's chain count
    changes nothing."""
    big = kr_.kernel_rates(pp_catalog(), pp_deck(tmp_path, 20000).spec, T)
    assert json.dumps(big, sort_keys=True) == json.dumps(pp["table"], sort_keys=True)


def test_progress_is_one_stderr_line_per_law_term(pp, capsys):
    table = kr_.kernel_rates(pp_catalog(), pp_deck().spec, T)
    err = capsys.readouterr().err
    n_terms = sum(len(k["templates"]) for k in table["kernels"].values())
    lines = [line for line in err.splitlines() if line.startswith("kernel_rates: law ")]
    assert len(lines) == n_terms
    assert all(re.search(r" \d+/%d " % n_terms, line) for line in lines)


TERTIARY_DONOR = "7313ebd9a7c7:54f630374459"      # H_Abstraction: a chain donor near an end -> tertiary radical


def test_the_far_law_check_does_not_read_which_interior_class_was_noted_first():
    """One length of an interior family carries two PP classes (CH2 and CH
    positions); its check must hold whichever class the length was labelled with
    first.  The donor here makes the tertiary radical at backbone atom 1 (since
    I-279 a product code carries its class: radical_interior:C, not the family)."""
    doc = json.loads(gzip.decompress((FIXTURES / "pp_tertiary_donor.json.gz").read_bytes()))
    original = kr_._Probes.label

    def neighbour(self, family, n, x=None, tau=None):
        """The label of an interior family of the declared type read at its second position."""
        xs = self.positions(family, n)
        if family not in kr_.INTERIOR_FAMILIES or len(xs) < 2 or x is not None or tau is not None:
            return original(self, family, n, x, tau)
        state, paths = self.melt()
        _seg, feature = self.put(state, paths[0], family, n, xs[1])
        return f"{family}:{kr_._label(state, feature)}"
    tables = [kr_.kernel_rates(doc, pp_deck().spec, T)]
    kr_._Probes.label = neighbour
    try:
        tables.append(kr_.kernel_rates(doc, pp_deck().spec, T))
    finally:
        kr_._Probes.label = original
    assert tables[0]["kernels"] == tables[1]["kernels"]
    # one term per acceptor class (since I-279 the ends of every chain type): the donor's product is the same in each
    donor = [law for t in terms(tables[0], TERTIARY_DONOR) if t["law"] for law in t["law"].values()
             if law["class"] == "chain"]
    assert donor and all(law["far"]["first"]["1"] == [[1, [["radical_interior:C", 1, 0, 0]]]] for law in donor)


def test_identifying_a_long_interior_radical_does_not_branch_over_hydrogens():
    """``_Probes.identify`` compares a molecule with the declared probe at every
    interior position; all but one are not isomorphic, with equal degree
    sequences.  VF2 on the explicit-hydrogen graphs branches over the
    interchangeable H of every CH2 and CH3 on each failed comparison (about 45 s
    each on 12 PP periods); with the hydrogen leaves folded it takes
    milliseconds."""
    import signal

    import numpy

    from rmgpy.kmc import melt_engine as me
    spec = pp_deck().spec
    engine = me.MeltEngine(kr_._single_chain_spec(spec, 4), me.PayloadCatalog([]), T,
                           numpy.random.default_rng(0))
    probes = kr_._Probes(spec, engine, 1, 2, 27)
    n = 12
    xs = probes.positions("radical_interior", n)
    for x in xs:
        probes.forms("radical_interior", n, x)
    state, paths = probes.melt()
    seg, _feature = probes.put(state, paths[0], "radical_interior", n, xs[-1])

    def too_slow(*_):
        raise TimeoutError("identify branched over the hydrogens")
    previous = signal.signal(signal.SIGALRM, too_slow)
    signal.alarm(20)
    try:
        label, length, positions = probes.identify(state, state.components.component_of(seg[0]))
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)
    assert (length, positions) == (n, [xs[-1]]) and label.startswith("radical_interior:")
