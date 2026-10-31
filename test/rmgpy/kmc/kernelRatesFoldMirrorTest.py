"""Red tests (I-288): a placed product folded to its nearer end on a chain type
that is its own reverse is named by the reading from that same end
(``rmgpy.kmc.kernel_rates._expand_placed``, ``_placed``).

A product code's class token is read from the end its distance is counted
from (an ``end_distance`` code's end, the head for ``offset`` and
``head_distance``, the tail for ``tail_distance``).  Expanded to another
length the position is folded to the nearer end; when the fold changes the
end, the token must be mirrored with it (``classify`` reads a measured product
the same way).  Before the fix the token was kept: on PP's c1~c1 and c2~c2 the
centre-bond alkene was named by its far-end reading at its near-end distance,
which is no molecule's name, at every length past the probes.

The structural check reads the molecule a placed product names off the
table's own pristine piece (``kernel_rates._piece``), independently of
``_expand_placed`` and ``classify``: its feature atoms' element and hydrogens
left, read from its nearer end (both readings at a tie).  The far law at the
centre of a long piece is also checked against the engine itself
(``kernelRatesStateModelTest.brute``), which the I-279b mirror test leaves out.
"""

from collections import Counter

import pytest

from kernelRatesMirrorTest import DONOR, donor_of, donor_terms, partner, pp  # noqa: F401 (pp is a fixture)
from kernelRatesStateModelTest import brute

from rmgpy.kmc import kernel_rates as kr_

PP_TYPES = {"c1~c1": {"start_offset": 0, "residue": 1, "symmetric": True, "declared": False},
            "c1~c2": {"start_offset": 0, "residue": 0, "symmetric": False, "declared": True},
            "c2~c2": {"start_offset": 1, "residue": 1, "symmetric": True, "declared": False}}
CH_C = "alkene_internal:CH=C/c1~c1"
C_CH = "alkene_internal:C=CH/c1~c1"


def expand(token, alpha, beta, gamma, where, L, k, types=PP_TYPES):
    (outcome,) = kr_._expand_placed([[1, [[token, alpha, beta, gamma, where]]]], L, k, 2, types)
    (prod,) = outcome
    return prod


# --------------------------------------------------------------------------
# the reader: a fold that changes the end mirrors the token
# --------------------------------------------------------------------------
def test_a_token_read_from_the_other_end_is_its_feature_reversed():
    assert kr_._mirror(CH_C) == C_CH
    assert kr_._mirror("alkene_internal:CH=C") == "alkene_internal:C=CH"
    assert kr_._mirror("radical_interior:CH/c2~c2") == "radical_interior:CH/c2~c2"


@pytest.mark.parametrize("n", [5, 12, 20])
def test_an_end_distance_folded_past_the_middle_is_read_from_the_nearer_end(n):
    """R6's unit case: a CH=C product coded at end distance n on 2n + 1 atoms (the centre bond read from the end
    it is counted from) is the C=CH reading at n - 1 from the other end."""
    assert expand(CH_C, 0, 0, 2 * n + 1, {"end_distance": n}, 2 * n + 1, 0) == (C_CH, n, "end_distance", n - 1)
    assert expand(CH_C, 0, 0, 2 * n + 1, {"end_distance": n - 1}, 2 * n + 1, 0) == (CH_C, n, "end_distance", n - 1)


@pytest.mark.parametrize("n", [5, 12, 20])
def test_an_offset_past_the_middle_of_a_whole_reactant_is_read_from_the_nearer_end(n):
    """A donor at the centre atom n of 2n + 1 atoms: the alkene one bond on (offset 0, from the head x = n) is the
    centre bond read from the tail; the one a bond back (offset -1, x = n - 1) keeps its head reading."""
    L = 2 * n + 1
    assert expand(CH_C, 1, 0, 0, {"offset": 0}, L, n) == (C_CH, n, "end_distance", n - 1)
    assert expand(CH_C, 1, 0, 0, {"offset": -1}, L, n) == (CH_C, n, "end_distance", n - 1)


def test_a_tie_keeps_the_forward_reading():
    """A feature exactly at the middle (lp - w - d == d, an own-reverse type of even residue) has two readings of
    one molecule; the forward one is kept, as classify keeps it."""
    types = {"c~c": {"start_offset": 0, "residue": 0, "symmetric": True, "declared": True}}
    token = "alkene_internal:CH=C"
    assert expand(token, 0, 0, 12, {"end_distance": 5}, 12, 0, types) == (token, 6, "end_distance", 5)
    assert expand(token, 1, 0, 0, {"offset": 0}, 12, 5, types) == (token, 6, "end_distance", 5)


def test_an_oriented_type_keeps_its_head_to_tail_token():
    token = "alkene_internal:CH=C"
    assert expand(token, 1, 0, 0, {"offset": 0}, 20, 15, PP_TYPES) == (token, 10, "tail_distance", 3)
    assert expand(token, 0, 0, 20, {"head_distance": 16}, 20, 0, PP_TYPES) == (token, 10, "tail_distance", 2)


def test_an_end_distance_on_an_oriented_type_is_refused_by_name():
    with pytest.raises(kr_.KernelRatesError, match="names no end"):
        expand("alkene_internal:CH=C", 0, 0, 20, {"end_distance": 3}, 20, 0, PP_TYPES)


# --------------------------------------------------------------------------
# the writer: a token is stored as read from the end its code counts from
# --------------------------------------------------------------------------
def fragment(index, lp, first, token):
    """A measured fragment of lp backbone atoms of c1~c1 holding its reactant's first (or last) atom, its alkene's
    head-side atom at fragment index ``index``, its token as classify reads it (from its nearer end)."""
    lo, hi = (0, lp - 1) if first else (5, 5 + lp - 1)
    place = kr_._Place(lo + index, lo, hi, first, not first, 2, min(index, lp - 2 - index), oriented=False)
    ends = ((first, not first, True),)
    return ("alkene_internal", token, (lp - 1) // 2, ends, place, lp, token)


@pytest.mark.parametrize("first", [True, False])
def test_a_cut_fragment_read_from_its_other_end_stores_the_cut_end_reading(first):
    """On 9 atoms: the feature 2 atoms from the reactant's end is 5 from the cut end.  classify reads it from the
    reactant's end; the code counts from the cut end, so it stores the cut-end reading, and expanding it back at
    the measured length gives classify's name."""
    lp = 9
    index = 2 if first else lp - 2 - 2
    prod = fragment(index, lp, first, CH_C)
    L, k = 30, (lp - 1 if first else 5)
    code = kr_._placed(prod, L, k, 2, False)
    assert code[0] == C_CH and code[-1] == ("end_distance", 5)
    where = {code[-1][0]: code[-1][1]}
    assert expand(code[0], *code[1:4], where, L, k) == (CH_C, 4, "end_distance", 2)


def test_a_cut_fragment_read_from_its_cut_end_is_stored_as_read():
    prod = fragment(5, 9, True, C_CH)       # 2 atoms from the cut end
    code = kr_._placed(prod, 30, 8, 2, False)
    assert code[0] == C_CH and code[-1] == ("end_distance", 2)


# --------------------------------------------------------------------------
# the table: every placed product is a molecule's canonical name
# --------------------------------------------------------------------------
def suffix(tau):
    return "" if tau == "c1~c2" else f"/{tau}"


_PIECES = {}


def piece(spec, t, n):
    """(backbone elements, hydrogens per backbone atom) of the pristine piece of chain type t and n periods."""
    key = (t["start_offset"], t["residue"], n)
    if key not in _PIECES:
        state, path = kr_._piece(spec, int(t["start_offset"]), int(t["residue"]), n, 1, 2)
        _PIECES[key] = ([state.atoms[a].element for a in path], [len(kr_._hydrogens(state, a)) for a in path])
    return _PIECES[key]


def names(spec, types, label, n, x):
    """The canonical names of the molecule with an interior feature of ``label``'s family and chain type at
    backbone atom x from the head of the pristine piece of n periods: its feature atoms' element and hydrogens
    left, read from its nearer end (both readings at a tie; none when x is past the middle of an own-reverse
    piece or the feature cannot sit there)."""
    family, tau = kr_.family_of(label), kr_.type_of(label) or "c1~c2"
    t, w = types[tau], kr_.WIDTH[family]
    elements, h = piece(spec, t, n)
    L = len(elements)
    if not 1 <= x <= L - w - 1 or (t["symmetric"] and x > L - w - x):
        return set()

    def read(idx):
        if any(h[j] < 1 for j in idx):
            return None
        return f"{family}:" + "=".join(elements[j] + ("H" if h[j] == 2 else f"H{h[j] - 1}" if h[j] > 2 else "")
                                       for j in idx) + suffix(tau)
    idx = list(range(x, x + w))
    got = {read(idx)}
    if t["symmetric"] and x == L - w - x:
        got.add(read(idx[::-1]))
    return got - {None}


def stated(law, near, n, x, types):
    """The positioned outcomes a participant law states on n periods (its feature at x when interior)."""
    try:
        if law["anchor"] == "interior":
            return kr_.predict_near(law, near, n, 2, x, types)
        if law["anchor"] == "end_distance":
            return kr_.predict_placed(law, n, 2, x, types)
        return kr_.predict_placed(law, n, 2, types=types) if law.get("positioned") else Counter()
    except kr_.KernelRatesError:
        return Counter()


def test_every_placed_mid_chain_product_is_a_canonical_name(pp):  # noqa: F811
    """Every law of the fixture table (with its near-end law), each interior participant at every backbone
    position its class sits on, on every length to twice the probes: each mid-chain product's label is a
    canonical name of the molecule its (n, kind, d) places on the pristine piece.  (A product placed at d = 0 is
    the producer's refusal, not a name.)"""
    table, spec, types = pp["table"], pp["spec"], pp["types"]
    checked, wrong = 0, Counter()
    for kern in table["kernels"].values():
        nslots = len(kern["classes"])
        for term in kern["templates"]:
            if term.get("law") is None:
                continue
            near = term.get("law_near_end") or {}
            slots = [(term["law"], near or None)] if nslots == 1 else \
                [(term["law"][f"S{p}"], near.get(f"S{p}")) for p in range(nslots)]
            for law, near_law in slots:
                tau = kr_.type_of(law["class"]) or "c1~c2"
                for n in range(1, 2 * pp["n_hi"] + 3):
                    L = 2 * n + int(types[tau]["residue"])
                    xs = [x for x in range(1, L) if law["class"] in names(spec, types, law["class"], n, x)] \
                        if law["family"] in kr_.INTERIOR_FAMILIES else [None]
                    for x in xs:
                        for outcome in stated(law, near_law, n, x, types):
                            for label, pn, kind, d in (p for p in outcome if p[0] != "join" and len(p) == 4):
                                if d < 1:
                                    continue
                                Lp = 2 * pn + int(types[kr_.type_of(label) or "c1~c2"]["residue"])
                                px = Lp - kr_.WIDTH[kr_.family_of(label)] - d if kind == "tail_distance" else d
                                checked += 1
                                if label not in names(spec, types, label, pn, px):
                                    wrong[(law["class"], label, kind, d, pn)] += 1
    assert checked
    assert not wrong, sorted(wrong)[:10]


@pytest.mark.parametrize("tau", ["c1~c1", "c2~c2"])
def test_the_donor_far_law_at_the_centre_equals_the_engine(pp, tau):  # noqa: F811
    """The donor on lengths past the probes at the two positions nearest the centre of its own-reverse piece (the
    ones the I-279b test leaves out): the far law's positioned products are the engine's."""
    label = f"{DONOR}/{tau}"
    law = donor_terms(pp["table"])[label]["law"]["S1"]
    ci, p = donor_of(pp, label)
    probes = pp["probes"]
    checked = 0
    for n in (pp["n_hi"] + 1, pp["n_hi"] + 2):
        L = 2 * n + 1
        for x in probes.positions("radical_interior", n, tau=tau):
            if not probes.carries(label, n, x) or L - 1 - 2 * x >= 4:
                continue
            want = brute(pp, ci, p, "radical_interior", n, x, tau=tau, partner=partner(pp))
            got = Counter({tuple(sorted(k, key=repr)): v
                           for k, v in kr_.predict_near(law, None, n, 2, x, types=pp["types"]).items()})
            assert got == want, (label, n, x)
            checked += 1
    assert checked
