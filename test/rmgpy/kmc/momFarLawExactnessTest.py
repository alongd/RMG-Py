"""Red tests (I-291): the MoM producer's far-law exactness test.

A far law (a move's A n + B, a pair transfer's c00 + c10 n + c01 m + c11 n m,
an interior scission's polynomial) must be what the producer's own per-length
rates say on the tail:

* it is checked against the exact rate at every checked length, with a
  tolerance relative to the rates at that length, which a large rate at
  another length of the term cannot widen;
* rates that are zero from some length on get the zero far branch, so the
  term's short part is a short-only channel (``record.short_lengths`` with
  ``far_law`` 0.0) and MoM carries no term for it;
* a far law negative anywhere on its tail is refused with an error naming the
  term;
* every far row and column rests on at least three checked lengths (one
  tested degree of freedom past the two that fix an affine law): the checked
  top is extended until it does, and a term still unresolved is refused and
  listed.

Every table here is synthetic (the generic toy linear polymer of the position
tests, or rates written as numbers): no polymer, term, template or pool name
enters a rule.
"""

import json

import pytest

import momPositionPoolsTest as T
from momPositionPoolsTest import table  # noqa: F401  (the toy table fixture)
from momProducerTest import BASE, LABEL, spec_of

from rmgpy.kmc import mom_producer as mp


def kernel_input(doc, kernels):
    return mp.kernel_input(doc, spec_of(BASE), kernels, temperature=700.0, phi=1.0, label=LABEL,
                           radical_label=T.RADICAL, alkene_label=T.ALKENE, interior_radical_label=T.MID,
                           interior_alkene_label=T.INTERNAL, volume_ratio=T.RATIO, far_law=True, positions=True)


def end_move(doc):
    """The toy table's end-radical move (one participant, anchored at the chain end, a positioned product):
    (kernel id, term)."""
    ((kid, term),) = [(k, t) for k, kernel in doc["kernels"].items() for t in kernel["templates"]
                      if len(kernel["classes"]) == 1 and (t.get("law") or {}).get("anchor") == "end"
                      and t["law"].get("positioned") and any(t["law"]["positioned"]["exact"].values())]
    return kid, term


def moves_of(terms, name):
    return [x for ts in terms["move"].values() for x in ts if name in x["template"].split("@")[0].split("+")]


def copy(doc):
    return json.loads(json.dumps(doc))


# --------------------------------------------------------------------------
# 1. per-length exactness: a large rate elsewhere in the term cannot mask a deviation
# --------------------------------------------------------------------------
def test_a_large_rate_at_another_length_cannot_mask_a_move_law_deviation():
    """Rates 1e10 at the shortest length, 1 at the next, 0 from there on: the line through the second and third
    lengths (7 - n) misses the exact zero rates by up to 8, inside 1e-9 of the term's largest rate (10) but not
    of the rate at those lengths.  The law is the zero law from the first zero."""
    r = {5: 1.0e10, 6: 1.0, **{n: 0.0 for n in range(7, 16)}}
    n0, c, _worst = mp._affine(r, 5, 15)
    assert (n0, c) == (7, [0.0, 0.0])


def test_a_large_rate_at_another_length_cannot_mask_a_pair_law_deviation():
    """The two-participant version: a row of 1e10 at the shortest length, a row of 1 at the next, zeros beyond."""
    F = {(n, m): 1.0e10 if n == 4 else 1.0 if n == 5 else 0.0 for n in range(4, 13) for m in range(4, 13)}
    n0, c, _worst = mp._bilinear(F, 4, 12)
    assert n0 == 6 and not any(c)


def test_a_polynomial_law_is_checked_per_length():
    """An interior-scission rate (degree 0) of 1e10 at the shortest length and 1 + 1e-6 n beyond: the constant
    through the top misses the longer lengths by up to 1e-5, inside 1e-9 of 1e10 but not of the rate there."""
    values = {4: 1.0e10, **{n: 1.0 + 1e-6 * n for n in range(5, 16)}}
    n0, _c = mp._polynomial_from(values, 0, 15)
    assert n0 is None


def test_a_constant_with_rounding_noise_stays_exact_per_length():
    """Per-length checking keeps the rounding noise of the table's sums (1e-15 relative) inside the law."""
    k = 3.096853570259764e-13
    F = {(n, m): k * (1.0 + 1e-15 * ((7 * n + 13 * m) % 5 - 2)) for n in range(4, 35) for m in range(4, 35)}
    n0, c, _worst = mp._bilinear(F, 4, 34)
    assert n0 == 4 and c[1:] == [0.0, 0.0, 0.0]
    n0, c = mp._polynomial_from({n: k * (1.0 + 1e-15 * (n % 3 - 1)) for n in range(4, 35)}, 0, 34)
    assert n0 == 4


# --------------------------------------------------------------------------
# 2. rates zero from some length on: a zero far branch, the short part a short-only channel
# --------------------------------------------------------------------------
def test_rates_zero_from_a_length_give_a_zero_far_branch_and_a_short_only_channel(table):
    """The end-radical move fires 1e10 times per chain at the shortest reachable length, once at the next, and
    never from there on (no far outcome).  Under a tolerance relative to the term's largest rate the producer
    fits the line through the last nonzero rate and the first zero, negative on its whole tail; per length it
    is the zero law from the first zero, so MoM carries no term and the record lists a short-only channel."""
    doc = copy(table)
    kid, term = end_move(doc)
    base = kernel_input(copy(table), [kid])["record"]["checks"]
    (first,) = {v["reachable_from_periods"] for k, v in base.items() if k.startswith(term["template"] + "@")}
    law = term["law"]
    product = next(rows for rows in law["positioned"]["exact"].values() if rows)[0][1][0]

    def at(n, mult):
        return [[mult, [[product[0], n, *product[2:]]]]]
    law["exact_to_n"] = first + 4
    law["positioned"]["exact"].update({str(first): at(first, 1.0e10), str(first + 1): at(first + 1, 1.0),
                                       **{str(n): [] for n in range(first + 2, first + 5)}})
    law["exact"].update({str(n): [] for n in range(first, first + 5)})
    law["far"]["outcomes"] = []
    law["positioned"]["far"]["outcomes"] = []
    terms = kernel_input(doc, [kid])
    assert not moves_of(terms, term["template"])
    (listed,) = [x for x in terms["record"]["short_lengths"] if x.get("far_law") == 0.0]
    assert listed["kernel"] == kid and listed["templates"] == [term["template"]]
    assert listed["law_from_periods"] == first + 2 and listed["last_nonzero_periods"] == first + 1


# --------------------------------------------------------------------------
# 3. a far law negative anywhere on its tail is refused, naming the term
# --------------------------------------------------------------------------
@pytest.mark.parametrize("c, n0, negative", [
    ([7.0, -1.0], 6, True),                   # 7 - n: negative from n = 8
    ([-1.0, 1.0], 0, True),                   # negative at its first length
    ([1.0, 0.0], 4, False),
    ([0.0, 2.0], 3, False),
    ([6.0, -1.0, 0.0, 0.0], 5, True),         # 6 - n
    ([100.0, 1.0, 1.0, -0.01], 4, True),      # positive on the checked lengths, negative at large n and m
    ([1.0, 0.5, 0.5, 0.1], 4, False),
    ([-37.0, 0.0, 0.0, 1.0], 6, True),        # -37 + n m: negative at (6, 6)
    ([-36.0, 0.0, 0.0, 1.0], 6, False),       # zero at (6, 6), positive beyond
    ([10.0, -1.0, 0.0, 1.0], 4, False),       # 10 + n (m - 1): its slope in n is m - 1 >= 3 on the tail
    ([10.0, -5.0, 0.0, 1.0], 4, True),        # its slope in n at m = 4 is -1
])
def test_a_far_law_is_negative_on_its_tail_where_its_least_value_there_is(c, n0, negative):
    why = mp._tail_negative(c, n0)
    assert (why is not None) == negative
    if negative:
        assert "negative" in why


def test_a_move_law_negative_on_its_tail_is_refused_naming_the_term(table, monkeypatch):
    """A fit that passed the per-length check on the checked lengths with a negative slope (the table's rates
    falling on every checked length): the producer refuses it with an error naming the kernel and template."""
    doc = copy(table)
    kid, term = end_move(doc)
    monkeypatch.setattr(mp, "_affine", lambda r, floor, top: (floor, [2.0 * top, -1.0], 0.0))
    with pytest.raises(mp.ProducerError) as err:
        kernel_input(doc, [kid])
    assert "negative" in str(err.value) and kid in str(err.value) and term["template"] in str(err.value)


# --------------------------------------------------------------------------
# 4. at least three checked lengths: the top is extended, an unresolved term refused and listed
# --------------------------------------------------------------------------
def test_a_move_law_on_two_checked_lengths_is_not_established():
    """Rates 1/n up to the second-last checked length, then on a line: the line rests on two lengths only."""
    r = {**{n: 1.0 / n for n in range(4, 14)}, 14: 5.0, 15: 6.0}
    assert mp._affine(r, 4, 15)[0] is None


def test_the_checked_top_extends_until_the_law_rests_on_three_lengths():
    """Rates 1/n below 14 and n from 14 on: at the checked top 15 no law rests on three lengths; one length more
    establishes n from 14, and the extended top is returned."""
    seen = []

    def values(t):
        seen.append(t)
        return {n: 1.0 / n if n < 14 else float(n) for n in range(4, t + 1)}
    n0, c, _worst, top = mp._fit_extending(mp._affine, values, 4, 15)
    assert (n0, c, top) == (14, [0.0, 1.0], 16) and seen == [15, 16]


def test_a_term_never_established_is_refused_at_the_extension_limit():
    """Rates 1/n on every length: no affine law on any range; the top extends to its limit and no law results."""
    n0, c, _worst, top = mp._fit_extending(mp._affine, lambda t: {n: 1.0 / n for n in range(4, t + 1)}, 4, 15)
    assert n0 is None and c is None and top == mp.CHECKED_TOP_LIMIT * 15


def test_a_term_unresolved_after_the_extension_is_refused_and_listed(table, monkeypatch):
    """A move with no law at any checked top: the producer extends the top to its limit, maps no term and lists
    the template with the top it was checked to."""
    doc = copy(table)
    kid, term = end_move(doc)
    tops = []

    def never(r, floor, top):
        tops.append(top)
        return None, None, 0.0
    monkeypatch.setattr(mp, "_affine", never)
    terms = kernel_input(doc, [kid])
    assert not moves_of(terms, term["template"])
    listed = [x for x in terms["record"]["unmapped"] if term["template"] in x["templates"]]
    assert listed and max(tops) == mp.CHECKED_TOP_LIMIT * min(tops)
    assert all(f"checked to {max(tops)} periods" in x["reason"] for x in listed)
