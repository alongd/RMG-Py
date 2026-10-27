"""The kernel table counts a transition two deduplicated templates share once (I-271).

A template that defers to another (``defers_to``) makes only the transitions
that template does not.  The kernel table reads the melt engine's index: a
unimolecular class's transitions, and a bimolecular class's sites per
participant, whose kernel term is a product of site counts.  The PP overlaps
the catalog deduplicates have no kernel term on the declared states (their
templates are kMC-only there), so the case is a PE toy template and a twin of
it deferring to it: every transition of the twin is the template's, and the
twin adds nothing to the table.
"""

import copy
import hashlib

import pytest

import melt_engine_support as S
from meltRunTest import PE_LIVE, toy_payloads

from rmgpy.kmc import kernel_rates as kr_
from rmgpy.kmc.window_matcher import AtomEdit, BondEdit, Template, WindowAtom
from rmgpy.polymer_spec import LinearHomopolymerSpec

T = 700.0


def _table(payloads):
    spec = LinearHomopolymerSpec.from_mapping(S.pe_mapping({10: 1}))
    return kr_.kernel_rates({"schema": "kmc-template-catalog/1", "payloads": payloads}, spec, T,
                            allow_unbound=True)


def _with_twin(payloads, name):
    """``payloads`` and a twin of template ``name`` deferring to it."""
    (original,) = [p for p in payloads if p["record_event_id"] == f"evt_{name}"]
    twin = copy.deepcopy(original)
    twin["record_event_id"] = f"evt_twin_{name}"
    (t,) = twin["templates"]
    t["key"] = hashlib.sha256(f"twin:{name}".encode()).hexdigest()
    t["defers_to"] = [{"record": original["record_event_id"], "key": original["templates"][0]["key"]}]
    return payloads + [twin], twin["record_event_id"][4:16] + ":" + t["key"][:12]


@pytest.mark.parametrize("name", ["h_abstraction", "homolysis"])
def test_a_template_deferring_every_transition_adds_nothing_to_the_kernel_table(name):
    payloads = toy_payloads([t for t in S.pe_templates() if t.name in PE_LIVE])
    base = _table(payloads)
    both, twin = _with_twin(payloads, name)
    table = _table(both)
    assert twin in [r["template"] for r in table["templates"]]
    # the twin is in no kernel: each of its transitions is its original's
    placed = {x["template"] for k in table["kernels"].values() for x in k["templates"]}
    assert twin not in placed
    reasons = {o["template"]: o["reason"] for o in table["kmc_only"]}
    assert reasons[twin] == {"h_abstraction": "deferred", "homolysis": "no_site_on_declared_states"}[name]
    # and every other template's row and term are the table without it
    assert [r for r in table["templates"] if r["template"] != twin] == base["templates"]
    assert {kid: [x for x in k["templates"] if x["template"] != twin] for kid, k in table["kernels"].items()} \
        == {kid: k["templates"] for kid, k in base["kernels"].items()}


def _narrow_h_abstraction(name, long_donor):
    """The toy H abstraction on a wider window: the acceptor's CH2 neighbour
    bonded to a further CH2 (no propyl end), and with ``long_donor`` the
    donor's too (no donor next to a chain's last CH2).  Its sites are the
    toy's at the same atoms in the same roles, on fewer states."""
    acceptor = [WindowAtom("C", 0, radical=1), WindowAtom("H", 0), WindowAtom("H", 0), WindowAtom("C", 0),
                WindowAtom("H", 0), WindowAtom("H", 0), WindowAtom("C", 0, "port", S.MID)]
    donor = [WindowAtom("C", 1), WindowAtom("H", 1), WindowAtom("H", 1), WindowAtom("C", 1, "port", S.MID)]
    bonds = [(0, 1, 1), (0, 2, 1), (0, 3, 1), (3, 4, 1), (3, 5, 1), (3, 6, 1), (7, 8, 1), (7, 9, 1), (7, 10, 1)]
    if long_donor:
        donor += [WindowAtom("C", 1), WindowAtom("H", 1), WindowAtom("H", 1), WindowAtom("C", 1, "port", S.MID)]
        bonds += [(7, 11, 1), (11, 12, 1), (11, 13, 1), (11, 14, 1)]
    else:
        donor += [WindowAtom("C", 1, "port", S.MID)]
        bonds += [(7, 11, 1)]
    return Template(name, "toy_h_abstraction", acceptor + donor, bonds,
                    [BondEdit(7, 8, 0), BondEdit(0, 8, 1), AtomEdit(0, radical=0), AtomEdit(7, radical=1)],
                    weight=S.PE_WEIGHTS["h_abstraction"])


def _deferring(long_donor):
    """The toy catalog with a narrower H abstraction, and the same with the toy's
    H abstraction deferring to it; the two loaded names."""
    payloads = toy_payloads([t for t in S.pe_templates() if t.name in PE_LIVE]
                            + [_narrow_h_abstraction("narrow_abstraction", long_donor)])
    deferring = copy.deepcopy(payloads)
    (mine,) = [p for p in deferring if p["record_event_id"] == "evt_h_abstraction"]
    (target,) = [p for p in deferring if p["record_event_id"] == "evt_narrow_abstraction"]
    mine["templates"][0]["defers_to"] = [{"record": target["record_event_id"],
                                          "key": target["templates"][0]["key"]}]
    name = lambda p: p["record_event_id"][4:16] + ":" + p["templates"][0]["key"][:12]
    return payloads, deferring, name(mine), name(target)


def _site_count(count, position):
    return count["profile"].get(position, [count["per_feature"]]) if "profile" in count else count


def test_a_term_counts_only_the_pairs_its_template_makes_itself():
    """The narrower template has every donor site of the toy's and its acceptor
    sites but the propyl end's: the toy's own transitions are its propyl-end
    acceptors with every donor, one product of site counts."""
    payloads, deferring, mine, target = _deferring(long_donor=False)
    base, table = _table(payloads), _table(deferring)
    assert table["templates"] == base["templates"]       # a row counts all its template's sites
    terms = lambda tab, name: {kid: x for kid, k in tab["kernels"].items() for x in k["templates"]
                               if x["template"] == name}
    assert {kid: x for kid, k in table["kernels"].items() for x in k["templates"] if x["template"] != mine} \
        == {kid: x for kid, k in base["kernels"].items() for x in k["templates"] if x["template"] != mine}
    full, own, narrow = terms(base, mine), terms(table, mine), terms(table, target)
    assert set(own) == set(full) == set(narrow) and own
    restricted = 0
    for kid, term in own.items():
        assert term["sites_deferred_to"] == [target]
        assert term["w_over_divisor"] == full[kid]["w_over_divisor"]
        for slot in ("S0", "S1"):
            if "per_feature" not in term[slot]:          # the donor's chain: every site is its own
                assert term[slot] == full[kid][slot]
                continue
            # the acceptor: its sites less the narrower template's, at every position
            positions = set(term[slot]["profile"]) | set(full[kid][slot]["profile"]) | \
                set(narrow[kid][slot]["profile"])
            for position in positions | {None}:
                (mine_n,), (all_n,), (its_n,) = (_site_count(c, position) for c in
                                                 (term[slot], full[kid][slot], narrow[kid][slot]))
                assert mine_n + its_n == all_n, position
            restricted += term[slot] != full[kid][slot]
        assert (term["law"] is None) == (full[kid]["law"] is None)
    assert restricted


def test_pairs_a_template_makes_itself_that_are_no_product_leave_it_kmc_only():
    """With the narrower donor too, the toy's own transitions are a propyl-end
    acceptor with any donor and any acceptor with a donor the narrower template
    misses: no single product of site counts."""
    _payloads, deferring, mine, _target = _deferring(long_donor=True)
    table = _table(deferring)
    assert {"template": mine, "family": "toy_h_abstraction", "reason": "deferred_not_a_product"} \
        in table["kmc_only"]
    assert mine not in {x["template"] for k in table["kernels"].values() for x in k["templates"]}
