"""The kernel table counts a transition two deduplicated templates share once (I-271).

A template that defers to another (``defers_to``) makes only the transitions
that template does not.  The kernel table reads the melt engine's index: a
unimolecular class's transitions, and a bimolecular class's sites per
participant, whose kernel term is a product of site counts.  The PP overlaps
the catalog deduplicates have no kernel term on the declared states (their
templates are kMC-only there), so the case is a PE toy template and a twin of
it deferring to it: every transition of the twin is the template's, and the
twin adds nothing to the table.  A deferring term's laws read its own sites
too: the end-distance and near-end laws of an interior participant (a
mid-chain radical) as much as its far law (I-275).
"""

import copy
import hashlib
import re

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


def _mid_abstraction(name, narrow=None):
    """A mid-chain radical (an interior participant) taking a mid CH2's H:
    the radical between two CH2 ports, or with ``narrow`` "far" a CH2 and a
    CH2 port on each side (no radical next to the CH2 beside a chain end),
    "near" a CH2 then the chain's CH3 on one side (that radical only)."""
    acceptor, bonds = [WindowAtom("C", 0, radical=1), WindowAtom("H", 0)], [(0, 1, 1)]
    if narrow is None:
        acceptor += [WindowAtom("C", 0, "port", S.MID), WindowAtom("C", 0, "port", S.MID)]
        bonds += [(0, 2, 1), (0, 3, 1)]
    elif narrow == "far":
        for o in (2, 6):
            acceptor += [WindowAtom("C", 0), WindowAtom("H", 0), WindowAtom("H", 0),
                         WindowAtom("C", 0, "port", S.MID)]
            bonds += [(0, o, 1), (o, o + 1, 1), (o, o + 2, 1), (o, o + 3, 1)]
    else:
        acceptor += [WindowAtom("C", 0), WindowAtom("H", 0), WindowAtom("H", 0), WindowAtom("C", 0),
                     WindowAtom("H", 0), WindowAtom("H", 0), WindowAtom("H", 0), WindowAtom("C", 0, "port", S.MID)]
        bonds += [(0, 2, 1), (2, 3, 1), (2, 4, 1), (2, 5, 1), (5, 6, 1), (5, 7, 1), (5, 8, 1), (0, 9, 1)]
    o = len(acceptor)
    donor = [WindowAtom("C", 1), WindowAtom("H", 1), WindowAtom("H", 1),
             WindowAtom("C", 1, "port", S.MID), WindowAtom("C", 1, "port", S.MID)]
    bonds += [(o, o + 1, 1), (o, o + 2, 1), (o, o + 3, 1), (o, o + 4, 1)]
    return Template(name, "toy_mid_abstraction", acceptor + donor, bonds,
                    [BondEdit(o, o + 1, 0), BondEdit(0, o + 1, 1), AtomEdit(0, radical=0), AtomEdit(o, radical=1)],
                    weight=S.PE_WEIGHTS["h_abstraction"])


def _deferring_term(narrow):
    """The term of a mid-chain H abstraction deferring to its ``narrow``
    version, and the slot of its mid-chain radical."""
    payloads = toy_payloads([t for t in S.pe_templates() if t.name in PE_LIVE]
                            + [_mid_abstraction("mid_abstraction"), _mid_abstraction("narrow_mid", narrow)])
    (mine,) = [p for p in payloads if p["record_event_id"] == "evt_mid_abstraction"]
    (target,) = [p for p in payloads if p["record_event_id"] == "evt_narrow_mid"]
    mine["templates"][0]["defers_to"] = [{"record": target["record_event_id"], "key": target["templates"][0]["key"]}]
    name = mine["record_event_id"][4:16] + ":" + mine["templates"][0]["key"][:12]
    (term,) = [x for k in _table(payloads)["kernels"].values() for x in k["templates"] if x["template"] == name]
    (slot,) = [s for s in ("S0", "S1") if "per_feature" in term[s]]
    return term, slot


def _near_end_sites(count):
    """{end distance: sites} of an interior count on its longest chain where they are not its far count."""
    keys = [re.fullmatch(r"chain_backbone_atoms=(\d+),end_distance=(\d+)", k) for k in count["profile"]]
    longest = max(int(m[1]) for m in keys if m)
    return {m[2]: _site_count(count, m[0])[0] for m in keys if m and int(m[1]) == longest}


def test_an_interior_participant_left_only_near_an_end_has_an_end_distance_law_of_its_own_sites():
    """The far-narrow version makes every mid-chain radical's transitions but
    the one beside the CH2 next to an end: the deferring term's radical sites
    are there only, and its law is the end-distance law of those sites, not a
    far law of every site."""
    term, slot = _deferring_term("far")
    assert term[slot]["per_feature"] == 0 and term["law"] is None and term["sites_deferred_to"]
    law = term["law_end_distance"][slot]
    assert law["anchor"] == "end_distance"
    assert {d: at["sites"] for d, at in law["at"].items()} \
        == {d: n for d, n in _near_end_sites(term[slot]).items() if n} == {"2": 1}


def test_a_near_end_law_leaves_out_the_sites_a_term_defers():
    """The near-narrow version makes the mid-chain radical's transitions
    beside the CH2 next to an end: the deferring term keeps its far sites,
    and its near-end law lists that end distance with none."""
    term, slot = _deferring_term("near")
    assert term[slot]["per_feature"] == 1 and term["sites_deferred_to"]
    own = _near_end_sites(term[slot])
    assert own == {"1": 0, "2": 0}
    assert {d: at["sites"] for d, at in term["law_near_end"][slot]["at"].items()} == own


def test_pairs_a_template_makes_itself_that_are_no_product_leave_it_kmc_only():
    """With the narrower donor too, the toy's own transitions are a propyl-end
    acceptor with any donor and any acceptor with a donor the narrower template
    misses: no single product of site counts."""
    _payloads, deferring, mine, _target = _deferring(long_donor=True)
    table = _table(deferring)
    assert {"template": mine, "family": "toy_h_abstraction", "reason": "deferred_not_a_product"} \
        in table["kmc_only"]
    assert mine not in {x["template"] for k in table["kernels"].values() for x in k["templates"]}
