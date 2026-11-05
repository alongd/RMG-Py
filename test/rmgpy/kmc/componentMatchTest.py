"""I-303: the exact component match of certification (``rewrite_problems``)
is decided in work bounded by the molecule, not by its symmetry.

The budget counts calls of ``state._order``, through which every bond-order
comparison of the pre-I-303 search goes, and aborts the search once it is
spent: a search that would run for hours fails here at once.  The fixed
search reads bond orders once per edge before it searches, so it is held to
a CPU-time bound as well (checked after the call, never by interruption).
``_reference_match`` is the search as it stood before I-303 (85e78cd70), the
oracle of the verdicts on cases small enough for it to finish."""
import random
import time

import pytest

from rmgpy.kmc import certificates as cm
from rmgpy.kmc import state
from rmgpy.molecule.molecule import Molecule

# polystyrene trimer, saturated ends: three phenyl rings on a chain
PS_TRIMER = "CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)C"
BUDGET = 200_000  # bond-order comparisons
CPU_SECONDS = 2.0


class _Spent(Exception):
    pass


def _graph(smiles):
    return state._parse_adjacency(Molecule().from_smiles(smiles).to_adjacency_list())


def _relabel(graph, perm):
    out = {}
    for i, node in graph.items():
        node = dict(node)
        node["edges"] = {perm[j]: o for j, o in node["edges"].items()}
        out[perm[i]] = node
    return out


def _permuted(graph, seed):
    keys = sorted(graph)
    shuffled = keys[:]
    random.Random(seed).shuffle(shuffled)
    perm = dict(zip(keys, shuffled))
    return _relabel(graph, perm), perm


def _double_edge_swap(graph, ab, cd):
    """(a, b), (c, d) -> (a, d), (c, b): every degree, signature and bond
    order is kept; the structure is not."""
    (a, b), (c, d) = ab, cd
    out = {i: dict(node, edges=dict(node["edges"])) for i, node in graph.items()}
    o1, o2 = out[a]["edges"].pop(b), out[c]["edges"].pop(d)
    del out[b]["edges"][a], out[d]["edges"][c]
    out[a]["edges"][d] = out[d]["edges"][a] = o1
    out[c]["edges"][b] = out[b]["edges"][c] = o2
    return out


def _ring_bonds(graph):
    """Two aromatic CH-CH bonds on different phenyl rings."""
    def ch(i):
        return graph[i]["element"] == "C" and sum(graph[j]["element"] == "H"
                                                  for j in graph[i]["edges"]) == 1
    bonds = sorted((i, j) for i in graph for j in graph[i]["edges"]
                   if i < j and ch(i) and ch(j))
    ring, grown = set(bonds[0]), True
    while grown:
        grown = False
        for i, j in bonds:
            if (i in ring) != (j in ring):
                ring |= {i, j}
                grown = True
    return bonds[0], next((i, j) for i, j in bonds if i not in ring and j not in ring)


def _within(budget, fn, *args):
    """fn(*args), raising _Spent past ``budget`` bond-order comparisons and
    failing if it took more than CPU_SECONDS."""
    calls = [0]
    original = state._order

    def counted(value):
        calls[0] += 1
        if calls[0] > budget:
            raise _Spent(f"more than {budget} bond-order comparisons")
        return original(value)

    state._order = counted
    start = time.process_time()
    try:
        out = fn(*args)
    finally:
        state._order = original
    assert time.process_time() - start < CPU_SECONDS
    return out


def _reference_match(graph, members, product, pins) -> bool:
    """The pre-I-303 search, verbatim (the oracle)."""
    _order = state._order
    if len(members) != len(product):
        return False
    pinned = set(pins.values())
    candidates = {}
    for src in members:
        node = graph[src]
        options = [pins[src]] if src in pins else [t for t in product if t not in pinned]
        candidates[src] = [t for t in options if t in product
                           and cm._runtime_signature(product[t]) == cm._runtime_signature(node)
                           and len(product[t]["edges"]) == len(node["edges"])]
        if not candidates[src]:
            return False
    order = sorted(members, key=lambda a: (len(candidates[a]), a))

    def visit(position, matched, used):
        if position == len(order):
            return True
        src = order[position]
        for target in candidates[src]:
            if target in used:
                continue
            if all(_order(product[target]["edges"].get(matched[n])) == _order(o)
                   for n, o in graph[src]["edges"].items() if n in matched):
                matched[src] = target
                used.add(target)
                if visit(position + 1, matched, used):
                    return True
                del matched[src]
                used.discard(target)
        return False

    return visit(0, {}, set())


def _rewired_trimer():
    product = _graph(PS_TRIMER)
    graph, _ = _permuted(_double_edge_swap(product, *_ring_bonds(product)), 7)
    return graph, product


def test_reference_search_blows_the_budget_on_the_rewired_trimer():
    """The reproduction: the pre-I-303 search does not decide this case
    within the budget (it would enumerate the rings' symmetries)."""
    graph, product = _rewired_trimer()
    with pytest.raises(_Spent):
        _within(BUDGET, _reference_match, graph, set(graph), product, {})


def test_symmetric_aromatic_mismatch_is_decided_within_a_budget():
    """Two phenyl rings of a PS trimer rewired into one twelve-ring: every
    degree and signature agrees and nothing is pinned. The match must be
    refused without enumerating the rings' symmetries."""
    graph, product = _rewired_trimer()
    assert _within(BUDGET, cm._exact_component_match, graph, set(graph), product, {}) is False


def test_symmetric_aromatic_match_is_found_within_a_budget():
    product = _graph(PS_TRIMER)
    graph, _ = _permuted(product, 11)
    assert _within(BUDGET, cm._exact_component_match, graph, set(graph), product, {}) is True


def test_heavy_pinned_hydrogens_free():
    """The shape rewrite_problems passes: every heavy atom pinned, hydrogens
    free.  A permuted copy matches, a rewired one does not."""
    product = _graph(PS_TRIMER)
    for source, expected in ((product, True),
                             (_double_edge_swap(product, *_ring_bonds(product)), False)):
        graph, perm = _permuted(source, 3)
        pins = {perm[i]: i for i in product if product[i]["element"] != "H"}
        assert _within(BUDGET, cm._exact_component_match, graph, set(graph), product,
                       pins) is expected


SMALL = ["c1ccccc1", "Cc1ccccc1", "CC(C)C", "C=CC=C", "C1CC1", "CC(=O)O", "[CH2]C=C",
         "C1=CC=CC1", "OCCO", "C#CC", "c1ccc2ccccc2c1", "CC(c1ccccc1)C", "C[CH]c1ccccc1",
         "C1CCC2CCCCC2C1"]


@pytest.mark.parametrize("smiles", SMALL)
def test_verdicts_equal_the_reference_search(smiles):
    """Wherever the pre-I-303 search finishes, the verdict is its verdict:
    permuted copies, degree-preserving rewires, partial, wrong and complete
    pins, sub-components of the wrong size."""
    rng = random.Random(smiles)
    product = _graph(smiles)
    heavy = [i for i in sorted(product) if product[i]["element"] != "H"]
    bonds = sorted((i, j) for i in product for j in product[i]["edges"] if i < j)
    cases = []
    for k in range(24):
        source = product
        if k % 3 == 1 and len(bonds) > 1:
            ab, cd = rng.sample(bonds, 2)
            if len({*ab, *cd}) == 4 and cd[1] not in product[ab[0]]["edges"] \
                    and ab[1] not in product[cd[0]]["edges"]:
                source = _double_edge_swap(product, ab, cd)
        graph, perm = _permuted(source, rng.random())
        chosen = heavy if k % 6 == 0 else rng.sample(heavy, rng.randint(0, len(heavy)))
        pins = {perm[i]: i for i in chosen}
        if k % 4 == 3 and len(heavy) > 1:
            a, b = rng.sample(heavy, 2)
            pins[perm[a]] = b  # a wrong pin
        cases.append((graph, set(graph), pins))
        if k % 5 == 4:
            cases.append((graph, set(list(graph)[:-1]), {}))
    verdicts = set()
    for graph, members, pins in cases:
        expected = _reference_match(graph, members, product, pins)
        assert cm._exact_component_match(graph, members, product, pins) is expected
        verdicts.add(expected)
    assert verdicts == {True, False}, "the cases exercise both verdicts"


def test_isolated_hydrogens_and_empty():
    h = {0: {"element": "H", "radical": 1, "charge": 0, "lone_pairs": 0,
             "implicit_hydrogens": 0, "edges": {}}}
    assert cm._exact_component_match(h, {0}, h, {}) is True
    assert cm._exact_component_match({}, set(), {}, {}) is True
    h2 = _graph("[H][H]")
    assert cm._exact_component_match(h2, set(h2), h2, {}) is True
    assert cm._exact_component_match(h, {0}, h2, {}) is False
