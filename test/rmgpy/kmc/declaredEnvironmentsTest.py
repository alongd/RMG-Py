"""The environment set of a linear copolymer declaration
(``rmgpy.kmc.declared_environments``): one chain period that holds every
k-residue word of the sequence support covers every site-class environment
a chain of the declaration can have, at k = 1 + ceil(2(r + 1) / n_min)."""

import random
from fractions import Fraction
from itertools import product

import pytest

from rmgpy.kmc import declared_environments as de
from rmgpy.kmc.declared_chain import ChainGeometry
from rmgpy.kmc.run import read_deck
from rmgpy.polymer_spec import LinkDeclaration, ResidueDeclaration, linear_junctions

from polymer_spec_v2_support import DECKS, EP_COMPOSITION, EP_ROWS, markov, q, v2_spec


def _ep_only(mapping):
    del mapping["residues"]["N"]
    mapping["sequence_support"] = {"E": ["E", "P"], "P": ["E", "P"]}
    mapping["initial_ensemble"]["composition"] = {"E": q(1, 2), "P": q(1, 2)}


def _environments(spec, period, whole_chain=False):
    return de.period_environments(dict(spec.residues), spec.junctions, spec.locality_policy,
                                  tuple(period), whole_chain=whole_chain)


@pytest.fixture(scope="module")
def epdm():
    return v2_spec()


@pytest.fixture(scope="module")
def epdm_environments(epdm):
    return de.sequence_environments(epdm)


def test_the_word_order_keeps_the_boundary_shell():
    assert de.word_order(2, 2) == 4
    assert de.word_order(2, 3) == 3
    assert de.word_order(0, 1) == 3


def test_the_epdm_word_order_is_four(epdm):
    assert de.spec_word_order(epdm) == 4


@pytest.mark.parametrize("k", [2, 3, 4])
def test_the_covering_period_holds_every_word_of_a_complete_support(k):
    support = {"E": ("E", "N", "P"), "N": ("E", "N", "P"), "P": ("E", "N", "P")}
    period = de.covering_period(support, k)
    cyclic = period * k
    words = {tuple(cyclic[i:i + k]) for i in range(len(period))}
    assert words == set(product("ENP", repeat=k)) and len(period) == 3 ** k


def test_the_covering_period_follows_a_restricted_support():
    support = {"E": ("E", "N", "P"), "N": ("E", "P"), "P": ("E", "N", "P")}
    period = de.covering_period(support, 4)
    cyclic = period * 4
    words = {tuple(cyclic[i:i + 4]) for i in range(len(period))}
    allowed = {w for w in product("ENP", repeat=4)
               if all(b in support[a] for a, b in zip(w, w[1:]))}
    assert words == allowed
    assert all(b in support[a] for a, b in zip(period, period[1:] + period[:1]))


def test_the_environment_set_equals_the_union_over_random_chains(epdm, epdm_environments):
    rng = random.Random(1)
    union = set()
    for _ in range(20):
        union |= set(_environments(epdm, [rng.choice("EPN") for _ in range(10)]))
    assert len(epdm_environments) == 123
    assert set(epdm_environments) == union


def test_two_residues_need_k_four_and_three_misses_three():
    spec = v2_spec(_ep_only)
    four = de.sequence_environments(spec)
    three = set(_environments(spec, de.covering_period(dict(spec.sequence_support), 3)))
    assert len(four) == 18
    assert three < set(four) and len(set(four) - three) == 3


def test_the_per_residue_union_misses_the_junction_environments(epdm, epdm_environments):
    per_residue = set()
    for name in "ENP":
        per_residue |= set(_environments(epdm, [name]))
    assert len(per_residue) == 13
    assert per_residue < set(epdm_environments)


@pytest.mark.parametrize("sequence", ["E", "P", "N", "EN", "ENP", "EEPP", "EEEPEPPP", "NN"])
def test_ball_keys_equal_whole_chain_keys(epdm, sequence):
    assert _environments(epdm, sequence) == _environments(epdm, sequence, whole_chain=True)


def test_every_witness_word_fits_in_k_residues(epdm, epdm_environments):
    k = de.spec_word_order(epdm)
    for occurrences in epdm_environments.values():
        for word, residue, atom in occurrences:
            assert 1 <= len(word) <= k
            assert residue in word
            assert atom in {a.key for a in dict(epdm.residues)[residue].graph.atoms}


def _lift(spec):
    """The v1 repeat unit as the one residue of a homopolymer support."""
    unit = spec.repeat_unit.as_mapping()
    ports = unit.pop("ports")
    residue = ResidueDeclaration.from_mapping({"graph": unit, "sites": {
        "head": dict(ports["head"], **{"class": "entry", "vacancy": "forbidden"}),
        "tail": dict(ports["tail"], **{"class": "exit", "vacancy": "forbidden"})}})
    residues = {"R": residue}
    links = {"j": LinkDeclaration.from_mapping(
        {"from": "exit", "to": "entry", "bond_order": ports["tail"]["bond_order"]})}
    return residues, linear_junctions(residues, links, {"R": ("R",)})


@pytest.mark.parametrize("name, classes", [
    ("pe_toy_deck.py", 1), ("pp_deck.py", 3), ("epdm_deck.py", 14), ("novolac_deck.py", 6)])
def test_a_one_residue_lift_of_each_v1_deck_has_its_atom_classes(name, classes):
    """Declared bond orders: a Kekule chain gives novolac 8 environments."""
    spec = read_deck(DECKS / name).spec
    residues, junctions = _lift(spec)
    found = de.period_environments(residues, junctions, spec.locality_policy, ("R",))
    assert len(found) == ChainGeometry.of(spec).heavy_classes == classes


def test_ideal_random_word_frequencies(epdm):
    frequencies = de.word_frequencies(epdm, 2)
    assert len(frequencies) == 9 and sum(frequencies.values()) == 1
    assert frequencies[("E", "P")] == Fraction(1, 2) * Fraction(2, 5)
    assert frequencies[("N", "N")] == Fraction(1, 100)


def test_alternating_word_frequencies():
    def alternating(mapping):
        mapping["initial_ensemble"]["composition"] = {"E": q(1, 2), "P": q(1, 2), "N": 0}
        mapping["initial_ensemble"]["sequence_statistics"] = {"model": "alternating", "cycle": ["E", "P"]}
    assert de.word_frequencies(v2_spec(alternating), 3) == {
        ("E", "P", "E"): Fraction(1, 2), ("P", "E", "P"): Fraction(1, 2)}


def test_markov_word_frequencies():
    frequencies = de.word_frequencies(v2_spec(markov(EP_ROWS, EP_COMPOSITION)), 3)
    assert sum(frequencies.values()) == 1
    assert frequencies[("E", "E", "P")] == Fraction(2, 3) * Fraction(3, 4) * Fraction(1, 4)
    assert all(word.count("N") == 0 for word in frequencies)
