"""Decks beyond PS and PE (``fixtures/decks``): polypropylene, an alternating
E-P-ENB approximation of EPDM and a linear ortho-ortho novolac, read through
the real deck reader and built by the generic declared-chain path."""

from pathlib import Path

import pytest

from rmgpy.kmc import emit
from rmgpy.kmc import windows as win
from rmgpy.kmc.run import read_deck

DECKS = Path(__file__).parent / "fixtures" / "decks"
DECK_NAMES = ("pp_deck.py", "epdm_deck.py", "novolac_deck.py")


def _heavy(graph_atoms):
    return sum(1 for atom in graph_atoms if atom.element != "H")


@pytest.mark.parametrize("name, label, repeat_heavy", [
    ("pp_deck.py", "polypropylene (CH2-CH(CH3))", 3),
    ("epdm_deck.py", "EPDM (alternating E-P-ENB repeat)", 14),
    ("novolac_deck.py", "linear ortho-ortho novolac (CH2-C6H2(OH))", 8),
])
def test_each_deck_loads_through_the_deck_reader(name, label, repeat_heavy):
    deck = read_deck(DECKS / name)
    assert deck.label == deck.spec.label == label
    assert _heavy(deck.spec.repeat_unit.atoms) == repeat_heavy
    assert set(deck.digests()) == {"chemistry_sha256", "declaration_sha256"}
    assert emit.match_provider(deck.spec.repeat_unit) is None


@pytest.mark.parametrize("name, dp", [(name, 1) for name in DECK_NAMES] + [("pp_deck.py", 3)])
def test_a_declared_chain_of_dp_repeat_units_holds_dp_repeat_units(name, dp):
    """``dp`` counts repeat units whatever the primitive period: a repeat of
    several backbone bonds whose period is the whole repeat is not built
    ``d_rep`` times over."""
    deck = read_deck(DECKS / name)
    declared = emit.declared_provider(deck.spec)
    build = win.build_spec(declared.spec_from_declaration(dp), "chain")
    heavy = sum(1 for atom in build.molecule.atoms if atom.element.symbol != "H")
    assert heavy == dp * _heavy(deck.spec.repeat_unit.atoms)


@pytest.mark.parametrize("dp", (1, 2, 3))
def test_a_period_shorter_than_the_repeat_still_counts_repeat_units(dp):
    """PE: the primitive period (CH2) is half the repeat (CH2-CH2)."""
    import melt_engine_support as S
    from rmgpy.polymer_spec import LinearHomopolymerSpec

    spec = LinearHomopolymerSpec.from_mapping(S.pe_mapping({10: 1}))
    declared = emit.declared_provider(spec)
    assert (declared.d_rep, declared.period) == (2, 1)
    build = win.build_spec(declared.spec_from_declaration(dp), "chain")
    assert sum(1 for atom in build.molecule.atoms if atom.element.symbol != "H") == 2 * dp
