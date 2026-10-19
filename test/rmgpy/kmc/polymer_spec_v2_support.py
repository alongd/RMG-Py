"""Shared builders for the schema-v2 tests: the v2 EPDM fixture declaration
(``fixtures/decks/epdm_v2_deck.py``) and mutations of it."""

import ast
import copy
from pathlib import Path

from rmgpy.polymer_spec import PolymerSpec

DECKS = Path(__file__).parent / "fixtures" / "decks"
V2_DECK = DECKS / "epdm_v2_deck.py"


def polymer_keyword(path, keyword):
    """The literal ``keyword`` argument of a deck's ``polymer(...)`` call."""
    tree = ast.parse(Path(path).read_text())
    call = next(node.value for node in tree.body if isinstance(node, ast.Expr)
                and isinstance(node.value, ast.Call) and getattr(node.value.func, "id", "") == "polymer")
    return ast.literal_eval(next(k.value for k in call.keywords if k.arg == keyword))


def v2_mapping():
    """The v2 EPDM fixture's declaration, as a fresh mutable mapping."""
    return copy.deepcopy(polymer_keyword(V2_DECK, "polymer_spec"))


def v2_spec(mutate=None):
    mapping = v2_mapping()
    if mutate is not None:
        mutate(mapping)
    return PolymerSpec.from_mapping(mapping)


def q(numerator, denominator):
    return {"numerator": numerator, "denominator": denominator}


# a test composition with N absent, and a Markov chain whose stationary
# distribution it is (E: 3/4 stay, P: 1/2 stay)
EP_COMPOSITION = {"E": q(2, 3), "P": q(1, 3), "N": q(0, 1)}
EP_ROWS = {"E": {"E": q(3, 4), "P": q(1, 4)}, "P": {"E": q(1, 2), "P": q(1, 2)}}


def markov(rows, composition=None):
    def mutate(mapping):
        if composition is not None:
            mapping["initial_ensemble"]["composition"] = composition
        mapping["initial_ensemble"]["sequence_statistics"] = {"model": "markov1", "transition": rows}
    return mutate
