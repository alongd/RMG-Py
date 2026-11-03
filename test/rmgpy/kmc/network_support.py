"""Shared builders for the network (schema v2, ``architecture: network``)
tests: the cured phenolic network fixture with an inert filler
(``fixtures/decks/phenolic_network_deck.py``) and mutations of it."""

import copy
from pathlib import Path

from polymer_spec_v2_support import polymer_keyword

DECKS = Path(__file__).parent / "fixtures" / "decks"
NETWORK_DECK = DECKS / "phenolic_network_deck.py"


def network_mapping():
    """The network fixture's declaration, as a fresh mutable mapping."""
    return copy.deepcopy(polymer_keyword(NETWORK_DECK, "polymer_spec"))


def network_spec(mutate=None):
    from rmgpy.polymer_spec import normalized_polymer_spec
    mapping = network_mapping()
    if mutate is not None:
        mutate(mapping)
    return normalized_polymer_spec(mapping)


def no_filler(mapping):
    mapping.pop("inert_phases", None)


def with_met(mapping):
    import melt_engine_support as S
    mapping["transport_model"] = {"met": S.pe_met()}


def counts(residues=None, links=None, seed=None):
    def mutate(mapping):
        ensemble = mapping["initial_ensemble"]
        if residues is not None:
            ensemble["representative_residue_counts"] = dict(residues)
        if links is not None:
            ensemble["link_counts"] = dict(links)
        if seed is not None:
            ensemble["seed"] = seed
    return mutate
