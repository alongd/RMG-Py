"""A provider constructed for a declaration stays registered for the life of the process
(``boundaries.PROVIDERS``, named by the declaration's digest), so that every later lookup
of that name resolves.  A lookup by repeat unit alone (``emit.match_provider`` without a
declaration) must not see it: what it returns would depend on which declarations the
process happened to build before (I-270)."""

import json
from pathlib import Path

import pytest

from rmgpy.kmc import boundaries as b
from rmgpy.kmc import emit
from rmgpy.kmc.declared_chain import DeclaredChainProvider
from rmgpy.kmc.run import read_deck
from rmgpy.polymer_spec import normalized_linear_homopolymer_spec

DECKS = Path(__file__).parent / "fixtures" / "decks"


@pytest.fixture
def registry():
    """The registry as it was, restored in place (the modules hold the dict itself)."""
    saved = dict(b.PROVIDERS)
    yield b.PROVIDERS
    b.PROVIDERS.clear()
    b.PROVIDERS.update(saved)


def _pp_and_another():
    """The PP deck's declaration and another declaration of the same repeat unit."""
    spec = normalized_linear_homopolymer_spec(read_deck(DECKS / "pp_deck.py").spec)
    mapping = json.loads(spec.normalized_json)
    mapping["label"] = "another declaration of the PP repeat unit"
    return spec, normalized_linear_homopolymer_spec(json.dumps(mapping))


def test_a_provider_built_for_a_declaration_matches_only_that_declaration(registry):
    spec, other = _pp_and_another()
    built = DeclaredChainProvider(other)
    assert registry[built.name] is built
    assert emit.match_provider(spec.repeat_unit) is None
    assert emit.match_provider(spec.repeat_unit, other) is built
    assert emit.declared_provider(spec) is not built


def test_a_roster_does_not_depend_on_the_declarations_built_before_it(registry):
    spec, other = _pp_and_another()
    DeclaredChainProvider(other)
    roster = emit.reference_systems(spec, emit.match_provider(spec.repeat_unit), dps=(1,),
                                    small=False, monomer=False, molecularity=frozenset((1,)))
    owners = {entry["chain"]["owner"] for entries in roster.systems.values()
              for entry in entries if "chain" in entry}
    systems, refusals = emit.build_reference_systems(roster)
    assert owners == {emit.declared_provider(spec).name}
    assert refusals == [] and len(systems) == len(roster.systems)
