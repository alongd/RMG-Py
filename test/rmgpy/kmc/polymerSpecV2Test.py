"""Declaration schema v2 (``PolymerSpec``): linear statistical copolymers.

Backward compatibility (every v1 deck keeps its bytes and digests), the
validation rules (one refusal per rule), the digests (composition and
statistics outside the chemistry digest, the sequence support inside it),
the deck reader and the named refusals of every v1 consumer."""

import copy
import json
from fractions import Fraction
from types import SimpleNamespace

import pytest

from rmgpy.kmc import run as kr
from rmgpy.polymer_spec import (
    LinearHomopolymerSpec,
    PolymerSpec,
    SchemaV2Refusal,
    normalized_polymer_spec,
)

from polymer_spec_v2_support import (
    DECKS, EP_COMPOSITION, EP_ROWS, V2_DECK, markov, polymer_keyword, q, v2_mapping, v2_spec)


# declaration / chemistry sha256 of the v1 decks before schema v2 existed (PE and PP declarations
# re-pinned when their transport_model gained cited values; their chemistry pins are the originals)
PINS = {
    "pe_toy_deck.py": ("fa994b4868a4392bccc9199e3bd3b7d308a510d94206d2f8efcc91eb3c69a4c6",
                       "0f20e6eade428a907610eb8a36696142dc1b5e6d4623e6fd96992e63d2a9f82c"),
    "pp_deck.py": ("6a4b6f421de35fe5afc86e521943a19d5b4386899e3c31129d342cdf8bf88f4d",
                   "8b978f4744f96c096087b29ec12c77d9fe515dbecfc3f8394216a5b085c43ad1"),
    "epdm_deck.py": ("9d5fa1dfc9708f0a57d1383cd6b4d96d7a97f085dac45ae5a7e5faeb8e46feac",
                     "25ca8f7d516f36885617b66d6076c71476296df9a29f6cd3a3b1cc631fdf1dcc"),
    "novolac_deck.py": ("614e0a16a4878c671b0e4dc21ba4128b0d9ae591ee05fde58b0c92c0759b955f",
                        "d5368fe1d64202d3cb3d0a9760763293d04b5f6a1ffe6bc548e43be3ee3e66c2"),
}
REFUSED = r"^refused:schema_v2: {} reads v1 declarations only$"


# --------------------------------------------------------------------------
# backward compatibility
# --------------------------------------------------------------------------
@pytest.mark.parametrize("name", sorted(PINS))
def test_v1_deck_digests_are_unchanged(name):
    deck = kr.read_deck(DECKS / name)
    assert type(deck.spec) is LinearHomopolymerSpec
    declaration, chemistry = PINS[name]
    assert deck.digests() == {"declaration_sha256": declaration, "chemistry_sha256": chemistry}


@pytest.mark.parametrize("name", sorted(PINS))
def test_v1_declarations_are_refused_by_the_v2_reader(name):
    """A v1 mapping never parses into ``PolymerSpec``."""
    with pytest.raises(ValueError, match="schema_version"):
        PolymerSpec.from_mapping(polymer_keyword(DECKS / name, "linear_homopolymer_spec"))


# --------------------------------------------------------------------------
# the v2 fixture and normalization
# --------------------------------------------------------------------------
def test_the_v2_fixture_reads_through_the_deck_reader():
    deck = kr.read_deck(V2_DECK)
    assert type(deck.spec) is PolymerSpec
    assert deck.label == deck.spec.label
    assert deck.spec.architecture == "linear"
    assert [name for name, _residue in deck.spec.residues] == ["E", "N", "P"]
    assert deck.spec.initial_ensemble.composition == (
        ("E", Fraction(1, 2)), ("N", Fraction(1, 10)), ("P", Fraction(2, 5)))
    assert set(deck.digests()) == {"chemistry_sha256", "declaration_sha256"}


def test_normalized_json_is_a_fixed_point():
    """Through the mapping and through the bytes, floats in policy blocks included."""
    def floats(mapping):
        mapping["diagnostics"] = {"note": 700.0, "fraction": 0.1}
    spec = v2_spec(floats)
    assert PolymerSpec.from_mapping(json.loads(spec.normalized_json)).normalized_json == spec.normalized_json
    assert normalized_polymer_spec(spec.normalized_json).normalized_json == spec.normalized_json
    assert normalized_polymer_spec(spec) is spec


def test_junctions_join_each_tail_site_to_the_next_head_site():
    spec = v2_spec()
    assert spec.junctions[("E", "P")] == ("e2", "p1", 1.0)
    assert spec.junctions[("N", "E")] == ("n3", "e1", 1.0)
    assert len(spec.junctions) == 9


# --------------------------------------------------------------------------
# validation: one refusal per rule
# --------------------------------------------------------------------------
def _composition(e, p, n):
    def mutate(mapping):
        mapping["initial_ensemble"]["composition"] = {"E": e, "P": p, "N": n}
    return mutate


def test_composition_must_sum_to_one():
    with pytest.raises(ValueError, match="sum exactly to one"):
        v2_spec(_composition(q(1, 2), q(2, 5), q(1, 5)))


def test_a_float_composition_is_refused():
    with pytest.raises(ValueError, match="exact rational"):
        v2_spec(_composition(0.5, q(2, 5), q(1, 10)))


def test_a_float_inside_a_rational_is_refused():
    with pytest.raises(ValueError, match="exact rational"):
        v2_spec(_composition({"numerator": 1.0, "denominator": 2}, q(2, 5), q(1, 10)))


def test_a_markov_chain_whose_stationary_distribution_is_the_composition_is_accepted():
    spec = v2_spec(markov(EP_ROWS, EP_COMPOSITION))
    assert spec.initial_ensemble.statistics.model == "markov1"


def test_a_markov_stationary_distribution_that_differs_from_the_composition_is_refused():
    rows = {"E": {"E": q(1, 2), "P": q(1, 2)}, "P": {"E": q(1, 2), "P": q(1, 2)}}
    with pytest.raises(ValueError, match="stationary"):
        v2_spec(markov(rows, EP_COMPOSITION))


def test_a_markov_row_that_does_not_sum_to_one_is_refused():
    rows = {"E": {"E": q(3, 4), "P": q(1, 2)}, "P": {"E": q(1, 2), "P": q(1, 2)}}
    with pytest.raises(ValueError, match="sums to"):
        v2_spec(markov(rows, EP_COMPOSITION))


def test_a_markov_entry_outside_the_support_is_refused():
    def mutate(mapping):
        markov(EP_ROWS, EP_COMPOSITION)(mapping)
        mapping["sequence_support"]["P"] = ["E", "N"]
    with pytest.raises(ValueError, match="outside sequence_support"):
        v2_spec(mutate)


def test_a_reducible_markov_chain_is_refused():
    rows = {"E": {"E": q(1, 1)}, "P": {"P": q(1, 1)}}
    with pytest.raises(ValueError, match="irreducible"):
        v2_spec(markov(rows, EP_COMPOSITION))


def _alternating(cycle, composition):
    def mutate(mapping):
        mapping["initial_ensemble"]["composition"] = composition
        mapping["initial_ensemble"]["sequence_statistics"] = {"model": "alternating", "cycle": cycle}
    return mutate


def test_an_alternating_cycle_with_its_own_composition_is_accepted():
    spec = v2_spec(_alternating(["E", "P"], {"E": q(1, 2), "P": q(1, 2), "N": 0}))
    assert spec.initial_ensemble.statistics.cycle == ("E", "P")


def test_an_alternating_cycle_with_an_inconsistent_composition_is_refused():
    with pytest.raises(ValueError, match="cycle"):
        v2_spec(_alternating(["E", "P"], {"E": q(2, 3), "P": q(1, 3), "N": 0}))


def test_ideal_random_needs_every_pair_of_present_residues_in_the_support():
    def mutate(mapping):
        mapping["sequence_support"]["N"] = ["E", "P"]
    with pytest.raises(ValueError, match="ideal_random"):
        v2_spec(mutate)


def test_a_site_class_missing_from_attachment_classes_is_refused():
    def mutate(mapping):
        mapping["locality_policy"]["bounds"]["attachment_classes"] = ["vinyl_head"]
    with pytest.raises(ValueError, match="attachment_classes"):
        v2_spec(mutate)


def test_a_sequence_model_outside_the_declared_vocabulary_is_refused():
    def mutate(mapping):
        mapping["locality_policy"]["bounds"]["topology_metadata_schema"]["sequence_models"] = ["markov1"]
    with pytest.raises(ValueError, match="topology_metadata_schema"):
        v2_spec(mutate)


def test_a_one_residue_declaration_is_refused_as_v1_expressible():
    def mutate(mapping):
        mapping["residues"] = {"E": mapping["residues"]["E"]}
        mapping["sequence_support"] = {"E": ["E"]}
        mapping["initial_ensemble"]["composition"] = {"E": 1}
    with pytest.raises(ValueError, match="declare it as v1"):
        v2_spec(mutate)


@pytest.mark.parametrize("architecture", ["network"])
def test_a_network_architecture_is_refused_by_name(architecture):
    """Networks are S7; a branched declaration has its own class and tests
    (``branchedSpecTest``)."""
    def mutate(mapping):
        mapping["architecture"] = architecture
        mapping["locality_policy"]["bounds"]["topology_metadata_schema"]["architecture"] = [architecture]
    with pytest.raises(ValueError, match=f"refused:architecture: {architecture}"):
        v2_spec(mutate)


def _ladder(mapping):
    """Two links between consecutive residues: a ladder (rule 2)."""
    carbon = {"element": "C", "hydrogens": 1, "electronic_state": {"radical_electrons": 0}}
    rung = {"graph": {"atoms": [dict(carbon, key="a"), dict(carbon, key="b")],
                      "bonds": [{"atoms": ["a", "b"], "order": 1}],
                      "electronic_state": {"multiplicity": 1}},
            "sites": {"in_a": {"atom": "a", "bond_order": 1, "class": "in_a", "vacancy": "forbidden"},
                      "out_a": {"atom": "a", "bond_order": 1, "class": "out_a", "vacancy": "forbidden"},
                      "in_b": {"atom": "b", "bond_order": 1, "class": "in_b", "vacancy": "forbidden"},
                      "out_b": {"atom": "b", "bond_order": 1, "class": "out_b", "vacancy": "forbidden"}}}
    mapping["residues"] = {"L": rung, "M": copy.deepcopy(rung)}
    mapping["links"] = {"la": {"from": "out_a", "to": "in_a", "bond_order": 1},
                        "lb": {"from": "out_b", "to": "in_b", "bond_order": 1}}
    mapping["sequence_support"] = {"L": ["L", "M"], "M": ["L", "M"]}
    mapping["locality_policy"]["bounds"]["attachment_classes"] = ["in_a", "in_b", "out_a", "out_b"]
    mapping["initial_ensemble"]["composition"] = {"L": q(1, 2), "M": q(1, 2)}


def test_a_link_that_puts_an_inter_residue_bond_on_a_ring_is_refused():
    with pytest.raises(ValueError, match="rule 2:"):
        v2_spec(_ladder)


def _conjugated(mapping):
    """Vinylene residues joined so one pi system crosses the link (rule 2b)."""
    carbon = {"element": "C", "hydrogens": 1, "electronic_state": {"radical_electrons": 0}}
    vinylene = {"graph": {"atoms": [dict(carbon, key="a"), dict(carbon, key="b")],
                          "bonds": [{"atoms": ["a", "b"], "order": 2}],
                          "electronic_state": {"multiplicity": 1}},
                "sites": {"head": {"atom": "a", "bond_order": 1, "class": "vinyl_head", "vacancy": "forbidden"},
                          "tail": {"atom": "b", "bond_order": 1, "class": "vinyl_tail", "vacancy": "forbidden"}}}
    mapping["residues"]["V"] = vinylene
    for successors in mapping["sequence_support"].values():
        successors.append("V")
    mapping["sequence_support"]["V"] = ["E", "P", "N", "V"]
    mapping["initial_ensemble"]["composition"]["V"] = 0


def test_a_link_inside_one_pi_system_is_refused():
    with pytest.raises(ValueError, match="rule 2b:"):
        v2_spec(_conjugated)


def test_a_disconnected_support_is_refused():
    def mutate(mapping):
        mapping["sequence_support"] = {"E": ["E", "P"], "P": ["P"], "N": ["N"]}
        mapping["initial_ensemble"]["composition"] = {"E": 0, "P": 1, "N": 0}
    with pytest.raises(ValueError, match="strongly connected"):
        v2_spec(mutate)


# --------------------------------------------------------------------------
# digests
# --------------------------------------------------------------------------
def test_composition_and_statistics_stay_outside_the_chemistry_digest():
    base = v2_spec()
    other = v2_spec(_composition(q(1, 3), q(1, 3), q(1, 3)))
    blocky = v2_spec(markov(EP_ROWS, EP_COMPOSITION))
    assert kr.chemistry_sha256(base) == kr.chemistry_sha256(other) == kr.chemistry_sha256(blocky)
    assert len({kr.declaration_sha256(s) for s in (base, other, blocky)}) == 3


def test_the_sequence_support_enters_the_chemistry_digest():
    composition = {"E": q(1, 2), "P": q(1, 2), "N": 0}
    full = v2_spec(_composition(*composition.values()))

    def fewer(mapping):
        _composition(*composition.values())(mapping)
        mapping["sequence_support"]["N"] = ["E", "P"]
    assert kr.chemistry_sha256(v2_spec(fewer)) != kr.chemistry_sha256(full)


# --------------------------------------------------------------------------
# the deck reader and the input file
# --------------------------------------------------------------------------
def test_a_deck_carrying_both_declarations_is_refused(tmp_path):
    v1 = polymer_keyword(DECKS / "pp_deck.py", "linear_homopolymer_spec")
    path = tmp_path / "both.py"
    path.write_text(f"polymer({v1['label']!r}, linear_homopolymer_spec={v1!r}, "
                    f"polymer_spec={v2_mapping()!r})\n")
    with pytest.raises(kr.RunError, match="both linear_homopolymer_spec and polymer_spec"):
        kr.read_deck(path)


def test_a_deck_carrying_neither_declaration_is_refused(tmp_path):
    path = tmp_path / "none.py"
    path.write_text("polymer('x', cutoff=4)\n")
    with pytest.raises(kr.RunError, match="carries no linear_homopolymer_spec or polymer_spec"):
        kr.read_deck(path)


def test_the_input_file_polymer_refuses_both_declarations():
    from rmgpy.exceptions import InputError
    from rmgpy.rmg import input as rmg_input

    v1 = polymer_keyword(DECKS / "pp_deck.py", "linear_homopolymer_spec")
    with pytest.raises(InputError, match="exactly one"):
        rmg_input.polymer(v1["label"], linear_homopolymer_spec=v1, polymer_spec=v2_mapping())


def test_the_input_file_polymer_refuses_a_v2_declaration_by_name():
    from rmgpy.rmg import input as rmg_input

    mapping = v2_mapping()
    with pytest.raises(SchemaV2Refusal, match=REFUSED.format("Polymer")):
        rmg_input.polymer(mapping["label"], polymer_spec=mapping)


# --------------------------------------------------------------------------
# named refusals at every v1 consumer
# --------------------------------------------------------------------------
def test_polymer_refuses_a_v2_declaration():
    from rmgpy.polymer import Polymer

    for keyword in ("polymer_spec", "linear_homopolymer_spec"):
        with pytest.raises(SchemaV2Refusal, match=REFUSED.format("Polymer")):
            Polymer(label="x", monomer=None, **{keyword: v2_spec()})


def test_the_shared_support_refuses_a_v2_declaration():
    from rmgpy.polymer_shared_support import _authoritative_spec

    pool = SimpleNamespace(label="x", linear_homopolymer_spec=v2_spec())
    with pytest.raises(SchemaV2Refusal,
                       match=REFUSED.format("polymer_shared_support._authoritative_spec")):
        _authoritative_spec(pool)


def test_chain_geometry_refuses_a_v2_declaration():
    from rmgpy.kmc.declared_chain import ChainGeometry

    with pytest.raises(SchemaV2Refusal, match=REFUSED.format(r"ChainGeometry\.of")):
        ChainGeometry.of(v2_spec())


def test_declared_proxies_refuse_a_v2_declaration():
    from rmgpy.kmc import declared_proxies

    for function in (declared_proxies.declared_proxy_set, declared_proxies.declared_structures):
        with pytest.raises(SchemaV2Refusal, match=REFUSED.format("declared_proxies")):
            function(v2_spec(), 3)


def test_emit_refuses_a_v2_declaration_in_a_legacy_roster():
    """``declared_provider`` and ``declared_reference_systems`` read v2
    (``sequenceProxiesTest``); the legacy DP roster does not."""
    from rmgpy.kmc import emit

    with pytest.raises(SchemaV2Refusal, match=REFUSED.format(r"emit\.reference_systems")):
        emit.reference_systems(v2_spec(), None)


def test_load_roster_refuses_a_v2_declaration(tmp_path):
    from rmgpy.kmc import emit

    path = tmp_path / "roster.json"
    path.write_text(json.dumps({"schema": emit.ROSTER_SCHEMA_V3, "polymer": json.loads(v2_spec().normalized_json),
                                "proxy_units": 3, "systems": {}}))
    with pytest.raises(SchemaV2Refusal, match=REFUSED.format(r"emit\.load_roster")):
        emit.load_roster(path)


def test_rationals_are_exact_after_normalization():
    """Composition and transition entries keep exact rationals in the bytes."""
    rows = {"E": {"E": q(6, 8), "P": q(1, 4)}, "P": {"E": q(1, 2), "P": q(2, 4)}}
    spec = v2_spec(markov(rows, EP_COMPOSITION))
    data = json.loads(spec.normalized_json)["initial_ensemble"]
    assert data["sequence_statistics"]["transition"]["E"]["E"] == {"numerator": 3, "denominator": 4}
    assert data["composition"]["N"] == {"numerator": 0, "denominator": 1}
    assert data["seed"] == 1 and data["sampler"] == "markov_chain_v1"
