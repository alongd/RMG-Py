"""Declaration schema v2, ``architecture: network`` (S7): a cured network
above the gel point, its initial state the declared ``random_bond_v1``
generator, fillers declared as ``inert_phases`` outside the chemistry digest;
the method of moments refuses a network by name (``refused:network_topology``)
and every v1, linear-v2 and branched declaration keeps its bytes."""

from fractions import Fraction

import pytest

from branched_support import BRANCHED_DECK, branched_spec
from network_support import DECKS, NETWORK_DECK, counts, network_mapping, network_spec, no_filler

from rmgpy.kmc import run as kr

# declaration / chemistry sha256 of every fixture deck at the S7 base (538658d6a)
PINS = {
    "ep_v2_deck.py": ("ddb53e4b1edd219abf5182d9502f8b84fb49df59997337226ed96e817bab550c",
                      "4dbceba4482bdeeffbf45acfe0de6d447538df7745218e89b10465c86559b51b"),
    "epdm_deck.py": ("9d5fa1dfc9708f0a57d1383cd6b4d96d7a97f085dac45ae5a7e5faeb8e46feac",
                     "25ca8f7d516f36885617b66d6076c71476296df9a29f6cd3a3b1cc631fdf1dcc"),
    "epdm_v2_deck.py": ("33664ce5b12206ed9f61d8d9ff1e0abc3728296840b0a07b10fbe32007a0a9f4",
                        "a0603a35c4e7abd2bf251f2fe2991b4e6a7e3a206c97392c77ca941887f23dc2"),
    "novolac_deck.py": ("614e0a16a4878c671b0e4dc21ba4128b0d9ae591ee05fde58b0c92c0759b955f",
                        "d5368fe1d64202d3cb3d0a9760763293d04b5f6a1ffe6bc548e43be3ee3e66c2"),
    "novolac_v2_deck.py": ("ff1167e9a31dd498d2e6df1c49869c1f594e339840b19a946a12ada607e9c2e2",
                           "1d937335a373be55c18e0451ff1bd583c0fcc18f715dcb8937cd022af8077b44"),
    "pe_toy_deck.py": ("fa994b4868a4392bccc9199e3bd3b7d308a510d94206d2f8efcc91eb3c69a4c6",
                       "0f20e6eade428a907610eb8a36696142dc1b5e6d4623e6fd96992e63d2a9f82c"),
    "pp_deck.py": ("6a4b6f421de35fe5afc86e521943a19d5b4386899e3c31129d342cdf8bf88f4d",
                   "8b978f4744f96c096087b29ec12c77d9fe515dbecfc3f8394216a5b085c43ad1"),
    "ps_deck.py": ("039b0665a904b19dfc8b8e05ff23cc37bb2a588591644cbf6e11ba1221cf9260",
                   "904b332bbbbf1dec6e792cbaa1fb8994eb3bbcfec298f0079108ec038a84331e"),
}


def _network_class():
    from rmgpy.polymer_spec import NetworkPolymerSpec
    return NetworkPolymerSpec


# --------------------------------------------------------------------------
# the declaration
# --------------------------------------------------------------------------
def test_the_network_fixture_reads_through_the_deck_reader():
    from rmgpy.polymer_spec import BranchedPolymerSpec
    deck = kr.read_deck(NETWORK_DECK)
    assert type(deck.spec) is _network_class()
    assert isinstance(deck.spec, BranchedPolymerSpec)
    assert deck.spec.architecture == "network"


def test_a_network_declaration_materializes():
    """Refused by name before S7; now one giant component plus sol."""
    from rmgpy.kmc.melt_state import materialize
    state = materialize(network_spec())
    assert len(state.components) >= 1


def test_the_normalized_json_is_a_fixed_point():
    from rmgpy.polymer_spec import normalized_polymer_spec
    spec = network_spec()
    assert normalized_polymer_spec(spec.normalized_json).normalized_json == spec.normalized_json


def test_counts_below_the_gel_point_are_refused_as_a_network():
    """One spelling per polymer: a sub-gel network is a branched declaration."""
    with pytest.raises(ValueError, match="refused:gel_point.*declare it branched"):
        network_spec(counts(residues={"Ph": 40, "M": 20}, links={"m_o": 30, "m_p": 10}))


def test_a_branched_declaration_above_the_gel_point_is_still_refused():
    with pytest.raises(ValueError, match="refused:gel_point"):
        branched_spec(counts(residues={"Ph": 40, "M": 36}, links={"m_o": 52, "m_p": 20}))


def test_the_architecture_must_be_in_the_topology_vocabulary():
    def branched_only(mapping):
        mapping["locality_policy"]["bounds"]["topology_metadata_schema"]["architecture"] = ["branched"]
    with pytest.raises(ValueError, match="architecture 'network' is not in"):
        network_spec(branched_only)


def test_inert_phases_parse_as_exact_rationals():
    spec = network_spec()
    phases = dict((p.name, p) for p in spec.inert_phases)
    assert phases["filler"].mass_fraction == Fraction(1, 2)
    assert dict(phases["filler"].elements) == {"C": 1}
    assert network_spec(no_filler).inert_phases == ()


@pytest.mark.parametrize("bad", [
    0.5,
    {"numerator": 1, "denominator": 1},
    {"numerator": 0, "denominator": 3},
    {"numerator": -1, "denominator": 3},
])
def test_an_inert_mass_fraction_outside_zero_one_or_inexact_is_refused(bad):
    def mutate(mapping):
        mapping["inert_phases"]["filler"]["mass_fraction"] = bad
    with pytest.raises(ValueError, match="inert"):
        network_spec(mutate)


def test_inert_mass_fractions_summing_to_one_are_refused():
    def two(mapping):
        mapping["inert_phases"]["glass"] = {"mass_fraction": {"numerator": 1, "denominator": 2},
                                            "elements": {"Si": 1, "O": 2}}
    with pytest.raises(ValueError, match="inert_phases.*sum"):
        network_spec(two)


def test_an_inert_phase_of_an_unknown_element_is_refused():
    def unknown(mapping):
        mapping["inert_phases"]["filler"]["elements"] = {"Xx": 1}
    with pytest.raises(ValueError, match="inert"):
        network_spec(unknown)


def test_a_branched_declaration_cannot_carry_inert_phases():
    from branched_support import branched_mapping
    from rmgpy.polymer_spec import normalized_polymer_spec
    mapping = branched_mapping()
    mapping["inert_phases"] = network_mapping()["inert_phases"]
    with pytest.raises(ValueError, match="unknown BranchedPolymerSpec fields"):
        normalized_polymer_spec(mapping)


def test_fillers_and_counts_stay_outside_the_chemistry_digest():
    base = network_spec()
    bare = network_spec(no_filler)
    other = network_spec(counts(residues={"Ph": 40, "M": 42}, links={"m_o": 60, "m_p": 24}, seed=11))
    assert kr.chemistry_sha256(bare) == kr.chemistry_sha256(base) == kr.chemistry_sha256(other)
    assert len({kr.declaration_sha256(s) for s in (base, bare, other)}) == 3


def test_a_network_and_a_branched_declaration_of_the_same_chemistry_differ_in_digest():
    """The architecture is chemistry: a network reads trees, a branched melt finite molecules."""
    assert kr.chemistry_sha256(network_spec(no_filler)) != kr.chemistry_sha256(branched_spec())


# --------------------------------------------------------------------------
# MoM refuses a network by name
# --------------------------------------------------------------------------
def test_polymer_refuses_a_network_declaration_by_name():
    from rmgpy.polymer import Polymer
    with pytest.raises(ValueError, match="refused:network_topology"):
        Polymer(label="phenolic", monomer=None, polymer_spec=network_spec())


def test_the_input_file_polymer_refuses_a_network_declaration_by_name():
    from rmgpy.rmg import input as rmg_input
    mapping = network_mapping()
    with pytest.raises(ValueError, match="refused:network_topology"):
        rmg_input.polymer(mapping["label"], polymer_spec=mapping)


def test_the_shared_pool_config_refuses_a_network_declaration_by_name():
    from types import SimpleNamespace
    from rmgpy.polymer_shared_support import _authoritative_spec
    with pytest.raises(ValueError, match="refused:network_topology"):
        _authoritative_spec(SimpleNamespace(label="p", linear_homopolymer_spec=network_spec()))


def test_the_melt_run_refuses_a_mom_comparison_of_a_network_deck(tmp_path):
    mom = tmp_path / "mom.json"
    mom.write_text("{}")
    with pytest.raises(kr.RunError, match="refused:network_topology"):
        kr.run(NETWORK_DECK, [], 700.0, [0.0], [1], tmp_path / "out", mom=mom)
    assert not (tmp_path / "out").exists()


def test_the_linear_melt_runner_refuses_a_network_deck_by_name(tmp_path):
    """The runner's observables read linear chains; a network run is
    ``melt_network.observed_run``."""
    with pytest.raises(kr.RunError, match="refused:network_observables"):
        kr.run(NETWORK_DECK, [], 700.0, [0.0], [1], tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_a_branched_declaration_is_still_refused_as_branched():
    from rmgpy.polymer import Polymer
    with pytest.raises(ValueError, match="refused:branched_mom"):
        Polymer(label="novolac", monomer=None, polymer_spec=branched_spec())


# --------------------------------------------------------------------------
# every earlier declaration keeps its bytes
# --------------------------------------------------------------------------
@pytest.mark.parametrize("name", sorted(PINS))
def test_earlier_deck_digests_are_byte_identical(name):
    spec = kr.read_deck(DECKS / name).spec
    assert (kr.declaration_sha256(spec), kr.chemistry_sha256(spec)) == PINS[name]


def test_the_branched_fixture_is_still_branched():
    from rmgpy.polymer_spec import BranchedPolymerSpec
    spec = kr.read_deck(BRANCHED_DECK).spec
    assert type(spec) is BranchedPolymerSpec
    assert "inert_phases" not in spec.as_mapping()
