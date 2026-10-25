"""Declaration schema v2, ``architecture: branched`` (S6): residues with
vacant-capable sites, the bridge its own residue, the declared minimum cycle
length and the ``random_bond_v1`` initial state.

The branched novolac test fixture parses (it is refused by architecture
before S6); each validation rule refuses by name what cannot be built; the
initial state stays outside the chemistry digest; MoM refuses a branched
declaration by name; v1 and linear v2 declarations keep their bytes."""

import pytest

from rmgpy.kmc import run as kr
from rmgpy.polymer_spec import LinearHomopolymerSpec, PolymerSpec

from branched_support import BRANCHED_DECK, branched_mapping, branched_spec, counts
from polymer_spec_v2_support import DECKS

# declaration / chemistry sha256 of the linear v2 decks at the S6 base (fa9f58f21)
LINEAR_V2_PINS = {
    "ep_v2_deck.py": ("ddb53e4b1edd219abf5182d9502f8b84fb49df59997337226ed96e817bab550c",
                      "4dbceba4482bdeeffbf45acfe0de6d447538df7745218e89b10465c86559b51b"),
    "epdm_v2_deck.py": ("33664ce5b12206ed9f61d8d9ff1e0abc3728296840b0a07b10fbe32007a0a9f4",
                        "a0603a35c4e7abd2bf251f2fe2991b4e6a7e3a206c97392c77ca941887f23dc2"),
}


def _branched_class():
    from rmgpy.polymer_spec import BranchedPolymerSpec
    return BranchedPolymerSpec


# --------------------------------------------------------------------------
# the declaration
# --------------------------------------------------------------------------
def test_the_branched_fixture_reads_through_the_deck_reader():
    deck = kr.read_deck(BRANCHED_DECK)
    assert type(deck.spec) is _branched_class()
    assert deck.spec.architecture == "branched"
    assert set(deck.digests()) == {"declaration_sha256", "chemistry_sha256"}


def test_normalized_json_is_a_fixed_point():
    from rmgpy.polymer_spec import normalized_polymer_spec
    spec = branched_spec()
    assert normalized_polymer_spec(spec.normalized_json).normalized_json == spec.normalized_json


def test_the_declared_minimum_cycle_and_its_default():
    """Absent, ``minimum_cycle_bonds`` is 2(r+1)+2 (8 at radius 2) and is
    written out, so stating the default changes no digest."""
    def drop(mapping):
        del mapping["minimum_cycle_bonds"]
    assert branched_spec().minimum_cycle_bonds == 8
    assert branched_spec(drop).normalized_json == branched_spec().normalized_json


def test_the_initial_state_stays_outside_the_chemistry_digest():
    base = kr.chemistry_sha256(branched_spec())
    other = branched_spec(counts(residues={"Ph": 30, "M": 12}, links={"m_o": 18, "m_p": 6}, seed=11))
    assert kr.chemistry_sha256(other) == base
    assert kr.declaration_sha256(other) != kr.declaration_sha256(branched_spec())


def test_the_links_and_the_minimum_cycle_enter_the_chemistry_digest():
    def para_only(mapping):
        del mapping["links"]["m_o"]
        del mapping["initial_ensemble"]["link_counts"]["m_o"]
        mapping["residues"]["Ph"]["sites"] = {"p4": mapping["residues"]["Ph"]["sites"]["p4"]}
        mapping["locality_policy"]["bounds"]["attachment_classes"] = ["methylene", "para"]
        mapping["initial_ensemble"]["representative_residue_counts"] = {"Ph": 20, "M": 5}
        mapping["initial_ensemble"]["link_counts"] = {"m_p": 10}

    def longer(mapping):
        mapping["minimum_cycle_bonds"] = 10
    base = kr.chemistry_sha256(branched_spec())
    assert kr.chemistry_sha256(branched_spec(longer)) != base
    assert kr.chemistry_sha256(branched_spec(para_only)) != base


# --------------------------------------------------------------------------
# refused by name
# --------------------------------------------------------------------------
def test_a_network_is_still_refused_by_name():
    def network(mapping):
        mapping["architecture"] = "network"
        mapping["locality_policy"]["bounds"]["topology_metadata_schema"]["architecture"] = ["network"]
    with pytest.raises(ValueError, match="refused:architecture: network"):
        branched_spec(network)


def test_a_minimum_cycle_below_the_tree_bound_is_refused():
    """Below 2(r+1)+2 bonds a cycle fits inside one environment: the
    environment set is no longer a tree set."""
    def short(mapping):
        mapping["minimum_cycle_bonds"] = 7
    with pytest.raises(ValueError, match="refused:minimum_cycle_bonds"):
        branched_spec(short)


def test_a_declaration_without_a_vacant_capable_site_is_linear():
    def forbidden(mapping):
        for site in mapping["residues"]["Ph"]["sites"].values():
            site["vacancy"] = "forbidden"
    with pytest.raises(ValueError, match="declare it as linear"):
        branched_spec(forbidden)


def test_one_class_both_forbidden_and_vacant_is_refused():
    def mixed(mapping):
        mapping["residues"]["Ph"]["sites"]["o6"]["vacancy"] = "forbidden"
    with pytest.raises(ValueError, match="refused:mixed_vacancy"):
        branched_spec(mixed)


def test_a_double_bond_site_is_refused():
    def double(mapping):
        mapping["residues"]["M"]["sites"]["a"]["bond_order"] = 2
        mapping["residues"]["M"]["graph"]["atoms"][0]["hydrogens"] = 1
    with pytest.raises(ValueError, match="refused:site_order"):
        branched_spec(double)


def test_link_counts_that_leave_a_bridge_site_vacant_are_refused():
    """Every methylene site is forbidden-vacant: 2 N_M link endpoints."""
    with pytest.raises(ValueError, match="refused:link_counts"):
        branched_spec(counts(links={"m_o": 30, "m_p": 9}))


def test_link_counts_beyond_the_sites_of_a_class_are_refused():
    with pytest.raises(ValueError, match="refused:link_counts"):
        branched_spec(counts(residues={"Ph": 10, "M": 20}, links={"m_o": 30, "m_p": 10}))


def test_counts_at_or_above_the_gel_point_are_refused():
    """Equal reactivity: alpha = 1/2 is Flory's gel point for three sites
    per ring; every site occupied is well above it."""
    with pytest.raises(ValueError, match="refused:gel_point"):
        branched_spec(counts(residues={"Ph": 40, "M": 60}, links={"m_o": 80, "m_p": 40}))
    with pytest.raises(ValueError, match="refused:gel_point"):
        branched_spec(counts(residues={"Ph": 40, "M": 30}, links={"m_o": 40, "m_p": 20}))
    branched_spec(counts(residues={"Ph": 40, "M": 29}, links={"m_o": 39, "m_p": 19}))


def test_a_link_inside_one_pi_system_is_refused():
    """Rule 2b on the residue pair each link joins."""
    def vinyl_bridge(mapping):
        mapping["residues"]["M"]["graph"]["atoms"][0]["hydrogens"] = 0
        mapping["residues"]["M"]["graph"]["atoms"].append(
            {"key": "n", "element": "C", "hydrogens": 2, "electronic_state": {"radical_electrons": 0}})
        mapping["residues"]["M"]["graph"]["bonds"] = [{"atoms": ["m", "n"], "order": 2}]
    with pytest.raises(ValueError, match="rule 2b"):
        branched_spec(vinyl_bridge)


# --------------------------------------------------------------------------
# MoM refuses a branched declaration by name
# --------------------------------------------------------------------------
def test_polymer_refuses_a_branched_declaration_by_name():
    from rmgpy.polymer import Polymer
    with pytest.raises(ValueError, match="refused:branched_mom"):
        Polymer(label="novolac", monomer=None, polymer_spec=branched_spec())


def test_the_input_file_polymer_refuses_a_branched_declaration_by_name():
    from rmgpy.rmg import input as rmg_input
    mapping = branched_mapping()
    with pytest.raises(ValueError, match="refused:branched_mom"):
        rmg_input.polymer(mapping["label"], polymer_spec=mapping)


def test_the_melt_run_refuses_a_mom_comparison_of_a_branched_deck(tmp_path):
    mom = tmp_path / "mom.json"
    mom.write_text("{}")
    with pytest.raises(kr.RunError, match="refused:branched_mom"):
        kr.run(BRANCHED_DECK, [], 700.0, [0.0], [1], tmp_path / "out", mom=mom)
    assert not (tmp_path / "out").exists()


def test_the_melt_run_of_a_branched_deck_is_refused_by_name(tmp_path):
    """The melt observables read linear chains (S6 builds the initial state
    and the declared linearity rule, not the observables)."""
    with pytest.raises(kr.RunError, match="refused:branched_observables"):
        kr.run(BRANCHED_DECK, [], 700.0, [0.0], [1], tmp_path / "out")
    assert not (tmp_path / "out").exists()


# --------------------------------------------------------------------------
# v1 and linear v2 keep their bytes
# --------------------------------------------------------------------------
@pytest.mark.parametrize("name", sorted(LINEAR_V2_PINS))
def test_linear_v2_deck_digests_are_unchanged(name):
    deck = kr.read_deck(DECKS / name)
    assert type(deck.spec) is PolymerSpec
    declaration, chemistry = LINEAR_V2_PINS[name]
    assert deck.digests() == {"declaration_sha256": declaration, "chemistry_sha256": chemistry}


@pytest.mark.parametrize("name", ["pe_toy_deck.py", "pp_deck.py", "epdm_deck.py", "novolac_deck.py"])
def test_v1_decks_are_still_v1(name):
    assert type(kr.read_deck(DECKS / name).spec) is LinearHomopolymerSpec
