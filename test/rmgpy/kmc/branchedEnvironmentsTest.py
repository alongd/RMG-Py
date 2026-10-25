"""The environment set of a branched declaration (S6).

The declaration is lifted to a linear alphabet of residue variants
(``declared_branches.branched_lift``): every site-bearing residue gives a
vacant (H cap), an occupied-in-chain and a branch word for each site, the
branch carrying a stub arm of ceil((r + 1) / n_min) residues whose atoms are
marked and never radical sites.  Coverage follows the S1 pattern: the set
equals an oracle that enumerates every neighbourhood of every residue from
the declaration mapping alone, and contains every environment of the melt
``random_bond_v1`` builds."""

from math import ceil

import pytest

from branched_support import Oracle, branched_mapping, branched_spec, state_keys


@pytest.fixture(scope="module")
def spec():
    return branched_spec()


@pytest.fixture(scope="module")
def lift(spec):
    from rmgpy.kmc.declared_branches import branched_lift
    return branched_lift(spec)


@pytest.fixture(scope="module")
def environments(spec):
    from rmgpy.kmc.declared_branches import branched_environments
    return branched_environments(spec)


@pytest.fixture(scope="module")
def oracle():
    return Oracle(branched_mapping())


@pytest.fixture(scope="module")
def exhaustive(oracle):
    return oracle.environments(2)


# --------------------------------------------------------------------------
# the variants
# --------------------------------------------------------------------------
def test_every_site_bearing_residue_gives_vacant_in_chain_and_branch_words(lift):
    """Ph: each site class (ortho, para) is an entry or exit (in chain),
    vacant (H) and a branch in some variant; M (two forbidden sites) is only
    ever in chain.  By class: the variants are kept up to isomorphism, so the
    mirror sites o2 and o6 share their variants."""
    roles = {}

    def add(residue, site, role):
        roles.setdefault((residue, lift.sites[residue][site].site_class), set()).add(role)
    for variant in lift.variants.values():
        for site in (variant.entry, variant.exit):
            add(variant.residue, site, "in_chain")
        for site, role in variant.roles:
            add(variant.residue, site, "vacant" if role == "H" else "branch")
    assert roles[("Ph", "ortho")] == roles[("Ph", "para")] == {"in_chain", "vacant", "branch"}
    assert roles[("M", "methylene")] == {"in_chain"}


def test_a_variant_name_is_its_residue_tree(lift):
    from rmgpy.kmc.declared_branches import parse_variant
    for name, variant in lift.variants.items():
        assert parse_variant(name) == variant


def test_the_stub_arm_is_ceil_r_plus_one_over_n_min_residues_deep(lift):
    """n_min = 1 (the bridge's one backbone atom), r = 2: three residues,
    padded to the first residue that can end an arm (a ring)."""
    assert lift.stub_depth == ceil((2 + 1) / 1) == 3
    arms = [role for v in lift.variants.values() for _s, role in v.roles if role != "H"]
    assert arms
    for arm in arms:
        assert len(arm) >= lift.stub_depth
        last_residue, _entry, last_exit = arm[-1]
        assert last_exit is None and last_residue == "Ph"


def test_stub_atoms_are_marked_and_never_radical_sites(spec, lift, environments):
    from rmgpy.kmc import declared_chain as dc
    from rmgpy.kmc.declared_sequence import SequenceChainSpec, sequence_provider
    branch = next(name for name, v in lift.variants.items() if any(r != "H" for _s, r in v.roles))
    stubs = lift.stub_keys(branch)
    assert stubs and all(key.startswith("~") for key in stubs)
    assert not any(key.startswith("~") for occ in environments.values() for *_rest, key in occ)
    provider = sequence_provider(spec)
    chain = SequenceChainSpec((branch,), edits=(dc.EditOp("remove_h", 0, min(stubs)),
                                                dc.EditOp("set_radical", 0, min(stubs), 1)))
    with pytest.raises(dc.ChainSpecError, match="refused:stub_site"):
        provider.build("stub", chain)


def test_each_stub_reads_as_a_stub_one_residue_deeper(spec, lift):
    """Every branch variant's own atoms have the same environment keys in one
    context with its arm and with every arm one residue deeper."""
    from rmgpy.kmc.declared_branches import branched_lift
    deeper = branched_lift(spec, stub_depth=lift.stub_depth + 1)
    compared = 0
    for name, variant in lift.variants.items():
        if all(role == "H" for _s, role in variant.roles):
            continue
        twins = [n for n, v in deeper.variants.items()
                 if (v.residue, v.entry, v.exit) == (variant.residue, variant.entry, variant.exit)
                 and all(r == "H" if r0 == "H" else (r != "H" and r[:lift.stub_depth] == r0[:lift.stub_depth])
                         for (_s, r0), (_t, r) in zip(variant.roles, v.roles))]
        assert twins, name
        context = lift.context(name)
        base = lift.word_keys(context.word, context.position)
        for twin in twins:
            word = context.word[:context.position] + (twin,) + context.word[context.position + 1:]
            assert deeper.word_keys(word, context.position) == base, (name, twin)
            compared += 1
    assert compared


# --------------------------------------------------------------------------
# coverage (the S1 pattern)
# --------------------------------------------------------------------------
def test_the_environment_set_equals_every_neighbourhood_of_the_declaration(environments, exhaustive):
    assert set(environments) == set(exhaustive)
    assert len(environments) == 66


def test_the_oracle_is_deep_enough(oracle):
    """One residue hop misses environments two hops show (the bridge reads the
    next ring's other sites): two hops are the oracle's depth."""
    assert len(oracle.environments(1)) < len(oracle.environments(2))


def test_the_random_melt_shows_only_environments_of_the_set(spec, environments, oracle):
    from rmgpy.kmc.melt_state import materialize
    found = state_keys(materialize(spec), oracle)
    assert found and set(found) <= set(environments)


def test_chain_words_alone_miss_ends_and_short_molecules(spec, environments):
    """The interior words see no ring with one or no occupied site, and no
    bridge between two end rings: each kind of word is needed."""
    from rmgpy.kmc.declared_branches import branched_environments
    interior = branched_environments(spec, kinds=("interior",))
    with_ends = branched_environments(spec, kinds=("interior", "end"))
    assert set(interior) < set(with_ends) < set(environments)
    assert {kind for occ in environments.values() for kind, *_rest in occ} == {"interior", "end", "short"}


def test_linear_v2_environment_sets_are_unchanged():
    from rmgpy.kmc import declared_environments as de
    from rmgpy.kmc import run as kr
    from polymer_spec_v2_support import DECKS
    assert len(de.sequence_environments(kr.read_deck(DECKS / "ep_v2_deck.py").spec)) == 18
    assert len(de.sequence_environments(kr.read_deck(DECKS / "epdm_v2_deck.py").spec)) == 123
