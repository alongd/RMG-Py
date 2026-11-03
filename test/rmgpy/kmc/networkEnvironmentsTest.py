"""The environment set of a network declaration (S7): rooted residue trees
over the occupied, vacant and branch alphabet, as deep as an environment
reads (``declared_networks.tree_environments``).  The set is finite because
``minimum_cycle_bonds`` keeps every neighbourhood within reach a tree.
Coverage follows the S1 pattern: the set equals the union over seeded
random networks above the gel point, equals an oracle that enumerates every
neighbourhood from the declaration mapping alone, and equals the S6 lift's
set, so the compile's alphabet covers every network environment."""

import pytest

from branched_support import Oracle, state_keys
from network_support import counts, network_mapping, network_spec


@pytest.fixture(scope="module")
def spec():
    return network_spec()


@pytest.fixture(scope="module")
def trees(spec):
    from rmgpy.kmc.declared_networks import tree_environments
    return tree_environments(spec)


@pytest.fixture(scope="module")
def oracle():
    return Oracle(network_mapping())


def test_the_tree_depth_is_derived_from_the_radius(spec):
    """An environment reads r + 1 bonds and the heavy bond orders of its
    boundary atoms, so a site is expanded when its atom lies within r + 1
    bonds of the root residue."""
    from rmgpy.kmc.declared_networks import site_reach
    assert site_reach(spec) == int(spec.locality_policy["radius"]) + 1


def test_the_tree_set_equals_every_neighbourhood_of_the_declaration(trees, oracle):
    exhaustive = oracle.environments(2)
    assert set(trees) == set(exhaustive)
    assert all(trees[k] == exhaustive[k] for k in trees)


def test_the_tree_set_is_deterministic(spec, trees):
    from rmgpy.kmc import declared_networks as dn
    dn._TREES.clear()
    assert dn.tree_environments(spec) == trees


def test_a_shallower_tree_misses_environments(spec, trees):
    from rmgpy.kmc.declared_networks import tree_environments
    assert set(tree_environments(spec, reach=1)) < set(trees)


def test_the_union_over_random_networks_equals_the_tree_set(spec, trees, oracle):
    """Seeded random networks above the gel point, at three degrees of cure,
    show every tree environment and nothing else but one: the bridge of an
    isolated p,p'-dimer (both rings otherwise vacant), a sol molecule too rare
    at 40 rings to be drawn in 24 networks (nor in 8 sub-gel melts); it is
    that molecule's environment exactly, built as the generator builds it."""
    from rmgpy.kmc import melt_branched as mb
    from rmgpy.kmc.melt_state import materialize
    seen = set()
    plans = [({"Ph": 40, "M": 36}, {"m_o": 52, "m_p": 20}),
             ({"Ph": 40, "M": 42}, {"m_o": 60, "m_p": 24}),
             ({"Ph": 40, "M": 48}, {"m_o": 64, "m_p": 32})]
    for residues, links in plans:
        for seed in range(1, 9):
            seen |= set(state_keys(materialize(network_spec(counts(residues, links, seed))), oracle))
    assert seen <= set(trees)
    dimer = state_keys(mb.build_molecules(spec, [(["M", "Ph", "Ph"], [(0, "a", 1, "p4"), (0, "b", 2, "p4")])]),
                       oracle)
    bridge = {key for key, roots in dimer.items() if ("M", "m") in roots}
    assert set(trees) - seen == bridge and len(bridge) == 1


def test_the_lift_covers_every_network_environment(spec, trees):
    """The compile reads the S6 lift's words; a network adds no environment
    the lift lacks, so nothing is missing (``network_environment_gap``)."""
    from rmgpy.kmc.declared_branches import branched_environments
    from rmgpy.kmc.declared_networks import network_environment_gap
    assert set(trees) == set(branched_environments(spec))
    assert network_environment_gap(spec) == ()


def test_a_gap_is_refused_by_name_in_the_compile(spec, monkeypatch):
    from rmgpy.kmc import declared_networks as dn
    from rmgpy.kmc import declared_proxies as dp
    monkeypatch.setattr(dn, "network_environment_gap", lambda _spec: ("missing-key",))
    with pytest.raises(ValueError, match="refused:network_environment"):
        dp.sequence_structures(spec, 3)
