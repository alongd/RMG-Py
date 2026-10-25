"""Compile structures and roster of a branched declaration (S6).

The compile runs the linear-v2 machinery on the lift (residue variants): its
chain provider builds the variants' chains under ``branched:<sha16>``; every
radical structure sits on a residue's own atom, never on a stub atom; the end
and short-molecule environments get their own radicals; pairs are over the
family-tree classes; the v4 roster's participants name variants, so each
carries its residue tree.  What the compile cannot build yet is refused by
name: branch growth in the depository exclusion (no growth sites)."""

import pytest

from branched_support import branched_spec


@pytest.fixture(scope="module")
def spec():
    return branched_spec()


@pytest.fixture(scope="module")
def structures(spec):
    from rmgpy.kmc.declared_proxies import sequence_structures
    return sequence_structures(spec, 1)


def test_the_provider_is_the_lifts_under_the_declarations_digest(spec):
    import hashlib
    from rmgpy.kmc.declared_branches import BranchedChainProvider
    from rmgpy.kmc.declared_sequence import sequence_provider
    provider = sequence_provider(spec)
    assert isinstance(provider, BranchedChainProvider)
    assert provider.name == "branched:" + hashlib.sha256(spec.normalized_json).hexdigest()[:16]


def test_the_growth_units_are_the_residues_tree_increments(spec):
    """A residue added anywhere in a tree adds its formula with every site
    occupied plus (sites - 2) H: Ph C6H3O + H, M CH2."""
    from rmgpy.kmc.declared_sequence import sequence_provider
    units = sequence_provider(spec).unit_formula()
    assert units == tuple(sorted({(("C", 6), ("H", 4), ("O", 1)), (("C", 1), ("H", 2))}))


def test_branch_growth_by_construction_is_refused_by_name(spec, structures):
    """The provider names no growth sites: a depository decision that needs
    the analogues built is refused:unknown, never decided from a linear
    insertion."""
    from rmgpy.kmc.declared_sequence import sequence_provider
    build = next(s.build for name, s in structures.items() if name.startswith("pristine:"))
    assert sequence_provider(spec).growth_sites(build.molecule) is None


def test_no_radical_structure_sits_on_a_stub_atom(spec, structures):
    from rmgpy.kmc.declared_branches import branched_lift
    lift = branched_lift(spec)
    radicals = 0
    for name, s in structures.items():
        if s.witness is None:
            continue
        chain = s.witness.chain
        for edit in chain.edits:
            if 0 <= edit.unit < len(chain.residues):
                assert edit.key not in lift.stub_keys(chain.residues[edit.unit]), name
        radicals += bool(s.build.molecule.get_radical_count())
    assert radicals


def test_every_hydrogen_bearing_environment_has_a_radical(spec, structures):
    """Interior, end and short-molecule classes alike (an end ring's vacant
    ortho site, the bridge between two end rings)."""
    from rmgpy.kmc.declared_branches import branched_environments, h_bearing_environments
    covered = {s.environment for s in structures.values() if s.environment is not None}
    assert h_bearing_environments(spec, branched_environments(spec)) <= covered
    assert any(name.startswith(("head_terminal_radical:", "tail_terminal_radical:"))
               for name in structures)


def test_the_pairs_are_over_family_tree_classes(spec):
    from rmgpy.kmc.declared_proxies import PAIR_RULE, sequence_proxy_set
    proxies, record = sequence_proxy_set(spec, 1, {})
    assert record["rule"] == PAIR_RULE == "family_tree_classes"
    assert proxies and all("+" not in p.site_type for p in proxies)


def test_the_roster_participants_carry_their_residue_tree(spec, tmp_path):
    from rmgpy.kmc import emit
    from rmgpy.kmc.declared_branches import parse_variant
    roster = emit.declared_reference_systems(spec, 1, pairs={"1": [], "3": []})
    assert roster.schema == emit.ROSTER_SCHEMA_V4
    loaded = emit.write_roster(roster, tmp_path / "roster.json")
    assert loaded.digests() == roster.digests()
    built, refusals = emit.build_reference_systems(loaded)
    assert refusals == [] and set(built) == set(roster.systems)
    trees = [parse_variant(name) for system in loaded.systems.values() for entry in system
             if "declared_sequence" in entry for name in entry["declared_sequence"]["residues"]]
    assert trees and any(any(role != "H" for _s, role in tree.roles) for tree in trees)
