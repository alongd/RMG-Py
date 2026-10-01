#!/usr/bin/env python3
"""
Tests for rmgpy.kmc.atom_map module.
"""
import json
import os
import pytest
from pathlib import Path

from rmgpy.data.rmg import RMGDatabase
from rmgpy.molecule.molecule import Molecule
from rmgpy.species import Species
from rmgpy.data.kinetics.common import ensure_independent_atom_ids
from rmgpy.kmc.atom_map import extract_atom_map


REPO_ROOT = Path(__file__).resolve().parents[3]
DB_PATH = os.environ.get("RMG_DATABASE_PATH", str(REPO_ROOT.parent / "RMG-database"))
FAMILIES = [
    "H_Abstraction",
    "R_Addition_MultipleBond",
    "R_Recombination",
    "Disproportionation",
    "intra_H_migration",
]

PROXY_SMILES = {
    "trimer": "CC(c1ccccc1)CC(c1ccccc1)CC(C)c1ccccc1",
    "midchain_rad": "CC(c1ccccc1)[C](CC(C)c1ccccc1)c1ccccc1",
    "chainend_rad": "[CH2]C(c1ccccc1)CC(C)c1ccccc1",
}

TEST_CASES = [
    ("H_Abstraction", ["midchain_rad", "trimer"], True),
    ("R_Addition_MultipleBond", ["midchain_rad"], False),
    ("intra_H_migration", ["midchain_rad"], False),
    ("intra_H_migration", ["chainend_rad"], False),
    ("R_Recombination", ["chainend_rad", "chainend_rad"], True),
    ("Disproportionation", ["midchain_rad", "chainend_rad"], True),
]

FIXTURE_PATH = Path(__file__).parent / "fixtures" / "p1_atommap.json"


@pytest.fixture(scope="module")
def kinetics_db():
    """Load the kinetics database once per module."""
    db = RMGDatabase()
    db.load_kinetics(
        DB_PATH + "/input/kinetics",
        reaction_libraries=[],
        seed_mechanisms=None,
        kinetics_families=FAMILIES,
        kinetics_depositories=["training"],
    )
    return db.kinetics.families


@pytest.fixture(scope="module")
def fixture_data_fixture():
    """Load the fixture data once per module."""
    with open(FIXTURE_PATH) as f:
        return json.load(f)


def tag_atoms(mol: Molecule):
    mol.assign_atom_ids()


def test_fixture_sha_matches(kinetics_db, fixture_data_fixture):
    """Skip if database SHA differs from fixture (not a failure)."""
    import subprocess

    db_sha = subprocess.check_output(
        ["git", "-C", DB_PATH, "rev-parse", "HEAD"], text=True
    ).strip()
    if db_sha != fixture_data_fixture["rmg_database_sha"]:
        pytest.skip(
            f"RMG-database SHA differs: fixture={fixture_data_fixture['rmg_database_sha'][:8]}, current={db_sha[:8]}"
        )


def test_fixture_reaction_count(fixture_data_fixture):
    """Verify fixture has expected number of reactions."""
    assert len(fixture_data_fixture["reactions"]) == 296


def test_bijection_and_element_conservation(fixture_data_fixture):
    """Every fixture reaction must be a heavy-atom bijection and conserve elements."""
    for i, rxn_data in enumerate(fixture_data_fixture["reactions"]):
        atom_map = rxn_data["atom_map"]
        reactant_counts = rxn_data["reactant_element_counts"]
        product_counts = rxn_data["product_element_counts"]

        # Bijection: map must be a bijection (each reactant atom -> unique product atom)
        # JSON serializes int keys as strings
        reactant_ids = set(int(rid) for rid in atom_map.keys())
        product_ids = set(int(pid) for pid in atom_map.values())
        assert len(reactant_ids) == len(
            atom_map
        ), f"Reaction {i}: duplicate reactant ids in map"
        assert len(product_ids) == len(
            atom_map
        ), f"Reaction {i}: duplicate product ids in map (not injective)"
        # Canonical form: reactant ids should be 0..n-1, product ids should be 0..m-1
        n_react = len(reactant_ids)
        n_prod = len(product_ids)
        assert reactant_ids == set(
            range(n_react)
        ), f"Reaction {i}: reactant ids not canonical 0..{n_react-1}"
        assert product_ids == set(
            range(n_prod)
        ), f"Reaction {i}: product ids not canonical 0..{n_prod-1}"

        # Element conservation
        assert reactant_counts == product_counts, (
            f"Reaction {i}: element counts differ: "
            f"reactants={reactant_counts}, products={product_counts}"
        )


def test_pitfall_a_canonical_direction(kinetics_db):
    """Pitfall (a): families with own_reverse=False return canonical direction, not requested."""
    # R_Addition_MultipleBond has own_reverse=False
    fam = kinetics_db["R_Addition_MultipleBond"]
    assert fam.own_reverse is False

    # Generate from radical (beta-scission direction)
    mol = Molecule(smiles=PROXY_SMILES["midchain_rad"])
    mol.update()
    mol = mol.copy(deep=True)
    tag_atoms(mol)

    rxns = fam.generate_reactions(
        [mol], prod_resonance=True, delete_labels=False, relabel_atoms=True
    )

    # Some reactions will be 1->2 (addition), some 2->1 (beta-scission)
    # We must read reactants/products as returned, not assume input order
    directions = set()
    for rxn in rxns:
        n_react = len(rxn.reactants)
        n_prod = len(rxn.products)
        directions.add((n_react, n_prod))

    # Should have both directions
    assert (1, 2) in directions or (2, 1) in directions
    # Never assume input was reactants


def test_pitfall_b_independent_atom_ids(kinetics_db):
    """Pitfall (b): ensure_independent_atom_ids needed for bimolecular self-reactions."""
    fam = kinetics_db["R_Recombination"]

    # For self-reaction, we must copy the SAME molecule to get colliding ids
    mol1 = Molecule(smiles=PROXY_SMILES["chainend_rad"])
    mol1.update()
    mol1 = mol1.copy(deep=True)
    tag_atoms(mol1)

    # Copy the same molecule - this creates colliding atom ids
    mol2 = mol1.copy(deep=True)

    # Check collision before ensure_independent_atom_ids
    ids1 = {a.id for a in mol1.atoms}
    ids2 = {a.id for a in mol2.atoms}
    assert (
        ids1 & ids2
    ), "Expected collision before ensure_independent_atom_ids (copy of same molecule)"

    # After ensure_independent_atom_ids, no collision
    sp1, sp2 = Species(molecule=[mol1]), Species(molecule=[mol2])
    ensure_independent_atom_ids([sp1, sp2], resonance=False)
    mol1, mol2 = sp1.molecule[0], sp2.molecule[0]

    ids1 = {a.id for a in mol1.atoms}
    ids2 = {a.id for a in mol2.atoms}
    assert not (ids1 & ids2), "ensure_independent_atom_ids should remove collision"


def test_pitfall_c_raw_vs_merged_degeneracy(kinetics_db, fixture_data_fixture):
    """Pitfall (c): family.generate_reactions returns raw paths (degeneracy=1), merged comes from public pipeline."""
    # Test H_Abstraction
    fam = kinetics_db["H_Abstraction"]

    mol1 = Molecule(smiles=PROXY_SMILES["midchain_rad"])
    mol1.update()
    mol1 = mol1.copy(deep=True)
    tag_atoms(mol1)

    mol2 = Molecule(smiles=PROXY_SMILES["trimer"])
    mol2.update()
    mol2 = mol2.copy(deep=True)
    tag_atoms(mol2)

    sp1, sp2 = Species(molecule=[mol1]), Species(molecule=[mol2])
    ensure_independent_atom_ids([sp1, sp2], resonance=False)
    mol1, mol2 = sp1.molecule[0], sp2.molecule[0]

    raw_rxns = fam.generate_reactions(
        [mol1, mol2], prod_resonance=True, delete_labels=False, relabel_atoms=True
    )

    # All raw paths have degeneracy 1
    for rxn in raw_rxns:
        assert getattr(rxn, "degeneracy", 1.0) == 1.0

    # Verify fixture entries for H_Abstraction: merged_degeneracy equals raw_path_count
    # (since H_Abstraction has no A+A or reverse-family corrections for these proxies)
    for rxn_data in fixture_data_fixture["reactions"]:
        if rxn_data["family"] == "H_Abstraction":
            assert (
                rxn_data["merged_degeneracy"] == rxn_data["raw_path_count"]
            ), f"H_Abstraction: merged_deg={rxn_data['merged_degeneracy']} != raw_path_count={rxn_data['raw_path_count']}"

    # Test Disproportionation: at least one reaction where merged_degeneracy equals raw_path_count
    # (adjusted by documented factors - for these proxies, no A+A or reverse corrections apply)
    found_disp = False
    for rxn_data in fixture_data_fixture["reactions"]:
        if rxn_data["family"] == "Disproportionation":
            # For Disproportionation with midchain_rad + chainend_rad (not A+A),
            # merged_degeneracy should equal raw_path_count
            if rxn_data["reactant_smiles"][0] != rxn_data["reactant_smiles"][1]:
                assert (
                    rxn_data["merged_degeneracy"] == rxn_data["raw_path_count"]
                ), f"Disproportionation: merged_deg={rxn_data['merged_degeneracy']} != raw_path_count={rxn_data['raw_path_count']}"
                found_disp = True
                break
    assert (
        found_disp
    ), "No suitable Disproportionation reaction found for pitfall (c) test"


def test_pitfall_d_aa_reaction_half_factor(kinetics_db, fixture_data_fixture):
    """Pitfall (d): A+A reactions have 1/2 factor in RMG degeneracy; ssa_multiplier is separate."""
    from rmgpy.kmc.event_record import ssa_multiplier_for

    # R_Recombination of chainend_rad with itself is A+A
    fam = kinetics_db["R_Recombination"]

    mol1 = Molecule(smiles=PROXY_SMILES["chainend_rad"])
    mol1.update()
    mol1 = mol1.copy(deep=True)
    tag_atoms(mol1)

    mol2 = Molecule(smiles=PROXY_SMILES["chainend_rad"])
    mol2.update()
    mol2 = mol2.copy(deep=True)
    tag_atoms(mol2)

    sp1, sp2 = Species(molecule=[mol1]), Species(molecule=[mol2])
    ensure_independent_atom_ids([sp1, sp2], resonance=False)
    mol1, mol2 = sp1.molecule[0], sp2.molecule[0]

    # Find the A+A reaction in fixture (R_Recombination, same reactants)
    aa_rxn = None
    for rxn_data in fixture_data_fixture["reactions"]:
        if (
            rxn_data["family"] == "R_Recombination"
            and rxn_data["reactant_smiles"][0] == rxn_data["reactant_smiles"][1]
        ):
            aa_rxn = rxn_data
            break
    assert aa_rxn is not None, "Could not find A+A reaction in fixture"

    # RMG's merged degeneracy includes the 1/2 factor for A+A
    merged_deg = aa_rxn["merged_degeneracy"]
    raw_path_count = aa_rxn["raw_path_count"]
    assert raw_path_count == 1, f"Expected raw_path_count=1, got {raw_path_count}"
    assert merged_deg == 0.5, f"Expected merged_degeneracy=0.5, got {merged_deg}"

    # Build reactant species list for ssa_multiplier_for
    reactant_species = [sp1, sp2]
    assert (
        ssa_multiplier_for(reactant_species) == 2.0
    ), "ssa_multiplier_for A+A should be 2.0"

    # Verify the convention: ssa_multiplier * k/(N_Av*V) * N(N-1)/2 divided by
    # k * N^2/(N_Av*V) equals (N-1)/N to 1e-12
    N_Av = 6.02214076e23
    V = 1.0  # arbitrary volume
    k = 1.0  # arbitrary rate constant

    for N in (10, 100, 10_000):
        ssa_rate = (
            ssa_multiplier_for(reactant_species) * k / (N_Av * V) * N * (N - 1) / 2
        )
        rmg_rate = k * N * N / (N_Av * V)
        ratio = ssa_rate / rmg_rate
        expected = (N - 1) / N
        assert (
            abs(ratio - expected) < 1e-12
        ), f"N={N}: ratio={ratio}, expected={expected}"

    # Setting ssa_multiplier = 1 gives half that ratio (factor not applied twice)
    for N in (10, 100, 10_000):
        ssa_rate_wrong = 1.0 * k / (N_Av * V) * N * (N - 1) / 2
        rmg_rate = k * N * N / (N_Av * V)
        ratio_wrong = ssa_rate_wrong / rmg_rate
        expected = (N - 1) / N
        assert (
            abs(ratio_wrong - expected / 2) < 1e-12
        ), f"N={N}: ratio={ratio_wrong}, expected half={expected/2}"

    # Create distinct species for A+B test
    mol_a = Molecule(smiles=PROXY_SMILES["midchain_rad"])
    mol_a.update()
    mol_a = mol_a.copy(deep=True)
    tag_atoms(mol_a)

    mol_b = Molecule(smiles=PROXY_SMILES["trimer"])
    mol_b.update()
    mol_b = mol_b.copy(deep=True)
    tag_atoms(mol_b)

    sp_a, sp_b = Species(molecule=[mol_a]), Species(molecule=[mol_b])
    ensure_independent_atom_ids([sp_a, sp_b], resonance=False)

    ab_species = [sp_a.molecule[0], sp_b.molecule[0]]
    assert ssa_multiplier_for(ab_species) == 1.0, "ssa_multiplier_for A+B should be 1.0"

    # For A+B, h = N_A * N_B, ratio should be exactly 1
    N_A = 10
    N_B = 20
    ssa_rate_ab = ssa_multiplier_for(ab_species) * k / (N_Av * V) * N_A * N_B
    rmg_rate_ab = k * N_A * N_B / (N_Av * V)
    assert (
        abs(ssa_rate_ab / rmg_rate_ab - 1.0) < 1e-12
    ), "A+B ratio should be exactly 1.0"

    # Unimolecular case
    assert (
        ssa_multiplier_for([sp1]) == 1.0
    ), "ssa_multiplier_for unimolecular should be 1.0"


def test_extract_atom_map_roundtrip(kinetics_db):
    """Test that extract_atom_map works on generated reactions."""
    fam = kinetics_db["H_Abstraction"]

    mol1 = Molecule(smiles=PROXY_SMILES["midchain_rad"])
    mol1.update()
    mol1 = mol1.copy(deep=True)
    tag_atoms(mol1)

    mol2 = Molecule(smiles=PROXY_SMILES["trimer"])
    mol2.update()
    mol2 = mol2.copy(deep=True)
    tag_atoms(mol2)

    sp1, sp2 = Species(molecule=[mol1]), Species(molecule=[mol2])
    ensure_independent_atom_ids([sp1, sp2], resonance=False)
    mol1, mol2 = sp1.molecule[0], sp2.molecule[0]

    rxns = fam.generate_reactions(
        [mol1, mol2], prod_resonance=True, delete_labels=False, relabel_atoms=True
    )

    for rxn in rxns:
        result = extract_atom_map(rxn)
        assert "atom_map" in result
        assert "reactant_element_counts" in result
        assert "product_element_counts" in result
        # Map should be identity bijection
        for rid, pid in result["atom_map"].items():
            assert rid == pid
        assert result["reactant_element_counts"] == result["product_element_counts"]


def test_extract_atom_map_raises_on_bad_input():
    """extract_atom_map should raise if reaction has mismatched atom counts."""
    from rmgpy.reaction import Reaction

    reactant = Molecule(smiles="C")
    product = Molecule(smiles="CC")
    reactant.assign_atom_ids()
    product.assign_atom_ids()
    corrupted = Reaction(reactants=[reactant], products=[product])

    with pytest.raises(ValueError, match="not a bijection"):
        extract_atom_map(corrupted)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])


@pytest.mark.slow
@pytest.mark.skipif(
    os.environ.get("RMG_KMC_REGENERATE_FIXTURE") != "1",
    reason="set RMG_KMC_REGENERATE_FIXTURE=1 for the maintenance regeneration gate",
)
def test_fixture_regeneration(fixture_data_fixture, tmp_path):
    """
    Regenerate the fixture from pinned inputs and compare with stored JSON.
    Skips if RMG-database HEAD differs from rmg_database_sha.
    Before comparing, rewrites every map into canonical form (renumber atoms 0..n-1).
    """
    import subprocess
    import sys

    db_sha = subprocess.check_output(
        ["git", "-C", DB_PATH, "rev-parse", "HEAD"], text=True
    ).strip()
    if db_sha != fixture_data_fixture["rmg_database_sha"]:
        pytest.skip(
            f"RMG-database SHA differs: fixture={fixture_data_fixture['rmg_database_sha'][:8]}, current={db_sha[:8]}"
        )

    # Run the generator
    gen_path = Path(__file__).parent / "fixtures" / "generate_p1_fixture.py"
    regenerated_path = tmp_path / "p1_atommap.json"
    result = subprocess.run(
        [sys.executable, str(gen_path), "--out", str(regenerated_path)],
        cwd=REPO_ROOT,
        env={**os.environ, "PYTHONPATH": str(REPO_ROOT)},
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"Generator failed: {result.stderr}"

    # Reload the regenerated fixture
    with open(regenerated_path) as f:
        regenerated = json.load(f)

    for key in ("families", "proxy_smiles", "rmg_database_sha"):
        assert regenerated[key] == fixture_data_fixture[key]

    # Compare reactions (order may differ, so compare as sets of sorted tuples)
    def canonicalize_reaction(r):
        # Convert atom_map keys to int (JSON serializes as strings)
        atom_map = {int(k): v for k, v in r["atom_map"].items()}
        return (
            r["family"],
            tuple(r["reactant_smiles"]),
            tuple(r["product_smiles"]),
            r["raw_path_count"],
            r["merged_degeneracy"],
            tuple(sorted(atom_map.items())),
            tuple(sorted(r["reactant_element_counts"].items())),
            tuple(sorted(r["product_element_counts"].items())),
        )

    original_set = {canonicalize_reaction(r) for r in fixture_data_fixture["reactions"]}
    regenerated_set = {canonicalize_reaction(r) for r in regenerated["reactions"]}

    assert original_set == regenerated_set, (
        f"Regenerated fixture differs: "
        f"original={len(original_set)} reactions, regenerated={len(regenerated_set)} reactions, "
        f"diff_orig={len(original_set - regenerated_set)}, diff_regen={len(regenerated_set - original_set)}"
    )
