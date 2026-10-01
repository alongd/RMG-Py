#!/usr/bin/env python3
"""
Generate the p1_atommap.json fixture from the pinned RMG-Py and RMG-database SHAs.
Run from the worktree root with PYTHONPATH=$PWD.

This generator uses RMG's public pipeline logic: it generates raw reactions from each
family, then uses find_degenerate_reactions to merge isomorphic paths and compute the
correct degeneracy (including A+A 1/2 factor and reverse-family corrections).
"""
import json
import argparse
import os
import subprocess
import sys
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from rmgpy.data.rmg import RMGDatabase
from rmgpy.molecule.molecule import Molecule
from rmgpy.species import Species
from rmgpy.data.kinetics.common import ensure_independent_atom_ids, find_degenerate_reactions, check_for_same_reactants
from rmgpy.kmc.atom_map import extract_atom_map
from rmgpy.reaction import Reaction
from rmgpy.data.kinetics.family import KineticsFamily


REPO_ROOT = Path(__file__).resolve().parents[4]
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


def get_git_sha(path: str) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", path, "rev-parse", "HEAD"], text=True
        ).strip()
    except subprocess.CalledProcessError:
        return "unknown"


def tag_atoms(mol: Molecule) -> dict[int, Molecule]:
    mol.assign_atom_ids()
    return {atom.id: atom for atom in mol.atoms}


def smiles(mol: Molecule) -> str:
    return mol.to_smiles()


def canonicalize_atom_map(reactants, products, atom_map):
    """
    Renumber atoms on each side to 0..n-1 in deterministic order:
    reactant index, then atom index within the molecule as returned.
    Returns new atom_map with canonical numbering.
    """
    # Build mapping from old reactant atom id -> new canonical id
    old_to_new_reactant = {}
    new_id = 0
    for mol in reactants:
        heavy_atoms = [a for a in mol.atoms if a.element.number != 1]
        heavy_atoms.sort(key=lambda a: a.id)  # deterministic order by original id
        for atom in heavy_atoms:
            old_to_new_reactant[atom.id] = new_id
            new_id += 1

    # Build mapping from old product atom id -> new canonical id
    old_to_new_product = {}
    new_id = 0
    for mol in products:
        heavy_atoms = [a for a in mol.atoms if a.element.number != 1]
        heavy_atoms.sort(key=lambda a: a.id)
        for atom in heavy_atoms:
            old_to_new_product[atom.id] = new_id
            new_id += 1

    # Rewrite atom_map with canonical ids
    canonical_map = {}
    for old_rid, old_pid in atom_map.items():
        new_rid = old_to_new_reactant[old_rid]
        new_pid = old_to_new_product[old_pid]
        canonical_map[new_rid] = new_pid

    return canonical_map


def find_degenerate_reactions_with_groups(rxn_list, same_reactants=None, kinetics_family=None):
    """
    Wrapper around find_degenerate_reactions that also returns the grouping.
    Returns (merged_reactions, groups) where groups is list of lists of raw reactions.
    """
    # This replicates the grouping logic from find_degenerate_reactions
    selected_rxns = rxn_list

    # We want to sort all the reactions into sublists composed of isomorphic reactions
    # with degenerate transition states
    sorted_rxns = []
    for rxn0 in selected_rxns:
        rxn0.ensure_species(save_order=False)
        if len(sorted_rxns) == 0:
            # This is the first reaction, so create a new sublist
            sorted_rxns.append([rxn0])
        else:
            # Loop through each sublist, which represents a unique reaction
            for sub_list in sorted_rxns:
                # Try to determine if the current rxn0 is identical or isomorphic to any reactions in the sublist
                isomorphic = False
                identical = False
                same_template = True
                for rxn in sub_list:
                    isomorphic = rxn0.is_isomorphic(rxn, check_identical=False, strict=False,
                                                    check_template_rxn_products=True, save_order=False)
                    if isomorphic:
                        identical = rxn0.is_isomorphic(rxn, check_identical=True, strict=False,
                                                       check_template_rxn_products=True, save_order=False)
                        if identical:
                            # An exact copy of rxn0 is already in our list, so we can move on
                            break
                        same_template = frozenset(rxn.template) == frozenset(rxn0.template)
                    else:
                        # This sublist contains a different product
                        break

                # Process the reaction depending on the results of the comparisons
                if identical:
                    # This reaction does not contribute to degeneracy
                    break
                elif isomorphic:
                    if same_template:
                        # We found the right sublist, and there is no identical reaction
                        # We should add rxn0 to the sublist as a degenerate rxn, and move on to the next rxn
                        sub_list.append(rxn0)
                        break
                    else:
                        # We found an isomorphic sublist, but the reaction templates are different
                        # We need to mark this as a duplicate and continue searching the remaining sublists
                        rxn0.duplicate = True
                        sub_list[0].duplicate = True
                        continue
                else:
                    # This is not an isomorphic sublist, so we need to continue searching the remaining sublists
                    continue
            else:
                # We did not break, which means that there was no isomorphic sublist, so create a new one
                sorted_rxns.append([rxn0])

    # Now collapse each sublist into a merged reaction with proper degeneracy
    merged_reactions = []
    for sub_list in sorted_rxns:
        # Collapse our sorted reaction list by taking one reaction from each sublist
        rxn = sub_list[0]
        # The degeneracy of each reaction is the number of reactions that were in the sublist
        rxn.degeneracy = sum([reaction0.degeneracy for reaction0 in sub_list])
        merged_reactions.append(rxn)

    # Apply same_reactants reduction (A+A factor) and reverse-family correction
    for rxn in merged_reactions:
        if rxn.is_forward:
            from rmgpy.data.kinetics.common import reduce_same_reactant_degeneracy
            reduce_same_reactant_degeneracy(rxn, same_reactants)
        else:
            # fix the degeneracy of (not ownReverse) reactions found in the backwards direction
            family = kinetics_family
            if not family.own_reverse:
                rxn.degeneracy = family.calculate_degeneracy(rxn, resonance=True)

    return merged_reactions, sorted_rxns


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out",
        type=Path,
        default=REPO_ROOT / "test/rmgpy/kmc/fixtures/p1_atommap.json",
    )
    args = parser.parse_args()

    print("Loading kinetics families...")
    db = RMGDatabase()
    db.load_kinetics(
        DB_PATH + "/input/kinetics",
        reaction_libraries=None,
        seed_mechanisms=None,
        kinetics_families=FAMILIES,
        kinetics_depositories=["training"],
    )
    families = db.kinetics.families

    rmgpy_sha = get_git_sha(str(REPO_ROOT))
    db_sha = get_git_sha(DB_PATH)

    print(f"RMG-Py SHA: {rmgpy_sha}")
    print(f"RMG-database SHA: {db_sha}")

    all_merged_data = []

    for fam_label, mol_names, bimolecular in TEST_CASES:
        fam = families[fam_label]
        reactant_mols = []
        for name in mol_names:
            mol = Molecule(smiles=PROXY_SMILES[name])
            mol.update()
            reactant_mols.append(mol.copy(deep=True))

        for mol in reactant_mols:
            tag_atoms(mol)

        if bimolecular:
            species_list = [Species(molecule=[m]) for m in reactant_mols]
            ensure_independent_atom_ids(species_list, resonance=False)
            reactant_mols = [sp.molecule[0] for sp in species_list]

        # Generate raw reactions from the family
        raw_rxns = fam.generate_reactions(
            reactant_mols, prod_resonance=True, delete_labels=False, relabel_atoms=True
        )
        print(f"  {fam_label}: {len(raw_rxns)} raw reactions")

        # Check for same reactants (A+A case)
        reactants_for_check = [Species(molecule=[m]) for m in reactant_mols]
        _, same_reactants = check_for_same_reactants(reactants_for_check)

        # Merge degenerate reactions using RMG's logic (this computes correct degeneracy)
        merged_rxns, groups = find_degenerate_reactions_with_groups(
            raw_rxns, same_reactants=same_reactants, kinetics_family=fam
        )
        print(f"  {fam_label}: {len(merged_rxns)} merged reactions")

        # For each merged reaction, create fixture entry
        for i, merged in enumerate(merged_rxns):
            raw_list = groups[i]
            raw_path_count = len(raw_list)

            # Extract molecules from Species objects (find_degenerate_reactions converts to Species)
            merged_reactants = [sp.molecule[0] for sp in merged.reactants]
            merged_products = [sp.molecule[0] for sp in merged.products]

            # Create a reaction with Molecules for extract_atom_map
            merged_mol = Reaction(reactants=merged_reactants, products=merged_products)

            # Extract atom map from the merged reaction
            atom_map_result = extract_atom_map(merged_mol)

            # Canonicalize atom map
            canon_map = canonicalize_atom_map(merged_reactants, merged_products, atom_map_result["atom_map"])

            all_merged_data.append({
                "family": fam_label,
                "reactant_smiles": [smiles(m) for m in merged_reactants],
                "product_smiles": [smiles(m) for m in merged_products],
                "raw_path_count": raw_path_count,
                "merged_degeneracy": merged.degeneracy,
                "atom_map": canon_map,
                "reactant_element_counts": atom_map_result["reactant_element_counts"],
                "product_element_counts": atom_map_result["product_element_counts"],
            })

    fixture = {
        "proxy_smiles": PROXY_SMILES,
        "families": FAMILIES,
        "rmgpy_sha": rmgpy_sha,
        "rmg_database_sha": db_sha,
        "reactions": all_merged_data,
    }

    out_path = args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(fixture, f, separators=(',', ':'), sort_keys=True)

    print(f"Fixture written to {out_path}")
    print(f"Total merged reactions: {len(all_merged_data)}")

    # Print degeneracy distribution
    from collections import Counter
    deg_dist = Counter(r["merged_degeneracy"] for r in all_merged_data)
    print("Merged degeneracy distribution:")
    for deg, count in sorted(deg_dist.items()):
        print(f"  {deg}: {count}")

    # Print raw_path_count distribution
    rpc_dist = Counter(r["raw_path_count"] for r in all_merged_data)
    print("raw_path_count distribution:")
    for rpc, count in sorted(rpc_dist.items()):
        print(f"  {rpc}: {count}")

    # Check fixture size
    import os
    size_kb = os.path.getsize(out_path) / 1024
    print(f"Fixture size: {size_kb:.1f} kB")


if __name__ == "__main__":
    main()
