from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from rmgpy.reaction import Reaction
    from rmgpy.molecule.molecule import Molecule
    from rmgpy.molecule.atom import Atom


def extract_atom_map(reaction: Reaction) -> dict:
    """
    Extract the heavy-atom bijection map from reactants to products for a reaction.

    Args:
        reaction: An RMG Reaction object with reactants and products that have
                  atom.id assigned (via Molecule.assign_atom_ids()).

    Returns:
        A dict with keys:
        - 'atom_map': dict mapping reactant atom id -> product atom id (heavy atoms only)
        - 'reactant_element_counts': dict element symbol -> count (incl. implicit H)
        - 'product_element_counts': dict element symbol -> count (incl. implicit H)

    Raises:
        ValueError: if the map is not a bijection or element counts differ.
    """
    reactants = reaction.reactants
    products = reaction.products

    reactant_heavy: dict[int, Atom] = {}
    for mol in reactants:
        for atom in mol.atoms:
            if atom.element.number != 1:
                if atom.id in reactant_heavy:
                    raise ValueError(f"Duplicate atom id {atom.id} in reactants")
                reactant_heavy[atom.id] = atom

    product_heavy: dict[int, Atom] = {}
    for mol in products:
        for atom in mol.atoms:
            if atom.element.number != 1:
                if atom.id in product_heavy:
                    raise ValueError(f"Duplicate atom id {atom.id} in products")
                product_heavy[atom.id] = atom

    if set(reactant_heavy.keys()) != set(product_heavy.keys()):
        missing_in_prod = set(reactant_heavy.keys()) - set(product_heavy.keys())
        missing_in_react = set(product_heavy.keys()) - set(reactant_heavy.keys())
        raise ValueError(
            f"Atom map is not a bijection: "
            f"{len(missing_in_prod)} reactant ids missing in products, "
            f"{len(missing_in_react)} product ids missing in reactants"
        )

    atom_map = {rid: pid for rid, pid in zip(sorted(reactant_heavy.keys()), sorted(product_heavy.keys()))}
    for rid in reactant_heavy:
        if atom_map[rid] != rid:
            raise ValueError(f"Atom id mismatch: reactant {rid} -> product {atom_map[rid]} (expected {rid})")

    def element_counts(mols: list[Molecule]) -> dict[str, int]:
        counts: dict[str, int] = {}
        for mol in mols:
            mol_counts = mol.get_element_count()
            for sym, cnt in mol_counts.items():
                counts[sym] = counts.get(sym, 0) + cnt
        return counts

    reactant_counts = element_counts(reactants)
    product_counts = element_counts(products)

    if reactant_counts != product_counts:
        raise ValueError(
            f"Element counts differ: reactants={reactant_counts}, products={product_counts}"
        )

    return {
        "atom_map": atom_map,
        "reactant_element_counts": reactant_counts,
        "product_element_counts": product_counts,
    }


def extract_atom_map_from_reaction(reaction: Reaction) -> dict:
    """
    Convenience wrapper that returns only the atom_map dict (reactant_id -> product_id).
    """
    return extract_atom_map(reaction)["atom_map"]
