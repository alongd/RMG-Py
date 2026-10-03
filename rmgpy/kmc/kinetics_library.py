"""One owner-selected PLP-SEC rate, scoped to an exact PS propagation rewrite."""

from __future__ import annotations

import copy
import math
from functools import lru_cache

from rmgpy import constants
from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.resonance import generate_aromatic_resonance_structure


_LIBRARY_ENTRY = {
    "library": "styrene_plpsec",
    "entry_id": "styrene_head_to_tail_propagation",
    "reaction": "R-CH2-CH*(Ph) + CH2=CHPh -> R-CH2-CH(Ph)-CH2-CH*(Ph)",
    "citation": {
        "authors": "Buback et al.",
        "title": "Critically evaluated rate coefficients for free-radical polymerization, 1. Propagation rate coefficient for styrene",
        "journal": "Macromolecular Chemistry and Physics",
        "year": 1995,
        "volume": 196,
        "pages": "3267-3280",
        "doi": "10.1002/macp.1995.021961016",
        "url": "https://doi.org/10.1002/macp.1995.021961016",
    },
    "method": "IUPAC benchmark PLP-SEC, bulk styrene, low conversion, ambient pressure",
    "measured_temperature_range_K": [261.15, 366.15],
    "arrhenius": {"A": 4.27e7, "A_units": "L/(mol*s)", "n": 0.0,
                  "Ea": 32.5, "Ea_units": "kJ/mol"},
    "output_units": "m^3/(mol*s)",
    "extrapolation_note": "Use at 600-800 K is an Arrhenius extrapolation beyond the measured range; not a high-temperature measurement.",
    "matching_scope": "Exact ordinary H-capped alternating PS chain (at least two repeat units), terminal CH*(Ph) attacking unsubstituted styrene CH2, retaining terminal CH*(Ph).",
    "sensitivity_alternative": "Displaced RMG family estimate is recorded per matched record; disable with use_plpsec_library=False or RMG_KMC_PLPSEC_LIBRARY=0.",
}


def load_plpsec_entry():
    """Return independent data so provenance cannot mutate the library."""
    return copy.deepcopy(_LIBRARY_ENTRY)


def plpsec_rate_table(temperatures):
    """Evaluate the measured Arrhenius expression in the compiler's SI units."""
    temperatures = list(temperatures)
    arrhenius = _LIBRARY_ENTRY["arrhenius"]
    return {
        "T": temperatures,
        "k": [arrhenius["A"] * 1e-3 * math.exp(-arrhenius["Ea"] * 1000 / (constants.R * t))
              for t in temperatures],
        "interpolation": "linear-ln-k",
        "extrapolation": "refuse",
    }


@lru_cache(maxsize=32)
def _reference(smiles):
    molecule = Molecule(smiles=smiles)
    aromatic = generate_aromatic_resonance_structure(molecule, copy=True, save_order=True)
    return aromatic[0] if aromatic else molecule


def _operation_key(operation):
    action = operation["action"]
    if action in {"break", "form"}:
        return action, tuple(sorted(operation["atoms"])), float(operation["order"])
    return action, operation.get("atom"), operation.get("value")


def matches_head_to_tail(record):
    """Check full graphs AND the five mapped edits, independent of templates.

    Serialized adjacency atom order is the compiler's rewrite index order.
    Isomorphism establishes the intact alternating backbone and phenyl groups;
    the edit check establishes which styrene carbon is attacked and excludes
    any accompanying chemistry. Reverse records fail the arity check.
    """
    data = vars(record) if not isinstance(record, dict) else record
    if (data["family"] != "R_Addition_MultipleBond" or data["arity"] != 2
            or len(data["reactant_graphs"]) != 2 or len(data["product_graphs"]) != 1):
        return False
    reactants = [Molecule().from_adjacency_list(graph) for graph in data["reactant_graphs"]]
    monomer_index = next((index for index, molecule in enumerate(reactants)
                          if molecule.is_isomorphic(_reference("C=Cc1ccccc1"), save_order=True)), None)
    if monomer_index is None:
        return False
    chain = reactants[1 - monomer_index]
    carbons = chain.get_element_count().get("C", 0)
    units, remainder = divmod(carbons, 8)
    if remainder or units < 2:
        return False
    chain_smiles = "CC(c1ccccc1)" * (units - 1) + "C[CH](c1ccccc1)"
    if not chain.is_isomorphic(_reference(chain_smiles), save_order=True):
        return False
    product = Molecule().from_adjacency_list(data["product_graphs"][0])
    if not product.is_isomorphic(_reference("CC(c1ccccc1)" * units + "C[CH](c1ccccc1)"), save_order=True):
        return False
    radical = next(atom for atom in chain.atoms if atom.radical_electrons)
    monomer = reactants[monomer_index]
    alkene = next(bond for bond in monomer.get_all_edges() if bond.is_double())
    tail = next(atom for atom in (alkene.atom1, alkene.atom2)
                if sum(neighbor.element.number == 1 for neighbor in atom.edges) == 2)
    head = alkene.atom2 if tail is alkene.atom1 else alkene.atom1
    atoms = [atom for molecule in reactants for atom in molecule.atoms]
    indices = {atom: index for index, atom in enumerate(atoms)}
    radical, tail, head = (indices[atom] for atom in (radical, tail, head))
    alkene_pair = tuple(sorted((tail, head)))
    expected = {
        ("break", alkene_pair, 2.0),
        ("form", alkene_pair, 1.0),
        ("form", tuple(sorted((radical, tail))), 1.0),
        ("set_radical", radical, 0),
        ("set_radical", head, 1),
    }
    return (len(data["bond_ops"]) == len(expected)
            and {_operation_key(operation) for operation in data["bond_ops"]} == expected)
