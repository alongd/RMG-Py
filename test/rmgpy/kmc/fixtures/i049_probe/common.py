"""Small structural helpers for the Phase 1 probes; never compile events."""

from collections import Counter
from functools import lru_cache

from rdkit import Chem

from rmgpy.kmc.state import _parse_adjacency


def benzylic_tail(units):
    """H-capped CH3–CH(Ph)–...–CH2–CH•(Ph); n=1 is ethylbenzyl."""
    assert units >= 1
    return "CC(c1ccccc1)" * (units - 1) + "C[CH](c1ccccc1)"


@lru_cache(maxsize=None)
def molecule(text):
    graph = _parse_adjacency(text)
    mol = Chem.RWMol()
    indices = {}
    for index, node in sorted(graph.items()):
        atom = Chem.Atom(node["element"])
        atom.SetNoImplicit(True)
        atom.SetNumRadicalElectrons(node["radical"])
        atom.SetFormalCharge(node["charge"])
        if any(float(order) == 1.5 for order in node["edges"].values()):
            atom.SetIsAromatic(True)
        indices[index] = mol.AddAtom(atom)
    orders = {1.0: Chem.BondType.SINGLE, 2.0: Chem.BondType.DOUBLE,
              3.0: Chem.BondType.TRIPLE, 1.5: Chem.BondType.AROMATIC}
    for index, node in graph.items():
        for other, order in node["edges"].items():
            if index < other:
                mol.AddBond(indices[index], indices[other], orders[float(order)])
    result = mol.GetMol()
    Chem.SanitizeMol(result)
    return result


@lru_cache(maxsize=None)
def describe(text):
    mol = molecule(text)
    radicals = []
    for atom in mol.GetAtoms():
        if not atom.GetNumRadicalElectrons():
            continue
        neighbors = list(atom.GetNeighbors())
        carbons = [other for other in neighbors if other.GetAtomicNum() == 6]
        aromatic = any(other.GetIsAromatic() for other in carbons)
        hydrogens = sum(other.GetAtomicNum() == 1 for other in neighbors)
        radicals.append({"index": atom.GetIdx(), "aromatic": atom.GetIsAromatic(),
                         "benzylic": aromatic, "H": hydrogens,
                         "C_neighbors": len(carbons),
                         "terminal_benzylic": not atom.GetIsAromatic() and aromatic
                         and hydrogens == 1 and len(carbons) == 2})
    return {"smiles": Chem.MolToSmiles(Chem.RemoveHs(mol)), "radicals": radicals,
            "formula": dict(Counter(a.GetSymbol() for a in mol.GetAtoms()))}


def canonical_smiles(smiles):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def terminal_benzylic(text):
    return any(radical["terminal_benzylic"] for radical in describe(text)["radicals"])
