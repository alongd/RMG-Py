#!/usr/bin/env python
"""Round-83 HIGH 2 probe: why the skeleton key cannot be refined to separate tautomers.

The wall recycle key is the standard InChI truncated before its charge (/q) and proton
(/p) layers, so a cation and its neutral share a key. Standard InChI DELIBERATELY merges
tautomers (mobile-H), so 2-pyridone and 2-hydroxypyridine key identically -- an energy
pick among them transmutes one into the other (round-79 HIGH-2 defect).

The obvious refinement -- generate a FixedH InChI, whose /f layer fixes the H positions
and so separates tautomers -- does NOT compose with the charge-truncation the key needs:
for a cation, InChI places /p and /f BEFORE /q and shifts the main formula, so truncating
at /q,/p to make the key charge-independent strips the /f layer that carries the tautomer
distinction. This probe prints the layer orders that show it. Conclusion: no finer key;
the fix is the multiplicity-requires-declaration rule (Option 1). Run with rmg_env on PATH,
PYTHONPATH=worktree, MPLCONFIGDIR set.
"""
from rmgpy.molecule import Molecule


def skel(inchi):
    cut = len(inchi)
    for sep in ('/q', '/p'):
        i = inchi.find(sep)
        if i != -1 and i < cut:
            cut = i
    return inchi[:cut]


def std(m):
    return m.to_inchi()


def fixedh(m):
    from rdkit import Chem
    return Chem.MolToInchi(m.to_rdkit_mol(remove_h=False, return_mapping=False),
                           options='/FixedH')


print("== standard InChI (the key today): tautomers MERGE, isomers SPLIT, Ar states UNIFY ==")
cases = {
    'Ar   (u0 p4)     ': Molecule().from_adjacency_list('1 Ar u0 p4 c0'),
    'Ar*  (u2 p3)     ': Molecule().from_adjacency_list('multiplicity 3\n1 Ar u2 p3 c0'),
    '2-pyridone       ': Molecule().from_smiles('O=c1cccc[nH]1'),
    '2-hydroxypyridine': Molecule().from_smiles('Oc1ccccn1'),
    'DME  (COC)       ': Molecule().from_smiles('COC'),
    'ethanol (CCO)    ': Molecule().from_smiles('CCO'),
}
for name, m in cases.items():
    print("  {0}  skel = {1}".format(name, skel(std(m))))
print("  -> pyridone and hydroxypyridine share a key (tautomer merge); DME != ethanol;"
      " Ar == Ar*")

print("\n== FixedH separates the tautomer NEUTRALS but its /f sits BEHIND /q,/p on IONS ==")
print("  neutral 2-pyridone       FixedH:", fixedh(Molecule().from_smiles('O=c1cccc[nH]1')))
print("  neutral 2-hydroxypyridine FixedH:", fixedh(Molecule().from_smiles('Oc1ccccn1')))
print("  cation  pyridinium (proton) FixedH:", fixedh(Molecule().from_smiles('O=c1cccc[nH+]1')))
print("  cation  formamide+ (radical) FixedH:", fixedh(Molecule().from_smiles('N[CH][O+]')))
print("  -> on the cations /p and /f precede /q and the main formula shifts, so truncating")
print("     at /q,/p to make the key charge-independent DROPS /f. FixedH cannot be the")
print("     ion-side key. No fourth projection: the fix is multiplicity -> declaration.")
