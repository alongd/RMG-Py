"""Consequence-1 axis: does a BOND-FREE neutral argon have two candidate atom types?"""
import rmgpy.molecule.atomtype as at
from rmgpy.molecule.atomtype import ATOMTYPES, get_atomtype
from rmgpy.molecule import Group, Molecule
from rmgpy.exceptions import AtomTypeError

print("Ar0s single=%s" % ATOMTYPES['Ar0s'].single)

# Every bond-free neutral argon the perceiver can be handed: p3 with u0/u1/u2, and p4 with u0.
for u in (0, 1, 2):
    adj = "1 Ar u%d p3 c0" % u
    try:
        m = Molecule().from_adjacency_list(adj, raise_charge_exception=False)
        m.update_atomtypes()
        print("  %-16s -> %s" % (adj, m.atoms[0].atomtype.label))
    except AtomTypeError as e:
        print("  %-16s -> AtomTypeError: %s" % (adj, str(e).splitlines()[0]))
    except Exception as e:
        print("  %-16s -> %s: %s" % (adj, type(e).__name__, str(e).splitlines()[0]))

adj = "1 Ar u0 p4 c0"
m = Molecule().from_adjacency_list(adj)
m.update_atomtypes()
print("  %-16s -> %s" % (adj, m.atoms[0].atomtype.label))

# Group-matching axis: does the bond-free Ar0s GROUP still exist as a matchable pattern?
for pat in ("1 Ar0s ux", "1 Ar0 ux"):
    g = Group().from_adjacency_list(pat)
    for mol_adj in ("1 Ar u0 p4 c0",):
        mm = Molecule().from_adjacency_list(mol_adj)
        mm.update_atomtypes()
        print("  '%s' matches '%s' : %s" % (mol_adj, pat, mm.is_subgraph_isomorphic(g)))

# Ar2+ both halves, re-asserted on this build
ar2p = Molecule().from_adjacency_list("1 Ar u1 p3 c0 {2,S}\n2 Ar u0 p3 c+1 {1,S}\n")
ar2p.update_atomtypes()
print("  Ar2+ types =", [a.atomtype.label for a in ar2p.atoms])
