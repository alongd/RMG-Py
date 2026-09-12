"""Probe the Ar0s sample-molecule path and the three controls. Asserts on values, never on absence."""
import logging
import rmgpy.molecule.atomtype as at
from rmgpy.molecule.atomtype import ATOMTYPES, get_atomtype
from rmgpy.molecule import Group, Molecule, Atom, Bond
from rmgpy.exceptions import UnexpectedChargeError

print("MODULE :", at.__file__)
print("Ar0s single=%s lone_pairs=%s charge=%s" % (
    ATOMTYPES['Ar0s'].single, ATOMTYPES['Ar0s'].lone_pairs, ATOMTYPES['Ar0s'].charge))
print()

# ---------- the sample-molecule path, step by step ----------
print("=== make_sample_molecule('1 Ar0s ux') ===")
g = Group().from_adjacency_list("1 Ar0s ux")
mg = g.copy(deep=True)
mg.pick_wildcards()
print("after pick_wildcards, radical_electrons =", mg.vertices[0].radical_electrons)
mg2 = mg.add_implicit_atoms_from_atomtype()
print("atoms after add_implicit_atoms_from_atomtype =", len(mg2.vertices))
sample_atom = mg2.vertices[0].make_sample_atom()
print("make_sample_atom -> symbol=%s u=%s p=%s c=%s" % (
    sample_atom.symbol, sample_atom.radical_electrons, sample_atom.lone_pairs, sample_atom.charge))
m = Molecule(atoms=[sample_atom])
sample_atom.update_charge()
print("after update_charge with 0 bonds -> charge =", sample_atom.charge)
try:
    result = g.make_sample_molecule()
    print("make_sample_molecule OK:\n" + result.to_adjacency_list())
except UnexpectedChargeError as e:
    print("make_sample_molecule RAISED UnexpectedChargeError")
    print("offending graph:\n" + e.graph.to_adjacency_list())
print()

# ---------- Control A: Ar2+ ----------
print("=== Control A: Ar2+ ===")
ar2p = Molecule().from_adjacency_list("""
1 Ar u1 p3 c0 {2,S}
2 Ar u0 p3 c+1 {1,S}
""")
ar2p.update_atomtypes()
types_a = [a.atomtype.label for a in ar2p.atoms]
print("Ar2+ atom types =", types_a)
print("Ar2+ net charge =", ar2p.get_net_charge())
print("Ar2+ SMILES-ish multiplicity =", ar2p.multiplicity)
assert types_a == ['Ar0s', 'Ar+'], types_a
print("CONTROL A PASS")
print()

# ---------- Control B: bare argon ----------
print("=== Control B: bare argon atom ===")
bare = Molecule().from_adjacency_list("1 Ar u0 p4 c0")
bare.update_atomtypes()
t_b = bare.atoms[0].atomtype.label
print("bare argon type =", t_b)
assert t_b == 'Ar0', t_b
# and check ambiguity directly through get_atomtype
direct = get_atomtype(bare.atoms[0], {})
print("get_atomtype(bare Ar) =", direct.label)
assert direct.label == 'Ar0', direct.label
print("CONTROL B PASS")
print()

# ---------- Control C: Ar+ / Ar++ ----------
print("=== Control C: Ar+ / Ar++ ===")
arp = Molecule().from_adjacency_list("1 Ar u1 p3 c+1")
arp.update_atomtypes()
print("Ar u1 p3 c+1 ->", arp.atoms[0].atomtype.label)
assert arp.atoms[0].atomtype.label == 'Ar+', arp.atoms[0].atomtype.label
arpp = Molecule().from_adjacency_list("1 Ar u0 p3 c+2")
arpp.update_atomtypes()
print("Ar u0 p3 c+2 ->", arpp.atoms[0].atomtype.label)
assert arpp.atoms[0].atomtype.label == 'Ar++', arpp.atoms[0].atomtype.label
print("CONTROL C PASS")
print()

# ---------- Group-matching surface ----------
print("=== Group matching surface ===")
for lbl in ['R', 'R!H', 'R!H!Val7', 'Rx', 'Rx!H', 'Ar']:
    print("  Ar0s is_specific_case_of(%-9s) = %s" % (
        lbl, ATOMTYPES['Ar0s'].is_specific_case_of(ATOMTYPES[lbl])))
# the Ar0s GroupAtom vs the bonded Ar0s atom in Ar2+
ga = Group().from_adjacency_list("1 Ar0s ux")
print("  Ar2+ subgraph-isomorphic to '1 Ar0s ux' :", ar2p.is_subgraph_isomorphic(ga))
gr = Group().from_adjacency_list("1 R ux")
print("  Ar2+ subgraph-isomorphic to '1 R ux'    :", ar2p.is_subgraph_isomorphic(gr))
g0 = Group().from_adjacency_list("1 Ar0 ux")
print("  bare Ar subgraph-isomorphic to '1 Ar0 ux' :", bare.is_subgraph_isomorphic(g0))
print("  bare Ar subgraph-isomorphic to '1 Ar0s ux':", bare.is_subgraph_isomorphic(ga))
print()
print("ALL CONTROLS DONE")
