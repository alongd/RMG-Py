"""I-222 baseline measurement: reproduce the brief's block before touching anything."""
import rmgpy.molecule.atomtype as at
import rmgpy.molecule.molecule as mol
from rmgpy.molecule import Molecule
from rmgpy.molecule.group import Group
from rmgpy.species import Species
from rmgpy.exceptions import AtomTypeError
from rmgpy import settings

print("database.directory =", settings["database.directory"])
print("loaded atomtype module __file__ =", at.__file__)
print("loaded molecule module __file__ =", mol.__file__)
print()

print("Ar0s.single       =", at.ATOMTYPES["Ar0s"].single)
print("Ar0s.lone_pairs   =", at.ATOMTYPES["Ar0s"].lone_pairs)
print("Ar0.lone_pairs    =", at.ATOMTYPES["Ar0"].lone_pairs)
print("existing Ar types =", [n for n in at.ATOMTYPES if n.startswith("Ar")])
print()

ADJ = "1 Ar u2 p3 c0"

try:
    m = Molecule().from_adjacency_list(ADJ)
    print(f"Molecule().from_adjacency_list({ADJ!r}) -> OK, atomtype = {m.atoms[0].atomtype.label}")
except AtomTypeError as e:
    print(f"Molecule().from_adjacency_list({ADJ!r}) -> RAISES AtomTypeError")
    print("   ", " ".join(str(e).split())[:400])

try:
    s = Species().from_adjacency_list(ADJ)
    print(f"Species().from_adjacency_list({ADJ!r}) -> OK, atomtype = "
          f"{s.molecule[0].atoms[0].atomtype.label}")
except AtomTypeError as e:
    print(f"Species().from_adjacency_list({ADJ!r}) -> RAISES AtomTypeError")
    print("   ", " ".join(str(e).split())[:400])
print()

# tolerant path
m = Molecule().from_adjacency_list(ADJ, raise_atomtype_exception=False)
m.update_atomtypes(log_species=False, raise_exception=False)
label = m.atoms[0].atomtype.label
print("tolerant update_atomtypes(log_species=False, raise_exception=False) -> atomtype =", label)

grp = Group().from_adjacency_list("1 *1 R u[1,2,3,4] px c[0,+1,+2,+3,+4]")
print('  and that argon matches "1 *1 R u[1,2,3,4] px c[0,+1,...]" ->',
      m.is_subgraph_isomorphic(grp))
