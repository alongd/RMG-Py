"""I-222 rework probes: HIGH 1 (generic templates act on Ar0e), HIGH 3 (Ar+ LOSE_CHARGE), smalls."""
import rmgpy.molecule.atomtype as at
from rmgpy.molecule import Atom, Bond, Molecule
from rmgpy.molecule.group import Group, GroupAtom
from rmgpy.molecule.atomtype import ATOMTYPES, get_atomtype
from rmgpy.exceptions import AtomTypeError, ActionError
from rmgpy.data.kinetics.family import ReactionRecipe
from rmgpy import settings

print("database.directory =", settings["database.directory"])
print("loaded atomtype module __file__ =", at.__file__)
print()


def perceive(atom, bonds=None):
    bonds = bonds or {}
    try:
        return get_atomtype(atom, bonds).label
    except AtomTypeError:
        return "<no type>"


print("=" * 78)
print("HIGH 1 -- what generic Ar / generic R declare, and what the molecule path does")
print("=" * 78)
for label in ("Ar", "R", "Ar0e"):
    a = ATOMTYPES[label]
    print(f"  ATOMTYPES[{label!r}].set_actions:")
    for action in ("increment_bond", "decrement_bond", "form_bond", "break_bond",
                   "increment_radical", "decrement_radical", "increment_lone_pair",
                   "decrement_lone_pair", "increment_charge", "decrement_charge"):
        print(f"      {action:22s} {[t.label for t in getattr(a, action)]}")
    print()

print("-- do generic groups match a concrete Ar0e molecule? --")
ar0e = Molecule().from_adjacency_list("1 Ar u2 p3 c0")
print("   Ar0e molecule atomtype:", ar0e.atoms[0].atomtype.label)
for spec in ("1 *1 R ux px cx",
             "1 *1 R u[2,3,4] px cx",
             "1 *1 Ar ux px cx",
             "1 *1 Rx ux px cx",
             "1 *1 R!H ux px cx"):
    g = Group().from_adjacency_list(spec)
    print(f"   {spec!r:28s} -> matches Ar0e: {ar0e.is_subgraph_isomorphic(g)}")

print()
print("-- GROUP path (GroupAtom.apply_action, which DOES read set_actions) --")
for label in ("R", "Ar", "Ar0e"):
    ga = GroupAtom(atomtype=[ATOMTYPES[label]])
    try:
        ga.apply_action(["GAIN_CHARGE", "*1", 1])
        print(f"   GroupAtom([{label}]) GAIN_CHARGE -> {[t.label for t in ga.atomtype]}")
    except ActionError as e:
        print(f"   GroupAtom([{label}]) GAIN_CHARGE -> ActionError: {' '.join(str(e).split())[:110]}")

print()
print("-- MOLECULE path (ReactionRecipe on a concrete Ar0e atom) --")
for recipe_actions in ([["GAIN_CHARGE", "*1", 1]],
                       [["LOSE_RADICAL", "*1", 1]],
                       [["LOSE_RADICAL", "*1", 2], ["GAIN_PAIR", "*1", 1]]):
    mol = Molecule().from_adjacency_list("1 *1 Ar u2 p3 c0")
    try:
        ReactionRecipe(recipe_actions).apply_forward(mol, unique=True)
        mol.update_atomtypes(log_species=False, raise_exception=False)
        a = mol.atoms[0]
        print(f"   {str(recipe_actions):48s} -> u{a.radical_electrons} p{a.lone_pairs} "
              f"c{a.charge:+d} -> {a.atomtype.label}")
    except Exception as e:
        print(f"   {str(recipe_actions):48s} -> {type(e).__name__}: {' '.join(str(e).split())[:90]}")

print()
print("-- the alkaline-family-shaped recipe on two metastable argons --")
two = Molecule().from_adjacency_list("1 *1 Ar u2 p3 c0\n2 *2 Ar u2 p3 c0")
print("   before:", [a.atomtype.label for a in two.atoms])
try:
    ReactionRecipe([["LOSE_RADICAL", "*1", 1], ["LOSE_RADICAL", "*2", 1],
                    ["FORM_BOND", "*1", 1, "*2"], ["GAIN_CHARGE", "*1", 1]]
                   ).apply_forward(two, unique=True)
    two.update_atomtypes(log_species=False, raise_exception=False)
    print("   after: ", [f"u{a.radical_electrons} p{a.lone_pairs} c{a.charge:+d} "
                         f"-> {a.atomtype.label}" for a in two.atoms],
          "net charge", two.get_net_charge())
except Exception as e:
    print(f"   {type(e).__name__}: {' '.join(str(e).split())[:200]}")

print()
print("=" * 78)
print("HIGH 3 -- what LOSE_CHARGE on Ar+ actually produces")
print("=" * 78)
print("  declared Ar+.decrement_charge =", [t.label for t in ATOMTYPES["Ar+"].decrement_charge])
print("  declared Ar0.increment_charge =", [t.label for t in ATOMTYPES["Ar0"].increment_charge])
print()

# bare Ar+ : u1 p3 c+1, no bonds
bare = Atom(element="Ar", radical_electrons=1, lone_pairs=3, charge=1)
print("  bare Ar+   u1 p3 c+1            perceives as", perceive(bare))
bare.apply_action(["LOSE_CHARGE", "*1", 1])
print("    LOSE_CHARGE -> "
      f"u{bare.radical_electrons} p{bare.lone_pairs} c{bare.charge:+d} -> {perceive(bare)}")

# bonded Ar+ : the cation half of Ar2+, u0 p3 c+1 with one single bond
cat = Atom(element="Ar", radical_electrons=0, lone_pairs=3, charge=1)
neu = Atom(element="Ar", radical_electrons=1, lone_pairs=3, charge=0)
b = Bond(cat, neu, order=1)
print("  bonded Ar+ u0 p3 c+1 (1 single) perceives as", perceive(cat, {neu: b}))
cat.apply_action(["LOSE_CHARGE", "*1", 1])
print("    LOSE_CHARGE -> "
      f"u{cat.radical_electrons} p{cat.lone_pairs} c{cat.charge:+d} 1 single -> {perceive(cat, {neu: b})}")

print()
print("  inverses (GAIN_CHARGE), to see what would have to close back:")
for label, kw, bonds in (("Ar0", dict(radical_electrons=0, lone_pairs=4, charge=0), False),
                         ("Ar0e", dict(radical_electrons=2, lone_pairs=3, charge=0), False),
                         ("Ar0s u1", dict(radical_electrons=1, lone_pairs=3, charge=0), True),
                         ("Ar0s u0", dict(radical_electrons=0, lone_pairs=3, charge=0), True)):
    a = Atom(element="Ar", **kw)
    bd = {}
    if bonds:
        p = Atom(element="Ar", radical_electrons=0, lone_pairs=3, charge=1)
        bd = {p: Bond(a, p, order=1)}
    before = perceive(a, bd)
    a.apply_action(["GAIN_CHARGE", "*1", 1])
    print(f"    {label:8s} ({before:5s}) GAIN_CHARGE -> "
          f"u{a.radical_electrons} p{a.lone_pairs} c{a.charge:+d} -> {perceive(a, bd)}")

print()
print("=" * 78)
print("SMALLER -- perception at p3 c0 across every u; group `Ar0e ux p3 c0`")
print("=" * 78)
for u in range(0, 5):
    a = Atom(element="Ar", radical_electrons=u, lone_pairs=3, charge=0)
    print(f"  direct perception of Ar u{u} p3 c0 (no bonds) -> {perceive(a)}")

print()
g = Group().from_adjacency_list("1 *1 Ar0e ux p3 c0")
print("  group `1 *1 Ar0e ux p3 c0` built; atomtypes =",
      [t.label for t in g.atoms[0].atomtype], "radical_electrons =", g.atoms[0].radical_electrons)
for u in range(0, 5):
    m = Molecule(atoms=[Atom(element="Ar", radical_electrons=u, lone_pairs=3, charge=0)])
    m.update_atomtypes(log_species=False, raise_exception=False)
    a = m.atoms[0]
    print(f"    molecule built directly at u{u}: after update -> u{a.radical_electrons} "
          f"p{a.lone_pairs} c{a.charge:+d} type={a.atomtype.label}; "
          f"matches group: {m.is_subgraph_isomorphic(g)}")

print()
print("  and the same group against the other argon molecules:")
for adj, name in (("1 Ar u0 p4 c0", "Ar0"), ("1 Ar u1 p3 c+1", "Ar+"), ("1 Ar u0 p3 c+2", "Ar++")):
    m = Molecule().from_adjacency_list(adj)
    print(f"    {name:4s} {adj!r:18s} type={m.atoms[0].atomtype.label} "
          f"matches: {m.is_subgraph_isomorphic(g)}")
