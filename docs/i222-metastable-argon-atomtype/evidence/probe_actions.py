"""I-222: which single actions, applied to Ar0e, land on a declared argon type?

Also measures which `u` states adjacency-list consistency admits at Ar p3 c0.
"""
import rmgpy.molecule.atomtype as at
from rmgpy.molecule import Atom, Molecule
from rmgpy.molecule.atomtype import get_atomtype
from rmgpy.exceptions import AtomTypeError
from rmgpy import settings

print("database.directory =", settings["database.directory"])
print("loaded atomtype module __file__ =", at.__file__)
print()


def describe(atom, bonds=None):
    bonds = bonds if bonds is not None else {}
    try:
        label = get_atomtype(atom, bonds).label
    except AtomTypeError:
        label = "<no type>"
    return (f"u{atom.radical_electrons} p{atom.lone_pairs} c{atom.charge:+d} "
            f"bonds={sorted(b.order for b in bonds.values())} -> {label}")


print("--- single actions applied to a concrete Ar0e atom (u2 p3 c0, no bonds) ---")
ACTIONS = [
    ("GAIN_PAIR (increment_lone_pair)", ["GAIN_PAIR", "*1", 1]),
    ("LOSE_PAIR (decrement_lone_pair)", ["LOSE_PAIR", "*1", 1]),
    ("GAIN_RADICAL (increment_radical)", ["GAIN_RADICAL", "*1", 1]),
    ("LOSE_RADICAL (decrement_radical)", ["LOSE_RADICAL", "*1", 1]),
    ("GAIN_CHARGE (increment_charge)", ["GAIN_CHARGE", "*1", 1]),
    ("LOSE_CHARGE (decrement_charge)", ["LOSE_CHARGE", "*1", 1]),
]
for name, action in ACTIONS:
    atom = Atom(element="Ar", radical_electrons=2, lone_pairs=3, charge=0)
    try:
        atom.apply_action(action)
        print(f"  {name:36s} {describe(atom)}")
    except Exception as e:  # ActionError etc.
        print(f"  {name:36s} REFUSED: {type(e).__name__}: {e}")

# FORM_BOND is a bond-level action; model it by giving the atom a single bond to an H.
from rmgpy.molecule import Bond
partner = Atom(element="H", radical_electrons=1, lone_pairs=0, charge=0)
atom = Atom(element="Ar", radical_electrons=2, lone_pairs=3, charge=0)
bond = Bond(atom, partner, order=1)
print(f"  {'FORM_BOND (single, to H)':36s} {describe(atom, {partner: bond})}")

print()
print("--- reverse direction: which declared argon types reach Ar0e in one action? ---")
for label in ("Ar0", "Ar0s", "Ar+", "Ar++"):
    a = at.ATOMTYPES[label]
    print(f"  {label:5s} single={a.single} lone_pairs={a.lone_pairs} charge={a.charge}")
    starts = {
        "Ar0": dict(radical_electrons=0, lone_pairs=4, charge=0),
        "Ar0s": dict(radical_electrons=1, lone_pairs=3, charge=0),
        "Ar+": dict(radical_electrons=1, lone_pairs=3, charge=1),
        "Ar++": dict(radical_electrons=0, lone_pairs=3, charge=2),
    }[label]
    for name, action in ACTIONS:
        atom = Atom(element="Ar", **starts)
        bonds = {}
        if label == "Ar0s":
            p = Atom(element="H", radical_electrons=1, lone_pairs=0, charge=0)
            bonds = {p: Bond(atom, p, order=1)}
        try:
            atom.apply_action(action)
            res = describe(atom, bonds)
        except Exception as e:
            res = f"REFUSED: {type(e).__name__}"
        flag = "  <== Ar0e" if res.endswith("-> Ar0e") else ""
        print(f"      {name:36s} {res}{flag}")

print()
print("--- which u states does the adjacency-list path admit at Ar p3 c0 (no bonds)? ---")
for u in range(0, 5):
    adj = f"1 Ar u{u} p3 c0"
    try:
        m = Molecule().from_adjacency_list(adj)
        print(f"  {adj!r:22s} -> OK, atomtype = {m.atoms[0].atomtype.label}")
    except Exception as e:
        print(f"  {adj!r:22s} -> {type(e).__name__}: {' '.join(str(e).split())[:150]}")
