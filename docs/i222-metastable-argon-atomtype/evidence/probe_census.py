"""I-222: a complete census of every argon atom type -- declared vs measured.

For each argon atom type: its declared feature ranges, its hierarchy links, every concrete atom
that perceives as it, its declared set_actions, and -- for every one of the ten actions applied to
every concrete representative -- what the resulting atom ACTUALLY perceives as. The last column is
the one that matters: it is the only place where "what the table says" and "what the code does" can
be compared.
"""
import itertools

import rmgpy.molecule.atomtype as at
from rmgpy.molecule import Atom, Bond, Molecule
from rmgpy.molecule.atomtype import ATOMTYPES, get_atomtype
from rmgpy.exceptions import AtomTypeError, ActionError
from rmgpy import settings

print("database.directory =", settings["database.directory"])
print("loaded atomtype module __file__ =", at.__file__)
print()

ARGON = ["Ar", "Ar0", "Ar0s", "Ar0e", "Ar+", "Ar++"]
FEATURES = ("single", "all_double", "r_double", "o_double", "s_double",
            "triple", "quadruple", "benzene", "lone_pairs", "charge")
ACTIONS = ("increment_bond", "decrement_bond", "form_bond", "break_bond",
           "increment_radical", "decrement_radical", "increment_lone_pair",
           "decrement_lone_pair", "increment_charge", "decrement_charge")
# action name -> the recipe primitive that drives it
PRIMITIVE = {
    "increment_bond": "CHANGE_BOND +1", "decrement_bond": "CHANGE_BOND -1",
    "form_bond": "FORM_BOND", "break_bond": "BREAK_BOND",
    "increment_radical": "GAIN_RADICAL", "decrement_radical": "LOSE_RADICAL",
    "increment_lone_pair": "GAIN_PAIR", "decrement_lone_pair": "LOSE_PAIR",
    "increment_charge": "GAIN_CHARGE", "decrement_charge": "LOSE_CHARGE",
}

# Concrete representatives: (nickname, radicals, lone_pairs, charge, n_single_bonds)
REPRESENTATIVES = [
    ("Ar0   bare",          0, 4, 0, 0),
    ("Ar0s  bonded u1",     1, 3, 0, 1),
    ("Ar0s  bonded u0",     0, 3, 0, 1),
    ("Ar0e  bare u2",       2, 3, 0, 0),
    ("Ar+   bare u1",       1, 3, 1, 0),
    ("Ar+   bonded u0",     0, 3, 1, 1),
    ("Ar++  bare u0",       0, 3, 2, 0),
    ("Ar++  bonded u1",     1, 3, 2, 1),
]


def build(radicals, pairs, charge, n_bonds):
    atom = Atom(element="Ar", radical_electrons=radicals, lone_pairs=pairs, charge=charge)
    bonds = {}
    for _ in range(n_bonds):
        partner = Atom(element="H", radical_electrons=1, lone_pairs=0, charge=0)
        bonds[partner] = Bond(atom, partner, order=1)
    return atom, bonds


def perceive(atom, bonds):
    try:
        return get_atomtype(atom, bonds).label
    except AtomTypeError:
        return "<none>"


print("=" * 96)
print("1. DECLARATIONS")
print("=" * 96)
hdr = f"{'type':6s} " + " ".join(f"{f[:9]:>9s}" for f in FEATURES)
print(hdr)
print("-" * len(hdr))
for label in ARGON:
    a = ATOMTYPES[label]
    row = f"{label:6s} " + " ".join(f"{str(getattr(a, f)):>9s}" for f in FEATURES)
    print(row)
print()
print("hierarchy links:")
for label in ARGON:
    a = ATOMTYPES[label]
    print(f"  {label:6s} generic={[t.label for t in a.generic]}")
    print(f"  {'':6s} specific={[t.label for t in a.specific]}")
print()

print("=" * 96)
print("2. WHICH CONCRETE ATOMS PERCEIVE AS WHICH TYPE")
print("=" * 96)
print(f"{'atom':34s} {'perceives as':12s}  constructible via adjacency list?")
print("-" * 96)
for nick, u, p, c, nb in REPRESENTATIVES:
    atom, bonds = build(u, p, c, nb)
    got = perceive(atom, bonds)
    if nb == 0:
        adj = f"1 Ar u{u} p{p} c{c:+d}".replace("c+0", "c0")
        try:
            Molecule().from_adjacency_list(adj)
            ok = "yes"
        except Exception as e:
            ok = f"NO ({type(e).__name__})"
    else:
        ok = "n/a (needs a partner)"
    print(f"  {nick:32s} {got:12s}  {ok}")

print()
print("  exhaustive sweep of bond-free argon, u0..u4 x p0..p4 x c-1..+3:")
print(f"  {'u':>2s} {'p':>2s} {'c':>3s}  {'perceives':10s} {'adjlist':10s}")
seen = {}
for u, p, c in itertools.product(range(5), range(5), range(-1, 4)):
    atom, bonds = build(u, p, c, 0)
    got = perceive(atom, bonds)
    if got == "<none>":
        continue
    adj = f"1 Ar u{u} p{p} c{c:+d}".replace("c+0", "c0")
    try:
        Molecule().from_adjacency_list(adj)
        ok = "builds"
    except Exception:
        ok = "refused"
    seen.setdefault(got, []).append((u, p, c, ok))
    print(f"  {u:2d} {p:2d} {c:+3d}  {got:10s} {ok:10s}")
print()
print("  summary -- how many (u,p,c) triples perceive as each type, and how many of those build:")
for label, rows in sorted(seen.items()):
    builds = sum(1 for r in rows if r[3] == "builds")
    print(f"    {label:6s} {len(rows):3d} triples perceive, {builds:3d} of them constructible")

print()
print("=" * 96)
print("3. ACTION GRAPH -- DECLARED vs MEASURED")
print("=" * 96)
for label in ARGON:
    a = ATOMTYPES[label]
    reps = [r for r in REPRESENTATIVES if r[0].split()[0] == label]
    print(f"\n  {label}")
    print(f"    {'action':20s} {'primitive':14s} {'DECLARED':18s} MEASURED (per representative)")
    print("    " + "-" * 90)
    for action in ACTIONS:
        declared = [t.label for t in getattr(a, action)]
        results = []
        for nick, u, p, c, nb in reps:
            atom, bonds = build(u, p, c, nb)
            try:
                if action == "form_bond":
                    partner = Atom(element="H", radical_electrons=1, lone_pairs=0, charge=0)
                    bonds[partner] = Bond(atom, partner, order=1)
                elif action == "break_bond":
                    if not bonds:
                        results.append(f"{nick.split()[1]}:n/a")
                        continue
                    bonds.pop(next(iter(bonds)))
                elif action == "increment_bond":
                    if not bonds:
                        results.append(f"{nick.split()[1]}:n/a")
                        continue
                    next(iter(bonds.values())).order += 1
                elif action == "decrement_bond":
                    if not bonds:
                        results.append(f"{nick.split()[1]}:n/a")
                        continue
                    next(iter(bonds.values())).order -= 1
                else:
                    prim = {"increment_radical": "GAIN_RADICAL", "decrement_radical": "LOSE_RADICAL",
                            "increment_lone_pair": "GAIN_PAIR", "decrement_lone_pair": "LOSE_PAIR",
                            "increment_charge": "GAIN_CHARGE", "decrement_charge": "LOSE_CHARGE"}[action]
                    atom.apply_action([prim, "*1", 1])
                results.append(f"{' '.join(nick.split()[1:]) or 'bare'}:{perceive(atom, bonds)}")
            except (ActionError, Exception) as e:
                results.append(f"{' '.join(nick.split()[1:]) or 'bare'}:{type(e).__name__}")
        flag = ""
        measured = {r.split(":", 1)[1] for r in results if not r.endswith(("n/a", "ActionError"))}
        measured.discard("<none>")
        if set(declared) != measured and (declared or measured):
            flag = "   <<< DISAGREES"
        print(f"    {action:20s} {PRIMITIVE[action]:14s} {str(declared):18s} {', '.join(results)}{flag}")

print()
print("=" * 96)
print("4. THE OPEN QUESTION -- Ar+ LOSE_CHARGE")
print("=" * 96)
print("  declared Ar+.decrement_charge :", [t.label for t in ATOMTYPES["Ar+"].decrement_charge])
for nick, u, p, c, nb in [r for r in REPRESENTATIVES if r[0].startswith("Ar+")]:
    atom, bonds = build(u, p, c, nb)
    atom.apply_action(["LOSE_CHARGE", "*1", 1])
    print(f"  measured {nick:20s} LOSE_CHARGE -> u{atom.radical_electrons} p{atom.lone_pairs} "
          f"c{atom.charge:+d} {len(bonds)} single -> {perceive(atom, bonds)}")
print()
print("  what each candidate inverse would have to declare back (GAIN_CHARGE):")
for nick, u, p, c, nb in REPRESENTATIVES:
    if c != 0:
        continue
    atom, bonds = build(u, p, c, nb)
    before = perceive(atom, bonds)
    atom.apply_action(["GAIN_CHARGE", "*1", 1])
    print(f"    {nick:20s} ({before:6s}) GAIN_CHARGE -> {perceive(atom, bonds)}")
