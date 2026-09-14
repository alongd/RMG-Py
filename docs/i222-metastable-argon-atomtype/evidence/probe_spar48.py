"""Measure the claims spar round 48 asks I-222 to WRITE DOWN (items 1, 2 and 3).

Nothing here changes code. It exists so that the sentences added to the declaration
comments and to the report are measurements rather than readings of the source.

Item 1 -- the four OLD argon edges commit 62f61558f replaced were each false. Re-measured
         here by applying the primitive to a concrete atom and re-perceiving it.
Item 2 -- Ar+.decrement_charge = ['Ar0e', 'Ar0s'] states a bond-conditioned result as an
         unconditional union. Measured: the two concrete answers; that
         GroupAtom._lose_charge keeps both and never looks at bond count; that the order
         it keeps them in is not the declared order, because it goes through set(); and
         that make_sample_atom and pick_wildcards read element [0] of that list.
Item 3 -- a group spelled `1 Ar0e ux p3 c0` matches non-metastable argon. (Already pinned
         by test_group_spelled_ar0e_matches_more_than_the_metastable_triplet; re-measured
         here so the declaration comment quotes a number it owns.)
"""

import os
import subprocess
import sys

import rmgpy.molecule.atomtype as at
from rmgpy.molecule import Atom, Molecule
from rmgpy.molecule.atomtype import ATOMTYPES, get_atomtype
from rmgpy.molecule.group import Group, GroupAtom


def rule(title):
    print("\n" + title)
    print("-" * len(title))


rule("0. what is loaded")
print("atomtype module:", at.__file__)
print("Ar+.decrement_charge as declared:", [t.label for t in ATOMTYPES["Ar+"].decrement_charge])

rule("1. the four edges 62f61558f replaced, re-measured on concrete atoms")


def apply_primitive(adjlist, primitive):
    """Apply the MOLECULE-path primitive to atom 0 of ``adjlist`` and re-perceive it.

    This is `Atom.apply_action` itself -- the same call the recipe makes -- not a hand
    mutation of the fields, so the answer is the one a running reaction would get.
    """
    from rmgpy.exceptions import ActionError, AtomTypeError

    molecule = Molecule().from_adjacency_list(adjlist)
    atom = molecule.atoms[0]
    try:
        atom.apply_action([primitive, "*1", 1])
    except ActionError as exc:
        return "REFUSED ({0})".format(exc)
    try:
        return get_atomtype(atom, atom.bonds).label
    except AtomTypeError:
        return "NO TYPE"


print("Ar0 LOSE_PAIR    -> {0:9s} (old edge said Ar+)".format(
    apply_primitive("1 Ar u0 p4 c0", "LOSE_PAIR")))
print("Ar0 GAIN_CHARGE  -> {0:9s} (old edge said Ar+)".format(
    apply_primitive("1 Ar u0 p4 c0", "GAIN_CHARGE")))
print("Ar+ GAIN_PAIR    -> {0:9s} (old edge said Ar0)".format(
    apply_primitive("1 Ar u1 p3 c+1", "GAIN_PAIR")))
print("Ar+ LOSE_CHARGE  -> {0:9s} (old edge said Ar0)".format(
    apply_primitive("1 Ar u1 p3 c+1", "LOSE_CHARGE")))

rule("2a. LOSE_CHARGE on a concrete Ar+ has TWO answers, chosen by bond count")
print("bare Ar+          -> ", apply_primitive("1 Ar u1 p3 c+1", "LOSE_CHARGE"))
dimer = Molecule().from_adjacency_list("1 Ar u0 p3 c+1 {2,S}\n2 Ar u1 p3 c0 {1,S}")
bonded = dimer.atoms[0]
bonded.apply_action(["LOSE_CHARGE", "*1", 1])
print("singly-bonded Ar+ -> ", get_atomtype(bonded, bonded.bonds).label)

rule("2b/2c. LOSE_CHARGE on a GroupAtom spelled Ar+, via the GROUP path")
ga = GroupAtom(atomtype=[ATOMTYPES["Ar+"]], charge=[1])
ga.apply_action(["LOSE_CHARGE", "*1", 1])
print("resulting atomtype list:", [t.label for t in ga.atomtype])
print("bond count was never consulted -- _lose_charge reads decrement_charge and nothing else")

# The order is whatever list(set(...)) yields. Re-measure it in fresh interpreters, where
# object identities -- and therefore the default hashes -- differ from run to run.
CHILD = (
    "from rmgpy.molecule.group import GroupAtom;"
    "from rmgpy.molecule.atomtype import ATOMTYPES;"
    "a=GroupAtom(atomtype=[ATOMTYPES['Ar+']], charge=[1]);"
    "a.apply_action(['LOSE_CHARGE','*1',1]);"
    "print([t.label for t in a.atomtype])"
)
for i in range(6):
    out = subprocess.run([sys.executable, "-c", CHILD], capture_output=True, text=True)
    print("fresh interpreter {0}: {1}".format(i, out.stdout.strip() or out.stderr.strip()[-200:]))

rule("2d. the consumers that read element [0] of that list")
# group.py is cythonized, so inspect.getsource cannot reach these; read the source file.
source = os.path.join(os.path.dirname(at.__file__), "group.py")
with open(source) as handle:
    lines = handle.read().splitlines()
for name in ("def make_sample_atom", "def pick_wildcards", "def _lose_charge"):
    start = next(n for n, line in enumerate(lines) if line.strip().startswith(name))
    body = lines[start:start + 60]
    end = next((n for n, line in enumerate(body[1:], 1) if line.startswith("    def ")), len(body))
    hits = [(start + n + 1, line.strip()) for n, line in enumerate(body[:end])
            if "atomtype[0]" in line or "list(set(" in line]
    print("{0:22s} line {1}: {2}".format(name[4:], start + 1, hits))

rule("3. a group spelled `Ar0e ux p3 c0`, against every u at p3 c0")
group = Group().from_adjacency_list("1 *1 Ar0e ux p3 c0")
for u in range(5):
    molecule = Molecule(atoms=[Atom(element="Ar", radical_electrons=u, lone_pairs=3, charge=0)])
    molecule.update_atomtypes(log_species=False, raise_exception=False)
    print("u{0} p3 c0 -> perceived {1:5s} matched by the group: {2}".format(
        u, molecule.atoms[0].atomtype.label, molecule.is_subgraph_isomorphic(group)))
