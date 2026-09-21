#!/usr/bin/env python3
"""
Reproduce both halves of the group/pairwise mismatch in ``mark_duplicate_reaction``.

The engine under test is chosen by ``sys.argv[1]``: the script chdirs into that
worktree *before* importing rmgpy, because this checkout's editable install
resolves ``import rmgpy`` from the current working directory. The loaded module's
path and a value-level fingerprint are printed first, so the arm can never be
confused with the other one.

Half A -- the pre-marked bypass. A pair whose electron placements differ but which
arrives with ``duplicate = True`` on both sides never reaches the refined branch at
all: the refinement lives in the ``else`` of ``if reaction1.duplicate and
reaction2.duplicate``.

Half B -- the cross-group clearing. Four legitimately pre-marked reactions, two
groups of two, where a cross-group comparison hits the mixed-pressure-dependence
un-marking branch and clears flags that belong to a *different* group. The deck
then carries two identical equations with no ``DUPLICATE`` line, which Chemkin
rejects.
"""

import os
import sys

WORKTREE = sys.argv[1]
os.chdir(WORKTREE)
sys.path.insert(0, WORKTREE)

import rmgpy.chemkin as chemkin_module  # noqa: E402
from rmgpy.chemkin import mark_duplicate_reaction, mark_duplicate_reactions, save_chemkin_file  # noqa: E402
from rmgpy.data.kinetics.library import LibraryReaction  # noqa: E402
from rmgpy.electron_balance import get_electron_placement_counts  # noqa: E402
from rmgpy.kinetics import Arrhenius, Chebyshev  # noqa: E402
from rmgpy.molecule import Molecule  # noqa: E402
from rmgpy.species import Species  # noqa: E402
from rmgpy.thermo import NASA, NASAPolynomial  # noqa: E402

DECLARED_TWO_SIDED_OWNER = "PlasmaElectronImpactIonization"
UNDECLARED_OWNER = "AnUndeclaredKineticsLibrary"


def make_species(label, index, molecule):
    spc = Species(label=label, molecule=[molecule])
    spc.index = index
    coeffs = [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, -11.7246]
    spc.thermo = NASA(
        polynomials=[
            NASAPolynomial(coeffs=coeffs, Tmin=(200, "K"), Tmax=(1000, "K")),
            NASAPolynomial(coeffs=coeffs, Tmin=(1000, "K"), Tmax=(6000, "K")),
        ],
        Tmin=(200, "K"), Tmax=(6000, "K"),
    )
    return spc


def species_trio():
    electron = make_species("e-", 1, Molecule().from_adjacency_list("1 e u1 p0 c-1"))
    li = make_species("Li", 2, Molecule().from_adjacency_list("multiplicity 2\n1 Li u1 p0 c0\n"))
    lip = make_species("Lip", 3, Molecule().from_adjacency_list("1 Li u0 p0 c+1\n"))
    return electron, li, lip


def library_reaction(reactants, products, owner, electrons, kinetics, duplicate=False):
    rxn = LibraryReaction(reactants=list(reactants), products=list(products),
                          library=owner, kinetics=kinetics, reversible=False,
                          duplicate=duplicate)
    rxn.electrons = electrons
    return rxn


def arrhenius(units, A=1.0e12):
    return Arrhenius(A=(A, units), n=0.0, Ea=(0.0, "kcal/mol"))


def chebyshev():
    return Chebyshev(coeffs=[[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                     kunits="s^-1", Tmin=(300, "K"), Tmax=(2000, "K"),
                     Pmin=(0.01, "bar"), Pmax=(100, "bar"))


def deck_entries(text):
    entries, in_reactions = [], False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("REACTIONS"):
            in_reactions = True
            continue
        if not in_reactions:
            continue
        if stripped == "END":
            break
        if stripped == "DUPLICATE":
            equation, _ = entries[-1]
            entries[-1] = (equation, True)
            continue
        if line[:1].strip() and ("=>" in line or "<=>" in line or "=" in line):
            entries.append((line.split()[0], False))
    return entries


def half_a():
    print("=" * 78)
    print("HALF A -- pre-marked bypass: the refinement sits in the else of the flag test")
    print("=" * 78)
    for arriving in [(True, True), (True, False), (False, True), (False, False)]:
        _electron, li, lip = species_trio()
        declared = library_reaction([li], [lip], DECLARED_TWO_SIDED_OWNER, 1,
                                    arrhenius("cm^3/(mol*s)"), duplicate=arriving[0])
        undeclared = library_reaction([li], [lip], UNDECLARED_OWNER, 1,
                                      arrhenius("s^-1", 2.0e12), duplicate=arriving[1])
        counts = (get_electron_placement_counts(declared), get_electron_placement_counts(undeclared))
        mark_duplicate_reaction(declared, [undeclared])
        flags = (declared.duplicate, undeclared.duplicate)
        verdict = "DUPLICATE pair emitted" if all(flags) else ("clean" if not any(flags) else "one-sided")
        print("  arrive {0!s:<14} counts {1} vs {2} -> flags {3!s:<14}  {4}".format(
            arriving, counts[0], counts[1], flags, verdict))
    print()


def half_a_deck(tmpdir):
    """
    The same pair as half A, taken all the way to a deck. This is the arm that
    matters: the pairwise function is an incremental hint and is left as it is,
    while the deck is what Chemkin reads.
    """
    print("=" * 78)
    print("HALF A (deck) -- the same pre-marked pair written through the real writer")
    print("=" * 78)
    for arriving in [(True, True), (False, False)]:
        electron, li, lip = species_trio()
        declared = library_reaction([li], [lip], DECLARED_TWO_SIDED_OWNER, 1,
                                    arrhenius("cm^3/(mol*s)"), duplicate=arriving[0])
        undeclared = library_reaction([li], [lip], UNDECLARED_OWNER, 1,
                                      arrhenius("s^-1", 2.0e12), duplicate=arriving[1])
        path = os.path.join(tmpdir, "half_a_{0}.inp".format(arriving[0]))
        save_chemkin_file(path, [electron, li, lip], [declared, undeclared], verbose=False)
        with open(path) as f:
            entries = deck_entries(f.read())
        print("  arrive {0!s:<14} -> deck {1}".format(
            arriving, ", ".join("{0} DUPLICATE={1}".format(e, m) for e, m in entries)))
        if all(m for _, m in entries):
            print("      >>> two different stoichiometries both stamped DUPLICATE"
                  "   <== CHEMKIN REJECTS THIS DECK")
    print()


def half_b(tmpdir):
    print("=" * 78)
    print("HALF B -- cross-group clearing: an un-marking branch clears another group's flags")
    print("=" * 78)
    electron, li, lip = species_trio()
    # Group A: two genuine duplicates of the declared (1, 2) channel, not pressure dependent.
    a1 = library_reaction([li], [lip], DECLARED_TWO_SIDED_OWNER, 1, arrhenius("cm^3/(mol*s)", 1.0e12), duplicate=True)
    a2 = library_reaction([li], [lip], DECLARED_TWO_SIDED_OWNER, 1, arrhenius("cm^3/(mol*s)", 2.0e12), duplicate=True)
    # Group B: two genuine duplicates of the net-derived (0, 1) channel, pressure dependent.
    b1 = library_reaction([li], [lip], UNDECLARED_OWNER, 1, chebyshev(), duplicate=True)
    b2 = library_reaction([li], [lip], "AnotherUndeclaredKineticsLibrary", 1, chebyshev(), duplicate=True)
    reactions = [a1, a2, b1, b2]
    names = ["A1 (1,2) non-pdep", "A2 (1,2) non-pdep", "B1 (0,1) pdep", "B2 (0,1) pdep"]

    print("  arriving flags: " + ", ".join(
        "{0}={1}".format(n, r.duplicate) for n, r in zip(names, reactions)))
    print("  placements:     " + ", ".join(
        "{0}={1}".format(n.split()[0], get_electron_placement_counts(r)) for n, r in zip(names, reactions)))

    mark_duplicate_reactions(reactions)
    print("  resulting flags: " + ", ".join(
        "{0}={1}".format(n, r.duplicate) for n, r in zip(names, reactions)))

    path = os.path.join(tmpdir, "half_b.inp")
    save_chemkin_file(path, [electron, li, lip], reactions, verbose=False)
    with open(path) as f:
        text = f.read()
    entries = deck_entries(text)
    print("  deck entries:")
    for equation, marked in entries:
        print("    {0:<40} DUPLICATE={1}".format(equation, marked))
    equations = [e for e, _ in entries]
    for equation in sorted(set(equations)):
        n = equations.count(equation)
        if n > 1:
            marks = [m for e, m in entries if e == equation]
            print("    >>> {0} identical entries of {1!r}, DUPLICATE marks {2}{3}".format(
                n, equation, marks,
                "   <== CHEMKIN REJECTS THIS DECK" if not all(marks) else ""))
        else:
            marks = [m for e, m in entries if e == equation]
            if any(marks):
                print("    >>> lone entry {0!r} carries DUPLICATE with no mate"
                      "   <== CHEMKIN REJECTS THIS DECK".format(equation))
    print()
    print("  REACTIONS block:")
    started = False
    for line in text.splitlines():
        if line.strip().startswith("REACTIONS"):
            started = True
        if started:
            print("    | " + line)
            if line.strip() == "END":
                break
    print()


def main():
    print("engine worktree : {0}".format(WORKTREE))
    print("chemkin module  : {0}".format(chemkin_module.__file__))
    doc = mark_duplicate_reaction.__doc__ or ""
    print("fingerprint     : 'Electrons are participants' in docstring = {0}".format(
        "Electrons are participants" in doc))
    print()
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        half_a()
        half_a_deck(tmpdir)
        half_b(tmpdir)


if __name__ == "__main__":
    main()
