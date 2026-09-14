#!/usr/bin/env python3

"""
Guard for the disposition of the seven atom types that commit dbd131221 restored
to ``rmgpy/molecule/atomtype.py`` -- Ar0, Ar+, Ar++, Li-, Na-, K- and N3dc --
after i159 established, by measurement, which of them are correct and needed.

Verdict encoded here (see docs/i159-atomtypes/report.md for the evidence):

* KEPT  -- ``Ar0`` / ``Ar+`` / ``Ar++`` (a faithful port of branch 99). They
  give the plasma campaign the neutral/singly/doubly-ionised argon distinction an
  argon discharge is built on, they type the only two argon signatures that occur
  in the RMG-database ``plasma`` branch, and their action graph closes both ways.

* DROPPED -- ``Li-`` / ``Na-`` / ``K-`` (removal commit b52045138's central
  reason -- no carried reaction, group or dictionary uses an alkali anion -- is
  still true: an exhaustive search of the database found zero uses) and ``N3dc``
  (zero database use, absent from mainline, unwired even in branch 99, and an
  isolated node that raises ActionError the moment a recipe touches it).

The *old* membership + no-crash version of this file PASSED against the broken tip
(dbd131221): there the four dropped types are still registered and the alkali-anion
action graph is a one-way sink, neither of which the old assertions could see. This
reworked file does NOT pass against that tip: its ``test_dropped_types_are_absent``
and closure assertions FAIL there (measured: 4 failed, 7 passed) and pass only once
the change set is reduced to the argon three.

An *incompletely* registered atom type does not announce itself as a test
failure: ``GroupAtom.make_sample_atom`` resolves an atom type to an element by
scanning ``allElements``; when no element matches, ``element`` stays ``None`` and
the ``mol.Atom(element=None)`` that follows can take the interpreter down rather
than raise. So the kept types are each exercised in their **own subprocess**: one
crash cannot hide the others.

The argon leaf set has since GROWN, and this file guards the new shape as well as
the i159 verdict:

* i218 added ``Ar0s``, the singly-bonded neutral -- the neutral half of Ar2+.
* i222 added ``Ar0e``, the bond-free neutral at three lone pairs -- metastable Ar*.

So "the argon three" is now the argon **five**, and three statements this file used
to make are no longer the whole truth. Each has been replaced by a test rather than
by an edited comment:

* the registry census (``test_argon_leaf_set_is_exactly_these_five``) fails if a
  sixth leaf appears or one of the five disappears, which the old membership
  assertions could not see;
* ``ARGON_DB_SIGNATURES`` still holds the only two argon species that occur in
  RMG-database (plasma) -- verified by grep at the time of writing, ``Ar u0 p4 c0``
  and ``Ar u1 p3 c+1``. Three of the five leaves therefore have NO concrete
  database species, which is now stated outright rather than left implied;
* the action-closure test covered only ``Ar``/``Ar0``/``Ar+``/``Ar++``, so the two
  newer leaves were outside it exactly while i222 was rewiring the charge edges
  through them. It now covers all five.
"""

import os
import subprocess
import sys

import pytest

from rmgpy.molecule import Molecule
from rmgpy.molecule.group import Group
from rmgpy.molecule.atomtype import ATOMTYPES, get_atomtype

# The three argon types that were kept by i159.
ARGON_KEPT = ["Ar0", "Ar+", "Ar++"]
# Added after i159: the singly-bonded neutral (i218) and the metastable (i222).
ARGON_ADDED = ["Ar0s", "Ar0e"]
# Every argon leaf that must be registered, in the order ATOMTYPES['Ar'].specific holds them.
ARGON_ALL = ["Ar0", "Ar0s", "Ar0e", "Ar+", "Ar++"]
# The four types dropped as unused / off-theme / isolated.
DROPPED = ["Li-", "Na-", "K-", "N3dc"]

# make_sample_atom expectations, per argon leaf: (charge, lone_pairs).
# Ar0s and Ar0e share a sample atom -- see test_added_argon_types_share_one_sample_atom.
ARGON_SAMPLE = {"Ar0": (0, 4), "Ar+": (1, 3), "Ar++": (2, 3), "Ar0s": (0, 3), "Ar0e": (0, 3)}

# The only two argon signatures that occur as concrete species in RMG-database
# (plasma): neutral ground-state argon and the Ar+ radical cation.
ARGON_DB_SIGNATURES = {
    "1 Ar u0 p4 c0": "Ar0",
    "1 Ar u1 p3 c+1": "Ar+",
}
# The leaves no database species reaches. Ar0s and Ar0e need a bonded Ar2+ and a
# metastable species respectively, both owned by other tickets; Ar++ has never had one.
ARGON_WITHOUT_DB_SPECIES = ["Ar0s", "Ar0e", "Ar++"]

# The inverse-action pairs that define both-ways action-graph closure.
_INVERSE_ACTIONS = [
    ("increment_bond", "decrement_bond"),
    ("form_bond", "break_bond"),
    ("increment_radical", "decrement_radical"),
    ("increment_lone_pair", "decrement_lone_pair"),
    ("increment_charge", "decrement_charge"),
]

# Run in a child interpreter: an incomplete registration may segfault rather than raise.
CHILD = """
import sys
from rmgpy.molecule.group import GroupAtom
from rmgpy.molecule.atomtype import ATOMTYPES

label = sys.argv[1]
if label not in ATOMTYPES:
    sys.exit("UNREGISTERED: {0} is not in ATOMTYPES".format(label))
atom = GroupAtom(atomtype=[ATOMTYPES[label]]).make_sample_atom()
assert atom is not None, "make_sample_atom returned None for {0}".format(label)
print("OK {0} -> element={1} charge={2} lone_pairs={3}".format(
    label, atom.symbol, atom.charge, atom.lone_pairs))
"""


def _run_child(label):
    """Exercise one atom type in a fresh interpreter loading this same source tree."""
    tree = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))))
    env = dict(os.environ)
    env["PYTHONPATH"] = tree + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        [sys.executable, "-c", CHILD, label],
        capture_output=True, text=True, timeout=300, cwd=tree, env=env,
    )


def _closure_violations(labels):
    """Return the both-ways closure violations restricted to ``labels``.

    A violation ``(atomtype, action, neighbour)`` means ``atomtype`` names
    ``neighbour`` under ``action`` but ``neighbour`` does not name ``atomtype``
    under the inverse action -- an edge in the action graph with no return edge.
    Restricting to ``labels`` keeps this test focused on the argon subgraph; the
    whole-table invariant lives in atomtypeTest.py.
    """
    def names(action_list):
        return [t.label for t in action_list]

    violations = set()
    for label in labels:
        atomtype = ATOMTYPES[label]
        for fwd, rev in _INVERSE_ACTIONS:
            for action, inverse in ((fwd, rev), (rev, fwd)):
                for neighbour in names(getattr(atomtype, action)):
                    if label not in names(getattr(ATOMTYPES[neighbour], inverse)):
                        violations.add((label, action, neighbour))
    return violations


@pytest.mark.parametrize("label", ARGON_KEPT)
def test_kept_argon_makes_correct_sample_atom(label):
    """Each kept argon type resolves to a real sample atom with the right charge
    and lone-pair count, in its own subprocess (an incomplete registration would
    segfault, not raise)."""
    proc = _run_child(label)
    assert proc.returncode == 0, (
        "atom type {0!r} failed: returncode={1}{2}\nstdout: {3}\nstderr: {4}".format(
            label,
            proc.returncode,
            " (killed by signal {0} -- a crash, i.e. an incomplete "
            "registration)".format(-proc.returncode) if proc.returncode < 0 else "",
            proc.stdout,
            proc.stderr,
        )
    )
    assert proc.stdout.startswith("OK "), proc.stdout
    charge, lone_pairs = ARGON_SAMPLE[label]
    assert "charge={0}".format(charge) in proc.stdout, proc.stdout
    assert "lone_pairs={0}".format(lone_pairs) in proc.stdout, proc.stdout


@pytest.mark.parametrize("label", DROPPED)
def test_dropped_types_are_absent(label):
    """The four unused / off-theme / isolated types must not be registered.

    Fails against tip dbd131221, where all four are present."""
    assert label not in ATOMTYPES, (
        "{0!r} is registered but has no consumer in the carried chemistry; "
        "i159 dropped it. See docs/i159-atomtypes/report.md.".format(label)
    )


@pytest.mark.parametrize("adjlist,expected", sorted(ARGON_DB_SIGNATURES.items()))
def test_concrete_argon_species_type_specifically(adjlist, expected):
    """The two argon signatures that occur in the database type to their specific
    argon atom type -- not to the lumped generic ``Ar`` -- now that ``Ar`` is out
    of ``nonSpecifics``."""
    mol = Molecule().from_adjacency_list(adjlist)
    atom = mol.atoms[0]
    atomtype = get_atomtype(atom, {b: bd for b, bd in atom.bonds.items()})
    assert atomtype.label == expected, (
        "{0!r} typed as {1}, expected {2}".format(adjlist, atomtype.label, expected)
    )


def test_generic_argon_group_matches_specific_argon_atom():
    """A group written with the generic ``Ar`` atom type still matches a concrete
    neutral argon molecule whose atom types as ``Ar0`` -- i.e. the generic/specific
    hierarchy carries the group match, so removing ``Ar`` from ``nonSpecifics``
    does not orphan any (hypothetical) generic-``Ar`` group."""
    group = Group().from_adjacency_list("1 Ar u0")
    neutral_argon = Molecule().from_adjacency_list("1 Ar u0 p4 c0")
    assert neutral_argon.is_subgraph_isomorphic(group)


def test_argon_action_graph_closes_both_ways():
    """The argon subgraph closes both ways: every declared action has its inverse
    declared back. This is the property the alkali anions violated (one-way sink)
    and the reason they, unlike argon, were not kept with only their forward actions.

    Covers all five leaves. It used to cover only Ar/Ar0/Ar+/Ar++, which left Ar0s
    and Ar0e outside the check exactly while i222 rewired the charge edges through
    them -- ``Ar+.decrement_charge`` now names both."""
    violations = _closure_violations(["Ar"] + ARGON_ALL)
    assert violations == set(), (
        "argon action graph is not both-ways closed: {0}".format(sorted(violations))
    )


def test_argon_leaf_set_is_exactly_these_five():
    """The registry census: exactly five argon leaves, each linked to generic ``Ar``
    in both directions, and generic ``Ar`` naming them in the order it resolves them.

    ``get_atomtype`` only ever returns entries of ``ATOMTYPES['Ar'].specific``, so a
    leaf missing from that list is unreachable however well it is declared, and a
    leaf present but unregistered is an AttributeError waiting to happen. A plain
    membership assertion per type cannot see a SIXTH leaf arriving; this can."""
    assert [t.label for t in ATOMTYPES["Ar"].specific] == ARGON_ALL
    for label in ARGON_ALL:
        assert label in ATOMTYPES, "{0} is not registered".format(label)
        assert ATOMTYPES["Ar"] in ATOMTYPES[label].generic, (
            "{0} does not name generic Ar".format(label)
        )
        assert ATOMTYPES[label].specific == [], (
            "{0} is a leaf and must have no specifics".format(label)
        )


@pytest.mark.parametrize("label", ARGON_ADDED)
def test_added_argon_types_make_correct_sample_atom(label):
    """The two leaves added after i159 resolve to a real sample atom too, each in its
    own subprocess, for the reason in the module docstring."""
    proc = _run_child(label)
    assert proc.returncode == 0, (
        "atom type {0!r} failed: returncode={1}\nstdout: {2}\nstderr: {3}".format(
            label, proc.returncode, proc.stdout, proc.stderr
        )
    )
    charge, lone_pairs = ARGON_SAMPLE[label]
    assert "charge={0}".format(charge) in proc.stdout, proc.stdout
    assert "lone_pairs={0}".format(lone_pairs) in proc.stdout, proc.stdout


def test_added_argon_types_share_one_sample_atom():
    """Ar0s and Ar0e are indistinguishable to the sample builder, and that is why
    both sit in ``EXPECTED_FAILING_ATOMTYPES`` in atomtypeTest.py.

    ``make_sample_atom`` takes the first entry of each feature list and has no rule
    for choosing ``u`` or for adding a single-bonded partner, so it builds the same
    ``Ar u0 p3 c0`` for both -- charge-inconsistent for either, which is what makes
    ``make_sample_molecule`` raise. The two types differ only in ``single``, the one
    feature the builder does not act on."""
    from rmgpy.molecule.group import GroupAtom

    samples = {}
    for label in ARGON_ADDED:
        atom = GroupAtom(atomtype=[ATOMTYPES[label]]).make_sample_atom()
        samples[label] = (atom.symbol, atom.radical_electrons, atom.lone_pairs, atom.charge)
    assert samples["Ar0s"] == samples["Ar0e"] == ("Ar", 0, 3, 0)
    assert ATOMTYPES["Ar0s"].single == [1] and ATOMTYPES["Ar0e"].single == [0]


@pytest.mark.parametrize("label", ARGON_ADDED)
def test_added_argon_types_cannot_make_a_sample_molecule(label):
    """The consequence of the above, pinned as behaviour: the sample atom is
    charge-inconsistent, so building a molecule from it raises. If this ever starts
    passing, ``make_sample_atom`` has learned to pick ``u`` and both labels should
    come out of ``EXPECTED_FAILING_ATOMTYPES``."""
    from rmgpy.exceptions import UnexpectedChargeError

    group = Group().from_adjacency_list("1 {0} ux".format(label))
    with pytest.raises(UnexpectedChargeError):
        group.make_sample_molecule()


@pytest.mark.parametrize("label", ARGON_WITHOUT_DB_SPECIES)
def test_argon_leaves_without_a_database_species(label):
    """Three of the five leaves are declared but unexercised by any database species.

    Stated outright so the gap is visible: ``ARGON_DB_SIGNATURES`` is the whole set
    of argon species in RMG-database (plasma), and it produces only Ar0 and Ar+.
    Whatever these three do is therefore pinned by unit tests alone -- the point the
    i222 report makes under "what this could not reach". When a metastable argon
    species lands, move ``Ar0e`` into ``ARGON_DB_SIGNATURES`` and delete it here."""
    typed = set()
    for adjlist in ARGON_DB_SIGNATURES:
        mol = Molecule().from_adjacency_list(adjlist)
        atom = mol.atoms[0]
        typed.add(get_atomtype(atom, {b: bd for b, bd in atom.bonds.items()}).label)
    assert label in ATOMTYPES, "{0} is not registered".format(label)
    assert label not in typed, (
        "{0} now has a database species -- move it into ARGON_DB_SIGNATURES".format(label)
    )
