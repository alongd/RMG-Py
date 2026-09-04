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
"""

import os
import subprocess
import sys

import pytest

from rmgpy.molecule import Molecule
from rmgpy.molecule.group import Group
from rmgpy.molecule.atomtype import ATOMTYPES, get_atomtype

# The three argon types that were kept.
ARGON_KEPT = ["Ar0", "Ar+", "Ar++"]
# The four types dropped as unused / off-theme / isolated.
DROPPED = ["Li-", "Na-", "K-", "N3dc"]

# make_sample_atom expectations for the kept argon types: (charge, lone_pairs).
ARGON_SAMPLE = {"Ar0": (0, 4), "Ar+": (1, 3), "Ar++": (2, 3)}

# The only two argon signatures that occur as concrete species in RMG-database
# (plasma): neutral ground-state argon and the Ar+ radical cation.
ARGON_DB_SIGNATURES = {
    "1 Ar u0 p4 c0": "Ar0",
    "1 Ar u1 p3 c+1": "Ar+",
}

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
    """The argon subgraph (generic ``Ar`` + Ar0/Ar+/Ar++) closes both ways: every
    declared action has its inverse declared back. This is the property the alkali
    anions violated (one-way sink) and the reason they, unlike argon, were not
    kept with only their forward actions."""
    violations = _closure_violations(["Ar", "Ar0", "Ar+", "Ar++"])
    assert violations == set(), (
        "argon action graph is not both-ways closed: {0}".format(sorted(violations))
    )
