#!/usr/bin/env python3

"""
Regression guard for the seven atom types that branch `plasma` was missing
relative to branch `99`: Ar0, Ar+, Ar++, Li-, Na-, K- and N3dc.

An *incompletely* registered atom type does not announce itself as a test
failure. ``GroupAtom.make_sample_atom`` resolves an atom type to an element by
scanning ``allElements`` for a type that either is ``ATOMTYPES[element_label]``
or appears in that element's ``.specific`` list; when neither holds, ``element``
stays ``None`` and the ``mol.Atom(element=None)`` that follows can take the
interpreter down rather than raise. A suite that dies mid-run reads as an
infrastructure problem rather than as a bug in the atom type roster, so each
atom type is exercised in its **own subprocess** here: one crash cannot hide
the other six.
"""

import os
import subprocess
import sys

import pytest

SEVEN = ["Ar0", "Ar+", "Ar++", "Li-", "Na-", "K-", "N3dc"]

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


@pytest.mark.parametrize("label", SEVEN)
def test_make_sample_atom_for_restored_atomtype(label):
    """Each of the seven resolves to a real sample atom, in its own subprocess."""
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
