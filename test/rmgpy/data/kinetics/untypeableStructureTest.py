#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@northeastern.edu) and the RMG Team            #
# (rmg_dev@mit.edu)                                                           #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL     #
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER  #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
Regression tests for I-262: when a reaction family recipe builds a product molecule that no atom
type can represent, generation must raise a *named* refusal (``UntypeableStructureError``) that
names the family, the matched top-level template, and the reacting species -- not the anonymous
``AtomTypeError`` that names only the impossible atom's shape.

The impossible-argon cells these tests drive live only in the plasma RMG-database, so the tests are
marked ``database``. They assert on the *mechanism* (class, hierarchy, named facts, round-trip),
not on exact database strings, so they survive an owner editing the templates.
"""

import os

import pytest

from rmgpy import settings
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import UntypeableStructureError
from rmgpy.exceptions import AtomTypeError
from rmgpy.molecule import Molecule


def _mol(adj):
    return Molecule().from_adjacency_list(adj, saturate_h=False)


# (family, [reactant adjacency lists]) -- two distinct families building distinct impossible argon.
CELLS = {
    "R_Recombination": [
        "1 Ar u0 p3 c+1 {2,S}\n2 Ar u1 p3 c0 {1,S}",  # Ar2+
        "1 H u1 p0 c0",                                # H
    ],
    "Birad_R_Recombination": [
        "1 Ar u0 p3 c+1 {2,S}\n2 Ar u1 p3 c0 {1,S}",  # Ar2+
        "1 Ar u2 p3 c0",                              # Ar(3P2) metastable
    ],
}


@pytest.mark.database
class TestUntypeableStructureRefusal:
    @classmethod
    def setup_class(cls):
        cls.database = KineticsDatabase()
        cls.database.load_families(
            os.path.join(settings["database.directory"], "kinetics", "families"),
            families=sorted(CELLS.keys()),
        )

    def _run(self, family):
        reactants = [_mol(a) for a in CELLS[family]]
        with pytest.raises(UntypeableStructureError) as exc_info:
            self.database.generate_reactions_from_families(reactants, only_families=[family])
        return exc_info.value

    def test_class_is_not_an_atomtypeerror_subclass(self):
        # A direct Exception subclass so no `except AtomTypeError` can silently swallow the refusal.
        assert issubclass(UntypeableStructureError, Exception)
        assert not issubclass(UntypeableStructureError, AtomTypeError)

    def test_refusal_names_family_and_chains_cause(self):
        exc = self._run("R_Recombination")
        assert exc.family == "R_Recombination"
        assert "R_Recombination" in str(exc)
        # The original anonymous AtomTypeError is chained, nothing lost.
        assert isinstance(exc.__cause__, AtomTypeError)
        assert "Unable to determine atom type" in str(exc)

    def test_two_cells_discriminate(self):
        # The whole point: two cells that gave the IDENTICAL anonymous AtomTypeError message now
        # name different families and different species.
        a = self._run("R_Recombination")
        b = self._run("Birad_R_Recombination")
        assert a.family != b.family
        assert str(a) != str(b)
        assert a.reactants != b.reactants

    def test_printed_reactant_adjlists_round_trip(self):
        exc = self._run("R_Recombination")
        assert exc.reactants  # the message carries the reacting species as adjacency lists
        for adj in exc.reactants:
            # Must re-parse without error -- a printed adj list that cannot be re-parsed is the
            # sibling-ticket defect this guards against.
            Molecule().from_adjacency_list(adj, saturate_h=False)
