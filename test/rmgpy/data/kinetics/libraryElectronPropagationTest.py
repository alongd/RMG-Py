#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
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
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
A charge-transfer reaction carries its electron count on its rate law
(``SurfaceChargeTransfer``, ``ArrheniusChargeTransfer``, ... all take ``electrons=``), not among
its reactants and products. ``Reaction.is_balanced`` folds that count into the net charges before
comparing them, so a charged reaction only balances if the count has been copied off the kinetics
onto the ``Reaction`` first.

``KineticsDepository.load`` does that copy. ``KineticsLibrary.load`` did not: it read only
``entry.item.electrons``, which for a library entry is zero, because ``KineticsLibrary.load_entry``
takes no electron argument. So a perfectly balanced charged library entry reached the balance check
declaring zero electrons and was rejected.

This went unseen because the two halves concealed each other. ``is_balanced`` computed the charge
balance but returned before comparing it, so nothing ever asked whether the count had been
propagated, and the loader was never required to feed the guard correctly. The comparison is now
live; this module locks the measurement that feeds it.

Two fixture libraries, same reaction, differing only in the declared count:

* ``charged-balanced`` -- ``electrons = -1``, which is exactly what closes the +1/0 charge gap.
  It must load, and its reaction must carry ``electrons == -1``.
* ``charged-unbalanced`` -- ``electrons = -2``, which does not close it. It must keep being
  rejected.

The second is the tripwire. Propagating the count and *weakening the balance check* both make the
first library load; only the propagation leaves the second rejected. Without it, "the charged entry
loads" would not distinguish "it is balanced" from "it was excused".
"""

import os.path

import pytest

from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.library import KineticsLibrary
from rmgpy.exceptions import DatabaseError

FIXTURE_ROOT = os.path.join(os.path.dirname(__file__), "library_electron_propagation_data")


class TestLibraryElectronPropagation:
    """``KineticsLibrary.load`` must copy an entry's declared electron count onto its reaction."""

    @staticmethod
    def _load_library(label):
        """Load one fixture library, returning it and whatever ``load`` raised (``None`` if it
        did not raise). The library object is returned either way: ``load`` mutates the reactions
        in place as it goes, so the entries are readable even after the balance check rejects one,
        which is what lets this test report the count the guard actually saw."""
        database = KineticsDatabase()
        library = KineticsLibrary(label=label)
        path = os.path.join(FIXTURE_ROOT, label, "reactions.py")
        try:
            library.load(path, database.local_context, database.global_context)
        except DatabaseError as exc:
            return library, exc
        return library, None

    def test_declared_electron_count_reaches_the_reaction(self):
        """The count on the rate law must be on the ``Reaction`` by the time the balance check
        runs. Pre-fix this fails reporting ``electrons=0`` -- the reaction was neutralised on the
        way to a guard that then, correctly, rejected it."""
        library, load_error = self._load_library("charged-balanced")
        entry = library.entries[1]
        assert entry.item.electrons == -1, (
            "reaction reached the balance check declaring electrons={0}, while its kinetics "
            "declares electrons={1}; load raised: {2}".format(
                entry.item.electrons, entry.data.electrons.value, load_error
            )
        )
        assert load_error is None, "balanced charged entry was rejected: {0}".format(load_error)

    def test_charged_library_loads_through_the_database_path(self):
        """The same thing through ``KineticsDatabase.load_libraries``, the call the production
        loader makes, so the fix is exercised where the shipped electrochemistry libraries hit
        it."""
        database = KineticsDatabase()
        database.load_libraries(FIXTURE_ROOT, libraries=["charged-balanced"])
        rxns = database.libraries["charged-balanced"].get_library_reactions()
        assert rxns, "no reactions came back from the library"
        assert all(rxn.electrons == -1 for rxn in rxns), [rxn.electrons for rxn in rxns]
        assert all(rxn.is_balanced() for rxn in rxns)

    def test_wrong_electron_count_is_still_rejected(self):
        """Tripwire: the balance check must stay live and un-narrowed. An entry whose declared
        count does not close its charge gap must still be refused, so that the sibling library's
        acceptance means "balanced" and not "excused".

        The rejection is asserted together with the propagated count, because the rejection alone
        is not evidence: pre-fix this library is refused too, but for the unrelated reason that its
        count never arrived. Only ``electrons == -2`` at the point of refusal shows the guard
        weighed the declared count and found it wanting."""
        library, load_error = self._load_library("charged-unbalanced")
        assert library.entries[1].item.electrons == -2, (
            "declared count did not reach the guard; rejection proves nothing"
        )
        assert load_error is not None, "entry off by one electron was accepted"
        assert "not balanced" in str(load_error), load_error
