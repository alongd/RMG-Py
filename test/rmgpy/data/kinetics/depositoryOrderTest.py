#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
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
The depository load order is inclusion order, and it decides nothing.

This file used to pin the opposite: that ``training`` was loaded first and therefore won
every reaction two depositories both carried.  Depository membership is provenance, not a
quality ranking, so that policy is gone and so are the tests that pinned it.  What is
pinned now is narrower and stronger:

1. the order is deterministic and is a property of the code, not of the directory listing;
2. a directory is a depository only if it is an immediate child of the family and holds a
   ``reactions.py`` -- non-depository directories are skipped rather than aborting the load,
   and descendants deeper than one level are not depositories at all; and
3. **the order does not choose a rate.**  Loading the same two depositories in either order
   produces the same refusal, naming the same candidates.  That third one is the load
   order's whole remaining contract, and it is the reason the first can be plain
   alphabetical without anybody having to think about it.

The refusal itself -- what it says and when -- lives in ``depositoryConflictTest.py``.

**Deliberate conflict, for whoever merges this.**  The unmerged branch
``i225-depositories-all`` carries ``test/rmgpy/data/kinetics/depositoryLoadingTest.py``,
whose class ``TestDepositoryOrderIsFilesystemOrder`` asserts that ``training`` comes before
``NIST`` and that the order is not ``sorted()``.  Both of those assertions now fail, on
purpose: the order IS ``sorted()``, and nothing downstream reads it.  When the branches
meet, delete that class rather than restoring the order it describes -- its subject no
longer exists.
"""

import os
import shutil

import pytest

from rmgpy import settings
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import KineticsFamily, order_depository_names
from rmgpy.exceptions import KineticsDepositoryConflictError
from rmgpy.reaction import Reaction
from rmgpy.species import Species

#: A real family small enough to load in a test that nevertheless carries two depository
#: directories on disk, so "order matters" and "order does not matter" are distinguishable.
TWO_DEPOSITORY_FAMILY = "Korcek_step1"

#: A second, also small, two-depository family that carries a reaction BOTH of its
#: depositories supply, in the same direction -- which is what is needed to show that the
#: load order no longer decides a rate.  ``Korcek_step1`` has two depositories but no
#: reaction in both, so it cannot demonstrate anything about selection.
CONTESTED_FAMILY = "HO2_Elimination_from_PeroxyRadical"
CONTESTED_REACTANTS = ['CCCO[O]']
CONTESTED_PRODUCTS = ['C=CC', '[O]O']


def _touch(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write('')


def _listdir_returning(ordering):
    """A drop-in ``os.listdir`` that hands back each directory listing in `ordering` order."""
    real_listdir = os.listdir

    def listdir(path='.'):
        return ordering(real_listdir(path))

    return listdir


class TestDepositoryOrderIsDeterministic:
    """
    The ordering itself, pinned without touching the database.

    These run on a synthetic family directory, so they pin what the code *decides* rather
    than what one particular checkout happens to contain.
    """

    def test_the_order_is_alphabetical(self, tmp_path):
        family = tmp_path / "Some_Family"
        for name in ["zzz", "NIST", "aaa", "training"]:
            _touch(str(family / name / "reactions.py"))

        assert order_depository_names(str(family)) == ["NIST", "aaa", "training", "zzz"]

    def test_no_depository_is_ranked_by_its_name(self, tmp_path):
        """
        ``training`` gets no special treatment; it sits where the alphabet puts it.  Pinned
        as its own assertion so that reintroducing a named priority -- for training, for
        NIST, for anything -- goes red here and names itself, rather than quietly changing
        which source supplies a rate.
        """
        family = tmp_path / "Some_Family"
        for name in ["NIST", "training"]:
            _touch(str(family / name / "reactions.py"))

        names = order_depository_names(str(family))

        assert names == sorted(names), (
            "the depository order is no longer the plain alphabet, which means something "
            "is ranking depositories again. Depository membership is provenance, not a "
            "quality ranking. Got {0}".format(names))

    @pytest.mark.parametrize("ordering", [sorted, lambda d: sorted(d, reverse=True)],
                             ids=["alphabetical", "reverse-alphabetical"])
    def test_order_does_not_depend_on_the_directory_listing(self, tmp_path, monkeypatch, ordering):
        """Hand the loader the directories in either order and get the same answer back."""
        import rmgpy.data.kinetics.family as family_module

        family = tmp_path / "Some_Family"
        for name in ["NIST", "training", "aaa"]:
            _touch(str(family / name / "reactions.py"))

        monkeypatch.setattr(family_module.os, 'listdir', _listdir_returning(ordering))

        assert order_depository_names(str(family)) == ["NIST", "aaa", "training"]


class TestDepositoryDiscoveryIsBounded:
    """
    A depository is an immediate child of the family that holds a ``reactions.py``.

    Two separate failures are pinned here.  A directory with no ``reactions.py`` used to
    abort the entire load-everything run -- the loader built the file path from the family
    directory rather than from the directory it was actually in, and never checked the file
    existed, so one stray ``__pycache__`` or scratch copy was enough.  And the search was
    ``os.walk``, which descends without limit, so any descendant holding a ``reactions.py``
    became a depository: an archive, a backup, a working copy left inside another
    depository.  ``reactions.py`` is executed on load, so that second one is not merely a
    way to pick up rates nobody meant to ship.
    """

    def test_a_directory_without_reactions_is_not_a_depository(self, tmp_path):
        family = tmp_path / "Some_Family"
        _touch(str(family / "training" / "reactions.py"))
        os.makedirs(str(family / "scratch"))
        os.makedirs(str(family / "__pycache__"))

        assert order_depository_names(str(family)) == ["training"]

    def test_a_file_that_is_not_a_directory_is_not_a_depository(self, tmp_path):
        family = tmp_path / "Some_Family"
        _touch(str(family / "training" / "reactions.py"))
        _touch(str(family / "groups.py"))

        assert order_depository_names(str(family)) == ["training"]

    def test_a_depository_nested_below_the_top_level_is_not_found(self, tmp_path):
        """
        The bound.  An archived copy sitting inside another depository holds a perfectly
        valid ``reactions.py``; it is still not a depository of this family, and loading it
        would both supply rates nobody asked for and execute a file nobody reviewed.
        """
        family = tmp_path / "Some_Family"
        _touch(str(family / "training" / "reactions.py"))
        _touch(str(family / "training" / "archive-2019" / "reactions.py"))
        _touch(str(family / "NIST" / "backup" / "old" / "reactions.py"))

        assert order_depository_names(str(family)) == ["training"]


@pytest.mark.database
class TestNonDepositoryDirectoriesDoNotBreakLoading:
    """The same two failures, end to end, against a real family copied out of the database."""

    @staticmethod
    def _copy_family(tmp_path):
        source = os.path.join(settings["database.directory"], "kinetics", "families",
                              TWO_DEPOSITORY_FAMILY)
        target = str(tmp_path / TWO_DEPOSITORY_FAMILY)
        shutil.copytree(source, target)
        return target

    @staticmethod
    def _load(path):
        database = KineticsDatabase()
        family = KineticsFamily(label=TWO_DEPOSITORY_FAMILY)
        family.load(path, database.local_context, database.global_context,
                    depository_labels='all')
        return family

    #: What a clean copy of the family must load, in order.
    EXPECTED = ["{0}/NIST".format(TWO_DEPOSITORY_FAMILY),
                "{0}/training".format(TWO_DEPOSITORY_FAMILY)]

    def test_the_family_under_test_really_has_two_depositories(self, tmp_path):
        """If the database loses this family's NIST directory, the assertions go vacuous."""
        path = self._copy_family(tmp_path)

        on_disk = sorted(d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d)))

        assert on_disk == ["NIST", "training"], (
            "{0} no longer carries exactly training and NIST; pick another family for "
            "these tests. Got {1}".format(TWO_DEPOSITORY_FAMILY, on_disk))

    def test_a_stray_top_level_directory_is_skipped(self, tmp_path):
        path = self._copy_family(tmp_path)
        os.makedirs(os.path.join(path, "scratch"))

        assert [d.label for d in self._load(path).depositories] == self.EXPECTED

    def test_a_directory_nested_inside_a_depository_is_skipped(self, tmp_path):
        path = self._copy_family(tmp_path)
        os.makedirs(os.path.join(path, "training", "notes"))

        assert [d.label for d in self._load(path).depositories] == self.EXPECTED

    def test_a_reactions_file_nested_inside_a_depository_is_not_loaded(self, tmp_path):
        """
        The bound, against real data: a copy of the NIST depository parked inside the
        training depository must not be loaded as a third depository.
        """
        path = self._copy_family(tmp_path)
        shutil.copytree(os.path.join(path, "NIST"), os.path.join(path, "training", "archive"))

        assert [d.label for d in self._load(path).depositories] == self.EXPECTED

    def test_pycache_under_a_family_is_skipped(self, tmp_path):
        """
        Not hypothetical: ``__pycache__`` already exists one level up, in the families
        directory itself, where ``KineticsDatabase.load_families`` discards it by name.
        Nothing stops one appearing under a family.
        """
        path = self._copy_family(tmp_path)
        os.makedirs(os.path.join(path, "__pycache__"))

        assert [d.label for d in self._load(path).depositories] == self.EXPECTED


@pytest.mark.database
class TestTheLoadOrderNoLongerDecidesARate:
    """
    The claim that makes the order safe to leave alphabetical.

    Before this change, reversing ``family.depositories`` changed k(1000 K) for
    ``C3H7O2 <=> C3H6 + HO2`` by a factor of 13 -- see
    ``docs/depository-order/logs/31-red-fixture-silent-choice.stdout.log``, taken on the tip
    that still ranked ``training`` first.  Now both orders refuse, and refuse identically.
    """

    @staticmethod
    def _load():
        database = KineticsDatabase()
        family = KineticsFamily(label=CONTESTED_FAMILY)
        family.load(os.path.join(settings["database.directory"], "kinetics", "families",
                                 CONTESTED_FAMILY),
                    database.local_context, database.global_context, depository_labels='all')
        return family

    @staticmethod
    def _ask(family):
        reaction = Reaction(reactants=[Species(smiles=s) for s in CONTESTED_REACTANTS],
                            products=[Species(smiles=s) for s in CONTESTED_PRODUCTS])
        return family.get_kinetics(
            reaction,
            template_labels=[g.label for g in family.forward_template.reactants],
            degeneracy=1, estimator='', return_all_kinetics=False)

    def test_both_load_orders_refuse_with_the_same_report(self):
        forwards = self._load()
        backwards = self._load()
        backwards.depositories.reverse()
        assert [d.label for d in forwards.depositories] != \
            [d.label for d in backwards.depositories], 'the two orders must differ'

        messages = []
        for family in (forwards, backwards):
            with pytest.raises(KineticsDepositoryConflictError) as caught:
                self._ask(family)
            messages.append(str(caught.value))

        assert messages[0] == messages[1], (
            "the conflict report changed when the load order did, so the load order is "
            "still carrying information it must not carry")

    def test_the_filesystem_listing_does_not_change_the_answer(self, monkeypatch):
        """
        The test that separates "the code decides" from "the filesystem decided and we got
        lucky" -- except that now neither decides, and what is pinned is only that the
        refusal is reached however the directory listing comes back.

        Deliberately says nothing about what the resulting order IS.  That belongs to
        ``TestDepositoryOrderIsDeterministic``, and asserting it here as well would make
        this test fail for a reason that has nothing to do with the refusal -- which is how
        a reader ends up believing the refusal depends on the order.
        """
        import rmgpy.data.kinetics.family as family_module

        monkeypatch.setattr(family_module.os, 'listdir',
                            _listdir_returning(lambda d: sorted(d, reverse=True)))

        with pytest.raises(KineticsDepositoryConflictError):
            self._ask(self._load())
