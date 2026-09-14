#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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
Tests for the two Cd/CO/CS/Cdd atom-type checks in databaseTest.py.

Both checks decide "is this a double bond?" from a :class:`Group` bond's ``order``. Until this
change they did so by comparing ``order[0]`` against the string ``"D"`` -- but a ``Group`` bond
stores ``order`` as a list of *numbers*, so the comparison was never true, ``num_of_d_bonds`` was
always ``0``, and neither check could fail for anything it was ever given. They passed every
kinetics family, thermo group, solvation group, statmech group and transport group in the database,
for free.

These tests drive the real check functions against synthesized groups, so the repaired predicate is
pinned both against a defect it must catch and against a correct group it must not complain about.
They do not load the database, and are fast.
"""

import importlib.util
import os
import types

import pytest

from rmgpy.molecule import Group


def _database_test_class():
    """
    Import ``databaseTest.py`` by path and return ``TestDatabase``.

    A plain relative import would also bring ``TestDatabase`` into this module's namespace, where
    pytest would collect it a second time and pay for its ``setup_class``, which loads the entire
    database. Importing by path under a private module name keeps that cost out of this file.
    """
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "databaseTest.py")
    spec = importlib.util.spec_from_file_location("_databaseTest_for_cdAtomTypeTest", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.TestDatabase


DatabaseChecks = _database_test_class()


#: A carbon double-bonded to an oxygen, typed ``Cd`` but *not* ``CO``. This is precisely the
#: authoring mistake both checks exist to catch.
MISUSED_CD = "1 Cd  u0 {2,D}\n2 O2d u0 {1,D}\n"

#: The same structure, typed correctly.
CORRECT_CO = "1 CO  u0 {2,D}\n2 O2d u0 {1,D}\n"


def make_entry(label, adjlist):
    return types.SimpleNamespace(label=label, item=Group().from_adjacency_list(adjlist), children=[])


def group_database(label, adjlist):
    """The minimum ``general_check_cd_atom_type`` touches: an object with an ``entries`` dict."""
    return types.SimpleNamespace(entries={label: make_entry(label, adjlist)})


def kinetics_harness(label, adjlist, family_name="Cd_Test"):
    """The minimum ``kinetics_check_cd_atom_type`` touches, reached through ``self.database``."""
    node = make_entry(label, adjlist)
    family = types.SimpleNamespace(
        own_reverse=True,  # empties the check's ignore list, so nothing is skipped
        forward_template=types.SimpleNamespace(products=[]),
        groups=types.SimpleNamespace(entries={label: node}),
    )
    harness = types.SimpleNamespace(
        database=types.SimpleNamespace(kinetics=types.SimpleNamespace(families={family_name: family}))
    )
    return harness, family_name


class TestCdAtomTypeCheckIsLive:
    """The repaired predicate must actually fire. Before the repair, every one of these passed."""

    def test_general_check_catches_a_misused_cd(self):
        with pytest.raises(ValueError):
            DatabaseChecks.general_check_cd_atom_type(None, "MisusedCd", group_database("bad", MISUSED_CD))

    def test_general_check_accepts_a_correctly_typed_group(self):
        assert DatabaseChecks.general_check_cd_atom_type(None, "CorrectCO", group_database("good", CORRECT_CO))

    def test_kinetics_check_catches_a_misused_cd(self):
        harness, family_name = kinetics_harness("bad", MISUSED_CD)
        with pytest.raises(ValueError):
            DatabaseChecks.kinetics_check_cd_atom_type(harness, family_name)

    def test_kinetics_check_accepts_a_correctly_typed_group(self):
        harness, family_name = kinetics_harness("good", CORRECT_CO)
        assert DatabaseChecks.kinetics_check_cd_atom_type(harness, family_name)


class TestCdAtomTypeCheckCannotRegress:
    """
    Pin the reason the checks were dead, so the string comparison cannot be reintroduced.

    The first test is a characterization of :class:`Group`, not of the checks: if ``order`` ever
    does become a list of bond-order letters, it fails and names the two checks that read the field
    numerically.
    """

    def test_group_bond_order_is_numeric_not_a_letter(self):
        group = Group().from_adjacency_list(MISUSED_CD)
        orders = [bond.order for atom in group.atoms for bond in atom.bonds.values()]
        assert orders, "the probe group has no bonds"
        for order in orders:
            assert len(order) == 1
            assert not isinstance(order[0], str), (
                "Group bond order became a string. Two checks in databaseTest.py "
                "(kinetics_check_cd_atom_type, general_check_cd_atom_type) count double bonds "
                "numerically and must be revisited."
            )
            assert abs(2 - order[0]) < 1e-7

    def test_the_old_string_comparison_would_have_been_dead(self):
        group = Group().from_adjacency_list(MISUSED_CD)
        atom = group.atoms[0]
        dead = sum([1 if x.order[0] == "D" and len(x.order) == 1 else 0 for x in atom.bonds.values()])
        live = sum([1 if len(x.order) == 1 and abs(2 - x.order[0]) < 1e-7 else 0 for x in atom.bonds.values()])
        assert dead == 0, "the pre-repair predicate is supposed to be dead on a real double bond"
        assert live == 1, "the repaired predicate must see the double bond it is given"
