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
This module contains unit tests of the rmgpy.molecule.atomtype module.
"""


import logging
import rmgpy.molecule
from rmgpy.molecule import Molecule
from rmgpy.molecule.atomtype import get_atomtype

import pytest


class TestAtomType:
    """
    Contains unit tests of the AtomType class.
    """

    def setup_class(self):
        self.atomtype = rmgpy.molecule.atomtype.ATOMTYPES["Cd"]

    def test_pickle(self):
        """
        Test that an AtomType object can be successfully pickled and
        unpickled with no loss of information.
        """
        import pickle

        atom_type = pickle.loads(pickle.dumps(self.atomtype))
        assert self.atomtype.label == atom_type.label
        assert len(self.atomtype.generic) == len(atom_type.generic)
        for item1, item2 in zip(self.atomtype.generic, atom_type.generic):
            assert item1.label == item2.label
        assert len(self.atomtype.specific) == len(atom_type.specific)
        for item1, item2 in zip(self.atomtype.specific, atom_type.specific):
            assert item1.label == item2.label
        assert len(self.atomtype.increment_bond) == len(atom_type.increment_bond)
        for item1, item2 in zip(self.atomtype.increment_bond, atom_type.increment_bond):
            assert item1.label == item2.label
        assert len(self.atomtype.decrement_bond) == len(atom_type.decrement_bond)
        for item1, item2 in zip(self.atomtype.decrement_bond, atom_type.decrement_bond):
            assert item1.label == item2.label
        assert len(self.atomtype.form_bond) == len(atom_type.form_bond)
        for item1, item2 in zip(self.atomtype.form_bond, atom_type.form_bond):
            assert item1.label == item2.label
        assert len(self.atomtype.break_bond) == len(atom_type.break_bond)
        for item1, item2 in zip(self.atomtype.break_bond, atom_type.break_bond):
            assert item1.label == item2.label
        assert len(self.atomtype.increment_radical) == len(atom_type.increment_radical)
        for item1, item2 in zip(self.atomtype.increment_radical, atom_type.increment_radical):
            assert item1.label == item2.label
        assert len(self.atomtype.decrement_radical) == len(atom_type.decrement_radical)
        for item1, item2 in zip(self.atomtype.decrement_radical, atom_type.decrement_radical):
            assert item1.label == item2.label
        for item1, item2 in zip(self.atomtype.increment_charge, atom_type.increment_charge):
            assert item1.label == item2.label
        assert len(self.atomtype.decrement_charge) == len(atom_type.decrement_charge)
        for item1, item2 in zip(self.atomtype.decrement_charge, atom_type.decrement_charge):
            assert item1.label == item2.label

    def test_output(self):
        """
        Test that we can reconstruct an AtomType object from its repr()
        with no loss of information.
        """
        namespace = {}
        exec(
            "atomtype = rmgpy.molecule.atomtype.ATOMTYPES[{0!r}]".format(self.atomtype.__repr__().split('"')[1]),
            globals(),
            namespace,
        )
        assert "atomtype" in namespace
        atomtype = namespace["atomtype"]
        assert self.atomtype.equivalent(atomtype)

    def test_equivalent(self):
        """
        Test the AtomType.equivalent() method.
        """
        assert self.atomtype.equivalent(rmgpy.molecule.atomtype.ATOMTYPES["Cd"])

    def test_is_specfic_case_of(self):
        """
        Test the AtomType.is_specific_case_of() method.
        """
        assert self.atomtype.is_specific_case_of(rmgpy.molecule.atomtype.ATOMTYPES["C"])

    def test_set_actions(self):
        """
        Test the AtomType.set_actions() method.
        """
        other = rmgpy.molecule.atomtype.AtomType("Test", generic=["R"], specific=[])
        other.set_actions(
            self.atomtype.increment_bond,
            self.atomtype.decrement_bond,
            self.atomtype.form_bond,
            self.atomtype.break_bond,
            self.atomtype.increment_radical,
            self.atomtype.decrement_radical,
            self.atomtype.increment_lone_pair,
            self.atomtype.decrement_lone_pair,
            self.atomtype.increment_charge,
            self.atomtype.decrement_charge,
        )
        assert self.atomtype.increment_bond == other.increment_bond
        assert self.atomtype.decrement_bond == other.decrement_bond
        assert self.atomtype.form_bond == other.form_bond
        assert self.atomtype.break_bond == other.break_bond
        assert self.atomtype.increment_radical == other.increment_radical
        assert self.atomtype.decrement_radical == other.decrement_radical
        assert self.atomtype.increment_charge == other.increment_charge
        assert self.atomtype.decrement_charge == other.decrement_charge

    """
    Currently RMG doesn't even detect aromaticity of furan or thiophene, so making
    them out of atom type group definitions is not an easy fix. Even if we did make a
    sample molecule out of it, we'd probably not realize it was right. Currently
    it tries to make a 6-membered ring (not realizing it should be 5-membered)
    and then gets the wrong number of pi electrons.
    """
    EXPECTED_FAILING_ATOMTYPES = ["O4b", "S4b"]

    def test_make_sample_molecule(self):
        """
        Test we can make a sample molecule for most atom type without crashing.
        """
        failed = []
        for name, atom_type in rmgpy.molecule.atomtype.ATOMTYPES.items():
            if name in self.EXPECTED_FAILING_ATOMTYPES:
                continue  # These are known to fail. See next WIP test
            adjlist = f"1 {name} ux"
            group = rmgpy.molecule.Group().from_adjacency_list(adjlist)
            try:
                result = group.make_sample_molecule()
                # logging.info(f"For {name} made\n{result.to_adjacency_list()}")
            except:
                logging.exception(f"Couldn't make sample molecule for atomType {name}")
                failed.append(name)
        assert len(failed) == 0, f"Couldn't make sample molecules for types {', '.join(failed)}"

    @pytest.mark.skip(reason="WIP")
    def test_make_sample_molecule_wip(self):
        """
        Test we can make a sample molecule for some failing atom types.
        """
        failed = []
        for name in self.EXPECTED_FAILING_ATOMTYPES:
            adjlist = f"1 {name} ux"
            group = rmgpy.molecule.Group().from_adjacency_list(adjlist)
            try:
                result = group.make_sample_molecule()
            except:
                logging.exception(f"Couldn't make sample molecule for atomType {name}")
                failed.append(name)
        assert not failed, f"Couldn't make sample molecules for types {', '.join(failed)}"

    @pytest.mark.skip(reason="WIP")
    def test_make_sample_molecule_right(self):
        """
        Test we can make the correct sample molecule for each atom type.
        """
        failed = []
        for name, atom_type in rmgpy.molecule.atomtype.ATOMTYPES.items():
            if name in self.EXPECTED_FAILING_ATOMTYPES:
                continue  # These are known to fail. See next WIP test
            adjlist = f"1 {name} ux"
            group = rmgpy.molecule.Group().from_adjacency_list(adjlist)
            try:
                result = group.make_sample_molecule()
                if not result.is_subgraph_isomorphic(group):
                    failed.append(name)
                    logging.error(
                        f"Sample molecule for {name} is not correct. "
                        f"Expected:\n{group.to_adjacency_list().strip()}\n"
                        f"Got:\n{result.to_adjacency_list().strip()}"
                    )
            except:
                logging.exception(f"Couldn't make sample molecule for atomType {name}")
                failed.append(name)
        assert not failed, f"Couldn't make correct sample molecules for types {', '.join(failed)}"


class TestGetAtomType:
    """
    Contains unit tests of the get_atomtype() method.
    """

    def setup_class(self):
        self.mol1 = Molecule().from_smiles("COC(=O)CC=C=CC#C")
        # self.mol2 = Molecule().from_smiles('c1ccccc1')
        # the from_smiles method currently Kekulizes, so to test Benzene we use from_adjacency_list
        self.mol2 = Molecule().from_adjacency_list(
            """1  C u0 p0 {2,B} {6,B} {7,S}
                                                    2  C u0 p0 {1,B} {3,B} {8,S}
                                                    3  C u0 p0 {2,B} {4,B} {9,S}
                                                    4  C u0 p0 {3,B} {5,B} {10,S}
                                                    5  C u0 p0 {4,B} {6,B} {11,S}
                                                    6  C u0 p0 {1,B} {5,B} {12,S}
                                                    7  H u0 p0 {1,S}
                                                    8  H u0 p0 {2,S}
                                                    9  H u0 p0 {3,S}
                                                    10 H u0 p0 {4,S}
                                                    11 H u0 p0 {5,S}
                                                    12 H u0 p0 {6,S}"""
        )
        self.mol3 = Molecule().from_smiles("[H]")
        self.mol4 = Molecule().from_smiles("O=[Si][Si][Si]=[Si]=[Si][Si]#[Si]SS=S")
        self.mol5 = Molecule().from_adjacency_list(
            """1 H u0 p0 {3,S}
                                                    2 H u0 p0 {3,S}
                                                    3 N u0 p0 c+1 {1,S} {2,S} {4,D}
                                                    4 N u0 p2 c-1 {3,D}"""
        )
        self.mol6 = Molecule().from_smiles("[Ar]")
        self.mol7 = Molecule().from_smiles("[He]")
        self.mol8 = Molecule().from_smiles("[Ne]")
        self.mol9 = Molecule().from_adjacency_list(
            """1 N u0 p1 {2,S} {3,S} {4,S}
                                                    2 H u0 p0 {1,S}
                                                    3 H u0 p0 {1,S}
                                                    4 H u0 p0 {1,S}"""
        )

        self.mol10 = Molecule().from_adjacency_list(
            """1 N u1 p1 {2,S} {3,S}
                                                     2 H u0 p0 {1,S}
                                                     3 H u0 p0 {1,S}"""
        )

        self.mol11 = Molecule().from_adjacency_list(
            """1 N u2 p1 {2,S}
                                                     2 H u0 p0 {1,S}"""
        )

        self.mol12 = Molecule().from_adjacency_list(
            """1 N u0 p1 {2,T}
                                                     2 C u1 p0 {1,T}"""
        )

        self.mol14 = Molecule().from_adjacency_list(
            """1 N u0 p2 c-1 {2,D}
                                                     2 N u0 p0 c+1 {1,D} {3,D}
                                                     3 O u0 p2 {2,D}"""
        )

        self.mol15 = Molecule().from_adjacency_list(
            """1 N u0 p1 c0 {2,T}
                                                     2 N u0 p0 c+1 {1,T} {3,S}
                                                     3 O u0 p3 c-1 {2,S}"""
        )

        self.mol16 = Molecule().from_adjacency_list(
            """1 N u0 p1 {2,D} {3,S}
                                                     2 O u0 p2 {1,D}
                                                     3 O u1 p2 {1,S}"""
        )

        self.mol17 = Molecule().from_adjacency_list(
            """1 N u1 p1 {2,D}
                                                     2 O u0 p2 {1,D}"""
        )

        self.mol18 = Molecule().from_adjacency_list(
            """1  N u0 p0 c+1 {2,B} {6,B} {7,S}
                                                     2  C u0 p0 {1,B} {3,B} {8,S}
                                                     3  C u0 p0 {2,B} {4,B} {9,S}
                                                     4  C u0 p0 {3,B} {5,B} {10,S}
                                                     5  C u0 p0 {4,B} {6,B} {11,S}
                                                     6  N u0 p1 {1,B} {5,B}
                                                     7  O u0 p3 c-1 {1,S}
                                                     8  H u0 p0 {2,S}
                                                     9  H u0 p0 {3,S}
                                                     10 H u0 p0 {4,S}
                                                     11 H u0 p0 {5,S}"""
        )

        self.mol19 = Molecule().from_smiles("C=S")

        self.mol20 = Molecule().from_smiles("[C-]#[O+]")

        self.mol21 = Molecule().from_adjacency_list(
            """1 S u0 p3 c-1 {2,S}
                                                     2 S u0 p2 c+1 {1,S}"""
        )

        self.mol22 = Molecule().from_adjacency_list("""1 S u0 p3 c0""")

        self.mol23 = Molecule().from_adjacency_list(
            """1 S u0 p2 c0 {2,S} {5,S}
                                                     2 S u0 p1 c+1 {1,S} {3,S} {4,S}
                                                     3 C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
                                                     4 O u0 p3 c-1 {2,S}
                                                     5 H u0 p0 c0 {1,S}
                                                     6 H u0 p0 c0 {3,S}
                                                     7 H u0 p0 c0 {3,S}
                                                     8 H u0 p0 c0 {3,S}"""
        )

        self.mol24 = Molecule().from_adjacency_list(
            """1 C u0 p0 c0 {2,D} {4,S} {5,S}
                                                     2 S u0 p2 c-1 {1,D} {3,S}
                                                     3 O u0 p2 c+1 {2,S}
                                                     4 H u0 p0 c0 {1,S}
                                                     5 H u0 p0 c0 {1,S}"""
        )

        self.mol25 = Molecule().from_adjacency_list(
            """1 S u0 p1 c0 {2,S} {5,S} {7,S} {8,S}
                                                     2 O u0 p2 c0 {1,S} {3,S}
                                                     3 S u0 p1 c0 {2,S} {4,S} {9,D}
                                                     4 O u0 p2 c0 {3,S} {6,S}
                                                     5 H u0 p0 c0 {1,S}
                                                     6 H u0 p0 c0 {4,S}
                                                     7 H u0 p0 c0 {1,S}
                                                     8 H u0 p0 c0 {1,S}
                                                     9 O u0 p2 c0 {3,D}"""
        )

        self.mol26 = Molecule().from_adjacency_list(
            """1 O u0 p3 c-1 {2,S}
                                                     2 S u0 p1 c+1 {1,S} {3,D}
                                                     3 O u0 p2 c0 {2,D}"""
        )

        # self.mol27 = Molecule().from_adjacency_list('''1 S u0 p1 c0 {2,B} {5,B}
        #                                              2 C u0 p0 c0 {1,B} {3,B} {6,S}
        #                                              3 C u0 p0 c0 {2,B} {4,B} {7,S}
        #                                              4 C u0 p0 c0 {3,B} {5,B} {8,S}
        #                                              5 C u0 p0 c0 {1,B} {4,B} {9,S}
        #                                              6 H u0 p0 c0 {2,S}
        #                                              7 H u0 p0 c0 {3,S}
        #                                              8 H u0 p0 c0 {4,S}
        #                                              9 H u0 p0 c0 {5,S}''')

        self.mol28 = Molecule().from_adjacency_list(
            """1  O u0 p2 c0 {2,D}
                                                     2  S u0 p1 c0 {1,D} {3,D}
                                                     3  C u0 p0 c0 {2,D} {4,S} {7,S}
                                                     4  C u0 p0 c0 {3,S} {5,T}
                                                     5  S u0 p1 c0 {4,T} {6,S}
                                                     6  S u0 p0 c0 {5,S} {8,S} {9,S} {10,S} {11,S} {12,S}
                                                     7  H u0 p0 c0 {3,S}
                                                     8  H u0 p0 c0 {6,S}
                                                     9  H u0 p0 c0 {6,S}
                                                     10 H u0 p0 c0 {6,S}
                                                     11 H u0 p0 c0 {6,S}
                                                     12 H u0 p0 c0 {6,S}"""
        )

        self.mol29 = Molecule().from_adjacency_list(
            """1 C u0 p1 c-1 {2,T}
                                                     2 S u0 p1 c+1 {1,T}"""
        )

        self.mol30 = Molecule().from_adjacency_list(
            """1 S u0 p0 c0 {2,D} {3,S} {4,S} {5,S} {6,S}
                                                     2 O u0 p2 c0 {1,D}
                                                     3 H u0 p0 c0 {1,S}
                                                     4 H u0 p0 c0 {1,S}
                                                     5 H u0 p0 c0 {1,S}
                                                     6 H u0 p0 c0 {1,S}"""
        )

        self.mol31 = Molecule().from_adjacency_list(
            """1 S u0 p0 c+1 {2,S} {3,D} {4,D}
                                                     2 O u0 p3 c-1 {1,S}
                                                     3 O u0 p2 c0 {1,D}
                                                     4 O u0 p2 c0 {1,D}"""
        )

        self.mol32 = Molecule().from_adjacency_list(
            """1 O u0 p2 c0 {2,D}
                                                     2 S u0 p0 c0 {1,D} {3,D} {4,S} {5,S}
                                                     3 O u0 p2 c0 {2,D}
                                                     4 O u0 p2 c0 {2,S} {6,S}
                                                     5 O u0 p2 c0 {2,S} {7,S}
                                                     6 H u0 p0 c0 {4,S}
                                                     7 H u0 p0 c0 {5,S}"""
        )

        self.mol33 = Molecule().from_adjacency_list(
            """1 O u0 p3 c-1 {2,S}
                                                     2 S u0 p0 c+1 {1,S} {3,D} {4,D}
                                                     3 O u0 p2 c0 {2,D}
                                                     4 O u0 p2 c0 {2,D}"""
        )

        self.mol34 = Molecule().from_adjacency_list(
            """1 O u0 p2 c0 {2,D}
                                                     2 S u0 p0 c0 {1,D} {3,D} {4,D}
                                                     3 O u0 p2 c0 {2,D}
                                                     4 O u0 p2 c0 {2,D}"""
        )

        self.mol35 = Molecule().from_adjacency_list(
            """1 S u0 p0 c0 {2,T} {3,S} {4,S} {5,S}
                                                     2 N u0 p1 c0 {1,T}
                                                     3 H u0 p0 c0 {1,S}
                                                     4 H u0 p0 c0 {1,S}
                                                     5 H u0 p0 c0 {1,S}"""
        )

        self.mol36 = Molecule().from_adjacency_list(
            """1 S u0 p0 c0 {2,T} {3,D} {4,S}
                                                     2 N u0 p1 c0 {1,T}
                                                     3 O u0 p2 c0 {1,D}
                                                     4 H u0 p0 c0 {1,S}"""
        )

        self.mol37 = Molecule().from_adjacency_list(
            """1 N u0 p1 c0 {2,T}
                                                     2 S u0 p0 c0 {1,T} {3,T}
                                                     3 N u0 p1 c0 {2,T}"""
        )

        self.mol38 = Molecule().from_smiles("O=S=O")

        self.mol39 = Molecule().from_adjacency_list(
            """1 N u0 p2 c-1 {2,S} {3,S}
                                                     2 H u0 p0 c0 {1,S}
                                                     3 N u0 p0 c+1 {1,S} {4,T}
                                                     4 C u0 p0 c0 {3,T} {5,S}
                                                     5 H u0 p0 c0 {4,S}"""
        )

        self.mol40 = Molecule().from_adjacency_list(
            """1 N u0 p0 c+1 {2,S} {3,T}
                                                     2 H u0 p0 c0 {1,S}
                                                     3 N u0 p0 c+1 {1,T} {4,S}
                                                     4 N u0 p3 c-2 {3,S}"""
        )

        self.mol41 = Molecule().from_adjacency_list(
            """1 N u0 p2 c0 {2,S}
                                                     2 H u0 p0 c0 {1,S}"""
        )

        self.mol42 = Molecule().from_adjacency_list(
            """1 N u0 p1 c0 {2,T}
                                                     2 N u0 p0 c+1 {1,T} {3,S}
                                                     3 S u0 p2 c-1 {2,S} {4,S} {5,S}
                                                     4 O u1 p2 c0 {3,S}
                                                     5 O u1 p2 c0 {3,S}"""
        )

        self.mol43 = Molecule().from_adjacency_list(
            """1 C u0 p1 c-1 {2,D} {3,S}
                                                     2 S u1 p0 c+1 {1,D} {4,S} {5,S}
                                                     3 H u0 p0 c0 {1,S}
                                                     4 H u0 p0 c0 {2,S}
                                                     5 H u0 p0 c0 {2,S}"""
        )

        self.mol44 = Molecule().from_adjacency_list("""1 O u0 p3 c0""")

        self.mol45 = Molecule().from_adjacency_list(
            """1 O u0 p2 c0 {2,S} {5,S}
                                                     2 N u0 p0 c+1 {1,S} {3,S} {4,D}
                                                     3 O u0 p3 c-1 {2,S}
                                                     4 O u0 p2 c0 {2,D}
                                                     5 H u0 p0 c0 {1,S}"""
        )

        self.mol49 = Molecule().from_adjacency_list(
            """1 O u0 p3 c-1 {2,S}
                                                     2 O u0 p1 c+1 {1,S} {3,S} {4,S}
                                                     3 H u0 p0 c0 {2,S}
                                                     4 S u0 p2 c0 {2,S} {5,S}
                                                     5 H u0 p0 c0 {4,S}"""
        )

        self.mol50 = Molecule().from_adjacency_list(
            """1 O u0 p3 c-1 {2,S}
                                                     2 O u0 p1 c+1 {1,S} {3,D}
                                                     3 C u0 p0 c0 {2,D} {4,S} {5,S}
                                                     4 H u0 p0 c0 {3,S}
                                                     5 H u0 p0 c0 {3,S}"""
        )

        self.mol51 = Molecule().from_adjacency_list(
            """1 O u0 p2 c0 {2,S} {7,S}
                                                     2 S u0 p0 c+1 {1,S} {3,S} {4,S} {5,S} {6,S}
                                                     3 H u0 p0 c0 {2,S}
                                                     4 H u0 p0 c0 {2,S}
                                                     5 H u0 p0 c0 {2,S}
                                                     6 O u0 p3 c-1 {2,S}
                                                     7 H u0 p0 c0 {1,S}"""
        )

        self.mol52 = Molecule().from_adjacency_list(
            """1  C u0 p0 c0 {2,D} {6,S} {8,S}
                                                     2  C u0 p0 c0 {1,D} {3,S} {9,S}
                                                     3  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
                                                     4  C u0 p0 c0 {3,S} {5,S} {6,S} {12,S}
                                                     5  O u0 p3 c-1 {4,S}
                                                     6  C u0 p0 c+1 {1,S} {4,S} {7,S}
                                                     7  H u0 p0 c0 {6,S}
                                                     8  H u0 p0 c0 {1,S}
                                                     9  H u0 p0 c0 {2,S}
                                                     10 H u0 p0 c0 {3,S}
                                                     11 H u0 p0 c0 {3,S}
                                                     12 H u0 p0 c0 {4,S}"""
        )

        self.mol53 = Molecule().from_adjacency_list(
            """1 N u0 p0 c-1 {2,D} {3,D} {4,D}
                                                     2 C u0 p0 c0 {1,D} {5,S} {6,S}
                                                     3 C u0 p0 c0 {1,D} {7,S} {8,S}
                                                     4 N u0 p0 c+1 {1,D} {9,S} {10,S}
                                                     5 H u0 p0 c0 {2,S}
                                                     6 H u0 p0 c0 {2,S}
                                                     7 H u0 p0 c0 {3,S}
                                                     8 H u0 p0 c0 {3,S}
                                                     9 H u0 p0 c0 {4,S}
                                                     10 H u0 p0 c0 {4,S}"""
        )

        self.mol54 = Molecule().from_adjacency_list(
            """1 C u0 p0 c+1 {2,S} {3,D}
                                                     2 H u0 p0 c0 {1,S}
                                                     3 C u0 p0 c0 {1,D} {4,D}
                                                     4 C u0 p1 c-1 {3,D} {5,S}
                                                     5 H u0 p0 c0 {4,S}"""
        )

        self.mol55 = Molecule().from_adjacency_list(
            """1  C u0 p0 c0 {2,B} {10,B} {11,S}
                                                     2  C u0 p0 c0 {1,B} {3,B} {12,S}
                                                     3  C u0 p0 c0 {2,B} {4,B} {13,S}
                                                     4  C u0 p0 c0 {3,B} {5,B} {9,B}
                                                     5  C u0 p0 c0 {4,B} {6,B} {14,S}
                                                     6  C u0 p0 c0 {5,B} {7,B} {15,S}
                                                     7  C u0 p0 c0 {6,B} {8,B} {16,S}
                                                     8  C u0 p0 c0 {7,B} {9,B} {17,S}
                                                     9  C u0 p0 c0 {4,B} {8,B} {10,B}
                                                     10 C u0 p0 c0 {1,B} {9,B} {18,S}
                                                     11 H u0 p0 c0 {1,S}
                                                     12 H u0 p0 c0 {2,S}
                                                     13 H u0 p0 c0 {3,S}
                                                     14 H u0 p0 c0 {5,S}
                                                     15 H u0 p0 c0 {6,S}
                                                     16 H u0 p0 c0 {7,S}
                                                     17 H u0 p0 c0 {8,S}
                                                     18 H u0 p0 c0 {10,S}"""
        )

        self.mol56 = Molecule().from_adjacency_list(
            """1 C u0 p1 c0 {2,S} {3,S}
                                                     2 H u0 p0 c0 {1,S}
                                                     3 H u0 p0 c0 {1,S}"""
        )

        self.mol57 = Molecule().from_adjacency_list(
            """1 C u0 p1 c-1 {2,S} {3,S} {4,S}
                                                     2 H u0 p0 c0 {1,S}
                                                     3 H u0 p0 c0 {1,S}
                                                     4 N u0 p0 c+1 {1,S} {5,T}
                                                     5 N u0 p1 c0 {4,T}"""
        )

        self.mol58 = Molecule().from_adjacency_list(
            """1 C u0 p1 c0 {2,D}
                                                     2 C u0 p0 c0 {1,D} {3,S} {4,S}
                                                     3 H u0 p0 c0 {2,S}
                                                     4 H u0 p0 c0 {2,S}"""
        )

        self.mol59 = Molecule().from_adjacency_list(
            """1 C u0 p1 c-1 {2,S} {3,D}
                                                     2 H u0 p0 c0 {1,S}
                                                     3 N u0 p0 c+1 {1,D} {4,D}
                                                     4 O u0 p2 c0 {3,D}"""
        )

        self.mol60 = Molecule().from_adjacency_list(
            """1 C u0 p0 c0 {2,D} {3,D}
                                                     2 C u0 p0 c+1 {1,D} {4,S}
                                                     3 C u0 p1 c-1 {1,D} {5,S}
                                                     4 H u0 p0 c0 {2,S}
                                                     5 H u0 p0 c0 {3,S}"""
        )

        self.mol64 = Molecule().from_adjacency_list(
            """1 N u0 p1 c0 {2,D} {4,S}
                                                     2 N u0 p0 c+1 {1,D} {3,D}
                                                     3 N u0 p2 c-1 {2,D}
                                                     4 H u0 p0 c0 {1,S}"""
        )

        self.mol69 = Molecule().from_adjacency_list(
            """1 N u0 p0 c+1 {2,T} {3,S}
                                                     2 S u0 p2 c-1 {1,T}
                                                     3 H u0 p0 c0 {1,S}"""
        )

        self.mol70 = Molecule().from_adjacency_list(
            """1 S u0 p0 c+1 {2,D} {3,T}
                                                     2 N u0 p2 c-1 {1,D}
                                                     3 N u0 p1 c0 {1,T}"""
        )

        # self.mol71 = Molecule().from_adjacency_list('''1 O u0 p1 c0 {2,B} {5,B}
        #                                              2 C u0 p0 c0 {1,B} {3,B} {6,S}
        #                                              3 C u0 p0 c0 {2,B} {4,B} {7,S}
        #                                              4 C u0 p0 c0 {3,B} {5,B} {8,S}
        #                                              5 C u0 p0 c0 {1,B} {4,B} {9,S}
        #                                              6 H u0 p0 c0 {2,S}
        #                                              7 H u0 p0 c0 {3,S}
        #                                              8 H u0 p0 c0 {4,S}
        #                                              9 H u0 p0 c0 {5,S}''')

        # self.mol72 = Molecule().from_adjacency_list('''1  N u0 p0 c0 {2,B} {6,B} {7,D}
        #                                              2  C u0 p0 {1,B} {3,B} {8,S}
        #                                              3  C u0 p0 {2,B} {4,B} {9,S}
        #                                              4  C u0 p0 {3,B} {5,B} {10,S}
        #                                              5  C u0 p0 {4,B} {6,B} {11,S}
        #                                              6  N u0 p1 {1,B} {5,B}
        #                                              7  O u0 p2 c0 {1,D}
        #                                              8  H u0 p0 {2,S}
        #                                              9  H u0 p0 {3,S}
        #                                              10 H u0 p0 {4,S}
        #                                              11 H u0 p0 {5,S}''')

        self.mol73 = Molecule().from_adjacency_list(
            """1 H  u0 p0 c0 {2,S}
                                                       2 Cl u0 p3 c0 {1,S}"""
        )

        self.mol74 = Molecule().from_adjacency_list(
            """1 H u0 p0 c0 {2,S}
                                                       2 I u0 p3 c0 {1,S}"""
        )

        self.mol75 = Molecule().from_adjacency_list(
            """1 H u0 p0 c0 {2,S}
                                                       2 F u0 p3 c0 {1,S}"""
        )

        self.mol75_anion = Molecule().from_adjacency_list("""1 F u0 p4 c-1""")

        self.mol75_cation = Molecule().from_adjacency_list("""1 F u0 p3 c+1""")

        self.mol_hydride = Molecule().from_adjacency_list("""1 H u0 p1 c-1""")

        self.mol_chloride = Molecule().from_adjacency_list("""1 Cl u0 p4 c-1""")

        self.mol_bromide = Molecule().from_adjacency_list("""1 Br u0 p4 c-1""")

        # SF6- is a doublet radical anion: 49 valence electrons, so the odd electron sits on S
        self.mol_sf6_anion = Molecule().from_adjacency_list(
            """multiplicity 2
               1 S u1 p0 c-1 {2,S} {3,S} {4,S} {5,S} {6,S} {7,S}
               2 F u0 p3 c0 {1,S}
               3 F u0 p3 c0 {1,S}
               4 F u0 p3 c0 {1,S}
               5 F u0 p3 c0 {1,S}
               6 F u0 p3 c0 {1,S}
               7 F u0 p3 c0 {1,S}"""
        )

        self.mol76 = Molecule().from_adjacency_list(
            """1 H u0 p0 c0 {2,S}
                                                       2 X u0 p0 c0 {1,S}"""
        )

        self.mol77 = Molecule().from_adjacency_list(
            """1 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
                                                       2 H u0 p0 c0 {1,S}
                                                       3 H u0 p0 c0 {1,S}
                                                       4 X u0 p0 c0
                                                       5 H u0 p0 c0 {1,S}
                                                       6 H u0 p0 c0 {1,S}"""
        )

        self.mol78 = Molecule().from_adjacency_list("""1 X u0 p0 c0""")

        self.mol79 = Molecule().from_adjacency_list(
            """1 H  u0 p0 c0 {2,S}
                                                       2 Br u0 p3 c0 {1,S}"""
        )

        self.mol80 = Molecule().from_adjacency_list(
            """1 P u0 p3 c-2 {2,S}
                                                       2 P u0 p0 c+1 {1,S} {3,T}
                                                       3 P u0 p0 c+1 {2,T} {4,S}
                                                       4 H u0 p0 c0 {3,S}"""
        )

        self.mol81 = Molecule().from_adjacency_list(
            """1 P u0 p2 c0 {2,S}
                                                       2 H u0 p0 c0 {1,S}"""
        )

        self.mol82 = Molecule().from_adjacency_list(
            """1 P u0 p2 c-1 {2,S} {3,S}
                                                       2 H u0 p0 c0 {1,S}
                                                       3 P u0 p0 c+1 {1,S} {4,T}
                                                       4 C u0 p0 c0 {3,T} {5,S}
                                                       5 H u0 p0 c0 {4,S}"""
        )

        self.mol83 = Molecule().from_adjacency_list(
            """1 H u0 p0 {3,S}
                                                       2 H u0 p0 {3,S}
                                                       3 P u0 p0 c+1 {1,S} {2,S} {4,D}
                                                       4 P u0 p2 c-1 {3,D}"""
        )

        self.mol84 = Molecule().from_adjacency_list(
            """1 P u0 p1 c0 {4,S} {7,S} {8,S}
                                                       2 P u0 p1 c0 {3,D} {4,S}
                                                       3 O u0 p2 c0 {2,D}
                                                       4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
                                                       5 H u0 p0 c0 {4,S}
                                                       6 H u0 p0 c0 {4,S}
                                                       7 H u0 p0 c0 {1,S}
                                                       8 H u0 p0 c0 {1,S}"""
        )

        self.mol85 = Molecule().from_adjacency_list(
            """1 P u0 p1 c0 {2,T}
                                                       2 C u0 p0 c0 {1,T} {3,S}
                                                       3 H u0 p0 c0 {2,S}"""
        )

        self.mol86 = Molecule().from_adjacency_list(
            """1  P u0 p0 c+1 {2,B} {6,B} {7,S}
                                                       2  C u0 p0 {1,B} {3,B} {8,S}
                                                       3  C u0 p0 {2,B} {4,B} {9,S}
                                                       4  C u0 p0 {3,B} {5,B} {10,S}
                                                       5  C u0 p0 {4,B} {6,B} {11,S}
                                                       6  P u0 p1 {1,B} {5,B}
                                                       7  O u0 p3 c-1 {1,S}
                                                       8  H u0 p0 {2,S}
                                                       9  H u0 p0 {3,S}
                                                       10 H u0 p0 {4,S}
                                                       11 H u0 p0 {5,S}"""
        )

        self.mol87 = Molecule().from_adjacency_list(
            """1 P  u0 p0 c0 {2,S} {3,S} {4,S} {5,S} {6,S}
                                                       2 Cl u0 p3 c0 {1,S}
                                                       3 Cl u0 p3 c0 {1,S}
                                                       4 Cl u0 p3 c0 {1,S}
                                                       5 Cl u0 p3 c0 {1,S}
                                                       6 Cl u0 p3 c0 {1,S}"""
        )

        self.mol88 = Molecule().from_adjacency_list(
            """1 P u0 p0 c+1 {2,S} {3,S} {4,S} {5,S}
                                                       2 O u0 p2 c0 {1,S} {6,S}
                                                       3 O u0 p3 c-1 {1,S}
                                                       4 H u0 p0 c0 {1,S}
                                                       5 H u0 p0 c0 {1,S}
                                                       6 H u0 p0 c0 {2,S}"""
        )

        self.mol89 = Molecule().from_adjacency_list(
            """1 P u0 p0 c0 {2,S} {3,S} {4,S} {5,D}
                                                       2 O u0 p2 c0 {1,S} {6,S}
                                                       3 O u0 p2 c0 {1,S} {7,S}
                                                       4 O u0 p2 c0 {1,S} {8,S}
                                                       5 O u0 p2 c0 {1,D}
                                                       6 H u0 p0 c0 {2,S}
                                                       7 H u0 p0 c0 {3,S}
                                                       8 H u0 p0 c0 {4,S}"""
        )

        self.mol90 = Molecule().from_adjacency_list(
            """1 P u0 p0 c0 {2,D} {3,D} {4,S}
                                                       2 O u0 p2 c0 {1,D}
                                                       3 O u0 p2 c0 {1,D}
                                                       4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
                                                       5 H u0 p0 c0 {4,S}
                                                       6 H u0 p0 c0 {4,S}
                                                       7 H u0 p0 c0 {4,S}"""
        )

        self.mol91 = Molecule().from_adjacency_list(
            """1 P u0 p0 c+1 {2,D} {3,D}
                                                       2 N u0 p2 c-1 {1,D}
                                                       3 C u0 p0 c0 {1,D} {4,S} {5,S}
                                                       4 H u0 p0 c0 {3,S}
                                                       5 H u0 p0 c0 {3,S}"""
        )

        self.mol92 = Molecule().from_adjacency_list(
            """1 P  u0 p0 c0 {2,T} {3,S} {4,S}
                                                       2 N  u0 p1 c0 {1,T}
                                                       3 Cl u0 p3 c0 {1,S}
                                                       4 Cl u0 p3 c0 {1,S}"""
        )

        self.mol93 = Molecule().from_adjacency_list(
            """1 P u0 p0 c0 {2,D} {3,T}
                                                       2 O u0 p2 c0 {1,D}
                                                       3 C u0 p0 c0 {1,T} {4,S}
                                                       4 H u0 p0 c0 {3,S}"""
        )

        self.mol94 = Molecule().from_adjacency_list(
            """1  P u0 p0 {2,B} {6,B} {7,D}
                                                       2  C u0 p0 {1,B} {3,B} {8,S}
                                                       3  C u0 p0 {2,B} {4,B} {9,S}
                                                       4  C u0 p0 {3,B} {5,B} {10,S}
                                                       5  C u0 p0 {4,B} {6,B} {11,S}
                                                       6  C u0 p0 {1,B} {5,B} {12,S}
                                                       7  S u0 p2 {1,D}
                                                       8  H u0 p0 {2,S}
                                                       9  H u0 p0 {3,S}
                                                       10 H u0 p0 {4,S}
                                                       11 H u0 p0 {5,S}
                                                       12 H u0 p0 {6,S}"""
        )

        self.mol95 = Molecule().from_adjacency_list(
            """1 C u0 p0 c+1 {2,T}
                                                       2 C u0 p1 c-1 {1,T}"""
        )

        self.mol96 = Molecule().from_adjacency_list('''1 Li u0 p0 c+1''')
        self.mol97 = Molecule().from_adjacency_list('''1 Li u1 p0 c0''')

        self.electron = Molecule().from_adjacency_list('''1 e u1 p0 c-1''')
        self.proton = Molecule().from_adjacency_list('''1 H u0 p0 c+1''')

    def atom_type(self, mol, atom_id):
        atom = mol.atoms[atom_id]
        atom_type = get_atomtype(atom, mol.get_bonds(atom))
        if atom_type is None:
            return atom_type
        else:
            return atom_type.label

    def test_hydrogen_type(self):
        """
        Test that get_atomtype() returns the hydrogen atom type.
        """
        assert self.atom_type(self.mol3, 0) == "H0"

    def test_hydride_type(self):
        """
        Test that get_atomtype() resolves the hydride anion.
        """
        assert self.atom_type(self.mol_hydride, 0) == "H-"
        assert self.mol_hydride.get_net_charge() == -1

    def test_hydride_survives_update(self):
        """
        Test that Molecule.update() leaves the hydride anion at -1.

        update_lone_pairs() used to force lone_pairs = 0 on every hydrogen regardless of charge,
        which stripped the anion back to neutral; update_charge() then re-derived +1 from the
        electron count. The inversion was silent, so this asserts on the charge, not on a raise.
        """
        mol = Molecule().from_adjacency_list("""1 H u0 p1 c-1""")
        mol.update()
        assert mol.get_net_charge() == -1
        assert mol.atoms[0].lone_pairs == 1
        assert mol.atoms[0].atomtype.label == "H-"

    def test_carbon_types(self):
        """
        Test that get_atomtype() returns appropriate carbon atom types.
        """
        assert self.atom_type(self.mol1, 0) == "Cs"
        assert self.atom_type(self.mol52, 5) == "Csc"
        assert self.atom_type(self.mol1, 5) == "Cd"
        assert self.atom_type(self.mol60, 1) == "Cdc"
        assert self.atom_type(self.mol1, 2) == "CO"
        assert self.atom_type(self.mol19, 0) == "CS"
        assert self.atom_type(self.mol1, 6) == "Cdd"
        assert self.atom_type(self.mol1, 9) == "Ct"
        assert self.atom_type(self.mol2, 0) == "Cb"
        assert self.atom_type(self.mol55, 3) == "Cbf"
        assert self.atom_type(self.mol56, 0) == "C2s"
        assert self.atom_type(self.mol57, 0) == "C2sc"
        assert self.atom_type(self.mol58, 0) == "C2d"
        assert self.atom_type(self.mol59, 0) == "C2dc"
        assert self.atom_type(self.mol60, 2) == "C2dc"
        assert self.atom_type(self.mol20, 0) == "C2tc"
        assert self.atom_type(self.mol29, 0) == "C2tc"  # todo: add in a ciq unit test?
        assert self.atom_type(self.mol95, 0) == "Ctc"

    def test_nitrogen_types(self):
        """
        Test that get_atomtype() returns appropriate nitrogen atom types.
        """
        assert self.atom_type(self.mol40, 3) == "N0sc"
        assert self.atom_type(self.mol41, 0) == "N1s"
        assert self.atom_type(self.mol39, 0) == "N1sc"
        assert self.atom_type(self.mol5, 3) == "N1dc"
        assert self.atom_type(self.mol9, 0) == "N3s"
        assert self.atom_type(self.mol10, 0) == "N3s"
        assert self.atom_type(self.mol11, 0) == "N3s"
        assert self.atom_type(self.mol16, 0) == "N3d"
        assert self.atom_type(self.mol17, 0) == "N3d"
        assert self.atom_type(self.mol12, 0) == "N3t"
        assert self.atom_type(self.mol18, 5) == "N3b"
        assert self.atom_type(self.mol5, 2) == "N5dc"
        assert self.atom_type(self.mol64, 1) == "N5ddc"
        assert self.atom_type(self.mol53, 0) == "N5dddc"
        assert self.atom_type(self.mol15, 1) == "N5tc"
        assert self.atom_type(self.mol39, 2) == "N5tc"
        assert self.atom_type(self.mol18, 0) == "N5b"

    def test_oxygen_types(self):
        """
        Test that get_atomtype() returns appropriate oxygen atom types.
        """
        assert self.atom_type(self.mol44, 0) == "Oa"
        assert self.atom_type(self.mol45, 2) == "O0sc"
        assert self.atom_type(self.mol49, 0) == "O0sc"
        assert self.atom_type(self.mol1, 1) == "O2s"
        assert self.atom_type(self.mol24, 2) == "O2sc"
        assert self.atom_type(self.mol1, 3) == "O2d"
        assert self.atom_type(self.mol49, 1) == "O4sc"
        assert self.atom_type(self.mol50, 1) == "O4dc"
        assert self.atom_type(self.mol20, 1) == "O4tc"

    def test_silicon_types(self):
        """
        Test that get_atomtype() returns appropriate silicon atom types.
        """
        assert self.atom_type(self.mol4, 2) == "Sis"
        assert self.atom_type(self.mol4, 1) == "SiO"
        assert self.atom_type(self.mol4, 5) == "Sid"
        assert self.atom_type(self.mol4, 4) == "Sidd"
        assert self.atom_type(self.mol4, 7) == "Sit"  # todo: add in Siq unit test?

    def test_phosphorus_types(self):
        """
        Test that get_atomtype() returns appropriate phosphorus atom types.
        """
        assert self.atom_type(self.mol80, 0) == "P0sc"
        assert self.atom_type(self.mol81, 0) == "P1s"
        assert self.atom_type(self.mol82, 0) == "P1sc"
        assert self.atom_type(self.mol83, 3) == "P1dc"
        assert self.atom_type(self.mol84, 0) == "P3s"
        assert self.atom_type(self.mol84, 1) == "P3d"
        assert self.atom_type(self.mol85, 0) == "P3t"
        assert self.atom_type(self.mol86, 5) == "P3b"
        assert self.atom_type(self.mol87, 0) == "P5s"
        assert self.atom_type(self.mol88, 0) == "P5sc"
        assert self.atom_type(self.mol89, 0) == "P5d"
        assert self.atom_type(self.mol90, 0) == "P5dd"
        assert self.atom_type(self.mol83, 2) == "P5dc"
        assert self.atom_type(self.mol91, 0) == "P5ddc"
        assert self.atom_type(self.mol92, 0) == "P5t"
        assert self.atom_type(self.mol93, 0) == "P5td"
        assert self.atom_type(self.mol80, 1) == "P5tc"
        assert self.atom_type(self.mol86, 0) == "P5b"
        assert self.atom_type(self.mol94, 0) == "P5bd"

    def test_sulfur_types(self):
        """
        Test that get_atomtype() returns appropriate sulfur atom types.
        """
        assert self.atom_type(self.mol22, 0) == "Sa"
        assert self.atom_type(self.mol21, 0) == "S0sc"
        assert self.atom_type(self.mol23, 0) == "S2s"
        assert self.atom_type(self.mol21, 1) == "S2sc"
        assert self.atom_type(self.mol42, 2) == "S2sc"
        assert self.atom_type(self.mol19, 1) == "S2d"
        assert self.atom_type(self.mol24, 1) == "S2dc"
        assert self.atom_type(self.mol69, 1) == "S2tc"
        assert self.atom_type(self.mol25, 0) == "S4s"
        assert self.atom_type(self.mol23, 1) == "S4sc"
        assert self.atom_type(self.mol25, 2) == "S4d"
        assert self.atom_type(self.mol28, 1) == "S4dd"
        assert self.atom_type(self.mol38, 1) == "S4dd"
        assert self.atom_type(self.mol26, 1) == "S4dc"
        assert self.atom_type(self.mol28, 4) == "S4t"
        assert self.atom_type(self.mol29, 1) == "S4tdc"
        assert self.atom_type(self.mol28, 5) == "S6s"
        assert self.atom_type(self.mol51, 1) == "S6sc"
        assert self.atom_type(self.mol30, 0) == "S6d"
        assert self.atom_type(self.mol32, 1) == "S6dd"
        assert self.atom_type(self.mol34, 1) == "S6ddd"
        assert self.atom_type(self.mol43, 1) == "S6dc"
        assert self.atom_type(self.mol31, 0) == "S6dc"
        assert self.atom_type(self.mol33, 1) == "S6dc"
        assert self.atom_type(self.mol35, 0) == "S6t"
        assert self.atom_type(self.mol36, 0) == "S6td"
        assert self.atom_type(self.mol37, 1) == "S6tt"
        assert self.atom_type(self.mol70, 0) == "S6tdc"

    def test_chlorine_types(self):
        """
        Test that get_atomtype() returns appropriate chlorine atom types.
        """
        assert self.atom_type(self.mol73, 1) == "Cl1s"
        assert self.atom_type(self.mol_chloride, 0) == "Cl0sc"

    def test_bromine_types(self):
        """
        Test that get_atomtype() returns appropriate bromine atom types.
        """
        assert self.atom_type(self.mol79, 1) == "Br1s"
        assert self.atom_type(self.mol_bromide, 0) == "Br0sc"

    def test_iodine_types(self):
        """
        Test that get_atomtype() returns appropriate iodine atom types.
        """
        assert self.atom_type(self.mol74, 1) == "I1s"

    def test_fluorine_types(self):
        """
        Test that get_atomtype() returns appropriate fluorine atom types.
        """
        assert self.atom_type(self.mol75, 1) == "F1s"
        assert self.atom_type(self.mol75_anion, 0) == "F0sc"
        assert self.atom_type(self.mol75_cation, 0) == "F1sc"

    def test_sf6_anion_types(self):
        """
        Test that get_atomtype() resolves every atom of the SF6 radical anion.
        """
        assert self.atom_type(self.mol_sf6_anion, 0) == "S6sc"
        for index in range(1, 7):
            assert self.atom_type(self.mol_sf6_anion, index) == "F1s"

    def test_lithium_types(self):
        """
        Test that get_atomtype() returns appropriate lithium atom types.
        """
        assert self.atom_type(self.mol96, 0) == "Li+"
        assert self.atom_type(self.mol97, 0) == "Li0"

    def test_other_types(self):
        """
        Test that get_atomtype() returns appropriate types for other misc inerts.
        """
        assert self.atom_type(self.mol6, 0) == "Ar"
        assert self.atom_type(self.mol7, 0) == "He"
        assert self.atom_type(self.mol8, 0) == "Ne"

    def test_occupied_surface_atom_type(self):
        """
        Test that get_atomtype() works for occupied surface sites and for regular atoms in the complex.
        """
        assert self.atom_type(self.mol76, 0) == "H0"
        assert self.atom_type(self.mol76, 1) == "Xo"

    def test_vacant_surface_site_atom_type(self):
        """
        Test that get_atomtype() works for vacant surface sites and for regular atoms in the complex.
        """
        assert self.atom_type(self.mol77, 0) == "Cs"
        assert self.atom_type(self.mol77, 1) == "H0"
        assert self.atom_type(self.mol77, 3) == "Xv"
        assert self.atom_type(self.mol78, 0) == "Xv"

    def test_electron(self):
        """
        Test that get_atomtype() returns the electron (e) atom type.
        """
        assert self.atom_type(self.electron, 0) == 'e'

    def test_proton(self):
        """
        Test that get_atomtype() returns the proton (H+) atom type.
        """
        assert self.atom_type(self.proton, 0) == 'H+'


class TestAnionRecipes:
    """
    Contains tests that a reaction recipe reaching a monatomic anion actually produces the anion.

    These deliberately assert on the net charge and the formula of the product rather than on an
    exception being raised. The defect they pin was silent: applying the recipe returned a species
    of the opposite charge with no exception and no log line, so a ``pytest.raises`` test would
    have passed against the bug and proved nothing.
    """

    def attach_electron(self, adjlist):
        """Apply the electron-attachment recipe (radical -> lone pair, charge -1) to `adjlist`."""
        from rmgpy.data.kinetics.family import ReactionRecipe

        mol = Molecule().from_adjacency_list(adjlist)
        product = mol.copy(deep=True)
        recipe = ReactionRecipe([["LOSE_RADICAL", "*1", 1], ["GAIN_PAIR", "*1", 1]])
        recipe.apply_forward(product, unique=True)
        product.update()
        return product

    def test_hydride_from_recipe(self):
        """A recipe reaching H- must yield a product of net charge -1, not the H+ it used to."""
        product = self.attach_electron("""1 *1 H u1 p0 c0""")
        assert product.get_net_charge() == -1
        assert product.get_formula() == "H"
        assert product.atoms[0].atomtype.label == "H-"

    def test_chloride_from_recipe(self):
        """A recipe reaching Cl- must yield a product of net charge -1."""
        product = self.attach_electron("""1 *1 Cl u1 p3 c0""")
        assert product.get_net_charge() == -1
        assert product.get_formula() == "Cl"
        assert product.atoms[0].atomtype.label == "Cl0sc"

    def test_bromide_from_recipe(self):
        """A recipe reaching Br- must yield a product of net charge -1."""
        product = self.attach_electron("""1 *1 Br u1 p3 c0""")
        assert product.get_net_charge() == -1
        assert product.get_formula() == "Br"
        assert product.atoms[0].atomtype.label == "Br0sc"


class TestGenericLonePairActions:
    """
    Contains tests that GAIN_PAIR is a legal action on the *generic* atomtype of every element
    whose anion this branch supports.

    A reaction family builds its product template by applying the recipe to the root group
    itself, so a recipe containing GAIN_PAIR can only ever match an atom whose generic atomtype
    permits the action. Adding the H-/F0sc/Cl0sc/Br0sc leaves does not by itself make attachment
    reachable: with an empty increment_lone_pair on generic H/F/Cl/Br, the root group cannot
    include those elements at all, so the leaf exists and the capability does not.
    """

    def test_generic_lone_pair_actions_are_wired(self):
        """The generic atomtypes of the supported anions must permit lone pair actions."""
        from rmgpy.molecule.atomtype import ATOMTYPES

        for label in ["H", "F", "Cl", "Br"]:
            assert ATOMTYPES[label].increment_lone_pair, f"{label} has no increment_lone_pair action"
            assert ATOMTYPES[label].decrement_lone_pair, f"{label} has no decrement_lone_pair action"

    def test_gain_pair_on_generic_root_group(self):
        """Applying the attachment recipe to a root group of each generic atomtype must not raise."""
        from rmgpy.data.kinetics.family import ReactionRecipe
        from rmgpy.molecule.group import Group

        for label in ["H", "F", "Cl", "Br"]:
            group = Group().from_adjacency_list(f"1 *1 {label} u1")
            recipe = ReactionRecipe([["LOSE_RADICAL", "*1", 1], ["GAIN_PAIR", "*1", 1]])
            recipe.apply_forward(group, unique=True)  # raises ActionError if the action is missing
            assert group.atoms[0].radical_electrons == [0]

    def test_gain_pair_on_halogen_radical_yields_anion(self):
        """The recipe applied to a real halogen radical must yield the anion, not the cation."""
        from rmgpy.data.kinetics.family import ReactionRecipe

        for adjlist, expected in (
            ("1 *1 Cl u1 p3 c0", "Cl0sc"),
            ("1 *1 Br u1 p3 c0", "Br0sc"),
            ("1 *1 F u1 p3 c0", "F0sc"),
            ("1 *1 H u1 p0 c0", "H-"),
        ):
            product = Molecule().from_adjacency_list(adjlist)
            ReactionRecipe([["LOSE_RADICAL", "*1", 1], ["GAIN_PAIR", "*1", 1]]).apply_forward(product, unique=True)
            product.update()
            assert product.get_net_charge() == -1
            assert product.atoms[0].atomtype.label == expected


class TestChargeClassParents:
    """
    Contains tests that the charge-class parent atomtypes stay group-matching labels.

    They are reachable from a group definition through the specific/generic lists, but must never
    win a get_atomtype lookup against their own leaves, because their feature ranges are only
    consulted there. A parent declared wider than its leaves resolves atoms the leaves reject.
    """

    def test_parents_are_reachable_as_group_labels(self):
        """Each parent must be linked to its leaves in both directions."""
        from rmgpy.molecule.atomtype import ATOMTYPES

        for parent, leaves in (("Om1", ["O0sc"]), ("Nm1", ["N1sc", "N1dc"]), ("Nm2", ["N0sc"])):
            for leaf in leaves:
                assert ATOMTYPES[leaf].is_specific_case_of(ATOMTYPES[parent])
                assert ATOMTYPES[leaf] in ATOMTYPES[parent].specific
                assert ATOMTYPES[parent] in ATOMTYPES[leaf].generic

    def test_parents_do_not_shadow_their_leaves(self):
        """get_atomtype must still return the leaf, not the charge-class parent."""
        for adjlist, expected in (
            # hydroxide [OH-]
            ("1 O u0 p3 c-1 {2,S}\n2 H u0 p0 c0 {1,S}", "O0sc"),
            # amide [NH2-]
            ("1 N u0 p2 c-1 {2,S} {3,S}\n2 H u0 p0 c0 {1,S}\n3 H u0 p0 c0 {1,S}", "N1sc"),
        ):
            mol = Molecule().from_adjacency_list(adjlist)
            assert mol.atoms[0].atomtype.label == expected

    def test_parents_admit_no_octet_violation(self):
        """
        No charge-class parent may resolve an atom carrying more than eight electrons.

        Om1 as originally written covered `O u1 p2 c-1` (nine electrons), which added a spurious
        third resonance structure to NO2 and a fourth and fifth to SO2.
        """
        from rmgpy.molecule.atomtype import ATOMTYPES

        valence = {"O": 6, "N": 5}
        for parent, symbol in (("Om1", "O"), ("Om2", "O"), ("Nm1", "N"), ("Nm2", "N"), ("Nm3", "N")):
            at = ATOMTYPES[parent]
            for single in at.single or [0]:
                for double in at.all_double or [0]:
                    for triple in at.triple or [0]:
                        for pairs in at.lone_pairs or [0]:
                            for charge in at.charge or [0]:
                                order = single + 2 * double + 3 * triple
                                radicals = valence[symbol] - order - 2 * pairs - charge
                                if radicals < 0:
                                    continue  # not a real atom
                                electrons = 2 * pairs + radicals + 2 * order
                                assert electrons <= 8, (
                                    f"{parent} resolves {symbol} u{radicals} p{pairs} c{charge:+d} "
                                    f"with {order} bond orders, carrying {electrons} electrons"
                                )
