################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

import unittest
import molecule

from .resonance import *
from .filtration import *
################################################################################


class FiltrationTest(unittest.TestCase):

    def basic_filtration_test(self):
        """Test that structures with higher octet deviation get filtered out"""
        adj1 = """
        multiplicity 2
        1 N u0 p1 c0 {2,D} {3,S}
        2 O u0 p2 c0 {1,D}
        3 O u1 p2 c0 {1,S}
        """
        adj2 = """
        multiplicity 2
        1 N u1 p1 c0 {2,S} {3,S}
        2 O u0 p2 c+1 {1,S}
        3 O u0 p3 c-1 {1,S}
        """
        adj3 = """
        multiplicity 2
        1 N u1 p0 c0 {2,T} {3,S}
        2 O u0 p1 c+1 {1,T}
        3 O u0 p3 c-1 {1,S}
        """

        mol1 = molecule.Molecule().fromAdjacencyList(adj1)
        mol2 = molecule.Molecule().fromAdjacencyList(adj2)
        mol3 = molecule.Molecule().fromAdjacencyList(adj3)

        mol_list = [mol1,mol2,mol3]
        octet_deviation_list = get_octet_deviation_list(mol_list)
        filtered_list = filter_structures(mol_list)

        self.assertEqual(octet_deviation_list,[1,3,3])
        self.assertEqual(len(filtered_list), 1)
        self.assertTrue(all([atom.charge == 0 for atom in filtered_list[0].vertices]))

    def penalty_for_N_val_9_test(self):
        """Test that N atoms with valance 9 get penalized in the octet deviation score"""
        adj = """
        multiplicity 2
        1 N u1 p0 c0 {2,S} {3,T}
        2 O u0 p2 c0 {1,S} {4,S}
        3 N u0 p1 c0 {1,T}
        4 H u0 p0 c0 {2,S}
        """
        mol = molecule.Molecule().fromAdjacencyList(adj)
        octet_deviation = get_octet_deviation(mol)
        self.assertEqual(octet_deviation, 2)

    def penalty_for_O4tc_test(self):
        """Test that an O4tc atomType with octet 8 gets penalized in the octet deviation score"""
        adj = """
        1 S u0 p1 c0 {2,S} {3,T}
        2 O u0 p3 c-1 {1,S}
        3 O u0 p1 c+1 {1,T}
        """
        mol = molecule.Molecule().fromAdjacencyList(adj)
        octet_deviation = get_octet_deviation(mol)
        self.assertEqual(octet_deviation, 1)
        self.assertEqual(mol.vertices[2].atomType.label, 'O4tc')

    def penalty_for_S_triple_S_test(self):
        """Test that an S#S substructure in a molecule gets penalized in the octet deviation score"""
        adj = """
        1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
        2  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
        3  S u0 p0 c0 {1,S} {4,T} {11,D}
        4  S u0 p1 c0 {2,S} {3,T}
        5  H u0 p0 c0 {1,S}
        6  H u0 p0 c0 {1,S}
        7  H u0 p0 c0 {1,S}
        8  H u0 p0 c0 {2,S}
        9  H u0 p0 c0 {2,S}
        10 H u0 p0 c0 {2,S}
        11 O u0 p2 c0 {3,D}
        """
        mol = molecule.Molecule().fromAdjacencyList(adj)
        octet_deviation = get_octet_deviation(mol)
        self.assertEqual(octet_deviation, 1)

    def penalty_for_fake_SO_test(self):
        """Test that an SO structure which has a multiplicity 3 and two radicals as the correct [S][O] structure, yet has
        both radicals on the same site gets penalized in the octet deviation score"""
        adj = """
        multiplicity 3
        1 O u0 p2 c0 {2,D}
        2 S u2 p1 c0 {1,D}
        """
        mol = molecule.Molecule().fromAdjacencyList(adj)
        octet_deviation = get_octet_deviation(mol)
        self.assertEqual(octet_deviation, 3)

    def is_OS_test(self):
        """Test that the is_OS() function works as expected"""
        smiles = "[O][S]=O"
        mol = Molecule().fromSMILES(smiles)
        self.assertEqual(is_OS(mol), 0)  # 0 - neither O2, S2, or SO

        smiles = "[S][S]"
        mol = Molecule().fromSMILES(smiles)
        self.assertEqual(is_OS(mol), 1)  # 1 - triplet ground state ([O.][O.], [S.][S.], or [S.][O.])

        smiles = "O=O"
        mol = Molecule().fromSMILES(smiles)
        self.assertEqual(is_OS(mol), 2)  # 2 - singlet excited state (O=O, S=S, or S=O)

        adj = """
        multiplicity 3
        1 O u0 p2 c0 {2,D}
        2 S u2 p1 c0 {1,D}
        """
        mol = molecule.Molecule().fromAdjacencyList(adj)
        self.assertEqual(is_OS(mol), 3)  # 3 - a O2/S2/SO structure which is neither case `1` or `2`
