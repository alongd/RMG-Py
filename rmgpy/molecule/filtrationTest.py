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

        filtered_list = filter_structures(mol_list)

        self.assertEqual(len(filtered_list), 1)
        self.assertTrue(all([atom.charge == 0 for atom in filtered_list[0].vertices]))




# tests specific species that get penalized, test OS and adjBirads, test the reactive attribute
























































