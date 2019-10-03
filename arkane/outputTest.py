#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
This module contains unit tests of the :mod:`arkane.output` module.
"""

import unittest

from rmgpy.statmech.conformer import Conformer
from rmgpy.statmech.translation import IdealGasTranslation
from rmgpy.statmech.rotation import NonlinearRotor
from rmgpy.statmech.vibration import HarmonicOscillator

from arkane.output import prettify

################################################################################


class OutputTest(unittest.TestCase):
    """
    Contains unit tests of Arkane's output module.
    """

    @classmethod
    def setUpClass(cls):
        """
        A method that is run before all unit tests in this class.
        """
        cls.maxDiff = None

    def test_prettify(self):
        """Test the prettify function"""
        # test a simple dictionary
        d = {'A': 5, 'B': 6, 'C': 20}
        pretty_d = prettify(d)
        self.assertEqual(pretty_d, "{'A': 5, 'B': 6, 'C': 20}")

        # test a complex object such as Conformer
        e0 = (237.219, 'kJ/mol')
        modes = [IdealGasTranslation(mass=(43.0296, 'amu')),
                 NonlinearRotor(inertia = ([4.19833, 48.1107, 52.309], 'amu*angstrom^2'), symmetry = 1),
                 HarmonicOscillator(frequencies = ([77.4703, 307.323, 478.716, 666.453, 958.485, 1090.88, 1254.11,
                                                    1584.84, 1792.08, 3529.53, 3559.75, 3710.93], 'cm^-1'))]
        conf = Conformer(E0=e0, modes=modes, spin_multiplicity=2)
        pretty_conf = prettify(conf.as_dict())
        print(pretty_conf)
        expected_conf = """{   'E0': {'class': 'ScalarQuantity', 'units': 'kJ/mol', 'value': 237.219},
    'class': 'Conformer',
    'modes': [   {   'class': 'IdealGasTranslation',
                     'mass': {   'class': 'ScalarQuantity',
                                 'units': 'amu',
                                 'value': 43.0296},
                     'quantum': False},
                 {   'class': 'NonlinearRotor',
                     'inertia': {   'class': 'ArrayQuantity',
                                    'units': 'amu*angstrom^2',
                                    'value': {   'class': 'np_array',
                                                 'object': [   4.19833,
                                                               48.1107,
                                                               52.309]}},
                     'quantum': False,
                     'rotationalConstant': {   'class': 'ArrayQuantity',
                                               'units': 'cm^-1',
                                               'value': {   'class': 'np_array',
                                                            'object': [   4.01531777000581,
                                                                          0.3503925125460343,
                                                                          0.3222701457368425]}},
                     'symmetry': 1},
                 {   'class': 'HarmonicOscillator',
                     'frequencies': {   'class': 'ArrayQuantity',
                                        'units': 'cm^-1',
                                        'value': {   'class': 'np_array',
                                                     'object': [   77.4703,
                                                                   307.323,
                                                                   478.716,
                                                                   666.453,
                                                                   958.485,
                                                                   1090.88,
                                                                   1254.11,
                                                                   1584.84,
                                                                   1792.08,
                                                                   3529.53,
                                                                   3559.75,
                                                                   3710.93]}},
                     'quantum': True}],
    'optical_isomers': 1,
    'spin_multiplicity': 2}"""
        self.assertEqual(pretty_conf, expected_conf)


################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
