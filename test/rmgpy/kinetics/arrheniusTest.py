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
This script contains unit tests of the :mod:`rmgpy.kinetics.arrhenius` module.
"""

import copy
import inspect
import math
import os
import pickle


import numpy as np

import rmgpy
import rmgpy.constants as constants
from rmgpy import settings
from rmgpy.kinetics.arrhenius import (
    Arrhenius,
    ArrheniusEP,
    ArrheniusBM,
    PDepArrhenius,
    MultiArrhenius,
    MultiPDepArrhenius,
    TwoTemperaturePlasma,
    ElectronCollisionPlasma,
    BadnellRRArrhenius,
    VoronovEIArrhenius,
)
from rmgpy.molecule.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial, ThermoData
import pytest


class TestArrhenius:
    """
    Contains unit tests of the :class:`Arrhenius` class.
    """

    def setup_method(self):
        self.A = 1.0e12
        self.n = 0.5
        self.Ea = 41.84
        self.T0 = 1.0
        self.Tmin = 300.0
        self.Tmax = 3000.0
        self.comment = "C2H6"
        self.arrhenius = Arrhenius(
            A=(self.A, "cm^3/(mol*s)"),
            n=self.n,
            Ea=(self.Ea, "kJ/mol"),
            T0=(self.T0, "K"),
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment=self.comment,
        )

    def test_a_factor(self):
        """
        Test that the Arrhenius A property was properly set.
        """
        assert abs(self.arrhenius.A.value_si * 1e6 - self.A) < 1e0

    def test_n(self):
        """
        Test that the Arrhenius n property was properly set.
        """
        assert round(abs(self.arrhenius.n.value_si - self.n), 6) == 0

    def test_ea(self):
        """
        Test that the Arrhenius Ea property was properly set.
        """
        assert round(abs(self.arrhenius.Ea.value_si * 0.001 - self.Ea), 6) == 0

    def test_temperature0(self):
        """
        Test that the Arrhenius T0 property was properly set.
        """
        assert round(abs(self.arrhenius.T0.value_si - self.T0), 6) == 0

    def test_temperature_min(self):
        """
        Test that the Arrhenius Tmin property was properly set.
        """
        assert round(abs(self.arrhenius.Tmin.value_si - self.Tmin), 6) == 0

    def test_temperature_max(self):
        """
        Test that the Arrhenius Tmax property was properly set.
        """
        assert round(abs(self.arrhenius.Tmax.value_si - self.Tmax), 6) == 0

    def test_comment(self):
        """
        Test that the Arrhenius comment property was properly set.
        """
        assert self.arrhenius.comment == self.comment

    def test_is_temperature_valid(self):
        """
        Test the Arrhenius.is_temperature_valid() method.
        """
        Tdata = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        validdata = np.array([False, True, True, True, True, True, True, True, True, True], bool)
        for T, valid in zip(Tdata, validdata):
            valid0 = self.arrhenius.is_temperature_valid(T)
            assert valid0 == valid

    def test_get_rate_coefficient(self):
        """
        Test the Arrhenius.get_rate_coefficient() method.
        """
        Tlist = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        kexplist = np.array(
            [
                1.6721e-4,
                6.8770e1,
                5.5803e3,
                5.2448e4,
                2.0632e5,
                5.2285e5,
                1.0281e6,
                1.7225e6,
                2.5912e6,
                3.6123e6,
            ]
        )
        for T, kexp in zip(Tlist, kexplist):
            kact = self.arrhenius.get_rate_coefficient(T)
            assert abs(kexp - kact) < 1e-4 * kexp

    def test_change_t0(self):
        """
        Test the Arrhenius.change_t0() method.
        """
        Tlist = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500])
        k0list = np.array([self.arrhenius.get_rate_coefficient(T) for T in Tlist])
        self.arrhenius.change_t0(300)
        assert self.arrhenius.T0.value_si == 300
        for T, kexp in zip(Tlist, k0list):
            kact = self.arrhenius.get_rate_coefficient(T)
            assert abs(kexp - kact) < 1e-6 * kexp

    def test_fit_to_data(self):
        """
        Test the Arrhenius.fit_to_data() method.
        """
        Tdata = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500])
        kdata = np.array([self.arrhenius.get_rate_coefficient(T) for T in Tdata])
        arrhenius = Arrhenius().fit_to_data(Tdata, kdata, kunits="m^3/(mol*s)")
        assert float(self.arrhenius.T0.value_si) == 1
        for T, k in zip(Tdata, kdata):
            assert abs(k - arrhenius.get_rate_coefficient(T)) < 1e-6 * k
        assert abs(arrhenius.A.value_si - self.arrhenius.A.value_si) < 1e0
        assert round(abs(arrhenius.n.value_si - self.arrhenius.n.value_si), 1) == 0, 4
        assert round(abs(arrhenius.Ea.value_si - self.arrhenius.Ea.value_si), 2) == 0
        assert round(abs(arrhenius.T0.value_si - self.arrhenius.T0.value_si), 4) == 0

    def test_fit_to_negative_data(self):
        """
        Test the Arrhenius.fit_to_data() method on negative rates
        """
        Tdata = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500])
        kdata = np.array([-1 * self.arrhenius.get_rate_coefficient(T) for T in Tdata])
        arrhenius = Arrhenius().fit_to_data(Tdata, kdata, kunits="m^3/(mol*s)")
        assert float(self.arrhenius.T0.value_si) == 1
        for T, k in zip(Tdata, kdata):
            assert abs(k - arrhenius.get_rate_coefficient(T)) < 1e-6 * abs(k)
        assert abs(arrhenius.A.value_si - -1 * self.arrhenius.A.value_si) < 1e0
        assert round(abs(arrhenius.n.value_si - self.arrhenius.n.value_si), 1) == 0, 4
        assert round(abs(arrhenius.Ea.value_si - self.arrhenius.Ea.value_si), 2) == 0
        assert round(abs(arrhenius.T0.value_si - self.arrhenius.T0.value_si), 4) == 0

    def test_pickle(self):
        """
        Test that an Arrhenius object can be pickled and unpickled with no loss
        of information.
        """
        import pickle

        arrhenius = pickle.loads(pickle.dumps(self.arrhenius, -1))
        assert abs(self.arrhenius.A.value - arrhenius.A.value) < 1e0
        assert self.arrhenius.A.units == arrhenius.A.units
        assert round(abs(self.arrhenius.n.value - arrhenius.n.value), 4) == 0
        assert round(abs(self.arrhenius.Ea.value - arrhenius.Ea.value), 4) == 0
        assert self.arrhenius.Ea.units == arrhenius.Ea.units
        assert round(abs(self.arrhenius.T0.value - arrhenius.T0.value), 4) == 0
        assert self.arrhenius.T0.units == arrhenius.T0.units
        assert round(abs(self.arrhenius.Tmin.value - arrhenius.Tmin.value), 4) == 0
        assert self.arrhenius.Tmin.units == arrhenius.Tmin.units
        assert round(abs(self.arrhenius.Tmax.value - arrhenius.Tmax.value), 4) == 0
        assert self.arrhenius.Tmax.units == arrhenius.Tmax.units
        assert self.arrhenius.comment == arrhenius.comment

    def test_repr(self):
        """
        Test that an Arrhenius object can be reconstructed from its repr()
        output with no loss of information.
        """
        namespace = {}
        exec("arrhenius = {0!r}".format(self.arrhenius), globals(), namespace)
        assert "arrhenius" in namespace
        arrhenius = namespace["arrhenius"]
        assert abs(self.arrhenius.A.value - arrhenius.A.value) < 1e0
        assert self.arrhenius.A.units == arrhenius.A.units
        assert round(abs(self.arrhenius.n.value - arrhenius.n.value), 4) == 0
        assert round(abs(self.arrhenius.Ea.value - arrhenius.Ea.value), 4) == 0
        assert self.arrhenius.Ea.units == arrhenius.Ea.units
        assert round(abs(self.arrhenius.T0.value - arrhenius.T0.value), 4) == 0
        assert self.arrhenius.T0.units == arrhenius.T0.units
        assert round(abs(self.arrhenius.Tmin.value - arrhenius.Tmin.value), 4) == 0
        assert self.arrhenius.Tmin.units == arrhenius.Tmin.units
        assert round(abs(self.arrhenius.Tmax.value - arrhenius.Tmax.value), 4) == 0
        assert self.arrhenius.Tmax.units == arrhenius.Tmax.units
        assert self.arrhenius.comment == arrhenius.comment

    def test_change_rate(self):
        """
        Test the Arrhenius.change_rate() method.
        """
        Tlist = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500])
        k0list = np.array([self.arrhenius.get_rate_coefficient(T) for T in Tlist])
        self.arrhenius.change_rate(2)
        for T, kexp in zip(Tlist, k0list):
            kact = self.arrhenius.get_rate_coefficient(T)
            assert abs(2 * kexp - kact) < 1e-6 * kexp

    def test_to_cantera_kinetics(self):
        """
        Test that the Arrhenius cantera object can be set properly within
        a cantera Reaction object
        """
        ctArrhenius = self.arrhenius.to_cantera_kinetics()
        assert round(abs(ctArrhenius.pre_exponential_factor - 1e9), 6) == 0
        assert round(abs(ctArrhenius.temperature_exponent - 0.5), 7) == 0
        assert round(abs(ctArrhenius.activation_energy - 41.84e6), 7) == 0

    def test_to_arrhenius_ep(self):
        """
        Tests that the Arrhenius object can be converted to ArrheniusEP
        """
        arr_rate = self.arrhenius.get_rate_coefficient(500)
        arr_ep = self.arrhenius.to_arrhenius_ep()
        arr_ep_rate = arr_ep.get_rate_coefficient(500, 10)  # the second number should not matter
        assert round(abs(arr_rate - arr_ep_rate), 7) == 0

    def test_to_arrhenius_ep_with_alpha_and_hrxn(self):
        """
        Tests that the Arrhenius object can be converted to ArrheniusEP given parameters
        """
        hrxn = 5
        arr_rate = self.arrhenius.get_rate_coefficient(500)
        arr_ep = self.arrhenius.to_arrhenius_ep(alpha=1, dHrxn=hrxn)
        assert round(abs(1.0 - arr_ep.alpha.value_si), 7) == 0
        arr_ep_rate = arr_ep.get_rate_coefficient(500, hrxn)
        assert round(abs(arr_rate - arr_ep_rate), 7) == 0

    def test_to_arrhenius_ep_throws_error_with_just_alpha(self):
        with pytest.raises(Exception):
            self.arrhenius.to_arrhenius_ep(alpha=1)


class TestArrheniusEP:
    """
    Contains unit tests of the :class:`ArrheniusEP` class.
    """

    def setup_method(self):
        self.A = 1.0e12
        self.n = 0.5
        self.alpha = 0.5
        self.E0 = 41.84
        self.Tmin = 300.0
        self.Tmax = 3000.0
        self.comment = "C2H6"
        self.arrhenius = ArrheniusEP(
            A=(self.A, "cm^3/(mol*s)"),
            n=self.n,
            alpha=self.alpha,
            E0=(self.E0, "kJ/mol"),
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment=self.comment,
        )

    def test_a_factor(self):
        """
        Test that the ArrheniusEP A property was properly set.
        """
        assert abs(self.arrhenius.A.value_si * 1e6 - self.A) < 1e0

    def test_n(self):
        """
        Test that the ArrheniusEP n property was properly set.
        """
        assert round(abs(self.arrhenius.n.value_si - self.n), 6) == 0

    def test_alpha(self):
        """
        Test that the ArrheniusEP alpha property was properly set.
        """
        assert round(abs(self.arrhenius.alpha.value_si - self.alpha), 6) == 0

    def test_e0(self):
        """
        Test that the ArrheniusEP E0 property was properly set.
        """
        assert round(abs(self.arrhenius.E0.value_si * 0.001 - self.E0), 6) == 0

    def test_temperature_min(self):
        """
        Test that the ArrheniusEP Tmin property was properly set.
        """
        assert round(abs(self.arrhenius.Tmin.value_si - self.Tmin), 6) == 0

    def test_temperature_max(self):
        """
        Test that the ArrheniusEP Tmax property was properly set.
        """
        assert round(abs(self.arrhenius.Tmax.value_si - self.Tmax), 6) == 0

    def test_comment(self):
        """
        Test that the ArrheniusEP comment property was properly set.
        """
        assert self.arrhenius.comment == self.comment

    def test_is_temperature_valid(self):
        """
        Test the ArrheniusEP.is_temperature_valid() method.
        """
        Tdata = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        validdata = np.array([False, True, True, True, True, True, True, True, True, True], bool)
        for T, valid in zip(Tdata, validdata):
            valid0 = self.arrhenius.is_temperature_valid(T)
            assert valid0 == valid

    def test_get_rate_coefficient(self):
        """
        Test the ArrheniusEP.get_rate_coefficient() method.
        """
        Tlist = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        kexplist = np.array(
            [
                1.6721e-4,
                6.8770e1,
                5.5803e3,
                5.2448e4,
                2.0632e5,
                5.2285e5,
                1.0281e6,
                1.7225e6,
                2.5912e6,
                3.6123e6,
            ]
        )
        for T, kexp in zip(Tlist, kexplist):
            kact = self.arrhenius.get_rate_coefficient(
                T,
            )
            assert abs(kexp - kact) < 1e-4 * kexp

    def test_pickle(self):
        """
        Test that an ArrheniusEP object can be pickled and unpickled with no loss
        of information.
        """
        import pickle

        arrhenius = pickle.loads(pickle.dumps(self.arrhenius, -1))
        assert abs(self.arrhenius.A.value - arrhenius.A.value) < 1e0
        assert self.arrhenius.A.units == arrhenius.A.units
        assert round(abs(self.arrhenius.n.value - arrhenius.n.value), 4) == 0
        assert round(abs(self.arrhenius.alpha.value - arrhenius.alpha.value), 4) == 0
        assert round(abs(self.arrhenius.E0.value - arrhenius.E0.value), 4) == 0
        assert self.arrhenius.E0.units == arrhenius.E0.units
        assert round(abs(self.arrhenius.Tmin.value - arrhenius.Tmin.value), 4) == 0
        assert self.arrhenius.Tmin.units == arrhenius.Tmin.units
        assert round(abs(self.arrhenius.Tmax.value - arrhenius.Tmax.value), 4) == 0
        assert self.arrhenius.Tmax.units == arrhenius.Tmax.units
        assert self.arrhenius.comment == arrhenius.comment

    def test_repr(self):
        """
        Test that an ArrheniusEP object can be reconstructed from its repr()
        output with no loss of information.
        """
        namespace = {}
        exec("arrhenius = {0!r}".format(self.arrhenius), globals(), namespace)
        assert "arrhenius" in namespace
        arrhenius = namespace["arrhenius"]
        assert abs(self.arrhenius.A.value - arrhenius.A.value) < 1e0
        assert self.arrhenius.A.units == arrhenius.A.units
        assert round(abs(self.arrhenius.n.value - arrhenius.n.value), 4) == 0
        assert round(abs(self.arrhenius.alpha.value - arrhenius.alpha.value), 4) == 0
        assert round(abs(self.arrhenius.E0.value - arrhenius.E0.value), 4) == 0
        assert self.arrhenius.E0.units == arrhenius.E0.units
        assert round(abs(self.arrhenius.Tmin.value - arrhenius.Tmin.value), 4) == 0
        assert self.arrhenius.Tmin.units == arrhenius.Tmin.units
        assert round(abs(self.arrhenius.Tmax.value - arrhenius.Tmax.value), 4) == 0
        assert self.arrhenius.Tmax.units == arrhenius.Tmax.units
        assert self.arrhenius.comment == arrhenius.comment

    def test_change_rate(self):
        """
        Test the ArrheniusEP.change_rate() method.
        """
        Tlist = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500])
        k0list = np.array([self.arrhenius.get_rate_coefficient(T) for T in Tlist])
        self.arrhenius.change_rate(2)
        for T, kexp in zip(Tlist, k0list):
            kact = self.arrhenius.get_rate_coefficient(T)
            assert abs(2 * kexp - kact) < 1e-6 * kexp


class TestArrheniusBM:
    """
    Contains unit tests of the :class:`ArrheniusBM` class.
    """

    def setup_method(self):
        self.A = 8.00037e12
        self.n = 0.391734
        self.w0 = 798000
        self.E0 = 116249.32617478925
        self.Tmin = 300.0
        self.Tmax = 2000.0
        self.comment = "rxn001084"
        self.arrhenius_bm = ArrheniusBM(
            A=(self.A, "s^-1"),
            n=self.n,
            w0=(self.w0, "J/mol"),
            E0=(self.E0, "J/mol"),
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment=self.comment,
        )

        self.rsmi = "NC(=NC=O)O"
        self.psmi = "O=CNC(=O)N"
        self.arrhenius = Arrhenius(
            A=(8.00037e12, "s^-1"),
            n=0.391734,
            Ea=(94.5149, "kJ/mol"),
            T0=(1, "K"),
            Tmin=(300, "K"),
            Tmax=(2000, "K"),
            comment="""Fitted to 50 data points; dA = *|/ 1.18377, dn = +|- 0.0223855, dEa = +|- 0.115431 kJ/mol""",
        )

        self.r_thermo = NASA(
            polynomials=[
                NASAPolynomial(
                    coeffs=[
                        3.90453,
                        0.0068491,
                        0.000125755,
                        -2.92973e-07,
                        2.12971e-10,
                        -45444.2,
                        10.0669,
                    ],
                    Tmin=(10, "K"),
                    Tmax=(433.425, "K"),
                ),
                NASAPolynomial(
                    coeffs=[
                        2.09778,
                        0.0367646,
                        -2.36023e-05,
                        7.24527e-09,
                        -8.51275e-13,
                        -45412,
                        15.8381,
                    ],
                    Tmin=(433.425, "K"),
                    Tmax=(3000, "K"),
                ),
            ],
            Tmin=(10, "K"),
            Tmax=(3000, "K"),
            E0=(-377.851, "kJ/mol"),
            Cp0=(33.2579, "J/(mol*K)"),
            CpInf=(232.805, "J/(mol*K)"),
            comment="""Thermo library: Spiekermann_refining_elementary_reactions""",
        )
        self.p_thermo = NASA(
            polynomials=[
                NASAPolynomial(
                    coeffs=[
                        3.88423,
                        0.00825528,
                        0.000133399,
                        -3.31802e-07,
                        2.52823e-10,
                        -51045.1,
                        10.3937,
                    ],
                    Tmin=(10, "K"),
                    Tmax=(428.701, "K"),
                ),
                NASAPolynomial(
                    coeffs=[
                        2.89294,
                        0.0351772,
                        -2.26349e-05,
                        7.00331e-09,
                        -8.2982e-13,
                        -51122.5,
                        12.4424,
                    ],
                    Tmin=(428.701, "K"),
                    Tmax=(3000, "K"),
                ),
            ],
            Tmin=(10, "K"),
            Tmax=(3000, "K"),
            E0=(-424.419, "kJ/mol"),
            Cp0=(33.2579, "J/(mol*K)"),
            CpInf=(232.805, "J/(mol*K)"),
            comment="""Thermo library: Spiekermann_refining_elementary_reactions""",
        )
        CF2 = Species().from_adjacency_list(
            """
            1 F u0 p3 c0 {2,S}
            2 C u0 p1 c0 {1,S} {3,S}
            3 F u0 p3 c0 {2,S}
            """
        )
        CF2.thermo = NASA(
            polynomials=[
                    NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,-24340.7,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), 
                    NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,-25190.9,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))
                ], 
                Tmin=(298,'K'), Tmax=(3000,'K'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'),
                comment="""Thermo library: halogens"""
            )
        C2H6 = Species(smiles="CC")
        C2H6.thermo = ThermoData(
            Tdata = ([300,400,500,600,800,1000,1500],'K'),
            Cpdata = ([12.565,15.512,18.421,21.059,25.487,28.964,34.591],'cal/(mol*K)','+|-',[0.8,1.1,1.3,1.4,1.5,1.5,1.2]),
            H298 = (-20.028,'kcal/mol','+|-',0.1),
            S298 = (54.726,'cal/(mol*K)','+|-',0.6),
            comment="""Thermo library: DFT_QCI_thermo"""
            )
        CH3CF2CH3 = Species(smiles="CC(F)(F)C")
        CH3CF2CH3.thermo = NASA(
            polynomials = [
                NASAPolynomial(coeffs=[3.89769,0.00706735,0.000140168,-3.37628e-07,2.51812e-10,-68682.1,8.74321], Tmin=(10,'K'), Tmax=(436.522,'K')),
                NASAPolynomial(coeffs=[2.78849,0.0356982,-2.16715e-05,6.45057e-09,-7.47989e-13,-68761.2,11.1597], Tmin=(436.522,'K'), Tmax=(3000,'K')),
            ],
            Tmin = (10,'K'), Tmax = (3000,'K'), Cp0 = (33.2579,'J/(mol*K)'), CpInf = (249.434,'J/(mol*K)'),
            comment="""Thermo library: CHOF_G4"""
            )
        kinetics = Arrhenius(A=(0.222791,'cm^3/(mol*s)'), n=3.59921, Ea=(320.496,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), Tmax=(2500,'K'), comment="""Training Rxn 54 for 1,2_Insertion_carbene""")
        self.reaction = Reaction(reactants=[CF2,C2H6],products=[CH3CF2CH3],kinetics=kinetics)
        self.reaction_w0 = 519000 # J/mol

    def test_a_factor(self):
        """
        Test that the ArrheniusBM A property was properly set.
        """
        assert abs(self.arrhenius_bm.A.value_si - self.A) < 1e0

    def test_n(self):
        """
        Test that the ArrheniusBM n property was properly set.
        """
        assert round(abs(self.arrhenius_bm.n.value_si - self.n), 6) == 0

    def test_w0(self):
        """
        Test that the ArrheniusBM w0 property was properly set.
        """
        assert round(abs(self.arrhenius_bm.w0.value_si - self.w0), 6) == 0

    def test_e0(self):
        """
        Test that the ArrheniusBM E0 property was properly set.
        """
        assert round(abs(self.arrhenius_bm.E0.value_si - self.E0), 6) == 0

    def test_temperature_min(self):
        """
        Test that the ArrheniusBM Tmin property was properly set.
        """
        assert round(abs(self.arrhenius_bm.Tmin.value_si - self.Tmin), 6) == 0

    def test_temperature_max(self):
        """
        Test that the ArrheniusBM Tmax property was properly set.
        """
        assert round(abs(self.arrhenius_bm.Tmax.value_si - self.Tmax), 6) == 0

    def test_is_temperature_valid(self):
        """
        Test the ArrheniusBM.is_temperature_valid() method.
        """
        Tdata = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        validdata = np.array([False, True, True, True, True, True, True, True, True, True], bool)
        for T, valid in zip(Tdata, validdata):
            valid0 = self.arrhenius_bm.is_temperature_valid(T)
            assert valid0 == valid

    def test_fit_to_data(self):
        """
        Test the ArrheniusBM.fit_to_data() method.
        """
        reactant = Molecule(smiles=self.rsmi)
        product = Molecule(smiles=self.psmi)
        reaction = Reaction(
            reactants=[
                Species(
                    molecule=[reactant],
                    thermo=self.r_thermo,
                )
            ],
            products=[Species(molecule=[product], thermo=self.p_thermo)],
            kinetics=self.arrhenius,
        )
        Tdata = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500])
        kdata = np.array([reaction.kinetics.get_rate_coefficient(T) for T in Tdata])
        arrhenius_bm = ArrheniusBM().fit_to_reactions([reaction], w0=self.w0)
        assert abs(arrhenius_bm.A.value_si - self.arrhenius_bm.A.value_si) < 1.5e1
        assert round(abs(arrhenius_bm.n.value_si - self.arrhenius_bm.n.value_si), 1) == 0, 4
        assert round(abs(arrhenius_bm.E0.value_si - self.arrhenius_bm.E0.value_si), 1) == 0
        arrhenius = arrhenius_bm.to_arrhenius(reaction.get_enthalpy_of_reaction(298))
        for T, k in zip(Tdata, kdata):
            assert abs(k - arrhenius.get_rate_coefficient(T)) < 1e-6 * k

        # A second check, with a different reaction
        arrhenius_bm = ArrheniusBM().fit_to_reactions([self.reaction], w0=self.reaction_w0)
        arrhenius = arrhenius_bm.to_arrhenius(self.reaction.get_enthalpy_of_reaction(298))
        kdata = np.array([self.reaction.kinetics.get_rate_coefficient(T) for T in Tdata])
        for T, k in zip(Tdata, kdata):
            assert abs(k - arrhenius.get_rate_coefficient(T)) < 1e-6 * k

    def test_get_activation_energy(self):
        """
        Test the ArrheniusBM.get_activation_energy() method.
        """
        Hrxn = -44000  # J/mol
        Ea = self.arrhenius_bm.get_activation_energy(Hrxn)
        w = self.w0
        E0 = self.E0
        Vp = 2 * w * (w + E0)/(w - E0)
        Ea_exp = (w + Hrxn/2) * (Vp - 2*w + Hrxn)**2 / (Vp*Vp - 4*w*w + Hrxn*Hrxn)

        assert abs(Ea - Ea_exp) < 1e1


class TestPDepArrhenius:
    """
    Contains unit tests of the :class:`PDepArrhenius` class.
    """

    def setup_method(self):
        self.arrhenius0 = Arrhenius(
            A=(1.0e6, "s^-1"),
            n=1.0,
            Ea=(10.0, "kJ/mol"),
            T0=(300.0, "K"),
            Tmin=(300.0, "K"),
            Tmax=(2000.0, "K"),
            comment="""This data is completely made up""",
        )
        self.arrhenius1 = Arrhenius(
            A=(1.0e12, "s^-1"),
            n=1.0,
            Ea=(20.0, "kJ/mol"),
            T0=(300.0, "K"),
            Tmin=(300.0, "K"),
            Tmax=(2000.0, "K"),
            comment="""This data is completely made up""",
        )
        self.pressures = np.array([0.1, 10.0])
        self.arrhenius = [self.arrhenius0, self.arrhenius1]
        self.Tmin = 300.0
        self.Tmax = 2000.0
        self.Pmin = 0.1
        self.Pmax = 10.0
        self.comment = """This data is completely made up"""
        self.kinetics = PDepArrhenius(
            pressures=(self.pressures, "bar"),
            arrhenius=self.arrhenius,
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            Pmin=(self.Pmin, "bar"),
            Pmax=(self.Pmax, "bar"),
            comment=self.comment,
        )

    def test_pressures(self):
        """
        Test that the PDepArrhenius pressures property was properly set.
        """
        assert len(self.kinetics.pressures.value_si) == 2
        for i in range(2):
            assert round(abs(self.kinetics.pressures.value_si[i] * 1e-5 - self.pressures[i]), 4) == 0

    def test_arrhenius(self):
        """
        Test that the PDepArrhenius arrhenius property was properly set.
        """
        assert len(self.kinetics.arrhenius) == 2
        for i in range(2):
            assert abs(self.kinetics.arrhenius[i].A.value - self.arrhenius[i].A.value) < 1e0
            assert self.kinetics.arrhenius[i].A.units == self.arrhenius[i].A.units
            assert round(abs(self.kinetics.arrhenius[i].n.value - self.arrhenius[i].n.value), 4) == 0
            assert round(abs(self.kinetics.arrhenius[i].Ea.value - self.arrhenius[i].Ea.value), 4) == 0
            assert self.kinetics.arrhenius[i].Ea.units == self.arrhenius[i].Ea.units
            assert round(abs(self.kinetics.arrhenius[i].T0.value - self.arrhenius[i].T0.value), 4) == 0
            assert self.kinetics.arrhenius[i].T0.units == self.arrhenius[i].T0.units
            assert round(abs(self.kinetics.arrhenius[i].Tmin.value - self.arrhenius[i].Tmin.value), 4) == 0
            assert self.kinetics.arrhenius[i].Tmin.units == self.arrhenius[i].Tmin.units
            assert round(abs(self.kinetics.arrhenius[i].Tmax.value - self.arrhenius[i].Tmax.value), 4) == 0
            assert self.kinetics.arrhenius[i].Tmax.units == self.arrhenius[i].Tmax.units
            assert self.kinetics.arrhenius[i].comment == self.arrhenius[i].comment

    def test_temperature_min(self):
        """
        Test that the PDepArrhenius Tmin property was properly set.
        """
        assert round(abs(self.kinetics.Tmin.value_si - self.Tmin), 6) == 0

    def test_temperature_max(self):
        """
        Test that the PDepArrhenius Tmax property was properly set.
        """
        assert round(abs(self.kinetics.Tmax.value_si - self.Tmax), 6) == 0

    def test_pressure_min(self):
        """
        Test that the PDepArrhenius Pmin property was properly set.
        """
        assert round(abs(self.kinetics.Pmin.value_si * 1e-5 - self.Pmin), 6) == 0

    def test_pressure_max(self):
        """
        Test that the PDepArrhenius Pmax property was properly set.
        """
        assert round(abs(self.kinetics.Pmax.value_si * 1e-5 - self.Pmax), 6) == 0

    def test_comment(self):
        """
        Test that the PDepArrhenius comment property was properly set.
        """
        assert self.kinetics.comment == self.comment

    def test_is_pressure_dependent(self):
        """
        Test the PDepArrhenius.is_pressure_dependent() method.
        """
        assert self.kinetics.is_pressure_dependent()

    def test_get_rate_coefficient(self):
        """
        Test the PDepArrhenius.get_rate_coefficient() method.
        """
        P = 1e4
        for T in [
            300,
            400,
            500,
            600,
            700,
            800,
            900,
            1000,
            1100,
            1200,
            1300,
            1400,
            1500,
        ]:
            k0 = self.kinetics.get_rate_coefficient(T, P)
            k1 = self.arrhenius0.get_rate_coefficient(T)
            assert abs(k0 - k1) < 1e-6 * k1
        P = 1e6
        for T in [
            300,
            400,
            500,
            600,
            700,
            800,
            900,
            1000,
            1100,
            1200,
            1300,
            1400,
            1500,
        ]:
            k0 = self.kinetics.get_rate_coefficient(T, P)
            k1 = self.arrhenius1.get_rate_coefficient(T)
            assert abs(k0 - k1) < 1e-6 * k1
        P = 1e5
        for T in [
            300,
            400,
            500,
            600,
            700,
            800,
            900,
            1000,
            1100,
            1200,
            1300,
            1400,
            1500,
        ]:
            k0 = self.kinetics.get_rate_coefficient(T, P)
            k1 = math.sqrt(self.arrhenius0.get_rate_coefficient(T) * self.arrhenius1.get_rate_coefficient(T))
            assert abs(k0 - k1) < 1e-6 * k1

    def test_fit_to_data(self):
        """
        Test the PDepArrhenius.fit_to_data() method.
        """
        Tdata = np.array(
            [300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500],
            float,
        )
        Pdata = np.array([1e4, 3e4, 1e5, 3e5, 1e6], float)
        kdata = np.zeros([len(Tdata), len(Pdata)], float)
        for t in range(len(Tdata)):
            for p in range(len(Pdata)):
                kdata[t, p] = self.kinetics.get_rate_coefficient(Tdata[t], Pdata[p])
        kinetics = PDepArrhenius().fit_to_data(Tdata, Pdata, kdata, kunits="s^-1")
        for t in range(len(Tdata)):
            for p in range(len(Pdata)):
                assert abs(kinetics.get_rate_coefficient(Tdata[t], Pdata[p]) - kdata[t, p]) < 1e-6 * kdata[t, p]

    def test_pickle(self):
        """
        Test that a PDepArrhenius object can be successfully pickled and
        unpickled with no loss of information.
        """
        import pickle

        kinetics = pickle.loads(pickle.dumps(self.kinetics, -1))
        Narrh = 2
        assert len(self.kinetics.pressures.value) == Narrh
        assert len(kinetics.pressures.value) == Narrh
        assert len(self.kinetics.arrhenius) == Narrh
        assert len(kinetics.arrhenius) == Narrh
        for i in range(Narrh):
            assert round(abs(self.kinetics.pressures.value[i] - kinetics.pressures.value[i]), 4) == 0
            assert abs(self.kinetics.arrhenius[i].A.value - kinetics.arrhenius[i].A.value) < 1e0
            assert self.kinetics.arrhenius[i].A.units == kinetics.arrhenius[i].A.units
            assert round(abs(self.kinetics.arrhenius[i].n.value - kinetics.arrhenius[i].n.value), 7) == 0
            assert round(abs(self.kinetics.arrhenius[i].T0.value - kinetics.arrhenius[i].T0.value), 4) == 0
            assert self.kinetics.arrhenius[i].T0.units == kinetics.arrhenius[i].T0.units
            assert round(abs(self.kinetics.arrhenius[i].Ea.value - kinetics.arrhenius[i].Ea.value), 4) == 0
            assert self.kinetics.arrhenius[i].Ea.units == kinetics.arrhenius[i].Ea.units
        assert round(abs(self.kinetics.Tmin.value - kinetics.Tmin.value), 4) == 0
        assert self.kinetics.Tmin.units == kinetics.Tmin.units
        assert round(abs(self.kinetics.Tmax.value - kinetics.Tmax.value), 4) == 0
        assert self.kinetics.Tmax.units == kinetics.Tmax.units
        assert round(abs(self.kinetics.Pmin.value - kinetics.Pmin.value), 4) == 0
        assert self.kinetics.Pmin.units == kinetics.Pmin.units
        assert round(abs(self.kinetics.Pmax.value - kinetics.Pmax.value), 4) == 0
        assert self.kinetics.Pmax.units == kinetics.Pmax.units
        assert self.kinetics.comment == kinetics.comment

    def test_repr(self):
        """
        Test that a PDepArrhenius object can be successfully reconstructed
        from its repr() output with no loss of information.
        """
        namespace = {}
        exec("kinetics = {0!r}".format(self.kinetics), globals(), namespace)
        assert "kinetics" in namespace
        kinetics = namespace["kinetics"]
        Narrh = 2
        assert len(self.kinetics.pressures.value) == Narrh
        assert len(kinetics.pressures.value) == Narrh
        assert len(self.kinetics.arrhenius) == Narrh
        assert len(kinetics.arrhenius) == Narrh
        for i in range(Narrh):
            assert round(abs(self.kinetics.pressures.value[i] - kinetics.pressures.value[i]), 4) == 0
            assert abs(self.kinetics.arrhenius[i].A.value - kinetics.arrhenius[i].A.value) < 1e0
            assert self.kinetics.arrhenius[i].A.units == kinetics.arrhenius[i].A.units
            assert round(abs(self.kinetics.arrhenius[i].n.value - kinetics.arrhenius[i].n.value), 7) == 0
            assert round(abs(self.kinetics.arrhenius[i].T0.value - kinetics.arrhenius[i].T0.value), 4) == 0
            assert self.kinetics.arrhenius[i].T0.units == kinetics.arrhenius[i].T0.units
            assert round(abs(self.kinetics.arrhenius[i].Ea.value - kinetics.arrhenius[i].Ea.value), 4) == 0
            assert self.kinetics.arrhenius[i].Ea.units == kinetics.arrhenius[i].Ea.units
        assert round(abs(self.kinetics.Tmin.value - kinetics.Tmin.value), 4) == 0
        assert self.kinetics.Tmin.units == kinetics.Tmin.units
        assert round(abs(self.kinetics.Tmax.value - kinetics.Tmax.value), 4) == 0
        assert self.kinetics.Tmax.units == kinetics.Tmax.units
        assert round(abs(self.kinetics.Pmin.value - kinetics.Pmin.value), 4) == 0
        assert self.kinetics.Pmin.units == kinetics.Pmin.units
        assert round(abs(self.kinetics.Pmax.value - kinetics.Pmax.value), 4) == 0
        assert self.kinetics.Pmax.units == kinetics.Pmax.units
        assert self.kinetics.comment == kinetics.comment

    def test_change_rate(self):
        """
        Test the PDepArrhenius.change_rate() method.
        """
        Tlist = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500])
        k0list = np.array([self.kinetics.get_rate_coefficient(T, 1e5) for T in Tlist])
        self.kinetics.change_rate(2)
        for T, kexp in zip(Tlist, k0list):
            kact = self.kinetics.get_rate_coefficient(T, 1e5)
            assert abs(2 * kexp - kact) < 1e-6 * kexp


class TestMultiArrhenius:
    """
    Contains unit tests of the :class:`MultiArrhenius` class.
    """

    def setup_method(self):
        self.Tmin = 350.0
        self.Tmax = 1500.0
        self.comment = "Comment"
        self.arrhenius = [
            Arrhenius(
                A=(9.3e-14, "cm^3/(molecule*s)"),
                n=0.0,
                Ea=(4740 * constants.R * 0.001, "kJ/mol"),
                T0=(1, "K"),
                Tmin=(self.Tmin, "K"),
                Tmax=(self.Tmax, "K"),
                comment=self.comment,
            ),
            Arrhenius(
                A=(1.4e-9, "cm^3/(molecule*s)"),
                n=0.0,
                Ea=(11200 * constants.R * 0.001, "kJ/mol"),
                T0=(1, "K"),
                Tmin=(self.Tmin, "K"),
                Tmax=(self.Tmax, "K"),
                comment=self.comment,
            ),
        ]
        self.kinetics = MultiArrhenius(
            arrhenius=self.arrhenius,
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment=self.comment,
        )
        self.single_kinetics = MultiArrhenius(
            arrhenius=self.arrhenius[:1],
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment=self.comment,
        )

    def test_arrhenius(self):
        """
        Test that the MultiArrhenius A property was properly set.
        """
        assert self.kinetics.arrhenius == self.arrhenius

    def test_temperature_min(self):
        """
        Test that the MultiArrhenius Tmin property was properly set.
        """
        assert round(abs(self.kinetics.Tmin.value_si - self.Tmin), 6) == 0

    def test_temperature_max(self):
        """
        Test that the MultiArrhenius Tmax property was properly set.
        """
        assert round(abs(self.kinetics.Tmax.value_si - self.Tmax), 6) == 0

    def test_comment(self):
        """
        Test that the MultiArrhenius comment property was properly set.
        """
        assert self.kinetics.comment == self.comment

    def test_is_temperature_valid(self):
        """
        Test the MultiArrhenius.is_temperature_valid() method.
        """
        Tdata = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        validdata = np.array([False, True, True, True, True, True, True, False, False, False], bool)
        for T, valid in zip(Tdata, validdata):
            valid0 = self.kinetics.is_temperature_valid(T)
            assert valid0 == valid

    def test_get_rate_coefficient(self):
        """
        Test the MultiArrhenius.get_rate_coefficient() method.
        """
        Tlist = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        kexplist = np.array(
            [
                2.85400e-06,
                4.00384e-01,
                2.73563e01,
                8.50699e02,
                1.20181e04,
                7.56312e04,
                2.84724e05,
                7.71702e05,
                1.67743e06,
                3.12290e06,
            ]
        )
        for T, kexp in zip(Tlist, kexplist):
            kact = self.kinetics.get_rate_coefficient(T)
            assert abs(kexp - kact) < 1e-4 * kexp

    def test_pickle(self):
        """
        Test that a MultiArrhenius object can be pickled and unpickled with no loss
        of information.
        """
        import pickle

        kinetics = pickle.loads(pickle.dumps(self.kinetics, -1))
        assert len(self.kinetics.arrhenius) == len(kinetics.arrhenius)
        for arrh0, arrh in zip(self.kinetics.arrhenius, kinetics.arrhenius):
            assert abs(arrh0.A.value - arrh.A.value) < 1e-18
            assert arrh0.A.units == arrh.A.units
            assert round(abs(arrh0.n.value - arrh.n.value), 4) == 0
            assert round(abs(arrh0.Ea.value - arrh.Ea.value), 4) == 0
            assert arrh0.Ea.units == arrh.Ea.units
            assert round(abs(arrh0.T0.value - arrh.T0.value), 4) == 0
            assert arrh0.T0.units == arrh.T0.units
        assert round(abs(self.kinetics.Tmin.value - kinetics.Tmin.value), 4) == 0
        assert self.kinetics.Tmin.units == kinetics.Tmin.units
        assert round(abs(self.kinetics.Tmax.value - kinetics.Tmax.value), 4) == 0
        assert self.kinetics.Tmax.units == kinetics.Tmax.units
        assert self.kinetics.comment == kinetics.comment

    def test_repr(self):
        """
        Test that a MultiArrhenius object can be reconstructed from its repr()
        output with no loss of information.
        """
        namespace = {}
        exec("kinetics = {0!r}".format(self.kinetics), globals(), namespace)
        assert "kinetics" in namespace
        kinetics = namespace["kinetics"]
        assert len(self.kinetics.arrhenius) == len(kinetics.arrhenius)
        for arrh0, arrh in zip(self.kinetics.arrhenius, kinetics.arrhenius):
            assert abs(arrh0.A.value - arrh.A.value) < 1e-18
            assert arrh0.A.units == arrh.A.units
            assert round(abs(arrh0.n.value - arrh.n.value), 4) == 0
            assert round(abs(arrh0.Ea.value - arrh.Ea.value), 4) == 0
            assert arrh0.Ea.units == arrh.Ea.units
            assert round(abs(arrh0.T0.value - arrh.T0.value), 4) == 0
            assert arrh0.T0.units == arrh.T0.units
        assert round(abs(self.kinetics.Tmin.value - kinetics.Tmin.value), 4) == 0
        assert self.kinetics.Tmin.units == kinetics.Tmin.units
        assert round(abs(self.kinetics.Tmax.value - kinetics.Tmax.value), 4) == 0
        assert self.kinetics.Tmax.units == kinetics.Tmax.units
        assert self.kinetics.comment == kinetics.comment

    def test_to_arrhenius(self):
        """
        Test that we can convert to an Arrhenius
        """
        answer = self.single_kinetics.arrhenius[0]
        fitted = self.single_kinetics.to_arrhenius()

        assert abs(fitted.A.value_si - answer.A.value_si) < 1e0
        assert round(abs(fitted.n.value_si - answer.n.value_si), 1) == 0, 4
        assert round(abs(fitted.Ea.value_si - answer.Ea.value_si), 2) == 0
        assert round(abs(fitted.T0.value_si - answer.T0.value_si), 4) == 0

    def test_to_arrhenius_temperature_range(self):
        """
        Test the to_arrhenius temperature range is set correctly.
        """
        answer = self.single_kinetics.arrhenius[0]
        fitted = self.single_kinetics.to_arrhenius(Tmin=800, Tmax=1200)
        assert round(abs(fitted.Tmin.value_si - 800.0), 7) == 0
        assert round(abs(fitted.Tmax.value_si - 1200.0), 7) == 0
        for T in [800, 1000, 1200]:
            assert round(abs(fitted.get_rate_coefficient(T) / answer.get_rate_coefficient(T) - 1.0), 7) == 0

    def test_to_arrhenius_multiple(self):
        """
        Test the to_arrhenius fitting multiple kinetics over a small range, see if we're within 5% at a few points
        """
        answer = self.kinetics
        fitted = self.kinetics.to_arrhenius(Tmin=800, Tmax=1200)
        assert round(abs(fitted.Tmin.value_si - 800.0), 7) == 0
        assert round(abs(fitted.Tmax.value_si - 1200.0), 7) == 0
        for T in [800, 1000, 1200]:
            assert abs(fitted.get_rate_coefficient(T) / answer.get_rate_coefficient(T) - 1.0) < 0.05

    def test_change_rate(self):
        """
        Test the MultiArrhenius.change_rate() method.
        """
        Tlist = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500])
        k0list = np.array([self.kinetics.get_rate_coefficient(T) for T in Tlist])
        self.kinetics.change_rate(2)
        for T, kexp in zip(Tlist, k0list):
            kact = self.kinetics.get_rate_coefficient(T)
            assert abs(2 * kexp - kact) < 1e-6 * kexp


class TestMultiPDepArrhenius:
    """
    Contains unit tests of the :class:`MultiPDepArrhenius` class.
    """

    def setup_method(self):
        self.Tmin = 350.0
        self.Tmax = 1500.0
        self.Pmin = 1e-1
        self.Pmax = 1e1
        self.pressures = np.array([1e-1, 1e1])
        self.comment = "CH3 + C2H6 <=> CH4 + C2H5 (Baulch 2005)"
        self.arrhenius = [
            PDepArrhenius(
                pressures=(self.pressures, "bar"),
                arrhenius=[
                    Arrhenius(
                        A=(9.3e-16, "cm^3/(molecule*s)"),
                        n=0.0,
                        Ea=(4740 * constants.R * 0.001, "kJ/mol"),
                        T0=(1, "K"),
                        Tmin=(self.Tmin, "K"),
                        Tmax=(self.Tmax, "K"),
                        comment=self.comment,
                    ),
                    Arrhenius(
                        A=(9.3e-14, "cm^3/(molecule*s)"),
                        n=0.0,
                        Ea=(4740 * constants.R * 0.001, "kJ/mol"),
                        T0=(1, "K"),
                        Tmin=(self.Tmin, "K"),
                        Tmax=(self.Tmax, "K"),
                        comment=self.comment,
                    ),
                ],
                Tmin=(self.Tmin, "K"),
                Tmax=(self.Tmax, "K"),
                Pmin=(self.Pmin, "bar"),
                Pmax=(self.Pmax, "bar"),
                comment=self.comment,
            ),
            PDepArrhenius(
                pressures=(self.pressures, "bar"),
                arrhenius=[
                    Arrhenius(
                        A=(1.4e-11, "cm^3/(molecule*s)"),
                        n=0.0,
                        Ea=(11200 * constants.R * 0.001, "kJ/mol"),
                        T0=(1, "K"),
                        Tmin=(self.Tmin, "K"),
                        Tmax=(self.Tmax, "K"),
                        comment=self.comment,
                    ),
                    Arrhenius(
                        A=(1.4e-9, "cm^3/(molecule*s)"),
                        n=0.0,
                        Ea=(11200 * constants.R * 0.001, "kJ/mol"),
                        T0=(1, "K"),
                        Tmin=(self.Tmin, "K"),
                        Tmax=(self.Tmax, "K"),
                        comment=self.comment,
                    ),
                ],
                Tmin=(self.Tmin, "K"),
                Tmax=(self.Tmax, "K"),
                Pmin=(self.Pmin, "bar"),
                Pmax=(self.Pmax, "bar"),
                comment=self.comment,
            ),
        ]
        self.kinetics = MultiPDepArrhenius(
            arrhenius=self.arrhenius,
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            Pmin=(self.Pmin, "bar"),
            Pmax=(self.Pmax, "bar"),
            comment=self.comment,
        )

    def test_arrhenius(self):
        """
        Test that the MultiPDepArrhenius arrhenius property was properly set.
        """
        assert self.kinetics.arrhenius == self.arrhenius

    def test_temperature_min(self):
        """
        Test that the MultiPDepArrhenius Tmin property was properly set.
        """
        assert round(abs(self.kinetics.Tmin.value_si - self.Tmin), 6) == 0

    def test_temperature_max(self):
        """
        Test that the MultiPDepArrhenius Tmax property was properly set.
        """
        assert round(abs(self.kinetics.Tmax.value_si - self.Tmax), 6) == 0

    def test_pressure_min(self):
        """
        Test that the MultiPDepArrhenius Pmin property was properly set.
        """
        assert round(abs(self.kinetics.Pmin.value_si * 1e-5 - self.Pmin), 6) == 0

    def test_pressure_max(self):
        """
        Test that the MultiPDepArrhenius Pmax property was properly set.
        """
        assert round(abs(self.kinetics.Pmax.value_si * 1e-5 - self.Pmax), 6) == 0

    def test_comment(self):
        """
        Test that the MultiPDepArrhenius comment property was properly set.
        """
        assert self.kinetics.comment == self.comment

    def test_is_temperature_valid(self):
        """
        Test the MultiPDepArrhenius.is_temperature_valid() method.
        """
        Tdata = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        validdata = np.array([False, True, True, True, True, True, True, False, False, False], bool)
        for T, valid in zip(Tdata, validdata):
            valid0 = self.kinetics.is_temperature_valid(T)
            assert valid0 == valid

    def test_is_pressure_valid(self):
        """
        Test the MultiPDepArrhenius.is_pressure_valid() method.
        """
        Pdata = np.array([1e3, 1e4, 1e5, 1e6, 1e7])
        validdata = np.array([False, True, True, True, False], bool)
        for P, valid in zip(Pdata, validdata):
            valid0 = self.kinetics.is_pressure_valid(P)
            assert valid0 == valid

    def test_get_rate_coefficient(self):
        """
        Test the MultiPDepArrhenius.get_rate_coefficient() method.
        """
        Tlist = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        Plist = np.array([1e4, 1e5, 1e6])
        kexplist = np.array(
            [
                [
                    2.85400e-08,
                    4.00384e-03,
                    2.73563e-01,
                    8.50699e00,
                    1.20181e02,
                    7.56312e02,
                    2.84724e03,
                    7.71702e03,
                    1.67743e04,
                    3.12290e04,
                ],
                [
                    2.85400e-07,
                    4.00384e-02,
                    2.73563e00,
                    8.50699e01,
                    1.20181e03,
                    7.56312e03,
                    2.84724e04,
                    7.71702e04,
                    1.67743e05,
                    3.12290e05,
                ],
                [
                    2.85400e-06,
                    4.00384e-01,
                    2.73563e01,
                    8.50699e02,
                    1.20181e04,
                    7.56312e04,
                    2.84724e05,
                    7.71702e05,
                    1.67743e06,
                    3.12290e06,
                ],
            ]
        ).T
        for i in range(Tlist.shape[0]):
            for j in range(Plist.shape[0]):
                kexp = kexplist[i, j]
                kact = self.kinetics.get_rate_coefficient(Tlist[i], Plist[j])
                assert abs(kexp - kact) < 1e-4 * kexp

    def test_get_rate_coefficient_diff_plist(self):
        """
        Test the MultiPDepArrhenius.get_rate_coefficient() when plists are different.
        """
        # modify the MultiPDepArrhenius object with an additional entry
        pressures = np.array([1e-1, 1e-1, 1e1])
        self.kinetics.arrhenius[0].pressures = (pressures, "bar")
        self.kinetics.arrhenius[0].arrhenius.insert(0, self.kinetics.arrhenius[0].arrhenius[0])

        Tlist = np.array([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000])
        Plist = np.array([1e4, 1e5, 1e6])
        kexplist = np.array(
            [
                [
                    2.85400e-08,
                    4.00384e-03,
                    2.73563e-01,
                    8.50699e00,
                    1.20181e02,
                    7.56312e02,
                    2.84724e03,
                    7.71702e03,
                    1.67743e04,
                    3.12290e04,
                ],
                [
                    2.85400e-07,
                    4.00384e-02,
                    2.73563e00,
                    8.50699e01,
                    1.20181e03,
                    7.56312e03,
                    2.84724e04,
                    7.71702e04,
                    1.67743e05,
                    3.12290e05,
                ],
                [
                    2.85400e-06,
                    4.00384e-01,
                    2.73563e01,
                    8.50699e02,
                    1.20181e04,
                    7.56312e04,
                    2.84724e05,
                    7.71702e05,
                    1.67743e06,
                    3.12290e06,
                ],
            ]
        ).T
        for i in range(Tlist.shape[0]):
            for j in range(Plist.shape[0]):
                kexp = kexplist[i, j]
                kact = self.kinetics.get_rate_coefficient(Tlist[i], Plist[j])
                assert abs(kexp - kact) < 1e-4 * kexp

    def test_pickle(self):
        """
        Test that a MultiPDepArrhenius object can be pickled and unpickled with
        no loss of information.
        """
        import pickle

        kinetics = pickle.loads(pickle.dumps(self.kinetics, -1))
        assert len(self.kinetics.arrhenius) == len(kinetics.arrhenius)
        assert round(abs(self.kinetics.Tmin.value - kinetics.Tmin.value), 4) == 0
        assert self.kinetics.Tmin.units == kinetics.Tmin.units
        assert round(abs(self.kinetics.Tmax.value - kinetics.Tmax.value), 4) == 0
        assert self.kinetics.Tmax.units == kinetics.Tmax.units
        assert self.kinetics.comment == kinetics.comment

    def test_repr(self):
        """
        Test that a MultiPDepArrhenius object can be reconstructed from its
        repr() output with no loss of information.
        """
        namespace = {}
        exec("kinetics = {0!r}".format(self.kinetics), globals(), namespace)
        assert "kinetics" in namespace
        kinetics = namespace["kinetics"]
        assert len(self.kinetics.arrhenius) == len(kinetics.arrhenius)
        assert round(abs(self.kinetics.Tmin.value - kinetics.Tmin.value), 4) == 0
        assert self.kinetics.Tmin.units == kinetics.Tmin.units
        assert round(abs(self.kinetics.Tmax.value - kinetics.Tmax.value), 4) == 0
        assert self.kinetics.Tmax.units == kinetics.Tmax.units
        assert self.kinetics.comment == kinetics.comment

    def test_change_rate(self):
        """
        Test the PDepMultiArrhenius.change_rate() method.
        """
        Tlist = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500])
        k0list = np.array([self.kinetics.get_rate_coefficient(T, 1e5) for T in Tlist])
        self.kinetics.change_rate(2)
        for T, kexp in zip(Tlist, k0list):
            kact = self.kinetics.get_rate_coefficient(T, 1e5)
            assert abs(2 * kexp - kact) < 1e-6 * kexp

    def test_generate_reverse_rate_coefficient(self):
        """
        Test ability to reverse a reaction rate.

        This is a real example from an imported chemkin file.
        """
        from rmgpy.species import Species
        from rmgpy.molecule import Molecule
        from rmgpy.data.kinetics import LibraryReaction
        from rmgpy.thermo import NASA, NASAPolynomial

        test_reaction = LibraryReaction(
            reactants=[
                Species(
                    label="C2H3",
                    thermo=NASA(
                        polynomials=[
                            NASAPolynomial(
                                coeffs=[
                                    3.12502,
                                    0.00235137,
                                    2.36803e-05,
                                    -3.35092e-08,
                                    1.39444e-11,
                                    34524.3,
                                    8.81538,
                                ],
                                Tmin=(200, "K"),
                                Tmax=(1000, "K"),
                            ),
                            NASAPolynomial(
                                coeffs=[
                                    4.37211,
                                    0.00746869,
                                    -2.64716e-06,
                                    4.22753e-10,
                                    -2.44958e-14,
                                    33805.2,
                                    0.428772,
                                ],
                                Tmin=(1000, "K"),
                                Tmax=(6000, "K"),
                            ),
                        ],
                        Tmin=(200, "K"),
                        Tmax=(6000, "K"),
                        E0=(285.696, "kJ/mol"),
                        Cp0=(33.2579, "J/mol/K"),
                        CpInf=(108.088, "J/mol/K"),
                        comment="""ATcT3E\nC2H3 <g> ATcT ver. 1.122, DHf298 = 296.91 ± 0.33 kJ/mol - fit JAN17""",
                    ),
                    molecule=[Molecule(smiles="[CH]=C")],
                    molecular_weight=(27.0452, "amu"),
                ),
                Species(
                    label="CH2O",
                    thermo=NASA(
                        polynomials=[
                            NASAPolynomial(
                                coeffs=[
                                    4.77187,
                                    -0.00976266,
                                    3.70122e-05,
                                    -3.76922e-08,
                                    1.31327e-11,
                                    -14379.8,
                                    0.696586,
                                ],
                                Tmin=(200, "K"),
                                Tmax=(1000, "K"),
                            ),
                            NASAPolynomial(
                                coeffs=[
                                    2.91333,
                                    0.0067004,
                                    -2.55521e-06,
                                    4.27795e-10,
                                    -2.44073e-14,
                                    -14462.2,
                                    7.43823,
                                ],
                                Tmin=(1000, "K"),
                                Tmax=(6000, "K"),
                            ),
                        ],
                        Tmin=(200, "K"),
                        Tmax=(6000, "K"),
                        E0=(-119.527, "kJ/mol"),
                        Cp0=(33.2579, "J/mol/K"),
                        CpInf=(83.1447, "J/mol/K"),
                        comment="""ATcT3E\nH2CO <g> ATcT ver. 1.122, DHf298 = -109.188 ± 0.099 kJ/mol - fit JAN17""",
                    ),
                    molecule=[Molecule(smiles="C=O")],
                    molecular_weight=(30.026, "amu"),
                ),
            ],
            products=[
                Species(
                    label="C2H4",
                    thermo=NASA(
                        polynomials=[
                            NASAPolynomial(
                                coeffs=[
                                    3.65151,
                                    -0.00535067,
                                    5.16486e-05,
                                    -6.36869e-08,
                                    2.50743e-11,
                                    5114.51,
                                    5.38561,
                                ],
                                Tmin=(200, "K"),
                                Tmax=(1000, "K"),
                            ),
                            NASAPolynomial(
                                coeffs=[
                                    4.14446,
                                    0.0102648,
                                    -3.61247e-06,
                                    5.74009e-10,
                                    -3.39296e-14,
                                    4190.59,
                                    -1.14778,
                                ],
                                Tmin=(1000, "K"),
                                Tmax=(6000, "K"),
                            ),
                        ],
                        Tmin=(200, "K"),
                        Tmax=(6000, "K"),
                        E0=(42.06, "kJ/mol"),
                        Cp0=(33.2579, "J/mol/K"),
                        CpInf=(133.032, "J/mol/K"),
                        comment="""ATcT3E\nC2H4 <g> ATcT ver. 1.122, DHf298 = 52.45 ± 0.13 kJ/mol - fit JAN17""",
                    ),
                    molecule=[Molecule(smiles="C=C")],
                    molecular_weight=(28.0532, "amu"),
                ),
                Species(
                    label="HCO",
                    thermo=NASA(
                        polynomials=[
                            NASAPolynomial(
                                coeffs=[
                                    3.97075,
                                    -0.00149122,
                                    9.54042e-06,
                                    -8.8272e-09,
                                    2.67645e-12,
                                    3842.03,
                                    4.4466,
                                ],
                                Tmin=(200, "K"),
                                Tmax=(1000, "K"),
                            ),
                            NASAPolynomial(
                                coeffs=[
                                    3.85781,
                                    0.00264114,
                                    -7.44177e-07,
                                    1.23313e-10,
                                    -8.88959e-15,
                                    3616.43,
                                    3.92451,
                                ],
                                Tmin=(1000, "K"),
                                Tmax=(6000, "K"),
                            ),
                        ],
                        Tmin=(200, "K"),
                        Tmax=(6000, "K"),
                        E0=(32.0237, "kJ/mol"),
                        Cp0=(33.2579, "J/mol/K"),
                        CpInf=(58.2013, "J/mol/K"),
                        comment="""HCO <g> ATcT ver. 1.122, DHf298 = 41.803 ± 0.099 kJ/mol - fit JAN17""",
                    ),
                    molecule=[Molecule(smiles="[CH]=O")],
                    molecular_weight=(29.018, "amu"),
                ),
            ],
            kinetics=MultiPDepArrhenius(
                arrhenius=[
                    PDepArrhenius(
                        pressures=([0.001, 0.01, 0.1, 1, 10, 100, 1000], "atm"),
                        arrhenius=[
                            Arrhenius(
                                A=(1.1e07, "cm^3/(mol*s)"),
                                n=1.09,
                                Ea=(1807, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(2.5e07, "cm^3/(mol*s)"),
                                n=0.993,
                                Ea=(1995, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(2.5e08, "cm^3/(mol*s)"),
                                n=0.704,
                                Ea=(2596, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(1.4e10, "cm^3/(mol*s)"),
                                n=0.209,
                                Ea=(3934, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(3.5e13, "cm^3/(mol*s)"),
                                n=-0.726,
                                Ea=(6944, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(3.3e14, "cm^3/(mol*s)"),
                                n=-0.866,
                                Ea=(10966, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(17, "cm^3/(mol*s)"),
                                n=3.17,
                                Ea=(9400, "cal/mol"),
                                T0=(1, "K"),
                            ),
                        ],
                    ),
                    PDepArrhenius(
                        pressures=([0.001, 0.01, 0.1, 1, 10, 100, 1000], "atm"),
                        arrhenius=[
                            Arrhenius(
                                A=(-2.3e16, "cm^3/(mol*s)"),
                                n=-1.269,
                                Ea=(20617, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(-5.2e16, "cm^3/(mol*s)"),
                                n=-1.366,
                                Ea=(20805, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(-1.5e18, "cm^3/(mol*s)"),
                                n=-1.769,
                                Ea=(22524, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(-8.5e19, "cm^3/(mol*s)"),
                                n=-2.264,
                                Ea=(23862, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(-4.4e23, "cm^3/(mol*s)"),
                                n=-3.278,
                                Ea=(27795, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(-4.2e24, "cm^3/(mol*s)"),
                                n=-3.418,
                                Ea=(31817, "cal/mol"),
                                T0=(1, "K"),
                            ),
                            Arrhenius(
                                A=(-2.1e11, "cm^3/(mol*s)"),
                                n=0.618,
                                Ea=(30251, "cal/mol"),
                                T0=(1, "K"),
                            ),
                        ],
                    ),
                ]
            ),
            duplicate=True,
        )
        test_reaction.generate_reverse_rate_coefficient()


################################################################################
#
#   Plasma kinetics: TwoTemperaturePlasma, ElectronCollisionPlasma,
#   BadnellRRArrhenius and VoronovEIArrhenius.
#
################################################################################


def _plasma_yaml_path(filename):
    """
    Locate one of the plasma parameter tables (``badnell.yaml``, ``voronov.yaml``),
    which live in RMG-database under ``input/kinetics/``.

    The configured database is tried first. Since an RMG-Py feature worktree is
    paired with a same-suffixed RMG-database worktree (``RMG-Py-<feature>`` next to
    ``RMG-database-<feature>``), that sibling checkout is tried second, so these
    tests do not depend on how ``database.directory`` happens to be configured on
    the machine running them.

    Returns the path to the table, or None if no checkout provides it.
    """
    candidates = [os.path.join(settings["database.directory"], "kinetics", filename)]

    source_root = os.path.dirname(os.path.dirname(os.path.abspath(rmgpy.__file__)))
    parent, name = os.path.split(source_root)
    if name.startswith("RMG-Py"):
        sibling = "RMG-database" + name[len("RMG-Py"):]
        candidates.append(os.path.join(parent, sibling, "input", "kinetics", filename))

    for path in candidates:
        if os.path.isfile(path):
            return path
    return None


def _require_plasma_yaml(filename):
    """As :func:`_plasma_yaml_path`, but skip the test when the table is absent."""
    path = _plasma_yaml_path(filename)
    if path is None:
        pytest.skip("{0} not found in any RMG-database checkout".format(filename))
    return path


def _alpha_rr_cm3_per_molecule_s(Te, A_cms_per_molecule, B, T0, T1, C=None, T2=None):
    """
    Helper: compute the Badnell alpha_RR in cm^3/(molecule*s) directly from parameters, for test expectations.
    """
    if Te <= 0 or T0 <= 0 or T1 <= 0:
        raise ValueError("Te, T0, T1 must be > 0")
    Bstar = B + (C * np.exp(-T2 / Te) if (C is not None and T2 is not None) else 0.0)
    s0 = np.sqrt(Te / T0)
    s1 = np.sqrt(Te / T1)
    denom = s0 * (1.0 + s0) ** (1.0 - Bstar) * (1.0 + s1) ** (1.0 + Bstar)
    return A_cms_per_molecule / denom  # cm^3/(molecule*s)


def _alpha_ei_cm3_per_molecule_s(T, A, P, X, K, dE_eV):
    """
    Voronov (1997) per-particle rate coefficient:
        <sigma v>(Te_eV) = A * [ U^K * exp(-U) ] / [ (1 + P*sqrt(U)) * (X + U) ]
    with U = dE_eV / Te_eV, and Te_eV = k_B_eV_per_K * T.
    """
    kB_eV_per_K = 8.617333262145e-5  # eV/K
    Te_eV = kB_eV_per_K * float(T)
    U = dE_eV / Te_eV
    return A * ((1 + P * np.sqrt(U)) / (X + U)) * (U**K) * np.exp(-U)


class TestTwoTemperaturePlasma:
    """
    Tests for the TwoTemperaturePlasma kinetics class.

    The underlying functional form is

        k(T, Te) = A * Te^n
                   * exp(-Ea_g / (R * T))
                   * exp(Ea_e * (Te - T) / (R * T * Te))

    In the one-temperature fallback interface `get_rate_coefficient(T)`,
    we expect k(T, Te=T).
    """

    def setup_method(self):
        """
        Initialize a TwoTemperaturePlasma instance for testing.
        """
        self.A = 1.5e-10          # cm^3/(molecule*s) (Target value for assertions)
        self.n = 0.5
        self.Ea_g = 10.0          # kJ/mol
        self.Ea_e = 50.0          # kJ/mol
        self.Tmin = 300.0
        self.Tmax = 3000.0
        self.comment = "Test 2T plasma"

        self.plasma = TwoTemperaturePlasma(
            A=(self.A, "cm^3/(molecule*s)"),
            n=self.n,
            Ea_g=(self.Ea_g, "kJ/mol"),
            Ea_e=(self.Ea_e, "kJ/mol"),
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment=self.comment,
        )

        self.plasma_2 = TwoTemperaturePlasma(
            A=(self.A, "cm^3/(molecule*s)"),
            n=self.n,
            Ea_g=(self.Ea_g, "kJ/mol"),
            Ea_e=(self.Ea_e, "kJ/mol"),
            T0=(300, "K"),
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment=self.comment,
        )

    def test_a_factor(self):
        """
        A should be stored as a RateCoefficient with correct SI value.
        We supplied A in cm^3/(molecule*s).
        RMG stores SI as m^3/(mol*s).
        """
        # Calculate expected SI value:
        # 1.5e-10 cm^3/s/molecule * (1e-6 m^3/cm^3) * (6.022e23 molecule/mol)
        expected_si = self.A * 1.0e-6 * constants.Na
        # Compare stored .value_si against expected_si
        assert np.isclose(self.plasma.A.value_si, expected_si, rtol=1e-6)

    def test_n(self):
        """Electron temperature exponent n should be stored correctly."""
        assert abs(self.plasma.n.value_si - self.n) < 1e-12

    def test_ea_g(self):
        """Gas activation energy Ea_g should be stored correctly in kJ/mol."""
        assert abs(self.plasma.Ea_g.value_si * 0.001 - self.Ea_g) < 1e-6

    def test_ea_e(self):
        """Electron activation energy Ea_e should be stored correctly in kJ/mol."""
        assert abs(self.plasma.Ea_e.value_si * 0.001 - self.Ea_e) < 1e-6

    def test_t0(self):
        """T0 should default to 1.0 K if not specified."""
        assert abs(self.plasma.T0.value_si - 1) < 1e-6
        assert abs(self.plasma_2.T0.value_si - 300) < 1e-6

    def test_temperature_min(self):
        """Tmin should be stored correctly."""
        assert abs(self.plasma.Tmin.value_si - self.Tmin) < 1e-6

    def test_temperature_max(self):
        """Tmax should be stored correctly."""
        assert abs(self.plasma.Tmax.value_si - self.Tmax) < 1e-6

    def test_comment(self):
        """Comment should be preserved."""
        assert self.plasma.comment == self.comment

    def test_uses_electron_temperature_flag(self):
        """TwoTemperaturePlasma should advertise Te-dependence."""
        assert self.plasma.uses_electron_temperature is True

    def test_is_temperature_valid(self):
        """
        KineticsModel.is_temperature_valid(T) should respect Tmin/Tmax.
        """
        assert not self.plasma.is_temperature_valid(200.0)
        assert self.plasma.is_temperature_valid(1000.0)
        assert not self.plasma.is_temperature_valid(4000.0)

    def test_reduces_to_single_temp_when_Te_equals_Tgas(self):
        """
        The standard interface get_rate_coefficient(T) should match
        get_rate_coefficient_two_temp(T, Te) when Te == T.
        """
        for T in [300.0, 800.0, 1500.0, 2500.0]:
            k_single = self.plasma.get_rate_coefficient(T)
            k_two = self.plasma.get_rate_coefficient_two_temp(T, T)
            # Allow a tiny numerical tolerance
            assert abs(k_single - k_two) <= 1e-12 * max(1.0, abs(k_single))

    def test_rate_increases_with_electron_temperature(self):
        """
        For positive Ea_e, increasing Te at fixed T should increase the rate.
        """
        T = 1000.0  # K
        k_low = self.plasma.get_rate_coefficient_two_temp(T, 3000.0)
        k_high = self.plasma.get_rate_coefficient_two_temp(T, 8000.0)
        assert k_high > k_low

    def test_change_rate(self):
        """
        change_rate(factor) should scale the rate coefficient by that factor.
        """
        T = 1000.0
        Te = 5000.0
        k1 = self.plasma.get_rate_coefficient_two_temp(T, Te)
        self.plasma.change_rate(2.0)
        k2 = self.plasma.get_rate_coefficient_two_temp(T, Te)
        assert abs(k2 / k1 - 2.0) < 1e-6

    def test_pickle(self):
        """
        The object should roundtrip through pickle with key parameters intact.
        """
        plasma2 = pickle.loads(pickle.dumps(self.plasma, -1))

        # Check A (SI units) - Comparing stored values directly
        assert np.isclose(plasma2.A.value_si, self.plasma.A.value_si, rtol=1e-12)

        # Check Ea
        assert np.isclose(plasma2.Ea_g.value_si, self.plasma.Ea_g.value_si, rtol=1e-12)

        # Check rate evaluation
        k1 = self.plasma.get_rate_coefficient_two_temp(1000, 2000)
        k2 = plasma2.get_rate_coefficient_two_temp(1000, 2000)
        assert np.isclose(k1, k2, rtol=1e-12)

    def test_repr(self):
        """
        repr() should at least return a string without error.
        (We don't require it to be eval-roundtrippable.)
        """
        s = repr(self.plasma)
        assert isinstance(s, str)
        assert "TwoTemperaturePlasma" in s

    def test_cantera_yaml_structure(self):
        """
        Verify that the converted Cantera object contains the correct data structure.
        """
        ct = pytest.importorskip("cantera")

        ct_rate = self.plasma.to_cantera_kinetics()
        assert isinstance(ct_rate, ct.TwoTempPlasmaRate)

        if hasattr(ct_rate, "input_data"):
            data = ct_rate.input_data

            # Cantera 3.0+ nests parameters under 'rate-constant'.
            # We use .get() to support both flat (older) and nested (newer) dicts.
            rc = data.get("rate-constant", data)

            # Check A (m^3/kmol/s)
            expected_A_kmol = self.plasma.A.value_si * 1000.0

            A_val = rc.get("A")
            # Handle Cantera Quantity objects if present
            if hasattr(A_val, "value"):
                A_val = A_val.value

            assert A_val is not None, f"Could not find 'A' in {data}"
            assert np.isclose(A_val, expected_A_kmol, rtol=1e-4)

            # Check b
            assert np.isclose(rc.get("b"), self.n)

            # Check Ea (J/kmol)
            expected_Ea_g_kmol = self.plasma.Ea_g.value_si * 1000.0
            expected_Ea_e_kmol = self.plasma.Ea_e.value_si * 1000.0

            Ea_g_val = rc.get("Ea-gas", 0.0)
            Ea_e_val = rc.get("Ea-electron", 0.0)

            if hasattr(Ea_g_val, "value"):
                Ea_g_val = Ea_g_val.value
            if hasattr(Ea_e_val, "value"):
                Ea_e_val = Ea_e_val.value

            assert np.isclose(Ea_g_val, expected_Ea_g_kmol, rtol=1e-4)
            assert np.isclose(Ea_e_val, expected_Ea_e_kmol, rtol=1e-4)


class TestElectronCollisionPlasma:
    """
    Tests for the ElectronCollisionPlasma kinetics model.

    We keep the tests fairly high-level so that they only depend on the
    public API of the class (attributes + rate interface), not on the
    detailed sigma(E) -> k(Te) implementation.
    """

    def setup_method(self):
        # Simple toy cross-section on an energy grid. The exact
        # numbers don't matter, only that they are positive and
        # not all equal so that k(Te) varies with Te.
        self.energy_grid = [1.0, 5.0, 10.0]           # arbitrary energies (eV)
        self.sigma_vals = [1.0e-20, 2.0e-20, 5.0e-21] # m^2 (example)

        self.Tmin = 3000.0       # K (electron temperature range)
        self.Tmax = 50000.0      # K
        self.comment = "e- + X -> X* test collision"

        self.plasma = ElectronCollisionPlasma(
            energies=(self.energy_grid, "eV/molecule"),
            sigma=(self.sigma_vals, "m^2"),
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment=self.comment,
        )

    def test_energy_grid(self):
        """Energy grid is stored and has the expected shape / ordering."""
        E = self.plasma.energies
        # Just check basic structure; don't depend on absolute units.
        # We assume Energy-like with value_si being a 1D array.
        assert hasattr(E, "value_si")
        assert len(E.value_si) == len(self.energy_grid)
        # Monotonic ordering should be preserved after unit conversion
        assert np.all(np.diff(E.value_si) > 0.0)

    def test_cross_section(self):
        """Cross section array is stored and matches the input ratios."""
        sigma = self.plasma.sigma
        assert hasattr(sigma, "value_si")
        assert len(sigma.value_si) == len(self.sigma_vals)

        # Ratios are independent of units, so this is robust
        ratio_input = self.sigma_vals[1] / self.sigma_vals[0]
        ratio_stored = sigma.value_si[1] / sigma.value_si[0]
        assert ratio_stored == pytest.approx(ratio_input)

    def test_temperature_min(self):
        assert self.plasma.Tmin is not None
        assert self.plasma.Tmin.value_si == pytest.approx(self.Tmin)

    def test_temperature_max(self):
        assert self.plasma.Tmax is not None
        assert self.plasma.Tmax.value_si == pytest.approx(self.Tmax)

    def test_comment(self):
        assert self.plasma.comment == self.comment

    def test_uses_electron_temperature_flag(self):
        """ElectronCollisionPlasma should advertise Te-dependence."""
        assert getattr(self.plasma, "uses_electron_temperature", False)

    def test_is_temperature_valid(self):
        """Validity check should use the Te bounds (Tmin, Tmax)."""
        assert not self.plasma.is_temperature_valid(self.Tmin - 1.0)
        assert self.plasma.is_temperature_valid(self.Tmin)
        assert self.plasma.is_temperature_valid((self.Tmin + self.Tmax) / 2.0)
        assert self.plasma.is_temperature_valid(self.Tmax)
        assert not self.plasma.is_temperature_valid(self.Tmax + 1.0)

    def test_get_rate_coefficient_alias(self):
        """
        get_rate_coefficient(T) should be a Te-only view and agree
        with the explicit electron-temperature method, if provided.
        """
        Te = 20000.0  # K
        # If the dedicated Te method exists, require equality.
        if hasattr(self.plasma, "get_rate_coefficient_electron_temp"):
            k1 = self.plasma.get_rate_coefficient(Te)
            k2 = self.plasma.get_rate_coefficient_electron_temp(Te)
            assert k1 == pytest.approx(k2)
        else:
            # At minimum, it must be callable and positive
            k = self.plasma.get_rate_coefficient(Te)
            assert k > 0.0

    def test_rate_increases_with_electron_temperature(self):
        """
        For a physically reasonable sigma(E) and EEDF model, the effective
        rate should increase with Te over a broad range.
        """
        Te_low = 10000.0
        Te_mid = 20000.0
        Te_high = 40000.0

        k_low = self.plasma.get_rate_coefficient(Te_low)
        k_mid = self.plasma.get_rate_coefficient(Te_mid)
        k_high = self.plasma.get_rate_coefficient(Te_high)

        assert k_low > 0.0
        assert k_mid > 0.0
        assert k_high > 0.0
        assert k_mid > k_low
        assert k_high > k_mid

    def test_change_rate(self):
        """
        Scaling the rate should uniformly scale k(Te) for any Te.
        This mirrors the behaviour of Arrhenius / TwoTemperaturePlasma.
        """
        Te = 15000.0
        k0 = self.plasma.get_rate_coefficient(Te)

        self.plasma.change_rate(3.0)
        k1 = self.plasma.get_rate_coefficient(Te)
        assert k1 == pytest.approx(3.0 * k0)

    def test_is_identical_to_distinguishes_cross_sections(self):
        """
        Two ElectronCollisionPlasma objects that share a validity window but carry
        different cross-section tables must not compare identical.

        ``is_identical_to`` is declared for this class in ``arrhenius.pxd``. Without
        its own implementation the call falls through to ``KineticsModel``'s, which
        compares only Tmin/Tmax -- so any two objects with the same window compare
        identical no matter what sigma(E) says.
        """
        def build(energies, sigma):
            return ElectronCollisionPlasma(
                energies=(energies, "eV/molecule"),
                sigma=(sigma, "m^2"),
                Tmin=(self.Tmin, "K"),
                Tmax=(self.Tmax, "K"),
                comment=self.comment,
            )

        same = build(self.energy_grid, self.sigma_vals)
        assert self.plasma.is_identical_to(same)

        # A completely different cross-section on the same energy grid.
        other_sigma = build(self.energy_grid, [9.0e-20, 8.0e-20, 7.0e-20])
        assert not self.plasma.is_identical_to(other_sigma)

        # The same cross-section values on a different energy grid.
        other_energies = build([2.0, 6.0, 11.0], self.sigma_vals)
        assert not self.plasma.is_identical_to(other_energies)

        # And a different kinetics type entirely.
        assert not self.plasma.is_identical_to(
            Arrhenius(A=(1.0e10, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kJ/mol"),
                      Tmin=(self.Tmin, "K"), Tmax=(self.Tmax, "K"))
        )

    def test_pickle(self):
        """ElectronCollisionPlasma should round-trip via pickle."""
        blob = pickle.dumps(self.plasma, protocol=-1)
        other = pickle.loads(blob)

        assert isinstance(other, ElectronCollisionPlasma)
        assert other.comment == self.plasma.comment

        # Energy / sigma arrays should come back intact
        assert np.allclose(other.energies.value_si, self.plasma.energies.value_si)
        assert np.allclose(other.sigma.value_si, self.plasma.sigma.value_si)

        # And the Te-only rate should match as well
        Te = 18000.0
        assert other.get_rate_coefficient(Te) == pytest.approx(
            self.plasma.get_rate_coefficient(Te)
        )

    def test_repr(self):
        """repr() should be informative and round-trippable in spirit."""
        s = repr(self.plasma)
        assert "ElectronCollisionPlasma" in s
        assert "comment=" in s
        assert self.comment in s

    def test_cantera_yaml_structure(self):
        """
        Verify that the converted Cantera object mimics the 'electron-collision-plasma' YAML.
        """
        pytest.importorskip("cantera")

        # 1. Convert
        ct_rate = self.plasma.to_cantera_kinetics()

        # 2. Check structure via input_data (Cantera 3.0+)
        if hasattr(ct_rate, "input_data"):
            data = ct_rate.input_data

            # The YAML type for this is usually 'electron-collision-plasma'
            # The fields are 'energy-levels' and 'cross-sections'

            assert "energy-levels" in data
            assert "cross-sections" in data

            energies_out = np.array(data["energy-levels"])
            sigma_out = np.array(data["cross-sections"])

            # 3. Validate Content
            # Cantera input_data usually returns the raw values passed to constructor.
            # We passed: energies in eV, sigma in m^2.

            # Energies (should match our input eV grid)
            # RMG stores J/mol -> we converted to eV for Cantera.
            # Factor: ~96485 J/mol per eV
            conversion_factor = 96485.33212
            expected_eV = self.plasma.energies.value_si / conversion_factor
            assert np.allclose(energies_out, expected_eV, rtol=1e-4)

            # Cross-sections (m^2)
            expected_sigma = self.plasma.sigma.value_si
            assert np.allclose(sigma_out, expected_sigma, rtol=1e-4)


class TestBadnellRRArrhenius:
    """
    Contains unit tests of the :class:`BadnellRRArrhenius` class.
    """

    def setup_method(self):
        # Primary test object uses per-molecule A and includes the C/T2 correction
        self.A_cms_per_molecule = 5.0e-12
        self.B = 0.10
        self.T0 = 350.0
        self.T1 = 2.30e6
        self.C = 0.50
        self.T2 = 5.00e4
        self.Tmin = 300.0
        self.Tmax = 5.0e4
        self.comment = "Na+ + e- -> Na + hv (test)"

        self.rr = BadnellRRArrhenius(
            A=(self.A_cms_per_molecule, "cm^3/(molecule*s)"),
            B=self.B,
            T0=(self.T0, "K"),
            T1=(self.T1, "K"),
            C=self.C,
            T2=(self.T2, "K"),
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment=self.comment,
        )

        # Secondary object: no C/T2, A provided per-mole
        self.A_cms_per_mol = 1.2e-7
        self.rr_molar = BadnellRRArrhenius(
            A=(self.A_cms_per_mol, "cm^3/(mol*s)"),
            B=self.B,
            T0=(self.T0, "K"),
            T1=(self.T1, "K"),
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment="per-mole variant",
        )

    # -------- property tests --------

    def test_a_factor_per_molecule_units(self):
        """
        RateCoefficient converts 'cm^3/(molecule*s)' -> 'm^3/(mol*s)' via (x 1e-6) and (x N_A).
        """
        expected_SI = self.A_cms_per_molecule * 1e-6 * constants.Na  # m^3/(mol*s)
        assert abs(self.rr.A.value_si - expected_SI) <= 1e-12 * max(1.0, expected_SI)

    def test_b_parameter(self):
        assert round(abs(self.rr.B.value_si - self.B), 12) == 0

    def test_t_params(self):
        assert round(abs(self.rr.T0.value_si - self.T0), 12) == 0
        assert round(abs(self.rr.T1.value_si - self.T1), 6) == 0

    def test_optional_cterms_present(self):
        assert self.rr.C is not None
        assert self.rr.T2 is not None
        assert round(abs(self.rr.C.value_si - self.C), 12) == 0
        assert round(abs(self.rr.T2.value_si - self.T2), 6) == 0

    def test_temperature_min_max(self):
        assert round(abs(self.rr.Tmin.value_si - self.Tmin), 6) == 0
        assert round(abs(self.rr.Tmax.value_si - self.Tmax), 6) == 0

    def test_comment(self):
        assert self.rr.comment == self.comment

    def test_electron_flags_on_direct_parameter_path(self):
        """Direct-parameter construction advertises Te- and ne-dependence."""
        assert self.rr.uses_electron_temperature is True
        assert self.rr.uses_electron_density is True

    def test_is_temperature_valid(self):
        Tdata = np.array([200, 400, 6000, 12000, 30000, 60000])
        validdata = np.array([False, True, True, True, True, False], dtype=bool)
        for T, valid in zip(Tdata, validdata):
            assert self.rr.is_temperature_valid(T) == valid

    # -------- evaluator tests --------

    def test_get_rate_coefficient_matches_formula_per_molecule_input(self):
        """
        Compare get_rate_coefficient (returns SI m^3/(mol*s)) to an independently computed expectation
        from the Badnell formula using per-molecule A, then converting with N_A and cm->m.
        """
        Tlist = np.array([800, 2000, 5000, 11600, 20000, 40000], dtype=float)
        for T in Tlist:
            alpha_cm = _alpha_rr_cm3_per_molecule_s(
                T, self.A_cms_per_molecule, self.B, self.T0, self.T1, self.C, self.T2
            )
            kexp_SI = alpha_cm * constants.Na * 1e-6  # m^3/(mol*s)
            kact = self.rr.get_rate_coefficient(T)
            assert abs(kexp_SI - kact) <= 1e-12 * max(1.0, kexp_SI)

    def test_get_rate_coefficient_matches_formula_per_mole_input(self):
        """
        When A is per-mole, the returned k should *not* multiply by N_A.
        """
        Tlist = np.array([800, 2000, 5000, 11600, 20000, 40000], dtype=float)
        for T in Tlist:
            # Expected: alpha_molar = (A_molar / denom). Convert cm^3->m^3 at end.
            B = self.B
            s0 = np.sqrt(T / self.T0)
            s1 = np.sqrt(T / self.T1)
            Bstar = B  # no C/T2 in this object
            denom = s0 * (1.0 + s0) ** (1.0 - Bstar) * (1.0 + s1) ** (1.0 + Bstar)
            alpha_molar_cm = self.A_cms_per_mol / denom  # cm^3/(mol*s)
            kexp_SI = alpha_molar_cm * 1e-6  # m^3/(mol*s)
            kact = self.rr_molar.get_rate_coefficient(T)
            assert abs(kexp_SI - kact) <= 1e-12 * max(1.0, kexp_SI)

    def test_change_rate(self):
        Tlist = np.array([1000, 5000, 15000, 30000], dtype=float)
        k0 = np.array([self.rr.get_rate_coefficient(T) for T in Tlist])
        self.rr.change_rate(2.0)
        k1 = np.array([self.rr.get_rate_coefficient(T) for T in Tlist])
        assert np.allclose(k1, 2.0 * k0, rtol=1e-12, atol=0.0)

    def test_repr(self):
        namespace = {}
        exec("rr_obj = {0!r}".format(self.rr), globals(), namespace)
        assert "rr_obj" in namespace
        rr2 = namespace["rr_obj"]
        # spot-check a few fields
        assert self.rr.A.units == rr2.A.units
        assert round(abs(self.rr.B.value - rr2.B.value), 12) == 0
        assert round(abs(self.rr.T0.value - rr2.T0.value), 12) == 0
        assert round(abs(self.rr.T1.value - rr2.T1.value), 6) == 0
        assert self.rr.comment == rr2.comment

    # -------- compatibility test against Arrhenius (power-law fit) --------

    def test_match_power_law_fit_over_narrow_window(self):
        """
        Over a narrow Te window, Badnell can be approximated by k = Apl * Te^n (Ea=0).
        Fit Arrhenius to a small set of points and verify reconstruction error is small.
        """
        # sample Te window
        Tdata = np.array([6000, 8000, 10000, 12000, 15000], dtype=float)
        kdata = np.array([self.rr.get_rate_coefficient(T) for T in Tdata])  # SI m^3/(mol*s)

        # Fit plain Arrhenius with Ea=0 (two-parameter fit): emulate by setting three_params=False
        arr = Arrhenius().fit_to_data(Tdata, kdata, kunits="m^3/(mol*s)", T0=1.0, three_params=True)

        # Reproduce rates within small error over the window
        for T, kexp in zip(Tdata, kdata):
            kfit = arr.get_rate_coefficient(T)
            assert abs(kfit - kexp) <= 5e-3 * kexp  # ~0.5% is plenty strict here

    def _expected_rate_SI_per_mol(self, entry, T):
        """
        Compute expected SI m^3/(mol*s) using Badnell (2006) form with:
        A in cm^3/(molecule*s) -> m^3/(mol*s) via 1e-6 * N_A
        """
        A = float(entry["A"])
        B = float(entry["B"])
        T0 = float(entry["T0"])
        T1 = float(entry["T1"])
        C = float(entry["C"]) if "C" in entry and entry["C"] is not None else None
        T2 = float(entry["T2"]) if "T2" in entry and entry["T2"] is not None else None

        N_A = 6.02214076e23
        A_SI = A * 1e-6 * N_A

        s0 = math.sqrt(T / T0)
        s1 = math.sqrt(T / T1)
        Bstar = B + (C * math.exp(-T2 / T) if (C is not None and T2 is not None) else 0.0)
        denom = s0 * (1.0 + s0) ** (1.0 - Bstar) * (1.0 + s1) ** (1.0 + Bstar)
        return A_SI / denom

    def test_electron_flags_on_yaml_path(self):
        """
        Building through the (Z, N) table path must still advertise Te- and
        ne-dependence.

        That path returns as soon as ``populate_from_yaml`` comes back, and
        ``populate_from_yaml`` sets neither flag. Cython zero-initialises a
        ``cdef public bint``, so a skipped assignment is observable as ``False``
        rather than as an AttributeError -- which is exactly how a table-built
        object came to claim it needs neither the electron temperature nor the
        electron density.
        """
        path = _require_plasma_yaml("badnell.yaml")

        m = BadnellRRArrhenius(Z=1, N=0, yaml_path_or_obj=path)
        assert m.uses_electron_temperature is True
        assert m.uses_electron_density is True

        # And for a row that carries the optional C/T2 terms.
        m2 = BadnellRRArrhenius(Z=2, N=1, yaml_path_or_obj=path)
        assert m2.uses_electron_temperature is True
        assert m2.uses_electron_density is True

    def test_init_from_yaml_Z1_N0_basic_loads_and_sets_window_and_rate(self):
        """
        Uses well-known hydrogen entry Z=1, N=0 from the real YAML.
        Asserts parameters, default Tmin/Tmax window, comment, and rate at a T.
        """
        m = BadnellRRArrhenius(Z=1, N=0, yaml_path_or_obj=_require_plasma_yaml("badnell.yaml"))

        # Parameters from file (Hydrogen Z=1,N=0)
        assert pytest.approx(m.B.value_si, rel=0, abs=1e-12) == 0.7472
        assert pytest.approx(m.T0.value_si, rel=0, abs=1e-12) == 2.965e0
        assert pytest.approx(m.T1.value_si, rel=0, abs=1e-3) == 7.001e5
        assert m.C is None
        assert m.T2 is None

        # Default validity window: z = Z - N = 1 -> [10, 1e7] K
        assert pytest.approx(m.Tmin.value_si, rel=0, abs=1e-12) == 10.0
        assert pytest.approx(m.Tmax.value_si, rel=0, abs=1e-3) == 1.0e7

        # Comment includes the tag
        assert "Z=1" in m.comment and "N=0" in m.comment

        # Rate check (numbers from the paper entry)
        T = 1.0e4
        entry = {"A": 8.318e-11, "B": 0.7472, "T0": 2.965e0, "T1": 7.001e5}
        k_expected = self._expected_rate_SI_per_mol(entry, T)
        assert np.isclose(m.get_rate_coefficient(T), k_expected)

    def test_init_from_yaml_Z2_N1_with_C_T2_and_rate(self):
        """
        Helium-like case that includes C and T2; verifies optional params and rate.
        """
        m = BadnellRRArrhenius(Z=2, N=1, yaml_path_or_obj=_require_plasma_yaml("badnell.yaml"))

        # Params present
        assert pytest.approx(m.B.value_si, abs=1e-12) == 0.6988
        assert m.C is not None and pytest.approx(m.C.value_si, abs=1e-12) == 8.29e-2
        assert m.T2 is not None and pytest.approx(m.T2.value_si, abs=1e-3) == 1.682e5

        # Rate at a temperature that exercises exp(-T2/T)
        T = 2.0e5
        entry = {"A": 5.235e-11, "B": 0.6988, "T0": 7.301e0, "T1": 4.475e6, "C": 8.29e-2, "T2": 1.682e5}
        k_expected = self._expected_rate_SI_per_mol(entry, T)
        assert np.isclose(m.get_rate_coefficient(T), k_expected)

    def test_init_from_yaml_comment_override_is_appended(self):
        m = BadnellRRArrhenius(Z=1, N=0, yaml_path_or_obj=_require_plasma_yaml("badnell.yaml"),
                               comment="custom note")
        assert "Badnell (2006)" in m.comment
        assert "Z=1" in m.comment and "N=0" in m.comment
        assert "custom note" in m.comment

    def test_init_from_yaml_default_path_resolution_works(self):
        """
        Ensure that omitting yaml_path_or_obj uses the default path via settings.

        Skipped when the configured database does not carry badnell.yaml -- the
        table ships with RMG-database, and `settings['database.directory']` may
        point at a checkout that predates it.
        """
        default_path = os.path.join(settings["database.directory"], "kinetics", "badnell.yaml")
        if not os.path.isfile(default_path):
            pytest.skip("configured database has no kinetics/badnell.yaml: {0}".format(default_path))

        m = BadnellRRArrhenius(Z=1, N=0)  # no yaml_path_or_obj
        assert pytest.approx(m.B.value_si, abs=1e-12) == 0.7472

    def test_init_from_yaml_blocks_Z_gt36_by_default(self):
        """
        __init__ calls populate_from_yaml without allow_Z_gt36=True, so Z>36 should raise.
        """
        with pytest.raises(ValueError):
            BadnellRRArrhenius(Z=80, N=6, yaml_path_or_obj=_require_plasma_yaml("badnell.yaml"))

    def test_init_from_yaml_early_return_overrides_manual_values(self):
        """
        If Z/N are provided, __init__ should populate from YAML and ignore manual A/B/T*.
        """
        m = BadnellRRArrhenius(
            A=(9.99e-10, "cm^3/(molecule*s)"),
            B=9.99,
            T0=(9.99, "K"),
            T1=(9.99, "K"),
            Z=1, N=0,
            yaml_path_or_obj=_require_plasma_yaml("badnell.yaml"),
        )
        # Values should be from YAML, not the manual ones above
        assert pytest.approx(m.B.value_si, abs=1e-12) == 0.7472
        assert pytest.approx(m.T0.value_si, abs=1e-12) == 2.965e0
        assert pytest.approx(m.T1.value_si, abs=1e-3) == 7.001e5

    def test_init_from_yaml_validates_integer_ZN(self):
        path = _require_plasma_yaml("badnell.yaml")
        with pytest.raises(TypeError):
            BadnellRRArrhenius(Z="not-int", N=0, yaml_path_or_obj=path)
        with pytest.raises(TypeError):
            BadnellRRArrhenius(Z=1, N="nope", yaml_path_or_obj=path)

    def test_to_two_temp_plasma_structure_and_flags(self):
        """
        BadnellRRArrhenius.to_two_temp_plasma should return a TwoTemperaturePlasma
        that uses Te and ne and inherits a sensible Tmin/Tmax window.
        """
        plasma = self.rr.to_two_temp_plasma()

        # Type and flags
        assert isinstance(plasma, TwoTemperaturePlasma)
        assert plasma.uses_electron_temperature is True

        # Tmin/Tmax should be based on the Badnell fit window
        assert pytest.approx(plasma.Tmin.value_si, rel=0, abs=1e-6) == self.rr.Tmin.value_si
        assert pytest.approx(plasma.Tmax.value_si, rel=0, abs=1e-3) == self.rr.Tmax.value_si

    def test_to_two_temp_plasma_matches_badnell_along_TeqTe(self):
        """
        Along T = Te, the TwoTemperaturePlasma mapping should reproduce the
        Badnell rate within a modest tolerance over a central Te window.
        """
        plasma = self.rr.to_two_temp_plasma()

        # Pick temperatures comfortably inside the validity window
        Tvals = np.array([2000.0, 8000.0, 20000.0], dtype=float)
        for T in Tvals:
            k_badnell = self.rr.get_rate_coefficient(T)
            # TwoTemperaturePlasma expects (T_gas, T_e)
            k_plasma = plasma.get_rate_coefficient(T, T)

            # The Badnell function (especially with C/T2 terms) has curvature
            # that a single Arrhenius form cannot fit to better than a few percent
            # over the wide range (300-50,000 K) used in this test object;
            # ~9% deviation is expected and acceptable for this approximation.
            assert np.isclose(k_plasma, k_badnell, rtol=0.15)

    def test_to_two_temp_plasma_Ea_partitioning(self):
        """
        For a k(Te)-only Badnell rate, the TwoTemperaturePlasma mapping should
        put all activation into the electron channel (Ea_g ~ 0, Ea_e ~ Ea_fit).
        """
        plasma = self.rr.to_two_temp_plasma()
        arr = self.rr.to_arrhenius()

        Ea_fit = arr.Ea.value_si

        # Gas activation channel "off"
        assert pytest.approx(plasma.Ea_g.value_si, rel=0, abs=1e-6) == 0.0
        # Electron activation follows the Arrhenius fit in Te
        assert pytest.approx(plasma.Ea_e.value_si, rel=1e-6, abs=0.0) == Ea_fit

    def test_to_cantera_kinetics_uses_two_temp_rate(self):
        """
        to_cantera_kinetics should go through the canonical TwoTemperaturePlasma
        mapping and return a Cantera TwoTempPlasmaRate (or equivalent type).
        """
        ct = pytest.importorskip("cantera")

        plasma = self.rr.to_two_temp_plasma()

        rate_from_plasma = plasma.to_cantera_kinetics()
        rate_from_badnell = self.rr.to_cantera_kinetics()

        # Both paths must produce the same underlying rate type
        assert isinstance(rate_from_badnell, type(rate_from_plasma))
        # Sanity check on the name for easier debugging
        assert "TwoTempPlasma" in rate_from_badnell.__class__.__name__
        assert isinstance(rate_from_badnell, ct.ReactionRate)


class TestVoronovEIArrhenius:
    """
    Contains unit tests of the :class:`VoronovEIArrhenius` class.
    """

    def setup_method(self):
        # Primary test object uses per-molecule A (constructed directly; no YAML).
        self.A_cms_per_molecule = 3.0e-8
        self.P = 0.25
        self.X = 1.80
        self.K = 0.42
        self.dE_eV = 5.14   # ~Na I first ionization (arbitrary for tests)
        self.Z = 11
        self.N = 11  # electrons before ionization (neutral Na)
        self.Tmin = 1.0e4
        self.Tmax = 2.0e5
        self.comment = "Na (I) -> Na+ + e- (test, per-molecule A)"

        self.ei = VoronovEIArrhenius(
            A=(self.A_cms_per_molecule, "cm^3/(molecule*s)"),
            P=self.P,
            X=self.X,
            K=self.K,
            dE=self.dE_eV,
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment=self.comment,
        )

        # Secondary object: A provided per-mole (constructed directly; no YAML).
        self.A_cms_per_mol = 1.0e-6
        self.ei_molar = VoronovEIArrhenius(
            A=(self.A_cms_per_mol, "cm^3/(mol*s)"),
            P=self.P,
            X=self.X,
            K=self.K,
            dE=self.dE_eV,
            Tmin=(self.Tmin, "K"),
            Tmax=(self.Tmax, "K"),
            comment="per-mole variant",
        )

        # YAML fixtures as dictionaries (schema-compatible)
        # Per-molecule A entry (Z=11, N=11)
        self.yaml_per_molecule = {
            "units": {"A": "cm^3/(molecule*s)", "dE": "eV", "Tmin": "eV", "Tmax": "eV"},
            "coefficients": [
                {
                    "Z": 11,
                    "entries": [
                        {
                            "N": 11,
                            "A": 3.0e-8,
                            "P": 0.25,
                            "X": 1.80,
                            "K": 0.42,
                            "dE": 5.14,
                            "Tmin": 0.010,   # eV
                            "Tmax": 12.0,    # eV
                            "species": "Na I",
                            "comment": "Na (I) -> Na+ + e- (YAML test)"
                        }
                    ]
                }
            ],
        }

        # Per-mole A entry (Z=11, N=10)
        self.yaml_per_mole = {
            "units": {"A": "cm^3/(mol*s)", "dE": "eV", "Tmin": "eV", "Tmax": "eV"},
            "coefficients": [
                {
                    "Z": 11,
                    "entries": [
                        {
                            "N": 10,
                            "A": 1.0e-6,
                            "P": 0.25,
                            "X": 1.80,
                            "K": 0.42,
                            "dE": 47.3,
                            "Tmin": 0.020,  # eV
                            "Tmax": 10.0,   # eV
                            "species": "Na II",
                            "comment": "Na (II) -> Na2+ + e- (YAML test)"
                        }
                    ]
                }
            ],
        }

    # -------- property tests --------

    def test_a_factor_per_molecule_units(self):
        """
        RateCoefficient converts 'cm^3/(molecule*s)' -> SI m^3/(mol*s) via (x 1e-6) and (x N_A).
        """
        expected_SI = self.A_cms_per_molecule * 1e-6 * constants.Na  # m^3/(mol*s)
        assert np.isclose(self.ei.A.value_si, expected_SI, rtol=0, atol=1e-12 * max(1.0, expected_SI))

    def test_dimensionless_params(self):
        assert round(abs(self.ei.P.value_si - self.P), 12) == 0
        assert round(abs(self.ei.X.value_si - self.X), 12) == 0
        assert round(abs(self.ei.K.value_si - self.K), 12) == 0

    def test_threshold_energy_stored_as_double(self):
        # dE_eV is stored as a plain double; accessor returns the numeric eV.
        assert round(abs(self.ei.dE_eV - self.dE_eV), 12) == 0

    def test_temperature_min_max(self):
        assert np.isclose(self.ei.Tmin.value_si, self.Tmin, rtol=0, atol=1e-9)
        assert np.isclose(self.ei.Tmax.value_si, self.Tmax, rtol=0, atol=1e-9)

    def test_comment_and_stage(self):
        # When constructed directly, comment should be preserved verbatim.
        assert self.comment in self.ei.comment

    def test_electron_flags_on_direct_parameter_path(self):
        """Direct-parameter construction advertises Te- and ne-dependence."""
        assert self.ei.uses_electron_temperature is True
        assert self.ei.uses_electron_density is True

    def test_is_temperature_valid(self):
        Tdata = np.array([200, 400, 1.0e4, 5.0e4, 2.1e5])
        validdata = np.array([False, False, True, True, False], dtype=bool)
        for T, valid in zip(Tdata, validdata):
            assert self.ei.is_temperature_valid(T) == valid

    # -------- evaluator tests --------

    def test_get_rate_coefficient_matches_formula_per_molecule_input(self):
        """
        Compare get_rate_coefficient (SI m^3/(mol*s)) with independent Voronov formula,
        using per-molecule A then converting with N_A and cm->m.
        """
        Tlist = np.array([800, 2000, 5000, 10000, 20000, 80000], dtype=float)
        for T in Tlist:
            alpha_cm = _alpha_ei_cm3_per_molecule_s(
                T, self.A_cms_per_molecule, self.P, self.X, self.K, self.dE_eV
            )  # cm^3/(molecule*s)
            kexp_SI = alpha_cm * constants.Na * 1e-6  # m^3/(mol*s)
            kact = self.ei.get_rate_coefficient(T)
            assert np.isclose(kact, kexp_SI, rtol=5e-4, atol=1e-3)

    def test_get_rate_coefficient_matches_formula_per_mole_input(self):
        """
        When A is per-mole, the returned k should *not* multiply by N_A.
        """
        Tlist = np.array([800, 2000, 5000, 10000, 20000, 80000], dtype=float)
        for T in Tlist:
            # Same functional form, but A already per mole -> only cm^3->m^3 at the end.
            alpha_molar_cm = _alpha_ei_cm3_per_molecule_s(
                T, self.A_cms_per_mol, self.P, self.X, self.K, self.dE_eV
            )  # cm^3/(mol*s)
            kexp_SI = alpha_molar_cm * 1e-6  # m^3/(mol*s)
            kact = self.ei_molar.get_rate_coefficient(T)
            assert np.isclose(kact, kexp_SI, rtol=5e-4, atol=1e-3)

    def test_change_rate(self):
        Tlist = np.array([1000, 5000, 15000, 30000], dtype=float)
        k0 = np.array([self.ei.get_rate_coefficient(T) for T in Tlist])
        self.ei.change_rate(2.5)
        k1 = np.array([self.ei.get_rate_coefficient(T) for T in Tlist])
        assert np.allclose(k1, 2.5 * k0, rtol=1e-12, atol=0.0)

    def test_repr_roundtrip(self):
        # The repr should contain key fields; we don't require eval roundtrip here.
        r = repr(self.ei)
        assert "VoronovEIArrhenius(" in r
        assert "dE=" in r and "eV" in r
        assert "comment" in r

    # -------- threshold-energy handling --------

    def test_missing_dE_raises_rather_than_fabricating_a_threshold(self):
        """
        dE is the ionisation threshold and the rate is exponential in it, so there
        is no defensible default: absent a value the direct-parameter path must
        fail loudly rather than quietly substitute one.
        """
        with pytest.raises(ValueError):
            VoronovEIArrhenius(
                A=(self.A_cms_per_molecule, "cm^3/(molecule*s)"),
                P=self.P,
                X=self.X,
                K=self.K,
            )

        # A bare constructor is the same case with everything omitted.
        with pytest.raises(ValueError):
            VoronovEIArrhenius()

        # dE=0 is a real (if degenerate) value, not an absent one, and must be kept.
        assert VoronovEIArrhenius(dE=0.0).dE_eV == 0.0

    # -------- YAML-driven tests (selection/queries & units) --------

    def test_electron_flags_on_yaml_path(self):
        """
        Building through the (Z, N) table path must still advertise Te- and
        ne-dependence -- see the sibling test on BadnellRRArrhenius for why a
        skipped assignment is observable as ``False``.
        """
        path = _require_plasma_yaml("voronov.yaml")

        m = VoronovEIArrhenius(Z=1, N=1, yaml_path_or_obj=path)
        assert m.uses_electron_temperature is True
        assert m.uses_electron_density is True

    def test_init_from_yaml_Z1_N1_loads_hydrogen_row(self):
        """
        The hydrogen row (Z=1, N=1) of the real table: H I -> H II at 13.6 eV.
        """
        m = VoronovEIArrhenius(Z=1, N=1, yaml_path_or_obj=_require_plasma_yaml("voronov.yaml"))

        assert pytest.approx(m.dE_eV, abs=1e-12) == 13.6
        assert pytest.approx(m.P.value_si, abs=1e-12) == 0.0
        assert pytest.approx(m.X.value_si, abs=1e-12) == 2.32e-1
        assert pytest.approx(m.K.value_si, abs=1e-12) == 0.39
        assert "Z=1" in m.comment and "N=1" in m.comment

        # A is tabulated per molecule; RMG stores SI molar.
        expected_SI_A = 2.910e-08 * 1e-6 * constants.Na
        assert np.isclose(m.A.value_si, expected_SI_A, rtol=1e-12)

        # Rate agrees with an independent evaluation of the Voronov formula.
        for T in (2.0e4, 1.0e5):
            alpha_cm = _alpha_ei_cm3_per_molecule_s(T, 2.910e-08, 0.0, 2.32e-1, 0.39, 13.6)
            assert np.isclose(m.get_rate_coefficient(T), alpha_cm * constants.Na * 1e-6, rtol=5e-4)

    def test_populate_from_yaml_selects_correct_entry_by_ZN_per_molecule(self):
        # dE is overwritten by populate_from_yaml; it is supplied here only because
        # the direct-parameter path refuses to invent a threshold.
        obj = VoronovEIArrhenius(dE=1.0)

        # YAML lists Tmin/Tmax in eV; compute the Kelvin values and pass with units
        kB = 8.617333262e-5  # eV/K
        expected_Tmin_K = 0.010 / kB
        expected_Tmax_K = 12.0 / kB

        # Pass explicit K-units to satisfy Quantity setter
        obj.populate_from_yaml(
            self.yaml_per_molecule,
            Z=11,
            N=11,
            Tmin=(expected_Tmin_K, "K"),
            Tmax=(expected_Tmax_K, "K"),
        )

        # A: cm^3/(molecule*s) -> m^3/(mol*s)
        expected_SI_A = 3.0e-8 * 1e-6 * constants.Na  # m^3/(mol*s)
        assert np.isclose(obj.A.value_si, expected_SI_A, rtol=0, atol=1e-12 * max(1.0, expected_SI_A))
        assert round(abs(obj.P.value_si - 0.25), 12) == 0
        assert round(abs(obj.X.value_si - 1.80), 12) == 0
        assert round(abs(obj.K.value_si - 0.42), 12) == 0
        assert round(abs(obj.dE_eV - 5.14), 12) == 0

        # Comment currently built as "Voronov ... Z=..., N=..."; ensure it reflects selection
        assert "Z=11" in obj.comment and "N=11" in obj.comment

        # Tmin/Tmax were provided in K; verify they are stored as Kelvin and match
        assert np.isclose(obj.Tmin.value_si, expected_Tmin_K, rtol=0, atol=1e-8 * expected_Tmin_K)
        assert np.isclose(obj.Tmax.value_si, expected_Tmax_K, rtol=0, atol=1e-8 * expected_Tmax_K)

    def test_populate_from_yaml_selects_correct_entry_by_ZN_per_mole(self):
        obj = VoronovEIArrhenius(dE=1.0)  # placeholder; overwritten below

        # YAML lists Tmin/Tmax in eV; convert to K and pass with units
        kB = 8.617333262e-5  # eV/K
        expected_Tmin_K = 0.020 / kB
        expected_Tmax_K = 10.0 / kB

        obj.populate_from_yaml(
            self.yaml_per_mole,
            Z=11,
            N=10,
            Tmin=(expected_Tmin_K, "K"),
            Tmax=(expected_Tmax_K, "K"),
        )

        # A per mole: only cm^3->m^3
        expected_SI_A = 1.0e-6 * 1e-6  # m^3/(mol*s)
        assert np.isclose(obj.A.value_si, expected_SI_A, rtol=0, atol=1e-12 * max(1.0, expected_SI_A))
        assert round(abs(obj.dE_eV - 47.3), 12) == 0

        # Comment contains selected stage
        assert "Z=11" in obj.comment and "N=10" in obj.comment

        # Tmin/Tmax were provided in K; verify stored values
        assert np.isclose(obj.Tmin.value_si, expected_Tmin_K, rtol=0, atol=1e-8 * expected_Tmin_K)
        assert np.isclose(obj.Tmax.value_si, expected_Tmax_K, rtol=0, atol=1e-8 * expected_Tmax_K)

    def test_populate_from_yaml_rate_matches_formula(self):
        """
        After populate_from_yaml, the evaluator should match a direct Voronov computation
        for the selected (Z,N) tuple across several temperatures.
        """
        obj = VoronovEIArrhenius(dE=1.0)  # placeholder; overwritten below

        # YAML lists Tmin/Tmax in eV; convert to K and pass with units
        kB = 8.617333262e-5  # eV/K
        expected_Tmin_K = 0.010 / kB
        expected_Tmax_K = 12.0 / kB

        obj.populate_from_yaml(
            self.yaml_per_molecule,
            Z=11,
            N=11,  # per-molecule case
            Tmin=(expected_Tmin_K, "K"),
            Tmax=(expected_Tmax_K, "K"),
        )

        A = 3.0e-8
        P = 0.25
        X = 1.80
        K = 0.42
        dE = 5.14
        Tlist = np.array([500, 2000, 10000, 40000, 120000], dtype=float)

        for T in Tlist:
            alpha_cm = _alpha_ei_cm3_per_molecule_s(T, A, P, X, K, dE)  # cm^3/(molecule*s)
            kexp_SI = alpha_cm * constants.Na * 1e-6  # m^3/(mol*s)
            kact = obj.get_rate_coefficient(T)
            assert np.isclose(kact, kexp_SI, rtol=5e-4, atol=1e-3)

    # -------- 2T-Plasma & Cantera conversion tests --------

    def test_to_two_temp_plasma_structure_and_flags(self):
        """
        Ensure to_two_temp_plasma returns a TwoTemperaturePlasma object
        with the correct temperature flags set.
        """
        plasma = self.ei.to_two_temp_plasma()

        assert isinstance(plasma, TwoTemperaturePlasma)
        assert plasma.uses_electron_temperature is True
        # Check window inheritance
        assert np.isclose(plasma.Tmin.value_si, self.Tmin, rtol=0, atol=1e-9)
        assert np.isclose(plasma.Tmax.value_si, self.Tmax, rtol=0, atol=1e-9)

    def test_to_two_temp_plasma_Ea_partitioning(self):
        """
        Verify that Voronov (electron-driven) kinetics are mapped such that:
        1. Ea_g is forced to 0.0 (gas T is a spectator).
        2. Ea_e contains the activation energy from the Arrhenius fit.
        """
        plasma = self.ei.to_two_temp_plasma()
        arr = self.ei.to_arrhenius()

        # 1. Gas activation must be zero
        assert pytest.approx(plasma.Ea_g.value_si, abs=1e-9) == 0.0

        # 2. Electron activation matches the fit
        assert pytest.approx(plasma.Ea_e.value_si, rel=1e-6) == arr.Ea.value_si

        # Sanity check: The fitted Ea should be roughly close to the threshold dE
        # dE = 5.14 eV ~ 495.9 kJ/mol.
        expected_Ea_J = self.dE_eV * 96485.332  # eV to J/mol
        assert np.isclose(plasma.Ea_e.value_si, expected_Ea_J, rtol=0.25)

    def test_to_cantera_kinetics_generates_correct_object(self):
        """
        Ensure to_cantera_kinetics produces a valid Cantera Rate object.
        """
        ct = pytest.importorskip("cantera")

        # This implicitly calls to_two_temp_plasma -> to_cantera_kinetics
        ct_rate = self.ei.to_cantera_kinetics()

        # Check type
        assert isinstance(ct_rate, ct.TwoTempPlasmaRate)


################################################################################
#
#   Cross-cutting checks on the plasma kinetics classes.
#
#   These deliberately test *properties* rather than individual call sites: each
#   of the three defects they were written for had a sibling on a code path the
#   original per-method tests never reached.
#
################################################################################


def _arrhenius_pyx_path():
    """
    Locate the Cython source of :mod:`rmgpy.kinetics.arrhenius`.

    The compiled extension sits next to its ``.pyx`` in an in-place build, which is
    how this repo is built. Returns None if the source is not on disk (e.g. an
    installed-wheel layout).
    """
    path = os.path.join(os.path.dirname(os.path.abspath(rmgpy.__file__)), "kinetics", "arrhenius.pyx")
    return path if os.path.isfile(path) else None


def _numeric_literal_arguments_of_calls_to(source, callee_names, line_range=None):
    """
    Find every call in `source` whose callee is one of `callee_names` and report the
    numeric literals passed directly as its arguments.

    Returns a list of ``(line, callee, [literals])`` for the call sites that pass at
    least one numeric literal at the top level of the argument list. Nested calls are
    not descended into: ``(1.0, "K")`` counts, ``foo(bar(2))`` reports only ``foo``'s
    own arguments.

    Tokenising rather than pattern-matching is what makes this reliable: a call
    spelled inside a string -- ``__repr__`` builds one -- never produces a NAME
    token, so it cannot be mistaken for a construction site.
    """
    import tokenize as _tokenize
    import io as _io

    tokens = list(_tokenize.generate_tokens(_io.StringIO(source).readline))
    hits = []

    for index, tok in enumerate(tokens):
        if tok.type != _tokenize.NAME or tok.string not in callee_names:
            continue
        if index + 1 >= len(tokens):
            continue
        opener = tokens[index + 1]
        if opener.type != _tokenize.OP or opener.string != "(":
            continue
        if line_range is not None and not (line_range[0] <= tok.start[0] <= line_range[1]):
            continue

        depth = 0
        literals = []
        for follower in tokens[index + 1:]:
            if follower.type == _tokenize.OP and follower.string in "([{":
                depth += 1
            elif follower.type == _tokenize.OP and follower.string in ")]}":
                depth -= 1
                if depth == 0:
                    break
            elif follower.type == _tokenize.NUMBER and depth in (1, 2):
                # depth 1 = a bare argument; depth 2 = inside a (value, "units") tuple.
                literals.append(follower.string)
        if literals:
            hits.append((tok.start[0], tok.string, literals))

    return hits


def _class_body_line_range(source, class_name):
    """Return the (first, last) 1-based line numbers of ``cdef class <class_name>``."""
    lines = source.splitlines()
    start = None
    for number, line in enumerate(lines, start=1):
        if line.startswith("cdef class "):
            if start is not None:
                return start, number - 1
            if line.startswith("cdef class {0}".format(class_name)):
                start = number
    assert start is not None, "class {0} not found in source".format(class_name)
    return start, len(lines)


class TestVoronovFactoryInventory:
    """
    Every way a :class:`VoronovEIArrhenius` can come into existence, and the rule
    each of them has to obey.

    ``dE`` is the ionisation threshold and the Voronov rate goes as
    ``exp(-dE / Te_eV)``, so a substituted threshold rescales every rate the object
    will ever produce, silently. There is no defensible default. The previous round
    of work enforced that on the direct-parameter constructor and left the sibling
    factory inventing one, which is why this reads as an inventory rather than as a
    test of one call.
    """

    # The complete set of entry points that can hand back a VoronovEIArrhenius.
    # Keep this list and the assertions below in step: `test_..._inventory_is_complete`
    # fails if a factory is added to the class without being accounted for here.
    DECLARED_FACTORIES = frozenset(
        {
            "__init__ (explicit parameters)",
            "__init__ (Z/N table lookup)",
            "from_yaml",
            "__reduce__ (pickle)",
        }
    )

    # Classmethods/staticmethods on the class that construct and return an instance.
    DECLARED_ALTERNATE_CONSTRUCTORS = frozenset({"from_yaml"})

    def setup_method(self):
        # A stand-in table whose threshold is nothing like any plausible fabricated
        # constant, so a substituted value would be unmistakable.
        self.table = {
            "units": {"A": "cm^3/(molecule*s)", "dE": "eV", "Tmin": "eV", "Tmax": "eV"},
            "coefficients": [
                {
                    "Z": 11,
                    "entries": [
                        {
                            "N": 11,
                            "A": 3.0e-8,
                            "P": 0.25,
                            "X": 1.80,
                            "K": 0.42,
                            "dE": 5.139076,
                            "Tmin": 0.010,
                            "Tmax": 12.0,
                        }
                    ],
                }
            ],
        }
        self.table_dE = 5.139076

    def test_voronov_factory_inventory_is_complete(self):
        """
        The enumeration above has to match the class. Anything callable on
        VoronovEIArrhenius that is bound to the class rather than an instance is an
        alternate constructor until proven otherwise, so a new one shows up here
        first and the rest of this class then has to grow a case for it.
        """
        alternate = set()
        for name in dir(VoronovEIArrhenius):
            if name.startswith("__"):
                continue
            attribute = inspect.getattr_static(VoronovEIArrhenius, name, None)
            if isinstance(attribute, (classmethod, staticmethod)):
                alternate.add(name)

        assert alternate == set(self.DECLARED_ALTERNATE_CONSTRUCTORS), (
            "the set of alternate constructors on VoronovEIArrhenius has changed; "
            "add the new one to DECLARED_FACTORIES and give it a threshold check"
        )
        # The declared inventory is those, plus the two __init__ branches and pickle.
        assert self.DECLARED_ALTERNATE_CONSTRUCTORS <= self.DECLARED_FACTORIES
        assert len(self.DECLARED_FACTORIES) == len(self.DECLARED_ALTERNATE_CONSTRUCTORS) + 3

    def test_voronov_factory_no_construction_site_invents_a_parameter(self):
        """
        No construction of a VoronovEIArrhenius anywhere in the module may pass a
        numeric literal.

        This is the check that catches the defect class rather than the instance.
        A fabricated threshold that is overwritten a line later by the table is
        invisible to every behavioural test -- the object never escapes carrying it --
        yet it is a physical constant invented in shipped source, one failed lookup
        away from being the value a mechanism runs on. The only place it is visible
        is where it is written, so that is where it gets asserted about.
        """
        path = _arrhenius_pyx_path()
        if path is None:
            pytest.skip("arrhenius.pyx not available next to the compiled module")
        with open(path) as source_file:
            source = source_file.read()

        line_range = _class_body_line_range(source, "VoronovEIArrhenius")
        hits = _numeric_literal_arguments_of_calls_to(
            source, {"VoronovEIArrhenius", "cls"}, line_range=line_range
        )

        assert hits == [], (
            "VoronovEIArrhenius is constructed with fabricated parameters at: "
            + "; ".join(
                "{0}:{1} {2}({3})".format(path, line, callee, ", ".join(literals))
                for line, callee, literals in hits
            )
        )

    def test_voronov_factory_direct_parameters_demands_an_explicit_threshold(self):
        """Factory 1: the explicit-parameter branch of __init__."""
        with pytest.raises(ValueError):
            VoronovEIArrhenius(A=(3.0e-8, "cm^3/(molecule*s)"), P=0.25, X=1.8, K=0.42)
        with pytest.raises(ValueError):
            VoronovEIArrhenius()

        # A supplied threshold is kept verbatim, including the degenerate zero.
        assert VoronovEIArrhenius(dE=13.598433).dE_eV == 13.598433
        assert VoronovEIArrhenius(dE=0.0).dE_eV == 0.0

    def test_voronov_factory_table_constructor_takes_the_threshold_from_the_table(self):
        """Factory 2: the (Z, N) branch of __init__."""
        obj = VoronovEIArrhenius(Z=11, N=11, yaml_path_or_obj=self.table)
        assert obj.dE_eV == self.table_dE

    def test_voronov_factory_from_yaml_takes_the_threshold_from_the_table(self):
        """Factory 3: the from_yaml classmethod."""
        obj = VoronovEIArrhenius.from_yaml(self.table, Z=11, N=11)
        assert obj.dE_eV == self.table_dE

        # The rest of the object has to be table-derived too, and carry the flags the
        # constructor sets -- from_yaml must not become a path that skips them.
        assert obj.P.value_si == 0.25
        assert obj.X.value_si == 1.80
        assert obj.K.value_si == 0.42
        assert obj.uses_electron_temperature is True
        assert obj.uses_electron_density is True
        assert "Z=11" in obj.comment and "N=11" in obj.comment

    def test_voronov_factory_table_paths_refuse_a_row_without_a_threshold(self):
        """
        Both table-driven factories: a row with no dE has to fail, not fall back.
        This is what a fabricated placeholder would eventually leak through.
        """
        table = copy.deepcopy(self.table)
        del table["coefficients"][0]["entries"][0]["dE"]

        with pytest.raises((KeyError, ValueError, TypeError)):
            VoronovEIArrhenius(Z=11, N=11, yaml_path_or_obj=table)
        with pytest.raises((KeyError, ValueError, TypeError)):
            VoronovEIArrhenius.from_yaml(table, Z=11, N=11)

    def test_voronov_factory_pickle_carries_the_threshold_across(self):
        """Factory 4: __reduce__, which reconstructs through __init__."""
        obj = VoronovEIArrhenius(
            A=(3.0e-8, "cm^3/(molecule*s)"), P=0.25, X=1.8, K=0.42, dE=13.598433,
            Tmin=(1.0e4, "K"), Tmax=(2.0e5, "K"),
        )
        restored = pickle.loads(pickle.dumps(obj, -1))
        assert restored.dE_eV == obj.dE_eV == 13.598433


class TestPlasmaKineticsReprRoundTrip:
    """
    ``repr`` of a kinetics object has to be loadable Python.

    RMG's Python database format persists an entry as ``repr(entry.data)`` and reads
    it back by evaluating it, so a repr that is not valid Python is silent corruption
    deferred to whoever reloads the file. Asserting the round-trip rather than the
    presence of a substring is deliberate: it tests the property that matters, so it
    catches the typo nobody thought to look for as well as the two that were found.
    """

    @staticmethod
    def _two_temperature_plasma():
        return TwoTemperaturePlasma(
            A=(1.5e-10, "cm^3/(molecule*s)"),
            n=0.5,
            Ea_g=(10.0, "kJ/mol"),
            Ea_e=(50.0, "kJ/mol"),
            T0=(347.25, "K"),
            Tmin=(300.0, "K"),
            Tmax=(3000.0, "K"),
            comment="round-trip fixture",
        )

    @staticmethod
    def _electron_collision_plasma():
        return ElectronCollisionPlasma(
            energies=([1.0, 5.0, 10.0], "eV/molecule"),
            sigma=([1.0e-20, 2.0e-20, 5.0e-21], "m^2"),
            Tmin=(300.0, "K"),
            Tmax=(3000.0, "K"),
            comment="round-trip fixture",
        )

    @staticmethod
    def _badnell_rr_arrhenius():
        return BadnellRRArrhenius(
            A=(1.0e-13, "cm^3/(molecule*s)"),
            B=0.7,
            T0=(100.0, "K"),
            T1=(1.0e5, "K"),
            C=0.3,
            T2=(2.0e5, "K"),
            Tmin=(300.0, "K"),
            Tmax=(3000.0, "K"),
            comment="round-trip fixture",
        )

    @staticmethod
    def _voronov_ei_arrhenius():
        return VoronovEIArrhenius(
            A=(3.0e-8, "cm^3/(molecule*s)"),
            P=0.25,
            X=1.80,
            K=0.42,
            dE=13.598433,
            Tmin=(1.0e4, "K"),
            Tmax=(2.0e5, "K"),
            comment="round-trip fixture",
        )

    # All four plasma kinetics classes, each with every optional field populated so
    # that the optional branches of __repr__ are exercised too.
    FIXTURES = [
        ("TwoTemperaturePlasma", _two_temperature_plasma),
        ("ElectronCollisionPlasma", _electron_collision_plasma),
        ("BadnellRRArrhenius", _badnell_rr_arrhenius),
        ("VoronovEIArrhenius", _voronov_ei_arrhenius),
    ]

    @staticmethod
    def _database_local_context():
        """
        The namespace an RMG database file is evaluated in: the kinetics classes by
        name, as `rmgpy.data.kinetics` supplies them.
        """
        import rmgpy.kinetics.arrhenius as arrhenius_module

        return dict(vars(arrhenius_module))

    @pytest.mark.parametrize("name,build", FIXTURES, ids=[name for name, _ in FIXTURES])
    def test_repr_roundtrip_is_loadable_python(self, name, build):
        original = build.__func__() if isinstance(build, staticmethod) else build()
        text = repr(original)
        assert text.startswith(name + "(")

        try:
            restored = eval(text, self._database_local_context())
        except SyntaxError as error:
            pytest.fail("{0}.__repr__ is not valid Python: {1}\n  {2}".format(name, error, text))
        except TypeError as error:
            pytest.fail(
                "{0}.__repr__ does not name the constructor's parameters: {1}\n  {2}".format(name, error, text)
            )

        assert isinstance(restored, type(original))
        assert original.is_identical_to(restored), (
            "{0} did not survive repr -> eval intact:\n  {1}".format(name, text)
        )
        assert restored.is_identical_to(original)

    def test_repr_roundtrip_is_loadable_python_preserves_reference_temperature(self):
        """
        `is_identical_to` compares quantities with a tolerance, so pin the two scalars
        whose exact value a lossy repr would quietly round away.
        """
        original = self._two_temperature_plasma()
        restored = eval(repr(original), self._database_local_context())
        assert restored.T0.value_si == original.T0.value_si == 347.25

    def test_repr_roundtrip_is_loadable_python_preserves_ionization_threshold(self):
        original = self._voronov_ei_arrhenius()
        restored = eval(repr(original), self._database_local_context())
        assert restored.dE_eV == original.dE_eV == 13.598433


class TestTwoTemperaturePlasmaIdentity:
    """
    ``is_identical_to`` has to answer False whenever the rate differs.

    ``KineticsModel.is_identical_to`` compares only Tmin/Tmax and does not even check
    the type, so every subclass owns the whole comparison. Parametrising over the
    rate-bearing fields is the point: the live defect was a missing T0, and the next
    one will be a different field added without a matching comparison.
    """

    BASE_KWARGS = dict(
        A=(1.5e-10, "cm^3/(molecule*s)"),
        n=0.5,
        Ea_g=(10.0, "kJ/mol"),
        Ea_e=(50.0, "kJ/mol"),
        T0=(300.0, "K"),
        Tmin=(300.0, "K"),
        Tmax=(3000.0, "K"),
        comment="identity fixture",
    )

    # Every constructor field that enters k(T, Te), with a perturbation big enough to
    # clear the tolerance in Quantity.equals.
    RATE_FIELDS = [
        ("A", (4.5e-10, "cm^3/(molecule*s)")),
        ("n", 1.5),
        ("Ea_g", (25.0, "kJ/mol")),
        ("Ea_e", (90.0, "kJ/mol")),
        ("T0", (1.0, "K")),
    ]

    @pytest.mark.parametrize("field,value", RATE_FIELDS, ids=[field for field, _ in RATE_FIELDS])
    def test_identity_is_sensitive_to_every_rate_bearing_field(self, field, value):
        base = TwoTemperaturePlasma(**self.BASE_KWARGS)
        perturbed_kwargs = dict(self.BASE_KWARGS)
        perturbed_kwargs[field] = value
        perturbed = TwoTemperaturePlasma(**perturbed_kwargs)

        # Guard the premise of the assertion below: if the perturbation did not move
        # the rate, the identity check has nothing to be sensitive to and this test
        # would pass for the wrong reason.
        k_base = base.get_rate_coefficient_two_temp(1000.0, 2000.0)
        k_perturbed = perturbed.get_rate_coefficient_two_temp(1000.0, 2000.0)
        assert not np.isclose(k_base, k_perturbed, rtol=1e-9), (
            "perturbing {0} did not change k, so this case proves nothing".format(field)
        )

        assert not base.is_identical_to(perturbed), (
            "TwoTemperaturePlasma.is_identical_to ignores {0}: two objects with "
            "different rates compare identical".format(field)
        )
        assert not perturbed.is_identical_to(base)

    def test_identity_is_sensitive_to_type(self):
        """A different class is never identical, whatever the base class compares."""
        base = TwoTemperaturePlasma(**self.BASE_KWARGS)
        other = Arrhenius(A=(1.5e-10, "cm^3/(molecule*s)"), n=0.5, Ea=(10.0, "kJ/mol"),
                          Tmin=(300.0, "K"), Tmax=(3000.0, "K"))
        assert not base.is_identical_to(other)

    def test_identity_holds_for_an_unperturbed_copy(self):
        """The other half of the property: equal parameters still compare identical."""
        base = TwoTemperaturePlasma(**self.BASE_KWARGS)
        twin = TwoTemperaturePlasma(**self.BASE_KWARGS)
        assert base.is_identical_to(twin)
        assert twin.is_identical_to(base)
        assert base.is_identical_to(pickle.loads(pickle.dumps(base, -1)))
