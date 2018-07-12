#!/usr/bin/env python
# -*- coding: utf-8 -*-

title = 'HSO2'

modelChemistry = "CCSD(T)-F12/cc-pVTZ-F12"
levelOfTheory = "CCSD(T)-F12a/cc-pVTZ-F12//B3LYP/6-311G(2d,d,p), rotors at B3LYP/6-311G(2d,pd)"
author = 'alongd'
frequencyScaleFactor = 0.975

useHinderedRotors = True
useBondCorrections = False

species('HSOO','HSOO.py')

species('HSO2','HSO2.py')

species(
    label = 'N2',
    structure = SMILES('N#N'),
    E0 = (-8.69489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0135,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(810.913,'J/mol'), sigma=(3.621,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.61263,-0.00100893,2.49898e-06,-1.43375e-09,2.58634e-13,-1051.1,2.6527], Tmin=(100,'K'), Tmax=(1817.04,'K')), NASAPolynomial(coeffs=[2.97592,0.00164138,-7.19708e-07,1.25375e-10,-7.91505e-15,-1025.86,5.53744], Tmin=(1817.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.69489,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: BurkeH2O2"""),
)

transitionState('TS1','TS1/TS1.py')

reaction(
    label='HSO2  = HSOO', reactants=['HSO2'], products=['HSOO'],
    transitionState = 'TS1',
    tunneling = 'Eckart',
)

kinetics(
label = 'HSO2  = HSOO',
Tmin = (300,'K'), Tmax = (3000,'K'), Tcount = 15,
)
