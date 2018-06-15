#!/usr/bin/env python
# encoding: utf-8

modelChemistry = "CBS-QB3"
useHinderedRotors = True
useBondCorrections = False
author = 'I.B. Modeling'

species('C2H4', 'ethene.py',
        structure=SMILES('C=C'),
        )

statmech('C2H4')

thermo('C2H4', 'NASA')
