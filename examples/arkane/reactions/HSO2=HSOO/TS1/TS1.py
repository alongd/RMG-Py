#!/usr/bin/env python
# -*- coding: utf-8 -*-

bonds = {'O=S': 2}

linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = {'CCSD(T)-F12/cc-pVTZ-F12': Log('sp.out')}

geometry = Log('freq.out')

frequencies = Log('freq.out')
