#!/usr/bin/env python
# encoding: utf-8

"""
This example loads all species data from the C2H6.yml file
To generate the .yml file, run the legacy input file in this example
"""

modelChemistry = "CBS-QB3"
useHinderedRotors = True
useBondCorrections = False

species('C2H6', 'SpeciesDatabase/C2H6.yml')

thermo('C2H6', 'NASA')
