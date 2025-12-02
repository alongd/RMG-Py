#!/usr/bin/env python
# encoding: utf-8

name = "Plasma_Radiative_Recombination/groups"
shortDesc = u""
longDesc = u"""
This RR family describes reactions of the sort:

A+ + e- => A  + hv
A  + e- => A- + hv

*1 is always the reactive atom, N and Z are determined by it.
"""

template(reactants=["A", "e-"], products=["A-"], ownReverse=False)

reverse = "photoionization"

reversible = False
custom_kinetics = True
allowChargedSpecies = True
product_electrons = 0

recipe(actions=[
    ['GAIN_RADICAL', '*1', 1],
    ['LOSE_CHARGE', '*1', 1],
])

entry(
    index = 0,
    label = "A",
    group =
"""
1 *1 R ux px c[0,+1,+2,+3,+4,+5,+6,+7,+8]
""",
    kinetics = None,
)

entry(
    index = 1,
    label = "e-",
    group =
"""
1 *2 e u1 p0 c-1
""",
    kinetics = None,
)

tree(
"""
L1: A
L1: e-
"""
)
