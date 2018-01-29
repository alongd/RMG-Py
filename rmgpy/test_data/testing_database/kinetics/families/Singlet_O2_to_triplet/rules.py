#!/usr/bin/env python
# encoding: utf-8

name = "Singlet_O2_to_triplet/rules"
shortDesc = u""
longDesc = u"""
"""

entry(
    index = 1,
    label = "singlet",
    kinetics = Arrhenius(A=(1e+08, 's^-1'), n=0, Ea=(0, 'kcal/mol')),
    rank = 1,
    shortDesc = u"""Default""",
    longDesc =
u"""

""",
)
