#!/usr/bin/env python
# encoding: utf-8

name = "Singlet_O2_to_triplet/training"
shortDesc = u"Kinetics used to train group additivity values"
longDesc = u"""

"""

entry(
    index = 1,
    label = "O2(S) <=> O2(T)",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.0e+11, 's^-1'), n=0, Ea=(0, 'kcal/mol')),
    rank = 3,
    shortDesc = u"""""",
    longDesc =
u"""
taken from:
FD.R. Kearns, Chem. Rev., 1971, 71(4), 395-427, doi: 10.1021/cr60272a004
(Table II, rough average of the different compounds) and adjusted to a first order reaction at 1 atm
""",
)
