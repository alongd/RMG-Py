#!/usr/bin/env python
# encoding: utf-8

name = "Plain_Recombination/training"
shortDesc = "Control fixture: the same training reaction in a family that declares no electrons"
longDesc = """
Identical to Electron_Carrying_Recombination/training except that the family's groups.py carries no
`electrons` declaration. It is the control: the reaction must load with `electrons = 0` and keep the
molecularity it has always had, so that the propagation cannot be mistaken for a blanket change.
"""

entry(
    index = 0,
    label = "CH3 + CH3 <=> C2H6",
    degeneracy = 1.0,
    kinetics = Arrhenius(A=(1.0e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'kcal/mol'), T0=(1, 'K')),
    rank = 5,
    shortDesc = u"""Fixture value, not a measurement""",
    longDesc =
u"""
Second order, two reactant particles, cm^3/(mol*s) -- the ordinary case.
""",
)
