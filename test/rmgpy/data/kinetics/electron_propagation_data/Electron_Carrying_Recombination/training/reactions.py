#!/usr/bin/env python
# encoding: utf-8

name = "Electron_Carrying_Recombination/training"
shortDesc = "Fixture training reactions for the electron-count propagation tests"
longDesc = """
A deliberately ordinary, non-plasma reaction -- methyl radical recombination -- in a family that
declares `electrons = -1`. The family is synthetic: nothing about methyl recombination requires an
electron. That is the point. The electron count has to survive the load boundary for any family
that declares it, not only for the plasma families that motivated the flag.

Three reactant particles collide (two methyls and the electron), so the rate is third order and the
A factor is given in cm^6/(mol^2*s), even though `reactants` holds only two species.
"""

entry(
    index = 0,
    label = "CH3 + CH3 <=> C2H6",
    degeneracy = 1.0,
    kinetics = Arrhenius(A=(1.0e+16, 'cm^6/(mol^2*s)'), n=0, Ea=(0, 'kcal/mol'), T0=(1, 'K')),
    rank = 5,
    shortDesc = u"""Fixture value, not a measurement""",
    longDesc =
u"""
The rate is arbitrary. Only its units matter here, and they matter because they encode the
molecularity that the units checker has to arrive at.
""",
)

entry(
    index = 1,
    label = "CH3 + CH3 <=> C2H6",
    degeneracy = 1.0,
    kinetics = SurfaceChargeTransfer(
        alpha = 0.5,
        A = (1.0e+13, 'm^3/(mol*s)'),
        n = 0,
        V0 = (0, 'V'),
        Ea = (0, 'eV/molecule'),
        electrons = 1,
    ),
    rank = 5,
    shortDesc = u"""Fixture value, not a measurement""",
    longDesc =
u"""
The same reaction again, but with kinetics that carry an electron count of their own, and one that
disagrees with the family's. Charge transfer kinetics have always supplied their own count; the
family's declaration is only a default, and must not overwrite them. Nothing here is physical --
the reaction is the fixture's methyl recombination and the kinetics are borrowed wholesale -- the
entry exists solely to pin the precedence between the two sources.
""",
)
