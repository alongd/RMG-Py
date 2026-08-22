#!/usr/bin/env python
# encoding: utf-8

name = "electron-precedence"
shortDesc = "Fixture depository: charge-transfer kinetics whose own electron count must win"
longDesc = """
Not chemistry. Three entries of the same synthetic methyl recombination, each with a charge-transfer
rate law that carries its own electron count -- one of each class that the old isinstance allowlist
in KineticsDepository.load missed (SurfaceChargeTransferBEP, ArrheniusChargeTransfer,
ArrheniusChargeTransferBM). Each declared count differs from the family default (-1) the test loads
the depository with, so the test can tell whether the data-borne count was kept or overwritten.
"""

entry(
    index = 0,
    label = "CH3 + CH3 <=> C2H6",
    kinetics = SurfaceChargeTransferBEP(
        A = (1.0e13, 'm^2/(mol*s)'), n = 0.0, E0 = (10.0, 'kJ/mol'),
        V0 = (0.0, 'V'), alpha = 0.5, electrons = 2,
    ),
    shortDesc = u"""SurfaceChargeTransferBEP carries electrons=2""",
)

entry(
    index = 1,
    label = "CH3 + CH3 <=> C2H6",
    kinetics = ArrheniusChargeTransfer(
        A = (1.0e13, 'm^3/(mol*s)'), n = 0.0, Ea = (10.0, 'kJ/mol'),
        V0 = (0.0, 'V'), alpha = 0.5, electrons = -3,
    ),
    shortDesc = u"""ArrheniusChargeTransfer carries electrons=-3""",
)

entry(
    index = 2,
    label = "CH3 + CH3 <=> C2H6",
    kinetics = ArrheniusChargeTransferBM(
        A = (1.0e13, 'm^3/(mol*s)'), n = 0.0, w0 = (1.0e5, 'J/mol'), E0 = (10.0, 'kJ/mol'),
        V0 = (0.0, 'V'), alpha = 0.5, electrons = 2,
    ),
    shortDesc = u"""ArrheniusChargeTransferBM carries electrons=2""",
)
