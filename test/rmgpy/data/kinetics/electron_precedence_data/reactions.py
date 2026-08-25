#!/usr/bin/env python
# encoding: utf-8

name = "electron-precedence"
shortDesc = "Fixture depository: charge-transfer kinetics whose own electron count must win"
longDesc = """
Three entries of the same two-electron cation recombination, `Li+ + H+ (+ 2e-) -> LiH`, each with a
charge-transfer rate law from a different class -- one of each that the old isinstance allowlist in
KineticsDepository.load missed (SurfaceChargeTransferBEP, ArrheniusChargeTransfer,
ArrheniusChargeTransferBM). Holding the reaction fixed and varying only the class is the point: the
class is the only thing that can explain a difference in outcome.

Every entry declares -2, which is the count the chemistry actually has (+2 on the reactant side,
neutral product) and which differs from the family default of -1 that the test loads the depository
with, so the test can tell whether the data-borne count was kept or overwritten. Because
`is_balanced` is live, an overwritten count is not a wrong number but a failed load.

The rates and their units are placeholders. Nothing here reads them.
"""

entry(
    index = 0,
    label = "Lip + Hp <=> LiH",
    kinetics = SurfaceChargeTransferBEP(
        A = (1.0e13, 'm^2/(mol*s)'), n = 0.0, E0 = (10.0, 'kJ/mol'),
        V0 = (0.0, 'V'), alpha = 0.5, electrons = -2,
    ),
    shortDesc = u"""SurfaceChargeTransferBEP carries electrons=-2""",
)

entry(
    index = 1,
    label = "Lip + Hp <=> LiH",
    kinetics = ArrheniusChargeTransfer(
        A = (1.0e13, 'm^3/(mol*s)'), n = 0.0, Ea = (10.0, 'kJ/mol'),
        V0 = (0.0, 'V'), alpha = 0.5, electrons = -2,
    ),
    shortDesc = u"""ArrheniusChargeTransfer carries electrons=-2""",
)

entry(
    index = 2,
    label = "Lip + Hp <=> LiH",
    kinetics = ArrheniusChargeTransferBM(
        A = (1.0e13, 'm^3/(mol*s)'), n = 0.0, w0 = (1.0e5, 'J/mol'), E0 = (10.0, 'kJ/mol'),
        V0 = (0.0, 'V'), alpha = 0.5, electrons = -2,
    ),
    shortDesc = u"""ArrheniusChargeTransferBM carries electrons=-2""",
)
