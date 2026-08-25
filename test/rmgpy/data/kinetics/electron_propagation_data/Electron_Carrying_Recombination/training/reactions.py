#!/usr/bin/env python
# encoding: utf-8

name = "Electron_Carrying_Recombination/training"
shortDesc = "Fixture training reactions for the electron-count propagation tests"
longDesc = """
Cation recombination, the smallest reaction whose declared electron count is actually true:
`Li+ + CH3 (+ e-) -> CH3Li`. The reactant side carries +1, the product side is neutral, and the
electron the family declares is what closes the gap. The dictionary names the cation `Lip` and the
proton `Hp` because a species label may not contain a `+`; the entry labels below are split on `+`
to find the reactants.

That the chemistry is real is what arms the fixture. `Reaction.is_balanced` folds `electrons` into
the net charges, so if the family's declaration fails to reach the loaded reaction the entry does
not merely carry the wrong count -- it fails the balance check and does not load at all. The rates
are placeholders; only their units carry meaning, and only for the molecularity assertion.

Three reactant particles collide in entry 0 (the cation, the radical and the electron), so the rate
is third order and the A factor is given in cm^6/(mol^2*s), even though `reactants` holds only two
species. That mismatch between `len(reactants)` and the true molecularity is the whole reason
`get_molecularity` exists.
"""

entry(
    index = 0,
    label = "Lip + CH3 <=> CH3Li",
    degeneracy = 1.0,
    kinetics = Arrhenius(A=(1.0e+16, 'cm^6/(mol^2*s)'), n=0, Ea=(0, 'kcal/mol'), T0=(1, 'K')),
    rank = 5,
    shortDesc = u"""Fixture value, not a measurement""",
    longDesc =
u"""
The kinetics carry no electron count of their own, so this entry can only get one from the family
declaration. It is the entry the propagation tests read.
""",
)

entry(
    index = 1,
    label = "Lip + Hp <=> LiH",
    degeneracy = 1.0,
    kinetics = ArrheniusChargeTransfer(
        A = (1.0e+13, 'm^3/(mol*s)'),
        n = 0,
        Ea = (0, 'kJ/mol'),
        V0 = (0, 'V'),
        alpha = 0.5,
        electrons = -2,
    ),
    rank = 5,
    shortDesc = u"""Fixture value, not a measurement""",
    longDesc =
u"""
A two-electron recombination in the same family: two cations, +2 on the reactant side, neutral
product, and a charge-transfer rate law that carries its own electron count of -2. Charge transfer
kinetics have always supplied their own count; the family's declaration is only a default and must
not overwrite them. The entry exists to pin that precedence, and it is self-arming in the same way
entry 0 is -- if the family's -1 won, the reaction would be one charge unit short of balanced and
the load would fail.
""",
)
