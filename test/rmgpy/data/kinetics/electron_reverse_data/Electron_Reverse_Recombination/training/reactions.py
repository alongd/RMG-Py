#!/usr/bin/env python
# encoding: utf-8

name = "Electron_Reverse_Recombination/training"
shortDesc = "Fixture training reaction stored backward, for the reverse electron-sign test"
longDesc = """
The same synthetic electron-carrying methyl recombination family as electron_propagation_data, but
its single training reaction is written in the dissociation direction -- C2H6 => CH3 + CH3 -- which
is backward relative to the family's recombination template. get_training_set therefore has to flip
it to the forward direction before returning it. The family declares electrons = -1, so the flipped
forward reaction must carry the negated count (+1); before the fix the flip built a bare Reaction()
that dropped the count to 0.
"""

entry(
    index = 0,
    label = "C2H6 <=> CH3 + CH3",
    degeneracy = 1.0,
    kinetics = Arrhenius(A=(1.0e+16, 'cm^6/(mol^2*s)'), n=0, Ea=(0, 'kcal/mol'), T0=(1, 'K')),
    rank = 5,
    shortDesc = u"""Fixture value, not a measurement""",
)
