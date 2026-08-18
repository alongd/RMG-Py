#!/usr/bin/env python
# encoding: utf-8

name = "Charge_Agnostic_Attachment/groups"
shortDesc = u"test fixture: the same family without the neutral-reactant declaration"
longDesc = u"""
Byte-for-byte the group tree, template, recipe and `allowChargedSpecies` of
`Neutral_Only_Attachment`, minus the `allowChargedReactants` line.

It is the control, and it is what makes the other fixture's result mean anything: the two
families differ in one declaration and nothing else, so if this one still generates
`O2- + e- => O2(2-)` and the other does not, the declaration is what refused it -- not the group,
not the recipe, not some incidental property of the anion.

It is also the regression guard on the default. `allowChargedReactants` left undeclared must
inherit `allowChargedSpecies`, or every existing family that legitimately consumes an ion -- the
`Cation_*` and `Surface_Proton_Electron_Reduction_*` families -- silently loses reactions.
"""

template(reactants=["Attacher"], products=["Anion"], ownReverse=False)

reversible = False

reactantNum = 1
productNum = 1

allowChargedSpecies = True

electrons = -1

recipe(actions=[
    ['LOSE_RADICAL', '*1', 1],
    ['GAIN_PAIR', '*1', 1],
])

entry(
    index = 1,
    label = "Attacher",
    group =
"""
1 *1 O u1 p2 c0
""",
    kinetics = None,
)
