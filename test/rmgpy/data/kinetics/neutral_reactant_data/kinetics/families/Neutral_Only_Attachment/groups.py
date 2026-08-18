#!/usr/bin/env python
# encoding: utf-8

name = "Neutral_Only_Attachment/groups"
shortDesc = u"test fixture: electron attachment restricted to neutral reactants"
longDesc = u"""
A test fixture, not chemistry. It exists to exercise `allowChargedReactants`, and it is
deliberately built the way a real family must not be: the root is one wide group,

    1 *1 O u1 p2 c0

which matches the radical oxygen of the *anion* `[O-][O]` just as readily as it matches either
oxygen of neutral `[O][O]`. That is the point. A group constrains the atom it is written on;
`c0` on *1 says nothing about the charge of the molecule *1 belongs to, because RMG groups match
subgraphs. So this root cannot express "the reactant is neutral" no matter how it is written, and
without an engine-side control the family generates `O2- + e- => O2(2-)`.

The committed `Plasma_Electron_Attachment` family avoids that today only by accident of shape --
its root is a narrow union that happens to exclude `O2-` -- which is exactly the fragile
arrangement this fixture refuses to reproduce. Widening the root here is what makes the engine
control, rather than the group shape, the only thing that can refuse the anion.

Paired with `Charge_Agnostic_Attachment`, which is this file minus the one declaration.
"""

template(reactants=["Attacher"], products=["Anion"], ownReverse=False)

reversible = False

reactantNum = 1
productNum = 1

# Charged *products* are the whole point of an attachment family; charged *reactants* are outside
# what it claims to model. Saying both with one flag is impossible, which is why the second one
# exists.
allowChargedSpecies = True
allowChargedReactants = False

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
