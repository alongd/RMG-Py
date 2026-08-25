#!/usr/bin/env python
# encoding: utf-8

name = "Electron_Reverse_Recombination/training"
shortDesc = "Fixture training reaction stored backward, for the reverse electron-sign test"
longDesc = """
The same cation recombination as electron_propagation_data, but its single training reaction is
written in the dissociation direction -- CH3Li -> Li+ + CH3, releasing an electron -- which is
backward relative to the family's recombination template. `get_training_set` therefore has to flip
it to the forward direction before returning it. The family declares `electrons = -1`, so the
flipped forward reaction must carry the negated stored count, -1; before the fix the flip built a
bare `Reaction()` that dropped the count to 0 entirely.

The entry has to carry its own electron count, and cannot lean on the family declaration the way
the forward fixture does. `KineticsDepository.load` stamps the family-*forward* count onto whatever
orientation is stored and then balances that, so a backward-stored entry handed the family's -1
comes out mismatched by 2*|electrons| and is rejected at load -- always, for any nonzero
declaration. Declaring +1 on the kinetics is what is physically true of the stored orientation: the
dissociation releases the electron the recombination consumed.
"""

entry(
    index = 0,
    label = "CH3Li <=> Lip + CH3",
    degeneracy = 1.0,
    kinetics = ArrheniusChargeTransfer(
        A = (1.0e+13, 's^-1'),
        n = 0,
        Ea = (0, 'kJ/mol'),
        V0 = (0, 'V'),
        alpha = 0.5,
        electrons = 1,
    ),
    rank = 5,
    shortDesc = u"""Fixture value, not a measurement""",
)
