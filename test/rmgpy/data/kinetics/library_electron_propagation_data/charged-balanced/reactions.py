#!/usr/bin/env python
# encoding: utf-8

name = "charged-balanced"
shortDesc = u"Fixture library: a charge-transfer entry that declares its electron count in its kinetics"
longDesc = u"""
Not chemistry to be used in a mechanism. This library exists so that a test can prove that
KineticsLibrary.load propagates the electron count from an entry's rate law onto the Reaction it
builds, the way KineticsDepository.load already does.

The single entry carries a net charge of +1 on the reactant side and 0 on the product side, and
its rate law declares electrons = -1. It therefore balances if and only if that declared count
reaches Reaction.is_balanced, which folds one consumed electron (charge -1) onto the reactant
side. The electron is left implicit, as the shipped electrochemistry libraries do.

The rate parameters are arbitrary but well-formed. See the sibling fixture library
"charged-unbalanced", which is the same reaction with a declared count that does NOT close the
charge; together they show that the entry here is accepted on its merits and not because the
balance check was weakened.
"""

entry(
    index = 1,
    label = "proton + CH2X => CH3X",
    reversible = False,
    kinetics = SurfaceChargeTransfer(
        A = (2.5e10, 'm^3/(mol*s)'),
        n = 0.0,
        Ea = (0.75, 'eV/molecule'),
        V0 = (0.0, 'V'),
        alpha = 0.5,
        electrons = -1,
        Tmin = (298, 'K'),
        Tmax = (3000, 'K'),
    ),
    shortDesc = u"""Balanced once the declared electron is counted""",
)
