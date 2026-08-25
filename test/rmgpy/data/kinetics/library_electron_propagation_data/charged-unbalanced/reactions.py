#!/usr/bin/env python
# encoding: utf-8

name = "charged-unbalanced"
shortDesc = u"Fixture library: a charge-transfer entry whose declared electron count does not close the charge"
longDesc = u"""
Not chemistry to be used in a mechanism. This is the sibling of the "charged-balanced" fixture
library: the same reaction, but the rate law declares electrons = -2 where the charges only close
at -1.

It exists as a tripwire. The correct fix for the missing electron propagation makes the
"charged-balanced" entry load; the wrong fix -- weakening, exempting or special-casing
Reaction.is_balanced so it stops asking about charge -- would make this entry load too. This
library must keep being rejected at load, so that "the charged entry loads" can only mean "it is
balanced", never "it was excused".
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
        electrons = -2,
        Tmin = (298, 'K'),
        Tmax = (3000, 'K'),
    ),
    shortDesc = u"""Off by one electron: must stay rejected""",
)
