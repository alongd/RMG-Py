#!/usr/bin/env python
# encoding: utf-8

name = "surface-charge-transfer-bep"
shortDesc = u"Fixture library exercising SurfaceChargeTransferBEP in the loader namespace"
longDesc = u"""
Not chemistry to be used in a mechanism. This library exists so that a test can prove that
SurfaceChargeTransferBEP can be named in a kinetics library entry and evaluated by the database
loader, which execs this file against KineticsDatabase.local_context. Before the loader knew the
class, this entry could not be evaluated at all. The rate parameters below are arbitrary but
well-formed; the test compares what loads against an object constructed in Python with the same
values. The electron is left implicit (the rate law carries it via electrons=-1), as the
electrochemistry libraries do, because Reaction.is_balanced treats the electron as a conserved
element and would reject an explicit form.
"""

entry(
    index = 1,
    label = "proton + CH2X => CH3X",
    reversible = False,
    kinetics = SurfaceChargeTransferBEP(
        A = (2.483e21, 'cm^3/(mol*s)'),
        n = 0.0,
        E0 = (10.0, 'kJ/mol'),
        V0 = (0.0, 'V'),
        alpha = 0.5,
        electrons = -1,
        Tmin = (300, 'K'),
        Tmax = (3000, 'K'),
    ),
    shortDesc = u"""SurfaceChargeTransferBEP in a library entry""",
)
