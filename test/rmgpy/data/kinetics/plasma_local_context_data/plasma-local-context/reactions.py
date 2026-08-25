#!/usr/bin/env python
# encoding: utf-8

name = "plasma-local-context"
shortDesc = u"Fixture library exercising the plasma kinetics classes in the loader namespace"
longDesc = u"""
Not chemistry to be used in a mechanism. This library exists so that a test can prove that each
of the four plasma kinetics classes can be named in a kinetics library entry and evaluated by
the database loader, which execs this file against KineticsDatabase.local_context. The rate
parameters below are arbitrary but well-formed; the test compares what loads against objects
constructed in Python with the same values.
"""

entry(
    index = 1,
    label = "e + H2 => e + H + H",
    reversible = False,
    kinetics = TwoTemperaturePlasma(
        A = (1.23e13, 'cm^3/(mol*s)'),
        n = 0.5,
        Ea_g = (12.0, 'kJ/mol'),
        Ea_e = (450.0, 'kJ/mol'),
        T0 = (300.0, 'K'),
        Tmin = (300, 'K'),
        Tmax = (5000, 'K'),
    ),
    shortDesc = u"""TwoTemperaturePlasma in a library entry""",
)

entry(
    index = 2,
    label = "e + N2 => e + N + N",
    reversible = False,
    kinetics = ElectronCollisionPlasma(
        energies = ([0.0, 5.0, 10.0, 20.0], 'eV/molecule'),
        sigma = ([0.0, 1.1e-21, 4.4e-21, 2.2e-21], 'm^2'),
        Tmin = (300, 'K'),
        Tmax = (5000, 'K'),
    ),
    shortDesc = u"""ElectronCollisionPlasma in a library entry""",
)

# Entries 3 and 4 change the number of free electrons, and entries 1 and 2 do not. That is the
# whole of why they are written differently, and it is not a stylistic choice.
#
# Reaction.is_balanced compares per-element atom counts over element_list, and element_list
# begins with the electron (rmgpy/molecule/element.py). An explicit electron species is therefore
# a CONSERVED element: it can pass through a reaction but it cannot be created or destroyed.
# Entries 1 and 2 destroy none, so they write the electron out on both sides. Entries 3 and 4
# each change the count by one, so that one electron cannot be a species -- it has to be the
# signed scalar that is_balanced folds into the charge comparison, which is what the `electrons`
# field on both rate laws now supplies and what KineticsLibrary.load copies onto the reaction.
#
# Entry 4 carries both: the incident electron survives electron-impact ionization, so it is a
# spectator and is written explicitly on each side, and only the surplus electron is the scalar.
# Dropping the explicit pair would still balance and still load, but the exported equation would
# be first order in H and rmgpy.electron_balance.check_electron_reactant_order refuses it.
entry(
    index = 3,
    label = "Heplus => He",
    reversible = False,
    kinetics = BadnellRRArrhenius(
        A = (3.44e-13, 'cm^3/(molecule*s)'),
        B = 0.7,
        T0 = (4.5, 'K'),
        T1 = (1.7e6, 'K'),
        C = 0.11,
        T2 = (3.1e5, 'K'),
        electrons = -1,
        Tmin = (300, 'K'),
        Tmax = (20000, 'K'),
    ),
    shortDesc = u"""BadnellRRArrhenius in a library entry""",
)

entry(
    index = 4,
    label = "e + H => proton + e",
    reversible = False,
    kinetics = VoronovEIArrhenius(
        A = (2.91e-8, 'cm^3/(molecule*s)'),
        P = 0.0,
        X = 0.232,
        K = 0.39,
        dE = (13.6, 'eV'),
        electrons = 1,
        Tmin = (300, 'K'),
        Tmax = (20000, 'K'),
    ),
    shortDesc = u"""VoronovEIArrhenius in a library entry""",
)
