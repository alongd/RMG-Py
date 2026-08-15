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

# Entries 3 and 4 leave the electron implicit, as the electrochemistry libraries do: the
# electron is supplied by the rate law (both classes set uses_electron_density) rather than
# written as a species. It is not a stylistic choice here -- Reaction.is_balanced treats the
# electron as a conserved element, so the explicit forms "e + proton => H" and
# "e + H => proton + e + e" are both rejected as unbalanced. See the module docstring of
# test/rmgpy/data/kinetics/plasmaLocalContextTest.py.
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
        Tmin = (300, 'K'),
        Tmax = (20000, 'K'),
    ),
    shortDesc = u"""BadnellRRArrhenius in a library entry""",
)

entry(
    index = 4,
    label = "H => proton",
    reversible = False,
    kinetics = VoronovEIArrhenius(
        A = (2.91e-8, 'cm^3/(molecule*s)'),
        P = 0.0,
        X = 0.232,
        K = 0.39,
        dE = (13.6, 'eV'),
        Tmin = (300, 'K'),
        Tmax = (20000, 'K'),
    ),
    shortDesc = u"""VoronovEIArrhenius in a library entry""",
)
