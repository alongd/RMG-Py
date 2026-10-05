#!/usr/bin/env python
# encoding: utf-8

name = "eedf-channel"
shortDesc = u"Fixture library exercising EEDFChannel persistence"
longDesc = u"""
Not chemistry to be used in a mechanism. This fixture proves that EEDFChannel
survives the kinetics-library execution context and save/load path.
"""

entry(
    index = 1,
    label = "Ar => Ar*",
    reversible = False,
    kinetics = EEDFChannel(
        process = "Ar -> Ar*",
        collision_set = "argon-lxcat-v1",
        side = "ine",
        Tmin = (300, "K"),
        Tmax = (5000, "K"),
        Pmin = (0.01, "bar"),
        Pmax = (10, "bar"),
        comment = "fixture EEDF channel",
    ),
    shortDesc = u"EEDF provider channel marker",
)
