#!/usr/bin/env python
# encoding: utf-8

name = "Electron_Reverse_Recombination/rules"
shortDesc = ""
longDesc = """
One rule at the root, so the family has a rate estimator at all. Nothing in the reverse-flip test
reads it; the number is a placeholder, not a measurement.
"""
entry(
    index=1,
    label="Root",
    kinetics=ArrheniusEP(
        A=(1e13, "cm^3/(mol*s)"),
        n=0,
        alpha=0,
        E0=(0, "kcal/mol"),
        Tmin=(300, "K"),
        Tmax=(1500, "K"),
    ),
    rank=0,
    shortDesc="""Default""",
)
