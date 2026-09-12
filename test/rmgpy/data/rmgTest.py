#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                  #
#                                                                             #
###############################################################################

"""
This module contains regression tests for :class:`rmgpy.data.rmg.RMGDatabase.load`,
specifically for the `surface` keyword argument.

`RMGDatabase.load()` used to call `self.load_thermo(...)` only inside `if surface:`,
so `load(..., surface=False)` -- an entirely reasonable thing to do for a gas-phase
run -- silently left `self.thermo` as `None`, with nothing raised at load time. See
docs/i217-thermo-surface-gate/report.md for the full history and the git archaeology
(regression introduced in commit 191253b486, which replaced an unconditional
`load_thermo(...)` call with the `if surface:`-gated one it inherited from surface-only
metal-database loading).

Both directions are exercised in the SAME process deliberately, in this order:
`surface=False` first, then `surface=True`. `RMGDatabase.__init__` sets the
`rmgpy.data.rmg.database` module global to `self` (`global database; database = self`),
so constructing a second `RMGDatabase` after the first does not corrupt anything this
test reads -- each `RMGDatabase` instance keeps its own `self.thermo` -- but the ordering
here is deliberate and documented rather than incidental: `surface=False` is checked
before anything else touches an `RMGDatabase.thermo` in this process.

Note on `db.thermo.groups` (a contradiction, discovered empirically, of the ticket
brief that dispatched this test): `ThermoDatabase.load_groups()` always loads the
`adsorptionPt111`/`adsorptionLi` group categories into `self.groups`, unconditionally
of the `surface` argument -- only `ThermoDatabase.load_surface()` (loading the
`MetalDatabase` into `self.surface`) is gated by `surface`. So the value that actually
distinguishes `surface=False` from `surface=True` is `db.thermo.surface` (`{}` vs.
`{'metal': MetalDatabase(...)}`), not the *presence* of a key in `db.thermo.groups`.
This test asserts on `db.thermo.surface`, which is the value that is actually gated.
"""

import os

import rmgpy
from rmgpy.data.rmg import RMGDatabase


def _load(surface):
    db = RMGDatabase()
    db.load(
        path=rmgpy.settings['database.directory'],
        thermo_libraries=['primaryThermoLibrary'],
        reaction_libraries=[],
        kinetics_families='none',
        kinetics_depositories=[],
        depository=False,
        solvation=False,
        surface=surface,
        testing=True,
    )
    return db


class RMGDatabaseLoadSurfaceTest:
    """
    Regression test for i217: `RMGDatabase.load(..., surface=False)` must still
    populate `self.thermo` -- thermo is not a surface concern.
    """

    def test_load_surface_false_populates_thermo_without_adsorption_groups_loaded(self):
        db = _load(surface=False)

        # The core defect: `db.thermo` must not be the sentinel `None` set by
        # `RMGDatabase.__init__`. Assert on the value, not on "did it raise".
        assert db.thermo is not None, (
            "RMGDatabase.load(surface=False) left db.thermo as None -- thermo was "
            "never loaded because load_thermo() was only called inside `if surface:`."
        )

        # Thermo groups and libraries are genuinely populated (not just a bare,
        # empty ThermoDatabase object).
        assert 'group' in db.thermo.groups
        assert len(db.thermo.groups['group'].entries) > 0
        assert 'primaryThermoLibrary' in db.thermo.libraries

        # The adsorption-groups switch (ThermoDatabase.load's `surface` argument,
        # forwarded from RMGDatabase.load_thermo's `surface` parameter) must still
        # be OFF: no metal database loaded into db.thermo.surface.
        assert db.thermo.surface == {}, (
            f"surface=False must not load the metal/adsorption surface database, "
            f"got db.thermo.surface = {db.thermo.surface!r}"
        )

    def test_load_surface_true_populates_thermo_with_adsorption_groups_loaded(self):
        db = _load(surface=True)

        assert db.thermo is not None
        assert 'group' in db.thermo.groups
        assert len(db.thermo.groups['group'].entries) > 0
        assert 'primaryThermoLibrary' in db.thermo.libraries

        # surface=True must load the metal database (today's default behaviour).
        assert 'metal' in db.thermo.surface
        assert db.thermo.surface['metal'] is not None
