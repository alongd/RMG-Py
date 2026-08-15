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
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
The four plasma kinetics classes are exported from ``rmgpy.kinetics``, which is not the same
thing as being loadable from a kinetics library. A library's ``reactions.py`` is evaluated by
``KineticsLibrary.load`` with ``exec(content, global_context, local_context)``, where
``local_context`` is the hand-maintained dict built in ``KineticsDatabase.__init__``. A class
missing from that dict is invisible to every library, and because ``__builtins__`` is set to
``None`` for that ``exec``, the failure does not even surface as a ``NameError`` -- the name
lookup falls through to the builtins mapping and raises
``TypeError: 'NoneType' object is not subscriptable`` at the entry's line.

These tests therefore go through ``KineticsDatabase.load_libraries``, the same call the
production loader makes, rather than eval'ing an entry by hand, and check that what comes back
carries the parameters the library file specified. Membership of ``local_context`` is
deliberately not asserted: it would pass while the load still failed for some other reason.

Two constraints on how a plasma reaction may be written surfaced while building the fixture, and
they shape the entries in ``plasma_local_context_data/plasma-local-context/reactions.py``:

1. ``Reaction.is_balanced`` counts the electron as a conserved element, so the physically honest
   forms of radiative recombination (``e + proton => H``) and electron-impact ionization
   (``e + H => proton + e + e``) are both rejected as unbalanced. Entries 3 and 4 leave the
   electron implicit instead, as the electrochemistry libraries do -- it is supplied by the rate
   law, which is what ``uses_electron_density`` on those two classes is for.
2. ``Reaction.is_balanced`` accumulates ``reactants_net_charge`` and ``products_net_charge`` and
   then returns without ever comparing them, so charge is not in fact enforced. That is what
   lets an electron-implicit entry such as ``H => proton`` load at all.
"""

import os.path

from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.kinetics import (
    BadnellRRArrhenius,
    ElectronCollisionPlasma,
    TwoTemperaturePlasma,
    VoronovEIArrhenius,
)


class TestPlasmaKineticsLocalContext:
    """
    Each of the four plasma kinetics classes must be usable by name in a kinetics library.
    """

    @classmethod
    def setup_class(cls):
        """
        Load the fixture library through the normal database path, once.
        """
        cls.library_root = os.path.join(os.path.dirname(__file__), "plasma_local_context_data")
        database = KineticsDatabase()
        database.load_libraries(cls.library_root, libraries=["plasma-local-context"])
        cls.library = database.libraries["plasma-local-context"]
        cls.entries = cls.library.entries

    def _loaded(self, index):
        """
        Return the kinetics object of the entry at `index`, checking the entry loaded at all.
        """
        assert index in self.entries, f"entry {index} did not load; got {sorted(self.entries)}"
        return self.entries[index].data

    def _assert_round_trips(self, index, expected):
        """
        The loaded kinetics must be the same class as, and identical to, `expected`, which the
        test constructs directly with the parameters the library file specifies.
        """
        loaded = self._loaded(index)
        assert isinstance(loaded, type(expected)), f"entry {index} loaded as {type(loaded).__name__}"
        # is_identical_to is this repo's equality idiom for kinetics models; none of them define
        # __eq__, so `==` would compare identity and fail for every class, plasma or not.
        assert loaded.is_identical_to(expected), f"entry {index}: {loaded!r} != {expected!r}"
        assert expected.is_identical_to(loaded), f"entry {index}: identity is not symmetric"
        assert repr(loaded) == repr(expected)

    def test_two_temperature_plasma_entry_loads(self):
        """A TwoTemperaturePlasma entry loads and keeps its parameters."""
        self._assert_round_trips(
            1,
            TwoTemperaturePlasma(
                A=(1.23e13, "cm^3/(mol*s)"),
                n=0.5,
                Ea_g=(12.0, "kJ/mol"),
                Ea_e=(450.0, "kJ/mol"),
                T0=(300.0, "K"),
                Tmin=(300, "K"),
                Tmax=(5000, "K"),
            ),
        )

    def test_electron_collision_plasma_entry_loads(self):
        """An ElectronCollisionPlasma entry loads and keeps its cross-section table."""
        self._assert_round_trips(
            2,
            ElectronCollisionPlasma(
                energies=([0.0, 5.0, 10.0, 20.0], "eV/molecule"),
                sigma=([0.0, 1.1e-21, 4.4e-21, 2.2e-21], "m^2"),
                Tmin=(300, "K"),
                Tmax=(5000, "K"),
            ),
        )

    def test_badnell_rr_arrhenius_entry_loads(self):
        """A BadnellRRArrhenius entry loads and keeps its parameters, optional C/T2 included."""
        self._assert_round_trips(
            3,
            BadnellRRArrhenius(
                A=(3.44e-13, "cm^3/(molecule*s)"),
                B=0.7,
                T0=(4.5, "K"),
                T1=(1.7e6, "K"),
                C=0.11,
                T2=(3.1e5, "K"),
                Tmin=(300, "K"),
                Tmax=(20000, "K"),
            ),
        )

    def test_voronov_ei_arrhenius_entry_loads(self):
        """A VoronovEIArrhenius entry loads and keeps its parameters, dE included."""
        self._assert_round_trips(
            4,
            VoronovEIArrhenius(
                A=(2.91e-8, "cm^3/(molecule*s)"),
                P=0.0,
                X=0.232,
                K=0.39,
                dE=(13.6, "eV"),
                Tmin=(300, "K"),
                Tmax=(20000, "K"),
            ),
        )

    def test_all_four_entries_loaded_as_library_reactions(self):
        """
        The library as a whole is usable: every entry became a LibraryReaction carrying plasma
        kinetics. A per-entry pass would not catch a library that loaded but produced nothing.
        """
        reactions = self.library.get_library_reactions()
        assert len(reactions) == 4, f"expected 4 library reactions, got {len(reactions)}"
        loaded_classes = {type(reaction.kinetics) for reaction in reactions}
        assert loaded_classes == {
            TwoTemperaturePlasma,
            ElectronCollisionPlasma,
            BadnellRRArrhenius,
            VoronovEIArrhenius,
        }, f"got {sorted(cls.__name__ for cls in loaded_classes)}"
