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
``SurfaceChargeTransferBEP`` is exported from ``rmgpy.kinetics`` but was missing from the
hand-maintained ``KineticsDatabase.local_context`` (and from the import list that feeds it), so a
kinetics library that declared that rate law could not be evaluated by ``KineticsLibrary.load``:
because ``__builtins__`` is set to ``None`` for that ``exec``, the missing name does not even
surface as a ``NameError`` but as ``TypeError: 'NoneType' object is not subscriptable``.

This test goes through ``KineticsDatabase.load_libraries`` -- the same call the production loader
makes -- and checks that the entry round-trips. Membership of ``local_context`` is deliberately not
asserted: it would pass while the load still failed for some other reason. The electron is left
implicit (carried by ``electrons=-1`` on the rate law), as the electrochemistry libraries do,
because ``Reaction.is_balanced`` counts the electron as a conserved element.
"""

import os.path

from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.kinetics import SurfaceChargeTransferBEP


class TestSurfaceChargeTransferBEPLocalContext:
    """
    ``SurfaceChargeTransferBEP`` must be usable by name in a kinetics library.
    """

    def _load_entries(self):
        """Load the fixture library through the normal database path."""
        library_root = os.path.join(os.path.dirname(__file__), "surface_charge_transfer_bep_data")
        database = KineticsDatabase()
        # Pre-fix this raises TypeError: 'NoneType' object is not subscriptable, because the class
        # is absent from local_context and __builtins__ is None during the entry's exec.
        database.load_libraries(library_root, libraries=["surface-charge-transfer-bep"])
        return database.libraries["surface-charge-transfer-bep"].entries

    def test_surface_charge_transfer_bep_entry_loads(self):
        """
        The entry declaring ``SurfaceChargeTransferBEP`` must load and round-trip its parameters.
        Pre-fix, the class was invisible to the loader and the entry failed to load at all.
        """
        entries = self._load_entries()
        assert 1 in entries, f"entry 1 did not load; got {sorted(entries)}"
        loaded = entries[1].data

        expected = SurfaceChargeTransferBEP(
            A=(2.483e21, 'cm^3/(mol*s)'),
            n=0.0,
            E0=(10.0, 'kJ/mol'),
            V0=(0.0, 'V'),
            alpha=0.5,
            electrons=-1,
            Tmin=(300, 'K'),
            Tmax=(3000, 'K'),
        )

        assert isinstance(loaded, SurfaceChargeTransferBEP), \
            f"entry 1 loaded as {type(loaded).__name__}"
        # is_identical_to is this repo's equality idiom for kinetics models.
        assert loaded.is_identical_to(expected), f"entry 1: {loaded!r} != {expected!r}"
        assert loaded.electrons.value_si == -1
