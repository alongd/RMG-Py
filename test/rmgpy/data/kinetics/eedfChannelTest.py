#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
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

from pathlib import Path

from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.kinetics import EEDFChannel


LIBRARY = "eedf-channel"
FIXTURE_ROOT = Path(__file__).with_name("eedf_channel_data")


def _load_library(root):
    database = KineticsDatabase()
    database.load_libraries(str(root), libraries=[LIBRARY])
    return database.libraries[LIBRARY]


def _assert_fixture_channel(channel):
    assert isinstance(channel, EEDFChannel)
    assert channel.process == "Ar -> Ar*"
    assert channel.collision_set == "argon-lxcat-v1"
    assert channel.side == "ine"
    assert channel.uses_eedf is True
    assert channel.Tmin.value_si == 300
    assert channel.Tmax.value_si == 5000
    assert channel.Pmin.value_si == 1000
    assert channel.Pmax.value_si == 1_000_000
    assert channel.comment == "fixture EEDF channel"


def test_eedf_channel_library_loads_from_database_execution_context():
    library = _load_library(FIXTURE_ROOT)

    _assert_fixture_channel(library.entries[1].data)
    reaction = library.entries[1].item
    assert reaction.reactants[0].molecule[0].electronic_state == ""
    assert reaction.products[0].molecule[0].electronic_state == "Ar(4s)"
    assert not reaction.reactants[0].is_isomorphic(reaction.products[0])


def test_eedf_channel_library_save_and_reload_round_trip(tmp_path):
    library = _load_library(FIXTURE_ROOT)
    output_root = tmp_path / "saved"
    output_library = output_root / LIBRARY
    output_library.mkdir(parents=True)
    library.save(str(output_library / "reactions.py"))

    reloaded = _load_library(output_root)

    reloaded_channel = next(iter(reloaded.entries.values())).data
    _assert_fixture_channel(reloaded_channel)
    assert repr(reloaded_channel) == repr(library.entries[1].data)
