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

from copy import deepcopy
import pickle

import pytest

from rmgpy.kinetics import Arrhenius, EEDFChannel, KineticsModel


def _channel_with_metadata():
    return EEDFChannel(
        "Ar -> Ar*",
        {"source": "LXCat"},
        "sup",
        Tmin=(300, "K"),
        Tmax=(5000, "K"),
        Pmin=(0.01, "bar"),
        Pmax=(10, "bar"),
        comment="provider-owned channel",
    )


def _assert_metadata_preserved(restored, original):
    assert restored.process == original.process
    assert restored.collision_set == original.collision_set
    assert restored.side == original.side
    assert restored.Tmin.equals(original.Tmin)
    assert restored.Tmax.equals(original.Tmax)
    assert restored.Pmin.equals(original.Pmin)
    assert restored.Pmax.equals(original.Pmax)
    assert restored.comment == original.comment
    assert restored.uses_eedf is True


def test_eedf_channel_fields_and_marker_flag():
    channel = EEDFChannel(process="Ar -> Ar*", collision_set="argon-v1", side="ine")

    assert isinstance(channel, KineticsModel)
    assert channel.process == "Ar -> Ar*"
    assert channel.collision_set == "argon-v1"
    assert channel.side == "ine"
    assert channel.uses_eedf is True


def test_eedf_channel_default_constructor_supports_native_record_census():
    channel = EEDFChannel()

    assert channel.process == ''
    assert channel.collision_set == ''
    assert channel.side == 'ine'
    with pytest.raises(NotImplementedError, match='EEDF context'):
        channel.get_rate_coefficient(300.)


def test_eedf_channel_internal_marker_state_is_not_python_writable():
    channel = EEDFChannel("ionization", "argon-v1", "ine")

    with pytest.raises(AttributeError):
        channel.uses_eedf = False
    with pytest.raises(AttributeError):
        channel._side = "sup"
    with pytest.raises(AttributeError):
        _ = channel._side

    channel.side = "sup"
    assert channel.side == "sup"
    with pytest.raises(ValueError, match="must be 'ine' or 'sup'"):
        channel.side = "min"


@pytest.mark.parametrize("side", ["ine", "sup"])
def test_eedf_channel_accepts_both_rate_sides(side):
    assert EEDFChannel("ionization", "argon-v1", side).side == side


@pytest.mark.parametrize("side", [None, "", "min", "INE"])
def test_eedf_channel_rejects_unknown_rate_side(side):
    with pytest.raises(ValueError, match="must be 'ine' or 'sup'"):
        EEDFChannel("ionization", "argon-v1", side)


def test_eedf_channel_refuses_standalone_evaluation():
    channel = EEDFChannel("ionization", "argon-v1", "sup")

    with pytest.raises(NotImplementedError, match="no standalone rate coefficient"):
        channel.get_rate_coefficient(1000.0)


def test_eedf_channel_repr_is_evaluable():
    channel = _channel_with_metadata()

    restored = eval(repr(channel), {"EEDFChannel": EEDFChannel})

    assert repr(restored) == repr(channel)
    _assert_metadata_preserved(restored, channel)


def test_eedf_channel_pickle_round_trip():
    channel = _channel_with_metadata()

    restored = pickle.loads(pickle.dumps(channel))

    assert type(restored) is EEDFChannel
    _assert_metadata_preserved(restored, channel)


def test_eedf_channel_deepcopy_preserves_metadata():
    channel = _channel_with_metadata()

    restored = deepcopy(channel)

    assert restored is not channel
    assert restored.collision_set is not channel.collision_set
    _assert_metadata_preserved(restored, channel)


def test_eedf_channel_identity_requires_matching_marker_fields():
    channel = EEDFChannel("Ar -> Ar*", "argon-v1", "ine")

    assert channel.is_identical_to(EEDFChannel("Ar -> Ar*", "argon-v1", "ine"))
    assert not channel.is_identical_to(EEDFChannel("Ar -> Ar+", "argon-v1", "ine"))
    assert not channel.is_identical_to(EEDFChannel("Ar -> Ar*", "argon-v2", "ine"))
    assert not channel.is_identical_to(EEDFChannel("Ar -> Ar*", "argon-v1", "sup"))
    assert not channel.is_identical_to(Arrhenius())
