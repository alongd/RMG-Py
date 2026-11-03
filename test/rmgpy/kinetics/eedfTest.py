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


def _branch(**kwargs):
    fields = dict(process='cascade', collision_set='argon-v1',
                  branch_product='Ar(1s5)', branch_yield=0.7420494699646644,
                  yield_table_sha256='a' * 64)
    fields.update(kwargs)
    return EEDFChannel(**fields)


@pytest.mark.parametrize('restore', [deepcopy, lambda c: pickle.loads(pickle.dumps(c)),
                                    lambda c: eval(repr(c), {'EEDFChannel': EEDFChannel})])
def test_branch_roundtrip(restore):
    original = _branch()
    restored = restore(original)
    assert restored.is_identical_to(original)
    assert restored.branch_product == 'Ar(1s5)'
    assert restored.branch_yield == 0.7420494699646644
    assert restored.yield_table_sha256 == 'a' * 64


def test_branch_native_record_preserves_metadata():
    from rmgpy.rmg.pdep import _KINETICS_DATA_FIELDS
    original = _branch(branch_yield=0.)
    # KineticsModel uses native descriptors and constructor/reduce records,
    # rather than the RMGObject.as_dict API used by species.
    fields = ('process', 'collision_set', 'side', 'branch_product',
              'branch_yield', 'yield_table_sha256')
    assert set(fields) <= set(_KINETICS_DATA_FIELDS[EEDFChannel])
    record = {name: getattr(original, name) for name in fields}
    restored = EEDFChannel(**record)
    assert restored.is_identical_to(original)
    assert restored.branch_yield == 0.


@pytest.mark.parametrize('kwargs', [
    {'branch_product': None}, {'branch_yield': None}, {'yield_table_sha256': None},
    {'branch_product': 'Ars'}, {'branch_yield': -1.}, {'branch_yield': 1.01},
    {'branch_yield': float('nan')}, {'branch_yield': float('inf')},
    {'branch_yield': True}, {'branch_yield': '0.5'}, {'branch_yield': 10 ** 1000},
    {'yield_table_sha256': 'bad'}, {'side': 'sup'},
])
def test_branch_invalid_metadata_refuses(kwargs):
    with pytest.raises(ValueError, match='branch|cascade'):
        _branch(**kwargs)


def test_branch_metadata_is_atomic_and_sup_assignment_refuses():
    channel = _branch()
    with pytest.raises(ValueError, match='branch|cascade'):
        channel.side = 'sup'
    assert channel.side == 'ine'
    for field, value in [('branch_yield', .5), ('branch_product', 'Ar(1s4)'),
                         ('yield_table_sha256', 'b' * 64)]:
        with pytest.raises(AttributeError):
            setattr(channel, field, value)


@pytest.mark.parametrize('kwargs', [
    {'branch_product': 'Ar(1s4)'}, {'branch_yield': .5}, {'yield_table_sha256': 'b' * 64},
])
def test_branch_identity_includes_source_and_product(kwargs):
    assert not _branch().is_identical_to(_branch(**kwargs))


def test_ordinary_repr_and_old_pickle_arguments_are_unchanged():
    channel = EEDFChannel('old', 'argon', 'ine', None, None, None, None, '')
    assert repr(channel) == "EEDFChannel(process='old', collision_set='argon', side='ine')"
    assert channel.__reduce__()[1] == ('old', 'argon', 'ine', None, None, None, None, '')
    assert channel.branch_yield is None
    assert not _branch().is_identical_to(EEDFChannel('cascade', 'argon-v1'))
