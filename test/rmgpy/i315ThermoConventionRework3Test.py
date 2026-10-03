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

"""Honest wrapper/helper mistakes and exact-byte pin regressions."""

import hashlib
from pathlib import Path

import pytest

from rmgpy.data.thermo import ThermoDatabase, ThermoLibrary
from rmgpy.thermo import ThermoData


DATA = Path(__file__).parent / 'test_data'
FIXTURES = DATA / 'i315_thermo_convention_rework3'


def _load(path, context):
    return ThermoLibrary(label=path.stem).load(
        str(path), context.local_context, context.global_context)


def test_helper_uses_own_thermo_data_wrapper_without_affecting_next_library():
    """The helper must apply its library's 1530 kJ/mol convention conversion."""
    context = ThermoDatabase()
    context.load_libraries(str(FIXTURES), libraries=[
        str(FIXTURES / 'helper_conversion.py'), str(FIXTURES / 'next_library.py')])
    converted = context.libraries['helper_conversion']
    subsequent = context.libraries['next_library']
    assert converted.thermo_convention == 'ion'
    # Check isolation before the assertion that reproduces the silent corruption.
    assert subsequent.entries['argon'].data.H298.value_si == 0
    assert subsequent.thermo_convention is None
    assert context.local_context['ThermoData'] is ThermoData
    assert converted.entries['proton'].data.H298.value_si == pytest.approx(1530000)


def test_global_thermo_convention_declaration_is_preserved():
    library = _load(FIXTURES / 'global_declaration.py', ThermoDatabase())
    assert library.thermo_convention == 'ion'
    assert library.entries['proton'].data.H298.value_si == pytest.approx(1530000)


def test_crlf_edit_invalidates_raw_byte_legacy_pin(tmp_path):
    reviewed = DATA / 'i315_thermo_convention/reviewed_legacy.py'
    raw = reviewed.read_bytes()
    assert b'\r\n' not in raw
    context = ThermoDatabase()
    clean = _load(reviewed, context)
    assert clean.thermo_convention == 'ion'
    assert clean._loaded_file_sha256 == hashlib.sha256(raw).hexdigest()
    crlf_raw = raw.replace(b'\n', b'\r\n')
    crlf = tmp_path / 'changed_line_endings.py'
    crlf.write_bytes(crlf_raw)
    changed = _load(crlf, context)
    assert changed._loaded_file_sha256 == hashlib.sha256(crlf_raw).hexdigest()
    assert changed._loaded_file_sha256 != clean._loaded_file_sha256
    assert changed.thermo_convention is None
    assert changed.entries['[Arp]'].data.H298.value_si == clean.entries['[Arp]'].data.H298.value_si
