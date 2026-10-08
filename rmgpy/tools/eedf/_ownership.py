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

"""Strict ownership conversion for public EEDF inputs and outputs."""

from types import MappingProxyType

import numpy as np


def _plain_numeric_dtype(dtype):
    """Return whether *dtype* contains only binary64-or-narrower numbers."""
    dtype = np.dtype(dtype)
    if dtype.metadata is not None:
        raise TypeError('NumPy dtype metadata is not supported')
    if dtype.hasobject or dtype.subdtype is not None:
        return False
    if dtype.fields is not None:
        return all(
            field[0].fields is None and _plain_numeric_dtype(field[0])
            for field in dtype.fields.values()
        )
    if dtype.kind == 'b':
        return dtype.itemsize <= 1
    if dtype.kind in ('i', 'u', 'f'):
        return dtype.itemsize <= 8
    if dtype.kind == 'c':
        return dtype.itemsize <= 16
    return False


def _detached_copy(value, name='value', *, coordinates=False):
    """Own an allow-listed value, refusing every other concrete type.

    Supported values are exact ``None``, ``bool``, ``int``, ``float``,
    ``complex`` and ``str``; exact ``dict``, ``list``, ``tuple`` and
    ``mappingproxy`` containers containing supported values; exact NumPy
    numeric/bool arrays (including numeric structured arrays without object
    fields); and exact binary64-or-narrower NumPy numeric/bool scalars.
    Container and scalar subclasses, object/string/datetime arrays, extended
    precision diagnostic values and all other types are refused.

    With ``coordinates=True``, a real NumPy scalar coordinate, including an
    extended-precision or zero-dimensional value, is instead normalized to a
    Python binary64 ``float``. Diagnostic values are always copied without
    coordinate normalization or refused.
    """
    value_type = type(value)
    if ((value_type is np.ndarray or isinstance(value, np.generic))
            and value.dtype.metadata is not None):
        raise TypeError(name + ' contains unsupported NumPy dtype metadata')
    if value is None or value_type in (bool, int, float, complex, str):
        return value
    if value_type in (dict, MappingProxyType):
        result = {}
        for key, item in value.items():
            copied_key = _detached_copy(
                key, name + ' key', coordinates=coordinates)
            copied_item = _detached_copy(
                item, '{}[{!r}]'.format(name, copied_key),
                coordinates=coordinates)
            try:
                result[copied_key] = copied_item
            except TypeError as error:
                raise TypeError(
                    name + ' contains an unsupported mapping key') from error
        return result
    if value_type is list:
        return [
            _detached_copy(item, name + ' item', coordinates=coordinates)
            for item in value
        ]
    if value_type is tuple:
        return tuple(
            _detached_copy(item, name + ' item', coordinates=coordinates)
            for item in value
        )
    if isinstance(value, np.generic) and value_type is value.dtype.type:
        dtype = value.dtype
        if coordinates and dtype.fields is None and dtype.kind in ('i', 'u', 'f'):
            return float(value)
        if _plain_numeric_dtype(dtype) and dtype.fields is None:
            if dtype.kind == 'b':
                return bool(value)
            if dtype.kind in ('i', 'u'):
                return int(value)
            if dtype.kind == 'f':
                return float(value)
            if dtype.kind == 'c':
                return complex(value)
    if value_type is np.ndarray:
        if (coordinates and value.ndim == 0 and value.dtype.fields is None
                and value.dtype.kind in ('i', 'u', 'f')):
            return float(value)
        if _plain_numeric_dtype(value.dtype):
            return np.array(value, copy=True)
    qualified_type = value_type.__module__ + '.' + value_type.__qualname__
    if value_type is np.ndarray:
        qualified_type += '[dtype={}]'.format(value.dtype)
    raise TypeError(
        '{} contains unsupported value type {}'.format(name, qualified_type))
