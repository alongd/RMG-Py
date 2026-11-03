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

"""Shared validation of sourced collision ownership and cascade yields.

The generated yield document owns the numbers. Consumers validate a byte
snapshot and never renormalize it or turn an unresolved family into chemistry.
"""

import copy
import hashlib
import math
from pathlib import Path
import re

import yaml

from rmgpy.tools.eedf.schema import FingerprintMismatch, SpecError


YIELD_SUM_ATOL = 1e-12
ARGON_PRODUCTS = ('Ar(1s5)', 'Ar(1s4)', 'Ar(1s3)', 'Ar(1s2)')


class _UniqueLoader(yaml.SafeLoader):
    """Refuse duplicate source keys instead of accepting the last value."""

    def construct_mapping(self, node, deep=False):
        result = {}
        for key_node, value_node in node.value:
            key = self.construct_object(key_node, deep=deep)
            if type(key) not in (str, int) or key in result:
                raise SpecError('duplicate or invalid cascade source key: ' + repr(key))
            result[key] = self.construct_object(value_node, deep=deep)
        return result


def validate_branch_yield(value):
    """Accept a finite plain numeric fraction, including an explicit zero."""
    if type(value) not in (int, float) or not 0 <= value <= 1 or not math.isfinite(value):
        raise SpecError('branch yield must be finite and between zero and one')
    return value


def validate_yield_vector(values):
    """Require all four sourced products and unity within an absolute 1e-12.

    This is the sole engine normalization rule. Return a detached copy with
    unchanged values; absent, null and zero yields have different meanings.
    """
    if not isinstance(values, dict) or set(values) != set(ARGON_PRODUCTS):
        raise SpecError('cascade yield products must contain all four argon 1s states')
    for value in values.values():
        validate_branch_yield(value)
    if abs(math.fsum(values.values()) - 1.) > YIELD_SUM_ATOL:
        raise SpecError('cascade yields must sum to one within 1e-12')
    return dict(values)


def _read_source(path, expected):
    """Hash and return the same bytes, even if the source path is replaced."""
    try:
        data = Path(path).read_bytes()
    except OSError as error:
        raise SpecError('cannot read cascade source: ' + str(error)) from error
    digest = hashlib.sha256(data).hexdigest()
    if expected is not None and expected != digest:
        raise FingerprintMismatch('cascade source hash: ' + str(path))
    return data, digest


def _validate_source_row(row):
    """Check source completeness without inferring missing transition data."""
    if not isinstance(row, dict):
        raise SpecError('invalid cascade source channel')
    if (type(row.get('id')) is not int or row['id'] < 1
            or not isinstance(row.get('process'), str) or not row['process'].strip()):
        raise SpecError('invalid cascade collision identity')
    threshold = row.get('threshold_eV')
    if type(threshold) not in (int, float) or not math.isfinite(threshold) or threshold < 0:
        raise SpecError('invalid cascade threshold')
    status = row.get('status')
    if status == 'unresolved':
        if 'product_yields' not in row or row['product_yields'] is not None:
            raise SpecError('unresolved cascade must have null yields')
    elif status == 'complete_for_E1':
        validate_yield_vector(row.get('product_yields'))
        transitions = row.get('transition_ids')
        if (row.get('missing_transitions') != [] or not isinstance(transitions, list)
                or not transitions or any(not isinstance(t, str) or not t for t in transitions)
                or len(set(transitions)) != len(transitions)):
            raise SpecError('complete cascade has missing or invalid transition provenance')
    else:
        raise SpecError('unsupported cascade source status')


def load_cascade_yields(path, radiative_source, *, yield_table_sha256=None, source_sha256=None):
    """Read a generated E1 yield document and verify its radiative source.

    Optional pins constrain the external bytes. Both files are read once;
    the returned snapshot contains their digests, never machine-local paths.
    Unresolved records may be carried as evidence, but cannot be requested
    as usable families through :func:`require_cascade_yields`.
    """
    data, digest = _read_source(path, yield_table_sha256)
    _, raw_digest = _read_source(radiative_source, source_sha256)
    try:
        document = yaml.load(data, Loader=_UniqueLoader)
    except (yaml.YAMLError, UnicodeError) as error:
        raise SpecError('malformed cascade yield YAML') from error
    if (not isinstance(document, dict)
            or set(document) != {'schema_version', 'radiative_model', 'source_sha256', 'channels'}
            or type(document['schema_version']) is not int or document['schema_version'] != 1
            or document['radiative_model'] != 'E1'):
        raise SpecError('cascade yield document schema')
    source_hash = document['source_sha256']
    if not isinstance(source_hash, str) or not re.fullmatch('[0-9a-f]{64}', source_hash):
        raise SpecError('invalid cascade radiative source hash')
    if source_hash != raw_digest:
        raise FingerprintMismatch('cascade radiative source')
    rows = document['channels']
    if not isinstance(rows, list) or not rows:
        raise SpecError('empty cascade source channels')
    ids, processes = set(), set()
    for row in rows:
        _validate_source_row(row)
        if row['id'] in ids or row['process'] in processes:
            raise SpecError('duplicate cascade collision identity')
        ids.add(row['id'])
        processes.add(row['process'])
    snapshot = copy.deepcopy(document)
    snapshot['yield_table_sha256'] = digest
    return snapshot


def require_cascade_yields(snapshot, channel_id, *, process=None, threshold_eV=None):
    """Return one complete sourced row, cross-checking its physical collision."""
    if type(channel_id) is not int or channel_id < 1:
        raise SpecError('invalid cascade channel id')
    rows = [row for row in snapshot['channels'] if row['id'] == channel_id]
    if len(rows) != 1:
        raise SpecError('missing or duplicate cascade source family: ' + str(channel_id))
    row = rows[0]
    _validate_source_row(row)
    if row['status'] != 'complete_for_E1':
        raise SpecError('unresolved cascade yields: ' + str(channel_id))
    if ((process is not None and row['process'] != process)
            or (threshold_eV is not None and row['threshold_eV'] != threshold_eV)):
        raise SpecError('cascade collision identity differs from sourced row: ' + str(channel_id))
    return copy.deepcopy(row)
