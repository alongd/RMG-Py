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

"""HDF5 storage with separate training and independent held-out groups."""
import json
import tempfile
from pathlib import Path
from urllib.parse import quote, unquote

import h5py
import numpy as np

from rmgpy.tools.eedf.schema import canonical_json, content_hash, file_hash, SpecError, FingerprintMismatch

NUMERIC = ('u', 'EN_Td', 'k_ine', 'k_sup', 'channel_power', 'attachment_energy_eV',
           'target_fractions', 'product_fractions', 'rate_floors', 'below_floor',
           'f0', 'converged', 'iteration_count', 'energy_eV', 'energy_edges_eV')
MAPPINGS = ('swarm', 'power_groups', 'gas_fractions', 'state_populations')


def storage_key(identifier):
    """Encode one HDF5 component reversibly, including empty and dot names."""
    if not isinstance(identifier, str):
        raise SpecError('storage identifiers must be strings')
    if identifier == '':
        return '%'
    if identifier in ('.', '..'):
        return identifier.replace('.', '%2E')
    return quote(identifier, safe="()[]'+* ")


def storage_items(group):
    """Read keys in the declared encoding; legacy simple keys stay readable."""
    encoded = group.attrs.get('identifier_encoding') == 'percent-v1'
    for key, value in group.items():
        identifier = ('' if key == '%' else unquote(key)) if encoded else key
        if encoded and storage_key(identifier) != key:
            raise SpecError('noncanonical HDF5 identifier')
        yield identifier, value


def write_rows(group, rows):
    """Preserve tensor shape, numeric fields and exact symbolic assumptions."""
    rows = np.asarray(rows, dtype=object)
    first = rows.flat[0]
    for name in NUMERIC:
        if name not in first:
            continue
        values = np.asarray([row[name] for row in rows.flat])
        shape = rows.shape + values.shape[1:]
        group.create_dataset(name, data=values.reshape(shape), compression='gzip')
    for name in MAPPINGS:
        child = group.create_group(name)
        child.attrs['identifier_encoding'] = 'percent-v1'
        for key in first[name]:
            values = [row[name][key] for row in rows.flat]
            if all(isinstance(v, (int, float)) for v in values):
                child.create_dataset(storage_key(key), data=np.asarray(values).reshape(rows.shape))
            else:
                raise SpecError('symbolic population/property must be resolved: ' + key)
    for name in ('setup', 'composition', 'branch_id'):
        group.attrs[name] = canonical_json([row.get(name) for row in rows.flat])


def validate_storage(h5):
    """Refuse storage that can change the meaning of a pinned file's bytes.

    Inspect links, layouts, filters and types before reading values. Generated
    tables use unique hard links, internal contiguous/chunked numeric datasets
    and numeric/string attributes. References and unsupported types are refused
    instead of followed; plugin filters are never invoked.
    """
    seen = set()
    builtin_filters = {1, 2, 3, 4, 5, 6}  # HDF5 deflate/shuffle/checksum/SZIP/NBIT/scaleoffset.

    def refuse(reason, obj):
        raise FingerprintMismatch('HDF5 storage: ' + reason + ' at ' + obj.name)

    def visit(obj):
        address = h5py.h5o.get_info(obj.id).addr
        if address in seen:
            refuse('hard-link alias or cycle', obj)
        seen.add(address)
        for name in obj.attrs:
            dtype = obj.attrs.get_id(name).dtype
            if dtype.kind not in 'biuf' and h5py.check_string_dtype(dtype) is None:
                refuse('unsupported attribute type: ' + name, obj)
        if isinstance(obj, h5py.Group):
            for name in obj:
                if not isinstance(obj.get(name, getlink=True), h5py.HardLink):
                    refuse('soft or external link: ' + name, obj)
                visit(obj[name])
        elif isinstance(obj, h5py.Dataset):
            properties = obj.id.get_create_plist()
            if properties.get_external_count():
                refuse('external raw data', obj)
            if properties.get_layout() not in (h5py.h5d.CONTIGUOUS, h5py.h5d.CHUNKED):
                refuse('unsupported dataset layout', obj)
            for i in range(properties.get_nfilters()):
                if properties.get_filter(i)[0] not in builtin_filters:
                    refuse('external plugin filter', obj)
            if obj.dtype.kind not in 'biuf':
                refuse('unsupported dataset type', obj)
        else:
            refuse('unsupported object type', obj)

    visit(h5)


def write_artifact(root, manifest, branches, held_out, *, extra_files=None):
    """Write an immutable content-addressed artifact; return its directory."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    for rows in list(branches.values()) + ([held_out] if held_out else []):
        for row in np.asarray(rows, dtype=object).flat:
            for name in ('energy_eV', 'energy_edges_eV'):
                if name not in row or not np.array_equal(row[name], manifest[name]):
                    raise SpecError('row energy grid differs bit for bit: ' + name)
    with tempfile.TemporaryDirectory(prefix='.writing-', dir=root) as temporary:
        destination = Path(temporary)
        manifest = dict(manifest)
        manifest['branches'] = list(branches)
        path = destination / 'table.h5'
        with h5py.File(path, 'w') as h5:
            h5.attrs['schema_version'] = manifest['schema_version']
            h5.attrs['fingerprint'] = manifest['fingerprint']
            h5.attrs['qualification_sha256'] = content_hash({
                key: manifest.get(key) for key in (
                    'accepted', 'held_out_verdicts', 'screen_results',
                    'branch_detection', 'A6b-source', 'A6b-runtime',
                    'T10 interpolation qualification')})
            h5.attrs['branch_certification_sha256'] = content_hash(manifest['branch_certification'])
            h5.attrs['channel_map_content_sha256'] = content_hash(manifest['channel_map'])
            h5.attrs['policy_sha256'] = content_hash({key: manifest[key] for key in ('tolerances', 'floors', 'envelopes')})
            h5.create_dataset('energy_eV', data=manifest['energy_eV'])
            h5.create_dataset('energy_edges_eV', data=manifest['energy_edges_eV'])
            axes = h5.create_group('axes')
            axes.attrs['identifier_encoding'] = 'percent-v1'
            for name, values in manifest['axes'].items():
                axes.create_dataset(storage_key(name), data=values)
            training = h5.create_group('branches')
            training.attrs['identifier_encoding'] = 'percent-v1'
            for branch_id, rows in branches.items():
                write_rows(training.create_group(storage_key(branch_id)), rows)
            independent = h5.create_group('held_out')
            if held_out:
                write_rows(independent, held_out)
            validate_storage(h5)
        manifest['artifact_sha256'] = file_hash(path)
        final = root / manifest['artifact_sha256']
        (destination / 'manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + '\n')
        for name, contents in (extra_files or {}).items():
            if name not in ('generation_spec.json', 'validation.md'):
                raise SpecError('unsupported artifact companion file')
            (destination / name).write_text(contents)
        destination.rename(final)
        return final
