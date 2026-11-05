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

"""Chemical dissociation rates covered by the existing N2 energy inventory.

The Cosby/Zipf event curves are not equivalent to sing-sup. Complete-set
energy coverage is an explicit approximation, never a new Boltzmann loss.
"""

import copy
import ast
import hashlib
from pathlib import Path
import re

import numpy as np

from rmgpy.tools.eedf.schema import SpecError, FingerprintMismatch, _content_hash


SOURCE_SHA256 = 'f2d2433c76c6e6667dc843069a53c0d7a7f78a5e8e183300ae967b1f7c63230a'
PROCESSES = (
    'e + N2(X,v=0) -> e + N(4S) + N(4S), Excitation',
    'e + N2(X,v=0) -> e + N(4S) + N(2D), Excitation',
)
APPROXIMATION = 'complete-N2-loss-inventory; event curves and thresholds differ from sing-sup'


def validate_rate_only_marker(role, side, branch_product):
    """Keep the approved chemical role forward-only and outside cascades."""
    if role not in ('collision', 'rate_only'):
        raise SpecError('unsupported EEDF owner role')
    if role == 'rate_only' and (side != 'ine' or branch_product is not None):
        raise SpecError('rate-only owner must be forward and cannot be a cascade')


def energy_schema(channels):
    """Describe chemical rates separately from the unchanged energy operator."""
    energy = [c for c in channels if c.get('role', 'collision') == 'collision']
    chemical = [c for c in channels if c.get('role') == 'rate_only']
    if len(energy) + len(chemical) != len(channels):
        raise SpecError('unsupported EEDF owner role')
    coverage = [c['description'] for c in energy if c['description'].startswith('e + N2(')]
    if chemical and not any('sing-sup' in name for name in coverage):
        raise SpecError('rate-only coverage requires the complete N2 inventory including sing-sup')
    seen = set()
    for channel in channels:
        from rmgpy.tools.eedf.ownership import iter_reaction_owners
        for owner in iter_reaction_owners(channel.get('reaction')):
            key = (owner['library'], owner['index'])
            if key in seen:
                raise SpecError('duplicate chemical owner in energy schema')
            seen.add(key)
    records = []
    for channel in chemical:
        digest = channel.get('source_sha256')
        if (channel['description'] not in PROCESSES or channel.get('classification') != 'A'
                or channel.get('kind') != 'excitation' or channel.get('reaction') is None
                or not isinstance(digest, str) or not re.fullmatch('[0-9a-f]{64}', digest)
                or '<->' in channel['description']):
            raise SpecError('invalid rate-only dissociation channel')
        try:
            marker = ast.parse(channel['reaction']['repr'], mode='eval').body
            if not isinstance(marker, ast.Call) or not isinstance(marker.func, ast.Name) or marker.func.id != 'EEDFChannel' or marker.args:
                raise ValueError('marker constructor')
            fields = {item.arg: ast.literal_eval(item.value) for item in marker.keywords}
            if (fields.get('role') != 'rate_only' or fields.get('side', 'ine') != 'ine'
                    or fields.get('process') != channel['description']
                    or not fields.get('collision_set') or fields.get('branch_product') is not None):
                raise ValueError('marker role/process')
        except (ValueError, SyntaxError, AttributeError, TypeError) as error:
            raise SpecError('rate-only owner differs from energy schema') from error
        records.append({name: copy.deepcopy(channel[name]) for name in
                        ('description', 'role', 'reaction', 'source_sha256', 'cross_section', 'threshold_eV')})
    return {'schema': 'eedf-chemical-energy/1',
            'energy_operator_sha256': _content_hash(energy),
            'chemical_rates': records,
            'coverage': {'approximation': APPROXIMATION, 'channels': coverage}}


def validate_energy_schema(channels, declaration):
    """Refuse missing, incomplete or stale chemical/energy metadata."""
    if declaration != energy_schema(channels):
        raise SpecError('chemical rates differ from energy schema or complete coverage')


def read_dissociation_source(path, expected=SOURCE_SHA256):
    """Select exactly the two approved blocks from the pinned extra file.

    Read/hash/parse one byte snapshot. The source uses 2N(4S); expand that
    spelling solely in the chemical process identity, retaining every datum.
    """
    data = Path(path).read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    if digest != expected or digest != SOURCE_SHA256:
        raise FingerprintMismatch('rate-only dissociation source')
    text = data.decode()
    result = []
    for process, threshold in zip(PROCESSES, (10., 12.)):
        source_process = process.replace('N(4S) + N(4S)', '2N(4S)')
        matches = [block for block in re.split(r'(?m)^EXCITATION\s*$', text)
                   if '[' + source_process + ']' in block]
        if len(matches) != 1:
            raise SpecError('missing or duplicate dissociation source branch')
        block = matches[0]
        match = re.search(r'PARAM\.:\s+E =\s*(\S+)\s+eV', block)
        if match is None or float(match[1]) != threshold:
            raise SpecError('dissociation source threshold')
        parts = re.split(r'(?m)^-{3,}\s*$', block)
        values = np.asarray([[float(x) for x in line.split()]
                             for line in parts[1].splitlines() if line.strip()])
        if values.ndim != 2 or values.shape[1] != 2 or len(values) < 2:
            raise SpecError('dissociation source data')
        result.append(dict(description=process, kind='excitation', threshold_eV=threshold,
                           cross_section={'energy_eV': values[:, 0].tolist(),
                                          'sigma_m2': values[:, 1].tolist()}))
    return result


def validate_spec_sources(spec, channels):
    """Keep rate sources out of the physical collision-file closure."""
    declaration = spec.get('rate_only')
    chemical = [c for c in channels if c.get('role') == 'rate_only']
    for item in spec['input_files'].values():
        if item['sha256'] == SOURCE_SHA256:
            raise SpecError('rate-only input cannot be a normal collision input')
    if not chemical:
        if declaration is not None:
            raise SpecError('rate-only declaration has no chemical rates')
        return
    if (not isinstance(declaration, dict) or set(declaration) != {'sources', 'energy_schema'}
            or set(declaration['sources']) != {SOURCE_SHA256}
            or len(chemical) != 2 or {c['description'] for c in chemical} != set(PROCESSES)):
        raise SpecError('rate-only source/branch coverage is incomplete')
    validate_energy_schema(channels, declaration['energy_schema'])
    sections = read_dissociation_source(declaration['sources'][SOURCE_SHA256], SOURCE_SHA256)
    for channel in chemical:
        source = next(s for s in sections if s['description'] == channel['description'])
        if channel['source_sha256'] != SOURCE_SHA256 or any(channel[k] != v for k, v in source.items()):
            raise FingerprintMismatch('rate-only channel differs from source')


def validate_artifact(channels, inputs):
    """Check rate metadata even for a table lacking physical population fields."""
    chemical = [c for c in channels if c.get('role') == 'rate_only']
    declaration = inputs.get('rate_only')
    known_processes = any(c.get('description') in PROCESSES for c in channels)
    if not chemical and declaration is None and not known_processes:
        return
    if (not chemical or not isinstance(declaration, dict)
            or set(declaration) != {'source_sha256', 'energy_schema'}
            or declaration['source_sha256'] != [SOURCE_SHA256]
            or len(chemical) != 2 or {c['description'] for c in chemical} != set(PROCESSES)
            or any(c.get('source_sha256') != SOURCE_SHA256 for c in chemical)):
        raise SpecError('rate-only artifact energy schema/source coverage')
    validate_energy_schema(channels, declaration['energy_schema'])


def check_zero_energy_rows(channels, group):
    """Never trust an artifact's claim that its chemical columns carry no heat."""
    indices = [i for i, c in enumerate(channels) if c.get('role') == 'rate_only']
    if indices and len(group):
        for name in ('channel_power', 'k_sup', 'product_fractions', 'attachment_energy_eV'):
            values = np.asarray(group[name])
            if np.any(values[..., indices] != 0.):
                raise FingerprintMismatch('rate-only column carries energy or an inverse: ' + name)


def refresh_spec(base_spec, database, source, destination):
    """Add the two sourced chemical columns without changing physical inputs."""
    import json
    import tempfile
    from rmgpy.tools.eedf.schema import load_spec, _file_hash
    from rmgpy.tools.eedf.channels import _load_channel_map, _validate_physical_map
    from rmgpy.tools.eedf.channel_map import _load_library

    spec = load_spec(base_spec)
    if 'rate_only' in spec:
        raise SpecError('rate-only owners already present')
    channels = _load_channel_map(spec['channel_map'])
    sections = read_dissociation_source(source)
    library_name = 'PlasmaAtomicNitrogen'
    with tempfile.TemporaryDirectory(prefix='dissociation-owner-') as directory:
        library, hashes = _load_library(Path(database), Path(directory), library_name)
        for index, section, label in zip((18, 19), sections, (
                'e- + N2 => e- + N + N', 'e- + N2 => e- + N + N(2D)')):
            entry = library.entries.get(index)
            if (entry is None or entry.label != label or entry.item.reversible
                    or not entry.item.is_balanced() or entry.data.role != 'rate_only'
                    or entry.data.side != 'ine' or entry.data.process != section['description']
                    or entry.data.collision_set != 'IST-Lisbon-N2-dissociation-extra'):
                raise SpecError('missing or mismatched dissociation library owner')
            channels.append(dict(section, classification='A', role='rate_only',
                                 source_sha256=SOURCE_SHA256,
                                 target_fraction='state:N2(X,v=0)', product_fraction=0.,
                                 sigma_max_m2=max(section['cross_section']['sigma_m2']),
                                 flux_group='N2-dissociation-chemical',
                                 reaction={'library': library_name, 'index': index, 'repr': repr(entry.data)}))
    spec['rate_only'] = {'sources': {SOURCE_SHA256: str(Path(source).resolve())},
                         'energy_schema': energy_schema(channels)}
    _validate_physical_map(channels, spec, {name: axis[0] for name, axis in spec['axes'].items() if name != 'u'})
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=False)
    channel_path = destination / 'channel-map.json'
    channel_path.write_text(json.dumps(channels, indent=2) + '\n')
    spec['channel_map'] = {'path': str(channel_path), 'sha256': _file_hash(channel_path)}
    spec['output_root'], spec['scratch_root'] = str(destination / 'tables'), str(destination / 'loki')
    path = destination / 'atomic-dissociation-spec.json'
    path.write_text(json.dumps(spec, indent=2) + '\n')
    load_spec(path)
    (destination / 'conversion-report.json').write_text(json.dumps({
        'source_sha256': SOURCE_SHA256, 'library_sha256': hashes,
        'energy_operator_sha256': spec['rate_only']['energy_schema']['energy_operator_sha256'],
        'chemical_owners': [18, 19], 'node_count_unchanged': True,
        'approximation': APPROXIMATION}, indent=2) + '\n')
    return path


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('base-spec', 'database', 'source', 'output-directory'):
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args()
    print(refresh_spec(args.base_spec, args.database, args.source, args.output_directory))
