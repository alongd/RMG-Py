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

"""Convert a frozen argon map without guessing chemistry for missing owners."""

import argparse
from collections import Counter, defaultdict
import copy
import hashlib
import json
from pathlib import Path
import re
import sys
import tempfile

import yaml

from rmgpy.kinetics import EEDFChannel
from rmgpy.exceptions import DatabaseError
from rmgpy.tools.eedf.channels import _load_channel_map, _validate_map, _validate_physical_map
from rmgpy.tools.eedf.schema import FingerprintMismatch, SpecError, load_spec
from rmgpy.tools.eedf.ownership import load_cascade_yields, make_branch_family


class MissingOwners(SpecError):
    """A complete diagnostic report, without a partial runnable map."""

    def __init__(self, report):
        self.report = report
        super().__init__('class-A channels lack owning PlasmaArgon reactions: ' +
                         ', '.join(str(row['id']) for row in report['missing_owners']))


class _UniqueLoader(yaml.SafeLoader):
    """Ambiguous YAML keys cannot silently replace classification metadata."""

    def construct_mapping(self, node, deep=False):
        result = {}
        for key_node, value_node in node.value:
            key = self.construct_object(key_node, deep=deep)
            if not isinstance(key, (str, int)) or key in result:
                raise SpecError('duplicate or invalid frozen YAML key: ' + repr(key))
            result[key] = self.construct_object(value_node, deep=deep)
        return result


def _read_pinned_bytes(path, digest):
    try:
        data = Path(path).read_bytes()
    except OSError as error:
        raise SpecError('cannot read pinned input: ' + str(error)) from error
    if hashlib.sha256(data).hexdigest() != digest:
        raise FingerprintMismatch('input hash: ' + str(path))
    return data


def _frozen_rows(frozen):
    """Check metadata completeness, leaving scientific class policy in the YAML."""
    if not isinstance(frozen, dict) or frozen.get('schema') != 'eedf-channel-map/2':
        raise SpecError('frozen channel map schema')
    if not isinstance(frozen.get('loki_files'), dict) or not frozen['loki_files']:
        raise SpecError('frozen channel map loki_files')
    gaps = frozen.get('class_a_gaps')
    if not isinstance(gaps, dict):
        raise SpecError('frozen channel map class_a_gaps')
    gap_status = {}
    for description, ids in gaps.items():
        match = re.fullmatch(r'\[([^]]+)\] .+', description) if isinstance(description, str) else None
        if not match or not isinstance(ids, list):
            raise SpecError('malformed class_a_gaps entry')
        for channel_id in ids:
            if type(channel_id) is not int or channel_id in gap_status:
                raise SpecError('duplicate/invalid class_a_gaps id')
            gap_status[channel_id] = match[1]
    channels = frozen.get('channels')
    if not isinstance(channels, list) or not channels:
        raise SpecError('empty frozen channel map')
    by_process, ids = {}, set()
    for row in channels:
        if not isinstance(row, dict):
            raise SpecError('frozen channel record')
        channel_id, process = row.get('id'), row.get('loki_process')
        if (type(channel_id) is not int or channel_id < 1 or channel_id in ids
                or not isinstance(process, str) or process in by_process
                or row.get('class') not in ('A', 'B', 'C', 'D')
                or not isinstance(row.get('lxcat_file'), str)
                or row.get('lxcat_file') not in frozen['loki_files']):
            raise SpecError('frozen channel identity/class/file')
        ids.add(channel_id)
        by_process[process] = row
        owner = row.get('rmg')
        if row['class'] == 'A':
            if (not isinstance(owner, dict) or owner.get('library') != 'PlasmaArgon'
                    or not isinstance(owner.get('reaction'), str)
                    or not owner['reaction']
                    or owner.get('status') not in
                    ('present', 'present_lumped', 'absent', 'lumped_internal')):
                raise SpecError('class-A reaction metadata: ' + str(channel_id))
            status = owner['status']
            if (status == 'present' and channel_id in gap_status or
                    status != 'present' and gap_status.get(channel_id) != status):
                raise SpecError('class_a_gaps disagrees with channel ' + str(channel_id))
            if status == 'present_lumped' and not owner.get('lumped'):
                raise SpecError('missing lumping declaration: ' + str(channel_id))
        elif owner is not None or channel_id in gap_status:
            raise SpecError('non-A channel has reaction metadata: ' + str(channel_id))
    if not set(gap_status) <= ids:
        raise SpecError('class_a_gaps contains unknown ids')
    counts = {label: sum(row['class'] == label for row in channels) for label in 'ABCD'}
    if frozen.get('class_counts') != counts:
        raise SpecError('frozen class_counts disagrees with channels')
    return by_process


def _owner_problem(row, library):
    owner = row['rmg']
    if owner['status'] in ('absent', 'lumped_internal'):
        return owner['status']
    index = owner.get('index')
    if type(index) is not int or index not in library.entries:
        return 'requested library entry is absent'
    from rmgpy.util import parse_reaction_equation
    try:
        left, right, reversible = parse_reaction_equation(owner['reaction'])
    except DatabaseError as error:
        raise SpecError('invalid requested reaction: ' + str(row['id'])) from error
    reaction = library.entries[index].item
    if (Counter(left.split(' + ')) != Counter(s.label for s in reaction.reactants)
            or Counter(right.split(' + ')) != Counter(s.label for s in reaction.products)
            or reversible != reaction.reversible):
        return 'requested reaction differs from loaded library entry'
    return None


def _convert_channel(physical, frozen, collision_set):
    """Keep each collision's physics, including distinct lumped-channel thresholds."""
    result = copy.deepcopy(physical)
    result['classification'] = frozen['class']
    result['reaction'] = None
    if frozen['class'] == 'A':
        owner = frozen['rmg']
        if (owner['status'] not in ('present', 'present_lumped')
                or type(owner.get('index')) is not int):
            raise SpecError('class-A channel has no owner: ' + str(frozen['id']))
        marker = EEDFChannel(process=frozen['loki_process'],
                             collision_set=collision_set, side='ine')
        result['reaction'] = {'library': owner['library'], 'index': owner['index'],
                              'repr': repr(marker)}
    return result


def _branch_family(row, library, collision_set, cascade_yields):
    """Discover declared markers, never infer branches from placeholder equations."""
    owners = []
    for index, entry in library.entries.items():
        marker = getattr(entry, 'data', None)
        if marker is None:
            marker = getattr(entry.item, 'kinetics', None)
        if getattr(marker, 'uses_eedf', False) and marker.process == row['loki_process']:
            owners.append((row['rmg']['library'], index, entry.item, marker))
    branched = [owner for owner in owners if owner[3].branch_product is not None]
    if not branched:
        sides = [owner[3].side for owner in owners]
        if len(sides) != len(set(sides)):
            raise SpecError('more than one unweighted EEDF owner would double-count one tabulated rate')
        return None
    if len(branched) != len(owners):
        raise SpecError('cascade family cannot share a collision with an unweighted owner')
    declared = row['rmg']
    if (declared['status'] in ('present', 'present_lumped')
            and declared.get('index') not in [owner[1] for owner in branched]):
        raise SpecError('cascade family overlaps the frozen ordinary owner')
    if cascade_yields is None:
        raise SpecError('cascade family requires a sourced yield table')
    return make_branch_family(cascade_yields, row['id'], row['loki_process'],
                              row['threshold_eV'], collision_set, branched)


def convert_channel_map(spec, frozen, library, collision_set, cascade_yields=None):
    """Convert complete physical coverage or report every unresolved selected A row."""
    if not isinstance(collision_set, str) or not collision_set.strip():
        raise SpecError('explicit collision_set required')
    by_process = _frozen_rows(frozen)
    channels = _load_channel_map(spec['channel_map'])
    _validate_physical_map(channels, spec, {
        name: values[0] for name, values in spec['axes'].items() if name != 'u'})
    selected_files = {}
    for name, entry in spec['input_files'].items():
        _read_pinned_bytes(entry['path'], entry['sha256'])
        if entry['kind'] == 'cross_section':
            filename = Path(name).name
            if filename in selected_files:
                raise SpecError('ambiguous cross-section basename: ' + filename)
            if frozen['loki_files'].get(filename) != entry['sha256']:
                raise FingerprintMismatch('frozen collision file: ' + filename)
            selected_files[filename] = entry['sha256']
    selected = [row for row in frozen['channels'] if row['lxcat_file'] in selected_files]
    if {row['loki_process'] for row in selected} != {row['description'] for row in channels}:
        raise SpecError('frozen map does not cover selected physical collisions')
    report = {'outcome': 'converted', 'selected_channel_count': len(selected),
              'class_counts': dict(sorted(Counter(row['class'] for row in selected).items())),
              'selected_ids': sorted(row['id'] for row in selected),
              'missing_owners': [], 'shared_owners': [], 'branch_families': []}
    shared = defaultdict(list)
    converted = []
    for physical in channels:
        row = by_process[physical['description']]
        if physical['threshold_eV'] != row.get('threshold_eV'):
            raise SpecError('frozen threshold differs from physical channel: ' + str(row['id']))
        if row['class'] == 'A':
            try:
                family = _branch_family(row, library, collision_set, cascade_yields)
                problem = None if family is not None else _owner_problem(row, library)
            except (SpecError, FingerprintMismatch) as error:
                family, problem = None, str(error)
            if problem:
                process = row['loki_process']
                report['missing_owners'].append({
                    'id': row['id'], 'process': process, 'classification': 'A',
                    'products': re.split(r' <-> | -> ', process)[-1].rsplit(', ', 1)[0],
                    'threshold_eV': row['threshold_eV'],
                    'class_a_gaps_status': row['rmg']['status'],
                    'requested_owner': copy.deepcopy(row['rmg']), 'reason': problem})
                continue
            if family is not None:
                channel = copy.deepcopy(physical)
                channel.update(classification='A', reaction=family)
                converted.append(channel)
                report['branch_families'].append({'id': row['id'], 'family': copy.deepcopy(family)})
                for owner in family['branches']:
                    shared[(owner['library'], owner['index'])].append(row['id'])
                continue
            shared[(row['rmg']['library'], row['rmg']['index'])].append(row['id'])
        converted.append(_convert_channel(physical, row, collision_set))
    report['missing_owners'].sort(key=lambda row: row['id'])
    report['branch_families'].sort(key=lambda row: row['id'])
    report['shared_owners'] = [
        {'library': library_name, 'index': index, 'ids': sorted(ids)}
        for (library_name, index), ids in sorted(shared.items()) if len(ids) > 1]
    if report['missing_owners']:
        report['outcome'] = 'refused'
        raise MissingOwners(report)
    _validate_map(converted)
    return converted, report


def _load_library(database, scratch):
    """Load a private, provenance-recorded copy through the existing RMG loader."""
    from rmgpy.data.kinetics.database import KineticsDatabase
    from rmgpy.data.kinetics.library import KineticsLibrary
    path = Path(database) / 'input/kinetics/libraries/PlasmaArgon'
    folder = scratch / 'PlasmaArgon'
    folder.mkdir()
    hashes = {}
    for name in ('reactions.py', 'dictionary.txt'):
        data = (path / name).read_bytes()
        (folder / name).write_bytes(data)
        hashes[name] = hashlib.sha256(data).hexdigest()
    # The native loader validates generation markers on the copied pair. Keep
    # the raw bytes (including markers and newlines) in both provenance and copy.
    contexts = KineticsDatabase()
    library = KineticsLibrary(label='PlasmaArgon')
    library.load(str(folder / 'reactions.py'), contexts.local_context, contexts.global_context)
    return library, hashes


def prepare_spec(base_spec, frozen_map, frozen_digest, database, collision_set, *,
                 yield_table=None, radiative_source=None, yield_table_sha256=None, source_sha256=None):
    """Snapshot external inputs; no solver is started and no source is written."""
    spec = load_spec(base_spec)
    frozen_bytes = _read_pinned_bytes(frozen_map, frozen_digest)
    try:
        frozen = yaml.load(frozen_bytes, Loader=_UniqueLoader)
    except (yaml.YAMLError, UnicodeError) as error:
        raise SpecError('malformed frozen channel map YAML') from error
    provenance = {'frozen_map_sha256': frozen_digest,
                  'base_channel_map': copy.deepcopy(spec['channel_map']),
                  'input_files': copy.deepcopy(spec['input_files']),
                  'collision_set': collision_set}
    with tempfile.TemporaryDirectory(prefix='eedf-map-') as directory:
        scratch = Path(directory)
        snapshot = copy.deepcopy(spec)
        # Validation parses these same pinned bytes, never a later source revision.
        for i, entry in enumerate(snapshot['input_files'].values()):
            data = _read_pinned_bytes(entry['path'], entry['sha256'])
            destination = scratch / (str(i) + Path(entry['path']).suffix)
            destination.write_bytes(data)
            entry['path'] = str(destination)
        data = _read_pinned_bytes(spec['channel_map']['path'], spec['channel_map']['sha256'])
        (scratch / 'channels.json').write_bytes(data)
        snapshot['channel_map']['path'] = str(scratch / 'channels.json')
        library, hashes = _load_library(database, scratch)
        provenance['library_sha256'] = hashes
        default_root = Path(database) / 'input/eedf'
        yield_path = Path(yield_table) if yield_table is not None else default_root / 'argon_cascade_yields.yaml'
        source_path = Path(radiative_source) if radiative_source is not None else default_root / 'argon_radiative_data.yaml'
        cascade = None
        # Old databases without any cascade documents remain valid for ordinary
        # conversion. Explicit pins, partial installs and requested families fail closed.
        if (yield_path.exists() or source_path.exists() or yield_table is not None
                or radiative_source is not None or yield_table_sha256 is not None or source_sha256 is not None):
            cascade = load_cascade_yields(yield_path, source_path,
                                          yield_table_sha256=yield_table_sha256, source_sha256=source_sha256)
            provenance['cascade_sources'] = {name: cascade[name] for name in
                                              ('yield_table_sha256', 'source_sha256', 'radiative_model')}
        try:
            converted, report = convert_channel_map(snapshot, frozen, library, collision_set, cascade)
        except MissingOwners as error:
            error.report['provenance'] = provenance
            raise
    report['provenance'] = provenance
    return spec, converted, report


def _write_outputs(destination, report, spec=None, channels=None):
    """Claim a fresh destination exclusively; the runnable spec is written last."""
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=False)
    files = {'conversion-report.json': report}
    if spec is not None:
        spec = copy.deepcopy(spec)
        data = json.dumps(channels, indent=2, allow_nan=False) + '\n'
        spec['channel_map'] = {'path': str(destination / 'channel-map.json'),
                               'sha256': hashlib.sha256(data.encode()).hexdigest()}
        spec['output_root'] = str(destination / 'tables')
        spec['scratch_root'] = str(destination / 'loki')
        files['channel-map.json'] = channels
        files['pure-ar-spec.json'] = spec
    # Serialize everything before any file is published. An I/O failure preserves
    # the claimed directory for inspection; a retry cannot overwrite it.
    serialized = {name: json.dumps(value, indent=2, allow_nan=False) + '\n'
                  for name, value in files.items()}
    for name, data in serialized.items():
        with (destination / name).open('x') as stream:
            stream.write(data)


def _main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ('base-spec', 'frozen-map', 'frozen-map-sha256', 'database',
                   'collision-set', 'output-directory'):
        parser.add_argument('--' + option, required=True)
    for option in ('yield-table', 'radiative-source', 'yield-table-sha256', 'source-sha256'):
        parser.add_argument('--' + option)
    args = parser.parse_args(argv)
    try:
        optional = {name: getattr(args, name) for name in
                    ('yield_table', 'radiative_source', 'yield_table_sha256', 'source_sha256')
                    if getattr(args, name) is not None}
        spec, channels, report = prepare_spec(
            args.base_spec, args.frozen_map, args.frozen_map_sha256,
            args.database, args.collision_set, **optional)
    except MissingOwners as error:
        report, spec, channels = error.report, None, None
        report['refusal'] = str(error)
    except (SpecError, FingerprintMismatch, DatabaseError, OSError, ValueError) as error:
        report = {'outcome': 'refused', 'missing_owners': [], 'refusal': str(error)}
        spec, channels = None, None
    try:
        _write_outputs(args.output_directory, report, spec, channels)
    except OSError as error:
        print('cannot publish conversion: ' + str(error), file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2))
    return 2 if report['outcome'] == 'refused' else 0


if __name__ == '__main__':
    sys.exit(_main())
