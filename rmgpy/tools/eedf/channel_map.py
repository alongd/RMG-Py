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

"""Convert gas-specific frozen channel maps without guessing missing owners."""

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
from scipy.constants import electron_mass

from rmgpy.kinetics import EEDFChannel
from rmgpy.exceptions import DatabaseError
from rmgpy.tools.eedf.channels import (
    _classic_sections, _gas_values, _load_channel_map, _validate_map, _validate_physical_map,
)
from rmgpy.tools.eedf.integrity import _physical_properties
from rmgpy.tools.eedf.schema import FingerprintMismatch, SpecError, load_spec
from rmgpy.tools.eedf.ownership import load_cascade_yields, make_branch_family


class MissingOwners(SpecError):
    """A complete diagnostic report, without a partial runnable map."""

    def __init__(self, report):
        self.report = report
        super().__init__('class-A channels lack owning configured library reactions: ' +
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


def _frozen_rows(frozen, owner_libraries=('PlasmaArgon',)):
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
            if (not isinstance(owner, dict) or owner.get('library') not in owner_libraries
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


def _owner_problem(row, libraries):
    owner = row['rmg']
    if owner['status'] in ('absent', 'lumped_internal'):
        return owner['status']
    index = owner.get('index')
    library = libraries if hasattr(libraries, 'entries') else libraries.get(owner['library'])
    if library is None or type(index) is not int or index not in library.entries:
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


def _derived_elastic(description):
    """Return LoKI's derived Elastic identity for a source Effective process."""
    suffix = ', Effective'
    return description[:-len(suffix)] + ', Elastic' if description.endswith(suffix) else None


def _target_gas(description):
    match = re.match(r'^e \+ ([A-Z][A-Za-z0-9]*)\(', description)
    if not match:
        raise SpecError('cannot identify collision target gas: ' + description)
    return match[1]


def _declared_gases(spec):
    raw = spec.get('gas_properties', {}).get('fraction')
    if not isinstance(raw, list) or not raw:
        raise SpecError('physical gas fractions are required')
    gases = []
    for entry in raw:
        match = re.fullmatch(r'\s*([A-Z][A-Za-z0-9]*)\s+=\s+.+', entry) if isinstance(entry, str) else None
        if not match or match[1] in gases:
            raise SpecError('invalid or duplicate physical gas fraction')
        gases.append(match[1])
    return set(gases)


def _upgrade_frozen_map(frozen, gas, spec):
    """Adapt the DB's generated N2/1 map into the converter's pinned /2 form."""
    if not isinstance(frozen, dict) or frozen.get('schema') == 'eedf-channel-map/2':
        return frozen
    if gas != 'N2' or frozen.get('schema') != 'eedf-channel-map/1':
        raise SpecError('unsupported frozen channel map schema for gas ' + gas)
    if (set(frozen) != {'schema', 'set', 'source_file', 'source_sha256', 'class_counts', 'channels'}
            or not isinstance(frozen['source_file'], str)
            or not isinstance(frozen['channels'], list)):
        raise SpecError('malformed generated N2 channel map')
    source_file = Path(frozen['source_file']).name
    candidates = [(name, entry) for name, entry in spec['input_files'].items()
                  if Path(name).name == source_file]
    if len(candidates) != 1 or candidates[0][1]['sha256'] != frozen['source_sha256']:
        raise FingerprintMismatch('N2 frozen map source file: ' + source_file)
    source_name, source_entry = candidates[0]
    sections = _classic_sections(source_entry['path'], source_entry['sha256'])
    ordered = sorted(frozen['channels'], key=lambda row: row.get('id', -1))
    if (len(ordered) != len(sections)
            or [row.get('id') for row in ordered] != list(range(1, len(ordered) + 1))):
        raise SpecError('generated N2 map does not cover source blocks in order')
    rows, counts = [], Counter()
    for row, section in zip(ordered, sections):
        if row.get('threshold_eV') != section['threshold_eV']:
            raise SpecError('generated N2 map threshold differs from source block ' + str(row.get('id')))
        label, index = row.get('owner_library'), row.get('owner_index')
        owner = None
        if row.get('class') == 'A':
            if not isinstance(label, str) or type(index) is not int or not isinstance(row.get('owner_label'), str):
                raise SpecError('generated N2 class-A owner is incomplete: ' + str(row.get('id')))
            owner = {'library': label, 'index': index, 'status': 'present',
                     'reaction': row['owner_label']}
        elif row.get('class') not in ('B', 'C', 'D'):
            raise SpecError('generated N2 channel classification: ' + str(row.get('id')))
        rows.append({'id': row['id'], 'loki_process': section['description'],
                     'lxcat_file': Path(source_name).name,
                     'threshold_eV': section['threshold_eV'],
                     'class': row['class'], 'rmg': owner})
        counts[row['class']] += 1
    loki_files = {Path(source_name).name: source_entry['sha256']}
    vib_entries = [(name, entry) for name, entry in spec['input_files'].items()
                   if Path(name).name == 'N2_vib_LXCat.txt']
    if vib_entries:
        if len(vib_entries) != 1:
            raise SpecError('ambiguous N2_vib_LXCat.txt source')
        vib_name, vib_entry = vib_entries[0]
        vib_sections = _classic_sections(vib_entry['path'], vib_entry['sha256'])
        loki_files[Path(vib_name).name] = vib_entry['sha256']
        for offset, section in enumerate(vib_sections, len(rows) + 1):
            rows.append({'id': offset, 'loki_process': section['description'],
                         'lxcat_file': Path(vib_name).name,
                         'threshold_eV': section['threshold_eV'], 'class': 'B', 'rmg': None})
            counts['B'] += 1
    return {'schema': 'eedf-channel-map/2', 'set': frozen['set'],
            'loki_files': loki_files, 'class_a_gaps': {},
            'class_counts': {label: counts[label] for label in 'ABCD'},
            'channels': rows}


def _diagnostic_d6():
    return {'set_by': 'owner', 'top_fraction_limit': .05,
            'top_leak_share_limit': .05, 'K_step': 5, 'status': 'diagnostic'}


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


def convert_channel_maps(spec, gas_maps, libraries, collision_set, cascade_yields=None,
                         library_by_gas=None):
    """Convert physical coverage using the owner library declared for each gas.

    ``gas_maps`` maps gas names to frozen maps. Libraries are keyed by their
    exact database labels. Processes and owner identities cannot be claimed by
    more than one gas map.
    """
    if not isinstance(collision_set, str) or not collision_set.strip():
        raise SpecError('explicit collision_set required')
    if (not isinstance(gas_maps, dict) or not gas_maps or not isinstance(libraries, dict) or not libraries
            or any(not isinstance(value, str) or not value or not hasattr(library, 'entries')
                   for value, library in libraries.items())
            or not isinstance(library_by_gas, dict) or set(library_by_gas) != set(gas_maps)):
        raise SpecError('gas maps and configured owner libraries are required')
    declared_gases = _declared_gases(spec)
    if set(gas_maps) != declared_gases:
        extra = sorted(set(gas_maps) - declared_gases)
        if extra:
            raise SpecError('channel map gas is not declared: ' + extra[0])
        missing = sorted(declared_gases - set(gas_maps))
        if missing:
            raise SpecError('missing channel map for gas ' + missing[0])
        raise SpecError('gas map configuration differs from physical gases')
    by_process, owner_gas, frozen_file_hashes, normalized_maps = {}, {}, {}, {}
    for gas, frozen in gas_maps.items():
        if not isinstance(gas, str) or not gas:
            raise SpecError('gas-map key must be a nonempty gas name')
        if library_by_gas[gas] not in libraries:
            raise SpecError('configured owner library is not loaded for gas ' + gas)
        frozen = _upgrade_frozen_map(frozen, gas, spec)
        normalized_maps[gas] = frozen
        if not isinstance(frozen, dict) or not isinstance(frozen.get('channels'), list):
            raise SpecError('invalid frozen channel map for gas ' + gas)
        for row in frozen['channels']:
            process = row.get('loki_process', '') if isinstance(row, dict) else ''
            if _target_gas(process) != gas:
                raise SpecError('channel target gas differs from map gas ' + gas)
            owner = row.get('rmg')
            if (row.get('class') == 'A' and isinstance(owner, dict)
                    and owner.get('library') != library_by_gas[gas]):
                raise SpecError('configured owner library mismatch for gas ' + gas)
        allowed = {library_by_gas[gas]}
        rows = _frozen_rows(frozen, allowed)
        for process, row in rows.items():
            if process in by_process:
                raise SpecError('two gas maps claim one LoKI process: ' + process)
            by_process[process] = row
            owner = row.get('rmg')
            if row['class'] == 'A' and owner.get('status') in ('present', 'present_lumped'):
                identity = (owner['library'], owner.get('index'))
                previous = owner_gas.get(identity)
                if previous is not None and previous != gas:
                    raise SpecError('two gases claim one owner: ' + repr(identity))
                owner_gas[identity] = gas
        for filename, digest in frozen['loki_files'].items():
            if filename in frozen_file_hashes:
                raise SpecError('two gas maps claim one collision filename: ' + filename)
            frozen_file_hashes[filename] = digest
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
            if frozen_file_hashes.get(filename) != entry['sha256']:
                raise FingerprintMismatch('frozen collision file: ' + filename)
            selected_files[filename] = entry['sha256']
    selected = [row for frozen in normalized_maps.values() for row in frozen['channels']
                if row['lxcat_file'] in selected_files]
    if {row['loki_process'] for row in selected} != {row['description'] for row in channels}:
        raise SpecError('frozen map does not cover selected physical collisions')
    report = {'outcome': 'converted', 'selected_channel_count': len(selected),
              'class_counts': dict(sorted(Counter(row['class'] for row in selected).items())),
              'selected_ids': sorted(row['id'] for row in selected),
              'missing_owners': [], 'shared_owners': [], 'branch_families': []}
    shared = defaultdict(list)
    gas_properties, _ = _physical_properties(spec)
    mass_entry = spec['input_files'].get(gas_properties.get('mass'))
    mass_values = (_gas_values(mass_entry['path'], mass_entry['sha256'])
                   if mass_entry is not None else {})
    converted = []
    for physical in channels:
        row = by_process[physical['description']]
        library = libraries.get(row['rmg']['library']) if row['class'] == 'A' else None
        if physical['threshold_eV'] != row.get('threshold_eV'):
            raise SpecError('frozen threshold differs from physical channel: ' + str(row['id']))
        if row['class'] == 'A':
            try:
                family = _branch_family(row, library, collision_set, cascade_yields)
                problem = None if family is not None else _owner_problem(row, libraries)
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
        channel = _convert_channel(physical, row, collision_set)
        converted.append(channel)
        if physical['kind'] == 'effective':
            derived_description = _derived_elastic(physical['description'])
            target = _target_gas(physical['description'])
            if target not in mass_values or mass_values[target] <= 0:
                raise SpecError('missing physical gas mass for derived Elastic ' + target)
            derived = copy.deepcopy(channel)
            derived.update(description=derived_description, flux_group=derived_description,
                           kind='elastic', classification='B', reaction=None,
                           mass_ratio=electron_mass / mass_values[target])
            converted.append(derived)
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


def convert_channel_map(spec, frozen, library, collision_set, cascade_yields=None):
    """Backward-compatible single-argon-map conversion."""
    return convert_channel_maps(spec, {'Ar': frozen}, {'PlasmaArgon': library},
                                collision_set, cascade_yields,
                                library_by_gas={'Ar': 'PlasmaArgon'})


def _load_library(database, scratch, label='PlasmaArgon'):
    """Load a private, provenance-recorded copy through the existing RMG loader."""
    from rmgpy.data.kinetics.database import KineticsDatabase
    from rmgpy.data.kinetics.library import KineticsLibrary
    if not isinstance(label, str) or not label or Path(label).name != label:
        raise SpecError('invalid configured kinetics library label')
    path = Path(database) / 'input/kinetics/libraries' / label
    folder = scratch / label
    folder.mkdir()
    hashes = {}
    for name in ('reactions.py', 'dictionary.txt'):
        data = (path / name).read_bytes()
        (folder / name).write_bytes(data)
        hashes[name] = hashlib.sha256(data).hexdigest()
    # The native loader validates generation markers on the copied pair. Keep
    # the raw bytes (including markers and newlines) in both provenance and copy.
    contexts = KineticsDatabase()
    library = KineticsLibrary(label=label)
    library.load(str(folder / 'reactions.py'), contexts.local_context, contexts.global_context)
    return library, hashes


def _cascade_snapshot(database, *, yield_table=None, radiative_source=None,
                      yield_table_sha256=None, source_sha256=None):
    """Load the legacy pinned cascade documents for either CLI conversion path."""
    default_root = Path(database) / 'input/eedf'
    yield_path = Path(yield_table) if yield_table is not None else default_root / 'argon_cascade_yields.yaml'
    source_path = Path(radiative_source) if radiative_source is not None else default_root / 'argon_radiative_data.yaml'
    if not (yield_path.exists() or source_path.exists() or yield_table is not None
            or radiative_source is not None or yield_table_sha256 is not None or source_sha256 is not None):
        return None, None
    snapshot = load_cascade_yields(yield_path, source_path,
                                   yield_table_sha256=yield_table_sha256,
                                   source_sha256=source_sha256)
    provenance = {name: snapshot[name] for name in
                  ('yield_table_sha256', 'source_sha256', 'radiative_model')}
    return snapshot, provenance


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
        # Old databases without documents remain valid for ordinary conversion;
        # explicit pins, partial installs and requested families fail closed.
        cascade, cascade_provenance = _cascade_snapshot(
            database, yield_table=yield_table, radiative_source=radiative_source,
            yield_table_sha256=yield_table_sha256, source_sha256=source_sha256)
        if cascade_provenance is not None:
            provenance['cascade_sources'] = cascade_provenance
        try:
            converted, report = convert_channel_map(snapshot, frozen, library, collision_set, cascade)
        except MissingOwners as error:
            error.report['provenance'] = provenance
            raise
    report['provenance'] = provenance
    return spec, converted, report


def build_mixture_spec(base_spec, ladder_manifest, bindings, *,
                       composition_axis='x_N2', tv_axis='Tv_K',
                       composition_values=(0., .5, 1.), tv_values=(1000., 2000.)):
    """Add explicit Ar/N2 composition and a coupled N2(X,v) Tv ladder."""
    from rmgpy.tools.eedf.vibrational import ladder_coupling
    spec = copy.deepcopy(base_spec)
    if (not isinstance(composition_axis, str) or not re.fullmatch(r'[A-Za-z_][A-Za-z_0-9]*', composition_axis)
            or not isinstance(tv_axis, str) or not re.fullmatch(r'[A-Za-z_][A-Za-z_0-9]*', tv_axis)
            or composition_axis == tv_axis or 'u' in (composition_axis, tv_axis)):
        raise SpecError('invalid mixture axis name')
    if (not composition_values or not tv_values or
            any(type(v) not in (int, float) for v in (*composition_values, *tv_values))):
        raise SpecError('mixture axes require numeric points')
    if (min(composition_values) < 0 or max(composition_values) > 1
            or sorted(set(composition_values)) != list(composition_values)
            or min(tv_values) <= 0 or sorted(set(tv_values)) != list(tv_values)):
        raise SpecError('mixture composition and Tv axes must be ordered physical values')
    ladder = ladder_coupling(ladder_manifest, bindings, axis=tv_axis,
                             parent='N2(X)', d6_policy=_diagnostic_d6())
    domain = ladder.get('source_domain')
    if domain != {'initial_levels': list(range(30)), 'final_v': 58}:
        raise SpecError('N2 ladder must retain source domain initial v=0..29, final v=0..58')
    spec['axes'] = {**spec['axes'], composition_axis: list(composition_values), tv_axis: list(tv_values)}
    spec['gas_properties'] = copy.deepcopy(spec['gas_properties'])
    spec['gas_properties']['fraction'] = [f'Ar = 1-{{{composition_axis}}}',
                                          f'N2 = {{{composition_axis}}}']
    populations = ['Ar(1S0) = 1', 'N2(X) = 1']
    zero_population_tail = {}
    for level in ladder['levels']:
        expression = (f'boltzmann:{tv_axis}:{level["v"]}'
                      if level['v'] <= ladder['K'] else '0')
        if level['v'] <= ladder['K']:
            populations.append(f'{level["loki_state"]} = {expression}')
        else:
            zero_population_tail[level['loki_state']] = 0.0
    ladder['zero_population_tail'] = zero_population_tail
    spec['state_properties'] = copy.deepcopy(spec['state_properties'])
    spec['state_properties']['population'] = populations
    child_prefix = ladder['parent'][:-1] + ',v='
    for name, field in (('energy', 'energy_eV'), ('statisticalWeight', 'statistical_weight')):
        # The base spec can carry a complete source ladder for accounting.
        # Remove those declarations before adding only solver-supported rows.
        entries = [entry for entry in spec['state_properties'].get(name, [])
                   if not entry.split(' = ', 1)[0].startswith(child_prefix)]
        entries.extend(f'{level["loki_state"]} = {level[field]}'
                       for level in ladder['levels'] if level['v'] <= ladder['K'])
        spec['state_properties'][name] = entries
    spec['vibrational_ladders'] = [ladder]
    return spec


def prepare_mixture_spec(base_spec, gas_maps, database, collision_set,
                         ladder_source_manifest, *, cascade_yields=None,
                         yield_table=None, radiative_source=None,
                         yield_table_sha256=None, source_sha256=None):
    """Prepare a multi-gas conversion from explicit (gas, map, hash, library) inputs."""
    spec = load_spec(base_spec)
    if not isinstance(gas_maps, (list, tuple)) or len(gas_maps) < 2:
        raise SpecError('mixture conversion requires at least two gas maps')
    frozen_by_gas, libraries_needed, map_hashes, library_hashes = {}, set(), {}, {}
    library_by_gas = {}
    for item in gas_maps:
        if not isinstance(item, dict) or set(item) != {'gas', 'path', 'sha256', 'library'}:
            raise SpecError('each gas map requires gas, path, sha256 and library')
        gas, label = item['gas'], item['library']
        if (not isinstance(gas, str) or not gas or gas in frozen_by_gas
                or not isinstance(label, str) or not label):
            raise SpecError('invalid or duplicate gas map identity')
        data = _read_pinned_bytes(item['path'], item['sha256'])
        try:
            frozen_by_gas[gas] = yaml.load(data, Loader=_UniqueLoader)
        except (yaml.YAMLError, UnicodeError) as error:
            raise SpecError('malformed frozen channel map YAML') from error
        libraries_needed.add(label)
        library_by_gas[gas] = label
        map_hashes[gas] = item['sha256']
    provenance = {'frozen_map_sha256': map_hashes,
                  'map_sets_by_gas': {gas: frozen.get('set')
                                      for gas, frozen in frozen_by_gas.items()},
                  'base_channel_map': copy.deepcopy(spec['channel_map']),
                  'input_files': copy.deepcopy(spec['input_files']),
                  'collision_set': collision_set}
    with tempfile.TemporaryDirectory(prefix='eedf-map-') as directory:
        scratch = Path(directory)
        snapshot = copy.deepcopy(spec)
        for i, entry in enumerate(snapshot['input_files'].values()):
            data = _read_pinned_bytes(entry['path'], entry['sha256'])
            destination = scratch / (str(i) + Path(entry['path']).suffix)
            destination.write_bytes(data)
            entry['path'] = str(destination)
        data = _read_pinned_bytes(spec['channel_map']['path'], spec['channel_map']['sha256'])
        (scratch / 'channels.json').write_bytes(data)
        snapshot['channel_map']['path'] = str(scratch / 'channels.json')
        loaded = {}
        for label in sorted(libraries_needed):
            loaded[label], library_hashes[label] = _load_library(database, scratch, label)
        cascade_provenance = None
        if cascade_yields is None:
            cascade_yields, cascade_provenance = _cascade_snapshot(
                database, yield_table=yield_table, radiative_source=radiative_source,
                yield_table_sha256=yield_table_sha256, source_sha256=source_sha256)
        else:
            cascade_provenance = {name: cascade_yields[name] for name in
                                  ('yield_table_sha256', 'source_sha256', 'radiative_model')}
        if cascade_provenance is not None:
            provenance['cascade_sources'] = cascade_provenance
        try:
            converted, report = convert_channel_maps(snapshot, frozen_by_gas, loaded,
                                                      collision_set, cascade_yields, library_by_gas)
        except MissingOwners as error:
            error.report['provenance'] = provenance
            raise
    provenance['library_labels'] = sorted(libraries_needed)
    provenance['library_sha256'] = library_hashes
    report['provenance'] = provenance
    from rmgpy.tools.eedf.n2_vibrational import generate
    generator_manifest = copy.deepcopy(ladder_source_manifest)
    generator_manifest['D6'] = _diagnostic_d6()
    generated = generate(generator_manifest, K=10)
    ladder_manifest = json.loads(generated['manifest.json'])
    bindings = [{'source_sha256': row['source']['sha256'],
                 'source_block': row['source_block'],
                 'description': (f"e + N2(X,v={row['initial_v']}) <-> "
                                 f"e + N2(X,v={row['final_v']}), Vibrational")}
                for row in ladder_manifest['transition_audit']]
    ladder_spec = build_mixture_spec(spec, ladder_manifest, bindings)
    return spec, converted, report, ladder_spec


def _write_outputs(destination, report, spec=None, channels=None, mixture_spec=None):
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
        if mixture_spec is not None:
            mixture_spec = copy.deepcopy(mixture_spec)
            mixture_spec['channel_map'] = copy.deepcopy(spec['channel_map'])
            mixture_spec['output_root'] = spec['output_root']
            mixture_spec['scratch_root'] = spec['scratch_root']
            files['ar-n2-spec.json'] = mixture_spec
        else:
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
    parser.add_argument('--base-spec', required=True)
    for option in ('frozen-map', 'frozen-map-sha256'):
        parser.add_argument('--' + option)
    parser.add_argument('--gas-map', nargs=4, action='append', metavar=('GAS', 'PATH', 'SHA256', 'LIBRARY'))
    parser.add_argument('--ladder-manifest')
    for option in ('database', 'collision-set', 'output-directory'):
        parser.add_argument('--' + option, required=True)
    for option in ('yield-table', 'radiative-source', 'yield-table-sha256', 'source-sha256'):
        parser.add_argument('--' + option)
    args = parser.parse_args(argv)
    try:
        optional = {name: getattr(args, name) for name in
                    ('yield_table', 'radiative_source', 'yield_table_sha256', 'source_sha256')
                    if getattr(args, name) is not None}
        mixture_spec = None
        if args.gas_map:
            if (not args.ladder_manifest
                    or args.frozen_map or args.frozen_map_sha256):
                raise SpecError('multi-gas conversion requires a ladder source manifest and no legacy frozen-map options')
            with open(args.ladder_manifest) as stream:
                ladder_manifest = json.load(stream)
            gas_maps = [{'gas': gas, 'path': path, 'sha256': digest, 'library': library}
                        for gas, path, digest, library in args.gas_map]
            spec, channels, report, mixture_spec = prepare_mixture_spec(
                args.base_spec, gas_maps, args.database, args.collision_set,
                ladder_manifest, **optional)
        else:
            if not args.frozen_map or not args.frozen_map_sha256:
                raise SpecError('legacy conversion requires --frozen-map and --frozen-map-sha256')
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
        _write_outputs(args.output_directory, report, spec, channels, mixture_spec)
    except OSError as error:
        print('cannot publish conversion: ' + str(error), file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2))
    return 2 if report['outcome'] == 'refused' else 0


if __name__ == '__main__':
    sys.exit(_main())
