"""Refresh an Ar or Ar/N2 generation spec for resolved argon electron owners.

The collision realization is produced by the database collector. This tool
preserves the supplied axes and non-argon data; it does not generate a table.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import tempfile

import yaml

from rmgpy.tools.eedf.channel_map import (
    _UniqueLoader, _cascade_snapshot, _derived_elastic, _load_library, _read_pinned_bytes,
    _write_outputs, convert_channel_maps,
)
from rmgpy.tools.eedf.schema import FingerprintMismatch, SpecError, load_spec


def refresh_spec(base_spec, database, realization, destination, *, max_parallel=5):
    """Publish a source-checked owner conversion into a new output directory."""
    base_spec, database, realization = map(lambda p: Path(p).resolve(), (base_spec, database, realization))
    raw = base_spec.read_bytes()
    original = json.loads(raw) if base_spec.suffix == '.json' else yaml.load(raw, Loader=_UniqueLoader)
    spec = copy.deepcopy(original)
    if type(max_parallel) is not int or not 1 <= max_parallel <= 5 or type(spec.get('max_parallel')) is not int:
        raise SpecError('max_parallel must be an integer in [1, 5]')
    spec['max_parallel'] = min(spec['max_parallel'], max_parallel)
    # Resolve against the original spec before moving its owned snapshot.
    for key in ('binary', 'cmake_cache', 'channel_map'):
        spec[key]['path'] = str((base_spec.parent / spec[key]['path']).resolve())
    for entry in spec['input_files'].values():
        entry['path'] = str((base_spec.parent / entry['path']).resolve())
    spec['shared_objects'] = {str((base_spec.parent / name).resolve()): digest
                              for name, digest in spec['shared_objects'].items()}
    gases = spec['arm']['gases']
    if set(gases) not in ({'Ar'}, {'Ar', 'N2'}):
        raise SpecError('argon refresh supports the supplied Ar or Ar/N2 spec')
    maps, map_hashes, libraries, library_hashes = {}, {}, {}, {}
    names = {'Ar': ('argon_channel_map.yaml', 'PlasmaArgon'),
             'N2': ('nitrogen_channel_map.yaml', 'PlasmaNitrogen')}
    for gas in gases:
        content = (database / 'input/eedf' / names[gas][0]).read_bytes()
        maps[gas] = yaml.load(content, Loader=_UniqueLoader)
        map_hashes[gas] = hashlib.sha256(content).hexdigest()
    declared = maps['Ar'].get('collision_realization')
    actual = json.loads((realization / 'realization.json').read_bytes())
    if declared != actual or actual.get('schema') != 'argon-reversible-realization/1':
        raise FingerprintMismatch('argon realization differs from frozen map')
    for filename, pins in actual['files'].items():
        path = realization / filename
        _read_pinned_bytes(path, pins['sha256'])
        matches = [key for key in spec['input_files'] if Path(key).name == filename]
        if len(matches) > 1:
            raise SpecError('ambiguous argon collision filename: ' + filename)
        key = matches[0] if matches else 'Argon/' + filename
        if matches and spec['input_files'][key]['sha256'] not in (pins['base_sha256'], pins['sha256']):
            raise FingerprintMismatch('base spec argon collision source: ' + filename)
        spec['input_files'][key] = {'path': str(path), 'sha256': pins['sha256'], 'kind': 'cross_section'}
    for key in ('excited_fraction_max', 'argon_elimination_electron_density_max_m3'):
        spec[key] = actual[key]
    states = maps['Ar']['states']
    resolved = ('Ar(1S0)', 'Ar(3P2)', 'Ar(3P1)', 'Ar(3P0)', 'Ar(1P1)')
    properties = spec['state_properties']
    for key in ('population', 'energy', 'statisticalWeight'):
        existing = properties.get(key, [])
        if not isinstance(existing, list):
            raise SpecError('argon refresh requires inline state properties')
        properties[key] = [line for line in existing if not line.startswith('Ar(')]
    for state in resolved:
        properties['population'].append(f'{state} = {1 if state == "Ar(1S0)" else 0}')
        properties['energy'].append(f'{state} = {states[state]["excitation_eV"]}')
        properties['statisticalWeight'].append(f'{state} = {states[state]["g"]}')
    if 'engine_state_map' in spec:
        present = {row['loki_state'] for row in spec['engine_state_map']}
        for state in resolved[1:]:
            if state not in present:
                term = state[3:-1]
                spec['engine_state_map'].append(dict(formula='Ar', electronic_state=term,
                    vibrational_level=-1, multiplicity=int(term[0]), loki_state=state,
                    statistical_weight=states[state]['g']))
    channels = json.loads(_read_pinned_bytes(spec['channel_map']['path'], spec['channel_map']['sha256']))
    # Rebuild the argon physical rows from their newly pinned collision files.
    derived_elastic = {_derived_elastic(row['description']) for row in channels if row['kind'] == 'effective'}
    retained = [row for row in channels if not row['description'].startswith('e + Ar(')
                and row['description'] not in derived_elastic]
    with tempfile.TemporaryDirectory(prefix='argon-4s-spec-') as directory:
        scratch = Path(directory)
        channel_bytes = (json.dumps(retained, indent=2) + '\n').encode()
        (scratch / 'channels.json').write_bytes(channel_bytes)
        spec['channel_map'] = {'path': str(scratch / 'channels.json'),
                               'sha256': hashlib.sha256(channel_bytes).hexdigest()}
        (scratch / 'spec.json').write_text(json.dumps(spec))
        spec = load_spec(scratch / 'spec.json')
        for gas in gases:
            label = names[gas][1]
            libraries[label], library_hashes[label] = _load_library(database, scratch, label)
        cascade, cascade_provenance = _cascade_snapshot(database)
        converted, report = convert_channel_maps(spec, maps, libraries, 'IST-Lisbon-Ar-2021', cascade,
                                                 library_by_gas={gas: names[gas][1] for gas in gases})
        report['provenance'] = {'base_spec_sha256': hashlib.sha256(raw).hexdigest(),
                                'frozen_map_sha256': map_hashes, 'library_sha256': library_hashes,
                                'cascade_sources': cascade_provenance, 'collision_realization': actual,
                                'max_parallel_before': original['max_parallel'],
                                'max_parallel_after': spec['max_parallel']}
        _write_outputs(destination, report, spec, converted, mixture_spec=spec if 'N2' in gases else None)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base-spec', required=True)
    parser.add_argument('--database', required=True, help='RMG-database repository root')
    parser.add_argument('--realization', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--max-parallel', type=int, default=5)
    args = parser.parse_args()
    report = refresh_spec(args.base_spec, args.database, args.realization, args.output,
                          max_parallel=args.max_parallel)
    print(json.dumps({'outcome': report['outcome'], 'selected_channels': report['selected_channel_count'],
                      'paired_owners': len(report['paired_owners']),
                      'branch_families': len(report['branch_families']),
                      'reviewed_inverse_exceptions': len(report['missing_owners'])}, indent=2))


if __name__ == '__main__':
    main()
