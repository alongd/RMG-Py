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

"""Solver build and physical-file closure for the pinned legacy setup."""
import copy
import json
import math
from decimal import Decimal
from pathlib import Path
import re
import subprocess

from rmgpy.tools.eedf.schema import FingerprintMismatch, SpecError, file_hash


def solver_environment(spec):
    """Use loader resolution independent of the invoking shell."""
    # No parent-shell values reach ldd or the solver, including future loader,
    # allocator, locale, OpenMP, and vendor-library switches.
    return {'PATH': '/usr/bin:/bin', 'LANG': 'C', 'LC_ALL': 'C', 'TZ': 'UTC',
            'OMP_NUM_THREADS': str(spec['omp_threads']), 'OMP_DYNAMIC': 'FALSE'}



def check_solver(spec):
    """Pin every shared object resolved from the immutable build by ldd."""
    if file_hash(spec['binary']['path']) != spec['binary']['sha256']:
        raise FingerprintMismatch('binary SHA-256')
    build = Path(spec['cmake_cache']['path']).resolve().parent
    pins = {str(Path(p).resolve()): digest for p, digest in spec.get('shared_objects', {}).items()}
    for path, digest in pins.items():
        if file_hash(path) != digest:
            raise FingerprintMismatch('shared object SHA-256: ' + path)
    completed = subprocess.run(['ldd', spec['binary']['path']], env=solver_environment(spec),
                               capture_output=True, text=True, check=False)
    if completed.returncode or 'not found' in completed.stdout:
        raise FingerprintMismatch('shared object resolution')
    resolved = set()
    for name in re.findall(r'(?:=>\s+|^\s*)(/\S+)', completed.stdout, re.M):
        path = Path(name).resolve()
        if build in path.parents:
            resolved.add(str(path))
    if not resolved or resolved != set(pins):
        raise FingerprintMismatch('shared objects from build do not match pins')
    return pins


def physical_properties(spec):
    """Close every property-file reference; return scratch-safe resolved values.

    Population files may include other files, but every leaf must be numeric.
    Symbolic population functions are refused rather than stored as assumptions.
    """
    registry = spec['input_files']
    by_path = {str(Path(entry['path']).resolve()): name for name, entry in registry.items()}
    for name, entry in registry.items():
        if legacy_string(name) != name:
            raise SpecError('ambiguous scratch input name: ' + name)
        if (Path(name).is_absolute() or '..' in Path(name).parts
                or any(c in name for c in ('%', '\n', '\r'))):
            raise SpecError('unsafe input_files name: ' + name)
        if file_hash(entry['path']) != entry['sha256']:
            raise FingerprintMismatch('input file ' + name)

    def reference(value):
        # OffSideToJSON strips % comments, then JSON-decodes quoted strings.
        value = legacy_string(value)
        if value in registry:
            return value
        if Path(value).is_absolute() and str(Path(value).resolve()) in by_path:
            return by_path[str(Path(value).resolve())]
        raise SpecError('unhashed physical file: ' + str(value))

    def entries(values, population=False, seen=()):
        if not isinstance(values, list):
            raise SpecError('legacy state properties must be lists')
        result = []
        for raw in values:
            value = legacy_string(raw)
            # LegacyToJSON.cpp:94 uses regex_search, NOT a test for '='.
            match = re.search(r'(\S+)\s+=\s+(\S+)\s*$', value)
            if match:
                if match.start() != 0:
                    raise SpecError('ambiguous state property assignment: ' + value)
                key, expr = match.groups()
                if not re.fullmatch(r'\{[A-Za-z_][A-Za-z_0-9]*\}|[+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+\-]?\d+)?', expr):
                    raise SpecError('symbolic state property must be resolved to numbers')
                result.append(neutral_state_key(key) + ' = ' + expr)
                continue
            name = reference(value)
            if name in seen:
                raise SpecError('recursive property file')
            path = Path(registry[name]['path'])
            if path.suffix == '.json':
                data = json.loads(path.read_text())
                if not isinstance(data, dict) or set(data) - {'states', 'files'}:
                    raise SpecError('unsupported JSON state property file')
                expanded = []
                for key, prop in data.get('states', {}).items():
                    if (not isinstance(prop, dict) or set(prop) != {'type', 'value'}
                            or prop['type'] != 'constant' or type(prop['value']) not in (int, float)):
                        raise SpecError('symbolic JSON state property must be resolved to numbers')
                    expanded.append(key + ' = ' + str(prop['value']))
                files = data.get('files', [])
                if not isinstance(files, list):
                    raise SpecError('JSON state property files must be a list')
                # The pinned loader applies inline states before file entries.
                if any(not isinstance(v, str) or '%' in v or '\n' in v or '\r' in v for v in files):
                    raise SpecError('ambiguous JSON file reference')
                expanded.extend(json.dumps(v) for v in files)
                result.extend(entries(expanded, population, seen + (name,)))
            else:
                expanded = []
                for raw_line in path.read_text().splitlines():
                    line = raw_line.split('%', 1)[0].strip()
                    if not line:
                        continue
                    # readLegacyStatePropertyFile searches for two trailing
                    # tokens first. Refuse any ignored prefix rather than
                    # resolving different bytes than the solver.
                    tokens = line.split()
                    if len(tokens) == 2:
                        expanded.append(tokens[0] + ' = ' + tokens[1])
                    elif len(tokens) == 1:
                        expanded.append(json.dumps(tokens[0]))
                    else:
                        raise SpecError('ambiguous legacy state property file line: ' + line)
                result.extend(entries(expanded, population, seen + (name,)))
        # Flatten to explicit constants.  LoKI resolves selector aliases to a
        # state before applying properties, so exact duplicate resolved values
        # are harmless while conflicting values would have precedence-dependent
        # behavior.  Keep one canonical entry and refuse only the latter.
        flattened = {}
        for entry in result:
            key, expr = entry.split(' = ', 1)
            value = (float(expr.replace('D', 'e').replace('d', 'e'))
                     if re.fullmatch(r'[+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+\-]?\d+)?', expr) else expr)
            if key in flattened:
                if flattened[key][0] != value:
                    raise SpecError('conflicting state property selector')
                continue
            flattened[key] = (value, expr)
        return [key + ' = ' + expr for key, (_, expr) in flattened.items()]

    gas = copy.deepcopy(spec['gas_properties'])
    for name, value in gas.items():
        if name == 'fraction':
            if not isinstance(value, list):
                raise SpecError('gas fractions must be a list')
        elif isinstance(value, str):
            gas[name] = reference(value)
        else:
            raise SpecError('unsupported gas property shape: ' + name)
    states = {name: entries(value, name == 'population') for name, value in spec['state_properties'].items()}
    if set(states) - {'energy', 'statisticalWeight', 'population'}:
        raise SpecError('unsupported state property key')
    allowed_options = {'isOn', 'eedfType', 'ionizationOperatorType', 'growthModelType',
                       'includeEECollisions', 'numerics', 'LXCatFiles',
                       'LXCatFilesExtra', 'effectiveCrossSectionPopulations'}
    if set(spec['solver_options']) - allowed_options:
        raise SpecError('unsupported solver option; file closure is not proven')
    for key in ('LXCatFiles', 'LXCatFilesExtra', 'effectiveCrossSectionPopulations'):
        if key in spec['solver_options']:
            values = spec['solver_options'][key]
            if not isinstance(values, list):
                raise SpecError(key + ' must be a list of hashed files')
            for filename in values:
                reference(filename)
            # These change the collision set/densities outside the stored map.
            # Refuse until a corresponding supported numeric identity exists.
            if key != 'LXCatFiles':
                raise SpecError(key + ' is unsupported for qualified tables')
    return gas, states


def legacy_string(value):
    """Decode a string as the pinned off-side setup parser does."""
    if not isinstance(value, str) or '\n' in value or '\r' in value:
        raise SpecError('legacy entry must be a single string')
    text = value.split('%', 1)[0].strip()
    try:
        decoded = json.loads(text)
    except ValueError:
        decoded = text
    if not isinstance(decoded, str) or not decoded or '\n' in decoded or '\r' in decoded:
        raise SpecError('ambiguous legacy string: ' + value)
    return decoded


def neutral_state_key(selector):
    """Resolve the supported subset of LoKI's propertyStateFromString.

    StateEntry.cpp:393-495 splits comma tokens with std::getline (which drops
    one trailing empty token), then consumes an optional charge before the
    electronic identifier. Gas::findState uses these fields, not the original
    spelling. Root, charge, wildcard and hierarchical selectors need a state
    tree and are deliberately refused for every state property.
    """
    outer = re.fullmatch(r'([A-Za-z][A-Za-z0-9]*)(\(.*\))?', selector)
    if not outer or not outer[2]:
        raise SpecError('state selector needs an explicit neutral electronic state: ' + selector)
    body = outer[2][1:-1]
    tokens = body.split(',')
    if body.endswith(','):
        tokens.pop()  # Match C++ std::getline, including repeated delimiters.
    charge, electronic = '', ''
    electronic_index = 0
    if tokens and re.fullmatch(r'\s*([-+]*)\s*', tokens[0]):
        charge = tokens[0].strip()
        electronic_index = 1
    if len(tokens) > electronic_index:
        electronic = tokens[electronic_index].strip()
    if (charge or not re.fullmatch(r'[^,()=\s*]+', electronic)
            or len(tokens) != electronic_index + 1):
        raise SpecError('unsupported state selector: ' + selector)
    return outer[1] + '(' + electronic + ')'


def resolved_properties(spec, coordinates):
    """Produce the very same explicit constants for setup and row metadata.

    Flattened numeric state selectors are deliberately limited to neutral
    electronic states. Wildcards, charge overrides and hierarchical selectors
    need LoKI's full state tree and are refused rather than guessed.
    """
    state = {name: item['reference'] for name, item in spec['envelopes'].items()}
    state.update(coordinates)
    gas, properties = physical_properties(spec)

    def numbers(entries, population=False):
        result = {}
        for entry in entries:
            match = re.fullmatch(r'(\S+)\s+=\s+(\S+)', legacy_string(entry))
            if not match:
                raise SpecError('ambiguous numeric property assignment')
            key, expr = match.groups()
            if population:
                key = neutral_state_key(key)
            if not population and not re.fullmatch(r'[A-Za-z][A-Za-z0-9]*', key):
                raise SpecError('ambiguous gas fraction identifier')
            try:
                value = float(expr.format_map(state).replace('D', 'e').replace('d', 'e'))
            except (ValueError, KeyError) as error:
                raise SpecError('symbolic population/fraction is unresolved') from error
            if (not math.isfinite(value) or value < 0 or value > 1 or key in result
                    or (value == 0 and math.copysign(1., value) < 0)):
                raise SpecError('invalid/duplicate population or fraction')
            result[key] = value
        return result

    fractions = numbers(gas['fraction'])
    populations = numbers(properties.get('population', []), population=True)
    # patchFractions does not understand signs or exponents. Exact decimal
    # spelling round-trips even tiny fractions without changing their bits.
    gas['fraction'] = [key + ' = ' + format(Decimal.from_float(value), 'f')
                       for key, value in fractions.items()]
    properties['population'] = [key + ' = ' + format(value, '.17g')
                                for key, value in populations.items()]
    for name in properties.keys() - {'population'}:
        expanded = []
        for entry in properties[name]:
            key, expr = entry.split(' = ', 1)
            try:
                value = float(expr.format_map(state).replace('D', 'e').replace('d', 'e'))
            except (ValueError, KeyError) as error:
                raise SpecError('unresolved state property') from error
            if not math.isfinite(value):
                raise SpecError('nonfinite state property')
            expanded.append(key + ' = ' + format(value, '.17g'))
        properties[name] = expanded
    return gas, properties, fractions, populations


def channel_fractions(channel_map, gas, populations, coordinates):
    """Check map fractions bit for bit against LoKI's neutral state densities.

    Gas.cpp:233-238 multiplies the electronic population by the neutral
    parent's density. The neutral/root factors default to one for nonzero gas
    fractions (Gas.cpp:275-286). Unspecified electronic states default to zero.
    Only unambiguous unit-stoichiometry channels are supported here.
    """
    targets, products = [], []

    def density(identifier):
        match = re.fullmatch(r'([A-Z][A-Za-z0-9]*)\([^,()=\s*]+\)', identifier)
        if (not match or match[1] not in gas
                or re.fullmatch(r'[-+]+', identifier.split('(', 1)[-1].rstrip(')'))):
            raise SpecError('channel population cannot be resolved: ' + identifier)
        return gas[match[1]] * populations.get(identifier, 0.)

    for channel in channel_map:
        match = re.fullmatch(r'e \+ (.+?) (?:<->|->) (.+), [A-Za-z]+', channel['description'])
        if not match:
            raise SpecError('channel population cannot be resolved from description')
        target = density(match[1])
        if channel['kind'] in ('excitation', 'vibrational', 'rotational'):
            product_match = re.fullmatch(r'e \+ (.+)', match[2])
            if not product_match:
                raise SpecError('superelastic population cannot be resolved')
            product = density(product_match[1])
        else:
            product = 0.  # No superelastic process for these collision types.
        for name, value in (('target_fraction', target), ('product_fraction', product)):
            try:
                supplied = float(str(channel[name]).format_map(coordinates))
            except (ValueError, KeyError) as error:
                raise SpecError('unresolved channel population') from error
            if supplied.hex() != value.hex():
                raise SpecError('channel ' + name + ' differs from LoKI population: ' + channel['description'])
        targets.append(target)
        products.append(product)
    return targets, products
