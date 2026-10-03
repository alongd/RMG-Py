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

"""Validate collision metadata against the pinned solver's physical inputs."""
from functools import lru_cache
import json
from pathlib import Path
import re

import numpy as np
from scipy.constants import electron_mass

from rmgpy.tools.eedf.integrity import physical_properties, resolved_properties
from rmgpy.tools.eedf.schema import SpecError


def validate_map(channels):
    """Check internal bounds and structural production classifications."""
    required = {'description', 'kind', 'classification', 'threshold_eV',
                'target_fraction', 'product_fraction', 'sigma_max_m2',
                'flux_group', 'reaction', 'cross_section'}
    optional = {'mass_ratio', 'opb_eV', 'statistical_weight_ratio'}
    if not channels:
        raise SpecError('empty channel map')
    for channel in channels:
        if (set(channel) - required - optional or required - set(channel)
                or channel['classification'] not in ('A', 'B', 'C', 'D')):
            raise SpecError('channel map fields/classification')
        if channel['kind'] not in ('elastic', 'excitation', 'vibrational', 'rotational', 'ionization', 'attachment'):
            raise SpecError('unsupported channel power kind')
        for name in ('threshold_eV', 'sigma_max_m2') + tuple(optional & set(channel)):
            value = channel[name]
            if type(value) not in (int, float) or not np.isfinite(value) or value < 0:
                raise SpecError('invalid channel ' + name)
        section = channel['cross_section']
        if set(section) != {'energy_eV', 'sigma_m2'}:
            raise SpecError('channel cross section fields')
        energy = np.asarray(section['energy_eV'], dtype=float)
        sigma = np.asarray(section['sigma_m2'], dtype=float)
        if (energy.ndim != 1 or len(energy) < 2 or sigma.shape != energy.shape
                or not np.all(np.isfinite(energy)) or not np.all(np.isfinite(sigma))
                or np.any(energy < 0) or np.any(np.diff(energy) <= 0) or np.any(sigma < 0)):
            raise SpecError('invalid channel cross section')
        if float(channel['sigma_max_m2']).hex() != float(np.max(sigma)).hex():
            raise SpecError('sigma_max_m2 differs from cross section')
        if channel['kind'] == 'elastic' and channel.get('mass_ratio', 0) <= 0:
            raise SpecError('elastic mass ratio')
        reaction = channel['reaction']
        if ((channel['classification'] == 'A') != (reaction is not None)
                or reaction is not None and set(reaction) != {'library', 'index', 'repr'}):
            raise SpecError('classification disagrees with mapped reaction identity')
        if not isinstance(channel['flux_group'], str) or not channel['flux_group']:
            raise SpecError('channel flux_group')


@lru_cache(maxsize=64)
def _classic_sections(path, digest):
    """Pinned EedfCollisions.cpp:934–1014, LookupTable.cpp:83–121 subset.

    The PARAM line and following bracketed process, rather than the LXCat
    title or PROCESS label, supply the actual collision and threshold.
    Unsupported stoichiometry/formats refuse rather than guess identities.
    Cache keys include the content digest checked by the caller.
    """
    if Path(path).suffix == '.json':
        raise SpecError('JSON collision metadata validation is unsupported')
    lines = Path(path).read_text().splitlines()
    result = []
    i = 0
    while i < len(lines):
        line = lines[i]
        i += 1
        if 'PARAM.:' not in line:
            continue
        threshold = re.search(r'E = +(\S*) +eV', line)
        threshold = float(threshold[1]) if threshold else 0.
        process = re.search(r'\[(.+?)(<->|->)(.+?), (\w+)\]', lines[i])
        i += 1
        if not process:
            raise SpecError('unsupported LXCat process metadata')
        lhs, arrow, rhs, kind = process.groups()
        left = re.split(r'\s+\+\s+', lhs.strip())
        right = re.split(r'\s+\+\s+', rhs.strip())
        heavy = [s for s in left if s != 'e']
        products = [s for s in right if s != 'e']
        if (len(left) != 2 or left.count('e') != 1 or len(heavy) != 1
                or not all(re.fullmatch(r'[A-Z][A-Za-z0-9]*\([^()]+\)', s) for s in heavy + products)):
            raise SpecError('unsupported LXCat collision identity')
        expected_electrons = 2 if kind == 'Ionization' else (0 if kind == 'Attachment' else 1)
        if right.count('e') != expected_electrons:
            raise SpecError('unsupported LXCat product stoichiometry')
        description = 'e + ' + heavy[0] + ' ' + arrow + ' '
        description += 'e + ' * expected_electrons + ' + '.join(products) + ', ' + kind
        while i < len(lines) and not lines[i].startswith('--'):
            i += 1
        i += 1
        values = []
        while i < len(lines) and not lines[i].startswith('--'):
            if not lines[i].startswith('#'):
                fields = lines[i].split()
                if len(fields) != 2:
                    raise SpecError('ambiguous LXCat data row')
                values.append([float(v) for v in fields])
            i += 1
        data = np.asarray(values)
        if data.ndim != 2 or data.shape[1] != 2:
            raise SpecError('missing LXCat cross section')
        result.append({'description': description, 'kind': kind.lower(),
                       'threshold_eV': threshold,
                       'cross_section': {'energy_eV': data[:, 0].tolist(), 'sigma_m2': data[:, 1].tolist()}})
    if not result:
        raise SpecError('missing LXCat collisions')
    return result


@lru_cache(maxsize=64)
def _gas_values(path, digest):
    """Read the supported numeric gas-property file grammar."""
    if Path(path).suffix == '.json':
        data = json.loads(Path(path).read_text())
    else:
        data = {}
        for line in Path(path).read_text().splitlines():
            fields = line.split('%', 1)[0].split()
            if fields:
                if len(fields) != 2 or fields[0] in data:
                    raise SpecError('ambiguous gas property file')
                data[fields[0]] = float(fields[1])
    if not isinstance(data, dict) or any(type(v) not in (int, float) or not np.isfinite(v) for v in data.values()):
        raise SpecError('unresolved gas property')
    return data


def validate_physical_map(channels, spec, coordinates):
    """Refuse verdict-changing metadata that differs from actual solver inputs."""
    validate_map(channels)
    gas, _ = physical_properties(spec)
    sections = {}
    for entry in spec['input_files'].values():
        if entry['kind'] == 'cross_section':
            for section in _classic_sections(entry['path'], entry['sha256']):
                name = section['description']
                if name in sections:
                    raise SpecError('duplicate physical collision')
                sections[name] = section
    if set(sections) != {c['description'] for c in channels} or len(sections) != len(channels):
        raise SpecError('channel map does not cover physical collision identities')
    _, states, _, _ = resolved_properties(spec, coordinates)
    weights = {v.split(' = ', 1)[0]: float(v.split(' = ', 1)[1])
               for v in states.get('statisticalWeight', [])}

    def gas_value(name, species, default=None):
        if name not in gas:
            if default is not None:
                return default
            raise SpecError('missing physical gas ' + name)
        entry = spec['input_files'][gas[name]]
        values = _gas_values(entry['path'], entry['sha256'])
        if species not in values and default is None:
            raise SpecError('missing physical gas ' + name + ' for ' + species)
        return values.get(species, default)

    for channel in channels:
        physical = sections[channel['description']]
        for field in ('kind', 'threshold_eV', 'cross_section'):
            if channel[field] != physical[field]:
                raise SpecError('channel ' + field + ' differs from hashed LoKI input')
        target = channel['description'].split(' ', 3)[2]
        species = target.split('(', 1)[0]
        if channel['kind'] == 'elastic':
            mass = gas_value('mass', species)
            if mass <= 0 or channel['mass_ratio'] != electron_mass / mass:
                raise SpecError('mass_ratio differs from hashed LoKI mass')
        if channel['kind'] == 'ionization':
            opb = gas_value('OPBParameter', species, -1.)
            if opb < 0:
                opb = channel['threshold_eV']
            if channel.get('opb_eV', channel['threshold_eV']) != opb:
                raise SpecError('opb_eV differs from hashed LoKI OPBParameter')
        if '<->' in channel['description']:
            product = channel['description'].split(' <-> e + ', 1)[-1].rsplit(', ', 1)[0]
            if target not in weights or product not in weights or min(weights[target], weights[product]) <= 0:
                raise SpecError('reverse statistical weights must be explicit positive setup constants')
            if channel.get('statistical_weight_ratio') != weights[target] / weights[product]:
                raise SpecError('statistical_weight_ratio differs from LoKI setup')
