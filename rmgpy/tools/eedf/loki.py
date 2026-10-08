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

"""Pinned subprocess driver and parsers for LoKI-B C++ text output.

The external GPL solver is executed, never linked. No study build is changed.
"""
import copy
import io
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
from scipy.constants import elementary_charge, electron_mass

from rmgpy.tools.eedf.channels import _validate_physical_map
from rmgpy.tools.eedf.schema import EEDFError, FingerprintMismatch, _file_hash
from rmgpy.tools.eedf.integrity import (
    _channel_fractions,
    _check_solver,
    _physical_properties,
    _resolved_properties,
    _solver_environment,
)


class LoKIError(EEDFError):
    """Missing, non-converged, malformed, or inconsistent solver output."""


def _finite_scalar(value, name):
    """Return one finite Python float without retaining a caller array/view."""
    try:
        array = np.asarray(value)
        if array.ndim != 0:
            raise ValueError
        result = float(array)
    except (TypeError, ValueError, OverflowError) as error:
        raise LoKIError(name + ' must be a finite scalar') from error
    if not np.isfinite(result):
        raise LoKIError(name + ' must be a finite scalar')
    return result


SWARM = {
    'Reduced electric field': 'EN_Td',
    'Reduced diffusion coefficient': 'diffusion_N',
    'Reduced mobility coefficient': 'mobility_N',
    'Drift velocity': 'drift_velocity_m_s',
    'Reduced Townsend coefficient': 'townsend_N',
    'Reduced attachment coefficient': 'attachment_N',
    'Reduced energy diffusion coefficient': 'energy_diffusion_N',
    'Reduced energy mobility': 'energy_mobility_N',
    'Mean energy': 'mean_energy_eV',
    'Characteristic energy': 'characteristic_energy_eV',
    'Electron temperature': 'Te_LoKI_eV',
}
POWER = {
    'Field': 'field', 'Elastic collisions (gain)': 'elastic_gain',
    'Elastic collisions (loss)': 'elastic_loss', 'CAR (gain)': 'car_gain',
    'CAR (loss)': 'car_loss', 'Excitation inelastic collisions': 'excitation_loss',
    'Excitation superelastic collisions': 'excitation_gain',
    'Vibrational inelastic collisions': 'vibrational_loss',
    'Vibrational superelastic collisions': 'vibrational_gain',
    'Rotational inelastic collisions': 'rotational_loss',
    'Rotational superelastic collisions': 'rotational_gain',
    'Ionization collisions': 'ionization', 'Attachment collisions': 'attachment',
    'Electron density growth': 'growth', 'Power Balance': 'balance',
    'Relative Power Balance': 'relative_balance_percent',
}


def _terms(path, names):
    result = {}
    pattern = re.compile(r'^\s*(.*?)\s*=\s*(\S+)')
    for line in Path(path).read_text().splitlines():
        match = pattern.match(line)
        if match and match[1] in names and names[match[1]] not in result:
            result[names[match[1]]] = _number(match[2])
    if set(result) != set(names.values()) or not np.all(np.isfinite(list(result.values()))):
        raise LoKIError('missing/nonfinite terms in ' + str(path))
    return result


def _number(value):
    try:
        number = float(value.replace('D', 'e').replace('d', 'e'))
    except ValueError as error:
        raise LoKIError('invalid numeric output: ' + value) from error
    if not np.isfinite(number):
        raise LoKIError('nonfinite numeric output')
    return number


def parse_output(folder, normalization_atol, energy_edges=None):
    """Parse the four per-solve output files, preserving LoKI's SI/eV units."""
    folder, normalization_atol, energy_edges = copy.deepcopy(
        (folder, normalization_atol, energy_edges))
    normalization_atol = _finite_scalar(
        normalization_atol, 'normalization_atol')
    if normalization_atol < 0:
        raise LoKIError('normalization_atol must be a finite scalar')
    if energy_edges is not None:
        energy_edges = np.array(energy_edges, dtype=float, copy=True)
        energy_edges.setflags(write=False)
    folder = Path(folder)
    swarm = _terms(folder / 'swarmParameters.txt', SWARM)
    power = _terms(folder / 'powerBalance.txt', POWER)
    rates, channels = [], []
    lines = (folder / 'rateCoefficients.txt').read_text().splitlines()
    header = lines[0].split()
    expected = ['Ine.R.Coeff.(m3/s)', 'Sup.R.Coeff.(m3/s)', 'Description']
    if header != expected:
        raise LoKIError('unknown or swapped rate header')
    columns = {name: index for index, name in enumerate(header)}
    for line_index, line in enumerate(lines[1:], 1):
        if line.strip().startswith('---'):
            footer = [part.strip() for part in lines[line_index:] if part.strip()]
            if len(footer) != 3 or footer[1] != '* Extra Rate Coefficients *' or not all(re.fullmatch('-+', part) for part in (footer[0], footer[2])):
                raise LoKIError('unsupported extra rate coefficient section')
            break
        if not line.strip():
            continue
        fields = line.split(maxsplit=2)
        if len(fields) != 3:
            raise LoKIError('malformed rate coefficient row')
        a = _number(fields[columns['Ine.R.Coeff.(m3/s)']])
        b = _number(fields[columns['Sup.R.Coeff.(m3/s)']])
        if len(fields) != 3 or not np.isfinite(a + b) or min(a, b) < 0:
            raise LoKIError('invalid rate coefficient')
        rates.append([a, b])
        channels.append(fields[columns['Description']])
    if not rates or len(set(channels)) != len(channels):
        raise LoKIError('missing/duplicate rate channels')
    lines = (folder / 'eedf.txt').read_text().splitlines()
    header = lines[0].split()
    if header != ['Energy(eV)', 'EEDF(eV^-(3/2))', 'Anisotropy(eV^-(3/2))']:
        raise LoKIError('unknown or swapped EEDF header')
    try:
        eedf = np.asarray([[_number(v) for v in line.split()] for line in lines[1:] if line.strip()])
    except ValueError as error:
        raise LoKIError('malformed EEDF output') from error
    if eedf.ndim != 2 or eedf.shape[1] != len(header):
        raise LoKIError('malformed EEDF output')
    columns = {name: index for index, name in enumerate(header)}
    energy, f0 = eedf[:, columns['Energy(eV)']], eedf[:, columns['EEDF(eV^-(3/2))']]
    if len(energy) < 2 or not np.all(np.isfinite(eedf)) or np.any(energy < 0) or np.any(np.diff(energy) <= 0):
        raise LoKIError('invalid EEDF energy grid')
    if energy_edges is None:
        # With the pinned solver's arithmetic cell centres and lower edge 0,
        # each next face is 2*centre - previous face; midpoint spacing is wrong
        # on nonuniform grids. A configured grid supplies exact faces instead.
        edges = np.empty(len(energy) + 1)
        edges[0] = 0.
        for i, centre in enumerate(energy):
            edges[i + 1] = 2 * centre - edges[i]
    else:
        edges = np.asarray(energy_edges, dtype=float)
    if len(edges) != len(energy) + 1 or not np.all(np.isfinite(edges)) or np.any(np.diff(edges) <= 0):
        raise LoKIError('invalid EEDF edge grid')
    if not np.allclose((edges[:-1]+edges[1:])/2, energy, rtol=normalization_atol, atol=normalization_atol):
        raise LoKIError('EEDF centres disagree with edge grid')
    if edges[0] < -normalization_atol or np.any(f0 < -normalization_atol * np.max(abs(f0))):
        raise LoKIError('negative EEDF')
    f0 = np.maximum(f0, 0)
    norm = np.sum(f0 * np.sqrt(energy) * np.diff(edges))
    if not np.isfinite(norm) or abs(norm - 1) > normalization_atol:
        raise LoKIError('EEDF normalization: ' + str(norm))
    rates = np.asarray(rates)
    return {'EN_Td': swarm.pop('EN_Td'), 'swarm': swarm,
            'power_groups': power, 'channels': channels,
            'k_ine': rates[:, 0], 'k_sup': rates[:, 1],
            'energy_eV': energy, 'energy_edges_eV': edges, 'f0': f0,
            'converged': True, 'iteration_count': -1}


def _format(value, state):
    if isinstance(value, str):
        return value.format_map(state)
    return str(value).lower() if isinstance(value, bool) else str(value)


def _setup_text(spec, coordinates, fields_Td, folder):
    """Write a complete legacy setup from explicit, fingerprinted properties."""
    spec, coordinates, fields_Td = copy.deepcopy(
        (spec, coordinates, fields_Td))
    state = {name: item['reference'] for name, item in spec['envelopes'].items()}
    state.update(coordinates)
    wc = copy.deepcopy(spec['working_conditions'])
    wc.update(reducedField='[' + ','.join(format(x, '.17g') for x in fields_Td) + ']',
              gasTemperature=state.get('Tg_K', spec['Tg_K']), gasPressure=state.get('P_Pa', spec['P_Pa']))
    options = copy.deepcopy(spec['solver_options'])
    options['isOn'] = True
    options['LXCatFiles'] = [name for name, entry in spec['input_files'].items() if entry['kind'] == 'cross_section']
    options['gasProperties'], options['stateProperties'], _, _ = _resolved_properties(spec, coordinates)
    data = {'workingConditions': wc, 'electronKinetics': options,
            'gui': {'isOn': False, 'refreshFrequency': 1},
            'output': {'isOn': True, 'dataFormat': 'txt', 'folder': folder,
                       'dataFiles': ['log', 'eedf', 'swarmParameters', 'rateCoefficients', 'powerBalance', 'lookUpTable']}}
    lines = []

    def emit(values, level=0):
        for name, value in values.items():
            if not re.fullmatch(r'[A-Za-z][A-Za-z_0-9]*', name):
                raise LoKIError('ambiguous setup key')
            indent = '  ' * level
            if isinstance(value, dict):
                lines.append(indent + name + ':')
                emit(value, level + 1)
            elif isinstance(value, list):
                lines.append(indent + name + ':')
                for item in value:
                    text = _format(item, state)
                    if any(c in text for c in ('%', '\n', '\r')):
                        raise LoKIError('ambiguous setup list entry')
                    lines.append(indent + '  - ' + text)
            else:
                text = _format(value, state)
                if '\n' in text or '\r' in text or not re.fullmatch(r'[A-Za-z][A-Za-z_0-9]*', name):
                    raise LoKIError('ambiguous setup scalar or key')
                lines.append(indent + name + ': ' + text)
    emit(data)
    return '\n'.join(lines) + '\n'


class LoKIDriver:
    """Execute an immutable binary in isolated scratch jobs, with bounded OMP."""

    def __init__(self, spec):
        spec = copy.deepcopy(spec)
        self.spec = spec
        _check_solver(spec)
        _physical_properties(spec)
        for name, item in spec['input_files'].items():
            if _file_hash(item['path']) != item['sha256']:
                raise FingerprintMismatch('input file ' + name)
        for key in ('channel_map', 'cmake_cache'):
            if _file_hash(spec[key]['path']) != spec[key]['sha256']:
                raise FingerprintMismatch(key)
        self.setups = {}

    def run(self, coordinates, fields_Td, job_name):
        """One ordered scan; returns a row per requested field in that order."""
        spec, coordinates, fields_Td, job_name = copy.deepcopy(
            (self.spec, coordinates, fields_Td, job_name))
        if not isinstance(coordinates, dict):
            raise LoKIError('coordinates must be a named mapping')
        coordinates = {
            str(name): _finite_scalar(value, str(name))
            for name, value in coordinates.items()
        }
        try:
            fields_Td = tuple(
                _finite_scalar(value, 'fields_Td') for value in fields_Td)
        except TypeError as error:
            raise LoKIError('fields_Td must be a finite sequence') from error
        if not fields_Td or any(value <= 0. for value in fields_Td):
            raise LoKIError('fields_Td must contain positive finite scalars')
        if not isinstance(job_name, str):
            raise LoKIError('job_name must be a string')
        job_name = str(job_name)
        _check_solver(spec)
        _physical_properties(spec)
        work = Path(spec['scratch_root']) / job_name
        work.mkdir(parents=True, exist_ok=False)
        for name, item in spec['input_files'].items():
            target = work / name
            if Path(name).is_absolute() or '..' in Path(name).parts:
                raise LoKIError('unsafe scratch input name')
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(item['path'], target)
            if _file_hash(target) != item['sha256']:
                raise FingerprintMismatch(name)
        setup = _setup_text(spec, coordinates, fields_Td, 'solve')
        (work / 'setup.in').write_text(setup)
        self.setups[job_name] = setup
        env = _solver_environment(spec)
        with (work / 'stdout.log').open('w') as stdout, (work / 'stderr.log').open('w') as stderr:
            try:
                completed = subprocess.run([spec['binary']['path'], 'setup.in'], cwd=work,
                                           env=env, stdin=subprocess.DEVNULL, stdout=stdout,
                                           stderr=stderr, timeout=spec['timeout_s'], check=False)
            except subprocess.TimeoutExpired as error:
                raise LoKIError('timed out; inspect ' + str(work)) from error
        diagnostics = (work / 'stderr.log').read_text()
        if completed.returncode or re.search(r'did not converge|not converged|Results might be incorrect', diagnostics, re.I):
            raise LoKIError('LoKI refused; inspect ' + str(work))
        grid = spec['solver_options']['numerics']['energyGrid']
        edges = np.linspace(0, grid['maxEnergy'], grid['cellNumber'] + 1)
        rows = [parse_output(p, spec['tolerances']['normalization_atol'], edges)
                for p in (work / 'output' / 'solve').glob('reducedField_*')]
        if len(rows) != len(fields_Td):
            raise LoKIError('LoKI output count; possible rounded folder collision')
        rows.sort(key=lambda row: row['EN_Td'])
        fields_sorted = sorted(fields_Td)
        for row, field in zip(rows, fields_sorted):
            if not np.isclose(row['EN_Td'], field, rtol=spec['tolerances']['fingerprint_rtol'], atol=0):
                raise LoKIError('output field mismatch')
            row['EN_Td'] = float(field)
            row['u'] = float(np.log(field))
            row['composition'] = dict(coordinates)
            row['setup'] = job_name
            row['converged'] = (abs(row['power_groups']['relative_balance_percent']) / 100 <=
                                spec['solver_options']['numerics']['maxPowerBalanceRelError'])
            if not row['converged']:
                raise LoKIError('power convergence')
        return rows if fields_Td[0] <= fields_Td[-1] else rows[::-1]


def _enrich_row(row, channel_map, spec):
    """Add collision powers, attachment moments, physical inputs and floor flags."""
    caller_row = row
    row, channel_map, spec = copy.deepcopy((row, channel_map, spec))
    if row['channels'] != [entry['description'] for entry in channel_map]:
        raise LoKIError('channel map does not exactly cover the solver collision set')
    _validate_physical_map(channel_map, spec, row['composition'])
    state = {name: env['reference'] for name, env in spec['envelopes'].items()}
    state.update(row['composition'])
    _, _, gas, populations = _resolved_properties(spec, row['composition'])
    targets, products = _channel_fractions(channel_map, gas, populations, state)
    row['gas_fractions'] = gas
    row['state_populations'] = populations
    row['channel_power'] = np.zeros(len(channel_map))
    row['attachment_energy_eV'] = np.zeros(len(channel_map))
    row['target_fractions'] = np.asarray(targets)
    row['product_fractions'] = np.asarray(products)
    energy, f0 = row['energy_eV'], row['f0']
    width = np.diff(row['energy_edges_eV'])
    row['rate_floors'] = np.zeros(len(channel_map))
    for i, channel in enumerate(channel_map):
        xtgt, xprod = targets[i], products[i]
        row['target_fractions'][i] = xtgt
        row['channel_power'][i] = channel['threshold_eV'] * (xtgt * row['k_ine'][i] - xprod * row['k_sup'][i])
        sigma_v = max(channel['cross_section']['sigma_m2']) * np.sqrt(2 * elementary_charge * energy[-1] / electron_mass)
        row['rate_floors'][i] = max(spec['floors']['rate_absolute'], spec['floors']['eedf_dynamic_range'] * sigma_v)
        if channel['kind'] == 'attachment':
            section = channel['cross_section']
            sigma_nodes = np.interp(row['energy_edges_eV'], section['energy_eV'], section['sigma_m2'], left=0, right=0)
            sigma_nodes[row['energy_edges_eV'] <= channel['threshold_eV']] = 0
            sigma = (sigma_nodes[:-1] + sigma_nodes[1:]) / 2
            integrand = f0 * energy * sigma * np.sqrt(2 * elementary_charge / electron_mass) * width
            rate = np.sum(integrand)
            if rate > 0:
                row['attachment_energy_eV'][i] = np.sum(integrand * energy) / rate
            row['channel_power'][i] = xtgt * row['k_ine'][i] * row['attachment_energy_eV'][i]
    if all('cross_section' in channel for channel in channel_map):
        from rmgpy.tools.eedf.moments import _distribution_moments
        physical = dict(row, product_fractions=np.asarray(products),
                        gas_temperature_K=state.get('Tg_K', spec['Tg_K']))
        row['channel_power'] = _distribution_moments(physical, channel_map)['channel_power']
    elastic = [i for i, channel in enumerate(channel_map) if channel['kind'] == 'elastic']
    if elastic and 'power_groups' in row:
        if len(elastic) > 1 and not all('cross_section' in channel and 'mass_ratio' in channel for channel in channel_map if channel['kind'] == 'elastic'):
            raise LoKIError('per-channel elastic partition needs cross sections and masses')
        if len(elastic) == 1:
            row['channel_power'][elastic[0]] = -row['power_groups']['elastic_gain'] - row['power_groups']['elastic_loss']
    if 'power_groups' in row:
        from rmgpy.tools.eedf.validation import _channel_power_checks
        checks = _channel_power_checks(row, {'channel_map': channel_map, 'tolerances': spec['tolerances'], 'floors': spec['floors']})
        if any(not check['passed'] for check in checks):
            raise LoKIError('per-channel power sum disagrees with LoKI groups: ' + str([c for c in checks if not c['passed']]))
    row['below_floor'] = np.stack([row['k_ine'] < row['rate_floors'], row['k_sup'] < row['rate_floors']], axis=-1)
    caller_row.clear()
    caller_row.update(copy.deepcopy(row))
    return caller_row
