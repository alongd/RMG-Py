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

"""Run and audit a fixed-mechanism plasma benchmark; external data remain in scratch.

The LoKI arm is independent chemistry plus a separately implemented stationary
budget, not validation of LoKI's native electron energy equation. This runner
exports declared inputs, never RMG residual code or precomputed steady roots.
"""
import argparse
import hashlib
import json
import logging
import math
from pathlib import Path
import subprocess
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = {"PlasmaArgonLegacy:" + str(i) for i in range(86, 94)} | {
    "PlasmaRadiativeRecombination:1"}


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError("Expected exactly one input declaration: " + old)
    return text.replace(old, new)


def migrated_deck(source, pressure, power):
    """Only the approved library/state migration and sweep coordinates change."""
    text = source.replace("PlasmaArgon", "PlasmaArgonLegacy")
    text = replace_once(text, "multiplicity 3\n1 Ar u2 p3 c0",
                        "multiplicity 3\nelectronicstate 3P2\n1 Ar u2 p3 c0")
    text = replace_once(text, "pressure=(5, 'torr')",
                        "pressure=({!r}, 'torr')".format(pressure))
    return replace_once(text, "'absorbedPower': (0.5, 'W')",
                        "'absorbedPower': ({!r}, 'W')".format(power))


def initialize(directory, pressure, power):
    from rmgpy.rmg.main import RMG
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    source = ROOT / 'examples/rmg/plasma_argon_energy_balance/input.py'
    deck = directory / 'input.py'
    deck.write_text(migrated_deck(source.read_text(), pressure, power))
    output = directory / 'rmg'
    output.mkdir(exist_ok=True)
    job = RMG(input_file=str(deck), output_directory=str(output))
    job.initialize()
    core = job.reaction_model.core.species
    reactions = job.reaction_model.core.reactions
    if len(core) != 4 or {s.label for s in core} != {'Ar', 'Ars', 'Arp', 'e-'}:
        raise ValueError('Migration did not produce exactly the four intended species')
    from rmgpy.species import Species
    untagged = Species().from_adjacency_list('multiplicity 3\n1 Ar u2 p3 c0')
    candidates = []
    for species in core:
        candidate = species.molecule[0].copy(deep=True)
        candidate.electronic_state = ''
        if candidate.is_isomorphic(untagged.molecule[0]):
            candidates.append(species)
    if len(candidates) != 1 or candidates[0].label != 'Ars':
        raise ValueError('Untagged Ars does not resolve uniquely')
    if candidates[0].molecule[0].electronic_state != '3P2':
        raise ValueError('Ars is not the approved 3P2 identity')
    ids = [str(r.library) + ':' + str(r.entry.index) for r in reactions]
    if len(ids) != 9 or set(ids) != EXPECTED:
        raise ValueError('Unexpected mechanism reaction inventory: ' + repr(ids))
    reactor = job.reaction_systems[0]
    reactor.initialize_model(core, reactions, [], [],
                             atol=job.simulator_settings_list[0].atol,
                             rtol=job.simulator_settings_list[0].rtol)
    return job, reactor, core, reactions


def independent_rate(law, te_k, tg, constants):
    """Number-density rate from exported parameters, independent of kinetics APIs."""
    if law['kind'] == 'ElectronCollisionPlasma':
        energy = np.asarray(law['energy_J_particle'])
        sigma = np.asarray(law['sigma_m2'])
        kt = constants['kB'] * te_k
        return math.sqrt(8 / (math.pi * constants['m_e'])) * kt ** -1.5 * np.trapz(
            sigma * energy * np.exp(-energy / kt), energy)
    t = te_k if law['kind'] == 'TwoTemperaturePlasma' else tg
    value = law['A_particle'] * (t / law['T0_K']) ** law['n']
    if law['kind'] == 'TwoTemperaturePlasma':
        return value * math.exp(-law['Ea_g_J_mol'] / (constants['R'] * tg)
                                - law['Ea_e_J_mol'] / constants['R'] * (1 / te_k - 1 / tg))
    return value * math.exp(-law['Ea_J_mol'] / (constants['R'] * tg))


def export_manifest(job, reactor, core, reactions, directory):
    import rmgpy.constants as c
    from rmgpy import settings
    from rmgpy.solver.plasma import PLASMA_LOSCHMIDT
    constants = {name: float(getattr(c, name)) for name in ('kB', 'R', 'Na', 'e', 'm_e')}
    views = reactor._resolve_electron_placements(reactions, core, [])
    energy = reactor.electron_energy_balance
    rows = []
    preflight = []
    for reaction, view in zip(reactions, views):
        key = str(reaction.library) + ':' + str(reaction.entry.index)
        kin = reaction.kinetics
        order = len(view.reactants)
        law = {'kind': type(kin).__name__}
        if law['kind'] == 'ElectronCollisionPlasma':
            law.update(energy_J_particle=(kin.energies.value_si / c.Na).tolist(),
                       sigma_m2=kin.sigma.value_si.tolist())
        else:
            law.update(A_particle=float(kin.A.value_si / c.Na ** (order - 1)),
                       n=float(kin.n.value_si), T0_K=float(kin.T0.value_si))
            names = ('Ea_g', 'Ea_e') if law['kind'] == 'TwoTemperaturePlasma' else ('Ea',)
            for name in names:
                law[name + '_J_mol'] = float(getattr(kin, name).value_si)
        epsilon = energy['electron_energies'].get(key, (0., 'eV'))[0]
        rows.append({'id': key, 'reactants': [s.label for s in view.reactants],
                     'products': [s.label for s in view.products], 'law': law,
                     'electron_energy_eV': epsilon})
        for theta in np.geomspace(0.2, 10, 81):
            te_k = theta * c.e / c.kB
            if hasattr(kin, 'get_rate_coefficient_two_temp'):
                rate = kin.get_rate_coefficient_two_temp(reactor.T.value_si, te_k)
            elif getattr(kin, 'uses_electron_temperature', False):
                rate = kin.get_rate_coefficient_electron_temp(te_k)
            else:
                rate = kin.get_rate_coefficient(reactor.T.value_si)
            actual = rate / c.Na ** (order - 1)
            expected = independent_rate(law, te_k, reactor.T.value_si, constants)
            relative = abs(expected / actual - 1) if actual else abs(expected)
            preflight.append({'id': key, 'Te_eV': float(theta), 'rmg': actual,
                              'independent': expected, 'relative_error': relative})
    if max(x['relative_error'] for x in preflight) > 1e-10:
        raise ValueError('Rate preflight failed')
    y0 = np.asarray(reactor.y0)
    v0 = reactor.compute_volume(y0)
    heavy = sum(y0[i] for i, s in enumerate(core) if not s.is_electron())
    record = {'schema': 1, 'mechanism': 'historical lumped four-species argon',
              'independence': 'LoKI chemistry plus a separately implemented budget',
              'runner_sha256': file_hash(__file__),
              'rmg_sha': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
              'database_sha': subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                cwd=settings['database.directory'], text=True).strip(),
              'constants': constants, 'reactions': rows,
              'species': [{'label': s.label, 'adjacency': s.molecule[0].to_adjacency_list(),
                           'particle_mass_kg': float(s.molecular_weight.value_si)} for s in core],
              'Ars_resolution': {'candidates': 1, 'label': 'Ars', 'electronic_state': '3P2'},
              'settings': {'Tg_K': float(reactor.T.value_si), 'pressure_Pa': float(reactor.P.value_si),
                  'power_W': float(reactor.absorbed_power_total),
                  'chamber_volume_m3': float(energy['chamber_volume'][0]),
                  'initial_heavy_density_m3': float(heavy * c.Na / v0),
                  'power_reference_heavy_density_m3': float(reactor.P.value_si / ((c.R / c.Na) * reactor.T.value_si)),
                  'initial_electron_density_m3': float(y0[reactor.electron_index] * c.Na / v0),
                  'initial_Te_K': float(y0[reactor.te_index]),
                  'lambda_m': float(reactor.diffusion_length.value_si),
                  'mobility_density_product': 0.0001535 * float(PLASMA_LOSCHMIDT),
                  'neutral_diffusion_density_product': 47e-4 * 133.32236842105263 / ((c.R/c.Na) * reactor.T.value_si),
                  'source_m3_s': 6.6e4,
                  'elastic': energy['elastic_collisions']['Ar']}}
    directory = Path(directory)
    write_json(directory / 'mechanism.json', record)
    write_json(directory / 'rate-preflight.json', preflight)
    logging.info('RATE_PREFLIGHT_PASS: %d values; max relative error %.3g',
                 len(preflight), max(x['relative_error'] for x in preflight))
    return record


def run_point(job, reactor, core, reactions, directory):
    import rmgpy.constants as c
    from rmgpy.solver.electronegative import manifest_values
    started = time.monotonic()
    simulation = reactor.simulate(core, reactions, [], [], [], [],
                                 model_settings=job.model_settings_list[0],
                                 simulator_settings=job.simulator_settings_list[0])
    y = np.asarray(reactor.y)
    volume = reactor.compute_volume(y)
    result = {'status': 'steady' if reactor.steady_state_reached else 'not_steady',
              'runtime_seconds': time.monotonic() - started,
              'physical_time_seconds': float(reactor.t),
              'accepted_state_records': len(reactor.energy_history),
              'Te_K': float(y[reactor.te_index]),
              'solver_terminated': bool(simulation[0]),
              'Te_eV': float(y[reactor.te_index] * (c.R / c.Na) / c.e),
              'densities_m3': {s.label: float(y[i] * c.Na / volume) for i, s in enumerate(core)},
              'amounts_mol': {s.label: float(y[i]) for i, s in enumerate(core)},
              'volume_m3': float(volume),
              'budget_W_inventory': manifest_values(reactor.energy_budget),
              'reaction_rates_m3_s': (np.asarray(reactor.core_reaction_rates) * c.Na).tolist(),
              'history': [[float(v) for v in row] for row in reactor.energy_history]}
    write_json(Path(directory) / 'rmg-result.json', result)
    if not reactor.steady_state_reached:
        raise RuntimeError('RMG did not reach steady state; result recorded')
    logging.info('RMG_CASE_COMPLETED: ne=%.9g Te=%.9g eV',
                 result['densities_m3']['e-'], result['Te_eV'])
    return result


def sweep_coordinates():
    ascending = [(p, w) for p in (1., 5., 20.) for w in (.01, .1, .5, 5.)]
    return [('ascending', ascending), ('descending', list(reversed(ascending)))]


def apply_seed(reactor, core, reactions, result, simulator):
    """Use the previous accepted composition/temperature as a continuation seed."""
    from rmgpy.quantity import Temperature
    amounts = result['amounts_mol']
    total = sum(amounts.values())
    reactor.initial_mole_fractions = {species: amounts[species.label] / total for species in core}
    reactor.Te = Temperature((result['Te_K'], 'K'))
    reactor.te_initial = result['Te_K']
    reactor.initialize_model(core, reactions, [], [], atol=simulator.atol, rtol=simulator.rtol)


def run_sweep(directory):
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    records = []
    for direction, coordinates in sweep_coordinates():
        seed = None
        for pressure, power in coordinates:
            label = '{}_p{:g}_w{:g}'.format(direction, pressure, power)
            output = directory / label
            row = {'point': label, 'pressure_torr': pressure, 'power_W': power,
                   'direction': direction, 'seed_from': records[-1]['point'] if seed else None}
            try:
                job, reactor, core, reactions = initialize(output, pressure, power)
                if seed:
                    apply_seed(reactor, core, reactions, seed, job.simulator_settings_list[0])
                export_manifest(job, reactor, core, reactions, output)
                write_json(output / 'continuation.json', row)
                result = run_point(job, reactor, core, reactions, output)
                row.update(status=result['status'], ne_m3=result['densities_m3']['e-'],
                           Te_eV=result['Te_eV'], physical_time_seconds=result['physical_time_seconds'])
                seed = result
            except Exception as error:
                logging.exception('Failed benchmark point %s', label)
                row.update(status='failed', error=repr(error))
                seed = None
            records.append(row)
            write_json(directory / 'sweep-index.json', records)
    if any(row['status'] != 'steady' for row in records):
        raise RuntimeError('One or more sweep points failed; see sweep-index.json')
    return records


def compare_point(rmg, loki):
    """Compare verified endpoints; a measured mismatch is a finding, not success."""
    if rmg.get('status') != 'steady':
        raise ValueError('RMG endpoint is not a demonstrated steady state')
    for name, tolerance in [('power_relative_residual', 1e-6),
                            ('charge_relative_residual', 1e-6),
                            ('eos_relative_residual', 1e-7),
                            ('particle_relative_residual', 1e-6)]:
        if name not in loki or not math.isfinite(loki[name]) or abs(loki[name]) >= tolerance:
            raise ValueError('LoKI endpoint failed ' + name)
    import rmgpy.constants as c
    budget = rmg['budget_W_inventory']
    te = budget['Te'] * (c.R / c.Na) / c.e
    volume = rmg['volume_m3']
    pairs = {
        'ne_m3': (rmg['densities_m3']['e-'], loki['ne_m3']),
        'Te_eV': (te, loki['Te_eV']),
        'Arp_m3': (rmg['densities_m3']['Arp'], loki['Arp_m3']),
        'Ars_m3': (rmg['densities_m3']['Ars'], loki['Ars_m3']),
        'absorbed_W_m3': (budget['P_abs'] / volume, loki['P_abs_W_m3']),
        'inelastic_W_m3': (budget['Q_inelastic'] / volume,
                          sum(loki['reaction_power_W_m3']) + sum(loki['consumed_electron_power_W_m3'])),
        'elastic_W_m3': (budget['Q_elastic'] / volume, loki['elastic_W_m3']),
        'wall_electron_W_m3': (budget['Q_wall_electron'] / volume, loki['wall_electron_W_m3']),
        'wall_ion_W_m3': (budget['Q_wall_ion'] / volume, loki['wall_ion_W_m3']),
        'wall_total_W_m3': ((budget['Q_wall_electron'] + budget['Q_wall_ion']) / volume,
                            loki['wall_electron_W_m3'] + loki['wall_ion_W_m3']),
    }
    rows = []
    for name, (a, b) in pairs.items():
        if not math.isfinite(a) or not math.isfinite(b) or b == 0:
            raise ValueError('Nonfinite/zero comparison quantity: ' + name)
        rows.append({'quantity': name, 'rmg': a, 'loki': b, 'relative_difference': a / b - 1})
    differences = {r['quantity']: abs(r['relative_difference']) for r in rows}
    return {'independence': 'LoKI chemistry plus a separately implemented budget',
            'acceptance': differences['ne_m3'] <= .02 and differences['Te_eV'] <= .01,
            'rows': rows, 'rmg_runtime_seconds': rmg['runtime_seconds'],
            'loki_runtime_seconds': loki['runtime_total_seconds'],
            'rmg_physical_settling_seconds': rmg['physical_time_seconds'],
            'loki_physical_settling_seconds': None,
            'loki_outer_iterations': loki['outer_iterations']}


def comparison_value(a, b):
    if not math.isfinite(a) or not math.isfinite(b):
        raise ValueError('Nonfinite comparison value')
    relative = a / b - 1 if b else (0. if a == 0 else None)
    return {'rmg': a, 'loki': b, 'absolute_difference': a - b,
            'relative_difference': relative}


def compare_reactions(rmg, loki, rmg_manifest, loki_manifest):
    """Align original laws by identity, keeping grouped LoKI powers separate."""
    r_sources = rmg_manifest['reactions']
    l_sources = loki_manifest['reactions']
    r_ids = [row['id'] for row in r_sources]
    l_ids = [row['id'] for row in l_sources]
    if len(set(r_ids)) != len(r_ids) or len(set(l_ids)) != len(l_ids):
        raise ValueError('Duplicate source-law identity')
    if {row['id']: row for row in r_sources} != {row['id']: row for row in l_sources}:
        raise ValueError('Source-law manifests differ')
    r_rates = rmg['reaction_rates_m3_s']
    l_rates = loki['source_rate_m3_s']
    r_power = rmg['budget_W_inventory']['Q_inelastic_by_reaction']
    l_power = loki['reaction_power_W_m3']
    l_consumed = loki['consumed_electron_power_W_m3']
    if any(len(values) != len(r_ids) for values in (r_rates, r_power)) or any(
            len(values) != len(l_ids) for values in (l_rates, l_power, l_consumed)):
        raise ValueError('Reaction vector length differs from source-law manifest')
    rows = []
    for i, identity in enumerate(r_ids):
        j = l_ids.index(identity)
        rate = comparison_value(r_rates[i], l_rates[j])
        power = comparison_value(r_power[i] / rmg['volume_m3'], l_power[j] + l_consumed[j])
        rows.append({'reaction_id': identity,
                     **{'rate_' + key: value for key, value in rate.items()},
                     **{'power_' + key: value for key, value in power.items()}})
    return rows


def comparison_markdown(result):
    lines = ['LoKI chemistry plus a separately implemented budget; historical lumped Case A.', '',
             '| Quantity | RMG | LoKI | Relative difference (RMG/LoKI − 1) |',
             '|---|---:|---:|---:|']
    for row in result['rows']:
        lines.append('| {quantity} | {rmg:.12g} | {loki:.12g} | {relative_difference:.6g} |'.format(**row))
    lines += ['', '| Reaction | RMG rate (m⁻³ s⁻¹) | LoKI rate (m⁻³ s⁻¹) | Relative difference | RMG power (W m⁻³) | LoKI power (W m⁻³) |',
              '|---|---:|---:|---:|---:|---:|---:|']
    for row in result['reaction_rows']:
        relative = row['rate_relative_difference']
        value = 'undefined (zero LoKI rate)' if relative is None else '{:.6g}'.format(relative)
        lines.append('| {} | {:.12g} | {:.12g} | {} | {:.12g} | {:.12g} |'.format(
            row['reaction_id'], row['rate_rmg'], row['rate_loki'], value,
            row['power_rmg'], row['power_loki']))
    lines += ['', 'Acceptance: **{}** (ne ≤ 2%, Te ≤ 1%).'.format('PASS' if result['acceptance'] else 'OUTSIDE TOLERANCE'),
              'Relative differences are fractions, not percentages; zero/zero reaction powers are recorded as equal.', '']
    return '\n'.join(lines)


def plot_comparison(result, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    names = ['elastic_W_m3', 'inelastic_W_m3', 'wall_electron_W_m3', 'wall_ion_W_m3']
    values = {row['quantity']: row for row in result['rows']}
    x = np.arange(len(names))
    axes[0].bar(x - .18, [values[name]['rmg'] for name in names], .36, label='RMG')
    axes[0].bar(x + .18, [values[name]['loki'] for name in names], .36, label='LoKI + independent budget')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(['Elastic', 'Inelastic', 'Electron wall', 'Ion wall'], rotation=25, ha='right')
    axes[0].set_yscale('log')
    axes[0].set_ylabel('Power loss (W m$^{-3}$)')
    axes[0].legend(fontsize=8)
    rates = result['reaction_rows']
    axes[1].barh(np.arange(len(rates)), [abs(row['rate_relative_difference']) if row['rate_relative_difference'] is not None else np.nan for row in rates])
    axes[1].set_yticks(np.arange(len(rates)))
    axes[1].set_yticklabels([row['reaction_id'].replace('PlasmaArgonLegacy:', 'Legacy ').replace('PlasmaRadiativeRecombination:', 'RR ') for row in rates])
    axes[1].set_xlabel('Absolute relative reaction-rate difference')
    axes[1].ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
    conditions = result['conditions']
    fig.suptitle('Historical lumped Case A: {:g} W, {:g} torr'.format(
        conditions['power_W'], conditions['pressure_torr']))
    fig.tight_layout()
    fig.savefig(output / 'comparison.png', dpi=180)
    fig.savefig(output / 'comparison.svg')
    plt.close(fig)


def write_comparison(rmg_path, loki_path, output, loki_manifest_path, plot=False):
    import csv
    rmg_path, loki_path = Path(rmg_path), Path(loki_path)
    rmg = json.loads(rmg_path.read_text())
    loki = json.loads(loki_path.read_text())
    rmg_manifest_path = rmg_path.parent / 'mechanism.json'
    result = compare_point(rmg, loki)
    settings = json.loads(rmg_manifest_path.read_text())['settings']
    result['conditions'] = {'power_W': settings['power_W'],
                            'pressure_torr': settings['pressure_Pa'] / 133.32236842105263}
    result['reaction_rows'] = compare_reactions(rmg, loki,
        json.loads(rmg_manifest_path.read_text()), json.loads(Path(loki_manifest_path).read_text()))
    result['input_sha256'] = {str(path): file_hash(path) for path in
        (rmg_path, loki_path, rmg_manifest_path, Path(loki_manifest_path))}
    result['script_sha256'] = file_hash(Path(__file__))
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    write_json(output / 'comparison.json', result)
    for name, rows in [('comparison.csv', result['rows']), ('reactions.csv', result['reaction_rows'])]:
        with (output / name).open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    (output / 'comparison.md').write_text(comparison_markdown(result))
    if plot:
        plot_comparison(result, output)
    logging.info('COMPARISON_RECORDED acceptance=%s', result['acceptance'])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['preflight', 'point', 'sweep', 'compare'])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--pressure', type=float, default=5.)
    parser.add_argument('--power', type=float, default=.5)
    parser.add_argument('--seed', type=Path, help='Previous accepted RMG result for continuation')
    parser.add_argument('--rtol', type=float, help='Explicit integration refinement')
    parser.add_argument('--atol', type=float, help='Explicit integration refinement')
    parser.add_argument('--rmg-result', type=Path)
    parser.add_argument('--loki-result', type=Path)
    parser.add_argument('--loki-manifest', type=Path)
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    if args.action == 'compare':
        if not args.rmg_result or not args.loki_result or not args.loki_manifest:
            parser.error('compare requires --rmg-result, --loki-result and --loki-manifest')
        write_comparison(args.rmg_result, args.loki_result, args.output, args.loki_manifest, args.plot)
        return
    if args.action == 'sweep':
        run_sweep(args.output)
        return
    job, reactor, core, reactions = initialize(args.output, args.pressure, args.power)
    simulator = job.simulator_settings_list[0]
    if args.rtol is not None:
        simulator.rtol = args.rtol
    if args.atol is not None:
        simulator.atol = args.atol
    if args.seed:
        seed = json.loads(args.seed.read_text())
        if seed.get('status') != 'steady':
            parser.error('Continuation requires a demonstrated steady seed')
        if 'Te_K' not in seed:
            seed['Te_K'] = seed['budget_W_inventory']['Te']
        apply_seed(reactor, core, reactions, seed, simulator)
        write_json(args.output / 'continuation.json', {'seed_path': str(args.seed.resolve()), 'seed_sha256': file_hash(args.seed)})
    export_manifest(job, reactor, core, reactions, args.output)
    if args.action == 'point':
        run_point(job, reactor, core, reactions, args.output)


if __name__ == '__main__':
    main()
