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

"""Summarize retained plasma benchmark outputs without rerunning either solver.

External source, cross sections and result files stay outside the repository.
A native iteration failure and a fixed-state power decomposition never count as
converged native-arm validation. The scratch layout is produced by the benchmark
runner and its externally staged chemistry adapter.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import tarfile


def load(path):
    return json.loads(Path(path).read_text())


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finite(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('Nonfinite benchmark value')
    return value


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def write_csv(path, rows):
    if not rows:
        raise ValueError('Empty result table')
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def compare_directions(rows):
    """Require complete pairs; failed rows never borrow the opposite endpoint."""
    indexed = {}
    for row in rows:
        key = (row['direction'], finite(row['pressure_torr']), finite(row['power_W']))
        if key[0] not in ('ascending', 'descending') or key in indexed:
            raise ValueError('Invalid or duplicate sweep coordinate')
        if row['status'] not in ('steady', 'failed', 'solver_refused', 'not_steady'):
            raise ValueError('Unknown sweep status')
        indexed[key] = row
    coordinates = sorted({key[1:] for key in indexed})
    if not coordinates:
        raise ValueError('Empty sweep')
    output = []
    for pressure, power in coordinates:
        try:
            up, down = [indexed[(direction, pressure, power)] for direction in ('ascending', 'descending')]
        except KeyError as error:
            raise ValueError('Missing continuation direction') from error
        pair = {'pressure_torr': pressure, 'power_W': power}
        for direction, row in [('ascending', up), ('descending', down)]:
            status = row['status']
            if status == 'failed' and row.get('error', '').startswith('PlasmaStateError('):
                status = 'solver_refused'
            pair[direction + '_status'] = status
            for name in ('ne_m3', 'Te_eV'):
                value = finite(row[name]) if row['status'] == 'steady' else None
                if value is not None and value <= 0:
                    raise ValueError('Nonpositive steady endpoint')
                pair[direction + '_' + name] = value
        for name in ('ne_m3', 'Te_eV'):
            a, b = pair['ascending_' + name], pair['descending_' + name]
            pair[name + '_relative_difference'] = abs(a / b - 1) if a is not None and b is not None else None
        output.append(pair)
    return output


def power_accounting(base, text, manifest, groups):
    """Account for archived channels at one fixed state, using source identities."""
    normalized = {}
    # Gas-specific breakdowns repeat mixture labels; only consume the mixture block.
    for line in text.split('*************************************', 1)[0].splitlines():
        match = re.match(r'\s*(.*?)\s*=\s*([+-]?[\d.]+e[+-]\d+)\s*\(eVm\^3s\^-1\)', line)
        if match:
            name, value = match[1], finite(match[2])
            if name in normalized and normalized[name] != value:
                raise ValueError('Conflicting native power labels')
            normalized[name] = value
    names = ('Field', 'Elastic collisions (net)', 'Excitation collisions (net)', 'Ionization collisions')
    if any(name not in normalized for name in names) or normalized['Field'] <= 0:
        raise ValueError('Missing or invalid native power channels')
    ids = [row['id'] for row in manifest['reactions']]
    reaction = [finite(value) for value in base['reaction_power_W_m3']]
    if len(set(ids)) != len(ids) or len(ids) != len(reaction):
        raise ValueError('Invalid source-law power mapping')
    if set(groups) != {'excitation', 'ionization'}:
        raise ValueError('Supply excitation and ionization source groups')
    owners = [name for group in groups.values() for name in group]
    if len(set(owners)) != len(owners) or any(name not in ids for name in owners):
        raise ValueError('Missing or shared power-group identity')
    by_id = dict(zip(ids, reaction))
    excitation = sum(by_id[name] for name in groups['excitation'])
    ionization = sum(by_id[name] for name in groups['ionization'])
    scale = finite(base['native_field_W_m3']) / normalized['Field']
    parts = [
        ('elastic', -normalized[names[1]] * scale, finite(base['elastic_W_m3'])),
        ('excitation', -normalized[names[2]] * scale, excitation),
        ('ionization', -normalized[names[3]] * scale, ionization),
        ('remaining signed chemistry', 0., sum(reaction) - excitation - ionization),
        ('consumed electron energy', 0., sum(finite(x) for x in base['consumed_electron_power_W_m3'])),
        ('electron wall', 0., finite(base['wall_electron_W_m3'])),
        ('ion sheath', 0., finite(base['wall_ion_W_m3'])),
    ]
    rows = [{'term': term, 'native_W_m3': a, 'matched_W_m3': b, 'native_minus_matched_W_m3': a-b}
            for term, a, b in parts]
    for column, reference in [('native_W_m3', 'native_field_W_m3'), ('matched_W_m3', 'loss_W_m3')]:
        if not math.isclose(sum(row[column] for row in rows), finite(base[reference]), rel_tol=1e-10, abs_tol=1e-9):
            raise ValueError('Power channels do not close: ' + column)
    return {'scope': 'archived matched-state accounting; not ablation-confirmed', 'ablation_confirmed': False,
            'rows': rows, 'native_field_W_m3': finite(base['native_field_W_m3']),
            'matched_absorbed_W_m3': finite(base['P_abs_W_m3']),
            'gap_to_matched_loss_W_m3': sum(row['native_minus_matched_W_m3'] for row in rows),
            'matched_loss_minus_absorbed_W_m3': finite(base['loss_W_m3'])-finite(base['P_abs_W_m3'])}


def parse_iterations(stdout, config):
    """Reconstruct printed iterates; never relabel them as accepted endpoints."""
    te = finite(config['workingConditions']['electronTemperature'])
    ne = finite(config['workingConditions']['electronDensity'])
    rows = []
    for line in stdout.splitlines():
        match = re.search(r'New neutrality cycle iteration \((\d+)\): relative error = ([+\-\deE.]+)', line)
        if match:
            rows.append({'iteration': int(match[1]), 'Te_eV': te, 'ne_used_m3': ne,
                         'relative_error': finite(match[2])})
        match = re.search(r'Updated electron density \(([+\-\deE.]+) m\^-3\)', line)
        if match:
            ne = finite(match[1])
        match = re.search(r'Updated electron temperature \(([+\-\deE.]+) eV\)', line)
        if match:
            te = finite(match[1])
    if not rows or [row['iteration'] for row in rows] != list(range(1, len(rows)+1)):
        raise ValueError('Missing or duplicated native iteration records')
    return rows


def verify_manifest(manifest_path, root):
    root = Path(root).resolve()
    entries = Path(manifest_path).read_text().splitlines()
    if not entries:
        raise ValueError('Empty hash manifest')
    seen = set()
    for line in entries:
        expected, name = line.split(None, 1)
        path = (root / name.strip()).resolve()
        if root not in path.parents or path in seen or digest(path) != expected:
            raise ValueError('Hash manifest mismatch: ' + name)
        seen.add(path)
    return len(seen)


def verify_archive(fetch):
    lines = (fetch / 'output.sha256').read_text().splitlines()
    if not lines:
        raise ValueError('Empty archived-output manifest')
    with tarfile.open(fetch / 'loki-output.tar.gz') as archive:
        members = [member for member in archive.getmembers() if member.isfile()]
        by_name = {member.name: member for member in members}
        if len(by_name) != len(members):
            raise ValueError('Duplicate archive members')
        names = []
        for line in lines:
            expected, name = line.split(None, 1)
            name = name.strip()
            with archive.extractfile(by_name[name]) as stream:
                actual = hashlib.sha256(stream.read()).hexdigest()
            if actual != expected or name in names:
                raise ValueError('Archived output hash mismatch')
            names.append(name)
    return names


def audit_native(fetch, stage):
    """Accept a non-convergence finding only with its failure and integrity evidence."""
    fetch, stage = Path(fetch), Path(stage)
    expected = {'job': 1, 'matlab-stage': 1, 'pack': 0, 'upstream-check': 0, 'adapter-check': 0}
    markers = {name: int((fetch / (name + '.exit')).read_text()) for name in expected}
    if markers != expected or int((fetch / 'N0/trial-01.exit').read_text()) != 1:
        raise ValueError('Unexpected native exit markers')
    counts = {name: verify_manifest(fetch / (name + '.sha256'), stage) for name in ('upstream', 'adapter')}
    outputs = verify_archive(fetch)
    config = load(fetch / 'N0/trial-01-config.json')
    cap = config['chemistry']['iterationSchemes']['neutralityMaxIterations']
    tolerance = config['chemistry']['iterationSchemes']['neutralityRelError']
    rows = parse_iterations((fetch / 'N0/trial-01.stdout').read_text(), config)
    stderr = ' '.join((fetch / 'N0/trial-01.stderr').read_text().split())
    if ('Chemistry/neutralityCycle' not in stderr or
            'Maximum number of iterations reached for the neutrality cycle without convergence.' not in stderr):
        raise ValueError('Failure is not native neutrality non-convergence')
    index = load(fetch / 'native-index.json')
    ids = [row['arm_id'] for row in index]
    expected_ids = [arm['id'] for arm in load(stage / 'native-modes.json')]
    if (ids != expected_ids or len(ids) != 17 or len(set(ids)) != 17 or ids[0] != 'N0' or index[0]['status'] != 'failed'
            or any(row['status'] != 'not_run' for row in index[1:])):
        raise ValueError('Native arm statuses disagree with failure claim')
    if len(rows) != cap or min(abs(row['relative_error']) for row in rows) <= tolerance:
        raise ValueError('Iteration history does not establish cap exhaustion')
    if not any(a['relative_error']*b['relative_error'] < 0 for a, b in zip(rows, rows[1:])):
        raise ValueError('No observed oscillation in native history')
    return {'status': 'native-non-convergent', 'behavior': 'oscillatory', 'iteration_cap': cap,
            'target_relative_error': tolerance, 'best_iteration': min(rows, key=lambda row: abs(row['relative_error'])),
            'last_iteration': rows[-1], 'sign_changes': sum(a['relative_error']*b['relative_error'] < 0
                                                          for a, b in zip(rows, rows[1:])),
            'other_arms_not_run': len(index)-1, 'ablation_confirmed': False,
            'exit_markers': markers, 'verified_hash_counts': counts, 'verified_archived_outputs': outputs,
            'rows': rows}


def render_markdown(result):
    lines = ['## RMG continuation sweep', '',
             'Direction consistency only; this is not a two-code sweep validation.', '',
             '| Pressure / torr | Power / W | Upward ne / m^-3 | Downward ne / m^-3 | Upward Te / eV | Downward Te / eV |',
             '|---|---|---|---|---|---|']
    def cell(row, direction, quantity):
        value = row[direction + '_' + quantity]
        return '{:.10g}'.format(value) if value is not None else row[direction + '_status']
    for row in result['directions']:
        lines.append('| {:g} | {:g} | {} | {} | {} | {} |'.format(
            row['pressure_torr'], row['power_W'], cell(row, 'ascending', 'ne_m3'), cell(row, 'descending', 'ne_m3'),
            cell(row, 'ascending', 'Te_eV'), cell(row, 'descending', 'Te_eV')))
    retry = result['retry']
    lines += ['', 'The upward refusal is retained separately from its downward endpoint. Original atol={original_atol:g}; '
              'retry atol={retry_atol:g}, rtol={retry_rtol:g}; retry status: **{retry_status}**.'.format(**retry), '',
              '## Archived-state power attribution — not ablation-confirmed', '',
              '| Term | Native / W m^-3 | Matched / W m^-3 | Native minus matched / W m^-3 |',
              '|---|---:|---:|---:|']
    for row in result['power']['rows']:
        lines.append('| {term} | {native_W_m3:.9g} | {matched_W_m3:.9g} | {native_minus_matched_W_m3:.9g} |'.format(**row))
    lines += ['', 'All terms are evaluated at the archived matched state. This attribution is **not ablation-confirmed**.', '',
              '## Native attempts', '', '| Cap | Best absolute residual | Final residual | Sign changes | Unrun arms |',
              '|---|---:|---:|---:|---:|']
    for native in result['native']:
        lines.append('| {} | {:.9g} | {:.9g} | {} | {} |'.format(native['iteration_cap'],
            abs(native['best_iteration']['relative_error']), native['last_iteration']['relative_error'],
            native['sign_changes'], native['other_arms_not_run']))
    lines += ['', 'Final verdict: **N0 native-non-convergent (oscillatory, {} iterations)**. '
              'No accepted native endpoint or ablation result exists; the other 16 arms did not run.'.format(result['native'][-1]['iteration_cap']), '']
    return '\n'.join(lines)


def plot_results(result, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True, layout='constrained')
    for native in result['native']:
        rows = native['rows']
        x = [row['iteration'] for row in rows]
        label = 'cap {}'.format(native['iteration_cap'])
        axes[0].semilogy(x, [abs(row['relative_error']) for row in rows], label=label)
        axes[1].plot(x, [row['Te_eV'] for row in rows], label=label)
    axes[0].axhline(result['native'][-1]['target_relative_error'], color='red', ls='--', label='target')
    axes[0].set_ylabel('Absolute neutrality residual')
    axes[0].set_title('Native N0 non-convergence; iterates are not accepted endpoints')
    axes[1].set_ylabel('Electron temperature / eV')
    axes[1].set_xlabel('Native neutrality iteration')
    for axis in axes:
        axis.legend()
        axis.grid(alpha=.2)
    for suffix in ('png', 'svg'):
        fig.savefig(output / ('native-nonconvergence.' + suffix), dpi=160)
    plt.close(fig)


def generate(scratch, output, groups, plot=False):
    scratch, output = Path(scratch).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    inputs = [scratch / 'case-a-rmg-sweep/sweep-index.json', scratch / 'case-a-low-power-retry/attempts.json',
              scratch / 'case-a-base-r4-fetch/loki-result.json', scratch / 'case-a-comparison/powerBalance.txt',
              scratch / 'zeus-case-a-base-r4/mechanism.json']
    sweep, retry, base = [load(path) for path in inputs[:3]]
    expected = {(direction, pressure, power) for direction in ('ascending', 'descending')
                for pressure in (1., 5., 20.) for power in (.01, .1, .5, 5.)}
    if len(sweep) != 24 or {(row['direction'], row['pressure_torr'], row['power_W']) for row in sweep} != expected:
        raise ValueError('Incomplete Case A sweep inventory')
    if retry['original_attempt'] not in sweep or retry['retry_status'] != 'solver_refused':
        raise ValueError('Retry record does not establish the accepted solver refusal')
    # These checks preserve the actual retry's provenance, not a newly rerun point.
    if digest(inputs[0]) != retry['original_sweep_index_sha256'] or digest(retry['seed_path']) != retry['seed_sha256']:
        raise ValueError('Original sweep or retry seed changed')
    result = {'script_sha256': digest(__file__), 'inputs': {str(path): digest(path) for path in inputs},
              'directions': compare_directions(sweep), 'retry': retry,
              'power': power_accounting(base, inputs[3].read_text(), load(inputs[4]), groups), 'native': []}
    for name, stage in [('native-ablations-fetch', 'zeus-native-ablations-r1'),
                        ('native-cap500-fetch', 'zeus-native-ablations-cap500-proposed')]:
        native = audit_native(scratch / name, scratch / stage)
        result['native'].append(native)
        write_csv(output / ('native-{}.csv'.format(native['iteration_cap'])), native['rows'])
    write_csv(output / 'direction-comparison.csv', result['directions'])
    write_csv(output / 'fixed-state-power.csv', result['power']['rows'])
    write_json(output / 'summary.json', result)
    (output / 'report-section.md').write_text(render_markdown(result))
    if plot:
        plot_results(result, output)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scratch', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--power-groups', type=Path,
                        help='JSON mapping of excitation/ionization to source-law identities')
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()
    groups = load(args.power_groups) if args.power_groups else {
        'excitation': ['PlasmaArgonLegacy:87'], 'ionization': ['PlasmaArgonLegacy:86', 'PlasmaArgonLegacy:88']}
    generate(args.scratch, args.output, groups, args.plot)


if __name__ == '__main__':
    main()
