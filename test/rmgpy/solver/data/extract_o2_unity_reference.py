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

"""Extract the I-314 rework-3 O2 unity-closure test oracle.

The source directory must contain ``rw_core.py``, ``sweep_sc.csv`` and
``materiality_rw.csv`` from the same I-314 rework-3 revision.  The script
re-solves the spatial references at the 18 CSV-listed unity operating states;
it does not make the resulting vectors available to production code.

Usage::

    python extract_o2_unity_reference.py SOURCE_DIR OUTPUT_JSON
"""

import argparse
import csv
import json
import multiprocessing
from pathlib import Path
import sys

import numpy as np


VARIANTS = (
    ('nominal', {}),
    ('mob200', {'ion_field_mob': 200.}),
    ('mob50', {'ion_field_mob': 50.}),
    ('vth05', {'vth_scale': 0.5}),
    ('vth2', {'vth_scale': 2.}),
    ('c50v05', {'ion_field_mob': 50., 'vth_scale': 0.5}),
    ('c50v2', {'ion_field_mob': 50., 'vth_scale': 2.}),
)


def rows_with_lines(path):
    with path.open(newline='') as handle:
        return [(line, dict(row))
                for line, row in enumerate(csv.DictReader(handle), 2)]


def number(row, key):
    return float(row[key])


def reference_is_valid(row):
    if row['status'] != 'ok' or not abs(number(row, 'J_rel')) < 1.e-6:
        return False
    if not number(row, 'f_hf') <= 0.01:
        return False
    for key in ('E_p_mob200', 'E_p_vth05', 'E_p_vth2'):
        if not np.isfinite(number(row, key)) or not number(row, key) <= 0.02:
            return False
    return True


def csv_variant_key(variant):
    return 'E_bracket_act' + ('' if variant == 'nominal' else '_' + variant)


def initialise_worker(source_dir):
    sys.path.insert(0, source_dir)
    global O, W, enref
    import o2_state as O
    import rw_core as W
    import enref


def extract_point(payload):
    sweep_line, row, materiality_line, materiality, confined_line, confined, confined_materiality_line, confined_materiality = payload
    state = O.state(number(row, 'p'), number(row, 'te_eff'), number(row, 'ne'),
                    row['case'], row['closure'])
    column, qn = W.solve_qn_alpha(state)
    poisson, nominal = W.solve_poisson_alpha(state, column, qn)
    references = {}
    resolved = nominal
    engine_frequencies = W.nu_ep(state, state['t_tr'])
    engine_cations = W.engine_cations(state, 1., state['t_tr'])
    engine = engine_frequencies * engine_cations
    for variant, options in VARIANTS:
        if options:
            resolved = W.poisson_resolve(poisson, state, nominal, **options)
        else:
            resolved = nominal
        reference = resolved['nu_ref'] * resolved['nbar_cat']
        error = enref.e_wall(engine, reference, np.ones_like(engine))
        csv_error = number(row, csv_variant_key(variant))
        if not np.isclose(error, csv_error, rtol=5.e-5, atol=5.e-7):
            raise RuntimeError(
                're-solved error does not reproduce sweep_sc.csv line {} {}: '
                '{} != {}'.format(sweep_line, variant, error, csv_error))
        references[variant] = {
            'reference_cation_wall_vector': reference.tolist(),
            'error': float(error),
        }

    valid = reference_is_valid(row)
    status = ('reference-qualified'
              if valid else
              'reference-unqualified (ion-heating validity), B3-pass numerically')
    return {
        'condition': {
            'pressure_o2_torr': number(row, 'p'),
            'absorbed_power_w': number(row, 'P'),
        },
        'source_rows': {
            'sweep_sc.csv': sweep_line,
            'materiality_rw.csv': materiality_line,
            'confined_sweep_sc.csv': confined_line,
            'confined_materiality_rw.csv': confined_materiality_line,
        },
        'state': {
            'electron_temperature_effective_eV': number(row, 'te_eff'),
            'mean_electron_energy_eV': number(row, 'eps'),
            'electron_density_m_3': number(row, 'ne'),
            'alpha': number(row, 'alpha'),
            'neutral_density_m_3': number(row, 'n_m'),
            'transport_energy_eV': number(row, 't_tr'),
        },
        'closure_state': [number(row, 'ne'), number(row, 'alpha') * number(row, 'ne')],
        'engine_cation_number_densities_m_3': engine_cations.tolist(),
        'reference_geometry_engine_frequencies_s_1': engine_frequencies.tolist(),
        'engine_unity_cation_wall_vector': engine.tolist(),
        'reference_status': status,
        'budgets': {
            budget: number(materiality, 'eps_uniform_' + budget)
            for budget in ('B1', 'B2', 'B3')
        },
        'reference_variants': references,
        'falsified_comparator': {
            'local_checks': {
                'minimum_confinement_ratio': number(confined, 'conf_min'),
                'attachment_metric': number(confined, 'att_tau'),
            },
            'actual_error': number(confined, 'E_simple_act'),
            'B3_budget': number(confined_materiality, 'eps_uniform_B3'),
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('source_dir', type=Path)
    parser.add_argument('output_json', type=Path)
    parser.add_argument('--workers', type=int, default=6)
    args = parser.parse_args()

    source = args.source_dir.resolve()
    sweep = rows_with_lines(source / 'sweep_sc.csv')
    materiality = rows_with_lines(source / 'materiality_rw.csv')
    mat_by_key = {
        (row['arm'], number(row, 'p'), number(row, 'P')): (line, row)
        for line, row in materiality
    }
    sweep_by_key = {
        (row['arm'], number(row, 'p'), number(row, 'P')): (line, row)
        for line, row in sweep
    }
    payloads = []
    for line, row in sweep:
        if row['arm'] != 'electropositiveBracket':
            continue
        key = (number(row, 'p'), number(row, 'P'))
        mat_line, mat = mat_by_key[('electropositiveBracket',) + key]
        confined_line, confined = sweep_by_key[('confinedAnion',) + key]
        confined_mat_line, confined_mat = mat_by_key[('confinedAnion',) + key]
        payloads.append((line, row, mat_line, mat, confined_line, confined,
                         confined_mat_line, confined_mat))
    if len(payloads) != 18:
        raise RuntimeError('expected 18 unity operating states, found {}'.format(len(payloads)))

    with multiprocessing.Pool(
            args.workers, initializer=initialise_worker,
            initargs=(str(source),)) as pool:
        points = list(pool.imap(extract_point, payloads))
    points.sort(key=lambda point: (
        point['condition']['pressure_o2_torr'],
        point['condition']['absorbed_power_w']))

    valid = [point for point in points
             if point['reference_status'] == 'reference-qualified']
    if len(valid) != 17:
        raise RuntimeError('expected 17 qualified reference states, found {}'.format(len(valid)))
    tally = {}
    for budget in ('B1', 'B2', 'B3'):
        tally[budget] = sum(
            max(variant['error'] for variant in point['reference_variants'].values())
            <= point['budgets'][budget]
            for point in valid)
    if tally != {'B1': 0, 'B2': 7, 'B3': 17}:
        raise RuntimeError('unexpected robust budget tallies: {!r}'.format(tally))

    output = {
        'schema_version': 2,
        'provenance': (
            'I-314 rework 3; rw_core.py frozen actual-vector error definition; '
            'sweep_sc.csv and materiality_rw.csv source rows recorded per point'),
        'cation_order': ['Ar+', 'O+', 'O2+'],
        'cation_charges': [1, 1, 1],
        'selected_closure': 'o2ReferenceQualifiedUnity',
        'legacy_source_arm': 'electropositiveBracket',
        'robust_valid_point_tallies': tally,
        'points': points,
    }
    args.output_json.write_text(json.dumps(output, indent=2, sort_keys=True) + '\n')


if __name__ == '__main__':
    main()
