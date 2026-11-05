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

"""Checks for independent benchmark rate evaluation and explicit input migration."""
import importlib.util
from pathlib import Path
import math

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('plasma_benchmark', ROOT / 'scripts/plasmaBenchmark.py')
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)


def test_migration_is_explicit_and_preserves_mechanism_numbers():
    text = (ROOT / 'examples/rmg/plasma_argon_energy_balance/input.py').read_text()
    migrated = benchmark.migrated_deck(text, 5, .5)
    expected = text.replace('PlasmaArgon', 'PlasmaArgonLegacy').replace(
        'multiplicity 3\n1 Ar u2 p3 c0', 'multiplicity 3\nelectronicstate 3P2\n1 Ar u2 p3 c0')
    assert migrated == expected
    assert "'PlasmaArgonLegacy:90': (-7.337097999999999, 'eV')" in migrated
    with pytest.raises(ValueError, match='exactly one'):
        benchmark.migrated_deck(text.replace('multiplicity 3', 'multiplicity 5'), 5, .5)


def test_maxwellian_integral_has_correct_particle_units():
    c = {'kB': 1.380649e-23, 'm_e': 9.1093837e-31}
    temperature = 10000.
    energy = np.linspace(0, 50 * c['kB'] * temperature, 100001)
    sigma = 2.5e-20
    law = {'kind': 'ElectronCollisionPlasma', 'energy_J_particle': energy,
           'sigma_m2': np.full_like(energy, sigma)}
    actual = benchmark.independent_rate(law, temperature, 298.15, c)
    exact = sigma * math.sqrt(8 * c['kB'] * temperature / (math.pi * c['m_e']))
    assert actual == pytest.approx(exact, rel=3e-8)


def test_two_temperature_equal_barriers_cancel_gas_temperature():
    c = {'R': 8.314472}
    law = {'kind': 'TwoTemperaturePlasma', 'A_particle': 5e-15,
           'n': .74, 'T0_K': 11604.51812, 'Ea_g_J_mol': 100000., 'Ea_e_J_mol': 100000.}
    te = 12000.
    expected = law['A_particle'] * (te / law['T0_K']) ** .74 * math.exp(-100000 / (c['R'] * te))
    for tg in [298.15, 1000, 10000]:
        assert benchmark.independent_rate(law, te, tg, c) == pytest.approx(expected, rel=1e-12)


def test_neutral_termolecular_law_uses_gas_temperature():
    law = {'kind': 'Arrhenius', 'A_particle': 1.1e-43, 'n': 0,
           'T0_K': 1., 'Ea_J_mol': 0.}
    assert benchmark.independent_rate(law, 11604., 298.15, {'R': 8.314472}) == 1.1e-43


def test_sweep_contains_every_point_in_both_directions():
    first, second = benchmark.sweep_coordinates()
    assert first[0] == 'ascending' and second[0] == 'descending'
    assert len(first[1]) == 12 and len(set(first[1])) == 12
    assert second[1] == list(reversed(first[1]))
    assert set(first[1]) == {(p, w) for p in (1., 5., 20.) for w in (.01, .1, .5, 5.)}


def test_comparison_refuses_partial_outputs():
    with pytest.raises(ValueError, match='demonstrated steady state'):
        benchmark.compare_point({'status': 'not_steady'}, {})
    with pytest.raises(ValueError, match='power_relative_residual'):
        benchmark.compare_point({'status': 'steady'}, {})


def test_comparison_reports_real_mismatch_instead_of_accepting_it():
    import rmgpy.constants as c
    rmg = {'status': 'steady', 'densities_m3': {'e-': 2e16, 'Arp': 2e16, 'Ars': 1e15},
           'volume_m3': 1., 'budget_W_inventory': {'Te': c.e / (c.R / c.Na),
           'P_abs': 10., 'Q_inelastic': 4., 'Q_elastic': 4.,
           'Q_wall_electron': 1., 'Q_wall_ion': 1.},
           'runtime_seconds': 1., 'physical_time_seconds': .1}
    loki = {'power_relative_residual': 0., 'charge_relative_residual': 0.,
            'eos_relative_residual': 0., 'particle_relative_residual': 0.,
            'ne_m3': 1e16, 'Te_eV': 1., 'Arp_m3': 1e16, 'Ars_m3': 1e15,
            'P_abs_W_m3': 10., 'reaction_power_W_m3': [4.],
            'consumed_electron_power_W_m3': [0.], 'elastic_W_m3': 4.,
            'wall_electron_W_m3': 1., 'wall_ion_W_m3': 1.,
            'runtime_total_seconds': 5., 'outer_iterations': 3}
    result = benchmark.compare_point(rmg, loki)
    assert result['acceptance'] is False
    assert result['rows'][0]['relative_difference'] == 1.
    assert result['loki_physical_settling_seconds'] is None
    loki['ne_m3'] = 2e16
    assert benchmark.compare_point(rmg, loki)['acceptance'] is True
    loki['particle_relative_residual'] = .01
    with pytest.raises(ValueError, match='particle_relative_residual'):
        benchmark.compare_point(rmg, loki)


def test_reaction_comparison_aligns_ids_and_preserves_signed_and_consumed_power():
    r_manifest = {'reactions': [{'id': 'superelastic'}, {'id': 'recombination'}]}
    l_manifest = {'reactions': list(reversed(r_manifest['reactions']))}
    rmg = {'volume_m3': 2., 'reaction_rates_m3_s': [3., 4.],
           'budget_W_inventory': {'Q_inelastic_by_reaction': [-10., 6.]}}
    loki = {'source_rate_m3_s': [4., 3.], 'reaction_power_W_m3': [0., -5.],
            'consumed_electron_power_W_m3': [3., 0.]}
    rows = benchmark.compare_reactions(rmg, loki, r_manifest, l_manifest)
    assert [row['reaction_id'] for row in rows] == ['superelastic', 'recombination']
    assert all(row['rate_relative_difference'] == 0 for row in rows)
    assert all(row['power_relative_difference'] == 0 for row in rows)
    assert rows[0]['power_rmg'] == -5.
    assert rows[1]['power_loki'] == 3.
    with pytest.raises(ValueError, match='manifests differ'):
        benchmark.compare_reactions(rmg, loki, r_manifest, {'reactions': [{'id': 'wrong'}]})
    loki['source_rate_m3_s'] = [7.]
    with pytest.raises(ValueError, match='vector length'):
        benchmark.compare_reactions(rmg, loki, r_manifest, l_manifest)


def test_zero_reference_never_hides_nonzero_reaction_difference():
    assert benchmark.comparison_value(0., 0.)['relative_difference'] == 0.
    difference = benchmark.comparison_value(3., 0.)
    assert difference['relative_difference'] is None
    assert difference['absolute_difference'] == 3.
    with pytest.raises(ValueError, match='Nonfinite'):
        benchmark.comparison_value(float('nan'), 0.)
