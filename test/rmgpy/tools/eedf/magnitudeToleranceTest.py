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

"""Hand-calculated impact budgets; no solver or generated table required."""
import copy
import json
from pathlib import Path

import numpy as np
import pytest

from rmgpy.tools.eedf import moments, validation
from rmgpy.tools.eedf.constants import BOLTZMANN
from rmgpy.tools.eedf.schema import (
    _content_hash, _row_inputs, _validate_magnitude_tolerances, SpecError,
)


@pytest.fixture
def comparison(monkeypatch):
    policy = json.loads((Path(__file__).parents[4] /
                         'scripts/eedf_field_drive/magnitude-tolerances.json').read_text())
    # N=1e20, ne_max=1e14, n_floor=1e8, t=10, x=.5:
    # k_atol=1e-6*1e8/(1e14*.5*1e20*10)=2e-33.
    policy['time_max_s'] = 10.
    tolerances = {name: {'rtol': .01, 'atol': 0.} for name in ('H1', 'H2', 'H3', 'H4')}
    tolerances.update(F0={'atol': 0.}, channel_power_sum={'rtol': .01, 'atol': 0.})
    tolerances['H4']['total_rtol'] = .01
    manifest = dict(
        tolerances=tolerances, floors=dict(absolute_flux_fraction=0., relative_power_share=0.,
                                          absolute_power_share=0.),
        channel_map=[{'kind': 'excitation'}, {'kind': 'elastic'}],
        row_inputs=dict(Tg_K=300., P_Pa=1e20*BOLTZMANN*300., magnitude_tolerances=policy))
    powers = {name: 0. for name in (
        'elastic_gain', 'elastic_loss', 'car_gain', 'car_loss', 'excitation_loss',
        'excitation_gain', 'vibrational_loss', 'vibrational_gain', 'rotational_loss',
        'rotational_gain', 'ionization', 'attachment', 'growth')}
    powers['field'] = 2e-20
    row = dict(energy_eV=np.array([1.]), energy_edges_eV=np.array([0., 2.]),
               swarm=dict(mean_energy_eV=1., townsend_N=0., drift_velocity_m_s=2.),
               k_ine=np.zeros(2), k_sup=np.zeros(2), target_fractions=np.array([.5, .5]),
               product_fractions=np.array([.25, .25]), below_floor=np.zeros((2, 2), dtype=bool),
               rate_floors=np.zeros(2), channel_power=np.zeros(2), power_groups=powers,
               f0=np.zeros(1), attachment_energy_eV=0.)
    # Fixed independent moment references. Differences may therefore occur in
    # BOTH H3 and consistency, and must be charged twice.
    expected = dict(mean_energy_eV=1., k_ine=np.zeros(2), k_sup=np.zeros(2),
                    channel_power=np.zeros(2))
    monkeypatch.setattr(moments, '_distribution_moments', lambda *args: expected)
    return copy.deepcopy(row), row, manifest, expected


def select(checks, rule, quantity, index=None):
    found = [c for c in checks if c['rule'] == rule and c['quantity'] == quantity
             and (index is None or c['index'] == index)]
    assert found
    return found[0]


def test_negligible_rate_passes_through_atol(comparison):
    predicted, direct, manifest, _ = comparison
    predicted['k_ine'][0] = 8e-34
    checks = validation._compare_row(predicted, direct, manifest)
    rate = select(checks, 'H3', 'k_ine:0')
    assert rate['allowed_error'] == pytest.approx(2e-33, rel=1e-14, abs=0)
    assert rate['absolute_error'] == 8e-34 and rate['passed']
    # Each of two constraints spends 4e-7 of inventory; total 8e-7.
    cumulative = select(checks, 'cumulative', 'species_fraction')
    assert cumulative['absolute_error'] == pytest.approx(8e-7, rel=1e-14, abs=0)
    assert all(c['passed'] for c in checks)


def test_invented_non_negligible_rate_fails_zero_reference(comparison):
    predicted, direct, manifest, _ = comparison
    predicted['k_ine'][0] = 3e-33
    c = select(validation._compare_row(predicted, direct, manifest), 'H3', 'k_ine:0')
    assert not c['passed'] and c['relative_error'] is None
    assert c['allowed_error'] == pytest.approx(2e-33, rel=1e-14, abs=0)


def test_cumulative_query_fails_when_every_individual_check_passes(comparison):
    predicted, direct, manifest, _ = comparison
    predicted['k_ine'][0] = 1.5e-33
    checks = validation._compare_row(predicted, direct, manifest)
    assert all(c['passed'] for c in checks if c['rule'] != 'cumulative')
    c = select(checks, 'cumulative', 'species_fraction')
    assert c['absolute_error'] == pytest.approx(1.5e-6, rel=1e-14, abs=0)
    assert not c['passed']


def test_cumulative_charges_excess_over_relative_allowance(comparison):
    predicted, direct, manifest, expected = comparison
    direct['k_ine'][0] = expected['k_ine'][0] = 1e-30
    predicted['k_ine'][0] = 1.0108e-30
    checks = validation._compare_row(predicted, direct, manifest)
    # d=1.08e-32, rtol*ref=1e-32, excess=8e-34, counted twice.
    c = select(checks, 'cumulative', 'species_fraction')
    assert c['absolute_error'] == pytest.approx(8e-7, rel=1e-12, abs=0)
    assert c['passed']


@pytest.mark.parametrize('floor_kind', ['rate', 'flux', 'absolute'])
def test_old_rate_allowance_passes_without_budget_charge(comparison, floor_kind):
    predicted, direct, manifest, expected = comparison
    direct['k_ine'][0] = expected['k_ine'][0] = 1e-30
    predicted['k_ine'][0] = 1.0108e-30
    if floor_kind == 'rate':
        predicted['rate_floors'][0] = direct['rate_floors'][0] = 2e-32
        direct['below_floor'][0, 0] = True
    elif floor_kind == 'flux':
        manifest['floors']['absolute_flux_fraction'] = .02
        # Isolate the H3 flux-floor comparison from the consistency check.
        expected['k_ine'][0] = predicted['k_ine'][0]
    else:
        manifest['tolerances']['H3']['atol'] = 2e-32
    old_manifest = copy.deepcopy(manifest)
    old_manifest['row_inputs'].pop('magnitude_tolerances')
    old = validation._compare_row(predicted, direct, old_manifest)
    new = validation._compare_row(predicted, direct, manifest)
    assert all(c['passed'] for c in old)
    assert all(c['passed'] for c in new)
    assert select(new, 'H3', 'k_ine:0') == select(old, 'H3', 'k_ine:0')
    assert select(new, 'cumulative', 'species_fraction')['absolute_error'] == 0.


def test_old_power_floor_passes_without_budget_charge(comparison):
    predicted, direct, manifest, _ = comparison
    manifest['floors'].update(absolute_power_share=1e-4, relative_power_share=1e-3)
    predicted['channel_power'][0] = 1e-24
    old_manifest = copy.deepcopy(manifest)
    old_manifest['row_inputs'].pop('magnitude_tolerances')
    old = validation._compare_row(predicted, direct, old_manifest)
    new = validation._compare_row(predicted, direct, manifest)
    # Old floor=2e-24 dominates impact floor=2e-29, including zero-rtol H4.
    for rule, quantity, index in [('H4', 'channel_power:0', []),
                                  ('consistency', 'f0.channel_power', [0])]:
        before = select(old, rule, quantity, index)
        assert before['passed']
        assert select(new, rule, quantity, index) == before
    assert select(new, 'cumulative', 'energy_fraction')['absolute_error'] == 0.


def test_new_allowance_charges_only_excess_over_old_atol_and_rtol(comparison):
    predicted, direct, manifest, expected = comparison
    direct['k_ine'][0] = expected['k_ine'][0] = 1e-30
    predicted['k_ine'][0] = 1.0108e-30
    predicted['rate_floors'][0] = direct['rate_floors'][0] = 4e-34
    old_manifest = copy.deepcopy(manifest)
    old_manifest['row_inputs'].pop('magnitude_tolerances')
    old = validation._compare_row(predicted, direct, old_manifest)
    new = validation._compare_row(predicted, direct, manifest)
    assert not select(old, 'H3', 'k_ine:0')['passed']
    assert all(c['passed'] for c in new)
    # d=1.08e-32, old allowance=1.04e-32; excess=4e-34.
    # Two independent constraints, each scaled by 5e26, spend 4e-7 total.
    c = select(new, 'cumulative', 'species_fraction')
    assert c['absolute_error'] == pytest.approx(4e-7, rel=1e-12, abs=0)
    assert select(new, 'H3', 'k_ine:0')['allowed_error'] == pytest.approx(1.2e-32, rel=1e-14, abs=0)


def test_townsend_and_power_formulas_and_energy_cumulative(comparison):
    predicted, direct, manifest, _ = comparison
    # alpha_atol=1e-6/(2*1e20*10)=5e-28; P_atol=1e-9*2e-20=2e-29.
    predicted['swarm']['townsend_N'] = 4e-28
    predicted['channel_power'][0] = 1.5e-29
    checks = validation._compare_row(predicted, direct, manifest)
    c = select(checks, 'H2', 'townsend_N')
    assert c['allowed_error'] == pytest.approx(5e-28, rel=1e-14, abs=0) and c['passed']
    c = select(checks, 'H4', 'channel_power:0')
    assert c['allowed_error'] == pytest.approx(2e-29, rel=1e-14, abs=0) and c['passed']
    assert select(checks, 'cumulative', 'electron_fraction')['absolute_error'] == pytest.approx(8e-7)
    c = select(checks, 'cumulative', 'energy_fraction')
    assert c['absolute_error'] == pytest.approx(1.5e-9, rel=1e-14, abs=0)
    assert not c['passed']


def test_subnormal_is_zero_before_comparison(comparison):
    predicted, direct, manifest, _ = comparison
    predicted['swarm']['townsend_N'] = 9.4e-323
    direct['swarm']['townsend_N'] = 5.924e-321
    c = select(validation._compare_row(predicted, direct, manifest), 'H2', 'townsend_N')
    assert c['absolute_error'] == 0 and c['relative_error'] is None and c['passed']


def test_larger_target_and_reverse_target_floor(comparison):
    predicted, direct, manifest, _ = comparison
    predicted['target_fractions'][0] = 1.
    direct['product_fractions'][0] = predicted['product_fractions'][0] = 0.
    checks = validation._compare_row(predicted, direct, manifest)
    assert select(checks, 'H3', 'k_ine:0')['allowed_error'] == pytest.approx(1e-33, rel=1e-14, abs=0)
    assert select(checks, 'H3', 'k_sup:0')['allowed_error'] == pytest.approx(1e-21, rel=1e-14, abs=0)


def test_elastic_and_aggregate_gates_remain_exact(comparison):
    predicted, direct, manifest, _ = comparison
    predicted['k_ine'][1] = 2e-34
    predicted['channel_power'][1] = 1e-30
    new = validation._compare_row(predicted, direct, manifest)
    manifest['row_inputs'].pop('magnitude_tolerances')
    old = validation._compare_row(predicted, direct, manifest)
    for c in old:
        if (c['rule'] in ('H1', 'F0', 'channel_power_sum') or
                c['quantity'] in ('k_ine:1', 'k_sup:1', 'channel_power:1') or
                (c['rule'] == 'consistency' and c['index'] == [1]) or
                c['quantity'].startswith('power.')):
            assert c in new
    assert not select(new, 'H3', 'k_ine:1')['passed']


@pytest.mark.parametrize('key', ['species_epsilon', 'electron_epsilon', 'energy_epsilon',
    'neutral_floor_fraction', 'electron_fraction_max', 'time_max_s',
    'target_fraction_floor', 'underflow_absolute', 'time_sources'])
def test_policy_is_validated_and_fingerprinted(comparison, key):
    _, _, manifest, _ = comparison
    policy = manifest['row_inputs']['magnitude_tolerances']
    _validate_magnitude_tolerances(policy)
    # _row_inputs needs only these explicit fields for identity assembly.
    spec = {k: {} for k in ('arm', 'Tg_K', 'P_Pa', 'loki_commit', 'axes', 'gas_properties',
        'state_properties', 'working_conditions', 'solver_options', 'screen',
        'envelopes', 'floors', 'tolerances', 'held_out', 'refinement')}
    spec['magnitude_tolerances'] = policy
    before = _content_hash(_row_inputs(spec))
    policy[key] = ['different horizon source'] if key == 'time_sources' else policy[key] * 2
    assert _content_hash(_row_inputs(spec)) != before
    policy[key] = [] if key == 'time_sources' else float('nan')
    with pytest.raises(SpecError, match='magnitude_tolerances'):
        _row_inputs(spec)


def test_query_aggregation_consumes_cumulative_failure(comparison):
    from types import SimpleNamespace
    predicted, direct, manifest, _ = comparison
    predicted['k_ine'][0] = 1.5e-33
    predicted['branch_id'] = 'branch_0'
    manifest['fingerprint'] = 'test-identity'
    predicted['swarm']['mobility_N'] = direct['swarm']['mobility_N'] = 1.
    table = SimpleNamespace(axis_names=('u', 'f'), manifest=manifest,
                            row=lambda *args: SimpleNamespace(_as_dict=lambda: predicted))
    result = validation._held_out_accuracy(table, direct, dict(u=0., f=.5),
        criteria={k: dict(rtol=.01, atol=0.) for k in ('mobility_N', 'field_power')},
        setup_fingerprint='test-identity')
    assert not result['passed']
    assert all(c['passed'] for c in result['checks'] if c['rule'] != 'cumulative')


def test_cutoff_is_strict_and_policy_absence_keeps_subnormals(comparison):
    predicted, direct, manifest, _ = comparison
    predicted['swarm']['townsend_N'] = 1e-300
    c = select(validation._compare_row(predicted, direct, manifest), 'H2', 'townsend_N')
    assert c['absolute_error'] == 1e-300
    predicted['swarm']['townsend_N'] = 5e-324
    manifest['row_inputs'].pop('magnitude_tolerances')
    c = select(validation._compare_row(predicted, direct, manifest), 'H2', 'townsend_N')
    assert c['absolute_error'] == 5e-324 and not c['passed']


def test_saved_c21_channel16_94_percent_error_still_fails(comparison):
    """c21 refine-1: u=1.0737477931164006, f=1.6666666666666667e-8."""
    predicted, direct, manifest, expected = comparison
    manifest['row_inputs'].update(P_Pa=666.61184, Tg_K=298.15)
    manifest['row_inputs']['magnitude_tolerances']['time_max_s'] = 600.
    manifest['tolerances']['H3']['rtol'] = .005
    predicted['target_fractions'][0] = direct['target_fractions'][0] = .9999999833333336
    predicted['k_ine'][0] = 2.685286677336356e-26
    direct['k_ine'][0] = expected['k_ine'][0] = 4.237540123631974e-25
    c = select(validation._compare_row(predicted, direct, manifest),
               'consistency', 'f0.k_ine', [0])
    # N=1.6194394159529672e23, k_atol=1.02916273e-38.
    # .005*4.237540123631974e-25 + k_atol = 2.118770061826279e-27.
    assert c['allowed_error'] == pytest.approx(2.118770061826279e-27, rel=1e-14, abs=0)
    assert c['absolute_error'] == pytest.approx(3.9690114558983384e-25, rel=1e-14, abs=0)
    assert c['relative_error'] == pytest.approx(.9366310029169751)
    assert not c['passed']


def test_saved_cell025_elastic_sum_error_still_fails(comparison):
    """Frozen predicted power scalars at u=-5.639666005844986, f=0."""
    predicted, direct, manifest, expected = comparison
    manifest['floors'].update(absolute_power_share=1e-4, relative_power_share=1e-3)
    manifest['tolerances']['channel_power_sum'] = dict(atol=1e-24, rtol=.005)
    for row in (predicted, direct):
        row['power_groups'].update(field=2.1505102859181662e-21,
                                  elastic_gain=1.5152498869282245e-21,
                                  elastic_loss=-3.665267010618173e-21)
        row['channel_power'][1] = 2.1504943003650167e-21
    expected['channel_power'][1] = predicted['channel_power'][1]
    checks = validation._compare_row(predicted, direct, manifest)
    c = select(checks, 'channel_power_sum', 'power.sum.elastic')
    # Group=2.1500171236899486e-21, error=4.771766750680394e-25.
    # allowed=(1e-4 + .005*1e-3)*field=2.2580358002140746e-25.
    assert c['absolute_error'] == pytest.approx(4.771766750680394e-25, rel=1e-14, abs=0)
    assert c['allowed_error'] == pytest.approx(2.2580358002140746e-25, rel=1e-14, abs=0)
    assert not c['passed']


def test_full_spec_admits_policy_and_rejects_missing_budget():
    import os
    from rmgpy.tools.eedf.schema import _validate_spec
    spec = json.loads(Path(os.environ['EEDF_REAL_SPEC']).read_text())
    spec['magnitude_tolerances'] = json.loads((Path(__file__).parents[4] /
                         'scripts/eedf_field_drive/magnitude-tolerances.json').read_text())
    _validate_spec(spec)
    del spec['magnitude_tolerances']['species_epsilon']
    with pytest.raises(SpecError, match='magnitude_tolerances'):
        _validate_spec(spec)


@pytest.mark.parametrize('coordinate', ['P_Pa', 'Tg_K'])
@pytest.mark.parametrize('source', ['axes', 'envelope', 'envelope_query'])
def test_query_density_controls_allowances_and_cumulative_charge(comparison, coordinate, source):
    predicted, direct, manifest, _ = comparison
    inputs = manifest['row_inputs']
    inputs['magnitude_tolerances']['time_max_s'] = 600.
    nominal = inputs[coordinate]
    query = nominal * (10 if coordinate == 'P_Pa' else .1)
    if source == 'axes':
        predicted['composition'] = direct['composition'] = {coordinate: query}
    else:
        inputs['envelopes'] = {coordinate: {'reference': query}}
        if source == 'envelope_query':
            inputs['envelopes'][coordinate]['reference'] = nominal
            direct['composition'] = {coordinate: query}
    # Reviewer's exact witness: nominal N=1e20 admits both errors, whereas
    # query N=1e21 rejects them individually before cumulative charging.
    predicted['k_ine'][0] = 1e-35
    predicted['swarm']['townsend_N'] = 5e-30
    checks = validation._compare_row(predicted, direct, manifest)
    rate = select(checks, 'H3', 'k_ine:0')
    townsend = select(checks, 'H2', 'townsend_N')
    assert rate['allowed_error'] == pytest.approx(1e-6 / 3e29, rel=1e-14, abs=0)
    assert townsend['allowed_error'] == pytest.approx(1e-6 / 1.2e24, rel=1e-14, abs=0)
    assert not rate['passed'] and not townsend['passed']
    # Reduce each error to below its query allowance: two independent rate
    # constraints still exceed the query's cumulative inventory limit.
    predicted['k_ine'][0] = 2e-36
    predicted['swarm']['townsend_N'] = 5e-31
    checks = validation._compare_row(predicted, direct, manifest)
    assert all(c['passed'] for c in checks if c['rule'] != 'cumulative')
    species = select(checks, 'cumulative', 'species_fraction')
    assert species['absolute_error'] == pytest.approx(1.2e-6, rel=1e-14, abs=0)
    assert not species['passed']
    electron = select(checks, 'cumulative', 'electron_fraction')
    assert electron['absolute_error'] == pytest.approx(6e-7, rel=1e-14, abs=0)


@pytest.mark.parametrize('zero_drift,zero_power', [(True, False), (False, True), (True, True)])
def test_zero_reference_scales_still_normalize_underflow(comparison, zero_drift, zero_power):
    predicted, direct, manifest, _ = comparison
    for row in (predicted, direct):
        if zero_drift:
            row['swarm']['drift_velocity_m_s'] = 0.
        if zero_power:
            row['power_groups']['field'] = 0.
    checks = validation._compare_row(predicted, direct, manifest)
    assert all(c['passed'] for c in checks)
    assert all(c['absolute_error'] == 0 for c in checks if c['rule'] == 'cumulative')
    # Zero impact allowance still admits values that both normalize to zero.
    if zero_drift:
        predicted['swarm']['townsend_N'] = 5e-324
    if zero_power:
        predicted['channel_power'][0] = 5e-324
    new = validation._compare_row(predicted, direct, manifest)
    old_manifest = copy.deepcopy(manifest)
    old_manifest['row_inputs'].pop('magnitude_tolerances')
    old = validation._compare_row(predicted, direct, old_manifest)
    keys = []
    if zero_drift:
        keys.append(('H2', 'townsend_N', []))
    if zero_power:
        keys.extend([('H4', 'channel_power:0', []), ('consistency', 'f0.channel_power', [0])])
    for rule, quantity, index in keys:
        assert not select(old, rule, quantity, index)['passed']
        check = select(new, rule, quantity, index)
        assert check['passed']
        assert check['absolute_error'] == check['allowed_error'] == 0.
    assert all(c['absolute_error'] == 0 for c in new if c['rule'] == 'cumulative')
    # Values at the cutoff are not normalized and gain no impact allowance.
    if zero_drift:
        predicted['swarm']['townsend_N'] = 1e-300
    if zero_power:
        predicted['channel_power'][0] = 1e-300
    boundary = validation._compare_row(predicted, direct, manifest)
    for rule, quantity, index in keys:
        assert not select(boundary, rule, quantity, index)['passed']


@pytest.mark.parametrize('value', [float('nan'), float('inf'), -1.])
@pytest.mark.parametrize('scale', ['drift', 'field', 'P_Pa', 'Tg_K'])
def test_invalid_reference_scales_still_raise(comparison, scale, value):
    predicted, direct, manifest, _ = comparison
    if scale == 'drift':
        direct['swarm']['drift_velocity_m_s'] = value
    elif scale == 'field':
        direct['power_groups']['field'] = value
    else:
        direct['composition'] = {scale: value}
    with pytest.raises(validation.EEDFError, match='reference scales'):
        validation._compare_row(predicted, direct, manifest)


def test_underflow_boundary_preserves_raw_old_power_pass(comparison):
    predicted, direct, manifest, expected = comparison
    predicted['channel_power'][0] = 1e-300
    direct['channel_power'][0] = expected['channel_power'][0] = 9.99e-301
    predicted['power_groups']['field'] = direct['power_groups']['field'] = 1e-300
    old_manifest = copy.deepcopy(manifest)
    old_manifest['row_inputs'].pop('magnitude_tolerances')
    old = validation._compare_row(predicted, direct, old_manifest)
    new = validation._compare_row(predicted, direct, manifest)
    for rule, quantity, index in [('H4', 'channel_power:0', []),
                                  ('consistency', 'f0.channel_power', [0])]:
        before = select(old, rule, quantity, index)
        assert before['passed']
        assert select(new, rule, quantity, index) == before
    assert select(new, 'cumulative', 'energy_fraction')['absolute_error'] == 0.


@pytest.mark.parametrize('kind', ['elastic', 'effective'])
def test_momentum_transfer_channels_keep_exact_old_checks(comparison, kind):
    predicted, direct, manifest, _ = comparison
    manifest['channel_map'][0]['kind'] = kind
    manifest['row_inputs']['magnitude_tolerances']['time_max_s'] = 600.
    predicted['k_ine'][0] = predicted['k_sup'][0] = 1e-35
    predicted['channel_power'][0] = 4e-30
    new = validation._compare_row(predicted, direct, manifest)
    manifest['row_inputs'].pop('magnitude_tolerances')
    old = validation._compare_row(predicted, direct, manifest)
    assert [c for c in new if c['rule'] != 'cumulative' and c['quantity'] != 'townsend_N'] == [
        c for c in old if c['quantity'] != 'townsend_N']
    assert not select(new, 'H3', 'k_ine:0')['passed']
    assert not select(new, 'H3', 'k_sup:0')['passed']
    assert not select(new, 'H4', 'channel_power:0')['passed']
    assert all(c['absolute_error'] == 0 for c in new if c['rule'] == 'cumulative')


def test_normalized_power_uses_larger_old_allowance(comparison):
    predicted, direct, manifest, expected = comparison
    predicted['channel_power'][0] = -9e-301
    direct['channel_power'][0] = expected['channel_power'][0] = 1e-300
    predicted['power_groups']['field'] = direct['power_groups']['field'] = 1e-300
    manifest['floors']['absolute_power_share'] = 1.
    manifest['tolerances']['channel_power_sum']['atol'] = 1e-300
    old_manifest = copy.deepcopy(manifest)
    old_manifest['row_inputs'].pop('magnitude_tolerances')
    old = validation._compare_row(predicted, direct, old_manifest)
    new = validation._compare_row(predicted, direct, manifest)
    for rule, quantity, index in [('H4', 'channel_power:0', []),
                                  ('consistency', 'f0.channel_power', [0])]:
        assert not select(old, rule, quantity, index)['passed']
        check = select(new, rule, quantity, index)
        assert check['passed']
        assert check['absolute_error'] == 1e-300
        assert check['allowed_error'] == pytest.approx(1.01e-300, rel=1e-14, abs=0)
    assert all(c['passed'] for c in new)
    assert select(new, 'cumulative', 'energy_fraction')['absolute_error'] == 0.


def test_normalized_power_charge_uses_normalized_old_reference(comparison):
    predicted, direct, manifest, expected = comparison
    predicted['channel_power'][0] = 1.1e-300
    direct['channel_power'][0] = expected['channel_power'][0] = 9e-301
    predicted['power_groups']['field'] = direct['power_groups']['field'] = 2e-290
    checks = validation._compare_row(predicted, direct, manifest)
    for rule, quantity, index in [('H4', 'channel_power:0', []),
                                  ('consistency', 'f0.channel_power', [0])]:
        check = select(checks, rule, quantity, index)
        assert check['passed']
        assert check['absolute_error'] == 1.1e-300
    # The normalized reference and its old relative allowance are zero.
    # Both constraints charge (1.1e-300 - 0) / 2e-290.
    cumulative = select(checks, 'cumulative', 'energy_fraction')
    assert cumulative['absolute_error'] == pytest.approx(1.1e-10, rel=1e-14, abs=0)
    assert cumulative['passed']
