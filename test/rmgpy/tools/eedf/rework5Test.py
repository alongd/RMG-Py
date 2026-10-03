#!/usr/bin/env python3

from pathlib import Path
import runpy

import numpy as np
import pytest

from rmgpy.tools.eedf.integrity import physical_properties, resolved_properties
from rmgpy.tools.eedf import moments, validation


helpers = runpy.run_path(str(Path(__file__).with_name('rework3Test.py')))
request_and_channels = helpers['request_and_channels']


def test_r5_regression_saved_real_forward_rate_with_inactive_target_qualifies(monkeypatch):
    """H3 keeps the relative scale of a saved real forward coefficient.

    The 12 values are the inactive-target channel added to the rework-4
    real-LoKI held-out output.  Their smallest relative interpolation error is
    1.34e-7; zero population must not reduce the allowance to 1e-40.
    """
    saved_rates = np.array([2.69258250e-21, 1.20303373e-21, 1.01689287e-22,
                            4.06742523e-22, 1.36608582e-25, 2.34630716e-26,
                            7.87889020e-26, 4.81939703e-27, 4.28326045e-27,
                            5.88909685e-27, 1.84170406e-27, 1.65182922e-27])
    powers = {'elastic_gain': 0., 'elastic_loss': 0., 'car_gain': 0., 'car_loss': 0.,
              'excitation_loss': 0., 'excitation_gain': 0., 'vibrational_loss': 0.,
              'vibrational_gain': 0., 'rotational_loss': 0., 'rotational_gain': 0.,
              'ionization': 0., 'attachment': 0., 'field': 1., 'growth': 0.}
    manifest = {'tolerances': {'H1': {'atol': 0., 'rtol': 1e-6},
                               'H2': {'atol': 0., 'rtol': 1e-6},
                               'H3': {'atol': 1e-40, 'rtol': 1e-6, 'total_rtol': 1e-6},
                               'H4': {'atol': 0., 'rtol': 1e-6, 'total_rtol': 1e-6},
                               'F0': {'atol': 0.},
                               'channel_power_sum': {'atol': 0., 'rtol': 1e-6}},
                'floors': {'rate_flux_fraction': 0., 'absolute_flux_fraction': 0.,
                           'relative_power_share': 0., 'absolute_power_share': 0.},
                'channel_map': [{'kind': 'excitation', 'description': 'e + Ar(3P2) -> Ar(4p), Excitation'}],
                'row_inputs': {'Tg_K': 300.}}
    monkeypatch.setattr(moments, 'distribution_moments',
                        lambda predicted, channels: {'mean_energy_eV': 1., 'k_ine': predicted['k_ine'],
                                                     'k_sup': predicted['k_sup'], 'channel_power': predicted['channel_power']})
    for rate in saved_rates:
        direct = {'energy_eV': np.array([1.]), 'energy_edges_eV': np.array([0., 2.]),
                  'swarm': {'mean_energy_eV': 1., 'townsend_N': 1.}, 'k_ine': np.array([rate]),
                  'k_sup': np.array([0.]), 'target_fractions': np.array([0.]),
                  'product_fractions': np.array([0.]), 'below_floor': np.array([[False, True]]),
                  'rate_floors': np.array([1e-40]), 'channel_power': np.array([0.]),
                  'power_groups': powers, 'composition': {'Tg_K': 300.},
                  'f0': np.array([0.]), 'attachment_energy_eV': 0.}
        predicted = dict(direct, k_ine=np.array([rate * 1.000000134]))
        checks = [check for check in validation.compare_row(predicted, direct, manifest)
                  if check['rule'] == 'H3' and check['quantity'] == 'k_ine:0']
        assert checks and checks[0]['passed'], checks
        predicted['k_ine'] *= 1.1
        checks = [check for check in validation.compare_row(predicted, direct, manifest)
                  if check['rule'] == 'H3' and check['quantity'] == 'k_ine:0']
        assert checks and not checks[0]['passed']


@pytest.mark.parametrize('property_name', ['energy', 'statisticalWeight'])
@pytest.mark.parametrize('alias', ['Ar(1S0,)', 'Ar(,1S0)', 'Ar(,1S0,)'])
def test_r5_regression_identical_selector_aliases_are_accepted(tmp_path, property_name, alias):
    spec, _ = request_and_channels(tmp_path)
    spec['state_properties'][property_name] = ['Ar(1S0) = 1', alias + ' = 1']
    _, properties = physical_properties(spec)
    assert properties[property_name] == ['Ar(1S0) = 1']
    _, resolved, _, _ = resolved_properties(spec, {})
    assert resolved[property_name] == ['Ar(1S0) = 1']
