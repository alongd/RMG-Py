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

"""Owner-ruling regressions for the named O2 unity wall closure."""

import json
import logging
import pickle
from pathlib import Path

import numpy as np
import pytest

import rmgpy.constants as constants
from plasmaElectronegativeWallTest import model, qualification
from rmgpy.exceptions import PlasmaStateError
from rmgpy.solver.electronegative import (
    closure_factor_gradient,
    qualify_envelope,
    wall_discrepancy,
)
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.species import Species
from rmgpy.thermo import ThermoData


UNITY = 'o2ReferenceQualifiedUnity'
ORACLE = Path(__file__).with_name('data') / 'o2_unity_reference.json'
logger = logging.getLogger(__name__)

REFERENCE_RADIUS = 0.05
REFERENCE_RADIAL_EIGENVALUE = (2.404825557695773 / REFERENCE_RADIUS) ** 2
REFERENCE_AXIAL_EIGENVALUE = REFERENCE_RADIAL_EIGENVALUE * 1.e-8
O2_REDUCED_MOBILITIES = {
    'Ar+': 1.535e-4,
    'O+': 3.2047284103791057e-4,
    'O2+': 2.5697516730204976e-4,
}


def _test_thermo():
    return ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([20.8] * 7, 'J/(mol*K)'),
        H298=(0., 'kJ/mol'), S298=(150., 'J/(mol*K)'))


def _reference_geometry_o2_engine(point):
    """Build the production operator at one frozen I-314 engine-side state."""
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    arp = Species(label='Ar+').from_adjacency_list('1 Ar u1 p3 c+1')
    oxygen = Species(label='O').from_adjacency_list('1 O u2 p2 c0')
    op = Species(label='O+').from_adjacency_list('1 O u3 p1 c+1')
    o2 = Species(label='O2').from_adjacency_list(
        'multiplicity 3\n1 O u1 p2 c0 {2,S}\n2 O u1 p2 c0 {1,S}')
    o2p = Species(label='O2+').from_adjacency_list(
        'multiplicity 2\n1 O u0 p2 c0 {2,D}\n2 O u1 p1 c+1 {1,D}')
    om = Species(label='O-').from_adjacency_list('1 O u1 p3 c-1')
    species = [electron, ar, arp, oxygen, op, o2, o2p, om]
    for item in species[1:]:
        item.thermo = _test_thermo()

    frozen = point['state']
    cation_densities = point['engine_cation_number_densities_m_3']
    initial_electron_density = sum(cation_densities)
    initial_densities = {
        electron: initial_electron_density,
        ar: frozen['neutral_density_m_3'],
        arp: cation_densities[0],
        oxygen: 0.,
        op: cation_densities[1],
        o2: 0.,
        o2p: cation_densities[2],
        om: 0.,
    }
    total_density = sum(initial_densities.values())
    amounts = {
        item: density / total_density
        for item, density in initial_densities.items()}
    electron_temperature = (
        frozen['transport_energy_eV'] * constants.e / constants.kB)
    anion_density = initial_electron_density - frozen['electron_density_m_3']
    pressure = constants.kB * (
        frozen['electron_density_m_3'] * electron_temperature
        + (frozen['neutral_density_m_3'] + initial_electron_density
           + anion_density) * 298.15)
    reactor = PlasmaReactor(
        (298.15, 'K'), (pressure, 'Pa'), amounts,
        (electron_temperature, 'K'),
        diffusion_length=(
            1. / np.sqrt(REFERENCE_RADIAL_EIGENVALUE + REFERENCE_AXIAL_EIGENVALUE),
            'm'),
        ion_reduced_mobilities={
            label: (mobility, 'm^2/(V*s)')
            for label, mobility in O2_REDUCED_MOBILITIES.items()},
        mobility_reference_density=(2.6867811e25, 'm^-3'),
        ambipolar_ion_temperature='gas',
        anion_reduced_mobilities={'O-': (1.5e-4, 'm^2/(V*s)')},
        wall_diffusion_components=(
            REFERENCE_RADIAL_EIGENVALUE, REFERENCE_AXIAL_EIGENVALUE),
        electronegative_wall_model=UNITY,
        electronegative_wall_geometry='fullFrequency',
        electronegative_wall_qualification=qualification(),
        wall_single_bath_approximation=True,
        thermo_source_assertions={
            'Ar+': 'ion', 'O+': 'ion', 'O2+': 'ion', 'O-': 'ion'})
    reactor.initialize_model(species, [], [], [])
    final_densities = dict(initial_densities)
    final_densities[electron] = frozen['electron_density_m_3']
    final_densities[om] = anion_density
    for item, density in final_densities.items():
        reactor.y[reactor.species_index[item]] = density / constants.Na
    return reactor, [arp, op, o2p]


def test_old_bracket_label_is_a_named_falsified_refusal():
    with pytest.raises(
            PlasmaStateError,
            match=r'electropositiveBracket.*FALSIFIED.*o2ReferenceQualifiedUnity'):
        model(closure='electropositiveBracket')


@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
def test_named_unity_closure_uses_the_existing_charge_balanced_operator(arm):
    reactor, _, _ = model(closure=UNITY, arm=arm, second_cation=True)
    y = reactor.y.copy()
    volume = reactor.compute_volume(y)
    frequencies = reactor.compute_ion_wall_frequencies(y, volume)
    components = reactor.compute_ion_wall_components(y, volume)
    cations = [j for j, charge in enumerate(reactor.species_charges) if charge > 0]
    anions = [j for j, charge in enumerate(reactor.species_charges) if charge < 0 and j != reactor.electron_index]

    assert reactor.electronegative_wall_factor(y) == (1., 1.)
    np.testing.assert_array_equal(frequencies[cations], components['total_ep'][cations])
    assert frequencies[cations[0]] != frequencies[cations[1]]
    assert all(reactor.wall_flux[j] == 0. for j in anions)
    cation_charge_flux = sum(
        reactor.species_charges[j] * frequencies[j] * y[j] for j in cations)
    assert -reactor.wall_flux[reactor.electron_index] == pytest.approx(cation_charge_flux)
    assert sum(charge * reactor.wall_flux[j]
               for j, charge in enumerate(reactor.species_charges)) == pytest.approx(0.)


def test_named_unity_qualification_and_manifest_record_owner_ruling():
    reactor, _, _ = model(closure=UNITY)
    result = qualify_envelope([(reactor, reactor.y.copy(), 7.)])
    assert result['samples'] == 1
    manifest = reactor.electronegative_wall_manifest()
    assert manifest['closure'] == UNITY
    assert manifest['claim_level'] == 'B3-QUALIFIED / MODEL-CONDITIONAL'
    assert manifest['claim_scope'] == (
        'closure qualification over the recorded I-314 O2 envelope only; '
        'the engine does not infer current-run envelope membership')
    assert manifest['qualifying_reference'] == (
        'plasma-pm3/reports/i314-en-wall-reference/rework3/report.md; '
        'I-314 rework 3, amended 2026-10-04')
    assert len(manifest['qualification_envelope']) == 18
    refused = [condition for condition in manifest['qualification_envelope']
               if condition['pressure_o2_torr'] == 0.05
               and condition['absorbed_power_w'] == 1.0]
    assert refused == [{
        'pressure_o2_torr': 0.05,
        'absorbed_power_w': 1.0,
        'reference_status': (
            'reference-unqualified (ion-heating validity), B3-pass numerically'),
    }]


def test_confined_anion_manifest_is_the_falsified_comparator():
    reactor, _, _ = model(closure='confinedAnion')
    manifest = reactor.electronegative_wall_manifest()
    assert manifest['scientific_status'] == 'FALSIFIED SIMPLIFIED CLOSURE'
    assert reactor.electronegative_wall_diagnostics['scientific_status'] == (
        'FALSIFIED SIMPLIFIED CLOSURE')
    assert manifest['closure'] == 'confinedAnion'


def test_named_unity_alpha_remains_dynamic_and_restart_preserves_selection():
    reactor, _, _ = model(closure=UNITY)
    first = reactor.y.copy()
    second = first.copy()
    second[4] = first[0]
    second[2] = first[0] + second[4]
    reactor.monitor_electronegative_wall(first, 1.)
    alpha_first = reactor.electronegative_wall_manifest()['alpha']
    reactor.monitor_electronegative_wall(second, 2.)
    manifest = reactor.electronegative_wall_manifest()
    assert alpha_first == pytest.approx(0.5)
    assert manifest['alpha'] == pytest.approx(1.)
    assert manifest['h'] == 1.
    rebuilt = pickle.loads(pickle.dumps(reactor))
    assert rebuilt.electronegative_wall_model == UNITY


def test_named_unity_residual_and_jacobian_share_the_wall_source():
    reactor, _, _ = model(closure=UNITY, second_cation=True)
    reactor.kf[:] = 0.
    y = reactor.y.copy()
    dydt = np.zeros_like(y)
    analytic = np.asarray(reactor.jacobian(0., y, dydt, 0.))
    numerical = np.zeros_like(analytic)
    columns = np.flatnonzero(y)
    for j in columns:
        delta = abs(y[j]) * 2.e-5
        plus, minus = y.copy(), y.copy()
        plus[j] += delta
        minus[j] -= delta
        numerical[:, j] = (
            np.asarray(reactor.residual(0., plus, dydt)[0])
            - np.asarray(reactor.residual(0., minus, dydt)[0])) / (2. * delta)
    np.testing.assert_allclose(
        analytic[:, columns], numerical[:, columns], rtol=2.e-6, atol=2.e-7)


def test_frozen_i314_oracle_replays_through_the_production_unity_operator():
    fixture = json.loads(ORACLE.read_text())
    assert fixture['selected_closure'] == UNITY
    assert fixture['legacy_source_arm'] == 'electropositiveBracket'
    charges = np.asarray(fixture['cation_charges'])
    tallies = {name: 0 for name in ('B1', 'B2', 'B3')}
    valid = 0
    logger.info('pressure_torr power_W reference_status max_error B1 B2 B3')
    for point in fixture['points']:
        state = np.asarray(point['closure_state'])
        h, factor, gradient = closure_factor_gradient(
            state, 0, [1], UNITY, 'fullFrequency', (1., 1.))
        assert h == factor == 1.
        np.testing.assert_array_equal(gradient, np.zeros_like(state))
        reactor, cation_species = _reference_geometry_o2_engine(point)
        assert reactor.electronegative_wall_model == UNITY
        volume = reactor.compute_volume(reactor.y)
        frequencies = reactor.compute_ion_wall_frequencies(reactor.y, volume)
        cation_indices = [reactor.species_index[item] for item in cation_species]
        cation_densities = np.asarray([
            reactor.y[index] * constants.Na / volume for index in cation_indices])
        np.testing.assert_allclose(
            cation_densities, point['engine_cation_number_densities_m_3'],
            rtol=1.e-7, atol=0.)
        np.testing.assert_allclose(
            frequencies[cation_indices],
            point['reference_geometry_engine_frequencies_s_1'],
            rtol=2.e-7, atol=0.)
        engine = frequencies[cation_indices] * cation_densities
        np.testing.assert_allclose(
            engine, point['engine_unity_cation_wall_vector'],
            rtol=3.e-7, atol=0.)
        errors = []
        for variant in point['reference_variants'].values():
            error = wall_discrepancy(
                engine, np.asarray(variant['reference_cation_wall_vector']), charges)
            assert error == pytest.approx(
                variant['error'], rel=2.e-6, abs=2.e-7)
            errors.append(error)
        robust = {name: max(errors) <= budget
                  for name, budget in point['budgets'].items()}
        if point['reference_status'] == 'reference-qualified':
            valid += 1
            for name, passed in robust.items():
                tallies[name] += int(passed)
        condition = point['condition']
        logger.info(
            '%s %s %s %s %s %s %s',
            condition['pressure_o2_torr'], condition['absorbed_power_w'],
            point['reference_status'], max(errors),
            *(robust[name] for name in ('B1', 'B2', 'B3')))
        assert all(row > 0 for row in point['source_rows'].values())
    assert valid == 17
    assert tallies == fixture['robust_valid_point_tallies'] == {
        'B1': 0, 'B2': 7, 'B3': 17}
    unqualified = [point for point in fixture['points']
                   if point['reference_status'] != 'reference-qualified']
    assert len(unqualified) == 1
    assert unqualified[0]['condition'] == {
        'absorbed_power_w': 1.0, 'pressure_o2_torr': 0.05}
    assert max(variant['error'] for variant in
               unqualified[0]['reference_variants'].values()) <= (
                   unqualified[0]['budgets']['B3'])


def test_falsified_comparator_oracle_fails_despite_passing_local_checks():
    fixture = json.loads(ORACLE.read_text())
    ratios = []
    for point in fixture['points']:
        comparator = point['falsified_comparator']
        assert comparator['local_checks']['attachment_metric'] > 1.
        assert comparator['local_checks']['minimum_confinement_ratio'] > 1.
        ratios.append(comparator['actual_error'] / comparator['B3_budget'])
    assert len(ratios) == 18
    assert min(ratios) == pytest.approx(4.147988621946723)
    assert max(ratios) == pytest.approx(18.134023001066044)
