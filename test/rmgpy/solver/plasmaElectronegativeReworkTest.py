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

"""Regressions for the five EN wall rework defects; references are test only."""
import copy
import csv
import json
from pathlib import Path

import numpy as np
import pytest
import rmgpy

from plasmaElectronegativeWallTest import model, qualification, fixture_reference
from rmgpy.exceptions import PlasmaStateError, ElectronegativeWallRegimeError
from rmgpy.solver.electronegative import qualify_envelope


@pytest.fixture(autouse=True)
def local_import():
    assert Path.cwd() in Path(rmgpy.__file__).resolve().parents
    import rmgpy.data.rmg as data_module
    saved = data_module.database
    data_module.database = None
    yield
    data_module.database = saved


@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
def test_subnormal_charged_state_refuses_before_infinite_jacobian(arm):
    r, _, _ = model(arm=arm)
    y = r.y.copy()
    y[0], y[4] = 1.e-310, 5.e-311
    y[2] = y[0] + y[4]
    old = copy.deepcopy(r.electronegative_wall_last_valid_state)
    with pytest.raises(ElectronegativeWallRegimeError, match='domain|normal|resolution'):
        r.check_wall_support(y)
    assert r.electronegative_wall_last_valid_state[0] == old[0]
    np.testing.assert_array_equal(r.electronegative_wall_last_valid_state[1], old[1])


@pytest.mark.parametrize('invalid', ['nan-electron', 'charge-imbalance', 'negative-neutral'])
@pytest.mark.parametrize('alpha', [0., 0.5])
def test_envelope_physical_validation_preserves_last_valid(invalid, alpha):
    r, _, _ = model(alpha=alpha)
    y = r.y.copy()
    if invalid == 'nan-electron':
        y[0] = np.nan
    elif invalid == 'charge-imbalance':
        y[2] *= 2.
    else:
        y[3] = -1.e-4
    old = copy.deepcopy(r.electronegative_wall_last_valid_state)
    manifest = r.electronegative_wall_manifest()
    with pytest.raises(PlasmaStateError):
        qualify_envelope([(r, y, 7.)])
    assert r.electronegative_wall_last_valid_state[0] == old[0]
    np.testing.assert_array_equal(r.electronegative_wall_last_valid_state[1], old[1])
    assert r.electronegative_wall_manifest() == manifest


def half_reference(context):
    result = fixture_reference(context)
    result['cation_flux'] *= 0.5
    return result


@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
@pytest.mark.parametrize('alpha', [0., 1.e-6, 0.5])
def test_public_energy_wall_jacobian_zero_small_and_gate_boundaries(arm, alpha):
    # Both C monitors sit exactly at their frozen discrepancy threshold.
    q = qualification(full_profile_reference=half_reference, geometry_reference=half_reference,
                      wall_threshold=2./3., geometry_threshold=2./3.)
    r, _, _ = model(alpha=alpha, arm=arm, energy=True, electronegative_wall_qualification=q)
    # Sample both sides of the A/B crossover at fixed state. Jacobians are
    # Newton operators, so the gate itself must not run here.
    transport = r.compute_anion_transport_frequencies(r.y, r.compute_volume(r.y))[4]
    ep = r.compute_ion_wall_components(r.y, r.compute_volume(r.y))['total_ep'][2]
    conc = r.y[3] / r.compute_volume(r.y)
    y = r.y.copy()
    for metric in (1.-1.e-6, 1., 1.+1.e-6):
        r.kf[:] = 0.
        wall = np.asarray(r.compute_electronegative_wall_jacobian(y))
        full = np.asarray(r.jacobian(0., y, np.zeros_like(y), 0.))
        numerical = np.zeros_like(wall)
        # Zero-anion face: electron/cation/neutral/Te derivatives stay on
        # the face. The anion derivative is discontinuous at the legacy switch
        # and is intentionally excluded there, as in the owner's alpha0 rule.
        columns = [j for j in range(len(y)) if alpha > 0. or j != 4]
        for j in columns:
            delta = abs(y[j]) * (1.e-3 if j in (1, 3) else 2.e-5)
            a, b = y.copy(), y.copy()
            a[j] += delta
            b[j] -= delta
            numerical[:, j] = (np.asarray(r.residual(0., a, np.zeros_like(y))[0]) -
                               np.asarray(r.residual(0., b, np.zeros_like(y))[0])) / (2.*delta)
        np.testing.assert_allclose(wall[:, columns], numerical[:, columns], rtol=3.e-5, atol=3.e-3)
        # The full energy Jacobian retains the legacy numerical remainder;
        # its small neutral columns lose several digits next to 1e13 charged columns.
        np.testing.assert_allclose(full[:, columns], numerical[:, columns], rtol=2.e-4, atol=3.e-3)
        # With the physical gate frequencies restored at the boundary, the
        # public full Jacobian must also match the actual chemistry+wall row.
        r.kf[0], r.kf[1] = metric * transport, metric * ep / conc
        full = np.asarray(r.jacobian(0., y, np.zeros_like(y), 0.))
        for j in columns:
            delta = abs(y[j]) * (1.e-4 if j in (1, 3) else 2.e-5)
            a, b = y.copy(), y.copy(); a[j] += delta; b[j] -= delta
            numerical[:, j] = (np.asarray(r.residual(0., a, np.zeros_like(y))[0]) -
                               np.asarray(r.residual(0., b, np.zeros_like(y))[0])) / (2.*delta)
        # The full energy Jacobian retains the legacy numerical remainder;
        # its small neutral columns lose several digits next to 1e13 charged columns.
        np.testing.assert_allclose(full[:, columns], numerical[:, columns], rtol=2.e-4, atol=3.e-3)


@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
@pytest.mark.parametrize('energy', [False, True])
def test_manifest_all_fields_follow_one_monitored_snapshot(arm, energy):
    r, _, _ = model(arm=arm, energy=energy)
    y = r.y.copy()
    y[4] = y[0]
    y[2] = y[0] + y[4]
    if energy:
        y[-1] *= 1.1
    r.monitor_electronegative_wall(y, 7.)
    record = r.electronegative_wall_manifest()
    assert record['time'] == 7.
    assert record['alpha'] == 1.
    assert record['h'] == 0.5
    np.testing.assert_array_equal(record['state'], y)
    assert record['electron_density'] == pytest.approx(y[0]*6.02214076e23/r.compute_volume(y))
    assert record['negative_ion_density'] == record['electron_density']
    assert record['floating_potential_e_over_kTe'] is not None
    flux = r.compute_ion_wall_frequencies(y, r.compute_volume(y))[2]*y[2]
    assert record['total_cation_wall_loss'] == pytest.approx(flux)
    assert record['wall_electron_energy_flux'] == pytest.approx(2.*8.314472*r.Te.value_si*flux, rel=1.e-5)
    if energy:
        assert record['energy_budget']['t'] == 7.
        assert record['energy_budget']['Te'] == y[-1]
    # Returning to alpha0 must not carry the previous sample's gate difference.
    y[4] = 0.; y[2] = y[0]
    r.monitor_electronegative_wall(y, 8.)
    record = r.electronegative_wall_manifest()
    assert record['time'] == 8.
    assert record['gates'] == 'not-evaluated-zero-anion'
    assert record['geometry_endpoint_wall_difference'] == 0.


@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
def test_real_deck_simulation_writes_wall_identity_gates_and_potential(tmp_path, arm):
    from rmgpy.rmg.main import RMG
    from rmgpy.rmg.input import read_input_file
    from rmgpy.rmg.listener import SimulationProfileWriter, SimulationProfilePlotter
    from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
    import rmgpy.data.rmg as data_module
    from rmgpy.kinetics import Arrhenius
    from rmgpy.reaction import Reaction
    from rmgpy.thermo import ThermoData
    q = qualification()
    for key in ('full_profile_reference', 'geometry_reference'):
        q[key] = fixture_reference.__module__ + ':fixture_reference'
    deck = tmp_path/'input.py'
    deck.write_text("""
database(thermoLibraries=[], reactionLibraries=[], seedMechanisms=[], kineticsFamilies='none')
species(label='e-', structure=adjacencyList('1 e u1 p0 c-1'))
species(label='Ar', structure=adjacencyList('1 Ar u0 p4 c0'))
species(label='Arp', structure=adjacencyList('multiplicity 2\\n1 Ar u1 p3 c+1'))
species(label='Cl', structure=adjacencyList('1 Cl u1 p3 c0'))
species(label='Cl-', structure=adjacencyList('1 Cl u0 p4 c-1'))
plasmaReactor(temperature=(298.15,'K'), pressure=(5,'torr'), electronTemperature=(35000.,'K'),
 initialMoleFractions={'Ar':.998996,'Cl':.001,'Arp':1.5e-6,'e-':1.e-6,'Cl-':.5e-6},
 chamberGeometry={'shape':'cylinder','radius':(5.,'cm'),'length':(30.,'cm')},
 ionReducedMobilities={'Arp':(1.535e-4,'m^2/(V*s)')},
 anionReducedMobilities={'Cl-':(1.5e-4,'m^2/(V*s)')},
 electronegativeWallModel='confinedAnion', electronegativeWallGeometry=%r,
 electronegativeWallQualification=%r, wallSingleBathApproximation=True,
 thermoSourceAssertions={'Arp': 'ion', 'Cl-': 'ion'}, terminationTime=(1.e-9,'s'))
simulator(atol=1.e-16, rtol=1.e-8)
model(toleranceMoveToCore=.1, toleranceInterruptSimulation=.1)
""" % (arm, q))
    job = RMG(); read_input_file(str(deck), job)
    r = job.reaction_systems[0]
    core = list(r.initial_mole_fractions)
    sp = {s.label:s for s in core}
    for s in core:
        if s.label != 'e-':
            s.thermo = ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'),
                Cpdata=([20.8]*7,'J/(mol*K)'), H298=(0.,'kJ/mol'), S298=(150.,'J/(mol*K)'))
    # Supplied synthetic mechanism: no database mutation or scientific validity claim.
    reactions = [Reaction(reactants=[sp['Cl-']], products=[sp['Cl'], sp['e-']],
        reversible=False, kinetics=Arrhenius(A=(1.e6,'s^-1'))),
        Reaction(reactants=[sp['e-'],sp['Cl']], products=[sp['Cl-']],
        reversible=False, kinetics=Arrhenius(A=(1.e9,'m^3/(mol*s)')))]
    data_module.database = None
    (tmp_path/'solver').mkdir()
    r.attach(SimulationProfileWriter(str(tmp_path), 0, core))
    r.attach(SimulationProfilePlotter(str(tmp_path), 0, core))
    result = r.simulate(core,reactions,[],[],[],[],
        model_settings=ModelSettings(tol_keep_in_edge=0.,tol_move_to_core=1.e5,tol_interrupt_simulation=1.e8),
        simulator_settings=SimulatorSettings())
    assert result[0]
    stem = tmp_path/'solver'/('simulation_1_%d' % len(core))
    manifest = json.loads(stem.with_suffix('.electronegative-wall.json').read_text())
    assert manifest['closure'] == 'confinedAnion'
    assert manifest['geometry_arm'] == arm
    assert set(manifest['gates']) == {'A','B','C_radial','C_geometry','C_full_profile'}
    assert np.isfinite(manifest['floating_potential_e_over_kTe'])
    with stem.with_suffix('.csv').open() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == len(r.snapshots)
    header = next(key for key in rows[-1] if key.startswith('Wall h'))
    assert 'confinedAnion' in header and arm in header
    assert float(rows[0][header]) != float(rows[-1][header])
    for row in rows:
        record = r.electronegative_wall_history[float(row['Time (s)'])]
        assert float(row[header]) == record['h']
        assert float(row['Wall gate A minimum conf']) == min(item['conf'] for item in record['gates']['A'].values())
    assert float(rows[-1]['Wall gate C full-profile error']) == 0.
    assert float(rows[-1]['Wall floating potential e/(kTe)']) == pytest.approx(manifest['floating_potential_e_over_kTe'])
    assert float(rows[-1]['Time (s)']) == manifest['time']
    assert stem.with_suffix('.png').stat().st_size > 1000
    # PNG metadata is read from the written production figure.
    from PIL import Image
    with Image.open(stem.with_suffix('.png')) as im:
        title = im.info['Title']
    assert 'confinedAnion' in title and arm in title and 'C_full_profile' in title
    assert 'phi=' in title


@pytest.mark.parametrize('entry', [(float('nan'),'m^2/(V*s)'), (0.,'m^2/(V*s)'), (1.,'m')])
def test_constructor_validates_unused_anion_mobility(entry):
    with pytest.raises(PlasmaStateError, match='unused'):
        model(anion_reduced_mobilities={'Cl-':(1.5e-4,'m^2/(V*s)'), 'unused':entry})


def test_zero_current_profile_records_unavailable_potential_and_plots(tmp_path):
    from rmgpy.rmg.listener import SimulationProfileWriter, SimulationProfilePlotter
    from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
    from rmgpy.solver.termination import TerminationTime
    r, core, _ = model(alpha=0., qualify=False, reactions=False,
                       ionisation_source=(1.e14, 'm^-3*s^-1'))
    r.termination = [TerminationTime((1.e-12, 's'))]
    for sp in core:
        if sp.get_net_charge() != 0:
            r.initial_mole_fractions[sp] = 0.
    (tmp_path/'solver').mkdir()
    r.attach(SimulationProfileWriter(str(tmp_path), 0, core))
    r.attach(SimulationProfilePlotter(str(tmp_path), 0, core))
    r.simulate(core,[],[],[],[],[],
        model_settings=ModelSettings(tol_keep_in_edge=0.,tol_move_to_core=1.e5,tol_interrupt_simulation=1.e8),
        simulator_settings=SimulatorSettings())
    stem = tmp_path/'solver'/('simulation_1_%d' % len(core))
    record = json.loads(stem.with_suffix('.electronegative-wall.json').read_text())
    assert r.electronegative_wall_history[0.]['floating_potential_e_over_kTe'] is None
    assert np.isfinite(record['floating_potential_e_over_kTe'])
    with stem.with_suffix('.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    assert np.isnan(float(rows[0]['Wall floating potential e/(kTe)']))
    assert stem.with_suffix('.png').stat().st_size > 1000
