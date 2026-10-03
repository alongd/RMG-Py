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

"""Second EN-wall rework: numerical scale and real infinite-confinement output."""
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
from rmgpy.solver.plasma import PlasmaReactor


@pytest.fixture(autouse=True)
def local_import():
    assert Path.cwd() in Path(rmgpy.__file__).resolve().parents
    import rmgpy.data.rmg as data_module
    saved = data_module.database
    data_module.database = None
    yield
    data_module.database = saved


@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
def test_energy_jacobians_refuse_outside_concentration_domain(arm):
    r, _, _ = model(arm=arm, energy=True)
    y = r.y.copy()
    y[0], y[4], y[2] = 1.e-200, 5.e-201, 1.5e-200
    old=copy.deepcopy(r.electronegative_wall_last_valid_state)
    with pytest.raises(ElectronegativeWallRegimeError, match='domain.*e-'):
        r.monitor_electronegative_wall(y, 0.)
    assert r.electronegative_wall_last_valid_state[0]==old[0]
    np.testing.assert_array_equal(r.electronegative_wall_last_valid_state[1],old[1])


@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
@pytest.mark.parametrize('energy', [False, True])
@pytest.mark.parametrize('scaling', ['charge-only', 'whole-inventory'])
@pytest.mark.parametrize('scale', [1.e-300, 1.e-250, 1.e-200, 1.e-160, 1.e-100,
                                  1.e-30, 1.e-6, 1., 1.e15, 1.e30])
def test_accepted_scale_sweep_has_finite_matching_jacobians(arm, energy, scaling, scale):
    r, _, _ = model(arm=arm, energy=energy)
    y = r.y.copy()
    if scaling == 'whole-inventory':
        y[:r.num_core_species] *= scale/y[0]
    y[0], y[4], y[2] = scale, 0.5*scale, 1.5*scale
    # Keep the chemical derivative representable at accepted inventory scales;
    # states outside the physical concentration/volume domain are refused.
    r.kf[0] = 1.e3
    try:
        r.monitor_electronegative_wall(y, 0.)
    except ElectronegativeWallRegimeError as error:
        assert 'domain:' in str(error)
        V=r.compute_volume(y)
        concentrations=y[:r.num_core_species]/V
        assert (V < 1.e-12 or V > 1.e3 or
                any(0. < c < 1.e-30 or c > (1.e3 if r.species_charges[j] else 1.e6)
                    for j,c in enumerate(concentrations)))
        return
    except PlasmaStateError as error:
        # Large charge-only states exceed the unchanged ionisation-degree gate.
        assert scaling == 'charge-only' and scale >= 1.
        assert 'ionisation' in str(error) or 'neutral heavy number density' in str(error)
        return
    residual = np.asarray(r.residual(0., y, np.zeros_like(y))[0])
    assert np.isfinite(residual).all()
    full = np.asarray(r.jacobian(0., y, np.zeros_like(y), 0.))
    wall = np.asarray(r.compute_electronegative_wall_jacobian(y))
    assert np.isfinite(full).all()
    assert np.isfinite(wall).all()
    # Differentiate the actual residual, preserving tiny component scales.
    fd = np.zeros_like(full)
    for j in range(len(y)):
        delta = y[j]*(1.e-3 if j in (1, 3) else 2.e-5)
        a, b = y.copy(), y.copy()
        a[j] += delta
        b[j] -= delta
        fd[:, j] = (np.asarray(r.residual(0., a, np.zeros_like(y))[0]) -
                    np.asarray(r.residual(0., b, np.zeros_like(y))[0]))/(2.*delta)
    assert np.isfinite(fd).all()
    # Normalize each derivative by its row and component scale so tolerances
    # cannot hide small derivatives at large inventories or large derivatives
    # at tiny inventories. Neutral finite differences incur cancellation.
    row_scale = np.maximum(np.max(np.abs(fd*y[None, :]), axis=1), np.abs(residual))
    row_scale[row_scale == 0.] = 1.
    np.testing.assert_allclose((full/row_scale[:, None])*y[None, :],
                               (fd/row_scale[:, None])*y[None, :], rtol=3.e-4, atol=3.e-6)
    # Public helper is the wall operator; compare it after removing chemistry.
    r.kf[:] = 0.
    wall_fd = np.zeros_like(wall)
    for j in range(len(y)):
        delta = y[j]*(1.e-3 if j in (1, 3) else 2.e-5)
        a, b = y.copy(), y.copy()
        a[j] += delta
        b[j] -= delta
        wall_fd[:, j] = (np.asarray(r.residual(0., a, np.zeros_like(y))[0]) -
                         np.asarray(r.residual(0., b, np.zeros_like(y))[0]))/(2.*delta)
    row_scale = np.max(np.abs(wall_fd*y[None, :]), axis=1)
    row_scale[row_scale == 0.] = 1.
    np.testing.assert_allclose((wall/row_scale[:, None])*y[None, :],
                               (wall_fd/row_scale[:, None])*y[None, :], rtol=3.e-4, atol=3.e-6)


@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
def test_unrepresentable_energy_derivative_refuses_without_last_valid_overwrite(arm):
    r, _, _ = model(arm=arm, energy=True)
    y = r.y.copy()
    y[0], y[4], y[2] = 1.e-300, 5.e-301, 1.5e-300
    old = copy.deepcopy(r.electronegative_wall_last_valid_state)
    with pytest.raises(ElectronegativeWallRegimeError, match='domain.*e-|energy Jacobian.*range'):
        r.check_wall_support(y)
    assert r.electronegative_wall_last_valid_state[0] == old[0]
    np.testing.assert_array_equal(r.electronegative_wall_last_valid_state[1], old[1])


@pytest.mark.parametrize('time', [float('nan'), float('inf'), -float('inf'), -1.])
def test_envelope_invalid_time_preserves_last_valid(time):
    r, _, _ = model()
    old = copy.deepcopy(r.electronegative_wall_last_valid_state)
    manifest = r.electronegative_wall_manifest()
    with pytest.raises(ElectronegativeWallRegimeError, match='time'):
        qualify_envelope([(r, r.y.copy(), time)])
    assert r.electronegative_wall_last_valid_state[0] == old[0]
    np.testing.assert_array_equal(r.electronegative_wall_last_valid_state[1], old[1])
    assert r.electronegative_wall_manifest() == manifest


def test_envelope_bracket_named_refusal_preserves_last_valid():
    r, _, _ = model(closure='electropositiveBracket')
    manifest = r.electronegative_wall_manifest()
    with pytest.raises(ElectronegativeWallRegimeError, match='electropositiveBracket'):
        qualify_envelope([(r, r.y.copy(), 7.)])
    assert r.electronegative_wall_manifest() == manifest


class CountDenseFixedTeReactor(PlasmaReactor):
    def _en_wall_linearization(self, y, volume):
        self.dense_calls = getattr(self, 'dense_calls', 0) + 1
        return super()._en_wall_linearization(y, volume)


def test_fixed_te_monitor_checks_full_jacobian_and_latches_scalar_potential():
    r, _, _ = model(reactor_cls=CountDenseFixedTeReactor)
    assert getattr(r, 'dense_calls', 0) == 1
    r.monitor_electronegative_wall(r.y, 7.)
    # Rework 5 requires one complete operator check per acceptance. The
    # scalar potential diagnostics must not construct another dense matrix.
    assert r.dense_calls == 2
    normal, _, _ = model()
    normal.monitor_electronegative_wall(normal.y, 7.)
    assert r.electronegative_wall_manifest()['floating_potential_e_over_kTe'] == pytest.approx(
        normal._en_wall_linearization(normal.y, normal.compute_volume(normal.y))[2], rel=1.e-14)


@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
@pytest.mark.parametrize('mixed', [False, True])
def test_real_deck_low_transport_positive_destruction_writes_plot(tmp_path, arm, mixed):
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
    # Low in-domain mobility keeps destruction dominant. This synthetic
    # input exercises the production deck/transport/output path; impossible
    # subnormal mobilities are refused by the rework-7 configuration tests.
    extra_species = "species(label='O', structure=adjacencyList('1 O u2 p2 c0'))\nspecies(label='O-', structure=adjacencyList('1 O u1 p3 c-1'))" if mixed else ''
    fractions = {'Ar':.998996,'Cl':.001,'Arp':1.5e-6,'e-':1.e-6,'Cl-':.5e-6}
    mobilities = {'Cl-':(1.5e-6,'m^2/(V*s)')}
    assertions = {'Arp': 'ion', 'Cl-': 'ion'}
    if mixed:
        fractions.update({'O':0., 'O-':1.e-7, 'Arp':1.6e-6})
        mobilities['O-'] = (1.5e-4,'m^2/(V*s)')
        assertions['O-'] = 'ion'
    deck.write_text("""
database(thermoLibraries=[], reactionLibraries=[], seedMechanisms=[], kineticsFamilies='none')
species(label='e-', structure=adjacencyList('1 e u1 p0 c-1'))
species(label='Ar', structure=adjacencyList('1 Ar u0 p4 c0'))
species(label='Arp', structure=adjacencyList('multiplicity 2\\n1 Ar u1 p3 c+1'))
species(label='Cl', structure=adjacencyList('1 Cl u1 p3 c0'))
species(label='Cl-', structure=adjacencyList('1 Cl u0 p4 c-1'))
%s
plasmaReactor(temperature=(298.15,'K'), pressure=(5,'torr'), electronTemperature=(35000.,'K'),
 initialMoleFractions=%r,
 chamberGeometry={'shape':'cylinder','radius':(5.,'cm'),'length':(30.,'cm')},
 ionReducedMobilities={'Arp':(1.535e-4,'m^2/(V*s)')},
 anionReducedMobilities=%r,
 electronegativeWallModel='confinedAnion', electronegativeWallGeometry=%r,
 electronegativeWallQualification=%r, wallSingleBathApproximation=True,
 thermoSourceAssertions=%r, terminationTime=(1.e-9,'s'))
simulator(atol=1.e-16, rtol=1.e-8)
model(toleranceMoveToCore=.1, toleranceInterruptSimulation=.1)
""" % (extra_species, fractions, mobilities, arm, q, assertions))
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
    if mixed:
        reactions.append(Reaction(reactants=[sp['O-']], products=[sp['O'], sp['e-']],
            reversible=False, kinetics=Arrhenius(A=(1.e6,'s^-1'))))
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
        assert float(row['Wall gate A minimum conf']) == min(float('inf') if item['conf'] == 'infinite' else item['conf']
            for item in record['gates']['A'].values())
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

    assert manifest['gates']['A']['Cl-']['conf'] > 1.e5
    assert manifest['gates']['A']['Cl-']['transport'] > 0.
    if mixed:
        assert np.isfinite(manifest['gates']['A']['O-']['conf'])
    for row in rows:
        conf = float(row['Wall gate A minimum conf'])
        assert np.isfinite(conf) and conf > 1.
    assert 'min conf=' in title
