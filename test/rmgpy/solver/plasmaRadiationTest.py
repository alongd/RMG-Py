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

"""Trapped radiation through the compiled PlasmaReactor boundary."""
import copy
import pickle
from types import SimpleNamespace
import numpy as np
import pytest
import rmgpy.constants as constants
import rmgpy.data.rmg as data_module
from rmgpy.data.base import Entry
from rmgpy.data.thermo import ThermoDatabase, ThermoLibrary
from rmgpy.exceptions import PlasmaStateError
from rmgpy.species import Species
from rmgpy.solver.plasma import PlasmaReactor
from plasmaEEDFCoreTest import thermo


@pytest.fixture
def radiation_model(monkeypatch):
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ground = Species(label='ground-renamed').from_adjacency_list('1 Ar u0 p4 c0')
    ion = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    r4 = Species(label='resonance-a').from_adjacency_list('multiplicity 3\nelectronicstate 3P1\n1 Ar u2 p3 c0')
    r2 = Species(label='resonance-b').from_adjacency_list('multiplicity 1\nelectronicstate 1P1\n1 Ar u0 p4 c0')
    core = [electron, ground, ion, r4, r2]
    library = ThermoLibrary(label='RadiationTest')
    library.thermo_convention = 'ion'
    for i, (spc, energy) in enumerate(zip(core[1:], [0., 15.76, 11.624, 11.828])):
        spc.thermo = thermo(energy)
        spc.thermo.Cp0 = (20., 'J/(mol*K)')
        spc.thermo.CpInf = (20., 'J/(mol*K)')
        library.entries[i] = Entry(index=i, label=spc.label, item=spc.molecule[0].copy(deep=True), data=copy.deepcopy(spc.thermo))
    database = ThermoDatabase()
    database.libraries = {library.label: library}
    database.library_order = [library.label]
    monkeypatch.setattr(data_module, 'database', SimpleNamespace(thermo=database, solvation=None))
    def build(model='holsteinWalsh', geometry=True, bath=None, energy=False, feed=False):
        opts = {} if model is None else {'resonance_radiation': {'model': model, 'geometryModel': 'infiniteCylinder'}}
        if geometry:
            opts['radiation_chamber_geometry'] = {'shape': 'cylinder', 'radius': .05, 'length': .30}
        species = list(core)
        fractions = dict(zip(core,[1e-6,.999978,1e-6,1e-5,1e-5]))
        if bath is not None:
            adjacency, amount = bath
            gas = Species(label='bath').from_adjacency_list(adjacency)
            gas.thermo = thermo()
            species.append(gas)
            fractions[gas] = amount
            fractions[ground] -= amount
        if energy:
            opts.update(diffusion_length=(.02,'m'), ion_reduced_mobility=(1e-4,'m^2/(V*s)'),
                wall_single_bath_approximation=True, wall_neutralization_products={'Ar+':'ground-renamed'},
                electron_energy_balance={'absorbed_power':(.5,'W'), 'chamber_volume':(.002,'m^3'),
                    'sheath':'floating_wall','electron_energies':{},
                    'elastic_collisions':{s.label:{'ignore':'isolate photon energy accounting'}
                                          for s in species if not s.get_net_charge()}})
        reactor = PlasmaReactor((298.15,'K'),(5,'torr'), fractions, (11604.5,'K'), termination=[], **opts)
        reactions = []
        if feed:
            from rmgpy.reaction import Reaction
            from rmgpy.kinetics import TwoTemperaturePlasma
            reactions = [Reaction(reactants=[electron,ground],products=[electron,upper],reversible=False,
                kinetics=TwoTemperaturePlasma(A=(k,'m^3/(molecule*s)'),n=0.,
                    Ea_g=(0.,'J/mol'),Ea_e=(0.,'J/mol'))) for upper,k in ((r4,1e-17),(r2,2e-17))]
        reactor.initialize_model(species,reactions,[],[])
        return reactor
    return core, build


def test_decay_returns_ground_and_publishes_only_accepted_states(radiation_model):
    core, build = radiation_model
    reactor = build()
    initial = copy.deepcopy(reactor.radiation_diagnostics)
    y = reactor.y0.copy()
    residual, status = reactor.residual(0., y, np.zeros_like(y))
    assert status == 1
    assert residual[3] < 0 and residual[4] < 0
    assert residual[1] == pytest.approx(-residual[3]-residual[4])
    assert residual[0] == 0 and residual[2] == 0
    assert reactor.radiation_diagnostics == initial
    trial = y.copy(); trial[1] *= .9
    reactor.residual(99., trial, np.zeros_like(y))
    assert reactor.radiation_diagnostics == initial
    reactor.advance(1e-6)
    assert reactor.radiation_diagnostics['time_s'] == reactor.t
    assert reactor.radiation_diagnostics['photon_power_W'] > 0


def test_complete_is_exactly_zero(radiation_model):
    core, build = radiation_model
    reactor = build('complete')
    residual, _ = reactor.residual(0., reactor.y0.copy(), np.zeros_like(reactor.y0))
    assert np.array_equal(residual, np.zeros_like(residual))
    assert all(row['g'] == 0. for row in reactor.radiation_diagnostics['lines'])
    reactor.advance(1e-6)
    assert reactor.radiation_diagnostics['photon_power_W'] == 0.


def test_resolved_resonance_requires_explicit_model(radiation_model):
    with pytest.raises(PlasmaStateError, match='resonanceRadiation'):
        radiation_model[1](None)


def test_complete_still_requires_geometry(radiation_model):
    with pytest.raises(PlasmaStateError, match='geometry'):
        radiation_model[1]('complete', geometry=False)


def test_jacobian_includes_opacity_and_volume_derivatives(radiation_model):
    reactor = radiation_model[1]()
    y = reactor.y0.copy()
    analytic = reactor.jacobian(0., y, np.zeros_like(y), 0.)
    for k in range(len(y)):
        h = max(abs(y[k])*1e-5,1e-14)
        yp, ym = y.copy(), y.copy(); yp[k]+=h; ym[k]-=h
        rp = reactor.residual(0.,yp,np.zeros_like(y))[0]
        rm = reactor.residual(0.,ym,np.zeros_like(y))[0]
        assert analytic[:,k] == pytest.approx((rp-rm)/(2*h), rel=2e-4, abs=2e-5)


def test_pickle_preserves_declaration(radiation_model):
    reactor = radiation_model[1]()
    clone = pickle.loads(pickle.dumps(reactor))
    assert clone.resonance_radiation == reactor.resonance_radiation
    assert clone.radiation_chamber_geometry == reactor.radiation_chamber_geometry


@pytest.mark.parametrize('bath', [('1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}', .20),
                                  ('multiplicity 3\n1 Ar u2 p3 c0', 1e-4)])
def test_admixture_uses_ground_argon_partial_density(radiation_model, bath):
    core, build = radiation_model
    reactor = build(bath=bath)
    row = reactor.radiation_diagnostics
    assert row['ground_density_m_3'] == pytest.approx(reactor.y0[1]*constants.Na/reactor.compute_volume(reactor.y0),rel=1e-6)
    assert all(0 < line['g'] < 1 for line in row['lines'])
    if bath[0].startswith('multiplicity'):
        residual,_ = reactor.residual(0.,reactor.y0.copy(),np.zeros_like(reactor.y0))
        assert residual[-1] == 0.  # No implicit metastable loss.


@pytest.mark.parametrize('bath', [('1 He u0 p1 c0',.01),
                                  ('1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}',.30)])
def test_unsupported_bath_refuses(radiation_model, bath):
    with pytest.raises(PlasmaStateError, match='radiation domain'):
        radiation_model[1](bath=bath)


def test_photons_do_not_cool_electrons_twice(radiation_model):
    build = radiation_model[1]
    trapped, complete = build(energy=True), build('complete',energy=True)
    y = trapped.y0.copy()
    rt,_ = trapped.residual(0.,y,np.zeros_like(y))
    rc,_ = complete.residual(0.,y,np.zeros_like(y))
    assert trapped.radiation_diagnostics['photon_power_W'] > 0
    assert rt[-1] == pytest.approx(rc[-1],rel=1e-12)
    assert rt[3] < rc[3] and rt[4] < rc[4]
    analytic = trapped.jacobian(0.,y,np.zeros_like(y),0.)
    for k in (1,3,4,len(y)-1):
        h = abs(y[k])*2e-5
        yp,ym = y.copy(),y.copy(); yp[k]+=h; ym[k]-=h
        fd = (trapped.residual(0.,yp,np.zeros_like(y))[0]-trapped.residual(0.,ym,np.zeros_like(y))[0])/(2*h)
        # The finite-difference energy row is cancellation-dominated for tiny
        # neutral perturbations; photon/electron separation is checked above.
        assert analytic[:-1,k] == pytest.approx(fd[:-1],rel=1e-3,abs=2e-5)


def test_fed_resonances_reach_steady_state_in_compiled_reactor(radiation_model):
    # A controlled kinetic integration benchmark; these constant excitation laws
    # are deliberately synthetic, not production EEDF evidence.
    reactor = radiation_model[1](feed=True)
    previous = None
    for time in (1e-3,2e-3,3e-3):
        reactor.advance(time)
        rows = reactor.radiation_diagnostics['lines']
        assert len(rows) == 2
        populations = np.array([row['density_m_3'] for row in rows])
        assert np.all(populations > 0)
        residual,_ = reactor.residual(time,reactor.y.copy(),np.zeros_like(reactor.y))
        for index,row in zip((3,4),rows):
            assert abs(residual[index])/row['loss_mol_s'] < 1e-7
        if previous is not None:
            assert populations == pytest.approx(previous,rel=1e-7)
        previous = populations
    complete = radiation_model[1]('complete',feed=True)
    complete.advance(3e-5)
    assert complete.radiation_diagnostics['photon_power_W'] == 0.
    assert all(row['g'] == row['loss_mol_s'] == 0. for row in complete.radiation_diagnostics['lines'])
    assert complete.y[3] > complete.y0[3] and complete.y[4] > complete.y0[4]


def test_declaration_views_cannot_diverge_from_live_operator(radiation_model):
    reactor = radiation_model[1]()
    declaration = reactor.resonance_radiation
    geometry = reactor.radiation_chamber_geometry
    declaration['model'] = 'complete'
    geometry['radius'] = 1.
    assert reactor.resonance_radiation['model'] == 'holsteinWalsh'
    assert reactor.radiation_chamber_geometry['radius'] == .05
    assert reactor.radiation_manifest()['declaration']['model'] == 'holsteinWalsh'


def test_compiled_accepted_trajectory_writes_native_radiation_json(radiation_model, tmp_path):
    import json
    from rmgpy.rmg.listener import SimulationProfileWriter
    core, build = radiation_model
    reactor = build(feed=True)
    reactor.snapshots = [[0.,reactor.compute_volume(reactor.y0)] + list(reactor.y0)]
    reactor.advance(1e-3)
    reactor.snapshots.append([reactor.t,reactor.compute_volume(reactor.y)] + list(reactor.y))
    (tmp_path/'solver').mkdir()
    SimulationProfileWriter(str(tmp_path),0,core).update(reactor)
    record = json.loads((tmp_path/'solver/simulation_1_5.radiation.json').read_text())
    assert [row['time_s'] for row in record['trajectory']] == [0.,1e-3]
    assert record['trajectory'][-1] == record['accepted']
    assert {row['state'] for row in record['accepted']['lines']} == {'3P1','1P1'}
    assert record['accepted']['photon_power_W'] > 0


def test_reinitialization_rebinds_resonance_indices(radiation_model):
    core, build = radiation_model
    reactor = build()
    reference = {row['state']: row['loss_mol_s'] for row in reactor.radiation_diagnostics['lines']}
    reordered = list(reversed(core))
    reactor.initialize_model(reordered,[],[],[])
    residual,_ = reactor.residual(0.,reactor.y0.copy(),np.zeros_like(reactor.y0))
    assert residual[reordered.index(core[3])] == pytest.approx(-reference['3P1'])
    assert residual[reordered.index(core[4])] == pytest.approx(-reference['1P1'])
    assert residual[reordered.index(core[0])] == 0.
    assert set(reactor.radiation_history) == {0.}


def test_nonfinite_resonance_trial_refuses_without_publishing(radiation_model):
    reactor = radiation_model[1]()
    accepted = copy.deepcopy(reactor.radiation_diagnostics)
    trial = reactor.y0.copy(); trial[3] = float('nan')
    with pytest.raises(PlasmaStateError,match='non-physical volume'):
        reactor.residual(0.,trial,np.zeros_like(trial))
    assert reactor.radiation_diagnostics == accepted
