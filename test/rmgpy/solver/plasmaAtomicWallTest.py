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

"""Atomic wall recombination: flux, energy and derivatives, independently checked."""
import copy
import pickle

import numpy as np
import pytest

from rmgpy import constants
from rmgpy.exceptions import PlasmaStateError
from rmgpy.species import Species
from rmgpy.thermo import ThermoData
from plasmaWallTest import _fixture_plasma_reactor, _argon_thermo


DN = 7.05616e20
RADIUS, LENGTH = .05, .30
COMPONENTS = ((2.405 / RADIUS)**2, (np.pi / LENGTH)**2)


def declaration(gamma=5e-5):
    return {'kind': 'homonuclearRecombination', 'product': 'N2',
            'diffusivity': (DN, '1/(m*s)'), 'recombinationProbability': gamma}


def reactor(gamma=5e-5, entry=None, geometry=True, atom=None, product=None,
            energy=False, blanc=False, nitrogen_enthalpy=472.68):
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    n = Species(label='N').from_adjacency_list(atom or 'multiplicity 4\n1 N u3 p1 c0')
    n2 = Species(label='N2').from_adjacency_list(
        product or '1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}')
    ar.thermo, arp.thermo = _argon_thermo(0), _argon_thermo(15.76)
    for sp, h, cp in [(n, nitrogen_enthalpy, 2.5), (n2, 0., 3.5)]:
        sp.thermo = ThermoData(Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
                              Cpdata=([cp * constants.R]*7, 'J/(mol*K)'),
                              H298=(h, 'kJ/mol'), S298=(150, 'J/(mol*K)'))
    core = [electron, ar, arp, n, n2]
    fractions = dict(zip(core, [1e-6, .899998, 1e-6, .01, .09]))
    entry = declaration(gamma) if entry is None else entry
    kwargs = dict(diffusion_length=(sum(COMPONENTS)**-.5, 'm'),
                  ion_reduced_mobility=(1.535e-4, 'm^2/(V*s)'),
                  wall_single_bath_approximation=True,
                  wall_neutralization_products={'Ar+': 'Ar'},
                  wall_neutral_diffusion={'N': entry})
    if geometry:
        kwargs.update(wall_chamber_geometry={'shape': 'cylinder', 'radius': RADIUS, 'length': LENGTH},
                      wall_diffusion_components=COMPONENTS)
    if blanc:
        del kwargs['ion_reduced_mobility']
        del kwargs['wall_single_bath_approximation']
        entry = dict(entry, diffusivity={'Ar': (DN, '1/(m*s)'),
                                        'N': (DN*.8, '1/(m*s)'),
                                        'N2': (DN*1.2, '1/(m*s)')})
        kwargs['wall_neutral_diffusion'] = {'N': entry}
        kwargs['ion_reduced_mobilities'] = {'Ar+': {'perBath': {
            name: (1.535e-4, 'm^2/(V*s)') for name in ['Ar', 'N', 'N2']}}}
    if energy:
        kwargs['electron_energy_balance'] = {
            'absorbed_power': (.5, 'W'),
            'chamber_volume': (np.pi*RADIUS**2*LENGTH, 'm^3'),
            'sheath': 'floating_wall', 'electron_energies': {},
            'elastic_collisions': {name: {'ignore': 'Synthetic wall-only derivative test'}
                                   for name in ['Ar', 'N', 'N2']}}
    r = _fixture_plasma_reactor((300, 'K'), (666.61184, 'Pa'), fractions, (30000, 'K'), **kwargs)
    r.initialize_model(core, [], [], [])
    y = np.array([fractions[s] for s in core] + ([30000.] if energy else []))
    return r, core, y


@pytest.mark.parametrize('gamma', [0., 1e-5, 5e-5, 1e-3, 1.])
def test_recombination_flux_atoms_and_enthalpy(gamma):
    r, core, y = reactor(gamma)
    volume = r.compute_volume(y)
    density = sum(y[[1, 3, 4]]) * constants.Na / volume
    speed = np.sqrt(8*constants.R*r.T.value_si/(np.pi*core[3].molecular_weight.value_si))
    vs = RADIUS*LENGTH/(2*(RADIUS+LENGTH))
    expected = (1/(density/(DN*sum(COMPONENTS)) + 2*vs*(2-gamma)/(gamma*speed))
                if gamma else 0.)
    assert r.compute_neutral_wall_frequencies(y, volume)[3] == pytest.approx(expected, rel=1e-12)
    residual = r.residual(0., y, np.zeros_like(y))[0]
    loss = expected*y[3]
    assert residual[3] == pytest.approx(-loss)
    assert residual[4] == pytest.approx(.5*loss)
    assert residual[3]+2*residual[4] == pytest.approx(0., abs=1e-15)
    assert residual[3]+residual[4] == pytest.approx(-.5*loss)
    r._latch_wall_diagnostics(y, volume, 0.)
    np.testing.assert_allclose(r.wall_flux[[3, 4]], [-loss, .5*loss])
    dh = core[3].get_enthalpy(300)-.5*core[4].get_enthalpy(300)
    assert r.wall_recombination_energy_flux == pytest.approx(loss*dh)
    assert 'available' in r.wall_energy_availability['wall_recombination_energy_flux']


@pytest.mark.parametrize('energy,blanc', [(False, False), (True, False), (False, True), (True, True)])
@pytest.mark.parametrize('gamma', [0., 5e-5, 1.])
def test_recombination_jacobian(energy, blanc, gamma):
    r, _, y = reactor(gamma, energy=energy, blanc=blanc)
    zero = np.zeros_like(y)
    analytic = np.array(r._jacobian(0., y, zero, 0.))
    fd = np.zeros_like(analytic)
    for col in range(len(y)):
        h = abs(y[col])*1e-4
        yp, ym = y.copy(), y.copy()
        yp[col] += h
        ym[col] -= h
        fd[:, col] = (r.residual(0., yp, zero)[0]-r.residual(0., ym, zero)[0])/(2*h)
    # Per-row scaling keeps the large electron-energy derivative from hiding
    # a wrong half-yield in the much smaller nitrogen rows.
    scale = np.maximum(np.max(np.abs(fd), axis=1), 1e-12)
    np.testing.assert_allclose(analytic/scale[:, None], fd/scale[:, None], rtol=2e-4, atol=1e-7)


def test_recombination_roundtrip():
    r, _, _ = reactor()
    for clone in [copy.deepcopy(r), pickle.loads(pickle.dumps(r))]:
        assert clone.wall_neutral_diffusion['N'] == declaration()
        assert clone.wall_chamber_geometry == r.wall_chamber_geometry


def test_ordinary_wall_still_refuses_atomic_molecular_product():
    with pytest.raises(PlasmaStateError, match='not the same nuclei and skeleton'):
        reactor(entry={'product': 'N2', 'diffusivity': (DN, '1/(m*s)')})


def test_recombination_requires_geometry():
    with pytest.raises(PlasmaStateError, match='chamber|cylinder'):
        reactor(geometry=False)


@pytest.mark.parametrize('change', [
    {'kind': 'arbitrary'}, {'recombinationProbability': None},
    {'recombinationProbability': True}, {'recombinationProbability': -1.},
    {'recombinationProbability': 1.1}, {'recombinationProbability': float('nan')},
    {'recombinationProbability': float('inf')}, {'yield': .5},
])
def test_invalid_recombination_declarations(change):
    with pytest.raises(PlasmaStateError):
        reactor(entry=dict(declaration(), **change))


def test_recombination_requires_declared_probability():
    entry = declaration()
    del entry['recombinationProbability']
    with pytest.raises(PlasmaStateError, match='recombinationProbability'):
        reactor(entry=entry)


def test_probability_does_not_change_ordinary_wall_kind():
    entry = declaration()
    del entry['kind']
    with pytest.raises(PlasmaStateError):
        reactor(entry=entry)


@pytest.mark.parametrize('atom,product', [
    ('1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}', None),
    (None, '1 N u3 p1 c0'),
    (None, '1 N u3 p1 c0\n2 N u3 p1 c0'),
    (None, '1 N u0 p1 c0 {2,T}\n2 C u1 p0 c0 {1,T}'),
    ('1 C u4 p0 c0', 'multiplicity 3\n1 C u1 p0 c0 i13 {2,T}\n2 C u1 p0 c0 i13 {1,T}'),
    (None, 'multiplicity 2\n1 N u0 p1 c0 {2,T}\n2 N u1 p0 c+1 {1,T}'),
])
def test_only_neutral_homonuclear_atom_to_diatomic_is_allowed(atom, product):
    with pytest.raises(PlasmaStateError):
        reactor(atom=atom, product=product)


def test_recombination_cannot_raise_material_enthalpy():
    with pytest.raises(PlasmaStateError, match='energy|enthalpy|H298'):
        reactor(nitrogen_enthalpy=-10.)
