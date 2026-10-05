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

"""Full EN energy derivatives, public state validation and accepted operators."""

from decimal import Decimal, localcontext
from itertools import product
from pathlib import Path
import copy
import numpy as np
import pytest
import rmgpy
from rmgpy import constants
from rmgpy.kinetics import Arrhenius
from rmgpy.reaction import Reaction
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.species import Species
from rmgpy.exceptions import ElectronegativeWallRegimeError, PlasmaStateError
from plasmaElectronegativeWallTest import model, FixtureEnergyReactor
from plasmaElectronegativeRework4Test import set_state, blanc_overrides


def decimal_operator(r, state, full=False, values=False):
    """Differentiate the complete constitutive law with 800-digit dual numbers.

    EOS, density floor, harmonic mobility, geometry, h, current and potential
    are independent equations. Only constant EP transport/mass prefactors are
    calibrated. Dual numbers avoid finite-difference artifacts at exact zeros.
    """
    D = lambda x: Decimal.from_float(float(x))
    n=len(state)
    class Dual:
        def __init__(self,value,gradient=None):
            self.value=value if isinstance(value,Decimal) else Decimal(value)
            self.gradient=gradient if gradient is not None else [Decimal(0)]*n
        def cast(self,value):
            return value if isinstance(value,Dual) else Dual(value)
        def __add__(self,other):
            other=self.cast(other)
            return Dual(self.value+other.value,[a+b for a,b in zip(self.gradient,other.gradient)])
        __radd__=__add__
        def __neg__(self):
            return Dual(-self.value,[-g for g in self.gradient])
        def __sub__(self,other):
            return self+-self.cast(other)
        def __rsub__(self,other):
            return self.cast(other)+-self
        def __mul__(self,other):
            other=self.cast(other)
            return Dual(self.value*other.value,[a*other.value+b*self.value for a,b in zip(self.gradient,other.gradient)])
        __rmul__=__mul__
        def __truediv__(self,other):
            other=self.cast(other)
            return Dual(self.value/other.value,[(a*other.value-b*self.value)/(other.value*other.value)
                for a,b in zip(self.gradient,other.gradient)])
        def __rtruediv__(self,other):
            return self.cast(other)/self
        def ln(self):
            return Dual(self.value.ln(),[g/self.value for g in self.gradient])
    with localcontext() as ctx:
        ctx.prec=800
        y=[]
        for j,value in enumerate(state):
            gradient=[Decimal(0)]*n;gradient[j]=Decimal(1)
            y.append(Dual(D(value),gradient))
        ns,ie=r.num_core_species,r.electron_index
        te0=D(r.Te.value_si)
        tg,pressure,gasR=D(r.T.value_si),D(r.P.value_si),D(constants.R)
        ions=[j for j in range(ns) if j!=ie and r.species_charges[j]>0]
        neutrals=[j for j in range(ns) if r.neutral_heavy_mask[j]]
        anions=[j for j in range(ns) if j!=ie and r.species_charges[j]<0]
        iontemp=tg if r.ambipolar_ion_temperature is not None else D(0)
        floor=D(r.wall_neutral_density_floor)/D(constants.Na)
        kr,kz=map(D,r.wall_diffusion_components)
        wr,wz=kr/(kr+kz),kz/(kr+kz)
        te=y[-1] if r.energy_balance else Dual(te0)
        volume=gasR*(tg*sum(y[j] for j in range(ns) if j!=ie)+te*y[ie])/pressure
        density=sum(y[j] for j in neutrals)/volume
        if density.value<=floor:
            density=Dual(floor)
        def mobility(j):
            if not r.wall_blanc:
                return Dual(1)
            row=list(map(D,r.wall_ion_bath_k0[j]))
            bath=[sum(y[k] for k in neutrals if r.wall_bath_group[k]==g) for g in range(len(row))]
            total=sum(bath[g] for g in range(len(row)) if row[g]>0)
            return total/sum(bath[g]/row[g] for g in range(len(row)) if row[g]>0)
        # Transport calibration also supports deliberate outside-domain trial chemistry.
        components=getattr(r, '_compute_ion_wall_components', r.compute_ion_wall_components)
        ep=components(state,r.compute_volume(state))['total_ep']
        k={j:mobility(j) for j in ions}
        pref={j:D(ep[j])*density.value/(te0+iontemp)/k[j].value for j in ions}
        species={index:sp for sp,index in r.species_index.items()}
        mass={j:D(np.exp(r.energy_ion_sheath_factor[j]-.5) if r.energy_balance else
            np.sqrt(species[j].molecular_weight.value_si/(2*np.pi*constants.m_e))) for j in ions}
        h=y[ie]/(y[ie]+sum(y[j] for j in anions))
        factor=Dual(1) if r.electronegative_wall_model=='o2ReferenceQualifiedUnity' else (
            h if r.electronegative_wall_geometry=='fullFrequency' else wr*h+wz)
        currents={j:pref[j]*k[j]*(te+iontemp)/density*y[j] for j in ions}
        gamma=sum(D(r.species_charges[j])*currents[j] for j in ions)
        plus=sum(D(r.species_charges[j])*y[j] for j in ions)
        mean=sum(mass[j]*currents[j] for j in ions)/gamma
        phi=mean.ln()+(y[ie]/plus).ln()-factor.ln()
        out=[Dual(0) for _ in range(n)]
        for j in ions:
            loss=factor*currents[j]
            out[j]-=loss
            out[ie]-=D(r.species_charges[j])*loss
            target=r.wall_recycle_target[j]
            if target>=0 and r.wall_recycling>0:
                out[target]+=D(r.wall_recycling*r.wall_recycle_multiplicity[j])*loss
        if r.energy_balance:
            out[-1]=-te*(D(1)+phi)*factor*gamma/y[ie]/D('1.5')
        if full:
            concentrations = [amount/volume for amount in y[:ns]]
            chemical_energy, electron_production = Dual(0), Dual(0)
            for j in range(r.num_core_reactions):
                left = [int(v) for v in r.reactant_indices[j] if v >= 0]
                right = [int(v) for v in r.product_indices[j] if v >= 0]
                for reactants, products, coefficient, sign in (
                        (left, right, r.kf[j], 1), (right, left, r.kb[j], -1)):
                    rate = Dual(D(coefficient))*volume
                    for index in reactants:
                        rate *= concentrations[index]
                    electron_production += (products.count(ie)-reactants.count(ie))*rate
                    if r.energy_participates[j]:
                        chemical_energy += sign*rate*(D(r.energy_threshold[j])/(D('1.5')*gasR)
                            +D(r.energy_electrons_consumed[j])*te)
            if r.energy_balance:
                out[-1] -= (chemical_energy+te*electron_production)/y[ie]
                out[-1] += D(r.absorbed_power_reactor)/(D('1.5')*gasR*y[ie])
                # Grid fixtures explicitly ignore elastic loss and have no source.
                assert len(r.energy_elastic_index) == 0
                assert r.ionisation_source.value_si == 0.
        if values:
            return np.array([v.value for v in out],float)
        return np.array([v.gradient for v in out],float)


@pytest.fixture(autouse=True)
def local_import():
    assert Path.cwd() in Path(rmgpy.__file__).resolve().parents
    import rmgpy.data.rmg as data
    saved, data.database = data.database, None
    yield
    data.database = saved


def attachment_model(threshold):
    r, species, reactions = model(energy=True)
    # Charge-conserving mutual neutralisation lets the physical monitor's
    # existing A/B/C interface gates pass at the review's trace-electron state.
    neutralisation = LibraryReaction(reactants=[species[2], species[4]],
        products=[species[1], species[3]], reversible=False,
        kinetics=Arrhenius(A=(1.e12, 'm^3/(mol*s)')), library='EN-test')
    r.electron_energy_balance['electron_energies']['EN-test:2'] = (threshold, 'eV')
    r._configure_energy_balance(copy.deepcopy(r.electron_energy_balance))
    r.initialize_model(species, reactions+[neutralisation], [], [])
    r.kf[0] = 0.
    y = set_state(r, [1.e-30, 1., 1.e-6, 1.e-3, 1.e-6], 3., volume=1.)
    return r, y


@pytest.mark.parametrize('threshold', [0., .3])
def test_attachment_full_energy_row_and_last_valid_are_accurate(threshold):
    r, y = attachment_model(threshold)
    expected = decimal_operator(r, y, full=True)[-1]
    expected_value = decimal_operator(r, y, full=True, values=True)[-1]
    r.monitor_electronegative_wall(y, 7.)
    assert r.electronegative_wall_last_valid_state[0] == 7.
    actual = r.jacobian(7., y, np.zeros_like(y), 0.)[-1]
    print('ATTACHMENT', threshold, 'actual', actual[0], 'Decimal', expected[0],
          'last_valid', r.electronegative_wall_last_valid_state[0])
    np.testing.assert_allclose(actual, expected, rtol=2.e-10, atol=0.)
    assert r.residual(7.,y,np.zeros_like(y))[0][-1] == pytest.approx(
        expected_value, rel=2.e-10, abs=0.)


@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
def test_prescribed_te_complete_jacobian_refuses_before_publication(arm):
    r, species, reactions = model(energy=False, arm=arm)
    n = Species(label='N').from_adjacency_list('1 N u3 p1 c0')
    ncl = Species(label='NCl').from_adjacency_list(
        'multiplicity 3\n1 N u2 p1 c0 {2,S}\n2 Cl u0 p3 c0 {1,S}')
    n.thermo = copy.deepcopy(species[3].thermo)
    ncl.thermo = copy.deepcopy(species[3].thermo)
    reaction = Reaction(reactants=[n, species[3]], products=[ncl], reversible=False,
                        kinetics=Arrhenius(A=(1., 'm^3/(mol*s)')))
    r.initial_mole_fractions.update({n:1.e-29, ncl:0.})
    r.initialize_model(species+[n, ncl], reactions+[reaction], [], [])
    # Start from a representable operator, then install the finite-rate
    # witness to test a later acceptance and preservation of last-valid.
    reaction.kinetics.A.value_si = 1.e308
    r.kf[-1] = 1.e308
    pressure = 55.8*101325./760.
    neutral = pressure/(constants.R*r.T.value_si)-2.-2.e-6-1.e-6*r.Te.value_si/r.T.value_si
    y = set_state(r, [1.e-6, neutral, 1.5e-6, 2., .5e-6, 1.e-30, 0.])
    old = copy.deepcopy(r.electronegative_wall_last_valid_state)
    residual = r.residual(7., y, np.zeros_like(y))[0]
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        operator = r.jacobian(7., y, np.zeros_like(y), 0.)
    assert np.isfinite(residual).all()
    assert np.isinf(operator[6, 5])
    print('PRESCRIBED', arm, 'residual_finite', True, 'J[NCl,N]', operator[6,5])
    with pytest.raises(PlasmaStateError, match='domain:.*rate coefficient'):
        r.monitor_electronegative_wall(y, 7.)
    assert r.electronegative_wall_last_valid_state[0] == old[0]
    np.testing.assert_array_equal(r.electronegative_wall_last_valid_state[1], old[1])


def test_reduced_reverse_energy_keeps_a_representable_product():
    r, y = attachment_model(1.e-300)
    # A finite directed rate can overflow after division by ne even when
    # the threshold-weighted energy and every final derivative are finite.
    # Supply the reverse coefficient explicitly to exercise that direction.
    r.kb[1] = 1.e300
    expected = decimal_operator(r, y, full=True)[-1]
    expected_value = decimal_operator(r, y, full=True, values=True)[-1]
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        actual = r.jacobian(0., y, np.zeros_like(y), 0.)
        value = r.residual(0., y, np.zeros_like(y))[0][-1]
    print('REDUCED_REVERSE', 'actual', actual[-1,0], 'Decimal', expected[0])
    assert np.isfinite(actual).all()
    np.testing.assert_allclose(actual[-1], expected, rtol=2.e-10, atol=0.)
    assert value == pytest.approx(expected_value, rel=2.e-10, abs=0.)


@pytest.mark.parametrize('api', ['compute_ion_wall_components', 'compute_ion_wall_frequencies'])
@pytest.mark.parametrize('energy', [False, True])
@pytest.mark.parametrize('volume', [1.e-13, 1.e4])
def test_public_wall_apis_refuse_volume(api, energy, volume):
    r, _, _ = model(energy=energy)
    y = r.y.copy()
    y[:r.num_core_species] *= volume/r.compute_volume(y)
    with pytest.raises(ElectronegativeWallRegimeError, match='volume'):
        getattr(r, api)(y, volume)


@pytest.mark.parametrize('api', ['compute_ion_wall_components', 'compute_ion_wall_frequencies'])
def test_public_wall_apis_synchronise_supplied_te(api):
    r, _, _ = model(energy=True)
    y = r.y.copy()
    y[-1] *= 2.
    old_te = r.Te.value_si
    r.residual(0., y, np.zeros_like(y))
    expected = getattr(r, api)(y, r.compute_volume(y))
    if isinstance(expected, dict):
        expected = expected['total_ep']
    r.Te.value_si = old_te
    actual = getattr(r, api)(y, r.compute_volume(y))
    if isinstance(actual, dict):
        actual = actual['total_ep']
    print('PUBLIC_TE', api, 'actual', actual[2], 'synchronised', expected[2])
    assert r.Te.value_si == y[-1]
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize('api', ['compute_ion_wall_components', 'compute_ion_wall_frequencies'])
@pytest.mark.parametrize('te_ev, held_ev', [(.05, 50.), (50., .05)])
def test_public_te_edges_validate_the_current_eos(api, te_ev, held_ev):
    r, _, _ = model(energy=True)
    y = set_state(r, [100.,1.e-30,100.,.001,1.e-30], te_ev, volume=1.)
    expected = getattr(r, api)(y, r.compute_volume(y))
    if isinstance(expected, dict):
        expected = expected['total_ep']
    r.Te.value_si = held_ev*constants.e*constants.Na/constants.R
    stale_volume = r.compute_volume(y)
    actual = getattr(r, api)(y, stale_volume)
    if isinstance(actual, dict):
        actual = actual['total_ep']
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize('blanc', [False, True])
@pytest.mark.parametrize('energy', [False, True])
@pytest.mark.parametrize('mixture', [False, True])
@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
@pytest.mark.parametrize('closure', ['confinedAnion', 'o2ReferenceQualifiedUnity'])
@pytest.mark.parametrize('ion_temperature', [None, 'gas'])
def test_independent_volume_unequal_cations_full_energy_grid(
        blanc, energy, mixture, arm, closure, ion_temperature):
    overrides = blanc_overrides(mixture) if blanc else {}
    if energy:
        labels = ('Ar','Cl','Ne') if mixture else ('Ar','Cl')
        overrides['electron_energy_balance'] = dict(absorbed_power=(0.,'W'),
            chamber_volume=(np.pi*.05**2*.30,'m^3'), sheath='floating_wall',
            elastic_collisions={label:{'ignore':'grid isolation'} for label in labels},
            electron_energies={'EN-test:1':(0.,'eV'),'EN-test:2':(.3,'eV')})
    r, _, _ = model(energy=energy, arm=arm, closure=closure, second_cation=mixture,
        ambipolar_ion_temperature=ion_temperature, **overrides)
    r.kf[0] = 0. # Isolate attachment; a detachment can hide its cancellation.
    worst, location, count = 0., None, 0
    for ce, cm, ratio, te, neutral, volume in product(
            [1.e-30,1.e-6,1.e3], [1.e-30,1.e-6,1.e3], [1.e-8,1.,1.e8],
            [.05,3.,50.], [1.e-30,1.,1.e6], [1.e-12,1.,1.e3]):
        # At the maximum total charge only an equal mixture would fit, so
        # omit that point rather than silently lose the unequal split.
        if (ce+cm)*(0.7 if mixture else 1.) > 1.e3:
            continue
        if mixture and .3*(ce+cm) < 1.e-30:
            continue # At this lower corner only an equal split is in-domain.
        trace = max(1.e-30, neutral*1.e-30) if blanc else .001
        concentrations = [ce, neutral, ce+cm, trace, cm]
        if mixture:
            concentrations[2] *= .3
            concentrations.extend([0., .7*(ce+cm)])
        y = set_state(r, concentrations, te, volume=volume, ratio=ratio)
        expected = decimal_operator(r, y, full=energy)
        actual = r.compute_electronegative_wall_jacobian(y)
        if energy:
            actual[-1] = r.jacobian(0., y, np.zeros_like(y), 0.)[-1]
        # Species rows compare wall terms; the last row includes chemistry.
        nonzero = expected != 0.
        error = float(np.max(np.abs((actual[nonzero]-expected[nonzero])/expected[nonzero]), initial=0.))
        if error > worst:
            worst, location = error, [ce,cm,ratio,te,neutral,volume]
        count += 1
        np.testing.assert_allclose(actual, expected, rtol=2.e-10, atol=0.)
    print('FULL_DOMAIN_SWEEP', blanc, energy, mixture, arm, closure, ion_temperature,
          'cases', count, 'worst', worst, 'state', location)
