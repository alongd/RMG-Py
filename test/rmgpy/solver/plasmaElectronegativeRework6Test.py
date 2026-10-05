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

"""Material cancellation, product and public-state regressions for EN walls."""
from decimal import Decimal, localcontext
from itertools import product, permutations
from pathlib import Path
import copy
import json
import numpy as np
import pytest
import rmgpy
from rmgpy import constants
from rmgpy.kinetics import Arrhenius, TwoTemperaturePlasma
from rmgpy.reaction import Reaction
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.species import Species
from rmgpy.exceptions import ElectronegativeWallRegimeError, PlasmaStateError
from plasmaElectronegativeWallTest import model, FixtureEnergyReactor
from plasmaElectronegativeRework4Test import set_state, blanc_overrides
from plasmaElectronegativeRework5Test import attachment_model
from plasmaEnergyBalanceTest import LL_AR_ELASTIC

def decimal_operator(r, state, full=False, values=False, precision=800):
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
        def exp(self):
            value = self.value.exp()
            return Dual(value, [g*value for g in self.gradient])
    with localcontext() as ctx:
        ctx.prec=precision
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
            inverse=sum(bath[g]/row[g] for g in range(len(row)) if row[g]>0)
            # Differentiate the harmonic mean from bath differences. At a
            # single populated bath its own column is mathematically zero;
            # finite Decimal division must not invent a 10**-precision entry.
            gradient = [Decimal(0)]*n
            for g in range(len(row)):
                if row[g] <= 0:
                    continue
                difference = sum(bath[h].value*(1/row[h]-1/row[g])
                    for h in range(len(row)) if row[h]>0)
                for column in range(n):
                    gradient[column] += difference/inverse.value**2*bath[g].gradient[column]
            return Dual(total.value/inverse.value, gradient)
        # Calibrate transport independently of deliberately out-of-domain
        # chemical caches; the public API now refuses those caches.
        components = getattr(r, '_compute_ion_wall_components', r.compute_ion_wall_components)
        ep=components(state,r.compute_volume(state))['total_ep']
        k={j:mobility(j) for j in ions}
        pref={j:D(ep[j])*density.value/(te0+iontemp)/k[j].value for j in ions}
        species={index:sp for sp,index in r.species_index.items()}
        mass={j:D(np.exp(r.energy_ion_sheath_factor[j]-.5) if r.energy_balance else
            np.sqrt(species[j].molecular_weight.value_si/(2*np.pi*constants.m_e))) for j in ions}
        h=y[ie]/(y[ie]+sum(y[j] for j in anions))
        factor=Dual(1) if r.electronegative_wall_model=='electropositiveBracket' else (
            h if r.electronegative_wall_geometry=='fullFrequency' else wr*h+wz)
        currents={j:pref[j]*k[j]*(te+iontemp)/density*y[j] for j in ions}
        gamma=sum(D(r.species_charges[j])*currents[j] for j in ions)
        plus=sum(D(r.species_charges[j])*y[j] for j in ions)
        numerator=sum(mass[j]*currents[j] for j in ions)
        mean_gradient = [Decimal(0)]*n
        for j in ions:
            difference = sum((mass[j]-mass[k])*currents[k].value for k in ions)
            for column in range(n):
                mean_gradient[column] += difference/gamma.value**2*currents[j].gradient[column]
        mean=Dual(numerator.value/gamma.value, mean_gradient)
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
            laws = dict(r.energy_te_rate_refresh or [])
            for j in range(r.num_core_reactions):
                left = [int(v) for v in r.reactant_indices[j] if v >= 0]
                right = [int(v) for v in r.product_indices[j] if v >= 0]
                for reactants, products, coefficient, sign in (
                        (left, right, r.kf[j], 1), (right, left, r.kb[j], -1)):
                    coefficient = Dual(D(coefficient))
                    if sign == 1 and j in laws:
                        kin = laws[j]
                        assert type(kin).__name__ == 'TwoTemperaturePlasma'
                        # Calibrate at the oracle's sample temperature, independently
                        # of the reactor cache. set_state changes public Te before
                        # the production operator refreshes its evaluated rates.
                        coefficient = Dual(D(kin.get_rate_coefficient_two_temp(
                            r.T.value_si, float(te0))))
                        exponent, activation = D(kin.n.value_si), D(kin.Ea_e.value_si)
                        if exponent == -1:
                            coefficient *= te0/te
                        elif exponent != 0:
                            coefficient *= ((te/te0).ln()*exponent).exp()
                        if activation != 0:
                            coefficient *= (activation/gasR*(1/te0-1/te)).exp()
                    rate = coefficient*volume
                    for index in reactants:
                        rate *= concentrations[index]
                    for index in reactants:
                        out[index] -= rate
                    for index in products:
                        out[index] += rate
                    electron_production += (products.count(ie)-reactants.count(ie))*rate
                    if r.energy_balance and r.energy_participates[j]:
                        chemical_energy += sign*rate*(D(r.energy_threshold[j])/(D('1.5')*gasR)
                            +D(r.energy_electrons_consumed[j])*te)
            ionisable = [j for j in neutrals if r.source_cation_target[j] >= 0]
            source = D(r.ionisation_source.value_si)/D(constants.Na)*volume
            total = sum(y[j] for j in ionisable)
            if not ionisable or total.value == 0:
                source = Dual(0)
            else:
                for j in ionisable:
                    rate = source*y[j]/total
                    out[j] -= rate
                    out[r.source_cation_target[j]] += rate
                    out[ie] += rate
            if r.energy_balance:
                out[-1] -= (chemical_energy+te*electron_production)/y[ie]
                out[-1] += D(r.absorbed_power_reactor)/(D('1.5')*gasR*y[ie])
                out[-1] -= te*source/y[ie]
                tev = te*gasR/D(constants.Na)/D(constants.e)
                lt = tev.ln()
                for p, partner in enumerate(r.energy_elastic_index):
                    a,b,c,d = map(D, r.energy_elastic_params[p])
                    k = a*(b*lt+c*lt*lt+d*lt*lt*lt).exp()
                    out[-1] -= 2*D(r.energy_elastic_mass_ratio[p])*k*D(constants.Na)*concentrations[partner]*(te-tg)
        if values == 'both':
            return np.array([v.value for v in out],float), np.array([v.gradient for v in out],float)
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


def monitor_with_parameter_domain(r, y, t):
    """Retain old trial-arithmetic witnesses, refuse their impossible parameters."""
    if np.any(r.kf > 1.e15) or np.any(r.kb > 1.e15):
        previous = copy.deepcopy(r.electronegative_wall_last_valid_state)
        with pytest.raises(PlasmaStateError, match='domain:.*rate coefficient'):
            r.monitor_electronegative_wall(y, t)
        assert r.electronegative_wall_last_valid_state[0] == previous[0]
        np.testing.assert_array_equal(r.electronegative_wall_last_valid_state[1], previous[1])
        return
    r.monitor_electronegative_wall(y, t)


def reinitialise(r, species, reactions):
    if r.energy_balance:
        config = copy.deepcopy(r.electron_energy_balance)
        config['elastic_collisions'].update({sp.label: {'ignore': 'test isolation'}
            for sp in species if sp.label not in config['elastic_collisions']
            and sp.get_net_charge() == 0})
        for i in range(len(reactions)):
            config['electron_energies'].setdefault('EN-test:'+str(i+1), (0., 'eV'))
        FixtureEnergyReactor.entries = {i+1: rx.kinetics for i, rx in enumerate(reactions)}
        r._configure_energy_balance(config)
    r.initialize_model(species, reactions, [], [])


def add_neutral_chemistry(energy=False, order=(0, 1, 2), reverse=False):
    r, sp, reactions = model(energy=energy)
    n = Species(label='N').from_adjacency_list('1 N u3 p1 c0')
    ncl = Species(label='NCl').from_adjacency_list(
        'multiplicity 3\n1 N u2 p1 c0 {2,S}\n2 Cl u0 p3 c0 {1,S}')
    n.thermo, ncl.thermo = copy.deepcopy(sp[3].thermo), copy.deepcopy(sp[3].thermo)
    sp += [n, ncl]
    r.initial_mole_fractions.update({n: 1.e-29, ncl: 0.})
    left = [sp[1], n, sp[3]]
    rx = LibraryReaction(reactants=[left[i] for i in order], products=[sp[1], ncl],
        reversible=False, kinetics=Arrhenius(A=(1., 'm^6/(mol^2*s)')), library='EN-test')
    if reverse:
        rx.reactants, rx.products = rx.products, rx.reactants
        rx.kinetics = Arrhenius(A=(1., 'm^3/(mol*s)'))
    reactions.append(rx)
    reinitialise(r, sp, reactions)
    r.kf[-1], r.kb[-1] = (0., 1.e308) if reverse else (1.e308, 0.)
    return r, sp, reactions


@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
@pytest.mark.parametrize('activation', [0., 1.e-20])
def test_te_energy_constant_product_cancels_before_accumulation(arm, activation):
    r, sp, reactions = model(energy=True, arm=arm)
    reactions.append(LibraryReaction(reactants=[sp[2], sp[4]], products=[sp[1], sp[3]],
        reversible=False, kinetics=Arrhenius(A=(1.e12, 'm^3/(mol*s)')), library='EN-test'))
    rx = LibraryReaction(reactants=[sp[0], sp[1]], products=[sp[2], sp[0], sp[0]],
        reversible=False, library='EN-test', kinetics=TwoTemperaturePlasma(
            A=(1., 'm^3/(mol*s)'), n=-1., Ea_g=(0., 'J/mol'),
            Ea_e=(activation, 'J/mol'), T0=(30000., 'K'), electrons=1))
    reactions.append(rx)
    reinitialise(r, sp, reactions)
    r.kf[0] = 0.
    rx.kinetics.A.value_si = 1.e20
    r.kf[-1] = r._evaluate_two_temperature_rate_coefficient(rx.kinetics)
    y = set_state(r, [1.e-30, 1., 1.e-6, .001, 1.e-6])
    expected = decimal_operator(r, y, full=True)
    monitor_with_parameter_domain(r, y, 7.)
    actual = r.jacobian(7., y, np.zeros_like(y), 0.)
    print('TE_ENERGY', arm, activation, 'actual', actual[-1,-1], 'Decimal', expected[-1,-1],
          'last_valid', r.electronegative_wall_last_valid_state[0])
    assert r.electronegative_wall_last_valid_state[0] == (0. if np.any(r.kf > 1.e15) or np.any(r.kb > 1.e15) else 7.)
    np.testing.assert_allclose(actual[-1], expected[-1], rtol=2.e-10, atol=0.)
    np.testing.assert_allclose(r.residual(7., y, np.zeros_like(y))[0],
        decimal_operator(r, y, full=True, values=True), rtol=2.e-10, atol=0.)


@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
def test_elastic_partner_keeps_other_pressure(arm):
    r, sp, reactions = model(energy=True, arm=arm)
    config = copy.deepcopy(r.electron_energy_balance)
    config['elastic_collisions']['Ar'] = copy.deepcopy(LL_AR_ELASTIC)
    r._configure_energy_balance(config)
    r.initialize_model(sp, reactions, [], [])
    r.kf[1] = 1.e35
    y = set_state(r, [1.e-30, 1., 2.e-30, 1.e-30, 1.e-30])
    expected = decimal_operator(r, y, full=True)[-1]
    monitor_with_parameter_domain(r, y, 7.)
    actual = r.jacobian(7., y, np.zeros_like(y), 0.)[-1]
    print('ELASTIC', arm, 'actual', actual[1], 'Decimal', expected[1], 'last_valid', 7.)
    np.testing.assert_allclose(actual, expected, rtol=2.e-10, atol=0.)


@pytest.mark.parametrize('energy', [False, True])
@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('order', list(permutations(range(3))))
def test_species_rate_is_reactant_order_independent(energy, reverse, order):
    r, _, _ = add_neutral_chemistry(energy, order, reverse)
    r.kf[1] = 1.e35 # Keep the existing Gate B qualified at trace Cl.
    y = set_state(r, [1.e-6, 1.e6, 1.5e-6, 1.e-30, .5e-6, 1.e-30, 0.])
    actual = r.residual(7., y, np.zeros_like(y))[0]
    print('PRODUCT', energy, reverse, order, 'NCl', actual[6], 'expected', 1.e254)
    assert np.isfinite(actual).all()
    assert actual[6] == pytest.approx(1.e254, rel=2.e-15)
    assert np.isfinite(r.jacobian(7., y, np.zeros_like(y), 0.)).all()
    monitor_with_parameter_domain(r, y, 7.)
    assert r.electronegative_wall_last_valid_state[0] == (0. if np.any(r.kf > 1.e15) or np.any(r.kb > 1.e15) else 7.)


@pytest.mark.parametrize('energy', [False, True])
def test_network_leak_product_preserves_finite_rate(energy):
    r, _, _ = add_neutral_chemistry(energy)
    r.kf[-1] = 0.
    r.network_indices = np.array([[1, 5, 3]], dtype=int)
    r.network_leak_coefficients = np.array([1.e308])
    r.network_leak_rates = np.zeros(1)
    y = set_state(r, [1.e-6, 1.e6, 1.5e-6, 1.e-30, .5e-6, 1.e-30, 0.])
    r.residual(0., y, np.zeros_like(y))
    print('NETWORK', energy, r.network_leak_rates[0])
    assert r.network_leak_rates[0] == pytest.approx(1.e254, rel=2.e-15)


@pytest.mark.parametrize('api', ['compute_anion_transport_data', 'compute_anion_transport_frequencies'])
@pytest.mark.parametrize('energy', [False, True])
@pytest.mark.parametrize('volume', [1.e-13, 1.e4])
def test_public_anion_apis_refuse_volume(api, energy, volume):
    r, _, _ = model(energy=energy)
    y = r.y.copy()
    y[:r.num_core_species] *= volume/r.compute_volume(y)
    with pytest.raises(ElectronegativeWallRegimeError, match='volume'):
        getattr(r, api)(y, volume)


@pytest.mark.parametrize('api', ['compute_anion_transport_data', 'compute_anion_transport_frequencies'])
@pytest.mark.parametrize('te_ev, held_ev', [(6., 3.), (.05, 50.), (50., .05)])
def test_public_anion_apis_synchronise_te(api, te_ev, held_ev):
    r, _, _ = model(energy=True)
    y = set_state(r, [1.e-6, 1., 1.5e-6, .001, .5e-6], te_ev)
    expected = getattr(r, api)(y, r.compute_volume(y))
    r.Te.value_si = held_ev*constants.e*constants.Na/constants.R
    actual = getattr(r, api)(y, r.compute_volume(y))
    print('ANION_TE', api, te_ev, 'actual', actual, 'expected', expected)
    assert r.Te.value_si == y[-1]
    assert actual == expected


@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('energy', [False, True])
def test_neutral_chemistry_keeps_other_pressure(energy, reverse):
    r, sp, reactions = add_neutral_chemistry(energy)
    reactions[-1].reactants = [sp[5], sp[3]]
    reactions[-1].products = [sp[6]]
    reactions[-1].kinetics = Arrhenius(A=(1., 'm^3/(mol*s)'))
    if reverse:
        reactions[-1].reactants, reactions[-1].products = reactions[-1].products, reactions[-1].reactants
        reactions[-1].kinetics = Arrhenius(A=(1., 's^-1'))
    reinitialise(r, sp, reactions)
    r.kf[-1], r.kb[-1] = (0., 1.e9) if reverse else (1.e9, 0.)
    y = set_state(r, [1.e-30, 1.e-30, 2.e-30, 1., 1.e-30, 1.e-30, 0.])
    expected = decimal_operator(r, y, full=True)[6]
    monitor_with_parameter_domain(r, y, 7.)
    actual = r.jacobian(7., y, np.zeros_like(y), 0.)[6]
    print('CHEM_EOS', energy, reverse, 'actual', actual[3], 'Decimal', expected[3], 'last_valid', 7.)
    np.testing.assert_allclose(actual, expected, rtol=2.e-10, atol=0.)


@pytest.mark.parametrize('energy', [False, True])
def test_source_partition_keeps_trace_pressure(energy):
    overrides = {}
    if energy:
        overrides['electron_energy_balance'] = dict(absorbed_power=(0., 'W'),
            chamber_volume=(np.pi*.05**2*.30, 'm^3'), sheath='floating_wall',
            elastic_collisions={label: {'ignore': 'test isolation'} for label in ('Ar', 'Cl', 'Ne')},
            electron_energies={'EN-test:1': (0., 'eV'), 'EN-test:2': (0., 'eV')})
    r, _, _ = model(energy=energy, second_cation=True, ionisation_source=(1.e15, 'm^-3*s^-1'), **overrides)
    y = set_state(r, [1.e-30, 1., 2.e-30, 1.e-30, 1.e-30, 1.e-30, 0.])
    expected = decimal_operator(r, y, full=True)[6,1]
    actual = r.jacobian(0., y, np.zeros_like(y), 0.)[6,1]
    print('SOURCE_PARTITION', energy, 'actual', actual, 'Decimal', expected)
    assert actual == pytest.approx(expected, rel=2.e-10, abs=0.)


@pytest.mark.parametrize('energy', [False, True])
def test_gate_attachment_multiplicity_uses_protected_product(energy):
    r, sp, reactions = model(energy=energy)
    cl2 = Species(label='Cl2').from_smiles('ClCl')
    cl2.thermo = copy.deepcopy(sp[3].thermo)
    r.initial_mole_fractions[cl2] = 1.e-29
    rx = LibraryReaction(reactants=[sp[0], sp[0], cl2], products=[sp[4], sp[4]],
        reversible=False, kinetics=Arrhenius(A=(1., 'm^6/(mol^2*s)')), library='EN-test')
    sp.append(cl2)
    reactions.append(rx)
    reinitialise(r, sp, reactions)
    r.kf[-1] = 1.e308
    y = set_state(r, [1.e-30, 1., 1.e-6, 1., 1.e-6, 1.e-30])
    monitor_with_parameter_domain(r, y, 7.)
    print('GATE_MULTIPLICITY', energy, 'last_valid', r.electronegative_wall_last_valid_state[0])
    assert r.electronegative_wall_last_valid_state[0] == (0. if np.any(r.kf > 1.e15) or np.any(r.kb > 1.e15) else 7.)


@pytest.mark.parametrize('energy', [False, True])
def test_total_electron_source_differentiates_total_without_partition_cancellation(energy):
    overrides = {}
    if energy:
        overrides['electron_energy_balance'] = dict(absorbed_power=(0., 'W'),
            chamber_volume=(np.pi*.05**2*.30, 'm^3'), sheath='floating_wall',
            elastic_collisions={label: {'ignore': 'test isolation'} for label in ('Ar', 'Cl', 'Ne')},
            electron_energies={'EN-test:1': (0., 'eV'), 'EN-test:2': (0., 'eV')})
    r, _, _ = model(energy=energy, second_cation=True, ionisation_source=(1.e6, 'm^-3*s^-1'), **overrides)
    y = set_state(r, [1.e-30, 1.e-30, 0., 1., 1.e-30, 0., 2.e-30], 3., volume=1.e-12)
    expected = decimal_operator(r, y, full=True)[0,5]
    monitor_with_parameter_domain(r, y, 7.)
    actual = r.jacobian(7., y, np.zeros_like(y), 0.)[0,5]
    print('TOTAL_ELECTRON_SOURCE', energy, 'actual', actual, 'Decimal', expected, 'last_valid', 7.)
    assert actual == pytest.approx(expected, rel=2.e-10, abs=0.)


def test_species_te_rate_and_eos_cancel_before_accumulation():
    r, sp, reactions = model(energy=True)
    rx = LibraryReaction(reactants=[sp[0], sp[1]], products=[sp[2], sp[0], sp[0]],
        reversible=False, library='EN-test', kinetics=TwoTemperaturePlasma(
            A=(1., 'm^3/(mol*s)'), n=1., Ea_g=(0., 'J/mol'),
            Ea_e=(0., 'J/mol'), T0=(30000., 'K'), electrons=1))
    reactions.append(rx)
    reinitialise(r, sp, reactions)
    rx.kinetics.A.value_si = 1.e290
    r.kf[-1] = r._evaluate_two_temperature_rate_coefficient(rx.kinetics)
    y = set_state(r, [1.e3, 1.e-30, 1.e-30, 1.e-30, 1.e-30])
    expected = decimal_operator(r, y, full=True)[2,-1]
    actual = r.jacobian(0., y, np.zeros_like(y), 0.)[2,-1]
    print('SPECIES_TE', 'actual', actual, 'Decimal', expected)
    assert actual == pytest.approx(expected, rel=2.e-10, abs=0.)


@pytest.mark.parametrize('blanc', [False, True])
@pytest.mark.parametrize('energy', [False, True])
@pytest.mark.parametrize('mixture', [False, True])
@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
@pytest.mark.parametrize('closure', ['confinedAnion', 'electropositiveBracket'])
@pytest.mark.parametrize('ion_temperature', [None, 'gas'])
def test_complete_population_faces_domain_sweep(
        blanc, energy, mixture, arm, closure, ion_temperature, record_property):
    overrides = blanc_overrides(mixture) if blanc else {}
    if energy:
        labels = ('Ar', 'Cl', 'Ne') if mixture else ('Ar', 'Cl')
        overrides['electron_energy_balance'] = dict(absorbed_power=(1.e-12, 'W'),
            chamber_volume=(np.pi*.05**2*.30, 'm^3'), sheath='floating_wall',
            elastic_collisions={label: (copy.deepcopy(LL_AR_ELASTIC) if label == 'Ar'
                else {'ignore': 'no elastic data in fixture'}) for label in labels},
            electron_energies={'EN-test:1': (0., 'eV'), 'EN-test:2': (.3, 'eV')})
        overrides['ionisation_source'] = (1.e6, 'm^-3*s^-1')
    r, sp, reactions = model(energy=energy, arm=arm, closure=closure,
        second_cation=mixture, ambipolar_ion_temperature=ion_temperature, **overrides)
    if energy:
        reactions.append(LibraryReaction(reactants=[sp[0], sp[1]], products=[sp[2], sp[0], sp[0]],
            reversible=False, library='EN-test', kinetics=TwoTemperaturePlasma(
                A=(1.e9, 'm^3/(mol*s)'), n=-1., Ea_g=(0., 'J/mol'),
                Ea_e=(0., 'J/mol'), T0=(30000., 'K'), electrons=1)))
        reinitialise(r, sp, reactions)
    r.kf[0] = 0.
    result = dict(configuration=dict(blanc=blanc, energy=energy, mixture=mixture,
        arm=arm, closure=closure, ion_temperature=ion_temperature), states=0,
        columns={}, energy_columns={}, residual_rows={}, population_faces={})
    labels = [spc.label for spc in sp]+(['Te'] if energy else [])

    def retain(bucket, label, actual, expected, location):
        nonzero = expected != 0.
        worst = float(np.max(np.abs((actual[nonzero]-expected[nonzero])/expected[nonzero]), initial=0.))
        entry = bucket.setdefault(label, dict(worst_relative_error=0., nonzero_entries=0,
            zero_entries=0, nonzero_at_exact_zero=0, location=None))
        entry['nonzero_entries'] += int(np.count_nonzero(nonzero))
        entry['zero_entries'] += int(np.count_nonzero(~nonzero))
        entry['nonzero_at_exact_zero'] += int(np.count_nonzero(actual[~nonzero]))
        if worst > entry['worst_relative_error']:
            entry['worst_relative_error'], entry['location'] = worst, location

    charges = [(1.e-30, 1.e-30), (1.e-30, 1.e3), (1.e3, 1.e-30), (1.e3, 1.e3), (1.e-6, 1.e-6)]
    for (ce, cm), ratio, te, neutral, volume in product(charges,
            [1.e-8, 1., 1.e8], [.05, 3., 50.], [0., 1.e-30, 1., 1.e6], [1.e-12, 1., 1.e3]):
        total = ce+cm
        fractions = sorted(set([0., .3, .5, 1.e-30/total, 1.])) if mixture else [1.]
        for fraction in fractions:
            populations = [total*fraction, total*(1-fraction)] if mixture else [total]
            if any(population > 1.e3 or 0. < population < 1.e-30*(1-4*np.finfo(float).eps)
                    for population in populations):
                continue
            trace = max(1.e-30, neutral*1.e-30) if blanc else .001
            concentrations = [ce, neutral, populations[0], trace, cm]
            if mixture:
                concentrations += [0., populations[1]]
            y = set_state(r, concentrations, te, volume=volume, ratio=ratio)
            # Refresh the Te-dependent coefficient at each supplied state.
            # set_state sets Te directly; force the ordinary accepted refresh.
            if energy:
                r.Te.value_si *= .99
            actual_value = r.residual(0., y, np.zeros_like(y))[0]
            actual = r.jacobian(0., y, np.zeros_like(y), 0.)
            expected_value, expected = decimal_operator(r, y, full=True, values='both', precision=200)
            location = dict(state=[ce, cm, ratio, te, neutral, volume], fraction=fraction)
            for column, label in enumerate(labels):
                retain(result['columns'], label, actual[:,column], expected[:,column], location)
                if energy:
                    retain(result['energy_columns'], label, actual[-1:,column], expected[-1:,column], location)
                retain(result['residual_rows'], label, actual_value[column:column+1],
                    expected_value[column:column+1], location)
            result['states'] += 1
            face = ('zero' if fraction in (0., 1.) else 'equal' if fraction == .5
                    else '30/70' if fraction == .3 else 'minority') if mixture else 'single'
            result['population_faces'][face] = result['population_faces'].get(face, 0)+1
            np.testing.assert_allclose(actual, expected, rtol=2.e-10, atol=0.)
            np.testing.assert_allclose(actual_value, expected_value, rtol=2.e-10, atol=0.)
    record_property('i313_sweep', json.dumps(result))
    print('COMPLETE_DOMAIN_SWEEP', json.dumps(result))
