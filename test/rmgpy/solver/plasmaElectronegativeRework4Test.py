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

"""Physical EN domain, exact chemistry and independently differentiated wall law."""
from decimal import Decimal, localcontext
from itertools import product
from pathlib import Path
import copy

import numpy as np
import pytest
import rmgpy
from rmgpy import constants
from rmgpy.exceptions import ElectronegativeWallRegimeError, PlasmaStateError
from rmgpy.kinetics import Arrhenius
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.species import Species
from plasmaElectronegativeWallTest import model, FixtureEnergyReactor


@pytest.fixture(autouse=True)
def local_import():
    assert Path.cwd() in Path(rmgpy.__file__).resolve().parents
    import rmgpy.data.rmg as data
    saved, data.database = data.database, None
    yield
    data.database = saved


def set_state(r, concentrations, te_ev=3., volume=1., ratio=1.):
    te = te_ev*constants.e*constants.Na/constants.R
    r.Te.value_si = te
    r.P.value_si = constants.R*((np.sum(concentrations)-concentrations[0])*r.T.value_si+te*concentrations[0])
    kr = 1.e3/(1.+ratio)
    r.wall_diffusion_components = (kr, kr*ratio)
    r.diffusion_length.value_si = 1./np.sqrt(sum(r.wall_diffusion_components))
    y = np.array(list(np.array(concentrations)*volume)+([te] if r.energy_balance else []))
    r.P.value_si=constants.R*((np.sum(y[:r.num_core_species])-y[0])*r.T.value_si+te*y[0])/volume
    return y


def blanc_overrides(mixture=False):
    unit='m^2/(V*s)'
    ar={'Ar':(1.535e-4,unit),'Cl':(3.07e-4,unit)}
    anion={'Ar':(1.5e-4,unit),'Cl':(3.e-4,unit)}
    ions={'Ar+':{'perBath':ar}}
    if mixture:
        ar['Ne']=(.7675e-4,unit)
        anion['Ne']=(.75e-4,unit)
        ions['Ne+']={'perBath':{'Ar':(3.07e-4,unit),'Cl':(1.535e-4,unit),'Ne':(1.535e-4,unit)}}
    return dict(ion_reduced_mobilities=ions,anion_reduced_mobilities={'Cl-':{'perBath':anion}},
                wall_single_bath_approximation=False)


def decimal_wall(r, state):
    """Differentiate the complete constitutive law with 450-digit dual numbers.

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
        ctx.prec=450
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
        ep=r.compute_ion_wall_components(state,r.compute_volume(state))['total_ep']
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
        return np.array([v.gradient for v in out],float)


@pytest.mark.parametrize('arm',['fullFrequency','radialOnly'])
def test_quadratic_trace_mass_action_is_exact(arm):
    r,sp,rxns=model(arm=arm,energy=True)
    n=Species(label='N').from_adjacency_list('1 N u3 p1 c0')
    n2=Species(label='N2').from_smiles('N#N')
    n.thermo=copy.deepcopy(sp[3].thermo); n2.thermo=copy.deepcopy(sp[3].thermo)
    rx=LibraryReaction(reactants=[n,n],products=[n2],reversible=False,
                       kinetics=Arrhenius(A=(1.e12,'m^3/(mol*s)')),library='EN-test')
    r.initial_mole_fractions[n]=1.e-29; r.initial_mole_fractions[n2]=0.
    r.electron_energy_balance['elastic_collisions'].update(N={'ignore':'test'},N2={'ignore':'test'})
    r.electron_energy_balance['electron_energies']['EN-test:3']=(0.,'eV')
    r._configure_energy_balance(copy.deepcopy(r.electron_energy_balance))
    FixtureEnergyReactor.entries[3]=rx.kinetics
    r.initialize_model(sp+[n,n2],rxns+[rx],[],[])
    y=r.y.copy(); y[5]=1.e-30
    v=r.compute_volume(y)
    expected=2.e12*y[5]/v-1.e12*(y[5]/v)**2*constants.R*r.T.value_si/r.P.value_si
    got=r.jacobian(0.,y,np.zeros_like(y),0.)[6,5]
    print('QUADRATIC_TRACE',arm,'actual',got,'exact',expected)
    assert got == pytest.approx(expected,rel=1.e-13,abs=0.), (got,expected)


@pytest.mark.parametrize('arm',['fullFrequency','radialOnly'])
def test_linear_trace_chemistry_remains_exact(arm):
    r,_,_=model(arm=arm,energy=True)
    y=r.y.copy(); y[1],y[3]=1.e-30,1.
    v=r.compute_volume(y)
    expected=r.kf[1]*y[0]*y[3]/v**2*constants.R*r.T.value_si/r.P.value_si
    wall=r._en_wall_linearization(y,v)[1][0,1]
    got=r.jacobian(0.,y,np.zeros_like(y),0.)[0,1]-wall
    print('LINEAR_TRACE',arm,'actual',got,'exact',expected)
    assert got == pytest.approx(expected,rel=1.e-13,abs=0.)


def test_small_axial_share_matches_decimal():
    r,_,_=model(energy=True,arm='radialOnly')
    y=set_state(r,[1.e-30,1.,1.e-6,0.001,1.e-6],ratio=1.e-8)
    expected=decimal_wall(r,y)
    actual=r.compute_electronegative_wall_jacobian(y)
    np.testing.assert_allclose(actual,expected,rtol=2.e-10,atol=0.)


def test_intermediate_products_preserve_domain_derivatives():
    r,_,_=model(energy=True)
    # Concentrations alone cannot bound floating-point intermediates: the
    # EOS also carries an inventory volume. Check its two physical edges,
    # then refuse an absurd normalization before publishing any derivative.
    for volume in (1.e-12,1.e3):
        y=set_state(r,[1.e-30,1.e6,2.e-30,0.001,1.e-30],volume=volume)
        expected=decimal_wall(r,y)
        actual=r.compute_electronegative_wall_jacobian(y)
        np.testing.assert_allclose(actual,expected,rtol=2.e-10,atol=0.)
    y=set_state(r,[1.e-30,1.e6,2.e-30,0.001,1.e-30],volume=1.e290)
    with pytest.raises(ElectronegativeWallRegimeError,match='volume'):
        r.compute_electronegative_wall_jacobian(y)


def test_monitor_dilution_guard_is_complete_between_checks():
    r,_,_=model(energy=True)
    old=copy.deepcopy(r.electronegative_wall_last_valid_state)
    # The old range heuristic ignores a large chemical electron-production
    # rate, whose energy dilution derivative is unrepresentable.
    y=set_state(r,[1.e-30,1.,1.e-6,0.001,1.e-6])
    r.kf[0]=1.e270
    residual=r.residual(1.,y,np.zeros_like(y))[0]
    with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
        operator=r.jacobian(1.,y,np.zeros_like(y),0.)
    assert np.isfinite(residual).all()
    assert not np.isfinite(operator).all()
    print('DILUTION_STATE','residual_finite',True,'Jacobian_finite',False)
    with pytest.raises(PlasmaStateError,match='domain:.*rate coefficient'):
        r.monitor_electronegative_wall(y,1.)
    assert r.electronegative_wall_last_valid_state[0] == old[0]
    np.testing.assert_array_equal(r.electronegative_wall_last_valid_state[1],old[1])


@pytest.mark.parametrize('frequency',[1.e-13,1.e-12,1.e12,1.e13])
def test_absolute_wall_frequency_domain_edges(frequency):
    r,_,_=model(energy=True)
    # Reach the frequency endpoints while keeping the live mobility,
    # concentration, lengths and eigenvalues inside the parameter domain.
    neutral = 1.e5 if frequency < 1. else 1.
    r.wall_ion_mobility_si[2] = 1.e-8 if frequency < 1. else 1.e-2
    y=set_state(r,[1.e-30,neutral,2.e-30,.001,1.e-30])
    base=r.compute_ion_wall_components(y,r.compute_volume(y))['total_ep'][2]
    scale=frequency/base
    r.wall_diffusion_components=tuple(k*scale for k in r.wall_diffusion_components)
    r.diffusion_length.value_si/=np.sqrt(scale)
    r.wall_chamber_geometry = None
    if frequency in (1.e-13,1.e13):
        old=copy.deepcopy(r.electronegative_wall_last_valid_state)
        with pytest.raises(PlasmaStateError,match='domain:.*(diffusion_length|EP wall frequency)'):
            r.monitor_electronegative_wall(y,1.)
        assert r.electronegative_wall_last_valid_state[0]==old[0]
        np.testing.assert_array_equal(r.electronegative_wall_last_valid_state[1],old[1])
        print('FREQUENCY_REFUSAL',frequency)
    else:
        actual=r.compute_electronegative_wall_jacobian(y)
        expected=decimal_wall(r,y)
        np.testing.assert_allclose(actual,expected,rtol=2.e-10,atol=0.)
        print('FREQUENCY_EDGE',frequency)


@pytest.mark.parametrize('variable',['electron','anion','neutral','Te','kz/kr','volume'])
@pytest.mark.parametrize('edge',['low','high'])
def test_outside_domain_refuses_without_publishing(variable,edge):
    r,_,_=model(energy=True)
    y=r.y.copy(); v=r.compute_volume(y)
    if variable in ('electron','anion','neutral'):
        j={'electron':0,'anion':4,'neutral':3}[variable]
        y[j]=v*(1.e-31 if edge=='low' else (1.e4 if j in (0,4) else 1.e7))
        if j == 0: y[2]=y[0]+y[4]
        if j == 4 and edge=='low': y[2]=y[0]+y[4]
        # Preserve the desired concentration by keeping the EOS pressure.
        r.P.value_si=constants.R*(r.T.value_si*sum(y[1:-1])+y[-1]*y[0])/v
        match={'electron':'e-','anion':'Cl-','neutral':'Cl'}[variable]
    elif variable=='Te':
        y[-1]=(0.04 if edge=='low' else 51.)*constants.e*constants.Na/constants.R
        match='Te'
    elif variable=='volume':
        y[:r.num_core_species]*=(1.e-13 if edge=='low' else 1.e4)/v
        match='volume'
    else:
        r.wall_diffusion_components=(1.e3,1.e3*(1.e-9 if edge=='low' else 1.e9))
        match='kz/kr'
    old=copy.deepcopy(r.electronegative_wall_last_valid_state)
    with pytest.raises(ElectronegativeWallRegimeError,match=match) as refusal:
        r.monitor_electronegative_wall(y,1.)
    print('DOMAIN_REFUSAL',variable,edge,str(refusal.value))
    assert r.electronegative_wall_last_valid_state[0]==old[0]
    np.testing.assert_array_equal(r.electronegative_wall_last_valid_state[1],old[1])


@pytest.mark.parametrize('blanc',[False,True])
@pytest.mark.parametrize('energy',[False,True])
@pytest.mark.parametrize('mixture',[False,True])
@pytest.mark.parametrize('arm',['fullFrequency','radialOnly'])
@pytest.mark.parametrize('closure',['confinedAnion','electropositiveBracket'])
@pytest.mark.parametrize('ion_temperature',[None,'gas'])
def test_domain_grid_every_wall_column(energy,mixture,arm,closure,ion_temperature,blanc):
    overrides=blanc_overrides(mixture) if blanc else {}
    if energy and mixture:
        overrides['electron_energy_balance']=dict(absorbed_power=(0.,'W'),
            chamber_volume=(np.pi*.05**2*.30,'m^3'),sheath='floating_wall',
            elastic_collisions={label:{'ignore':'test-only wall isolation'} for label in ('Ar','Cl','Ne')},
            electron_energies={'EN-test:1':(0.,'eV'),'EN-test:2':(0.,'eV')})
    r,_,_=model(energy=energy,arm=arm,closure=closure,
                ambipolar_ion_temperature=ion_temperature,second_cation=mixture,**overrides)
    # Cross both independent charge-density edges, geometry and temperature
    # edges, with interior points. Neutral extremes exercise the transport
    # floor as well as charged-pressure cancellation. Zero Cl is also covered.
    states=list(product([1.e-30,1.e-6,1.e3],[1.e-30,1.e-6,1.e3],
                        [1.e-8,1.,1.e8],[0.05,3.,50.],[1.e-30,1.,1.e6]))
    worst=0.; location=None; tested=0
    for index,(ce,cm,ratio,te,neutral) in enumerate(states):
        if not mixture and ce+cm > 1.e3:
            continue  # No charge-balanced single-cation state lies in-domain.
        tested+=1
        # A trace second bath resolves the dominant-bath mobility gradient.
        trace=max(1.e-30,neutral*1.e-30) if blanc else 0.
        concentrations=[ce,neutral,ce+cm,trace,cm]
        if mixture:
            concentrations[2]*=0.5
            concentrations.extend([0.,0.5*(ce+cm)])
        volume=[1.e-12,1.,1.e3][index%3]
        y=set_state(r,concentrations,te,volume=volume,ratio=ratio)
        expected=decimal_wall(r,y)
        actual=r.compute_electronegative_wall_jacobian(y)
        assert np.isfinite(actual).all()
        nz=expected!=0.
        errors=np.abs((actual[nz]-expected[nz])/expected[nz])
        error=float(np.max(errors,initial=0.))
        if error>worst: worst,location=error,[ce,cm,ratio,te,neutral,volume]
        np.testing.assert_allclose(actual,expected,rtol=2.e-10,atol=0.)
    print('DOMAIN_SWEEP',energy,mixture,arm,closure,ion_temperature,blanc,'cases',tested,'worst',worst,'state',location)


@pytest.mark.parametrize('energy',[False,True])
@pytest.mark.parametrize('arm',['fullFrequency','radialOnly'])
def test_trace_blanc_bath_keeps_representable_derivatives(energy,arm):
    r,_,_=model(energy=energy,arm=arm,**blanc_overrides())
    y=set_state(r,[1.e-30,1.,2.e-30,1.e-30,1.e-30])
    expected=decimal_wall(r,y)
    actual=r.compute_electronegative_wall_jacobian(y)
    print('BLANC_TRACE',energy,arm,'actual',actual[2,1],'Decimal',expected[2,1])
    np.testing.assert_allclose(actual,expected,rtol=2.e-10,atol=0.)


@pytest.mark.parametrize('te_ev',[0.05,3.,50.])
def test_analytic_electron_rate_slopes(te_ev):
    from rmgpy.kinetics import (TwoTemperaturePlasma, ElectronCollisionPlasma,
                               BadnellRRArrhenius, VoronovEIArrhenius)
    from rmgpy.solver.electronegative import electron_rate_derivative
    kinetics=[
        TwoTemperaturePlasma(A=(1.e9,'m^3/(mol*s)'),n=.5,Ea_g=(10.,'kJ/mol'),Ea_e=(50.,'kJ/mol')),
        ElectronCollisionPlasma(energies=([1.,5.,10.],'eV/molecule'),sigma=([1.e-20,2.e-20,5.e-21],'m^2')),
        BadnellRRArrhenius(A=(1.e-13,'cm^3/(molecule*s)'),B=.7,T0=(100.,'K'),T1=(1.e5,'K'),C=.3,T2=(2.e5,'K')),
        VoronovEIArrhenius(A=(3.e-8,'cm^3/(molecule*s)'),P=.25,X=1.8,K=.42,dE=13.598433)]
    te=te_ev*constants.e*constants.Na/constants.R
    for kin in kinetics:
        evaluate=(lambda v:kin.get_rate_coefficient_two_temp(298.15,v)) if hasattr(kin,'get_rate_coefficient_two_temp') else kin.get_rate_coefficient_electron_temp
        step=te*2.e-6
        expected=(evaluate(te+step)-evaluate(te-step))/(2.*step)
        assert electron_rate_derivative(kin,te,evaluate(te)) == pytest.approx(expected,rel=5.e-8,abs=0.)


@pytest.mark.parametrize('arm',['fullFrequency','radialOnly'])
@pytest.mark.parametrize('closure',['confinedAnion','electropositiveBracket'])
@pytest.mark.parametrize('quasineutral',[False,True])
def test_full_energy_operator_with_power_elastic_source_and_te_chemistry(arm,closure,quasineutral):
    from rmgpy.kinetics import TwoTemperaturePlasma
    from plasmaEnergyBalanceTest import LL_AR_ELASTIC
    r,sp,rxns=model(arm=arm,energy=True,closure=closure,quasineutral_electron=quasineutral,
                    ionisation_source=(1.e15,'m^-3*s^-1'))
    config=copy.deepcopy(r.electron_energy_balance)
    config['absorbed_power']=(1.e-4,'W')
    config['elastic_collisions']['Ar']=LL_AR_ELASTIC
    config['electron_energies']={'EN-test:1':(-.2,'eV'),'EN-test:2':(.3,'eV')}
    rxns[1].kinetics=TwoTemperaturePlasma(A=(1.e9,'m^3/(mol*s)'),n=.5,
        Ea_g=(5.,'kJ/mol'),Ea_e=(5.,'kJ/mol'),T0=(3.e4,'K'))
    FixtureEnergyReactor.entries[2]=rxns[1].kinetics
    r._configure_energy_balance(config)
    r.initialize_model(sp,rxns,[],[])
    y=r.y.copy()
    actual=r.jacobian(0.,y,np.zeros_like(y),0.)
    expected=np.zeros_like(actual)
    for column in range(len(y)):
        step=y[column]*2.e-5
        a,b=y.copy(),y.copy();a[column]+=step;b[column]-=step
        expected[:,column]=(r.residual(0.,a,np.zeros_like(y))[0]-r.residual(0.,b,np.zeros_like(y))[0])/(2.*step)
    scale=np.max(np.abs(expected*y[None,:]),axis=1)
    scale[scale==0.]=1.
    np.testing.assert_allclose(actual*y[None,:]/scale[:,None],expected*y[None,:]/scale[:,None],rtol=2.e-5,atol=2.e-7)


@pytest.mark.parametrize('arm',['fullFrequency','radialOnly'])
def test_energy_diagnostics_get_potential_without_dense_jacobian(monkeypatch,arm):
    r,_,_=model(arm=arm,energy=True)
    y=r.y.copy()
    expected=r._en_wall_linearization(y,r.compute_volume(y))[2]
    def forbidden(*args,**kwargs):
        raise AssertionError('diagnostics built a dense wall Jacobian')
    monkeypatch.setattr(FixtureEnergyReactor,'_en_wall_linearization',forbidden)
    r._latch_wall_diagnostics(y,r.compute_volume(y),1.)
    assert r.electronegative_wall_diagnostics['floating_potential_e_over_kTe'] == pytest.approx(expected,rel=1.e-14)
