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


"""Accepted-state gate census and Decimal Jacobian review witnesses."""
import copy
import math
import pytest
from pathlib import Path
from decimal import Decimal, localcontext
import numpy as np
import rmgpy.data.rmg as data
from rmgpy import constants
from rmgpy.kinetics import Arrhenius, TwoTemperaturePlasma
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.species import Species
from plasmaElectronegativeWallTest import model
from plasmaElectronegativeRework4Test import set_state
from plasmaElectronegativeRework6Test import reinitialise, decimal_operator


@pytest.fixture(autouse=True)
def isolated_database():
    previous=data.database
    data.database=None
    yield
    data.database=previous

def rx(left, right, law):
    return LibraryReaction(reactants=left, products=right, reversible=False,
        library='EN-test', kinetics=law)

def k(value, units):
    return Arrhenius(A=(value, units))

def check_domain(r, y):
    r._check_plasma_rate_parameters()
    r._check_en_transport_parameters()
    r._check_en_wall_domain(y)
    assert np.isfinite(r.kf).all() and np.all((r.kf >= 0)&(r.kf <= 1.e15))

def test_catalyst_jacobians():
    for energy in (False, True):
        r, sp, reactions = model(energy=energy)
        n = Species(label='N').from_adjacency_list('1 N u3 p1 c0')
        ncl = Species(label='NCl').from_adjacency_list(
            'multiplicity 3\n1 N u2 p1 c0 {2,S}\n2 Cl u0 p3 c0 {1,S}')
        n.thermo, ncl.thermo = copy.deepcopy(sp[3].thermo), copy.deepcopy(sp[3].thermo)
        sp.extend((n, ncl))
        r.initial_mole_fractions.update({n:1e-29, ncl:0.})
        reactions.extend((rx([sp[0],sp[2]], [sp[1]], k(2.,'m^3/(mol*s)')),
            rx([sp[1],n,sp[3]], [sp[1],ncl], k(1e15,'m^6/(mol^2*s)'))))
        reinitialise(r, sp, reactions)
        r.wall_recycling = 0.
        y = set_state(r, [.3e-6,1.,1.3e-6,1.,1e-6,1e6,0.])
        check_domain(r,y)
        actual = r.jacobian(7.,y,np.zeros_like(y),0.)
        reference = decimal_operator(r,y,full=True,precision=200)
        r.monitor_electronegative_wall(y,7.)
        assert r.electronegative_wall_last_valid_state[0] == 7.
        np.testing.assert_allclose(actual[1,[0,2]], reference[1,[0,2]], rtol=2e-10, atol=0.)
        assert reference[1,0] > 2.59e-6 and reference[1,2] > 5.99e-7
        print('P1_CATALYST_JACOBIAN', energy, 'published_t',7.,
            'J[Ar,e]',actual[1,0],reference[1,0],
            'J[Ar,Ar+]',actual[1,2],reference[1,2],flush=True)

def test_independent_flux_conditioning_is_documented():
    for energy in (False,True):
        r,sp,reactions = model(energy=energy)
        a = Species(label='NCl3').from_adjacency_list(
            'multiplicity 3\n1 N u2 p1 c0 {2,S}\n2 Cl u0 p3 c0 {1,S}')
        b = Species(label='NCl1').from_adjacency_list(
            '1 N u0 p2 c0 {2,S}\n2 Cl u0 p3 c0 {1,S}')
        a.thermo,b.thermo = copy.deepcopy(sp[3].thermo),copy.deepcopy(sp[3].thermo)
        sp.extend((a,b));r.initial_mole_fractions.update({a:1e-8,b:1e-8})
        reactions.extend((rx([a],[b],k(1e15,'s^-1')),
            rx([b],[a],k(.001,'s^-1')),
            rx([b],[a],Arrhenius(A=(1e15,'s^-1'),T0=(2.,'K')))))
        reinitialise(r,sp,reactions)
        y = set_state(r,[1e-30,1.,2e-30,.001,1e-30,1.,1.])
        check_domain(r,y)
        actual = r.residual(7.,y,np.zeros_like(y))[0]
        reference = decimal_operator(r,y,full=True,values=True,precision=200)
        # Independent amount-rate arithmetic: unary V*(k*y/V) = k*y.
        with localcontext() as ctx:
            ctx.prec=200; D=lambda x:Decimal.from_float(float(x))
            exact = -D(r.kf[-3])*D(y[5])+D(r.kf[-2])*D(y[6])+D(r.kf[-1])*D(y[6])
        r.monitor_electronegative_wall(y,7.)
        assert r.electronegative_wall_last_valid_state[0] == 7.
        assert actual[5] == actual[6] == 0.
        assert reference[5] == float(exact) == .001 and reference[6] == -.001
        assert "opposing 1e15" in (Path(__file__).resolve().parents[3]/"documentation/source/users/rmg/input.rst").read_text()
        print('P1_FLUX_ACCUMULATION',energy,'published_t',7.,
            'NCl3 mol/s',actual[5],reference[5],
            'NCl1 mol/s',actual[6],reference[6],flush=True)

def test_signed_energy_jacobian():
    r,sp,reactions = model(energy=True)
    reactions[0].kinetics.A.value_si = 0.
    reactions.extend((rx([sp[2],sp[4]],[sp[1],sp[3]],k(1e12,'m^3/(mol*s)')),
        rx([sp[0],sp[1]],[sp[2],sp[0],sp[0]],TwoTemperaturePlasma(
            A=(1e14,'m^3/(mol*s)'),n=1.,Ea_g=(0.,'J/mol'),Ea_e=(0.,'J/mol'),
            T0=(30000.,'K'),electrons=1))))
    reinitialise(r,sp,reactions)
    y=set_state(r,[1e-30,1e6,1e-6,.001,1e-6])
    r.kf[-1]=r.evaluate_two_temperature_rate_coefficient(reactions[-1].kinetics)
    r.energy_threshold[-1]=-868368.1002340637
    check_domain(r,y)
    actual=r.jacobian(7.,y,np.zeros_like(y),0.)[-1,-1]
    reference=decimal_operator(r,y,full=True,precision=200)[-1,-1]
    r.monitor_electronegative_wall(y,7.)
    assert r.electronegative_wall_last_valid_state[0] == 7.
    np.testing.assert_allclose(actual,reference,rtol=2e-10,atol=0.)
    print('P1_SIGNED_ENERGY_JACOBIAN','published_t',7.,'k',r.kf[-1],
        'J[Te,Te]',actual,reference,flush=True)
    for precision in (200,800):
        with localcontext() as ctx:
            ctx.prec=precision;D=lambda x:Decimal.from_float(float(x))
            te,ne,tg,p,gasR=map(D,(r.Te.value_si,y[0],r.T.value_si,r.P.value_si,constants.R))
            volume=gasR*(tg*sum(D(v) for v in y[1:r.num_core_species])+te*ne)/p
            theta=D(r.energy_threshold[-1])/(D(1.5)*gasR)
            coefficient=D(r.kf[-1]);dv=gasR*ne/p
            chemical=-coefficient*D(y[1])/volume*(theta/te+2)
            chemical+=coefficient*D(y[1])/volume*(theta+te)*dv/volume
            wall=r.compute_electronegative_wall_jacobian(y)[-1,-1]
            independently=float(chemical+D(wall))
            assert independently == reference
            print('INDEPENDENT_ENERGY_REFERENCE',precision,independently,flush=True)


from rmgpy import constants
from rmgpy.kinetics import TwoTemperaturePlasma
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.species import Species
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.exceptions import PlasmaStateError
from plasmaEnergyBalanceTest import _thermo
EVK=constants.e*constants.Na/constants.R
law=TwoTemperaturePlasma(A=(1e16*math.exp(20.),'m^3/(mol*s)'),n=-20.,
    Ea_g=(40.,'eV/molecule'),Ea_e=(40.,'eV/molecule'),T0=(2*EVK,'K'),electrons=1)
class Reactor(PlasmaReactor):
    def _declared_entry_kinetics(self,library,index):return law
    def check_wall_support(self,y,accepted=True):
        self.checked.append((self.t,float(y[-1])/EVK))
        return super().check_wall_support(y,accepted)
    def __init__(self,*args,**kw):
        self.checked=[];self.evaluations=[];super().__init__(*args,**kw)
    def residual(self,t,y,dydt,*args,**kw):
        result=super().residual(t,y,dydt,*args,**kw)
        self.evaluations.append((float(t),float(y[-1])/EVK,float(self.kf[0])))
        return result

def build_ep_cap_crossing():
    e=Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar=Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    ap=Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    ar.thermo,ap.thermo=_thermo(0.),_thermo(15.76)
    reaction=LibraryReaction(reactants=[e,ar],products=[ar,e],reversible=False,
        library='EN-test',kinetics=law)
    r=Reactor((298.15,'K'),(5*101325/760,'Pa'),{e:1e-6,ar:1.-2e-6,ap:1e-6},
        (.5*EVK,'K'),diffusion_length=(.02,'m'),
        ion_reduced_mobility=(1.535e-4,'m^2/(V*s)'),thermo_source_assertions={'Ar+': 'ion'},
        electron_energy_balance=dict(absorbed_power=(1000.,'W'),chamber_volume=(.002,'m^3'),
            sheath='floating_wall',elastic_collisions={'Ar':{'ignore':'isolate cap enforcement'}},
            electron_energies={'EN-test:1':(0.,'eV')}))
    r.initialize_model([e,ar,ap],[reaction],[],[])
    return r



def test_legacy_advance_checks_internal_accepted_states():
    r=build_ep_cap_crossing()
    previous=copy.deepcopy(r.energy_budget)
    flux=r.wall_flux.copy()
    tout=1.5*constants.R*r.y[0]*(6.-.5)*EVK/r.absorbed_power_reactor
    with pytest.raises(PlasmaStateError, match=r'rate coefficient.*2697707.*1e\+15') as refusal:
        r.advance(tout)
    assert r.t < tout
    assert r.energy_budget == previous
    np.testing.assert_array_equal(r.wall_flux,flux)
    print('ADVANCE_INTERNAL_REFUSAL',str(refusal.value),'time',r.t,'published',r.energy_budget['t'])


PUBLIC_STATE_APIS = (
    'electronegative_wall_factor', 'compute_mixture_reduced_mobilities',
    'compute_ion_wall_frequencies', 'compute_ion_wall_components',
    'compute_neutral_wall_frequencies', 'compute_anion_transport_data',
    'compute_anion_transport_frequencies', 'compute_reference_reaction_data',
    'compute_electronegative_wall_jacobian', 'check_wall_support',
    'monitor_electronegative_wall', 'electronegative_wall_manifest',
    'get_non_chemical_char_rate','steady_state_relaxation_time',
    'steady_state_external_armed','steady_state_external_residual','advance_same_time')


def call_public(r,name,y):
    if name=='advance_same_time':
        return r.advance(r.t)
    if name=='get_non_chemical_char_rate':
        return r.get_non_chemical_char_rate()
    if name in ('steady_state_relaxation_time','steady_state_external_armed'):
        return getattr(r,name)(7.,y)
    if name=='steady_state_external_residual':
        return r.steady_state_external_residual(7.,y,6.,y)
    if name=='electronegative_wall_manifest':
        return r.electronegative_wall_manifest()
    if name=='monitor_electronegative_wall':
        return r.monitor_electronegative_wall(y,7.)
    if name in ('compute_ion_wall_frequencies','compute_ion_wall_components',
                'compute_neutral_wall_frequencies','compute_anion_transport_data',
                'compute_anion_transport_frequencies'):
        return getattr(r,name)(y,r.compute_volume(y))
    return getattr(r,name)(y)


@pytest.mark.parametrize('name',PUBLIC_STATE_APIS)
@pytest.mark.parametrize('mutation',('mobility','eigenvalue','rate'))
def test_public_state_gate_census(name,mutation):
    r,_,_=model()
    y=r.y.copy()
    previous=copy.deepcopy(r.electronegative_wall_last_valid_state)
    diag=copy.deepcopy(r.electronegative_wall_diagnostics)
    if mutation=='mobility':
        r.wall_ion_mobility_si[2]=1e-9
        match=r'mobility.*1e-09.*1e-08'
    elif mutation=='eigenvalue':
        r.wall_diffusion_components=(1e-11,1e-6)
        r.diffusion_length.value_si=1./np.sqrt(sum(r.wall_diffusion_components))
        r.wall_chamber_geometry=None
        match=r'radial.*1e-11.*1e-10'
    else:
        r.kf[0]=1e16
        match=r'rate coefficient.*1e\+16.*1e\+15'
    with pytest.raises(PlasmaStateError,match=match) as refusal:
        call_public(r,name,y)
    assert previous[0]==r.electronegative_wall_last_valid_state[0]
    np.testing.assert_array_equal(previous[1],r.electronegative_wall_last_valid_state[1])
    assert r.electronegative_wall_diagnostics==diag
    print('PUBLIC_GATE',name,mutation,str(refusal.value),'last_valid',previous[0])


def test_gate_static_census():
    """Acceptance/publication sinks cannot lose their domain-gate call."""
    import re
    source=(Path(__file__).resolve().parents[3]/'rmgpy/solver/plasma.pyx').read_text()
    blocks=dict((m.group(1),m.group(2)) for m in re.finditer(
        r'^    (?:def|cpdef|cdef) (?:[\w.]+ )?(\w+)\([^\n]*\n(.*?)(?=^    (?:def|cpdef|cdef|@)|\Z)',
        source,re.M|re.S))
    direct=('step','advance','check_wall_support','monitor_electronegative_wall',
            '_evaluate_electronegative_wall_regime','_latch_wall_diagnostics',
            '_latch_energy_budget','electronegative_wall_manifest',
            'compute_reference_reaction_data','compute_mixture_reduced_mobilities',
            'electronegative_wall_factor','compute_electronegative_wall_jacobian',
            'evaluate_two_temperature_rate_coefficient','extinction_persistence_time',
            'get_non_chemical_char_rate','steady_state_external_residual',
            'steady_state_external_armed','steady_state_relaxation_time',
            '_record_electronegative_wall_output','_update_terminal_state',
            '_check_en_complete_operator')
    prepared=('compute_nu_wall','compute_ion_wall_components','compute_ion_wall_frequencies',
              'compute_neutral_wall_frequencies','compute_anion_transport_data')
    for name in direct:
        assert '_check_accepted_plasma_domain(' in blocks[name],name
    for name in prepared:
        assert '_prepare_public_wall_state(' in blocks[name],name
    assert {name for name in blocks if name.startswith(('compute_', 'evaluate_'))} == {
        'compute_volume', *prepared, 'compute_reference_reaction_data',
        'compute_mixture_reduced_mobilities','compute_anion_transport_frequencies',
        'compute_electronegative_wall_jacobian','evaluate_two_temperature_rate_coefficient'}
    assert '_check_accepted_plasma_domain(' in blocks['_prepare_public_wall_state']
    assert 'self.compute_anion_transport_data(' in blocks['compute_anion_transport_frequencies']
    assert 'ReactionSystem.advance(' not in blocks['advance']
    assert 'ReactionSystem.step(' not in blocks['advance']
    assert 'self._native_solver_step(' in blocks['advance']
    assert 'self.step(' in blocks['advance']
    print('GATE_CENSUS',len(direct)+len(prepared)+1,'acceptance/publication/public paths')


def test_public_rate_evaluator_checks_evaluated_value():
    r,_,_=model(energy=True)
    previous=copy.deepcopy(r.electronegative_wall_last_valid_state)
    law=TwoTemperaturePlasma(A=(1e16,'m^3/(mol*s)'),n=0.,Ea_g=(0.,'J/mol'),
        Ea_e=(0.,'J/mol'),T0=(30000.,'K'),electrons=1)
    with pytest.raises(PlasmaStateError,match=r'evaluated rate coefficient.*1e\+16.*1e\+15') as refusal:
        r.evaluate_two_temperature_rate_coefficient(law)
    assert previous[0]==r.electronegative_wall_last_valid_state[0]
    np.testing.assert_array_equal(previous[1],r.electronegative_wall_last_valid_state[1])
    print('PUBLIC_RATE_EVALUATOR',str(refusal.value))


@pytest.mark.parametrize('quantity,value,key',[
    ('mobility_T_factor',1e-5,'mobility_temperature_factor'),
    ('mobility_reference_density',1e19,'mobility_reference_density_m_3'),
    ('wall_recycling',1.1,'wall_recycling'),
    ('max_ionisation_degree',1e7,'max_ionisation_degree')])
def test_other_live_transport_parameters(quantity,value,key):
    r,_,_=model()
    previous=copy.deepcopy(r.electronegative_wall_last_valid_state)
    setattr(r,quantity,value)
    with pytest.raises(PlasmaStateError,match=key):
        r.compute_reference_reaction_data(r.y.copy())
    assert previous[0]==r.electronegative_wall_last_valid_state[0]
    np.testing.assert_array_equal(previous[1],r.electronegative_wall_last_valid_state[1])


@pytest.mark.parametrize('component',('mobility','diffusion'))
@pytest.mark.parametrize('value',(0.,1e-11,float('nan'),float('inf')))
def test_transport_cache_values_are_checked_before_blanc_arithmetic(component,value):
    from plasmaElectronegativeRework4Test import blanc_overrides
    overrides=blanc_overrides(False)
    overrides['wall_neutral_diffusion']={
        'Cl*':{'product':'Cl','diffusivity':{
            'Ar':(1e20,'1/(m*s)'), 'Cl':(1e20,'1/(m*s)')}}}
    r,sp,reactions=model(**overrides)
    metastable=Species(label='Cl*').from_adjacency_list('1 Cl u1 p3 c0')
    metastable.thermo=copy.deepcopy(sp[3].thermo)
    metastable.thermo.H298.value_si += 1000.
    sp.append(metastable)
    r.initial_mole_fractions[metastable]=0.
    r.wall_bath_lumping['Cl*']='Cl'
    reinitialise(r,sp,reactions)
    previous=copy.deepcopy(r.electronegative_wall_last_valid_state)
    if component=='mobility':
        r.wall_ion_bath_k0[2,0]=value
        match='cation bath mobility'
    else:
        r.wall_neutral_bath_dn[5,0]=value
        match=r'neutral bath D\*N'
    with pytest.raises(PlasmaStateError,match=match):
        r.compute_reference_reaction_data(r.y.copy())
    assert previous[0]==r.electronegative_wall_last_valid_state[0]
    np.testing.assert_array_equal(previous[1],r.electronegative_wall_last_valid_state[1])
