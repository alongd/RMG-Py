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

"""Public-interface tests for the declared finite-cylinder extension.

All reference callables in this module are test-only fixtures. They qualify
interface behavior, never the oxygen mechanism or a production reference.
"""
import copy
import pickle
from pathlib import Path

import numpy as np
import pytest

import rmgpy
import rmgpy.constants as constants
import rmgpy.data.rmg as data_module
from rmgpy.kinetics import Arrhenius
from rmgpy.reaction import Reaction
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.exceptions import PlasmaStateError
try:
    from rmgpy.exceptions import ElectronegativeWallRegimeError
except ImportError:
    # The base has no named regime exception; its missing input surface must
    # still make these behavior tests RED rather than stopping collection.
    ElectronegativeWallRegimeError = PlasmaStateError
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.species import Species
from rmgpy.thermo import ThermoData

RADIUS, LENGTH = 0.05, 0.30
KR, KZ = (2.405 / RADIUS)**2, (np.pi / LENGTH)**2


@pytest.fixture(autouse=True)
def tested_tree():
    assert Path.cwd() in Path(rmgpy.__file__).resolve().parents
    saved = data_module.database
    data_module.database = None
    yield
    data_module.database = saved


def fixture_reference(context):
    # Deliberately synthetic: the monitor must report its provenance.
    return dict(cation_flux=context['simple_cation_flux'].copy(), in_domain=True,
                independent=True, spans_transition=True, provenance='TEST ONLY')


def qualification(**overrides):
    value = dict(full_profile_reference=fixture_reference, wall_threshold=0.01,
                 geometry_reference=fixture_reference, geometry_threshold=0.01,
                 regime=dict(collisional=True, unmagnetised=True,
                             confined_anions=True, homogeneous_profile=True,
                             no_boundary_ionisation=True, isotropic_diffusion=True,
                             compatible_boundaries=True))
    value.update(overrides)
    return value


class FixtureEnergyReactor(PlasmaReactor):
    entries = {}

    def _declared_entry_kinetics(self, library, index):
        assert library == 'EN-test'
        return self.entries[index]


class FixtureTransportReactor(PlasmaReactor):
    # Test-only public transport seam for the mandated degenerate cases.
    prescribed = None

    def compute_anion_transport_frequencies(self,y,V):
        actual = super().compute_anion_transport_frequencies(y,V)
        if self.prescribed is None:
            return actual
        return {j:self.prescribed for j in actual}


def model(alpha=0.5, arm='fullFrequency', closure='confinedAnion',
          qualify=True, energy=False, reactions=True, **overrides):
    cls = overrides.pop('reactor_cls',PlasmaReactor)
    second = overrides.pop('second_anion',False)
    second_cation = overrides.pop('second_cation',False)
    dianion = overrides.pop('dianion',False)
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    arp = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    cl = Species(label='Cl').from_adjacency_list('1 Cl u1 p3 c0')
    clm = Species(label='Cl-').from_adjacency_list('1 O u0 p4 c-2' if dianion else '1 Cl u0 p4 c-1')
    for spc in (ar, arp, cl, clm):
        spc.thermo = ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'),
                              Cpdata=([20.8]*7,'J/(mol*K)'),
                              H298=(15.76 if spc is arp else 0.,'kJ/mol'), S298=(150.,'J/(mol*K)'))
    ne = 1.e-6
    species = [electron, ar, arp, cl, clm]
    amounts = {electron:ne, arp:ne*(1.+alpha), clm:ne*alpha,
               ar:0.999-2.*ne*(1.+alpha), cl:0.001}
    rxns = []
    if reactions:
        # Exact gross chemical frequencies, with both destruction and attachment;
        # chemistry is intentionally artificial and carries no production data.
        rxns = [Reaction(reactants=[clm],products=[cl,electron],reversible=False,
                         kinetics=Arrhenius(A=(1.e6,'s^-1'))),
                Reaction(reactants=[electron,cl],products=[clm],reversible=False,
                         kinetics=Arrhenius(A=(1.e9,'m^3/(mol*s)')))]
    if second:
        om = Species(label='O-').from_adjacency_list('1 O u1 p3 c-1')
        om.thermo = copy.deepcopy(clm.thermo)
        species.append(om)
        amounts[om] = ne*0.1
        amounts[arp] += ne*0.1
        amounts[ar] -= ne*0.2
        rxns.append(Reaction(reactants=[om],products=[cl,electron],reversible=False,
                             kinetics=Arrhenius(A=(1.,'s^-1'))))
        # Stoichiometry is deliberately not used as oxygen chemistry data.
        # Replace the neutral product with O to conserve atoms.
        o = Species(label='O').from_adjacency_list('1 O u2 p2 c0')
        o.thermo = copy.deepcopy(cl.thermo); species.append(o); amounts[o]=0.
        rxns[-1].products = [o,electron]
    if second_cation:
        ne0 = Species(label='Ne').from_adjacency_list('1 Ne u0 p4 c0')
        nep = Species(label='Ne+').from_adjacency_list('multiplicity 2\n1 Ne u1 p3 c+1')
        ne0.thermo = copy.deepcopy(ar.thermo); nep.thermo = copy.deepcopy(arp.thermo)
        species.extend([ne0,nep]);amounts[ne0]=0.;amounts[nep]=ne*0.3
        amounts[arp] -= ne*0.3
    kw = dict(diffusion_length=(1./np.sqrt(KR+KZ),'m'),
              ion_reduced_mobilities={'Ar+':(1.535e-4,'m^2/(V*s)')},
              anion_reduced_mobilities={'Cl-':(1.5e-4,'m^2/(V*s)')},
              wall_diffusion_components=(KR,KZ),
              wall_chamber_geometry=dict(shape='cylinder', radius=RADIUS, length=LENGTH),
              electronegative_wall_model=closure, electronegative_wall_geometry=arm,
              electronegative_wall_qualification=qualification() if qualify else None,
              wall_single_bath_approximation=True, thermo_source_assertions=['Ar+','Cl-'])
    if second_cation:
        kw['ion_reduced_mobilities']['Ne+'] = (3.07e-4,'m^2/(V*s)')
        kw['thermo_source_assertions'].append('Ne+')
    if second:
        kw['thermo_source_assertions'].append('O-')
        kw['anion_reduced_mobilities']['O-'] = (1.5e-4,'m^2/(V*s)')
    if energy:
        cls = FixtureEnergyReactor
        rxns = [LibraryReaction(reactants=r.reactants,products=r.products,
                kinetics=r.kinetics,reversible=False,library='EN-test') for r in rxns]
        cls.entries = {i+1:r.kinetics for i,r in enumerate(rxns)}
        kw['electron_energy_balance'] = dict(absorbed_power=(0.,'W'),
            chamber_volume=(np.pi*RADIUS**2*LENGTH,'m^3'),sheath='floating_wall',
            elastic_collisions={'Ar':{'ignore':'test-only wall isolation'},
                                'Cl':{'ignore':'test-only wall isolation'}},
            electron_energies={'EN-test:'+str(i+1):(0.,'eV') for i in range(len(rxns))})
    kw.update(overrides)
    r = cls((298.15,'K'),(5.*101325./760.,'Pa'),amounts,
                      (3./8.617333262e-5,'K'),**kw)
    r.initialize_model(species,rxns,[],[])
    return r,species,rxns


@pytest.mark.parametrize('arm',['fullFrequency','radialOnly'])
@pytest.mark.parametrize('alpha,expected_h',[(0.5,2./3.),(0.6,0.625)])
def test_dynamic_cation_and_charge_balanced_wall_sources(arm,alpha,expected_h):
    r,spcs,_ = model(alpha=alpha,arm=arm)
    y = r.y.copy(); v = r.compute_volume(y)
    components = r.compute_ion_wall_components(y,v)
    ep = components['radial_ep'][2]+components['axial_ep'][2]
    expected = expected_h*ep if arm=='fullFrequency' else expected_h*components['radial_ep'][2]+components['axial_ep'][2]
    nu = r.compute_ion_wall_frequencies(y,v)
    assert nu[2] == pytest.approx(expected,rel=1e-14)
    assert r.electronegative_wall_diagnostics['h'] == expected_h
    assert r.wall_flux[4] == 0.
    assert r.wall_flux[0] == r.wall_flux[2]
    assert sum(s.get_net_charge()*r.wall_flux[i] for i,s in enumerate(spcs)) == 0.
    y[4] *= 2.
    assert r.compute_ion_wall_frequencies(y,r.compute_volume(y))[2] < nu[2]


def test_unqualified_finite_anions_refuse_instead_of_receiving_a_default():
    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError,match='regime|reference|frozen'):
        model(qualify=False)


@pytest.mark.parametrize('arm',['fullFrequency','radialOnly'])
def test_analytic_jacobian_at_half_electronegativity(arm):
    # No chemistry: isolates the actual wall Jacobian without cancellation by
    # a much larger chemical destruction derivative.
    r,_,_ = model(arm=arm,closure='electropositiveBracket')
    r.kf[:] = 0.
    r.electronegative_wall_model = 'confinedAnion'
    y = r.y.copy(); dydt = np.zeros_like(y)
    analytic = np.asarray(r.jacobian(0.,y,dydt,0.))
    numerical = np.zeros_like(analytic)
    for j in range(len(y)):
        delta = abs(y[j])*2.e-5
        plus,minus = y.copy(),y.copy(); plus[j]+=delta; minus[j]-=delta
        numerical[:,j] = (np.asarray(r.residual(0.,plus,dydt)[0])-
                          np.asarray(r.residual(0.,minus,dydt)[0]))/(2.*delta)
    np.testing.assert_allclose(analytic,numerical,rtol=1.e-6,atol=2.e-7)


@pytest.mark.parametrize('arm',['fullFrequency','radialOnly'])
def test_pickle_reinitialization_preserves_closure_and_geometry(arm):
    r,_,_ = model(arm=arm)
    rebuilt = pickle.loads(pickle.dumps(r))
    assert rebuilt.electronegative_wall_model == 'confinedAnion'
    assert rebuilt.electronegative_wall_geometry == arm
    assert rebuilt.wall_diffusion_components == (KR,KZ)
    spcs=list(rebuilt.initial_mole_fractions)
    # Reinitialization recalculates monitoring at the restart physical state;
    # missing chemistry must refuse Gate A instead of bypassing qualification.
    from rmgpy.exceptions import PlasmaStateError
    with pytest.raises(PlasmaStateError,match='A:'):
        rebuilt.initialize_model(spcs,[],[],[])


@pytest.mark.parametrize('arm',['fullFrequency','radialOnly'])
def test_input_writer_preserves_named_arm_decomposition_and_monitor(arm):
    from rmgpy.rmg.input import _format_plasma_wall
    r,_,_ = model(arm=arm)
    text = ''.join(_format_plasma_wall(r))
    assert 'electronegativeWallModel' in text
    assert 'electronegativeWallGeometry' in text
    assert 'wallDiffusionComponents' in text
    assert 'electronegativeWallQualification' in text
    # Parse the writer's actual kwargs: malformed reprs of callables fail here.
    parsed = eval('dict('+text+')')
    assert parsed['electronegativeWallGeometry'] == arm
    assert parsed['electronegativeWallQualification']['full_profile_reference'].endswith(':fixture_reference')


@pytest.mark.parametrize('arm',['fullFrequency','radialOnly'])
def test_bracket_has_unit_factor_and_zero_anion_flux(arm):
    r,_,_ = model(arm=arm,closure='electropositiveBracket')
    assert r.electronegative_wall_factor(r.y) == (1.,1.)
    assert r.wall_flux[4] == 0.
    assert r.wall_flux[0] == r.wall_flux[2]


@pytest.mark.parametrize('multiplier',[0.,0.999,1.])
def test_gate_a_rejects_at_and_below_per_anion_crossover(multiplier):
    r,_,_ = model()
    transport = r.compute_anion_transport_frequencies(r.y,r.compute_volume(r.y))[4]
    r.kf[0] = multiplier*transport
    with pytest.raises(ElectronegativeWallRegimeError,match='A:'):
        r.monitor_electronegative_wall(r.y,0.)


@pytest.mark.parametrize('multiplier',[0.,1.])
def test_gate_b_is_independent_of_large_confinement(multiplier):
    r,_,_ = model()
    ep = r.compute_ion_wall_components(r.y,r.compute_volume(r.y))['total_ep'][2]
    conc_cl = r.y[3]/r.compute_volume(r.y)
    r.kf[1] = multiplier*ep/conc_cl
    with pytest.raises(ElectronegativeWallRegimeError,match='B:'):
        r.monitor_electronegative_wall(r.y,0.)


@pytest.mark.parametrize('key',['full_profile_reference','geometry_reference','wall_threshold','geometry_threshold'])
def test_missing_reference_or_frozen_threshold_refuses(key):
    with pytest.raises(ElectronegativeWallRegimeError,match='C-'):
        model(electronegative_wall_qualification=qualification(**{key:None}))


def biased_reference(context):
    result = fixture_reference(context)
    result['cation_flux'] *= 2.
    return result


def nonfinite_reference(context):
    result = fixture_reference(context)
    result['cation_flux'][:] = np.nan
    return result


def outside_reference(context):
    result = fixture_reference(context)
    result['in_domain'] = False
    return result


def circular_reference(context):
    result = fixture_reference(context)
    result['independent'] = False
    return result


def incomplete_reference(context):
    result = fixture_reference(context)
    result['spans_transition'] = False
    return result


@pytest.mark.parametrize('reference',[biased_reference,nonfinite_reference,outside_reference,circular_reference,incomplete_reference])
def test_full_profile_reference_error_domain_and_independence_are_hard_gates(reference):
    with pytest.raises(ElectronegativeWallRegimeError,match='C-full-profile'):
        model(electronegative_wall_qualification=qualification(full_profile_reference=reference))


def test_low_confinement_anion_cannot_be_averaged_away():
    with pytest.raises(ElectronegativeWallRegimeError,match='A: O-'):
        model(second_anion=True)


@pytest.mark.parametrize('which',[0,1])
def test_nonfinite_chemical_or_attachment_frequency_refuses(which):
    r,_,_ = model()
    r.kf[which] = np.nan
    with pytest.raises(ElectronegativeWallRegimeError,match='A:|B:'):
        r.monitor_electronegative_wall(r.y,0.)


def test_zero_transport_positive_destruction_is_infinite_confinement():
    FixtureTransportReactor.prescribed = 0.
    try:
        r,_,_ = model(reactor_cls=FixtureTransportReactor)
        assert r.electronegative_wall_diagnostics['gates']['A']['Cl-']['conf'] == float('inf')
        r.kf[0] = 0.
        with pytest.raises(ElectronegativeWallRegimeError,match='both destruction and transport zero'):
            r.monitor_electronegative_wall(r.y,0.)
    finally:
        FixtureTransportReactor.prescribed = None


@pytest.mark.parametrize('transport',[np.nan,np.inf,-1.])
def test_nonfinite_or_negative_transport_refuses(transport):
    FixtureTransportReactor.prescribed = transport
    try:
        with pytest.raises(ElectronegativeWallRegimeError,match='A:'):
            model(reactor_cls=FixtureTransportReactor)
    finally:
        FixtureTransportReactor.prescribed = None


def fail_after_initial_state(context):
    result = fixture_reference(context)
    result['in_domain'] = context['time'] == 0.
    return result


def test_first_accepted_violating_state_terminates_and_preserves_last_valid_state():
    r,_,_ = model(electronegative_wall_qualification=qualification(full_profile_reference=fail_after_initial_state))
    saved = r.y.copy(); flux = r.wall_flux.copy()
    # Reference is a monitor only; Newton residuals at positive time do not call it.
    r.residual(1.,saved,np.zeros_like(saved))
    with pytest.raises(ElectronegativeWallRegimeError,match='C-full-profile'):
        r.advance(1.e-8)
    valid_t,valid_y = r.electronegative_wall_last_valid_state
    assert valid_t == 0.
    np.testing.assert_array_equal(valid_y,saved)
    np.testing.assert_array_equal(r.wall_flux,flux)
    assert r.wall_diagnostics_time == 0.


def test_thresholds_are_frozen_and_diagnostics_have_every_gate():
    r,_,_ = model()
    external = r.electronegative_wall_qualification
    external['wall_threshold'] = 2.
    assert r.electronegative_wall_qualification['wall_threshold'] == 0.01
    manifest = r.electronegative_wall_manifest()
    assert manifest['closure'] == 'confinedAnion'
    assert manifest['geometry_arm'] == 'fullFrequency'
    assert manifest['scientific_status'] == 'FINITE-CYLINDER GEOMETRIC EXTENSION'
    assert set(manifest['gates']) == {'A','B','C_radial','C_geometry','C_full_profile'}
    assert manifest['gates']['A']['Cl-']['destruction_channels']
    assert manifest['f_z'] == pytest.approx(0.04525379687210524)


@pytest.mark.parametrize('arm',['fullFrequency','radialOnly'])
def test_no_anion_state_bypasses_all_qualification(arm):
    r,_,_ = model(alpha=0.,arm=arm,qualify=False,reactions=False)
    assert r.electronegative_wall_factor(r.y) == (1.,1.)
    assert r.electronegative_wall_diagnostics['gates'] == 'not-evaluated-zero-anion'


def test_no_implicit_closure_and_no_scalar_or_missing_mobility():
    with pytest.raises(PlasmaStateError,match='explicitly selected'):
        model(closure=None)
    with pytest.raises(PlasmaStateError,match='map-mode'):
        model(ion_reduced_mobilities=None,ion_reduced_mobility=(1.535e-4,'m^2/(V*s)'))
    with pytest.raises(PlasmaStateError,match='no declared mobility'):
        model(anion_reduced_mobilities={'other':(1.e-4,'m^2/(V*s)')})


@pytest.mark.parametrize('arm',['fullFrequency','radialOnly'])
@pytest.mark.parametrize('gas_temperature',[None,'gas'])
def test_energy_mode_analytic_wall_and_full_jacobian_at_half_alpha(arm,gas_temperature):
    r,_,_ = model(arm=arm,energy=True,ambipolar_ion_temperature=gas_temperature,
                    closure='electropositiveBracket')
    r.kf[:] = 0.
    r.electronegative_wall_model = 'confinedAnion'
    # Disable the synthetic chemistry after valid initialization to isolate
    # the actual wall residual and its analytic Te/phi contribution.
    y=r.y.copy(); dydt=np.zeros_like(y)
    analytic=np.asarray(r.jacobian(0.,y,dydt,0.))
    numerical=np.zeros_like(analytic)
    for j in range(len(y)):
        # A larger neutral step resolves the small EOS derivative of the
        # energy source beside its 1e7 K/s value, without relaxing accuracy.
        delta=abs(y[j])*(1.e-3 if j in (1,3) else 2.e-5)
        a,b=y.copy(),y.copy(); a[j]+=delta;b[j]-=delta
        numerical[:,j]=(np.asarray(r.residual(0.,a,dydt)[0])-np.asarray(r.residual(0.,b,dydt)[0]))/(2.*delta)
    np.testing.assert_allclose(analytic,numerical,rtol=2.e-5,atol=2.e-3)
    r.residual(0.,y,dydt)
    assert np.all(np.isfinite(r.compute_electronegative_wall_jacobian(y)))


@pytest.mark.parametrize('arm',['fullFrequency','radialOnly'])
def test_multiple_cations_keep_individual_mobility_and_charge_balance(arm):
    r,species,_=model(arm=arm,second_cation=True)
    frequencies=r.compute_ion_wall_frequencies(r.y,r.compute_volume(r.y))
    assert frequencies[-1]/frequencies[2] == pytest.approx(2.,rel=1.e-14)
    assert r.wall_flux[0] == r.wall_flux[2]+r.wall_flux[-1]
    assert np.dot(r.species_charges,r.wall_flux) == pytest.approx(0.,abs=1.e-18)


def test_multiply_charged_anion_refuses_before_transport():
    with pytest.raises(PlasmaStateError,match='multiply charged anions'):
        model(dianion=True)


def conditional_bias(context):
    return fixture_reference(context) if context['time']==0. else biased_reference(context)


def test_just_above_crossover_does_not_hide_large_reference_error():
    r,_,_=model(electronegative_wall_qualification=qualification(full_profile_reference=conditional_bias))
    transport=r.compute_anion_transport_frequencies(r.y,r.compute_volume(r.y))[4]
    r.kf[0]=transport*1.001
    with pytest.raises(ElectronegativeWallRegimeError,match='C-full-profile') as error:
        r.monitor_electronegative_wall(r.y,1.)
    assert error.value.gate == 'C-full-profile'
    assert error.value.time == 1.
    np.testing.assert_array_equal(error.value.state,r.y)


def alpha_dependent_reference(context):
    result=fixture_reference(context)
    result['cation_flux'] *= 1.+context['alpha']
    return result


def test_equal_confinement_does_not_fix_reference_error():
    FixtureTransportReactor.prescribed=100.
    try:
        records=[]
        for alpha in (0.5,0.6):
            r,_,_=model(alpha=alpha,reactor_cls=FixtureTransportReactor,
                electronegative_wall_qualification=qualification(full_profile_reference=alpha_dependent_reference,wall_threshold=1.))
            records.append(r.electronegative_wall_diagnostics['gates'])
        assert records[0]['A']['Cl-']['conf'] == records[1]['A']['Cl-']['conf']
        assert records[0]['C_full_profile']['error'] != records[1]['C_full_profile']['error']
    finally:
        FixtureTransportReactor.prescribed=None


def test_nonfinite_ratio_from_finite_frequencies_refuses():
    FixtureTransportReactor.prescribed=1.e-309
    try:
        with pytest.raises(ElectronegativeWallRegimeError,match='nonfinite confinement'):
            model(reactor_cls=FixtureTransportReactor)
    finally:
        FixtureTransportReactor.prescribed=None


def test_synthetic_envelope_reports_gate_extrema_locations_and_margins():
    from rmgpy.solver.electronegative import qualify_envelope
    a,_,_=model(alpha=0.5);b,_,_=model(alpha=0.6)
    result=qualify_envelope([(a,a.y,0.),(b,b.y,0.)])
    assert result['samples'] == 2
    assert set(result['extrema']) == {'min_conf','min_attachment_metric','max_full_profile_error','max_geometry_error'}
    for record in result['extrema'].values():
        assert record['margin'] > 0.
        assert record['sample'] in (1,2)
        assert record['state']


def test_actual_input_reference_restart_and_manifest_json():
    import json
    from rmgpy.solver.electronegative import serializable_qualification
    r,_,_=model(electronegative_wall_qualification=serializable_qualification(qualification()))
    assert r.electronegative_wall_diagnostics['gates']['C_full_profile']['error']==0.
    json.dumps(r.electronegative_wall_manifest(),allow_nan=False)


@pytest.mark.parametrize('arm',['fullFrequency','radialOnly'])
def test_floating_potential_rebalances_reduced_total_flux(arm):
    r,species,_=model(arm=arm,energy=True)
    y=r.y.copy();ne=y[0];plus=y[2]
    phi=r.electronegative_wall_manifest()['floating_potential_e_over_kTe']
    # Independent Maxwellian/Bohm current normalisation, with the explicit
    # electron/cation density ratio and actual total collection frequency.
    mass=species[2].molecular_weight.value_si
    phi_ep=np.log(np.sqrt(mass/(2.*np.pi*constants.m_e)))
    ep=r.compute_ion_wall_components(y,r.compute_volume(y))['total_ep'][2]*plus
    actual=-r.wall_flux[0]
    expected=phi_ep+np.log(ne/plus)+np.log(ep/actual)
    assert phi == pytest.approx(expected,rel=1.e-14)
    assert actual*np.exp(phi) == pytest.approx(ep*np.exp(phi_ep)*ne/plus,rel=1.e-14)


@pytest.mark.parametrize('which,gate',[(0,'A:'),(1,'B:')])
def test_comparator_keeps_confined_anion_physical_gates(which,gate):
    r,_,_=model(closure='electropositiveBracket')
    r.kf[which]=0.
    with pytest.raises(ElectronegativeWallRegimeError,match=gate):
        r.monitor_electronegative_wall(r.y,0.)


@pytest.mark.parametrize('closure', ['confinedAnion', 'electropositiveBracket'])
@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
def test_both_models_refuse_at_first_accepted_violation(closure, arm):
    class AcceptedMonitorReactor(PlasmaReactor):
        def monitor_electronegative_wall(self, state, time):
            if time > 0.:
                self.accepted_check_times.append(time)
                # Test-only change at the first accepted state makes the real
                # gross-destruction Gate A fail; Newton trials cannot do this.
                self.kf[0] = 0.
            return super().monitor_electronegative_wall(state, time)

    r, _, _ = model(closure=closure, arm=arm, reactor_cls=AcceptedMonitorReactor)
    r.accepted_check_times = []
    initial = r.y.copy()
    with pytest.raises(ElectronegativeWallRegimeError, match='A:'):
        r.advance(1.e-8)
    assert len(r.accepted_check_times) == 1
    assert 0. < r.accepted_check_times[0] < 1.e-8
    valid_time, valid_state = r.electronegative_wall_last_valid_state
    assert valid_time == 0.
    np.testing.assert_array_equal(valid_state, initial)


def test_wall_discrepancy_cannot_cancel_cation_errors_and_weights_charge():
    from rmgpy.solver.electronegative import wall_discrepancy
    # Opposing errors cancel in a net sum, but actual source error cannot.
    # Nonuniform charges independently pin the charge weighting.
    assert wall_discrepancy([5., 2.], [1., 4.], [1., 2.]) == pytest.approx(8./9.)
    assert wall_discrepancy([0., 0.], [0., 0.], [1., 2.]) == 0.


def test_reference_accepts_exact_frozen_threshold_and_refuses_above():
    from rmgpy.solver.electronegative import check_reference, wall_discrepancy
    context = dict(simple_cation_flux=np.array([1.]), cation_charges=np.array([1.]))
    def reference(context):
        return dict(cation_flux=np.array([0.5]), in_domain=True,
                    independent=True, provenance='TEST ONLY threshold fixture')
    threshold = wall_discrepancy([1.], [0.5], [1.])
    assert check_reference('C-test', reference, threshold, context)['error'] == threshold
    with pytest.raises(ElectronegativeWallRegimeError, match='exceeds'):
        check_reference('C-test', reference, np.nextafter(threshold, 0.), context)



def test_canceling_invalid_anion_samples_cannot_bypass_zero_branch():
    from rmgpy.solver.electronegative import qualify_envelope
    FixtureTransportReactor.prescribed = 0.1
    try:
        r, _, _ = model(second_anion=True, reactor_cls=FixtureTransportReactor)
        saved = r.electronegative_wall_last_valid_state[1].copy()
        invalid = r.y.copy()
        invalid[4], invalid[5] = 1.e-7, -1.e-7
        with pytest.raises(ElectronegativeWallRegimeError, match='invalid anion|domain.*O-'):
            qualify_envelope([(r, invalid, 1.)])
        np.testing.assert_array_equal(r.electronegative_wall_last_valid_state[1], saved)
    finally:
        FixtureTransportReactor.prescribed = None


def test_supplied_energy_sample_synchronizes_temperature_before_monitoring():
    r, _, _ = model(energy=True)
    initial_metric = r.electronegative_wall_diagnostics['gates']['B']['metrics'][0]
    initial_volume = r.compute_volume(r.y)
    state = r.y.copy()
    state[r.te_index] *= 2.
    r.monitor_electronegative_wall(state, 1.)
    assert r.Te.value_si == state[r.te_index]
    # Attachment is proportional to 1/V, and ambipolar diffusion to Te*V.
    expected = initial_metric / 2. * (initial_volume / r.compute_volume(state))**2
    assert r.electronegative_wall_diagnostics['gates']['B']['metrics'][0] == pytest.approx(expected)


def test_varied_envelope_pins_every_extremum_and_location():
    from rmgpy.solver.electronegative import qualify_envelope
    q = qualification(full_profile_reference=alpha_dependent_reference,
                      geometry_reference=alpha_dependent_reference,
                      wall_threshold=1., geometry_threshold=1.)
    a, _, _ = model(alpha=0.5, electronegative_wall_qualification=q)
    b, _, _ = model(alpha=0.6, electronegative_wall_qualification=q)
    b.kf[0] *= 2.
    result = qualify_envelope([(a, a.y, 1.), (b, b.y, 2.)])
    gate_sets = [r.electronegative_wall_diagnostics['gates'] for r in (a, b)]
    values = dict(
        min_conf=[g['A']['Cl-']['conf'] for g in gate_sets],
        min_attachment_metric=[g['B']['metrics'][0] for g in gate_sets],
        max_full_profile_error=[g['C_full_profile']['error'] for g in gate_sets],
        max_geometry_error=[g['C_geometry']['error'] for g in gate_sets])
    for key, samples in values.items():
        expected = min(samples) if key.startswith('min') else max(samples)
        index = samples.index(expected)
        record = result['extrema'][key]
        assert record['value'] == expected
        assert record['sample'] == index + 1
        assert record['time'] == index + 1.
        assert record['margin'] == (expected - 1. if key.startswith('min') else 1. - expected)


@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
def test_zero_electron_newton_trial_keeps_right_jacobian_coupling(arm):
    r, _, _ = model(arm=arm)
    r.kf[:] = 0.
    state = r.y.copy()
    state[r.electron_index] = 0.
    step = r.y[r.electron_index] * 1.e-6
    shifted = state.copy()
    shifted[r.electron_index] = step
    analytic = np.asarray(r.jacobian(0., state, np.zeros_like(state), 0.))[:, r.electron_index]
    numerical = (np.asarray(r.residual(0., shifted, np.zeros_like(state))[0]) -
                 np.asarray(r.residual(0., state, np.zeros_like(state))[0])) / step
    np.testing.assert_allclose(analytic, numerical, rtol=3.e-6, atol=2.e-7)


@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
def test_energy_newton_trial_below_temperature_floor_has_zero_te_derivative(arm):
    r, _, _ = model(arm=arm, energy=True)
    r.kf[:] = 0.
    state = r.y.copy()
    state[r.te_index] = 0.25 * r.T.value_si
    analytic = np.asarray(r.jacobian(0., state, np.zeros_like(state), 0.))[:, r.te_index]
    np.testing.assert_array_equal(analytic, np.zeros_like(analytic))


def test_reference_context_exposes_declared_geometry_transport_and_accepted_rates():
    records = []
    def adapter(context):
        records.append(context)
        result = fixture_reference(context)
        result.update(alpha_ref=context['alpha'], eigenvalues={'s': 1., 'a': 1.},
                      convergence={'converged': True, 'iterations': 2})
        return result
    r, _, _ = model(electronegative_wall_qualification=qualification(full_profile_reference=adapter))
    # A reference must never consume the last Newton scratch rates.
    r.core_reaction_rates[:] = np.nan
    r.monitor_electronegative_wall(r.y, 1.)
    context = records[-1]
    assert context['schema_version'] == 1
    assert context['chamber_radius'] == RADIUS
    assert context['chamber_length'] == LENGTH
    assert context['geometry_provenance'] == 'declared-cylinder'
    assert context['anion_indices'] == (4,)
    transport = context['anion_transport'][4]
    assert transport['reduced_mobility'] == 1.5e-4
    assert transport['mobility'] == pytest.approx(1.5e-4 * r.mobility_reference_density / context['neutral_number_density'])
    assert context['cation_transport'][2]['reduced_mobility'] == 1.535e-4
    rows = context['reactions']
    assert [row['classification'] for row in rows] == ['anion_neutral', 'electron_attachment']
    for row in rows:
        expected = row['coefficient_si'] * np.prod(context['concentrations_mol_per_volume'][list(row['reactants'])])
        assert row['rate_mol_per_volume_s'] == expected
        assert row['rate_mol_per_s'] == context['volume'] * expected
        assert row['reference_domain_reason'] is None
    diagnostic = r.electronegative_wall_manifest()['gates']['C_full_profile']['reference_diagnostics']
    assert diagnostic['alpha_ref'] == pytest.approx(0.5)
    assert diagnostic['eigenvalues'] == {'s': 1., 'a': 1.}
    assert diagnostic['convergence']['converged'] is True


@pytest.mark.parametrize('left,right,expected', [
    ((0, 1), (0, 0, 2), 'electron_ionisation'),
    ((0, 1), (3,), 'electron_attachment'),
    ((3, 1), (0, 1), 'anion_neutral'),
    ((3, 0), (0, 0, 1), 'anion_electron'),
    ((2, 0), (1,), 'cation_electron'),
    ((2, 3), (1, 1), 'cation_anion'),
    ((2, 1), (2, 1), 'cation_neutral'),
    ((0, 1), (0, 1), 'neutral'),
    ((3, 3), (1, 1), 'unsupported'),
    ((4, 0), (1,), 'unsupported'),
])
def test_reference_reaction_classes_use_actual_charged_reactants(left, right, expected):
    from rmgpy.solver.electronegative import classify_charged_reaction
    classification, reason = classify_charged_reaction(left, right, [-1, 0, 1, -1, 2], 0)
    assert classification == expected
    assert (reason is not None) == (expected == 'unsupported')
    assert classify_charged_reaction(left, right, [-1, 0, 1, -1, 2], 0, True) == ('unsupported', 'charged third body')


def test_eigenvalues_alone_do_not_claim_a_declared_chamber_radius():
    records = []
    def adapter(context):
        records.append(context)
        return fixture_reference(context)
    model(wall_chamber_geometry=None,
          electronegative_wall_qualification=qualification(full_profile_reference=adapter))
    assert records[0]['chamber_radius'] is None
    assert records[0]['geometry_provenance'] == 'eigenvalues-only'
    assert records[0]['cylinder_equivalent_radius'] == RADIUS


def test_declared_cylinder_must_match_the_transport_eigenvalues():
    with pytest.raises(PlasmaStateError, match='dimensions must match'):
        model(wall_chamber_geometry=dict(shape='cylinder', radius=2.*RADIUS, length=LENGTH))



def test_reference_cannot_redefine_the_production_vector_during_comparison():
    from rmgpy.solver.electronegative import check_reference
    context = dict(simple_cation_flux=np.array([1.]), cation_charges=np.array([1.]))
    def mutating_reference(context):
        context['simple_cation_flux'][:] = 0.5
        return dict(cation_flux=np.array([0.5]), in_domain=True,
                    independent=True, provenance='TEST ONLY mutating fixture')
    with pytest.raises(ElectronegativeWallRegimeError, match='exceeds'):
        check_reference('C-test', mutating_reference, 0.01, context)
    np.testing.assert_array_equal(context['simple_cation_flux'], [1.])



@pytest.mark.parametrize('arm', ['fullFrequency', 'radialOnly'])
@pytest.mark.parametrize('alpha', [1.e-1, 1.e-3, 1.e-7, 1.e-12])
def test_vanishing_anion_population_approaches_electropositive_path(alpha, arm):
    r, _, _ = model(alpha=alpha, arm=arm)
    components = r.compute_ion_wall_components(r.y, r.compute_volume(r.y))
    active = r.compute_ion_wall_frequencies(r.y, r.compute_volume(r.y))[2]
    weight = 1. if arm == 'fullFrequency' else KR / (KR + KZ)
    expected_difference = weight * alpha / (1. + alpha)
    assert 1. - active / components['total_ep'][2] == pytest.approx(expected_difference, abs=1.e-15)
    assert r.wall_flux[4] == 0.
    assert r.wall_flux[0] == r.wall_flux[2]


def test_zero_electrons_with_finite_anions_refuse_without_an_alpha_floor():
    r, _, _ = model()
    state = r.y.copy()
    state[0] = 0.
    state[2] = state[4]  # Charge-consistent finite ions, no electrons.
    with pytest.raises(ElectronegativeWallRegimeError, match='electron density undefined'):
        r.monitor_electronegative_wall(state, 1.)


def test_no_charged_particles_have_zero_charged_wall_sources():
    r, _, _ = model()
    r.kf[:] = 0.
    state = r.y.copy()
    state[[0, 2, 4]] = 0.
    r.monitor_electronegative_wall(state, 1.)
    assert r.electronegative_wall_diagnostics['gates'] == 'not-evaluated-zero-anion'
    residual = np.asarray(r.residual(1., state, np.zeros_like(state))[0])
    np.testing.assert_array_equal(residual[[0, 2, 4]], [0., 0., 0.])
