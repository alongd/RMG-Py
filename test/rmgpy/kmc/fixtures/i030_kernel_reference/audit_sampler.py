"""Additional independent lattice and continuum checks of the live oracle."""
import importlib.util
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import erfc

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('rouse_oracle_audit', ROOT/'reference/run.py')
oracle = importlib.util.module_from_spec(spec)
spec.loader.exec_module(oracle)


def direct_laplacian_checks(cfg):
    output = []
    for case in cfg['cases']:
        mapped = oracle.mapped_parameters(cfg, case)
        for pair in cfg['pair_classes']:
            covariance, local_d = 0., mapped['D_relative_reduced']
            kinds = {'end/end':('end','end'), 'end/mid':('end','mid'), 'mid/mid':('mid','mid')}[pair]
            for key, kind in zip(('units_i','units_j'), kinds):
                n = mapped[key]
                if n == 1:
                    continue
                lap = np.zeros((n,n))
                for bond in range(n-1):
                    lap[bond,bond] += 1
                    lap[bond+1,bond+1] += 1
                    lap[bond,bond+1] -= 1
                    lap[bond+1,bond] -= 1
                eigen, vec = np.linalg.eigh(lap)
                b2 = (cfg['chains'][str(n)]['bond_rms_m']/mapped['scale_m'])**2
                site = 0 if kind == 'end' else n//2-1
                covariance += float(np.sum(vec[site,1:]**2*b2/(3*eigen[1:])))
                local_d += cfg['D0_m2_s']/mapped['diffusion_unit_m2_s']*float(np.sum(vec[site,1:]**2))
                assert math.isclose(float(np.sum(b2/eigen[1:])/n), (cfg['chains'][str(n)]['Rg_m']/mapped['scale_m'])**2,rel_tol=1e-12)
            rates, weights, var = oracle.mode_coefficients(cfg,mapped,pair)
            assert math.isclose(covariance, float(np.sum(weights**2*var)), abs_tol=1e-12, rel_tol=1e-12)
            assert math.isclose(local_d,mapped['D_local_reduced'],rel_tol=1e-12)
            output.append({'case':case['name'],'pair_class':pair,'independent_covariance':covariance,
                           'local_diffusion_reduced':local_d,'pass':True})
    return output


def continuum_sphere_check(cfg):
    """Infinite Brownian absorbing sphere, with a bounded cube truncation error.

    P(hit by t | r>a)=a/r*erfc((r-a)/sqrt(4Dt)). Its integral over the
    outside of the sink is 4*pi*a*D*t+8*sqrt(pi)*a^2*sqrt(D*t).
    Replacing the cube by R3 can only overcount the r>=L/2 tail.
    """
    settings = dict(cfg['bead_crosscheck'], N=1, replicas=65536,
                    steps=800, dt_reduced=.0001, box_over_scale=6., seed=209001)
    own = dict(cfg, bead_crosscheck=settings)
    contacts = oracle.mode_grid_contact(own,'end/end',settings['seed'])
    a, diffusion, side, elapsed = 1.,2.,6.,settings['steps']*settings['dt_reduced']
    volume = side**3-4*math.pi*a**3/3
    infinite = (4*math.pi*a*diffusion*elapsed+8*math.sqrt(math.pi)*a*a*math.sqrt(diffusion*elapsed))/volume
    tail = quad(lambda r:4*math.pi*a*r*erfc((r-a)/math.sqrt(4*diffusion*elapsed)),side/2,np.inf,epsabs=1e-13)[0]/volume
    probability = contacts/settings['replicas']
    se = math.sqrt(probability*(1-probability)/settings['replicas'])
    # This is a continuum diagnostic, with statistics and cube-tail only.
    passed = infinite-tail-5*se <= probability <= infinite+5*se
    assert passed, 'exact Brownian continuum contact cross-check failed'
    return {'replicas':settings['replicas'],'contacts':contacts,'probability':probability,'se':se,
            'analytic_infinite_probability':infinite,'cube_tail_upper_bound':tail,
            'dt_reduced':settings['dt_reduced'],'elapsed_reduced':elapsed,'pass':passed}


def stationary_variance_checks(cfg):
    output = []
    for ci, case in enumerate(cfg['cases']):
        mapped = oracle.mapped_parameters(cfg, case)
        for pi, pair in enumerate(cfg['pair_classes']):
            rates, weights, variances = oracle.mode_coefficients(cfg, mapped, pair)
            if not len(rates):
                continue
            rng = np.random.Generator(np.random.PCG64(210001+100*ci+pi))
            q = rng.normal(size=(4096,len(rates),3))*np.sqrt(variances)[None,:,None]
            dt = np.full(4096,.001)
            theory = float(np.sum(weights**2*variances))
            for step in range(201):
                if step in (0,1,10,100,200):
                    coordinates = oracle.projected_site(q,weights)[:,0]
                    observed, se = oracle.mean_se(coordinates**2)
                    passed = abs(observed-theory) <= 6*se
                    output.append({'case':case['name'],'pair_class':pair,'step':step,
                                   'observed_variance':observed,'se':se,'theory':theory,'pass':passed})
                if step < 200:
                    q = oracle.advance_modes(q,dt,rates,variances,rng)
    assert all(row['pass'] for row in output), 'actual sampled stationary variance failed'
    return output


def adaptive_bead_contact_checks(cfg):
    setup = cfg['bead_crosscheck']
    n = setup['N']
    mapped = oracle.mapped_parameters(cfg,{'name':'adaptive_small_N','units_i':n,'units_j':n})
    elapsed = setup['steps']*setup['dt_reduced']
    run = {'name':'adaptive_contact_crosscheck','replicas':setup['replicas'],
           'box_over_scale':setup['box_over_scale'],'dt_reduced':setup['dt_reduced'],
           'max_dt_reduced':5*setup['dt_reduced'],'safety':10.,'window':[elapsed/2,elapsed]}
    output = []
    for index,pair in enumerate(cfg['pair_classes']):
        bead_count = oracle.independent_bead_contact(cfg,pair,setup['seed']+index)
        sampled = oracle.simulate(run,cfg,mapped,pair,212001+index)
        adaptive_count = setup['replicas']-sampled['survivors'][-1]
        pb,pa = bead_count/setup['replicas'],adaptive_count/setup['replicas']
        se = math.sqrt((pb*(1-pb)+pa*(1-pa))/setup['replicas'])
        euler_scale = 3*n/(cfg['chains'][str(n)]['bond_rms_m']/mapped['scale_m'])**2*4*setup['dt_reduced']
        bound = setup['z']*se+euler_scale*max(pb,pa)
        passed = abs(pb-pa) <= bound
        output.append({'pair_class':pair,'bead_contacts':bead_count,'actual_adaptive_contacts':adaptive_count,
                       'replicas':setup['replicas'],'bead_probability':pb,'actual_adaptive_probability':pa,
                       'combined_se':se,'euler_relative_scale':euler_scale,'bound':bound,'pass':passed,
                       'actual_hit_times_sha256':sampled['hit_times_sha256']})
    assert all(row['pass'] for row in output), 'actual adaptive trajectory / bead contacts failed'
    return output


def run_audits(cfg):
    return {'laplacian':direct_laplacian_checks(cfg),'continuum_sphere':continuum_sphere_check(cfg),
            'stationary_variance':stationary_variance_checks(cfg),
            'adaptive_bead_contact':adaptive_bead_contact_checks(cfg)}
