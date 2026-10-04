"""Equilibrium-state guard regression and joint-uncertainty arithmetic."""
import importlib.util
import os
from pathlib import Path
import sys

for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import numpy as np
import pytest

sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location('guard_oracle', Path(__file__).parent/'reference/run.py')
oracle = importlib.util.module_from_spec(spec)
spec.loader.exec_module(oracle)
cfg,_ = oracle.read_parameters()


@pytest.mark.parametrize('case',cfg['cases'],ids=lambda case:case['name'])
@pytest.mark.parametrize('pair',cfg['pair_classes'])
def test_equilibrium_guard(case,pair):
    mapped=oracle.mapped_parameters(cfg,case)
    rates,weights,var=oracle.mode_coefficients(cfg,mapped,pair)
    rng=np.random.Generator(np.random.PCG64(230001))
    count=4096
    q=rng.normal(size=(count,len(rates),3))*np.sqrt(var)[None,:,None]
    gap=10**rng.uniform(-4,2,count)
    for dt_min,max_dt,safety in ((1e-4,.2,6.),(6.25e-6,.05,8.),(6.25e-6,.025,10.)):
        proposed=np.minimum(max_dt,(gap/safety)**2/(2*mapped['D_local_reduced']))
        guarded=oracle.guard_observation_step(q,rates,weights,gap,proposed,dt_min,safety)
        mean=oracle.projected_site(q*np.expm1(-guarded[:,None]*rates)[:,:,None],weights)
        assert np.all((np.linalg.norm(mean,axis=1)<=gap/safety)|(guarded==dt_min))
        assert np.all(guarded<=np.maximum(proposed,dt_min))


def test_old_linear_shrink_violates_equilibrium_guard():
    mapped=oracle.mapped_parameters(cfg,next(case for case in cfg['cases'] if case['name']=='N16'))
    rates,weights,var=oracle.mode_coefficients(cfg,mapped,'mid/mid')
    rng=np.random.Generator(np.random.PCG64(230001))
    q=rng.normal(size=(4096,len(rates),3))*np.sqrt(var)[None,:,None]
    gap=10**rng.uniform(-4,2,4096)
    dt=np.minimum(.2,(gap/6)**2/(2*mapped['D_local_reduced']))
    shift=oracle.projected_site(q*np.expm1(-dt[:,None]*rates)[:,:,None],weights)
    dt=np.maximum(dt*np.minimum(1.,gap/(6*np.maximum(np.linalg.norm(shift,axis=1),1e-30))),1e-4)
    actual=oracle.projected_site(q*np.expm1(-dt[:,None]*rates)[:,:,None],weights)
    assert np.any((dt>1e-4)&(np.linalg.norm(actual,axis=1)>gap/6))


def test_linear_forms_represent_observed_contrasts():
    logs=np.random.Generator(np.random.PCG64(230002)).normal(size=(1000,8))
    explicit=(abs(logs[:,0]-logs[:,1])+abs(logs[:,1]-logs[:,2])+
              abs(logs[:,2]-logs[:,3])+np.maximum(abs(logs[:,2]-logs[:,4]),
              abs(logs[:,4]-logs[:,5]))+abs(logs[:,6]-logs[:,7]))
    forms=oracle.audit_linear_forms()
    assert forms.shape==(64,8)
    np.testing.assert_allclose((logs@forms.T).max(axis=1),explicit,rtol=1e-14,atol=1e-14)


def test_whole_half_covariance_is_preserved():
    mapped=oracle.mapped_parameters(cfg,cfg['cases'][0])
    hit=np.concatenate([np.full(1000,20.),np.full(800,30.),np.full(700,36.),np.full(7500,np.inf)])
    full=oracle.hazard(hit,20.,[24.,40.],mapped,cfg)
    halves=[oracle.hazard(hit,20.,[24.,32.],mapped,cfg),oracle.hazard(hit,20.,[32.,40.],mapped,cfg)]
    own={name:{'long':dict(full),'late_halves':halves} for name in ('base','step','contact','adaptive','box','box2')}
    logs,cov=oracle.audit_log_covariance(own)
    assert cov[5,6]>0 and cov[5,7]>0 and cov[6,7]==0
    assert np.linalg.eigvalsh(cov).min()>-1e-14
    result=oracle.joint_numerical_uncertainty(own,4.)
    assert result['audit_upper_log']>=result['observed_audit_contrasts_log']
    assert result['joint_normal_quantile']>4.
