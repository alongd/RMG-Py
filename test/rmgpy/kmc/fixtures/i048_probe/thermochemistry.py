"""Fixed-sequence oligomer increments with explicit reference and components.

Command: PYTHONPATH=$PWD rmg_env python .../i048_probe/thermochemistry.py
Command: same, --available-only (partial scientific evidence, no completion claim)
"""
from __future__ import annotations

import argparse
import math
from functools import lru_cache
from pathlib import Path
import sys

from common import SCRATCH, PLAN, TEMPERATURES, bootstrap, load, save

import numpy as np
from scipy.special import logsumexp
from scipy.optimize import brentq

SPECIES, SEQUENCES=bootstrap()
# This file intentionally has the same descriptive name as the old module;
# load the scientific implementation explicitly to avoid import shadowing.
import importlib.util
spec=importlib.util.spec_from_file_location('i043_thermal',Path(__file__).resolve().parent.parent/'i043_probe/thermochemistry.py')
old=importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
R=old.R
REFERENCE={'cumene':-1.,'n_propylbenzene':-1.,'ethylbenzene':1.,'ethane':1.}


def chemical_ancestor(name,index):
    while True:
        record=load(SCRATCH/'xtb_checks'/name/f'{index:04d}'/'result.json')
        parent=record.get('derived_by_reflection_from',record.get('derived_by_permutation_from'))
        if parent is None:
            return index
        index=parent


class Ensemble(old.Ensemble):
    @lru_cache(maxsize=256)
    def hs(self,temperature,level,internal=False):
        data=[old.conformer_thermal(record,temperature,self.symmetry) for record in self.records]
        energy=self.energy[level]
        low=float(min(energy))
        logq=np.array([row['log_q'] for row in data])
        thermal=np.array([row['H_thermal_J_mol'] for row in data])
        entropy=np.array([row['S_J_mol_K'] for row in data])
        external=np.array([row['S_translation_J_mol_K']+row['S_external_rotation_J_mol_K'] for row in data])
        if internal:
            logq -= external/R-4.
            thermal -= 4*R*temperature
            entropy -= external
        logw=logq-(energy-low)/(R*temperature)
        p=np.exp(logw-logsumexp(logw))
        total_h=energy+thermal
        h=float(p@total_h)
        mixing=-R*float(p@np.log(np.maximum(p,1e-300)))
        s=float(p@entropy+mixing)
        sc={label:float(p@np.array([row[key] for row in data])) for label,key in (
            ('translation','S_translation_J_mol_K'),('external_rotation','S_external_rotation_J_mol_K'),
            ('non_torsional_vibration','S_non_torsional_vibration_J_mol_K'),('coupled_rotors','S_coupled_rotor_J_mol_K'))}
        if internal:
            sc['translation']=sc['external_rotation']=0.
        sc['basin_mixing']=mixing
        hc={'electronic':float(p@energy),'translation':0. if internal else 2.5*R*temperature,
            'external_rotation':0. if internal else 1.5*R*temperature,
            'non_torsional_vibration':float(p@np.array([row['H_thermal_J_mol']-row['H_tor_J_mol']-4*R*temperature for row in data])),
            'coupled_rotors':float(p@np.array([row['H_tor_J_mol'] for row in data])), 'basin_mixing':0.}
        if abs(sum(hc.values())-h)>1e-4 or abs(sum(sc.values())-s)>1e-7:
            raise AssertionError('component sum mismatch')
        variances=np.zeros(3)
        for pi,hi,si,row in zip(p,total_h,entropy,data):
            influence_h=pi*(row['H_influence']+(hi-h)*row['log_influence'])
            influence_s=pi*(row['S_influence']+(si-R*math.log(max(pi,1e-300))-s)*row['log_influence'])
            influence_g=-R*temperature*pi*row['log_influence']
            variances += np.array([np.var(i,ddof=1)/len(i) for i in (influence_h,influence_s,influence_g)])
        return {'T_K':temperature,'H_J_mol':h,'S_J_mol_K':s,'H_components_J_mol':hc,
                'S_components_J_mol_K':sc,'MC_errors_H_S_G':np.sqrt(variances).tolist(),
                'probabilities':p.tolist(),'minimum_ESS':min(row['effective_samples'] for row in data)}


def make_ensembles(baseline,new_samples=2048,available_only=False):
    ensembles={}
    for name in SPECIES:
        prior=name in PLAN['prior_samples']
        counts={int(i):n for i,n in PLAN['prior_basin_overrides'].get(name,{}).items()}
        proposal='correlated' if name.startswith('ps') or name.startswith('diphenylpentane') or name.startswith('triphenylheptane') else 'diagonal'
        try:
            ensemble=Ensemble(SCRATCH,name,baseline['species'][name],
                PLAN['prior_samples'].get(name,new_samples),PLAN['seed'],proposal,counts)
            ensemble.energy['gfn2']=np.array([load(SCRATCH/'xtb_checks'/name/f'{i:04d}'/'result.json')['energy_Eh'] for i in ensemble.indices])*old.EH_MOL
            roots=[chemical_ancestor(name,i) for i in ensemble.indices]
            # Apply the already-declared larger-chain energy rule to the full
            # shorter-chain evidence as a diagnostic only. Production energies
            # and the fixed one-ring H reference remain unchanged.
            for level in ('pbe','blyp'):
                if name.startswith('ps') or name in REFERENCE:
                    ensemble.energy[level+'_sparse']=ensemble.energy[level].copy()
                    continue
                chemical=sorted(set(roots))
                gfn={i:load(SCRATCH/'xtb_checks'/name/f'{i:04d}'/'result.json')['energy_Eh'] for i in chemical}
                selected=sorted(chemical,key=gfn.get)[:PLAN['DFT_chemical_minima_per_sequence']]
                quantum={i:load(SCRATCH/'composite'/name/f'{i:04d}'/level/'result.json')['energy_hartree'] for i in selected}
                offset=quantum[selected[0]]-gfn[selected[0]]
                ensemble.energy[level+'_sparse']=np.array([quantum[i] if i in quantum else gfn[i]+offset for i in roots])*old.EH_MOL
            ensembles[name]=ensemble
        except FileNotFoundError:
            if not available_only:
                raise
    return ensembles


def sum_rows(rows,coefficients):
    result={'H_J_mol':0.,'S_J_mol_K':0.,'H_components_J_mol':{},'S_components_J_mol_K':{}}
    variance=np.zeros(3)
    for name,c in coefficients.items():
        row=rows[name]
        for key in ('H_J_mol','S_J_mol_K'):
            result[key]+=c*row[key]
        for key in ('H_components_J_mol','S_components_J_mol_K'):
            for part,value in row[key].items():
                result[key][part]=result[key].get(part,0.)+c*value
        variance += c*c*np.array(row['MC_errors_H_S_G'])**2
    result['MC_errors_H_S_G']=np.sqrt(variance).tolist()
    return result


def gav_orientation(baseline,ensemble,name):
    row=baseline['species'][name]
    end_sigma=ensemble.symmetry['nonoptical_RMG_sigma']/ensemble.symmetry['internal_rotor_divisor']
    return R*math.log(end_sigma)-R*math.log(2)*row['optical_atom_half_factors']


def evaluate_step(ensembles,baseline,n,t,level):
    weights={name:-r['weight'] for name,r in SEQUENCES[str(n)].items()}
    weights.update({name:r['weight'] for name,r in SEQUENCES[str(n+1)].items()})
    rows={name:e.hs(t,level) for name,e in ensembles.items() if name in weights or name in REFERENCE}
    internal={name:ensembles[name].hs(t,level,True) for name in weights}
    chain=sum_rows(rows,weights)
    local_chain=sum_rows(internal,weights)
    ref=sum_rows(rows,REFERENCE)
    raw_cycle=sum_rows(rows,dict(REFERENCE,**weights))
    local_rows=dict(rows)
    local_rows.update(internal)
    local_cycle=sum_rows(local_rows,dict(REFERENCE,**weights))
    # Balanced additive source vectors are zero here, so H/Cp are identically
    # zero and S is constant. Do not extrapolate a nonzero-source case.
    model=baseline['increments'][str(n)]
    if model['balanced_source_vector']:
        raise ValueError('nonzero GAV source vector needs continuous thermo evaluator')
    gh=model['balanced'][0]['H_J_mol']
    gs=model['balanced'][0]['S_J_mol_K']
    gav_shift=sum(c*gav_orientation(baseline,ensembles[name],name) for name,c in weights.items())
    qm_shift=sum(c*R*math.log(ensembles[name].symmetry['external_divisor']) for name,c in weights.items())
    normalized_gav=gs+gav_shift
    corrected_raw={'H_J_mol':raw_cycle['H_J_mol']-gh,
                   'S_J_mol_K':raw_cycle['S_J_mol_K']+qm_shift-normalized_gav}
    corrected_local={'H_J_mol':local_chain['H_J_mol']+ref['H_J_mol']-gh,
                     'S_J_mol_K':local_chain['S_J_mol_K']+ref['S_J_mol_K']-normalized_gav}
    stripped={'H_J_mol':chain['H_J_mol']+ref['H_J_mol']-gh,
              'S_J_mol_K':chain['S_J_mol_K']-chain['S_components_J_mol_K']['translation']
               -chain['S_components_J_mol_K']['external_rotation']+ref['S_J_mol_K']-normalized_gav}
    all_local=dict(internal)
    all_local.update({name:ensembles[name].hs(t,level,True) for name in REFERENCE})
    all_cycle=sum_rows(all_local,dict(REFERENCE,**weights))
    all_removed={'H_J_mol':all_cycle['H_J_mol']-gh,'S_J_mol_K':all_cycle['S_J_mol_K']-normalized_gav}
    local_variance=np.array(local_chain['MC_errors_H_S_G'])**2+np.array(ref['MC_errors_H_S_G'])**2
    lookup={length:{bits:name for name,r in classes.items() for bits in r['assignments']} for length,classes in SEQUENCES.items()}
    channels=[]
    for bits,child in lookup[str(n+1)].items():
        parent=lookup[str(n)][bits[:-1]]
        ga_shift=gav_orientation(baseline,ensembles[child],child)-gav_orientation(baseline,ensembles[parent],parent)
        raw_model_h=sum(c*baseline['species'][name]['thermo'][0]['H_J_mol'] for name,c in dict(REFERENCE,**{parent:-1.,child:1.}).items())
        raw_model_s=sum(c*baseline['species'][name]['thermo'][0]['S_J_mol_K'] for name,c in dict(REFERENCE,**{parent:-1.,child:1.}).items())
        channels.append({'assignment':bits,'parent':parent,'child':child,
            'delta_H_J_mol':internal[child]['H_J_mol']-internal[parent]['H_J_mol']+ref['H_J_mol']-raw_model_h,
            'delta_S_J_mol_K':internal[child]['S_J_mol_K']-internal[parent]['S_J_mol_K']+ref['S_J_mol_K']-raw_model_s-ga_shift})
    spread=[float(np.std([r[k] for r in channels])) for k in ('delta_H_J_mol','delta_S_J_mol_K')]
    if abs(np.mean([r['delta_H_J_mol'] for r in channels])-corrected_local['H_J_mol'])>1e-4 or abs(np.mean([r['delta_S_J_mol_K'] for r in channels])-corrected_local['S_J_mol_K'])>1e-7:
        raise AssertionError('oriented channel average differs from class-weighted increment')
    # Absolute entropy needs no electronic zero. Anchor formation H at 298 K
    # and retain chain Cp thereafter; this preserves the H/S derivative identity.
    dense=baseline['dense_GAV_increments'][str(n)]
    grid=np.array([r['T_K'] for r in dense])
    gav_h=float(np.interp(t,grid,[r['H_J_mol'] for r in dense]))
    gav_s=float(np.interp(t,grid,[r['S_J_mol_K'] for r in dense]))+gav_shift
    ref298=sum_rows({name:ensembles[name].hs(298.15,level) for name in REFERENCE},REFERENCE)
    gav298=model['raw'][0]['H_J_mol']
    h_anchor=ref298['H_J_mol']-gh+gav298
    direct_raw={'H_J_mol':chain['H_J_mol']+h_anchor-gav_h,
                'S_J_mol_K':chain['S_J_mol_K']+qm_shift-gav_s}
    direct_local={'H_J_mol':local_chain['H_J_mol']+h_anchor-gav_h,
                  'S_J_mol_K':local_chain['S_J_mol_K']-gav_s}
    direct_var=np.array(local_chain['MC_errors_H_S_G'])**2
    direct_var[0]+=ref298['MC_errors_H_S_G'][0]**2
    direct_var[2]+=ref298['MC_errors_H_S_G'][0]**2
    channel_shift_h=direct_local['H_J_mol']-corrected_local['H_J_mol']
    channel_shift_s=direct_local['S_J_mol_K']-corrected_local['S_J_mol_K']
    direct_channels=[dict(r,delta_H_J_mol=r['delta_H_J_mol']+channel_shift_h,
                         delta_S_J_mol_K=r['delta_S_J_mol_K']+channel_shift_s) for r in channels]
    return {'T_K':t,'chain_raw':chain,'chain_local_partition':local_chain,'balanced_reference_side':ref,
            'raw_balanced_cycle':raw_cycle,'local_balanced_cycle':local_cycle,
            'raw_delta':direct_raw,'local_delta':direct_local,
            'balanced_raw_delta':corrected_raw,'balanced_local_delta':corrected_local,
            'GAV_oriented_increment':{'H_J_mol':gav_h,'S_J_mol_K':gav_s},
            'formation_H_anchor_J_mol':h_anchor,
            'gas_weight_component_subtraction_delta':stripped,'all_cycle_external_removal_delta':all_removed,
            'QM_chain_orientation_shift_S_J_mol_K':qm_shift,'GAV_chain_orientation_shift_S_J_mol_K':gav_shift,
            'local_MC_errors_H_S_G':np.sqrt(direct_var).tolist(),
            'raw_MC_errors_H_S_G':np.sqrt(np.array(chain['MC_errors_H_S_G'])**2+
                np.array([ref298['MC_errors_H_S_G'][0]**2,0.,ref298['MC_errors_H_S_G'][0]**2])).tolist(),
            'balanced_local_MC_errors_H_S_G':np.sqrt(local_variance).tolist(),
            'sequence_population_spread_H_S':spread,'oriented_channels':direct_channels}


def compute(available_only=False,new_samples=2048,roots=True):
    baseline=load(SCRATCH/'baseline/baseline.json')
    ensembles=make_ensembles(baseline,new_samples,available_only)
    result={'available_species':sorted(ensembles),'complete':set(ensembles)==set(SPECIES),
            'new_samples':new_samples,'sequence_catalogue':SEQUENCES,'steps':{},'molecules':{},
            'convergence_comparison':{}}
    for name,e in ensembles.items():
        approximated=[]
        offset_spread={}
        for level in ('pbe','blyp'):
            original=[i for i in e.indices if i<10000]
            offsets=[]
            for i in original:
                quantum=load(SCRATCH/'composite'/name/f'{i:04d}'/level/'result.json')
                if not quantum.get('approximation'):
                    offsets.append(quantum['energy_hartree']-load(SCRATCH/'xtb_checks'/name/f'{i:04d}'/'result.json')['energy_Eh'])
            offset_spread[level]=(max(offsets)-min(offsets))*old.EH_MOL if offsets else 0.
        for i in e.indices:
            i=chemical_ancestor(name,i)
            approximated.append(bool(load(SCRATCH/'composite'/name/f'{i:04d}'/'pbe/result.json').get('approximation')))
        result['molecules'][name]={'indices':e.indices,'symmetry':e.symmetry,
            'lowest_ESS_at_298':e.hs(298.15,'pbe')['minimum_ESS'],
            'basin_298_internal_diagnostics':[
                {'index':index,
                 'ESS':old.conformer_thermal(record,298.15,e.symmetry)['effective_samples'],
                 'probability':{level:e.hs(298.15,level,True)['probabilities'][position]
                                for level in ('pbe','blyp')}}
                for position,(index,record) in enumerate(zip(e.indices,e.records))],
            'DFT_minus_GFN2_offset_range_J_mol':offset_spread,
            'approximate_internal_population':{level:[{'T_K':t,'fraction':float(np.array(e.hs(t,level,True)['probabilities'])@np.array(approximated))} for t in TEMPERATURES] for level in ('pbe','blyp')}}
    for n in range(2,5):
        needed=set(SEQUENCES[str(n)])|set(SEQUENCES[str(n+1)])|set(REFERENCE)
        if not needed<=set(ensembles):
            continue
        result['steps'][str(n)]={}
        for level in ('pbe','blyp','gfn2','pbe_sparse','blyp_sparse'):
            rows=[evaluate_step(ensembles,baseline,n,t,level) for t in TEMPERATURES]
            result['steps'][str(n)][level]={'rows':rows}
            if roots and level in ('pbe','blyp'):
                grid=baseline['propagation']
                tx=np.array([r['T_K'] for r in grid])
                hx=np.array([r['H_J_mol'] for r in grid])
                sx=np.array([r['S_J_mol_K'] for r in grid])
                values={}
                for kind in ('raw','local'):
                    def free(t):
                        row=evaluate_step(ensembles,baseline,n,t,level)
                        dh=row[kind+'_delta']['H_J_mol']
                        ds=row[kind+'_delta']['S_J_mol_K']
                        return np.interp(t,tx,hx)+dh-t*(np.interp(t,tx,sx)+ds)-R*t*math.log(1000*R*t/100000)
                    coarse=np.linspace(298.15,800.,26)
                    brackets=[(a,b) for a,b in zip(coarse[:-1],coarse[1:]) if free(a)*free(b)<0]
                    values[kind]=[]
                    for a,b in brackets:
                        root=brentq(free,a,b,xtol=1e-6)
                        row=evaluate_step(ensembles,baseline,n,root,level)
                        slope=(free(root+.01)-free(root-.01))/.02
                        spread=np.std([c['delta_H_J_mol']-root*c['delta_S_J_mol_K'] for c in row['oriented_channels']])
                        values[kind].append({'Tc_K':float(root),'residual_J_mol':float(free(root)),
                            'MC_SE_K':row[kind+'_MC_errors_H_S_G'][2]/abs(slope),
                            'sequence_local_spread_linearized_K':float(spread/abs(slope))})
                result['steps'][str(n)][level]['finite_fragment_conditional_Tc']=values
    if {'3','4'}<=set(result['steps']):
        coefficients={name:r['weight'] for name,r in SEQUENCES['5'].items()}
        coefficients.update({name:-2*r['weight'] for name,r in SEQUENCES['4'].items()})
        coefficients.update({name:r['weight'] for name,r in SEQUENCES['3'].items()})
        lookup={n:{bits:name for name,r in classes.items() for bits in r['assignments']} for n,classes in SEQUENCES.items()}
        for level in ('pbe','blyp','gfn2','pbe_sparse','blyp_sparse'):
            compared=[]
            for j,t in enumerate(TEMPERATURES):
                before=result['steps']['3'][level]['rows'][j]['local_delta']
                after=result['steps']['4'][level]['rows'][j]['local_delta']
                rows={name:ensembles[name].hs(t,level,True) for name in coefficients}
                difference=sum_rows(rows,coefficients)
                # The shared H anchor and identical oriented GAV increments
                # cancel. Use covariance-aware coefficients, especially -2
                # for the common n=4 ensembles, rather than independent errors.
                dh=after['H_J_mol']-before['H_J_mol']
                ds=after['S_J_mol_K']-before['S_J_mol_K']
                if abs(dh-difference['H_J_mol'])>1e-3 or abs(ds-difference['S_J_mol_K'])>1e-6:
                    raise AssertionError('terminal comparison reference cancellation mismatch')
                channels=[]
                for bits,name5 in lookup['5'].items():
                    name4=lookup['4'][bits[:-1]]
                    name3=lookup['3'][bits[:-2]]
                    channels.append([rows[name5][key]-2*rows[name4][key]+rows[name3][key] for key in ('H_J_mol','S_J_mol_K')])
                compared.append({'T_K':t,'delta_H_change_J_mol':dh,'delta_S_change_J_mol_K':ds,
                    'MC_SE_H_S_G':difference['MC_errors_H_S_G'],
                    'sequence_population_SD_H_S':np.std(channels,axis=0).tolist()})
            result['convergence_comparison'][level]=compared
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--available-only',action='store_true')
    parser.add_argument('--samples',type=int,default=2048)
    parser.add_argument('--no-roots',action='store_true')
    args=parser.parse_args()
    if load(SCRATCH/'method_plan_i048.json')!=PLAN:
        raise AssertionError('method plan not declared')
    result=compute(args.available_only,args.samples,not args.no_roots)
    path=SCRATCH/('increments_partial.json' if args.available_only else f'increments_{args.samples}.json')
    save(path,result)
    if args.available_only:
        save(SCRATCH/f'increments_available_{args.samples}.json',result)
    for n,levels in result['steps'].items():
        for level,record in levels.items():
            row=record['rows'][0]
            print('I048 %s->%d %s local delta H298 %.6f kJ/mol, S298 %.6f J/mol/K; MC %.3f kJ/mol %.3f J/mol/K; sequence SD %.3f %.3f'%(
                n,int(n)+1,level,row['local_delta']['H_J_mol']/1000,row['local_delta']['S_J_mol_K'],
                row['local_MC_errors_H_S_G'][0]/1000,row['local_MC_errors_H_S_G'][1],
                row['sequence_population_spread_H_S'][0]/1000,row['sequence_population_spread_H_S'][1]),flush=True)
    print('I048 thermochemistry coverage: '+str(len(result['available_species']))+'/'+str(len(SPECIES))+' species; complete='+str(result['complete']),flush=True)


if __name__=='__main__':
    main()
