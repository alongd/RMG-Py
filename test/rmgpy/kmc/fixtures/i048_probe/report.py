"""Render measured tables from saved, reproducible scientific source evidence.

Command: rmg_env python .../i048_probe/report.py --write-report
Command: rmg_env python .../i048_probe/report.py (check report)
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
import re
import time

from common import SCRATCH,HERE,PLAN,CPUS as PLAN_CPUS,bootstrap,load,save,digest,production_deadline

SPECIES,SEQUENCES=bootstrap()
REPORT=HERE.parent/'I048_oligomer_series.md'
START='<!-- BEGIN I048:measured -->'
END='<!-- END I048:measured -->'


def cost_snapshot():
    start=load(SCRATCH/'started.json')['unix_s']
    jobs=[]
    for path in SCRATCH.rglob('job*.json'):
        if path.name=='job.json' and list(path.parent.glob('job_*.json')):
            continue
        record=load(path)
        if record.get('started_unix_s',0)<start:
            continue
        resources=Path(record.get('resources_path',str(path.parent/'resources.txt')))
        text=resources.read_text() if resources.exists() else ''
        def field(pattern):
            found=re.search(pattern,text)
            return float(found.group(1)) if found else 0.
        jobs.append({'path':str(path),'sha256':digest(path),'label':record['label'],
                     'elapsed_s':record['elapsed_s'],'exit_code':record['exit_code'],
                     'user_CPU_s':field(r'User time \(seconds\):\s*([\d.]+)'),
                     'system_CPU_s':field(r'System time \(seconds\):\s*([\d.]+)'),
                     'max_RSS_kB':field(r'Maximum resident set size \(kbytes\):\s*([\d.]+)'),
                     'resources_path':str(resources),'resources_sha256':digest(resources) if resources.exists() else None})
    leaf=[j for j in jobs if not j['label'].startswith(('production stage ','replay '))]
    observations=[__import__('json').loads(line) for line in (SCRATCH/'resource_snapshots.jsonl').read_text().splitlines()]
    for observation in observations:
        if observation['RSS_kB']>16*1024**2:
            raise AssertionError('observed resource cap exceeded')
        for process in observation['processes']:
            if any(not set(mask)<=set(PLAN_CPUS) for mask in process['thread_affinities']):
                raise AssertionError('observed affinity cap exceeded')
    cutoff=time.time()
    extension_path=SCRATCH/'extension_authorization.json'
    extension=load(extension_path) if extension_path.exists() else None
    return {'cutoff_unix_s':cutoff,'elapsed_wall_s':cutoff-start,
            'recorded_jobs':jobs,'CPU_s':sum(j['user_CPU_s']+j['system_CPU_s'] for j in leaf),
            'max_single_job_RSS_kB':max([j['max_RSS_kB'] for j in jobs] or [0.]),
            'resource_observations':len(observations),
            'last_resource_observation_unix_s':observations[-1]['unix_s'] if observations else None,
            'max_observed_aggregate_RSS_kB':max([o['RSS_kB'] for o in observations] or [0.]),
            'owner_extension':extension,
            'owner_extension_sha256':digest(extension_path) if extension else None,
            'extension_CPU_s':sum(j['user_CPU_s']+j['system_CPU_s'] for j in leaf
                                 if extension and load(Path(j['path']))['started_unix_s']>=extension['recorded_unix_s']),
            'successful_snapshot_masks_within_declared_cores':True,
            'development_execution_exception':load(SCRATCH/'development_thread_cap_exception.json'),
            'development_execution_exception_sha256':digest(SCRATCH/'development_thread_cap_exception.json'),
            'workflow_queue_restart':load(SCRATCH/'workflow_queue_restart.json') if (SCRATCH/'workflow_queue_restart.json').exists() else None,
            'workflow_queue_restart_sha256':digest(SCRATCH/'workflow_queue_restart.json') if (SCRATCH/'workflow_queue_restart.json').exists() else None,
            'timeout_guard_holds':[{'path':str(p),'sha256':digest(p),'record':load(p)}
                                   for p in sorted((SCRATCH/'ensembles').glob('ps*/timeout_hold.json'))],
            'timeout_queue_restart':{'path':str(SCRATCH/'workflow_timeout_queue_restart.json'),
                                     'sha256':digest(SCRATCH/'workflow_timeout_queue_restart.json'),
                                     'record':load(SCRATCH/'workflow_timeout_queue_restart.json')}
                                    if (SCRATCH/'workflow_timeout_queue_restart.json').exists() else None,
            'captured_resource_violations':[{'path':str(p),'sha256':digest(p),'record':load(p)} for p in sorted(SCRATCH.glob('resource_violation_*.json'))],
            'scope':'recorded leaf scientific jobs through this cutoff; enclosing supervisors excluded to avoid double counting CPU; build and small authoring/audit processes excluded; copied prior results incurred no new quantum cost'}


def render(data,base,baseline):
    lines=[]
    def out(text=''):
        lines.append(text)
    out('Pinned baseline Tc: **%.6f K** at 1 mol/L. Database commit `%s`; snapshot `%s`, %d allowlisted files.'%(
        baseline['continuous_Tc_K'],baseline['database_SHA'],baseline['snapshot']['sha256'],baseline['snapshot']['files']))
    out()
    out('### Sequence and molecular coverage')
    out()
    missing=sorted(set(SPECIES)-set(data['available_species']))
    out('Thermochemical coverage: **%d/%d species**. Missing completed thermochemistry: %s.'%(
        len(data['available_species']),len(SPECIES),', '.join(missing) or 'none'))
    out('An increment is calculated only when every frozen class at both lengths and every reference is complete. Missing classes are never omitted and the remaining weights are never renormalized.')
    new_names=[name for name in SPECIES if name.startswith('ps')]
    searches=sum((SCRATCH/'ensembles'/name/'completed.json').exists() for name in new_names)
    pools=sum((SCRATCH/'composite'/name/'selection.json').exists() for name in new_names)
    out('New-case source coverage: %d/%d completed searches and %d/%d selected minimum pools. The table retains available search and minimum evidence even when a case lacks completed thermochemistry.'%(
        searches,len(new_names),pools,len(new_names)))
    out()
    out('| n | Class | Frozen atactic weight | Oriented assignments | CREST candidates | Within 12 kJ/mol | Checked candidates | Rotor wells | Minimum ESS at 298 K | Largest-population basin ESS: PBE / BLYP |')
    out('| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |')
    for n,classes in SEQUENCES.items():
        for name,row in classes.items():
            if name not in data['molecules']:
                rec_path=SCRATCH/'ensembles'/name/'completed.json'
                pool_path=SCRATCH/'composite'/name/'selection.json'
                reduction_path=SCRATCH/'composite'/name/'candidate_reduction.json'
                rec=load(rec_path) if rec_path.exists() else None
                pool=load(pool_path) if pool_path.exists() else None
                checked=(sum((SCRATCH/'xtb_checks'/name/f'{i:04d}'/'result.json').exists()
                             for i in load(reduction_path)['retained_indices'])
                         if reduction_path.exists() else 'unavailable')
                values=[n,name,'%.6f'%row['weight'],len(row['assignments']),
                        rec['frames'] if rec else 'unavailable',
                        len(rec['selected_indices_within_12_kJ']) if rec else 'unavailable',
                        checked,len(pool['unique_minima_indices']) if pool else 'unavailable',
                        'unavailable','unavailable']
                out('| '+' | '.join(map(str,values))+' |')
                continue
            rec=load(SCRATCH/'ensembles'/name/'completed.json')
            pool=load(SCRATCH/'composite'/name/'selection.json')
            checked=(len(load(SCRATCH/'composite'/name/'candidate_reduction.json')['retained_indices'])
                     if name.startswith('ps') else len(pool['checked']))
            diagnostics=data['molecules'][name]['basin_298_internal_diagnostics']
            dominant=[max(diagnostics,key=lambda r:r['probability'][level])['ESS']
                      for level in ('pbe','blyp')]
            out('| %s | %s | %.6f | %d | %d | %d | %d | %d | %.2f | %.2f / %.2f |'%(
                n,name,row['weight'],len(row['assignments']),rec['frames'],
                len(rec['selected_indices_within_12_kJ']),checked,len(pool['unique_minima_indices']),
                data['molecules'][name]['lowest_ESS_at_298'],*dominant))
    out()
    out(('All stereochemical classes are covered through n=5. ' if data['complete'] else
         'The declared full n=5 series is incomplete; the 4→5 atactic increment is unavailable. ')+
        'Complete lengths are averaged with their frozen assignment weights, without equilibrium diastereomer mixing entropy. The population spread below measures variation among oriented addition channels; it is not a standard error on these exact weighted means.')
    out()
    out('Monte Carlo standard errors condition on the retained wells and declared partition model. They do not quantify omitted conformers, candidate truncation, rigidity of the coupled-rotor potential, or the electronic fallback. Sequence population spread is reported separately. These quantities alone therefore do not provide a total uncertainty or prove physical convergence; low effective sample sizes also limit the linearized Monte Carlo error estimate.')
    out('The minimum effective sample size can belong to a basin of small probability. The largest-population basin column and the per-basin internal probabilities saved in JSON provide the population context; none of these diagnostics changes the declared quadrature.')
    out()
    out('### Increment and residual table at 298.15 K')
    out()
    out('H increments use the one-ring formation-enthalpy anchor at 298 K described above, followed by the chain thermal increment. S increments are absolute molecular differences; GAV S and the raw comparison receive matching oriented-end normalization. The component tables retain the original molecular raw sums. Local values come from the internal chain partition. Balanced-cycle residuals are separately retained in JSON and in the component diagnostics.')
    out()
    out('| Step | Level | GAV ΔH (kJ/mol) | GAV oriented ΔS (J/mol/K) | Calibrated raw ΔH | Raw oriented ΔS | Calibrated local ΔH | Internal ΔS | Raw δH | Raw δS | Local δH | Local δS | Local MC SE H/S | Sequence population SD H/S |')
    out('| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |')
    for n,levels in data['steps'].items():
        ga=baseline['increments'][n]['raw'][0]
        for level in ('pbe','blyp'):
            r=levels[level]['rows'][0]
            raw=r['raw_delta']; local=r['local_delta']
            err=r['local_MC_errors_H_S_G']; spread=r['sequence_population_spread_H_S']
            out('| %s→%d | %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f / %.3f | %.3f / %.3f |'%(
                n,int(n)+1,level,ga['H_J_mol']/1000,ga['S_J_mol_K']+r['GAV_chain_orientation_shift_S_J_mol_K'],
                (ga['H_J_mol']+raw['H_J_mol'])/1000,r['chain_raw']['S_J_mol_K']+r['QM_chain_orientation_shift_S_J_mol_K'],
                (ga['H_J_mol']+local['H_J_mol'])/1000,r['chain_local_partition']['S_J_mol_K'],
                raw['H_J_mol']/1000,raw['S_J_mol_K'],local['H_J_mol']/1000,local['S_J_mol_K'],
                err[0]/1000,err[1],spread[0]/1000,spread[1]))
    out()
    out('### Temperature-dependent raw and local corrections')
    out()
    out('Only new n=4 and 5 quadratures change in the 1024/2048 comparison; the reused shorter ensembles retain their original production counts.')
    out()
    out('| Step | Level | T (K) | GAV ΔH / oriented ΔS | Calibrated raw ΔH / local ΔH | Raw δH (kJ/mol) | Raw δS (J/mol/K) | Local δH | Local δS | MC SE H/S/G (kJ/mol, J/mol/K, kJ/mol) | Sequence SD H/S | Change from 1024 to 2048 points: local H/S |')
    out('| --- | --- | ---: | --- | --- | ---: | ---: | ---: | ---: | --- | --- | --- |')
    for n,levels in data['steps'].items():
        for level in ('pbe','blyp'):
            for j,r in enumerate(levels[level]['rows']):
                b=base['steps'][n][level]['rows'][j]
                raw=r['raw_delta']; local=r['local_delta']; err=r['local_MC_errors_H_S_G']; sd=r['sequence_population_spread_H_S']
                ga=r['GAV_oriented_increment']
                out('| %s→%d | %s | %.2f | %.3f / %.3f | %.3f / %.3f | %.3f | %.3f | %.3f | %.3f | %.3f / %.3f / %.3f | %.3f / %.3f | %.3f / %.3f |'%(
                    n,int(n)+1,level,r['T_K'],ga['H_J_mol']/1000,ga['S_J_mol_K'],
                    (ga['H_J_mol']+raw['H_J_mol'])/1000,(ga['H_J_mol']+local['H_J_mol'])/1000,
                    raw['H_J_mol']/1000,raw['S_J_mol_K'],local['H_J_mol']/1000,local['S_J_mol_K'],
                    err[0]/1000,err[1],err[2]/1000,sd[0]/1000,sd[1],
                    (local['H_J_mol']-b['local_delta']['H_J_mol'])/1000,local['S_J_mol_K']-b['local_delta']['S_J_mol_K']))
    out()
    out('### Entropy component increments of the chains')
    out()
    out('These are chain increments before the balanced-reference subtraction. Local partitions remove external factors before reweighting conformer basins; their internal components need not equal the gas-weight component subtraction. No class-mixing bit is added.')
    out()
    out('| Step | Level | T (K) | Partition | Translation | External rotation | Non-torsional vibration | Coupled rotors | Basin mixing | Sum (J/mol/K) |')
    out('| --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |')
    for n,levels in data['steps'].items():
        for level in ('pbe','blyp'):
            for r in levels[level]['rows']:
                for label,key in (('raw','chain_raw'),('local','chain_local_partition')):
                    row=r[key]; c=row['S_components_J_mol_K']
                    out('| %s→%d | %s | %.2f | %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |'%(
                        n,int(n)+1,level,r['T_K'],label,c['translation'],c['external_rotation'],
                        c['non_torsional_vibration'],c['coupled_rotors'],c['basin_mixing'],row['S_J_mol_K']))
    out()
    out('### Enthalpy component increments of the chains')
    out()
    out('The electronic column includes the constant formation-enthalpy reference determined by the balanced cycle at 298 K. The uncalibrated electronic increment and all balanced-cycle components are also preserved in JSON. Translation and rotation each have zero H increment for one nonlinear chain growing into another. Basin mixing contributes entropy rather than a separate enthalpy; electronic ensemble averaging includes population energy shifts.')
    out()
    out('| Step | Level | T (K) | Partition | Electronic + formation anchor | Translation | External rotation | Non-torsional vibration | Coupled rotors | Basin mixing | Sum (kJ/mol) |')
    out('| --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |')
    for n,levels in data['steps'].items():
        for level in ('pbe','blyp'):
            for r in levels[level]['rows']:
                for label,key in (('raw','chain_raw'),('local','chain_local_partition')):
                    row=r[key];c=row['H_components_J_mol']
                    out('| %s→%d | %s | %.2f | %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |'%(
                        n,int(n)+1,level,r['T_K'],label,(c['electronic']+r['formation_H_anchor_J_mol'])/1000,c['translation']/1000,
                        c['external_rotation']/1000,c['non_torsional_vibration']/1000,c['coupled_rotors']/1000,c['basin_mixing']/1000,
                        (row['H_J_mol']+r['formation_H_anchor_J_mol'])/1000))
    out()
    out('### Balanced-reference diagnostics and finite-fragment conditional ceilings')
    out()
    out('| Step | Level | Chain translation at 298 | Chain rotation at 298 | Reference-side translation | Reference-side rotation | Balanced raw δS | Balanced local δS | Component-subtracted balanced δS | All-cycle internal δS | Raw conditional Tc (K) | Local conditional Tc (K) |')
    out('| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |')
    for n,levels in data['steps'].items():
        for level in ('pbe','blyp'):
            r=levels[level]['rows'][0]; c=r['chain_raw']['S_components_J_mol_K'];ref=r['balanced_reference_side']['S_components_J_mol_K']
            roots=levels[level].get('finite_fragment_conditional_Tc',{'raw':[],'local':[]})
            def fmt(values,local=False):
                return ', '.join('%.3f ± %.3f (MC)'%(v['Tc_K'],v['MC_SE_K'])+
                    ('; channel SD %.3f K (linearized)'%v['sequence_local_spread_linearized_K'] if local else '')
                    for v in values) or 'none in 298–800 K'
            out('| %s→%d | %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %s | %s |'%(
                n,int(n)+1,level,c['translation'],c['external_rotation'],ref['translation'],ref['external_rotation'],
                r['balanced_raw_delta']['S_J_mol_K'],r['balanced_local_delta']['S_J_mol_K'],
                r['gas_weight_component_subtraction_delta']['S_J_mol_K'],r['all_cycle_external_removal_delta']['S_J_mol_K'],fmt(roots['raw']),fmt(roots['local'],True)))
    out()
    out('### Functional spread and length dependence')
    out()
    out('| Step | T (K) | PBE−BLYP local δH (kJ/mol) | PBE−BLYP local δS (J/mol/K) |')
    out('| --- | ---: | ---: | ---: |')
    for n,levels in data['steps'].items():
        for p,b in zip(levels['pbe']['rows'],levels['blyp']['rows']):
            out('| %s→%d | %.2f | %.3f | %.3f |'%(n,int(n)+1,p['T_K'],
                (p['local_delta']['H_J_mol']-b['local_delta']['H_J_mol'])/1000,
                p['local_delta']['S_J_mol_K']-b['local_delta']['S_J_mol_K']))
    out()
    out('| Level | T (K) | Terminal change: (4→5)−(3→4), local δH (kJ/mol) | Terminal local δS change (J/mol/K) | Covariance-aware MC SE H/S | Paired sequence population SD H/S |')
    out('| --- | ---: | ---: | ---: | --- | --- |')
    for level in ('pbe','blyp','gfn2','pbe_sparse','blyp_sparse'):
        for row in data['convergence_comparison'].get(level,[]):
            err=row['MC_SE_H_S_G']; sd=row['sequence_population_SD_H_S']
            out('| %s | %.2f | %.3f | %.3f | %.3f / %.3f | %.3f / %.3f |'%(level,row['T_K'],
                row['delta_H_change_J_mol']/1000,row['delta_S_change_J_mol_K'],err[0]/1000,err[1],sd[0]/1000,sd[1]))
    out()
    out('### Convergence verdict and transfer')
    out()
    if missing:
        out('The completed classes support only the fully covered increments listed above. Missing n=5 classes prevent the terminal length comparison. No converged long-chain correction or corresponding long-chain Tc is established. Finite-fragment conditional roots, when computed, are comparisons under an assumed transfer and cannot replace the missing convergence evidence.')
    else:
        out('The terminal length changes are listed with covariance-aware Monte Carlo errors and paired sequence spread. Three finite increments, sparse electronic re-ranking, candidate truncation and one doubled quadrature do not establish an asymptotic plateau with a bounded total error. No validated long-chain correction or Tc is adopted; the terminal 4→5 result and its conditional roots remain finite-fragment estimates.')
    if {'2','3'}<=set(data['steps']):
        changes=[]
        for level in ('pbe','blyp'):
            before=data['steps']['2'][level]['rows'][0]['local_delta']
            after=data['steps']['3'][level]['rows'][0]['local_delta']
            sparse_before=data['steps']['2'][level+'_sparse']['rows'][0]['local_delta']
            sparse_after=data['steps']['3'][level+'_sparse']['rows'][0]['local_delta']
            changes.append('%s: production length change %.3f kJ/mol and %.3f J/mol/K; uniform sparse diagnostic %.3f kJ/mol and %.3f J/mol/K'%(
                level.upper(),(after['H_J_mol']-before['H_J_mol'])/1000,
                after['S_J_mol_K']-before['S_J_mol_K'],
                (sparse_after['H_J_mol']-sparse_before['H_J_mol'])/1000,
                sparse_after['S_J_mol_K']-sparse_before['S_J_mol_K']))
        out('At 298.15 K, comparing (3→4) with (2→3): '+ '; '.join(changes)+'. The enthalpy trend changes sign under the already-declared energy-treatment diagnostic. The length trend therefore cannot be interpreted independently of that approximation boundary; neither diagnostic replaces the production result.')
    out('The internal chain increment is the candidate for transfer because whole-chain translation and external rotation do not accompany local growth of a macroscopic chain. Removing these factors before basin reweighting differs from subtracting gas-weighted component totals. The reported correction also contains intrarepeat QM-versus-GAV differences; it is not an isolated adjacent-phenyl pair interaction.')
    out('No n=6 calculation was allocated: required n≤5 searches, both frozen quadratures, composite energies and verification had priority within the original 48-hour cap. The owner-authorized extension explicitly prohibits n=6. This budget decision uses measured computational cost and outstanding required work, without inspecting a ceiling temperature.' if (SCRATCH/'extension_authorization.json').exists() else 'No n=6 calculation was allocated: required n≤5 searches, both frozen quadratures, composite energies and verification had priority within the original 48-hour cap. This budget decision uses measured computational cost and outstanding required work, without inspecting a ceiling temperature.')
    out()
    out('### Electronic approximation sensitivity')
    out()
    out('The pure-GFN2 comparison changes electronic energies only, retaining the same geometries, Hessians and quadratures. It is a sensitivity diagnostic, not an error bar on the composite. The range of DFT−GFN2 offsets among the few calculated minima tests the constant-offset assumption locally; it cannot bound offsets for omitted conformers.')
    out()
    out('| Class | PBE / BLYP calculated offset range (kJ/mol) | Approximate internal population at 298 K: PBE / BLYP | At 800 K: PBE / BLYP |')
    out('| --- | --- | --- | --- |')
    for name,molecule in data['molecules'].items():
        if not name.startswith('ps'):
            continue
        ranges=molecule['DFT_minus_GFN2_offset_range_J_mol']; pop=molecule['approximate_internal_population']
        out('| %s | %.3f / %.3f | %.6f / %.6f | %.6f / %.6f |'%(name,ranges['pbe']/1000,ranges['blyp']/1000,
            pop['pbe'][0]['fraction'],pop['blyp'][0]['fraction'],pop['pbe'][-1]['fraction'],pop['blyp'][-1]['fraction']))
    out()
    out('| Step | T (K) | Pure-GFN2 local δH / δS (kJ/mol, J/mol/K) | PBE−GFN2 local δH / δS | BLYP−GFN2 local δH / δS |')
    out('| --- | ---: | --- | --- | --- |')
    for n,levels in data['steps'].items():
        for j,g in enumerate(levels['gfn2']['rows']):
            gv=g['local_delta']; p=levels['pbe']['rows'][j]['local_delta']; b=levels['blyp']['rows'][j]['local_delta']
            out('| %s→%d | %.2f | %.3f / %.3f | %.3f / %.3f | %.3f / %.3f |'%(n,int(n)+1,g['T_K'],
                gv['H_J_mol']/1000,gv['S_J_mol_K'],(p['H_J_mol']-gv['H_J_mol'])/1000,p['S_J_mol_K']-gv['S_J_mol_K'],
                (b['H_J_mol']-gv['H_J_mol'])/1000,b['S_J_mol_K']-gv['S_J_mol_K']))
    out()
    out('The following diagnostic applies the same lowest-three-plus-offset rule to the complete n=2 and 3 electronic evidence, keeping the production one-ring reference fixed. The larger-chain energies already obey that rule. Production minus this uniform sparse diagnostic therefore measures sensitivity to the energy-treatment boundary; it does not change the reported production increments or select another method.')
    out()
    out('| Step | Level | T (K) | Production−uniform sparse local δH (kJ/mol) | Production−uniform sparse local δS (J/mol/K) |')
    out('| --- | --- | ---: | ---: | ---: |')
    for n,levels in data['steps'].items():
        for level in ('pbe','blyp'):
            for j,p in enumerate(levels[level]['rows']):
                s=levels[level+'_sparse']['rows'][j]['local_delta']; q=p['local_delta']
                out('| %s→%d | %s | %.2f | %.3f | %.3f |'%(n,int(n)+1,level,p['T_K'],
                    (q['H_J_mol']-s['H_J_mol'])/1000,q['S_J_mol_K']-s['S_J_mol_K']))
    out()
    out('### Literature cross-check')
    out()
    lit=load(SCRATCH/'literature_crosscheck.json')
    out('**Complete copies of the two early RIS papers could not be fetched.** The reached primary indexed excerpts support the specific parameters below; a full RIS entropy-per-dyad reconstruction was not reproduced. No literature parameter was fitted.')
    out()
    y=lit['Yoon_1975']
    out('Yoon, Sundararajan and Flory give rounded relative statistical-weight prefactors 0.8 and 1.3. Their well-shape entropy terms, calculated here as R ln(prefactor), are %.6f and %.6f J/mol/K. These small local terms are not an absolute per-dyad entropy or a GAV residual. [Original study, p. 781](%s).'%(y['R_ln_eta_prefactor_J_mol_K'],y['R_ln_omega_prefactor_J_mol_K'],y['indexed_primary_pdf']))
    out()
    out('Relative RIS weights determine state-probability ratios. Multiplying every transfer-matrix weight by the same temperature-independent constant leaves those probabilities unchanged but adds R ln(constant) per step to the partition-derived entropy. An absolute intrawell reference and matching basin definitions are therefore needed to compare them with the present continuous coupled-rotor and whole QM-minus-GAV increment. A discrete-state Shannon entropy alone does not supply that reference. This is a mathematical reference limitation, independent of the unavailable full-paper copies.')
    out()
    w=lit['Williams_Flory_1969']
    out('Williams and Flory report a relative conformer preference of −700 cal/mol with an opposing entropy difference of −1.4 cal/mol/K: %.6f kJ/mol and %.6f J/mol/K after unit conversion. Phenyl rotational restriction has a negative relative entropy contribution. This relative-conformer quantity does not establish the sign or magnitude of an entire chain-increment correction against GAV. [Original study, discussion after eq. 28](%s).'%(w['relative_H_kJ_mol'],w['relative_S_J_mol_K'],w['indexed_primary_pdf']))
    out()
    k=lit['Khare_Paulaitis_1994']
    out('Khare and Paulaitis explicitly study coupled phenyl/backbone motions in polystyrene hexamers, supporting the use of coupled torsions. The accessible abstract supplies no absolute entropy-per-dyad comparator. The literature therefore supports the mechanism of coupling and relative entropy penalties, while agreement in sign and magnitude with the computed QM-minus-GAV chain increment is not established from the sources reached. [Primary university record and abstract](%s).'%k['accessible_primary_abstract'])
    out()
    out('### Recorded computation cost')
    out()
    cost=load(SCRATCH/'cost_snapshot.json')
    if cost['timeout_guard_holds']:
        out('An eight-hour per-search timeout was an execution limit added by this worker, rather than the dispatched total-wall limit. Only each recorded GNU timeout wrapper was held while its unchanged CREST child continued. Original-deadline guardians preserve the actual wrapper status and scientific-child evidence separately. Guard outcomes: %s. Successful acceptance requires time-v exit 0, normal CREST termination and source hashes; an unfinished child at the original cap is not accepted. A controlled sleep process reproduced a successful child with wrapper status 124. Preparation-supervisor reload evidence and prior logs remain in scratch. Neither a wrapper hold nor a lane recovery extends the original 48-hour deadline.'%('; '.join('%s: scientific completed=%s, wrapper exit=%s'%(item['record']['name'],item['record']['scientific_completed'],item['record'].get('wrapper_exit_code','pending')) for item in cost['timeout_guard_holds'])))
        out()
    out('Through the recorded cutoff: %.3f h wall since preparation; %.3f CPU h in %d recorded jobs; maximum single-job RSS %.3f GiB. %s.'%(
        cost['elapsed_wall_s']/3600,cost['CPU_s']/3600,len(cost['recorded_jobs']),
        cost['max_single_job_RSS_kB']/1024**2,cost['scope']))
    deadline=load(SCRATCH/'started.json')['unix_s']+48*3600
    out('The original scientific deadline is 2026-10-05 12:16:40 UTC. The cost cutoff is %.3f h %s that deadline. Report rendering and cached numerical audits after the cap do not allocate additional quantum time.'%(
        abs(cost['cutoff_unix_s']-deadline)/3600,'before' if cost['cutoff_unix_s']<deadline else 'after'))
    if cost.get('owner_extension'):
        extension=cost['owner_extension']
        out('The owner subsequently authorized a hard 16-hour extension from ruling time, with absolute deadline %s (epoch %d), to finish n=5 from local baseline `%s` (22/25 cases). This authorization supersedes the original production cutoff for the resumed jobs; the original start and scientific method declaration are unchanged. The two interrupted searches restarted with the installed CREST 3.0.2 GFN2/quick/6 kcal/mol/four-thread settings. A requested native restart did not recover their checkpoint stages: actual streams show new metadynamics, so these are fresh same-settings restarts. Prior trees and launch receipts are preserved. No per-case eight-hour timer or n=6 allocation is used. Extension leaf scientific CPU recorded since its initialization: %.3f h. Cached rendering/audits after the new cap allocate no quantum time.'%(extension['deadline_UTC'],extension['deadline_unix_s'],extension['baseline_SHA'],cost['extension_CPU_s']/3600))
    out()
    out('%d successful resource snapshots observed an aggregate calculation RSS peak of %.3f GiB with threads confined to the eight declared physical cores. Monitoring exceptions and a gap are described below; these periodic observations are not a continuous peak-memory or affinity proof. Production DFT uses three lanes on 3/3/2 distinct physical cores, native OpenMP bounded by each lane, and one thread in each BLAS pool. The workspace setting is 5000 MB per job; observed RSS is measured separately.'%(
        cost['resource_observations'],cost['max_observed_aggregate_RSS_kB']/1024**2))
    out()
    out('Completed cases received minima checks and electronic calculations while remaining searches continued on the same eight pinned physical cores. The primary pipeline prepared one case at a time. A supplemental producer later used a second preparation lane, with minimum checks on core 10 and DFT on cores 6, 8 and 10. Per-case process locks serialize minimum and electronic writes. The supplemental producer finishes its active case when all searches complete; the bulk electronic stage waits on its process lock before launching the declared three lanes. The already-active ps4_0000 was excluded from supplemental scheduling because it began before the case locks were installed. This scheduling overlap did not change candidate selection, energy levels, sequence weights or sampling counts.')
    out('The final bulk electronic stage distributes independent conformer/functional single points across those same three lanes, including when only one case remains. It holds the case locks, retains the frozen lowest-three selection, and applies the existing offset and symmetry-image rules only after all direct points complete. The saved task-to-lane matrix and producing jobs are audited; a development check also compares cached energy-file hashes before and after this scheduling refactor.')
    out()
    if cost['workflow_queue_restart'] is not None:
        out('The primary queue was restarted after its preparation CLI was confirmed to have no calculation children and to be waiting on the supplemental-owned ps4_0110 minimum lock. The interrupted waiting-stage exit 143, timing record and original logs were preserved; scientific producers continued. The restarted queue checks case ownership and its children return status 75 if another producer acquired the case. Such a case remains pending. `workflow_queue_restart.json` records this change; the original wall deadline is retained.')
        out()
    out('The first prepared tetramer also received its declared rotor integrations on core 6 during the searches. A per-case process lock protects these integrations and later cached bulk replay from overlapping writes. This added no physical core or scientific reduction.')
    out()
    out('A later early-rotor producer integrated prepared cases on core 12 during the searches. It finishes its active case when all searches complete, then releases a process lock required by the bulk eight-lane rotor stage. Both declared sampling counts are retained.')
    out()
    layouts={}
    for path in (SCRATCH/'composite').glob('ps*/checks/execution.json'):
        execution=load(path);nlanes=len(execution['lane_CPUs'])
        layouts[nlanes]=layouts.get(nlanes,0)+1
    if layouts:
        out('Later primary minimum preparation divided independent candidates between two single-core lanes; bulk preparation assigned the same eight cores among cases still missing their pools. The original supplemental process retained its single check lane. Recorded stage layouts: '+', '.join('%d case(s) with %d lane(s)'%(count,nlanes) for nlanes,count in sorted(layouts.items()))+'. Each candidate batch and its core is recorded in `checks/execution.json`; the Verifier checks the complete unchanged candidate set and single-thread execution metadata. Cached candidates retain their original computation metadata. No selection, threshold, Hessian treatment or energy level changes.')
        out()
    out('One development partial-thermochemistry replay omitted explicit thread environment settings; its actual BLAS use was not measured. The monitor stopped during that replay on an affinity escape whose argv/mask were not captured. Its first restart identified an unpinned tee logger, which was corrected. Quantum-job environments inspected during this period had the caps set. The partial analysis was repeated with caps, and analysis entry points now set them before scientific imports. `development_thread_cap_exception.json` records the exceptions and monitoring gap; the final report uses capped reproductions.')
    return '\n'.join(lines).strip()


PREFACE='''# I048 — polystyrene oligomer per-unit increment probe

This report distinguishes numerical replay of a declared approximate partition
from physical convergence and transfer to a radical polymer chain. No correction
is adopted. Only probe scripts and this report are changed.
The [reproduction commands](i048_probe/README.md) include the full
`verify_results.py --replay` source audit and numerical replay.

## Declared reductions and reference definitions

The previous probe's methyl cap convention is retained:
`CH3–[CH(Ph)–CH2]_(n−1)–CH(Ph)–CH3`, formula `C(8n+1)H(8n+4)`.
This differs from literal hydrogen-capped `H–(CH2–CHPh)n–H` by the same
terminal carbon at every length; successive increments still add C8H8.
The existing n=2 and 3 molecular evidence and the four one-ring reference
molecules are copied from the prior named scratch directory, with byte hashes
and source provenance. The n=4 and 5 inputs and searches are new.

All distinct stereo classes through n=5 are enumerated from all oriented
binary assignments and weighted as a frozen unbiased atactic chain.
Enantiomer equivalence does not imply an equilibrium mixture of sequence
classes. External rotational symmetry, local optical factors, global mirror
partners and chain configuration entropy are kept separate. Existing RMG
propagation configuration entropy is retained once.

The original reduction declaration is retained in
`method_plan_before_local_clarification.json`; `method_plan_i048.json` contains
the same reductions plus the recorded partition/reference clarifications.
After replaying the earlier short-fragment cycle, the entropy comparison was
clarified to use its absolute increment and the H reference was held constant
at 298 K. These clarifications preserve H/S consistency and precede any Tc
calculation; they do not change sampling, candidate choices, or energy levels.
The reductions declared before increment calculations
are: CREST/GFN2 `--quick`, 6 kcal/mol search window; up to 32 candidates within
12 kJ/mol for each new sequence, chosen by GFN2 energy; tight extreme xTB
optimization and Hessians, symmetry deduplication and proper/reflection orbit
completion; 1024 and 2048 correlated-importance points per basin; DFT on the
lowest three chemical minima for each larger sequence. Other retained minima
use GFN2 energies plus the functional offset of the lowest GFN2 conformer.
The declared levels are PBE-D3(BJ) and BLYP-D3(BJ)/def2-SVP//GFN2, density
fitting, grid level 3, SCF tolerance 1e-9 Eh, no optional three-body D3 term.
No reduction or ranking uses Tc or a thermochemical comparator.
The same already-declared sparse energy rule was subsequently evaluated on
the complete short-chain evidence as a diagnostic of the treatment boundary.
It does not replace any production energy or increment; all retained prior
and new source evidence is unchanged.

The exact prior fixed-angle coupled-rotor model is reused. Every acyclic
heavy-atom single bond is a torsional coordinate. Projection removes that
subspace before non-torsional harmonic vibrations are counted. The rigid
torsional Hessian supplies the Pitzer–Gwinn harmonic reference. Periodic
Voronoi basins partition the full torsional domain once, with methyl and
phenyl fundamental periods. Severe contacts below 0.45 Å have zero weight
and remain in the Monte Carlo denominator; other unresolved SCC failures
stop the producer. The stochastic CREST search is not proved converged.

## Reference and the meaning of local

Electronic energies of different atom counts cannot be subtracted from GAV
formation enthalpies directly. Following the prior probe, use the balanced
cycle `C_n + cumene + n-propylbenzene → C_(n+1) + ethylbenzene + ethane`.
Its additive source vector cancels at every length. The balanced QM-minus-GAV
H residual at 298.15 K anchors the formation-enthalpy increment to the pinned
GAV one-ring reference; the anchor is held constant and the chain thermal
increment supplies its temperature dependence. Absolute entropy increments
need no electronic-energy reference and are compared directly to oriented
GAV. This preserves the H/S derivative identity. The temperature-dependent
balanced-cycle H/S residuals are retained as separate diagnostics reproducing
the prior probe's comparison. This is not an independently determined absolute
formation enthalpy or an explicit styrene-addition calculation.

The raw chain partition includes ideal-gas translation and external rotation
at 1 bar. The primary local partition omits both factors before basin
probabilities are calculated. Thus it describes internal chain conformations,
rather than retaining an inertia bias from gas rotational weights. The gas
one-ring reference used for the H anchor remains fixed. Component subtraction at the gas weights
and removal from every reference molecule are reported as diagnostics.
GAV chain entropies are normalized to an oriented fixed sequence by undoing
finite-end symmetry and removing its constitutional optical factors; no new
R ln 2 is inserted. JSON preserves every molecular component and channel.
The identical additive increment at each length defines the group model's
per-unit baseline. The local residual compares that fixed increment with
the internal QM chain increment; the raw residual retains the measured
finite-chain external factors.
The total residual also contains errors in intrarepeat vibrations and torsions.
The balanced-cycle diagnostics retain the original GAV cycle reference; their
different quantum reference partitions describe separate comparisons. The
direct local increment relative to the fixed GAV per-unit baseline supplies
the conditional chain-transfer calculation.

Whole-chain translation and rotation change as mass and moments of inertia
change; their increments tend to zero in a long-chain limit. It is therefore
the chain external increment that can be removed to estimate that limit.
The prior cycle's negative translation/rotation terms largely belong to its
small-molecule reference. Deleting those terms changes the reference, rather
than merely removing a finite-chain artifact. This is an evidenced problem
with the motivating interpretation in the dispatch.

The conditional gas roots solve
`ΔH_p(T)+δH(T) − T[ΔS_p(T)+δS(T)] − RT ln(c0 RT/p0) = 0`,
with c0=1000 mol/m³ and p0=100000 Pa. Only the measured 298–800 K domain is
used. The pinned propagation H/S are evaluated on a one-kelvin grid and
interpolated for root finding. Each finite-fragment root is a comparison
under its transfer assumption, not a validated long-chain ceiling.

## Measured evidence

'''


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true')
    parser.add_argument('--refresh-cost',action='store_true')
    parser.add_argument('--available-only',action='store_true',help='render incomplete coverage explicitly; does not satisfy the full Verifier')
    args=parser.parse_args()
    stem='increments_available_' if args.available_only else 'increments_'
    data=load(SCRATCH/(stem+'2048.json'))
    base=load(SCRATCH/(stem+'1024.json'))
    baseline=load(SCRATCH/'baseline/baseline.json')
    if not data['complete'] and not args.available_only:
        raise AssertionError('cannot render a complete report from incomplete series')
    if args.refresh_cost or not (SCRATCH/'cost_snapshot.json').exists():
        save(SCRATCH/'cost_snapshot.json',cost_snapshot())
    rendered=render(data,base,baseline)
    if args.write_report:
        text=PREFACE+START+'\n'+rendered+'\n'+END+'\n'
        REPORT.write_text(text)
        print('I048 report written from measured evidence',flush=True)
    else:
        text=REPORT.read_text().split(START,1)[1].split(END,1)[0].strip()
        if text!=rendered:
            raise AssertionError('report differs from sources')
        print('I048 report numbers match measured evidence',flush=True)


if __name__=='__main__':
    main()
