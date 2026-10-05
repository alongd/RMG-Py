"""Audit source evidence independently and replay the committed producers.

Command: PYTHONPATH=$PWD rmg_env python .../i048_probe/verify_results.py --replay
No fresh-search completeness or physical convergence follows from this audit.
"""
from __future__ import annotations

import argparse
import importlib.util
import math
import os
from pathlib import Path
import re
import sys
import time

from common import (SCRATCH,HERE,LEGACY,PLAN,CPUS,RMG_PYTHON,DFT_PYTHON,
                    TEMPERATURES,bootstrap,digest,load,save)

import numpy as np
from rdkit import Chem

SPECIES,SEQUENCES=bootstrap()


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


thermal=module('i048_thermal',HERE/'thermochemistry.py')
legacy=module('i043_audit',LEGACY/'verify_results.py')
from rotors import make_molecule,angles,wrap,rotated,XTB,quadrature_filename,EH_J,NA,rotor_definitions,projected_modes
from summarize import stereo_check
from run_ensembles import read_xyz
R=thermal.R


def close(a,b,tolerance=1e-6):
    if not math.isclose(a,b,rel_tol=1e-10,abs_tol=tolerance):
        raise AssertionError('independent verification mismatch: %r != %r'%(a,b))


def compare(a,b,path='root'):
    if isinstance(a,dict):
        if set(a)!=set(b):
            raise AssertionError('key mismatch at '+path)
        for key in a:
            compare(a[key],b[key],path+'.'+key)
    elif isinstance(a,list):
        if len(a)!=len(b):
            raise AssertionError('length mismatch at '+path)
        for i,(x,y) in enumerate(zip(a,b)):
            compare(x,y,path+'.'+str(i))
    elif isinstance(a,(int,float)) and not isinstance(a,bool):
        close(a,b,1e-5)
    elif a!=b:
        raise AssertionError('value mismatch at '+path)


def audit_symmetry(source,minimum,smiles):
    key=('derived_by_reflection_from' if 'derived_by_reflection_from' in minimum else 'derived_by_permutation_from')
    if key not in minimum:
        return
    reflection=key=='derived_by_reflection_from'
    parent=source.parent/f"{minimum[key]:04d}"
    before=np.array([[float(v) for v in line.split()[1:]] for line in read_xyz(parent/'xtbopt.xyz')[0]['atoms']])
    actual=np.array([[float(v) for v in line.split()[1:]] for line in read_xyz(source/'xtbopt.xyz')[0]['atoms']])
    mapping=minimum['reflection_atom_permutation' if reflection else 'permutation_atom_order']
    expected=before.copy()
    if reflection:
        expected[:,0]*=-1
    if not np.allclose(actual,expected[mapping],rtol=0,atol=1e-6):
        raise AssertionError('symmetry geometry transform mismatch')
    molecule=Chem.AddHs(Chem.MolFromSmiles(smiles))
    for atom in molecule.GetAtoms():
        if atom.GetDegree()!=4:
            continue
        order=sorted(a.GetIdx() for a in atom.GetNeighbors())
        p=before[order]; q=actual[order]
        if np.linalg.det(p[:3]-p[3])*np.linalg.det(q[:3]-q[3])<=0:
            raise AssertionError('symmetry transform changed signed tetrahedral orientation')
    def hessian(path):
        return np.array([float(v) for line in path.read_text().splitlines() if not line.startswith('$') for v in line.split()]).reshape(3*len(before),3*len(before))
    if minimum['parent_hessian_sha256']!=digest(parent/'hessian'):
        raise AssertionError('symmetry Hessian parent stale')
    permutation=[3*i+a for i in mapping for a in range(3)]
    expected=hessian(parent/'hessian')[np.ix_(permutation,permutation)]
    if reflection:
        signs=np.tile([-1.,1.,1.],len(before))
        expected*=signs[:,None]*signs[None,:]
    if not np.allclose(expected,hessian(source/'hessian'),rtol=0,atol=1e-12):
        raise AssertionError('symmetry Hessian transform mismatch')


def audit_rotor_orbit(name,ensemble,manifest):
    """Check reflection/proper closure independently of the saved image list."""
    molecules=[make_molecule(SPECIES[name],SCRATCH/'xtb_checks'/name/f'{i:04d}'/'xtbopt.xyz')
               for i in ensemble.indices]
    definitions=rotor_definitions(molecules[0])
    periods=np.array([2*math.pi/r['symmetry'] for r in definitions])
    centers=np.array([wrap(angles(m,definitions),periods) for m in molecules])
    expected=Chem.AddHs(Chem.MolFromSmiles(SPECIES[name]))
    for molecule in molecules:
        xyz=np.array(molecule.GetConformer().GetPositions())
        tetra=[sorted(a.GetIdx() for a in atom.GetNeighbors())
               for atom in molecule.GetAtoms() if atom.GetDegree()==4]
        def volume(coordinates,order):
            points=coordinates[order]
            return np.linalg.det(points[:3]-points[3])
        volumes=[volume(xyz,order) for order in tetra]
        def check_image(coordinates,mapping,tolerance):
            image=Chem.Mol(molecule)
            for i,point in enumerate(coordinates[list(mapping)]):
                image.GetConformer().SetAtomPosition(i,point)
            point=wrap(angles(image,definitions),periods)
            distance=np.min(np.max(np.abs(wrap(point-centers,periods)),axis=1))
            if distance>tolerance:
                raise AssertionError('symmetry-related rotor well omitted: '+name)
        if manifest[name]['achiral']:
            reflected=xyz.copy();reflected[:,0]*=-1
            image=Chem.Mol(molecule)
            for i,point in enumerate(reflected):
                image.GetConformer().SetAtomPosition(i,point)
            Chem.AssignStereochemistryFrom3D(image,replaceExistingTags=True)
            matches=image.GetSubstructMatches(expected,useChirality=True,uniquify=False,maxMatches=1000000)
            if len(matches)>=1000000:
                raise AssertionError('reflection enumeration truncated')
            mapping=next((p for p in matches if all(v*volume(reflected,[p[i] for i in order])>0
                         for v,order in zip(volumes,tetra))),None)
            if mapping is None:
                raise AssertionError('achiral reflection changed tetrahedral labels')
            check_image(reflected,mapping,.10)
        if ensemble.symmetry['external_divisor']>1:
            matches=molecule.GetSubstructMatches(molecule,useChirality=True,uniquify=False,maxMatches=1000000)
            if len(matches)>=1000000:
                raise AssertionError('proper permutation enumeration truncated')
            for mapping in matches:
                if all(v*volume(xyz,[mapping[i] for i in order])>0 for v,order in zip(volumes,tetra)):
                    check_image(xyz,mapping,.050001)


def audit_partial_sources(complete_names):
    """Check report rows with search/minimum evidence but no thermochemistry."""
    from run_series import search_succeeded
    rows=[]
    for name in sorted(set(SPECIES)-set(complete_names)):
        directory=SCRATCH/'ensembles'/name
        path=directory/'completed.json'
        if not path.exists():
            continue
        search=load(path)
        frames=read_xyz(directory/'crest_conformers.xyz')
        energies=[float(frame['comment'].split()[0]) for frame in frames]
        if len(frames)!=search['frames'] or digest(directory/'crest_conformers.xyz')!=search['ensemble_sha256']:
            raise AssertionError('partial search frame count/hash differs')
        compare(energies,search['energy_hartree'])
        selected=[i for i,e in enumerate(energies) if (e-min(energies))*2625.499639<12.]
        compare(selected,search['selected_indices_within_12_kJ'])
        for frame in frames:
            stereo_check(SPECIES[name],frame)
        job=load(directory/'job.json')
        if not search_succeeded(name,job):
            raise AssertionError('partial search scientific child failed')
        receipt=directory/'timeout_hold.json'
        if receipt.exists():
            proof=load(receipt)
            if proof['job_sha256']!=digest(directory/'job.json') or proof['wrapper_exit_code']!=job['exit_code']:
                raise AssertionError('partial held-wrapper receipt differs')
        row={'name':name,'frames':len(frames),'within_12_kJ':len(selected),
             'completed_sha256':digest(path),'checked_candidates':None,'rotor_wells':None}
        reduction_path=SCRATCH/'composite'/name/'candidate_reduction.json'
        if reduction_path.exists():
            reduction=load(reduction_path)
            compare(reduction['retained_indices'],sorted(selected,key=lambda i:energies[i])[:PLAN['new_candidate_cap']])
            row['checked_candidates']=sum((SCRATCH/'xtb_checks'/name/f'{i:04d}'/'result.json').exists()
                                          for i in reduction['retained_indices'])
        pool_path=SCRATCH/'composite'/name/'selection.json'
        if pool_path.exists():
            indices=load(pool_path)['unique_minima_indices']
            if len(indices)!=len(set(indices)):
                raise AssertionError('partial minimum pool repeats a well')
            for index in indices:
                source=SCRATCH/'xtb_checks'/name/f'{index:04d}'
                minimum=load(source/'result.json')
                if digest(source/'xtbopt.xyz')!=minimum['xyz_sha256'] or not minimum['minimum_confirmed'] or min(minimum['frequencies_cm1'])<-1:
                    raise AssertionError('partial minimum is changed or unstable')
                stereo_check(SPECIES[name],read_xyz(source/'xtbopt.xyz')[0])
                audit_symmetry(source,minimum,SPECIES[name])
                legacy.verify_frequency_units(source,SPECIES[name],minimum)
            row.update(rotor_wells=len(indices),selection_sha256=digest(pool_path))
        rows.append(row)
    save(SCRATCH/'verification/partial_source_audit.json',rows)
    if rows:
        print('I048 partial search/minimum source rows reproduced: '+', '.join(row['name'] for row in rows)+'; no missing thermochemistry supplied',flush=True)


def audit_sources(ensembles):
    from run_series import search_succeeded
    rotor_execution=SCRATCH/'pipeline/bulk_rotors_execution.json'
    if rotor_execution.exists():
        execution=load(rotor_execution)
        cpus=execution['lane_CPUs']
        if cpus!=list(CPUS) or execution['science_changes'] is not False:
            raise AssertionError('bulk rotor allocation differs')
        keys=[(t['name'],t['index'],t['samples']) for t in execution['tasks']]
        if len(keys)!=len(set(keys)):
            raise AssertionError('bulk rotor basin duplicated')
        for task in execution['tasks']:
            name=task['name']
            if name not in ensembles:
                continue
            selection_path=SCRATCH/'composite'/name/'selection.json'
            if digest(selection_path)!=task['selection_sha256'] or task['index'] not in load(selection_path)['unique_minima_indices']:
                raise AssertionError('bulk rotor partition changed')
            if task['samples'] not in PLAN['new_rotor_counts']:
                raise AssertionError('bulk rotor sampling count changed')
            job=load(SCRATCH/'rotor_jobs'/name/'basins'/str(task['samples'])/f"{task['index']:04d}"/'job.json')
            if job['exit_code']!=0 or job['cpus']!=[cpus[task['lane']]]:
                raise AssertionError('bulk rotor did not complete on its lane')
    execution_path=SCRATCH/'pipeline/bulk_electronic_execution.json'
    if execution_path.exists():
        execution=load(execution_path)
        lanes=execution['lane_CPUs']
        flattened=[cpu for lane in lanes for cpu in lane]
        if len(lanes)!=3 or sorted(flattened)!=list(CPUS) or execution['science_changes'] is not False:
            raise AssertionError('bulk single-point allocation differs')
        keys=[(t['name'],t['index'],t['level']) for t in execution['tasks']]
        if len(keys)!=len(set(keys)):
            raise AssertionError('bulk single point duplicated')
        for task in execution['tasks']:
            name=task['name']
            if name not in ensembles:
                continue
            selected=load(SCRATCH/'composite'/name/'electronic_reduction.json')['DFT_indices']
            if task['index'] not in selected or task['level'] not in ('pbe','blyp'):
                raise AssertionError('bulk scheduling changed electronic selection')
            job=load(SCRATCH/'composite'/name/f"{task['index']:04d}"/task['level']/'job.json')
            if job['exit_code']!=0 or job['cpus']!=lanes[task['lane']]:
                raise AssertionError('bulk single point did not complete on its lane')
    manifest=load(SCRATCH/'species.json')
    count=points_count=0
    streams={}
    for name,ensemble in ensembles.items():
        directory=SCRATCH/'ensembles'/name
        search=load(directory/'completed.json')
        if digest(directory/'crest_conformers.xyz')!=search['ensemble_sha256']:
            raise AssertionError('search hash mismatch')
        frames=read_xyz(directory/'crest_conformers.xyz')
        if len(frames)!=search['frames']:
            raise AssertionError('search frame count mismatch')
        energies=[float(f['comment'].split()[0]) for f in frames]
        compare(energies,search['energy_hartree'])
        selection=[i for i,e in enumerate(energies) if (e-min(energies))*2625.499639<12.]
        if selection!=search['selected_indices_within_12_kJ']:
            raise AssertionError('candidate energy rule mismatch')
        for frame in frames:
            stereo_check(SPECIES[name],frame)
        if name.startswith('ps'):
            job=load(directory/'job.json')
            if not search_succeeded(name,job):
                raise AssertionError('search scientific child did not succeed: '+name)
            if (directory/'timeout_hold.json').exists():
                proof=load(directory/'timeout_hold.json')
                if proof['job_sha256']!=digest(directory/'job.json') or proof['wrapper_exit_code']!=job['exit_code']:
                    raise AssertionError('held timeout status was not preserved: '+name)
            reduction=load(SCRATCH/'composite'/name/'candidate_reduction.json')
            if reduction['retained_indices']!=sorted(selection,key=lambda i:energies[i])[:32]:
                raise AssertionError('predeclared reduction rule mismatch')
            execution_path=SCRATCH/'composite'/name/'checks/execution.json'
            if execution_path.exists():
                execution=load(execution_path)
                cpus=execution['lane_CPUs'];batches=execution['candidate_batches']
                if len(cpus)!=len(batches) or len(set(cpus))!=len(cpus) or not set(cpus)<=set(CPUS):
                    raise AssertionError('minimum lanes exceed their allocation')
                if sorted(i for batch in batches for i in batch)!=sorted(reduction['retained_indices']):
                    raise AssertionError('minimum lanes changed or duplicated the candidate set')
                for batch in batches:
                    for i in batch:
                        candidate=load(SCRATCH/'xtb_checks'/name/f'{i:04d}'/'result.json')
                        adapter=candidate['execution_adapter']
                        if adapter['OMP_NUM_THREADS']!=1 or adapter['BLAS_threads']!=1 or len(adapter['cpu_affinity'])!=1 or not set(adapter['cpu_affinity'])<=set(CPUS):
                            raise AssertionError('candidate check exceeded its thread/core allocation')
                        command=candidate['command']
                        if command[command.index('--parallel')+1]!='1':
                            raise AssertionError('candidate xTB native thread count differs')
            chemical=[i for i in ensemble.indices if i<10000]
            gfn={i:load(SCRATCH/'xtb_checks'/name/f'{i:04d}'/'result.json')['energy_Eh'] for i in chemical}
            selected=sorted(chemical,key=gfn.get)[:PLAN['DFT_chemical_minima_per_sequence']]
            compare(load(SCRATCH/'composite'/name/'electronic_reduction.json'),
                    {'DFT_indices':selected,'GFN2_indices':sorted(set(chemical)-set(selected))})
            for i in chemical:
                for level in ('pbe','blyp'):
                    quantum=load(SCRATCH/'composite'/name/f'{i:04d}'/level/'result.json')
                    if bool(quantum.get('approximation'))!=(i not in selected):
                        raise AssertionError('direct DFT/fallback assignment differs from declared rule')
        audit_rotor_orbit(name,ensemble,manifest)
        for index,record in zip(ensemble.indices,ensemble.records):
            if record['seed'] in streams:
                raise AssertionError('MC variance assumes independent basins but RNG seeds collide: '+str((streams[record['seed']],(name,index))))
            streams[record['seed']]=(name,index)
            if name.startswith('ps') and record['samples']==2048:
                lower=load(SCRATCH/'rotors'/name/f'{index:04d}'/quadrature_filename(1024,PLAN['seed'],'correlated'))
                if lower['samples']!=1024:
                    raise AssertionError('declared lower sampling count differs')
                for key in ('seed','input_sha256','definitions','periods_rad','modes','basin_centers_rad','unique_minima_indices','coordinates_A','importance_parameters','importance_proposal','hard_core_cutoff_A'):
                    compare(lower[key],record[key])
                compare(lower['in_basin'],record['in_basin'][:1024])
                compare(lower['proposal_densities_rad_minus_d'],record['proposal_densities_rad_minus_d'][:1024])
                for a,b in zip(lower['energies_above_minimum_J_mol'],record['energies_above_minimum_J_mol'][:1024]):
                    if (a is None)!=(b is None):
                        raise AssertionError('same-stream finite-energy prefix differs')
                    if a is not None:
                        close(a,b,.1)
            source=SCRATCH/'xtb_checks'/name/f'{index:04d}'
            minimum=load(source/'result.json')
            if digest(source/'xtbopt.xyz')!=minimum['xyz_sha256']:
                raise AssertionError('minimum geometry hash mismatch')
            if not minimum['minimum_confirmed'] or min(minimum['frequencies_cm1'])<-1:
                raise AssertionError('unstable minimum included')
            stereo_check(SPECIES[name],read_xyz(source/'xtbopt.xyz')[0])
            audit_symmetry(source,minimum,SPECIES[name])
            legacy.verify_frequency_units(source,SPECIES[name],minimum)
            if record['input_sha256']!=digest(source/'xtbopt.xyz'):
                raise AssertionError('rotor geometry is stale')
            molecule=make_molecule(SPECIES[name],source/'xtbopt.xyz')
            definitions=rotor_definitions(molecule)
            if record['definitions']!=definitions or len(record['modes']['non_torsional_cm1'])!=3*molecule.GetNumAtoms()-6-len(definitions):
                raise AssertionError('rotor/vibration mode accounting mismatch')
            if not np.array_equal(record['coordinates_A'],molecule.GetConformer().GetPositions()):
                raise AssertionError('rotor coordinates differ from their minimum')
            compare(record['periods_rad'],[2*math.pi/r['symmetry'] for r in definitions])
            numbers=[float(v) for line in (source/'hessian').read_text().splitlines()
                     if not line.startswith('$') for v in line.split()]
            modes=projected_modes(molecule,definitions,np.array(numbers).reshape((3*molecule.GetNumAtoms(),)*2))
            for key in ('non_torsional_cm1','torsional_cm1','inertia_amu_A2','curvature_J_mol_rad2','masses_amu'):
                if not np.allclose(modes[key],record['modes'][key],rtol=1e-10,atol=1e-8):
                    raise AssertionError('projected modes differ from the source Hessian: '+key)
            # Independent Pitzer–Gwinn Gaussian-well harmonic limit.
            from scipy.linalg import eigvalsh
            from rotors import KB,HPLANCK,AMU_KG
            inertia=np.array(record['modes']['inertia_amu_A2'])*AMU_KG*1e-20
            curvature=np.array(record['modes']['curvature_J_mol_rad2'])
            d=len(definitions); temperature=500.
            x=HPLANCK*np.sqrt(eigvalsh(curvature/NA,inertia))/(2*math.pi*KB*temperature)
            gaussian=d/2*math.log(2*math.pi*R*temperature)-.5*np.linalg.slogdet(curvature)[1]
            classical=d/2*math.log(2*math.pi*KB*temperature)-d*math.log(HPLANCK)+.5*np.linalg.slogdet(inertia)[1]+gaussian
            pg=float(np.sum(np.log(x)-x/2-np.log(-np.expm1(-x))))
            oscillator=float(np.sum(-x/2-np.log(-np.expm1(-x))))
            close(classical+pg,oscillator,1e-8)
            points=legacy.verify_proposal(record,name,record.get('importance_proposal','diagonal'),ensemble.indices)
            centers=np.array(record['basin_centers_rad'])
            expected=np.array([wrap(angles(make_molecule(SPECIES[name],SCRATCH/'xtb_checks'/name/f'{i:04d}'/'xtbopt.xyz'),definitions),np.array(record['periods_rad'])) for i in ensemble.indices])
            if not np.allclose(centers,expected,rtol=0,atol=1e-10):
                raise AssertionError('basin centers mismatch')
            for level in ('pbe','blyp'):
                energy_path=SCRATCH/'composite'/name/f'{index:04d}'/level/'result.json'
                quantum=load(energy_path)
                if quantum['input_sha256']!=digest(source/'xtbopt.xyz'):
                    raise AssertionError('electronic geometry is stale')
                cursor=quantum
                cursor_min=minimum
                while 'derived_by_reflection_from' in cursor or 'derived_by_permutation_from' in cursor:
                    reflection='derived_by_reflection_from' in cursor
                    key='derived_by_reflection_from' if reflection else 'derived_by_permutation_from'
                    parent=cursor[key]
                    parent_source=SCRATCH/'xtb_checks'/name/f'{parent:04d}'
                    parent_path=SCRATCH/'composite'/name/f'{parent:04d}'/level/'result.json'
                    parent_energy=load(parent_path)
                    close(cursor['energy_hartree'],parent_energy['energy_hartree'])
                    if cursor['parent_result_sha256']!=digest(parent_path) or cursor_min['parent_hessian_sha256']!=digest(parent_source/'hessian'):
                        raise AssertionError('symmetry parent evidence stale')
                    cursor=parent_energy
                    cursor_min=load(parent_source/'result.json')
                if cursor.get('approximation'):
                    anchor_path=SCRATCH/'composite'/name/f"{cursor['anchor_index']:04d}"/level/'result.json'
                    if cursor['anchor_sha256']!=digest(anchor_path):
                        raise AssertionError('hybrid energy anchor stale')
                    if name.startswith('ps'):
                        if cursor['anchor_index']!=selected[0]:
                            raise AssertionError('fallback anchor is not lowest GFN2 chemical minimum')
                        close(cursor['offset_Eh'],load(anchor_path)['energy_hartree']-gfn[selected[0]],1e-10)
                    close(cursor['energy_hartree'],cursor_min['energy_Eh']+cursor['offset_Eh'],1e-10)
                else:
                    if name.startswith('ps'):
                        pools=cursor['threadpools']
                        if cursor['BLAS_threads']!=1 or any(p['num_threads']!=1 for p in pools if p['user_api']=='blas'):
                            raise AssertionError('new electronic result did not enforce one BLAS thread')
                        if not set(cursor['cpu_affinity'])<=set(CPUS) or cursor['native_OpenMP_threads']>len(cursor['cpu_affinity']):
                            raise AssertionError('new electronic result exceeded its CPU allocation')
                    path=SCRATCH/'composite'/name/f'{index:04d}'/level/'stdout.log'
                    # Exact symmetry images audit their real ancestor's log.
                    if quantum!=cursor:
                        parent=quantum.get('derived_by_reflection_from',quantum.get('derived_by_permutation_from'))
                        while True:
                            candidate=load(SCRATCH/'composite'/name/f'{parent:04d}'/level/'result.json')
                            next_parent=candidate.get('derived_by_reflection_from',candidate.get('derived_by_permutation_from'))
                            if next_parent is None:
                                break
                            parent=next_parent
                        path=SCRATCH/'composite'/name/f'{parent:04d}'/level/'stdout.log'
                    matches=re.findall(r'converged SCF energy\s*=\s*([-+\d.Ee]+)',path.read_text())
                    if not matches:
                        raise AssertionError('SCF source log missing')
                    close(float(matches[-1]),cursor['energy_hartree'],5e-9)
            count+=1
            points_count+=record['samples']
        for level in ('pbe','blyp'):
            for internal in (False,True):
                for t in (298.15,500.,800.):
                    row=ensemble.hs(t,level,internal)
                    plus=ensemble.hs(t+.01,level,internal)
                    minus=ensemble.hs(t-.01,level,internal)
                    dh=(plus['H_J_mol']-minus['H_J_mol'])/.02
                    ds=(plus['S_J_mol_K']-minus['S_J_mol_K'])/.02
                    close(dh,t*ds,0.05)
                    # Independent partition identity, including basin mixing.
                    p=np.array(row['probabilities'])
                    logparts=[]
                    for rec in ensemble.records:
                        item=thermal.old.conformer_thermal(rec,t,ensemble.symmetry)
                        logpart=item['log_q']
                        if internal:
                            logpart-=(item['S_translation_J_mol_K']+item['S_external_rotation_J_mol_K'])/R-4
                        logparts.append(logpart)
                    low=min(ensemble.energy[level])
                    logq=__import__('scipy').special.logsumexp(np.array(logparts)-(ensemble.energy[level]-low)/(R*t))
                    close(row['S_J_mol_K'],R*logq+(row['H_J_mol']-low)/t,1e-5)
    print('I048 source audit: %d stable wells, %d quadrature points; spectra, projected modes, stereo, symmetry orbit, energy logs, proposals and H/S identities verified'%(count,points_count),flush=True)


def fresh_checks(ensembles):
    from run_series import run_job
    if list((SCRATCH/'ensembles').glob('ps*/timeout_hold.json')):
        run_job([RMG_PYTHON,str(HERE/'hold_timeout.py'),'--self-check'],
                SCRATCH/'verification/timeout_protocol_cli',(CPUS[-1],),'fresh timeout-wrapper protocol')
    baseline_saved=load(SCRATCH/'baseline/baseline.json')
    run_job([RMG_PYTHON,str(HERE/'baseline.py')],SCRATCH/'verification/baseline',CPUS,'fresh pinned baseline')
    compare(baseline_saved,load(SCRATCH/'baseline/baseline.json'))
    for name in ('ethane','ps4_0000'):
        index=(load(SCRATCH/'composite'/name/'electronic_reduction.json')['DFT_indices'][0]
               if name.startswith('ps') else min(i for i in ensembles[name].indices if i<10000))
        source=SCRATCH/'xtb_checks'/name/f'{index:04d}'/'xtbopt.xyz'
        target=SCRATCH/'verification'/name/'pbe/result.json'
        run_job([DFT_PYTHON,str(HERE/'single_point.py'),'--xyz',str(source),'--output',str(target),'--level','pbe'],
                target.parent,CPUS[:4],'fresh DFT '+name)
        close(load(target)['energy_hartree'],load(SCRATCH/'composite'/name/f'{index:04d}'/'pbe/result.json')['energy_hartree'],.02/2625.499639)
    name='ps5_00000'
    ensemble=ensembles[name]
    record=ensemble.records[0]
    index=record['index']
    points=legacy.verify_proposal(record,name,'correlated',ensemble.indices)
    molecule=make_molecule(SPECIES[name],SCRATCH/'xtb_checks'/name/f'{index:04d}'/'xtbopt.xyz')
    definitions=record['definitions']
    initial=angles(molecule,definitions)
    center=np.array(record['basin_centers_rad'])[0]
    periods=np.array(record['periods_rad'])
    calculator=XTB(molecule,molecule.GetConformer().GetPositions())
    evidence=[]
    try:
        for i,(inside,value) in enumerate(zip(record['in_basin'],record['energies_above_minimum_J_mol'])):
            if not inside or value is None:
                continue
            xyz=rotated(molecule,definitions,initial+wrap(points[i]-center,periods))
            fresh=(calculator.energy(xyz)-record['reference_energy_Eh'])*EH_J*NA
            close(fresh,value,0.1)
            evidence.append({'sample':i,'saved_J_mol':value,'fresh_J_mol':fresh})
            if len(evidence)==8:
                break
    finally:
        calculator.close()
    if len(evidence)<8:
        raise AssertionError('not enough finite pentamer rotor points')
    save(SCRATCH/'verification/fresh_pentamer_rotor.json',{'name':name,'index':index,'points':evidence})
    print('I048 fresh checks: pinned baseline, ethane and tetramer PBE energies, eight pentamer torsional energies reproduced',flush=True)


def audit_cost():
    cost=load(SCRATCH/'cost_snapshot.json')
    guard_evidence=cost['timeout_guard_holds']+[cost['timeout_queue_restart']] if cost['timeout_queue_restart'] else cost['timeout_guard_holds']
    for evidence in guard_evidence:
        path=Path(evidence['path'])
        if digest(path)!=evidence['sha256']:
            raise AssertionError('timeout-guard execution evidence changed')
        compare(load(path),evidence['record'])
    if cost['workflow_queue_restart'] is not None:
        restart=SCRATCH/'workflow_queue_restart.json'
        if digest(restart)!=cost['workflow_queue_restart_sha256']:
            raise AssertionError('queue-restart evidence changed')
        compare(load(restart),cost['workflow_queue_restart'])
        interrupted=load(Path(cost['workflow_queue_restart']['preserved_interrupted_job']))
        if interrupted['exit_code']!=143 or not interrupted['label'].startswith('production stage early minima ps4_0110'):
            raise AssertionError('interrupted waiting-stage record differs')
    exception=SCRATCH/'development_thread_cap_exception.json'
    if digest(exception)!=cost['development_execution_exception_sha256']:
        raise AssertionError('execution-exception source changed')
    compare(load(exception),cost['development_execution_exception'])
    for violation in cost['captured_resource_violations']:
        if digest(Path(violation['path']))!=violation['sha256']:
            raise AssertionError('resource-violation source changed')
        compare(load(Path(violation['path'])),violation['record'])
    close(cost['elapsed_wall_s'],cost['cutoff_unix_s']-load(SCRATCH/'started.json')['unix_s'])
    total=0.
    peaks=[]
    for job in cost['recorded_jobs']:
        path=Path(job['path']); resources=Path(job['resources_path'])
        if digest(path)!=job['sha256'] or digest(resources)!=job['resources_sha256']:
            raise AssertionError('recorded cost source changed')
        record=load(path)
        close(job['elapsed_s'],record['elapsed_s'])
        if record['exit_code']!=job['exit_code']:
            raise AssertionError('job exit status differs')
        text=resources.read_text()
        def field(pattern):
            found=re.search(pattern,text)
            return float(found.group(1)) if found else 0.
        user=field(r'User time \(seconds\):\s*([\d.]+)')
        system=field(r'System time \(seconds\):\s*([\d.]+)')
        rss=field(r'Maximum resident set size \(kbytes\):\s*([\d.]+)')
        close(user,job['user_CPU_s']); close(system,job['system_CPU_s']); close(rss,job['max_RSS_kB'])
        peaks.append(rss)
        if not job['label'].startswith(('production stage ','replay ')):
            total+=user+system
    close(total,cost['CPU_s']); close(max(peaks or [0.]),cost['max_single_job_RSS_kB'])
    observations=[__import__('json').loads(line) for line in (SCRATCH/'resource_snapshots.jsonl').read_text().splitlines()]
    observations=[o for o in observations if o['unix_s']<=cost['last_resource_observation_unix_s']]
    if len(observations)!=cost['resource_observations']:
        raise AssertionError('cost resource observation count differs')
    close(max([o['RSS_kB'] for o in observations] or [0.]),cost['max_observed_aggregate_RSS_kB'])
    for o in observations:
        if sum(p['RSS_kB'] for p in o['processes'])!=o['RSS_kB'] or o['RSS_kB']>16*1024**2:
            raise AssertionError('resource sum or cap differs')
        if any(not set(mask)<=set(CPUS) for p in o['processes'] for mask in p['thread_affinities']):
            raise AssertionError('observed thread escaped core allocation')
    print('I048 recorded job cost, aggregate RSS observations and thread affinity reproduced',flush=True)


def main():
    invoked_at=time.time()
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay',action='store_true')
    parser.add_argument('--sources-only',action='store_true',help='development audit; does not satisfy dispatch Verifier')
    parser.add_argument('--available-only',action='store_true',help='audit and reproduce an explicitly incomplete report; does not satisfy the full dispatch Verifier')
    args=parser.parse_args()
    if load(SCRATCH/'method_plan_i048.json')!=PLAN:
        raise AssertionError('declared method changed')
    initial=load(SCRATCH/'method_plan_before_local_clarification.json')
    for key in ('candidate_window_kJ_mol','new_candidate_cap','energy','DFT_chemical_minima_per_sequence','DFT','rotors','new_rotor_counts','seed','prior_samples','prior_basin_overrides'):
        compare(initial[key],PLAN[key])
    for name,value in load(SCRATCH/'legacy_sources.json').items():
        if digest(HERE.parent/name)!=value:
            raise AssertionError('legacy implementation changed')
    for item in load(SCRATCH/'reuse.json'):
        if digest(Path(item['source']))!=item['sha256'] or digest(Path(item['target']))!=item['sha256']:
            raise AssertionError('read-only prior reuse changed')
    if load(SCRATCH/'sequences.json')!=SEQUENCES or [len(SEQUENCES[str(n)]) for n in range(2,6)]!=[2,3,6,10]:
        raise AssertionError('stereo enumeration mismatch')
    for n,classes in SEQUENCES.items():
        close(sum(row['weight'] for row in classes.values()),1.,1e-12)
        if sum(len(row['assignments']) for row in classes.values())!=2**int(n):
            raise AssertionError('oriented sequences missing')
    baseline=load(SCRATCH/'baseline/baseline.json')
    ensembles=thermal.make_ensembles(baseline,available_only=args.sources_only or args.available_only)
    audit_sources(ensembles)
    audit_partial_sources(ensembles)
    if args.sources_only:
        print('I048 DEVELOPMENT SOURCE AUDIT ONLY: full dispatch completion not verified',flush=True)
        return
    stem='increments_available_' if args.available_only else 'increments_'
    saved=load(SCRATCH/(stem+'2048.json'))
    if not args.available_only and (not saved['complete'] or set(saved['steps'])!={'2','3','4'}):
        raise AssertionError('required series incomplete')
    if set(saved['available_species'])!=set(ensembles):
        raise AssertionError('saved thermal coverage differs from available sources')
    if args.replay:
        # Completed replay must remain runnable after the production deadline.
        # The author's own replay stays within its original wall allocation;
        # a later manager replay gets a new, bounded verification invocation.
        production_end=load(SCRATCH/'started.json')['unix_s']+48*3600
        deadline=invoked_at+2*3600
        if invoked_at<production_end:
            deadline=min(deadline,production_end)
        if os.environ.get('I048_REPLAY_DEADLINE_UNIX'):
            deadline=min(deadline,float(os.environ['I048_REPLAY_DEADLINE_UNIX']))
        os.environ['I048_REPLAY_DEADLINE_UNIX']=str(deadline)
        from run_series import run_job
        for stage in ('search','minima','electronic','rotors'):
            command=[RMG_PYTHON,str(HERE/'run_series.py'),'--stage',stage]
            if args.available_only:
                command.extend(['--species',*sorted(name for name in ensembles if name.startswith('ps'))])
            run_job(command,SCRATCH/'verification'/stage,CPUS,'replay '+stage)
        fresh_checks(ensembles)
    reproduced=thermal.compute(available_only=args.available_only)
    compare(saved,reproduced)
    for n,levels in saved['steps'].items():
        for level in ('pbe','blyp'):
            row=thermal.evaluate_step(ensembles,baseline,int(n),500.,level)
            close(row['GAV_chain_orientation_shift_S_J_mol_K'],-R*math.log(2),1e-8)
            plus=thermal.evaluate_step(ensembles,baseline,int(n),500.01,level)
            minus=thermal.evaluate_step(ensembles,baseline,int(n),499.99,level)
            for kind in ('raw_delta','local_delta'):
                dh=(plus[kind]['H_J_mol']-minus[kind]['H_J_mol'])/.02
                ds=(plus[kind]['S_J_mol_K']-minus[kind]['S_J_mol_K'])/.02
                close(dh,500*ds,.1)
    lit=module('i048_literature',HERE/'literature.py')
    compare(lit.DATA,load(SCRATCH/'literature_crosscheck.json'))
    base=load(SCRATCH/(stem+'1024.json'))
    compare(base,thermal.compute(available_only=args.available_only,new_samples=1024,roots=False))
    audit_cost()
    report=module('i048_report',HERE/'report.py')
    document=report.REPORT.read_text()
    if document.split(report.START,1)[0]!=report.PREFACE or document.split(report.END,1)[1].strip():
        raise AssertionError('report method narrative or unverified suffix differs from the reproducible template')
    measured=document.split(report.START,1)[1].split(report.END,1)[0].strip()
    if measured!=report.render(saved,base,baseline):
        raise AssertionError('report numbers differ from replayed sources')
    print(('I048 INCOMPLETE REPORT ONLY: %d/%d species; full dispatch completion NOT verified'%(len(ensembles),len(SPECIES)) if args.available_only else 'I048 full series report and 1024/2048 quadrature comparison reproduced; no physical convergence or fresh-search identity asserted'),flush=True)


if __name__=='__main__':
    main()
