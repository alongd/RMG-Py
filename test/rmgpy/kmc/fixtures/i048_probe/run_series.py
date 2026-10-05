"""Prepare and execute the resource-bounded I048 calculation stages.

Commands (rmg_env, repository root, PYTHONPATH=$PWD):
python .../i048_probe/run_series.py --prepare
python .../i048_probe/run_series.py --stage search
python .../i048_probe/run_series.py --stage minima
python .../i048_probe/run_series.py --stage electronic
python .../i048_probe/run_series.py --stage rotors
Each subprocess persists both streams and its timing/resource record.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from threading import Lock
from contextlib import contextmanager,ExitStack
import datetime
import fcntl
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time

from common import (SCRATCH, PRIOR, HERE, LEGACY, RMG_PYTHON, DFT_PYTHON,
                    CPUS, PLAN, bootstrap, numerical_environment, save, load, digest)

SPECIES, SEQUENCES = bootstrap()


class StageBusy(Exception):
    """Another producer owns this case; no calculation has been started."""


def prepare():
    from rdkit import Chem
    from rdkit.Chem import AllChem
    path = SCRATCH/'method_plan_i048.json'
    if path.exists() and load(path) != PLAN:
        raise ValueError('method plan changed; preserve evidence and explicitly revise before increments')
    save(path, PLAN)
    if not (SCRATCH/'started.json').exists():
        save(SCRATCH/'started.json', {'unix_s':time.time(), 'UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()})
    manifest = {}
    for name,smiles in SPECIES.items():
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        params = AllChem.ETKDGv3()
        params.randomSeed = 43002
        params.numThreads = 1
        if AllChem.EmbedMolecule(mol,params) or AllChem.MMFFOptimizeMolecule(mol,maxIters=2000):
            raise RuntimeError('starting embedding/optimization failed: '+name)
        xyz = Chem.MolToXYZBlock(mol)
        target = SCRATCH/'ensembles'/name/'input.xyz'
        target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists() and target.read_text()!=xyz:
            raise AssertionError('starting input changed: '+name)
        target.write_text(xyz)
        mirror = Chem.Mol(mol)
        for a in mirror.GetAtoms():
            if a.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED:
                a.InvertChirality()
        canon=Chem.MolToSmiles(Chem.RemoveHs(mol))
        mirrored=Chem.MolToSmiles(Chem.RemoveHs(mirror))
        manifest[name]={'smiles':smiles,'canonical_smiles':canon,'mirror_smiles':mirrored,
                        'achiral':canon==mirrored,'natoms':mol.GetNumAtoms(),
                        'input_sha256':digest(target)}
    save(SCRATCH/'species.json', manifest)
    save(SCRATCH/'sequences.json', SEQUENCES)
    sources = []
    for name in PLAN['prior_samples']:
        for stage in ('ensembles','xtb_checks','composite','rotors'):
            source=PRIOR/stage/name
            target=SCRATCH/stage/name
            if source.exists():
                if stage=='ensembles' and digest(source/'input.xyz')!=digest(target/'input.xyz'):
                    raise AssertionError('reuse input mismatch: '+name)
                if not (SCRATCH/'reuse.json').exists():
                    shutil.copytree(source,target,dirs_exist_ok=True)
                for f in sorted(source.rglob('*')):
                    if f.is_file():
                        sources.append({'source':str(f),'target':str(target/f.relative_to(source)),
                                        'sha256':digest(f)})
    if (SCRATCH/'reuse.json').exists():
        if load(SCRATCH/'reuse.json')!=sources:
            raise AssertionError('prior evidence changed')
    else:
        save(SCRATCH/'reuse.json',sources)
    save(SCRATCH/'legacy_sources.json', {str(p.relative_to(HERE.parent)):digest(p) for p in sorted(LEGACY.glob('*.py'))})
    print('I048 method fixed before increments; 2/3/6/10 stereo classes; prior molecular evidence copied read-only',flush=True)


def run_job(command, directory, cpus, label, timeout_s=None):
    directory.mkdir(parents=True,exist_ok=True)
    env=numerical_environment()
    affinity=','.join(map(str,cpus))
    if os.environ.get('I048_REPLAY_DEADLINE_UNIX'):
        remaining=float(os.environ['I048_REPLAY_DEADLINE_UNIX'])-time.time()
    else:
        remaining=48*3600-(time.time()-load(SCRATCH/'started.json')['unix_s'])
    if remaining<=31:
        raise TimeoutError('48 hour total wall budget exhausted')
    # Reserve the wrapper's termination grace inside the original allocation.
    timeout_s=min(timeout_s or int(remaining),int(remaining)-31)
    run_id=str(int(time.time()*1e6))
    resources=directory/('resources_'+run_id+'.txt')
    wrapped=['timeout','--kill-after=30',str(timeout_s),'taskset','-c',affinity,
             '/usr/bin/time','-v','-o',str(resources),*command]
    shell=shlex.join(wrapped)+' > >(taskset -c '+affinity+' tee -a stdout.log) 2> >(taskset -c '+affinity+' tee -a stderr.log >&2)'
    started=time.time()
    print('I048 START '+label+': '+shlex.join(command),flush=True)
    result=subprocess.run(['taskset','-c',affinity,'bash','-o','pipefail','-c',shell],cwd=directory,env=env)
    record={'label':label,'command':command,'cpus':list(cpus),'BLAS_threads':1,
            'started_unix_s':started,'elapsed_s':time.time()-started,'exit_code':result.returncode,
            'resources_path':str(resources)}
    save(directory/('job_'+run_id+'.json'),record)
    save(directory/'job.json',record)
    if result.returncode:
        raise RuntimeError('job failed: '+label+'; exit '+str(result.returncode))
    print('I048 DONE '+label+' %.1f s'%record['elapsed_s'],flush=True)
    return record


def search_succeeded(name,record):
    """Distinguish a preserved timeout-wrapper status from the CREST exit."""
    if record['exit_code']==0:
        return True
    directory=SCRATCH/'ensembles'/name
    receipt=directory/'timeout_hold.json'
    if record['exit_code']!=124 or not receipt.exists():
        return False
    proof=load(receipt)
    resources=Path(proof['resources_path'])
    deadline=load(SCRATCH/'started.json')['unix_s']+48*3600
    expected=['timeout','--kill-after=30','28800','taskset','-c',
              ','.join(map(str,record['cpus'])),'/usr/bin/time','-v','-o',str(resources),*record['command']]
    return (proof.get('scientific_completed') is True
            and proof['science_changes'] is False
            and proof['command']==record['command']
            and proof['wrapper_command']==expected
            and resources.parent==directory
            and proof['production_deadline_unix_s']==deadline
            and proof['child_completed_unix_s']<deadline
            and proof['stopped_unix_s']<proof['original_timeout_unix_s']
            and digest(resources)==proof['resources_sha256']
            and re.search(r'Exit status:\s*0\s*$',resources.read_text()) is not None
            and digest(directory/'stdout.log')==proof['stdout_sha256']
            and 'CREST terminated normally.' in (directory/'stdout.log').read_text())


def finalize_search(name,record):
    from run_ensembles import CREST,read_xyz
    directory=SCRATCH/'ensembles'/name
    if not search_succeeded(name,record):
        raise RuntimeError('CREST did not complete successfully: '+name)
    frames=read_xyz(directory/'crest_conformers.xyz')
    energies=[float(f['comment'].split()[0]) for f in frames]
    selected=[i for i,e in enumerate(energies) if (e-min(energies))*2625.499639<12.]
    save(directory/'completed.json',{'command':record['command'],'elapsed_s':record['elapsed_s'],'frames':len(frames),
                 'energy_hartree':energies,'selected_indices_within_12_kJ':selected,
                 'ensemble_sha256':digest(directory/'crest_conformers.xyz')})


def search(name,cpus):
    from run_ensembles import CREST
    directory=SCRATCH/'ensembles'/name
    if (directory/'completed.json').exists():
        return
    env_path=Path(CREST).parent
    command=[str(env_path/'crest'),'input.xyz','--gfn2','--quick','--ewin','6','--T',str(len(cpus))]
    previous=os.environ.get('PATH','')
    os.environ['PATH']=str(env_path)+':'+previous if not previous.startswith(str(env_path)+':') else previous
    # The dispatch limits total wall time, rather than eight hours per search.
    record=run_job(command,directory,cpus,'search '+name)
    finalize_search(name,record)


@contextmanager
def case_lock(name,stage,wait=True):
    directory=SCRATCH/'composite'/name
    directory.mkdir(parents=True,exist_ok=True)
    with (directory/(stage+'.lock')).open('a') as lock:
        try:
            fcntl.flock(lock,fcntl.LOCK_EX|(0 if wait else fcntl.LOCK_NB))
        except BlockingIOError:
            raise StageBusy(name+' '+stage) from None
        yield


def minima(name,cpu,wait=True):
    with case_lock(name,'minima',wait):
        _minima(name,cpu)


def _minima(name,cpu):
    directory=SCRATCH/'composite'/name
    if (directory/'selection.json').exists():
        return
    rec=load(SCRATCH/'ensembles'/name/'completed.json')
    chosen=sorted(rec['selected_indices_within_12_kJ'],key=lambda i:rec['energy_hartree'][i])[:PLAN['new_candidate_cap']]
    save(directory/'candidate_reduction.json',{'all_within_12_kJ':len(rec['selected_indices_within_12_kJ']),
        'retained_indices':chosen,'omitted':len(rec['selected_indices_within_12_kJ'])-len(chosen),
        'rule':'first 32 by GFN2 energy, fixed before increments'})
    cpus=(cpu,) if isinstance(cpu,int) else tuple(cpu)
    batches=[chosen[i::len(cpus)] for i in range(len(cpus))]
    save(directory/'checks/execution.json',{'lane_CPUs':list(cpus),'candidate_batches':batches,
         'reason':'independent candidate checks; unchanged candidate set and one native/BLAS thread per xTB process'})
    def check_lane(i):
        if not batches[i]:
            return
        target=directory/'checks' if len(cpus)==1 else directory/'checks'/('lane_'+str(i))
        run_job([RMG_PYTHON,str(HERE/'legacy_cli.py'),'check_xtb','--species',name,'--indices',*map(str,batches[i])],
                target,(cpus[i],),'minima '+name+' lane '+str(i))
    with ThreadPoolExecutor(max_workers=len(cpus)) as pool:
        for result in pool.map(check_lane,range(len(cpus))):
            pass
    # Selection and symmetry closure run in a separate single-core producer;
    # avoiding global json monkeypatches across supervisor threads.
    run_job([RMG_PYTHON,str(HERE/'select_pool.py'),name],directory/'select',(cpus[0],),'selection '+name)


def electronic(name,cpus,wait=True):
    with case_lock(name,'electronic',wait):
        _electronic(name,cpus)


def _electronic(name,cpus):
    inputs=electronic_inputs(name)
    if not missing_electronic(name,inputs):
        return
    for i in inputs[-1]:
        for level in ('pbe','blyp'):
            direct_electronic(name,i,level,cpus)
    finish_electronic(name,inputs)


def electronic_inputs(name):
    pool=load(SCRATCH/'composite'/name/'selection.json')
    chemical=[i for i in pool['unique_minima_indices'] if i<10000]
    energies={i:load(SCRATCH/'xtb_checks'/name/f'{i:04d}'/'result.json')['energy_Eh'] for i in chemical}
    selected=sorted(chemical,key=energies.get)[:PLAN['DFT_chemical_minima_per_sequence']]
    save(SCRATCH/'composite'/name/'electronic_reduction.json',{'DFT_indices':selected,'GFN2_indices':sorted(set(chemical)-set(selected))})
    return pool,chemical,energies,selected


def direct_electronic(name,i,level,cpus):
    directory=SCRATCH/'composite'/name/f'{i:04d}'/level
    if not (directory/'result.json').exists():
        run_job([DFT_PYTHON,str(HERE/'single_point.py'),'--xyz',str(SCRATCH/'xtb_checks'/name/f'{i:04d}'/'xtbopt.xyz'),
                 '--output',str(directory/'result.json'),'--level',level],directory,cpus,'DFT '+name+' '+str(i)+' '+level)


def missing_electronic(name,inputs):
    return any(not (SCRATCH/'composite'/name/f'{i:04d}'/level/'result.json').exists()
               for i in inputs[0]['unique_minima_indices'] for level in ('pbe','blyp'))


def finish_electronic(name,inputs):
    from complete_symmetry import materialize_energy
    pool,chemical,energies,selected=inputs
    for level in ('pbe','blyp'):
        anchor=selected[0]
        anchor_path=SCRATCH/'composite'/name/f'{anchor:04d}'/level/'result.json'
        offset=load(anchor_path)['energy_hartree']-energies[anchor]
        for i in chemical:
            if i in selected:
                continue
            source=SCRATCH/'xtb_checks'/name/f'{i:04d}'/'xtbopt.xyz'
            save(SCRATCH/'composite'/name/f'{i:04d}'/level/'result.json',
                 {'energy_hartree':energies[i]+offset,'input_sha256':digest(source),
                  'approximation':'GFN2 plus functional offset from lowest GFN2 conformer',
                  'anchor_index':anchor,'anchor_sha256':digest(anchor_path),'offset_Eh':offset})
    for i in pool['unique_minima_indices']:
        for level in ('pbe','blyp'):
            materialize_energy(SCRATCH,name,i,level)


def bulk_electronic(names,slots):
    """Use the same three lanes for independent points, including one case."""
    with ExitStack() as locks:
        for name in names:
            locks.enter_context(case_lock(name,'electronic'))
        inputs={name:electronic_inputs(name) for name in names}
        needs_finish={name for name in names if missing_electronic(name,inputs[name])}
        tasks=[(name,i,level) for name in names for i in inputs[name][-1]
               for level in ('pbe','blyp')
               if not (SCRATCH/'composite'/name/f'{i:04d}'/level/'result.json').exists()]
        if tasks:
            execution={'lane_CPUs':[list(cpus) for cpus in slots],
                       'tasks':[{'name':name,'index':i,'level':level,'lane':j%len(slots)}
                                for j,(name,i,level) in enumerate(tasks)],
                       'science_changes':False,'unix_s':time.time()}
            save(SCRATCH/'pipeline'/('bulk_electronic_'+str(int(time.time()*1e6))+'.json'),execution)
            save(SCRATCH/'pipeline/bulk_electronic_execution.json',execution)
        def lane(j):
            for name,i,level in tasks[j::len(slots)]:
                direct_electronic(name,i,level,slots[j])
        with ThreadPoolExecutor(max_workers=len(slots)) as pool:
            for _ in pool.map(lane,range(len(slots))):
                pass
        for name in names:
            if name in needs_finish:
                finish_electronic(name,inputs[name])
    return len(tasks)


def rotor(name,cpu):
    directory=SCRATCH/'rotor_jobs'/name
    directory.mkdir(parents=True,exist_ok=True)
    # A prepared case may be integrated before the bulk stage. The same
    # process lock also protects a concurrent replay from duplicate writers.
    with (directory/'producer.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        for samples in PLAN['new_rotor_counts']:
            run_job([RMG_PYTHON,str(HERE/'legacy_cli.py'),'rotors','--species',name,
                     '--samples',str(samples),'--proposal','correlated'],
                    directory/str(samples),(cpu,),'rotors '+name+' '+str(samples))


def bulk_rotors(names,cpus):
    """Spread independent frozen basins over the same eight single-core lanes."""
    from rotors import quadrature_filename
    with ExitStack() as locks:
        tasks=[]
        pending=[]
        for name in names:
            directory=SCRATCH/'rotor_jobs'/name
            directory.mkdir(parents=True,exist_ok=True)
            lock=locks.enter_context((directory/'producer.lock').open('a'))
            fcntl.flock(lock,fcntl.LOCK_EX)
            indices=load(SCRATCH/'composite'/name/'selection.json')['unique_minima_indices']
            missing=[(name,i,count) for i in indices for count in PLAN['new_rotor_counts']
                     if not (SCRATCH/'rotors'/name/f'{i:04d}'/quadrature_filename(count,PLAN['seed'],'correlated')).exists()]
            if missing:
                run_job([RMG_PYTHON,str(HERE/'rotor_worker.py'),name,'--seal'],
                        directory/'partition',(cpus[0],),'rotor partition '+name)
                # Closure must precede enumerating tasks, including its images.
                indices=load(SCRATCH/'composite'/name/'selection.json')['unique_minima_indices']
                tasks.extend((name,i,count) for i in indices for count in PLAN['new_rotor_counts']
                             if not (SCRATCH/'rotors'/name/f'{i:04d}'/quadrature_filename(count,PLAN['seed'],'correlated')).exists())
                pending.append(name)
        if tasks:
            execution={'lane_CPUs':list(cpus),'science_changes':False,'unix_s':time.time(),
                       'tasks':[{'name':name,'index':i,'samples':count,'lane':None,
                                 'selection_sha256':digest(SCRATCH/'composite'/name/'selection.json')}
                                for name,i,count in tasks]}
            archive=SCRATCH/'pipeline'/('bulk_rotors_'+str(int(time.time()*1e6))+'.json')
            save(archive,execution)
            save(SCRATCH/'pipeline/bulk_rotors_execution.json',execution)
        next_task=[0]
        assignment_lock=Lock()
        def lane(j):
            while True:
                # A lane that finishes a basin can take the next fixed task.
                # Samples and seeds depend on the basin, never on its lane.
                with assignment_lock:
                    position=next_task[0]
                    if position==len(tasks):
                        return
                    next_task[0]+=1
                    execution['tasks'][position]['lane']=j
                    save(archive,execution)
                    save(SCRATCH/'pipeline/bulk_rotors_execution.json',execution)
                name,i,count=tasks[position]
                run_job([RMG_PYTHON,str(HERE/'rotor_worker.py'),name,'--index',str(i),'--samples',str(count)],
                        SCRATCH/'rotor_jobs'/name/'basins'/str(count)/f'{i:04d}',(cpus[j],),
                        'rotor basin '+name+' '+str(i)+' '+str(count))
        with ThreadPoolExecutor(max_workers=len(cpus)) as pool:
            for _ in pool.map(lane,range(len(cpus))):
                pass
        # Canonical case CLIs verify their entire saved pools and preserve the
        # existing positive per-count job receipts. Cached data are not rewritten.
        for name in names:
            for count in PLAN['new_rotor_counts']:
                directory=SCRATCH/'rotor_jobs'/name/str(count)
                if name not in pending and (directory/'job.json').exists() and load(directory/'job.json')['exit_code']==0:
                    continue
                run_job([RMG_PYTHON,str(HERE/'legacy_cli.py'),'rotors','--species',name,
                         '--samples',str(count),'--proposal','correlated'],directory,(cpus[0],),
                        'rotors '+name+' '+str(count))
    return len(tasks)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare',action='store_true')
    parser.add_argument('--stage',choices=('search','minima','electronic','rotors'))
    parser.add_argument('--species',nargs='+',choices=tuple(name for name in SPECIES if name.startswith('ps')))
    parser.add_argument('--search-lane',type=int,choices=(0,1),help='existing four-core lane for one search')
    parser.add_argument('--if-idle',action='store_true',help='one preparation case: return 75 if owned by another producer')
    args=parser.parse_args()
    if args.if_idle and (args.stage not in ('minima','electronic') or not args.species or len(args.species)!=1):
        parser.error('--if-idle requires one minimum/electronic preparation case')
    if args.search_lane is not None and (args.stage!='search' or not args.species or len(args.species)!=1):
        parser.error('--search-lane requires exactly one search case')
    os.sched_setaffinity(0,CPUS)
    if args.prepare:
        prepare()
        return
    if load(SCRATCH/'method_plan_i048.json')!=PLAN:
        raise AssertionError('method declaration missing/changed')
    names=args.species or [name for name in SPECIES if name.startswith('ps')]
    if args.stage=='search':
        workers=2
        fn=search
        slots=[CPUS[:4],CPUS[4:]]
        if args.search_lane is not None:
            slots=[slots[args.search_lane]]
            workers=1
    elif args.stage=='electronic':
        slots=[tuple(slot) for slot in PLAN['DFT_execution_lanes']]
        workers=len(slots)
        fn=electronic
    elif args.stage=='minima':
        names=[name for name in names if not (SCRATCH/'composite'/name/'selection.json').exists()]
        if not names:
            print('I048 stage complete: minima (cached)',flush=True)
            return
        workers=min(len(CPUS),len(names))
        slots=[tuple(CPUS[i::workers]) for i in range(workers)]
        if args.species:
            # Early preparation shares the allocation with ongoing searches.
            # Bulk preparation can use all eight lanes after searches finish.
            slots=[slot[:2] for slot in slots]
        fn=minima
    elif args.stage=='rotors':
        workers=8
        fn=rotor
        slots=list(CPUS)
    else:
        parser.error('choose --prepare or --stage')
    if args.if_idle:
        from functools import partial
        fn=partial(fn,wait=False)
    # A fixed lane owns each affinity slot, including its descendants/loggers.
    def lane(index):
        for name in names[index::workers]:
            fn(name,slots[index])
    (SCRATCH/'pipeline').mkdir(parents=True,exist_ok=True)
    with (SCRATCH/'pipeline/supplemental.lock').open('a') as extra_lock:
        if args.stage=='electronic' and not args.species:
            print('I048 bulk electronic stage waits for supplemental producer',flush=True)
            fcntl.flock(extra_lock,fcntl.LOCK_SH)
            points=bulk_electronic(names,slots)
            print('I048 stage complete: electronic (%d independent single points)'%points,flush=True)
            return
        with (SCRATCH/'pipeline/early_rotors.lock').open('a') as rotor_lock:
            if args.stage=='rotors' and not args.species:
                print('I048 bulk rotor stage waits for early rotor producer',flush=True)
                fcntl.flock(rotor_lock,fcntl.LOCK_SH)
                points=bulk_rotors(names,CPUS)
                print('I048 stage complete: rotors (%d independent basin integrations)'%points,flush=True)
                return
            with ThreadPoolExecutor(max_workers=workers) as pool:
                for result in pool.map(lane,range(workers)):
                    pass
    print('I048 stage complete: '+args.stage,flush=True)


if __name__=='__main__':
    try:
        main()
    except StageBusy as error:
        print('I048 preparation already owned: '+str(error),flush=True)
        sys.exit(75)
