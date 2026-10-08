#!/home/alon/anaconda3/envs/rmg_env/bin/python
"""Authoring only: add independent floor cohorts after a qualified Rouse run.

The normal oracle and Verifier always generate every trajectory afresh. This
driver preserves actual old input/program hashes instead of relabelling them.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
if hasattr(os,'sched_getaffinity'):
    os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:8]))

import argparse
from concurrent.futures import ProcessPoolExecutor
import importlib.util
import json
import multiprocessing
from pathlib import Path
import resource
import sys
import time

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
import checks
spec=importlib.util.spec_from_file_location('floor_extension_oracle',ROOT/'reference/run.py')
oracle=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=oracle
spec.loader.exec_module(oracle)


def versions():
    return {'python':oracle.platform.python_version(),'numpy':oracle.np.__version__,
            'scipy':oracle.scipy.__version__,'rng':'numpy.PCG64',
            'platform':oracle.platform.system(),'machine':oracle.platform.machine()}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--previous',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workers',type=int,choices=range(1,9),default=8)
    args=parser.parse_args()
    begin=time.monotonic()
    cfg,digest=oracle.read_parameters()
    old_bytes=(ROOT/'reference/provenance/parameters_04.json').read_bytes()
    old_source=(ROOT/'reference/provenance/measurement_04.py').read_bytes()
    old=json.loads(old_bytes)
    source=(ROOT/'reference/run.py').read_bytes()
    driver=Path(__file__).read_bytes()
    previous=json.loads(args.previous.read_text())
    assert previous['parameters_sha256']==oracle.fingerprint(old_bytes)
    assert previous['program_sha256']==oracle.fingerprint(old_source)
    assert previous['versions']==versions()
    checks.assert_generation_inputs_equivalent(old,cfg,old_source,source)
    assert len(previous['stage_simulation'])==96*old['precision_stages']
    assert all(row['qualified'] for row in previous['qualification'] if row['case']!='floor')
    args.output.mkdir(parents=True,exist_ok=False)
    checkpoint=args.output/'ensembles';checkpoint.mkdir()
    (args.output/'parsed_inputs.json').write_text(json.dumps({
        'parameters_sha256':digest,'program_sha256':oracle.fingerprint(source),
        'driver_program_sha256':oracle.fingerprint(driver),'parameters':cfg},indent=2)+'\n')
    case_index,floor=next((i,case) for i,case in enumerate(cfg['cases']) if case['name']=='floor')
    mapped=oracle.mapped_parameters(cfg,floor)
    jobs=[(stage,run,cfg,mapped,'end/end',run['seed']+100*case_index+100000*stage)
          for stage in range(old['precision_stages'],floor['precision_stages'])
          for run in cfg['runs']]
    print('AUTHORING FLOOR EXTENSION: archived input/source hashes, versions, physical inputs and generator ASTs PASS.',flush=True)
    rows=list(previous['stage_simulation'])
    with ProcessPoolExecutor(max_workers=args.workers,mp_context=multiprocessing.get_context('spawn')) as pool:
        futures=[(job[0],pool.submit(oracle.run_one,*job[1:])) for job in jobs]
        for index,(stage,future) in enumerate(futures):
            row=dict(future.result(),precision_stage=stage)
            rows.append(row)
            pending=checkpoint/f'{index:03d}.pending'
            pending.write_text(json.dumps(row,indent=2,allow_nan=False)+'\n')
            pending.rename(checkpoint/f'{index:03d}.json')
    result=oracle.result_from_stages(rows,cfg)
    for old_row in previous['qualification']:
        if old_row['case']!='floor':
            assert old_row==next(row for row in result['qualification']
                                if (row['case'],row['pair_class'])==(old_row['case'],old_row['pair_class']))
    result['parameters_sha256']=digest
    result['program_sha256']=oracle.fingerprint(source)
    result['versions']=versions()
    result['generation_provenance']=[
        {'precision_stage':stage,'parameters_sha256':previous['parameters_sha256'],
         'program_sha256':previous['program_sha256'],
         'driver_program_sha256':previous['program_sha256'],
         'cases':[case['name'] for case in old['cases']]}
        for stage in range(old['precision_stages'])]
    result['generation_provenance'] += [
        {'precision_stage':stage,'parameters_sha256':digest,
         'program_sha256':oracle.fingerprint(source),
         'driver_program_sha256':oracle.fingerprint(driver),'cases':['floor']}
        for stage in range(old['precision_stages'],floor['precision_stages'])]
    assert (ROOT/'reference/run.py').read_bytes()==source
    assert Path(__file__).read_bytes()==driver
    checks.assert_generation_history(cfg,result,versions())
    (args.output/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    (args.output/'results.md').write_text(oracle.markdown(result))
    children=resource.getrusage(resource.RUSAGE_CHILDREN)
    own=resource.getrusage(resource.RUSAGE_SELF)
    runtime={'wall_seconds':time.monotonic()-begin,'child_cpu_seconds':children.ru_utime+children.ru_stime,
             'parent_cpu_seconds':own.ru_utime+own.ru_stime,'child_max_rss_MiB':children.ru_maxrss/1024,
             'parent_max_rss_MiB':own.ru_maxrss/1024,'workers':args.workers,'threads_per_worker':1,
             'cpu_affinity':sorted(os.sched_getaffinity(0)),'new_physical_ensembles':len(jobs)}
    (args.output/'runtime.json').write_text(json.dumps(runtime,indent=2)+'\n')
    print('FLOOR EXTENSION COMPUTED:',json.dumps(runtime),flush=True)


if __name__=='__main__':
    main()
