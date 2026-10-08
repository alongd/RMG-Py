"""Integrate prepared cases during searches within the existing core set.

Command: rmg_env python .../i048_probe/early_rotors.py
Command: same, --basin-parallel (independent basins on the same eight cores)
Both declared sample counts are retained; the bulk rotor stage waits on a lock.
"""
import argparse
import fcntl
import os
import time

from common import SCRATCH,CPUS,bootstrap,load,save
from run_series import rotor,bulk_rotors


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--basin-parallel',action='store_true')
    args=parser.parse_args()
    os.sched_setaffinity(0,CPUS)
    species,_=bootstrap()
    names=[name for name in species if name.startswith('ps')]
    (SCRATCH/'pipeline').mkdir(parents=True,exist_ok=True)
    with (SCRATCH/'pipeline/early_rotors.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        execution_path=SCRATCH/'early_rotor_execution.json'
        if execution_path.exists():
            save(SCRATCH/('early_rotor_execution_'+str(int(time.time()*1e6))+'.json'),load(execution_path))
        save(execution_path,{
            'started_unix_s':time.time(),'physical_core_set':list(CPUS),
            'rotor_core':12,'science_changes':False,
            'basin_parallel':args.basin_parallel,
            'basin_lane_CPUs':list(CPUS) if args.basin_parallel else [12],
            'reason':'first prepared tetramer required 1230.2+2356.0 seconds for the two declared quadratures; spread this work during searches',
            'bulk_rotors_wait_on_this_process_lock':True})
        while True:
            if all((SCRATCH/'ensembles'/name/'completed.json').exists() for name in names):
                break
            if time.time()-load(SCRATCH/'started.json')['unix_s']>=48*3600:
                raise TimeoutError('total wall budget exhausted')
            ready=[]
            for name in names:
                prepared=any((SCRATCH/'pipeline'/stage/name/'completed.json').exists()
                             for stage in ('early','supplemental'))
                result=SCRATCH/'rotor_jobs'/name/'2048/job.json'
                if prepared and not result.exists():
                    ready.append(name)
                elif prepared and load(result)['exit_code']:
                    raise RuntimeError('rotor integration failed; see '+str(result))
            if not ready:
                time.sleep(15)
                continue
            print('I048 early rotor preparation '+ready[0],flush=True)
            if args.basin_parallel:
                bulk_rotors([ready[0]],CPUS)
            else:
                rotor(ready[0],12)
        print('I048 early rotor producer finished; bulk rotor lock released',flush=True)


if __name__=='__main__':
    main()
