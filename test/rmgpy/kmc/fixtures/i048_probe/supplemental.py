"""Prepare completed cases on a second lane during ongoing searches.

Command: rmg_env python .../i048_probe/supplemental.py
Same eight physical cores and declared chemistry as the primary pipeline.
"""
import fcntl
import os
import time

from common import SCRATCH,CPUS,bootstrap,load,save
from run_series import minima,electronic,search_succeeded


def main():
    os.sched_setaffinity(0,CPUS)
    species,_=bootstrap()
    names=[name for name in species if name.startswith('ps')]
    (SCRATCH/'pipeline').mkdir(parents=True,exist_ok=True)
    with (SCRATCH/'pipeline/supplemental.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        save(SCRATCH/'supplemental_execution.json',{
            'started_unix_s':time.time(),'physical_core_set':list(CPUS),
            'minima_core':10,'DFT_cores':[6,8,10],
            'excluded_case':'ps4_0000: already active before case locks were installed',
            'reason':'completed four-core searches averaged 2.0555 and 2.3170 CPU cores; use a second preparation lane within the same physical allocation',
            'science_changes':False,'bulk_DFT_waits_on_this_process_lock':True})
        while True:
            complete=[name for name in names if (SCRATCH/'ensembles'/name/'completed.json').exists()]
            if len(complete)==len(names):
                break
            if time.time()-load(SCRATCH/'started.json')['unix_s']>=48*3600:
                raise TimeoutError('total wall budget exhausted')
            for name in names:
                path=SCRATCH/'ensembles'/name/'job.json'
                if path.exists() and not search_succeeded(name,load(path)):
                    raise RuntimeError('search failed; see '+str(path))
            ready=[name for name in reversed(names) if name in complete and name!='ps4_0000'
                   and not (SCRATCH/'pipeline/early'/name/'completed.json').exists()
                   and not (SCRATCH/'pipeline/supplemental'/name/'completed.json').exists()]
            if not ready:
                time.sleep(15)
                continue
            name=ready[0]
            print('I048 supplemental preparation '+name,flush=True)
            minima(name,10)
            electronic(name,(6,8,10))
            save(SCRATCH/'pipeline/supplemental'/name/'completed.json',{
                'name':name,'unix_s':time.time(),'stages':['minima','electronic']})
        print('I048 supplemental producer finished; bulk DFT lock released',flush=True)


if __name__=='__main__':
    main()
