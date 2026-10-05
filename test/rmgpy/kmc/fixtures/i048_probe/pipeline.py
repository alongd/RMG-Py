"""Finish the declared producers after all searches have completed.

Command: rmg_env python .../i048_probe/pipeline.py --after-search
Command: same, --overlap-electronic (prepare completed cases during searches)
Does not author conclusions, calculate Tc, or claim Verifier completion.
"""
import argparse
import fcntl
import time
from common import SCRATCH,HERE,CPUS,RMG_PYTHON,bootstrap,load,save
from run_series import run_job,search_succeeded

species,_=bootstrap()
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--after-search',action='store_true',required=True)
parser.add_argument('--overlap-electronic',action='store_true')
args=parser.parse_args()
names=[name for name in species if name.startswith('ps')]


def preparing_elsewhere(name):
    directory=SCRATCH/'composite'/name
    directory.mkdir(parents=True,exist_ok=True)
    for stage in ('minima','electronic'):
        with (directory/(stage+'.lock')).open('a') as lock:
            try:
                fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            except BlockingIOError:
                return True
    return False


last=None
while True:
    complete=[name for name in names if (SCRATCH/'ensembles'/name/'completed.json').exists()]
    if complete!=last:
        print('I048 pipeline: %d/%d new searches complete'%(len(complete),len(names)),flush=True)
        last=complete
    if len(complete)==len(names):
        break
    for name in names:
        path=SCRATCH/'ensembles'/name/'job.json'
        if path.exists() and not search_succeeded(name,load(path)):
            raise RuntimeError('search failed; see '+str(path))
    if time.time()-load(SCRATCH/'started.json')['unix_s']>=48*3600:
        raise TimeoutError('total wall budget exhausted while waiting for searches')
    if args.overlap_electronic:
        ready=[name for name in names if name in complete
               and not (SCRATCH/'pipeline/early'/name/'completed.json').exists()
               and not preparing_elsewhere(name)]
        if ready:
            # Prepare one idle case at a time, allowing another producer to
            # retain its case while this lane handles other available work.
            name=ready[0]
            finished=True
            for stage in ('minima','electronic'):
                directory=SCRATCH/'pipeline/early'/name/stage
                try:
                    run_job([RMG_PYTHON,str(HERE/'run_series.py'),'--stage',stage,'--species',name,'--if-idle'],
                            directory,CPUS,'production stage early '+stage+' '+name)
                except RuntimeError:
                    if load(directory/'job.json')['exit_code']!=75:
                        raise
                    print('I048 preparation owned elsewhere; retaining pending case '+name,flush=True)
                    finished=False
                    break
            if not finished:
                continue
            save(SCRATCH/'pipeline/early'/name/'completed.json',{'name':name,'unix_s':time.time(),'stages':['minima','electronic']})
            continue
    time.sleep(15)
for stage in ('minima','electronic','rotors'):
    run_job([RMG_PYTHON,str(HERE/'run_series.py'),'--stage',stage],
            SCRATCH/'pipeline'/stage,CPUS,'production stage '+stage)
print('I048 declared molecular producers complete; thermochemistry and full Verifier remain',flush=True)
