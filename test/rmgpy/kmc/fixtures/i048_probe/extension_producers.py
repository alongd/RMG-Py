"""Finish the owner-authorized pentamer extension without changing science.

Command: rmg_env python .../i048_probe/extension_producers.py
The two fresh searches with frozen settings and first electronic batch launch first.
Only one three-lane electronic batch and one eight-lane rotor batch run at once.
"""
from concurrent.futures import ThreadPoolExecutor
import fcntl
import os
import time

from common import SCRATCH,HERE,CPUS,RMG_PYTHON,load,save,production_deadline
from run_series import run_job

NAMES=('ps5_01001','ps5_01010','ps5_01110')


def stage(name,kind):
    command=[RMG_PYTHON,str(HERE/'run_series.py'),'--stage',kind,'--species',name]
    if kind in ('electronic','rotors'):
        command.append('--bulk-selected')
    result=run_job(command,SCRATCH/'pipeline/extension'/name/kind,CPUS,
                   'production stage extension '+kind+' '+name)
    save(SCRATCH/'pipeline/extension'/name/(kind+'_completed.json'),result)
    return name


def main():
    os.sched_setaffinity(0,CPUS)
    directory=SCRATCH/'pipeline/extension'
    directory.mkdir(parents=True,exist_ok=True)
    with (directory/'producer.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        states={NAMES[0]:'external electronic',NAMES[1]:'search',NAMES[2]:'search'}
        futures={}
        with ThreadPoolExecutor(max_workers=2) as minima, ThreadPoolExecutor(max_workers=1) as electronic, ThreadPoolExecutor(max_workers=1) as rotors:
            while time.time()<production_deadline():
                for name,(kind,future) in list(futures.items()):
                    if not future.done():continue
                    future.result()
                    del futures[name]
                    states[name]={'minima':'electronic ready','electronic':'rotors ready','rotors':'done'}[kind]
                    if kind=='electronic':
                        save(SCRATCH/'pipeline/early'/name/'completed.json',{'name':name,'unix_s':time.time(),'stages':['minima','electronic'],'owner_extension':True})
                if states[NAMES[0]]=='external electronic':
                    stdout=SCRATCH/'logs/extension-electronic-ps5_01001.stdout.log'
                    stderr=SCRATCH/'logs/extension-electronic-ps5_01001.stderr.log'
                    if stdout.exists() and 'I048 stage complete: electronic' in stdout.read_text()[-1000:]:
                        states[NAMES[0]]='rotors ready'
                        save(SCRATCH/'pipeline/early'/NAMES[0]/'completed.json',{'name':NAMES[0],'unix_s':time.time(),'stages':['minima','electronic'],'owner_extension':True})
                    elif stderr.exists() and 'Traceback' in stderr.read_text():
                        raise RuntimeError('initial extension electronic batch failed; inspect persisted streams')
                for name in NAMES[1:]:
                    if states[name]=='search' and (SCRATCH/'ensembles'/name/'completed.json').exists():
                        states[name]='minima running'
                        futures[name]=('minima',minima.submit(stage,name,'minima'))
                if states[NAMES[0]]!='external electronic' and not any(kind=='electronic' for kind,_ in futures.values()):
                    ready=next((name for name in NAMES if states[name]=='electronic ready'),None)
                    if ready:
                        states[ready]='electronic running'
                        futures[ready]=('electronic',electronic.submit(stage,ready,'electronic'))
                if not any(kind=='rotors' for kind,_ in futures.values()):
                    ready=next((name for name in NAMES if states[name]=='rotors ready'),None)
                    if ready:
                        states[ready]='rotors running'
                        futures[ready]=('rotors',rotors.submit(stage,ready,'rotors'))
                save(directory/'state.json',{'unix_s':time.time(),'states':states,'deadline_unix_s':production_deadline(),'science_changes':False})
                if all(state=='done' for state in states.values()):
                    save(directory/'completed.json',{'unix_s':time.time(),'states':states,'science_changes':False})
                    print('I048 owner extension molecular producers complete; thermochemistry and full Verifier remain',flush=True)
                    return
                time.sleep(min(15,max(.1,production_deadline()-time.time())))
        raise TimeoutError('owner extension global deadline exhausted')


if __name__=='__main__':
    main()
