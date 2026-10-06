"""Close the authorized probe when producers finish or the global cap ends.

Command: PYTHONPATH=$PWD rmg_env python .../i048_probe/finish_series.py
Persists each command's streams, runs the full Verifier, and commits locally.
Available-only numerical audits never stand in for the full Verifier.
"""
from __future__ import annotations
import fcntl
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import time
from common import SCRATCH,HERE,CPUS,RMG_PYTHON,load,save,production_deadline


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extension',action='store_true')
    args=parser.parse_args()
    os.sched_setaffinity(0,CPUS)
    deadline=production_deadline()
    os.environ['I048_REPLAY_DEADLINE_UNIX']=str(deadline)
    (SCRATCH/'logs').mkdir(exist_ok=True)
    records=[]
    affinity=','.join(map(str,CPUS))
    prefix='extension-' if args.extension else ''
    def run(arguments,stem,quantum=False):
        stem=prefix+stem
        command=[RMG_PYTHON,str(HERE/arguments[0]),*arguments[1:]]
        if quantum and time.time()<deadline-31:
            command=['timeout','--kill-after=30',str(int(deadline-time.time())-31),*command]
        stdout=SCRATCH/'logs'/(stem+'.stdout.log')
        stderr=SCRATCH/'logs'/(stem+'.stderr.log')
        shell=shlex.join(['taskset','-c',affinity,*command])+' > >(taskset -c '+affinity+' tee -a '+shlex.quote(str(stdout))+') 2> >(taskset -c '+affinity+' tee -a '+shlex.quote(str(stderr))+' >&2)'
        print('I048 FINALIZER '+stem,flush=True)
        status=subprocess.run(['bash','-o','pipefail','-c',shell]).returncode
        records.append({'command':command,'exit_code':status,'stdout':str(stdout),'stderr':str(stderr),'finished_unix_s':time.time()})
        save(SCRATCH/'verification'/(prefix+'final_commands.json'),records)
        return status
    def required(arguments,stem):
        status=run(arguments,stem)
        if status: raise RuntimeError('finalization failed: '+stem+'; exit '+str(status))
    with (SCRATCH/'finalization.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        while time.time()<deadline:
            path=SCRATCH/'logs/pipeline-budgeted.stdout.log'
            ready=(SCRATCH/'pipeline/extension/completed.json').exists() if args.extension else (path.exists() and 'I048 declared molecular producers complete; thermochemistry and full Verifier remain' in path.read_text()[-1500:])
            if ready: break
            time.sleep(min(30,max(.1,deadline-time.time())))
        if time.time()>=deadline:
            # Let capped supervisors persist their actual failed-job receipts.
            time.sleep(10)
        required(['thermochemistry.py','--available-only','--samples','1024','--no-roots'],'final-base-thermo')
        required(['thermochemistry.py','--available-only','--samples','2048'],'final-thermo')
        data=load(SCRATCH/'increments_available_2048.json')
        complete=data['complete']
        if complete:
            for count in (1024,2048):
                save(SCRATCH/f'increments_{count}.json',load(SCRATCH/f'increments_available_{count}.json'))
        report=['report.py']+([] if complete else ['--available-only'])
        required(report+['--refresh-cost','--write-report'],'final-report-before-verifier')
        full_status=run(['verify_results.py','--replay'],'final-full-verifier',quantum=True)
        required(report+['--refresh-cost','--write-report'],'final-report')
        audit=['verify_results.py']+([] if complete else ['--available-only'])
        audit_status=run(audit,'final-numerical-audit')
        if audit_status in (124,137,143,-9,-15) and time.time()>=deadline-1:
            # The hard guardian also covers the verifier's direct quantum
            # checks. Retry an interrupted cached audit after that guardian
            # has ended; this retry cannot start any quantum producer.
            time.sleep(max(0,deadline+5-time.time()))
            required(report+['--refresh-cost','--write-report'],'final-report-after-cap')
            required(audit,'final-numerical-audit-after-cap')
        elif audit_status:
            raise RuntimeError('final numerical audit failed; exit '+str(audit_status))
        required(report,'final-report-check')
        save(SCRATCH/('monitor_stop_extension.json' if args.extension else 'monitor_stop.json'),{'unix_s':time.time(),'reason':'scientific work and author verification ended'})
        branch=subprocess.check_output(['git','branch','--show-current'],text=True).strip()
        if branch!='i048-oligomer-series': raise RuntimeError('unexpected branch: '+branch)
        paths=[HERE.parent/'I048_oligomer_series.md',HERE/'README.md',*sorted(HERE.glob('*.py'))]
        relative=[str(path.relative_to(Path.cwd())) for path in paths]
        staged=subprocess.check_output(['git','diff','--cached','--name-only'],text=True).splitlines()
        if set(staged)-set(relative): raise RuntimeError('unrelated staged changes prevent the local commit')
        subprocess.run(['git','add',*relative],check=True)
        subprocess.run(['git','diff','--cached','--check'],check=True)
        subject='kmc: extend pentamer thermochemistry probe' if args.extension else 'kmc: probe oligomer thermochemistry increments'
        subprocess.run(['git','commit','-m',subject],check=True)
        sha=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
        summary={'SHA':sha,'complete':complete,'available_species':data['available_species'],
                 'full_verifier_exit_code':full_status,'numerical_audit_exit_code':0,
                 'deadline_unix_s':deadline,'finished_unix_s':time.time(),'commands':records}
        save(SCRATCH/'verification'/(prefix+'final_summary.json'),summary)
        with Path('/tmp/i048-live-worker-state.md').open('a') as state:
            proof='verification/'+prefix+'final_summary.json and logs/'+prefix+'final-*'
            state.write('\nFINALIZER DONE: SHA '+sha+'; coverage '+str(len(data['available_species']))+'/25; full --replay exit '+str(full_status)+'; final numerical audit/report check exit0. Proof '+proof+'. No push. Root must inspect actual summary/logs/gitstatus and report exact three sections.\n')
        print('I048 FINALIZER DONE '+sha+' coverage '+str(len(data['available_species']))+'/25 full_verifier_exit='+str(full_status),flush=True)


if __name__=='__main__':
    main()
