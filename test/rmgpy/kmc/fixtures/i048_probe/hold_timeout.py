"""Preserve an active search through an unnecessarily short wrapper timer.

Command: rmg_env python .../i048_probe/hold_timeout.py NAME TIMEOUT_PID
Only the worker-owned GNU timeout is held. CREST and xTB keep running.
The original 48-hour deadline remains enforced. The wrapper's actual status
124 is retained; acceptance requires the timed scientific child's exit 0.
"""
import argparse
import os
from pathlib import Path
import re
import signal
import subprocess
import tempfile
import time

import psutil
from common import SCRATCH,CPUS,PLAN,digest,load,save
from run_series import finalize_search


def self_check():
    """Reproduce wrapper 124 and successful timed child independently."""
    parent=SCRATCH/'verification/timeout_protocol'
    parent.mkdir(parents=True,exist_ok=True)
    directory=Path(tempfile.mkdtemp(prefix='probe_',dir=parent))
    resources=directory/'resources.txt'
    command=['timeout','--kill-after=30','2','/usr/bin/time','-v','-o',str(resources),'sleep','6']
    started=time.time()
    with (directory/'stdout.log').open('w') as stdout,(directory/'stderr.log').open('w') as stderr:
        child=subprocess.Popen(command,stdout=stdout,stderr=stderr)
        time.sleep(.25)
        assert time.time()<started+2
        os.kill(child.pid,signal.SIGSTOP)
        try:
            while not resources.exists() or 'Exit status:' not in resources.read_text():
                if time.time()-started>20:
                    raise TimeoutError('controlled child did not finish')
                time.sleep(.1)
        finally:
            os.kill(child.pid,signal.SIGCONT)
        status=child.wait(timeout=5)
    assert status==124 and re.search(r'Exit status:\s*0\s*$',resources.read_text())
    save(directory/'proof.json',{'command':command,'wrapper_exit_code':status,
                                'timed_child_exit_code':0,'elapsed_s':time.time()-started,
                                'resources_sha256':digest(resources)})
    print('I048 controlled timeout hold: wrapper 124; timed child 0',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('name',nargs='?')
    parser.add_argument('pid',nargs='?',type=int)
    parser.add_argument('--self-check',action='store_true')
    args=parser.parse_args()
    if args.self_check:
        self_check()
        return
    if not args.name or args.pid is None:
        parser.error('NAME and TIMEOUT_PID are required for an active guardian')
    directory=SCRATCH/'ensembles'/args.name
    wrapper=psutil.Process(args.pid)
    command=wrapper.cmdline()
    assert command[:3]==['timeout','--kill-after=30','28800']
    assert wrapper.cwd()==str(directory) and wrapper.uids().real==os.getuid()
    assert set(wrapper.cpu_affinity())<=set(CPUS)
    crest=command.index('/home/alon/anaconda3/envs/crest_env/bin/crest')
    scientific=command[crest:]
    assert scientific==[scientific[0],'input.xyz','--gfn2','--quick','--ewin','6','--T','4']
    assert '--T' in scientific and PLAN['budget']['cores']==8
    resources=Path(command[command.index('-o')+1])
    assert resources.parent==directory
    assert 'Exit status:' not in (resources.read_text() if resources.exists() else '')
    children=wrapper.children(recursive=True)
    assert any(p.cmdline()==scientific for p in children)
    deadline=load(SCRATCH/'started.json')['unix_s']+48*3600
    proof={'name':args.name,'timeout_PID':args.pid,'wrapper_command':command,
           'command':scientific,'resources_path':str(resources),
           'original_timeout_unix_s':wrapper.create_time()+28800,
           'production_deadline_unix_s':deadline,'scientific_completed':False,
           'science_changes':False,'reason':'eight-hour per-search timer is shorter than the dispatched total budget'}
    assert time.time()<proof['original_timeout_unix_s']
    os.kill(args.pid,signal.SIGSTOP)
    proof['stopped_unix_s']=time.time()
    save(directory/'timeout_hold.json',proof)
    print('I048 held only timeout wrapper for '+args.name,flush=True)
    try:
        while time.time()<deadline-31:
            if resources.exists():
                text=resources.read_text()
                if re.search(r'Exit status:\s*0\s*$',text):
                    stdout=directory/'stdout.log'
                    if 'CREST terminated normally.' not in stdout.read_text():
                        time.sleep(1)
                        continue
                    proof.update(scientific_completed=True,child_completed_unix_s=time.time(),
                                 resources_sha256=digest(resources),stdout_sha256=digest(stdout))
                    save(directory/'timeout_hold.json',proof)
                    break
                if 'Exit status:' in text:
                    raise RuntimeError('held search exited unsuccessfully')
            time.sleep(1)
    finally:
        os.kill(args.pid,signal.SIGCONT)
        proof['resumed_unix_s']=time.time()
        save(directory/'timeout_hold.json',proof)
    if not proof['scientific_completed']:
        raise TimeoutError('original production deadline reached; timeout wrapper resumed')
    for _ in range(20):
        job=directory/'job.json'
        if job.exists():
            record=load(job)
            proof.update(wrapper_exit_code=record['exit_code'],job_sha256=digest(job))
            save(directory/'timeout_hold.json',proof)
            finalize_search(args.name,record)
            print('I048 preserved wrapper status %d; CREST exit 0; search finalized %s'%(
                record['exit_code'],args.name),flush=True)
            return
        time.sleep(1)
    raise RuntimeError('original search supervisor did not record the wrapper status')


if __name__=='__main__':
    main()
