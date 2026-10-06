"""Enforce the original cap on this user-owned probe process tree.

Command: rmg_env python .../i048_probe/deadline_guard.py (start BEFORE the cap)
No debugger, privileged operation, or signal to unrelated processes is used.
"""
import json,os,pathlib,signal,time,psutil
from common import production_deadline
root=pathlib.Path('/home/alon/runs/i048-oligomer-series')
deadline=production_deadline()
scripts={'run_series.py','pipeline.py','supplemental.py','early_rotors.py','hold_timeout.py','lane1_recovery.py','verify_results.py'}
tracked={};events=[]
def collect():
    for p in psutil.process_iter(['pid','uids','cmdline']):
        try:
            a=p.info['cmdline'] or []
            if p.info['uids'].real!=os.getuid() or len(a)<2 or '/i048_probe/' not in a[1] or pathlib.Path(a[1]).name not in scripts: continue
            if pathlib.Path(a[1]).name=='verify_results.py': tracked[p.pid]=p.create_time()
            for child in p.children(recursive=True):
                if child.uids().real==os.getuid(): tracked[child.pid]=child.create_time()
        except psutil.Error: pass
while time.time()<deadline-31:
    time.sleep(min(30,max(.1,deadline-31-time.time())))
collect()
for pid,born in list(tracked.items()):
    try:
        p=psutil.Process(pid)
        if p.create_time()!=born or p.name()!='timeout': continue
        p.send_signal(signal.SIGTERM)
        events.append({'unix_s':time.time(),'pid':pid,'signal':'SIGTERM','command':p.cmdline()})
    except psutil.Error: pass
while time.time()<deadline-.5:
    collect();time.sleep(min(1,max(.01,deadline-.5-time.time())))
collect()
for pid,born in list(tracked.items()):
    try:
        p=psutil.Process(pid)
        if p.create_time()!=born or p.name() in ('tee','bash'): continue
        command=p.cmdline();p.kill()
        events.append({'unix_s':time.time(),'pid':pid,'signal':'SIGKILL','command':command})
    except psutil.Error: pass
(root/'deadline_guard.json').write_text(json.dumps({'deadline_unix_s':deadline,'events':events,'completed_unix_s':time.time()},indent=2)+'\n')
print('I048 original deadline enforced on own calculation descendants',flush=True)
