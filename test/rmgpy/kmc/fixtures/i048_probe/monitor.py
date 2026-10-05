"""Read-only affinity and RSS snapshots for this probe's calculation processes.

Command: rmg_env python .../i048_probe/monitor.py --once
Command: same without --once; snapshots append every 30 s, within the wall cap.
"""
import argparse
import os
from pathlib import Path
import time
from common import SCRATCH,CPUS,save,load


def snapshot():
    processes={}
    roots=set()
    for directory in Path('/proc').iterdir():
        if not directory.name.isdigit():
            continue
        try:
            status={key:value.strip() for line in (directory/'status').read_text().splitlines()
                    if ':' in line for key,value in [line.split(':',1)]}
            pid=int(directory.name)
            processes[pid]=(directory,status)
            if status['Name'] in ('python','python3'):
                argv=(directory/'cmdline').read_bytes().decode(errors='replace').split('\0')
                # Identify actual probe script invocations, not a compiler or
                # read-only interpreter with a probe path in another argument.
                if len(argv)>1 and '/i048_probe/' in argv[1]:
                    roots.add(pid)
        except (OSError,KeyError):
            continue
    selected=set(roots)
    changed=True
    while changed:
        changed=False
        for pid,(_,status) in processes.items():
            if int(status.get('PPid','0')) in selected and pid not in selected:
                selected.add(pid)
                changed=True
    rows=[]
    for pid in selected:
        directory,status=processes[pid]
        affinities=[]
        try:
            for thread in (directory/'task').iterdir():
                mask=sorted(os.sched_getaffinity(int(thread.name)))
                if not set(mask)<=set(CPUS):
                    save(SCRATCH/('resource_violation_%d.json'%int(time.time()*1e6)),
                         {'unix_s':time.time(),'pid':pid,'thread':int(thread.name),
                          'argv':(directory/'cmdline').read_bytes().decode(errors='replace').split('\0'),
                          'affinity':mask,'declared_cores':list(CPUS)})
                    raise RuntimeError('probe thread escaped its eight-core set: '+str(pid))
                affinities.append(mask)
        except (FileNotFoundError,ProcessLookupError):
            continue
        rows.append({'pid':pid,'name':status['Name'],'RSS_kB':int(status.get('VmRSS','0 kB').split()[0]),
                     'thread_affinities':affinities})
    total=sum(row['RSS_kB'] for row in rows)
    if total>16*1024**2:
        raise RuntimeError('observed calculation RSS exceeded 16 GiB')
    return {'unix_s':time.time(),'RSS_kB':total,'processes':rows}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--once',action='store_true')
    args=parser.parse_args()
    os.sched_setaffinity(0,CPUS)
    while True:
        result=snapshot()
        with (SCRATCH/'resource_snapshots.jsonl').open('a') as handle:
            import json
            handle.write(json.dumps(result)+'\n')
        print('I048 resource snapshot: %.3f GiB, %d calculation processes; all observed threads confined to declared cores'%(
            result['RSS_kB']/1024**2,len(result['processes'])),flush=True)
        if args.once or (SCRATCH/'monitor_stop.json').exists() or time.time()-load(SCRATCH/'started.json')['unix_s']>=48*3600:
            break
        time.sleep(30)
