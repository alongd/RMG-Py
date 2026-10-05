"""Reuse I043 scientific CLIs with I048 paths and inherited CPU affinity.

Command: rmg_env python .../i048_probe/legacy_cli.py check_xtb --species NAME --indices I
Command: rmg_env python .../i048_probe/legacy_cli.py rotors --species NAME --samples 2048 --proposal correlated
"""
import sys
import importlib
import subprocess
import os
from common import SCRATCH,load,save
from common import bootstrap

bootstrap()
name = sys.argv.pop(1)
if name not in ('check_xtb','rotors','complete_symmetry'):
    raise ValueError('only specified scientific producers may be replayed')
module = importlib.import_module(name)
original_run = subprocess.run
def confined_run(command, *args, **kwargs):
    if isinstance(command, list):
        command = [part.replace('--parallel 8','--parallel 1') if isinstance(part,str) else part for part in command]
    if 'env' in kwargs:
        kwargs['env'] = dict(kwargs['env'], OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    return original_run(command,*args,**kwargs)
subprocess.run = confined_run
module.main()
if name=='check_xtb':
    # The unchanged producer records its original parallel flag. Record the
    # actual adapter execution rather than leaving a misleading command.
    arguments=sys.argv[1:]
    start=arguments.index('--species')+1
    names=[]
    for argument in arguments[start:]:
        if argument.startswith('--'):
            break
        names.append(argument)
    indices=None
    if '--indices' in arguments:
        indices=set()
        for argument in arguments[arguments.index('--indices')+1:]:
            if argument.startswith('--'):
                break
            indices.add(int(argument))
    for species in names:
        if not species.startswith('ps'):
            continue
        stage='xtb_pilots' if '--starting-input' in arguments else 'xtb_checks'
        for path in (SCRATCH/stage/species).glob('*/result.json'):
            # Concurrent lanes own disjoint candidate records. Never annotate
            # another lane's output or race its atomic result-file write.
            if indices is not None and int(path.parent.name) not in indices:
                continue
            record=load(path)
            if record.get('execution_adapter'):
                # A cached candidate keeps the actual computation's metadata.
                continue
            command=record.get('command',[])
            if '--parallel' in command:
                command[command.index('--parallel')+1]='1'
                record['execution_adapter']={'OMP_NUM_THREADS':1,'BLAS_threads':1,
                                             'cpu_affinity':sorted(os.sched_getaffinity(0))}
                save(path,record)
