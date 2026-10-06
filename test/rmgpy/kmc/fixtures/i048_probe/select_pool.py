"""Select checked candidates and close the I043 rotor symmetry orbit.

Command: rmg_env python .../i048_probe/select_pool.py NAME
"""
import json
import sys
from common import SCRATCH, bootstrap, load

bootstrap()
import run_composite
name=sys.argv[1]
reduction=load(SCRATCH/'composite'/name/'candidate_reduction.json')
search=load(SCRATCH/'ensembles'/name/'completed.json')
original=json.loads
def narrowed(text,*args,**kwargs):
    data=original(text,*args,**kwargs)
    if isinstance(data,dict) and data.get('ensemble_sha256')==search['ensemble_sha256']:
        data=dict(data,selected_indices_within_12_kJ=reduction['retained_indices'])
    return data
run_composite.json.loads=narrowed
pool=run_composite.select_minima(SCRATCH,name)
print('I048 '+name+': '+str(len(pool['unique_minima_indices']))+' stable rotor-space wells',flush=True)
