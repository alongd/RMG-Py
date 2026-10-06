"""Measure current tetramer composite cost without calculating increments.

Command: rmg_env python .../i048_probe/cost_pilot.py
Command: same, --native-replay (compare the current capped execution setting)
Pilot thermochemistry is not promoted into the ensemble or used to rank it.
"""
from common import SCRATCH,HERE,RMG_PYTHON,DFT_PYTHON,CPUS
from common import load
from run_series import run_job
import argparse

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--native-replay',action='store_true')
args=parser.parse_args()

name='ps4_0000'
if not args.native_replay:
    run_job([RMG_PYTHON,str(HERE/'legacy_cli.py'),'check_xtb','--species',name,'--starting-input'],
            SCRATCH/'pilot'/name/'checks',CPUS[4:],'cost pilot tetramer minimum')
source=SCRATCH/'xtb_pilots'/name/'0000/xtbopt.xyz'
target=SCRATCH/'pilot'/name/('pbe_native_replay' if args.native_replay else 'pbe')/'result.json'
if not target.exists():
    run_job([DFT_PYTHON,str(HERE/'single_point.py'),'--xyz',str(source),
             '--output',str(target),'--level','pbe'],target.parent,CPUS[4:],'cost pilot tetramer PBE'+(' native replay' if args.native_replay else ''))
if args.native_replay:
    first=load(SCRATCH/'pilot'/name/'pbe/result.json')
    replay=load(target)
    difference=(replay['energy_hartree']-first['energy_hartree'])*2625.499639
    if abs(difference)>.02:
        raise AssertionError('execution-setting energy difference exceeds 0.02 kJ/mol')
    print('I048 native-thread pilot energy difference %.9f kJ/mol'%difference,flush=True)
print('I048 cost pilot complete; no increment or Tc computed',flush=True)
