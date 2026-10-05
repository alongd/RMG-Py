"""Integrate one frozen basin with the unchanged I043 quadrature.

Command: rmg_env python .../i048_probe/rotor_worker.py NAME --seal
Command: same, NAME --index I --samples 2048
The seal serializes symmetry closure; workers never rewrite the shared pool.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from common import SCRATCH,PLAN,bootstrap,digest,load,save


def run(name,index=None,samples=None,seal=False,scratch=SCRATCH):
    bootstrap()
    import rotors
    import complete_symmetry
    selection_path=scratch/'composite'/name/'selection.json'
    receipt=selection_path.parent/'rotor_partition.json'
    if seal:
        # Empty work selection still constructs the entire orbit and partition.
        # No XTB calculator is initialized by this call.
        rotors.integrate(scratch,name,PLAN['new_rotor_counts'][0],PLAN['seed'],
                         proposal_kind='correlated',only_indices=[])
        selection=load(selection_path)
        inputs={str(i):{filename:digest(scratch/'xtb_checks'/name/f'{i:04d}'/filename)
                         for filename in ('xtbopt.xyz','hessian','result.json')}
                for i in selection['unique_minima_indices']}
        save(receipt,{'selection_sha256':digest(selection_path),'inputs':inputs,
                      'science_changes':False})
        return
    proof=load(receipt)
    if digest(selection_path)!=proof['selection_sha256'] or proof['science_changes'] is not False:
        raise AssertionError('sealed rotor selection changed')
    selection=load(selection_path)
    if index not in selection['unique_minima_indices']:
        raise AssertionError('worker index is outside the full partition')
    for i,inputs in proof['inputs'].items():
        for filename,expected in inputs.items():
            if digest(scratch/'xtb_checks'/name/f'{int(i):04d}'/filename)!=expected:
                raise AssertionError('sealed rotor geometry/Hessian changed')
    original_augment=complete_symmetry.augment
    original_write=Path.write_text
    def frozen_augment(root,species,actual):
        if root!=scratch or species!=name or actual!=selection:
            raise AssertionError('worker tried to change the rotor partition')
        return actual
    def frozen_write(path,text,*args,**kwargs):
        if path==selection_path:
            if json.loads(text)!=selection or digest(path)!=proof['selection_sha256']:
                raise AssertionError('worker tried to rewrite a changed selection')
            return len(text)
        return original_write(path,text,*args,**kwargs)
    # These adapters affect this isolated worker process only. All integration
    # steps, global centers, projected modes, draws and XTB calls remain I043.
    complete_symmetry.augment=frozen_augment
    Path.write_text=frozen_write
    try:
        rotors.integrate(scratch,name,samples,PLAN['seed'],proposal_kind='correlated',only_indices=[index])
    finally:
        complete_symmetry.augment=original_augment
        Path.write_text=original_write
    if digest(selection_path)!=proof['selection_sha256']:
        raise AssertionError('worker altered the frozen selection')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('name')
    parser.add_argument('--seal',action='store_true')
    parser.add_argument('--index',type=int)
    parser.add_argument('--samples',type=int,choices=PLAN['new_rotor_counts'])
    args=parser.parse_args()
    if not args.seal and (args.index is None or args.samples is None):
        parser.error('a worker needs --index and --samples')
    run(args.name,args.index,args.samples,args.seal)


if __name__=='__main__':
    main()
