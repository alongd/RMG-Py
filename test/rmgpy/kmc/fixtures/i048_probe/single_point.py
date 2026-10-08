"""Declared composite single point, one BLAS thread, 5 GB PySCF workspace.

Command: pyscf_env python .../i048_probe/single_point.py --xyz FILE --output FILE --level pbe
"""
import argparse
import importlib.util
import os
import tempfile
import time
from common import SCRATCH, PLAN, digest, save
from pyscf import dft,gto,lib,__version__

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--xyz',type=__import__('pathlib').Path,required=True)
parser.add_argument('--output',type=__import__('pathlib').Path,required=True)
parser.add_argument('--level',choices=('pbe','blyp'),required=True)
args=parser.parse_args()
if SCRATCH not in args.output.resolve().parents:
    raise ValueError('output outside dispatched scratch')
args.output.parent.mkdir(parents=True,exist_ok=True)
temporary=args.output.parent/'tmp'
temporary.mkdir(exist_ok=True)
tempfile.tempdir=lib.param.TMPDIR=str(temporary)
lib.num_threads(min(PLAN['DFT_native_OpenMP_threads'],len(os.sched_getaffinity(0))))
control_path='/home/alon/anaconda3/envs/rmg_env/lib/python3.9/site-packages/threadpoolctl.py'
control_spec=importlib.util.spec_from_file_location('i048_blas_control',control_path)
control=importlib.util.module_from_spec(control_spec)
control_spec.loader.exec_module(control)
lines=args.xyz.read_text().splitlines()
mol=gto.M(atom='\n'.join(lines[2:2+int(lines[0])]),basis='def2-svp',verbose=4)
mol.max_memory=PLAN['DFT']['max_memory_MB_per_job']
calculation=dft.RKS(mol,xc=args.level).density_fit()
calculation.disp='d3bj'
calculation.grids.level=3
calculation.conv_tol=1e-9
calculation.max_cycle=100
start=time.monotonic()
with control.threadpool_limits(limits=1,user_api='blas'):
    energy=float(calculation.kernel())
    threadpools=control.threadpool_info()
    if any(pool['num_threads']!=1 for pool in threadpools if pool['user_api']=='blas'):
        raise RuntimeError('BLAS thread cap was not enforced')
if not calculation.converged:
    raise RuntimeError('SCF did not converge')
save(args.output,{'pyscf_version':__version__,'level':args.level+'-D3(BJ)/def2-SVP//GFN2-xTB',
     'input_sha256':digest(args.xyz),'energy_hartree':energy,
     'dispersion_hartree':float(calculation.get_dispersion()),'electronic_s':time.monotonic()-start,
     'atoms':[mol.atom_symbol(i) for i in range(mol.natm)],'coordinates_angstrom':mol.atom_coords(unit='Angstrom').tolist(),
     'composite_single_point':True,'BLAS_threads':1,'cpu_affinity':sorted(os.sched_getaffinity(0)),
     'native_OpenMP_threads':lib.num_threads(),'threadpool_controller_version':control.__version__,
     'threadpools':[{'user_api':pool['user_api'],'internal_api':pool['internal_api'],
                    'num_threads':pool['num_threads']} for pool in threadpools]})
print('I048 fresh '+args.level+' energy '+str(args.output),flush=True)
