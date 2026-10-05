"""Fixed I048 method and adapters to the committed I043 scientific routines.

Command: rmg_env python test/rmgpy/kmc/fixtures/i048_probe/run_series.py --prepare
"""
from __future__ import annotations

import hashlib
import itertools
import json
import os
from pathlib import Path
import sys

# Apply the numerical thread cap before any scientific-library imports, even
# when a replay command omitted the documented shell environment.
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'

SCRATCH = Path('/home/alon/runs/i048-oligomer-series')
PRIOR = Path('/home/alon/runs/i043-diphenyl-nonadditivity')
HERE = Path(__file__).resolve().parent
LEGACY = HERE.parent / 'i043_probe'
RMG_PYTHON = '/home/alon/anaconda3/envs/rmg_env/bin/python'
DFT_PYTHON = '/home/alon/anaconda3/envs/pyscf_env/bin/python'
CPUS = (0, 2, 4, 6, 8, 10, 12, 14)
DB_SHA = '4a12d36fcdc193ede82c8d1ab5c1653495d445bc'
TEMPERATURES = (298.15, 300., 400., 500., 600., 700., 800.)
PLAN = {
    'caps': 'I043 methyl caps: CH3-[CH(Ph)-CH2]_(n-1)-CH(Ph)-CH3; C(8n+1)H(8n+4)',
    'sequences': 'all distinct classes modulo enantiomerism and end reversal for n=2,3,4,5; frozen unbiased atactic weights from all 2**n oriented assignments',
    'search': 'CREST GFN2 --quick --ewin 6; stochastic; RDKit seed 43002',
    'candidate_window_kJ_mol': 12.,
    'new_candidate_cap': 32,
    'reduction': 'n=4,5: first 32 candidates in increasing GFN2 energy within 12 kJ/mol, tight extreme xTB minima/Hessians, symmetry deduplication; count omitted candidates. No claim of search convergence.',
    'energy': 'PBE-D3(BJ) and BLYP-D3(BJ)/def2-SVP//GFN2; n<=3 use prior complete composites. n=4,5 DFT lowest three unique chemical minima; others GFN2 with offset from lowest GFN2 conformer for each functional. DFT re-ranking and approximation sensitivity reported.',
    'DFT_chemical_minima_per_sequence': 3,
    'DFT': {'density_fitting': True, 'grid': 3, 'scf_Eh': 1e-9, 'three_body_dispersion': False, 'max_memory_MB_per_job': 5000},
    'DFT_native_OpenMP_threads': 4,
    'DFT_execution_lanes': [[0,2,4],[6,8,10],[12,14]],
    'rotors': 'I043 full-dimensional fixed-angle potential, disjoint periodic Voronoi basins, projected non-torsional Hessian and rigid-curvature Pitzer-Gwinn factor; correlated importance proposal, steric cutoff 0.45 A',
    'new_rotor_counts': [1024, 2048],
    'seed': 43003,
    'prior_samples': {'diphenylpentane_meso':4096, 'diphenylpentane_racemo':4096, 'triphenylheptane_iso':4096, 'triphenylheptane_syndio':4096, 'triphenylheptane_hetero':4096, 'ethane':4096, 'ethylbenzene':4096, 'n_propylbenzene':4096, 'cumene':4096},
    'prior_basin_overrides': {'triphenylheptane_iso': {'16':16384}, 'triphenylheptane_hetero': {'22':32768,'28':32768,'31':32768,'46':32768}},
    'local_definition': 'primary local: omit chain translation/rotation from each basin partition and reweight basins; retain the fixed balanced one-ring reference unchanged. Also show component subtraction at gas basin weights and all-cycle external removal as diagnostics.',
    'formation_reference': 'use the I043 balanced one-ring surrogate to anchor H increment at 298.15 K, then use chain thermal increments at other T. Absolute chain entropy increments need no electronic reference and are compared directly with oriented GAV. Balanced-cycle entropy residuals are also reported.',
    'configuration': 'fixed sequences, no equilibrium diastereomer mixture; remove artificial finite-cap end symmetry; existing propagation R ln 2 retained once',
    'convergence': 'compare successive increments at every T; MC errors and population sequence spread reported separately; no fit/extrapolated limit if sampling/truncation or length dependence prevents it',
    'budget': {'cores':8,'BLAS_threads_per_process':1,'total_RAM_GB':16,'wall_hours':48,'n6':'optional, not in initial allocation'},
    'selection': 'no choice uses Tc or comparator agreement',
}


def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def load(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def catalogue():
    from rdkit import Chem
    sys.path.insert(0, str(LEGACY))
    from run_ensembles import SPECIES as old_species
    species = {name: smiles for name,smiles in old_species.items() if name != 'diphenylpropane' and not name.startswith('ps')}
    by_key = {}
    def identity(smiles):
        mol = Chem.MolFromSmiles(smiles)
        mirror = Chem.Mol(mol)
        for atom in mirror.GetAtoms():
            if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                atom.InvertChirality()
        return min(Chem.MolToSmiles(mol), Chem.MolToSmiles(mirror))
    for name,smiles in species.items():
        by_key[identity(smiles)] = name
    sequences = {}
    for n in range(2,6):
        classes = {}
        for bits in itertools.product((0,1), repeat=n):
            smiles = 'C'+''.join('[C'+('@' if b else '@@')+'H](c1ccccc1)C' for b in bits)
            key = identity(smiles)
            name = by_key.get(key, 'ps%d_%s' % (n,''.join(map(str,bits))))
            by_key[key] = name
            species.setdefault(name, smiles)
            row = classes.setdefault(name, {'n':n, 'weight':0., 'assignments':[], 'canonical_mod_mirror':key})
            row['weight'] += 1/2**n
            row['assignments'].append(''.join(map(str,bits)))
        sequences[str(n)] = classes
    return species, sequences


def bootstrap():
    """Bind unchanged I043 routines to the new species and scratch root."""
    species, sequences = catalogue()
    import run_ensembles
    run_ensembles.SPECIES = species
    run_ensembles.SCRATCH = SCRATCH
    return species, sequences


def numerical_environment():
    return dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1', OMP_MAX_ACTIVE_LEVELS='1',
                PYTHONPATH=str(HERE.parents[4]), MPLCONFIGDIR=str(SCRATCH/'mpl'))
