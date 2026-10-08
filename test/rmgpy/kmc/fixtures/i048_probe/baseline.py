"""Read-only git-show materialization and pinned GAV oligomer increments.

Command: PYTHONPATH=$PWD rmg_env python .../i048_probe/baseline.py
"""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

from common import SCRATCH, DB_SHA, TEMPERATURES, bootstrap, save

import numpy as np
from scipy.optimize import brentq

SPECIES, SEQUENCES = bootstrap()
from rmgpy import constants
from rmgpy.data.rmg import RMGDatabase
from rmgpy.molecule import Molecule
from rmgpy.species import Species

FIXTURES=Path(__file__).resolve().parent.parent
spec=importlib.util.spec_from_file_location('i039_reference',FIXTURES/'i039_probe/run_probe.py')
i039=importlib.util.module_from_spec(spec)
spec.loader.exec_module(i039)


def build():
    if i039.prior.DATABASE_SHA!=DB_SHA:
        raise AssertionError('database commit mismatch')
    snapshot=i039.prior.snapshot_database(Path('/home/alon/Code/RMG-database'),SCRATCH/'baseline/database')
    db=RMGDatabase()
    db.load_thermo(str(SCRATCH/'baseline/database/input/thermo'),thermo_libraries=['primaryThermoLibrary'])
    db.load_kinetics(str(SCRATCH/'baseline/database/input/kinetics'),reaction_libraries=[],
                     kinetics_families=list(i039.prior.PS_FAMILY_CANDIDATES),kinetics_depositories=['training'])
    propagation,compiled=i039.compile_pair(db,3)
    for item in propagation.reactants+propagation.products:
        item.thermo=db.thermo.get_thermo_data(item)
    molecules={}
    result={'snapshot':snapshot,'database_SHA':DB_SHA,'R':constants.R,'p0_Pa':100000.,'c0_mol_m3':1000.,
            'continuous_Tc_K':brentq(lambda t:math.log(propagation.get_equilibrium_constant(t)*1000),250,1800),
            'species':{},'averages':{},'increments':{},'dense_GAV_increments':{},
            'propagation': [{'T_K':float(t),'H_J_mol':propagation.get_enthalpy_of_reaction(t),
                             'S_J_mol_K':propagation.get_entropy_of_reaction(t)}
                            for t in sorted(set([*np.arange(250.,1801.),*TEMPERATURES]))]}
    for name,smiles in SPECIES.items():
        item=Species(molecule=[Molecule(smiles=smiles)])
        item.generate_resonance_structures()
        item.thermo=db.thermo.get_thermo_data(item)
        molecules[name]=item
        decomposition=i039.decompose(db.thermo,item)
        result['species'][name]={
            'adjacency':item.molecule[0].to_adjacency_list(),'sigma':float(item.get_symmetry_number()),
            'optical_atom_half_factors':decomposition['optical_atom_half_factors'],
            'source_weights':decomposition['source_weights'],'thermo_comment':item.thermo.comment,
            'thermo':[{'T_K':t,'H_J_mol':item.thermo.get_enthalpy(t),'S_J_mol_K':item.thermo.get_entropy(t)} for t in TEMPERATURES]}
    for n,classes in SEQUENCES.items():
        result['averages'][n]=[{'T_K':t,
            'H_J_mol':sum(row['weight']*molecules[name].thermo.get_enthalpy(t) for name,row in classes.items()),
            'S_J_mol_K':sum(row['weight']*molecules[name].thermo.get_entropy(t) for name,row in classes.items())} for t in TEMPERATURES]
    reference={'cumene':-1.,'n_propylbenzene':-1.,'ethylbenzene':1.,'ethane':1.}
    for n in range(2,5):
        coefficients=dict(reference)
        coefficients.update({name:-row['weight'] for name,row in SEQUENCES[str(n)].items()})
        coefficients.update({name:row['weight'] for name,row in SEQUENCES[str(n+1)].items()})
        vector={}
        for name,c in coefficients.items():
            for source,w in result['species'][name]['source_weights'].items():
                vector[source]=vector.get(source,0.)+c*w
        vector={source:w for source,w in vector.items() if abs(w)>1e-10}
        result['increments'][str(n)]={'balanced_coefficients':coefficients,'balanced_source_vector':vector,
            'balanced':[{'T_K':t,'H_J_mol':sum(c*molecules[name].thermo.get_enthalpy(t) for name,c in coefficients.items()),
                         'S_J_mol_K':sum(c*molecules[name].thermo.get_entropy(t) for name,c in coefficients.items())} for t in TEMPERATURES],
            'raw':[{'T_K':t,
                'H_J_mol':result['averages'][str(n+1)][j]['H_J_mol']-result['averages'][str(n)][j]['H_J_mol'],
                'S_J_mol_K':result['averages'][str(n+1)][j]['S_J_mol_K']-result['averages'][str(n)][j]['S_J_mol_K']} for j,t in enumerate(TEMPERATURES)]}
        chain_weights={name:-r['weight'] for name,r in SEQUENCES[str(n)].items()}
        chain_weights.update({name:r['weight'] for name,r in SEQUENCES[str(n+1)].items()})
        result['dense_GAV_increments'][str(n)]=[{'T_K':row['T_K'],
            'H_J_mol':sum(c*molecules[name].thermo.get_enthalpy(row['T_K']) for name,c in chain_weights.items()),
            'S_J_mol_K':sum(c*molecules[name].thermo.get_entropy(row['T_K']) for name,c in chain_weights.items())}
            for row in result['propagation']]
    return result,propagation,molecules


def main():
    result,_,_=build()
    save(SCRATCH/'baseline/baseline.json',result)
    print('I048 pinned snapshot: '+str(result['snapshot']['files'])+' files; '+result['snapshot']['sha256'])
    print('I048 baseline gas Tc at 1 mol/L: %.6f K'%result['continuous_Tc_K'])
    for n,row in result['increments'].items():
        h=row['raw'][0]['H_J_mol']/1000
        s=row['raw'][0]['S_J_mol_K']
        print('I048 GAV %s->%d: H298 %.6f kJ/mol; S298 %.6f J/mol/K; balanced source vector %s'%(n,int(n)+1,h,s,row['balanced_source_vector']))


if __name__=='__main__':
    main()
