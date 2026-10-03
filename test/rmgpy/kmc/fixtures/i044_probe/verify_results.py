"""Independently verify artifact, pinned sources, thermodynamics and report.

Command: see I044_kp_benchmark.md. Reads no network and compiles no event set.
"""

from __future__ import annotations

import argparse
import ast
from bisect import bisect_left
import hashlib
import json
import math
from pathlib import Path
import subprocess
from types import SimpleNamespace

from scipy.optimize import brentq

from rmgpy.data.rmg import RMGDatabase
from rmgpy.kmc.ssa import record_propensity
from rmgpy.kmc.state import AtomRef, Site
from rmgpy.molecule.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species

from render_tables import check_report
from run_probe import ARTIFACT, DATABASE, DATABASE_SHA, LITERATURE, ROOT, digest


def close(a, b, relative=3e-10, absolute=1e-10):
    assert math.isclose(a, b, rel_tol=relative, abs_tol=absolute), (a, b)


def interpolate(table, t):
    grid, values = table['T'], table['k']
    if not grid[0] <= t <= grid[-1]:
        return None
    i = bisect_left(grid, t)
    if grid[i] == t:
        return values[i]
    f = (t - grid[i-1]) / (grid[i] - grid[i-1])
    return math.exp((1-f) * math.log(values[i-1]) + f * math.log(values[i]))


def reaction(db, record):
    sides = []
    for key in ('reactant_graphs', 'product_graphs'):
        species = [Species(molecule=[Molecule().from_adjacency_list(a)]) for a in record[key]]
        for s in species:
            s.thermo = db.thermo.get_thermo_data(s)
        sides.append(species)
    return Reaction(reactants=sides[0], products=sides[1])


def kc(rxn, t, gas_constant):
    # Explicit independent pressure-standard to dimensional concentration Kc.
    dh = sum(s.thermo.get_enthalpy(t) for s in rxn.products) - sum(s.thermo.get_enthalpy(t) for s in rxn.reactants)
    ds = sum(s.thermo.get_entropy(t) for s in rxn.products) - sum(s.thermo.get_entropy(t) for s in rxn.reactants)
    dn = len(rxn.products) - len(rxn.reactants)
    return math.exp(-(dh - t*ds) / (gas_constant*t)) * (100000.0 / (gas_constant*t))**dn


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    parser.add_argument('--scratch', type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.results.read_text())
    assert data['database_sha'] == DATABASE_SHA
    assert data['artifact_sha256'] == digest(ARTIFACT) == ARTIFACT.stem
    artifact = json.loads(ARTIFACT.read_text())
    assert data['record_count'] == len(artifact['records']) == 14998
    assert data['artifact_rmgpy_sha'] == artifact['provenance']['rmgpy_sha']
    assert data['literature'] == LITERATURE
    assert data['source_sha256']['rmgpy/kmc/compiler.py'] == artifact['provenance']['compiler_sha256']
    for path, sha in data['source_sha256'].items():
        assert digest(ROOT / path) == sha, path
    snapshot = args.scratch / 'database'
    paths = subprocess.check_output(
        ['git', '-C', str(DATABASE), 'ls-tree', '-r', '--name-only', DATABASE_SHA, '--',
         *data['snapshot']['prefixes']]).decode().splitlines()
    h = hashlib.sha256()
    for path in paths:
        assert 'catalog' not in Path(path).parts
        blob = subprocess.check_output(['git', '-C', str(DATABASE), 'show', f'{DATABASE_SHA}:{path}'])
        assert (snapshot / path).read_bytes() == blob, path
        h.update(path.encode() + b'\0' + blob)
    assert len(paths) == data['snapshot']['files'] == 65
    assert h.hexdigest() == data['snapshot']['sha256']
    tree = ast.parse((snapshot / 'input/kinetics/families/R_Addition_MultipleBond/rules.py').read_text())
    entries = [n.value for n in tree.body if isinstance(n, ast.Expr) and isinstance(n.value, ast.Call)
               and isinstance(n.value.func, ast.Name) and n.value.func.id == 'entry']
    assert len(entries) == 1
    assert data['parameters']['rule_count'] == len(entries)
    values = {k.arg: k.value for k in entries[0].keywords}
    assert ast.literal_eval(values['label']) == data['parameters']['rule_label'] == 'R_R;YJ'
    assert ast.literal_eval(values['index']) == data['parameters']['rule_index'] == 3000
    assert ast.literal_eval(values['rank']) == data['parameters']['rank'] == 0
    params = {k.arg: ast.literal_eval(k.value) for k in values['kinetics'].keywords}
    assert params['A'][1] == 'cm^3/(mol*s)' and params['E0'][1] == 'kcal/mol'
    a, ea = params['A'][0] * 1e-6, params['E0'][0] * 4184.0
    close(a, data['parameters']['A_m3_mol_s'])
    close(ea, data['parameters']['Ea_J_mol'])
    assert params['n'] == data['parameters']['n'] == 0
    assert params['alpha'] == data['parameters']['alpha'] == 0
    assert params['Tmin'][1] == params['Tmax'][1] == 'K'
    close(params['Tmin'][0], data['parameters']['Tmin_K'])
    close(params['Tmax'][0], data['parameters']['Tmax_K'])
    training_tree = ast.parse((snapshot / 'input/kinetics/families/R_Addition_MultipleBond/training/reactions.py').read_text())
    training_entries = [n.value for n in training_tree.body
                        if isinstance(n, ast.Expr) and isinstance(n.value, ast.Call)
                        and isinstance(n.value.func, ast.Name) and n.value.func.id == 'entry']
    assert len(training_entries) == data['parameters']['training_entry_count'] == 2962
    gas_constant = data['R_J_mol_K']
    close(gas_constant, 8.314472)
    forward = lambda t: a * math.exp(-ea/(gas_constant*t))
    benchmark = lambda t: LITERATURE['iupac']['A_L_mol_s'] * math.exp(-LITERATURE['iupac']['Ea_J_mol']/(gas_constant*t))
    db = RMGDatabase()
    db.load_thermo(str(snapshot / 'input/thermo'), thermo_libraries=['primaryThermoLibrary'], depository=True)
    assert len(data['pairs']) == len(artifact.get('ps_primary_end_ceiling_pairs', artifact['ps_ceiling_pairs'])) == 2
    checked = 0
    for pair_data in data['pairs']:
        pair = pair_data['pair']
        assert pair in artifact.get('ps_primary_end_ceiling_pairs', artifact['ps_ceiling_pairs'])
        prop, dep = pair_data['propagation'], pair_data['depropagation']
        for record in (prop, dep):
            assert record == next(r for r in artifact['records'] if r['event_id'] == record['event_id'])
            assert record['degeneracy'] == record['ssa_multiplier'] == record['raw_path_degeneracy'] == 1
        rxn = reaction(db, prop)
        assert len(prop['k_table']['T']) == len(prop['k_table']['k']) == len(dep['k_table']['k']) == len(prop['thermo_provenance']['equilibrium_constant_table']['Kc']) == pair_data['grid_points_checked'] == 29
        for t, kf, kr, kt in zip(prop['k_table']['T'], prop['k_table']['k'], dep['k_table']['k'],
                                 prop['thermo_provenance']['equilibrium_constant_table']['Kc']):
            close(forward(t), kf)
            close(kc(rxn, t, gas_constant), kt)
            close(forward(t)/kc(rxn, t, gas_constant), kr)
            checked += 1
    assert checked == 58
    first = data['pairs'][0]
    prop, dep = first['propagation'], first['depropagation']
    rxn = reaction(db, prop)
    for row in data['comparison']:
        t = row['T_K']
        close(row['rule_L_mol_s'], forward(t)*1000)
        close(row['benchmark_L_mol_s'], benchmark(t))
        close(row['rule_over_benchmark'], forward(t)*1000/benchmark(t))
        close(row['A_factor'], a*1000/LITERATURE['iupac']['A_L_mol_s'])
        close(row['Ea_factor'], math.exp((LITERATURE['iupac']['Ea_J_mol']-ea)/(gas_constant*t)))
        k = interpolate(prop['k_table'], t)
        if k is None:
            assert all(row[key] is None for key in ('ssa_L_mol_s', 'ssa_over_benchmark', 'interpolation_factor'))
        else:
            close(row['ssa_L_mol_s'], k*1000)
            close(row['ssa_over_benchmark'], k*1000/benchmark(t))
            close(row['interpolation_factor'], k/forward(t))
    for row in data['reverse']:
        t = row['T_K']
        kf, kr, kct = interpolate(prop['k_table'], t), interpolate(dep['k_table'], t), kc(rxn, t, gas_constant)
        close(row['kf_L_mol_s'], kf*1000)
        close(row['kr_s_1'], kr)
        close(row['Kc_m3_mol'], kct)
        close(row['M_eq_SSA_mol_L'], kr/kf/1000)
        close(row['M_eq_thermo_mol_L'], 1/kct/1000)
        close(row['kr_with_benchmark_same_Kc_s_1'], benchmark(t)/1000/kct)
        close(row['forward_per_site_at_1M_s_1'], kf*1000)
        close(row['reverse_over_forward_at_1M'], kr/kf/1000)
    for row in data['equilibrium']:
        datum = LITERATURE[row['key']]
        t, ce = row['T_K'], 1/kc(rxn, row['T_K'], gas_constant)/1000
        assert t == datum['T_K'] and row['literature_M_mol_L'] == datum['M_mol_L']
        close(row['model_M_mol_L'], ce)
        close(row['model_over_literature'], ce/datum['M_mol_L'])
        close(row['delta_G_required_J_mol'], gas_constant*t*math.log(datum['M_mol_L']/ce))
        close(row['model_rule_kr_s_1'], forward(t)/kc(rxn,t,gas_constant))
        close(row['conditional_radical_kr_from_benchmark_s_1'], benchmark(t)*datum['M_mol_L'])
        close(row['model_rule_over_conditional_kr'], row['model_rule_kr_s_1']/row['conditional_radical_kr_from_benchmark_s_1'])
        kr = interpolate(dep['k_table'], t)
        if kr is None:
            assert row['ssa_kr_s_1'] is None
        else:
            close(row['ssa_kr_s_1'], kr)
    root = brentq(lambda t: math.log(kc(rxn,t,gas_constant)*1000), 600, 800)
    close(root, data['Tc_thermo_K'])
    grid_root = brentq(lambda t: math.log(interpolate(prop['k_table'],t)*1000/interpolate(dep['k_table'],t)), 700, 725)
    close(grid_root, data['Tc_artifact_grid_K'])
    control_data = data['structure_control']
    s = [Species(molecule=[Molecule(smiles=smiles)]) for smiles in control_data['species_smiles']]
    for item in s:
        item.generate_resonance_structures()
        item.thermo = db.thermo.get_thermo_data(item)
    control = Reaction(reactants=s[:2], products=s[2:])
    ct = brentq(lambda t: math.log(kc(control,t,gas_constant)*1000),600,800)
    close(ct,control_data['Tc_K'])
    close(ct-root,control_data['Tc_change_K'])
    close(control.get_enthalpy_of_reaction(298.15)-rxn.get_enthalpy_of_reaction(298.15),control_data['H_change_J_mol_at_298_15'])
    close(control.get_entropy_of_reaction(298.15)-rxn.get_entropy_of_reaction(298.15),control_data['S_change_J_mol_K_at_298_15'])
    # Exercise the public propensity function with supplied eligible candidates.
    # This checks count/volume normalization, not production site discovery.
    sites = [Site('probe', str(i), 0, 1, {1:{'atom_ref':AtomRef(str(i),str(i),0)}}) for i in range(5)]
    state = SimpleNamespace(components={str(i):i for i in range(5)})
    candidates = tuple((r,m) for r in sites[:2] for m in sites[2:])
    index = SimpleNamespace(candidates=lambda event_id: candidates)
    volume = 1e-21
    close(record_propensity(state,index,prop,700.0,volume),forward(700)*6/(6.02214076e23*volume))
    index = SimpleNamespace(candidates=lambda event_id: tuple((s,) for s in sites[:2]))
    close(record_propensity(state,index,dep,700.0,volume),interpolate(dep['k_table'],700)*2)
    print('I044 artifact, 65 pinned database files and product-source hashes verified')
    print('I044 both pairs: default Arrhenius, dimensional Kc and reverse rates independently verified')
    print('I044 benchmark ratios, A/Ea factors, equilibrium comparisons and crossings independently verified')
    print('I044 public SSA propensity count/volume normalization verified')
    check_report(data)


if __name__ == '__main__':
    main()
