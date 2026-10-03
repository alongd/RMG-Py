"""Audit an existing event set; never compile events or regenerate rate rules.

Command: see I044_kp_benchmark.md (run with --scratch and --output).
Literature inputs below reproduce published transcriptions and arithmetic;
they do not reproduce historical experiments. No network is needed to rerun.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path

from scipy.optimize import brentq

from rmgpy import constants
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.kmc.ssa import record_rate

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
ARTIFACT = Path('/home/alon/runs/polymer/i038-bench/head-artifact/'
                '0883a1292e20708a17ee8dc7960cf5648b354e0c24ce3373dee83813cea85c9d.json')
DATABASE = Path('/home/alon/Code/RMG-database')
DATABASE_SHA = '4a12d36fcdc193ede82c8d1ab5c1653495d445bc'
spec = importlib.util.spec_from_file_location('i034', HERE.parent / 'i034_probe/run_probe.py')
prior = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)

# Rounded, commonly cited benchmark in IUPAC 2019 Table 1. The original
# 1995 report specifies -12 to 93 C; the later summary gives -12 to 90 C.
LITERATURE = {
    'iupac': {
        'A_L_mol_s': 4.27e7, 'Ea_J_mol': 32500.0,
        'range_K': [261.15, 363.15], 'original_upper_K': 366.15,
        'assumed_relative_k_uncertainty_1995': 0.10,
        'assumed_T_uncertainty_K_1995': 0.5,
        'source': 'https://doi.org/10.1515/pac-2018-1108',
        'original': 'https://pure.tue.nl/ws/portalfiles/portal/1333387/617494.pdf',
        'uncertainty_meaning': '1995 EVM input assumptions; NOT independent A/Ea confidence limits',
    },
    'equilibrium_110C': {
        'T_K': 383.15, 'M_mol_L': 1.2e-4,
        'source': 'https://doi.org/10.1002/(SICI)1099-0518(20000615)38:12<2137::AID-POLA20>3.0.CO;2-D',
        'original': 'https://doi.org/10.1002/pol.1962.1205816633',
        'status': 'Ivin 2000 quotes Bywater/Worsfold 1962; exact solvent assignment not retrieved',
        'phase': 'liquid hydrocarbon solution; original studied benzene and cyclohexane',
        'mechanism': 'living anionic, butyllithium; not a radical reverse-rate measurement',
        'uncertainty': None,
    },
    'equilibrium_25C': {
        'T_K': 298.15, 'M_mol_L': 1e-6,
        'source': 'https://doi.org/10.1002/(SICI)1099-0518(20000615)38:12<2137::AID-POLA20>3.0.CO;2-D',
        'status': 'approximate literature estimate in review; not established here as a direct measurement',
        'phase': 'condensed polymer/monomer reference; exact medium not specified in retrieved excerpt',
        'uncertainty': None,
    },
    'ESR_1992': {
        'range_K': [273.15, 403.15], 'Ea_J_mol': 39700.0,
        'source': 'https://doi.org/10.1007/BF00944835',
        'A_L_mol_s': None, 'status': 'publisher abstract; no numerical high-T kp points retrieved',
    },
    'Bywater_Worsfold_1962': {
        'range_K': [373.15, 423.15],
        'source': 'https://doi.org/10.1002/pol.1962.1205816633',
        'retrieved_abstract': 'https://www.researchgate.net/publication/225109668_Thermodynamics_of_polymerization_with_special_emphasis_on_living_polymers',
        'phase': 'benzene and cyclohexane solution; butyllithium initiation',
        'status': 'original abstract mirrored/indexed; full concentration tables not retrieved',
    },
    'EPR_2004': {
        'T_K': 393.15, 'diffusion_onset_conversion_approx': 0.8,
        'source': 'https://doi.org/10.1002/macp.200300148',
        'status': 'publisher abstract; no numerical high-T kp points retrieved',
    },
    'chain_length_2002': {
        'range_K': [298.15, 343.15], 'observed_variation_percent': [25.0, 35.0],
        'extrapolated_reduction_percent': [40.0, 60.0], 'half_change_chain_length_order': 100,
        'source': 'https://doi.org/10.1021/ma011215b',
    },
    'DMF_2000': {
        'T_K': 313.15, 'M_mol_L': 1.0, 'kp_over_bulk_approx': 0.75,
        'source': 'https://doi.org/10.1016/S0014-3057(00)00021-5',
    },
    'reanalysis_2022': {
        'independent_studies': 81, 'sd_ln_k_25C': 0.08,
        'sd_Ea_J_mol': 1400.0, 'correlation': 0.04,
        'source': 'https://doi.org/10.1039/D2PY00147K',
        'status': 'pooled interlaboratory errors from abstract; revised styrene parameters not retrieved',
    },
}
TEMPERATURES = (261.15, 273.15, 298.15, 300.0, 323.15, 350.0,
                363.15, 366.15, 393.15, 403.15, 600.0, 650.0, 700.0, 750.0, 800.0)
SOURCE_FILES = ('rmgpy/kmc/compiler.py', 'rmgpy/kmc/ssa.py', 'rmgpy/kmc/met.py',
                'rmgpy/kmc/event_record.py', 'rmgpy/data/kinetics/family.py',
                'rmgpy/data/kinetics/rules.py', 'rmgpy/reaction.py',
                'test/rmgpy/kmc/fixtures/i034_probe/run_probe.py')


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def close(a, b, relative=2e-10):
    assert math.isclose(a, b, rel_tol=relative, abs_tol=1e-12), (a, b)


def benchmark(t):
    data = LITERATURE['iupac']
    return data['A_L_mol_s'] * math.exp(-data['Ea_J_mol'] / (constants.R * t))


def safe_output(path):
    path = path.resolve()
    if (path.is_relative_to(ROOT) or path.is_relative_to(DATABASE)
            or path.is_relative_to(Path('/home/alon/Code/polymers')) or 'catalog' in path.parts):
        raise ValueError('output must be in external, non-excluded scratch')
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scratch', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    scratch, output = safe_output(args.scratch), safe_output(args.output)
    print('[I044] materializing pinned allowlist; no event compilation', flush=True)
    snapshot = prior.snapshot_database(DATABASE, scratch / 'database')
    db = RMGDatabase()
    db.load_kinetics(str(scratch / 'database/input/kinetics'), reaction_libraries=[],
                     seed_mechanisms=None, kinetics_families=['R_Addition_MultipleBond'],
                     kinetics_depositories=['training'])
    db.load_thermo(str(scratch / 'database/input/thermo'),
                   thermo_libraries=['primaryThermoLibrary'], depository=True)
    artifact = json.loads(ARTIFACT.read_text())
    assert len(artifact['records']) == 14998
    assert artifact['provenance']['rmg_database_sha'] == DATABASE_SHA
    family = db.kinetics.families['R_Addition_MultipleBond']
    rule_list = [entry for entries in family.rules.entries.values() for entry in entries]
    assert len(rule_list) == 1
    default = rule_list[0]
    training = family.get_training_depository()
    pairs = []
    for pair in artifact.get('ps_primary_end_ceiling_pairs', artifact['ps_ceiling_pairs']):
        prop = next(r for r in artifact['records'] if r['event_id'] == pair['propagation_event_id'])
        dep = next(r for r in artifact['records'] if r['event_id'] == pair['depropagation_event_id'])
        reaction = prior.reaction_from_record(prop)
        estimate = TemplateReaction(reactants=reaction.reactants, products=reaction.products,
                                    family=prop['family'], template=prop['template'].split(';'),
                                    degeneracy=prop['raw_path_degeneracy'])
        kinetics, source, matched, direction = family.get_kinetics(
            estimate, estimate.template, degeneracy=estimate.degeneracy, return_all_kinetics=False)
        assert source == 'rate rules' and matched is None
        for species in reaction.reactants + reaction.products:
            species.thermo = db.thermo.get_thermo_data(species)
        kc_table = prop['thermo_provenance']['equilibrium_constant_table']
        for t, kf, kr, kc in zip(prop['k_table']['T'], prop['k_table']['k'],
                                dep['k_table']['k'], kc_table['Kc']):
            close(kinetics.get_rate_coefficient(t), kf)
            close(default.data.get_rate_coefficient(t), kf)
            close(reaction.get_equilibrium_constant(t, type='Kc'), kc)
            close(kf / kc, kr)
        for field in ('degeneracy', 'raw_path_degeneracy', 'ssa_multiplier'):
            assert prop[field] == dep[field] == 1.0
        pairs.append({'pair': pair, 'propagation': prop, 'depropagation': dep,
                      'reactants_smiles': prior.smiles(reaction.reactants),
                      'products_smiles': prior.smiles(reaction.products),
                      'kinetics_repr': repr(kinetics), 'kinetics_comment': kinetics.comment,
                      'grid_points_checked': len(prop['k_table']['T']),
                      'thermo_species': [{'adjacency': s.molecule[0].to_adjacency_list(),
                                         'thermo_comment': s.thermo.comment}
                                        for s in reaction.reactants + reaction.products]})
    first = pairs[0]
    prop, dep = first['propagation'], first['depropagation']
    rxn = prior.reaction_from_record(prop)
    for species in rxn.reactants + rxn.products:
        species.thermo = db.thermo.get_thermo_data(species)
    parameters = {'A_m3_mol_s': default.data.A.value_si, 'n': default.data.n.value_si,
                  'alpha': default.data.alpha.value_si, 'Ea_J_mol': default.data.E0.value_si,
                  'Tmin_K': default.data.Tmin.value_si, 'Tmax_K': default.data.Tmax.value_si,
                  'rule_index': default.index, 'rule_label': default.label,
                  'rank': default.rank, 'short_desc': default.short_desc,
                  'rule_count': len(rule_list), 'training_entry_count': len(training.entries)}
    rows = []
    for t in TEMPERATURES:
        k_rule = default.data.get_rate_coefficient(t) * 1000.0
        kp = benchmark(t)
        try:
            k_ssa = record_rate(prop, t) * 1000.0
        except ValueError:
            k_ssa = None
        a_factor = parameters['A_m3_mol_s'] * 1000.0 / LITERATURE['iupac']['A_L_mol_s']
        ea_factor = math.exp((LITERATURE['iupac']['Ea_J_mol'] - parameters['Ea_J_mol'])
                             / (constants.R * t))
        close(k_rule / kp, a_factor * ea_factor)
        rows.append({'T_K': t, 'rule_L_mol_s': k_rule, 'ssa_L_mol_s': k_ssa,
                     'benchmark_L_mol_s': kp, 'rule_over_benchmark': k_rule / kp,
                     'ssa_over_benchmark': None if k_ssa is None else k_ssa / kp,
                     'A_factor': a_factor, 'Ea_factor': ea_factor,
                     'interpolation_factor': None if k_ssa is None else k_ssa / k_rule})
    reverse = []
    for t in (300.0, 350.0, 383.15, 600.0, 650.0, 700.0, 750.0, 800.0):
        kf, kr = record_rate(prop, t), record_rate(dep, t)
        kc = rxn.get_equilibrium_constant(t, type='Kc')
        reverse.append({'T_K': t, 'kf_L_mol_s': kf * 1000.0, 'kr_s_1': kr,
                        'Kc_m3_mol': kc, 'M_eq_SSA_mol_L': kr / kf / 1000.0,
                        'M_eq_thermo_mol_L': 1.0 / kc / 1000.0,
                        'kr_with_benchmark_same_Kc_s_1': benchmark(t) / 1000.0 / kc,
                        'forward_per_site_at_1M_s_1': kf * 1000.0,
                        'reverse_over_forward_at_1M': kr / (kf * 1000.0)})
    equilibrium = []
    for key in ('equilibrium_110C', 'equilibrium_25C'):
        datum = LITERATURE[key]
        t, measured = datum['T_K'], datum['M_mol_L']
        kc = rxn.get_equilibrium_constant(t, type='Kc')
        ce = 1.0 / kc / 1000.0
        kf_rule = default.data.get_rate_coefficient(t)
        equilibrium.append({'key': key, 'T_K': t, 'literature_M_mol_L': measured,
                            'model_M_mol_L': ce, 'model_over_literature': ce / measured,
                            'delta_G_required_J_mol': constants.R * t * math.log(measured / ce),
                            'model_rule_kr_s_1': kf_rule / kc,
                            'conditional_radical_kr_from_benchmark_s_1': benchmark(t) * measured,
                            'model_rule_over_conditional_kr': (kf_rule / kc) / (benchmark(t) * measured),
                            'ssa_kr_s_1': record_rate(dep, t) if t >= 300 else None})
    root = brentq(lambda t: math.log(rxn.get_equilibrium_constant(t) * 1000.0), 600.0, 800.0)
    # I039's structure-only benzylic control; no family generation or compilation.
    molecules = ('CC(c1ccccc1)C[CH](c1ccccc1)', 'C=Cc1ccccc1',
                 'CC(c1ccccc1)CC(c1ccccc1)C[CH](c1ccccc1)')
    species = [prior.Species(molecule=[prior.Molecule(smiles=s)]) for s in molecules]
    for item in species:
        item.molecule[0].update()
        item.generate_resonance_structures()
        item.thermo = db.thermo.get_thermo_data(item)
    control = prior.Reaction(reactants=species[:2], products=species[2:])
    control_root = brentq(lambda t: math.log(control.get_equilibrium_constant(t) * 1000.0), 600.0, 800.0)
    structure = {'species_smiles': molecules,
                 'species_adjacencies': [s.molecule[0].to_adjacency_list() for s in species],
                 'Tc_K': control_root,
                 'H_change_J_mol_at_298_15': control.get_enthalpy_of_reaction(298.15) - rxn.get_enthalpy_of_reaction(298.15),
                 'S_change_J_mol_K_at_298_15': control.get_entropy_of_reaction(298.15) - rxn.get_entropy_of_reaction(298.15),
                 'Tc_change_K': control_root - root}
    for field in ('propagation', 'depropagation'):
        for a, b in zip(pairs[0][field]['k_table']['k'], pairs[1][field]['k_table']['k']):
            close(a, b)
    result = {'artifact': str(ARTIFACT), 'artifact_sha256': digest(ARTIFACT),
              'artifact_rmgpy_sha': artifact['provenance']['rmgpy_sha'],
              'database_sha': DATABASE_SHA, 'snapshot': snapshot,
              'record_count': len(artifact['records']), 'R_J_mol_K': constants.R,
              'parameters': parameters, 'pairs': pairs, 'literature': LITERATURE,
              'comparison': rows, 'reverse': reverse, 'equilibrium': equilibrium,
              'structure_control': structure,
              'Tc_artifact_grid_K': artifact['ps_ceiling_temperature_K'], 'Tc_thermo_K': root,
              'source_sha256': {p: digest(ROOT / p) for p in SOURCE_FILES}}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print('I044 two pairs: 58 forward rates, 58 Kc values and 58 reverse rates reproduced')
    print('I044 generic default rule recovered; no matched training source or specific styrene rule')
    print('I044 results: ' + str(output))


if __name__ == '__main__':
    main()
