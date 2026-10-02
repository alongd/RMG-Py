"""Generate or check the report's numerical blocks; command in the report."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

REPORT = Path(__file__).resolve().parent.parent / 'I044_kp_benchmark.md'


def table(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |',
                      '| ' + ' | '.join(['---'] * len(headers)) + ' |',
                      *['| ' + ' | '.join(map(str, row)) + ' |' for row in rows]])


def number(value):
    return 'refused' if value is None else f'{value:.6g}'


def blocks(result):
    p, lit = result['parameters'], result['literature']
    blocks = {}
    blocks['parameters'] = table(
        ['Quantity', 'Recovered compiled source', 'Common IUPAC bulk benchmark'],
        [['A (L/mol/s)', number(p['A_m3_mol_s'] * 1000), number(lit['iupac']['A_L_mol_s'])],
         ['Ea (kJ/mol)', f"{p['Ea_J_mol']/1000:.3f}", f"{lit['iupac']['Ea_J_mol']/1000:.3f}"],
         ['Temperature exponent n', number(p['n']), '0'],
         ['Source range (K)', f"{p['Tmin_K']:.2f}–{p['Tmax_K']:.2f}",
          '–'.join(f'{t:.2f}' for t in lit['iupac']['range_K'])],
         ['Source rank / matched training entry', f"{p['rank']} / none", 'PLP-SEC multi-laboratory benchmark'],
         ['Loaded addition rules / training entries', f"{p['rule_count']} / {p['training_entry_count']}", 'not applicable']])
    blocks['provenance'] = table(
        ['Quantity', 'Reproduced value'],
        [['Artifact records', result['record_count']], ['Artifact SHA256', result['artifact_sha256']],
         ['Pinned snapshot files', result['snapshot']['files']],
         ['Pinned snapshot SHA256', result['snapshot']['sha256']],
         ['R (J/mol/K)', f"{result['R_J_mol_K']:.8f}"],
         ['Artifact-grid crossing at 1 M (K)', f"{result['Tc_artifact_grid_K']:.6f}"],
         ['Continuous thermo crossing at 1 M (K)', f"{result['Tc_thermo_K']:.6f}"]])
    blocks['uncertainty'] = table(
        ['Source', 'Stated uncertainty', 'Meaning'],
        [['Original IUPAC EVM fit (1995)',
          f"±{lit['iupac']['assumed_T_uncertainty_K_1995']:.1f} K, ±{100*lit['iupac']['assumed_relative_k_uncertainty_1995']:.0f}% k_p; 95% joint A/Ea region",
          'Assumed input measurement errors; A and Ea are strongly correlated'],
         ['Interlaboratory reanalysis (2022)',
          f"SD ln(k_p at 25°C) = {lit['reanalysis_2022']['sd_ln_k_25C']:.2f}; SD Ea = {lit['reanalysis_2022']['sd_Ea_J_mol']/1000:.1f} kJ/mol; correlation = {lit['reanalysis_2022']['correlation']:.2f}; {lit['reanalysis_2022']['independent_studies']} studies",
          'Pooled error across monomers; not a styrene-specific extrapolation confidence band']])
    for name, chosen in [('measured', result['comparison'][:8]), ('high_temperature', result['comparison'][8:])]:
        blocks[name] = table(
            ['T (K)', 'Recovered rule (L/mol/s)', 'SSA table (L/mol/s)',
             'IUPAC Arrhenius (L/mol/s)', 'Rule / IUPAC', 'SSA / IUPAC'],
            [[f"{row['T_K']:.2f}", *[number(row[k]) for k in
              ('rule_L_mol_s', 'ssa_L_mol_s', 'benchmark_L_mol_s',
               'rule_over_benchmark', 'ssa_over_benchmark')]] for row in chosen])
    blocks['contributions'] = table(
        ['T (K)', 'A-factor ratio', 'Ea exponential ratio', 'SSA interpolation / rule'],
        [[f"{row['T_K']:.2f}", number(row['A_factor']), number(row['Ea_factor']),
          number(row['interpolation_factor'])] for row in result['comparison']
         if row['T_K'] in (261.15, 300.0, 363.15, 600.0, 700.0, 800.0)])
    blocks['reverse'] = table(
        ['T (K)', 'Compiled k_rev (s⁻¹)', 'Gas Kc (m³/mol)',
         '[M]eq from thermo (mol/L)', '[M]eq from SSA (mol/L)',
         'IUPAC k_p / same gas Kc (s⁻¹)'],
        [[f"{row['T_K']:.2f}", *[number(row[k]) for k in
          ('kr_s_1', 'Kc_m3_mol', 'M_eq_thermo_mol_L', 'M_eq_SSA_mol_L',
           'kr_with_benchmark_same_Kc_s_1')]] for row in result['reverse']])
    blocks['equilibrium'] = table(
        ['T (K)', 'Literature [M]eq (mol/L)', 'Gas-model [M]eq (mol/L)',
         'Model / literature', 'Conditional ΔG addition (kJ/mol)'],
        [[f"{row['T_K']:.2f}", number(row['literature_M_mol_L']),
          number(row['model_M_mol_L']), number(row['model_over_literature']),
          f"{row['delta_G_required_J_mol']/1000:+.6f}"] for row in result['equilibrium']])
    blocks['conditional_reverse'] = table(
        ['T (K)', 'Model rule k_rev (s⁻¹)', 'SSA k_rev (s⁻¹)',
         'Conditional k_p,IUPAC × literature [M]eq (s⁻¹)', 'Model rule / conditional'],
        [[f"{row['T_K']:.2f}", *[number(row[k]) for k in
          ('model_rule_kr_s_1', 'ssa_kr_s_1', 'conditional_radical_kr_from_benchmark_s_1',
           'model_rule_over_conditional_kr')]] for row in result['equilibrium']])
    blocks['other_evidence'] = table(
        ['Study', 'Retrieved numerical evidence', 'Retrieval / interpretation'],
        [['Yamada et al. 1992 (ESR)',
          f"{lit['ESR_1992']['range_K'][0]:.2f}–{lit['ESR_1992']['range_K'][1]:.2f} K; Ea = {lit['ESR_1992']['Ea_J_mol']/1000:.1f} kJ/mol",
          'Publisher abstract; A and high-T rate values not retrieved'],
         ['Zetterlund et al. 2004 (EPR)', f"{lit['EPR_2004']['T_K']:.2f} K; diffusion effects near {100*lit['EPR_2004']['diffusion_onset_conversion_approx']:.0f}% conversion",
          'Publisher abstract; no digitized rate values'],
         ['Olaj et al. 2002 (PLP chain length)',
          f"{lit['chain_length_2002']['range_K'][0]:.2f}–{lit['chain_length_2002']['range_K'][1]:.2f} K; observed variation {lit['chain_length_2002']['observed_variation_percent'][0]:.0f}–{lit['chain_length_2002']['observed_variation_percent'][1]:.0f}%; extrapolated reduction {lit['chain_length_2002']['extrapolated_reduction_percent'][0]:.0f}–{lit['chain_length_2002']['extrapolated_reduction_percent'][1]:.0f}%; half-change chain length order {lit['chain_length_2002']['half_change_chain_length_order']}",
          'Abstract; short/infinite-chain limits depend on modeling function'],
         ['DMF solution PLP (2000)', f"{lit['DMF_2000']['T_K']:.2f} K, {lit['DMF_2000']['M_mol_L']:.2f} M; k_p ≈ {lit['DMF_2000']['kp_over_bulk_approx']:.2f} × bulk",
          'Acetonitrile also reduced k_p, with a minimum at intermediate dilution']])
    structure = result['structure_control']
    blocks['structure_control'] = table(
        ['I039 benzylic end control', 'Reproduced value'],
        [['Continuous Tc at 1 M (K)', f"{structure['Tc_K']:.6f}"],
         ['ΔH change at 298.15 K (J/mol)', f"{structure['H_change_J_mol_at_298_15']:+.6f}"],
         ['ΔS change at 298.15 K (J/mol/K)', f"{structure['S_change_J_mol_K_at_298_15']:+.6f}"],
         ['Tc change from compiled primary end (K)', f"{structure['Tc_change_K']:+.6f}"]])
    return blocks


def check_report(result, update=False):
    text = REPORT.read_text()
    for label, block in blocks(result).items():
        begin, end = f'<!-- BEGIN I044:{label} -->', f'<!-- END I044:{label} -->'
        assert text.count(begin) == text.count(end) == 1, label
        start, finish = text.index(begin) + len(begin), text.index(end)
        if update:
            text = text[:start] + '\n' + block + '\n' + text[finish:]
        else:
            assert text[start:finish].strip() == block, f'numeric report block differs: {label}'
    if update:
        REPORT.write_text(text)
    print(f"I044 all {len(blocks(result))} numeric report blocks {'rendered' if update else 'verified'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    parser.add_argument('--update', action='store_true')
    args = parser.parse_args()
    check_report(json.loads(args.results.read_text()), update=args.update)


if __name__ == '__main__':
    main()
