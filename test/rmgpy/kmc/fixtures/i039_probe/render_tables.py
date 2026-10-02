"""Render every numeric report block from the reproduced JSON (no network)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPORT = Path(__file__).resolve().parent.parent / "I039_tc_gap_probe.md"


def table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |",
                       "| " + " | ".join(["---"] * len(headers)) + " |",
                       *["| " + " | ".join(map(str, row)) + " |" for row in rows]])


def blocks(result):
    sizes = result["sizes"]
    gas = sizes["3"]["continuous_K"]
    solv = result["baselines"]["toluene"]
    output = {}
    output["baselines"] = table(
        ["Model at 1 mol/L", "Compiler-grid Tc (K)", "Continuous Tc (K)", "Continuous shift (K)"],
        [["gas", f'{result["baselines"]["gas_grid_K"]:.6f}', f"{gas:.6f}", "0.000000"],
         ["toluene / linear H−TS", f'{solv["tabulated_K"]:.6f}',
          f'{solv["continuous_K"]:.6f}', f'{solv["continuous_K"]-gas:+.6f}']])
    output["sizes"] = table(
        ["L", "Chain reaction", "Grid Tc (K)", "Continuous Tc (K)", "Same-state 1 M Tc (K)",
         "1-bar monomer Tc (K)", "ΔH at anchor, gas (kJ/mol)", "ΔS at anchor, gas (J/mol/K)"],
        [[n, f"P{int(n)-1}• + styrene ⇌ P{n}•", f'{x["rates"]["ceiling_tabulated_K"]:.6f}',
          f'{x["continuous_K"]:.6f}', f'{x["same_state_1M_K"]:.6f}',
          f'{x["one_bar_monomer_root_K"]:.6f}',
          f'{x["reaction_rows"][0]["pressure_HS"][0]/1000:.6f}',
          f'{x["reaction_rows"][0]["pressure_HS"][1]:.6f}'] for n, x in sizes.items()])
    output["species"] = table(
        ["L", "Side", "Species (SMILES)", "σ including optical", "Optical half factors", "Thermo source"],
        [[n, "reactant" if sign < 0 else "product", f'`{s["smiles"]}`',
          f'{s["symmetry_number_including_optical"]:.6f}', s["optical_atom_half_factors"],
          s["thermo_comment"]] for n, x in sizes.items() for sign, s in zip(x["signs"], x["species"])])
    output["thermo"] = table(
        ["L", "T (K)", "ΔH gas, 1 bar (kJ/mol)", "ΔS gas, 1 bar (J/mol/K)",
         "ΔH gas, 1 M (kJ/mol)", "ΔS gas, 1 M (J/mol/K)", "ΔCp gas, 1 bar (J/mol/K)"],
        [[n, f'{row["T_K"]:.2f}', f'{row["pressure_HS"][0]/1000:.6f}',
          f'{row["pressure_HS"][1]:.6f}', f'{row["concentration_HS"][0]/1000:.6f}',
          f'{row["concentration_HS"][1]:.6f}', f'{row["delta_Cp_pressure_J_mol_K"]:.6f}']
         for n, x in sizes.items() for row in x["reaction_rows"]])
    groups = sizes["3"]["reaction_rows"]
    output["groups"] = table(
        ["Net source", "Net weight", "ΔH at anchor (kJ/mol)", "ΔS at anchor (J/mol/K)",
         "ΔH at 700 K (kJ/mol)", "ΔS at 700 K (J/mol/K)"],
        [[label, sum(sign * s["source_weights"].get(label, 0)
                     for sign, s in zip(sizes["3"]["signs"], sizes["3"]["species"]))
          if not label.startswith(("symmetry:", "optical:")) else "applied per species",
          f'{values[0]/1000:.6f}', f'{values[1]:.6f}',
          f'{groups[3]["net_source_terms"][label][0]/1000:.6f}',
          f'{groups[3]["net_source_terms"][label][1]:.6f}']
         for label, values in sorted(groups[0]["net_source_terms"].items())])
    output["solvation"] = table(
        ["Toluene LSER term", "Reaction descriptor / intercept count", "ΔΔH (kJ/mol)",
         "ΔΔS (J/mol/K)", "ΔΔG at anchor 298 K (kJ/mol)"],
        [[label, f'{values["descriptor_delta"]:.8f}', f'{values["ddH_J_mol"]/1000:.6f}',
          f'{values["ddS_J_mol_K"]:.6f}', f'{values["ddG298_J_mol"]/1000:.6f}']
         for label, values in result["solvation"]["terms"].items()]
        + [["Sum", "—", f'{result["solvation"]["ddH_J_mol"]/1000:.6f}',
            f'{result["solvation"]["ddS_J_mol_K"]:.6f}',
            f'{sum(v["ddG298_J_mol"] for v in result["solvation"]["terms"].values())/1000:.6f}']])
    output["solvation_sources"] = table(
        ["Species", "S", "B", "E", "L", "A", "Source"],
        [[f'`{key}`', *[f'{value["descriptors"][field]:.6f}' for field in ("S", "B", "E", "L", "A")],
          value["comment"]] for key, value in result["solvation"]["species"].items()])
    output["radicals"] = table(
        ["Radical", "Matched node", "Lookup result", "Descriptors minus saturated"],
        [[f'`{row["smiles"]}`', row["matched_node"], row.get("correction_error", "correction found"),
          ", ".join(f'{field}={value:+.6f}' for field, value in sorted(row["radical_minus_saturated"].items()))]
         for row in result["solvation"]["radical_audit"]])
    output["solvation_groups"] = table(
        ["Net solute descriptor source", "Weight", "ΔΔH (kJ/mol)", "ΔΔS (J/mol/K)"],
        [[key, row["weight"], f'{row["ddH_J_mol"]/1000:.6f}', f'{row["ddS_J_mol_K"]:.6f}']
         for key, row in sorted(result["solvation"]["descriptor_source_terms"].items())]
        + [["LSER intercept (molecule-count change)", -1,
            f'{result["solvation"]["terms"]["intercept"]["ddH_J_mol"]/1000:.6f}',
            f'{result["solvation"]["terms"]["intercept"]["ddS_J_mol_K"]:.6f}']])
    control = result["chain_end_control"]
    output["chain_end"] = table(
        ["Structural control", "Continuous Tc (K)", "ΔH at anchor, 1 bar (kJ/mol)",
         "ΔS at anchor, 1 bar (J/mol/K)", "Shift from compiled primary-end pair (K)"],
        [["terminal benzylic radical", f'{control["continuous_K"]:.6f}',
          f'{control["pressure_HS_at_anchor"][0]/1000:.6f}',
          f'{control["pressure_HS_at_anchor"][1]:.6f}', f'{control["continuous_K"]-gas:+.6f}']])
    output["corrected_thermo"] = table(
        ["T (K)", "Toluene-corrected ΔH, 1 M (kJ/mol)", "Toluene-corrected ΔS, 1 M (J/mol/K)", "ΔΔGsolv (kJ/mol)"],
        [[f'{row["T_K"]:.2f}', f'{(row["concentration_HS"][0]+result["solvation"]["ddH_J_mol"])/1000:.6f}',
          f'{row["concentration_HS"][1]+result["solvation"]["ddS_J_mol_K"]:.6f}',
          f'{(result["solvation"]["ddH_J_mol"]-row["T_K"]*result["solvation"]["ddS_J_mol_K"])/1000:.6f}']
         for row in sizes["3"]["reaction_rows"]])
    output["solvation_controls"] = table(
        ["Unapplied solvation component control", "Continuous Tc (K)", "Shift from gas (K)"],
        [[label.replace("_", " "), f"{value:.6f}", f"{value-gas:+.6f}"]
         for label, value in result["solvation"]["continuous_controls_K"].items()])
    literature = result["literature"]
    output["literature"] = table(
        ["Reference", "ΔHp (kJ/mol)", "ΔSp (J/mol/K)", "Temperature", "State / uncertainty"],
        [["Roberts et al.: solid PS", f'{literature["Roberts_solid"]["H_J_mol"]/1000:.2f}', "not measured here",
          f'{literature["Roberts_solid"]["T_K"]:.2f} K', "pure liquid styrene → solid PS; published ±0.66 kJ/mol"],
         ["Roberts et al.: solution", f'{literature["Roberts_solution"]["H_J_mol"]/1000:.2f}', "not measured here",
          f'{literature["Roberts_solution"]["T_K"]:.2f} K',
          f'PS in styrene, {literature["Roberts_solution"]["polymer_weight_percent"]:.1f} wt%; published ±0.69 kJ/mol'],
         ["Warfield & Petree", "not a formation enthalpy datum", f'{literature["Warfield"]["S_J_mol_K"]:.5f}',
          f'{literature["Warfield"]["T_K"]:.2f} K',
          f'condensed-state third-law estimate; entropy loss {literature["Warfield"]["entropy_loss_cal_mol_K"]:.2f} cal/mol/K; uncertainty not supplied'],
         ["Odian Table 3-15", f'{literature["Odian"]["H_J_mol"]/1000:.0f}',
          f'{literature["Odian"]["S_J_mol_K"]:.0f}', "not specified in table footnote", "liquid-monomer ΔH; 1 M monomer ΔS; uncertainties not supplied"],
         ["Odian constant-H/S calculation", "—", "—", f'{literature["Odian"]["H_J_mol"]/literature["Odian"]["S_J_mol_K"]:.6f} K',
          "1 M ceiling; calculation, not a measured ceiling"],
         ["Cowie Table 3.6", f'{literature["Cowie"]["H_J_mol"]/1000:.1f}', "not specified", f'{literature["Cowie"]["pure_liquid_Tc_K"]:.0f} K',
          "pure liquid monomer ceiling; uncertainty not supplied"],
         ["Kirk-Othmer liquid-monomer column", f'{literature["Kirk_Othmer"]["liquid_H_J_mol"]/1000:.1f}',
          f'{literature["Kirk_Othmer"]["liquid_S_J_mol_K"]:.1f}', f'{literature["Kirk_Othmer"]["T_K"]:.2f} K',
          f'liquid monomer / condensed polymer reference, {literature["Kirk_Othmer"]["P_Pa"]/1000:.1f} kPa; tabulated Tc={literature["Kirk_Othmer"]["liquid_Tc_C"]:.0f} °C; no error bars'],
         ["Kirk-Othmer gas-monomer column", f'{literature["Kirk_Othmer"]["gas_H_J_mol"]/1000:.1f}',
          f'{literature["Kirk_Othmer"]["gas_S_J_mol_K"]:.1f}', f'{literature["Kirk_Othmer"]["T_K"]:.2f} K',
          f'gas monomer / condensed polymer reference, {literature["Kirk_Othmer"]["P_Pa"]/1000:.1f} kPa; tabulated Tc={literature["Kirk_Othmer"]["gas_Tc_C"]:.0f} °C; no error bars'],
         ["Dispatch's literature comparator", "not supplied", "not supplied", f'{literature["contract_target_K"]:.0f} K',
          "consistent with rounded Kirk-Othmer liquid-monomer value; source identity not supplied in dispatch"]])
    comparison = result["literature_convention_comparison"]
    output["literature_conventions"] = table(
        ["Recomputation / comparison", "Temperature or difference (K)"],
        [["Compiled literature liquid H/S ratio (frozen anchor)", f'{comparison["liquid_constant_HS_K"]:.6f}'],
         ["Compiled literature gas H/S ratio at source pressure", f'{comparison["gas_constant_HS_at_source_pressure_K"]:.6f}'],
         ["Compiled literature gas reference at 1 bar", f'{comparison["gas_constant_HS_at_1bar_K"]:.6f}'],
         ["Compiled literature gas reference, monomer at 1 M", f'{comparison["gas_constant_HS_at_1M_K"]:.6f}'],
         ["Monomer-condition shift in literature gas case: 1 bar → 1 M", f'{comparison["gas_constant_HS_at_1M_K"]-comparison["gas_constant_HS_at_1bar_K"]:+.6f}'],
         ["Model minus compiled literature on gas-monomer 1 M convention", f'{gas-comparison["gas_constant_HS_at_1M_K"]:+.6f}'],
         ["Model monomer-condition shift: 1 bar → 1 M", f'{gas-sizes["3"]["one_bar_monomer_root_K"]:+.6f}']])
    central = result["conditional_accounting"]
    labels = {"convention": "Convention (same physical state)",
              "enthalpy_error_conditional": "Enthalpy error (conditional allocation)",
              "entropy_error_conditional": "Entropy error (conditional allocation)",
              "proxy_size_observed": "Proxy size (L=5 → L=3, observed model shift)",
              "solvation_including_grid_definition": "Solvation (plus grid interpolation)",
              "residual_physics_not_represented_or_unresolved_reference": "Residual: physics not represented / unresolved reference"}
    output["attribution"] = table(
        ["Contribution to compiled toluene Tc minus dispatch comparator", "ΔT (K)",
         "Rounding-only sensitivity (K)", "Scientific uncertainty"],
        [[labels[name], f'{value:+.6f}', f'±{central["rounding_only_half_ranges_K"][name]:.6f}',
          "exact coordinate identity" if name == "convention" else "not identifiable; no confidence interval"]
         for name in labels for value in [central["rows_K"][name]]]
        + [["Sum", f'{central["sum_K"]:+.6f}', "correlated contributions; do not add error bars", "not identifiable"]])
    output["counterfactuals"] = table(
        ["Anchor H", "Anchor S", "Continuous Tc (K)"],
        [["model" if key[0] == "1" else "Odian (conditional anchor)",
          "model" if key[1] == "1" else "Odian (conditional anchor)", f'{value:.6f}']
         for key, value in central["roots_K"].items()])
    output["definitions"] = table(
        ["Quantity", "Value"],
        [["R (J/mol/K)", f'{result["constants"]["R"]:.8f}'],
         ["Gas reference pressure (Pa)", f'{result["constants"]["P0_Pa"]:.0f}'],
         ["Monomer concentration / concentration reference (mol/m³)", f'{result["constants"]["C0_mol_m3"]:.0f}'],
         ["Conditional thermo anchor (K)", f'{result["constants"]["T0_K"]:.2f}'],
         ["Toluene molecular critical temperature (K)", f'{result["solvation"]["critical_temperature_K"]:.6f}'],
         ["Gas grid minus continuous crossing (K)", f'{result["baselines"]["gas_grid_K"]-gas:+.6f}'],
         ["Toluene grid minus continuous crossing (K)", f'{solv["tabulated_K"]-solv["continuous_K"]:+.6f}'],
         ["Model L=4 minus L=5 continuous crossing (K)", f'{sizes["4"]["continuous_K"]-sizes["5"]["continuous_K"]:+.6f}'],
         ["Unapplied optical-term-suppressed control, continuous Tc (K)", f'{result["optical_term_suppressed_control_K"]:.6f}'],
         ["Persisted pinned database files", result["provenance"]["snapshot"]["files"]]])
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("--report", type=Path, default=REPORT)
    parser.add_argument("--update", action="store_true")
    args = parser.parse_args()
    rendered = blocks(json.loads(args.results.read_text()))
    report = args.report.read_text()
    for name, body in rendered.items():
        begin, end = f"<!-- BEGIN I039:{name} -->", f"<!-- END I039:{name} -->"
        assert report.count(begin) == report.count(end) == 1, name
        prefix, rest = report.split(begin)
        actual, suffix = rest.split(end)
        if args.update:
            report = prefix + begin + "\n" + body + "\n" + end + suffix
        else:
            assert actual == "\n" + body + "\n", f"report block differs: {name}"
    if args.update:
        args.report.write_text(report)
        print(f"I039 rendered {len(rendered)} numeric report blocks")
    else:
        print(f"I039 all {len(rendered)} numeric report blocks match reproduced output")


if __name__ == "__main__":
    main()
