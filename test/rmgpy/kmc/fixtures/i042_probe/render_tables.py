"""Render/check every numeric report table from the saved probe results.

Command: PYTHONPATH=$PWD /home/alon/anaconda3/envs/rmg_env/bin/python
test/rmgpy/kmc/fixtures/i042_probe/render_tables.py
/home/alon/runs/i042-melt-reference/reproduce/results.json
Add --update only when authoring the report.
"""

from __future__ import annotations
import argparse
import json
from pathlib import Path
import models as m

REPORT = Path(__file__).resolve().parent.parent / "I042_melt_reference.md"


def number(value, places=3):
    if value is None:
        return "unavailable"
    if abs(value) < 0.5*10**(-places):
        value = 0.0
    return f"{value:.{places}f}"


def table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"]*len(headers)) + " |",
                      *("| " + " | ".join(map(str, row)) + " |" for row in rows)])


def blocks(r):
    out = {}
    d = r["definitions"]
    out["definitions"] = table(["Quantity", "Value"], [
        ["Gas constant (J/mol/K)", number(d["R"], 6)], ["Avogadro constant (1/mol)", f'{d["NA"]:.8e}'],
        ["Styrene-equivalent molar mass (g/mol)", number(d["MW_g_mol"], 6)],
        ["Fixed total mass density (kg/m³)", number(d["density_kg_m3"])],
        ["Conserved repeat-equivalent concentration (mol/L)", number(d["repeat_c_mol_L"], 6)],
        ["Volume per repeat equivalent (m³)", f'{d["repeat_volume_m3"]:.8e}'],
        ["Illustrative carrier DP", number(d["DP"], 6)],
        ["Concentration reference (mol/m³)", number(d["c0_mol_m3"])],
        ["Short / long radical surrogate masses (g/mol)", " / ".join(number(x, 6) for x in d["radical_masses_g_mol"])],
        ["Pinned database files", r["provenance"]["snapshot"]["files"]],
        ["RMG molecular solvent entries", r["solvent_inventory_count"]],
        ["Polystyrene-labelled solvent entries", len(r["polymer_labels"])]])
    main_rows = [["Gas baseline", "0 / 0 / 0", "0", "0", number(r["gas"]["continuous_K"], 6), "0", "gas reference"]]
    for label, pc in r["PC"].items():
        mid = pc["rows"][1]
        main_rows.append([label, " / ".join(number(x["G"]/1000) for x in pc["rows"]),
                          number(mid["H_path"]/1000), number(mid["S_path"]),
                          number(pc["roots_K"][0], 6), number(pc["shifts_K"][0], 6), "radical / T extrapolation"])
    for label, v in r["solvents"].items():
        if "H" not in v:
            continue
        main_rows.append(["RMG linear " + label, " / ".join(number(x["G"]/1000) for x in v["rows"]),
                          number(v["H"]/1000), number(v["S"]), number(v["roots_K"][0], 6),
                          number(v["shifts_K"][0], 6), "298 K anchor extrapolation"])
    out["main"] = table(["Model at 1 M", "ΔΔG at 600 / 700 / 800 K (kJ/mol)", "ΔΔH at 700 K (kJ/mol)",
                         "ΔΔS at 700 K (J/mol/K)", "Continuous Tc (K)", "Shift from gas (K)", "Status"], main_rows)
    out["gas"] = table(["T (K)", "Gas ΔH, 1 M (kJ/mol)", "Gas ΔS, 1 M (J/mol/K)", "Gas ΔG, 1 M (kJ/mol)",
                        "Gas equilibrium monomer (mol/L)"],
                       [[number(v["T_K"], 0), number(v["H"]/1000), number(v["S"]), number(v["G"]/1000),
                         number(e["gas_equilibrium_c_mol_L"], 6)] for v, e in zip(r["gas"]["rows"], r["equilibrium"])])
    out["pc_parameters"] = table(["Set", "m(styrene)", "σ(styrene), Å", "ε/k(styrene), K", "m(PS)/M, mol/g",
                                   "σ(PS), Å", "ε/k(PS), K", "k(styrene,PS)"],
                                  [[label, *(number(v[key], 6) for key in ("monomer_m", "monomer_sigma", "monomer_epsilon",
                                     "polymer_m_per_mass", "polymer_sigma", "polymer_epsilon", "kij"))]
                                   for label, p in r["PC"].items() for v in [p["parameters"]]])
    pc_rows, corrected, isobar, controls = [], [], [], []
    for label, p in r["PC"].items():
        for v in p["rows"]:
            pc_rows.append([label, number(v["T_K"], 0), number(v["G"]/1000), number(v["H_path"]/1000), number(v["S_path"]),
                            number(v["H_pressure"]/1000), number(v["S_pressure"]), number(v["P_Pa"]/1e6),
                            number(v["monomer_concentration_activity"], 6), number(v["pair_equilibrium_multiplier"], 6)])
            corrected.append([label, number(v["T_K"], 0), number(v["corrected_H_path"]/1000), number(v["corrected_S_path"]),
                              number(v["corrected_H_pressure"]/1000), number(v["corrected_S_pressure"])])
        for v in p["isobar_1bar"]:
            isobar.append([label, number(v["T_K"], 0), number(v["density_kg_m3"]), number(v["G"]/1000),
                           number(v["H_pressure"]/1000), number(v["S_pressure"])])
        for v in p["controls"]:
            controls.append([label, number(v["density_kg_m3"], 0), number(v["DP"], 3), number(v["G700"]/1000),
                             ", ".join(number(x, 6) for x in v["roots_K"]) or "no root in scan"])
    out["pc_transfer"] = table(["Set", "T (K)", "ΔΔG (kJ/mol)", "H along isochore (kJ/mol)", "S along isochore (J/mol/K)",
                                "H along isobar (kJ/mol)", "S along isobar (J/mol/K)", "EOS P (MPa)",
                                "Monomer γ relative to ideal gas at equal c", "Kc multiplier for pair"], pc_rows)
    out["corrected"] = table(["Set", "T (K)", "Effective ΔH, isochore (kJ/mol)", "Effective ΔS, isochore (J/mol/K)",
                              "Effective ΔH, isobar (kJ/mol)", "Effective ΔS, isobar (J/mol/K)"], corrected)
    out["isobar"] = table(["Set", "T (K)", "EOS density at 1 bar (kg/m³)", "ΔΔG (kJ/mol)",
                           "Isobar ΔΔH (kJ/mol)", "Isobar ΔΔS (J/mol/K)"], isobar)
    out["controls"] = table(["Set", "Fixed density (kg/m³)", "Carrier DP", "ΔΔG at 700 K (kJ/mol)", "Tc (K)"], controls)
    out["concentrations"] = table(["Free fraction of conserved equivalents", "Monomer (mol/L)", "Gas Tc (K)", "Gas shift (K)",
                                   "Effective ΔS concentration shift (J/mol/K)", "FH χ=0 activity", "UNIFAC activity at 700 K"],
                                  [[number(v["fraction_free_styrene"], 6), number(v["c_mol_L"], 6), number(v["gas_Tc_K"], 6),
                                    number(v["gas_shift_K"], 6), number(v["concentration_S_shift"]), number(v["FH_athermal_activity"], 6),
                                    number(v["UNIFAC_activity700"], 6)] for v in r["concentrations"]])
    mix_rows = []
    for label, rows in (("FH χ=0 control", r["FH_athermal"]), ("Original UNIFAC", r["UNIFAC"])):
        for v in rows:
            mix_rows.append([label, number(v["T_K"], 0), number(v["G_mixing_1M"]/1000), number(v["H_mixing_1M"]/1000),
                             number(v["S_mixing_1M"]), number(v["monomer_activity"], 6),
                             number(v.get("effective_chi", 0), 6), "reference transfer missing"])
    out["mixing"] = table(["Mixing model", "T (K)", "ΔmixG, 1 M (kJ/mol)", "ΔmixH (kJ/mol)", "ΔmixS (J/mol/K)",
                           "Pure-monomer activity", "χ equivalent", "Total transfer / Tc"], mix_rows)
    out["fh_sensitivity"] = table(["T (K)", "∂ΔmixG/∂χ (kJ/mol)", "∂ΔmixH/∂Bχ (J/mol/K)"],
                                  [[number(v["T_K"], 0), number(v["dG_dchi_J_mol"]/1000), number(v["dH_dchi_b_J_mol_K"])]
                                   for v in r["FH_athermal"]])
    rm_rows = []
    for label, v in r["solvents"].items():
        roots = v.get("kfactor_roots_K")
        if "kfactor" not in v:
            rm_rows.append([label, number(v["critical_K"]), "all", "missing enthalpy coefficients", "unavailable"])
            continue
        for row in v["kfactor"]:
            if "error" in row:
                desc = "above/equal solvent critical T" if "critical" in row["error"] else "no CoolProp fluid name"
            else:
                desc = " / ".join(number(row[k]/(1000 if k != "S" else 1)) for k in ("G", "H", "S"))
            rm_rows.append([label, number(v["critical_K"]), number(row["T_K"], 0), desc,
                            "unavailable" if roots is None else (", ".join(number(x, 6) for x in roots) or "no root in liquid domain")])
    out["rmg_validity"] = table(["Solvent", "Molecular critical T (K)", "T (K)", "K-factor G / H / S (kJ/mol / kJ/mol / J/mol/K) or failure",
                                 "K-factor Tc (K)"], rm_rows)
    fields = ("s_g", "b_g", "e_g", "l_g", "a_g", "c_g", "s_h", "b_h", "e_h", "l_h", "a_h", "c_h")
    out["rmg_parameters"] = table(["Solvent", *fields], [[label, *(number(v["coefficients"][f], 6) for f in fields)]
                                                           for label, v in r["solvents"].items()])
    out["unifac_parameters"] = table(["Subgroup", "R", "Q", "Main group"],
                                     [[group, number(rg, 4), number(qg, 4), main+1]
                                      for group, rg, qg, main in zip(m.SUBGROUPS, m.UR, m.UQ, m.MAIN)])
    out["unifac_interactions"] = table(["Main group i", "a(i,1), K", "a(i,2), K", "a(i,3), K", "a(i,4), K"],
                                       [[i+1, *(number(x, 3) for x in row)] for i, row in enumerate(m.INTERACTION)])
    return out


def verify_or_update(results, report=REPORT, update=False):
    text = report.read_text()
    values = blocks(results)
    for name, expected in values.items():
        start, end = f"<!-- BEGIN I042:{name} -->", f"<!-- END I042:{name} -->"
        before, rest = text.split(start)
        current, after = rest.split(end)
        if update:
            text = before+start+"\n"+expected+"\n"+end+after
        elif current.strip() != expected:
            raise AssertionError("report table differs: " + name)
    if update:
        report.write_text(text)
    print(f"I042 all {len(values)} numeric report blocks " + ("updated" if update else "verified"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("--report", type=Path, default=REPORT)
    parser.add_argument("--update", action="store_true")
    args = parser.parse_args()
    verify_or_update(json.loads(args.results.read_text()), args.report, args.update)


if __name__ == "__main__":
    main()
