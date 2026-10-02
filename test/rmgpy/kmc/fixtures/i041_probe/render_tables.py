"""Render/check all numerical report blocks; command is in the report."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPORT = Path(__file__).resolve().parent.parent / "I041_config_entropy.md"


def table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |",
                       "| " + " | ".join(["---"] * len(headers)) + " |",
                       *["| " + " | ".join(map(str, row)) + " |" for row in rows]])


def blocks(result):
    sizes = result["sizes"]
    output = {}
    output["constants"] = table(["Quantity", "Value"], [
        ["R (J/mol/K)", f'{result["constants"]["R"]:.8f}'],
        ["Gas reference pressure (Pa)", f'{result["constants"]["P0"]:.0f}'],
        ["Monomer concentration (mol/m³)", f'{result["constants"]["C0"]:.0f}'],
        ["Pinned database files", result["provenance"]["snapshot"]["files"]],
        ["Python / RDKit", result["versions"]["python"] + " / " + result["versions"]["rdkit"]]])
    output["species"] = table(
        ["L", "Side", "SMILES", "Atom product", "Bond product", "Axis", "Cyclic",
         "RMG σ", "Half factors m", "σ without optical", "Enumerated stereoisomers g"],
        [[units, "reactant" if sign < 0 else "product", f'`{item["smiles"]}`',
          *[f'{item["symmetry"][key]:.6f}' for key in (
              "atom_product", "bond_product", "axis_factor", "cyclic_factor", "sigma")],
          item["symmetry"]["centres"], f'{item["symmetry"]["sigma_nonoptical"]:.6f}',
          item["stereo"]["count"]]
         for units, size in sizes.items() for sign, item in zip(size["signs"], size["species"])])
    output["species_entropy"] = table(
        ["L=3 side", "S without optical symmetry (J/mol/K)",
         "S optical (J/mol/K)", "Total symmetry S (J/mol/K)"],
        [["reactant chain" if i == 0 else "styrene" if i == 1 else "product chain",
          f'{item["symmetry"]["S_nonoptical"]:.6f}',
          f'{item["symmetry"]["S_optical"]:.6f}',
          f'{item["symmetry"]["S_optical"] + item["symmetry"]["S_nonoptical"]:.6f}']
         for i, item in enumerate(sizes["3"]["species"])])
    output["main"] = table(
        ["T (K)", "Source ΔH, 1 bar (kJ/mol)", "Source ΔS (J/mol/K)",
         "Nonoptical ΔS", "Stereo ΔS", "1 M conversion ΔH (kJ/mol)",
         "1 M conversion ΔS", "Total ΔH, 1 M (kJ/mol)", "Total ΔS, 1 M (J/mol/K)"],
        [[f'{row["T"]:.0f}', f'{row["H_source"]/1000:.6f}', f'{row["S_source"]:.6f}',
          f'{row["S_symmetry"]:.6f}', f'{row["S_stereo"]:.6f}',
          f'{row["H_convention"]/1000:.6f}', f'{row["S_convention"]:.6f}',
          f'{row["H_c"]/1000:.6f}', f'{row["S_c"]:.6f}'] for row in sizes["3"]["rows"]])
    output["optical_energy"] = table(
        ["T (K)", "Stereo ΔH (kJ/mol)", "Stereo −TΔS (kJ/mol)",
         "Kc with / without term", "Reverse with / without term, fixed total forward", "Kc (m³/mol)"],
        [[f'{row["T"]:.0f}', "0.000000", f'{row["G_stereo"]/1000:.6f}',
          f'{row["Kc_optical_multiplier"]:.6f}',
          f'{1/row["Kc_optical_multiplier"]:.6f}', f'{row["Kc"]:.9g}']
         for row in sizes["3"]["rows"]])
    nominal = sizes["3"]["cases"][0]["Tc_continuous"]
    output["cases"] = table(
        ["Convention/control", "ΔΔH (kJ/mol)", "ΔΔS (J/mol/K)", "Kc / baseline",
         "Reverse / baseline at fixed forward", "Continuous gas Tc (K)",
         "Grid gas Tc (K)", "Continuous ΔTc (K)"],
        [[item["name"], f'{item["shift_H"]/1000:.6f}', f'{item["shift_S"]:.6f}',
          f'{item["K_multiplier"]:.6f}', f'{item["reverse_multiplier_at_fixed_total_forward"]:.6f}',
          f'{item["Tc_continuous"]:.6f}', f'{item["Tc_grid"]:.6f}',
          f'{item["Tc_continuous"]-nominal:+.6f}'] for item in sizes["3"]["cases"]])
    output["sizes"] = table(
        ["Product L", "Net half factors", "Nonoptical ΔS (J/mol/K)",
         "Stereo ΔS (J/mol/K)", "Continuous gas Tc (K)", "Grid gas Tc (K)"],
        [[units, item["delta_centres"], f'{item["delta_S_nonoptical"]:.6f}',
          f'{item["delta_S_stereo"]:.6f}', f'{item["cases"][0]["Tc_continuous"]:.6f}',
          f'{item["cases"][0]["Tc_grid"]:.6f}'] for units, item in sizes.items()]
        + [["benzylic end control", result["benzylic_control"]["delta_centres"],
            f'{result["benzylic_control"]["delta_S_nonoptical"]:.6f}',
            f'{result["benzylic_control"]["delta_S_stereo"]:.6f}',
            f'{result["benzylic_control"]["Tc_continuous"]:.6f}', "not compiled"]])
    output["families"] = table(
        ["Selected forward reaction", "Constitutional SMILES: reactants → products", "Δm",
         "RMG stereo ΔS (J/mol/K)", "RMG stereo K multiplier", "Enumerated g-product/g-reactant",
         "Nonoptical ΔS (J/mol/K)", "Full gas ΔH at 700 K (kJ/mol)", "Full gas ΔS at 700 K (J/mol/K)"],
        [[name, " + ".join(f'`{s["smiles"]}`' for sign, s in zip(item["signs"], item["species"]) if sign < 0)
          + " → " + " + ".join(f'`{s["smiles"]}`' for sign, s in zip(item["signs"], item["species"]) if sign > 0),
          item["delta_centres"], f'{item["delta_S_stereo"]:.6f}',
          f'{item["rows"][1]["Kc_optical_multiplier"]:.6f}',
          f'{stereo_ratio(item):.6f}', f'{item["delta_S_nonoptical"]:.6f}',
          f'{item["rows"][1]["H_p"]/1000:.6f}', f'{item["rows"][1]["S_p"]:.6f}']
         for name, item in result["families"].items()])
    output["general"] = table(
        ["Δm", "Stereo ΔS (J/mol/K)", "K multiplier", "−TΔS at 600 K (kJ/mol)",
         "at 700 K", "at 800 K"],
        [[row["delta_centres"], f'{row["S"]:.6f}', f'{row["K_multiplier"]:.6f}',
          *[f'{g/1000:.6f}' for g in row["G_at_T"]]] for row in result["general_stereo_terms"]])
    return output


def stereo_ratio(item):
    ratio = 1.0
    for sign, species in zip(item["signs"], item["species"]):
        ratio *= species["stereo"]["count"]**sign
    return ratio


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("--report", type=Path, default=REPORT)
    parser.add_argument("--update", action="store_true")
    args = parser.parse_args()
    rendered = blocks(json.loads(args.results.read_text()))
    report = args.report.read_text()
    for name, body in rendered.items():
        begin, end = f"<!-- BEGIN I041:{name} -->", f"<!-- END I041:{name} -->"
        assert report.count(begin) == report.count(end) == 1, name
        prefix, rest = report.split(begin)
        actual, suffix = rest.split(end)
        if args.update:
            report = prefix + begin + "\n" + body + "\n" + end + suffix
        else:
            assert actual == "\n" + body + "\n", f"report block differs: {name}"
    if args.update:
        args.report.write_text(report)
        print(f"I041 rendered {len(rendered)} numeric report blocks")
    else:
        print(f"I041 all {len(rendered)} numeric report blocks match reproduced output")


if __name__ == "__main__":
    main()
