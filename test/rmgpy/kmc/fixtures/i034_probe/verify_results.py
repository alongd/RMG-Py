"""Independently check the probe and every generated numeric block in its report.

Command: PYTHONPATH=$PWD /home/alon/anaconda3/envs/rmg_env/bin/python
test/rmgpy/kmc/fixtures/i034_probe/verify_results.py /tmp/i034-reproduce/results.json
Use --snapshot when run_probe.py used a nondefault scratch directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from CoolProp.CoolProp import PropsSI
from scipy.optimize import brentq

from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.solvation import SoluteData
from rmgpy.molecule.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species


REPORT = Path(__file__).resolve().parent.parent / "I034_melt_standin_probe.md"
REACTION_LABELS = {"homolysis": "homolysis", "H_baseline": "H (prior styrene pair)",
                   "H_chain": "H (chain-to-chain)", "propagation": "propagation"}


def close(actual, expected, relative=1e-9, absolute=1e-8):
    if not math.isclose(actual, expected, rel_tol=relative, abs_tol=absolute):
        raise AssertionError(f"{actual} != {expected}")


def delta(reaction, values):
    return sum(values[key] for key in reaction["products"]) - sum(values[key] for key in reaction["reactants"])


def lser(descriptors, coefficients):
    terms = (("S", "s"), ("B", "b"), ("E", "e"), ("L", "l"), ("A", "a"))
    log_partition = coefficients["c_g"] + sum(descriptors[field] * coefficients[f"{prefix}_g"] for field, prefix in terms)
    gibbs = -8.314 * 298.0 * 2.303 * log_partition
    enthalpy = 1000.0 * (coefficients["c_h"] + sum(descriptors[field] * coefficients[f"{prefix}_h"] for field, prefix in terms))
    return gibbs, enthalpy, (enthalpy - gibbs) / 298.0


def independent_kfactor(anchor, solvent, temperature, gas_constant):
    critical = PropsSI("T_critical", solvent)
    transition = 0.75 * critical
    critical_density = PropsSI("rhomolar_critical", solvent)
    def liquid(temp):
        return PropsSI("Dmolar", "T", temp, "Q", 0, solvent)
    def vapor(temp):
        return PropsSI("Dmolar", "T", temp, "Q", 1, solvent)
    def reduced_log(temp):
        gibbs = anchor["H298_J_mol"] - temp * anchor["S298_J_mol_K"]
        return (gibbs / (gas_constant * temp) + math.log(liquid(temp) / vapor(temp))) * temp / critical
    def basis(temp):
        reduced = temp / critical
        return [1.0, (1.0 - reduced) ** 0.355, math.exp(1.0 - reduced) * reduced ** 0.59]
    def derivative(temp):
        reduced = temp / critical
        return [0.0, -0.355 / critical * (1.0 - reduced) ** -0.645,
                math.exp(1.0 - reduced) / critical * (0.59 * reduced ** -0.41 - reduced ** 0.59)]
    matrix = np.array([
        [*basis(298.0), 0.0], [*derivative(298.0), 0.0],
        [*basis(transition), -(liquid(transition) / critical_density - 1.0)],
        [*derivative(transition), -(liquid(transition + 1.0) - liquid(transition)) / critical_density],
    ])
    rhs = np.array([reduced_log(298.0), reduced_log(299.0) - reduced_log(298.0), 0.0, 0.0])
    parameters = np.linalg.solve(matrix, rhs)
    if temperature < transition:
        log_factor = np.dot(parameters[:3], basis(temperature)) * critical / temperature
    else:
        log_factor = parameters[3] * (liquid(temperature) / critical_density - 1.0) * critical / temperature
    return gas_constant * temperature * (log_factor + math.log(vapor(temperature) / liquid(temperature)))


def verify(result, snapshot):
    root = Path(__file__).resolve().parents[5]
    assert result["provenance"]["database_sha"] == "4a12d36fcdc193ede82c8d1ab5c1653495d445bc"
    assert result["provenance"]["repository_base_sha"] == "ff8e75ce48b0ff4e5ce8a3a28b985932f784149c"
    for relative, expected in result["provenance"]["code_sha256"].items():
        assert hashlib.sha256((root / relative).read_bytes()).hexdigest() == expected
    close(result["rates"]["Kc_association"][1], 1.63104e14, relative=6e-6)
    close(result["rates"]["homolysis_s-1"][1], 2.02492e-9, relative=6e-6, absolute=1e-20)
    close(result["rates"]["H_Kc"][1], 1.25830e-4, relative=6e-6)
    close(result["rates"]["ceiling_tabulated_K"], 710.2487, absolute=5e-5)
    assert result["solute_failures"] == []
    database = RMGDatabase()
    database.load_thermo(str(snapshot / "input/thermo"), thermo_libraries=["primaryThermoLibrary"], depository=True)
    database.load_solvation(str(snapshot / "input/solvation"))
    solvent_entries = database.solvation.libraries["solvent"].entries
    assert {row["label"] for row in result["inventory"]["rows"]} == set(solvent_entries)
    for row in result["inventory"]["rows"]:
        data = solvent_entries[row["label"]].data
        for suffix in ("g", "h"):
            fields = [f"{prefix}_{suffix}" for prefix in ("s", "b", "e", "l", "a", "c")]
            assert row[f"missing_{suffix}"] == [field for field in fields if getattr(data, field) is None]
        assert row["coolprop"] == data.name_in_coolprop
        if row["coolprop"]:
            close(row["Tc_K"], PropsSI("T_critical", row["coolprop"]))
    species = {}
    for key, entry in result["species"].items():
        item = Species(molecule=[Molecule().from_adjacency_list(entry["adjacency"])])
        item.thermo = database.thermo.get_thermo_data(item)
        for force_groups, field in ((False, "descriptors"), (True, "group_descriptors")):
            copy = item.copy(deep=True)
            if force_groups:
                solute = database.solvation.get_solute_data_from_groups(copy)
                solute.set_mcgowan_volume(copy)
            else:
                solute = database.solvation.get_solute_data(copy)
            for descriptor, expected in entry[field].items():
                close(getattr(solute, descriptor), expected)
        species[key] = item
    for audit in result["radical_audit"]:
        molecule = species[audit["smiles"]].molecule[0]
        radical_atom = next(atom for atom in molecule.atoms if atom.radical_electrons)
        node = database.solvation.groups["radical"].descend_tree(molecule, {"*": radical_atom}, None)
        assert node.label == audit["matched_node"] == "R_rad"
        try:
            database.solvation._add_group_solute_data(SoluteData(S=0.0, B=0.0, E=0.0, L=0.0, A=0.0),
                database.solvation.groups["radical"], molecule, {"*": radical_atom})
        except KeyError as error:
            assert str(error) == audit["correction_error"]
        else:
            raise AssertionError("expected absent carbon-radical correction")
        saturated = molecule.copy(deep=True)
        saturated.saturate_radicals()
        stable = database.solvation.get_solute_data(Species(molecule=[saturated]))
        for field, expected in audit["radical_minus_saturated"].items():
            close(result["species"][audit["smiles"]]["descriptors"][field] - getattr(stable, field), expected)
            if field != "V":
                close(expected, 0.0, absolute=1e-12)
    reactions = {name: Reaction(reactants=[species[key] for key in pair["reactants"]],
                                products=[species[key] for key in pair["products"]])
                 for name, pair in result["reactions"].items()}
    gas_constant = result["provenance"]["R_J_mol_K"]
    checked = 0
    for index, row in enumerate(result["gas"]):
        temperature = row["T_K"]
        for name, values in row["reactions"].items():
            close(reactions[name].get_equilibrium_constant(temperature, type="Kc"), values["Kc"])
            close(reactions[name].get_free_energy_of_reaction(temperature), values["dG_J_mol"])
            expected = math.exp(-values["dG_J_mol"] / (gas_constant * temperature))
            expected *= (1e5 / (gas_constant * temperature)) ** result["reactions"][name]["delta_n"]
            close(values["Kc"], expected)
        close(1.0 / row["reactions"]["homolysis"]["Kc"], result["rates"]["Kc_association"][index])
        close(row["reactions"]["H_baseline"]["Kc"], result["rates"]["H_Kc"][index])
        close(result["rates"]["H_forward_m3_mol_s"][index] / result["rates"]["H_Kc"][index],
              result["rates"]["H_reverse_m3_mol_s"][index])
    for label, candidate in result["candidates"].items():
        source = database.solvation.get_solvent_data(label)
        for field, expected in candidate["coefficients"].items():
            actual = getattr(source, field)
            assert (actual is None) == (expected is None)
            if expected is not None:
                close(actual, expected)
        if not candidate["missing"]:
            for key, entry in result["species"].items():
                gibbs, enthalpy, entropy = lser(entry["descriptors"], candidate["coefficients"])
                anchor = candidate["anchors"][key]
                close(anchor["G298_J_mol"], gibbs)
                close(anchor["H298_J_mol"], enthalpy)
                close(anchor["S298_J_mol_K"], entropy)
        for model, rows in candidate["models"].items():
            for row in rows:
                if "error" in row:
                    if candidate["missing"] or not candidate["coolprop"]:
                        assert row["error_type"] == "DatabaseError"
                    else:
                        assert row["T_K"] >= candidate["Tc_K"]
                        assert row["error_type"] == "InputError"
                    continue
                temperature = row["T_K"]
                fresh = {}
                for key, anchor in candidate["anchors"].items():
                    if model == "linear_HS":
                        fresh[key] = anchor["H298_J_mol"] - temperature * anchor["S298_J_mol_K"]
                    else:
                        fresh[key] = independent_kfactor(anchor, candidate["coolprop"], temperature, gas_constant)
                    close(fresh[key], row["species_solvation_J_mol"][key], relative=1e-8, absolute=1e-5)
                for name, values in row["reactions"].items():
                    close(values["ddG_J_mol"], delta(result["reactions"][name], fresh), relative=1e-8, absolute=1e-5)
                    close(values["ratio"], math.exp(-values["ddG_J_mol"] / (gas_constant * temperature)))
                    close(values["melt_Kc"], values["gas_Kc"] * values["ratio"])
                    close(values["gas_Kc"], reactions[name].get_equilibrium_constant(temperature, type="Kc"))
                    checked += 1
        if "linear_ceiling" in candidate:
            ceiling = candidate["linear_ceiling"]
            def residual(temperature):
                values = {key: anchor["H298_J_mol"] - temperature * anchor["S298_J_mol_K"]
                          for key, anchor in candidate["anchors"].items()}
                shift = delta(result["reactions"]["propagation"], values)
                return math.log(reactions["propagation"].get_equilibrium_constant(temperature, type="Kc")
                                * result["monomer_concentration_mol_m3"]) - shift / (gas_constant * temperature)
            if ceiling["continuous_K"] is not None:
                close(ceiling["continuous_K"], brentq(residual, 600.0, 800.0), absolute=1e-7)
                close(residual(ceiling["continuous_K"]), 0.0, absolute=1e-8)
            residuals = [residual(temperature) for temperature in ceiling["temperatures_K"]]
            for actual, expected in zip(ceiling["rate_residuals"], residuals):
                close(actual, expected, absolute=1e-8)
            crossing = None
            for index in range(len(residuals) - 1):
                if residuals[index] * residuals[index + 1] < 0:
                    fraction = -residuals[index] / (residuals[index + 1] - residuals[index])
                    crossing = ceiling["temperatures_K"][index] + fraction * (
                        ceiling["temperatures_K"][index + 1] - ceiling["temperatures_K"][index])
                    break
            if crossing is None:
                assert ceiling["tabulated_K"] is None
            else:
                close(crossing, ceiling["tabulated_K"])
        for name, observation in candidate["boundaries"].items():
            if candidate["missing"] or not candidate["coolprop"] or name == "linear_at_800":
                continue
            if name in ("at_Tc", "above_Tc"):
                assert observation["error_type"] == "InputError"
            elif name == "below_triple":
                assert observation["error_type"] == "DatabaseError"
            else:
                assert "error" not in observation
    close(result["rates"]["ceiling_continuous_K"], brentq(lambda temperature: math.log(
        reactions["propagation"].get_equilibrium_constant(temperature, type="Kc")
        * result["monomer_concentration_mol_m3"]), 600.0, 800.0))
    return checked


def table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |",
                      *["| " + " | ".join(str(value) for value in row) + " |" for row in rows]])


def render_blocks(result):
    blocks = {}
    rows = result["inventory"]["rows"]
    categories = {
        "Complete Abraham + Mintz": [row for row in rows if not row["missing_g"] and not row["missing_h"]],
        "Abraham only": [row for row in rows if not row["missing_g"] and row["missing_h"]],
        "Mintz only": [row for row in rows if row["missing_g"] and not row["missing_h"]],
    }
    inventory_parts = [f"Pinned solvent entries: **{len(rows)}**. Solvent-library lines mentioning polymer/polystyrene/melt: **{len(result['inventory']['keyword_lines_in_solvent_library'])}**."]
    for heading, group in categories.items():
        inventory_parts.append(f"**{heading} ({len(group)}):** " + ", ".join(f"`{row['label']}`" for row in group) + ".")
    supported = [row for row in rows if not row["missing_g"] and not row["missing_h"] and row["coolprop"]]
    assert all(row["model_error_298"] is None for row in supported)
    inventory_parts.append("**Complete LSERs + mapped CoolProp identity** (the independent reference-solute call at the anchor succeeds for every row):\n\n" + table(
        ["Solvent", "CoolProp identity", "Tc (K)"], [[row["label"], row["coolprop"], f"{row['Tc_K']:.5f}"] for row in supported]))
    blocks["inventory"] = "\n\n".join(inventory_parts)
    blocks["provenance"] = table(["Quantity", "Reproduced value"], [
        ["Python", result["provenance"]["python"]], ["CoolProp", result["provenance"]["CoolProp"]],
        ["R (J/mol/K)", f"{result['provenance']['R_J_mol_K']:.12g}"],
        ["Pinned snapshot files", result["provenance"]["snapshot"]["files"]],
        ["Snapshot SHA256", f"`{result['provenance']['snapshot']['sha256']}`"],
        ["Prior probe script SHA256", f"`{result['provenance']['prior_script_sha256']}`"],
        ["Monomer concentration (mol/m³)", f"{result['monomer_concentration_mol_m3']:.1f}"],
        ["Default ceiling grid (K)", f"{result['rates']['propagation_grid']['T'][0]:.0f}–{result['rates']['propagation_grid']['T'][-1]:.0f}, step {result['rates']['propagation_grid']['T'][1] - result['rates']['propagation_grid']['T'][0]:.0f}"],
    ])
    blocks["reactions"] = "\n\n".join(f"**{REACTION_LABELS[name]}:** `{' + '.join(pair['reactants'])} -> {' + '.join(pair['products'])}` (Δn = {pair['delta_n']:+d})."
                                                    for name, pair in result["reactions"].items())
    blocks["baseline"] = table(
        ["T (K)", "Kc,association (m³/mol)", "k_homolysis (s⁻¹)", "Kc,H prior", "Kc,H chain", "Kc,propagation (m³/mol)"],
        [[f"{row['T_K']:.0f}", f"{result['rates']['Kc_association'][index]:.6e}",
          f"{result['rates']['homolysis_s-1'][index]:.6e}", f"{row['reactions']['H_baseline']['Kc']:.6e}",
          f"{row['reactions']['H_chain']['Kc']:.6e}", f"{row['reactions']['propagation']['Kc']:.6e}"] for index, row in enumerate(result["gas"])])
    blocks["gas_magnitudes"] = table(["Reaction", "T (K)", "ΔrGgas (kJ/mol)", "Kc,gas (SI concentration convention)"],
        [[REACTION_LABELS[name], f"{row['T_K']:.0f}", f"{values['dG_J_mol']/1000:.6f}", f"{values['Kc']:.6e}"]
         for row in result["gas"] for name, values in row["reactions"].items()])
    solute_rows = []
    for key, entry in result["species"].items():
        descriptors = entry["descriptors"]
        origin = "library" if "Solute library:" in entry["comment"] else "group additivity"
        if "radical(" in entry["comment"]:
            origin += " + radical HBI"
        solute_rows.append([f"`{key}`", *[f"{descriptors[field]:.6f}" for field in ("S", "B", "E", "L", "A", "V")], origin])
    blocks["solutes"] = f"Selected unique species: **{len(result['species'])}**; failures: **{len(result['solute_failures'])}**. Default descriptor lookup and forced group additivity both succeed on every selected species.\n\n" + table(
        ["Species (SMILES)", "S", "B", "E", "L", "A", "V", "Default origin"], solute_rows)
    blocks["radicals"] = f"Selected carbon-radical species lacking a tabulated radical correction: **{len(result['radical_audit'])}**.\n\n" + table(
        ["Radical (SMILES)", "Saturated analogue", "Matched node", "Explicit correction lookup", "max abs(ΔS,ΔB,ΔE,ΔL,ΔA)", "ΔV"],
        [[f"`{row['smiles']}`", f"`{row['saturated_smiles']}`", row["matched_node"], "KeyError: no parent with data",
          f"{max(abs(row['radical_minus_saturated'][field]) for field in ('S','B','E','L','A')):.3e}", f"{row['radical_minus_saturated']['V']:.6f}"]
         for row in result["radical_audit"]])
    blocks["magnitude_comparison"] = table(["Toluene / linear H−TS, reaction", "T (K)", "abs(ΔΔGsolv) / abs(ΔrGgas)", "Kc,melt / Kc,gas"],
        [[REACTION_LABELS[name], f"{row['T_K']:.0f}",
          f"{abs(values['ddG_J_mol'] / result['gas'][index]['reactions'][name]['dG_J_mol']):.6f}", f"{values['ratio']:.6e}"]
         for index, row in enumerate(result["candidates"]["toluene"]["models"]["linear_HS"])
         for name, values in row["reactions"].items()])
    candidate_rows = []
    for label, candidate in result["candidates"].items():
        critical = f"{candidate['Tc_K']:.5f}" if candidate["Tc_K"] else f"{result['hexadecane_literature_Tc']['value_K']:.0f} ± {result['hexadecane_literature_Tc']['uncertainty_K']:.0f} [NIST]"
        statuses = [row.get("error_type", "computed") for row in candidate["models"]["Kfactor"]]
        candidate_rows.append([label, "G only" if candidate["missing"] else "G + H", critical,
                               candidate["coolprop"] or "none", "/".join(statuses), "unavailable" if candidate["missing"] else "explicit extrapolation at all requested T"])
    blocks["candidates"] = table(["Stand-in", "LSERs", "Tc (K)", "CoolProp mapping", "K-factor at 600/700/800 K", "Linear H−TS at 600–800 K"], candidate_rows)
    anchors = []
    for label, candidate in result["candidates"].items():
        if not candidate["missing"]:
            for name, pair in result["reactions"].items():
                anchors.append([label, REACTION_LABELS[name],
                                f"{delta(pair, {key: value['G298_J_mol'] for key, value in candidate['anchors'].items()})/1000:.6f}",
                                f"{delta(pair, {key: value['H298_J_mol'] for key, value in candidate['anchors'].items()})/1000:.6f}",
                                f"{delta(pair, {key: value['S298_J_mol_K'] for key, value in candidate['anchors'].items()}):.6f}"])
    blocks["anchors"] = table(["Stand-in", "Reaction", "ΔΔG298 (kJ/mol)", "ΔΔH298 (kJ/mol)", "ΔΔS298 (J/mol/K)"], anchors)
    shifts = []
    kfactor = []
    for label, candidate in result["candidates"].items():
        for model, rows in candidate["models"].items():
            for row in rows:
                if "error" in row:
                    continue
                target = shifts if model == "linear_HS" else kfactor
                for name, values in row["reactions"].items():
                    gas_row = result["gas"][result["temperatures_K"].index(row["T_K"])]
                    fraction = abs(values["ddG_J_mol"] / gas_row["reactions"][name]["dG_J_mol"])
                    target.append([label, REACTION_LABELS[name], f"{row['T_K']:.0f}",
                                   f"{values['ddG_J_mol']/1000:+.6f}", f"{values['ratio']:.6e}",
                                   f"{values['melt_Kc']:.6e}", f"{fraction:.6f}"])
    headers = ["Stand-in", "Reaction", "T (K)", "ΔΔGsolv (kJ/mol)", "Kc,melt/Kc,gas", "Kc,melt (SI)", "abs(ΔΔGsolv)/abs(ΔrGgas)"]
    blocks["linear_shifts"] = table(headers, shifts)
    blocks["kfactor_shifts"] = table(headers, kfactor)
    ceiling_rows = [["gas", f"{result['rates']['ceiling_tabulated_K']:.6f}", f"{result['rates']['ceiling_continuous_K']:.6f}", "0.000000"]]
    for label, candidate in result["candidates"].items():
        if "linear_ceiling" in candidate:
            ceiling = candidate["linear_ceiling"]
            ceiling_rows.append([label + " / linear H−TS", f"{ceiling['tabulated_K']:.6f}" if ceiling['tabulated_K'] else "no crossing",
                                 f"{ceiling['continuous_K']:.6f}" if ceiling['continuous_K'] else "no crossing",
                                 f"{ceiling['continuous_K'] - result['rates']['ceiling_continuous_K']:+.6f}" if ceiling['continuous_K'] else "N/A"])
    blocks["ceilings"] = table(["Reference / model", "Compiler-grid crossing (K)", "Continuous Kc[M]=1 root (K)", "Continuous shift from gas (K)"], ceiling_rows)
    boundaries = []
    for label, candidate in result["candidates"].items():
        for name, observation in candidate["boundaries"].items():
            status = observation.get("error_type", f"Gsolv={observation.get('gibbs_J_mol', 0.0)/1000:+.6f} kJ/mol")
            boundaries.append([label, name, f"{observation['T_K']:.5f}", status])
    blocks["boundaries"] = f"Fixed reference solute: `{result['boundary_reference_smiles']}`.\n\n" + table(["Stand-in", "Call", "T (K)", "Observed result (one fixed reference solute)"], boundaries)
    g_only = result["candidates"]["ethylbenzene"]["G298_reactions_J_mol"]
    blocks["ethylbenzene_anchor"] = table(["Reaction", "Ethylbenzene ΔΔG298 (kJ/mol)"],
                                         [[REACTION_LABELS[name], f"{value/1000:.6f}"] for name, value in g_only.items()])
    return blocks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("--snapshot", type=Path, default=Path("/tmp/i034-reproduce/database"))
    parser.add_argument("--report", type=Path, default=REPORT)
    parser.add_argument("--update-report-blocks", action="store_true")
    args = parser.parse_args()
    result = json.loads(args.results.read_text())
    count = verify(result, args.snapshot)
    blocks = render_blocks(result)
    report = args.report.read_text()
    for name, text in blocks.items():
        begin = f"<!-- BEGIN I034:{name} -->"
        end = f"<!-- END I034:{name} -->"
        assert report.count(begin) == report.count(end) == 1, f"missing/duplicated block {name}"
        prefix, remainder = report.split(begin)
        actual, suffix = remainder.split(end)
        if args.update_report_blocks:
            report = prefix + begin + "\n" + text + "\n" + end + suffix
        else:
            assert actual == "\n" + text + "\n", f"report numbers differ in {name}"
    if args.update_report_blocks:
        args.report.write_text(report)
    print("I-034 gas baseline, descriptors, independent LSER/K-factor calculations and ceilings verified")
    print(f"I-034 {count} corrected reaction-temperature rows and {len(blocks)} report blocks reproduced")


if __name__ == "__main__":
    main()
