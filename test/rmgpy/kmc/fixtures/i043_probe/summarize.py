"""Audit sampled structures and render measured I043 tables.

No candidate list is promoted to a thermochemical ensemble without refinement,
minima validation, and rotor/conformer accounting. Commands are in the report.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from rdkit import Chem

from run_ensembles import SPECIES, PLAN, read_xyz, allowed_scratch, limit_cores

REPORT = Path(__file__).resolve().parent.parent / "I043_diphenyl_nonadditivity.md"
START = "<!-- BEGIN I043:measured -->"
END = "<!-- END I043:measured -->"


def stereo_check(smiles, frame):
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if molecule.GetNumAtoms() != len(frame["atoms"]):
        raise AssertionError("atom count changed during search")
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for index, line in enumerate(frame["atoms"]):
        symbol, x, y, z = line.split()
        if symbol != molecule.GetAtomWithIdx(index).GetSymbol():
            raise AssertionError("atom ordering changed")
        conformer.SetAtomPosition(index, [float(x), float(y), float(z)])
    conformer.Set3D(True)
    molecule.AddConformer(conformer)
    Chem.AssignStereochemistryFrom3D(molecule, replaceExistingTags=True)
    actual = Chem.MolToSmiles(Chem.RemoveHs(molecule))
    expected = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    if actual != expected:
        raise AssertionError(f"stereo changed: {actual} != {expected}")
    return actual


def audit(scratch):
    manifest = json.loads((scratch / "species.json").read_text())
    rows = []
    for name, smiles in SPECIES.items():
        directory = scratch / "ensembles" / name
        marker = directory / "completed.json"
        if not marker.exists():
            rows.append({"name": name, "atoms": manifest[name]["natoms"], "complete": False})
            continue
        record = json.loads(marker.read_text())
        frames = read_xyz(directory / "crest_conformers.xyz")
        digest = hashlib.sha256((directory / "crest_conformers.xyz").read_bytes()).hexdigest()
        if digest != record["ensemble_sha256"] or len(frames) != record["frames"]:
            raise AssertionError(f"ensemble bytes changed: {name}")
        energy = [float(frame["comment"].split()[0]) for frame in frames]
        if energy != record["energy_hartree"]:
            raise AssertionError("stored energies differ from XYZ")
        selected = [i for i, value in enumerate(energy)
                    if (value - min(energy)) * 2625.499639 < PLAN["refinement_window_kJ_mol"]]
        if selected != record["selected_indices_within_12_kJ"]:
            raise AssertionError("energy selection changed")
        for frame in frames:
            stereo_check(smiles, frame)
        rows.append({"name": name, "atoms": manifest[name]["natoms"], "complete": True,
                     "candidates": len(frames), "within_12_kJ": len(selected),
                     "minimum_energy_Eh": min(energy), "search_s": record["elapsed_s"],
                     "stereo_checked": len(frames), "ensemble_sha256": digest})
        selection = scratch / "composite" / name / "selection.json"
        if selection.exists():
            selected = json.loads(selection.read_text())
            indices = selected["unique_minima_indices"]
            rows[-1]["unique_minima"] = len(indices)
            rows[-1]["chemical_minima"] = selected.get("chemical_minima_count",len(indices))
            rows[-1]["permutation_partners"] = len(selected.get("permutation_partners",[]))
            rows[-1]["reflection_partners"] = len(selected.get("reflection_partners",[]))
            rows[-1]["minimum_frequency_cm1"] = min(
                min(v for v in json.loads((scratch / "xtb_checks" / name / f"{i:04d}" / "result.json").read_text())["frequencies_cm1"]
                    if abs(v) > 1e-6)
                for i in indices)
    pilots = {}
    for path in sorted((scratch / "pilot").glob("*.json")):
        data = json.loads(path.read_text())
        if not data.get("frequencies_complete"):
            continue
        n = len(data["atoms"])
        molecule_name = path.stem.rsplit("_",1)[0]
        if molecule_name in SPECIES:
            frame = {"atoms": [symbol+" "+" ".join(str(v) for v in point)
                                for symbol,point in zip(data["atoms"],data["coordinates_angstrom"])]}
            stereo_check(SPECIES[molecule_name],frame)
        if len(data["frequencies_cm1_real"]) != 3*n-6:
            raise AssertionError("wrong number of vibrational modes")
        gradient = np.array(data["gradient_Eh_bohr"])
        if np.max(np.linalg.norm(gradient, axis=1)) > 3e-4:
            raise AssertionError("pilot is not a converged geometry")
        imaginary = float(max(data["frequencies_cm1_imag"]))
        pilots[path.stem] = {
            "level": data["level"], "atoms": n, "energy_Eh": data["energy_hartree"],
            "minimum_frequency_cm1": min(data["frequencies_cm1_real"]),
            "maximum_imaginary_cm1": imaginary,
            "optimization_s": data["optimization_s"], "hessian_s": data["hessian_s"],
            "minimum_confirmed": imaginary <= 20,
        }
    baseline = json.loads((scratch / "baseline/baseline.json").read_text())
    if baseline["snapshot"]["sha256"] != "61e2c7357bb4945dce6e6a1a39eccf216e945a01771e461391e7975165ea29eb":
        raise AssertionError("pinned database snapshot changed")
    result = {"candidates": rows, "pilots": pilots, "baseline": baseline}
    reflection = scratch / "reflection_validation/diphenylpropane_pbe.json"
    if reflection.exists():
        fresh = json.loads(reflection.read_text())
        parent = json.loads((scratch / "composite/diphenylpropane/0000/pbe/result.json").read_text())
        difference = (fresh["energy_hartree"]-parent["energy_hartree"])*2625.499639
        if abs(difference) > 0.02:
            raise AssertionError("fresh reflection energy check failed")
        result["reflection_validation"] = {"level": "PBE", "energy_difference_kJ_mol": difference}
    thread_probe = scratch / "thread_probe/iso_0000_pbe_blas1.json"
    if thread_probe.exists():
        fresh = json.loads(thread_probe.read_text())
        parent = json.loads((scratch / "composite/triphenylheptane_iso/0000/pbe/result.json").read_text())
        difference = (fresh["energy_hartree"]-parent["energy_hartree"])*2625.499639
        if abs(difference) > 0.02 or fresh["input_sha256"] != parent["input_sha256"]:
            raise AssertionError("fresh trimer single-point check failed")
        result["trimer_energy_validation"] = {"energy_difference_kJ_mol": difference}
    thermochemistry = scratch / "thermochemistry.json"
    if thermochemistry.exists():
        result["thermochemistry"] = json.loads(thermochemistry.read_text())
    result["boundary_probes"]={}
    for path in sorted((scratch/"hard_core_probe").glob("*.json")):
        probe=json.loads(path.read_text())
        source=scratch/"xtb_checks"/probe["name"]/f"{probe['index']:04d}"/"xtbopt.xyz"
        if probe["source_input_sha256"]!=hashlib.sha256(source.read_bytes()).hexdigest():
            raise AssertionError("steric boundary probe geometry changed")
        result["boundary_probes"][path.stem]=probe
    return result


def render(data):
    lines = ["| Molecule / stereo sequence | Atoms | CREST candidates | Within 12 kJ/mol | Stereo checked | Chemical minima incl. mirrors | Rotor-space wells | Mirror partners | Proper-permutation images | Lowest nonzero frequency (cm⁻¹) |",
             "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in data["candidates"]:
        columns = [str(row[key]) for key in ("candidates", "within_12_kJ", "stereo_checked")] if row["complete"] else ["pending"]*3
        minima = str(row.get("unique_minima", "pending"))
        chemical = str(row.get("chemical_minima","pending"))
        permutations = str(row.get("permutation_partners","pending"))
        frequency = f"{row['minimum_frequency_cm1']:.3f}" if "minimum_frequency_cm1" in row else "pending"
        mirrors = str(row.get("reflection_partners", "pending"))
        lines.append(f"| {row['name']} | {row['atoms']} | " + " | ".join(columns) + f" | {chemical} | {minima} | {mirrors} | {permutations} | {frequency} |")
    lines += ["", "| DFT pilot | Level | Energy (Eh) | Lowest real frequency (cm⁻¹) | Largest imaginary frequency (cm⁻¹) |",
              "| --- | --- | ---: | ---: | ---: |"]
    for name, row in data["pilots"].items():
        lines.append(f"| {name} | {row['level']} | {row['energy_Eh']:.9f} | {row['minimum_frequency_cm1']:.3f} | {row['maximum_imaginary_cm1']:.3f} |")
    if "reflection_validation" in data:
        lines += ["", f"Fresh PBE single-point reflection check, diphenylpropane: energy difference {data['reflection_validation']['energy_difference_kJ_mol']:.6f} kJ/mol (acceptance tolerance 0.02 kJ/mol)."]
    if "trimer_energy_validation" in data:
        lines += ["", f"Fresh PBE trimer single-point check with BLAS threading changed: energy difference {data['trimer_energy_validation']['energy_difference_kJ_mol']:.9f} kJ/mol (acceptance tolerance 0.02 kJ/mol)."]
    baseline = data["baseline"]
    lines += ["", "| Reproduced RMG quantity | Value |", "| --- | ---: |",
              f"| Gas continuous Tc at 1 mol/L (K) | {baseline['continuous_Tc_K']:.6f} |",
              f"| Compiled grid Tc (K) | {baseline['compiled_grid_Tc_K']:.6f} |",
              f"| Propagation ΔH° at 298.15 K (kJ/mol) | {baseline['propagation'][0]['H_J_mol']/1000:.6f} |",
              f"| Propagation ΔS° at 298.15 K (J/mol/K) | {baseline['propagation'][0]['S_J_mol_K']:.6f} |"]
    for name, cycle in baseline["cycles"].items():
        row = cycle["thermo"][0]
        lines += [f"| {name}, ΔH° at 298.15 K (kJ/mol) | {row['H_J_mol']/1000:.6f} |",
                  f"| {name}, ΔS° at 298.15 K (J/mol/K) | {row['S_J_mol_K']:.6f} |"]
    if data["boundary_probes"]:
        lines += ["", "Steric boundary controls evaluate the uncut GFN2 Hamiltonian at retained points nearest the 0.45 Å boundary. These individual points do not bound the whole excluded region.", "",
                  "| Molecule / well | Side | Closest atoms | Distance (Å) | V above minimum (kJ/mol) | ln Boltzmann factor at 800 K | Status |",
                  "| --- | --- | --- | ---: | ---: | ---: | --- |"]
        for probe in data["boundary_probes"].values():
            for point in probe["points"]:
                energy=f"{point['energy_above_minimum_J_mol']/1000:.3f}" if point['status']=='finite' else '—'
                log_factor=f"{point['log_Boltzmann_factor_at_800_K']:.3f}" if point['status']=='finite' else '—'
                lines.append(f"| {probe['name']} / {probe['index']} | {point['side']} | {'–'.join(point['pair_symbols'])} | {point['minimum_distance_A']:.6f} | {energy} | {log_factor} | {point['status']} |")
    if "thermochemistry" not in data:
        return "\n".join(lines)
    thermo = data["thermochemistry"]
    lines += ["", "| Molecule | PBE electronic order, low to high | BLYP electronic order, low to high | PBE minimum E (Eh) | BLYP minimum E (Eh) |",
              "| --- | --- | --- | ---: | ---: |"]
    for name, species in thermo["species"].items():
        orders = species["electronic_order_indices"]
        energies = species["minimum_electronic_energy_hartree"]
        lines.append(f"| {name} | {', '.join(map(str,orders['pbe']))} | {', '.join(map(str,orders['blyp']))} | {energies['pbe']:.9f} | {energies['blyp']:.9f} |")
    lines += ["", "### Conditional effect on compiled propagation (main table)", "",
              "H and S anchors below are at 298.15 K. Tc uses the temperature-dependent correction, one-bar gas thermo and 1 mol/L monomer. Each row is a separate transfer hypothesis; rows are not added together.", "",
              "| Finding | Electronic level | δΔH° (kJ/mol) | δΔS° (J/mol/K) | Corrected ΔH° (kJ/mol) | Corrected ΔS° (J/mol/K) | Gas Tc (K) | Rotor MC SE of Tc (K) |",
              "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for name, levels in thermo["conditional_corrections"].items():
        for level, row in levels.items():
            roots = row["roots_298_to_800_K"]
            tc = ", ".join(f"{value:.3f}" for value in roots) or "no root in 298–800"
            se = ", ".join(f"{value:.3f}" for value in row["root_MC_standard_errors_K"]) or "—"
            lines.append(f"| {name} | {level.upper()} | {row['delta_H298_kJ_mol']:.3f} | {row['delta_S298_J_mol_K']:.3f} | {row['propagation_H298_kJ_mol']:.3f} | {row['propagation_S298_J_mol_K']:.3f} | {tc} | {se} |")
    lines += ["", "| Finding | Two-method δΔH° range (kJ/mol) | Two-method δΔS° range (J/mol/K) | Two-method Tc range (K) |",
              "| --- | ---: | ---: | ---: |"]
    for name, spread in thermo["method_spread"].items():
        def interval(key):
            return f"{spread[key]['minimum']:.3f} to {spread[key]['maximum']:.3f}" if key in spread else "not bracketed"
        lines.append(f"| {name} | {interval('delta_H298_kJ_mol')} | {interval('delta_S298_J_mol_K')} | {interval('Tc_K')} |")
    lines += ["", "### Balanced cycle thermochemistry", "",
              "| Cycle | Level | T (K) | ΔH° (kJ/mol) | ΔS° (J/mol/K) | Rotor MC SE H (kJ/mol) | Rotor MC SE S (J/mol/K) |",
              "| --- | --- | ---: | ---: | ---: | ---: | ---: |"]
    for name, cycle in thermo["cycles"].items():
        for level, rows in cycle["levels"].items():
            for row in rows:
                error = row["MC_standard_errors_H_S_G"]
                lines.append(f"| {name} | {level.upper()} | {row['T_K']:.2f} | {row['H_J_mol']/1000:.3f} | {row['S_J_mol_K']:.3f} | {error[0]/1000:.3f} | {error[1]:.3f} |")
    lines += ["", "Cycle entropy components below are in J/mol/K. Basin mixing includes coordinate images whose multiplicity is canceled by external symmetry. Translation/rotation are capped-molecule contributions and are not automatically transferable local interactions.", "",
              "| Cycle | Level | T (K) | Translation | External rotation | Non-torsional vibrations | Coupled rotors | Basin mixing | Total raw ΔS° |",
              "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for name,diagnostic in thermo["cycles"].items():
        for level,rows in diagnostic["levels"].items():
            for row in (rows[0],rows[-1]):
                parts=row["entropy_components"]
                values=[parts[key] for key in ("S_translation_J_mol_K","S_external_rotation_J_mol_K",
                    "S_non_torsional_vibration_J_mol_K","S_coupled_rotor_J_mol_K","S_basin_mixing_J_mol_K")]
                lines.append(f"| {name} | {level.upper()} | {row['T_K']:.2f} | "+" | ".join(f"{v:.3f}" for v in values)+f" | {row['S_J_mol_K']:.3f} |")
    if "quadrature_comparison" in thermo:
        comparison=thermo["quadrature_comparison"]
        lines += ["", f"Diphenylpropane quadrature comparison: {comparison['points'][0]} versus {comparison['points'][1]} points per basin on the same final partition. Other species are unchanged; these are first-cycle changes, not new chemistry.", "",
                  "| Level | T (K) | ΔH° change (kJ/mol) | ΔS° change (J/mol/K) | Higher-count MC SE H (kJ/mol) | Higher-count MC SE S (J/mol/K) |",
                  "| --- | ---: | ---: | ---: | ---: | ---: |"]
        for level,rows in comparison["levels"].items():
            for low,high in zip(rows,thermo["cycles"]["diphenylpropane_isodesmic"]["levels"][level]):
                error=high["MC_standard_errors_H_S_G"]
                lines.append(f"| {level.upper()} | {high['T_K']:.2f} | {(high['H_J_mol']-low['H_J_mol'])/1000:.3f} | {high['S_J_mol_K']-low['S_J_mol_K']:.3f} | {error[0]/1000:.3f} | {error[1]:.3f} |")
    lines += ["", "| Conditional capped-fragment channel | Level | ΔH°298 (kJ/mol) | ΔS°298 (J/mol/K) | Atactic weight |",
              "| --- | --- | ---: | ---: | ---: |"]
    for name, channel in thermo["stereo_channel_residuals"].items():
        for level,rows in channel["levels"].items():
            lines.append(f"| {name} | {level.upper()} | {rows[0]['H_J_mol']/1000:.3f} | {rows[0]['S_J_mol_K']:.3f} | {channel['unconditional_atactic_weight']:.2f} |")
    lines += ["", "The inherited corrected-G3MP2 comparison is +7.15 ± 5.58 kJ/mol and 30.4 J/mol/K at 298.15 K.", "",
              "| DPP cycle level | H minus inherited G3MP2 (kJ/mol) | S minus inherited G3MP2 (J/mol/K) | S minus pinned RMG (J/mol/K) |",
              "| --- | ---: | ---: | ---: |"]
    for level, rows in thermo["cycles"]["diphenylpropane_isodesmic"]["levels"].items():
        row = rows[0]
        rmg_s = baseline["cycles"]["diphenylpropane_isodesmic"]["thermo"][0]["S_J_mol_K"]
        lines.append(f"| {level.upper()} | {row['H_J_mol']/1000-7.15:.3f} | {row['S_J_mol_K']-30.4:.3f} | {row['S_J_mol_K']-rmg_s:.3f} |")
    lines += ["", "### Molecular ensemble H(T), S(T)", "",
              f"H* is H minus the tabulated minimum electronic energy of that molecule at that electronic level; it includes zero-point and thermal energy and is not a formation enthalpy. Full H (J/mol) = minimum E (Eh) × {thermo['hartree_J_mol']:.9f} + 1000 H*. Reaction cycles above use full electronic plus thermal H, not differences between H* columns.", "",
              "Srep uses a specified enantiomer when the molecule is chiral. Srac includes its equal mirror population at the same total gas pressure; it equals Srep for achiral molecules. The fixed-sequence cycles use Srep; including Srac consistently leaves the weighted added-unit residual unchanged.", "",
              "| Molecule | T (K) | PBE H* (kJ/mol) | PBE Srep (J/mol/K) | PBE Srac (J/mol/K) | BLYP H* (kJ/mol) | BLYP Srep (J/mol/K) | BLYP Srac (J/mol/K) |",
              "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for name, species in thermo["species"].items():
        for pbe, blyp in zip(species["levels"]["pbe"], species["levels"]["blyp"]):
            addition = species["racemic_gas_entropy_addition_J_mol_K"]
            lines.append(f"| {name} | {pbe['T_K']:.2f} | {pbe['H_above_electronic_minimum_kJ_mol']:.3f} | {pbe['S_J_mol_K']:.3f} | {pbe['S_J_mol_K']+addition:.3f} | {blyp['H_above_electronic_minimum_kJ_mol']:.3f} | {blyp['S_J_mol_K']:.3f} | {blyp['S_J_mol_K']+addition:.3f} |")
    lines += ["", "| Molecule | Points per basin | Importance proposal | External symmetry divisor | Internal rotor divisor | PBE conformer S298 (J/mol/K) | BLYP conformer S298 (J/mol/K) | Lowest basin ESS at 298 K, PBE/BLYP |",
              "| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: |"]
    for name, species in thermo["species"].items():
        a, b = (species["levels"][level][0] for level in ("pbe", "blyp"))
        sym = species["symmetry"]
        overrides=thermo['samples_per_basin_overrides'][name]
        counts=[int(overrides.get(str(i),thermo['samples_per_species'][name])) for i in species['indices']]
        count_text=str(min(counts)) if min(counts)==max(counts) else f"{min(counts)}–{max(counts)}"
        lines.append(f"| {name} | {count_text} | {thermo['importance_proposals'][name]} | {sym['external_divisor']:.0f} | {sym['internal_rotor_divisor']} | {a['S_conformer_J_mol_K']:.3f} | {b['S_conformer_J_mol_K']:.3f} | {a['minimum_effective_samples']:.1f}/{b['minimum_effective_samples']:.1f} |")
    allocated=[(n,i,count) for n,values in thermo['samples_per_basin_overrides'].items() for i,count in values.items()]
    lines += ["", "| Steric zero-weight points among selected proposal draws | Count | All draws |", "| --- | ---: | ---: |"]
    lines += [f"| {n} | {s['hard_core_zero_weight_points']} | {s['total_quadrature_points']} |" for n,s in thermo['species'].items()]
    if allocated:
        lines += ["", "| Extra point allocation | Well index | Points |", "| --- | ---: | ---: |"]
        lines += [f"| {n} | {i} | {count} |" for n,i,count in allocated]
    if 'precision_comparison' in thermo:
        comparison=thermo['precision_comparison']
        lines += ["", f"Same-proposal convergence: {comparison['base_points']} points in every well versus the selected allocations. The physical partitions and random-stream prefixes are identical.", "",
                  "| Cycle | Level | Maximum absolute change in H over 298–800 K (kJ/mol) | In S (J/mol/K) | In G (kJ/mol) |",
                  "| --- | --- | ---: | ---: | ---: |"]
        for name,levels in comparison['cycles'].items():
            for level,base_rows in levels.items():
                selected=thermo['cycles'][name]['levels'][level]
                changes=[(abs(a['H_J_mol']-b['H_J_mol'])/1000,abs(a['S_J_mol_K']-b['S_J_mol_K']),
                          abs(a['H_J_mol']-b['H_J_mol']-a['T_K']*(a['S_J_mol_K']-b['S_J_mol_K']))/1000)
                         for a,b in zip(selected,base_rows)]
                maximum=np.max(changes,axis=0)
                lines.append(f"| {name} | {level} | {maximum[0]:.3f} | {maximum[1]:.3f} | {maximum[2]:.3f} |")
    stereo = thermo["stereo_accounting"]
    lines += ["", f"R ln 2 = {stereo['R_ln_2_J_mol_K']:.6f} J/mol/K remains in the original compiler once.",
              f"The global mirror contribution excluded from both fixed-sequence averages would change the added-unit entropy by {stereo['global_mirror_increment_if_included_J_mol_K']:.6f} J/mol/K if included consistently.",
              f"The QM capped molecules' external end-symmetry increment is {stereo['external_end_symmetry_increment_J_mol_K']:.6f} J/mol/K; the GAV end increment is {stereo['GAV_external_end_symmetry_increment_J_mol_K']:.6f}, and its local optical increment is {stereo['GAV_local_optical_increment_J_mol_K']:.6f}. Both sides are normalized before transfer. Their difference adds {stereo['matched_entropy_offset_J_mol_K']:.6f} to the raw QM-minus-GAV entropy residual.",
              f"The excluded class-mixing increment is {stereo['excluded_class_mixing_entropy_increment_J_mol_K']:.6f} J/mol/K. QM end symmetry + class mixing + global mirror contributions give {stereo['end_plus_class_plus_mirror_increment_J_mol_K']:.6f} J/mol/K. GAV end + local optical increments give the same one new-centre bit. No independent R ln 2 is added to the compiler or the fixed-sequence averages."]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--write-report", action="store_true")
    args = parser.parse_args()
    limit_cores()
    scratch = allowed_scratch(args.scratch)
    data = audit(scratch)
    rendered = render(data)
    (scratch / "audit.json").write_text(json.dumps(data, indent=2) + "\n")
    current = REPORT.read_text()
    prefix, body = current.split(START, 1)
    old, suffix = body.split(END, 1)
    if args.write_report:
        REPORT.write_text(prefix + START + "\n" + rendered + "\n" + END + suffix)
    elif old.strip() != rendered:
        raise AssertionError("report differs from reproduced tables")
    print("I043 ensemble hashes, energy selections, all sampled stereo assignments and pilot stationarity verified")
    print("I043 pinned RMG baseline and measured report tables verified")


if __name__ == "__main__":
    main()
