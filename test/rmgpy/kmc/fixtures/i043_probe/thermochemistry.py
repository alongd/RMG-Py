"""Evaluate the composite conformer/rotor partition and diagnostic cycles.

Run with rmg_env. Rotors are integrated on a fixed-geometry torsional manifold;
non-torsional vibrations use projected Hessians and a Pitzer-Gwinn correction.
The report must identify these approximations and the conditional chain transfer.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.linalg import eigvalsh
from scipy.special import logsumexp
from rdkit import Chem

from run_ensembles import SPECIES, SCRATCH, allowed_scratch
from rotors import R, KB, HPLANCK, AMU_KG, EH_J, NA, quadrature_filename

EH_MOL = EH_J*NA
THETA_CM1 = HPLANCK*29979245800./KB
TEMPERATURES = (298.15, 300., 400., 500., 600., 700., 800.)
CYCLE_COEFFICIENTS = {
    "diphenylpropane_isodesmic": {"diphenylpropane": -1, "ethane": -1,
        "ethylbenzene": 1, "n_propylbenzene": 1},
    "atactic_added_unit": {"diphenylpentane_meso": -0.5, "diphenylpentane_racemo": -0.5,
        "triphenylheptane_iso": 0.25, "triphenylheptane_syndio": 0.25,
        "triphenylheptane_hetero": 0.5, "cumene": -1, "n_propylbenzene": -1,
        "ethylbenzene": 1, "ethane": 1},
}


def external_symmetry(smiles, baseline, definitions, coordinates):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    counts = []
    for chirality in (False, True):
        matches = mol.GetSubstructMatches(mol, uniquify=False, useChirality=chirality, maxMatches=1000000)
        if len(matches) == 1000000:
            raise AssertionError("automorphism enumeration was truncated")
        counts.append(len(matches))
    # Constitutional automorphisms include reflections and tetrahedral ligand
    # exchanges that no proper rotation/internal torsion can realize. Filter
    # using signed tetrahedron volumes, including formally achiral CH2/CH3
    # atoms. A stereo/constitutional ratio alone incorrectly halves the
    # trimers' already orientation-restricted graph factor.
    xyz = np.array(coordinates)
    tetrahedra = [sorted(a.GetIdx() for a in atom.GetNeighbors())
                  for atom in mol.GetAtoms() if atom.GetDegree() == 4]
    def volume(indices):
        points = xyz[indices]
        return float(np.linalg.det(points[:3]-points[3]))
    volumes = [volume(indices) for indices in tetrahedra]
    if any(abs(value) < 1e-4 for value in volumes):
        raise AssertionError("degenerate tetrahedral geometry in symmetry audit")
    proper = sum(all(v*volume([match[i] for i in indices]) > 0
                     for v,indices in zip(volumes,tetrahedra))
                 for match in mol.GetSubstructMatches(mol,uniquify=False,useChirality=True,maxMatches=1000000))
    non_optical = baseline["sigma"]*2**baseline["optical_atom_half_factors"]
    internal = math.prod(r["symmetry"] for r in definitions)
    external = proper/internal
    if external < 1-1e-10 or abs(external-round(external)) > 1e-10:
        raise AssertionError("nonintegral fixed-stereo external symmetry")
    return {"external_divisor": external, "internal_rotor_divisor": internal,
            "nonoptical_RMG_sigma": non_optical,
            "constitutional_automorphisms": counts[0], "stereo_automorphisms": counts[1],
            "orientation_preserving_automorphisms": proper}


def conformer_thermal(record, temperature, symmetry):
    modes = record["modes"]
    d = len(record["definitions"])
    inertia = np.array(modes["inertia_amu_A2"])*AMU_KG*1e-20
    v = np.array([value if value is not None else 0.0 for value in record["energies_above_minimum_J_mol"]])
    mask = np.array([value is not None for value in record["energies_above_minimum_J_mol"]])
    density = np.array(record["proposal_densities_rad_minus_d"])
    log_weight = -v/(R*temperature)-np.log(density)
    log_weight[~mask] = -np.inf
    largest = np.max(log_weight)
    weight = np.exp(log_weight-largest)
    if np.sum(weight) == 0:
        raise AssertionError("rotor integral is empty")
    normalized = weight/np.sum(weight)
    integral_log = float(largest + math.log(np.mean(weight)))
    mean_v = float(normalized@v)
    mean_weight = np.mean(weight)
    relative_weight = weight/mean_weight
    log_influence = relative_weight-1
    h_influence = (v-mean_v)*relative_weight
    s_influence = R*log_influence + h_influence/temperature
    effective_samples = float(np.sum(weight)**2/np.sum(weight**2))
    log_qtor = d/2*math.log(2*math.pi*KB*temperature)-d*math.log(HPLANCK)
    log_qtor += 0.5*np.linalg.slogdet(inertia)[1] + integral_log
    curvature = np.array(modes["curvature_J_mol_rad2"])/NA
    eigenvalues = eigvalsh(curvature,inertia)
    if np.min(eigenvalues) <= 0:
        raise AssertionError("rigid torsional curvature is unstable")
    # The Pitzer-Gwinn reference must match the fixed-angle classical surface,
    # rather than the Schur complement that permits non-torsional relaxation.
    xt = HPLANCK*np.sqrt(eigenvalues)/(2*math.pi*KB*temperature)
    # log[x/(2 sinh(x/2))], evaluated without large-temperature cancellation.
    correction = float(np.sum(np.log(xt)-xt/2-np.log(-np.expm1(-xt))))
    correction_h = R*temperature*float(np.sum(-1+xt/2/np.tanh(xt/2)))
    log_qtor += correction
    htor = d/2*R*temperature + mean_v + correction_h

    xv = np.array(modes["non_torsional_cm1"])*THETA_CM1/temperature
    log_qvib = float(np.sum(-xv/2-np.log(-np.expm1(-xv))))
    hvib = R*temperature*float(np.sum(xv*(0.5+1/np.expm1(xv))))
    masses = np.array(modes["masses_amu"])*AMU_KG
    coordinates = np.array(record["coordinates_A"])*1e-10
    centered = coordinates-np.average(coordinates, axis=0, weights=masses)
    tensor = np.eye(3)*np.sum(masses*np.sum(centered**2, axis=1))
    tensor -= np.einsum("i,ij,ik->jk", masses, centered, centered)
    moments = np.linalg.eigvalsh(tensor)
    theta = HPLANCK**2/(8*math.pi**2*moments*KB)
    log_qrot = 0.5*math.log(math.pi)+1.5*math.log(temperature)
    log_qrot -= math.log(symmetry["external_divisor"])+0.5*np.sum(np.log(theta))
    log_qtrans = 1.5*math.log(2*math.pi*sum(masses)*KB*temperature/HPLANCK**2)
    log_qtrans += math.log(KB*temperature/100000.)
    h = float(htor+hvib+4*R*temperature)
    log_q = float(log_qtor+log_qvib+log_qrot+log_qtrans)
    s = R*log_q+h/temperature
    return {"H_thermal_J_mol": h, "S_J_mol_K": s, "log_q": log_q,
            "effective_samples": effective_samples,
            "log_qtor": log_qtor, "H_tor_J_mol": htor,
            "S_translation_J_mol_K": R*log_qtrans+2.5*R,
            "S_external_rotation_J_mol_K": R*log_qrot+1.5*R,
            "S_non_torsional_vibration_J_mol_K": R*log_qvib+hvib/temperature,
            "S_coupled_rotor_J_mol_K": R*log_qtor+htor/temperature,
            "log_influence": log_influence, "H_influence": h_influence,
            "S_influence": s_influence}


class Ensemble:
    def __init__(self, scratch, name, baseline, samples, seed, proposal_kind="diagonal",samples_by_index=None):
        self.name = name
        selection = json.loads((scratch / "composite" / name / "selection.json").read_text())
        self.indices = selection["unique_minima_indices"]
        counts=samples_by_index or {}
        self.records = [json.loads((scratch / "rotors" / name / f"{i:04d}" / quadrature_filename(counts.get(i,samples),seed,proposal_kind)).read_text())
                        for i in self.indices]
        self.energy = {}
        for level in ("pbe", "blyp"):
            self.energy[level] = np.array([json.loads((scratch / "composite" / name / f"{i:04d}" / level / "result.json").read_text())["energy_hartree"]
                                         for i in self.indices])*EH_MOL
        self.symmetry = external_symmetry(SPECIES[name], baseline,
            self.records[0]["definitions"], self.records[0]["coordinates_A"])

    def hs(self, temperature, level):
        data = [conformer_thermal(record, temperature, self.symmetry) for record in self.records]
        energy = self.energy[level]
        low = float(np.min(energy))
        log_weights = np.array([row["log_q"] for row in data])-(energy-low)/(R*temperature)
        probabilities = np.exp(log_weights-logsumexp(log_weights))
        total_h = np.array([row["H_thermal_J_mol"] for row in data])+energy
        entropy = np.array([row["S_J_mol_K"] for row in data])
        conformational_s = -R*float(np.sum(probabilities*np.log(np.maximum(probabilities, 1e-300))))
        components = {key:float(probabilities@np.array([row[key] for row in data]))
                      for key in ("S_translation_J_mol_K","S_external_rotation_J_mol_K",
                                  "S_non_torsional_vibration_J_mol_K","S_coupled_rotor_J_mol_K")}
        h = float(probabilities@total_h)
        s = float(probabilities@entropy + conformational_s)
        if abs(sum(components.values())+conformational_s-s)>1e-8:
            raise AssertionError("entropy components do not sum to the ensemble entropy")
        variances = np.zeros(3)
        basin_variances=[]
        for p, hi, si, logp, row in zip(probabilities, total_h, entropy, np.log(np.maximum(probabilities, 1e-300)), data):
            influence_h = p*(row["H_influence"]+(hi-h)*row["log_influence"])
            influence_s = p*(row["S_influence"]+(si-R*logp-s)*row["log_influence"])
            influence_g = -R*temperature*p*row["log_influence"]
            n = len(influence_h)
            contribution=np.array([np.var(influence_h, ddof=1), np.var(influence_s, ddof=1),
                                   np.var(influence_g, ddof=1)])/n
            variances += contribution
            basin_variances.append(contribution.tolist())
        return {"T_K": temperature, "H_J_mol": h, "S_J_mol_K": s,
                "H_above_electronic_minimum_kJ_mol": (h-low)/1000,
                "S_conformer_J_mol_K": conformational_s,
                "entropy_components": components,
                "probabilities": probabilities.tolist(),
                "minimum_effective_samples": min(row["effective_samples"] for row in data),
                "MC_basin_variances_H_S_G":basin_variances,
                "MC_standard_errors_H_S_G": np.sqrt(variances).tolist()}


def cycle(ensembles, coefficients, temperature, level):
    h = s = 0.
    variance = np.zeros(3)
    components={"S_basin_mixing_J_mol_K":0.0}
    for name, coefficient in coefficients.items():
        row = ensembles[name].hs(temperature, level)
        h += coefficient*row["H_J_mol"]
        s += coefficient*row["S_J_mol_K"]
        variance += coefficient**2*np.array(row["MC_standard_errors_H_S_G"])**2
        for key,value in row["entropy_components"].items():
            components[key]=components.get(key,0.0)+coefficient*value
        components["S_basin_mixing_J_mol_K"]+=coefficient*row["S_conformer_J_mol_K"]
    if abs(sum(components.values())-s)>1e-8:
        raise AssertionError("cycle entropy components do not sum")
    return {"T_K": temperature, "H_J_mol": h, "S_J_mol_K": s,
            "entropy_components":components,
            "MC_standard_errors_H_S_G": np.sqrt(variance).tolist()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch", type=Path, default=SCRATCH)
    parser.add_argument("--species", choices=tuple(SPECIES), nargs="+")
    parser.add_argument("--samples", type=int, default=4096)
    parser.add_argument("--samples-for",action="append",default=[],metavar="MOLECULE=N",
                        help="increase quadrature for specified molecules without changing the others")
    parser.add_argument("--seed", type=int, default=43003)
    parser.add_argument("--proposal-for",action="append",default=[],metavar="MOLECULE=PROPOSAL",
                        help="select diagonal or correlated importance sampling for specified molecules")
    parser.add_argument("--samples-in",action="append",default=[],metavar="MOLECULE:INDEX=N",
                        help="allocate extra points to individual high-variance wells")
    parser.add_argument("--preflight", action="store_true", help="check rotor precision without computing corrected Tc")
    args = parser.parse_args()
    if args.samples<64:
        raise ValueError("quadrature needs at least 64 points per well")
    physical = {}
    for cpu in sorted(os.sched_getaffinity(0)):
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        core = (topology / "physical_package_id").read_text(), (topology / "core_id").read_text()
        physical.setdefault(core, cpu)
    os.sched_setaffinity(0, sorted(physical.values())[:8])
    scratch = allowed_scratch(args.scratch)
    baseline = json.loads((scratch / "baseline/baseline.json").read_text())
    manifest = json.loads((scratch / "species.json").read_text())
    names = args.species or SPECIES
    sample_counts = {name: args.samples for name in names}
    for option in args.samples_for:
        name,count = option.split("=",1)
        if name not in sample_counts or int(count) < args.samples:
            raise ValueError("invalid per-molecule quadrature count")
        sample_counts[name] = int(count)
    proposals=dict.fromkeys(names,"diagonal")
    for option in args.proposal_for:
        name,proposal=option.split("=",1)
        if name not in proposals or proposal not in ("diagonal","correlated"):
            raise ValueError("invalid per-molecule importance proposal")
        proposals[name]=proposal
    basin_counts={name:{} for name in names}
    for option in args.samples_in:
        well,count=option.split("=",1);name,index=well.split(":",1)
        if name not in basin_counts or int(count)<sample_counts[name]:
            raise ValueError("invalid individual-basin quadrature count")
        basin_counts[name][int(index)]=int(count)
    ensembles = {name: Ensemble(scratch, name, baseline["species"][name], sample_counts[name], args.seed,proposals[name],basin_counts[name]) for name in names}
    if any(set(basin_counts[name])-set(ensembles[name].indices) for name in names):
        raise ValueError("point allocation names an absent well")
    result = {"p0_Pa": 100000., "R_J_mol_K": R, "hartree_J_mol": EH_MOL,
              "samples_per_basin": args.samples, "quadrature_base_seed": args.seed,
              "samples_per_species": sample_counts,
              "importance_proposals":proposals,
              "samples_per_basin_overrides":basin_counts,
              "state": "fixed stereosequence, conformers at equilibrium; no diastereomer mixing term",
              "species": {}, "cycles": {}, "conditional_corrections": {}}
    for name, ensemble in ensembles.items():
        result["species"][name] = {"indices": ensemble.indices, "symmetry": ensemble.symmetry,
            "hard_core_zero_weight_points":sum(sum(inside and value is None for inside,value in
                zip(record["in_basin"],record["energies_above_minimum_J_mol"])) for record in ensemble.records),
            "total_quadrature_points":sum(record["samples"] for record in ensemble.records),
            "minimum_electronic_energy_hartree": {level: float(np.min(energy)/EH_MOL)
                                                 for level,energy in ensemble.energy.items()},
            "racemic_gas_entropy_addition_J_mol_K": 0.0 if manifest[name]["achiral"] else R*math.log(2),
            "electronic_order_indices": {level: [ensemble.indices[i] for i in np.argsort(ensemble.energy[level])]
                                         for level in ("pbe","blyp")},
            "levels": {level: [ensemble.hs(t, level) for t in TEMPERATURES] for level in ("pbe", "blyp")}}
        anchor = result["species"][name]["levels"]["pbe"][0]
        print(f"I043 {name}: PBE-composite S298={anchor['S_J_mol_K']:.3f} J/mol/K; "
              f"thermal H298 above electronic minimum={anchor['H_above_electronic_minimum_kJ_mol']:.3f} kJ/mol; "
              f"minimum effective samples={anchor['minimum_effective_samples']:.1f}", flush=True)
    if set(ensembles) == set(SPECIES):
        sequence_weights = {
            "dyads": {"diphenylpentane_meso": 0.5, "diphenylpentane_racemo": 0.5},
            "triads": {"triphenylheptane_iso": 0.25, "triphenylheptane_syndio": 0.25,
                       "triphenylheptane_hetero": 0.5},
        }
        mirror_s = {side: R*sum(w*math.log(1 if manifest[n]["achiral"] else 2)
                               for n,w in weights.items()) for side,weights in sequence_weights.items()}
        external_s = {side: -R*sum(w*math.log(ensembles[n].symmetry["external_divisor"])
                                  for n,w in weights.items()) for side,weights in sequence_weights.items()}
        class_s = {side: -R*sum(w*math.log(w) for w in weights.values())
                   for side,weights in sequence_weights.items()}
        end_increment = external_s["triads"]-external_s["dyads"]
        class_increment = class_s["triads"]-class_s["dyads"]
        configuration_increment = end_increment+class_increment+mirror_s["triads"]-mirror_s["dyads"]
        if abs(configuration_increment-R*math.log(2)) > 1e-8:
            raise AssertionError("finite-fragment stereo counting does not recover one new-centre bit")
        gav_external_s = {side: -R*sum(w*math.log(
            ensembles[n].symmetry["nonoptical_RMG_sigma"]/
            ensembles[n].symmetry["internal_rotor_divisor"])
            for n,w in weights.items()) for side,weights in sequence_weights.items()}
        gav_optical_s = {side: R*math.log(2)*sum(
            w*baseline["species"][n]["optical_atom_half_factors"]
            for n,w in weights.items()) for side,weights in sequence_weights.items()}
        gav_end_increment = gav_external_s["triads"]-gav_external_s["dyads"]
        gav_optical_increment = gav_optical_s["triads"]-gav_optical_s["dyads"]
        if abs(gav_end_increment+gav_optical_increment-R*math.log(2)) > 1e-8:
            raise AssertionError("capped GAV bookkeeping does not contain one new-centre bit")
        entropy_offset = gav_end_increment+gav_optical_increment-end_increment
        if abs(entropy_offset-class_increment) > 1e-8:
            raise AssertionError("matched QM/GAV normalization disagrees with class counting")
        result["stereo_accounting"] = {
            "weights": sequence_weights, "R_ln_2_J_mol_K": R*math.log(2),
            "excluded_global_mirror_entropy_J_mol_K": mirror_s,
            "global_mirror_increment_if_included_J_mol_K": mirror_s["triads"]-mirror_s["dyads"],
            "external_end_symmetry_increment_J_mol_K": end_increment,
            "excluded_class_mixing_entropy_increment_J_mol_K": class_increment,
            "end_plus_class_plus_mirror_increment_J_mol_K": configuration_increment,
            "GAV_external_end_symmetry_increment_J_mol_K": gav_end_increment,
            "GAV_local_optical_increment_J_mol_K": gav_optical_increment,
            "matched_entropy_offset_J_mol_K": entropy_offset,
            "new_centre_configurational_term": "unchanged RMG term retained once; no additional term applied",
        }
        coefficients = CYCLE_COEFFICIENTS
        for name, values in coefficients.items():
            result["cycles"][name] = {"coefficients": values,
                "levels": {level: [cycle(ensembles, values, t, level) for t in TEMPERATURES] for level in ("pbe", "blyp")}}
        increased=[n for n in names if sample_counts[n]!=args.samples or basin_counts[n]]
        if increased:
            base=dict(ensembles)
            for n in increased:
                base[n]=Ensemble(scratch,n,baseline["species"][n],args.samples,args.seed,proposals[n])
            result["precision_comparison"]={"base_points":args.samples,"increased_species":increased,
                "cycles":{name:{level:[cycle(base,values,t,level) for t in TEMPERATURES]
                    for level in ("pbe","blyp")} for name,values in coefficients.items()}}
        result["stereo_channel_residuals"] = {}
        if sample_counts["diphenylpropane"]!=args.samples:
            alternate=dict(ensembles)
            alternate["diphenylpropane"]=Ensemble(scratch,"diphenylpropane",baseline["species"]["diphenylpropane"],args.samples,args.seed)
            result["quadrature_comparison"]={"species":"diphenylpropane","points":[args.samples,sample_counts["diphenylpropane"]],
                "levels":{level:[cycle(alternate,coefficients["diphenylpropane_isodesmic"],t,level)
                                  for t in TEMPERATURES] for level in ("pbe","blyp")}}
        for dyad,triad in (("meso","iso"),("meso","hetero"),("racemo","syndio"),("racemo","hetero")):
            values = {"diphenylpentane_"+dyad: -1, "triphenylheptane_"+triad: 1,
                      "cumene": -1, "n_propylbenzene": -1, "ethylbenzene": 1, "ethane": 1}
            result["stereo_channel_residuals"][dyad+"_to_"+triad] = {
                "unconditional_atactic_weight": 0.25, "coefficients": values,
                "levels": {level: [cycle(ensembles,values,t,level) for t in TEMPERATURES] for level in ("pbe","blyp")}}
        for level in ("pbe","blyp"):
            for i,t in enumerate(TEMPERATURES):
                for key in ("H_J_mol","S_J_mol_K"):
                    averaged = sum(0.25*channel["levels"][level][i][key]
                                   for channel in result["stereo_channel_residuals"].values())
                    if abs(averaged-result["cycles"]["atactic_added_unit"]["levels"][level][i][key]) > 1e-5:
                        raise AssertionError("conditional channel average differs from atactic residual")
        targets = {"H_standard_error_kJ_mol": 0.5, "S_standard_error_J_mol_K": 1.0,
                   "G_standard_error_kJ_mol": 0.5}
        failures = []
        for name, data in result["cycles"].items():
            for level, rows in data["levels"].items():
                for row in rows:
                    errors = row["MC_standard_errors_H_S_G"]
                    if errors[0] > 500 or errors[1] > 1 or errors[2] > 500:
                        failures.append({"cycle": name, "level": level, "T_K": row["T_K"],
                                         "errors_H_S_G": errors})
        result["numerical_preflight"] = {"targets": targets, "passed": not failures,
                                          "failures": failures, "corrected_Tc_computed": False}
        if args.preflight:
            (scratch / "numerical_preflight.json").write_text(json.dumps(result,indent=2)+"\n")
            print(f"I043 numerical preflight: {'PASS' if not failures else 'MORE QUADRATURE NEEDED'}; {len(failures)} target misses; no corrected Tc computed",flush=True)
            return
        if failures:
            raise RuntimeError("numerical precision targets missed; run --preflight and increase quadrature before computing Tc")
        from rmgpy.data.rmg import RMGDatabase
        from rmgpy.molecule import Molecule
        from rmgpy.reaction import Reaction
        from rmgpy.species import Species
        db = RMGDatabase()
        db.load_thermo(str(scratch / "baseline/database/input/thermo"), thermo_libraries=["primaryThermoLibrary"])
        graphs = baseline["propagation_graphs"]
        sides = {}
        for side in ("reactants", "products"):
            sides[side] = [Species(molecule=[Molecule().from_adjacency_list(adj)]) for adj in graphs[side]]
            for item in sides[side]:
                item.thermo = db.thermo.get_thermo_data(item)
        propagation = Reaction(**sides)
        for name, values in coefficients.items():
            result["conditional_corrections"][name] = {}
            for level in ("pbe", "blyp"):
                def delta(t):
                    measured = cycle(ensembles, values, t, level)
                    if name == "diphenylpropane_isodesmic":
                        model = baseline["cycles"][name]["thermo"][0]
                        return -(measured["H_J_mol"]-model["H_J_mol"]), -(measured["S_J_mol_K"]-model["S_J_mol_K"])
                    # Normalize both capped QM and GAV sides to an oriented
                    # fixed sequence. GAV already represents the new-centre bit
                    # through end/optical bookkeeping; the compiler keeps it.
                    model = baseline["cycles"][name]["thermo"][0]
                    return measured["H_J_mol"]-model["H_J_mol"], measured["S_J_mol_K"]-model["S_J_mol_K"]+entropy_offset
                def free_energy(t):
                    dh, ds = delta(t)
                    return propagation.get_enthalpy_of_reaction(t)+dh-t*(propagation.get_entropy_of_reaction(t)+ds)-R*t*math.log(1000*R*t/100000)
                roots = []
                grid = np.linspace(298.15, 800, 102)
                for a, b in zip(grid[:-1], grid[1:]):
                    if free_energy(a)*free_energy(b) < 0:
                        roots.append(float(brentq(free_energy, a, b, xtol=1e-7)))
                dh, ds = delta(298.15)
                row = {"delta_H298_kJ_mol": dh/1000, "delta_S298_J_mol_K": ds,
                       "roots_298_to_800_K": roots,
                       "propagation_H298_kJ_mol": (baseline["propagation"][0]["H_J_mol"]+dh)/1000,
                       "propagation_S298_J_mol_K": baseline["propagation"][0]["S_J_mol_K"]+ds,
                       "QM_removed_end_symmetry_entropy_J_mol_K": 0.0 if name=="diphenylpropane_isodesmic" else end_increment,
                       "GAV_removed_end_symmetry_entropy_J_mol_K": 0.0 if name=="diphenylpropane_isodesmic" else gav_end_increment,
                       "matched_entropy_offset_J_mol_K": 0.0 if name=="diphenylpropane_isodesmic" else entropy_offset,
                       "conditional_transfer": "closed-shell local model; both QM and GAV normalized to oriented fixed sequences; RMG R ln 2 retained once"}
                measured = cycle(ensembles, values, 298.15, level)
                row["MC_standard_errors_delta_H298_kJ_S298_J"] = [
                    measured["MC_standard_errors_H_S_G"][0]/1000,
                    measured["MC_standard_errors_H_S_G"][1]]
                row["root_MC_standard_errors_K"] = []
                for root in roots:
                    slope = (free_energy(root+0.01)-free_energy(root-0.01))/0.02
                    root_cycle = cycle(ensembles, values, root, level)
                    row["root_MC_standard_errors_K"].append(
                        root_cycle["MC_standard_errors_H_S_G"][2]/abs(slope))
                result["conditional_corrections"][name][level] = row
                print(f"I043 {name} {level}: delta H298={dh/1000:.3f} kJ/mol; delta S298={ds:.3f} J/mol/K; roots={roots}", flush=True)
        result["method_spread"] = {}
        for name, levels in result["conditional_corrections"].items():
            spread = {}
            for key in ("delta_H298_kJ_mol", "delta_S298_J_mol_K"):
                values = [levels[level][key] for level in ("pbe", "blyp")]
                spread[key] = {"minimum": min(values), "maximum": max(values),
                               "half_range": (max(values)-min(values))/2}
            roots = [levels[level]["roots_298_to_800_K"] for level in ("pbe", "blyp")]
            if all(len(value) == 1 for value in roots):
                values = [value[0] for value in roots]
                spread["Tc_K"] = {"minimum": min(values), "maximum": max(values),
                                  "half_range": (max(values)-min(values))/2}
            result["method_spread"][name] = spread
        result["numerical_preflight"]["corrected_Tc_computed"]=True
    suffix = "thermochemistry.json" if args.species is None else "thermochemistry_partial.json"
    (scratch / suffix).write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
