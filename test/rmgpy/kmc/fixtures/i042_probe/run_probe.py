"""Generate I042 evidence using only allowlisted thermo/solvation database paths.

Command (from repository root; persist both streams as in the report):
PYTHONPATH=$PWD /home/alon/anaconda3/envs/rmg_env/bin/python
test/rmgpy/kmc/fixtures/i042_probe/run_probe.py
--scratch /home/alon/runs/i042-melt-reference/reproduce
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import shutil
import sys

import numpy as np
from scipy.optimize import brentq
from rmgpy import constants
from rmgpy.data.rmg import RMGDatabase
from rmgpy.kmc.compiler import ps_proxy_set

import models as m

HERE = Path(__file__).resolve().parent
FIXTURES = HERE.parent
REPO = HERE.parents[4]
spec = importlib.util.spec_from_file_location("i034", FIXTURES / "i034_probe/run_probe.py")
prior = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)
TEMPERATURES = (600.0, 700.0, 800.0)


def progress(message):
    print("[I042] " + message, file=sys.stderr, flush=True)


def crossing(function, lower=300., upper=1000.):
    grid = np.linspace(lower, upper, 71)
    roots = []
    for a, b in zip(grid[:-1], grid[1:]):
        if function(a)*function(b) < 0:
            roots.append(float(brentq(function, a, b, xtol=1e-8)))
    return roots


def gas_hs(reaction, t):
    h = reaction.get_enthalpy_of_reaction(t) + m.R*t
    s = reaction.get_entropy_of_reaction(t) + m.R*(1+math.log(m.C0*m.R*t/1e5))
    return float(h), float(s)


def safe_scratch(path):
    path = path.resolve()
    if path.is_relative_to(REPO) or path.is_relative_to(Path("/home/alon/Code/RMG-database")):
        raise ValueError("output must be outside repositories")
    if "catalog" in path.parts or path.is_relative_to(Path("/home/alon/Code/polymers")):
        raise ValueError("excluded path")
    return path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scratch", type=Path, required=True)
    args = parser.parse_args()
    scratch = safe_scratch(args.scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    snapshot_path = scratch / "database"
    progress("materializing pinned allowlist with git show")
    snapshot = prior.snapshot_database(Path("/home/alon/Code/RMG-database"), snapshot_path)
    db = RMGDatabase()
    db.load_kinetics(str(snapshot_path / "input/kinetics"), reaction_libraries=[],
                     seed_mechanisms=None, kinetics_families=list(prior.PS_FAMILY_CANDIDATES),
                     kinetics_depositories=["training"])
    db.load_thermo(str(snapshot_path / "input/thermo"),
                   thermo_libraries=["primaryThermoLibrary"], depository=True)
    db.load_solvation(str(snapshot_path / "input/solvation"))
    progress("reproducing the existing I034/I039 propagation pair")
    reactions, rates = prior.gas_baseline(db, ps_proxy_set(3))
    reaction = reactions["propagation"]
    entries, failures = prior.estimate_species(db, {"propagation": reaction})
    assert not failures
    assert constants.R == m.R and constants.Na == m.NA
    signs = [-1]*len(reaction.reactants)+[1]*len(reaction.products)
    participants = reaction.reactants+reaction.products
    species = [{"smiles": s.molecule[0].to_smiles(),
                "adjacency": s.molecule[0].to_adjacency_list(),
                "formula": s.molecule[0].get_element_count(), "sign": sign,
                "solute": entries[s.molecule[0].to_smiles()]["descriptors"]}
               for s, sign in zip(participants, signs)]
    radical_species = [s for s in species if s["formula"]["C"] != 8]
    masses = [s["formula"]["C"]*12.011+s["formula"]["H"]*1.008 for s in radical_species]
    assert math.isclose(masses[1]-masses[0], m.MW, abs_tol=1e-10)
    gas_g = lambda t: gas_hs(reaction, t)[0]-t*gas_hs(reaction, t)[1]
    gas_root = crossing(gas_g)[0]
    progress("PC-SAFT residual chemical potentials at dispatch volume")
    pc = {}
    for label in m.PC_SETS:
        state = m.pc_state(label, masses)
        rows = []
        for t in TEMPERATURES:
            row = m.pc_transfer(t, state)
            h, s = gas_hs(reaction, t)
            row.update(T_K=t, corrected_H_path=h+row["H_path"],
                       corrected_S_path=s+row["S_path"],
                       corrected_H_pressure=h+row["H_pressure"],
                       corrected_S_pressure=s+row["S_pressure"])
            row["monomer_concentration_activity"] = math.exp(row["mu_J_mol"][1]/(m.R*t))
            row["pair_equilibrium_multiplier"] = math.exp(-row["G"]/(m.R*t))
            # Local composition stability of the carrier/monomer background.
            c = state[0]
            jac = np.zeros((2, 2))
            for j in range(2):
                step = c[j]*1e-4
                def mu_at(value):
                    trial = c.copy(); trial[j] = value
                    return m.chemical_potentials(t, trial, *state[1:])[:2]
                jac[:, j] = m.derivative(mu_at, c[j], step)
            jac += np.diag(1/c[:2])
            row["stability_min_eigenvalue_m3_mol"] = float(np.linalg.eigvalsh(jac).min())
            rows.append(row)
        roots = crossing(lambda t: gas_g(t)+m.pc_values(t, state)["G"])
        # The isobar is a state control, not the fixed-volume kMC calculation.
        pressure_rows = []
        for t in TEMPERATURES:
            density = brentq(lambda rho: m.pc_values(t, m.pc_state(label, masses, density=rho))["P_Pa"]-1e5,
                             600, 1150)
            row = m.pc_transfer(t, m.pc_state(label, masses, density=density))
            pressure_rows.append({"T_K": t, "density_kg_m3": density,
                                  "G": row["G"], "H_pressure": row["H_pressure"],
                                  "S_pressure": row["S_pressure"], "P_Pa": row["P_Pa"]})
        # Controls selected by state variables/chain size, never by ceiling distance.
        controls = []
        for density, dp in ((900., m.DP), (1050., 100.), (1050., 10000.)):
            control_state = m.pc_state(label, masses, density=density, dp=dp)
            roots_control = crossing(lambda t: gas_g(t)+m.pc_values(t, control_state)["G"])
            controls.append({"density_kg_m3": density, "DP": dp,
                             "G700": m.pc_values(700, control_state)["G"], "roots_K": roots_control})
        pc[label] = {"parameters": m.PC_SETS[label], "rows": rows, "roots_K": roots,
                     "shifts_K": [r-gas_root for r in roots], "isobar_1bar": pressure_rows,
                     "controls": controls}
    progress("offline UNIFAC and FH activity diagnostics")
    us = m.unifac_state()
    unifac = []
    for t in TEMPERATURES:
        row = m.unifac_values(t, us)
        row["S_mixing_1M"] = -m.derivative(lambda u: m.unifac_values(u, us)["G_mixing_1M"], t, 0.05)
        row["H_mixing_1M"] = row["G_mixing_1M"]+t*row["S_mixing_1M"]
        phi = m.C0/m.CTOT
        row["effective_chi"] = (math.log(row["monomer_activity"]/phi)-(1-1/m.DP)*(1-phi))/(1-phi)**2
        row["T_K"] = t
        unifac.append(row)
    phi = m.C0/m.CTOT
    fh_coefficient = math.log(3/2)+math.log(m.CTOT/m.C0)-1
    fh = [{"T_K": t, "G_mixing_1M": m.R*t*fh_coefficient,
           "H_mixing_1M": 0., "S_mixing_1M": -m.R*fh_coefficient,
           "monomer_activity": m.fh_activity(t, m.C0),
           "dG_dchi_J_mol": m.R*t*(2*phi-1),
           "dH_dchi_b_J_mol_K": m.R*(2*phi-1)} for t in TEMPERATURES]
    concentrations = []
    for f in (0.001, 0.01, m.C0/m.CTOT, 0.1, 0.5):
        c = f*m.CTOT
        root = crossing(lambda t: gas_g(t)-m.R*t*math.log(c/m.C0))[0]
        concentrations.append({"fraction_free_styrene": f, "c_mol_L": c/1000,
                               "gas_Tc_K": root, "gas_shift_K": root-gas_root,
                               "FH_athermal_activity": m.fh_activity(700, c),
                               "UNIFAC_activity700": m.unifac_values(700, m.unifac_state(c))["monomer_activity"],
                               "concentration_S_shift": m.R*math.log(c/m.C0)})
    equilibrium = [{"T_K": t, "gas_equilibrium_c_mol_L": 1/reaction.get_equilibrium_constant(t)/1000}
                   for t in TEMPERATURES]
    progress("RMG alternative molecular solvents and validity guards")
    solvents = {}
    for label in prior.CANDIDATES:
        solvent = db.solvation.get_solvent_data(label)
        coefficients = {field: getattr(solvent, field) for field in (*prior.FIELDS_G, *prior.FIELDS_H)}
        if any(value is None for value in coefficients.values()):
            solvents[label] = {"status": "missing LSER coefficients", "coefficients": coefficients,
                              "critical_K": prior.get_critical_temperature(solvent.name_in_coolprop)
                              if solvent.name_in_coolprop else None}
            continue
        anchors = {key: db.solvation.get_solvation_correction(e["solute"], solvent) for key, e in entries.items()}
        dh = prior.reaction_delta(reaction, {key: v.enthalpy for key, v in anchors.items()})
        ds = prior.reaction_delta(reaction, {key: v.entropy for key, v in anchors.items()})
        roots = crossing(lambda t: gas_g(t)+dh-t*ds)
        kfactor = []
        for t in TEMPERATURES:
            try:
                def g(u):
                    return prior.reaction_delta(reaction, {key: db.solvation.get_T_dep_solvation_energy_from_LSER_298(
                        e["solute"], solvent, u)[0] for key, e in entries.items()})
                dg = g(t)
                entropy = -m.derivative(g, t, 0.05)
                kfactor.append({"T_K": t, "G": dg, "S": entropy, "H": dg+t*entropy})
            except Exception as error:
                kfactor.append({"T_K": t, "error": type(error).__name__+": "+str(error)})
        critical = prior.get_critical_temperature(solvent.name_in_coolprop) if solvent.name_in_coolprop else None
        kfactor_roots = crossing(lambda t: gas_g(t)+g(t), 300, critical-0.1) if critical else None
        solvents[label] = {"H": dh, "S": ds, "rows": [{"T_K": t, "G": dh-t*ds} for t in TEMPERATURES],
                          "roots_K": roots, "shifts_K": [r-gas_root for r in roots],
                          "critical_K": critical, "kfactor_roots_K": kfactor_roots,
                          "kfactor": kfactor,
                          "coefficients": coefficients}
    inventory = prior.inventory(db.solvation, snapshot_path)
    hash_paths = ["rmgpy/kmc/compiler.py", "rmgpy/kmc/ssa.py", "rmgpy/data/solvation.py", "rmgpy/constants.py",
                  "test/rmgpy/kmc/fixtures/i034_probe/run_probe.py", "test/rmgpy/kmc/fixtures/i042_probe/models.py",
                  "test/rmgpy/kmc/fixtures/i042_probe/run_probe.py"]
    result = {"provenance": {"database_sha": prior.DATABASE_SHA, "snapshot": snapshot,
                             "source_sha256": {p: hashlib.sha256((REPO/p).read_bytes()).hexdigest() for p in hash_paths}},
              "definitions": {"R": m.R, "NA": m.NA, "MW_g_mol": m.MW, "density_kg_m3": m.RHO,
                              "repeat_c_mol_L": m.CTOT/1000, "repeat_volume_m3": 1/(m.NA*m.CTOT),
                              "DP": m.DP, "c0_mol_m3": m.C0, "radical_masses_g_mol": masses},
              "species": species, "rates": {"propagation": rates["propagation_grid"],
                                            "depropagation": rates["depropagation_grid"]},
              "gas": {"continuous_K": gas_root,
                      "rows": [{"T_K": t, "H": gas_hs(reaction, t)[0], "S": gas_hs(reaction, t)[1], "G": gas_g(t)}
                               for t in TEMPERATURES]},
              "PC": pc, "UNIFAC": unifac, "FH_athermal": fh,
              "concentrations": concentrations, "equilibrium": equilibrium,
              "solvents": solvents, "radical_audit": prior.audit_radicals(db.solvation, entries),
              "solvent_inventory_count": len(inventory["rows"]),
              "polymer_labels": [r["label"] for r in inventory["rows"] if "polystyrene" in r["label"].lower()],
              "availability": {x: shutil.which(x) for x in ("amspython", "COSMOtherm", "cosmors")}}
    result["literature_domains"] = {
        "PC_A_styrene_vapor_pressure_fit_K": [-30+273.15, 363+273.15],
        "PC_A_Sty_PS_kij_calibration_K": 65+273.15,
        "Gornert_Sadowski_ternary_T_K": 338.15,
        "Gornert_Sadowski_ternary_P_MPa": [10., 15.],
        "Gornert_Sadowski_ternary_carrier_kg_mol": [6., 105.],
    }
    output = scratch / "results.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print("I042 results written to " + str(output))
    for label, values in pc.items():
        print(label + " continuous roots K: " + str(values["roots_K"]))


if __name__ == "__main__":
    main()
