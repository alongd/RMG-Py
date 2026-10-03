"""Reproduce the pinned RMG reference and balanced diagnostic cycles.

Run with rmg_env; --scratch is constrained to the dispatch scratch directory.
Reuse I034's allowlist materializer and I039's propagation selection/conventions.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
from pathlib import Path

from scipy.optimize import brentq
from rmgpy import constants
from rmgpy.data.rmg import RMGDatabase
from rmgpy.molecule import Molecule
from rmgpy.species import Species

from run_ensembles import SPECIES, allowed_scratch

FIXTURES = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("i039", FIXTURES / "i039_probe/run_probe.py")
i039 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(i039)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch", type=Path, required=True)
    args = parser.parse_args()
    physical = {}
    for cpu in sorted(os.sched_getaffinity(0)):
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        core = (topology / "physical_package_id").read_text(), (topology / "core_id").read_text()
        physical.setdefault(core, cpu)
    os.sched_setaffinity(0, sorted(physical.values())[:8])
    scratch = allowed_scratch(args.scratch)
    snapshot = i039.prior.snapshot_database(Path("/home/alon/Code/RMG-database"), scratch / "database")
    db = RMGDatabase()
    db.load_thermo(str(scratch / "database/input/thermo"), thermo_libraries=["primaryThermoLibrary"])
    db.load_kinetics(str(scratch / "database/input/kinetics"), reaction_libraries=[],
                    kinetics_families=list(i039.prior.PS_FAMILY_CANDIDATES),
                    kinetics_depositories=["training"])
    propagation, compiled = i039.compile_pair(db, 3)
    for item in propagation.reactants + propagation.products:
        item.thermo = db.thermo.get_thermo_data(item)
    temperatures = [298.15, 400., 500., 600., 700., 800.]
    molecules = {}
    for name, smiles in SPECIES.items():
        item = Species(molecule=[Molecule(smiles=smiles)])
        item.generate_resonance_structures()
        item.thermo = db.thermo.get_thermo_data(item)
        molecules[name] = item
    cycles = {
        "diphenylpropane_isodesmic": {"diphenylpropane": -1, "ethane": -1,
                                      "ethylbenzene": 1, "n_propylbenzene": 1},
        "meso_iso_increment": {"diphenylpentane_meso": -1, "cumene": -1,
                               "n_propylbenzene": -1, "triphenylheptane_iso": 1,
                               "ethylbenzene": 1, "ethane": 1},
        "atactic_added_unit": {"diphenylpentane_meso": -0.5, "diphenylpentane_racemo": -0.5,
            "triphenylheptane_iso": 0.25, "triphenylheptane_syndio": 0.25,
            "triphenylheptane_hetero": 0.5, "cumene": -1, "n_propylbenzene": -1,
            "ethylbenzene": 1, "ethane": 1},
    }
    data = {"snapshot": snapshot, "R": constants.R, "p0_Pa": 100000.,
            "c0_mol_m3": 1000., "compiled_grid_Tc_K": compiled["ceiling_tabulated_K"],
            "continuous_Tc_K": brentq(lambda t: math.log(propagation.get_equilibrium_constant(t)*1000), 250, 1800),
            "propagation": [{"T_K": t, "H_J_mol": propagation.get_enthalpy_of_reaction(t),
                             "S_J_mol_K": propagation.get_entropy_of_reaction(t)} for t in temperatures],
            "propagation_graphs": {"reactants": [s.molecule[0].to_adjacency_list() for s in propagation.reactants],
                                   "products": [s.molecule[0].to_adjacency_list() for s in propagation.products]},
            "species": {}, "cycles": {}}
    for name, item in molecules.items():
        decomposition = i039.decompose(db.thermo, item)
        data["species"][name] = {
            "adjacency": item.molecule[0].to_adjacency_list(),
            "sigma": float(item.get_symmetry_number()),
            "optical_atom_half_factors": decomposition["optical_atom_half_factors"],
            "source_weights": decomposition["source_weights"],
            "thermo_comment": item.thermo.comment,
            "thermo": [{"T_K": t, "H_J_mol": item.thermo.get_enthalpy(t),
                        "S_J_mol_K": item.thermo.get_entropy(t)} for t in temperatures]}
    for name, coefficients in cycles.items():
        vector = {}
        for species, coefficient in coefficients.items():
            for group, weight in data["species"][species]["source_weights"].items():
                vector[group] = vector.get(group, 0) + coefficient*weight
        vector = {key: value for key, value in vector.items() if value}
        if vector:
            raise AssertionError(f"nonadditivity cycle did not cancel: {vector}")
        data["cycles"][name] = {"coefficients": coefficients, "source_vector": vector,
            "thermo": [{"T_K": t,
                        "H_J_mol": sum(c*molecules[s].thermo.get_enthalpy(t) for s,c in coefficients.items()),
                        "S_J_mol_K": sum(c*molecules[s].thermo.get_entropy(t) for s,c in coefficients.items())}
                       for t in temperatures]}
    (scratch / "baseline.json").write_text(json.dumps(data, indent=2) + "\n")
    print(f"I043 pinned snapshot: {snapshot['files']} files, {snapshot['sha256']}")
    print(f"I043 compiled baseline: Tc={data['continuous_Tc_K']:.6f} K; "
          f"H298={data['propagation'][0]['H_J_mol']/1000:.6f} kJ/mol; "
          f"S298={data['propagation'][0]['S_J_mol_K']:.6f} J/mol/K")
    for name, cycle in data["cycles"].items():
        row = cycle["thermo"][0]
        print(f"I043 {name}: zero net additive source vector; "
              f"H298={row['H_J_mol']/1000:.6f} kJ/mol; S298={row['S_J_mol_K']:.6f} J/mol/K")


if __name__ == "__main__":
    main()
