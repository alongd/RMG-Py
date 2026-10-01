"""Run the diagnosis-only melt stand-in probe against a read-only database snapshot.

Command: PYTHONPATH=$PWD /home/alon/anaconda3/envs/rmg_env/bin/python
test/rmgpy/kmc/fixtures/i034_probe/run_probe.py --output /tmp/i034-reproduce/results.json
Persist stdout and stderr separately as shown in the accompanying report.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import subprocess
import sys
from pathlib import Path

import CoolProp
from CoolProp.CoolProp import PropsSI
from scipy.optimize import brentq

from rmgpy import constants
from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.solvation import SoluteData, get_critical_temperature
from rmgpy.kmc.compiler import (
    DEFAULT_T_GRID,
    PS_FAMILY_CANDIDATES,
    PS_PROXY_UNITS,
    EventSetCompiler,
    _orient_to_proxy,
    _reverse_view,
    ceiling_temperature,
    ps_proxy_set,
)
from rmgpy.molecule.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species


DATABASE_SHA = "4a12d36fcdc193ede82c8d1ab5c1653495d445bc"
BASE_SHA = "ff8e75ce48b0ff4e5ce8a3a28b985932f784149c"
PRIOR_SHA = "413887456db0502e6a452da74d0ed4d85127bb53"
TEMPERATURES = (600.0, 700.0, 800.0)
CANDIDATES = ("benzene", "toluene", "ethylbenzene", "dodecane", "hexadecane")
FIELDS_G = tuple(f"{prefix}_g" for prefix in ("s", "b", "e", "l", "a", "c"))
FIELDS_H = tuple(f"{prefix}_h" for prefix in ("s", "b", "e", "l", "a", "c"))
SOLUTE_FIELDS = ("S", "B", "E", "L", "A", "V")
MONOMER_CONCENTRATION = 1000.0
HEXADECANE_TC = {"value_K": 722.0, "uncertainty_K": 4.0, "source": "NIST SRD 69, C544763, phase-change table, average of nine values; accessed 2026-09-30"}


def progress(message):
    print(f"[I-034] {message}", file=sys.stderr, flush=True)


def git_output(repository, *arguments):
    return subprocess.check_output(["git", "-C", str(repository), *arguments])


def snapshot_database(database, scratch):
    prefixes = [
        "input/solvation",
        "input/thermo/groups",
        "input/thermo/depository",
        "input/thermo/libraries/primaryThermoLibrary.py",
        "input/kinetics/families/recommended.py",
        *[f"input/kinetics/families/{family}" for family in PS_FAMILY_CANDIDATES],
    ]
    listing = git_output(database, "ls-tree", "-r", "--name-only", DATABASE_SHA, "--", *prefixes)
    paths = listing.decode().splitlines()
    digest = hashlib.sha256()
    for relative in paths:
        if "catalog" in Path(relative).parts:
            raise RuntimeError("refusing a catalog path")
        content = git_output(database, "show", f"{DATABASE_SHA}:{relative}")
        target = scratch / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(content)
        digest.update(relative.encode() + b"\0" + content)
    return {"files": len(paths), "sha256": digest.hexdigest(), "prefixes": prefixes}


def smiles(species):
    return [item.molecule[0].to_smiles() for item in species]


def generate(database, reactants, family):
    return database.kinetics.generate_reactions_from_families(
        [item.copy(deep=True) for item in reactants], only_families=[family], resonance=True
    )


def reaction_from_record(record):
    return Reaction(
        reactants=[Species(molecule=[Molecule().from_adjacency_list(graph)]) for graph in record["reactant_graphs"]],
        products=[Species(molecule=[Molecule().from_adjacency_list(graph)]) for graph in record["product_graphs"]],
    )


def benzylic_radical(species):
    for molecule in species.molecule:
        for atom in molecule.atoms:
            if atom.radical_electrons and not any(bond.is_benzene() for bond in atom.bonds.values()):
                if any(any(bond.is_benzene() for bond in neighbor.bonds.values()) for neighbor in atom.bonds):
                    return True
    return False


def require_close(actual, expected, tolerance=6e-6):
    if not math.isclose(actual, expected, rel_tol=tolerance, abs_tol=1e-10):
        raise AssertionError(f"gas baseline mismatch: {actual} != {expected}")


def gas_baseline(database, proxies):
    by_site = {proxy.site_type: proxy for proxy in proxies}
    pristine = by_site["pristine"]
    dissociations = generate(database, pristine.reactants, "R_Recombination")
    target = {"[CH](CCc1ccccc1)c1ccccc1", "[CH2]C(C)c1ccccc1"}
    association = next(reaction for reaction in dissociations if set(smiles(reaction.reactants)) == target)
    homolysis, orientation = _orient_to_proxy(pristine, association)
    abstractions = generate(database, by_site["end_radical+styrene"].reactants, "H_Abstraction")
    abstraction = next(reaction for reaction in abstractions if "C=Cc1[c]cccc1" in smiles(reaction.products))
    chain_pair = (by_site["end_radical"].reactants[0], pristine.reactants[0])
    chain_abstractions = generate(database, chain_pair, "H_Abstraction")
    donor_smiles = pristine.reactants[0].molecule[0].to_smiles()
    chain_choices = [
        reaction for reaction in chain_abstractions
        if any(benzylic_radical(species) for species in reaction.products)
        and donor_smiles not in smiles(reaction.products)
    ]
    chain_abstraction = min(chain_choices, key=lambda reaction: tuple(sorted(smiles(reaction.products))))
    compiler = EventSetCompiler(database.kinetics, (), PS_FAMILY_CANDIDATES,
                                thermo_database=database.thermo, temperature_grid=TEMPERATURES)
    homolysis_rates, homolysis_source = compiler._rate_table(homolysis)
    h_forward, _ = compiler._rate_table(abstraction)
    h_reverse, h_source = compiler._rate_table(_reverse_view(abstraction))
    if homolysis_rates is None or h_forward is None or h_reverse is None:
        raise AssertionError("baseline rate estimation failed")
    require_close(homolysis_source["equilibrium_constant_table"]["Kc"][1], 1.63104e14)
    require_close(homolysis_rates["k"][1], 2.02492e-9)
    require_close(h_source["equilibrium_constant_table"]["Kc"][1], 1.25830e-4)
    progress("gas homolysis and H-abstraction baselines reproduced; compiling the addition-only ceiling oracle")
    ceiling_proxies = tuple(proxy for proxy in proxies if proxy.site_type in ("end_radical", "end_radical+styrene"))
    active, excluded, cache = EventSetCompiler.discover_family_reactions(
        database.kinetics, ceiling_proxies, ("R_Addition_MultipleBond",), family_universe=["R_Addition_MultipleBond"]
    )
    artifact = EventSetCompiler(database.kinetics, ceiling_proxies, active, excluded_families=excluded,
                               thermo_database=database.thermo, reaction_cache=cache,
                               rmg_database_sha=DATABASE_SHA).compile()
    require_close(artifact["ps_ceiling_temperature_K"], 710.2487, tolerance=1e-7)
    selected = artifact["ps_ceiling_pairs"][0]
    prop_record = next(record for record in artifact["records"] if record["event_id"] == selected["propagation_event_id"])
    dep_record = next(record for record in artifact["records"] if record["event_id"] == selected["depropagation_event_id"])
    propagation = reaction_from_record(prop_record)
    rates = {
        "homolysis_s-1": homolysis_rates["k"],
        "Kc_association": homolysis_source["equilibrium_constant_table"]["Kc"],
        "H_forward_m3_mol_s": h_forward["k"],
        "H_reverse_m3_mol_s": h_reverse["k"],
        "H_Kc": h_source["equilibrium_constant_table"]["Kc"],
        "ceiling_tabulated_K": artifact["ps_ceiling_temperature_K"],
        "ceiling_pair_count": len(artifact["ps_ceiling_pairs"]),
        "ceiling_pair": selected,
        "propagation_grid": prop_record["k_table"],
        "depropagation_grid": dep_record["k_table"],
        "orientation": orientation,
        "forced_pristine_dissociations": len(dissociations),
        "chain_abstraction_count": len(chain_abstractions),
    }
    return {"homolysis": Reaction(reactants=homolysis.reactants, products=homolysis.products), "H_baseline": abstraction,
            "H_chain": chain_abstraction, "propagation": propagation}, rates


def estimate_species(database, reactions):
    entries = {}
    failures = []
    for name, reaction in reactions.items():
        for species in reaction.reactants + reaction.products:
            species.thermo = database.thermo.get_thermo_data(species)
            key = species.molecule[0].to_smiles()
            if key in entries:
                continue
            try:
                solute = database.solvation.get_solute_data(species)
                group_species = species.copy(deep=True)
                pure_groups = database.solvation.get_solute_data_from_groups(group_species)
                pure_groups.set_mcgowan_volume(group_species)
                descriptors = {field: float(getattr(solute, field)) for field in SOLUTE_FIELDS}
                group_descriptors = {field: float(getattr(pure_groups, field)) for field in SOLUTE_FIELDS}
                if not all(math.isfinite(value) for value in descriptors.values()):
                    raise ValueError("nonfinite descriptor")
                entries[key] = {"smiles": key, "adjacency": species.molecule[0].to_adjacency_list(remove_h=False),
                                "descriptors": descriptors, "comment": solute.comment,
                                "group_descriptors": group_descriptors, "group_comment": pure_groups.comment,
                                "thermo_comment": species.thermo.comment, "solute": solute}
            except Exception as error:
                failures.append({"reaction": name, "smiles": key, "error": str(error)})
    return entries, failures


def inventory(solvation, scratch):
    rows = []
    reference = SoluteData(S=0.5, B=0.0, E=0.5, L=2.0, A=0.0, V=1.0)
    for label, entry in sorted(solvation.libraries["solvent"].entries.items()):
        data = entry.data
        missing_g = [field for field in FIELDS_G if getattr(data, field) is None]
        missing_h = [field for field in FIELDS_H if getattr(data, field) is None]
        critical = None
        model_error = None
        if data.name_in_coolprop:
            critical = get_critical_temperature(data.name_in_coolprop)
        if not missing_g and not missing_h and data.name_in_coolprop:
            try:
                solvation.get_T_dep_solvation_energy_from_LSER_298(reference, data, 298.0)
            except Exception as error:
                model_error = f"{type(error).__name__}: {error}"
        rows.append({"label": label, "missing_g": missing_g, "missing_h": missing_h,
                     "coolprop": data.name_in_coolprop, "Tc_K": critical, "model_error_298": model_error})
    solvent_text = (scratch / "input/solvation/libraries/solvent.py").read_text()
    keyword_lines = [{"line": number, "text": line.strip()} for number, line in enumerate(solvent_text.splitlines(), 1)
                     if any(word in line.lower() for word in ("polymer", "polystyrene", "melt"))]
    return {"rows": rows, "keyword_lines_in_solvent_library": keyword_lines}


def audit_radicals(solvation, entries):
    rows = []
    for key, entry in sorted(entries.items()):
        molecule = Molecule().from_adjacency_list(entry["adjacency"])
        if not molecule.is_radical():
            continue
        saturated = molecule.copy(deep=True)
        saturated.saturate_radicals()
        stable = solvation.get_solute_data(Species(molecule=[saturated]))
        differences = {field: entry["descriptors"][field] - float(getattr(stable, field))
                       for field in SOLUTE_FIELDS}
        for atom in molecule.atoms:
            if not atom.radical_electrons:
                continue
            node = solvation.groups["radical"].descend_tree(molecule, {"*": atom}, None)
            observation = {"smiles": key, "saturated_smiles": saturated.to_smiles(),
                           "matched_node": node.label if node else None,
                           "radical_minus_saturated": differences}
            try:
                correction = solvation._add_group_solute_data(
                    SoluteData(S=0.0, B=0.0, E=0.0, L=0.0, A=0.0),
                    solvation.groups["radical"], molecule, {"*": atom})
                observation["correction"] = {field: float(getattr(correction, field)) for field in SOLUTE_FIELDS if field != "V"}
            except KeyError as error:
                observation["correction_error"] = str(error)
            rows.append(observation)
    return rows


def collect_boundaries(solvation, solvent_data, solute):
    observations = {}
    if solvent_data.name_in_coolprop:
        critical = get_critical_temperature(solvent_data.name_in_coolprop)
        triple = PropsSI("Ttriple", solvent_data.name_in_coolprop)
        tests = {"below_triple": triple - 1.0, "below_298": 297.0,
                 "near_Tc": critical - 0.01, "at_Tc": critical, "above_Tc": critical + 1.0}
    else:
        tests = {"298": 298.0}
    for label, temperature in tests.items():
        try:
            gibbs, kfactor, henry = solvation.get_T_dep_solvation_energy_from_LSER_298(solute, solvent_data, temperature)
            observations[label] = {"T_K": temperature, "gibbs_J_mol": gibbs, "Kfactor": kfactor, "henry_Pa_m3_mol": henry}
        except Exception as error:
            observations[label] = {"T_K": temperature, "error_type": type(error).__name__, "error": str(error)}
    try:
        correction = solvation.get_solvation_correction(solute, solvent_data)
        observations["linear_at_800"] = {"T_K": 800.0, "gibbs_J_mol": correction.enthalpy - 800.0 * correction.entropy}
    except Exception as error:
        observations["linear_at_800"] = {"T_K": 800.0, "error_type": type(error).__name__, "error": str(error)}
    return observations


def reaction_delta(reaction, values):
    return (sum(values[species.molecule[0].to_smiles()] for species in reaction.products)
            - sum(values[species.molecule[0].to_smiles()] for species in reaction.reactants))


def correction_table(solvation, solvent_data, entries, reactions, model, temperatures):
    results = []
    for temperature in temperatures:
        try:
            values = {}
            for key, entry in entries.items():
                if model == "linear_HS":
                    correction = solvation.get_solvation_correction(entry["solute"], solvent_data)
                    values[key] = correction.enthalpy - temperature * correction.entropy
                else:
                    values[key] = solvation.get_T_dep_solvation_energy_from_LSER_298(entry["solute"], solvent_data, temperature)[0]
            rows = {}
            for name, reaction in reactions.items():
                delta = reaction_delta(reaction, values)
                ratio = math.exp(-delta / (constants.R * temperature))
                gas_kc = reaction.get_equilibrium_constant(temperature, type="Kc")
                rows[name] = {"ddG_J_mol": delta, "ratio": ratio, "gas_Kc": gas_kc, "melt_Kc": gas_kc * ratio}
            results.append({"T_K": temperature, "species_solvation_J_mol": values, "reactions": rows})
        except Exception as error:
            results.append({"T_K": temperature, "error_type": type(error).__name__, "error": str(error)})
    return results


def crossing(temperatures, residuals):
    for index in range(len(temperatures) - 1):
        if residuals[index] == 0.0:
            return temperatures[index]
        if residuals[index] * residuals[index + 1] < 0:
            fraction = -residuals[index] / (residuals[index + 1] - residuals[index])
            return temperatures[index] + fraction * (temperatures[index + 1] - temperatures[index])
    return None


def ceiling_results(solvation, solvent_data, entries, propagation, rates, model, temperatures):
    data = correction_table(solvation, solvent_data, entries, {"propagation": propagation}, model, temperatures)
    if any("error" in row for row in data):
        return {"status": "unavailable on full ceiling grid", "errors": [row for row in data if "error" in row]}
    shifts = [row["reactions"]["propagation"]["ddG_J_mol"] for row in data]
    residuals = [math.log(prop * MONOMER_CONCENTRATION / dep) - shift / (constants.R * temperature)
                 for prop, dep, shift, temperature in zip(rates["propagation_grid"]["k"],
                                                        rates["depropagation_grid"]["k"], shifts, temperatures)]
    table_root = crossing(temperatures, residuals)
    def residual(temperature):
        row = correction_table(solvation, solvent_data, entries, {"propagation": propagation}, model, [temperature])[0]
        return math.log(row["reactions"]["propagation"]["melt_Kc"] * MONOMER_CONCENTRATION)
    exact_root = brentq(residual, 600.0, 800.0) if residual(600.0) * residual(800.0) < 0 else None
    return {"status": "computed", "tabulated_K": table_root, "continuous_K": exact_root,
            "temperatures_K": list(temperatures), "rate_residuals": residuals,
            "continuous_root_residual": residual(exact_root) if exact_root else None}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", type=Path, default=Path("/home/alon/Code/RMG-database"))
    parser.add_argument("--scratch", type=Path, default=Path("/tmp/i034-reproduce/database"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    repository = Path(__file__).resolve().parents[5]
    subprocess.check_call(["git", "-C", str(repository), "merge-base", "--is-ancestor", BASE_SHA, "HEAD"])
    args.scratch = args.scratch.resolve()
    args.output = args.output.resolve()
    if args.scratch.is_relative_to(args.database.resolve()) or args.output.is_relative_to(args.database.resolve()):
        raise ValueError("scratch/output must not be inside the database")
    progress("materializing only pinned thermo, solvation and five-family files with git show")
    snapshot = snapshot_database(args.database, args.scratch)
    progress("loading pinned database; no rate-tree regeneration")
    database = RMGDatabase()
    database.load_kinetics(str(args.scratch / "input/kinetics"), reaction_libraries=[], seed_mechanisms=None,
                           kinetics_families=list(PS_FAMILY_CANDIDATES), kinetics_depositories=["training"])
    database.load_thermo(str(args.scratch / "input/thermo"), thermo_libraries=["primaryThermoLibrary"], depository=True)
    database.load_solvation(str(args.scratch / "input/solvation"))
    reactions, rates = gas_baseline(database, ps_proxy_set(PS_PROXY_UNITS))
    entries, failures = estimate_species(database, reactions)
    if failures:
        raise AssertionError(f"solute failures: {failures}")
    progress("all requested gas baselines reproduced before any solvation correction")
    gas_rows = []
    for temperature in TEMPERATURES:
        gas_rows.append({"T_K": temperature, "reactions": {
            name: {"Kc": reaction.get_equilibrium_constant(temperature, type="Kc"),
                   "dG_J_mol": reaction.get_free_energy_of_reaction(temperature)}
            for name, reaction in reactions.items()}})
    prop = reactions["propagation"]
    rates["ceiling_continuous_K"] = brentq(lambda temperature: math.log(
        prop.get_equilibrium_constant(temperature, type="Kc") * MONOMER_CONCENTRATION), 600.0, 800.0)
    inventory_result = inventory(database.solvation, args.scratch)
    candidate_results = {}
    for label in CANDIDATES:
        progress(f"solvation stand-in: {label}")
        solvent = database.solvation.get_solvent_data(label)
        missing = [field for field in FIELDS_G + FIELDS_H if getattr(solvent, field) is None]
        candidate = {"coefficients": {field: getattr(solvent, field) for field in FIELDS_G + FIELDS_H},
                     "missing": missing, "coolprop": solvent.name_in_coolprop,
                     "Tc_K": get_critical_temperature(solvent.name_in_coolprop) if solvent.name_in_coolprop else None,
                     "boundaries": collect_boundaries(database.solvation, solvent, next(iter(entries.values()))["solute"]),
                     "G298_reactions_J_mol": {name: reaction_delta(reaction, {
                         key: database.solvation.calc_g(entry["solute"], solvent) for key, entry in entries.items()})
                         for name, reaction in reactions.items()}, "models": {}}
        if not missing:
            candidate["anchors"] = {key: {
                "G298_J_mol": correction.gibbs, "H298_J_mol": correction.enthalpy, "S298_J_mol_K": correction.entropy}
                for key, entry in entries.items()
                for correction in [database.solvation.get_solvation_correction(entry["solute"], solvent)]}
            candidate["models"]["linear_HS"] = correction_table(database.solvation, solvent, entries, reactions,
                                                                  "linear_HS", TEMPERATURES)
            candidate["linear_ceiling"] = ceiling_results(database.solvation, solvent, entries, prop, rates,
                                                            "linear_HS", DEFAULT_T_GRID)
        candidate["models"]["Kfactor"] = correction_table(database.solvation, solvent, entries, reactions,
                                                           "Kfactor", TEMPERATURES)
        candidate_results[label] = candidate
    result = {
        "provenance": {"repository_base_sha": BASE_SHA, "database_sha": DATABASE_SHA,
                       "prior_probe_sha": PRIOR_SHA, "prior_script_sha256": hashlib.sha256(git_output(
                           repository, "show", f"{PRIOR_SHA}:test/rmgpy/kmc/fixtures/i032_probe/run_probe.py")).hexdigest(),
                       "python": platform.python_version(), "CoolProp": CoolProp.__version__,
                       "R_J_mol_K": constants.R, "snapshot": snapshot,
                       "code_sha256": {relative: hashlib.sha256((repository / relative).read_bytes()).hexdigest()
                                       for relative in ("rmgpy/data/solvation.py", "rmgpy/thermo/thermoengine.py", "rmgpy/kmc/compiler.py")}},
        "temperatures_K": list(TEMPERATURES), "monomer_concentration_mol_m3": MONOMER_CONCENTRATION,
        "boundary_reference_smiles": next(iter(entries)),
        "rates": rates, "gas": gas_rows, "inventory": inventory_result, "candidates": candidate_results,
        "hexadecane_literature_Tc": HEXADECANE_TC,
        "species": {key: {field: value for field, value in entry.items() if field != "solute"}
                    for key, entry in sorted(entries.items())},
        "solute_failures": failures,
        "radical_audit": audit_radicals(database.solvation, entries),
        "reactions": {name: {"reactants": smiles(reaction.reactants), "products": smiles(reaction.products),
                             "delta_n": len(reaction.products) - len(reaction.reactants)}
                      for name, reaction in reactions.items()},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n")
    print(f"I-034 gas baselines reproduced; {len(entries)} selected solutes, {len(failures)} failures")
    print(f"I-034 wrote {args.output}")


if __name__ == "__main__":
    main()
