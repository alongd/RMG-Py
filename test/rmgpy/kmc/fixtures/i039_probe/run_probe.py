"""Diagnose PS equilibrium thermo; all mutations are to scratch or probe output.

Run from the repository root with the Python and logging command in the report.
Pinned database materialization and baseline selection reuse the committed I034
probe, so this probe follows precisely the same library and reaction choices.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys

from scipy.optimize import brentq

from rmgpy import constants
from rmgpy.data.rmg import RMGDatabase
from rmgpy.kmc.compiler import DEFAULT_T_GRID, EventSetCompiler, ps_proxy_set
from rmgpy.molecule.symmetry import calculate_atom_symmetry_number
from rmgpy.thermo import ThermoData

FIXTURES = Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location("i034_baseline", FIXTURES / "i034_probe/run_probe.py")
prior = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(prior)


T0 = 298.15
TEMPERATURES = (T0, 600.0, 668.0, 700.0, 800.0)
P0 = 100000.0
C0 = 1000.0
# Literature inputs, not fitted parameters. Full citations and qualifications
# are in the report. Odian's table does not specify a reference temperature;
# using T0 below is explicitly a CONDITIONAL accounting exercise.
LITERATURE = {
    "contract_target_K": 668.0,
    "Odian": {"H_J_mol": -73000.0, "S_J_mol_K": -104.0,
              "H_rounding_J_mol": 500.0, "S_rounding_J_mol_K": 0.5,
              "monomer_standard": "1 mol/L; condensed polymer",
              "reference_temperature_K": None,
              "source": "Odian, Principles of Polymerization, 4th ed., Table 3-15"},
    "Roberts_solid": {"T_K": 298.15, "H_J_mol": -69790.0,
                      "published_uncertainty_J_mol": 660.0},
    "Roberts_solution": {"T_K": 298.15, "H_J_mol": -73380.0,
                         "published_uncertainty_J_mol": 690.0,
                         "polymer_weight_percent": 6.9},
    "Warfield": {"T_K": 298.16, "entropy_loss_cal_mol_K": 26.69,
                 "S_J_mol_K": -26.69 * 4.184,
                 "uncertainty": None},
    "Cowie": {"pure_liquid_Tc_K": 583.0, "H_J_mol": -68500.0},
    "Kirk_Othmer": {"T_K": 298.15, "P_Pa": 101300.0,
                    "liquid_H_J_mol": -69900.0, "liquid_S_J_mol_K": -104.6,
                    "liquid_Tc_C": 395.0,
                    "gas_H_J_mol": -113400.0, "gas_S_J_mol_K": -212.1,
                    "gas_Tc_C": 262.0,
                    "source": "Kirk-Othmer, Xylylene Polymers, Table 2, p. 420; liquid/gas monomer reference"},
}


def progress(message):
    print("[I039] " + message, file=sys.stderr, flush=True)


def close(a, b, absolute=2e-6):
    if not math.isclose(a, b, rel_tol=2e-10, abs_tol=absolute):
        raise AssertionError(f"decomposition mismatch: {a} != {b}")


def hs(data, temperature):
    return [float(data.get_enthalpy(temperature)), float(data.get_entropy(temperature))]


def resolved_data(database, entry):
    while isinstance(entry.data, str):
        entry = database.entries[entry.data]
    return entry.data


def decompose(thermo_database, species):
    source = thermo_database.extract_source_from_comments(species)
    contributions = []
    if "Library" in source:
        library = thermo_database.libraries[source["Library"]]
        mol = species.molecule[0].copy(deep=True)
        mol.saturate_radicals()
        saturated = prior.Species(molecule=[mol])
        hit = thermo_database.get_thermo_data_from_library(saturated, library)
        if hit is None:
            raise AssertionError("source library does not contain saturated structure")
        contributions.append(("library:" + source["Library"] + ":" + hit[2].label, 1, hit[0]))
    for kind, entries in source.get("GAV", {}).items():
        for entry, weight in entries:
            if weight:
                contributions.append((kind + ":" + entry.label, weight,
                                      resolved_data(thermo_database.groups[kind], entry)))
    radical_count = species.molecule[0].get_radical_count()
    if radical_count and "GAV" in source:
        contributions.append(("HBI:subtract_H_atoms", radical_count,
                              ThermoData(Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
                                         Cpdata=([0.0] * 7, "J/(mol*K)"),
                                         H298=(-52.103 * 4184, "J/mol"), S298=(0.0, "J/(mol*K)"))))
    sigma = float(species.get_symmetry_number())
    hybrid = species.get_resonance_hybrid()
    optical_sites = sum(
        calculate_atom_symmetry_number(hybrid, atom) == 0.5
        for atom in hybrid.atoms if not hybrid.is_atom_in_cycle(atom)
    )
    symmetry_applied = "Library" not in source or radical_count > 0
    # Library HBI removes saturated symmetry internally; expose it separately.
    if "Library" in source and radical_count:
        saturated_sigma = float(saturated.get_symmetry_number())
        contributions.append(("HBI:undo_saturated_library_symmetry", 1,
                              ThermoData(Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
                                         Cpdata=([0.0] * 7, "J/(mol*K)"), H298=(0.0, "J/mol"),
                                         S298=(constants.R * math.log(saturated_sigma), "J/(mol*K)"))))
    rows = []
    for temperature in TEMPERATURES:
        terms = {label: [weight * value for value in hs(data, temperature)]
                 for label, weight, data in contributions}
        if symmetry_applied:
            terms["symmetry:without_optical_factor"] = [0.0, -constants.R * math.log(sigma * 2**optical_sites)]
            terms["optical:atom_half_factors"] = [0.0, constants.R * optical_sites * math.log(2)]
        actual = hs(species.thermo, temperature)
        total = [sum(term[i] for term in terms.values()) for i in range(2)]
        for a, b in zip(actual, total):
            close(a, b)
        rows.append({"T_K": temperature, "actual_HS": actual, "terms": terms})
    return {"smiles": species.molecule[0].to_smiles(),
            "adjacency": species.molecule[0].to_adjacency_list(),
            "thermo_comment": species.thermo.comment,
            "symmetry_number_including_optical": sigma,
            "optical_atom_half_factors": optical_sites,
            "symmetry_separately_applied": symmetry_applied,
            "source_weights": {label: weight for label, weight, data in contributions},
            "rows": rows}


def compile_pair(database, units):
    proxies = tuple(p for p in ps_proxy_set(units)
                    if p.site_type in ("end_radical", "end_radical+styrene"))
    active, excluded, cache = EventSetCompiler.discover_family_reactions(
        database.kinetics, proxies, ("R_Addition_MultipleBond",),
        family_universe=["R_Addition_MultipleBond"])
    artifact = EventSetCompiler(database.kinetics, proxies, active, excluded_families=excluded,
                               thermo_database=database.thermo, reaction_cache=cache,
                               rmg_database_sha=prior.DATABASE_SHA).compile()
    pair = artifact["ps_ceiling_pairs"][0]
    prop = next(r for r in artifact["records"] if r["event_id"] == pair["propagation_event_id"])
    dep = next(r for r in artifact["records"] if r["event_id"] == pair["depropagation_event_id"])
    return prior.reaction_from_record(prop), {
        "ceiling_tabulated_K": artifact["ps_ceiling_temperature_K"],
        "ceiling_pair_count": len(artifact["ps_ceiling_pairs"]),
        "propagation_grid": prop["k_table"], "depropagation_grid": dep["k_table"]}


def pressure_hs(reaction, temperature):
    return [float(reaction.get_enthalpy_of_reaction(temperature)),
            float(reaction.get_entropy_of_reaction(temperature))]


def concentration_hs(reaction, temperature):
    h, s = pressure_hs(reaction, temperature)
    nu = len(reaction.products) - len(reaction.reactants)
    alpha = C0 * constants.R * temperature / P0
    # Derivatives of Gc = Gp + nu RT ln(c0 RT / p0), at fixed c0.
    return [h - nu * constants.R * temperature,
            s - nu * constants.R * (math.log(alpha) + 1.0)]


def chain_end_audit(database):
    """Check a benzylic terminal radical, selected by structure, not by Tc."""
    molecules = ["CC(c1ccccc1)C[CH](c1ccccc1)", "C=Cc1ccccc1",
                 "CC(c1ccccc1)CC(c1ccccc1)C[CH](c1ccccc1)"]
    species = [prior.Species(molecule=[prior.Molecule(smiles=s)]) for s in molecules]
    for item in species:
        item.molecule[0].update()
        item.generate_resonance_structures()
        item.thermo = database.thermo.get_thermo_data(item)
    reaction = prior.Reaction(reactants=species[:2], products=species[2:])
    return {"reactants": prior.smiles(reaction.reactants), "products": prior.smiles(reaction.products),
            "species": [decompose(database.thermo, item) for item in species],
            "continuous_K": root(lambda t: math.log(reaction.get_equilibrium_constant(t) * C0)),
            "pressure_HS_at_anchor": pressure_hs(reaction, T0),
            "purpose": "terminal benzylic propagation; structure-only control, not a correction selection"}


def root(function):
    return float(brentq(function, 250.0, 1800.0, xtol=1e-9))


def accounting(reaction, gas3, gas5, solvated, lit_h, lit_s):
    """Symmetric H/S counterfactual allocation; not a causal error estimate."""
    base_h, base_s = concentration_hs(reaction, T0)
    def free_energy(t, model_h, model_s):
        h, s = concentration_hs(reaction, t)
        if not model_h:
            h += lit_h - base_h
        if not model_s:
            s += lit_s - base_s
        return h - t * s
    roots = {f"{int(h)}{int(s)}": root(lambda t: free_energy(t, h, s))
             for h in (False, True) for s in (False, True)}
    h_effect = ((roots["10"] - roots["00"]) + (roots["11"] - roots["01"])) / 2
    s_effect = ((roots["01"] - roots["00"]) + (roots["11"] - roots["10"])) / 2
    rows = {
        "convention": 0.0,
        "enthalpy_error_conditional": h_effect,
        "entropy_error_conditional": s_effect,
        "proxy_size_observed": gas3 - gas5,
        "solvation_including_grid_definition": solvated - gas3,
        "residual_physics_not_represented_or_unresolved_reference": roots["00"] - LITERATURE["contract_target_K"],
    }
    close(sum(rows.values()), solvated - LITERATURE["contract_target_K"])
    return {"roots_K": roots, "rows_K": rows, "sum_K": sum(rows.values())}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", type=Path, default=Path("/home/alon/Code/RMG-database"))
    parser.add_argument("--scratch", type=Path, default=Path("/tmp/i039-reproduce/database"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    repository = FIXTURES.parents[3]
    for target in (args.scratch, args.output):
        target = target.resolve()
        if target.is_relative_to(args.database.resolve()) or target.is_relative_to(repository):
            raise ValueError("scratch and output must be outside repositories")
        if "catalog" in target.parts or target.is_relative_to(Path("/home/alon/Code/polymers")):
            raise ValueError("excluded path")
    progress("materializing pinned database with git show; no tree regeneration")
    snapshot = prior.snapshot_database(args.database, args.scratch)
    database = RMGDatabase()
    database.load_kinetics(str(args.scratch / "input/kinetics"), reaction_libraries=[],
                           seed_mechanisms=None, kinetics_families=list(prior.PS_FAMILY_CANDIDATES),
                           kinetics_depositories=["training"])
    database.load_thermo(str(args.scratch / "input/thermo"),
                         thermo_libraries=["primaryThermoLibrary"], depository=True)
    database.load_solvation(str(args.scratch / "input/solvation"))
    progress("reproducing I032 and I034 propagation selection and rates first")
    baseline_reactions, baseline_rates = prior.gas_baseline(database, ps_proxy_set(3))
    entries, failures = prior.estimate_species(database, {"propagation": baseline_reactions["propagation"]})
    assert not failures
    prop3 = baseline_reactions["propagation"]
    solvent = database.solvation.get_solvent_data("toluene")
    solvent_ceiling = prior.ceiling_results(database.solvation, solvent, entries, prop3,
                                           baseline_rates, "linear_HS", DEFAULT_T_GRID)
    close(solvent_ceiling["continuous_K"], 781.324909, absolute=1e-6)
    close(solvent_ceiling["tabulated_K"], 781.491601, absolute=1e-6)
    results = {}
    reactions = {}
    for units in (3, 4, 5):
        progress(f"thermo decomposition and compiled ceiling for L={units}")
        reaction, rates = (prop3, baseline_rates) if units == 3 else compile_pair(database, units)
        for species in reaction.reactants + reaction.products:
            if species.thermo is None:
                species.thermo = database.thermo.get_thermo_data(species)
        continuous = root(lambda t: math.log(reaction.get_equilibrium_constant(t, type="Kc") * C0))
        convention_root = root(lambda t: concentration_hs(reaction, t)[0]
                               - t * concentration_hs(reaction, t)[1])
        close(continuous, convention_root)
        sources = [decompose(database.thermo, species)
                   for species in reaction.reactants + reaction.products]
        signs = [-1] * len(reaction.reactants) + [1] * len(reaction.products)
        rows = []
        for i, temperature in enumerate(TEMPERATURES):
            net = {}
            for sign, species_source in zip(signs, sources):
                for label, values in species_source["rows"][i]["terms"].items():
                    net.setdefault(label, [0.0, 0.0])
                    net[label] = [a + sign * b for a, b in zip(net[label], values)]
            rows.append({"T_K": temperature, "pressure_HS": pressure_hs(reaction, temperature),
                         "concentration_HS": concentration_hs(reaction, temperature),
                         "net_source_terms": net,
                         "delta_Cp_pressure_J_mol_K": float(sum(s.thermo.get_heat_capacity(temperature)
                                                                 for s in reaction.products)
                                                           - sum(s.thermo.get_heat_capacity(temperature)
                                                                 for s in reaction.reactants))})
        reactions[units] = reaction
        results[str(units)] = {"rates": rates, "continuous_K": continuous,
                               "same_state_1M_K": convention_root,
                               "one_bar_monomer_root_K": root(lambda t: reaction.get_free_energy_of_reaction(t)),
                               "species": sources, "reaction_rows": rows,
                               "reactants": prior.smiles(reaction.reactants),
                               "products": prior.smiles(reaction.products), "signs": signs}
    progress("decomposing toluene LSER terms and radical lookup failures")
    fields = (("S", "s"), ("B", "b"), ("E", "e"), ("L", "l"), ("A", "a"))
    descriptors = {name: prior.reaction_delta(prop3, {key: e["descriptors"][name]
                                                   for key, e in entries.items()})
                   for name, _ in fields}
    nu = len(prop3.products) - len(prop3.reactants)
    terms = {}
    for name, prefix in (*fields, ("intercept", "c")):
        weight = descriptors[name] if name != "intercept" else nu
        dh = 1000 * weight * getattr(solvent, prefix + "_h")
        dg = -8.314 * 298 * 2.303 * weight * getattr(solvent, prefix + "_g")
        terms[name] = {"descriptor_delta": weight, "ddH_J_mol": dh,
                       "ddG298_J_mol": dg, "ddS_J_mol_K": (dh - dg) / 298}
    anchors = {key: {"H": c.enthalpy, "S": c.entropy, "G": c.gibbs}
               for key, e in entries.items()
               for c in [database.solvation.get_solvation_correction(e["solute"], solvent)]}
    dh = prior.reaction_delta(prop3, {key: a["H"] for key, a in anchors.items()})
    ds = prior.reaction_delta(prop3, {key: a["S"] for key, a in anchors.items()})
    close(sum(v["ddH_J_mol"] for v in terms.values()), dh)
    close(sum(v["ddS_J_mol_K"] for v in terms.values()), ds)
    descriptor_sources = {}
    for sign, species in [(-1, s) for s in prop3.reactants] + [(1, s) for s in prop3.products]:
        entry = entries[species.molecule[0].to_smiles()]
        comment = entry["comment"].replace(" + ", " +").replace(" - ", " -")
        for token in comment.split():
            weight = -1 if token.startswith("-") else 1
            token = token.lstrip("+-")
            for kind, groups in database.solvation.groups.items():
                if token.startswith(kind + "(") and token.endswith(")"):
                    label = token[len(kind) + 1:-1]
                    key = kind + ":" + label
                    data = resolved_data(groups, groups.entries[label])
                    contribution = descriptor_sources.setdefault(key, {"weight": 0, "descriptors": {f: 0.0 for f, _ in fields}})
                    contribution["weight"] += sign * weight
                    for f, _ in fields:
                        contribution["descriptors"][f] += sign * weight * getattr(data, f)
                    break
    for f, _ in fields:
        close(sum(row["descriptors"][f] for row in descriptor_sources.values()), descriptors[f])
    for row in descriptor_sources.values():
        row["ddH_J_mol"] = 1000 * sum(row["descriptors"][f] * getattr(solvent, p + "_h") for f, p in fields)
        row["ddG298_J_mol"] = -8.314 * 298 * 2.303 * sum(row["descriptors"][f] * getattr(solvent, p + "_g") for f, p in fields)
        row["ddS_J_mol_K"] = (row["ddH_J_mol"] - row["ddG298_J_mol"]) / 298
    structural_control = chain_end_audit(database)
    solvent_controls = {}
    for label, control_h, control_s in (
            ("enthalpy_only", dh, 0.0), ("entropy_only", 0.0, ds),
            ("both", dh, ds),
            ("intercept_only", terms["intercept"]["ddH_J_mol"], terms["intercept"]["ddS_J_mol_K"]),
            ("descriptor_terms_only", dh - terms["intercept"]["ddH_J_mol"], ds - terms["intercept"]["ddS_J_mol_K"])):
        solvent_controls[label] = root(lambda t: math.log(prop3.get_equilibrium_constant(t) * C0)
                                      - (control_h - t * control_s) / (constants.R * t))
    optical_delta = results["3"]["reaction_rows"][0]["net_source_terms"]["optical:atom_half_factors"][1]
    optical_control = root(lambda t: math.log(prop3.get_equilibrium_constant(t) * C0) - optical_delta / constants.R)
    # The compilation's gas MONOMER -> condensed polymer standard, keeping its
    # H/S fixed, normalized to our gas pressure and concentration conventions.
    compiled_lit = LITERATURE["Kirk_Othmer"]
    lit_gas_s_bar = compiled_lit["gas_S_J_mol_K"] + constants.R * math.log(P0 / compiled_lit["P_Pa"])
    lit_gas_h = compiled_lit["gas_H_J_mol"]
    comparison = {
        "liquid_constant_HS_K": compiled_lit["liquid_H_J_mol"] / compiled_lit["liquid_S_J_mol_K"],
        "gas_constant_HS_at_source_pressure_K": compiled_lit["gas_H_J_mol"] / compiled_lit["gas_S_J_mol_K"],
        "gas_constant_HS_at_1bar_K": lit_gas_h / lit_gas_s_bar,
        "gas_constant_HS_at_1M_K": root(lambda t: lit_gas_h - t * lit_gas_s_bar
                                      - constants.R * t * math.log(C0 * constants.R * t / P0)),
        "gas_entropy_1bar_J_mol_K": lit_gas_s_bar,
        "full_common_liquid_convention": "unidentified: requires monomer chemical potential/activity and condensed-polymer increments",
    }
    central = accounting(reactions[5], results["3"]["continuous_K"], results["5"]["continuous_K"],
                         solvent_ceiling["tabulated_K"], -73000, -104)
    sensitivity = [accounting(reactions[5], results["3"]["continuous_K"], results["5"]["continuous_K"],
                              solvent_ceiling["tabulated_K"], h, s)
                   for h in (-73500, -72500) for s in (-104.5, -103.5)]
    for name, value in central["rows_K"].items():
        central.setdefault("rounding_only_half_ranges_K", {})[name] = max(
            abs(case["rows_K"][name] - value) for case in sensitivity)
    output = {
        "provenance": {"database_sha": prior.DATABASE_SHA, "snapshot": snapshot,
                       "repository_sha_before_probe_commit": subprocess.check_output(
                           ["git", "-C", str(repository), "rev-parse", "HEAD"], text=True).strip(),
                       "source_sha256": {path: hashlib.sha256((repository / path).read_bytes()).hexdigest()
                                         for path in ("rmgpy/kmc/compiler.py", "rmgpy/data/thermo.py",
                                                      "rmgpy/data/solvation.py", "rmgpy/reaction.py",
                                                      "rmgpy/molecule/symmetry.py")}},
        "constants": {"R": constants.R, "P0_Pa": P0, "C0_mol_m3": C0, "T0_K": T0},
        "literature": LITERATURE, "baselines": {"gas_grid_K": baseline_rates["ceiling_tabulated_K"],
                                                  "toluene": solvent_ceiling},
        "sizes": results,
        "chain_end_control": structural_control,
        "literature_convention_comparison": comparison,
        "optical_term_suppressed_control_K": optical_control,
        "solvation": {"terms": terms, "ddH_J_mol": dh, "ddS_J_mol_K": ds,
                      "descriptor_source_terms": descriptor_sources,
                      "continuous_controls_K": solvent_controls,
                      "anchors": anchors, "critical_temperature_K": prior.get_critical_temperature(solvent.name_in_coolprop),
                      "species": {key: {k: v for k, v in e.items() if k != "solute"}
                                  for key, e in entries.items()},
                      "radical_audit": prior.audit_radicals(database.solvation, entries)},
        "conditional_accounting": central,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(f"I039 reproduced propagation baselines; decomposed L=3,4,5; wrote {args.output}")


if __name__ == "__main__":
    main()
