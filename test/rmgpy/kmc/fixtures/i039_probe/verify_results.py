"""Re-estimate species and independently reproduce the probe's numerical claims.

Run with results JSON; --snapshot defaults to /tmp/i039-reproduce/database.
This script intentionally does not import the generating probe's calculations.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from scipy.optimize import brentq

from rmgpy import constants
from rmgpy.data.rmg import RMGDatabase
from rmgpy.molecule.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species

from render_tables import REPORT, blocks


def close(actual, expected, absolute=3e-6, relative=1e-9):
    if not math.isclose(actual, expected, abs_tol=absolute, rel_tol=relative):
        raise AssertionError(f"{actual} != {expected}")


def crossing(temperatures, residuals):
    for t1, t2, f1, f2 in zip(temperatures, temperatures[1:], residuals, residuals[1:]):
        if f1 == 0:
            return t1
        if f1 * f2 < 0:
            return t1 - f1 * (t2 - t1) / (f2 - f1)
    raise AssertionError("unbracketed grid crossing")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("--snapshot", type=Path, default=Path("/tmp/i039-reproduce/database"))
    parser.add_argument("--report", type=Path, default=REPORT)
    args = parser.parse_args()
    result = json.loads(args.results.read_text())
    root = Path(__file__).resolve().parents[5]
    assert result["provenance"]["database_sha"] == "4a12d36fcdc193ede82c8d1ab5c1653495d445bc"
    digest = hashlib.sha256()
    paths = sorted(path for prefix in result["provenance"]["snapshot"]["prefixes"]
                   for path in ([args.snapshot / prefix] if (args.snapshot / prefix).is_file()
                                else (args.snapshot / prefix).rglob("*")) if path.is_file())
    for path in paths:
        relative = path.relative_to(args.snapshot).as_posix()
        assert "catalog" not in path.parts
        digest.update(relative.encode() + b"\0" + path.read_bytes())
    assert digest.hexdigest() == result["provenance"]["snapshot"]["sha256"]
    for path, expected in result["provenance"]["source_sha256"].items():
        assert hashlib.sha256((root / path).read_bytes()).hexdigest() == expected
    db = RMGDatabase()
    db.load_thermo(str(args.snapshot / "input/thermo"), thermo_libraries=["primaryThermoLibrary"], depository=True)
    db.load_solvation(str(args.snapshot / "input/solvation"))
    reactions = {}
    for units, probe in result["sizes"].items():
        participants = []
        for source in probe["species"]:
            species = Species(molecule=[Molecule().from_adjacency_list(source["adjacency"])])
            species.thermo = db.thermo.get_thermo_data(species)
            close(species.get_symmetry_number(), source["symmetry_number_including_optical"])
            for row in source["rows"]:
                t = row["T_K"]
                actual = [species.thermo.get_enthalpy(t), species.thermo.get_entropy(t)]
                for i, value in enumerate(actual):
                    close(value, row["actual_HS"][i])
                    close(value, sum(term[i] for term in row["terms"].values()))
            participants.append(species)
        reaction = Reaction(reactants=participants[:2], products=participants[2:])
        assert len(reaction.reactants) == 2 and len(reaction.products) == 1
        assert reaction.reactants[1].molecule[0].get_element_count()["C"] == 8
        assert reaction.products[0].molecule[0].get_element_count()["C"] == 8 * int(units)
        for row in probe["reaction_rows"]:
            t = row["T_K"]
            h = sum(s.thermo.get_enthalpy(t) for s in reaction.products) - sum(s.thermo.get_enthalpy(t) for s in reaction.reactants)
            s = sum(s.thermo.get_entropy(t) for s in reaction.products) - sum(s.thermo.get_entropy(t) for s in reaction.reactants)
            close(h, row["pressure_HS"][0])
            close(s, row["pressure_HS"][1])
            close(h + constants.R * t, row["concentration_HS"][0])
            close(s + constants.R * (1 + math.log(1000 * constants.R * t / 1e5)), row["concentration_HS"][1])
            close(h, sum(term[0] for term in row["net_source_terms"].values()))
            close(s, sum(term[1] for term in row["net_source_terms"].values()))
            kp = math.exp(-(h - t * s) / (constants.R * t))
            kc = kp * constants.R * t / 1e5
            close(kc, reaction.get_equilibrium_constant(t, type="Kc"))
        exact = brentq(lambda t: reaction.get_free_energy_of_reaction(t)
                       - constants.R * t * math.log(1000 * constants.R * t / 1e5), 250, 1800)
        close(exact, probe["continuous_K"])
        close(exact, probe["same_state_1M_K"])
        p_root = brentq(lambda t: reaction.get_free_energy_of_reaction(t), 250, 1800)
        close(p_root, probe["one_bar_monomer_root_K"])
        forward, reverse = probe["rates"]["propagation_grid"], probe["rates"]["depropagation_grid"]
        assert forward["T"] == reverse["T"]
        residuals = [math.log(kp * 1000 / kd) for kp, kd in zip(forward["k"], reverse["k"])]
        for t, kp, kd in zip(forward["T"], forward["k"], reverse["k"]):
            close(kp / kd, reaction.get_equilibrium_constant(t), absolute=1e-12)
        close(crossing(forward["T"], residuals), probe["rates"]["ceiling_tabulated_K"])
        reactions[units] = reaction
    close(result["baselines"]["gas_grid_K"], 710.248673, absolute=1e-6)
    close(result["sizes"]["3"]["continuous_K"], 710.020462, absolute=1e-6)
    # Independent LSER algebra, including the species-count intercept.
    solvent = db.solvation.get_solvent_data("toluene")
    net_h, net_s = 0.0, 0.0
    prop = reactions["3"]
    for sign, species in [( -1, s) for s in prop.reactants] + [(1, s) for s in prop.products]:
        descriptors = db.solvation.get_solute_data(species)
        dh = 1000 * (solvent.c_h + sum(getattr(descriptors, f) * getattr(solvent, p + "_h")
                                     for f, p in (("S", "s"), ("B", "b"), ("E", "e"), ("L", "l"), ("A", "a"))))
        dg = -8.314 * 298 * 2.303 * (solvent.c_g + sum(getattr(descriptors, f) * getattr(solvent, p + "_g")
                                                       for f, p in (("S", "s"), ("B", "b"), ("E", "e"), ("L", "l"), ("A", "a"))))
        net_h += sign * dh
        net_s += sign * (dh - dg) / 298
    close(net_h, result["solvation"]["ddH_J_mol"])
    close(net_s, result["solvation"]["ddS_J_mol_K"])
    for name, term in result["solvation"]["terms"].items():
        prefix = {"S": "s", "B": "b", "E": "e", "L": "l", "A": "a", "intercept": "c"}[name]
        weight = term["descriptor_delta"]
        close(term["ddH_J_mol"], 1000 * weight * getattr(solvent, prefix + "_h"))
        close(term["ddG298_J_mol"], -8.314 * 298 * 2.303 * weight * getattr(solvent, prefix + "_g"))
        close(term["ddS_J_mol_K"], (term["ddH_J_mol"] - term["ddG298_J_mol"]) / 298)
    for name, term in result["solvation"]["descriptor_source_terms"].items():
        kind, label = name.split(":", 1)
        entry = db.solvation.groups[kind].entries[label]
        while isinstance(entry.data, str):
            entry = db.solvation.groups[kind].entries[entry.data]
        for field, value in term["descriptors"].items():
            close(value, term["weight"] * getattr(entry.data, field))
        h = 1000 * sum(value * getattr(solvent, field.lower() + "_h") for field, value in term["descriptors"].items())
        g = -8.314 * 298 * 2.303 * sum(value * getattr(solvent, field.lower() + "_g") for field, value in term["descriptors"].items())
        close(h, term["ddH_J_mol"])
        close(g, term["ddG298_J_mol"])
        close((h-g)/298, term["ddS_J_mol_K"])
    corrected_root = brentq(lambda t: math.log(prop.get_equilibrium_constant(t) * 1000)
                            - (net_h - t * net_s) / (constants.R * t), 600, 800)
    close(corrected_root, result["baselines"]["toluene"]["continuous_K"])
    close(corrected_root, 781.324909, absolute=1e-6)
    forward = result["sizes"]["3"]["rates"]["propagation_grid"]
    residuals = [math.log(prop.get_equilibrium_constant(t) * 1000)
                 - (net_h - t * net_s) / (constants.R * t) for t in forward["T"]]
    corrected_grid = crossing(forward["T"], residuals)
    close(corrected_grid, result["baselines"]["toluene"]["tabulated_K"])
    intercept = result["solvation"]["terms"]["intercept"]
    for name, h, s in (("enthalpy_only", net_h, 0), ("entropy_only", 0, net_s),
                       ("both", net_h, net_s),
                       ("intercept_only", intercept["ddH_J_mol"], intercept["ddS_J_mol_K"]),
                       ("descriptor_terms_only", net_h-intercept["ddH_J_mol"], net_s-intercept["ddS_J_mol_K"])):
        t = brentq(lambda temp: math.log(prop.get_equilibrium_constant(temp) * 1000)
                   - (h-temp*s)/(constants.R*temp), 250, 1800)
        close(t, result["solvation"]["continuous_controls_K"][name])
    optical_entropy = result["sizes"]["3"]["reaction_rows"][0]["net_source_terms"]["optical:atom_half_factors"][1]
    t = brentq(lambda temp: math.log(prop.get_equilibrium_constant(temp) * 1000)
               - optical_entropy / constants.R, 250, 1800)
    close(t, result["optical_term_suppressed_control_K"])
    control_species = []
    for source in result["chain_end_control"]["species"]:
        species = Species(molecule=[Molecule().from_adjacency_list(source["adjacency"])])
        species.generate_resonance_structures()
        species.thermo = db.thermo.get_thermo_data(species)
        for row in source["rows"]:
            close(species.thermo.get_enthalpy(row["T_K"]), row["actual_HS"][0])
            close(species.thermo.get_entropy(row["T_K"]), row["actual_HS"][1])
        control_species.append(species)
    control_reaction = Reaction(reactants=control_species[:2], products=control_species[2:])
    close(brentq(lambda temp: math.log(control_reaction.get_equilibrium_constant(temp) * 1000), 250, 1800),
          result["chain_end_control"]["continuous_K"])
    compiled_lit = result["literature"]["Kirk_Othmer"]
    h, s = compiled_lit["gas_H_J_mol"], compiled_lit["gas_S_J_mol_K"]
    s += constants.R * math.log(1e5 / compiled_lit["P_Pa"])
    comparison = result["literature_convention_comparison"]
    close(comparison["gas_constant_HS_at_1bar_K"], h/s)
    close(comparison["gas_constant_HS_at_1M_K"], brentq(lambda temp: h-temp*s
                                                       - constants.R*temp*math.log(1000*constants.R*temp/1e5), 250, 1800))
    # Independently solve all H/S allocation endpoints.
    chain = reactions["5"]
    def c_hs(t):
        return (chain.get_enthalpy_of_reaction(t) + constants.R * t,
                chain.get_entropy_of_reaction(t) + constants.R * (1 + math.log(1000 * constants.R * t / 1e5)))
    h0, s0 = c_hs(298.15)
    roots = {}
    for key in ("00", "01", "10", "11"):
        def balance(t):
            h, s = c_hs(t)
            return h + (-73000 - h0 if key[0] == "0" else 0) - t * (s + (-104 - s0 if key[1] == "0" else 0))
        roots[key] = brentq(balance, 250, 1800)
        close(roots[key], result["conditional_accounting"]["roots_K"][key])
    rows = result["conditional_accounting"]["rows_K"]
    close(rows["enthalpy_error_conditional"], ((roots["10"]-roots["00"])+(roots["11"]-roots["01"]))/2)
    close(rows["entropy_error_conditional"], ((roots["01"]-roots["00"])+(roots["11"]-roots["10"]))/2)
    close(sum(rows.values()), corrected_grid - 668)
    # Every numeric block must agree with freshly reproduced data.
    report = args.report.read_text()
    rendered = blocks(result)
    for name, text in rendered.items():
        begin, end = f"<!-- BEGIN I039:{name} -->", f"<!-- END I039:{name} -->"
        assert report.count(begin) == report.count(end) == 1
        actual = report.split(begin)[1].split(end)[0]
        assert actual == "\n" + text + "\n", f"report mismatch: {name}"
    print("I039 pinned snapshot and product-source hashes verified")
    print("I039 propagation baselines, L=3,4,5 thermo and equilibrium roots independently reproduced")
    print("I039 toluene LSER and conditional additive gap accounting independently reproduced")
    print(f"I039 all {len(rendered)} numeric report blocks verified")


if __name__ == "__main__":
    main()
