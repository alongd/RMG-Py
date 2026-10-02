"""Reproduce stereo/symmetry bookkeeping without changing product thermochemistry.

Command and logged invocation: ../I041_config_entropy.md. No literature ceiling
or polymerization H/S value is an input. Only the pinned I034 materializer and
baseline selectors are reused; no database tree is regenerated.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from scipy.optimize import brentq

from rmgpy import constants
from rmgpy.data.rmg import RMGDatabase
from rmgpy.kmc.compiler import EventSetCompiler, ps_proxy_set
from rmgpy.molecule.symmetry import (
    calculate_atom_symmetry_number, calculate_axis_symmetry_number,
    calculate_bond_symmetry_number, calculate_cyclic_symmetry_number,
)
from rmgpy.thermo import ThermoData


FIXTURES = Path(__file__).resolve().parent.parent
REPOSITORY = FIXTURES.parents[3]
spec = importlib.util.spec_from_file_location("i034_baseline", FIXTURES / "i034_probe/run_probe.py")
prior = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)
TEMPERATURES = (600.0, 700.0, 800.0)
C0, P0 = 1000.0, 100000.0
SOURCE_PATHS = (
    "rmgpy/molecule/symmetry.py", "rmgpy/molecule/molecule.py", "rmgpy/species.py",
    "rmgpy/data/thermo.py", "rmgpy/reaction.py", "rmgpy/kmc/compiler.py",
    "rmgpy/statmech/conformer.pyx", "test/rmgpy/kmc/fixtures/i034_probe/run_probe.py",
)


def close(actual, expected, absolute=3e-6):
    if not math.isclose(actual, expected, rel_tol=2e-10, abs_tol=absolute):
        raise AssertionError(f"{actual} != {expected}")


def progress(message):
    print("[I041] " + message, file=sys.stderr, flush=True)


def stereo_enumeration(smiles):
    """Independently enumerate tetrahedral assignments, retaining enantiomers.

    Canonical isomeric SMILES deduplicates graph-equivalent assignments, including
    meso structures. No embeddings, conformer energies, or assumed 2**n count.
    """
    # Explicit hydrogens are needed for RDKit to distinguish CH2-radical and
    # CH3 ligands at an adjacent tetrahedral centre.
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    options = StereoEnumerationOptions(onlyUnassigned=True, unique=True, maxIsomers=4096)
    structures = sorted({Chem.MolToSmiles(Chem.RemoveHs(item), isomericSmiles=True)
                         for item in EnumerateStereoisomers(molecule, options=options)})
    return {"count": len(structures), "isomeric_smiles": structures,
            "potential_centres": len(Chem.FindMolChiralCenters(
                molecule, includeUnassigned=True, useLegacyImplementation=False))}


def symmetry_audit(species):
    sigma = float(species.get_symmetry_number())
    hybrid = species.get_resonance_hybrid()
    # Stable indices in the *hybrid*, whose numerical connectivity is saved.
    atoms = [{"index": i + 1, "element": atom.element.symbol,
              "radicals": atom.radical_electrons,
              "factor": float(calculate_atom_symmetry_number(hybrid, atom))}
             for i, atom in enumerate(hybrid.atoms) if not hybrid.is_atom_in_cycle(atom)]
    bonds = []
    for i, atom in enumerate(hybrid.atoms):
        for neighbor in list(atom.edges):
            j = hybrid.atoms.index(neighbor)
            if i < j and not hybrid.is_bond_in_cycle(atom.edges[neighbor]):
                bonds.append({"indices": [i + 1, j + 1], "factor": float(
                    calculate_bond_symmetry_number(hybrid, atom, neighbor))})
    atom_product = math.prod(row["factor"] for row in atoms)
    bond_product = math.prod(row["factor"] for row in bonds)
    axis = float(calculate_axis_symmetry_number(hybrid))
    cyclic = float(calculate_cyclic_symmetry_number(hybrid)) if hybrid.is_cyclic() else 1.0
    close(atom_product * bond_product * axis * cyclic, sigma)
    centres = sum(row["factor"] == 0.5 for row in atoms)
    nonoptical = sigma * 2**centres
    return {"sigma": sigma, "centres": centres, "sigma_nonoptical": nonoptical,
            "S_nonoptical": -constants.R * math.log(nonoptical),
            "S_optical": centres * constants.R * math.log(2),
            "atom_product": atom_product, "bond_product": bond_product,
            "axis_factor": axis, "cyclic_factor": cyclic,
            "atoms": atoms, "bonds": bonds,
            "hybrid_atoms": [{"element": atom.element.symbol,
                              "radicals": atom.radical_electrons} for atom in hybrid.atoms],
            "hybrid_bonds": [{"indices": [i + 1, hybrid.atoms.index(neighbor) + 1],
                              "order": float(bond.order)}
                             for i, atom in enumerate(hybrid.atoms)
                             for neighbor, bond in atom.edges.items()
                             if i < hybrid.atoms.index(neighbor)]}


def resolve(database, entry):
    while isinstance(entry.data, str):
        entry = database.entries[entry.data]
    return entry.data


def source_terms(database, species):
    """Expand the actual source comment; require the baseline GAV route.

    Selected participants use no library, QM or fitted oligomer thermochemistry.
    HBI subtracts hydrogen-atom enthalpy but no hydrogen-atom entropy.
    """
    source = database.extract_source_from_comments(species)
    if set(source) != {"GAV"}:
        raise AssertionError(f"unexpected thermo source: {source}")
    terms = []
    for kind, entries in source["GAV"].items():
        for entry, weight in entries:
            if weight:
                terms.append((kind + ":" + entry.label, weight,
                              resolve(database.groups[kind], entry)))
    radicals = species.molecule[0].get_radical_count()
    if radicals:
        terms.append(("HBI:subtract_H_atoms", radicals,
                      ThermoData(Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
                                 Cpdata=([0.0] * 7, "J/(mol*K)"),
                                 H298=(-52.103 * 4184, "J/mol"), S298=(0, "J/(mol*K)"))))
    return [{"label": label, "weight": weight,
             "HS": [[weight * float(data.get_enthalpy(t)),
                     weight * float(data.get_entropy(t))] for t in TEMPERATURES]}
            for label, weight, data in terms]


def species_audit(database, species):
    if species.thermo is None:
        species.thermo = database.get_thermo_data(species)
    symmetry = symmetry_audit(species)
    terms = source_terms(database, species)
    rows = []
    for i, t in enumerate(TEMPERATURES):
        h, s = float(species.thermo.get_enthalpy(t)), float(species.thermo.get_entropy(t))
        group_h = sum(term["HS"][i][0] for term in terms)
        group_s = sum(term["HS"][i][1] for term in terms)
        close(h, group_h)
        close(s, group_s + symmetry["S_nonoptical"] + symmetry["S_optical"])
        rows.append({"T": t, "H": h, "S": s, "H_source": group_h, "S_source": group_s})
    smiles = species.molecule[0].to_smiles()
    return {"smiles": smiles, "adjacency": species.molecule[0].to_adjacency_list(),
            "comment": species.thermo.comment, "symmetry": symmetry, "sources": terms,
            "stereo": stereo_enumeration(smiles), "rows": rows}


def reaction_audit(database, reaction):
    # Keep participant/report ordering independent of family-generation sets.
    key = lambda s: (-s.molecule[0].get_element_count()["C"], s.molecule[0].to_smiles())
    reaction.reactants.sort(key=key)
    reaction.products.sort(key=key)
    species = [species_audit(database, item) for item in reaction.reactants + reaction.products]
    signs = [-1] * len(reaction.reactants) + [1] * len(reaction.products)
    delta_n = sum(sign * item["symmetry"]["centres"] for sign, item in zip(signs, species))
    ds_opt = delta_n * constants.R * math.log(2)
    ds_sym = sum(sign * item["symmetry"]["S_nonoptical"] for sign, item in zip(signs, species))
    rows = []
    for i, t in enumerate(TEMPERATURES):
        h, s = float(reaction.get_enthalpy_of_reaction(t)), float(reaction.get_entropy_of_reaction(t))
        h_source = sum(sign * item["rows"][i]["H_source"] for sign, item in zip(signs, species))
        s_source = sum(sign * item["rows"][i]["S_source"] for sign, item in zip(signs, species))
        close(h, h_source)
        close(s, s_source + ds_opt + ds_sym)
        nu = len(reaction.products) - len(reaction.reactants)
        h_convention = -nu * constants.R * t
        s_convention = -nu * constants.R * (1 + math.log(C0 * constants.R * t / P0))
        rows.append({"T": t, "H_p": h, "S_p": s, "H_source": h_source,
                     "S_source": s_source, "S_symmetry": ds_sym, "S_stereo": ds_opt,
                     "H_convention": h_convention, "S_convention": s_convention,
                     "H_c": h + h_convention, "S_c": s + s_convention,
                     "G_stereo": -t * ds_opt, "Kc": float(reaction.get_equilibrium_constant(t)),
                     "Kc_optical_multiplier": 2.0**delta_n})
    return {"species": species, "signs": signs, "delta_centres": delta_n,
            "delta_S_stereo": ds_opt, "delta_S_nonoptical": ds_sym, "rows": rows}


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
        "propagation_grid": prop["k_table"], "depropagation_grid": dep["k_table"]}


def cases(reaction, rates, ds):
    output = []
    for name, shift in (("atactic constitution lump (retain once)", 0.0),
                        ("one prescribed stereo channel / suppression control", -ds),
                        ("extra R ln 2 / double-counting control", ds)):
        root = brentq(lambda t: math.log(reaction.get_equilibrium_constant(t) * C0)
                      + shift / constants.R, 250, 1800, xtol=1e-9)
        forward, reverse = rates["propagation_grid"], rates["depropagation_grid"]
        residuals = [math.log(kf * C0 / kr) + shift / constants.R
                     for kf, kr in zip(forward["k"], reverse["k"])]
        crossing = None
        for i in range(len(residuals) - 1):
            if residuals[i] * residuals[i + 1] <= 0:
                a, b = forward["T"][i:i + 2]
                crossing = a - residuals[i] * (b - a) / (residuals[i + 1] - residuals[i])
                break
        if crossing is None:
            raise AssertionError("no compiler-grid crossing")
        output.append({"name": name, "shift_H": 0.0, "shift_S": shift,
                       "K_multiplier": math.exp(shift / constants.R),
                       "reverse_multiplier_at_fixed_total_forward": math.exp(-shift / constants.R),
                       "Tc_continuous": float(root), "Tc_grid": float(crossing)})
    close(output[0]["Tc_grid"], rates["ceiling_tabulated_K"])
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", type=Path, default=Path("/home/alon/Code/RMG-database"))
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    for path in (args.scratch, args.output):
        path = path.resolve()
        if path.is_relative_to(args.database.resolve()) or path.is_relative_to(REPOSITORY):
            raise ValueError("outputs must be outside repositories")
        if "catalog" in path.parts or path.is_relative_to(Path("/home/alon/Code/polymers")):
            raise ValueError("excluded output path")
    progress("materializing only the pinned database allowlist with git show")
    snapshot = prior.snapshot_database(args.database, args.scratch)
    database = RMGDatabase()
    database.load_kinetics(str(args.scratch / "input/kinetics"), reaction_libraries=[],
                           seed_mechanisms=None, kinetics_families=list(prior.PS_FAMILY_CANDIDATES),
                           kinetics_depositories=["training"])
    database.load_thermo(str(args.scratch / "input/thermo"),
                         thermo_libraries=["primaryThermoLibrary"], depository=True)
    progress("reproducing the existing gas baseline, then auditing L=3,4,5")
    reactions, rates = prior.gas_baseline(database, ps_proxy_set(3))
    sizes = {}
    for units in (3, 4, 5):
        reaction, grid = ((reactions["propagation"], rates) if units == 3
                          else compile_pair(database, units))
        item = reaction_audit(database.thermo, reaction)
        item["rates"] = {key: grid[key] for key in (
            "ceiling_tabulated_K", "propagation_grid", "depropagation_grid")}
        item["cases"] = cases(reaction, grid, item["delta_S_stereo"])
        sizes[str(units)] = item
    progress("auditing selected H-abstraction and recombination/homolysis graphs")
    families = {name: reaction_audit(database.thermo, reactions[name])
                for name in ("H_baseline", "H_chain", "homolysis")}
    # Select abstraction *at* an asymmetric backbone CH(Ph), in addition to
    # I034's terminal-benzylic control. The structural predicate requires the
    # product radical carbon to have three carbon ligands, including phenyl.
    by_site = {proxy.site_type: proxy for proxy in ps_proxy_set(3)}
    pair = (by_site["end_radical"].reactants[0], by_site["pristine"].reactants[0])
    candidates = prior.generate(database, pair, "H_Abstraction")
    choices = [reaction for reaction in candidates if any(
        atom.radical_electrons == 1
        and len([neighbor for neighbor in atom.edges if neighbor.element.symbol == "C"]) == 3
        and any(any(bond.is_benzene() for bond in neighbor.edges.values()) for neighbor in atom.edges)
        for species in reaction.products for atom in species.molecule[0].atoms)]
    selected = min(choices, key=lambda reaction: tuple(sorted(prior.smiles(reaction.products))))
    families["H_stereocentre"] = reaction_audit(database.thermo, selected)
    # Same structural end control as I039; selected by connectivity, never Tc.
    smiles = ("CC(c1ccccc1)C[CH](c1ccccc1)", "C=Cc1ccccc1",
              "CC(c1ccccc1)CC(c1ccccc1)C[CH](c1ccccc1)")
    participants = [prior.Species(molecule=[prior.Molecule(smiles=s)]) for s in smiles]
    for item in participants:
        item.generate_resonance_structures()
    benzylic = prior.Reaction(reactants=participants[:2], products=participants[2:])
    control = reaction_audit(database.thermo, benzylic)
    control["Tc_continuous"] = float(brentq(
        lambda t: math.log(benzylic.get_equilibrium_constant(t) * C0), 250, 1800))
    result = {"provenance": {"database_sha": prior.DATABASE_SHA, "snapshot": snapshot,
                             "source_sha256": {p: hashlib.sha256((REPOSITORY / p).read_bytes()).hexdigest()
                                               for p in SOURCE_PATHS},
                             "base_sha": "0bbfe8ea33781b34ec6a29d12512aaa8f73f1e45"},
              "constants": {"R": constants.R, "P0": P0, "C0": C0},
              "sizes": sizes, "families": families, "benzylic_control": control,
              "versions": {"python": sys.version.split()[0], "rdkit": Chem.rdBase.rdkitVersion},
              "general_stereo_terms": [{"delta_centres": n, "S": n * constants.R * math.log(2),
                                         "K_multiplier": 2.0**n,
                                         "G_at_T": [-t * n * constants.R * math.log(2)
                                                    for t in TEMPERATURES]}
                                        for n in (-2, -1, 0, 1, 2)]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"I041 wrote {args.output}")
    print("I041 no comparator, literature H/S, solvation, or excluded dataset is a numerical input")


if __name__ == "__main__":
    main()
