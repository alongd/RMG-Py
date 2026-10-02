"""Independently reconstruct all measured numbers and check the report.

Command is in ../I041_config_entropy.md. Checks chemistry, stereo enumeration,
source decomposition, rate/Kc ratios, roots, and a microchannel balance model.
It does not claim to validate tacticity-independent molecular free energies.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess

from rdkit import Chem
from scipy.optimize import brentq

from rmgpy import constants
from rmgpy.data.rmg import RMGDatabase
from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.symmetry import (
    calculate_atom_symmetry_number, calculate_axis_symmetry_number,
    calculate_bond_symmetry_number, calculate_cyclic_symmetry_number,
)
from rmgpy.reaction import Reaction
from rmgpy.species import Species

from render_tables import REPORT, blocks


REPOSITORY = Path(__file__).resolve().parents[5]
DATABASE_SHA = "4a12d36fcdc193ede82c8d1ab5c1653495d445bc"


def close(actual, expected, absolute=3e-6):
    assert math.isclose(actual, expected, rel_tol=2e-10, abs_tol=absolute), (actual, expected)


def brute_force_stereo(smiles):
    """Cartesian tetrahedral assignments, independent of EnumerateStereoisomers."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    positions = [item.centeredOn for item in Chem.FindPotentialStereo(mol)
                 if item.type == Chem.StereoType.Atom_Tetrahedral]
    isomers = set()
    for tags in itertools.product((Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                                   Chem.ChiralType.CHI_TETRAHEDRAL_CCW), repeat=len(positions)):
        item = Chem.Mol(mol)
        for index, tag in zip(positions, tags):
            item.GetAtomWithIdx(index).SetChiralTag(tag)
        isomers.add(Chem.MolToSmiles(Chem.RemoveHs(item), isomericSmiles=True))
    return sorted(isomers)


def check_species(database, saved):
    item = Species(molecule=[Molecule().from_adjacency_list(saved["adjacency"])])
    item.thermo = database.get_thermo_data(item)
    assert item.molecule[0].to_smiles() == saved["smiles"]
    assert item.thermo.comment == saved["comment"]
    symmetry = saved["symmetry"]
    hybrid = item.get_resonance_hybrid()
    factors = [float(calculate_atom_symmetry_number(hybrid, atom))
               for atom in hybrid.atoms if not hybrid.is_atom_in_cycle(atom)]
    atom_product = math.prod(factors)
    # Atom ordering is allowed to change, so compare the multiset of factors.
    assert sorted(factors) == sorted(row["factor"] for row in symmetry["atoms"])
    bonds = []
    for i, atom in enumerate(hybrid.atoms):
        for neighbor in list(atom.edges):
            if i < hybrid.atoms.index(neighbor) and not hybrid.is_bond_in_cycle(atom.edges[neighbor]):
                bonds.append(float(calculate_bond_symmetry_number(hybrid, atom, neighbor)))
    assert sorted(bonds) == sorted(row["factor"] for row in symmetry["bonds"])
    bond_product = math.prod(bonds)
    axis = float(calculate_axis_symmetry_number(hybrid))
    cyclic = float(calculate_cyclic_symmetry_number(hybrid)) if hybrid.is_cyclic() else 1.0
    n = factors.count(0.5)
    for actual, key in ((item.get_symmetry_number(), "sigma"), (atom_product, "atom_product"),
                        (bond_product, "bond_product"), (axis, "axis_factor"),
                        (cyclic, "cyclic_factor"), (n, "centres")):
        close(actual, symmetry[key])
    close(atom_product * bond_product * axis * cyclic, symmetry["sigma"])
    close(item.get_symmetry_number() * 2**n, symmetry["sigma_nonoptical"])
    close(-constants.R * math.log(symmetry["sigma_nonoptical"]), symmetry["S_nonoptical"])
    close(constants.R * n * math.log(2), symmetry["S_optical"])
    assert brute_force_stereo(saved["smiles"]) == saved["stereo"]["isomeric_smiles"]
    assert len(saved["stereo"]["isomeric_smiles"]) == saved["stereo"]["count"]
    # Calling the uncorrected group estimator independently checks that no
    # additional optical/symmetry term is concealed in the source sum.
    raw = database.get_thermo_data_from_groups(item.copy(deep=True))
    extracted = database.extract_source_from_comments(item)
    expected_weights = {}
    for kind, entries in extracted["GAV"].items():
        for entry, weight in entries:
            if weight:
                expected_weights[kind + ":" + entry.label] = weight
    radicals = item.molecule[0].get_radical_count()
    if radicals:
        expected_weights["HBI:subtract_H_atoms"] = radicals
    assert {term["label"]: term["weight"] for term in saved["sources"]} == expected_weights
    for i, row in enumerate(saved["rows"]):
        t = row["T"]
        close(item.thermo.get_enthalpy(t), row["H"])
        close(item.thermo.get_entropy(t), row["S"])
        close(raw.get_enthalpy(t), row["H_source"])
        close(raw.get_entropy(t), row["S_source"])
        close(sum(term["HS"][i][0] for term in saved["sources"]), row["H_source"])
        close(sum(term["HS"][i][1] for term in saved["sources"]), row["S_source"])
        close(row["S_source"] - constants.R * math.log(symmetry["sigma"]), row["S"])
        for term in saved["sources"]:
            if term["label"] == "HBI:subtract_H_atoms":
                close(term["HS"][i][0], -52.103 * 4184 * radicals)
                close(term["HS"][i][1], 0)
                continue
            kind, label = term["label"].split(":", 1)
            data = database.groups[kind].entries[label].data
            while isinstance(data, str):
                data = database.groups[kind].entries[data].data
            close(data.get_enthalpy(t) * term["weight"], term["HS"][i][0])
            close(data.get_entropy(t) * term["weight"], term["HS"][i][1])
    return item


def check_reaction(database, saved):
    participants = [check_species(database, row) for row in saved["species"]]
    nr = saved["signs"].count(-1)
    assert saved["signs"] == [-1] * nr + [1] * (len(participants) - nr)
    reaction = Reaction(reactants=participants[:nr], products=participants[nr:])
    n = sum(sign * row["symmetry"]["centres"] for sign, row in zip(saved["signs"], saved["species"]))
    ds_stereo = constants.R * n * math.log(2)
    ds_sym = -constants.R * sum(sign * math.log(row["symmetry"]["sigma_nonoptical"])
                                for sign, row in zip(saved["signs"], saved["species"]))
    close(n, saved["delta_centres"])
    close(ds_stereo, saved["delta_S_stereo"])
    close(ds_sym, saved["delta_S_nonoptical"])
    nu = len(reaction.products) - nr
    for i, row in enumerate(saved["rows"]):
        t = row["T"]
        h, s = reaction.get_enthalpy_of_reaction(t), reaction.get_entropy_of_reaction(t)
        hc = h - nu * constants.R * t
        sc = s - nu * constants.R * (1 + math.log(1000 * constants.R * t / 1e5))
        values = {"H_p": h, "S_p": s, "H_source": h,
                  "S_source": s - ds_stereo - ds_sym,
                  "S_stereo": ds_stereo, "S_symmetry": ds_sym,
                  "H_c": hc, "S_c": sc, "H_convention": hc - h, "S_convention": sc - s,
                  "G_stereo": -t * ds_stereo, "Kc_optical_multiplier": 2.0**n}
        for key, value in values.items():
            close(value, row[key])
        kc = math.exp(-(hc - t * sc) / (constants.R * t)) * 1000**nu
        close(kc, row["Kc"], absolute=1e-12)
        close(kc, reaction.get_equilibrium_constant(t), absolute=1e-12)
    if "rates" in saved:
        rates = saved["rates"]
        fwd, rev = rates["propagation_grid"], rates["depropagation_grid"]
        assert fwd["T"] == rev["T"]
        for t, kf, kr in zip(fwd["T"], fwd["k"], rev["k"]):
            close(kf / kr, reaction.get_equilibrium_constant(t), absolute=1e-12)
        assert n == 1
        assert [row["stereo"]["count"] for row in saved["species"]] == [
            2**saved["species"][0]["symmetry"]["centres"], 1,
            2**saved["species"][2]["symmetry"]["centres"]]
        for case, shift in zip(saved["cases"], (0.0, -ds_stereo, ds_stereo)):
            close(case["shift_H"], 0)
            close(case["shift_S"], shift)
            close(case["K_multiplier"], math.exp(shift / constants.R))
            close(case["reverse_multiplier_at_fixed_total_forward"], math.exp(-shift / constants.R))
            root = brentq(lambda t: reaction.get_free_energy_of_reaction(t) - t * shift
                          - constants.R * t * math.log(1000 * constants.R * t / 1e5),
                          250, 1800, xtol=1e-9)
            close(root, case["Tc_continuous"])
            y = [math.log(reaction.get_equilibrium_constant(t) * 1000) + shift / constants.R
                 for t in fwd["T"]]
            for i in range(len(y) - 1):
                if y[i] * y[i + 1] <= 0:
                    crossing = fwd["T"][i] - y[i] * (fwd["T"][i + 1] - fwd["T"][i]) / (y[i + 1] - y[i])
                    close(crossing, case["Tc_grid"])
                    break
            else:
                raise AssertionError("missing grid crossing")
        close(rates["ceiling_tabulated_K"], saved["cases"][0]["Tc_grid"])
    return reaction


def microchannel_check(kc_lump):
    """Fixed, even nonuniform prefixes: two daughters and one reverse each.

    Arbitrary normalized parent populations suffice because rates do not depend
    on the prefix under the stated ideal model. Both steady-state fluxes and
    arbitrary population derivatives are checked, not just state counts.
    """
    kmicro = kc_lump / 2
    for m in (1, 2, 3, 4):
        histories = list(itertools.product((0, 1), repeat=m))
        total_weight = sum(range(1, len(histories) + 1))
        parents = {history: (i + 1) / total_weight for i, history in enumerate(histories)}
        children = {history + (s,): kmicro * concentration * 1000
                    for history, concentration in parents.items() for s in (0, 1)}
        kbranch, kreverse = 1.0, 1.0 / kmicro
        for daughter, concentration in children.items():
            close(kbranch * parents[daughter[:-1]] * 1000, kreverse * concentration)
        close(sum(children.values()) / (sum(parents.values()) * 1000), kc_lump)
        # A non-equilibrium product population still closes exactly when summed.
        arbitrary_children = {history: concentration * (i + 1) / len(children)
                              for i, (history, concentration) in enumerate(children.items())}
        summed_flux = sum(kbranch * parents[history[:-1]] * 1000 - kreverse * concentration
                          for history, concentration in arbitrary_children.items())
        coarse_flux = 2 * kbranch * sum(parents.values()) * 1000 - kreverse * sum(arbitrary_children.values())
        close(summed_flux, coarse_flux)
        close((2 * kbranch) / kreverse, kc_lump)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--database", type=Path, default=Path("/home/alon/Code/RMG-database"))
    parser.add_argument("--report", type=Path, default=REPORT)
    args = parser.parse_args()
    result = json.loads(args.results.read_text())
    assert result["provenance"]["database_sha"] == DATABASE_SHA
    close(result["constants"]["R"], constants.R)
    close(result["constants"]["P0"], 1e5)
    close(result["constants"]["C0"], 1000)
    paths = sorted(p for p in args.snapshot.rglob("*")
                   if p.is_file() and "__pycache__" not in p.parts)
    allowed = subprocess.check_output([
        "git", "-C", str(args.database), "ls-tree", "-r", "--name-only", DATABASE_SHA, "--",
        *result["provenance"]["snapshot"]["prefixes"]], text=True).splitlines()
    assert [p.relative_to(args.snapshot).as_posix() for p in paths] == allowed
    assert len(paths) == result["provenance"]["snapshot"]["files"]
    digest = hashlib.sha256()
    for path in paths:
        relative = path.relative_to(args.snapshot).as_posix()
        assert "catalog" not in path.parts
        expected = subprocess.check_output(["git", "-C", str(args.database), "show", f"{DATABASE_SHA}:{relative}"])
        assert path.read_bytes() == expected
        digest.update(relative.encode() + b"\0" + expected)
    assert digest.hexdigest() == result["provenance"]["snapshot"]["sha256"]
    for path, expected in result["provenance"]["source_sha256"].items():
        assert hashlib.sha256((REPOSITORY / path).read_bytes()).hexdigest() == expected
    db = RMGDatabase()
    db.load_thermo(str(args.snapshot / "input/thermo"), thermo_libraries=["primaryThermoLibrary"], depository=True)
    reactions = {units: check_reaction(db.thermo, item) for units, item in result["sizes"].items()}
    for item in result["families"].values():
        check_reaction(db.thermo, item)
    benzylic = check_reaction(db.thermo, result["benzylic_control"])
    close(brentq(lambda t: math.log(benzylic.get_equilibrium_constant(t) * 1000), 250, 1800),
          result["benzylic_control"]["Tc_continuous"])
    for row in result["general_stereo_terms"]:
        ds = row["delta_centres"] * constants.R * math.log(2)
        close(row["S"], ds)
        close(row["K_multiplier"], 2.0**row["delta_centres"])
        for t, g in zip((600, 700, 800), row["G_at_T"]):
            close(g, -t * ds)
    microchannel_check(reactions["3"].get_equilibrium_constant(700))
    # Reproduction anchors from the antecedent report, never selection targets.
    close(result["sizes"]["3"]["cases"][0]["Tc_continuous"], 710.020462, absolute=1e-6)
    close(result["sizes"]["3"]["cases"][1]["Tc_continuous"], 672.121955, absolute=1e-6)
    report = args.report.read_text()
    rendered = blocks(result)
    for name, body in rendered.items():
        begin, end = f"<!-- BEGIN I041:{name} -->", f"<!-- END I041:{name} -->"
        assert report.count(begin) == report.count(end) == 1
        assert report.split(begin)[1].split(end)[0] == "\n" + body + "\n", name
    print("I041 pinned git-show snapshot and product-source hashes verified")
    print("I041 species thermo, source terms, symmetry primitives and independent stereo enumeration verified")
    print("I041 L=3,4,5 rate/Kc ratios, concentration conversion and all gas Tc controls reproduced")
    print("I041 selected H-abstraction, homolysis and benzylic-end numbers reproduced")
    print("I041 fixed-prefix microchannel detailed balance and exact population lumping verified")
    print(f"I041 all {len(rendered)} numeric report blocks verified")


if __name__ == "__main__":
    main()
