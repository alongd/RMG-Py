"""Run: PYTHONPATH=$PWD python test/rmgpy/kmc/fixtures/i049_probe/rates.py

Only two radical+styrene calls, with product filtering; no event compilation.
Verify every loaded database file against git show at the dispatch pin.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import time

from common import benzylic_tail, canonical_smiles, describe
from rmgpy import constants
from rmgpy.data.rmg import RMGDatabase
from rmgpy.kmc.compiler import _rate_rule_source, prepare_rate_rules
from rmgpy.molecule.molecule import Molecule
from rmgpy.species import Species

PIN = "4a12d36fcdc193ede82c8d1ab5c1653495d445bc"
DEFAULT_DATABASE = "/home/alon/runs/i046-rules-from-training/database"
REPOSITORY = "/home/alon/Code/RMG-database"


def verify_snapshot(database):
    roots = [database / "input/thermo", database / "input/kinetics/families/R_Addition_MultipleBond"]
    files = [database / "input/kinetics/families/recommended.py"]
    files += [p for root in roots for p in sorted(root.rglob("*"))
              if p.is_file() and p.suffix in {".py", ".txt"}]
    verified = []
    for path in files:
        if not path.is_file():
            continue
        relative = str(path.relative_to(database))
        expected = subprocess.check_output(["git", "-C", REPOSITORY, "show", f"{PIN}:{relative}"])
        actual = path.read_bytes()
        assert actual == expected, relative
        verified.append({"path": relative, "sha256": hashlib.sha256(actual).hexdigest()})
    assert verified
    return verified


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, default=Path(DEFAULT_DATABASE))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    start = time.monotonic()
    files = verify_snapshot(args.database)
    print(f"verified {len(files)} pinned database files; loading one family", flush=True)
    db = RMGDatabase()
    db.load_kinetics(str(args.database / "input/kinetics"), reaction_libraries=[],
                     seed_mechanisms=None, kinetics_families=["R_Addition_MultipleBond"],
                     kinetics_depositories=["training"])
    db.load_thermo(str(args.database / "input/thermo"),
                   thermo_libraries=["primaryThermoLibrary"], depository=True)
    preparation = prepare_rate_rules(db.kinetics, db.thermo, verbose=True)
    family = db.kinetics.families["R_Addition_MultipleBond"]
    rows = []
    for units in (2, 3):
        reactants = [Species(molecule=[Molecule(smiles=s)])
                     for s in (benzylic_tail(units), "C=Cc1ccccc1")]
        target = Molecule(smiles=benzylic_tail(units + 1))
        reactions = db.kinetics.generate_reactions_from_families(
            reactants, products=[target], only_families=[family.label], resonance=True)
        assert len(reactions) == 1, (units, len(reactions))
        reaction = reactions[0]
        assert len(reaction.products) == 1
        product_smiles = [describe(m.to_adjacency_list(remove_h=False))["smiles"]
                          for m in reaction.products[0].molecule]
        # Family generation aromatizes benzene; the raw builder uses Kekule bonds.
        assert canonical_smiles(benzylic_tail(units + 1)) in product_smiles, product_smiles
        kinetics, source, entry, forward = family.get_kinetics(
            reaction, template_labels=reaction.template, degeneracy=reaction.degeneracy,
            return_all_kinetics=False)
        assert forward, "Unexpected reverse estimate: report separately before changing orientation"
        reaction.kinetics = kinetics
        training_origin = re.search(r"From training reaction (\d+)", kinetics.comment)
        assert training_origin, kinetics.comment
        original = family.get_training_depository().entries[int(training_origin.group(1))]
        data = {"units": units, "reactants": [s.molecule[0].to_smiles() for s in reaction.reactants],
                "product": target.to_smiles(), "template": reaction.template,
                "degeneracy": reaction.degeneracy, "source": source,
                "entry": str(entry) if entry else None,
                "rate_source": _rate_rule_source(family, reaction),
                "training_origin": {"index": original.index, "label": original.label,
                                    "rank": original.rank, "short_desc": original.short_desc,
                                    "long_desc": original.long_desc,
                                    "kinetics": repr(original.data)},
                "kinetics": repr(kinetics), "grid": []}
        for temperature in (600., 700., 800.):
            # Exactly the compiler's _rate_table call; no selection/conversion change.
            rate = kinetics.get_rate_coefficient(temperature) * 1000
            benchmark = 4.27e7 * math.exp(-32500 / (constants.R * temperature))
            assert rate > 0 and math.isfinite(rate)
            data["grid"].append({"T_K": temperature, "RMG_L_mol_s": rate,
                                 "IUPAC_L_mol_s": benchmark, "ratio": rate / benchmark})
        rows.append(data)
        print(json.dumps(data, indent=2), flush=True)
    result = {"database_sha": PIN, "verified_files": files,
              "preparation": preparation, "reactions": rows,
              "elapsed_s": time.monotonic() - start}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS: two head-to-tail reactions estimated after pinned training preparation", flush=True)


if __name__ == "__main__":
    main()
