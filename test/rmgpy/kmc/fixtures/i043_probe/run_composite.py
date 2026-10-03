"""Resume tight minimum checks and both composite energy levels serially.

The report gives the command. The wall-time cap includes waiting for searches.
This creates molecular electronic/frequency evidence, not rotor-corrected thermo.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import time

from rdkit import Chem
from rdkit.Chem import rdMolAlign

from run_ensembles import SPECIES, SCRATCH, allowed_scratch, read_xyz
from summarize import stereo_check
from complete_symmetry import augment, materialize_energy

RMG_PYTHON = "/home/alon/anaconda3/envs/rmg_env/bin/python"
QUANTUM_PYTHON = "/home/alon/anaconda3/envs/pyscf_env/bin/python"
HERE = Path(__file__).resolve().parent


def heavy_molecule(smiles, frame):
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for index, line in enumerate(frame["atoms"]):
        conformer.SetAtomPosition(index, [float(v) for v in line.split()[1:]])
    molecule.AddConformer(conformer)
    return Chem.RemoveHs(molecule)


def select_minima(scratch, name):
    source = scratch / "ensembles" / name
    ensemble = json.loads((source / "completed.json").read_text())
    selected = []
    checked = []
    for index in ensemble["selected_indices_within_12_kJ"]:
        directory = scratch / "xtb_checks" / name / f"{index:04d}"
        record = json.loads((directory / "result.json").read_text())
        if not record["minimum_confirmed"]:
            continue
        frame = read_xyz(directory / "xtbopt.xyz")[0]
        stereo_check(SPECIES[name], frame)
        molecule = heavy_molecule(SPECIES[name], frame)
        duplicate = None
        for earlier in selected:
            if abs(record["energy_Eh"] - earlier["energy_Eh"]) * 2625.499639 < 0.05:
                if rdMolAlign.GetBestRMS(molecule, earlier["molecule"], maxMatches=10000) < 0.10:
                    duplicate = earlier["index"]
                    break
        checked.append({"index": index, "duplicate_of": duplicate, "energy_Eh": record["energy_Eh"]})
        if duplicate is None:
            selected.append({"index": index, "energy_Eh": record["energy_Eh"], "molecule": molecule})
    if not selected:
        raise RuntimeError(f"no stable minima for {name}")
    result = {"name": name, "checked": checked, "unique_minima_indices": [row["index"] for row in selected],
              "deduplication": {"heavy_atom_RMS_A": 0.10, "energy_kJ_mol": 0.05,
                                "methyl_rotamers": "not separate conformers; require rotor treatment"}}
    result = augment(scratch,name,result)
    target = scratch / "composite" / name
    target.mkdir(parents=True, exist_ok=True)
    (target / "selection.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch", type=Path, default=SCRATCH)
    parser.add_argument("--hours", type=float, default=20.)
    parser.add_argument("--species", nargs="+", choices=tuple(SPECIES))
    parser.add_argument("--minima-only", action="store_true",
                        help="prepare validated minimum selections while electronic jobs run separately")
    args = parser.parse_args()
    scratch = allowed_scratch(args.scratch)
    physical = {}
    for cpu in sorted(os.sched_getaffinity(0)):
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        core = (topology / "physical_package_id").read_text(), (topology / "core_id").read_text()
        physical.setdefault(core, cpu)
    os.sched_setaffinity(0, sorted(physical.values())[:8])
    env = dict(os.environ, OMP_NUM_THREADS="8", OPENBLAS_NUM_THREADS="8", MKL_NUM_THREADS="8",
               OMP_MAX_ACTIVE_LEVELS="1")
    started = time.monotonic()

    def remaining():
        seconds = int(args.hours*3600 - (time.monotonic() - started))
        if seconds <= 0:
            raise TimeoutError("composite wall-time budget exhausted")
        return seconds

    def run(command, directory):
        directory.mkdir(parents=True, exist_ok=True)
        command = ["timeout", "--kill-after=30", str(remaining()), *command]
        shell = shlex.join(command) + " > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)"
        subprocess.run(["bash", "-o", "pipefail", "-c", shell], cwd=directory, env=env, check=True)

    for name in args.species or SPECIES:
        marker = scratch / "ensembles" / name / "completed.json"
        print(f"I043 composite awaiting completed search: {name}", flush=True)
        while not marker.exists():
            time.sleep(min(30, remaining()))
        run([RMG_PYTHON, str(HERE / "check_xtb.py"), "--scratch", str(scratch), "--species", name],
            scratch / "composite" / name / "checks")
        selection = select_minima(scratch, name)
        print(f"I043 composite {name}: {len(selection['unique_minima_indices'])} unique checked minima", flush=True)
        if args.minima_only:
            continue
        for index in selection["unique_minima_indices"]:
            source = scratch / "xtb_checks" / name / f"{index:04d}" / "xtbopt.xyz"
            for level in ("pbe", "blyp"):
                if materialize_energy(scratch,name,index,level):
                    continue
                directory = scratch / "composite" / name / f"{index:04d}" / level
                output = directory / "result.json"
                if output.exists():
                    cached = json.loads(output.read_text())
                    if cached["input_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest():
                        continue
                print(f"I043 composite calculating {name} {index} {level}", flush=True)
                run([QUANTUM_PYTHON, str(HERE / "refine.py"), "--xyz", str(source),
                     "--output", str(output), "--level", level, "--single-point"], directory)
        print(f"I043 composite completed both electronic levels for {name}", flush=True)
    print("I043 composite completed requested minimum selections" if args.minima_only else
          "I043 composite completed all requested molecular electronic/frequency calculations", flush=True)


if __name__ == "__main__":
    main()
