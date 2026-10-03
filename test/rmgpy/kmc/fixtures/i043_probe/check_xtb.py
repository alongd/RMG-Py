"""Tightly optimize and frequency-check CREST candidates with GFN2-xTB.

This supplies the geometry/Hessian component of the chosen composite; it
does not claim to be a DFT Hessian.
Run with rmg_env, using the report's command and persisted streams.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess

import numpy as np
from scipy.linalg import null_space
from rdkit import Chem

from run_ensembles import SPECIES, SCRATCH, allowed_scratch, read_xyz
from rotors import make_molecule, rotor_definitions, projected_modes

XTB = "/home/alon/anaconda3/envs/crest_env/bin/xtb"


def frequencies(path):
    modes = []
    for line in path.read_text().splitlines():
        match = re.match(r"\s*\d+\s+(?:[a-zA-Z]+\s+)?([-+\d.]+)\s+", line)
        if match:
            modes.append(float(match.group(1)))
    return modes


def unstable_displacement(smiles, directory):
    frame = read_xyz(directory / "xtbopt.xyz")[0]
    xyz = np.array([[float(v) for v in line.split()[1:]] for line in frame["atoms"]])
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    mass = np.array([a.GetMass() for a in mol.GetAtoms()])
    root_mass = np.repeat(np.sqrt(mass), 3)
    values = [float(v) for line in (directory / "hessian").read_text().splitlines()
              if not line.startswith("$") for v in line.split()]
    hessian = np.array(values).reshape((3*len(xyz),)*2)/root_mass[:,None]/root_mass[None,:]
    center = xyz-np.average(xyz,axis=0,weights=mass)
    rigid=[]
    for axis in np.eye(3):
        rigid.append((np.tile(axis,(len(xyz),1))*np.sqrt(mass[:,None])).ravel())
        rigid.append((np.cross(center,axis)*np.sqrt(mass[:,None])).ravel())
    basis = null_space(np.array(rigid))
    eigenvalues, eigenvectors=np.linalg.eigh(basis.T@hessian@basis)
    displacement=(basis@eigenvectors[:,0]/root_mass).reshape((-1,3))
    displacement/=np.max(np.linalg.norm(displacement,axis=1))
    return xyz, displacement


def projected_minimum(smiles, directory):
    mol=make_molecule(smiles, directory / "xtbopt.xyz")
    values=[float(v) for line in (directory / "hessian").read_text().splitlines()
            if not line.startswith("$") for v in line.split()]
    hessian=np.array(values).reshape((3*mol.GetNumAtoms(),)*2)
    try:
        projected_modes(mol,rotor_definitions(mol),hessian)
        return True
    except AssertionError:
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch", type=Path, default=SCRATCH)
    parser.add_argument("--species", choices=tuple(SPECIES), nargs="+", required=True)
    parser.add_argument("--indices", type=int, nargs="+")
    parser.add_argument("--starting-input", action="store_true",
                        help="pilot on the prepared input instead of a completed CREST ensemble")
    args = parser.parse_args()
    scratch = allowed_scratch(args.scratch)
    physical = {}
    for cpu in sorted(os.sched_getaffinity(0)):
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        core = (topology / "physical_package_id").read_text(), (topology / "core_id").read_text()
        physical.setdefault(core, cpu)
    os.sched_setaffinity(0, sorted(physical.values())[:8])
    env = dict(os.environ, OMP_NUM_THREADS="8", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               OMP_MAX_ACTIVE_LEVELS="1")
    for name in args.species:
        source = scratch / "ensembles" / name
        if args.starting_input:
            frames = read_xyz(source / "input.xyz")
            indices = [0]
        else:
            candidates = json.loads((source / "completed.json").read_text())
            frames = read_xyz(source / "crest_conformers.xyz")
            indices = args.indices if args.indices is not None else candidates["selected_indices_within_12_kJ"]
        for index in indices:
            frame = frames[index]
            directory = scratch / ("xtb_pilots" if args.starting_input else "xtb_checks") / name / f"{index:04d}"
            directory.mkdir(parents=True, exist_ok=True)
            if (directory / "result.json").exists():
                old = json.loads((directory / "result.json").read_text())
                if min(old["frequencies_cm1"]) >= -1.0 and projected_minimum(SPECIES[name], directory):
                    continue
            xyz = str(len(frame["atoms"])) + "\n" + frame["comment"] + "\n" + "\n".join(frame["atoms"]) + "\n"
            (directory / "input.xyz").write_text(xyz)
            command = [XTB, "input.xyz", "--gfn", "2", "--ohess", "extreme", "--acc", "0.1", "--parallel", "8"]
            shell = shlex.join(command) + " > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)"
            subprocess.run(["bash", "-o", "pipefail", "-c", shell], cwd=directory, env=env, check=True)
            modes = frequencies(directory / "vibspectrum")
            retry=0
            while (min(modes) < -1.0 or not projected_minimum(SPECIES[name], directory)) and retry < 4:
                coordinates, displacement=unstable_displacement(SPECIES[name], directory)
                amplitude=0.20*(retry+1)
                coordinates += amplitude*displacement
                perturbed = str(len(frame["atoms"])) + "\nnegative-mode displacement\n"
                for line, row in zip(frame["atoms"], coordinates):
                    perturbed += line.split()[0]+" "+" ".join(f"{value:.12f}" for value in row)+"\n"
                (directory / f"negative_mode_retry_{retry}.xyz").write_text(perturbed)
                (directory / "input.xyz").write_text(perturbed)
                subprocess.run(["bash", "-o", "pipefail", "-c", shell], cwd=directory, env=env, check=True)
                modes=frequencies(directory / "vibspectrum")
                retry+=1
            if len(modes) != 3*len(frame["atoms"]):
                raise AssertionError("wrong number of xTB modes")
            optimized = read_xyz(directory / "xtbopt.xyz")[0]
            match = re.search(r"energy:\s*([-+\d.Ee]+)", optimized["comment"])
            if match is None:
                raise ValueError("xTB energy missing from optimized geometry")
            energy = float(match.group(1))
            record = {"name": name, "index": index, "command": command,
                      "input_sha256": hashlib.sha256(xyz.encode()).hexdigest(),
                      "energy_Eh": energy, "frequencies_cm1": modes,
                      "minimum_confirmed": min(modes) >= -1.0 and projected_minimum(SPECIES[name], directory),
                      "negative_mode_retries": retry,
                      "xyz_sha256": hashlib.sha256((directory / "xtbopt.xyz").read_bytes()).hexdigest()}
            (directory / "result.json").write_text(json.dumps(record, indent=2) + "\n")
            print(f"I043 xTB check {name} {index}: lowest frequency {min(modes):.2f} cm-1; "
                  f"minimum={record['minimum_confirmed']}", flush=True)


if __name__ == "__main__":
    main()
