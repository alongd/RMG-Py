"""Generate the dispatched GFN2-xTB ensembles; commands are in the report.

Run with rmg_env (RDKit). All external jobs run serially with eight threads.
No thermochemistry or ceiling-temperature comparator selects conformers.
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

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

CREST = "/home/alon/anaconda3/envs/crest_env/bin/crest"
SCRATCH = Path("/home/alon/runs/i043-diphenyl-nonadditivity")
SPECIES = {
    "ethane": "CC",
    "ethylbenzene": "CCc1ccccc1",
    "n_propylbenzene": "CCCc1ccccc1",
    "diphenylpropane": "c1ccccc1CCCc1ccccc1",
    "diphenylpentane_meso": "C[C@H](c1ccccc1)C[C@H](c1ccccc1)C",
    "diphenylpentane_racemo": "C[C@H](c1ccccc1)C[C@@H](c1ccccc1)C",
    "triphenylheptane_iso": "C[C@H](c1ccccc1)C[C@H](c1ccccc1)C[C@H](c1ccccc1)C",
    "triphenylheptane_syndio": "C[C@H](c1ccccc1)C[C@@H](c1ccccc1)C[C@H](c1ccccc1)C",
    "triphenylheptane_hetero": "C[C@H](c1ccccc1)C[C@H](c1ccccc1)C[C@@H](c1ccccc1)C",
    # Supporting one-ring reference: balances the trimer-minus-dimer increment.
    "cumene": "CC(C)c1ccccc1",
}
PLAN = {
    "search": "CREST 3.0.2, GFN2-xTB, --quick, 6 kcal/mol ensemble window",
    "refinement_levels": ["PBE-D3(BJ)/def2-SVP//GFN2-xTB", "BLYP-D3(BJ)/def2-SVP//GFN2-xTB"],
    "refinement_window_kJ_mol": 12.0,
    "frequencies": "GFN2-xTB minimum Hessians; audited vtight minima retained, new/repaired minima extreme; separate DFT pilots",
    "minima_acceptance": "no mode below -1 cm-1; positive projected torsion/vibration Hessian",
    "sampling_seed": 43002,
    "cores": 8,
    "selection": "energy window only; no Tc or external comparator used",
    "state": "ideal gas, 100000 Pa; monomer 1000 mol/m3; I039 conventions",
    "rotors": "must be resolved before reporting ensemble H/S or corrected Tc",
    "stereo": "fixed random dyads; meso/racemo 1/2 each; mm,rr,mr/rm 1/4,1/4,1/2; retain R ln 2 once",
}


def allowed_scratch(path):
    path = path.resolve()
    if path != SCRATCH and SCRATCH not in path.parents:
        raise ValueError("scratch must be inside the dispatched scratch directory")
    path.mkdir(parents=True, exist_ok=True)
    return path


def limit_cores():
    physical = {}
    for cpu in sorted(os.sched_getaffinity(0)):
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        key = (topology / "physical_package_id").read_text(), (topology / "core_id").read_text()
        physical.setdefault(key,cpu)
    os.sched_setaffinity(0,sorted(physical.values())[:8])


def read_xyz(path):
    """Read concatenated CREST XYZ frames, retaining the energy comments."""
    lines = path.read_text().splitlines()
    frames = []
    offset = 0
    while offset < len(lines):
        if not lines[offset].strip():
            offset += 1
            continue
        n = int(lines[offset])
        comment = lines[offset + 1]
        atoms = lines[offset + 2:offset + 2 + n]
        if len(atoms) != n or any(len(line.split()) != 4 for line in atoms):
            raise ValueError(f"incomplete XYZ frame in {path}")
        frames.append({"comment": comment, "atoms": atoms})
        offset += n + 2
    return frames


def prepare(scratch, update_plan=False):
    plan = dict(PLAN, rdkit_version=rdBase.rdkitVersion)
    plan_path = scratch / "method_plan.json"
    if plan_path.exists() and json.loads(plan_path.read_text()) != plan:
        if not update_plan:
            raise ValueError("existing method plan differs; explicit --update-plan required")
        archive = scratch / "method_plan_initial.json"
        if not archive.exists():
            archive.write_bytes(plan_path.read_bytes())
    plan_path.write_text(json.dumps(plan, indent=2) + "\n")
    manifest = {}
    for name, smiles in SPECIES.items():
        directory = scratch / "ensembles" / name
        directory.mkdir(parents=True, exist_ok=True)
        molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = PLAN["sampling_seed"]
        parameters.numThreads = 1
        if AllChem.EmbedMolecule(molecule, parameters) != 0:
            raise RuntimeError(f"embedding failed: {name}")
        if AllChem.MMFFOptimizeMolecule(molecule, maxIters=2000) != 0:
            raise RuntimeError(f"starting geometry optimization failed: {name}")
        xyz = Chem.MolToXYZBlock(molecule)
        target = directory / "input.xyz"
        if target.exists() and target.read_text() != xyz:
            raise ValueError(f"starting geometry changed: {name}")
        target.write_text(xyz)
        mirror = Chem.Mol(molecule)
        for atom in mirror.GetAtoms():
            if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                atom.InvertChirality()
        canonical = Chem.MolToSmiles(Chem.RemoveHs(molecule))
        mirror_smiles = Chem.MolToSmiles(Chem.RemoveHs(mirror))
        centers = Chem.FindMolChiralCenters(molecule, includeUnassigned=True)
        manifest[name] = {
            "smiles": smiles, "canonical_smiles": canonical,
            "mirror_smiles": mirror_smiles, "achiral": canonical == mirror_smiles,
            "meso": bool(centers) and canonical == mirror_smiles,
            "natoms": molecule.GetNumAtoms(),
            "centers": centers,
            "input_sha256": hashlib.sha256(xyz.encode()).hexdigest(),
        }
    (scratch / "species.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch", type=Path, default=SCRATCH)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--update-plan", action="store_true")
    parser.add_argument("--species", choices=tuple(SPECIES), nargs="+")
    parser.add_argument("--hours", type=float, default=20.0)
    args = parser.parse_args()
    # Limit the actual CPU set as well as thread variables: nested BLAS/OpenMP
    # runtimes in external tools otherwise can oversubscribe --T.
    physical = {}
    for cpu in sorted(os.sched_getaffinity(0)):
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        core = (topology / "physical_package_id").read_text(), (topology / "core_id").read_text()
        physical.setdefault(core, cpu)
    os.sched_setaffinity(0, sorted(physical.values())[:8])
    scratch = allowed_scratch(args.scratch)
    prepare(scratch, args.update_plan)
    if args.prepare_only:
        print(f"I043 prepared {len(SPECIES)} input geometries and fixed method plan", flush=True)
        return
    env = dict(os.environ, OMP_NUM_THREADS="8", MKL_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", OMP_MAX_ACTIVE_LEVELS="1")
    env["PATH"] = str(Path(CREST).parent) + ":" + env.get("PATH", "")
    started = time.monotonic()
    for name in args.species or SPECIES:
        directory = scratch / "ensembles" / name
        marker = directory / "completed.json"
        if marker.exists():
            print(f"I043 already completed {name}", flush=True)
            continue
        remaining = int(args.hours * 3600 - (time.monotonic() - started))
        if remaining <= 0:
            raise TimeoutError("ensemble wall-time budget exhausted")
        command = ["timeout", "--kill-after=30", str(remaining), CREST,
                   "input.xyz", "--gfn2", "--quick", "--ewin", "6", "--T", "8"]
        launched = time.monotonic()
        print(f"I043 start {name}: {shlex.join(command)}", flush=True)
        shell_command = shlex.join(command) + " > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)"
        result = subprocess.run(["bash", "-o", "pipefail", "-c", shell_command],
                                cwd=directory, env=env)
        elapsed = time.monotonic() - launched
        if result.returncode != 0:
            raise RuntimeError(f"CREST failed for {name}: exit {result.returncode}")
        frames = read_xyz(directory / "crest_conformers.xyz")
        energies = [float(frame["comment"].split()[0]) for frame in frames]
        low = min(energies)
        selected = [i for i, e in enumerate(energies) if (e - low) * 2625.499639 < 12.0]
        marker.write_text(json.dumps({"command": command, "elapsed_s": elapsed,
            "frames": len(frames), "energy_hartree": energies,
            "selected_indices_within_12_kJ": selected,
            "ensemble_sha256": hashlib.sha256((directory / "crest_conformers.xyz").read_bytes()).hexdigest()},
            indent=2) + "\n")
        print(f"I043 finished {name}: {len(frames)} conformers; {len(selected)} within 12 kJ/mol; {elapsed:.1f} s", flush=True)


if __name__ == "__main__":
    main()
