"""Replay producing CLIs and independently check I043 outputs against the report.

The replay uses retained candidate pools and completed electronic/quadrature
jobs; it does not claim a fresh stochastic search gives an identical ensemble.
The pinned baseline and an ethane DFT optimization/Hessian are recalculated.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import subprocess

import numpy as np
from rdkit import Chem

from run_ensembles import SPECIES, SCRATCH, allowed_scratch, read_xyz
from summarize import audit, render, REPORT, START, END, stereo_check
from thermochemistry import Ensemble, R, TEMPERATURES, CYCLE_COEFFICIENTS
from rotors import KB,HPLANCK,AMU_KG,NA,quadrature_filename,HARD_CORE_A
from scipy.linalg import eigvalsh,null_space
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
RMG_PYTHON = "/home/alon/anaconda3/envs/rmg_env/bin/python"
QUANTUM_PYTHON = "/home/alon/anaconda3/envs/pyscf_env/bin/python"


def close(actual, expected, absolute=1e-5):
    if not math.isclose(actual, expected, rel_tol=1e-10, abs_tol=absolute):
        raise AssertionError(f"numeric mismatch: {actual} != {expected}")


def verify_proposal(record,name,kind,indices=None):
    """Independently reconstruct draws, normalized densities and basin masks."""
    periods=np.array(record["periods_rad"])
    centers=np.array(record["basin_centers_rad"])
    position=(record.get("unique_minima_indices") or indices).index(record["index"])
    curvature=np.array(record["modes"]["curvature_J_mol_rad2"])
    deviation=np.minimum(np.sqrt(R*550/np.diag(curvature)),periods/2)
    regression=np.zeros_like(curvature)
    conditional=deviation
    fractions=np.array([{3:0.8,2:0.7,1:0.3}[r["symmetry"]] for r in record["definitions"]])
    if kind=="correlated":
        covariance=R*550*np.linalg.inv(curvature)
        scale=np.minimum(periods/(2*np.sqrt(np.diag(covariance))),1)
        covariance*=scale[:,None]*scale[None,:]
        conditional=np.sqrt(np.diag(covariance)).copy()
        for i in range(1,len(periods)):
            regression[i,:i]=np.linalg.solve(covariance[:i,:i],covariance[:i,i])
            conditional[i]=math.sqrt(covariance[i,i]-regression[i,:i]@covariance[:i,i])
        for key,value in (("regression",regression),("conditional_deviation",conditional),
                          ("covariance",covariance),("local_fractions",fractions)):
            if not np.allclose(value,record["importance_parameters"][key],rtol=1e-10,atol=1e-10):
                raise AssertionError("importance conditional parameters disagree")
    if record.get("importance_proposal","diagonal")!=kind:
        raise AssertionError("quadrature proposal provenance mismatch")
    seed=43003+record["index"]+int(hashlib.sha256(name.encode()).hexdigest()[:8],16)%1000000
    if kind=="correlated":seed+=2000000
    if record["seed"]!=seed:raise AssertionError("quadrature seed mismatch")
    rng=np.random.default_rng(seed)
    def periodic(x):return (x+periods/2)%periods-periods/2
    points=[]
    for _ in range(record["samples"]):
        choice=rng.random()
        if kind=="diagonal":
            point=(periodic(centers[position]+rng.normal(size=len(periods))*deviation)
                   if choice<0.8 else rng.uniform(-periods/2,periods/2))
        elif choice<0.5:
            noise=rng.normal(size=len(periods));delta=np.zeros_like(periods)
            for i in range(len(periods)):
                delta[i]=(regression[i]@delta+noise[i]*conditional[i]+periods[i]/2)%periods[i]-periods[i]/2
            point=periodic(centers[position]+delta)
        elif choice<0.9:
            local=rng.normal(size=len(periods))*deviation
            broad=rng.uniform(-periods/2,periods/2)
            delta=np.where(rng.random(len(periods))<fractions,local,broad)
            point=periodic(centers[position]+delta)
        else:point=rng.uniform(-periods/2,periods/2)
        points.append(point)
    points=np.array(points);delta=periodic(points-centers[position])
    def factors(x,sigma):
        return sum(norm.pdf(x+shift*periods,scale=sigma) for shift in range(-8,9))
    uniform=1/np.prod(periods)
    density=0.8*np.prod(factors(delta,deviation),axis=1)+0.2*uniform
    if kind=="correlated":
        residual=periodic(delta-delta@regression.T)
        density=0.5*np.prod(factors(residual,conditional),axis=1)
        density+=0.4*np.prod(fractions*factors(delta,deviation)+(1-fractions)/periods,axis=1)+0.1*uniform
    if not np.allclose(density,record["proposal_densities_rad_minus_d"],rtol=1e-9,atol=1e-15):
        raise AssertionError("importance densities differ from independent reconstruction")
    nearest=np.argmin(np.sum((periodic(points[:,None,:]-centers[None,:,:])/periods)**2,axis=2),axis=1)
    if not np.array_equal(nearest==position,record["in_basin"]):
        raise AssertionError("quadrature Voronoi mask differs")
    excluded=[i for i,(inside,value) in enumerate(zip(record["in_basin"],record["energies_above_minimum_J_mol"]))
              if inside and value is None]
    if record.get("hard_core_cutoff_A",HARD_CORE_A)!=HARD_CORE_A:
        raise AssertionError("quadrature steric cutoff differs")
    if record.get("hard_core_rejections",len(excluded))!=len(excluded):
        raise AssertionError("steric rejection count differs")
    if excluded:
        from rotors import angles,rotated
        from scipy.spatial.distance import pdist
        molecule=Chem.AddHs(Chem.MolFromSmiles(SPECIES[name]));conformer=Chem.Conformer(molecule.GetNumAtoms())
        for i,point in enumerate(record["coordinates_A"]):conformer.SetAtomPosition(i,point)
        conformer.Set3D(True);molecule.AddConformer(conformer)
        initial=angles(molecule,record["definitions"])
        for i in excluded:
            coordinates=rotated(molecule,record["definitions"],initial+periodic(points[i]-centers[position]))
            if np.min(pdist(coordinates))>=HARD_CORE_A+1e-9:
                raise AssertionError("an accepted point without a finite energy is not a steric cutoff point")
    for period,sigma in zip(periods,conditional):
        # Each conditional factor integrates to one for any conditioning mean.
        for mean in (0.,0.37*period):
            shifts=np.arange(-8,9)*period
            integral=float(np.sum(norm.cdf((period/2-mean+shifts)/sigma)-
                                  norm.cdf((-period/2-mean+shifts)/sigma)))
            close(integral,1.,1e-8)
    return points


def verify_frequency_units(source,smiles,minimum):
    """Rebuild the full Cartesian spectrum and compare with the xTB CLI."""
    from rotors import make_molecule,WAVENUMBER
    molecule=make_molecule(smiles,source/"xtbopt.xyz");n=molecule.GetNumAtoms()
    masses=np.array([a.GetMass() for a in molecule.GetAtoms()]);root=np.repeat(np.sqrt(masses),3)
    h=np.array([float(v) for line in (source/"hessian").read_text().splitlines()
                if not line.startswith("$") for v in line.split()]).reshape(3*n,3*n)
    xyz=molecule.GetConformer().GetPositions();xyz-=np.average(xyz,axis=0,weights=masses)
    rigid=[]
    for axis in np.eye(3):
        rigid.append((np.tile(axis,(n,1))*np.sqrt(masses[:,None])).ravel())
        rigid.append((np.cross(xyz,axis)*np.sqrt(masses[:,None])).ravel())
    complement=null_space(np.array(rigid))
    eigenvalues=np.linalg.eigvalsh(complement.T@(h/root[:,None]/root[None,:])@complement)
    if min(eigenvalues)<=0:raise AssertionError("full projected Cartesian Hessian is unstable")
    spectrum=np.sqrt(eigenvalues)*WAVENUMBER
    reference=np.sort(minimum["frequencies_cm1"])[-(3*n-6):]
    # CLI and RDKit standard atomic masses differ slightly; observed maximum
    # differences in ethane/DPP/trimer calibrations are below 0.1 cm^-1.
    if len(reference)!=3*n-6 or np.max(np.abs(spectrum-reference))>0.5:
        raise AssertionError("Cartesian Hessian spectrum disagrees with xTB CLI units/masses")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch", type=Path, default=SCRATCH)
    parser.add_argument("--replay", action="store_true")
    parser.add_argument("--partial", action="store_true", help="development check; not a complete dispatch verifier")
    parser.add_argument("--samples", type=int, default=4096)
    args = parser.parse_args()
    scratch = allowed_scratch(args.scratch)
    sample_counts = dict.fromkeys(SPECIES,args.samples)
    proposals=dict.fromkeys(SPECIES,"diagonal")
    basin_counts={name:{} for name in SPECIES}
    if not args.partial:
        saved = json.loads((scratch / "thermochemistry.json").read_text())
        if saved["samples_per_basin"] != args.samples:
            raise ValueError("verifier base count differs from reported thermochemistry")
        sample_counts.update(saved["samples_per_species"])
        proposals.update(saved["importance_proposals"])
        basin_counts={n:{int(i):c for i,c in values.items()} for n,values in saved["samples_per_basin_overrides"].items()}
    physical = {}
    for cpu in sorted(os.sched_getaffinity(0)):
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        core = (topology / "physical_package_id").read_text(), (topology / "core_id").read_text()
        physical.setdefault(core, cpu)
    os.sched_setaffinity(0, sorted(physical.values())[:8])
    verification = scratch / "verification"
    verification.mkdir(exist_ok=True)

    def run(python, name, arguments):
        directory = verification / name.replace(".py", "")
        directory.mkdir(exist_ok=True)
        command = [python, str(HERE / name), *arguments]
        print("I043 replay:",shlex.join(command),flush=True)
        shell = shlex.join(command) + " > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)"
        subprocess.run(["bash", "-o", "pipefail", "-c", shell], cwd=directory,
                       env=dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1"), check=True)

    if args.replay:
        if args.partial:
            raise ValueError("partial checks cannot be called a full replay")
        run(RMG_PYTHON, "run_ensembles.py", ["--scratch", str(scratch), "--prepare-only"])
        run(RMG_PYTHON, "run_ensembles.py", ["--scratch", str(scratch)])
        run(RMG_PYTHON, "baseline.py", ["--scratch", str(verification / "baseline")])
        run(RMG_PYTHON, "check_xtb.py", ["--scratch", str(scratch), "--species", *SPECIES])
        run(RMG_PYTHON, "run_composite.py", ["--scratch", str(scratch), "--hours", "2"])
        run(RMG_PYTHON, "complete_symmetry.py", ["--scratch", str(scratch)])
        for name in SPECIES:
            indices=json.loads((scratch/"composite"/name/"selection.json").read_text())["unique_minima_indices"]
            counts={i:basin_counts[name].get(i,sample_counts[name]) for i in indices}
            for count in sorted(set(counts.values())):
                wells=[i for i in indices if counts[i]==count]
                run(RMG_PYTHON,"rotors.py",["--scratch",str(scratch),"--samples",str(count),
                    "--proposal",proposals[name],"--species",name,"--indices",*map(str,wells)])
            if sample_counts[name]!=args.samples or basin_counts[name]:
                run(RMG_PYTHON,"rotors.py",["--scratch",str(scratch),"--samples",str(args.samples),
                    "--proposal",proposals[name],"--species",name])
        quad_path=scratch / "rotors/ethane/0000" / quadrature_filename(sample_counts['ethane'],43003,proposals['ethane'])
        original_quad=json.loads(quad_path.read_text())
        run(RMG_PYTHON, "rotors.py", ["--scratch", str(scratch), "--samples", str(sample_counts['ethane']),
            "--species", "ethane", "--proposal",proposals['ethane'],"--force"])
        fresh_quad=json.loads(quad_path.read_text())
        if not np.allclose(original_quad["energies_above_minimum_J_mol"],fresh_quad["energies_above_minimum_J_mol"],atol=0.02,rtol=0):
            raise AssertionError("fresh ethane rotor energies differ")
        for path in sorted((scratch/"hard_core_probe").glob("*.json")):
            probe=json.loads(path.read_text())
            run(RMG_PYTHON,"rotors.py",["--scratch",str(scratch),"--boundary-probe",
                "--samples",str(probe["samples"]),"--seed",str(probe["seed"]),"--proposal",probe["proposal"],
                "--species",probe["name"],"--indices",str(probe["index"])])
            reproduced=json.loads(path.read_text())
            for a,b in zip(probe["points"],reproduced["points"]):
                if a['sample']!=b['sample'] or a['status']!=b['status']:
                    raise AssertionError("boundary probe sample/status differs")
                if a['status']=='finite':close(a['energy_above_minimum_J_mol'],b['energy_above_minimum_J_mol'],0.05)
        thermal_args=["--scratch",str(scratch),"--samples",str(args.samples)]
        for name,count in sample_counts.items():
            if count!=args.samples:thermal_args += ["--samples-for",f"{name}={count}"]
            if proposals[name]!="diagonal":thermal_args += ["--proposal-for",f"{name}={proposals[name]}"]
            for index,npoints in basin_counts[name].items():
                thermal_args += ["--samples-in",f"{name}:{index}={npoints}"]
        run(RMG_PYTHON, "thermochemistry.py", thermal_args)
        run(QUANTUM_PYTHON, "refine.py", ["--xyz", str(scratch / "ensembles/ethane/input.xyz"),
            "--output", str(verification / "ethane_pbe.json"), "--level", "pbe"])
        original = json.loads((scratch / "pilot/ethane_pbe.json").read_text())
        fresh = json.loads((verification / "ethane_pbe.json").read_text())
        close(fresh["energy_hartree"], original["energy_hartree"], 5e-7)
        if np.max(np.abs(np.array(fresh["frequencies_cm1_real"])-original["frequencies_cm1_real"])) > 0.2:
            raise AssertionError("fresh ethane DFT frequencies differ")
        run(QUANTUM_PYTHON, "refine.py", ["--xyz", str(scratch / "xtb_checks/diphenylpropane/10000/xtbopt.xyz"),
            "--output", str(verification / "reflection_pbe.json"), "--level", "pbe", "--single-point"])
        reflected = json.loads((verification / "reflection_pbe.json").read_text())
        parent = json.loads((scratch / "composite/diphenylpropane/0000/pbe/result.json").read_text())
        close(reflected["energy_hartree"], parent["energy_hartree"], 0.02/2625.499639)
        run(QUANTUM_PYTHON, "refine.py", ["--xyz", str(scratch / "xtb_checks/triphenylheptane_iso/0000/xtbopt.xyz"),
            "--output", str(verification / "trimer_pbe.json"), "--level", "pbe", "--single-point"])
        fresh_trimer=json.loads((verification / "trimer_pbe.json").read_text())
        parent_trimer=json.loads((scratch / "composite/triphenylheptane_iso/0000/pbe/result.json").read_text())
        close(fresh_trimer["energy_hartree"],parent_trimer["energy_hartree"],0.02/2625.499639)
        run(RMG_PYTHON, "summarize.py", ["--scratch", str(scratch)])

    data = audit(scratch)
    text = REPORT.read_text().split(START, 1)[1].split(END, 1)[0].strip()
    if text != render(data):
        raise AssertionError("measured report block differs")
    baseline = json.loads((scratch / "baseline/baseline.json").read_text())
    from scipy.spatial.distance import pdist
    from rotors import EH_J
    for probe in data['boundary_probes'].values():
        close(probe['hard_core_cutoff_A'],HARD_CORE_A)
        for point in probe['points']:
            close(float(min(pdist(np.array(point['coordinates_A'])))),point['minimum_distance_A'],1e-9)
            if (point['minimum_distance_A']<HARD_CORE_A)!=(point['side']=='below'):
                raise AssertionError("boundary probe side differs")
            if point['status']=='finite':
                close((point['energy_Eh']-point['reference_energy_Eh'])*EH_J*NA,point['energy_above_minimum_J_mol'],0.01)
                close(-point['energy_above_minimum_J_mol']/(R*800),point['log_Boltzmann_factor_at_800_K'])
    if args.replay:
        fresh = json.loads((verification / "baseline/baseline.json").read_text())
        if fresh["snapshot"] != baseline["snapshot"]:
            raise AssertionError("fresh pinned snapshot differs")
        for key in ("continuous_Tc_K", "compiled_grid_Tc_K"):
            close(fresh[key], baseline[key])
        for a,b in zip(fresh["propagation"],baseline["propagation"]):
            close(a["H_J_mol"],b["H_J_mol"])
            close(a["S_J_mol_K"],b["S_J_mol_K"])
        if fresh["propagation_graphs"] != baseline["propagation_graphs"]:
            raise AssertionError("fresh compiled propagation graphs differ")
        for name, original in baseline["cycles"].items():
            calculated=fresh["cycles"][name]
            if calculated["coefficients"] != original["coefficients"] or calculated["source_vector"] != original["source_vector"]:
                raise AssertionError("fresh additive cycle coefficients or source vector differ")
            if len(calculated["thermo"]) != len(original["thermo"]):
                raise AssertionError("fresh additive cycle temperature grid differs")
            for a,b in zip(calculated["thermo"],original["thermo"]):
                for key in ("T_K","H_J_mol","S_J_mol_K"):close(a[key],b[key])
        for name, original in baseline["species"].items():
            for key in ("adjacency","sigma","optical_atom_half_factors","source_weights"):
                if fresh["species"][name][key] != original[key]:
                    raise AssertionError("fresh molecular GAV or symmetry provenance differs")
    completed=[]
    ensembles={}
    for name in SPECIES:
        selection_path=scratch / "composite" / name / "selection.json"
        if not selection_path.exists():
            if args.partial:continue
            raise AssertionError(f"missing minimum selection: {name}")
        selection=json.loads(selection_path.read_text())
        for index in selection["unique_minima_indices"]:
            source=scratch / "xtb_checks" / name / f"{index:04d}"
            minimum=json.loads((source / "result.json").read_text())
            if not minimum["minimum_confirmed"] or min(minimum["frequencies_cm1"]) < -1:
                raise AssertionError("unstable candidate used as a minimum")
            verify_frequency_units(source,SPECIES[name],minimum)
            stereo_check(SPECIES[name],read_xyz(source / "xtbopt.xyz")[0])
            for level in ("pbe","blyp"):
                target=scratch / "composite" / name / f"{index:04d}" / level
                output=target / "result.json"
                if not output.exists():
                    if args.partial:continue
                    raise AssertionError(f"missing composite calculation {name} {index} {level}")
                quantum=json.loads(output.read_text())
                if quantum["input_sha256"] != hashlib.sha256((source / "xtbopt.xyz").read_bytes()).hexdigest():
                    raise AssertionError("electronic energy used an obsolete geometry")
                cursor_quantum,cursor_min,cursor_source=quantum,minimum,source
                while "derived_by_reflection_from" in cursor_quantum or "derived_by_permutation_from" in cursor_quantum:
                    reflection="derived_by_reflection_from" in cursor_quantum
                    key="derived_by_reflection_from" if reflection else "derived_by_permutation_from"
                    parent=cursor_quantum[key]
                    parent_path=scratch / "composite" / name / f"{parent:04d}" / level / "result.json"
                    parent_result=json.loads(parent_path.read_text())
                    close(cursor_quantum["energy_hartree"],parent_result["energy_hartree"])
                    if cursor_quantum["parent_result_sha256"] != hashlib.sha256(parent_path.read_bytes()).hexdigest():
                        raise AssertionError("obsolete symmetry-parent electronic energy")
                    parent_min=scratch / "xtb_checks" / name / f"{parent:04d}"
                    coordinates=np.array([[float(v) for v in line.split()[1:]] for line in read_xyz(parent_min / "xtbopt.xyz")[0]["atoms"]])
                    if reflection:coordinates[:,0]*=-1
                    mapping=cursor_min["reflection_atom_permutation" if reflection else "permutation_atom_order"]
                    coordinates=coordinates[mapping]
                    actual=np.array([[float(v) for v in line.split()[1:]] for line in read_xyz(cursor_source / "xtbopt.xyz")[0]["atoms"]])
                    if not np.allclose(actual,coordinates,rtol=0,atol=1e-6):
                        raise AssertionError("invalid symmetry-image geometry")
                    graph=Chem.AddHs(Chem.MolFromSmiles(SPECIES[name]))
                    before=np.array([[float(v) for v in line.split()[1:]] for line in read_xyz(parent_min/"xtbopt.xyz")[0]["atoms"]])
                    for atom in graph.GetAtoms():
                        if atom.GetDegree()!=4:continue
                        order=sorted(a.GetIdx() for a in atom.GetNeighbors())
                        a,b=before[order],actual[order]
                        if np.linalg.det(a[:3]-a[3])*np.linalg.det(b[:3]-b[3])<=0:
                            raise AssertionError("symmetry image changed the labeled tetrahedral orientation")
                    if cursor_min["parent_hessian_sha256"] != hashlib.sha256((parent_min/"hessian").read_bytes()).hexdigest():
                        raise AssertionError("obsolete symmetry-parent Hessian")
                    def read_hessian(path):
                        values=[float(v) for line in path.read_text().splitlines() if not line.startswith("$") for v in line.split()]
                        return np.array(values).reshape((3*len(coordinates),)*2)
                    permutation=[3*i+a for i in mapping for a in range(3)]
                    transformed=read_hessian(parent_min/"hessian")[np.ix_(permutation,permutation)]
                    if reflection:
                        signs=np.tile([-1,1,1],len(coordinates))
                        transformed*=signs[:,None]*signs[None,:]
                    if not np.allclose(transformed,read_hessian(cursor_source/"hessian"),rtol=0,atol=1e-12):
                        raise AssertionError("invalid transformed symmetry-image Hessian")
                    target=parent_path.parent
                    cursor_quantum=parent_result
                    cursor_source=parent_min
                    cursor_min=json.loads((parent_min/"result.json").read_text())
                values=re.findall(r"converged SCF energy =\s*([-+\d.Ee]+)",(target / "stdout.log").read_text())
                if not values:raise AssertionError("SCF convergence output missing")
                close(float(values[-1]),quantum["energy_hartree"],5e-10)
        if all((scratch / "rotors" / name / f"{i:04d}" / quadrature_filename(basin_counts[name].get(i,sample_counts[name]),43003,proposals[name])).exists()
               and all((scratch / "composite" / name / f"{i:04d}" / level / "result.json").exists() for level in ("pbe","blyp"))
               for i in selection["unique_minima_indices"]):
            ensemble=Ensemble(scratch,name,baseline["species"][name],sample_counts[name],43003,proposals[name],basin_counts[name])
            ensembles[name]=ensemble
            from rotors import make_molecule, rotor_definitions, projected_modes, angles, wrap
            molecules=[make_molecule(SPECIES[name],scratch / "xtb_checks" / name / f"{i:04d}" / "xtbopt.xyz")
                       for i in selection["unique_minima_indices"]]
            definitions=rotor_definitions(molecules[0])
            periods=np.array([2*math.pi/r["symmetry"] for r in definitions])
            centers=[wrap(angles(m,definitions),periods).tolist() for m in molecules]
            manifest=json.loads((scratch/"species.json").read_text())
            if manifest[name]["achiral"]:
                expected=Chem.AddHs(Chem.MolFromSmiles(SPECIES[name]))
                for molecule in molecules:
                    xyz=np.array(molecule.GetConformer().GetPositions());reflected=xyz.copy();reflected[:,0]*=-1
                    image=Chem.Mol(molecule)
                    for i,point in enumerate(reflected):image.GetConformer().SetAtomPosition(i,point)
                    Chem.AssignStereochemistryFrom3D(image,replaceExistingTags=True)
                    tetra=[sorted(a.GetIdx() for a in t.GetNeighbors()) for t in molecule.GetAtoms() if t.GetDegree()==4]
                    def volume(points,order):
                        p=points[order]
                        return np.linalg.det(p[:3]-p[3])
                    volumes=[volume(xyz,order) for order in tetra]
                    matches=image.GetSubstructMatches(expected,useChirality=True,uniquify=False,maxMatches=1000000)
                    mapping=next((p for p in matches if all(v*volume(reflected,[p[i] for i in order])>0 for v,order in zip(volumes,tetra))),None)
                    if mapping is None:raise AssertionError("achiral reflection could not preserve tetrahedral labels")
                    for i,point in enumerate(reflected[list(mapping)]):image.GetConformer().SetAtomPosition(i,point)
                    point=wrap(angles(image,definitions),periods)
                    distance=np.min(np.max(np.abs(wrap(point-np.array(centers),periods)),axis=1))
                    if distance>0.10:raise AssertionError("reflection-related rotor well omitted")
            if ensemble.symmetry["external_divisor"]>1:
                for molecule in molecules:
                    xyz=np.array(molecule.GetConformer().GetPositions())
                    tetra=[sorted(a.GetIdx() for a in t.GetNeighbors()) for t in molecule.GetAtoms() if t.GetDegree()==4]
                    def vol(order):
                        points=xyz[order]
                        return np.linalg.det(points[:3]-points[3])
                    volumes=[vol(order) for order in tetra]
                    for perm in molecule.GetSubstructMatches(molecule,uniquify=False,useChirality=True,maxMatches=1000000):
                        if not all(v*vol([perm[i] for i in order])>0 for v,order in zip(volumes,tetra)):continue
                        image=Chem.Mol(molecule)
                        for i,point in enumerate(xyz[list(perm)]):image.GetConformer().SetAtomPosition(i,point)
                        point=wrap(angles(image,definitions),periods)
                        distance=np.min(np.max(np.abs(wrap(point-np.array(centers),periods)),axis=1))
                        if distance>0.050001:raise AssertionError("proper-permutation rotor well omitted")
            for i,molecule,record in zip(selection["unique_minima_indices"],molecules,ensemble.records):
                source=scratch / "xtb_checks" / name / f"{i:04d}"
                if record["input_sha256"] != hashlib.sha256((source / "xtbopt.xyz").read_bytes()).hexdigest():
                    raise AssertionError("obsolete geometry in rotor integral")
                if not np.array_equal(record["coordinates_A"],molecule.GetConformer().GetPositions()):
                    raise AssertionError("rotor coordinates do not match their source geometry")
                if record["basin_centers_rad"] != centers:
                    raise AssertionError("obsolete Voronoi partition in rotor integral")
                verify_proposal(record,name,proposals[name],ensemble.indices)
                numbers=[float(v) for line in (source / "hessian").read_text().splitlines()
                         if not line.startswith("$") for v in line.split()]
                modes=projected_modes(molecule,definitions,np.array(numbers).reshape((3*molecule.GetNumAtoms(),)*2))
                if len(modes["non_torsional_cm1"]) != 3*molecule.GetNumAtoms()-6-len(definitions) or len(modes["torsional_cm1"]) != len(definitions):
                    raise AssertionError("torsional and vibrational degrees of freedom were counted incorrectly")
                for key in ("non_torsional_cm1","torsional_cm1","inertia_amu_A2"):
                    if not np.allclose(modes[key],record["modes"][key],rtol=1e-10,atol=1e-8):
                        raise AssertionError("obsolete Hessian in rotor integral")
                inertia=np.array(record["modes"]["inertia_amu_A2"])*AMU_KG*1e-20
                curvature=np.array(record["modes"]["curvature_J_mol_rad2"])/NA
                frequencies=np.sqrt(eigvalsh(curvature,inertia))/(2*math.pi)
                d=len(definitions)
                for t in (298.15,600.,800.):
                    x=HPLANCK*frequencies/(KB*t)
                    log_classical=d*math.log(2*math.pi*KB*t/HPLANCK)
                    log_classical+=(np.linalg.slogdet(inertia)[1]-np.linalg.slogdet(curvature)[1])/2
                    close(log_classical,float(-np.log(x).sum()),1e-9)
                    log_pg=float((np.log(x)-x/2-np.log(-np.expm1(-x))).sum())
                    quantum_harmonic=float((-x/2-np.log(-np.expm1(-x))).sum())
                    close(log_classical+log_pg,quantum_harmonic,1e-9)
            for level in ("pbe","blyp"):
                for t in (298.15,600.,800.):
                    row=ensemble.hs(t,level)
                    left=ensemble.hs(t-0.02,level)
                    right=ensemble.hs(t+0.02,level)
                    g_left=left["H_J_mol"]-(t-0.02)*left["S_J_mol_K"]
                    g_right=right["H_J_mol"]-(t+0.02)*right["S_J_mol_K"]
                    close(-(g_right-g_left)/0.04,row["S_J_mol_K"],0.004)
            completed.append(name)
    if not args.partial and len(completed)!=len(SPECIES):
        raise AssertionError("thermochemistry inputs incomplete")
    if not args.partial:
        thermo=json.loads((scratch / "thermochemistry.json").read_text())
        weights={"dyads":{"diphenylpentane_meso":0.5,"diphenylpentane_racemo":0.5},
                 "triads":{"triphenylheptane_iso":0.25,"triphenylheptane_syndio":0.25,"triphenylheptane_hetero":0.5}}
        external={side:-R*sum(w*math.log(ensembles[n].symmetry["external_divisor"]) for n,w in values.items())
                  for side,values in weights.items()}
        mixing={side:-R*sum(w*math.log(w) for w in values.values()) for side,values in weights.items()}
        stereo=thermo["stereo_accounting"]
        close(external["triads"]-external["dyads"],stereo["external_end_symmetry_increment_J_mol_K"])
        close(mixing["triads"]-mixing["dyads"],stereo["excluded_class_mixing_entropy_increment_J_mol_K"])
        close(stereo["end_plus_class_plus_mirror_increment_J_mol_K"],R*math.log(2))
        gav_end = {side:-R*sum(w*math.log(baseline["species"][n]["sigma"]*
            2**baseline["species"][n]["optical_atom_half_factors"]/
            ensembles[n].symmetry["internal_rotor_divisor"]) for n,w in values.items())
            for side,values in weights.items()}
        gav_opt = {side:R*math.log(2)*sum(w*baseline["species"][n]["optical_atom_half_factors"]
            for n,w in values.items()) for side,values in weights.items()}
        close(gav_end["triads"]-gav_end["dyads"],stereo["GAV_external_end_symmetry_increment_J_mol_K"])
        close(gav_opt["triads"]-gav_opt["dyads"],stereo["GAV_local_optical_increment_J_mol_K"])
        entropy_offset=(gav_end["triads"]-gav_end["dyads"]+gav_opt["triads"]-gav_opt["dyads"]
                        -external["triads"]+external["dyads"])
        close(entropy_offset,stereo["matched_entropy_offset_J_mol_K"])
        close(entropy_offset,mixing["triads"]-mixing["dyads"])
        for name,ensemble in ensembles.items():
            count=sum(sum(inside and value is None for inside,value in zip(record['in_basin'],record['energies_above_minimum_J_mol'])) for record in ensemble.records)
            close(count,thermo['species'][name]['hard_core_zero_weight_points'],0)
            close(sum(record['samples'] for record in ensemble.records),thermo['species'][name]['total_quadrature_points'],0)
            for level in ("pbe","blyp"):
                for t,record in zip(TEMPERATURES,thermo["species"][name]["levels"][level]):
                    calculated=ensemble.hs(t,level)
                    for key in ("H_J_mol","S_J_mol_K","H_above_electronic_minimum_kJ_mol",
                                "S_conformer_J_mol_K","minimum_effective_samples"):
                        close(calculated[key],record[key])
                    for key,value in calculated["entropy_components"].items():close(value,record["entropy_components"][key])
                    close(sum(record["entropy_components"].values())+record["S_conformer_J_mol_K"],record["S_J_mol_K"])
                    if not np.allclose(calculated["MC_standard_errors_H_S_G"],record["MC_standard_errors_H_S_G"],rtol=1e-10,atol=1e-8):
                        raise AssertionError("molecular quadrature errors differ")
        for name,coefficients in CYCLE_COEFFICIENTS.items():
            # Independently check atom balance and sum tabulated molecular H/S.
            atom_balance={}
            for species,coefficient in coefficients.items():
                molecule=Chem.AddHs(Chem.MolFromSmiles(SPECIES[species]))
                for atom in molecule.GetAtoms():
                    key=atom.GetSymbol()
                    atom_balance[key]=atom_balance.get(key,0)+coefficient
            if any(abs(v)>1e-10 for v in atom_balance.values()):
                raise AssertionError("unbalanced diagnostic cycle")
            for level in ("pbe","blyp"):
                rows=thermo["cycles"][name]["levels"][level]
                for i,row in enumerate(rows):
                    for key in ("H_J_mol","S_J_mol_K"):
                        value=sum(c*thermo["species"][s]["levels"][level][i][key]
                                  for s,c in coefficients.items())
                        close(value,row[key])
                    close(sum(row["entropy_components"].values()),row["S_J_mol_K"])
                    translation=1.5*R*sum(c*math.log(sum(r["modes"]["masses_amu"]))
                        for species,c in coefficients.items() for r in ensembles[species].records[:1])
                    close(translation,row["entropy_components"]["S_translation_J_mol_K"])
                correction=thermo["conditional_corrections"][name][level]
                h,s=rows[0]["H_J_mol"],rows[0]["S_J_mol_K"]
                if name=="diphenylpropane_isodesmic":
                    model=baseline["cycles"][name]["thermo"][0]
                    h,s=-(h-model["H_J_mol"]),-(s-model["S_J_mol_K"])
                else:
                    model=baseline["cycles"][name]["thermo"][0]
                    h,s=h-model["H_J_mol"],s-model["S_J_mol_K"]+entropy_offset
                close(h/1000,correction["delta_H298_kJ_mol"])
                close(s,correction["delta_S298_J_mol_K"])
                close((baseline["propagation"][0]["H_J_mol"]+h)/1000,correction["propagation_H298_kJ_mol"])
                close(baseline["propagation"][0]["S_J_mol_K"]+s,correction["propagation_S298_J_mol_K"])
        if "quadrature_comparison" in thermo:
            comparison=thermo["quadrature_comparison"]
            alternative=Ensemble(scratch,"diphenylpropane",baseline["species"]["diphenylpropane"],comparison["points"][0],43003)
            for low,high in zip(alternative.records,ensembles["diphenylpropane"].records):
                if low["input_sha256"]!=high["input_sha256"] or low["basin_centers_rad"]!=high["basin_centers_rad"]:
                    raise AssertionError("quadrature comparison uses different coordinate partitions")
                a,b=low["energies_above_minimum_J_mol"],high["energies_above_minimum_J_mol"][:low["samples"]]
                if [v is None for v in a]!=[v is None for v in b] or any(abs(x-y)>0.05 for x,y in zip(a,b) if x is not None):
                    raise AssertionError("same-seed quadrature prefixes differ")
            for level,rows in comparison["levels"].items():
                for row in rows:
                    data={name:(alternative if name=="diphenylpropane" else ensemble).hs(row["T_K"],level)
                          for name,ensemble in ensembles.items()}
                    for key in ("H_J_mol","S_J_mol_K"):
                        value=sum(c*data[name][key] for name,c in CYCLE_COEFFICIENTS["diphenylpropane_isodesmic"].items())
                        close(value,row[key])
            print("I043 same-partition quadrature comparison and deterministic energy prefixes verified")
        if "precision_comparison" in thermo:
            comparison=thermo["precision_comparison"];base=dict(ensembles)
            for name in comparison["increased_species"]:
                base[name]=Ensemble(scratch,name,baseline["species"][name],comparison["base_points"],43003,proposals[name])
                for low,high in zip(base[name].records,ensembles[name].records):
                    if low["input_sha256"]!=high["input_sha256"] or low["basin_centers_rad"]!=high["basin_centers_rad"]:
                        raise AssertionError("precision comparison changed physical partitions")
                    a,b=low["energies_above_minimum_J_mol"],high["energies_above_minimum_J_mol"][:low["samples"]]
                    if [v is None for v in a]!=[v is None for v in b] or any(abs(x-y)>0.05 for x,y in zip(a,b) if x is not None):
                        raise AssertionError("precision comparison random stream prefixes differ")
            for name,levels in comparison["cycles"].items():
                for level,rows in levels.items():
                    for row in rows:
                        data={n:e.hs(row["T_K"],level) for n,e in base.items()}
                        for key in ("H_J_mol","S_J_mol_K"):
                            close(sum(c*data[n][key] for n,c in CYCLE_COEFFICIENTS[name].items()),row[key])
            print("I043 all increased-well allocations and same-proposal convergence reproduced")
        if not thermo["numerical_preflight"]["passed"]:
            raise AssertionError("reported thermochemistry failed numerical precision preflight")
        for data in thermo["cycles"].values():
            for rows in data["levels"].values():
                for row in rows:
                    h,s,g=row["MC_standard_errors_H_S_G"]
                    if h>500 or s>1 or g>500:
                        raise AssertionError("reported cycle misses a predeclared precision target")
        from rmgpy.data.rmg import RMGDatabase
        from rmgpy.molecule import Molecule
        from rmgpy.reaction import Reaction
        from rmgpy.species import Species
        db=RMGDatabase()
        db.load_thermo(str(scratch / "baseline/database/input/thermo"),thermo_libraries=["primaryThermoLibrary"])
        sides={side:[Species(molecule=[Molecule().from_adjacency_list(adj)])
                     for adj in baseline["propagation_graphs"][side]] for side in ("reactants","products")}
        for side in sides.values():
            for species in side:species.thermo=db.thermo.get_thermo_data(species)
        propagation=Reaction(**sides)
        for name,coefficients in CYCLE_COEFFICIENTS.items():
            for level in ("pbe","blyp"):
                for root in thermo["conditional_corrections"][name][level]["roots_298_to_800_K"]:
                    values={s:e.hs(root,level) for s,e in ensembles.items()}
                    h=sum(c*values[s]["H_J_mol"] for s,c in coefficients.items())
                    s=sum(c*values[name]["S_J_mol_K"] for name,c in coefficients.items())
                    if name=="diphenylpropane_isodesmic":
                        model=baseline["cycles"][name]["thermo"][0]
                        h,s=-(h-model["H_J_mol"]),-(s-model["S_J_mol_K"])
                    else:
                        model=baseline["cycles"][name]["thermo"][0]
                        h,s=h-model["H_J_mol"],s-model["S_J_mol_K"]+entropy_offset
                    residual=propagation.get_enthalpy_of_reaction(root)+h-root*(propagation.get_entropy_of_reaction(root)+s)-R*root*math.log(1000*R*root/100000)
                    close(residual,0,1e-4)
        print("I043 all molecular tables, balanced cycle sums, correction signs and corrected Tc residuals verified")
    print(f"I043 {'PARTIAL development check' if args.partial else 'verification'}: retained pool hashes/stereo and composite source energies verified")
    print(f"I043 H/S thermodynamic derivative identities independently verified for {len(completed)} complete ensembles")
    print("I043 rotor/vibration mode counts and the Pitzer-Gwinn harmonic reference limit independently verified")
    print("I043 full Cartesian Hessian spectra agree with xTB CLI frequencies; importance draws/densities/basin masks independently reconstructed")
    if args.replay:
        print("I043 every producing CLI replayed; pinned baseline, ethane DFT optimization/Hessian, ethane rotor quadrature, a DPP reflection energy and a trimer energy freshly recalculated")
        print("I043 retained-pool report reproduction passed; fresh stochastic ensemble identity is not asserted")


if __name__ == "__main__":
    main()
