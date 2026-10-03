"""Coupled torsional quadrature on disjoint conformer basins, with GFN2-xTB.

Commands and approximation limits belong in the report. Cartesian Hessian
projection removes the torsional harmonic subspace before integrating it.
Only one Voronoi basin accepts a torsion point, so wells are not counted twice.
"""
from __future__ import annotations

import argparse
import ctypes as ct
import hashlib
import json
import math
import os
from pathlib import Path
import time

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from scipy.linalg import null_space
from scipy.spatial.distance import pdist

from run_ensembles import SPECIES, SCRATCH, allowed_scratch, read_xyz

BOHR_A = 0.529177210903
EH_J = 4.3597447222071e-18
AMU_KG = 1.66053906660e-27
HBAR = 1.0545718176461565e-34
HPLANCK = 2*math.pi*HBAR
KB = 1.380649e-23
R = 8.314472
NA = R/KB  # use the dispatch's R consistently, including the energy conversion
WAVENUMBER = math.sqrt(EH_J/(AMU_KG*(BOHR_A*1e-10)**2))/(2*math.pi*29979245800.)
HARD_CORE_A = 0.45


class XTB:
    """Small wrapper of the installed, filesystem-documented xTB C API."""
    def __init__(self, molecule, coordinates):
        self.template = Chem.Mol(molecule)
        self.retries = 0
        self.library = ct.CDLL("/home/alon/anaconda3/envs/crest_env/lib/libxtb.so")
        lib = self.library
        pointer = ct.c_void_p
        specifications = {
            "xtb_newEnvironment": (pointer, []),
            "xtb_newCalculator": (pointer, []),
            "xtb_newResults": (pointer, []),
            "xtb_setVerbosity": (None, [pointer, ct.c_int]),
            "xtb_newMolecule": (pointer, [pointer, pointer, pointer, pointer, pointer, pointer, pointer, pointer]),
            "xtb_loadGFN2xTB": (None, [pointer, pointer, pointer, ct.c_char_p]),
            "xtb_updateMolecule": (None, [pointer, pointer, pointer, pointer]),
            "xtb_singlepoint": (None, [pointer, pointer, pointer, pointer]),
            "xtb_getEnergy": (None, [pointer, pointer, pointer]),
            "xtb_checkEnvironment": (ct.c_int, [pointer]),
            "xtb_getError": (None, [pointer, pointer, pointer]),
            "xtb_setMaxIter": (None, [pointer, pointer, ct.c_int]),
            "xtb_setAccuracy": (None, [pointer, pointer, ct.c_double]),
            "xtb_delEnvironment": (None, [pointer]),
            "xtb_delMolecule": (None, [pointer]),
            "xtb_delCalculator": (None, [pointer]),
            "xtb_delResults": (None, [pointer]),
        }
        for name, (result, arguments) in specifications.items():
            getattr(lib, name).restype = result
            getattr(lib, name).argtypes = arguments
        self.env = ct.c_void_p(lib.xtb_newEnvironment())
        lib.xtb_setVerbosity(self.env, 0)
        numbers = np.array([atom.GetAtomicNum() for atom in molecule.GetAtoms()], dtype=np.int32)
        coordinates = np.array(coordinates/BOHR_A, dtype=np.float64, order="C")
        n = ct.c_int(len(numbers))
        charge = ct.c_double(0)
        spin = ct.c_int(0)
        self.mol = ct.c_void_p(lib.xtb_newMolecule(self.env, ct.byref(n), numbers.ctypes.data,
            coordinates.ctypes.data, ct.byref(charge), ct.byref(spin), None, None))
        self.calc = ct.c_void_p(lib.xtb_newCalculator())
        self.res = ct.c_void_p(lib.xtb_newResults())
        lib.xtb_loadGFN2xTB(self.env, self.mol, self.calc, None)
        lib.xtb_setMaxIter(self.env, self.calc, 500)
        lib.xtb_setAccuracy(self.env, self.calc, 0.1)
        self.check()

    def check(self):
        if self.library.xtb_checkEnvironment(self.env):
            message = ct.create_string_buffer(2048)
            size = ct.c_int(2048)
            self.library.xtb_getError(self.env, message, ct.byref(size))
            raise RuntimeError(message.value.decode())

    def singlepoint(self, coordinates):
        coordinates = np.array(coordinates/BOHR_A, dtype=np.float64, order="C")
        self.library.xtb_updateMolecule(self.env, self.mol, coordinates.ctypes.data, None)
        self.library.xtb_singlepoint(self.env, self.mol, self.calc, self.res)
        self.check()
        energy = ct.c_double()
        self.library.xtb_getEnergy(self.env, self.res, ct.byref(energy))
        self.check()
        if not math.isfinite(energy.value):
            raise RuntimeError("xTB returned a non-finite electronic energy")
        return energy.value

    def energy(self, coordinates,apply_hard_core=True):
        if apply_hard_core and np.min(pdist(coordinates)) < HARD_CORE_A:
            return math.inf
        try:
            return self.singlepoint(coordinates)
        except RuntimeError as error:
            if "iterator did not converge" not in str(error):
                raise
        template = self.template
        attempts = self.retries
        for accuracy in (0.1,1.0):
            self.close()
            self.__init__(template,coordinates)
            attempts += 1
            self.retries = attempts
            self.library.xtb_setMaxIter(self.env,self.calc,2000)
            self.library.xtb_setAccuracy(self.env,self.calc,accuracy)
            try:
                energy = self.singlepoint(coordinates)
                if accuracy != 0.1:
                    self.library.xtb_setAccuracy(self.env,self.calc,0.1)
                    energy = self.singlepoint(coordinates)
                return energy
            except RuntimeError as error:
                if "iterator did not converge" not in str(error):
                    raise
        raise RuntimeError("xTB SCC still failed after fresh electronic restarts")

    def close(self):
        for name, handle in (("Results", self.res), ("Calculator", self.calc),
                             ("Molecule", self.mol), ("Environment", self.env)):
            getattr(self.library, "xtb_del" + name)(ct.byref(handle))


def make_molecule(smiles, xyz):
    frame = read_xyz(xyz)[0]
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    conformer = Chem.Conformer(mol.GetNumAtoms())
    for i, line in enumerate(frame["atoms"]):
        conformer.SetAtomPosition(i, [float(v) for v in line.split()[1:]])
    conformer.Set3D(True)
    mol.AddConformer(conformer)
    return mol


def rotor_definitions(mol):
    rotors = []
    for bond in mol.GetBonds():
        if bond.IsInRing() or bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        j, k = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if mol.GetAtomWithIdx(j).GetAtomicNum() == 1 or mol.GetAtomWithIdx(k).GetAtomicNum() == 1:
            continue
        left = sorted([a.GetIdx() for a in mol.GetAtomWithIdx(j).GetNeighbors() if a.GetIdx() != k],
                      key=lambda i: (-mol.GetAtomWithIdx(i).GetAtomicNum(), i))
        right = sorted([a.GetIdx() for a in mol.GetAtomWithIdx(k).GetNeighbors() if a.GetIdx() != j],
                       key=lambda i: (-mol.GetAtomWithIdx(i).GetAtomicNum(), i))
        if not left or not right:
            continue
        ends = [mol.GetAtomWithIdx(j), mol.GetAtomWithIdx(k)]
        methyl = any(sum(a.GetAtomicNum() == 1 for a in end.GetNeighbors()) == 3 for end in ends)
        phenyl = any(end.GetIsAromatic() for end in ends)
        symmetry = 3 if methyl else 2 if phenyl else 1
        rotors.append({"atoms": [left[0], j, k, right[0]], "symmetry": symmetry,
                       "kind": "methyl" if methyl else "phenyl" if phenyl else "backbone"})
    return rotors


def wrap(angle, periods):
    return (angle + periods/2) % periods - periods/2


def angles(mol, definitions):
    conf = mol.GetConformer()
    return np.array([rdMolTransforms.GetDihedralRad(conf, *r["atoms"]) for r in definitions])


def rotated(mol, definitions, values):
    clone = Chem.Mol(mol)
    conf = clone.GetConformer()
    for rotor, value in zip(definitions, values):
        rdMolTransforms.SetDihedralRad(conf, *rotor["atoms"], float(value))
    return np.array(conf.GetPositions())


def projected_modes(mol, definitions, hessian):
    coordinates = np.array(mol.GetConformer().GetPositions())
    masses = np.array([atom.GetMass() for atom in mol.GetAtoms()])
    root_mass = np.repeat(np.sqrt(masses), 3)
    mw = hessian/root_mass[:, None]/root_mass[None, :]
    centered = coordinates-np.average(coordinates, axis=0, weights=masses)
    rigid = []
    for axis in np.eye(3):
        rigid.append((np.tile(axis, (len(masses), 1))*np.sqrt(masses[:, None])).ravel())
        rigid.append((np.cross(centered, axis)*np.sqrt(masses[:, None])).ravel())
    rigid = np.linalg.qr(np.array(rigid).T)[0]
    initial = angles(mol, definitions)
    tangents = []
    for index in range(len(definitions)):
        shifted = initial.copy()
        shifted[index] += 1e-4
        derivative = (rotated(mol, definitions, shifted)-coordinates)/1e-4
        tangent = derivative.ravel()*root_mass
        tangent -= rigid@(rigid.T@tangent)
        tangents.append(tangent)
    tangents = np.array(tangents).T
    inertia = tangents.T@tangents
    if np.min(np.linalg.eigvalsh(inertia)) <= 1e-8:
        raise AssertionError("dependent torsional coordinates")
    basis = np.linalg.qr(tangents)[0]
    complement = null_space(np.column_stack([rigid, basis]).T)
    hv = complement.T@mw@complement
    coupling = basis.T@mw@complement
    ht = basis.T@mw@basis - coupling@np.linalg.solve(hv, coupling.T)
    non_tor = np.linalg.eigvalsh(hv)
    tor = np.linalg.eigvalsh(ht)
    if min(np.min(non_tor), np.min(tor)) <= 0:
        raise AssertionError("projected rotor/vibration Hessian is not positive definite")
    # Rigid torsion curvature controls the sampling proposal and Pitzer-Gwinn
    # reference. The Schur complement diagnoses stability; the orthogonal
    # complement supplies the non-torsional harmonic modes.
    curvature = tangents.T@mw@tangents/(BOHR_A**2)*EH_J*NA
    return {"non_torsional_cm1": (np.sqrt(non_tor)*WAVENUMBER).tolist(),
            "torsional_cm1": (np.sqrt(tor)*WAVENUMBER).tolist(),
            "inertia_amu_A2": inertia.tolist(), "curvature_J_mol_rad2": curvature.tolist(),
            "masses_amu": masses.tolist()}


def wrapped_normal_factors(displacement, deviation, periods):
    terms = np.zeros_like(displacement)
    for shift in range(-8, 9):
        terms += np.exp(-0.5*((displacement + shift*periods)/deviation)**2)/(math.sqrt(2*math.pi)*deviation)
    return terms


def wrapped_normal_pdf(displacement,deviation,periods):
    return float(np.prod(wrapped_normal_factors(displacement,deviation,periods)))


def quadrature_filename(samples,seed,proposal_kind="diagonal"):
    suffix="" if proposal_kind=="diagonal" else "_correlated"
    return f"quadrature_{samples}_{seed}{suffix}.json"


def importance_parameters(curvature,periods,symmetries,proposal_kind):
    parameters={"independent_deviation":np.minimum(np.sqrt(R*550/np.diag(curvature)),periods/2)}
    if proposal_kind=="diagonal":return parameters
    covariance=R*550*np.linalg.inv(curvature)
    scale=np.minimum(periods/(2*np.sqrt(np.diag(covariance))),1)
    covariance=covariance*scale[:,None]*scale[None,:]
    cholesky=np.linalg.cholesky((covariance+covariance.T)/2)
    regression=np.zeros_like(covariance)
    for i in range(1,len(periods)):
        regression[i,:i]=np.linalg.solve(cholesky[:i,:i].T,cholesky[i,:i])
    parameters.update(covariance=covariance,regression=regression,
                      conditional_deviation=np.diag(cholesky),
                      local_fractions=np.array([{3:0.8,2:0.7,1:0.3}[s] for s in symmetries]))
    return parameters


def sample_point(rng,center,periods,parameters,proposal_kind):
    choice=rng.random()
    if proposal_kind=="diagonal":
        if choice<0.8:return wrap(center+rng.normal(size=len(periods))*parameters["independent_deviation"],periods)
        return rng.uniform(-periods/2,periods/2)
    if choice<0.5:
        noise=rng.normal(size=len(periods));displacement=np.zeros_like(periods)
        for i in range(len(periods)):
            mean=parameters["regression"][i]@displacement
            displacement[i]=wrap(mean+noise[i]*parameters["conditional_deviation"][i],periods[i])
        return wrap(center+displacement,periods)
    if choice<0.9:
        normal=rng.normal(size=len(periods))*parameters["independent_deviation"]
        uniform=rng.uniform(-periods/2,periods/2)
        displacement=np.where(rng.random(len(periods))<parameters["local_fractions"],normal,uniform)
        return wrap(center+displacement,periods)
    return rng.uniform(-periods/2,periods/2)


def point_density(displacement,periods,parameters,proposal_kind):
    uniform=1/np.prod(periods)
    if proposal_kind=="diagonal":
        return 0.8*wrapped_normal_pdf(displacement,parameters["independent_deviation"],periods)+0.2*uniform
    residual=wrap(displacement-parameters["regression"]@displacement,periods)
    correlated=wrapped_normal_pdf(residual,parameters["conditional_deviation"],periods)
    marginal=wrapped_normal_factors(displacement,parameters["independent_deviation"],periods)
    fractions=parameters["local_fractions"]
    product=float(np.prod(fractions*marginal+(1-fractions)/periods))
    # Every wrapped conditional is normalized on its period. The triangular
    # dependency makes their product normalized by successive integration.
    return 0.5*correlated+0.4*product+0.1*uniform


def integrate(scratch, name, samples, seed, force=False,proposal_kind="diagonal",only_indices=None):
    selection_path = scratch / "composite" / name / "selection.json"
    selection = json.loads(selection_path.read_text())
    from complete_symmetry import augment
    selection = augment(scratch,name,selection)
    selection_path.write_text(json.dumps(selection,indent=2)+"\n")
    indices = selection["unique_minima_indices"]
    molecules = [make_molecule(SPECIES[name], scratch / "xtb_checks" / name / f"{i:04d}" / "xtbopt.xyz") for i in indices]
    definitions = rotor_definitions(molecules[0])
    periods = np.array([2*math.pi/r["symmetry"] for r in definitions])
    centers = np.array([wrap(angles(m, definitions), periods) for m in molecules])
    results = []
    for position, (index, mol) in enumerate(zip(indices, molecules)):
        if only_indices is not None and index not in only_indices:continue
        directory = scratch / "rotors" / name / f"{index:04d}"
        directory.mkdir(parents=True, exist_ok=True)
        target = directory / quadrature_filename(samples,seed,proposal_kind)
        if target.exists() and not force:
            cached = json.loads(target.read_text())
            current_source=scratch/"xtb_checks"/name/f"{index:04d}"/"xtbopt.xyz"
            if cached["basin_centers_rad"] == centers.tolist() and cached["input_sha256"]==hashlib.sha256(current_source.read_bytes()).hexdigest() and np.array_equal(cached["coordinates_A"],mol.GetConformer().GetPositions()):
                results.append(cached)
                continue
        source = scratch / "xtb_checks" / name / f"{index:04d}"
        if not np.array_equal(mol.GetConformer().GetPositions(),make_molecule(SPECIES[name],source/"xtbopt.xyz").GetConformer().GetPositions()):
            raise RuntimeError("input geometry changed after basin construction; restart integration")
        geometry_digest=hashlib.sha256((source/"xtbopt.xyz").read_bytes()).hexdigest()
        numbers = [float(value) for line in (source / "hessian").read_text().splitlines()
                   if not line.startswith("$") for value in line.split()]
        hessian = np.array(numbers).reshape((3*mol.GetNumAtoms(),)*2)
        if np.max(np.abs(hessian-hessian.T)) > 1e-6:
            raise AssertionError("xTB Hessian is not symmetric")
        modes = projected_modes(mol, definitions, hessian)
        curvature = np.array(modes["curvature_J_mol_rad2"])
        parameters=importance_parameters(curvature,periods,[r["symmetry"] for r in definitions],proposal_kind)
        molecule_seed = seed+index+int(hashlib.sha256(name.encode()).hexdigest()[:8], 16)%1000000+(2000000 if proposal_kind=="correlated" else 0)
        rng = np.random.default_rng(molecule_seed)
        reference = np.array(mol.GetConformer().GetPositions())
        calculator = XTB(mol, reference)
        energy0 = calculator.energy(reference)
        xtb_record = json.loads((source / "result.json").read_text())
        if abs(energy0-xtb_record["energy_Eh"])*2625.499639 > 0.05:
            raise AssertionError("xTB C API and command-line reference energies differ")
        energy, density, accepted = [], [], []
        initial = angles(mol, definitions)
        start = time.monotonic()
        for sample in range(samples):
            point=sample_point(rng,centers[position],periods,parameters,proposal_kind)
            distances = np.sum((wrap(point-centers, periods)/periods)**2, axis=1)
            belongs = int(np.argmin(distances)) == position
            displacement = wrap(point-centers[position], periods)
            proposal = point_density(displacement,periods,parameters,proposal_kind)
            if belongs:
                coordinates = rotated(mol, definitions, initial+wrap(point-centers[position], periods))
                try:
                    value = (calculator.energy(coordinates)-energy0)*EH_J*NA
                except RuntimeError:
                    xyz = str(mol.GetNumAtoms())+f"\nfailed quadrature sample {sample}\n"
                    xyz += "\n".join(atom.GetSymbol()+" "+" ".join(f"{v:.12f}" for v in point)
                                      for atom,point in zip(mol.GetAtoms(),coordinates))+"\n"
                    (directory / "failed_scf.xyz").write_text(xyz)
                    raise
            else:
                value = math.inf
            energy.append(value if math.isfinite(value) else None)
            density.append(proposal)
            accepted.append(belongs)
        retries = calculator.retries
        calculator.close()
        if geometry_digest!=hashlib.sha256((source/"xtbopt.xyz").read_bytes()).hexdigest():
            raise RuntimeError("input geometry changed during rotor integration; result not saved")
        result = {"name": name, "index": index, "seed": molecule_seed, "samples": samples,
                  "definitions": definitions, "periods_rad": periods.tolist(),
                  "basin_centers_rad": centers.tolist(), "energies_above_minimum_J_mol": energy,
                  "unique_minima_indices": indices,
                  "proposal_densities_rad_minus_d": density, "in_basin": accepted,
                  "modes": modes, "reference_energy_Eh": energy0,
                  "coordinates_A": reference.tolist(), "elapsed_s": time.monotonic()-start,
                  "input_sha256": geometry_digest}
        result["fresh_SCC_restarts"] = retries
        result["importance_proposal"]=proposal_kind
        result["hard_core_cutoff_A"]=HARD_CORE_A
        result["hard_core_rejections"]=sum(inside and value is None for inside,value in zip(accepted,energy))
        result["importance_parameters"]={key:value.tolist() for key,value in parameters.items()}
        target.write_text(json.dumps(result, indent=2) + "\n")
        print(f"I043 coupled rotors {name} {index}: {len(periods)} torsions, {samples} samples, "
              f"{result['elapsed_s']:.1f} s", flush=True)
        results.append(result)
    return results


def probe_boundary(scratch,name,index,samples,seed,proposal_kind):
    """Evaluate one retained point on either side of the steric boundary."""
    path=scratch/"rotors"/name/f"{index:04d}"/quadrature_filename(samples,seed,proposal_kind)
    record=json.loads(path.read_text())
    source=scratch/"xtb_checks"/name/f"{index:04d}"/"xtbopt.xyz"
    if record["input_sha256"]!=hashlib.sha256(source.read_bytes()).hexdigest():
        raise RuntimeError("boundary probe input is obsolete")
    molecule=make_molecule(SPECIES[name],source);definitions=record["definitions"]
    periods=np.array(record["periods_rad"]);centers=np.array(record["basin_centers_rad"])
    ids=json.loads((scratch/"composite"/name/"selection.json").read_text())["unique_minima_indices"]
    position=ids.index(index);initial=angles(molecule,definitions)
    parameters=importance_parameters(np.array(record["modes"]["curvature_J_mol_rad2"]),periods,
                                     [r["symmetry"] for r in definitions],proposal_kind)
    rng=np.random.default_rng(record["seed"]);candidates={"below":None,"above":None}
    reference=np.array(molecule.GetConformer().GetPositions())
    for sample,inside in enumerate(record["in_basin"]):
        point=sample_point(rng,centers[position],periods,parameters,proposal_kind)
        if not inside:continue
        coordinates=rotated(molecule,definitions,initial+wrap(point-centers[position],periods))
        distances=pdist(coordinates);distance=float(min(distances))
        side="below" if distance<HARD_CORE_A else "above"
        candidate=candidates[side]
        if candidate is None or abs(distance-HARD_CORE_A)<abs(candidate["minimum_distance_A"]-HARD_CORE_A):
            pairs=list(zip(*np.triu_indices(len(coordinates),1)));pair=pairs[int(np.argmin(distances))]
            candidates[side]={"sample":sample,"minimum_distance_A":distance,
                "pair_atoms":[int(i) for i in pair],"pair_symbols":[molecule.GetAtomWithIdx(int(i)).GetSymbol() for i in pair],
                "coordinates_A":coordinates.tolist()}
    rows=[]
    for side,candidate in candidates.items():
        if candidate is None:continue
        calculator=XTB(molecule,reference)
        try:
            energy0=calculator.energy(reference)
            energy=calculator.energy(np.array(candidate["coordinates_A"]),apply_hard_core=False)
            difference=(energy-energy0)*EH_J*NA
            log_factor=-difference/(R*800)
            candidate.update(side=side,status="finite",energy_above_minimum_J_mol=difference,
                             reference_energy_Eh=energy0,energy_Eh=energy,
                             log_Boltzmann_factor_at_800_K=log_factor,
                             Boltzmann_factor_at_800_K=(math.exp(log_factor) if log_factor>-745 else 0.0)
                             if log_factor<700 else None)
        except RuntimeError as error:
            candidate.update(side=side,status="SCC_failed",error=str(error))
        finally:calculator.close()
        rows.append(candidate)
        print(f"I043 boundary probe {name} {index} {side}: distance={candidate['minimum_distance_A']:.6f} A, status={candidate['status']}",flush=True)
    result={"name":name,"index":index,"samples":samples,"seed":seed,"proposal":proposal_kind,
            "source_input_sha256":record["input_sha256"],"hard_core_cutoff_A":HARD_CORE_A,"points":rows,
            "scope":"representative retained points near the boundary; not a global bound on the excluded region"}
    directory=scratch/"hard_core_probe";directory.mkdir(exist_ok=True)
    (directory/f"{name}_{index}_{proposal_kind}.json").write_text(json.dumps(result,indent=2)+"\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch", type=Path, default=SCRATCH)
    parser.add_argument("--species", choices=tuple(SPECIES), nargs="+", required=True)
    parser.add_argument("--samples", type=int, default=4096)
    parser.add_argument("--seed", type=int, default=43003)
    parser.add_argument("--proposal",choices=("diagonal","correlated"),default="diagonal")
    parser.add_argument("--indices",type=int,nargs="+",help="pilot selected wells on the unchanged full partition")
    parser.add_argument("--boundary-probe",action="store_true",help="evaluate retained points immediately below/above the steric cutoff")
    parser.add_argument("--wait", action="store_true", help="wait for each tight-minimum selection to be ready")
    parser.add_argument("--force", action="store_true", help="recompute quadrature instead of reusing matching inputs")
    args = parser.parse_args()
    if args.samples < 64:
        raise ValueError("quadrature needs at least 64 samples")
    scratch = allowed_scratch(args.scratch)
    physical = {}
    for cpu in sorted(os.sched_getaffinity(0)):
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        core = (topology / "physical_package_id").read_text(), (topology / "core_id").read_text()
        physical.setdefault(core, cpu)
    os.sched_setaffinity(0, sorted(physical.values())[:8])
    for name in args.species:
        if args.boundary_probe:
            if not args.indices:raise ValueError("boundary probes require explicit well indices")
            for index in args.indices:probe_boundary(scratch,name,index,args.samples,args.seed,args.proposal)
            continue
        marker = scratch / "composite" / name / "selection.json"
        if args.wait:
            print(f"I043 coupled rotors awaiting tight-minimum selection: {name}", flush=True)
            while not marker.exists():
                time.sleep(30)
        integrate(scratch, name, args.samples, args.seed,args.force,args.proposal,args.indices)


if __name__ == "__main__":
    main()
