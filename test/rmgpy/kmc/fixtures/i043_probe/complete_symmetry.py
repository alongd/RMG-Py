"""Close rotor-space pools under reflection and proper atom permutations.

Commands are in the report. Exact spatial reflection/permutation transforms
the Cartesian Hessian; its frequencies and electronic energy are unchanged.
Chiral stereoisomers are not supplemented with their opposite enantiomer.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign

from run_ensembles import SPECIES, SCRATCH, allowed_scratch, read_xyz, limit_cores
from rotors import make_molecule, rotor_definitions, angles, wrap
from summarize import stereo_check


def materialize_energy(scratch, name, index, level):
    directory = scratch / "xtb_checks" / name / f"{index:04d}"
    record = json.loads((directory / "result.json").read_text())
    kind = "derived_by_reflection_from" if "derived_by_reflection_from" in record else "derived_by_permutation_from"
    parent = record.get(kind)
    if parent is None:
        return False
    source = scratch / "composite" / name / f"{parent:04d}" / level / "result.json"
    if not source.exists():
        return False
    result = json.loads(source.read_text())
    result.pop("derived_by_reflection_from",None)
    result.pop("derived_by_permutation_from",None)
    frame = read_xyz(directory / "xtbopt.xyz")[0]
    result.update(input_sha256=hashlib.sha256((directory / "xtbopt.xyz").read_bytes()).hexdigest(),
                  coordinates_angstrom=[[float(v) for v in line.split()[1:]] for line in frame["atoms"]],
                  parent_result_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    result[kind] = parent
    result.pop("electronic_s", None)
    target = scratch / "composite" / name / f"{index:04d}" / level
    target.mkdir(parents=True, exist_ok=True)
    (target / "result.json").write_text(json.dumps(result, indent=2)+"\n")
    return True


def permutation_closure(scratch, name, selection):
    """Restore label-space wells removed by graph-symmetry RMS deduplication.

    The full torsion domain is divided by the external symmetry factor. Its
    proper-permutation images must therefore be represented before quadrature.
    Internal methyl/phenyl periods already quotient their own permutations.
    """
    indices = list(selection["unique_minima_indices"])
    molecules = [make_molecule(SPECIES[name],scratch/"xtb_checks"/name/f"{i:04d}"/"xtbopt.xyz") for i in indices]
    definitions = rotor_definitions(molecules[0])
    periods = np.array([2*math.pi/r["symmetry"] for r in definitions])
    centers = [wrap(angles(m,definitions),periods) for m in molecules]
    additions = []
    for index,molecule in zip(indices,molecules):
        xyz = np.array(molecule.GetConformer().GetPositions())
        tetrahedra = [sorted(a.GetIdx() for a in atom.GetNeighbors())
                      for atom in molecule.GetAtoms() if atom.GetDegree()==4]
        def volume(order):
            points=xyz[order]
            return float(np.linalg.det(points[:3]-points[3]))
        volumes=[volume(order) for order in tetrahedra]
        matches=molecule.GetSubstructMatches(molecule,uniquify=False,useChirality=True,maxMatches=1000000)
        if len(matches)==1000000:
            raise AssertionError("permutation enumeration truncated")
        proper=[p for p in matches if all(v*volume([p[i] for i in order])>0
                                          for v,order in zip(volumes,tetrahedra))]
        if len(proper)/math.prod(r["symmetry"] for r in definitions)<=1:
            continue
        for mapping in proper:
            image=Chem.Mol(molecule)
            coordinates=xyz[list(mapping)]
            for i,point in enumerate(coordinates):
                image.GetConformer().SetAtomPosition(i,point)
            center=wrap(angles(image,definitions),periods)
            if any(np.max(np.abs(wrap(center-old,periods)))<=0.05 for old in centers):
                continue
            centers.append(center)
            image_index=100000+len(additions)
            source=scratch/"xtb_checks"/name/f"{index:04d}"
            target=scratch/"xtb_checks"/name/f"{image_index:04d}"
            target.mkdir(parents=True,exist_ok=True)
            block=Chem.MolToXYZBlock(image)
            (target/"xtbopt.xyz").write_text(block)
            stereo_check(SPECIES[name],read_xyz(target/"xtbopt.xyz")[0])
            values=[float(v) for line in (source/"hessian").read_text().splitlines()
                    if not line.startswith("$") for v in line.split()]
            hessian=np.array(values).reshape((3*molecule.GetNumAtoms(),)*2)
            permutation=np.array([3*i+a for i in mapping for a in range(3)])
            transformed=hessian[np.ix_(permutation,permutation)]
            if not np.allclose(np.linalg.eigvalsh(hessian),np.linalg.eigvalsh(transformed),atol=1e-10):
                raise AssertionError("atom permutation changed Hessian eigenvalues")
            (target/"hessian").write_text("$hessian\n"+"\n".join(" ".join(f"{v:.15e}" for v in row)
                                                               for row in transformed)+"\n$end\n")
            record=json.loads((source/"result.json").read_text())
            record.pop("derived_by_reflection_from",None)
            record.pop("reflection_atom_permutation",None)
            record.update(derived_by_permutation_from=index,permutation_atom_order=list(mapping),
                          xyz_sha256=hashlib.sha256(block.encode()).hexdigest(),
                          parent_hessian_sha256=hashlib.sha256((source/"hessian").read_bytes()).hexdigest(),
                          command=["complete_symmetry.py","exact proper atom permutation"],elapsed_s=0.0)
            (target/"result.json").write_text(json.dumps(record,indent=2)+"\n")
            additions.append({"index":image_index,"parent":index})
    selection["chemical_minima_count"]=len(indices)
    selection["unique_minima_indices"]=indices+[a["index"] for a in additions]
    selection["permutation_partners"]=additions
    selection["angular_orbit_tolerance_rad"]=0.05
    return selection


def augment(scratch, name, selection):
    manifest = json.loads((scratch / "species.json").read_text())
    originals = [i for i in selection["unique_minima_indices"] if i < 10000]
    molecules = [make_molecule(SPECIES[name], scratch / "xtb_checks" / name / f"{i:04d}" / "xtbopt.xyz")
                 for i in originals]
    retained = [Chem.RemoveHs(molecule) for molecule in molecules]
    additions = []
    for index, molecule in zip(originals, molecules):
        if not manifest[name]["achiral"]:
            continue
        reflected = Chem.Mol(molecule)
        original_coordinates = np.array(molecule.GetConformer().GetPositions())
        coordinates = original_coordinates.copy()
        coordinates[:, 0] *= -1
        for i,point in enumerate(coordinates):
            reflected.GetConformer().SetAtomPosition(i,point)
        Chem.AssignStereochemistryFrom3D(reflected, replaceExistingTags=True)
        expected = Chem.AddHs(Chem.MolFromSmiles(SPECIES[name]))
        tetrahedra=[sorted(a.GetIdx() for a in atom.GetNeighbors())
                    for atom in molecule.GetAtoms() if atom.GetDegree()==4]
        def volume(points,order):
            xyz=points[order]
            return np.linalg.det(xyz[:3]-xyz[3])
        volumes=[volume(original_coordinates,order) for order in tetrahedra]
        matches=reflected.GetSubstructMatches(expected,useChirality=True,uniquify=False,maxMatches=1000000)
        if len(matches)==1000000:raise AssertionError("reflection mapping enumeration truncated")
        mapping=next((p for p in matches if all(v*volume(coordinates,[p[i] for i in order])>0
                      for v,order in zip(volumes,tetrahedra))),None)
        if not mapping:
            raise AssertionError("reflection changed an allegedly achiral stereoisomer")
        reordered = Chem.Mol(molecule)
        coordinates = coordinates[list(mapping)]
        for i,point in enumerate(coordinates):
            reordered.GetConformer().SetAtomPosition(i,point)
        heavy = Chem.RemoveHs(reordered)
        if any(rdMolAlign.GetBestRMS(Chem.Mol(heavy), earlier, maxMatches=10000) < 0.10 for earlier in retained):
            continue
        retained.append(heavy)
        mirror_index = 10000+index
        source = scratch / "xtb_checks" / name / f"{index:04d}"
        target = scratch / "xtb_checks" / name / f"{mirror_index:04d}"
        target.mkdir(parents=True,exist_ok=True)
        xyz = Chem.MolToXYZBlock(reordered)
        (target / "xtbopt.xyz").write_text(xyz)
        stereo_check(SPECIES[name],read_xyz(target / "xtbopt.xyz")[0])
        values = [float(v) for line in (source / "hessian").read_text().splitlines()
                  if not line.startswith("$") for v in line.split()]
        hessian = np.array(values).reshape((3*molecule.GetNumAtoms(),)*2)
        permutation = np.array([3*i+a for i in mapping for a in range(3)])
        signs = np.tile([-1,1,1],molecule.GetNumAtoms())
        transformed = hessian[np.ix_(permutation,permutation)]*signs[:,None]*signs[None,:]
        if not np.allclose(np.linalg.eigvalsh(hessian),np.linalg.eigvalsh(transformed),atol=1e-10):
            raise AssertionError("reflection changed Hessian eigenvalues")
        (target / "hessian").write_text("$hessian\n"+"\n".join(" ".join(f"{v:.15e}" for v in row)
                                                                  for row in transformed)+"\n$end\n")
        record = json.loads((source / "result.json").read_text())
        record.update(derived_by_reflection_from=index, reflection_atom_permutation=list(mapping),
                      xyz_sha256=hashlib.sha256(xyz.encode()).hexdigest(),
                      parent_hessian_sha256=hashlib.sha256((source / "hessian").read_bytes()).hexdigest(),
                      command=["complete_symmetry.py", "exact reflection and atom permutation"],
                      elapsed_s=0.0)
        (target / "result.json").write_text(json.dumps(record,indent=2)+"\n")
        additions.append({"index": mirror_index, "parent": index})
    selection["unique_minima_indices"] = originals+[row["index"] for row in additions]
    selection["reflection_partners"] = additions
    selection["reflection_mapping"]="preserve signed tetrahedral volumes, including CH2 and CH3 labels"
    selection = permutation_closure(scratch,name,selection)
    for index in selection["unique_minima_indices"]:
        for level in ("pbe","blyp"):
            materialize_energy(scratch,name,index,level)
    return selection


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch",type=Path,default=SCRATCH)
    parser.add_argument("--species",choices=tuple(SPECIES),nargs="+")
    args = parser.parse_args()
    limit_cores()
    scratch = allowed_scratch(args.scratch)
    for name in args.species or SPECIES:
        path = scratch / "composite" / name / "selection.json"
        selection = augment(scratch,name,json.loads(path.read_text()))
        path.write_text(json.dumps(selection,indent=2)+"\n")
        print(f"I043 {name}: {len(selection['reflection_partners'])} reflection and "
              f"{len(selection['permutation_partners'])} proper-permutation rotor images added")


if __name__ == "__main__":
    main()
