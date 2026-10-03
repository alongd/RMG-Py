"""Optimize one XYZ conformer and calculate its Hessian with pyscf_env.

Commands and fixed levels are documented in I043_diphenyl_nonadditivity.md.
Results preserve energies, gradient, Cartesian Hessian, normal modes and timings.
An optimization checkpoint is written before the more expensive Hessian.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import time
import tempfile

import numpy as np
from pyscf import dft, gto, lib, __version__
from pyscf.geomopt.geometric_solver import kernel
from pyscf.hessian import thermo

SCRATCH = Path("/home/alon/runs/i043-diphenyl-nonadditivity")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xyz", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--level", choices=("pbe", "blyp"), required=True)
    parser.add_argument("--opt-only", action="store_true")
    parser.add_argument("--single-point", action="store_true",
                        help="composite electronic energy at an already validated xTB minimum")
    args = parser.parse_args()
    physical = {}
    for cpu in sorted(os.sched_getaffinity(0)):
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        core = (topology / "physical_package_id").read_text(), (topology / "core_id").read_text()
        physical.setdefault(core, cpu)
    os.sched_setaffinity(0, sorted(physical.values())[:8])
    target = args.output.resolve()
    if SCRATCH not in target.parents:
        raise ValueError("output must be under dispatched scratch")
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.parent / "tmp"
    temporary.mkdir(exist_ok=True)
    tempfile.tempdir = str(temporary)
    lib.param.TMPDIR = str(temporary)
    lib.num_threads(8)
    lines = args.xyz.read_text().splitlines()
    n = int(lines[0])
    molecule = gto.M(atom="\n".join(lines[2:2+n]), basis="def2-svp", verbose=4)
    molecule.max_memory = 12000

    def mean_field(mol):
        calculation = dft.RKS(mol, xc=args.level).density_fit()
        calculation.disp = "d3bj"
        calculation.grids.level = 3
        calculation.conv_tol = 1e-9
        calculation.max_cycle = 100
        return calculation

    start = time.monotonic()
    if args.single_point:
        calculation = mean_field(molecule)
        energy = float(calculation.kernel())
        if not calculation.converged:
            raise RuntimeError("composite SCF did not converge")
        result = {"pyscf_version": __version__,
                  "level": args.level + "-D3(BJ)/def2-SVP//GFN2-xTB",
                  "input_sha256": hashlib.sha256(args.xyz.read_bytes()).hexdigest(),
                  "energy_hartree": energy,
                  "dispersion_hartree": float(calculation.get_dispersion()),
                  "electronic_s": time.monotonic()-start,
                  "atoms": [molecule.atom_symbol(i) for i in range(molecule.natm)],
                  "coordinates_angstrom": molecule.atom_coords(unit="Angstrom").tolist(),
                  "masses_amu": molecule.atom_mass_list().tolist(),
                  "frequencies_complete": False, "composite_single_point": True}
        target.write_text(json.dumps(result, indent=2) + "\n")
        print(f"I043 composite energy {args.xyz} {args.level}: {result['electronic_s']:.1f} s", flush=True)
        return
    converged, optimized = kernel(mean_field(molecule), maxsteps=150,
                                 convergence_energy=1e-6,
                                 convergence_grms=1e-4,
                                 convergence_gmax=3e-4,
                                 convergence_drms=6e-4,
                                 convergence_dmax=9e-4)
    if not converged:
        raise RuntimeError("DFT geometry did not converge")
    optimized.tofile(str(target.with_suffix(".xyz")), format="xyz")
    calculation = mean_field(optimized)
    energy = float(calculation.kernel())
    if not calculation.converged:
        raise RuntimeError("final SCF did not converge")
    gradient = calculation.nuc_grad_method().kernel()
    result = {"pyscf_version": __version__, "level": args.level + "-D3(BJ)/def2-SVP",
              "input_sha256": hashlib.sha256(args.xyz.read_bytes()).hexdigest(),
              "energy_hartree": energy, "dispersion_hartree": float(calculation.get_dispersion()),
              "gradient_Eh_bohr": gradient.tolist(),
              "optimization_s": time.monotonic()-start,
              "atoms": [optimized.atom_symbol(i) for i in range(optimized.natm)],
              "coordinates_angstrom": optimized.atom_coords(unit="Angstrom").tolist(),
              "masses_amu": optimized.atom_mass_list().tolist(),
              "frequencies_complete": False}
    target.write_text(json.dumps(result, indent=2) + "\n")
    if args.opt_only:
        print(f"I043 optimized {args.xyz} {args.level}: {result['optimization_s']:.1f} s", flush=True)
        return
    start = time.monotonic()
    hessian = calculation.Hessian().kernel()
    modes = thermo.harmonic_analysis(optimized, hessian)
    frequencies = modes["freq_wavenumber"]
    result.update(hessian_Eh_bohr2=hessian.tolist(),
                  frequencies_cm1_real=np.real(frequencies).tolist(),
                  frequencies_cm1_imag=np.imag(frequencies).tolist(),
                  normal_modes=modes["norm_mode"].tolist(),
                  hessian_s=time.monotonic()-start,
                  frequencies_complete=True)
    target.write_text(json.dumps(result, indent=2) + "\n")
    if np.max(np.abs(np.imag(frequencies))) > 20:
        raise RuntimeError("refined geometry has a significant imaginary frequency")
    print(f"I043 refined minimum {args.xyz} {args.level}: optimization {result['optimization_s']:.1f} s, Hessian {result['hessian_s']:.1f} s", flush=True)


if __name__ == "__main__":
    main()
