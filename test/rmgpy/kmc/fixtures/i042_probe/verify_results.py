"""Independently check source snapshot, numbers, thermodynamic identities and report.

Command: MPLCONFIGDIR=/home/alon/runs/i042-melt-reference/mpl PYTHONPATH=$PWD
/home/alon/anaconda3/envs/rmg_env/bin/python
test/rmgpy/kmc/fixtures/i042_probe/verify_results.py
/home/alon/runs/i042-melt-reference/reproduce/results.json
--snapshot /home/alon/runs/i042-melt-reference/reproduce/database

No generator import. PC-SAFT is checked against external component regressions
and independent finite-difference derivatives; this is not experimental validation.
"""

from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

import numpy as np
from scipy.optimize import brentq
from rmgpy.data.rmg import RMGDatabase
from rmgpy.molecule.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species

import models as m
from render_tables import verify_or_update


def close(actual, expected, atol=3e-6, rtol=2e-9):
    np.testing.assert_allclose(actual, expected, atol=atol, rtol=rtol)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("--snapshot", type=Path, required=True)
    args = parser.parse_args()
    r = json.loads(args.results.read_text())
    repo = Path(__file__).resolve().parents[5]
    sha = r["provenance"]["database_sha"]
    assert sha == "4a12d36fcdc193ede82c8d1ab5c1653495d445bc"
    digest = hashlib.sha256()
    prefixes = r["provenance"]["snapshot"]["prefixes"]
    paths = subprocess.check_output(["git", "-C", "/home/alon/Code/RMG-database",
                                     "ls-tree", "-r", "--name-only", sha, "--", *prefixes]).decode().splitlines()
    assert len(paths) == r["provenance"]["snapshot"]["files"]
    for p in paths:
        assert "catalog" not in Path(p).parts
        expected = subprocess.check_output(["git", "-C", "/home/alon/Code/RMG-database", "show", sha+":"+p])
        assert (args.snapshot / p).read_bytes() == expected
        digest.update(p.encode()+b"\0"+expected)
    assert digest.hexdigest() == r["provenance"]["snapshot"]["sha256"]
    for p, expected in r["provenance"]["source_sha256"].items():
        assert hashlib.sha256((repo / p).read_bytes()).hexdigest() == expected
    db = RMGDatabase()
    db.load_thermo(str(args.snapshot / "input/thermo"), thermo_libraries=["primaryThermoLibrary"], depository=True)
    db.load_solvation(str(args.snapshot / "input/solvation"))
    participants = [Species(molecule=[Molecule().from_adjacency_list(v["adjacency"])]) for v in r["species"]]
    for s, saved in zip(participants, r["species"]):
        assert s.molecule[0].get_element_count() == saved["formula"]
        s.thermo = db.thermo.get_thermo_data(s)
    reaction = Reaction(reactants=participants[:2], products=participants[2:])
    # Independently form concentration-standard G from pressure-standard H/S.
    def gc(t):
        hp = sum(s.thermo.get_enthalpy(t) for s in reaction.products)-sum(s.thermo.get_enthalpy(t) for s in reaction.reactants)
        sp = sum(s.thermo.get_entropy(t) for s in reaction.products)-sum(s.thermo.get_entropy(t) for s in reaction.reactants)
        return hp-t*sp-m.R*t*math.log(m.C0*m.R*t/1e5)
    root = brentq(gc, 600, 800)
    close(root, r["gas"]["continuous_K"])
    for v, e in zip(r["gas"]["rows"], r["equilibrium"]):
        t = v["T_K"]
        close(gc(t), v["G"])
        close(gc(t), v["H"]-t*v["S"])
        close(v["H"], reaction.get_enthalpy_of_reaction(t)+m.R*t)
        close(v["S"], reaction.get_entropy_of_reaction(t)+m.R*(math.log(m.C0*m.R*t/1e5)+1))
        close(math.exp(gc(t)/(m.R*t)), e["gas_equilibrium_c_mol_L"])
    for t, kp, kd in zip(r["rates"]["propagation"]["T"], r["rates"]["propagation"]["k"], r["rates"]["depropagation"]["k"]):
        close(kp/kd, math.exp(-gc(t)/(m.R*t))/m.C0, atol=1e-12)
    for v in r["concentrations"]:
        c = v["fraction_free_styrene"]*m.RHO*1000/m.MW
        close(c/1000, v["c_mol_L"])
        expected = brentq(lambda t: gc(t)-m.R*t*math.log(c/m.C0), 400, 950)
        close(expected, v["gas_Tc_K"])
        close(expected-root, v["gas_shift_K"])
        close(m.R*math.log(c/m.C0), v["concentration_S_shift"])
        close(m.fh_activity(700, c), v["FH_athermal_activity"])
        close(m.unifac_values(700, m.unifac_state(c))["monomer_activity"], v["UNIFAC_activity700"])
    close(m.CTOT/1000, r["definitions"]["repeat_c_mol_L"])
    close(1/(m.NA*m.CTOT), r["definitions"]["repeat_volume_m3"], atol=0)
    print("I042 pinned database, source hashes, gas thermo, rates and volume/concentration roots verified")

    m.self_check()
    masses = r["definitions"]["radical_masses_g_mol"]
    for label, pc in r["PC"].items():
        assert pc["parameters"] == m.PC_SETS[label]
        state = m.pc_state(label, masses)
        c = state[0]
        for v in pc["rows"]:
            t = v["T_K"]
            mu_fd = []
            for i in range(4):
                step = 1e-3
                samples = []
                for offset in (-2, -1, 1, 2):
                    trial = c.copy(); trial[i] += offset*step
                    samples.append(m.helmholtz_density(t, trial, *state[1:]))
                mu_fd.append((samples[0]-8*samples[1]+8*samples[2]-samples[3])/(12*step))
            close(np.array(mu_fd)*m.R*t, v["mu_J_mol"], atol=0.01, rtol=1e-8)
            # Reaction-direction derivative, rather than subtracting chemical potentials.
            direction = np.array([0., -1., -1., 1.])
            step = 1e-3
            dg = m.R*t*(m.helmholtz_density(t, c+step*direction, *state[1:])
                       -m.helmholtz_density(t, c-step*direction, *state[1:]))/(2*step)
            close(dg, v["G"], atol=0.005)
            gt = (m.pc_values(t+0.1, state)["G"]-m.pc_values(t-0.1, state)["G"])/0.2
            close(-gt, v["S_path"], atol=2e-5)
            close(v["G"]-t*gt, v["H_path"], atol=0.02)
            # Independent pressure derivative of A at fixed amounts (scale=1/V).
            eps = 1e-5
            fplus = m.helmholtz_density(t, c/(1+eps), *state[1:])*(1+eps)
            fminus = m.helmholtz_density(t, c/(1-eps), *state[1:])*(1-eps)
            pressure = m.R*t*(sum(c)-(fplus-fminus)/(2*eps))
            close(pressure, v["P_Pa"], atol=1, rtol=1e-8)
            # Hold this state's pressure and composition fixed at neighboring T.
            isobar_g = []
            for u in (t-0.1, t+0.1):
                scale = brentq(lambda scale: m.pc_values(u, (c*scale, *state[1:]))["P_Pa"]-v["P_Pa"], 0.99, 1.01)
                isobar_g.append(m.pc_values(u, (c*scale, *state[1:]))["G"])
            sp = -(isobar_g[1]-isobar_g[0])/0.2
            close(sp, v["S_pressure"], atol=2e-5)
            close(v["G"]+t*sp, v["H_pressure"], atol=0.02)
            close(v["monomer_concentration_activity"], math.exp(v["mu_J_mol"][1]/(m.R*t)))
            close(v["pair_equilibrium_multiplier"], math.exp(-v["G"]/(m.R*t)))
            hgas = reaction.get_enthalpy_of_reaction(t)+m.R*t
            sgas = (hgas-gc(t))/t
            for path in ("path", "pressure"):
                close(hgas+v["H_"+path], v["corrected_H_"+path])
                close(sgas+v["S_"+path], v["corrected_S_"+path])
            assert v["dP_dlogrho_Pa"] > 0 and v["stability_min_eigenvalue_m3_mol"] > 0
        expected = brentq(lambda t: gc(t)+m.pc_values(t, state)["G"], 600, 900)
        close(expected, pc["roots_K"][0])
        close(expected-root, pc["shifts_K"][0])
        for v in pc["isobar_1bar"]:
            state_p = m.pc_state(label, masses, density=v["density_kg_m3"])
            actual = m.pc_transfer(v["T_K"], state_p)
            close(actual["P_Pa"], 1e5, atol=0.1)
            for key in ("G", "H_pressure", "S_pressure"):
                close(actual[key], v[key])
        for v in pc["controls"]:
            control = m.pc_state(label, masses, density=v["density_kg_m3"], dp=v["DP"])
            close(m.pc_values(700, control)["G"], v["G700"])
            for tc in v["roots_K"]:
                close(gc(tc)+m.pc_values(tc, control)["G"], 0, atol=1e-5)
    print("I042 PC-SAFT transfer, isochoric/isobaric derivatives, pressure and roots independently checked")

    for v in r["UNIFAC"]:
        t = v["T_K"]
        state = m.unifac_state()
        actual = m.unifac_values(t, state)
        for key in actual:
            close(actual[key], v[key])
        sp = -(m.unifac_values(t+0.1, state)["G_mixing_1M"]-m.unifac_values(t-0.1, state)["G_mixing_1M"])/0.2
        close(sp, v["S_mixing_1M"], atol=2e-6)
        close(v["G_mixing_1M"]+t*sp, v["H_mixing_1M"], atol=0.002)
        # Gibbs-Duhem consistency: derivative of total excess G produces ln(gamma).
        groups, x, _ = state
        def excess(amounts):
            return np.dot(amounts, m.unifac_log_gamma(t, groups, amounts/sum(amounts)))
        for i in range(4):
            h = 1e-6
            plus, minus = x.copy(), x.copy()
            plus[i] += h; minus[i] -= h
            close((excess(plus)-excess(minus))/(2*h), v["log_gamma"][i], atol=2e-6)
        phi = m.C0/m.CTOT
        close(m.fh_activity(t, m.C0, chi_a=v["effective_chi"]), v["monomer_activity"])
    for v in r["FH_athermal"]:
        coefficient = math.log(1.5)+math.log(m.CTOT/m.C0)-1
        close(m.R*v["T_K"]*coefficient, v["G_mixing_1M"])
        close(-m.R*coefficient, v["S_mixing_1M"])
        close(v["G_mixing_1M"]+v["T_K"]*v["S_mixing_1M"], v["H_mixing_1M"])
        close(m.R*v["T_K"]*(2*m.C0/m.CTOT-1), v["dG_dchi_J_mol"])
    # Independently differentiate a multicomponent FH lattice free energy.
    # These arbitrary nonzero chi coefficients exercise algebra, not chemistry.
    c = np.array([(m.CTOT-m.C0)/m.DP, m.C0, 1e-3, 1e-3])
    sizes = np.array([m.DP, 1., 2., 3.])
    chi = 0.123+47./700
    def fh_total(amounts):
        lattice = np.dot(amounts, sizes)
        fractions = amounts*sizes/lattice
        return np.dot(amounts, np.log(fractions))+chi*amounts[1]*(lattice-amounts[1])/lattice
    mu = []
    for i in range(4):
        trial = c.astype(complex); trial[i] += 1e-20j
        mu.append(fh_total(trial).imag/1e-20-math.log(c[i]/m.C0))
    lattice = np.dot(c, sizes)
    close(mu[3]-mu[2]-mu[1], math.log(1.5)+math.log(lattice/m.C0)-1+chi*(2*c[1]/lattice-1))
    for label, v in r["solvents"].items():
        solvent = db.solvation.get_solvent_data(label)
        for field, value in v["coefficients"].items():
            assert getattr(solvent, field) == value
        if "H" not in v:
            assert any(value is None for value in v["coefficients"].values())
            continue
        dh = ds = 0.
        for species, saved in zip(participants, r["species"]):
            solute = db.solvation.get_solute_data(species)
            h = 1000*(solvent.c_h+sum(getattr(solute, a)*getattr(solvent, b+"_h")
                                    for a, b in (("S", "s"), ("B", "b"), ("E", "e"), ("L", "l"), ("A", "a"))))
            g = -8.314*298*2.303*(solvent.c_g+sum(getattr(solute, a)*getattr(solvent, b+"_g")
                                              for a, b in (("S", "s"), ("B", "b"), ("E", "e"), ("L", "l"), ("A", "a"))))
            dh += saved["sign"]*h; ds += saved["sign"]*(h-g)/298
        close(dh, v["H"]); close(ds, v["S"])
        for row in v["rows"]:
            close(dh-row["T_K"]*ds, row["G"])
        tc = brentq(lambda t: gc(t)+dh-t*ds, 600, 900)
        close(tc, v["roots_K"][0]); close(tc-root, v["shifts_K"][0])
        for row in v["kfactor"]:
            def g(u):
                return sum(saved["sign"]*db.solvation.get_T_dep_solvation_energy_from_LSER_298(
                    db.solvation.get_solute_data(species), solvent, u)[0] for species, saved in zip(participants, r["species"]))
            if "error" in row:
                try:
                    g(row["T_K"])
                except Exception as error:
                    assert type(error).__name__+": "+str(error) == row["error"]
                else:
                    raise AssertionError("expected a validity guard")
            else:
                dg = g(row["T_K"])
                entropy = -(g(row["T_K"]+0.05)-g(row["T_K"]-0.05))/0.1
                close(dg, row["G"])
                close(entropy, row["S"], atol=0.0001)
                close(dg+row["T_K"]*entropy, row["H"], atol=0.1)
        if v["critical_K"] is not None:
            # Independently scan the valid molecular-liquid domain for sign crossings.
            grid = np.linspace(300, v["critical_K"]-0.1, 501)
            residuals = [gc(t)+g(t) for t in grid]
            expected = [brentq(lambda t: gc(t)+g(t), a, b) for a, b, x, y in
                        zip(grid[:-1], grid[1:], residuals[:-1], residuals[1:]) if x*y < 0]
            close(expected, v["kfactor_roots_K"])
    print("I042 FH/UNIFAC mixing quantities, LSER parameters and K-factor validity guards verified")
    verify_or_update(r)


if __name__ == "__main__":
    main()
