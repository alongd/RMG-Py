"""Reestimate I040 from the pinned database and verify its entire report.

See I040_gas_thermo_bench.md for commands, literature provenance, and limits.
No network, database mutation, rate-tree generation, or QM calculation occurs.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

from scipy.integrate import quad
from scipy.optimize import brentq

from rmgpy import constants
from rmgpy.data.rmg import RMGDatabase
from rmgpy.kmc.compiler import ps_proxy_set

FIXTURES = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("i039", FIXTURES / "i039_probe/run_probe.py")
i039 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(i039)
prior = i039.prior
R, T0, P0, C0 = constants.R, 298.15, 100000.0, 1000.0
TEMPS = [T0, 300., 400., 500., 600., 700., 800., 1000., 1500.]

SMILES = {
    "styrene": "C=Cc1ccccc1",
    "ethylbenzene": "CCc1ccccc1",
    "cumene": "CC(C)c1ccccc1",
    "1,3-diphenylpropane": "c1ccccc1CCCc1ccccc1",
    "1,3-diphenylbutane": "CC(CCc1ccccc1)c1ccccc1",
    "2,4-diphenylpentane": "CC(c1ccccc1)CC(C)c1ccccc1",
    "1-phenylethyl radical": "C[CH]c1ccccc1",
    "cumyl radical": "C[C](C)c1ccccc1",
    "ethane": "CC",
    "n-propylbenzene": "CCCc1ccccc1",
    "P1 primary radical": "[CH2]Cc1ccccc1",
    "P2 primary radical": "[CH2]C(CCc1ccccc1)c1ccccc1",
    "P3 primary radical": "[CH2]C(CC(CCc1ccccc1)c1ccccc1)c1ccccc1",
    "B2 benzylic radical": "CC(c1ccccc1)C[CH]c1ccccc1",
    "B3 benzylic radical": "CC(c1ccccc1)CC(c1ccccc1)C[CH]c1ccccc1",
}

# Literal, independently published inputs, never adjusted using a Tc target.
# H/uH in kJ/mol; S/uS and Cp in J/mol/K. None means NOT supplied.
REFERENCES = {
    "styrene": dict(H=148.59, uH=.51, S=345.1, uS=2.1, T=T0,
                    Hsource="A-sty", Ssource="N-sty",
                    Cp=[120.19, 120.94, 159.79, 192.59, 219., 240.4, 258., 285.2, 325.2],
                    CpT=TEMPS, Cpsource="N-sty", uCp=None),
    "ethylbenzene": dict(H=29.98, uH=.52, S=360.6, uS=.5, T=T0,
                         Hsource="A-ethyl", Ssource="N-ethyl",
                         Cp=[127.40, 128.19, 169.95, 206.58, 236.75, 261.51, 282.08, 314.04, 361.27],
                         CpT=TEMPS, Cpsource="N-ethyl", uCp=None),
    "cumene": dict(H=3.9, uH=1.1, S=386.53, uS=None, T=T0,
                   Hsource="N-cum", Ssource="N-cum"),
    "n-propylbenzene": dict(H=7.82, uH=.84, S=397.86, uS=None, T=T0,
                           Hsource="N-nprop", Ssource="N-nprop"),
    "ethane": dict(H=-84.05, uH=.12, S=229.16, uS=.10, T=T0,
                   Hsource="A-index", Ssource="N-ethane"),
    "1,3-diphenylpropane": dict(H=114.7, uH=4.1, S=498.9, uS=None, T=T0,
                               Hsource="V24", Ssource="V24"),
    "1-phenylethyl radical": dict(H=174.9, uH=None, S=359.0, uS=None, T=298.,
                                Hsource="I17", Ssource="I17", raw_G4_H=178.1,
                                Cp=[124.6, 164.6, 199., 227.1, 269.1, 298.8, 343.1],
                                CpT=[300., 400., 500., 600., 800., 1000., 1500.],
                                Cpsource="I17", uCp=None),
}
# Gas C-H dissociation inputs from the primary JACS Table 3 image. The
# molecular radical references derived below assume the usual 298.15 K gas
# convention; the SI pressure/rotor details were not retrievable. They remain
# explicitly conditional rather than being silently promoted to measured data.
H_ATOM = dict(H=217.998, uH=.006, S=114.717, uS=.002)
BOND_REFERENCES = {
    "ethylbenzene": dict(radical="1-phenylethyl radical", CBS_BDE_kcal=87.6,
                         CBS_BDFE_kcal=79.1, Luo_BDE_kcal=85.4),
    "cumene": dict(radical="cumyl radical", CBS_BDE_kcal=86.9,
                   CBS_BDFE_kcal=77.4, Luo_BDE_kcal=83.2),
}
_cum = BOND_REFERENCES["cumene"]
REFERENCES["cumyl radical"] = dict(
    H=REFERENCES["cumene"]["H"] + _cum["CBS_BDE_kcal"] * 4.184 - H_ATOM["H"],
    uH=None,
    S=REFERENCES["cumene"]["S"] + (_cum["CBS_BDE_kcal"] - _cum["CBS_BDFE_kcal"]) * 4184 / T0 - H_ATOM["S"],
    uS=None, T=T0, Hsource="J21", Ssource="J21", conditional_derived=True)
SOURCES = {
    "A-sty": ("ATcT v1.140 (2024), phenylethene, species 312; H298, not H0",
              "https://atct.anl.gov/Thermochemical%20Data/version%201.140/species/?species_number=312"),
    "A-ethyl": ("ATcT v1.140 (2024), ethylbenzene, species 395; H298",
                "https://atct.anl.gov/Thermochemical%20Data/version%201.140/species/?species_number=395"),
    "A-index": ("ATcT v1.140 (2024), ethane, H298",
                "https://atct.anl.gov/Thermochemical%20Data/version%201.140/"),
    "N-sty": ("NIST WebBook, styrene: Pitzer et al. 1946 S; TRC 1997 recommended gas Cp",
              "https://webbook.nist.gov/cgi/cbook.cgi?ID=C100425&Mask=1"),
    "N-ethyl": ("NIST WebBook, ethylbenzene: Miller 1978 S; TRC 1997 recommended gas Cp",
                "https://webbook.nist.gov/cgi/cbook.cgi?ID=C100414&Mask=1"),
    "N-cum": ("NIST WebBook, cumene: Prosen et al. 1945 H; Kishimoto et al. 1973 S",
              "https://webbook.nist.gov/cgi/cbook.cgi?ID=C98828&Mask=1"),
    "N-nprop": ("NIST WebBook, n-propylbenzene: Prosen et al. 1945/1946 H; Messerly et al. 1965 S",
                "https://webbook.nist.gov/cgi/cbook.cgi?ID=C103651&Mask=1"),
    "N-ethane": ("NIST CCCBDB, ethane, Gurvich compilation S298",
                 "https://cccbdb.nist.gov/exp2x.asp?casno=74840&charge=0"),
    "V24": ("Verevkin et al., Oxygen 2024, 4, 266–285, DOI 10.3390/oxygen4030015; Tables 5,7: corrected G3MP2",
            "https://www.mdpi.com/2673-9801/4/3/15"),
    "I17": ("Ince et al., AIChE J. 2017, DOI 10.1002/aic.15588; SI Table S1, T6/53: G4/BAC H, molecular S, Cp",
            "https://biblio.ugent.be/publication/8525511/file/8525521.pdf"),
    "J21": ("Salamone et al., JACS 2021, 143, 11759–11776, DOI 10.1021/jacs.1c05566; Table 3, rows 11,12, gas (RO)CBS-QB3 BDE/BDFE",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC8343544/#tbl3"),
    "N-H": ("NIST WebBook, atomic hydrogen, CODATA 1984 H298 and one-bar S298",
            "https://webbook.nist.gov/cgi/cbook.cgi?ID=C12385136&Mask=1"),
}

REACTIONS = {
    "closed ethyl": {"ethylbenzene": -1, "styrene": -1, "1,3-diphenylbutane": 1},
    "closed cumene": {"cumene": -1, "styrene": -1, "2,4-diphenylpentane": 1},
    "smallest benzylic": {"1-phenylethyl radical": -1, "styrene": -1, "B2 benzylic radical": 1},
    "smallest primary": {"P1 primary radical": -1, "styrene": -1, "P2 primary radical": 1},
    "compiled P2 to P3": {"P2 primary radical": -1, "styrene": -1, "P3 primary radical": 1},
    "benzylic B2 to B3": {"B2 benzylic radical": -1, "styrene": -1, "B3 benzylic radical": 1},
    "balanced group surrogate": {"ethylbenzene": -1, "ethane": -1, "styrene": -1,
                                 "cumene": 1, "n-propylbenzene": 1},
    "two-ring nonadditivity": {"1,3-diphenylpropane": -1, "ethane": -1,
                              "ethylbenzene": 1, "n-propylbenzene": 1},
}


def close(a, b, tol=2e-5):
    if not math.isclose(a, b, rel_tol=2e-10, abs_tol=tol):
        raise AssertionError(f"{a} != {b}")


def thermo(species, t):
    return [float(species.thermo.get_enthalpy(t)),
            float(species.thermo.get_entropy(t)),
            float(species.thermo.get_heat_capacity(t))]


def linear_combination(coeff, mapping):
    return [sum(weight * mapping[name][j] for name, weight in coeff.items()) for j in range(3)]


def net_weights(coeff, decompositions):
    out = Counter()
    for name, weight in coeff.items():
        for label, count in decompositions[name]["source_weights"].items():
            out[label] += weight * count
    return dict(sorted((k, v) for k, v in out.items() if v))


def bisection(function, bounds):
    lo, hi = bounds
    flo, fhi = function(lo), function(hi)
    if flo * fhi > 0:
        return None
    for _ in range(80):
        mid = (lo + hi) / 2
        fm = function(mid)
        if flo * fm <= 0:
            hi = mid
        else:
            lo, flo = mid, fm
    return (lo + hi) / 2


def root(reaction, dh=0., ds=0., cp_correction=None, bounds=(250., 1800.)):
    def changes(t):
        return [dh, ds] if cp_correction is None else cp_correction(t)
    def f(t):
        h, s = i039.pressure_hs(reaction, t)
        hh, ss = changes(t)
        return h + hh - t * (s + ss) - R * t * math.log(C0 * R * t / P0)
    result = bisection(f, bounds)
    if result is None:
        return None
    # A second algorithm and RMG's Kc verify the root and convention.
    close(result, brentq(f, *bounds, xtol=1e-9), 1e-7)
    hh, ss = changes(result)
    close(math.log(reaction.get_equilibrium_constant(result, type="Kc") * C0)
          - (hh - result * ss) / (R * result), 0., 2e-8)
    return result


def interval(reaction, dh, ds, uh, us):
    if uh is None or us is None:
        return None
    vals = [root(reaction, dh + a * uh, ds + b * us) for a in (-1, 1) for b in (-1, 1)]
    return [min(vals), max(vals)] if all(v is not None for v in vals) else None


def interpolated_cp(reference, t):
    ts, cs = reference["CpT"], reference["Cp"]
    if not ts[0] <= t <= ts[-1]:
        raise ValueError("reference Cp extrapolation prohibited")
    for j in range(len(ts) - 1):
        if t <= ts[j + 1]:
            return cs[j] + (cs[j + 1] - cs[j]) * (t - ts[j]) / (ts[j + 1] - ts[j])
    return cs[-1]


def estimate(database):
    species, details = {}, {}
    for name, smiles in SMILES.items():
        item = prior.Species(molecule=[prior.Molecule(smiles=smiles)])
        item.molecule[0].update()
        item.generate_resonance_structures()
        item.thermo = database.thermo.get_thermo_data(item)
        species[name] = item
        details[name] = i039.decompose(database.thermo, item)
        details[name]["HSCp"] = [thermo(item, t) for t in TEMPS]
        details[name]["source_metadata"] = {}
        for label in details[name]["source_weights"]:
            kind, key = label.split(":", 1)
            if kind not in database.thermo.groups:
                continue
            group = database.thermo.groups[kind]
            entry = group.entries[key]
            aliases = [key]
            while isinstance(entry.data, str):
                entry = group.entries[entry.data]
                aliases.append(entry.label)
            details[name]["source_metadata"][label] = dict(
                resolved_aliases=aliases, shortDesc=entry.short_desc, longDesc=entry.long_desc,
                H298=str(entry.data.H298), S298=str(entry.data.S298),
                Cpdata=str(entry.data.Cpdata))
    return species, details


def calculate(database):
    print("[I040] estimating species and compiling I034 baseline", file=sys.stderr, flush=True)
    species, details = estimate(database)
    baseline, rates = prior.gas_baseline(database, ps_proxy_set(3))
    compiled = baseline["propagation"]
    for item in compiled.reactants + compiled.products:
        item.thermo = database.thermo.get_thermo_data(item)
    pcoeff = REACTIONS["compiled P2 to P3"]
    expected = prior.Reaction(reactants=[species[k] for k, v in pcoeff.items() if v == -1],
                              products=[species[k] for k, v in pcoeff.items() if v == 1])
    assert sorted(prior.smiles(compiled.reactants)) == sorted(prior.smiles(expected.reactants))
    assert sorted(prior.smiles(compiled.products)) == sorted(prior.smiles(expected.products))
    reactions = {}
    model_at_anchor = {k: thermo(v, T0) for k, v in species.items()}
    for label, coeff in REACTIONS.items():
        # Atom balance from molecular graphs, including the synthetic surrogate.
        elements = Counter()
        for name, weight in coeff.items():
            for atom in species[name].molecule[0].atoms:
                elements[atom.element.symbol] += weight
        assert all(v == 0 for v in elements.values()), (label, elements)
        hsc = linear_combination(coeff, model_at_anchor)
        reactions[label] = dict(coefficients=coeff, HSCp=hsc, groups=net_weights(coeff, details))
        ref = {k: [v["H"] * 1000, v["S"], 0.] for k, v in REFERENCES.items() if v["T"] == T0}
        if all(k in ref for k in coeff):
            reactions[label]["reference_HS"] = linear_combination(coeff, ref)[:2]
            reactions[label]["reference_H_bound"] = sum(abs(v) * REFERENCES[k]["uH"] * 1000
                                                       for k, v in coeff.items())
        else:
            reactions[label]["missing_reference"] = [k for k in coeff if k not in ref]
    for t in TEMPS:
        actual = i039.pressure_hs(compiled, t)
        wanted = linear_combination(pcoeff, {k: thermo(v, t) for k, v in species.items()})
        close(actual[0], wanted[0]); close(actual[1], wanted[1])
    surrogate = reactions["balanced group surrogate"]
    assert surrogate["groups"] == reactions["compiled P2 to P3"]["groups"], (
        surrogate["groups"], reactions["compiled P2 to P3"]["groups"])
    assert not reactions["two-ring nonadditivity"]["groups"]
    entropy_offset = reactions["compiled P2 to P3"]["HSCp"][1] - surrogate["HSCp"][1]
    close(reactions["compiled P2 to P3"]["HSCp"][0], surrogate["HSCp"][0])
    for t in TEMPS:
        vec = {k: thermo(v, t) for k, v in species.items()}
        a, b = linear_combination(pcoeff, vec), linear_combination(surrogate["coefficients"], vec)
        close(a[0], b[0]); close(a[1] - b[1], entropy_offset); close(a[2], b[2])
    tc = root(compiled)
    comparisons = {}
    for name, ref in REFERENCES.items():
        h, s, cp = thermo(species[name], ref["T"])
        comparisons[name] = dict(T=ref["T"], model_H=h / 1000., model_S=s,
                                 ref_minus_model_H=ref["H"] - h / 1000.,
                                 ref_minus_model_S=ref["S"] - s)
    bond_cycles = {}
    for parent, ref in BOND_REFERENCES.items():
        radical = ref["radical"]
        atom = prior.Species(molecule=[prior.Molecule(smiles="[H]")])
        atom.thermo = database.thermo.get_thermo_data(atom)
        ph, ps, _ = thermo(species[parent], T0)
        rh, rs, _ = thermo(species[radical], T0)
        ah, ass, _ = thermo(atom, T0)
        # direct dissociation comparison does not depend on parent Hf anchors
        bond_cycles[parent] = dict(
            model_BDE_kJ=(rh + ah - ph) / 1000,
            model_BDFE_kJ=(rh + ah - ph - T0 * (rs + ass - ps)) / 1000,
            CBS_BDE_kJ=ref["CBS_BDE_kcal"] * 4.184,
            CBS_BDFE_kJ=ref["CBS_BDFE_kcal"] * 4.184,
            derived_CBS_radical_H=REFERENCES[parent]["H"] + ref["CBS_BDE_kcal"] * 4.184 - H_ATOM["H"],
            derived_CBS_radical_S=REFERENCES[parent]["S"] + (ref["CBS_BDE_kcal"] - ref["CBS_BDFE_kcal"]) * 4184 / T0 - H_ATOM["S"],
            derived_Luo_radical_H=REFERENCES[parent]["H"] + ref["Luo_BDE_kcal"] * 4.184 - H_ATOM["H"],
        )
    effects = []
    for name, coefficient in surrogate["coefficients"].items():
        comp, ref = comparisons[name], REFERENCES[name]
        dh, ds = coefficient * comp["ref_minus_model_H"] * 1000, coefficient * comp["ref_minus_model_S"]
        uh = abs(coefficient) * ref["uH"] * 1000 if ref["uH"] is not None else None
        us = abs(coefficient) * ref["uS"] if ref["uS"] is not None else None
        for axis in ("H", "S"):
            hh, ss = (dh, 0.) if axis == "H" else (0., ds)
            hu, su = (uh, 0.) if axis == "H" else (0., us)
            effects.append(dict(name=name, axis=axis, dH=hh, dS=ss,
                                Tc=root(compiled, hh, ss), interval=interval(compiled, hh, ss, hu, su),
                                assumption="direct monomer substitution; other species fixed" if name == "styrene"
                                else "one surrogate input substitution; exact GAV identity transferred"))
    dh = surrogate["reference_HS"][0] - surrogate["HSCp"][0]
    ds = surrogate["reference_HS"][1] - surrogate["HSCp"][1]
    uh = surrogate["reference_H_bound"]
    # The covariance is available only for this pair; this is an illustration,
    # NOT a combined confidence interval across heterogeneous sources.
    rho = .433
    illustrative_rss = math.sqrt(sum(REFERENCES[k]["uH"]**2 for k in surrogate["coefficients"])
                                 + 2 * rho * REFERENCES["styrene"]["uH"] * REFERENCES["ethylbenzene"]["uH"])
    combined = dict(dH=dh, dS=ds, H_bound=uh, illustrative_H_rss_kJ=illustrative_rss,
                    Tc_H=root(compiled, dh, 0.), Tc_S=root(compiled, 0., ds),
                    Tc_HS=root(compiled, dh, ds),
                    Tc_H_bound=interval(compiled, dh, 0., uh, 0.),
                    entropy_offset=entropy_offset, S_bound=None)
    # Separate anchor-H, anchor-S, and Cp contributions. Styrene is direct;
    # ethylbenzene transfers only through the conditional group surrogate.
    def cp_change(name, t):
        reference = REFERENCES[name]
        if not T0 <= t <= 1500.:
            raise ValueError("Cp transfer outside the published temperature range")
        def dc(u):
            return surrogate["coefficients"][name] * (
                interpolated_cp(reference, u) - species[name].thermo.get_heat_capacity(u))
        points = [v for v in reference["CpT"] if T0 < v < t]
        return [quad(dc, T0, t, points=points, epsabs=1e-6)[0],
                quad(lambda u: dc(u) / u, T0, t, points=points, epsabs=1e-8)[0]]
    cp_effects = {}
    full_input_tc = {}
    for name in ("styrene", "ethylbenzene"):
        for axis in ("H", "S", "HS"):
            def correction(t, axis=axis, name=name):
                h, s = cp_change(name, t)
                return [h if "H" in axis else 0., s if "S" in axis else 0.]
            candidate = root(compiled, cp_correction=correction, bounds=(T0, 1500.))
            assert T0 < candidate < 1500.
            cp_effects[name + ":" + axis] = dict(Tc=candidate, change_at_baseline_Tc=correction(tc), uncertainty=None)
        def full_input(t, name=name):
            h, s = cp_change(name, t)
            coefficient = surrogate["coefficients"][name]
            return [h + coefficient * comparisons[name]["ref_minus_model_H"] * 1000,
                    s + coefficient * comparisons[name]["ref_minus_model_S"]]
        full_input_tc[name] = root(compiled, cp_correction=full_input, bounds=(T0, 1500.))
    def known_cp_surrogate(t):
        corrections = [cp_change(name, t) for name in ("styrene", "ethylbenzene")]
        return [dh + sum(c[0] for c in corrections), ds + sum(c[1] for c in corrections)]
    known_cp_tc = root(compiled, cp_correction=known_cp_surrogate, bounds=(T0, 1500.))
    step = 1e-3
    def f(t):
        h, s = i039.pressure_hs(compiled, t)
        return h - t * s - R * t * math.log(C0 * R * t / P0)
    slope = (f(tc + step) - f(tc - step)) / (2 * step)
    _, ss = i039.pressure_hs(compiled, tc)
    close(slope, -ss - R * (math.log(C0 * R * tc / P0) + 1), 1e-4)
    # Planned core-hours, NOT measured timings. Quotas are explicit and revisable.
    costs = [
        ("small-species G4, optimization/frequency included", 7 * 3, 600.),
        ("large-species DFT optimization/frequency", 7 * 12, 80.),
        ("large-species DLPNO-CCSD(T1)/CBS single-point pair", 7 * 12, 500.),
        ("small-species 24-point relaxed rotor scans", 7 * 3 * 24, .5),
        ("large-species 24-point relaxed rotor scans", 7 * 6 * 24, 2.),
        ("coupled-rotor 12-by-12 grids", 4 * 12 * 12, 2.),
        ("small-species conformer search", 7, 20.),
        ("large-species conformer search", 7, 100.),
        ("canonical G4 cross-checks on large representative conformers", 3 * 2, 1500.),
    ]
    return dict(species=details, comparisons=comparisons, reactions=reactions, bond_cycles=bond_cycles,
                baseline=dict(Tc=tc, grid_Tc=rates["ceiling_tabulated_K"],
                              HSCp=reactions["compiled P2 to P3"]["HSCp"],
                              slope=slope, sensitivity_H_kJ=-1000 / slope, sensitivity_S=tc / slope,
                              reactants=prior.smiles(compiled.reactants), products=prior.smiles(compiled.products)),
                combined=combined, effects=effects, cp_effects=cp_effects, full_input_Tc=full_input_tc,
                known_cp_surrogate_Tc=known_cp_tc,
                pressure_S_offset_per_species=R * math.log(101325 / P0),
                costs=[dict(task=k, jobs=n, core_hours_per_job=c, core_hours=n*c) for k, n, c in costs])


def table(headers, rows):
    def cell(value):
        return str(value).replace("|", "\\|").replace("\n", " ")
    return "\n".join(["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |",
                      *["| " + " | ".join(map(cell, row)) + " |" for row in rows]])


def measured(value, uncertainty, digits=2):
    return f"{value:.{digits}f}" + (f" ± {uncertainty:.{digits}f}" if uncertainty is not None else " (u missing)")


def render(d):
    b, comb = d["baseline"], d["combined"]
    out = ["# I040 — independent gas thermochemistry benchmark\n",
           "Independent small-molecule data constrain RMG's additive propagation increment, but do not establish the real oligomer gas thermochemistry. The balanced surrogate's enthalpy discrepancy lies within the sum of its heterogeneous source uncertainty magnitudes. Entropy, the complete heat-capacity increment, and two-ring nonadditivity remain less securely constrained. The surrogate is an exact identity inside this pinned additive model; transferring its experimental error to an oligomer remains conditional. No correction is selected, and no result is ranked against a ceiling-temperature comparator.\n",
           "## Reproduction and verifier\n",
           "Run from `/home/alon/Code/RMG-Py-kmc-i040-gas-thermo-bench`, branch `i040-gas-thermo-bench`, base `0bbfe8ea33781b34ec6a29d12512aaa8f73f1e45`. Use the supplied environment. Build products are ignored. Both streams must be persisted.\n",
           "```bash\nmkdir -p /home/alon/runs/i040-gas-thermo-bench/build\nPYTHONPATH=$PWD /home/alon/anaconda3/envs/rmg_env/bin/python setup.py build_ext --inplace > >(tee -a /home/alon/runs/i040-gas-thermo-bench/build/stdout.log) 2> >(tee -a /home/alon/runs/i040-gas-thermo-bench/build/stderr.log >&2)\nmkdir -p /home/alon/runs/i040-gas-thermo-bench/verify\nPYTHONPATH=$PWD /home/alon/anaconda3/envs/rmg_env/bin/python test/rmgpy/kmc/fixtures/i040_probe/run_probe.py --scratch /home/alon/runs/i040-gas-thermo-bench/verify > >(tee -a /home/alon/runs/i040-gas-thermo-bench/verify/stdout.log) 2> >(tee -a /home/alon/runs/i040-gas-thermo-bench/verify/stderr.log >&2)\n```\n",
           "There is one committed script. Every invocation materializes the pinned allowlist via I034's `git show` helper, reestimates all species, recompiles the baseline, verifies reaction atom balance, compares exact net source vectors and thermo over the Cp grid, checks roots by bisection, Brent and RMG Kc, and compares the complete generated report byte for byte. It writes `results.json` under scratch. `--write-report` is the authoring mode; omit it for verification. Literature values are transcribed inputs in the script, with original citations below; rerunning verifies their arithmetic and model comparisons, not new measurements or a live website scrape.\n",
           f"Database `{prior.DATABASE_SHA}`: {d['snapshot']['files']} allowlisted files, SHA256 `{d['snapshot']['sha256']}`. No database mutation, rate-tree generation, prohibited dataset access, product-code changes, or QM runs.\n",
           "## State, structures, and baseline\n",
           f"All RMG tables use {T0:.2f} K, gas ideal standard pressure {P0:.0f} Pa. Concentration is {C0:.0f} mol/m³. H is formation enthalpy or reaction enthalpy in kJ/mol, S and Cp in J/mol/K. The Ince radical comparison alone uses its published {REFERENCES['1-phenylethyl radical']['T']:.2f} K anchor; its RMG comparison is evaluated at that same temperature.\n",
           f"Fresh compiled baseline: ΔH° = {b['HSCp'][0]/1000:.6f}; ΔS° = {b['HSCp'][1]:.6f}; ΔCp°({T0:.2f}) = {b['HSCp'][2]:.6f}. Continuous gas Tc = **{b['Tc']:.6f} K**, compiler-grid Tc = {b['grid_Tc']:.6f} K. The selected end radical is primary, not benzylic.\n",
           "```text\n" + " + ".join(b["reactants"]) + " -> " + " + ".join(b["products"]) + "\n```\n",
           "Tc solves `H°(T) - T S°(T) - RT ln(c0 RT/p0) = 0`, equivalently `Kc(T)c0 = 1`. Thus `Hc = H°+RT` and `Sc = S°+R[ln(c0 RT/p0)+1]` for this association. Gas entropy errors are never inserted into a liquid-monomer standard.\n",
           f"Implicit sensitivities at the uncorrected root are dTc/d(δH) = {b['sensitivity_H_kJ']:.6f} K/(kJ/mol), dTc/d(δS) = {b['sensitivity_S']:.6f} K/(J/mol/K); actual shifts below use nonlinear roots, not this linear approximation.\n",
           "## Main independent-reference table\n",
           "δ means reference minus RMG. Quoted uncertainties are source magnitudes, not RMG confidence intervals. A blank reference means no defensible independent value was retrieved in this probe. ATcT is the explicitly versioned archive retrieved, not asserted to be today's latest network.\n"]
    rows = []
    for name in SMILES:
        h, s, _ = d["species"][name]["HSCp"][0]
        if name in REFERENCES:
            ref, comparison = REFERENCES[name], d["comparisons"][name]
            mark = "†" if ref["T"] != T0 else "‡" if ref.get("conditional_derived") else ""
            rows.append([name + mark, f"{comparison['model_H']:.3f}", measured(ref['H'], ref['uH']),
                         f"{comparison['ref_minus_model_H']:+.3f}", f"{comparison['model_S']:.3f}",
                         measured(ref['S'], ref['uS']), f"{comparison['ref_minus_model_S']:+.3f}",
                         f"[{ref['Hsource']}]({SOURCES[ref['Hsource']][1]})/"
                         f"[{ref['Ssource']}]({SOURCES[ref['Ssource']][1]})"])
        else:
            rows.append([name, f"{h/1000:.3f}", "—", "unidentified", f"{s:.3f}", "—", "unidentified", "missing"])
    out += [table(["Species", "RMG H", "Reference H ± u", "δH", "RMG S", "Reference S ± u", "δS", "Sources"], rows), "\n",
            "† The displayed RMG radical row is at the reference anchor, whereas the Cp table below also exposes the standard RMG anchor. The radical enthalpy uses G4 with published bond-additivity correction; the raw G4 value is also retained below. No uncertainty is supplied for that individual calculation. Its molecular entropy, rather than the intrinsic group entropy, is used.\n",
            "V24's diphenylpropane result is corrected **G3MP2**, with a published expanded uncertainty (95%); it does not meet the requested G4/CBS/W1 tier. Its gas entropy has no supplied uncertainty. This entry is supporting evidence, not a substitute for a higher-level independent benchmark. The RMG benzylic correction's CBS-QB3 fit already includes phenylethyl chemistry; reusing that training fit as independent evidence would be circular. I17 provides a separate G4 calculation.\n",
            f"The radical raw G4 H is {REFERENCES['1-phenylethyl radical']['raw_G4_H']:.1f} kJ/mol, versus BAC H {REFERENCES['1-phenylethyl radical']['H']:.1f} kJ/mol. Their difference is a method variant, not a confidence interval.\n",
            "‡ Cumyl H/S are **derived conditionally**, not tabulated absolute values: `Hf(radical)=Hf(parent)+BDE−Hf(H)` and `S(radical)=S(parent)+(BDE−BDFE)/T−S(H)`. These use the independent NIST parent and atomic-H anchors with J21's gas CBS values, assuming their dissociation thermochemistry is at the nominal anchor with matching gas pressure conventions. The accessible main paper verifies that the BDE/BDFE difference largely reflects gas H-atom entropy, but its SI temperature/rotor details were not retrieved. Its molecule-specific method uncertainties are absent; the W1BD benchmark mean error is not a species uncertainty. These derived entropies do not establish a converged hindered-rotor/conformer entropy. [J21](https://pmc.ncbi.nlm.nih.gov/articles/PMC8343544/#tbl3), [N-H](https://webbook.nist.gov/cgi/cbook.cgi?ID=C12385136&Mask=1).\n",
            "The source entropy standard pressures are not explicit on every accessible compilation page. Published S values are treated as nominal one-bar values; this is an unresolved convention uncertainty for legacy entries. If every entry instead used one atmosphere, each S would increase by "
            f"{d['pressure_S_offset_per_species']:.6f} J/mol/K on conversion to one bar. The surrogate has Δν = −1, so its reference ΔS would decrease by that amount; mixed-source pressure conventions require source-specific resolution. No such normalization is silently applied.\n",
            "## Heat capacities\n",
            "Each cell is RMG Cp at the column temperature. NIST recommended Cp arrays are comparison inputs, not independently measured values with assigned uncertainties.\n",
            table(["Species"] + [f"{t:g} K" for t in TEMPS],
                  [[name] + [f"{v[2]:.3f}" for v in item["HSCp"]] for name, item in d["species"].items()]), "\n"]
    for name, ref in REFERENCES.items():
        if "Cp" not in ref:
            continue
        model = d["species"][name]
        vals = [float(ref["Cp"][j]) for j in range(len(ref["CpT"]))]
        # RMG Cp on extra reference knots is saved independently in main().
        out += [f"Reference Cp and δCp for **{name}** ({ref['Cpsource']}; uCp missing):\n",
                table(["T / K", "reference Cp", "RMG Cp", "δCp"],
                      [[f"{t:g}", f"{v:.3f}", f"{m:.3f}", f"{v-m:+.3f}"]
                       for t, v, m in zip(ref["CpT"], vals, model["reference_knot_Cp"])]), "\n"]
    out += ["No retrieved gas Cp reference for cumene, n-propylbenzene, the three diphenyl compounds or cumyl is represented as a zero error. These missing functions prevent a fully independent high-temperature propagation prediction.\n",
            "## Independent radical bond cycles\n",
            "J21's primary table was visually read from its archived image. The CBS columns below are calculated gas dissociation quantities; the Luo column is a compilation reproduced by the authors, without an uncertainty on these rows. It is retained as a conflicting comparison, not chosen by agreement with RMG. H/S derived from different parent anchors are not the internally consistent absolute CBS species values. All cycle values assume the nominal anchor; the entropy-derived column is conditional as above.\n",
            table(["Parent → radical + H", "RMG BDE / kJ", "CBS BDE / kJ", "RMG BDFE / kJ", "CBS BDFE / kJ", "Derived CBS radical H", "Derived CBS radical S", "Derived Luo radical H"],
                  [[parent, f"{row['model_BDE_kJ']:.6f}", f"{row['CBS_BDE_kJ']:.6f}",
                    f"{row['model_BDFE_kJ']:.6f}", f"{row['CBS_BDFE_kJ']:.6f}",
                    f"{row['derived_CBS_radical_H']:.6f}", f"{row['derived_CBS_radical_S']:.6f}",
                    f"{row['derived_Luo_radical_H']:.6f}"] for parent, row in d["bond_cycles"].items()]), "\n",
            "The cumyl enthalpy method spread is substantial and cannot be settled by assigning a fabricated uncertainty or by choosing whichever radical value agrees with a desired Tc. Benzyl_T still has coefficient zero in the compiled step. Its possible species-level error therefore does not shift that step when confined to this group; a correction to shared saturated groups remains a different inference.\n",
            "## Model reactions and group cancellation\n",
            "Reaction coefficients below are products positive, reactants negative. All are checked for atom balance. Only the two complete small-molecule reference reactions have independent numeric ΔH/ΔS. Missing product-radical or oligomer data prevents independent values for the direct addition reactions.\n"]
    rows = []
    for name, reaction in d["reactions"].items():
        ref = reaction.get("reference_HS")
        rows.append([name, f"{reaction['HSCp'][0]/1000:.6f}", f"{reaction['HSCp'][1]:.6f}",
                     measured(ref[0]/1000, reaction["reference_H_bound"]/1000) if ref else "missing: " + ", ".join(reaction["missing_reference"]),
                     f"{ref[1]:.6f} (u incomplete)" if ref else "unidentified"])
    out += [table(["Reaction", "RMG ΔH", "RMG ΔS", "Reference ΔH ± source-magnitude bound", "Reference ΔS"], rows), "\n"]
    for name, reaction in d["reactions"].items():
        out += [f"**{name}:** `" + "; ".join(f"{v:+d} {k}" for k, v in reaction["coefficients"].items()) + "`\n",
                "Net source vector: `" + json.dumps(reaction["groups"], sort_keys=True) + "`.\n"]
    out += [f"The balanced surrogate is **ethylbenzene + ethane + styrene → cumene + n-propylbenzene**. Its non-symmetry source vector, ΔH(T), and ΔCp(T) exactly equal the compiled increment on every reported temperature. Compiled ΔS minus surrogate ΔS is a temperature-independent {comb['entropy_offset']:+.6f} J/mol/K, explicitly retained. Group identity does not prove equality of real-molecule nonadditivity, conformer populations, or stereochemical ensembles.\n",
            "The two-ring nonadditivity reaction is **diphenylpropane + ethane → ethylbenzene + n-propylbenzene**. Its entire additive source vector vanishes, so its nonzero supporting reference enthalpy suggests physics invisible to those groups, subject to the lower-level diphenylpropane calculation. It cannot be assigned uniquely to a propagation group or used as a unique correction. Diphenylbutane/pentane stereochemistry is unspecified in these RMG graphs; RMG's local optical factors are not an independently established meso/racemic equilibrium ensemble.\n",
            "The smallest primary addition is explicitly checked separately: its radical environment changes and it need not reproduce the saturated long-chain increment. Benzylic self-similar propagation cancels the same Benzyl_S correction at both ends. The compiled primary propagation cancels Isobutyl at both ends. Consequently an error confined to either radical correction gives **δΔH = 0, δΔS = 0, δTc = 0** for its self-similar step. The cumyl Benzyl_T correction has coefficient zero in the compiled step. Whole-molecule radical discrepancies do not identify errors in the individual HBI group; library, saturated skeleton, symmetry and Cp contributions must be separated first.\n",
            "## Conditional effects on compiled propagation and Tc\n",
            "Rows change one published surrogate input at a time, keeping the other inputs and RMG ΔCp fixed. The styrene rows are direct one-species replacements in the actual compiled reaction. The other rows transfer through the proved additive identity and therefore require absence of additional oligomer nonadditivity. Their signs follow stoichiometry. They are controls, not independent additive uncertainties; do not add the combined row again. Tc intervals are nonlinear endpoint envelopes from the published input uncertainty magnitude alone. They exclude model error and omitted Cp/nonadditivity/stereo/pressure uncertainties. An unavailable bound is printed explicitly.\n"]
    rows = []
    for row in d["effects"]:
        bounds = row["interval"]
        rows.append([row["name"], row["axis"], f"{row['dH']/1000:+.6f}", f"{row['dS']:+.6f}",
                     f"{row['Tc']:.6f}", f"{row['Tc']-b['Tc']:+.6f}",
                     f"[{bounds[0]:.3f}, {bounds[1]:.3f}]" if bounds else "u missing"])
    out += [table(["Input replaced", "Axis", "δΔH", "δΔS", "Tc / K", "δTc / K", "Input-bound Tc / K"], rows), "\n",
            f"Combining the complete surrogate references gives δΔH = **{comb['dH']/1000:+.6f} kJ/mol**, δΔS = **{comb['dS']:+.6f} J/mol/K**, retaining the model symmetry offset. Enthalpy alone gives Tc {comb['Tc_H']:.6f} K; entropy alone {comb['Tc_S']:.6f} K; both {comb['Tc_HS']:.6f} K. This is a small-molecule transfer exercise with unchanged RMG ΔCp, not an independent oligomer Tc estimate.\n",
            f"The sum of published H uncertainty magnitudes is ±{comb['H_bound']/1000:.3f} kJ/mol, giving the enthalpy-only Tc envelope [{comb['Tc_H_bound'][0]:.3f}, {comb['Tc_H_bound'][1]:.3f}] K. ATcT H uncertainties use its network confidence convention; legacy NIST magnitudes are not all harmonized confidence intervals. This triangle-inequality envelope is conditional on interpreting each supplied magnitude as a bound, not a guaranteed statistical interval. The ATcT styrene–ethylbenzene correlation is +{.433:.3f}; with both coefficients negative its covariance contribution is positive. An illustrative RSS using that pair and treating all other sources as uncorrelated is {comb['illustrative_H_rss_kJ']:.6f} kJ/mol, **not** a defensible combined confidence interval. The S envelope is unavailable because cumene/n-propylbenzene lack reported uS.\n",
            "For the styrene and ethylbenzene Cp discrepancies, the code integrates piecewise-linear published Cp minus RMG Cp from the anchor, separately for δH(T) and δS(T), then applies their negative surrogate coefficients. Styrene transfers directly; ethylbenzene remains conditional on the group identity. This is additional to the constant anchor errors. Separate H-only or S-only Cp rows are accounting controls; the thermodynamically consistent correction includes both. No Cp uncertainty envelope can be supplied.\n",
            table(["Input: Cp contribution", "δΔH at baseline Tc / kJ/mol", "δΔS at baseline Tc", "Tc / K", "δTc / K"],
                  [[axis, f"{row['change_at_baseline_Tc'][0]/1000:+.6f}", f"{row['change_at_baseline_Tc'][1]:+.6f}",
                    f"{row['Tc']:.6f}", f"{row['Tc']-b['Tc']:+.6f}"] for axis, row in d["cp_effects"].items()]), "\n",
            f"H anchor + S anchor + Cp correction for styrene alone gives Tc {d['full_input_Tc']['styrene']:.6f} K; for ethylbenzene alone, via the surrogate, {d['full_input_Tc']['ethylbenzene']:.6f} K. Combining all surrogate anchor errors with only the known styrene/ethylbenzene Cp discrepancies gives {d['known_cp_surrogate_Tc']:.6f} K. All of these controls retain the other species' RMG Cp functions; a complete independent ΔCp remains unavailable.\n",
            "Phenylethyl's benchmark discrepancy can test Benzyl_S only after subtracting the independently assessed saturated parent and consistent entropy conventions. Even such an identified Benzyl_S correction cancels in self-similar benzylic propagation, and Benzyl_S is absent from the compiled primary propagation. Diphenylpropane's whole-molecule residual has no unique transfer coefficient; diphenylbutane, diphenylpentane and the actual oligomer radicals have no complete independent references here. Cumyl has a conditional cycle reference with unresolved method/entropy uncertainties and no unique shared-group attribution. For those whole-molecule residuals δΔH, δΔS and δTc are **unidentified**, rather than set to zero.\n",
            "## Exact species provenance\n",
            "These are the actual comments, source weights, aliases, graph symmetries and optical-site counts used by RMG. Each H/S decomposition is checked against the estimated species at the I039 decomposition temperatures. Full source entry descriptions and quantity uncertainty strings are also saved in scratch results.json.\n"]
    for name, row in d["species"].items():
        out += [f"**{name}** — `{row['smiles']}`; σ = {row['symmetry_number_including_optical']:g}; optical half-factor sites = {row['optical_atom_half_factors']}.\n",
                "```text\n" + row["thermo_comment"] + "\n" + json.dumps(row["source_weights"], sort_keys=True) + "\n```\n"]
    metadata = {}
    for item in d["species"].values():
        metadata.update(item["source_metadata"])
    out += ["Resolved source entries (quoted quantities are database inputs, not newly estimated error bars):\n",
            table(["Matched entry", "Resolved data entry", "H298", "S298", "Attribution"],
                  [[label, " → ".join(meta["resolved_aliases"]), meta["H298"], meta["S298"],
                    meta["shortDesc"].replace("\n", " ").replace("|", "/")]
                   for label, meta in sorted(metadata.items())]), "\n",
            "Shared Benson/Stein group input uncertainties and radical fit method estimates have unknown covariance and validation scope. They cannot be treated as independent species uncertainties or combined with the external errors as if each molecule were an independent measurement. Common groups cancel before any error propagation.\n",
            "## Literature gaps and costed QM plan — description only\n",
            "First calculate the actual gas compounds and reaction increments, without using any ceiling-temperature target for selection or calibration. Small-species set: styrene, ethylbenzene, cumene, n-propylbenzene, ethane, 1-phenylethyl radical and cumyl radical. Large-species set: the three diphenyl compounds, P1, P2, P3 primary radicals and B2 benzylic radical. These cover the requested species, both smallest additions and the actual compiled increment. B3 can be added if benzylic length convergence is needed; it is not budgeted initially.\n",
            "Use CREST/GFN2-xTB for conformer discovery (both doublet and singlet as appropriate), followed by dispersion-aware ωB97X-D/def2-TZVP optimization/frequencies and connectivity checks. Enumerate stereoisomers before conformers, explicitly including meso and racemic diphenylpentane; report each stereoisomer and a defined ensemble, rather than silently averaging across fixed tacticities. Include electronic spin degeneracy, standard pressure, rotational symmetry, and enantiomer degeneracy consistently. Deduplicate minima, retain low-energy families and expand the energy window until ensemble H/S/Cp converge. Initial conformer quotas in the budget are starting estimates, not proof of convergence.\n",
            "Use canonical G4 for the small set, with raw atomization and optional BAC results reported separately. For the large set use TightPNO DLPNO-CCSD(T1) single points at def2-TZVPP/QZVPP with documented CBS extrapolation, plus canonical G4 cross-checks on representative large conformers. Check unrestricted spin contamination and open-shell diagnostics; if multi-reference character or PNO/basis sensitivity exceeds the planned uncertainty budget, escalate the affected pair rather than forcing an increment. Reaction energies should be computed directly from consistently treated species and isodesmic cycles; external ATcT/NIST anchors test method bias.\n",
            "Treat methyl, phenyl and backbone torsions with relaxed periodic scans at the DFT level and replace corresponding harmonic modes by hindered-rotor partition functions. Test coupled backbone/phenyl torsions with two-dimensional grids on representative minima, and avoid double counting conformers and rotor states. Compute Boltzmann-weighted H/S/Cp over the relevant gas-temperature range. Convergence targets: sub-kJ/mol reaction-enthalpy changes and about one J/mol/K reaction-entropy changes under conformer/rotor/basis expansion; these are proposed numerical goals, not achieved uncertainties. Use the observed convergence and cross-method/cycle residuals for a covariance-aware uncertainty budget.\n",
            "Budget below is **estimated serial core-hours**, not measured wall time, and includes no QM executed in this probe. A single core-hour is one CPU core for one hour; parallel speedup/memory limits need a pilot. The costed subset has finite scan/conformer quotas, so the budget is conditional on convergence.\n",
            table(["Planned task", "Jobs", "Estimated core-h/job", "Estimated core-h"],
                  [[row["task"], row["jobs"], f"{row['core_hours_per_job']:g}", f"{row['core_hours']:g}"] for row in d["costs"]]), "\n"]
    total = sum(row["core_hours"] for row in d["costs"])
    out += [f"Estimated initial total: **{total:,.0f} core-hours**. A factor-of-three planning range is {total/3:,.0f}–{total*3:,.0f} core-hours. At ideal use of {32} cores this corresponds to {total/32:.1f} aggregate node-hours; it is not an elapsed-time prediction. Run a small fixed pilot to replace unit costs before authorizing this QM campaign. No QM run is authorized by this report.\n",
            "## Limits of the dispatch premise\n",
            "The stated gas ceiling is reproduced, but a request to propagate *each molecular benchmark error* uniquely into propagation overstates identifiability. Shared groups, cancelling radical corrections, unknown two-ring nonadditivity and unspecified stereochemistry prevent that inference. The tables show direct replacements, explicitly conditional surrogate transfers, exact zero coefficients, and unresolved effects separately. No complete independent reference exists in the retrieved evidence for the actual P2/P3 reaction. Its total uncertainty and an independently validated gas Tc remain unsettled.\n",
            "A quoted liquid-monomer ceiling is not an independent validation of an ideal-gas molecular benchmark. The source-pressure ambiguity and legacy uncertainty conventions also prevent calling all supplied S/H values a common confidence interval. None of these limitations requires changing the task's implementation scope; the report and verifier are complete, while the scientific gaps require the described additional calculations or better references.\n",
            "## Primary references\n"]
    for key, (title, url) in SOURCES.items():
        out.append(f"- **{key}**: [{title}]({url}).\n")
    out += ["\nV24 accessible primary-paper copy: [PDF](https://pdfs.semanticscholar.org/d7d0/6da8ff5248dc2515cbc4ccd8b8cc19c9d56b.pdf). Its gas H is corrected G3MP2 and the gas S uncertainty is not the uncertainty of vaporization entropy. Ince's group-fit residual statistics are not uncertainty bounds on a specific G4 species.\n",
            "Legacy diphenylpropane combustion compilations contain conflicting derived liquid values and lack a complete vetted gas conversion/uncertainty chain in the retrieved record; they were not substituted for a gas reference. Missing references are retrieval/validation gaps in this probe, not assertions that no reference can exist. Pyrolysis results and prohibited datasets were not used.\n"]
    return "\n".join(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--database", type=Path, default=Path("/home/alon/Code/RMG-database"))
    parser.add_argument("--write-report", action="store_true")
    args = parser.parse_args()
    args.scratch = args.scratch.resolve()
    if args.database.resolve() != Path("/home/alon/Code/RMG-database").resolve():
        raise ValueError("only the dispatch-named read-only database may be accessed")
    allowed = Path("/home/alon/runs/i040-gas-thermo-bench")
    if not args.scratch.is_relative_to(allowed) or "catalog" in args.scratch.parts:
        raise ValueError("scratch must be under dispatch's allowed run directory")
    args.scratch.mkdir(parents=True, exist_ok=True)
    snapshot = prior.snapshot_database(args.database, args.scratch / "database")
    database = RMGDatabase()
    database.load_kinetics(str(args.scratch / "database/input/kinetics"), reaction_libraries=[],
                           seed_mechanisms=None, kinetics_families=list(prior.PS_FAMILY_CANDIDATES),
                           kinetics_depositories=["training"])
    database.load_thermo(str(args.scratch / "database/input/thermo"),
                         thermo_libraries=["primaryThermoLibrary"], depository=True)
    results = calculate(database)
    for name, ref in REFERENCES.items():
        if "CpT" in ref:
            species = prior.Species(molecule=[prior.Molecule().from_adjacency_list(results["species"][name]["adjacency"])])
            species.generate_resonance_structures()
            species.thermo = database.thermo.get_thermo_data(species)
            results["species"][name]["reference_knot_Cp"] = [thermo(species, t)[2] for t in ref["CpT"]]
    results["snapshot"] = snapshot
    results["inputs"] = REFERENCES
    results["bond_inputs"] = BOND_REFERENCES
    results["H_atom_input"] = H_ATOM
    results["sources"] = SOURCES
    results["script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (args.scratch / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    report = render(results)
    target = FIXTURES / "I040_gas_thermo_bench.md"
    if args.write_report:
        target.write_text(report)
        print("WROTE " + str(target))
    else:
        if target.read_text() != report:
            raise AssertionError("report differs from fresh independently reproduced calculation")
        print("PASS: pinned snapshot, species H/S/Cp and sources, atom balances, source-vector identities,")
        print("      compiler baseline, nonlinear roots by three formulations, reference arithmetic,")
        print("      separate H/S/Cp transfers, uncertainty envelopes, QM budget, entire report byte comparison")
    print(f"baseline Tc={results['baseline']['Tc']:.6f} K; "
          f"surrogate dH={results['combined']['dH']/1000:+.6f} kJ/mol; "
          f"dS={results['combined']['dS']:+.6f} J/mol/K")


if __name__ == "__main__":
    main()
