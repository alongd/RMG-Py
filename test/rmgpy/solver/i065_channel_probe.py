#!/usr/bin/env python3
"""I-065 two-arm channel probe: the chain-end RELEASE derivatives, printed as
hex floats so the fixed and unfixed builds can be differenced bit-for-bit.

Runs identically on BOTH builds -- it imports nothing that exists only on one
of them -- so the SAME file measures the defect (unfixed) and its absence
(fixed). Emits one JSON object on stdout:

    {"build": <path to the loaded polymer extension>,
     "rows": [{"channel": ..., "mu": [...], "dmu0": hex, "dmu1": hex,
               "dmu2": hex, "gas": hex}, ...]}

Every number is recomputed from the real HybridPolymerSystem residual; nothing
is restated.
"""

import contextlib
import io
import json
import math
import sys

import numpy as np

from rmgpy.molecule import Molecule
from rmgpy.species import Species
import rmgpy.solver.polymer as _polymer_mod
from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig

GAS, MU0, MU1, MU2 = 1, 3, 4, 5
R_GAS = 8.314
KDEP_TRIPLET = dict(A=1.0e13, n=0.5, Ea=1.2e5)
QSSA_CHANNEL = dict(
    initiation=dict(A=1.0e13, n=0.0, Ea=1.5e5),
    depropagation=dict(A=1.0e13, n=0.0, Ea=1.0e5),
    termination=dict(A=1.0e8, n=0.0, Ea=0.0),
    efficiency=1.0, monomer_yield=1.0)

# States: two inside the realizable cone (mu1 >= mu0 >= 0, where the fix must
# be a bit-for-bit no-op) and three outside it (where it must stop the
# fabrication). (0.099, 0.7378, 16.4973) is poly_104/poly_105's own t=0 state.
CONE_STATES = [(1.0, 5.0, 30.0), (0.4, 0.4, 0.4), (0.25, 0.3, 0.5),
               (0.099, 0.7378430665, 16.4973451743), (1.0, 1.0, 1.0)]
OFF_CONE_STATES = [(0.5, 0.0, 0.0), (0.5, 0.2, 0.15), (0.116, 0.0139, 0.05)]

# Randomized sweep of FULLY realizable states -- mu1 >= mu0 > 0 AND
# mu0*mu2 >= mu1^2, i.e. both halves of the cone, spanning eight decades of
# mu0 and mean DP from 1 to ~1e4. Fixed seed so the verifier is reproducible.
# Five hand-picked states can be a lucky sample; this cannot.
def _random_cone_states(n=150, seed=20260903):
    rng = np.random.default_rng(seed)
    out = []
    while len(out) < n:
        mu0 = float(10.0 ** rng.uniform(-6.0, 2.0))
        mean = float(10.0 ** rng.uniform(0.0, 4.0))       # DP >= 1
        pdi = float(1.0 + 10.0 ** rng.uniform(-4.0, 1.0))  # PDI > 1
        mu1 = mu0 * mean
        mu2 = mu1 * mean * pdi
        if mu1 >= mu0 > 0.0 and mu0 * mu2 >= mu1 * mu1 and np.isfinite(mu2):
            out.append((mu0, mu1, mu2))
    return out


def _spc(smiles, label):
    s = Species(molecule=[Molecule().from_smiles(smiles)])
    s.label = label
    return s


def _system(moments, k_unzip=0.0, k_depropagation=None, radical_qssa=None,
            T=800.0):
    inert = _spc("N#N", "N2")
    core = [inert, _spc("C=CC", "monomer_gas"),
            _spc("[CH2]CC", "poly"), _spc("CCO", "poly_mu0"),
            _spc("CC=O", "poly_mu1"), _spc("CC#N", "poly_mu2")]
    mask = np.array([s.label in ("N2", "monomer_gas") for s in core], dtype=bool)
    pools = [PolymerPoolConfig(
        label="poly", xs=2, explicit_dp_to_species_index={},
        mu_indices=(MU0, MU1, MU2), monomer_poly_index=GAS,
        monomer_mw_g_mol=42.08, k_scission=0.0, k_unzip=k_unzip,
        radical_qssa_unzip=radical_qssa, k_depropagation=k_depropagation,
        tail_kinetics=None)]
    rs = HybridPolymerSystem(
        T=T, P=1.0e5, initial_mole_fractions={inert: 1.0}, V_poly=1.0,
        polymer_pools=pools, mass_transfer=[], gas_species_mask=mask,
        constant_gas_volume=False,
        initial_polymer_moments={"poly": tuple(moments)}, termination=[])
    with contextlib.redirect_stdout(io.StringIO()):
        rs.initialize_model(core, [], [], [])
    return rs


def _row(channel, mu, **kw):
    rs = _system(mu, **kw)
    dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
    return {"channel": channel, "mu": list(mu),
            "dmu0": float(dn[MU0]).hex(), "dmu1": float(dn[MU1]).hex(),
            "dmu2": float(dn[MU2]).hex(), "gas": float(dn[GAS]).hex()}


def kdep_k(T=800.0):
    return (KDEP_TRIPLET["A"] * T ** KDEP_TRIPLET["n"]
            * math.exp(-KDEP_TRIPLET["Ea"] / (R_GAS * T)))


def main():
    rows = []
    cone = CONE_STATES + _random_cone_states()
    for mu in cone + OFF_CONE_STATES:
        rows.append(_row("legacy_k_unzip", mu, k_unzip=0.7))
        rows.append(_row("k_depropagation", mu, k_depropagation=KDEP_TRIPLET))
        rows.append(_row("radical_qssa_unzip", mu, radical_qssa=QSSA_CHANNEL))
    json.dump({"build": _polymer_mod.__file__,
               "k_dep_800K": kdep_k(),
               "cone_states": [list(s) for s in cone],
               "off_cone_states": [list(s) for s in OFF_CONE_STATES],
               "rows": rows}, sys.stdout)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
