#!/usr/bin/env python3
"""I-065 check 4b: run the ACCEPTED-state census over real DASPK trajectories
and report how often each half of the realizability check fires.

The census hook (_assert_pool_moments_accepted, reached from
_phase_gate_flux_census) is what base.pyx calls once per accepted snapshot
during simulate(). This replays that hook against the DASPK integrator on the
production decks' own pool parameters and t=0 moments, plus a sweep around
them, which reaches thousands of accepted states in seconds where a full RMG
generation reaches them in hours.

It is a COMPLEMENT to running the decks, not a substitute: it exercises the
polymer moment subsystem with the deck's kinetics, not the deck's full
generated mechanism. Prints a count per arm and exits 0 always -- the caller
decides what a nonzero count means.
"""

import contextlib
import io
import json
import logging
import sys

import numpy as np

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig

MU0, MU1, MU2 = 3, 4, 5


def _spc(smiles, label):
    s = Species(molecule=[Molecule().from_smiles(smiles)])
    s.label = label
    return s


class _Grab(logging.Handler):
    def __init__(self):
        super().__init__()
        self.cone = []
        self.var = []

    def emit(self, record):
        m = record.getMessage()
        if "CONE CENSUS" in m:
            self.cone.append(m)
        elif "VARIANCE CENSUS" in m:
            self.var.append(m)


def _system(moments, k_unzip, k_scission, T, monomer_mw=42.08):
    inert = _spc("N#N", "N2")
    core = [inert, _spc("C=CC", "monomer_gas"),
            _spc("[CH2]CC", "poly"), _spc("CCO", "poly_mu0"),
            _spc("CC=O", "poly_mu1"), _spc("CC#N", "poly_mu2")]
    mask = np.array([s.label in ("N2", "monomer_gas") for s in core], dtype=bool)
    pools = [PolymerPoolConfig(
        label="poly", xs=2, explicit_dp_to_species_index={},
        mu_indices=(MU0, MU1, MU2), monomer_poly_index=1,
        monomer_mw_g_mol=monomer_mw, k_scission=k_scission, k_unzip=k_unzip,
        tail_kinetics=None)]
    rs = HybridPolymerSystem(
        T=T, P=1.0e5, initial_mole_fractions={inert: 1.0}, V_poly=1.0,
        polymer_pools=pools, mass_transfer=[], gas_species_mask=mask,
        constant_gas_volume=False,
        initial_polymer_moments={"poly": tuple(moments)}, termination=[])
    with contextlib.redirect_stdout(io.StringIO()):
        rs.initialize_model(core, [], [], [])
    return rs


def run_arm(name, moments, k_unzip, k_scission, T, t_end, n_steps=400):
    """Step DASPK to t_end and invoke the accepted-state census at every
    accepted snapshot, exactly as base.pyx's simulate() does."""
    rs = _system(moments, k_unzip, k_scission, T)
    grab = _Grab()
    root = logging.getLogger()
    root.addHandler(grab)
    old = root.level
    root.setLevel(logging.WARNING)
    accepted = 0
    worst_gap = float("inf")
    worst_var = float("inf")
    err = None
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            for i in range(1, n_steps + 1):
                t = t_end * i / n_steps
                try:
                    rs.advance(t)
                except Exception as exc:            # integrator gave up
                    err = f"{type(exc).__name__}: {exc}"[:160]
                    break
                accepted += 1
                m0, m1, m2 = rs.y[MU0], rs.y[MU1], rs.y[MU2]
                worst_gap = min(worst_gap, float(m1 - m0))
                worst_var = min(worst_var, float(m0 * m2 - m1 * m1))
                try:
                    rs._phase_gate_flux_census([], [], [], 1.0, 0.0)
                except ValueError as exc:       # the r81 accepted-state raise
                    err = f"r81 raise at t={t:g}: {exc}"[:200]
                    break
    finally:
        root.removeHandler(grab)
        root.setLevel(old)
    return {"arm": name, "moments": list(moments), "k_unzip": k_unzip,
            "k_scission": k_scission, "T": T, "t_end": t_end,
            "accepted_states": accepted, "cone_hits": len(grab.cone),
            "variance_hits": len(grab.var),
            "min_mu1_minus_mu0": worst_gap,
            "min_mu0mu2_minus_mu1sq": worst_var,
            "final": [float(rs.y[MU0]), float(rs.y[MU1]), float(rs.y[MU2])],
            "error": err}


# poly_104/poly_105 (phenol-formaldehyde, 1100 K): mu0 = 0.099 mol chains,
# Mn = 1000 g/mol over the C8H9O repeat unit -> the moments below; the deck
# runs k_unzip = 100, k_scission = 1, terminationTime = 100 s.
POLY105 = (0.099, 0.7378430665, 16.4973451743)
ARMS = [
    ("poly_105 (k_u=100, k_s=1, 1100 K)", POLY105, 100.0, 1.0, 1100.0, 100.0),
    ("poly_105 short window", POLY105, 100.0, 1.0, 1100.0, 1.0),
    # examples/rmg/polystyrene: initialMoles PS = 0.01 mol chains, Mn = 5000,
    # Mw = 6000 over the styrene repeat unit (104.15 g/mol) ->
    # DPn = 48.007, DPw = 57.609, mu1 = mu0*DPn, mu2 = mu1*DPw.
    ("polystyrene example (k_u=0, k_s=0)",
     (0.01, 0.01 * 5000.0 / 104.15, 0.01 * 5000.0 / 104.15 * 6000.0 / 104.15),
     0.0, 0.0, 1000.0, 0.1),
    # sweep around the production point: the k_unzip/k_scission ratio is what
    # sets whether the pool reaches the cone boundary at all (I-055 measured
    # the negative-mu1 threshold at k_scission/4)
    ("sweep k_u=0.1 k_s=1", POLY105, 0.1, 1.0, 1100.0, 100.0),
    ("sweep k_u=0.25 k_s=1", POLY105, 0.25, 1.0, 1100.0, 100.0),
    ("sweep k_u=1 k_s=1", POLY105, 1.0, 1.0, 1100.0, 100.0),
    ("sweep k_u=10 k_s=0", POLY105, 10.0, 0.0, 1100.0, 100.0),
    ("sweep k_u=1000 k_s=1", POLY105, 1000.0, 1.0, 1100.0, 10.0),
]


def main():
    # DASPK writes its own diagnostics to fd 1 from Fortran, underneath
    # Python's redirect_stdout, so the JSON goes to a named file rather than
    # to stdout -- otherwise the two interleave and the result is unparseable.
    if len(sys.argv) != 2:
        sys.stderr.write("usage: i065_census_replay.py <out.json>\n")
        return 2
    out = [run_arm(*a) for a in ARMS]
    with open(sys.argv[1], "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
