"""Cross-method rate identity of the catalog pool kernels: MoM's events
against the kinetic Monte Carlo's propensities on the same explicit states.

The rate-identity tests of polymerKernelsTest sum over chain pairs with the
solver's own ordered-pair normalization, so a convention both share (a
missing unordered 1/2, a reverse coefficient off its Kc) cannot fail them.
These take the other side from the kMC code instead: the fixture
``fixtures/i230_cross_method_propensities.json`` is written by
``test/rmgpy/kmc/momCrossMethodTest.py`` of the kMC checkout (which pins the
same sha256 and regenerates it): the melt engine's events per second on
explicit chains of the in-tree PE catalog at 700 K, the MoM terms its
producer writes for those kernels (volume ratio 1, V_poly = the engine's
volume), and Kc(700 K) of the D1 record's thermo.

Each case puts its chains (count / N_A mol of chains of n periods) in MoM's
pools and integrates nothing.  MoM's events per second of the case's term
are its term law's on those chains (the producer's record, every case), and
where the chains lie past each participant's exact-from length they are the
engine's mean field: per template k S0 S1 / N_A V when both participants
come from one pool -- the engine's k (S0 S1 - X) / N_A V plus the X
same-molecule site pairs it excludes (recombination and D1 of two radicals:
exactly 2) -- and the engine's events with two pools.  A radical of 4
periods (``d1_two_template_short``) lies below the D1 law's exact length:
there MoM integrates its law, not the mean field (the law-domain error).
"""

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest

from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig
from rmgpy.species import Species
from rmgpy.molecule import Molecule

FIXTURE = Path(__file__).parent / "fixtures" / "i230_cross_method_propensities.json"
CROSS_METHOD_SHA256 = "761840f5562848876eee3fcd9effa111e1623c54975dcc18af18fb7e63f8ebd4"
POOLS = ("P", "Re", "Ut", "Rm", "Ui")


@pytest.fixture(scope="module")
def doc():
    assert hashlib.sha256(FIXTURE.read_bytes()).hexdigest() == CROSS_METHOD_SHA256
    return json.loads(FIXTURE.read_text())


def _spc(smiles, label):
    return Species(label=label, molecule=[Molecule(smiles=smiles)])


def _term(doc, template):
    for kind in ("recombination",):
        for t in doc["terms"][kind]:
            if t["template"] == template:
                return kind, "Re", t
    for owner, ts in doc["terms"]["transfer"].items():
        for t in ts:
            if t["template"] == template:
                return "transfer", owner, t
    raise KeyError(template)


def mom_events(doc, case):
    """MoM's events per second of the case's term alone on its chains."""
    kind, owner, term = _term(doc, case["term"])
    ppu, na, volume = doc["periods_per_unit"], doc["avogadro"], doc["volume_m3"]
    core = [_spc("N#N", "N2")]
    pools = []
    for i, label in enumerate(POOLS):
        stem = "C" * (i + 1)
        core += [_spc("C" * (i + 4), label), _spc(stem + "O", f"{label}_mu0"), _spc(stem + "N", f"{label}_mu1"),
                 _spc(stem + "S", f"{label}_mu2")]
        kernels = {}
        if label == owner:
            kernels = {"kernel_recombination": [term]} if kind == "recombination" else {"kernel_transfer": [term]}
        pools.append(PolymerPoolConfig(label=label, xs=2, explicit_dp_to_species_index={},
                                       mu_indices=(2 + 4 * i, 3 + 4 * i, 4 + 4 * i), **kernels))
    moments = {label: tuple(math.fsum(count / na * (int(n) / ppu) ** k
                                      for n, count in case["chains"].get(label, {}).items()) for k in range(3))
               for label in POOLS}
    mask = np.array([s.label == "N2" for s in core], dtype=bool)
    rs = HybridPolymerSystem(T=doc["temperature_K"], P=1.0e5, initial_mole_fractions={core[0]: 1.0},
                             V_poly=volume, polymer_pools=pools, mass_transfer=[], gas_species_mask=mask,
                             constant_gas_volume=False, initial_polymer_moments=moments, termination=[])
    rs.initialize_model(core, [], [], [])
    dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
    i = POOLS.index(owner)
    lost = -dn[2 + 4 * i]                                # chains of the owner pool per second (mol/s)
    both_here = kind == "recombination" or term["partner"] == owner
    return (lost / 2.0 if both_here else lost) * na


def _ids(doc_cases):
    return [f"{c['name']}:{c['term'][:12]}" for c in doc_cases]


CASE_NAMES = ("recombination_like", "recombination_unlike", "d1_like", "d1_unlike", "reverse_d1", "h_transfer",
              "d1_two_template_short")


@pytest.mark.parametrize("name", CASE_NAMES)
def test_moms_events_are_its_term_laws_on_the_same_chains(doc, name):
    cases = [c for c in doc["cases"] if c["name"] == name]
    assert cases
    for case in cases:
        assert mom_events(doc, case) == pytest.approx(case["mom_law_events_per_s"], rel=1e-9), case["term"]


@pytest.mark.parametrize("name", CASE_NAMES)
def test_moms_events_are_the_melt_engines_mean_field_where_the_law_is_exact(doc, name):
    for case in [c for c in doc["cases"] if c["name"] == name]:
        if case["law_exact"]:
            assert mom_events(doc, case) == pytest.approx(case["mean_field_events_per_s"], rel=1e-9), case["term"]
        else:
            assert mom_events(doc, case) != pytest.approx(case["mean_field_events_per_s"], rel=1e-6)


def test_a_pair_from_one_pool_differs_from_the_engine_by_its_same_molecule_pairs_only(doc):
    """Two radicals: the engine has one unordered pair, the mean field four
    ordered site pairs of which two are a radical with itself -- the factor 2
    is the finite system, not a convention (it vanishes as 1/N); summed per
    template, with each template's own weight and divisor."""
    for case in doc["cases"]:
        if case["name"] in ("recombination_like", "d1_like"):
            term = [t for t in case["templates"] if t["in_term"]]
            mean_field = sum(t["k"] * t["site_pairs"] for t in term)
            engine = sum(t["k"] * (t["site_pairs"] - t["same_molecule_pairs"]) for t in term)
            assert mean_field / engine == 2.0


def test_the_forward_over_reverse_rate_coefficient_is_kc(doc):
    """D1 (R + R -> P + U) and its thermo-derived reverse: MoM's events per
    volume over [R]^2 and over [U][P] are k_f and k_r, and k_f / k_r is Kc of
    the D1 record's thermo at C0 = P0 / R T (the kMC's single convention)."""
    db = doc["detailed_balance"]
    na, volume = doc["avogadro"], doc["volume_m3"]
    (fwd,) = [c for c in doc["cases"] if c["name"] == "d1_like" and c["term"] == db["forward_term"]]
    (rev,) = [c for c in doc["cases"] if c["name"] == "reverse_d1" and c["term"] == db["reverse_term"]]

    def conc(case, pool):
        return sum(case["chains"][pool].values()) / na / volume
    kf = mom_events(doc, fwd) / na / volume / conc(fwd, "Re") ** 2
    kr = mom_events(doc, rev) / na / volume / (conc(rev, "Ut") * conc(rev, "P"))
    assert kf / kr == pytest.approx(db["kc_forward"], rel=1e-9)
