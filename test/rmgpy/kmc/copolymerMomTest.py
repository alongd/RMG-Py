"""Red tests (I-240, slice S5 of the copolymer design): the method of moments
for a linear statistical copolymer as one pool.

The pool's unit is one residue, of mean molar mass Mbar0 = sum x_i M_i
(exact rationals); its moments count residues.  Its kernels come from the
kMC catalog's templates, each template's sites at a residue read on that
residue's window of the chain and weighted by the window's frequency in the
declared sequence statistics (``word_frequencies``).  Chemistry whose rate
depends on which residues are there makes the composition drift, which one
pool cannot follow: it is reported by name (``limited:composition_drift``),
never corrected.

The chemistry is synthetic: the E/P test declaration with methyl end caps
(so every chain end looks alike to a four-atom window) and generated toy
homolysis templates, one per four-atom backbone window (read either way),
cutting the window's middle bond.  With one weight for every window the
catalog cuts every backbone bond whose window holds no cap, whatever the
residues (non-selective); a weight only on the windows whose cut bond
touches P's CH is P-only scission (residue-selective).  The references are
exact at the sequence level: every sequence of the initial length,
materialized, its sites counted by the melt engine's own index.
"""

import copy
import dataclasses
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy
import pytest

import melt_engine_support as S
from meltRunTest import toy_payloads, write_catalog
from meltSequenceTest import write_v2_deck
from polymer_spec_v2_support import DECKS, polymer_keyword, q

from rmgpy.kmc import kernel_rates as K
from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import melt_sequence as msq
from rmgpy.kmc import melt_state as ms
from rmgpy.kmc import run as kr
from rmgpy.kmc.window_matcher import AtomEdit, BondEdit, Template, WindowAtom
from rmgpy.polymer_spec import PolymerSpec

T = 700.0
W = 0.1
COMPOSITION = {"E": q(2, 3), "P": q(1, 3)}
# a Markov chain whose stationary distribution is (2/3, 1/3)
STICKY = {"model": "markov1", "transition": {"E": {"E": q(3, 4), "P": q(1, 4)}, "P": {"E": q(1, 2), "P": q(1, 2)}}}
METHYL_CAP = {"atoms": [{"key": "m", "element": "C", "hydrogens": 3, "electronic_state": {"radical_electrons": 0}}],
              "bonds": [], "ports": {"cap": {"atom": "m", "bond_order": 1}},
              "electronic_state": {"multiplicity": 1}}
HOMOLYSIS = "homolysis:substrate=chain"


def copolymer_pool():
    from rmgpy import copolymer_pool
    return copolymer_pool


def copolymer_kernels():
    from rmgpy.kmc import copolymer_kernels
    return copolymer_kernels


# --------------------------------------------------------------------------
# the declaration and the synthetic catalog
# --------------------------------------------------------------------------
def ep_mapping(counts=None, *, statistics=None, caps=True):
    mapping = copy.deepcopy(polymer_keyword(DECKS / "ep_v2_deck.py", "polymer_spec"))
    ensemble = mapping["initial_ensemble"]
    ensemble["dp_counts"] = dict(counts or {8: 30})
    ensemble["composition"] = dict(COMPOSITION)
    ensemble["composition_relative_tolerance"] = q(1, 1)
    if statistics is not None:
        ensemble["sequence_statistics"] = statistics
    mapping["transport_model"] = {"met": S.pe_met()}
    if caps:
        mapping["end_groups"] = {"head": copy.deepcopy(METHYL_CAP), "tail": copy.deepcopy(METHYL_CAP)}
    return mapping


def ep_spec(counts=None, **kw):
    return PolymerSpec.from_mapping(ep_mapping(counts, **kw))


CH2_PORT = ((1.0, "C", False), (1.0, "H", False), (1.0, "H", False))
CH_PORT = ((1.0, "C", False), (1.0, "C", False), (1.0, "H", False))
METHYL_PORT = ((1.0, "H", False),) * 3


def _core(kind):
    if kind == "CH2":
        return [WindowAtom("C"), WindowAtom("H"), WindowAtom("H")]
    return [WindowAtom("C"), WindowAtom("H"), WindowAtom("C", 0, "port", METHYL_PORT)]


def window_cut(left, b, c, right, weight):
    """Homolysis of the middle bond of a backbone window left-b-c-right (CH2
    or CH(CH3) carbons; the ends are ports)."""
    atoms = [WindowAtom("C", 0, "port", CH2_PORT if left == "CH2" else CH_PORT)]
    bonds = []
    ib = len(atoms)
    atoms += _core(b)
    bonds += [(ib, ib + 1, 1), (ib, ib + 2, 1), (0, ib, 1)]
    ic = len(atoms)
    atoms += _core(c)
    bonds += [(ic, ic + 1, 1), (ic, ic + 2, 1), (ib, ic, 1)]
    atoms.append(WindowAtom("C", 0, "port", CH2_PORT if right == "CH2" else CH_PORT))
    bonds.append((ic, len(atoms) - 1, 1))
    return Template(f"cut_{left}_{b}_{c}_{right}", "toy_homolysis", atoms, bonds,
                    [BondEdit(ib, ic, 0), AtomEdit(ib, radical=1), AtomEdit(ic, radical=1)], weight=weight)


def cut_templates(every=W, p_extra=0.0):
    """One template per backbone window read either way: weight ``every``,
    plus ``p_extra`` when the cut bond touches P's CH."""
    out = []
    for left, b, c, right in itertools.product(("CH2", "CH"), repeat=4):
        if (b, c) == ("CH", "CH") or (right, c, b, left) < (left, b, c, right):
            continue
        weight = every + (p_extra if "CH" in (b, c) else 0.0)
        if weight:
            out.append(window_cut(left, b, c, right, weight))
    return out


def catalog_doc(tmp_path, spec, templates, name="ep_cut.json"):
    path = write_catalog(tmp_path / name, toy_payloads(templates), tag="ep-cut")
    doc = json.loads(path.read_text())
    doc["provenance"]["chemistry_sha256"] = kr.chemistry_sha256(spec)
    path.write_text(json.dumps(doc, sort_keys=True, indent=1) + "\n")
    return doc, path


class _Catalog:
    def __init__(self, templates):
        self.ts = list(templates)

    def templates(self, _temperature):
        return list(self.ts)


# --------------------------------------------------------------------------
# the exact sequence-level reference
# --------------------------------------------------------------------------
def every_sequence(spec, n):
    """(sequence, probability) of every sequence of n residues: pi(s1) T(s1, s2) ..."""
    x = dict(spec.initial_ensemble.composition)
    stats = spec.initial_ensemble.statistics
    rows = ({i: dict(x) for i in x} if stats.model == "ideal_random" else
            {i: dict(row) for i, row in stats.transition})
    out = []
    for seq in itertools.product(sorted(r for r, v in x.items() if v), repeat=n):
        p = x[seq[0]]
        for a, b in zip(seq, seq[1:]):
            p *= rows[a].get(b, Fraction(0))
        if p:
            out.append((seq, p))
    return out


def backbone_path(state, cid):
    """The chain's backbone atoms from the head end."""
    origin = state.origin
    bb = {a for a in state.components.atoms(cid) if origin[a].kind == "chain" and origin[a].is_backbone}
    nbrs = {a: [b for b in state.atoms[a].edges if b in bb] for a in bb}
    head = min((a for a in bb if len(nbrs[a]) <= 1),
               key=lambda a: int(origin[a].segment.split()[1]))
    path, prev = [head], None
    while len(path) < len(bb):
        nxt = [b for b in nbrs[path[-1]] if b != prev]
        prev = path[-1]
        path.append(nxt[0])
    return path


def engine_cuts(spec, templates, seqs):
    """Per sequence: [(template, cut period q)] of every site, q the head-side
    backbone atom's position (1-based): fragments q and N - q periods."""
    state = ms.materialize_sequences(spec, seqs)
    engine = me.MeltEngine(spec, _Catalog(templates), T, numpy.random.default_rng(0), state=state)
    out = []
    for cid in sorted(state.components.ids()):
        position = {a: i + 1 for i, a in enumerate(backbone_path(state, cid))}
        entry = engine._scan([cid], state=state)[cid]
        sites = []
        for ci, (allowed, _leak) in K._uni_all(entry).items():
            for tr in allowed:
                (cut,) = [e for e in tr.edits if isinstance(e, ms.SetBond)]
                sites.append((engine.templates[ci].name, min(position[cut.a], position[cut.b])))
        out.append(sites)
    return out


def exact_moments(spec, templates, n, chains_mol, times, nbb=2):
    """Chain and end-radical pool moments (residues) of a monodisperse pool of
    n residues under the catalog's homolysis alone, exactly: each sequence
    decays at its own rate c(s); its cuts fire in proportion to their rates."""
    seqs = every_sequence(spec, n)
    weight = {t.name: t.weight for t in templates}
    cuts = [[(weight[name], qq) for name, qq in sites]
            for sites in engine_cuts(spec, templates, [s for s, _p in seqs])]
    N = n * nbb
    rows = []
    for t in times:
        chain, radical = [0.0] * 3, [0.0] * 3
        for (_seq, p), sites in zip(seqs, cuts):
            c = math.fsum(w for w, _q in sites)
            alive = math.exp(-c * t)
            for k in range(3):
                chain[k] += chains_mol * float(p) * alive * n ** k
            if c:
                fired = (1.0 - alive) / c
                for w, qq in sites:
                    for k in range(3):
                        radical[k] += chains_mol * float(p) * w * fired * ((qq / nbb) ** k + ((N - qq) / nbb) ** k)
        rows.append((chain, radical))
    return rows


# --------------------------------------------------------------------------
# the one-pool MoM run (the polymer-mom solver, the producer's terms)
# --------------------------------------------------------------------------
# The kMC line's own solver has no kernel_homolysis: the solve runs in a
# polymer-mom checkout (KMC_MOM_CHECKOUT), in a subprocess that imports rmgpy
# from it alone.  Only numbers cross: the seed's mass and moments and the
# producer's terms in, the moment rows out.
MOM_CHECKOUT_ENV = "KMC_MOM_CHECKOUT"
MOM_RUN = """
import json, sys
import numpy
import rmgpy
print("RMGPY " + rmgpy.__file__)
from rmgpy.molecule import Molecule
from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig
from rmgpy.species import Species
job = json.loads(sys.stdin.read())
def spc(smiles, label):
    return Species(label=label, molecule=[Molecule(smiles=smiles)])
core = [spc("N#N", "N2"),
        spc("CCCC", "P"), spc("CO", "P_mu0"), spc("C=O", "P_mu1"), spc("C#N", "P_mu2"),
        spc("[CH2]CCC", "Re"), spc("CCO", "Re_mu0"), spc("CC=O", "Re_mu1"), spc("CC#N", "Re_mu2")]
pools = [PolymerPoolConfig(label="P", xs=2, explicit_dp_to_species_index={}, mu_indices=(2, 3, 4),
                           monomer_mw_g_mol=job["mbar"], kernel_homolysis=job["terms"]),
         PolymerPoolConfig(label="Re", xs=2, explicit_dp_to_species_index={}, mu_indices=(6, 7, 8),
                           monomer_mw_g_mol=job["mbar"])]
mask = numpy.array([s.label == "N2" for s in core], dtype=bool)
rs = HybridPolymerSystem(T=job["T"], P=1.0e5, initial_mole_fractions={core[0]: 1.0}, V_poly=1.0,
                         polymer_pools=pools, mass_transfer=[], gas_species_mask=mask,
                         constant_gas_volume=False,
                         initial_polymer_moments={"P": tuple(job["mu"]), "Re": (0.0, 0.0, 0.0)},
                         termination=[])
rs.initialize_model(core, [], [], [])
rows = []
for t in job["times"]:
    if t > 0.0:
        rs.advance(t)
    y = numpy.asarray(rs.y)
    rows.append((list(map(float, y[2:5])), list(map(float, y[6:9]))))
print("ROWS " + json.dumps(rows))
"""


def mom_run(spec, terms, times):
    import os
    import subprocess
    import sys
    checkout = os.environ.get(MOM_CHECKOUT_ENV)
    if not checkout:
        pytest.skip(f"{MOM_CHECKOUT_ENV} names no polymer-mom checkout: the kMC line's solver has no "
                    "kernel_homolysis, so the one-pool MoM solve cannot run here")
    seed = copolymer_pool().copolymer_pool_seed(spec)
    job = {"mbar": float(seed.mean_residue_molar_mass), "mu": list(seed.moments_mol), "terms": terms,
           "T": T, "times": list(times)}
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    env["PYTHONPATH"] = checkout
    done = subprocess.run([sys.executable, "-c", MOM_RUN], input=json.dumps(job), capture_output=True,
                          text=True, cwd=checkout, env=env, timeout=600)
    assert done.returncode == 0, done.stderr[-4000:]
    (where,) = [x for x in done.stdout.splitlines() if x.startswith("RMGPY ")]
    assert os.path.realpath(where[len("RMGPY "):]).startswith(os.path.realpath(checkout) + os.sep), where
    (line,) = [x for x in done.stdout.splitlines() if x.startswith("ROWS ")]
    return seed, [tuple(r) for r in json.loads(line[len("ROWS "):])]


@pytest.fixture(scope="module")
def tables(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("tables")
    spec = ep_spec()
    out = {}
    for name, templates in (("uniform", cut_templates()), ("p_only", cut_templates(every=0.0, p_extra=W)),
                            ("mixed", cut_templates(p_extra=W))):
        doc, _path = catalog_doc(tmp, spec, templates, f"{name}.json")
        out[name] = (templates, copolymer_kernels().copolymer_kernel_rates(doc, spec, T))
    return spec, out


# --------------------------------------------------------------------------
# the adapter: one pool, mean residue mass
# --------------------------------------------------------------------------
def test_the_seed_is_one_pool_of_mean_residue_mass():
    spec = ep_spec({6: 10, 9: 5})
    seed = copolymer_pool().copolymer_pool_seed(spec)
    masses = {name: msq.residue_molar_mass(spec, name) for name, _r in spec.residues}
    assert seed.mean_residue_molar_mass == Fraction(2, 3) * masses["E"] + Fraction(1, 3) * masses["P"]
    assert seed.mean_residue_molar_mass == msq.mean_residue_molar_mass(spec)
    assert isinstance(seed.mean_residue_molar_mass, Fraction)
    assert seed.backbone_atoms_per_residue == 2
    mu0, mu1, mu2 = seed.moments_mol
    assert mu1 / mu0 == pytest.approx((6 * 10 + 9 * 5) / 15, rel=1e-14)
    assert mu2 / mu0 == pytest.approx((36 * 10 + 81 * 5) / 15, rel=1e-14)
    assert mu1 * float(seed.mean_residue_molar_mass) + mu0 * float(seed.end_molar_mass) == pytest.approx(
        spec.initial_mass_kg * 1000.0, rel=1e-14)


@pytest.mark.parametrize("architecture, name", [("network", "refused:network_topology"),
                                                ("branched", "refused:branched_copolymer")])
def test_networks_and_branched_copolymers_stay_refused_by_name(architecture, name):
    spec = dataclasses.replace(ep_spec(), architecture=architecture)
    with pytest.raises(ValueError, match=name):
        copolymer_pool().copolymer_pool_seed(spec)


def _polymer_calls(text):
    """{label: {keyword: literal}} of every polymer(...) call of a MoM input."""
    import ast
    out = {}
    for node in ast.walk(ast.parse(text)):
        if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "polymer":
            kw = {k.arg: ast.literal_eval(k.value) for k in node.keywords}
            out[kw["label"]] = kw
    return out


def test_the_mom_deck_carries_one_copolymer_pool_in_legacy_arguments(kmc_run, tmp_path):
    """polymer-mom's polymer() takes legacy arguments only: the copolymer
    pool is its proxy unit and caps, the seed's moments, and the mean
    residue's mass and exact element counts (unit_molar_mass, unit_atoms);
    the end-radical pool is born empty with the same unit.  No declaration
    crosses to the MoM checkout."""
    from momProducerTest import CHECKOUT
    from rmgpy.kmc import mom_producer as mp
    from rmgpy.molecule import Molecule
    deck, _out, spec, doc = kmc_run
    table = copolymer_kernels().copolymer_kernel_rates(doc, spec, T)
    path = tmp_path / "copolymer_kernels.json"
    path.write_text(json.dumps(table))
    (tmp_path / "db" / "input" / "kinetics" / "families").mkdir(parents=True)
    assert mp.produce(deck, checkout=CHECKOUT, temperature=T, end=1.0, snapshots=2, families=(), catalogs=(),
                      out=tmp_path / "mom", database=str(tmp_path / "db"), kernel_rates=str(path),
                      kernels=[HOMOLYSIS])["status"] == "written"
    text = (tmp_path / "mom" / "rmg" / "input.py").read_text()
    assert "polymer_spec" not in text
    seed = copolymer_pool().copolymer_pool_seed(spec)
    label = mp.mom_label(spec.label)
    calls = _polymer_calls(text)
    assert sorted(calls) == [label, label + "_re"]
    chain, radical = calls[label], calls[label + "_re"]
    for pool in (chain, radical):
        assert pool["unit_molar_mass"] == float(seed.mean_residue_molar_mass)
        assert {k: Fraction(v) for k, v in pool["unit_atoms"].items()} == {"C": Fraction(7, 3),
                                                                           "H": Fraction(14, 3)}
        assert Molecule().from_adjacency_list(pool["monomer"]).is_isomorphic(
            copolymer_pool().proxy_unit(spec))
        assert pool["end_groups"] == ["[CH3]", "[CH3]"]
    assert chain["moments"] == list(seed.moments_mol)
    assert chain["initial_mass"] == float(seed.mean_residue_molar_mass) * seed.moments_mol[1] / 1000.0
    assert chain["kernel_homolysis"]
    assert radical["moments"] == [0.0, 0.0, 0.0] and radical["initial_mass"] == 0.0
    assert radical["end_radical_site"] == "primary"


# --------------------------------------------------------------------------
# the kernels: per-residue sites weighted by the word frequencies
# --------------------------------------------------------------------------
@pytest.mark.parametrize("statistics", [None, STICKY])
@pytest.mark.parametrize("which", ["uniform", "p_only"])
def test_the_weighting_rule_is_the_expectation_over_every_sequence(tmp_path, which, statistics):
    """Each template's sites a*N + b (N periods) are its exact mean count over
    every sequence of the chain, each at its probability."""
    spec = ep_spec(statistics=statistics)
    templates = cut_templates() if which == "uniform" else cut_templates(every=0.0, p_extra=W)
    doc, _path = catalog_doc(tmp_path, spec, templates)
    table = copolymer_kernels().copolymer_kernel_rates(doc, spec, T)
    rows = {r["template"]: r for r in table["templates"]}
    templates = list(me.PayloadCatalog(doc["payloads"]).templates(T))     # the names the catalog loads
    assert set(rows) == {t.name for t in templates}
    for n in (4, 7):
        seqs = every_sequence(spec, n)
        cuts = engine_cuts(spec, templates, [s for s, _p in seqs])
        for t in templates:
            mean = sum((p * sum(1 for name, _q in sites if name == t.name)
                        for (_s, p), sites in zip(seqs, cuts)), Fraction(0))
            a, b = Fraction(rows[t.name]["sites"]["a"]), Fraction(rows[t.name]["sites"]["b"])
            assert a * (2 * n) + b == mean, (t.name, n, a, b)


def test_a_uniform_catalog_is_not_residue_selective(tables):
    _spec, out = tables
    kernel = out["uniform"][1]["kernels"][HOMOLYSIS]
    assert kernel["residue_selective"] is False
    assert kernel["mom"]["min_fragment"] == 2
    assert not [lim for lim in out["uniform"][1]["limitations"] if lim["name"] == "limited:composition_drift"]


def test_a_non_selective_catalog_is_exact_as_one_pool(tables):
    spec, out = tables
    templates, table = out["uniform"]
    terms = copolymer_kernels().mom_homolysis_terms(table, phi=1.0, daughter="Re")
    times = [0.0, 0.5, 2.0, 6.0]
    seed, mom = mom_run(spec, terms, times)
    (n, _count), = seed.dp_counts
    exact = exact_moments(spec, templates, n, seed.moments_mol[0], times)
    for (chain, radical), (chain_x, radical_x) in zip(mom, exact):
        assert chain == pytest.approx(chain_x, rel=1e-6)
        assert radical == pytest.approx(radical_x, rel=1e-6, abs=1e-12 * seed.moments_mol[1])


def test_a_residue_selective_catalog_carries_the_composition_drift_limitation(tables):
    """Every bond at W plus P's bonds at W again: E residues react at 2W, P
    residues at 4W; the selective part is the propensity above the slowest
    residue window, (a - a_floor) mu1 over a mu1 + b mu0, a = W(1 + x_P),
    a_floor = W per period, b = -3 W (1 + x_P)."""
    _spec, out = tables
    _templates, table = out["mixed"]
    kernel = table["kernels"][HOMOLYSIS]
    assert kernel["residue_selective"] is True
    assert kernel["per_residue_rate"] == pytest.approx({"E": 2 * W, "P": 4 * W}, rel=1e-14)
    (limitation,) = [lim for lim in table["limitations"] if lim["name"] == "limited:composition_drift"]
    assert limitation["kernels"] == [HOMOLYSIS]
    assert set(limitation["templates"]) == {t["template"] for t in kernel["templates"]}
    drift = copolymer_kernels().composition_drift(table, [(0.0, (1.0, 8.0)), (1.0, (2.0, 8.0))])
    assert drift["name"] == "limited:composition_drift"
    # mu in residues: mu1 periods = 2 mu1; fraction = (W/3) 2 mu1 / ((4W/3) 2 mu1 - 4 W mu0)
    assert [row["selective_fraction"] for row in drift["series"]] == pytest.approx(
        [(16 / 3) / (64 / 3 - 4), (16 / 3) / (64 / 3 - 8)], rel=1e-12)


def test_the_selective_gap_is_reported_not_corrected(tables):
    """P-only scission: one pool decays every chain at the mean rate; the
    exact sequence-level decay is slower (a chain with few P lives longer).
    MoM integrates its own law; the difference is the method's."""
    spec, out = tables
    templates, table = out["p_only"]
    terms = copolymer_kernels().mom_homolysis_terms(table, phi=1.0, daughter="Re")
    times = [0.0, 20.0]
    seed, mom = mom_run(spec, terms, times)
    (n, _count), = seed.dp_counts
    exact = exact_moments(spec, templates, n, seed.moments_mol[0], times)
    assert mom[0][0] == pytest.approx(exact[0][0], rel=1e-12)
    assert exact[1][0][0] > 1.05 * mom[1][0][0]
    assert table["kernels"][HOMOLYSIS]["residue_selective"] is True


def test_a_feature_kernel_on_a_copolymer_is_refused_by_name(tmp_path):
    """A radical end on a copolymer sits on one residue or another; one
    end-radical pool cannot carry which."""
    spec = ep_spec()
    templates = cut_templates() + [t for t in S.pe_templates() if t.name == "recombination"]
    doc, _path = catalog_doc(tmp_path, spec, templates)
    table = copolymer_kernels().copolymer_kernel_rates(doc, spec, T)
    (reason,) = [r["reason"] for r in table["refused"] if r["template"].startswith("recombinatio")]
    assert reason.startswith("refused:residue_context")
    assert HOMOLYSIS in table["kernels"]


# --------------------------------------------------------------------------
# the comparison: composition is a kMC observable
# --------------------------------------------------------------------------
TIMES = [0.0, 1.0, 2.0]


@pytest.fixture(scope="module")
def kmc_run(tmp_path_factory):
    """A seeded kMC run of the copolymer under the mixed (selective) catalog."""
    tmp = tmp_path_factory.mktemp("kmc")
    mapping = ep_mapping({8: 6})
    spec = PolymerSpec.from_mapping(mapping)
    deck = write_v2_deck(tmp / "ep.py", mapping)
    doc, catalog = catalog_doc(tmp, spec, cut_templates(p_extra=W))
    kr.run(deck, [catalog], T, TIMES, [1, 2], tmp / "out", max_events=60)
    return deck, tmp / "out", spec, doc


def test_the_kmc_rows_of_a_copolymer_carry_its_composition(kmc_run):
    """Refused today (refused:schema_v2_mom; the rows read one repeat unit)."""
    from rmgpy.kmc import comparison_rows as cr
    _deck, out, _spec, _doc = kmc_run
    doc = cr.kmc_rows(out)
    rows = {(r["time_s"], r["row"]): r for r in doc["rows"]}
    for t in TIMES:
        values = [rows[(t, f"composition_fraction.{r}")] for r in ("E", "P")]
        assert all(v["kind"] == "exact" for v in values)
        assert sum(v["value"] for v in values) == pytest.approx(1.0, rel=1e-12)
        assert rows[(t, "whole_mu0_mol")]["kind"] == "exact"
    assert doc["geometry"] == {"periods_per_declared_unit": 2, "backbone_atoms_per_period": 1,
                               "cutoff_periods": 8}


def copolymer_mom_file(kmc_run, tmp_path):
    """The producer's kernel-mode copolymer deck, converted from a fake pool
    table at Mbar0: (mom.json path, its deck, the seed, the pool label, the
    radical pool's chain moments at t)."""
    from momProducerTest import CHECKOUT, ZERO_GAS, write_run
    from rmgpy.kmc import mom_producer as mp
    deck, out, spec, doc = kmc_run
    table = copolymer_kernels().copolymer_kernel_rates(doc, spec, T)
    path = tmp_path / "copolymer_kernels.json"
    path.write_text(json.dumps(table))
    (tmp_path / "db" / "input" / "kinetics" / "families").mkdir(parents=True)
    args = dict(checkout=CHECKOUT, temperature=T, end=TIMES[-1], snapshots=len(TIMES), families=(),
                catalogs=(), out=tmp_path / "mom", database=str(tmp_path / "db"), kernel_rates=str(path),
                kernels=[HOMOLYSIS])
    assert mp.produce(deck, **args)["status"] == "written"
    label = mp.mom_label(spec.label)
    seed = copolymer_pool().copolymer_pool_seed(spec)
    mu = seed.moments_mol
    d = 1e-6
    # the pool table as a copolymer run writes it: the mean residue's mass (PolymerPool.to_config)
    write_run(tmp_path / "mom", times=[0.0, 1.0, 2.0, 2.5],
              pools=[(label, lambda t: (mu[0] - d * t, mu[1] - 8 * d * t, mu[2] - 64 * d * t), 0.0),
                     (label + "_re", lambda t: (2 * d * t, 8 * d * t, 40 * d * t), 0.0)],
              gas=ZERO_GAS, rows=[], end_groups=("[CH3]", "[CH3]"), label=label,
              mass_override=float(seed.mean_residue_molar_mass))
    assert mp.produce(deck, mode="convert", **args)["status"] == "converted"
    return tmp_path / "mom" / "mom.json", deck, seed, label


def test_the_mom_rows_of_a_copolymer_run_withhold_composition_by_its_limitation(kmc_run, tmp_path):
    """The producer's kernel mode writes the copolymer as one pool from the
    copolymer kernel table; its file carries the limitation, and its rows
    withhold the composition the kMC rows report, naming it."""
    from rmgpy.kmc import comparison_rows as cr
    path, deck, _seed, _label = copolymer_mom_file(kmc_run, tmp_path)
    out = kmc_run[1]
    mom = cr.mom_rows(path, deck)
    assert "limited:composition_drift" in mom["producer"]["flags"]
    rows = {(r["time_s"], r["row"]): r for r in mom["rows"]}
    for t in TIMES:
        assert rows[(t, "composition_fraction.P")]["kind"] == "withheld"
        assert "limited:composition_drift" in rows[(t, "composition_fraction.P")]["reason"]
        assert rows[(t, "whole_mu0_mol")]["kind"] == "exact"
    kmc = cr.kmc_rows(out)
    mom["binding"]["pass"] = kmc["binding"]["pass"]           # the comparison under test is kMC vs MoM rows
    joined = {(r["time_s"], r["row"]): r for r in cr.join(kmc, mom)}
    assert joined[(1.0, "composition_fraction.P")]["status"] == "withheld_b"
    assert joined[(0.0, "whole_mu0_mol")]["status"] == "compared"


def test_the_classification_free_rows_of_a_copolymer_mom_file(kmc_run, tmp_path):
    """I-171's rows on a copolymer MoM file.  The pool counts residues of the
    mean residue mass Mbar0, so a chain of n residues weighs anywhere from
    n M_E + e to n M_P + e: the booked inventory Mbar0 mu1 + e_p mu0 is exact
    (MoM's own booking, the seed's mass here, compared with kMC's total),
    mass_outside_certified_g is the inventory with nothing certified, and
    certification and the low-DP bounds of an occupied pool are withheld by
    COPOLYMER_REACH; an empty pool still bounds nothing at [0, 0]."""
    from rmgpy.kmc import comparison_rows as cr
    from rmgpy.kmc import mom_producer as mp
    path, deck, seed, label = copolymer_mom_file(kmc_run, tmp_path)
    mom = cr.mom_rows(path, deck)
    rows = {(r["time_s"], r["row"]): r for r in mom["rows"]}
    booked = (float(seed.mean_residue_molar_mass) * seed.moments_mol[1]
              + float(seed.end_molar_mass) * seed.moments_mol[0])
    species = sorted({name for (_t, name) in rows if name.startswith("certified_species_mol.")})
    assert species                                        # the deck's declared species reach the rows
    for t in TIMES:
        inventory = rows[(t, "inventory_mass_g")]
        assert inventory["kind"] == "exact"
        assert inventory["value"] == pytest.approx(booked, rel=1e-12)
        outside = rows[(t, "mass_outside_certified_g")]
        assert outside["kind"] == "exact" and outside["certified"] == []
        assert outside["value"] == pytest.approx(inventory["value"], rel=1e-15)
        for name in species:
            assert rows[(t, name)]["kind"] == "withheld"
            assert mp.COPOLYMER_REACH in rows[(t, name)]["reason"]
        main = rows[(t, f"low_dp_mass_bounds_g.{label}")]
        assert main["kind"] == "withheld" and mp.COPOLYMER_REACH in main["reason"]
        assert rows[(t, "low_dp_mass_bounds_g")]["kind"] == "withheld"
    empty = rows[(0.0, f"low_dp_mass_bounds_g.{label}_re")]   # the end-radical pool is born empty
    assert (empty["kind"], empty["lower"], empty["upper"]) == ("interval", 0.0, 0.0)
    kmc = cr.kmc_rows(kmc_run[1])
    mom["binding"]["pass"] = kmc["binding"]["pass"]
    joined = {(r["time_s"], r["row"]): r for r in cr.join(kmc, mom)}
    for t in TIMES:
        assert joined[(t, "inventory_mass_g")]["status"] == "compared"
        assert joined[(t, "low_dp_mass_bounds_g")]["status"] == "withheld_b"
