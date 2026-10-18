"""Red tests: the method-of-moments producer for the runner's ``--mom`` file
(``rmgpy.kmc.mom_producer``, I-169).  No RMG run: the MoM side is a
hand-built solver profile, pool table, species dictionary and annotated
mechanism with known answers.

Expected values are derived here without the producer: the moments of the
declared ensemble by hand, the ledger-booked condensed mass written out, the
piecewise-linear interpolant of samples of t^2 in closed form, and the
profile column of a long moment-dummy label by the MoM writer's rule.  The
producer must export a key only under kMC's definition of it: at t = 0 every
exported quantity equals the kMC observation of the same deck.
"""

import csv
import hashlib
import json

import numpy as np
import pytest

import melt_engine_support as S
from meltRunTest import pe_catalog_file, write_deck, write_catalog, toy_payloads

from rmgpy import constants
from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import melt_observables as mo
from rmgpy.kmc import mom_producer as mp
from rmgpy.kmc import run as kr
from rmgpy.molecule import Molecule
from rmgpy.polymer_spec import LinearHomopolymerSpec

CHECKOUT = "/nonexistent/mom-checkout"   # never run here; only recorded
PE_FAMILIES = ["toy_beta_scission", "toy_branching", "toy_disproportionation", "toy_h_abstraction",
               "toy_homolysis", "toy_recombination"]
PE_LIVE = ("homolysis", "recombination", "h_abstraction", "disproportionation", "beta_scission",
           "mid_end_recombination", "intra_mid_end")
LABEL = "toylinearp"                      # mom_label of "toy linear polymer (CH2-CH2)"
DAUGHTER = LABEL + "_rad_primary_end"     # a spawned pool with a long label
BASE = {10: 1, 20: 1}

GOLDEN = '''polymer(
    label='toylinearp',
    monomer="""multiplicity 3
1 *2 C u1 p0 c0 {2,S} {3,S} {4,S}
2 *1 C u1 p0 c0 {1,S} {5,S} {6,S}
3    H u0 p0 c0 {1,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {2,S}
6    H u0 p0 c0 {2,S}
""",
    end_groups=['[H]', '[H]'],
    cutoff=4,
    moments=[2.3651094350938804, 35.47664152640821, 591.2773587734702],
    initial_mass=0.9952320765335197,
)
phase = polymer_phase(label='toylinearp_melt', species=['toylinearp'], solvent='toylinearp', density=(850.0, 'kg/m^3'))
hybridPolymerReactor(
    temperature=(700.0, 'K'),
    pressure=(1.0, 'bar'),
    initialMoles={'toylinearp': 2.3651094350938804, 'N2': 1.0},
    polymerPhase=phase,
    terminationTime=(40.0, 's'),
    constant_gas_volume=False,
)
'''

ADJ = {
    "N2": "1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}\n",
    "H": "multiplicity 2\n1 H u1 p0 c0\n",
    "H2": "1 H u0 p0 c0 {2,S}\n2 H u0 p0 c0 {1,S}\n",
    "dummy": "1 Ne u0 p4 c0\n",
}


def spec_of(counts, **kw):
    return LinearHomopolymerSpec.from_mapping(S.pe_mapping(counts, **kw))


def masses(counts=BASE):
    spec = spec_of(counts)
    return spec.repeat_unit.composition.molar_mass_g_mol(), spec.end_composition.molar_mass_g_mol()


def mw(smiles):
    return Molecule(smiles=smiles).get_molecular_weight() * 1000.0


def make_setup(tmp_path, counts=BASE, **kw):
    deck = write_deck(tmp_path / "pe_deck.py", S.pe_mapping(counts, **kw))
    catalog = pe_catalog_file(tmp_path)
    (tmp_path / "db" / "input" / "kinetics" / "families").mkdir(parents=True, exist_ok=True)
    return {"deck": deck, "catalog": catalog, "db": str(tmp_path / "db"), "tmp": tmp_path,
            "spec": spec_of(counts, **kw)}


@pytest.fixture
def setup(tmp_path):
    return make_setup(tmp_path)


def produce(setup, out, **kw):
    args = dict(checkout=CHECKOUT, temperature=700.0, end=40.0, snapshots=9, families=PE_FAMILIES,
                catalogs=[setup["catalog"]], out=out, database=setup["db"])
    args.update(kw)
    return mp.produce(setup["deck"], **args)


# --------------------------------------------------------------------------
# a hand-built MoM run
# --------------------------------------------------------------------------
def write_run(out, *, times, pools, gas=(), admitted=None, label=LABEL, end_groups=("[H]", "[H]"),
              mass_override=None, duplicate=None, extra_columns=()):
    """A MoM output directory under ``out/rmg``.

    ``pools``: [(label, f(t) -> (mu0, mu1, mu2), defect_g_mol)]; ``gas``:
    [(identifier, adjacency, f(t) -> mol)]; ``admitted``: {family: (rows,
    admitted rows)} of pool-coupled reactions (default one admitted
    toy_homolysis row)."""
    rmg = out / "rmg"
    for sub in ("solver", "chemkin"):
        (rmg / sub).mkdir(parents=True, exist_ok=True)
    a, _e = masses()
    moment_ids = {p: [mp.moment_column(f"{p}_mu{k}") for k in range(3)] for p, _f, _d in pools}
    proxy = f"{label}(3)"
    header = ["Time (s)", "Volume (m^3)", "N2"] + [c for p, _f, _d in pools for c in moment_ids[p]] + \
        [g[0] for g in gas] + [proxy] + list(extra_columns)
    rows = []
    for t in times:
        row = [t, 1.0, 1.0]
        for _p, f, _d in pools:
            row += list(f(t))
        row += [g[2](t) for g in gas] + [7.0] + [0.0] * len(extra_columns)
        rows.append(row)
    if duplicate is not None:
        k, delta = duplicate
        bad = list(rows[k])
        bad[3] += delta
        rows.insert(k + 1, bad)
    with open(rmg / "solver" / "simulation_1_8.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)
    with open(rmg / "solver" / "simulation_1_7.csv", "w", newline="") as handle:   # an earlier core: ignored
        csv.writer(handle).writerow(["Time (s)", "Volume (m^3)"])
    dictionary = [f"N2\n{ADJ['N2']}"]
    for p, _f, _d in pools:
        dictionary += [f"{c}\n{ADJ['dummy']}" for c in moment_ids[p]]
    dictionary += [f"{g[0]}\n{g[1]}" for g in gas]
    dictionary.append(f"{proxy}\n" + Molecule(smiles="CCCCCC").to_adjacency_list())
    (rmg / "chemkin" / "species_dictionary.txt").write_text("\n".join(dictionary) + "\n")
    pool_rows = [{"label": p, "monomer_mw_g_mol": a if mass_override is None else mass_override,
                  "chain_mass_defect_g_mol": d, "end_groups": list(end_groups),
                  "bookkeeping_species": moment_ids[p] + ([proxy] if p == label else [])}
                 for p, _f, d in pools]
    admitted = {"toy_homolysis": (1, 1)} if admitted is None else admitted
    reactions, annotated, k = [], ["REACTIONS    KCAL/MOLE   MOLES", ""], 0
    for family, (n, ok) in sorted(admitted.items()):
        for j in range(n):
            row = {"id": f"r{k}", "cantera": {"index": k, "equation": f"X{k} <=> {proxy}"}, "refused": j >= ok}
            if j >= ok:
                row["refused_reason"] = "conduit-deferred"
            reactions.append(row)
            annotated += [f"! Reaction index: Chemkin #{k + 1}; RMG #{k + 1}", f"! Template reaction: {family}",
                          f"X{k}<=>{proxy}    1.0 0.0 0.0", ""]
            k += 1
    (rmg / "chemkin" / "polymer_pools.json").write_text(json.dumps(
        {"schema_version": "2.0", "pools": pool_rows, "reactions": reactions,
         "conventions": {"stale_topology": False}}))
    (rmg / "chemkin" / "chem_annotated.inp").write_text("\n".join(annotated + ["END", ""]))


def declared(spec):
    return [float(x) for x in spec.macroscopic_moments_mol]


ZERO_GAS = [("H(1)", ADJ["H"], lambda t: 0.0), ("H2(2)", ADJ["H2"], lambda t: 0.0)]


def constant_run(out, spec, times=(0.0, 10.0, 20.0, 40.0, 45.0), **kw):
    mu = declared(spec)
    write_run(out, times=list(times), pools=[(LABEL, lambda t: mu, 0.0)], gas=ZERO_GAS, **kw)


PROFILE_T = [0.0, 10.0, 20.0, 40.0, 45.0]       # the solver's last step passes the end time


def interpolant(f, t):
    """The piecewise-linear interpolant of f sampled at PROFILE_T."""
    for lo, hi in zip(PROFILE_T, PROFILE_T[1:]):
        if lo <= t <= hi:
            return f(lo) + (f(hi) - f(lo)) * (t - lo) / (hi - lo)
    raise AssertionError(t)


def moving_run(out, spec):
    """Main pool loses mass to H2 and chains to a spawned daughter pool;
    total mass closes exactly (linear in the columns, so also between rows)."""
    mu = declared(spec)
    a, _e = masses()
    m_h2 = mw("[H][H]")
    d, h = 1e-3, 2e-3
    r = 5 * d + m_h2 * h / a

    def main(t):
        return mu[0] * (1 + 1e-3 * t * t), mu[1] - r * t * t, mu[2] * (1 + 1e-3 * t * t)

    def daughter(t):
        return d * t * t, 5 * d * t * t, 40 * d * t * t

    write_run(out, times=PROFILE_T, pools=[(LABEL, main, 0.0), (DAUGHTER, daughter, 0.0)],
              gas=[("H(1)", ADJ["H"], lambda t: 0.0), ("H2(2)", ADJ["H2"], lambda t: h * t * t)])
    return main, daughter, h, m_h2


# --------------------------------------------------------------------------
# deck -> input file
# --------------------------------------------------------------------------
def test_deck_gives_the_golden_mom_input(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    text = (out / "rmg" / "input.py").read_text()
    start = text.index("polymer(\n")
    assert text[start:text.index("model(")] == GOLDEN
    a, e = masses()
    pair = (10 * a + e) + (20 * a + e)
    mu0 = 2 * 1000.0 / pair
    assert [2.3651094350938804, 35.47664152640821, 591.2773587734702] == pytest.approx(
        [mu0, 15 * mu0, 250 * mu0], rel=1e-14)
    assert 0.9952320765335197 == pytest.approx(a * 15 * mu0 / 1000.0, rel=1e-14)
    assert f"kineticsFamilies={PE_FAMILIES!r}," in text
    assert "reactionLibraries=[]," in text and "seedMechanisms=[]," in text
    for kernel in ("k_scission", "k_unzip", "k_homolysis", "radical_qssa_unzip", "side_group_homolysis",
                   "explicit_dp=True"):
        assert kernel not in text
    for name, smiles in (("H", "[H]"), ("H2", "[H][H]")):
        assert f"species(label={name!r}, reactive=True, structure=SMILES({smiles!r}))" in text
    assert "species(label='N2', reactive=False, structure=SMILES('N#N'))" in text
    assert "saveSimulationProfiles=True" in text
    assert (out / "rmg" / "rmgrc").read_text() == f"database.directory : {setup['db']}/input\n"


def test_conventions_are_written_and_name_every_setting(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    doc = json.loads((out / "conventions.json").read_text())
    names = {c["name"]: c for c in doc["conventions"]}
    for name in ("label", "legacy_translation", "initial_moments", "mass_booking", "mass_closure",
                 "tail_keys", "volume", "declared_species", "bath_gas", "pressure", "families",
                 "family_check", "database", "kernels_off", "no_cage", "model_settings", "profile",
                 "time_interpolation", "species_constraints", "phase_kinetics", "inactive_chemistry",
                 "initial_inventory"):
        assert name in names and names[name]["statement"], name
    assert names["time_interpolation"]["value"] == "linear in time"
    assert names["label"]["value"] == LABEL
    assert set(doc["quantity_sources"]) | set(doc["not_produced"]) == mp.runner_quantity_keys()


# --------------------------------------------------------------------------
# item 1: a key is exported only under kMC's definition of it
# --------------------------------------------------------------------------
def test_tail_keys_are_not_produced():
    for key in ("mu0_mol", "mu1_mol", "mu2_mol", "Mn_g_mol", "Mw_g_mol"):
        assert key in mp.NOT_PRODUCED and key not in mp.QUANTITY_SOURCES
    assert set(mp.QUANTITY_SOURCES) == {"polymer_condensed_mass_g", "polymer_conversion",
                                        "polymer_derived_molecules_mol", "species_mol"}
    assert set(mp.QUANTITY_SOURCES).isdisjoint(mp.NOT_PRODUCED)
    assert set(mp.QUANTITY_SOURCES) | set(mp.NOT_PRODUCED) == mp.runner_quantity_keys()


@pytest.mark.parametrize("counts", [{2: 1, 8: 1}, BASE])
def test_every_exported_quantity_equals_the_kmc_value_of_the_same_deck_at_t0(tmp_path, counts):
    """Spar's fixture (DP2 + DP8 with xs = 4: the kMC tail holds one chain, the
    MoM pool both): no exported row can disagree at t = 0 by definition."""
    s = make_setup(tmp_path, counts)
    out = tmp_path / "out"
    produce(s, out)
    constant_run(out, s["spec"])
    produce(s, out, mode="convert")
    quantities = json.loads((out / "mom.json").read_text())[0]["quantities"]
    engine = me.MeltEngine(s["spec"], S.ToyCatalog(names=PE_LIVE), 700.0, np.random.default_rng(1))
    obs = mo.MeltObserver(s["spec"], xs=4, identity=engine.identity).observe(engine)
    rows = [r for r in mo.side_by_side(obs, quantities) if r["mom"] is not None]
    assert {r["quantity"] for r in rows} >= {"polymer_condensed_mass_g", "polymer_derived_molecules_mol",
                                             "polymer_conversion", "species_mol.H", "species_mol.H2"}
    assert all(r["agreement"] == "equal" for r in rows), [r for r in rows if r["agreement"] != "equal"]


# --------------------------------------------------------------------------
# item 2: the condensed mass is MoM's ledger booking
# --------------------------------------------------------------------------
def test_scission_does_not_invent_end_group_mass(tmp_path):
    """PE DP6 -> two DP3 (chain count doubles, repeat units kept): the mass
    stays 1000 g (adding e per chain would give 1011.8 g)."""
    s = make_setup(tmp_path, {6: 1})
    out = tmp_path / "out"
    produce(s, out, end=10.0, snapshots=2)
    mu = declared(s["spec"])

    def split(t):
        return mu[0] * (1 + t / 10.0), mu[1], mu[2] * (1 - 0.5 * t / 10.0)

    write_run(out, times=[0.0, 10.0], pools=[(LABEL, split, 0.0)], gas=ZERO_GAS)
    produce(s, out, end=10.0, snapshots=2, mode="convert")
    entries = json.loads((out / "mom.json").read_text())
    a, e = masses()
    assert 1000.0 * (1 + e / (6 * a + e)) == pytest.approx(1011.8, abs=0.1)     # the invented value
    for entry in entries:
        assert entry["quantities"]["polymer_condensed_mass_g"] == pytest.approx(1000.0, rel=1e-12)
        assert entry["quantities"]["polymer_conversion"] == pytest.approx(0.0, abs=1e-12)
    assert entries[1]["quantities"]["polymer_derived_molecules_mol"] == pytest.approx(2 * mu[0], rel=1e-12)


def test_profile_converts_to_the_mom_file_at_the_observation_times(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    main, daughter, h, m_h2 = moving_run(out, setup["spec"])
    produce(setup, out, mode="convert")
    entries = json.loads((out / "mom.json").read_text())
    times = kr.observation_times(end=40.0, snapshots=9)
    assert [x["time_s"] for x in entries] == times
    for entry in entries:
        t, q = entry["time_s"], entry["quantities"]
        mu0 = interpolant(lambda s: main(s)[0], t) + interpolant(lambda s: daughter(s)[0], t)
        mass = 1000.0 - m_h2 * interpolant(lambda s: h * s * s, t)
        if t == 25.0:   # closed form between the 20 s and 40 s samples of t^2: (20 + 40) 25 - 800 = 700
            assert interpolant(lambda s: s * s, t) == pytest.approx(700.0, rel=1e-15)
        assert q["polymer_derived_molecules_mol"] == pytest.approx(mu0, rel=1e-12)
        assert q["polymer_condensed_mass_g"] == pytest.approx(mass, rel=1e-12)
        assert q["polymer_conversion"] == pytest.approx(1.0 - mass / 1000.0, rel=1e-9, abs=1e-15)
        assert q["species_mol"]["H"] == 0.0
        assert q["species_mol"]["H2"] == pytest.approx(interpolant(lambda s: h * s * s, t), rel=1e-12)
        assert not set(q) & set(mp.NOT_PRODUCED)
    manifest = json.loads((out / "producer.json").read_text())
    assert manifest["read"]["profile"].endswith("simulation_1_8.csv")
    assert manifest["read"]["pool_labels"] == [LABEL, DAUGHTER]


def test_chain_mass_defect_enters_the_booked_mass(setup):
    """A feature pool whose chains lost H2 (defect = M(H2) per chain): chains
    move to it while H2 is released; the booked mass is 1000 - M(H2) y t^2."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    mu = declared(setup["spec"])
    m_h2 = mw("[H][H]")
    y = 1e-3

    def main(t):
        return mu[0] - y * t * t, mu[1] - 15 * y * t * t, mu[2] - 250 * y * t * t

    def feature(t):
        return y * t * t, 15 * y * t * t, 250 * y * t * t

    write_run(out, times=PROFILE_T, pools=[(LABEL, main, 0.0), (LABEL + "_mod", feature, m_h2)],
              gas=[("H(1)", ADJ["H"], lambda t: 0.0), ("H2(2)", ADJ["H2"], lambda t: y * t * t)])
    produce(setup, out, mode="convert")
    for entry in json.loads((out / "mom.json").read_text()):
        t = entry["time_s"]
        assert entry["quantities"]["polymer_condensed_mass_g"] == pytest.approx(
            1000.0 - m_h2 * interpolant(lambda s: y * s * s, t), rel=1e-12)


def test_a_booking_that_does_not_close_exports_no_mass(setup):
    """A capped DP3 chain released to the gas phase: MoM's ledger books its end
    groups as condensed while the gas species carries them too, so the booked
    total grows by e per released chain: mass and conversion are not exported
    where the booking does not close."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    mu = declared(setup["spec"])
    x = 1e-2

    def main(t):
        return mu[0] - x * t * t, mu[1] - 3 * x * t * t, mu[2] - 9 * x * t * t

    hexane = Molecule(smiles="CCCCCC").to_adjacency_list()
    write_run(out, times=PROFILE_T, pools=[(LABEL, main, 0.0)],
              gas=ZERO_GAS + [("C6H14(4)", hexane, lambda t: x * t * t)])
    produce(setup, out, mode="convert")
    entries = json.loads((out / "mom.json").read_text())
    assert "polymer_condensed_mass_g" in entries[0]["quantities"]
    for entry in entries[1:]:
        assert "polymer_condensed_mass_g" not in entry["quantities"]
        assert "polymer_conversion" not in entry["quantities"]
        assert "polymer_derived_molecules_mol" in entry["quantities"]
    diagnostics = json.loads((out / "mom_diagnostics.json").read_text())
    _a, e = masses()
    assert diagnostics[-1]["mass_closure_residual_g"] == pytest.approx(e * x * 40.0 ** 2, rel=1e-6)


# --------------------------------------------------------------------------
# item 3: inactive chemistry
# --------------------------------------------------------------------------
def test_inactive_chemistry_is_refused(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], admitted={"toy_homolysis": (17, 0), "R_Recombination": (2, 2)})
    with pytest.raises(mp.ProducerError, match="inactive_chemistry"):
        produce(setup, out, mode="convert")
    assert not (out / "mom.json").exists()
    manifest = json.loads((out / "producer.json").read_text())
    assert manifest["status"] == "refused:inactive_chemistry"
    census = manifest["census"]
    assert census["families"]["toy_homolysis"] == {"rows": 17, "admitted": 0}
    assert census["families"]["R_Recombination"] == {"rows": 2, "admitted": 2}   # not a requested family
    assert census["admitted_in_requested_families"] == 0


def test_no_pool_reaction_at_all_is_inactive(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], admitted={})
    with pytest.raises(mp.ProducerError, match="inactive_chemistry"):
        produce(setup, out, mode="convert")


def test_inactive_chemistry_allowed_is_flagged(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], admitted={"toy_homolysis": (3, 0)})
    produce(setup, out, mode="convert", allow_inactive=True)
    assert (out / "mom.json").exists()
    manifest = json.loads((out / "producer.json").read_text())
    conv = json.loads((out / "conventions.json").read_text())
    for doc in (manifest, conv):
        assert doc["flags"]["inactive_chemistry"]["admitted_in_requested_families"] == 0
    assert manifest["status"] == "converted"


def test_admitted_chemistry_converts_unflagged(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], admitted={"toy_homolysis": (3, 1)})
    manifest = produce(setup, out, mode="convert")
    assert manifest["census"]["admitted_in_requested_families"] == 1
    assert "inactive_chemistry" not in manifest["flags"]


# --------------------------------------------------------------------------
# item 4: phase kinetics
# --------------------------------------------------------------------------
def test_phase_kinetics_ratio_is_recorded_and_flagged(setup):
    out = setup["tmp"] / "out"
    manifest = produce(setup, out)
    conv = json.loads((out / "conventions.json").read_text())
    names = {c["name"]: c["value"] for c in conv["conventions"]}
    v_kmc = 1.0 / 850.0                                 # 1 kg dry polymer at 850 kg/m^3
    v_gas = 1.0 * constants.R * 700.0 / 1.0e5          # 1 mol bath at 700 K, 1 bar
    ratio = v_kmc / v_gas
    kinetics = names["phase_kinetics"]
    assert set(kinetics["ratio_by_species"]) == {"H", "H2"}
    for value in kinetics["ratio_by_species"].values():
        assert value == pytest.approx(ratio, rel=1e-12)
    for doc in (manifest, conv):
        assert doc["flags"]["kinetics_not_comparable"]["ratio_by_species"]["H"] == pytest.approx(ratio, rel=1e-12)


# --------------------------------------------------------------------------
# item 5: species naming
# --------------------------------------------------------------------------
def test_a_formula_named_species_column_resolves(tmp_path):
    """RMG writes monomer_product as C2H4(4) (label + index > 16 characters):
    it is found by its structure in the species dictionary."""
    s = make_setup(tmp_path, BASE, monomer=S.ETHYLENE_GRAPH)
    out = tmp_path / "out"
    produce(s, out)
    mu = declared(s["spec"])
    ethylene = Molecule(smiles="C=C").to_adjacency_list()
    write_run(out, times=PROFILE_T, pools=[(LABEL, lambda t: mu, 0.0)],
              gas=ZERO_GAS + [("C2H4(4)", ethylene, lambda t: 1e-3 * t)])
    produce(s, out, mode="convert")
    q = json.loads((out / "mom.json").read_text())[-1]["quantities"]
    assert q["species_mol"]["monomer_product"] == pytest.approx(0.04, rel=1e-12)
    assert set(q["species_mol"]) == {"monomer_product", "H", "H2"}


def test_a_species_label_colliding_with_the_pool_label_is_refused(tmp_path):
    s = make_setup(tmp_path, BASE, small={LABEL: S.METHYL})
    with pytest.raises(mp.ProducerError, match="collide"):
        produce(s, tmp_path / "out")
    assert not (tmp_path / "out" / "rmg" / "input.py").exists()


# --------------------------------------------------------------------------
# items 6, 7 and the profile's own checks
# --------------------------------------------------------------------------
def test_a_wrong_initial_inventory_is_refused(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    mu = declared(setup["spec"])
    write_run(out, times=PROFILE_T, pools=[(LABEL, lambda t: [0.5 * x for x in mu], 0.0)], gas=ZERO_GAS)
    with pytest.raises(mp.ProducerError, match="initial inventory"):
        produce(setup, out, mode="convert")


def test_conflicting_duplicate_times_are_refused_anywhere(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    mu = declared(setup["spec"])
    write_run(out, times=PROFILE_T, pools=[(LABEL, lambda t: mu, 0.0)], gas=ZERO_GAS,
              duplicate=(4, 1e-3))                     # at 45 s, past every requested time
    with pytest.raises(mp.ProducerError, match="duplicate"):
        produce(setup, out, mode="convert")


def test_a_repeat_mass_that_differs_from_the_deck_is_refused(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    a, _e = masses()
    constant_run(out, setup["spec"], mass_override=a * (1 + 1e-6))
    with pytest.raises(mp.ProducerError, match="repeat-unit mass"):
        produce(setup, out, mode="convert")


def test_a_first_row_not_at_zero_is_refused(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], times=(-1.0, 0.0, 10.0, 20.0, 40.0, 45.0))
    with pytest.raises(mp.ProducerError, match="first row"):
        produce(setup, out, mode="convert")


def test_a_time_on_a_profile_row_takes_the_row_exactly(setup):
    out = setup["tmp"] / "out"
    produce(setup, out, end=20.0, snapshots=3)
    main, daughter, _h, _m = moving_run(out, setup["spec"])
    produce(setup, out, end=20.0, snapshots=3, mode="convert")
    entries = json.loads((out / "mom.json").read_text())
    assert [x["time_s"] for x in entries] == [0.0, 10.0, 20.0]
    for entry in entries:
        t = entry["time_s"]
        assert entry["quantities"]["polymer_derived_molecules_mol"] == main(t)[0] + daughter(t)[0]


def test_a_time_outside_the_profile_is_refused(setup):
    out = setup["tmp"] / "out"
    produce(setup, out, end=50.0, snapshots=3)
    constant_run(out, setup["spec"])
    with pytest.raises(mp.ProducerError, match="outside the profile"):
        produce(setup, out, end=50.0, snapshots=3, mode="convert")


def test_moment_columns_of_an_unknown_pool_are_refused(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], extra_columns=("other_mu0",))
    with pytest.raises(mp.ProducerError, match="belong to no pool"):
        produce(setup, out, mode="convert")


def test_moment_column_follows_the_mom_writer_rule():
    assert mp.moment_column("PS_mu0") == "PS_mu0"
    assert mp.moment_column("toylinearp_mu2") == "toylinearp_mu2"            # <= 16: verbatim
    long = DAUGHTER + "_mu1"
    assert mp.moment_column(long) == "toyli-" + hashlib.md5(long.encode()).hexdigest()[:6] + "_mu1"
    assert len(mp.moment_column(long)) == 16


# --------------------------------------------------------------------------
# quantity keys
# --------------------------------------------------------------------------
def test_quantity_keys_are_the_runners():
    keys = mp.runner_quantity_keys()
    assert keys == {"mu0_mol", "mu1_mol", "mu2_mol", "end_inventory_mol", "explicit_dp_mol",
                    "polymer_condensed_mass_g", "Mn_g_mol", "Mw_g_mol", "polymer_conversion", "released_mol",
                    "species_mol", "polymer_derived_molecules_mol"}
    mp.check_quantity_keys({"polymer_derived_molecules_mol": 1.0, "species_mol": {"H": 0.0}}, ["H", "H2"])
    with pytest.raises(mp.ProducerError, match="not keys of"):
        mp.check_quantity_keys({"polymer_derived_molecules_mol": 1.0, "Mn": 1.0}, [])
    with pytest.raises(mp.ProducerError, match="not produced"):
        mp.check_quantity_keys({"mu0_mol": 1.0}, [])
    with pytest.raises(mp.ProducerError, match="not declared"):
        mp.check_quantity_keys({"species_mol": {"CH4": 1.0}}, ["H", "H2"])


# --------------------------------------------------------------------------
# families
# --------------------------------------------------------------------------
def test_a_family_mismatch_is_refused_unless_allowed_and_then_recorded(setup):
    with pytest.raises(mp.ProducerError, match="differ from the catalogs' provenance families"):
        produce(setup, setup["tmp"] / "refused", families=["R_Recombination"])
    assert not (setup["tmp"] / "refused" / "rmg" / "input.py").exists()
    out = setup["tmp"] / "allowed"
    manifest = produce(setup, out, families=["R_Recombination"], allow_family_mismatch=True)
    record = manifest["family_check"]
    assert record["match"] is False and record["allowed_mismatch"] is True
    assert record["only_in_families"] == ["R_Recombination"] and record["only_in_catalogs"] == PE_FAMILIES
    conv = {c["name"]: c["value"] for c in json.loads((out / "conventions.json").read_text())["conventions"]}
    assert conv["family_check"] == record
    assert "kineticsFamilies=['R_Recombination']," in (out / "rmg" / "input.py").read_text()
    match = produce(setup, setup["tmp"] / "match", families=list(reversed(PE_FAMILIES)))
    assert match["family_check"]["match"] is True and "allowed_mismatch" not in match["family_check"]


def test_two_catalogs_give_the_union_of_their_families(setup):
    tmp = setup["tmp"]
    first = write_catalog(tmp / "a.json", toy_payloads([t for t in S.pe_templates() if t.name == "homolysis"]),
                          tag="a")
    second = write_catalog(tmp / "b.json", toy_payloads([t for t in S.pe_templates()
                                                         if t.name == "recombination"]), tag="b")
    manifest = produce(setup, tmp / "two", catalogs=[first, second],
                       families=["toy_homolysis", "toy_recombination"])
    assert manifest["family_check"]["catalog_families"] == ["toy_homolysis", "toy_recombination"]


# --------------------------------------------------------------------------
# round trip and the command line
# --------------------------------------------------------------------------
def test_a_produced_mom_file_is_accepted_by_the_runner(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    moving_run(out, setup["spec"])
    produce(setup, out, mode="convert")
    times = kr.observation_times(end=40.0, snapshots=9)
    mom, record = kr.load_mom(out / "mom.json", times)
    assert record["times"] == len(times) and all(q is not None for q in mom)


def test_command_line_writes_and_refuses(setup, capsys):
    out = setup["tmp"] / "cli"
    argv = [str(setup["deck"]), "--mom-checkout", CHECKOUT, "--temperature", "700", "--end", "40",
            "--snapshots", "9", "--families", *PE_FAMILIES, "--catalog", str(setup["catalog"]),
            "--out", str(out), "--database", setup["db"], "--write-only"]
    assert mp.main(argv) == 0
    assert (out / "rmg" / "input.py").is_file() and (out / "conventions.json").is_file()
    assert mp.main(argv) == 2                                  # a non-empty output directory
    assert "not empty" in capsys.readouterr().err
    bad = argv[:argv.index("--families") + 1] + ["R_Recombination"] + argv[argv.index("--catalog"):]
    bad[bad.index("--out") + 1] = str(setup["tmp"] / "cli2")
    assert mp.main(bad) == 2
    assert "--allow-family-mismatch" in capsys.readouterr().err
    constant_run(out, setup["spec"], admitted={"toy_homolysis": (1, 0)})
    convert = [x for x in argv if x != "--write-only"] + ["--convert"]
    assert mp.main(convert) == 2
    assert "inactive_chemistry" in capsys.readouterr().err
    assert mp.main(convert + ["--allow-inactive"]) == 0


def test_species_constraints_keep_the_declared_species_admissible(setup):
    out = setup["tmp"] / "out"
    manifest = produce(setup, out, species_constraints=["maximumRadicalElectrons=0"])
    text = (out / "rmg" / "input.py").read_text()
    assert "generatedSpeciesConstraints(allowed=['input species'], maximumRadicalElectrons=0)" in text
    assert manifest["facts"]["constraints"] == {"allowed": ["input species"], "maximumRadicalElectrons": 0}
    with pytest.raises(mp.ProducerError, match="KEY=VALUE"):
        produce(setup, setup["tmp"] / "bad", species_constraints=["maximumRadicalElectrons"])
