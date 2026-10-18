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


def make_setup(tmp_path, counts=BASE, *, edit=None, **kw):
    """``edit``: a function that changes the declaration mapping in place."""
    mapping = S.pe_mapping(counts, **kw)
    if edit is not None:
        edit(mapping)
    deck = write_deck(tmp_path / "pe_deck.py", mapping)
    catalog = pe_catalog_file(tmp_path)
    (tmp_path / "db" / "input" / "kinetics" / "families").mkdir(parents=True, exist_ok=True)
    return {"deck": deck, "catalog": catalog, "db": str(tmp_path / "db"), "tmp": tmp_path,
            "spec": LinearHomopolymerSpec.from_mapping(mapping)}


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
def proxy_of(label=LABEL):
    return f"{label}(3)"


def family_rows(admitted, label=LABEL):
    """{family: (rows, admitted rows)} -> pool-coupled reaction rows, each
    with its own stoichiometry (proxy <=> proxy + X<k>)."""
    rows = []
    for family, (n, ok) in sorted(admitted.items()):
        for j in range(n):
            rows.append({"family": family, "admitted": j < ok, "reactants": [proxy_of(label)],
                         "products": [proxy_of(label), f"X{len(rows)}"], "reversible": True})
    return rows


def write_run(out, *, times, pools, gas=(), admitted=None, rows=None, label=LABEL, end_groups=("[H]", "[H]"),
              mass_override=None, duplicate=None, extra_columns=(), n2=lambda t: 1.0, gas_reactions=()):
    """A MoM output directory under ``out/rmg``.

    ``pools``: [(label, f(t) -> (mu0, mu1, mu2), defect_g_mol)]; ``gas``:
    [(identifier, adjacency, f(t) -> mol)]; ``admitted``: {family: (rows,
    admitted rows)} of pool-coupled reactions (default one admitted
    toy_homolysis row); ``rows``: the reaction rows themselves (dicts with
    family, admitted, reactants, products, reversible and optionally
    cantera_null / expanded: no Cantera index / a MultiArrhenius the Chemkin
    writer expands without its family annotation; A: the rate coefficient's
    A, default 1; no_chemkin: the row has no Chemkin entry);
    ``gas_reactions``: Chemkin-only reactions (reactants, products,
    reversible, A)."""
    rmg = out / "rmg"
    for sub in ("solver", "chemkin"):
        (rmg / sub).mkdir(parents=True, exist_ok=True)
    a, _e = masses()
    moment_ids = {p: [mp.moment_column(f"{p}_mu{k}") for k in range(3)] for p, _f, _d in pools}
    proxy = f"{label}(3)"
    header = ["Time (s)", "Volume (m^3)", "N2"] + [c for p, _f, _d in pools for c in moment_ids[p]] + \
        [g[0] for g in gas] + [proxy] + list(extra_columns)
    table = []
    for t in times:
        row = [t, 1.0, n2(t)]
        for _p, f, _d in pools:
            row += list(f(t))
        row += [g[2](t) for g in gas] + [7.0] + [0.0] * len(extra_columns)
        table.append(row)
    if duplicate is not None:
        k, delta = duplicate
        bad = list(table[k])
        bad[3] += delta
        table.insert(k + 1, bad)
    with open(rmg / "solver" / "simulation_1_8.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(table)
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
    if rows is None:
        rows = family_rows({"toy_homolysis": (1, 1)} if admitted is None else admitted, label)
    reactions, annotated, k, c = [], ["REACTIONS    KCAL/MOLE   MOLES", ""], 0, 0
    pool_labels = {p for p, _f, _d in pools}

    def is_proxy(x):
        return x.rsplit("(", 1)[0] in pool_labels and x.endswith(")")

    for r in rows:
        arrow = "<=>" if r["reversible"] else "=>"
        equation = " + ".join(r["reactants"]) + f" {arrow} " + " + ".join(r["products"])
        if r.get("cantera_null"):
            row = {"id": f"d{r['family']}:{equation}:0", "cantera": None}
        else:
            row = {"id": f"r{k}", "cantera": {"index": k, "equation": equation}}
            k += 1
        big_a = r.get("A", 1.0)
        row.update({"kinetics": {"A": big_a, "n": 0.0, "Ea": 0.0, "reversible": r["reversible"]},
                    "reactants": list(r["reactants"]), "products": list(r["products"]),
                    "proxy_reactants": [x for x in r["reactants"] if is_proxy(x)],
                    "proxy_products": [x for x in r["products"] if is_proxy(x)],
                    "refused": not r["admitted"]})
        if not r["admitted"]:
            row["refused_reason"] = "conduit-deferred"
        reactions.append(row)
        line = "+".join(r["reactants"]) + arrow + "+".join(r["products"])
        for _ in range(0 if r.get("no_chemkin") else 2 if r.get("expanded") else 1):
            c += 1
            annotated.append(f"! Reaction index: Chemkin #{c}; RMG #{len(reactions)}")
            if not r.get("expanded"):
                annotated.append(f"! Template reaction: {r['family']}")
            annotated += [f"{line:<51} {big_a} 0.0 0.0"] + (["DUPLICATE"] if r.get("expanded") else []) + [""]
    for j, g in enumerate(gas_reactions):
        c += 1
        line = "+".join(g["reactants"]) + ("<=>" if g["reversible"] else "=>") + "+".join(g["products"])
        annotated += [f"! Reaction index: Chemkin #{c}; RMG #{len(reactions) + 1 + j}",
                      "! Template reaction: toy_h_abstraction", f"{line:<51} {g['A']} 0.0 0.0", ""]
    (rmg / "chemkin" / "polymer_pools.json").write_text(json.dumps(
        {"schema_version": "2.0", "pools": pool_rows, "reactions": reactions,
         "conventions": {"stale_topology": False}}))
    (rmg / "chemkin" / "chem_annotated.inp").write_text("\n".join(annotated + ["END", ""]))


def mom_entries(out):
    """The entries of a produced --mom file (header form)."""
    doc = json.loads((out / "mom.json").read_text())
    assert set(doc) == {"header", "entries"}
    return doc["entries"]


WHOLE = ("whole_mu0_mol", "whole_mu1_mol", "whole_mu2_mol", "whole_Mn_g_mol", "whole_Mw_g_mol")


def declared(spec):
    return [float(x) for x in spec.macroscopic_moments_mol]


ZERO_GAS = [("H(1)", ADJ["H"], lambda t: 0.0), ("H2(2)", ADJ["H2"], lambda t: 0.0)]


def drift(mu):
    """Pool moments whose chain count grows (scission-like) at fixed mu1: the
    profile moves, the booked mass (ledger) stays 1000 g."""
    return lambda t: (mu[0] * (1 + 1e-3 * t), mu[1], mu[2])


def constant_run(out, spec, times=(0.0, 10.0, 20.0, 40.0, 45.0), flat=False, **kw):
    """A run whose mass stays put; ``flat``: nothing in the profile moves."""
    mu = declared(spec)
    write_run(out, times=list(times), pools=[(LABEL, (lambda t: mu) if flat else drift(mu), 0.0)],
              gas=ZERO_GAS, **kw)


PROFILE_T = [0.0, 10.0, 20.0, 40.0, 45.0]       # the solver's last step passes the end time


def interpolant(f, t):
    """The piecewise-linear interpolant of f sampled at PROFILE_T."""
    for lo, hi in zip(PROFILE_T, PROFILE_T[1:]):
        if lo <= t <= hi:
            return f(lo) + (f(hi) - f(lo)) * (t - lo) / (hi - lo)
    raise AssertionError(t)


def moving_run(out, spec, **kw):
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
              gas=[("H(1)", ADJ["H"], lambda t: 0.0), ("H2(2)", ADJ["H2"], lambda t: h * t * t)], **kw)
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
                                        "polymer_derived_molecules_mol", "species_mol", *WHOLE}
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
    quantities = mom_entries(out)[0]["quantities"]
    engine = me.MeltEngine(s["spec"], S.ToyCatalog(names=PE_LIVE), 700.0, np.random.default_rng(1))
    obs = mo.MeltObserver(s["spec"], xs=4, identity=engine.identity).observe(engine)
    rows = [r for r in mo.side_by_side(obs, quantities) if r["mom"] is not None]
    assert {r["quantity"] for r in rows} >= {"polymer_condensed_mass_g", "polymer_derived_molecules_mol",
                                             "polymer_conversion", "species_mol.H", "species_mol.H2", *WHOLE}
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
    entries = mom_entries(out)
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
    entries = mom_entries(out)
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
    for entry in mom_entries(out):
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
    entries = mom_entries(out)
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
    q = mom_entries(out)[-1]["quantities"]
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
    entries = mom_entries(out)
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
                    "species_mol", "polymer_derived_molecules_mol", *WHOLE}
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


# --------------------------------------------------------------------------
# rework 2 (r127)
# --------------------------------------------------------------------------
def test_gas_oligomers_are_booked_on_the_side_kmc_puts_them(tmp_path):
    """Spar's fixture: PE DP6 -> DP3 (kept in the pool) + hexene in MoM's gas
    phase.  kMC counts the undeclared hexene (n = 3 >= 1) as polymer-derived:
    the mass stays 1000 g, the molecules double, the conversion is 0."""
    s = make_setup(tmp_path, {6: 1})
    out = tmp_path / "out"
    produce(s, out, end=10.0, snapshots=2)
    u, mu1, mu2 = declared(s["spec"])

    def pool(t):
        f = t / 10.0
        return u, mu1 - 3 * f * u, mu2 * (1 - f) + 9 * f * u

    hexene = Molecule(smiles="C=CCCCC").to_adjacency_list()
    write_run(out, times=[0.0, 10.0], pools=[(LABEL, pool, 0.0)],
              gas=ZERO_GAS + [("C6H12(4)", hexene, lambda t: u * t / 10.0)])
    produce(s, out, end=10.0, snapshots=2, mode="convert")
    q = mom_entries(out)[-1]["quantities"]
    assert q["polymer_condensed_mass_g"] == pytest.approx(1000.0, rel=1e-12)
    assert q["polymer_conversion"] == pytest.approx(0.0, abs=1e-12)
    assert q["polymer_derived_molecules_mol"] == pytest.approx(2 * u, rel=1e-12)
    manifest = json.loads((out / "producer.json").read_text())
    assert manifest["read"]["species_classes"]["C6H12(4)"] == "polymer"


def _methyl_capped(mapping):
    mapping["end_groups"]["head"] = {"atoms": [S._atom("m", hydrogens=3)], "bonds": [],
                                     "ports": {"cap": {"atom": "m", "bond_order": 1}},
                                     "electronic_state": {"multiplicity": 1}}


def test_a_species_kmc_cannot_classify_exports_no_polymer_rows(tmp_path):
    """Methyl-capped PE: carbon is on the backbone and in an end group, so
    whether a released propane carries a whole repeat unit's backbone (kMC:
    polymer-derived) depends on ancestry MoM does not keep.  The polymer rows
    are not exported where it is present, even though the booking closes."""
    s = make_setup(tmp_path, BASE, edit=_methyl_capped)
    out = tmp_path / "out"
    produce(s, out)
    mu = declared(s["spec"])
    a, _e = masses()
    m3 = mw("CCC")
    x = 1e-4

    def pool(t):
        return mu[0], mu[1] - m3 * x * t / a, mu[2]

    propane = Molecule(smiles="CCC").to_adjacency_list()
    write_run(out, times=PROFILE_T, pools=[(LABEL, pool, 0.0)], end_groups=("[CH3]", "[H]"),
              gas=ZERO_GAS + [("C3H8(4)", propane, lambda t: x * t)])
    produce(s, out, mode="convert")
    entries = mom_entries(out)
    assert {"polymer_condensed_mass_g", "polymer_derived_molecules_mol", *WHOLE} <= set(entries[0]["quantities"])
    diagnostics = json.loads((out / "mom_diagnostics.json").read_text())
    for entry, diag in list(zip(entries, diagnostics))[1:]:
        for key in ("polymer_condensed_mass_g", "polymer_conversion", "polymer_derived_molecules_mol", *WHOLE):
            assert key not in entry["quantities"]
        assert diag["mass_closure_residual_g"] == pytest.approx(0.0, abs=1e-9)
        assert diag["unclassified_species"] == ["C3H8(4)"]


def test_zero_flux_chemistry_is_inactive(setup):
    """17 refused homolysis rows and one admitted irreversible H abstraction
    that needs H, which nothing admitted makes: inactive.  With one homolysis
    row admitted (it makes H) both rows are live."""
    proxy = proxy_of()
    hab = {"family": "toy_h_abstraction", "admitted": True, "reactants": ["H(1)", proxy],
           "products": ["H2(2)", proxy], "reversible": False}

    def homolysis(admitted):
        return {"family": "toy_homolysis", "admitted": admitted, "reactants": [proxy],
                "products": [proxy, "H(1)"], "reversible": False}

    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], rows=[homolysis(False)] * 17 + [hab])
    with pytest.raises(mp.ProducerError, match="inactive_chemistry"):
        produce(setup, out, mode="convert")
    census = json.loads((out / "producer.json").read_text())["census"]
    assert census["admitted_in_requested_families"] == 1 and census["live_rows"] == []
    live = setup["tmp"] / "live"
    produce(setup, live)
    constant_run(live, setup["spec"], rows=[homolysis(True)] + [homolysis(False)] * 16 + [hab])
    manifest = produce(setup, live, mode="convert")
    assert "inactive_chemistry" not in manifest["flags"]
    assert manifest["census"]["live_rows"] == ["r0", "r17"]


def test_a_cantera_null_row_is_counted_by_its_serialized_family(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    proxy = proxy_of()
    constant_run(out, setup["spec"], rows=[{"family": "toy_homolysis", "admitted": True, "reactants": [proxy],
                                            "products": [proxy, "X0"], "reversible": True,
                                            "cantera_null": True}])
    manifest = produce(setup, out, mode="convert")
    assert manifest["census"]["families"]["toy_homolysis"] == {"rows": 1, "admitted": 1}


def test_an_expanded_multi_arrhenius_row_is_counted_by_name(setup):
    """The Chemkin writer expands a template reaction's MultiArrhenius into
    plain entries without the family: counted by its equation, as a row of
    the requested families (the input admits no other source)."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    proxy = proxy_of()
    constant_run(out, setup["spec"], rows=[{"family": "toy_homolysis", "admitted": True, "reactants": [proxy],
                                            "products": [proxy, "X0"], "reversible": True, "expanded": True}])
    manifest = produce(setup, out, mode="convert")
    census = manifest["census"]
    assert census["unattributed"] == {f"{proxy} <=> {proxy} + X0": {"rows": 1, "admitted": 1}}
    assert census["admitted_in_requested_families"] == 1 and "inactive_chemistry" not in manifest["flags"]


def test_the_bath_does_not_dilute_the_mass_closure(tmp_path):
    """1 ug of PE in 1 mol N2: a declared butane appears while the pool keeps
    its mass (the exported polymer mass would be 3.59 % high).  Against the
    bath it closes to 1e-9; against the kMC inventory it does not."""
    def tiny(mapping):
        mapping["initial_mass_kg"] = 1.0e-9

    s = make_setup(tmp_path, BASE, edit=tiny, small={"butane": S.BUTANE})
    out = tmp_path / "out"
    produce(s, out)
    mu = declared(s["spec"])
    b = 0.0359e-6 / mw("CCCC")
    butane = Molecule(smiles="CCCC").to_adjacency_list()
    write_run(out, times=PROFILE_T, pools=[(LABEL, lambda t: mu, 0.0)],
              gas=ZERO_GAS + [("C4H10(4)", butane, lambda t: b * t / 40.0)])
    produce(s, out, mode="convert")
    entries = mom_entries(out)
    assert "polymer_condensed_mass_g" in entries[0]["quantities"]
    for entry in entries[1:]:
        assert "polymer_condensed_mass_g" not in entry["quantities"]
    assert entries[-1]["quantities"]["species_mol"]["butane"] == pytest.approx(b, rel=1e-12)


def test_a_bath_that_does_not_close_exports_no_mass(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    mu = declared(setup["spec"])
    write_run(out, times=PROFILE_T, pools=[(LABEL, lambda t: mu, 0.0)], gas=ZERO_GAS,
              n2=lambda t: 1.0 - 1e-3 * t)
    produce(setup, out, mode="convert")
    entries = mom_entries(out)
    diagnostics = json.loads((out / "mom_diagnostics.json").read_text())
    for entry, diag in list(zip(entries, diagnostics))[1:]:
        assert "polymer_condensed_mass_g" not in entry["quantities"]
        assert diag["bath_closure_residual_g"] == pytest.approx(-1e-3 * entry["time_s"] * mw("N#N"), rel=1e-9)


def test_the_initial_mobile_inventory_is_checked(tmp_path):
    s = make_setup(tmp_path, BASE, small={"ethylene": S.ETHYLENE}, mobile={"ethylene": 3})
    out = tmp_path / "out"
    manifest = produce(s, out)
    amount = manifest["facts"]["initial_moles"]["ethylene"]
    assert amount > 0
    mu = declared(s["spec"])
    ethylene = Molecule(smiles="C=C").to_adjacency_list()
    write_run(out, times=PROFILE_T, pools=[(LABEL, lambda t: mu, 0.0)],
              gas=ZERO_GAS + [("C2H4(4)", ethylene, lambda t: 0.0)])
    with pytest.raises(mp.ProducerError, match="initial inventory"):
        produce(s, out, mode="convert")
    write_run(out, times=PROFILE_T, pools=[(LABEL, drift(mu), 0.0)],
              gas=ZERO_GAS + [("C2H4(4)", ethylene, lambda t: amount)])
    q = produce(s, out, mode="convert") and mom_entries(out)[-1]["quantities"]
    assert q["species_mol"]["ethylene"] == amount


def test_a_resonance_form_of_a_declared_species_resolves(tmp_path):
    """Declared [CH2]C=CC; RMG's dictionary holds the form C=C[CH]C."""
    s = make_setup(tmp_path, BASE, small={"crotyl": S.CROTYL})
    out = tmp_path / "out"
    produce(s, out)
    mu = declared(s["spec"])
    other = Molecule(smiles="C=C[CH]C")
    assert not other.is_isomorphic(Molecule(smiles="[CH2]C=CC"))
    write_run(out, times=PROFILE_T, pools=[(LABEL, lambda t: mu, 0.0)],
              gas=ZERO_GAS + [("C4H7(4)", other.to_adjacency_list(), lambda t: 1e-6 * t)])
    produce(s, out, mode="convert")
    q = mom_entries(out)[-1]["quantities"]
    assert q["species_mol"]["crotyl"] == pytest.approx(4e-5, rel=1e-12)


def test_initial_moments_off_by_a_part_per_million_are_refused(setup):
    """Between INITIAL_RELATIVE (1e-9) and a loose bound such as 0.1."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    mu = declared(setup["spec"])
    write_run(out, times=PROFILE_T, pools=[(LABEL, lambda t: [x * (1 + 1e-6) for x in mu], 0.0)], gas=ZERO_GAS)
    with pytest.raises(mp.ProducerError, match="initial inventory"):
        produce(setup, out, mode="convert")


def test_a_doubled_initial_mu2_is_refused(setup):
    """mu0 and mu1 (hence the booked mass) right, mu2 doubled: only the moment
    check sees it."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    mu = declared(setup["spec"])
    write_run(out, times=PROFILE_T, pools=[(LABEL, lambda t: [mu[0], mu[1], 2 * mu[2]], 0.0)], gas=ZERO_GAS)
    with pytest.raises(mp.ProducerError, match="initial inventory"):
        produce(setup, out, mode="convert")


def test_an_inert_species_kmc_cannot_hold_blocks_nothing_while_absent(setup):
    """RMG puts Ar in every core; the kMC melt has no argon: unclassifiable,
    harmless at zero amount."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    mu = declared(setup["spec"])
    write_run(out, times=PROFILE_T, pools=[(LABEL, drift(mu), 0.0)],
              gas=ZERO_GAS + [("Ar", "1 Ar u0 p4 c0\n", lambda t: 0.0)])
    manifest = produce(setup, out, mode="convert")
    assert manifest["read"]["species_classes"]["Ar"] == "undetermined"
    q = mom_entries(out)[-1]["quantities"]
    assert q["polymer_condensed_mass_g"] == pytest.approx(1000.0, rel=1e-12)


# --------------------------------------------------------------------------
# rework 3 (r129)
# --------------------------------------------------------------------------
def _refused_census(setup, out):
    with pytest.raises(mp.ProducerError, match="inactive_chemistry"):
        produce(setup, out, mode="convert")
    return json.loads((out / "producer.json").read_text())["census"]


def test_a_zero_rate_row_is_dead(setup):
    """An admitted homolysis with A = 0 on a moving profile: not live."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    proxy = proxy_of()
    constant_run(out, setup["spec"], rows=[{"family": "toy_homolysis", "admitted": True, "reactants": [proxy],
                                            "products": [proxy, "X0"], "reversible": True, "A": 0.0}])
    census = _refused_census(setup, out)
    assert census["live_rows"] == [] and census["zero_rate_rows"] == ["r0"]


def test_a_zero_rate_gas_reaction_unlocks_nothing(setup):
    """The admitted H abstraction needs H; the only maker of H is a gas
    reaction from the bath: with A = 0 it makes nothing, with A = 1 it does."""
    proxy = proxy_of()
    hab = [{"family": "toy_h_abstraction", "admitted": True, "reactants": ["H(1)", proxy],
            "products": ["H2(2)", proxy], "reversible": False}]
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], rows=hab,
                 gas_reactions=[{"reactants": ["N2"], "products": ["H(1)"], "reversible": False, "A": 0.0}])
    assert _refused_census(setup, out)["live_rows"] == []
    live = setup["tmp"] / "live"
    produce(setup, live)
    constant_run(live, setup["spec"], rows=hab,
                 gas_reactions=[{"reactants": ["N2"], "products": ["H(1)"], "reversible": False, "A": 1.0}])
    assert produce(setup, live, mode="convert")["census"]["live_rows"] == ["r0"]


def test_a_flat_profile_is_inactive_whatever_the_census(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], flat=True)
    census = _refused_census(setup, out)
    assert census["live_rows"] == ["r0"] and census["profile_changes"] is False


def test_a_reversible_row_runs_from_its_products(setup):
    """Survivor reversible-not-run: products present at t = 0, reactants absent."""
    proxy = proxy_of()

    def row(reversible):
        return [{"family": "toy_homolysis", "admitted": True, "reactants": ["X5"], "products": [proxy],
                 "reversible": reversible}]

    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], rows=row(True))
    assert produce(setup, out, mode="convert")["census"]["live_rows"] == ["r0"]
    one_way = setup["tmp"] / "one_way"
    produce(setup, one_way)
    constant_run(one_way, setup["spec"], rows=row(False))
    assert _refused_census(setup, one_way)["live_rows"] == []


def test_a_proxy_of_an_empty_pool_is_absent(setup):
    """Survivor pool-mu0-zero: the daughter pool has mu0 = 0 at t = 0; a row
    that needs its proxy is not live."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    proxy = f"{DAUGHTER}(5)"
    moving_run(out, setup["spec"], rows=[{"family": "toy_homolysis", "admitted": True, "reactants": [proxy],
                                          "products": [proxy, "X0"], "reversible": False}])
    assert _refused_census(setup, out)["live_rows"] == []


def test_duplicate_stoichiometries_keep_their_own_families(setup):
    """Spar's fixture: one refused toy_homolysis row and one admitted
    toy_beta_scission row with the same species."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    proxy = proxy_of()
    rows = [{"family": family, "admitted": ok, "reactants": [proxy], "products": [proxy, "X0"],
             "reversible": True} for family, ok in (("toy_homolysis", False), ("toy_beta_scission", True))]
    constant_run(out, setup["spec"], rows=rows)
    census = produce(setup, out, mode="convert")["census"]
    assert census["families"] == {"toy_homolysis": {"rows": 1, "admitted": 0},
                                  "toy_beta_scission": {"rows": 1, "admitted": 1}}


def test_rows_without_a_one_to_one_chemkin_match_are_refused(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    proxy = proxy_of()
    rows = [{"family": "toy_homolysis", "admitted": True, "reactants": [proxy], "products": [proxy, "X0"],
             "reversible": True, "no_chemkin": j == 1} for j in range(2)]
    constant_run(out, setup["spec"], rows=rows)
    with pytest.raises(mp.ProducerError, match="one-to-one"):
        produce(setup, out, mode="convert")


def test_toluene_in_a_ps_deck_is_undetermined_and_withholds_the_polymer_rows(tmp_path):
    """Pins the conservative bracket (r129 item 3): PS carbon is backbone and
    ring, so a released toluene is light or polymer-derived by ancestry MoM
    does not keep.  A later, tighter bracket must change this test."""
    from meltStateTest import ps_mapping
    mapping = ps_mapping({10: 1, 20: 1}, mobile={})
    mapping["small_species"] = {}
    mapping["transport_model"] = {"met": S.ps_met()}
    deck = write_deck(tmp_path / "ps_deck.py", mapping)
    spec = LinearHomopolymerSpec.from_mapping(mapping)
    s = {"deck": deck, "catalog": pe_catalog_file(tmp_path), "db": str(tmp_path / "db"), "tmp": tmp_path,
         "spec": spec}
    (tmp_path / "db" / "input" / "kinetics" / "families").mkdir(parents=True, exist_ok=True)
    out = tmp_path / "out"
    produce(s, out)
    a = spec.repeat_unit.composition.molar_mass_g_mol()
    mu = declared(spec)
    m7, x = mw("Cc1ccccc1"), 1e-4
    toluene = Molecule(smiles="Cc1ccccc1").to_adjacency_list()
    write_run(out, times=PROFILE_T, label="PS", pools=[("PS", lambda t: (mu[0], mu[1] - m7 * x * t / a, mu[2]),
                                                        0.0)],
              mass_override=a, gas=ZERO_GAS + [("C7H8(4)", toluene, lambda t: x * t)])
    manifest = produce(s, out, mode="convert")
    assert manifest["read"]["species_classes"]["C7H8(4)"] == "undetermined"
    entries = mom_entries(out)
    diagnostics = json.loads((out / "mom_diagnostics.json").read_text())
    for entry, diag in list(zip(entries, diagnostics))[1:]:
        for key in ("polymer_condensed_mass_g", "polymer_conversion", "polymer_derived_molecules_mol"):
            assert key not in entry["quantities"]
        assert diag["unclassified_species"] == ["C7H8(4)"]
        assert diag["mass_closure_residual_g"] == pytest.approx(0.0, abs=1e-9)


# --------------------------------------------------------------------------
# I-170: whole-ensemble moments and the header form
# --------------------------------------------------------------------------
def test_whole_ensemble_keys_are_the_pools_plus_polymer_class_species(tmp_path):
    """Pool moments plus the hexene (kMC class polymer, n = 3): mu_k + x 3^k;
    Mn / Mw of those moments with the deck's a and e (capped chains)."""
    from rmgpy.polymer_mass import observables_from_moments
    s = make_setup(tmp_path, {6: 1})
    out = tmp_path / "out"
    produce(s, out, end=10.0, snapshots=2)
    u, mu1, mu2 = declared(s["spec"])

    def pool(t):
        f = t / 10.0
        return u, mu1 - 3 * f * u, mu2 * (1 - f) + 9 * f * u

    hexene = Molecule(smiles="C=CCCCC").to_adjacency_list()
    write_run(out, times=[0.0, 10.0], pools=[(LABEL, pool, 0.0)],
              gas=ZERO_GAS + [("C6H12(4)", hexene, lambda t: u * t / 10.0)])
    produce(s, out, end=10.0, snapshots=2, mode="convert")
    q = mom_entries(out)[-1]["quantities"]
    whole = [u + u, mu1 - 3 * u + 3 * u, 9 * u + 9 * u]
    assert [q["whole_mu0_mol"], q["whole_mu1_mol"], q["whole_mu2_mol"]] == pytest.approx(whole, rel=1e-12)
    a, e = masses()
    conv = observables_from_moments(whole, a, e)
    assert q["whole_Mn_g_mol"] == pytest.approx(conv.mn_g_mol, rel=1e-12)
    assert q["whole_Mw_g_mol"] == pytest.approx(conv.mw_g_mol, rel=1e-12)


def test_the_mom_file_carries_the_flags_and_the_conventions_digest(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], admitted={"toy_homolysis": (3, 0)})
    manifest = produce(setup, out, mode="convert", allow_inactive=True)
    doc = json.loads((out / "mom.json").read_text())
    header = doc["header"]
    assert sorted(header["flags"]) == ["inactive_chemistry", "kinetics_not_comparable"]
    assert header["flags"] == manifest["flags"]
    assert header["conventions_sha256"] == hashlib.sha256((out / "conventions.json").read_bytes()).hexdigest()
    assert header["producer"]["schema"] == mp.SCHEMA and header["producer"]["deck_sha256"] == manifest["deck"]["sha256"]
    mom, record = kr.load_mom(out / "mom.json", kr.observation_times(end=40.0, snapshots=9))
    assert record["header"] == header and all(q is not None for q in mom)
