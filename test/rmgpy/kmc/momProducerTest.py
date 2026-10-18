"""Red tests: the method-of-moments producer for the runner's ``--mom`` file
(``rmgpy.kmc.mom_producer``, I-169).  No RMG run: the MoM side is a
hand-built solver profile and pool table with known answers.

Expected values are derived here without the producer: the moments of the
declared ensemble by hand, the capped-chain Mn / Mw / mass written out, the
piecewise-linear interpolant of samples of t^2 in closed form, and the
profile column of a long moment-dummy label by the MoM writer's rule.
"""

import csv
import hashlib
import json
import math
from pathlib import Path

import pytest

import melt_engine_support as S
from meltRunTest import pe_catalog_file, write_deck, write_catalog, toy_payloads

from rmgpy.kmc import mom_producer as mp
from rmgpy.kmc import run as kr
from rmgpy.polymer_spec import LinearHomopolymerSpec

CHECKOUT = "/nonexistent/mom-checkout"   # never run here; only recorded
PE_FAMILIES = ["toy_beta_scission", "toy_branching", "toy_disproportionation", "toy_h_abstraction",
               "toy_homolysis", "toy_recombination"]
LABEL = "toylinearp"                      # mom_label of "toy linear polymer (CH2-CH2)"
DAUGHTER = LABEL + "_rad_primary_end"     # a spawned pool with a long label

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


def masses():
    spec = LinearHomopolymerSpec.from_mapping(S.pe_mapping({10: 1, 20: 1}))
    return spec.repeat_unit.composition.molar_mass_g_mol(), spec.end_composition.molar_mass_g_mol()


@pytest.fixture
def setup(tmp_path):
    deck = write_deck(tmp_path / "pe_deck.py", S.pe_mapping({10: 1, 20: 1}))
    catalog = pe_catalog_file(tmp_path)
    (tmp_path / "db" / "input" / "kinetics" / "families").mkdir(parents=True)
    return {"deck": deck, "catalog": catalog, "db": str(tmp_path / "db"), "tmp": tmp_path}


def produce(setup, out, **kw):
    args = dict(checkout=CHECKOUT, temperature=700.0, end=40.0, snapshots=9, families=PE_FAMILIES,
                catalogs=[setup["catalog"]], out=out, database=setup["db"])
    args.update(kw)
    return mp.produce(setup["deck"], **args)


# --------------------------------------------------------------------------
# deck -> input file
# --------------------------------------------------------------------------
def test_deck_gives_the_golden_mom_input(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    text = (out / "rmg" / "input.py").read_text()
    start = text.index("polymer(\n")
    assert text[start:text.index("model(")] == GOLDEN
    # the golden moments are the declared ensemble's, by hand: one DP10 and one DP20 chain, 1 kg with ends
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
    for name in ("label", "legacy_translation", "initial_moments", "mass_basis", "capped_chains",
                 "moments_all_dp", "volume", "declared_species", "bath_gas", "pressure", "families",
                 "family_check", "database", "kernels_off", "no_cage", "model_settings", "profile",
                 "time_interpolation", "species_constraints"):
        assert name in names and names[name]["statement"], name
    assert names["time_interpolation"]["value"] == "linear in time"
    assert names["label"]["value"] == LABEL
    assert set(doc["not_produced"]) == {"end_inventory_mol", "explicit_dp_mol", "released_mol"}
    assert set(doc["quantity_sources"]) <= mp.runner_quantity_keys()
    assert set(doc["quantity_sources"]) | set(doc["not_produced"]) == mp.runner_quantity_keys()


# --------------------------------------------------------------------------
# profile -> --mom file
# --------------------------------------------------------------------------
PROFILE_T = [0.0, 10.0, 20.0, 40.0, 45.0]       # the solver's last step passes the end time


def write_profile(out, *, times=PROFILE_T, pools=(LABEL, DAUGHTER), extra_columns=()):
    """mu0 of the main pool samples 1 + (t/10)^2 (x0.5 for mu0), the daughter
    t^2/100, H 1e-3 t^2, H2 2e-3 t^2: known samples, interpolated in closed form."""
    rmg = out / "rmg"
    (rmg / "solver").mkdir(parents=True, exist_ok=True)
    (rmg / "chemkin").mkdir(parents=True, exist_ok=True)
    a, _e = masses()
    columns = {LABEL: [LABEL + "_mu0", LABEL + "_mu1", LABEL + "_mu2"],
               DAUGHTER: ["toyli-" + hashlib.md5((DAUGHTER + f"_mu{k}").encode()).hexdigest()[:6] + f"_mu{k}"
                          for k in range(3)]}
    header = ["Time (s)", "Volume (m^3)", "N2"] + [c for p in pools for c in columns[p]] + \
             ["H(1)", "H2(2)", LABEL + "(3)"] + list(extra_columns)
    rows = []
    for t in times:
        main = 0.5 * (1.0 + (t / 10.0) ** 2)
        daughter = t * t / 100.0
        values = {LABEL: [main, 30 * main, 900 * main], DAUGHTER: [daughter, 5 * daughter, 40 * daughter]}
        rows.append([t, 1.0, 1.0] + [v for p in pools for v in values[p]] + [1e-3 * t * t, 2e-3 * t * t, 7.0]
                    + [0.0] * len(extra_columns))
    with open(rmg / "solver" / "simulation_1_8.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)
    with open(rmg / "solver" / "simulation_1_7.csv", "w", newline="") as handle:   # an earlier core: ignored
        csv.writer(handle).writerow(["Time (s)", "Volume (m^3)"])
    pools_doc = {"schema_version": "2.0", "pools": [{"label": p, "monomer_mw_g_mol": a} for p in pools]}
    (rmg / "chemkin" / "polymer_pools.json").write_text(json.dumps(pools_doc))


def interpolant(f, t):
    """The piecewise-linear interpolant of f sampled at PROFILE_T."""
    for lo, hi in zip(PROFILE_T, PROFILE_T[1:]):
        if lo <= t <= hi:
            return f(lo) + (f(hi) - f(lo)) * (t - lo) / (hi - lo)
    raise AssertionError(t)


def test_profile_converts_to_the_mom_file_at_the_observation_times(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    write_profile(out)
    produce(setup, out, mode="convert")
    entries = json.loads((out / "mom.json").read_text())
    times = kr.observation_times(end=40.0, snapshots=9)
    assert [x["time_s"] for x in entries] == times == [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0]
    a, e = masses()
    mass0 = None
    for entry in entries:
        t, q = entry["time_s"], entry["quantities"]
        main = interpolant(lambda s: 0.5 * (1.0 + (s / 10.0) ** 2), t)
        daughter = interpolant(lambda s: s * s / 100.0, t)
        mu0, mu1, mu2 = main + daughter, 30 * main + 5 * daughter, 900 * main + 40 * daughter
        # closed form at 25 s: between the 20 s and 40 s samples of t^2, (20 + 40) 25 - 800 = 700
        if t == 25.0:
            assert daughter == pytest.approx(7.0, rel=1e-15)
        mass = a * mu1 + e * mu0
        mass0 = mass if mass0 is None else mass0
        expected = {"mu0_mol": mu0, "mu1_mol": mu1, "mu2_mol": mu2,
                    "Mn_g_mol": mass / mu0,
                    "Mw_g_mol": (a * a * mu2 + 2 * a * e * mu1 + e * e * mu0) / mass,
                    "polymer_condensed_mass_g": mass, "polymer_conversion": 1.0 - mass / mass0,
                    "polymer_derived_molecules_mol": mu0}
        for key, value in expected.items():
            assert q[key] == pytest.approx(value, rel=1e-12, abs=1e-15), (t, key)
        assert q["species_mol"] == pytest.approx({"H": interpolant(lambda s: 1e-3 * s * s, t),
                                                  "H2": interpolant(lambda s: 2e-3 * s * s, t)}, rel=1e-12)
        assert not set(q) & set(mp.NOT_PRODUCED)
    manifest = json.loads((out / "producer.json").read_text())
    assert manifest["read"]["profile"].endswith("simulation_1_8.csv")
    assert manifest["read"]["pool_labels"] == [LABEL, DAUGHTER]


def test_a_time_on_a_profile_row_takes_the_row_exactly(setup):
    out = setup["tmp"] / "out"
    produce(setup, out, end=20.0, snapshots=3)
    write_profile(out)
    produce(setup, out, end=20.0, snapshots=3, mode="convert")
    entries = json.loads((out / "mom.json").read_text())
    assert [x["time_s"] for x in entries] == [0.0, 10.0, 20.0]
    for entry, main in zip(entries, (0.5, 1.0, 2.5)):
        daughter = entry["time_s"] ** 2 / 100.0
        assert entry["quantities"]["mu0_mol"] == main + daughter        # bit for bit: no interpolation


def test_a_time_outside_the_profile_is_refused(setup):
    out = setup["tmp"] / "out"
    produce(setup, out, end=50.0, snapshots=3)
    write_profile(out)
    with pytest.raises(mp.ProducerError, match="outside the profile"):
        produce(setup, out, end=50.0, snapshots=3, mode="convert")
    late = setup["tmp"] / "late"
    produce(setup, late)
    write_profile(late, times=[1.0, 10.0, 20.0, 45.0])
    with pytest.raises(mp.ProducerError, match="outside the profile|first row"):
        produce(setup, late, mode="convert")


def test_moment_columns_of_an_unknown_pool_are_refused(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    write_profile(out, extra_columns=("other_mu0",))
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
    mp.check_quantity_keys({"mu0_mol": 1.0, "species_mol": {"H": 0.0}}, ["H", "H2"])
    with pytest.raises(mp.ProducerError, match="not keys of"):
        mp.check_quantity_keys({"mu0_mol": 1.0, "Mn": 1.0}, [])
    with pytest.raises(mp.ProducerError, match="not produced"):
        mp.check_quantity_keys({"end_inventory_mol": 1.0}, [])
    with pytest.raises(mp.ProducerError, match="not declared"):
        mp.check_quantity_keys({"species_mol": {"CH4": 1.0}}, ["H", "H2"])


def test_produced_keys_are_declared_and_partition_the_runners():
    assert set(mp.QUANTITY_SOURCES).isdisjoint(mp.NOT_PRODUCED)
    assert set(mp.QUANTITY_SOURCES) | set(mp.NOT_PRODUCED) == mp.runner_quantity_keys()


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
    # matching families, given in any order
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
    write_profile(out)
    produce(setup, out, mode="convert")
    times = kr.observation_times(end=40.0, snapshots=9)
    mom, record = kr.load_mom(out / "mom.json", times)
    assert record["times"] == len(times) and all(q is not None for q in mom)
    assert mom[2]["mu0_mol"] == json.loads((out / "mom.json").read_text())[2]["quantities"]["mu0_mol"]


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


def test_species_constraints_keep_the_declared_species_admissible(setup):
    out = setup["tmp"] / "out"
    manifest = produce(setup, out, species_constraints=["maximumRadicalElectrons=0"])
    text = (out / "rmg" / "input.py").read_text()
    assert "generatedSpeciesConstraints(allowed=['input species'], maximumRadicalElectrons=0)" in text
    assert manifest["facts"]["constraints"] == {"allowed": ["input species"], "maximumRadicalElectrons": 0}
    with pytest.raises(mp.ProducerError, match="KEY=VALUE"):
        produce(setup, setup["tmp"] / "bad", species_constraints=["maximumRadicalElectrons"])
