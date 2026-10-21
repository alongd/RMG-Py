"""Red tests (I-214, S2): the kernel references -- an exact length-resolved
population balance (V0) and a length-level stochastic simulation (V1) of the
kernel set a kernel-rates table states (``rmgpy.kmc.kernel_reference``).

The null tests of the shared-chemistry design (``/home/alon/runs/i204/
design.md`` §4): recombination only and phi = 0 are flat; random scission is
Montroll and Simha's solution; PE homolysis alone is its exact linear master
equation, with the cut law read from the melt engine's own transitions on
the declared chain.  The stochastic run agrees with the balance within 3
standard errors on a bilinear toy and on PE homolysis; a depository override
is honoured at its exact sizes.
"""

import base64
import json
import math
import zlib
from collections import Counter
from pathlib import Path

import numpy
import pytest

import melt_engine_support as S
from meltRunTest import PE_LIVE, pe_catalog_file, toy_payloads, write_catalog

from rmgpy.kmc import kernel_rates as K
from rmgpy.kmc import kernel_reference as R
from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import run as kr
from rmgpy.kmc.melt_observables import WHOLE_ENSEMBLE_KEYS
from rmgpy.kmc.melt_state import materialize

FIXTURES = Path(__file__).parent / "fixtures"
T = 700.0
HOMOLYSIS = "homolysis:substrate=chain"
RECOMBINATION = "recombination:radical=radical_end:CH2+radical=radical_end:CH2"
BETA = "beta_scission:substrate=radical_interior:CH"


def toy_table(tmp_path, payloads=None):
    path = pe_catalog_file(tmp_path) if payloads is None else \
        write_catalog(tmp_path / "pe_ovr.json", payloads, tag="pe-ovr", bind=S.pe_mapping({1: 1}))
    doc = json.loads(Path(path).read_text())
    return K.kernel_rates(doc, S.pe_spec({10: 1}), T)


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    return toy_table(tmp_path_factory.mktemp("toy"))


@pytest.fixture(scope="module")
def pe():
    text = (FIXTURES / "i213_pe_catalog.json.zlib.b64").read_bytes().replace(b"\n", b"")
    doc = json.loads(zlib.decompress(base64.b64decode(text)))
    deck = kr.read_deck(FIXTURES / "i213_pe_toy_deck.py")
    return {"doc": doc, "deck": deck, "table": K.kernel_rates(doc, deck.spec, T)}


def flat(rows):
    first = {k: rows[0][k] for k in WHOLE_ENSEMBLE_KEYS}
    return all({k: r[k] for k in WHOLE_ENSEMBLE_KEYS} == first for r in rows)


def within(comparison, z=3.0):
    return [r for r in comparison if not abs(r["z"]) <= z]


# --------------------------------------------------------------------------
# null 1 and 2: nothing moves
# --------------------------------------------------------------------------
def test_recombination_alone_without_radicals_is_flat(toy):
    model = R.KernelModel(toy, S.pe_spec({10: 50}), kernels=[RECOMBINATION])
    times = [0.0, 1.0, 10.0]
    assert flat(R.population_balance(model, times))
    rows, events = R.stochastic(model, times, 0)
    assert events == 0 and flat(rows)


def test_phi_zero_is_flat_for_the_whole_kernel_set(toy):
    """Every radical-creating term is cage-thinned to nothing: the dead
    chains never react (to the bit, both references)."""
    model = R.KernelModel(toy, S.pe_spec({10: 50}), phi=0.0)
    assert len(model.bis) >= 3 and model.uni
    times = [0.0, 0.5, 5.0]
    pbe = R.population_balance(model, times)
    assert flat(pbe) and all(numpy.array_equal(r["counts"], pbe[0]["counts"]) for r in pbe)
    rows, events = R.stochastic(model, times, 3)
    assert events == 0 and flat(rows)


def test_phi_scales_the_radical_creating_rates(toy):
    full = R.KernelModel(toy, S.pe_spec({10: 50}), kernels=[HOMOLYSIS])
    half = R.KernelModel(toy, S.pe_spec({10: 50}), kernels=[HOMOLYSIS], phi=0.5)
    i = full.index("chain", 20)
    assert sum(r for r, _t in half.uni[i]) == pytest.approx(0.5 * sum(r for r, _t in full.uni[i]), rel=1e-15)


# --------------------------------------------------------------------------
# null 3: random scission, analytic
# --------------------------------------------------------------------------
def scission_table(spec, w):
    """A hand-built table: one chain kernel cutting every backbone bond into
    two chains at w per bond (no catalog behind it)."""
    piece = [["chain", 0, 1, 1], ["chain", 1, -1, -1]]
    law = {"class": "chain", "family": "chain", "anchor": "chain", "exact": {"1": []}, "exact_to_n": 1,
           "far": {"edge": 1, "bulk": [[1, piece]], "first": {"0": [[1, piece]]}, "last": {"0": []}}}
    e = spec.end_composition.molar_mass_g_mol()
    return {"schema": K.SCHEMA, "temperature_K": T, "catalog": {"chemistry_sha256": kr.chemistry_sha256(spec)},
            "moment_unit": {"periods_per_declared_unit": 2, "backbone_atoms_per_period": 1,
                            "molar_mass_g_mol": spec.repeat_unit.composition.molar_mass_g_mol() / 2,
                            "composition": {"C": "1", "H": "2"}},
            "families": {"chain": {"label_far": "chain", "min_n": 1, "labels": {},
                                   "molar_mass_offset_g_mol": e, "composition_offset": {"H": "2"},
                                   "declared_species": {}}},
            "kernels": {HOMOLYSIS: {"type": "homolysis", "molecularity": 1, "classes": ["chain"],
                                    "templates": [{"template": "scission", "phi_thinned": False,
                                                   "law": law, "overrides": []}]}},
            "kmc_only": [], "templates": [{"template": "scission", "weight": w}]}


def montroll_simha(L, n, kt, n0):
    """Chains of n bonds-plus-one periods from n0 chains of L periods cut at
    random, per-bond rate k, at time t (kt = k*t)."""
    p = 1.0 - math.exp(-kt)
    if n == L:
        return n0 * math.exp(-kt * (L - 1))
    return n0 * math.exp(-kt * (n - 1)) * (2.0 * p + (L - n - 1) * p * p)


def test_random_scission_is_montroll_simha_in_both_references():
    spec, w, n0 = S.pe_spec({10: 40}), 0.7, 40
    model = R.KernelModel(scission_table(spec, w), spec)
    times = [0.0, 0.05, 0.2, 1.0]
    pbe = R.population_balance(model, times)
    exact_rows = []
    for t, row in zip(times, pbe):
        exact = numpy.zeros(model.size)
        for n in range(1, 21):
            exact[model.index("chain", n)] = montroll_simha(20, n, w * t, n0)
        assert numpy.allclose(row["counts"], exact, rtol=1e-7, atol=1e-9 * n0)
        # mu0 of every chain (bonds broken = mu0 - mu0(0)) is mu1 - (mu1 - mu0(0)) e^{-kt}
        assert row["counts"].sum() == pytest.approx(20 * n0 - 19 * n0 * math.exp(-w * t), rel=1e-8)
        exact_rows.append(model.observables(exact))
    ens = R.ensemble(model, times, list(range(48)))
    fake_pbe = [dict(r, time_s=t) for r, t in zip(exact_rows, times)]
    assert not within(R.comparison(fake_pbe, ens))


# --------------------------------------------------------------------------
# null 4: PE homolysis alone, the exact linear master equation
# --------------------------------------------------------------------------
def engine_cuts(pe, units):
    """The fragment lengths (backbone atoms) of every homolysis the melt
    engine allows on the declared chain of ``units`` units."""
    from rmgpy.kmc.kernel_rates import _single_chain_spec
    spec = _single_chain_spec(pe["deck"].spec, units)
    engine = me.MeltEngine(spec, me.PayloadCatalog(pe["doc"]["payloads"]), T, numpy.random.default_rng(0))
    (ci,) = [i for i, t in enumerate(engine.templates) if t.name.startswith("5a88d9573d8b")]
    state = engine.state
    (cid,) = state.components.ids()
    allowed = engine._scan([cid])[cid].uni[ci][0]
    view = me._LocalView(state, [cid])
    cuts = []
    for tr in allowed:
        tx, cids = engine._begin_apply(tr, view)
        cuts.append(tuple(sorted(state.components.backbone_count(c) for c in cids)))
        tx.rollback()
    return sorted(cuts)


def test_pe_homolysis_cut_law_is_the_engines(pe):
    """The model's channels on the declared 60-period chain are the melt
    engine's own transitions on it: 53 cuts, fragments 4..56."""
    model = R.KernelModel(pe["table"], pe["deck"].spec, kernels=[HOMOLYSIS])
    w = next(r["weight"] for r in pe["table"]["templates"] if r["template"].startswith("5a88d9573d8b"))
    cuts = engine_cuts(pe, 30)
    assert len(cuts) == 53 and cuts[0] == (4, 56)
    channels = model.uni[model.index("chain", 60)]
    got = {tuple(sorted(model.state(t)[1] for t in targets)): rate / w for rate, targets in channels}
    want = {c: float(k) for c, k in Counter(cuts).items()}
    assert set(got) == set(want) and got == pytest.approx(want, rel=1e-12)


def test_pe_homolysis_alone_is_its_exact_master_equation(pe):
    model = R.KernelModel(pe["table"], pe["deck"].spec, kernels=[HOMOLYSIS])
    w = next(r["weight"] for r in pe["table"]["templates"] if r["template"].startswith("5a88d9573d8b"))
    times = [0.0, 1.0e7, 3.0e7, 6.0e7]
    pbe = R.population_balance(model, times)
    for t, row in zip(times, pbe):
        decay = math.exp(-53.0 * w * t)
        assert row["counts"][model.index("chain", 60)] == pytest.approx(100.0 * decay, rel=1e-7)
        for m in range(1, 61):
            want = 2.0 * 100.0 * (1.0 - decay) / 53.0 if 4 <= m <= 56 else 0.0
            assert row["counts"][model.index("radical_end", m)] == pytest.approx(want, rel=1e-7, abs=1e-9)
    ens = R.ensemble(model, times, list(range(32)))
    assert not within(R.comparison(pbe, ens))


# --------------------------------------------------------------------------
# V1 against V0 on a bilinear toy
# --------------------------------------------------------------------------
def test_the_stochastic_run_agrees_with_the_balance_on_the_toy(toy):
    model = R.KernelModel(toy, S.pe_spec({10: 400}), n_max=120)
    assert {b.template for b in model.bis} and any(b.join for b in model.bis)
    times = [0.0, 0.1, 0.25, 0.5]
    pbe = R.population_balance(model, times)
    ens = R.ensemble(model, times, list(range(16)))
    assert min(ens["events"]) > 500
    assert not within(R.comparison(pbe, ens))


def test_a_symmetric_recombination_counts_each_pair_once(toy):
    """A radicals of 6 periods: the balance makes chains of 12 at
    (w/2) A^2 / (N_A V), the stochastic run's propensity is (w/2) A (A-1) /
    (N_A V) (the melt engine's divisor and same-molecule term)."""
    model = R.KernelModel(toy, S.pe_spec({10: 50}), kernels=[RECOMBINATION])
    w = S.PE_WEIGHTS["recombination"]
    N = numpy.zeros(model.size)
    A = 10.0
    N[model.index("radical_end", 6)] = A
    bal = R._Balance(model)
    d = bal.rhs(0.0, numpy.concatenate([N, numpy.zeros(3 * len(model.families))]))
    assert d[model.index("chain", 12)] == pytest.approx(0.5 * w * A * A / model.nav, rel=1e-12)
    assert R.propensities(model, N)[1] == pytest.approx(0.5 * w * A * (A - 1) / model.nav, rel=1e-12)


def test_the_jacobian_is_the_derivative_of_the_balance(toy):
    model = R.KernelModel(toy, S.pe_spec({10: 30}), n_max=50)
    bal = R._Balance(model)
    rng = numpy.random.default_rng(5)
    y = numpy.concatenate([rng.random(model.size) * 3, numpy.zeros(3 * len(model.families))])
    J = bal.jac(0.0, y)
    h = 1e-6
    for col in rng.choice(model.size, 25, replace=False):
        e = numpy.zeros_like(y)
        e[col] = h
        fd = (bal.rhs(0.0, y + e) - bal.rhs(0.0, y - e)) / (2 * h)
        assert numpy.allclose(J[:, col], fd, rtol=1e-6, atol=1e-9 * numpy.abs(fd).max())


# --------------------------------------------------------------------------
# depository overrides at their exact sizes
# --------------------------------------------------------------------------
def _rows_json(smiles):
    rows = S.smiles_rows(smiles)
    return [[a, el, iso, ch, rad, lp, [[n, o] for n, o in edges]] for a, (el, iso, ch, rad, lp, edges)
            in sorted(rows.items())]


def override_payloads(spec):
    """The toy catalog with depository overrides: pentyl + heptyl ->
    n-dodecane at 10x the recombination weight, and 4-dodecyl -> 1-pentene +
    heptyl at 10x the beta-scission weight."""
    payloads = toy_payloads([t for t in S.pe_templates() if t.name in PE_LIVE])
    extra = {"recombination": (["[CH2]CCCC", "[CH2]CCCCCC"], ["CCCCCCCCCCCC"], 10 * S.PE_WEIGHTS["recombination"]),
             "beta_scission": (["CCC[CH]CCCCCCCC"], ["C=CCCC", "[CH2]CCCCCC"], 10 * S.PE_WEIGHTS["beta_scission"])}
    for p in payloads:
        for t in p["templates"]:
            if t["name"] in extra:
                r, q, w = extra[t["name"]]
                t["overrides"] = [{"reaction": {"reactants": _rows_json(r), "products": _rows_json(q)},
                                   "weight": {"T": [300.0, 1000.0], "w": [w, w]}}]
    return payloads


@pytest.fixture(scope="module")
def ovr(tmp_path_factory):
    return toy_table(tmp_path_factory.mktemp("ovr"), override_payloads(S.pe_spec({1: 1})))


def term(table, kernel):
    (t,) = table["kernels"][kernel]["templates"]
    return t


def test_the_table_carries_the_overrides_at_their_sizes(ovr):
    (rec,) = term(ovr, RECOMBINATION)["overrides"]
    assert rec["reactants"] == [["radical_end:CH2", 5, None], ["radical_end:CH2", 7, None]]
    assert rec["products"] == [["chain", 12]] and rec["multiplicity"] == 2
    assert rec["weight"] == pytest.approx(10 * S.PE_WEIGHTS["recombination"])
    (beta,) = term(ovr, BETA)["overrides"]
    assert beta["reactants"] == [["radical_interior:CH", 12, [3]]]
    assert beta["products"] == [["alkene_terminal:CH2=CH", 5], ["radical_end:CH2", 7]]
    assert beta["multiplicity"] == 1 and ovr["overrides_unplaced"] == []


def test_a_bimolecular_override_takes_its_weight_at_its_sizes_only(ovr):
    """Pentyl (A) and heptyl (B): dodecane forms at w_o A B / (N_A V); two
    heptyls make tetradecane at the template's (w/2) B^2 / (N_A V)."""
    model = R.KernelModel(ovr, S.pe_spec({10: 50}), kernels=[RECOMBINATION])
    w = S.PE_WEIGHTS["recombination"]
    N = numpy.zeros(model.size)
    A, B = 7.0, 11.0
    N[model.index("radical_end", 5)] = A
    N[model.index("radical_end", 7)] = B
    d = R._Balance(model).rhs(0.0, numpy.concatenate([N, numpy.zeros(3 * len(model.families))]))
    assert d[model.index("chain", 12)] == pytest.approx(10 * w * A * B / model.nav, rel=1e-12)
    assert d[model.index("chain", 14)] == pytest.approx(0.5 * w * B * B / model.nav, rel=1e-12)
    assert d[model.index("chain", 10)] == pytest.approx(0.5 * w * A * A / model.nav, rel=1e-12)


def test_a_unimolecular_override_takes_its_weight_at_its_size_only(ovr):
    """4-dodecyl is 2 of the 18 placed beta scissions of a 12-period
    interior radical (x = 3 and its mirror, one per outcome over 9 valid
    positions): 1-pentene + heptyl at w_o * 2/9; 1-hexene + hexyl, one
    period off, at w * 2/9."""
    model = R.KernelModel(ovr, S.pe_spec({10: 50}), kernels=[BETA])
    w = S.PE_WEIGHTS["beta_scission"]
    N = numpy.zeros(model.size)
    N[model.index("radical_interior", 12)] = 1.0
    d = R._Balance(model).rhs(0.0, numpy.concatenate([N, numpy.zeros(3 * len(model.families))]))
    assert d[model.index("alkene_terminal", 5)] == pytest.approx(10 * w * 2 / 9, rel=1e-12)
    assert d[model.index("alkene_terminal", 6)] == pytest.approx(w * 2 / 9, rel=1e-12)
    N[:] = 0.0
    N[model.index("radical_interior", 13)] = 1.0
    d = R._Balance(model).rhs(0.0, numpy.concatenate([N, numpy.zeros(3 * len(model.families))]))
    assert d[model.index("alkene_terminal", 5)] == pytest.approx(w * 2 / 10, rel=1e-12)


def test_the_stochastic_run_honours_an_override(ovr):
    """Chains of 12 periods cut into pentyl + heptyl (one homolysis channel
    of seven) and recombine: both references at the override's weight."""
    model = R.KernelModel(ovr, S.pe_spec({6: 300}), kernels=[HOMOLYSIS, RECOMBINATION], n_max=48)
    assert len(model.pair_overrides) == 1 and model.pair_overrides[0]["delta"] > 0
    times = [0.0, 0.2, 0.6]
    pbe = R.population_balance(model, times)
    ens = R.ensemble(model, times, list(range(16)))
    assert not within(R.comparison(pbe, ens))


# --------------------------------------------------------------------------
# what the references refuse, list and write
# --------------------------------------------------------------------------
def test_inert_and_unmodelled_terms_are_listed_not_dropped(pe):
    model = R.KernelModel(pe["table"], pe["deck"].spec)
    terms = sum(len(k["templates"]) for k in pe["table"]["kernels"].values())
    assert len(model.inert) + len(model.not_modelled) + len(model.uni_terms) + len(model.bis) == terms
    assert [x["template"][:12] for x in model.not_modelled] == ["55593e44e6ca"]
    assert model.unbalanced == []       # every measured law conserves each element


def test_a_table_of_another_chemistry_is_refused(toy):
    other = dict(toy, catalog={"chemistry_sha256": "0" * 64})
    with pytest.raises(R.KernelReferenceError, match="chemistry digest"):
        R.KernelModel(other, S.pe_spec({10: 5}))


def test_the_command_line_writes_both_time_series_beside_the_table(toy, tmp_path, capsys):
    path = tmp_path / "toy.kernel_rates.json"
    path.write_text(json.dumps(toy))
    deck = tmp_path / "deck.py"
    from meltRunTest import write_deck
    write_deck(deck, S.pe_mapping({10: 60}))
    code = R.main([str(deck), str(path), "--end", "0.2", "--snapshots", "3", "--seeds", "4", "--n-max", "80"])
    assert code == 0
    doc = json.loads((tmp_path / "toy.kernel_rates.reference.json").read_text())
    assert doc["schema"] == R.SCHEMA and len(doc["pbe"]) == 3 and len(doc["ssa"]["rows"]) == 3
    assert {r["quantity"] for r in doc["comparison"]} == set(WHOLE_ENSEMBLE_KEYS)
    assert "written" in capsys.readouterr().out


# --------------------------------------------------------------------------
# the law-domain diagnostic of a MoM kernel set
# --------------------------------------------------------------------------
def test_law_domain_fractions_count_each_familys_chains_below_each_validity_length(toy):
    """A MoM term's site law miscounts below its validity length (homolysis,
    min_fragment 4: a chain of 4 periods under sites n - 7 counts -3); the
    balance knows the distribution, so it reports the fraction of each
    participant family's chains there."""
    model = R.KernelModel(toy, S.pe_spec({10: 1}), n_max=20)
    counts = numpy.zeros(model.size)
    counts[model.index("chain", 4)] = 1.0
    counts[model.index("chain", 20)] = 1.0
    over = numpy.zeros(3 * len(model.families))
    validity = [{"kernel": "k", "template": "t", "families": ["chain", "radical_end"], "valid_from_periods": 7}]
    assert R.law_domain_fractions(model, counts, over, validity) == [
        {"kernel": "k", "template": "t", "family": "chain", "valid_from_periods": 7, "mu0_fraction_below": 0.5},
        {"kernel": "k", "template": "t", "family": "radical_end", "valid_from_periods": 7,
         "mu0_fraction_below": None}]
    over[3 * model.families.index("chain")] = 2.0          # chains past the grid are above any length
    (chain, _end) = R.law_domain_fractions(model, counts, over, validity)
    assert chain["mu0_fraction_below"] == 0.25
    validity[0]["valid_from_periods"] = 4                  # strictly shorter only
    assert R.law_domain_fractions(model, counts, over, validity)[0]["mu0_fraction_below"] == 0.0


def test_the_validity_lengths_are_read_from_a_producer_manifest(tmp_path):
    entry = {"kernel": "k", "template": "t", "families": ["chain"], "valid_from_periods": 5}
    path = tmp_path / "producer.json"
    path.write_text(json.dumps({"facts": {"kernels": {"record": {"validity": [entry]}}}}))
    assert R.read_validity(path) == [entry]
    path.write_text(json.dumps([entry]))
    assert R.read_validity(path) == [entry]
    path.write_text(json.dumps({"facts": {}}))
    with pytest.raises(R.KernelReferenceError, match="validity"):
        R.read_validity(path)
