"""The MoM producer's kernel mode (slice S5 a + b of the shared-chemistry
design): a kernel-only MoM input written from a kernel-rates table.

The table fixture is the toy PE catalog's (``meltRunTest.pe_catalog_file``)
at 700 K from ``rmgpy.kmc.kernel_rates`` at kmc-i188 db7660c23
(``kmc-kernel-rates/3``; homolysis a = 1, b = -5 in CH2 periods, fragments
3..n - 3; end recombination per feature 1, zero at 2 periods).  Families are
off; the homolysis and end-recombination kernels carry the table's
coefficients, per template, onto a dead-chain pool and an explicit
end-radical pool.
"""

import ast
import json
from pathlib import Path

import pytest

from momProducerTest import (BASE, CHECKOUT, LABEL, ZERO_GAS, constant_run, make_setup, mom_entries, spec_of,
                             write_run)

from rmgpy.kmc import mom_producer as mp
from rmgpy.kmc import run as kr

TABLE = Path(__file__).parent / "fixtures" / "i219_toy_kernel_rates_700K.json"
H = "homolysis:substrate=chain"
C = "recombination:radical=radical_end:CH2+radical=radical_end:CH2"
RADICAL = LABEL + "_re"


def table():
    return json.loads(TABLE.read_text())


def row(doc, kernel):
    (name,) = [t["template"] for t in doc["kernels"][kernel]["templates"]]
    return next(r for r in doc["templates"] if r["template"] == name)


def terms(doc=None, kernels=(H, C), phi=1.0, counts=BASE, ratio=1.0):
    return mp.kernel_input(table() if doc is None else doc, spec_of(counts), list(kernels), temperature=700.0,
                           phi=phi, label=LABEL, radical_label=RADICAL, volume_ratio=ratio)


def produce_kernels(setup, out, **kw):
    args = dict(checkout=CHECKOUT, temperature=700.0, end=40.0, snapshots=9, families=(), catalogs=(),
                out=out, database=setup["db"], kernel_rates=str(TABLE), kernels=[H, C])
    args.update(kw)
    return mp.produce(setup["deck"], **args)


def polymers(out):
    tree = ast.parse((out / "rmg" / "input.py").read_text())
    calls = [s.value for s in tree.body if isinstance(s, ast.Expr) and isinstance(s.value, ast.Call)
             and getattr(s.value.func, "id", None) == "polymer"]
    out = {}
    for call in calls:
        kwargs = {k.arg: ast.literal_eval(k.value) for k in call.keywords}
        out[kwargs["label"]] = kwargs
    return out


@pytest.fixture
def setup(tmp_path):
    return make_setup(tmp_path)


# --------------------------------------------------------------------------
# the terms
# --------------------------------------------------------------------------
def test_the_terms_are_the_tables_coefficients():
    doc = table()
    out = terms(doc)
    h, c = row(doc, H), row(doc, C)
    assert out["homolysis"] == [{"template": h["template"], "w": h["weight"], "phi": 1.0, "a": 1.0, "b": -5.0,
                                 "min_fragment": 3, "periods_per_unit": 2, "daughter": RADICAL}]
    assert out["recombination"] == [{"template": c["template"], "k_pair": c["weight"] / c["divisor"],
                                     "sites_per_chain": 1.0, "periods_per_unit": 2, "product": LABEL}]


def test_phi_thins_the_radical_creating_templates_only():
    out = terms(phi=0.25)
    assert out["homolysis"][0]["phi"] == 0.25
    assert out["recombination"][0]["k_pair"] == terms()["recombination"][0]["k_pair"]


def test_the_pair_coefficient_is_converted_to_moms_volume():
    """MoM's concentration is n / V_poly, kMC's n / V: k_pair = (w / divisor) V_poly / V."""
    c = row(table(), C)
    out = terms(ratio=0.97)
    assert out["recombination"][0]["k_pair"] == pytest.approx(c["weight"] / c["divisor"] * 0.97, rel=1e-15)
    assert out["record"]["volume_ratio_mom_over_kmc"] == 0.97


def test_a_single_kernel_is_written_alone():
    assert terms(kernels=[C])["homolysis"] == [] and len(terms(kernels=[C])["recombination"]) == 1
    assert terms(kernels=[H])["recombination"] == [] and len(terms(kernels=[H])["homolysis"]) == 1


# --------------------------------------------------------------------------
# refusals
# --------------------------------------------------------------------------
def test_a_kernel_moms_kernels_do_not_carry_is_refused():
    with pytest.raises(mp.ProducerError, match="no MoM kernel"):
        terms(kernels=[H, C, "beta_scission:substrate=radical_interior:CH"])


def test_an_unknown_kernel_is_refused():
    with pytest.raises(mp.ProducerError, match="not in the table"):
        terms(kernels=["homolysis:substrate=nothing"])


def test_a_table_of_another_chemistry_or_temperature_or_schema_is_refused():
    doc = table()
    doc["catalog"]["chemistry_sha256"] = "0" * 64
    with pytest.raises(mp.ProducerError, match="chemistry digest"):
        terms(doc)
    with pytest.raises(mp.ProducerError, match="temperature"):
        mp.kernel_input(table(), spec_of(BASE), [H, C], temperature=701.0, phi=1.0, label=LABEL,
                        radical_label=RADICAL, volume_ratio=1.0)
    doc = table()
    doc["schema"] = "kmc-kernel-rates/1"
    with pytest.raises(mp.ProducerError, match="kmc-kernel-rates"):
        terms(doc)


def test_a_homolysis_law_that_is_not_the_uniform_cut_is_refused():
    doc = table()
    law = next(t for t in doc["kernels"][H]["templates"])["law"]
    outcomes = law["exact"]["12"]
    outcomes[0][0] += 1                     # one position counted twice
    with pytest.raises(mp.ProducerError, match="not the uniform cut"):
        terms(doc)


def test_sites_off_the_cut_positions_are_refused():
    doc = table()
    for r in doc["templates"]:
        if r["template"] == row(doc, H)["template"]:
            r["participants"][0]["chain"]["b"] = -4.0
    with pytest.raises(mp.ProducerError, match="b = a"):
        terms(doc)


def test_a_recombination_profile_on_a_fragment_length_is_refused():
    """Zero sites on a radical length homolysis makes: per_feature is then not
    every reachable radical's count."""
    doc = table()
    for r in doc["templates"]:
        if r["template"] == row(doc, C)["template"]:
            for p in r["participants"]:
                p["radical_end:CH2"]["profile"]["fragment_backbone_atoms=4"] = [0]
    with pytest.raises(mp.ProducerError, match="profile"):
        terms(doc)


def test_initial_chains_shorter_than_the_sites_allow_are_refused():
    """a*n + b < 0 at n = 4 periods (2 units): the affine sites are not the count."""
    with pytest.raises(mp.ProducerError, match="shorter than"):
        terms(counts={2: 1, 20: 1})


# --------------------------------------------------------------------------
# the input and its run
# --------------------------------------------------------------------------
def test_kernel_mode_writes_a_kernel_only_input(setup):
    out = setup["tmp"] / "out"
    manifest = produce_kernels(setup, out)
    text = (out / "rmg" / "input.py").read_text()
    assert "kineticsFamilies='none'" in text
    pools = polymers(out)
    assert set(pools) == {LABEL, RADICAL}
    record = manifest["facts"]["kernels"]
    assert pools[LABEL]["kernel_homolysis"] == record["homolysis"]
    assert pools[RADICAL]["kernel_recombination"] == record["recombination"]
    assert pools[RADICAL]["end_radical_site"] == "primary"
    assert pools[RADICAL]["moments"] == [0.0, 0.0, 0.0] and pools[RADICAL]["initial_mass"] == 0.0
    assert pools[RADICAL]["end_groups"] == pools[LABEL]["end_groups"]
    assert record["record"]["volume_ratio_mom_over_kmc"] == manifest["facts"]["volume_ratio_mom_over_kmc"]
    assert record["record"]["table_sha256"] == mp._sha256(TABLE.read_bytes())
    names = {c["name"] for c in json.loads((out / "conventions.json").read_text())["conventions"]}
    assert {"kernels", "kernel_phi", "kernel_volume"} <= names


def test_kernel_mode_writes_no_families(setup):
    with pytest.raises(mp.ProducerError, match="kernel mode writes no families"):
        produce_kernels(setup, setup["tmp"] / "out", families=["toy_homolysis"])


def test_kernel_mode_converts_a_moving_run_without_a_family_census(setup):
    out = setup["tmp"] / "out"
    produce_kernels(setup, out)
    mu = [float(x) for x in setup["spec"].macroscopic_moments_mol]
    d = 1e-3

    def chains(t):
        return mu[0] - d * t, mu[1] - 3 * d * t, mu[2] - 40 * d * t

    def radicals(t):
        return 2 * d * t, 3 * d * t, 40 * d * t
    write_run(out, times=[0.0, 10.0, 20.0, 40.0, 45.0], pools=[(LABEL, chains, 0.0), (RADICAL, radicals, 0.0)],
              gas=ZERO_GAS, rows=[])
    manifest = produce_kernels(setup, out, mode="convert")
    assert manifest["status"] == "converted"
    assert "inactive_chemistry" not in manifest["flags"]
    assert manifest["census"]["mode"] == "kernels" and manifest["census"]["profile_changes"]
    assert len(mom_entries(out)) == len(kr.observation_times(end=40.0, snapshots=9))


def test_kernel_mode_refuses_a_flat_run_as_inactive(setup):
    out = setup["tmp"] / "out"
    produce_kernels(setup, out)
    constant_run(out, setup["spec"], flat=True, rows=[])
    with pytest.raises(mp.ProducerError, match="inactive"):
        produce_kernels(setup, out, mode="convert")
