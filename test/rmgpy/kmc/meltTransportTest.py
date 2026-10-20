"""Declared melt transport with its sources, and the well-mixed bound beside it.

``transport_model["sources"]`` (beside ``met``) names a citation for every MET
field; it is outside ``chemistry_sha256`` like the rest of ``transport_model``.
``run(..., also_well_mixed=True)`` runs the same seeds a second time with the
segment diffusivity at ``run.WELL_MIXED_D0_M2_S`` (every escape factor exactly
1, the encounter limit off) and writes both passes, labelled, under one
``run.json``."""

import json
from pathlib import Path

import pytest

import melt_engine_support as S
from meltRunTest import lines, pe_catalog_file, pe_deck, write_deck

from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import run as kr
from rmgpy.polymer_spec import LinearHomopolymerSpec

DECKS = Path(__file__).parent / "fixtures" / "decks"
CITED_DECKS = ("pe_toy_deck.py", "pp_deck.py", "ps_deck.py")
MET_FIELDS = ("segment_diffusivity", "chain_scaling", "mean_square_end_to_end_per_molar_mass_m2_mol_per_g",
              "contact_radius_m", "spin_factor")


def cited(mapping, sources=None):
    out = json.loads(json.dumps(mapping))
    out["transport_model"]["sources"] = sources if sources is not None else {
        field: [f"a citation for {field}"] for field in MET_FIELDS}
    return out


# --------------------------------------------------------------------------
# sources
# --------------------------------------------------------------------------
def test_sources_stay_outside_the_chemistry_digest():
    """The catalogs bind on ``chemistry_sha256``: citing (or re-citing) the
    transport must not unbind them; the declaration digest does move."""
    bare = LinearHomopolymerSpec.from_mapping(S.pe_mapping({10: 1}))
    one = LinearHomopolymerSpec.from_mapping(cited(S.pe_mapping({10: 1})))
    other = LinearHomopolymerSpec.from_mapping(cited(S.pe_mapping({10: 1}), {
        field: ["another citation"] for field in MET_FIELDS}))
    assert kr.chemistry_sha256(bare) == kr.chemistry_sha256(one) == kr.chemistry_sha256(other)
    assert len({kr.declaration_sha256(s) for s in (bare, one, other)}) == 3


def test_sources_must_cite_every_met_field_and_nothing_else(tmp_path):
    missing = cited(S.pe_mapping({10: 1}))
    del missing["transport_model"]["sources"]["spin_factor"]
    with pytest.raises(kr.RunError, match="spin_factor"):
        kr.read_deck(write_deck(tmp_path / "missing.py", missing))
    unknown = cited(S.pe_mapping({10: 1}))
    unknown["transport_model"]["sources"]["viscosity"] = ["somebody"]
    with pytest.raises(kr.RunError, match="viscosity"):
        kr.read_deck(write_deck(tmp_path / "unknown.py", unknown))
    empty = cited(S.pe_mapping({10: 1}), {field: [] for field in MET_FIELDS})
    with pytest.raises(kr.RunError, match="citation"):
        kr.read_deck(write_deck(tmp_path / "empty.py", empty))
    assert kr.transport_sources(kr.read_deck(write_deck(tmp_path / "ok.py", cited(S.pe_mapping({10: 1})))).spec) \
        == {field: [f"a citation for {field}"] for field in MET_FIELDS}
    # a deck without sources still runs (the synthetic test decks)
    assert kr.transport_sources(kr.read_deck(pe_deck(tmp_path)).spec) is None


@pytest.mark.parametrize("name", CITED_DECKS)
def test_each_polymer_deck_declares_cited_transport(name):
    """PE, PP and PS declare their MET block with a citation per field, and the
    law is finite and positive over 600-800 K."""
    deck = kr.read_deck(DECKS / name)
    sources = kr.transport_sources(deck.spec)
    assert sources is not None and set(sources) == set(MET_FIELDS)
    assert all(sources[f] and all(isinstance(c, str) and len(c) > 20 for c in sources[f]) for f in MET_FIELDS)
    law = me.TransportLaw.from_spec(deck.spec)
    law.check_range(600.0, 800.0)
    # not the placeholder that froze the PE runs (cage escape ~1e-9)
    assert law.d0(700.0) > 1.0e-12


# --------------------------------------------------------------------------
# the well-mixed bound
# --------------------------------------------------------------------------
def test_the_well_mixed_spec_changes_only_the_segment_diffusivity():
    spec = LinearHomopolymerSpec.from_mapping(cited(S.pe_mapping({10: 1})))
    bound = kr.well_mixed_spec(spec, 700.0)
    assert kr.chemistry_sha256(bound) == kr.chemistry_sha256(spec)
    met, declared = bound.transport_model["met"], spec.transport_model["met"]
    assert met["segment_diffusivity"] == {"law": "arrhenius", "D_ref_m2_s": kr.WELL_MIXED_D0_M2_S,
                                          "T_ref_K": 700.0, "Ea_J_mol": 0.0}
    assert met["chain_scaling"] == {"law": "rouse"}
    assert {k: v for k, v in met.items() if k not in ("segment_diffusivity", "chain_scaling")} == \
        {k: v for k, v in declared.items() if k not in ("segment_diffusivity", "chain_scaling")}
    assert "sources" not in bound.transport_model
    # every escape factor is exactly one: the encounter limit is off
    law = me.TransportLaw.from_spec(bound)
    for n1, n2 in ((1, 1), (1, 10 ** 6), (10 ** 6, 10 ** 6)):
        assert me.escape_factor(1.0e30, law.k_diffusion(700.0, n1, n2)) == 1.0


def test_also_well_mixed_writes_both_passes_labelled(tmp_path):
    """The declared deck carries a slow transport (cage escape ~0); the bound
    pass runs the same seeds at phi = 1 into ``well_mixed/`` and both are in
    one run.json, each observation labelled with its pass."""
    deck = write_deck(tmp_path / "pe.py", cited(S.pe_mapping({10: 2, 20: 2})))
    catalog = pe_catalog_file(tmp_path)
    out = tmp_path / "out"
    manifest = kr.run(deck, [catalog], 700.0, [0.0, 1.0, 3.0], [1, 2], out, also_well_mixed=True)
    on_disk = json.loads((out / "run.json").read_text())
    assert on_disk["status"] == manifest["status"] == "complete"
    assert on_disk["transport"]["passes"] == ["declared", "well_mixed_bound"]
    assert on_disk["transport"]["declared"]["met"] == json.loads(
        kr.read_deck(deck).spec.normalized_json)["transport_model"]["met"]
    assert on_disk["transport"]["declared"]["sources"] == kr.transport_sources(kr.read_deck(deck).spec)
    bound_met = on_disk["transport"]["well_mixed_bound"]["met"]
    assert bound_met["segment_diffusivity"]["D_ref_m2_s"] == kr.WELL_MIXED_D0_M2_S
    wm = on_disk["well_mixed"]
    assert wm["status"] == "complete" and wm["comparison"] == "well_mixed/comparison.csv"
    assert [r["seed"] for r in wm["runs"]] == [r["seed"] for r in on_disk["runs"]] == [1, 2]
    declared_escape, bound_escape = [], []
    for r in on_disk["runs"]:
        rows = lines(out / r["path"])
        assert rows and all(x["transport"] == "declared" for x in rows)
        declared_escape.append(r["report"]["thinning"])
    for r in wm["runs"]:
        assert r["path"] == f"well_mixed/seed-{r['seed']}/observables.jsonl"
        rows = lines(out / r["path"])
        assert rows and all(x["transport"] == "well_mixed_bound" for x in rows)
        bound_escape.append(r["report"]["thinning"])
    assert (out / "comparison.csv").exists() and (out / "well_mixed" / "comparison.csv").exists()
    assert wm["leak"]["pooled"] == me.pooled_leak_fraction(wm["leak"]["per_seed"])
    # the bound pass accepts every thinned proposal; the declared one does not
    assert sum(t["proposals"] for t in bound_escape) > 0
    for t in bound_escape:
        assert t["accepted"] == t["proposals"]
        assert all(c["escape_probability"] in (None, 1.0) for c in t["classes"])
    assert sum(t["proposals"] for t in declared_escape) > 0
    escapes = [t["escape_probability"] for t in declared_escape if t["proposals"]]
    assert escapes and all(e < 1.0 for e in escapes), escapes


def test_without_the_flag_the_run_is_unchanged(tmp_path):
    deck, catalog = pe_deck(tmp_path), pe_catalog_file(tmp_path)
    manifest = kr.run(deck, [catalog], 700.0, [0.0, 1.0], [1], tmp_path / "plain")
    assert "well_mixed" not in manifest and "transport" not in manifest
    assert not (tmp_path / "plain" / "well_mixed").exists()
    assert all("transport" not in x for x in lines(tmp_path / "plain" / "seed-1" / "observables.jsonl"))


def test_cli_also_well_mixed(tmp_path):
    deck = write_deck(tmp_path / "pe.py", cited(S.pe_mapping({10: 1, 20: 1})))
    catalog = pe_catalog_file(tmp_path)
    out = tmp_path / "cli"
    argv = [str(deck), str(catalog), "--temperature", "700", "--end", "2", "--snapshots", "2",
            "--seeds", "3", "--out", str(out), "--also-well-mixed"]
    assert kr.main(argv) == 0
    manifest = json.loads((out / "run.json").read_text())
    assert manifest["well_mixed"]["status"] == "complete"
    assert (out / "well_mixed" / "seed-3" / "observables.jsonl").exists()
