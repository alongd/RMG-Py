#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""I-284: the converter classes a released monomer the way kMC classes it.

A depropagation releases the declared repeat unit with its head-tail bond
closed (propene from PP, styrene from PS).  kMC gives that molecule the
ancestry of the unit it was cut from: the atoms on the unit's head-to-tail
path are backbone (``melt_state.monomer_backbone_keys``' map), so B_c = b and
``MeltObserver.category`` puts it in the polymer class at DP 1.  Where carbon
is both backbone and side group, the element bracket alone cannot decide, and
the converter withheld every whole-ensemble quantity while the monomer was
present.  A species that is not the closed unit (toluene from PS, propane
from PP) stays a named refusal.
"""

import json

import pytest

import melt_engine_support as S
from meltRunTest import pe_catalog_file, write_deck
from meltStateTest import ps_mapping
from momProducerTest import LABEL, PROFILE_T, WHOLE, ZERO_GAS, declared, make_setup, mom_entries, mw, produce, \
    write_run

from rmgpy.kmc import mom_producer as mp
from rmgpy.molecule import Molecule
from rmgpy.polymer_spec import LinearHomopolymerSpec


def _pp(mapping):
    """CH2-CH(CH3), H caps: the methyl carbon is a side group, c1 and c2 the backbone."""
    mapping["label"] = "PP"
    mapping["repeat_unit"] = {
        "atoms": [S._atom("c1", hydrogens=2), S._atom("c2", hydrogens=1), S._atom("c3", hydrogens=3)],
        "bonds": [{"atoms": ["c1", "c2"], "order": 1}, {"atoms": ["c2", "c3"], "order": 1}],
        "ports": {"head": {"atom": "c1", "bond_order": 1}, "tail": {"atom": "c2", "bond_order": 1}},
        "electronic_state": {"multiplicity": 1},
    }


def _ps_setup(tmp_path):
    mapping = ps_mapping({10: 1, 20: 1}, mobile={})
    mapping["small_species"] = {}
    mapping["transport_model"] = {"met": S.ps_met()}
    deck = write_deck(tmp_path / "ps_deck.py", mapping)
    (tmp_path / "db" / "input" / "kinetics" / "families").mkdir(parents=True, exist_ok=True)
    return {"deck": deck, "catalog": pe_catalog_file(tmp_path), "db": str(tmp_path / "db"), "tmp": tmp_path,
            "spec": LinearHomopolymerSpec.from_mapping(mapping)}


def _setup(tmp_path, deck):
    return {"pe": lambda: make_setup(tmp_path), "pp": lambda: make_setup(tmp_path, edit=_pp),
            "ps": lambda: _ps_setup(tmp_path)}[deck]()


def _release(s, out, label, smiles, x):
    """A run whose one pool loses x t mol of ``smiles`` (column M(4)) to the gas: the mass closes."""
    spec = s["spec"]
    a = spec.repeat_unit.composition.molar_mass_g_mol()
    mu = declared(spec)
    m = mw(smiles)
    produce(s, out)
    write_run(out, times=PROFILE_T, label=label,
              pools=[(label, lambda t: (mu[0], mu[1] - m * x * t / a, mu[2]), 0.0)], mass_override=a,
              gas=ZERO_GAS + [("M(4)", Molecule(smiles=smiles).to_adjacency_list(), lambda t: x * t)])
    return produce(s, out, mode="convert")


@pytest.mark.parametrize("deck, label, smiles, rule", [("pp", "PP", "C=CC", "closed_unit"),
                                                       ("ps", "PS", "C=Cc1ccccc1", "closed_unit"),
                                                       ("pe", LABEL, "C=C", "element")])
def test_a_released_monomer_is_polymer_class_and_the_whole_rows_stand(tmp_path, deck, label, smiles, rule):
    """The monomer is a polymer-derived molecule at its mass DP (M - e) / a,
    counted beside the pool's chains; PE's ethylene was already decided by
    element (every carbon is backbone) and keeps that rule."""
    s = _setup(tmp_path, deck)
    out = tmp_path / "out"
    x = 1e-4
    manifest = _release(s, out, label, smiles, x)
    assert manifest["read"]["species_classes"]["M(4)"] == "polymer"
    assert manifest["read"]["species_class_rules"]["M(4)"] == rule
    assert manifest["read"]["species_class_rules"]["H(1)"] == "element"
    spec = s["spec"]
    a, e = spec.repeat_unit.composition.molar_mass_g_mol(), spec.end_composition.molar_mass_g_mol()
    mu = declared(spec)
    n_m = (mw(smiles) - e) / a
    entries = mom_entries(out)
    diagnostics = json.loads((out / "mom_diagnostics.json").read_text())
    assert len(entries) == 9
    for entry, diag in zip(entries, diagnostics):
        t, q = entry["time_s"], entry["quantities"]
        assert diag["unclassified_species"] == []
        assert diag["polymer_withheld"] is None
        assert set(WHOLE) <= set(q)
        assert q["whole_mu0_mol"] == pytest.approx(mu[0] + x * t, rel=1e-12)
        assert q["polymer_derived_molecules_mol"] == pytest.approx(mu[0] + x * t, rel=1e-12)
        # pool chains at mass DP (pool mu1 a + e mu0 - e mu0) / a = pool mu1, plus the monomer's n_m
        assert q["whole_mu1_mol"] == pytest.approx(mu[1] - mw(smiles) * x * t / a + n_m * x * t, rel=1e-12)
        assert diag["mass_closure_residual_g"] == pytest.approx(0.0, abs=1e-9)
        assert q["polymer_condensed_mass_g"] == pytest.approx(spec.initial_mass_kg * 1000.0, rel=1e-12)


@pytest.mark.parametrize("deck, label, smiles", [("pp", "PP", "CCC"), ("ps", "PS", "Cc1ccccc1")])
def test_a_species_that_is_not_the_closed_unit_stays_a_named_refusal(tmp_path, deck, label, smiles):
    """Propane from PP and toluene from PS: they are not the declared unit
    closed head-tail, so kMC's class depends on ancestry MoM does not keep."""
    s = _setup(tmp_path, deck)
    out = tmp_path / "out"
    manifest = _release(s, out, label, smiles, 1e-4)
    assert manifest["read"]["species_classes"]["M(4)"] == "undetermined"
    assert manifest["read"]["species_class_rules"]["M(4)"] == "undetermined"
    entries = mom_entries(out)
    diagnostics = json.loads((out / "mom_diagnostics.json").read_text())
    for entry, diag in list(zip(entries, diagnostics))[1:]:
        assert not set(WHOLE) & set(entry["quantities"])
        assert diag["unclassified_species"] == ["M(4)"]
        assert "unclassifiable species present: ['M(4)']" in diag["polymer_withheld"]


def test_the_class_rules_are_named():
    assert set(mp.SPECIES_CLASS_RULES) == {"element", "closed_unit", "undetermined"}
    assert all(isinstance(v, str) and v for v in mp.SPECIES_CLASS_RULES.values())
