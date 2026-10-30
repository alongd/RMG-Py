#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""I-286: the converter books what kernel-mode pools hold, and a depropagation's residue is declared.

In kernel mode a pool's chain weighs M_u n + e_p, e_p its kernel-rates family's
offset, and every kernel term conserves sum_p (M_u mu1_p + e_p mu0_p): a cut
that is not on a unit boundary moves mass from units into end groups.  The
booking kept each pool's end groups at their t = 0 count (the H_p = 0 ledger
of the legacy kernels, which mint chain ends from nothing), so a born-empty
pool with e_p != 0 booked none, and the PP full table lost 2.5 g of 1000 by
1e4 s.  A pool whose end mass is its family's books e_p mu0_p(t) and no
ledger.

MoM's k_depropagation drains the chains its last unzip leaves with no unit;
what is left is the pool's end groups.  The producer names that residue for
the engine (``depropagation_residue``): the pool's cap whose mass is e_p (the
cap away from the radical end), and the converter classes it by kMC's
ancestry of that cap (its backbone atoms, B < b: light).
"""

import json

import pytest

from momKernelProducerTest import RADICAL, produce_kernels, table
from momProducerTest import ADJ, LABEL, make_setup, mom_entries, write_run
from momReleasedMonomerTest import _pp

from rmgpy.kmc import mom_producer as mp
from rmgpy.molecule import Molecule

H_MASS = Molecule(smiles="[H]").get_molecular_weight() * 1000.0


@pytest.fixture
def setup(tmp_path):
    return make_setup(tmp_path)


def test_kernel_mode_books_each_pools_end_groups_at_its_current_chain_count(setup):
    """Chains cut into end radicals (2 x one H cap from two H caps), then
    radicals drained at DP 0 whose H caps go to H: the system's mass is
    constant, the radical pool's end groups fall below their t = 0 count (0),
    and the mass closes only when each pool books e_p mu0_p(t).  The whole
    mu2 and Mw stand: every chain still weighs a n + e_p."""
    out = setup["tmp"] / "out"
    produce_kernels(setup, out)
    mu = [float(x) for x in setup["spec"].macroscopic_moments_mol]
    d, x = 1e-3, 4e-4

    def chains(t):
        return mu[0] - d * t, mu[1] - 3 * d * t, mu[2] - 40 * d * t

    def radicals(t):
        return 2 * d * t - x * t, 3 * d * t, 40 * d * t
    write_run(out, times=[0.0, 10.0, 20.0, 40.0, 45.0], pools=[(LABEL, chains, 0.0), (RADICAL, radicals, 0.0)],
              gas=[("H(1)", ADJ["H"], lambda t: x * t), ("H2(2)", ADJ["H2"], lambda t: 0.0)], rows=[])
    manifest = produce_kernels(setup, out, mode="convert")
    offsets = {f: v["molar_mass_offset_g_mol"] for f, v in table()["families"].items()}
    a, e = manifest["facts"]["a_g_mol"], manifest["facts"]["e_g_mol"]
    s = (offsets["radical_end"] - e) / a
    diagnostics = json.loads((out / "mom_diagnostics.json").read_text())
    for entry, diag in list(zip(mom_entries(out), diagnostics))[1:]:
        q, t = entry["quantities"], entry["time_s"]
        c, r = chains(t), radicals(t)
        assert diag["mass_closure_residual_g"] == pytest.approx(0.0, abs=1e-9)
        assert q["polymer_condensed_mass_g"] == pytest.approx(
            a * (c[1] + r[1]) + offsets["chain"] * c[0] + offsets["radical_end"] * r[0], rel=1e-12)
        assert q["whole_mu2_mol"] == pytest.approx(c[2] + r[2] + 2 * s * r[1] + s * s * r[0], rel=1e-9)
        assert "whole_Mw_g_mol" in q


def test_the_radical_pools_depropagation_names_its_residue(tmp_path):
    """The PE toy's end beta scission: the radical pool (one H cap) leaves
    an H atom when its last unit unzips."""
    from endClassPoolsTest import FIXTURES, K, carried, kr, spec
    from momKernelProducerTest import polymers
    from momProducerTest import CHECKOUT
    payloads = []
    for name in ("i261_pe_toy_recombination_payloads.json", "i262_pe_toy_addition_payloads.json"):
        payloads += json.loads((FIXTURES / name).read_text())["payloads"]
    doc = K.kernel_rates({"schema": "kmc-template-catalog/1", "payloads": payloads,
                          "provenance": {"chemistry_sha256": kr.chemistry_sha256(spec())}}, spec(), 700.0)
    setup = make_setup(tmp_path)
    path = tmp_path / "pe.kernel_rates.json"
    path.write_text(json.dumps(doc))
    out = tmp_path / "out"
    manifest = mp.produce(setup["deck"], checkout=CHECKOUT, temperature=700.0, end=40.0, snapshots=9, families=(),
                          catalogs=(), out=out, database=setup["db"], kernel_rates=str(path), kernels=carried(doc))
    pools = polymers(out)
    assert pools[RADICAL]["depropagation_residue"] == "[H]"
    assert manifest["facts"]["kernels"]["record"]["depropagation_residue"] == {
        RADICAL: {"smiles": "[H]", "backbone_atoms": 0, "molar_mass_g_mol": pytest.approx(H_MASS, rel=1e-12)}}
    assert all("depropagation_residue" not in p for label, p in pools.items() if label != RADICAL)


def test_a_chain_types_caps_carry_their_backbone_atoms(tmp_path):
    """PP (CH2-CH(CH3), b = 2 backbone atoms per unit): the c2~c2 type
    starts one atom into the period, so its head cap is c2 with its methyl
    (ethyl, one backbone atom); the c1~c1 type's tail cap is c1 (methyl, one
    backbone atom); the deck's H caps carry none."""
    spec = make_setup(tmp_path, edit=_pp)["spec"]

    def caps(s, r):
        return [(m.to_smiles(), b) for m, b in mp._type_caps(spec, {"start_offset": s, "residue": r}, 2, 1)]
    assert caps(1, 1) == [("C[CH2]", 1), ("[H]", 0)]
    assert caps(0, 1) == [("[H]", 0), ("[CH3]", 1)]
    assert caps(0, 0) == [("[H]", 0), ("[H]", 0)]


@pytest.mark.parametrize("caps, e, want", [
    ([("C[CH2]", 1), ("[H]", 0)], 29.061135952360928, ("C[CH2]", 1)),
    ([("[H]", 0), ("[H]", 0)], H_MASS, ("[H]", 0)),
    ([("[H]", 0), ("[CH3]", 1)], 15.034553478471935, ("[CH3]", 1)),
])
def test_the_residue_is_the_cap_that_weighs_the_pools_end_mass(caps, e, want):
    got = mp._depropagation_residue("P", [(Molecule(smiles=s), b) for s, b in caps], e)
    assert (got[0].to_smiles(), got[1]) == want


@pytest.mark.parametrize("caps, e", [([("[H]", 0), ("[H]", 0)], 2.0159), ([("C[CH2]", 1), ("C[CH2]", 0)], 29.0611)])
def test_a_residue_no_cap_or_two_ancestries_weigh_is_refused(caps, e):
    with pytest.raises(mp.ProducerError, match="residue"):
        mp._depropagation_residue("P", [(Molecule(smiles=s), b) for s, b in caps], e)


def test_a_declared_residue_takes_its_caps_ancestry(tmp_path):
    """Ethyl from a PP c2c2 radical: the element bracket cannot class it (C
    is backbone and methyl), the residue rule gives kMC's class at B = 1 <
    b = 2: light.  Without the declaration, or for another molecule, it stays
    undetermined."""
    spec = make_setup(tmp_path, edit=_pp)["spec"]
    ethyl, propane = Molecule(smiles="C[CH2]"), Molecule(smiles="CCC")
    plain = mp.KmcClasses(spec, 4)
    declared = mp.KmcClasses(spec, 4, residues=[(Molecule(smiles="C[CH2]"), 1)])
    assert plain.classify_with_rule(ethyl) == ("undetermined", "undetermined")
    assert declared.classify_with_rule(ethyl) == ("light", "end_residue")
    assert declared.classify_with_rule(propane) == ("undetermined", "undetermined")
    assert declared.classify_with_rule(Molecule(smiles="[H]")) == ("species", "element")   # the deck declares H
