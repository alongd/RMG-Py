"""Bookkeeping guards of the catalog-driven pool kernels (``kernel_homolysis``,
``kernel_recombination``, ``kernel_transfer``, ``kernel_move``).

The kernels move chains between pools; nothing about the move may be lost
or counted twice.  A kernel is refused, by name, where it cannot be
integrated honestly:

* on (or into) a side-group homolysis carrier: the carrier's removable-X
  inventory Z is seeded once from its mu1 and drains on its own, so a
  kernel that moves the carrier's chains strands or double-counts it;
* between pools of different repeat units (atoms or molar mass): the moment
  move keeps the unit count, not the atoms;
* on a pool connection a live generated reaction row also carries, or with
  one template listed twice: the same events would be counted twice.

A generated bond fission cutting a pool whose homolysis kernel carries that
channel never reaches that check: the row itself is refused (superseded by
channel identity) and the kernel runs. A row of another channel on the
kernel's connection still meets it.

What pure moments cannot represent is reported, never patched: the solver
records, per pool, every accepted step whose moments leave the realizable
set (mu1 >= mu0 >= 0, mu0*mu2 >= mu1^2), and every kernel term's validity
length (the shortest chain, in periods, its site law counts correctly),
declared or, for homolysis, derived (2*min_fragment - 1).
"""

import dataclasses
import json
import math

import numpy as np
import pytest

from rmgpy.kinetics import Arrhenius
from rmgpy.reaction import Reaction
from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig

from polymerKernelsTest import (PPU, _spc, bb_term, h_term, moments, r_term, rates5, system, system5, t_term)
from solverPolymerTest import _SGH_PARENT_MW, _sgh_core_and_pools, _sgh_system, _two_pool_species

C2H4 = (("C", 2), ("H", 4))
CH2 = (("C", 1), ("H", 2))


def _pool(label, i, **kw):
    return PolymerPoolConfig(label=label, xs=2, explicit_dp_to_species_index={},
                             mu_indices=(2 + 4 * i, 3 + 4 * i, 4 + 4 * i), **kw)


def _species(labels):
    core = [_spc("N#N", "N2")]
    for i, label in enumerate(labels):
        stem = "C" * (i + 1)
        core += [_spc("C" * (i + 4), label), _spc(stem + "O", f"{label}_mu0"), _spc(stem + "N", f"{label}_mu1"),
                 _spc(stem + "S", f"{label}_mu2")]
    return core


def _system(pools, chains, rxns=()):
    core = _species([p.label for p in pools])
    mask = np.array([s.label == "N2" for s in core], dtype=bool)
    rs = HybridPolymerSystem(T=700.0, P=1.0e5, initial_mole_fractions={core[0]: 1.0}, V_poly=1.0,
                             polymer_pools=pools, mass_transfer=[], gas_species_mask=mask,
                             constant_gas_volume=False,
                             initial_polymer_moments={p.label: moments(chains.get(p.label, {})) for p in pools},
                             termination=[])
    rs.initialize_model(core, list(rxns), [], [])
    return rs


def move(product="B", **kw):
    term = dict(template="move_t", k=1.0, a=0.0, b=1.0, periods_per_unit=PPU, product=product)
    term.update(kw)
    return term


# --------------------------------------------------------------------------
# side-group homolysis carriers
# --------------------------------------------------------------------------
def _sgh(carrier=None, other=None):
    """The SGH fixture's carrier 'poly' (1, 5, 30) beside an equal-MW pool Q;
    ``carrier`` / ``other`` are each pool's kernel fields."""
    core, mask, pools, mom = _sgh_core_and_pools()
    base = len(core)
    core += [_spc("CCCCC", "Q"), _spc("CCO", "Q_mu0"), _spc("CC=O", "Q_mu1"), _spc("CC#N", "Q_mu2")]
    mask = np.append(mask, [False] * 4)
    parent = dataclasses.replace(pools[0], **(carrier or {}))
    q = PolymerPoolConfig(label="Q", xs=2, explicit_dp_to_species_index={},
                          mu_indices=(base + 1, base + 2, base + 3), monomer_mw_g_mol=_SGH_PARENT_MW,
                          **(other or {}))
    mom["Q"] = (1.0, 5.0, 30.0)
    return _sgh_system(core, mask, [parent, q], mom)


@pytest.mark.parametrize("kernel, terms", [
    ("kernel_move", [move(product="Q")]),
    ("kernel_homolysis", [h_term(daughter="Q")]),
    ("kernel_recombination", [r_term(product="Q")]),
    ("kernel_transfer", [t_term(partner="poly", product="Q", partner_product="Q")]),
])
def test_a_catalog_kernel_on_an_sgh_carrier_is_refused(kernel, terms):
    """The review's case: a constant 1/s move off the carrier (Z = 9 of 10
    sites) gained 799 g/s of condensed + gas mass."""
    with pytest.raises(ValueError, match=rf"SGH-v2 carrier pool 'poly'.*{kernel}"):
        _sgh(carrier={kernel: terms})


@pytest.mark.parametrize("kernel, terms", [
    ("kernel_move", [move(product="poly")]),
    ("kernel_homolysis", [h_term(daughter="poly")]),
    ("kernel_recombination", [r_term(product="poly")]),
    ("kernel_transfer", [t_term(partner="Q", product="poly", partner_product="poly")]),
    ("kernel_transfer", [t_term(partner="poly", product="poly", partner_product="Q")]),
])
def test_a_catalog_kernel_into_an_sgh_carrier_is_refused(kernel, terms):
    with pytest.raises(ValueError, match=rf"SGH-v2 carrier pool 'poly'.*{kernel}.*'Q'"):
        _sgh(other={kernel: terms})


def test_a_carrier_beside_another_pool_still_builds():
    rs = _sgh()
    assert rs.sgh_enabled[0] == 1 and rs.sgh_enabled[1] == 0


# --------------------------------------------------------------------------
# repeat units across a connection
# --------------------------------------------------------------------------
def _units(kind, a_unit, b_unit, a_mw=0.0, b_mw=0.0):
    """Pools A (source) and B (target) with the given units; the kernel of
    ``kind`` connects them."""
    kernels = {"move": dict(kernel_move=[move()]),
               "homolysis": dict(kernel_homolysis=[h_term(daughter="B")]),
               "recombination": dict(kernel_recombination=[r_term(product="B")]),
               "transfer_product": dict(kernel_transfer=[t_term(partner="A", product="B", partner_product="A2")]),
               "transfer_partner": dict(kernel_transfer=[t_term(partner="B", product="A2", partner_product="A2")])}
    pools = [_pool("A", 0, repeat_unit_atoms=a_unit, monomer_mw_g_mol=a_mw, **kernels[kind]),
             _pool("B", 1, repeat_unit_atoms=b_unit, monomer_mw_g_mol=b_mw),
             _pool("A2", 2, repeat_unit_atoms=a_unit, monomer_mw_g_mol=a_mw)]
    return _system(pools, {"A": {40: 1.0}, "B": {40: 1.0}})


KINDS = ("move", "homolysis", "recombination", "transfer_product", "transfer_partner")


@pytest.mark.parametrize("kind", KINDS)
def test_a_kernel_between_pools_of_different_repeat_units_is_refused(kind):
    """The review's case: a C2 -> C1 move lost 10 mol C/s."""
    with pytest.raises(ValueError, match=r"repeat unit.*'A2?'.*'B'|repeat unit.*'B'.*'A2?'"):
        _units(kind, C2H4, CH2)


@pytest.mark.parametrize("kind", KINDS)
def test_a_kernel_between_pools_of_different_molar_mass_is_refused(kind):
    with pytest.raises(ValueError, match=r"repeat unit.*'A2?'.*'B'|repeat unit.*'B'.*'A2?'"):
        _units(kind, None, None, 28.054, 14.027)


def test_a_unit_declared_on_one_side_only_is_refused():
    with pytest.raises(ValueError, match=r"repeat unit.*'A'.*'B'"):
        _units("move", C2H4, None)


@pytest.mark.parametrize("kind", KINDS)
def test_a_kernel_between_pools_of_one_repeat_unit_builds(kind):
    rs = _units(kind, C2H4, C2H4, 28.054, 28.054)
    assert rs.polymer_pools[0].label == "A"


def test_the_kernel_test_pools_declare_no_unit_and_still_build():
    rs = system({60: 1.0}, {}, homolysis=[h_term()])
    assert rs.polymer_pools[0].repeat_unit_atoms is None


# --------------------------------------------------------------------------
# one event, one count
# --------------------------------------------------------------------------
def _migration(kernel=None, kernel_b=None, archetype=2, fission=False):
    """A -> B (a chain move), or with ``fission`` a generated bond fission
    A -> B + G stamped as generation stamps it."""
    sp, core, mask = _two_pool_species()
    if fission:
        # the daughter end keeps the fission's radical (a balanced row, I-255)
        sp["B"].molecule = _spc("[CH2]CCCC", "B").molecule
    rxn = Reaction(reactants=[sp["A"]], products=[sp["B"], sp["G"]] if fission else [sp["B"]],
                   kinetics=Arrhenius(A=(2.0, "1/s"), n=0.0, Ea=(0.0, "kcal/mol"), T0=(298.15, "K")),
                   reversible=False)
    rxn.polymer_flux_archetype = archetype
    if fission:
        rxn.polymer_reaction_channel = "bond_fission"
    rxn.is_end_group_reaction = True
    pa = PolymerPoolConfig(label="A", xs=2, explicit_dp_to_species_index={}, mu_indices=(1, 2, 3),
                           **(kernel or {}))
    pb = PolymerPoolConfig(label="B", xs=2, explicit_dp_to_species_index={}, mu_indices=(5, 6, 7),
                           **(kernel_b or {}))
    rs = HybridPolymerSystem(T=800.0, P=1.0e5, initial_mole_fractions={core[8]: 0.0}, V_poly=1.0,
                             polymer_pools=[pa, pb], mass_transfer=[], gas_species_mask=mask.copy(),
                             constant_gas_volume=False,
                             initial_polymer_moments={"A": (1.0, 10.0, 100.0), "B": (0.0, 0.0, 0.0)},
                             termination=[])
    rs.initialize_model(core, [rxn], [], [])
    return rs


def test_a_generated_row_on_a_kernels_connection_is_refused():
    """The review's case: an end-group migration A -> B at 2/s gives
    (-2, -20, -200) on A; the same move as a kernel doubled it."""
    rs = _migration()
    dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
    assert list(dn[1:4]) == pytest.approx([-2.0, -20.0, -200.0])
    with pytest.raises(ValueError, match=r"kernel_move.*move_t.*'A'.*'B'.*row 0"):
        _migration(kernel=dict(kernel_move=[move(product="B", k=2.0, b=1.0, periods_per_unit=1)]))


def test_a_fission_row_beside_a_homolysis_kernel_into_its_daughter_is_superseded():
    """The double-count ruling: the generated bond fission is refused by name
    (superseded, zero flux) and counted; the kernel stays."""
    rs = _migration(kernel=dict(kernel_homolysis=[h_term(daughter="B")]), archetype=3, fission=True)
    assert rs.reaction_refused[0] == 1
    assert [(c["pool"], c["kernels"], c["reason"]) for c in rs.kernel_supersession_census] == [
        ("A", ["kernel_homolysis"], "superseded:backbone_bond_fission")]
    assert len(rs.kernel_homolysis_terms[0]) == 1


def test_a_scission_row_of_another_channel_on_a_homolysis_kernels_connection_is_refused():
    """A chain move stamped scission is no bond fission: the channel rule
    leaves it live, and the count-once guard refuses the build."""
    with pytest.raises(ValueError, match=r"kernel_homolysis.*homolysis_t.*'A'.*'B'.*row 0"):
        _migration(kernel=dict(kernel_homolysis=[h_term(daughter="B")]), archetype=3)


def test_an_irreversible_row_and_the_reverse_kernel_are_not_an_overlap():
    rs = _migration(kernel_b=dict(kernel_move=[move(product="A")]))
    assert len(rs.kernel_move_terms[1]) == 1


@pytest.mark.parametrize("kernel", ["kernel_move", "kernel_transfer"])
def test_one_template_on_two_pools_is_refused(kernel):
    """A transfer listed on both of its pools fires on each."""
    terms = {"kernel_move": lambda this, other: [move(product=other, template="same_t")],
             "kernel_transfer": lambda this, other: [t_term(template="same_t", partner=other, product=other,
                                                            partner_product=this)]}[kernel]
    pools = [_pool("A", 0, **{kernel: terms("A", "B")}), _pool("B", 1, **{kernel: terms("B", "A")})]
    with pytest.raises(ValueError, match=rf"{kernel}.*'same_t'.*'A'.*'B'"):
        _system(pools, {"A": {40: 1.0}, "B": {40: 1.0}})


# --------------------------------------------------------------------------
# realizability census
# --------------------------------------------------------------------------
def _cone_exit():
    """The review's case: Re {18: 0.99, 218: 0.01} periods moves at n - 18
    per chain into Rm, whose moment flux has mu0*mu2 - mu1^2 < 0."""
    pools = [_pool("Re", 0, kernel_move=[move(product="Rm", a=1.0, b=-18.0)]), _pool("Rm", 1)]
    return _system(pools, {"Re": {18: 0.99, 218: 0.01}})


def test_the_kernel_keeps_its_flux_outside_the_realizable_set():
    """Intrinsic to pure moments: the daughter flux is not a distribution's
    (no change to the kernel)."""
    rs = _cone_exit()
    dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
    f = dn[6:9]
    assert list(f[:2]) == pytest.approx([0.01 * 200, 0.01 * 200 * 218 / PPU], rel=1e-12)
    assert f[0] * f[2] - f[1] ** 2 < 0.0


def test_the_census_counts_every_accepted_step_outside_the_realizable_set():
    rs = _cone_exit()
    rs.advance(1e-3)
    census = rs.realizability_census
    rm, re = census["Rm"], census["Re"]
    y = np.asarray(rs.y)
    assert y[6] * y[8] - y[7] ** 2 < 0.0
    assert rm["accepted_steps"] == re["accepted_steps"] >= 1
    assert rm["steps_outside"] >= 1
    assert 0.0 < rm["first_outside_t"] <= rm["last_outside_t"] == rs.t
    assert rm["worst_variance_deficit"] < 0.0
    assert rm["excursions"][0][0] == rm["first_outside_t"] and rm["excursions"][-1][1] is None
    assert re["steps_outside"] == 0 and re["first_outside_t"] is None and re["excursions"] == []


def test_the_census_counts_each_accepted_state_once():
    rs = _cone_exit()
    rs.advance(1e-3)
    n = rs.realizability_census["Rm"]["accepted_steps"]
    rs.record_realizability_census()        # the same accepted state again
    assert rs.realizability_census["Rm"]["accepted_steps"] == n
    rs.advance(2e-3)
    assert rs.realizability_census["Rm"]["accepted_steps"] > n


def test_a_realizable_run_records_nothing_outside():
    rs = system5({"Re": {40: 1.0, 60: 0.5}}, move={"Re": [bb_term()]})
    rs.advance(1e-2)
    assert all(c["steps_outside"] == 0 for c in rs.realizability_census.values())
    assert rs.realizability_census["Rm"]["accepted_steps"] >= 1


def test_the_census_summary_counts_the_pools_and_steps_outside():
    rs = _cone_exit()
    rs.advance(1e-3)
    summary = rs.realizability_summary()
    assert summary["pools_outside"] == ["Rm"]
    assert summary["steps_outside"] == {"Re": 0, "Rm": rs.realizability_census["Rm"]["steps_outside"]}
    assert summary["accepted_steps"] == rs.realizability_census["Rm"]["accepted_steps"]
    json.dumps(summary)


def test_the_profile_writer_records_the_census_beside_the_profile(tmp_path):
    from rmgpy.rmg.listener import SimulationProfileWriter
    rs = _cone_exit()
    rs.advance(1e-3)
    rs.snapshots = []
    (tmp_path / "solver").mkdir()
    SimulationProfileWriter(str(tmp_path), 0, _species(["Re", "Rm"])).update(rs)
    doc = json.loads((tmp_path / "solver" / "realizability_1_9.json").read_text())
    assert doc["summary"] == json.loads(json.dumps(rs.realizability_summary()))
    assert doc["pools"]["Rm"]["steps_outside"] == rs.realizability_census["Rm"]["steps_outside"]


# --------------------------------------------------------------------------
# validity lengths
# --------------------------------------------------------------------------
def test_homolysis_below_its_validity_length_is_reported_not_corrected():
    """The review's case: min_fragment 4, sites n - 7 on {4: 1, 20: 1}
    periods gives 1.50 events/s against the exact 1.95: the 4-period chain
    counts -3 sites.  Intrinsic; the term's validity length (7 periods) is
    reported."""
    rs = system({4: 1.0, 20: 1.0}, {}, homolysis=[h_term()])
    dn = rs.residual(0.0, rs.y, np.zeros_like(rs.y))[0]
    assert -dn[2] == pytest.approx(0.3 * 0.5 * (13 - 3))
    assert rs.kernel_validity == [{"pool": "P", "kernel": "kernel_homolysis", "template": "homolysis_t",
                                   "valid_from_periods": 2 * 4 - 1, "source": "derived"}]


def test_a_declared_validity_length_is_carried_and_checked():
    rs = system5({"Re": {40: 1.0}}, move={"Re": [bb_term(valid_from_periods=10)]})
    assert rs.kernel_validity == [{"pool": "Re", "kernel": "kernel_move", "template": "backbiting_t",
                                   "valid_from_periods": 10, "source": "declared"}]
    with pytest.raises(ValueError, match=r"valid_from_periods"):
        system5({"Re": {40: 1.0}}, move={"Re": [bb_term(valid_from_periods=0)]})
    with pytest.raises(ValueError, match=r"valid_from_periods=5.*2\*min_fragment - 1 = 7"):
        system({40: 1.0}, {}, homolysis=[h_term(valid_from_periods=5)])


def test_an_undeclared_validity_length_is_reported_as_undeclared():
    rs = system5({"Re": {40: 1.0}}, move={"Re": [bb_term()]})
    assert rs.kernel_validity == [{"pool": "Re", "kernel": "kernel_move", "template": "backbiting_t",
                                   "valid_from_periods": None, "source": "undeclared"}]
