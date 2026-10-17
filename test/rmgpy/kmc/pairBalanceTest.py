"""Reciprocal record pairs are detailed-balanced by construction (I-148).

Convention under test (the one the template-emission cycle check applies):
``S_f / S_r = (ssa_f / ssa_r) * Kc`` with ``S = ssa_multiplier * k`` and ``Kc``
from the reference thermo provider, in the primary record's direction.  DB-free:
the provider is a stand-in with a closed-form Kc.
"""

import math
from dataclasses import replace

import pytest

import rmgpy.kmc.compiler as compiler_module
from rmgpy.kmc.compiler import EventSetCompiler
from rmgpy.kmc.event_record import EventRecord
from rmgpy.kmc.reference_thermo import ReferenceThermoResult, ThermoUnavailable
from rmgpy.molecule.molecule import Molecule

GRID = (600.0, 800.0, 1000.0)


def _adjacency(smiles):
    return Molecule(smiles=smiles).to_adjacency_list()


class _Provider:
    """Closed-form Kc(T) = 5 * exp(2000/T); records the directions asked for."""

    provenance = {"reference_thermo": "RMG gas-phase Kc", "stand_in": True}

    def __init__(self, fail=False):
        self.fail = fail
        self.asked = []

    def evaluate(self, reaction, temperatures):
        self.asked.append(
            (
                [s.molecule[0].to_smiles() for s in reaction.reactants],
                [s.molecule[0].to_smiles() for s in reaction.products],
            )
        )
        if self.fail:
            raise ThermoUnavailable("thermo unavailable for species [CH3]: no entry")
        return ReferenceThermoResult(
            equilibrium_constants=tuple(5.0 * math.exp(2000.0 / t) for t in temperatures),
            reaction_enthalpy_298_J_per_mol=-1.0e5,
            species_thermo_assignments=(),
        )


def _kc(temperature):
    return 5.0 * math.exp(2000.0 / temperature)


def _pair(reactants, products, k_primary, k_partner, *, ssa_primary=1.0, ssa_partner=1.0,
          degeneracy_partner=1.0):
    primary = EventRecord(
        family="X", template="t", arity=len(reactants), ssa_multiplier=ssa_primary,
        k_table={"T": list(GRID), "k": list(k_primary)}, rate_units="u",
        reactant_graphs=[_adjacency(s) for s in reactants],
        product_graphs=[_adjacency(s) for s in products],
        rate_source={"kind": "RMG family estimate", "available": True},
        reverse_of="evt_" + "0" * 64,
    )
    partner = EventRecord(
        family="X", template="t", arity=len(products), ssa_multiplier=ssa_partner,
        degeneracy=degeneracy_partner, orientation="reversed",
        k_table={"T": list(GRID), "k": list(k_partner)}, rate_units="u",
        reactant_graphs=[_adjacency(s) for s in products],
        product_graphs=[_adjacency(s) for s in reactants],
        rate_source={"kind": "RMG family estimate", "available": True,
                     "reference": "reverse from RMG gas-phase Kc",
                     "equilibrium_constant_table": {"T": list(GRID), "Kc": [1.0, 1.0, 1.0]}},
        reverse_of="evt_" + "0" * 64,
    )
    return primary, partner


def _balance(primary, partner, provider):
    return compiler_module.balance_reciprocal_pair(primary, partner, provider, GRID)


def _convention_error(primary, partner):
    worst = 0.0
    for i, t in enumerate(primary.k_table["T"]):
        s_f = primary.ssa_multiplier * primary.k_table["k"][i]
        s_r = partner.ssa_multiplier * partner.k_table["k"][i]
        want = primary.ssa_multiplier / partner.ssa_multiplier * _kc(t)
        worst = max(worst, abs(s_f / s_r / want - 1.0))
    return worst


def test_two_separate_estimates_with_degeneracies_2_and_1_come_out_balanced():
    k_f = [1.0e7, 2.0e7, 3.0e7]
    # the partner is its own family estimate, degeneracy 2 against 1: a flat factor 2
    k_r = [2.0 * k / _kc(t) for k, t in zip(k_f, GRID)]
    primary, partner = _pair(["[CH3]", "CC"], ["CCCC"], k_f, k_r, degeneracy_partner=2.0)
    primary2, partner2 = _balance(primary, partner, _Provider())
    assert primary2.k_table["k"] == k_f
    assert _convention_error(primary2, partner2) <= 1e-12
    assert partner2.k_table["k"] == pytest.approx([k / _kc(t) for k, t in zip(k_f, GRID)],
                                                  rel=1e-14)
    assert partner2.status == "enabled"


def test_rate_rule_forward_with_reverse_from_a_different_thermo_path_is_balanced():
    k_f = [4.0e6, 5.0e6, 6.0e6]
    # the stale reverse used a thermo whose Kc is 1.2 percent lower than the reference
    k_r = [k / (0.988 * _kc(t)) for k, t in zip(k_f, GRID)]
    primary, partner = _pair(["[CH3]", "CC"], ["CCCC"], k_f, k_r)
    primary = replace(primary, rate_source={"kind": "R-009 archived RMG rate rule",
                                            "available": True}, event_id="")
    primary2, partner2 = _balance(primary, partner, _Provider())
    assert primary2.k_table["k"] == k_f
    assert _convention_error(primary2, partner2) <= 1e-12
    # provenance: which side is primary, the convention, the reference thermo
    assert primary2.rate_source["pair_balance"]["role"] == "primary"
    note = partner2.rate_source["pair_balance"]
    assert note["role"] == "derived"
    assert "ssa" in note["convention"]
    assert partner2.thermo_provenance["reference_thermo"] == "RMG gas-phase Kc"
    assert partner2.thermo_provenance["equilibrium_constant_table"]["Kc"] == pytest.approx(
        [_kc(t) for t in GRID], rel=1e-14)
    # the stale compiler-thermo Kc is not carried on the derived record
    assert "equilibrium_constant_table" not in partner2.rate_source


def test_identical_reactant_pair_follows_the_ssa_convention():
    k_f = [1.0e7, 2.0e7, 3.0e7]
    primary, partner = _pair(["[CH3]", "[CH3]"], ["CC"], k_f, [1.0, 1.0, 1.0],
                             ssa_primary=2.0, ssa_partner=1.0)
    primary2, partner2 = _balance(primary, partner, _Provider())
    assert _convention_error(primary2, partner2) <= 1e-12
    # and in the other primary direction (unimolecular primary, ssa 2 on the partner)
    back, fwd = _pair(["CC"], ["[CH3]", "[CH3]"], k_f, [1.0, 1.0, 1.0],
                      ssa_primary=1.0, ssa_partner=2.0)
    back2, fwd2 = _balance(back, fwd, _Provider())
    assert _convention_error(back2, fwd2) <= 1e-12


def test_pair_with_a_participant_lacking_thermo_is_refused_with_a_reason():
    primary, partner = _pair(["[CH3]", "CC"], ["CCCC"], [1.0e7] * 3, [1.0e6] * 3)
    primary2, partner2 = _balance(primary, partner, _Provider(fail=True))
    for record in (primary2, partner2):
        assert record.status == "refused"
        assert "cannot be balanced" in record.status_reason
        assert "thermo unavailable for species [CH3]" in record.status_reason
    assert partner2.k_table is None


def test_derivation_is_asked_in_the_primary_direction():
    provider = _Provider()
    primary, partner = _pair(["[CH3]", "CC"], ["CCCC"], [1.0e7] * 3, [1.0e6] * 3)
    _balance(primary, partner, provider)
    assert provider.asked == [(["[CH3]", "CC"], ["CCCC"])]


def _junction_compiler(provider):
    return EventSetCompiler(None, [], ["R_Recombination"], reference_thermo_provider=provider,
                            temperature_grid=GRID, rmgpy_sha="t", rmg_database_sha="t")


def test_archived_para_pair_emitted_by_the_compiler_is_balanced():
    records = _junction_compiler(_Provider())._para_junction_records({})
    forward, reverse = records
    assert forward.k_table is not None and reverse.k_table is not None
    assert _convention_error(forward, reverse) <= 1e-12
    assert forward.rate_source["pair_balance"]["role"] == "primary"
    assert reverse.rate_source["pair_balance"]["role"] == "derived"


def test_archived_ortho_pairs_emitted_by_the_compiler_are_balanced_without_the_factor_2():
    records = _junction_compiler(_Provider())._ortho_junction_records({})
    assert len(records) == 4
    for forward, reverse in (records[0:2], records[2:4]):
        assert _convention_error(forward, reverse) <= 1e-12
        assert forward.degeneracy == 1.0


def test_finite_huge_kc_does_not_overflow_the_ssa_scaling():
    class Huge(_Provider):
        def evaluate(self, reaction, temperatures):
            return ReferenceThermoResult(
                equilibrium_constants=tuple(1e308 for _ in temperatures),
                reaction_enthalpy_298_J_per_mol=-1.0e5,
                species_thermo_assignments=(),
            )

    primary, partner = _pair(["[CH3]", "[CH3]"], ["CC"], [1e7] * 3, [1.0] * 3,
                             ssa_primary=2.0, ssa_partner=1.0)
    primary2, partner2 = _balance(primary, partner, Huge())
    assert partner2.status == "enabled"
    assert partner2.k_table["k"] == pytest.approx([1e7 / 1e308] * 3, rel=1e-14)
