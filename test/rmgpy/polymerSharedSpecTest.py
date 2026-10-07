"""Verifier for the shared linear-homopolymer Step-2 contract."""

import json
from copy import deepcopy
from types import SimpleNamespace

import pytest

from rmgpy import constants
from rmgpy.polymer import Polymer
from rmgpy.polymer_mass import (
    Composition,
    EvolvedPoolMass,
    ExplicitInventory,
    InventoryTerm,
    MassVolumeLedger,
    moments_from_observables,
    observables_from_moments,
    subtract_explicit_moments,
)
from rmgpy.polymer_spec import LinearHomopolymerSpec
from rmgpy.quantity import Quantity
from rmgpy.rmg.polymer_input import PolymerPhase


def _atom(key, element="C", hydrogens=0, isotope=None, radicals=0):
    value = {
        "key": key,
        "element": element,
        "hydrogens": hydrogens,
        "electronic_state": {"radical_electrons": radicals},
    }
    if isotope is not None:
        value["isotope"] = isotope
    return value


def _graph(atoms, bonds=(), ports=None):
    return {
        "atoms": atoms,
        "bonds": [{"atoms": list(pair), "order": order} for pair, order in bonds],
        "ports": ports or {},
    }


def _hydrogen_end(port_name):
    return _graph(
        [_atom("h", "H")],
        ports={port_name: {"atom": "h", "bond_order": 1}},
    )


def _pe_repeat():
    return _graph(
        [_atom("c1", hydrogens=2), _atom("c2", hydrogens=2)],
        [(('c1', 'c2'), 1)],
        {
            "head": {"atom": "c1", "bond_order": 1},
            "tail": {"atom": "c2", "bond_order": 1},
        },
    )


def _pp_repeat():
    return _graph(
        [_atom("c1", hydrogens=2), _atom("c2", hydrogens=1), _atom("m", hydrogens=3)],
        [(('c1', 'c2'), 1), (('c2', 'm'), 1)],
        {
            "head": {"atom": "c1", "bond_order": 1},
            "tail": {"atom": "c2", "bond_order": 1},
        },
    )


def _ps_repeat():
    atoms = [_atom("c1", hydrogens=2), _atom("c2", hydrogens=1)]
    atoms += [_atom(f"p{i}", hydrogens=0 if i == 1 else 1) for i in range(1, 7)]
    bonds = [(('c1', 'c2'), 1), (('c2', 'p1'), 1)]
    bonds += [((f"p{i}", f"p{i + 1}"), 1.5) for i in range(1, 6)]
    bonds.append((('p1', 'p6'), 1.5))
    return _graph(
        atoms,
        bonds,
        {
            "head": {"atom": "c1", "bond_order": 1},
            "tail": {"atom": "c2", "bond_order": 1},
        },
    )


def _alkene_product(repeat):
    product = deepcopy(repeat)
    product["ports"] = {}
    product["bonds"][0]["order"] = 2
    return product


def _spec_mapping(label="PE", repeat=None, ensemble=None, monomer_product=None):
    repeat = repeat or _pe_repeat()
    value = {
        "schema_version": 1,
        "label": label,
        "repeat_unit": repeat,
        "end_groups": {
            "head": _hydrogen_end("cap"),
            "tail": _hydrogen_end("cap"),
        },
        "small_species": {},
        "initial_ensemble": ensemble or {
            "dp_counts": {2: 1, 4: 2},
            "rounding": "largest_remainder_dp_ascending_v1",
        },
        "initial_mass_kg": 1.25,
        "density_kg_m3": 900.0,
        "locality_policy": {
            "radius": 2,
            "bounds": {"elements": ["C", "H"], "maximum_heavy_atom_degree": 4},
        },
        "family_manifest": {},
        "library_manifest": {},
        "transport_model": {},
        "diagnostics": {},
    }
    if monomer_product is not None:
        value["monomer_product"] = monomer_product
    return value


def _spec(**kwargs):
    return LinearHomopolymerSpec.from_mapping(_spec_mapping(**kwargs))


def test_n3_counts_canonicalize_to_reduced_thirds_reproducibly():
    value = _spec_mapping(ensemble={"dp_counts": {1: 1, 2: 1, 3: 1}})
    first = LinearHomopolymerSpec.from_mapping(value)
    second = LinearHomopolymerSpec.from_mapping(json.loads(first.normalized_json))
    fractions = first.as_mapping()["initial_ensemble"]["dp_number_fractions"]
    assert fractions == {
        "1": {"numerator": 1, "denominator": 3},
        "2": {"numerator": 1, "denominator": 3},
        "3": {"numerator": 1, "denominator": 3},
    }
    assert first.normalized_json == second.normalized_json


def test_largest_remainder_ties_break_by_ascending_dp():
    spec = _spec(ensemble={
        "dp_number_fractions": {
            1: {"numerator": 1, "denominator": 4},
            2: {"numerator": 1, "denominator": 4},
            3: {"numerator": 1, "denominator": 2},
        },
        "representative_chain_count": 3,
        "rounding": "largest_remainder_dp_ascending_v1",
    })
    assert spec.initial_ensemble.counts == {1: 1, 2: 1, 3: 1}


@pytest.mark.parametrize(
    "label,repeat,formula",
    [("PS", _ps_repeat(), {"C": 8, "H": 8}),
     ("PE", _pe_repeat(), {"C": 2, "H": 4}),
     ("PP", _pp_repeat(), {"C": 3, "H": 6})],
)
def test_ps_pe_pp_use_the_same_generic_declaration_path(label, repeat, formula):
    spec = _spec(label=label, repeat=repeat, monomer_product=_alkene_product(repeat))
    assert spec.repeat_unit.composition == Composition.from_mapping(formula)
    assert spec.monomer_product is not None
    assert type(spec).__name__ == "LinearHomopolymerSpec"


def test_ports_are_non_atomic_and_physical_radicals_are_separate():
    value = _spec_mapping()
    value["repeat_unit"]["atoms"][0]["electronic_state"]["radical_electrons"] = 1
    spec = LinearHomopolymerSpec.from_mapping(value)
    assert spec.repeat_unit.composition == Composition.from_mapping({"C": 2, "H": 4})
    assert spec.repeat_unit.atoms[0].radical_electrons == 1
    assert set(dict(spec.repeat_unit.ports)) == {"head", "tail"}


def test_end_sensitive_forward_inverse_moment_equations_round_trip():
    moments = (2.0, 5.0, 13.0)
    observed = observables_from_moments(moments, 100.0, 10.0)
    assert observed.mass_g == pytest.approx(520.0)
    assert observed.mn_g_mol == pytest.approx(260.0)
    assert observed.mw_g_mol == pytest.approx(269.61538461538464)
    assert moments_from_observables(
        observed.mass_g, observed.mn_g_mol, observed.mw_g_mol, 100.0, 10.0
    ) == pytest.approx(moments)


def test_ps_formula_uses_end_corrected_mass_and_realized_observables():
    spec = _spec(label="PS", repeat=_ps_repeat())
    counts = spec.initial_ensemble.counts
    expected_mn = sum(
        count * (dp * Composition.from_mapping({"C": 8, "H": 8}).molar_mass_g_mol()
                 + Composition.from_mapping({"H": 2}).molar_mass_g_mol())
        for dp, count in counts.items()
    ) / sum(counts.values())
    assert spec.ensemble_observables.mn_g_mol == pytest.approx(expected_mn)
    round_trip = observables_from_moments(
        spec.macroscopic_moments_mol,
        spec.repeat_unit.composition.molar_mass_g_mol(),
        spec.end_composition.molar_mass_g_mol(),
    )
    assert round_trip.mass_g == pytest.approx(spec.initial_mass_kg * 1000.0)


def test_explicit_chains_subtract_all_moments_and_are_counted_once():
    total = (3.0, 8.0, 22.0)  # DP 2, 2, 4
    tail = subtract_explicit_moments(total, {2: 1.0})
    assert tail == pytest.approx((2.0, 6.0, 18.0))
    repeat = Composition.from_mapping({"C": 2, "H": 4})
    ends = Composition.from_mapping({"H": 2})
    pool = EvolvedPoolMass(
        "PE", repeat, tail,
        end_inventory=(InventoryTerm("ends", ends, 2.0, True),),
    )
    explicit = ExplicitInventory("dp2", repeat.scaled(2).plus(ends), 1.0)
    expected = observables_from_moments(
        total, repeat.molar_mass_g_mol(), ends.molar_mass_g_mol()
    ).mass_g
    assert pool.mass_g + explicit.mass_g == pytest.approx(expected)


def test_kmc_and_mom_consume_byte_identical_json_and_close_volume():
    spec = _spec()
    kmc = spec.for_kmc()
    mom = spec.for_method_of_moments()
    assert kmc.normalized_json == mom.normalized_json == spec.normalized_json
    assert kmc.moments_mol == pytest.approx(mom.moments_mol)
    assert kmc.volume_m3 == pytest.approx(
        kmc.representative_polymer_mass_kg / spec.density_kg_m3
    )
    assert kmc.scale_factor * kmc.representative_polymer_mass_kg == pytest.approx(
        spec.initial_mass_kg
    )


def test_mixture_density_requires_and_counts_explicit_mobile_inventory():
    value = _spec_mapping()
    value["small_species"] = {
        "water": {
            "graph": _graph([_atom("o", "O", hydrogens=2)]),
            "disposition": "mobile",
        }
    }
    value["initial_mobile_counts"] = {"water": 5}
    with pytest.raises(ValueError, match="whole_initial_mixture"):
        LinearHomopolymerSpec.from_mapping(value)
    value["density_basis"] = "whole_initial_mixture"
    spec = LinearHomopolymerSpec.from_mapping(value)
    initialized = spec.initialization()
    expected_mobile = (
        5 * Composition.from_mapping({"H": 2, "O": 1}).molar_mass_kg_mol()
        / constants.Na
    )
    assert initialized.representative_mobile_mass_kg == pytest.approx(expected_mobile)
    assert initialized.volume_m3 == pytest.approx(
        (initialized.representative_polymer_mass_kg + expected_mobile)
        / spec.density_kg_m3
    )


def test_correction5_physical_isotope_mass_and_unknown_rejection_delta():
    labelled = Composition.from_mapping({("C", 13): 1, "H": 4})
    assert labelled.legacy_nominal_mass_g_mol() == 16.0
    assert labelled.molar_mass_g_mol() == pytest.approx(17.035238, rel=2.0e-6)
    unknown = Composition.from_mapping({"Qq": 1})
    assert unknown.legacy_nominal_mass_g_mol() == 0.0
    with pytest.raises(ValueError, match="unknown element/isotope"):
        unknown.molar_mass_g_mol()


def test_mom_adapter_uses_end_corrected_shared_mass_without_legacy_fallback():
    spec = _spec()
    monomer, ends = spec.legacy_molecules()
    pool = Polymer(
        label=spec.label,
        monomer=monomer,
        end_groups=ends,
        moments=spec.macroscopic_moments_mol,
        initial_mass=spec.initial_mass_kg,
        linear_homopolymer_spec=spec,
    )
    assert pool.monomer_mw_g_mol == pytest.approx(
        spec.repeat_unit.composition.molar_mass_g_mol()
    )
    assert pool.end_group_mw_g_mol == pytest.approx(spec.end_composition.molar_mass_g_mol())
    assert (pool.Mn, pool.Mw) == pytest.approx(
        (spec.ensemble_observables.mn_g_mol, spec.ensemble_observables.mw_g_mol)
    )
    assert pool.linear_homopolymer_spec.normalized_json == spec.normalized_json


def test_mom_phase_volume_delegates_to_shared_end_corrected_mass():
    spec = _spec()
    moments = spec.macroscopic_moments_mol
    proxy = SimpleNamespace(linear_homopolymer_spec=spec)
    pool = SimpleNamespace(
        label=spec.label,
        proxy_species=proxy,
        explicit_map={},
        monomer=object(),
    )
    phase = PolymerPhase(
        density=Quantity((spec.density_kg_m3, "kg/m^3")),
        initial_moments={spec.label: moments},
        initial_explicit={},
        pools=[pool],
    )
    assert phase.calculate_volume() == pytest.approx(
        spec.initial_mass_kg / spec.density_kg_m3
    )


def test_combined_depropagation_credits_monomer_and_leaves_loss_ledger_unchanged():
    """Combined update pin: c220aa83:polymer.pyx:7304-7315."""
    repeat = Composition.from_mapping({"C": 2, "H": 4})
    ends = Composition.from_mapping({"H": 2})
    pool = EvolvedPoolMass(
        "PE", repeat, (1.0, 2.0, 4.0),
        end_inventory=(InventoryTerm("ends", ends, 1.0, True),),
    )
    before = MassVolumeLedger((pool,), (), 900.0)
    updated, monomer = pool.depropagate(
        0.25, "released-monomer", mu2_debit=0.75
    )
    after = MassVolumeLedger((updated,), (monomer,), 900.0)
    assert updated.losses == pool.losses
    assert before.mass_kg == pytest.approx(after.mass_kg)
    assert before.volume_m3 == pytest.approx(after.volume_m3)


def test_combined_complete_bundle_transfer_has_no_duplicate_loss():
    """Combined update pin: c220aa83:polymer.pyx:5739 and explicit handshake."""
    repeat = Composition.from_mapping({"C": 2, "H": 4})
    ends = Composition.from_mapping({"H": 2})
    stationary_loss = InventoryTerm(
        "stationary-Br", Composition.from_mapping({"Br": 1}), 0.1
    )
    per_chain_defect = InventoryTerm(
        "per-chain-H-defect", Composition.from_mapping({"H": -1}), 2.0, True
    )
    pool = EvolvedPoolMass(
        "PE", repeat, (2.0, 5.0, 13.0),
        end_inventory=(InventoryTerm("ends", ends, 2.0, True),),
        defects=(per_chain_defect,),
        losses=(stationary_loss,),
    )
    before = MassVolumeLedger((pool,), (), 900.0)
    updated, explicit = pool.transfer_complete_chain_bundle(2, 0.4, "dp2")
    after = MassVolumeLedger((updated,), (explicit,), 900.0)
    assert updated.losses == pool.losses
    assert updated.defects[0].amount_mol == pytest.approx(1.6)
    assert explicit.composition == repeat.scaled(2).plus(ends).plus(
        per_chain_defect.composition
    )
    assert before.mass_kg == pytest.approx(after.mass_kg)
    assert before.volume_m3 == pytest.approx(after.volume_m3)


def test_combined_defect_pool_and_disjoint_stationary_depletion_close_once():
    """Combined update pins: c220aa83:polymer.pyx:1180 and 5739."""
    repeat = Composition.from_mapping({"C": 2, "H": 4, "Br": 1})
    ends = Composition.from_mapping({"H": 2})
    defect = InventoryTerm(
        "per-chain-H-defect", Composition.from_mapping({"H": -1}), 2.0, True
    )
    pool = EvolvedPoolMass(
        "defect-bearing", repeat, (2.0, 4.0, 8.0),
        end_inventory=(InventoryTerm("ends", ends, 2.0, True),),
        defects=(defect,),
    )
    before = MassVolumeLedger((pool,), (), 1100.0)
    updated, bromine = pool.deplete_stationary_inventory(
        "stationary-Br-depletion",
        Composition.from_mapping({"Br": 1}),
        0.5,
        "Br-pool",
    )
    after = MassVolumeLedger((updated,), (bromine,), 1100.0)
    assert updated.defects == pool.defects
    assert {term.provenance_id for term in updated.defects}.isdisjoint(
        term.provenance_id for term in updated.losses
    )
    assert before.mass_kg == pytest.approx(after.mass_kg)
    assert before.volume_m3 == pytest.approx(after.volume_m3)


def test_per_pool_and_per_chain_defects_compose_with_distinct_losses():
    repeat = Composition.from_mapping({"C": 2, "H": 4})
    pool = EvolvedPoolMass(
        "mixed-defects", repeat, (2.0, 5.0, 13.0),
        defects=(
            InventoryTerm("per-chain", Composition.from_mapping({"H": -1}), 2.0, True),
            InventoryTerm("per-pool", Composition.from_mapping({"O": 1}), 0.2, False),
        ),
        losses=(InventoryTerm("stationary-loss", Composition.from_mapping({"H": 1}), 0.1),),
    )
    expected = (
        5.0 * repeat.molar_mass_g_mol()
        - 2.0 * Composition.from_mapping({"H": 1}).molar_mass_g_mol()
        + 0.2 * Composition.from_mapping({"O": 1}).molar_mass_g_mol()
        - 0.1 * Composition.from_mapping({"H": 1}).molar_mass_g_mol()
    )
    assert pool.mass_g == pytest.approx(expected)


def test_unknown_schema_and_malformed_ports_are_rejected():
    future = _spec_mapping()
    future["schema_version"] = 2
    with pytest.raises(ValueError, match="unknown.*schema_version"):
        LinearHomopolymerSpec.from_mapping(future)
    malformed = _spec_mapping()
    malformed["repeat_unit"]["ports"].pop("tail")
    with pytest.raises(ValueError, match="requires exactly ports"):
        LinearHomopolymerSpec.from_mapping(malformed)
