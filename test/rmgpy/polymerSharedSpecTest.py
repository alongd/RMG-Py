"""Verifier for the shared linear-homopolymer Step-2 contract."""

import json
from copy import deepcopy
from types import SimpleNamespace

import numpy as np
import pytest

from rmgpy import constants
from rmgpy.molecule import Molecule
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
from rmgpy.polymer_spec import LinearHomopolymerSpec, normalized_linear_homopolymer_spec
from rmgpy.quantity import Quantity
from rmgpy.rmg.polymer_input import (
    PolymerPhase,
    PolymerPhaseBlueprint,
    compile_polymer_phase,
)
from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig
from rmgpy.species import Species


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
            "bounds": {
                "elements": ["C", "H"],
                "isotopes": [-1],
                "formal_charges": [0],
                "radical_electrons": [0, 1],
                "lone_pairs": [0],
                "bond_orders": [1, 1.5, 2, 3],
                "spin_multiplicities": [1, 2, 3],
                "maximum_heavy_atom_degree": 4,
                "maximum_explicit_hydrogens_per_heavy_atom": 4,
                "maximum_total_hydrogens_on_free_small_species": 12,
                "maximum_family_arity": 2,
                "maximum_root_tuple_size": 4,
                "maximum_free_small_species_heavy_atoms": 8,
                "attachment_classes": [],
                "topology_metadata_schema": {},
            },
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


def _compiled_spec_phase(spec, *, density=(900.0, "kg/m^3"), chain_moles=None,
                         extra_species=()):
    monomer, ends = spec.legacy_molecules()
    pool = Polymer(
        label=spec.label,
        monomer=monomer,
        end_groups=ends,
        moments=spec.macroscopic_moments_mol,
        initial_mass=spec.initial_mass_kg,
        linear_homopolymer_spec=spec,
    )
    dummies = [Species(label=f"{spec.label}_mu{index}") for index in range(3)]
    species_dict = {pool.label: pool}
    species_dict.update({species.label: species for species in dummies})
    species_dict.update({species.label: species for species in extra_species})
    labels = [pool.label]
    if density is None:
        blueprint = PolymerPhaseBlueprint("phase", labels, "")
    else:
        blueprint = PolymerPhaseBlueprint("phase", labels, "", density=density)
    initial_moles = {
        pool: spec.macroscopic_moments_mol[0] if chain_moles is None else chain_moles
    }
    phase = compile_polymer_phase(blueprint, initial_moles, species_dict)
    return pool, dummies, phase, initial_moles, species_dict


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


def test_compile_phase_uses_realized_integer_ensemble_moments():
    spec = _spec(ensemble={
        "dp_number_fractions": {
            1: {"numerator": 1, "denominator": 4},
            2: {"numerator": 1, "denominator": 4},
            3: {"numerator": 1, "denominator": 2},
        },
        "representative_chain_count": 3,
        "rounding": "largest_remainder_dp_ascending_v1",
    })
    _pool, _dummies, phase, _initial, _species = _compiled_spec_phase(spec)
    assert phase.initial_moments[spec.label] == pytest.approx(
        spec.macroscopic_moments_mol
    )


def test_shared_spec_rejects_conflicting_density_and_chain_amount():
    spec = _spec()
    with pytest.raises(ValueError, match="density.*shared"):
        _compiled_spec_phase(spec, density=(1000.0, "kg/m^3"))
    with pytest.raises(ValueError, match="chain amount.*shared"):
        _compiled_spec_phase(
            spec,
            density=(spec.density_kg_m3, "kg/m^3"),
            chain_moles=2.0 * spec.macroscopic_moments_mol[0],
        )


def test_shared_spec_materializes_density_chain_and_mobile_amounts():
    value = _spec_mapping()
    value["small_species"] = {
        "water": {
            "graph": _graph([_atom("o", "O", hydrogens=2)]),
            "disposition": "mobile",
        }
    }
    value["initial_mobile_counts"] = {"water": 5}
    value["density_basis"] = "whole_initial_mixture"
    value["locality_policy"]["bounds"]["elements"].append("O")
    spec = LinearHomopolymerSpec.from_mapping(value)
    water = Species(label="water", molecule=[Molecule().from_smiles("O")])
    pool, _dummies, phase, initial_moles, _species = _compiled_spec_phase(
        spec, density=None, extra_species=(water,)
    )
    expected_mobile_mol = (
        spec.initialization().scale_factor * dict(spec.initial_mobile_counts)["water"]
        / constants.Na
    )
    assert phase.density.value_si == pytest.approx(spec.density_kg_m3)
    assert initial_moles[pool] == pytest.approx(spec.macroscopic_moments_mol[0])
    assert initial_moles[water] == pytest.approx(expected_mobile_mol)
    assert phase.initial_explicit[water] == pytest.approx(expected_mobile_mol)


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
    value["repeat_unit"]["atoms"][0]["hydrogens"] = 1
    value["repeat_unit"]["atoms"][0]["electronic_state"]["radical_electrons"] = 1
    spec = LinearHomopolymerSpec.from_mapping(value)
    assert spec.repeat_unit.composition == Composition.from_mapping({"C": 2, "H": 3})
    assert spec.repeat_unit.atoms[0].radical_electrons == 1
    assert set(dict(spec.repeat_unit.ports)) == {"head", "tail"}


def test_capped_dimer_rejects_inconsistent_physical_electronic_state():
    value = _spec_mapping()
    value["repeat_unit"]["atoms"][0]["electronic_state"]["radical_electrons"] = 1
    with pytest.raises(ValueError, match="electronic state"):
        LinearHomopolymerSpec.from_mapping(value)


def test_end_sensitive_forward_inverse_moment_equations_round_trip():
    moments = (2.0, 5.0, 13.0)
    observed = observables_from_moments(moments, 100.0, 10.0)
    assert observed.mass_g == pytest.approx(520.0)
    assert observed.mn_g_mol == pytest.approx(260.0)
    assert observed.mw_g_mol == pytest.approx(269.61538461538464)
    assert moments_from_observables(
        observed.mass_g, observed.mn_g_mol, observed.mw_g_mol, 100.0, 10.0
    ) == pytest.approx(moments)


def test_moment_validation_rejects_nonfinite_small_scale_and_bad_residuals():
    with pytest.raises(ValueError, match="finite"):
        observables_from_moments((float("nan"), float("nan"), float("nan")), 100.0, 10.0)
    with pytest.raises(ValueError, match="mu0\*mu2"):
        observables_from_moments((1.0e-9, 2.0e-9, 2.0e-9), 100.0, 10.0)
    with pytest.raises(ValueError, match="mu0\*mu2"):
        subtract_explicit_moments((2.0, 4.0, 8.0), {1: 0.5})


def test_depropagation_requires_actual_second_moment_debit():
    pool = EvolvedPoolMass(
        "PE", Composition.from_mapping({"C": 2, "H": 4}), (1.0, 2.0, 4.0)
    )
    with pytest.raises(ValueError, match="mu2_debit.*explicit"):
        pool.depropagate(1.0, "monomer")


def test_locality_policy_requires_complete_known_finite_bounds():
    missing = _spec_mapping()
    missing["locality_policy"] = {"radius": 2}
    with pytest.raises(ValueError, match="bounds"):
        LinearHomopolymerSpec.from_mapping(missing)
    unknown = _spec_mapping()
    unknown["locality_policy"]["bounds"]["mystery_bound"] = 3
    with pytest.raises(ValueError, match="unknown.*bound"):
        LinearHomopolymerSpec.from_mapping(unknown)


def test_dp_aliases_and_duplicate_json_keys_are_rejected():
    aliases = _spec_mapping(ensemble={
        "dp_number_fractions": {
            1: {"numerator": 1, "denominator": 2},
            "1": {"numerator": 1, "denominator": 2},
        },
        "representative_chain_count": 2,
        "rounding": "largest_remainder_dp_ascending_v1",
    })
    with pytest.raises(ValueError, match="duplicate DP"):
        LinearHomopolymerSpec.from_mapping(aliases)
    with pytest.raises(ValueError, match="duplicate JSON key"):
        normalized_linear_homopolymer_spec(
            '{"schema_version":1,"schema_version":1}'
        )


def test_decimal_fraction_tokens_remain_exact_and_isotopes_require_integers():
    value = _spec_mapping(ensemble={
        "dp_number_fractions": {
            1: {"numerator": 1, "denominator": 2},
            2: {"numerator": 1, "denominator": 2},
        },
        "representative_chain_count": 10,
        "rounding": "largest_remainder_dp_ascending_v1",
    })
    payload = json.dumps(value)
    payload = payload.replace(
        '{"numerator": 1, "denominator": 2}',
        '0.123456789012345678901',
        1,
    ).replace(
        '{"numerator": 1, "denominator": 2}',
        '0.876543210987654321099',
        1,
    )
    spec = normalized_linear_homopolymer_spec(payload)
    assert spec.initial_ensemble.fractions[0][1].numerator == 123456789012345678901
    assert spec.initial_ensemble.fractions[0][1].denominator == 10 ** 21
    with pytest.raises(ValueError, match="isotope.*integer"):
        Composition.from_mapping({("C", 13.9): 1})


def test_stationary_depletion_rejects_absent_elemental_inventory():
    pool = EvolvedPoolMass(
        "PE",
        Composition.from_mapping({"C": 2, "H": 4}),
        (1.0, 2.0, 4.0),
        end_inventory=(
            InventoryTerm("ends", Composition.from_mapping({"H": 2}), 1.0, True),
        ),
    )
    with pytest.raises(ValueError, match="available elemental inventory"):
        pool.deplete_stationary_inventory(
            "nonexistent-Br", Composition.from_mapping({"Br": 1}), 0.1, "Br-pool"
        )


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


def test_compiled_pool_carries_end_mass_into_live_solver_accounting():
    spec = _spec()
    pool, dummies, _phase, _initial, _species = _compiled_spec_phase(spec)
    compiled_pool = _phase.pools[0]
    core = [pool] + dummies
    config = compiled_pool.to_config({species: index for index, species in enumerate(core)})
    assert config.end_group_mw_g_mol == pytest.approx(
        spec.end_composition.molar_mass_g_mol()
    )
    mu0, mu1, _mu2 = spec.macroscopic_moments_mol
    assert config.condensed_mass_g(mu0, mu1) == pytest.approx(
        spec.initial_mass_kg * 1000.0
    )


def test_paired_depropagation_solver_and_ledger_close_mass_and_volume():
    """Production residual pair for pinned polymer.pyx:7304-7315."""
    repeat = Composition.from_mapping({"C": 2, "H": 4})
    ends = Composition.from_mapping({"H": 2})
    repeat_mass = repeat.molar_mass_g_mol()
    end_mass = ends.molar_mass_g_mol()
    core = [
        Species(label="PE", molecule=[Molecule().from_smiles("CCCC")]),
        Species(label="monomer", molecule=[Molecule().from_smiles("C=C")]),
        Species(label="PE_mu0"),
        Species(label="PE_mu1"),
        Species(label="PE_mu2"),
    ]
    config = PolymerPoolConfig(
        label="PE", xs=2, explicit_dp_to_species_index={}, mu_indices=(2, 3, 4),
        monomer_poly_index=1, monomer_mw_g_mol=repeat_mass,
        end_group_mw_g_mol=end_mass, k_unzip=1.0,
    )
    solver = HybridPolymerSystem(
        T=800.0, P=1.0e5, initial_mole_fractions={core[1]: 0.0}, V_poly=1.0,
        polymer_pools=[config], mass_transfer=[],
        gas_species_mask=np.array([False, True, False, False, False]),
        constant_gas_volume=False,
        initial_polymer_moments={"PE": (1.0, 2.0, 4.0)}, termination=[],
    )
    solver.initialize_model(core, [], [], [])
    assert solver.get_total_polymer_condensed_mass_g() == pytest.approx(
        config.condensed_mass_g(1.0, 2.0)
    )
    rates = np.asarray(solver.residual(0.0, solver.y, np.zeros_like(solver.y))[0])
    assert rates[[2, 3, 4]] == pytest.approx((0.0, -1.0, -3.0))
    assert rates[1] == pytest.approx(1.0)
    assert config.condensed_mass_g(rates[2], rates[3]) + rates[1] * repeat_mass == pytest.approx(0.0)

    ledger_pool = EvolvedPoolMass(
        "PE", repeat, (1.0, 2.0, 4.0),
        end_inventory=(InventoryTerm("ends", ends, 1.0, True),),
    )
    before = MassVolumeLedger((ledger_pool,), (), 900.0)
    updated, monomer = ledger_pool.depropagate(
        1.0, "released-monomer", mu2_debit=3.0
    )
    after = MassVolumeLedger((updated,), (monomer,), 900.0)
    assert updated.moments == pytest.approx((1.0, 1.0, 1.0))
    assert updated.losses == ledger_pool.losses
    assert monomer.amount_mol == pytest.approx(rates[1])
    assert before.mass_kg == pytest.approx(after.mass_kg)
    assert before.volume_m3 == pytest.approx(after.volume_m3)


def test_paired_complete_bundle_adapter_and_ledger_count_exact_graph_once():
    """Production phase-adapter pair for pinned polymer.pyx:5739."""
    spec = _spec()
    repeat = Composition.from_mapping({"C": 2, "H": 4})
    ends = Composition.from_mapping({"H": 2})
    per_chain_defect = InventoryTerm(
        "transferred-H-defect", Composition.from_mapping({"H": -1}), 0.4, True
    )
    ledger_pool = EvolvedPoolMass(
        "PE", repeat, (2.0, 5.0, 13.0),
        end_inventory=(InventoryTerm("ends", ends, 2.0, True),),
        defects=(per_chain_defect,),
    )
    before = MassVolumeLedger((ledger_pool,), (), 900.0)
    updated, explicit = ledger_pool.transfer_complete_chain_bundle(2, 0.4, "dp2")
    after = MassVolumeLedger((updated,), (explicit,), 900.0)
    assert updated.moments == pytest.approx((1.6, 4.2, 11.4))
    assert updated.losses == ledger_pool.losses
    assert updated.defects[0].amount_mol == pytest.approx(0.0)
    assert explicit.composition == repeat.scaled(2).plus(ends).plus(
        per_chain_defect.composition
    )
    assert before.mass_kg == pytest.approx(after.mass_kg)
    assert before.volume_m3 == pytest.approx(after.volume_m3)

    explicit_species = Species(
        label="defective_dp2", molecule=[Molecule().from_smiles("[CH2]CCC")]
    )
    proxy = SimpleNamespace(linear_homopolymer_spec=spec)
    adapter_pool = SimpleNamespace(
        label=spec.label, proxy_species=proxy, explicit_map={2: explicit_species},
        monomer=object(),
    )
    phase = PolymerPhase(
        density=Quantity((900.0, "kg/m^3")),
        initial_moments={spec.label: (2.0, 5.0, 13.0)},
        initial_explicit={explicit_species: 0.4}, pools=[adapter_pool],
    )
    tail = subtract_explicit_moments((2.0, 5.0, 13.0), {2: 0.4})
    production_tail = PolymerPoolConfig(
        label="PE", xs=2, explicit_dp_to_species_index={}, mu_indices=(0, 1, 2),
        monomer_mw_g_mol=repeat.molar_mass_g_mol(),
        end_group_mw_g_mol=ends.molar_mass_g_mol(),
    )
    explicit_mass_g = (
        0.4 * explicit_species.molecular_weight.value_si * constants.Na * 1000.0
    )
    expected_mass_g = production_tail.condensed_mass_g(tail[0], tail[1]) + explicit_mass_g
    assert phase.calculate_volume() == pytest.approx(expected_mass_g / 1000.0 / 900.0)
    assert phase.calculate_volume() == pytest.approx(after.volume_m3)


def test_paired_solver_mass_law_and_disjoint_depletion_close_exactly_once():
    """Production mass-law pair for pinned polymer.pyx:1180 and 5739."""
    repeat = Composition.from_mapping({"C": 2, "H": 4, "Br": 1})
    ends = Composition.from_mapping({"H": 2})
    defect = InventoryTerm(
        "per-chain-H-defect", Composition.from_mapping({"H": -1}), 2.0, True
    )
    ledger_pool = EvolvedPoolMass(
        "defect-bearing", repeat, (2.0, 4.0, 8.0),
        end_inventory=(InventoryTerm("ends", ends, 2.0, True),),
        defects=(defect,),
    )
    solver_pool = PolymerPoolConfig(
        label="defect-bearing", xs=2, explicit_dp_to_species_index={},
        mu_indices=(0, 1, 2), monomer_mw_g_mol=repeat.molar_mass_g_mol(),
        end_group_mw_g_mol=ends.molar_mass_g_mol(),
        chain_mass_defect_g_mol=Composition.from_mapping({"H": 1}).molar_mass_g_mol(),
    )
    before = MassVolumeLedger((ledger_pool,), (), 1100.0)
    assert solver_pool.condensed_mass_g(2.0, 4.0) == pytest.approx(
        ledger_pool.mass_g
    )
    updated, bromine = ledger_pool.deplete_stationary_inventory(
        "stationary-Br-depletion",
        Composition.from_mapping({"Br": 1}),
        0.5,
        "Br-pool",
    )
    after = MassVolumeLedger((updated,), (bromine,), 1100.0)
    remaining_solver_mass = solver_pool.condensed_mass_g(
        2.0, 4.0, stationary_depletion_g=bromine.mass_g
    )
    assert remaining_solver_mass == pytest.approx(updated.mass_g)
    assert remaining_solver_mass + bromine.mass_g == pytest.approx(
        solver_pool.condensed_mass_g(2.0, 4.0)
    )
    assert updated.defects == ledger_pool.defects
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
