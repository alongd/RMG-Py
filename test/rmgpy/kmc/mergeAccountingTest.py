"""Acceptance tests for observational RMG merge accounting."""

import copy
from fractions import Fraction
import os
from pathlib import Path

import pytest

from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import KineticsFamily
from rmgpy.kmc.merge_accounting import (
    MergeAccountingError,
    RecipeEditV2,
    RootBindingV2,
    RootedPathV2,
    _assert_no_reaction_aliases,
    _corrected_divisor,
    _root_weights,
    _trace_base_degeneracy,
    generate_reactions_with_merge_accounting,
)
from rmgpy.molecule.molecule import Molecule


DATABASE = Path(os.environ["RMG_DATABASE_PATH"])


@pytest.fixture(scope="module")
def kinetics_database():
    database = KineticsDatabase()
    database.load_families(
        str(DATABASE / "input/kinetics/families"),
        families=[
            "H_Abstraction",
            "R_Recombination",
            "R_Addition_MultipleBond",
        ],
    )
    return database


def _records(database, smiles):
    result = generate_reactions_with_merge_accounting(
        database,
        [Molecule(smiles=value) for value in smiles],
        only_families=["R_Recombination"],
        resonance=False,
    )
    assert result.reactions
    assert len(result.reactions) == len(result.records)
    return result.records


def test_full_species_identity_controls_identical_reactant_divisor(
    kinetics_database,
):
    identical = _records(kinetics_database, ["[CH3]", "[CH3]"])
    distinct_same_site_class = _records(kinetics_database, ["[CH3]", "C[CH2]"])

    assert {record.same_reactants for record in identical} == {2}
    assert {record.identical_reactant_divisor for record in identical} == {2.0}
    assert all(record.full_species_reactants_identical for record in identical)

    assert {record.same_reactants for record in distinct_same_site_class} == {
        0
    }
    assert {
        record.identical_reactant_divisor
        for record in distinct_same_site_class
    } == {1.0}
    assert not any(
        record.full_species_reactants_identical
        for record in distinct_same_site_class
    )


def test_root_weights_close_exactly_and_preserve_distinct_roots(
    kinetics_database,
):
    identical = _records(kinetics_database, ["[CH3]", "[CH3]"])[0]

    # RMG generates both swapped embeddings, then rejects one as an exact
    # duplicate of the same mapped root.
    assert len(identical.raw_paths) == 2
    assert (
        sum(
            path.contributes_to_pre_correction_sum
            for path in identical.raw_paths
        )
        == 1
    )
    assert len(identical.root_weights) == 1
    assert (
        sum(weight.weight for weight in identical.root_weights)
        == identical.G_final
    )

    records = generate_reactions_with_merge_accounting(
        kinetics_database,
        [Molecule(smiles="[CH3]"), Molecule(smiles="CCC")],
        only_families=["H_Abstraction"],
        resonance=False,
    ).records

    assert records
    assert any(len(record.root_weights) > 1 for record in records)
    for record in records:
        assert len({weight.root for weight in record.root_weights}) == len(
            record.root_weights
        )
        assert (
            sum(weight.weight for weight in record.root_weights)
            == record.G_final
        )
        assert {
            atom_id
            for graph in record.reactant_graphs
            for atom_id in graph.atom_ids
        } == {
            atom_id
            for graph in record.product_graphs
            for atom_id in graph.atom_ids
        }


def test_inconsistent_root_decomposition_is_refused():
    root = (RootBindingV2("*1", participant_index=0, atom_id=1),)
    recipe = (RecipeEditV2("GAIN_RADICAL", ("*1", 1)),)
    paths = (
        RootedPathV2("raw-000000", 1.0, root, recipe),
        RootedPathV2("raw-000001", 1.0, root, recipe),
    )

    with pytest.raises(MergeAccountingError, match="does not match G_final"):
        _root_weights(paths, divisor=1.0, final_degeneracy=1.0)


def test_divisor_is_derived_from_rmg_branch_not_observed_ratio():
    class DistinctReactant:
        def is_isomorphic(self, other):
            return False

    class ReactionStub:
        is_forward = True
        reactants = [DistinctReactant(), DistinctReactant()]

    with pytest.raises(MergeAccountingError, match="derived divisor"):
        _corrected_divisor(
            ReactionStub(),
            same_reactants=0,
            pre_correction_degeneracy=2.0,
            final_degeneracy=1.0,
        )


def test_role_distinct_paths_with_same_atom_set_remain_distinct_roots():
    recipe = (RecipeEditV2("BREAK_BOND", ("*1", 1, "*2")),)
    paths = (
        RootedPathV2(
            "raw-000000",
            1.0,
            (
                RootBindingV2("*1", 0, 1),
                RootBindingV2("*2", 0, 2),
            ),
            recipe,
        ),
        RootedPathV2(
            "raw-000001",
            1.0,
            (
                RootBindingV2("*1", 0, 2),
                RootBindingV2("*2", 0, 1),
            ),
            recipe,
        ),
    )

    weights = _root_weights(paths, divisor=1, final_degeneracy=2)

    assert len(weights) == 2
    assert weights[0].root != weights[1].root


def test_fractional_root_allocation_closes_exactly():
    recipe = (RecipeEditV2("GAIN_RADICAL", ("*1", 1)),)
    paths = tuple(
        RootedPathV2(
            "raw-{0:06d}".format(index),
            1.0,
            (RootBindingV2("*1", 0, index),),
            recipe,
        )
        for index in range(6)
    )

    weights = _root_weights(paths, divisor=6, final_degeneracy=1.0)

    assert {weight.weight for weight in weights} == {Fraction(1, 6)}
    assert sum(weight.weight for weight in weights) == Fraction(1, 1)


def test_repeated_reference_raw_paths_are_refused():
    reaction = object()

    with pytest.raises(
        MergeAccountingError, match="aliases an earlier occurrence"
    ):
        _assert_no_reaction_aliases([reaction, reaction])


def test_reverse_recalculation_refuses_unrooted_special_implementations():
    class ChargeTransferReaction:
        def is_charge_transfer_reaction(self):
            return True

    with pytest.raises(MergeAccountingError, match="charge-transfer"):
        _trace_base_degeneracy(
            object(), ChargeTransferReaction(), resonance=True
        )

    class OverriddenFamily(KineticsFamily):
        def calculate_degeneracy(self, reaction, resonance=True):
            return 999

    class OrdinaryReaction:
        def is_charge_transfer_reaction(self):
            return False

    with pytest.raises(MergeAccountingError, match="override"):
        _trace_base_degeneracy(
            OverriddenFamily(label="override"),
            OrdinaryReaction(),
            resonance=True,
        )


def _reaction_signature(reaction):
    def side(participants):
        return tuple(
            (
                participant.molecule[0].to_adjacency_list(remove_h=False),
                tuple(atom.id for atom in participant.molecule[0].atoms),
            )
            for participant in participants
        )

    reverse = getattr(reaction, "reverse", None)
    reverse_signature = (
        None
        if reverse is None
        else (
            side(reverse.reactants),
            side(reverse.products),
            tuple(reverse.template),
            reverse.degeneracy,
            reverse.family,
            reverse.is_forward,
            repr(reverse.kinetics),
        )
    )
    return (
        side(reaction.reactants),
        side(reaction.products),
        tuple(reaction.template),
        reaction.degeneracy,
        reaction.family,
        reaction.is_forward,
        repr(reaction.kinetics),
        reverse_signature,
    )


def _cold_reaction_signature(reaction):
    def side(participants):
        return tuple(
            participant.molecule[0].to_adjacency_list(remove_h=False)
            for participant in participants
        )

    return (
        side(reaction.reactants),
        side(reaction.products),
        tuple(reaction.template),
        reaction.degeneracy,
        reaction.family,
        reaction.is_forward,
    )


def _load_database(*families):
    database = KineticsDatabase()
    database.load_families(
        str(DATABASE / "input/kinetics/families"), families=list(families)
    )
    return database


@pytest.mark.parametrize("traced_first", [True, False])
def test_cold_tracing_parity_in_both_orders_with_resonance_and_failure(
    traced_first,
):
    order = (True, False) if traced_first else (False, True)
    results = {}
    for trace in order:
        database = _load_database("H_Abstraction")
        results[trace] = generate_reactions_with_merge_accounting(
            database,
            [Molecule(smiles="[CH2]C=C"), Molecule(smiles="C")],
            only_families=["H_Abstraction"],
            resonance=True,
            trace=trace,
        )

    assert [
        _cold_reaction_signature(item) for item in results[True].reactions
    ] == [_cold_reaction_signature(item) for item in results[False].reactions]
    assert results[True].records
    assert results[False].records == ()

    failures = {}
    for trace in order:
        database = _load_database("R_Recombination")
        with pytest.raises(AttributeError) as error:
            generate_reactions_with_merge_accounting(
                database,
                [Molecule(smiles="[CH3]"), Molecule(smiles="[CH3]")],
                products=[object()],
                only_families=["R_Recombination"],
                resonance=True,
                trace=trace,
            )
        failures[trace] = str(error.value)
    assert failures[True] == failures[False]


def test_resonant_bimolecular_paths_map_into_serialized_graphs(
    kinetics_database,
):
    result = generate_reactions_with_merge_accounting(
        kinetics_database,
        [Molecule(smiles="[CH2]C=C"), Molecule(smiles="C")],
        only_families=["H_Abstraction"],
        resonance=True,
    )

    assert result.records
    assert any(len(record.raw_paths) > 1 for record in result.records)
    for record in result.records:
        target_atoms = {
            (graph.participant_index, atom_index, atom_id)
            for graph in record.reactant_graphs
            for atom_index, atom_id in enumerate(graph.atom_ids)
        }
        for path in record.raw_paths:
            mapped_targets = {
                (
                    item.serialized_participant_index,
                    item.serialized_atom_index,
                    item.serialized_atom_id,
                )
                for item in path.atom_map
            }
            assert mapped_targets == target_atoms
            assert all(
                (item.raw_element, item.raw_isotope)
                == (item.serialized_element, item.serialized_isotope)
                for item in path.atom_map
            )
        assert (
            sum(weight.weight for weight in record.root_weights)
            == record.G_final
        )


def test_tracing_on_and_off_matches_native_generation_and_failures(
    kinetics_database,
):
    cases = (
        ("R_Recombination", ["[CH3]", "C[CH2]"]),
        ("R_Addition_MultipleBond", ["[CH2]CCC"]),
        ("H_Abstraction", ["[CH3]", "CCC"]),
    )
    for family, smiles in cases:
        reactants = [Molecule(smiles=value) for value in smiles]
        for molecule in reactants:
            molecule.assign_atom_ids()

        traced = generate_reactions_with_merge_accounting(
            kinetics_database,
            copy.deepcopy(reactants),
            only_families=[family],
            resonance=False,
            trace=True,
        )
        untraced = generate_reactions_with_merge_accounting(
            kinetics_database,
            copy.deepcopy(reactants),
            only_families=[family],
            resonance=False,
            trace=False,
        )

        assert [_reaction_signature(item) for item in traced.reactions] == [
            _reaction_signature(item) for item in untraced.reactions
        ]
        assert traced.records
        assert untraced.records == ()

    too_many = [Molecule(smiles="[CH3]") for _ in range(4)]
    failures = []
    for trace in (False, True):
        with pytest.raises(ValueError) as error:
            generate_reactions_with_merge_accounting(
                kinetics_database,
                copy.deepcopy(too_many),
                only_families=["R_Recombination"],
                resonance=False,
                trace=trace,
            )
        failures.append(str(error.value))
    assert failures[0] == failures[1]


def test_backward_non_own_reverse_uses_only_recalculated_rooted_paths(
    kinetics_database,
):
    result = generate_reactions_with_merge_accounting(
        kinetics_database,
        [Molecule(smiles="[CH2]CCC")],
        only_families=["R_Addition_MultipleBond"],
        resonance=False,
    )

    assert result.records
    assert all(record.direction == "reverse" for record in result.records)
    for record in result.records:
        assert not record.family_own_reverse
        assert record.reverse_degeneracy_recalculated
        assert record.recalculated_rooted_paths
        recalculated = {
            path.path_id for path in record.recalculated_rooted_paths
        }
        superseded = {path.path_id for path in record.raw_paths}
        allocated = {
            path_id
            for weight in record.root_weights
            for path_id in weight.contributing_path_ids
        }
        assert allocated <= recalculated
        assert allocated.isdisjoint(superseded)
        assert (
            sum(weight.weight for weight in record.root_weights)
            == record.G_final
        )


def test_small_ps_proxy_witness_is_accounted_without_special_case_logic(
    kinetics_database,
):
    from rmgpy.kmc.compiler import ps_proxy_set

    proxy = next(
        item
        for item in ps_proxy_set(3)
        if item.site_type == "end_radical+end_radical"
    )
    result = generate_reactions_with_merge_accounting(
        kinetics_database,
        copy.deepcopy(list(proxy.reactants)),
        only_families=["R_Recombination"],
        resonance=False,
    )

    assert result.records
    assert all(
        record.schema == "RMGMergeAccountingV2" for record in result.records
    )
    assert all(
        sum(weight.weight for weight in record.root_weights) == record.G_final
        for record in result.records
    )


def test_c_to_two_a_keeps_kc_degeneracy_and_unordered_counting_independent(
    kinetics_database,
):
    from scipy.constants import Avogadro

    from rmgpy.data.kinetics.family import TemplateReaction
    from rmgpy.kmc.event_record import ssa_multiplier_for
    from rmgpy.kinetics.arrhenius import Arrhenius
    from rmgpy.species import Species
    from rmgpy.thermo import ThermoData

    def thermo(enthalpy, entropy):
        return ThermoData(
            Tdata=([300, 400, 600, 800, 1000], "K"),
            Cpdata=([30, 30, 30, 30, 30], "J/(mol*K)"),
            H298=(enthalpy, "kJ/mol"),
            S298=(entropy, "J/(mol*K)"),
            Cp0=(30, "J/(mol*K)"),
            CpInf=(30, "J/(mol*K)"),
        )

    species_a = Species(
        label="A", molecule=[Molecule(smiles="[CH3]")], thermo=thermo(20, 150)
    )
    second_a = species_a.copy(deep=True)
    second_a.thermo = species_a.thermo
    species_c = Species(
        label="C", molecule=[Molecule(smiles="CC")], thermo=thermo(-50, 200)
    )
    source_record = _records(kinetics_database, ["[CH3]", "[CH3]"])[0]
    source_trace_degeneracy = source_record.G_final
    assert source_trace_degeneracy == sum(
        weight.weight for weight in source_record.root_weights
    )
    reaction = TemplateReaction(
        reactants=[species_c],
        products=[species_a, second_a],
        kinetics=Arrhenius(A=(10.0, "s^-1"), n=0, Ea=(0, "kJ/mol")),
        degeneracy=source_trace_degeneracy,
        is_forward=True,
    )

    temperature = 600.0
    volume = 2.5e-21
    k_forward = reaction.kinetics.get_rate_coefficient(temperature)
    equilibrium_constant = reaction.get_equilibrium_constant(
        temperature, type="Kc"
    )
    k_reverse = k_forward / equilibrium_constant

    forward_candidates = 1
    assert ssa_multiplier_for([species_c]) == 1.0
    assert k_forward * forward_candidates == k_forward

    unordered_pairs = 2 * (2 - 1) // 2
    assert unordered_pairs == 1
    assert ssa_multiplier_for([species_a, second_a]) == 2.0
    reverse_propensity = (
        ssa_multiplier_for([species_a, second_a])
        * k_reverse
        * unordered_pairs
        / (Avogadro * volume)
    )
    assert reverse_propensity == pytest.approx(
        2.0 * k_reverse / (Avogadro * volume)
    )
    assert k_reverse == k_forward / equilibrium_constant
    assert k_reverse != k_forward / (
        equilibrium_constant * source_trace_degeneracy
    )
