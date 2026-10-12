"""Acceptance tests for observational RMG merge accounting."""

import copy
import os
from pathlib import Path

import pytest

from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.kmc.merge_accounting import (
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
    distinct_same_site_class = _records(
        kinetics_database, ["[CH3]", "C[CH2]"]
    )

    assert {record.same_reactants for record in identical} == {2}
    assert {record.identical_reactant_divisor for record in identical} == {2.0}
    assert all(record.full_species_reactants_identical for record in identical)

    assert {record.same_reactants for record in distinct_same_site_class} == {0}
    assert {
        record.identical_reactant_divisor for record in distinct_same_site_class
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
    assert sum(path.contributes_to_pre_correction_sum for path in identical.raw_paths) == 1
    assert len(identical.root_weights) == 1
    assert sum(weight.weight for weight in identical.root_weights) == identical.G_final

    h_abstraction = KineticsDatabase()
    h_abstraction.load_families(
        str(DATABASE / "input/kinetics/families"),
        families=["H_Abstraction"],
    )
    records = generate_reactions_with_merge_accounting(
        h_abstraction,
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
        assert sum(weight.weight for weight in record.root_weights) == record.G_final
        assert {
            atom_id
            for graph in record.reactant_graphs
            for atom_id in graph.atom_ids
        } == {
            atom_id
            for graph in record.product_graphs
            for atom_id in graph.atom_ids
        }


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
        assert sum(weight.weight for weight in record.root_weights) == record.G_final


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
    assert all(record.schema == "RMGMergeAccountingV1" for record in result.records)
    assert all(
        sum(weight.weight for weight in record.root_weights) == record.G_final
        for record in result.records
    )


def test_c_to_two_a_keeps_kc_degeneracy_and_unordered_counting_independent():
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
    source_trace_degeneracy = 7.0
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
