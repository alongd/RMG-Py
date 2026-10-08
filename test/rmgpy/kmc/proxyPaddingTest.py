from types import SimpleNamespace

import pytest

from rmgpy.kmc.proxy_padding import (
    BoundaryPort,
    PaddingLimitExceeded,
    RepeatUnitGraph,
    heavy_atom_distance,
    pad_molecule,
    pad_reaction_witness,
)
from rmgpy.molecule.molecule import Bond, Molecule


def _heavy_atoms(molecule):
    return [atom for atom in molecule.atoms if atom.element.number != 1]


def test_pe_artificial_boundary_is_padded_to_reacting_atom_distance():
    molecule = Molecule(smiles="CC")
    molecule.assign_atom_ids()
    backbone = _heavy_atoms(molecule)
    result = pad_molecule(
        molecule,
        [BoundaryPort(backbone[-1].id, "tail")],
        RepeatUnitGraph.from_smiles(
            "PE", "CC", head_atom_index=0, tail_atom_index=1
        ),
        reacting_atom_ids=[backbone[0].id],
        min_distance=5,
    )

    assert len(_heavy_atoms(result.molecule)) == 6
    assert result.extensions == 2
    assert heavy_atom_distance(
        result.molecule, backbone[0].id, result.boundaries[0].atom_id
    ) == 5


@pytest.mark.parametrize(
    "name,smiles,tail_index,expected_heavy_atoms",
    [
        ("PE", "CC", 1, 6),
        ("PP", "CC(C)", 1, 9),
        ("PS", "CC(c1ccccc1)", 1, 24),
    ],
)
def test_repeat_extension_interface_is_polymer_generic(
    name, smiles, tail_index, expected_heavy_atoms
):
    molecule = Molecule(smiles=smiles)
    molecule.assign_atom_ids()
    heavy = _heavy_atoms(molecule)
    result = pad_molecule(
        molecule,
        [BoundaryPort(heavy[tail_index].id, "tail")],
        RepeatUnitGraph.from_smiles(
            name, smiles, head_atom_index=0, tail_atom_index=tail_index
        ),
        reacting_atom_ids=[heavy[0].id],
        min_distance=5,
    )

    assert len(_heavy_atoms(result.molecule)) == expected_heavy_atoms
    assert heavy_atom_distance(
        result.molecule, heavy[0].id, result.boundaries[0].atom_id
    ) == 5


def test_physical_end_is_untouched_while_artificial_end_is_padded():
    molecule = Molecule(smiles="CC")
    molecule.assign_atom_ids()
    head, tail = _heavy_atoms(molecule)
    result = pad_molecule(
        molecule,
        [
            BoundaryPort(head.id, "head", kind="physical"),
            BoundaryPort(tail.id, "tail"),
        ],
        RepeatUnitGraph.from_smiles(
            "PE", "CC", head_atom_index=0, tail_atom_index=1
        ),
        reacting_atom_ids=[head.id],
        min_distance=5,
    )

    assert result.boundaries[0].atom_id == head.id
    assert result.boundaries[1].atom_id != tail.id
    assert len(_heavy_atoms(result.molecule)) == 6


def test_padding_is_idempotent_at_the_completed_boundary():
    molecule = Molecule(smiles="CC")
    molecule.assign_atom_ids()
    head, tail = _heavy_atoms(molecule)
    repeat = RepeatUnitGraph.from_smiles(
        "PE", "CC", head_atom_index=0, tail_atom_index=1
    )
    first = pad_molecule(
        molecule,
        [BoundaryPort(tail.id, "tail")],
        repeat,
        reacting_atom_ids=[head.id],
        min_distance=5,
    )
    second = pad_molecule(
        first.molecule,
        first.boundaries,
        repeat,
        reacting_atom_ids=[head.id],
        min_distance=5,
    )

    assert first.extensions == 2
    assert second.extensions == 0
    assert second.molecule.to_adjacency_list() == first.molecule.to_adjacency_list()


def test_reaction_padding_preserves_map_degeneracy_and_spectator():
    active = Molecule(smiles="CC")
    spectator = Molecule(smiles="[CH2]C(C)")
    active.assign_atom_ids()
    spectator.assign_atom_ids()
    active_head, active_tail = _heavy_atoms(active)
    spectator_tail = _heavy_atoms(spectator)[1]
    reaction = SimpleNamespace(
        reactants=[active, spectator],
        products=[active.copy(deep=True), spectator.copy(deep=True)],
        degeneracy=3.0,
        template=["root"],
    )
    spectator_before = spectator.to_adjacency_list()
    original_ids = {
        atom.id
        for participant in reaction.reactants
        for atom in _heavy_atoms(participant)
    }
    witness = pad_reaction_witness(
        reaction,
        [
            BoundaryPort(active_tail.id, "tail", participant_index=0),
            BoundaryPort(spectator_tail.id, "tail", participant_index=1),
        ],
        RepeatUnitGraph.from_smiles(
            "PE", "CC", head_atom_index=0, tail_atom_index=1
        ),
        min_distance=5,
        reacting_atom_ids=[active_head.id],
    )
    reactant_ids = {
        atom.id
        for participant in witness.reaction.reactants
        for atom in _heavy_atoms(participant)
    }
    product_ids = {
        atom.id
        for participant in witness.reaction.products
        for atom in _heavy_atoms(participant)
    }

    assert witness.extensions == 2
    assert original_ids < reactant_ids
    assert reactant_ids == product_ids
    assert witness.reaction.degeneracy == 3.0
    assert witness.reaction.reactants[1].to_adjacency_list() == spectator_before


def test_ambiguous_boundary_is_reported_not_guessed():
    molecule = Molecule(smiles="CC")
    molecule.assign_atom_ids()
    head, tail = _heavy_atoms(molecule)
    reaction = SimpleNamespace(
        reactants=[molecule],
        products=[molecule.copy(deep=True)],
        degeneracy=1.0,
    )
    repeat = RepeatUnitGraph.from_smiles(
        "PE", "CC", head_atom_index=0, tail_atom_index=1
    )
    completed = pad_reaction_witness(
        reaction,
        [BoundaryPort(tail.id, "tail")],
        repeat,
        min_distance=3,
        reacting_atom_ids=[head.id],
    )

    assert completed.extensions == 1
    with pytest.raises(ValueError, match="exactly one participant"):
        pad_reaction_witness(
            reaction,
            [BoundaryPort(999999, "tail")],
            repeat,
            min_distance=3,
            reacting_atom_ids=[head.id],
        )


def test_product_side_distance_also_triggers_bijective_padding():
    radical = Molecule(smiles="[CH3]")
    chain = Molecule(smiles="[CH2]C")
    for atom_id, atom in enumerate(radical.atoms + chain.atoms, 1):
        atom.id = atom_id
    radical_atom = _heavy_atoms(radical)[0]
    chain_head, chain_tail = _heavy_atoms(chain)
    product = radical.copy(deep=True).merge(chain.copy(deep=True))
    product_by_id = {atom.id: atom for atom in product.atoms}
    product.add_bond(
        Bond(product_by_id[radical_atom.id], product_by_id[chain_head.id], 1)
    )
    product_by_id[radical_atom.id].radical_electrons = 0
    product_by_id[chain_head.id].radical_electrons = 0
    product.update(sort_atoms=False)
    reaction = SimpleNamespace(
        reactants=[radical, chain],
        products=[product],
        degeneracy=1.0,
    )

    witness = pad_reaction_witness(
        reaction,
        [BoundaryPort(chain_tail.id, "tail", participant_index=1)],
        RepeatUnitGraph.from_smiles(
            "PE", "CC", head_atom_index=0, tail_atom_index=1
        ),
        min_distance=4,
        reacting_atom_ids=[radical_atom.id],
    )

    assert witness.extensions == 1
    assert heavy_atom_distance(
        witness.reaction.products[0],
        radical_atom.id,
        witness.boundaries[0].atom_id,
    ) == 4


def test_reaction_padding_fails_at_the_finite_graph_ceiling():
    molecule = Molecule(smiles="CC")
    molecule.assign_atom_ids()
    head, tail = _heavy_atoms(molecule)
    reaction = SimpleNamespace(
        reactants=[molecule],
        products=[molecule.copy(deep=True)],
        degeneracy=1.0,
    )

    with pytest.raises(PaddingLimitExceeded, match="maximum graph size"):
        pad_reaction_witness(
            reaction,
            [BoundaryPort(tail.id, "tail")],
            RepeatUnitGraph.from_smiles(
                "PE", "CC", head_atom_index=0, tail_atom_index=1
            ),
            min_distance=7,
            reacting_atom_ids=[head.id],
            max_heavy_atoms=4,
        )
