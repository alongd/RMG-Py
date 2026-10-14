"""Compiled artifact bytes depend on chemistry, not on generator hidden state.

Unit tests pin each structural field the compiler derives from one generated
reaction (graphs, atom map, bond rewrite, thermo source graph) against the
three kinds of hidden state the I-122 probe found reaching the artifact:

* raw atom ids, which come from RMG's process-global counter and wrap from
  32767 to -32768;
* the input atom order of a participant (and with it RMG's traversal order);
* which member of a degenerate set RMG keeps as the representative path.

The database-backed tests repeat the probe on a reduced real compile, each
perturbation in its own fresh interpreter.
"""

import copy
import hashlib
import json
import os
import random
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from rmgpy.kmc.compiler import (
    _bond_ops,
    _canonical_adjacency,
    _canonical_atom_map,
    _graph_adjacencies,
)
from rmgpy.molecule.molecule import Bond, Molecule
from rmgpy.species import Species

HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[2]
DATABASE = Path(os.environ.get("RMG_DATABASE_PATH", "/nonexistent"))
needs_database = pytest.mark.skipif(
    not (DATABASE / "input/kinetics/families").is_dir(),
    reason="RMG_DATABASE_PATH is not a database checkout",
)


# --- hand-built reactions ------------------------------------------------------


def _with_ids(molecule, identifiers):
    for atom, identifier in zip(molecule.atoms, identifiers):
        atom.id = identifier
    return molecule


def _by_id(molecule, identifier):
    return next(atom for atom in molecule.atoms if atom.id == identifier)


def _radical_carbon(molecule):
    return next(atom for atom in molecule.atoms if atom.radical_electrons)


def _recombination(first_smiles, second_smiles, first_ids, second_ids):
    """``A* + B* -> A-B`` with conserved atom ids and RMG-style *1/*2 labels."""
    first = _with_ids(Molecule(smiles=first_smiles), first_ids)
    second = _with_ids(Molecule(smiles=second_smiles), second_ids)
    _radical_carbon(first).label = "*1"
    _radical_carbon(second).label = "*2"
    product = first.copy(deep=True).merge(second.copy(deep=True))
    left = _by_id(product, _radical_carbon(first).id)
    right = _by_id(product, _radical_carbon(second).id)
    product.add_bond(Bond(left, right, 1))
    left.radical_electrons = right.radical_electrons = 0
    product.update_multiplicity()
    return SimpleNamespace(
        reactants=[Species(molecule=[first]), Species(molecule=[second])],
        products=[Species(molecule=[product])],
        family="R_Recombination",
    )


def _addition(first_ids, second_ids, attacked=0):
    """``CH3* + C=C -> CH3CH2CH2*`` with a cross-participant bond formation.

    ``attacked`` picks which (symmetry-equivalent) ethene carbon the methyl
    radical adds to: two representatives of one degenerate reaction.
    """
    methyl = _with_ids(Molecule(smiles="[CH3]"), first_ids)
    ethene = _with_ids(Molecule(smiles="C=C"), second_ids)
    carbon = _radical_carbon(methyl)
    ethene_carbons = sorted(
        (atom for atom in ethene.atoms if atom.is_carbon()), key=lambda atom: atom.id
    )
    target, partner = (
        ethene_carbons if attacked == 0 else list(reversed(ethene_carbons))
    )
    carbon.label, target.label, partner.label = "*3", "*1", "*2"
    product = methyl.copy(deep=True).merge(ethene.copy(deep=True))
    attacked_atom = _by_id(product, target.id)
    other = _by_id(product, partner.id)
    product.add_bond(Bond(_by_id(product, carbon.id), attacked_atom, 1))
    attacked_atom.edges[other].order = 1
    _by_id(product, carbon.id).radical_electrons = 0
    other.radical_electrons = 1
    product.update_multiplicity()
    return SimpleNamespace(
        reactants=[Species(molecule=[methyl]), Species(molecule=[ethene])],
        products=[Species(molecule=[product])],
        family="R_Addition_MultipleBond",
    )


def _structure(reaction):
    return {
        "reactant_graphs": _graph_adjacencies(reaction.reactants),
        "product_graphs": _graph_adjacencies(reaction.products),
        "atom_map": _canonical_atom_map(reaction),
        "bond_ops": _bond_ops(reaction),
    }


def _shuffled(reaction, seed):
    """The same reaction with every participant's atom list reordered."""
    rng = random.Random(seed)
    clone = copy.deepcopy(reaction)
    for participant in list(clone.reactants) + list(clone.products):
        for molecule in participant.molecule:
            atoms = list(molecule.atoms)
            rng.shuffle(atoms)
            molecule.atoms = atoms
    return clone


# --- raw atom ids ----------------------------------------------------------------


def test_bond_ops_do_not_depend_on_raw_atom_ids_or_the_counter_wrap():
    # ethene's ids above methyl's, then below it, as when RMG's counter wraps
    # from 32767 to -32768 between the two participants.
    ordered = _addition((100, 101, 102, 103), (200, 201, 202, 203, 204, 205))
    wrapped = _addition((32764, 32765, 32766, 32767), (-32768, -32767, -32766, -32765, -32764, -32763))
    assert _bond_ops(ordered) == _bond_ops(wrapped)


def test_bond_ops_do_not_depend_on_ids_wrapping_inside_one_participant():
    ordered = _recombination("[CH2]c1ccccc1", "[CH3]", range(1, 16), range(20, 24))
    wrapped = _recombination(
        "[CH2]c1ccccc1", "[CH3]", list(range(32760, 32768)) + list(range(-32768, -32761)), range(20, 24)
    )
    assert _bond_ops(ordered) == _bond_ops(wrapped)


def test_bond_ops_list_bonds_by_artifact_index_not_by_atom_id():
    reaction = _addition((300, 301, 302, 303), (-5, -4, -3, -2, -1, 0))
    pairs = [tuple(op["atoms"]) for op in _bond_ops(reaction) if "atoms" in op]
    assert all(first < second for first, second in pairs)
    assert pairs == sorted(pairs)


# --- input atom order / representative path ----------------------------------------


def _toluene_labelled_at_ortho(which, aromatic):
    molecule = Molecule(smiles="Cc1ccccc1")
    if aromatic:
        from rmgpy.molecule.resonance import generate_aromatic_resonance_structure

        molecule = generate_aromatic_resonance_structure(molecule)[0]
    methyl = next(atom for atom in molecule.atoms if atom.is_carbon() and not molecule.is_atom_in_cycle(atom))
    ipso = next(atom for atom in methyl.edges if atom.is_carbon())
    orthos = sorted(
        (atom for atom in ipso.edges if atom.is_carbon() and atom is not methyl),
        key=lambda atom: molecule.atoms.index(atom),
    )
    orthos[which].label = "*1"
    return molecule


@pytest.mark.parametrize("aromatic", [True, False])
@pytest.mark.parametrize("which", [0, 1])
def test_labelled_canonical_adjacency_is_independent_of_atom_order(which, aromatic):
    expected = _canonical_adjacency(_toluene_labelled_at_ortho(which, aromatic))
    if aromatic:
        # Either ortho carbon of the aromatic ring is the same labelled graph.
        assert _canonical_adjacency(_toluene_labelled_at_ortho(1 - which, True)) == expected
    for seed in range(12):
        clone = _toluene_labelled_at_ortho(which, aromatic)
        atoms = list(clone.atoms)
        random.Random(seed).shuffle(atoms)
        clone.atoms = atoms
        assert _canonical_adjacency(clone) == expected, seed


def test_kekule_canonical_adjacency_is_independent_of_atom_order():
    """RDKit aromatizes a Kekule ring; the two mirror drawings must not tie."""
    molecule = Molecule(smiles="Cc1ccccc1")
    molecule.kekulize()
    expected = _canonical_adjacency(molecule)
    for seed in range(12):
        clone = molecule.copy(deep=True)
        atoms = list(clone.atoms)
        random.Random(seed).shuffle(atoms)
        clone.atoms = atoms
        assert _canonical_adjacency(clone) == expected, seed


@pytest.mark.parametrize(
    "first, second",
    [
        ("[CH2]c1ccccc1", "[CH3]"),
        ("[CH2]c1ccccc1", "[CH2]c1ccccc1"),
        ("C[CH]c1ccccc1", "[CH2]CCc1ccccc1"),
    ],
)
def test_structural_fields_are_independent_of_participant_atom_order(first, second):
    reaction = _recombination(first, second, range(1, 100), range(200, 300))
    expected = _structure(reaction)
    for seed in range(12):
        assert _structure(_shuffled(reaction, seed)) == expected, seed


def test_structural_fields_are_the_same_for_every_symmetric_representative():
    """Adding to either ethene carbon is one degenerate reaction, two paths."""
    ids = (1, 2, 3, 4), (10, 11, 12, 13, 14, 15)
    assert _structure(_addition(*ids, attacked=1)) == _structure(_addition(*ids, attacked=0))


def _hydrogen_transfer(moving, labelled):
    """``CH4 + CH3* -> CH3* + CH4`` moving methane hydrogen ``moving``."""
    methane = _with_ids(Molecule(smiles="C"), (1, 2, 3, 4, 5))
    methyl = _with_ids(Molecule(smiles="[CH3]"), (10, 11, 12, 13))
    donor = next(atom for atom in methane.atoms if atom.is_carbon())
    acceptor = _radical_carbon(methyl)
    hydrogen = [atom for atom in methane.atoms if atom.is_hydrogen()][moving]
    donor.label, acceptor.label = "*1", "*3"
    if labelled:
        hydrogen.label = "*2"
    left = methane.copy(deep=True)
    right = methyl.copy(deep=True)
    moved = _by_id(left, hydrogen.id)
    left.remove_bond(moved.edges[_by_id(left, donor.id)])
    left.remove_atom(moved)
    right.add_atom(moved)
    right.add_bond(Bond(_by_id(right, acceptor.id), moved, 1))
    _by_id(left, donor.id).radical_electrons = 1
    _by_id(right, acceptor.id).radical_electrons = 0
    left.update_multiplicity()
    right.update_multiplicity()
    return SimpleNamespace(
        reactants=[Species(molecule=[methane]), Species(molecule=[methyl])],
        products=[Species(molecule=[left]), Species(molecule=[right])],
        family="H_Abstraction",
    )


@pytest.mark.parametrize("labelled", [True, False])
def test_the_moving_hydrogen_is_not_chosen_by_input_order(labelled):
    """Any of methane's four hydrogens is the same reaction, labelled or not."""
    expected = _structure(_hydrogen_transfer(0, labelled))
    for moving in range(1, 4):
        assert _structure(_hydrogen_transfer(moving, labelled)) == expected, moving
    for seed in range(8):
        assert _structure(_shuffled(_hydrogen_transfer(seed % 4, labelled), seed)) == expected


# --- thermo provenance --------------------------------------------------------------


class _ThermoDatabase:
    def get_thermo_data(self, species):
        from rmgpy.thermo import ThermoData

        return ThermoData(
            Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
            Cpdata=([1, 2, 3, 4, 5, 6, 7], "J/(mol*K)"),
            H298=(len(species.molecule[0].atoms), "kJ/mol"),
            S298=(1, "J/(mol*K)"),
            comment="Thermo group additivity estimation: test",
        )


def test_thermo_source_graph_is_independent_of_atom_order():
    from rmgpy.kmc.reference_thermo import SharedThermoAssignment

    reaction = _recombination("C[CH]c1ccccc1", "[CH3]", range(1, 100), range(200, 300))
    expected = SharedThermoAssignment(_ThermoDatabase(), None).assign_reaction(reaction)
    for seed in range(6):
        actual = SharedThermoAssignment(_ThermoDatabase(), None).assign_reaction(
            _shuffled(reaction, seed)
        )
        assert actual == expected, seed


def test_thermo_source_graph_ignores_which_kekule_forms_the_form_set_holds():
    """RMG's filtered resonance set depends on the starting form (I-122 probe).

    The same radical, reached by two degenerate paths, carried these two form
    sets: they differ only in which radical form the Kekule phenyl rides on.
    """
    from rmgpy.kmc.reference_thermo import SharedThermoAssignment

    shared = [
        "[CH2]C(CCC1=CC=CC[CH]1)c1ccccc1",
        "[CH2]C(CCC1=CCC=C[CH]1)c1ccccc1",
        "[CH2]C(CCC1=CC[CH]C=C1)c1ccccc1",
    ]
    first = shared + ["[CH2]C(CCC1=CC=CC[CH]1)C1=CC=CC=C1"]
    second = list(reversed(shared)) + ["[CH2]C(CCC1=CC[CH]C=C1)C1=CC=CC=C1"]

    def assignment(forms):
        species = Species(molecule=[Molecule(smiles=smiles) for smiles in forms])
        reaction = SimpleNamespace(reactants=[species], products=[])
        return SharedThermoAssignment(_ThermoDatabase(), None).assign_reaction(reaction)

    assert assignment(first)[0]["source_graph_sha256"] == assignment(second)[0][
        "source_graph_sha256"
    ]


# --- reduced real compile, one fresh interpreter per perturbation -----------------

DRIVER = HERE / "determinism_probe_driver.py"
# A reduced fixture: every PS family on the end radical and on its styrene
# addition pair (about a minute per cold compile).
FAMILIES = (
    "Disproportionation,H_Abstraction,R_Addition_MultipleBond,"
    "R_Recombination,intra_H_migration"
)
PROXIES = "end_radical,end_radical+styrene"


def _compile(tmp_path, name, perturbation, seed="0"):
    output = tmp_path / name
    output.mkdir()
    environment = dict(os.environ)
    environment.update(
        PYTHONHASHSEED=seed,
        PYTHONPATH=str(REPOSITORY) + os.pathsep + environment.get("PYTHONPATH", ""),
    )
    environment.pop("RMG_DATABASE_SHA", None)
    completed = subprocess.run(
        [
            sys.executable, str(DRIVER), str(DATABASE), str(output), perturbation,
            "--families", FAMILIES, "--proxies", PROXIES,
        ],
        cwd=REPOSITORY, env=environment, capture_output=True, text=True,
    )
    (output / "stdout.log").write_text(completed.stdout)
    (output / "stderr.log").write_text(completed.stderr)
    assert completed.returncode == 0, completed.stderr[-3000:]
    artifact = next(output.glob("*.json"))
    return artifact.read_bytes(), completed.stdout


@pytest.fixture(scope="module")
def baseline(tmp_path_factory):
    if not (DATABASE / "input/kinetics/families").is_dir():
        pytest.skip("RMG_DATABASE_PATH is not a database checkout")
    return _compile(tmp_path_factory.mktemp("baseline"), "baseline", "baseline")


def _record_key(record):
    return (
        record["family"],
        record["template"],
        tuple(record["reactant_graphs"]),
        tuple(record["product_graphs"]),
        record["orientation"],
    )


def _differences(expected: bytes, actual: bytes) -> str:
    """Summarize two artifacts by chemistry-keyed records, not by position."""
    first = {_record_key(record): record for record in json.loads(expected)["records"]}
    second = {_record_key(record): record for record in json.loads(actual)["records"]}
    fields = {}
    for key in set(first) & set(second):
        for name in sorted(set(first[key]) | set(second[key])):
            if first[key].get(name) != second[key].get(name):
                fields[name] = fields.get(name, 0) + 1
    return (
        f"sha {hashlib.sha256(expected).hexdigest()[:12]} != "
        f"{hashlib.sha256(actual).hexdigest()[:12]}; records only in one: "
        f"{len(set(first) ^ set(second))}; differing fields of shared records {fields}"
    )


@needs_database
@pytest.mark.parametrize(
    "perturbation, seed",
    [
        # The counter wraps between the end radical and styrene ids.
        ("counter_at_proxies:32467", "0"),
        # The counter wraps inside the short end radical.
        ("counter_at_proxies:32554", "0"),
        # The counter starts mid-range, so a later generation crosses the wrap.
        ("counter_start:12345", "0"),
        # RMG's own traversal order: participant atom lists reordered.
        ("permute_atoms", "0"),
        ("reverse_atoms", "0"),
        # A different string-hash seed reorders RMG's set iteration.
        ("baseline", "1"),
        # The generator's output list reordered (compiler candidate order).
        ("shuffle_output", "0"),
        ("reverse_reactants", "0"),
        # Allocator state.
        ("heap_noise", "0"),
    ],
)
def test_reduced_compile_bytes_ignore_generator_hidden_state(
    baseline, tmp_path, perturbation, seed
):
    expected, _ = baseline
    actual, stdout = _compile(tmp_path, "perturbed", perturbation, seed)
    assert actual == expected, f"{perturbation} (PYTHONHASHSEED={seed}): " + _differences(
        expected, actual
    )


@needs_database
def test_reshuffled_degenerate_paths_change_only_what_rmg_itself_returns(baseline, tmp_path):
    """RMG's degeneracy collapse keeps the first path of each degenerate set.

    Shuffling its input changes RMG's own output: which template a mixed-
    template duplicate set reports, and whether a representative keeps its
    reaction labels.  Every record whose chemistry RMG still reports the same
    (family, template, graphs, direction) must be byte-identical: the
    compiler adds no dependence on which path was kept.
    """
    expected, _ = baseline
    actual, _ = _compile(tmp_path, "perturbed", "shuffle_degenerate")

    def keyed(payload):
        records = {}
        for record in json.loads(payload)["records"]:
            records.setdefault(_record_key(record), []).append(
                {name: value for name, value in record.items() if name != "canonical_index"}
            )
        return records

    first, second = keyed(expected), keyed(actual)
    shared = set(first) & set(second)
    assert len(shared) >= 0.9 * len(first)
    differing = [key[:2] + key[4:] for key in shared if first[key] != second[key]]
    assert not differing, differing
