"""Fast contract tests for the offline event-set compiler."""

import copy
import hashlib
import subprocess
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest

import rmgpy.kmc.compiler as compiler_module
from rmgpy.kmc.compiler import (
    EventSetCompiler,
    PS_FAMILY_FILTER_REASON,
    SiteProxy,
    _archived_ortho_reaction,
    _canonical_adjacency,
    _graph_list_key,
    _implicit_h,
    _molecule,
    _ortho_junction_proxy,
    _orient_to_proxy,
    _reflect_ortho_reaction,
    ceiling_temperature,
    ps_proxy_set,
    short_ps_molecule_catalogue,
    validate_artifact,
)


class _Element:
    number = 6


class _Atom:
    element = _Element()
    radical_electrons = 0
    implicit_hydrogens = 0

    def __init__(self, identifier):
        self.id = identifier


class _Molecule:
    def __init__(self):
        self.atoms = [_Atom(1), _Atom(2)]

    def copy(self, deep=True):
        return _Molecule()

    def get_element_count(self):
        return {"C": 2, "H": 6}

    def get_all_edges(self):
        return []

    def to_adjacency_list(self, remove_h=False):
        return "1 C u0\n2 C u0"


class _Species:
    def __init__(self, molecule):
        self.molecule = [molecule]


class _Kinetics:
    units = "s^-1"

    def get_rate_coefficient(self, temperature):
        return temperature / 100.0


class _Reaction:
    family = "H_Abstraction"
    template = []
    degeneracy = 2.0
    kinetics = _Kinetics()

    def __init__(self, molecule):
        self.reactants = [_Species(molecule)]
        self.products = [_Species(molecule.copy(deep=True))]


class _Database:
    def __init__(self, reaction):
        self.reaction = reaction
        self.calls = []
        self.families = {
            "H_Abstraction": object(),
            "Loaded_But_Filtered": object(),
            "No_Firing_Path": object(),
        }

    def generate_reactions_from_families(
        self, reactants, products=None, only_families=None, resonance=True
    ):
        self.calls.append((only_families, resonance))
        return [self.reaction] if self.reaction.family in only_families else []


def _compiler():
    molecule = _Molecule()
    database = _Database(_Reaction(molecule))
    return (
        EventSetCompiler(
            database,
            [SiteProxy("interior", [_Species(molecule)])],
            ["H_Abstraction"],
            temperature_grid=[600.0, 700.0],
        ),
        database,
    )


def test_uses_public_family_pipeline_and_emits_f1_fields():
    compiler, database = _compiler()
    artifact = compiler.compile()
    expected_rmgpy_sha = subprocess.check_output(
        ["git", "-C", str(compiler.rmgpy_path), "rev-parse", "HEAD"], text=True
    ).strip()
    assert (
        artifact["provenance"]["compiler_sha256"]
        == hashlib.sha256(Path(compiler_module.__file__).read_bytes()).hexdigest()
    )
    assert artifact["provenance"]["rmgpy_sha"] == expected_rmgpy_sha
    assert database.calls == [(["H_Abstraction"], True)]
    record = artifact["records"][0]
    assert record["site_type"] == "interior"
    assert record["raw_path_degeneracy"] == 2.0
    assert record["ssa_multiplier"] == 1.0
    assert record["k_table"]["interpolation"] == "linear-ln-k"
    assert record["k_table"]["extrapolation"] == "refuse"
    assert record["formula_delta"] == {"C": 0, "H": 0}
    assert record["provenance"]["family_list_sha256"]
    assert "archived_j_para_rate" not in artifact["provenance"]
    assert "archived_j_para_rate" not in record["provenance"]


def test_family_filter_accounts_for_every_loaded_family_and_keeps_only_firing():
    molecule = _Molecule()
    database = _Database(_Reaction(molecule))
    proxy = SiteProxy("featured", [_Species(molecule)])

    active, excluded, reactions = EventSetCompiler.discover_family_reactions(
        database,
        [proxy],
        candidate_families=["H_Abstraction", "No_Firing_Path"],
        family_universe=[*database.families, "Filesystem_Only"],
    )

    assert active == ["H_Abstraction"]
    assert set(active) | set(excluded) == {
        *database.families,
        "Filesystem_Only",
    }
    assert excluded["Loaded_But_Filtered"] == PS_FAMILY_FILTER_REASON
    assert excluded["Filesystem_Only"] == PS_FAMILY_FILTER_REASON
    assert "bounded L=3 PS proxy set" in excluded["No_Firing_Path"]
    assert [reaction.family for reaction in reactions["featured"]] == active
    assert database.calls == [
        (["H_Abstraction", "No_Firing_Path"], True),
    ]


def test_orientation_uses_rmg_direction_for_nonisomorphic_resonance_form():
    from rmgpy.molecule.molecule import Molecule
    from rmgpy.species import Species

    proxy = SiteProxy(
        "end_radical",
        [Species(molecule=[Molecule(smiles="[CH2]C=C")])],
    )
    reaction = SimpleNamespace(
        reactants=[Species(molecule=[Molecule(smiles="C=C")])],
        products=[Species(molecule=[Molecule(smiles="[CH3]")])],
        family="R_Addition_MultipleBond",
        template=[],
        degeneracy=1.0,
        reversible=True,
        is_forward=False,
    )

    oriented, orientation = _orient_to_proxy(proxy, reaction)

    assert orientation == "reversed"
    assert oriented.reactants is reaction.products
    assert oriented.products is reaction.reactants


def test_implicit_h_count_accepts_real_atoms_without_optional_attribute():
    from rmgpy.molecule.molecule import Molecule

    molecule = Molecule(smiles="CC")
    assert _implicit_h([molecule]) == sum(
        getattr(atom, "implicit_hydrogens", 0) for atom in molecule.atoms
    )


def test_junction_proxy_places_the_capture_radical_on_a_ring_atom():
    proxies = ps_proxy_set(3)
    proxy = next(proxy for proxy in proxies if proxy.site_type == "junction_radical")
    molecule = proxy.reactants[0].molecule[0]
    radical_atoms = [atom for atom in molecule.atoms if atom.radical_electrons]

    assert len(radical_atoms) == 1
    assert molecule.is_atom_in_cycle(radical_atoms[0])
    assert proxy.metadata["family_candidates"] == []
    pair = next(
        proxy for proxy in proxies if proxy.site_type == "junction_radical+end_radical"
    )
    assert pair.metadata["family_candidates"] == ["R_Recombination"]


def test_ortho_reflection_swaps_the_persistent_attacked_atom_id():
    proxy = _ortho_junction_proxy()
    representative = _archived_ortho_reaction(proxy)
    reflected = _reflect_ortho_reaction(representative)
    representative_secondary = _molecule(representative.reactants[0])
    reflected_secondary = _molecule(reflected.reactants[0])
    representative_attacked = next(
        atom for atom in representative_secondary.atoms if atom.radical_electrons
    )
    reflected_attacked = next(
        atom for atom in reflected_secondary.atoms if atom.radical_electrons
    )

    assert representative_attacked.id != reflected_attacked.id
    assert representative_secondary.is_isomorphic(reflected_secondary)
    compiler = EventSetCompiler(
        None,
        [],
        [],
        rmgpy_sha="test-rmgpy",
        rmg_database_sha="test-database",
    )
    record = compiler._record(
        replace(
            proxy,
            metadata={**proxy.metadata, "junction_kind": "J_ortho_S6"},
        ),
        reflected,
        {},
    )
    assert [operation["action"] for operation in record.bond_ops] == [
        "form",
        "set_radical",
        "set_radical",
    ]


def test_symmetry_equivalent_representative_uses_canonical_total_order():
    from rmgpy.molecule.molecule import Molecule

    first = Molecule(smiles="Cc1ccccc1")
    reflected = first.copy(deep=True)
    reflected.atoms.reverse()
    left = SimpleNamespace(molecule=[first, reflected])
    right = SimpleNamespace(molecule=[reflected, first])

    assert _canonical_adjacency(_molecule(left)) == _canonical_adjacency(
        _molecule(right)
    )


def test_opposite_kekule_localisations_collapse_to_one_representative():
    from rmgpy.molecule.molecule import Molecule

    first = Molecule(smiles="c1ccccc1")
    first.kekulize()
    opposite = first.copy(deep=True)
    ring = opposite.get_smallest_set_of_smallest_rings()[0]
    for bond in opposite.get_edges_in_cycle(ring):
        bond.order = 2.0 if bond.order == 1.0 else 1.0
    opposite.update_atomtypes(log_species=False)

    first_graph = _canonical_adjacency(_molecule(SimpleNamespace(molecule=[first])))
    opposite_graph = _canonical_adjacency(
        _molecule(SimpleNamespace(molecule=[opposite]))
    )
    assert first_graph == opposite_graph
    assert ",B}" in first_graph


def test_artifact_is_byte_deterministic_and_content_addressed(tmp_path):
    compiler, database = _compiler()
    first, _ = compiler.write_artifact(tmp_path)
    second, _ = compiler.write_artifact(tmp_path)
    assert first == second
    assert len(database.calls) == 1
    payload = first.read_bytes()
    assert first.stem == hashlib.sha256(payload).hexdigest()

    corrupted = copy.deepcopy(compiler.compile())
    corrupted["records"][0]["atom_map"] = {1: 1, 2: 1}
    with pytest.raises(ValueError, match="bijection"):
        validate_artifact(corrupted)
    assert compiler.compile()["records"][0]["atom_map"] != {1: 1, 2: 1}


def test_artifact_identity_changes_with_each_declared_input(tmp_path):
    compiler, _ = _compiler()
    baseline, _ = compiler.write_artifact(tmp_path / "base")

    molecule = _Molecule()
    database = _Database(_Reaction(molecule))
    changed_family = EventSetCompiler(
        database,
        [SiteProxy("interior", [_Species(molecule)])],
        ["Different_Family"],
        temperature_grid=[600.0, 700.0],
    )
    changed_proxy = EventSetCompiler(
        database,
        [SiteProxy("interior", [_Species(molecule)], metadata={"revision": 2})],
        ["H_Abstraction"],
        temperature_grid=[600.0, 700.0],
    )
    changed_database = EventSetCompiler(
        database,
        [SiteProxy("interior", [_Species(molecule)])],
        ["H_Abstraction"],
        temperature_grid=[600.0, 700.0],
        rmg_database_sha="different-database-sha",
    )

    changed = {
        changed_family.write_artifact(tmp_path / "family")[0].stem,
        changed_proxy.write_artifact(tmp_path / "proxy")[0].stem,
        changed_database.write_artifact(tmp_path / "database")[0].stem,
    }
    assert baseline.stem not in changed
    assert len(changed) == 3


def test_c7_catalogue_size_is_derived_and_corruption_is_rejected():
    catalogue = short_ps_molecule_catalogue(2)
    assert catalogue["max_units"] == 5
    assert catalogue["size"] == len(catalogue["molecules"])
    assert catalogue["termination_bound"] >= catalogue["size"]
    assert any(
        molecule["repeat_units"] == 5
        and molecule["state"] == "closed_shell_linear"
        and molecule["formula"] == {"C": 40, "H": 42}
        for molecule in catalogue["molecules"]
    )

    compiler, _ = _compiler()
    corrupted = compiler.compile()
    corrupted["short_molecule_catalogue"]["size"] += 1
    with pytest.raises(ValueError, match="catalogue size"):
        validate_artifact(corrupted)


def test_c9_ceiling_temperature_uses_compiled_rate_crossing():
    propagation = {
        "k_table": {"T": [500.0, 600.0], "k": [1.0, 100.0]},
    }
    depropagation = {
        "k_table": {"T": [500.0, 600.0], "k": [10.0, 10.0]},
    }
    assert ceiling_temperature(propagation, depropagation, 1.0) == pytest.approx(550.0)

    corrupted = copy.deepcopy(depropagation)
    corrupted["k_table"]["T"] = [500.0, 601.0]
    assert ceiling_temperature(propagation, corrupted, 1.0) is None


def test_ceiling_pair_graph_key_is_canonical_and_structure_sensitive():
    from rmgpy.molecule.molecule import Molecule

    radical = Molecule(smiles="[CH2]CC")
    reordered = radical.copy(deep=True)
    reordered.atoms.reverse()
    alkene = Molecule(smiles="C=CC")

    assert _graph_list_key([radical.to_adjacency_list(remove_h=False)]) == (
        _graph_list_key([reordered.to_adjacency_list(remove_h=False)])
    )
    assert _graph_list_key([radical.to_adjacency_list(remove_h=False)]) != (
        _graph_list_key([alkene.to_adjacency_list(remove_h=False)])
    )


def test_c1_rewrite_application_reproduces_product_and_detects_corruption():
    from rmgpy.kmc.compiler import apply_record
    from rmgpy.molecule.molecule import Molecule

    ethane = Molecule(smiles="CC")
    products = apply_record(
        {
            "atom_map": {0: 0, 1: 1},
            "bond_ops": [
                {"action": "break", "atoms": [0, 1], "order": "1.0"},
                {"action": "set_radical", "atom": 0, "value": 1},
                {"action": "set_radical", "atom": 1, "value": 1},
            ],
        },
        [ethane],
    )
    methyl = Molecule(smiles="[CH3]")
    assert len(products) == 2
    assert all(product.is_isomorphic(methyl) for product in products)

    bad_record = {
        "atom_map": {0: 0, 1: 1},
        "bond_ops": [{"action": "set_radical", "atom": 0, "value": 1}],
    }
    with pytest.raises(ValueError, match="invalid molecular graph"):
        apply_record(bad_record, [ethane])
