"""Fast contract tests for the offline event-set compiler."""

import copy
import hashlib
import os
import subprocess
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.kmc.barrier_e0 import FixedBBarrierE0Provider
import rmgpy.kmc.compiler as compiler_module
from rmgpy.kmc.compiler import (
    DEFAULT_T_GRID,
    EventSetCompiler,
    ARCHIVED_J_PARA_RATE_PROVENANCE,
    PS_FAMILY_FILTER_REASON,
    SiteProxy,
    _archived_ortho_reaction,
    _canonical_adjacency,
    _graph_list_key,
    _implicit_h,
    _molecule,
    _ortho_junction_proxy,
    _orient_to_proxy,
    _proxy_fingerprint,
    _production_padded_root,
    _reflect_ortho_reaction,
    _remap_witness_projection,
    ceiling_temperature,
    compiler_source_hash,
    ps_proxy_set,
    short_ps_molecule_catalogue,
    _validate_thermo_assignment_consistency,
    validate_artifact,
)


def test_proxy_fingerprint_is_independent_of_process_global_atom_ids():
    first = ps_proxy_set(3)[0]
    second = copy.deepcopy(first)
    offset = 100000
    second = replace(
        second,
        artificial_boundaries=tuple(
            replace(boundary, atom_id=boundary.atom_id + offset)
            for boundary in second.artificial_boundaries
        ),
    )
    for participant in second.reactants:
        for molecule in participant.molecule:
            for atom in molecule.atoms:
                atom.id += offset

    assert _proxy_fingerprint(first, include_padding=True) == _proxy_fingerprint(
        second, include_padding=True
    )


def test_padding_is_off_by_default_and_requires_explicit_k():
    instance = EventSetCompiler(None, [], [])
    assert instance.proxy_padding_distance is None
    assert "proxy_boundary_padding" not in instance.compile()["provenance"]


def test_padding_summary_is_published_only_when_enabled():
    disabled = EventSetCompiler(None, [], []).compile()
    enabled = EventSetCompiler(None, [], [], proxy_padding_distance=3).compile()

    assert "proxy_padding_summary" not in disabled
    assert enabled["proxy_padding_summary"] == {
        "padded": 0,
        "not_paddable": 0,
        "excluded": 0,
    }


def test_padded_root_degeneracy_is_derived_by_production_generation():
    from rmgpy.data.kinetics.family import TemplateReaction
    from rmgpy.molecule.molecule import Molecule
    from rmgpy.species import Species

    reactant = Species(molecule=[Molecule(smiles="[CH3]")])
    product = Species(molecule=[Molecule(smiles="[CH3]")])
    expected = TemplateReaction(
        reactants=[reactant], products=[product], family="fake", degeneracy=1.0
    )
    derived = copy.deepcopy(expected)
    derived.degeneracy = 7.0

    class Database:
        def generate_reactions_from_families(self, reactants, **kwargs):
            assert kwargs["products"] and kwargs["only_families"] == ["fake"]
            assert kwargs["resonance"] is False
            return [derived]

    matched, generated_count = _production_padded_root(Database(), expected)
    assert generated_count == 1
    assert matched.degeneracy == 7.0


def test_projection_is_remapped_to_regenerated_participant_order():
    from rmgpy.molecule.molecule import Molecule

    methane = Molecule(smiles="C")
    ethane = Molecule(smiles="CC")
    product = Molecule(smiles="CCC")
    expected = SimpleNamespace(
        reactants=[methane, ethane], products=[product]
    )
    regenerated = SimpleNamespace(
        reactants=[ethane.copy(deep=True), methane.copy(deep=True)],
        products=[product.copy(deep=True)],
    )
    projection = {
        "reactants": [
            {
                "participant_index": 0,
                "witness_participant_index": 0,
                "executable_atom_index": 0,
                "witness_atom_index": 0,
            },
            {
                "participant_index": 1,
                "witness_participant_index": 1,
                "executable_atom_index": 0,
                "witness_atom_index": 0,
            },
            {
                "participant_index": 1,
                "witness_participant_index": 1,
                "executable_atom_index": 1,
                "witness_atom_index": 1,
            },
        ],
        "products": [
            {
                "participant_index": 0,
                "witness_participant_index": 0,
                "executable_atom_index": index,
                "witness_atom_index": index,
            }
            for index in range(3)
        ],
    }

    remapped = _remap_witness_projection(projection, expected, regenerated)

    assert {
        item["witness_participant_index"]
        for item in remapped["reactants"]
        if item["participant_index"] == 0
    } == {1}
    assert {
        item["witness_participant_index"]
        for item in remapped["reactants"]
        if item["participant_index"] == 1
    } == {0}


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


def test_archived_j_para_provenance_keeps_extraction_and_compile_sha_distinct():
    assert ARCHIVED_J_PARA_RATE_PROVENANCE["extraction_database_sha"] == (
        "4a12d36fcdc193ede82c8d1ab5c1653495d445bc"
    )
    assert ARCHIVED_J_PARA_RATE_PROVENANCE["compile_database_sha"] == (
        "cd86d4e1c187a132109e16cd86f624ed9fb217df"
    )
    assert (
        ARCHIVED_J_PARA_RATE_PROVENANCE["extraction_database_sha"]
        != ARCHIVED_J_PARA_RATE_PROVENANCE["compile_database_sha"]
    )


def test_barrier_e0_provider_is_off_by_default_and_explicit_when_enabled():
    disabled, _ = _compiler()
    disabled_artifact = disabled.compile()
    assert "barrier_e0_provider" not in disabled_artifact["provenance"]
    assert "barrier_e0_provider" not in disabled_artifact["inputs"]

    enabled, _ = _compiler()
    enabled.barrier_e0_provider = FixedBBarrierE0Provider(900.0)
    artifact = enabled.compile()
    assert artifact["provenance"]["barrier_e0_provider"] == {
        "enabled": True,
        "name": "fixed-b-wilhoit",
        "version": "1",
        "B_K": 900.0,
        "fit_temperature_grid": "participant ThermoData.Tdata",
        "fit_weights": "uniform least squares",
    }
    assert artifact["inputs"]["barrier_e0_provider"] == artifact["provenance"][
        "barrier_e0_provider"
    ]


def test_disabled_barrier_e0_provider_preserves_base_artifact_bytes(
    monkeypatch, tmp_path
):
    monkeypatch.setattr(
        compiler_module,
        "_LOADED_COMPILER_HASH",
        "17a155a4026a058e2c4933eb58a5cfde79af0274496e0242a4e2b5642d796ac2",
    )
    monkeypatch.setattr(
        compiler_module,
        "compiler_source_hash",
        lambda: "48c1f8b73e311f92d4169734c92b37f7bcc984ee01e2532f0da83bbc573c833d",
    )
    compiler, _ = _compiler()
    compiler.rmgpy_sha = "1e51f1f8f33b81ce66baf3c8b52b84f15ab9f974"

    path, artifact = compiler.write_artifact(tmp_path)

    assert "barrier_e0_provider" not in artifact["provenance"]
    assert "barrier_e0_provider" not in artifact["inputs"]
    # Moves only with PERSISTENT_CARBENE_POLICY_VERSION (v1 reproduces bab79a76...).
    assert path.stem == (
        "bd10d88703512ae17d66d1f543fe6cf4e3186f8f7a45398f0732dc49d33348a4"
    )
    assert hashlib.sha256(path.read_bytes()).hexdigest() == path.stem


def test_default_off_artifact_identity_rejects_unpadded_field_mutant(
    monkeypatch, tmp_path
):
    monkeypatch.setattr(
        compiler_module,
        "_LOADED_COMPILER_HASH",
        "17a155a4026a058e2c4933eb58a5cfde79af0274496e0242a4e2b5642d796ac2",
    )
    monkeypatch.setattr(
        compiler_module,
        "compiler_source_hash",
        lambda: "48c1f8b73e311f92d4169734c92b37f7bcc984ee01e2532f0da83bbc573c833d",
    )
    compiler, _ = _compiler()
    compiler.rmgpy_sha = "1e51f1f8f33b81ce66baf3c8b52b84f15ab9f974"
    compiler.proxies = (
        replace(compiler.proxies[0], site_type="identity-field-mutant"),
    )

    path, _ = compiler.write_artifact(tmp_path)

    with pytest.raises(AssertionError):
        assert path.stem == (
            "bab79a76501a9dd64015caf08e16d236c1c6fb4d3ddedc5fd917dfcd17f75337"
        )


def test_changing_barrier_e0_provider_invalidates_in_memory_artifact():
    compiler, _ = _compiler()
    disabled = compiler.compile()

    compiler.barrier_e0_provider = FixedBBarrierE0Provider(900.0)
    enabled = compiler.compile()

    assert "barrier_e0_provider" not in disabled["provenance"]
    assert "barrier_e0_provider" not in disabled["inputs"]
    assert enabled["provenance"]["barrier_e0_provider"]["B_K"] == 900.0
    assert enabled is not disabled


def test_compiler_source_hash_includes_barrier_e0_provider():
    source_root = Path(compiler_module.__file__).parent
    expected = hashlib.sha256(
        b"".join(
            (source_root / filename).read_bytes()
            for filename in (
                "compiler.py",
                "reference_thermo.py",
                "event_record.py",
                "atom_map.py",
                "kinetics_library.py",
                "database_provenance.py",
                "proxy_padding.py",
                "barrier_e0.py",
            )
        )
    ).hexdigest()
    assert compiler_source_hash() == expected


@pytest.mark.parametrize("kinetics_type", ["ArrheniusBM", "ArrheniusEP"])
def test_compiled_enthalpy_dependent_tree_rate_uses_reaction_enthalpy(
    kinetics_type,
):
    from rmgpy.kinetics.arrhenius import ArrheniusBM, ArrheniusEP

    class TreeFamily:
        auto_generated = True

        def get_kinetics(self, *args, **kwargs):
            entry = SimpleNamespace(rank=7)
            return kinetics, "rate rules", entry, True

    class ThermoDatabase:
        def get_thermo_data(self, species):
            return SimpleNamespace(source="pinned test thermo")

    class ReferenceThermo:
        provenance = {"reference_thermo": "test gas-phase Kc"}

        def equilibrium_constants(self, reaction, temperatures):
            return [1.0 for _ in temperatures]

    class EnthalpyDependentReaction(_Reaction):
        family = "R_Recombination"
        template = ["Root"]
        kinetics = None
        reversible = True
        is_forward = True
        dHrxn298 = -100000.0

        def get_enthalpy_of_reaction(self, temperature):
            assert temperature == 298
            assert all(
                species.thermo is not None
                for species in self.reactants + self.products
            )
            return self.dHrxn298

        def fix_barrier_height(self):
            intrinsic_barrier = self.kinetics.E0.value_si
            self.kinetics = self.kinetics.to_arrhenius(
                self.get_enthalpy_of_reaction(298)
            )
            if (
                self.kinetics.Ea.value_si < 0.0
                and self.kinetics.Ea.value_si < intrinsic_barrier
            ):
                self.kinetics.Ea.value_si = min(0.0, intrinsic_barrier)

    if kinetics_type == "ArrheniusBM":
        kinetics = ArrheniusBM(
            A=(1.0e6, "s^-1"),
            n=0.0,
            w0=(500.0, "kJ/mol"),
            E0=(20.0, "kJ/mol"),
        )
    else:
        kinetics = ArrheniusEP(
            A=(1.0e6, "s^-1"),
            n=0.0,
            alpha=0.5,
            E0=(20.0, "kJ/mol"),
        )
    molecule = _Molecule()
    reaction = EnthalpyDependentReaction(molecule)
    database = _Database(reaction)
    database.families[reaction.family] = TreeFamily()
    temperatures = [600.0, 800.0, 1000.0]
    artifact = EventSetCompiler(
        database,
        [SiteProxy("interior", [_Species(molecule)])],
        [reaction.family],
        temperature_grid=temperatures,
        thermo_database=ThermoDatabase(),
        reference_thermo_provider=ReferenceThermo(),
    ).compile()

    record = next(
        record
        for record in artifact["records"]
        if record["rate_source"]["kind"] == "RMG family estimate"
    )
    model_generation_reaction = copy.deepcopy(reaction)
    model_generation_reaction.kinetics = copy.deepcopy(kinetics)
    for species in (
        model_generation_reaction.reactants + model_generation_reaction.products
    ):
        species.thermo = ThermoDatabase().get_thermo_data(species)
    model_generation_reaction.fix_barrier_height()
    expected = [
        model_generation_reaction.kinetics.get_rate_coefficient(t)
        for t in temperatures
    ]
    zero_enthalpy = [kinetics.get_rate_coefficient(t) for t in temperatures]

    assert record["k_table"]["k"] == pytest.approx(expected, rel=1e-12)
    assert record["k_table"]["k"] != pytest.approx(zero_enthalpy, rel=1e-3)
    if kinetics_type == "ArrheniusBM":
        # Corrected R_Recombination BM pin; the former dHrxn=0 values were
        # [18150.19, 49449.39, 90225.39] on this temperature grid.
        assert record["k_table"]["k"] == pytest.approx(
            [1.0e6, 1.0e6, 1.0e6], rel=1e-12
        )
    conversion = record["rate_source"]["kinetics_conversion"]
    assert {
        key: value
        for key, value in conversion.items()
        if key != "species_thermo_assignments"
    } == {
        "input_model": kinetics_type,
        "method": "to_arrhenius(reaction.get_enthalpy_of_reaction(298))",
        "post_conversion": "reaction.fix_barrier_height()",
        "output_model": "Arrhenius",
        "reaction_enthalpy_J_per_mol": reaction.dHrxn298,
        "activation_energy_J_per_mol": (
            model_generation_reaction.kinetics.Ea.value_si
        ),
    }
    assert len(conversion["species_thermo_assignments"]) == 2
    split_thermo = copy.deepcopy(artifact)
    split_record = next(
        record
        for record in split_thermo["records"]
        if record["rate_source"]["kind"] == "RMG family estimate"
    )
    split_record["rate_source"]["kinetics_conversion"][
        "species_thermo_assignments"
    ][0]["chemical_identity_sha256"] = "deliberately-split-assignment"
    split_record["thermo_provenance"]["species_thermo_assignments"] = copy.deepcopy(
        conversion["species_thermo_assignments"]
    )
    with pytest.raises(ValueError, match="different thermo assignments"):
        _validate_thermo_assignment_consistency(split_thermo["records"])
    reverse = next(
        record
        for record in artifact["records"]
        if record["rate_source"]["kind"] == "reference-thermo reverse"
    )
    assert (
        reverse["rate_source"]["forward_kinetics_conversion"]
        == record["rate_source"]["kinetics_conversion"]
    )


def test_shared_featured_assignment_is_call_order_independent_and_closes_kc():
    from rmgpy.data.kinetics.family import TemplateReaction
    from rmgpy.kinetics.arrhenius import ArrheniusBM
    from rmgpy.molecule.molecule import Molecule
    from rmgpy.species import Species
    from rmgpy.thermo import ThermoData

    class ThermoDatabase:
        library_order = ["primaryThermoLibrary"]

        def get_thermo_data(self, species):
            molecule = species.molecule[0]
            is_h_atom = len(molecule.atoms) == 1 and molecule.atoms[0].is_hydrogen()
            is_radical = molecule.get_radical_count() > 0
            if is_h_atom:
                enthalpy = 218.0
            elif is_radical:
                enthalpy = -80.0 if len(species.molecule) > 1 else 120.0
            else:
                enthalpy = -40.0
            return ThermoData(
                Tdata=([300, 400, 600, 800, 1000], "K"),
                Cpdata=([30, 30, 30, 30, 30], "J/(mol*K)"),
                H298=(enthalpy, "kJ/mol"),
                S298=(100, "J/(mol*K)"),
                E0=(enthalpy, "kJ/mol"),
                comment=(
                    "Thermo library: primaryThermoLibrary"
                    if is_h_atom
                    else "Thermo group additivity estimation: pinned test"
                ),
            )

    def make_reaction(preassigned=False):
        from rmgpy.molecule.molecule import Bond

        toluene = Molecule(smiles="Cc1ccccc1")
        for identifier, atom in enumerate(toluene.atoms, start=1):
            atom.id = identifier
        hydrogen_atom = Molecule(smiles="[H]")
        hydrogen_atom.atoms[0].id = len(toluene.atoms) + 1
        benzyl = toluene.copy(deep=True)
        benzylic_atom = next(
            atom
            for atom in benzyl.atoms
            if atom.element.symbol == "C"
            and not benzyl.is_atom_in_cycle(atom)
        )
        hydrogen = next(
            atom for atom in benzylic_atom.edges if atom.is_hydrogen()
        )
        benzyl.remove_bond(benzyl.get_bond(benzylic_atom, hydrogen))
        benzyl.remove_atom(hydrogen)
        benzylic_atom.radical_electrons = 1
        benzyl.update(sort_atoms=False)
        h2_atoms = [hydrogen_atom.atoms[0].copy(), hydrogen.copy()]
        h2_atoms[0].radical_electrons = 0
        h2_atoms[0].charge = 0
        h2 = Molecule(atoms=h2_atoms, multiplicity=1)
        h2.add_bond(Bond(h2_atoms[0], h2_atoms[1], order=1))
        h2.update(sort_atoms=False)
        reaction = TemplateReaction(
            reactants=[
                Species(label="H", molecule=[hydrogen_atom]),
                Species(label="toluene", molecule=[toluene]),
            ],
            products=[
                Species(label="H2", molecule=[h2]),
                Species(label="benzyl", molecule=[benzyl]),
            ],
            kinetics=ArrheniusBM(
                A=(1.0e6, "m^3/(mol*s)"),
                n=0.0,
                w0=(500.0, "kJ/mol"),
                E0=(20.0, "kJ/mol"),
            ),
            reversible=True,
            degeneracy=1.0,
            family="H_Abstraction",
            template=["featured-root"],
            is_forward=True,
        )
        if preassigned:
            reaction.products[1].thermo = ThermoData(
                Tdata=([300, 400, 600, 800, 1000], "K"),
                Cpdata=([30, 30, 30, 30, 30], "J/(mol*K)"),
                H298=(999.0, "kJ/mol"),
                S298=(100, "J/(mol*K)"),
                comment="stale call-order thermo",
            )
        return reaction

    def compile_reaction(reaction, *, kc_first=False):
        original_graphs = [
            species.molecule[0].to_adjacency_list(remove_h=False)
            for species in reaction.reactants + reaction.products
        ]
        kinetics_database = SimpleNamespace(
            families={"H_Abstraction": SimpleNamespace(auto_generated=True)}
        )
        compiler = EventSetCompiler(
            kinetics_database,
            [SiteProxy("featured", reaction.reactants)],
            ["H_Abstraction"],
            thermo_database=ThermoDatabase(),
            rmg_database_sha="test-only-database-commit",
            reaction_cache={"featured": [reaction]},
        )
        if kc_first:
            compiler.reference_thermo_provider.evaluate(
                copy.deepcopy(reaction), DEFAULT_T_GRID
            )
        artifact = compiler.compile()
        assert original_graphs == [
            species.molecule[0].to_adjacency_list(remove_h=False)
            for species in reaction.reactants + reaction.products
        ]
        forward = next(
            item
            for item in artifact["records"]
            if item["rate_source"]["kind"] == "RMG family estimate"
        )
        reverse_records = [
            item
            for item in artifact["records"]
            if item["rate_source"]["kind"] == "reference-thermo reverse"
        ]
        assert reverse_records, [
            (item["status"], item["status_reason"], item["rate_source"])
            for item in artifact["records"]
        ]
        reverse = reverse_records[0]
        return forward, reverse

    clean_forward, clean_reverse = compile_reaction(make_reaction())
    kc_first_forward, kc_first_reverse = compile_reaction(
        make_reaction(), kc_first=True
    )
    stale_forward, stale_reverse = compile_reaction(make_reaction(preassigned=True))
    assert kc_first_forward["k_table"] == clean_forward["k_table"]
    assert kc_first_reverse["k_table"] == clean_reverse["k_table"]
    assert stale_forward["k_table"] == clean_forward["k_table"]
    assert stale_reverse["k_table"] == clean_reverse["k_table"]
    assert clean_forward["k_table"]["T"] == list(DEFAULT_T_GRID)

    conversion = clean_forward["rate_source"]["kinetics_conversion"]
    thermo = clean_forward["thermo_provenance"]
    assert conversion["reaction_enthalpy_J_per_mol"] == pytest.approx(
        thermo["reaction_enthalpy_298_J_per_mol"]
    )
    for forward_rate, reverse_rate, equilibrium_constant in zip(
        clean_forward["k_table"]["k"],
        clean_reverse["k_table"]["k"],
        thermo["equilibrium_constant_table"]["Kc"],
    ):
        assert forward_rate / reverse_rate == pytest.approx(
            equilibrium_constant, rel=1e-12
        )

    old_direct = make_reaction()
    old_database = ThermoDatabase()
    for species in old_direct.reactants + old_direct.products:
        species.thermo = old_database.get_thermo_data(species)
    assert abs(
        conversion["reaction_enthalpy_J_per_mol"]
        - old_direct.get_enthalpy_of_reaction(298)
    ) > 40000.0
    assert any(
        item["thermo_source"] == {
            "kind": "library",
            "library": "primaryThermoLibrary",
        }
        for item in conversion["species_thermo_assignments"]
    )


def test_compiled_real_bm_tree_node_matches_rmg_model_generation():
    """A real family tree node must not inherit the BM dHrxn=0 default."""
    from rmgpy.kinetics.arrhenius import ArrheniusBM
    from rmgpy.species import Species

    database_path = os.path.join(
        settings["test_data.directory"], "testing_database"
    )
    database = RMGDatabase()
    database.load_kinetics(
        os.path.join(database_path, "kinetics"),
        reaction_libraries=[],
        seed_mechanisms=None,
        kinetics_families=["Disproportionation"],
        kinetics_depositories=["training"],
    )
    database.load_thermo(
        os.path.join(database_path, "thermo"),
        thermo_libraries=["primaryThermoLibrary"],
        depository=True,
    )
    family = database.kinetics.families["Disproportionation"]
    # Force the normal tree estimate instead of this test database's exact
    # ethyl-disproportionation training match.
    family.depositories = []
    reactants = [Species(label="ethyl", smiles="C[CH2]") for _ in range(2)]
    generated = database.kinetics.generate_reactions_from_families(
        reactants,
        only_families=["Disproportionation"],
        resonance=True,
    )
    assert len(generated) == 1
    reaction = generated[0]
    kinetics, source, entry, estimated_forward = family.get_kinetics(
        reaction,
        template_labels=reaction.template,
        degeneracy=reaction.degeneracy,
        return_all_kinetics=False,
    )
    assert isinstance(kinetics, ArrheniusBM)
    assert source == "rate rules"
    assert entry.label == reaction.template[0]
    assert estimated_forward

    for species in reaction.reactants + reaction.products:
        species.thermo = database.thermo.get_thermo_data(species)
    reaction_enthalpy = reaction.get_enthalpy_of_reaction(298)
    temperatures = [600.0, 800.0, 1000.0]
    expected_kinetics = kinetics.to_arrhenius(reaction_enthalpy)
    expected_rates = [
        expected_kinetics.get_rate_coefficient(temperature)
        for temperature in temperatures
    ]
    model_reaction = copy.deepcopy(reaction)
    model_reaction.kinetics = copy.deepcopy(kinetics)
    model_reaction.fix_barrier_height()
    assert [
        model_reaction.kinetics.get_rate_coefficient(temperature)
        for temperature in temperatures
    ] == pytest.approx(expected_rates, rel=1e-12)

    artifact = EventSetCompiler(
        database.kinetics,
        [SiteProxy("real_disproportionation", reactants)],
        ["Disproportionation"],
        temperature_grid=temperatures,
        thermo_database=database.thermo,
        reaction_cache={"real_disproportionation": generated},
    ).compile()
    forward = next(
        record
        for record in artifact["records"]
        if record["rate_source"]["kind"] == "RMG family estimate"
    )
    reverse = next(
        record
        for record in artifact["records"]
        if record["rate_source"]["kind"] == "reference-thermo reverse"
    )

    # Re-pinned from the invalid dHrxn=0 path to the real reaction enthalpy.
    assert forward["k_table"]["k"] == pytest.approx(
        [7752417.767623479, 6930886.275599589, 6354101.468792986],
        rel=1e-12,
    )
    assert forward["k_table"]["k"] == pytest.approx(expected_rates, rel=1e-12)
    assert forward["k_table"]["k"] != pytest.approx(
        [kinetics.get_rate_coefficient(t) for t in temperatures], rel=1e-3
    )
    conversion = forward["rate_source"]["kinetics_conversion"]
    assert conversion["input_model"] == "ArrheniusBM"
    assert conversion["post_conversion"] == "reaction.fix_barrier_height()"
    assert conversion["activation_energy_J_per_mol"] == pytest.approx(0.0)
    assert conversion["reaction_enthalpy_J_per_mol"] == pytest.approx(
        -272269.616, abs=1e-6
    )
    assert reverse["rate_source"]["forward_kinetics_conversion"] == conversion


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
        (["H_Abstraction"], True),
        (["No_Firing_Path"], True),
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
    assert proxy.metadata["family_candidates"] == list(
        compiler_module.PS_FAMILY_CANDIDATES
    )
    pair = next(
        proxy for proxy in proxies if proxy.site_type == "junction_radical+end_radical"
    )
    assert pair.metadata["family_candidates"] == list(
        compiler_module.PS_FAMILY_CANDIDATES
    )


def test_pair_graph_cache_is_scoped_and_cleared_after_failure(monkeypatch):
    from rmgpy.molecule.molecule import Molecule
    from rmgpy.species import Species

    species = Species(molecule=[Molecule(smiles="[CH3]")])
    compiler = EventSetCompiler(None, [], [])
    selected = []
    ordered = []
    original_order = compiler_module._molecule_total_order

    def observed_order(molecule):
        ordered.append(molecule)
        return original_order(molecule)

    monkeypatch.setattr(compiler_module, "_molecule_total_order", observed_order)

    def fail_pair(*args):
        selected.append(compiler_module._molecule(species))
        assert compiler_module._molecule(species) is selected[0]
        raise ValueError("induced pair failure")

    monkeypatch.setattr(compiler, "_build_linked_family_pair", fail_pair)
    with pytest.raises(ValueError, match="induced pair failure"):
        compiler._linked_family_pair(None, None, {})
    assert compiler_module._PAIR_MOLECULE_CACHE.get() is None
    assert len(ordered) == 1
    compiler_module._molecule(species)
    assert len(ordered) == 2


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
