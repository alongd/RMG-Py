"""Identities, dependency manifests, snapshot binding and the request-binding envelope."""

import copy
import io
import json
import os
import subprocess
import sys
import textwrap
from pathlib import Path
from types import SimpleNamespace

import pytest

from rmgpy.kmc import compiler, generation_identity as gi
from rmgpy.kmc.compile import CompileRequest
from rmgpy.molecule.molecule import Molecule
from rmgpy.species import Species

DATABASE = Path(os.environ.get("RMG_DATABASE_PATH", "/home/alon/runs/phase2b/database"))
FAMILIES = ["H_Abstraction", "R_Recombination", "intra_H_migration"]
ALL_FAMILIES = ["H_Abstraction", "R_Recombination", "intra_H_migration", "Disproportionation"]


@pytest.fixture(autouse=True)
def stable_seed(monkeypatch):
    monkeypatch.setenv("PYTHONHASHSEED", "0")
    monkeypatch.delenv("RMG_KMC_PLPSEC_LIBRARY", raising=False)
    monkeypatch.delenv("RMG_DATABASE_SHA", raising=False)
    gi.reset_caches()


@pytest.fixture(scope="module")
def kinetics():
    if not (DATABASE / "input/kinetics/families").is_dir():
        pytest.skip("pinned RMG database is not available")
    from rmgpy.data.rmg import RMGDatabase

    database = RMGDatabase()
    database.load_kinetics(
        str(DATABASE / "input/kinetics"),
        reaction_libraries=[],
        seed_mechanisms=None,
        kinetics_families=ALL_FAMILIES,
        kinetics_depositories=["training"],
    )
    return database.kinetics


def species(smiles, label="", index=-1):
    item = Species(molecule=[Molecule(smiles=smiles)], label=label, index=index)
    for number, atom in enumerate(item.molecule[0].atoms, 1):
        atom.id = number
    return item


def request_sha(kinetics_database, reactants, **options):
    return gi.generation_request_identity(kinetics_database, reactants, **options)


def pair():
    return [species("C[CH]c1ccccc1"), species("CC(C)c1ccccc1")]


# ---------------------------------------------------------------------------
# Canonical encoding
# ---------------------------------------------------------------------------


def test_canonical_bytes_are_key_order_independent_and_type_exact():
    assert gi.canonical_bytes({"b": 1, "a": [1, 2]}) == gi.canonical_bytes({"a": [1, 2], "b": 1})
    assert gi.canonical_bytes(True) != gi.canonical_bytes(1)
    assert gi.canonical_bytes(1) != gi.canonical_bytes(1.0)
    assert gi.canonical_bytes((1, 2)) == gi.canonical_bytes([1, 2])


def test_floats_are_encoded_exactly_through_float_hex():
    assert float.hex(0.1).encode() in gi.canonical_bytes(0.1)
    assert gi.canonical_bytes(0.1) != gi.canonical_bytes(0.1 + 2**-56)
    assert gi.canonical_bytes(0.0) != gi.canonical_bytes(-0.0)
    for bad in (float("nan"), float("inf")):
        with pytest.raises(gi.IdentityError):
            gi.canonical_bytes(bad)


def test_sha_is_full_length_and_domain_separated():
    a = gi.sha256_canonical("generation-request", {"x": 1})
    b = gi.sha256_canonical("compilation-request", {"x": 1})
    assert len(a) == 64 and a != b
    assert a == gi.sha256_canonical("generation-request", {"x": 1})


def test_non_json_values_are_refused():
    with pytest.raises(gi.IdentityError):
        gi.canonical_bytes({"x": object()})
    with pytest.raises(gi.IdentityError):
        gi.canonical_bytes({1: "non-string key"})


# ---------------------------------------------------------------------------
# Hash seed
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("value", [None, "random", "", "abc", "-1", "4294967296"])
def test_unstable_hash_seed_refuses_every_identity(monkeypatch, kinetics, value):
    if value is None:
        monkeypatch.delenv("PYTHONHASHSEED", raising=False)
    else:
        monkeypatch.setenv("PYTHONHASHSEED", value)
    gi.reset_caches()
    with pytest.raises(gi.UnstableHashSeedError):
        request_sha(kinetics, pair())
    with pytest.raises(gi.UnstableHashSeedError):
        gi.generation_dependency_manifest()
    with pytest.raises(gi.UnstableHashSeedError):
        gi.compilation_dependency_manifest()
    with pytest.raises(gi.UnstableHashSeedError):
        gi.compilation_request_identity(CompileRequest("/x", "/y", ("a",)))


@pytest.mark.parametrize("value", ["0", "1", "12345", "4294967295"])
def test_numeric_hash_seed_is_accepted(monkeypatch, value):
    monkeypatch.setenv("PYTHONHASHSEED", value)
    assert gi.require_stable_hash_seed() == value


# ---------------------------------------------------------------------------
# Concrete generation request
# ---------------------------------------------------------------------------


def test_generation_identity_is_deterministic_and_does_not_mutate_inputs(kinetics):
    reactants = pair()
    before = [m.to_adjacency_list() for s in reactants for m in s.molecule]
    first = request_sha(kinetics, reactants)
    assert first == request_sha(kinetics, reactants)
    assert first == request_sha(kinetics, pair())
    assert before == [m.to_adjacency_list() for s in reactants for m in s.molecule]
    assert len(first) == 64


def test_product_filters_hash_graphs_not_labels(kinetics):
    one = species("CCc1ccccc1", label="target")
    two = species("CC(C)c1ccccc1", label="target")
    assert one.label == two.label and not one.molecule[0].is_isomorphic(two.molecule[0])
    assert request_sha(kinetics, pair(), products=[one]) != request_sha(kinetics, pair(), products=[two])


def test_product_filter_none_differs_from_empty_and_from_populated(kinetics):
    shas = {
        request_sha(kinetics, pair(), products=None),
        request_sha(kinetics, pair(), products=[]),
        request_sha(kinetics, pair(), products=[species("CCc1ccccc1")]),
    }
    assert len(shas) == 3


def test_effective_family_list_is_resolved_never_none(kinetics):
    everything = list(kinetics.families)
    assert request_sha(kinetics, pair(), only_families=None) == request_sha(
        kinetics, pair(), only_families=everything
    )
    assert request_sha(kinetics, pair(), only_families=None) == request_sha(
        kinetics, pair(), only_families=list(reversed(everything))
    )
    assert request_sha(kinetics, pair(), only_families=["H_Abstraction"]) != request_sha(
        kinetics, pair(), only_families=None
    )
    assert request_sha(kinetics, pair(), only_families=["H_Abstraction"]) != request_sha(
        kinetics, pair(), only_families=["R_Recombination"]
    )


def test_every_call_argument_matters(kinetics):
    base = request_sha(kinetics, pair(), only_families=["H_Abstraction"])
    assert base != request_sha(kinetics, pair(), only_families=["H_Abstraction"], resonance=False)
    assert base != request_sha(kinetics, pair()[:1], only_families=["H_Abstraction"])
    assert base != request_sha(kinetics, list(reversed(pair())), only_families=["H_Abstraction"])


def mutated(edit):
    reactants = pair()
    edit(reactants)
    return reactants


@pytest.mark.parametrize(
    "name, edit",
    {
        "reactive flag": lambda r: setattr(r[0].molecule[0], "reactive", False),
        "atom id": lambda r: setattr(r[0].molecule[0].atoms[3], "id", 999),
        "atom label": lambda r: setattr(r[0].molecule[0].atoms[0], "label", "*1"),
        "atom order": lambda r: r[0].molecule[0].atoms.reverse(),
        "charge": lambda r: setattr(r[0].molecule[0].atoms[0], "charge", 1),
        "radical electrons": lambda r: setattr(r[0].molecule[0].atoms[0], "radical_electrons", 2),
        "lone pairs": lambda r: setattr(r[0].molecule[0].atoms[0], "lone_pairs", 1),
        "isotope": lambda r: setattr(
            r[0].molecule[0].atoms[0], "element", __import__("rmgpy.molecule.element", fromlist=["x"]).get_element("C", 13)
        ),
        "bond order": lambda r: setattr(next(iter(r[0].molecule[0].atoms[0].edges.values())), "order", 2.0),
        "multiplicity": lambda r: setattr(r[0].molecule[0], "multiplicity", 4),
        "species label": lambda r: setattr(r[0], "label", "named"),
        "species index": lambda r: setattr(r[0], "index", 7),
        "atom props": lambda r: r[0].molecule[0].atoms[0].props.update({"extra": True}),
        "added resonance structure": lambda r: (
            r[0].molecule.append(Molecule(smiles="C[C](C)c1ccccc1")),
        ),
        "resonance order": lambda r: (
            r[0].molecule.append(Molecule(smiles="C[CH]c1ccccc1")),
            r[0].molecule.reverse(),
        ),
    }.items(),
)
def test_observable_input_state_changes_the_identity(kinetics, name, edit):
    baseline = request_sha(kinetics, pair())
    assert request_sha(kinetics, mutated(edit)) != baseline, name


def test_resonance_order_is_observable():
    first = species("C[CH]c1ccccc1")
    first.generate_resonance_structures()
    assert len(first.molecule) > 2
    second = copy.deepcopy(first)
    second.molecule.reverse()
    assert gi.participant_state(first) != gi.participant_state(second)


def test_effective_global_generation_policy_is_hashed(kinetics, monkeypatch):
    from rmgpy.data.kinetics import family as family_module
    from rmgpy.rmg import input as rmg_input

    baseline = request_sha(kinetics, pair())
    monkeypatch.setattr(family_module, "TOLERATE_MISSING_REVERSE", True)
    tolerant = request_sha(kinetics, pair())
    assert tolerant != baseline
    monkeypatch.setattr(family_module, "TOLERATE_MISSING_REVERSE", 5)
    assert request_sha(kinetics, pair()) not in (baseline, tolerant)
    monkeypatch.setattr(family_module, "TOLERATE_MISSING_REVERSE", False)
    assert request_sha(kinetics, pair()) == baseline
    monkeypatch.setattr(family_module, "DROPPED_REVERSE_LOG_PATH", "/tmp/somewhere.jsonl")
    assert request_sha(kinetics, pair()) != baseline
    monkeypatch.setattr(family_module, "DROPPED_REVERSE_LOG_PATH", None)

    monkeypatch.setattr(
        rmg_input, "rmg",
        SimpleNamespace(species_constraints={"maximumCarbonAtoms": 12}, polymer_constraints=None),
        raising=False,
    )
    capped = request_sha(kinetics, pair())
    assert capped != baseline
    monkeypatch.setattr(
        rmg_input, "rmg",
        SimpleNamespace(species_constraints={"maximumCarbonAtoms": 13}, polymer_constraints=None),
        raising=False,
    )
    assert request_sha(kinetics, pair()) != capped
    monkeypatch.setattr(
        rmg_input, "rmg",
        SimpleNamespace(
            species_constraints={"maximumCarbonAtoms": 13}, polymer_constraints={"polymerSizeThreshold": 9}
        ),
        raising=False,
    )
    assert request_sha(kinetics, pair()) not in (baseline, capped)


def test_generation_policy_description_lists_the_global_switches(kinetics):
    policy = gi.generation_policy()
    assert {"tolerate_missing_reverse", "dropped_reverse_log_path", "species_constraints",
            "polymer_constraints"} <= set(policy)


def edit_group(family):
    entry = next(iter(family.groups.entries.values()))
    entry.short_desc = entry.short_desc + " edited"


def edit_group_graph(family):
    entry = next(e for e in family.groups.entries.values() if hasattr(e.item, "atoms") and e.item.atoms)
    entry.item.atoms[0].label = "*99"


def edit_rule(family):
    entry = next(e for rules in family.rules.entries.values() for e in rules)
    entry.rank = (entry.rank or 0) + 1


def edit_recipe(family):
    family.forward_recipe.actions.append(["CHANGE_BOND", "*1", 0, "*2"])


def edit_forbidden(family):
    from rmgpy.data.base import ForbiddenStructures

    family.forbidden = family.forbidden or ForbiddenStructures()
    family.forbidden.load_entry(label="fixture-forbidden", group="""
1 *1 R u0 {2,S}
2 *2 R u0 {1,S}
""")


@pytest.mark.parametrize("editor", [edit_group, edit_group_graph, edit_rule, edit_recipe, edit_forbidden])
def test_editing_a_selected_family_changes_the_key_but_an_unselected_one_does_not(kinetics, editor):
    selected = request_sha(kinetics, pair(), only_families=["H_Abstraction"])
    unselected_other = request_sha(kinetics, pair(), only_families=["R_Recombination"])
    family = kinetics.families["H_Abstraction"]
    backup = copy.deepcopy(family)
    try:
        editor(family)
        assert request_sha(kinetics, pair(), only_families=["H_Abstraction"]) != selected
        assert request_sha(kinetics, pair(), only_families=["R_Recombination"]) == unselected_other
        # a request that selects every family sees the edit as well
        assert gi.family_content_fingerprint(family) != gi.family_content_fingerprint(backup)
    finally:
        kinetics.families["H_Abstraction"] = backup
    assert request_sha(kinetics, pair(), only_families=["H_Abstraction"]) == selected


def test_family_fingerprint_is_stable_and_distinguishes_families(kinetics):
    a = gi.family_content_fingerprint(kinetics.families["H_Abstraction"])
    assert a == gi.family_content_fingerprint(kinetics.families["H_Abstraction"])
    assert a != gi.family_content_fingerprint(kinetics.families["R_Recombination"])
    assert len(a) == 64


def test_generation_identity_binds_the_dependency_manifest(kinetics, monkeypatch):
    baseline = request_sha(kinetics, pair())
    real = gi.generation_dependency_manifest()
    forged = copy.deepcopy(real)
    forged["python"]["version"] = "forged"
    monkeypatch.setattr(gi, "generation_dependency_manifest", lambda: forged)
    assert request_sha(kinetics, pair()) != baseline


# ---------------------------------------------------------------------------
# Dependency manifests
# ---------------------------------------------------------------------------


def test_manifests_pin_interpreter_libraries_and_module_bytes():
    manifest = gi.generation_dependency_manifest()
    assert manifest["python"]["version"] == sys.version
    assert manifest["python"]["implementation"]
    assert manifest["libraries"]["numpy"]["version"]
    assert manifest["libraries"]["scipy"]["version"]
    assert manifest["libraries"]["rdkit"]["version"] and "build" in manifest["libraries"]["rdkit"]
    assert manifest["hash_seed"] == "0"
    modules = manifest["modules"]
    assert "rmgpy.data.kinetics.family" in modules and "rmgpy.molecule.molecule" in modules
    extension = modules["rmgpy.molecule.molecule"]
    assert any(name.endswith(".so") for name in extension), "extension module bytes are pinned"
    assert all(len(digest) == 64 for files in modules.values() for digest in files.values())
    assert "git" in manifest


def test_compilation_manifest_extends_the_generation_manifest():
    generation = gi.generation_dependency_manifest()
    compilation = gi.compilation_dependency_manifest()
    assert set(generation["modules"]) < set(compilation["modules"])
    assert "rmgpy.kmc.compiler" in compilation["modules"] and "rmgpy.kmc.compiler" not in generation["modules"]
    assert gi.manifest_sha(generation) != gi.manifest_sha(compilation)


def test_every_named_module_resolves_to_a_file():
    for name in (*gi.GENERATION_MODULES, *gi.COMPILATION_EXTRA_MODULES):
        assert gi.module_files(name), name


def test_module_manifest_changes_with_source_bytes(tmp_path, monkeypatch):
    package = tmp_path / "fixturepkg_i107"
    package.mkdir()
    (package / "__init__.py").write_text("")
    (package / "leaf.py").write_text("X = 1\n")
    monkeypatch.syspath_prepend(str(tmp_path))
    before = gi.module_manifest(("fixturepkg_i107.leaf",))
    assert gi.module_manifest(("fixturepkg_i107.leaf",)) == before
    (package / "leaf.py").write_text("X = 2\n")
    assert gi.module_manifest(("fixturepkg_i107.leaf",)) != before
    with pytest.raises(gi.IdentityError):
        gi.module_manifest(("fixturepkg_i107.missing",))


def test_manifest_is_located_from_the_package_not_the_working_directory(tmp_path):
    code = textwrap.dedent(
        """
        import json, os
        os.chdir(%r)
        from rmgpy.kmc import generation_identity as gi
        print(json.dumps([gi.manifest_sha(gi.generation_dependency_manifest()),
                          gi.manifest_sha(gi.compilation_dependency_manifest())]))
        """
    ) % str(tmp_path)
    environment = {**os.environ, "PYTHONHASHSEED": "0",
                   "PYTHONPATH": os.pathsep.join(filter(None, [str(Path(compiler.__file__).parents[2]),
                                                                os.environ.get("PYTHONPATH", "")]))}
    output = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, cwd=tmp_path, env=environment, check=True
    ).stdout.strip().splitlines()[-1]
    here = [gi.manifest_sha(gi.generation_dependency_manifest()),
            gi.manifest_sha(gi.compilation_dependency_manifest())]
    assert json.loads(output) == here


def test_git_state_is_absent_outside_a_checkout_and_present_inside(tmp_path):
    assert gi.git_state(tmp_path) == {"sha": None, "dirty": None}
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    subprocess.run(["git", "-c", "user.email=a@b", "-c", "user.name=n", "commit", "-q", "--allow-empty", "-m", "x"],
                   cwd=tmp_path, check=True)
    state = gi.git_state(tmp_path)
    assert len(state["sha"]) == 40 and state["dirty"] is False
    (tmp_path / "tracked.txt").write_text("a")
    subprocess.run(["git", "add", "tracked.txt"], cwd=tmp_path, check=True)
    assert gi.git_state(tmp_path)["dirty"] is True


def test_unpinnable_items_are_documented():
    assert "native libraries" in gi.UNPINNED_ITEMS.lower() or "native" in gi.UNPINNED_ITEMS.lower()
    assert "drift" in gi.UNPINNED_ITEMS.lower()


# ---------------------------------------------------------------------------
# Compilation request
# ---------------------------------------------------------------------------


def make_database(root: Path, extra: str = "") -> Path:
    (root / "input/kinetics/families/H_Abstraction").mkdir(parents=True)
    (root / "input/kinetics/families/H_Abstraction/groups.py").write_text("name = 'h' " + extra + "\n")
    (root / "input/kinetics/families/Other").mkdir(parents=True)
    (root / "input/kinetics/families/Other/groups.py").write_text("name = 'o'\n")
    (root / "input/thermo").mkdir(parents=True)
    (root / "input/thermo/library.py").write_text("x = 1\n")
    return root


def compile_request(database, output="/out/a", **options):
    return CompileRequest(
        database_path=str(database), output_dir=str(output), families=("H_Abstraction",), **options
    )


def compile_sha(request):
    return gi.compilation_request_identity(request)


def test_output_dir_and_database_location_do_not_enter_the_compilation_identity(tmp_path):
    first = make_database(tmp_path / "a")
    second = make_database(tmp_path / "b")
    base = compile_sha(compile_request(first, output="/out/one"))
    assert compile_sha(compile_request(first, output="/somewhere/else")) == base
    assert compile_sha(compile_request(second, output="/out/one")) == base
    description = gi.compilation_request_description(compile_request(first))
    assert "output_dir" not in json.dumps(description) and str(first) not in json.dumps(description)


def test_defaults_are_resolved_before_hashing(tmp_path, monkeypatch):
    database = make_database(tmp_path / "db")
    explicit = compile_request(
        database,
        proxy_units=compiler.PS_PROXY_UNITS,
        family_universe=("H_Abstraction", "Other"),
        database_sha="f" * 40,
        temperature_grid=tuple(compiler.DEFAULT_T_GRID),
    )
    monkeypatch.setenv("RMG_DATABASE_SHA", "f" * 40)
    implicit = compile_request(database)
    assert compile_sha(explicit) == compile_sha(implicit)
    monkeypatch.setenv("RMG_DATABASE_SHA", "e" * 40)
    assert compile_sha(explicit) != compile_sha(compile_request(database))
    assert compile_sha(compile_request(database, proxy_units=2)) != compile_sha(implicit)
    assert compile_sha(compile_request(database, family_universe=("H_Abstraction",))) != compile_sha(implicit)


@pytest.mark.parametrize(
    "options",
    [
        {"barrier_e0_fixed_b": 0.5},
        {"temperature_grid": (600.0, 700.0)},
        {"proxy_units": 5},
        {"database_sha": "a" * 40},
    ],
)
def test_every_request_field_is_hashed(tmp_path, options):
    database = make_database(tmp_path / "db")
    assert compile_sha(compile_request(database)) != compile_sha(compile_request(database, **options))


def test_request_families_are_hashed(tmp_path):
    database = make_database(tmp_path / "db")
    other = CompileRequest(str(database), "/out", families=("H_Abstraction", "Other"))
    assert compile_sha(compile_request(database)) != compile_sha(other)


def test_compile_environment_switches_are_hashed(tmp_path, monkeypatch):
    database = make_database(tmp_path / "db")
    base = compile_sha(compile_request(database))
    monkeypatch.setenv("RMG_KMC_PLPSEC_LIBRARY", "1")
    assert compile_sha(compile_request(database)) == base  # "1" is the default
    monkeypatch.setenv("RMG_KMC_PLPSEC_LIBRARY", "0")
    assert compile_sha(compile_request(database)) != base
    monkeypatch.setenv("RMG_KMC_PLPSEC_LIBRARY", "2")
    with pytest.raises(ValueError):
        compile_sha(compile_request(database))
    monkeypatch.setenv("RMG_KMC_PLPSEC_LIBRARY", "0")
    assert gi.compile_environment_switches() == {"RMG_KMC_PLPSEC_LIBRARY": "0"}


def test_database_content_edits_additions_and_removals_change_the_identity(tmp_path):
    database = make_database(tmp_path / "db")
    base = compile_sha(compile_request(database))
    groups = database / "input/kinetics/families/H_Abstraction/groups.py"
    original = groups.read_text()
    groups.write_text(original + "# edit\n")
    assert compile_sha(compile_request(database)) != base
    groups.write_text(original)
    assert compile_sha(compile_request(database)) == base
    added = database / "input/thermo/extra.py"
    added.write_text("y = 2\n")
    assert compile_sha(compile_request(database)) != base
    added.unlink()
    assert compile_sha(compile_request(database)) == base


def git(database, *arguments):
    subprocess.run(
        ["git", "-c", "user.email=a@b", "-c", "user.name=n", *arguments], cwd=database, check=True,
        capture_output=True,
    )


def test_database_git_identity_fields_embedded_in_artifacts_are_hashed(tmp_path):
    database = make_database(tmp_path / "db")
    base = compile_sha(compile_request(database))
    git(database, "init", "-q")
    git(database, "add", "-A")
    git(database, "commit", "-q", "-m", "one")
    clean = compile_sha(compile_request(database))
    assert clean != base  # git presence, revision and cleanliness are embedded fields
    (database / "input/thermo/library.py").write_text("x = 1\n")  # identical bytes: still clean
    assert compile_sha(compile_request(database)) == clean
    git(database, "commit", "-q", "--allow-empty", "-m", "two")
    revised = compile_sha(compile_request(database))
    assert revised != clean  # revision changed, content did not
    (database / "untracked.txt").write_text("dirty")
    assert compile_sha(compile_request(database)) != revised  # cleanliness changed, content did not


def test_producer_metadata_embedded_in_artifacts_is_hashed(tmp_path, monkeypatch):
    database = make_database(tmp_path / "db")
    base = compile_sha(compile_request(database))
    with monkeypatch.context() as patch:
        patch.setattr(compiler, "_LOADED_SOURCE_HASH", "0" * 64)
        assert compile_sha(compile_request(database)) != base
    with monkeypatch.context() as patch:
        patch.setattr(compiler, "_LOADED_COMPILER_HASH", "0" * 64)
        assert compile_sha(compile_request(database)) != base
    assert compile_sha(compile_request(database)) == base
    with monkeypatch.context() as patch:
        patch.setattr(compiler, "_git_sha", lambda path: "d" * 40)
        assert compile_sha(compile_request(database)) != base


def test_compilation_identity_binds_the_compilation_manifest(tmp_path, monkeypatch):
    database = make_database(tmp_path / "db")
    base = compile_sha(compile_request(database))
    forged = copy.deepcopy(gi.compilation_dependency_manifest())
    forged["libraries"]["numpy"]["version"] = "0.0"
    monkeypatch.setattr(gi, "compilation_dependency_manifest", lambda: forged)
    assert compile_sha(compile_request(database)) != base


def test_published_content_sha_is_the_artifact_id():
    artifact = {"b": [1, 2.5, None], "a": {"z": "é"}}
    buffer = io.BytesIO()
    assert compiler._write_canonical_json(artifact, buffer) == gi.published_content_sha(artifact)
    assert gi.published_content_sha(artifact) == compiler.sha256_json(artifact)


# ---------------------------------------------------------------------------
# Snapshot binding
# ---------------------------------------------------------------------------


def test_snapshot_manifest_lists_the_input_files_and_verifies(tmp_path):
    database = make_database(tmp_path / "db")
    manifest = gi.database_snapshot_manifest(database)
    assert set(manifest["files"]) == {
        "input/kinetics/families/H_Abstraction/groups.py",
        "input/kinetics/families/Other/groups.py",
        "input/thermo/library.py",
    }
    assert all(len(value) == 64 for value in manifest["files"].values())
    gi.verify_database_snapshot(manifest, database)
    relocated = tmp_path / "relocated"
    subprocess.run(["cp", "-r", str(database), str(relocated)], check=True)
    gi.verify_database_snapshot(manifest, relocated)


@pytest.mark.parametrize(
    "tamper",
    ["edit", "replace", "replace_with_symlink", "add_kinetics", "add_thermo", "remove", "swap_directory"],
)
def test_snapshot_verification_detects_every_tamper(tmp_path, tamper):
    database = make_database(tmp_path / "db")
    manifest = gi.database_snapshot_manifest(database)
    target = database / "input/kinetics/families/H_Abstraction/groups.py"
    if tamper == "edit":
        target.write_text(target.read_text() + "# x\n")
    elif tamper == "replace":
        replacement = target.with_name("replacement.tmp")
        replacement.write_text("name = 'evil'\n")
        os.replace(replacement, target)
    elif tamper == "replace_with_symlink":
        evil = tmp_path / "evil.py"
        evil.write_text("name = 'h' \n")
        target.unlink()
        target.symlink_to(evil)
    elif tamper == "add_kinetics":
        (database / "input/kinetics/families/H_Abstraction/rules.py").write_text("x = 1\n")
    elif tamper == "add_thermo":
        (database / "input/thermo/extra.py").write_text("x = 1\n")
    elif tamper == "remove":
        target.unlink()
    else:
        moved = database / "input/kinetics/families/H_Abstraction_old"
        (database / "input/kinetics/families/H_Abstraction").rename(moved)
        make_database_again = database / "input/kinetics/families/H_Abstraction"
        make_database_again.mkdir()
        (make_database_again / "groups.py").write_text("name = 'h' \nplanted = 1\n")
    with pytest.raises(gi.SnapshotMismatchError):
        gi.verify_database_snapshot(manifest, database)


def test_snapshot_manifest_hash_is_tamper_evident(tmp_path):
    database = make_database(tmp_path / "db")
    manifest = gi.database_snapshot_manifest(database)
    forged = copy.deepcopy(manifest)
    forged["files"]["input/thermo/library.py"] = "0" * 64
    with pytest.raises(gi.SnapshotMismatchError):
        gi.verify_database_snapshot(forged, database)
    forged = copy.deepcopy(manifest)
    forged["manifest_sha256"] = "0" * 64
    with pytest.raises(gi.SnapshotMismatchError):
        gi.verify_database_snapshot(forged, database)


def test_snapshot_manifest_refuses_missing_and_empty_databases(tmp_path):
    with pytest.raises(gi.IdentityError):
        gi.database_snapshot_manifest(tmp_path / "missing")
    (tmp_path / "empty").mkdir()
    with pytest.raises(gi.IdentityError):
        gi.database_snapshot_manifest(tmp_path / "empty")


# ---------------------------------------------------------------------------
# Request-binding envelope
# ---------------------------------------------------------------------------

REQUEST_A = "a" * 64
REQUEST_B = "b" * 64
PRODUCER = "c" * 64


def test_envelope_round_trips_its_payload():
    envelope = gi.build_envelope(REQUEST_A, PRODUCER, b'{"x":1}')
    assert set(envelope) == {"schema_version", "request_sha", "producer_identity", "payload_sha256", "payload"}
    assert gi.verify_envelope(envelope, request_sha=REQUEST_A, producer_identity=PRODUCER) == b'{"x":1}'
    again = gi.envelope_from_bytes(gi.envelope_to_bytes(envelope))
    assert gi.verify_envelope(again, request_sha=REQUEST_A) == b'{"x":1}'


def test_payload_bound_to_request_a_fails_against_request_b():
    envelope = gi.build_envelope(REQUEST_A, PRODUCER, b"payload")
    with pytest.raises(gi.EnvelopeError):
        gi.verify_envelope(envelope, request_sha=REQUEST_B)


def test_envelope_rejects_tampering_and_malformed_shapes():
    envelope = gi.build_envelope(REQUEST_A, PRODUCER, b"payload")
    cases = {
        "payload": dict(envelope, payload="tampered"),
        "payload digest": dict(envelope, payload_sha256="0" * 64),
        "schema": dict(envelope, schema_version=99),
        "producer": dict(envelope, producer_identity="d" * 64),
        "extra key": dict(envelope, extra=1),
        "missing key": {k: v for k, v in envelope.items() if k != "payload"},
        "bad request type": dict(envelope, request_sha=5),
        "short request sha": dict(envelope, request_sha="abc"),
        "non-dict": [envelope],
    }
    for name, candidate in cases.items():
        with pytest.raises(gi.EnvelopeError):
            gi.verify_envelope(candidate, request_sha=REQUEST_A, producer_identity=PRODUCER)
    assert name


def test_envelope_bytes_reject_duplicate_keys_and_garbage():
    good = gi.envelope_to_bytes(gi.build_envelope(REQUEST_A, PRODUCER, b"p"))
    duplicated = good.replace(b'"schema_version":1', b'"schema_version":1,"schema_version":1')
    assert duplicated != good
    for bad in (duplicated, b"\xff", b"[]", b"", b"{"):
        with pytest.raises(gi.EnvelopeError):
            gi.envelope_from_bytes(bad)


def test_envelope_builder_requires_full_hex_shas_and_bytes():
    with pytest.raises(gi.EnvelopeError):
        gi.build_envelope("short", PRODUCER, b"x")
    with pytest.raises(gi.EnvelopeError):
        gi.build_envelope(REQUEST_A, PRODUCER, "not bytes")  # type: ignore[arg-type]


def test_producer_identities_come_from_the_manifests():
    assert gi.generation_producer_identity() == gi.manifest_sha(gi.generation_dependency_manifest())
    assert gi.compilation_producer_identity() == gi.manifest_sha(gi.compilation_dependency_manifest())


def test_rules_holding_reactions_after_rate_rule_preparation_are_fingerprintable(kinetics):
    from rmgpy.reaction import Reaction

    family = copy.deepcopy(kinetics.families["R_Recombination"])
    entry = next(e for rules in family.rules.entries.values() for e in rules)
    entry.item = Reaction(reactants=[species("C[CH2]")], products=[species("CC")], reversible=True)
    first = gi.family_content_fingerprint(family)
    assert first == gi.family_content_fingerprint(family)
    entry.item.products = [species("CCC")]
    assert gi.family_content_fingerprint(family) != first
