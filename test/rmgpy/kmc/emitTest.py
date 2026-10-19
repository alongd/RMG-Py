"""Tests for rmgpy/kmc/emit.py (I-155): a compiled event-set artifact plus a
roster of declared reference systems becomes one gate-checked template catalog.

Database-free.  The reference systems are declared as a ``LinearHomopolymerSpec``
mapping plus chain variants; ``emit`` resolves the repeat unit to a registered
provider by graph isomorphism.  Assembly runs on stored row fixtures
(``fixtures/i155_rows_dp8.jsonl``: rows of the I-138 acceptance run on the closed
DP8 reference, template maps dropped) with stub services; the database-bound
balance functions (``derive_weights``, ``check_cycles``) are replaced where a
test needs a chosen verdict.
"""

import ast
import copy
import inspect
import json
import re
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from rmgpy.kmc import boundaries as b
from rmgpy.kmc import declared_chain as dc
from rmgpy.kmc import emit
from rmgpy.kmc import template_emission as te
from rmgpy.kmc import windows as win
from rmgpy.polymer_spec import LinearHomopolymerSpec

HERE = Path(__file__).resolve().parent
ROWS = HERE / "fixtures" / "i155_rows_dp8.jsonl"
S = b.PsChainSpec


# --------------------------------------------------------------------------
# declared polymers
# --------------------------------------------------------------------------
def _atom(key, element="C", hydrogens=0, radicals=0):
    return {"key": key, "element": element, "hydrogens": hydrogens,
            "electronic_state": {"radical_electrons": radicals}}


def _cap():
    return {"atoms": [_atom("h", "H", 0)], "bonds": [],
            "ports": {"cap": {"atom": "h", "bond_order": 1}},
            "electronic_state": {"multiplicity": 1}}


def polymer_mapping(kind):
    """A ``LinearHomopolymerSpec`` mapping: ``ps`` (styrene repeat), ``pe``
    (ethylene), ``pp`` (propylene: no provider declares it)."""
    ring = [_atom(f"p{i}", hydrogens=0 if i == 1 else 1) for i in range(1, 7)]
    ring_bonds = [{"atoms": [f"p{i}", f"p{i % 6 + 1}"], "order": 1.5} for i in range(1, 7)]
    if kind == "ps":
        atoms = [_atom("c1", hydrogens=2), _atom("c2", hydrogens=1)] + ring
        bonds = [{"atoms": ["c1", "c2"], "order": 1}, {"atoms": ["c2", "p1"], "order": 1}] + ring_bonds
    elif kind == "pe":
        atoms = [_atom("c1", hydrogens=2), _atom("c2", hydrogens=2)]
        bonds = [{"atoms": ["c1", "c2"], "order": 1}]
    elif kind == "pp":
        atoms = [_atom("c1", hydrogens=2), _atom("c2", hydrogens=1), _atom("c3", hydrogens=3)]
        bonds = [{"atoms": ["c1", "c2"], "order": 1}, {"atoms": ["c2", "c3"], "order": 1}]
    else:
        raise ValueError(kind)
    repeat = {"atoms": atoms, "bonds": bonds,
              "ports": {"head": {"atom": "c1", "bond_order": 1},
                        "tail": {"atom": "c2", "bond_order": 1}},
              "electronic_state": {"multiplicity": 1}}
    return {
        "schema_version": 1, "label": kind, "repeat_unit": repeat,
        "end_groups": {"head": _cap(), "tail": _cap()},
        "small_species": {},
        "initial_ensemble": {"dp_counts": {"8": 1}},
        "initial_mobile_counts": {}, "density_basis": "whole_initial_mixture",
        "initial_mass_kg": 1.25, "density_kg_m3": 1050.0,
        "locality_policy": {"radius": 2, "bounds": {
            "elements": ["C", "H"], "isotopes": [-1], "formal_charges": [0],
            "radical_electrons": [0, 1], "lone_pairs": [0],
            "bond_orders": [1, 1.5, 2, 3], "spin_multiplicities": [1, 2, 3],
            "maximum_heavy_atom_degree": 4,
            "maximum_explicit_hydrogens_per_heavy_atom": 4,
            "maximum_total_hydrogens_on_free_small_species": 12,
            "maximum_family_arity": 2, "maximum_root_tuple_size": 4,
            "maximum_free_small_species_heavy_atoms": 8,
            "attachment_classes": [], "topology_metadata_schema": {}}},
        "family_manifest": {}, "library_manifest": {}, "transport_model": {}, "diagnostics": {},
    }


def chain(dp, head=None, tail=None, radical_unit=None, owner=None):
    value = {"dp": dp, "head": head, "tail": tail, "radical_unit": radical_unit}
    if owner is not None:
        value["owner"] = owner
    return {"chain": value}


def write_roster(tmp_path, kind, systems, name="roster.json"):
    path = tmp_path / name
    default_owner = {"ps": "polystyrene", "pe": "linear_alkane",
                     "pp": "missing-provider"}[kind]
    systems = copy.deepcopy(systems)
    for entries in systems.values():
        for entry in entries:
            if "chain" in entry:
                entry["chain"].setdefault("owner", default_owner)
    path.write_text(json.dumps({"schema": emit.ROSTER_SCHEMA, "polymer": polymer_mapping(kind),
                                "systems": systems}))
    return path


def systems_of(tmp_path, kind, systems):
    roster = emit.load_roster(write_roster(tmp_path, kind, systems))
    return emit.build_reference_systems(roster)


def generated_of(tmp_path, kind, **kwargs):
    spec = LinearHomopolymerSpec.from_mapping(polymer_mapping(kind))
    provider = emit.match_provider(spec.repeat_unit)
    generated = emit.reference_systems(spec, provider, **kwargs)
    path = tmp_path / f"{kind}.generated.roster.json"
    written = emit.write_roster(generated, path)
    systems, refusals = emit.build_reference_systems(written)
    assert refusals == []
    return generated, written, systems


def generated_participants(systems):
    return [participants[0] for name, participants in systems.items()
            if name.startswith("single:")]


def contains_species(participants, expected):
    built = win.build_spec(expected, "legacy")
    return any(isomorphic(participant.build("generated"), built)
               for participant in participants)


def isomorphic(build_a, build_b):
    return build_a.molecule.is_isomorphic(build_b.molecule)


# --------------------------------------------------------------------------
# item 1: reference systems from a declaration
# --------------------------------------------------------------------------
def test_ps_roster_chains_use_one_declared_provider(tmp_path):
    systems, refusals = systems_of(tmp_path, "ps", {
        "dp8_closed": [chain(8)],
        "dp4+pe": [chain(4), chain(1, tail="benzylic_radical")],
        "pe2_identical": [chain(2, tail="benzylic_radical"), chain(2, tail="benzylic_radical")]})
    assert refusals == []
    specs = [participant.spec for participants in systems.values() for participant in participants]
    assert all(isinstance(spec, win.DeclaredSpec) for spec in specs)
    assert len({spec.provider for spec in specs}) == 1
    assert specs[0].provider.startswith("declared:")


def test_ps_declaration_molecules_are_isomorphic_to_the_harness_ones(tmp_path):
    systems, _ = systems_of(tmp_path, "ps", {
        "a": [chain(8)], "b": [chain(4), chain(1, tail="benzylic_radical")],
        "c": [chain(2, tail="benzylic_radical")]})
    pairs = [(systems["a"][0], S(8)), (systems["b"][0], S(4)),
             (systems["b"][1], S(1, tail_end="benzylic_radical")),
             (systems["c"][0], S(2, tail_end="benzylic_radical"))]
    for participant, spec in pairs:
        assert isomorphic(participant.build("p0"), win.build_spec(spec, "p0"))


def test_pe_declaration_is_the_n_alkane(tmp_path):
    systems, refusals = systems_of(tmp_path, "pe", {"c16": [chain(8)]})
    assert refusals == []
    (participant,) = systems["c16"]
    assert isinstance(participant.spec, win.DeclaredSpec)
    assert participant.spec.provider.startswith("declared:")
    build = participant.build("p0")
    assert isomorphic(build, win.build_spec(("alkane", 16), "p0"))
    assert build.molecule.get_formula() == "C16H34"


def test_pe_declaration_with_ends_is_the_alkane_end_spec(tmp_path):
    systems, _ = systems_of(tmp_path, "pe", {"r": [chain(3, head="primary_radical", tail="methyl")]})
    participant = systems["r"][0]
    assert isinstance(participant.spec, win.DeclaredSpec)
    assert isomorphic(participant.build("declared"),
                      win.build_spec(("alkane", 6, "primary_radical", "methyl", None)))


def test_unregistered_vocabulary_owner_is_refused_with_a_row(tmp_path):
    systems, refusals = systems_of(tmp_path, "pp", {"pp8": [chain(8)]})
    assert systems == {} or systems.get("pp8") in (None, [])
    assert [r["status"] for r in refusals] == ["refused:owner"]
    assert refusals[0]["system"] == "pp8" and refusals[0]["reason"]


def test_unknown_end_variant_is_refused_never_dropped(tmp_path):
    systems, refusals = systems_of(tmp_path, "ps", {"x": [chain(4, tail="bogus")]})
    assert [r["status"] for r in refusals] == ["refused:variant"]
    assert "x" not in systems


def test_interior_radical_on_the_alkane_is_built_open_shell(tmp_path):
    systems, refusals = systems_of(tmp_path, "pe", {"r": [chain(5, radical_unit=2)]})
    assert refusals == []
    participant = systems["r"][0]
    assert isinstance(participant.spec, win.DeclaredSpec)
    molecule = participant.build("r").molecule
    assert molecule.get_formula() == "C10H21"
    assert molecule.get_radical_count() == 1


def test_ps_interior_radical_is_the_radical_unit_spec(tmp_path):
    systems, _ = systems_of(tmp_path, "ps", {"r": [chain(5, radical_unit=2)]})
    participant = systems["r"][0]
    assert isinstance(participant.spec, win.DeclaredSpec)
    assert isomorphic(participant.build("declared"),
                      win.build_spec(S(5, radical_units=(2,))))


def test_a_refused_participant_refuses_its_whole_system(tmp_path):
    systems, refusals = systems_of(tmp_path, "ps", {
        "mixed": [chain(4), chain(2, tail="bogus")], "ok": [chain(2)]})
    assert "mixed" not in systems and "ok" in systems
    assert len(refusals) == 1


def test_explicit_small_molecules(tmp_path):
    systems, refusals = systems_of(tmp_path, "ps", {"s": [
        {"smiles": "C"},
        {"adjacency_list": "1 H u1 p0 c0"}]})
    assert refusals == []
    methane, h_atom = systems["s"]
    assert methane.build("p0").molecule.get_formula() == "CH4"
    assert h_atom.build("p1").molecule.get_formula() == "H"
    assert methane.build("p0").metadata.provider == "small"


def test_roster_sha_is_the_file_hash(tmp_path):
    import hashlib
    path = write_roster(tmp_path, "ps", {"a": [chain(2)]})
    assert emit.load_roster(path).sha256 == hashlib.sha256(path.read_bytes()).hexdigest()


def test_generated_roster_counts_are_pinned_per_declaration(tmp_path):
    ps, ps_written, ps_systems = generated_of(tmp_path, "ps", dps=(1, 2))
    pe, pe_written, pe_systems = generated_of(tmp_path, "pe", dps=(1, 2))
    assert (len(ps.systems), len(ps_systems)) == (440, 440)
    assert (len(pe.systems), len(pe_systems)) == (173, 173)
    assert ps.sha256 == ps_written.sha256
    assert pe.sha256 == pe_written.sha256


def test_generated_pp_roster_loads_and_builds_in_a_fresh_process(tmp_path):
    spec = LinearHomopolymerSpec.from_mapping(polymer_mapping("pp"))
    roster = emit.reference_systems(spec, emit.match_provider(spec.repeat_unit), dps=(1, 2),
                                    small=False, monomer=False, molecularity={1})
    path = tmp_path / "pp.generated.roster.json"
    emit.write_roster(roster, path)
    script = """
import json
import sys
from pathlib import Path
from rmgpy.kmc import boundaries, emit

path = Path(sys.argv[1])
data = json.loads(path.read_text())
owners = {entry["chain"]["owner"] for entries in data["systems"].values()
          for entry in entries if "chain" in entry}
assert owners and not (owners & set(boundaries.PROVIDERS))
roster = emit.load_roster(path)
assert owners <= set(boundaries.PROVIDERS)
systems, refusals = emit.build_reference_systems(roster)
assert not refusals
assert len(systems) == len(roster.systems)
for name, participants in systems.items():
    for index, participant in enumerate(participants):
        participant.build(f"{name}:{index}")
"""
    result = subprocess.run([sys.executable, "-c", script, str(path)],
                            cwd=Path(__file__).parents[3], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("kind", ("ps", "pe"))
def test_generated_participants_round_trip_the_exact_vocabulary_owner(tmp_path, kind):
    spec = LinearHomopolymerSpec.from_mapping(polymer_mapping(kind))
    declared = emit.declared_provider(spec)
    roster = emit.reference_systems(spec, declared, dps=(1, 2, 3, 5, 8), molecularity={1})
    expected = {}
    for name, (entry,) in roster.systems.items():
        if "chain" in entry:
            value = entry["chain"]
            declared_spec = declared.spec_from_declaration(
                value["dp"], value.get("head"), value.get("tail"),
                value.get("radical_unit"))
            build = win.build_spec(declared_spec, f"generated:{name}")
        else:
            (kind_, value), = entry.items()
            build = emit.explicit_build(f"generated:{name}", **{kind_: value})
        expected[name] = build.molecule.to_adjacency_list()

    written = emit.write_roster(roster, tmp_path / f"{kind}.owner-round-trip.json")
    systems, refusals = emit.build_reference_systems(written)
    assert refusals == []
    assert set(systems) == set(expected)
    for name, participants in systems.items():
        assert len(participants) == 1
        assert participants[0].build(name).molecule.to_adjacency_list() == expected[name]


def test_roster_refuses_a_vocabulary_owner_for_another_repeat_unit(tmp_path):
    spec = LinearHomopolymerSpec.from_mapping(polymer_mapping("ps"))
    roster = emit.reference_systems(spec, emit.declared_provider(spec), dps=(2,),
                                    small=False, monomer=False, molecularity={1})
    data = json.loads(emit.roster_bytes(roster))
    name = "single:chain:dp2:interior-radical"
    data["systems"] = {name: data["systems"][name]}
    data["systems"][name][0]["chain"]["owner"] = "linear_alkane"
    path = tmp_path / "wrong-owner.json"
    path.write_text(json.dumps(data))

    systems, refusals = emit.build_reference_systems(emit.load_roster(path))
    assert systems == {}
    assert [row["status"] for row in refusals] == ["refused:owner"]
    assert "linear_alkane" in refusals[0]["reason"]


@pytest.mark.parametrize("kind", ("ps", "pe"))
def test_every_generated_system_uses_only_its_declaration_chain_provider(tmp_path, kind):
    _roster, _written, systems = generated_of(tmp_path, kind, dps=(1, 2, 3, 5, 8))
    chain_providers = set()
    for name, participants in systems.items():
        builds = [participant.build(f"{name}:{index}")
                  for index, participant in enumerate(participants)]
        win._builds_policy(builds)
        chain_providers.update(build.metadata.provider for build in builds
                               if build.metadata.provider != "small")
    assert len(chain_providers) == 1
    assert next(iter(chain_providers)).startswith("declared:")


def test_ps_generated_roster_contains_the_recombination_legacy_oracle(tmp_path):
    _roster, _written, systems = generated_of(tmp_path, "ps", dps=(1, 2))
    participants = generated_participants(systems)
    expected = []
    for kind in win.PARTICIPANT_KINDS:
        for dp in (1, 2):
            expected.extend(win.kind_specs(kind, dp))
    expected.extend(make() for _name, make in win.PROBE_CANDIDATES)
    missing = [repr(spec) for spec in expected if not contains_species(participants, spec)]
    assert missing == []


def test_generated_radical_pairs_have_both_orders_and_identical_pairs(tmp_path):
    _roster, _written, systems = generated_of(tmp_path, "ps", dps=(1,))
    pairs = {name: participants for name, participants in systems.items()
             if name.startswith("pair:radical:")}
    unequal = next(participants for participants in pairs.values()
                   if not isomorphic(participants[0].build("a"), participants[1].build("b")))
    assert any(isomorphic(participants[0].build("a"), unequal[1].build("b")) and
               isomorphic(participants[1].build("b"), unequal[0].build("a"))
               for participants in pairs.values())
    assert any(isomorphic(participants[0].build("a"), participants[1].build("b"))
               for participants in pairs.values())


def test_generated_radicals_pair_with_every_chain_role_in_both_orders():
    spec = LinearHomopolymerSpec.from_mapping(polymer_mapping("ps"))
    roster = emit.reference_systems(spec, emit.match_provider(spec.repeat_unit), dps=(1,))
    singles = {name.removeprefix("single:"): entries[0]
               for name, entries in roster.systems.items()
               if name.startswith("single:chain:")}
    radicals = {label: entry for label, entry in singles.items()
                if "radical" in label}
    assert radicals
    for radical_label, radical in radicals.items():
        for chain_label, chain_entry in singles.items():
            assert any(entries == [radical, chain_entry]
                       for entries in roster.systems.values()), (radical_label, chain_label, "left")
            assert any(entries == [chain_entry, radical]
                       for entries in roster.systems.values()), (radical_label, chain_label, "right")


def test_declared_provider_builds_the_declarations_caps_through_emit(tmp_path):
    mapping = polymer_mapping("pe")
    mapping["end_groups"]["head"] = {
        "atoms": [_atom("methyl", hydrogens=3)], "bonds": [],
        "ports": {"cap": {"atom": "methyl", "bond_order": 1}},
        "electronic_state": {"multiplicity": 1}}
    spec = LinearHomopolymerSpec.from_mapping(mapping)
    roster = emit.reference_systems(spec, emit.match_provider(spec.repeat_unit), dps=(3,),
                                    small=False, monomer=False, molecularity={1})
    systems, refusals = emit.build_reference_systems(
        emit.write_roster(roster, tmp_path / "declared-pe.json"))
    assert refusals == []
    molecule = systems["single:chain:dp3:declared"][0].build("declared").molecule
    assert molecule.get_formula() == "C7H16"


def test_generated_roster_skips_an_unbuildable_interior_radical(monkeypatch):
    spec = LinearHomopolymerSpec.from_mapping(polymer_mapping("pe"))
    provider = emit.match_provider(spec.repeat_unit)
    original = provider.spec_for_declared_provider

    def without_interior(declared, dp, head=None, tail=None, interior_radical=None):
        if interior_radical is not None:
            raise b.ExtensionRefused("no hydrogen-bearing interior backbone atom")
        return original(declared, dp, head, tail, interior_radical)

    monkeypatch.setattr(provider, "spec_for_declared_provider", without_interior)
    roster = emit.reference_systems(spec, provider, dps=(2,), small=False, monomer=False,
                                    molecularity={1})
    assert not any("interior-radical" in name for name in roster.systems)


def test_generated_explicit_species_are_valid_window_specs(tmp_path):
    mapping = polymer_mapping("pe")
    mapping["small_species"] = {
        "methane": {"graph": {"atoms": [_atom("c", hydrogens=4)], "bonds": [],
                                "ports": {}, "electronic_state": {"multiplicity": 1}},
                    "disposition": "mobile"}}
    spec = LinearHomopolymerSpec.from_mapping(mapping)
    roster = emit.reference_systems(spec, emit.match_provider(spec.repeat_unit), dps=(1,))
    systems, refusals = emit.build_reference_systems(
        emit.write_roster(roster, tmp_path / "small-window.json"))
    assert refusals == []
    for name, formula in (("single:small:H", "H"), ("single:small:methane", "CH4")):
        participant = systems[name][0]
        (rebuilt,) = win.builds_of([participant.spec])
        assert rebuilt.molecule.get_formula() == formula


def test_declared_small_species_multiplicity_survives_roster_round_trip(tmp_path):
    mapping = polymer_mapping("pe")
    atoms = [_atom("a", hydrogens=2, radicals=1),
             _atom("b", hydrogens=2, radicals=1)]
    mapping["small_species"] = {
        name: {"graph": {"atoms": atoms, "bonds": [{"atoms": ["a", "b"], "order": 1}],
                         "ports": {}, "electronic_state": {"multiplicity": multiplicity}},
               "disposition": "mobile"}
        for name, multiplicity in (("singlet", 1), ("triplet", 3))}
    spec = LinearHomopolymerSpec.from_mapping(mapping)
    roster = emit.reference_systems(spec, emit.match_provider(spec.repeat_unit), dps=(1,))
    systems, refusals = emit.build_reference_systems(
        emit.write_roster(roster, tmp_path / "multiplicity.json"))
    assert refusals == []
    assert systems["single:small:singlet"][0].build("s").molecule.multiplicity == 1
    assert systems["single:small:triplet"][0].build("t").molecule.multiplicity == 3


def test_generated_interior_radical_is_not_a_terminal_backbone_atom(tmp_path):
    _roster, _written, systems = generated_of(tmp_path, "pe", dps=(1, 2))
    assert "single:chain:dp1:interior-radical" not in systems
    molecule = systems["single:chain:dp2:interior-radical"][0].build("r").molecule
    (radical,) = [atom for atom in molecule.atoms if atom.radical_electrons]
    assert sum(not atom.is_hydrogen() for atom in radical.bonds) == 2


def test_mismatched_explicit_monomer_is_a_named_refusal():
    mapping = polymer_mapping("pe")
    mapping["monomer_product"] = {
        "atoms": [_atom("a", hydrogens=3), _atom("b", hydrogens=3)],
        "bonds": [{"atoms": ["a", "b"], "order": 1}], "ports": {},
        "electronic_state": {"multiplicity": 1}}
    spec = LinearHomopolymerSpec.from_mapping(mapping)
    with pytest.raises(emit.RosterRefusal, match="refused:monomer"):
        emit.reference_systems(spec, emit.match_provider(spec.repeat_unit), dps=(1,))


def test_declared_small_species_cannot_shadow_builtin_h():
    mapping = polymer_mapping("pe")
    mapping["small_species"] = {
        "H": {"graph": {"atoms": [_atom("h", "H", radicals=1)], "bonds": [],
                          "ports": {}, "electronic_state": {"multiplicity": 2}},
              "disposition": "mobile"}}
    spec = LinearHomopolymerSpec.from_mapping(mapping)
    with pytest.raises(emit.RosterRefusal, match="refused:small_name_collision"):
        emit.reference_systems(spec, emit.match_provider(spec.repeat_unit), dps=(1,))


def test_pe_generated_roster_has_the_same_core_roles_without_ps_provider(tmp_path, monkeypatch):
    def ps_reached(*_args, **_kwargs):
        raise AssertionError("the PE path reached the PS provider")
    monkeypatch.setattr(b.PROVIDERS["polystyrene"], "spec_from_declaration", ps_reached)
    _roster, _written, systems = generated_of(tmp_path, "pe", dps=(1, 2))
    names = set(systems)
    for role in ("capped", "declared", "end-radical", "interior-radical",
                 "monomer", "small:H"):
        assert any(role in name for name in names), role


def test_declared_small_species_are_singletons_and_pair_with_each_chain(tmp_path):
    mapping = polymer_mapping("pe")
    mapping["small_species"] = {
        "methane": {"graph": {"atoms": [_atom("c", hydrogens=4)], "bonds": [],
                                "ports": {}, "electronic_state": {"multiplicity": 1}},
                    "disposition": "mobile"}}
    spec = LinearHomopolymerSpec.from_mapping(mapping)
    roster = emit.reference_systems(spec, emit.match_provider(spec.repeat_unit), dps=(1,))
    path = tmp_path / "pe-small.roster.json"
    systems, refusals = emit.build_reference_systems(emit.write_roster(roster, path))
    assert refusals == []
    assert "single:small:methane" in systems
    chain_names = [name for name in systems if name.startswith("single:chain:")]
    pairs = [name for name in systems if name.startswith("pair:small:methane:")]
    assert len(pairs) == 2 * len(chain_names)
    methane = systems["single:small:methane"][0].build("methane").molecule
    assert methane.get_formula() == "CH4"


def test_no_chemistry_specific_code_in_emit():
    text = Path(emit.__file__).read_text()
    found = re.findall(r"polystyrene|PsChainSpec|alkane|polyethylene|styrene", text, flags=re.IGNORECASE)
    found += re.findall(r"\bPS\b|\bPE\b", text)
    assert not found, found


# --------------------------------------------------------------------------
# harness equivalence of the moved pure functions
# --------------------------------------------------------------------------
def _harness():
    sys.path.insert(0, str(HERE))
    import template_emission_acceptance as acc
    return acc


@pytest.mark.parametrize("name", ["record_key", "bind", "molecularity"])
def test_moved_harness_functions_are_verbatim(name):
    def body(function):
        node = ast.parse(textwrap.dedent(inspect.getsource(function))).body[0]
        statements = node.body
        if isinstance(statements[0], ast.Expr) and isinstance(statements[0].value, ast.Constant):
            statements = statements[1:]       # the docstring
        return ast.dump(node.args), [ast.dump(x) for x in statements]
    assert body(getattr(_harness(), name)) == body(getattr(emit, name))


# --------------------------------------------------------------------------
# record selection
# --------------------------------------------------------------------------
def test_chunks_partition_the_sorted_records():
    records = [{"event_id": f"evt_{k:03d}", "family": "F"} for k in (5, 1, 9, 3, 7, 2, 8)]
    parts = [emit.select_records(records, f"{k}/3") for k in range(3)]
    ids = sorted(r["event_id"] for part in parts for r in part)
    assert ids == sorted(r["event_id"] for r in records)
    assert all(not ({r["event_id"] for r in a} & {r["event_id"] for r in c})
               for i, a in enumerate(parts) for c in parts[i + 1:])
    assert [r["event_id"] for r in parts[0]] == ["evt_001", "evt_005", "evt_009"]


# --------------------------------------------------------------------------
# assembly
# --------------------------------------------------------------------------
class Services:
    counter = None
    kc = None
    database_content_sha256 = "h" * 64


def ok_midpoint(row, worst=0.0):
    """The row's per-template midpoint entries, as the emit step stores them."""
    return {t["key"]: {"status": "ok", "worst": worst, "worst_T": 650.0, "n": 28}
            for t in (row.get("payload") or {}).get("templates", [])
            if t["status"] == "emitted" and t.get("weight")}


def fixture_rows():
    rows = [json.loads(line) for line in ROWS.read_text().splitlines()]
    rows.sort(key=lambda r: r["event_id"])
    for row in rows:
        if row.get("payload"):
            row["midpoint"] = ok_midpoint(row)
    return rows


def record_id(row):
    return (row.get("payload") or {}).get("record_event_id") or row["event_id"] + "0" * 44


def header(chunk, rows, **over):
    out = {"kind": "header", "schema": emit.ROWS_SCHEMA, "artifact_id": "art",
           "artifact_sha256": "a" * 64, "artifact_records": 6, "families": ["R_Recombination"],
           "database_path": "/db", "database_head": "d" * 40, "code_sha": "c" * 40,
           "database_content_sha256": "h" * 64, "code_sources_sha256": "s" * 64,
           "code_extensions_sha256": "x" * 64,
           "roster_sha256": "r" * 64, "temperatures": [600.0, 700.0, 800.0], "chunk": chunk,
           "selection": "all", "assigned": [record_id(r) for r in rows], "roster_refusals": [],
           "systems": {"dp8_closed": ["chain dp=8"]}}
    out.update(over)
    return out


def write_rows(path, head, rows):
    with open(path, "w") as handle:
        handle.write(json.dumps(head) + "\n")
        for row in rows:
            handle.write(json.dumps(dict(row, kind="row")) + "\n")
    return path


def shards(tmp_path, n=2, rows=None, **over):
    rows = rows or fixture_rows()
    paths = []
    for k in range(n):
        mine = rows[k::n]
        paths.append(write_rows(tmp_path / f"rows.{k}of{n}.jsonl", header(f"{k}/{n}", mine, **over), mine))
    return paths


ASSEMBLER = {"database_content_sha256": "h" * 64, "code_sources_sha256": "s" * 64,
             "code_extensions_sha256": "x" * 64}


@pytest.fixture(autouse=True)
def assembler_identity(monkeypatch):
    """The assembling process's own database and code digests: those of the
    fixture shards unless a test sets others."""
    current = dict(ASSEMBLER)
    monkeypatch.setattr(emit, "assembler_identity", lambda database: dict(current), raising=False)
    return current


@pytest.fixture(autouse=True)
def balance_stubs(monkeypatch):
    """The database-bound balance functions: identity unless a test sets a verdict."""
    calls = {"derive": 0, "cycles": 0}

    def derive(payloads, counter):
        calls["derive"] += 1
        return payloads

    def cycles(payloads, counter, kc, complete=True, **kw):
        calls["cycles"] += 1
        return payloads, []
    monkeypatch.setattr(te, "derive_weights", derive)
    monkeypatch.setattr(te, "check_cycles", cycles)
    return calls


def run(paths, **kw):
    return emit.assemble([str(p) for p in paths], services=Services(), **kw)


def test_catalog_schema_and_provenance(tmp_path):
    catalog = run(shards(tmp_path))
    assert catalog["schema"] == "kmc-template-catalog/1"
    assert set(catalog) == {"schema", "provenance", "payloads", "refusals", "gates",
                            "completeness", "complete"}
    prov = catalog["provenance"]
    assert prov["artifact_id"] == "art" and prov["artifact_sha256"] == "a" * 64
    assert prov["database_path"] == "/db" and prov["database_head"] == "d" * 40
    assert prov["code_sha"] == "c" * 40 and prov["roster_sha256"] == "r" * 64
    assert prov["families"] == ["R_Recombination"]
    assert prov["temperatures"] == [600.0, 700.0, 800.0]
    assert set(catalog["gates"]) >= {"collisions", "cycles", "balance", "failures"}
    assert set(catalog["gates"]["collisions"]) >= {"keys", "merged", "refused"}
    assert set(catalog["gates"]["cycles"]) >= {"pairs", "refused"}
    assert catalog["complete"] is True and catalog["gates"]["failures"] == []


def test_payloads_hold_exactly_the_fields_loaded_templates_reads(tmp_path):
    catalog = run(shards(tmp_path))
    assert catalog["payloads"]
    for payload in catalog["payloads"]:
        assert set(payload) == {"record_event_id", "family", "reverse_of", "templates"}
        for t in payload["templates"]:
            assert t["status"] == "emitted"
            assert set(t) <= {"key", "name", "family", "status", "direction", "executable",
                              "weight", "derived"}


def test_refusals_and_completeness(tmp_path):
    catalog = run(shards(tmp_path))
    refused = {r["record_event_id"][:24] for r in catalog["refusals"]}
    assert refused == {r["event_id"] for r in fixture_rows() if r["status"] == "refused"}
    assert all(r["status"] and r["reason"] is not None for r in catalog["refusals"])
    totals = {s: sum(c[s] for c in catalog["completeness"].values())
              for s in ("emitted", "partial", "refused", "error")}
    assert totals == {"emitted": 4, "partial": 0, "refused": 2, "error": 0}
    assert all(k.startswith("R_Recombination|") for k in catalog["completeness"])


def test_two_shards_equal_one_unsharded_row_set(tmp_path):
    rows = fixture_rows()
    one = write_rows(tmp_path / "all.jsonl", header(None, rows), rows)
    first = emit.dumps(run([one]))
    two = emit.dumps(run(shards(tmp_path, 2)))
    three = emit.dumps(run(list(reversed(shards(tmp_path, 3)))))
    assert first == two == three


def test_dropped_shard_is_incomplete_and_listed(tmp_path):
    paths = shards(tmp_path, 3)
    catalog = run(paths[:2])
    assert catalog["complete"] is False
    assert any("missing shard" in f and "2/3" in f for f in catalog["gates"]["failures"])


def test_duplicate_shard_is_listed(tmp_path):
    paths = shards(tmp_path, 2)
    catalog = run(paths + [paths[0]])
    assert catalog["complete"] is False
    assert any("duplicate" in f for f in catalog["gates"]["failures"])


def test_mixed_artifacts_are_listed(tmp_path):
    rows = fixture_rows()
    a = write_rows(tmp_path / "a.jsonl", header("0/2", rows[0::2]), rows[0::2])
    c = write_rows(tmp_path / "c.jsonl", header("1/2", rows[1::2], artifact_sha256="b" * 64), rows[1::2])
    catalog = run([a, c])
    assert catalog["complete"] is False
    assert any("artifact_sha256" in f for f in catalog["gates"]["failures"])


def test_error_row_is_incomplete_with_the_failure_listed(tmp_path):
    rows = fixture_rows()
    bad = {"event_id": "evt_ffffffffffffffffffff", "family": "R_Recombination",
           "orientation": "as_generated", "status": "error", "refusal": "error:exception",
           "reason": "Traceback: boom", "seconds": 0.1}
    rows = rows + [bad]
    bad_header = header("0/1", rows, artifact_records=7)
    path = write_rows(tmp_path / "e.jsonl", bad_header, rows)
    catalog = run([path])
    assert catalog["complete"] is False
    assert any("error" in f and "boom" in f for f in catalog["gates"]["failures"])
    assert any(r["status"] == "error" for r in catalog["refusals"])


def test_injected_collision_is_incomplete_with_the_failure_listed(tmp_path):
    rows = fixture_rows()
    emitted = [r for r in rows if r["status"] == "emitted"]
    donor, host = emitted[0], emitted[1]
    clone = copy.deepcopy(donor["payload"]["templates"][0])
    clone["weight"] = dict(clone["weight"], w=[w * 1.0e3 for w in clone["weight"]["w"]])
    host["payload"]["templates"] = host["payload"]["templates"] + [clone]
    host["payload"] = te.finish(host["payload"])
    path = write_rows(tmp_path / "c.jsonl", header(None, rows), rows)
    catalog = run([path])
    assert catalog["gates"]["collisions"]["refused"] >= 1
    assert catalog["complete"] is False
    assert any("collision" in f for f in catalog["gates"]["failures"])


def test_injected_cycle_failure_is_incomplete_with_the_failure_listed(tmp_path, monkeypatch):
    def cycles(payloads, counter, kc, complete=True, **kw):
        return payloads, [{"forward": "evt_" + "a" * 20, "reverse": "evt_" + "b" * 20,
                           "status": "refused:cycle", "reason": "ratio off by 1e3"}]
    monkeypatch.setattr(te, "check_cycles", cycles)
    catalog = run(shards(tmp_path))
    assert catalog["gates"]["cycles"]["pairs"] == 1 and catalog["gates"]["cycles"]["refused"] == 1
    assert catalog["complete"] is False
    assert any("cycle" in f and "ratio off" in f for f in catalog["gates"]["failures"])


def test_assemble_runs_derive_weights_then_check_cycles(tmp_path, balance_stubs):
    run(shards(tmp_path))
    assert balance_stubs == {"derive": 1, "cycles": 1}


def test_check_cycles_receives_completeness_of_coverage(tmp_path, monkeypatch):
    seen = []

    def cycles(payloads, counter, kc, complete=True, **kw):
        seen.append(complete)
        return payloads, []
    monkeypatch.setattr(te, "check_cycles", cycles)
    run(shards(tmp_path, 2))
    run(shards(tmp_path, 3)[:2])
    assert seen == [True, False]


def test_roster_refusal_makes_the_catalog_incomplete(tmp_path):
    refusal = {"kind": "roster_refusal", "system": "pp8", "participant": 0,
               "status": "refused:no_provider", "reason": "no provider declares this repeat unit"}
    catalog = run(shards(tmp_path, roster_refusals=[refusal]))
    assert catalog["complete"] is False
    assert any("refused:no_provider" in f and "pp8" in f for f in catalog["gates"]["failures"])


def test_rows_without_a_shard_header_are_incomplete(tmp_path):
    path = tmp_path / "bare.jsonl"
    path.write_text(ROWS.read_text())
    catalog = run([path])
    assert catalog["complete"] is False
    assert any("header" in f for f in catalog["gates"]["failures"])


def test_round_trip_loads_at_every_provenance_temperature(tmp_path):
    catalog = run(shards(tmp_path))
    expected = sum(1 for r in fixture_rows() if r["payload"]
                   for t in r["payload"]["templates"] if t["status"] == "emitted")
    # the fixture's two rows that emit one executable template are merged into one
    merged = catalog["gates"]["collisions"]["copies"] - catalog["gates"]["collisions"]["keys"]
    assert expected > 0 and merged >= 1
    expected -= merged
    for T in catalog["provenance"]["temperatures"]:
        assert len(te.loaded_templates(catalog["payloads"], T)) == expected
    assert catalog["gates"]["round_trip"]["templates"] == expected
    assert catalog["gates"]["round_trip"]["ok"] is True


def test_round_trip_failure_is_a_gate(tmp_path, monkeypatch):
    def broken(payloads, temperature, **kw):
        raise te.EmissionError("orphaned derived template")
    monkeypatch.setattr(te, "loaded_templates", broken)
    catalog = run(shards(tmp_path))
    assert catalog["complete"] is False
    assert any("round trip" in f for f in catalog["gates"]["failures"])


# --------------------------------------------------------------------------
# CLI end to end
# --------------------------------------------------------------------------
def test_cli_assemble_end_to_end(tmp_path, monkeypatch):
    monkeypatch.setattr(emit, "load_services", lambda database, families: Services())
    paths = [str(p) for p in shards(tmp_path, 2)]
    out_two, out_one = tmp_path / "two.json", tmp_path / "one.json"
    assert emit.main(["assemble", *paths, "--out", str(out_two)]) == 0
    rows = fixture_rows()
    one = write_rows(tmp_path / "all.jsonl", header(None, rows), rows)
    assert emit.main(["assemble", str(one), "--out", str(out_one)]) == 0
    assert out_two.read_bytes() == out_one.read_bytes()
    catalog = json.loads(out_two.read_text())
    assert catalog["schema"] == "kmc-template-catalog/1" and catalog["complete"] is True
    assert te.loaded_templates(catalog["payloads"], 700.0)


def test_cli_incomplete_catalog_is_written_and_exits_nonzero(tmp_path, monkeypatch):
    monkeypatch.setattr(emit, "load_services", lambda database, families: Services())
    paths = [str(p) for p in shards(tmp_path, 2)][:1]
    out = tmp_path / "partial.json"
    assert emit.main(["assemble", *paths, "--out", str(out)]) == 1
    catalog = json.loads(out.read_text())
    assert catalog["complete"] is False and catalog["gates"]["failures"]


def test_module_runs_as_a_script():
    import subprocess
    done = subprocess.run([sys.executable, "-m", "rmgpy.kmc.emit", "assemble", "--help"],
                          capture_output=True, text=True, timeout=120)
    assert done.returncode == 0 and "--out" in done.stdout


# --------------------------------------------------------------------------
# rework 1: the midpoint gate (computed per record at emission)
# --------------------------------------------------------------------------
def _emitted_row(rows, index=0):
    return [r for r in rows if r["status"] == "emitted"][index]


def test_midpoint_gate_aggregates_the_rows(tmp_path):
    catalog = run(shards(tmp_path))
    mid = catalog["gates"]["midpoint"]
    assert mid["evaluated"] is True and mid["threshold"] == 0.01
    assert mid["templates"] == catalog["gates"]["round_trip"]["templates"]
    assert mid["worst"] == 0.0 and mid["failing"] == [] and mid["unverified"] == []
    assert catalog["complete"] is True


def test_injected_midpoint_failure_is_incomplete_with_the_failure_listed(tmp_path):
    rows = fixture_rows()
    row = _emitted_row(rows)
    key = next(iter(row["midpoint"]))
    row["midpoint"][key] = {"status": "ok", "worst": 0.05, "worst_T": 650.0, "n": 28}
    path = write_rows(tmp_path / "m.jsonl", header(None, rows), rows)
    catalog = run([path])
    mid = catalog["gates"]["midpoint"]
    assert mid["worst"] == 0.05 and len(mid["failing"]) == 1
    assert catalog["complete"] is False
    assert any("midpoint" in f and "650" in f for f in catalog["gates"]["failures"])


def test_error_at_the_threshold_passes(tmp_path):
    rows = fixture_rows()
    row = _emitted_row(rows)
    key = next(iter(row["midpoint"]))
    row["midpoint"][key]["worst"] = 0.01
    assert run([write_rows(tmp_path / "m.jsonl", header(None, rows), rows)])["complete"] is True


def test_missing_midpoint_entry_is_unverified_and_incomplete(tmp_path):
    rows = fixture_rows()
    row = _emitted_row(rows)
    row["midpoint"].pop(next(iter(row["midpoint"])))
    catalog = run([write_rows(tmp_path / "m.jsonl", header(None, rows), rows)])
    mid = catalog["gates"]["midpoint"]
    assert mid["evaluated"] is False and len(mid["unverified"]) == 1
    assert catalog["complete"] is False
    assert any("midpoint check not evaluable" in f for f in catalog["gates"]["failures"])


def test_unavailable_midpoint_entry_is_unverified_and_incomplete(tmp_path):
    rows = fixture_rows()
    row = _emitted_row(rows)
    key = next(iter(row["midpoint"]))
    row["midpoint"][key] = {"status": "unavailable", "reason": "no emitted witness"}
    catalog = run([write_rows(tmp_path / "m.jsonl", header(None, rows), rows)])
    assert catalog["gates"]["midpoint"]["evaluated"] is False
    assert catalog["complete"] is False
    assert any("no emitted witness" in f for f in catalog["gates"]["failures"])


def test_a_derived_template_is_not_midpoint_gated(tmp_path, monkeypatch):
    def derive(payloads, counter):
        out = copy.deepcopy(payloads)
        out[0]["templates"][0]["derived"] = {"primary": {"record": "x", "key": "y"}}
        return out
    monkeypatch.setattr(te, "derive_weights", derive)
    rows = fixture_rows()
    for row in rows:
        row["midpoint"] = {}          # no entry at all for any template
    catalog = run([write_rows(tmp_path / "m.jsonl", header(None, rows), rows)])
    mid = catalog["gates"]["midpoint"]
    assert mid["unavailable"]["derived"] == 1
    assert len(mid["unverified"]) == catalog["gates"]["round_trip"]["templates"] - 1


class FakeRates:
    """Stand-in for ``ProductionRates``: a fixed smooth k(T) on the requested grid."""
    instances = []

    def __init__(self, database, families, *, database_path=None, database_sha=None,
                 temperature_grid=None, **kw):
        self.grid = tuple(temperature_grid)
        FakeRates.instances.append(self)

    def __call__(self, request):
        import math
        from types import SimpleNamespace
        return SimpleNamespace(problems=[], temperatures=self.grid, ssa_multiplier=1.0,
                               k=[1.0e7 * math.exp(-400.0 / t) for t in self.grid])


def _one_row_with_maps():
    return json.loads((HERE / "fixtures" / "i155_row_with_maps.json").read_text())


def test_midpoint_equals_the_harness_midpoint_table(monkeypatch):
    acc = _harness()
    row = _one_row_with_maps()
    payload = row["payload"]
    record = {"event_id": payload["record_event_id"], "family": "R_Recombination",
              "template": "T", "orientation": "as_generated"}
    monkeypatch.setattr(te, "ProductionRates", FakeRates)
    monkeypatch.setitem(acc._S, "by_id", {record["event_id"]: record})
    monkeypatch.setitem(acc._S, "db", None)
    monkeypatch.setitem(acc._S, "records", [record])
    table = acc.midpoint_table([payload])
    entries = emit.midpoint_entries(record, payload, emit.RatesByGrid(None, ["R_Recombination"], "/db"))
    (entry,) = entries.values()
    assert table["templates"] == 1 and entry["status"] == "ok"
    assert entry["worst"] == table["worst"] and entry["worst_T"] == table["worst_T"]
    assert table["worst"] > 0.0


def test_midpoint_without_a_witness_is_unavailable_not_skipped(monkeypatch):
    monkeypatch.setattr(te, "ProductionRates", FakeRates)
    payload = copy.deepcopy(_one_row_with_maps()["payload"])
    for t in payload["templates"]:
        t["maps"] = {}
    record = {"event_id": payload["record_event_id"], "family": "R_Recombination",
              "template": "T", "orientation": "as_generated"}
    entries = emit.midpoint_entries(record, payload, emit.RatesByGrid(None, ["R_Recombination"], "/db"))
    (entry,) = entries.values()
    assert entry["status"] == "unavailable" and "witness" in entry["reason"]


def test_midpoint_regenerates_on_the_templates_own_generation_side(monkeypatch):
    """The emitted witness is stored in template direction; the rate is
    regenerated on the side its member was generated on, which is the
    template's own orientation (members of both orientations mix in one
    record), not one flag for the whole payload."""
    seen = []

    class Recording(FakeRates):
        def __call__(self, request):
            seen.append(request)
            return super().__call__(request)
    monkeypatch.setattr(te, "ProductionRates", Recording)
    payload = copy.deepcopy(_one_row_with_maps()["payload"])
    payload.pop("member_reverse", None)
    record = {"event_id": payload["record_event_id"], "family": "R_Recombination",
              "template": "T", "orientation": "as_generated"}
    for flag in (False, True):
        seen.clear()
        for t in payload["templates"]:
            t["member_reverse"] = flag
        entries = emit.midpoint_entries(record, payload,
                                        emit.RatesByGrid(None, ["R_Recombination"], "/db"))
        (t,) = [t for t in payload["templates"] if t["key"] in entries]
        wit = t["maps"]["witnesses"]["emitted"]
        r_rows, p_rows = te.rows_from_json(wit["r"]), te.rows_from_json(wit["p"])
        assert entries[t["key"]]["status"] == "ok"
        assert [dict(seen[0].generation), dict(seen[0].products)] == \
            ([p_rows, r_rows] if flag else [r_rows, p_rows])


def test_midpoint_rate_failure_is_unavailable_with_the_reason(monkeypatch):
    class Failing(FakeRates):
        def __call__(self, request):
            raise te.EmissionError("no rate")
    monkeypatch.setattr(te, "ProductionRates", Failing)
    payload = _one_row_with_maps()["payload"]
    record = {"event_id": payload["record_event_id"], "family": "R_Recombination",
              "template": "T", "orientation": "as_generated"}
    (entry,) = emit.midpoint_entries(record, payload, emit.RatesByGrid(None, [], "/db")).values()
    assert entry["status"] == "unavailable" and "no rate" in entry["reason"]


# --------------------------------------------------------------------------
# rework 1: coverage is pinned into the cycle check
# --------------------------------------------------------------------------
def test_a_missing_shard_alone_makes_the_cycle_check_partial(tmp_path, monkeypatch):
    """Every record is assigned and has a row, the artifact count matches, no
    header disagrees: the only defect is that shard 1/2 never arrived."""
    seen = []

    def cycles(payloads, counter, kc, complete=True, **kw):
        seen.append(complete)
        status = "ok" if complete else "not_evaluable"
        return payloads, [{"forward": "evt_" + "a" * 20, "reverse": "evt_" + "b" * 20,
                           "status": status, "reason": "reverse record absent"}]
    monkeypatch.setattr(te, "check_cycles", cycles)
    rows = fixture_rows()
    path = write_rows(tmp_path / "only0.jsonl", header("0/2", rows), rows)
    catalog = run([path])
    assert seen == [False]
    assert catalog["gates"]["cycles"]["by_status"] == {"not_evaluable": 1}
    assert catalog["gates"]["coverage"]["complete"] is False
    assert catalog["complete"] is False
    assert any("missing shard 1/2" in f for f in catalog["gates"]["failures"])
    full = run(shards(tmp_path, 2))
    assert seen[-1] is True and full["gates"]["cycles"]["by_status"] == {"ok": 1}


def test_the_catalog_lists_every_key_that_entered_collision_resolution(tmp_path):
    """I-162 S-g rework 1: a catalog set compares its members' entry keys, so
    a catalog lists every executable key that entered collision resolution,
    the refused ones included (a key refused here still collides with the same
    key in another family's catalog under one joint assemble)."""
    rows = fixture_rows()
    emitted = [r for r in rows if r["status"] == "emitted"]
    donor, host = emitted[0], emitted[1]
    clone = copy.deepcopy(donor["payload"]["templates"][0])
    clone["weight"] = dict(clone["weight"], w=[w * 1.0e3 for w in clone["weight"]["w"]])
    host["payload"]["templates"] = host["payload"]["templates"] + [clone]
    host["payload"] = te.finish(host["payload"])
    want = sorted({t["key"] for r in rows for t in (r.get("payload") or {}).get("templates", [])
                   if t.get("status") == "emitted"})
    catalog = run([write_rows(tmp_path / "c.jsonl", header(None, rows), rows)])
    assert catalog["gates"]["collisions"]["refused"] >= 1
    assert clone["key"] in want
    assert catalog["gates"]["collisions"]["entry_keys"] == want
    loaded = {t["key"] for p in catalog["payloads"] for t in p["templates"]}
    assert clone["key"] not in loaded


def test_the_catalog_carries_the_content_identities_of_its_inputs(tmp_path):
    catalog = run(shards(tmp_path, database_content_sha256="h" * 64, code_sources_sha256="s" * 64))
    assert catalog["provenance"]["database_content_sha256"] == "h" * 64
    assert catalog["provenance"]["code_sources_sha256"] == "s" * 64
    mixed = shards(tmp_path, n=2, code_sources_sha256="s" * 64)
    head = json.loads(mixed[1].read_text().splitlines()[0])
    head["code_sources_sha256"] = "t" * 64
    lines = mixed[1].read_text().splitlines()
    mixed[1].write_text("\n".join([json.dumps(head)] + lines[1:]) + "\n")
    with pytest.raises(emit.AssembleRefused, match="code_sources_sha256"):
        run(mixed)


# --------------------------------------------------------------------------
# I-162 S-g rework 2: assemble derives with its own code and database
# --------------------------------------------------------------------------
@pytest.mark.parametrize("field_", ["database_content_sha256", "code_sources_sha256",
                                    "code_extensions_sha256"])
def test_assemble_refuses_shards_emitted_with_other_code_or_database(tmp_path, field_,
                                                                     assembler_identity):
    """Derivation and the cycle check run with the assembler's own code and
    database: shards emitted with others are refused by field name, and the
    catalog would otherwise carry their digests."""
    assembler_identity[field_] = "z" * 64
    with pytest.raises(emit.AssembleRefused) as caught:
        run(shards(tmp_path))
    assert field_ in str(caught.value) and "z" * 16 in str(caught.value)


@pytest.mark.parametrize("field_", ["database_content_sha256", "code_sources_sha256",
                                    "code_extensions_sha256"])
def test_assemble_refuses_shards_without_the_digests(tmp_path, field_):
    with pytest.raises(emit.AssembleRefused, match=field_):
        run(shards(tmp_path, **{field_: None}))


def test_cli_refusal_writes_no_catalog_and_exits_nonzero(tmp_path, monkeypatch, assembler_identity):
    monkeypatch.setattr(emit, "load_services", lambda database, families: Services())
    assembler_identity["code_sources_sha256"] = "x" * 64
    out = tmp_path / "refused.json"
    assert emit.main(["assemble", *[str(p) for p in shards(tmp_path)], "--out", str(out)]) == 1
    assert not out.exists()


def test_the_code_digest_covers_every_rmgpy_source(tmp_path):
    """Every .py, .pyx and .pxd file under rmgpy/ (paths included), not only
    rmgpy/kmc: a change to rmgpy/reaction.py changes the chemistry."""
    root = tmp_path / "rmgpy"
    for rel, text in {"reaction.py": "a = 1\n", "kmc/emit.py": "b = 1\n", "molecule/x.pyx": "c\n",
                      "molecule/x.pxd": "d\n", "molecule/x.c": "generated\n"}.items():
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        (root / rel).write_text(text)
    base = emit.code_sources_sha256(root)
    for rel, changes in (("reaction.py", True), ("molecule/x.pyx", True), ("molecule/x.pxd", True),
                         ("kmc/emit.py", True), ("molecule/x.c", False)):
        before = (root / rel).read_text()
        (root / rel).write_text(before + "# changed\n")
        assert (emit.code_sources_sha256(root) != base) is changes, rel
        (root / rel).write_text(before)
    (root / "reaction.py").rename(root / "reaction2.py")              # the paths are digested
    assert emit.code_sources_sha256(root) != base
    import rmgpy
    assert emit.code_sources_sha256() == emit.code_sources_sha256(Path(rmgpy.__file__).parent)



# --------------------------------------------------------------------------
# I-162 S-g rework 3: the compiled extensions and the supplied services
# --------------------------------------------------------------------------
def test_the_extension_digest_reads_the_compiled_modules(tmp_path):
    """Edited sources next to extensions built from another revision: the
    source digest moves, the binaries do not, so the binaries are digested
    too (every extension module file under rmgpy/, paths included); a loaded
    rmgpy extension from outside that tree makes the digest unknown."""
    from types import SimpleNamespace
    root = tmp_path / "rmgpy"
    (root / "molecule").mkdir(parents=True)
    so = root / "molecule" / "molecule.cpython-39-x86_64-linux-gnu.so"
    so.write_bytes(b"old build")
    (root / "molecule" / "molecule.pyx").write_text("cdef int a\n")
    loaded = [SimpleNamespace(__name__="rmgpy.molecule.molecule", __file__=str(so))]
    base = emit.code_extensions_sha256(root, modules=loaded)
    assert base and len(base) == 64
    (root / "molecule" / "molecule.pyx").write_text("cdef int b\n")     # a source edit alone
    assert emit.code_extensions_sha256(root, modules=loaded) == base
    so.write_bytes(b"new build")
    assert emit.code_extensions_sha256(root, modules=loaded) != base
    stray = tmp_path / "elsewhere" / "kinetics.cpython-39-x86_64-linux-gnu.so"
    stray.parent.mkdir()
    stray.write_bytes(b"other")
    outside = loaded + [SimpleNamespace(__name__="rmgpy.kinetics.kinetics", __file__=str(stray))]
    assert emit.code_extensions_sha256(root, modules=outside) is None
    assert emit.code_extensions_sha256() is not None                 # this process, this build


def test_the_identity_of_this_process_covers_the_extensions():
    """The shard header and the assembler compare the same three digests."""
    assert set(emit.IDENTITY_FIELDS) == {"database_content_sha256", "code_sources_sha256",
                                         "code_extensions_sha256"}
    assert set(emit.IDENTITY_FIELDS) <= set(emit.HEADER_MATCH)


@pytest.mark.parametrize("digest", ["y" * 64, None])
def test_assemble_refuses_services_built_from_another_database(tmp_path, digest):
    """Services carry the digest of the database they were built from; the
    balance check would otherwise run on database B under provenance A."""
    services = Services()
    services.database_content_sha256 = digest
    with pytest.raises(emit.AssembleRefused, match="services"):
        emit.assemble([str(p) for p in shards(tmp_path)], services=services)


def test_load_services_records_its_database_digest(monkeypatch):
    from types import SimpleNamespace
    monkeypatch.setattr(emit, "load_database", lambda families, database: SimpleNamespace(
        kinetics=SimpleNamespace(families={}), thermo=None))
    monkeypatch.setattr(emit, "database_content_sha256", lambda database: "d" * 64)
    monkeypatch.setattr(emit.te, "ThermoKc", lambda provider: None)
    import rmgpy.kmc.reference_thermo as rt
    monkeypatch.setattr(rt, "GasPhaseRMGReferenceThermo", lambda thermo, label: None)
    assert emit.load_services("/db", []).database_content_sha256 == "d" * 64


# --------------------------------------------------------------------------
# I-188: a catalog with partial payloads merges; a raising gate refuses it
# --------------------------------------------------------------------------
def i188_rows():
    """Rows of four records of a full emit (``i188_reciprocal_partial_payloads``):
    two reciprocal record pairs sharing two base templates, two records
    partial with a base template refused at emission (no key)."""
    payloads = json.loads((HERE / "fixtures" / "i188_reciprocal_partial_payloads.json").read_text())
    rows = [{"event_id": p["record_event_id"][:24], "family": p["family"],
             "orientation": p["orientation"], "status": p["status"], "refusal": None,
             "payload": p} for p in payloads]
    for row in rows:
        row["midpoint"] = ok_midpoint(row)
    return rows


def test_partial_payloads_with_keyless_refused_templates_assemble_merged(tmp_path, balance_stubs):
    rows = i188_rows()
    assert any(t["status"] == "refused" and "key" not in t
               for r in rows for t in r["payload"]["templates"])
    catalog = run(shards(tmp_path, n=1, rows=rows, artifact_records=len(rows)))
    failures = catalog["gates"]["failures"]
    assert not [f for f in failures if f.startswith("catalog consistency")], failures
    collisions = catalog["gates"]["collisions"]
    assert collisions["keys"] == collisions["merged"] == 2 and collisions["refused"] == 0
    keys = [t["key"] for p in catalog["payloads"] for t in p["templates"]]
    assert keys and len(keys) == len(set(keys))
    assert balance_stubs["derive"] == 1


def test_a_raising_collision_resolution_refuses_the_catalog(tmp_path, monkeypatch, balance_stubs):
    """No silent fallback to the unmerged payloads: the error is named and
    nothing loads (no balance check runs on unmerged copies)."""
    def split(payloads):
        raise te.ReciprocalSplitError("reciprocal split: test")
    monkeypatch.setattr(te, "resolve_collisions", split)
    catalog = run(shards(tmp_path))
    failures = catalog["gates"]["failures"]
    assert any(f.startswith("catalog consistency: ReciprocalSplitError: reciprocal split: test")
               for f in failures), failures
    assert catalog["payloads"] == [] and catalog["complete"] is False
    assert balance_stubs["derive"] == balance_stubs["cycles"] == 0


# -- I-195: an artifact compiled from another declaration is never emitted against this roster


class _Reached(Exception):
    """Raised by the first expensive step: the binding let the run through."""


def _digests(kind):
    from rmgpy.kmc.run import chemistry_sha256, declaration_sha256
    from rmgpy.polymer_spec import normalized_linear_homopolymer_spec
    spec = normalized_linear_homopolymer_spec(polymer_mapping(kind))
    return {"chemistry_sha256": chemistry_sha256(spec), "declaration_sha256": declaration_sha256(spec)}


def _bound_artifact(tmp_path, proxy_source):
    """An artifact whose provenance carries ``proxy_source`` (None: a legacy
    compile without one)."""
    provenance = {} if proxy_source is None else {"proxy_source": proxy_source}
    path = tmp_path / "artifact.json"
    path.write_text(json.dumps({"provenance": provenance, "records": []}))
    return path


def _declared_source(kind):
    return {"kind": "declaration", "label": kind, "deck_sha256": "deck", **_digests(kind),
            "proxy_units": 3, "proxy_set": "declared"}


def _emit_cli(tmp_path, monkeypatch, artifact, *flags, roster=None):
    def reached(_roster):
        raise _Reached
    roster = roster or write_roster(tmp_path, "pe", {"a": [chain(2)]})
    monkeypatch.setattr(emit, "build_reference_systems", reached)
    return emit.main([str(artifact), str(roster), "--database", "db",
                      "--out", str(tmp_path / "rows.jsonl"), *flags])


def test_a_roster_carries_the_digests_of_the_declaration_it_was_built_from(tmp_path):
    from meltRunTest import pe_deck
    from rmgpy.kmc.run import read_deck

    deck = read_deck(pe_deck(tmp_path))
    generated = emit.reference_systems(deck.spec, emit.match_provider(deck.spec.repeat_unit),
                                       dps=(1, 2))
    written = emit.write_roster(generated, tmp_path / "roster.json")
    assert generated.digests() == written.digests() == deck.digests()
    assert emit.load_roster(write_roster(tmp_path, "pe", {"a": [chain(2)]})).digests() == _digests("pe")


def test_emit_refuses_a_roster_of_another_declaration_before_any_expensive_step(
        tmp_path, monkeypatch, capsys):
    artifact = _bound_artifact(tmp_path, _declared_source("ps"))
    assert _emit_cli(tmp_path, monkeypatch, artifact) == 1
    out = capsys.readouterr().out
    assert "refused:declaration" in out and "declaration_sha256" in out
    assert not (tmp_path / "rows.jsonl").exists()


def test_emit_refuses_another_declaration_even_with_the_legacy_flag(tmp_path, monkeypatch, capsys):
    artifact = _bound_artifact(tmp_path, _declared_source("ps"))
    assert _emit_cli(tmp_path, monkeypatch, artifact, "--legacy-artifact") == 1
    assert "refused:declaration" in capsys.readouterr().out


def _v3_roster(tmp_path, kind, units):
    from rmgpy.polymer_spec import normalized_linear_homopolymer_spec
    spec = normalized_linear_homopolymer_spec(polymer_mapping(kind))
    path = tmp_path / f"{kind}.v3.json"
    emit.write_roster(emit.declared_reference_systems(spec, units), path)
    return path


def test_emit_accepts_the_roster_of_the_artifacts_declaration(tmp_path, monkeypatch):
    artifact = _bound_artifact(tmp_path, _declared_source("pe"))
    with pytest.raises(_Reached):
        _emit_cli(tmp_path, monkeypatch, artifact, roster=_v3_roster(tmp_path, "pe", 3))


@pytest.mark.parametrize("units", (None, 2))
def test_a_declared_artifact_refuses_a_roster_not_built_from_its_compile_structures(
        tmp_path, monkeypatch, capsys, units):
    """A v2 roster, or the v3 roster of other proxy units, is refused by
    name before anything is built."""
    artifact = _bound_artifact(tmp_path, _declared_source("pe"))
    roster = None if units is None else _v3_roster(tmp_path, "pe", units)
    assert _emit_cli(tmp_path, monkeypatch, artifact, roster=roster) == 1
    out = capsys.readouterr().out
    assert "refused:roster" in out and "kmc-reference-roster/3" in out


def test_a_legacy_artifact_refuses_a_v3_roster(tmp_path, monkeypatch, capsys):
    artifact = _bound_artifact(tmp_path, None)
    assert _emit_cli(tmp_path, monkeypatch, artifact, "--legacy-artifact",
                     roster=_v3_roster(tmp_path, "pe", 2)) == 1
    assert "refused:roster" in capsys.readouterr().out


@pytest.mark.parametrize("source", [None, {"kind": "built-in", "proxy_set": "built-in",
                                           "proxy_units": None, "reason": "no deck was given"}])
def test_emit_refuses_a_legacy_artifact_unless_flagged(tmp_path, monkeypatch, capsys, source):
    artifact = _bound_artifact(tmp_path, source)
    assert _emit_cli(tmp_path, monkeypatch, artifact) == 1
    out = capsys.readouterr().out
    assert "refused:unbound_artifact" in out and "--legacy-artifact" in out
    with pytest.raises(_Reached):
        _emit_cli(tmp_path, monkeypatch, artifact, "--legacy-artifact")


@pytest.mark.parametrize("key", ["chemistry_sha256", "declaration_sha256"])
def test_each_carried_digest_binds(key):
    source = dict(_declared_source("pe"), **{key: "0" * 64})
    with pytest.raises(emit.ArtifactBindingRefusal, match=rf"refused:declaration: .*{key}"):
        emit.bind_artifact({"provenance": {"proxy_source": source}}, _digests("pe"), "the roster")
    assert emit.bind_artifact({"provenance": {"proxy_source": _declared_source("pe")}},
                              _digests("pe"), "the roster") == ("chemistry_sha256",
                                                                "declaration_sha256")


# -- the catalog carries the declaration it was emitted for, so the runner binds it to the deck


class _StubEmitter:
    """An emitter with no records: enough for the shard header."""

    def __init__(self, artifact, database, systems, **kwargs):
        self.artifact, self.database, self.systems = Path(artifact), database, systems
        self.records, self.state = [], {}

    def init(self):
        self.state.update(artifact_records=0, families=["R_Recombination"])


def _emitted_header(tmp_path, monkeypatch, artifact, *flags, roster):
    monkeypatch.setattr(emit, "build_reference_systems", lambda roster: ({}, []))
    monkeypatch.setattr(emit, "Emitter", _StubEmitter)
    monkeypatch.setattr(emit, "database_content_sha256", lambda path: "h" * 64)
    out = tmp_path / "rows.jsonl"
    assert emit.main([str(artifact), str(roster), "--database", str(tmp_path),
                      "--out", str(out), *flags]) == 0
    return json.loads(out.read_text().splitlines()[0])


def test_a_bound_shard_header_carries_the_chemistry_digest_only(tmp_path, monkeypatch):
    artifact = _bound_artifact(tmp_path, _declared_source("pe"))
    head = _emitted_header(tmp_path, monkeypatch, artifact, roster=_v3_roster(tmp_path, "pe", 3))
    assert head["chemistry_sha256"] == _digests("pe")["chemistry_sha256"]
    assert "declaration_sha256" not in head


def test_a_legacy_shard_header_carries_no_declaration_digest(tmp_path, monkeypatch):
    artifact = _bound_artifact(tmp_path, None)
    head = _emitted_header(tmp_path, monkeypatch, artifact, "--legacy-artifact",
                           roster=write_roster(tmp_path, "pe", {"a": [chain(2)]}))
    assert not set(emit.BINDING_KEYS) & set(head)


def test_the_catalog_of_bound_shards_binds_to_the_deck_of_its_declaration(tmp_path):
    from meltRunTest import pe_deck
    from rmgpy.kmc import run as kmc_run

    deck = kmc_run.read_deck(pe_deck(tmp_path))
    catalog = run(shards(tmp_path, chemistry_sha256=deck.digests()["chemistry_sha256"]))
    assert catalog["provenance"]["chemistry_sha256"] == deck.digests()["chemistry_sha256"]
    assert catalog["complete"] is True
    assert kmc_run._bound_by("catalog.json", catalog["provenance"], deck, False) == \
        ("chemistry_sha256",)


def test_the_catalog_of_unbound_shards_carries_no_declaration_digest(tmp_path):
    catalog = run(shards(tmp_path))
    assert not set(emit.BINDING_KEYS) & set(catalog["provenance"])


def test_shards_bound_to_different_declarations_are_a_listed_failure(tmp_path):
    paths = shards(tmp_path, chemistry_sha256=_digests("pe")["chemistry_sha256"])
    other = json.loads(paths[1].read_text().splitlines()[0])
    other["chemistry_sha256"] = _digests("ps")["chemistry_sha256"]
    lines = paths[1].read_text().splitlines()
    paths[1].write_text("\n".join([json.dumps(other)] + lines[1:]) + "\n")
    catalog = run(paths)
    assert catalog["complete"] is False
    assert any("chemistry_sha256" in failure for failure in catalog["gates"]["failures"])


def test_the_catalog_binds_to_every_deck_of_its_chemistry(tmp_path, monkeypatch):
    """Emission reads only the declaration's chemistry, so a catalog emitted
    for one population and transport block runs on a deck that changes them."""
    from meltRunTest import S, pe_deck, write_deck
    from rmgpy.kmc import run as kmc_run

    emitted_for = kmc_run.read_deck(pe_deck(tmp_path))
    other = S.pe_mapping({30: 4}, d_ref=3.0e-17)
    other["initial_mass_kg"], other["density_kg_m3"] = 2.0, 900.0
    deck = kmc_run.read_deck(write_deck(tmp_path / "other_deck.py", other))
    assert deck.digests()["declaration_sha256"] != emitted_for.digests()["declaration_sha256"]
    artifact = _bound_artifact(tmp_path, dict(_declared_source("pe"), **emitted_for.digests()))
    roster = emit.write_roster(emit.declared_reference_systems(emitted_for.spec, 3),
                               tmp_path / "roster.json")
    head = _emitted_header(tmp_path, monkeypatch, artifact, roster=roster.path)
    catalog = run(shards(tmp_path, **{key: head[key] for key in emit.BINDING_KEYS if key in head}))
    assert kmc_run._bound_by("catalog.json", catalog["provenance"], deck, False) == \
        ("chemistry_sha256",)
