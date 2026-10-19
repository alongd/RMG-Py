"""The compile proxies, pair classes and roster of a linear copolymer
declaration (schema v2): the residue-sequence chain provider
(``declared_sequence``), the structures and proxies built on the witness
words of ``sequence_environments`` (``declared_proxies.sequence_*``), pairs
over the family-tree classes of each participant (``pair_classes``) and the
v4 roster built from the same structures (``emit``).

The family-tree tests need the RMG database (``rmg_moves_harness.DATABASE``)."""

import importlib.util
import itertools
import json
from dataclasses import asdict, replace
from pathlib import Path
from types import SimpleNamespace

import pytest

import rmg_moves_harness as harness
from rmgpy.kmc import boundaries as b
from rmgpy.kmc import declared_environments as de
from rmgpy.kmc import emit
from rmgpy.kmc import windows as win
from rmgpy.kmc.declared_chain import ChainSpecError, EditOp, site_class_policy
from rmgpy.kmc.declared_proxies import (
    sequence_declared_structures,
    sequence_proxy_set,
    sequence_structures,
)
from rmgpy.kmc.declared_sequence import SequenceChainSpec, sequence_provider
from rmgpy.kmc.pair_classes import class_pairs, tree_fingerprint
from rmgpy.kmc.run import read_deck
from rmgpy.polymer_spec import (
    LinkDeclaration,
    ResidueDeclaration,
    SchemaV2Refusal,
    linear_junctions,
)

from polymer_spec_v2_support import DECKS, q, v2_spec

EP_DECK = DECKS / "ep_v2_deck.py"
SCRIPT = Path(__file__).parents[3] / "scripts" / "check-kmc-reference-mapping.py"
needs_database = pytest.mark.skipif(not harness.DATABASE.is_dir(),
                                    reason=f"needs the RMG database at {harness.DATABASE}")


@pytest.fixture(scope="module")
def ep():
    return read_deck(EP_DECK).spec


@pytest.fixture(scope="module")
def provider(ep):
    return sequence_provider(ep)


@pytest.fixture(scope="module")
def structures(ep):
    return sequence_structures(ep, 3)


def _heavy(molecule):
    return [atom for atom in molecule.atoms if not atom.is_hydrogen()]


def _key(spec, build, unit, key):
    """The site-class key of the atom (unit, key) of a build, under the
    environment set's size-independent policy."""
    from rmgpy.kmc.site_class import extract_site_class

    residue = dict(spec.residues)[build.metadata.residues[unit]]
    position = sorted(atom.key for atom in residue.graph.atoms).index(key)
    atom_id = next(a for a, u, p in build.metadata.coordinates if (u, p) == (unit, position))
    atom = next(atom for atom in build.molecule.atoms if atom.id == atom_id)
    policy = site_class_policy(spec.locality_policy, de.MAX_HEAVY_ATOMS)
    return extract_site_class(build.molecule, [("site", atom)], policy)[0]


# -- the residue-sequence chain provider ------------------------------------------


def test_a_sequence_chain_builds_its_residues_in_order(provider):
    build = provider.build("p", SequenceChainSpec(("E", "P", "E")))
    assert len(_heavy(build.molecule)) == 2 + 3 + 2
    assert build.molecule.get_formula() == "C7H16"
    assert build.metadata.residues == ("E", "P", "E")
    assert build.metadata.provider == provider.name
    units = sorted({(u, p) for _a, u, p in build.metadata.coordinates})
    assert units == [(0, 0), (0, 1), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1)]
    ports = {p.port_id: p for p in build.metadata.ports}
    assert {p.tag for p in ports.values()} == {"backbone_head", "backbone_tail"}
    assert all(p.artificial and len(p.cap_atom_ids) == 1 for p in ports.values())
    # the head port sits on the first residue's entry atom, the tail on the last's exit
    assert build.metadata.coordinate(ports["head"].atom_id) == (0, 0)
    assert build.metadata.coordinate(ports["tail"].atom_id) == (2, 1)


def test_a_step_outside_the_support_is_refused():
    spec = v2_spec()
    support = dict(spec.sequence_support)
    support["E"] = ("N", "P")
    restricted = sequence_provider(replace(spec, sequence_support=tuple(sorted(support.items()))))
    restricted.build("ok", SequenceChainSpec(("E", "P", "E")))
    with pytest.raises(ChainSpecError, match="not in the sequence support"):
        restricted.build("bad", SequenceChainSpec(("E", "E")))
    with pytest.raises(ChainSpecError, match="not a declared residue"):
        restricted.build("bad", SequenceChainSpec(("E", "X")))


def test_edits_address_a_residue_atom_by_unit_and_key(provider):
    chain = SequenceChainSpec(("E", "P", "E"), edits=(
        EditOp("remove_h", 1, "p2"), EditOp("set_radical", 1, "p2", 1)))
    build = provider.build("r", chain)
    radicals = [a for a in build.molecule.atoms if a.radical_electrons]
    assert len(radicals) == 1
    assert build.metadata.coordinate(radicals[0].id) == (1, 1)
    with pytest.raises(ChainSpecError):
        provider.build("x", replace(chain, edits=(EditOp("remove_h", 1, "e1"),)))


def test_radical_and_declared_ends(provider):
    radical = provider.build("h", SequenceChainSpec(("P", "E"), head="radical"))
    head = radical.metadata.port("head")
    assert not head.artificial and not head.cap_atom_ids
    assert next(a for a in radical.molecule.atoms if a.id == head.atom_id).radical_electrons == 1
    declared = provider.build("d", SequenceChainSpec(("P", "E"), tail="declared"))
    assert not declared.metadata.port("tail").artificial
    assert declared.molecule.get_formula() == "C5H12"


def test_a_sequence_chain_is_not_extended_past_an_artificial_port(provider):
    port, coordinate = provider.probe_ports()[0]
    with pytest.raises(b.ExtensionRefused, match="refused:sequence_extension"):
        provider.fragment(port, coordinate)
    spec = win.DeclaredSpec(provider.name, SequenceChainSpec(("E", "P")))
    with pytest.raises(b.ExtensionRefused, match="refused:sequence_extension"):
        win.longer_specs([spec], [(0, "tail")])


def test_the_windows_read_a_sequence_chain_as_a_declared_chain(provider):
    spec = win.DeclaredSpec(provider.name, SequenceChainSpec(("E", "P", "E")))
    build = win.build_spec(spec, "p")
    _molecules, participants = b.witness_reactants([build])
    policy = win.policy_for(SimpleNamespace(participants=tuple(participants)))
    assert policy.name == win.declared_interim({"C", "H"}).name
    with pytest.raises(win.WindowError, match="refused:legacy_spec"):
        win.witness_spec(win.DeclaredSpec(provider.name, object()))


def test_the_growth_unit_is_the_common_formula_of_the_residues(provider):
    assert provider.unit_formula() == (("C", 1), ("H", 2))
    with pytest.raises(ValueError, match="refused:growth_unit"):
        sequence_provider(v2_spec()).unit_formula()


# -- the compile structures --------------------------------------------------------


def _homopolymer_keys(spec):
    keys = set()
    for name, _residue in spec.residues:
        keys |= set(de.period_environments(dict(spec.residues), spec.junctions,
                                           spec.locality_policy, (name,)))
    return keys


def test_every_h_bearing_environment_has_an_interior_radical_junctions_included(ep, structures):
    environments = de.sequence_environments(ep)
    residues = dict(ep.residues)
    h_bearing = {key for key, occurrences in environments.items()
                 if any(atom.hydrogens for atom in residues[occurrences[0][1]].graph.atoms
                        if atom.key == occurrences[0][2])}
    radicals = {name: s for name, s in structures.items() if name.startswith("interior_radical:")}
    assert len(h_bearing) == 18
    assert {s.environment for s in radicals.values()} == h_bearing
    assert len(radicals) == 18
    junctions = {s.environment for s in radicals.values() if s.junction}
    assert junctions == h_bearing - _homopolymer_keys(ep) and len(junctions) == 14


def test_each_interior_radical_realizes_its_class_on_its_witness_word(ep, structures):
    k = de.spec_word_order(ep)
    for name, s in structures.items():
        if not name.startswith("interior_radical:"):
            continue
        chain = s.witness.chain
        edit = next(e for e in chain.edits if e.kind == "set_radical")
        pristine = sequence_provider(ep).build("p", replace(chain, edits=()))
        assert _key(ep, pristine, edit.unit, edit.key) == s.environment, name
        # the shortest witness word, the radical's residue on the centre residue
        assert len(s.word) <= k and len(chain.residues) >= len(s.word)
        left, right = edit.unit, len(chain.residues) - 1 - edit.unit
        assert abs(left - right) <= 1, name


def test_the_pristine_chains_are_the_witness_chains_unedited(ep, structures):
    pristine = [s for name, s in structures.items() if name.startswith("pristine:")]
    chains = {s.witness.chain for s in pristine}
    assert all(not chain.edits for chain in chains)
    words = {replace(s.witness.chain, edits=()).residues for name, s in structures.items()
             if name.startswith("interior_radical:")}
    assert {chain.residues for chain in chains} <= words and pristine


def test_end_radicals_one_per_end_environment(ep, structures):
    for side in ("head", "tail"):
        ends = {name: s for name, s in structures.items()
                if name.startswith(f"{side}_end_radical:")}
        assert ends, side
        assert len({s.environment for s in ends.values()}) == len(ends)
        for s in ends.values():
            assert getattr(s.witness.chain, side) == "radical"


def test_one_monomer_per_monomer_product(ep):
    assert not [name for name in sequence_structures(ep, 3) if name.startswith("monomer")]

    def monomers(mapping):
        del mapping["residues"]["N"]
        mapping["sequence_support"] = {"E": ["E", "P"], "P": ["E", "P"]}
        mapping["initial_ensemble"]["composition"] = {"E": q(1, 2), "P": q(1, 2)}
        def graph(hydrogens, orders):
            keys = "abc"[:len(hydrogens)]
            return {"atoms": [{"key": k, "element": "C", "hydrogens": h,
                               "electronic_state": {"radical_electrons": 0}}
                              for k, h in zip(keys, hydrogens)],
                    "bonds": [{"atoms": [x, y], "order": o}
                              for x, y, o in zip(keys, keys[1:], orders)],
                    "electronic_state": {"multiplicity": 1}}
        mapping["monomer_products"] = {"E": graph((2, 2), (2,)), "P": graph((2, 1, 3), (2, 1))}
    found = sequence_structures(v2_spec(monomers), 3)
    names = sorted(name for name in found if name.startswith("monomer"))
    assert names == ["monomer:E", "monomer:P"]
    assert found["monomer:P"].build.molecule.get_formula() == "C3H6"
    assert found["monomer:P"].witness is None and not found["monomer:P"].build.metadata.ports


def test_the_full_context_copies_hold_one_more_residue_at_each_artificial_end(ep):
    structures, _proxies = sequence_declared_structures(ep, 3, {"3": [], "5": []})
    for name, (witness, _build) in structures.items():
        if not name.endswith("@3") or witness is None:
            continue
        wide = structures[name[:-2] + "@5"][0]
        artificial = sum(getattr(witness.chain, side) == "artificial" for side in ("head", "tail"))
        assert len(wide.chain.residues) == len(witness.chain.residues) + artificial, name


def test_the_proxy_units_bound_the_chain_length(ep):
    for units in (1, 6):
        for name, s in sequence_structures(ep, units).items():
            if s.witness is not None:
                assert len(s.witness.chain.residues) >= units, name


# -- pairs over family-tree classes (Q5(a)) ---------------------------------------------


def _lift(path):
    """A v1 deck as a one-residue v2 declaration, built through the internal
    constructor (a v1-expressible deck is refused by the reader)."""
    v1 = read_deck(path).spec
    unit = v1.repeat_unit.as_mapping()
    ports = unit.pop("ports")
    residue = ResidueDeclaration.from_mapping({"graph": unit, "sites": {
        "head": dict(ports["head"], **{"class": "entry", "vacancy": "forbidden"}),
        "tail": dict(ports["tail"], **{"class": "exit", "vacancy": "forbidden"})}})
    links = {"j": LinkDeclaration.from_mapping(
        {"from": "exit", "to": "entry", "bond_order": ports["tail"]["bond_order"]})}
    junctions = linear_junctions({"R": residue}, links, {"R": ("R",)})
    template = read_deck(EP_DECK).spec
    return replace(template, label=v1.label, residues=(("R", residue),), links=tuple(links.items()),
                   sequence_support=(("R", ("R",)),), head_end_group=v1.head_end_group,
                   tail_end_group=v1.tail_end_group, monomer_products=(),
                   locality_policy=v1.locality_policy, junction_table=tuple(junctions.items()))


def _templates(family, first, second):
    return {tuple(sorted(getattr(entry, "label", entry) for entry in reaction.template))
            for reaction in family.generate_reactions([first.copy(deep=True),
                                                       second.copy(deep=True)])}


@needs_database
@pytest.mark.parametrize("family_name", ["R_Recombination", "H_Abstraction"])
def test_pairs_over_tree_classes_reach_every_all_pairs_template_on_pp(family_name):
    """The tripwire: on PP (one residue, where all pairs are cheap) the pair
    proxies over family-tree classes reach exactly the templates of all
    pairs, and every pair of two classes reaches only templates of the pair
    that stands for them."""
    family = harness.families()[family_name]
    lifted = _lift(DECKS / "pp_deck.py")
    molecules = {name: s.build.molecule for name, s in sequence_structures(lifted, 3).items()}
    pairs, classes = class_pairs({family_name: family}, molecules)
    names = list(molecules)
    every = {}
    for first, second in itertools.combinations_with_replacement(names, 2):
        every[(first, second)] = _templates(family, molecules[first], molecules[second])
    reached = set().union(*(every[pair] for pair in pairs))
    assert reached == set().union(*every.values()) and reached
    group = {member: index for index, members in enumerate(classes) for member in members}
    stands = {frozenset((group[a], group[c])): every[(a, c)] for a, c in pairs}
    for (first, second), found in every.items():
        key = frozenset((group[first], group[second]))
        assert found <= stands.get(key, set()), (first, second)
    assert len(pairs) < len(every)


@needs_database
def test_equivalent_sites_share_a_fingerprint_and_distinct_ones_do_not(provider):
    family = harness.families()["R_Recombination"]

    def radical(residues, unit, key):
        return provider.build("r", SequenceChainSpec(residues, edits=(
            EditOp("remove_h", unit, key), EditOp("set_radical", unit, key, 1)))).molecule
    secondary = tree_fingerprint(family, radical(("E", "E", "E"), 1, "e1"))
    assert secondary == tree_fingerprint(family, radical(("E", "E", "E"), 1, "e2"))
    assert secondary != tree_fingerprint(family, radical(("E", "P", "E"), 1, "p2"))
    assert tree_fingerprint(family, provider.build("p", SequenceChainSpec(("E", "P"))).molecule) \
        == frozenset()


# -- the compile -----------------------------------------------------------------------


@needs_database
def test_the_sequence_proxy_set_records_its_pairs(ep):
    families = {"R_Recombination": harness.families()["R_Recombination"]}
    proxies, record = sequence_proxy_set(ep, 3, families)
    site_types = {p.site_type for p in proxies}
    assert set(record["pairs"]) == {"3", "5"} and record["rule"] == "family_tree_classes"
    for depth, suffix in (("3", ""), ("5", "@5")):
        for first, second in record["pairs"][depth]:
            assert f"{first}+{second}{suffix}" in site_types
    singles = {t for t in site_types if "+" not in t and "@" not in t}
    assert {f"interior_radical:{i}" for i in range(18)} <= singles
    for proxy in proxies:
        assert len(proxy.boundary_builds) == len(proxy.reactants)
        assert proxy.repeat_unit is None
        assert all(port.source == "construction" for port in proxy.artificial_boundaries)
        assert proxy.metadata["declaration"].startswith("sequence:")
    # far fewer pairs than all pairs of the radicals
    radicals = [t for t in singles if "radical" in t]
    assert len(record["pairs"]["3"]) < len(radicals) * (len(radicals) + 1) // 2


class _Recorder:
    seen = {}

    @staticmethod
    def discover_family_reactions(kinetics, proxies, candidates, **kwargs):
        _Recorder.seen["proxies"] = proxies
        return [], {}, {}

    def __init__(self, kinetics, proxies, active, **kwargs):
        _Recorder.seen["kwargs"] = kwargs

    def write_artifact(self, directory):
        return Path(directory) / "artifact.json", {}


@needs_database
def test_compile_with_a_v2_deck_builds_its_sequence_proxies(monkeypatch, tmp_path):
    from rmgpy.kmc import compile as compile_module

    families = {"R_Recombination": harness.families()["R_Recombination"]}
    database = SimpleNamespace(
        load_kinetics=lambda *a, **k: None, load_thermo=lambda *a, **k: None,
        kinetics=SimpleNamespace(families=families, generate_reactions_from_families=None),
        thermo=object())
    monkeypatch.setattr(compile_module, "RMGDatabase", lambda: database)
    monkeypatch.setattr(compile_module, "prepare_rate_rules", lambda *a, **k: None)
    monkeypatch.setattr(compile_module, "resolve_database_declaration", lambda sha: "sha")
    monkeypatch.setattr(compile_module, "EventSetCompiler", _Recorder)
    request = compile_module.CompileRequest(
        database_path=str(tmp_path), output_dir=str(tmp_path / "out"),
        families=("R_Recombination",), family_universe=("R_Recombination",), deck=str(EP_DECK),
        proxy_units=2)
    compile_module.compile_artifact(request)
    deck = read_deck(EP_DECK)
    source = _Recorder.seen["kwargs"]["proxy_source"]
    assert source["kind"] == "declaration" and source["proxy_units"] == 2
    assert source["declaration_sha256"] == deck.digests()["declaration_sha256"]
    assert source["proxy_set"].endswith("sequence_proxy_set")
    proxies, record = sequence_proxy_set(deck.spec, 2, families)
    assert source["pairs"] == record["pairs"]
    assert [p.site_type for p in _Recorder.seen["proxies"]] == [p.site_type for p in proxies]


# -- the v4 roster ---------------------------------------------------------------------

PAIRS = {"3": [["interior_radical:0", "interior_radical:1"]],
         "5": [["interior_radical:0", "interior_radical:0"]]}


def _artifact(spec, pairs=PAIRS, units=3):
    from rmgpy.kmc.run import chemistry_sha256, declaration_sha256
    return {"provenance": {"proxy_source": {
        "kind": "declaration", "proxy_units": units, "pairs": pairs,
        "chemistry_sha256": chemistry_sha256(spec), "declaration_sha256": declaration_sha256(spec)}},
        "records": []}


def test_the_v4_roster_is_built_from_the_compile_structures_and_pairs(ep, tmp_path):
    roster = emit.declared_reference_systems(ep, 3, PAIRS)
    assert roster.schema == emit.ROSTER_SCHEMA_V4 and roster.pairs == PAIRS
    structures, proxies = sequence_declared_structures(ep, 3, PAIRS)
    entry = {name: ({"declared_sequence": asdict(witness.chain)} if witness is not None
                    else None) for name, (witness, _build) in structures.items()}
    assert "pair:interior_radical:0+interior_radical:1:left" in roster.systems
    assert "pair:interior_radical:0+interior_radical:0@5:left" in roster.systems
    for site_type, names in proxies.items():
        kind = "single" if len(names) == 1 else "pair"
        system = roster.systems[f"{kind}:{site_type}" + ("" if kind == "single" else ":left")]
        assert system == [entry[name] for name in names], site_type
    loaded = emit.write_roster(roster, tmp_path / "roster.json")
    assert loaded.schema == emit.ROSTER_SCHEMA_V4 and loaded.pairs == PAIRS
    assert loaded.digests() == read_deck(EP_DECK).digests()
    systems, refusals = emit.build_reference_systems(loaded)
    assert not refusals and len(systems) == len(roster.systems)


def test_a_v2_roster_needs_the_compiles_pairs(ep):
    with pytest.raises(emit.RosterRefusal, match="refused:roster"):
        emit.declared_reference_systems(ep, 3)


def test_load_roster_refuses_a_declaration_of_the_other_schema(ep, tmp_path):
    v1 = read_deck(DECKS / "pp_deck.py").spec
    path = tmp_path / "roster.json"
    path.write_text(json.dumps({"schema": emit.ROSTER_SCHEMA_V4, "polymer": v1.as_mapping(),
                                "proxy_units": 3, "pairs": PAIRS, "systems": {}}))
    with pytest.raises(emit.RosterRefusal, match="refused:roster: a kmc-reference-roster/4 roster "
                                                 "carries a schema-v2 declaration"):
        emit.load_roster(path)
    path.write_text(json.dumps({"schema": emit.ROSTER_SCHEMA_V3,
                                "polymer": json.loads(ep.normalized_json),
                                "proxy_units": 3, "systems": {}}))
    with pytest.raises(SchemaV2Refusal, match=r"emit\.load_roster reads v1 declarations only"):
        emit.load_roster(path)


def test_a_v2_artifact_takes_the_v4_roster_of_its_own_pairs(ep):
    artifact = _artifact(ep)
    roster = emit.declared_reference_systems(ep, 3, PAIRS)
    emit.check_roster_for_artifact(roster, artifact)
    with pytest.raises(emit.RosterRefusal, match="pairs"):
        emit.check_roster_for_artifact(roster, _artifact(ep, {"3": [], "5": []}))
    v1 = read_deck(DECKS / "pp_deck.py").spec
    with pytest.raises(emit.RosterRefusal, match="refused:roster"):
        emit.check_roster_for_artifact(emit.declared_reference_systems(v1, 3), artifact)


def test_the_mapping_check_builds_the_v4_roster_from_the_artifact(ep):
    spec = importlib.util.spec_from_file_location("check_kmc_reference_mapping_v2", SCRIPT)
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
    roster = checker.reference_roster(_artifact(ep), ep, (1, 2, 3, 8))
    assert roster.sha256 == emit.declared_reference_systems(ep, 3, PAIRS).sha256


def test_a_v2_artifact_is_certified_under_its_declaration(ep):
    rules = emit.certification_rules(_artifact(ep), ep)
    assert rules["repeat"] == (("C", 1), ("H", 2))
    assert rules["domain"] is win.declared_interim({"C", "H"}).domain


def test_the_declared_provider_of_a_v2_declaration_is_its_sequence_provider(ep, provider):
    assert emit.declared_provider(ep) is provider
    assert provider.name.startswith("sequence:") and b.PROVIDERS[provider.name] is provider
