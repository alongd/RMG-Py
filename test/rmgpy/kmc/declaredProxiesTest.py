"""The compile proxies of a declared polymer (``rmgpy.kmc.declared_proxies``)
and the compile CLI's ``--deck`` path."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

import melt_engine_support as S
from meltRunTest import pe_deck, ps_deck, ps_deck_mapping
from rmgpy.kmc import compile as compile_module
from rmgpy.kmc.compiler import PS_FAMILY_CANDIDATES, EventSetCompiler, ps_proxy_set
from rmgpy.polymer_spec import LinearHomopolymerSpec, normalized_linear_homopolymer_spec


def _spec(mapping):
    return normalized_linear_homopolymer_spec(LinearHomopolymerSpec.from_mapping(mapping))


def _declared(spec, units=3):
    from rmgpy.kmc.declared_proxies import declared_proxy_set
    return declared_proxy_set(spec, units)


def _molecule(participant):
    return participant.molecule[0] if hasattr(participant, "molecule") else participant


def _ports(build):
    """(tag, artificial, cap count) of a participant build's ports."""
    return sorted((p.tag, p.artificial, len(p.cap_atom_ids)) for p in build.metadata.ports)


# -- a declared polymer's own proxies ------------------------------------------


def test_the_pe_declaration_gives_pe_proxies_without_an_aromatic_ring():
    spec = _spec(S.pe_mapping({10: 1}))
    proxies = _declared(spec)
    names = {p.site_type for p in proxies}
    assert {"pristine", "interior_radical", "doubly_featured", "head_end_radical",
            "head_end_radical+pristine", "head_end_radical+monomer",
            "pristine@5", "head_end_radical+head_end_radical@5"} <= names
    # the PE chain's two end radicals are one structure; it has no side group
    assert not any("tail_end_radical" in name or "pendant" in name for name in names)
    for proxy in proxies:
        assert len(proxy.boundary_builds) == len(proxy.reactants)
        for participant in proxy.reactants:
            molecule = _molecule(participant)
            assert {atom.element.symbol for atom in molecule.atoms} <= {"C", "H"}
            assert not any(molecule.is_atom_in_cycle(atom) for atom in molecule.atoms), \
                proxy.site_type
        assert all(b.source == "construction" for b in proxy.artificial_boundaries)
        assert proxy.repeat_unit is not None or not proxy.artificial_boundaries


def test_proxies_need_a_short_variant_of_at_least_one_unit():
    with pytest.raises(ValueError, match="at least two"):
        _declared(_spec(S.pe_mapping({10: 1})), units=1)


# -- PS parity -------------------------------------------------------------------

# the built-in set's names -> the declared set's (generic) names
RENAMED = {"pristine": "pristine", "interior_radical": "interior_radical",
           "doubly_featured": "doubly_featured", "end_radical": "head_end_radical",
           "benzylic_end_radical": "tail_end_radical", "junction_radical": "pendant_radical",
           "styrene": "monomer"}
# declared pairs the built-in set does not try (every single radical meets
# every single radical and the monomer)
ADDED_PAIRS = {"interior_radical+interior_radical", "interior_radical+head_end_radical",
               "interior_radical+tail_end_radical", "interior_radical+pendant_radical",
               "pendant_radical+pendant_radical", "interior_radical+monomer",
               "pendant_radical+monomer"}


def _renamed(site_type):
    name, _, depth = site_type.partition("@")
    parts = [RENAMED[part] for part in name.split("+")]
    order = ["interior_radical", "head_end_radical", "tail_end_radical", "pendant_radical",
             "pristine", "monomer", "doubly_featured"]
    return "+".join(sorted(parts, key=order.index)) + (f"@{depth}" if depth else "")


def _same_participant(legacy_build, declared_build):
    """Isomorphic molecules whose port atoms correspond tag for tag, with the
    same port kinds and caps."""
    if _ports(legacy_build) != _ports(declared_build):
        return False
    a, b = legacy_build.molecule, declared_build.molecule
    by_id = {atom.id: atom for atom in b.atoms}
    a_by_id = {atom.id: atom for atom in a.atoms}
    pins = {}
    for port in legacy_build.metadata.ports:
        other = next(p for p in declared_build.metadata.ports if p.tag == port.tag)
        pins[a_by_id[port.atom_id]] = by_id[other.atom_id]
    return bool(a.find_isomorphism(b, initial_map=pins, strict=False)) if pins else \
        a.is_isomorphic(b, strict=False)


def _participants_match(legacy, declared):
    remaining = list(declared.boundary_builds)
    for build in legacy.boundary_builds:
        match = next((d for d in remaining if _same_participant(build, d)), None)
        if match is None:
            return False
        remaining.remove(match)
    return not remaining


def test_the_ps_declaration_reproduces_the_builtin_set_structure_by_structure():
    """Every built-in proxy has a declared proxy of the same participants
    (molecules, port atoms, port kinds and caps), named generically; the
    differences are exactly the ones listed here."""
    legacy = {p.site_type: p for p in ps_proxy_set()}
    declared = {p.site_type: p for p in _declared(_spec(ps_deck_mapping({6: 3})))}
    unmatched = {}
    for name, proxy in legacy.items():
        other = declared.get(_renamed(name))
        assert other is not None, (name, _renamed(name))
        if not _participants_match(proxy, other):
            unmatched[name] = _renamed(name)
    # at five units the built-in set pairs the FULL head end radical with the
    # monomer (its short variant only at three units); the declared set
    # always pairs the short variant
    assert unmatched == {"end_radical+styrene@5": "head_end_radical+monomer@5"}
    extra = set(declared) - {_renamed(name) for name in legacy}
    assert extra == ADDED_PAIRS | {f"{name}@5" for name in ADDED_PAIRS}
    # the padding repeat unit is the same chemistry
    unit = next(p.repeat_unit for p in declared.values() if p.repeat_unit is not None)
    assert unit.molecule.is_isomorphic(legacy["pristine"].repeat_unit.molecule)


# -- family candidates -------------------------------------------------------------


class _Families:
    def __init__(self, names):
        self.families = {name: object() for name in names}

    def generate_reactions_from_families(self, *args, **kwargs):
        raise AssertionError("generation bypassed the hook")


def test_a_family_outside_the_builtin_candidates_is_tried_and_recorded_when_silent():
    family = "1,2_shiftC"
    assert family not in PS_FAMILY_CANDIDATES
    calls = []

    def provider(reactants, products=None, only_families=None, resonance=True):
        calls.append(tuple(only_families))
        return []

    proxies = _declared(_spec(S.pe_mapping({10: 1})))
    active, excluded, _cache = EventSetCompiler.discover_family_reactions(
        _Families([family]), proxies, [family], family_universe=[family],
        generation_provider=provider, proxy_set_name="declared proxy set")
    assert calls and set(calls) == {(family,)}
    assert len(calls) == len(proxies)
    assert active == []
    assert "no reaction" in excluded[family] and "PS" not in excluded[family]


# -- the compile CLI ------------------------------------------------------------------


class _Recorder:
    """EventSetCompiler stand-in: records discovery and construction."""
    seen = {}

    @staticmethod
    def discover_family_reactions(kinetics, proxies, candidates, **kwargs):
        _Recorder.seen["proxies"] = proxies
        _Recorder.seen["candidates"] = tuple(candidates)
        _Recorder.seen["discover"] = kwargs
        return [], {}, {}

    def __init__(self, kinetics, proxies, active, **kwargs):
        _Recorder.seen["kwargs"] = kwargs

    def write_artifact(self, directory):
        return Path(directory) / "artifact.json", {}


def _compile(monkeypatch, tmp_path, deck=None, families=("R_Recombination",)):
    database = SimpleNamespace(
        load_kinetics=lambda *a, **k: None, load_thermo=lambda *a, **k: None,
        kinetics=SimpleNamespace(families={}, generate_reactions_from_families=None),
        thermo=object())
    monkeypatch.setattr(compile_module, "RMGDatabase", lambda: database)
    monkeypatch.setattr(compile_module, "prepare_rate_rules", lambda *a, **k: None)
    monkeypatch.setattr(compile_module, "resolve_database_declaration", lambda sha: "sha")
    monkeypatch.setattr(compile_module, "EventSetCompiler", _Recorder)
    _Recorder.seen = {}
    request = compile_module.CompileRequest(
        database_path=str(tmp_path), output_dir=str(tmp_path / "out"), families=families,
        family_universe=families, deck=None if deck is None else str(deck))
    compile_module.compile_artifact(request)
    return _Recorder.seen


def test_compile_with_a_deck_uses_its_declared_proxies_and_names_its_digest(monkeypatch,
                                                                           tmp_path):
    from rmgpy.kmc.run import read_deck

    deck_path = pe_deck(tmp_path)
    seen = _compile(monkeypatch, tmp_path, deck_path, families=("R_Recombination", "1,2_shiftC"))
    deck = read_deck(deck_path)
    names = {p.site_type for p in seen["proxies"]}
    assert "head_end_radical" in names and "pristine" in names
    for proxy in seen["proxies"]:
        for participant in proxy.reactants:
            molecule = _molecule(participant)
            assert not any(molecule.is_atom_in_cycle(atom) for atom in molecule.atoms)
    assert seen["candidates"] == ("R_Recombination", "1,2_shiftC")
    assert seen["kwargs"]["family_candidates"] == ("R_Recombination", "1,2_shiftC")
    assert seen["discover"]["proxy_set_name"] == "declared proxy set"
    assert "PS" not in seen["kwargs"]["family_filter_reason"]
    source = seen["kwargs"]["proxy_source"]
    assert source["kind"] == "declaration"
    assert source["declaration_sha256"] == deck.digests()["declaration_sha256"]
    assert source["chemistry_sha256"] == deck.digests()["chemistry_sha256"]
    assert source["deck_sha256"] == deck.sha256
    assert source["proxy_units"] == 3


def test_compile_without_a_deck_keeps_the_builtin_set_and_says_so(monkeypatch, tmp_path):
    seen = _compile(monkeypatch, tmp_path)
    assert [p.site_type for p in seen["proxies"]] == [p.site_type for p in ps_proxy_set()]
    assert seen["kwargs"]["proxy_source"]["kind"] == "built-in"
    assert "family_candidates" not in seen["kwargs"]
    assert "proxy_set_name" not in seen["discover"]


def test_the_cli_passes_the_deck(monkeypatch, tmp_path):
    families = tmp_path / "families.json"
    families.write_text(json.dumps(["R_Recombination"]))
    deck_path = ps_deck(tmp_path)
    requests = []
    monkeypatch.setattr(compile_module, "compile_in_child",
                        lambda request: requests.append(request) or (tmp_path / "a.json", "a"))
    assert compile_module.main([str(tmp_path), str(tmp_path / "out"), "--families",
                                str(families), "--deck", str(deck_path)]) == 0
    assert requests[0].deck == str(deck_path.resolve())


def test_the_artifact_provenance_names_the_proxy_source(tmp_path):
    from compilerTest import _Database, _Molecule, _Reaction, _Species
    from rmgpy.kmc.compiler import SiteProxy

    molecule = _Molecule()
    source = {"kind": "declaration", "declaration_sha256": "d" * 64}
    compiler = EventSetCompiler(_Database(_Reaction(molecule)),
                                [SiteProxy("interior", [_Species(molecule)])],
                                ["H_Abstraction"], temperature_grid=[600.0, 700.0],
                                proxy_source=source)
    _path, artifact = compiler.write_artifact(tmp_path)
    assert artifact["provenance"]["proxy_source"] == source
