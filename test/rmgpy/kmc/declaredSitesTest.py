"""A declared polymer's records are certified under its own species rules,
and the compile and the roster offer every distinct radical site
(``declared_proxies.radical_sites``, ``emit.declared_reference_systems``)."""

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

import melt_engine_support as S
from rmgpy.kmc import emit
from rmgpy.kmc.run import read_deck
from rmgpy.polymer_spec import LinearHomopolymerSpec, normalized_linear_homopolymer_spec

DECKS = Path(__file__).parent / "fixtures" / "decks"
SCRIPT = Path(__file__).parents[3] / "scripts" / "check-kmc-reference-mapping.py"
_LOADER = importlib.util.spec_from_file_location("check_kmc_reference_mapping_sites", SCRIPT)
CHECKER = importlib.util.module_from_spec(_LOADER)
_LOADER.loader.exec_module(CHECKER)


def _deck(name):
    return read_deck(DECKS / name).spec


def _pe():
    return normalized_linear_homopolymer_spec(LinearHomopolymerSpec.from_mapping(S.pe_mapping({10: 1})))


def _ps():
    from meltRunTest import ps_deck_mapping
    return normalized_linear_homopolymer_spec(LinearHomopolymerSpec.from_mapping(ps_deck_mapping({6: 3})))


def _declared_artifact(spec, units=2):
    from rmgpy.kmc.run import chemistry_sha256, declaration_sha256
    return {"provenance": {"proxy_source": {
        "kind": "declaration", "label": spec.label, "deck_sha256": "deck",
        "chemistry_sha256": chemistry_sha256(spec), "declaration_sha256": declaration_sha256(spec),
        "proxy_units": units, "proxy_set": "declared"}}, "records": []}


def _singles(spec, units=2):
    from rmgpy.kmc.declared_proxies import declared_proxy_set
    return [p for p in declared_proxy_set(spec, units) if len(p.reactants) == 1]


def _radical_atoms(molecule):
    return [atom for atom in molecule.atoms if atom.radical_electrons]


# -- gap A: certification under the declaration's own species rules ------------------


def test_a_declared_artifact_is_certified_under_its_declarations_rules():
    pp, novolac = _deck("pp_deck.py"), _deck("novolac_deck.py")
    rules = emit.certification_rules(_declared_artifact(pp), pp)
    assert rules["domain"].name == "declared_interim_v1"
    assert rules["repeat"] == emit.declared_provider(pp).unit_formula()
    # the depository exclusion grows the witness at the declared chain's sites
    assert rules["growth"].name == emit.declared_provider(pp).name
    assert rules["growth"].unit_formula() == rules["repeat"]
    assert emit.certification_rules(_declared_artifact(novolac), novolac)["domain"].name == \
        "declared_interim_v1[C,H,O]"


def test_a_legacy_artifact_keeps_the_v1_rules():
    assert emit.certification_rules({"records": []}, _deck("pp_deck.py")) == {}
    built_in = {"provenance": {"proxy_source": {"kind": "built-in"}}, "records": []}
    assert emit.certification_rules(built_in, _deck("pp_deck.py")) == {}


def test_the_declared_rules_accept_what_the_ps_rules_refuse():
    from rmgpy.molecule import Molecule
    from rmgpy.kmc import domain as dom

    for name, smiles in (("pp_deck.py", "C[C](C)CC(C)C"), ("novolac_deck.py", "Cc1ccccc1[O]")):
        spec = _deck(name)
        domain = emit.certification_rules(_declared_artifact(spec), spec)["domain"]
        molecule = Molecule(smiles=smiles)
        assert domain.check(molecule).accepted, (name, smiles)
        assert not dom.check(molecule).accepted, (name, smiles)


def test_the_checker_certifies_a_declared_artifact_under_its_rules(tmp_path, monkeypatch):
    from meltRunTest import pe_deck
    from rmgpy.kmc.event_record import EventRecord

    deck = pe_deck(tmp_path)
    spec = read_deck(deck).spec
    artifact = tmp_path / "artifact.json"
    document = _declared_artifact(spec)
    document["records"] = [{"event_id": "e0", "status": "enabled", "site_type": "pristine",
                            "orientation": "as_generated", "family": "family",
                            "template": "t", "reactant_graphs": ["g"], "product_graphs": ["g"],
                            "rate_source": {"kind": "RMG family estimate"}}]
    artifact.write_text(json.dumps(document))
    seen = {}

    def certifier(*_args, **kwargs):
        seen.update(kwargs)
        return object()

    def certify(records, _families, _certifier, required=True):
        return [dict(r, rate_certificate={"sha256": "s"}) for r in records], \
            {"s": {"certified": True, "gaps": []}}, {}
    monkeypatch.setattr(CHECKER.emit, "load_database", lambda *_args:
                        SimpleNamespace(kinetics=SimpleNamespace(families={})))
    monkeypatch.setattr(CHECKER.cm, "Certifier", certifier)
    monkeypatch.setattr(CHECKER.cm, "certify_records", certify)
    monkeypatch.setattr(EventRecord, "from_dict", staticmethod(dict))
    monkeypatch.setattr(CHECKER.emit, "reference_cases", lambda *_args: ({}, {}))
    monkeypatch.setattr(CHECKER.emit, "record_key", lambda *_args: ("k",))
    args = SimpleNamespace(artifact=artifact, deck=deck, database="db", old_rows=None,
                           roster_out=None, dps=(1, 2, 3), limit=None, legacy_artifact=False)
    CHECKER.check(args)
    assert seen["domain"].name == "declared_interim_v1"
    assert seen["repeat"] == emit.declared_provider(spec).unit_formula()


def test_the_emitter_certifies_a_declared_artifact_under_its_rules(tmp_path, monkeypatch):
    spec = _pe()
    artifact = tmp_path / "artifact.json"
    artifact.write_text(json.dumps(_declared_artifact(spec)))
    seen = {}

    class Stop(Exception):
        pass

    def certifier(*_args, **kwargs):
        seen.update(kwargs)
        raise Stop
    monkeypatch.setattr(emit, "load_database", lambda *_args: SimpleNamespace(
        kinetics=SimpleNamespace(families={}), thermo=None))
    from rmgpy.kmc import certificates as cm
    monkeypatch.setattr(cm, "Certifier", certifier)
    emitter = emit.Emitter(artifact, "db", {}, polymer=spec)
    with pytest.raises(Stop):
        emitter.init()
    assert seen["domain"].name == "declared_interim_v1"


# -- gap C: one structure per class of hydrogen-bearing atoms ------------------------


def test_one_radical_site_per_class_of_hydrogen_bearing_atoms():
    from rmgpy.kmc.declared_proxies import radical_sites

    counts = {name: len(radical_sites(emit.declared_provider(spec)))
              for name, spec in (("pe", _pe()), ("ps", _ps()), ("pp", _deck("pp_deck.py")),
                                 ("novolac", _deck("novolac_deck.py")))}
    # PE: CH2; PS: CH2, benzylic CH, ortho, meta, para; PP: CH2, CH, CH3;
    # novolac: the bridge CH2, ring meta (two equivalent), ring para, OH
    assert counts == {"pe": 1, "ps": 5, "pp": 3, "novolac": 4}


def test_pp_gets_its_tertiary_radical():
    tertiary = [p.site_type for p in _singles(_deck("pp_deck.py"))
                for atom in _radical_atoms(p.reactants[0].molecule[0])
                if sum(1 for other in atom.bonds if other.is_carbon()) == 3]
    assert tertiary


def test_novolac_gets_its_phenoxy_radical():
    phenoxy = [p.site_type for p in _singles(_deck("novolac_deck.py"))
               for atom in _radical_atoms(p.reactants[0].molecule[0]) if atom.is_oxygen()]
    assert phenoxy


# -- gap B: the roster is built from the compile structures --------------------------


def _same_structures(proxy_builds, participants, name):
    remaining = [p.build(name).molecule for p in participants]
    for build in proxy_builds:
        match = next((m for m in remaining if m.is_isomorphic(build.molecule)), None)
        if match is None:
            return False
        remaining.remove(match)
    return not remaining


def _proxies_without_a_system(spec, roster, units=2):
    from rmgpy.kmc.declared_proxies import declared_proxy_set

    systems, refusals = emit.build_reference_systems(roster)
    assert refusals == []
    return sorted(p.site_type for p in declared_proxy_set(spec, units)
                  if not any(_same_structures(p.boundary_builds, participants, name)
                             for name, participants in systems.items()))


@pytest.mark.parametrize("which", ("pe", "pp"))
def test_every_compile_structure_has_a_roster_system(which):
    """Singles and pairs, doubly-featured and pendant radicals included."""
    spec = _pe() if which == "pe" else _deck("pp_deck.py")
    roster = emit.declared_reference_systems(spec, 2)
    assert _proxies_without_a_system(spec, roster) == []


def test_the_checker_builds_a_declared_artifacts_roster_from_its_compile_structures():
    spec = _pe()
    roster = CHECKER.reference_roster(_declared_artifact(spec, units=2), spec, (1, 2, 3))
    assert _proxies_without_a_system(spec, roster) == []
    assert any(name.startswith("single:doubly_featured") for name in roster.systems)


def test_the_checker_keeps_the_generated_roster_for_a_legacy_artifact():
    spec = _pe()
    roster = CHECKER.reference_roster({"records": [{"site_type": "pristine@4"}]}, spec, (1, 2))
    assert "single:chain:dp4:capped" in roster.systems


def test_a_v3_roster_round_trips(tmp_path):
    spec = _deck("pp_deck.py")
    roster = emit.declared_reference_systems(spec, 2)
    written = emit.write_roster(roster, tmp_path / "v3.json")
    assert json.loads((tmp_path / "v3.json").read_text())["schema"] == emit.ROSTER_SCHEMA_V3
    assert written.digests() == roster.digests()
    built, refusals = emit.build_reference_systems(written)
    assert refusals == [] and set(built) == set(roster.systems)
    assert _proxies_without_a_system(spec, written) == []


def test_a_v2_roster_still_loads_and_an_unknown_schema_is_refused_by_name(tmp_path):
    v2 = emit.write_roster(emit.reference_systems(_pe(), None, dps=(1, 2)), tmp_path / "v2.json")
    assert emit.build_reference_systems(v2)[1] == []
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"schema": "kmc-reference-roster/9", "polymer": {}, "systems": {}}))
    with pytest.raises(emit.RosterRefusal, match="refused:roster"):
        emit.load_roster(bad)
