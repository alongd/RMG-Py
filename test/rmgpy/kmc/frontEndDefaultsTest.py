"""No polystyrene default remains in the kMC front end: an input that names no
polymer is refused, and the PS / linear-alkane providers are reached only by
name (``CompileRequest.legacy_provider``, the explicit owner lookup of
``windows.build_spec``)."""

import ast
from pathlib import Path

import pytest

from rmgpy.kmc import boundaries as b
from rmgpy.kmc import compile as compile_module
from rmgpy.kmc import windows as win
from rmgpy.kmc.compile import CompileRequest
from rmgpy.kmc.compiler import EventSetCompiler, ps_proxy_set, ps_repeat_unit

REPOSITORY = Path(__file__).resolve().parents[3]
DECKS = REPOSITORY / "test/rmgpy/kmc/fixtures/decks"


def request(**options):
    return CompileRequest("/db", "/out", ("R_Recombination",),
                          family_universe=("R_Recombination",), **options)


# -- the compile front end -------------------------------------------------------


def test_a_request_naming_no_polymer_is_refused_with_both_ways_out():
    with pytest.raises(ValueError, match="refused:no_polymer") as error:
        compile_module._proxy_inputs(request())
    assert "deck" in str(error.value)
    assert "legacy_provider='polystyrene'" in str(error.value)


def test_the_refusal_comes_before_the_database_loads(monkeypatch, tmp_path):
    def loads(*args, **kwargs):
        raise AssertionError("the database was loaded before the request was refused")

    monkeypatch.setattr(compile_module, "RMGDatabase", loads)
    with pytest.raises(ValueError, match="refused:no_polymer"):
        compile_module.compile_artifact(
            CompileRequest(str(tmp_path), str(tmp_path / "out"), ("R_Recombination",),
                           family_universe=("R_Recombination",)))


def test_the_named_reference_provider_gives_the_built_in_set_unchanged():
    proxies, source, options, discovery = compile_module._proxy_inputs(
        request(legacy_provider="polystyrene", proxy_units=3))
    assert [p.site_type for p in proxies] == [p.site_type for p in ps_proxy_set(3)]
    assert source == {"kind": "built-in",
                      "proxy_set": "rmgpy.kmc.compiler.ps_proxy_set",
                      "proxy_units": 3, "reason": "no deck was given"}
    assert options == {} and discovery == {}


@pytest.mark.parametrize("value", ["", "polystyrene ", "Polystyrene", "linear_alkane", "pp", 3,
                                   True, ("polystyrene",)])
def test_any_other_legacy_provider_is_refused_never_defaulted(value):
    with pytest.raises(ValueError, match="refused:no_polymer"):
        compile_module._proxy_inputs(request(legacy_provider=value))


def test_a_deck_and_a_legacy_provider_together_are_refused():
    with pytest.raises(ValueError, match="refused:no_polymer"):
        compile_module._proxy_inputs(
            request(deck=str(DECKS / "pp_deck.py"), legacy_provider="polystyrene"))


def test_the_legacy_provider_is_not_part_of_the_request_identity(tmp_path, monkeypatch):
    from generationIdentityTest import compile_request, make_database
    from rmgpy.kmc import generation_identity as gi
    monkeypatch.setenv("PYTHONHASHSEED", "0")
    database = make_database(tmp_path / "db")
    named = compile_request(database, legacy_provider="polystyrene")
    assert gi.compilation_request_identity(named) == \
        gi.compilation_request_identity(compile_request(database))


def test_the_cli_names_the_legacy_provider(tmp_path):
    import argparse
    families = tmp_path / "families.json"
    families.write_text('["F"]')
    parser = argparse.ArgumentParser()
    parser.add_argument("database")
    compile_module.add_request_arguments(parser)
    base = ["/db", "--families", str(families)]
    named = compile_module.request_from_arguments(
        parser, parser.parse_args(base + ["--legacy-provider", "polystyrene"]), "/out")
    assert named.legacy_provider == "polystyrene" and named.deck is None
    bare = compile_module.request_from_arguments(parser, parser.parse_args(base), "/out")
    assert bare.legacy_provider is None


def test_a_work_directory_begun_before_the_field_matches_the_same_request():
    from rmgpy.kmc.compile_shards import _request_identity
    named = _request_identity(request(legacy_provider="polystyrene"))
    assert "legacy_provider" not in named
    assert named == _request_identity(request())


# -- the witness builder -----------------------------------------------------------


@pytest.mark.parametrize("spec", [b.PsChainSpec(3), ("alkane", 6), ("small", "H2")])
def test_build_spec_builds_a_spec_its_named_owner_claims(spec):
    assert win.build_spec(spec).molecule is not None


@pytest.mark.parametrize("spec", [None, 3, "polystyrene", ("styrene", 3), ("alkane",),
                                  {"units": 3}, object()])
def test_build_spec_refuses_a_bare_spec_that_names_no_provider(spec):
    with pytest.raises(win.WindowError, match="refused:legacy_spec"):
        win.build_spec(spec)


def test_build_spec_has_no_polystyrene_fallback():
    source = (REPOSITORY / "rmgpy/kmc/windows.py").read_text()
    tree = ast.parse(source)
    function = next(n for n in ast.walk(tree)
                    if isinstance(n, ast.FunctionDef) and n.name == "build_spec")
    names = {c.value for c in ast.walk(function)
             if isinstance(c, ast.Constant) and isinstance(c.value, str)}
    assert "polystyrene" not in names     # the alkane branch is its own named owner


def test_the_spec_owners_derive_from_the_one_reference_list():
    assert win.SPEC_OWNERS == b.LEGACY_PROVIDERS + (win.SMALL_PROVIDER,)


def test_provider_defaults_are_explicit_arguments():
    import inspect
    for function in (win.SearchClassifier.__init__, win._fragment_facts):
        parameter = inspect.signature(function).parameters["providers"]
        assert parameter.default is inspect.Parameter.empty, function.__qualname__
    facts = inspect.signature(win.certificate_support).parameters["fragment_facts"]
    assert facts.default is inspect.Parameter.empty
    reference = inspect.signature(win.observation_providers).parameters["reference"]
    assert reference.default is inspect.Parameter.empty and reference.kind is reference.KEYWORD_ONLY


def test_pair_roster_is_not_a_front_end_function():
    assert not hasattr(win, "pair_roster")


def test_the_identitys_default_proxy_units_is_the_declared_default():
    # generation_identity hashes the PS constant for an unset proxy_units; a
    # declared compile defaults to DEFAULT_PROXY_UNITS: the same number
    from rmgpy.kmc import compiler
    from rmgpy.kmc.declared_proxies import DEFAULT_PROXY_UNITS
    assert compiler.PS_PROXY_UNITS == DEFAULT_PROXY_UNITS


# -- padding provenance --------------------------------------------------------------


def padding_unit(proxies):
    compiler = EventSetCompiler(None, proxies, [], proxy_padding_distance=3)
    return compiler._compile_provenance()[1]["proxy_boundary_padding"]["repeat_unit"]


def test_padding_provenance_records_the_declared_polymers_repeat_unit_not_ps():
    from rmgpy.kmc.declared_proxies import declared_proxy_set
    from rmgpy.kmc.run import read_deck
    proxies = declared_proxy_set(read_deck(DECKS / "pp_deck.py").spec, 3)
    unit = padding_unit(proxies)
    assert unit != ps_repeat_unit().provenance()
    assert unit["name"] != "polystyrene"
    assert "c1ccccc1" not in unit["adjacency_list"] and "{7,B}" not in unit["adjacency_list"]


def test_padding_provenance_of_the_ps_set_is_the_ps_unit():
    proxies = tuple(p for p in ps_proxy_set(3) if p.repeat_unit is not None)
    assert proxies
    assert padding_unit(proxies) == ps_repeat_unit().provenance()


def test_padding_provenance_with_no_polymer_names_none():
    assert padding_unit(()) is None


# -- every caller that wants PS asks for it ---------------------------------------------


def deckless_request_calls(path):
    tree = ast.parse(path.read_text())
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            name = getattr(node.func, "attr", getattr(node.func, "id", None))
            if name == "CompileRequest":
                keywords = {k.arg for k in node.keywords}
                if None in keywords:        # **options: cannot be read statically
                    continue
                if not keywords & {"deck", "legacy_provider"}:
                    yield node.lineno


def test_every_compile_that_wants_ps_names_it():
    # Scripts that only hash a request (generationIdentityTest) never compile it.
    hash_only = {"generationIdentityTest.py", "frontEndDefaultsTest.py"}
    offenders = [f"{p.name}:{line}"
                 for p in sorted((REPOSITORY / "test/rmgpy/kmc").glob("*.py"))
                 if p.name not in hash_only
                 for line in deckless_request_calls(p)]
    assert offenders == []
