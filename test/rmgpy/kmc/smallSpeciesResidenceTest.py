"""Where each small species lives, declared once and honoured by both methods.

``small_species.<name>.residence`` (optional) is ``"melt"``, the default and
the canonical form (absent from the normalized declaration, so stating it
moves no digest); ``"volatile"`` is a later extension and refused.  H and H2
may be declared (disposition ``mobile``, their built-in graphs): the engine
treats them as it treats the built-ins.  The MoM producer forwards every
species' residence into ``polymer_phase(..., residence=...)``: the declared
species, H and H2 and every other reactive species (``'*'``) melt, the bath
gas."""

import ast
import json

import numpy
import pytest

import melt_engine_support as S
from meltRunTest import pe_catalog_file
from momProducerTest import make_setup, produce
from polymer_spec_v2_support import v2_mapping

from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import melt_observables as mo
from rmgpy.kmc import mom_producer as mp
from rmgpy.kmc import run as kr
from rmgpy.kmc.melt_state import MeltState, Node, Origin
from rmgpy.polymer_spec import LinearHomopolymerSpec, PolymerSpec

H_GRAPH = {"atoms": [S._atom("h", "H", radical=1)], "bonds": [], "electronic_state": {"multiplicity": 2}}
H2_GRAPH = {"atoms": [S._atom("h1", "H"), S._atom("h2", "H")], "bonds": [{"atoms": ["h1", "h2"], "order": 1}],
            "electronic_state": {"multiplicity": 1}}


def methyl(**extra):
    return dict(S.METHYL, **extra)


def spec(small):
    return LinearHomopolymerSpec.from_mapping(S.pe_mapping({10: 1, 20: 1}, small=small))


# --------------------------------------------------------------------------
# the declaration
# --------------------------------------------------------------------------
def test_melt_residence_is_the_canonical_default():
    bare, stated = spec({"methyl": methyl()}), spec({"methyl": methyl(residence="melt")})
    assert kr.declaration_sha256(bare) == kr.declaration_sha256(stated)
    assert kr.chemistry_sha256(bare) == kr.chemistry_sha256(stated)
    assert b"residence" not in stated.normalized_json
    assert stated.residence("methyl") == "melt" == bare.residence("methyl")


def test_volatile_and_unknown_residences_are_refused():
    with pytest.raises(ValueError, match="volatile.*later extension"):
        spec({"methyl": methyl(residence="volatile")})
    with pytest.raises(ValueError, match="residence"):
        spec({"methyl": methyl(residence="gas")})
    with pytest.raises(ValueError, match="graph and disposition"):
        spec({"methyl": methyl(colour="blue")})


def test_v2_declarations_take_the_same_residence_rule():
    mapping = v2_mapping()
    mapping["small_species"] = {"methyl": methyl(residence="melt")}
    stated = PolymerSpec.from_mapping(mapping)
    mapping["small_species"] = {"methyl": methyl()}
    bare = PolymerSpec.from_mapping(mapping)
    assert stated.normalized_json == bare.normalized_json and stated.residence("methyl") == "melt"
    mapping["small_species"] = {"methyl": methyl(residence="volatile")}
    with pytest.raises(ValueError, match="volatile"):
        PolymerSpec.from_mapping(mapping)


def test_h_and_h2_are_declarable_as_the_builtins(tmp_path):
    declared = spec({"H": {"graph": H_GRAPH, "disposition": "mobile", "residence": "melt"},
                     "H2": {"graph": H2_GRAPH, "disposition": "mobile"}})
    bare = spec({})
    # a different declaration: the chemistry digest moves (small_species is chemistry)
    assert kr.chemistry_sha256(declared) != kr.chemistry_sha256(bare)
    assert declared.residence("H") == declared.residence("H2") == "melt" == bare.residence("H")
    identity = me.SpeciesIdentity(declared)
    names = mo.declared_species_names(declared, identity)
    assert sorted(names) == ["H", "H2"] and len(names) == len(set(names))
    h = MeltState({0: Node("H", None, 0, 1, 0, {})}, {0: Origin("mobile", 0, "H", "H", False)})
    assert identity.identify(h, h.components.ids()[0]) == "H"
    # the engine loads it on a catalog in which H takes no part (as a built-in may)
    catalog, _records = kr.load_catalogs([pe_catalog_file(tmp_path)], 700.0)
    engine = me.MeltEngine(declared, catalog, 700.0, numpy.random.default_rng(0))
    assert engine.declared_mobile_participation == {}


def test_h_declared_otherwise_than_the_builtin_is_refused():
    with pytest.raises(me.DeclarationError, match="built-in"):
        me.SpeciesIdentity(spec({"H": {"graph": H2_GRAPH, "disposition": "mobile"}}))
    with pytest.raises(me.DeclarationError, match="built-in"):
        me.SpeciesIdentity(spec({"H2": {"graph": H2_GRAPH, "disposition": "inert"}}))
    with pytest.raises(me.DeclarationError, match="monomer_product"):
        me.SpeciesIdentity(spec({"monomer_product": methyl()}))


# --------------------------------------------------------------------------
# the producer
# --------------------------------------------------------------------------
def _phase_call(text):
    tree = ast.parse(text)
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "polymer_phase":
            return {k.arg: ast.literal_eval(k.value) for k in node.keywords}
    raise AssertionError("no polymer_phase call")


def test_the_producer_forwards_every_residence(tmp_path):
    s = make_setup(tmp_path, small={"methyl": methyl()}, mobile=None)
    out = tmp_path / "out"
    produce(s, out)
    phase = _phase_call((out / "rmg" / "input.py").read_text())
    assert phase["residence"] == {"*": "melt", "H": "melt", "H2": "melt", "N2": "gas", "methyl": "melt"}
    conv = json.loads((out / "conventions.json").read_text())
    names = {c["name"]: c["value"] for c in conv["conventions"]}
    assert names["residence"] == phase["residence"]
    facts = json.loads((out / "producer.json").read_text())["facts"]
    assert facts["residence"] == phase["residence"]


def test_phase_kinetics_compare_melt_residents_with_the_melt_volume(tmp_path):
    """MoM gives a melt-resident species n / V_poly (its polymer volume: the
    pool mass a mu1 over the density); kMC n / V.  The ratio is V / V_poly."""
    s = make_setup(tmp_path)
    out = tmp_path / "out"
    produce(s, out)
    facts = json.loads((out / "producer.json").read_text())["facts"]
    kinetics = facts["phase_kinetics"]
    a, mu1 = facts["a_g_mol"], facts["moments_mol"][1]
    v_poly = a * mu1 / 1000.0 / 850.0
    v_kmc = 1.0 / 850.0
    assert set(kinetics["ratio_by_species"]) == {"H", "H2"}
    for value in kinetics["ratio_by_species"].values():
        assert value == pytest.approx(v_kmc / v_poly, rel=1e-12)
    assert kinetics["V_poly_m3"] == pytest.approx(v_poly, rel=1e-12)
