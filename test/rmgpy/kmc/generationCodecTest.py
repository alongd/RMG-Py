"""Codec for generated reactions: lossless on real and populated fixtures, bounded on hostile bytes."""

import copy
import json
import os
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.kinetics.arrhenius import Arrhenius, ArrheniusBM, ArrheniusEP
from rmgpy.kinetics.uncertainties import RateUncertainty
from rmgpy.kmc import compiler, generation_codec
from rmgpy.kmc.generation_codec import (
    DEFAULT_LIMITS,
    DecodeError,
    DecodeLimits,
    EncodeError,
    decode_generated_reactions,
    encode_generated_reactions,
    roundtrip_generation_provider_factory,
)
from rmgpy.molecule.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.thermo.thermodata import ThermoData

DATABASE = Path(os.environ.get("RMG_DATABASE_PATH", "/home/alon/runs/phase2b/database"))
REAL_FAMILIES = {
    "H_Abstraction": (["C[CH]c1ccccc1", "CC(C)c1ccccc1"], 15),
    "R_Recombination": (["C[CH]c1ccccc1", "C[CH]c1ccccc1"], 6),
    "intra_H_migration": (["CC(C)C[CH]c1ccccc1"], 20),
    # The pinned database has no plain beta-scission family; these two are the
    # scission families that fire on small radical inputs.
    "1,4_Linear_birad_scission": (["[CH2]CC[CH2]"], 1),
    "Cyclopentadiene_scission": (["C1=CC=CC1"], 1),
}


# ---------------------------------------------------------------------------
# Structural comparison that also checks sharing
# ---------------------------------------------------------------------------


class Pairing:
    """Bijection between original and decoded instances; sharing must match."""

    def __init__(self):
        self.forward = {}
        self.backward = {}

    def check(self, original, decoded):
        """True if this pair was already compared (so recursion can stop)."""
        a, b = id(original), id(decoded)
        if a in self.forward or b in self.backward:
            assert self.forward.get(a) == b and self.backward.get(b) == a, "sharing differs"
            return True
        self.forward[a], self.backward[b] = b, a
        return False


def same_quantity(a, b):
    if a is None or b is None:
        assert a is b
        return
    assert type(a) is type(b)
    assert a.units == b.units and a.uncertainty_type == b.uncertainty_type
    assert np.array_equal(a.value_si, b.value_si)
    assert np.array_equal(a.uncertainty_si, b.uncertainty_si)


def same_atom(pairing, a, b):
    if pairing.check(a, b):
        return
    for name in ("radical_electrons", "charge", "label", "lone_pairs", "site", "morphology", "id",
                 "connectivity1", "connectivity2", "connectivity3", "sorting_label", "terminal", "ignore"):
        assert getattr(a, name) == getattr(b, name), name
    assert a.element is b.element
    assert a.atomtype is b.atomtype
    assert a.props == b.props
    assert len(a.edges) == len(b.edges)
    for (na, ba), (nb, bb) in zip(a.edges.items(), b.edges.items()):  # edge order matters
        same_atom(pairing, na, nb)
        if not pairing.check(ba, bb):
            assert ba.order == bb.order and type(ba.order) is type(bb.order)
            same_atom(pairing, ba.vertex1, bb.vertex1)
            same_atom(pairing, ba.vertex2, bb.vertex2)


def same_molecule(pairing, a, b):
    if pairing.check(a, b):
        return
    assert type(a) is type(b)
    for name in ("multiplicity", "symmetry_number", "reactive", "is_polymer_proxy", "props", "metal", "facet"):
        assert getattr(a, name) == getattr(b, name), name
    assert len(a.atoms) == len(b.atoms)
    for x, y in zip(a.atoms, b.atoms):  # atom order matters
        same_atom(pairing, x, y)


def same_thermo(pairing, a, b):
    if a is None or b is None:
        assert a is b
        return
    if pairing.check(a, b):
        return
    assert type(a) is type(b)
    for name in ("_Tdata", "_Cpdata", "_H298", "_S298", "_Cp0", "_CpInf", "_Tmin", "_Tmax", "_E0"):
        same_quantity(getattr(a, name), getattr(b, name))
    assert a.label == b.label and a.comment == b.comment


def same_species(pairing, a, b):
    if pairing.check(a, b):
        return
    assert type(a) is type(b)
    for name in ("index", "label", "reactive", "props", "aug_inchi", "symmetry_number", "is_solvent",
                 "is_moment_dummy", "_is_polymer_proxy", "creation_iteration", "explicitly_allowed"):
        assert getattr(a, name) == getattr(b, name), name
    same_quantity(a._molecular_weight, b._molecular_weight)
    same_thermo(pairing, a.thermo, b.thermo)
    assert len(a.molecule) == len(b.molecule)
    for x, y in zip(a.molecule, b.molecule):  # resonance order matters
        same_molecule(pairing, x, y)


def same_kinetics(pairing, a, b):
    if a is None or b is None:
        assert a is b
        return
    if pairing.check(a, b):
        return
    assert type(a) is type(b)
    assert a.comment == b.comment
    for name in ("_A", "_n", "_Ea", "_T0", "_alpha", "_w0", "_E0", "_Tmin", "_Tmax", "_Pmin", "_Pmax"):
        if hasattr(a, name):
            same_quantity(getattr(a, name), getattr(b, name))
    ua, ub = a.uncertainty, b.uncertainty
    if ua is None or ub is None:
        assert ua is ub
    else:
        for name in ("mu", "var", "Tref", "N", "data_mean", "correlation"):
            assert getattr(ua, name) == getattr(ub, name), name
    assert a.solute is None and b.solute is None


def same_reaction(pairing, a, b):
    if pairing.check(a, b):
        return
    assert type(a) is type(b)
    for name in ("index", "label", "reversible", "duplicate", "degeneracy", "electrons", "_protons", "allow_pdep_route",
                 "elementary_high_p", "comment", "is_forward", "allow_max_rate_violation", "is_end_group_reaction", "rank"):
        assert getattr(a, name) == getattr(b, name), name
    assert getattr(a, "__dict__", {}).keys() == getattr(b, "__dict__", {}).keys()
    for side in ("reactants", "products"):
        assert len(getattr(a, side)) == len(getattr(b, side))
        for x, y in zip(getattr(a, side), getattr(b, side)):
            same_species(pairing, x, y)
    same_kinetics(pairing, a.kinetics, b.kinetics)
    if a.pairs is None or b.pairs is None:
        assert a.pairs is b.pairs
    else:
        assert len(a.pairs) == len(b.pairs)
        for pa, pb in zip(a.pairs, b.pairs):
            assert type(pa) is type(pb) and len(pa) == len(pb)
            for x, y in zip(pa, pb):
                same_species(pairing, x, y)
    if isinstance(a, TemplateReaction):
        assert a.family == b.family and a.template == b.template and a.estimator == b.estimator
        if a.reverse is None or b.reverse is None:
            assert a.reverse is b.reverse
        else:
            same_reaction(pairing, a.reverse, b.reverse)


def assert_same(original, decoded):
    assert len(original) == len(decoded)
    pairing = Pairing()
    for a, b in zip(original, decoded):
        same_reaction(pairing, a, b)


def roundtrip(reactions, limits=DEFAULT_LIMITS):
    payload = encode_generated_reactions(reactions)
    assert payload == encode_generated_reactions(reactions), "encoding must be deterministic"
    return payload, decode_generated_reactions(payload, limits=limits)


# ---------------------------------------------------------------------------
# Compiler-consumer equality
# ---------------------------------------------------------------------------


def consumer_view(reaction):
    token = compiler._PAIR_MOLECULE_CACHE.set({})
    try:
        participants = reaction.reactants + reaction.products
        proxy_same = SimpleNamespace(reactants=reaction.reactants)
        proxy_reverse = SimpleNamespace(reactants=reaction.products)
        _, forward_orientation = compiler._orient_to_proxy(proxy_same, reaction)
        reversed_reaction, reverse_orientation = compiler._orient_to_proxy(proxy_reverse, reaction)
        return {
            "atom_map": compiler._canonical_atom_map(reaction),
            "bond_ops": compiler._bond_ops(reaction),
            "molecules": [
                (
                    compiler._molecule(item).to_adjacency_list(),
                    [atom.id for atom in compiler._molecule(item).atoms],
                )
                for item in participants
            ],
            "orientation": (forward_orientation, reverse_orientation),
            "reversed_sides": (len(reversed_reaction.reactants), len(reversed_reaction.products)),
            "heavy": compiler._formula(participants),
        }
    finally:
        compiler._PAIR_MOLECULE_CACHE.reset(token)


def assert_consumers_agree(original, decoded):
    for a, b in zip(original, decoded):
        assert consumer_view(a) == consumer_view(b)


# ---------------------------------------------------------------------------
# Real generator fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def real_generation():
    if not (DATABASE / "input/kinetics/families").is_dir():
        pytest.skip("pinned RMG database is not available")
    from rmgpy.data.rmg import RMGDatabase

    database = RMGDatabase()
    database.load_kinetics(
        str(DATABASE / "input/kinetics"),
        reaction_libraries=[],
        seed_mechanisms=None,
        kinetics_families=list(REAL_FAMILIES),
        kinetics_depositories=["training"],
    )
    database.load_thermo(
        str(DATABASE / "input/thermo"), thermo_libraries=["primaryThermoLibrary"], depository=True
    )
    compiler.prepare_rate_rules(database.kinetics, database.thermo, verbose=False)
    results = {}
    for family, (smiles, expected) in REAL_FAMILIES.items():
        results[family] = database.kinetics.generate_reactions_from_families(
            [Species(molecule=[Molecule(smiles=item)]) for item in smiles], only_families=[family]
        )
        assert len(results[family]) == expected, family
    return SimpleNamespace(database=database, results=results)


@pytest.mark.parametrize("family", list(REAL_FAMILIES))
def test_real_generator_output_round_trips_exactly(real_generation, family):
    reactions = real_generation.results[family]
    assert reactions and all(isinstance(item, TemplateReaction) for item in reactions)
    _, decoded = roundtrip(reactions)
    assert_same(reactions, decoded)
    assert_consumers_agree(reactions, decoded)


def test_real_output_exercises_the_features_the_codec_claims(real_generation):
    every = [item for reactions in real_generation.results.values() for item in reactions]
    species = [s for r in every for s in r.reactants + r.products]
    molecules = [m for s in species for m in s.molecule]
    atoms = [a for m in molecules for a in m.atoms]
    assert any(len(s.molecule) > 1 for s in species), "multi-resonance species"
    assert any(a.label for a in atoms), "kept atom labels"
    assert any(a.props for a in atoms), "atom props"
    assert any(r.reverse is not None for r in every), "reverse links"
    assert any(r.pairs for r in every), "pairs"
    assert len({id(m) for m in molecules}) < len(molecules), "shared molecule instances"
    assert any(r.degeneracy != 1 for r in every), "non-unit degeneracy"
    assert any(r.template for r in every), "templates"
    assert any(not r.is_forward for r in every) or all(r.is_forward for r in every)


def test_decoded_real_output_keeps_distinct_isomorphic_instances_distinct(real_generation):
    reactions = real_generation.results["R_Recombination"]
    _, decoded = roundtrip(reactions)
    species = [s for r in decoded for s in r.reactants + r.products]
    isomorphic = [
        (x, y)
        for i, x in enumerate(species)
        for y in species[i + 1:]
        if x is not y and x.molecule[0].is_isomorphic(y.molecule[0])
    ]
    assert isomorphic, "fixture must contain distinct isomorphic species"
    original = [s for r in reactions for s in r.reactants + r.products]
    assert len({id(s) for s in species}) == len({id(s) for s in original})


# ---------------------------------------------------------------------------
# Deliberately populated fixtures
# ---------------------------------------------------------------------------


def species_of(smiles, resonance=False, label=""):
    species = Species(molecule=[Molecule(smiles=smiles)], label=label)
    if resonance:
        species.generate_resonance_structures()
    for number, atom in enumerate(species.molecule[0].atoms, 1):
        atom.id = number
    for molecule in species.molecule[1:]:
        for atom, source in zip(molecule.atoms, species.molecule[0].atoms):
            atom.id = source.id
    return species


def make_kinetics():
    uncertainty = RateUncertainty(mu=0.25, var=0.5, Tref=1000.0, N=7, correlation="rule-7")
    uncertainty.data_mean = 0.125
    return [
        Arrhenius(
            A=(1.5e6, "cm^3/(mol*s)", "*|/", 3.0), n=(0.5, "", "+|-", 0.1), Ea=(10.2, "kcal/mol", "+|-", 1.0),
            T0=(300, "K"), Tmin=(300, "K"), Tmax=(2000, "K"), Pmin=(1e4, "Pa"), Pmax=(1e7, "Pa"),
            uncertainty=uncertainty, comment="Arrhenius with every field",
        ),
        ArrheniusEP(A=(1e9, "s^-1"), n=0.3, alpha=0.4, E0=(55, "kJ/mol", "+|-", 2.0), Tmin=(500, "K"),
                    uncertainty=copy.copy(uncertainty), comment="EP rate rule"),
        ArrheniusBM(A=(2e3, "m^3/(mol*s)"), n=0.1, w0=(8e5, "J/mol"), E0=(61, "kJ/mol"),
                    Tmax=(1500, "K"), comment="Blowers-Masel"),
    ]


def make_thermo():
    return ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
        Cpdata=([60.2599, 68.0494, 74.7775, 80.9311, 90.8846, 97.4337, 105.393], "J/(mol*K)", "+|-", 1.5),
        H298=(-47.1, "kJ/mol", "+|-", 2.0), S298=(269.5, "J/(mol*K)"), Cp0=(33.2, "J/(mol*K)"),
        CpInf=(103.9, "J/(mol*K)"), Tmin=(300, "K"), Tmax=(2000, "K"), E0=(-40.0, "kJ/mol"),
        label="fixture", comment="hand populated thermo",
    )


def populated_reactions():
    benzyl = species_of("C[CH]c1ccccc1", resonance=True, label="benzyl")
    ethylbenzene = species_of("CCc1ccccc1", label="EB")
    cumene = species_of("CC(C)c1ccccc1", label="cumene")
    benzyl.molecule[0].atoms[0].label = "*1"
    benzyl.thermo = make_thermo()
    benzyl.props["kmc_frozen_thermo"] = True
    kinetics = make_kinetics()
    forward = TemplateReaction(
        reactants=[benzyl, ethylbenzene], products=[cumene, species_of("C[CH]c1ccccc1")],
        kinetics=kinetics[0], degeneracy=2.0, family="H_Abstraction", template=["a", "b"], estimator="rate rules",
        is_forward=True, reversible=True, duplicate=True, electrons=0,
    )
    backward = TemplateReaction(
        reactants=[cumene, benzyl], products=[ethylbenzene, benzyl],
        kinetics=kinetics[1], degeneracy=1.0, family="H_Abstraction", template=["c"], estimator="rate rules",
        is_forward=False,
    )
    forward.reverse, backward.reverse = backward, forward  # a cycle
    forward.pairs = [[benzyl, cumene], (ethylbenzene, forward.products[1])]
    backward.pairs = [(cumene, ethylbenzene)]
    other = TemplateReaction(
        reactants=[benzyl], products=[ethylbenzene], kinetics=kinetics[2], family="intra_H_migration",
        template=["z"], is_forward=True,
    )
    plain = Reaction(reactants=[benzyl], products=[cumene], kinetics=None, label="plain", comment="not templated",
                     index=4, rank=3)
    return [forward, backward, other, plain]


def test_populated_fixture_round_trips_and_covers_every_claimed_feature():
    reactions = populated_reactions()
    forward, backward, other, plain = reactions
    benzyl = forward.reactants[0]
    # features are really present before the round trip
    assert len(benzyl.molecule) > 1
    assert benzyl.molecule[0].atoms[0].label == "*1"
    assert forward.reverse is backward and backward.reverse is forward
    assert isinstance(forward.pairs[0], list) and isinstance(forward.pairs[1], tuple)
    assert benzyl is backward.reactants[1] is other.reactants[0] is plain.reactants[0]
    assert benzyl.thermo is not None and benzyl.props["kmc_frozen_thermo"] is True
    assert {type(r.kinetics).__name__ for r in reactions[:3]} == {"Arrhenius", "ArrheniusEP", "ArrheniusBM"}
    assert forward.kinetics.uncertainty.correlation == "rule-7" and forward.kinetics.comment
    assert forward.kinetics._A.uncertainty_type == "*|/"
    assert not isinstance(plain, TemplateReaction)
    _, decoded = roundtrip(reactions)
    assert_same(reactions, decoded)
    assert decoded[0].reactants[0] is decoded[1].reactants[1] is decoded[2].reactants[0]
    assert decoded[0].reverse is decoded[1] and decoded[1].reverse is decoded[0]


def test_kinetics_and_thermo_values_are_bit_exact():
    reactions = populated_reactions()
    _, decoded = roundtrip(reactions)
    for original, copy_ in zip(reactions, decoded):
        if original.kinetics is not None:
            assert original.kinetics._A.value_si == copy_.kinetics._A.value_si
            assert original.kinetics.comment == copy_.kinetics.comment
            assert (original.kinetics.uncertainty is None) == (copy_.kinetics.uncertainty is None)
            assert original.get_rate_coefficient(800.0) == copy_.get_rate_coefficient(800.0) if hasattr(original, "get_rate_coefficient") and not isinstance(original.kinetics, (ArrheniusEP, ArrheniusBM)) else True
    thermo = reactions[0].reactants[0].thermo
    again = decoded[0].reactants[0].thermo
    assert thermo.get_enthalpy(800.0) == again.get_enthalpy(800.0)
    assert thermo.get_entropy(800.0) == again.get_entropy(800.0)


def test_mapped_reaction_compiles_to_identical_records_after_round_trip():
    from benzylicEndTest import SmallReferenceThermo, mapped_addition

    grid = [600.0, 700.0, 800.0]

    def compile_pair(reaction):
        instance = compiler.EventSetCompiler(
            None, [], [], temperature_grid=grid, reference_thermo_provider=SmallReferenceThermo()
        )
        proxy = compiler.SiteProxy(
            "benzylic_end_radical+styrene", reaction.reactants,
            participant_site_types=("benzylic_end_radical", "styrene"),
        )
        records = instance._linked_family_pair(proxy, reaction, {})[0]
        return [record.to_dict() for record in records]

    reaction, _ = mapped_addition(2)
    reaction.kinetics = make_kinetics()[0]
    reaction.kinetics.comment = "Multiplied by reaction path degeneracy 1.0"
    reaction.template = ["Cds-HH_Cds-CbH", "CsJ-CbCsH"]
    assert reaction.kinetics.uncertainty is not None
    _, decoded = roundtrip([reaction])
    assert_same([reaction], decoded)
    assert_consumers_agree([reaction], decoded)
    original_records = compile_pair(copy.deepcopy(reaction))
    decoded_records = compile_pair(decoded[0])
    assert len(original_records) == 2
    assert original_records == decoded_records
    assert original_records[0]["thermo_provenance"] and original_records[0]["k_table"]


def test_shared_versus_isomorphic_species_are_not_merged_on_decode():
    first, second = species_of("C[CH]c1ccccc1"), species_of("C[CH]c1ccccc1")
    assert first is not second and first.molecule[0].is_isomorphic(second.molecule[0])
    reactions = [
        TemplateReaction(reactants=[first], products=[second], family="x", template=["t"], is_forward=True),
        TemplateReaction(reactants=[first, second], products=[first], family="x", template=["t"], is_forward=True),
    ]
    _, decoded = roundtrip(reactions)
    assert decoded[0].reactants[0] is decoded[1].reactants[0] is decoded[1].products[0]
    assert decoded[0].products[0] is decoded[1].reactants[1]
    assert decoded[0].reactants[0] is not decoded[0].products[0]


def test_resonance_order_labels_and_atom_correspondence_are_preserved():
    species = species_of("C[CH]c1ccccc1", resonance=True)
    assert len(species.molecule) > 2
    first_atoms = [id(a) for a in species.molecule[0].atoms]
    species.molecule.reverse()  # a non-canonical order must survive
    species.molecule[1].atoms[2].label = "*7"
    reaction = TemplateReaction(reactants=[species], products=[species], family="x", template=["t"], is_forward=True)
    _, decoded = roundtrip([reaction])
    got = decoded[0].reactants[0]
    assert [m.to_adjacency_list() for m in got.molecule] == [m.to_adjacency_list() for m in species.molecule]
    assert got.molecule[1].atoms[2].label == "*7"
    assert [a.id for a in got.molecule[0].atoms] == [a.id for a in species.molecule[0].atoms]
    assert first_atoms  # keep the id list used above alive


def test_pairs_are_references_not_rebuilt_by_matching():
    species = species_of("C[CH]c1ccccc1")
    twin = species_of("C[CH]c1ccccc1")
    reaction = TemplateReaction(reactants=[species], products=[twin], family="x", template=["t"], is_forward=True)
    reaction.pairs = [(species, twin), (twin, species)]
    _, decoded = roundtrip([reaction])
    (a, b), (c, d) = decoded[0].pairs
    assert a is decoded[0].reactants[0] and b is decoded[0].products[0]
    assert c is b and d is a


def test_reverse_is_dropped_only_when_absent():
    reaction = populated_reactions()[0]
    assert reaction.reverse is not None
    _, decoded = roundtrip([reaction])
    assert decoded[0].reverse is not None and decoded[0].reverse.reverse is decoded[0]
    reaction.reverse = None
    _, decoded = roundtrip([reaction])
    assert decoded[0].reverse is None


# ---------------------------------------------------------------------------
# Encode-side refusals
# ---------------------------------------------------------------------------


def simple_reaction():
    species = species_of("C[CH]c1ccccc1")
    return TemplateReaction(reactants=[species], products=[species_of("CCc1ccccc1")], family="x",
                            template=["t"], is_forward=True)


def corrupt_conformer(reaction):
    reaction.reactants[0].conformer = object()


UNSUPPORTED = {
    "species conformer": lambda r: setattr(r.reactants[0], "conformer", __import__("rmgpy.statmech.conformer", fromlist=["Conformer"]).Conformer()),
    "species transport data": lambda r: setattr(r.reactants[0], "transport_data", object()),
    "species props": lambda r: r.reactants[0].props.update({"other": 1}),
    "molecule props": lambda r: r.reactants[0].molecule[0].props.update({"x": 1}),
    "molecule metal": lambda r: setattr(r.reactants[0].molecule[0], "metal", "Pt"),
    "atom coords": lambda r: setattr(r.reactants[0].molecule[0].atoms[0], "coords", np.array([0.0, 0.0, 0.0])),
    "atom props nested": lambda r: r.reactants[0].molecule[0].atoms[0].props.update({"x": [1]}),
    "specific collider": lambda r: setattr(r, "specific_collider", species_of("[He]")),
    "network kinetics": lambda r: setattr(r, "network_kinetics", Arrhenius(A=(1, "s^-1"))),
    "polymer stamp": lambda r: setattr(r, "polymer_chip_units", 2),
    "rate cache": lambda r: r.k_effective_cache.update({"x": 1.0}),
    "extra reaction attribute": lambda r: setattr(r, "surprise", 1),
    "labeled atoms": lambda r: r.labeled_atoms["reactants"].update({"*1": r.reactants[0].molecule[0].atoms[0]}),
    "other kinetics class": lambda r: setattr(r, "kinetics", __import__("rmgpy.kinetics", fromlist=["MultiArrhenius"]).MultiArrhenius(arrhenius=[Arrhenius(A=(1, "s^-1"))])),
    "kinetics solute": lambda r: setattr(r, "kinetics", Arrhenius(A=(1, "s^-1"), solute=object())),
    "unsupported thermo": lambda r: setattr(r.reactants[0], "thermo", __import__("rmgpy.thermo", fromlist=["NASA"]).NASA(polynomials=[], Tmin=(300, "K"), Tmax=(1000, "K"))),
    "non-finite degeneracy": lambda r: setattr(r, "_degeneracy", float("inf")),
    "molecule subclass": lambda r: setattr(r.reactants[0], "molecule", [__import__("rmgpy.molecule.fragment", fromlist=["Fragment"]).Fragment(smiles="CC")]),
}


@pytest.mark.parametrize("name", sorted(UNSUPPORTED))
def test_unsupported_populated_fields_raise_a_typed_error(name):
    reaction = simple_reaction()
    UNSUPPORTED[name](reaction)
    with pytest.raises(EncodeError):
        encode_generated_reactions([reaction])


def test_unsupported_top_level_types_are_refused():
    with pytest.raises(EncodeError):
        encode_generated_reactions((simple_reaction(),))
    with pytest.raises(EncodeError):
        encode_generated_reactions([object()])


# ---------------------------------------------------------------------------
# Hostile and malformed payloads
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def good_payload():
    return encode_generated_reactions(populated_reactions())


def tweak(payload, mutate):
    document = json.loads(payload)
    mutate(document)
    return json.dumps(document, sort_keys=True, separators=(",", ":")).encode()


def never_builds(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("an RMG object was constructed before validation finished")

    monkeypatch.setattr(generation_codec._Builder, "__init__", forbidden)


def test_good_payload_decodes(good_payload):
    assert len(decode_generated_reactions(good_payload, limits=DEFAULT_LIMITS)) == 4


def test_byte_cap_is_checked_before_any_parsing(good_payload, monkeypatch):
    def explode(*args, **kwargs):
        raise AssertionError("parsed before the byte cap")

    monkeypatch.setattr(generation_codec.json, "loads", explode)
    monkeypatch.setattr(generation_codec, "_preflight", explode)
    with pytest.raises(DecodeError):
        decode_generated_reactions(good_payload, limits=replace(DEFAULT_LIMITS, max_bytes=len(good_payload) - 1))


def test_preflight_runs_before_json_loads(monkeypatch):
    def explode(*args, **kwargs):
        raise AssertionError("json.loads touched a payload that preflight must reject")

    monkeypatch.setattr(generation_codec.json, "loads", explode)
    with pytest.raises(DecodeError):
        decode_generated_reactions(b"[" * 5000, limits=DEFAULT_LIMITS)
    with pytest.raises(DecodeError):
        decode_generated_reactions(b"[1," * 100 + b"1" + b"]" * 100, limits=replace(DEFAULT_LIMITS, max_tokens=50))


@pytest.mark.parametrize(
    "payload",
    [
        b"",
        b"[",
        b'{"a":',
        b'{"a":"unterminated}',
        b"]}",
        b"[]]",
        b'{"schema":"kmc.generated-reactions/1"}',
        b"\xff\xfe{}",
        b'{"a":NaN}',
        b'{"a":Infinity}',
        b'{"a":-Infinity}',
        b'{"a":1e400}',
        b'{"a":-1e400}',
        b'{"a":123456789012345678901234567890}',
        b'{"a":1,"a":2}',
        b"[" * 100 + b"]" * 100,
    ],
)
def test_malformed_bytes_raise_the_one_typed_error(payload, monkeypatch):
    never_builds(monkeypatch)
    with pytest.raises(DecodeError):
        decode_generated_reactions(payload, limits=DEFAULT_LIMITS)


def test_duplicate_keys_are_rejected_even_when_values_agree(good_payload, monkeypatch):
    never_builds(monkeypatch)
    text = good_payload.decode()
    duplicated = text.replace('"schema":"kmc.generated-reactions/1"', '"schema":"kmc.generated-reactions/1","schema":"kmc.generated-reactions/1"')
    assert duplicated != text
    with pytest.raises(DecodeError, match="duplicate"):
        decode_generated_reactions(duplicated.encode(), limits=DEFAULT_LIMITS)


def test_out_of_range_float_inside_a_valid_document_is_rejected(good_payload, monkeypatch):
    never_builds(monkeypatch)
    text = good_payload.decode()
    marker = '"degeneracy":2.0'
    assert marker in text
    with pytest.raises(DecodeError):
        decode_generated_reactions(text.replace(marker, '"degeneracy":1e400').encode(), limits=DEFAULT_LIMITS)
    with pytest.raises(DecodeError):
        decode_generated_reactions(text.replace(marker, '"degeneracy":NaN').encode(), limits=DEFAULT_LIMITS)


def first_edge(document):
    for index, atom in enumerate(document["atoms"]):
        if atom["edges"]:
            return index, atom["edges"][0]
    raise AssertionError("no bonded atom")


STRUCTURAL = {
    "unknown schema": lambda d: d.__setitem__("schema", "other/1"),
    "extra top-level key": lambda d: d.__setitem__("extra", []),
    "missing table": lambda d: d.pop("thermo"),
    "result reference out of range": lambda d: d["result"].append(10**6),
    "negative reference": lambda d: d["result"].append(-1),
    "boolean as reference": lambda d: d["result"].append(True),
    "reactant reference out of range": lambda d: d["reactions"][0]["reactants"].append(10**6),
    "species molecule reference out of range": lambda d: d["species"][0]["molecules"].append(10**6),
    "molecule atom reference out of range": lambda d: d["molecules"][0]["atoms"].append(10**6),
    "reverse reference out of range": lambda d: d["reactions"][0].__setitem__("reverse", 10**6),
    "kinetics reference out of range": lambda d: d["reactions"][0].__setitem__("kinetics", 10**6),
    "thermo reference out of range": lambda d: d["species"][0].__setitem__("thermo", 10**6),
    "pair endpoint out of range": lambda d: d["reactions"][0]["pairs"][0]["items"].__setitem__(0, 10**6),
    "pair with three endpoints": lambda d: d["reactions"][0]["pairs"][0]["items"].append(0),
    "float for int": lambda d: d["reactions"][0].__setitem__("index", 1.0),
    "int for float": lambda d: d["reactions"][0].__setitem__("degeneracy", 2),
    "string for bool": lambda d: d["reactions"][0].__setitem__("reversible", "yes"),
    "unknown reaction class": lambda d: d["reactions"][0].__setitem__("class", "os.system"),
    "unknown kinetics class": lambda d: d["kinetics"][0].__setitem__("class", "Marcus"),
    "unknown element": lambda d: d["atoms"][0].__setitem__("element", "Zz"),
    "unknown atom type": lambda d: d["atoms"][0].__setitem__("atomtype", "nope"),
    "bad unit string": lambda d: d["kinetics"][0]["A"].__setitem__("units", "__import__('os')"),
    "unresolvable unit string": lambda d: d["kinetics"][0]["A"].__setitem__("units", "zorkmid_per_flarg"),
    "overlong string": lambda d: d["reactions"][0].__setitem__("comment", "x" * 100000),
    "bad uncertainty type": lambda d: d["kinetics"][0]["A"].__setitem__("uncertainty_type", "??"),
    "unsupported species props": lambda d: d["species"][0]["props"].__setitem__("evil", True),
    "nested atom props": lambda d: d["atoms"][0]["props"].__setitem__("x", {"y": 1}),
    "atom bonds to itself": lambda d: d["atoms"][first_edge(d)[0]]["edges"].append([first_edge(d)[0], first_edge(d)[1][1]]),
    "duplicate neighbour": lambda d: d["atoms"][first_edge(d)[0]]["edges"].append(list(first_edge(d)[1])),
    "asymmetric bond": lambda d: d["atoms"][first_edge(d)[1][0]]["edges"].clear(),
    "bond with bad order": lambda d: d["bonds"][first_edge(d)[1][1]].__setitem__("order", 99),
    "self bond record": lambda d: d["bonds"][0].__setitem__("v2", d["bonds"][0]["v1"]),
    "duplicate bond record": lambda d: d["bonds"].append(dict(d["bonds"][0])),
    "bond joins other atoms": lambda d: d["bonds"][0].__setitem__("v1", (d["bonds"][0]["v1"] + 1) % len(d["atoms"])),
    "molecule repeats an atom": lambda d: d["molecules"][0]["atoms"].append(d["molecules"][0]["atoms"][0]),
    "empty molecule": lambda d: d["molecules"][0]["atoms"].clear(),
    "bond escapes its molecule": lambda d: d["molecules"][0]["atoms"].pop(),
    "short out of range": lambda d: d["atoms"][0].__setitem__("charge", 2**20),
    "int out of range": lambda d: d["reactions"][0].__setitem__("index", 2**40),
    "empty thermo array": lambda d: (d["thermo"][0]["Tdata"]["value_si"].clear(), d["thermo"][0]["Tdata"]["uncertainty_si"].clear()),
    "ragged thermo array": lambda d: d["thermo"][0]["Tdata"]["uncertainty_si"].pop(),
}


@pytest.mark.parametrize("name", sorted(STRUCTURAL))
def test_structural_corruption_is_rejected_before_any_object_is_built(good_payload, monkeypatch, name):
    never_builds(monkeypatch)
    corrupted = tweak(good_payload, STRUCTURAL[name])
    with pytest.raises(DecodeError):
        decode_generated_reactions(corrupted, limits=DEFAULT_LIMITS)


@pytest.mark.parametrize(
    "limit, value",
    [
        ("max_reactions", 3),
        ("max_species", 3),
        ("max_molecules", 3),
        ("max_atoms", 30),
        ("max_bonds", 30),
        ("max_kinetics", 2),
        ("max_thermo", 0),
        ("max_references", 100),
        ("max_depth", 3),
        ("max_string_length", 5),
    ],
)
def test_cumulative_count_limits_apply_across_the_reference_graph(good_payload, monkeypatch, limit, value):
    never_builds(monkeypatch)
    with pytest.raises(DecodeError):
        decode_generated_reactions(good_payload, limits=replace(DEFAULT_LIMITS, **{limit: value}))


def test_limits_argument_is_required_keyword(good_payload):
    with pytest.raises(TypeError):
        decode_generated_reactions(good_payload)  # type: ignore[call-arg]


# ---------------------------------------------------------------------------
# Parity provider (fast variant; the full PS compile is the integration proof)
# ---------------------------------------------------------------------------


def test_roundtrip_provider_returns_equal_results_and_new_objects(real_generation):
    database = real_generation.database
    captured = database.kinetics.generate_reactions_from_families(
        [species_of("C[CH]c1ccccc1"), species_of("CC(C)c1ccccc1")], only_families=["H_Abstraction"]
    )
    calls = []

    def fake_generate(reactants, products=None, only_families=None, resonance=True):
        calls.append((products, only_families, resonance))
        return captured

    provider = roundtrip_generation_provider_factory()(fake_generate)
    wrapped = provider([], only_families=["H_Abstraction"])
    assert calls == [(None, ["H_Abstraction"], True)]
    assert captured and len(captured) == len(wrapped)
    assert all(a is not b for a, b in zip(captured, wrapped))
    assert_same(captured, wrapped)
    assert_consumers_agree(captured, wrapped)


def test_discovery_through_the_roundtrip_provider_matches_cold_generation(real_generation):
    """Small-family-set variant of the full-compile parity proof (generation + discovery only)."""
    kinetics_database = real_generation.database.kinetics
    proxies = compiler.ps_proxy_set()[:3]
    cold = compiler.EventSetCompiler.discover_family_reactions(
        kinetics_database, proxies, ("R_Recombination",), family_universe=("R_Recombination",)
    )
    provider = roundtrip_generation_provider_factory()(kinetics_database.generate_reactions_from_families)
    wrapped = compiler.EventSetCompiler.discover_family_reactions(
        kinetics_database, proxies, ("R_Recombination",), family_universe=("R_Recombination",),
        generation_provider=provider,
    )
    assert cold[:2] == wrapped[:2]
    assert set(cold[2]) == set(wrapped[2])
    total = 0
    for site in cold[2]:
        assert len(cold[2][site]) == len(wrapped[2][site])
        total += len(cold[2][site])
        assert_same(list(cold[2][site]), list(wrapped[2][site]))
        assert_consumers_agree(cold[2][site], wrapped[2][site])
    assert total > 0
