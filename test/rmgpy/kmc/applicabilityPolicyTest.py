"""Behavioral contract for the persistent-carbene applicability policy."""

import copy
import importlib.util
import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

import rmgpy.kmc.compiler as compiler_module
from rmgpy.data.kinetics.common import ensure_independent_atom_ids
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.kmc.compiler import (
    EventSetCompiler,
    PERSISTENT_CARBENE_POLICY_VERSION,
    SiteProxy,
    apply_persistent_carbene_policy,
    classify_persistent_carbene,
    duplicate_transition_groups,
    _apply_record,
    _canonical_adjacency,
    _mapped_reaction_u2_roots,
    _reaction_from_record,
)
from rmgpy.kmc.event_record import EventRecord
from rmgpy.molecule.molecule import Molecule
from rmgpy.species import Species


FIXTURES = Path(__file__).with_name("fixtures")
REPO_ROOT = Path(__file__).resolve().parents[3]
DB_PATH = Path(
    os.environ.get("RMG_DATABASE_PATH", str(REPO_ROOT.parent / "RMG-database"))
)

_I084_FIXTURE_SPEC = importlib.util.spec_from_file_location(
    "i084_applicability_record", FIXTURES / "i084_applicability_record.py"
)
_I084_FIXTURE = importlib.util.module_from_spec(_I084_FIXTURE_SPEC)
_I084_FIXTURE_SPEC.loader.exec_module(_I084_FIXTURE)
_I084_QUINOID_SPEC = importlib.util.spec_from_file_location(
    "i084_quinoid_product", FIXTURES / "i084_quinoid_product.py"
)
_I084_QUINOID = importlib.util.module_from_spec(_I084_QUINOID_SPEC)
_I084_QUINOID_SPEC.loader.exec_module(_I084_QUINOID)


def _h_abstraction_record(event_id="evt_forward", reverse_of="evt_reverse"):
    return {
        "event_id": event_id,
        "reverse_of": reverse_of,
        "family": "H_Abstraction",
        "template": "synthetic",
        "status": "enabled",
        "rate_source": {"kind": "RMG family estimate", "available": True},
        "reactant_graphs": [
            """multiplicity 2
1 *1 C u1 p0 c0 {2,S} {3,S} {4,S}
2    H u0 p0 c0 {1,S}
3    H u0 p0 c0 {1,S}
4    H u0 p0 c0 {1,S}
""",
            """multiplicity 2
1 *2 N u1 p1 c0 {2,S} {3,S}
2    H u0 p0 c0 {1,S}
3    H u0 p0 c0 {1,S}
""",
        ],
        "product_graphs": [
            """multiplicity 3
1 *1 C u2 p0 c0 {2,S} {3,S}
2    H u0 p0 c0 {1,S}
3    H u0 p0 c0 {1,S}
""",
            """1 *2 N u0 p1 c0 {2,S} {3,S} {4,S}
2    H u0 p0 c0 {1,S}
3    H u0 p0 c0 {1,S}
4    H u0 p0 c0 {1,S}
""",
        ],
        "atom_map": {0: 0, 1: 1},
        "bond_ops": [
            {"action": "break", "atoms": [0, 3], "order": "1.0"},
            {"action": "form", "atoms": [4, 3], "order": "1.0"},
            {"action": "set_radical", "atom": 0, "value": 2},
            {"action": "set_radical", "atom": 4, "value": 0},
        ],
    }


def _pair():
    forward = _h_abstraction_record()
    reverse = copy.deepcopy(forward)
    reverse.update(
        event_id="evt_reverse",
        reverse_of="evt_forward",
        reactant_graphs=copy.deepcopy(forward["product_graphs"]),
        product_graphs=copy.deepcopy(forward["reactant_graphs"]),
        bond_ops=[
            {"action": "break", "atoms": [3, 6], "order": "1.0"},
            {"action": "form", "atoms": [0, 6], "order": "1.0"},
            {"action": "set_radical", "atom": 0, "value": 1},
            {"action": "set_radical", "atom": 3, "value": 1},
        ],
    )
    return [forward, reverse]


def _complete_inverse_record(record):
    """Build an inverse using the forward rewrite's atom identities."""
    reactants = [
        Molecule().from_adjacency_list(graph)
        for graph in record["reactant_graphs"]
    ]
    original_atoms = [
        atom for molecule in reactants for atom in molecule.atoms
    ]
    index_map = {}
    reverse_reactants = None
    for atom_index in range(len(original_atoms)):
        products, tracked = _apply_record(record, reactants, atom_index)
        product_atoms = [
            atom for molecule in products for atom in molecule.atoms
        ]
        index_map[atom_index] = product_atoms.index(tracked)
        if reverse_reactants is None:
            reverse_reactants = products

    old_values = {
        "set_radical": "radical_electrons",
        "set_charge": "charge",
        "set_lone_pairs": "lone_pairs",
        "set_implicit_hydrogens": "implicit_hydrogens",
    }
    inverse_operations = []
    for operation in reversed(record["bond_ops"]):
        inverse = copy.deepcopy(operation)
        action = operation["action"]
        if action in {"break", "form", "change"}:
            inverse["atoms"] = [
                index_map[index] for index in operation["atoms"]
            ]
            if action == "break":
                inverse["action"] = "form"
            elif action == "form":
                inverse["action"] = "break"
        elif action in old_values:
            inverse["atom"] = index_map[operation["atom"]]
            inverse["value"] = int(
                getattr(original_atoms[operation["atom"]], old_values[action])
            )
        inverse_operations.append(inverse)

    reverse_atoms = [
        atom for molecule in reverse_reactants for atom in molecule.atoms
    ]
    original_heavy = [
        atom for atom in original_atoms if atom.element.number != 1
    ]
    reverse_heavy = [
        atom for atom in reverse_atoms if atom.element.number != 1
    ]
    inverse_atom_map = {}
    for original_heavy_index, atom in enumerate(original_heavy):
        original_index = original_atoms.index(atom)
        reverse_atom = reverse_atoms[index_map[original_index]]
        reverse_heavy_index = reverse_heavy.index(reverse_atom)
        inverse_atom_map[reverse_heavy_index] = original_heavy_index

    inverse = copy.deepcopy(record)
    inverse.update(
        event_id=record["reverse_of"],
        reverse_of=record["event_id"],
        reactant_graphs=[
            molecule.to_adjacency_list(remove_h=False)
            for molecule in reverse_reactants
        ],
        product_graphs=copy.deepcopy(record["reactant_graphs"]),
        atom_map=inverse_atom_map,
        bond_ops=inverse_operations,
    )
    return inverse


def test_rewrite_tracks_the_copied_root_object_through_product_rebuilding():
    record = _h_abstraction_record()
    participants = [
        Molecule().from_adjacency_list(graph)
        for graph in record["reactant_graphs"]
    ]
    products, tracked_root = _apply_record(record, participants, 0)

    assert tracked_root in [
        atom for product in products for atom in product.atoms
    ]
    assert tracked_root.radical_electrons == 2


def _quinoid_record():
    return {
        "event_id": "evt_quinoid_forward",
        "reverse_of": "evt_quinoid_reverse",
        "status": "enabled",
        "rate_source": {"available": True},
        "reactant_graphs": [
            _canonical_adjacency(Molecule().from_adjacency_list(graph))
            for graph in _I084_FIXTURE.REACTANT_GRAPHS
        ],
        "product_graphs": [_I084_QUINOID.CAPTURED_PRODUCT_GRAPH],
        "template": "R7HJ_2;C_rad_out_OneDe/Cs;Cb_H_out",
        "atom_map": _I084_QUINOID.ATOM_MAP,
        "bond_ops": _I084_QUINOID.BOND_OPS,
    }


def test_real_quinoid_intra_h_migration_mapping_uses_generated_resonance_form():
    database = KineticsDatabase()
    database.load_families(
        str(DB_PATH / "input/kinetics/families"),
        families=["intra_H_migration"],
    )
    record = _quinoid_record()
    reaction = _reaction_from_record(record)
    roots = _mapped_reaction_u2_roots(
        database.families["intra_H_migration"], reaction, record
    )
    assert len(roots) == 1
    assert roots[0]["recipe_label"] == "*2"
    assert roots[0]["family_forward_role"] == "product"
    assert roots[0]["attribution"] == "verified"
    touched = compiler_module._touched_atom_indices(record["bond_ops"])
    assert roots[0]["reactant_atom_index"] in touched
    decision = classify_persistent_carbene(record, roots[0], None)
    assert decision["disposition"] == "not-applicable"
    assert decision["mapped_root"]["mapping_verified"] is True
    assert (
        decision["mapped_root"]["persistent_neutral_divalent_carbon"]
        is False
    )
    stored_atoms = [
        atom
        for graph in record["reactant_graphs"]
        for atom in Molecule().from_adjacency_list(graph).atoms
    ]
    untouched_carbon = next(
        index
        for index, atom in enumerate(stored_atoms)
        if atom.element.number == 6 and index not in touched
    )
    misplaced = {
        **roots[0],
        "reactant_atom_index": untouched_carbon,
        "mapping_verified": False,
        "mapping_error": "misplaced root",
    }
    misplaced_decision = classify_persistent_carbene(record, misplaced, None)
    assert misplaced_decision["disposition"] == (
        "retained-unresolved-applicability"
    )
    assert misplaced_decision["mapped_root"]["mapping_verified"] is False
    assert (
        misplaced_decision["mapped_root"]["persistent_neutral_divalent_carbon"]
        is None
    )

    reverse = _complete_inverse_record(record)
    published, refusal = apply_persistent_carbene_policy(
        (record, reverse), roots[0], None
    )
    assert refusal is None
    assert len(published) == 2


def test_repeated_root_mapping_reuses_family_reaction_generation():
    database = KineticsDatabase()
    database.load_families(
        str(DB_PATH / "input/kinetics/families"),
        families=["intra_H_migration"],
    )
    family = database.families["intra_H_migration"]
    record = _quinoid_record()
    cache = {}
    stats = {}

    cache_token = compiler_module._APPLICABILITY_GENERATION_CACHE.set(cache)
    stats_token = compiler_module._APPLICABILITY_GENERATION_STATS.set(stats)
    try:
        first = _mapped_reaction_u2_roots(
            family, _reaction_from_record(record), record
        )
        second = _mapped_reaction_u2_roots(
            family, _reaction_from_record(record), record
        )
    finally:
        compiler_module._APPLICABILITY_GENERATION_CACHE.reset(cache_token)
        compiler_module._APPLICABILITY_GENERATION_STATS.reset(stats_token)

    assert second == first
    assert len(cache) == 1
    assert stats[family.label]["requests"] == 2
    assert stats[family.label]["generation_calls"] == 1
    assert stats[family.label]["cache_hits"] == 1
    assert stats[family.label]["seconds"] >= 0.0


def _root():
    return {
        "reactant_atom_index": 0,
        "recipe_label": "*1",
        "family_forward_role": "product",
    }


def _source(label="*1", role="product", persistent=True, source="training 629"):
    return [
        {
            "source": source,
            "mapped_roots": [
                {
                    "recipe_label": label,
                    "family_forward_role": role,
                    "persistent_neutral_divalent_carbon": persistent,
                }
            ],
        }
    ]


def test_unsupported_persistent_carbene_refuses_both_directions():
    published, refusal = apply_persistent_carbene_policy(
        _pair(), _root(), _source(persistent=False)
    )
    assert published == []
    assert refusal["policy_version"] == PERSISTENT_CARBENE_POLICY_VERSION
    assert refusal["record_ids"] == ["evt_forward", "evt_reverse"]
    assert (
        refusal["reason"]
        == "no selected rate contributor supports persistent u2 in the mapped role"
    )
    assert refusal["mapped_root"]["recipe_label"] == "*1"
    assert refusal["rate_source"] == _source(persistent=False)


def test_prepublication_pair_is_refused_atomically_before_links_exist():
    pair = _pair()
    for record in pair:
        record["reverse_of"] = None
    published, refusal = apply_persistent_carbene_policy(
        pair, _root(), _source(persistent=False)
    )
    assert published == []
    assert refusal["record_ids"] == ["evt_forward", "evt_reverse"]


def test_compiler_applies_prepublication_refusal_before_rate_evaluation(
    monkeypatch,
):
    pair = _pair()
    structural = [
        EventRecord(
            reactant_graphs=record["reactant_graphs"],
            product_graphs=record["product_graphs"],
            atom_map=record["atom_map"],
            bond_ops=record["bond_ops"],
            rate_source={"available": True},
        )
        for record in pair
    ]
    family = SimpleNamespace(auto_generated=False)
    compiler = object.__new__(EventSetCompiler)
    compiler.kinetics_database = SimpleNamespace(families={"H_Abstraction": family})
    compiler._applicability_source_cache = {}
    compiler._applicability_refusals = []
    compiler._direction_proxy = lambda proxy, reaction: proxy
    compiler._record = lambda *args, **kwargs: structural.pop(0)
    compiler._rate_table = lambda reaction: pytest.fail(
        "rate evaluation ran before applicability refusal"
    )
    monkeypatch.setattr(compiler_module, "_reverse_view", lambda reaction: reaction)
    monkeypatch.setattr(
        compiler_module,
        "_mapped_reaction_u2_roots",
        lambda family, reaction, record: [_root()],
    )
    monkeypatch.setattr(
        compiler_module,
        "_training_source_domain",
        lambda family, reaction, cache: _source(persistent=False),
    )
    reaction = SimpleNamespace(
        family="H_Abstraction",
        kinetics=object(),
        template=["Root"],
    )
    published, refused = compiler._build_linked_family_pair(
        SiteProxy("synthetic", []), reaction, {}
    )
    assert published == []
    assert refused is None
    assert len(compiler._applicability_refusals) == 1
    assert compiler._applicability_refusals[0]["record_ids"]


@pytest.mark.parametrize("corrupt_direction", [0, 1])
def test_compiler_refuses_structurally_corrupt_pair_before_publication(
    monkeypatch, corrupt_direction
):
    pair = _pair()
    pair[corrupt_direction]["atom_map"] = {0: 1, 1: 0}
    structural = [
        EventRecord(
            reactant_graphs=record["reactant_graphs"],
            product_graphs=record["product_graphs"],
            atom_map=record["atom_map"],
            bond_ops=record["bond_ops"],
            rate_source={"available": True},
        )
        for record in pair
    ]
    family = SimpleNamespace(auto_generated=False)
    compiler = object.__new__(EventSetCompiler)
    compiler.kinetics_database = SimpleNamespace(families={"H_Abstraction": family})
    compiler._applicability_source_cache = {}
    compiler._applicability_refusals = []
    compiler._direction_proxy = lambda proxy, reaction: proxy
    compiler._record = lambda *args, **kwargs: structural.pop(0)
    compiler._rate_table = lambda reaction: pytest.fail(
        "rate evaluation ran for a structurally corrupt pair"
    )
    monkeypatch.setattr(compiler_module, "_reverse_view", lambda reaction: reaction)
    monkeypatch.setattr(
        compiler_module,
        "_mapped_reaction_u2_roots",
        lambda family, reaction, record: [_root()],
    )
    monkeypatch.setattr(
        compiler_module,
        "_training_source_domain",
        lambda family, reaction, cache: _source(),
    )
    reaction = SimpleNamespace(
        family="H_Abstraction",
        kinetics=object(),
        template=["Root"],
    )

    published, refused = compiler._build_linked_family_pair(
        SiteProxy("synthetic", []), reaction, {}
    )

    assert published == []
    assert refused is None
    assert len(compiler._applicability_refusals) == 1
    refusal = compiler._applicability_refusals[0]
    assert refusal["disposition"] == "refused-structural-inconsistency"
    assert refusal["reason"].startswith("structural-inconsistency:")
    assert len(refusal["record_ids"]) == 2


def test_supported_same_role_is_retained_with_transfer_caveat():
    published, refusal = apply_persistent_carbene_policy(_pair(), _root(), _source())
    assert refusal is None
    assert len(published) == 2
    assert {record["reverse_of"] for record in published} == {
        "evt_forward",
        "evt_reverse",
    }
    annotations = [record["rate_source"]["applicability"] for record in published]
    assert all(
        item["disposition"] == "retained-supported-transfer" for item in annotations
    )
    assert all(
        item["caveat"] == "same reacting role; donor/acceptor environment may differ"
        for item in annotations
    )


def test_unserialized_source_domain_is_retained_as_unresolved():
    published, refusal = apply_persistent_carbene_policy(_pair(), _root(), None)
    assert refusal is None
    assert len(published) == 2
    assert all(
        record["rate_source"]["applicability"]["disposition"]
        == "retained-unresolved-applicability"
        for record in published
    )


def test_u2_contributor_in_wrong_recipe_role_does_not_authorize_transfer():
    published, refusal = apply_persistent_carbene_policy(
        _pair(), _root(), _source(label="*2", role="reactant")
    )
    assert published == []
    assert refusal["disposition"] == "refused-unsupported-transfer"


def test_reverse_stored_training_source_flips_to_family_forward_role(monkeypatch):
    class Entry:
        index = 629
        label = "reverse training"
        item = SimpleNamespace(reactants=["forward"], products=["reverse"])

    class Family:
        auto_generated = False
        label = "H_Abstraction"

        @staticmethod
        def extract_source_from_comments(reaction):
            return True, ("Training", Entry(), True)

    monkeypatch.setattr(
        compiler_module,
        "_mapped_reaction_u2_roots",
        lambda family, reaction: [
            {
                "recipe_label": "*1",
                "family_forward_role": (
                    "reactant" if reaction.reactants == ["reverse"] else "product"
                ),
                "persistent_neutral_divalent_carbon": True,
            }
        ],
    )
    domain = compiler_module._training_source_domain(Family(), object(), {})
    assert domain[0]["source"]["stored_reverse"] is True
    assert domain[0]["mapped_roots"][0]["family_forward_role"] == "reactant"


@pytest.mark.parametrize(
    "permute_reactants,permute_products",
    [(False, False), (True, False), (False, True), (True, True)],
)
def test_real_training_629_forward_and_reverse_bind_the_same_supported_source(
    permute_reactants, permute_products
):
    database = KineticsDatabase()
    database.load_families(
        str(DB_PATH / "input/kinetics/families"),
        families=["H_Abstraction"],
    )
    family = database.families["H_Abstraction"]
    entry = family.get_training_depository().entries[629]
    reactants = copy.deepcopy(entry.item.reactants)
    ensure_independent_atom_ids(reactants, resonance=False)
    generated = family.generate_reactions(
        [species.molecule for species in reactants],
        prod_resonance=True,
        delete_labels=False,
        relabel_atoms=True,
    )
    generated = next(
        reaction
        for reaction in generated
        if any(
            atom.element.number == 6 and atom.radical_electrons == 2
            for product in reaction.products
            for atom in product.atoms
        )
    )
    base = TemplateReaction(
        reactants=[Species(molecule=[molecule]) for molecule in generated.reactants],
        products=[Species(molecule=[molecule]) for molecule in generated.products],
        family=family.label,
        template=generated.template,
        degeneracy=generated.degeneracy,
    )
    # A record estimated in the CH2 + NH3 direction carries that direction's
    # family template, exactly as the compiler stores ``estimate.template``.
    reverse_template = next(
        reaction.template
        for reaction in family.generate_reactions(
            [molecule.copy(deep=True) for molecule in generated.products],
            prod_resonance=True,
        )
        if all(
            any(product.is_isomorphic(reactant) for reactant in generated.reactants)
            for product in reaction.products
        )
    )
    dispositions = []
    recipe_labels = []
    roles = []
    complete_records = []
    mapped_roots = []
    source_domains = []
    for reverse in (False, True):
        reaction = TemplateReaction(
            reactants=copy.deepcopy(base.products if reverse else base.reactants),
            products=copy.deepcopy(base.reactants if reverse else base.products),
            kinetics=copy.deepcopy(entry.data),
            family=family.label,
            template=reverse_template if reverse else base.template,
            degeneracy=base.degeneracy,
        )
        if permute_reactants:
            reaction.reactants.reverse()
        if permute_products:
            reaction.products.reverse()
        reaction.kinetics.comment = f"Matched reaction 629 {entry.label} in training"
        compiler = object.__new__(EventSetCompiler)
        compiler.reference_thermo_provider = SimpleNamespace(provenance={})
        proxy = SiteProxy(
            "training-629",
            reaction.reactants,
            participant_site_types=["training-629"] * len(reaction.reactants),
        )
        record = compiler._record(
            proxy,
            reaction,
            {},
            rate_override=(None, {"available": True}),
        )
        record = record.to_dict()
        assert sum(
            atom.element.number != 1
            for graph in record["reactant_graphs"]
            for atom in Molecule().from_adjacency_list(graph).atoms
        ) == 2
        for side in ("reactant_graphs", "product_graphs"):
            record[side] = [
                graph.replace("*1", "*swap")
                .replace("*3", "*1")
                .replace("*swap", "*3")
                for graph in record[side]
            ]
        roots = compiler_module._mapped_reaction_u2_roots(family, reaction, record)
        assert len(roots) == 1
        complete_records.append(record)
        mapped_roots.append(roots[0])
        stored_atoms = [
            atom
            for graph in record["reactant_graphs"]
            for atom in Molecule().from_adjacency_list(graph).atoms
        ]
        assert (
            stored_atoms[roots[0]["reactant_atom_index"]].label
            != roots[0]["recipe_label"]
        )
        source_domain = compiler_module._training_source_domain(
            family, reaction, {}
        )
        source_domains.append(source_domain)
        decision = classify_persistent_carbene(record, roots[0], source_domain)
        dispositions.append(decision["disposition"])
        recipe_labels.append(decision["mapped_root"]["recipe_label"])
        roles.append(decision["mapped_root"]["family_forward_role"])
        if not reverse:
            wrong_rewrite = copy.deepcopy(record)
            radical_ops = [
                operation
                for operation in wrong_rewrite["bond_ops"]
                if operation["action"] == "set_radical"
            ]
            assert sorted(operation["value"] for operation in radical_ops) == [0, 2]
            for operation in radical_ops:
                operation["value"] = 2 - operation["value"]
            wrong_roots = compiler_module._mapped_reaction_u2_roots(
                family, reaction, wrong_rewrite
            )
            assert len(wrong_roots) == 1
            wrong_rewrite.update(
                event_id="evt_corrupt_forward",
                reverse_of="evt_corrupt_reverse",
            )
            wrong_reverse = copy.deepcopy(wrong_rewrite)
            wrong_reverse.update(
                event_id="evt_corrupt_reverse",
                reverse_of="evt_corrupt_forward",
            )
            wrong_published, wrong_decision = apply_persistent_carbene_policy(
                (wrong_rewrite, wrong_reverse),
                wrong_roots[0],
                compiler_module._training_source_domain(family, reaction, {}),
            )
            assert wrong_published == []
            assert wrong_decision["disposition"] == ("refused-structural-inconsistency")
            assert wrong_decision["record_ids"] == [
                "evt_corrupt_forward",
                "evt_corrupt_reverse",
            ]
            assert wrong_decision["mapped_root"]["mapping_verified"] is False
            assert wrong_decision["reason"].startswith("structural-inconsistency:")

    assert dispositions == [
        "retained-supported-transfer",
        "retained-supported-transfer",
    ]
    assert recipe_labels == ["*1", "*3"]
    assert roles == ["product", "reactant"]
    complete_records[0]["reverse_of"] = complete_records[1]["event_id"]
    complete_records[1]["reverse_of"] = complete_records[0]["event_id"]
    published, refusal = apply_persistent_carbene_policy(
        complete_records,
        mapped_roots[0],
        source_domains[0],
    )
    assert refusal is None
    assert len(published) == 2


def test_mixed_rule_training_source_is_retained_unresolved_as_a_whole(
    monkeypatch,
):
    class Entry:
        index = 629
        label = "training"
        item = object()

    class Rule:
        index = 12
        label = "averaged rule"

    class Family:
        auto_generated = False
        label = "H_Abstraction"

        @staticmethod
        def extract_source_from_comments(reaction):
            return False, (
                "Rate Rules",
                {
                    "training": [(Rule(), Entry(), 0.75)],
                    "rules": [(Rule(), 0.25)],
                },
            )

    monkeypatch.setattr(
        compiler_module,
        "_mapped_reaction_u2_roots",
        lambda family, reaction: [
            {
                "recipe_label": "*1",
                "family_forward_role": "product",
                "persistent_neutral_divalent_carbon": True,
            }
        ],
    )
    domain = compiler_module._training_source_domain(Family(), object(), {})
    assert [item["source"]["kind"] for item in domain] == ["training", "rule"]
    assert domain[1]["mapped_roots"] is None
    published, refusal = apply_persistent_carbene_policy(_pair(), _root(), domain)
    assert refusal is None
    assert len(published) == 2
    assert published[0]["rate_source"]["applicability"]["disposition"] == (
        "retained-unresolved-applicability"
    )


def test_nine_form_delocalised_witness_has_no_carbene_reacting_centre():
    witness = json.loads((FIXTURES / "i078_resonance_witness.json").read_text())
    database = KineticsDatabase()
    database.load_families(
        str(DB_PATH / "input/kinetics/families"),
        families=["H_Abstraction"],
    )
    family = database.families["H_Abstraction"]
    reactants = compiler_module._molecules_from_graphs(witness["reactant_graphs"])
    compiler_module._tag_atom_ids(reactants)
    tracked, _ = _apply_record(witness, reactants)
    attributions, failure = compiler_module._verified_attributions(
        family, witness, reactants, tracked, {"*1", "*2", "*3"}, {}, {}
    )
    # The stored product draws u2 on atom 74 only as one of nine resonance
    # forms; the verified recipe centre is *1 = 96, which never carries u2.
    assert failure is None
    assert attributions == {((25, "*3"), (72, "*2"), (96, "*1"))}
    assert (
        _mapped_reaction_u2_roots(family, _reaction_from_record(witness), witness)
        == []
    )
    assert (
        witness["event_id"]
        == "evt_453febdaa6658ad2fa22f6375f21b3c46d0eb8ac92c534f69fb7dfafb72665b5"
    )
    stored_root = compiler_module._stored_u2_root_fallback(witness)
    assert [root["reactant_atom_index"] for root in stored_root] == [74]
    decision = classify_persistent_carbene(witness, stored_root[0], _source())
    assert decision["disposition"] == "not-applicable"
    assert decision["mapped_root"]["product_rewrite_verified"] is True
    assert decision["mapped_root"]["resonance_form_count"] == 9
    assert decision["mapped_root"]["persistent_neutral_divalent_carbon"] is False
    reverse = _complete_inverse_record(witness)
    published, refusal = apply_persistent_carbene_policy(
        (witness, reverse), stored_root[0], _source()
    )
    assert refusal is None
    assert len(published) == 2
    assert "applicability" not in published[0]["rate_source"]


def test_wrong_recipe_root_is_retained_only_as_unresolved():
    wrong_root = {
        "reactant_atom_index": 1,
        "recipe_label": "*1",
        "family_forward_role": "product",
    }
    decision = classify_persistent_carbene(_pair()[0], wrong_root, _source())
    assert decision["disposition"] == "retained-unresolved-applicability"
    assert decision["mapped_root"]["mapping_verified"] is False
    assert "not touched" in decision["reason"]


@pytest.mark.parametrize(
    "root_index, mapping_error",
    [
        (-1, "no family-generated candidate matched the stored template"),
        (0, "no template-matched candidate reproduces the stored rewrite"),
    ],
)
def test_unverified_root_is_retained_only_when_every_role_is_supported(
    root_index, mapping_error
):
    root = {
        **_root(),
        "reactant_atom_index": root_index,
        "mapping_verified": False,
        "mapping_error": mapping_error,
    }
    decision = classify_persistent_carbene(_pair()[0], root, _source())
    assert decision["disposition"] == "retained-unresolved-applicability"
    assert decision["reason"] == mapping_error
    assert decision["attribution"] == "unverified"


@pytest.mark.parametrize("root_index", [-1, 0])
@pytest.mark.parametrize("admissible_roles", [None, "declared"])
def test_unverified_root_cannot_hide_an_unsupported_role(root_index, admissible_roles):
    root = {
        **_root(),
        "reactant_atom_index": root_index,
        "mapping_verified": False,
        "mapping_error": "no family-generated candidate matched the stored template",
    }
    if admissible_roles is None:
        root.update(
            recipe_label="__runtime_recipe_root_unresolved__",
            admissible_roles=None,
        )
    published, refusal = apply_persistent_carbene_policy(
        _pair(), root, _source(persistent=False)
    )
    assert published == []
    assert refusal["disposition"] == "refused-unsupported-transfer"
    assert refusal["reason"].startswith("stored-transition attribution unverified")
    assert refusal["record_ids"] == ["evt_forward", "evt_reverse"]


def test_unverified_any_role_is_refused_even_when_one_role_is_supported():
    root = {
        **_root(),
        "recipe_label": "__runtime_recipe_root_unresolved__",
        "admissible_roles": None,
        "mapping_verified": False,
        "mapping_error": "no family-generated candidate matched the stored template",
    }
    published, refusal = apply_persistent_carbene_policy(_pair(), root, _source())
    assert published == []
    assert refusal["disposition"] == "refused-unsupported-transfer"


@pytest.mark.parametrize("corrupt_direction", [0, 1])
def test_unbound_family_root_does_not_hide_an_invalid_stored_rewrite(
    corrupt_direction,
):
    pair = _pair()
    radical_operation = next(
        operation
        for operation in pair[corrupt_direction]["bond_ops"]
        if operation["action"] == "set_radical"
    )
    radical_operation["value"] += 1
    published, decision = apply_persistent_carbene_policy(
        pair,
        {
            **_root(),
            "reactant_atom_index": -1,
            "mapping_verified": False,
            "mapping_error": "family root is unbound",
        },
        _source(),
    )
    assert published == []
    assert decision["disposition"] == "refused-structural-inconsistency"
    assert decision["record_ids"] == ["evt_forward", "evt_reverse"]
    assert decision["reason"].startswith("structural-inconsistency:")
    assert "stored rewrite" in decision["reason"]


@pytest.mark.parametrize(
    "root_index, label, role",
    [(0, "*1", "product"), (4, "*2", "reactant")],
)
def test_corrupt_atom_map_refuses_both_directions(root_index, label, role):
    record = _pair()[0]
    record["atom_map"] = {0: 1, 1: 0}
    root = {
        "reactant_atom_index": root_index,
        "recipe_label": label,
        "family_forward_role": role,
    }
    pair = _pair()
    pair[0] = record
    published, refusal = apply_persistent_carbene_policy(
        pair, root, _source(label=label, role=role)
    )
    assert published == []
    assert refusal["disposition"] == "refused-structural-inconsistency"
    assert refusal["record_ids"] == ["evt_forward", "evt_reverse"]
    assert refusal["mapped_root"]["mapping_verified"] is False
    assert refusal["reason"].startswith("structural-inconsistency:")
    assert "element" in refusal["reason"]


def test_persistent_u2_with_remote_spectator_radical_is_classified_per_centre():
    record = {
        "reactant_graphs": [
            """multiplicity 3
1 *1 C u1 p0 c0 {2,S} {3,S} {9,S}
2    H u0 p0 c0 {1,S}
3    C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4    H u0 p0 c0 {3,S}
5    H u0 p0 c0 {3,S}
6    C u1 p0 c0 {3,S} {7,S} {8,S}
7    H u0 p0 c0 {6,S}
8    H u0 p0 c0 {6,S}
9    H u0 p0 c0 {1,S}
""",
            """multiplicity 2
1 *2 N u1 p1 c0 {2,S} {3,S}
2    H u0 p0 c0 {1,S}
3    H u0 p0 c0 {1,S}
""",
        ],
        "product_graphs": [
            """multiplicity 4
1 *1 C u2 p0 c0 {2,S} {3,S}
2    H u0 p0 c0 {1,S}
3    C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4    H u0 p0 c0 {3,S}
5    H u0 p0 c0 {3,S}
6    C u1 p0 c0 {3,S} {7,S} {8,S}
7    H u0 p0 c0 {6,S}
8    H u0 p0 c0 {6,S}
""",
            """1 *2 N u0 p1 c0 {2,S} {3,S} {4,S}
2    H u0 p0 c0 {1,S}
3    H u0 p0 c0 {1,S}
4    H u0 p0 c0 {1,S}
""",
        ],
        "atom_map": {0: 0, 1: 1, 2: 2, 3: 3},
        "bond_ops": [
            {"action": "break", "atoms": [0, 8], "order": "1.0"},
            {"action": "form", "atoms": [9, 8], "order": "1.0"},
            {"action": "set_radical", "atom": 0, "value": 2},
            {"action": "set_radical", "atom": 9, "value": 0},
        ],
    }
    decision = classify_persistent_carbene(record, _root(), _source())
    assert decision["disposition"] == "retained-supported-transfer"
    assert decision["mapped_root"]["persistent_neutral_divalent_carbon"] is True
    assert decision["mapped_root"]["resonance_form_count"] >= 1


def test_disproportionation_fallback_root_can_never_be_supported():
    fallback = compiler_module._stored_u2_root_fallback(_pair()[0])
    assert len(fallback) == 1
    decision = classify_persistent_carbene(
        _pair()[0],
        fallback[0],
        _source(label="__runtime_recipe_root_unresolved__", role="product"),
    )
    assert decision["disposition"] == "refused-unsupported-transfer"
    assert "family relabeling failed" in decision["reason"]
    uncalibrated = classify_persistent_carbene(_pair()[0], fallback[0], None)
    assert uncalibrated["disposition"] == "retained-unresolved-applicability"
    assert "family relabeling failed" in uncalibrated["reason"]


def test_ordinary_radical_channel_is_byte_unchanged():
    pair = _pair()
    before = json.dumps(pair, sort_keys=True)
    ordinary_root = {
        "reactant_atom_index": 4,
        "recipe_label": "*2",
        "family_forward_role": "product",
    }
    published, refusal = apply_persistent_carbene_policy(pair, ordinary_root, _source())
    assert refusal is None
    assert json.dumps(published, sort_keys=True) == before


def test_retained_pair_remains_reciprocal():
    published, _ = apply_persistent_carbene_policy(_pair(), _root(), _source())
    assert len(published) == 2
    by_id = {record["event_id"]: record for record in published}
    assert all(
        by_id[record["reverse_of"]]["reverse_of"] == record["event_id"]
        for record in published
    )


def _transition_pair(prefix, junction_kind=None, rate=1.0):
    forward, reverse = _pair()
    forward["event_id"], reverse["event_id"] = f"{prefix}_f", f"{prefix}_r"
    forward["reverse_of"], reverse["reverse_of"] = (
        reverse["event_id"],
        forward["event_id"],
    )
    for record in (forward, reverse):
        record["k_table"] = {"T": [1000.0], "k": [rate]}
        record["junction_ops"] = None
    if junction_kind:
        label = junction_kind.rsplit("_", 1)[-1]
        forward["junction_ops"] = [
            {
                "action": "create",
                "junction_kind": junction_kind,
                "attacked_atom_label": label,
                "chosen_resonance_localisation": f"ortho:{label}",
                "reflection": None if label == "S7" else {"S6": "S7", "S7": "S6"},
                "reverse_event_handle": reverse["event_id"],
            }
        ]
        reverse["junction_ops"] = [
            {
                **forward["junction_ops"][0],
                "action": "dissociate",
                "reverse_event_handle": forward["event_id"],
            }
        ]
    return [forward, reverse]


def test_s6_s7_split_rates_and_reverse_histories_are_not_duplicates():
    records = _transition_pair("s6", "J_ortho_S6", 0.5) + _transition_pair(
        "s7", "J_ortho_S7", 0.5
    )
    assert sum(record["k_table"]["k"][0] for record in records[::2]) == pytest.approx(
        1.0
    )
    assert duplicate_transition_groups(records) == []
    by_id = {record["event_id"]: record for record in records}
    assert all(
        by_id[record["reverse_of"]]["reverse_of"] == record["event_id"]
        for record in records
    )


def test_synthetic_duplicated_reciprocal_pair_is_detected():
    records = _transition_pair("first", rate=1.0) + _transition_pair(
        "duplicate", rate=2.0
    )
    assert duplicate_transition_groups(records) == [
        [("duplicate_f", "duplicate_r"), ("first_f", "first_r")]
    ]


# --- I-088: verified stored-transition attribution -------------------------

_CENSUS = json.loads((FIXTURES / "i088_carbene_census_records.json").read_text())


def _pairs(records):
    by_id = {record["event_id"]: record for record in records}
    return [
        (record, by_id[record["reverse_of"]])
        for record in records
        if record["rate_source"]["kind"] == "RMG family estimate"
    ]


@pytest.fixture(scope="module")
def census_database():
    from rmgpy.data.rmg import RMGDatabase

    database = RMGDatabase()
    database.load_kinetics(
        str(DB_PATH / "input/kinetics"),
        reaction_libraries=[],
        seed_mechanisms=None,
        kinetics_families=["H_Abstraction", "Disproportionation", "R_Recombination"],
        kinetics_depositories=["training"],
    )
    database.load_thermo(
        str(DB_PATH / "input/thermo"),
        thermo_libraries=["primaryThermoLibrary"],
        depository=True,
    )
    compiler_module.prepare_rate_rules(
        database.kinetics, database.thermo, verbose=True
    )
    return database


def _census_decision(database, forward, reverse):
    """Run the production attribution, source domain and policy on one pair."""
    from rmgpy.kinetics import Arrhenius

    family = database.kinetics.families[forward["family"]]
    roots = _mapped_reaction_u2_roots(
        family, _reaction_from_record(forward), forward
    )
    estimate = _reaction_from_record(forward)
    estimate.kinetics = Arrhenius(
        A=(1.0, "s^-1"),
        n=0.0,
        Ea=(0.0, "J/mol"),
        T0=(1.0, "K"),
        comment=forward["rate_source"].get("comment", ""),
    )
    domain = compiler_module._training_source_domain(family, estimate, {})
    results = [
        apply_persistent_carbene_policy((forward, reverse), root, domain)
        for root in roots
    ]
    return roots, domain, results


@pytest.mark.parametrize(
    "pair",
    _pairs(_CENSUS["refused_pairs"]),
    ids=lambda pair: pair[0]["template"] + ":" + pair[0]["event_id"][4:12],
)
def test_census_u05_u07_pairs_are_refused_on_verified_roots(census_database, pair):
    forward, reverse = pair
    roots, domain, results = _census_decision(census_database, forward, reverse)
    assert len(roots) == 1
    assert roots[0]["attribution"] == "verified"
    assert roots[0]["recipe_label"] == "*1"
    assert roots[0]["family_forward_role"] == "product"
    assert roots[0]["attribution_count"] == 1
    assert domain and all(
        contributor["source"]["kind"] == "training" for contributor in domain
    )
    (published, refusal), = results
    assert published == []
    assert refusal["disposition"] == "refused-unsupported-transfer"
    assert refusal["reason"] == (
        "no selected rate contributor supports persistent u2 in the mapped role"
    )
    assert refusal["record_ids"] == [forward["event_id"], reverse["event_id"]]
    assert refusal["mapped_root"]["persistent_neutral_divalent_carbon"] is True


@pytest.mark.parametrize(
    "pair",
    _pairs(_CENSUS["u08_pairs"]),
    ids=lambda pair: pair[0]["event_id"][4:12],
)
def test_census_u08_training_629_pairs_are_supported(census_database, pair):
    forward, reverse = pair
    roots, domain, results = _census_decision(census_database, forward, reverse)
    assert [
        (root["recipe_label"], root["family_forward_role"], root["attribution"])
        for root in roots
    ] == [("*1", "product", "verified")]
    assert [contributor["source"]["entry_index"] for contributor in domain] == [629]
    assert domain[0]["mapped_roots"][0]["attribution_orbits"] == 1
    (published, refusal), = results
    assert refusal is None
    decisions = [record["rate_source"]["applicability"] for record in published]
    assert [decision["disposition"] for decision in decisions] == [
        "retained-supported-transfer"
    ] * 2
    assert decisions[0]["caveat"] == (
        "same reacting role; donor/acceptor environment may differ"
    )


def test_real_disproportionation_shard_pair_is_unresolved_only_by_source_domain(
    census_database,
):
    (forward, reverse), = _pairs(_CENSUS["disproportionation_pair"])
    roots, domain, results = _census_decision(census_database, forward, reverse)
    assert domain is None
    assert len(roots) == 1
    root = roots[0]
    assert (root["recipe_label"], root["family_forward_role"]) == ("*1", "reactant")
    assert root["attribution"] == "verified"
    # Two resonance-equivalent *3 choices, related by the ring mirror that
    # preserves the whole mapped transition.
    assert root["attribution_count"] == 2
    assert root["attribution_orbits"] == 1
    (published, refusal), = results
    assert refusal is None
    decision = published[0]["rate_source"]["applicability"]
    assert decision["disposition"] == "retained-unresolved-applicability"
    assert decision["reason"] == (
        "selected rate has no serialized calibrated contributor domain"
    )
    assert decision["mapped_root"]["mapping_verified"] is True
    assert decision["mapped_root"]["persistent_neutral_divalent_carbon"] is True


def test_real_recombination_binds_both_repeated_star_labels(census_database):
    family = census_database.kinetics.families["R_Recombination"]
    # A remote spectator carbene rides on the radical end that recombines.
    carbene_radical = Molecule().from_adjacency_list(
        """multiplicity 4
1 C u2 p0 c0 {2,S} {4,S}
2 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
3 C u1 p0 c0 {2,S} {7,S} {8,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
"""
    )
    methyl = Molecule().from_smiles("[CH3]")
    participants = [Species(molecule=[carbene_radical]), Species(molecule=[methyl])]
    ensure_independent_atom_ids(participants, resonance=False)
    generated = [
        reaction
        for reaction in family.generate_reactions(
            [species.molecule for species in participants], prod_resonance=True
        )
        if len(reaction.products) == 1
    ]
    assert generated
    reaction = TemplateReaction(
        reactants=[Species(molecule=[molecule]) for molecule in generated[0].reactants],
        products=[Species(molecule=[molecule]) for molecule in generated[0].products],
        family=family.label,
        template=generated[0].template,
        degeneracy=generated[0].degeneracy,
    )
    compiler = object.__new__(EventSetCompiler)
    compiler.reference_thermo_provider = SimpleNamespace(provenance={})
    proxy = SiteProxy(
        "recombination",
        reaction.reactants,
        participant_site_types=["recombination"] * len(reaction.reactants),
    )
    record = compiler._record(
        proxy, reaction, {}, rate_override=(None, {"available": True})
    ).to_dict()
    reactants = compiler_module._molecules_from_graphs(record["reactant_graphs"])
    compiler_module._tag_atom_ids(reactants)
    tracked, _ = _apply_record(record, reactants)
    attributions, failure = compiler_module._verified_attributions(
        family, record, reactants, tracked, {"*"}, {}, {}
    )
    assert failure is None
    (attribution,) = attributions
    atoms = [atom for molecule in reactants for atom in molecule.atoms]
    # Both centres carry the repeated recipe label; neither is the carbene.
    assert [label for _, label in attribution] == ["*", "*"]
    assert sorted(atoms[index].radical_electrons for index, _ in attribution) == [1, 1]
    assert {index for index, _ in attribution} <= compiler_module._touched_atom_indices(
        record["bond_ops"]
    )
    assert any(atom.radical_electrons == 2 for atom in atoms)
    assert _mapped_reaction_u2_roots(family, _reaction_from_record(record), record) == []


def test_wrong_template_is_never_upgraded_by_the_stored_rewrite():
    database = KineticsDatabase()
    database.load_families(
        str(DB_PATH / "input/kinetics/families"),
        families=["H_Abstraction"],
    )
    family = database.families["H_Abstraction"]
    record = _pair()[0]
    roots = _mapped_reaction_u2_roots(family, _reaction_from_record(record), record)
    # The rewrite alone fixes *1/*2/*3, but no generated candidate carries the
    # stored template, so nothing is verified and no role is admitted.
    assert [
        (root["reactant_atom_index"], root["attribution"], root["admissible_roles"])
        for root in roots
    ] == [(0, "unverified", None)]
    assert roots[0]["mapping_error"] == (
        "no family-generated candidate matched the stored template"
    )
    published, refusal = apply_persistent_carbene_policy(_pair(), roots[0], _source())
    assert published == []
    assert refusal["disposition"] == "refused-unsupported-transfer"


def _disagreeing_attributions():
    return [((0, "*1"), (3, "*2"), (4, "*3")), ((0, "*3"), (3, "*2"), (4, "*1"))]


def test_ambiguity_is_declared_when_orbits_disagree_on_the_u2_role():
    record = _pair()[0]
    reactants = compiler_module._molecules_from_graphs(record["reactant_graphs"])
    products = compiler_module._molecules_from_graphs(record["product_graphs"])
    attributions = _disagreeing_attributions()
    roots = compiler_module._u2_roots_from_attributions(
        record, reactants, products, set(attributions), attributions
    )
    assert len(roots) == 1
    root = roots[0]
    assert root["attribution"] == "ambiguous"
    assert root["admissible_roles"] == [
        {"recipe_label": "*1", "family_forward_role": "product"},
        {"recipe_label": "*3", "family_forward_role": "product"},
    ]
    published, refusal = apply_persistent_carbene_policy(_pair(), root, _source())
    assert published == []
    assert refusal["disposition"] == "refused-unsupported-transfer"
    assert refusal["reason"].startswith("ambiguous stored-transition attribution")
    assert "*3/product" in refusal["reason"]
    both = _source()
    both[0]["mapped_roots"].append(
        {
            "recipe_label": "*3",
            "family_forward_role": "product",
            "persistent_neutral_divalent_carbon": True,
        }
    )
    published, refusal = apply_persistent_carbene_policy(_pair(), root, both)
    assert refusal is None
    decision = published[0]["rate_source"]["applicability"]
    assert decision["disposition"] == "retained-supported-transfer"
    assert decision["attribution"] == "ambiguous"
    assert "attribution_caveat" in decision


def test_disagreement_on_non_u2_labels_is_not_ambiguity():
    record = _pair()[0]
    reactants = compiler_module._molecules_from_graphs(record["reactant_graphs"])
    products = compiler_module._molecules_from_graphs(record["product_graphs"])
    attributions = [((0, "*1"), (1, "*2"), (4, "*3")), ((0, "*1"), (3, "*2"), (4, "*3"))]
    roots = compiler_module._u2_roots_from_attributions(
        record, reactants, products, set(attributions), attributions
    )
    assert [(root["recipe_label"], root["attribution"]) for root in roots] == [
        ("*1", "verified")
    ]


def test_symmetric_recombination_attributions_collapse_to_one_orbit():
    methyl = """multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
"""
    record = {
        "reactant_graphs": [methyl, methyl],
        "product_graphs": [],
        "atom_map": {0: 0, 1: 1},
        "bond_ops": [
            {"action": "form", "atoms": [0, 4], "order": "1.0"},
            {"action": "set_radical", "atom": 0, "value": 0},
            {"action": "set_radical", "atom": 4, "value": 0},
        ],
    }
    reactants = compiler_module._molecules_from_graphs(record["reactant_graphs"])
    compiler_module._tag_atom_ids(reactants)
    tracked, _ = _apply_record(record, reactants)
    swapped = [((0, "*1"), (4, "*2")), ((0, "*2"), (4, "*1"))]
    assert compiler_module._collapse_automorphic_attributions(
        swapped, reactants, tracked
    ) == [swapped[0]]
    # Moving a role onto an atom the transition does not touch is no symmetry.
    moved = [((0, "*1"), (4, "*2")), ((1, "*1"), (4, "*2"))]
    assert compiler_module._collapse_automorphic_attributions(
        moved, reactants, tracked
    ) == sorted(moved)


def test_generation_cache_is_a_bounded_lru(monkeypatch):
    calls = []
    family = SimpleNamespace(
        label="fake",
        _generate_reactions=lambda forms, **kwargs: calls.append(forms) or [],
    )
    monkeypatch.setattr(compiler_module, "APPLICABILITY_GENERATION_CACHE_LIMIT", 2)
    cache, stats = {}, {}
    for key in ("a", "b", "a", "c", "b"):
        compiler_module._generated_forward_transitions(
            family, {"reactant_graphs": [key]}, [[key]], cache, stats
        )
    assert list(cache) == [("fake", ("c",)), ("fake", ("b",))]
    assert stats["fake"]["requests"] == 5
    assert stats["fake"]["cache_hits"] == 1
    assert stats["fake"]["generation_calls"] == 4
    assert stats["fake"]["evictions"] == 2


def test_unresolved_ledger_entries_carry_the_published_final_ids(monkeypatch):
    pair = _pair()
    structural = [
        EventRecord(
            reactant_graphs=record["reactant_graphs"],
            product_graphs=record["product_graphs"],
            atom_map=record["atom_map"],
            bond_ops=record["bond_ops"],
            rate_source={"available": True},
        )
        for record in pair
    ]
    structural_ids = {record.event_id for record in structural}
    grid = compiler_module.DEFAULT_T_GRID
    compiler = object.__new__(EventSetCompiler)
    compiler.kinetics_database = SimpleNamespace(
        families={"Disproportionation": SimpleNamespace(auto_generated=True)}
    )
    compiler._applicability_source_cache = {}
    compiler._applicability_refusals = []
    compiler.temperature_grid = grid
    compiler.use_plpsec_library = False
    compiler.reference_thermo_provider = SimpleNamespace(
        provenance={},
        equilibrium_constants=lambda reaction, temperatures: [1.0] * len(temperatures),
    )
    compiler._direction_proxy = lambda proxy, reaction: proxy
    compiler._record = lambda *args, **kwargs: structural.pop(0)
    compiler._rate_table = lambda reaction: (
        {"T": list(grid), "k": [1.0] * len(grid)},
        {"kind": "RMG family estimate", "available": True, "units": "m^3/(mol*s)"},
    )
    monkeypatch.setattr(compiler_module, "_reverse_view", lambda reaction: reaction)
    monkeypatch.setattr(
        compiler_module,
        "_mapped_reaction_u2_roots",
        lambda family, reaction, record: [_root()],
    )
    monkeypatch.setattr(
        compiler_module, "_training_source_domain", lambda family, reaction, cache: None
    )
    reaction = SimpleNamespace(
        family="Disproportionation",
        kinetics=object(),
        template=["Root"],
        reactants=[],
        products=[],
    )
    published, _ = compiler._build_linked_family_pair(
        SiteProxy("synthetic", []), reaction, {}
    )
    assert len(published) == 2
    (entry,) = compiler._applicability_refusals
    assert entry["disposition"] == "retained-unresolved-applicability"
    assert entry["record_ids"] == [record.event_id for record in published]
    assert not structural_ids & set(entry["record_ids"])
    assert all(
        record.rate_source["applicability"]["disposition"]
        == "retained-unresolved-applicability"
        for record in published
    )
