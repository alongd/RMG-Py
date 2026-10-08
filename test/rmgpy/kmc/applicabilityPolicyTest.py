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


def test_real_quinoid_intra_h_migration_mapping_uses_generated_resonance_form():
    database = KineticsDatabase()
    database.load_families(
        str(DB_PATH / "input/kinetics/families"),
        families=["intra_H_migration"],
    )
    record = {
        "reactant_graphs": _I084_FIXTURE.REACTANT_GRAPHS,
        "product_graphs": [_I084_QUINOID.PRODUCT_GRAPH],
        "template": "R7HJ_2;C_rad_out_OneDe/Cs;Cb_H_out",
    }
    reaction = _reaction_from_record(record)
    roots = _mapped_reaction_u2_roots(
        database.families["intra_H_migration"], reaction, record
    )
    assert roots
    assert roots[0]["recipe_label"] == "*2"
    assert roots[0]["family_forward_role"] == "product"
    assert roots[0]["persistent_neutral_divalent_carbon"] is False


def _root():
    return {
        "reactant_atom_index": 0,
        "recipe_label": "*1",
        "family_forward_role": "product",
    }


def _source(label="*1", role="product", persistent=True, source="training 629"):
    return [{
        "source": source,
        "mapped_roots": [{
            "recipe_label": label,
            "family_forward_role": role,
            "persistent_neutral_divalent_carbon": persistent,
        }],
    }]


def test_unsupported_persistent_carbene_refuses_both_directions():
    published, refusal = apply_persistent_carbene_policy(
        _pair(), _root(), _source(persistent=False)
    )
    assert published == []
    assert refusal["policy_version"] == PERSISTENT_CARBENE_POLICY_VERSION
    assert refusal["record_ids"] == ["evt_forward", "evt_reverse"]
    assert refusal["reason"] == "no selected rate contributor supports persistent u2 in the mapped role"
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
        compiler_module, "_mapped_reaction_u2_roots",
        lambda family, reaction, record: [_root()],
    )
    monkeypatch.setattr(
        compiler_module, "_training_source_domain",
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
        compiler_module, "_mapped_reaction_u2_roots",
        lambda family, reaction, record: [_root()],
    )
    monkeypatch.setattr(
        compiler_module, "_training_source_domain",
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
    published, refusal = apply_persistent_carbene_policy(
        _pair(), _root(), _source()
    )
    assert refusal is None
    assert len(published) == 2
    assert {record["reverse_of"] for record in published} == {
        "evt_forward", "evt_reverse"
    }
    annotations = [record["rate_source"]["applicability"] for record in published]
    assert all(item["disposition"] == "retained-supported-transfer" for item in annotations)
    assert all(item["caveat"] == "same reacting role; donor/acceptor environment may differ" for item in annotations)


def test_unserialized_source_domain_is_retained_as_unresolved():
    published, refusal = apply_persistent_carbene_policy(
        _pair(), _root(), None
    )
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
        lambda family, reaction: [{
            "recipe_label": "*1",
            "family_forward_role": (
                "reactant" if reaction.reactants == ["reverse"] else "product"
            ),
            "persistent_neutral_divalent_carbon": True,
        }],
    )
    domain = compiler_module._training_source_domain(Family(), object(), {})
    assert domain[0]["source"]["stored_reverse"] is True
    assert domain[0]["mapped_roots"][0]["family_forward_role"] == "reactant"


def test_real_training_629_forward_and_reverse_bind_the_same_supported_source():
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
    dispositions = []
    recipe_labels = []
    roles = []
    for reverse in (False, True):
        reaction = TemplateReaction(
            reactants=copy.deepcopy(base.products if reverse else base.reactants),
            products=copy.deepcopy(base.reactants if reverse else base.products),
            kinetics=copy.deepcopy(entry.data),
            family=family.label,
            template=base.template,
            degeneracy=base.degeneracy,
        )
        reaction.kinetics.comment = (
            f"Matched reaction 629 {entry.label} in training"
        )
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
        for side in ("reactant_graphs", "product_graphs"):
            record[side] = [
                graph.replace("*1", "*swap").replace("*3", "*1").replace(
                    "*swap", "*3"
                )
                for graph in record[side]
            ]
        roots = compiler_module._mapped_reaction_u2_roots(
            family, reaction, record
        )
        assert len(roots) == 1
        stored_atoms = [
            atom
            for graph in record["reactant_graphs"]
            for atom in Molecule().from_adjacency_list(graph).atoms
        ]
        assert stored_atoms[roots[0]["reactant_atom_index"]].label != roots[0][
            "recipe_label"
        ]
        decision = classify_persistent_carbene(
            record,
            roots[0],
            compiler_module._training_source_domain(family, reaction, {}),
        )
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
            assert wrong_decision["disposition"] == (
                "refused-structural-inconsistency"
            )
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
            return False, ("Rate Rules", {
                "training": [(Rule(), Entry(), 0.75)],
                "rules": [(Rule(), 0.25)],
            })

    monkeypatch.setattr(
        compiler_module,
        "_mapped_reaction_u2_roots",
        lambda family, reaction: [{
            "recipe_label": "*1",
            "family_forward_role": "product",
            "persistent_neutral_divalent_carbon": True,
        }],
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


def test_nine_form_delocalised_witness_stays_enabled_and_identity_bound():
    witness = json.loads(
        (FIXTURES / "i078_resonance_witness.json").read_text()
    )
    decision = classify_persistent_carbene(
        witness,
        {"reactant_atom_index": 74, "recipe_label": "*3", "family_forward_role": "product"},
        _source(persistent=False),
    )
    assert witness["event_id"] == "evt_453febdaa6658ad2fa22f6375f21b3c46d0eb8ac92c534f69fb7dfafb72665b5"
    assert decision["disposition"] == "not-applicable"
    assert decision["mapped_root"]["product_rewrite_verified"] is True
    assert decision["mapped_root"]["resonance_form_count"] == 9
    assert decision["mapped_root"]["persistent_neutral_divalent_carbon"] is False


def test_wrong_recipe_root_is_retained_only_as_unresolved():
    wrong_root = {
        "reactant_atom_index": 1,
        "recipe_label": "*1",
        "family_forward_role": "product",
    }
    decision = classify_persistent_carbene(
        _pair()[0], wrong_root, _source()
    )
    assert decision["disposition"] == "retained-unresolved-applicability"
    assert decision["mapped_root"]["mapping_verified"] is False
    assert "not touched" in decision["reason"]


@pytest.mark.parametrize(
    "family_candidates, preferred, reason",
    [
        ({1, 2}, {1, 2}, "ambiguous"),
        ({1}, {2}, "does not touch"),
        (set(), set(), "ambiguous"),
    ],
)
def test_ambiguous_or_unbound_family_root_becomes_unresolved_evidence(
    family_candidates, preferred, reason
):
    index, mapping_error = compiler_module._select_family_root_candidate(
        family_candidates, preferred, "*1"
    )
    assert index == min(family_candidates or preferred or {-1})
    assert reason in mapping_error
    decision = classify_persistent_carbene(
        _pair()[0],
        {
            **_root(),
            "reactant_atom_index": index,
            "mapping_verified": False,
            "mapping_error": mapping_error,
        },
        _source(),
    )
    assert decision["disposition"] == "retained-unresolved-applicability"


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
    assert decision["disposition"] == "retained-unresolved-applicability"
    assert "family relabeling failed" in decision["reason"]


def test_ordinary_radical_channel_is_byte_unchanged():
    pair = _pair()
    before = json.dumps(pair, sort_keys=True)
    ordinary_root = {
        "reactant_atom_index": 4,
        "recipe_label": "*2",
        "family_forward_role": "product",
    }
    published, refusal = apply_persistent_carbene_policy(
        pair, ordinary_root, _source()
    )
    assert refusal is None
    assert json.dumps(published, sort_keys=True) == before


def test_retained_pair_remains_reciprocal():
    published, _ = apply_persistent_carbene_policy(_pair(), _root(), _source())
    assert len(published) == 2
    by_id = {record["event_id"]: record for record in published}
    assert all(by_id[record["reverse_of"]]["reverse_of"] == record["event_id"] for record in published)


def _transition_pair(prefix, junction_kind=None, rate=1.0):
    forward, reverse = _pair()
    forward["event_id"], reverse["event_id"] = f"{prefix}_f", f"{prefix}_r"
    forward["reverse_of"], reverse["reverse_of"] = reverse["event_id"], forward["event_id"]
    for record in (forward, reverse):
        record["k_table"] = {"T": [1000.0], "k": [rate]}
        record["junction_ops"] = None
    if junction_kind:
        label = junction_kind.rsplit("_", 1)[-1]
        forward["junction_ops"] = [{
            "action": "create", "junction_kind": junction_kind,
            "attacked_atom_label": label,
            "chosen_resonance_localisation": f"ortho:{label}",
            "reflection": None if label == "S7" else {"S6": "S7", "S7": "S6"},
            "reverse_event_handle": reverse["event_id"],
        }]
        reverse["junction_ops"] = [{
            **forward["junction_ops"][0], "action": "dissociate",
            "reverse_event_handle": forward["event_id"],
        }]
    return [forward, reverse]


def test_s6_s7_split_rates_and_reverse_histories_are_not_duplicates():
    records = _transition_pair("s6", "J_ortho_S6", 0.5) + _transition_pair(
        "s7", "J_ortho_S7", 0.5
    )
    assert sum(record["k_table"]["k"][0] for record in records[::2]) == pytest.approx(1.0)
    assert duplicate_transition_groups(records) == []
    by_id = {record["event_id"]: record for record in records}
    assert all(by_id[record["reverse_of"]]["reverse_of"] == record["event_id"] for record in records)


def test_synthetic_duplicated_reciprocal_pair_is_detected():
    records = _transition_pair("first", rate=1.0) + _transition_pair(
        "duplicate", rate=2.0
    )
    assert duplicate_transition_groups(records) == [
        [("duplicate_f", "duplicate_r"), ("first_f", "first_r")]
    ]
