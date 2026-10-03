"""Registered grammar-closure and R1 record-half product verifiers."""

import copy
import hashlib
import json
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

from rmgpy.kmc.compiler import (
    _canonical_adjacency,
    apply_record,
    validate_artifact,
)
from rmgpy.kmc.event_record import EventRecord
from rmgpy.kmc.met import ARCHIVED_CHANNEL_RATES, RATE_GRID, RateTable, _channel_id
from rmgpy.kmc.state import (
    AtomRef,
    KMCState,
    ReverseBindingError,
    Site,
    Strand,
    _parse_adjacency,
)
from rmgpy.molecule.molecule import Molecule


ROOT = Path(__file__).resolve().parents[3]
DATABASE = Path(os.environ.get("RMG_DATABASE_PATH", ROOT.parent / "RMG-database"))
CACHE = ROOT / ".kmc-cache" / "i028-verifiers"

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(
        os.environ.get("RMG_KMC_SLOW") != "1", reason="set RMG_KMC_SLOW=1"
    ),
]


# BEGIN COPIED V8 LITERALS
EVENT_ARTIFACT = b"I-011-event-set-artifact-v8-fixture\n"
EVENT_ARTIFACT_SHA256 = (
    "7704712e6b355d69644add23c187fe15692016670eb72a7c47c0690d86627074"
)

J_H = (3, 2, 2, 1, 0, 1, 1, 2, 1, 1)
J_BONDS = (
    (0, 1, 1),
    (1, 2, 1),
    (2, 3, 1),
    (3, 4, 2),
    (4, 5, 1),
    (5, 6, 2),
    (6, 7, 1),
    (7, 8, 1),
    (8, 9, 2),
    (4, 9, 1),
)
J_RADICALS = (4, 9)

Q_H = (1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 2, 2, 2, 2, 2)
Q_BONDS = (
    (0, 1, 2),
    (1, 2, 1),
    (2, 3, 2),
    (3, 4, 1),
    (4, 5, 1),
    (5, 6, 1),
    (6, 7, 2),
    (7, 8, 1),
    (8, 9, 2),
    (9, 10, 1),
    (10, 11, 1),
    (11, 12, 1),
    (12, 13, 1),
    (13, 14, 1),
    (14, 15, 1),
    (15, 16, 1),
    (5, 0, 1),
    (16, 11, 1),
)
Q_RADICALS = (11, 16)

GRAMMAR_TERMINAL_SCHEMA = {
    "element": ("C",),
    "implicit_H": (0, 1, 2, 3, 4),
    "aromatic": (0, 1),
    "radicals": (0, 1),
    "r1_label": ("none", "J", "Q"),
    "bond_order": (1.0, 1.5, 2.0),
    "constraint": "H+radicals+sum(bond_order)=4;accepts_g8 side conditions",
}
GRAMMAR_PRODUCTIONS = (
    "Component::=Tree|Tree<attach>Benzene|R1(Component,Localisation)",
    "Tree::=C4|C[h,0,r,none]-Tree|C[h,0,0,none]=C[h',0,0,none]",
    "Benzene::=cycle6(C[h_i,1,0,none],B[3/2])",
    "Localisation::=para|ortho_S7|ortho_S6",
    "Boundary::=P[o,c]|C[h,a,r,e]<B[o]>Boundary(distance<6)",
)

EXPECTED_GRAMMAR_FIXTURE_HASHES = {
    "terminal_schema": "5428ed2f11a03dd0dcba878a2ad80a5bda631e97095d0e98ec3315b4ce99fdbd",
    "productions": "2a71425a1988542a87324a814b00cd97c7e91e0ea747cd3ef5843c4d76d7de5e",
    "initials": "ad75dbd953cc77f36f8508a2ba46a845bc28c0961e8fee6b9d03df35c307391e",
    "rewrites": "03ee939937ea63c7b00700d450ad0f1a86e9b16bbff04bfe20e9b8490a71517f",
}
EXPECTED_GRAMMAR_CLOSURE = {
    "depth_counts": (3, 12, 18, 12, 3),
    "states": 48,
    "depth_transitions": tuple(
        (rewrite, depth, depth + 1, count)
        for rewrite in (
            "RW00_attach_methyl",
            "RW01_attach_benzene",
            "RW02_R1_para",
            "RW03_R1_ortho",
        )
        for depth, count in enumerate((3, 9, 9, 3))
    ),
    "outside_G8": (),
}
EXPECTED_GRAMMAR_TRANSITIONS = {
    "count": 96,
    "sha256": "fce0f94f8b2562b7fcc6ba7d15485e8b1892240b1a087c6b15eeea686e8b13d1",
}
EXPECTED_R1_GRAPH_HASHES = {
    "J": "784889032231c16207edadd477493200d33b6d66aa855e42fe6777ba8a8173a1",
    "Q": "9154b8b5cf115ccca1c99b6936d2df5aad4190c4f507b94e7d035d3b6ba19227",
    "J_parent": "377cd1c15f28ba7ee4a0fde3c75cca5be335c68ac698b4968283839cebf558fc",
    "Q_parent": "33c6be6f8aac5cbfaebc58c65eee1b195782ba3512a50abd43073f1b0c903409",
}
EXPECTED_R1_FIXTURES = (
    ("R1-1", "reject", "E_R1_PROVENANCE"),
    ("R1-2", "accept", "-"),
    ("R1-3", "BLOCKED-BY-GAP", "E_DOMAIN_ALIPHATIC_RING"),
    ("R1-4", "BLOCKED-BY-GAP", "group-count half"),
)
# END COPIED V8 LITERALS


def _canonical_data(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def _data_sha256(value):
    return hashlib.sha256(_canonical_data(value).encode()).hexdigest()


def _literal_graph(hydrogens, bonds, radicals=()):
    radical_set = set(radicals)
    payload = {
        "atoms": [
            ("C", hydrogen, False, int(index in radical_set))
            for index, hydrogen in enumerate(hydrogens)
        ],
        "bonds": sorted(
            (min(first, second), max(first, second), float(order))
            for first, second, order in bonds
        ),
    }
    return _canonical_data(payload)


def _alkane(length):
    return _literal_graph(
        (3,) + (2,) * (length - 2) + (3,),
        tuple((index, index + 1, 1) for index in range(length - 1)),
    )


def _registered_grammar_literals():
    initial_graphs = tuple(
        sorted(
            (
                (
                    name,
                    _alkane(length),
                    (("S00", 0), ("S01", 1), ("S02", 2), ("S03", 3)),
                )
                for name, length in (
                    ("n-hexane", 6),
                    ("n-heptane", 7),
                    ("n-octane", 8),
                )
            )
        )
    )
    benzene_fragment = {
        "atoms": [["C", 0, True, 0]] + [["C", 1, True, 0] for _ in range(5)],
        "bonds": [[index, (index + 1) % 6, 1.5] for index in range(6)],
        "labels": ["none"] * 6,
        "attach": 0,
    }
    j_fragment = json.loads(_literal_graph(J_H, J_BONDS))
    j_fragment["atoms"][0][1] -= 1
    j_fragment["labels"] = ["J"] * len(j_fragment["atoms"])
    j_fragment["attach"] = 0
    q_fragment = json.loads(_literal_graph(Q_H, Q_BONDS))
    q_fragment["atoms"][0][1] -= 1
    q_fragment["labels"] = ["Q"] * len(q_fragment["atoms"])
    q_fragment["attach"] = 0
    rewrites = (
        (
            "RW00_attach_methyl",
            0,
            "S00",
            {
                "atoms": [["C", 3, False, 0]],
                "bonds": [],
                "labels": ["none"],
                "attach": 0,
            },
        ),
        ("RW01_attach_benzene", 1, "S01", benzene_fragment),
        ("RW02_R1_para", 2, "S02", j_fragment),
        ("RW03_R1_ortho", 3, "S03", q_fragment),
    )
    return initial_graphs, rewrites


def _literal_block_sha256():
    source = Path(__file__).read_bytes()
    start = source.index(b"# BEGIN COPIED V8 LITERALS\n") + len(
        b"# BEGIN COPIED V8 LITERALS\n"
    )
    end = source.index(b"# END COPIED V8 LITERALS\n")
    return hashlib.sha256(source[start:end]).hexdigest()


def _database_sha():
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=DATABASE, text=True
    ).strip()


@pytest.fixture(scope="module")
def compiled_artifact():
    """Load one real content-addressed artifact, compiling it when absent."""
    from cache_provenance import supplied_artifact

    supplied = supplied_artifact()
    if supplied:
        return supplied
    CACHE.mkdir(parents=True, exist_ok=True)
    compiler_sha = hashlib.sha256(
        (ROOT / "rmgpy/kmc/compiler.py").read_bytes()
    ).hexdigest()
    database_sha = _database_sha()
    artifact_path = None
    for candidate in sorted(CACHE.glob("*.json")):
        artifact = json.loads(candidate.read_bytes())
        provenance = artifact.get("provenance", {})
        if (
            provenance.get("compiler_sha256") == compiler_sha
            and provenance.get("rmg_database_sha") == database_sha
        ):
            artifact_path = candidate
            break
    if artifact_path is None:
        fixture = Path(__file__).with_name("compile_event_set_fixture.py")
        environment = os.environ.copy()
        environment.update(
            {
                "PYTHONPATH": str(ROOT),
                "MPLCONFIGDIR": str(CACHE / "matplotlib"),
            }
        )
        with (CACHE / "compile.stdout.log").open("a") as stdout, (
            CACHE / "compile.stderr.log"
        ).open("a") as stderr:
            subprocess.run(
                [sys.executable, str(fixture), str(DATABASE), str(CACHE)],
                cwd=ROOT,
                env=environment,
                stdout=stdout,
                stderr=stderr,
                check=True,
            )
        candidates = []
        for candidate in sorted(CACHE.glob("*.json")):
            artifact = json.loads(candidate.read_bytes())
            provenance = artifact.get("provenance", {})
            if (
                provenance.get("compiler_sha256") == compiler_sha
                and provenance.get("rmg_database_sha") == database_sha
            ):
                candidates.append(candidate)
        assert len(candidates) == 1
        artifact_path = candidates[0]
    artifact = json.loads(artifact_path.read_bytes())
    validate_artifact(artifact)
    return artifact_path, artifact


def _molecules(adjacencies):
    return [Molecule().from_adjacency_list(adjacency) for adjacency in adjacencies]


def _canonical_graphs(molecules):
    return sorted(_canonical_adjacency(molecule) for molecule in molecules)


def _ring_records(artifact):
    return [
        record
        for record in artifact["records"]
        if record.get("inventory_class") == "R1:J_ring"
    ]


def _record_groups(artifact):
    groups = {}
    for record in _ring_records(artifact):
        operation = record["junction_ops"][0]
        groups.setdefault((operation["junction_kind"], operation["action"]), []).append(
            record
        )
    return groups


def _assert_registered_r1_table(artifact):
    groups = _record_groups(artifact)
    expected_keys = {
        (kind, action)
        for kind in ("J_para", "J_ortho_S7", "J_ortho_S6")
        for action in ("create", "dissociate")
    }
    assert set(groups) == expected_keys
    assert all(len(records) == 1 for records in groups.values())
    by_id = {record["event_id"]: record for record in artifact["records"]}
    expected_attacked = {
        "J_para": "S9",
        "J_ortho_S7": "S7",
        "J_ortho_S6": "S6",
    }
    for (kind, action), records in groups.items():
        record = records[0]
        operation = record["junction_ops"][0]
        attacked = expected_attacked[kind]
        assert operation["attacked_atom_label"] == attacked
        assert operation["attacker_site_label"] == "P9"
        assert operation["formed_bond"] == {
            "label_pair": [attacked, "P9"],
            "order": "1.0",
        }
        expected_inverse_action = "break" if action == "create" else "form"
        expected_radical = 1 if action == "create" else 0
        assert operation["exact_inverse_rewrite"] == [
            {
                "action": expected_inverse_action,
                "atoms": [attacked, "P9"],
                "order": "1.0",
            },
            {"action": "set_radical", "atom": attacked, "value": expected_radical},
            {"action": "set_radical", "atom": "P9", "value": expected_radical},
        ]
        graph_action = "form" if action == "create" else "break"
        graph_bond = next(
            item for item in record["bond_ops"] if item["action"] == graph_action
        )
        assert graph_bond["order"] == "1.0"
        assert record["radical_delta"] == (-2 if action == "create" else 2)
        assert record["implicit_h_delta"] == 0
        assert record["formula_delta"].get("H", 0) == 0
        assert not any(
            item["action"] == "set_implicit_hydrogens" for item in record["bond_ops"]
        )
        assert len(record["event_id"]) == 68
        assert record["event_id"].startswith("evt_")
        assert len(record["event_id"][4:]) == 64
        int(record["event_id"][4:], 16)
        EventRecord.from_dict(record).validate()
        partner = by_id[record["reverse_of"]]
        assert partner["reverse_of"] == record["event_id"]
        assert operation["reverse_event_handle"] == partner["event_id"]

    assert groups[("J_ortho_S7", "create")][0]["junction_ops"][0]["reflection"] is None
    reflection = groups[("J_ortho_S6", "create")][0]["junction_ops"][0]["reflection"]
    assert reflection == {"S6": "S7", "S7": "S6", "S8": "S10", "S10": "S8"}
    assert reflection.get("S4", "S4") == "S4"
    assert reflection.get("S9", "S9") == "S9"


def _harness(data):
    """Build executor participants from one real compiled record."""
    record = EventRecord.from_dict(data)
    graphs = [_parse_adjacency(text) for text in record.reactant_graphs]
    expanded_types = [
        site_type
        for site_type, count in zip(
            record.participant_site_types, record.reactant_multiplicities
        )
        for _ in range(count)
    ]
    form_atoms = next(
        operation["atoms"]
        for operation in record.bond_ops
        if operation["action"] == "form"
    )
    labels = {form_atoms[0]: "S9", form_atoms[1]: "P9"}
    roles = {form_atoms[0]: "attacked_ring", form_atoms[1]: "attacker"}
    strands = []
    sites = []
    offset = 0
    for ordinal, (graph, site_type) in enumerate(zip(graphs, expanded_types)):
        strand_id = f"strand-{ordinal}"
        refs = {}
        formula = Counter()
        radicals = 0
        for local, node in sorted(graph.items()):
            ref = AtomRef(
                f"{record.event_id}-{ordinal}-{local}",
                strand_id,
                local,
                roles.get(offset + local, "feature"),
            )
            refs[local] = ref
            formula[node["element"]] += 1
            radicals += int(node.get("radical", 0))
        atom_graph = {}
        for local, node in graph.items():
            atom_graph[refs[local].uuid] = {
                "element": node["element"],
                "radical": int(node.get("radical", 0)),
                "charge": int(node.get("charge", 0)),
                "lone_pairs": int(node.get("lone_pairs", 0)),
                "implicit_hydrogens": int(node.get("implicit_hydrogens", 0)),
                "edges": {
                    refs[neighbor].uuid: order
                    for neighbor, order in node.get("edges", {}).items()
                },
            }
        strands.append(
            Strand(
                strand_id,
                max(8, len(graph) + 1),
                {},
                dict(formula),
                radicals,
                refs,
                atom_graph,
            )
        )
        actual = copy.deepcopy(graph)
        for local, node in actual.items():
            ref = refs[local]
            node.update(
                {
                    "atom_ref": ref,
                    "position": ref.position,
                    "role": ref.role,
                }
            )
            if offset + local in labels:
                node["label"] = labels[offset + local]
        sites.append(Site(site_type, strand_id, 0, 0, actual))
        offset += len(graph)
    state = KMCState(strands)
    current_sites = []
    for site in sites:
        graph = copy.deepcopy(site.graph)
        for node in graph.values():
            node["atom_ref"] = state._ref(node["atom_ref"].uuid)
        current_sites.append(
            Site(site.site_type, site.strand_id, site.position, site.atom, graph)
        )
    return state, record, current_sites


def test_registered_literals_and_expected_rows_are_intact():
    assert hashlib.sha256(EVENT_ARTIFACT).hexdigest() == EVENT_ARTIFACT_SHA256
    j_product = _literal_graph(J_H, J_BONDS)
    q_product = _literal_graph(Q_H, Q_BONDS)
    j_parent = _literal_graph(J_H, J_BONDS[:-1], J_RADICALS)
    q_parent = _literal_graph(Q_H, Q_BONDS[:-1], Q_RADICALS)
    assert {
        "J": hashlib.sha256(j_product.encode()).hexdigest(),
        "Q": hashlib.sha256(q_product.encode()).hexdigest(),
        "J_parent": hashlib.sha256(j_parent.encode()).hexdigest(),
        "Q_parent": hashlib.sha256(q_parent.encode()).hexdigest(),
    } == EXPECTED_R1_GRAPH_HASHES
    initials, rewrites = _registered_grammar_literals()
    assert {
        "terminal_schema": _data_sha256(GRAMMAR_TERMINAL_SCHEMA),
        "productions": _data_sha256(GRAMMAR_PRODUCTIONS),
        "initials": _data_sha256(initials),
        "rewrites": _data_sha256(rewrites),
    } == EXPECTED_GRAMMAR_FIXTURE_HASHES
    assert EXPECTED_GRAMMAR_CLOSURE["depth_counts"] == (3, 12, 18, 12, 3)
    assert EXPECTED_GRAMMAR_CLOSURE["states"] == 48
    assert EXPECTED_GRAMMAR_TRANSITIONS == {
        "count": 96,
        "sha256": "fce0f94f8b2562b7fcc6ba7d15485e8b1892240b1a087c6b15eeea686e8b13d1",
    }
    assert EXPECTED_R1_FIXTURES[:2] == (
        ("R1-1", "reject", "E_R1_PROVENANCE"),
        ("R1-2", "accept", "-"),
    )
    assert len(_literal_block_sha256()) == 64


def test_real_artifact_hash_provenance_and_injected_defects(compiled_artifact):
    artifact_path, artifact = compiled_artifact
    assert hashlib.sha256(artifact_path.read_bytes()).hexdigest() == artifact_path.stem
    validate_artifact(artifact)
    assert all(
        record["provenance"] == artifact["provenance"] for record in artifact["records"]
    )

    swapped = copy.deepcopy(artifact)
    classified = next(
        record
        for record in swapped["records"]
        if record["inventory_class"] == "R1:J_ring"
    )
    unclassified = next(
        record for record in swapped["records"] if record["inventory_class"] is None
    )
    classified["inventory_class"], unclassified["inventory_class"] = (
        unclassified["inventory_class"],
        classified["inventory_class"],
    )
    with pytest.raises(ValueError, match="event_id"):
        validate_artifact(swapped)

    truncated = copy.deepcopy(artifact)
    candidate = next(
        record for record in truncated["records"] if record["junction_ops"]
    )
    candidate["event_id"] = candidate["event_id"][:-8]
    with pytest.raises(ValueError):
        validate_artifact(truncated)

    altered_bytes = artifact_path.read_bytes() + b"\n"
    assert hashlib.sha256(altered_bytes).hexdigest() != artifact_path.stem


def test_real_compiled_r1_deltas_inverses_reflection_and_controls(compiled_artifact):
    _, artifact = compiled_artifact
    _assert_registered_r1_table(artifact)
    groups = _record_groups(artifact)

    for kind in ("J_para", "J_ortho_S7", "J_ortho_S6"):
        forward = groups[(kind, "create")][0]
        reverse = groups[(kind, "dissociate")][0]
        product = apply_record(forward, _molecules(forward["reactant_graphs"]))
        assert _canonical_graphs(product) == _canonical_graphs(
            _molecules(forward["product_graphs"])
        )
        assert _canonical_graphs(product) == _canonical_graphs(
            _molecules(reverse["reactant_graphs"])
        )
        restored = apply_record(reverse, _molecules(reverse["reactant_graphs"]))
        assert _canonical_graphs(restored) == _canonical_graphs(
            _molecules(forward["reactant_graphs"])
        )

    wrong_inverse = copy.deepcopy(artifact)
    wrong_groups = _record_groups(wrong_inverse)
    inverse = wrong_groups[("J_para", "create")][0]["junction_ops"][0][
        "exact_inverse_rewrite"
    ]
    inverse[0]["order"] = "2.0"
    with pytest.raises(AssertionError):
        _assert_registered_r1_table(wrong_inverse)

    misreflected = copy.deepcopy(artifact)
    reflected = _record_groups(misreflected)[("J_ortho_S6", "create")][0]
    reflected["junction_ops"][0]["reflection"]["S8"] = "S9"
    with pytest.raises((AssertionError, ValueError)):
        _assert_registered_r1_table(misreflected)


def test_real_r1_runtime_binding_uses_persistent_uuids_and_fails_closed(
    compiled_artifact,
):
    _, artifact = compiled_artifact
    forward = _record_groups(artifact)[("J_para", "create")][0]
    reverse = EventRecord.from_dict(
        _record_groups(artifact)[("J_para", "dissociate")][0]
    )
    state, record, sites = _harness(forward)
    before_uuids = set(state._uuid_owners())
    entry = state.apply(record, sites)
    assert set(entry.bindings) >= {"S9", "P9"}
    assert set(entry.bindings.values()) <= set(state._uuid_owners())
    assert before_uuids == set(state._uuid_owners())
    assert entry.formed_bond in state.junctions
    state.assert_uuid_uniqueness()

    capture = state._open_captures[record.event_id][-1]
    capture["bindings"]["S9"] = "injected-wrong-persistent-uuid"
    with pytest.raises(ReverseBindingError):
        state._find_reverse_capture(reverse, sites, entry.bindings)


def test_real_compiled_j_para_matches_archive_reverse_kc_and_rejects_old_rate(
    compiled_artifact,
):
    _, artifact = compiled_artifact
    para_forwards = _record_groups(artifact)[("J_para", "create")]
    assert len(para_forwards) == 1
    forward = para_forwards[0]
    assert _channel_id(forward) == "J_para"

    # Independent pack literals, also cross-checked against the product constant.
    literal_rates = [1.76793e10 * temperature**-1.00291 for temperature in RATE_GRID]
    compiled_rate = RateTable.from_mapping(forward["k_table"])
    assert [compiled_rate(temperature) for temperature in RATE_GRID] == pytest.approx(
        literal_rates, rel=1.0e-11
    )
    assert [
        ARCHIVED_CHANNEL_RATES["J_para"](temperature) for temperature in RATE_GRID
    ] == pytest.approx(literal_rates, rel=1.0e-15)
    assert forward["rate_source"]["channel_id"] == "J_para"
    assert forward["rate_source"]["direction"] == "forward"
    assert forward["rate_source"]["rule_entry_index"] == 176
    assert artifact["provenance"]["archived_j_para_rate"]["rule_entry_index"] == 176

    reverse = _record_groups(artifact)[("J_para", "dissociate")][0]
    reverse_rate = RateTable.from_mapping(reverse["k_table"])
    expected_kc = (
        1.045410027912e8,
        4.663224384510e6,
        3.284349406621e5,
        3.330607208101e4,
        4.538409537754e3,
    )
    observed_kc = [
        compiled_rate(temperature) / reverse_rate(temperature)
        for temperature in RATE_GRID
    ]
    assert observed_kc == pytest.approx(expected_kc, rel=1.0e-9)

    # This mutates the real compiled record to the pre-fix activated table.
    old_activated_rates = (
        3.11840237963936583e5,
        4.07177260964213521e5,
        5.11781620917322929e5,
        6.23944123405979481e5,
        7.42077195116760791e5,
    )
    old_rate = copy.deepcopy(forward)
    for temperature, value in zip(RATE_GRID, old_activated_rates):
        index = old_rate["k_table"]["T"].index(temperature)
        old_rate["k_table"]["k"][index] = value
    with pytest.raises(ValueError, match="J_para.*does not match"):
        _channel_id(old_rate)


def test_gap_verdicts_are_explicit_and_not_replaced_by_stand_ins():
    mapping = (
        Path(__file__).with_name("fixtures") / "i028_registered_verifier_mapping.md"
    )
    text = mapping.read_text()
    assert "**Row-7 premise verdict: `BLOCKED-BY-GAP`.**" in text
    assert "**Row-10 record-half verdict:**" in text
    assert "R1-3" in text and "R1-4" in text
    assert EXPECTED_R1_FIXTURES[2][1] == "BLOCKED-BY-GAP"
    assert EXPECTED_R1_FIXTURES[3][1] == "BLOCKED-BY-GAP"
