"""Slow real-RMG acceptance tests for the PS event-set compiler."""

import copy
import hashlib
import json
import multiprocessing
import os
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import pytest

from rmgpy.data.rmg import RMGDatabase
from rmgpy.kmc.compiler import (
    PS_FAMILY_CANDIDATES,
    PS_FAMILY_FILTER_REASON,
    PS_PROXY_UNITS,
    _canonical_adjacency,
    apply_record,
    canonical_json_bytes,
    ceiling_temperature,
    ps_proxy_set,
    short_ps_molecule_catalogue,
    validate_artifact,
)
from rmgpy.molecule.molecule import Molecule


REPO_ROOT = Path(__file__).resolve().parents[3]
DB_PATH = os.environ.get("RMG_DATABASE_PATH", str(REPO_ROOT.parent / "RMG-database"))
CROSS_PROCESS_TIMEOUT_SECONDS = 14 * 60
PRE_FIX_UNCHANGED_RECORD_COUNT = 312
PRE_FIX_UNCHANGED_RECORDS_SHA256 = (
    "7f3616d39972ce530fc5b2b1748a36f20a4452ac77bd71ba0661b7a5316b2a7d"
)
ORTHO_NODES = (
    ("P1", "C", 0, 1),
    ("P2", "C", 0, 3),
    ("P3", "C", 0, 0),
    ("P4", "C", 0, 1),
    ("P5", "C", 0, 1),
    ("P6", "C", 0, 1),
    ("P7", "C", 0, 1),
    ("P8", "C", 0, 1),
    ("P9", "C", 1, 2),
    ("S1", "C", 0, 2),
    ("S10", "C", 0, 1),
    ("S2", "C", 0, 2),
    ("S3", "C", 0, 3),
    ("S4", "C", 0, 0),
    ("S5", "C", 0, 1),
    ("S6", "C", 0, 1),
    ("S7", "C", 1, 1),
    ("S8", "C", 0, 1),
    ("S9", "C", 0, 1),
)
ORTHO_EDGES = (
    ("P1", "P2", 1.0),
    ("P1", "P3", 1.0),
    ("P1", "P9", 1.0),
    ("P3", "P4", 1.5),
    ("P3", "P5", 1.5),
    ("P4", "P6", 1.5),
    ("P5", "P8", 1.5),
    ("P6", "P7", 1.5),
    ("P7", "P8", 1.5),
    ("S1", "S2", 1.0),
    ("S1", "S3", 1.0),
    ("S10", "S7", 1.0),
    ("S10", "S9", 2.0),
    ("S2", "S5", 1.0),
    ("S4", "S5", 2.0),
    ("S4", "S6", 1.0),
    ("S4", "S7", 1.0),
    ("S6", "S8", 2.0),
    ("S8", "S9", 1.0),
)
pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(
        os.environ.get("RMG_KMC_SLOW") != "1", reason="set RMG_KMC_SLOW=1"
    ),
]


def _archive_molecule(prefixes, *, product=False, reflected=False):
    """Build one exact explicit-H graph from the archived mapped tuples."""
    selected = [node for node in ORTHO_NODES if node[0][0] in prefixes]
    reflection = {"S6": "S7", "S7": "S6", "S8": "S10", "S10": "S8"} if reflected else {}
    selected = [
        (reflection.get(label, label), element, radical, hydrogens)
        for label, element, radical, hydrogens in selected
    ]
    labels = {label for label, *_ in selected}
    attacked_label = "S6" if reflected else "S7"
    nodes = [
        (
            label,
            element,
            0 if product and label in {"P9", attacked_label} else radical,
            hydrogens,
        )
        for label, element, radical, hydrogens in selected
    ]
    edges = [
        (reflection.get(first, first), reflection.get(second, second), order)
        for first, second, order in ORTHO_EDGES
        if reflection.get(first, first) in labels
        and reflection.get(second, second) in labels
    ]
    if product:
        edges.append(("P9", attacked_label, 1.0))

    indices = {label: index for index, (label, *_rest) in enumerate(nodes, 1)}
    neighbors = {label: [] for label in labels}
    order_code = {1.0: "S", 1.5: "B", 2.0: "D"}
    for first, second, order in edges:
        neighbors[first].append((indices[second], order_code[order]))
        neighbors[second].append((indices[first], order_code[order]))

    lines = []
    next_index = len(nodes) + 1
    hydrogen_lines = []
    for label, element, radical, hydrogen_count in nodes:
        atom_index = indices[label]
        atom_bonds = list(neighbors[label])
        for _ in range(hydrogen_count):
            atom_bonds.append((next_index, "S"))
            hydrogen_lines.append(f"{next_index} H u0 p0 c0 {{{atom_index},S}}")
            next_index += 1
        bonds = " ".join(
            f"{{{neighbor},{order}}}" for neighbor, order in sorted(atom_bonds)
        )
        lines.append(
            f"{atom_index} *{label} {element} u{radical} p0 c0 {bonds}".rstrip()
        )
    multiplicity = 1 + sum(node[2] for node in nodes)
    adjacency = "\n".join([f"multiplicity {multiplicity}", *lines, *hydrogen_lines])
    molecule = Molecule().from_adjacency_list(adjacency)
    by_label = {atom.label: atom for atom in molecule.atoms if atom.label}
    for atom in molecule.atoms:
        atom.label = ""
    return molecule, by_label


def _canonical_graphs(molecules):
    return sorted(_canonical_adjacency(molecule) for molecule in molecules)


def _junction_records(artifact, *, action=None):
    records = [
        record
        for record in artifact["records"]
        if record.get("junction_ops")
        and record["junction_ops"][0]["junction_kind"].startswith("J_ortho")
    ]
    if action is not None:
        records = [
            record
            for record in records
            if record["junction_ops"][0]["action"] == action
        ]
    return records


def _normalized_base_record(record):
    normalized = copy.deepcopy(record)
    for field in ("canonical_index", "event_id", "provenance", "reverse_of"):
        normalized.pop(field, None)
    for operation in normalized.get("junction_ops") or ():
        operation.pop("reverse_event_handle", None)
    return normalized


def _assert_base_record_invariance(records):
    normalized = sorted(
        canonical_json_bytes(_normalized_base_record(record)).decode("ascii")
        for record in records
    )
    assert len(normalized) == PRE_FIX_UNCHANGED_RECORD_COUNT
    assert hashlib.sha256(canonical_json_bytes(normalized)).hexdigest() == (
        PRE_FIX_UNCHANGED_RECORDS_SHA256
    )


def _family_universe():
    root = Path(DB_PATH) / "input/kinetics/families"
    return sorted(
        path.name for path in root.iterdir() if (path / "groups.py").is_file()
    )


def _molecules(adjacencies):
    return [Molecule().from_adjacency_list(adjacency) for adjacency in adjacencies]


def _same_molecule_multiset(first, second):
    if len(first) != len(second):
        return False
    remaining = list(second)
    for molecule in first:
        for index, candidate in enumerate(remaining):
            if molecule.is_isomorphic(candidate):
                del remaining[index]
                break
        else:
            return False
    return True


def _record_reproduces_product(record):
    try:
        actual = apply_record(record, _molecules(record["reactant_graphs"]))
    except ValueError:
        return False
    return _same_molecule_multiset(actual, _molecules(record["product_graphs"]))


def _template(reaction):
    return ";".join(
        getattr(item, "label", str(item)) for item in (reaction.template or [])
    )


def _reaction_key(reaction, proxy):
    return (
        proxy.site_type,
        reaction.family,
        _template(reaction),
    )


def _record_key(record):
    return (
        record["site_type"],
        record["family"],
        record["template"],
    )


def _generate_independent_oracle(connection, kinetics_database, proxies):
    """Generate a fresh public-pipeline oracle in a forked process."""
    try:
        degeneracies = defaultdict(float)
        keys = set()
        families_by_site = {}
        for proxy in proxies:
            families = tuple(proxy.metadata["family_candidates"])
            generated = (
                kinetics_database.generate_reactions_from_families(
                    [participant.copy(deep=True) for participant in proxy.reactants],
                    only_families=list(families),
                    resonance=True,
                )
                if families
                else []
            )
            firing = {reaction.family for reaction in generated}
            if firing != set(families):
                raise AssertionError(
                    f"non-firing family routed to {proxy.site_type}: "
                    f"{sorted(set(families) - firing)}"
                )
            families_by_site[proxy.site_type] = families
            for reaction in generated:
                key = _reaction_key(reaction, proxy)
                keys.add(key)
                degeneracies[key] += float(reaction.degeneracy)
        connection.send(
            {
                "degeneracies": dict(degeneracies),
                "keys": keys,
                "families_by_site": families_by_site,
            }
        )
    finally:
        connection.close()


@pytest.fixture(scope="module")
def rmg_database():
    database = RMGDatabase()
    database.load_kinetics(
        DB_PATH + "/input/kinetics",
        reaction_libraries=[],
        seed_mechanisms=None,
        kinetics_families=list(PS_FAMILY_CANDIDATES),
        kinetics_depositories=["training"],
    )
    database.load_thermo(
        DB_PATH + "/input/thermo",
        thermo_libraries=["primaryThermoLibrary"],
        depository=True,
    )
    return database


@pytest.fixture(scope="module")
def compilation(rmg_database, tmp_path_factory):
    proxies = ps_proxy_set(PS_PROXY_UNITS)
    context = multiprocessing.get_context("fork")
    receive_oracle, send_oracle = context.Pipe(duplex=False)
    oracle_process = context.Process(
        target=_generate_independent_oracle,
        args=(send_oracle, rmg_database.kinetics, ps_proxy_set(PS_PROXY_UNITS)),
    )
    oracle_process.start()
    send_oracle.close()
    deadline = time.monotonic() + CROSS_PROCESS_TIMEOUT_SECONDS
    fixture_script = Path(__file__).with_name("compile_event_set_fixture.py")
    seeded_runs = []
    for seed in ("0", "4242"):
        output = tmp_path_factory.mktemp(f"hashseed-{seed}")
        environment = os.environ.copy()
        environment.update(
            {
                "PYTHONHASHSEED": seed,
                "PYTHONPATH": str(REPO_ROOT),
                "MPLCONFIGDIR": str(output / "matplotlib"),
            }
        )
        stdout = (output / "compile.stdout.log").open("w")
        stderr = (output / "compile.stderr.log").open("w")
        process = subprocess.Popen(
            [sys.executable, str(fixture_script), DB_PATH, str(output)],
            cwd=REPO_ROOT,
            env=environment,
            stdout=stdout,
            stderr=stderr,
        )
        seeded_runs.append((seed, output, process, stdout, stderr))
    try:
        for seed, output, process, _, _ in seeded_runs:
            try:
                process.wait(timeout=max(0.0, deadline - time.monotonic()))
            except subprocess.TimeoutExpired:
                process.terminate()
                process.wait(timeout=10)
                raise TimeoutError(
                    f"PYTHONHASHSEED={seed} compilation exceeded 14 minutes; "
                    f"see {output}"
                )
        if not receive_oracle.poll(max(0.0, deadline - time.monotonic())):
            raise TimeoutError("independent RMG oracle exceeded the suite budget")
        oracle = receive_oracle.recv()
    finally:
        for _, _, process, stdout, stderr in seeded_runs:
            if process.poll() is None:
                process.terminate()
                process.wait(timeout=10)
            stdout.close()
            stderr.close()
        receive_oracle.close()
        oracle_process.join(timeout=5)
        if oracle_process.is_alive():
            oracle_process.terminate()
            oracle_process.join()
    assert oracle_process.exitcode == 0
    paths = []
    for seed, output, process, _, _ in seeded_runs:
        assert process.returncode == 0, (
            f"PYTHONHASHSEED={seed} compilation failed; "
            f"see {output / 'compile.stderr.log'}"
        )
        artifacts = sorted(output.glob("*.json"))
        assert len(artifacts) == 1
        paths.append(artifacts[0])
    artifact = json.loads(paths[0].read_bytes())
    active = artifact["families"]
    excluded = artifact["excluded_families"]
    return (
        None,
        artifact,
        proxies,
        None,
        active,
        excluded,
        paths[0],
        oracle,
        tuple(paths),
    )


def test_family_enumeration_covers_every_database_family(compilation, rmg_database):
    """Every database family is active or has a documented exclusion reason."""
    _, artifact, _, _, active, excluded, _, _, _ = compilation
    universe = set(_family_universe())
    assert set(rmg_database.kinetics.families) == set(PS_FAMILY_CANDIDATES)
    assert set(active) | set(excluded) == universe
    assert set(active) == set(artifact["families"])
    assert all(reason for reason in excluded.values())
    assert set(active) <= set(PS_FAMILY_CANDIDATES)
    assert artifact["inputs"]["family_filter"]["candidates"] == sorted(
        PS_FAMILY_CANDIDATES
    )
    assert all(
        reason == PS_FAMILY_FILTER_REASON
        for family, reason in excluded.items()
        if family not in PS_FAMILY_CANDIDATES
    )

    corrupted_active = active[:-1]
    assert set(corrupted_active) | set(excluded) != universe


def test_cross_process_hash_seed_determinism(compilation):
    """Independent interpreters must emit identical canonical artifact bytes."""
    *_, paths = compilation
    seed_zero, seed_4242 = paths
    zero_bytes = seed_zero.read_bytes()
    other_bytes = seed_4242.read_bytes()
    assert zero_bytes == other_bytes
    assert seed_zero.name == seed_4242.name
    assert seed_zero.stem == hashlib.sha256(zero_bytes).hexdigest()


def test_j_ortho_archive_agreement_exact_reverses_and_controls(compilation):
    """Both ortho localisations reproduce the archived graph and invert exactly."""
    _, artifact, _, _, _, _, _, _, _ = compilation
    forward_records = _junction_records(artifact, action="create")
    assert {
        record["junction_ops"][0]["attacked_atom_label"] for record in forward_records
    } == {"S6", "S7"}

    primary, _ = _archive_molecule({"P"})
    secondary, secondary_labels = _archive_molecule({"S"})
    expected_product, _ = _archive_molecule({"P", "S"}, product=True)
    expected_reflected_product, _ = _archive_molecule(
        {"P", "S"}, product=True, reflected=True
    )
    records_by_id = {record["event_id"]: record for record in artifact["records"]}

    for forward in forward_records:
        expected = (
            expected_product
            if forward["junction_ops"][0]["attacked_atom_label"] == "S7"
            else expected_reflected_product
        )
        actual_product = apply_record(forward, _molecules(forward["reactant_graphs"]))
        assert _canonical_graphs(actual_product) == _canonical_graphs([expected])

        reverse = records_by_id[forward["reverse_of"]]
        assert reverse["reverse_of"] == forward["event_id"]
        assert _canonical_graphs(actual_product) == _canonical_graphs(
            _molecules(reverse["reactant_graphs"])
        )
        restored = apply_record(reverse, _molecules(reverse["reactant_graphs"]))
        assert _canonical_graphs(restored) == _canonical_graphs([primary, secondary])
        assert forward["radical_delta"] == -2
        assert reverse["radical_delta"] == 2
        for record in (forward, reverse):
            assert record["element_delta"] == {"C": 0, "H": 0}
            assert record["formula_delta"] == {"C": 0, "H": 0}
            assert record["implicit_h_delta"] == 0
            assert record["mass_delta"] == 0

    s7 = next(
        record
        for record in forward_records
        if record["junction_ops"][0]["attacked_atom_label"] == "S7"
    )
    corrupted = copy.deepcopy(s7)
    reactants = _molecules(corrupted["reactant_graphs"])
    secondary_actual = next(
        molecule for molecule in reactants if len(molecule.atoms) == 23
    )
    ordered_atoms = [atom for molecule in reactants for atom in molecule.atoms]
    isomorphism = secondary.find_isomorphism(secondary_actual)[0]
    s8_actual = isomorphism[secondary_labels["*S8"]]
    s8_index = ordered_atoms.index(s8_actual)
    formed = next(op for op in corrupted["bond_ops"] if op["action"] == "form")
    attacked_index = next(
        index
        for index in formed["atoms"]
        if ordered_atoms[index].radical_electrons
        and secondary_actual.is_atom_in_cycle(ordered_atoms[index])
    )
    formed["atoms"][formed["atoms"].index(attacked_index)] = s8_index
    with pytest.raises(ValueError):
        apply_record(corrupted, reactants)

    s6 = next(
        record
        for record in forward_records
        if record["junction_ops"][0]["attacked_atom_label"] == "S6"
    )
    missing_s6_forward = copy.deepcopy(artifact)
    missing_s6_forward["records"] = [
        record
        for record in missing_s6_forward["records"]
        if record["event_id"] != s6["event_id"]
    ]
    with pytest.raises(ValueError, match="S6"):
        validate_artifact(missing_s6_forward)

    missing_s6_reverse = copy.deepcopy(artifact)
    missing_s6_reverse["records"] = [
        record
        for record in missing_s6_reverse["records"]
        if record["event_id"] != s6["reverse_of"]
    ]
    with pytest.raises(ValueError, match="S6"):
        validate_artifact(missing_s6_reverse)


def test_j_ortho_rate_counts_source_degeneracy_once(compilation):
    """The two localisations partition, rather than duplicate, the g=2 rate."""
    _, artifact, _, _, _, _, _, _, _ = compilation
    forward_records = _junction_records(artifact, action="create")
    assert len(forward_records) == 2
    temperatures = forward_records[0]["k_table"]["T"]
    expected = [3.53586e10 * temperature**-1.00291 for temperature in temperatures]
    summed = [
        sum(record["k_table"]["k"][index] for record in forward_records)
        for index in range(len(temperatures))
    ]
    assert summed == pytest.approx(expected, rel=1e-6)
    assert [2.0 * value for value in summed] != pytest.approx(expected, rel=1e-6)
    assert sum(record["degeneracy"] for record in forward_records) == 2
    reverse_records = _junction_records(artifact, action="dissociate")
    index_600 = temperatures.index(600.0)
    assert all(
        record["k_table"]["k"][index_600] == pytest.approx(8.801895488669665, rel=1e-11)
        for record in reverse_records
    )


def test_archived_junctions_preserve_every_other_pre_fix_record(compilation):
    """Every non-J_para record matches the pre-fix artifact canonical content."""
    _, artifact, _, _, _, _, _, _, _ = compilation
    base_records = [
        record
        for record in artifact["records"]
        if not (
            record.get("junction_ops")
            and record["junction_ops"][0].get("junction_kind") == "J_para"
        )
    ]
    _assert_base_record_invariance(base_records)

    corrupted = copy.deepcopy(base_records)
    candidate = next(record for record in corrupted if record["bond_ops"])
    candidate["bond_ops"] = candidate["bond_ops"][1:]
    with pytest.raises(AssertionError):
        _assert_base_record_invariance(corrupted)


def test_c1_every_rewrite_reproduces_rmg_product(compilation):
    """C1: every rewrite is graph-complete, canonical, and corruption-sensitive."""
    _, artifact, _, _, _, _, _, _, _ = compilation
    assert artifact["records"]
    for record in artifact["records"]:
        size = len(record["atom_map"])
        assert set(map(int, record["atom_map"])) == set(range(size))
        assert set(record["atom_map"].values()) == set(range(size))
        assert _record_reproduces_product(record), record["event_id"]

    candidate = next(record for record in artifact["records"] if record["bond_ops"])
    corrupted = copy.deepcopy(candidate)
    corrupted["bond_ops"] = corrupted["bond_ops"][1:]
    assert not _record_reproduces_product(corrupted)


def test_c2_per_site_degeneracy_matches_whole_molecule_pipeline(compilation):
    """C2: per-site sums equal the independent public whole-molecule oracle."""
    _, artifact, proxies, _, _, _, _, oracle, _ = compilation
    generated_j_para = {
        key: value
        for key, value in oracle["degeneracies"].items()
        if key[:2] == ("junction_radical+end_radical", "R_Recombination")
    }
    assert len(generated_j_para) == 1
    expected = {
        key: value
        for key, value in oracle["degeneracies"].items()
        if key not in generated_j_para
    }
    proxy_by_site = {proxy.site_type: proxy for proxy in proxies}
    actual = defaultdict(float)
    for record in artifact["records"]:
        # Both approved junction channels replace the L=3 proxy estimate with
        # their archived, labelled R-009 reactions.
        if record.get("inventory_class") != "R1:J_ring":
            actual[_record_key(record)] += record["raw_path_degeneracy"]
    assert dict(actual) == pytest.approx(expected)
    required_cases = {
        "pristine",
        "interior_radical",
        "doubly_featured",
        "end_radical",
        "junction_radical",
        "end_radical+end_radical",
    }
    assert required_cases <= set(proxy_by_site)
    identical_pair_records = [
        record
        for record in artifact["records"]
        if record["site_type"] == "end_radical+end_radical" and record["arity"] == 2
    ]
    assert identical_pair_records
    assert all(record["ssa_multiplier"] == 2.0 for record in identical_pair_records)
    assert all(
        record["reactant_multiplicities"] == [2] for record in identical_pair_records
    )
    assert {proxy.metadata["context_class"] for proxy in proxies} == {
        "pristine",
        "featured",
        "end-proximal",
        "junction",
    }
    assert all(proxy.metadata["proxy_units"] == PS_PROXY_UNITS for proxy in proxies)
    assert oracle["families_by_site"] == {
        proxy.site_type: tuple(proxy.metadata["family_candidates"]) for proxy in proxies
    }

    corrupted = dict(actual)
    candidate_key = next(iter(corrupted))
    corrupted[candidate_key] += 1.0
    assert corrupted != pytest.approx(expected)


def test_c3_bounded_l3_generation_has_no_missing_family_template(compilation):
    """C3: independent bounded L=3 generation has no uncompiled template."""
    _, artifact, _, _, _, _, _, oracle, _ = compilation
    compiled = {
        _record_key(record)
        for record in artifact["records"]
        if record["site_type"] != "J_ring"
    }
    generated_j_para = {
        key
        for key in oracle["keys"]
        if key[:2] == ("junction_radical+end_radical", "R_Recombination")
    }
    assert len(generated_j_para) == 1
    required = oracle["keys"] - generated_j_para
    assert required <= compiled, sorted(required - compiled)[:10]
    archived_j_para = [
        record
        for record in artifact["records"]
        if record.get("junction_ops")
        and record["junction_ops"][0].get("junction_kind") == "J_para"
        and record["junction_ops"][0].get("action") == "create"
    ]
    assert len(archived_j_para) == 1
    assert archived_j_para[0]["rate_source"]["rule_entry_index"] == 176

    corrupted = set(compiled)
    corrupted.remove(next(iter(required)))
    assert not required <= corrupted


def test_c7_c9_inventory_and_real_artifact(compilation):
    """Check C7/C9 inventory, real ceiling, and ring reverse records."""
    _, artifact, _, _, _, _, artifact_path, _, _ = compilation
    catalogue = short_ps_molecule_catalogue(artifact["inputs"]["span_radius"])
    assert catalogue == artifact["short_molecule_catalogue"]
    assert catalogue["terminates"]
    assert catalogue["size"] == len(catalogue["molecules"])
    assert artifact["irreversible_pairs"]
    assert artifact["ps_ceiling_pairs"]
    pair = artifact["ps_ceiling_pairs"][0]
    records_by_id = {record["event_id"]: record for record in artifact["records"]}
    propagation = records_by_id[pair["propagation_event_id"]]
    depropagation = records_by_id[pair["depropagation_event_id"]]
    recomputed = ceiling_temperature(
        propagation,
        depropagation,
        pair["monomer_concentration_mol_m3"],
    )
    assert recomputed == pair["temperature_K"]
    assert artifact["ps_ceiling_temperature_K"] == recomputed
    ring_records = [
        record
        for record in artifact["records"]
        if record["inventory_class"] == "R1:J_ring"
    ]
    para_records = [
        record
        for record in ring_records
        if record["junction_ops"][0]["junction_kind"] == "J_para"
    ]
    assert ring_records
    assert para_records
    assert all(record["reverse_of"] for record in ring_records)
    assert all(len(record["event_id"]) == 68 for record in artifact["records"])
    for record in para_records:
        operation = record["junction_ops"][0]
        assert operation["junction_kind"] == "J_para"
        assert operation["attacker_site_label"] == "P9"
        assert operation["attacked_atom_label"] == "S9"
        assert operation["formed_bond"] == {
            "label_pair": ["S9", "P9"],
            "order": "1.0",
        }
        assert operation["chosen_resonance_localisation"] == "para:S9"
        assert operation["reverse_event_handle"] == record["reverse_of"]
        assert records_by_id[record["reverse_of"]]["reverse_of"] == record["event_id"]
        inverse = operation["exact_inverse_rewrite"]
        assert inverse
        for rewrite in inverse:
            references = rewrite.get("atoms", [rewrite.get("atom")])
            assert all(isinstance(reference, str) for reference in references)

    tampered_values = {
        "junction_kind": "J_ortho_S7",
        "attacker_site_label": "P8",
        "attacked_atom_label": "S7",
        "formed_bond": {"label_pair": ["S7", "P9"], "order": "1.0"},
        "chosen_resonance_localisation": "ortho:S7",
        "reverse_event_handle": "evt_" + "f" * 64,
        "exact_inverse_rewrite": [],
    }
    candidate_id = para_records[0]["event_id"]
    for field, value in tampered_values.items():
        corrupted = copy.deepcopy(artifact)
        corrupted_by_id = {item["event_id"]: item for item in corrupted["records"]}
        corrupted_by_id[candidate_id]["junction_ops"][0][field] = value
        with pytest.raises(ValueError):
            validate_artifact(corrupted)

    corrupted = copy.deepcopy(artifact)
    unclassified = next(
        item for item in corrupted["records"] if item["inventory_class"] is None
    )
    para = next(
        item for item in corrupted["records"] if item["inventory_class"] == "R1:J_ring"
    )
    unclassified["inventory_class"], para["inventory_class"] = (
        para["inventory_class"],
        unclassified["inventory_class"],
    )
    with pytest.raises(ValueError, match="event_id"):
        validate_artifact(corrupted)
    assert artifact_path.stem == hashlib.sha256(artifact_path.read_bytes()).hexdigest()

    corrupted = copy.deepcopy(artifact)
    corrupted["short_molecule_catalogue"]["size"] += 1
    with pytest.raises(ValueError, match="catalogue size"):
        validate_artifact(corrupted)
    corrupted_depropagation = copy.deepcopy(depropagation)
    corrupted_depropagation["k_table"]["k"] = [
        1.0e300 for _ in corrupted_depropagation["k_table"]["k"]
    ]
    assert (
        ceiling_temperature(
            propagation,
            corrupted_depropagation,
            pair["monomer_concentration_mol_m3"],
        )
        is None
    )
