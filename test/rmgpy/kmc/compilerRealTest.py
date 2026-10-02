"""Slow real-RMG acceptance tests for the PS event-set compiler."""

import copy
import base64
import zlib
from dataclasses import replace
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
    EventSetCompiler,
    PS_FAMILY_CANDIDATES,
    PS_FAMILY_FILTER_REASON,
    PS_PROXY_UNITS,
    _canonical_adjacency,
    compiler_source_hash,
    apply_record,
    canonical_json_bytes,
    ceiling_temperature,
    ps_proxy_set,
    short_ps_molecule_catalogue,
    validate_artifact,
)
from rmgpy.molecule.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from completeness_oracle import (
    independent_c3_oracle,
    oracle_cache_key,
    reaction_graph_key,
)
from cache_provenance import generator_code_unchanged


REPO_ROOT = Path(__file__).resolve().parents[3]
DB_PATH = os.environ.get("RMG_DATABASE_PATH", str(REPO_ROOT.parent / "RMG-database"))
CROSS_PROCESS_TIMEOUT_SECONDS = 8 * 60 * 60
ARCHIVED_PACK_EXEMPTION_REASON = (
    "owner-approved transport and cage pack rates and Kc (R-009 v12)"
)
ARCHIVED_PACK_EVENT_IDS = (
    "evt_16eef4f25e21b91eaf71e9e99b1442421d71c971962030170571e48edb781e90",
    "evt_230f8f7e0de2176cc15eb832cdc83ad350eebf271ebb1f08541c6d9830d945c7",
    "evt_3700408a3c4cfbf91f2d9ae151a6c513349c28c1ced68257c7ce9836f3dd2a04",
    "evt_a511ddbb6b960c85a8fc91105770e0210d740deedbd47898e801a04283469baf",
    "evt_a86efb350529124c22adf46d1b420f7fb22d87eb73e13cd1cc915b5688de649c",
    "evt_a88efbdc0f390f31d6c4b1fb9e50e52ebc45ca262c2009e1f317a2854a578b9f",
)
PARA_PACK_KC = (
    1.045410027912e8,
    4.663224384510e6,
    3.284349406621e5,
    3.330607208101e4,
    4.538409537754e3,
)
ORTHO_PACK_KC = (
    6571782.521596457,
    377346.21958396066,
    32995.400975202465,
    4036.333967679853,
    648.2373622040287,
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


def _pre_change_artifact():
    path = Path(__file__).with_name("fixtures") / "i035_pre_change_artifact.zlib.b64"
    return json.loads(zlib.decompress(base64.b64decode(path.read_bytes())))


def _assert_base_record_invariance(artifact):
    baseline = _pre_change_artifact()
    paired_graphs = {
        reaction_graph_key(
            record["family"],
            _molecules(record["reactant_graphs"]),
            _molecules(record["product_graphs"]),
        )
        for record in artifact["records"]
        if record["reverse_of"] and record["inventory_class"] != "R1:J_ring"
    }
    for record in baseline["records"]:
        if record["inventory_class"] != "R1:J_ring":
            assert (
                reaction_graph_key(
                    record["family"],
                    _molecules(record["reactant_graphs"]),
                    _molecules(record["product_graphs"]),
                )
                in paired_graphs
            )
    baseline_records = [
        record
        for record in baseline["records"]
        if record["inventory_class"] == "R1:J_ring"
    ]
    actual = [
        record
        for record in artifact["records"]
        if record["inventory_class"] == "R1:J_ring"
    ]

    def normalized(records):
        return sorted(
            canonical_json_bytes(_normalized_base_record(record)) for record in records
        )

    assert normalized(actual) == normalized(baseline_records)


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
        cache = REPO_ROOT / ".kmc-cache/c3-oracle"
        key = oracle_cache_key(REPO_ROOT, DB_PATH, 1)
        destination = cache / (key + ".json")
        current = key.split("-")[0]
        if not destination.is_file():
            suffix = key[len(current) :] + ".json"
            for previous in sorted(cache.glob("*" + suffix)):
                origin = previous.name.split("-")[0]
                if generator_code_unchanged(REPO_ROOT, origin, current):
                    data = json.loads(previous.read_bytes())
                    data["cache_provenance"] = {
                        "independent_generation_cache": str(previous),
                        "generator_origin_commit": origin,
                        "validated_current_commit": current,
                        "database_commit": key.split("-")[1],
                        "identical_oracle_and_input_source_hash": key.split("-")[2],
                    }
                    destination.write_text(json.dumps(data, sort_keys=True))
                    break
        full = independent_c3_oracle(kinetics_database, REPO_ROOT, DB_PATH, 1, cache)
        connection.send(
            {
                "degeneracies": {
                    tuple(item["key"]): item["value"] for item in full["degeneracies"]
                },
                "keys": {tuple(key) for key in full["keys"]},
                "counts": full["counts"],
                "reactions": full["reactions"],
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
def compilation(rmg_database):
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
    cache_key = "-".join(
        [
            subprocess.check_output(
                ["git", "-C", str(path), "rev-parse", "HEAD"], text=True
            ).strip()
            for path in (REPO_ROOT, DB_PATH)
        ]
        + [compiler_source_hash()]
    )
    cache_root = (
        Path(os.environ.get("RMG_KMC_CACHE_ROOT", str(REPO_ROOT / ".kmc-cache")))
        / "hash-seeded-compiler"
        / cache_key
    )
    seeded_runs = []
    for seed in ("0", "4242"):
        output = cache_root / f"hashseed-{seed}"
        output.mkdir(parents=True, exist_ok=True)
        artifacts = list(output.glob("*.json"))
        if len(artifacts) == 1:
            payload = artifacts[0].read_bytes()
            assert artifacts[0].stem == hashlib.sha256(payload).hexdigest()
            cached = json.loads(payload)
            assert (
                cached["provenance"]["compiler_sources_sha256"]
                == compiler_source_hash()
            )
            seeded_runs.append((seed, output, None, None, None))
            continue
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
            if process is None:
                continue
            try:
                process.wait(timeout=max(0.0, deadline - time.monotonic()))
            except subprocess.TimeoutExpired:
                process.terminate()
                process.wait(timeout=10)
                raise TimeoutError(
                    f"PYTHONHASHSEED={seed} compilation exceeded {CROSS_PROCESS_TIMEOUT_SECONDS} seconds; "
                    f"see {output}"
                )
        if not receive_oracle.poll(max(0.0, deadline - time.monotonic())):
            raise TimeoutError(
                "independent full-molecule RMG oracle exceeded the suite budget"
            )
        oracle = receive_oracle.recv()
    finally:
        for _, _, process, stdout, stderr in seeded_runs:
            if process is not None and process.poll() is None:
                process.terminate()
                process.wait(timeout=10)
            if stdout is not None:
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
        assert process is None or process.returncode == 0, (
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
    """The unpaired-change exemption cannot conceal a changed archived record."""
    _, artifact, _, _, _, _, _, _, _ = compilation
    _assert_base_record_invariance(artifact)
    corrupted = copy.deepcopy(artifact)
    candidate = next(
        record
        for record in corrupted["records"]
        if record["inventory_class"] == "R1:J_ring" and record["bond_ops"]
    )
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
    generated_j_para = {_record_key(entry) for entry in artifact["excluded_channels"]}
    expected = {
        key: value
        for key, value in oracle["degeneracies"].items()
        if key not in generated_j_para
    }
    proxy_by_site = {proxy.site_type: proxy for proxy in proxies}

    def observed_degeneracies(compiled_artifact):
        totals = defaultdict(float)
        by_id = {record["event_id"]: record for record in compiled_artifact["records"]}
        selected_sources = set()
        for discovery in compiled_artifact["discovery"]:
            assert discovery["event_id"] in by_id
            target = by_id[discovery["event_id"]]
            pair = frozenset((target["event_id"], target["reverse_of"]))
            if pair not in selected_sources:
                assert target["raw_path_degeneracy"] == discovery["raw_path_degeneracy"]
                selected_sources.add(pair)
            if discovery["proxy_units"] == 5:
                totals[_record_key(discovery)] += discovery["raw_path_degeneracy"]
        return totals

    actual = observed_degeneracies(artifact)
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
    assert all(
        record["ssa_multiplier"] in {1.0, 2.0} for record in identical_pair_records
    )
    assert all(
        sum(record["reactant_multiplicities"]) == 2 for record in identical_pair_records
    )
    assert {proxy.metadata["context_class"] for proxy in proxies} == {
        "pristine",
        "featured",
        "end-proximal",
        "junction",
    }
    assert all(
        proxy.metadata["proxy_units"] in {PS_PROXY_UNITS, 5} for proxy in proxies
    )

    corrupted = copy.deepcopy(artifact)
    candidate = next(
        record
        for record in corrupted["records"]
        if record["rate_source"]["kind"] == "RMG family estimate"
        and any(
            discovery["event_id"] == record["event_id"]
            and discovery["proxy_site_type"] == record["site_type"]
            for discovery in corrupted["discovery"]
        )
    )
    candidate["raw_path_degeneracy"] += 1.0
    with pytest.raises(AssertionError):
        assert dict(observed_degeneracies(corrupted)) == pytest.approx(expected)


def test_c3_full_2r_plus_3_generation_covers_all_keys_and_128_probe_misses(compilation):
    """C3: full five-unit molecules, all families, and closed-chain donors."""
    _, artifact, _, _, _, _, _, oracle, _ = compilation
    records = {record["event_id"] for record in artifact["records"]}
    compiled = {
        _record_key(entry)
        for entry in artifact["discovery"]
        if entry["event_id"] in records
    }
    exclusions = {
        _record_key(entry)
        for entry in artifact["excluded_channels"]
        if entry["reason"]
        == "bounded J_ring capture replaced by immutable R-009 archived para/ortho pairs"
    }
    required = oracle["keys"] - exclusions
    assert required <= compiled, sorted(required - compiled)[:10]
    baseline = {_record_key(record) for record in _pre_change_artifact()["records"]}
    legacy_sites = {
        "pristine",
        "interior_radical",
        "doubly_featured",
        "end_radical",
        "junction_radical",
        "end_radical+end_radical",
        "junction_radical+end_radical",
        "end_radical+styrene",
    }
    current_missing = {
        key for key in oracle["keys"] if key[0] in legacy_sites
    } - baseline
    archived_replacement = current_missing & exclusions
    assert {key[:2] for key in archived_replacement} == {
        ("junction_radical+end_radical", "R_Recombination")
    }
    probe_missing = current_missing - archived_replacement
    assert len(probe_missing) == 128
    assert probe_missing <= compiled | exclusions

    def graph_keys(records):
        return {
            reaction_graph_key(
                record["family"],
                _molecules(record["reactant_graphs"]),
                _molecules(record["product_graphs"]),
            )
            for record in records
        }

    required_graphs = {
        (family, tuple(tuple(side) for side in sides))
        for site, family, template, sides in oracle["reactions"]
        if (site, family, template) not in exclusions
    }
    actual_graphs = graph_keys(artifact["records"])
    assert required_graphs <= actual_graphs, sorted(required_graphs - actual_graphs)[
        :10
    ]
    print(
        f"C3: full 5-unit oracle keys={len(oracle['keys'])}; coverage={len(required & compiled)}/{len(required)}; probe misses accounted={len(probe_missing & (compiled | exclusions))}/128 (compiled={len(probe_missing & compiled)}, explicitly excluded={len(probe_missing & exclusions)}); exclusions={len(exclusions)}; additional post-probe archived replacement={len(archived_replacement)}"
    )
    archived_j_para = [
        record
        for record in artifact["records"]
        if record.get("junction_ops")
        and record["junction_ops"][0].get("junction_kind") == "J_para"
        and record["junction_ops"][0].get("action") == "create"
    ]
    assert len(archived_j_para) == 1
    assert archived_j_para[0]["rate_source"]["rule_entry_index"] == 176

    corrupted = copy.deepcopy(artifact)
    candidate_key = next(iter(required))
    corrupted["discovery"] = [
        entry for entry in corrupted["discovery"] if _record_key(entry) != candidate_key
    ]
    corrupted_keys = {_record_key(entry) for entry in corrupted["discovery"]}
    assert not required <= corrupted_keys
    missing_pair = next(iter(required_graphs))
    corrupted_records = [
        record
        for record in artifact["records"]
        if reaction_graph_key(
            record["family"],
            _molecules(record["reactant_graphs"]),
            _molecules(record["product_graphs"]),
        )
        != missing_pair
    ]
    assert not required_graphs <= graph_keys(corrupted_records)


def _pack_exemptions(artifact):
    """Match six frozen event IDs exactly, allowing only provenance/handle changes."""
    baseline = {
        record["event_id"]: record for record in _pre_change_artifact()["records"]
    }
    permitted = {
        canonical_json_bytes(_normalized_base_record(baseline[event_id])): event_id
        for event_id in ARCHIVED_PACK_EVENT_IDS
    }
    exemptions = {}
    for record in artifact["records"]:
        content = canonical_json_bytes(_normalized_base_record(record))
        if content in permitted:
            assert record["event_id"] not in exemptions
            exemptions[record["event_id"]] = {
                "reason": ARCHIVED_PACK_EXEMPTION_REASON,
                "pack_event_id": permitted[content],
            }
    assert len(exemptions) == len(ARCHIVED_PACK_EVENT_IDS)
    assert {entry["pack_event_id"] for entry in exemptions.values()} == set(
        ARCHIVED_PACK_EVENT_IDS
    )
    return exemptions


def _assert_pack_exemption(record, exemption):
    baseline = {
        record["event_id"]: record for record in _pre_change_artifact()["records"]
    }
    assert exemption["reason"] == ARCHIVED_PACK_EXEMPTION_REASON
    assert exemption["pack_event_id"] in ARCHIVED_PACK_EVENT_IDS
    assert _normalized_base_record(record) == _normalized_base_record(
        baseline[exemption["pack_event_id"]]
    )


def test_all_pairs_have_exact_graphs_maps_degeneracies_and_detailed_balance(
    compilation, rmg_database
):
    """Re-derive Kc from RMG thermo, not from the compiler's provider or tables."""
    artifact = compilation[1]
    exemptions = _pack_exemptions(artifact)
    by_id = {record["event_id"]: record for record in artifact["records"]}
    thermo_cache = {}
    checked = set()
    paired_graphs = set()
    maximum_error = 0.0

    def thermochemical_species(graph):
        molecule = Molecule().from_adjacency_list(graph)
        key = molecule.to_smiles()
        if key not in thermo_cache:
            species = Species(molecule=[molecule])
            species.generate_resonance_structures()
            species.thermo = rmg_database.thermo.get_thermo_data(species)
            thermo_cache[key] = species
        return thermo_cache[key]

    for record in artifact["records"]:
        if not record["reverse_of"] or record["event_id"] in checked:
            continue
        reverse = by_id[record["reverse_of"]]
        graph_pair = reaction_graph_key(
            record["family"],
            _molecules(record["reactant_graphs"]),
            _molecules(record["product_graphs"]),
        )
        assert reverse["reverse_of"] == record["event_id"]
        assert sorted(record["reactant_graphs"]) == sorted(reverse["product_graphs"])
        assert sorted(record["product_graphs"]) == sorted(reverse["reactant_graphs"])
        assert {
            str(value): int(key) for key, value in record["atom_map"].items()
        } == reverse["atom_map"]
        if record["event_id"] in exemptions:
            _assert_pack_exemption(record, exemptions[record["event_id"]])
            _assert_pack_exemption(reverse, exemptions[reverse["event_id"]])
            checked.update((record["event_id"], reverse["event_id"]))
            continue
        assert reverse["event_id"] not in exemptions
        assert graph_pair not in paired_graphs, graph_pair
        paired_graphs.add(graph_pair)
        assert record["raw_path_degeneracy"] == reverse["raw_path_degeneracy"]
        assert record["degeneracy"] == reverse["degeneracy"]
        direct = [
            candidate
            for candidate in (record, reverse)
            if candidate["rate_source"]["kind"] == "RMG family estimate"
        ]
        assert len(direct) == 1
        forward = direct[0]
        partner = reverse if forward is record else record
        assert partner["rate_source"]["kind"] == "reference-thermo reverse"
        assert forward["thermo_provenance"]["reference_thermo"] == "RMG gas-phase Kc"
        assert (
            forward["thermo_provenance"]["rmg_database_sha"]
            == "4a12d36fcdc193ede82c8d1ab5c1653495d445bc"
        )
        reaction = Reaction(
            reactants=[
                thermochemical_species(graph) for graph in forward["reactant_graphs"]
            ],
            products=[
                thermochemical_species(graph) for graph in forward["product_graphs"]
            ],
        )
        constants = [
            reaction.get_equilibrium_constant(temperature, type="Kc")
            for temperature in (600.0, 700.0, 800.0)
        ]
        for temperature, constant in zip((600.0, 700.0, 800.0), constants):
            index = forward["k_table"]["T"].index(temperature)
            ratio = forward["k_table"]["k"][index] / partner["k_table"]["k"][index]
            error = abs(ratio / constant - 1.0)
            maximum_error = max(maximum_error, error)
            assert error <= 1e-9, (forward["event_id"], temperature, error)
        checked.update((record["event_id"], reverse["event_id"]))
    assert checked
    print(
        f"DETAILED-BALANCE: {(len(checked) - len(exemptions)) // 2} generic pairs at 600/700/800 K; max relative error={maximum_error:.12g}; exemptions={json.dumps(exemptions, sort_keys=True)}"
    )


def test_exact_id_pack_exemptions_preserve_kc_and_reject_non_ring_pair(compilation):
    artifact = compilation[1]
    by_id = {record["event_id"]: record for record in artifact["records"]}
    exemptions = _pack_exemptions(artifact)
    forward_records = [
        by_id[event_id] for event_id in exemptions if by_id[event_id]["arity"] == 2
    ]
    assert len(forward_records) == 3
    for forward in forward_records:
        reverse = by_id[forward["reverse_of"]]
        ortho = forward["junction_ops"][0]["junction_kind"] != "J_para"
        observed = []
        for temperature in (600.0, 650.0, 700.0, 750.0, 800.0):
            index = forward["k_table"]["T"].index(temperature)
            observed.append(
                forward["k_table"]["k"][index]
                / reverse["k_table"]["k"][index]
                * (2.0 if ortho else 1.0)
            )
        assert observed == pytest.approx(
            ORTHO_PACK_KC if ortho else PARA_PACK_KC, rel=1e-9
        )
        assert (forward["degeneracy"], reverse["degeneracy"]) == (
            (1.0, 2.0) if ortho else (1.0, 1.0)
        )
    non_ring = next(
        record
        for record in artifact["records"]
        if record["event_id"] not in exemptions and record["reverse_of"]
    )
    assert non_ring["inventory_class"] != "R1:J_ring"
    disguised = copy.deepcopy(non_ring)
    disguised["inventory_class"] = "R1:J_ring"
    disguised["event_id"] = forward_records[0]["event_id"]
    with pytest.raises(AssertionError):
        _assert_pack_exemption(disguised, exemptions[disguised["event_id"]])


def test_real_frontier_recombination_stays_refused_and_out_of_normal_total(
    rmg_database,
):
    from compile_event_set_fixture import (
        dump_generated_reactions,
        load_generated_reactions,
    )
    from rmgpy.kmc.met import TRANSPORT_ARMS, _is_forward_met_record, compile_bulk_table
    from rmgpy.kmc.ssa import SiteIndex, channel_propensities
    import stateTest as state_oracle

    proxy = next(
        proxy
        for proxy in ps_proxy_set(3)
        if proxy.site_type == "end_radical+end_radical"
    )
    frontier = replace(proxy, frontier=True)
    generated = rmg_database.kinetics.generate_reactions_from_families(
        [species.copy(deep=True) for species in frontier.reactants],
        only_families=["R_Recombination"],
        resonance=True,
    )
    recovered = load_generated_reactions(dump_generated_reactions(generated))
    assert [
        atom.id
        for reaction in generated
        for species in reaction.reactants + reaction.products
        for molecule in species.molecule
        for atom in molecule.atoms
    ] == [
        atom.id
        for reaction in recovered
        for species in reaction.reactants + reaction.products
        for molecule in species.molecule
        for atom in molecule.atoms
    ]
    artifact = EventSetCompiler(
        rmg_database.kinetics,
        [frontier],
        ["R_Recombination"],
        thermo_database=rmg_database.thermo,
        database_path=DB_PATH,
        reaction_cache={frontier.site_type: recovered},
    ).compile()
    real = [
        record
        for record in artifact["records"]
        if record["inventory_class"] != "R1:J_ring"
    ]
    assert real
    assert all(record["status"] == "refused" for record in real)
    assert all(not _is_forward_met_record(record, "R0") for record in real)
    table = compile_bulk_table(
        artifact["records"], TRANSPORT_ARMS["A0_REF_H_CROSS_NEc"], "R0"
    )
    assert not table.channels
    record = next(record for record in real if record["arity"] == 2)
    state, _, _ = state_oracle._harness(record)
    index = SiteIndex(state, real)
    report = channel_propensities(state, index, real, 700.0, 1e-24)
    assert report.total_enabled == 0.0
    assert report.total_refused > 0.0
    mutated = copy.deepcopy(record)
    mutated["status"] = "irreversible"
    assert _is_forward_met_record(mutated, "R0")


def test_c7_c9_inventory_and_real_artifact(compilation):
    """Check C7/C9 inventory, real ceiling, and ring reverse records."""
    _, artifact, _, _, _, _, artifact_path, _, _ = compilation
    catalogue = short_ps_molecule_catalogue(artifact["inputs"]["span_radius"])
    assert catalogue == artifact["short_molecule_catalogue"]
    assert catalogue["terminates"]
    assert catalogue["size"] == len(catalogue["molecules"])
    irreversible = {entry["event_id"] for entry in artifact["irreversible_pairs"]}
    assert all(
        record["reverse_of"]
        or record["event_id"] in irreversible
        or record["status"] == "refused"
        for record in artifact["records"]
    )
    assert all(
        entry["reason"].startswith("thermo unavailable for species ")
        or entry["reason"].startswith("Kc unavailable for species ")
        or entry["reason"] == "RMG reaction is declared irreversible"
        for entry in artifact["irreversible_pairs"]
    )
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
    literature = {
        "temperature_K": 668.15,
        "range_K": [583.15, 668.15],
        "reference": "Kinetic Phenomena in Mechanochemical Depolymerization of Poly(styrene), ACS Sustainable Chemistry & Engineering, DOI:10.1021/acssuschemeng.3c05296, section 4.1 reports 310-395 C",
        "reference_state": "literature range (upper endpoint quoted); compiled value uses gas-phase Kc and the stated monomer concentration, not the same solvent/activity",
    }
    linked_pairs = (
        sum(bool(record["reverse_of"]) for record in artifact["records"]) // 2
    )
    one_way_records = sum(not record["reverse_of"] for record in artifact["records"])
    print(
        f"C9: PS gas-reference ceiling={recomputed:.6f} K at [styrene]={pair['monomer_concentration_mol_m3']} mol/m3; literature={literature['temperature_K']} K ({literature['reference']}); {literature['reference_state']}; irreversible={len(irreversible)}/{linked_pairs + one_way_records} reaction pairs"
    )
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
