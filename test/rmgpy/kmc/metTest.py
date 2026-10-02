"""Contract tests for the melt encounter-limited termination kernel."""

import math
import os
import hashlib
import json
import subprocess
import sys
from dataclasses import replace
from itertools import permutations
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import rmgpy.kmc.met as met_module
from rmgpy.data.rmg import RMGDatabase
from rmgpy.kmc.compiler import validate_artifact
from rmgpy.kmc.met import (
    CageBalanceError,
    CoverageClass,
    CoverageError,
    ESCAPED_RADICAL_PAIR,
    FAMILY_CLASSIFICATION,
    FamilyClassification,
    FrozenMetricVerdict,
    MappedGraph,
    R0_DISPROPORTIONATION_PRODUCTS,
    TRANSPORT_ARMS,
    RateTable,
    TransportArm,
    UnsupportedTopologyError,
    apply_compiled_graph_rewrite,
    all_arm_conjunction,
    audit_met_coverage,
    collins_kimball,
    compare_compiled_tables,
    compiled_graph_postcondition_failures,
    compile_bulk_table,
    compile_cage_table,
    compile_r0_fallback_stress,
    diffusion_rate,
    draw_ortho_attack_site,
    enumerate_radical_decreasing_families,
    enumerate_radical_decreasing_library_reactions,
    graph_postcondition_failures,
    reconcile_cage_balance,
    record_activation_rate,
    rmg_reference_equilibrium_constant,
    solve_cage,
    validate_met_topology,
)
from rmgpy.kmc.event_record import EventRecord


TEMPERATURES = (600.0, 650.0, 700.0, 750.0, 800.0)
REPO_ROOT = Path(__file__).resolve().parents[3]
REAL_DATABASE_PATH = Path(
    os.environ.get("RMG_DATABASE_PATH", "/home/alon/Code/RMG-database")
)
REAL_CACHE_ROOT = Path(
    os.environ.get(
        "RMG_KMC_CACHE_ROOT", str(Path(__file__).with_name(".real-event-cache"))
    )
)
ARCHIVED_D0 = {
    "H": (
        6.636772327361e-9,
        1.512876971830e-8,
        2.755363626148e-8,
        4.346438783277e-8,
        6.218869377200e-8,
    ),
    "U": (
        2.888355626300e-9,
        5.194737616032e-9,
        8.113610273553e-9,
        1.151642486291e-8,
        1.527388428419e-8,
    ),
    "B": (
        2.667092489713e-8,
        6.959334540015e-8,
        1.583409029823e-7,
        3.228606625277e-7,
        6.022258564396e-7,
    ),
}

ARCHIVED_RATES = {
    "D1": (
        2.85002820943072380e5,
        3.27835038939649821e5,
        3.68520160926779383e5,
        4.06775334046629316e5,
        4.42476593677184777e5,
    ),
    "D2": (
        2.73305346453879378e5,
        2.99018641244220489e5,
        3.22000156835251546e5,
        3.42442600462969684e5,
        3.60563723405797558e5,
    ),
    "D3": (
        1.55901020305904839e6,
        1.49186594037869922e6,
        1.43228097483286751e6,
        1.37894978137219185e6,
        1.33086085842475668e6,
    ),
    "D4": (2.9e6, 2.9e6, 2.9e6, 2.9e6, 2.9e6),
    "D5": (
        2.34449584058339105e-5,
        1.34358601442878577e-4,
        6.00018871714558841e-4,
        2.19488154723838967e-3,
        6.82733939244908102e-3,
    ),
    "null": (
        7.03220236239306442e6,
        6.72546468392766826e6,
        6.45341129656218272e6,
        6.21003669605328795e6,
        5.99069032098846417e6,
    ),
    "J_para": (
        2.89220714646790028e7,
        2.66910790121819414e7,
        2.47792290474214815e7,
        2.31226376634507552e7,
        2.16734020085316598e7,
    ),
    "J_ortho": (
        5.78441429293580055e7,
        5.33821580243638828e7,
        4.95584580948429629e7,
        4.62452753269015104e7,
        4.33468040170633197e7,
    ),
}
ARCHIVED_KC = {
    "J_para": (
        1.04541002791187644e8,
        4.66322438450995088e6,
        3.28434940662126872e5,
        3.33060720810081402e4,
        4.53840953775354592e3,
    ),
    "J_ortho": (
        6.57178252159645688e6,
        3.77346219583960657e5,
        3.29954009752024649e4,
        4.03633396767985278e3,
        6.48237362204028727e2,
    ),
}
EXPECTED_FATES = (
    1.689288479229e-1,
    6.846389753251e-3,
    6.565390887287e-3,
    3.745082748346e-2,
    6.966432900113e-2,
    5.631990675179e-13,
    7.105442149514e-1,
)
EXPECTED_MEANS = {
    "G": 2.402218241418e-8,
    "J_para": 2.511303038812,
    "J_ortho": 1.578685585201e-1,
    "absorption_time": 2.669171621354,
    "cycles": 2.084313829562,
}
SSA_FATE_TOLERANCES = (
    6.209652555666e-3,
    1.366580783966e-3,
    1.338431730996e-3,
    3.146577687350e-3,
    4.219114901730e-3,
    1.243733266766e-8,
    7.515934005542e-3,
)
SSA_MEAN_TOLERANCES = {
    "G": 3.981151894680e-10,
    "J_para": 8.196631986289e-2,
    "J_ortho": 4.086265213050e-3,
    "absorption_time": 8.338438019610e-2,
    "cycles": 4.202010399341e-2,
}


def _git_head(path):
    return subprocess.check_output(
        ["git", "-C", str(path), "rev-parse", "HEAD"], text=True
    ).strip()


def real_slow(test_function):
    test_function = pytest.mark.slow(test_function)
    return pytest.mark.skipif(
        os.environ.get("RMG_KMC_SLOW") != "1" or not REAL_DATABASE_PATH.is_dir(),
        reason="set RMG_KMC_SLOW=1 and provide RMG_DATABASE_PATH",
    )(test_function)


@pytest.fixture(scope="session")
def real_ps_inputs():
    """Compile/cache one validated real PS artifact and its real kinetics DB."""
    from rmgpy.kmc.compiler import compiler_source_hash

    cache_key = f"{_git_head(REPO_ROOT)}-{_git_head(REAL_DATABASE_PATH)}-{compiler_source_hash()}"
    cache_dir = REAL_CACHE_ROOT / cache_key
    cache_dir.mkdir(parents=True, exist_ok=True)
    artifacts = sorted(cache_dir.glob("*.json"))
    if not artifacts:
        fixture_script = Path(__file__).with_name("compile_event_set_fixture.py")
        stdout_path = cache_dir / "compile.stdout.log"
        stderr_path = cache_dir / "compile.stderr.log"
        environment = os.environ.copy()
        environment.update(
            {
                "PYTHONPATH": str(REPO_ROOT),
                "PYTHONHASHSEED": "0",
                "MPLCONFIGDIR": str(cache_dir / "matplotlib"),
            }
        )
        with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
            subprocess.run(
                [
                    sys.executable,
                    str(fixture_script),
                    str(REAL_DATABASE_PATH),
                    str(cache_dir),
                ],
                cwd=REPO_ROOT,
                env=environment,
                stdout=stdout,
                stderr=stderr,
                check=True,
                timeout=4 * 60 * 60,
            )
        artifacts = sorted(cache_dir.glob("*.json"))
    assert len(artifacts) == 1
    artifact_path = artifacts[0]
    payload = artifact_path.read_bytes()
    assert hashlib.sha256(payload).hexdigest() == artifact_path.stem
    artifact = json.loads(payload)
    validate_artifact(artifact)

    database = RMGDatabase()
    database.load_kinetics(
        str(REAL_DATABASE_PATH / "input/kinetics"),
        reaction_libraries=None,
        seed_mechanisms=None,
        kinetics_families=["all"],
        kinetics_depositories=[],
    )
    return SimpleNamespace(
        artifact=artifact,
        records=tuple(artifact["records"]),
        database=database,
        artifact_path=artifact_path,
        cache_key=cache_key,
    )


def _table(rates):
    return {
        "T": list(TEMPERATURES),
        "k": list(rates),
        "interpolation": "linear-ln-k",
        "extrapolation": "refuse",
    }


def _record_bond_ops(channel_id):
    if channel_id != "D2":
        return [{"channel": channel_id, "action": "archived"}]
    return [
        {"action": "break", "atoms": [0, 1], "order": "1.0"},
        {"action": "form", "atoms": [0, 1], "order": "2.0"},
        {"action": "set_implicit_hydrogens", "atom": 0, "value": 0},
        {"action": "set_radical", "atom": 1, "value": 0},
        {"action": "set_implicit_hydrogens", "atom": 2, "value": 2},
        {"action": "set_radical", "atom": 2, "value": 0},
    ]


def _compiled_records():
    records = []
    channels = [
        (channel_id, channel_id, None)
        for channel_id in ("null", "D1", "D2", "D3", "D4", "D5", "J_para")
    ] + [
        ("J_ortho", "J_ortho", "S7"),
        ("J_ortho-S6", "J_ortho", "S6"),
    ]
    for event_name, channel_id, ortho_label in channels:
        ring = channel_id.startswith("J_")
        reverse_id = f"reverse-{event_name}" if ring else None
        forward_rates = (
            tuple(rate / 2.0 for rate in ARCHIVED_RATES[channel_id])
            if channel_id == "J_ortho"
            else ARCHIVED_RATES[channel_id]
        )
        records.append(
            SimpleNamespace(
                event_id=f"forward-{event_name}",
                family=(
                    "R_Recombination"
                    if channel_id in {"null", "J_para", "J_ortho"}
                    else "Disproportionation"
                ),
                arity=2,
                radical_delta=-2,
                k_table=_table(forward_rates),
                ssa_multiplier=1.0,
                inventory_class="R1:J_ring" if ring else None,
                reverse_of=reverse_id,
                rate_source={},
                atom_map={0: 0, 1: 1},
                bond_ops=_record_bond_ops(channel_id),
                junction_ops=(
                    [
                        {
                            "action": "create",
                            "junction_kind": (
                                "J_para"
                                if channel_id == "J_para"
                                else f"J_ortho_{ortho_label}"
                            ),
                            "attacked_atom_label": (
                                "S9" if channel_id == "J_para" else ortho_label
                            ),
                            "attacker_site_label": "P9",
                            "formed_bond": {
                                "label_pair": [
                                    "S9" if channel_id == "J_para" else ortho_label,
                                    "P9",
                                ],
                                "order": "1.0",
                            },
                            "chosen_resonance_localisation": (
                                "para:S9"
                                if channel_id == "J_para"
                                else f"ortho:{ortho_label}"
                            ),
                            "exact_inverse_rewrite": [
                                {
                                    "action": "break",
                                    "atoms": [
                                        "S9" if channel_id == "J_para" else ortho_label,
                                        "P9",
                                    ],
                                    "order": "1.0",
                                },
                                {
                                    "action": "set_radical",
                                    "atom": (
                                        "S9" if channel_id == "J_para" else ortho_label
                                    ),
                                    "value": 1,
                                },
                                {
                                    "action": "set_radical",
                                    "atom": "P9",
                                    "value": 1,
                                },
                            ],
                            "aromaticity": "dearomatised",
                            "eligibility_refresh": True,
                        }
                    ]
                    if ring
                    else None
                ),
                product_graphs=[f"archived-{channel_id}-product"],
                status="enabled",
            )
        )
        if ring:
            reverse_rates = tuple(
                forward / equilibrium
                for forward, equilibrium in zip(
                    ARCHIVED_RATES[channel_id], ARCHIVED_KC[channel_id]
                )
            )
            records.append(
                SimpleNamespace(
                    event_id=reverse_id,
                    family="R_Recombination",
                    arity=1,
                    radical_delta=2,
                    k_table=_table(reverse_rates),
                    ssa_multiplier=1.0,
                    inventory_class="R1:J_ring",
                    reverse_of=f"forward-{event_name}",
                    rate_source={},
                    atom_map={0: 0, 1: 1},
                    bond_ops=[{"channel": channel_id, "action": "inverse"}],
                    junction_ops=[
                        {
                            "action": "dissociate",
                            "junction_kind": (
                                "J_para"
                                if channel_id == "J_para"
                                else f"J_ortho_{ortho_label}"
                            ),
                            "attacked_atom_label": (
                                "S9" if channel_id == "J_para" else ortho_label
                            ),
                            "attacker_site_label": "P9",
                            "formed_bond": {
                                "label_pair": [
                                    "S9" if channel_id == "J_para" else ortho_label,
                                    "P9",
                                ],
                                "order": "1.0",
                            },
                            "aromaticity": "restore",
                            "eligibility_refresh": True,
                        }
                    ],
                    product_graphs=[f"archived-{channel_id}-pair"],
                    status="enabled",
                )
            )
    return tuple(records)


def _replace_channel(table, channel_id, **changes):
    channels = tuple(
        replace(channel, **changes) if channel.channel_id == channel_id else channel
        for channel in table.channels
    )
    return replace(table, channels=channels)


def _without_channel(table, channel_id):
    return replace(
        table,
        channels=tuple(
            channel for channel in table.channels if channel.channel_id != channel_id
        ),
    )


@pytest.mark.parametrize("arm", TRANSPORT_ARMS.values(), ids=TRANSPORT_ARMS)
def test_transport_arms_round_trip_archived_d0(arm):
    actual = tuple(arm.d0(temperature) for temperature in TEMPERATURES)
    assert actual == pytest.approx(ARCHIVED_D0[arm.package], rel=1e-12)


def test_transport_arms_have_frozen_identifiers_and_no_default():
    assert tuple(TRANSPORT_ARMS) == (
        "A0_REF_H_CROSS_NEc",
        "A1_H_ROUSE",
        "A2_H_CROSS_NElo",
        "A3_H_CROSS_NEhi",
        "A4_LU_ROUSE",
        "A5_LU_CROSS_NElo",
        "A6_LU_CROSS_NEhi",
        "A7_LB_ROUSE",
        "A8_LB_CROSS_NElo",
        "A9_LB_CROSS_NEhi",
    )
    assert all(arm.arm_id for arm in TRANSPORT_ARMS.values())


def test_collins_kimball_and_identical_pair_convention():
    table = RateTable(
        temperatures=(600.0, 700.0),
        rates=(5.0, 10.0),
    )
    record = type(
        "Record",
        (),
        {"k_table": table.to_mapping(), "ssa_multiplier": 2.0},
    )()
    k_act = record_activation_rate(record, 700.0)
    assert k_act == 20.0
    assert collins_kimball(k_act, 30.0) == pytest.approx(12.0)
    assert collins_kimball(math.inf, 30.0) == 30.0


def test_diffusion_rate_uses_injected_arm_and_pair_class():
    arm = TRANSPORT_ARMS["A5_LU_CROSS_NElo"]
    expected = (
        4.0
        * math.pi
        * 6.02214076e23
        * (arm.chain_diffusivity(700.0, 1000.0) * 2.0)
        * max(6.766081442101e-10, 2.0 * math.sqrt(0.434e-20 * 104.15 * 1000 / 6.0))
    )
    assert diffusion_rate(
        arm, 700.0, 1000, 1000, "end/end", spin_factor=1.0
    ) == pytest.approx(expected)
    with pytest.raises(ValueError, match="pair class"):
        diffusion_rate(arm, 700.0, 10, 10, "unknown", spin_factor=1.0)


def test_met_coverage_accounts_for_every_report_family_and_compiled_record():
    records = (
        SimpleNamespace(
            event_id="recombination",
            family="R_Recombination",
            radical_delta=-2,
            arity=2,
            status="enabled",
            rate_source={},
        ),
        SimpleNamespace(
            event_id="intra",
            family="Intra_Disproportionation",
            radical_delta=-2,
            arity=1,
            status="enabled",
            rate_source={},
        ),
        SimpleNamespace(
            event_id="held",
            family="R_Addition_MultipleBond_Disprop",
            radical_delta=-2,
            arity=2,
            status="refused",
            rate_source={},
        ),
    )
    report = audit_met_coverage(records, FAMILY_CLASSIFICATION)
    assert set(report.classified_families) == set(FAMILY_CLASSIFICATION)
    assert report.admitted_event_ids == ("recombination",)
    assert (
        report.classified_families["R_Recombination"].category
        is CoverageClass.BIMolecular_MET
    )
    assert (
        report.classified_families["Intra_Disproportionation"].category
        is CoverageClass.UNIMOLECULAR
    )


def test_met_coverage_dropped_family_control_and_library_audit():
    dropped = dict(FAMILY_CLASSIFICATION)
    dropped.pop("Disproportionation-Y")
    with pytest.raises(CoverageError, match="Disproportionation-Y"):
        audit_met_coverage((), FAMILY_CLASSIFICATION, classification=dropped)
    library_record = SimpleNamespace(
        event_id="library",
        family="",
        radical_delta=-2,
        arity=2,
        status="enabled",
        rate_source={"library_reaction": "lib:pair-1"},
    )
    with pytest.raises(CoverageError, match="unclassified.*library"):
        audit_met_coverage((library_record,), (), {})
    with pytest.raises(CoverageError, match="held out.*library"):
        audit_met_coverage(
            (library_record,),
            (),
            {"lib:pair-1": FAMILY_CLASSIFICATION["R_Addition_MultipleBond_Disprop"]},
        )


def test_met_coverage_enumerates_loaded_database_recipes_and_library_entries():
    losing_recipe = SimpleNamespace(
        actions=[("LOSE_RADICAL", "*1", 1), ("LOSE_RADICAL", "*2", 1)]
    )
    database = SimpleNamespace(
        families={
            "R_Recombination": SimpleNamespace(forward_recipe=losing_recipe),
            "New_Termination": SimpleNamespace(forward_recipe=losing_recipe),
        },
        libraries={},
    )
    with pytest.raises(CoverageError, match="New_Termination"):
        audit_met_coverage((), database)

    def species(radicals):
        atoms = [SimpleNamespace(radical_electrons=1) for _ in range(radicals)]
        return SimpleNamespace(molecule=[SimpleNamespace(atoms=atoms)])

    reaction = SimpleNamespace(
        reactants=[species(1), species(1)], products=[species(0)]
    )
    database.families.pop("New_Termination")
    database.libraries = {
        "loaded": SimpleNamespace(entries={"pair": SimpleNamespace(item=reaction)})
    }
    with pytest.raises(CoverageError, match="loaded:pair"):
        audit_met_coverage((), database)


@real_slow
def test_met_coverage_on_real_database_and_compiled_ps_artifact(real_ps_inputs):
    kinetics = real_ps_inputs.database.kinetics
    radical_decreasing_libraries = set(
        enumerate_radical_decreasing_library_reactions(kinetics)
    )
    assert radical_decreasing_libraries
    library_classification = {
        reaction_id: FamilyClassification(
            CoverageClass.OUT_OF_SCOPE,
            "global kinetics-library reaction is absent from the scoped PS proxy set",
        )
        for reaction_id in radical_decreasing_libraries
    }
    report = audit_met_coverage(
        real_ps_inputs.records,
        kinetics,
        library_classification,
    )
    radical_decreasing = set(enumerate_radical_decreasing_families(kinetics))
    assert set(report.classified_families) == radical_decreasing
    assert set(report.library_reactions) == radical_decreasing_libraries
    admitted = set(report.admitted_event_ids)
    assert admitted
    assert all(
        record["event_id"] in admitted
        for record in real_ps_inputs.records
        if record["radical_delta"] < 0
        and record["arity"] == 2
        and record["family"] in {"R_Recombination", "Disproportionation"}
        and record["status"] == "enabled"
    )

    dropped_name = sorted(radical_decreasing)[0]
    dropped = dict(FAMILY_CLASSIFICATION)
    dropped.pop(dropped_name)
    with pytest.raises(CoverageError, match=dropped_name):
        audit_met_coverage(
            real_ps_inputs.records,
            kinetics,
            library_classification,
            classification=dropped,
        )


def test_met_topology_rejects_branched_and_gel_components_with_typed_error():
    validate_met_topology({"topology": "linear", "branch_points": 0, "is_gel": False})
    with pytest.raises(UnsupportedTopologyError, match="branched"):
        validate_met_topology({"topology": "branched", "branch_points": 1})
    with pytest.raises(UnsupportedTopologyError, match="gel"):
        validate_met_topology({"topology": "linear", "is_gel": True})


def test_r0_r1_are_one_shared_inventory_rule_on_separate_compile_paths():
    records = _compiled_records()
    arm = TRANSPORT_ARMS["A4_LU_ROUSE"]
    cage_r0 = compile_cage_table(records, arm, "R0")
    bulk_r0 = compile_bulk_table(records, arm, "R0")
    cage_r1 = compile_cage_table(records, arm, "R1")
    bulk_r1 = compile_bulk_table(records, arm, "R1")
    assert cage_r0 is not bulk_r0
    assert cage_r0.channels is not bulk_r0.channels
    assert set(cage_r0.indexed_channels()) == {"null", "D1", "D2", "D3", "D4", "D5"}
    assert set(cage_r1.indexed_channels()) == {
        "null",
        "D1",
        "D2",
        "D3",
        "D4",
        "D5",
        "J_para",
        "J_ortho",
    }
    assert compare_compiled_tables(cage_r0, bulk_r0) == ()
    assert compare_compiled_tables(cage_r1, bulk_r1) == ()


def test_ortho_declared_channel_rejects_para_summed_rate():
    records = []
    partitioned_para = _table(tuple(rate / 2.0 for rate in ARCHIVED_RATES["J_para"]))
    for record in _compiled_records():
        operation = record.junction_ops[0] if record.junction_ops else {}
        if operation.get("action") == "create" and operation.get(
            "junction_kind", ""
        ).startswith("J_ortho"):
            record = SimpleNamespace(
                **{
                    **vars(record),
                    "k_table": partitioned_para,
                    "rate_source": {"channel_id": "J_ortho"},
                }
            )
        records.append(record)

    with pytest.raises(ValueError, match="J_ortho.*disagrees|disagrees.*J_ortho"):
        compile_cage_table(records, TRANSPORT_ARMS["A4_LU_ROUSE"], "R1")


def test_ortho_declared_channel_rejects_unrecognized_summed_rate():
    records = []
    partitioned_invalid = _table((1.0, 1.0, 1.0, 1.0, 1.0))
    for record in _compiled_records():
        operation = record.junction_ops[0] if record.junction_ops else {}
        if operation.get("action") == "create" and operation.get(
            "junction_kind", ""
        ).startswith("J_ortho"):
            record = SimpleNamespace(
                **{
                    **vars(record),
                    "k_table": partitioned_invalid,
                    "rate_source": {"channel_id": "J_ortho"},
                }
            )
        records.append(record)

    with pytest.raises(ValueError, match="J_ortho.*does not match"):
        compile_cage_table(records, TRANSPORT_ARMS["A4_LU_ROUSE"], "R1")


def test_para_declared_channel_rejects_unrecognized_rate():
    records = []
    for record in _compiled_records():
        operation = record.junction_ops[0] if record.junction_ops else {}
        if (
            operation.get("action") == "create"
            and operation.get("junction_kind") == "J_para"
        ):
            record = SimpleNamespace(
                **{
                    **vars(record),
                    "k_table": _table((1.0, 1.0, 1.0, 1.0, 1.0)),
                    "rate_source": {"channel_id": "J_para"},
                }
            )
        records.append(record)

    with pytest.raises(ValueError, match="J_para.*does not match"):
        compile_cage_table(records, TRANSPORT_ARMS["A4_LU_ROUSE"], "R1")


def test_inventory_domain_and_refused_record_admission_are_enforced():
    records = list(_compiled_records())
    with pytest.raises(ValueError, match="inventory must"):
        compile_cage_table(records, TRANSPORT_ARMS["A4_LU_ROUSE"], "R2")
    refused = SimpleNamespace(
        **{
            **vars(records[0]),
            "event_id": "refused-extra",
            "status": "refused",
            "rate_source": {"channel_id": "refused-extra"},
        }
    )
    table = compile_bulk_table(records + [refused], TRANSPORT_ARMS["A4_LU_ROUSE"], "R1")
    assert "refused-extra" not in table.indexed_channels()


def test_bulk_table_uses_compiled_activation_rates_and_injected_transport():
    records = _compiled_records()
    arm = TRANSPORT_ARMS["A1_H_ROUSE"]
    bulk = compile_bulk_table(records, arm, "R0")
    k_act = sum(
        ARCHIVED_RATES[channel][2] for channel in ("null", "D1", "D2", "D3", "D4", "D5")
    )
    k_diffusion = diffusion_rate(arm, 700.0, 100.0, 100.0, "end/mid", spin_factor=0.25)
    assert bulk.bulk_rate(
        700.0, 100.0, 100.0, "end/mid", spin_factor=0.25
    ) == pytest.approx(collins_kimball(k_act, k_diffusion))


@real_slow
def test_real_r0_tables_use_compiled_rmg_rates_for_every_transport_arm(
    real_ps_inputs,
):
    records_by_id = {record["event_id"]: record for record in real_ps_inputs.records}
    expected_event_ids = {
        record["event_id"]
        for record in real_ps_inputs.records
        if record["radical_delta"] < 0
        and record["arity"] == 2
        and record["family"] in {"R_Recombination", "Disproportionation"}
        and record["inventory_class"] != "R1:J_ring"
        and record["status"] != "refused"
        and record["k_table"]
    }
    assert expected_event_ids
    for arm in TRANSPORT_ARMS.values():
        cage = compile_cage_table(real_ps_inputs.records, arm, "R0")
        bulk = compile_bulk_table(real_ps_inputs.records, arm, "R0")
        for table in (cage, bulk):
            assert {
                channel.event_id for channel in table.channels
            } == expected_event_ids
            for channel in table.channels:
                record = records_by_id[channel.event_id]
                for temperature in TEMPERATURES:
                    assert channel.forward_rate(temperature) == pytest.approx(
                        record_activation_rate(record, temperature), rel=1.0e-12
                    )
        assert [channel.event_id for channel in cage.channels] == [
            channel.event_id for channel in bulk.channels
        ]

    identical_pair = next(
        record
        for record in real_ps_inputs.records
        if record["site_type"] == "end_radical+end_radical"
        and record["family"] == "R_Recombination"
        and record["radical_delta"] < 0
    )
    assert identical_pair["raw_path_degeneracy"] == 0.5
    assert identical_pair["ssa_multiplier"] == 2.0


@real_slow
def test_real_r1_tables_compile_ortho_for_every_transport_arm(real_ps_inputs):
    error_type = met_module.MissingMETChannelError
    for arm in TRANSPORT_ARMS.values():
        for compile_table in (compile_cage_table, compile_bulk_table):
            table = compile_table(real_ps_inputs.records, arm, "R1")
            assert "J_ortho" in table.indexed_channels()

    ortho_forwards = {
        record["junction_ops"][0]["attacked_atom_label"]: record
        for record in real_ps_inputs.records
        if record.get("junction_ops")
        and record["junction_ops"][0].get("action") == "create"
        and record["junction_ops"][0].get("junction_kind", "").startswith("J_ortho")
    }
    assert set(ortho_forwards) == {"S6", "S7"}
    arm = TRANSPORT_ARMS["A4_LU_ROUSE"]
    for missing_label, missing_record in sorted(ortho_forwards.items()):
        incomplete = tuple(
            record
            for record in real_ps_inputs.records
            if record["event_id"] != missing_record["event_id"]
        )
        for compile_table in (compile_cage_table, compile_bulk_table):
            with pytest.raises(error_type, match="J_ortho") as error:
                compile_table(incomplete, arm, "R1")
            assert error.value.channel_id == "J_ortho", missing_label


def test_exact_cage_observables_and_fate_semantics_match_preregistered_targets():
    cage = compile_cage_table(_compiled_records(), TRANSPORT_ARMS["A4_LU_ROUSE"], "R1")
    observed = solve_cage(cage, 600.0)
    assert observed.terminal_sequence == (
        "null",
        "D1",
        "D2",
        "D3",
        "D4",
        "D5",
        "escape",
    )
    for actual, expected in zip(observed.terminal_probabilities, EXPECTED_FATES):
        assert abs(math.log(actual / expected)) <= 1.0e-10
    for state in ("G", "J_para", "J_ortho"):
        assert observed.occupancies[state] == pytest.approx(
            EXPECTED_MEANS[state], rel=1e-10
        )
    assert observed.mean_absorption_time == pytest.approx(
        EXPECTED_MEANS["absorption_time"], rel=1e-10
    )
    assert observed.expected_cycles == pytest.approx(
        EXPECTED_MEANS["cycles"], rel=1e-10
    )
    assert observed.f_capture == pytest.approx(0.675778778, rel=1e-8)
    assert observed.f_redissociate == 1.0
    assert observed.f_escape == observed.terminal_probabilities[-1]
    assert observed.f_permanent == pytest.approx(
        sum(observed.terminal_probabilities[1:6])
    )


def test_continuous_rmg_reference_reverse_evaluator_matches_ramp_archive():
    cage = compile_cage_table(_compiled_records(), TRANSPORT_ARMS["A4_LU_ROUSE"], "R1")
    channels = cage.indexed_channels()
    targets = {
        600.0: (1.045410027912e8, 2.766576816032e-1, 8.801895488670),
        625.0: (2.071026575002e7, 1.340489342565, 3.739927269518e1),
        700.0: (3.284349406621e5, 7.544638520331e1, 1.501980780051e3),
        725.0: (1.004078742947e5, 2.382515322651e2, 4.304634388454e3),
        800.0: (4.538409537754e3, 4.775550075029e3, 6.686872208304e4),
        825.0: (1.840057291734e3, 1.142070438999e4, 1.483911699937e5),
        850.0: (7.882824828745e2, 2.587257501442e4, 3.132908096697e5),
    }
    for temperature, (para_kc, para_kdiss, ortho_kdiss) in targets.items():
        assert rmg_reference_equilibrium_constant(
            temperature, "J_para"
        ) == pytest.approx(para_kc, rel=1e-11)
        assert channels["J_para"].dissociation_rate(temperature) == pytest.approx(
            para_kdiss, rel=1e-11
        )
        assert channels["J_ortho"].dissociation_rate(temperature) == pytest.approx(
            ortho_kdiss, rel=1e-11
        )


def test_r0_omits_capture_and_has_no_cycles():
    cage = compile_cage_table(_compiled_records(), TRANSPORT_ARMS["A4_LU_ROUSE"], "R0")
    observed = solve_cage(cage, 600.0)
    assert observed.f_capture == 0.0
    assert observed.f_redissociate == 0.0
    assert observed.expected_cycles == 0.0
    assert set(observed.occupancies) == {"G"}


def test_terminal_j_negative_control_reports_required_fate_mass():
    cage = compile_cage_table(_compiled_records(), TRANSPORT_ARMS["A4_LU_ROUSE"], "R1")
    with pytest.raises(CageBalanceError) as caught:
        solve_cage(cage, 600.0, dissociation_scale=0.0)
    assert caught.value.required_terminal_fate_mass == pytest.approx(
        0.3242212223722, rel=1e-12
    )


PAIR_NODES = {
    "P1": ("C", 0, 1),
    "P2": ("C", 0, 3),
    "P3": ("C", 0, 0),
    "P4": ("C", 0, 1),
    "P5": ("C", 0, 1),
    "P6": ("C", 0, 1),
    "P7": ("C", 0, 1),
    "P8": ("C", 0, 1),
    "P9": ("C", 1, 2),
    "S1": ("C", 0, 2),
    "S2": ("C", 0, 2),
    "S3": ("C", 0, 3),
    "S4": ("C", 0, 0),
    "S5": ("C", 0, 1),
    "S6": ("C", 0, 1),
    "S7": ("C", 0, 1),
    "S8": ("C", 0, 1),
    "S9": ("C", 1, 1),
    "S10": ("C", 0, 1),
}
PAIR_EDGES = {
    ("P1", "P2"): 1.0,
    ("P1", "P3"): 1.0,
    ("P1", "P9"): 1.0,
    ("P3", "P4"): 1.5,
    ("P3", "P5"): 1.5,
    ("P4", "P6"): 1.5,
    ("P5", "P8"): 1.5,
    ("P6", "P7"): 1.5,
    ("P7", "P8"): 1.5,
    ("S1", "S2"): 1.0,
    ("S1", "S3"): 1.0,
    ("S2", "S5"): 1.0,
    ("S4", "S5"): 2.0,
    ("S4", "S6"): 1.0,
    ("S4", "S7"): 1.0,
    ("S6", "S8"): 2.0,
    ("S8", "S9"): 1.0,
    ("S9", "S10"): 1.0,
    ("S7", "S10"): 2.0,
}
PARA_J_NODES = {
    "P1": ("C", 0, 1),
    "P2": ("C", 0, 3),
    "P3": ("C", 0, 0),
    "P4": ("C", 0, 1),
    "P5": ("C", 0, 1),
    "P6": ("C", 0, 1),
    "P7": ("C", 0, 1),
    "P8": ("C", 0, 1),
    "P9": ("C", 0, 2),
    "S1": ("C", 0, 2),
    "S2": ("C", 0, 2),
    "S3": ("C", 0, 3),
    "S4": ("C", 0, 0),
    "S5": ("C", 0, 1),
    "S6": ("C", 0, 1),
    "S7": ("C", 0, 1),
    "S8": ("C", 0, 1),
    "S9": ("C", 0, 1),
    "S10": ("C", 0, 1),
}
PARA_J_EDGES = {
    ("P1", "P2"): 1.0,
    ("P1", "P3"): 1.0,
    ("P1", "P9"): 1.0,
    ("P3", "P4"): 1.5,
    ("P3", "P5"): 1.5,
    ("P4", "P6"): 1.5,
    ("P5", "P8"): 1.5,
    ("P6", "P7"): 1.5,
    ("P7", "P8"): 1.5,
    ("P9", "S9"): 1.0,
    ("S1", "S2"): 1.0,
    ("S1", "S3"): 1.0,
    ("S2", "S5"): 1.0,
    ("S4", "S5"): 2.0,
    ("S4", "S6"): 1.0,
    ("S4", "S7"): 1.0,
    ("S6", "S8"): 2.0,
    ("S8", "S9"): 1.0,
    ("S9", "S10"): 1.0,
    ("S7", "S10"): 2.0,
}
ORTHO_J_NODES = {
    "P1": ("C", 0, 1),
    "P2": ("C", 0, 3),
    "P3": ("C", 0, 0),
    "P4": ("C", 0, 1),
    "P5": ("C", 0, 1),
    "P6": ("C", 0, 1),
    "P7": ("C", 0, 1),
    "P8": ("C", 0, 1),
    "P9": ("C", 0, 2),
    "S1": ("C", 0, 2),
    "S2": ("C", 0, 2),
    "S3": ("C", 0, 3),
    "S4": ("C", 0, 0),
    "S5": ("C", 0, 1),
    "S6": ("C", 0, 1),
    "S7": ("C", 0, 1),
    "S8": ("C", 0, 1),
    "S9": ("C", 0, 1),
    "S10": ("C", 0, 1),
}
ORTHO_J_EDGES = {
    ("P1", "P2"): 1.0,
    ("P1", "P3"): 1.0,
    ("P1", "P9"): 1.0,
    ("P3", "P4"): 1.5,
    ("P3", "P5"): 1.5,
    ("P4", "P6"): 1.5,
    ("P5", "P8"): 1.5,
    ("P6", "P7"): 1.5,
    ("P7", "P8"): 1.5,
    ("P9", "S7"): 1.0,
    ("S1", "S2"): 1.0,
    ("S1", "S3"): 1.0,
    ("S2", "S5"): 1.0,
    ("S4", "S5"): 2.0,
    ("S4", "S6"): 1.0,
    ("S4", "S7"): 1.0,
    ("S6", "S8"): 2.0,
    ("S8", "S9"): 1.0,
    ("S9", "S10"): 2.0,
    ("S7", "S10"): 1.0,
}
D2_NODES = {
    "P1": ("C", 0, 0),
    "P2": ("C", 0, 3),
    "P3": ("C", 0, 0),
    "P4": ("C", 0, 1),
    "P5": ("C", 0, 1),
    "P6": ("C", 0, 1),
    "P7": ("C", 0, 1),
    "P8": ("C", 0, 1),
    "P9": ("C", 0, 2),
    "S1": ("C", 0, 2),
    "S2": ("C", 0, 2),
    "S3": ("C", 0, 3),
    "S4": ("C", 0, 0),
    "S5": ("C", 0, 1),
    "S6": ("C", 0, 1),
    "S7": ("C", 0, 1),
    "S8": ("C", 0, 1),
    "S9": ("C", 0, 2),
    "S10": ("C", 0, 1),
}
D2_EDGES = {
    ("P1", "P2"): 1.0,
    ("P1", "P3"): 1.0,
    ("P1", "P9"): 2.0,
    ("P3", "P4"): 1.5,
    ("P3", "P5"): 1.5,
    ("P4", "P6"): 1.5,
    ("P5", "P8"): 1.5,
    ("P6", "P7"): 1.5,
    ("P7", "P8"): 1.5,
    ("S1", "S2"): 1.0,
    ("S1", "S3"): 1.0,
    ("S2", "S5"): 1.0,
    ("S4", "S5"): 2.0,
    ("S4", "S6"): 1.0,
    ("S4", "S7"): 1.0,
    ("S6", "S8"): 2.0,
    ("S8", "S9"): 1.0,
    ("S9", "S10"): 1.0,
    ("S7", "S10"): 2.0,
}


def _ring_mirror_from_topology(edges):
    ring = ("S4", "S6", "S8", "S9", "S10", "S7")
    ring_edges = {frozenset((a, b)) for a, b in edges if a in ring and b in ring}
    candidates = []
    for image in permutations(ring):
        mapping = dict(zip(ring, image))
        if mapping["S4"] != "S4" or mapping["S9"] != "S9" or mapping["S7"] != "S6":
            continue
        mapped = {
            frozenset((mapping[a], mapping[b]))
            for edge in ring_edges
            for a, b in (tuple(edge),)
        }
        if mapped == ring_edges:
            candidates.append(mapping)
    assert len(candidates) == 1
    return candidates[0]


def test_mapped_para_and_ortho_capture_inverse_and_d2_routing():
    records = _compiled_records()
    by_event = {record.event_id: record for record in records}
    para_record = by_event["forward-J_para"]
    ortho_record = by_event["forward-J_ortho"]
    d2_record = by_event["forward-D2"]
    para_pair = MappedGraph(PAIR_NODES, PAIR_EDGES)
    captured = apply_compiled_graph_rewrite(para_record, para_pair)
    expected_para = MappedGraph(PARA_J_NODES, PARA_J_EDGES)
    assert (
        graph_postcondition_failures("J_para_capture_graph", captured, expected_para)
        == ()
    )
    restored = apply_compiled_graph_rewrite(para_record, captured, inverse=True)
    assert (
        graph_postcondition_failures("J_para_dissociation_graph", restored, para_pair)
        == ()
    )

    ortho_nodes = dict(PAIR_NODES)
    ortho_nodes["S9"] = ("C", 0, 1)
    ortho_nodes["S7"] = ("C", 1, 1)
    ortho_edges = dict(PAIR_EDGES)
    ortho_edges[("S7", "S10")] = 1.0
    ortho_edges[("S9", "S10")] = 2.0
    ortho_pair = MappedGraph(ortho_nodes, ortho_edges)
    ortho_s7 = apply_compiled_graph_rewrite(
        ortho_record, ortho_pair, attacked_atom_id="S7"
    )
    expected_ortho_s7 = MappedGraph(ORTHO_J_NODES, ORTHO_J_EDGES)
    assert (
        graph_postcondition_failures(
            "J_ortho_S7_capture_graph", ortho_s7, expected_ortho_s7
        )
        == ()
    )
    ortho_s7_restored = apply_compiled_graph_rewrite(
        ortho_record, ortho_s7, attacked_atom_id="S7", inverse=True
    )
    assert (
        graph_postcondition_failures(
            "J_ortho_S7_dissociation_graph", ortho_s7_restored, ortho_pair
        )
        == ()
    )
    reflection = _ring_mirror_from_topology(ortho_edges)
    ortho_s6_pair = ortho_pair.relabel(reflection)
    ortho_s6 = apply_compiled_graph_rewrite(
        ortho_record, ortho_s6_pair, attacked_atom_id="S6"
    )
    expected_ortho_s6 = expected_ortho_s7.relabel(reflection)
    assert (
        graph_postcondition_failures(
            "J_ortho_S6_capture_graph", ortho_s6, expected_ortho_s6
        )
        == ()
    )
    ortho_s6_restored = apply_compiled_graph_rewrite(
        ortho_record, ortho_s6, attacked_atom_id="S6", inverse=True
    )
    assert (
        graph_postcondition_failures(
            "J_ortho_S6_dissociation_graph", ortho_s6_restored, ortho_s6_pair
        )
        == ()
    )

    d2 = apply_compiled_graph_rewrite(
        d2_record, para_pair, atom_labels={0: "P1", 1: "P9", 2: "S9"}
    )
    assert (
        graph_postcondition_failures(
            "D2_product_graph", d2, MappedGraph(D2_NODES, D2_EDGES)
        )
        == ()
    )


def test_real_event_record_indices_and_product_graph_postcondition_are_consumed():
    record = EventRecord(
        atom_map={0: 0, 1: 1},
        bond_ops=[
            {"action": "form", "atoms": [0, 1], "order": "1.0"},
            {"action": "set_radical", "atom": 0, "value": 0},
            {"action": "set_radical", "atom": 1, "value": 0},
        ],
        product_graphs=["1 C u0 p0 c0 {2,S}\n2 C u0 p0 c0 {1,S}\n"],
    )
    pair = MappedGraph({"secondary": ("C", 1, 3), "primary": ("C", 1, 3)}, {})
    product = apply_compiled_graph_rewrite(
        record, pair, atom_labels={0: "secondary", 1: "primary"}
    )
    assert (
        compiled_graph_postcondition_failures(
            "EventRecord_product_graph",
            record,
            product,
            {0: "secondary", 1: "primary"},
        )
        == ()
    )
    corrupted = replace(
        record,
        product_graphs=["1 C u1 p0 c0\n", "1 C u1 p0 c0\n"],
        event_id="",
    )
    assert compiled_graph_postcondition_failures(
        "EventRecord_product_graph",
        corrupted,
        product,
        {0: "secondary", 1: "primary"},
    ) == ("EventRecord_product_graph",)


def test_wrong_para_mapping_and_corrupt_ortho_bond_are_observable():
    records = {record.event_id: record for record in _compiled_records()}
    pair = MappedGraph(PAIR_NODES, PAIR_EDGES)
    expected = apply_compiled_graph_rewrite(records["forward-J_para"], pair)
    wrong_operation = dict(records["forward-J_para"].junction_ops[0])
    wrong_operation["attacked_atom_label"] = "S7"
    wrong_record = SimpleNamespace(
        **{**vars(records["forward-J_para"]), "junction_ops": [wrong_operation]}
    )
    wrong = apply_compiled_graph_rewrite(wrong_record, pair)
    assert graph_postcondition_failures("J_para_capture_graph", wrong, expected) == (
        "J_para_capture_graph",
    )
    ortho_nodes = dict(PAIR_NODES)
    ortho_nodes["S9"] = ("C", 0, 1)
    ortho_nodes["S7"] = ("C", 1, 1)
    ortho_edges = dict(PAIR_EDGES)
    ortho_edges[("S7", "S10")] = 1.0
    ortho_edges[("S9", "S10")] = 2.0
    ortho_pair = MappedGraph(ortho_nodes, ortho_edges)
    ortho = apply_compiled_graph_rewrite(records["forward-J_ortho"], ortho_pair)
    corrupt_operation = dict(records["forward-J_ortho"].junction_ops[0])
    corrupt_operation["formed_bond"] = {
        **corrupt_operation["formed_bond"],
        "order": "2.0",
    }
    corrupt_record = SimpleNamespace(
        **{
            **vars(records["forward-J_ortho"]),
            "junction_ops": [corrupt_operation],
        }
    )
    corrupt = apply_compiled_graph_rewrite(corrupt_record, ortho_pair)
    assert graph_postcondition_failures("J_ortho_product_graph", corrupt, ortho) == (
        "J_ortho_product_graph",
    )


@pytest.mark.parametrize("inventory", ("R0", "R1"))
def test_shared_inventory_behavior_matches_for_all_transport_arms(inventory):
    records = _compiled_records()
    for arm in TRANSPORT_ARMS.values():
        cage = compile_cage_table(records, arm, inventory)
        bulk = compile_bulk_table(records, arm, inventory)
        assert compare_compiled_tables(cage, bulk) == ()


def test_shared_inventory_negative_controls_have_exact_named_failures():
    records = _compiled_records()
    arm = TRANSPORT_ARMS["A5_LU_CROSS_NElo"]
    cage = compile_cage_table(records, arm, "R1")
    bulk = compile_bulk_table(records, arm, "R1")
    bulk_r0 = compile_bulk_table(records, arm, "R0")
    assert compare_compiled_tables(cage, bulk_r0) == (
        "missing_bulk_channel:J_ortho",
        "missing_bulk_channel:J_para",
    )

    original_para = bulk.indexed_channels()["J_para"]
    changed_rate = RateTable(
        original_para.rate.temperatures,
        tuple(value * 2.0 for value in original_para.rate.rates),
    )
    assert compare_compiled_tables(
        cage, _replace_channel(bulk, "J_para", rate=changed_rate)
    ) == ("rate_object:J_para",)
    assert compare_compiled_tables(
        cage, _replace_channel(bulk, "J_ortho", reversible=False)
    ) == ("reversibility:J_ortho",)
    assert compare_compiled_tables(
        cage,
        _replace_channel(
            bulk,
            "J_para",
            rewrite_graph=original_para.rewrite_graph + (("corrupt",),),
        ),
    ) == ("rewrite_graph:J_para",)
    restricted_cage = _replace_channel(
        cage, "J_para", restrictions=frozenset({"cage_only_contact"})
    )
    assert compare_compiled_tables(restricted_cage, bulk) == (
        "unlisted_restriction:J_para:cage_only_contact",
    )


def test_shared_inventory_allow_list_nonfinite_missing_and_transport_controls():
    records = _compiled_records()
    arm = TRANSPORT_ARMS["A5_LU_CROSS_NElo"]
    cage = compile_cage_table(records, arm, "R1")
    bulk = compile_bulk_table(records, arm, "R1")
    missing = _without_channel(bulk, "J_para")
    assert (
        compare_compiled_tables(
            cage,
            missing,
            allowed_restrictions={
                ("bulk", "J_para", "unavailable"): "independent control justification"
            },
        )
        == ()
    )
    assert compare_compiled_tables(cage, missing) == ("missing_bulk_channel:J_para",)

    corrupt_rate = object.__new__(RateTable)
    object.__setattr__(corrupt_rate, "temperatures", tuple(TEMPERATURES))
    object.__setattr__(corrupt_rate, "rates", (math.nan,) * len(TEMPERATURES))
    assert compare_compiled_tables(
        cage, _replace_channel(bulk, "J_para", rate=corrupt_rate)
    ) == ("rate_object:J_para",)

    miswired = TransportArm(
        arm.arm_id,
        arm.package,
        arm.scaling,
        TRANSPORT_ARMS["A6_LU_CROSS_NEhi"].n_e,
    )
    assert compare_compiled_tables(cage, replace(bulk, arm=miswired)) == (
        "transport_Ne:A5_LU_CROSS_NElo",
    )


def test_duplicate_channel_is_rejected_before_comparison():
    records = _compiled_records()
    arm = TRANSPORT_ARMS["A4_LU_ROUSE"]
    cage = compile_cage_table(records, arm, "R1")
    bulk = compile_bulk_table(records, arm, "R1")
    duplicate = replace(bulk, channels=bulk.channels + (bulk.channels[0],))
    with pytest.raises(ValueError, match="duplicate channel"):
        compare_compiled_tables(cage, duplicate)


def test_binding_cage_reconciliation_and_all_arm_conjunction():
    r1 = reconcile_cage_balance("R1", False)
    assert r1.status == "B4-unresolved"
    assert r1.fallback_grid == ()

    required = reconcile_cage_balance("R0", False)
    assert required.status == "fallback-stress-required"
    assert required.fallback_grid == (0.0, 0.25, 0.5, 0.75, 1.0)
    rates = dict(zip(R0_DISPROPORTIONATION_PRODUCTS, (1.0, 3.0)))
    stress = compile_r0_fallback_stress(12.0, rates)
    assert tuple(point.f_const for point in stress) == required.fallback_grid
    assert all(point.k_fiss == 12.0 for point in stress)
    assert all(point.escaped_pair_topology == ESCAPED_RADICAL_PAIR for point in stress)
    assert stress[0].escaped_pair_probability == 0.0
    assert stress[-1].escaped_pair_probability == 1.0
    assert stress[2].disproportionation_probabilities == pytest.approx(
        {
            R0_DISPROPORTIONATION_PRODUCTS[0]: 0.125,
            R0_DISPROPORTIONATION_PRODUCTS[1]: 0.375,
        }
    )
    with pytest.raises(ValueError, match="exactly the two archived"):
        compile_r0_fallback_stress(12.0, {"placeholder": 1.0})

    verdict = FrozenMetricVerdict(
        True,
        sign="positive",
        direction="increasing",
        ranking=("transport", "baseline"),
        interpretation="robust",
    )
    executed = []

    def run_passing(point):
        executed.append(point)
        return verdict

    assert (
        reconcile_cage_balance(
            "R0",
            False,
            fallback_k_fiss=12.0,
            disproportionation_rates=rates,
            run_fallback_case=run_passing,
        ).status
        == "fallback-stress"
    )
    assert tuple(point.f_const for point in executed) == required.fallback_grid

    def run_discordant(point):
        if point.f_const == 0.5:
            return replace(verdict, direction="decreasing")
        return verdict

    assert (
        reconcile_cage_balance(
            "R0",
            False,
            fallback_k_fiss=12.0,
            disproportionation_rates=rates,
            run_fallback_case=run_discordant,
        ).status
        == "transport-nonidentified"
    )

    all_pass = {arm_id: verdict for arm_id in TRANSPORT_ARMS}
    assert all_arm_conjunction(all_pass)
    all_pass["A9_LB_CROSS_NEhi"] = replace(verdict, ranking=("baseline", "transport"))
    assert not all_arm_conjunction(all_pass)
    all_pass.pop("A9_LB_CROSS_NEhi")
    with pytest.raises(ValueError, match="all-arm conjunction"):
        all_arm_conjunction(all_pass)


def _ssa_observables(cage, dissociation_scale=1.0):
    channels = cage.indexed_channels()
    terminal_names = ("null", "D1", "D2", "D3", "D4", "D5")
    terminal_rates = [channels[name].forward_rate(600.0) for name in terminal_names]
    terminal_rates.append(
        4.0 * math.pi * 6.02214076e23 * 6.766081442101e-10 * 2.0 * cage.arm.d0(600.0)
    )
    capture_names = ("J_para", "J_ortho")
    capture_rates = [channels[name].forward_rate(600.0) for name in capture_names]
    if isinstance(dissociation_scale, dict):
        scales = [dissociation_scale.get(name, 1.0) for name in capture_names]
    else:
        scales = [dissociation_scale, dissociation_scale]
    reverse_rates = [
        channels[name].dissociation_rate(600.0) * scale
        for name, scale in zip(capture_names, scales)
    ]
    g_rates = np.asarray(terminal_rates + capture_rates)
    probabilities = g_rates / g_rates.sum()
    rng = np.random.default_rng(270927)
    fate_counts = np.zeros(7, dtype=int)
    dwell_sums = np.zeros(3)
    cycles = 0
    for _ in range(131072):
        while True:
            dwell_sums[0] += rng.exponential(1.0 / g_rates.sum())
            event = int(rng.choice(len(g_rates), p=probabilities))
            if event < 7:
                fate_counts[event] += 1
                break
            state = event - 7 + 1
            cycles += 1
            dwell_sums[state] += rng.exponential(1.0 / reverse_rates[state - 1])
    means = {
        "G": dwell_sums[0] / 131072,
        "J_para": dwell_sums[1] / 131072,
        "J_ortho": dwell_sums[2] / 131072,
        "absorption_time": dwell_sums.sum() / 131072,
        "cycles": cycles / 131072,
    }
    return fate_counts / 131072, means


def _ssa_failures(fates, means):
    failures = []
    for label, actual, expected, tolerance in zip(
        ("null", "D1", "D2", "D3", "D4", "D5", "escape"),
        fates,
        EXPECTED_FATES,
        SSA_FATE_TOLERANCES,
    ):
        if abs(actual - expected) > tolerance:
            failures.append(label)
    for label in ("G", "J_para", "J_ortho", "absorption_time", "cycles"):
        if abs(means[label] - EXPECTED_MEANS[label]) > SSA_MEAN_TOLERANCES[label]:
            failures.append(label)
    return tuple(failures)


@pytest.mark.skipif(
    os.environ.get("RMG_KMC_SLOW") != "1",
    reason="set RMG_KMC_SLOW=1 for the pre-registered 131072-replica cage verifier",
)
def test_slow_cage_ssa_rows_rate_controls_and_ortho_draw():
    cage = compile_cage_table(_compiled_records(), TRANSPORT_ARMS["A4_LU_ROUSE"], "R1")
    fates, means = _ssa_observables(cage)
    assert _ssa_failures(fates, means) == ()

    fates, means = _ssa_observables(cage, 10.0)
    assert _ssa_failures(fates, means) == (
        "J_para",
        "J_ortho",
        "absorption_time",
    )

    fates, means = _ssa_observables(cage, {"J_para": 10.0})
    assert _ssa_failures(fates, means) == ("J_para", "absorption_time")

    rng = np.random.default_rng(270927)
    draws = [draw_ortho_attack_site(rng) for _ in range(131072)]
    assert set(draws) == {"S6", "S7"}
    frequency = draws.count("S6") / 131072
    assert abs(frequency - 0.5) <= 0.008286407592
    assert cage.indexed_channels()["J_ortho"].forward_rate(600.0) == pytest.approx(
        ARCHIVED_RATES["J_ortho"][0], rel=1e-11
    )


@pytest.mark.skip(
    reason="MET-KERNEL-REFERENCE deferred: no executable Rouse/reptation reference defined in the termination report; manager tracking"
)
def test_met_kernel_reference():
    pass
