"""Verification of the frozen R-010 v13 volatility law."""

from collections import defaultdict
from dataclasses import FrozenInstanceError
from functools import lru_cache
import hashlib
import math
import os
from pathlib import Path
import subprocess
import sys

import pytest
from rdkit import Chem

from rmgpy.molecule import Molecule
import rmgpy.kmc.volatility as volatility
from rmgpy.kmc.compiler import short_ps_molecule_catalogue
from rmgpy.kmc.volatility import evaluate


REPORT = Path("/home/alon/Code/polymer-pm/reports/R-010_v13_2026-09-28_volatility-numerics-prereg.md")
REPORT_SHA256 = "10994e5189ec9e91d750d2acf7eb9230f3d74ab9c342d42232e2ba85aba0caf3"
SOURCE_ROOT = Path("/home/alon/Code/polymer-pm/reports/R-010_sources")
LOG10P_ATOL = 5e-12
TBTC_ATOL = 5e-9
PC_RTOL = 5e-12


def _molecule(smiles):
    return Molecule().from_smiles(smiles)


def _alkane(n):
    return _molecule("C" * n)


def _radical_alkane(n):
    return _molecule("[CH2]" + "C" * (n - 1))


def _read_pressure_rows():
    data = REPORT.read_bytes()
    assert hashlib.sha256(data).hexdigest() == REPORT_SHA256
    lines = data.decode().splitlines()
    header = "id|M_g_mol|n_ar|T_K|outcome|log10P_bar|Tb_K|Tc_K|Pc_bar|groups|regime|flags"
    start = lines.index(header)
    rows = {}
    for line in lines[start + 1 :]:
        if line.startswith("R1_FIXTURES|"):
            break
        fields = line.split("|")
        rows[fields[0]] = fields
    assert len(rows) == 32
    return rows


def _fixture(case_id):
    fixed = {
        "G01M": "C",
        "A01": "C=Cc1ccccc1",
        "A02": "C=C(C)c1ccccc1",
        "G03": "Cc1ccccc1",
        "G08D": "c1ccccc1Cc2ccccc2",
        "G127": "Cc1ccccc1C",
        "G128": "Cc1cccc(C)c1",
        "G59": "c1ccccc1C=Cc2ccccc2",
        "G130": "CC(C)(C)c1ccccc1",
        "G130X": "C=CC(C)(C)c1ccccc1",
        "G131": "CC(C)C(C)C",
        "G132": "CC(C)(C)C(C)C",
        "G133": "CC(C)(C)C(C)(C)C",
        "CA-": "C=Cc1ccccc1",
        "CA+": "C=Cc1ccccc1",
    }
    if case_id in fixed:
        return _molecule(fixed[case_id]), {}
    alkane_sizes = {
        "S01": 10,
        "S02": 14,
        "S03": 20,
        "S04": 28,
        "S05": 35,
        "S06": 36,
        "X01": 37,
        "X02": 50,
        "TB-": 20,
        "TB0": 20,
        "TB+": 20,
        "C-": 10,
        "C+": 10,
        "B36": 36,
        "O01": 107,
        "OC+": 107,
    }
    if case_id == "RCO+":
        return _radical_alkane(107), {"off_above_1500": True}
    kwargs = {"off_above_1500": True} if case_id in {"O01", "OC+"} else {}
    return _alkane(alkane_sizes[case_id]), kwargs


def _optional_float(text):
    return None if text == "-" else float(text)


def _groups(text):
    if not text:
        return ()
    return tuple((int(group), int(count)) for group, count in (item.split(":") for item in text.split(",")))


@lru_cache(maxsize=1)
def _source_smarts_patterns():
    patterns = {}
    for filename in ("nannoolal_primary.smarts", "nannoolal_secondary.smarts"):
        patterns.update(
            {
                int(group): Chem.MolFromSmarts(smarts)
                for group, smarts in (
                    line.split(None, 1)
                    for line in (SOURCE_ROOT / filename).read_text().splitlines()
                    if line.strip()
                )
            }
        )
    return patterns


def _source_smarts_counts(structure, *, include_group8_variants=True):
    smiles = structure.to_smiles() if isinstance(structure, Molecule) else structure
    molecule = Chem.MolFromSmiles(smiles)
    patterns = _source_smarts_patterns()
    matches = defaultdict(
        int,
        {
            group: len(molecule.GetSubstructMatches(pattern, uniquify=True))
            for group, pattern in patterns.items()
        },
    )
    counts = {
        1: matches[1] + matches[7000] + 2 * matches[7001],
        3: matches[3],
        4: matches[4],
        5: matches[5],
        6: matches[6],
        8: matches[8] + (matches[801] + matches[802] if include_group8_variants else 0),
        9: matches[9],
        10: matches[10],
        11: matches[11],
        15: matches[15],
        16: matches[16],
        58: matches[58] + matches[5801] + matches[5802],
        59: matches[59] + matches[5901] + matches[5902] + matches[5903],
        61: matches[61] + matches[6101] + matches[6102],
        62: sum(matches[group] for group in (62, 6201, 6202, 6203, 6204, 6205, 7006)),
        88: sum(matches[group] for group in (88, 7006, 8801, 8802, 8803, 8804, 8805, 8806, 8807, 8808, 8809)),
        89: sum(matches[group] for group in (89, 8901, 8902, 8903, 8904, 8905, 8906, 8907, 8908, 8909,
                                             8910, 8911, 8912, 8913, 8914, 8915, 8916, 8917, 8918, 8919)),
        127: int(matches[127] == 1 and matches[128] == 0 and matches[129] == 0),
        128: int(matches[127] == 0 and matches[128] in {1, 3} and matches[129] == 0),
        129: int(matches[127] == 0 and matches[128] == 0 and matches[129] == 2),
        130: matches[130],
        131: matches[131] + matches[1311] + matches[1312] + matches[1313],
        132: matches[132] + matches[1321],
        133: matches[133],
    }
    return tuple(sorted((group, count) for group, count in counts.items() if count))


def _closed_shell_ps_catalogue():
    catalogue = short_ps_molecule_catalogue(radius=3, feature_cap=2)
    return [
        (
            f"PS-{entry['repeat_units']}",
            Molecule().from_adjacency_list(entry["adjacency_list"]),
        )
        for entry in catalogue["molecules"]
        if entry["state"] == "closed_shell_linear"
    ]


def _production_differential_corpus():
    corpus = _closed_shell_ps_catalogue()
    corpus.extend((f"n-C{size}", _alkane(size)) for size in range(1, 41))
    corpus.extend(
        (f"{size}-alkylbenzene", _molecule("C" * size + "c1ccccc1"))
        for size in range(1, 11)
    )
    corpus.extend(
        (
            ("1,2-dimethylbenzene", _molecule("Cc1ccccc1C")),
            ("1,3-dimethylbenzene", _molecule("Cc1cccc(C)c1")),
            ("1,4-dimethylbenzene", _molecule("Cc1ccc(C)cc1")),
            ("1,3,5-trimethylbenzene", _molecule("Cc1cc(C)cc(C)c1")),
        )
    )
    return corpus


def _assert_row(case_id, row=None):
    row = row or _read_pressure_rows()[case_id]
    molecule, kwargs = _fixture(case_id)
    result = evaluate(molecule, float(row[3]), **kwargs)
    assert result.outcome == row[4]
    assert result.error == "-"
    assert result.M == float(row[1])
    assert result.n_ar == int(row[2])
    assert result.groups == _groups(row[9])
    assert result.regime == row[10]
    assert result.flags == (() if row[11] == "-" else tuple(row[11].split(",")))
    expected_logp = _optional_float(row[5])
    expected_tb = _optional_float(row[6])
    expected_tc = _optional_float(row[7])
    expected_pc = _optional_float(row[8])
    if expected_logp is None:
        assert result.log10P is None
    else:
        assert result.log10P == pytest.approx(expected_logp, abs=LOG10P_ATOL)
    for actual, expected in ((result.Tb, expected_tb), (result.Tc, expected_tc)):
        if expected is None:
            assert actual is None
        else:
            assert actual == pytest.approx(expected, abs=TBTC_ATOL)
    if expected_pc is None:
        assert result.Pc is None
    else:
        assert result.Pc == pytest.approx(expected_pc, rel=PC_RTOL)


def test_golden_pressure_table():
    for case_id, row in _read_pressure_rows().items():
        _assert_row(case_id, row)


@pytest.mark.parametrize(
    "case_id",
    ("TB-", "TB0", "TB+", "C-", "C+", "CA-", "CA+", "B36", "O01", "OC+", "RCO+"),
)
def test_boundary_and_precedence_rows(case_id):
    _assert_row(case_id)


def test_extrapolated_mass_depends_only_on_mass():
    light_radical = evaluate(_radical_alkane(10), 600.0)
    closed_c37 = evaluate(_alkane(37), 800.0)
    assert light_radical.outcome == "retained_radical"
    assert "extrapolated_M" not in light_radical.flags
    assert closed_c37.outcome == "P_sat"
    assert "extrapolated_M" in closed_c37.flags


@pytest.mark.parametrize(
    "temperature,option,code",
    (
        ("600", False, "E_TYPE_T"),
        (math.inf, False, "E_NONFINITE_T"),
        (math.nan, False, "E_NONFINITE_T"),
        (0.0, False, "E_RANGE_T"),
        (600.0, "false", "E_TYPE_OPTION"),
    ),
)
def test_input_errors(temperature, option, code):
    assert evaluate(_alkane(10), temperature, off_above_1500=option).error == code


@pytest.mark.parametrize(
    "smiles,code",
    (
        ("c1ccccc1-c2ccccc2", "E_DOMAIN_FUSED_BIARYL"),
        ("C1CCCCC1", "E_DOMAIN_ALIPHATIC_RING"),
        ("C=CC=C", "E_DOMAIN_CONJUGATED_ALKENE"),
        ("C=C=C", "E_DOMAIN_CUMULATED_ALKENE"),
        ("CC.CC", "E_DOMAIN_DISCONNECTED"),
    ),
)
def test_graph_domain_errors(smiles, code):
    assert evaluate(_molecule(smiles), 550.0).error == code


def test_psat_over_pc_is_last_rejection():
    isooctane = _molecule("CC(C)CC(C)(C)C")
    Tc = evaluate(isooctane, 300.0).Tc
    assert evaluate(isooctane, Tc - 1e-6).error == "E_PSAT_GT_PC"


def test_critical_group_domain_error(monkeypatch):
    monkeypatch.setattr(volatility, "_TC", {group: -abs(value) for group, value in volatility._TC.items()})
    assert evaluate(_alkane(10), 550.0).error == "E_CRITICAL_CORRELATION_DOMAIN"


def test_retained_routing_and_frozen_record():
    j_result = evaluate(_molecule("C1CCCCC1"), 550.0, r1_class="J_para")
    assert j_result.outcome == "retained_J_ring"
    assert j_result.groups == ()
    with pytest.raises(FrozenInstanceError):
        j_result.outcome = "changed"
    radical_j = evaluate(_molecule("[CH]1CCCCC1"), 550.0, r1_class="J_ortho")
    assert radical_j.outcome == "retained_radical"
    radical_biaryl = evaluate(_molecule("[c]1ccccc1-c2ccccc2"), 550.0)
    assert radical_biaryl.outcome == "retained_radical"
    j_biaryl = evaluate(_molecule("c1ccccc1-c2ccccc2"), 550.0, r1_class="J_para")
    assert j_biaryl.outcome == "retained_J_ring"


def test_source_smarts_finite_sample_counts():
    samples = {
        "C": ((1, 1),),
        "CCCCCCCCCC": ((1, 2), (4, 8)),
        "c1ccccc1": ((15, 6),),
        "Cc1ccccc1C": ((3, 2), (15, 4), (16, 2), (127, 1)),
        "c1ccccc1Cc2ccccc2": ((8, 2), (15, 10), (16, 2)),
        "C1=CC=CCC1C=CC=CCC2CCCCC2": ((4, 1), (9, 6), (10, 2), (88, 1), (89, 1)),
    }
    for smiles, expected in samples.items():
        molecule = _molecule(smiles)
        r1_class = "quinoid_free" if smiles.startswith("C1=") else None
        compiled = volatility._compile_graph(molecule, r1_class)
        assert compiled[2] == _source_smarts_counts(smiles)
        assert compiled[2] == expected


@pytest.mark.skipif(os.environ.get("RMG_KMC_SLOW") != "1", reason="set RMG_KMC_SLOW=1")
def test_source_smarts_production_differential():
    corpus = _production_differential_corpus()
    exclusions = []
    compared = 0
    for label, molecule in corpus:
        production = volatility._compile_graph(molecule, None)[2]
        source = _source_smarts_counts(molecule)
        production_groups = {group for group, _ in production}
        source_groups = {group for group, _ in source}
        unsupported = tuple(sorted((production_groups | source_groups) & {11, 62}))
        if unsupported:
            exclusions.append((label, unsupported))
            continue
        assert production == source, label
        compared += 1
    assert len(corpus) == 61
    assert compared + len(exclusions) == len(corpus)
    print(
        "PSAT_SMARTS_DIFFERENTIAL|"
        f"requested={len(corpus)}|compared={compared}|excluded={len(exclusions)}|"
        f"exclusions={exclusions or '-'}"
    )


@pytest.mark.skipif(os.environ.get("RMG_KMC_SLOW") != "1", reason="set RMG_KMC_SLOW=1")
def test_source_smarts_wrong_mapping_control():
    mismatches = []
    catalogue = _closed_shell_ps_catalogue()
    for label, molecule in catalogue:
        production = volatility._compile_graph(molecule, None)[2]
        deliberately_wrong = _source_smarts_counts(
            molecule,
            include_group8_variants=False,
        )
        if production != deliberately_wrong:
            mismatches.append(label)
    assert mismatches
    print(
        "PSAT_SMARTS_CONTROL|"
        f"catalogue={len(catalogue)}|mismatches={len(mismatches)}|structures={mismatches}"
    )


def test_atom_order_does_not_change_output():
    molecule = _molecule("CC(C)(C)c1ccccc1")
    expected = evaluate(molecule, 350.0)
    molecule.atoms.reverse()
    assert evaluate(molecule, 350.0) == expected


def test_coefficient_mutation_is_killed(monkeypatch):
    row = _read_pressure_rows()["S01"]
    changed = dict(volatility._TC)
    changed[1] += 1e-3
    monkeypatch.setattr(volatility, "_TC", changed)
    with pytest.raises(AssertionError):
        _assert_row("S01", row)


def test_tb_branch_mutation_is_killed(monkeypatch):
    row = _read_pressure_rows()["TB-"]
    monkeypatch.setattr(volatility, "_below_tb", lambda T, Tb: T >= Tb)
    with pytest.raises(AssertionError):
        _assert_row("TB-", row)


@pytest.mark.skipif(os.environ.get("RMG_KMC_SLOW") != "1", reason="set RMG_KMC_SLOW=1")
def test_python_hash_seed_determinism():
    script = (
        "from rmgpy.molecule import Molecule; "
        "from rmgpy.kmc.volatility import evaluate; "
        "print(repr(evaluate(Molecule().from_smiles('CC(C)(C)c1ccccc1'),350.0)))"
    )
    outputs = []
    for seed in ("0", "1", "8675309"):
        environment = os.environ.copy()
        environment["PYTHONHASHSEED"] = seed
        outputs.append(subprocess.check_output([sys.executable, "-c", script], text=True, env=environment))
    assert len(set(outputs)) == 1


# RMG Molecule inputs contain no external cached mass, aromatic-count, group,
# or radical-count fields.  Therefore ERRORS rows E02-E04, E06-E38 are not
# expressible through this API; E01, E05 and E39-E48 are exercised above.
