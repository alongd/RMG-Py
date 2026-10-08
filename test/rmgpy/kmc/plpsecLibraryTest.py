"""Small mapped-inventory verification of the local propagation library."""

import hashlib
from pathlib import Path
import math
from types import SimpleNamespace

import pytest

from rmgpy import constants
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.kinetics.arrhenius import Arrhenius
from rmgpy.kmc import compiler
from rmgpy.kmc.kinetics_library import load_plpsec_entry, matches_head_to_tail, plpsec_rate_table
from rmgpy.kmc.reference_thermo import ReferenceThermoResult
from rmgpy.molecule.molecule import Bond, Molecule
from rmgpy.species import Species
from benzylicEndTest import SmallReferenceThermo, mapped_addition, small_inventory


GRID = [600.0, 700.0, 800.0]


def _instance(*, use_plpsec_library=None, database=None):
    return compiler.EventSetCompiler(
        database, [], [], temperature_grid=GRID,
        reference_thermo_provider=SmallReferenceThermo(),
        use_plpsec_library=use_plpsec_library,
    )


def _proxy(reaction):
    # Deliberately misleading declarations/templates cannot set matching scope.
    return compiler.SiteProxy("benzylic_end_radical+styrene", reaction.reactants,
                              participant_site_types=("benzylic_end_radical", "styrene"))


def _pair(reaction, **options):
    return _instance(**options)._linked_family_pair(_proxy(reaction), reaction, {})[0]


def _ring_addition():
    chain = Molecule(smiles="CC(c1ccccc1)C[CH](c1ccccc1)")
    styrene = Molecule(smiles="C=Cc1ccccc1")
    for identifier, atom in enumerate(chain.atoms + styrene.atoms, 1):
        atom.id = identifier
    radical = next(atom for atom in chain.atoms if atom.radical_electrons)
    bond = next(bond for bond in styrene.get_all_edges() if bond.is_double()
                and styrene.is_atom_in_cycle(bond.atom1)
                and any(atom.element.number == 1 for atom in bond.atom1.edges))
    attacked, new_radical = bond.atom1, bond.atom2
    product = chain.copy(deep=True).merge(styrene.copy(deep=True))
    by_id = {atom.id: atom for atom in product.atoms}
    product.add_bond(Bond(by_id[radical.id], by_id[attacked.id], 1))
    product.get_bond(by_id[attacked.id], by_id[new_radical.id]).order = 1
    by_id[radical.id].radical_electrons = 0
    by_id[new_radical.id].radical_electrons = 1
    product.update(sort_atoms=False)
    reaction = TemplateReaction(
        reactants=[Species(molecule=[chain]), Species(molecule=[styrene])],
        products=[Species(molecule=[product])],
        kinetics=Arrhenius(A=(1, "m^3/(mol*s)")), reversible=True,
    )
    reaction.family = "R_Addition_MultipleBond"
    reaction.template = ["same_as_propagation"]
    reaction.is_forward = True
    return reaction


@pytest.mark.parametrize("units", [2, 3, 4, 6])
@pytest.mark.parametrize("reverse_participants", [False, True])
def test_match_uses_exact_graphs_and_rewrite(units, reverse_participants):
    reaction, _ = mapped_addition(units)
    if reverse_participants:
        reaction.reactants.reverse()
    reaction.template = ["unrelated_template"]
    reaction.degeneracy = 3
    record = _instance()._record(_proxy(reaction), reaction, {})
    assert matches_head_to_tail(record)
    # A correct product graph with an incorrect mapped rewrite must be rejected.
    corrupted = record.to_dict()
    edit = next(op for op in corrupted["bond_ops"] if op["action"] == "set_radical")
    edit["value"] = 2
    assert not matches_head_to_tail(corrupted)
    forward, reverse = _pair(reaction)
    assert forward.rate_source["kind"] == "kMC kinetics library"
    assert not matches_head_to_tail(reverse)
    assert forward.k_table["k"] == pytest.approx(plpsec_rate_table(GRID)["k"])
    assert forward.raw_path_degeneracy == 3  # No extra degeneracy multiplier on measured kp.
    assert all(not record.proxy_padding for record in (forward, reverse))
    assert all(
        not record.rate_witness_reactant_graphs for record in (forward, reverse)
    )


@pytest.mark.parametrize("lookalike", ["primary", "head_to_head", "ring"])
def test_nonmatching_additions_keep_both_rmg_rate_tables(lookalike):
    reaction = (_ring_addition() if lookalike == "ring" else
                mapped_addition(2, primary=lookalike == "primary", opposite=lookalike == "head_to_head")[0])
    reaction.template = ["same_as_propagation"]
    on, off = _pair(reaction), _pair(reaction, use_plpsec_library=False)
    for actual, expected in zip(on, off):
        assert not matches_head_to_tail(actual)
        assert actual.k_table == expected.k_table
        assert actual.rate_source == expected.rate_source
        assert actual.event_id == expected.event_id


@pytest.mark.parametrize("family", ["H_Abstraction", "intra_H_migration", "Disproportionation", "R_Recombination"])
def test_all_other_families_are_untouched(family):
    reaction, _ = mapped_addition(2)
    reaction.family = family
    on, off = _pair(reaction), _pair(reaction, use_plpsec_library=False)
    assert [record.to_dict() for record in on] == [record.to_dict() for record in off]


def test_units_arrhenius_detailed_balance_and_provenance():
    reaction, _ = mapped_addition(2)
    reaction.template = ["Cds-HH_Cds-CbH", "CsJ-CbCsH"]
    reaction.kinetics = Arrhenius(A=(926, "cm^3/(mol*s)"), n=2.41,
                                 Ea=(31547.4, "J/mol"), comment="sensitivity estimate")
    rule = SimpleNamespace(index=3253, label="mock sensitivity rule", rank=10, short_desc="mock rule")
    family = SimpleNamespace(
        auto_generated=False,
        extract_source_from_comments=lambda reaction: (False, ("R_Addition_MultipleBond", {
            "exact": True, "rules": [(rule, 1.0)], "training": [],
        })),
    )
    instance = _instance()
    instance.kinetics_database = SimpleNamespace(families={"R_Addition_MultipleBond": family})
    forward, reverse = instance._linked_family_pair(_proxy(reaction), reaction, {})[0]
    expected_L = [4.27e7 * math.exp(-32500 / (constants.R * t)) for t in GRID]
    assert [rate * 1000 for rate in forward.k_table["k"]] == pytest.approx(expected_L, rel=1e-12)
    assert expected_L == pytest.approx([63257.201921, 160434.678657, 322432.976207], rel=1e-9)
    assert forward.rate_units == "m^3/(mol*s)" and reverse.rate_units == "s^-1"
    assert [kf / kr for kf, kr in zip(forward.k_table["k"], reverse.k_table["k"])] == pytest.approx(
        forward.thermo_provenance["equilibrium_constant_table"]["Kc"], rel=1e-12)
    assert forward.reverse_of == reverse.event_id and reverse.reverse_of == forward.event_id
    source = forward.rate_source
    assert source["library_entry"] == load_plpsec_entry()["entry_id"]
    assert source["citation"]["doi"] == "10.1002/macp.1995.021961016"
    assert source["entry"]["measured_temperature_range_K"] == [261.15, 366.15]
    assert "600-800" in source["entry"]["extrapolation_note"]
    displaced = source["replaced_rmg_estimate"]
    assert displaced["template"] == "Cds-HH_Cds-CbH;CsJ-CbCsH"
    assert displaced["source"]["kind"] == "RMG family estimate"
    assert displaced["source"]["comment"] == "sensitivity estimate"
    assert displaced["source"]["source"] == "attached"
    assert displaced["source"]["rules"][0]["entry"]["index"] == 3253
    assert displaced["source"]["rules"][0]["weight"] == 1.0
    assert displaced["units"] == "m^3/(mol*s)"
    assert displaced["k_table"] == displaced["propagation_k_table"]
    assert displaced["k_table"]["k"] == pytest.approx([reaction.kinetics.get_rate_coefficient(t) for t in GRID])
    assert [k * 1000 / measured for k, measured in zip(displaced["k_table"]["k"], expected_L)] == pytest.approx(
        [0.130153877, 0.183633454, 0.248214469], rel=1e-8)
    assert reverse.rate_source["forward_rate_source"] == source
    assert reverse.rate_source["reference"] == "k_library / Kc"


def test_estimated_reverse_direction_is_normalized_to_propagation():
    reaction, _ = mapped_addition(2)
    reaction.kinetics = None

    class ReverseEstimatingFamily:
        def get_kinetics(self, *args, **kwargs):
            return Arrhenius(A=(2, "s^-1")), "reverse rule", None, False

    class ReferenceThermo:
        provenance = {"reference_thermo": "mock reciprocal Kc"}

        def evaluate(self, reaction, temperatures):
            assignments = tuple(
                {
                    "role": role,
                    "index": index,
                    "source_label": f"source_{role}_{index}",
                }
                for role, species_list in (
                    ("reactant", reaction.reactants),
                    ("product", reaction.products),
                )
                for index, _ in enumerate(species_list)
            )
            return ReferenceThermoResult(
                tuple([0.25] * len(temperatures)),
                12345.0,
                assignments,
            )

    database = SimpleNamespace(families={"R_Addition_MultipleBond": ReverseEstimatingFamily()})
    instance = _instance(database=database)
    instance.reference_thermo_provider = ReferenceThermo()
    records, initiating = instance._linked_family_pair(_proxy(reaction), reaction, {})
    forward, reverse = records
    assert matches_head_to_tail(forward)
    assert initiating is forward
    assert forward.k_table == plpsec_rate_table(GRID)
    assert reverse.rate_units == "s^-1"
    thermo = forward.thermo_provenance
    assert thermo["equilibrium_constant_table"]["Kc"] == [4.0] * 3
    assert thermo["reaction_enthalpy_298_J_per_mol"] == -12345.0
    assert [
        (assignment["role"], assignment["source_label"])
        for assignment in thermo["species_thermo_assignments"]
    ] == [
        ("product", "source_reactant_0"),
        ("reactant", "source_product_0"),
        ("reactant", "source_product_1"),
    ]
    displaced = forward.rate_source["replaced_rmg_estimate"]
    assert displaced["source"]["family_template_direction"] == "reverse"
    assert displaced["units"] == "s^-1"
    assert displaced["k_table"]["k"] == [2.0] * 3
    assert displaced["propagation_k_table"]["k"] == [8.0] * 3
    assert [kf / kr for kf, kr in zip(forward.k_table["k"], reverse.k_table["k"])] == pytest.approx([4.0] * 3)


@pytest.mark.parametrize("switch", ["constructor", "environment"])
def test_disable_switch_restores_rmg_estimate(monkeypatch, switch):
    reaction, _ = mapped_addition(2)
    monkeypatch.setenv("RMG_KMC_PLPSEC_LIBRARY", "0" if switch == "environment" else "1")
    instance = _instance(use_plpsec_library=False if switch == "constructor" else None)
    forward, reverse = instance._linked_family_pair(_proxy(reaction), reaction, {})[0]
    assert forward.rate_source["kind"] == "RMG family estimate"
    assert forward.k_table["k"] == [1.0] * 3
    assert reverse.rate_source["reference"] == "k_family / Kc"
    assert instance.compile()["provenance"]["kinetics_libraries"]["styrene_plpsec"]["enabled"] is False


def test_constructor_on_overrides_environment_off(monkeypatch):
    monkeypatch.setenv("RMG_KMC_PLPSEC_LIBRARY", "0")
    assert _instance(use_plpsec_library=True).use_plpsec_library


def test_invalid_environment_switch_fails(monkeypatch):
    monkeypatch.setenv("RMG_KMC_PLPSEC_LIBRARY", "no")
    with pytest.raises(ValueError, match="must be 0 or 1"):
        _instance()


def assert_compiled_plpsec_library(artifact):
    """Independent phase-2b census: exactly the two declared ordinary anchors.

    Expected graphs and numerical rates are spelled out independently of the
    production matcher/loader. This is also called by the existing default
    benzylic-anchor phase-2b acceptance check, never a new compile fixture.
    """
    block = artifact["provenance"]["kinetics_libraries"]["styrene_plpsec"]
    assert block["enabled"] is True
    assert block["entry_sha256"] == compiler.sha256_json(block["entry"])
    expected_ids = set()
    by_id = {record["event_id"]: record for record in artifact["records"]}
    for units in (2, 4):
        reactants = compiler._graph_adjacencies([
            Molecule(smiles="CC(c1ccccc1)" * (units - 1) + "C[CH](c1ccccc1)"),
            Molecule(smiles="C=Cc1ccccc1"),
        ])
        products = compiler._graph_adjacencies([
            Molecule(smiles="CC(c1ccccc1)" * units + "C[CH](c1ccccc1)"),
        ])
        matches = [record for record in artifact["records"]
                   if record["family"] == "R_Addition_MultipleBond"
                   and compiler._graph_lists_isomorphic(record["reactant_graphs"], reactants)
                   and compiler._graph_lists_isomorphic(record["product_graphs"], products)]
        assert len(matches) == 1, "phase-2b compile required: ordinary head-to-tail record absent/duplicated"
        forward = matches[0]
        expected_ids.add(forward["event_id"])
        source = forward["rate_source"]
        assert source["kind"] == "kMC kinetics library"
        assert source["entry"] == block["entry"]
        assert source["library_entry"] == "styrene_head_to_tail_propagation"
        assert source["citation"]["doi"] == "10.1002/macp.1995.021961016"
        displaced = source["replaced_rmg_estimate"]
        assert displaced["template"] == forward["template"]
        assert displaced["source"]["kind"] == "RMG family estimate"
        assert len(displaced["k_table"]["k"]) == len(forward["k_table"]["T"])
        assert displaced["propagation_k_table"]["T"] == forward["k_table"]["T"]
        expected_rates = [4.27e4 * math.exp(-32500 / (constants.R * t)) for t in forward["k_table"]["T"]]
        assert forward["k_table"]["k"] == pytest.approx(expected_rates, rel=1e-12)
        reverse = by_id[forward["reverse_of"]]
        assert reverse["reverse_of"] == forward["event_id"]
        assert reverse["rate_source"]["kind"] == "reference-thermo reverse"
        assert reverse["rate_source"]["forward_rate_source"] == source
        assert forward["rate_units"] == "m^3/(mol*s)" and reverse["rate_units"] == "s^-1"
        assert [kf / kr for kf, kr in zip(forward["k_table"]["k"], reverse["k_table"]["k"])] == pytest.approx(
            forward["thermo_provenance"]["equilibrium_constant_table"]["Kc"], rel=1e-12)
    observed_ids = {record["event_id"] for record in artifact["records"]
                    if record["rate_source"].get("kind") == "kMC kinetics library"}
    assert observed_ids == expected_ids, "library matched-record count/scope must be exactly the two head-to-tail anchors"
    return len(observed_ids)


def test_postcompile_census_and_rate_check_accept_small_inventory():
    artifact = small_inventory()
    assert assert_compiled_plpsec_library(artifact) == 2
    assert compiler.compiler_source_hash() == hashlib.sha256(b"".join(
        (Path(compiler.__file__).parent / name).read_bytes()
        for name in (
            "compiler.py",
            "reference_thermo.py",
            "event_record.py",
            "atom_map.py",
            "kinetics_library.py",
            "database_provenance.py",
            "proxy_padding.py",
            "barrier_e0.py",
        )
    )).hexdigest()


@pytest.mark.parametrize("mutation", ["rate", "extra_match"])
def test_postcompile_check_rejects_wrong_rate_or_scope(mutation):
    artifact = small_inventory()
    direct = [record for record in artifact["records"] if record["rate_source"]["kind"] == "kMC kinetics library"]
    if mutation == "rate":
        direct[0]["k_table"]["k"][0] *= 2
    else:
        other = next(record for record in artifact["records"] if record["arity"] == 2 and record not in direct)
        other["rate_source"]["kind"] = "kMC kinetics library"
    with pytest.raises(AssertionError):
        assert_compiled_plpsec_library(artifact)


def test_library_records_are_excluded_while_unvalidated_estimates_are_not_paddable(
    monkeypatch,
):
    monkeypatch.setenv("RMG_KMC_PLPSEC_LIBRARY", "1")
    on = small_inventory(proxy_padding_distance=3)
    monkeypatch.setenv("RMG_KMC_PLPSEC_LIBRARY", "0")
    off = small_inventory(proxy_padding_distance=3)
    on_records = [record for record in on["records"] if matches_head_to_tail(record)]
    off_records = [record for record in off["records"] if matches_head_to_tail(record)]
    assert on_records and off_records
    assert {
        record["proxy_padding"]["status"] for record in on_records
    } == {"excluded_unchanged"}
    assert {
        record["proxy_padding"]["status"] for record in off_records
    } == {"not_paddable"}
    on_anchor = next(
        record
        for record in on["records"]
        if record["event_id"] == on["ps_ceiling_anchor_event_id"]
    )
    off_anchor = next(
        record
        for record in off["records"]
        if record["event_id"] == off["ps_ceiling_anchor_event_id"]
    )
    assert on_anchor["k_table"]["k"] != off_anchor["k_table"]["k"]


def test_padding_off_identity_covers_library_archived_and_ordinary_records(
    monkeypatch,
):
    proxies = compiler.ps_proxy_set()
    cache = {proxy.site_type: [] for proxy in proxies}
    for size, site, primary in (
        (2, "benzylic_end_radical+styrene", False),
        (4, "benzylic_end_radical+styrene@5", False),
        (2, "end_radical+styrene", True),
        (5, "end_radical+styrene@5", True),
    ):
        cache[site] = [mapped_addition(size, primary=primary)[0]]
    cache["benzylic_end_radical+styrene"].append(
        mapped_addition(2, opposite=True)[0]
    )
    cache["end_radical@5"] = [mapped_addition(4, primary=True)[0]]
    monkeypatch.setattr(compiler, "_LOADED_COMPILER_HASH", "8a" * 32)
    monkeypatch.setattr(compiler, "compiler_source_hash", lambda: "5c" * 32)

    artifact = compiler.EventSetCompiler(
        None,
        proxies,
        ["R_Addition_MultipleBond", "R_Recombination"],
        reaction_cache=cache,
        temperature_grid=[500, 550, 600],
        rmgpy_sha="1e51f1f8f33b81ce66baf3c8b52b84f15ab9f974",
        rmg_database_sha="cd86d4e1c187a132109e16cd86f624ed9fb217df",
        reference_thermo_provider=SmallReferenceThermo(),
        ceiling_monomer_concentration_mol_m3=1,
        use_plpsec_library=True,
    ).compile()
    records = artifact["records"]
    library = [
        record
        for record in records
        if record["rate_source"].get("kind") == "kMC kinetics library"
    ]
    archived = [
        record for record in records if record.get("inventory_class") == "R1:J_ring"
    ]
    ordinary = [
        record for record in records if record not in library and record not in archived
    ]

    assert len(library) == 2
    assert sorted(
        (
            record["junction_ops"][0]["junction_kind"],
            record["junction_ops"][0]["action"],
        )
        for record in archived
    ) == [
        ("J_ortho_S6", "create"),
        ("J_ortho_S6", "dissociate"),
        ("J_ortho_S7", "create"),
        ("J_ortho_S7", "dissociate"),
        ("J_para", "create"),
        ("J_para", "dissociate"),
    ]
    assert ordinary
    assert all("proxy_padding" not in record for record in records)
    assert hashlib.sha256(compiler.canonical_json_bytes(artifact)).hexdigest() == (
        "82168c4f90b8179e7aa49bbf4f385e42cee409f23143594006406db96e60fa75"
    )
