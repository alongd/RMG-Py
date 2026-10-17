"""Tests for rmgpy/kmc/declared_chain.py (I-159): every chain witness from a
polymer declaration alone, with parity to the legacy PS and alkane providers.

Database-free: no RMG database, no reaction generation."""

import json
import os
import random
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from rmgpy.kmc import boundaries as b
from rmgpy.kmc import declared_chain as dc
from rmgpy.molecule.molecule import Molecule
from rmgpy.polymer_spec import LinearHomopolymerSpec

import declared_parity_support as sup
from meltStateTest import ps_mapping

TEST_DIR = Path(__file__).resolve().parent
ROOT = TEST_DIR.parents[2]
E = dc.EditOp
CEILING = 120


# --------------------------------------------------------------------------
# declarations (test side; no polymer is named in the product code)
# --------------------------------------------------------------------------
def _atom(key, element="C", hydrogens=0, lone_pairs=0):
    atom = {"key": key, "element": element, "hydrogens": hydrogens,
            "electronic_state": {"radical_electrons": 0}}
    if lone_pairs:
        atom["lone_pairs"] = lone_pairs
    return atom


def _bond(x, y, order=1):
    return {"atoms": [x, y], "order": order}


def _declaration(atoms, bonds, head, tail, elements=("C", "H"), lone_pairs=(0,), order=1):
    mapping = ps_mapping()
    mapping["label"] = "test declaration"
    mapping["repeat_unit"] = {"atoms": atoms, "bonds": bonds,
                              "ports": {"head": {"atom": head, "bond_order": order},
                                        "tail": {"atom": tail, "bond_order": order}},
                              "electronic_state": {"multiplicity": 1}}
    mapping["small_species"] = {}
    mapping["initial_mobile_counts"] = {}
    bounds = mapping["locality_policy"]["bounds"]
    bounds["elements"] = list(elements)
    bounds["lone_pairs"] = list(lone_pairs)
    return mapping


def pp_mapping():
    return _declaration([_atom("c1", hydrogens=2), _atom("c2", hydrogens=1),
                         _atom("c3", hydrogens=3)],
                        [_bond("c1", "c2"), _bond("c2", "c3")], "c1", "c2")


def pmma_mapping():
    """[*]CC([*])(C)C(=O)OC, ester and carbonyl O with two lone pairs."""
    return _declaration(
        [_atom("c1", hydrogens=2), _atom("c2"), _atom("c3", hydrogens=3), _atom("c4"),
         _atom("o1", "O", lone_pairs=2), _atom("o2", "O", lone_pairs=2),
         _atom("c5", hydrogens=3)],
        [_bond("c1", "c2"), _bond("c2", "c3"), _bond("c2", "c4"), _bond("c4", "o1", 2),
         _bond("c4", "o2"), _bond("o2", "c5")], "c1", "c2",
        elements=("C", "H", "O"), lone_pairs=(0, 2))


def para_mapping():
    """CH2-C6H4 joined at the para position."""
    ring = [_atom(f"r{i}", hydrogens=0 if i in (1, 4) else 1) for i in range(1, 7)]
    bonds = [_bond(f"r{i}", f"r{i % 6 + 1}", 1.5) for i in range(1, 7)]
    return _declaration([_atom("c", hydrogens=2)] + ring, [_bond("c", "r1")] + bonds, "c", "r4")


def polyacetylene_mapping():
    return _declaration([_atom("c1", hydrogens=1), _atom("c2", hydrogens=1)],
                        [_bond("c1", "c2", 2)], "c1", "c2")


def ppv_mapping():
    ring = [_atom(f"r{i}", hydrogens=0 if i in (1, 4) else 1) for i in range(1, 7)]
    bonds = [_bond(f"r{i}", f"r{i % 6 + 1}", 1.5) for i in range(1, 7)]
    return _declaration(ring + [_atom("v1", hydrogens=1), _atom("v2", hydrogens=1)],
                        bonds + [_bond("r4", "v1"), _bond("v1", "v2", 2)], "r1", "v2")


def spec_of(mapping):
    return LinearHomopolymerSpec.from_mapping(mapping)


@pytest.fixture(scope="module")
def ps():
    return sup.ps_context()


@pytest.fixture(scope="module")
def pe():
    return dc.DeclaredChainProvider(sup.pe_declaration())


@pytest.fixture(scope="module")
def pp():
    return dc.DeclaredChainProvider(spec_of(pp_mapping()))


@pytest.fixture(scope="module")
def pmma():
    return dc.DeclaredChainProvider(spec_of(pmma_mapping()))


def rows_of(build):
    return dict(b.rows_of_molecule(build.molecule))


def heavy_formula(build):
    counts = {}
    for row in rows_of(build).values():
        counts[row[0]] = counts.get(row[0], 0) + 1
    return counts


def at(build, unit, key, provider):
    for a, u, k in build.metadata.coordinates:
        if u == unit and provider.key_of(k)[0] == key:
            return a
    raise KeyError((unit, key))


# --------------------------------------------------------------------------
# registration and the PS position map
# --------------------------------------------------------------------------
def test_provider_is_registered_under_the_declaration_digest(ps):
    provider = ps.provider
    assert provider.name.startswith("declared:") and len(provider.name) == len("declared:") + 16
    assert b.PROVIDERS[provider.name] is provider
    build = provider.build("x", dc.ChainSpec(2))
    assert build.metadata.provider == provider.name


def test_ps_position_map_is_derived_and_keeps_roles(ps):
    geometry = ps.provider.geometry
    assert ps.keys[0] == geometry.head_key and ps.keys[1] == geometry.tail_key
    ring = set(ps.keys[2:])
    assert len(ring) == 6 and not ring & {ps.keys[0], ps.keys[1]}
    # para (legacy position 5) is the ring atom farthest from the ipso
    assert geometry.backbone_offsets.get(ps.keys[5]) is None
    assert sorted(ps.pi) == list(range(8))


# --------------------------------------------------------------------------
# PS parity
# --------------------------------------------------------------------------
CATALOG = sup.ps_catalog()


@pytest.mark.parametrize("group", sorted(CATALOG))
def test_ps_parity_with_the_legacy_provider(ps, group):
    mismatches, refused = sup.ps_parity(ps, CATALOG[group])
    assert not mismatches, mismatches[:5]
    assert len(CATALOG[group]) - len(refused) > 0


def test_ps_catalog_covers_every_spec_field():
    specs = [s for group in CATALOG.values() for s in group]
    fields = ("radical_units", "radical_end", "benzylic_end", "ring_radical_unit", "monomer",
              "head_end", "tail_end", "head_shift", "tail_shift", "secondary_radical_units",
              "pendant_loss_units")
    for field in fields:
        if field == "ring_radical_unit":
            assert any(s.ring_radical_unit is not None for s in specs)
        else:
            assert any(getattr(s, field) not in (None, False, ()) for s in specs), field


def test_ps_radical_ends_are_physical_and_cap_free(ps):
    legacy = b.ps_build("end_radical", 3)
    declared = sup.declared_ps_build(ps, b.ps_structure_specs(3)["end_radical"])
    head = declared.metadata.port("head")
    assert head.physicality == "physical" and head.cap_atom_ids == ()
    assert rows_of(declared)[head.atom_id][3] == 1
    assert sup.builds_match(legacy, declared, ps.coord_map)


@pytest.mark.parametrize("chain", [dc.ChainSpec(3), dc.ChainSpec(2, head="radical"),
                                   dc.ChainSpec(3, tail_phase=1, tail="declared")])
def test_caps_are_the_last_port_hydrogens_created(ps, chain):
    build = ps.provider.build("c", chain)
    r = rows_of(build)
    for port in build.metadata.ports:
        if not port.artificial:
            assert port.cap_atom_ids == ()
            continue
        port_h = [a for a in build.construction_order
                  if r[a][0] == "H" and r[a][5][0][0] == port.atom_id]
        assert port.cap_atom_ids == (port_h[-1],)


def test_ps_monomer_has_no_ports(ps):
    monomer = ps.provider.monomer("m")
    assert monomer.metadata.ports == ()
    assert heavy_formula(monomer) == {"C": 8, "H": 8}


def test_ps_fragment_parity_at_head_and_tail(ps):
    legacy = b.ps_build("pristine", 4)
    declared = ps.provider.build("p", dc.ChainSpec(4))
    for port_id in ("head", "tail"):
        lp, dp = legacy.metadata.port(port_id), declared.metadata.port(port_id)
        lc = legacy.metadata.coordinate(lp.atom_id)
        dcoord = declared.metadata.coordinate(dp.atom_id)
        assert ps.coord_map(lc) == dcoord
        assert sup.fragments_match(ps.legacy.fragment(lp, lc), ps.provider.fragment(dp, dcoord),
                                   ps.coord_map), port_id


def test_ps_fragment_refuses_a_phase_cut_port(ps):
    build = ps.provider.build("s", dc.ChainSpec(3, head_phase=ps.phase))
    port = build.metadata.port("head")
    with pytest.raises(b.ExtensionRefused, match="phase cut"):
        ps.provider.fragment(port, build.metadata.coordinate(port.atom_id))


DP4 = ("pristine", "interior_radical", "junction_radical")


@pytest.mark.parametrize("name", DP4)
def test_ps_extend_vector_parity(ps, name):
    spec = b.ps_structure_specs(4)[name]
    legacy = sup.cut_witness(ps.legacy.build("l", spec), (1, 1), (2, 0))
    declared = sup.cut_witness(sup.declared_ps_build(ps, spec),
                               ps.coord_map((1, 1)), ps.coord_map((2, 0)))
    assert sup.witnesses_match(legacy, declared, ps.coord_map)
    lres = {r.vector: r for r in b.enumerate_vectors(legacy, 2, ceiling=CEILING)}
    dres = {r.vector: r for r in b.enumerate_vectors(declared, 2, ceiling=CEILING)}
    assert sorted(lres) == sorted(dres) and len(lres) == 9
    for vector, lr in lres.items():
        dr = dres[vector]
        assert lr.status == dr.status, vector
        if lr.status == "ok":
            assert sup.witnesses_match(lr.witness, dr.witness, ps.coord_map), vector
            ext, _record = b.extend_vector(declared, vector, ceiling=CEILING)
            assert sup.witnesses_match(lr.witness, ext, ps.coord_map), vector
    assert b.vector_duplicates(list(lres.values())) == \
        b.vector_duplicates(list(dres.values()))


# --------------------------------------------------------------------------
# alkane parity (the PE declaration, primitive period one CH2)
# --------------------------------------------------------------------------
ALKANE_ENDS = (None, "methyl", "primary_radical")


def _alkane_key(pe):
    (key,) = [k for k in pe.geometry.keys]
    return key


@pytest.mark.parametrize("n", range(2, 26))
def test_alkane_parity(pe, n):
    key = _alkane_key(pe)
    for head in ALKANE_ENDS:
        for tail in ALKANE_ENDS:
            spec = ("alkane", n, head, tail)
            legacy = sup.legacy_alkane_build(spec)
            declared = pe.build("d", sup.legacy_alkane_to_chain(spec, key))
            assert sup.builds_match(legacy, declared), spec
    plain = ("alkane", n)
    assert sup.builds_match(sup.legacy_alkane_build(plain),
                            pe.build("d", sup.legacy_alkane_to_chain(plain, key)))


@pytest.mark.parametrize("n", range(2, 26))
def test_alkane_interior_radical_parity(pe, n):
    key = _alkane_key(pe)
    sites = sorted({0, 1, n // 2, n - 2, n - 1} & set(range(n)))
    for r in sites:
        for head, tail in ((None, None), ("methyl", None), (None, "primary_radical"),
                           ("primary_radical", "methyl")):
            spec = ("alkane", n, head, tail, r)
            legacy = sup.legacy_alkane_build(spec)
            declared = pe.build("d", sup.legacy_alkane_to_chain(spec, key))
            assert sup.builds_match(legacy, declared), spec


def test_c25r11_interior_radical(pe):
    spec = ("alkane", 25, None, None, 11)
    declared = pe.build("c25r11", sup.legacy_alkane_to_chain(spec, _alkane_key(pe)))
    assert sup.builds_match(sup.legacy_alkane_build(spec), declared)
    assert heavy_formula(declared) == {"C": 25, "H": 51}


def test_alkane_fragment_parity(pe):
    legacy = b.PROVIDERS["linear_alkane"].build(6)
    declared = pe.build("d", dc.ChainSpec(6))
    for port_id in ("head", "tail"):
        lp, dp = legacy.metadata.port(port_id), declared.metadata.port(port_id)
        coord = legacy.metadata.coordinate(lp.atom_id)
        assert declared.metadata.coordinate(dp.atom_id) == coord
        assert sup.fragments_match(b.PROVIDERS["linear_alkane"].fragment(lp, coord),
                                   pe.fragment(dp, coord)), port_id


def test_c16_extended_by_one_period_is_c17(pe):
    c16 = sup.identity_witness(pe.build("c16", dc.ChainSpec(16)))
    c17 = sup.identity_witness(pe.build("c17", dc.ChainSpec(17)))
    port = c16.port_keys().index((0, "tail"))
    vector = [0] * len(c16.port_keys())
    vector[port] = 1
    extended, _record = b.extend_vector(c16, tuple(vector), ceiling=CEILING)
    assert sup.witnesses_match(c17, extended)
    legacy = sup.identity_witness(b.PROVIDERS["linear_alkane"].build(16))
    legacy_ext, _ = b.extend_vector(legacy, tuple(vector), ceiling=CEILING)
    assert sup.witnesses_match(legacy_ext, extended)


# --------------------------------------------------------------------------
# geometry
# --------------------------------------------------------------------------
def test_geometry_ps(ps):
    g = ps.provider.geometry
    assert (g.n_bb, g.d_rep, g.period) == (2, 2, 2)
    assert len(g.phases) == 1 and g.bridge_ok and g.pi_ok
    assert g.heavy_classes == 6
    assert not g.reflection_symmetric
    # ortho pair and meta pair merge under the ring flip
    sizes = sorted(len(positions) for _key, positions in g.atom_classes)
    assert sizes == [1, 1, 1, 1, 2, 2]


def test_geometry_pe(pe):
    g = pe.geometry
    assert (g.n_bb, g.d_rep, g.period) == (2, 2, 1)
    assert g.reflection_symmetric
    assert len(g.keys) == 1 and g.heavy_classes == 1 and g.phases == ()


def test_geometry_pp(pp):
    g = pp.geometry
    assert (g.n_bb, g.d_rep, g.period) == (2, 2, 2)
    assert g.heavy_classes == 3 and len(g.phases) == 1 and not g.reflection_symmetric
    junction = pp.build("hh", dc.ChainSpec(2, orientation=(1, -1), head="declared",
                                           tail="declared"))
    assert junction.molecule.is_isomorphic(Molecule(smiles="CC(C)C(C)C"))
    c2 = pp.geometry.tail_key
    a, c = at(junction, 0, c2, pp), at(junction, 1, c2, pp)
    assert any(n == c for n, _o in rows_of(junction)[a][5])
    straight = pp.build("ht", dc.ChainSpec(2, head="declared", tail="declared"))
    assert straight.molecule.is_isomorphic(Molecule(smiles="CCCC(C)C"))


def test_geometry_pmma(pmma):
    g = pmma.geometry
    assert (g.n_bb, g.d_rep, g.period) == (2, 2, 2)
    assert g.unit_formula == {"C": 5, "H": 8, "O": 2}
    assert g.largest_pi == 2 and g.pi_ok and g.bridge_ok
    dp3 = pmma.build("dp3", dc.ChainSpec(3, head="declared", tail="declared"))
    molecule = dp3.molecule.copy(deep=True)
    molecule.update()
    assert dp3.molecule.is_isomorphic(
        Molecule(smiles="CC(C)(C(=O)OC)CC(C)(C(=O)OC)CC(C)(C(=O)OC)"))
    oxygens = [row for row in rows_of(dp3).values() if row[0] == "O"]
    assert len(oxygens) == 6 and all(row[4] == 2 for row in oxygens)
    open_ = pmma.build("o", dc.ChainSpec(3))
    port = open_.metadata.port("tail")
    fragment = pmma.fragment(port, open_.metadata.coordinate(port.atom_id))
    assert sorted(atom[4] for atom in fragment.atoms if atom[0] == "O") == [2, 2]


def test_geometry_para_ring_backbone():
    provider = dc.DeclaredChainProvider(spec_of(para_mapping()))
    g = provider.geometry
    # n_bb = |melt_state._backbone_keys| (5 path atoms); d_rep = c..c' = 5
    assert (g.n_bb, g.d_rep, g.period) == (5, 5, 5)
    assert len(g.phases) == 1
    p = g.positions
    assert provider.backbone_index(0, p["r2"]) == provider.backbone_index(0, p["r6"]) == 2
    assert provider.backbone_index(0, p["r3"]) == provider.backbone_index(0, p["r5"]) == 3
    assert provider.backbone_index(1, p["c"]) == 5
    assert provider.backbone_index(1, p["r4"]) == 9


@pytest.mark.parametrize("mapping", [polyacetylene_mapping, ppv_mapping])
def test_conjugated_backbones_are_refused_by_rule_2b(mapping):
    spec = spec_of(mapping())
    with pytest.raises(dc.DeclarationRefused) as info:
        dc.ChainGeometry.of(spec)
    assert info.value.rule == "2b"
    assert info.value.component and "crosses" in str(info.value)


def _ladder_unit(spiro=False):
    """Two joins per residue pair.  Ladder: a-b rung, joins c->a' and d->b'.
    Spiro: both joins leave c (c->a', c->b'), closing the ring a'-c'-d'-b'."""
    from rmgpy.polymer_spec import GraphDeclaration
    if spiro:
        atoms = [_atom("a", hydrogens=2), _atom("b", hydrogens=2), _atom("c"),
                 _atom("d", hydrogens=2)]
        bonds = [_bond("a", "c"), _bond("b", "d"), _bond("c", "d")]
        tail_2 = "c"
    else:
        atoms = [_atom("a", hydrogens=1), _atom("b", hydrogens=1), _atom("c", hydrogens=2),
                 _atom("d", hydrogens=2)]
        bonds = [_bond("a", "b"), _bond("a", "c"), _bond("b", "d")]
        tail_2 = "d"
    ports = {"head": {"atom": "a", "bond_order": 1}, "tail": {"atom": "c", "bond_order": 1},
             "head:2": {"atom": "b", "bond_order": 1},
             "tail:2": {"atom": tail_2, "bond_order": 1}}
    return GraphDeclaration.from_mapping({"atoms": atoms, "bonds": bonds, "ports": ports,
                                          "electronic_state": {"multiplicity": 1}})


@pytest.mark.parametrize("spiro", [False, True])
def test_ladder_and_spiro_units_are_refused_by_rule_2(spiro):
    with pytest.raises(dc.DeclarationRefused,
                       match="inter-residue bond 0:c-1:a lies on a cycle") as info:
        dc.ChainGeometry.of_unit(_ladder_unit(spiro))
    assert info.value.rule == "2" and info.value.component


def test_single_join_units_pass_rule_2(pp):
    assert pp.geometry.bridge_ok


# --------------------------------------------------------------------------
# edits
# --------------------------------------------------------------------------
def test_pp_edits(pp):
    c1, c2, c3 = "c1", "c2", "c3"
    tertiary = pp.build("t", dc.ChainSpec(3, edits=(E("remove_h", 1, c2),
                                                    E("set_radical", 1, c2, 1))))
    row = rows_of(tertiary)[at(tertiary, 1, c2, pp)]
    assert row[3] == 1 and sum(1 for n, _o in row[5] if rows_of(tertiary)[n][0] == "H") == 0
    methyl = pp.build("m", dc.ChainSpec(3, edits=(E("remove_h", 1, c3),
                                                  E("set_radical", 1, c3, 1))))
    assert rows_of(methyl)[at(methyl, 1, c3, pp)][3] == 1
    # add_h: a radical tail capped back is the declared tail
    capped = pp.build("c", dc.ChainSpec(3, tail="radical", edits=(
        E("add_h", 2, c2), E("set_radical", 2, c2, 0))))
    declared = pp.build("d", dc.ChainSpec(3, tail="declared"))
    tail = capped.metadata.port("tail")
    assert tail.physicality == "physical"
    assert capped.molecule.is_isomorphic(declared.molecule)
    # an end alkene on the PP backbone (set_bond_order across units)
    alkene = pp.build("v", dc.ChainSpec(2, head="declared", tail="declared", edits=(
        E("remove_h", 0, c2), E("remove_h", 1, c1),
        E("set_bond_order", 0, c2, 2, other=(1, c1)))))
    assert alkene.molecule.is_isomorphic(Molecule(smiles="CC(C)=CCC"))


def test_pmma_edits(pmma):
    end = pmma.build("me", dc.ChainSpec(3, head="declared", tail="radical", edits=(
        E("remove_h", 2, "c1"), E("set_bond_order", 2, "c1", 2, other=(2, "c2")),
        E("set_radical", 2, "c2", 0))))
    assert end.molecule.is_isomorphic(
        Molecule(smiles="CC(C)(C(=O)OC)CC(C)(C(=O)OC)C=C(C)C(=O)OC"))
    loss = pmma.build("acyl", dc.ChainSpec(3, head="declared", tail="declared", edits=(
        E("delete_atom", 1, "o2"), E("delete_atom", 1, "c5"), E("set_radical", 1, "c4", 1))))
    assert loss.molecule.is_isomorphic(
        Molecule(smiles="CC(C)(C(=O)OC)CC(C)([C]=O)CC(C)(C(=O)OC)"))
    assert heavy_formula(loss)["O"] == 5
    positions = {(u, k) for _a, u, k in loss.metadata.coordinates}
    p = pmma.geometry.positions
    assert (1, p["o2"]) not in positions and (1, p["c5"]) not in positions


@pytest.mark.parametrize("edit, message", [
    (E("set_radical", 3, "c2", 1), "unit 3 is outside the chain"),
    (E("set_radical", -1, "c2", 1), "unit -1 is outside the chain"),
    (E("remove_h", 1, "zz"), "key 'zz' is not an atom of the period"),
    (E("set_bond_order", 1, "c1", 2, other=(1, "c3")), "are not bonded"),
    (E("set_bond_order", 1, "c1", 2, other=(5, "c2")), "unit 5 is outside the chain"),
    (E("set_radical", 1, "c2", 1), "its valence needs"),
    (E("remove_h", 1, "c2", 2), "removable hydrogens"),
    (E("delete_atom", 0, "c1"), "carries a port"),
    (E("delete_atom", 2, "c2"), "carries a port"),
])
def test_edit_refusals(pp, edit, message):
    """Each refusal is pinned by class and message: a guard that is removed
    must not be silently replaced by a later check."""
    with pytest.raises(dc.EditError, match=message):
        pp.build("x", dc.ChainSpec(3, edits=(edit,)))


def test_unbonded_bond_edit_is_refused_even_when_valence_would_fit(pp):
    # c3 of unit 0 and c3 of unit 1 are both methyls: removing one H on each
    # and bonding them would satisfy valence, so only the bonded check refuses
    with pytest.raises(dc.EditError, match="are not bonded"):
        pp.build("x", dc.ChainSpec(2, edits=(
            E("remove_h", 0, "c3"), E("remove_h", 1, "c3"),
            E("set_bond_order", 0, "c3", 1, other=(1, "c3")))))


def test_valence_is_checked_after_all_edits(pp):
    # remove_h then set_radical is legal only as a pair
    with pytest.raises(dc.EditError):
        pp.build("x", dc.ChainSpec(3, edits=(E("remove_h", 1, "c2"),)))
    pp.build("ok", dc.ChainSpec(3, edits=(E("remove_h", 1, "c2"), E("set_radical", 1, "c2", 1))))


def test_edit_on_trimmed_atom_is_refused(ps):
    with pytest.raises(dc.EditError):
        ps.provider.build("x", dc.ChainSpec(2, head_phase=ps.phase, edits=(
            E("remove_h", 0, ps.keys[0]), E("set_radical", 0, ps.keys[0], 1))))


@pytest.mark.parametrize("bad", [
    {"units": 2, "colour": "red"},
    {"units": 0},
    {"units": 2, "head": "bogus"},
    {"units": 2, "orientation": (1,)},
    {"units": 2, "orientation": (1, 2)},
    {"units": 2, "head_phase": 7},
    {"units": 2, "edits": ({"kind": "explode", "unit": 0, "key": "c1"},)},
])
def test_chain_spec_refusals(pp, bad):
    with pytest.raises(dc.ChainSpecError):
        pp.build("x", dc.ChainSpec.from_mapping(bad))


def test_heavy_declared_end_groups_carry_end_coordinates():
    mapping = pp_mapping()
    methyl = {"atoms": [_atom("m", hydrogens=3)], "bonds": [],
              "ports": {"cap": {"atom": "m", "bond_order": 1}},
              "electronic_state": {"multiplicity": 1}}
    hydroxyl = {"atoms": [_atom("o", "O", hydrogens=1, lone_pairs=2)], "bonds": [],
                "ports": {"cap": {"atom": "o", "bond_order": 1}},
                "electronic_state": {"multiplicity": 1}}
    mapping["end_groups"] = {"head": methyl, "tail": hydroxyl}
    mapping["locality_policy"]["bounds"]["elements"] = ["C", "H", "O"]
    mapping["locality_policy"]["bounds"]["lone_pairs"] = [0, 2]
    provider = dc.DeclaredChainProvider(spec_of(mapping))
    build = provider.build("e", dc.ChainSpec(2, head="declared", tail="declared"))
    assert build.molecule.is_isomorphic(Molecule(smiles="CCC(C)CC(C)O"))
    coords = {a: (u, k) for a, u, k in build.metadata.coordinates}
    r = rows_of(build)
    assert coords[0] == (-1, 0) and r[0][0] == "C"          # head group first
    (oxygen,) = [a for a, row in r.items() if row[0] == "O"]
    assert coords[oxygen] == (2, 0)
    assert all(p.physicality == "physical" and not p.cap_atom_ids
               for p in build.metadata.ports)


def test_explicit_hydrogen_atoms_build_like_hydrogen_counts(pp):
    mapping = pp_mapping()
    atoms = mapping["repeat_unit"]["atoms"]
    atoms[2]["hydrogens"] = 0
    atoms += [_atom(f"h{i}", "H") for i in range(3)]
    mapping["repeat_unit"]["bonds"] += [_bond("c3", f"h{i}") for i in range(3)]
    explicit = dc.DeclaredChainProvider(spec_of(mapping))
    for chain in (dc.ChainSpec(3), dc.ChainSpec(2, head="declared", tail="radical")):
        assert sup.builds_match(pp.build("a", chain), explicit.build("b", chain))
    open_ = explicit.build("o", dc.ChainSpec(2))
    port = open_.metadata.port("tail")
    assert sup.fragments_match(
        pp.fragment(pp.build("o", dc.ChainSpec(2)).metadata.port("tail"), (1, 1)),
        explicit.fragment(port, open_.metadata.coordinate(port.atom_id)),
        lambda c: c)


def test_declared_end_edit_equals_the_declared_build(pp):
    open_ = sup.identity_witness(pp.build("o", dc.ChainSpec(3)))
    built, record = b.reconstruct_end(open_, (0, "tail"), "declared")
    assert built.port((0, "tail")).physicality == "physical"
    direct = sup.identity_witness(pp.build("d", dc.ChainSpec(3, tail="declared")))
    assert sup.witnesses_match(direct, built)
    assert b.end_chemistries(open_, (0, "head")) == ("declared",)


# --------------------------------------------------------------------------
# determinism
# --------------------------------------------------------------------------
def _shuffled(spec, seed):
    rng = random.Random(seed)
    repeat = spec.repeat_unit
    atoms, bonds = list(repeat.atoms), list(repeat.bonds)
    rng.shuffle(atoms)
    rng.shuffle(bonds)
    if [a.key for a in atoms] == sorted(a.key for a in atoms):
        atoms.reverse()
    return replace(spec, repeat_unit=replace(repeat, atoms=tuple(atoms), bonds=tuple(bonds)))


@pytest.mark.parametrize("make", [sup.ps_declaration, sup.pe_declaration,
                                  lambda: spec_of(pmma_mapping())])
def test_shuffled_declaration_gives_identical_builds(make):
    spec = make()
    reference = sup.fingerprint(spec)
    for seed in (1, 2):
        shuffled = _shuffled(spec, seed)
        assert [a.key for a in shuffled.repeat_unit.atoms] != \
            [a.key for a in spec.repeat_unit.atoms]
        assert sup.fingerprint(shuffled) == reference
    # the mapping-level shuffle is normalized by the declaration parser
    mapping = json.loads(spec.normalized_json)
    mapping["repeat_unit"]["atoms"].reverse()
    assert sup.fingerprint(LinearHomopolymerSpec.from_mapping(mapping)) == reference


def test_identical_across_hash_seeds():
    code = ("import sys; sys.path.insert(0, %r); import declared_parity_support as s; "
            "print(s.fingerprint())" % str(TEST_DIR))
    out = []
    for seed in ("0", "1"):
        env = dict(os.environ, PYTHONHASHSEED=seed,
                   MPLCONFIGDIR=os.environ.get("TMPDIR", "/tmp"))
        proc = subprocess.run([sys.executable, "-c", code], cwd=str(ROOT), env=env,
                              capture_output=True, text=True, timeout=600)
        assert proc.returncode == 0, proc.stderr[-2000:]
        out.append(proc.stdout.strip().splitlines()[-1])
    assert out[0] == out[1] == sup.fingerprint()


# --------------------------------------------------------------------------
# rework 1 (review r101)
# --------------------------------------------------------------------------
def _radical_atom(key, hydrogens):
    atom = _atom(key, hydrogens=hydrogens)
    atom["electronic_state"] = {"radical_electrons": 1}
    return atom


def _singlet_diradical_monomer():
    from melt_engine_support import pe_mapping
    singlet = {"atoms": [_radical_atom("a", 2), _radical_atom("b", 2)],
               "bonds": [_bond("a", "b")], "electronic_state": {"multiplicity": 1}}
    return pe_mapping({10: 1}, transport=False, monomer=singlet)


def pendant_diradical_mapping(multiplicity=3):
    """c1H2-c2(c3H2.)(c4H2.): two pendant radicals per unit."""
    mapping = _declaration([_atom("c1", hydrogens=2), _atom("c2"), _radical_atom("c3", 2),
                            _radical_atom("c4", 2)],
                           [_bond("c1", "c2"), _bond("c2", "c3"), _bond("c2", "c4")], "c1", "c2")
    mapping["repeat_unit"]["electronic_state"] = {"multiplicity": multiplicity}
    return mapping


def _pendant_radical_singlet():
    """-CH2-CH(CH2.)- declared as a singlet."""
    mapping = _declaration([_atom("c1", hydrogens=2), _atom("c2", hydrogens=1),
                            _radical_atom("c3", 2)],
                           [_bond("c1", "c2"), _bond("c2", "c3")], "c1", "c2")
    mapping["repeat_unit"]["electronic_state"] = {"multiplicity": 1}
    return mapping


def _radical_end_group():
    mapping = pp_mapping()
    mapping["end_groups"]["head"] = {
        "atoms": [_radical_atom("m", 2)], "bonds": [],
        "ports": {"cap": {"atom": "m", "bond_order": 1}},
        "electronic_state": {"multiplicity": 2}}
    return mapping


@pytest.mark.parametrize("make, where", [
    (_singlet_diradical_monomer, "monomer_product"),
    (pendant_diradical_mapping, "repeat unit"),
    (_pendant_radical_singlet, "repeat unit"),
    (_radical_end_group, "head end group"),
])
def test_open_shell_declarations_are_refused(make, where):
    """S1 supports closed-shell declarations only (no declared radical
    electrons, multiplicity 1): witness rows carry no multiplicity."""
    spec = spec_of(make())
    with pytest.raises(dc.DeclarationRefused, match=f"{where} is open-shell") as info:
        dc.DeclaredChainProvider(spec)
    assert info.value.rule == "closed_shell"


def test_end_and_edit_radicals_keep_rmg_default_multiplicity(pp):
    closed = pp.build("c", dc.ChainSpec(3, head="declared", tail="declared"))
    assert closed.molecule.multiplicity == 1
    three = pp.build("r", dc.ChainSpec(3, head="radical", tail="radical", edits=(
        E("remove_h", 1, "c2"), E("set_radical", 1, "c2", 1))))
    assert three.molecule.multiplicity == 4
    monomer_free = dc.DeclaredChainProvider(sup.ps_declaration()).monomer("m")
    assert monomer_free.molecule.multiplicity == 1


def _backbone_walk(build, provider):
    """backbone_index of the atoms on the built head -> tail path."""
    import networkx as nx
    r = rows_of(build)
    graph = nx.Graph()
    for a, row in r.items():
        if row[0] != "H":
            for c, _o in row[5]:
                if r[c][0] != "H":
                    graph.add_edge(a, c)
    head, tail = build.metadata.port("head").atom_id, build.metadata.port("tail").atom_id
    coords = {a: (u, k) for a, u, k in build.metadata.coordinates}
    path = nx.shortest_path(graph, head, tail) if head != tail else [head]
    return [provider.backbone_index(*coords[a]) for a in path]


@pytest.mark.parametrize("orientation", [(1, -1), (-1, 1), (1, -1, 1), (-1, -1), (1, 1, -1, -1)])
def test_backbone_index_is_consecutive_along_reversed_chains(pp, orientation):
    build = pp.build("o", dc.ChainSpec(len(orientation), orientation=orientation))
    walk = _backbone_walk(build, pp)
    assert walk == list(range(walk[0], walk[0] + len(walk))), walk
    assert walk[0] == 0


def test_backbone_index_on_reversed_ps(ps):
    build = ps.provider.build("o", dc.ChainSpec(3, orientation=(1, -1, 1)))
    walk = _backbone_walk(build, ps.provider)
    assert walk == list(range(6))


def test_extension_at_reversed_ends(pp):
    # tail of a reversed last period: the chain continues reversed
    hh = sup.identity_witness(pp.build("hh", dc.ChainSpec(2, orientation=(1, -1))))
    extended, _ = b.extend_vector(hh, (0, 1), ceiling=CEILING)
    direct = sup.identity_witness(pp.build("d", dc.ChainSpec(3, orientation=(1, -1, -1))))
    assert sup.witnesses_match(direct, extended)
    # head of a reversed first period
    tt = sup.identity_witness(pp.build("tt", dc.ChainSpec(2, orientation=(-1, 1))))
    extended, _ = b.extend_vector(tt, (1, 0), ceiling=CEILING)
    direct = sup.identity_witness(pp.build("d", dc.ChainSpec(3, orientation=(-1, -1, 1))))
    assert sup.witnesses_match(direct, extended, lambda c: (c[0] - 1, c[1]))
    # both junction kinds survive two extensions at both ends
    for orientation in ((1, -1), (-1, 1)):
        w = sup.identity_witness(pp.build("j", dc.ChainSpec(2, orientation=orientation)))
        extended, _ = b.extend_vector(w, (2, 2), ceiling=CEILING)
        assert extended.heavy_count() == 6 * 3


def pe_explicit_mapping():
    """PE with c1 carrying implicit H2 and c2 two explicit H atoms."""
    from melt_engine_support import pe_mapping
    mapping = pe_mapping({10: 1}, transport=False)
    repeat = mapping["repeat_unit"]
    repeat["atoms"] = [_atom("c1", hydrogens=2), _atom("c2"), _atom("h1", "H"), _atom("h2", "H")]
    repeat["bonds"] = [_bond("c1", "c2"), _bond("c2", "h1"), _bond("c2", "h2")]
    return mapping


def test_hydrogen_representation_does_not_change_the_period(pe):
    explicit = dc.DeclaredChainProvider(spec_of(pe_explicit_mapping()))
    assert explicit.geometry.period == 1
    c3 = explicit.build("c3", dc.ChainSpec(3))
    assert heavy_formula(c3) == {"C": 3, "H": 8}
    assert sup.builds_match(pe.build("c3", dc.ChainSpec(3)), c3)
    w = sup.identity_witness(c3)
    extended, _ = b.extend_vector(w, (0, 1), ceiling=CEILING)
    assert extended.heavy_count() == 4


def test_aromatic_representation_does_not_change_the_period():
    """CH2-C6H4-CH2-C6H4 (para), one ring declared aromatic (1.5), the other
    Kekule: the primitive period is still CH2-C6H4."""
    ring_a = [_atom(f"a{i}", hydrogens=0 if i in (1, 4) else 1) for i in range(1, 7)]
    ring_b = [_atom(f"b{i}", hydrogens=0 if i in (1, 4) else 1) for i in range(1, 7)]
    bonds = [_bond(f"a{i}", f"a{i % 6 + 1}", 1.5) for i in range(1, 7)]
    bonds += [_bond(f"b{i}", f"b{i % 6 + 1}", 2 if i % 2 else 1) for i in range(1, 7)]
    bonds += [_bond("m1", "a1"), _bond("a4", "m2"), _bond("m2", "b1")]
    mapping = _declaration([_atom("m1", hydrogens=2), _atom("m2", hydrogens=2)] + ring_a + ring_b,
                           bonds, "m1", "b4")
    provider = dc.DeclaredChainProvider(spec_of(mapping))
    assert (provider.geometry.d_rep, provider.geometry.period) == (10, 5)
    para = dc.DeclaredChainProvider(spec_of(para_mapping()))
    assert provider.build("x", dc.ChainSpec(3, head="declared", tail="declared")).molecule \
        .is_isomorphic(para.build("y", dc.ChainSpec(3, head="declared", tail="declared")).molecule)


def _cyclobutadiene_mapping(alternate=True):
    """m1-[c4 ring a, 1,2-attached, attachment edge double]-m2-[c4 ring b,
    attachment edge single (alternate) or double]: two cyclobutadiene
    patterns that differ in real bond orders (not aromatic)."""
    atoms = [_atom("m1", hydrogens=2), _atom("m2", hydrogens=2)]
    for ring in "ab":
        atoms += [_atom(f"{ring}1"), _atom(f"{ring}2"), _atom(f"{ring}3", hydrogens=1),
                  _atom(f"{ring}4", hydrogens=1)]
    bonds = [_bond("m1", "a1"), _bond("a2", "m2"), _bond("m2", "b1"),
             _bond("a1", "a2", 2), _bond("a2", "a3"), _bond("a3", "a4", 2), _bond("a4", "a1")]
    if alternate:
        bonds += [_bond("b1", "b2"), _bond("b2", "b3", 2), _bond("b3", "b4"),
                  _bond("b4", "b1", 2)]
    else:
        bonds += [_bond("b1", "b2", 2), _bond("b2", "b3"), _bond("b3", "b4", 2),
                  _bond("b4", "b1")]
    return _declaration(atoms, bonds, "m1", "b2")


def test_non_aromatic_rings_keep_their_bond_orders_in_period_detection():
    alternating = dc.DeclaredChainProvider(spec_of(_cyclobutadiene_mapping()))
    g = alternating.geometry
    assert (g.d_rep, g.period) == (6, 6)
    uniform = dc.DeclaredChainProvider(spec_of(_cyclobutadiene_mapping(alternate=False)))
    assert (uniform.geometry.d_rep, uniform.geometry.period) == (6, 3)
    ends = dict(head="declared", tail="declared")
    first = alternating.build("a", dc.ChainSpec(2, **ends))
    second = uniform.build("u", dc.ChainSpec(4, **ends))
    assert heavy_formula(first) == heavy_formula(second) == {"C": 20, "H": 18}
    assert not first.molecule.is_isomorphic(second.molecule)


def ends_mapping():
    mapping = pp_mapping()
    methyl = {"atoms": [_atom("m", hydrogens=3)], "bonds": [],
              "ports": {"cap": {"atom": "m", "bond_order": 1}},
              "electronic_state": {"multiplicity": 1}}
    hydroxyl = {"atoms": [_atom("o", "O", hydrogens=1, lone_pairs=2)], "bonds": [],
                "ports": {"cap": {"atom": "o", "bond_order": 1}},
                "electronic_state": {"multiplicity": 1}}
    mapping["end_groups"] = {"head": methyl, "tail": hydroxyl}
    mapping["locality_policy"]["bounds"]["elements"] = ["C", "H", "O"]
    mapping["locality_policy"]["bounds"]["lone_pairs"] = [0, 2]
    return mapping


def test_heavy_declared_ends_are_reconstructed():
    provider = dc.DeclaredChainProvider(spec_of(ends_mapping()))
    open_ = sup.identity_witness(provider.build("o", dc.ChainSpec(2)))
    built, record = b.reconstruct_end(open_, (0, "tail"), "declared")
    built, record2 = b.reconstruct_end(built, (0, "head"), "declared")
    both = dc.ChainSpec(2, head="declared", tail="declared")
    direct = sup.identity_witness(provider.build("d", both))
    assert sup.witnesses_match(direct, built)
    assert [a[1][0] for a in record.fresh_atoms] == ["O", "H"]
    assert [a[1][0] for a in record2.fresh_atoms] == ["C", "H", "H", "H"]
    one = b.reconstruct_end(open_, (0, "head"), "declared")[0]
    half = sup.identity_witness(provider.build("h", dc.ChainSpec(2, head="declared")))
    assert sup.witnesses_match(half, one)


def cyclopentane_mapping():
    """A cyclopentane-1,3-diyl backbone: head r1, tail r3; r2 is the short
    backbone side, r4-r5 the long one."""
    atoms = [_atom("r1", hydrogens=1), _atom("r2", hydrogens=2), _atom("r3", hydrogens=1),
             _atom("r4", hydrogens=2), _atom("r5", hydrogens=2)]
    bonds = [_bond("r1", "r2"), _bond("r2", "r3"), _bond("r3", "r4"), _bond("r4", "r5"),
             _bond("r5", "r1")]
    return _declaration(atoms, bonds, "r1", "r3")


def test_backbone_atoms_follow_deletions():
    provider = dc.DeclaredChainProvider(spec_of(cyclopentane_mapping()))
    assert provider.n_bb == 3
    pristine = dc.ChainSpec(2)
    assert provider.backbone_atoms(pristine) == 6
    opened = dc.ChainSpec(2, edits=(E("delete_atom", 1, "r2"), E("set_radical", 1, "r1", 1),
                                    E("set_radical", 1, "r3", 1)))
    provider.build("x", opened)
    # the built head -> tail path now runs r1-r5-r4-r3 in unit 1
    assert provider.backbone_atoms(opened) == 7


def _built_path(build):
    import networkx as nx
    r = rows_of(build)
    graph = nx.Graph()
    for a, row in r.items():
        if row[0] != "H":
            graph.add_node(a)
            for c, _o in row[5]:
                if r[c][0] != "H":
                    graph.add_edge(a, c)
    return nx.shortest_path(graph, build.metadata.port("head").atom_id,
                            build.metadata.port("tail").atom_id)


def test_backbone_index_follows_the_built_path_after_a_deletion():
    provider = dc.DeclaredChainProvider(spec_of(cyclopentane_mapping()))
    opened = provider.build("x", dc.ChainSpec(2, edits=(
        E("delete_atom", 1, "r2"), E("set_radical", 1, "r1", 1), E("set_radical", 1, "r3", 1))))
    coords = {a: (u, k) for a, u, k in opened.metadata.coordinates}
    path = _built_path(opened)
    indices = provider.backbone_indices(opened)
    walk = [indices.get(a) for a in path]
    assert walk == list(range(7)), walk
    assert [provider.backbone_index(*coords[a], build=opened) for a in path] == walk
    # off-path atoms of the opened unit carry no index
    p = provider.geometry.positions
    assert indices.get(at(opened, 0, "r4", provider)) is None
    # a pristine build agrees with the declaration-only index
    pristine = provider.build("p", dc.ChainSpec(2))
    pc = {a: (u, k) for a, u, k in pristine.metadata.coordinates}
    built = provider.backbone_indices(pristine)
    for a in _built_path(pristine):
        assert built[a] == provider.backbone_index(*pc[a])
    assert p["r2"] == 1


def test_built_backbone_index_starts_at_a_phase_cut(ps):
    build = ps.provider.build("s", dc.ChainSpec(3, head_phase=ps.phase))
    path = _built_path(build)
    indices = ps.provider.backbone_indices(build)
    assert [indices[a] for a in path] == list(range(1, 6))


def test_overlapping_legacy_flags_translate(ps):
    for spec in (b.PsChainSpec(3, head_end="methyl", tail_end="benzyl", radical_units=(1,),
                               pendant_loss_units=(1,)),
                 b.PsChainSpec(3, benzylic_end=True, radical_units=(2,)),
                 b.PsChainSpec(3, radical_units=(1, 1))):
        mismatches, refused = sup.ps_parity(ps, [spec])
        assert not mismatches and not refused, (spec, mismatches, refused)


def test_overlap_catalog_is_large():
    assert len(CATALOG["overlap"]) >= 250
