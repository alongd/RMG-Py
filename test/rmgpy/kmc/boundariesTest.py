"""Tests for rmgpy/kmc/boundaries.py: construction-carried boundary identity,
exact explicit-H extension maps, and per-port extension vectors.

Fast part only (a few small RMG generations; no kinetics).  The full PS
acceptance sweep is ``boundaries_acceptance.py`` (not collected).
"""

import hashlib
import os
import random
import subprocess
import sys
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest

from rmgpy.kmc import boundaries as b
from rmgpy.kmc import rmg_moves as rm
from rmgpy.molecule.molecule import Molecule

import rmg_moves_harness as harness

TEST_DIR = Path(__file__).resolve().parent
ROOT = TEST_DIR.parents[2]

pytestmark = pytest.mark.skipif(not harness.DATABASE.is_dir(),
                                reason=f"needs the RMG database at {harness.DATABASE}")


@pytest.fixture(scope="module")
def fam():
    return harness.families()


def kept_witnesses(fam, builds, family, direction="forward"):
    molecules, participants = b.witness_reactants(builds)
    forms = [rm.reactant_forms(m) for m in molecules]
    moves = rm.generate_moves(fam[family], molecules, direction)
    return [b.witness_from_move(m, molecules, participants, forms=forms)
            for m in moves if m.status == "kept"]


def alkane(n):
    return b.PROVIDERS["linear_alkane"].build(n)


def homolysis(fam, n, bond):
    """The C(n) homolysis witness of carbon bond ``bond`` (ids along the chain)."""
    found = [w for w in kept_witnesses(fam, [alkane(n)], "R_Recombination", "reverse")
             if set(b.edits(w)[1][0][0]) == set(bond) and len(b.edits(w)[1]) == 1]
    assert found, f"no homolysis of {bond}"
    return found[0]


def rows(side):
    return {a: row for a, row in side}


def side_without(side, atom_ids):
    atoms = {a: row[:5] for a, row in side if a not in atom_ids}
    bonds = {}
    for a, row in side:
        for c, o in row[5]:
            if a not in atom_ids and c not in atom_ids:
                bonds[(min(a, c), max(a, c))] = o
    return b._side_from(atoms, bonds)


def relabel_side(side, mapping):
    atoms = {mapping.get(a, a): row[:5] for a, row in side}
    bonds = {}
    for a, row in side:
        for c, o in row[5]:
            x, y = mapping.get(a, a), mapping.get(c, c)
            bonds[(min(x, y), max(x, y))] = o
    return b._side_from(atoms, bonds)


def hydrogens_of(side, atom_id):
    r = rows(side)
    return sorted(n for n, _o in r[atom_id][5] if r[n][0] == "H")


# --------------------------------------------------------------------------
# item 1: construction-carried identity
# --------------------------------------------------------------------------
EXPECTED_PORTS = {
    "pristine": {"head": "artificial", "tail": "artificial"},
    "interior_radical": {"head": "artificial", "tail": "artificial"},
    "doubly_featured": {"head": "artificial", "tail": "artificial"},
    "end_radical": {"head": "physical", "tail": "artificial"},
    "end_radical_short": {"head": "physical", "tail": "artificial"},
    "benzylic_end_radical": {"head": "artificial", "tail": "physical"},
    "benzylic_end_radical_short": {"head": "artificial", "tail": "physical"},
    "junction_radical": {"head": "artificial", "tail": "artificial"},
    "styrene": {},
}


@pytest.mark.parametrize("units", [3, 5])
def test_ps_builds_carry_ports_caps_and_repeat_coordinates(units):
    for name, expected in EXPECTED_PORTS.items():
        build = b.ps_build(name, units)
        meta = build.metadata
        r = dict(b.rows_of_molecule(build.molecule))
        assert meta.provider == "polystyrene"
        assert {p.port_id: p.physicality for p in meta.ports} == expected
        tags = {p.port_id: p.tag for p in meta.ports}
        assert tags == {k: f"backbone_{k}" for k in expected}
        for port in meta.ports:
            if port.artificial:
                (cap,) = port.cap_atom_ids
                assert r[cap][0] == "H" and r[cap][5] == ((port.atom_id, 1.0),)
                # the declared cap is the LAST H of the port added at construction
                port_h = [a for a in build.construction_order
                          if r[a][0] == "H" and r[a][5][0][0] == port.atom_id]
                assert cap == port_h[-1]
            else:
                assert port.cap_atom_ids == ()
        heavy = sorted(a for a, row in r.items() if row[0] != "H")
        assert sorted(a for a, _u, _k in meta.coordinates) == heavy
        coords = {a: (u, k) for a, u, k in meta.coordinates}
        # role sanity from the graph (test side only): position 1 bonds the
        # ring ipso (position 2); ring positions are in a cycle.
        by_role = {v: a for a, v in coords.items()}
        n_units = 1 if name == "styrene" else len(heavy) // 8
        for u in range(n_units):
            assert any(c == by_role[(u, 2)] for c, _o in r[by_role[(u, 1)]][5])
            assert build.molecule.is_atom_in_cycle(
                next(x for x in build.molecule.atoms if x.id == by_role[(u, 5)]))
        if "head" in expected:
            assert coords[meta.port("head").atom_id] == (0, 0)
            assert coords[meta.port("tail").atom_id] == (n_units - 1, 1)


def test_ps_builds_reproduce_the_existing_proxy_smiles():
    from rmgpy.kmc.compiler import benzylic_ps_end_smiles, linear_ps_smiles

    for units in (3, 5):
        legacy = {
            "pristine": linear_ps_smiles(units),
            "interior_radical": linear_ps_smiles(units, radical_unit=units // 2),
            "doubly_featured": linear_ps_smiles(units, radical_units=(0, 1)),
            "end_radical": linear_ps_smiles(units, radical_end=True),
            "end_radical_short": linear_ps_smiles(units - 1, radical_end=True),
            "benzylic_end_radical": benzylic_ps_end_smiles(units),
            "benzylic_end_radical_short": benzylic_ps_end_smiles(units - 1),
            "junction_radical": linear_ps_smiles(units).replace("c1ccccc1", "c1cc[c]cc1", 1),
            "styrene": "C=Cc1ccccc1",
        }
        specs = b.ps_structure_specs(units)
        assert set(specs) == set(legacy)
        for name, smiles in legacy.items():
            assert b.PolystyreneProvider.smiles(specs[name]) == smiles


def test_ps_proxy_set_content_is_unchanged_and_carries_builds():
    from rmgpy.kmc import compiler as cm

    for proxy in cm.ps_proxy_set(3):
        assert len(proxy.boundary_builds) == len(proxy.reactants)
        legacy_ports = []
        for index, (species, build) in enumerate(zip(proxy.reactants, proxy.boundary_builds)):
            molecule = species.molecule[0]
            assert build.molecule is molecule
            legacy = Molecule(smiles=b.PolystyreneProvider.smiles(
                b.ps_structure_specs(proxy.metadata["proxy_units"])[build.name]))
            assert legacy.is_isomorphic(molecule, strict=True)
            legacy_ports += list(cm._ps_artificial_ports(build.name, molecule, index))
        # the metadata-sourced ports equal the old inference on every PS proxy
        assert [(p.atom_id, p.orientation, p.participant_index, p.kind)
                for p in proxy.artificial_boundaries] == [
            (p.atom_id, p.orientation, p.participant_index, p.kind) for p in legacy_ports]
        assert all(p.source == "construction" and len(p.cap_atom_ids) == 1
                   for p in proxy.artificial_boundaries)
        assert all(p.source == "inferred" for p in legacy_ports)


def test_compiler_orientation_comes_from_metadata_not_phenyl_attachment():
    from rmgpy.kmc import compiler as cm

    build = alkane(6)
    from_metadata = cm._ps_artificial_ports("C6", build.molecule, 0, metadata=build.metadata)
    assert sorted((p.atom_id, p.orientation) for p in from_metadata) == [(0, "head"), (5, "tail")]
    assert {p.cap_atom_ids for p in from_metadata} == {
        port.cap_atom_ids for port in build.metadata.ports}
    # the legacy inference (no phenyl anywhere) calls both ends "head"
    inferred = cm._ps_artificial_ports("C6", build.molecule, 0)
    assert sorted(p.orientation for p in inferred) == ["head", "head"]


def test_alkane_build_is_explicit_unbranched_with_methyl_ends():
    build = alkane(8)
    r = dict(b.rows_of_molecule(build.molecule))
    carbons = [a for a in sorted(r) if r[a][0] == "C"]
    assert carbons == list(range(8))
    for c in carbons:
        heavy = [n for n, _o in r[c][5] if r[n][0] == "C"]
        assert sorted(heavy) == sorted(x for x in (c - 1, c + 1) if 0 <= x < 8)
    assert build.metadata.port("head") == b.Port("head", 0, (hydrogens_of(tuple(r.items()), 0)[-1],),
                                                  "backbone_head", "artificial")
    assert build.metadata.port("tail").atom_id == 7
    assert build.metadata.coordinates == tuple((i, i, 0) for i in range(8))
    assert build.molecule.to_smiles() == "CCCCCCCC"


@pytest.mark.parametrize("name,family,direction", [
    ("pristine", "R_Recombination", "reverse"),
    ("interior_radical", "intra_H_migration", "forward"),
    ("junction_radical", "intra_H_migration", "forward"),
    ("end_radical", "intra_H_migration", "forward"),
    ("end_radical", "R_Addition_MultipleBond", "reverse"),
    ("benzylic_end_radical", "intra_H_migration", "forward"),
])
def test_metadata_survives_rmg_generation(fam, name, family, direction):
    build = b.ps_build(name, 3)
    witnesses = kept_witnesses(fam, [build], family, direction)
    assert witnesses
    for w in witnesses:
        meta = w.participants[0].metadata
        r, p = w.reactant_rows(), w.product_rows()
        for port in meta.ports:
            assert r[port.atom_id][0] == p[port.atom_id][0] == "C"
            for cap in port.cap_atom_ids:
                assert r[cap][0] == p[cap][0] == "H"
        assert {a for a, _u, _k in meta.coordinates} == {a for a, row in r.items() if row[0] != "H"}


def test_metadata_survives_bimolecular_generation_with_disjoint_ids(fam):
    builds = [b.ps_build("end_radical", 3), b.ps_build("end_radical", 3)]
    witnesses = kept_witnesses(fam, builds, "R_Recombination")
    assert witnesses
    w = witnesses[0]
    assert [p.index for p in w.participants] == [0, 1]
    ids0, ids1 = (set(p.atom_ids) for p in w.participants)
    assert not ids0 & ids1
    assert w.port_keys() == ((0, "head"), (0, "tail"), (1, "head"), (1, "tail"))
    assert b.vector_ports(w) == ((0, "tail"), (1, "tail"))


def test_lost_or_changed_metadata_atom_is_reported_with_the_family(fam):
    build = b.ps_build("interior_radical", 3)
    molecules, participants = b.witness_reactants([build])
    moves = [m for m in rm.generate_moves(fam["intra_H_migration"], molecules) if m.status == "kept"]
    move = moves[0]
    cap = participants[0].metadata.port("tail").cap_atom_ids[0]
    port = participants[0].metadata.port("tail").atom_id
    products = tuple((a, row) for a, row in move.products if a != cap)
    products = tuple((a, row[:5] + (tuple(e for e in row[5] if e[0] != cap),))
                     for a, row in products)
    with pytest.raises(b.MetadataLostError, match="intra_H_migration"):
        b.witness_from_move(replace(move, products=products), molecules, participants)
    relabelled = tuple((a + 10_000 if a == port else a,
                        row[:5] + (tuple((n + 10_000 if n == port else n, o) for n, o in row[5]),))
                       for a, row in move.products)
    with pytest.raises(b.MetadataLostError, match="intra_H_migration"):
        b.witness_from_move(replace(move, products=relabelled), molecules, participants)


# --------------------------------------------------------------------------
# item 3: exact maps
# --------------------------------------------------------------------------
def test_single_extension_map_is_exact_and_deletes_the_declared_cap(fam):
    w = homolysis(fam, 8, (3, 4))
    for key in b.vector_ports(w):
        port = w.port(key)
        others = [h for h in hydrogens_of(w.reactants, port.atom_id) if h not in port.cap_atom_ids]
        # the declared cap is NOT the first neighbouring H
        assert min(hydrogens_of(w.reactants, port.atom_id)) not in port.cap_atom_ids
        ext, rec = b.extend(w, key)
        assert rec.deleted == port.cap_atom_ids
        assert [s[0] for s in rec.surviving] == sorted(set(w.reactant_rows()) - set(port.cap_atom_ids))
        for side in (ext.reactants, ext.products):
            ids = set(rows(side))
            assert not ids & set(port.cap_atom_ids)
            assert set(others) <= ids
            assert rows(side)[port.atom_id][5] == tuple(sorted(
                [e for e in rows(w.reactants)[port.atom_id][5] if e[0] not in port.cap_atom_ids]
                + [(rec.steps[0].fresh_ids[0], 1.0)]))
        new = ext.port(key)
        assert new.atom_id in rec.fresh_ids and set(new.cap_atom_ids) <= set(rec.fresh_ids)
        assert new.tag == port.tag and new.artificial
        assert ext.participants[0].metadata.coordinate(new.atom_id) == (
            (-1, 0) if key[1] == "head" else (8, 0))
        b.verify_extension(w, rec, ext)


def test_ps_extension_at_head_and_tail(fam):
    w = kept_witnesses(fam, [b.ps_build("pristine", 3)], "R_Recombination", "reverse")[0]
    for key in b.vector_ports(w):
        if w.port(key).atom_id in b.edit_support(w):
            continue
        ext, rec = b.extend(w, key)
        heavy = [a for a, attrs in rec.fresh_atoms if attrs[0] != "H"]
        assert len(heavy) == 8 and len(rec.fresh_atoms) == 17
        unit = ext.participants[0].metadata.coordinate(ext.port(key).atom_id)
        assert unit == ((-1, 0) if key[1] == "head" else (3, 1))
        ext2, rec2 = b.extend(ext, key)
        b.verify_extension(ext, rec2, ext2)
        b.verify_extension(w, b.compose(rec, rec2), ext2)


def _corrupt_cases(w, key):
    ext, rec = b.extend(w, key)
    port = w.port(key)
    cap = port.cap_atom_ids[0]
    other_h = [h for h in hydrogens_of(w.reactants, port.atom_id) if h != cap][0]
    survivor_h = hydrogens_of(w.reactants, 5)[0]
    far_h = hydrogens_of(w.reactants, 6)[0]
    fresh0 = rec.steps[0].fresh_ids[0]
    neighbour = next(n for n, _o in rows(w.reactants)[port.atom_id][5]
                     if rows(w.reactants)[n][0] == "C")
    cases = {}
    # a dropped H (both sides, so the sides still agree)
    cases["dropped_h"] = (rec, replace(ext, reactants=side_without(ext.reactants, {survivor_h}),
                                       products=side_without(ext.products, {survivor_h})))
    # a duplicated id: a fresh atom reuses a surviving id
    dup_atoms = tuple(((survivor_h if a == rec.fresh_atoms[-1][0] else a), attrs)
                      for a, attrs in rec.fresh_atoms)
    cases["duplicated_id"] = (replace(rec, fresh_atoms=dup_atoms), ext)
    cases["duplicated_row"] = (rec, replace(ext, reactants=ext.reactants + (ext.reactants[0],)))
    # the cap deleted on one side only
    one_side = b._apply_side(w.products, set(), rec.fresh_atoms, rec.fresh_bonds, "corrupt")
    cases["cap_one_side"] = (rec, replace(ext, products=one_side))
    # the wrong H deleted on BOTH sides (record and graph agree with each other)
    wrong_step = replace(rec.steps[0], deleted=(other_h,))
    wrong = replace(rec, steps=(wrong_step,), deleted=(other_h,),
                    surviving=tuple(s for s in b.identity_map(w).surviving if s[0] != other_h))
    cases["wrong_h_both_sides"] = (wrong, b.apply_map(w, wrong))
    cases["wrong_h_graph_only"] = (rec, b.apply_map(w, wrong))
    # two surviving H ids exchanged (both sides, then one side)
    swap = {survivor_h: far_h, far_h: survivor_h}
    cases["h_swapped_both"] = (rec, replace(ext, reactants=relabel_side(ext.reactants, swap),
                                            products=relabel_side(ext.products, swap)))
    cases["h_swapped_one_side"] = (rec, replace(ext, products=relabel_side(ext.products, swap)))
    # a fresh atom bonded to the wrong atom (record and graph agree)
    bad_bonds = tuple(sorted((min(neighbour, fresh0), max(neighbour, fresh0), o)
                             if {x, y} == {port.atom_id, fresh0} else (x, y, o)
                             for x, y, o in rec.fresh_bonds))
    wrong_bond = replace(rec, fresh_bonds=bad_bonds)
    cases["fresh_wrong_atom"] = (wrong_bond, b.apply_map(w, replace(
        wrong_bond, deleted=wrong_bond.deleted)))
    cases["fresh_wrong_atom_graph_only"] = (rec, b.apply_map(w, wrong_bond))
    return ext, rec, cases


CORRUPTIONS = ["dropped_h", "duplicated_id", "duplicated_row", "cap_one_side",
               "wrong_h_both_sides", "wrong_h_graph_only", "h_swapped_both",
               "h_swapped_one_side", "fresh_wrong_atom", "fresh_wrong_atom_graph_only"]


@pytest.mark.parametrize("corruption", CORRUPTIONS)
def test_corrupted_extension_is_refused(fam, corruption):
    w = homolysis(fam, 8, (3, 4))
    ext, rec, cases = _corrupt_cases(w, (0, "head"))
    b.verify_extension(w, rec, ext)  # the honest one passes
    bad_rec, bad_ext = cases[corruption]
    with pytest.raises(b.MapCheckError):
        b.verify_extension(w, bad_rec, bad_ext)


def test_extension_deletes_the_cap_on_both_sides(fam):
    w = homolysis(fam, 8, (2, 3))
    ext, rec = b.extend(w, (0, "tail"))
    cap = w.port((0, "tail")).cap_atom_ids[0]
    assert cap not in rows(ext.reactants) and cap not in rows(ext.products)
    assert set(rows(ext.reactants)) == set(rows(ext.products))


def test_edit_support_atoms_cannot_be_deleted(fam):
    # an H abstraction / migration moving the declared cap of a port
    build = b.ps_build("end_radical", 3)
    witnesses = kept_witnesses(fam, [build], "intra_H_migration")
    tail = (0, "tail")
    moved_cap = [w for w in witnesses if w.port(tail).cap_atom_ids[0] in b.edit_support(w)]
    assert moved_cap, "no migration moves the tail cap"
    w = moved_cap[0]
    with pytest.raises(b.ExtensionRefused, match="reacting atom"):
        b.extend(w, tail)
    # forcing the record past the port check still cannot delete a support atom
    clean = next(x for x in witnesses if x.port(tail).atom_id not in b.edit_support(x))
    _ext, rec = b.extend(clean, tail)
    with pytest.raises(b.MapCheckError, match="edit-support"):
        b.apply_map(w, rec)
    with pytest.raises(b.MapCheckError, match="edit-support"):
        b.verify_extension(w, rec, _ext)
    with pytest.raises(b.MapCheckError, match="edit-support"):
        b.plan_extension(w, tail, b.PROVIDERS["polystyrene"], support=frozenset(
            w.port(tail).cap_atom_ids))


def test_reacting_and_physical_ports_are_refused(fam):
    w = homolysis(fam, 8, (0, 1))
    with pytest.raises(b.ExtensionRefused, match="reacting"):
        b.extend(w, (0, "head"))
    b.extend(w, (0, "tail"))
    end = kept_witnesses(fam, [b.ps_build("end_radical", 3)], "intra_H_migration")[0]
    with pytest.raises(b.ExtensionRefused, match="physical"):
        b.extend(end, (0, "head"))


def test_composition_equals_sequential_application(fam):
    w = homolysis(fam, 16, (7, 8))
    e1, r1 = b.extend(w, (0, "head"))
    e2, r2 = b.extend(e1, (0, "tail"))
    e3, r3 = b.extend(e2, (0, "head"))
    composite = b.compose(b.compose(r1, r2), r3)
    assert b.apply_map(w, composite) == e3
    assert b.compose(r1, b.compose(r2, r3)) == composite
    b.verify_extension(w, composite, e3)
    # the second head step deleted the first head step's fresh cap
    assert r3.deleted == r1.steps[0].new_port.cap_atom_ids
    assert set(composite.deleted) == {w.port(k).cap_atom_ids[0] for k in b.vector_ports(w)}
    assert len(composite.fresh_atoms) == 3 * 4 - 1


class StarEtherProvider(b.BoundaryProvider):
    """Toy chemistry: three ports, an O-CH2 repeat (heteroatom continuation)
    and a two-atom OH cap group."""

    name = "toy_star_ether"

    def build(self):
        from rmgpy.kmc.rmg_moves import molecule_from_rows

        # central CH(-CH2-OH)3 style: C0 with three arms C1, C2, C3 each with an
        # O-H cap group (O at 4..6, H at 7..9); ring of H on carbons
        rows_ = {}
        def add(a, el, edges):
            rows_[a] = (el, None, 0, 0, 2 if el == "O" else 0, tuple(sorted(edges)))
        add(0, "C", [(1, 1.0), (2, 1.0), (3, 1.0), (10, 1.0)])
        hid = 11
        for arm in (1, 2, 3):
            o = arm + 3
            add(arm, "C", [(0, 1.0), (o, 1.0), (hid, 1.0), (hid + 1, 1.0)])
            add(hid, "H", [(arm, 1.0)])
            add(hid + 1, "H", [(arm, 1.0)])
            hid += 2
            add(o, "O", [(arm, 1.0), (o + 3, 1.0)])
            add(o + 3, "H", [(o, 1.0)])
        add(10, "H", [(0, 1.0)])
        molecule = molecule_from_rows(rows_)
        ports = tuple(b.Port(f"arm{arm}", arm, (arm + 3, arm + 6), "other", "artificial")
                      for arm in (1, 2, 3))
        meta = b.BoundaryMetadata(self.name, ports, ((0, 0, 0), (1, 1, 0), (2, 2, 0), (3, 3, 0),
                                                     (4, 1, 1), (5, 2, 1), (6, 3, 1)))
        return b.ProxyBuild("star", molecule, meta, tuple(sorted(rows_)))

    def fragment(self, port, coordinate):
        unit, _ = coordinate
        # O (0) - CH2 (1) + new cap O (2) - H (3); H on CH2: 4, 5
        atoms = (("O", None, 0, 0, 2), ("C", None, 0, 0, 0), ("O", None, 0, 0, 2),
                 ("H", None, 0, 0, 0), ("H", None, 0, 0, 0), ("H", None, 0, 0, 0))
        bonds = ((0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (1, 4, 1.0), (1, 5, 1.0))
        return b.Fragment(atoms, bonds, ((0, 1.0),), 1, (2, 3),
                          ((0, unit + 100, 0), (1, unit + 100, 1), (2, unit + 100, 2)))


def test_core_assumes_no_two_ports_no_carbon_continuation_no_h_cap():
    provider = StarEtherProvider()
    build = provider.build()
    molecules, participants = b.witness_reactants([build])
    side = b.rows_of_molecule(molecules[0])
    w = b.make_witness(side, side, participants)
    providers = {provider.name: provider}
    assert len(b.vector_ports(w)) == 3
    for result in b.enumerate_vectors(w, 1, providers, ceiling=50):
        assert result.status == "ok"
        b.verify_extension(w, result.record, result.witness, providers)
    ext, rec = b.extend_vector(w, (2, 0, 1), providers, ceiling=50)
    assert len(rec.deleted) == 4  # two two-atom cap groups
    assert sum(1 for _a, attrs in rec.fresh_atoms if attrs[0] == "O") == 3 + 2


# --------------------------------------------------------------------------
# item 4: extension vectors
# --------------------------------------------------------------------------
def test_vector_ports_are_ordered_by_participant_then_port_id(fam):
    builds = [b.ps_build("pristine", 3), b.ps_build("end_radical", 3)]
    w = kept_witnesses(fam, builds, "H_Abstraction")[0]
    assert b.vector_ports(w) == ((0, "head"), (0, "tail"), (1, "tail"))


def test_enumeration_box_statuses_and_ceiling(fam):
    w = homolysis(fam, 8, (3, 4))
    results = b.enumerate_vectors(w, 2)
    assert [r.vector for r in results] == [(i, j) for i in range(3) for j in range(3)]
    assert all(r.status == "ok" for r in results)
    capped = b.enumerate_vectors(w, (2, 2), ceiling=10)
    status = {r.vector: r.status for r in capped}
    assert status[(1, 1)] == "ok" and status[(2, 0)] == "ok"
    assert status[(2, 1)] == "ceiling_exceeded" and status[(2, 2)] == "ceiling_exceeded"
    with pytest.raises(b.CeilingExceeded):
        b.extend_vector(w, (2, 1), ceiling=10)
    terminal = homolysis(fam, 8, (0, 1))
    res = {r.vector: r for r in b.enumerate_vectors(terminal, 1)}
    assert res[(1, 0)].status == "refused" and "reacting" in res[(1, 0)].reason
    assert res[(0, 1)].status == "ok"


def test_default_ceiling_is_the_compiler_constant(fam):
    from rmgpy.kmc.compiler import PROXY_PADDING_MAX_HEAVY_ATOMS

    w = homolysis(fam, 8, (3, 4))
    assert b._ceiling(None) == PROXY_PADDING_MAX_HEAVY_ATOMS
    big = 1 + (PROXY_PADDING_MAX_HEAVY_ATOMS - 8) // 2
    assert b.evaluate_vector(w, (big, big)).status == "ceiling_exceeded"
    assert b.evaluate_vector(w, (big - 1, big - 1), verify=False,
                             check_order=False).status == "ok"
    with pytest.raises(b.CeilingExceeded):
        b.extend_vector(w, (big, big), verify=False)


def test_extension_order_does_not_matter(fam):
    w = homolysis(fam, 8, (3, 4))
    a, ra = b.extend_vector(w, (2, 1), order=(0, 1))
    c, rc = b.extend_vector(w, (2, 1), order=(1, 0))
    assert a == c and ra == rc
    builds = [b.ps_build("pristine", 3), b.ps_build("end_radical", 3)]
    pair = next(x for x in kept_witnesses(fam, builds, "H_Abstraction")
                if not {x.port(k).atom_id for k in b.vector_ports(x)} & b.edit_support(x))
    orders = [(0, 1, 2), (2, 1, 0), (1, 2, 0)]
    out = [b.extend_vector(pair, (1, 2, 1), order=o) for o in orders]
    assert all(x == out[0] for x in out)


def _unmarked(witness):
    g = b.marked_graph(witness)
    for _n, data in g.nodes(data=True):
        data["label"] = data["label"].rsplit(",", 1)[0]
    return g


def test_marked_distinctness(fam):
    from networkx.algorithms import isomorphism

    w = homolysis(fam, 8, (3, 4))  # symmetric: C0..C3 | C4..C7
    results = {r.vector: r for r in b.enumerate_vectors(w, 1)}
    left, right = results[(1, 0)].witness, results[(0, 1)].witness
    # isomorphic only by moving marked atoms (the two reacting carbons, the ports)
    gl, gr = _unmarked(left), _unmarked(right)
    assert isomorphism.GraphMatcher(
        gl, gr, node_match=lambda x, y: x["label"] == y["label"],
        edge_match=lambda x, y: x["label"] == y["label"]).is_isomorphic()
    assert not b.marked_isomorphic(left, right)
    assert b.vector_duplicates(tuple(results.values())) == ()
    # isomorphic preserving every mark: unmarked H ids permuted -> duplicate
    h5 = hydrogens_of(left.reactants, 5)
    h6 = hydrogens_of(left.reactants, 6)
    swap = {h5[0]: h6[0], h6[0]: h5[0]}
    twin = replace(left, reactants=relabel_side(left.reactants, swap),
                   products=relabel_side(left.products, swap))
    assert twin != left
    assert b.marked_isomorphic(left, twin)
    assert b.marked_duplicates([("left", left), ("right", right), ("twin", twin)]) == (
        ("left", "twin"),)


# --------------------------------------------------------------------------
# item 5: compiler wiring (padding path reads the declared cap)
# --------------------------------------------------------------------------
def test_padding_path_removes_the_declared_cap_on_both_sides():
    from rmgpy.kmc.proxy_padding import BoundaryPort, RepeatUnitGraph, pad_reaction_witness

    molecule = Molecule(smiles="CCC")
    molecule.assign_atom_ids()
    heavy = [a for a in molecule.atoms if not a.is_hydrogen()]
    head, tail = heavy[0], heavy[-1]
    tail_h = sorted(n.id for n in tail.edges if n.is_hydrogen())
    first_h = next(n.id for n in tail.edges if n.is_hydrogen())
    declared = max(tail_h) if first_h != max(tail_h) else min(tail_h)
    reaction = SimpleNamespace(reactants=[molecule], products=[molecule.copy(deep=True)],
                               degeneracy=1.0)
    repeat = RepeatUnitGraph.from_smiles("PE", "CC", head_atom_index=0, tail_atom_index=1)
    padded = pad_reaction_witness(
        reaction, [BoundaryPort(tail.id, "tail", cap_atom_ids=(declared,), source="construction")],
        repeat, min_distance=4, reacting_atom_ids=[head.id])
    for side in (padded.reaction.reactants, padded.reaction.products):
        ids = {a.id for a in side[0].atoms}
        assert declared not in ids
        assert set(tail_h) - {declared} <= ids
    port = padded.boundaries[0]
    assert port.source == "construction" and len(port.cap_atom_ids) == 1
    legacy = pad_reaction_witness(reaction, [BoundaryPort(tail.id, "tail")], repeat,
                                  min_distance=4, reacting_atom_ids=[head.id])
    assert first_h not in {a.id for a in legacy.reaction.reactants[0].atoms}
    assert legacy.boundaries[0].source == "inferred"


# --------------------------------------------------------------------------
# item 6: determinism
# --------------------------------------------------------------------------
def digest(fam):
    parts = []
    for build in (alkane(8), b.ps_build("pristine", 3), b.ps_build("end_radical", 3)):
        parts.append(repr(build.canonical_metadata()))
        witnesses = kept_witnesses(fam, [build], "R_Recombination", "reverse")[:6]
        for w in witnesses:
            for r in b.enumerate_vectors(w, 1):
                parts.append(repr((r.vector, r.status, r.witness, r.record)))
    return hashlib.sha256("\n".join(parts).encode()).hexdigest()


def test_output_is_identical_across_processes_and_hash_seeds(fam):
    code = ("import sys; sys.path.insert(0, %r); import boundariesTest as t, "
            "rmg_moves_harness as h; print(t.digest(h.families()))" % str(TEST_DIR))
    out = []
    for seed in ("0", "12345"):
        env = dict(os.environ, PYTHONHASHSEED=seed, MPLCONFIGDIR=os.environ.get(
            "TMPDIR", "/tmp"))
        proc = subprocess.run([sys.executable, "-c", code], cwd=str(ROOT), env=env,
                              capture_output=True, text=True, timeout=600)
        assert proc.returncode == 0, proc.stderr[-2000:]
        out.append(proc.stdout.strip().splitlines()[-1])
    assert out[0] == out[1] == digest(fam)


def test_permuted_input_atom_orders_give_the_same_witnesses(fam):
    build = b.ps_build("end_radical", 3)
    reference = kept_witnesses(fam, [build], "intra_H_migration")
    rng = random.Random(7)
    for _ in range(2):
        molecule = build.molecule.copy(deep=True)
        for atom, original in zip(molecule.atoms, build.molecule.atoms):
            atom.id = original.id
        atoms = list(molecule.atoms)
        rng.shuffle(atoms)
        molecule.atoms = atoms
        permuted = replace(build, molecule=molecule)
        again = kept_witnesses(fam, [permuted], "intra_H_migration")
        assert sorted(map(repr, again)) == sorted(map(repr, reference))
    w = reference[0]
    keys = [k for k in b.vector_ports(w) if w.port(k).atom_id not in b.edit_support(w)]
    if keys:
        assert b.extend(w, keys[0]) == b.extend(w, keys[0])
