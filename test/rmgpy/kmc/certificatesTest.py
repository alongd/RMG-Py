"""Tests for rmgpy/kmc/certificates.py (complete rate-selection certificates,
record binding, recertification) and rmgpy/kmc/domain.py (the declared v1 PS
domain).

Fast part: root checks, predicate semantics and boundary signatures, the
unanchored (list-label) predicates of R_Recombination, depository search and
whole-species policy, every member and mapping certified, the EventRecord
reference, attachment, recertification through a renumbered extension (also
for a reversed record), the domain checker.
Slow part (``RMG_KMC_SLOW=1``): the acceptance tables on the PS references
(DP3, DP8 head radical, DP4 + 1-phenylethyl) and on the closed PS DP8
reference compile (certificates attached at compile time).
"""

import copy
import json
import os
import time
from collections import Counter, defaultdict
from dataclasses import replace
from pathlib import Path

import pytest

from rmgpy.kmc import certificates as cm
from rmgpy.kmc import channels as chm
from rmgpy.kmc import domain as dom
from rmgpy.kmc import rmg_moves as rm
from rmgpy.kmc.event_record import EventRecord

import channelsTest as ct
import rmg_moves_harness as harness

SLOW = os.environ.get("RMG_KMC_SLOW") == "1"
TEST_DIR = Path(__file__).resolve().parent

pytestmark = pytest.mark.skipif(not harness.DATABASE.is_dir(),
                                reason=f"needs the RMG database at {harness.DATABASE}")

mol = ct.mol
only = ct.only


@pytest.fixture(scope="module")
def db():
    return ct.database()


@pytest.fixture(scope="module")
def fam(db):
    return db.kinetics.families


@pytest.fixture(scope="module")
def certifier(fam):
    return cm.Certifier(fam)


def witness_atoms(channel_set):
    return {a.id: a for m in channel_set.reactants for a in m.atoms}


# --------------------------------------------------------------------------
# item 1: complete certificate
# --------------------------------------------------------------------------
def test_root_checks_semantics_and_boundary_signatures(fam, certifier):
    """CH3* + propane, methyl H (leaf C/H3/Cs\\H2\\Cs): every descent starts
    with its ROOT check; every group predicate is anchored, has a radius and a
    ball, and every ball atom has a row with its atom type and its FULL bond
    list (bonds leaving the ball included).  Node semantics carry atom types
    and props; the radius unit, H handling and global guards are stated."""
    family = fam["H_Abstraction"]
    cs = chm.build_channels(family, [mol("[CH3]", 0), mol("CCC", 10)])
    ch = only(cs.find(pathway=("C/H3/Cs\\H2\\Cs", "C_methyl")))
    cert = certifier.certify_channel(family, cs, ch)
    assert cert["complete"], cert["gaps"]
    assert "ordinary bonds" in cert["radius_unit"] and "explicit" in cert["h_handling"]
    assert {g["guard"] for g in cert["global_guards"]} >= {
        "multiplicity_convention", "depository_whole_species", "inRing",
        "group_multiplicity", "group_radical_count", "group_element_count"}
    witness = witness_atoms(cs)
    (c,) = cert["channels"]
    leaving = 0
    tops = set()
    for mem in c["members"]:
        for mp in mem["mappings"]:
            descents = cert["descents"][mp["descent"]]
            for d in descents:
                if d.get("skipped"):
                    continue
                assert d["root"]["node"] == d["top"]
                if d["leaf"] is None:   # RMG tried this top on the other reactant
                    assert d["root"]["result"] == "failed" and d["steps"] == []
                    continue
                tops.add(d["top"])
                assert d["root"]["result"] == "matched"
                table = {row[0]: row for row in d["atoms"]}
                for p in cm.iter_descent_predicates(d):
                    if p.get("logic"):
                        continue
                    if p["anchoring"] == "labels_missing":
                        assert p["result"] == "failed"
                        continue
                    assert p["anchoring"] == "anchored" and p["radius"] is not None
                    ball = {i for i, _d in p["ball"]}
                    assert set(p["anchors"].values()) <= ball
                    for i in ball:
                        row = table[i]
                        atom = witness[i]
                        assert row[2] == atom.atomtype.label
                        assert sorted(n for n, _o in row[8]) == sorted(b.id for b in atom.bonds)
                        leaving += sum(1 for n, _o in row[8] if n not in ball)
    assert tops == {"X_H", "Y_rad_birad_trirad_quadrad"} or len(tops) == 2
    assert leaving > 0
    semantics = cert["nodes"]["C/H3/Cs\\H2\\Cs"]
    assert semantics["kind"] == "group"
    assert all(a["atomtypes"] for a in semantics["atoms"])
    # the failed child reads the far CH3 hydrogens (outside the leaf's own atoms)
    failed = [p for d in cert["descents"].values() for x in d
              for p in cm.iter_descent_predicates(x) if p["node"] == "C/H3/Cs\\H2\\Cs|O"]
    assert failed and all(p["result"] == "failed" for p in failed)
    assert all({18, 19, 20} <= {i for i, _d in p["ball"]} for p in failed)
    print("CERT CH3*+propane", len(cm.canonical_json(cert)), "bytes")


def test_unanchored_list_label_predicates_search_the_whole_structure(fam, certifier):
    """R_Recombination labels both reacting atoms '*': RMG's
    match_node_to_structure builds no initial map for a list label, so every
    predicate searches the whole merged structure.  The certificate says so
    (unanchored), lists every candidate atom (each radical centre of the
    form, the remote ones included) and its ball covers them all."""
    family = fam["R_Recombination"]
    radical, _closed, h_atom = ct.dp3_molecules()
    seen_remote = False
    for cs in (chm.build_channels(family, [radical, h_atom]),
               chm.build_channels(family, ct.remote_benzylic_pair())):
        for ch in cs.channels:
            cert = certifier.certify_channel(family, cs, ch)
            assert cert["complete"], cert["gaps"]
            for d in (x for v in cert["descents"].values() for x in v):
                assert d["structure"] == "merged reactants"
                labelled = set(d["labelled"]["*"])
                for p in cm.iter_descent_predicates(d):
                    if p.get("logic"):
                        continue
                    assert p["anchoring"] == "unanchored", p["node"]
                    assert p["scan"].startswith("whole structure")
                    cands = set(p["candidates"]["*"])
                    assert cands <= {i for i, _d in p["ball"]}
                    if p["result"] == "matched":
                        assert p["embeddings"]["count"] >= 1
                        assert set(p["embeddings"]["images"]) <= {i for i, _d in p["ball"]}
                    seen_remote = seen_remote or bool(cands - labelled)
    assert seen_remote


def test_logic_nodes_record_rmg_evaluation_order(fam, certifier):
    """A LogicOr top (H_Abstraction Y_rad...): components in RMG's order,
    stopping at the first match, and the recorded verdict is RMG's."""
    family = fam["H_Abstraction"]
    cs = chm.build_channels(family, [mol("[CH3]", 0), mol("CC", 10)])
    cert = certifier.certify_channel(family, cs, only(cs.channels))
    logic = [p for v in cert["descents"].values() for d in v
             for p in cm.iter_descent_predicates(d) if p.get("logic")]
    assert logic
    for p in logic:
        if p["logic"] == "LogicOr" and p["result"] == "matched" and not p["invert"]:
            assert p["components"][-1]["result"] == "matched"
            assert all(c["result"] == "failed" for c in p["components"][:-1])


# --------------------------------------------------------------------------
# item 2: depository / training policy
# --------------------------------------------------------------------------
def test_depository_search_recorded_whole_species(fam, certifier):
    """CH3* + C2H6: the training depository is searched (every entry, by whole
    reaction isomorphism), entry 215 matches and is selected; the certificate
    is whole_species and RMG's own get_kinetics agrees."""
    family = fam["H_Abstraction"]
    cs = chm.build_channels(family, [mol("[CH3]", 0), mol("CC", 10)])
    cert = certifier.certify_channel(family, cs, only(cs.channels))
    (c,) = cert["channels"]
    dep = c["depository"]
    (search,) = dep["searches"]
    assert search["depository"].endswith("training") and search["entries_searched"] > 1000
    assert 215 in [m["index"] for m in search["matches"]]
    assert dep["selected"]["index"] == 215 and "Reaction.is_isomorphic" in dep["identified_by"]
    assert c["source_kind"] == "training_exact_match" and c["whole_species"]
    assert c["support_scope"] == cm.WHOLE_SPECIES and "depository_exclusion" not in c
    assert c["rmg_check"]["agrees"] and cert["complete"]


def test_rule_estimated_channel_has_a_sound_depository_exclusion(fam, certifier):
    """DP3 benzylic radical + closed chain: the search is recorded with no
    match, and every long-chain analogue is excluded by size (a 24-heavy-atom
    chain against a largest training species of 22)."""
    family = fam["H_Abstraction"]
    radical, closed, _h = ct.dp3_molecules()
    cs = chm.build_channels(family, [radical, closed])
    ch = cs.channels[0]
    cert = certifier.certify_channel(family, cs, ch)
    (c,) = cert["channels"]
    assert c["depository"]["searches"][0]["matches"] == []
    assert c["depository"]["selected"] is None and not c["whole_species"]
    ex = c["depository_exclusion"]
    assert ex["verdict"] == "excluded:size", ex
    assert ex["depository_max_heavy_atoms"] < ex["witness_max_heavy_atoms"]
    assert c["rmg_check"]["agrees"]


def test_depository_exclusion_by_formula_and_refusal(fam):
    """Small witnesses: the formula argument excludes CH3* + propane-like
    analogues only when no entry has a grown formula; otherwise refused with
    the candidates (never assumed)."""
    family = fam["H_Abstraction"]
    ex = cm.depository_exclusion(family, [mol("[CH3]", 0), mol("CC", 10)],
                                 [mol("C", 20), mol("C[CH2]", 30)])
    assert ex["verdict"] in ("excluded:formula", "refused:unknown")
    if ex["verdict"] == "refused:unknown":
        assert ex["candidates"]
    print("EXCLUSION CH3*+C2H6", json.dumps(ex)[:400])


# --------------------------------------------------------------------------
# item 1: every member, every mapping
# --------------------------------------------------------------------------
def test_every_member_and_mapping_is_certified(fam, certifier):
    """A channel with several members (CH3* + propane, six methyl H) and a
    member with several mappings (the benzylic recombination found from three
    resonance forms): each member and each mapping has its own descent."""
    family = fam["H_Abstraction"]
    cs = chm.build_channels(family, [mol("[CH3]", 0), mol("CCC", 10)])
    ch = only(cs.find(pathway=("C/H3/Cs\\H2\\Cs", "C_methyl")))
    assert ch.size == 6
    cert = certifier.certify_channel(family, cs, ch)
    (c,) = cert["channels"]
    assert len(c["members"]) == ch.size
    assert [len(m["mappings"]) for m in c["members"]] == list(ch.mappings)
    donors = {next(a for _s, lab, a in map(tuple, mp["labels"]) if lab == "*2")
              for m in c["members"] for mp in m["mappings"]}
    assert len(donors) == 6
    data = json.loads((TEST_DIR / "fixtures" /
                       "i134_benzylic_recombination_witness.json").read_text())
    reactants = ct._witness(data["reactant_graphs"])
    rfam = fam["R_Recombination"]
    rcs = chm.build_channels(rfam, reactants)
    (rch,) = rcs.find(pathway=(data["template"],))
    rcert = certifier.certify_channel(rfam, rcs, rch)
    assert rcert["complete"], rcert["gaps"]
    (member,) = rcert["channels"][0]["members"]
    assert len(member["mappings"]) == rch.mappings[0] > 1
    assert len({mp["descent"] for mp in member["mappings"]}) > 1


# --------------------------------------------------------------------------
# item 3: record binding
# --------------------------------------------------------------------------
def _heavy_positions(molecules):
    """Original atom id -> position among the heavy atoms of the concatenated
    participants (the EventRecord ``atom_map`` index)."""
    heavy = [a.id for m in molecules for a in m.atoms if not a.is_hydrogen()]
    return {a: k for k, a in enumerate(heavy)}


def production_rate_source(family, cs, move):
    """The rate source the compiler serializes for an RMG family estimate
    (compiler._rate_table + _build_linked_family_pair): RMG's get_kinetics
    on the move's template reaction, its source/entry/rank, and for a
    hand-built tree compiler._rate_rule_source (rules, training, weights)."""
    from rmgpy.kmc import compiler as cmp
    rxn, _r, _p = chm.template_reaction(family, cs, move)
    kinetics, source, entry, _fwd = family.get_kinetics(
        rxn, template_labels=list(move.template), degeneracy=1, return_all_kinetics=False)
    rxn.kinetics = kinetics
    out = {"kind": "RMG family estimate", "available": True, "source": str(source),
           "entry": str(entry) if entry else None, "rank": getattr(entry, "rank", None),
           "family_template_direction": "forward"}
    if not family.auto_generated:
        out.update(cmp._rate_rule_source(family, rxn))
    return out


def _record_pair(fam, family_name, reactants, pathway=None, products=None, member=0):
    """A forward (as_generated, RMG family estimate) and reversed
    (reference-thermo reverse) EventRecord of one channel member on
    ``reactants``, shaped as the compiler publishes them: production rate
    source, heavy-atom ``atom_map``, reciprocal ``reverse_of`` links."""
    from rmgpy.kmc.compiler import LINK_PLACEHOLDER
    family = fam[family_name]
    cs = chm.build_channels(family, reactants)
    chs = [c for c in cs.channels if (pathway is None or c.pathway == tuple(pathway))]
    ch = chs[0]
    move = ch.representatives[member]
    prods = [rm.molecule_from_rows(p) for p in rm.split_rows(move.product_rows())]
    rg = [m.to_adjacency_list() for m in reactants]
    pg = [m.to_adjacency_list() for m in prods]
    rpos, ppos = _heavy_positions(reactants), _heavy_positions(prods)
    amap = {rpos[a]: ppos[a] for a in rpos}
    common = dict(family=family_name, template=";".join(ch.pathway))
    forward = EventRecord(**common, arity=len(reactants),
                          participant_site_types=["s"] * len(reactants),
                          reactant_multiplicities=[1] * len(reactants),
                          orientation="as_generated",
                          rate_source=production_rate_source(family, cs, move),
                          reactant_graphs=rg, product_graphs=pg, atom_map=amap,
                          reverse_of=LINK_PLACEHOLDER)
    reverse = EventRecord(**common, arity=len(prods), participant_site_types=["s"] * len(prods),
                          reactant_multiplicities=[1] * len(prods), orientation="reversed",
                          rate_source={"kind": "reference-thermo reverse", "available": True,
                                       "reference": "k_family / Kc"},
                          reactant_graphs=pg, product_graphs=rg,
                          atom_map={v: k for k, v in amap.items()},
                          reverse_of=forward.event_id)
    forward = replace(forward, reverse_of=reverse.event_id)
    return forward, reverse, cs, ch


def test_event_record_certificate_reference_is_additive():
    """Old records (no field) load and keep their event ids; the reference is
    serialized, round-trips, is outside the event hash and is validated."""
    r = EventRecord(family="H_Abstraction", template="a;b", arity=0)
    d = r.to_dict()
    assert "rate_certificate" not in d
    assert EventRecord.from_dict(dict(d)).event_id == r.event_id
    ref = {"schema": cm.REFERENCE_SCHEMA, "sha256": "0" * 64, "store": "artifact.rate_certificates",
           "owner_event_id": r.event_id, "template_side": "reactant_graphs", "padded": False,
           "complete": True, "eligible": True, "bound": True, "certified": True,
           "source_kinds": ["rate_rule"], "whole_species": False}
    r2 = replace(r, rate_certificate=ref)
    assert r2.event_id == r.event_id
    r2.validate()
    d2 = r2.to_dict()
    assert d2["rate_certificate"] == ref
    back = EventRecord.from_json(r2.to_json())
    assert back == r2 and back.event_id == r.event_id
    with pytest.raises(ValueError):
        replace(r, rate_certificate={"schema": cm.REFERENCE_SCHEMA, "sha256": "xyz"}).validate()


def test_attach_certificates_from_the_records_rate_source(fam):
    """attach_rate_certificates certifies the as_generated record on its
    reactants and the reversed record on its products (the family-template
    reactants); references resolve into the store; event ids unchanged; a
    library record is listed as not an RMG selection."""
    forward, reverse, _cs, _ch = _record_pair(fam, "H_Abstraction",
                                              [mol("[CH3]", 0), mol("CC", 10)])
    library = replace(forward, rate_source={"kind": "kMC kinetics library", "available": True},
                      event_id="")
    records, store = cm.attach_rate_certificates([forward, reverse, library], fam)
    assert [r.event_id for r in records] == [forward.event_id, reverse.event_id, library.event_id]
    f, r, lib = (store[x.rate_certificate["sha256"]] for x in records)
    assert f["template_side"] == "reactant_graphs" and r["template_side"] == "product_graphs"
    assert f["complete"] and r["complete"], (f["gaps"], r["gaps"])
    assert f["role"] == "rate_source" and r["role"] == "reverse_of_rate_source"
    assert records[0].rate_certificate["whole_species"]
    assert cm.selection_view(f["certificate"]) == cm.selection_view(r["certificate"])
    assert lib["source_kind"] == "not_rmg_selection" and not lib["complete"]
    for x in records:
        x.validate()
        assert cm.certificate_sha256(store[x.rate_certificate["sha256"]]) == \
            x.rate_certificate["sha256"]


def test_compiler_attaches_certificates_when_asked():
    """EventSetCompiler(certify_rates=True) attaches a reference to every
    record and publishes the store; ids and records are those of a compile
    without certificates."""
    plain, _db, _p = ct._fake_compiler(["pathway_A", "pathway_B"])
    base = plain.compile()
    certified, _db, _p = ct._fake_compiler(["pathway_A", "pathway_B"])
    certified.certify_rates = True
    art = certified.compile()
    assert [r["event_id"] for r in art["records"]] == [r["event_id"] for r in base["records"]]
    assert "rate_certificates" not in base
    for r in art["records"]:
        assert r["rate_certificate"]["sha256"] in art["rate_certificates"]


# --------------------------------------------------------------------------
# item 4: recertification
# --------------------------------------------------------------------------
def test_recertify_through_a_renumbered_extension(fam, certifier):
    """1-phenylethyl* abstracting a head-CH3 H of closed PS DP3: the record is
    recertified on the witness extended by one repeat unit at the far chain
    end, every id renumbered.  Compared through the atom map the selection is
    the same (raw ids would differ)."""
    closed = mol(ct.ps(3), 17)
    pe = mol("C[CH]c1ccccc1", 0)
    family = fam["H_Abstraction"]
    cs = chm.build_channels(family, [pe, closed])
    head = next(c for c in cs.channels if c.pathway[0].startswith("C/H3/Cs")
                and c.pathway[1].startswith("C_rad/H/CbCs"))
    forward, _rev, _cs, _ch = _record_pair(fam, "H_Abstraction", [pe, closed], head.pathway)
    rc = certifier.certify_record(forward)
    assert rc["complete"], rc["gaps"]
    _side, mols, _others = certifier.template_side(forward)
    reacting = {a for ch in rc["certificate"]["channels"] for m in ch["members"]
                for mp in m["mappings"] for _s, _l, a in mp["labels"]}
    ext, amap, info = cm.extend_witness(mols, reacting)
    assert any(i["extended"] for i in info)
    assert all(v != k for k, v in amap.items())
    res = certifier.recertify(forward, ext, amap, certificate=rc)
    assert res["certificate"]["complete"], res["certificate"]["gaps"]
    assert res["same_selection"], res["flips"]
    assert res["flips"] == []
    assert cm.selection_view(rc["certificate"]) != cm.selection_view(res["certificate"])


def test_recertify_reports_a_changed_selection(fam, certifier):
    """CH3* + propane, secondary H (leaf C/H2/Cs\\H3/Cs\\H3), recertified on
    CH3* + butane made by replacing a far methyl H (not a reacting atom) by a
    methyl: ids kept through the map, but a predicate of the donor descent
    flips and the leaf changes, so the selection is not the same."""
    family = fam["H_Abstraction"]
    propane = mol("CCC", 10)
    cs = chm.build_channels(family, [mol("[CH3]", 0), propane])
    ch = only(cs.find(pathway=("C/H2/Cs\\H3/Cs\\H3", "C_methyl")))
    cert = certifier.certify_channel(family, cs, ch)
    assert cert["complete"], cert["gaps"]
    rows = rm.molecule_rows(propane)
    rows.pop(18)
    rows[12] = rows[12][:5] + (((11, 1.0), (19, 1.0), (20, 1.0), (21, 1.0)),)
    rows[21] = ("C", None, 0, 0, 0, ((12, 1.0), (22, 1.0), (23, 1.0), (24, 1.0)))
    for h in (22, 23, 24):
        rows[h] = ("H", None, 0, 0, 0, ((21, 1.0),))
    butane = rm.molecule_from_rows(rows)
    amap = {a.id: a.id for m in cs.reactants for a in m.atoms if a.id != 18}
    res = certifier.recertify_certificate(family, cert, [mol("[CH3]", 0), butane], amap)
    assert res["certificate"]["complete"], res["certificate"]["gaps"]
    assert not res["same_selection"] and res["changed_mappings"]
    assert res["flips"]
    print("RECERT changed:", sorted({(f["node"], f["old"], f["new"]) for f in res["flips"]}))


def test_recertify_a_reversed_record_on_its_template_reactants(fam, certifier):
    """A reversed record is certified on its products (the family-template
    reactants); recertification takes the extension of that side."""
    pe = mol("C[CH]c1ccccc1", 0)
    closed = mol(ct.ps(3), 17)
    family = fam["H_Abstraction"]
    cs = chm.build_channels(family, [pe, closed])
    head = next(c for c in cs.channels if c.pathway[0].startswith("C/H3/Cs")
                and c.pathway[1].startswith("C_rad/H/CbCs"))
    _f, reverse, _cs, _ch = _record_pair(fam, "H_Abstraction", [pe, closed], head.pathway)
    rc = certifier.certify_record(reverse)
    assert rc["template_side"] == "product_graphs" and rc["complete"], rc["gaps"]
    _side, mols, _o = certifier.template_side(reverse)
    reacting = {a for ch in rc["certificate"]["channels"] for m in ch["members"]
                for mp in m["mappings"] for _s, _l, a in mp["labels"]}
    ext, amap, _info = cm.extend_witness(mols, reacting)
    res = certifier.recertify(reverse, ext, amap, certificate=rc)
    assert res["same_selection"], res["flips"]


# --------------------------------------------------------------------------
# required supported cases
# --------------------------------------------------------------------------
def _ring(mol_):
    return {a.id for a in mol_.atoms if mol_.is_atom_in_cycle(a)}


def test_required_supported_cases(fam, certifier):
    """CH3* + C2H6; H abstraction by a benzylic radical from a closed chain;
    chain-end recombination of two chains; beta-scission of a mid-chain
    radical: each is spin supported, every species is in the v1 domain, and
    its certificate is complete (beta-scission on the family-template side,
    its products)."""
    rows = []
    # CH3* + C2H6
    family = fam["H_Abstraction"]
    cs = chm.build_channels(family, [mol("[CH3]", 0), mol("CC", 10)])
    rows.append(("CH3* + C2H6", family, cs, only(cs.channels)))
    # benzylic radical abstracting from a closed chain (head CH3 donor)
    pe, closed = mol("C[CH]c1ccccc1", 0), mol(ct.ps(3), 100)
    cs = chm.build_channels(family, [pe, closed])
    ch = next(c for c in cs.channels if c.pathway[0].startswith("C/H3/Cs")
              and c.pathway[1].startswith("C_rad/H/CbCs"))
    rows.append(("benzylic H abstraction", family, cs, ch))
    # chain-end recombination of two chains
    rfam = fam["R_Recombination"]
    a = mol("CC(c1ccccc1)C[CH]c1ccccc1", 0)
    b = mol("[CH2]C(c1ccccc1)CC(c1ccccc1)C", 100)
    cs = chm.build_channels(rfam, [a, b])
    rings = _ring(a) | _ring(b)
    ch = next(c for c in cs.channels
              if all(x not in rings for _s, _l, x in c.representatives[0].labels))
    rows.append(("chain-end recombination", rfam, cs, ch))
    for name, fam_, cs, ch in rows:
        assert ch.support == chm.SUPPORTED, (name, ch.support)
        cert = certifier.certify_channel(fam_, cs, ch)
        assert cert["complete"], (name, cert["gaps"])
        assert all(v["domain"]["accepted"] for v in cert["species"].values()), name
        assert all(v["spin"]["status"] == "supported" for v in cert["species"].values()), name
        print("REQUIRED", name, [c["source_kind"] for c in cert["channels"]],
              [c["support_scope"] for c in cert["channels"]])
    # beta-scission of a mid-chain benzylic radical (R_Addition_MultipleBond reverse)
    afam = fam["R_Addition_MultipleBond"]
    radical = mol("CC(c1ccccc1)C[C](c1ccccc1)CC(c1ccccc1)C", 0)
    cs = chm.build_channels(afam, [radical], "reverse")
    ring = _ring(radical)
    # C-C scission of the backbone (C-H scission makes an internal alkene,
    # outside the declared domain: reported in the acceptance tables)
    backbone = [c for c in cs.channels
                if all(x not in ring for _s, _l, x in c.representatives[0].labels)
                and "[H]" not in c.products]
    assert backbone
    for ch in backbone:
        assert ch.support == chm.SUPPORTED
        for cert, mols in certifier.certify_reverse_channel(afam, cs, ch):
            assert cert["complete"], cert["gaps"]
            assert all(dom.check(m).accepted for m in mols)
        print("REQUIRED beta-scission", ch.pathway, ch.products,
              [c["source_kind"] for c in cert["channels"]])


# --------------------------------------------------------------------------
# item 5: the declared v1 domain
# --------------------------------------------------------------------------
def test_domain_accepts_the_ps_references_and_declared_ends():
    radical, closed, h_atom = ct.dp3_molecules()
    accepted = {
        "DP3 radical": radical, "DP3 closed": closed, "H": h_atom,
        "DP8 head radical": mol("[CH2]C(c1ccccc1)" + ct.ps(7), 0),
        "DP4": mol(ct.ps(4), 0), "1-phenylethyl": mol("C[CH]c1ccccc1", 0),
        "DP8 closed": mol(ct.ps(8), 0), "CH3": mol("[CH3]", 0), "C2H6": mol("CC", 0),
        "CH4": mol("C", 0), "H2": mol("[H][H]", 0),
        "vinylidene end": mol("C=C(c1ccccc1)CC(c1ccccc1)C", 0),
        "styryl end": mol("c1ccccc1C=CCC(c1ccccc1)C", 0),
        "Kekule DP2": mol("CC(C1=CC=CC=C1)CC(C1=CC=CC=C1)", 0),
        "two far benzylic radicals": mol("C[C](c1ccccc1)" + ct.ps(2) + "C[C](c1ccccc1)C", 0),
        "secondary backbone radical": mol("CC(c1ccccc1)[CH]C(c1ccccc1)C", 0),
    }
    for name, m in accepted.items():
        v = dom.check(m)
        assert v.accepted, (name, v.reason)
        print("DOMAIN", name, v.kind, v.details.get("ends"), v.details.get("radical_sites"))
    assert dom.check(accepted["vinylidene end"]).details["ends"] == ["methyl", "vinylidene"]
    assert "styryl" in dom.check(accepted["styryl end"]).details["ends"]


def test_domain_refuses_everything_else():
    refused = {
        "branched": "CC(C)(c1ccccc1)CC",
        "aryl radical on a chain": "CC(c1ccc[c]c1)CC(c1ccccc1)C",
        "phenyl radical": "[c]1ccccc1",
        "oxygen": "CC(O)C",
        "internal unsaturation": "CC(c1ccccc1)=CC(c1ccccc1)C",
        "radicals two bonds apart": "C[C](c1ccccc1)C[CH]C",
        "three radicals": "[CH2]C(c1ccccc1)" + ct.ps(3) + "C[C](c1ccccc1)" + ct.ps(3) +
                          "C[CH]c1ccccc1",
        "allylic radical": "C=C(c1ccccc1)[CH]C(c1ccccc1)C",
        "naphthyl": "CC(c1ccc2ccccc2c1)C",
        "carbene": "[CH2]",
        "two phenyls on one carbon": "CC(c1ccccc1)(c1ccccc1)C",
    }
    for name, smiles in refused.items():
        v = dom.check(mol(smiles, 0))
        assert not v.accepted and v.reason, name
        print("REFUSED", name, "|", v.reason)


# --------------------------------------------------------------------------
# rework 1: ownership, certification mode, production binding, the hang,
# cache semantics, nonvacuous recertification, mapping identity, support,
# averaged contributors, domain eligibility, bounded caches
# --------------------------------------------------------------------------
def _certified_fake_artifact():
    compiler, _db, _p = ct._fake_compiler(["pathway_A", "pathway_B"])
    compiler.certify_rates = True
    return compiler.compile()


def test_certificate_references_are_owned_by_their_record():
    """Validation binds every stored certificate to its owning record and
    derives the reference summary from the stored certificate: swapped
    references, a summary contradicting the certificate, a corrupted body and
    a re-hashed body naming another owner each refuse the artifact."""
    from rmgpy.kmc.compiler import validate_artifact
    art = _certified_fake_artifact()
    validate_artifact(art)
    first, second = art["records"][0], art["records"][1]
    assert first["rate_certificate"]["sha256"] != second["rate_certificate"]["sha256"]
    assert first["rate_certificate"]["owner_event_id"] == first["event_id"]
    swapped = copy.deepcopy(art)
    r0, r1 = swapped["records"][0], swapped["records"][1]
    r0["rate_certificate"], r1["rate_certificate"] = r1["rate_certificate"], r0["rate_certificate"]
    with pytest.raises(ValueError, match="owner"):
        validate_artifact(swapped)
    lying = copy.deepcopy(art)
    ref = lying["records"][0]["rate_certificate"]
    ref["whole_species"] = not ref["whole_species"]
    with pytest.raises(ValueError, match="summary"):
        validate_artifact(lying)
    corrupt = copy.deepcopy(art)
    sha = corrupt["records"][0]["rate_certificate"]["sha256"]
    corrupt["rate_certificates"][sha]["gaps"].append("tampered after hashing")
    with pytest.raises(ValueError, match="digest"):
        validate_artifact(corrupt)
    stolen = copy.deepcopy(art)
    ref = stolen["records"][0]["rate_certificate"]
    body = stolen["rate_certificates"].pop(ref["sha256"])
    body["owner"]["event_id"] = stolen["records"][1]["event_id"]
    stolen["rate_certificates"][cm.certificate_sha256(body)] = body
    ref.update(cm.certificate_reference(body))
    with pytest.raises(ValueError, match="owner"):
        validate_artifact(stolen)


def test_certification_mode_requires_coverage_of_every_record():
    """A certified compile declares ``rate_certification`` (mode required,
    coverage every record); a removed reference, references without a
    declared mode, and an unreferenced stored certificate are refused.  An
    uncertified compile declares nothing and validates."""
    from rmgpy.kmc.compiler import validate_artifact
    art = _certified_fake_artifact()
    assert art["rate_certification"]["mode"] == "required"
    assert art["rate_certification"]["coverage"] == "every record"
    removed = copy.deepcopy(art)
    del removed["records"][1]["rate_certificate"]
    with pytest.raises(ValueError, match="coverage"):
        validate_artifact(removed)
    undeclared = copy.deepcopy(art)
    del undeclared["rate_certification"]
    with pytest.raises(ValueError, match="certification mode"):
        validate_artifact(undeclared)
    orphan = copy.deepcopy(art)
    body = {"orphan": True}
    orphan["rate_certificates"][cm.certificate_sha256(body)] = body
    with pytest.raises(ValueError, match="unreferenced"):
        validate_artifact(orphan)
    plain, _db, _p = ct._fake_compiler(["pathway_A"])
    base = plain.compile()
    assert "rate_certification" not in base
    validate_artifact(base)


def test_record_certificate_binds_the_production_selection(fam, certifier):
    """1-phenylethyl* abstracting a head-CH3 H of closed PS DP3 (an averaged
    rule: its training contributor and weight; depository excluded by size).
    The record certificate names its owner, the members its heavy-atom map
    selects and RMG's production selection, which must equal the record's
    serialized rate source.  A changed weight, a changed source or entry, or
    an atom map of another event refuses the binding (the rate is never
    silently re-estimated)."""
    reactants = [mol("C[CH]c1ccccc1", 0), mol(ct.ps(3), 17)]
    forward, reverse, _cs, ch = _record_pair(fam, "H_Abstraction", reactants,
                                             ("C/H3/Cs", "C_rad/H/CbCs"))
    rc = certifier.certify_record(forward)
    assert rc["owner"]["event_id"] == forward.event_id
    assert rc["binding"]["status"] == "bound", rc["binding"]
    assert rc["complete"] and rc["eligible"] and rc["certified"], rc["gaps"]
    # the heavy-atom map pins the donor carbon; the three methyl H members
    # (equivalent up to H identity) all correspond to the record's event
    assert 1 <= len(rc["binding"]["members"]) <= ch.size
    assert forward.rate_source["training"] or forward.rate_source["rules"]
    tampered = copy.deepcopy(forward.rate_source)
    (tampered["training"] or tampered["rules"])[0]["weight"] += 0.25
    bad = certifier.certify_record(replace(forward, rate_source=tampered, event_id=""))
    assert bad["binding"]["status"] == "refused" and not bad["certified"]
    assert any("weight" in m or "training" in m or "rules" in m
               for m in bad["binding"]["mismatches"]), bad["binding"]
    moved = dict(forward.rate_source, source="Kinetics Depository H_Abstraction/training",
                 entry="C3H8 + CH3")
    other_source = certifier.certify_record(replace(forward, rate_source=moved, event_id=""))
    assert other_source["binding"]["status"] == "refused"
    amap = dict(forward.atom_map)
    amap[0], amap[1] = amap[1], amap[0]     # the methyl radical carbon ends up in the propyl
    foreign = certifier.certify_record(replace(forward, atom_map=amap, event_id=""))
    assert foreign["binding"]["status"] == "refused" and not foreign["certified"]
    assert any("member" in m for m in foreign["binding"]["mismatches"]), foreign["binding"]
    owned = certifier.certify_record(reverse, rate_owner=forward)
    assert owned["binding"]["status"] == "bound", owned["binding"]
    assert owned["owner"]["rate_owner_event_id"] == forward.event_id
    unowned = certifier.certify_record(reverse)
    assert unowned["binding"]["status"] == "refused" and not unowned["certified"]


def test_library_reverse_rates_are_not_rmg_selections(fam, certifier):
    """A reference-thermo reverse of a library rate (k_library / Kc) has no
    RMG provenance: neither by its own rate source nor through its
    rate-owning record."""
    from rmgpy.kmc.compiler import LINK_PLACEHOLDER
    forward, reverse, _cs, _ch = _record_pair(fam, "H_Abstraction", [mol("[CH3]", 0), mol("CC", 10)])
    library_reverse = replace(reverse, rate_source={
        "kind": "reference-thermo reverse", "available": True, "reference": "k_library / Kc",
        "forward_rate_source": {"kind": "kMC kinetics library", "available": True}}, event_id="")
    rc = certifier.certify_record(library_reverse)
    assert rc["certificate"] is None and rc["source_kind"] == "not_rmg_selection"
    assert any("library" in g for g in rc["gaps"])
    library_forward = replace(forward, rate_source={"kind": "kMC kinetics library",
                                                    "available": True},
                              reverse_of=LINK_PLACEHOLDER, event_id="")
    linked = replace(reverse, reverse_of=library_forward.event_id)
    rc = certifier.certify_record(linked, rate_owner=library_forward)
    assert rc["certificate"] is None and rc["source_kind"] == "not_rmg_selection"


def test_tree_selection_refuses_a_nonterminating_walk_before_calling_rmg(monkeypatch):
    """The certifier's own tree selection checks RMG's walk termination
    BEFORE get_kinetics_for_template (whose loop would spin forever)."""
    from types import SimpleNamespace
    calls = []

    class Family:
        auto_generated = True
        label = "fake"

        def retrieve_template(self, labels):
            return labels

        def get_kinetics_for_template(self, *args, **kwargs):
            calls.append(args)
            raise AssertionError("RMG's non-terminating walk was entered")

    leaf = SimpleNamespace(label="leaf")
    monkeypatch.setattr(chm, "autogenerated_walk", lambda family, pathway: {
        "matched": leaf, "supplying": leaf, "rmg_terminates": False,
        "walk": [{"from": "leaf", "parent": "root", "moved": False}]})
    with pytest.raises(chm.UnsupportedSelection, match="does not terminate"):
        cm.Certifier({})._tree_selection(Family(), SimpleNamespace(pathway=("leaf",)))
    assert calls == []


def test_descent_cache_hits_carry_predicate_semantics(fam):
    """A second certificate built from cached descents carries the semantics
    of every node its descents reference; completeness refuses a missing
    one; the cache key keeps the template order of the participants."""
    certifier = cm.Certifier(fam)
    family = fam["H_Abstraction"]
    cs = chm.build_channels(family, [mol("[CH3]", 0), mol("CCC", 10)])
    ch = only(cs.find(pathway=("C/H3/Cs\\H2\\Cs", "C_methyl")))
    first = certifier.certify_channel(family, cs, ch)
    second = certifier.certify_channel(family, cs, ch)
    for cert in (first, second):
        assert cert["complete"], cert["gaps"]
        referenced = {p["node"] for v in cert["descents"].values() for d in v
                      for p in cm.iter_descent_predicates(d)}
        assert referenced and referenced <= set(cert["nodes"]), referenced - set(cert["nodes"])
    assert second["nodes"] == first["nodes"]
    broken = copy.deepcopy(second)
    broken["nodes"].pop("C/H3/Cs\\H2\\Cs")
    assert any("semantics" in g for g in cm.completeness_gaps(broken))
    a, b = mol("[CH3]", 0), mol("CCC", 10)
    assert certifier._descent_key(family, [a, b]) != certifier._descent_key(family, [b, a])


def test_recertification_is_never_vacuous(fam, certifier):
    """No certified mapping, an incomplete or a corrupt-map recertification
    is refused (same_selection False with the reason), never vacuously
    true; the identity map of a complete certificate is a fixed point."""
    family = fam["H_Abstraction"]
    reactants = [mol("[CH3]", 0), mol("CC", 10)]
    ident = {a.id: a.id for m in reactants for a in m.atoms}
    empty = certifier.certify(family, reactants, pathway=("no such node",))
    assert not empty["complete"] and empty["channels"] == []
    res = certifier.recertify_certificate(family, empty, reactants, ident)
    assert res["same_selection"] is False and res["refusals"]
    cs = chm.build_channels(family, reactants)
    cert = certifier.certify_channel(family, cs, only(cs.channels))
    ok = certifier.recertify_certificate(family, cert, list(cs.reactants), ident)
    assert ok["same_selection"] and ok["fixed_point"], ok
    incomplete = copy.deepcopy(cert)
    incomplete["complete"], incomplete["gaps"] = False, ["incomplete on purpose"]
    res = certifier.recertify_certificate(family, incomplete, list(cs.reactants), ident)
    assert res["same_selection"] is False and any("complete" in r for r in res["refusals"])
    carbon = next(a.id for a in reactants[1].atoms if a.is_carbon())
    hydrogen = next(a.id for a in reactants[1].atoms if a.is_hydrogen())
    swapped = dict(ident)
    swapped[carbon], swapped[hydrogen] = ident[hydrogen], ident[carbon]
    res = certifier.recertify_certificate(family, cert, list(cs.reactants), swapped)
    assert res["same_selection"] is False and any("atom map" in r for r in res["refusals"])
    dropped = {k: v for k, v in ident.items() if k != carbon}
    res = certifier.recertify_certificate(family, cert, list(cs.reactants), dropped)
    assert res["same_selection"] is False and any("atom map" in r for r in res["refusals"])


def test_mapping_identity_keeps_template_slots(fam, certifier):
    """CH3* + CH3*: the two opposite '*' assignments are two mappings of one
    member; the mapping view keeps both (slot, mapped atoms, edits) instead
    of collapsing them to one label set."""
    family = fam["R_Recombination"]
    cs = chm.build_channels(family, [mol("[CH3]", 0), mol("[CH3]", 10)])
    ch = only(cs.channels)
    cert = certifier.certify_channel(family, cs, ch)
    assert cert["complete"], cert["gaps"]
    (member,) = cert["channels"][0]["members"]
    assert len(member["mappings"]) == 2
    assert len(cm.mapping_view(cert)) == 2
    swap = {a.id: (a.id + 10 if a.id < 10 else a.id - 10) for m in cs.reactants for a in m.atoms}
    res = certifier.recertify_certificate(family, cert, list(cs.reactants), swap)
    assert res["same_selection"], res


def test_equal_selection_is_not_equal_support(fam, certifier):
    """A nested logic component that changes while the top verdict stays
    changes the selection signature.  A recertified selection whose
    transported support differs (here: one ball row gains a bond) is the
    same selection but not a fixed point."""
    d1 = {"top": "T", "structure": "s", "leaf": "T", "steps": [],
          "root": {"node": "T", "result": "matched", "logic": "LogicOr", "invert": False,
                   "components": [{"node": "A", "result": "failed", "logic": None},
                                  {"node": "B", "result": "matched", "logic": None}]}}
    d2 = copy.deepcopy(d1)
    d2["root"]["components"] = [{"node": "A", "result": "matched", "logic": None}]
    assert cm.descent_signature([d1]) != cm.descent_signature([d2])
    closed, pe = mol(ct.ps(3), 17), mol("C[CH]c1ccccc1", 0)
    family = fam["H_Abstraction"]
    cs = chm.build_channels(family, [pe, closed])
    head = next(c for c in cs.channels if c.pathway[0].startswith("C/H3/Cs")
                and c.pathway[1].startswith("C_rad/H/CbCs"))
    forward, _rev, _cs, _ch = _record_pair(fam, "H_Abstraction", [pe, closed], head.pathway)
    rc = certifier.certify_record(forward)
    _side, mols, _others = certifier.template_side(forward)
    ext, amap, _info = cm.extend_witness(mols, _reacting(rc["certificate"]))
    res = certifier.recertify(forward, ext, amap, certificate=rc)
    assert res["same_selection"], res
    assert res["fixed_point"] == (res["same_selection"] and res["same_support"])
    print("SUPPORT renumbered extension same_support", res["same_support"],
          res.get("support_changes", [])[:3])
    new = copy.deepcopy(res["certificate"])
    d = next(x for v in new["descents"].values() for x in v if x.get("atoms"))
    d["atoms"][0][8].append([987654321, 1.0])
    out = cm.compare_certificates(rc["certificate"], new, amap)
    assert out["same_selection"] and not out["same_support"] and not out["fixed_point"]


def test_averaged_rule_contributors_are_the_production_extraction(fam, certifier):
    """An averaged rule (C/H2/NonDeC;C_rad/H/CbCs) is certified with its
    underlying rules and training entries and their weights, exactly as the
    compiler extracts them from RMG's estimate (compiler._rate_rule_source)."""
    from rmgpy.kmc import compiler as cmp
    family = fam["H_Abstraction"]
    cs = chm.build_channels(family, [mol("C[CH]c1ccccc1", 0), mol(ct.ps(3), 17)])
    ch = next(c for c in cs.channels if c.pathway == ("C/H2/NonDeC", "C_rad/H/CbCs"))
    cert = certifier.certify_channel(family, cs, ch)
    (c,) = cert["channels"]
    assert c["source_kind"] == "averaged"
    rxn, _r, _p = chm.template_reaction(family, cs, ch.representatives[0])
    rxn.kinetics = family.get_kinetics(rxn, template_labels=list(ch.pathway), degeneracy=1,
                                       return_all_kinetics=False)[0]
    expected = cmp._rate_rule_source(family, rxn)
    got = c["selection"]["sources"]
    assert got["rules"] == expected["rules"] and got["training"] == expected["training"]
    assert got["exact"] == expected["exact"]
    assert len(got["rules"]) + len(got["training"]) >= 2


def test_domain_checks_the_multiplicity_convention_before_normalizing(fam, certifier):
    """A separated biradical given with multiplicity one violates the declared
    maximum-multiplicity convention: refused by the domain, and a certificate
    of it is not eligible (it would otherwise be rebuilt at maximum
    multiplicity, possibly changing the selection)."""
    m = mol("C[C](c1ccccc1)" + ct.ps(2) + "C[C](c1ccccc1)C", 0)
    assert dom.check(m).accepted
    singlet = m.copy(deep=True)
    singlet.multiplicity = 1
    v = dom.check(singlet)
    assert not v.accepted and "multiplicity" in v.reason
    cert = certifier.certify(fam["R_Recombination"], [singlet, mol("[H]", 500)])
    assert cert["input_refusals"] and not cert["eligible"]


def test_domain_styryl_end_checks_both_alkene_carbons():
    """The declared styryl end is PhCH=CH-: PhCH=C(Ph)-CH2-CH3 is refused;
    the declared styryl and vinylidene ends and styrene stay accepted."""
    v = dom.check(mol("c1ccccc1C=C(c1ccccc1)CC", 0))
    assert not v.accepted, v
    assert "styryl" in dom.check(mol("c1ccccc1C=CCC(c1ccccc1)C", 0)).details["ends"]
    assert dom.check(mol("C=C(c1ccccc1)CC(c1ccccc1)C", 0)).accepted
    assert dom.check(mol("C=Cc1ccccc1", 0)).accepted


def test_domain_refusal_is_a_channel_disposition(fam, certifier):
    """CH3* + ethylbenzene: the ring-H-loss channels (Cb_H) are spin/pathway
    supported but their aryl-radical product is outside the v1 domain, so
    the certificate's channel eligibility is refused, the channel is not
    publishable and the certificate is not eligible; recertification refuses
    it.  The chain channels stay eligible."""
    family = fam["H_Abstraction"]
    cs = chm.build_channels(family, [mol("[CH3]", 0), mol("CCc1ccccc1", 10)])
    ring = [c for c in cs.channels if c.pathway[0] == "Cb_H"]
    chain = [c for c in cs.channels if c.pathway[0] != "Cb_H"]
    assert ring and chain
    for ch in ring:
        assert ch.support == chm.SUPPORTED
        cert = certifier.certify_channel(family, cs, ch)
        (c,) = cert["channels"]
        assert c["eligibility"]["disposition"] == "refused:domain", c["eligibility"]
        assert not c["publishable"] and not cert["eligible"]
    for ch in chain:
        cert = certifier.certify_channel(family, cs, ch)
        assert cert["eligible"] and cert["channels"][0]["publishable"], \
            cert["channels"][0]["eligibility"]
    ident = {a.id: a.id for m in cs.reactants for a in m.atoms}
    cert = certifier.certify_channel(family, cs, ring[0])
    res = certifier.recertify_certificate(family, cert, list(cs.reactants), ident)
    assert res["same_selection"] is False and any("domain" in r for r in res["refusals"])


def test_certifier_caches_are_bounded(fam):
    """Every certifier cache is bounded; resonance forms live with their
    channel set and leave with it."""
    certifier = cm.Certifier(fam, max_channel_sets=2, max_entries=8)
    family = fam["H_Abstraction"]
    for smiles in ("CC", "CCC", "CCCC", "CC(C)C"):
        cs = certifier.channel_set(family, [mol("[CH3]", 0), mol(smiles, 10)])
        certifier.certify_channel(family, cs, cs.channels[0])
    sizes = certifier.cache_sizes()
    assert sizes["channel_sets"] <= 2, sizes
    assert sizes["forms"] <= 2, sizes
    assert all(v <= 8 for k, v in sizes.items() if k not in ("channel_sets", "forms")), sizes


# --------------------------------------------------------------------------
# acceptance (slow): the PS references and the closed DP8 compile
# --------------------------------------------------------------------------
def _old_spin(mols):
    return "unsupported:spin" if any(chm._radical_count(m) >= 2 for m in mols) else "supported"


def _species_of(channel_set, ch):
    move = ch.representatives[0]
    return list(channel_set.reactants) + [rm.molecule_from_rows(p)
                                          for p in rm.split_rows(move.product_rows())]


def _reacting(cert):
    return {a for ch in cert["channels"] for m in ch.get("members", [])
            for mp in m["mappings"] for _s, _l, a in mp["labels"]}


def recert_row(certifier, family, cert, mols):
    """Recertification on +1 repeat unit: outcome ``same`` / ``changed`` /
    ``unmapped`` (a reacting atom was removed by the extension) / ``not
    extended``, keyed by the channel support, with the bond distance from the
    extension point to the nearest reacting atom."""
    support = ",".join(sorted({c.get("support", "?") for c in cert["channels"]})) or "none"
    ext, amap, info = cm.extend_witness(mols, _reacting(cert))
    if not any(i["extended"] for i in info):
        return {"extended": False, "outcome": f"{support}:not extended"}
    distance = min((i["distance_to_protected"] for i in info
                    if i.get("extended") and i.get("distance_to_protected") is not None),
                   default=None)
    new = certifier.recertify_certificate(family, cert, ext, amap)
    if new.get("certificate") is None:
        outcome = "unmapped"
    else:
        outcome = "same" if new["same_selection"] else "changed"
    return {"extended": True, "same": new["same_selection"], "distance": distance,
            "outcome": f"{support}:{outcome}",
            "flips": sorted({f["node"] for f in new["flips"]}),
            "new_complete": (new["certificate"] or {}).get("complete")}


# Rate selection runs in the family TEMPLATE direction only, so the forward
# channel sets certify what RMG reads on the references; reverse channels are
# certified on their products (certify_reverse_channel) for the unimolecular
# decompositions (homolysis, beta-scission).  Disproportionation reverse
# (alkene + H donor) is not certified on the references: its template-side
# witnesses are radical pairs that are not reference species, and its mappings
# on Kekule ring "alkenes" of two PS chains run > 10 min per pair.
REFERENCE_SETS = {"H_Abstraction": ("forward",), "Disproportionation": ("forward",),
                  "R_Recombination": ("forward", "reverse"),
                  "R_Addition_MultipleBond": ("forward", "reverse"),
                  "intra_H_migration": ("forward",)}


def reference_table(certifier, fam, label, mols, out):
    t0 = time.perf_counter()
    kinds, whole, incomplete = Counter(), set(), []
    spin_before, spin_after = Counter(), Counter()
    domain, recert, flips, ambiguous = {}, Counter(), Counter(), []
    changed_distance = Counter()
    n = 0
    keys = chm.SpeciesKeys()
    for name, directions in sorted(REFERENCE_SETS.items()):
        for direction in directions:
            family = fam[name]
            groups = [[m] for m in mols] + [[mols[i], mols[j]] for i in range(len(mols))
                                            for j in range(i + 1, len(mols))]
            for group in groups:
                t1 = time.perf_counter()
                cs = chm.build_channels(family, group, direction, species=keys)
                if not cs.channels:
                    continue
                for ch in cs.channels:
                    n += 1
                    species = _species_of(cs, ch)
                    spin_before[_old_spin(species)] += 1
                    spin_after[ch.support] += 1
                    if ch.support == chm.PATHWAY_AMBIGUOUS:
                        ambiguous.append((name, direction, ch.pathway, ch.ambiguous_with))
                    for m in species:
                        v = dom.check(m)
                        domain[keys(m)] = (v.accepted, v.reason)
                    if direction == "forward":
                        certs = [(certifier.certify_channel(family, cs, ch), list(cs.reactants))]
                    else:
                        certs = certifier.certify_reverse_channel(family, cs, ch)
                    for cert, tmols in certs:
                        for c in cert["channels"]:
                            kinds[c.get("source_kind")] += 1
                            if c.get("whole_species"):
                                whole.add((name, ";".join(c["pathway"])))
                        if not cert["complete"]:
                            incomplete.append((name, direction, ch.pathway, cert["gaps"][:3]))
                        row = recert_row(certifier, family, cert, tmols)
                        recert[row["outcome"]] += 1
                        if row["outcome"].endswith(":changed"):
                            changed_distance[(row["outcome"], row["distance"])] += 1
                            for node in row["flips"]:
                                flips[(name, node)] += 1
                print(f"PROGRESS {label} {name} {direction} {[len(m.atoms) for m in group]} "
                      f"channels {len(cs.channels)} {time.perf_counter() - t1:.0f} s", flush=True)
    refused = {k: r for k, (ok, r) in domain.items() if not ok}
    table = {"reference": label, "channels": n, "certificate_source_kinds": dict(kinds),
             "incomplete": len(incomplete), "whole_species": sorted(whole),
             "spin_before": dict(spin_before), "spin_after": dict(spin_after),
             "pathway_ambiguous": [list(map(str, x)) for x in ambiguous],
             "domain_accepted": len(domain) - len(refused),
             "domain_refused": sorted([k, r] for k, r in refused.items()),
             "recertification": dict(recert),
             "changed_by_distance": sorted([o, d, c] for (o, d), c in changed_distance.items()),
             "flips": sorted([f, node, c] for (f, node), c in flips.items()),
             "incomplete_detail": [list(map(str, x)) for x in incomplete],
             "seconds": round(time.perf_counter() - t0)}
    out.append(table)
    print("ACCEPTTABLE", json.dumps(table), flush=True)
    return table


@pytest.mark.slow
@pytest.mark.skipif(not SLOW, reason="set RMG_KMC_SLOW=1")
@pytest.mark.parametrize("label", ["DP3", "DP8 head radical", "DP4 + 1-phenylethyl"])
def test_reference_acceptance_tables(fam, certifier, label):
    """Per PS reference: every channel certified (forward channels on the
    reference, reverse channels on each member's family-template reactants),
    completeness, source kinds, whole-species list, domain verdict of every
    species, spin verdicts before/after the separated-radical policy, and
    recertification on +1 repeat unit.  The reference species themselves
    must be in the domain (gate)."""
    mols = ct.ps_references()[label]
    tables = []
    table = reference_table(certifier, fam, label, mols, tables)
    out = os.environ.get("RMG_KMC_I135_TABLES")
    if out:
        with open(out, "a") as stream:
            stream.write(json.dumps(table) + "\n")
    for m in mols:
        assert dom.check(m).accepted, (label, dom.check(m).reason)
    assert table["incomplete"] == 0


@pytest.mark.slow
@pytest.mark.skipif(not SLOW, reason="set RMG_KMC_SLOW=1")
def test_dp8_closed_compile_certificates(fam, certifier, tmp_path):
    """The closed PS DP8 reference compile with certify_rates: 132 records
    (unchanged by the duplicate-validation fix), every record carries a
    reference into the store; per record completeness, source kinds,
    whole-species records, domain and spin verdicts, recertification on +1
    repeat of the template side."""
    t0 = time.perf_counter()
    artifact = ct.compile_dp8_closed(tmp_path, certify_rates=True)
    t1 = time.perf_counter()
    records = artifact["records"]
    store = artifact["rate_certificates"]
    print(f"DP8CERT compile+certify {t1 - t0:.0f} s records {len(records)}")
    assert len(records) == 132
    kinds, whole, incomplete, recert = Counter(), [], [], Counter()
    domain, spin_before, spin_after, flips = {}, Counter(), Counter(), Counter()
    by_site = Counter()
    for record in records:
        ref = record["rate_certificate"]
        rc = store[ref["sha256"]]
        assert cm.certificate_sha256(rc) == ref["sha256"]
        by_site[(record["site_type"], rc.get("source_kind") or "rmg", rc["complete"])] += 1
        cert = rc.get("certificate")
        if cert is None:
            kinds[rc.get("source_kind")] += 1
            continue
        if record["site_type"] == "ref_dp8_closed" and not rc["complete"]:
            incomplete.append((record["event_id"][:16], rc["gaps"][:3]))
        for c in cert["channels"]:
            kinds[c.get("source_kind")] += 1
            if c.get("whole_species"):
                whole.append(record["event_id"][:16])
            spin_after[c["support"]] += 1
        _side, mols, others = certifier.template_side(record)
        spin_before[_old_spin(mols + others)] += 1
        keys = chm.SpeciesKeys()
        for m in mols + others:
            v = dom.check(m)
            domain[keys(m)] = (v.accepted, v.reason)
        row = recert_row(certifier, fam[record["family"]], cert, mols)
        recert[row["outcome"]] += 1
        if row["outcome"].endswith(":changed"):
            for node in row["flips"]:
                flips[node] += 1
    t2 = time.perf_counter()
    print(f"DP8CERT table {t2 - t1:.0f} s")
    print("DP8CERT by site/kind/complete", dict(by_site))
    print("DP8CERT source kinds", dict(kinds))
    print("DP8CERT whole_species", whole)
    print("DP8CERT spin before", dict(spin_before), "after", dict(spin_after))
    refused = {k: r for k, (ok, r) in domain.items() if not ok}
    print(f"DP8CERT domain {len(domain) - len(refused)} accepted {len(refused)} refused")
    for k, r in sorted(refused.items()):
        print("DP8CERT refused", k, "|", r)
    print("DP8CERT recertification", dict(recert))
    for node, c in sorted(flips.items()):
        print("DP8CERT flip", node, c)
    for x in incomplete:
        print("DP8CERT INCOMPLETE", x)
    out = os.environ.get("RMG_KMC_I135_TABLES")
    if out:
        with open(out, "a") as stream:
            stream.write(json.dumps({
                "reference": "DP8 closed compile", "records": len(records),
                "by_site_kind_complete": [[list(k), v] for k, v in sorted(by_site.items())],
                "certificate_source_kinds": dict(kinds), "whole_species_records": whole,
                "spin_before": dict(spin_before), "spin_after": dict(spin_after),
                "domain_accepted": len(domain) - len(refused),
                "domain_refused": sorted([k, r] for k, r in refused.items()),
                "recertification": dict(recert),
                "flips": sorted([n, c] for n, c in flips.items()),
                "incomplete": [list(map(str, x)) for x in incomplete],
                "compile_and_certify_s": round(t1 - t0), "table_s": round(t2 - t1),
                "store_bytes": len(cm.canonical_json(store))}) + "\n")
    assert not incomplete
    assert dom.check(ct.mol(ct.ps(8), 0)).accepted
