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
# depository exclusion decided by construction: the declared chain grows the
# witness at its own growth sites and every analogue is compared by
# Reaction.is_isomorphic (formula equality only selects the candidates)
# --------------------------------------------------------------------------
BACKBONE_SCISSION = ("Cds-HH_Cds-CsH", "CsJ-CsHH")
CH2 = (("C", 1), ("H", 2))
# sha256 of json.dumps(verdict, sort_keys=True) at the base (a8c5a9f91), no
# growth provider: /home/alon/runs/i208/probe/pin_probe.py
LEGACY_HABS_DIGEST = "500c30e0a7e0965425372065d7a4b04b9b012005ee3cf2d52705f756cf49ca1f"
LEGACY_SCISSION_PS_DIGEST = "093bf239be4d5a90d789cb907ed70e7bdba966c36eb90bbd1000c5f986b5b007"
LEGACY_SCISSION_CH2_DIGEST = "c0d8b5f2b849ca4f0aa72ed4df0d115d9534d3b74108d127880b0334032477ba"


@pytest.fixture(scope="module")
def pe_chain():
    import declared_parity_support as sup
    from rmgpy.kmc import declared_chain as dc
    return dc.DeclaredChainProvider(sup.pe_declaration())


def _digest(value):
    import hashlib
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def _backbone_scission(family):
    """The backbone beta-scission of a PE mid-chain radical in the family
    template direction (1-pentene + n-propyl -> oct-4-yl, the compiled PE
    witness): reactants and products with shared atom ids."""
    cs = chm.build_channels(family, [mol("C=CCCC", 0), mol("[CH2]CC", 20)])
    for ch in cs.find(pathway=BACKBONE_SCISSION):
        rep = ch.representatives[0]
        products = [rm.molecule_from_rows(p) for p in rm.split_rows(rep.product_rows())]
        if [p.to_smiles() for p in products] == ["CCC[CH]CCCC"]:
            return list(cs.reactants), products, cs, ch
    raise AssertionError("no backbone beta-scission channel on the witness")


def _training(family, indices):
    (training,) = [d for d in family.depositories if d.label.endswith("training")]
    found = [e for e in training.entries.values() if int(e.index) in indices]
    assert sorted(int(e.index) for e in found) == sorted(indices)
    return found


def _stub_family(family, entries, tag="stub"):
    """``family`` reduced to one depository holding only ``entries``."""
    from types import SimpleNamespace
    depository = SimpleNamespace(label=f"{family.label}/{tag}",
                                 entries={int(e.index): e for e in entries})
    return SimpleNamespace(label=f"{family.label}|{tag}", depositories=[depository])


def _reversed(entry):
    """The entry stored the other way round (products -> reactants)."""
    from types import SimpleNamespace
    from rmgpy.reaction import Reaction
    return SimpleNamespace(index=int(entry.index), label=entry.label, rank=entry.rank,
                           item=Reaction(reactants=list(entry.item.products),
                                         products=list(entry.item.reactants)))


def test_a_declared_chain_names_where_a_period_inserts(pe_chain):
    """The declared chain's growth sites are its capped artificial ends: a whole
    terminal period whose join is capped by hydrogen, read from the declared
    geometry.  Both methyl ends of oct-4-yl, the methyl end of n-propyl and of
    1-pentene; none on a species without such an end (methyl, ethylene, H; the
    radical and vinyl ends are physical).  Each site extends by the provider's
    own fragment.  A legacy provider names none: None, it cannot say."""
    from rmgpy.kmc import boundaries as b
    counts = {}
    for smiles in ("CCC[CH]CCCC", "[CH2]CC", "C=CCCC", "[CH3]", "C=C", "[H]"):
        molecule = mol(smiles, 0)
        by_id = {a.id: a for a in molecule.atoms}
        sites = pe_chain.growth_sites(molecule)
        counts[smiles] = len(sites)
        for port, coordinate in sites:
            atom = by_id[port.atom_id]
            assert atom.is_carbon() and atom.radical_electrons == 0
            assert sum(1 for n in atom.bonds if n.is_hydrogen()) == 3
            assert len(port.cap_atom_ids) == 1 and by_id[port.cap_atom_ids[0]] in atom.bonds
            assert port.artificial
            pe_chain.fragment(port, coordinate)     # the declared extension path
    assert counts == {"CCC[CH]CCCC": 2, "[CH2]CC": 1, "C=CCCC": 1, "[CH3]": 0, "C=C": 0,
                      "[H]": 0}
    for provider in b.LEGACY_PROVIDERS:
        assert provider.growth_sites(mol("CCC[CH]CCCC", 0)) is None


def test_growth_never_touches_a_reacting_atom(fam, pe_chain):
    """H abstraction from a methyl end of propane: the reacting methyl is no
    growth site of the reaction (its H leaves); the other methyl end is.  The
    methyl radical has no chain end."""
    family = fam["H_Abstraction"]
    cs = chm.build_channels(family, [mol("CCC", 0), mol("[CH3]", 20)])
    ch = only(cs.find(pathway=("C/H3/Cs\\H2\\Cs", "C_methyl")))
    for move in [m for m in cs.moves if m.template is not None
                 and tuple(m.template) == ch.pathway and m.products is not None]:
        products = [rm.molecule_from_rows(p) for p in rm.split_rows(move.product_rows())]
        sites = cm.reaction_growth_sites(pe_chain, list(cs.reactants), products)
        reacting = {a for a in rm.molecule_rows(cs.reactants[0])
                    if rm.molecule_rows(cs.reactants[0])[a] !=
                    {k: v for p in products for k, v in rm.molecule_rows(p).items()}.get(a)}
        assert len(sites) == 1
        ((species, port, _coordinate),) = sites
        assert species == 0 and port.atom_id not in reacting
        assert not set(port.cap_atom_ids) & reacting


def test_pe_backbone_beta_scission_certifies_except_at_its_matched_sizes(fam, pe_chain):
    """The formula argument refuses the PE backbone beta-scission (nine
    candidates, ``refused:unknown``).  Grown at its declared sites, the
    witness IS training 2915-2919 (n-alkyl + 1-alkene -> a linear dodecyl
    radical) at the sizes listed, and none of the formula-only candidates
    (456, 636 and the branched 653, 654).  Those are finite sizes RMG rates
    from the depository: ``excluded:except_matched`` names them, each match
    carrying the entry's kinetics as RMG reads them and the analogue itself,
    and the channel certifies."""
    from rmgpy.molecule import Molecule
    family = fam["R_Addition_MultipleBond"]
    reactants, products, cs, ch = _backbone_scission(family)
    legacy = cm.depository_exclusion(family, reactants, products, repeat=CH2)
    assert legacy["verdict"] == "refused:unknown"
    assert sorted(c["index"] for c in legacy["candidates"]) == \
        [456, 636, 653, 654, 2915, 2916, 2917, 2918, 2919]
    ex = cm.depository_exclusion(family, reactants, products, repeat=CH2, growth=pe_chain)
    assert ex["verdict"] == "excluded:except_matched", ex
    sizes = {m["index"]: (m["orientation"], m["species_heavy_atoms"]) for m in ex["matches"]}
    assert sizes == {2915: ("forward", {"reactants": [5, 7], "products": [12]}),
                     2916: ("forward", {"reactants": [3, 9], "products": [12]}),
                     2917: ("forward", {"reactants": [6, 6], "products": [12]}),
                     2918: ("forward", {"reactants": [4, 8], "products": [12]}),
                     2919: ("forward", {"reactants": [5, 7], "products": [12]})}, ex["matches"]
    label = "R_Addition_MultipleBond/training"
    assert ex["matched"] == [f"{label}#{i}" for i in sorted(sizes)]
    entries = {int(e.index): e for e in _training(family, tuple(sizes))}
    for m in ex["matches"]:
        entry = entries[m["index"]]
        assert m["kinetics"] == repr(entry.data) and m["rank"] == entry.rank
        # the analogue round-trips from its adjacency lists to the listed species
        for side in ("reactants", "products"):
            got = sorted(Molecule().from_adjacency_list(a).to_smiles()
                         for a in m["reaction"][side])
            assert got == sorted(m[side])
    assert cm.exclusion_matches(ex) == ex["matches"] and not cm.exclusion_refused(ex)
    assert ex["analogues_compared"] >= 5
    # the certificate certifies it (no refused exclusion left)
    certifier = cm.Certifier(fam, repeat=CH2, growth=pe_chain)
    cert = certifier.certify_channel(family, cs, ch)
    (c,) = cert["channels"]
    assert c["depository_exclusion"]["verdict"] == "excluded:except_matched"
    assert c["depository_exclusion"]["matches"] == ex["matches"]
    assert cert["depository_exclusion_refused"] == []


def test_an_undecided_exclusion_still_refuses(fam, pe_chain):
    """Only a decided exclusion certifies: the formula candidates without a
    growth provider, or with one that cannot say where a unit inserts, stay
    ``refused:unknown`` and carry no matches to override."""
    from rmgpy.kmc import boundaries as b
    family = fam["R_Addition_MultipleBond"]
    reactants, products, _cs, _ch = _backbone_scission(family)
    alkane = b.PROVIDERS[b.LinearAlkaneProvider.name]
    for ex in (cm.depository_exclusion(family, reactants, products, repeat=CH2),
               cm.depository_exclusion(family, reactants, products, repeat=CH2, growth=alkane)):
        assert ex["verdict"] == "refused:unknown"
        assert cm.exclusion_refused(ex) and cm.exclusion_matches(ex) == []


def test_formula_only_candidates_are_excluded_by_isomorphism(fam, pe_chain):
    """The same witness against a depository of only its formula-only
    candidates (456, 636, the branched 653 and 654): the formula argument
    refuses it; no grown analogue is any of them, so it is excluded by
    isomorphism and certifies under the rule estimate like
    ``excluded:formula``."""
    family = fam["R_Addition_MultipleBond"]
    reactants, products, _cs, _ch = _backbone_scission(family)
    stub = _stub_family(family, _training(family, (456, 636, 653, 654)))
    assert cm.depository_exclusion(stub, reactants, products,
                                   repeat=CH2)["verdict"] == "refused:unknown"
    ex = cm.depository_exclusion(stub, reactants, products, repeat=CH2, growth=pe_chain)
    assert ex["verdict"] == "excluded:isomorphism", ex
    assert ex["matches"] == [] and ex["analogues_compared"] > 0
    assert not cm.exclusion_refused(ex)


def test_a_reversed_entry_is_matched_in_reverse(fam, pe_chain):
    """An entry stored products -> reactants matches the grown analogue only
    in the reverse orientation."""
    family = fam["R_Addition_MultipleBond"]
    reactants, products, _cs, _ch = _backbone_scission(family)
    stub = _stub_family(family, [_reversed(e) for e in _training(family, (2916,))])
    ex = cm.depository_exclusion(stub, reactants, products, repeat=CH2, growth=pe_chain)
    assert ex["verdict"] == "excluded:except_matched", ex
    ((index, orientation),) = {(m["index"], m["orientation"]) for m in ex["matches"]}
    assert (index, orientation) == (2916, "reverse")


@pytest.mark.parametrize("index, grown", [(2915, "[CH2]CCCCCC"), (2916, "C=CCCCCCCC")])
def test_every_declared_growth_site_is_grown(fam, pe_chain, index, grown):
    """2915 needs the radical grown (n-heptyl + 1-pentene), 2916 the alkene
    (n-propyl + 1-nonene): growth at one fixed position finds at most one."""
    family = fam["R_Addition_MultipleBond"]
    reactants, products, _cs, _ch = _backbone_scission(family)
    stub = _stub_family(family, _training(family, (index,)))
    ex = cm.depository_exclusion(stub, reactants, products, repeat=CH2, growth=pe_chain)
    assert ex["verdict"] == "excluded:except_matched", ex
    assert ex["matched"] == [f"{stub.depositories[0].label}#{index}"]
    assert any(grown in m["reactants"] for m in ex["matches"]), ex["matches"]


def test_a_provider_without_growth_sites_keeps_the_formula_refusal(fam, pe_chain):
    """No growth provider, or one that names no sites (the legacy linear
    alkane): the verdict is the formula argument's, byte-identical to the base
    (the PS rules included); a provider that cannot say refuses as unknown and
    says why."""
    from rmgpy.kmc import boundaries as b
    habs = cm.depository_exclusion(fam["H_Abstraction"], [mol("[CH3]", 0), mol("CC", 10)],
                                   [mol("C", 20), mol("C[CH2]", 30)])
    assert _digest(habs) == LEGACY_HABS_DIGEST
    family = fam["R_Addition_MultipleBond"]
    reactants, products, _cs, _ch = _backbone_scission(family)
    assert _digest(cm.depository_exclusion(family, reactants, products)) == \
        LEGACY_SCISSION_PS_DIGEST
    legacy = cm.depository_exclusion(family, reactants, products, repeat=CH2)
    assert _digest(legacy) == LEGACY_SCISSION_CH2_DIGEST
    alkane = b.PROVIDERS[b.LinearAlkaneProvider.name]
    ex = cm.depository_exclusion(family, reactants, products, repeat=CH2, growth=alkane)
    assert ex["verdict"] == "refused:unknown"
    assert ex["candidates"] == legacy["candidates"]
    assert "names no growth sites" in ex["argument"]
    assert cm.exclusion_refused(ex) and cm.exclusion_refused(legacy)


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
    on the move's template reaction, its source/entry/rank/uncertainty, the
    selection direction, and for a hand-built tree compiler._rate_rule_source
    on the reaction AFTER production's normalization (swapped when RMG
    selected the rate for the opposite direction).  Returns (rate source,
    is_forward)."""
    from rmgpy.kmc import compiler as cmp
    rxn, _r, _p = chm.template_reaction(family, cs, move)
    kinetics, source, entry, fwd = family.get_kinetics(
        rxn, template_labels=list(move.template), degeneracy=1, return_all_kinetics=False)
    rxn.kinetics = kinetics
    if not fwd:
        rxn.reactants, rxn.products = rxn.products, rxn.reactants
    out = {"kind": "RMG family estimate", "available": True, "source": str(source),
           "entry": str(entry) if entry else None, "rank": getattr(entry, "rank", None),
           "uncertainty": str(getattr(kinetics, "uncertainty", None) or "UNKNOWN"),
           "family_template_direction": "forward" if fwd else "reverse"}
    if not family.auto_generated:
        out.update(cmp._rate_rule_source(family, rxn))
    return out, bool(fwd)


def _production_record(family_name, template, reaction, orientation, rate_source, reverse_of):
    """An EventRecord shaped exactly as the compiler writes it: canonical
    explicit-H graphs, canonical heavy-atom map and the executable bond_ops
    of ``reaction`` (reactants and products sharing atom ids)."""
    from rmgpy.kmc import compiler as cmp
    ops = cmp._bond_ops(reaction)
    return EventRecord(family=family_name, template=template, arity=len(reaction.reactants),
                       participant_site_types=["s"] * len(reaction.reactants),
                       reactant_multiplicities=[1] * len(reaction.reactants),
                       orientation=orientation, rate_source=rate_source,
                       reactant_graphs=cmp._graph_adjacencies(reaction.reactants),
                       product_graphs=cmp._graph_adjacencies(reaction.products),
                       atom_map=cmp._canonical_atom_map(reaction),
                       bond_ops=ops, reverse_of=reverse_of,
                       # as production: the set_* subset of the bond_ops
                       feature_ops=[op for op in ops if op["action"].startswith("set_")])


def _record_pair(fam, family_name, reactants, pathway=None, products=None, member=0):
    """The rate-owning (RMG family estimate) and the linked (reference-thermo
    reverse) EventRecord of one channel member on ``reactants``, shaped as the
    compiler publishes them: production rate source, canonical graphs, heavy
    ``atom_map``, ``bond_ops``, reciprocal ``reverse_of`` links.  As in
    production the rate owner is the record in the direction RMG selected the
    rate for: as_generated when RMG's get_kinetics is forward, else the
    reversed record.  Returns (rate owner, partner, channel set, channel)."""
    from types import SimpleNamespace
    from rmgpy.kmc.compiler import LINK_PLACEHOLDER
    family = fam[family_name]
    cs = chm.build_channels(family, reactants)
    chs = [c for c in cs.channels if (pathway is None or c.pathway == tuple(pathway))]
    ch = chs[0]
    move = ch.representatives[member]
    prods = [rm.molecule_from_rows(p) for p in rm.split_rows(move.product_rows())]
    as_generated = SimpleNamespace(reactants=list(reactants), products=prods)
    reversed_view = SimpleNamespace(reactants=prods, products=list(reactants))
    rate_source, fwd = production_rate_source(family, cs, move)
    template = ";".join(ch.pathway)
    reverse_rs = {"kind": "reference-thermo reverse", "available": True,
                  "reference": "k_family / Kc"}
    owner_view, owner_orientation, partner_view, partner_orientation = (
        (as_generated, "as_generated", reversed_view, "reversed") if fwd else
        (reversed_view, "reversed", as_generated, "as_generated"))
    owner = _production_record(family_name, template, owner_view, owner_orientation,
                               rate_source, LINK_PLACEHOLDER)
    partner = _production_record(family_name, template, partner_view, partner_orientation,
                                 reverse_rs, owner.event_id)
    owner = replace(owner, reverse_of=partner.event_id)
    return owner, partner, cs, ch


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
    record and publishes the store.  Certification is required: a record of
    the required population (enabled, an RMG family selection) whose
    certification failed (here: the fake family has no certifiable tree) is
    never published enabled; it is refused by status, with its reason, and
    otherwise equals the record of a compile without certificates."""
    plain, _db, _p = ct._fake_compiler(["pathway_A", "pathway_B"])
    base = plain.compile()
    certified, _db, _p = ct._fake_compiler(["pathway_A", "pathway_B"])
    certified.certify_rates = True
    art = certified.compile()
    assert "rate_certificates" not in base
    assert len(art["records"]) == len(base["records"])
    store = art["rate_certificates"]
    for r in art["records"]:
        assert r["rate_certificate"]["sha256"] in store
        rc = store[r["rate_certificate"]["sha256"]]
        assert not rc["certified"]
        assert r["status"] != "enabled", r["status"]
    refused = [r for r in art["records"] if r["status"] == "refused"]
    assert refused and all("rate certification" in r["status_reason"] for r in refused)
    plain_by_template = {(r["template"], r["orientation"]): r for r in base["records"]}
    for r in art["records"]:
        twin = plain_by_template[(r["template"], r["orientation"])]
        same = {k: v for k, v in r.items() if k not in (
            "event_id", "status", "status_reason", "rate_certificate", "reverse_of",
            "canonical_index")}
        assert same == {k: v for k, v in twin.items() if k in same}


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
    forward, reverse, _cs, _ch = _record_pair(fam, "H_Abstraction", [pe, closed], head.pathway)
    rc = certifier.certify_record(reverse, rate_owner=forward)
    assert rc["template_side"] == "product_graphs" and rc["complete"], rc["gaps"]
    assert rc["certified"], rc["binding"]
    _side, mols, _o = certifier.template_side(reverse)
    reacting = {a for ch in rc["certificate"]["channels"] for m in ch["members"]
                for mp in m["mappings"] for _s, _l, a in mp["labels"]}
    ext, amap, _info = cm.extend_witness(mols, reacting)
    res = certifier.recertify(reverse, ext, amap, certificate=rc, rate_owner=forward)
    assert res["same_selection"], (res["flips"], res["refusals"])


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
    """(artifact, families): the families the compile certified with are the
    database validation re-derives every certificate against."""
    compiler, db, _p = ct._fake_compiler(["pathway_A", "pathway_B"])
    compiler.certify_rates = True
    return compiler.compile(), db.families


def test_certificate_references_are_owned_by_their_record():
    """Validation binds every stored certificate to its owning record and
    derives the reference summary from the stored certificate: swapped
    references, a summary contradicting the certificate, a corrupted body and
    a re-hashed body naming another owner each refuse the artifact."""
    from rmgpy.kmc.compiler import validate_artifact
    art, families = _certified_fake_artifact()
    validate_artifact(art, families=families)
    first, second = art["records"][0], art["records"][1]
    assert first["rate_certificate"]["sha256"] != second["rate_certificate"]["sha256"]
    assert first["rate_certificate"]["owner_event_id"] == first["event_id"]
    swapped = copy.deepcopy(art)
    r0, r1 = swapped["records"][0], swapped["records"][1]
    r0["rate_certificate"], r1["rate_certificate"] = r1["rate_certificate"], r0["rate_certificate"]
    with pytest.raises(ValueError, match="owner"):
        validate_artifact(swapped, families=families)
    lying = copy.deepcopy(art)
    ref = lying["records"][0]["rate_certificate"]
    ref["whole_species"] = not ref["whole_species"]
    with pytest.raises(ValueError, match="summary"):
        validate_artifact(lying, families=families)
    corrupt = copy.deepcopy(art)
    sha = corrupt["records"][0]["rate_certificate"]["sha256"]
    corrupt["rate_certificates"][sha]["gaps"].append("tampered after hashing")
    with pytest.raises(ValueError, match="digest"):
        validate_artifact(corrupt, families=families)
    stolen = copy.deepcopy(art)
    ref = stolen["records"][0]["rate_certificate"]
    body = stolen["rate_certificates"].pop(ref["sha256"])
    body["owner"]["event_id"] = stolen["records"][1]["event_id"]
    stolen["rate_certificates"][cm.certificate_sha256(body)] = body
    ref.update(cm.certificate_reference(body))
    with pytest.raises(ValueError, match="owner"):
        validate_artifact(stolen, families=families)


def test_certification_mode_requires_coverage_of_every_record():
    """A certified compile declares ``rate_certification`` (mode required,
    coverage every record); a removed reference, references without a
    declared mode, and an unreferenced stored certificate are refused.  An
    uncertified compile declares nothing and validates."""
    from rmgpy.kmc.compiler import validate_artifact
    art, families = _certified_fake_artifact()
    assert art["rate_certification"]["mode"] == "required"
    assert art["rate_certification"]["coverage"] == "every record"
    removed = copy.deepcopy(art)
    del removed["records"][1]["rate_certificate"]
    with pytest.raises(ValueError, match="coverage"):
        validate_artifact(removed, families=families)
    undeclared = copy.deepcopy(art)
    del undeclared["rate_certification"]
    with pytest.raises(ValueError, match="certification mode"):
        validate_artifact(undeclared, families=families)
    orphan = copy.deepcopy(art)
    body = {"orphan": True}
    orphan["rate_certificates"][cm.certificate_sha256(body)] = body
    with pytest.raises(ValueError, match="unreferenced"):
        validate_artifact(orphan, families=families)
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
    assert any("member" in m or "rewrite" in m for m in foreign["binding"]["mismatches"]), \
        foreign["binding"]
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
# rework 2: forged claims, executable rewrite, reverse-selected training,
# supplier identity, record-level recertification, fixed points, transport,
# hydrogen-only records, required population, acceptance gate
# --------------------------------------------------------------------------
def _mini_artifact(records, fam, required=True):
    """A certified record set as the compiler publishes it (certification
    gate included), as the dict validate_certification reads."""
    records, store, _replacements = cm.certify_records(records, fam, required=required)
    return {"records": [r.to_dict() for r in records], "rate_certificates": store,
            "rate_certification": cm.certification_declaration()}


def _forge(art, index, body_changes=None, record_changes=None, reown=True, rebind=False):
    """A forger's edit of record ``index`` and its stored certificate: the
    record fields are changed (event id recomputed, links kept), the body is
    edited, its owner block regenerated from the changed record (``rebind``:
    and its binding evidence), re-hashed, and the reference regenerated from
    the forged body."""
    art = copy.deepcopy(art)
    data = art["records"][index]
    ref = data.pop("rate_certificate")
    body = art["rate_certificates"].pop(ref["sha256"])
    old_id = data["event_id"]
    if record_changes:
        data.update(record_changes)
        data["event_id"] = ""
        data = EventRecord.from_dict(data).to_dict()
        for k, other in enumerate(art["records"]):
            if k == index or other.get("reverse_of") != old_id:
                continue
            other["reverse_of"] = data["event_id"]
            oref = other["rate_certificate"]
            obody = art["rate_certificates"].pop(oref["sha256"])
            if obody["owner"].get("rate_owner_event_id") == old_id:
                obody["owner"]["rate_owner_event_id"] = data["event_id"]
            other["rate_certificate"] = cm.certificate_reference(obody)
            art["rate_certificates"][cm.certificate_sha256(obody)] = obody
    for key, value in (body_changes or {}).items():
        body[key] = value
    if rebind:
        body["binding"] = dict(body["binding"], evidence=cm.binding_evidence(data))
    if reown:
        owner = cm.record_owner(data)
        owner.update({k: body["owner"].get(k) for k in ("rate_owner_event_id",
                                                        "rate_source_sha256")})
        if body["owner"].get("rate_owner_event_id") == old_id:
            owner["rate_owner_event_id"] = data["event_id"]
        body["owner"] = owner
        body["event_id"] = data["event_id"]
    data["rate_certificate"] = cm.certificate_reference(body)
    art["rate_certificates"][cm.certificate_sha256(body)] = body
    art["records"][index] = data
    return art


def test_forged_certificate_claims_are_refused(fam):
    """Rehashed wrapper claims do not validate: the validator checks the
    schema, the inner ownership and witness consistency and the binding
    coverage, and recomputes every flag from the stored derivation
    (certified=True over a refused binding, certificate None, a bound
    binding without members, a bogus schema)."""
    owner, partner, _cs, _ch = _record_pair(fam, "H_Abstraction",
                                            [mol("[CH3]", 0), mol("CC", 10)])
    art = _mini_artifact([owner, partner], fam)
    cm.validate_certification(art, families=fam)
    k = next(i for i, r in enumerate(art["records"]) if r["event_id"] == owner.event_id)
    body = art["rate_certificates"][art["records"][k]["rate_certificate"]["sha256"]]
    assert body["certified"] and body["binding"]["status"] == "bound"
    refused = dict(body["binding"], status="refused", mismatches=["source: wrong rate source"])
    with pytest.raises(ValueError, match="derived"):
        cm.validate_certification(_forge(art, k, {"binding": refused}), families=fam)
    with pytest.raises(ValueError, match="schema"):
        cm.validate_certification(_forge(art, k, {"schema": "bogus/9"}), families=fam)
    with pytest.raises(ValueError, match="derived|certificate"):
        cm.validate_certification(_forge(art, k, {"certificate": None}), families=fam)
    hollow = dict(body["binding"], members=[], mismatches=[])
    with pytest.raises(ValueError, match="binding"):
        cm.validate_certification(_forge(art, k, {"binding": hollow}), families=fam)
    lying_witness = copy.deepcopy(body["certificate"])
    lying_witness["witness"]["atoms"][0][1] = "O"
    with pytest.raises(ValueError, match="witness"):
        cm.validate_certification(_forge(art, k, {"certificate": lying_witness}), families=fam)


def test_owner_binding_covers_the_rewrite_and_degeneracy(fam):
    """A record whose rewrite (bond_ops) or degeneracy changed, its event id
    recomputed and the certificate re-owned and re-hashed, is refused: the
    binding evidence names the rewrite and degeneracy it verified."""
    owner, partner, _cs, _ch = _record_pair(fam, "H_Abstraction",
                                            [mol("[CH3]", 0), mol("CC", 10)])
    art = _mini_artifact([owner, partner], fam)
    k = next(i for i, r in enumerate(art["records"]) if r["event_id"] == owner.event_id)
    with pytest.raises(ValueError, match="rewrite|bond_ops"):
        cm.validate_certification(_forge(art, k, record_changes={"bond_ops": []}), families=fam)
    with pytest.raises(ValueError, match="degeneracy"):
        cm.validate_certification(_forge(art, k, record_changes={"degeneracy": 7.0}),
                                  families=fam)


def test_required_population_must_be_certified(fam):
    """With certification required, a firing record (enabled or
    irreversible) whose rate is an RMG family selection must carry a
    certified certificate; the certification gate refuses a failed record
    and its linked partner by status (with the reason) instead."""
    owner, partner, _cs, _ch = _record_pair(fam, "H_Abstraction",
                                            [mol("[CH3]", 0), mol("CC", 10)])
    tampered = copy.deepcopy(owner.rate_source)
    tampered["entry"] = "not the entry RMG selected"
    bad = replace(owner, rate_source=tampered, event_id="")
    bad_partner = replace(partner, reverse_of=bad.event_id, event_id="")
    bad = replace(bad, reverse_of=bad_partner.event_id)
    loose = _mini_artifact([bad, bad_partner], fam, required=False)
    with pytest.raises(ValueError, match="required"):
        cm.validate_certification(loose, families=fam)
    gated = _mini_artifact([bad, bad_partner], fam, required=True)
    cm.validate_certification(gated, families=fam)
    assert {r["status"] for r in gated["records"]} == {"refused"}
    assert all("rate certification" in r["status_reason"] for r in gated["records"])
    by_id = {r["event_id"]: r for r in gated["records"]}
    for r in gated["records"]:
        assert by_id[r["reverse_of"]]["reverse_of"] == r["event_id"]


def test_binding_replays_the_executable_rewrite(fam, certifier):
    """Binding replays the record's explicit-H bond_ops against its declared
    products under a COMPLETE canonical heavy-atom map: cleared or foreign
    operations and a one-pin map for two carbons refuse; an equivalent
    methyl H (same products under the heavy map) still binds."""
    owner, _p, _cs, _ch = _record_pair(fam, "H_Abstraction", [mol("[CH3]", 0), mol("CC", 10)])
    rc = certifier.certify_record(owner)
    assert rc["certified"], rc["binding"]
    cleared = certifier.certify_record(replace(owner, bond_ops=[], event_id=""))
    assert cleared["binding"]["status"] == "refused" and not cleared["certified"]
    assert any("rewrite" in m for m in cleared["binding"]["mismatches"]), cleared["binding"]
    # the abstracted H moved to an equivalent H of the same methyl: same
    # products under the heavy map
    ops = copy.deepcopy(owner.bond_ops)
    broken = next(op for op in ops if op["action"] == "break")
    graphs = cm.witness_molecules(owner.reactant_graphs)
    atoms = {a.id: a for m in graphs for a in m.atoms}
    carbon, hydrogen = (broken["atoms"] if atoms[broken["atoms"][0]].is_carbon()
                        else broken["atoms"][::-1])
    twin = next(n.id for n in atoms[carbon].edges if n.is_hydrogen() and n.id != hydrogen)
    for op in ops:
        if "atoms" in op:
            op["atoms"] = [twin if a == hydrogen else a for a in op["atoms"]]
        if op.get("atom") == hydrogen:
            op["atom"] = twin
    equivalent = certifier.certify_record(replace(owner, bond_ops=ops, event_id=""))
    assert equivalent["binding"]["status"] == "bound", equivalent["binding"]
    rec, _p2, _cs2, _ch2 = _record_pair(fam, "R_Recombination", [mol("[CH3]", 0), mol("[CH3]", 10)])
    assert len(rec.atom_map) == 2
    one_pin = dict(list(rec.atom_map.items())[:1])
    partial = certifier.certify_record(replace(rec, atom_map=one_pin, event_id=""))
    assert partial["binding"]["status"] == "refused" and not partial["certified"]
    assert any("heavy" in m for m in partial["binding"]["mismatches"]), partial["binding"]


def test_reverse_selected_training_rates_bind(fam, certifier):
    """CH2CH3* + CH4 on (C_methane, C_rad/H2/Cs\\H3): RMG selects training
    entry C2H6 + CH3 <=> C2H5 + CH4 for the OPPOSITE direction
    (is_forward False).  Production swaps the reaction before extracting
    its contributors (training reverse False) and the rate owner is the
    reversed record; the certifier mirrors that normalization and binds."""
    owner, partner, _cs, ch = _record_pair(fam, "H_Abstraction", [mol("[CH2]C", 0), mol("C", 10)],
                                           ("C_methane", "C_rad/H2/Cs\\H3"))
    assert owner.orientation == "reversed"
    assert owner.rate_source["family_template_direction"] == "reverse"
    assert owner.rate_source["training"][0]["reverse"] is False
    rc = certifier.certify_record(owner)
    assert rc["binding"]["status"] == "bound", rc["binding"]
    assert rc["certified"], rc["gaps"]
    linked = certifier.certify_record(partner, rate_owner=owner)
    assert linked["binding"]["status"] == "bound", linked["binding"]
    flipped = dict(owner.rate_source, family_template_direction="forward")
    wrong = certifier.certify_record(replace(owner, rate_source=flipped, event_id=""))
    assert any("direction" in m for m in wrong["binding"]["mismatches"]), wrong["binding"]


def test_autogenerated_supplier_identity_is_bound(fam, certifier):
    """An auto-generated tree rate is bound to its SUPPLYING node (RMG's
    'Estimated from node' and the walk), which production serializes as the
    rate's uncertainty correlation; a stated supplier that differs refuses
    even when source, entry and rank agree (a non-exact selection returns
    entry None, so nothing else names the supplier)."""
    owner, _p, _cs, ch = _record_pair(fam, "R_Recombination",
                                      [mol("[CH3]", 0), mol("[CH2]CC", 10)])
    rc = certifier.certify_record(owner)
    assert rc["binding"]["status"] == "bound", rc["binding"]
    (c,) = [c for c in rc["certificate"]["channels"]]
    supplier = c["selection"]["supplying_node"]
    assert rc["binding"]["supplier"] == supplier
    assert f"correlation='{supplier}'" in owner.rate_source["uncertainty"]
    lying = dict(owner.rate_source, uncertainty=owner.rate_source["uncertainty"].replace(
        f"correlation='{supplier}'", "correlation='Root'"))
    bad = certifier.certify_record(replace(owner, rate_source=lying, event_id=""))
    assert any("supplier" in m for m in bad["binding"]["mismatches"]), bad["binding"]
    # a non-exact selection: entry None, the supplier is named only by RMG's
    # comment and the record's uncertainty correlation
    rmg = {"source": "rate rules", "entry": None, "rank": None, "is_forward": True,
           "sources": None, "supplier": "Root_parent"}
    prs = {"available": True, "source": "rate rules", "entry": None, "rank": None,
           "family_template_direction": "forward",
           "uncertainty": "RateUncertainty(mu=0.0, var=1.0, Tref=1000.0, N=1, "
                          "data_mean=0.0, correlation='Root_parent',)"}
    assert cm.production_mismatches(prs, rmg, "leaf", supplier="Root_parent") == []
    other = dict(prs, uncertainty=prs["uncertainty"].replace("Root_parent", "Root_other"))
    assert any("supplier" in m for m in cm.production_mismatches(
        other, rmg, "leaf", supplier="Root_parent"))
    assert any("supplier" in m for m in cm.production_mismatches(
        prs, rmg, "leaf", supplier="Root_walked_elsewhere"))
    undirected = {k: v for k, v in prs.items() if k != "family_template_direction"}
    assert any("direction" in m for m in cm.production_mismatches(
        undirected, rmg, "leaf", supplier="Root_parent"))


def test_record_recertification_requires_a_certified_owned_wrapper(fam, certifier):
    """Record-level recertification refuses a wrapper that is not certified
    (binding refused), and a wrapper owned by another record, instead of
    unwrapping its inner certificate."""
    closed, pe = mol(ct.ps(3), 17), mol("C[CH]c1ccccc1", 0)
    family = fam["H_Abstraction"]
    cs = chm.build_channels(family, [pe, closed])
    head = next(c for c in cs.channels if c.pathway[0].startswith("C/H3/Cs")
                and c.pathway[1].startswith("C_rad/H/CbCs"))
    owner, partner, _cs, _ch = _record_pair(fam, "H_Abstraction", [pe, closed], head.pathway)
    rc = certifier.certify_record(owner)
    assert rc["certified"]
    _side, mols, _o = certifier.template_side(owner)
    ext, amap, _info = cm.extend_witness(mols, _reacting(rc["certificate"]))
    ok = certifier.recertify(owner, ext, amap, certificate=rc)
    assert ok["same_selection"] and not ok["refusals"], ok["refusals"]
    uncertified = dict(rc, certified=False, binding=dict(
        rc["binding"], status="refused", mismatches=["source: wrong rate source"]))
    res = certifier.recertify(owner, ext, amap, certificate=uncertified)
    assert res["same_selection"] is False and not res["fixed_point"]
    assert any("certified" in r for r in res["refusals"]), res["refusals"]
    other = certifier.certify_record(partner, rate_owner=owner)
    res = certifier.recertify(owner, ext, amap, certificate=other)
    assert res["same_selection"] is False and any("owner" in r for r in res["refusals"])
    res = certifier.recertify(owner, ext, amap, certificate=rc["certificate"])
    assert res["same_selection"] is False and res["refusals"]


def test_fixed_point_reports_added_mappings(fam, certifier):
    """A mapping that appears on the extension within the recertified cohort
    is reported and is not a fixed point."""
    family = fam["H_Abstraction"]
    reactants = [mol("[CH3]", 0), mol("CC", 10)]
    cs = chm.build_channels(family, reactants)
    cert = certifier.certify_channel(family, cs, only(cs.channels))
    ident = {a.id: a.id for m in cs.reactants for a in m.atoms}
    assert cm.compare_certificates(cert, cert, ident)["fixed_point"]
    grown = copy.deepcopy(cert)
    member = grown["channels"][0]["members"][0]
    extra = copy.deepcopy(member["mappings"][0])
    extra["labels"] = [[s, lab, 999 if lab == "*2" else a] for s, lab, a in extra["labels"]]
    member["mappings"].append(extra)
    out = cm.compare_certificates(cert, grown, ident)
    assert out["added_mappings"] and not out["fixed_point"], out


def test_equivalent_hydrogen_renumbering_is_a_fixed_point(fam, certifier):
    """Swapping the ids of two equivalent methyl H (an automorphism of the
    witness) keeps the selection AND the transported support: embeddings
    are compared independently of the id order."""
    family = fam["H_Abstraction"]
    cs = chm.build_channels(family, [mol("[CH3]", 0), mol("CC", 10)])
    cert = certifier.certify_channel(family, cs, only(cs.channels))
    ethane = cs.reactants[1]
    carbon = next(a for a in ethane.atoms if a.is_carbon())
    h1, h2 = sorted(n.id for n in carbon.edges if n.is_hydrogen())[:2]
    swap = {a.id: a.id for m in cs.reactants for a in m.atoms}
    swap[h1], swap[h2] = h2, h1
    res = certifier.recertify_certificate(family, cert, list(cs.reactants), swap)
    assert res["same_selection"], res
    assert res["same_support"] and res["fixed_point"], res["support_changes"][:2]


def test_transport_refuses_isotope_and_attribute_changes(fam, certifier):
    """C -> 13C (or any change of a surviving atom's isotope, charge, lone
    pairs, radical electrons or explicit-H bonds) is not a transport."""
    family = fam["H_Abstraction"]
    cs = chm.build_channels(family, [mol("[CH3]", 0), mol("CC", 10)])
    cert = certifier.certify_channel(family, cs, only(cs.channels))
    ident = {a.id: a.id for m in cs.reactants for a in m.atoms}
    rows = rm.molecule_rows(cs.reactants[1])
    carbon = next(i for i, r in sorted(rows.items()) if r[0] == "C")
    rows[carbon] = (rows[carbon][0], 13) + tuple(rows[carbon][2:])
    labelled = rm.molecule_from_rows(rows)
    assert any(rm._isotope(a) == 13 for a in labelled.atoms)
    problems = cm.transport_problems(cert, [cs.reactants[0], labelled], ident)
    assert any("isotope" in p for p in problems), problems
    res = certifier.recertify_certificate(family, cert, [cs.reactants[0], labelled], ident)
    assert res["same_selection"] is False and any("isotope" in r for r in res["refusals"])


def test_hydrogen_only_records_bind(fam, certifier):
    """H* + H* -> H2: both sides have zero heavy atoms, so the complete
    canonical heavy map is empty; the binding accepts it and validates the
    explicit-H rewrite instead."""
    owner, partner, _cs, _ch = _record_pair(fam, "R_Recombination",
                                            [mol("[H]", 0), mol("[H]", 10)])
    assert owner.atom_map == {}
    rc = certifier.certify_record(owner)
    assert rc["binding"]["status"] == "bound", rc["binding"]
    assert rc["certified"], rc["gaps"]
    broken = certifier.certify_record(replace(owner, bond_ops=[], event_id=""))
    assert broken["binding"]["status"] == "refused"


# --------------------------------------------------------------------------
# rework 3: validation re-derives (never trusts stored assertions), the
# runtime operation stream, the rate owner of a supplied certificate,
# acceptance of the in-domain supported population
# --------------------------------------------------------------------------
def _owner_index(art, record):
    return next(i for i, r in enumerate(art["records"]) if r["event_id"] == record.event_id)


def _stored(art, k):
    return art["rate_certificates"][art["records"][k]["rate_certificate"]["sha256"]]


def _retemplated(body, labels):
    """The inner certificate with every channel pathway and mapping leaf
    renamed to ``labels`` (self-consistent: completeness and flags hold)."""
    cert = copy.deepcopy(body["certificate"])
    for ch in cert["channels"]:
        ch["pathway"] = list(labels)
        for mem in ch.get("members", []):
            for mp in mem["mappings"]:
                mp["leaves"] = list(labels)
    return cert


def test_fabricated_certificates_are_never_valid(fam):
    """Validation never turns stored assertions into a verdict: a certificate
    is valid only if its derivation is re-established against the loaded RMG
    database and the record.  A self-consistent fabrication (nonexistent
    template, no operations, no node semantics, an empty descent, invented
    bound members) is refused; so are a nonexistent template alone and stored
    descents that are not chains of the family tree.  Without a database the
    verdict is unverified, which refuses publication in required mode."""
    owner, partner, _cs, _ch = _record_pair(fam, "H_Abstraction",
                                            [mol("[CH3]", 0), mol("C", 10)])
    art = _mini_artifact([owner, partner], fam)
    k = _owner_index(art, owner)
    body = _stored(art, k)
    assert body["certified"]
    # the probe's fabrication, every stored claim consistent with every other
    fake = ("Nonexistent_A", "Nonexistent_B")
    cert = _retemplated(body, fake)
    cert["nodes"] = {}
    cert["descents"] = {"d": []}
    for ch in cert["channels"]:
        for mem in ch["members"]:
            for mp in mem["mappings"]:
                mp["descent"] = "d"
    binding = dict(body["binding"], members=[{"channel": 0, "member": 0, "mappings": [0]}])
    fabricated = _forge(art, k, {"certificate": cert, "template": ";".join(fake),
                                 "binding": binding},
                        record_changes={"template": ";".join(fake), "bond_ops": [],
                                        "feature_ops": []}, rebind=True)
    with pytest.raises(ValueError, match="unverified"):
        cm.validate_certification(fabricated)
    with pytest.raises(ValueError, match="unverified"):
        cm.validate_certification(art)
    cm.validate_certification(art, families=fam)
    with pytest.raises(ValueError):
        cm.validate_certification(fabricated, families=fam)
    # a nonexistent template, everything else authentic
    renamed = _forge(art, k, {"certificate": _retemplated(body, fake),
                              "template": ";".join(fake)},
                     record_changes={"template": ";".join(fake)})
    with pytest.raises(ValueError, match="does not exist"):
        cm.validate_certification(renamed, families=fam)
    # stored descents that are not chains of the family tree
    empty = copy.deepcopy(body["certificate"])
    empty["descents"] = {"d": []}
    for ch in empty["channels"]:
        for mem in ch["members"]:
            for mp in mem["mappings"]:
                mp["descent"] = "d"
    with pytest.raises(ValueError, match="family tree"):
        cm.validate_certification(_forge(art, k, {"certificate": empty}), families=fam)
    unresolved = copy.deepcopy(body["certificate"])
    group = next(lab for lab, sem in unresolved["nodes"].items() if sem["kind"] == "group")
    unresolved["nodes"][group]["atoms"][0]["atomtypes"] = ["O2s"]
    with pytest.raises(ValueError, match="family tree"):
        cm.validate_certification(_forge(art, k, {"certificate": unresolved}), families=fam)


def test_rewrite_replay_is_the_runtime_operation_stream(fam, certifier):
    """The binding replays the complete executable operation stream the
    runtime applies (bond_ops AND feature_ops, every field: radical, charge,
    lone pairs, implicit hydrogens) with the runtime's semantics and compares
    it to the declared products BEFORE any chemistry normalisation: a lone
    pair Molecule.update() would erase, a feature operation the compiler
    replay ignored and an implicit-hydrogen count are problems."""
    owner, _p, _cs, _ch = _record_pair(fam, "H_Abstraction", [mol("[CH3]", 0), mol("CC", 10)])
    assert owner.feature_ops
    assert cm.rewrite_problems(owner) == []
    graphs = cm.witness_molecules(owner.reactant_graphs)
    atoms = [a for m in graphs for a in m.atoms]
    touched = {i for op in owner.bond_ops for i in op.get("atoms", [op.get("atom")])}
    spectator = next(i for i, a in enumerate(atoms) if a.is_hydrogen() and i not in touched)
    extra = {"action": "set_radical", "atom": spectator, "value": 1}
    featured = replace(owner, feature_ops=list(owner.feature_ops) + [extra], event_id="")
    assert cm.rewrite_problems(featured), "feature_ops are executed by the runtime"
    bad = certifier.certify_record(featured)
    assert bad["binding"]["status"] == "refused" and not bad["certified"]
    carbon = next(i for i, a in enumerate(atoms) if a.is_carbon() and i not in touched)
    implicit = {"action": "set_implicit_hydrogens", "atom": carbon, "value": 1}
    assert cm.rewrite_problems(replace(owner, feature_ops=list(owner.feature_ops) + [implicit],
                                       event_id=""))
    hh, _p2, _cs2, _ch2 = _record_pair(fam, "R_Recombination", [mol("[H]", 0), mol("[H]", 10)])
    assert cm.rewrite_problems(hh) == []
    lone = list(hh.bond_ops) + [{"action": "set_lone_pairs", "atom": 0, "value": 1}]
    paired = replace(hh, bond_ops=lone, event_id="")
    problems = cm.rewrite_problems(paired)
    assert problems, "Molecule.update() would erase the requested lone pair"
    assert not certifier.certify_record(paired)["certified"]
    # declared products carrying the same lone pair: the replay equals them,
    # but chemistry normalisation would change the requested value
    declared = [g.replace("1 H u0 p0 c0", "1 H u0 p1 c0", 1) for g in hh.product_graphs]
    assert declared != list(hh.product_graphs)
    both = replace(hh, bond_ops=lone, product_graphs=declared, event_id="")
    assert any("normalisation" in p for p in cm.rewrite_problems(both)), \
        cm.rewrite_problems(both)


def test_recertification_checks_the_rate_owner(fam, certifier):
    """A supplied record certificate undergoes the artifact's ownership and
    derivation validation against the SUPPLIED rate owner (one code path):
    after the owner's rate source changed, the relinked reverse record keeps
    its event id, but its old certificate (bound to the previous owner) is
    refused; without the rate owner it cannot be checked and is refused."""
    owner, partner, _cs, _ch = _record_pair(fam, "H_Abstraction",
                                            [mol("[CH3]", 0), mol("CC", 10)])
    rc = certifier.certify_record(partner, rate_owner=owner)
    assert rc["certified"], rc["gaps"]
    _side, mols, _o = certifier.template_side(partner)
    ident = {a.id: a.id for m in mols for a in m.atoms}
    ok = certifier.recertify(partner, mols, ident, certificate=rc, rate_owner=owner)
    assert not ok["refusals"] and ok["same_selection"], ok["refusals"]
    moved = dict(owner.rate_source, rank=(owner.rate_source.get("rank") or 0) + 1)
    new_owner = replace(owner, rate_source=moved, event_id="")
    relinked = replace(partner, reverse_of=new_owner.event_id, event_id="")
    assert relinked.event_id == partner.event_id      # links are outside the event hash
    new_owner = replace(new_owner, reverse_of=relinked.event_id)
    res = certifier.recertify(relinked, mols, ident, certificate=rc, rate_owner=new_owner)
    assert res["certificate"] is None and not res["same_selection"]
    assert any("rate-owning record" in r for r in res["refusals"]), res["refusals"]
    assert cm.record_certificate_refusals(relinked, rc, rate_owner=new_owner)
    # a rate-source change the selection does not read (RMG's comment): only
    # the rate-owner check sees that the certificate is the previous owner's
    commented = dict(owner.rate_source, comment="re-estimated elsewhere")
    other_owner = replace(owner, rate_source=commented, event_id="")
    relinked = replace(partner, reverse_of=other_owner.event_id, event_id="")
    other_owner = replace(other_owner, reverse_of=relinked.event_id)
    res = certifier.recertify(relinked, mols, ident, certificate=rc, rate_owner=other_owner)
    assert res["certificate"] is None, res["refusals"]
    assert any("rate-owning record" in r or "rate source" in r for r in res["refusals"]), \
        res["refusals"]
    unchecked = certifier.recertify(partner, mols, ident, certificate=rc)
    assert unchecked["certificate"] is None
    assert any("rate-owning record" in r for r in unchecked["refusals"]), unchecked["refusals"]


def _synthetic(support="supported", disposition="eligible", eligible=True):
    """A complete certificate with one channel and a binding bound to it."""
    cert = {"complete": True, "eligible": eligible, "depository_exclusion_refused": [],
            "channels": [{"pathway": ["a", "b"], "support": support,
                          "eligibility": {"disposition": disposition}}]}
    binding = {"status": "bound", "members": [{"channel": 0, "member": 0, "mappings": [0]}],
               "mismatches": []}
    return cert, binding


def test_certified_requires_eligibility_and_publishable_bound_channels():
    """``certified`` needs the certificate domain-eligible (an input refusal
    makes it ineligible even when the bound channel is publishable) and every
    BOUND channel publishable (supported spin/pathway status and domain
    eligible): a pathway-ambiguous bound channel blocks certification."""
    cert, binding = _synthetic()
    assert cm.derived_certified(cert, binding)
    ineligible, _b = _synthetic(eligible=False)
    assert not cm.derived_certified(ineligible, binding)
    ambiguous, _b = _synthetic(support=chm.PATHWAY_AMBIGUOUS)
    assert not cm.derived_certified(ambiguous, binding)
    refused, _b = _synthetic(disposition="refused:domain")
    assert not cm.derived_certified(refused, binding)


def test_recertification_refuses_an_incomplete_or_ineligible_extension(fam, certifier,
                                                                        monkeypatch):
    """CH3* + propane (secondary H) recertified on CH3* + propan-1-ol (a far
    methyl H replaced by OH; every certified atom transported): the extended
    certificate is outside the declared domain, so recertification is
    refused (``extension:``) instead of compared; an extended certificate
    that is incomplete is refused the same way."""
    family = fam["H_Abstraction"]
    propane = mol("CCC", 10)
    cs = chm.build_channels(family, [mol("[CH3]", 0), propane])
    ch = only(cs.find(pathway=("C/H2/Cs\\H3/Cs\\H3", "C_methyl")))
    cert = certifier.certify_channel(family, cs, ch)
    assert cert["complete"] and cert["eligible"], cert["gaps"]
    ident = {a.id: a.id for m in cs.reactants for a in m.atoms}
    same = certifier.recertify_certificate(family, cert, list(cs.reactants), ident)
    assert same["fixed_point"] and not same["refusals"], same["refusals"]
    rows = rm.molecule_rows(propane)
    rows.pop(18)
    rows[12] = rows[12][:5] + (((11, 1.0), (19, 1.0), (20, 1.0), (21, 1.0)),)
    rows[21] = ("O", None, 0, 0, 0, ((12, 1.0), (22, 1.0)))
    rows[22] = ("H", None, 0, 0, 0, ((21, 1.0),))
    propanol = rm.molecule_from_rows(rows)
    amap = {a.id: a.id for m in cs.reactants for a in m.atoms if a.id != 18}
    res = certifier.recertify_certificate(family, cert, [mol("[CH3]", 0), propanol], amap)
    assert res["certificate"] is not None and not res["certificate"]["eligible"]
    assert not res["same_selection"] and not res["fixed_point"]
    assert any(r.startswith("extension:") and "domain" in r for r in res["refusals"]), \
        res["refusals"]
    real = certifier.certify

    def incomplete(*args, **kwargs):
        out = dict(real(*args, **kwargs))
        out.update(complete=False, gaps=["forced gap of the extended certificate"])
        return out
    monkeypatch.setattr(certifier, "certify", incomplete)
    res = certifier.recertify_certificate(family, cert, list(cs.reactants), ident)
    assert not res["same_selection"] and not res["fixed_point"]
    assert any(r.startswith("extension:") and "complete" in r for r in res["refusals"]), \
        res["refusals"]


def test_acceptance_gate_fails_on_supported_selection_changes():
    """The reference acceptance gate fails on an incomplete certificate, on
    any supported channel whose recertification changed its selection, and
    on ANY non-success of the in-domain supported population (support
    ``supported`` from the channel enumeration, every species accepted by
    ``domain.check``: both computed independently of the certificate):
    a refusal or an unmapped extension there fails.  Out-of-domain or
    unsupported refusals and witnesses without a +1 analogue are accounted,
    not failed; a row without its domain token is unaccounted."""
    ok = {"incomplete": 0, "recertification": {"supported;in-domain:same": 4,
                                               "supported;out-of-domain:refused": 100,
                                               "supported;in-domain:not extended": 2,
                                               "unsupported:spin;in-domain:changed": 3,
                                               "unsupported:spin;out-of-domain:refused": 5,
                                               "pathway_ambiguous;in-domain:refused": 2}}
    assert acceptance_failures(ok) == []
    assert acceptance_failures(dict(ok, incomplete=1))
    for bad in ("supported;in-domain:refused", "supported;in-domain:unmapped",
                "pathway_ambiguous,supported;in-domain:refused"):
        refused = dict(ok, recertification=dict(ok["recertification"], **{bad: 1}))
        assert any("in-domain supported" in f for f in acceptance_failures(refused)), bad
    changed = dict(ok, recertification=dict(ok["recertification"],
                                            **{"supported;out-of-domain:changed": 1}))
    assert any("supported" in f for f in acceptance_failures(changed))
    mixed = dict(ok, recertification={"pathway_ambiguous,supported;out-of-domain:changed": 1})
    assert acceptance_failures(mixed)
    grown = dict(ok, recertification={"supported;in-domain:cohort changed": 1})
    assert acceptance_failures(grown)
    assert acceptance_failures(dict(ok, recertification={"supported;in-domain:vanished": 1}))
    assert acceptance_failures(dict(ok, recertification={"supported:refused": 1}))


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


RECERT_OUTCOMES = ("same", "changed", "cohort changed", "unmapped", "refused", "not extended")


DOMAIN_TOKENS = ("in-domain", "out-of-domain")


def acceptance_failures(table) -> list:
    """Why a reference acceptance table fails the gate (empty: passes).

    Rows are keyed ``<support>;<domain>:<outcome>``: ``support`` the channel
    support tokens (``chm`` enumeration), ``domain`` ``in-domain`` when every
    species of the channel is accepted by ``domain.check`` (both independent
    of the certificate's own verdicts).  Fails on an incomplete certificate;
    on a SUPPORTED channel whose recertification on +1 repeat unit changed
    its selection or its cohort (an added mapping); on any outcome of the
    in-domain supported population other than success (``same``): a
    refusal (``refused``) or an extension that removed a reacting atom
    (``unmapped``) fails there.  A witness without a +1 analogue (``not
    extended``: no chain end to grow) has nothing to recertify and is
    outside that population.  Out-of-domain and unsupported refusals are
    accounted, not failed; an outcome outside ``RECERT_OUTCOMES`` or a row
    without its domain token is unaccounted (fails)."""
    out = []
    if table["incomplete"]:
        out.append(f"{table['incomplete']} incomplete certificates")
    for key, count in sorted(table.get("recertification", {}).items()):
        left, _sep, outcome = key.rpartition(":")
        support, _sep2, domain = left.rpartition(";")
        supported = "supported" in support.split(",")
        if outcome not in RECERT_OUTCOMES or domain not in DOMAIN_TOKENS:
            out.append(f"unaccounted recertification outcome {key!r} ({count})")
        elif supported and domain == "in-domain" and outcome not in ("same", "not extended"):
            out.append(f"{count} in-domain supported channel(s) not recertified on +1 repeat "
                       f"unit ({key})")
        elif outcome in ("changed", "cohort changed") and supported:
            out.append(f"{count} supported channel(s) changed on +1 repeat unit ({key})")
    return out


def recert_row(certifier, family, cert, mols, species):
    """Recertification on +1 repeat unit: outcome ``same`` / ``changed`` /
    ``cohort changed`` (an added mapping) / ``unmapped`` (a reacting atom was
    removed by the extension) / ``refused`` (recertification refused, e.g. a
    domain-ineligible channel; the reason is returned) / ``not extended``,
    keyed by the channel support and the domain token of ``species`` (every
    species of the channel, judged by ``domain.check``, not by the
    certificate), with the bond distance from the extension point to the
    nearest reacting atom."""
    support = ",".join(sorted({c.get("support", "?") for c in cert["channels"]})) or "none"
    in_domain = bool(species) and all(dom.check(m).accepted for m in species)
    support = f"{support};{'in-domain' if in_domain else 'out-of-domain'}"
    ext, amap, info = cm.extend_witness(mols, _reacting(cert))
    if not any(i["extended"] for i in info):
        return {"extended": False, "outcome": f"{support}:not extended"}
    distance = min((i["distance_to_protected"] for i in info
                    if i.get("extended") and i.get("distance_to_protected") is not None),
                   default=None)
    new = certifier.recertify_certificate(family, cert, ext, amap)
    if new.get("certificate") is None:
        refusals = new.get("refusals") or []
        outcome = "unmapped" if refusals and all(r.startswith("atom map")
                                                 for r in refusals) else "refused"
    elif not new["same_selection"]:
        outcome = "changed"
    else:
        outcome = "same" if not new.get("added_mappings") else "cohort changed"
    return {"extended": True, "same": new["same_selection"], "distance": distance,
            "outcome": f"{support}:{outcome}",
            "flips": sorted({f["node"] for f in new["flips"]}),
            "refusal": (new.get("refusals") or [None])[0],
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
    refusals = Counter()
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
                        row = recert_row(certifier, family, cert, tmols, species)
                        recert[row["outcome"]] += 1
                        if row.get("refusal"):
                            refusals[f"{row['outcome']}: {row['refusal'][:120]}"] += 1
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
             "recertification_refusals": dict(refusals),
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
    assert not acceptance_failures(table), acceptance_failures(table)


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
    statuses = Counter((r["status"], store[r["rate_certificate"]["sha256"]]["certified"])
                       for r in records)
    print("DP8CERT status/certified", dict(statuses))
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
        row = recert_row(certifier, fam[record["family"]], cert, mols, mols + others)
        recert[row["outcome"]] += 1
        if row.get("refusal"):
            print("DP8CERT recert refusal", row["outcome"], row["refusal"][:160])
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
                "status_certified": [[list(k), v] for k, v in sorted(statuses.items())],
                "flips": sorted([n, c] for n, c in flips.items()),
                "incomplete": [list(map(str, x)) for x in incomplete],
                "compile_and_certify_s": round(t1 - t0), "table_s": round(t2 - t1),
                "store_bytes": len(cm.canonical_json(store))}) + "\n")
    assert not incomplete
    gate = acceptance_failures({"incomplete": len(incomplete), "recertification": recert})
    assert not gate, gate
    # certification required: no record of the required population (a firing
    # RMG family selection) is uncertified; library-derived records are not in it
    by_id = {r["event_id"]: r for r in records}
    required = [r for r in records if cm.requires_certification(r, by_id)]
    assert required and all(store[r["rate_certificate"]["sha256"]]["certified"]
                            for r in required)
    assert dom.check(ct.mol(ct.ps(8), 0)).accepted


# --------------------------------------------------------------------------
# hardening: the rate owner's direction, and guards each proven by a test
# --------------------------------------------------------------------------
def _opposite_owner_pair(fam):
    """The reversed record of [CH3] + CC given the forward record's rate source
    (training entry 215, family_template_direction forward) and a forward rate
    table, with a reciprocal linked partner: a rate owner that runs opposite to
    the kinetics RMG selected."""
    from rmgpy.kmc.compiler import LINK_PLACEHOLDER
    forward, reverse, _cs, _ch = _record_pair(fam, "H_Abstraction",
                                              [mol("[CH3]", 0), mol("CC", 10)])
    assert forward.orientation == "as_generated" and reverse.orientation == "reversed"
    assert forward.rate_source["entry"].startswith("C2H6 + CH3")
    assert forward.rate_source["family_template_direction"] == "forward"
    k_table = {"700": [122.8126400424058]}
    owner = replace(reverse, rate_source=copy.deepcopy(forward.rate_source), reverse_of=None,
                    k_table=k_table, event_id="")
    owner = replace(owner, reverse_of=LINK_PLACEHOLDER, event_id="")
    partner = replace(forward, rate_source=copy.deepcopy(reverse.rate_source),
                      reverse_of=owner.event_id, k_table=None, event_id="")
    owner = replace(owner, reverse_of=partner.event_id)
    return owner, partner


def test_rate_owner_must_run_in_the_direction_rmg_selected(fam, certifier):
    """The stored source direction was compared with RMG's selection but not
    with the orientation of the record that owns the rate: the reversed record
    of [CH3] + CC carrying the FORWARD record's rate source (training entry
    215, direction forward) was bound and certified, and a reciprocal pair
    built on it validated.  The owner is the record itself for a family
    estimate and its reverse_of partner for a thermo reverse; the owner must be
    as_generated exactly when RMG selected the rate forward."""
    owner, partner = _opposite_owner_pair(fam)
    rc = certifier.certify_record(owner)
    assert rc["binding"]["status"] == "refused" and not rc["certified"], rc["binding"]
    assert any("direction: rate owner runs opposite to selected kinetics" in m
               for m in rc["binding"]["mismatches"]), rc["binding"]
    linked = certifier.certify_record(partner, rate_owner=owner)
    assert linked["binding"]["status"] == "refused" and not linked["certified"], \
        linked["binding"]
    assert any("direction: rate owner runs opposite to selected kinetics" in m
               for m in linked["binding"]["mismatches"]), linked["binding"]
    art = _mini_artifact([owner, partner], fam)
    cm.validate_certification(art, families=fam)
    assert [r["status"] for r in art["records"]] == ["refused", "refused"], \
        [r["status"] for r in art["records"]]
    # a rate owner running the selected way still binds: both orientations
    for reactants, pathway, selected in (
            ([mol("[CH3]", 0), mol("CC", 10)], None, "as_generated"),
            ([mol("[CH2]C", 0), mol("C", 10)], ("C_methane", "C_rad/H2/Cs\\H3"), "reversed")):
        good, other, _cs, _ch = _record_pair(fam, "H_Abstraction", reactants, pathway)
        assert good.orientation == selected
        assert certifier.certify_record(good)["certified"]
        assert certifier.certify_record(other, rate_owner=good)["certified"]


def _forged_binding(art, k, certificate=None, binding=None, **flags):
    """``art`` with record ``k``'s stored certificate edited (ownership and
    digests regenerated; the edit is internally consistent)."""
    changes = dict(flags)
    if certificate is not None:
        changes["certificate"] = certificate
    if binding is not None:
        changes["binding"] = binding
    return _forge(art, k, changes)


def test_stored_certificate_must_equal_its_rederivation(fam):
    """A stored certificate whose descents and tree path are genuine and whose
    claims are consistent with each other (ownership and digests regenerated),
    but which is not what the production binding re-derives from the record
    on the loaded database, is refused as not re-derived: fewer binding
    members, another selected source rank, and the flags complete / eligible
    (certified follows) each differ from the re-derivation."""
    owner, partner, _cs, _ch = _record_pair(fam, "H_Abstraction",
                                            [mol("[CH3]", 0), mol("CC", 10)])
    art = _mini_artifact([owner, partner], fam)
    cm.validate_certification(art, families=fam)
    k = next(i for i, r in enumerate(art["records"]) if r["event_id"] == owner.event_id)
    body = art["rate_certificates"][art["records"][k]["rate_certificate"]["sha256"]]
    assert len(body["binding"]["members"]) > 1

    members = copy.deepcopy(body["binding"])
    members["members"] = members["members"][:-1]
    rank = copy.deepcopy(body["certificate"])
    rank["channels"][0]["selection"]["training_entry"]["rank"] += 1
    incomplete = copy.deepcopy(body["certificate"])
    incomplete["gaps"] = ["an extra gap the re-derivation does not have"]
    incomplete["complete"] = False
    ineligible = copy.deepcopy(body["certificate"])
    channel = ineligible["channels"][0]
    channel["eligibility"]["disposition"] = "refused"
    channel["publishable"] = cm._publishable(channel)
    ineligible["eligible"] = False
    for label, forged in (
            ("binding members", _forged_binding(art, k, binding=members)),
            ("selected source rank", _forged_binding(art, k, certificate=rank)),
            ("complete flag", _forged_binding(art, k, certificate=incomplete, complete=False,
                                              certified=False)),
            ("eligible flag", _forged_binding(art, k, certificate=ineligible, eligible=False,
                                              certified=False))):
        with pytest.raises(ValueError, match="not re-derived"):
            cm.validate_certification(forged, families=fam)
        print("refused as not re-derived:", label)


def test_rewrite_compares_the_declared_lone_pairs(fam):
    """The runtime signature of a node includes its lone pairs: a declared
    product whose lone-pair count differs from the replay (a value
    Molecule.update() leaves alone on the replayed graph) is a rewrite
    problem."""
    hh, _p, _cs, _ch = _record_pair(fam, "R_Recombination", [mol("[H]", 0), mol("[H]", 10)])
    assert cm.rewrite_problems(hh) == []
    declared = [g.replace("1 H u0 p0 c0", "1 H u0 p1 c0", 1) for g in hh.product_graphs]
    assert declared != list(hh.product_graphs)
    problems = cm.rewrite_problems(replace(hh, product_graphs=declared, event_id=""))
    assert any("do not reproduce the declared products" in p for p in problems), problems


def test_rewrite_compares_the_declared_bond_count(fam):
    """A declared product with a bond the replay does not have is a rewrite
    problem (the candidate filter compares edge counts: the exact bond check
    only walks the replay's own bonds, so an extra declared bond between
    already matched atoms would pass without it); a missing bond is one too."""
    owner, _p, _cs, _ch = _record_pair(fam, "H_Abstraction", [mol("[CH3]", 0), mol("CC", 10)])
    assert cm.rewrite_problems(owner) == []
    methane = owner.product_graphs[0]
    assert "1 H u0 p0 c0 {5,S}" in methane
    extra = methane.replace("1 H u0 p0 c0 {5,S}", "1 H u0 p0 c0 {5,S} {2,S}", 1)
    problems = cm.rewrite_problems(replace(owner, product_graphs=[extra, owner.product_graphs[1]],
                                           event_id=""))
    assert any("do not reproduce the declared products" in p for p in problems), problems
    missing = methane.replace("1 H u0 p0 c0 {5,S}", "1 H u0 p0 c0", 1).replace(
        "{1,S} {2,S} {3,S} {4,S}", "{2,S} {3,S} {4,S}", 1)
    assert missing != methane
    problems = cm.rewrite_problems(replace(owner, product_graphs=[missing, owner.product_graphs[1]],
                                           event_id=""))
    assert any("do not reproduce the declared products" in p for p in problems), problems


# intra_H_migration's R5HJ_2 as the database writes it: its file order is not
# the order RMG's matcher sorts it into (by connectivity) on first use
R5HJ_2 = """
1 *1 R!H u1 {2,[S,D,T,B]}
2 *4 R!H u0 {1,[S,D,T,B]} {3,[S,D,T,B]}
3 *6 R!H u1 {2,[S,D,T,B]} {4,[S,D,T,B]}
4 *5 R!H u0 {3,[S,D,T,B]} {5,[S,D,T,B]}
5 *2 R!H u0 {4,[S,D,T,B]} {6,S}
6 *3 H   u0 {5,S}
"""


def _one_node_family(adjacency):
    from types import SimpleNamespace
    from rmgpy.data.base import Entry
    from rmgpy.molecule.group import Group
    entry = Entry(label="N", item=Group().from_adjacency_list(adjacency))
    groups = SimpleNamespace(entries={"N": entry}, top=[entry])
    return SimpleNamespace(label="toy", groups=groups), entry


def test_node_semantics_resolve_whatever_order_the_matcher_left_the_atoms_in():
    """A node's semantics recorded before RMG's matcher sorted its atoms
    resolve against the same node after (and the other way round): the
    stored claim is the group, not the in-memory atom order."""
    family, entry = _one_node_family(R5HJ_2)
    stored = cm._entry_semantics(entry)
    entry.item.sort_vertices()
    assert cm._entry_semantics(entry) != stored        # the order did change
    assert cm.tree_problems(family, {"nodes": {"N": stored}}) == []
    family, entry = _one_node_family(R5HJ_2)
    entry.item.atoms.reverse()
    assert cm.tree_problems(family, {"nodes": {"N": stored}}) == []


def test_node_semantics_of_another_group_still_do_not_resolve():
    """A different group (one radical moved) is still a stored claim that
    does not resolve, in any atom order."""
    family, entry = _one_node_family(R5HJ_2)
    other = R5HJ_2.replace("3 *6 R!H u1", "3 *6 R!H u0").replace("4 *5 R!H u0", "4 *5 R!H u1")
    _f, moved = _one_node_family(other)
    stored = cm._entry_semantics(moved)
    for _ in range(2):
        assert any("do not resolve" in p for p in cm.tree_problems(family, {"nodes": {"N": stored}}))
        entry.item.sort_vertices()
