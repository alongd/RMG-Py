"""Tests for rmgpy/kmc/channels.py: kinetic channels of RMG records.

Fast part: whole-channel normalization on small molecules, member collapse,
filtered moves (unchanged_species, forbidden) on the PS references, the
rate-selection certificate, spin scope, bounded caches, determinism, the
compiler pathway key, matcher agreement on the DP3 positive control.
Slow part (``RMG_KMC_SLOW=1``): matcher agreement on the S2 acceptance states
DP8 and DP12, and the closed PS DP8 homolysis reference compile.
"""

import json
import math
import os
import subprocess
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import pytest

from rmgpy.kmc import channels as chm
from rmgpy.kmc import rmg_moves as rm
from rmgpy.molecule import Molecule

import rmg_moves_harness as harness

TEST_DIR = Path(__file__).resolve().parent
SLOW = os.environ.get("RMG_KMC_SLOW") == "1"
TEMPERATURES = (600.0, 700.0, 800.0)

pytestmark = pytest.mark.skipif(not harness.DATABASE.is_dir(),
                                reason=f"needs the RMG database at {harness.DATABASE}")

_DB = {}


def database():
    """The S2 families with RMG's training rules prepared (as the compiler does)
    and the thermo database (for Kc)."""
    if "db" not in _DB:
        from rmgpy.data.rmg import RMGDatabase
        from rmgpy.kmc.compiler import prepare_rate_rules
        path = str(harness.DATABASE)
        db = RMGDatabase()
        db.load_kinetics(path + "/input/kinetics", reaction_libraries=[], seed_mechanisms=None,
                         kinetics_families=harness.FAMILY_NAMES, kinetics_depositories=["training"])
        db.load_thermo(path + "/input/thermo", thermo_libraries=["primaryThermoLibrary"],
                       depository=True)
        prepare_rate_rules(db.kinetics, db.thermo, verbose=False)
        _DB["db"] = db
    return _DB["db"]


@pytest.fixture(scope="module")
def db():
    return database()


@pytest.fixture(scope="module")
def fam(db):
    return db.kinetics.families


def mol(smiles, start):
    """A Molecule from SMILES with atom ids start, start + 1, ... in RMG atom order."""
    m = Molecule().from_smiles(smiles)
    for k, atom in enumerate(m.atoms):
        atom.id = start + k
    return m


def only(channels):
    assert len(channels) == 1, channels
    return channels[0]


def rule_k(family, channel, degeneracy, temperature):
    """Rate-rule / tree estimate (no training depository) at ``degeneracy``."""
    template = family.retrieve_template(list(channel.pathway))
    kinetics, _entry = family.get_kinetics_for_template(template, degeneracy)
    try:
        return float(kinetics.get_rate_coefficient(temperature, 0.0))
    except TypeError:
        return float(kinetics.get_rate_coefficient(temperature))


def assert_whole_channel(channel, k_record, ssa):
    """Exact (Fraction) and float sums of the member weights."""
    exact = chm.channel_rate_sum(channel, k_record, ssa, exact=True)
    assert exact == Fraction(ssa) * Fraction(k_record)
    weights = chm.member_weights(channel, k_record, ssa)
    assert len(weights) == channel.size
    assert set(weights.values()) == {ssa * k_record / channel.size}
    assert math.isclose(math.fsum(weights.values()), ssa * k_record, rel_tol=1e-15)


# --------------------------------------------------------------------------
# item 2: whole-channel normalization
# --------------------------------------------------------------------------
NORMALIZATION_TABLE = []


def test_six_fold_overcount_ch3_ethane(fam, db):
    """CH3* + C2H6 -> CH4 + C2H5*: ONE RMG record (degeneracy 6, k_record =
    6 k_node) and six mapped products.  Each transition carries k_record / 6;
    normalizing by the transitions whose MAPPED product equals one fixed
    product (always 1) would give each the whole k_record and the channel
    6 k_record."""
    family = fam["H_Abstraction"]
    reactants = [mol("[CH3]", 0), mol("CC", 10)]
    cs = chm.build_channels(family, reactants)
    ch = only(cs.channels)
    assert ch.size == 6 and ch.mappings == (1,) * 6 and not ch.identical_reactants
    assert ch.ssa_multiplier == 1.0 and ch.degeneracy_convention == 6.0
    (rmg,) = rm.rmg_degeneracy(db.kinetics, reactants, "H_Abstraction")
    assert rmg.degeneracy == ch.degeneracy_convention == 6.0
    for T in TEMPERATURES:
        k_node = rule_k(family, ch, 1, T)
        k_record = rule_k(family, ch, rmg.degeneracy, T)
        assert math.isclose(k_record, 6 * k_node, rel_tol=1e-12)
        assert_whole_channel(ch, k_record, 1.0)
        for member, w in chm.member_weights(ch, k_record).items():
            assert math.isclose(w, k_node, rel_tol=1e-12)
        NORMALIZATION_TABLE.append(("CH3* + C2H6", ch.size, rmg.degeneracy, 1.0, T, k_record,
                                    chm.channel_rate_sum(ch, k_record)))


def test_a_plus_a_recombination_ch3_ch3(fam, db):
    """CH3* + CH3* -> C2H6.  Two Molecule objects of one species; RMG tries both
    slot orders and both give the same mapped product: T_ch = 1 (2 mappings).
    Pair convention: unordered identical pair, N(N-1)/2 pairs, SSA multiplier
    2.  RMG's degeneracy is 1/2 (k_record = k_node / 2), so the one transition
    of one pair carries w = 2 * k_record / 1 = k_node and N(N-1)/2 pairs give
    sum ~ k_record N^2 (RMG's rate law k_record [CH3]^2 at large N)."""
    family = fam["R_Recombination"]
    reactants = [mol("[CH3]", 0), mol("[CH3]", 10)]
    cs = chm.build_channels(family, reactants)
    ch = only(cs.channels)
    assert (ch.size, ch.mappings, ch.identical_reactants) == (1, (2,), True)
    assert ch.ssa_multiplier == 2.0 and ch.pair_convention == "unordered-identical-pair N(N-1)/2"
    (rmg,) = rm.rmg_degeneracy(db.kinetics, reactants, "R_Recombination")
    assert rmg.degeneracy == ch.degeneracy_convention == 0.5
    for T in TEMPERATURES:
        k_node = rule_k(family, ch, 1, T)
        k_record = rule_k(family, ch, rmg.degeneracy, T)
        assert math.isclose(k_record, 0.5 * k_node, rel_tol=1e-12)
        assert_whole_channel(ch, k_record, 2.0)
        assert math.isclose(chm.transition_weight(ch, k_record), k_node, rel_tol=1e-12)
        NORMALIZATION_TABLE.append(("CH3* + CH3*", ch.size, rmg.degeneracy, 2.0, T, k_record,
                                    chm.channel_rate_sum(ch, k_record)))


def test_ethyl_ethyl_disproportionation(fam, db):
    """Et* + Et* -> C2H4 + C2H6 (A + A, ordered family): (acceptor A, donor B)
    and (acceptor B, donor A) are different mapped products, 3 beta H each:
    T_ch = 6; RMG degeneracy 3 = 6 x 1/2; SSA multiplier 2."""
    family = fam["Disproportionation"]
    reactants = [mol("C[CH2]", 0), mol("C[CH2]", 10)]
    cs = chm.build_channels(family, reactants)
    ch = only(cs.channels)
    assert (ch.size, ch.mappings, ch.identical_reactants) == (6, (1,) * 6, True)
    acceptors = Counter(m.order for m in ch.representatives)
    assert acceptors == {(0, 1): 3, (1, 0): 3}
    (rmg,) = rm.rmg_degeneracy(db.kinetics, reactants, "Disproportionation")
    assert rmg.degeneracy == ch.degeneracy_convention == 3.0
    for T in TEMPERATURES:
        k_record = rule_k(family, ch, rmg.degeneracy, T)
        assert_whole_channel(ch, k_record, 2.0)
        NORMALIZATION_TABLE.append(("Et* + Et*", ch.size, rmg.degeneracy, 2.0, T, k_record,
                                    chm.channel_rate_sum(ch, k_record)))


def reverse_k(db, family, cs, ch, T):
    """A homolysis record's rate: the recombination (template direction) rate at
    its forward degeneracy over Kc from the thermo database."""
    from rmgpy.reaction import Reaction
    from rmgpy.species import Species
    forward_deg = chm.template_direction_degeneracy(family, cs, ch)
    k_f = rule_k(family, ch, forward_deg, T)
    move = ch.representatives[0]
    species = {}

    def spc(molecule):
        s = Species(molecule=[molecule.copy(deep=True)])
        s.generate_resonance_structures()
        s.thermo = db.thermo.get_thermo_data(s)
        return s
    reactant = spc(cs.reactants[0])
    products = [spc(rm.molecule_from_rows(p)) for p in rm.split_rows(move.product_rows())]
    rxn = Reaction(reactants=products, products=[reactant])
    return k_f / rxn.get_equilibrium_constant(T), forward_deg


def test_hexane_homolysis_channels(fam, db):
    """n-C6H14 bond dissociation (R_Recombination reverse).  Every bond is
    mapped in both '*' orientations (2 mappings per member).  1-hexyl + H*:
    the six terminal C-H bonds, T_ch = 6; C3-C4: 2 x n-propyl (A + A products),
    T_ch = 1.  Unimolecular: SSA multiplier 1.  RMG's degeneracy is the
    recombination (template) direction's: 1 and 1/2."""
    family = fam["R_Recombination"]
    cs = chm.build_channels(family, [mol("CCCCCC", 0)], "reverse")
    hexyl = only(cs.find(products=("[CH2]CCCCC", "[H]")))
    propyl = only(cs.find(products=("[CH2]CC", "[CH2]CC")))
    assert (hexyl.size, set(hexyl.mappings)) == (6, {2})
    assert (propyl.size, propyl.mappings) == (1, (2,))
    assert hexyl.ssa_multiplier == propyl.ssa_multiplier == 1.0
    assert all(ch.pair_convention == "ordered-unimolecular" for ch in cs.channels)
    assert sum(ch.size for ch in cs.channels) == 6 + 4 + 4 + 2 + 1 + 2
    for ch, expected_forward in ((hexyl, 1.0), (propyl, 0.5)):
        for T in TEMPERATURES:
            k_record, forward_deg = reverse_k(db, family, cs, ch, T)
            assert forward_deg == expected_forward
            assert k_record > 0 and math.isfinite(k_record)
            assert_whole_channel(ch, k_record, 1.0)
            NORMALIZATION_TABLE.append((f"n-C6H14 -> {' + '.join(ch.products)}", ch.size,
                                        forward_deg, 1.0, T, k_record,
                                        chm.channel_rate_sum(ch, k_record)))


def test_resonance_donor_transfer_through_two_forms_is_one_member(fam, db):
    """CH3* + pent-1-en-3-yl (two forms).  Abstraction of a terminal methyl H is
    found from both donor forms (2 mappings) and is one member per H: T_ch = 3
    for that channel.  Its product is a biradical, so the channel is also
    ``unsupported:spin`` (item 7)."""
    family = fam["H_Abstraction"]
    reactants = [mol("[CH3]", 0), mol("CC[CH]C=C", 10)]
    cs = chm.build_channels(family, reactants)
    remote = [ch for ch in cs.channels if all(m.atom("*1") == 10 for m in ch.representatives)]
    ch = only(remote)
    assert (ch.size, ch.mappings) == (3, (2, 2, 2))
    for member in ch.members:
        forms = {m.forms for m in cs.moves if m.status == "kept" and rm.product_key(m) == member}
        assert forms == {(0, 0), (0, 1)}
    assert ch.support == chm.SPIN_UNSUPPORTED
    for T in TEMPERATURES:
        k_record = rule_k(family, ch, ch.degeneracy_convention, T)
        assert_whole_channel(ch, k_record, 1.0)
        NORMALIZATION_TABLE.append(("CH3* + pent-1-en-3-yl (CH3 site)", ch.size,
                                    ch.degeneracy_convention, 1.0, T, k_record,
                                    chm.channel_rate_sum(ch, k_record)))


def test_print_normalization_table():
    print("NORMTABLE case | T_ch | RMG degeneracy | ssa | T | k_record | sum w | sum/(ssa k)")
    for case, t_ch, deg, ssa, T, k, total in NORMALIZATION_TABLE:
        print(f"NORMTABLE {case} | {t_ch} | {deg} | {ssa} | {T:.0f} | {k:.6e} | {total:.6e} | "
              f"{total / (ssa * k):.17g}")


# --------------------------------------------------------------------------
# item 1: how raw mappings collapse into members
# --------------------------------------------------------------------------
def test_slot_orders_collapse_for_a_symmetric_family_not_an_ordered_one(fam):
    """A + B: CH3* + Et* recombination, both slot orders, one mapped product
    (T_ch = 1, 2 mappings).  Disproportionation CH3* + Et*: only (acceptor
    CH3*, donor Et*) has a donor H next to the radical: 3 members, 1 mapping
    each; acceptor and donor are never folded."""
    rc = chm.build_channels(fam["R_Recombination"], [mol("[CH3]", 0), mol("C[CH2]", 10)])
    ch = only(rc.channels)
    assert (ch.size, ch.mappings, ch.identical_reactants) == (1, (2,), False)
    assert {m.order for m in rc.moves} == {(0, 1), (1, 0)}
    ds = chm.build_channels(fam["Disproportionation"], [mol("[CH3]", 0), mol("C[CH2]", 10)])
    ch = only(ds.channels)
    assert (ch.size, ch.mappings) == (3, (1, 1, 1))


def test_template_labels_are_part_of_the_channel_key(fam):
    """Two pathways with the same reactant and product species are two
    channels: the six CH3* + C2H6 moves, three of them relabelled to another
    template, give two channels of three members (one record each)."""
    from dataclasses import replace
    reactants = [mol("[CH3]", 0), mol("CC", 10)]
    moves = rm.generate_moves(fam["H_Abstraction"], reactants)
    relabelled = [replace(m, template=("C/H3/Cs", "C_methyl")) if m.atom("*1") == 11 else m
                  for m in moves]
    cs = chm.channels_from_moves("H_Abstraction", "forward", reactants, relabelled)
    assert sorted((ch.pathway, ch.size) for ch in cs.channels) == [
        (("C/H3/Cs", "C_methyl"), 3), (("C/H3/Cs\\H3", "C_methyl"), 3)]
    assert len({ch.products for ch in cs.channels}) == 1


def test_raw_mappings_are_diagnostics_only(fam):
    """Benzyl + ethane: 24 mappings (4 acceptor forms x 6 donor H) are 24
    different mapped products in 4 channels (CH2, ortho, ortho', para acceptor
    sites give different product species): toluene's channel has 6 members."""
    cs = chm.build_channels(fam["H_Abstraction"], [mol("[CH2]c1ccccc1", 0), mol("CC", 100)])
    assert sum(ch.size for ch in cs.channels) == 24
    assert all(set(ch.mappings) == {1} for ch in cs.channels)
    toluene = [ch for ch in cs.channels if any("C1=CC=CC=C1" in p.upper() or "c1ccccc1" in p
                                              for p in ch.products)]
    assert toluene and all(ch.size == 6 for ch in toluene)


# --------------------------------------------------------------------------
# item 4: pathway retention and padded-rate provenance in the compiler
# --------------------------------------------------------------------------
def _fake_compiler(templates, padding=None):
    """compilerTest's fake pipeline: one proxy, one family returning one
    reaction per template label, all with identical species sides."""
    from types import SimpleNamespace
    import compilerTest as ct
    from rmgpy.kmc.compiler import EventSetCompiler, SiteProxy

    class Molecule(ct._Molecule):
        # a parseable graph, so the duplicate-pair branch can compare sides
        def copy(self, deep=True):
            return Molecule()

        def to_adjacency_list(self, remove_h=False):
            return "1 H u0 p0 c0 {2,S}\n2 H u0 p0 c0 {1,S}"

    molecule = Molecule()
    reactions = []
    for label in templates:
        reaction = ct._Reaction(molecule)
        reaction.template = [SimpleNamespace(label=label)]
        reactions.append(reaction)

    class Database(ct._Database):
        def generate_reactions_from_families(self, reactants, products=None,
                                             only_families=None, resonance=True):
            self.calls.append((only_families, resonance))
            return list(reactions) if "H_Abstraction" in only_families else []

    database = Database(reactions[0])
    kwargs = {} if padding is None else {"proxy_padding_distance": padding}
    proxy = SiteProxy("interior", [ct._Species(molecule)])
    return EventSetCompiler(database, [proxy], ["H_Abstraction"],
                            temperature_grid=[600.0, 700.0], **kwargs), database, proxy


def test_compiler_keeps_two_pathways_with_identical_species_sides():
    """One family, two templates (pathways), identical reactant and product
    sides: the pair key carries the template, so both pairs are compiled and
    both survive to the artifact (the fake reaction yields one record per
    pair; without the template in the key the second pathway is dropped)."""
    compiler, _db, _proxy = _fake_compiler(["pathway_A", "pathway_B"])
    artifact = compiler.compile()
    templates = Counter(r["template"] for r in artifact["records"])
    assert templates == {"pathway_A": 1, "pathway_B": 1}
    single, _db, _proxy = _fake_compiler(["pathway_A"])
    assert len(single.compile()["records"]) == 1
    discovery = Counter(d["template"] for d in artifact["discovery"])
    assert discovery == {"pathway_A": 1, "pathway_B": 1}


def test_padded_rate_provenance_comes_from_the_padded_rate_reaction(monkeypatch):
    """When a padded witness supplies the rate, the rate-rule provenance
    (RMG's comment and the rules it names) is read from that padded rate
    reaction, not from the unpadded estimate."""
    import copy
    from types import SimpleNamespace
    import compilerTest as ct
    import rmgpy.kmc.compiler as cm

    compiler, database, proxy = _fake_compiler(["pathway_A"], padding=2)
    proxy = cm.SiteProxy("interior", proxy.reactants, artificial_boundaries=("port",),
                         repeat_unit="unit")
    compiler.proxies = [proxy]
    seen = []

    class Family:
        auto_generated = False
        label = "H_Abstraction"

        def get_kinetics(self, reaction, template_labels, degeneracy, return_all_kinetics):
            kinetics = copy.deepcopy(ct._Kinetics())
            kinetics.comment = "padded rule comment"
            return kinetics, "rate rules", None, True

        def extract_source_from_comments(self, reaction):
            seen.append(reaction.kinetics.comment)
            return False, ["H_Abstraction", {"exact": True, "rules": [], "training": []}]

    database.families = {"H_Abstraction": Family()}
    database.reaction.kinetics = copy.deepcopy(ct._Kinetics())
    database.reaction.kinetics.comment = "unpadded estimate comment"

    class Padded:
        reaction = copy.deepcopy(database.reaction)

        def provenance(self, *args):
            return {"status": "padded", "reacting_atoms": {"reactants": [], "products": []},
                    "artificial_boundaries": [],
                    "executable_to_witness_projection": {"reactants": [], "products": []}}

    monkeypatch.setattr(cm, "_map_proxy_boundaries_to_reaction", lambda p, r: p)
    monkeypatch.setattr(cm, "pad_reaction_witness", lambda *a, **k: Padded())
    monkeypatch.setattr(cm, "_production_padded_root",
                        lambda db_, reaction, provider: (copy.deepcopy(reaction), 1))
    monkeypatch.setattr(cm, "_remap_witness_projection", lambda projection, *a: projection)
    monkeypatch.setattr(cm, "_remap_padding_references", lambda metadata, *a: metadata)
    monkeypatch.setattr(cm, "matches_head_to_tail", lambda record: False)
    records, _initiating = compiler._build_linked_family_pair(proxy, database.reaction, {})
    assert seen == ["padded rule comment"]
    forward = next(r for r in records if r.rate_source.get("comment"))
    assert forward.rate_source["comment"] == "padded rule comment"
    assert forward.proxy_padding["status"] == "padded"


# --------------------------------------------------------------------------
# PS references (DP3 positive control, DP8 head radical, DP4 + 1-phenylethyl)
# --------------------------------------------------------------------------
def ps(n):
    return "".join("CC(c1ccccc1)" for _ in range(n))


def dp3_molecules():
    import rmgMovesTest as rmt
    mols = rm.molecules_from_state(rmt.dp3_state())
    radical = next(m for m in mols if len(m.atoms) == 49)
    closed = next(m for m in mols if len(m.atoms) == 50)
    h_atom = next(m for m in mols if len(m.atoms) == 1)
    return radical, closed, h_atom


def ps_references():
    radical, closed, h_atom = dp3_molecules()
    return {
        "DP3": [radical, closed, h_atom],
        "DP8 head radical": [mol("[CH2]C(c1ccccc1)" + ps(7), 0)],
        "DP4 + 1-phenylethyl": [mol(ps(4), 0), mol("C[CH]c1ccccc1", 1000)],
    }


# every direction RMG generates for each S2 family (own-reverse families: forward)
DIRECTIONS = {"H_Abstraction": ("forward",), "Disproportionation": ("forward", "reverse"),
              "R_Recombination": ("forward", "reverse"),
              "R_Addition_MultipleBond": ("forward", "reverse"),
              "intra_H_migration": ("forward",)}


def reference_channel_sets(fam, mols, families=DIRECTIONS, species=None):
    species = species or chm.SpeciesKeys()
    out = []
    for name, directions in sorted(families.items()):
        for direction in directions:
            for m in mols:
                cs = chm.build_channels(fam[name], [m], direction, species=species)
                if cs.moves:
                    out.append(cs)
            for i in range(len(mols)):
                for j in range(i + 1, len(mols)):
                    cs = chm.build_channels(fam[name], [mols[i], mols[j]], direction,
                                            species=species)
                    if cs.moves:
                        out.append(cs)
    return out


# --------------------------------------------------------------------------
# item 5: unchanged_species moves
# --------------------------------------------------------------------------
def test_unchanged_species_moves_carry_template_labels(fam):
    """rmg_moves reads the template before RMG's same-species refusal (the
    status stays ``unchanged_species``); CH3* + CH4 identity moves read
    (C_methane, C_methyl)."""
    moves = rm.generate_moves(fam["H_Abstraction"], [mol("[CH3]", 0), mol("C", 10)])
    assert {m.status for m in moves} == {"unchanged_species"}
    assert {m.template for m in moves} == {("C_methane", "C_methyl")}
    assert all(m.reason == "products are the reactants" for m in moves)


def test_dp3_chain_to_chain_transfer_gets_the_kept_weight(fam):
    """DP3 positive control: the benzylic radical of chain 0 taking a tail-CH2
    H of chain 1 swaps the two species (RMG: unchanged_species, 2 moves).  The
    synthetic channel has both moves as members (T_ch = 2), its kinetics come
    from RMG's family kinetics for (C/H2/CbCs, C_rad/H/CbCs) with provenance
    synthetic:unchanged_species, and each move's weight is finite, positive and
    equal to the weight of a kept move of the same template and node (the same
    transfer from a DP4 chain, where the species differ)."""
    import windowMatcherTest as wmt
    family = fam["H_Abstraction"]
    radical, closed, _h = dp3_molecules()
    cs = chm.build_channels(family, [radical, closed])
    (syn,) = [ch for ch in cs.channels if ch.kind == "synthetic:unchanged_species"]
    assert syn.pathway == ("C/H2/CbCs", "C_rad/H/CbCs") and syn.size == 2
    assert sorted(syn.reactants) == sorted(syn.products)
    assert {m.status for m in syn.representatives} == {"unchanged_species"}
    syn_kin = chm.channel_kinetics(family, cs, syn)
    assert syn_kin.provenance == "synthetic:unchanged_species"
    state = wmt.merge(wmt.materialize(wmt.ps_spec({3: 1})), wmt.materialize(wmt.ps_spec({4: 1})))
    wmt.remove_h(state, wmt.find(state, 0, "unit 3", "c2"))
    mols = rm.molecules_from_state(state)
    r3 = next(m for m in mols if m.multiplicity == 2 and len(m.atoms) > 1)
    c4 = next(m for m in mols if len(m.atoms) == 66)
    kept_cs = chm.build_channels(family, [r3, c4])
    kept = only(kept_cs.find(pathway=syn.pathway))
    assert kept.kind == "kept" and kept.size == 2
    kept_kin = chm.channel_kinetics(family, kept_cs, kept)
    assert kept_kin.provenance == "rmg:family"
    syn_cert = chm.selection_certificate(family, cs, syn)
    kept_cert = chm.selection_certificate(family, kept_cs, kept)
    assert syn_cert["supplying_templates"] == kept_cert["supplying_templates"] == \
        [list(syn.pathway)]
    for T in TEMPERATURES:
        w_syn = chm.member_weights(syn, syn_kin.rate(T))
        w_kept = chm.transition_weight(kept, kept_kin.rate(T))
        assert len(w_syn) == 2
        for w in w_syn.values():
            assert math.isfinite(w) and w > 0
            assert w == w_kept


# --------------------------------------------------------------------------
# item 5: forbidden moves
# --------------------------------------------------------------------------
FORBIDDEN_TABLE = []


def reference_move_sets(fam, mols):
    """(family, direction, reactants, moves) for every S2 family and direction
    over every molecule and every pair of the reference."""
    for name, directions in sorted(DIRECTIONS.items()):
        for direction in directions:
            subsets = [[m] for m in mols] + [[mols[i], mols[j]] for i in range(len(mols))
                                             for j in range(i + 1, len(mols))]
            for reactants in subsets:
                moves = rm.generate_moves(fam[name], reactants, direction)
                if moves:
                    yield name, direction, reactants, moves


def test_forbidden_disposition_table_on_the_ps_references(fam):
    """Every forbidden rule that fires on the PS references, for every S2
    family and direction (moves only; channels are built where a rule fires).
    Only H_Abstraction fires, only on DP3 (the free H atom taking an H from the
    benzylic radical chain).  Each rule is read on the labelled resonance form
    it fires on, within ``radius`` bonds of the reacting atoms; the form is a
    form of the benzyl subsystem, within ``resonance_radius`` (5) of them.
    Each excluded outcome has the sigma skeleton of a kept Disproportionation
    outcome on the same two molecules: H2 + the alkene / quinoid diene."""
    rows = []
    for label, mols in ps_references().items():
        species = chm.SpeciesKeys()
        for name, direction, reactants, moves in reference_move_sets(fam, mols):
            if not any(m.status == "forbidden" for m in moves):
                continue
            cs = chm.channels_from_moves(name, direction, reactants, moves, species=species)
            covering = {f"{other} ({d})": chm.build_channels(fam[other], reactants, d,
                                                             species=species)
                        for other in ("Disproportionation", "R_Recombination")
                        for d in DIRECTIONS[other] if (other, d) != (name, direction)
                        and len(reactants) == 2 and d == "forward"}
            for row in chm.forbidden_rule_rows(fam[name], cs, covering):
                rows.append({"reference": label, "direction": direction, **row})
    FORBIDDEN_TABLE.extend(rows)
    print("FORBIDDEN reference | family | rule | where | moves | radius | resonance radius | "
          "decidable | excluded | allowed by quotient | covered by | disposition")
    for r in rows:
        print("FORBIDDEN", r["reference"], "|", r["family"], "|", r["rule"], "|", r["where"], "|",
              r["moves"], "|", r["radius"], "|", r["resonance_radius"], "|",
              r["decidable_locally"], "|", r["excluded_outcomes"], "|",
              r["allowed_by_quotient"], "|", r["covered_by"], "|", r["disposition"])
    assert {(r["reference"], r["family"]) for r in rows} == {("DP3", "H_Abstraction")}
    assert {r["rule"] for r in rows} == {"disprop1_base_case", "disprop2", "disprop4",
                                         "disprop5"}
    for r in rows:
        assert r["decidable_locally"] and r["radius"] is not None
        assert r["disposition"] == "exclude"
        assert r["uncovered_excluded"] == 0
        assert r["covered_by"] == "Disproportionation (forward)"
    assert {r["rule"]: r["radius"] for r in rows} == {
        "disprop1_base_case": 1, "disprop2": 3, "disprop4": 5, "disprop5": 5}
    assert {r["resonance_radius"] for r in rows} == {5}


def test_forbidden_is_applied_after_quotienting_mappings(fam):
    """DP3, H* + benzylic radical chain.  Two outcomes (a ring H next to the
    radical in a quinoid form) have one kept mapping and three forbidden ones:
    they are allowed (members of an H_Abstraction channel).  Two outcomes (the
    beta-CH2 H) are forbidden in every form: excluded, disposition exclude."""
    radical, _closed, h_atom = dp3_molecules()
    cs = chm.build_channels(fam["H_Abstraction"], [radical, h_atom])
    by_outcome = defaultdict(Counter)
    for m in cs.moves:
        if m.products is not None:
            by_outcome[rm.product_key(m)][m.status] += 1
    mixed = {o for o, c in by_outcome.items() if c["forbidden"] and c["kept"]}
    pure = {o for o, c in by_outcome.items() if c["forbidden"] and not c["kept"]}
    assert len(mixed) == 2 and len(pure) == 2
    members = {o for ch in cs.channels for o in ch.members}
    assert mixed <= members
    assert {e.outcome for e in cs.excluded} == pure
    assert {e.disposition for e in cs.excluded} == {"exclude"}
    quotiented = {o: q for ch in cs.channels for o, q in zip(ch.members, ch.quotiented_forbidden)}
    assert {quotiented[o] for o in mixed} == {3}


# --------------------------------------------------------------------------
# item 7: spin scope
# --------------------------------------------------------------------------
def test_spin_scope_marks_biradical_channels(fam):
    """A channel with a reactant or product of two or more unpaired electrons
    is unsupported:spin (MeltState has no multiplicity).  On DP3: every
    H_Abstraction channel of H* + the radical chain (biradical products);
    none of the closed chain + radical chain channels.  A biradical reactant
    (but-1,4-diyl) marks its channels too."""
    radical, closed, h_atom = dp3_molecules()
    ha = chm.build_channels(fam["H_Abstraction"], [radical, h_atom])
    assert ha.channels and {ch.support for ch in ha.channels} == {chm.SPIN_UNSUPPORTED}
    pair = chm.build_channels(fam["H_Abstraction"], [radical, closed])
    assert {ch.support for ch in pair.channels} == {"supported"}
    birad = chm.build_channels(fam["H_Abstraction"], [mol("[CH3]", 0), mol("[CH2]CC[CH2]", 10)])
    assert birad.channels and {ch.support for ch in birad.channels} == {chm.SPIN_UNSUPPORTED}
    listed = chm.spin_unsupported([ha, pair, birad])
    assert len(listed) == len(ha.channels) + len(birad.channels)
    for key in listed:
        print("SPIN", key[0], key[1], ";".join(key[2]), "|", " + ".join(key[3]), "->",
              " + ".join(key[4]))


@pytest.mark.slow
@pytest.mark.skipif(not SLOW, reason="set RMG_KMC_SLOW=1")
def test_spin_unsupported_channels_on_the_ps_references(fam):
    """Listed: every unsupported:spin channel of the PS references."""
    total = 0
    for label, mols in ps_references().items():
        sets = reference_channel_sets(fam, mols)
        listed = chm.spin_unsupported(sets)
        total += len(listed)
        n = sum(len(cs.channels) for cs in sets)
        print(f"SPINREF {label}: {len(listed)} of {n} channels unsupported:spin")
        for key in listed:
            print("SPINREF", label, "|", key[0], key[1], "|", " + ".join(key[3]), "->",
                  " + ".join(key[4]))
    assert total > 0


# --------------------------------------------------------------------------
# item 6: rate-selection certificate
# --------------------------------------------------------------------------
CERTIFICATE_KINDS = Counter()


def check_descent(cert, channel):
    assert [d["leaf"] for d in cert["descent"]] == list(channel.pathway)
    for d in cert["descent"]:
        steps = d["steps"]
        for a, b in zip(steps, steps[1:]):
            assert a["choice"] == b["node"]
            chosen = [c for c in a["children"] if c["node"] == a["choice"]]
            assert len(chosen) == 1 and chosen[0]["result"] == "matched"
            assert chosen[0]["mapping"] is not None
        assert steps[-1]["choice"] is None or steps[-1]["rule"].startswith("no match")
        assert all(c["result"] == "failed" for c in steps[-1]["children"])
        for step in steps:
            for c in step["children"]:
                assert c["result"] in ("matched", "failed") and c["labelled_atoms"]
                assert c["scope"] and c["radius"] is not None


def test_certificate_rate_rule_family_h_abstraction(fam):
    """H_Abstraction (hand-built trees, rate rules) on the DP3 pair: every
    channel's certificate is supported, serializable and deterministic; the
    descent records every child predicate tried (matched / failed) mapped onto
    the witness atom ids, and its leaves are the move's template.  No exact
    training match occurs on a PS chain."""
    family = fam["H_Abstraction"]
    radical, closed, _h = dp3_molecules()
    cs = chm.build_channels(family, [radical, closed])
    kinds = Counter()
    failed = 0
    for ch in cs.channels:
        cert = chm.selection_certificate(family, cs, ch)
        assert cert["source_kind"] in ("rate_rule", "rate_rule_from_training", "averaged"), cert
        check_descent(cert, ch)
        text = chm.certificate_json(cert)
        assert chm.certificate_json(json.loads(text)) == text
        assert chm.certificate_json(chm.selection_certificate(family, cs, ch)) == text
        assert cert["contributors"] and all(c["rule"] for c in cert["contributors"])
        if cert["source_kind"] == "averaged":
            assert len(cert["supplying_templates"]) > 1 or \
                cert["contributors"][0]["rule"]["origin"] == "averaged rule"
        kinds[cert["source_kind"]] += 1
        failed += sum(c["result"] == "failed" for d in cert["descent"] for s in d["steps"]
                      for c in s["children"])
    CERTIFICATE_KINDS.update({("DP3 radical + closed", "H_Abstraction", k): v
                              for k, v in kinds.items()})
    assert "training_exact_match" not in kinds and failed > 0
    print("CERTKINDS DP3 H_Abstraction", dict(kinds))


def test_certificate_autogenerated_tree_r_recombination(fam):
    """R_Recombination (auto-generated tree) on the DP3 benzylic radical + H*:
    source kind autogenerated_tree, one merged-reactant descent, the supplying
    node equals the matched leaf.  RMG walks to a parent only through rule
    entries' parent links; in this database load no R_Recombination rule entry
    has one, so the walk is empty (recorded, not assumed)."""
    family = fam["R_Recombination"]
    radical, _closed, h_atom = dp3_molecules()
    cs = chm.build_channels(family, [radical, h_atom])
    assert len(cs.channels) == 3
    for ch in cs.channels:
        cert = chm.selection_certificate(family, cs, ch)
        assert cert["source_kind"] == "autogenerated_tree", cert
        check_descent(cert, ch)
        (d,) = cert["descent"]
        assert d["structure"] == "merged reactants"
        assert cert["matched_node"] == cert["supplying_node"] == ch.pathway[0]
        assert cert["parent_walk"] == [] and cert["returned_entry_is_none"] is False
    assert not any(e.parent for v in family.rules.entries.values() for e in v)


def test_certificate_identifies_exact_training_matches(fam):
    """CH3* + C2H6 matches H_Abstraction training reaction 215 by RMG's whole
    species isomorphism (Reaction.is_isomorphic), checked before any rule; the
    certificate says so and still records the descent."""
    family = fam["H_Abstraction"]
    cs = chm.build_channels(family, [mol("[CH3]", 0), mol("CC", 10)])
    ch = only(cs.channels)
    cert = chm.selection_certificate(family, cs, ch)
    assert cert["source_kind"] == "training_exact_match"
    assert cert["training_entry"]["index"] == 215
    assert "Reaction.is_isomorphic" in cert["training_entry"]["identified_by"]
    check_descent(cert, ch)
    print("CERTIFICATE CH3*+C2H6", chm.certificate_json(cert))


def _leaf_children(cert):
    return {c["node"]: c for d in cert["descent"] for c in d["steps"][-1]["children"]}


def test_failed_child_predicate_outside_the_leaf_changes_the_certificate(fam):
    """CH3* + propane, methyl H: leaf C/H3/Cs\\H2\\Cs, whose group maps 8 atoms
    (not the far CH3's hydrogens).  Its child C/H3/Cs\\H2\\Cs|O FAILED and
    reads those hydrogens (radius 3).  Replacing one of them by OH (propan-1-ol,
    same methyl H) is outside the leaf's atoms but inside that failed
    predicate: the predicate now matches and becomes the leaf, so the
    certificate (and RMG's rate node) changes."""
    family = fam["H_Abstraction"]
    leaf = "C/H3/Cs\\H2\\Cs"
    child = "C/H3/Cs\\H2\\Cs|O"
    base = chm.build_channels(family, [mol("[CH3]", 0), mol("CCC", 10)])
    ch = only(base.find(pathway=(leaf, "C_methyl")))
    cert = chm.selection_certificate(family, base, ch)
    donor = next(d for d in cert["descent"] if d["leaf"] == leaf)
    parent_step = donor["steps"][-2]
    leaf_map = next(c for c in parent_step["children"] if c["node"] == leaf)["mapping"]
    leaf_atoms = {a for _g, a in leaf_map}
    failed = _leaf_children(cert)[child]
    assert failed["result"] == "failed"
    outside = {row[0] for row in failed["scope"]} - leaf_atoms
    far_h = {18, 19, 20}  # the hydrogens of propane's far carbon (atom 12)
    assert far_h <= outside
    alcohol = chm.build_channels(family, [mol("[CH3]", 0), mol("CCCO", 10)])
    changed = only([c for c in alcohol.channels if c.pathway[1] == "C_methyl"
                    and all(m.atom("*1") == 10 for m in c.representatives)])
    new = chm.selection_certificate(family, alcohol, changed)
    assert changed.pathway[0] == child
    hit = next(c for d in new["descent"] for s in d["steps"] for c in s["children"]
               if c["node"] == child)
    assert hit["result"] == "matched"
    assert chm.certificate_json(new) != chm.certificate_json(cert)
    print("CERTIFICATE CH3*+propane", chm.certificate_json(cert))


def test_unsupported_selection_is_explicit(fam):
    """A selection the certificate cannot state is ``unsupported`` with a
    reason, never partial: here a channel whose template disagrees with the
    recorded descent."""
    from dataclasses import replace
    family = fam["H_Abstraction"]
    cs = chm.build_channels(family, [mol("[CH3]", 0), mol("CC", 10)])
    ch = replace(only(cs.channels), pathway=("C/H3/Cs", "C_methyl"))
    cert = chm.selection_certificate(family, cs, ch)
    assert cert["source_kind"] == "unsupported" and "descent" in cert["reason"]
    assert set(cert) == {"version", "family", "direction", "pathway", "channel_kind",
                         "reactants", "products", "source_kind", "reason"}


# --------------------------------------------------------------------------
# item 8: bounded caches
# --------------------------------------------------------------------------
def test_product_key_cache_is_bounded(monkeypatch):
    """Keying many distinct witnesses never grows the product-key cache past
    its bound (also the module default), and eviction never changes a key."""
    keys = rm.ProductKeys(max_entries=8)
    fresh = rm.ProductKeys(cache=False)
    witnesses = [rm.molecule_rows(mol("C" * n + "[CH2]", 100 * n)) for n in range(1, 31)]
    for rows in witnesses + witnesses[:5]:
        assert keys(rows) == fresh(rows)
        assert len(keys._known) <= 8
    assert rm._DEFAULT_KEYS.max_entries == rm.DEFAULT_KEY_CACHE_ENTRIES
    monkeypatch.setattr(rm._DEFAULT_KEYS, "max_entries", 4)
    rm._DEFAULT_KEYS._known.clear()
    for rows in witnesses:
        rm.product_key(rows)
        assert len(rm._DEFAULT_KEYS._known) <= 4
    with pytest.raises(rm.RmgMovesError):
        rm.ProductKeys(max_entries=0)


def test_species_keys_are_scoped_per_build(fam):
    """Channel builds keep species identity in a per-call SpeciesKeys, not in
    module state."""
    import rmgpy.kmc.channels as module
    before = {k: v for k, v in vars(module).items() if isinstance(v, (dict, list, set))}
    chm.build_channels(fam["H_Abstraction"], [mol("[CH3]", 0), mol("CC", 10)])
    after = {k: v for k, v in vars(module).items() if isinstance(v, (dict, list, set))}
    assert {k: len(v) for k, v in before.items()} == {k: len(v) for k, v in after.items()}


# --------------------------------------------------------------------------
# item 9: determinism
# --------------------------------------------------------------------------
FINGERPRINT = """
import hashlib, sys, random
sys.path.insert(0, {test_dir!r})
import channelsTest as t
from rmgpy.kmc import channels as chm
fam = t.database().kinetics.families
cases = [("H_Abstraction", ["[CH3]", "CC"], "forward"),
         ("Disproportionation", ["C[CH2]", "C[CH2]"], "forward"),
         ("R_Recombination", ["CCCCCC"], "reverse"),
         ("H_Abstraction", ["[CH3]", "CC[CH]C=C"], "forward")]
out = []
for name, smiles, direction in cases:
    out.append(t.channel_fingerprint(chm.build_channels(
        fam[name], [t.mol(s, 100 * k) for k, s in enumerate(smiles)], direction)))
radical, closed, _h = t.dp3_molecules()
cs = chm.build_channels(fam["H_Abstraction"], [radical, closed])
out.append(t.channel_fingerprint(cs))
out.append([chm.certificate_json(chm.selection_certificate(fam["H_Abstraction"], cs, ch))
            for ch in cs.channels[:5]])
print(hashlib.sha256(repr(out).encode()).hexdigest())
"""


def channel_fingerprint(cs):
    return [(ch.key, ch.kind, ch.members, ch.mappings, ch.quotiented_forbidden, ch.support)
            for ch in cs.channels] + [(e.outcome, e.disposition) for e in cs.excluded]


def permuted(molecule, seed):
    """The same molecule (same atom ids) with its atoms in another order."""
    import random
    out = molecule.copy(deep=True)
    for atom, original in zip(out.atoms, molecule.atoms):
        atom.id = original.id
    atoms = list(out.atoms)
    random.Random(seed).shuffle(atoms)
    out.atoms = atoms
    return out


def test_channels_identical_across_processes_hash_seeds_and_atom_orders(fam):
    code = FINGERPRINT.format(test_dir=str(TEST_DIR))
    digests = set()
    for seed in ("0", "4242"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                             env=env, cwd=str(TEST_DIR.parents[2]), timeout=600)
        assert out.returncode == 0, out.stderr[-2000:]
        digests.add(out.stdout.strip().splitlines()[-1])
    assert len(digests) == 1
    family = fam["H_Abstraction"]
    radical, closed, _h = dp3_molecules()
    reference = channel_fingerprint(chm.build_channels(family, [radical, closed]))
    for seed in (1, 2):
        shuffled = channel_fingerprint(chm.build_channels(
            family, [permuted(radical, seed), permuted(closed, seed + 10)]))
        assert shuffled == reference
    hexane = mol("CCCCCC", 0)
    reference = channel_fingerprint(chm.build_channels(fam["R_Recombination"], [hexane],
                                                       "reverse"))
    assert channel_fingerprint(chm.build_channels(fam["R_Recombination"], [permuted(hexane, 3)],
                                                  "reverse")) == reference


# --------------------------------------------------------------------------
# item 3: matcher agreement per channel
# --------------------------------------------------------------------------
def channel_agreement(fam, state, case, label):
    """The S2 template set's transitions on ``state``, each placed in the RMG
    channel(s) that own its mapped outcome, against the channel members from
    rmg_moves (the case's reacting-atom scope as ``select``).  Compared as
    sets of (direction, pathway, mapped outcome); one transition per member.
    An outcome that is a member of several channels (RMG reads a different
    tree leaf from another resonance form of a remote radical centre) is
    counted in each and reported; only biradical (unsupported:spin) channels
    may do that.

    Returns the table row, the (template class, pathway) table and the
    one-sided items.  A matcher outcome RMG excludes (every mapping forbidden)
    is classified by its disposition, never silently matched."""
    fam_name, direction, arity, _core, scope = harness.CASES[case]
    family = fam[fam_name]
    mols = rm.molecules_from_state(state)
    g = harness.Graph(mols)

    def select(c):
        return scope(g, c.labels) is None
    species = chm.SpeciesKeys()
    t0 = time.perf_counter()
    sets = []
    if arity == 1:
        for m in mols:
            if m.multiplicity > 1:
                sets.append(chm.build_channels(family, [m], direction, select=select,
                                               species=species))
    else:
        for i in range(len(mols)):
            for j in range(i + 1, len(mols)):
                if mols[i].multiplicity == 1 and mols[j].multiplicity == 1:
                    continue
                sets.append(chm.build_channels(family, [mols[i], mols[j]], direction,
                                               select=select, species=species))
    t1 = time.perf_counter()
    owner = defaultdict(list)
    for cs in sets:
        for ch in cs.channels:
            for o in ch.members:
                owner[o].append(ch)
    excluded = {e.outcome: e for cs in sets for e in cs.excluded}
    rmg_side = {(direction, ch.pathway, o) for o, chs in owner.items() for ch in chs}
    mat = harness.matcher_side(state, case)
    t2 = time.perf_counter()
    matcher_side, class_table, one_sided = Counter(), Counter(), []
    excluded_hits = Counter()
    multi = []
    for o, transitions in mat["by_key"].items():
        chs = owner.get(o, [])
        if chs:
            # an outcome in several channels (pathways) carries a share of each
            if len(chs) > 1:
                multi.append((o, chs))
            for ch in chs:
                for tr in transitions:
                    matcher_side[(direction, ch.pathway, o)] += 1
                    class_table[(tr.template, ch.pathway)] += 1
        elif o in excluded:
            e = excluded[o]
            excluded_hits[(e.disposition, e.rules)] += len(transitions)
        else:
            one_sided.append(("matcher-only", o, [tr.template for tr in transitions]))
    for key in sorted(rmg_side - set(matcher_side)):
        one_sided.append(("rmg-only", key[2], list(key[1])))
    duplicates = sum(1 for v in matcher_side.values() if v != 1)
    row = {"state": label, "case": case, "rmg_family": f"{fam_name} ({direction})",
           "channels": sum(len(cs.channels) for cs in sets),
           "rmg_members": len(rmg_side), "matcher_transitions": len(mat["transitions"]),
           "matched": len(rmg_side & set(matcher_side)), "duplicates": duplicates,
           "matcher_only": sum(1 for x in one_sided if x[0] == "matcher-only"),
           "multi_pathway_outcomes": len(multi),
           "multi_pathway_outcomes_spin_supported": sum(
               1 for _o, chs in multi if any(ch.support != chm.SPIN_UNSUPPORTED for ch in chs)),
           "rmg_only": sum(1 for x in one_sided if x[0] == "rmg-only"),
           "matcher_excluded_by_disposition": dict(
               (f"{d}: {'+'.join(sorted({r for _w, r in rules}))}", n)
               for (d, rules), n in sorted(excluded_hits.items())),
           "synthetic_members": sum(ch.size for cs in sets for ch in cs.channels
                                    if ch.kind != "kept"),
           "spin_unsupported_channels": sum(1 for cs in sets for ch in cs.channels
                                            if ch.support == chm.SPIN_UNSUPPORTED),
           "classes": len({c for c, _p in class_table}),
           "classes_spanning_several_pathways": sum(
               1 for c in {c for c, _p in class_table}
               if len({p for cc, p in class_table if cc == c}) > 1),
           "pathways_fed_by_several_classes": sum(
               1 for p in {p for _c, p in class_table}
               if len({c for c, pp in class_table if pp == p}) > 1),
           "rmg_s": round(t1 - t0, 1), "matcher_s": round(t2 - t1, 1)}
    return row, class_table, one_sided


def assert_agreement(row, one_sided):
    assert row["duplicates"] == 0, one_sided[:5]
    assert row["multi_pathway_outcomes_spin_supported"] == 0
    assert row["rmg_only"] == 0 and row["matcher_only"] == 0, one_sided[:5]
    assert set(k.split(":")[0] for k in row["matcher_excluded_by_disposition"]) <= {"exclude"}
    assert row["matched"] == row["rmg_members"]


def print_agreement(row, class_table):
    print("AGREEROW", json.dumps(row))
    for (cls, pathway), n in sorted(class_table.items()):
        print("AGREECLASS", row["state"], row["case"], cls, "->", ";".join(pathway), n)


def test_matcher_agreement_dp3_positive_control(fam):
    """DP3 positive control (benzylic radical chain + closed chain + H*), the
    H-abstraction S2 set: every matcher transition lies in exactly one RMG
    channel, once, and every in-scope member is hit, including the two
    chain-to-chain transfers RMG files as unchanged_species (members of the
    synthetic channel).  The other S2 sets are printed only: DP3 is not an S2
    acceptance state (its 1,5 H shift from the chain-end CH3 has no window
    among the S2 witnesses DP6/7/9)."""
    import rmgMovesTest as rmt
    state = rmt.dp3_state()
    for case in harness.CONTRACT_CASES:
        row, table, one_sided = channel_agreement(fam, state, case, "DP3")
        print_agreement(row, table)
        for x in one_sided:
            print("AGREEMIS DP3", case, x[0], x[2])
        if case == "h_abstraction_benzylic":
            assert_agreement(row, one_sided)
            assert (row["matched"], row["synthetic_members"]) == (11, 2)


# (direction, pathway, outcome) members matched.  Chain-end recombination: 5
# outcomes, 2 of them members of two pathways each (the biradical chain's
# remote benzylic centre in an aromatic or Kekule/quinoid form changes the
# R_Recombination tree leaf), all unsupported:spin: 7 members.
AGREEMENT_EXPECTED = {
    8: {"h_abstraction_benzylic": 260, "beta_scission_midchain": 6,
        "recombination_chain_end": 7, "disproportionation": 43, "h_shift": 5,
        "h_atom_recombination": 10},
    12: {"h_abstraction_benzylic": 410, "beta_scission_midchain": 6,
         "recombination_chain_end": 7, "disproportionation": 43, "h_shift": 6,
         "h_atom_recombination": 10},
}
MULTI_PATHWAY_EXPECTED = {"recombination_chain_end": 2}


@pytest.mark.slow
@pytest.mark.skipif(not SLOW, reason="set RMG_KMC_SLOW=1")
@pytest.mark.parametrize("case", harness.CONTRACT_CASES)
@pytest.mark.parametrize("dp", [8, 12])
def test_matcher_agreement_s2_states(fam, case, dp):
    """S2 acceptance states DP8 and DP12, the six S2 template sets."""
    import windowMatcherTest as wmt
    state = wmt.state_at(dp)
    row, table, one_sided = channel_agreement(fam, state, case, f"DP{dp}")
    print_agreement(row, table)
    for x in one_sided:
        print("AGREEMIS", row["state"], case, x[0], x[2])
    assert_agreement(row, one_sided)
    assert row["matched"] == AGREEMENT_EXPECTED[dp][case]
    assert row["multi_pathway_outcomes"] == MULTI_PATHWAY_EXPECTED.get(case, 0)


# --------------------------------------------------------------------------
# acceptance: the closed PS DP8 homolysis reference compile
# --------------------------------------------------------------------------
def compile_dp8_closed(out_dir):
    """The I-127 closed PS DP8 reference compile (R_Recombination; ~50 s)."""
    from rmgpy.kmc import compile as compile_module
    from rmgpy.kmc.compile import CompileRequest, compile_artifact
    from rmgpy.kmc.compiler import SiteProxy
    from rmgpy.species import Species

    def proxy_set(*_a, **_k):
        m = Molecule(smiles=ps(8))
        m.update()
        m.assign_atom_ids()
        return (SiteProxy("ref_dp8_closed", (Species(molecule=[m]),),
                          participant_site_types=("ref0",),
                          metadata={"centred": True, "context_class": "reference",
                                    "proxy_units": 0, "family_candidates": ["R_Recombination"]}),)
    original = compile_module.ps_proxy_set
    compile_module.ps_proxy_set = proxy_set
    try:
        path, _ = compile_artifact(CompileRequest(
            database_path=str(harness.DATABASE), output_dir=str(out_dir),
            families=("R_Recombination",), family_universe=("R_Recombination",),
            database_sha="i134-reference"))
    finally:
        compile_module.ps_proxy_set = original
    return json.loads(Path(path).read_bytes())


def _witness(graphs, offset=0):
    mols = []
    for graph in graphs:
        m = Molecule().from_adjacency_list(graph)
        for atom in m.atoms:
            atom.id = offset
            offset += 1
        mols.append(m)
    return mols


def dp8_reference_table(fam, artifact):
    """Per record: the channel on the record's own witness, T_ch, and the sum of
    member weights against ssa_multiplier * k_record at 600/700/800 K."""
    from rmgpy.kmc.met import RateTable
    rows = []
    cache = {}
    for record in artifact["records"]:
        family = fam[record["family"]]
        reactant_graphs = record.get("rate_witness_reactant_graphs") or record["reactant_graphs"]
        product_graphs = record.get("rate_witness_product_graphs") or record["product_graphs"]
        species = chm.SpeciesKeys()
        reactants = _witness(reactant_graphs)
        want_r = tuple(sorted(species(m) for m in reactants))
        want_p = tuple(sorted(species(m) for m in _witness(product_graphs)))
        pathway = tuple(record["template"].split(";")) if record["template"] else ()
        found, direction, cs = [], None, None
        for d in ("forward", "reverse"):
            key = (tuple(reactant_graphs), d)
            if key not in cache:
                try:
                    cache[key] = chm.build_channels(family, reactants, d, species=species)
                except rm.RmgMovesError:
                    cache[key] = None
            cs = cache[key]
            if cs is None:
                continue
            found = [ch for ch in cs.channels if ch.pathway == pathway
                     and ch.reactants == want_r and ch.products == want_p]
            if found:
                direction = d
                break
        row = {"event_id": record["event_id"][:16], "site_type": record["site_type"],
               "arity": record["arity"], "orientation": record["orientation"],
               "degeneracy": record["degeneracy"], "ssa": record["ssa_multiplier"],
               "pair_convention": record["reactant_pair_convention"],
               "template": record["template"][-50:], "status": record["status"]}
        if len(found) != 1:
            row.update({"T_ch": None, "reason": f"{len(found)} channels on the record's witness "
                                                f"with its template and species"})
            rows.append(row)
            continue
        ch = found[0]
        table = RateTable.from_mapping(record["k_table"]) if record["k_table"] else None
        sums = []
        for T in TEMPERATURES:
            k = table(T) if table else 0.0
            total = chm.channel_rate_sum(ch, k, record["ssa_multiplier"])
            target = record["ssa_multiplier"] * k
            sums.append(abs(total - target) / target if target else abs(total))
        template_deg = chm.template_direction_degeneracy(family, cs, ch, species)
        if record["degeneracy"] == ch.size:
            reason = ""
        elif direction == "reverse" and record["degeneracy"] == template_deg:
            reason = (f"degeneracy is RMG's in the template (recombination) direction "
                      f"({template_deg}); T_ch counts the {ch.size} bonds cut")
        elif ch.identical_reactants and record["degeneracy"] == ch.size / 2:
            reason = "A + A: RMG's degeneracy holds the 1/2 (T_ch / 2)"
        else:
            reason = "UNEXPLAINED"
        row.update({"direction": direction, "T_ch": ch.size, "mappings": sum(ch.mappings),
                    "channel_ssa": ch.ssa_multiplier, "channel_pair": ch.pair_convention,
                    "template_direction_degeneracy": template_deg,
                    "max_rel_error": max(sums), "reason": reason, "support": ch.support})
        rows.append(row)
    return rows


@pytest.mark.slow
@pytest.mark.skipif(not SLOW, reason="set RMG_KMC_SLOW=1")
def test_dp8_closed_reference_compile_channels(fam, tmp_path):
    """Every record of the closed PS DP8 homolysis reference compile: T_ch from
    rmg_moves on the record's own witness and sum(w) = ssa * k_record at
    600/700/800 K to 1e-12; records whose degeneracy differs from T_ch are
    listed with the reason."""
    t0 = time.perf_counter()
    artifact = compile_dp8_closed(tmp_path)
    t1 = time.perf_counter()
    rows = dp8_reference_table(fam, artifact)
    t2 = time.perf_counter()
    print(f"DP8REF compile {t1 - t0:.1f} s, table {t2 - t1:.1f} s, records {len(rows)}")
    for r in rows:
        print("DP8REF", json.dumps(r))
    located = [r for r in rows if r["T_ch"] is not None]
    for r in located:
        assert r["max_rel_error"] <= 1e-12, r
        assert r["ssa"] == r["channel_ssa"] and r["pair_convention"] == r["channel_pair"], r
        assert r["reason"] != "UNEXPLAINED", r
    missing = [r for r in rows if r["T_ch"] is None]
    print("DP8REF located", len(located), "missing", len(missing),
          Counter(r["site_type"] for r in missing))
    assert all(r["site_type"] != "ref_dp8_closed" for r in missing)


def test_benzylic_recombination_from_two_kekule_forms_is_one_member(fam):
    """2-phenylethyl* + the DP7 benzylic tail radical (the witness of a DP8
    reference record, as compiled): recombination at the benzylic carbon is found from the
    aromatic form and from both Kekule forms of its ring.  RMG returns a
    product written with a Kekule ring as a REACTIVE form beside the aromatic
    one, so the least reactive signature of each product's own forms depends
    on the Kekule pattern; product_key re-roots at the most aromatic form, and
    the outcome is one member (T_ch = 1, not 2)."""
    from rmgpy.molecule.resonance import generate_resonance_structures
    data = json.loads((TEST_DIR / "fixtures" /
                       "i134_benzylic_recombination_witness.json").read_text())
    reactants = _witness(data["reactant_graphs"])
    cs = chm.build_channels(fam["R_Recombination"], reactants)
    (ch,) = cs.find(pathway=(data["template"],))
    assert ch.size == 1 and ch.mappings[0] > 1
    bond = {a for _s, _l, a in ch.representatives[0].labels}
    moves = [m for m in cs.moves if m.status == "kept"
             and {a for _s, _l, a in m.labels} == bond]
    naive = set()
    for m in moves:
        forms = generate_resonance_structures(rm.molecule_from_rows(m.product_rows()),
                                              keep_isomorphic=True, filter_structures=True)
        naive.add(min(rm._signature(rm.molecule_rows(f)) for f in forms if f.reactive))
    assert len(naive) > 1
    assert len({rm.ProductKeys(cache=False)(m) for m in moves}) == 1
