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
