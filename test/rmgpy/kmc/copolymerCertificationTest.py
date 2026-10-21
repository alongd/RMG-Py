"""Certification of a linear statistical copolymer (schema v2): the growth
unit of the depository exclusion is the set of residue formulas, the formula
argument covers every non-negative integer combination of them, and the
exclusion is decided by construction at the copolymer's capped residue ends,
inserting every residue word the sequence support allows there.

A single formula certifies a copolymer wrongly: one residue's formula misses
the analogues grown by the others (a spurious ``excluded:formula``), and the
gcd of the residue formulas admits sizes no residue word reaches.  The v1
homopolymer certificates (PE and PP declared chains, the PS default) are
pinned byte-identical.

Needs the RMG database (``rmg_moves_harness.DATABASE``)."""

import hashlib
import json
from types import SimpleNamespace

import pytest

import certificatesTest as ctest
import channelsTest as ct
import declared_parity_support as sup
import rmg_moves_harness as harness
from polymer_spec_v2_support import DECKS, polymer_keyword
from rmgpy.kmc import boundaries as b
from rmgpy.kmc import certificates as cm
from rmgpy.kmc import channels as chm
from rmgpy.kmc import declared_chain as dc
from rmgpy.kmc import emit
from rmgpy.kmc import windows as win
from rmgpy.kmc.declared_sequence import SequenceChainSpec, sequence_provider
from rmgpy.kmc.run import read_deck
from rmgpy.molecule import Molecule
from rmgpy.polymer_spec import PolymerSpec

pytestmark = pytest.mark.skipif(not harness.DATABASE.is_dir(),
                                reason=f"needs the RMG database at {harness.DATABASE}")

mol = ct.mol
E = (("C", 2), ("H", 4))
P = (("C", 3), ("H", 6))
CH2 = (("C", 1), ("H", 2))
TRAINING = "R_Addition_MultipleBond/training"
# sha256 of json.dumps(value, sort_keys=True) at e8df4cde7 (before the residue
# set): /home/alon/runs/pk-i234/probe/p1_premises.py and p2_pp_pins.py
V1_PINS = {
    "pe_exclusion": "3ef25f75531c3b309934ee9ac325d194203538d57fc93c7eb5476aaed6728077",
    "pe_certificate": "e38cb5d0cd45489c742e8a4065aeeba0bfcbdcae5f6235a62184e289cdc49648",
    "ps_certificate": "aef477748c91c64a0900ce3a84cfd963ca5f1bc08614779732b9570043208af8",
    "pp_habs_certificates": "79f3a318820f356b76b4ed900b79ea2aff6d12490d9db4be918e71121f6c4c61",
    "pp_radd_certificates": "8d6825ad78cd5878c6c4d4bc80638e9a794f139131316d133396b5ba0489b3af",
}


@pytest.fixture(scope="module")
def fam():
    return ct.database().kinetics.families


@pytest.fixture(scope="module")
def radd(fam):
    return fam["R_Addition_MultipleBond"]


@pytest.fixture(scope="module")
def ep():
    return sequence_provider(read_deck(DECKS / "ep_v2_deck.py").spec)


@pytest.fixture(scope="module")
def scission(radd):
    """1-pentene + n-propyl -> oct-4-yl (the PE backbone beta-scission
    witness, template direction): a chain piece of E/P words too."""
    return ctest._backbone_scission(radd)


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def _entry(index, reactants, product):
    """A depository entry built from SMILES (independently of any provider)."""
    from rmgpy.reaction import Reaction
    from rmgpy.species import Species

    def species(smiles):
        return Species(molecule=[Molecule(smiles=smiles)])

    return SimpleNamespace(index=index, label=f"constructed-{index}", rank=None, data=None,
                           item=Reaction(reactants=[species(s) for s in reactants],
                                         products=[species(product)]))


def _matched_indices(ex):
    return sorted({m["index"] for m in ex.get("matches", ())})


# -- the growth unit is the set of residue formulas -----------------------------------


def test_the_growth_unit_of_a_copolymer_is_its_residue_formulas(ep):
    """One formula per residue, sorted; EPDM (E, P and ENB share no common
    unit) has one too, it is no longer refused."""
    assert ep.unit_formula() == (E, P)
    epdm = sequence_provider(read_deck(DECKS / "epdm_v2_deck.py").spec)
    assert epdm.unit_formula() == (E, P, (("C", 9), ("H", 12)))
    assert cm.growth_units(ep.unit_formula()) == (E, P)
    assert cm.growth_units(CH2) == (CH2,) and cm.growth_units(dict(CH2)) == (CH2,)


def test_a_witness_of_two_providers_is_still_refused(ep):
    """A copolymer is ONE provider: ``growth_unit`` reads its residue set; two
    chain providers in one witness stay refused."""
    pe = dc.DeclaredChainProvider(sup.pe_declaration())
    copolymer = ep.build("a", SequenceChainSpec(("E", "P", "E")))
    _mols, parts = b.witness_reactants([copolymer])
    assert win.growth_unit(SimpleNamespace(participants=tuple(parts))) == (E, P)
    _mols, parts = b.witness_reactants([copolymer, pe.build("b", dc.ChainSpec(3))])
    with pytest.raises(win.WindowError, match="two chain providers"):
        win.growth_unit(SimpleNamespace(participants=tuple(parts)))


def test_an_epdm_artifact_is_certified_under_its_residue_set():
    spec = read_deck(DECKS / "epdm_v2_deck.py").spec
    from rmgpy.kmc.run import chemistry_sha256, declaration_sha256
    artifact = {"provenance": {"proxy_source": {
        "kind": "declaration", "proxy_units": 1, "pairs": {},
        "chemistry_sha256": chemistry_sha256(spec), "declaration_sha256": declaration_sha256(spec)}}}
    rules = emit.certification_rules(artifact, spec)
    assert rules["repeat"] == (E, P, (("C", 9), ("H", 12)))
    assert rules["growth"] is emit.declared_provider(spec)


# -- the spurious proofs of a single formula ------------------------------------------


def test_one_residue_formula_proves_a_spurious_exclusion_the_residue_set_does_not(
        radd, ep, scission):
    """On the E/P witness 1-pentene + n-propyl -> oct-4-yl: P's formula alone
    (C3H6) excludes every depository entry by formula, yet training 2915, 2916
    and 2919 ARE its E/P analogues (n-alkyl + 1-alkene -> a linear dodecyl,
    grown by E residues).  Under the residue set the formula argument keeps
    them, and construction at the declared residue ends matches exactly those
    three (2917 and 2918 need one carbon more on one species: no residue
    word, although the gcd unit CH2 admits them)."""
    reactants, products, _cs, _ch = scission
    one = cm.depository_exclusion(radd, reactants, products, repeat=P)
    assert one["verdict"] == "excluded:formula"          # the spurious proof
    gcd = cm.depository_exclusion(radd, reactants, products, repeat=CH2)
    assert {2917, 2918} <= {c["index"] for c in gcd["candidates"]}
    ex = cm.depository_exclusion(radd, reactants, products, repeat=ep.unit_formula(), growth=ep)
    assert ex["verdict"] == cm.EXCEPT_MATCHED, ex
    assert ex["repeat_units"] == [dict(E), dict(P)] and "repeat_unit" not in ex
    assert not {2917, 2918} & {c["index"] for c in ex["candidates"]}
    assert ex["matched"] == [f"{TRAINING}#{i}" for i in (2915, 2916, 2919)]
    for m in ex["matches"]:
        assert m["species_heavy_atoms"]["products"] == [12]
        assert all(set(p["residues"]) == {"E"} for p in m["periods"] if p is not None), m
    assert not cm.exclusion_refused(ex)


def test_one_residue_formula_misses_an_analogue_grown_by_another_residue(radd, ep, scission):
    """A depository entry that is a P-grown analogue (4-methylpentyl, n-propyl
    with a P residue prepended at its capped end, + 1-pentene), built from
    SMILES: E's formula alone (C2H4) excludes it by formula; the residue set
    matches it by construction."""
    reactants, products, _cs, _ch = scission
    stub = ctest._stub_family(radd, [_entry(900001, ["C=CCCC", "[CH2]CCC(C)C"],
                                            "CC(C)CCCC[CH]CCC")])
    assert cm.depository_exclusion(stub, reactants, products,
                                   repeat=E)["verdict"] == "excluded:formula"
    ex = cm.depository_exclusion(stub, reactants, products, repeat=ep.unit_formula(), growth=ep)
    assert ex["verdict"] == cm.EXCEPT_MATCHED, ex
    assert _matched_indices(ex) == [900001]
    (match,) = ex["matches"]
    assert sorted(match["reactants"]) == sorted(Molecule(smiles=x).to_smiles()
                                                for x in ("C=CCCC", "[CH2]CCC(C)C"))
    (grown,) = [p for p in match["periods"] if p is not None]
    assert grown == {"side": "head", "end": "E", "residues": ["P"]}


def test_the_formula_argument_is_integer_feasibility_over_the_residues(radd, ep, scission):
    """n-butyl + 1-pentene -> nonan-4-yl is the witness plus one CH2: no
    combination of E and P residues, so excluded by formula under the
    residue set (the gcd unit CH2 cannot exclude it)."""
    reactants, products, _cs, _ch = scission
    stub = ctest._stub_family(radd, [_entry(900002, ["C=CCCC", "[CH2]CCC"], "CCCCC[CH]CCC")])
    assert cm.depository_exclusion(stub, reactants, products,
                                   repeat=CH2)["verdict"] == "refused:unknown"
    ex = cm.depository_exclusion(stub, reactants, products, repeat=ep.unit_formula(), growth=ep)
    assert ex["verdict"] == "excluded:formula", ex


def test_a_residue_set_keeps_the_formula_of_a_species_without_its_elements(fam, ep):
    """The copolymer form of the H atom check (``certificatesTest``): H +
    oct-4-yl -> octane, grown by E residues at the methyl ends, IS training
    160-162, 171 and 172."""
    reactants, products = ctest.h_recombination(fam["R_Recombination"])
    ex = cm.depository_exclusion(fam["R_Recombination"], reactants, products,
                                 repeat=ep.unit_formula(), growth=ep)
    assert ex["verdict"] == cm.EXCEPT_MATCHED, ex
    assert ex["matched"] == ctest.H_RECOMBINATION_MATCHES


def _alternating():
    mapping = polymer_keyword(DECKS / "ep_v2_deck.py", "polymer_spec")
    mapping["sequence_support"] = {"E": ["P"], "P": ["E"]}
    mapping["initial_ensemble"]["sequence_statistics"] = {"model": "alternating",
                                                          "cycle": ["E", "P"]}
    return sequence_provider(PolymerSpec.from_mapping(mapping))


def test_construction_inserts_only_words_the_support_allows(radd, ep, scission):
    """n-pentyl (n-propyl with an E prepended) + 1-pentene: an analogue of the
    statistical E/P copolymer, not of the alternating one, where an E residue
    is never next to another E: same formula argument, decided by
    construction both ways."""
    reactants, products, _cs, _ch = scission
    stub = ctest._stub_family(radd, [_entry(900003, ["C=CCCC", "[CH2]CCCC"], "CCCCCC[CH]CCC")])
    statistical = cm.depository_exclusion(stub, reactants, products,
                                          repeat=ep.unit_formula(), growth=ep)
    assert statistical["verdict"] == cm.EXCEPT_MATCHED, statistical
    alternating = _alternating()
    ex = cm.depository_exclusion(stub, reactants, products, repeat=alternating.unit_formula(),
                                 growth=alternating)
    # the formula argument admits it (every combination of residues); no word
    # of the alternating support realizes it, so nothing is even built
    assert ex["verdict"] == "excluded:isomorphism", ex
    assert ex["candidates"] == statistical["candidates"] and ex["analogues_compared"] == 0


def test_a_copolymer_names_every_parse_of_its_capped_residue_ends(ep):
    """A growth site is a capped terminal residue matched exactly (its atoms
    with every neighbour, the join neighbour by element and bond order), one
    site per end atom carrying every residue it can be read as and the side
    it ends: n-propyl's methyl end is an E end both ways; 1-pentene adds a P
    tail end (the CH2 next to its methyl); a radical end, a vinyl end and a
    species without a whole residue are no site."""
    found = {}
    for smiles in ("[CH2]CC", "C=CCCC", "CCC[CH]CCCC", "[CH3]", "C=C", "[H]"):
        molecule = mol(smiles, 0)
        by_id = {a.id: a for a in molecule.atoms}
        sites = ep.growth_sites(molecule)
        found[smiles] = sorted(sorted((side, residue) for side, residue, _caps in coordinate)
                               for _port, coordinate in sites)
        for port, coordinate in sites:
            atom = by_id[port.atom_id]
            assert atom.radical_electrons == 0 and port.artificial
            for _side, _residue, caps in coordinate:
                assert caps and all(by_id[c] in atom.bonds and by_id[c].is_hydrogen()
                                    for c in caps)
    e_end = [("head", "E"), ("tail", "E")]
    assert found == {"[CH2]CC": [e_end], "C=CCCC": [e_end, [("tail", "P")]],
                     "CCC[CH]CCCC": [e_end, e_end, [("tail", "P")], [("tail", "P")]],
                     "[CH3]": [], "C=C": [], "[H]": []}


def test_the_insertions_are_the_support_words_within_the_budget(ep):
    """At an E end: appended words start with a successor of E, prepended
    words end with a predecessor of E (chain order); every word within the
    heavy-atom budget, its formula the sum of its residues."""
    (port, coordinate), = ep.growth_sites(mol("[CH2]CC", 0))
    options = list(ep.growth_options(coordinate, 6))
    words = sorted((side, residue, "".join(word)) for ((side, residue, _c), word), _f, _h in options)
    assert words == sorted((side, "E", w) for side in ("head", "tail")
                           for w in ("E", "P", "EE", "EP", "PE", "PP", "EEE"))
    for ((_side, _residue, _caps), word), formula, heavy in options:
        assert heavy == sum({"E": 2, "P": 3}[r] for r in word) <= 6
        assert formula == {"C": heavy, "H": 2 * heavy}
    alternating = _alternating()
    (port, coordinate), = alternating.growth_sites(mol("[CH2]CC", 0))
    words = sorted((side, "".join(word)) for ((side, _r, _c), word), _f, _h
                   in alternating.growth_options(coordinate, 8))
    assert words == [("head", "EP"), ("head", "P"), ("head", "PEP"),
                     ("tail", "P"), ("tail", "PE"), ("tail", "PEP")]


# -- the v1 homopolymers stay byte-identical -------------------------------------------


def test_v1_homopolymer_certification_is_byte_identical(fam, radd, scission):
    """PE (declared chain, decided by construction: training 2915-2919
    matched), PP (declared chain, formula argument with C3H6) and the PS
    default (C8H8): exclusions and certificates as at e8df4cde7."""
    reactants, products, cs, ch = scission
    pe = dc.DeclaredChainProvider(sup.pe_declaration())
    pp = dc.DeclaredChainProvider(read_deck(DECKS / "pp_deck.py").spec)
    got = {"pe_exclusion": _digest(cm.depository_exclusion(radd, reactants, products, repeat=CH2,
                                                           growth=pe)),
           "pe_certificate": _digest(cm.Certifier(fam, repeat=CH2, growth=pe)
                                     .certify_channel(radd, cs, ch)),
           "ps_certificate": _digest(cm.Certifier(fam).certify_channel(radd, cs, ch))}
    for key, family, witness in (
            ("pp_habs_certificates", fam["H_Abstraction"],
             [mol("CC(C)CC(C)CCC", 0), mol("[CH3]", 40)]),
            ("pp_radd_certificates", radd, [mol("C=C(C)CC(C)C", 0), mol("C[CH]CC(C)C", 30)])):
        pcs = chm.build_channels(family, witness)
        certifier = cm.Certifier(fam, repeat=pp.unit_formula(), growth=pp)
        got[key] = _digest([certifier.certify_channel(family, pcs, c)
                            for c in sorted(pcs.channels, key=lambda c: (c.pathway, c.products))])
    assert got == V1_PINS
