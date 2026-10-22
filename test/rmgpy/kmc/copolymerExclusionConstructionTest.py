"""The construction behind the depository exclusion (``certificates
.decide_by_construction``) on a linear copolymer (schema v2), and on a v1
homopolymer where the same construction runs:

- one atom-mapped witness is excluded alike in both orientations: its growth
  sites are read on both sides of the reaction;
- a residue whose entry and exit are one atom grows at both ends of a
  singleton at once;
- a capped end is read as a residue only when its species parses from that
  end under the declared sequence support (never an element-only join);
- the residue words are built only for the formulas of a candidate entry,
  within that entry's species sizes.

Needs the RMG database (``rmg_moves_harness.DATABASE``)."""

import subprocess
import sys
import textwrap
import time
from pathlib import Path

import pytest

import certificatesTest as ctest
import channelsTest as ct
import copolymerCertificationTest as cct
import rmg_moves_harness as harness
from polymer_spec_v2_support import DECKS, polymer_keyword
from rmgpy.kmc import certificates as cm
from rmgpy.kmc import channels as chm
from rmgpy.kmc import declared_chain as dc
from rmgpy.kmc import rmg_moves as rm
from rmgpy.kmc.declared_sequence import sequence_provider
from rmgpy.kmc.run import read_deck
from rmgpy.polymer_spec import PolymerSpec

pytestmark = pytest.mark.skipif(not harness.DATABASE.is_dir(),
                                reason=f"needs the RMG database at {harness.DATABASE}")

mol = ct.mol


@pytest.fixture(scope="module")
def fam():
    return ct.database().kinetics.families


@pytest.fixture(scope="module")
def rrec(fam):
    return fam["R_Recombination"]


@pytest.fixture(scope="module")
def ep():
    return sequence_provider(read_deck(DECKS / "ep_v2_deck.py").spec)


def _recombination(family, radical):
    """[H] + ``radical`` -> its alkane, atom-mapped (shared ids)."""
    cs = chm.build_channels(family, [mol("[H]", 0), mol(radical, 10)])
    (ch,) = cs.channels
    return list(cs.reactants), [rm.molecule_from_rows(p) for p in
                                rm.split_rows(ch.representatives[0].product_rows())]


def _orientations(family, reactants, products, provider):
    """The exclusion of the witness and of its reverse, each with fresh caches."""
    return {name: cm.depository_exclusion(family, a, b, repeat=provider.unit_formula(),
                                          growth=provider)
            for name, a, b in (("forward", reactants, products), ("reverse", products, reactants))}


def _outcome(ex):
    return (ex["verdict"], sorted({(m["depository"], m["index"]) for m in ex.get("matches", ())}),
            sorted(site["atom"] for site in ex.get("growth_sites", ())))


def test_one_mapped_witness_is_excluded_alike_in_both_orientations(rrec, ep):
    """[H] + ethyl <=> ethane: ethyl's methyl end is no whole E (its other
    atom is the radical), ethane's is, and the reaction leaves that end
    alone.  The end is a growth site of the reaction read on either side, so
    both orientations build n-butyl + [H] <=> n-butane (a constructed entry)
    and match it; on the real depository [H] + isopropyl <=> propane does
    the same.  A v1 homopolymer too: [H] + 2-methylbutyl <=> 2-methylbutane
    under the PP declaration names the same growth sites both ways."""
    reactants, products = _recombination(rrec, "[CH2]C")
    stub = ctest._stub_family(rrec, [cct._entry(900101, ["[H]", "[CH2]CCC"], "CCCC")])
    both = _orientations(stub, reactants, products, ep)
    assert _outcome(both["forward"]) == _outcome(both["reverse"]), both
    assert both["forward"]["verdict"] == cm.EXCEPT_MATCHED, both["forward"]
    reactants, products = _recombination(rrec, "C[CH]C")
    both = _orientations(rrec, reactants, products, ep)
    assert _outcome(both["forward"]) == _outcome(both["reverse"]), both
    assert both["forward"]["verdict"] == cm.EXCEPT_MATCHED
    pp = dc.DeclaredChainProvider(read_deck(DECKS / "pp_deck.py").spec)
    reactants, products = _recombination(rrec, "[CH2]C(C)CC")
    both = _orientations(rrec, reactants, products, pp)
    assert _outcome(both["forward"]) == _outcome(both["reverse"]), both
    assert _outcome(both["forward"])[2], both


def _x_and_e():
    """E plus X = -CH(CH3)- (entry and exit on its methine), every step
    allowed."""
    mapping = polymer_keyword(DECKS / "ep_v2_deck.py", "polymer_spec")
    x = {"graph": {"atoms": [{"key": "x1", "element": "C", "hydrogens": 1,
                              "electronic_state": {"radical_electrons": 0}},
                             {"key": "x2", "element": "C", "hydrogens": 3,
                              "electronic_state": {"radical_electrons": 0}}],
                   "bonds": [{"atoms": ["x1", "x2"], "order": 1}],
                   "electronic_state": {"multiplicity": 1}},
         "sites": {"head": {"atom": "x1", "bond_order": 1, "class": "vinyl_head",
                            "vacancy": "forbidden"},
                   "tail": {"atom": "x1", "bond_order": 1, "class": "vinyl_tail",
                            "vacancy": "forbidden"}}}
    mapping["residues"] = {"E": mapping["residues"]["E"], "X": x}
    mapping["sequence_support"] = {"E": ["E", "X"], "X": ["E", "X"]}
    mapping["initial_ensemble"]["composition"] = {"E": {"numerator": 1, "denominator": 2},
                                                  "X": {"numerator": 1, "denominator": 2}}
    return sequence_provider(PolymerSpec.from_mapping(mapping))


def test_a_residue_with_one_backbone_atom_grows_at_both_ends_of_a_singleton(rrec):
    """Ethane read as one X (its methine capped twice): the long-chain
    analogues of CC <=> [H] + [CH2]C (the radical on the pendant carbon)
    grow at the head and the tail of that one atom at once.  XXX with the
    middle pendant radical, CCC(C)CC <=> [H] + [CH2]C(CC)CC, is one: built
    and matched, its insertion naming both words."""
    x = _x_and_e()
    reactants, products = _recombination(rrec, "[CH2]C")
    stub = ctest._stub_family(rrec, [cct._entry(900102, ["[H]", "[CH2]C(CC)CC"], "CCC(C)CC")])
    for name, ex in _orientations(stub, reactants, products, x).items():
        assert ex["verdict"] == cm.EXCEPT_MATCHED, (name, ex)
        assert [m["periods"] for m in ex["matches"]] == [
            [{"side": "both", "end": "X", "residues": [["X"], ["X"]]}]], (name, ex["matches"])
    (_port, coordinate), = [s for s in x.growth_sites(mol("CC", 0)) if s[0].atom_id == 0]
    (caps,) = [caps for side, residue, caps in coordinate if (side, residue) == ("both", "X")]
    assert len(set(caps)) == 2


def test_a_reading_needs_a_parse_of_its_species_in_the_declared_support(rrec, ep):
    """Under an alternating E/P support n-pentyl is E then P (the radical on
    E's head): its methyl at the far end is P's pendant, and read as an E
    tail its species does not parse (the residue before that E would be a P
    whose methine has no methyl).  So it is no E end, n-octyl is no analogue
    and an [H] + n-octyl <=> n-octane entry is not matched; P's tail is still
    read.  Under the statistical support the same methyl does parse as an E
    end (n-pentyl as E E and one atom of a residue cut by a scission), and
    the entry is matched."""
    alternating = cct._alternating()
    readings = {port.atom_id: sorted((side, residue) for side, residue, _c in coordinate)
                for port, coordinate in alternating.growth_sites(mol("[CH2]CCCC", 0))}
    assert readings == {3: [("tail", "P")]}, readings
    reactants, products = _recombination(rrec, "[CH2]CCCC")
    stub = ctest._stub_family(rrec, [cct._entry(900104, ["[H]", "[CH2]CCCCCCC"], "CCCCCCCC")])
    for name, ex in _orientations(stub, reactants, products, alternating).items():
        assert ex["verdict"] == "excluded:isomorphism", (name, ex)
    for name, ex in _orientations(stub, reactants, products, ep).items():
        assert ex["verdict"] == cm.EXCEPT_MATCHED, (name, ex)


BOUNDED = textwrap.dedent("""
    import json, resource, sys, time
    from types import SimpleNamespace
    sys.path.insert(0, {tests!r})
    from polymer_spec_v2_support import DECKS
    from rmgpy.kmc import certificates as cm
    from rmgpy.kmc import rmg_moves as rm
    from rmgpy.kmc.declared_sequence import sequence_provider
    from rmgpy.kmc.run import read_deck
    from rmgpy.molecule import Molecule
    from rmgpy.reaction import Reaction
    from rmgpy.species import Species

    def entry(index, reactants, products):
        def sp(smiles):
            return Species(molecule=[Molecule(smiles=smiles)])
        return SimpleNamespace(index=index, label=f"constructed-{{index}}", rank=None, data=None,
                               item=Reaction(reactants=[sp(s) for s in reactants],
                                             products=[sp(s) for s in products]))

    ethane = Molecule(smiles="CC")
    for k, atom in enumerate(ethane.atoms):
        atom.id = k
    rows = rm.molecule_rows(ethane)
    carbon, h = 0, next(a for a, _o in rows[0][5] if rows[a][0] == "H")
    radical = {{a: row for a, row in rows.items() if a != h}}
    radical[carbon] = rows[carbon][:3] + (1,) + rows[carbon][4:5] + \\
        (tuple(e for e in rows[carbon][5] if e[0] != h),)
    reactants = [rm.molecule_from_rows({{h: ("H", None, 0, 1, 0, ())}}),
                 rm.molecule_from_rows(radical)]
    products = [ethane]
    ep = sequence_provider(read_deck(DECKS / "ep_v2_deck.py").spec)
    family = SimpleNamespace(label="bounded", depositories=[SimpleNamespace(
        label="bounded/stub", entries={{900105: entry(900105, ["[H]", "[CH2]CCC"], ["CCCC"]),
                                       900106: entry(900106, ["C=C" * 30], ["C=C" * 30])}})])
    out = {{}}
    for name, a, b in (("forward", reactants, products), ("reverse", products, reactants)):
        t0 = time.time()
        ex = cm.depository_exclusion(family, a, b, repeat=ep.unit_formula(), growth=ep)
        out[name] = [ex["depository_max_heavy_atoms"], ex["verdict"],
                     sorted(m["index"] for m in ex.get("matches", ())), time.time() - t0]
    print(json.dumps(out))
""")


def test_the_words_are_built_only_for_a_candidate_within_its_sizes(tmp_path):
    """E/P (every step allowed) with an unrelated unsaturated C60 entry in
    the depository: H_max is 60, so every residue word up to 58 heavy atoms
    is a growth option (about 1.5e7 per reading).  Only the n-butane entry
    has a grown formula: the words are built for its formula and within its
    species sizes, so the exclusion finishes in seconds within 3 GB of
    address space (a subprocess) and matches it both ways."""
    tests = str(Path(__file__).parent)
    script = tmp_path / "bounded.py"
    script.write_text(BOUNDED.format(tests=tests))
    limit = 3 * 1024 ** 3

    def cap():
        import resource
        resource.setrlimit(resource.RLIMIT_AS, (limit, limit))

    t0 = time.time()
    done = subprocess.run([sys.executable, str(script)], cwd=tmp_path, capture_output=True,
                          text=True, timeout=180, preexec_fn=cap)
    assert done.returncode == 0, (time.time() - t0, done.stderr[-2000:])
    import json
    out = json.loads(done.stdout.strip().splitlines()[-1])
    for name, (hmax, verdict, matched, seconds) in out.items():
        assert hmax == 60, out
        assert (verdict, matched) == (cm.EXCEPT_MATCHED, [900105]), out
        assert seconds < 30, out
