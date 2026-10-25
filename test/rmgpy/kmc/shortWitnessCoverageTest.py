"""The compile witnesses every pristine chain shorter than its proxies
(I-261), so a two-fragment event whose fragments both end inside the window
has a template, and the PE toy catalog joins every radical pair and cuts
every chain that its homolysis makes from the shortest fragment.

A template's window closes on a side only where the witness molecule ends
inside it.  The shortest pristine witness was the proxy itself (3 repeat
units, n-hexane for PE), so no template cut n-butane or n-pentane and none
joined ethyl + ethyl or ethyl + propyl, which RMG recombines like any
other pair.

Ethyl + ethyl and n-butane were the last gaps (I-262): that closed context
occurs on n-butane alone, which RMG rates from a training entry (a
whole-species match), and no rule template makes it.  The depository
reaction is now its own template (a closed context), and so is every
whole species the depository rates whose cut leaves a methyl
(``METHYL``: ethane, propane, n-decane), the only methyl cuts: elsewhere a
reacting atom on an artificial port stays a named refusal.  Primary
radical + ethylene is the open case: every witness up to octyl + C2H4 is a
training molecule, so the member grows to the first size the depository
does not match and the rule template formed there carries the training
sizes as overrides, in both directions.

``fixtures/i261_pe_toy_recombination_payloads.json`` holds the PE toy
catalog's R_Recombination payloads whose templates break or form a C-C bond
(``fixtures/decks/pe_toy_deck.py``, emitted on the stage-A artifact of the
code under test); ``fixtures/i262_pe_toy_addition_payloads.json`` its
R_Addition_MultipleBond payloads of the same kind;
``fixtures/i261_pe_toy_split_kernels_700K.json`` the homolysis and
end-recombination kernels of that catalog's kernel-rates table at 700 K.
Sizes count CH2 periods.
"""

import json
import os
from pathlib import Path

import numpy
import pytest

import melt_engine_support as S
from momProducerTest import LABEL

from rmgpy.kmc import emit
from rmgpy.kmc import kernel_rates as K
from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import mom_producer as mp
from rmgpy.kmc import run as kr
from rmgpy.kmc.declared_chain import ChainSpec
from rmgpy.molecule import Molecule
from rmgpy.polymer_spec import LinearHomopolymerSpec, normalized_linear_homopolymer_spec

FIXTURES = Path(__file__).parent / "fixtures"
DECKS = FIXTURES / "decks"
PAYLOADS = FIXTURES / "i261_pe_toy_recombination_payloads.json"
ADDITION = FIXTURES / "i262_pe_toy_addition_payloads.json"
TABLE = FIXTURES / "i261_pe_toy_split_kernels_700K.json"
H = "homolysis:substrate=chain"
C = "recombination:radical=radical_end:CH2+radical=radical_end:CH2"
LONGEST = 12
# ethyl: the shortest fragment a chain's homolysis leaves (module doc)
SHORTEST = 2
# n-alkane -> its methyl cuts (left fragment sizes): whole species RMG rates
# from its depository (module doc)
METHYL = {2: (1,), 3: (1, 2), 10: (1, 9)}
# the shortest end radical witnessed with ethylene (``pair:head_end_radical_short+monomer``)
SHORTEST_ADDING = 4
# the depository reactions of the PE toy that join or cut a C-C bond and
# that no loaded template carried before I-262 (both directions each); the
# training pairs of n-decane and n-dodecane already rode a longer member's
# rule template
DEPOSITORY_ONLY = (
    "C2H5 + C2H5 <=> C4H10", "CH3 + CH3 <=> C2H6", "CH3 + C2H5 <=> C3H8", "CH3 + C9H19 <=> C10H22",
    "C4H9-1 + C2H4 <=> C6H13-10", "C5H11-10 + C2H4 <=> C7H15-9", "C6H13-11 + C2H4 <=> C8H17-7",
    "C7H15-10 + C2H4 <=> C9H19-4", "C8H17-8 + C2H4 <=> C10H21-2", "C10H21-3 + C2H4 <=> C12H25")
DATABASE = Path(os.environ.get("RMG_DATABASE_PATH", "/nonexistent"))
needs_database = pytest.mark.skipif(not (DATABASE / "input/kinetics/families").is_dir(),
                                    reason="RMG_DATABASE_PATH is not a database checkout")


def _pe():
    return normalized_linear_homopolymer_spec(LinearHomopolymerSpec.from_mapping(S.pe_mapping({10: 1})))


def _alkane(n):
    return Molecule(smiles="C" * n)


def _radical(n):
    return "[CH3]" if n == 1 else "[CH2]" + "C" * (n - 1)


# -- the compile's witnesses ---------------------------------------------------


@pytest.mark.parametrize("deck, units", [("pe_toy_deck.py", 3), ("pp_deck.py", 3), ("ps_deck.py", 3)])
def test_the_declared_proxies_witness_a_pristine_chain_of_every_shorter_period_count(deck, units):
    from rmgpy.kmc.declared_proxies import _Builder, declared_proxy_set
    spec = kr.read_deck(DECKS / deck).spec
    declared = _Builder.of(spec, units).declared
    singles = [p for p in declared_proxy_set(spec, units) if len(p.reactants) == 1]
    for k in range(1, units * declared.periods_per_repeat):
        chain = declared.build("expected", ChainSpec(k)).molecule
        found = [p for p in singles if p.reactants[0].molecule[0].is_isomorphic(chain, strict=False)]
        assert len(found) == 1, (deck, k, [p.site_type for p in found])
        (proxy,) = found
        assert proxy.participant_site_types == ("pristine",)
        assert proxy.metadata["context_class"] == "pristine"
        # artificial ends, as every other proxy: a reacting end atom stays a cut point
        assert sorted(port.orientation for port in proxy.artificial_boundaries) == ["head", "tail"]


def test_the_pe_short_witnesses_are_the_n_alkanes_below_n_hexane():
    from rmgpy.kmc.declared_proxies import declared_proxy_set
    singles = [p.reactants[0].molecule[0] for p in declared_proxy_set(_pe(), 3) if len(p.reactants) == 1]
    for n in range(1, 6):
        assert sum(m.is_isomorphic(_alkane(n), strict=False) for m in singles) == 1, n


def test_the_roster_offers_each_short_witness_as_a_system_of_its_own():
    roster = emit.declared_reference_systems(_pe(), 3)
    chains = sorted(entry["declared_chain"]["units"] for (entry,) in
                    (v for k, v in roster.systems.items() if k.startswith("single:")) if
                    "declared_chain" in entry and not entry["declared_chain"]["edits"] and
                    entry["declared_chain"]["head"] == entry["declared_chain"]["tail"] == "artificial")
    # 1..5 periods, the 3-unit pristine chain (6) and its 5-unit copy (10)
    assert chains == [1, 2, 3, 4, 5, 6, 10]


# -- the PE toy catalog ----------------------------------------------------------


def _engine(path):
    payloads = json.loads(path.read_text())["payloads"]
    spec = kr.read_deck(DECKS / "pe_toy_deck.py").spec
    return me.MeltEngine(K._single_chain_spec(spec, 4), me.PayloadCatalog(payloads), 700.0,
                         numpy.random.default_rng(0))


@pytest.fixture(scope="module")
def engine():
    return _engine(PAYLOADS)


@pytest.fixture(scope="module")
def adding():
    return _engine(ADDITION)


def _cuts(engine, n):
    """left fragment size -> times a C-C bond of n-alkane n is cut."""
    state = me._rows_state(S.smiles_rows(["C" * n]))
    (cid,) = list(state.components.ids())
    entry = engine._scan([cid], state=state)[cid]
    out = {}
    # a matched size fires its depository override in place of the template (uni_o)
    for allowed, _leak in list(entry.uni.values()) + list(entry.uni_o.values()):
        for transition in allowed:
            for edit in transition.edits:
                if type(edit).__name__ == "SetBond" and edit.new == 0.0 and max(edit.a, edit.b) < n:
                    left = min(edit.a, edit.b) + 1
                    out[left] = out.get(left, 0) + 1
    return out


def _pair_sites(engine, a_smiles, b_smiles):
    """Unordered pair sites, summed over the bimolecular templates (each
    over its divisor), of molecules ``a_smiles`` and ``b_smiles``."""
    state = me._rows_state(S.smiles_rows([a_smiles, b_smiles]))
    a, b = list(state.components.ids())
    entries = engine._scan([a, b], state=state)
    total = 0.0
    for cls in engine.classes:
        if not cls.bimolecular:
            continue
        s = {(cid, p): sum(len(v) for (ci, q, _k), v in entries[cid].sites.items() if ci == cls.index and q == p)
             for cid in (a, b) for p in (0, 1)}
        total += (s[(a, 0)] * s[(b, 1)] + s[(b, 0)] * s[(a, 1)]) / cls.divisor
    return total


def _scissions(engine, n):
    """(rule, override) counts of the C-C bonds of n-alkyl radical n that
    the loaded templates break: at the rule weight (uni) and at a depository
    override (uni_o)."""
    rows = S.smiles_rows([_radical(n)])
    state = me._rows_state(rows)
    (cid,) = list(state.components.ids())
    entry = engine._scan([cid], state=state)[cid]

    def count(items):
        return sum(1 for allowed, _leak in items for tr in allowed
                   if any(type(ed).__name__ == "SetBond" and ed.new == 0.0 and
                          rows[ed.a][0] == rows[ed.b][0] == "C" for ed in tr.edits))
    return count(entry.uni.values()), count(entry.uni_o.values())


def test_every_chain_is_cut_once_at_every_bond_from_the_shortest_fragment(engine):
    """Every C-C bond whose fragments are at least ethyl, from n-butane on,
    and the methyl cuts of the depository's whole species, nothing else."""
    wrong = {}
    for n in range(2, LONGEST + 1):
        want = {left: 1 for left in range(SHORTEST, n - SHORTEST + 1)}
        want.update({left: 1 for left in METHYL.get(n, ())})
        got = _cuts(engine, n)
        if got != want:
            wrong[n] = got
    assert wrong == {}, f"chains cut otherwise: {wrong}"


def test_every_radical_pair_from_the_shortest_fragment_is_joined_once(engine):
    """Every pair of n-alkyl radicals from ethyl + ethyl; methyl only with
    the partners of the depository's whole species."""
    sizes = range(SHORTEST, 10)
    wrong = {(i, j): _pair_sites(engine, _radical(i), _radical(j)) for i in sizes for j in sizes if i <= j}
    wrong = {pair: got for pair, got in wrong.items() if got != pytest.approx(1.0)}
    assert wrong == {}, f"pairs joined otherwise: {wrong}"
    with_methyl = {n - left for n, lefts in METHYL.items() for left in lefts[:1]}
    got = {j: _pair_sites(engine, _radical(1), _radical(j)) for j in range(1, 10)}
    assert got == {j: pytest.approx(1.0 if j in with_methyl else 0.0) for j in range(1, 10)}


def test_every_end_radical_from_butyl_adds_ethylene_and_releases_it_back(adding):
    """Primary radical + C2H4 (both ends of the alkene) from butyl on, and
    the beta-scission back from hexyl on: each n-alkyl radical breaks one C-C
    bond, ethylene off its radical end."""
    far = _pair_sites(adding, _radical(LONGEST), "C=C")
    assert far > 0.0
    got = {i: _pair_sites(adding, _radical(i), "C=C") for i in range(1, LONGEST + 1)}
    assert got == {i: pytest.approx(far if i >= SHORTEST_ADDING else 0.0) for i in got}
    for n in range(SHORTEST_ADDING + 2, LONGEST + 3):
        assert sum(_scissions(adding, n)) == 1, (n, _scissions(adding, n))


@needs_database
def test_each_whole_species_cut_is_a_reaction_rmg_rates_only_from_its_depository():
    """The closed contexts (n-butane, and the methyl cuts) occur on their own
    molecules only, and RMG rates those from a training entry (a
    whole-species match): the reason each is its own template."""
    from rmgpy.kmc import certificates as cm
    db = emit.load_database(["R_Recombination"], DATABASE)
    family = db.kinetics.families["R_Recombination"]
    for n, lefts in sorted({4: (SHORTEST,), **METHYL}.items()):
        radicals = [Molecule(smiles=_radical(k)) for k in (lefts[0], n - lefts[0])]
        selected = cm.depository_search(family, radicals, [_alkane(n)])["selected"]
        assert selected is not None and selected["depository"].endswith("training"), n


@needs_database
def test_an_ethylene_scission_takes_its_override_exactly_where_rmg_rates_it_from_its_depository(adding):
    """At a size RMG rates from training the scission fires the depository
    override in place of the rule template; at every other size the rule
    weight."""
    from rmgpy.kmc import certificates as cm
    db = emit.load_database(["R_Addition_MultipleBond"], DATABASE)
    family = db.kinetics.families["R_Addition_MultipleBond"]
    for n in range(SHORTEST_ADDING + 2, LONGEST + 3):
        selected = cm.depository_search(family, [Molecule(smiles=_radical(n - 2)), Molecule(smiles="C=C")],
                                        [Molecule(smiles=_radical(n))])["selected"]
        rule, override = _scissions(adding, n)
        assert (rule, override) == ((0, 1) if selected is not None else (1, 0)), (n, selected)


@pytest.mark.parametrize("label", DEPOSITORY_ONLY)
def test_every_depository_only_reaction_rides_a_template_in_both_directions(label):
    """Each depository reaction no loaded template carried (its sizes are
    all training molecules, or it is a whole species) is an override of a
    loaded template, forward and reverse: RMG's depository rate on the
    catalog, never a silent skip."""
    carried = {(o["entry"]["label"], o["direction"])
               for path in (PAYLOADS, ADDITION) for p in json.loads(path.read_text())["payloads"]
               for t in p["templates"] for o in t.get("overrides") or ()}
    assert {d for name, d in carried if name == label} == {"forward", "reverse"}, label


def _methyl_cut(term) -> bool:
    return any(cls == "radical_end:CH3" for outs in term["law"]["exact"].values()
               for _m, prods in outs for cls, _n in prods)


def test_grouped_recombination_is_refused_on_the_whole_unit_kernel_states():
    """MoM's kernel mode on homolysis + end recombination: (2, 2) is joined
    (no pair-site refusal any more).  The methyl fragments of the
    depository's whole species are size departures (I-266: they cut single
    sizes only); this table, built on whole-unit states (before I-266), then
    refuses on the kernel-rates states: homolysis's law holds from 6 periods
    (no state of n-pentane), ethyl + ethyl makes 4.  On the period grid it
    maps (``endClassPoolsTest``)."""
    doc = json.loads(TABLE.read_text())
    spec = kr.read_deck(DECKS / "pe_toy_deck.py").spec
    terms = doc["kernels"][H]["templates"]
    assert 0 < sum(map(_methyl_cut, terms)) < len(terms)
    with pytest.raises(mp.ProducerError, match=r"radicals from 2 periods join into chains of 4, shorter than 6, "
                                               r"where the homolysis kernel's law holds$"):
        mp.kernel_input(doc, spec, [H, C], temperature=700.0, phi=1.0, label=LABEL,
                        radical_label=LABEL + "_re", volume_ratio=1.0)
