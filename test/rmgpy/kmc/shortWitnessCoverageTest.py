"""The compile witnesses every pristine chain shorter than its proxies
(I-261), so a two-fragment event whose fragments both end inside the window
has a template, and the PE toy catalog joins every radical pair and cuts
every chain that its homolysis makes, but the named gaps.

A template's window closes on a side only where the witness molecule ends
inside it.  The shortest pristine witness was the proxy itself (3 repeat
units, n-hexane for PE), so no template cut n-butane or n-pentane and none
joined ethyl + ethyl or ethyl + propyl, which RMG recombines like any
other pair.

The named gaps (``UNJOINED``, ``UNCUT``): ethyl + ethyl and n-butane.  That
closed context occurs on n-butane alone, which RMG rates from a training
entry (a whole-species match): emission refuses a whole-species match (its
rate rides as an override on the rule template that makes it), and no rule
template makes this one.  A fix that covers them must shrink the sets.

``fixtures/i261_pe_toy_recombination_payloads.json`` holds the PE toy
catalog's R_Recombination payloads whose templates break or form a C-C bond
(``fixtures/decks/pe_toy_deck.py``, emitted on the stage-A artifact of the
code under test); ``fixtures/i261_pe_toy_split_kernels_700K.json`` the
homolysis and end-recombination kernels of that catalog's kernel-rates table
at 700 K.  Sizes count CH2 periods.
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
TABLE = FIXTURES / "i261_pe_toy_split_kernels_700K.json"
H = "homolysis:substrate=chain"
C = "recombination:radical=radical_end:CH2+radical=radical_end:CH2"
LONGEST = 12
# a closed context whose only molecules RMG rates from a depository (module doc)
UNJOINED = {(2, 2)}
UNCUT = {4}
DATABASE = Path(os.environ.get("RMG_DATABASE_PATH", "/nonexistent"))
needs_database = pytest.mark.skipif(not (DATABASE / "input/kinetics/families").is_dir(),
                                    reason="RMG_DATABASE_PATH is not a database checkout")


def _pe():
    return normalized_linear_homopolymer_spec(LinearHomopolymerSpec.from_mapping(S.pe_mapping({10: 1})))


def _alkane(n):
    return Molecule(smiles="C" * n)


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


@pytest.fixture(scope="module")
def engine():
    payloads = json.loads(PAYLOADS.read_text())["payloads"]
    spec = kr.read_deck(DECKS / "pe_toy_deck.py").spec
    return me.MeltEngine(K._single_chain_spec(spec, 4), me.PayloadCatalog(payloads), 700.0,
                         numpy.random.default_rng(0))


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


def _pair_sites(engine, i, j):
    """Unordered pair sites, summed over the recombination templates (each
    over its divisor), of n-alkyl radicals of i and j periods."""
    state = me._rows_state(S.smiles_rows(["[CH2]" + "C" * (i - 1), "[CH2]" + "C" * (j - 1)]))
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


def _shortest_fragment(engine):
    return min(left for n in range(2, LONGEST + 1) for left in _cuts(engine, n))


def test_every_chain_is_cut_once_at_every_bond_from_the_shortest_fragment_but_the_named_gaps(engine):
    e = _shortest_fragment(engine)
    assert e >= 1
    wrong = {n: _cuts(engine, n) for n in range(2 * e, LONGEST + 1)}
    wrong = {n: got for n, got in wrong.items() if got != {left: 1 for left in range(e, n - e + 1)}}
    assert wrong == {n: {} for n in UNCUT}, f"shortest fragment {e}; chains cut otherwise: {wrong}"


def test_every_radical_pair_from_the_shortest_fragment_is_joined_once_but_the_named_gaps(engine):
    e = _shortest_fragment(engine)
    sizes = range(e, 9)
    wrong = {(i, j): _pair_sites(engine, i, j) for i in sizes for j in sizes if i <= j}
    wrong = {pair: got for pair, got in wrong.items() if got != pytest.approx(1.0)}
    assert wrong == {pair: 0.0 for pair in UNJOINED}, f"shortest fragment {e}; pairs joined otherwise: {wrong}"


@needs_database
def test_each_named_gap_is_a_reaction_rmg_rates_only_from_its_depository():
    """The closed context of each gap occurs on its own molecules only, and
    RMG rates those from a training entry (a whole-species match), so no
    rule template exists for an override to ride on."""
    from rmgpy.kmc import certificates as cm
    db = emit.load_database(["R_Recombination"], DATABASE)
    family = db.kinetics.families["R_Recombination"]
    for (i, j), n in zip(sorted(UNJOINED), sorted(UNCUT)):
        assert i + j == n
        radicals = [Molecule(smiles="[CH2]" + "C" * (k - 1)) for k in (i, j)]
        selected = cm.depository_search(family, radicals, [_alkane(n)])["selected"]
        assert selected is not None and selected["depository"].endswith("training"), (i, j)


def test_grouped_recombination_is_refused_on_the_named_gaps_only():
    """MoM's kernel mode with homolysis making ethyl radicals (min_fragment
    2): the recombination law fails on (2, 2) alone; (2, 3) is joined."""
    doc = json.loads(TABLE.read_text())
    with pytest.raises(mp.ProducerError, match=r"pair sites \(2, 2\): 0 \(0: no template joins the pair\)$"):
        mp.kernel_input(doc, kr.read_deck(DECKS / "pe_toy_deck.py").spec, [H, C], temperature=700.0, phi=1.0,
                        label=LABEL, radical_label=LABEL + "_re", volume_ratio=1.0)
