"""Emission speed without output change (I-235).

Three repeated computations of template emission, each pinned on a real
case with a bound and checked against the base code's own outputs:

* ``window_matcher.match`` evaluated every embedding of a window although
  an outcome is decided by the edited atoms' images, the touched components
  and the subsystem forms: embeddings that differ only elsewhere
  (interchangeable hydrogens, context atoms) repeat it.  The corpus is every
  emitted template of a PE and a PP H_Abstraction sample on its own emitted
  witness, with ``repr`` of the transitions and classes at the base code.
* ``certificates.depository_exclusion`` scanned the list of grown side
  formulas once per depository entry and orientation.  The corpus is its
  inputs met while emitting a PE and a PP record, with the base's results.
* Emission's certifier kept 8 channel sets, fewer than an artifact's
  certification revisits, so it rebuilt most of them.
"""

import gzip
import json
import random
import time
from pathlib import Path

import pytest

from rmgpy.kmc import certificates as cm
from rmgpy.kmc import channels as chm
from rmgpy.kmc import emit
from rmgpy.kmc import template_emission as te
from rmgpy.kmc import window_matcher as wm
from rmgpy.molecule import Molecule

import rmg_moves_harness as harness
from channelsTest import mol

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def _corpus(path):
    with gzip.open(path, "rt") as handle:
        return json.load(handle)["cases"]


MATCH_CASES = _corpus(FIXTURES / "window_matcher" / "emitted_witnesses.json.gz")
EXCLUSION_CASES = [c for name in ("pe_h_abstraction", "pp_h_abstraction")
                   for c in _corpus(FIXTURES / "depository_exclusion" / f"{name}.json.gz")]


def _matcher_case(case):
    template = te.matcher_from_payload(case["template"], 700.0, weight=1.0)
    return template, te.melt_state(te.rows_from_json(case["rows"]))


def _cpu(fn, repeat=1):
    start = time.process_time()
    for _ in range(repeat):
        fn()
    return time.process_time() - start


# ---------------------------------------------------------------------------
# window_matcher.match
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("case", MATCH_CASES, ids=lambda c: f"{c['record']}-{c['template']['name'][-8:]}")
def test_match_gives_the_base_transitions_and_classes(case):
    template, state = _matcher_case(case)
    result = wm.match([template], state)
    assert repr(result.transitions) == case["transitions"]
    assert repr(result.classes) == case["classes"]


def test_match_evaluates_each_edit_site_and_form_choice_once(monkeypatch):
    """The slowest emitted witness: 40960 embeddings of 10 transitions at the
    base (two hydrogens per CH2 of the window, every swap an embedding)."""
    case = max(MATCH_CASES, key=lambda c: c["base_evaluated"])
    template, state = _matcher_case(case)
    calls = []
    real = wm._evaluate

    def counting(view, tmpl, index, mapping, comps, choice):
        calls.append((index, comps, choice, tuple(mapping[t] for t in tmpl.edited_all)))
        return real(view, tmpl, index, mapping, comps, choice)

    monkeypatch.setattr(wm, "_evaluate", counting)
    result = wm.match([template], state)
    assert len(calls) == len(set(calls)) == case["distinct_sites"]
    assert result.classes[template.name].evaluated == case["base_evaluated"]


def test_match_on_the_slowest_emitted_witness_is_within_its_time_bound():
    case = max(MATCH_CASES, key=lambda c: c["base_cpu_s"])
    template, state = _matcher_case(case)
    wm.match([template], state)   # warm the view
    assert _cpu(lambda: wm.match([template], te.melt_state(te.rows_from_json(case["rows"]))),
                repeat=3) < 0.9


# ---------------------------------------------------------------------------
# certificates.depository_exclusion
# ---------------------------------------------------------------------------
pytestmark_db = pytest.mark.skipif(not harness.DATABASE.is_dir(),
                                   reason=f"needs the RMG database at {harness.DATABASE}")


@pytest.fixture(scope="module")
def families():
    return emit.load_database(["H_Abstraction"], harness.DATABASE).kinetics.families


def _exclusion_args(case):
    reactants = [Molecule().from_adjacency_list(a) for a in case["reactants"]]
    products = [Molecule().from_adjacency_list(a) for a in case["products"]]
    return reactants, products, tuple(tuple(x) for x in case["repeat"])


@pytestmark_db
@pytest.mark.parametrize("case", EXCLUSION_CASES, ids=lambda c: str(c["base_cpu_s"]))
def test_depository_exclusion_gives_the_base_result(families, case):
    reactants, products, repeat = _exclusion_args(case)
    out = cm.depository_exclusion(families[case["family"]], reactants, products, cache={},
                                  repeat=repeat)
    assert json.dumps(out, sort_keys=True, default=str) == case["result"]


def _list_scan(rows, rs, ps):
    """The base's scan: for every entry and orientation, the growths of each
    side with the entry's formulas, listed."""
    out = []
    for label, index, er, ep in rows:
        for a, b in ((er, ep), (ep, er)):
            hit_r = [g for g, f in rs if f == a]
            hit_p = [g for g, f in ps if f == b]
            if hit_r and hit_p and (any(hit_r) or any(hit_p)):
                out.append({"depository": label, "index": index})
                break
    return out


@pytest.mark.parametrize("seed", range(40))
def test_formula_candidates_are_the_list_scans(seed):
    rng = random.Random(seed)
    formulas = [((("C", c), ("H", h)),) for c in range(1, 4) for h in range(0, 5)]
    formulas += [tuple(sorted(rng.sample(formulas, 2))) for _ in range(6)]

    def side():
        return {(rng.random() < 0.5, rng.choice(formulas)) for _ in range(rng.randrange(0, 12))}

    rows = [("training", k, rng.choice(formulas), rng.choice(formulas))
            for k in range(rng.randrange(1, 60))]
    rs, ps = side(), side()
    assert cm.formula_candidates(rows, rs, ps) == _list_scan(rows, rs, ps)


@pytestmark_db
def test_depository_exclusion_of_the_slowest_witness_is_within_its_time_bound(families):
    """A PE witness: CH2 growth up to the 22 heavy atoms of the largest
    depository species gives hundreds of grown side formulas per side; the
    base scanned them for every entry and orientation (about 70 ms a call)."""
    case = max((c for c in EXCLUSION_CASES if c["repeat"] == [["C", 1], ["H", 2]]),
               key=lambda c: c["base_cpu_s"])
    reactants, products, repeat = _exclusion_args(case)
    family = families[case["family"]]
    cache = {}
    cm.depository_exclusion(family, reactants, products, cache=cache, repeat=repeat)
    assert _cpu(lambda: cm.depository_exclusion(family, reactants, products, cache=cache,
                                                repeat=repeat), repeat=20) < 0.4


# ---------------------------------------------------------------------------
# emission's certifier: channel sets
# ---------------------------------------------------------------------------
@pytestmark_db
def test_emission_certifier_builds_each_channel_set_of_a_round_robin_once(families,
                                                                        monkeypatch):
    """Twelve reactant pairs certified twice in turn (an artifact's records
    revisit their reactant sets after other records): each channel set is
    built once."""
    built = []
    real = chm.build_channels

    def counting(family, roots, direction, **kw):
        built.append(1)
        return real(family, roots, direction, **kw)

    monkeypatch.setattr(chm, "build_channels", counting)
    certifier = emit.make_certifier(families, {})
    family = families["H_Abstraction"]
    donors = ["C", "CC", "CCC", "CC(C)C"]
    radicals = ["[CH3]", "C[CH2]", "CC[CH2]"]
    pairs = [(d, r) for d in donors for r in radicals]
    assert len(pairs) == 12
    for _ in range(2):
        for d, r in pairs:
            certifier.channel_set(family, [mol(d, 0), mol(r, 100)])
    assert len(built) == len(pairs)
