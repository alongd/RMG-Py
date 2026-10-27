"""The melt engine makes a transition two deduplicated templates share once (I-271).

``resolve_overlaps`` deduplicates a class overlap nested neither way at one
weight: both sides load, and the reversed side's templates defer to the
as_generated ones (``defers_to``, ``window_matcher.Template.defers``).  The
engine indexes per component: a unimolecular class's transitions, and a
bimolecular class's sites per participant, paired across components.  Its
counts must be the deferred transition set's (``loaded_transitions``): a
deferring class indexes only the transitions no class it defers to makes.
The fixtures are the two PP components the assembly tests use: an
H_Abstraction one (bimolecular) and an intra_H_migration one (unimolecular).
"""

import gzip
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import template_emission as te

import melt_engine_support as S

FIXTURES = Path(__file__).resolve().parent / "fixtures"
CASES = ("i251/neither_nested_cascade_payloads.json.gz", "i249/neither_nested_overlap.json.gz")
T = 700.0


def _deduplicated(name):
    with gzip.open(FIXTURES / name, "rt") as handle:
        doc = json.load(handle)
    payloads = doc["payloads"] if isinstance(doc, dict) else doc
    merged, _ = te.resolve_collisions(payloads)
    out, table = te.resolve_overlaps(merged)
    (row,) = [r for r in table if r["disposition"] == "deduplicated"]
    return out, row


def _states(out):
    """Every stored witness state of the fixture (both sides), once."""
    states, seen = [], set()
    for p in out:
        for t in p["templates"]:
            for wit in ((t.get("maps") or {}).get("witnesses") or {}).values():
                for side in ("r", "p"):
                    rows = te.rows_from_json(wit[side])
                    sig = te.rm._signature(rows)
                    if sig not in seen:
                        seen.add(sig)
                        states.append(rows)
    return states


def engine_counts(engine) -> Counter:
    """Transitions per template the engine's index holds (allowed and leak)."""
    out = Counter()
    for cls in engine.classes:
        name = engine.templates[cls.index].name
        if not cls.bimolecular:
            out[name] += engine.uni_allowed[cls.index] + engine.uni_leak[cls.index]
            continue
        rows = [(k0, k1) for (ci, k0, k1) in engine.channels if ci == cls.index] if \
            not cls.met_forward else \
            [(k0, k1) for (ci, p, k0) in engine.S if ci == cls.index and p == 0
             for (cj, q, k1) in engine.S if cj == ci and q == 1]
        pairs = sum(engine.S[(cls.index, 0, k0)] * engine.S[(cls.index, 1, k1)]
                    - engine.X[(cls.index, k0, k1)] for k0, k1 in rows)
        out[name] += pairs / cls.divisor
    return +out


@pytest.mark.parametrize("name", CASES)
def test_the_engine_indexes_each_shared_transition_once(name):
    out, row = _deduplicated(name)
    catalog = me.PayloadCatalog(out)
    shared = 0
    for rows in _states(out):
        engine = me.MeltEngine(S.ps_spec({2: 1}), catalog, T, np.random.default_rng(0),
                               state=te.melt_state(rows))
        made = te.loaded_transitions(engine.templates, engine.state)
        alone = Counter(tr.identity[1:] for tr in te.loaded_transitions(
            [te.wm.Template(t.name, t.family, t.atoms, list(t.bonds), t.edits, weight=t.weight)
             for t in engine.templates], engine.state))
        shared += sum(1 for n in alone.values() if n > 1)
        assert engine_counts(engine) == +Counter(tr.template for tr in made)
    assert shared


def test_events_on_shared_sites_keep_the_index_a_full_recompute():
    out, row = _deduplicated(CASES[0])
    catalog = me.PayloadCatalog(out)
    fired = 0
    for rows in _states(out):
        engine = me.MeltEngine(S.ps_spec({2: 1}), catalog, T, np.random.default_rng(3),
                               state=te.melt_state(rows))

        def check(engine, _event):
            engine.check_consistency()
            assert engine.snapshot() == engine.fresh_snapshot()
        fired += engine.run(max_events=6, check=check).events
    assert fired
