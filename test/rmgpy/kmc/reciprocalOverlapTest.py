"""Loaded templates that make one transition twice (I-210).

The fixture is the six H_Abstraction payloads of the three PE reciprocal pairs
the I-199 PE catalog refused (``fixtures/i210_pe_exchange_payloads.zlib.b64``,
copied unchanged from ``pe.rows.jsonl``): a secondary radical and a secondary
C-H exchange their H, both sites at least five carbons from a chain end.  The
balance checks run with the hydrogen-leaf quotient (``loaded_transitions``),
which these tests hold equal to the matcher's full enumeration.
"""

import base64
import json
import zlib
from pathlib import Path

import pytest

from rmgpy.kmc import template_emission as te
from rmgpy.kmc import window_matcher as wm

HERE = Path(__file__).resolve().parent
DATABASE = Path("/home/alon/runs/phase2b/database")
EXCHANGE = ("3e484af0e1", "53b80ec7cc", "a7562b8ed1", "247abac61b")


def exchange_payloads():
    path = HERE / "fixtures" / "i210_pe_exchange_payloads.zlib.b64"
    return json.loads(zlib.decompress(base64.b64decode(path.read_bytes())))


def bucket(payloads, keys=EXCHANGE):
    """The loaded H_Abstraction templates whose key starts with one of ``keys``,
    as matcher templates named by their key prefix."""
    loaded = te._Loaded(payloads)
    entries = [t for t in loaded.by[("H_Abstraction", "forward")] if t["key"][:10] in keys]
    return [te.matcher_from_payload(dict(t, family="H_Abstraction"), 0.0, weight=1.0,
                                    name=t["key"][:10]) for t in entries]


def witness_state(payloads, key, side):
    t = next(t for p in payloads for t in p["templates"] if t["key"].startswith(key))
    return te.melt_state(te.rows_from_json(t["maps"]["witnesses"]["emitted"][side]))


def transition_fields(tr):
    """Everything a balance check reads of a transition (not the embedding
    count, nor the mapping of unedited hydrogen leaves)."""
    return (tr.template, tr.family, tr.class_index, tr.edit_key, tr.components, tr.edits,
            tr.form_bonds, tr.form_radicals, tr.identity, tr.delta, tr.weight,
            tr.state_revision, tr.state_hash)


@pytest.mark.parametrize("key,side", [("3e484af0e1", "r"), ("53b80ec7cc", "r")])
def test_loaded_transitions_equal_the_matchers_on_pe_witnesses(key, side):
    payloads = te.resolve_collisions(exchange_payloads())[0]
    templates = bucket(payloads)
    state = witness_state(payloads, key, side)
    full = wm.match(templates, state).transitions
    quotient = te.loaded_transitions(templates, state)
    assert len(full) > 0
    assert [transition_fields(x) for x in quotient] == [transition_fields(x) for x in full]
    assert sum(x.embeddings for x in quotient) < sum(x.embeddings for x in full)
