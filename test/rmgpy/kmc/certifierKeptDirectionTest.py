"""Rate certification binds a family that is its own reverse in the direction RMG keeps
(I-254).

The compiler records the direction RMG's model keeps for such a reaction
(``CoreEdgeReactionModel.generate_kinetics``: an exact entry, then the lower rank, then
the exergonic direction at 298 K): when the reverse is kept, the rate is the reverse's
estimate (its template, entry and direction) and the reversed record owns it.  The
certifier must bind the pair to that same selection; binding the generated direction's
estimate instead refuses every kept-reverse pair, and an emitted catalog then reports
complete while those channels are absent.

An exactly thermoneutral reaction (dG298 = 0) is RMG's own tie: ``gibbs_is_positive =
G298 > -1e-8`` counts it endergonic, so the reverse of the generated direction is kept.
Compiler and certifier decide it through one function.
"""

import json
import math
import os
from pathlib import Path

import pytest

from rmgpy.kmc import compile as compile_module

REPOSITORY = Path(__file__).resolve().parents[3]
DECKS = REPOSITORY / "test/rmgpy/kmc/fixtures/decks"
DATABASE = Path(os.environ.get("RMG_DATABASE_PATH", "/nonexistent"))
needs_database = pytest.mark.skipif(
    not (DATABASE / "input/kinetics/families").is_dir(),
    reason="RMG_DATABASE_PATH is not a database checkout",
)
FAMILIES = ("H_Abstraction", "intra_H_migration")
SITES = ("interior_radical:0", "head_end_radical", "interior_radical:0+pristine", "doubly_featured",
         "interior_radical:0@5")


@pytest.fixture(scope="module")
def emitted(tmp_path_factory):
    """The PE toy deck's own-reverse records on a few sites, compiled, then certified as
    emission certifies them (``Emitter.init``)."""
    from rmgpy.kmc import emit
    from rmgpy.kmc.run import read_deck
    from rmgpy.kmc import compiler
    deck = DECKS / "pe_toy_deck.py"
    original = compile_module._proxy_inputs
    decide = getattr(compiler, "kept_own_reverse_kinetics", None)
    calls = {"compile": 0, "certify": 0}
    phase = ["compile"]

    def subset(request):
        proxies, source, options, discovery = original(request)
        return [p for p in proxies if p.site_type in SITES], source, options, discovery

    def counted(*args, **kwargs):
        calls[phase[0]] += 1
        return decide(*args, **kwargs)

    compile_module._proxy_inputs = subset
    if decide is not None:
        compiler.kept_own_reverse_kinetics = counted
    try:
        path, _ = compile_module.compile_artifact(compile_module.CompileRequest(
            database_path=str(DATABASE), output_dir=str(tmp_path_factory.mktemp("art")),
            families=FAMILIES, deck=str(deck), database_sha="kept-direction-binding"))
        phase[0] = "certify"
        emitter = emit.Emitter(path, DATABASE, {}, families=FAMILIES,
                               polymer=read_deck(deck).spec)
        emitter.init()
    finally:
        compile_module._proxy_inputs = original
        if decide is not None:
            compiler.kept_own_reverse_kinetics = decide
    return json.loads(Path(path).read_text())["records"], emitter.state, calls


def _pairs(records):
    """(as_generated record, rate-owning record) of every enabled own-reverse pair."""
    by_id = {r["event_id"]: r for r in records}
    out = []
    for record in records:
        if record["orientation"] != "as_generated" or record["status"] != "enabled":
            continue
        twin = by_id[record["reverse_of"]]
        owner = record if record["rate_source"]["kind"] == "RMG family estimate" else twin
        out.append((record, owner))
    return out


def _certified(state, records):
    """Original event id -> (status, reason) after emission's certification: a record
    that fails is republished refused under a new id (same reaction, same orientation)."""
    def key(r):
        return (r["family"], r["orientation"], r["template"], tuple(r["reactant_graphs"]),
                tuple(r["product_graphs"]))
    after = {key(r): (r["status"], r.get("status_reason", "")) for r in state["records"]}
    return {r["event_id"]: after[key(r)] for r in records}


DOMAIN = "not publishable"     # a channel outside the declared domain: refused either way


@needs_database
def test_kept_reverse_and_forward_kept_pairs_both_bind(emitted):
    records, state, _calls = emitted
    pairs = _pairs(records)
    kept = {"as_generated": 0, "reversed": 0}
    for _record, owner in pairs:
        kept[owner["orientation"]] += 1
    # the subset holds pairs RMG keeps in either direction
    assert kept["as_generated"] > 0 and kept["reversed"] > 0, kept
    refused = [(r["family"], r["template"], r["status_reason"][:160])
               for r in state["records"]
               if r["status"] == "refused" and "binding:" in r.get("status_reason", "")]
    assert not refused, refused
    # every pair certifies, or is refused by its domain only (whichever direction is kept)
    certified = _certified(state, records)
    bound = {"as_generated": 0, "reversed": 0}
    for record, owner in pairs:
        statuses = [certified[r["event_id"]] for r in (record, owner)]
        for status, reason in statuses:
            assert status == "enabled" or DOMAIN in reason, \
                (record["template"], owner["orientation"], reason[:160])
        if all(status == "enabled" for status, _reason in statuses):
            bound[owner["orientation"]] += 1
    assert bound["as_generated"] > 0 and bound["reversed"] > 0, bound


@needs_database
def test_an_exact_tie_keeps_the_reverse_of_the_generated_direction(emitted):
    """dG298 = 0 (dH298 = 0 and Kc = 1): RMG counts it endergonic, so the reverse is kept;
    the compiled pair says so and the certifier binds it."""
    records, state, _calls = emitted
    ties = []
    for record, owner in _pairs(records):
        provenance = record.get("thermo_provenance") or \
            next(r for r in records if r["event_id"] == record["reverse_of"])["thermo_provenance"]
        kc = provenance["equilibrium_constant_table"]["Kc"]
        if provenance["reaction_enthalpy_298_J_per_mol"] == 0 and \
                all(abs(math.log(k)) < 1e-9 for k in kc):
            ties.append((record, owner))
    assert ties
    certified = _certified(state, records)
    bound = 0
    for record, owner in ties:
        rs = owner["rate_source"]
        assert owner["orientation"] == "reversed", (record["template"], owner["orientation"])
        assert rs["family_template_direction"] == "reverse"
        statuses = [certified[r["event_id"]] for r in (record, owner)]
        for status, reason in statuses:
            assert status == "enabled" or DOMAIN in reason, (record["template"], reason[:160])
        bound += all(status == "enabled" for status, _reason in statuses)
    assert bound > 0


@needs_database
def test_compiler_and_certifier_decide_through_one_function(emitted):
    """The compiler and rate certification both call
    ``compiler.kept_own_reverse_kinetics`` (looked up at call time), never a copy."""
    _records, _state, calls = emitted
    assert calls["compile"] > 0 and calls["certify"] > 0, calls
