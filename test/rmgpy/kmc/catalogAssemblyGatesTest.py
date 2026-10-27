"""The gates a polymer catalog with secondary, tertiary and resonance sites met at
assembly (I-249), one root cause each.

* A reaction of a family that is its own reverse has one kinetics in RMG: the model
  (``CoreEdgeReactionModel.generate_kinetics``) estimates both directions and keeps one
  (an exact entry, then the lower rank, then the exergonic direction).  A compiled
  pair must carry RMG's choice; estimating the generated direction instead gives the
  two records RMG generates for one reaction (one per direction) two kinetics, which
  made one transition twice with weights decades apart (class overlaps, collisions).
* Two records of opposite generated direction can partition one region of transitions
  twice.  The partition contained in the other is the one covered, whichever side
  generated it; when neither contains the other the overlap stays refused.
* The two ends of an allyl radical are one orbit of the species' symmetry, related
  only across its resonance forms; the channel's orbit check compares each transition
  in the form it was matched in.
"""

import gzip
import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from rmgpy.kmc import compile as compile_module
from rmgpy.kmc import template_emission as te

REPOSITORY = Path(__file__).resolve().parents[3]
DECKS = REPOSITORY / "test/rmgpy/kmc/fixtures/decks"
FIXTURES = REPOSITORY / "test/rmgpy/kmc/fixtures/i249"
DATABASE = Path(os.environ.get("RMG_DATABASE_PATH", "/nonexistent"))
needs_database = pytest.mark.skipif(
    not (DATABASE / "input/kinetics/families").is_dir(),
    reason="RMG_DATABASE_PATH is not a database checkout",
)
FAMILIES = ("R_Recombination", "H_Abstraction", "R_Addition_MultipleBond",
            "intra_H_migration", "Disproportionation")
SITES = ("interior_radical:1", "pristine", "tail_end_radical+pendant_radical:0",
         "head_end_radical_short+monomer")


def _fixture(name):
    with gzip.open(FIXTURES / name, "rt") as handle:
        return json.load(handle)


# --- own-reverse families: the direction RMG keeps ----------------------------------


def _species(graphs, thermo):
    from rmgpy.molecule import Molecule
    from rmgpy.species import Species
    out = []
    for graph in graphs:
        species = Species(molecule=[Molecule().from_adjacency_list(graph)])
        species.generate_resonance_structures()
        species.thermo = thermo.get_thermo_data(species)
        out.append(species)
    return out


@needs_database
def test_own_reverse_records_carry_the_direction_rmg_keeps(tmp_path, monkeypatch):
    from rmgpy.data.rmg import get_db
    from rmgpy.rmg.model import CoreEdgeReactionModel

    original = compile_module._proxy_inputs

    def subset(request):
        proxies, source, options, discovery = original(request)
        return [p for p in proxies if p.site_type in SITES], source, options, discovery

    monkeypatch.setattr(compile_module, "_proxy_inputs", subset)
    path, _ = compile_module.compile_artifact(compile_module.CompileRequest(
        database_path=str(DATABASE), output_dir=str(tmp_path), families=FAMILIES,
        deck=str(DECKS / "pp_deck.py"), database_sha="own-reverse-direction"))
    records = json.loads(Path(path).read_text())["records"]
    by_id = {record["event_id"]: record for record in records}
    kinetics, thermo = get_db("kinetics"), get_db("thermo")
    model = SimpleNamespace(kinetics_estimator="rate rules", verbose_comments=True)
    kept = {True: 0, False: 0}
    wrong = []
    for record in records:
        family = kinetics.families[record["family"]]
        if not family.own_reverse or record["orientation"] != "as_generated" \
                or record["status"] != "enabled":
            continue
        generated = [reaction for reaction in kinetics.generate_reactions_from_families(
            _species(record["reactant_graphs"], thermo),
            products=_species(record["product_graphs"], thermo),
            only_families=[record["family"]]) if reaction.reverse is not None]
        if len(generated) != 1:
            continue
        reaction = generated[0]
        for species in reaction.reactants + reaction.products:
            species.thermo = thermo.get_thermo_data(species)
        # RMG's own decision, on the record's molecules
        _kinetics, _source, _entry, is_forward = \
            CoreEdgeReactionModel.generate_kinetics(model, reaction)
        twin = by_id[record["reverse_of"]]
        estimated, derived = (record, twin) if is_forward else (twin, record)
        template = ";".join(reaction.template if is_forward else reaction.reverse.template)
        kept[is_forward] += 1
        if (estimated["rate_source"]["kind"], derived["rate_source"]["kind"]) != \
                ("RMG family estimate", "reference-thermo reverse") or \
                f"[{template}]" not in estimated["rate_source"]["comment"] or \
                estimated["rate_source"]["family_template_direction"] != \
                ("forward" if is_forward else "reverse"):
            wrong.append((record["event_id"][:16], record["template"],
                          "RMG keeps " + ("forward" if is_forward else "reverse " + template),
                          estimated["rate_source"]["kind"]))
    # the reduced set holds own-reverse reactions RMG keeps in either direction
    assert kept[True] > 0 and kept[False] > 0, kept
    assert not wrong, wrong


# --- class overlaps: the contained partition is covered -------------------------------


def test_a_nested_overlap_covers_the_contained_partition():
    payloads = _fixture("nested_overlap.json.gz")["payloads"]
    merged, collisions = te.resolve_collisions(payloads)
    assert all(c["disposition"] == "merged" for c in collisions)
    out, table = te.resolve_overlaps(merged)
    (row,) = table
    # the as_generated templates make no transition the reversed ones do not, but not
    # the converse: the as_generated partition is the contained one
    assert row["log10_diameter"] < te.LOCALITY_TOLERANCE_LOG10
    assert row["disposition"] == "covered", row["reason"]
    orientation = {p["record_event_id"]: p["orientation"] for p in payloads}
    assert {orientation[c["record"]] for c in row["covered"]} == {"as_generated"}
    assert {orientation[c["record"]] for c in row["covered_by"]} == {"reversed"}
    loaded = [t for p in out for t in p["templates"] if t["status"] == "emitted"]
    assert loaded and all(orientation[p["record_event_id"]] == "reversed"
                          for p in out for t in p["templates"] if t["status"] == "emitted")


def test_an_overlap_nested_neither_way_loads_both_sides_once():
    """A record and its reversed twin whose windows each miss end placements the other
    covers: covering either partition would lose transitions, loading both as they are
    makes the shared ones twice.  At one weight both load, the reversed side deferring
    the shared transitions to the as_generated one, so each is made once and none is
    lost (I-271)."""
    from collections import Counter
    from rmgpy.kmc import window_matcher as wm
    payloads = _fixture("neither_nested_overlap.json.gz")["payloads"]
    merged, _ = te.resolve_collisions(payloads)
    out, table = te.resolve_overlaps(merged)
    (row,) = table
    assert row["log10_diameter"] < te.LOCALITY_TOLERANCE_LOG10
    assert row["disposition"] == "deduplicated", row.get("reason")
    assert te.overlap_gates(table) == []
    orientation = {p["record_event_id"]: p["orientation"] for p in out}
    for p in out:
        for t in p["templates"]:
            assert t["status"] == "emitted"
            defers = [(d["record"], d["key"]) for d in t.get("defers_to", ())]
            assert bool(defers) == (orientation[p["record_event_id"]] == "reversed")
            assert all(orientation[r] == "as_generated" for r, _k in defers)
    loaded = te.loaded_templates(out, 700.0)
    alone = [wm.Template(t.name, t.family, t.atoms, list(t.bonds), t.edits, weight=t.weight)
             for t in loaded]
    shared = 0
    for p in out:
        for t in p["templates"]:
            for wit in t["maps"]["witnesses"].values():
                for side in ("r", "p"):
                    state = te.melt_state(te.rows_from_json(wit[side]))
                    made = Counter(tr.identity[1:] for tr in te.loaded_transitions(loaded, state))
                    each = Counter(tr.identity[1:] for tr in te.loaded_transitions(alone, state))
                    assert set(made) == set(each) and set(made.values()) <= {1}
                    shared += sum(1 for n in each.values() if n > 1)
    assert shared


# --- one orbit across resonance forms -------------------------------------------------


def test_the_two_ends_of_an_allyl_radical_are_one_orbit():
    doc = _fixture("allyl_orbit_witness.json.gz")
    template = doc["template"]
    witness = template["maps"]["witnesses"][doc["label"]]
    r, p = te.rows_from_json(witness["r"]), te.rows_from_json(witness["p"])
    matcher = te.matcher_from_payload(dict(template, family=doc["family"]), 0.0,
                                      weight=1.0, name="L0")
    transitions = te.loaded_transitions([matcher], te.melt_state(r))
    # the witness's own channel: the hydrogen leaves the carbon that ends as the radical
    (donor,) = [a for a in p if p[a][3] == 1 and r[a][3] == 0 and p[a][0] == "C"]
    channel = [tr for tr in transitions
               if any(getattr(e, "atom", None) == donor and e.new == 1 for e in tr.edits)]
    ends = {a for tr in channel for a, value in tr.form_radicals if value == 1}
    assert len(channel) == 4 and len(ends) == 2      # two hydrogens to either allyl end
    assert te.channel_orbit_problem(r, channel) is None
    # different donors stay apart
    other = [tr for tr in transitions if tr not in channel]
    assert te.channel_orbit_problem(r, [channel[0]] + other) is not None
