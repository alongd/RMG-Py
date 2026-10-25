"""A reaction RMG rates from its training depository reaches the catalog
even where no rule template makes the same context (I-262).

Emission refuses a whole-species depository match: such a reaction rode only
as an override on a rule template that makes it, and assembly skipped it
without a row when none did.  On the PE toy that removed ethyl + ethyl
recombination and the n-butane cut, and every primary radical + ethylene
template in both directions (butyl to decyl + C2H4 are all training
entries).  The fix follows RMG:

* open context (artificial ports left): the member grows to the first size
  no depository entry matches (``windows.depository_start``), the rule
  template forms there and the sizes passed over are its override sites;
* closed context (every end physical, the whole species): the depository
  reaction is its own template (``windows.CLOSED_DEPOSITORY``);
* a depository reaction no loaded template carries is a named row of the
  override table, never a silent skip.
"""

import os
from pathlib import Path

import pytest

from templateEmissionTest import (DepositoryCounter, ToyDepository, _c13_match, _certificate,
                                  _depository_rates, _res, _toy_pair, toy_entry, toy_kc)

from rmgpy.kmc import boundaries as b
from rmgpy.kmc import emit
from rmgpy.kmc import run as kr
from rmgpy.kmc import template_emission as te
from rmgpy.kmc import windows as win

DECKS = Path(__file__).parent / "fixtures" / "decks"
DATABASE = Path(os.environ.get("RMG_DATABASE_PATH", "/nonexistent"))
needs_database = pytest.mark.skipif(not (DATABASE / "input/kinetics/families").is_dir(),
                                    reason="RMG_DATABASE_PATH is not a database checkout")


# -- assembly: no silent skip -------------------------------------------------


def _match(entry):
    return {"depository": "toy_cc/training", "index": entry.index, "label": entry.label,
            "orientation": "forward", "rank": 3, "kinetics": "None",
            "reaction": {side: [spc.molecule[0].to_adjacency_list() for spc in
                                getattr(entry.item, side)] for side in ("reactants", "products")}}


def test_a_depository_reaction_no_loaded_template_carries_is_a_named_row():
    """A certified match no loaded template makes (methyl + methyl: the C12
    pair's windows read two carbons per side) is a ``refused:no_template``
    row naming its source; a depository entry of the family that no payload
    names and no template makes is counted; the matched C13 still rides."""
    entry, match = _c13_match()
    methyl = toy_entry(11, ["[CH3]", "[CH3]"], ["CC"])
    other = toy_entry(12, ["[CH3]", "[H]"], ["C"])
    counter = DepositoryCounter([ToyDepository([entry, methyl, other])])
    payloads, counter = _toy_pair(counter=counter, certificate=_certificate(match, _match(methyl)))
    resolved = _res(payloads, counter)
    out, table, failures = te.depository_overrides(
        resolved, counter=counter, rates=_depository_rates(counter, None), kc=toy_kc)
    assert failures == [], failures
    (named,) = [r for r in table if r["status"] == "refused:no_template"]
    assert {s.get("index") for s in named["sources"]} == {11}
    assert all(not s.get("enumerated") for s in named["sources"])
    (counted,) = [r for r in table if r["status"] == "not_made:enumerated"]
    assert counted["reactions"] == 1
    assert sorted(r["status"] for r in table if r.get("direction")) == ["ok", "ok"]
    assert all(t.get("overrides") for p in out for t in p["templates"] if t["base"])


# -- windows on the PE toy declaration ------------------------------------------


@pytest.fixture(scope="module")
def pe():
    from rmgpy.kmc import channels as chm
    spec = kr.read_deck(DECKS / "pe_toy_deck.py").spec
    systems, refusals = emit.build_reference_systems(emit.declared_reference_systems(spec, 3))
    assert not refusals
    db = emit.load_database(["R_Addition_MultipleBond", "R_Recombination", "intra_H_migration"],
                            DATABASE)
    fam = db.kinetics.families
    declared = emit.declared_provider(spec)
    policy = win.declared_interim(declared.declared_elements())
    certifier = emit.make_certifier(fam, {"domain": win.certificate_domain(policy),
                                          "repeat": declared.unit_formula(),
                                          "growth": declared}, emit.reference_thermo(db))

    def case(system, family, template_side, direction=None, other_side=None):
        cases, _stats = emit.reference_cases(certifier, {family: fam[family]},
                                             {system: systems[system]}, chm.SpeciesKeys())
        for found in cases.values():
            for c, _specs in found:
                mols = win.template_molecules(c.witness, c.direction)
                other = sorted(m.to_smiles() for m in te.molecules_of(dict(c.witness.product_rows())))
                if sorted(m.to_smiles() for m in mols) == sorted(template_side) and \
                        direction in (None, c.direction) and other_side in (None, other):
                    return c
        raise AssertionError(f"no {family} member on {system} with {template_side}")
    return fam, certifier, case


def _run(pe, family, case, **kw):
    fam, certifier, _case = pe
    keep = []
    row = win.run_member(certifier, fam[family], case,
                         bound=win.span_bound(win.family_constants(fam[family])), keep=keep, **kw)
    return row, keep[0]


@needs_database
def test_an_open_member_grows_past_the_training_sizes_to_its_rule_template(pe):
    """Butyl + C2H4: pentyl, hexyl, heptyl and octyl + C2H4 are training
    entries too, so the member grows one period a step to nonyl, the first
    size RMG rates by rule; the rule template formed there makes butyl +
    C2H4, every size passed over is named, and the port budget counts from
    nonyl."""
    fam, _certifier, case = pe
    family = "R_Addition_MultipleBond"
    member = case("pair:head_end_radical_short+monomer:left", family, ["C=C", "[CH2]CCC"],
                  "forward")
    row, fp = _run(pe, family, member, contexts=False)
    assert row["status"] == "fixed_point", (row["status"], row["reason"])
    assert row["depository"] == {"matched": 5, "closed": False}
    assert [w.heavy_count() for w in fp.matched] == [6, 7, 8, 9, 10]
    assert all(win.depository_match(fam[family], w) is not None for w in fp.matched)
    start = fp.matched[-1]
    for p in b.vector_ports(start):
        start = b.extend(start, p)[0]
    assert win.depository_match(fam[family], start) is None
    cut = te.cut_template(fp.witness, fp.window.atoms, family=family, pathway=member.pathway,
                          direction=member.direction)
    r, p = te.template_sides(member.witness)
    count, _inside, _ch = te.count_on([cut.matcher(1.0, "t")], r, p, family=fam[family],
                                      direction=member.direction, pathway=tuple(member.pathway),
                                      counter=te.RmgCounter(fam))
    assert count.matcher and count.equal and count.target_found


@needs_database
def test_ethyl_plus_ethyl_is_a_closed_context_of_its_own(pe):
    """n-Butane -> ethyl + ethyl: RMG rates it from training (C2H5 + C2H5 <=>
    C4H10); the member grows to propyl + propyl for its rule window, and the
    context whose two ends are both physical is the whole species: a closed
    fixed point whose window is every atom."""
    fam, _certifier, case = pe
    family = "R_Recombination"
    member = case("single:pristine~4", family, ["C[CH2]", "C[CH2]"])
    out = []
    row, fp = _run(pe, family, member, contexts=True, contexts_out=out)
    assert row["status"] == "fixed_point", (row["status"], row["reason"])
    assert row["depository"] == {"matched": 1, "closed": False}
    (results, _merges), = out
    closed = [r for r in results if r.fixed_point is not None and r.fixed_point.closed]
    assert len(closed) == 1
    (c,) = closed
    assert c.reason == win.CLOSED_DEPOSITORY
    assert sorted(e for _p, e in c.assignment) == ["declared", "declared"]
    assert c.window.atoms == frozenset(dict(c.witness.reactant_rows()))
    assert win.depository_match(fam[family], c.witness)["label"] == "C2H5 + C2H5 <=> C4H10"


@needs_database
def test_a_reacting_atom_on_a_port_keeps_only_its_closed_context(pe):
    """Methyl + methyl -> ethane: the reacting carbons are the ports, so no
    growth exists and the open base is refused as the cut-point event it is;
    the context whose two ends are declared is ethane whole, the depository
    reaction's own template, and the member is a fixed point on it alone."""
    fam, _certifier, case = pe
    family = "R_Recombination"
    member = case("single:pristine~2", family, ["[CH3]", "[CH3]"])
    row, _fp = _run(pe, family, member, contexts=False, covers=False)
    assert row["status"] == "refused:reacting_port", (row["status"], row["reason"])
    assert "depository size 1" in row["reason"]
    out = []
    row, fp = _run(pe, family, member, contexts=True, covers=False, contexts_out=out)
    assert row["status"] == "fixed_point", (row["status"], row["reason"])
    assert fp.status == "refused:reacting_port"
    (results, merges), = out
    assert merges == () and results[0].fixed_point is fp
    (c,) = results[1:]
    assert c.fixed_point.closed and c.window.atoms == frozenset(dict(c.witness.reactant_rows()))
    assert win.depository_match(fam[family], c.witness)["label"] == "CH3 + CH3 <=> C2H6"


@needs_database
def test_a_reacting_atom_a_declared_end_would_delete_stays_a_named_refusal(pe):
    """H + methyl -> methane: the reacting H stands on a port, and a declared
    end there deletes it, so no context closes: the member stays refused by
    its open growth, every closed context tried is named."""
    _fam, _certifier, case = pe
    family = "R_Recombination"
    member = case("single:pristine~1", family, ["[CH3]", "[H]"])
    row, _fp = _run(pe, family, member, contexts=True, covers=False)
    assert row["status"] == "refused:reacting_port", (row["status"], row["reason"])
    assert row["closed_contexts"]
    assert all(r["status"] == "refused:end" and not r["closed"] for r in row["closed_contexts"])


@needs_database
def test_an_own_reverse_family_at_a_depository_size_is_a_named_refusal(pe):
    """n-Hexyl -> 2-hexyl (intra_H_migration, one template both ways): RMG
    rates it from training and no depository override exists for an
    own-reverse family, so the member is refused by name and keeps no closed
    context: the rule template of a longer member makes n-hexyl's shift too
    (the assembly overlap gate refuses a second template on it)."""
    fam, _certifier, case = pe
    family = "intra_H_migration"
    member = case("single:head_end_radical", family, ["[CH2]CCCCC"], other_side=["C[CH]CCCC"])
    assert win.depository_match(fam[family], member.witness) is not None
    out = []
    row, fp = _run(pe, family, member, contexts=True, covers=False, contexts_out=out)
    assert fp.status == row["status"] == "refused:own_reverse_depository", (row["status"], row["reason"])
    assert fp.matched == () and "closed_contexts" not in row
    assert all(r.fixed_point is None or not r.fixed_point.closed for results, _m in out for r in results)