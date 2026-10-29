"""Red tests (I-283): the monomer read from either end is one species, and a
propagation template's sites on it are counted once.

``kernel_rates`` reads every cut of a declared chain from both sides, so on a
/4 table an asymmetric monomer (propene) is a state of two terminal-alkene
classes, one per end (``alkene_terminal:CH2=CH`` and ``CH=CH2``).  Each
reading counts the same sites of the same molecule, and the table has one
addition kernel per reading, both carrying the same template.  The melt engine
(kMC) fires a template over the sites it indexes on the molecule, which are
the sites of one reading, not of both.  The producer wrote both kernels' terms
onto the radical pool, which MoM refuses ("kernel_propagation has a duplicate
template") and which, summed, would count every addition twice.  Now the
template's addition is one term, counted once, and ``record.propagation_shared``
names the kernels that state it and the one whose term carries it.

The two kernels are the full polypropylene table's monomer additions onto
``radical_end:CH`` (``fixtures/i283/pp_monomer_addition_kernels.json``), set into
the I-279 fixture table.
"""

import copy
import json
from pathlib import Path

import pytest

from momPartialKernelTest import pp_table, pp_terms   # noqa: F401  (pp_table is a fixture)

from rmgpy.kmc import mom_producer as mp

FIXTURES = Path(__file__).parent / "fixtures"
ADDITIONS = json.loads((FIXTURES / "i283" / "pp_monomer_addition_kernels.json").read_text())
A = "addition:alkene=alkene_terminal:CH2=CH+radical=radical_end:CH"
B = "addition:alkene=alkene_terminal:CH=CH2+radical=radical_end:CH"
SHARED = "a85416d1286c:6a4b867414f0"


def with_additions(doc, kernels=(A, B)):
    """The fixture table with the monomer-addition kernels ``kernels`` and their templates' rows."""
    doc = copy.deepcopy(doc)
    for kid in kernels:
        doc["kernels"][kid] = copy.deepcopy(ADDITIONS["kernels"][kid])
    for row in ADDITIONS["rows"]:
        doc["templates"].append(dict(copy.deepcopy(row), kernels=[k for k in row["kernels"] if k in kernels]))
    return doc


def propagation(out):
    (pool,) = out["propagation"]
    return pool, out["propagation"][pool]


@pytest.mark.parametrize("kernel", [A, B])
def test_one_reading_alone_is_one_propagation_term(pp_table, kernel):
    pool, terms = propagation(pp_terms(with_additions(pp_table, [kernel])))
    (term,) = terms
    assert term["template"] == f"{SHARED}@{pool}"


def test_both_readings_of_the_monomer_give_the_one_term_of_either_alone(pp_table):
    """MoM refuses a pool whose kernel_propagation repeats a template; summed, the two would count it twice."""
    alone = propagation(pp_terms(with_additions(pp_table, [A])))
    both = propagation(pp_terms(with_additions(pp_table)))
    assert both == alone
    templates = [t["template"] for t in both[1]]
    assert len(templates) == len(set(templates))


def test_the_record_names_the_kernels_that_state_the_shared_template(pp_table):
    out = pp_terms(with_additions(pp_table))
    pool, (term,) = propagation(out)
    (shared,) = out["record"]["propagation_shared"]
    assert shared["pool"] == pool and shared["template"] == SHARED
    assert shared["kernels"] == [A, B] and shared["counted_in"] == A
    assert shared["monomer_classes"] == ["alkene_terminal:CH2=CH", "alkene_terminal:CH=CH2"]
    assert shared["k_pair"] == term["k_pair"]
    assert "once" in shared["reason"]


def test_a_kernel_whose_template_another_reading_carries_is_not_omitting_it(pp_table):
    out = pp_terms(with_additions(pp_table))
    refused = "dc87561126f2:29d83aca351f"
    assert [p for p in out["record"]["partial_kernels"] if p["kernel"] in (A, B)] == [
        {"kernel": kid, "templates": 2, "mapped": [SHARED],
         "omitted": [{"template": refused, "reason": next(t["law_refused"] for t in ADDITIONS["kernels"][kid]
                                                          ["templates"] if t["template"] == refused)}]}
        for kid in (A, B)]


def test_readings_that_count_the_template_differently_are_refused(pp_table):
    """Two readings of one molecule count the same sites; a difference means they are not one species."""
    doc = with_additions(pp_table)
    (term,) = [t for t in doc["kernels"][B]["templates"] if t["template"] == SHARED]
    term["S0"]["profile"] = {k: [2] for k in term["S0"]["profile"]}
    term["law"]["S0"]["exact"] = {n: [[2 * m, o] for m, o in x] for n, x in term["law"]["S0"]["exact"].items()}
    with pytest.raises(mp.ProducerError, match="counts its sites on the monomer"):
        pp_terms(doc)


def test_a_table_without_shared_templates_records_none(pp_table):
    assert pp_terms(pp_table)["record"]["propagation_shared"] == []
