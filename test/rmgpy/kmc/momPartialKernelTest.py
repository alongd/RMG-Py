"""Red tests (I-282): a kernel the table gives a law on some templates and a
named refusal on others maps its lawful templates and lists the others.

A kernel's rate is the sum over its templates, and the table gives each
template its own law or its own refusal (``law_refused`` / ``law_inert``).
Homolysis and end recombination refused the whole kernel on one template
without a law; on the full polypropylene table that removed every radical
source (each homolysis kernel carries a template whose product is read as
another orientation).  Now the kernel carries its lawful templates exactly as
it would carry them alone, every omitted template is listed under
``record.inert`` with the table's reason, and ``record.partial_kernels``
names and counts what each partly mapped kernel omits.

The refused templates are copies of the kernel's own lawful template with
the law removed and a refusal reason: the lawful mapping must not move.
"""

import copy
import gzip
import json
from pathlib import Path

import pytest

from momKernelProducerTest import C, H, RADICAL, table, terms
from momProducerTest import BASE, LABEL, spec_of

from rmgpy.kmc import kernel_rates as kr_
from rmgpy.kmc import mom_producer as mp
from rmgpy.kmc import run as kr

REFUSED = "ffffffffffff:000000000001"
REFUSAL = "a product is not the declared radical_end of 1 periods"
FIXTURES = Path(__file__).parent / "fixtures"
PP_DECK = FIXTURES / "decks" / "pp_deck.py"
PP_LABELS = dict(label="PP", radical_label="PP_re", alkene_label="PP_ut", interior_radical_label="PP_rm",
                 interior_alkene_label="PP_ui")


def with_refused_copy(doc, kernel, name=REFUSED, reason=REFUSAL):
    """The table with a second template on ``kernel``: its first template's row and term, the law removed."""
    doc = copy.deepcopy(doc)
    term = copy.deepcopy(doc["kernels"][kernel]["templates"][0])
    row = copy.deepcopy(next(r for r in doc["templates"] if r["template"] == term["template"]))
    term.update(template=name, law=None, law_refused=reason)
    term.pop("law_inert", None)
    row["template"] = name
    doc["kernels"][kernel]["templates"].append(term)
    doc["templates"].append(row)
    return doc


def listed(out, kernel, name):
    return [x for x in out["record"]["inert"] if x["kernel"] == kernel and name in x["templates"]]


# --------------------------------------------------------------------------
# kmc-kernel-rates/3: homolysis and end recombination
# --------------------------------------------------------------------------
@pytest.mark.parametrize("kernel, key", [(H, "homolysis"), (C, "recombination")])
def test_a_refused_template_beside_a_lawful_one_leaves_the_lawful_mapping(kernel, key):
    alone = terms(table())
    out = terms(with_refused_copy(table(), kernel))
    assert out[key] == alone[key]
    assert out["record"]["laws"] == alone["record"]["laws"]
    (entry,) = listed(out, kernel, REFUSED)
    assert entry["reason"] == REFUSAL


@pytest.mark.parametrize("kernel", [H, C])
def test_a_partly_mapped_kernel_names_and_counts_what_it_omits(kernel):
    out = terms(with_refused_copy(table(), kernel))
    (part,) = out["record"]["partial_kernels"]
    lawful = table()["kernels"][kernel]["templates"][0]["template"]
    assert part == {"kernel": kernel, "templates": 2, "mapped": [lawful],
                    "omitted": [{"template": REFUSED, "reason": REFUSAL}]}


def test_a_kernel_mapped_whole_is_not_a_partial_kernel():
    assert terms(table())["record"]["partial_kernels"] == []


def test_a_homolysis_kernel_without_a_lawful_template_maps_nothing_and_is_listed():
    doc = table()
    (term,) = doc["kernels"][H]["templates"]
    term.update(law=None, law_refused=REFUSAL)
    out = mp.kernel_input(doc, spec_of(BASE), [H, C], temperature=700.0, phi=1.0, label=LABEL,
                          radical_label=RADICAL, volume_ratio=1.0)
    assert out["homolysis"] == []
    (entry,) = listed(out, H, term["template"])
    assert entry["reason"] == REFUSAL
    assert out["record"]["partial_kernels"] == []


# --------------------------------------------------------------------------
# kmc-kernel-rates/4 (chain types): homolysis per chain pool
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def pp_table():
    doc = json.loads(gzip.decompress((FIXTURES / "i279" / "pp_state_model_catalog.json.gz").read_bytes()))
    return kr_.kernel_rates(doc, kr.read_deck(PP_DECK).spec, 700.0)


def pp_terms(doc):
    return mp.kernel_input(doc, kr.read_deck(PP_DECK).spec, sorted(doc["kernels"]), temperature=700.0, phi=1.0,
                           volume_ratio=1.0, far_law=True, positions=True, **PP_LABELS)


def test_a_typed_homolysis_kernel_maps_its_lawful_templates_per_chain_pool(pp_table):
    kernel = "homolysis:substrate=chain"
    alone = pp_terms(pp_table)
    out = pp_terms(with_refused_copy(pp_table, kernel))
    assert out["homolysis_by_pool"] == alone["homolysis_by_pool"]
    (entry,) = listed(out, kernel, REFUSED)
    assert entry["reason"] == REFUSAL
    # the fixture's own partly mapped kernels (a backbite template without a far-field site) are unchanged
    (part,) = [p for p in out["record"]["partial_kernels"] if p["kernel"] == kernel]
    assert part["templates"] == 2
    assert part["omitted"] == [{"template": REFUSED, "reason": REFUSAL}]
    assert [p for p in out["record"]["partial_kernels"] if p["kernel"] != kernel] == \
        alone["record"]["partial_kernels"]


# --------------------------------------------------------------------------
# kmc-kernel-rates/4: a join template at single sizes beside the far-field join
# --------------------------------------------------------------------------
JOIN = "recombination:radical=radical_end:CH+radical=radical_end:CH2"
SINGLE = json.loads((FIXTURES / "i282" / "pp_single_size_join_template.json").read_text())


def with_single_size_join(doc):
    """The fixture table with the full table's single-size template on the end recombination."""
    doc = copy.deepcopy(doc)
    doc["kernels"][JOIN]["templates"].append(copy.deepcopy(SINGLE["term"]))
    doc["templates"].append(dict(copy.deepcopy(SINGLE["row"]), kernels=[JOIN]))
    return doc


def test_a_single_size_join_template_is_a_size_departure_beside_the_far_field_join(pp_table):
    alone = pp_terms(pp_table)
    out = pp_terms(with_single_size_join(pp_table))
    assert out["joins"] == alone["joins"]
    name = SINGLE["term"]["template"]
    (dep,) = [d for d in out["record"]["size_departures"] if d["kernel"] == JOIN]
    assert dep["template"] == name
    assert "single sizes" in dep["reason"]
    (part,) = [p for p in out["record"]["partial_kernels"] if p["kernel"] == JOIN]
    assert part["mapped"] == [t["template"] for t in pp_table["kernels"][JOIN]["templates"]]
    assert part["omitted"] == [{"template": name, "reason": dep["reason"]}]


# --------------------------------------------------------------------------
# kmc-kernel-rates/4: a law the table refuses at one reachable length and position
# --------------------------------------------------------------------------
MOVE = "intra_h_transfer:substrate=radical_interior:CH"
SHORT = "1cfafb5c5066:e34cd8d621f9"            # its end-distance law is exact below 6 periods
OFF = "a product with more than one feature, or a feature off the backbone"


def test_a_law_refused_at_a_reachable_length_lists_its_template_unmapped(pp_table):
    """Homolysis makes radicals of 2 periods, so a move's end-distance law is read at 5 periods, where the
    full table refuses some entries: the template is listed with the table's reason, the kernel's others map."""
    doc = copy.deepcopy(pp_table)
    (term,) = [t for t in doc["kernels"][MOVE]["templates"] if t["template"] == SHORT]
    (x,) = term["law_end_distance"]["exact"]["5"]
    term["law_end_distance"]["exact"]["5"][x] = {"refused": OFF}
    out = pp_terms(doc)
    (entry,) = [u for u in out["record"]["unmapped"] if u["kernel"] == MOVE and SHORT in u["templates"]]
    assert OFF in entry["reason"] and "5 periods" in entry["reason"]
    # the template also sits on the /c2~c2 kernel: this kernel's laws only
    labels = [law["template"] for law in out["record"]["laws"] + out["record"]["laws_outside_attribution"]
              if law["kernel"] == MOVE]
    assert not any(SHORT in label for label in labels)
    assert any(label.startswith("1cfafb5c5066:704f98848600") for label in labels)
    (part,) = [p for p in out["record"]["partial_kernels"] if p["kernel"] == MOVE]
    assert SHORT in [o["template"] for o in part["omitted"]]
