"""The kernel table counts a transition two deduplicated templates share once (I-271).

A template that defers to another (``defers_to``) makes only the transitions
that template does not.  The kernel table reads the melt engine's index: a
unimolecular class's transitions, and a bimolecular class's sites per
participant, whose kernel term is a product of site counts.  The PP overlaps
the catalog deduplicates have no kernel term on the declared states (their
templates are kMC-only there), so the case is a PE toy template and a twin of
it deferring to it: every transition of the twin is the template's, and the
twin adds nothing to the table.
"""

import copy
import hashlib

import pytest

import melt_engine_support as S
from meltRunTest import PE_LIVE, toy_payloads

from rmgpy.kmc import kernel_rates as kr_
from rmgpy.polymer_spec import LinearHomopolymerSpec

T = 700.0


def _table(payloads):
    spec = LinearHomopolymerSpec.from_mapping(S.pe_mapping({10: 1}))
    return kr_.kernel_rates({"schema": "kmc-template-catalog/1", "payloads": payloads}, spec, T,
                            allow_unbound=True)


def _with_twin(payloads, name):
    """``payloads`` and a twin of template ``name`` deferring to it."""
    (original,) = [p for p in payloads if p["record_event_id"] == f"evt_{name}"]
    twin = copy.deepcopy(original)
    twin["record_event_id"] = f"evt_twin_{name}"
    (t,) = twin["templates"]
    t["key"] = hashlib.sha256(f"twin:{name}".encode()).hexdigest()
    t["defers_to"] = [{"record": original["record_event_id"], "key": original["templates"][0]["key"]}]
    return payloads + [twin], twin["record_event_id"][4:16] + ":" + t["key"][:12]


@pytest.mark.parametrize("name", ["h_abstraction", "homolysis"])
def test_a_template_deferring_every_transition_adds_nothing_to_the_kernel_table(name):
    payloads = toy_payloads([t for t in S.pe_templates() if t.name in PE_LIVE])
    base = _table(payloads)
    both, twin = _with_twin(payloads, name)
    table = _table(both)
    assert twin in [r["template"] for r in table["templates"]]
    # the twin is in no kernel: each of its transitions is its original's
    placed = {x["template"] for k in table["kernels"].values() for x in k["templates"]}
    assert twin not in placed
    assert twin in {o["template"] for o in table["kmc_only"]}
    # and every other template's row and term are the table without it
    assert [r for r in table["templates"] if r["template"] != twin] == base["templates"]
    assert {kid: [x for x in k["templates"] if x["template"] != twin] for kid, k in table["kernels"].items()} \
        == {kid: k["templates"] for kid, k in base["kernels"].items()}
