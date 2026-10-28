"""Red tests (I-279b): an interior law on a chain type that is its own reverse
is read off the piece's mirror, and a law the state model cannot represent is
refused by name instead of aborting the table (``rmgpy.kmc.kernel_rates``).

The fixture is the final PP catalog's disproportionation of two mid-chain
radicals (payload ``159646410c91`` with its reciprocal twin), bound to the PP
deck.  On a c1~c1 or c2~c2 piece (odd backbone, its own reverse) positions are
folded to the nearer end and an internal alkene is read from its nearer end.
A CH donor gives up a hydrogen of either backbone neighbour, so far from the
ends it makes one alkene read C=CH and one read CH=C.  At the centre atom of
the piece the two alkenes are mirror images: both are read from their own
nearer end, so both read C=CH.  The full PP build read that centre state as a
far-field state and stopped the whole table at law 103/1573.

Every expectation is read independently of the table builder with the
state-model tests' brute force (``kernelRatesStateModelTest.brute``).
"""

import gzip
import json
from collections import Counter
from pathlib import Path

import numpy
import pytest

import melt_engine_support as S
from kernelRatesStateModelTest import brute
from meltRunTest import toy_payloads, write_catalog

from rmgpy.kmc import kernel_rates as kr_
from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import run as kr
from rmgpy.polymer_spec import LinearHomopolymerSpec

T = 700.0
FIXTURES = Path(__file__).parent / "fixtures"
DISPROPORTIONATION = "159646410c91:82b4b3b95a92"    # tertiary acceptor takes a backbone CH2 H of a CH donor
DONOR = "radical_interior:CH"
ACCEPTOR = "radical_interior:C"


@pytest.fixture(scope="module")
def pp():
    doc = json.loads(gzip.decompress((FIXTURES / "i279" / "pp_disproportionation_catalog.json.gz").read_bytes()))
    spec = kr.read_deck(FIXTURES / "decks" / "pp_deck.py").spec
    table = kr_.kernel_rates(doc, spec, T)
    n_hi = table["laws_exact_to_n"]
    engine = me.MeltEngine(kr_._single_chain_spec(spec, n_hi), me.PayloadCatalog(doc["payloads"]), T,
                           numpy.random.default_rng(0))
    probes = kr_._Probes(spec, engine, 1, 2, 2 * n_hi + 3)
    return {"table": table, "engine": engine, "probes": probes, "spec": spec, "n_hi": n_hi,
            "types": table.get("chain_types")}


def donor_terms(table):
    """{donor class: term} of the disproportionation's terms with the declared type's tertiary acceptor."""
    out = {}
    for k in table["kernels"].values():
        for t in k["templates"]:
            if t["template"] == DISPROPORTIONATION and k["classes"][0] == ACCEPTOR:
                out[k["classes"][1]] = t
    return out


def donor_of(pp, label):
    """(class index, engine participant) of the disproportionation's donor of class ``label``."""
    (ci,) = [i for i, t in enumerate(pp["engine"].templates) if t.name == DISPROPORTIONATION]
    (row,) = [r for r in pp["table"]["templates"] if r["template"] == DISPROPORTIONATION]
    (p,) = [q for q, classes in enumerate(row["participants"]) if label in classes]
    return ci, p


def partner(pp):
    """A tertiary radical far from the ends of a declared chain: the acceptor of every brute-force firing."""
    n = pp["n_hi"]
    probes = pp["probes"]
    xs = [x for x in probes.positions("radical_interior", n) if probes.carries(ACCEPTOR, n, x)]
    return ("radical_interior", n, xs[len(xs) // 2], None)


def labels(counter):
    return Counter({tuple(sorted(p[0] for p in prods)): m for prods, m in counter.items()})


# --------------------------------------------------------------------------
# the mirror of an own-reverse piece is not a far-field state
# --------------------------------------------------------------------------
def test_the_centre_of_an_own_reverse_piece_folds_both_alkenes_onto_one_reading(pp):
    """The engine itself (no table): a c1~c1 donor at the centre of 33 backbone atoms makes two C=CH readings,
    the donor two atoms off the centre one C=CH and one CH=C."""
    ci, p = donor_of(pp, f"{DONOR}/c1~c1")
    centre = brute(pp, ci, p, "radical_interior", 16, 16, tau="c1~c1", partner=partner(pp))
    off = brute(pp, ci, p, "radical_interior", 16, 14, tau="c1~c1", partner=partner(pp))
    assert labels(centre) == Counter({("alkene_internal:C=CH/c1~c1",): 2})
    assert labels(off) == Counter({("alkene_internal:C=CH/c1~c1",): 1, ("alkene_internal:CH=C/c1~c1",): 1})


def test_the_donor_has_a_law_on_every_chain_type(pp):
    got = donor_terms(pp["table"])
    assert sorted(got) == [DONOR, f"{DONOR}/c1~c1", f"{DONOR}/c2~c2"]
    for label, t in got.items():
        assert t["law"] is not None, (label, t.get("law_refused"))


@pytest.mark.parametrize("tau", [None, "c1~c1", "c2~c2"])
def test_the_donor_far_law_equals_a_brute_force_enumeration_off_the_mirror(pp, tau):
    """Every donor position on two lengths past the probes, at least a quarter of the piece from each end and,
    on an own-reverse piece, two atoms or more off its centre."""
    label = DONOR if tau is None else f"{DONOR}/{tau}"
    law = donor_terms(pp["table"])[label]["law"]["S1"]
    assert law["class"] == label
    ci, p = donor_of(pp, label)
    probes = pp["probes"]
    checked = 0
    for n in (pp["n_hi"] + 1, pp["n_hi"] + 2):
        L = 2 * n + (0 if tau is None else 1)
        for x in probes.positions("radical_interior", n, tau=tau):
            if not probes.carries(label, n, x) or min(x, L - 1 - x) < L // 4:
                continue
            if tau is not None and L - 1 - 2 * x < 4:
                continue
            want = brute(pp, ci, p, "radical_interior", n, x, tau=tau, partner=partner(pp))
            got = Counter({tuple(sorted(k, key=repr)): v
                           for k, v in kr_.predict_near(law, None, n, 2, x, types=pp["types"]).items()})
            assert got == want, (label, n, x)
            checked += 1
    assert checked


def test_every_term_has_a_law_or_says_why(pp):
    for k in pp["table"]["kernels"].values():
        for t in k["templates"]:
            assert t["law"] is not None or t.get("law_refused") or t.get("law_inert"), t


# --------------------------------------------------------------------------
# a law the state model cannot represent is refused by name
# --------------------------------------------------------------------------
def test_a_law_the_model_cannot_represent_is_refused_by_name_not_raised(monkeypatch, tmp_path):
    """The PE toy with every relative product code made to depend on its anchor's parity: an interior law's
    neighbouring states then differ and a chain law's end effects reach the middle of the probe.  Each such
    term is refused with its reason; the table is still built."""
    original = kr_._relative

    def parity(product, L, k, nbb):
        token, alpha, beta, gamma = original(product, L, k, nbb)
        return token, alpha, beta, gamma + nbb * (k % 2)

    monkeypatch.setattr(kr_, "_relative", parity)
    path = write_catalog(tmp_path / "catalog.json", toy_payloads(S.pe_templates()), tag="mirror-pe",
                         bind=S.pe_mapping({1: 1}))
    table = kr_.kernel_rates(json.loads(path.read_text()), LinearHomopolymerSpec.from_mapping(S.pe_mapping({10: 1})),
                             T)
    got = {kid: t for kid, k in table["kernels"].items() for t in k["templates"]}
    beta = got["beta_scission:substrate=radical_interior:CH"]
    assert beta["law"] is None and "differ between neighbouring positions" in beta["law_refused"]
    homolysis = got["homolysis:substrate=chain"]
    assert homolysis["law"] is None and homolysis["law_refused"]
    for kid, t in got.items():
        assert t["law"] is not None or t.get("law_refused") or t.get("law_inert"), kid
