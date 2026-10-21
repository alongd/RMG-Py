"""Cross-method rate identity (I-230): the melt engine's per-event propensity
on explicit states, against the MoM terms the producer writes for the same
kernels, so that polymer-mom's own test (its
``test/rmgpy/solver/polymerCrossMethodTest.py``, a separate checkout) can
check its rates against the kMC code rather than against sums that share its
normalization.

The chemistry is the in-tree PE catalog (``i213_pe_catalog.json.zlib.b64``,
emitted templates with thermo-derived reverses) on its deck, at 700 K, in a
melt volume of 1e-24 m3 (V_poly the same: volume ratio 1).  Every chain is at
least the terms' validity length (60 periods).  The document
``cross_method_document`` builds is the fixture
``fixtures/i230_cross_method_propensities.json``; the MoM checkout carries a
byte-identical copy (``CROSS_METHOD_SHA256`` pins both).

Per case: the explicit chains (pool -> {periods: count}), the kernel's MoM
term, the engine's events per second over the term's templates (well-mixed:
a MET pair's intrinsic k_allowed, before the cage factor), and, for a kernel
whose two participants come from one pool, the engine's site pairs S0 S1 and
its same-molecule pairs X (excluded from kMC events, present in the mean
field): MoM's events are the engine's times S0 S1 / (S0 S1 - X).

Detailed balance: the D1 term's reverse (the thermo-derived reverse of its
template) carries its record's thermo; Kc(700 K) of R + R -> P + U from it,
at the declared reference state C0 = P0 / R T (``ThermoData-Kc/1``), must be
the MoM forward coefficient over the reverse one.
"""

import base64
import hashlib
import json
import math
import zlib
from pathlib import Path

import numpy as np
import pytest

import melt_engine_support as S

from rmgpy import constants
from rmgpy.kmc import kernel_rates as K
from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import mom_producer as mp
from rmgpy.kmc import run as kr
from rmgpy.thermo import ThermoData

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE = FIXTURES / "i230_cross_method_propensities.json"
CROSS_METHOD_SHA256 = "8d07044b69583c68a08968d78cc502c98e5811016c1c4c3b6f5dbf2037c01255"
T = 700.0
VOLUME = 1.0e-24
C = "recombination:radical=radical_end:CH2+radical=radical_end:CH2"
D1 = "disproportionation:acceptor=radical_end:CH2+donor=radical_end:CH2"
DR = "reverse_disproportionation:alkene=alkene_terminal:CH2=CH+donor=chain"
T1 = "h_transfer:acceptor=radical_end:CH2+donor=chain"
LABELS = dict(label="P", radical_label="Re", alkene_label="Ut", interior_radical_label="Rm",
              interior_alkene_label="Ui")


def radical(n):
    return "[CH2]" + "C" * (n - 1)


def alkane(n):
    return "C" * n


def alkene(n):
    return "C=C" + "C" * (n - 2)


# case: (name, kernel, MoM pool owning the term, the term's products / partner, explicit molecules by pool)
CASES = (
    ("recombination_like", C, "Re", {"Re": [radical(60), radical(60)]}),
    ("recombination_unlike", C, "Re", {"Re": [radical(60), radical(66)]}),
    ("d1_like", D1, "Re", {"Re": [radical(60), radical(60)]}),
    ("d1_unlike", D1, "Re", {"Re": [radical(60), radical(66)]}),
    ("reverse_d1", DR, "Ut", {"Ut": [alkene(64)], "P": [alkane(70)]}),
    ("h_transfer", T1, "Re", {"Re": [radical(60)], "P": [alkane(70)]}),
)


def pe():
    text = (FIXTURES / "i213_pe_catalog.json.zlib.b64").read_bytes().replace(b"\n", b"")
    doc = json.loads(zlib.decompress(base64.b64decode(text)))
    deck = kr.read_deck(FIXTURES / "i213_pe_toy_deck.py")
    return doc, deck


def per_template(engine):
    """Events per second of each loaded template on the engine's state: a MET
    pair's intrinsic rate (k_allowed / N_A V, no cage factor) and every
    bimolecular channel's propensity."""
    names = [c.name for c in engine.classes]
    out = {}
    for pr in engine.met_pair_rates().values():
        for row in pr.classes:
            out[names[row[0]]] = out.get(names[row[0]], 0.0) + engine.weights[row[0]] * row[1] / engine.nav
    for key, a in engine.propensities().channels:
        if key[0] == "bi":
            out[names[key[1]]] = out.get(names[key[1]], 0.0) + a
    return out


def site_pairs(engine, templates):
    """(S0 S1, X) of the templates: all ordered site pairs of the state and
    those inside one molecule, from the engine's index entries."""
    names = [c.name for c in engine.classes]
    s0 = s1 = x = 0
    for entry in engine.entries.values():
        n = [0, 0]
        for (ci, p, _kind), sites in entry.shape.sites.items():
            if names[ci] in templates:
                n[p] += len(sites)
        s0, s1, x = s0 + n[0], s1 + n[1], x + n[0] * n[1]
    return s0 * s1, x


def kc(derived, temperature):
    """Kc of the forward reaction (reactants -> products of the record's
    thermo) at C0 = P0 / R T, from its species' ThermoData."""
    thermo = derived["thermo"]
    assert thermo["kind"] == "ThermoData-Kc/1"

    def g(entry):
        td = ThermoData(Tdata=(entry["Tdata"], "K"), Cpdata=(entry["Cpdata"], "J/(mol*K)"),
                        H298=(entry["H298"], "J/mol"), S298=(entry["S298"], "J/(mol*K)"),
                        Cp0=(entry["Cp0"], "J/(mol*K)"), CpInf=(entry["CpInf"], "J/(mol*K)"))
        return td.get_free_energy(temperature)
    dg = math.fsum(g(e) for e in thermo["products"]) - math.fsum(g(e) for e in thermo["reactants"])
    dn = len(thermo["products"]) - len(thermo["reactants"])
    return math.exp(-dg / (constants.R * temperature)) * (thermo["P0"] / (constants.R * temperature)) ** dn


def cross_method_document():
    doc, deck = pe()
    table = K.kernel_rates(doc, deck.spec, T)
    terms = mp.kernel_input(table, deck.spec, [C, D1, DR, T1], temperature=T, phi=1.0, volume_ratio=1.0,
                            **LABELS)
    by_owner = {"Re": terms["recombination"] + terms["transfer"].get("Re", []), "Ut": terms["transfer"]["Ut"]}
    kernel_of = {}
    for kid in (C, D1, DR, T1):
        for t in table["kernels"][kid]["templates"]:
            kernel_of[t["template"]] = kid
    cases = []
    for name, kernel, owner, molecules in CASES:
        smiles = [s for pool in sorted(molecules) for s in molecules[pool]]
        engine = me.MeltEngine(deck.spec, me.PayloadCatalog(doc["payloads"]), T, np.random.default_rng(0),
                               state=S.witness_state(smiles), volume_m3=VOLUME)
        rates = per_template(engine)
        for term in by_owner[owner]:
            templates = term["template"].split("+")
            if kernel_of[templates[0]] != kernel:
                continue
            pairs, same = site_pairs(engine, set(templates))
            chains = {pool: {} for pool in molecules}
            for pool, ss in molecules.items():
                for s in ss:
                    n = sum(1 for ch in s if ch == "C")
                    chains[pool][str(n)] = chains[pool].get(str(n), 0) + 1
            cases.append({"name": name, "kernel": kernel, "owner": owner, "term": term["template"],
                          "chains": chains, "kmc_events_per_s": math.fsum(rates.get(t, 0.0) for t in templates),
                          "site_pairs": pairs, "same_molecule_pairs": same,
                          "one_pool": len(molecules) == 1})
    # D1's term template and its thermo-derived reverse, by record
    payloads = {p["record_event_id"]: p for p in doc["payloads"]}
    (d1_term,) = [t for t in terms["transfer"]["Re"] if kernel_of[t["template"].split("+")[0]] == D1]
    forward = d1_term["template"]
    record = next(p for p in doc["payloads"] if p["record_event_id"][4:16] == forward.split(":")[0])
    reverse = payloads[record["reverse_of"]]
    derived = reverse["templates"][0]["derived"]
    reverse_name = f"{reverse['record_event_id'][4:16]}:{reverse['templates'][0]['key'][:12]}"
    (dr_term,) = [t for t in terms["transfer"]["Ut"] if t["template"] == reverse_name]
    return {"schema": "kmc-mom-cross-method/1",
            "generator": "test/rmgpy/kmc/momCrossMethodTest.py cross_method_document",
            "catalog": "i213_pe_catalog.json.zlib.b64", "deck": "i213_pe_toy_deck.py", "temperature_K": T,
            "volume_m3": VOLUME, "avogadro": constants.Na,
            "periods_per_unit": table["moment_unit"]["periods_per_declared_unit"],
            "terms": {"recombination": terms["recombination"], "transfer": terms["transfer"]},
            "cases": cases,
            "detailed_balance": {"forward_term": forward, "reverse_term": dr_term["template"],
                                 "kc_forward": kc(derived, T), "reference_state": derived["thermo"]["kind"]}}


@pytest.fixture(scope="module")
def document():
    return cross_method_document()


def test_the_fixture_is_what_the_kmc_code_computes(document):
    assert json.loads(FIXTURE.read_text()) == document
    assert hashlib.sha256(FIXTURE.read_bytes()).hexdigest() == CROSS_METHOD_SHA256


def test_every_case_has_events_and_its_term(document):
    names = {c["name"] for c in document["cases"]}
    assert names == {name for name, *_ in CASES}
    for case in document["cases"]:
        assert case["kmc_events_per_s"] > 0.0
    like = [c for c in document["cases"] if c["one_pool"]]
    assert all(c["same_molecule_pairs"] > 0 for c in like)          # the mean field's extra pairs


def test_the_moms_forward_over_reverse_coefficient_is_kc(document):
    """R + R -> P + U at 700 K: the MoM D1 coefficient (ordered pairs, as
    RMG's k for identical reactants) over its reverse's is Kc of the record's
    thermo at C0 = P0 / R T."""
    db = document["detailed_balance"]
    terms = {t["template"]: t for ts in document["terms"]["transfer"].values() for t in ts}
    kf, kr_ = terms[db["forward_term"]]["k_pair"], terms[db["reverse_term"]]["k_pair"]
    assert "sites" not in terms[db["forward_term"]] and "sites" not in terms[db["reverse_term"]]
    assert kf / kr_ == pytest.approx(db["kc_forward"], rel=1e-9)


def write_fixture():
    """Regenerate the fixture (by hand, after a deliberate change)."""
    FIXTURE.write_text(json.dumps(cross_method_document(), indent=1, sort_keys=True) + "\n")
    return hashlib.sha256(FIXTURE.read_bytes()).hexdigest()
