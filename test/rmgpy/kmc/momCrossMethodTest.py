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
a MET pair's intrinsic k_allowed, before the cage factor), per template its
k = w / divisor, site pairs S0 S1 and same-molecule pairs X (the engine
counts k (S0 S1 - X) / N_A V) when both participants come from one pool,
the mean field sum_T k S0 S1 / N_A V there and the engine's events with two
pools (MoM's
events where its law is exact; summed per template, never over templates'
pooled sites, which invents cross-template pairs), and the MoM law's own
events on those chains (the producer's term law, ``record.laws``), which is
what MoM integrates everywhere.  ``law_exact`` says whether the chains lie
where the law is exact (then law events = mean field).

A chain carrying two radicals (``two_radical_chain``, not a MoM state): the
catalog has no intramolecular radical-radical template (its families:
Disproportionation, H_Abstraction, R_Recombination, intra_H_migration), and
the engine excludes same-molecule pairs, so its two radicals never meet.

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
CROSS_METHOD_SHA256 = "761840f5562848876eee3fcd9effa111e1623c54975dcc18af18fb7e63f8ebd4"
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
    ("d1_two_template_short", D1, "Re", {"Re": [radical(4), radical(60)]}),
)
DIRADICAL = "[CH2]" + "C" * 20 + "[CH]" + "C" * 38           # end radical and interior radical, 60 C
INTRAMOLECULAR_RADICAL_FAMILIES = ("Birad_recombination", "Intra_Disproportionation", "Intra_R_Add_Endocyclic",
                                   "Intra_R_Add_Exocyclic", "Birad_R_Recombination")


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


def site_pairs(engine, template):
    """(S0 S1, X) of one template: all ordered site pairs of the state and
    those inside one molecule, from the engine's index entries."""
    names = [c.name for c in engine.classes]
    s0 = s1 = x = 0
    for entry in engine.entries.values():
        n = [0, 0]
        for (ci, p, _kind), sites in entry.shape.sites.items():
            if names[ci] == template:
                n[p] += len(sites)
        s0, s1, x = s0 + n[0], s1 + n[1], x + n[0] * n[1]
    return s0 * s1, x


def site_totals(engine, template):
    """[S0, S1] of one template over the state."""
    names = [c.name for c in engine.classes]
    out = [0, 0]
    for entry in engine.entries.values():
        for (ci, p, _kind), sites in entry.shape.sites.items():
            if names[ci] == template:
                out[p] += len(sites)
    return out


def template_k(engine, template):
    """w / divisor of one template in the engine (no override on these states)."""
    (cls,) = [c for c in engine.classes if c.name == template]
    return engine.weights[cls.index] / cls.divisor


def law_events(law, chains, pools_of, nav):
    """Events per second of a MoM term's law (the producer's record.laws: per
    pair over N_A V, lengths in periods) on explicit chains, ordered pairs."""
    c = law["coefficients"]
    f0, f1 = (pools_of[f] for f in law["families"])
    total = 0.0
    for n, a in chains.get(f0, {}).items():
        for m, b in chains.get(f1, {}).items():
            n_, m_ = int(n), int(m)
            total += (c[0] + c[1] * n_ + c[2] * m_ + c[3] * n_ * m_) * a * b
    return total / nav


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
    pools_of = {"chain": "P", "radical_end": "Re", "alkene_terminal": "Ut", "radical_interior": "Rm",
                "alkene_internal": "Ui"}
    laws = {law["template"]: law for law in terms["record"]["laws"]}

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
            chains = {pool: {} for pool in molecules}
            for pool, ss in molecules.items():
                for s in ss:
                    n = sum(1 for ch in s if ch == "C")
                    chains[pool][str(n)] = chains[pool].get(str(n), 0) + 1
            # every template of the kernel (the term may carry fewer: one inert on the reachable lengths)
            per = [{"template": t["template"], "in_term": t["template"] in templates,
                    "kmc_events_per_s": rates.get(t["template"], 0.0)} for t in table["kernels"][kernel]["templates"]]
            if len(molecules) == 1:
                for entry in per:
                    t = entry["template"]
                    pairs, same = site_pairs(engine, t)
                    entry.update(k=template_k(engine, t), site_pairs=pairs, same_molecule_pairs=same,
                                 sites=site_totals(engine, t))
                mean_field = math.fsum(e["k"] * e["site_pairs"] for e in per if e["in_term"]) / engine.nav
            else:
                # two pools: no molecule is in both, so the engine's distinct pairs are MoM's pool pairs (a
                # participant's own sites of the other slot pair only with itself, which neither counts)
                mean_field = math.fsum(e["kmc_events_per_s"] for e in per if e["in_term"])
            law = laws[term["template"]]
            shortest = {}
            for v in [v for v in terms["record"]["validity"] if v["template"] == term["template"]]:
                for p in v["participants"]:
                    pool = pools_of[p["family"]]
                    shortest[pool] = max(shortest.get(pool, 0), p["exact_from_periods"])
            cases.append({"name": name, "kernel": kernel, "owner": owner, "term": term["template"],
                          "chains": chains, "kmc_events_per_s": math.fsum(rates.get(t, 0.0) for t in templates),
                          "templates": per,
                          "mean_field_events_per_s": mean_field,
                          "mom_law_events_per_s": law_events(law, chains, pools_of, engine.nav),
                          "law_exact": all(int(n) >= shortest[pool] for pool, ns in chains.items() for n in ns),
                          "one_pool": len(molecules) == 1})
    # the two-radical chain: every template's sites inside the one molecule, and what the engine offers it
    engine = me.MeltEngine(deck.spec, me.PayloadCatalog(doc["payloads"]), T, np.random.default_rng(0),
                           state=S.witness_state([DIRADICAL]), volume_m3=VOLUME)
    names = [c.name for c in engine.classes]
    channels = engine.propensities().channels
    own = []
    for kid in (C, D1, "disproportionation:acceptor=radical_end:CH2+donor=radical_interior:CH",
                "disproportionation:acceptor=radical_interior:CH+donor=radical_end:CH2"):
        for t in table["kernels"][kid]["templates"]:
            pairs, same = site_pairs(engine, t["template"])
            if same:
                own.append({"kernel": kid, "template": t["template"], "same_molecule_pairs": same,
                            "mean_field_events_per_s": template_k(engine, t["template"]) * pairs / engine.nav})
    two_radical = {"smiles": DIRADICAL, "pairs_inside_the_molecule": own,
                   "bimolecular_events_per_s": math.fsum(a for key, a in channels if key[0] in ("bi", "met")),
                   "unimolecular_templates": sorted({names[key[1]] for key, a in channels
                                                     if key[0] == "uni" and a > 0.0}),
                   "catalog_families": sorted({p["family"] for p in doc["payloads"]})}
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
            "cases": cases, "two_radical_chain": two_radical,
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
    assert all(c["mean_field_events_per_s"] == c["kmc_events_per_s"] for c in document["cases"] if not c["one_pool"])
    assert all(sum(t["same_molecule_pairs"] for t in c["templates"]) > 0 for c in like)   # the mean field's extra


def test_each_templates_events_are_its_k_times_its_distinct_pairs(document):
    """The engine's own form, per template, with its weight and divisor:
    k (S0 S1 - X) / N_A V; the mean field drops X, template by template."""
    nav = document["avogadro"] * document["volume_m3"]
    for case in [c for c in document["cases"] if c["one_pool"]]:
        for t in case["templates"]:
            assert t["kmc_events_per_s"] == pytest.approx(t["k"] * (t["site_pairs"] - t["same_molecule_pairs"]) / nav,
                                                          rel=1e-12), (case["name"], t["template"])
        assert case["mean_field_events_per_s"] == pytest.approx(
            math.fsum(t["k"] * t["site_pairs"] for t in case["templates"] if t["in_term"]) / nav, rel=1e-12)


def test_pooling_sites_across_templates_invents_pairs(document):
    """Radicals of 4 and 60 periods under D1's two templates, whose sites
    differ by length: summing the sites over the templates before pairing
    them counts template-A sites against template-B sites, pairs no event
    has.  The per-template mean field is the one MoM's law must match.  (On
    this deck's lengths, from 60 periods, one of the two is inert -- no site
    on a long radical -- and MoM's term carries the other; the 4-period
    radical is outside the reachable lengths, so both fire there.)"""
    (case,) = [c for c in document["cases"] if c["name"] == "d1_two_template_short"]
    assert len(case["templates"]) == 2 and case["one_pool"]
    assert all(t["kmc_events_per_s"] > 0.0 for t in case["templates"])
    pooled = math.prod(sum(t["sites"][p] for t in case["templates"]) for p in (0, 1))
    per_template = sum(t["site_pairs"] for t in case["templates"])
    assert pooled != per_template


def test_the_mom_law_is_the_mean_field_where_it_is_exact(document):
    """MoM integrates its term law (the producer's record) on every length;
    where the chains lie past each participant's exact-from length it is
    the mean field of the engine's templates.  The two-template case has a
    4-period radical below the D1 law's exact length: there the law and the
    mean field differ by what the law-domain attribution reports."""
    for case in document["cases"]:
        if case["law_exact"]:
            assert case["mom_law_events_per_s"] == pytest.approx(case["mean_field_events_per_s"], rel=1e-9), \
                case["name"]
    (short,) = [c for c in document["cases"] if c["name"] == "d1_two_template_short"]
    assert not short["law_exact"]
    assert short["mom_law_events_per_s"] != pytest.approx(short["mean_field_events_per_s"], rel=1e-6)
    assert all(c["law_exact"] for c in document["cases"] if c["name"] != "d1_two_template_short")


def test_a_two_radical_chain_has_no_channel_between_its_radicals(document):
    """One chain with an end and an interior radical: its templates count
    pairs inside it (the mean field would react them), the engine excludes
    same-molecule pairs and offers it no bimolecular event alone, and no
    catalog family joins or disproportionates two radicals of one molecule:
    only its unimolecular H shifts and homolysis remain.  A missing channel
    (kMC only; no MoM pool holds a two-radical chain)."""
    two = document["two_radical_chain"]
    assert two["pairs_inside_the_molecule"]
    assert all(p["mean_field_events_per_s"] > 0.0 for p in two["pairs_inside_the_molecule"])
    assert two["bimolecular_events_per_s"] == 0.0
    assert not set(two["catalog_families"]) & set(INTRAMOLECULAR_RADICAL_FAMILIES)
    doc, _deck = pe()
    family = {p["record_event_id"][4:16]: p["family"] for p in doc["payloads"]}
    assert {family[t.split(":")[0]] for t in two["unimolecular_templates"]} <= {"intra_H_migration",
                                                                               "R_Recombination"}


def test_the_moms_forward_over_reverse_coefficient_is_kc(document):
    """R + R -> P + U at 700 K: the MoM D1 coefficient (ordered pairs, as
    RMG's k for identical reactants) over its reverse's is Kc of the record's
    thermo at C0 = P0 / R T."""
    db = document["detailed_balance"]
    terms = {t["template"]: t for ts in document["terms"]["transfer"].values() for t in ts}
    kf, kr_ = terms[db["forward_term"]]["k_pair"], terms[db["reverse_term"]]["k_pair"]
    assert "sites" not in terms[db["forward_term"]] and "sites" not in terms[db["reverse_term"]]
    assert kf / kr_ == pytest.approx(db["kc_forward"], rel=1e-9)


def test_homolysis_over_recombination_is_kc_at_the_declared_reference_state():
    """Detailed balance with a change in mole number (dn = -1), which the D1
    test (dn = 0) cannot see: R + R -> P at 700 K and its thermo-derived
    reverse, the homolysis template, as the catalog emits them.  MoM's
    recombination coefficient per site pair (k_pair c^2, ordered pairs) over
    its homolysis rate per cut position (w a) is Kc of the recombination
    record's thermo at C0 = P0 / R T: an error in C0 moves the ratio by
    (C0' / C0)^-1."""
    doc, deck = pe()
    table = K.kernel_rates(doc, deck.spec, T)
    (homolysis,) = [k for k, v in table["kernels"].items() if v["type"] == "homolysis"]
    terms = mp.kernel_input(table, deck.spec, [homolysis, C], temperature=T, phi=1.0, volume_ratio=1.0, **LABELS)
    (h,), (c,) = terms["homolysis"], terms["recombination"]
    payload = {p["record_event_id"][4:16]: p for p in doc["payloads"]}
    hp, cp = payload[h["template"].split(":")[0]], payload[c["template"].split(":")[0]]
    assert hp["reverse_of"] == cp["record_event_id"]               # the homolysis is the derived reverse
    derived = hp["templates"][0]["derived"]
    assert len(derived["thermo"]["reactants"]) == 2 and len(derived["thermo"]["products"]) == 1
    assert c["k_pair"] * c["sites_per_chain"] ** 2 / (h["w"] * h["a"]) == pytest.approx(kc(derived, T), rel=1e-9)


def write_fixture():
    """Regenerate the fixture (by hand, after a deliberate change)."""
    FIXTURE.write_text(json.dumps(cross_method_document(), indent=1, sort_keys=True) + "\n")
    return hashlib.sha256(FIXTURE.read_bytes()).hexdigest()
