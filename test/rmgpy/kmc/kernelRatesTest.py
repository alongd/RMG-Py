"""Red tests (I-213, S1): MoM kernel rates from a kMC template catalog
(``rmgpy.kmc.kernel_rates``).

The PE fixture is the i207 catalog (64 loaded templates; payloads unchanged,
``refusals``/``gates``/``completeness`` dropped) with the PE toy deck it is
bound to and the well-mixed run's template counts.  Expected values come from
the I-204 design (``/home/alon/runs/i204/design.md`` §1.4, §2): the backbone
homolysis has 2n - 7 sites per chain of n C2 units, i.e. n - 7 per CH2
period; the template classes and their firings are the design's table.  The
toy catalog (``melt_engine_support.pe_templates``) adds a beta scission and
two templates that make branches.
"""

import base64
import json
import math
import zlib
from pathlib import Path

import pytest

import melt_engine_support as S
from meltRunTest import pe_catalog_file

from rmgpy.kmc import kernel_rates as kr_
from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import run as kr
from rmgpy.polymer_spec import LinearHomopolymerSpec

FIXTURES = Path(__file__).parent / "fixtures"
T = 700.0
HOMOLYSIS = "5a88d9573d8b:f5338ad0db2e"
RECOMBINATION = "1af8fada657a:e977ec702e62"


def pe_catalog():
    text = (FIXTURES / "i213_pe_catalog.json.zlib.b64").read_bytes().replace(b"\n", b"")
    return json.loads(zlib.decompress(base64.b64decode(text)))


@pytest.fixture(scope="module")
def pe():
    doc = pe_catalog()
    deck = kr.read_deck(FIXTURES / "i213_pe_toy_deck.py")
    run = json.loads((FIXTURES / "i213_pe_wellmixed_counts.json").read_text())
    table = kr_.kernel_rates(doc, deck.spec, T, run=run)
    weights = {t.name: t.weight for t in me.PayloadCatalog(doc["payloads"]).templates(T)}
    return {"doc": doc, "deck": deck, "table": table, "weights": weights}


def placed(table):
    return {k: sorted(x["template"] for x in v["templates"]) for k, v in table["kernels"].items()}


def toy_name(name):
    """A toy template's loaded name: its record id (``evt_<name>``) [4:16], then its key."""
    return lambda loaded: loaded.split(":")[0] == name[:12]


# --------------------------------------------------------------------------
# the PE tables of the design
# --------------------------------------------------------------------------
def test_the_moment_unit_is_the_primitive_period(pe):
    unit = pe["table"]["moment_unit"]
    assert unit["periods_per_declared_unit"] == 2
    ch2 = pe["deck"].spec.repeat_unit.composition.molar_mass_g_mol() / 2
    assert unit["molar_mass_g_mol"] == pytest.approx(ch2, rel=1e-12)


def test_backbone_homolysis_has_n_minus_7_sites_per_chain_of_n_periods(pe):
    kernel = pe["table"]["kernels"]["homolysis:substrate=chain"]
    w = pe["weights"][HOMOLYSIS]
    assert [x["template"] for x in kernel["templates"]] == [HOMOLYSIS]
    assert kernel["k_mu1"] == pytest.approx(w, rel=1e-12)
    assert kernel["k_mu0"] == pytest.approx(-7.0 * w, rel=1e-12)
    row = next(r for r in pe["table"]["templates"] if r["template"] == HOMOLYSIS)
    assert row["participants"][0]["chain"]["a"] == 1.0 and row["participants"][0]["chain"]["b"] == -7.0
    assert row["phi_thinned"] is True


def test_every_loaded_template_is_placed_or_listed(pe):
    names = set(pe["weights"])
    assert len(names) == 64
    in_kernels = {t for ts in placed(pe["table"]).values() for t in ts}
    listed = [e["template"] for e in pe["table"]["kmc_only"]]
    assert len(listed) == len(set(listed))
    assert in_kernels | set(listed) == names
    assert not in_kernels & set(listed)


def test_the_kernel_classes_are_the_designs(pe):
    counts = {k: len(v) for k, v in placed(pe["table"]).items()}
    assert counts == {
        "homolysis:substrate=chain": 1,
        "recombination:radical=radical_end:CH2+radical=radical_end:CH2": 1,
        "disproportionation:acceptor=radical_end:CH2+donor=radical_end:CH2": 2,
        "disproportionation:acceptor=radical_end:CH2+donor=radical_interior:CH": 5,
        "disproportionation:acceptor=radical_interior:CH+donor=radical_end:CH2": 2,
        "disproportionation:acceptor=radical_interior:CH+donor=radical_interior:CH": 7,
        "reverse_disproportionation:alkene=alkene_internal:CH=CH+donor=chain": 12,
        "reverse_disproportionation:alkene=alkene_terminal:CH2=CH+donor=chain": 4,
        "h_transfer:acceptor=radical_end:CH2+donor=chain": 4,
        "h_transfer:acceptor=radical_interior:CH+donor=chain": 4,
        "h_transfer:acceptor=radical_end:CH2+donor=alkene_terminal:CH2=CH2": 1,
        "intra_h_transfer:substrate=radical_end:CH2": 3,
        "intra_h_transfer:substrate=radical_interior:CH": 13,
    }


def test_the_kmc_only_templates_and_their_reasons(pe):
    got = sorted((e["template"], e["reason"]) for e in pe["table"]["kmc_only"])
    assert got == [("29c6bc22a5c6:4e4311c6c096", "no_site_on_declared_states"),   # homolysis at a branch point
                   ("ae2879012a1f:b09ba07d22c5", "branching"),                    # end + mid recombination
                   ("afa69acb7ca3:f0706b832fa9", "no_site_on_declared_states"),   # homolysis between branch points
                   ("d1cf2b57ab48:b326207b6cc5", "small_species_participant"),    # vinyl radical acceptor
                   ("efa40220a975:c6072eed23bf", "branching")]                    # mid + mid recombination


def test_a_bimolecular_coefficient_divides_by_the_role_swap(pe):
    kernel = pe["table"]["kernels"]["recombination:radical=radical_end:CH2+radical=radical_end:CH2"]
    (term,) = kernel["templates"]
    assert term["w_over_divisor"] == pytest.approx(pe["weights"][RECOMBINATION] / 2.0, rel=1e-12)
    assert term["S0"]["per_feature"] == term["S1"]["per_feature"] == 1


def test_the_hazard_shares_are_the_designs(pe):
    h = pe["table"]["hazard"]
    assert h["events"] == 23133 and h["net_events"] == 1175
    assert h["leak_fraction"] == pytest.approx(0.17096817018036953, rel=1e-12)
    assert h["kmc_only_share_gross"] == 0.0
    by = h["by_kernel_gross"]
    assert by["homolysis:substrate=chain"] == 1062
    assert by["recombination:radical=radical_end:CH2+radical=radical_end:CH2"] == 759
    assert by["disproportionation:acceptor=radical_interior:CH+donor=radical_interior:CH"] == 9946
    assert by["reverse_disproportionation:alkene=alkene_internal:CH=CH+donor=chain"] == 9954
    assert sum(by.values()) == 23133


def test_the_command_line_writes_the_table_beside_the_catalog(tmp_path, capsys):
    catalog = tmp_path / "pe.catalog.json"
    catalog.write_text(json.dumps(pe_catalog()))
    code = kr_.main([str(FIXTURES / "i213_pe_toy_deck.py"), str(catalog), "--temperature", "700"])
    assert code == 0
    table = json.loads((tmp_path / "pe.catalog.kernel_rates.json").read_text())
    assert table["schema"] == kr_.SCHEMA and "homolysis:substrate=chain" in table["kernels"]
    assert "written" in capsys.readouterr().out


def test_a_catalog_of_another_chemistry_is_refused(tmp_path, capsys):
    doc = pe_catalog()
    doc["provenance"]["chemistry_sha256"] = "0" * 64
    catalog = tmp_path / "other.json"
    catalog.write_text(json.dumps(doc))
    assert kr_.main([str(FIXTURES / "i213_pe_toy_deck.py"), str(catalog), "--temperature", "700"]) == 2
    assert "chemistry digest" in capsys.readouterr().err


# --------------------------------------------------------------------------
# the toy catalog: beta scission, branching and a template of no kernel
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    path = pe_catalog_file(tmp_path_factory.mktemp("toy"))
    doc = json.loads(Path(path).read_text())
    spec = LinearHomopolymerSpec.from_mapping(S.pe_mapping({10: 1}))
    return kr_.kernel_rates(doc, spec, T)


def test_beta_scission_of_a_mid_chain_radical_is_its_own_kernel(toy):
    (name,) = placed(toy)["beta_scission:substrate=radical_interior:CH"]
    assert toy_name("beta_scission")(name)


def test_a_template_of_no_kernel_type_is_listed_not_dropped(toy):
    reasons = {e["template"].split(":")[0]: e["reason"] for e in toy["kmc_only"]}
    assert reasons["intra_mid_end"[:12]] == "no_kernel"            # intramolecular radical coupling
    assert reasons["mid_end_recombination"[:12]] == "branching"


def test_the_toy_homolysis_counts_two_ch2_neighbours_on_each_side(toy):
    (row,) = [r for r in toy["templates"] if toy_name("homolysis")(r["template"])]
    chain = row["participants"][0]["chain"]
    assert (chain["a"], chain["b"]) == (1.0, -5.0)       # ports are interior CH2: bonds 3..N-3
    assert not math.isnan(row["weight"])


def test_a_copolymer_declaration_is_refused():
    with pytest.raises(kr_.KernelRatesError, match="linear homopolymer"):
        kr_.kernel_rates({"schema": "kmc-template-catalog/1", "payloads": []}, object(), T)


# --------------------------------------------------------------------------
# I-214: product laws, families and depository overrides
# --------------------------------------------------------------------------
def law_sites(rows):
    return sum(m for m, _p in rows)


def test_the_pe_homolysis_law_cuts_four_periods_from_either_end(pe):
    (term,) = pe["table"]["kernels"]["homolysis:substrate=chain"]["templates"]
    law = term["law"]
    assert law["anchor"] == "chain" and law["far"]["edge"] == 4
    assert all(law["exact"][str(n)] == [] for n in range(1, 8))
    assert law["exact"]["8"] == [[1, [["radical_end:CH2", 4], ["radical_end:CH2", 4]]]]
    assert law["far"]["bulk"] == [[1, [["radical_end", 0, 1, 1], ["radical_end", 1, -1, -1]]]]


def test_every_law_has_the_tables_site_count(pe):
    """A law's multiplicities are the counted sites: the chain polynomial on
    every length it holds, an end feature's profile or far count on every
    fragment, an interior feature's far count."""
    nbb = pe["table"]["moment_unit"]["backbone_atoms_per_period"]
    checked = 0
    for kernel in pe["table"]["kernels"].values():
        for term in kernel["templates"]:
            if term["law"] is None:
                continue
            if kernel["molecularity"] == 1:
                pairs = [(term["law"], next(r for r in pe["table"]["templates"]
                                            if r["template"] == term["template"])["participants"][0])]
            else:
                pairs = [(term["law"]["S0"], {term["law"]["S0"]["class"]: term["S0"]}),
                         (term["law"]["S1"], {term["law"]["S1"]["class"]: term["S1"]})]
            for law, counts in pairs:
                c = counts[law["class"]]
                if law["anchor"] == "chain":
                    for n, rows in law["exact"].items():
                        if int(n) >= c["exact_from_n"]:
                            assert law_sites(rows) == c["a"] * int(n) + c["b"]
                elif law["anchor"] == "end":
                    for n, rows in law["exact"].items():
                        want = c["profile"].get(f"fragment_backbone_atoms={int(n) * nbb}", [c["per_feature"]])
                        assert [law_sites(rows)] == want, (term["template"], n)
                else:
                    assert law_sites(law["far"]["outcomes"]) == c["per_feature"]
                checked += 1
    assert checked >= 30


def test_terms_without_a_law_say_why(pe):
    terms = [t for k in pe["table"]["kernels"].values() for t in k["templates"]]
    inert = [t for t in terms if t.get("law_inert")]
    refused = [t for t in terms if t.get("law_refused")]
    assert len(terms) == 59 and len(inert) == 34
    assert [t["template"][:12] for t in refused] == ["55593e44e6ca"]       # ethylene donor -> vinyl radical
    assert "more than one feature" in refused[0]["law_refused"]


def test_the_families_of_the_declared_states(toy):
    fam = toy["families"]
    assert fam["radical_end"]["labels"] == {"1": "radical_end:CH3"}
    assert fam["alkene_terminal"]["labels"] == {"2": "alkene_terminal:CH2=CH2"}
    assert fam["alkene_terminal"]["min_n"] == 2 and fam["radical_interior"]["min_n"] == 3
    h = 1.00797
    assert fam["chain"]["molar_mass_offset_g_mol"] == pytest.approx(2 * h, rel=1e-3)
    assert fam["radical_end"]["molar_mass_offset_g_mol"] == pytest.approx(h, rel=1e-3)
    assert fam["alkene_internal"]["molar_mass_offset_g_mol"] == pytest.approx(0.0, abs=1e-9)


def test_an_override_does_not_change_the_counted_sites(tmp_path):
    """The depository's 4-dodecyl and pentyl/heptyl are the template's sites
    at another weight: the counts are those of the catalog without them."""
    from kernelReferenceTest import override_payloads, toy_table
    (tmp_path / "o").mkdir()
    (tmp_path / "p").mkdir()
    with_o = toy_table(tmp_path / "o", override_payloads(S.pe_spec({1: 1})))
    without = toy_table(tmp_path / "p")
    for kid in ("beta_scission:substrate=radical_interior:CH",
                "recombination:radical=radical_end:CH2+radical=radical_end:CH2"):
        a, b = with_o["kernels"][kid]["templates"][0], without["kernels"][kid]["templates"][0]
        assert {k: a[k] for k in ("S0", "S1", "k_per_feature") if k in a} == \
            {k: b[k] for k in ("S0", "S1", "k_per_feature") if k in b}
        assert a["overrides"] and not b["overrides"]
