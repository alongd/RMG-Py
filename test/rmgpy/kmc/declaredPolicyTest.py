"""I-165 S2b: an explicit chemistry policy per witness (domain, certificate
eligibility, recertification, the ring-edit generation filter), correction
C-1 (the depository growth unit follows the chain provider), and the
legacy/declared member parity harness (``declared_windows_parity``).

Database-free unless a test says [DB: <family>]; a database test loads that
one kinetics family through ``one_family.family_db`` (``_only_family``
evicts any other first, so one process never holds two)."""

import gc
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import one_family  # noqa: E402

DATA = HERE / "data"
needs_db = pytest.mark.skipif(not one_family.available(),
                              reason=f"needs the RMG database at {one_family.DB}")


def _only_family(name):
    """``family_db(name)`` after evicting every other cached family."""
    for other in [k for k in one_family._DB if k != name]:
        del one_family._DB[other]
    gc.collect()
    return one_family.family_db(name)


def _pp():
    from declaredChainTest import pp_mapping, spec_of
    from rmgpy.kmc import declared_chain as dc
    return dc.DeclaredChainProvider(spec_of(pp_mapping()))


def _pmma():
    from declaredChainTest import pmma_mapping, spec_of
    from rmgpy.kmc import declared_chain as dc
    return dc.DeclaredChainProvider(spec_of(pmma_mapping()))


def _pp_specs(units=4):
    from rmgpy.kmc import declared_chain as dc
    from rmgpy.kmc import windows as win
    return [win.DeclaredSpec(_pp().name, dc.ChainSpec(units, head="artificial",
                                                       tail="artificial")), ("small", "H")]


def _witness_of(specs):
    """A participants-only witness of provider builds (policy resolution
    reads nothing else)."""
    from rmgpy.kmc import boundaries as b
    from rmgpy.kmc import windows as win
    _molecules, participants = b.witness_reactants(win.builds_of(specs))
    return SimpleNamespace(participants=participants)


# --------------------------------------------------------------------------
# 1. the comparator (database-free, fabricated rows)
# --------------------------------------------------------------------------
KEY = "ref|Fam|forward"


def _row(pair, status="fixed_point", **changes):
    row = {"key": KEY, "pair": pair, "run": True, "status": status, "family": "Fam",
           "window": "w-" + pair, "view": "v-" + pair, "span": {"0": [0, 3]},
           "contexts": {"(0,'head'):declared": "fixed_point"}, "deferred": []}
    row.update(changes)
    return row


def _expected(members=2, guarded=(), deferred=0):
    return {KEY: {"members": members, "guarded": sorted(guarded), "deferred": deferred}}


GUARD = "refused:whole_species_search: Node_X label *3 candidates [7, 9]"


@pytest.fixture(scope="module")
def parity():
    import declared_windows_parity as parity
    return parity


def test_comparator_equal_rows_pass(parity):
    rows = [_row("a"), _row("b")]
    report = parity.compare(rows, [dict(r) for r in rows], _expected())
    assert report.ok and report.failures == [] and report.exit_code == 0


@pytest.mark.parametrize("side", ["legacy", "declared"])
def test_comparator_an_unpaired_member_fails(parity, side):
    both = [_row("a"), _row("b")]
    one = [_row("a")]
    legacy, declared = (both, one) if side == "legacy" else (one, both)
    report = parity.compare(legacy, declared, _expected())
    assert not report.ok and report.exit_code == 1
    assert any("without a partner" in f and "b" in f for f in report.failures)


@pytest.mark.parametrize("field,value", [("status", "refused:domain"), ("window", "other"),
                                         ("view", "other"), ("span", {"0": [0, 4]}),
                                         ("contexts", {"(0,'head'):declared": "ceiling"})])
def test_comparator_a_disagreeing_pair_fails(parity, field, value):
    report = parity.compare([_row("a"), _row("b")], [_row("a"), _row("b", **{field: value})],
                            _expected())
    assert not report.ok
    assert any("b" in f and field in f for f in report.failures), report.failures


def test_comparator_contexts_compare_only_end_classes_on_both_sides(parity):
    legacy = _row("a", contexts={"(0,'head'):declared": "fixed_point",
                                 "(0,'tail'):radical": "fixed_point"})
    declared = _row("a", contexts={"(0,'head'):declared": "fixed_point"})
    report = parity.compare([legacy, _row("b")], [declared, _row("b")], _expected())
    assert report.ok, report.failures


def test_comparator_a_one_sided_guard_fails(parity):
    report = parity.compare([_row("a", status=GUARD), _row("b")], [_row("a"), _row("b")],
                            _expected(guarded=["a"]))
    assert not report.ok
    assert any("guarded on one side" in f for f in report.failures)


def test_comparator_a_both_sided_guard_is_listed_not_failed(parity):
    report = parity.compare([_row("a", status=GUARD, window=None, view=None), _row("b")],
                            [_row("a", status=GUARD, window=None, view=None), _row("b")],
                            _expected(guarded=["a"]))
    assert report.ok, report.failures
    listed = report.listed["guarded"]
    assert [(g["pair"], g["family"], g["node"], g["label"], g["candidates"]) for g in listed] == \
        [("a", "Fam", "Node_X", "*3", "[7, 9]")]


def test_comparator_a_guarded_member_is_never_compared_as_equal(parity):
    """A guard on both sides is listed; it is not a passing comparison of two
    windows (the guard row's window is never read)."""
    report = parity.compare([_row("a", status=GUARD), _row("b")],
                            [_row("a", status=GUARD, window="different"), _row("b")],
                            _expected(guarded=["a"]))
    assert report.ok
    assert report.compared == 1        # only b


@pytest.mark.parametrize("members", [1, 3])
def test_comparator_a_member_count_off_expected_fails(parity, members):
    rows = [_row("a"), _row("b")]
    report = parity.compare(rows, [dict(r) for r in rows], _expected(members=members))
    assert not report.ok
    assert any("member count" in f for f in report.failures)


def test_comparator_zero_members_fails(parity):
    report = parity.compare([], [], _expected(members=0))
    assert not report.ok
    assert any("zero members" in f for f in report.failures)


def test_comparator_a_reference_missing_from_expected_fails(parity):
    rows = [_row("a"), _row("b")]
    report = parity.compare(rows, [dict(r) for r in rows], {})
    assert not report.ok


def test_comparator_a_guarded_set_off_expected_fails(parity):
    rows = [_row("a", status=GUARD), _row("b")]
    report = parity.compare(rows, [dict(r) for r in rows], _expected(guarded=["b"]))
    assert not report.ok
    assert any("guarded set" in f for f in report.failures)
    report = parity.compare(rows, [dict(r) for r in rows], _expected(guarded=[]))
    assert not report.ok


def test_comparator_a_deferred_context_is_listed_not_failed(parity):
    legacy = _row("a", deferred=["(0,'head'):vinylidene"])
    report = parity.compare([legacy, _row("b")], [_row("a"), _row("b")],
                            _expected(deferred=1))
    assert report.ok, report.failures
    assert report.listed["deferred"] == [{"key": KEY, "pair": "a",
                                          "context": "(0,'head'):vinylidene"}]


def test_comparator_rows_not_run_count_but_are_not_compared(parity):
    """``--limit``: every member is a row (the count is total); members past
    the limit are not run on either side and are not compared, and the
    guarded set is checked on the members run."""
    legacy = [_row("a"), _row("b", run=False, status="not_run")]
    declared = [_row("a"), _row("b", run=False, status="not_run")]
    report = parity.compare(legacy, declared, _expected(guarded=["b"]))
    assert report.ok, report.failures
    assert report.compared == 1


# --------------------------------------------------------------------------
# 2. the interim domain (database-free)
# --------------------------------------------------------------------------
def _mol(smiles):
    from rmgpy.molecule import Molecule
    m = Molecule().from_smiles(smiles)
    for k, atom in enumerate(m.atoms):
        atom.id = k
    return m


@pytest.mark.parametrize("name,smiles", [
    ("PP DP3 closed", "CCC(C)CC(C)CC(C)C"),
    ("PP DP3 secondary radical", "CCC(C)C[C](C)CC(C)C"),
    ("H", "[H]"),
    ("H2", "[H][H]"),
    ("PMMA DP2 closed", "COC(=O)C(C)(C)CC(C)(C)C(=O)OC"),
])
def test_interim_check_accepts(name, smiles):
    from rmgpy.kmc import domain as dom
    _pp(), _pmma()          # registered declarations (PMMA declares O)
    verdict = dom.declared_interim_check(_mol(smiles))
    assert verdict.accepted, (name, verdict.reason)
    assert verdict.kind == "declared_interim"
    sites = [a.id for a in _mol(smiles).atoms if a.radical_electrons]
    assert sorted(verdict.details["radical_sites"]) == sites


def _h3():
    from rmgpy.molecule import Atom, Bond, Molecule
    atoms = [Atom(element="H", radical_electrons=0, lone_pairs=0) for _ in range(3)]
    m = Molecule(atoms=atoms)
    m.add_bond(Bond(atoms[0], atoms[1], order=1))
    atoms[2].radical_electrons = 1
    m.multiplicity = 2
    return m


def _o_one_lone_pair():
    m = _mol("CO")
    next(a for a in m.atoms if a.symbol == "O").lone_pairs = 1
    return m


def _one_site_triplet():
    m = _mol("CC[CH2]")
    m.multiplicity = 3
    return m


@pytest.mark.parametrize("name,make,reason", [
    ("carbene (u2)", lambda: _mol("[CH2]"), "u2"),
    ("charged", lambda: _mol("C[CH2+]"), "charged"),
    ("hydrogen-only H3", _h3, "hydrogen-only"),
    ("O with one lone pair", _o_one_lone_pair, "lone pair"),
    ("one-site triplet", _one_site_triplet, "multiplicity convention"),
])
def test_interim_check_refuses(name, make, reason):
    from rmgpy.kmc import domain as dom
    _pp(), _pmma()
    verdict = dom.declared_interim_check(make())
    assert not verdict.accepted, name
    assert reason in verdict.reason, (name, verdict.reason)


def test_interim_check_refuses_an_undeclared_element():
    from rmgpy.kmc import domain as dom
    _pp()
    verdict = dom.declared_interim_check(_mol("CCCl"))
    assert not verdict.accepted and "element" in verdict.reason


def test_interim_check_does_not_apply_the_ps_shape_rules():
    """Branched backbones (PP), backbone radicals two bonds apart and
    pendant-free chains are accepted: the interim widening is the point."""
    from rmgpy.kmc import domain as dom
    for smiles in ("CC(C)(C)CC", "C[CH]C[CH]C"):
        assert dom.declared_interim_check(_mol(smiles)).accepted, smiles
        assert not dom.check(_mol(smiles)).accepted, smiles


def test_check_is_unchanged_on_the_domain_test_molecules():
    """``dom.check`` on the inputs of the existing domain tests equals the
    S2a output frozen in the red commit (``capture_domain_before.py``)."""
    import capture_domain_before as cap
    frozen = json.loads((DATA / "domain_check_before.json").read_text())
    assert cap.verdicts() == frozen


def test_universal_check_is_the_prefix_of_check():
    from rmgpy.kmc import domain as dom
    assert dom.universal_check(_mol("CC(O)C"), ("C", "H")).reason == \
        "an element other than C and H"
    assert dom.universal_check(_mol("CCCC"), ("C", "H")) is None
    assert dom.universal_check(_mol("[H]"), ("C", "H")).kind == "small:H"


# --------------------------------------------------------------------------
# 3. policy resolution (database-free)
# --------------------------------------------------------------------------
def test_policy_for_legacy_witnesses():
    from rmgpy.kmc import boundaries as b
    from rmgpy.kmc import windows as win
    for specs in ([b.PsChainSpec(3), ("small", "H")], [("alkane", 8, None, None)],
                  [("small", "H"), ("small", "H2")]):
        assert win.policy_for(_witness_of(specs)) is win.LEGACY_V1, specs


def test_policy_for_a_declared_witness():
    from rmgpy.kmc import windows as win
    assert win.policy_for(_witness_of(_pp_specs())) is win.DECLARED_INTERIM


def test_policy_for_refuses_a_mixed_witness():
    from rmgpy.kmc import boundaries as b
    from rmgpy.kmc import windows as win
    with pytest.raises(win.WindowError, match="policy"):
        win.policy_for(_witness_of([b.PsChainSpec(2), _pp_specs()[0]]))


def test_policy_objects():
    from rmgpy.kmc import domain as dom
    from rmgpy.kmc import windows as win
    assert (win.LEGACY_V1.name, win.LEGACY_V1.ring_edits) == ("legacy_v1", "exclude")
    assert (win.DECLARED_INTERIM.name, win.DECLARED_INTERIM.ring_edits) == \
        ("declared_interim_v1", "exclude")
    m = _mol("CC(C)(C)CC")
    assert win.LEGACY_V1.domain.check(m) == dom.check(m)
    assert win.DECLARED_INTERIM.domain.check(m) == dom.declared_interim_check(m)
    assert win.DECLARED_INTERIM.domain.spin_verdict(m) == dom.spin_verdict(m)
    with pytest.raises(Exception):
        win.LEGACY_V1.name = "x"


def test_an_explicit_policy_wins(monkeypatch):
    """``species_domain_refusal(witness, policy=)``: the explicit policy is
    used, never the resolved one."""
    from rmgpy.kmc import windows as win
    calls = []

    def spy(molecule):
        calls.append(molecule)
        return win.LEGACY_V1.domain.check(molecule)
    policy = win.ChemistryPolicy("spy", win.PolicyDomain("spy", spy, win.LEGACY_V1.domain
                                                         .spin_verdict), "exclude")
    import declared_parity_support as sup
    from rmgpy.kmc import boundaries as b
    witness = sup.identity_witness(win.builds_of([b.PsChainSpec(2)])[0])
    assert win.species_domain_refusal(witness, policy=policy) is None
    assert calls
    resolved = win.species_domain_refusal(witness)
    assert resolved is None and len(calls) == 2          # one species, both sides; unchanged


# --------------------------------------------------------------------------
# 4-7. [DB: H_Abstraction] eligibility, recertification, PP + H, ring edits
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def habs():
    if not one_family.available():
        pytest.skip(f"needs the RMG database at {one_family.DB}")
    from rmgpy.kmc import certificates as cm
    db = _only_family("H_Abstraction")
    fams = db.kinetics.families
    return fams["H_Abstraction"], cm.Certifier(fams, max_channel_sets=16)


@pytest.fixture(scope="module")
def pp_cases(habs):
    from rmgpy.kmc import windows as win
    family, certifier = habs
    stats = {}
    cases = win.member_cases(certifier, family, win.builds_of(_pp_specs()), "forward",
                             label="PP4+H", stats=stats)
    return cases, stats


def _pp_certificate(habs, case, **kwargs):
    from rmgpy.kmc import windows as win
    family, certifier = habs
    return certifier.certify(family, win.template_molecules(case.witness, "forward"),
                             member_labels=case.wanted, **kwargs)


@needs_db
def test_eligibility_follows_the_declared_policy(habs, pp_cases):
    from rmgpy.kmc import windows as win
    case = pp_cases[0][0]
    cert = _pp_certificate(habs, case, domain=win.DECLARED_INTERIM.domain)
    assert cert["eligible"] is True
    assert cert["domain_policy"] == "declared_interim_v1"
    legacy = _pp_certificate(habs, case)
    assert legacy["eligible"] is False and "domain_policy" not in legacy
    assert {(ch.get("eligibility") or {}).get("disposition") for ch in legacy["channels"]} == \
        {"refused:domain"}


@needs_db
def test_a_legacy_certificate_without_a_domain_is_byte_identical(habs):
    import capture_c1_before as cap
    family, certifier = habs
    frozen = json.loads((DATA / "legacy_ps_dp3_h_certificate.json").read_text())
    live = json.loads(cap.canonical(cap.ps_dp3_h_certificate(certifier, family)))
    assert live["member"] == frozen["member"]
    assert cap.canonical(live["certificate"]) == cap.canonical(frozen["certificate"])


def _longer(case):
    """The member's witness one period longer at its first artificial port,
    and the identity map of the surviving atoms."""
    from rmgpy.kmc import boundaries as b
    from rmgpy.kmc import windows as win
    port = next(p for p in b.vector_ports(case.witness) if case.witness.port(p).artificial)
    longer, _record = b.extend(case.witness, port)
    return win.template_molecules(longer, "forward")


@needs_db
def test_recertification_follows_the_policy(habs, pp_cases):
    from rmgpy.kmc import windows as win
    family, certifier = habs
    for case in pp_cases[0]:
        cert = _pp_certificate(habs, case, domain=win.DECLARED_INTERIM.domain)
        mols = _longer(case)
        present = {a.id for m in mols for a in m.atoms}
        amap = {row[0]: row[0] for row in cert["witness"]["atoms"] if row[0] in present}
        same = certifier.recertify_certificate(family, cert, mols, amap,
                                               domain=win.DECLARED_INTERIM.domain)
        if same.get("same_selection"):
            break
    else:
        pytest.fail("no PP member recertifies with the same selection one period longer")
    assert not same["refusals"]
    changed = certifier.recertify_certificate(family, cert, mols, amap)
    assert not changed.get("same_selection")
    assert any("refused:domain_policy_changed" in r for r in changed["refusals"])


@pytest.fixture(scope="module")
def pp_runs(habs, pp_cases):
    """Every PP + H member through ``run_member(contexts=True)`` with the
    provider's span bound (the declared policy is resolved)."""
    from rmgpy.kmc import windows as win
    family, certifier = habs
    bound = win.span_bound(win.family_constants(family), _pp())
    out = []
    for case in pp_cases[0]:
        out.append((case, win.run_member(certifier, family, case, bound=bound, contexts=True,
                                         specs=_pp_specs())))
    return bound, out


@needs_db
def test_pp_plus_h_is_not_vacuous(habs, pp_cases, pp_runs):
    from rmgpy.kmc import windows as win
    cases, stats = pp_cases
    expected = json.loads((DATA / "declared_parity_expected.json").read_text())
    assert len(cases) > 0
    assert len(cases) == expected["pp4+h|H_Abstraction|forward"]["members"]
    assert stats["ring_change_mappings"] == 0
    bound, runs = pp_runs
    statuses = [r["status"] for _c, r in runs]
    assert not [s for s in statuses if s.startswith("refused:domain")], statuses
    assert not [s for s in statuses if "provider" in s or "roster" in s or "legacy_spec" in s]
    extended = [(c, r) for c, r in runs if r["status"] == "fixed_point"
                and any(d > 0 for _p, d in r["depths"])]
    assert extended, statuses
    case, res = extended[0]
    assert res["span_bound_held"]
    reached = [ctx for _c, r in runs for ctx in r.get("contexts") or ()
               if any(end in ("declared", "radical") for _p, end in ctx["assignment"])
               and not ctx["status"].startswith("refused")]
    assert reached, [r.get("contexts") for _c, r in runs][:3]
    assert all(r.get("seconds") is not None for _c, r in runs)


@needs_db
def test_pp_plus_h_runs_under_the_declared_policy(habs, pp_cases):
    """The certificate path sees the declared policy: a run_member forced to
    LEGACY_V1 is refused by the domain, the resolved one is not."""
    from rmgpy.kmc import windows as win
    family, certifier = habs
    case = pp_cases[0][0]
    legacy = win.run_member(certifier, family, case, contexts=False, policy=win.LEGACY_V1)
    assert legacy["status"] == "refused:domain"


@needs_db
def test_ring_edit_policy(habs):
    from rmgpy.kmc import boundaries as b
    from rmgpy.kmc import windows as win
    family, certifier = habs
    builds = win.builds_of([b.PsChainSpec(3), ("small", "H")])

    def labels(cases):
        return sorted(c.label for c in cases)
    default = win.member_cases(certifier, family, builds, "forward")
    assert labels(win.member_cases(certifier, family, builds, "forward",
                                   policy=win.LEGACY_V1)) == labels(default)
    admit = win.ChemistryPolicy("admit_test", win.LEGACY_V1.domain, "admit")
    stats = {}
    admitted = win.member_cases(certifier, family, builds, "forward", policy=admit)
    every = win.member_cases(certifier, family, builds, "forward", skip_ring_changes=False)
    assert labels(admitted) == labels(every)
    excluded = win.member_cases(certifier, family, builds, "forward", stats=stats)
    assert stats["ring_change_mappings"] > 0
    assert sum(len(c.moves) for c in admitted) > sum(len(c.moves) for c in excluded)


# --------------------------------------------------------------------------
# 8-9. [DB: R_Recombination] C-1 and local parity
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def recombination():
    if not one_family.available():
        pytest.skip(f"needs the RMG database at {one_family.DB}")
    return _only_family("R_Recombination")


@pytest.fixture(scope="module")
def c16_parity(recombination, parity):
    return parity.run_reference("c16", "R_Recombination", "reverse", limit=10)


@needs_db
def test_c1_before_is_frozen_and_after_is_live(c16_parity):
    frozen = json.loads((DATA / "c1_before.json").read_text())["members"]
    rows = c16_parity["corrections"]["C-1"]["rows"]
    assert rows
    for row in rows:
        assert row["before"] == frozen[row["member"]], row["member"]
        assert all(ch["repeat_unit"] == {"C": 1, "H": 2} for ch in row["after"]["channels"])
        assert all(ch["repeat_unit"] == {"C": 8, "H": 8} for ch in row["before"]["channels"])
    changed = [r for r in rows if [c["verdict"] for c in r["before"]["channels"]] !=
               [c["verdict"] for c in r["after"]["channels"]]]
    summary = c16_parity["corrections"]["C-1"]["summary"]
    if not changed:
        assert "no depository_exclusion verdict changes" in summary
        assert "before" in summary and "after" in summary


@needs_db
@pytest.mark.parametrize("reference,limit", [("c16", 10), ("dp8_closed", 20)])
def test_local_parity(recombination, parity, c16_parity, reference, limit):
    result = c16_parity if reference == "c16" else \
        parity.run_reference(reference, "R_Recombination", "reverse", limit=limit)
    report = result["report"]
    assert report["failures"] == [], report["failures"]
    key = f"{reference}|R_Recombination|reverse"
    expected = json.loads((DATA / "declared_parity_expected.json").read_text())[key]
    assert result["members"]["legacy"] == result["members"]["declared"] == \
        expected["members"] > 0
