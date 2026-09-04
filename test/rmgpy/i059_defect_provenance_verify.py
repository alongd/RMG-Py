#!/usr/bin/env python
"""
I-059 verifier -- serialized evidence for ``chain_mass_defect_g_mol``.

Recomputes, from a REAL run artifact and from live producer/consumer code,
that:

1. the consumer's defect-without-provenance census warning REPRODUCES at the
   pinned pre-fix consumer commit against a real pre-fix artifact, with its
   text captured verbatim;
2. it is GONE after the fix FOR THE RIGHT REASON -- proved two-armed, by
   stripping the evidence back out and asserting the warning returns, and by
   forging the evidence four different ways and asserting each is REJECTED;
3. the serialized ``chain_mass_defect_g_mol`` is byte-identical before and
   after (the new producer's value for the same pool, against the value the
   PRE-FIX producer wrote into the on-disk artifact);
4. the mass-curve bias the defect produces once the daughters hold non-zero
   moments is measured and printed, in grams and as a fraction;
5. artifacts written BEFORE this change still load through the post-fix
   consumer, with their behaviour asserted rather than assumed;
6. the RMG and TA test suites are run and their pass/fail sets reported.

Exits 0 only if every check passes.

Portability: the RMG tree under test is derived from this script's own
location, never hardcoded, so the verifier keeps measuring the tree it ships
in after this lands on mainline. Base commits are pinned by SHA, not computed
against mainline, so they keep meaning the same thing forever.

Environment overrides (all optional):
    I059_TA_ROOT     path to the TA repo            (default ~/Code/TA)
    I059_SIDECAR     path to a real pre-fix sidecar (default: newest usable
                     one under ~/runs/RMG/poly_*/chemkin/polymer_pools.json)
    I059_CONDA_SH    conda profile script     (default ~/anaconda3/etc/profile.d/conda.sh)
    I059_TA_ENV      TA conda env name              (default ta_env_py314)
    I059_QUICK=1     scope check 6's RMG arm to the polymer suites instead of
                     the full unit suite (the full suite is the default)

Run in the RMG conda env (rmg_env); the consumer arms are shelled out to the
TA env, which carries a different Python and Cantera on purpose.
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# --- the tree under test: derived from this file, never hardcoded ----------
RMG_ROOT = Path(__file__).resolve().parents[2]
TA_ROOT = Path(os.environ.get(
    "I059_TA_ROOT", os.path.expanduser("~/Code/TA"))).resolve()

# --- PINNED pre-fix bases. By SHA on purpose: computing these against
#     mainline would make them mean "whatever is current", i.e. nothing, the
#     moment this change lands.
#     Moved from 9e7e0c4d5 to d78c7211f when this branch was rebased onto the
#     polymer mainline that carries the i061 DASPK moment error-weight floor.
#     RMG_BASE_FAILURES was measured against the OLD base on DASSL binaries and
#     is therefore stale on BOTH counts; see its docstring.
RMG_BASE_SHA = "d78c7211f"     # RMG-Py polymer mainline, pre-I-059
TA_BASE_SHA = "e29642b"        # TA main, pre-I-059

#: The solver backend is chosen at COMPILE time from ``rmgpy/solver/
#: settings.pxi``, which is GITIGNORED -- so it is invisible to every
#: tracked-content check, and a worktree whose extensions are symlinked to
#: another checkout silently inherits THAT checkout's backend. poly_106 rev 1
#: was integrated on DASSL for exactly this reason while every tracked source
#: compared equal. Any suite comparison must therefore hold the backend fixed
#: across the two trees, or backend differences surface as phantom new
#: failures.
REQUIRED_BACKEND = "pydas.daspk"

CONDA_SH = os.environ.get(
    "I059_CONDA_SH", os.path.expanduser("~/anaconda3/etc/profile.d/conda.sh"))
TA_ENV = os.environ.get("I059_TA_ENV", "ta_env_py314")

PROV_KEY = "chain_mass_defect_provenance"
DEFECT_KEY = "chain_mass_defect_g_mol"
#: The ORIGINAL defect-without-provenance census warning, matched on the
#: phrase unique to it. Its family marker ("side-group mass-contract census")
#: is shared with the I-059 unresolvable-parent warning on purpose -- they
#: are the same census -- so matching the family would conflate the two.
CENSUS_MARKER = "WITHOUT the side_group_homolysis spawn provenance"
UNRESOLVED_MARKER = "is NOT present in this artifact"

FAILURES = []
NOTES = []

#: STALE -- DO NOT TRUST UNTIL RE-MEASURED. Retained only so the re-measured
#: set can be diffed against it.
#:
#: This set was measured at the OLD base 9e7e0c4d5 (28 failed, 3472 passed, 40
#: skipped, 39 deselected, 1:33:52) on binaries that were SYMLINKS into
#: /home/alon/Code/RMG-Py -- i.e. on DASSL. Two things invalidate it:
#:
#:   1. the base moved to d78c7211f, which carries the i061 DASPK moment
#:      error-weight floor and whatever else landed with it;
#:   2. the measurement was taken on the wrong solver backend, so several of
#:      the solver cone/realizability and variance-census entries below may be
#:      DASSL artefacts rather than real pre-existing failures.
#:
#: The bar itself is unchanged and still correct: the suite is not green at the
#: base, so demanding green would be a claim this ticket cannot make and did
#: not cause; the bar is that the after-set introduces nothing NEW. But both
#: sides of that comparison must now be re-measured at d78c7211f with
#: REQUIRED_BACKEND, in trees that each built their own extensions.
RMG_BASE_FAILURES_STALE_9E7E0C4D5 = frozenset((
    "test.rmgpy.data.kinetics.familyTest.TestGenerateReactions::test_beta_scission_generates_polymer_fragments",
    "test.rmgpy.data.kinetics.familyTest.TestGenerateReactions::test_generate_reactions_retains_polymer_identity",
    "test.rmgpy.data.kinetics.familyTest.TestGenerateReactions::test_h_abs_reaction_generation_with_polymer_input",
    "test.rmgpy.solver.solverPolymerJacobianTest.TestScopedJacobianPoly102::test_v2_v3_crash_window_twin_and_episode_invariance",
    "test.rmgpy.solver.solverPolymerTest.TestAcceptedStateVarianceCensus::test_variance_census_fires_on_an_unrealizable_second_moment",
    "test.rmgpy.solver.solverPolymerTest.TestAcceptedStateVarianceCensus::test_variance_census_is_warn_once_and_independent_of_the_cone_half",
    "test.rmgpy.solver.solverPolymerTest.TestHybridPolymerReactor::test_cross_pool_reverse_flux_vanishes_continuously",
    "test.rmgpy.solver.solverPolymerTest.TestHybridPolymerReactor::test_cross_pool_ve_detailed_balance_in_depletion_band",
    "test.rmgpy.solver.solverPolymerTest.TestHybridPolymerReactor::test_legacy_state_layout_and_rhs_golden_frozen",
    "test.rmgpy.solver.solverPolymerTest.TestHybridPolymerReactor::test_qssa_moment_signature",
    "test.rmgpy.solver.solverPolymerTest.TestHybridPolymerReactor::test_radical_qssa_unzip_footprint_confined_to_signature",
    "test.rmgpy.solver.solverPolymerTest.TestHybridPolymerReactor::test_weaklink_u_is_massless",
    "test.rmgpy.solver.solverPolymerTest.TestLegacyUnzipRealizability::test_accepted_state_cone_census_fires_below_the_r81_threshold",
    "test.rmgpy.solver.solverPolymerTest.TestLegacyUnzipRealizability::test_qssa_channel_carries_the_same_chain_debit",
    "test.rmgpy.solver.solverPolymerTest.TestLegacyUnzipRealizability::test_unzip_all_dp1_boundary_is_invariant",
    "test.rmgpy.solver.solverPolymerTest.TestLegacyUnzipRealizability::test_unzip_law_carries_the_chain_termination_debit",
    "test.rmgpy.solver.solverPolymerTest.TestLegacyUnzipRealizability::test_unzip_only_pool_stays_in_the_cone[0.05]",
    "test.rmgpy.solver.solverPolymerTest.TestLegacyUnzipRealizability::test_unzip_only_pool_stays_in_the_cone[1.0]",
    "test.rmgpy.solver.solverPolymerTest.TestLegacyUnzipRealizability::test_unzip_only_pool_stays_in_the_cone[100.0]",
    "test.rmgpy.solver.solverPolymerTest.TestLegacyUnzipRealizability::test_unzip_with_scission_stays_in_the_cone",
    "test.rmgpy.solver.solverPolymerTest.TestReleaseAvailabilityGate::test_deprop_sibling_emits_exactly_zero_without_repeat_units",
    "test.rmgpy.solver.solverPolymerTest.TestReleaseAvailabilityGate::test_gate_is_exactly_one_on_the_cone_and_exactly_zero_with_no_units",
    "test.rmgpy.solver.solverPolymerTest.TestReleaseAvailabilityGate::test_legacy_unzip_fabricated_gas_without_repeat_units",
    "test.rmgpy.solver.solverPolymerTest.TestReleaseAvailabilityGate::test_no_gas_is_emitted_beyond_the_units_the_pool_held",
    "test.rmgpy.solver.solverPolymerTest.TestReleaseAvailabilityGate::test_qssa_sibling_is_already_gated_by_its_initiation_term",
    "test.rmgpy.solver.solverPolymerTest.TestRtolNearFloorConviction::test_rtol_1e4_near_floor_conviction_canary",
    "test.rmgpy.tools.simulateTest.SimulateTest::test_liquid",
    "test.rmgpy.tools.simulateTest.SimulateTest::test_minimal",
))


def _solver_backend(tree):
    """Which integrator did ``tree``'s built extensions actually link?

    Read off the emitted binary, not off ``settings.pxi``: the flag is what
    the source SAYS and the symbol is what the compiler DID, and poly_106
    rev 1 is the case where the two disagreed. Also refuses a tree whose
    extensions are symlinks out to another checkout -- that is the mechanism
    by which a foreign, gitignored build flag wins silently.
    """
    tree = Path(tree)
    sos = sorted(tree.glob("rmgpy/solver/base*.so"))
    if not sos:
        return None, "no built base*.so"
    so = sos[0]
    if so.is_symlink():
        return None, f"SYMLINK to {os.readlink(so)}"
    out = subprocess.run(["strings", str(so)], capture_output=True,
                         text=True).stdout
    found = sorted({m for m in ("pydas.daspk", "pydas.dassl") if m in out})
    if len(found) != 1:
        return None, f"ambiguous backend symbols: {found}"
    return found[0], str(so)


def _junit_failures(xml_path):
    """The set of ``classname::name`` ids that failed or errored."""
    import xml.etree.ElementTree as ET
    out = set()
    for tc in ET.parse(str(xml_path)).iter("testcase"):
        if tc.find("failure") is not None or tc.find("error") is not None:
            out.add(f"{tc.get('classname')}::{tc.get('name')}")
    return out


def check(label, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        for line in str(detail).splitlines():
            print(f"         {line}")
    if not ok:
        FAILURES.append(label)
    return ok


def die(msg):
    print(f"\nVERIFIER CANNOT RUN: {msg}")
    sys.exit(2)


# ==========================================================================
# artifact discovery -- a REAL completed run's sidecar, never a fixture
# ==========================================================================
def find_pre_fix_sidecar():
    env = os.environ.get("I059_SIDECAR")
    cands = ([Path(env)] if env else
             sorted(Path(os.path.expanduser("~/runs/RMG")).glob(
                 "poly_*/chemkin/polymer_pools.json"),
                 key=lambda p: p.stat().st_mtime, reverse=True))
    for p in cands:
        if not p.is_file():
            continue
        try:
            art = json.loads(p.read_text())
        except Exception:
            continue
        pools = [x for x in art.get("pools", []) if isinstance(x, dict)]
        defect = [x for x in pools if DEFECT_KEY in x]
        # "pre-fix" means: carries spawned daughters WITH a defect and
        # WITHOUT the evidence block. That is the shape the warning is about.
        if defect and all(PROV_KEY not in x for x in defect):
            return p, art
    die("no real pre-fix polymer sidecar found (need one with defect-"
        "carrying pools and no %s block). Set I059_SIDECAR." % PROV_KEY)


# ==========================================================================
# TA arm: load a sidecar in the TA env, report warnings / exception
# ==========================================================================
TA_CHILD = r'''
import json, sys, warnings
from ta.mechanism import _load_sidecar
out = {"warnings": [], "error": None}
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    try:
        _load_sidecar(sys.argv[1])
    except Exception as exc:
        out["error"] = f"{type(exc).__name__}: {exc}"
    out["warnings"] = [str(x.message) for x in w]
print("@@@" + json.dumps(out))
'''


def ta_load(sidecar_path, ta_dir):
    """Load ``sidecar_path`` through the TA loader living in ``ta_dir``."""
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as fh:
        fh.write(TA_CHILD)
        child = fh.name
    try:
        cmd = (f"source {CONDA_SH} && conda activate {TA_ENV} && "
               f"cd {ta_dir} && python -W always {child} {sidecar_path}")
        res = subprocess.run(["bash", "-c", cmd], capture_output=True,
                             text=True)
        for line in res.stdout.splitlines():
            if line.startswith("@@@"):
                return json.loads(line[3:])
        die(f"TA arm produced no result.\nstdout:\n{res.stdout}\n"
            f"stderr:\n{res.stderr}")
    finally:
        os.unlink(child)


def census_warnings(result):
    return [w for w in result["warnings"] if CENSUS_MARKER in w]


# ==========================================================================
def main():
    sidecar_path, pre_art = find_pre_fix_sidecar()
    pools = [p for p in pre_art["pools"] if isinstance(p, dict)]
    by_label = {p["label"]: p for p in pools}
    defect_pools = [p for p in pools if DEFECT_KEY in p]
    print(f"artifact       : {sidecar_path}")
    print(f"rmg_commit     : {pre_art.get('rmg_commit')}")
    print(f"schema_version : {pre_art.get('schema_version')}")
    print(f"defect pools   : {[p['label'] for p in defect_pools]}")
    print(f"RMG tree       : {RMG_ROOT}")
    print(f"TA tree        : {TA_ROOT}")
    print(f"pinned bases   : RMG {RMG_BASE_SHA} / TA {TA_BASE_SHA}")
    print("check-3 'before': the defect value the PRE-FIX producer wrote "
          "into the artifact above")
    print(f"                  (its own rmg_commit "
          f"{pre_art.get('rmg_commit')}, a pre-{RMG_BASE_SHA} build)")
    print()

    tmp = Path(tempfile.mkdtemp(prefix="i059_"))
    ta_base_wt = tmp / "ta_base"
    try:
        # ------------------------------------------------------------------
        # CHECK 1 -- the warning reproduces at the PINNED PRE-FIX consumer
        # against the untouched real artifact.
        # ------------------------------------------------------------------
        r = subprocess.run(
            ["git", "-C", str(TA_ROOT), "worktree", "add", "--detach",
             str(ta_base_wt), TA_BASE_SHA],
            capture_output=True, text=True)
        if r.returncode != 0:
            die(f"could not check out pinned pre-fix TA {TA_BASE_SHA}:\n"
                f"{r.stderr}")
        pre = ta_load(sidecar_path, ta_base_wt)
        pre_census = census_warnings(pre)
        ok1 = check(
            "1. census warning REPRODUCES at pre-fix consumer "
            f"({TA_BASE_SHA}) on the real pre-fix artifact",
            pre["error"] is None and len(pre_census) == len(defect_pools),
            f"error={pre['error']}  census warnings="
            f"{len(pre_census)} (expected {len(defect_pools)})")
        print("\n--- reproduced warning, verbatim (1 of "
              f"{len(pre_census)}) ---")
        print(pre_census[0] if pre_census else "<none>")
        print("--- end ---\n")

        # ------------------------------------------------------------------
        # Build the POST-FIX artifact. The evidence block is taken from the
        # LIVE PRODUCER, not hand-written: each daughter is re-spawned
        # through the same Polymer path the run used, and the block the
        # serializer emits for it is spliced onto the real pool entry.
        # ------------------------------------------------------------------
        sys.path.insert(0, str(RMG_ROOT))
        from rmgpy.molecule import Molecule
        from rmgpy.polymer import Polymer, _serialize_pool_for_sidecar

        post_art = json.loads(json.dumps(pre_art))
        post_by_label = {p["label"]: p for p in post_art["pools"]
                         if isinstance(p, dict)}
        value_rows = []
        for entry in defect_pools:
            parent_label = entry.get("parent_pool")
            parent_entry = by_label.get(parent_label)
            if parent_entry is None:
                die(f"pool {entry['label']!r} names parent "
                    f"{parent_label!r}, absent from the artifact")
            mono = Molecule().from_adjacency_list(
                parent_entry["monomer_adj_list"])
            egs = [Molecule().from_smiles(s)
                   for s in parent_entry["end_groups"]]
            parent = Polymer(label=parent_label, monomer=mono, end_groups=egs,
                             cutoff=parent_entry["cutoff"],
                             Mn=parent_entry["mn_g_mol"],
                             Mw=parent_entry["mw_g_mol"],
                             initial_mass=parent_entry["initial_mass_g"])
            feat = Molecule().from_smiles(entry["feature_monomers_smiles"][0])
            daughter = parent._born_at_zero_mod_daughter(
                feat, source=entry["spawn_event_metadata"]["source"])
            fresh = _serialize_pool_for_sidecar(daughter)
            value_rows.append((entry["label"], entry[DEFECT_KEY],
                               fresh.get(DEFECT_KEY),
                               sorted(set(fresh) - set(entry))))
            post_by_label[entry["label"]][PROV_KEY] = fresh[PROV_KEY]

        # --------------------------------------------------------------
        # CHECK 3 -- the VALUE is unchanged, byte for byte. The "before" is
        # the number the PRE-FIX producer actually wrote into this artifact;
        # the "after" is what the post-fix producer emits for the same pool.
        # --------------------------------------------------------------
        same_value = all(
            json.dumps(old) == json.dumps(new) for _, old, new, _ in value_rows)
        added_keys = sorted({k for _, _, _, ks in value_rows for k in ks})
        ok3 = check(
            "3. serialized chain_mass_defect_g_mol is byte-identical "
            "before/after, and the ONLY key the producer adds is the "
            "evidence block",
            same_value and added_keys == [PROV_KEY],
            "\n".join(f"{lbl}: pre-fix artifact={json.dumps(old)} "
                      f"post-fix producer={json.dumps(new)} "
                      f"added_keys={ks}"
                      for lbl, old, new, ks in value_rows))

        # --------------------------------------------------------------
        # CHECK 2 -- the warning is gone AFTER, for the right reason.
        # --------------------------------------------------------------
        post_path = tmp / "post.json"
        post_path.write_text(json.dumps(post_art, indent=1))
        post = ta_load(post_path, TA_ROOT)
        ok2a = check(
            "2a. post-fix consumer + post-fix artifact: NO census warning",
            post["error"] is None and not census_warnings(post),
            f"error={post['error']}  warnings={post['warnings']}")

        # arm two: strip the evidence back out -> the warning must RETURN.
        stripped = json.loads(json.dumps(post_art))
        for p in stripped["pools"]:
            p.pop(PROV_KEY, None)
        strip_path = tmp / "stripped.json"
        strip_path.write_text(json.dumps(stripped, indent=1))
        st = ta_load(strip_path, TA_ROOT)
        ok2b = check(
            "2b. strip the evidence back out -> the census warning RETURNS "
            "(so 2a is evidence-gated, not a skipped check)",
            st["error"] is None
            and len(census_warnings(st)) == len(defect_pools),
            f"error={st['error']}  census warnings="
            f"{len(census_warnings(st))} (expected {len(defect_pools)})")

        # arms three..six: forged evidence must be REJECTED, not accepted.
        target = defect_pools[0]["label"]
        forgeries = {}

        f1 = json.loads(json.dumps(post_art))
        for p in f1["pools"]:
            if p["label"] == target:
                p[PROV_KEY]["shed_g_mol"] = p[PROV_KEY]["shed_g_mol"] * 2.0
        forgeries["2c. shed_g_mol doubled (arithmetic no longer closes)"] = f1

        f2 = json.loads(json.dumps(post_art))
        for p in f2["pools"]:
            if p["label"] == target:
                p[PROV_KEY]["source"] = "side_group_homolysis"
        forgeries["2d. source contradicts the pool's spawn_event_metadata"] = f2

        f3 = json.loads(json.dumps(post_art))
        for p in f3["pools"]:
            if p["label"] == target:
                p[PROV_KEY]["inherited_g_mol"] = 0.5
                p[PROV_KEY]["shed_g_mol"] = p[DEFECT_KEY] - 0.5
        forgeries["2e. lineage does not close (parent's defect != inherited)"] = f3

        f4 = json.loads(json.dumps(post_art))
        for p in f4["pools"]:
            if p["label"] == target:
                p[PROV_KEY]["extra"] = 1
        forgeries["2f. extra key smuggled into the evidence block"] = f4

        ok2rest = True
        for lbl, art in forgeries.items():
            fp = tmp / (lbl.split(".")[0] + ".json")
            fp.write_text(json.dumps(art, indent=1))
            res = ta_load(fp, TA_ROOT)
            got = res["error"] or ""
            ok2rest &= check(lbl + " -> REJECTED",
                             got.startswith("ValueError"), got or "<no error>")

        # 2h -- the COPY shape, which is the trap in this design. RMG's
        # Polymer.copy() carries the defect and its evidence but NOT
        # spawn_metadata, so a copied defect-bearing pool serializes with the
        # {'source': 'input'} absence sentinel beside a block naming the
        # channel that created the defect upstream. A cross-pin demanding
        # those two agree would hard-REJECT a real producer artifact -- this
        # arm exists because the first version of the check did exactly that.
        f6 = json.loads(json.dumps(post_art))
        for p in f6["pools"]:
            if p["label"] == target:
                p["spawn_event_metadata"] = {"source": "input"}
        fp = tmp / "copyshape.json"
        fp.write_text(json.dumps(f6, indent=1))
        res = ta_load(fp, TA_ROOT)
        ok2rest &= check(
            "2h. the Polymer.copy() shape (defect + evidence, no spawn "
            "event of its own) is ACCEPTED, not rejected",
            res["error"] is None and not census_warnings(res),
            f"error={res['error']}  census warnings="
            f"{len(census_warnings(res))}")

        # 2g -- the third state: the block's own arithmetic closes, but the
        # parent it names is not in the artifact, so the inherited-mass half
        # cannot be checked. Must NOT be accepted (that would grant silence
        # to half-evidence) and must NOT raise (a pruned parent is a real
        # artifact shape) -- it warns, and stays in the census.
        f5 = json.loads(json.dumps(post_art))
        for p in f5["pools"]:
            if p["label"] == target:
                p[PROV_KEY]["parent_pool"] = "no_such_pool_in_this_artifact"
        fp = tmp / "unresolved.json"
        fp.write_text(json.dumps(f5, indent=1))
        res = ta_load(fp, TA_ROOT)
        unres = [w for w in res["warnings"] if UNRESOLVED_MARKER in w]
        ok2rest &= check(
            "2g. evidence whose named parent is ABSENT -> not accepted, not "
            "raised: warned, and still counted in the census",
            res["error"] is None and len(unres) == 1
            and len(census_warnings(res)) == 1,
            f"error={res['error']}  unresolved warnings={len(unres)}  "
            f"census warnings={len(census_warnings(res))} (expect 1 and 1: "
            f"only the tampered pool loses its evidence)")

        # --------------------------------------------------------------
        # CHECK 5 -- old artifacts still load through the POST-FIX consumer,
        # with behaviour asserted, not assumed.
        # --------------------------------------------------------------
        old = ta_load(sidecar_path, TA_ROOT)
        old_census = census_warnings(old)
        ok5 = check(
            "5. PRE-fix artifact through the POST-fix consumer: loads, no "
            "exception, and still warns once per unevidenced defect pool",
            old["error"] is None
            and len(old_census) == len(defect_pools),
            f"error={old['error']}  census warnings={len(old_census)} "
            f"(expected {len(defect_pools)}); intended behaviour = WARN, "
            f"never fail: absence of evidence is legal, it just confirms "
            f"nothing")
        ok5 &= check(
            "5b. the pre-fix warning text is unchanged by this ticket "
            "(pre-fix consumer vs post-fix consumer, same artifact)",
            sorted(old_census) == sorted(pre_census),
            "identical" if sorted(old_census) == sorted(pre_census)
            else "TEXT DIFFERS -- the warning was modified")

        # --------------------------------------------------------------
        # CHECK 7 -- the OTHER direction of compatibility, and the one the
        # no-schema-bump decision actually rests on: an OLD consumer reading
        # a NEW artifact. Check 5 measured old-artifact/new-consumer; that is
        # the easy direction. This one was originally argued rather than
        # measured ("an older reader ignores the key and degrades to today's
        # warning"), and arguing it is not good enough for the claim that
        # holds up a compatibility decision.
        # --------------------------------------------------------------
        oldc = ta_load(post_path, ta_base_wt)
        oldc_census = census_warnings(oldc)
        ok7 = check(
            f"7a. NEW artifact through the OLD consumer ({TA_BASE_SHA}): "
            "loads, does not raise, and degrades to exactly today's census "
            "warning",
            oldc["error"] is None
            and len(oldc_census) == len(defect_pools),
            f"error={oldc['error']}  census warnings={len(oldc_census)} "
            f"(expected {len(defect_pools)})")
        # ...and the degradation is INERT, not merely survivable: the old
        # consumer must behave identically with and without the new key.
        ok7 &= check(
            "7b. the new key is INERT at the old consumer -- its warning set "
            "on the NEW artifact is identical to its warning set on the OLD "
            "one (ignored, not merely tolerated)",
            sorted(oldc["warnings"]) == sorted(pre["warnings"]),
            "identical" if sorted(oldc["warnings"]) == sorted(pre["warnings"])
            else ("DIFFERS\n old-artifact: "
                  f"{sorted(pre['warnings'])}\n new-artifact: "
                  f"{sorted(oldc['warnings'])}"))

        # 7c -- the boundary of 7a/7b, and it is NOT where the no-bump
        # rationale assumed. TA holds a composition-carrying (copolymer) pool
        # to a CLOSED key set (_COMPOSITION_CARRIER_POOL_KEYS) and raises on
        # any key outside it, with a message that says in as many words that
        # the producer must bump the schema minor. So on such a pool the new
        # key is fatal to an old consumer, not inert. Measured both ways.
        # The stamp must move to 3.2 with the block: TA elects the copolymer
        # rung on PRESENCE of a composition block and refuses an older stamp
        # that carries one, so a 2.5-stamped probe is rejected by the version
        # rule before it ever reaches the key allowlist. (Found by this check
        # failing for that reason -- the earlier guard is doing its job.)
        comp = json.loads(json.dumps(post_art))
        comp["schema_version"] = "3.2"
        for p in comp["pools"]:
            if p["label"] == target:
                p["composition"] = {}
        fp = tmp / "composition_carrier.json"
        fp.write_text(json.dumps(comp, indent=1))
        old_comp = ta_load(fp, ta_base_wt)
        comp_ctl = json.loads(json.dumps(comp))
        for p in comp_ctl["pools"]:
            if p["label"] == target:
                p.pop(PROV_KEY, None)
        fp2 = tmp / "composition_carrier_control.json"
        fp2.write_text(json.dumps(comp_ctl, indent=1))
        old_comp_ctl = ta_load(fp2, ta_base_wt)
        blames_key = PROV_KEY in (old_comp["error"] or "")
        ctl_blames_key = PROV_KEY in (old_comp_ctl["error"] or "")
        ok7 &= check(
            "7c. BOUNDARY: on a composition-carrying pool the old consumer "
            "REFUSES the new key outright (closed carrier key set) -- and "
            "the same artifact without the key is not refused for that "
            "reason, so the refusal is attributable to this ticket's key",
            blames_key and not ctl_blames_key,
            f"with key   : {old_comp['error']}\n"
            f"without key: {old_comp_ctl['error']}")
        # ...and the post-fix consumer does not refuse for that reason.
        new_comp = ta_load(fp, TA_ROOT)
        ok7 &= check(
            "7d. the post-fix consumer does NOT refuse the same pool over "
            "this key (it is in the carrier allowlist)",
            PROV_KEY not in (new_comp["error"] or ""),
            f"error={new_comp['error']}")

        # 7e -- so is 7c a live break or a future landmine? Measured against
        # every real artifact on disk: if no producer output has ever carried
        # a composition block, no existing consumer can hit 7c today.
        all_sidecars = sorted(Path(os.path.expanduser("~/runs/RMG")).glob(
            "poly_*/chemkin/polymer_pools.json"))
        with_comp = []
        for sc in all_sidecars:
            try:
                a = json.loads(sc.read_text())
            except Exception:
                continue
            if any(isinstance(p, dict) and "composition" in p
                   for p in a.get("pools", [])):
                with_comp.append(str(sc))
        ok7 &= check(
            "7e. 7c is a FUTURE landmine, not a live break: no real producer "
            "artifact on this box carries a composition block, so no pool "
            "can carry both it and this key today",
            not with_comp,
            f"scanned {len(all_sidecars)} real sidecars under ~/runs/RMG; "
            f"carrying a composition block: {with_comp or 'none'}")

        # --------------------------------------------------------------
        # CHECK 4 -- the bias, measured. The daughters hold zero moments in
        # every artifact today, so mu0*defect is exactly zero and the
        # correction contributes nothing. Construct the case where it does:
        # move the parent's whole chain population into the daughters at
        # constant (mu0, mu1) -- the producer's own intact-chain-transfer
        # contract -- and read off what the defect term then subtracts.
        # --------------------------------------------------------------
        zero_now = all(
            not any(float(m) for m in (p.get("moments") or [0, 0, 0]))
            for p in defect_pools)
        check("4a. the correction currently multiplies ZERO (every "
              "defect-carrying pool holds mu0 = mu1 = mu2 = 0 today)",
              zero_now,
              "; ".join(f"{p['label']}: moments={p.get('moments')}"
                        for p in defect_pools))

        parent_label = defect_pools[0]["parent_pool"]
        par = by_label[parent_label]
        mu0_p, mu1_p, mu2_p = [float(x) for x in par["moments"]]
        MW = float(par["monomer_mw_g_mol"])
        n = len(defect_pools)
        share = 1.0 / n
        naive = 0.0        # sum mu1*MW           (defect-blind)
        aware = 0.0        # sum mu1*MW - mu0*d   (NORMATIVE, defect-aware)
        for p in defect_pools:
            mu0, mu1 = mu0_p * share, mu1_p * share
            d = float(p[DEFECT_KEY])
            naive += mu1 * MW
            aware += mu1 * MW - mu0 * d
        bias_g = naive - aware
        dp = mu1_p / mu0_p if mu0_p else float("nan")
        d0 = float(defect_pools[0][DEFECT_KEY])
        print()
        print("  --- CHECK 4: measured mass-curve bias -----------------")
        print(f"  parent pool          : {parent_label}")
        print(f"  parent moments       : mu0={mu0_p!r} mu1={mu1_p!r} mol")
        print(f"  monomer_mw_g_mol     : {MW!r}")
        print(f"  per-chain defect     : {d0!r} g/mol")
        print(f"  chains moved         : all {mu0_p!r} mol, split {n} ways, "
              f"at constant (mu0, mu1)")
        print(f"  defect-blind  mass   : {naive!r} g")
        print(f"  defect-aware  mass   : {aware!r} g")
        print(f"  BIAS                 : {bias_g!r} g "
              f"= {100.0 * bias_g / naive:.6f} % of the condensed mass")
        print(f"  initial sample mass  : {par['initial_mass_g']!r} g "
              f"-> {100.0 * bias_g / float(par['initial_mass_g']):.6f} % "
              f"of the TGA basis")
        print(f"  number-average DP    : {dp:.4f}  (bias = (1/DP)*(d/MW))")
        print(f"  ceiling as DP -> 1   : "
              f"{100.0 * d0 / MW:.6f} % of the condensed mass")
        print("  -------------------------------------------------------")
        ok4 = check(
            "4b. the bias is real and non-zero once the daughters hold "
            "moments (this is the number that makes the ticket matter)",
            bias_g > 0.0, f"bias = {bias_g!r} g")

        # --------------------------------------------------------------
        # CHECK 6 -- suites.
        # --------------------------------------------------------------
        if os.environ.get("I059_SKIP_SUITES") == "1":
            # Escape hatch for iterating on checks 1-5. NOT a pass: it FAILS
            # the check, loudly, so a run that skipped the suites can never
            # be mistaken for one that ran them.
            print("\n  --- CHECK 6: test suites ------------------------------")
            check("6. RMG and TA suites both green", False,
                  "SKIPPED by I059_SKIP_SUITES=1 -- this run does NOT "
                  "verify the suites")
            return 1 if FAILURES else 0
        print("\n  --- CHECK 6: test suites ------------------------------")
        print("  The RMG unit suite is NOT green at the pinned base, so the "
              "bar is 'introduces no NEW failure', not")
        print("  'green' -- the latter would be a claim this ticket cannot "
              "make and did not cause. BOTH sides are now")
        print("  measured live: the before-set was previously a constant, and "
              "that constant was invalidated twice over")
        print("  (the base moved to d78c7211f, and it had been measured on "
              "DASSL binaries -- see RMG_BASE_FAILURES_STALE).")
        print("  Both trees are asserted to be on the SAME backend first, "
              "because the backend comes from a gitignored")
        print("  settings.pxi and a mismatch would surface as phantom new "
              "failures.")

        # Backend parity FIRST -- an unequal comparison is worse than none.
        base_tree = Path(os.environ.get(
            "I059_RMG_BASE_TREE",
            os.path.expanduser("~/Code/RMG-Py-i059-baseline")))
        ok6 = True
        for name, tree, want_sha in (("tree under test", RMG_ROOT, None),
                                     ("pinned base", base_tree,
                                      RMG_BASE_SHA)):
            bknd, why = _solver_backend(tree)
            head = subprocess.run(
                ["git", "-C", str(tree), "rev-parse", "HEAD"],
                capture_output=True, text=True).stdout.strip()
            sha_ok = want_sha is None or head.startswith(want_sha)
            ok6 &= check(
                f"6-pre. {name} ({tree.name}): own extensions, backend "
                f"{REQUIRED_BACKEND}"
                + (f", HEAD == {want_sha}" if want_sha else ""),
                bknd == REQUIRED_BACKEND and sha_ok,
                f"backend={bknd} ({why})  HEAD={head[:9]}")
        if not ok6:
            check("6. no NEW RMG failure vs the pinned base, and the TA "
                  "suite is green", False,
                  "NOT RUN -- the two trees are not comparable (see 6-pre). "
                  "Running them anyway would produce a number that looks "
                  "like evidence and is not.")
            return 1

        rmg_xml = tmp / "rmg_after.xml"
        base_xml = tmp / "rmg_base.xml"
        env_name = os.environ.get('CONDA_DEFAULT_ENV', 'rmg_env')

        def _suite(tree, xml):
            cmd = (f"source {CONDA_SH} && conda activate {env_name} && "
                   f"cd {tree} && python -m pytest "
                   f'-m "not functional and not database" '
                   f"-q --tb=no -p no:cacheprovider --junitxml={xml}")
            r = subprocess.run(["bash", "-c", cmd], capture_output=True,
                               text=True)
            tail = [x for x in r.stdout.strip().splitlines() if x][-1:]
            return r, (tail[0] if tail else "<no output>")

        base_run, base_tail = _suite(base_tree, base_xml)
        print(f"\n  RMG BASE  ({RMG_BASE_SHA}, {base_tree.name}): {base_tail}")
        before = _junit_failures(base_xml)
        rmg, rmg_tail = _suite(RMG_ROOT, rmg_xml)
        print(f"  RMG AFTER ({RMG_ROOT.name}): {rmg_tail}")
        after = _junit_failures(rmg_xml)
        new_fail = sorted(after - before)
        fixed = sorted(before - after)
        drift = sorted(before ^ RMG_BASE_FAILURES_STALE_9E7E0C4D5)
        print(f"  base failures: {len(before)}   after: {len(after)}")
        print(f"  (informational) base-set drift vs the stale 9e7e0c4d5/DASSL "
              f"constant: {len(drift)} entries differ")
        print(f"  failures now: {len(after)}   NEW vs {RMG_BASE_SHA}: "
              f"{len(new_fail)}   no longer failing: {len(fixed)}")
        for f in new_fail:
            print(f"    NEW FAILURE: {f}")
        for f in fixed:
            print(f"    no longer failing: {f}")

        ta_cmd = (f"source {CONDA_SH} && conda activate {TA_ENV} && "
                  f"cd {TA_ROOT} && python -m pytest -q --tb=no "
                  f"-p no:cacheprovider")
        ta = subprocess.run(["bash", "-c", ta_cmd], capture_output=True,
                            text=True)
        ta_tail = [x for x in ta.stdout.strip().splitlines() if x][-1:]
        print(f"  TA  (full suite, must be GREEN -- it is at "
              f"{TA_BASE_SHA}): "
              f"{ta_tail[0] if ta_tail else '<no output>'}")
        print("  -------------------------------------------------------")
        ok6 = check(
            "6. no NEW RMG failure vs the pinned pre-fix base, and the TA "
            "suite is green",
            not new_fail and ta.returncode == 0,
            f"new RMG failures={new_fail}  ta exit={ta.returncode}")

    finally:
        subprocess.run(["git", "-C", str(TA_ROOT), "worktree", "remove",
                        "--force", str(ta_base_wt)],
                       capture_output=True, text=True)
        shutil.rmtree(tmp, ignore_errors=True)

    print()
    if FAILURES:
        print(f"VERIFIER FAILED -- {len(FAILURES)} check(s):")
        for f in FAILURES:
            print(f"  - {f}")
        return 1
    print("VERIFIER PASSED -- all checks green.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
