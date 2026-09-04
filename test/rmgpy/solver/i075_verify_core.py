"""i075 verifier -- python half.

Recomputes, on the live build, the checks that need a running solver:

  3. a malformed / zero RMG_STALL_TRACE_EVERY cannot affect a run with the
     diagnostic OFF (and the pre-fix statement, taken verbatim out of git, is
     shown to raise on exactly those values);
  4. the open-descriptor count is BOUNDED across a sequence of simulations with
     the diagnostic ON (measured from /proc/self/fd, not argued);
  5. the i061 floor's own behaviour is unchanged by the i075 edits -- physical
     species untouched, the healthy trajectory still agrees to its recorded
     magnitude, the near-floor canary still BITWISE identical;
  7. the moment-slot scope guard fires on a slot that carries a real molar
     amount, and stays silent on a legitimate pool.

Every number is recomputed here. Nothing is restated from a report.
Exits 0 only when all of them pass.
"""
import json
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TREE = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, TREE)
sys.path.insert(0, HERE)

import numpy as np  # noqa: E402

FAILS = []
# The commit that introduced the stall trace. Pinned by sha, never by branch
# name, so this keeps working after the branch is merged and deleted.
DIAG_COMMIT = "ee73b7890c79354c549c742bb0208cce4d797940"
# The recorded magnitude of the healthy-trajectory difference from the i061
# closure ("unchanged to 1.6e-11 rel"). The bar below is the closure's own
# acceptance bar, 1e-8 (four decades under the deck rtol of 1e-4); the recorded
# magnitude is printed alongside so a drift is visible even while passing.
I061_RECORDED_MAX_REL = 1.58e-11
HEALTHY_REL_BAR = 1.0e-8


def check(name, ok, detail=""):
    print("%-4s %-42s %s" % ("PASS" if ok else "FAIL", name, detail))
    if not ok:
        FAILS.append(name)


def _run(args, env=None, cwd=None, timeout=1800):
    e = dict(os.environ)
    if env:
        e.update({k: v for k, v in env.items() if v is not None})
        for k, v in env.items():
            if v is None:
                e.pop(k, None)
    return subprocess.run(args, env=e, cwd=cwd or TREE, capture_output=True,
                          text=True, timeout=timeout)


print("tree under test : %s" % TREE)
print("python          : %s" % sys.executable)

# ===========================================================================
print("\n--- 3. a bad RMG_STALL_TRACE_EVERY cannot reach a diagnostic-OFF run")
# ===========================================================================
# (a) MEASURED on the pre-fix source. The pre-fix setup line is extracted
#     verbatim from the commit that introduced it and executed, so "it would
#     have failed" is a measurement of the historical code rather than a claim
#     about it. This does NOT re-run the pre-fix compiled extension -- see the
#     limitation printed below.
old_src = subprocess.run(
    ["git", "-C", TREE, "show", "%s:rmgpy/solver/base.pyx" % DIAG_COMMIT],
    capture_output=True, text=True)
m = re.search(r"^\s*(_stall_every = int\(os\.environ\.get\("
              r"'RMG_STALL_TRACE_EVERY', '2000'\)\))\s*$",
              old_src.stdout, re.M)
check("3a_prefix_line_recovered", bool(m),
      "pre-fix setup line recovered from %s: %r"
      % (DIAG_COMMIT[:9], m.group(1) if m else None))

if m:
    stmt = m.group(1)
    outcomes = {}
    for bad in ("banana", "", "0"):
        ns = {"os": os}
        os.environ["RMG_STALL_TRACE_EVERY"] = bad
        os.environ.pop("RMG_STALL_TRACE", None)   # diagnostic OFF
        try:
            exec(stmt, ns)
            # a parsed value still divides at `_stall_iter % _stall_every`
            try:
                1 % ns["_stall_every"]
                outcomes[bad] = "ok"
            except ZeroDivisionError as exc:
                outcomes[bad] = type(exc).__name__
        except Exception as exc:                  # noqa: BLE001
            outcomes[bad] = type(exc).__name__
        os.environ.pop("RMG_STALL_TRACE_EVERY", None)
    check("3b_prefix_would_have_failed",
          outcomes.get("banana") == "ValueError"
          and outcomes.get("") == "ValueError"
          and outcomes.get("0") == "ZeroDivisionError",
          "pre-fix outcomes with the diagnostic OFF: %r" % outcomes)

# (b) LIVE on the built extension: the same values must leave a diagnostic-off
#     run bit-for-bit where it was. Compared against a control run with the
#     variable absent entirely.
arm = os.path.join(HERE, "i075_ab_arm.py")
tmp = tempfile.mkdtemp(prefix="i075_verify.")
live = {}
for tag, val in (("control", None), ("malformed", "banana"),
                 ("zero", "0"), ("empty", "")):
    out = os.path.join(tmp, "gate_%s.json" % tag)
    r = _run([sys.executable, arm, "100.0", out],
             env={"RMG_STALL_TRACE_EVERY": val, "RMG_STALL_TRACE": None},
             cwd="/tmp")
    live[tag] = (r.returncode, out if r.returncode == 0 else r.stderr[-400:])
ok_rc = all(v[0] == 0 for v in live.values())
check("3c_live_runs_all_succeed", ok_rc,
      "returncodes: %r" % {k: v[0] for k, v in live.items()})
if ok_rc:
    ctl = json.load(open(live["control"][1]))
    same = {}
    for tag in ("malformed", "zero", "empty"):
        d = json.load(open(live[tag][1]))
        same[tag] = (d["healthy_traj"] == ctl["healthy_traj"]
                     and d["canary_y"] == ctl["canary_y"]
                     and d["atol_array"] == ctl["atol_array"])
    check("3d_live_diagnostic_off_unaffected", all(same.values()),
          "trajectory+canary+atol bitwise equal to the control run: %r" % same)
print("    LIMIT: 3b executes the pre-fix STATEMENT recovered from git, not a "
      "rebuilt pre-fix\n           extension. It establishes that the parse "
      "raised on those values where it\n           stood; it does not re-run "
      "the old binary.")

# ===========================================================================
print("\n--- 4. descriptor count is bounded with the diagnostic ON")
# ===========================================================================
FD_PROBE = r'''
import os, sys, glob
sys.path.insert(0, %(tree)r)
sys.path.insert(0, %(here)r)
import numpy as np
from rmgpy.reaction import Reaction
from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig
from rmgpy.solver.base import TerminationTime
from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
from solverPolymerTest import _two_pool_species, _KIN

TRACE = %(trace)r
os.environ['RMG_STALL_TRACE'] = TRACE
os.environ['RMG_STALL_TRACE_EVERY'] = '1'


def open_trace_fds():
    n = 0
    for fd in glob.glob('/proc/self/fd/*'):
        try:
            if os.readlink(fd) == TRACE:
                n += 1
        except OSError:
            pass
    return n


def one_simulation():
    sp, core, mask = _two_pool_species()
    feed = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
    feed.polymer_flux_archetype = 2
    fb = Reaction(reactants=[sp["B"]], products=[sp["B"], sp["G"]], **_KIN)
    mk = lambda l, mi: PolymerPoolConfig(
        label=l, xs=2, explicit_dp_to_species_index={}, mu_indices=mi,
        monomer_poly_index=None, k_scission=0.0, k_unzip=0.0,
        tail_kinetics=None)
    rs = HybridPolymerSystem(
        T=800.0, P=1.0e5, initial_mole_fractions={core[8]: 0.0}, V_poly=1.0,
        polymer_pools=[mk("A", (1, 2, 3)), mk("B", (5, 6, 7))],
        mass_transfer=[], gas_species_mask=mask.copy(),
        constant_gas_volume=False,
        initial_polymer_moments={"A": (1.0, 5.0, 30.0),
                                 "B": (0.8, 4.0, 24.0)},
        termination=[TerminationTime((1.0, "s"))],
        allow_default_prospective_edge=True, allow_unstamped_proxy_rows=True)
    ms = ModelSettings(tol_keep_in_edge=0.0, tol_move_to_core=1.0e8,
                       tol_interrupt_simulation=1.0e8)
    rs.simulate(list(core), [feed, fb], [], [], [], [],
                model_settings=ms, simulator_settings=SimulatorSettings())


counts = []
for _ in range(%(n)d):
    one_simulation()
    counts.append(open_trace_fds())
print('FDCOUNTS', counts, 'wrote=%%d' %% os.path.getsize(TRACE))
'''

trace_path = os.path.join(tmp, "stall_trace.log")
probe = FD_PROBE % {"tree": TREE, "here": HERE, "trace": trace_path, "n": 6}
probe_file = os.path.join(tmp, "fd_probe.py")
with open(probe_file, "w") as fh:
    fh.write(probe)
r = _run([sys.executable, probe_file], cwd="/tmp")
fm = re.search(r"FDCOUNTS (\[[^\]]*\]) wrote=(\d+)", r.stdout)
if not fm:
    check("4a_fd_probe_ran", False, "probe produced no FDCOUNTS: %s"
          % (r.stderr[-500:] or r.stdout[-500:]))
else:
    counts = json.loads(fm.group(1))
    wrote = int(fm.group(2))
    check("4a_fd_probe_ran", wrote > 0,
          "%d simulations with the diagnostic ON wrote %d bytes of trace"
          % (len(counts), wrote))
    check("4b_descriptors_bounded", max(counts) == 0,
          "open descriptors on the trace file after each simulation: %r "
          "(bounded: the handle never outlives one write)" % counts)

# The brief's premise for defect 2 was a growing descriptor count. Measure what
# the OLD shape actually did, so the fix is reported against a measured
# baseline rather than an assumed one.
LEAK_PROBE = r'''
import os, glob, sys
p = %(trace)r
def n_open():
    return sum(1 for fd in glob.glob('/proc/self/fd/*')
               if (os.path.islink(fd) and os.path.realpath(fd) == p))
def old_shape():
    fh = open(p, 'a', buffering=1)      # opened per call, never closed
    fh.write('x\n')
    return
counts = []
for _ in range(%(n)d):
    old_shape()
    counts.append(n_open())
print('OLDCOUNTS', counts)
'''
leak_file = os.path.join(tmp, "leak_probe.py")
with open(leak_file, "w") as fh:
    fh.write(LEAK_PROBE % {"trace": os.path.join(tmp, "old_shape.log"),
                           "n": 6})
r2 = _run([sys.executable, leak_file], cwd="/tmp")
om = re.search(r"OLDCOUNTS (\[[^\]]*\])", r2.stdout)
old_counts = json.loads(om.group(1)) if om else None
print("    measured, for the record: the pre-fix SHAPE (open per call, no "
      "close) leaves\n           %r open descriptors after 6 calls. CPython "
      "refcounting reclaims the\n           handle at frame exit, so the "
      "brief's 'leaks one descriptor per simulation'\n           is not what "
      "happens; the fix makes the lifetime EXPLICIT rather than\n           "
      "dependent on the interpreter's collection strategy." % (old_counts,))

# ===========================================================================
print("\n--- 5. the i061 floor's own behaviour is unchanged by the i075 edits")
# ===========================================================================
on_json = os.path.join(tmp, "arm_on.json")
off_json = os.path.join(tmp, "arm_off.json")
r_on = _run([sys.executable, arm, "100.0", on_json], cwd="/tmp")
r_off = _run([sys.executable, arm, "0.0", off_json], cwd="/tmp")
if r_on.returncode or r_off.returncode:
    check("5_arms_ran", False, "on rc=%d off rc=%d\n%s\n%s"
          % (r_on.returncode, r_off.returncode, r_on.stderr[-600:],
             r_off.stderr[-600:]))
else:
    on = json.load(open(on_json))
    off = json.load(open(off_json))
    check("5_arms_ran", True, "floor-ON and floor-OFF arms both completed")
    check("5a_off_arm_is_really_off", off["floored_slots"] == [],
          "K=0 floors nothing, so the OFF arm is the pre-fix array "
          "(module global is patchable: proven, not assumed)")

    inc = on["include_mask"]
    n_core = on["num_core_species"]
    floored = on["floored_slots"]
    phys = [i for i in floored if i < n_core and inc[i]]
    appended = [i for i in floored if i >= n_core]
    check("5b_physical_species_untouched",
          bool(floored) and not phys and not appended,
          "floored=%s ; physical floored=%s ; appended U/Z floored=%s"
          % (floored, phys, appended))
    check("5c_floored_are_moment_coords",
          all((i in on["mu_indices"]) or (not inc[i]) for i in floored),
          "every floored slot is a pool mu index or an is_moment_dummy slot "
          "(mu_indices=%s)" % on["mu_indices"])
    kept = [i for i in range(len(inc)) if i >= n_core or inc[i]]
    check("5d_atol_bitwise_equal_off_moment_slots",
          all(on["atol_array"][i] == off["atol_array"][i] for i in kept),
          "%d non-moment slots bitwise identical between the arms"
          % len(kept))

    maxrel, worst = 0.0, None
    for si, (ron, roff) in enumerate(zip(on["healthy_traj"],
                                         off["healthy_traj"])):
        for j, (x, y) in enumerate(zip(ron, roff)):
            fx, fy = float.fromhex(x), float.fromhex(y)
            if fx != 0.0:
                rel = abs(fx - fy) / abs(fx)
                if rel > maxrel:
                    maxrel, worst = rel, (si, j)
    check("5e_healthy_traj_unchanged", maxrel < HEALTHY_REL_BAR,
          "max rel diff = %.3e (i061 recorded %.2e; bar %.0e = 4 decades "
          "under the deck rtol 1e-4); worst at checkpoint/component %r"
          % (maxrel, I061_RECORDED_MAX_REL, HEALTHY_REL_BAR, worst))

    check("5f_canary_bit_identical",
          on["canary_y"] == off["canary_y"]
          and on["canary_dn"] == off["canary_dn"],
          "near-floor canary y and dn BITWISE identical across the arms "
          "(y[2]=%r)" % float.fromhex(on["canary_y2"]))
    check("5g_chem_consumers_bit_identical",
          on["chem_atol_array"] == off["chem_atol_array"]
          and on["pool_mu_floors"] == off["pool_mu_floors"]
          and on["softclamp_lam"] == off["softclamp_lam"]
          and on["jac_wt_atol"] == off["jac_wt_atol"],
          "_chem_atol_array, _pool_mu_floors, _softclamp_lam and "
          "_jac_wt_atol all bitwise identical")

# ===========================================================================
print("\n--- 7. the moment-slot scope guard")
# ===========================================================================
GUARD_PROBE = r'''
import sys, json
sys.path.insert(0, %(tree)r)
sys.path.insert(0, %(here)r)
import numpy as np
from rmgpy.reaction import Reaction
from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig
from solverPolymerTest import _two_pool_species, _spc, _KIN

MK = lambda l, mi: PolymerPoolConfig(
    label=l, xs=2, explicit_dp_to_species_index={}, mu_indices=mi,
    monomer_poly_index=None, k_scission=0.0, k_unzip=0.0, tail_kinetics=None)


def build(mutate):
    sp, core, mask = _two_pool_species()
    core = list(core)
    mask = list(mask)
    feed = Reaction(reactants=[sp["A"]], products=[sp["B"]], **_KIN)
    feed.polymer_flux_archetype = 2
    fb = Reaction(reactants=[sp["B"]], products=[sp["B"], sp["G"]], **_KIN)
    rxns = [feed, fb]
    imf = {core[8]: 0.0}
    mutate(sp, core, mask, rxns, imf)
    rs = HybridPolymerSystem(
        T=800.0, P=1.0e5, initial_mole_fractions=imf, V_poly=1.0,
        polymer_pools=[MK("A", (1, 2, 3)), MK("B", (5, 6, 7))],
        mass_transfer=[], gas_species_mask=np.array(mask, dtype=bool),
        constant_gas_volume=False,
        initial_polymer_moments={"A": (1.0, 5.0, 30.0),
                                 "B": (0.8, 4.0, 24.0)},
        termination=[], allow_default_prospective_edge=True,
        allow_unstamped_proxy_rows=True)
    rs.initialize_model(core, rxns, [], [], atol=1e-14, rtol=1e-4)
    return rs


def clean(sp, core, mask, rxns, imf):
    return None


# Every arm below drives the is_moment_dummy ARM of the mask complement, which
# is the one nothing upstream validates. The pool-mu_indices arm is already
# closed by validate_configuration (moment isolation: an index in reaction
# stoichiometry is refused there) and by _apply_pool_phase_overrides (mu slots
# are force-condensed), so it cannot be driven into these states at all -- see
# the note the driver prints.
def flag_in_network(sp, core, mask, rxns, imf):
    """G is a real, reactive gas species that carries the flag by mistake."""
    sp["G"].is_moment_dummy = True


def flag_deck_loading(sp, core, mask, rxns, imf):
    """A real species with a deck mole fraction that carries the flag."""
    extra = _spc("CCCCCCCC", "DILUENT")
    extra.is_moment_dummy = True
    core.append(extra)
    mask.append(False)          # condensed, and in no reaction
    imf[extra] = 0.0            # ... but the deck loads it


def flag_gas(sp, core, mask, rxns, imf):
    """A real GAS species, in no reaction and with no deck loading, that
    carries the flag."""
    extra = _spc("CCCCCCCCC", "INERT_GAS")
    extra.is_moment_dummy = True
    core.append(extra)
    mask.append(True)


out = {}
for name, mut in (("clean", clean), ("flag_in_network", flag_in_network),
                  ("flag_deck_loading", flag_deck_loading),
                  ("flag_gas", flag_gas)):
    try:
        build(mut)
        out[name] = None
    except Exception as e:
        out[name] = "%%s: %%s" %% (type(e).__name__, str(e))
print("GUARD" + json.dumps(out))
'''
guard_file = os.path.join(tmp, "guard_probe.py")
with open(guard_file, "w") as fh:
    fh.write(GUARD_PROBE % {"tree": TREE, "here": HERE})
r3 = _run([sys.executable, guard_file], cwd="/tmp")
gm = re.search(r"^GUARD(\{.*\})$", r3.stdout, re.M)
if not gm:
    check("7_guard_probe_ran", False,
          "probe produced no GUARD line: %s" % (r3.stderr[-800:] or r3.stdout[-500:]))
else:
    g = json.loads(gm.group(1))
    check("7a_guard_silent_on_a_legitimate_pool", g["clean"] is None,
          "an ordinary two-pool model initialises with no guard error")
    for arm_name, want in (("flag_in_network", "reaction network"),
                           ("flag_deck_loading", "initial MOLE FRACTION"),
                           ("flag_gas", "classified GAS")):
        msg = g[arm_name] or ""
        check("7b_guard_fires_%s" % arm_name,
              msg.startswith("ValueError:")
              and "moment error-weight floor" in msg and want in msg,
              (msg[:150] + "...") if msg else "NO ERROR RAISED")
    print("    The guard's three arms are exercised through the is_moment_dummy "
          "route, which\n           is the arm nothing upstream validates. The "
          "pool-mu_indices route CANNOT be\n           driven into these states: "
          "validate_configuration's moment-isolation check\n           refuses a "
          "mu index that appears in reaction stoichiometry, and\n           "
          "_apply_pool_phase_overrides force-condenses every mu slot before the "
          "gas\n           mask is read. That is a stronger result than the guard "
          "itself -- measured\n           in this session, not assumed.")

# ===========================================================================
print("\n" + "=" * 72)
if FAILS:
    print("i075 verifier core: FAILURES -> %s" % FAILS)
else:
    print("i075 verifier core: ALL CHECKS PASSED")
print("=" * 72)
sys.exit(1 if FAILS else 0)
