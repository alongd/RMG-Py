# I-195 — the reverse-rate wall on the irreversible argon ionisation

The 5 torr argon deck builds a clean mechanism (core `e-`, `Ar`, `Ar+`, one reaction) and dies at
final validation: `check_model()` computes a **reverse** collision-limit and reverse rate
coefficient for `Ar + e- => Arp + e- + e-`, which is declared `reversible = False`, and
`ElectronCollisionPlasma` has no reverse-rate branch, so `generate_reverse_rate_coefficient()`
raises `ReactionError`. This ticket **measured first**, took the fork honestly ((A) fails its bar),
and — after the manager reviewed the measurement and **explicitly lifted the check-weakening gate for
fix (C) only** — implemented **fix (C), graceful degradation**. The argon deck now completes.

Runtime: built from source (`python setup.py build_ext --inplace -j 8`, `BUILD_EXIT=0`). Database:
`RMG-database-i193-plasmaargon` (carries `PlasmaArgon`). Every run below imported this worktree's own
rmgpy, forced ahead of the shared `easy-install.pth` pin at the primary checkout via `PYTHONPATH`.

**Which commit produced which number.** Two distinct provenances, kept separate throughout:
- **Baseline / diagnosis (§ Q1–Q4, and the pre-fix wall)** — the *unfixed* branch tip `3b479a638`
  (Version 4.0.0). These are the measurements that establish the bug and the blast radius. The Q1
  instrumented run and the pre-fix argon run both report `git HEAD = 3b479a638`.
- **After-fix validation (§ implementation, run, verification)** — the *fixed* working tree
  (committed as the fix on top of the measurement commit). The RED/GREEN test, the `84 passed` suite,
  and the completing argon run were produced here, after rebuilding `reaction.py`. The post-fix argon
  run reports the fix commit's SHA, not `3b479a638`.

## The four questions, answered by measurement

### Q1 — is `self.reversible` False *at the point `check_model()` sees it*? — **YES**

Instrumented the real run: a `sitecustomize.py` wrapped `RMG.check_model` (pure Python) to log every
core reaction's flags at the exact point the collision-limit loop iterates them
(`main.py:1912`), before delegating to the original. Captured, verbatim, immediately before the crash:

```
Q1-PROBE core rxn='Ar(2) + e-(1) => Arp(3) + e-(1) + e-(1)' reversible=False kinetics=ElectronCollisionPlasma nreact=2 nprod=3
...
File "rmgpy/reaction.py", line 2016, in ...check_collision_limit_violation
File "rmgpy/reaction.py", line 1614, in ...generate_reverse_rate_coefficient
rmgpy.exceptions.ReactionError: Unexpected kinetics type <class 'rmgpy.kinetics.arrhenius.ElectronCollisionPlasma'>
```

The runtime core reaction is `reversible=False`, 2 reactants, 3 products — a library flag confirmed
on the runtime object, not merely in the file. The reverse branch (`reaction.py:2010`, guarded by
`len(self.products) >= 2`, **not** `self.reversible`) is entered because it has 3 products.

### Q2 — is `generate_reverse_rate_coefficient()` meaningless for an irreversible reaction? — **NO. It is a thermodynamic quantity, `k_f/K_eq`, independent of the `reversible` flag.**

`generate_reverse_rate_coefficient()` computes `k_r(T) = k_f(T) / K_eq(T)`, where `K_eq` comes from
species thermochemistry (`reaction.py:1528,1557,...`). It never reads `self.reversible`. Proven on
`reactionTest`'s `reaction2` (acetyl + O2 => acetylperoxy, Arrhenius, full thermo), flipping the flag:

```
reversible=True   kr(1000K)=8.318728e+07   kf/Keq=8.416289e+07   type=Arrhenius
reversible=False  kr(1000K)=8.318728e+07   kf/Keq=8.416289e+07   type=Arrhenius
```

Byte-identical. The `reversible` flag is a **modelling directive** (does the solver integrate the
reverse reaction?), not a physical claim that `k_r = 0`. Microscopic reversibility guarantees
`k_r = k_f/K_eq > 0` at any finite T. A **reverse** collision-limit violation therefore flags that
the specified **forward** rate is thermodynamically inconsistent — it implies an unphysically fast
reverse — which is actionable data-quality information **regardless of the `reversible` flag**.
The reverse check is not meaningless for irreversible reactions.

### Q3 — how many irreversible reactions reach the reverse branch? — **1964** (blast radius)

Measured over **69 kinetics libraries / 20 599 reactions** in `RMG-database-i193-plasmaargon`
(every `libraries/*/reactions.py`, loaded through the real `KineticsDatabase.load_libraries` path):

| bucket | count |
|---|---|
| reactions with `len(products) >= 2` (enter reverse branch) | 17 586 |
| of those, `reversible = False` (**blast radius**) | **1964** |
| ├─ reverse-**supported** kinetics — silently reverse-checked today, no crash | 1941 |
| └─ reverse-**unsupported** kinetics — would crash the run | 23 |
| &nbsp;&nbsp;&nbsp;&nbsp;`TwoTemperaturePlasma` (PlasmaAir) | 19 |
| &nbsp;&nbsp;&nbsp;&nbsp;`ElectronCollisionPlasma` (PlasmaAir 3 + PlasmaArgon 1) | 4 |

So **1941 ordinary irreversible reactions already receive a silent reverse collision-limit check
today** (exactly as the contract predicted for irreversible Arrhenius with ≥2 products), and Q2
proves that check is legitimate. Only **23** plasma-rate-law reactions crash, all confined to
kinetics types with no reverse implementation. (Top contributors to the 1964: CH3Cl 1156,
CurranPentane 369, CF2BrCl 134, Narayanaswamy 93, Chernov 45.)

### Q4 — is the forward check affected by the change under consideration? — **NO**

The forward branch (`reaction.py:2003-2009`, `2017-2022`) uses `calculate_coll_limit(reverse=False)`
and `get_rate_coefficient` — it never calls `generate_reverse_rate_coefficient`. Only the reverse
branch (`2010-2016`, crash at `2016`) does. Any fix scoped to the reverse call leaves the forward
check byte-for-byte unchanged.

## The verdict — (A) does NOT clear its bar

**(A) Guard the reverse branch with `self.reversible`.** Its bar: *show the suppressed reverse
checks on irreversible reactions are meaningless.* **They are not** (Q2): the reverse rate is a
flag-independent thermodynamic quantity, and a reverse violation flags a thermodynamically
inconsistent forward rate. (A) would suppress **1941 legitimate, currently-running** reverse checks
across ordinary combustion chemistry to dodge **23** plasma crashes. On this campaign the burden of
proof is on weakening a refusal; (A) cannot carry it. **(A) fails its bar, and is not implemented.**

**(B) Give `ElectronCollisionPlasma` reverse-rate support** — **GATED to the owner, forbidden, and
not attempted.** It is also the physically hard fix: the forward rate is evaluated at the electron
temperature `Te` while `K_eq` is a gas-temperature quantity, so the reverse of a two-temperature
plasma reaction is *not* simply `k_f/K_eq(T_gas)`. Detailed balance across two temperatures is a
research question — which is exactly why it is the owner's.

**(C) Degrade gracefully** — the collision-limit check should still run, but when a kinetics type has
no reverse implementation it should **warn and skip that one reaction's reverse direction** rather
than raise and abort the whole `check_model()`. This is strictly better than (A): it preserves all
1941 legitimate reverse checks and only skips the 23 genuinely un-computable ones (whose reverse rate
cannot be formed without the gated physics of B). It converts an uncaught crash — which today
produces *no* verdict for *any* reaction once it fires mid-loop — into a completed check plus a
logged warning. **Approved by the manager and implemented.**

## The manager ruling and why (C)'s gate does not bind

The manager accepted the refutation of (A), verified independently that
`generate_reverse_rate_coefficient` never references `self.reversible`, and **lifted the
check-weakening gate for (C) specifically**, recording why: the crash fires *mid-loop* inside
`check_model()`, so today a mechanism containing one plasma reaction gets **zero** collision-limit
verdicts for **any** of its reactions. (C) therefore makes the check fire **more** often — 1941
completed verdicts plus a warning, against 0 today — the opposite of the suppression the gate was
written to stop. The two-temperature detailed-balance question (fix B) is filed to the owner as the
real path; it was **not** attempted here.

## Implementation of fix (C) — dedicated-exception form (rework D-035)

The first form used a capability predicate that inspected only the **outer** kinetics object. But
`generate_reverse_rate_coefficient` **recurses**: its `MultiArrhenius` / `PDepArrhenius` /
`MultiPDepArrhenius` branches build a temporary `Reaction` per child and call themselves
(`reaction.py:~1591,1600,1609`). A `MultiArrhenius([ElectronCollisionPlasma])` reports a *supported*
outer type, is not skipped, and the recursive child call re-raises the exact `ReactionError` that had
been removed — the mid-loop crash returns. This shape exists in the repo already
(`plasmaExportTest.py:724`, `i167CanteraExportPathTest.py:441,459`:
`MultiArrhenius(arrhenius=[VoronovEIArrhenius(...)])`). The predicate approach is therefore wrong in
reach and is replaced.

All in `rmgpy/reaction.py` (+ one class in `rmgpy/exceptions.py`); cythonized module rebuilt,
`BUILD_EXIT=0`, before any test result was trusted:

- **A dedicated exception** `UnsupportedReverseRateError(ReactionError)` (in `rmgpy/exceptions.py`),
  raised **only** from `generate_reverse_rate_coefficient`'s final unsupported branch. The isinstance
  chain above that branch is the single source of truth for what is supported; this branch is its
  complement, so the two cannot drift, and the hand-maintained `_REVERSE_RATE_SUPPORTED_KINETICS`
  tuple and predicate are **deleted**.
- **Caught narrowly, handles recursion for free.** `check_collision_limit_violation` catches **only**
  `UnsupportedReverseRateError` around the reverse rate computation. Because the wrapper branches
  recurse, a child raising it propagates unchanged to the outermost caller — so the wrapped case is
  covered without any special handling. The dedicated subclass **honours the earlier objection to a
  broad catch**: it cannot be confused with the generic `ReactionError` that `get_equilibrium_constant`
  raises (`'Got equilibrium constant of 0'`), so a real defect on a supported reaction is never
  masked. The instinct against `try/except ReactionError` was right; the answer was one subclass away.
- **Loud, honest skip** (requirements 1, 2). One `logging.warning` names the reaction, names the outer
  kinetics class, adds "(or a component of it)" for the wrapped case, and says the reverse check was
  skipped and why. It **no longer claims the forward direction was checked** — that claim was false
  whenever `len(reactants) < 2` or `calculate_coll_limit(reverse=False)` raised `ValueError`. The
  forward block is byte-for-byte unchanged (Q4).

**Named limitation, accepted knowingly — the real breadth.** The rule, not a list, is the definition:
the reverse-direction collision-limit check is skipped for **any** kinetics type the
`generate_reverse_rate_coefficient` isinstance chain does not name — i.e. anything that reaches its
final `else`. It is **not only the plasma rate laws**. Verified against the chain at this commit
(`issubclass` of every `rmgpy.kinetics` class against the handled bases), the else-bound set is
`ArrheniusBM`, `ArrheniusEP`, `ArrheniusChargeTransferBM`, `Marcus`, `BadnellRRArrhenius`,
`VoronovEIArrhenius`, `PDepKineticsData`, `SurfaceArrheniusBEP`, `StickingCoefficientBEP`,
`SurfaceChargeTransferBEP`, `ElectronCollisionPlasma`, `TwoTemperaturePlasma` (plus the abstract base
`PDepKineticsModel`) — and anything added later without a reverse branch. This list is an illustration
verified at this commit, not the definition; treat the rule above as authoritative. For such a
reaction the **reverse-direction** check does not run; a forward check still runs when it applies. And
for some of these types, reaching `check_model` at all may itself signal an upstream pipeline defect (a
rate law that should have been converted or handled before final validation) — this change **converts
that hard failure into a warning**. The campaign accepts this knowingly; it is written here so it is
not discovered later. The real fix is (B): a physically-correct reverse (`k_f` at `Te`, `K_eq` at
`T_gas`), filed to the owner.

**Observed, not touched (pre-existing, outside this ticket; filed separately by the manager):**
1. `check_collision_limit_violation`'s `except ValueError: continue` (reverse `calculate_coll_limit`)
   drops a condition while the violator report still indexes `conditions[i]`, so a reported T/P can be
   silently wrong; a *forward* `ValueError` also `continue`s past that condition's reverse check.
   Both predate this branch and are untouched here.
2. `Reaction.get_rate_coefficient` special-cases `ArrheniusChargeTransfer` because its second argument
   is electrode potential, not pressure (`reaction.py:~1109`); the reverse loop calls the generated
   kinetics as `get_rate_coefficient(T, P)`, passing pascals where volts are expected. Pre-existing,
   untouched.

## The argon re-run (after fix C)

`python rmg.py input.py`, true exit captured from the interpreter itself:

```
PYTHON_EXIT=0
```

- **`MODEL GENERATION COMPLETED`: 1 occurrence** — the run finishes cleanly.
- The loud skip fired verbatim (honest text — no forward claim; the offending type is surfaced from
  the propagated exception, which for a wrapper names the child):
  `Warning: Reaction Ar(2) + e-(1) => Arp(3) + e-(1) + e-(1): skipping the reverse-direction
  collision-limit check. Its kinetics type ElectronCollisionPlasma (or a component of it) has no
  reverse-rate implementation in generate_reverse_rate_coefficient() (Reverse rate coefficient
  generation is not implemented for kinetics type <class
  'rmgpy.kinetics.arrhenius.ElectronCollisionPlasma'>), so no reverse rate coefficient could be formed
  to compare against the collision limit. The reverse direction is left unchecked for this reaction.`
- `No collision rate violators found in the model's core.` (forward check ran; no violators log written.)
- Final core: **3 species / 1 reaction** (`e-(1)`, `Ar(2)`, `Arp(3)`; `Ar + e- => Arp + e- + e-`).
- Final edge: **2 species / 2 reactions** (`[Li](4)`, `[Lip](5)` — inert lithium, not ours).
- `chemkin/chem.inp`: `SPECIES {e-(1), Ar(2), Arp(3)}`; the one reaction with its `TDEP/e-(1)/`
  annotation.

No further wall is reached — the completion line prints and the process exits 0. The argon campaign
run is complete.

### For the record: the pre-fix behaviour (I-194 wall)

Before the fix (same deck, stale `.so`): `PYTHON_EXIT=1`, **0×** `MODEL GENERATION COMPLETED`,
crash at `reaction.py:2016→1614` — identical to the I-194 record. That is the wall this ticket removed.

## Verification

| # | check | result |
|---|---|---|
| 1 | Q1–Q4 answered with measurements | **PASS** — above, each with file:line or command+output |
| 2 | verdict stated with blast-radius number | **PASS** — (A) fails; blast radius 1964 (1941 legit / 23 crash) |
| 3 | fix + RED/GREEN test pinning direct AND wrapped cases | **PASS** — `TestReverseCollisionLimitDegradation`: RED (predicate `.so`) `2 failed, 1 passed` — the wrapped `MultiArrhenius([ElectronCollisionPlasma])` raised `ReactionError: Unexpected kinetics type ElectronCollisionPlasma` *through the recursion*; GREEN (rebuilt) `3 passed` |
| 4 | argon re-run reported exactly | **PASS** — `PYTHON_EXIT=0`, `MODEL GENERATION COMPLETED` ×1, reworked skip warning fired, core 3/1, edge 2/2, chem.inp above |
| 5 | `pytest test/rmgpy/reactionTest.py` green, baseline + after | **PASS** — baseline `81 passed, 2 skipped`; after `84 passed, 2 skipped` (+3 new) |
| 6 | one suite per process | **PASS** — reactionTest run alone |

Rework (D-035) requirements: **HIGH (recursion)** — **FIXED** via dedicated
`UnsupportedReverseRateError`, caught narrowly; wrapped-child case now warns and continues (pinned by
`test_wrapped_unsupported_child_skipped_loudly_not_raised`, RED→GREEN). **(1) warning no longer lies**
— **DONE** (reports only the reverse skip, asserts nothing about the forward direction). **(2) real
breadth stated** — **DONE** (every unlisted ≥2-product type, plus the pipeline-defect-to-warning note).
**Observed-not-touched** (a: `ValueError`/condition-index drift; b: `ArrheniusChargeTransfer`
potential-vs-pressure) — **recorded, untouched**.

## Durable findings

- `check_collision_limit_violation` checks **both directions** by design (`main.py:1910` comment),
  gated on reactant/product count, **never** on `reversible`. The reverse rate it forms is
  `k_f/K_eq`, a thermodynamic quantity independent of the `reversible` flag (Q2). A reverse
  violation on an *irreversible* reaction is a real signal about the *forward* rate's
  thermodynamic consistency — not noise.
- 1964 irreversible ≥2-product reactions across the database enter the reverse branch; 1941 are
  silently, legitimately reverse-checked today; only 23 (all `TwoTemperaturePlasma` /
  `ElectronCollisionPlasma`, in PlasmaAir + PlasmaArgon) crash for lack of a reverse implementation.
- The physically-correct reverse for a two-temperature plasma rate is **not** `k_f/K_eq(T_gas)`
  (forward at `Te`, `K_eq` at `T_gas`); this is the owner's gated (B) work.
- Shared-env gotcha: `easy-install.pth` pins `import rmgpy` at the **primary** checkout
  `/home/alon/Code/RMG-Py`. A worktree run must lead `PYTHONPATH` with its own path, or an early
  `sitecustomize` import will cache the primary's rmgpy before the script dir is on `sys.path`.
