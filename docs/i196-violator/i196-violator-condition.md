# I-196 — the collision-limit violator report can name a condition the rate was never evaluated at

`Reaction.check_collision_limit_violation` (rmgpy/reaction.py) reported each violation with a
temperature/pressure string recovered by index (`conditions[i]`) while enumerating the *shortened*
`kf_list` / `kr_list`. `conditions` is never shortened when a condition is skipped, so a real
ratio could be attributed to a condition it did not come from. This ticket fixes the attribution by
carrying the retained condition **alongside** its value, and — deliberately — decouples the forward
failure from the reverse block so a computable reverse check is no longer suppressed when the
forward one cannot be computed.

Base: `i196-violator-condition` rebased onto `plasma` @ `b12ef4c43` (`git rev-parse HEAD^` ==
`b12ef4c43bd0c62f6267e34a76d49ee39ee1bfa6`), zero commits behind `plasma`. `b12ef4c43` is a
docs-only superset of the original `4cc872eb4` (it adds only `docs/i194-ar5torr-plasma-lineage/`, the
argon deck and its record — `git diff --stat 4cc872eb4 b12ef4c43` touches two docs files, no code),
so every measurement below is identical on either and the reaction.py `.so` did not need rebuilding
across the rebase. Built from source (`python setup.py build_ext --inplace -j 8`, `BUILD_EXIT=0`)
after editing the cythonized `reaction.py`. Database: `RMG-database-i193-plasmaargon` @
`ae0474a949de5c6c44796565d44ea6cebd5ca95a` (carries `PlasmaArgon`; SHA as reported by the deck run).
Runs imported this worktree's own rmgpy via `PYTHONPATH`. `/proc/loadavg` before the timed deck run:
`3.97 3.55 2.87`.

## The three questions, answered by measurement

### Q1 — what raises `ValueError` here, and how often?

The only `ValueError` source in the loop is `get_mean_sigma_and_epsilon`, at **reaction.py:2096**:

```python
if any([x == 0 for x in sigmas + epsilons]):
    raise ValueError
```

It fires when any participating species' Lennard-Jones `sigma` or `epsilon` is exactly 0 — i.e. a
species with missing/degenerate transport data (`get_transport_data()` → `generate_transport_data()`
returned zeros; in the plasma lineage an electron or a bare ion has no meaningful LJ transport).
`calculate_coll_limit` (reaction.py:2049) reaches it via `self.get_mean_sigma_and_epsilon(reverse)`
at line 2055.

**Crucially, this trigger is condition-INDEPENDENT.** `get_mean_sigma_and_epsilon` reads only
species transport data; `calculate_coll_limit` uses `temp` only in arithmetic that cannot raise
(`math.sqrt`/`math.exp` of positive quantities) and never touches pressure. So for a given reaction
it raises for **every** condition or **none** — it can never skip the first condition and keep a
later one. Consequently the two skip lists are all-or-nothing (`len(kf_list) ∈ {0, len(conditions)}`,
same for `kr_list`), and the *natural* transport `ValueError` **cannot** partially desynchronise a
list in the current code. Reachable in ordinary plasma use (a zero-transport species is common), but
only ever as a whole-reaction skip, not a per-condition one.

The literal contract scenario ("raises on the first condition, succeeds on a later one") therefore
requires a *condition-dependent* skip, which the natural trigger is not. The index-based attribution
is nonetheless wrong the instant any condition is skipped; the regression test induces a
condition-dependent skip (below) to expose it, and it is exactly this future fragility the fix
removes.

### Q2 — is the forward `continue` suppressing the reverse block intended? **No — a second defect. Decoupled.**

The old forward branch `continue`d on `ValueError`, which skipped the *rest of the loop body*,
including the reverse block. But the reverse collision limit depends on the **products'** transport
data, wholly independent of the reactants'. A reaction whose reactants lack transport data but whose
products have it had its perfectly computable reverse check silently dropped — and, unlike Q1's
attribution bug, this **is** reachable through the natural trigger. I decoupled them: the forward
`ValueError` now skips only the forward direction and falls through to the reverse block. This makes
the reverse check fire **more** often, never less — in scope per the ticket, and pinned by a test.

### Q3 — does any caller depend on the list lengths? **No.**

The sole consumer is `rmgpy/rmg/main.py:1916-1918`:

```python
violator_list = rxn.check_collision_limit_violation(t_min=..., t_max=..., p_min=..., p_max=...)
if violator_list:
    violators.extend(violator_list)
```

It consumes the flat list of `[reaction, direction, ratio, condition]` entries. The return shape is
unchanged (same 4-element entries, same list-of-lists). The internal `kf_list`/`conditions` lengths
never leaked out; nothing downstream depends on them.

## The fix

`rmgpy/reaction.py`, `check_collision_limit_violation`:

- Add parallel `conditions_f` / `conditions_r`, appended in lockstep with
  `kf_list`/`collision_limit_f` and `kr_list`/`collision_limit_r`. The report loops now format
  `conditions_f[i]` / `conditions_r[i]` instead of `conditions[i]`, so each violation carries the
  condition its ratio was computed at.
- Forward `except ValueError: continue` → `pass`, with the value appended in the `else` branch, so
  a forward failure no longer suppresses the reverse block (Q2 decoupling).

The reverse block's own `ValueError` `continue` and the `UnsupportedReverseRateError` handling
(landed by I-195) are untouched; both remain condition-independent and cannot leave a partially
populated `kr_list` (they latch/skip on the first reverse iteration).

## Verifier evidence

**1. Attribution test — RED before / GREEN after** (`test/rmgpy/reactionTest.py::TestCollisionLimitAttribution`).
Induces a condition-dependent skip via a `Reaction` subclass overriding the `cpdef`
`calculate_coll_limit` (dispatched through the extension-type vtable, verified). The violation's
ratio is computed at the surviving condition `t_max=2000 K`; the report must say so.

RED (stale pre-fix `.so`):
```
E       AssertionError: 300.0 K, 1.0 bar
E       assert '300.0 K, 1.0 bar' == '2000.0 K, 1.0 bar'
FAILED ...::test_forward_violation_names_the_condition_its_ratio_came_from
FAILED ...::test_reverse_still_checked_when_forward_uncomputable   (AssertionError: [])
1 passed  (test_no_skip_ratios_and_conditions_unchanged — correct in old code too)
```
GREEN (after rebuild): `6 passed, 83 deselected` (3 new + the 3 I-195 reverse-degradation tests).

**2. Forward + reverse ratios unchanged for the no-skip case.**
`test_no_skip_ratios_and_conditions_unchanged` asserts, at both `300 K` and `2000 K`, that the
reported forward ratio equals `get_rate_coefficient(T,P) / calculate_coll_limit(T, False)` and the
reverse ratio equals `reverse_kin.get_rate_coefficient(T,P) / calculate_coll_limit(T, True)`, to
1e-9, and that the condition strings map identity. GREEN — and this test passes on the pre-fix `.so`
too (no skip ⇒ nothing to misattribute), confirming the fix is a no-op when no condition is skipped.

**3. No regression — `test/rmgpy/reactionTest.py`.**
- Baseline (pre-existing tests, new class deselected): `84 passed, 2 skipped`.
- After (full file): `87 passed, 2 skipped` (= 84 + 3 new). One suite per process.

**4. 5 torr argon deck.** Ran the canonical deck `docs/i194-ar5torr-plasma-lineage/input.py`
(already on `plasma`, the single source of truth for this configuration) **unmodified** — no second
copy is kept in this ticket, so there is nothing to drift. It was executed from a scratch run
directory holding only a generated `rmgrc` (`database.directory` → `RMG-database-i193-plasmaargon/input`)
and RMG's own output artifacts; `PYTHONPATH=<this worktree> python rmg.py <canonical input.py>`, exit
captured from the interpreter itself:
```
Warning: Reaction Ar(2) + e-(1) => Arp(3) + e-(1) + e-(1): skipping the reverse-direction
collision-limit check. Its kinetics type ElectronCollisionPlasma ... left unchecked for this reaction.
No collision rate violators found in the model's core.
MODEL GENERATION COMPLETED
The final model core has 3 species and 1 reactions
PYTHON_EXIT=0
```
Final core `e-(1)`, `Ar(2)`, `Arp(3)`; `Ar + e- => Arp + e- + e-`. The `[Li]`/`[Lip]` edge species
are the inert lithium-library pair noted in the ticket, not ours. My changed function ran (the skip
warning and the "no violators" line are emitted by it) and completed cleanly.

## Named limitations

- **The decoupling can surface product-side transport-generation failures the old forward `continue`
  was masking.** When a reaction has >= 2 reactants whose forward collision limit is uncomputable,
  the old code `continue`d past the reverse block entirely; the reverse block never called
  `calculate_coll_limit(reverse=True)` (hence never `get_transport_data()` on the products) for that
  reaction. It now does. In a normal RMG run the transport database is fully loaded before
  `check_model()` runs, so this is strictly *increased checking* (a computable reverse limit that was
  being dropped now gets computed), not a new failure mode. But it is a real behaviour change on the
  reverse side: a product whose transport data must be *generated* (rather than looked up) is now
  exercised in this path where it was not before, so any latent failure in that generation would now
  be reached here rather than stay masked. Stated for completeness; no such failure was observed on
  the argon deck or the suite.
- **A caller passing a negative temperature directly can still get a condition-dependent `ValueError`
  via `math.sqrt`** in `calculate_coll_limit` (the `sqrt(8 * pi * kB * T * ...)` term). This is not a
  normal RMG path — `conditions` is built from `t_min`/`t_max`, both positive — and it is not this
  ticket's concern (not the transport `ValueError` of Q1). Recorded so it is not re-derived: it is
  the one way the skip predicate could become genuinely condition-dependent through the natural code,
  and the attribution fix makes even that case report the correct surviving condition.

## Non-goals honoured

- The electric-potential defect on the generated reverse kinetics object's `.get_rate_coefficient`
  (gated to owner) was left exactly as found.
- No change to `rmgpy/kinetics/`, `rmgpy/molecule/`, `rmgpy/data/`, or RMG-database.
- `check_collision_limit_violation` still checks both directions — the reverse now fires in strictly
  more cases (forward-uncomputable reactions), never fewer.
