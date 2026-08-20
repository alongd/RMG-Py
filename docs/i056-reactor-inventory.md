# I-056 — Archived Plasma Reactor: Port Inventory

**Report-only probe.** No `plasma.pyx` was created; nothing under `rmgpy/` was edited. This
document inventories what the archived reactor on branch `99` would require to re-land on today's
tree, so the eventual port ticket starts from a map instead of a blind `git cherry-pick`.

- Comparison basis: branch `99` vs today's HEAD `7e4fc58ed` (branch `i056-reactor-probe`).
  merge-base is `0e4889230`, but `99` is the *original full plasma branch* and the campaign has
  re-applied its other pieces (charge-transfer Arrhenius, plasma atomtypes, etc.) ticket by
  ticket. So `diff merge-base..99` is mostly already-landed work; the semantically correct diff
  for *the reactor* is `git diff HEAD 99 -- rmgpy/solver/`.
- Archived artifacts read in full: `99:rmgpy/solver/plasma.pyx` (1214 lines),
  `99:test/rmgpy/solver/plasmaTest.py` (896 lines), the `base.pyx`/`base.pxd`/`__init__.py`/
  `setup.py` deltas.
- Scope boundary (from the ticket): the reactor solver and its solver-infra hooks **only**. The
  input DSL builder, reactor-network integration, mechanism generation, benchmarking, and export
  are **out of scope** and are noted where they bound the reactor's contract, not inventoried.
  **But note (Q3/Q4): part of the electron machinery `plasma.pyx` depends on lives in the RMG
  driver `rmgpy/rmg/main.py`, not the input DSL — that is a hard port dependency, not sugar.**
- This inventory was hardened by an adversarial Codex review (persisted at
  `~/agents/adversarial/alongd-RMG-Py-i056-reactor/`, round 1). Every finding below marked
  *(adv.)* was raised there and then re-verified against source before inclusion.

---

## Category 1 — New files, re-appliable essentially as-is

| Artifact | On 99 | Notes |
|---|---|---|
| `rmgpy/solver/plasma.pyx` (1214 lines) | new | The reactor. Absent from HEAD (`git cat-file -e HEAD:rmgpy/solver/plasma.pyx` → absent). Structurally a near-verbatim clone of `simple.pyx` (`SimpleReactor`) with electron/`Te` machinery grafted in. Method roster below. |
| `test/rmgpy/solver/plasmaTest.py` (896 lines) | new | Imports `from rmgpy.solver.plasma import PlasmaReactor`. Classes: `PlasmaReactorTest` (`test_solve`, `test_collider_model`, `test_specific_collider_model`), `TestPlasmaReactorTe` (`test_Te_opt_in_uses_Te_for_kinetics`, `test_default_uses_T_when_not_opted_in`, `test_Te_override_via_conditions`, `test_volume_uses_T_not_Te`). Helper `TeProbeKinetics(Arrhenius)`. This is the natural port Verifier. |

`plasma.pyx` method roster (all present in the archived file):
`__init__`, `__reduce__`, `convert_initial_keys_to_species_objects`, `get_const_spc_indices`,
`get_electron_species_index` (new vs simple), `initialize_model`, `_temperature_for_kinetics`
(new), `calculate_effective_pressure`, `generate_rate_coefficients`,
`get_threshold_rate_constants`, `set_colliders`, `set_initial_conditions`, `residual`, `jacobian`.

**No `plasma.pxd` exists on 99 — and none is needed.** RMG solver reactors (`simple.pyx` has no
`simple.pxd`) declare their `cdef public` attributes *inline* in the class body. `plasma.pyx`
follows that pattern (attributes declared inline at lines 55–118). The generic CLAUDE.md ".pxd
per cdef class" rule does not apply to this module; do **not** author a `plasma.pxd`.

---

## Category 2 — Edits to existing infrastructure required by the port

These are the load-bearing ones. Split into **carry** (genuinely plasma-required) and **DO NOT
CARRY** (mainline drift `99` predates — carrying reverts a later bugfix).

### 2a. Carry — required, small, clean

| File | Change | Old (HEAD) → New (needed) |
|---|---|---|
| `rmgpy/solver/__init__.py` | register the class | add `from rmgpy.solver.plasma import PlasmaReactor` after the `SimpleReactor` import (line 32). One line. |
| `setup.py` `ext_modules` | build wiring | HEAD's solver list (lines 131–136) has **no** `plasma.pyx`; 99 lists it at line 134. Add `'rmgpy/solver/plasma.pyx'`. **Without this the `.so` is silently never built** (CLAUDE.md gotcha) and `import` fails at runtime. |
| `rmgpy/solver/base.pyx` | electron-temperature condition hook | in the `set_conditions` block (~line 208 on HEAD): add `if 'Te' in keys and hasattr(self, 'Te'): self.Te = Quantity(conditions['Te'])`. `hasattr`-guarded, so it is inert for non-plasma reactors. `self.Te` is **not** declared in `base.pxd`; `PlasmaReactor` declares `cdef public ScalarQuantity Te` itself. Clean. |

### 2b. DO NOT CARRY — `99` base edits that are stale drift, not plasma requirements

| File | Hunk on 99 | Why not to carry |
|---|---|---|
| `base.pxd` / `base.pyx` | replaces `if DASPK == 1: from pydas.daspk … else: from pydas.dassl …` with a **hardcoded** `from pydas.dassl cimport DASSL as DASx` | Removes the `settings.pxi`-driven compile-time DASPK/DASSL selection. HEAD's branch is correct; carrying this **breaks DASPK builds**. *(adv. — provenance:* the hardcode is side-branch commit `45c2e1c194 "Compile with the DASSL solver"`, unrelated to plasma logic.) Keep HEAD's version verbatim. |
| `base.pyx` | removes the `if not self.constant_volume:` guard around `RTP = R*T/P` in the sensitivity path (~line 713) | HEAD guards this; `99` computes `RTP` unconditionally. For a constant-volume reactor `P` can be unset → this is a **regression for Simple/Liquid/Surface constant-V sensitivity**. *(adv. — provenance:* HEAD's guard is the later explicit bugfix `6694323183 "Don't compute RTP if constant volume reactor"`. `99` predates it.) Keep HEAD's guard. |
| all touched solver files | copyright header `2002-2026` → `2002-2023` | Pure age artifact. Ignore. |

### 2c. Flagged for judgment — do not carry without justification

| File | Hunk on 99 | Concern |
|---|---|---|
| `base.pyx` | `from rmgpy.exceptions import NetworkError` + broadens the step-error catch from `except DASxError` to `except (DASxError, NetworkError)` (~line 732) | This **widens what the solver swallows** during `step()` — logs `logging.error` and continues instead of raising. That is exactly the error-suppression class the campaign exists to remove. It is not obviously plasma-required (plasma networks don't uniquely raise `NetworkError`). Treat as suppression-candidate: require an explicit justification and a live-path test, or drop it. |

---

## Category 3 — Single-temperature assumptions that conflict with the two-temperature design

The 2-T EOS itself is settled and correct in two of three places. The conflicts are where the
archived code still reasons at a **single** temperature and thereby contradicts its own `Te` path:

1. **`jacobian` uses the single-T EOS while `residual` uses the 2-T EOS — internal inconsistency
   (strongest item).** `jacobian` (lines 749, 751) computes `V = R*T*Σy/P` and
   `Ctot = P/(R*T)` with **no `Te`, no electron term**, whereas `residual` (lines 577–588) and
   `set_initial_conditions` (lines 477–484) use `V = (R/P)(N_heavy·T + N_e·Te)`. When `Te ≠ T`
   the analytic Jacobian no longer matches the residual's volume, degrading Newton convergence
   and corrupting sensitivity coefficients. The 2-T port is half-done here. Must be reconciled.

2. **Reverse rates / equilibrium constant are always evaluated at heavy-gas `T`.**
   `generate_rate_coefficients` (line 342): `self.Keq[j] = rxn.get_equilibrium_constant(self.T.value_si)`
   and `kb = kf/Keq`, even for reactions whose *forward* `k` was evaluated at `Te` (via
   `_temperature_for_kinetics`). Detailed balance in a 2-T plasma is subtle; pairing a `Te`
   forward rate with a `T` equilibrium constant is a modeling decision that needs an explicit
   ruling, not a silent default. The same `Keq(T)` reuse recurs in the `residual` pdep-collider
   recompute (lines 542, 559).

3. **Filter threshold uses heavy-gas `T`.** `get_threshold_rate_constants` (line 353):
   `unimolecular_threshold_rate_constant = 2.08366122e10 * self.T.value_si` (≈ kB·T/h). Electron-
   driven unimolecular processes would be thresholded at the wrong temperature. Minor, but a
   single-T assumption to record.

4. **Design intent locked by tests.** `TestPlasmaReactorTe.test_volume_uses_T_not_Te` asserts the
   EOS volume tracks `T` (the electron term `N_e·Te` is negligible at realistic ionization
   fractions), and `test_Te_opt_in_uses_Te_for_kinetics` asserts `Te` is used for kinetics only
   when `kinetics.uses_electron_temperature` is set. So the intended split is: **volume/EOS from
   `T`; rate constants optionally from `Te`.** Items 1–2 are where the archived code violates that
   split. The port must decide them, not inherit them.

5. **The whole `Te`-in-EOS path is largely DORMANT in `plasma.pyx` alone *(adv.)*.** The 2-T EOS
   only differs from single-T when `N_e = y0[electron_index] > 0`. But `plasma.pyx` **never
   converts `electron_density` into electron moles** — it stores the density, may append an
   electron species (with `y0[e] = 0` unless a mole fraction is separately supplied), and stops.
   The density→moles conversion lives in `rmgpy/rmg/main.py` (see Q3), which HEAD lacks. So in a
   `plasma.pyx`-only port, `N_e ≡ 0`, the EOS reduces to single-T everywhere, and item 1's
   inconsistency is masked — which is exactly why the archived `test_volume_uses_T_not_Te` passes.
   The inconsistency is real and must be fixed, but it is currently hidden by a deeper gap.

---

## Category 4 — Error-suppression and silent-corruption (do NOT carry blind)

**Exception-swallowing in `plasma.pyx`: none.** Rigorous check — full 1214-line read plus
`grep -nE "try:|except|pass|\.warning|\.error|getattr\(|= *0|continue"`:

- **Zero `try`/`except` blocks.** No exception swallowing anywhere in the reactor.
- Two **defensive `raise`s** (good, keep): `initialize_model` raises `ValueError` when
  `electron_density` is set but no electron species exists (lines 266–268);
  `set_initial_conditions` raises `ValueError` on non-positive mole-fraction **sum** (line 447).
- `getattr(kin, 'uses_electron_temperature', False)` (line 290) is **feature-detection for the
  opt-in `Te` path**, not error suppression.

**The only exception-swallowing in the archived diff is the `base.pyx` `NetworkError` catch (2c).**

But "no try/except" is not "no silent corruption." The reactor has several **silent behavior
changes / latent bugs** that must be fixed or consciously accepted during the port, not copied
blind *(adv.)*:

1. **Reverse rate for specific-collider pdep reactions is mis-scoped.** `residual` line 559,
   `kr[j] = kf[j] / equilibrium_constants[j]`, is **dedented out of the `for i` loop**, so it runs
   once after the loop against the *last* `j` only. With more than one specific-collider pdep
   reaction, every earlier reaction keeps a **stale `kb`**; if the last one took the
   `else: kf[j]=0` branch, `kr` is computed from a zeroed `kf`. Indentation bug — `simple.pyx` is
   the correct reference.
2. **Negative individual mole fractions pass.** The `X_sum > 0` check (line 447) guards only the
   *sum*; a component with a negative mole fraction paired with larger positives still yields a
   positive sum and produces a **negative species mole** in `y0`. `simple.pyx` has the same shape,
   but it is a silent-corruption risk to record.
3. **`electron_density` silently produces zero electron concentration** — see Category 3 item 5 /
   Q3. A nonzero density can leave `N_e = 0` with no error.
4. **`kf[j] = 0` for a missing specific collider** (line 558) is silent **model deletion** of that
   reaction, not merely harmless edge handling.

---

## Category 5 — Stale API calls (verified against HEAD)

Every non-trivial external call in `plasma.pyx` was checked against today's tree. **Result: the
API surface the reactor needs already exists on HEAD** — the charge-transfer Arrhenius tickets
already landed the plasma kinetics classes. Details:

| Call in plasma.pyx | HEAD status | Evidence |
|---|---|---|
| `spc.is_electron()` (line 229) | **present** | `rmgpy/species.py:558`, `rmgpy/molecule/molecule.py:363,1278` |
| `kin.uses_electron_temperature` (line 290) | **present** | `cdef public bint uses_electron_temperature` on `TwoTemperaturePlasma` (arrhenius.pyx:361), `ElectronCollisionPlasma` (:604), **and also `BadnellRRArrhenius` (:860) and `VoronovEIArrhenius` (:1451)** — the inventory's earlier "two classes" undercounted *(adv.)*. |
| electron species via `from_adjacency_list('1 e u1 p0 c-1')` (line 234) | **element present; parse to verify at build** | `e` is a real element (`element.py:178`, valences/valence_electrons include `'e'`). The exact `1 e u1 p0 c-1` adjlist should be round-tripped once at port-build time; nothing suggests it is stale, but it is the one item to smoke-test rather than assume. |
| `kin.get_effective_collider_efficiencies(core_species)` (line 383) | **present** | `cpdef` on `rmgpy/kinetics/model.pxd:76`; used identically by `simple.pyx:315`, `mbSampled.pyx:259` |
| base methods: `set_initial_conditions`, `set_initial_derivative`, `compute_network_variables`, `set_initial_reaction_thresholds`, `compute_rate_derivative`, `get_species_index`, `initialize_solver` | **all present** | grep of `rmgpy/solver/base.pyx` |

No *missing* calls found, but there is a **semantic API mismatch that a naive port would ship as a
scientific bug** *(adv., verified)*:

> **`_temperature_for_kinetics` routes `Te` through the wrong entry point.**
> `plasma.pyx` calls `rxn.get_rate_coefficient(Tkin, Peff)` with `Tkin = Te` for
> electron-temperature reactions (lines 338–339, 541, 556). But `Reaction.get_rate_coefficient(self,
> T, P=0, ...)` (reaction.py:912) has **no `Te` parameter**, and HEAD's
> `TwoTemperaturePlasma.get_rate_coefficient(T)` is defined as `return get_rate_coefficient_two_temp(T, T)`
> (arrhenius.pyx:451) — i.e. it means `k(T, Te=T)`. So passing `Tkin = Te` evaluates **`k(Te, Te)`**,
> collapsing the gas-temperature activation the 2-T rate law exists to keep separate. The real
> two-temperature entry is `get_rate_coefficient_two_temp(T, Te)` (arrhenius.pyx:428), which the
> `Reaction` wrapper does not currently expose. **The port must plumb `(T, Te)` through
> `Reaction`/kinetics to `get_rate_coefficient_two_temp`, not reuse the single-argument
> `_temperature_for_kinetics` shortcut.** This is the one place the ticket's "interfaces have since
> moved" narrative is genuinely true — HEAD grew a two-temperature kinetics API that the archived
> reactor predates and does not call.

Everything else the reactor calls (`is_electron`, `get_effective_collider_efficiencies`, electron
element `e`, all base methods) is present and unchanged. The genuine friction is: this kinetics
mismatch, the build-wiring (2a), the base-drift you must *avoid* carrying (2b/2c), and the
out-of-`plasma.pyx` electron machinery (Q3).

---

## The four design questions

### Q1 — What does the state vector hold? Is there an electron energy balance?
The DAE state vector `y` is **moles of core species only** (electron included as an ordinary
gas-phase species at `electron_index`). There is **no electron-energy equation and no gas-energy
equation** — the reactor is isothermal in *both* `T` and `Te`. `Te` is an externally imposed
constant (or a 2-element range across sims via `Te_range`), never integrated. Volume is not a
state variable either; it is recomputed algebraically from the EOS every `residual` call
(lines 577–588). So: species moles in, everything thermal held fixed.

### Q2 — Equation of state
Two-temperature ideal gas: `P·V = R·(N_heavy·T + N_e·Te)` (`set_initial_conditions` 477–484,
`residual` 577–588), degenerating to `P·V = R·T·N_total` when there is no electron species or no
`Te`. **Settled; not to be revisited** — flagged only that `jacobian` fails to use it (Category 3
item 1).

### Q3 — Zero-electron limit, and the electron-density dependency
The zero-electron *limit* is a clean graceful degradation: with no electron and no `Te`, both EOS
branches fall back to standard single-T ideal gas and the reactor behaves as a `SimpleReactor`.

But the **non-zero** electron path is only half-present in `plasma.pyx` *(adv., verified — this
corrects the first draft's "handled coherently")*: `get_electron_species_index` (224–238) locates
or synthesizes a 1-atom electron species via `from_adjacency_list('1 e u1 p0 c-1')`, but
**`plasma.pyx` never converts `electron_density`/`ne_range` into an electron mole fraction or into
`y0[electron_index]`.** That conversion — number-density → mole-fraction range, then injection —
lives in `rmgpy/rmg/main.py._initialize_electrons_from_density()` (`99:main.py:769`, called at
`:517`), and `Te_range` sampling lives in `RMG_Memory`/the driver (`99:main.py:892, 2358`). **HEAD
has none of this** (`grep electron_density rmgpy/rmg/main.py` → empty). Consequences for a
`plasma.pyx`-only port:

- `electron_density` becomes an **effectively dead constructor parameter**: it is stored, guards
  `get_electron_species_index`'s synthesis, and is checked by the `initialize_model` raise — but
  never reaches `N_e`, so the 2-T EOS silently degrades to single-T (this is what dormant-izes
  Category 3 item 1). The archived tests never assert electron amount/density, so they don't catch
  it.
- **`Te_range` and `ne_range` are dead inside `plasma.pyx`** — assigned in `__init__`, never
  sampled by any reactor method. A ranged plasma reactor initializes with `self.Te is None` unless
  an external `conditions` dict supplies a sampled `Te`.

So the port's true dependency set includes the `main.py` electron/`Te` driver logic, **not just**
the input DSL. That is a scoping correction, not a footnote.

### Q4 — Construction path
On `99` the path is `plasmaReactor {...}` input block → `plasma_reactor(...)` builder in
`rmgpy/rmg/input.py` (line 443, with `electronTemperature`, `electron_density`, `T`/`P`
scalars-or-ranges) → `PlasmaReactor.__init__`, **plus** the driver-side
`_initialize_electrons_from_density()` and `Te` sampling in `rmgpy/rmg/main.py` (Q3). HEAD has
neither the `input.py` builder nor the `main.py` machinery. The `input.py` builder is defensible as
a later input-DSL ticket, but the `main.py` electron/`Te` logic is **not** input sugar — without it
the reactor's `electron_density`/`Te`-range features do nothing (Q3). The reactor's own `__init__`
signature is the stable contract boundary:
`PlasmaReactor(T, P, initial_mole_fractions, electron_density, n_sims=1, termination=None,
sensitive_species=None, sensitivity_threshold=1e-3, sens_conditions=None, const_spc_names=None,
Te=None)`.

---

## Additional port landmines (from adversarial review, verified)

- **`__reduce__` drops the ranges.** `plasma.pyx` `__reduce__` (177–189) passes `self.T`, `self.P`,
  `self.Te` (each `None` when a range was given) and only falls back for the electron arg
  (`ne_arg`). So an unpickled ranged plasma reactor comes back as `T=None/P=None/Te=None` —
  `Trange`, `Prange`, `Te_range` are lost. RMG pickles reaction systems for multiprocessing; fix
  `__reduce__` to preserve all three ranges before relying on ranged plasma runs.
- **Synthetic electron vs object-keyed indexing.** `get_electron_species_index` **mutates the
  caller's `core_species` list** (appends an electron) *before* `ReactionSystem.initialize_model`.
  The base solver indexes species by **object identity** (`species_index[spec]`, base.pyx ~:425),
  not by isomorphism. Any electron-consuming reaction that references a *different* electron
  `Species` instance will not map to the synthesized one → `KeyError`/mis-index. Safe only if all
  electron reactions already share that exact object, or there are none.
- **`get_const_spc_indices` lost its duplicate guard.** HEAD's `SimpleReactor.get_const_spc_indices`
  has an `else: return` when `const_spc_indices` is already populated (simple.pyx:172–188); the
  archived `plasma.pyx` version (210–222) omits it, so a second `initialize_model` call **appends
  duplicate indices**. Restore the guard when porting.
- **Electron adjacency `1 e u1 p0 c-1` is UNPROVEN by execution.** Element `e` and `is_electron()`
  exist on HEAD (Category 5), but the exact adjlist string was not round-tripped; smoke-test it at
  port-build time rather than assume.

---

## What this probe could NOT reach
- **No build/compile was performed.** `plasma.pyx` was not cythonized against HEAD; the one
  parse-level unknown is the electron adjacency string (Category 5). Everything else is verified
  by source inspection, not execution.
- **Numerical correctness of the Jacobian/residual pair was not run** — the inconsistency in
  Category 3 item 1 is established by reading the two EOS expressions, not by a convergence test.
- **The detailed-balance question (Category 3 item 2) is a modeling decision, not a code fact** —
  flagged for a human/scientific ruling, not resolvable by inventory.
