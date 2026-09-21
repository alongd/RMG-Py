# I-246 — the charged-particle wall boundary operator

Branch `i246-ambipolar-wall-operator`, cut from `plasma` at `311818121`.
Envelope, observables and conditioning threshold: [`envelope.md`](envelope.md), committed at
`14042eebb` **before** any sweep was written or run — `git show --stat 14042eebb` contains that
file and nothing else.

Nothing here was pushed, merged, or opened as a PR. No file under `rmgpy/data/` and nothing in
`/home/alon/Code/RMG-database-plasma` was modified.

---

## The answer to the question this ticket exists to ask

> **Which quantities does the wall operator determine robustly, and which remain free or
> ill-conditioned until the discharge power balance is closed?**

**Robustly determined — by geometry and sourced transport data alone, with no reference to any
electron density:**

| quantity | evidence |
|---|---|
| `ν_wall = D_a/Λ²`, the common ambipolar loss frequency | matches an independent closed-form re-derivation to 0 relative error; 185.166 s⁻¹ at the nominal condition |
| its geometry scaling | halving `Λ²` doubles `ν_wall` to 15 digits; `ν_wall·Λ²` constant to 2.2e-16 across R = 5–200 mm |
| `Δν = ν_ion − ν_wall`, and the **sign** of it | classified at all 588 grid points, both sides sampled (347 / 241) |
| the **sustainment boundary**, `Te_thr(p·Λ)` | located by bisection at 42 (p, R) pairs; **collapses onto a similarity curve in p·Λ** |
| the **extinction decay rate** | equals `\|Δν\|` to 7.0e-6 relative, and is independent of the initial electron amount across 4 decades (spread 4.5e-5) |
| that volume recombination is irrelevant | `ν_RR/ν_wall ≤ 2.6e-12` and `ν_3b/ν_wall ≤ 2.4e-24` over every in-support point |

**Not determined — and the ticket's central finding:**

**The absolute electron density.** Not because the model cannot reach a realistic one — it can —
but because at a realistic one it is hopelessly conditioned. At 5 torr, R = 0.05 m:

```
n_e (m^-3)     Te needed (eV)   Te - Te_thr (eV)     C_wall = |dln n_e / dln nu_wall|
  1e16            0.852518          2.0e-07                    2.37e+05
  1e17            0.852520          2.0e-06                    2.37e+04
  1e18            0.852538          2.0e-05                    2.37e+03
```

A real argon glow discharge (n_e ~ 1e16–1e17 m⁻³) sits **2e-7 eV above the sustainment
threshold**, at `C_wall ≈ 2.4e5`. The ion mobility this rests on is known to ±3 % — the accuracy
of the Ellis/McDaniel/Albritton compilation itself — and 3 % on `ν_wall` there moves `n_e` by
`exp(0.03 × 2.4e5)`. That is not an error bar; it is the absence of a prediction.

The operationally useful form of the same statement: **`Te` would have to be known to a few parts
in 10⁷ of an eV for the electron density to be determined, and `Te` is prescribed in this model,
not solved for.** No wall operator can supply it. Closing the discharge power balance — the next
milestone, out of scope here — is what would, because it would determine `Te` instead of
prescribing it, and the distance to the boundary is the only thing the sensitivity depends on.

**So the gap does not close by fitting, and this is the result the brief asked for.** It was
reachable only because every wall parameter came from an independent source; a free wall
coefficient would have reproduced any of the rows in that table.

---

## 1. The opening probe, and where it contradicts the brief

`probe_opening.py`, every number taken through the compiled `PlasmaReactor.compute_volume`, never
a re-implementation. Log: `logs/probe-opening.stdout.log`.

| brief's claim | measured | verdict |
|---|---|---|
| electron carries 99.15 % of the two-temperature EOS sum at f = 1 | 99.1509 % | **confirmed** |
| recycling wall suppresses `n_e` response by ~1/117 | 1/117.8 | **confirmed** |
| pumping wall suppresses it by ~1/96,800 | **exactly 0** | **refuted** |
| the same event moves `n_e` and the neutral fraction orders apart | 4.34e-6 vs 5.11e+5 | **confirmed** |

With `f` the ionisation fraction (`N_h = N₀` is conserved because `Ar → Ar⁺` preserves heavy
count):

```
S_recycling(f) = Tgas / (Tgas + f*Te)
S_pumping(f)   = 1 - f*(Tgas + Te)/(Tgas + f*Te)
```

`S_pumping(1) = 0` identically and to all orders: a pumping wall removes heavy and electron moles
in exactly the ratio that scales the EOS sum and `N_e` together, so at complete ionisation it
cannot move `n_e` at all. 96,800 is a zero crossing, not a model property — `1/S_pumping` runs from
~119 at f→0 to ∞ at f→1 and passes 96,800 at f = 0.99878.

**Consequence, and it inverts the brief's (b):** the recycling-vs-pumping choice is not a bounded
~800× lever. It is a near-exact cancellation whose ratio is set by `f`, unbounded near the
operating point. That makes the no-fitting rule *more* load-bearing, not less — an unbounded knob
is worse than an 800× one.

### Reconciliation of the two parametrisations

A second derivation (using `f_φ = N_e·Te/S`, the electron's share of the EOS sum, rather than the
ionisation fraction) gives `1 − f_φ` and `1 − f_φ(Tg+Te)/Te`. These are **the same two functions**:
substituting `f_φ = f·Te/(Tg + f·Te)` reproduces the expressions above identically, and the zero
at `f_φ = Te/(Tg+Te) = 0.9915254` maps to `f = 1` exactly. The parametrisations differ only in
where the singular point sits — a domain edge in one, `f = 1` in the other — which is why
independent recomputations of the residue land on different figures (96,800, 39,000, 116,000).
**The residue is infinitely sensitive to the rounding of `f_φ` near its ceiling**, which is a
stronger reason to discard the figure than any of the three values.

### Two further probe findings

**The wall residual is `−ν·y_i` with `V` cancelling exactly.** `res = (dC/dt)·V` with `C = y/V`, so
a first-order loss contributes `−ν·C_i·V = −ν·y_i`. The brief's warning that the wall Jacobian
"also perturbs V" is right, but **not for the reason given**: the `C = y/V` conversion cancels
algebraically (and, written as `−ν·y_i`, bit-for-bit — written the other way it is off by an ulp,
measured). The `dV/dy` coupling is real only because the ion mobility carries `1/n_gas`. It enters
through the **mobility**, not the concentration conversion.

**The steady-state conditioning has a closed form.** `C_wall = (1−f)/(2f)`, and
`C_Te = C_wall·|1 − dln k_iz/dln Te|`. The two are proportional with a constant ratio, so **`Te`
and `ν_wall` are exactly degenerate directions** — any error in the wall coefficient is
indistinguishable from an error in `Te`. That is the mechanism by which a free wall coefficient
can reproduce any target density, stated as algebra rather than as a warning. The factor of 2 in
the denominator is the composition-dependent mobility partly self-correcting; a constant-mobility
model would give `(1−f)/f`.

---

## 2. What was built

`rmgpy/solver/plasma.pyx`, plus the `plasmaReactor(...)` keywords in `rmgpy/rmg/input.py` and
their documentation in `documentation/source/users/rmg/input.rst`.

```
mu_i    = ionReducedMobility * mobilityReferenceDensity / n_neutral
D_a     = mu_i * k_B*Te/e                      (Te >> Ti limit of D_i(1+Te/Ti))
nu_wall = D_a / Lambda**2
loss_i  = nu_wall * n_i                        for EVERY charged species
```

Against the six constraints:

1. **Reactor operator, never chemistry.** The wall terms are added to `res`, deliberately *not*
   folded into `core_species_rates` — that array is the gas-phase chemistry diagnostic that model
   enlargement and the rate-ratio termination compare edge fluxes against, and a transport term is
   not a reaction flux. Keeping them apart is the one place the separation can be enforced rather
   than asserted. No kinetics entry was added anywhere.
2. **Geometry is an input.** `chamberGeometry` takes a named shape and its dimensions (cylinder /
   sphere / slab, lowest diffusion eigenmode) or a stated `diffusionLength`. Giving both is
   refused — two sources of truth for one number.
3. **One common `ν_wall` for all charged species.** Not a lifetime each. Net charge then decays at
   a rate proportional to the net charge itself, which is zero in a quasineutral gas — so the
   zero-net-wall-current (floating wall) rule holds **by construction**, measured at exactly `0.0`.
4. **The electron can be made algebraic** (`quasineutralElectron`), carried on a DAE
   charge-conservation row with no `dydt` term. Independently switchable from the wall on purpose,
   so the conditioning claim is a measurement against a control rather than an assertion.
5. **The seed is a declared mechanism.** `ionisationSource`, a volumetric external pair source,
   defaulted to nothing and documented against the cosmic-ray background (2–10 pairs cm⁻³ s⁻¹ at
   1 atm, density-scaled). It is what creates the sub-threshold branch.
6. Retained scope items are covered in §5's verifier table.

**Recycling.** `wallRecycling` (γ ∈ [0,1], default 1.0). Each ion's neutral counterpart is resolved
once at initialization by heavy composition; an ion with no counterpart, or an **ambiguous** one
(two neutral argon states in the core, say ground and metastable), is refused rather than guessed
at — recycling into a coin-flip would silently move population between species with different
chemistry.

---

## 3. Two defects found by measurement, both fixed

### 3.1 The domain check was in the wrong place, and did not work

The first implementation raised `PlasmaStateError` from inside `compute_nu_wall` when the
ionisation degree exceeded the validity ceiling. Measured behaviour:

```
 All zero matrix    (x10)
 DASPK--  ITERATION MATRIX IS SINGULAR.
SystemError: <built-in method format of str object> returned a result with an error set
```

**A Newton trial state is not a physical state.** The solver evaluates the residual far outside the
domain on its way to a converged step; a Python exception raised inside the Fortran callback is not
propagated, DASPK continues with a corrupted iteration matrix, and the run dies with something
unreadable — the exact opposite of the loud failure the check exists to give.

Fixed by moving the check to the **accepted-step boundary**: `compute_nu_wall` is now total (a
numerical floor on the neutral population keeps a wild trial state finite so the solver can reject
it on its own error test), and `check_wall_support()` runs on the initial composition and after
every accepted step, via `cpdef` overrides of both `advance` and `step`. `cpdef`, matching the base
class, because a `def` override of a `cpdef` method is bypassed by C-level callers — and `step` is
the one `ReactionSystem.simulate` actually calls.

Measured after the fix (`logs/confirm-branches.stdout.log`), at Te = 1.5 eV, above threshold:

| configuration | how it ended |
|---|---|
| wall, integrated electron | `PlasmaStateError: ionisation degree n_e/n_neutral = 0.003428 exceeds the declared validity ceiling 0.001` |
| wall, quasineutral electron | same, α agreeing to **13 significant figures** |
| wall + cosmic-ray source | same |
| **control: no wall at all** | ran to t = 1 s, α = 2.18e34 — nothing stops the growth, as it must |

### 3.2 The algebraic electron was a 175× conditioning trap

Constraint 4's stated purpose is to remove a stiff direction from the Jacobian. Measured, it did
the opposite: `cond₂` of the iteration matrix was **5.39e5 with the electron integrated and 9.41e7
with the algebraic row** — a 175× penalty from the substitution meant to help.

The cause is scaling, not correctness. The charge row's entries are the charges themselves, exactly
±1, sitting among differential rows of order 1e2–1e5; DASPK factors with partial pivoting and **no
equilibration**. Scaling one row of a DAE residual by a positive constant is exact — it changes the
factored matrix and nothing about the solution — and it removes all of it:

| | raw κ | error-weight scaled κ | equilibrated κ |
|---|---|---|---|
| integrated electron | 5.39e5 | 3.59e5 | 2.54e5 |
| algebraic, unscaled row | 9.41e7 | 4.70e7 | 2.62 |
| **algebraic, row scaled** | **6.64e5** | **4.11e5** | **2.62** |

`PlasmaReactor` now sets `charge_row_scale` once at initialization, from the largest magnitude
among the differential rows — what an equilibration step would use, and no chosen constant. The
dependence is flat over three decades of the multiplier.

**Honest reading:** with its row scaled, the algebraic electron is **conditioning-neutral** in the
matrix the solver actually factors, and ~1e5 better in the matrix an equilibrating solver would.
It is also independently correct (the 13-figure agreement above). It does **nothing** for the 784×,
which is not a linear-algebra quantity at all — see §4.

---

## 4. The 784×, relocated

The governing ruling measured that a 0.01 % perturbation in `Te` moved `n_e` by ~784×. Taken by
bracketing **along** the branch, never across the boundary, as `envelope.md` predeclared:

```
Te_thr - Te (eV)     factor on n_e per 0.01% in Te      C_Te
    1.0e-01                 1.000108                  1.08e+00
    1.0e-02                 1.007619                  7.59e+01
    1.0e-03                 1.092106                  8.81e+02
    1.0e-04                 6.769855                  1.91e+04
```

**The 784× is not a property of the model, the solver, or the Jacobian.** It is a reading taken at
a particular distance from the sustainment boundary, and it can be made any number at all by moving
along the branch: ~1 far from the boundary, divergent at it. The transferable statement is a
**distance**, not a factor.

This is why §3.2's result is not a disappointment. No repacking of the unknowns can make a quantity
determined that the equations do not determine.

---

## 5. The sweep

588 grid points: Te ∈ [0.5, 4.0] eV × p ∈ [0.5, 50] torr × R ∈ [5, 200] mm. `sweep.csv`,
`logs/sweep.stdout.log`. Classified by the rules frozen in `envelope.md`:

| class | count | share |
|---|---|---|
| volume-growth (steady state out of support) | 344 | 58.5 % |
| extinction | 239 | 40.6 % |
| near-threshold / ill-conditioned | 5 | 0.9 % |

**Both sides sampled: 347 points with Δν > 0, 241 with Δν < 0.**

### The similarity law — the strongest result in the sweep

`ν_ion ∝ p` and `ν_wall ∝ 1/(p·Λ²)`, so the boundary should depend on `p` and `Λ` only through
their product. It does, and the collapse was not fitted:

| p·Λ (torr·cm) | Te_thr (eV) | p (torr) | R (m) |
|---|---|---|---|
| 0.2077 | 1.52390 | 0.50 | 0.010 |
| 0.2079 | 1.52368 | 1.00 | 0.005 |
| 0.4142 | 1.33495 | 0.50 | 0.020 |
| 0.4154 | 1.33427 | 1.00 | 0.010 |
| 0.4157 | 1.33410 | 2.00 | 0.005 |
| 0.8285 | 1.18806 | 1.00 | 0.020 |
| 0.8308 | 1.18753 | 2.00 | 0.010 |

Two conditions differing by a factor of two in pressure *and* radius agree on the threshold
temperature to **1.4e-4 relative**. The residual scatter is explained: `Λ` also carries the axial
term `(π/L)²`, which does not scale with R, so `p·Λ` is not an exact invariant for a finite
cylinder.

### The in-support window

The constant-pressure stationary state is `n_Ar* = sqrt(K/k_iz)` — a **square root**, not the
linear `ν_wall/k_iz` a constant-mobility model gives. Above threshold it runs quickly out of the
ion-neutral transport model's domain. The band in `Te` between the boundary and the ceiling is
narrow: **relative width 2.0e-3 to 1.7e-2** across the envelope, with `C_wall` = 7–20 at the top of
it — inside the predeclared limit of 1e2.

> **A correction to my own arithmetic, made by the measurement.** I had argued the in-support and
> well-conditioned regions were provably incompatible. They are not: the ionisation *fraction* `f`
> and the ionisation *degree* `α` differ by the temperature ratio,
> `α = f/[(1 + Te/Tgas)(1−f)]`, so at `Te/Tgas ≈ 33` an in-support `α = 1e-3` means `f ≈ 0.033`,
> not `f ≈ 1e-3`. Reading one as the other overstates `C_wall` by ~34× and made a usable window
> look unusable. Verified against the sweep's own `α` to 0 relative error.

### The sub-threshold branch

`n_e = S_ext/(ν_wall − ν_ion)`, which exists only where the wall wins. At 5 torr, R = 0.05 m it runs
from 1.4e3 m⁻³ far below threshold to 1.1e5 m⁻³ near it. **Its width is exactly the cosmic-ray
source interval, a factor of 5.0 at every point, and that is not narrowed anywhere** — the absolute
sub-threshold density is known no better than the background ionisation rate is.

---

## 6. Verifier

| # | requirement | evidence |
|---|---|---|
| 1 | correct dimensions **and population normalisation** | `ν_wall` matches an independent closed form to 0 rel; `−(wall flux)/N_e = ν_wall` exactly; the per-*neutral* normalisation is a factor 1e-6 away, asserted so the trap stays live |
| 2 | doubling `A_eff/V` doubles the loss frequency | halving `Λ²` gives ratio **2.000000000000000** |
| 3 | geometry affects flux per the transport model | `ν_wall·Λ²` constant to 2.2e-16 over R = 5–200 mm |
| 4 | current-balance rule obeyed | net wall charge flux **exactly 0.0** on a neutral state; negative control: non-neutral gives non-zero, and the ratio is `−ν_wall` |
| 5 | no spurious charge drift | `d(net charge)/dt` from the full residual = **0.0** against term magnitudes of 94 mol/s |
| 6 | gas+wall argon conserved under the declared model | heavy loss = `(1−γ)·ν·N_ion` exactly, at γ = 1.0, 0.5, 0.0 |
| 7 | **zero wall area recovers the pre-change equations bit-for-bit** | a second extension built from `git show 311818121:rmgpy/solver/plasma.pyx`; residual, Jacobian, `jacobian_matrix`, `y0`, `dydt0` and the EOS **identical** across 7 ionisation degrees × 3 `cj` values with non-zero `dydt`; negative control confirms a walled reactor does differ |
| 8 | `ν_wall > ν_ion` decays to extinction | `N_e(1 s)/N_e(0) = 1.1e-9`; decay rate = `\|Δν\|` to 7.0e-6 |
| 9 | no negative amounts near extinction | min `N_e` = 1.1e-18 mol, min `N_Ar+` = 1.1e-18 mol, both positive |
| 10 | residual and Jacobian identical; analytic vs FD away from the boundary | **step scan**, not one step: best agreement 9.9e-7 (integrated), 2.0e-7 (quasineutral), 9.9e-7 (+source), 1.0e-6 (γ=0.5); **no-wall control arm gives the same 9.9e-7**; negative control (γ 0.5→0.6 between analytic and FD) fails at 1.6e-1 |
| 11 | sweeps reproduce the crossover, both sides | 347 / 241; similarity collapse to 1.4e-4 |
| 12 | **no target electron density in the implementation or criteria** | the wall code introduces exactly two numeric constants: the Loschmidt density (a unit convention — the density the tabulated mobility is normalised to) and a dimensionless ionisation-degree ceiling. Asserted in the harness. Every transport parameter's source is in `envelope.md` §1 |
| 13 | ordinary reactors unchanged | see §7 |
| 14 | every new test has a **demonstrated red state** | **16 of 16**, in `logs/red-states.log`: test name, the exact diff hunk of the breakage, the captured failing output, and the revert confirmed by checksum against `bb4ce4dfb37a74bba012795cff779180` |
| 15 | `git status` clean; both streams captured | all runs logged to `logs/*.stdout.log` + `logs/*.stderr.log` |

Verifier 10 deserves a note. At `h/|y| = 1e-7` **every** configuration disagrees at ~1e-3 —
including the pre-existing no-wall path — because `res[0]` is O(1) while the perturbation changes
it by O(1e-14). A single step size proves nothing about a Jacobian; the control arm is what
attributes the number.

Two notes on verifier 14, both of which the red-state work turned up:

- **`test_diffusion_length_and_mobility_must_be_declared_together`** fails when its guard is
  removed, but with an `AttributeError` rather than the `PlasmaStateError` it expects. The test
  correctly rejects that — a `pytest.raises(PlasmaStateError)` is not satisfied by an
  `AttributeError` — so it is a valid demonstration, and it incidentally shows what the guard is
  worth: without it the code crashes anyway, just unreadably. Turning a crash into a message is
  exactly the guard's job.
- The cycles are only meaningful because each one rebuilds. `plasma.pyx` is cythonized; a
  breakage without a rebuild runs the old `.so`, the test passes, and a red state gets logged that
  never happened. Every entry records the rebuild in both directions.

---

## 7. Suites

```
$ pytest test/rmgpy/solver/ test/rmgpy/rmg/ test/rmgpy/kinetics/ -p no:cacheprovider --no-cov -q
678 passed, 2 skipped, 11 warnings in 157.99s
```

678 = the 658 that pass on the branch base, plus the 20 new wall test IDs. Zero failures, zero
errors.

An earlier run of the same three suites reported **one additional error**, and it is worth
recording because it is environmental rather than a flake: `mainTest.TestMain`'s *teardown* tries
to `shutil.rmtree` `RMG-database-plasma/input/kinetics/libraries/testSeed`, i.e. the test harness
writes into the database directory, which is read-only for this ticket. It did not recur on the
final run, it is unrelated to this change, and the database was not written to either way.

Note for anyone reproducing: a fresh worktree has **no `rmgrc`**, and without it these three suites
report 26 failures and 21 errors that are entirely database-resolution, nothing to do with the
code. `cp rmgrc.template rmgrc` first.

---

## 8. What is NOT claimed

Self-consistent `Te`; absorbed power; discharge current; plasma potential; sheath structure;
absolute steady-state `n_e`; agreement with any glow-discharge measurement. No electron-energy or
power balance was implemented — that is the next milestone and was out of scope by instruction.

Known limitations, stated rather than discovered later:

- **Single-ion ambipolar.** One common `ν_wall` is what constraint 3 asks for and it conserves
  charge exactly for any number of ion species, but a genuinely multi-ion ambipolar model needs
  per-ion mobilities under a common field. Not implemented.
- **No sheath-area or electrode-area factor**, and no `h`-factor. The model is diffusion-limited
  loss to a single enclosing surface. Adding either would be adding a free parameter, which is
  what the ticket prohibits.
- **`D_a = μ_i·kTe/e`** drops `(1 + Ti/Te)` — a 0.85 % underestimate at these ratios, carried
  inside the stated interval rather than corrected away.
- **Unmagnetised, no flow, closed batch**, infinite residence time.
- **`ν_RR` and `ν_3b` use order-of-magnitude recombination coefficients.** Their only role is to
  establish that volume recombination cannot compete with the wall in the supported regime
  (≤ 2.6e-12), and no conclusion rests on their precision.
- **Pre-existing, unrelated, not fixed here:** `ReactionSystem.set_initial_derivative` computes
  `dydt0 = -residual(t0, y0, 0)`, and every RMG reactor writes `delta = res - dydt`, so `dydt0`
  comes out as **minus** the true initial derivative. Project-wide, in `base.pyx`, and outside this
  ticket. Probed explicitly (`logs/probe-dae-row.stdout.log`, D4): DASPK is insensitive to it, here
  and for the algebraic row.

---

## 9. Reproducing

```bash
export PATH=/home/alon/anaconda3/envs/rmg_env/bin:$PATH
cp rmgrc.template rmgrc                      # a fresh worktree has none
python setup.py build_ext --inplace -j 8     # plasma.pyx is cythonized

cd docs/i246-ambipolar-wall-operator
export PYTHONPATH=$(git -C .. rev-parse --show-toplevel):.
python probe_opening.py          # the opening probe, against the compiled EOS
python probe_dae_row.py          # can the solver integrate the algebraic row at all
python verify_operator.py        # verifier checks 1-6, 8-12, exit 0 only if all pass
bash build_and_verify_zero_wall.sh   # check 7, bit-for-bit against a pre-change build
python sweep.py                  # the envelope sweep -> sweep.csv
python confirm_branches.py       # branches through the real integrator
python conditioning_ab.py        # integrated vs algebraic electron, and the 784x
python glow_discharge_gap.py     # where a real discharge sits, and its conditioning
```
