# I-199 — Justifying the stopping time of the 5 torr argon deck

**Deck:** `docs/i194-ar5torr-plasma-lineage/input.py`
**Question:** the deck stops at `terminationTime = (1e-3, 's')`, a placeholder its own
header admits is "a substitution" for the specified *"integrated to quasi-steady state"*.
Given a model that can only ionise, what stopping time is defensible, and what does the
number at that time mean?

**Answer, in one line:** there is no quasi-steady state to integrate to — the model has no
stationary state — and no finite exhaustion time either, since ionisation completes only
asymptotically; so the defensible stop is a **tolerance cutoff near complete ionisation,
1×10⁻⁴ s** (residual neutral fraction `1−f < 10⁻⁸`, and every core quantity stationary to ≥4
significant figures), and the run there represents **an all-but-fully-ionised snapshot, not a
steady state.** The deck is changed from 1×10⁻³ s to 1×10⁻⁴ s with that justification.

All numbers below are computed with RMG's own compiled machinery on this worktree
(RMG-Py `bb0a5d8e7`, RMG-database-plasma `3c1c57be2`); the trajectory is the actual solver
output of the deck. Load average at run time was 4.49 (`/proc/loadavg`).

---

## Q1 — The ionisation frequency, computed

`ElectronCollisionPlasma.integrate_rate_coefficient(Te)`
(`rmgpy/kinetics/arrhenius.pyx:570`) is the Maxwellian average the run actually uses:

```
k(Te) = √(8/(π mₑ)) · (k_B Te)^(-3/2) · Na · ∫ σ(E) · E · exp(-E/k_B Te) dE
```

integrating the `PlasmaArgon` 51-point Golyatina2021 cross-section table against a
Maxwellian EEDF at the deck's `Te = 34813.5 K` (= 3 eV). Evaluated directly:

| quantity | value |
|---|---|
| `k_iz(3 eV)` | **8.299×10⁷ m³ mol⁻¹ s⁻¹** = 1.378×10⁻¹⁰ cm³ s⁻¹ |
| `n_Ar` at 5 torr (666.612 Pa) / 298.15 K = P/(k_B T) | **1.619×10²³ m⁻³** |
| `nu_iz = k_iz · n_Ar / Na` (electron multiplication frequency) | **2.232×10⁷ s⁻¹** |
| e-folding time `τ = 1/nu_iz` | **4.481×10⁻⁸ s (44.8 ns)** |

**Units.** `integrate_rate_coefficient` multiplies by `Na` before returning
(`arrhenius.pyx:693`), so `k_iz` is **molar** (m³ mol⁻¹ s⁻¹). The multiplication frequency is
therefore `nu_iz = k_iz · n_Ar / Na = k_iz · P/(R·T_gas)`, not `k_iz · n_Ar` (which is
dimensionally wrong by a factor `Na`); numerically 8.299×10⁷ · 1.619×10²³ / 6.022×10²³ =
2.232×10⁷ s⁻¹.

(Sensitivity: `k_iz` climbs steeply with Te — 2.08×10⁻¹⁵ cm³/s at 1 eV, 7.89×10⁻¹² at 2 eV,
1.38×10⁻¹⁰ at 3 eV, 6.01×10⁻¹⁰ at 4 eV — so τ ranges from ~3 ms at 1 eV to ~10 ns at 4 eV.
The 3 eV working point is on the steep part of the curve.)

`nu_iz` is only the *initial* small-signal growth rate (argon undepleted, no volume change);
it is the slope of the full nonlinear trajectory at `f→0`. As ionisation proceeds the growth
slows sharply — the reactor is constant-pressure and each ionisation adds a hot electron, so
the volume expands and dilutes the concentrations — see the governing equation in Q3. The
final approach to complete ionisation runs on a much slower clock,
`kP/(R(T_gas+Te)) = 1.90×10⁵ s⁻¹` (τ_late = 5.28 µs).

## Q2 — How many e-foldings 1×10⁻³ s is

**1×10⁻³ s / 44.8 ns = 2.23×10⁴ e-folding times** (measured against the *initial* clock).
That framing is a lower bound on the real timescale, because the growth slows as ionisation
proceeds (Q3): the run reaches near-complete ionisation not "in well under a microsecond" but
at ~1×10⁻⁴ s. Either way the 1×10⁻³ s stopping time is **an order of magnitude past the point
at which the model has effectively nothing left to ionise** (residual neutral fraction < 10⁻⁸;
Q3 table). This is a finding about the current deck, not a detail: the placeholder sits inside
a static, essentially fully-ionised plateau for ~90% of the run.

## Q3 — What the electron density actually does (the run's own trajectory)

Trajectory saved to `solver/simulation_1_3.csv` (890 points; `saveSimulationProfiles=True`).
Key points (mole numbers normalised to 1 mol argon at t=0; `f` = ionised fraction
`Ar⁺/(Ar+Ar⁺)`; `n_e` = electron number density):

| time (s) | f (ionised) | n_e (m⁻³) | note |
|---|---|---|---|
| 0 | 6.18×10⁻⁸ | 1.00×10¹⁶ | seed |
| 1×10⁻⁶ | 0.069 | 1.22×10²¹ | avalanche |
| 4.45×10⁻⁶ | 0.50 | ~7×10²⁰ | half-ionised |
| 8.79×10⁻⁶ | ~0.80 | ~1.1×10²¹ | (rate-ratio would fire here — Q4) |
| 2.51×10⁻⁵ | 0.99 | 1.36×10²¹ | |
| 1×10⁻⁴ | 1−f = 6.8×10⁻⁹ | 1.375×10²¹ | near-complete ionisation |
| 1×10⁻³ | 1−f → 0 (mass-action asymptote) | 1.375×10²¹ | static (identical to 1×10⁻⁴) |

**n_e grows by a nonlinear avalanche and then saturates as the neutral argon is
asymptotically consumed — it does not grow without bound, and there is no dynamic balance.**

*Governing equation.* The seeding is charge-neutral (`Ne ≈ N_Ar⁺`), and the heavy count is
conserved (`N_Ar + N_Ar⁺ = N_h0`), so a single variable — the ionised fraction
`f = N_Ar⁺/N_h0` — closes the system. But this reactor is **not** fixed-volume: it is
constant-pressure and two-temperature, `V = R(N_heavy·T_gas + Ne·Te)/P`
(`plasma.pyx:831`), with concentrations `C = N/V` (`plasma.pyx:931`). Substituting
`N_Ar = N_h0(1−f)`, `Ne = N_h0 f`, and `V = (R N_h0/P)(T_gas + Te·f)` into the extensive rate
`dN_Ar⁺/dt = k·C_Ar·C_e·V` (with `k` molar) gives, exactly,

```
df/dt = (kP/R) · f(1−f) / (T_gas + Te·f)          [Riccati form, NOT a logistic]
```

Only in the (false) fixed-volume limit `Te = T_gas = const, V = const` would this collapse to
a logistic `d[e⁻]/dt = k[e⁻](n_Ar−[e⁻])`. Here the extra `1/(T_gas + Te·f)` factor — the
volume dilution, with `Te/T_gas ≈ 117` — dominates the dynamics as soon as `f` is even a
percent. Its consequences, all matched by the run:

- **Initial slope** (`f→0`): `df/dt → (kP/(R·T_gas))·f = nu_iz·f`, so the early e-folding time
  is the `τ = 44.8 ns` of Q1.
- **Final approach** (`f→1`): `df/dt → (kP/(R(T_gas+Te)))·(1−f)`, so `1−f` decays with the
  ~117× slower constant `τ_late = R(T_gas+Te)/(kP) = 5.28 µs`. **f approaches 1
  asymptotically; there is no finite exhaustion time.**
- **Exact trajectory:** `T_gas·ln f − (T_gas+Te)·ln(1−f) = (kP/R)·t + const`. This puts
  `f = 0.5` at 4.37 µs (run: 4.45 µs) — not the ~0.5 µs a naive `ln(n_Ar/n_e0)/nu_iz` estimate
  gives, precisely because of the dilution term.

*Plateau value.* The saturated electron density is not `n_Ar`. It is **1.375×10²¹ m⁻³**, set
by the same constant-pressure two-temperature constraint: at full ionisation the gas is Ar⁺
(298.15 K) plus an equal number of electrons (34813.5 K), and holding total pressure fixed
gives

```
n_e,sat = n_Ar · T_gas/(T_gas + Te) = 1.619×10²³ · 298.15/35111.65 = 1.375×10²¹ m⁻³
```

which matches the run to three figures. The hot electrons carry ~99% of the pressure at
saturation, so complete ionisation expands the gas ×117.8. **After ~1×10⁻⁴ s every core
quantity is static to ≥4 significant figures for the remaining ~90% of the 1×10⁻³ s run.**

## Q4 — What `terminationRateRatio` does here (and a prior-art record it refutes)

`terminationRateRatio` is RMG's nearest thing to a steady-state test. Mechanism
(`rmgpy/solver/base.pyx`):

- line **873**: `char_rate = sqrt(Σ core_species_rates²)` — the L2 norm of the *net*
  core-species production/consumption rates (mol m⁻³ s⁻¹).
- lines **875–876**: `max_char_rate` = running maximum of `char_rate` over the integration.
- lines **1426–1428**: fires when `max_char_rate != 0.0 and char_rate/max_char_rate < ratio`.

For this deck the single reaction `Ar + e⁻ → Ar⁺ + 2e⁻` makes the three net rates
`(−R, +R, +R)` with `R = k·C_Ar·C_e` (concentration-space), so `char_rate = √3·R`.
`char_rate` is thus in **concentration** units and inherits the volume dilution of Q3.
Writing `R` in the ionised fraction `f` (using `C_Ar, C_e ∝ (1−f), f` over `V ∝ (T_gas+Te·f)`):

```
char_rate ∝ f(1−f) / (T_gas + Te·f)²
```

The `(T_gas+Te·f)²` denominator (volume-expansion, *squared* because both concentrations
dilute) is what governs the shape. `char_rate` therefore **peaks very early**, at
`f* = T_gas/(2T_gas + Te) = 0.0084` (≈ 0.84 % ionised), and then falls monotonically for the
rest of the avalanche — **not** near half-ionisation. `max_char_rate` is set at that ~0.84 %
point; from there `char_rate/max_char_rate` falls through `ratio` while dilution steepens.

**Measured:** with `terminationRateRatio=0.01` the run fires at
**t = 8.79×10⁻⁶ s, char_rate/max_char_rate = 0.00938** (`stdout.log`:
*"reached target termination RateRatio: 0.009377…"*), where `f ≈ 0.80`. At that point the
argon is only ~80 % consumed, so **the trigger is dilution-dominated, not reactant depletion:**
between `f*` and `f = 0.80` the numerator `f(1−f)` actually *rises* ~19×, but the squared
dilution denominator rises ~2270×, netting the ~100× drop that trips the ratio.

> **Prior-art correction.** An earlier ticket on this campaign recorded that for this deck
> the rate-ratio criterion "correctly never fires, because the model has no sink." **That
> record is wrong.** It fires, at 8.79 µs. But the reason it fires is not the reason one might
> first reach for either: `char_rate` declines because the constant-pressure volume expansion
> dilutes both reactant concentrations (∝ `1/(T_gas+Te·f)²`), long before the argon is
> depleted — reactant depletion would only supply the `(1−f)` factor, still 0.2 at the
> trigger. The absence of a *product* sink (Ar⁺/e⁻ recombination) is irrelevant to whether
> `char_rate` declines; the two-temperature dilution does the job. Verified against the solver
> source (`char_rate` at base.pyx:873; `V`, `C=y/V` at plasma.pyx:831/931) and by the run.

But note **what** it detects: it fires mid-avalanche, at ~80% ionisation, on the *rate*
collapse — the composition is still moving toward its terminal value at that instant. So
`terminationRateRatio` is not a reliable "no longer changing" test for a fuel-limited system
like this one: it triggers early, on the rate downslope, and its `"reached … RateRatio"` log
line reads misleadingly like convergence. That is why the deck uses a derived
`terminationTime` (below) rather than this criterion.

---

## Deliverable — the stopping time and its justification

**There is no defensible steady-state stopping time, because the model has no stationary
state — and there is no finite exhaustion time either, because mass action approaches
complete ionisation only asymptotically** (`1−f ∝ e^{−t/τ_late}`, Q3). So the stopping time
cannot be "the instant the argon runs out"; it can only be a **cutoff at a stated tolerance**,
justified as what it is. The deck keeps 1×10⁻⁴ s, defined by two criteria the run satisfies:

```
terminationTime = (1e-4, 's')
```

- **Composition tolerance:** the residual neutral fraction `1−f` falls below **10⁻⁸** at
  **t = 9.80×10⁻⁵ s** (run: `1−f = 6.8×10⁻⁹` at 1×10⁻⁴ s). Choosing ε = 10⁻⁸ is a stated
  convention, not a derivation; at 1×10⁻⁴ s the model is >99.999999 % ionised.
- **Stationarity tolerance:** every core quantity is constant to **≥4 significant figures**
  well before the cutoff — `n_e` reaches 4-significant-figure stability (relative change to its
  asymptote < 10⁻⁴) at **t = 2.45×10⁻⁵ s**, and at 1×10⁻⁴ s it differs from its asymptote by
  1.3×10⁻¹⁰. The independent diagnostic run to 1×10⁻³ s (a decade beyond) confirms nothing
  changes thereafter.

So 1×10⁻⁴ s is **a numerical cutoff at which the model is ε-ionised (ε = 10⁻⁸) and every core
quantity is already stationary to 4+ significant figures — not an instant of physical
completion.** For scale it is ≈ 2.2×10³ initial e-folding times of `τ = 44.8 ns`, but the
justification is the two tolerances above, not that count. The previous 1×10⁻³ s let the run
finish and is not "wrong"; it is ~10× past this cutoff and was mislabelled a "QSS substitute",
which — there being no quasi-steady state — it is not.

## Limitations (for direct quotation into a report)

This argon model resolves ionisation and nothing else. Its single reaction, electron-impact
ionisation of neutral argon, is irreversible and has no counterpart: there is no
recombination and no transport channel that removes electrons or Ar⁺. Consequently the model
has **no stationary state** — the ionisation runs asymptotically to complete ionisation of the
argon, and any quantity the run reports is a snapshot of that transient or of the trivial
essentially-fully-ionised end state it approaches once the neutral argon is all but consumed
(here below one part in 10⁸ by ~1×10⁻⁴ s). No electron density, ionisation degree, or plasma
composition drawn from this
model should be read as a physical steady state of a 5 torr argon discharge; a real discharge
reaches a stationary ionisation degree of order 10⁻⁷–10⁻⁴, orders of magnitude below the
complete ionisation this sink-free model runs to.

**The missing electron loss is wall transport, not volume recombination.** At 5 torr in a
centimetre-scale tube the dominant electron sink is ambipolar diffusion to the wall, whose
frequency exceeds that of radiative recombination by roughly seven orders of magnitude.
Adding a volume-recombination reaction (`Ar⁺ + e⁻ → Ar`) would produce a stationary balance,
but at an ionisation degree wrong by orders of magnitude — a converged-looking number that
means nothing. The physically correct closure is a wall-loss operator, which is a later
milestone and is deliberately absent here; this is why the run is stopped at a stated,
computed transient time rather than being driven to a fictitious equilibrium.

Two further errors in the ionisation source term itself are known, and they act in
**opposite directions.** First, the mechanism carries only *direct* ground-state ionisation
and omits *stepwise* ionisation through the Ar(4s) metastable manifold, which at a few torr
and a few eV is a significant and often dominant channel; its absence **underpredicts** the
true ionisation rate. Second, the rate is computed by averaging the cross-section over a
**hardcoded Maxwellian** electron energy distribution, whereas a real low-pressure argon
discharge has a depleted high-energy tail (Druyvesteyn-like); because the ionisation
cross-section lives entirely above the 15.76 eV threshold, the Maxwellian tail
**overpredicts** the population of ionising electrons and hence the rate. **The model does
not resolve which of these opposing errors dominates**, so the computed `k_iz`, the avalanche
timescale, and the cutoff time should be treated as order-of-magnitude, not quantitative.

---

## Reproduction

```
conda activate rmg_env
cd docs/i194-ar5torr-plasma-lineage        # rmgrc pins RMG-database-plasma/input
python /path/to/rmg.py input.py > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)
```

`k_iz` and `nu_iz` are computed by instantiating `ElectronCollisionPlasma` with the
`PlasmaArgon` table and calling `integrate_rate_coefficient(34813.5)`; the trajectory is
`solver/simulation_1_3.csv` from a run with `saveSimulationProfiles=True`; the rate-ratio
firing time is from a run with `terminationRateRatio=0.01` added.
