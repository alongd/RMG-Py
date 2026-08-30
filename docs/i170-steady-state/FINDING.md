# I-170 — what a plasma trajectory actually does, and the criterion that follows from it

Measured 2026-08-30 on branch `i170-steady-state`, cut from `plasma` at `4f18dc389`, against
RMG-database-plasma at `053d4808`. Nobody on this campaign had plotted a plasma trajectory
before this; the numbers below are the first.

Artifacts in this directory:

| file | what it is |
|---|---|
| `trajectory-argon.csv` | the I-164 5 torr argon deck's final simulation, 51 solver steps |
| `trajectory-lithium.csv` | the lithium-seeded design case, 498 solver steps |
| `design-case-lithium.py` | the deck that produced it |
| `emulate.py` | replays a saved profile through candidate criteria |

---

## 1. The argon deck's trajectory is not "nearly trivial". It is exactly constant.

The contract warned the argon deck might be boring. It is worse than boring, and the exact
form matters for what follows.

`trajectory-argon.csv` holds 51 rows spanning `t = 0` to `1.0995e-3 s`. **Every mole fraction
in every row is bit-identical to its `t = 0` value**, as is the volume:

```
Time (s),Volume (m^3),He,Ne,N2,e-(1),Ar(2),Arp(3),[Li](4),[Lip](5)
0.0,                3.7187721422257316,0,0,0,6.175165367910961e-08,0.9999998764966928,6.175165350358922e-08,0,0
1e-15,              3.7187721422257316,0,0,0,6.175165367910961e-08,0.9999998764966928,6.175165350358922e-08,0,0
...
0.001099511627775,  3.7187721422257316,0,0,0,6.175165367910961e-08,0.9999998764966928,6.175165350358922e-08,0,0
```

The cause is I-164's B2: the plasma libraries carry lithium and only lithium, so an argon deck
generates no argon reaction at all. The final core has 8 species and 2 reactions, and both
reactions are lithium's, on species present at exactly zero (I-164 B3).

**Consequence for this ticket.** A tolerance cannot be argued from this trajectory: the residual
is identically zero at every step, so *every* tolerance terminates at the first step and none is
distinguishable from any other. Fitting a criterion to it would produce a number with no
information in it. Three further core species — He, Ne and N2, injected by RMG and never asked
for (I-164 B1) — sit at exactly `0.0`, which is also why the residual needs a floor: their
relative change is `0/0`.

## 2. A design case that does relax: lithium-seeded argon

Lithium is the one element for which the plasma libraries carry **both** directions of the
ionisation balance — `[Li] => [Lip]` (Voronov electron-impact ionisation) and `[Lip] => [Li]`
(Badnell radiative recombination). An alkali-seeded noble gas is a standard laboratory
configuration, so seeding 100 ppm of lithium into the I-164 argon deck, at identical conditions
(5 torr, 298.15 K gas, 3 eV electrons, `n_e = 1e16 m^-3`, charge-neutral by the same hand
arithmetic), gives a real case rather than a contrivance.

It relaxes over 498 solver steps:

| | t = 0 | t = 112 s |
|---|---|---|
| Li | 1.0000e-04 | 6.6692e-10 |
| Li+ | 6.1752e-08 | 1.0006e-04 |
| e- | 6.1752e-08 | 1.0006e-04 |
| V (m³) | 3.7187721422 | 3.7621938015 |

Lithium ionises essentially completely; the electron and cation densities rise by three orders
of magnitude; the volume moves 1.17% purely through the two-temperature equation of state, with
no chemistry in it — which is why the residual below is taken on **mole fractions** and not on
moles or concentrations.

## 3. The residual, and the finding that forced its shape

**Definition.** For every core species holding more moles than the simulator's own `atol`:

```
R = max_i | ln x_i(t_k) − ln x_i(t_{k−1}) | / ( ln t_k − ln t_{k−1} )
    = max_i | d ln x_i / d ln t |
```

Dimensionless, independent of the time unit, and independent of the solver's step size. Its
physical reading is the one that matters: since `|d ln x_i/dt| = 1/τ_i`, **`R = t / τ_min` is the
elapsed integration time expressed in units of the fastest chemistry still running.**

**On the lithium trajectory** R rises from ~1e-10, peaks at **15.5 at t = 7.89e-5 s**, then falls
off a cliff:

| t (s) | R | slowest species |
|---|---|---|
| 3.0e-15 | 2.830e-10 | Li+ |
| 9.75e-08 | 1.880e-02 | Li+ |
| 7.89e-05 | 1.548e+01 | Li |
| 1.00e-04 | 3.454e+00 | Li |
| 1.10e-04 | 5.730e-01 | Li |
| 1.24e-04 | 3.328e-02 | Li |
| 1.52e-04 | 1.355e-04 | Li |
| ≥ 1.84e-04 | ≤ 1.0e-07 | Li |

### 3.1 The residual alone is not a criterion — and this is the finding

Because `R = t/τ_min`, a **small R has two causes**: the system has converged, or the integration
has not yet reached its own chemistry. Both are real, and early in any run it is always the
second.

Measured, on the deck above:

```
max over species, NO arming   fires at t = 3.0000e-15 s   (R = 2.830e-10, step 2 of 497)
max over species, armed       fires at t = 1.8418e-04 s   (R = 1.036e-07, step 468 of 497)
```

**An unarmed criterion terminates 6e10 too early, at step 2 of 497, and reports success.** It is
not tunable away: excluding `R = 2.83e-10` needs a tolerance below the deck's own `rtol` of
`1e-8`, i.e. below the integrator's noise floor, at which point the criterion is chasing round-off.
Sweeping the tolerance down to `1e-10` makes the armed and unarmed variants coincide — on the
*wrong* answer.

So the criterion **arms**: it may fire only after `R` has reached 1, i.e. after the integration
has run longer than the fastest relaxation time in the system. On this trajectory `R` peaks at
15.5, so arming is comfortable, not marginal.

### 3.2 The tolerance is not load-bearing — which is what makes a default safe

Eight orders of magnitude of tolerance move the termination point by a factor of 2.9 in time and
by 25 of 497 solver steps:

| tolerance | fires at t (s) | step |
|---|---|---|
| 1e-2 | 1.3063e-04 | 447 |
| 1e-3 | 1.4267e-04 | 455 |
| 1e-4 | 1.5379e-04 | 461 |
| 1e-5 | 1.6750e-04 | 465 |
| **1e-6** | **1.8418e-04** | **468** |
| 1e-7 | 1.9752e-04 | 469 |
| 1e-8 | 2.2420e-04 | 470 |
| 1e-9 | 2.7756e-04 | 471 |
| 1e-10 | 3.8429e-04 | 472 |

`1e-6` sits in the middle of that plateau and is `100 × rtol` for the deck's own `rtol = 1e-8` —
tied to a number the deck already states, and two decades clear of the integrator's noise. The
composition at the `1e-2` and `1e-6` stops differs by 3.0e-4 relative in the slowest species.

### 3.3 All core species, not the electron alone — and the reason is counter-intuitive

**The electron settles first and the depleting neutral settles last.** The species carrying the
residual switches from Li+ to Li at `t = 3.6e-5 s` and stays Li for the rest of the run.

An electron-only criterion fires at `t = 1.1673e-04 s`, 1.58× earlier. At that instant:

| species | relative distance from its final value |
|---|---|
| e- | 3.667e-08 |
| Li | **5.502e-03** |

Five orders of magnitude apart. Watching the electron — the obvious choice for a plasma, and the
one a modeller would reach for — would declare convergence with the neutral still 0.55% out.

### 3.4 The norm: an honest caveat

L∞ (max over species) and L2 fire at **exactly the same step** on this trajectory, because one
species dominates the residual throughout. **The data does not separate them.** L∞ is chosen for
conservatism and because it names the offending species in the log, not because it was measured
to be better.

## 4. The three ways a run can end, and why only one of them is a steady state

Verifier 2 is the one that matters: a criterion that cannot fire is how a job runs for a week.
The implementation makes all three outcomes explicit, and end-to-end runs of each are recorded
below.

**(a) Converged.** Armed, then `R < tolerance` for `window` consecutive steps.

```
At time 2.2420e-04 s, reached steady state: residual 4.8449e-09 below tolerance
1.0000e-06 for 3 consecutive steps (slowest species: Li(3)).
```

(The first sub-tolerance step is 1.8418e-04 s as predicted offline; the window of 3 costs the two
following steps.)

**(b) Did not converge.** The backstop `terminationTime` fires first. `terminationSteadyState`
*requires* a companion `terminationTime` — refused at parse time otherwise — so this path always
exists and the run cannot integrate without limit.

```
At time 1.0000e-05 s, reached target termination time.
Warning: Terminated on the backstop terminationTime at 1.0270e-05 s WITHOUT reaching steady
state: residual 2.1340e+00 against tolerance 1.0000e-06 (0 consecutive steps below it, 3
needed; the integration passed the fastest relaxation time in the system). This result is NOT
a steady-state result.
```

It is a **warning, not an exception**, deliberately: `simulate()` is called once per model
enlargement during generation, and not reaching steady state mid-generation is normal. Raising
would kill every plasma run.

**(c) No flux, never armed.** The system carries no net rate anywhere, so the composition cannot
change and there is nothing left to integrate — but it never started, so it has demonstrated
nothing. This is the argon deck:

```
Warning: At time 0.0000e+00 s, terminating terminationSteadyState WITHOUT a steady state: the
system carries no flux, so the composition cannot change and there is nothing further to
integrate. NO STEADY STATE WAS DEMONSTRATED -- this model never started, which is not the same
result as a stationary ionisation balance, and must not be reported as one. There is no
residual to quote: the criterion was never evaluated, because a composition that cannot move
carries no information about whether it has settled.
```

`steady_state_reached` stays **False** on this path. The zero-rate condition itself is diagnosed
upstream by the solver; this message is only its termination consequence, so the two do not
double-report.

### 4.1 A frozen composition after a real transient still counts, and the latch is what says so

A first draft special-cased the no-flux state on *gross* rates: net zero with gross non-zero was
read as "production exactly balances loss" and reported as converged. Running the lithium deck at
tolerance `1e-30` showed why that was wrong — it fired anyway, at `t = 3.586e-3 s`, **bypassing
both the user's tolerance and the window**. A criterion that can fire without evaluating the
tolerance it was given is the shape this ticket exists to remove.

The arming latch subsumes the distinction with no extra machinery, and the gross-rate special
case was deleted. A run that went through a transient and then froze arrives at the no-flux state
*armed*, with a residual of exactly zero, which clears any tolerance and fills the window through
the ordinary path:

```
At time 3.5861e-03 s, reached steady state: residual 0.0000e+00 below tolerance 1.0000e-30
for 3 consecutive steps (slowest species: e-(1)).
```

Three consecutive steps, tolerance respected — not a bypass. A run that never armed gets outcome
(c).

**Corollary worth recording:** the lithium deck cannot be made to never-converge by tightening
the tolerance, because its composition reaches exactly zero net rate in floating point. To
demonstrate outcome (b) end-to-end the backstop has to land inside the transient (`1e-5 s`,
where `R ≈ 2.1`).

## 5. What this could not reach

- **The argon deck still has no argon chemistry.** Its steady state remains undemonstrated and
  undemonstrable until I-164's B2 is fixed (Voronov Z=18 and Badnell entries added to the two
  plasma libraries). This ticket makes that state *reportable as unresolved* rather than
  disguised as a converged 1e-3 s run; it does not resolve it.
- **The criterion is only wired to `plasmaReactor`.** The solver-side code in `base.pyx` is
  generic and would serve `simpleReactor` unchanged, but no verifier here exercises that.
- **One trajectory, one chemistry.** Every number above comes from a two-state ionisation balance
  with a single slow species. A deck with several comparably slow species, or with genuine
  oscillation, would test the L∞-versus-L2 question this data cannot settle.
