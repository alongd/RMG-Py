# I-246 — support envelope, observables, and conditioning threshold

**This file is written and committed BEFORE the sweep is run.** Its point is that the
classification rule and the conditioning threshold cannot have been chosen after seeing which
points failed. The commit that adds this file contains no sweep output; check its diff.

---

## 1. The support envelope

Everything the wall-transport formula is claimed to cover, and nothing else. A grid point outside
any of these is out of support and the reactor is expected to refuse rather than extrapolate.

| quantity | value / range | how it is fixed |
|---|---|---|
| gas | pure argon | deck (`docs/i194-ar5torr-plasma-lineage/input.py`) |
| buffer gas | none | deck |
| composition | Ar / Ar⁺ / e⁻ only | deck |
| pressure | **nominal 5 torr; swept 0.5 – 50 torr** | nominal from the deck; range chosen to bracket the crossover, both sides |
| gas temperature `Tgas` | 298.15 K, constant | deck; no gas-energy balance in scope |
| electron temperature `Te` | **swept 0.5 – 4.0 eV** (5800 – 46 400 K) | engineering range for an argon positive column; the deck's 3 eV sits inside it. **Prescribed, never solved for** — closing it is the next milestone |
| geometry | right circular cylinder, **R = 0.05 m, L = 0.30 m** nominal; R swept 0.005 – 0.20 m | declared chamber, not calibrated |
| diffusion length `Λ` | `1/Λ² = (2.405/R)² + (π/L)²` | lowest-order diffusion eigenmode of a finite cylinder (Schottky positive-column standard). Λ(nominal) = 2.031 cm |
| wall area / electrode area | not separately resolved | the model is a diffusion-limited ambipolar loss to a single enclosing surface; there is no sheath-area or electrode-area factor, and none is fitted. Stated as a limitation, not a parameter |
| wall material | unresolved **on purpose** | see γ below: for a noble gas the recycling coefficient is material-independent to the accuracy claimed here |
| electrical boundary | **floating wall, zero net wall current** | this is the declared current-balance rule; it is what one common ν for all charged species enforces |
| magnetic field | **none** (unmagnetised) | declared; a magnetised wall loss needs an anisotropic D and is out of scope |
| flow / residence time | **none — closed batch**, infinite residence time | property of `PlasmaReactor`; there is no inlet, outlet or makeup stream |
| ionisation degree | **α = n_e/n_neutral ≤ 1e-3** | validity ceiling of the ion-*neutral* mobility model. Above it Coulomb ion-ion momentum transfer is no longer negligible against the Ar⁺/Ar resonant charge-transfer cross section (~5e-19 m²), and the 1/n_neutral scaling is wrong. The reactor **hard-fails** above it |

### Transport data — one independent source each

| parameter | value | source |
|---|---|---|
| reduced zero-field mobility `μ₀(Ar⁺ in Ar)` | **1.535 × 10⁻⁴ m²/(V·s)**, i.e. 1.535 cm²/(V·s), **±3 %** | Ellis, McDaniel & Albritton, *Transport properties of gaseous ions over a wide energy range*, At. Data Nucl. Data Tables **17** (1976) 177. Low-E/N limit. The ±3 % is the compilation's own stated accuracy |
| reference density `N₀` | 2.6867811 × 10²⁵ m⁻³ | CODATA Loschmidt constant (273.15 K, 101 325 Pa). A **unit convention**: it is the density the tabulated μ₀ is normalised to |
| ambipolar diffusivity `D_a` | `μ_i · k_B T_e / e` | `D_a = D_i(1 + T_e/T_i)` in the `T_e ≫ T_i` limit, via the Einstein relation. **Drops a factor `(1 + T_i/T_e)`** — at Te/Tgas = 117 that is a 0.85 % underestimate of ν_wall, carried inside the interval below, not corrected away |
| ionisation rate `k_iz(Te)` | Voronov fit, Ar (Z=18, N=18): A = 5.99e-8 cm³/(molecule·s), P = 1, X = 0.136, K = 0.26, ΔE = 15.8 eV | Voronov, At. Data Nucl. Data Tables **65** (1997) 1. Read from the repository's own `voronov.yaml` through the compiled `VoronovEIArrhenius`, not re-implemented here |
| recycling coefficient `γ` | **1.0** | Ar⁺ incident on any practical wall is Auger-neutralised with probability ≈ 1, and argon does not chemisorb, so the neutralised atom re-enters the gas. Material-independent for a noble gas at these energies. `γ` remains a declared input in [0,1]; **1.0 is a sourced physical statement, not a fitted value** |
| external pair source `S_ext` | **1.3 × 10⁴ – 6.6 × 10⁴ m⁻³ s⁻¹** at 5 torr | cosmic-ray / ambient-radiation background ionisation, 2–10 ion pairs cm⁻³ s⁻¹ at 1 atm, scaled by gas density. **Carried as the interval it is and not narrowed**; the absolute sub-threshold electron density inherits the full factor-5 width of this interval |

### Total uncertainty carried on ν_wall

`ν_wall = μ₀ N₀ (k_B T_e/e) / (n_neutral Λ²)`. Independent contributions: μ₀ ±3 %, the dropped
`(1 + T_i/T_e)` −0.85 %. **ν_wall is therefore reported as a central value with a −0.85 %/+3 %
interval**, and every quantity derived from it inherits it. No term of this interval was obtained
by comparing anything to an electron density.

---

## 2. Predeclared observables

The opening probe measured that the same physical event moves `n_e` by 4.3 ppm and the neutral
fraction by 5.1e5 — four orders apart on one event — because the two are normalised to
populations at opposite ends of their range. So **every claim below is labelled with the
observable it is stated in**, and claims stated in different observables are never compared.

| # | observable | units | role |
|---|---|---|---|
| O1 | `ν_wall`, `ν_ion = k_iz(Te)·n_Ar`, `ν_RR`, `Δν = ν_ion − ν_wall` | s⁻¹ | **primary.** Determined by geometry and transport alone; this is what the operator is judged on |
| O2 | ionisation degree `α = n_e / n_neutral` | — | **primary.** The quantity the wall actually acts on, and the one the validity ceiling is stated in |
| O3 | neutral mole fraction `x_Ar` | — | **primary.** The heavy-composition observable; volume-independent |
| O4 | absolute electron density `n_e` | m⁻³ | **reported, never an acceptance criterion.** No threshold, fixture or constant in the implementation or the tests is stated in it |
| O5 | net gas charge `Σ z_j N_j`, and gas+wall heavy inventory | mol | **conservation.** Must not drift |
| O6 | `C_wall = |dln n_e/dln ν_wall|`, `C_Te = |dln n_e/dln Te|`, `C_Ar = |dln n_e/dln n_Ar|` | — | **conditioning** |

---

## 3. Predeclared classification and conditioning threshold

Fixed here, before any grid point is evaluated.

A grid point is classified by these rules, applied in order:

1. **out-of-support** — if the state violates any envelope row in §1 (in practice: α > 1e-3).
2. **near-threshold / ill-conditioned** — if **either**
   - `|Δν| / max(ν_ion, ν_wall) < 1e-2`  (the balance is a cancellation of two numbers agreeing
     to better than 1 %), **or**
   - any of `C_wall`, `C_Te`, `C_Ar` exceeds **1e2**.
3. **volume-growth** — `Δν > 0` and not near-threshold.
4. **extinction** — `Δν < 0` and not near-threshold.

**Why 1e2 and 1e-2, chosen now and not later.** `C_wall = (1−f)/f` in the ionisation fraction
`f`, so `C_wall > 1e2` is exactly `f < ~1 %`: the threshold states that a 1 % relative error in
the wall coefficient must not move `n_e` by more than a factor e. It is a statement about
admissible error propagation, fixed by the accuracy of the transport data in §1 (±3 % on μ₀),
not by which points turn out to fail. The 1 % on `|Δν|` is the same statement on the balance
itself: ν_ion and ν_wall each carry ≳1 % uncertainty, so a difference smaller than that is not
resolvable by this model at all and must not be reported as a sign.

**Finite differences across the boundary are forbidden.** The conditioning numbers are taken by
**bracketing/continuation along each branch**, never by a one-sided difference that steps across
the extinction boundary. Where a branch terminates, that is reported as a branch terminus, not as
a derivative.

**A single Newton solve returning a positive density is not evidence.** Each reported steady state
must be reached from at least two distinct initial conditions, or shown to be the unique root of
the scalar balance on its branch.

---

## 4. What is NOT claimed, stated before the results exist

Self-consistent `Te`; absorbed power; discharge current; plasma potential; sheath structure;
absolute steady-state `n_e`; agreement with any glow-discharge measurement. The model is five to
seven orders of magnitude from a real argon glow discharge and **the gap is not to be closed by
choosing a wall coefficient**. If it cannot be closed with independently sourced parameters, that
is the result.
