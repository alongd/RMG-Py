# I-274 — Electron energy balance: Te solved from absorbed power, n_e from it (M8-B)

> **ENGINEERING INTERMEDIATE.** The discharge closure is a *specified absorbed power*: a
> global-model closure that says nothing about how a real source couples its power. A DC glow
> sustained by secondary emission would need a discharge-current + sheath circuit closure,
> which is not built. The EEDF is Maxwellian. Every P_abs value below is **illustrative, not
> measured**, and none was chosen to match a density.

Engine: branch `electron-energy-balance` (base `89cb58890`). Database:
`RMG-database-plasma @ 0d9c5bc24`, read only. Runs, logs (stdout + stderr per arm), scripts and
report-ready artifacts: `/home/alon/runs/i274-energy-balance-20260924-002658/`
(`results.md`, `budget_0p5W.md`, `summary.md`, `fig_ne_te_vs_power.png|svg`, `figure_caption.md`;
`lxcat/` holds the LXCat source and fit, `arms/` the production LXCat arms and `arms-LL/` the L&L sensitivity set).

## What was built

`PlasmaReactor(electron_energy_balance=...)` / input keyword `electronEnergyBalance`. With it, Te
is one extra DAE state after the core species:

    d(3/2 N_e k_B Te)/dt = P_abs - Q_inelastic - Q_elastic - Q_wall - Q_flow

| term | form | source |
|---|---|---|
| P_abs | absorbedPower / chamber volume (cylinder/sphere geometry, or `chamberVolume`) | deck (illustrative) |
| Q_inelastic | sum_j r_j eps_j, eps_j **declared per library entry** (`'PlasmaArgon:86': (15.76, 'eV')`); electron-consuming reactions also remove 3/2 k_B Te per electron | declaration; thermo dH logged as a cross-check only |
| Q_elastic | 3 (m_e/M) K_m(Te) n_Ar n_e k_B (Te - Tg), K_m = A Te^n exp(b ln²Te + c ln³Te); production A = 1.8017e-14 m³/s, n = 1.54096, b = −0.019608, c = −0.057221 | **LXCat Phelps** e+Ar EFFECTIVE σ_m (www.lxcat.net, Phelps database, retrieved 2026-09-23; Yamabe, Buckman & Phelps, PRA 27, 1345 (1983), rev. 1997), Maxwell-averaged and fitted over 0.5–3 eV (finding 2). Sensitivity case: Lieberman & Lichtenberg 2005 Table 3.3 (2.336e-14 Te^1.609 exp(0.0618 ln² − 0.1171 ln³)) |
| Q_wall | 2 k_B Te per lost electron + k_B Te (1/2 + 1/2 ln(M/2π m_e)) = 5.18 k_B Te per lost Ar⁺ (Bohm presheath + floating sheath) | L&L 2005 §10.2; built from the **same** wall-loss array the species rows apply |
| Q_flow | 0: the reactor is a closed batch | — |

**Per-reaction energies are declared, not inferred (PM correction, 2026-09-24).** PlasmaArgon 91,
`Ars + e- => Ar + e-`, is a collapsed proxy for m→r mixing followed by radiation. Its thermo dH is
-11.548 eV, but the electron actually pays the 3P2→3P1 gap. The FULL-deck declaration (NIST ASD
levels):

| entry | reaction | declared (eV) | thermo dH cross-check (eV) |
|---|---|---|---|
| 86 | Ar + e → Ar⁺ + 2e | +15.7596 | +15.7597 |
| 87 | Ar + e → Ars + e | +11.5484 | +11.5484 |
| 88 | Ars + e → Ar⁺ + 2e | +4.2113 | +4.2114 |
| 89 | Ars + e → Ar + e (superelastic) | −11.5484 | −11.5484 |
| 90 | Ars + Ars → Ar + Ar⁺ + e (pooling) | −7.3371, **credited to the ejected electron** | −7.3370 |
| 91 | Ars + e → Ar + e (m→r mixing proxy) | **+0.0752** (3P1 − 3P2) | −11.5484 |
| RR:1 | Ar⁺ + e → Ar (radiative) | 0, plus 3/2 k_B Te automatically | −15.7597 |

An electron reaction in the core without a declaration is refused, and the refusal names its
thermo value. Te-dependent rate coefficients are re-evaluated from the Te state. The Jacobian in
this mode is a central difference of the residual (the analytic one assumes a fixed Te). In this
mode only, the charged species' atol is 1e-12× smaller, so N_e stays resolved down to the
source floor, because the Te row divides by N_e. A power budget is latched at accepted states,
with dU_e/dt taken from the solver's own derivative. `discharge_state()` returns
`self-sustained` when ν_ionisation ≥ (1 − 1e-3) ν_loss, `source-supported` when the external
source makes up the deficit, and `extinct` otherwise (round 3; it was a ≥ ½ cutoff). The M8-A field `wall_ion_energy_flux`,
left `declared-absent`, is now filled with the sheath term (`available-floating-wall-sheath-model`).

## Verifier results (contract §11 clause numbers)

| # | item | result | status |
|---|---|---|---|
| 1 | **Cl. 1 + 9**: Te solved; two initial Te converge | Te0 = 0.6, 1.0, 2.0 eV → 0.89881340 eV all three; n_e 1.17344391e16 m⁻³ agree to 1e-9 | **closed** |
| 2 | Te vs particle-balance Te* (tolerance 0.5 meV) | 0.5 W: 0.89881 eV vs Te*_FULL 0.89886. Across 0.001–10 W Te follows Te*(n_e) of report7: 0.89979 @2.3e14 (0.8998), 0.89960 @1.2e15 (0.8996), 0.89881 @1.2e16 (0.8989), 0.89903 @1.2e17 (0.8990). Energy equation does **not** overturn particle balance. | **closed** |
| 3 | **Cl. 2 + 3**: n_e vs P_abs, no wall retune; hand check | 13 powers, 0.001–10 W: n_e 2.33e13 → 2.34e17 m⁻³, log-log slope **1.0006**. Global formula at 0.5 W (LXCat K_el): ν_w 55.470 s⁻¹, E_c 2025.5 eV, E_e 1.798 eV, E_i 4.655 eV → n_e(hand) 1.1734e16 vs solver 1.1734e16 (2.3e-6, RMG's Ar mass) | **closed** |
| 4 | **Cl. 4**: radius ±50 % at 0.5 W | R 2.5 / 5 / 7.5 cm → Te 0.9755 / 0.8988 / 0.8610 eV (smaller chamber, larger wall loss, higher Te*); n_e 2.93e16 / 1.17e16 / 6.35e15 | **closed** |
| 5 | **Cl. 5**: extinction vs sustained | P_abs = 0 (LXCat K_el, default atol): Te within 1 % of Tg by 1.9 ms; the run **exits 0 with termination = 'extinct' at t = 1.28 s**: the extinct condition held from 2.1 ms for 1.28 s against the required 1.26 s (2 × the slower of 1/ν_loss and the energy-relaxation time; round 3, finding 4). No refusal in stderr. Every tested P_abs, 0.001–10 W, is `self-sustained` and never terminates; the cold-start test does not stop while heating. Energy-off FULL stays byte-identical to `mainline-FULL` (`round3-FULL-off/`) | **closed** |
| 6 | **Cl. 6 + 10**: ±20 % sensitivity at 0.5 W | ionisation ×0.8/×1.2: Te +11.6/−9.3 meV, n_e −6.0/+5.0 %. D_a ×0.8/×1.2: Te −11.6/+9.7 meV, n_e +6.4/−5.3 %. P_abs ×0.8/×1.2: Te +0.07/−0.04 meV, n_e −20.0/+20.0 %. K_el ×0.8/×1.2: Te −0.04/+0.04 meV, n_e +18.8/−13.7 %. **Joint 2⁴ corner sweep** (PM ruling 3, LXCat K_el baseline, 16 runs, all steady): Te **0.8781–0.9205 eV** (−20.7/+21.7 meV), n_e **7.25e15–1.88e16 m⁻³** (×0.62/×1.60 of 1.17e16). Each corner is within 1.3 % of the product of the one-at-a-time factors. The same sweep on L&L (`arms-LL/`) gives ×0.62/×1.60 about 9.55e15 (`results.md`) | Cl. 10 **closed for these four inputs**. Not covered: the Maxwellian-EEDF systematic |
| 7 | **Cl. 7**: budget closure | 0.5 W steady state: \|P − ΣQ − dU/dt\|/P = **1.2e-12**, dU/dt/P = 1.5e-12. Table in `budget_0p5W.md`: elastic 79.0 %, excitation (87) 19.9 %, ionisation 0.76 %, wall ions 0.23 %, wall electrons 0.09 % (L&L: elastic 82.9 %, closure 6.9e-14) | **closed** |
| 8 | **Cl. 8**: wall energy uses the particle flux | Q_wall,e / latched `wall_electron_energy_flux` = 1.000000000000000. Both are built from the same `wall_loss_rates` evaluation (unit test pins it at 1e-15) | **closed** |
| 9 | Regression, energy off | FULL arm vs `mainline-FULL`: final state, rates, wall flux, V, t and every trace snapshot **bit-identical**, re-checked after round 3 (`round3-FULL-off/`). nu_wall re-measured **15.954491 s⁻¹** (`round3-nu-wall/`). The five report7 arms refused on an Ars negative are refused again with report7's own message and a byte-identical trace (`round3-report7/`): with Te prescribed the wall check is the pre-I-274 code. Tests `test/rmgpy/solver/` + `inputTest.py`: **398 passed, 1 skipped** before I-274 → **467 passed, 1 skipped** after round 3 | **closed** |
| 10 | Red-first unit tests | `test/rmgpy/solver/plasmaEnergyBalanceTest.py`: 23 tests run red before implementation (`red/stdout.log`), including a hand value for each loss term and the energy-conservation row. The declared-energy and input-keyword tests were added with the PM correction and written alongside the code, not run red first | closed, with that exception |
| 11 | Both streams captured | every arm dir has `stdout.log` + `stderr.log`, as do the suites, analysis and remeasure | **closed** |

Clauses 11–12 (experimental observable) are out of scope: no observable is ratified.

## Entry 91: Te with the mixing proxy on and off (0.5 W)

| case | Te (eV) | n_e (m⁻³) |
|---|---|---|
| 91 declared +0.0752 eV (production) | 0.89881 | 1.1734e16 |
| 91 removed | 0.89721 (−1.60 meV) | 1.1845e16 (+0.9 %) |
| 91 at its enthalpy −11.548 eV (ruled out; counterfactual) | 0.89879 (−0.02 meV) | 1.2844e16 (+9.5 %) |

Production (LXCat K_el) Te is 0.898813 eV. On the L&L sensitivity case the same three cases gave
0.89888 / 0.89769 / 0.89885 eV and n_e +0.6 % / +6.5 %. The PM expected that crediting 91 with its
enthalpy would move Te. **It does not.** Te is pinned by particle balance, so the spurious heating
shows up only in n_e: +9.5 % at 0.5 W on LXCat, +6.5 % on L&L. The elastic loss is 79 % of P_abs
on LXCat and 83 % on L&L, so on LXCat the same spurious heating is a larger share of the loss it
offsets, and n_e moves more. Removing 91 moves Te by −1.6 meV, as report7 found under prescribed Te
(−1.3 meV).

## Findings

1. **Te is set by particle balance, n_e by power, exactly as the contract predicts.** Te spans
   1.1 meV (0.89872–0.89984 eV) over four decades of P_abs, tracking Te*(n_e), and n_e is linear in P_abs. The ill-posed n_e of the
   prescribed-Te model (21× over 0.9 meV) becomes well-posed. The price is that n_e now carries
   the full uncertainty of P_abs (d ln n_e/d ln P = 1.00) and ~0.25 of the ionisation rate and D_a.
2. **The elastic rate was 29 % too high; production now uses LXCat (PM ruling 1).** Elastic loss
   dominates the energy cost at 5 torr (2025 eV per lost pair, 79 % of the power), so the absolute
   n_e depends most on K_el at 0.9 eV. The Maxwell average of the LXCat Phelps EFFECTIVE σ_m
   (fetched by the PM; `lxcat/SOURCE.md`) is K_m(0.9 eV) = **1.5321e-14 m³/s**, reproduced
   independently here (`scripts/lxcat_fit.py`: own parser and quadrature, 1.53213e-14 against the PM's
   `km.py` 1.53212e-14). The L&L fit gives 1.9734e-14 there, so LXCat/L&L = **0.776**: 0.758–0.781
   over 0.5–1 eV and 0.70 at 3 eV. That is over the 10 % threshold, so production switched. The
   EFFECTIVE set adds inelastic momentum transfer above 11.5 eV, 2.4e-4 of the integral at 0.9 eV.
   The engine takes the four-coefficient form, fitted by least squares on ln K over 0.5–3 eV
   (`lxcat/fit.json`, written by this session): max error 0.7 % over the range, ≤ 0.09 % over the
   operating 0.878–0.921 eV. **n_e rises 22.8 % at fixed P_abs** (22.7–22.9 % across the sweep),
   matching the expected d ln n_e/d ln K_el × ln 0.776 = −0.83 × −0.254 = +21 to +24 %. Te moves by
   ≤ 0.07 meV. L&L is kept as the sensitivity case (`arms-LL/`, section in `results.md`). Limit:
   the four-coefficient form cannot follow the Ramsauer minimum, so below 0.5 eV the fit is
   extrapolated. At Tg it gives 8.1e-16 against the true 2.3e-15 (L&L: 4.6e-14, 20× too high). This
   matters only for P = 0 relaxation, where it sets Tg − Te (finding 4). A 0.025–3 eV fit would be
   6 % off at 0.9 eV, so it was not used.
3. **RMG carries argon at 39.8775 g/mol** (`rmgpy.molecule.element`), not the standard 39.948
   (Ar-40 is 39.962). Through m_e/M this shifts elastic loss, and hence n_e, by 0.15 %. The hand
   check with the standard mass is off by exactly that (1.47e-3); with RMG's mass it closes to
   2.5e-6. Not changed here (out of scope). Worth a ticket.
4. **Extinction run: closed as a terminal state.** At P_abs = 0 the state is labelled `extinct`
   within 0.2 ms (L&L; 1.9 ms on LXCat). Before the fix the FULL run exited 1 at 10 ms: once Te ≈ Tg, Ars decays to ~0 and
   its accepted value is solver noise (−2e-32 mol even at atol 1e-30), which the M8-A
   `check_wall_support` refused as a negative population. Ruling: a **neutral** population in
   [−atol_j, 0), with atol_j that species' own absolute tolerance, is clamped to exactly 0 in the
   published accepted state; anything more negative, any non-finite value and any **charged**
   negative stay refused. The clamp acts on the state the model reads and publishes at the
   accepted step, not on DASPK's own history (which carries the same noise inside its tolerance).
   Six tests in `plasmaWallTest.py`; the two clamp tests ran red first (`clamp-red/`). On the L&L
   K_el the P0 arm (default atol 1e-16) then runs to steady state at the hand floor
   6.6e4/1.5856 = 4.16e4 m⁻³ (`arms-LL/P0`). The pre-clamp run is kept in `arms/P0-preclamp/`.

   **On the LXCat production rate P0 was refused again, correctly under the ruling.** The Ars
   accepted value overshoots its own atol: −1.02e-16 at atol 1e-16 (t = 2.65 s), −1.54e-20 at 1e-20
   (1.62 s), −2.73e-30 at 1e-30 (0.84 s), i.e. 1.02–2.7× atol (`arms/P0-atol16-refused`,
   `arms/P0-atol20-refused`, `arms/P0-atol30`). DASPK bounds only the weighted RMS error over the
   neq = 5 components, so one component can reach about √5 × its own atol.

   **PM ruling: keep the 1× bound; extinction is a terminal state** (round 2, reworked in round 3).
   Widening the bound would tune a tolerance to make one run pass. In energy-balance mode only, the
   extinct condition — `discharge_state()` `extinct`, ν_iz/ν_loss < 1e-6 and |Te − Tg| ≤ 1 K
   (`PLASMA_EXTINCT_RATIO`, `PLASMA_EXTINCT_TE_BAND_K`) — must hold at every accepted step over
   **physical time** for `PLASMA_EXTINCT_PERSIST_MULTIPLE` = 2 times the slower of 1/ν_loss and the
   elastic energy-relaxation time 1/(3 Σ (m_e/M) K_m n). Calls at a time not after the last counted
   one never count, and a miss restarts the clock. It is allowed only at P_abs = 0 or once the
   discharge has been `self-sustained`. The round-2 rule counted 10 accepted steps, and the round-3
   red run showed it stopping a powered cold start (Te0 = Tg + 0.5 K) during its heating transient.
   On success the reactor records `energy_terminal` (termination, start, end, duration, required
   persistence, Te, Tg, n_e, the three frequencies and the state vector). `ReactionSystem.simulate`
   asks the `terminal_state()` hook after each accepted step (base default None) and stops normally.
   On the production deck P0 (LXCat, default atol) exits 0 as `extinct` at **t = 1.28 s**. The
   condition held from 2.1 ms, against a required 1.26 s (2/ν_loss at Tg; the energy-relaxation time
   there is ~0.2 ms). n_e at the stop is 1.28e15 m⁻³, and there is no refusal. The decay to the
   4.16e4 m⁻³ floor, where the state would read `source-supported`, is not integrated; the figure
   marks that floor as a hand value. Tg − Te in the tail is 0.49 K (LXCat fit) against 0.011 K
   (L&L): wall-lost electrons carry 2kTe against a mean 3/2 kTe, which cools the rest, and elastic
   heating from the gas restores them, so the deficit scales as 1/K_el(Tg). The ratio of the two
   rates at Tg (57) predicts 0.63 K.

   **report7 refusals: unchanged.** In round 2 the clamp also ran with Te prescribed, and the five
   f0d9c arms refused on an Ars negative (FULL-91OFF, SWEEP-0.80, SWEEP-0.80-91OFF, SWEEP-0.85,
   SWEEP-0.85-91OFF) instead failed one check later at the same step, on a negative electron amount
   (`clamp/report7/`). Round 3 gates the clamp to energy mode at accepted steps, so with Te prescribed
   the check is the pre-I-274 code. All five are again refused with report7's own Ars message at
   the same step, with byte-identical traces (`round3-report7/`).

5. **This closure has no minimum sustaining power** over the tested range. Every tested P_abs,
   0.001–10 W, sustains a discharge at Te*
   with n_e ∝ P_abs, down to 0.001 W (2.3e13 m⁻³), until n_e would meet the source floor near
   ~2e-12 W. A real discharge's extinction threshold (coupling-mode change, sheath/circuit
   limits) is outside a specified-power model. Contract §10's caveat bites here.
6. **Numerics, in declared settings only.** (a) Charged-species atol ×1e-12 in energy mode (engine).
   (b) P_abs ≤ 0.02 W arms use deck `atol = 1e-22`; at 0.05 and 0.5 W that setting changes Te and
   n_e by < 1e-9. P_abs = 0.001 and 0.01 W start from n_e = 1e13/1e14 m⁻³ rather than 1e16: from 1e16
   they must decay ~50× first, and hit finding 4 on the way. The first attempts are kept in
   `failed-first-attempt/` and `failed-second-attempt/`. Steady state (terminationSteadyState 1e-6)
   is reached at ~0.1 s, against ~26 s under prescribed Te: the power closure turns the marginal
   n_e mode into a damped one.
7. **Sheath model scope.** Floating wall, collisionless Bohm sheath, single ion. At 5 torr the bulk
   is collisional, so the ion's ambipolar-field energy loss in the bulk is not counted. The wall
   terms together are 0.32 % of P_abs here, so this is immaterial at 5 torr. It would matter at
   low pressure.

8. **What moves Te and what moves n_e (corner sweep).** Te is set by the ratio of ionisation to
   ambipolar loss. Moving ionisation and D_a in opposite directions shifts Te by ±21 meV, and moving
   them together cancels (≤0.3 meV). K_el and P_abs barely touch Te (≤0.1 meV) but carry n_e:
   d ln n_e/d ln P_abs = +1.00 and d ln n_e/d ln K_el ≈ −0.79 (LXCat; −0.83 on L&L). The elastic
   term is 79 % of the power, so it enters E_c with that weight. Ionisation and D_a each move n_e by
   only about ∓0.28 per unit log. The response is nearly separable: every corner is within 1.3 % of
   the product of the one-at-a-time factors. The n_e envelope, ×0.62 to ×1.60 on either K_el
   baseline, is therefore dominated by P_abs and K_el.

## Round 3 (Codex round 19 on aa3dbf83a; PM brief `i274-round3.md`)

The physics bookkeeping held. Seven boundary fixes, each with tests run red on aa3dbf83a
(`round3-red/`) and green after (`round3-green/`):

1. **Energy-off regression.** The neutral-noise clamp runs only with the energy balance and only
   at accepted steps (`check_wall_support(y, accepted=True)`). The initial composition is refused,
   never repaired. The energy-off test that blessed the clamp is replaced by one asserting the old
   refusal. Energy-off FULL is byte-identical and the report7 arms are refused as before (row 9).
2. **Extinction by physical time**, eligible only at P_abs = 0 or after self-sustainment (finding 4).
   Tested: identical and backward times never count, a miss restarts the clock, a powered run that
   was never sustained cannot stop, and a powered cold start runs to its backstop.
3. **Three discharge states** from `classify_discharge(ν_iz, ν_source, ν_loss)`. Every boundary is
   tested, including the old ½ cutoff, which now reads `extinct`. The P = 0 toy decays as `extinct`
   and settles on the source-held floor as `source-supported`.
4. **Zero electrons.** Energy mode refuses N_e ≤ 0 at initialisation ("the electron temperature is
   undefined"), even with a source declared. The 1e-300 floor in the Te row is now a named
   unreachable-state error. The Jacobian goes one-sided within one step of zero, so it never
   evaluates a negative amount; a spy test on the residual pins that. No run in either set reached
   the unreachable error.
5. **P_abs is a total power.** The reactor is a fixed-inventory image of the chamber, scaled once
   by V_ref = N_heavy(t0) R Tg / P, and heated with the constant P_abs V_ref / V_chamber (W), not
   (P_abs/V_chamber) V(state). **Electron pdV work under the constant-pressure EOS is neglected.**
   Bound: V(state)/V_ref − 1 = n_e k Te / P, at most **5.1e-5** over every accepted state of every
   production arm (10 W, `results.md`). The round-2 → round-3 shift in n_e is the same size (5.05e-5
   at 10 W, 2.5e-6 at 0.5 W). The hand check now uses the same definition and closes to 2.3e-6.
6. **Elastic partners complete.** Every core neutral declares an elastic rate or
   `{'ignore': '<reason>'}`, and an undeclared neutral is refused. The production deck declares Ars
   ignored: it is ~5e-8 of the gas, so its elastic loss is < 1e-7 of Ar's.
7. **Initial Te** below Tg/2, the evaluation floor, is refused.

**Physics declared.**
- **Ion wall energy.** The 5.18 kTe ion wall energy is a collisionless floating sheath, while at
  5 torr the sheath is collisional, so the true value lies between 0 and that. A 0.5 W rerun with
  the term at 0 (`arms/ION0`, `ION_SHEATH_SCALE=0`) gives ΔTe = −0.001 meV and Δn_e = +0.23 %, the
  size of the term (0.23 % of P_abs). Immaterial at 5 torr.
- **Source electrons.** Electrons from the external ionisation source enter with **zero energy**.
  This is a declaration of the model, in the code and the docs.

**Tests added:**
- pooling (entry 90's shape, a declared −7.34 eV credit, no 3/2 kTe charge);
- the production deck, added as `examples/rmg/plasma_argon_energy_balance/input.py`, with its LXCat
  coefficients, the Ars ignore and the −7.3371 / +0.0752 / −11.5484 / +15.7596 eV declarations
  pinned in `inputTest.py`;
- every Jacobian column (electron, ion, ground, metastable, Te) and the mass matrix at cj = 3.7e3
  against finite differences of the residual in y and in dy/dt;
- a budget closure recomputed from the rates, the declared energies, the hand elastic fit and the
  wall scratch. Its sabotage twin, the engine's ion sheath term zeroed, fails it. The tolerance is
  5e-6 of P: the hand/engine eV-constant vintage differs by 1.0e-6 in Te, which the elastic term
  (∝ ~Te^1.6) carries into the sum.

All arms of both sets were rerun on the round-3 engine (`arms/`, `arms-LL/`; the round-2 arms are
kept in `arms-r2/`, `arms-LL-r2/`), and the results, budget, figure, caption, summary and both
corner sweeps were regenerated. At the precision quoted above, only the initial-Te row's n_e
moved.

## What is and isn't predicted

See `summary.md` in the runs directory (five lines, report-ready). In short: Te at 5 torr is a
model result, within the Maxwellian-EEDF bias i258 already flagged. n_e is predicted
**conditionally on P_abs** and is not a prediction of any measured discharge until the
absorbed power and discharge mode are ratified and an observable is chosen (clauses 11–12).

## Reproduce

    R=/home/alon/runs/i274-energy-balance-20260924-002658
    $R/scripts/run_arm.sh P0.5 --power 0.5          # one arm, LXCat K_el (--elastic ll: L&L)
    python $R/scripts/lxcat_fit.py 0.5 3.0           # K_m(Te) from lxcat/Phelps_Ar_effective.txt -> lxcat/fit.json
    python $R/scripts/analyse.py                    # tables, figure, summary
    PYTHONPATH=$PWD pytest --no-cov test/rmgpy/solver/plasmaEnergyBalanceTest.py
