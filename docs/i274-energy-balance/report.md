# I-274 — Electron energy balance: Te solved from absorbed power, n_e from it (M8-B)

> **ENGINEERING INTERMEDIATE.** The discharge closure is a *specified absorbed power*: a
> global-model closure that says nothing about how a real source couples its power. A DC glow
> sustained by secondary emission would need a discharge-current + sheath circuit closure,
> which is not built. The EEDF is Maxwellian. Every P_abs value below is **illustrative, not
> measured**, and none was chosen to match a density.

Engine: branch `electron-energy-balance` (base `89cb58890`). Database:
`RMG-database-plasma @ 0d9c5bc24`, read only. Runs, logs (stdout + stderr per arm), scripts and
report-ready artifacts: `/home/alon/runs/i274-energy-balance-20260924-002658/`
(`results.md`, `budget_0p5W.md`, `summary.md`, `fig_ne_te_vs_power.png|svg`, `figure_caption.md`).

## What was built

`PlasmaReactor(electron_energy_balance=...)` / input keyword `electronEnergyBalance`. With it, Te
is one extra DAE state after the core species:

    d(3/2 N_e k_B Te)/dt = P_abs - Q_inelastic - Q_elastic - Q_wall - Q_flow

| term | form | source |
|---|---|---|
| P_abs | absorbedPower / chamber volume (cylinder/sphere geometry, or `chamberVolume`) | deck (illustrative) |
| Q_inelastic | sum_j r_j eps_j, eps_j **declared per library entry** (`'PlasmaArgon:86': (15.76, 'eV')`); electron-consuming reactions also remove 3/2 k_B Te per electron | declaration; thermo dH logged as a cross-check only |
| Q_elastic | 3 (m_e/M) K_m(Te) n_Ar n_e k_B (Te - Tg), K_m = 2.336e-14 Te^1.609 exp(0.0618 ln²Te - 0.1171 ln³Te) m³/s | Lieberman & Lichtenberg 2005, Table 3.3 (argon elastic) |
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
with dU_e/dt taken from the solver's own derivative. `discharge_state()` returns `sustained`
when ν_ionisation ≥ ½ ν_loss and `extinct` otherwise. The M8-A field `wall_ion_energy_flux`,
left `declared-absent`, is now filled with the sheath term (`available-floating-wall-sheath-model`).

## Verifier results (contract §11 clause numbers)

| # | item | result | status |
|---|---|---|---|
| 1 | **Cl. 1 + 9**: Te solved; two initial Te converge | Te0 = 0.6, 1.0, 2.0 eV → 0.89887502 eV all three; n_e 9.5538224e15 m⁻³ agree to 1e-9 | **closed** |
| 2 | Te vs particle-balance Te* (tolerance 0.5 meV) | 0.5 W: 0.89888 eV vs Te*_FULL 0.89886. Across 0.001–10 W Te follows Te*(n_e) of report7: 0.89980 @1.9e14 (0.8998), 0.89964 @9.5e14 (0.8996), 0.89888 @9.6e15 (0.8989), 0.89898 @9.5e16 (0.8990). Energy equation does **not** overturn particle balance. | **closed** |
| 3 | **Cl. 2 + 3**: n_e vs P_abs, no wall retune; hand check | 13 powers, 0.001–10 W: n_e 1.90e13 → 1.91e17 m⁻³, log-log slope **1.0006**. Global formula at 0.5 W: ν_w 55.474 s⁻¹, E_c 2489 eV, E_e 1.798 eV, E_i 4.655 eV → n_e(hand) 9.5538e15 vs solver 9.5538e15 (2.5e-6) | **closed** |
| 4 | **Cl. 4**: radius ±50 % at 0.5 W | R 2.5 / 5 / 7.5 cm → Te* 0.9754 / 0.8989 / 0.8611 eV (smaller chamber, larger wall loss, higher Te*); n_e 2.49e16 / 9.55e15 / 5.09e15 | **closed** |
| 5 | **Cl. 5**: extinction vs sustained | P_abs = 0: Te reaches Tg within 0.2 ms; state `extinct` (ν_iz/ν_loss = 1e-48). Toy-model test integrates P = 0 to the source floor. Every P_abs > 0 sustains. After PM ruling (2) the FULL run completes: exit 0 at t = 28.2 s (steady), Te = 298.14 K, n_e = 4.1626e4 m⁻³ vs hand floor 4.1624e4 (3.8e-5) | **closed** (finding 4) |
| 6 | **Cl. 6 + 10**: ±20 % sensitivity at 0.5 W | ionisation ×0.8/×1.2: Te +11.6/−9.3 meV, n_e −5.6/+4.6 %. D_a ×0.8/×1.2: Te −11.6/+9.7 meV, n_e +5.9/−4.8 %. P_abs ×0.8/×1.2: Te +0.08/−0.06 meV, n_e −20.0/+20.0 %. | reported. Cl. 10 **open**: one-at-a-time only, no joint propagation and no EEDF systematic |
| 7 | **Cl. 7**: budget closure | 0.5 W steady state: \|P − ΣQ − dU/dt\|/P = **6.9e-14**, dU/dt/P = 1.9e-12. Table in `budget_0p5W.md`: elastic 82.9 %, excitation (87) 16.2 %, ionisation 0.62 %, wall ions 0.19 %, wall electrons 0.07 % | **closed** |
| 8 | **Cl. 8**: wall energy uses the particle flux | Q_wall,e / latched `wall_electron_energy_flux` = 1.000000000000000. Both are built from the same `wall_loss_rates` evaluation (unit test pins it at 1e-15) | **closed** |
| 9 | Regression, energy off | FULL arm vs `mainline-FULL`: final state, rates, wall flux, V, t and every trace snapshot **bit-identical**. nu_wall re-measured **15.954491 s⁻¹** (this worktree's .so). Tests `test/rmgpy/solver/` + `inputTest.py`: **398 passed, 1 skipped** before → **428 passed, 1 skipped** after; **434 / 1** after the clamp. Re-checked after the clamp: FULL-off still bit-identical to `mainline-FULL` (`clamp/FULL-off/`) | **closed** |
| 10 | Red-first unit tests | `test/rmgpy/solver/plasmaEnergyBalanceTest.py`: 23 tests run red before implementation (`red/stdout.log`), including a hand value for each loss term and the energy-conservation row. The declared-energy and input-keyword tests were added with the PM correction and written alongside the code, not run red first | closed, with that exception |
| 11 | Both streams captured | every arm dir has `stdout.log` + `stderr.log`, as do the suites, analysis and remeasure | **closed** |

Clauses 11–12 (experimental observable) are out of scope: no observable is ratified.

## Entry 91: Te with the mixing proxy on and off (0.5 W)

| case | Te (eV) | n_e (m⁻³) |
|---|---|---|
| 91 declared +0.0752 eV (production) | 0.89888 | 9.554e15 |
| 91 removed | 0.89769 (−1.19 meV) | 9.616e15 (+0.6 %) |
| 91 at its enthalpy −11.548 eV (ruled out; counterfactual) | 0.89885 (−0.03 meV) | 1.018e16 (+6.5 %) |

The PM expected that crediting 91 with its enthalpy would move Te. **It does not.** Te is pinned
by particle balance, so the spurious heating shows up only in n_e, as +6.5 % at 0.5 W. That
error is smaller than the declaration decision might suggest because elastic loss carries 83 %
of the power at 5 torr. Removing 91 moves Te by −1.19 meV, as report7 found under prescribed Te
(−1.3 meV).

## Findings

1. **Te is set by particle balance, n_e by power, exactly as the contract predicts.** Te spans
   1.1 meV (0.89873–0.89984 eV) over four decades of P_abs, tracking Te*(n_e), and n_e is linear in P_abs. The ill-posed n_e of the
   prescribed-Te model (21× over 0.9 meV) becomes well-posed. The price is that n_e now carries
   the full uncertainty of P_abs (d ln n_e/d ln P = 1.00) and ~0.25 of the ionisation rate and D_a.
2. **At 5 torr, elastic loss dominates the energy cost:** 2489 eV per lost electron–ion pair,
   83 % of the power. The absolute n_e therefore depends most on K_el at 0.9 eV. The L&L fit is
   an analytic fit; an LXCat cross-section cross-check was not possible because the fetch was
   refused by the sandbox permission layer. **Open:** check K_el(0.9 eV) against a tabulated
   momentum-transfer cross section, and whether 0.9 eV lies inside the fit's stated range.
3. **RMG carries argon at 39.8775 g/mol** (`rmgpy.molecule.element`), not the standard 39.948
   (Ar-40 is 39.962). Through m_e/M this shifts elastic loss, and hence n_e, by 0.15 %. The hand
   check with the standard mass is off by exactly that (1.47e-3); with RMG's mass it closes to
   2.5e-6. Not changed here (out of scope). Worth a ticket.
4. **Extinction run: fixed by PM ruling (2).** At P_abs = 0 the state is labelled `extinct`
   within 0.2 ms. Before the fix the FULL run exited 1 at 10 ms: once Te ≈ Tg, Ars decays to ~0 and
   its accepted value is solver noise (−2e-32 mol even at atol 1e-30), which the M8-A
   `check_wall_support` refused as a negative population. Ruling: a **neutral** population in
   [−atol_j, 0), with atol_j that species' own absolute tolerance, is clamped to exactly 0 in the
   published accepted state; anything more negative, any non-finite value and any **charged**
   negative stay refused. The clamp acts on the state the model reads and publishes at the
   accepted step, not on DASPK's own history (which carries the same noise inside its tolerance).
   Six tests in `plasmaWallTest.py`; the two clamp tests ran red first (`clamp-red/`). The P0 arm
   (default atol 1e-16) now runs to steady state at the hand floor 6.6e4/1.5856 = 4.16e4 m⁻³.
   The pre-clamp run is kept in `arms/P0-preclamp/`.

   **report7 refusals do not change outcome.** The five f0d9c arms refused on an Ars negative
   (FULL-91OFF, SWEEP-0.80, SWEEP-0.80-91OFF, SWEEP-0.85, SWEEP-0.85-91OFF; Ars −1.2e-19 to
   −1.0e-18 mol, atol 1e-16) were rerun with this engine (`clamp/report7/`). Each still exits 1 at
   the **same accepted step** (same t to all digits, trace byte-identical). The Ars negative is now
   clamped, but that same state also carries a **negative electron population** (−4.6e-19 to
   −9.5e-18 mol), which is charged and stays refused. The neutral check simply ran first and
   named Ars. The real limit of those sub-threshold prescribed-Te arms is that the electrons are
   below atol, which the energy-mode charged-atol ×1e-12 addresses and prescribed-Te mode does
   not. They are unchanged here.
5. **This closure has no minimum sustaining power.** Every P_abs > 0 sustains a discharge at Te*
   with n_e ∝ P_abs, down to 0.001 W (1.9e13 m⁻³), until n_e would meet the source floor near
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
   terms together are 0.26 % of P_abs here, so this is immaterial at 5 torr. It would matter at
   low pressure.

## What is and isn't predicted

See `summary.md` in the runs directory (five lines, report-ready). In short: Te at 5 torr is a
model result, within the Maxwellian-EEDF bias i258 already flagged. n_e is predicted
**conditionally on P_abs** and is not a prediction of any measured discharge until the
absorbed power and discharge mode are ratified and an observable is chosen (clauses 11–12).

## Reproduce

    R=/home/alon/runs/i274-energy-balance-20260924-002658
    $R/scripts/run_arm.sh P0.5 --power 0.5          # one arm (see results.md for the rest)
    python $R/scripts/analyse.py                    # tables, figure, summary
    PYTHONPATH=$PWD pytest --no-cov test/rmgpy/solver/plasmaEnergyBalanceTest.py
