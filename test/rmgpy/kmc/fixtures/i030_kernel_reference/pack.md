# MET-KERNEL-REFERENCE: executable Rouse first-contact reference

Pack: `test/rmgpy/kmc/fixtures/i030_kernel_reference/pack.md`.
Read-only target: `rmgpy/kmc/met.py` in the enclosing repository checkout.
The skipped target is `test_met_kernel_reference` in the explicitly named `test/rmgpy/kmc/metTest.py`.
Source HEAD at original inspection: `d027e21d0bb48bb0e624835815acfdd8ee96c3b4`.
The fixture was copied to `/home/alon/Code/RMG-Py-kmc-i030-kernel-reference`
on branch `i030-met-kernel-reference`, based on `2798f5ce4e34e651943d8c7a78e30bfe7aa3a227`.
Its target module is byte-for-byte identical to the originally verified module.
This pack does not change either file. It accesses no prohibited data or validation outcomes.

## What this test means

The executable quantitative oracle is a pair of independent, equilibrium Gaussian Rouse chains, each with 16 beads, one reactive site on each chain, and irreversible first microscopic contact. Ends are bead 0; the central site is bead 7 (indices start at zero). Three separate ensembles measure end/end, end/mid and mid/mid. Exact Gaussian mode dynamics are projected onto the reactive-site separation and sampled using positive circulant embedding; this has the same finite-dimensional distribution as evolving every bead. The only trajectory discretization is contact detection. A boundary correction and separate timestep and box audits qualify the measured long-time oracle. An entangled extension supplies ideal-tube asymptotic exponents and normalized rate ratios; it does not invent site-specific numerical prefactors.

This is a build test of a deliberately weak long-time scaling approximation. It is not a demonstration that the present kernel resolves first contact accurately at all times or predicts nascent-chain conformations. The transient diagnostic remains visible even when the long-time assertion passes.

## Literal parameters and mapping to the actual API

The physical model, nine long-time/audit ensembles and acceptance policy were written before the first simulation or candidate comparison; their unchanged snapshot is `reference/parameters_initial.json`. The active `reference/parameters.json` adds one short-time precision allocation, declared after the base long-time rates but before any transient comparison. No physical parameter or tolerance changed. No parameter is fitted to an encounter result or experiment. The geometry comes only from constants in the named source, and the friction comes only from its selected transport arm.

| Parameter | Literal specification and source mapping |
|---|---|
| Temperature | 700 K; passed as `temperature` to `TransportArm.chain_diffusivity` and `diffusion_rate` |
| Rouse transport | `TRANSPORT_ARMS["A1_H_ROUSE"]`; H package, Rouse scaling, no entanglement cutoff in this arm |
| Polymer lengths | `units_i = units_j = 16`; 16 Gaussian beads are a coarse-grained realization of these 16 transport units, not an atomistic chain |
| Static size | `Rg² = PS_C_R2 * PS_M0 * 16 / 6`; source constants `PS_C_R2 = 0.434e-20`, `PS_M0 = 104.15` |
| RMS bond | `b² = 6*N*Rg²/(N²-1)`; this makes the exact discrete-chain `Rg` equal the kernel's `Rg` |
| Bond energy | `U = (3*k_B*T/(2*b²))*sum(|r[n+1]-r[n]|²)`; zero preferred bond vector, harmonic spring constant `3*k_B*T/b²` |
| Bead friction | `zeta = k_B*T/(N*D_CM)`; bead diffusivity `k_B*T/zeta = N*D_CM = arm.d0(T)` |
| Long-time COM diffusion | `D_CM = arm.chain_diffusivity(700,16)`; relative COM diffusion is `2*D_CM` |
| **Physical site contact** | `a = SIGMA_CONTACT = 6.766081442101e-10 m` (the dispatch's σ₀); first `|r_site_i-r_site_j| ≤ a` kills the pair |
| Candidate effective radius | `sigma_ij = max(SIGMA_CONTACT,2*Rg)`; this is the candidate's effective coil capture radius, not the microscopic reference sink |
| Chemistry and spin | Perfect absorption, no activation barrier, no cage, `spin_factor=1`; compare `diffusion_rate`, not activation-limited `bulk_rate` |
| Spatial initial condition | Reactive-site separation uniformly distributed in a periodic cube, conditioned outside the physical sphere `a`; equilibrium independent internal chain coordinates |
| Conformational initial condition | Equilibrium Gaussian chains. Time zero activates one reactive site on each existing chain. It does not create collapsed, stretched or scission-correlated chains |
| Box/concentration | One unlike pair per cube; base/refinement `L=10*Rg`, size audit `L=14*Rg`; each replica is independent, not an interacting many-chain melt; site density for each species is `1/L³`, molar concentration `1/(N_A*L³)` |
| Reduced units | Length `Rg`, time `t0=Rg²/D_CM`, rate `D_CM*Rg`; molar rate unit `N_A*D_CM*Rg` |
| Time grid | End time `8*t0`; bin edges `0, .02, .05, .1, .2, .4, .8, 1.6, 3.2, 4, 6, 8` times `t0` |
| Base ensemble | 4096 replicas/class, step `0.0001*t0`, seeds 103001, 103002, 103003 in class order end/end, end/mid, mid/mid |
| Step refinement | 2048 replicas/class, step `0.000025*t0`, seeds 103101, 103102, 103103; this run is the final oracle regardless of measured agreement |
| Box refinement | 4096 replicas/class, `L=14*Rg`, step `0.0001*t0`, seeds 103201, 103202, 103203 |
| Dedicated transient ensemble | 32768 replicas/class, `L=6*Rg`, step `0.000025*t0`, stop at `0.2*t0`, seeds 103301, 103302, 103303; independent pairs at number density `1/(216*Rg³)`; all four transient bins are reported |
| Random stream | NumPy `Generator(PCG64(seed))`; independent seeds across both classes and audit runs |
| Compute | Three worker processes at most, one FFT/BLAS thread each; `--workers 1` gives identical serial results; actual peak memory and elapsed times are recorded in the logs |

The actual physical values, concentrations, all survival values, all rate estimates and errors are preserved in `reference/results/results.json`; the tables below give the parameter values and primary observables. API units are m³ mol⁻¹ s⁻¹: `diffusion_rate` includes `N_A`. This reference's per-pair volume coefficient has units m³ s⁻¹. There is no factor of two for identical-reactant consumption: this is an unlike-labelled A–B pair coefficient, including for end/end. An initial serial run was interrupted after a completed base end/end ensemble to switch to three workers for runtime. Its partial log remains in `stdout.log`; the subsequent completed long-time ensembles plus the declared transient supplement supply the pack numbers. The Verifier repeats all twelve ensembles together. No model parameter or tolerance changed.

## Reference dynamics and first-passage numerics

For a chain with free ends, use the orthonormal modes

```
B_np = sqrt(2/N)*cos(pi*p*(n+1/2)/N), p=1,...,N-1
lambda_p = 4*sin²(pi*p/(2*N))
omega_p = (3*D_bead/b²)*lambda_p
Var(q_p,axis) = b²/(3*lambda_p)
```

Each `q_p` is an independent stationary Ornstein–Uhlenbeck process. For sites `n_i,n_j` on independent chains, the relative internal coordinate has per-axis covariance

```
C_ij(t) = sum_p [(B_n_i,p²+B_n_j,p²)*b²/(3*lambda_p)]*exp(-omega_p*|t|).
```

The separation path is `r(t)=r(0)+W_rel_CM(t)+q_rel(t)-q_rel(0)`. The three axes are independent. `W_rel_CM` has per-axis variance `4*D_CM*t`. The code samples the stationary Gaussian process with exactly this covariance on the observation grid using a circulant embedding of at least twice the trajectory length. It checks that the spectrum is nonnegative and that its inverse FFT recovers the specified covariance. An independently constructed bond Laplacian eigensystem checks the mode weights and exact discrete-chain `Rg`. Full paths are needed: neither MSD matching nor a memoryless site random walk is a first-passage reference.

The Rouse relaxation time is the slowest mode time `tau_R=1/omega_1`. The transient reporting window ends at `0.2*t0`, approximately `tau_R`. Long-time extraction uses the predeclared `[4,8]*t0` window; the two halves `[4,6]` and `[6,8]` provide a plateau audit. This definition gives a finite-run approximation to the asymptote, with a stated numerical audit, not an exact infinite-time measurement.

At arbitrarily short time the site-pair diffusion coefficient is `D_local=2*N*D_CM`, for all site classes. Discrete observations miss some between-step contacts. The stopping surface is therefore shifted outward, to

```
a_num = a + 0.5825971579390107*sqrt(2*D_local*dt).
```

This is the inward-domain boundary correction for the allowed region outside the absorbing sphere. It approximates continuum first contact; it is not a second physical contact radius. The physical initial exclusion uses `a`, not `a_num`, and no contact is counted at time zero. The sign follows the killed-diffusion boundary-shift prescription of [Gobet and Menozzi](https://arxiv.org/abs/0706.4042): discrete killing requires a smaller allowed domain. Applying this prescription to exact OU sampling is an approximation whose timestep audit is required. The quoted correction constant and the observed step difference are not fitted to the kernel.

Periodic minimum-image separation is used for contact only; the ideal Gaussian chain remains unwrapped. The long-time finite-cube hazard has an enhancement caused by periodic returns. Use the leading small-target periodic Green-function correction:

```
k_inf = k_box / (1 + 2.837297479*k_box/(4*pi*D_rel*L)), D_rel=2*D_CM.
dk_inf/dk_box = 1/(1 + 2.837297479*k_box/(4*pi*D_rel*L))².
```

This is a matched far-field correction, not an exact reduction of polymer first passage to an absorbing Brownian sphere. Its independent larger-box audit is essential. The dimensionless regular part of the cubic Laplace Green function can be calculated, without a polymer model, using an Ewald split with `alpha=2`:

```
sum_(n≠0) erfc(alpha*|n|)/|n|
+ sum_(m≠0) exp(-pi²*|m|²/alpha²)/(pi*|m|²)
- pi/alpha² - 2*alpha/sqrt(pi) = -2.837297479...
```

Remaining finite-volume, late-tail and contact-resolution errors are diagnosed, not absorbed into a calibrated capture radius.

## Observable and extraction from both sides

Save the first-contact time `T` of each replica, censoring at `8*t0`. Estimate `S(t)=P(T>t)` by the surviving fraction. The one-standard-error binomial error is `sqrt(S*(1-S)/M)`. For edges `u<v`, with survivor counts `n_u,n_v`, extract the finite-volume bin rate as

```
k_bin = L³*log(n_u/n_v)/(v-u)
SE(k_bin) = L³/(v-u)*sqrt((n_u-n_v)/(n_u*n_v)).
```

These formulas use conditional survival; nested survivor counts are correlated, and treating their errors as independent would overstate this error. The late-window error is propagated through the periodic correction with its derivative. In the dilute, infinite-volume limit the encounter-volume observable is `L³*(1-S(t))` and `k(t)` is its derivative. Finite-cell hazard and encounter-volume derivative differ after appreciable depletion; the transient table explicitly reports the former and does not call it an exact infinite-volume coefficient. Reaction is irreversible, so contacts after the first do not count. There is no competing partner and no replenishment.

For the present kMC transport law, call

```python
met_module.diffusion_rate(
    met_module.TRANSPORT_ARMS["A1_H_ROUSE"], 700.0, 16.0, 16.0,
    pair_class, spin_factor=1.0,
)
```

Divide the returned rate by `N_A` to obtain a per-pair volume coefficient. A single well-mixed pair at this concentration has `S_kmc(t)=exp[-k_D*t/(N_A*L³)]`. The bin hazard is constant, exactly `k_D/N_A`. No stochastic kMC run is needed to estimate an analytically known exponential rate. If using an SSA engine instead, use one A–B pair, propensity `k_D/(N_A*L³)`, one irreversible encounter event, the same survival estimator, and additional SSA error bars. Do not insert the microscopic radius into the candidate API: it has no radius or age argument.

## Acceptance policy declared before comparisons

All policies below are encoded in `reference/parameters.json` and apply independently to all three classes.

1. **Reference numerical qualification.** Step-refined versus base corrected long rates must agree within `4*sqrt(SE_step²+SE_base²)+0.10*k_base`. Larger versus base corrected long rates must agree within `4*sqrt(SE_box²+SE_base²)+0.15*k_base`. The base late-half rates must agree within `4*sqrt(SE_first²+SE_second²)+0.10*k_base`. Nonzero late events and survivors are required. Failure invalidates numerical acceptance; it does not justify loosening the kernel criterion. The allowances budget residual contact error, leading periodic-correction error and finite relaxation-window error, respectively. They are engineering accuracy targets checked against measurements, not rigorous bounds on unobserved systematic bias.
2. **Long-time kernel comparison.** Use the step-refined, corrected rate `k_ref` and its statistical standard error `s_ref`. Require `abs(log(k_candidate/k_ref)) <= log(2)+4*s_ref/k_ref`, after numerical qualification. With vanishing sampling error this allows a factor two in either direction. Four standard errors guard against chance failures across the three comparisons; it is an approximate Monte Carlo confidence allowance. A factor two is a declared scale-level model tolerance: compact exploration gives `k~D_CM*Rg`, while replacing a fluctuating absorbing site by a spherical coil sink leaves an uncontrolled order-one, site-dependent prefactor. It is **not** a derived physical confidence interval or proof that all those prefactors lie within two. A tighter prefactor-accuracy claim would need a different owner-approved specification; passing this loose test establishes only scale consistency for this stated Rouse model.
3. **Transient regime.** Report survival, every bin rate and uncertainty, and the ratio to the constant candidate. Use the dedicated short-time ensemble for this diagnostic; the large-box bins remain visible as well. A diagnostic is outside its band when `abs(k_candidate-k_bin)>0.20*k_bin+4*SE(k_bin)`. This creates no build failure, expected-failure decorator, or acceptance of a specific transient shape. It explicitly records mismatches. The 20% allowance recognizes finite concentration and boundary/discrete-chain approximations; the very short time bin includes local-bead diffusion and the boundary correction and is not a fit to a Rouse power law. No transient exponent is estimated from a handful of noisy bins. Increasing the short-time sample size was an allocation for precision, using a smaller volume and existing fine step, not a change to this policy.

The numerical qualification and factor-two comparison are separable assertions. Statistical uncertainty alone cannot validate the kernel's effective-radius prefactor. None of these tolerances were adjusted after viewing comparisons.

## Entangled extension: executable asymptotes and what they do not determine

Use `TRANSPORT_ARMS["A0_REF_H_CROSS_NEc"]`, `T=700 K`, `units_i=units_j=2048`, with literal `Ne=14800/104.15`, `b²=PS_C_R2*PS_M0`, and `D0=arm.d0(T)`. This arm documents `D_CM=D0*Ne/N²` above `Ne`.

The local bead friction is `zeta=k_B*T/D0=k_B*T*Ne/(N²*D_CM)`, equal to the H-package friction in the Rouse parameter table. The microscopic contact remains `SIGMA_CONTACT`, with perfect absorption and spin factor one. The continuum bond mapping is `b²=6*Rg²/N`; it differs from the finite-chain correction used in the Rouse simulation. These inputs establish the scaling conventions, not a determined first-contact amplitude.

```
ta=b²/D0, te=Ne²*ta, tau_R=N²*ta, tau_d=N³*ta/Ne,
r_e=b*sqrt(Ne), r_b=r_e*(N/Ne)^(1/4), R=b*sqrt(N).
x(t)~b*(t/ta)^(1/4)                 [ta << t << te]
x(t)~r_e*(t/te)^(1/8)               [te << t << tau_R]
x(t)~r_b*(t/tau_R)^(1/4)            [tau_R << t << tau_d]
x(t)~R*(t/tau_d)^(1/2)              [tau_d << t]
```

`R=sqrt(6)*Rg` in this convention; `R` must not silently be substituted for `Rg`. The times above are scaling time conventions, not exact Rouse eigenmode times or experimentally known entanglement times. They depend only on the named source inputs. No box, seed or integration step applies to this analytic extension. For high reactivity and dilute reactive sites, compact exploration gives the asymptotic rate exponents

| Pair | pre-tube | tube breathing | coherent reptation | relaxed long time |
|---|---:|---:|---:|---|
| end/end | −1/4 | −5/8 | −1/4 | constant ∝ `D_CM*Rg` |
| end/mid | −1/4 | −5/8 | −1/4 | constant ∝ `D_CM*Rg` |
| mid/mid | −1/4 | −5/8 | −1/4 | constant ∝ `D_CM*Rg` |

The end-functionalized-chain result is supported by [O'Shaughnessy and Vavylonis, equations (2), (23) and (26)](https://arxiv.org/pdf/cond-mat/9805331). The extension to mixed and central sites is an explicitly stated representative-monomer ideal-tube assumption about exponents, **not** a sourced calculation of their distinct prefactors. This also omits free-end tube effects. The source's approximate closure determines scaling, not precise encounter amplitudes.

This extension is executable in `tube_asymptotes()` in `reference/run.py`. It prints the source-mapped times and exponent windows for each class into JSON. Within a regime, ratios such as `k(t2)/k(t1)=(t2/t1)^q` eliminate the unknown amplitude. The candidate has ratio one at all times. For entangled lengths the long-time length ratio is `(N2/N1)^(-3/2)`; the candidate gives this ratio exactly. Its putative dimensionless amplitude is `16*pi` for equal chains and spin one. The literature does not make that amplitude an exact reference value. Therefore this branch validates only asymptotic length scaling and reports the missing time dependence; it has no numerical absolute-rate error bar and no independent entangled amplitude assertion.

For a literal analytic survival observable in a single regime, write `k(t)=k_star*(t/t_star)^q`, where each pair's `k_star` is unidentified. A dilute well-mixed pair then has `log[S(t2)/S(t1)] = -k_star*t_star*((t2/t_star)^(q+1)-(t1/t_star)^(q+1))/(V*(q+1))`. This symbolic expression and normalized rate ratios are the available asymptotic observables. They cannot supply absolute survival values until the independent amplitude is known.

An exact, class-resolved entangled first-passage oracle would require a specified tube/repton model with free ends, tube renewal and contour-length fluctuations, its own converged simulation and uncertainties. That is an adoption limitation, not evidence supplied by this Rouse simulation. The pack must not label the tube scaling table as a measured reptation simulation or claim an exact prefactor reference.

## Actual measurements

<!-- BEGIN REPRODUCED NUMBERS -->

All ± entries below are one Monte Carlo standard error; systematic/model errors are separate.

| Quantity | Value |
|---|---:|
| `N` | 16 |
| `temperature_K` | 700 |
| `Rg_m` | 1.09789009772e-09 |
| `bond_rms_m` | 6.73634613241e-10 |
| `sigma0_m` | 6.7660814421e-10 |
| `sigma_effective_m` | 2.19578019544e-09 |
| `D0_m2_s` | 2.75536362615e-08 |
| `D_CM_m2_s` | 1.72210226634e-09 |
| `bead_friction_kg_s` | 3.50753813699e-13 |
| `spring_constant_N_m` | 0.0638930748073 |
| `time_unit_s` | 6.99936751855e-10 |
| `tau_R_reduced` | 0.204091899843 |
| `tau_R_s` | 1.42851421456e-10 |
| `capture_reduced` | 0.616280396022 |
| `rate_unit_m3_s` | 1.89067902547e-18 |
| `rate_unit_m3_mol_s` | 1138593.52234 |
| `bond_squared_reduced` | 0.376470588235 |
| `entanglement_threshold_units` | 142.102736438 |
| `mid_bead_index_zero_based` | 7 |

| Run | Pair | seed | microscopic numerical radius / Rg | L / Rg | replicas | dt / t0 | late events | raw k_box / (D_CM Rg) | corrected k∞ / (D_CM Rg) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| base | end/end | 103001 | 0.66288817 | 10 | 4096 | 0.0001 | 786 | 85.03684 ± 3.0478072 | 43.386099 ± 0.79336762 |
| base | end/mid | 103002 | 0.66288817 | 10 | 4096 | 0.0001 | 746 | 72.568661 ± 2.6662653 | 39.889425 ± 0.80560227 |
| base | mid/mid | 103003 | 0.66288817 | 10 | 4096 | 0.0001 | 672 | 58.172124 ± 2.2491044 | 35.112838 ± 0.81942956 |
| step | end/end | 103101 | 0.63958428 | 10 | 2048 | 2.5e-05 | 430 | 95.269524 ± 4.6221544 | 45.901492 ± 1.0729761 |
| step | end/mid | 103102 | 0.63958428 | 10 | 2048 | 2.5e-05 | 370 | 72.77522 ± 3.7967762 | 39.951756 ± 1.1442468 |
| step | mid/mid | 103103 | 0.63958428 | 10 | 2048 | 2.5e-05 | 312 | 53.629374 ± 3.0419903 | 33.40488 ± 1.1802451 |
| box | end/end | 103201 | 0.66288817 | 14 | 4096 | 0.0001 | 375 | 76.041149 ± 3.9287587 | 47.13752 ± 1.5097028 |
| box | end/mid | 103202 | 0.66288817 | 14 | 4096 | 0.0001 | 332 | 65.10223 ± 3.57429 | 42.690882 ± 1.5369804 |
| box | mid/mid | 103203 | 0.66288817 | 14 | 4096 | 0.0001 | 269 | 50.936343 ± 3.1063571 | 36.106171 ± 1.5608413 |

Final long-time oracle uses the step-refined run, as declared before simulation.

| Pair | Reference k∞ (m³ mol⁻¹ s⁻¹) | Candidate k (m³ mol⁻¹ s⁻¹) | Candidate / reference | numerical audits | factor-two scale comparison |
|---|---:|---:|---:|---|---|
| end/end | 52263142 ± 1221683.7 | 57231953 | 1.095073 | PASS | PASS |
| end/mid | 45488811 ± 1302831.9 | 57231953 | 1.2581545 | PASS | PASS |
| mid/mid | 38034580 ± 1343819.5 | 57231953 | 1.5047347 | PASS | PASS |

Step-refined end/end; L/Rg=10, replicas=2048, seed=103101, dt/t0=2.5e-05. Time is in t0; bin k is a finite-box hazard coefficient.

| t1 | t2 | S(t2) ± SE | k_bin / (D_CM Rg) ± SE | k_bin (m³ mol⁻¹ s⁻¹) ± SE |
|---:|---:|---:|---:|---:|
| 0 | 0.02 | 0.98925781 ± 0.0022779078 | 540.01507 ± 115.13216 | 6.1485766e+08 ± 1.3108873e+08 |
| 0.02 | 0.05 | 0.98388672 ± 0.0027822719 | 181.47368 ± 54.716441 | 2.0662476e+08 ± 62299785 |
| 0.05 | 0.1 | 0.97753906 ± 0.0032742816 | 129.45029 ± 35.903113 | 1.4739126e+08 ± 40879052 |
| 0.1 | 0.2 | 0.96191406 ± 0.0042294655 | 161.13138 ± 28.484581 | 1.8346315e+08 ± 32432360 |
| 0.2 | 0.4 | 0.93945312 ± 0.005270095 | 118.13595 ± 17.418601 | 1.3450883e+08 ± 19832706 |
| 0.4 | 0.8 | 0.89550781 ± 0.0067594541 | 119.76745 ± 12.625805 | 1.3636644e+08 ± 14375659 |
| 0.8 | 1.6 | 0.82519531 ± 0.0083924727 | 102.21356 ± 8.5201696 | 1.1637969e+08 ± 9701009.9 |
| 1.6 | 3.2 | 0.71630859 ± 0.0099611205 | 88.443144 ± 5.9275338 | 1.0070079e+08 ± 6749051.6 |
| 3.2 | 4 | 0.66259766 ± 0.010448021 | 97.428898 ± 9.2918326 | 1.1093191e+08 ± 10579620 |
| 4 | 6 | 0.54882812 ± 0.010995734 | 94.191315 ± 6.1798066 | 1.0724562e+08 ± 7036287.8 |
| 6 | 8 | 0.45263672 ± 0.010998862 | 96.347732 ± 6.875117 | 1.097009e+08 ± 7827963.7 |

Step-refined end/mid; L/Rg=10, replicas=2048, seed=103102, dt/t0=2.5e-05. Time is in t0; bin k is a finite-box hazard coefficient.

| t1 | t2 | S(t2) ± SE | k_bin / (D_CM Rg) ± SE | k_bin (m³ mol⁻¹ s⁻¹) ± SE |
|---:|---:|---:|---:|---:|
| 0 | 0.02 | 0.98828125 ± 0.0023780224 | 589.39779 ± 120.31102 | 6.710845e+08 ± 1.3698534e+08 |
| 0.02 | 0.05 | 0.98291016 ± 0.0028639207 | 181.65349 ± 54.770656 | 2.0682949e+08 ± 62361515 |
| 0.05 | 0.1 | 0.97363281 ± 0.0035404994 | 189.6695 ± 43.513331 | 2.1595647e+08 ± 49543997 |
| 0.1 | 0.2 | 0.96142578 ± 0.0042554106 | 126.16872 ± 25.233911 | 1.4365489e+08 ± 28731168 |
| 0.2 | 0.4 | 0.94287109 ± 0.0051284856 | 97.438981 ± 15.80694 | 1.1094339e+08 ± 17997680 |
| 0.4 | 0.8 | 0.91308594 ± 0.0062249501 | 80.248931 ± 10.275264 | 91370914 ± 11699349 |
| 0.8 | 1.6 | 0.85986328 ± 0.0076705357 | 75.070752 ± 7.1915546 | 85475072 ± 8188257.5 |
| 1.6 | 3.2 | 0.75732422 ± 0.0094730355 | 79.363716 ± 5.4802969 | 90363013 ± 6239830.6 |
| 3.2 | 4 | 0.71533203 ± 0.0099714464 | 71.305802 ± 7.690145 | 81188324 ± 8755949.3 |
| 4 | 6 | 0.61865234 ± 0.010732945 | 72.601671 ± 5.1641099 | 82663792 ± 5879822.1 |
| 6 | 8 | 0.53466797 ± 0.011021954 | 72.948769 ± 5.5672261 | 83058996 ± 6338807.6 |

Step-refined mid/mid; L/Rg=10, replicas=2048, seed=103103, dt/t0=2.5e-05. Time is in t0; bin k is a finite-box hazard coefficient.

| t1 | t2 | S(t2) ± SE | k_bin / (D_CM Rg) ± SE | k_bin (m³ mol⁻¹ s⁻¹) ± SE |
|---:|---:|---:|---:|---:|
| 0 | 0.02 | 0.98876953 ± 0.0023285282 | 564.70033 ± 117.74878 | 6.4296414e+08 ± 1.34068e+08 |
| 0.02 | 0.05 | 0.98535156 ± 0.0026547662 | 115.42595 ± 43.626932 | 1.3142324e+08 ± 49673342 |
| 0.05 | 0.1 | 0.97509766 ± 0.0034433343 | 209.21735 ± 45.655177 | 2.3821352e+08 ± 51982689 |
| 0.1 | 0.2 | 0.96533203 ± 0.0040423841 | 100.6551 ± 22.50726 | 1.1460525e+08 ± 25626621 |
| 0.2 | 0.4 | 0.953125 ± 0.0046706852 | 63.630281 ± 12.726142 | 72449026 ± 14489903 |
| 0.4 | 0.8 | 0.93164062 ± 0.0055764559 | 56.997287 ± 8.5928505 | 64896742 ± 9783763.9 |
| 0.8 | 1.6 | 0.89160156 ± 0.0068696078 | 54.909739 ± 6.0642538 | 62519873 ± 6904720.1 |
| 1.6 | 3.2 | 0.82226562 ± 0.0084474729 | 50.597416 ± 4.2472005 | 57609891 ± 4835834.9 |
| 3.2 | 4 | 0.7890625 ± 0.0090150393 | 51.522445 ± 6.2484565 | 58663122 ± 7114452 |
| 4 | 6 | 0.70117188 ± 0.010114816 | 59.046245 ± 4.403605 | 67229672 ± 5013916.2 |
| 6 | 8 | 0.63671875 ± 0.010627481 | 48.212504 ± 4.1979857 | 54894444 ± 4779799.3 |

Dedicated transient end/end; L/Rg=6, replicas=32768, seed=103301, dt/t0=2.5e-05. Time is in t0; bin k is a finite-box hazard coefficient.

| t1 | t2 | S(t2) ± SE | k_bin / (D_CM Rg) ± SE | k_bin (m³ mol⁻¹ s⁻¹) ± SE |
|---:|---:|---:|---:|---:|
| 0 | 0.02 | 0.96325684 ± 0.0010392843 | 404.30015 ± 11.652417 | 4.6033353e+08 ± 13267366 |
| 0.02 | 0.05 | 0.93295288 ± 0.0013816402 | 230.15076 ± 7.3039246 | 2.6204817e+08 ± 8316201.2 |
| 0.05 | 0.1 | 0.8888855 ± 0.0017361343 | 209.02868 ± 5.5012913 | 2.379987e+08 ± 6263734.7 |
| 0.1 | 0.2 | 0.81985474 ± 0.002123024 | 174.61711 ± 3.6724742 | 1.9881791e+08 ± 4181455.3 |

Dedicated transient end/mid; L/Rg=6, replicas=32768, seed=103302, dt/t0=2.5e-05. Time is in t0; bin k is a finite-box hazard coefficient.

| t1 | t2 | S(t2) ± SE | k_bin / (D_CM Rg) ± SE | k_bin (m³ mol⁻¹ s⁻¹) ± SE |
|---:|---:|---:|---:|---:|
| 0 | 0.02 | 0.96447754 ± 0.0010225219 | 390.62233 ± 11.449967 | 4.4476006e+08 ± 13036859 |
| 0.02 | 0.05 | 0.93676758 ± 0.0013445002 | 209.88966 ± 6.9656797 | 2.3897901e+08 ± 7931077.8 |
| 0.05 | 0.1 | 0.89883423 ± 0.0016658337 | 178.57403 ± 5.0653989 | 2.0332323e+08 ± 5767430.3 |
| 0.1 | 0.2 | 0.84124756 ± 0.0020188179 | 143.01931 ± 3.2929721 | 1.6284086e+08 ± 3749356.7 |

Dedicated transient mid/mid; L/Rg=6, replicas=32768, seed=103303, dt/t0=2.5e-05. Time is in t0; bin k is a finite-box hazard coefficient.

| t1 | t2 | S(t2) ± SE | k_bin / (D_CM Rg) ± SE | k_bin (m³ mol⁻¹ s⁻¹) ± SE |
|---:|---:|---:|---:|---:|
| 0 | 0.02 | 0.96569824 ± 0.0010054349 | 376.96182 ± 11.244399 | 4.2920629e+08 ± 12802800 |
| 0.02 | 0.05 | 0.94158936 ± 0.0012955429 | 182.03152 ± 6.4765633 | 2.0725991e+08 ± 7374173.1 |
| 0.05 | 0.1 | 0.90960693 ± 0.0015840521 | 149.28488 ± 4.6116549 | 1.699748e+08 ± 5250800.4 |
| 0.1 | 0.2 | 0.86398315 ± 0.001893756 | 111.15207 ± 2.8750466 | 1.2655703e+08 ± 3273509.5 |

| Pair | Transient t1 | t2 | Candidate / finite-box reference | 20% + 4 SE diagnostic |
|---|---:|---:|---:|---|
| end/end | 0 | 0.02 | 0.12432714 | outside band (report-only) |
| end/end | 0.02 | 0.05 | 0.21840242 | outside band (report-only) |
| end/end | 0.05 | 0.1 | 0.2404717 | outside band (report-only) |
| end/end | 0.1 | 0.2 | 0.28786115 | outside band (report-only) |
| end/mid | 0 | 0.02 | 0.12868051 | outside band (report-only) |
| end/mid | 0.02 | 0.05 | 0.23948527 | outside band (report-only) |
| end/mid | 0.05 | 0.1 | 0.2814826 | outside band (report-only) |
| end/mid | 0.1 | 0.2 | 0.35145941 | outside band (report-only) |
| mid/mid | 0 | 0.02 | 0.1333437 | outside band (report-only) |
| mid/mid | 0.02 | 0.05 | 0.27613615 | outside band (report-only) |
| mid/mid | 0.05 | 0.1 | 0.33670845 | outside band (report-only) |
| mid/mid | 0.1 | 0.2 | 0.45222262 | outside band (report-only) |

| Entangled theory parameter | Value |
|---|---:|
| `units` | 2048 |
| `Ne_units` | 142.102736438 |
| `D_CM_m2_s` | 9.33515336888e-13 |
| `Rg_m` | 1.24212085295e-08 |
| `bond_rms_m` | 6.72317633266e-10 |
| `ta_s` | 1.64047676216e-11 |
| `te_s` | 3.31264551809e-07 |
| `tau_R_s` | 6.88065824544e-05 |
| `tau_d_s` | 0.000991647904882 |

Entangled entries are scaling predictions with unidentified prefactors, not measured rates or error bars.

<!-- END REPRODUCED NUMBERS -->

## Replacement pytest body (text only)

This replaces the skip and empty body after adoption. It keeps the real transport API and arms; it does not copy the candidate formula into a mock. The external reference is intentionally slow and is regenerated before assertions. The fixture path defaults to the sibling `fixtures/i030_kernel_reference` directory when adopted into `metTest.py`; the wrapper provisions `MET_KERNEL_REFERENCE_PACK` for its isolated temporary test. The scaffold imports `os` along with the other standard-library imports used below. A reference audit failure must fail the test. The entangled assertions test ratios and never claim amplitude agreement.

```python
def test_met_kernel_reference(tmp_path):
    pack_root = Path(os.environ.get(
        "MET_KERNEL_REFERENCE_PACK",
        Path(__file__).resolve().parent / "fixtures/i030_kernel_reference",
    ))
    reference_output = tmp_path / "reference"
    # Build a fresh oracle for these inputs; do not assert the old candidate's fingerprint.
    subprocess.run(
        [sys.executable, str(pack_root / "reference/run.py"), "--output", str(reference_output)],
        cwd=pack_root, check=True,
    )
    cfg = json.loads((pack_root / "reference/parameters.json").read_text())
    data = json.loads((reference_output / "results.json").read_text())
    assert data["parameters_sha256"] == hashlib.sha256(
        (pack_root / "reference/parameters.json").read_bytes()
    ).hexdigest()
    assert data["met_source_sha256"] == hashlib.sha256(
        Path(met_module.__file__).read_bytes()
    ).hexdigest()
    assert data["independent_checks"]["FFT_covariance_check"] == "PASS"
    arm = met_module.TRANSPORT_ARMS[cfg["rouse_arm"]]
    mapped = data["mapped"]
    comparisons = {row["pair_class"]: row for row in data["comparison"]}
    assert set(comparisons) == {"end/end", "end/mid", "mid/mid"}
    for pair_class, row in comparisons.items():
        assert all(check["pass"] for check in row["numerical_checks"].values()), row
        k = met_module.diffusion_rate(
            arm, cfg["temperature_K"], mapped["N"], mapped["N"],
            pair_class, spin_factor=1.0,
        ) / mapped["rate_unit_m3_mol_s"]
        kref = row["reference_reduced"]
        se = row["reference_se_reduced"]
        assert abs(math.log(k/kref)) <= math.log(cfg["long_model_factor"]) + 4*se/kref
        # Report every transient diagnostic, including disagreement; no xfail or assertion.
        print(pair_class, "transient/report-only", row["transient_reporting_only"])
    tube_arm = met_module.TRANSPORT_ARMS[cfg["entangled_arm"]]
    n1, n2 = 2048.0, 4096.0
    assert n1 > tube_arm.n_e
    for pair_class in comparisons:
        k1 = met_module.diffusion_rate(tube_arm, 700.0, n1, n1, pair_class, spin_factor=1.0)
        k2 = met_module.diffusion_rate(tube_arm, 700.0, n2, n2, pair_class, spin_factor=1.0)
        assert math.isclose(k2/k1, (n2/n1)**(-1.5), rel_tol=1e-12)
        theory = next(row for row in data["tube"]["rows"] if row["pair_class"] == pair_class)
        for regime in theory["transient_regimes"]:
            # Normalized observable: two times separated by a factor four in one regime.
            # Amplitudes cancel. This is a scaling prediction, with no Monte Carlo error.
            ratio_reference = 4.0**regime["k_exponent"]
            print(pair_class, regime["name"], "k(4t)/k(t): theory", ratio_reference,
                  "candidate", 1.0, "report-only")
```

## Commands, reproducibility and verification evidence

Python is `/home/alon/anaconda3/envs/rmg_env/bin/python`. The observed environment has NumPy 1.26.4 and SciPy 1.13.1. No additional dependencies, RMG compiled extensions, mechanism generation or database loading are required. `reference/` contains one runnable script, `run.py`, plus its frozen JSON configuration and output artifacts. There is no hidden build step. To run:

```bash
cd /home/alon/Code/RMG-Py-kmc-i030-kernel-reference/test/rmgpy/kmc/fixtures/i030_kernel_reference
/home/alon/anaconda3/envs/rmg_env/bin/python -u reference/run.py > >(tee -a /home/alon/runs/i030-met-kernel-reference/fixture-reference-stdout.log) 2> >(tee -a /home/alon/runs/i030-met-kernel-reference/fixture-reference-stderr.log >&2)
/home/alon/anaconda3/envs/rmg_env/bin/python -u reference/run.py --verify-again > >(tee -a /home/alon/runs/i030-met-kernel-reference/fixture-reference-verifier-stdout.log) 2> >(tee -a /home/alon/runs/i030-met-kernel-reference/fixture-reference-verifier-stderr.log >&2)
```

The first measured set was assembled from nine completed long-time/audit ensembles, then the predeclared precision extension using `reference/run.py --supplement-transient` (with the same two-stream log redirection). That option checks the original policy and model snapshot, adds only the three short-time ensembles and recomputes the diagnostic tables. A default run or `--verify-again` always generates all twelve; neither depends on partial results or this assembly shortcut.

The fixture includes only the small frozen baseline `reference/results/results.json` and `results.md` (65,814 bytes total), the script, both parameter files, this pack and its verifier. Original logs and redundant reproduced outputs remain in the run directory. Both scripts resolve `rmgpy/kmc/met.py` from their enclosing checkout. `reference/run.py` defaults to writing fresh results under `/home/alon/runs/i030-met-kernel-reference/fixture-verification/reference`; `--output` can select another output directory. `--verify-again` compares every trajectory against the fixture's frozen baseline and writes a separate `reproduced/` directory under its output directory. It requires exact JSON agreement (including first-contact-index hashes), and requires the freshly rendered numeric table to occur verbatim in this pack. Wall times and memory measurements live in logs, not in the deterministic comparison. Different NumPy/SciPy builds may change the random-to-FFT path at contact boundaries; an exact-reproduction failure then needs investigation rather than silently accepting changed numbers.

`verify_pack.py --assemble` embeds the generated numeric table into this document. `verify_pack.py` extracts the quoted pytest body, adds only imports of the real `met.py`, and executes it under an isolated pytest configuration in a temporary directory. It then requires exact agreement of the new results with the saved report and requires the new numeric table to occur verbatim in this pack. The pytest body itself compares the current candidate against a fresh reference, so an adopted build test need not reproduce an old candidate's fingerprint or comparison results. The wrapper's additional equality check verifies this frozen research report. Both checks execute without editing the worktree. The worker's final verification command is:

```bash
/home/alon/anaconda3/envs/rmg_env/bin/python -u verify_pack.py --output /home/alon/runs/i030-met-kernel-reference/fixture-verification > >(tee -a /home/alon/runs/i030-met-kernel-reference/fixture-verifier-stdout.log) 2> >(tee -a /home/alon/runs/i030-met-kernel-reference/fixture-verifier-stderr.log >&2)
```

The wrapper defaults to the same external output directory. Its pytest temporary directory,
fresh simulation results and retained `reproduced/results.json` and `results.md` all live
there; a verification run writes no files into the fixture or repository.

For pytest execution of the quoted body, use an isolated configuration (`-c /dev/null`, `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1`) from the pack directory. The shared repository's pytest defaults try to erase its coverage file; that operation was refused by the read-only filesystem during an initial `pytest --version` check. No repository pytest suite or shared coverage artifact is changed by the verification reported here.

The original worker executed the quoted pytest body through `verify_pack.py`, with both streams persisted in `/home/alon/runs/i030-met-kernel-reference/verifier-stdout.log` and `verifier-stderr.log`. The complete run exited 0. Actual terminal output before fixture relocation:

```text
1 passed in 1448.50s (0:24:08)
VERIFIER PASS: all 12 first-passage runs reproduced exactly; every generated pack number matches.
PACK VERIFIER PASS: reproduced table and proposed pytest assertions passed.
```

The three step audits, three box audits and three plateau audits passed. All three long-time scale assertions passed. All twelve dedicated transient bins were outside the predeclared diagnostic band; those differences were reported and did not fail the build. The candidate/reference transient ratios range from 0.12432714 to 0.45222262, as tabulated above. This confirms the stated limitation for the declared initial condition and finite-volume observable; it does not validate a nascent-chain conformation model.

The simulation used three workers with one FFT/BLAS thread each. Actual recorded memory was `PARENT_PEAK_RSS_MiB=56.0` and `CHILD_MAX_RSS_MiB=115.1`. Both simulation stderr logs are empty. The reference result and reproduced result JSON objects match exactly, including every contact-index hash. The numerical table remains an exact copy of the freshly reproduced table.

The reproduced target-source SHA-256 is `1a8176c2e49c3fcfb5909c73bf01ad8d7230f2e5ad724e7461984a7b8388a7c1`. The frozen active parameter-file SHA-256 is `0e99129828cd966e7acb3db73f2c202466941caa618ef9eae7874a2a0e1262fe`. The source-mapped values and all uncertainties are in the saved JSON. The proposed body was executed with the actual source module, not a mock. No repository test suite or product edit was performed.

The relocated fixture verifier was run from the fixture directory using the command above. Both streams were saved outside the repository in `fixture-verifier-stdout.log` and `fixture-verifier-stderr.log` under `/home/alon/runs/i030-met-kernel-reference/`. The complete run exited 0; stderr was empty. Actual output:

```text
1 passed in 1800.87s (0:30:00)
VERIFIER PASS: all 12 first-passage runs reproduced exactly; every generated pack number matches.
PACK VERIFIER PASS: reproduced table and proposed pytest assertions passed.
```

All twelve results, including contact-index hashes, were identical to the frozen baseline. Fresh retained outputs are in `/home/alon/runs/i030-met-kernel-reference/fixture-verification/reproduced/`. The rerun recorded `PARENT_PEAK_RSS_MiB=56.1`, `CHILD_MAX_RSS_MiB=114.4`, three workers and one FFT/BLAS thread per worker.

## Limitations and contract issues supported by evidence

- **Excluded volume and hydrodynamics.** This is a screened, phantom, free-draining Gaussian-chain model for an unentangled melt. It has no incompressibility field, hard cores, hydrodynamic interaction, chemistry-dependent bead friction or polydispersity. Bonds can cross. Its validity does not extend to unscreened dilute good-solvent chains. The bead-spring model is the standard [Rouse model, J. Chem. Phys. 21, 1272–1280 (1953)](https://doi.org/10.1063/1.1699180), applied here to screened melt chains.
- **Finite chain length.** Sixteen beads and the finite `a/Rg` give limited scale separation. The transient window includes crossover from local bead diffusion to Rouse motion. It cannot establish a universal power-law exponent. Exact finite-chain `Rg` mapping prevents an additional continuum-size error; it does not make a short coarse-grained chain atomistic.
- **Finite time and volume.** The primary long-time observable is a corrected late-window hazard, not a direct infinite-domain asymptote. Periodic absorption, depletion and finite-window effects are approximated and audited with the frozen budgets. These audits can be statistically weak when contact counts are small; four-error overlap is not evidence of zero systematic bias. Full raw counts and errors are available for adversarial review.
- **Chain birth.** The input lacks a birth conformation or initial separation correlation. Equilibrium activation is a literal chosen initial condition. Scission daughters, a stretched new chain, or a geminate pair could have different transient survival. This test cannot specify those from `sigma0`, `Rg` and `D_CM` alone.
- **Entanglement.** The cutoff `Ne` is the transport arm's documented input, not an empirically inferred onset. The Rouse simulation is below it. The ideal-tube branch supplies exact asymptotic exponents within its scaling assumptions and unknown class amplitudes; it is not an independent entangled first-contact simulation. Constraint release and contour-length fluctuations are omitted.
- **The contract overstates what is established at long time.** In `met.py:429–455`, the implementation chooses `max(sigma0,2*Rg)` and applies the same formula to all pair classes. Compact-exploration theory supplies an order-one coil-scale radius, not proof that the prefactor is exactly `2*Rg` or class independent. A statistical tolerance cannot convert that uncontrolled approximation into an exact theorem. The factor-two tolerance is an explicit, weak engineering definition for this build test.
- **The kernel cannot express the proposed transient dependence.** Its `diffusion_rate` signature has arm, temperature, two lengths, class and spin only, and its docstring explicitly says no class-specific multiplier is defined. It has no age/history input. This pack measures and reports that limitation; no product change is implied.
- **Microscopic versus effective capture must be specified.** Using `sigma_ij=2*Rg` as a site-contact distance inside a flexible-chain simulation would include coil-scale internal exploration as well as an already coarse-grained sink. This pack chooses the literal microscopic `sigma0` and compares its emergent long-time capture to the candidate's effective sphere. The dispatch did not settle that distinction. Both radii and the choice are explicit here; the choice was made before results.
- **Exact entangled amplitude data are absent.** The allowed scaling theory cannot by itself supply an absolute entangled rate with a statistical error for all three classes. A future adoption must either accept the explicitly limited asymptotic branch or commission that additional simulation. The current outputs do not certify an entangled amplitude test.

Adversarial review and project-owner adoption approval remain outstanding. No gate was needed to produce this pack, and no report-back pane was named.
