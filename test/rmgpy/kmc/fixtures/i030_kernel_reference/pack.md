# MET-KERNEL-REFERENCE: independent Rouse reference, rework round 1

This pack is a research fixture in `test/rmgpy/kmc/fixtures/i030_kernel_reference/`.
It changes no product code or repository tests. The target is the enclosing
checkout's `rmgpy/kmc/met.py`; the work branch is `i030-met-kernel-reference`,
starting at `52aa9addb488e08cf309560e5dbd3413d3a9a818`. Scratch, logs and fresh
verifier outputs are under `/home/alon/runs/i030-met-kernel-reference/rework/`.
No prohibited dataset, catalog, or pyrolysis literature is used.

## Meaning and adoption status

The independent observable is irreversible first microscopic contact of two
stationary, free-draining Gaussian Rouse chains in a periodic cube. The three
reactive-site classes are end/end, end/mid and mid/mid. Chains start at equilibrium,
with uniform initial site separation conditioned outside the physical contact
sphere. This is an unlike-labelled pair coefficient, including for end/end.
It measures neither a nascent-chain conformation nor a validation outcome.

Reproduction, numerical qualification, and candidate acceptance are distinct.
A reproduced report can contain a failed qualification or a rejected candidate.
The proposed scientific test below asserts all three; its real pytest outcome
is retained by the verifier. **Reproduction success is not adoption success.**
The internal `adoption_pass` flag records that proposed test's numerical
criteria. Owner approval and the second adversarial review are still required
before adoption.
The current target has one identical prefactor for all three classes. Consequently,
a swap of its end/end and mid/mid labels is observationally identical to the
unmodified target. No sound test can kill that no-op after passing the baseline.
The pack reports that separate label-only diagnostic as an expected survivor.
The review's required perturbation instead reverses the **reference-derived**
site prefactors in the candidate: multiply end/end by `C_mid/C_end` and mid/mid
by `C_end/C_mid`. This substantive in-memory source mutation is included in the
required set, with a passing unmodified literal-rule control. An independent
oracle-backed ordering witness is also reported, separately from rate-budget
discrimination and its numerical qualification.

## Frozen physical inputs and independence

`reference/parameters.json` is the sole physical-input source. `reference/run.py`
never imports, loads, evaluates, hashes or reads `met.py`. Target loading and
comparison live in `checks.py`, outside `reference/`. Parameter bytes are read
once, SHA-256 hashed, and parsed before any simulation. The fingerprint travels
with that in-memory snapshot. The oracle program hash and Python/NumPy/SciPy,
RNG, platform and machine versions are recorded; the verifier also fingerprints
all checking programs and records pytest's version. A program edit during a
run is detected. Exact reproduction is conditional on the recorded software
versions and numerical platform; it is not a portable bitwise guarantee.

The physical literals are benchmark conventions frozen from the named target
snapshot, **not newly established experimental polymer parameters**. Their
provenance and the independently evaluated H/WLF expression are cited in the
parameter JSON. Changing the current target cannot change these literals. The
transport comparison separately requires `D0`, every `D_CM`, the microscopic
radius, the static `Rg` coefficient, every length's kernel `Rg`, and SI constants
to match their literals to relative tolerance `2e-12`. This addresses numerical
independence, not epistemic independence of the original material convention.

The temperature is 700 K, `D0 = 2.755363626147875e-8 m2/s`, and microscopic
contact radius `sigma0 = 6.766081442101e-10 m`. The frozen size coefficient is
`Rg^2/N = 0.434e-20*104.15/6 m2`. `kB` and Avogadro's constant are SI defining
constants ([BIPM SI Brochure](https://www.bipm.org/en/publications/si-brochure)).
For each finite chain the JSON contains literal `D_CM`, `Rg`, and RMS bond
lengths satisfying

```
D_CM(N) = D0/N
Rg^2(N) = b_N^2*(N^2-1)/(6*N)
zeta = kB*T/D0
U = 3*kB*T/(2*b_N^2) * sum(|r[n+1]-r[n]|^2)
```

The finite-chain bond mapping changes slightly with N to make the exact discrete
size match the frozen kernel's static-size convention. It is not an atomistic
chain model or a fixed-bond asymptotic sequence. The underlying dynamics are
those of [Rouse, J. Chem. Phys. 21, 1272 (1953)](https://doi.org/10.1063/1.1699180).
No independent first-contact amplitude follows just from
[coil-scale reaction theory](https://arxiv.org/abs/cond-mat/9805331).

Coverage includes equal N=4,8,16, both 4/16 and 16/4 orientations, and N=1/1.
The swapped unequal orientations matter for end/mid. Together they distinguish
`min(i,j)` from `max`, either argument or their mean. The N=1 limit is a point
Brownian bead with physical Rg=0; its nonzero **kernel** static-size convention
is recorded separately. `sigma0 > 2*kernel_Rg(1)`, so the candidate's floor is
active. All three site labels coincide for a single bead; only one physical
ensemble per setting is generated and explicitly reused for those three labels.
N=4,8,16 are below the frozen transport cutoff 14800/104.15; only the Rouse arm
is compared. No entangled first-contact amplitude is claimed.

## Dynamics, sampler, and contact detector

For each chain, orthonormal free-end modes have

```
B_np = sqrt(2/N)*cos(pi*p*(n+1/2)/N)
lambda_p = 4*sin^2(pi*p/(2*N))
omega_p = 3*D0*lambda_p/b_N^2
Var(q_p,axis) = b_N^2/(3*lambda_p)
q_p(t+dt) = exp(-omega_p*dt)*q_p(t)
            + sqrt(Var(q_p)*[1-exp(-2*omega_p*dt)])*Z
```

Both chains' independent modes are evolved explicitly; the relative COM has
diffusivity `D_i+D_j`. End is bead 0 and mid is bead `N//2-1`. Their site
separation is relative COM plus the two weighted internal coordinates. COM
and modes remain unwrapped; minimum-image separation is used for contact only.
Exact conditional OU transitions avoid spring-integration error. Full state
is preserved: matching an MSD alone is not a first-passage oracle.

Observation intervals adapt to the current distance from the sink. With `g`
the distance outside the physical sink, the interval is bounded by
`[g/(safety*sqrt(2*D_local))]^2`, a maximum interval, and an exact conditional
mean-displacement guard; `D_local=2*D0`. The guard reduces the proposed interval
when its mean site displacement exceeds `g/safety`. It is then clipped below
by the declared minimum contact step. A pair is absorbed at the first observed
site separation below

```
a_numerical = sigma0 + 0.5825971579390107*sqrt(2*D_local*dt_min).
```

The physical initial exclusion uses `sigma0`, not this numerical surface. The
shift's sign follows the reduced allowed-domain prescription of
[Gobet–Menozzi](https://arxiv.org/abs/0706.4042). Their Euler theorem does not prove
an error bound for this exact-OU/adaptive implementation. The independent contact
checks and the separately refined minimum and far-step intervals therefore
matter. The **adaptive** audit halves the final maximum interval and increases
safety from 8 to 10 with the same minimum contact step. The step and contact
runs reduce the minimum interval by factors of four. These changes are not
fits to candidate rates.

The actual `advance_modes()` used by the trajectories is used to sample measured
stationary variances and covariances at several lags in all nontrivial cases and classes, with independent
replicas and Monte Carlo SE. An independently assembled bead bond Laplacian
checks the mode covariance, local diffusivity and exact finite-chain Rg.
A separate Euler implementation evolves every bead's force and noise at small
N=4, starting from independent Gaussian bonds. Its first-contact probability
is compared both with exact mode evolution on a fixed observation grid and
with the **actual adaptive `simulate()` trajectory sampler**,
using independent seeds, 5 combined SE and a separately stated Euler weak-error
scale. The detector's continuum behavior is additionally checked in the N=1
limit against the exact absorbing Brownian-sphere hitting probability integrated
over initial positions. Its omitted cube-exterior contribution is bounded by
the spherical `r>=L/2` tail. Real mutations doubling the trajectory sampler's
OU noise or projected internal coordinates must fail measured variance or
covariance after their unmodified controls pass.
The bead cross-check checks small-N fixed-grid and adaptive contact distributions; it
does not itself prove all adaptive large-N continuum errors are bounded.

## Rate extraction and empirical error budget

For late window `[u,v]`, surviving counts `n_u,n_v` give

```
k_box = L^3/(v-u)*log(n_u/n_v)
SE(k_box) = L^3/(v-u)*sqrt((n_u-n_v)/(n_u*n_v))
k_inf = k_box/(1 + 2.837297479*k_box/(4*pi*D_relative*L))
SE(k_inf) = SE(k_box)/(1 + 2.837297479*k_box/(4*pi*D_relative*L))^2
```

The conditional survival formula accounts for nested counts. The periodic
Green-function correction is a leading far-field approximation, not an exact
polymer reduction. Its cubic Ewald constant is recomputed independently. Three
boxes (20,28,36 times the case's length unit) and later windows scaled by relative
COM diffusion test its remaining bias. Every late window begins after at least
four relative-COM spatial mixing times; the corresponding values are tabulated.
The actual correction fractions are shown; the final largest-box correction is
capped at 25% during qualification,
against the rejected report's 38–52%. Statistical overlap does not establish
that this approximation is exact.

The **largest-box, contact-refined** run (`box2`) is the primary oracle. Its own
two late halves provide the plateau audit. All three boxes' halves are reported.
With `k,s` the
corrected rate and one SE, define the upper observed log contrast between two
runs as

```
U(A,B) = abs(log(k_A/k_B)) + 4*hypot(s_A/k_A, s_B/k_B)
```

The additive log error budget is

```
statistics = 4*s_box2/k_box2
step       = U(base,step)
contact    = U(step,contact)
adaptive   = U(contact,adaptive)
box        = max(U(contact,box),U(box,box2))
plateau    = U(box2_half1,box2_half2)
B_total    = statistics+step+contact+adaptive+box+plateau
```

A candidate falls within the numerical budget iff
`abs(log(k_candidate/k_reference)) <= B_total`. **There is no factor-two model
allowance.** The measured budgets include sampling uncertainty of the audits;
this is conservative and double-counts some statistical uncertainty. They are
empirical envelopes at sampled settings, not rigorous bounds on the unobserved
continuum/infinite-volume/infinite-time remainder or simultaneous confidence
intervals over every reported row. Passing the error-budget comparison alone
cannot justify a stronger mathematical statement.

Qualification also requires nonzero events and survivors, <=25% final correction,
>=4 mixing times, and individual log ceilings in the parameter JSON: stat .16,
step/contact/adaptive .12 each, box .16, plateau .20, total .55. The total cap
is a factor `exp(.55)=1.73325`, strictly below two: any qualified oracle-backed
rate control multiplied by 0.5,2 or4 must fail. A row exceeding any ceiling is
**unqualified**, not accepted with a looser tolerance. These ceilings are accuracy
targets; the comparison tolerance is the measured budget, not the ceiling.
The preliminary allocation was stopped before candidate comparison to increase
counts and batch size; its exact source, inputs and partial logs remain outside
the fixture. That preliminary allocation supplies no final measurements.
The final allocation uses three independent seed stages for each Rouse case
and eight for the point-bead floor, pooling integer survivors before rate
estimation. The three-stage floor exceeded the fixed step and total ceilings;
five additional floor stages reduce sampling uncertainty at little cost. The
ceilings and all physical/per-run settings stay fixed. `extend_floor.py` is an
authoring driver: it checks the completed three-stage input/program hashes,
dependency versions, unchanged physical inputs, and identical ASTs for every
generator and numerical-check function before adding the five independent
floor cohorts. The original whole program and parsed JSON bytes are archived
in `reference/provenance/`. Non-floor reference means, errors and budgets must
remain exactly equal to the completed three-stage report.

The normal reference CLI and the full Verifier generate all 318 physical
ensembles from scratch; neither reads authoring checkpoints. Each stage retains
its actual input, reference-program and driver hashes in `generation_provenance`.
The root hashes identify the final analysis allocation and reference program.
Authoring and fresh verification have genuinely different generation histories;
their histories are separately checked, retained and hashed. Every scientific
JSON value, event count, trajectory hash and table must reproduce exactly.
Only actual generation-history metadata is exempt from equality; it is never
rewritten to pretend that old cohorts ran under later input/program bytes.
Per-ensemble JSON records and the parsed-input snapshot are written to the
unique run directory as it progresses. Those records preserve work after an
interruption; they are never read as simulation inputs by the Verifier.
An earlier complete trajectory allocation failed at final JSON serialization
of a NumPy boolean; its rounded logs supply no final reference measurements.

## Reproduced reference measurements

<!-- BEGIN REPRODUCED NUMBERS -->

All +/- values are one Monte Carlo standard error. Log budgets are conservative empirical diagnostics, not rigorous continuum-error bounds.

| Case | Pair | Reference / (D_unit length_unit) | Reference (m3 mol-1 s-1) | stat | step | contact | adaptive | box | plateau | Total log | Factor | Qualified |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| N4 | end/end | 51.7866107 +/- 0.4787583 | 117927799 +/- 1090222 | 0.03698 | 0.05302 | 0.04462 | 0.04258 | 0.05839 | 0.09622 | 0.33182 | 1.39350 | YES |
| N4 | end/mid | 49.144749 +/- 0.4713032 | 111911786 +/- 1073246 | 0.03836 | 0.05713 | 0.04752 | 0.05221 | 0.05943 | 0.09533 | 0.34997 | 1.41903 | YES |
| N4 | mid/mid | 45.9045185 +/- 0.461476 | 104533175 +/- 1050867 | 0.04021 | 0.05706 | 0.04491 | 0.05090 | 0.07245 | 0.08357 | 0.34910 | 1.41779 | YES |
| N8 | end/end | 50.1519856 +/- 0.4735122 | 80755449.4 +/- 762456.2 | 0.03777 | 0.07012 | 0.04463 | 0.04359 | 0.06863 | 0.07792 | 0.34266 | 1.40868 | YES |
| N8 | end/mid | 45.2449363 +/- 0.4588953 | 72854048.1 +/- 738919.8 | 0.04057 | 0.06470 | 0.06424 | 0.06038 | 0.07111 | 0.08744 | 0.38844 | 1.47468 | YES |
| N8 | mid/mid | 39.7007607 +/- 0.4399938 | 63926736.6 +/- 708484.3 | 0.04433 | 0.07467 | 0.05151 | 0.05279 | 0.07893 | 0.09946 | 0.40169 | 1.49435 | YES |
| N16 | end/end | 47.4820351 +/- 0.4665417 | 54062737.6 +/- 531201.4 | 0.03930 | 0.07020 | 0.04330 | 0.04853 | 0.05470 | 0.10909 | 0.36512 | 1.44069 | YES |
| N16 | end/mid | 41.9034732 +/- 0.4479025 | 47711023.2 +/- 509978.9 | 0.04276 | 0.06528 | 0.05192 | 0.05130 | 0.06731 | 0.08903 | 0.36759 | 1.44426 | YES |
| N16 | mid/mid | 34.6808847 +/- 0.4195076 | 39487430.7 +/- 477648.6 | 0.04838 | 0.08557 | 0.06040 | 0.06366 | 0.08296 | 0.11687 | 0.45785 | 1.58067 | YES |
| N4_N16 | end/end | 40.5070045 +/- 0.3174643 | 92242025.9 +/- 722925.5 | 0.03135 | 0.04345 | 0.03469 | 0.04785 | 0.06380 | 0.08923 | 0.31037 | 1.36393 | YES |
| N4_N16 | end/mid | 36.6709804 +/- 0.3092471 | 83506681.4 +/- 704213.6 | 0.03373 | 0.04522 | 0.04066 | 0.03847 | 0.06225 | 0.07486 | 0.29519 | 1.34338 | YES |
| N4_N16 | mid/mid | 33.5150814 +/- 0.3025236 | 76320109.2 +/- 688902.9 | 0.03611 | 0.06501 | 0.05192 | 0.04830 | 0.05491 | 0.10902 | 0.36527 | 1.44091 | YES |
| N16_N4 | end/end | 41.063395 +/- 0.3181848 | 93509031.1 +/- 724566.3 | 0.03099 | 0.04401 | 0.04100 | 0.03276 | 0.06033 | 0.07242 | 0.28152 | 1.32514 | YES |
| N16_N4 | end/mid | 38.516392 +/- 0.3135843 | 87709028.8 +/- 714090.1 | 0.03257 | 0.04120 | 0.03891 | 0.04337 | 0.06229 | 0.08150 | 0.29984 | 1.34964 | YES |
| N16_N4 | mid/mid | 34.2030772 +/- 0.3035654 | 77886804.3 +/- 691275.2 | 0.03550 | 0.05009 | 0.04990 | 0.03790 | 0.05884 | 0.07399 | 0.30622 | 1.35828 | YES |
| floor | end/end | 25.2228164 +/- 0.2271163 | 283178726 +/- 2549854 | 0.03602 | 0.06040 | 0.04209 | 0.04919 | 0.06587 | 0.09298 | 0.34654 | 1.41417 | YES |
| floor | end/mid | 25.2228164 +/- 0.2271163 | 283178726 +/- 2549854 | 0.03602 | 0.06040 | 0.04209 | 0.04919 | 0.06587 | 0.09298 | 0.34654 | 1.41417 | YES |
| floor | mid/mid | 25.2228164 +/- 0.2271163 | 283178726 +/- 2549854 | 0.03602 | 0.06040 | 0.04209 | 0.04919 | 0.06587 | 0.09298 | 0.34654 | 1.41417 | YES |

| Case | Pair | Run | Events | k_box +/- SE | k_infinite +/- SE | Correction removed | Late half 1 | Late half 2 | Mixing times at start |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| N4 | end/end | base | 5192 | 71.5603 +/- 0.993975 | 50.97138 +/- 0.504294 | 0.28771 | 49.93027 +/- 0.70136 | 51.99559 +/- 0.72381 | 4.73741 |
| N4 | end/end | step | 10366 | 71.108915 +/- 0.699012 | 50.741953 +/- 0.355935 | 0.28642 | 51.23528 +/- 0.49396 | 50.24475 +/- 0.51285 | 4.73741 |
| N4 | end/end | contact | 10415 | 71.612023 +/- 0.702308 | 50.997617 +/- 0.356168 | 0.28786 | 51.14488 +/- 0.49449 | 50.85001 +/- 0.51283 | 4.73741 |
| N4 | end/end | box | 8678 | 64.785271 +/- 0.695596 | 51.367713 +/- 0.437306 | 0.20711 | 51.53823 +/- 0.61319 | 51.1969 +/- 0.62369 | 4.83409 |
| N4 | end/end | box2 | 8210 | 61.82721 +/- 0.682403 | 51.786611 +/- 0.478758 | 0.16240 | 52.36121 +/- 0.67506 | 51.20953 +/- 0.67909 | 4.87388 |
| N4 | end/end | adaptive | 10393 | 71.307987 +/- 0.700061 | 50.843239 +/- 0.355898 | 0.28699 | 50.63559 +/- 0.49442 | 51.05021 +/- 0.51194 | 4.73741 |
| N4 | end/mid | base | 5001 | 67.28025 +/- 0.952109 | 48.761871 +/- 0.500118 | 0.27524 | 48.04287 +/- 0.6957 | 49.47291 +/- 0.71804 | 4.73741 |
| N4 | end/mid | step | 9918 | 66.655151 +/- 0.669797 | 48.432681 +/- 0.353633 | 0.27338 | 48.60144 +/- 0.49164 | 48.26348 +/- 0.50853 | 4.73741 |
| N4 | end/mid | contact | 10005 | 67.241149 +/- 0.67275 | 48.741329 +/- 0.353491 | 0.27513 | 49.29321 +/- 0.49116 | 48.18467 +/- 0.50875 | 4.73741 |
| N4 | end/mid | box | 8235 | 60.739818 +/- 0.669455 | 48.791105 +/- 0.431972 | 0.19672 | 49.03977 +/- 0.60625 | 48.54182 +/- 0.61556 | 4.83409 |
| N4 | end/mid | box2 | 7781 | 58.098489 +/- 0.658682 | 49.144749 +/- 0.471303 | 0.15411 | 49.60079 +/- 0.66461 | 48.68716 +/- 0.66844 | 4.87388 |
| N4 | end/mid | adaptive | 9836 | 66.237969 +/- 0.668367 | 48.212044 +/- 0.354089 | 0.27214 | 48.36313 +/- 0.49234 | 48.0606 +/- 0.50911 | 4.73741 |
| N4 | mid/mid | base | 4610 | 60.697227 +/- 0.89451 | 45.208281 +/- 0.49623 | 0.25518 | 44.87233 +/- 0.691 | 45.54253 +/- 0.71215 | 4.73741 |
| N4 | mid/mid | step | 9293 | 60.97282 +/- 0.632889 | 45.36099 +/- 0.350284 | 0.25605 | 45.01783 +/- 0.48774 | 45.70237 +/- 0.50272 | 4.73741 |
| N4 | mid/mid | contact | 9285 | 61.073462 +/- 0.634208 | 45.416669 +/- 0.350717 | 0.25636 | 44.95502 +/- 0.48832 | 45.8751 +/- 0.50329 | 4.73741 |
| N4 | mid/mid | box | 7587 | 55.141475 +/- 0.633154 | 45.112015 +/- 0.423776 | 0.18189 | 45.55222 +/- 0.59578 | 44.66989 +/- 0.60286 | 4.83409 |
| N4 | mid/mid | box2 | 7252 | 53.623772 +/- 0.629728 | 45.904518 +/- 0.461476 | 0.14395 | 45.97669 +/- 0.64988 | 45.83231 +/- 0.65537 | 4.87388 |
| N4 | mid/mid | adaptive | 9211 | 60.495482 +/- 0.630717 | 45.096268 +/- 0.350485 | 0.25455 | 45.48409 +/- 0.48813 | 44.70616 +/- 0.50323 | 4.73741 |
| N8 | end/end | base | 4905 | 65.97605 +/- 0.942718 | 48.073134 +/- 0.500512 | 0.27135 | 47.7957 +/- 0.6962 | 48.34938 +/- 0.71906 | 4.73741 |
| N8 | end/end | step | 10044 | 67.766335 +/- 0.676695 | 49.016691 +/- 0.354041 | 0.27668 | 49.4668 +/- 0.4919 | 48.5634 +/- 0.50954 | 4.73741 |
| N8 | end/end | contact | 9975 | 67.423758 +/- 0.675594 | 48.837207 +/- 0.354455 | 0.27567 | 48.57714 +/- 0.49288 | 49.09622 +/- 0.50939 | 4.73741 |
| N8 | end/end | box | 8338 | 61.535589 +/- 0.674027 | 49.303264 +/- 0.432689 | 0.19878 | 49.0935 +/- 0.60639 | 49.51259 +/- 0.61734 | 4.83409 |
| N8 | end/end | box2 | 7968 | 59.511455 +/- 0.666739 | 50.151986 +/- 0.473512 | 0.15727 | 50.092 +/- 0.66604 | 50.21194 +/- 0.67323 | 4.87388 |
| N8 | end/end | adaptive | 9951 | 67.192816 +/- 0.674087 | 48.715928 +/- 0.354334 | 0.27498 | 48.6335 +/- 0.49265 | 48.79825 +/- 0.50937 | 4.73741 |
| N8 | end/mid | base | 4517 | 59.041025 +/- 0.878984 | 44.283059 +/- 0.49448 | 0.24996 | 44.73561 +/- 0.68902 | 43.82741 +/- 0.70966 | 4.73741 |
| N8 | end/mid | step | 8944 | 58.272298 +/- 0.616513 | 43.849193 +/- 0.349093 | 0.24751 | 44.31096 +/- 0.48657 | 43.38421 +/- 0.50088 | 4.73741 |
| N8 | end/mid | contact | 9154 | 59.815394 +/- 0.625556 | 44.717263 +/- 0.349616 | 0.25241 | 44.65989 +/- 0.48697 | 44.77458 +/- 0.50175 | 4.73741 |
| N8 | end/mid | box | 7494 | 54.297065 +/- 0.627312 | 44.545263 +/- 0.422215 | 0.17960 | 44.42225 +/- 0.59236 | 44.66813 +/- 0.60178 | 4.83409 |
| N8 | end/mid | box2 | 7159 | 52.725879 +/- 0.623191 | 45.244936 +/- 0.458895 | 0.14188 | 45.38739 +/- 0.64656 | 45.10233 +/- 0.65139 | 4.87388 |
| N8 | end/mid | adaptive | 8981 | 58.566247 +/- 0.618348 | 44.01543 +/- 0.34926 | 0.24845 | 43.96412 +/- 0.48663 | 44.0667 +/- 0.5011 | 4.73741 |
| N8 | mid/mid | base | 3997 | 50.312183 +/- 0.79614 | 39.184146 +/- 0.482907 | 0.22118 | 39.34036 +/- 0.67453 | 39.02757 +/- 0.69131 | 4.73741 |
| N8 | mid/mid | step | 8150 | 51.264066 +/- 0.5681 | 39.759115 +/- 0.341721 | 0.22443 | 39.21648 +/- 0.47647 | 40.2975 +/- 0.48972 | 4.73741 |
| N8 | mid/mid | contact | 8123 | 51.077538 +/- 0.566971 | 39.646824 +/- 0.3416 | 0.22379 | 39.28056 +/- 0.4765 | 40.01114 +/- 0.48944 | 4.73741 |
| N8 | mid/mid | box | 6515 | 46.265209 +/- 0.573249 | 38.991846 +/- 0.407176 | 0.15721 | 38.88394 +/- 0.57182 | 39.09964 +/- 0.5798 | 4.83409 |
| N8 | mid/mid | box2 | 6241 | 45.346267 +/- 0.574027 | 39.700761 +/- 0.439994 | 0.12450 | 39.91489 +/- 0.62079 | 39.4863 +/- 0.6237 | 4.87388 |
| N8 | mid/mid | adaptive | 8151 | 51.34957 +/- 0.569013 | 39.810528 +/- 0.342014 | 0.22472 | 40.17248 +/- 0.47778 | 39.44665 +/- 0.48962 | 4.73741 |
| N16 | end/end | base | 4745 | 63.215217 +/- 0.918317 | 46.590507 +/- 0.498821 | 0.26299 | 47.46209 +/- 0.69395 | 45.70712 +/- 0.71717 | 4.73741 |
| N16 | end/end | step | 9708 | 64.788069 +/- 0.658012 | 47.43931 +/- 0.352794 | 0.26778 | 47.72388 +/- 0.4907 | 47.15349 +/- 0.50715 | 4.73741 |
| N16 | end/end | contact | 9696 | 64.680961 +/- 0.657329 | 47.381858 +/- 0.35274 | 0.26745 | 47.45539 +/- 0.4907 | 47.30825 +/- 0.5069 | 4.73741 |
| N16 | end/end | box | 8005 | 58.819186 +/- 0.657526 | 47.544042 +/- 0.429603 | 0.19169 | 47.65115 +/- 0.60286 | 47.43682 +/- 0.61223 | 4.83409 |
| N16 | end/end | box2 | 7504 | 55.788956 +/- 0.644063 | 47.482035 +/- 0.466542 | 0.14890 | 48.2031 +/- 0.65902 | 46.75711 +/- 0.66055 | 4.87388 |
| N16 | end/end | adaptive | 9778 | 65.262655 +/- 0.660462 | 47.693261 +/- 0.352722 | 0.26921 | 47.37527 +/- 0.49071 | 48.0097 +/- 0.50664 | 4.73741 |
| N16 | end/mid | base | 4290 | 54.936994 +/- 0.839179 | 41.933478 +/- 0.48893 | 0.23670 | 42.29236 +/- 0.68218 | 41.57268 +/- 0.70076 | 4.73741 |
| N16 | end/mid | step | 8503 | 54.363203 +/- 0.589838 | 41.598342 +/- 0.345362 | 0.23481 | 41.1021 +/- 0.48138 | 42.09097 +/- 0.49511 | 4.73741 |
| N16 | end/mid | contact | 8432 | 54.022127 +/- 0.588596 | 41.398341 +/- 0.345653 | 0.23368 | 41.5123 +/- 0.48226 | 41.28419 +/- 0.49535 | 4.73741 |
| N16 | end/mid | box | 6977 | 49.888911 +/- 0.597343 | 41.534434 +/- 0.41403 | 0.16746 | 41.47429 +/- 0.58135 | 41.59455 +/- 0.58967 | 4.83409 |
| N16 | end/mid | box2 | 6604 | 48.242832 +/- 0.593675 | 41.903473 +/- 0.447903 | 0.13141 | 41.97716 +/- 0.63111 | 41.82975 +/- 0.63574 | 4.87388 |
| N16 | end/mid | adaptive | 8520 | 54.320882 +/- 0.58879 | 41.573558 +/- 0.344875 | 0.23467 | 41.54609 +/- 0.48104 | 41.60101 +/- 0.49431 | 4.73741 |
| N16 | mid/mid | base | 3529 | 43.171247 +/- 0.726949 | 34.712352 +/- 0.469983 | 0.19594 | 34.41779 +/- 0.65661 | 35.0057 +/- 0.67246 | 4.73741 |
| N16 | mid/mid | step | 7246 | 44.240883 +/- 0.519896 | 35.400547 +/- 0.332881 | 0.19982 | 35.48552 +/- 0.46568 | 35.31547 +/- 0.47581 | 4.73741 |
| N16 | mid/mid | contact | 7303 | 44.650194 +/- 0.522657 | 35.662138 +/- 0.333415 | 0.20130 | 35.56058 +/- 0.46605 | 35.76355 +/- 0.4769 | 4.73741 |
| N16 | mid/mid | box | 5902 | 41.129792 +/- 0.535419 | 35.279403 +/- 0.393934 | 0.14224 | 35.06206 +/- 0.55308 | 35.4963 +/- 0.56107 | 4.83409 |
| N16 | mid/mid | box2 | 5429 | 38.91289 +/- 0.528137 | 34.680885 +/- 0.419508 | 0.10876 | 34.33191 +/- 0.58937 | 35.02901 +/- 0.59711 | 4.87388 |
| N16 | mid/mid | adaptive | 7201 | 44.066736 +/- 0.519463 | 35.288956 +/- 0.333128 | 0.19919 | 35.32963 +/- 0.46597 | 35.24826 +/- 0.47621 | 4.73741 |
| N4_N16 | end/end | base | 6369 | 62.088381 +/- 0.779272 | 39.781249 +/- 0.319908 | 0.35928 | 39.7207 +/- 0.44125 | 39.84169 +/- 0.46324 | 4.73741 |
| N4_N16 | end/end | step | 12758 | 61.698922 +/- 0.547131 | 39.621006 +/- 0.225625 | 0.35783 | 39.82979 +/- 0.31062 | 39.41099 +/- 0.3275 | 4.73741 |
| N4_N16 | end/end | contact | 12809 | 61.94285 +/- 0.548207 | 39.721455 +/- 0.22543 | 0.35874 | 39.4873 +/- 0.31136 | 39.95407 +/- 0.32588 | 4.73741 |
| N4_N16 | end/end | box | 11085 | 55.44751 +/- 0.526847 | 40.83954 +/- 0.285813 | 0.26346 | 41.01916 +/- 0.39923 | 40.65935 +/- 0.40916 | 4.83409 |
| N4_N16 | end/end | box2 | 10338 | 50.839759 +/- 0.500082 | 40.507004 +/- 0.317464 | 0.20324 | 41.04179 +/- 0.44654 | 39.96859 +/- 0.45141 | 4.87388 |
| N4_N16 | end/end | adaptive | 13054 | 63.518872 +/- 0.556901 | 40.363675 +/- 0.224882 | 0.36454 | 40.4251 +/- 0.3097 | 40.30215 +/- 0.3262 | 4.73741 |
| N4_N16 | end/mid | base | 5829 | 52.812004 +/- 0.692552 | 35.757083 +/- 0.317476 | 0.32294 | 35.83549 +/- 0.43922 | 35.67851 +/- 0.45861 | 4.73741 |
| N4_N16 | end/mid | step | 11672 | 52.948052 +/- 0.490678 | 35.819397 +/- 0.22456 | 0.32350 | 36.0472 +/- 0.31041 | 35.5902 +/- 0.32474 | 4.73741 |
| N4_N16 | end/mid | contact | 11749 | 53.364659 +/- 0.492925 | 36.009575 +/- 0.224445 | 0.32522 | 36.14131 +/- 0.31035 | 35.87737 +/- 0.32441 | 4.73741 |
| N4_N16 | end/mid | box | 9753 | 47.017143 +/- 0.476222 | 36.075248 +/- 0.28036 | 0.23272 | 35.76363 +/- 0.39218 | 36.38524 +/- 0.40064 | 4.83409 |
| N4_N16 | end/mid | box2 | 9365 | 44.939628 +/- 0.464429 | 36.67098 +/- 0.309247 | 0.18399 | 36.5352 +/- 0.43436 | 36.80653 +/- 0.44028 | 4.87388 |
| N4_N16 | end/mid | adaptive | 11727 | 53.114893 +/- 0.491072 | 35.895675 +/- 0.224283 | 0.32419 | 36.05485 +/- 0.31012 | 35.73582 +/- 0.32421 | 4.73741 |
| N4_N16 | mid/mid | base | 5412 | 47.327918 +/- 0.643952 | 33.155865 +/- 0.316038 | 0.29944 | 33.12659 +/- 0.43845 | 33.18512 +/- 0.45526 | 4.73741 |
| N4_N16 | mid/mid | step | 11003 | 48.603982 +/- 0.463825 | 33.777115 +/- 0.224004 | 0.30505 | 33.57031 +/- 0.31082 | 33.98282 +/- 0.32252 | 4.73741 |
| N4_N16 | mid/mid | contact | 10840 | 47.627163 +/- 0.457889 | 33.302451 +/- 0.223874 | 0.30077 | 33.48094 +/- 0.31032 | 33.12313 +/- 0.32287 | 4.73741 |
| N4_N16 | mid/mid | box | 9117 | 43.096499 +/- 0.45146 | 33.721422 +/- 0.276405 | 0.21754 | 33.8001 +/- 0.38728 | 33.64265 +/- 0.3945 | 4.83409 |
| N4_N16 | mid/mid | box2 | 8494 | 40.290299 +/- 0.4372 | 33.515081 +/- 0.302524 | 0.16816 | 32.8962 +/- 0.42383 | 34.12938 +/- 0.43165 | 4.87388 |
| N4_N16 | mid/mid | adaptive | 11017 | 48.350367 +/- 0.461106 | 33.654436 +/- 0.223402 | 0.30395 | 33.47758 +/- 0.30997 | 33.83049 +/- 0.32168 | 4.73741 |
| N16_N4 | end/end | base | 6431 | 62.479748 +/- 0.780411 | 39.94155 +/- 0.318929 | 0.36073 | 40.26567 +/- 0.43852 | 39.61443 +/- 0.4636 | 4.73741 |
| N16_N4 | end/end | step | 13000 | 62.970976 +/- 0.553227 | 40.141732 +/- 0.224809 | 0.36254 | 40.45141 +/- 0.30907 | 39.82931 +/- 0.32682 | 4.73741 |
| N16_N4 | end/end | contact | 12828 | 62.077261 +/- 0.548993 | 39.776683 +/- 0.225403 | 0.35924 | 39.63207 +/- 0.3111 | 39.92071 +/- 0.32613 | 4.73741 |
| N16_N4 | end/end | box | 11062 | 55.28992 +/- 0.525894 | 40.753984 +/- 0.285724 | 0.26290 | 40.65225 +/- 0.39918 | 40.85554 +/- 0.40888 | 4.83409 |
| N16_N4 | end/end | box2 | 10502 | 51.719289 +/- 0.504748 | 41.063395 +/- 0.318185 | 0.20603 | 40.84888 +/- 0.44649 | 41.27733 +/- 0.4534 | 4.87388 |
| N16_N4 | end/end | adaptive | 12791 | 62.012165 +/- 0.549208 | 39.749946 +/- 0.225661 | 0.35900 | 39.28233 +/- 0.31221 | 40.21148 +/- 0.32549 | 4.73741 |
| N16_N4 | end/mid | base | 6186 | 58.327079 +/- 0.742669 | 38.202797 +/- 0.3186 | 0.34502 | 38.1686 +/- 0.44004 | 38.23697 +/- 0.46082 | 4.73741 |
| N16_N4 | end/mid | step | 12374 | 58.296871 +/- 0.524831 | 38.189836 +/- 0.225229 | 0.34491 | 37.96108 +/- 0.31148 | 38.41716 +/- 0.32523 | 4.73741 |
| N16_N4 | end/mid | contact | 12291 | 57.81487 +/- 0.522234 | 37.982396 +/- 0.225398 | 0.34303 | 38.04928 +/- 0.31117 | 37.91539 +/- 0.32623 | 4.73741 |
| N16_N4 | end/mid | box | 10581 | 51.990548 +/- 0.505604 | 38.932829 +/- 0.283527 | 0.25116 | 38.77958 +/- 0.39636 | 39.08567 +/- 0.40547 | 4.83409 |
| N16_N4 | end/mid | box2 | 9821 | 47.742887 +/- 0.481815 | 38.516392 +/- 0.313584 | 0.19325 | 38.20056 +/- 0.44002 | 38.83099 +/- 0.44684 | 4.87388 |
| N16_N4 | end/mid | adaptive | 12478 | 58.705392 +/- 0.526313 | 38.364728 +/- 0.224777 | 0.34649 | 37.96523 +/- 0.31116 | 38.75986 +/- 0.32413 | 4.73741 |
| N16_N4 | mid/mid | base | 5378 | 47.216112 +/- 0.644455 | 33.100954 +/- 0.316733 | 0.29895 | 33.16549 +/- 0.43929 | 33.03631 +/- 0.45645 | 4.73741 |
| N16_N4 | mid/mid | step | 10821 | 47.438112 +/- 0.456468 | 33.209908 +/- 0.223713 | 0.29993 | 33.54101 +/- 0.30995 | 32.87595 +/- 0.32286 | 4.73741 |
| N16_N4 | mid/mid | contact | 10957 | 48.261209 +/- 0.461513 | 33.611216 +/- 0.223849 | 0.30356 | 33.74823 +/- 0.31024 | 33.47371 +/- 0.32286 | 4.73741 |
| N16_N4 | mid/mid | box | 9273 | 43.855316 +/- 0.455532 | 34.184234 +/- 0.276774 | 0.22052 | 34.15312 +/- 0.38762 | 34.21533 +/- 0.39518 | 4.83409 |
| N16_N4 | mid/mid | box2 | 8713 | 41.288717 +/- 0.442369 | 34.203077 +/- 0.303565 | 0.17161 | 34.15203 +/- 0.42674 | 34.25409 +/- 0.43185 | 4.87388 |
| N16_N4 | mid/mid | adaptive | 10953 | 48.246246 +/- 0.461454 | 33.603958 +/- 0.223863 | 0.30349 | 33.45066 +/- 0.3106 | 33.75665 +/- 0.32237 | 4.73741 |
| floor | end/end | base | 6670 | 28.716808 +/- 0.351668 | 24.711227 +/- 0.260405 | 0.13949 | 24.50128 +/- 0.36478 | 24.9206 +/- 0.37168 | 4.73741 |
| floor | end/end | step | 13487 | 29.015373 +/- 0.24988 | 24.93199 +/- 0.184497 | 0.14073 | 24.80377 +/- 0.25867 | 25.06 +/- 0.26313 | 4.73741 |
| floor | end/end | contact | 13494 | 29.007407 +/- 0.249747 | 24.926108 +/- 0.184413 | 0.14070 | 25.00592 +/- 0.25912 | 24.84621 +/- 0.26248 | 4.73741 |
| floor | end/end | box | 10969 | 27.58056 +/- 0.263352 | 24.820488 +/- 0.21328 | 0.10007 | 24.47444 +/- 0.29906 | 25.16547 +/- 0.30413 | 4.83409 |
| floor | end/end | box2 | 10460 | 27.3892 +/- 0.267806 | 25.222816 +/- 0.227116 | 0.07910 | 25.48658 +/- 0.32168 | 24.95858 +/- 0.32069 | 4.87388 |
| floor | end/end | adaptive | 13585 | 29.259446 +/- 0.251072 | 25.111986 +/- 0.184939 | 0.14175 | 24.82632 +/- 0.25884 | 25.39658 +/- 0.26417 | 4.73741 |
| floor | end/mid | base | 6670 | 28.716808 +/- 0.351668 | 24.711227 +/- 0.260405 | 0.13949 | 24.50128 +/- 0.36478 | 24.9206 +/- 0.37168 | 4.73741 |
| floor | end/mid | step | 13487 | 29.015373 +/- 0.24988 | 24.93199 +/- 0.184497 | 0.14073 | 24.80377 +/- 0.25867 | 25.06 +/- 0.26313 | 4.73741 |
| floor | end/mid | contact | 13494 | 29.007407 +/- 0.249747 | 24.926108 +/- 0.184413 | 0.14070 | 25.00592 +/- 0.25912 | 24.84621 +/- 0.26248 | 4.73741 |
| floor | end/mid | box | 10969 | 27.58056 +/- 0.263352 | 24.820488 +/- 0.21328 | 0.10007 | 24.47444 +/- 0.29906 | 25.16547 +/- 0.30413 | 4.83409 |
| floor | end/mid | box2 | 10460 | 27.3892 +/- 0.267806 | 25.222816 +/- 0.227116 | 0.07910 | 25.48658 +/- 0.32168 | 24.95858 +/- 0.32069 | 4.87388 |
| floor | end/mid | adaptive | 13585 | 29.259446 +/- 0.251072 | 25.111986 +/- 0.184939 | 0.14175 | 24.82632 +/- 0.25884 | 25.39658 +/- 0.26417 | 4.73741 |
| floor | mid/mid | base | 6670 | 28.716808 +/- 0.351668 | 24.711227 +/- 0.260405 | 0.13949 | 24.50128 +/- 0.36478 | 24.9206 +/- 0.37168 | 4.73741 |
| floor | mid/mid | step | 13487 | 29.015373 +/- 0.24988 | 24.93199 +/- 0.184497 | 0.14073 | 24.80377 +/- 0.25867 | 25.06 +/- 0.26313 | 4.73741 |
| floor | mid/mid | contact | 13494 | 29.007407 +/- 0.249747 | 24.926108 +/- 0.184413 | 0.14070 | 25.00592 +/- 0.25912 | 24.84621 +/- 0.26248 | 4.73741 |
| floor | mid/mid | box | 10969 | 27.58056 +/- 0.263352 | 24.820488 +/- 0.21328 | 0.10007 | 24.47444 +/- 0.29906 | 25.16547 +/- 0.30413 | 4.83409 |
| floor | mid/mid | box2 | 10460 | 27.3892 +/- 0.267806 | 25.222816 +/- 0.227116 | 0.07910 | 25.48658 +/- 0.32168 | 24.95858 +/- 0.32069 | 4.87388 |
| floor | mid/mid | adaptive | 13585 | 29.259446 +/- 0.251072 | 25.111986 +/- 0.184939 | 0.14175 | 24.82632 +/- 0.25884 | 25.39658 +/- 0.26417 | 4.73741 |

| Pair | Bead contacts | Mode contacts | Replicas | Combined SE of probabilities | Allowed difference | Pass |
|---|---:|---:|---:|---:|---:|---|
| end/end | 4173 | 4093 | 32768 | 0.002593692 | 0.01373256 | True |
| end/mid | 4072 | 4035 | 32768 | 0.002572203 | 0.01360662 | True |
| mid/mid | 3932 | 3878 | 32768 | 0.002531162 | 0.01337578 | True |

Sampled covariance: 75 checks; maximum absolute discrepancy / SE = 3.0215409 (limit 6).
Independent Ewald constant: 2.83729747948.

<!-- END REPRODUCED NUMBERS -->

## Actual candidate comparison

Target rates include Avogadro's constant and have units m3 mol-1 s-1. The
reference volume rate is multiplied by the literal constant once. There is no
identical-reactant consumption factor. The candidate's analytically exponential
single-pair survival makes a separate kMC sampling run unnecessary.

<!-- BEGIN CANDIDATE NUMBERS -->

| Case | Pair | Candidate (m3 mol-1 s-1) | Candidate / reference | Within numerical budget | Reference qualified | Accepted |
|---|---|---:|---:|---|---|---|
| N4 | end/end | 114463905 | 0.970627 | True | True | True |
| N4 | end/mid | 114463905 | 1.0228047 | True | True | True |
| N4 | mid/mid | 114463905 | 1.0950008 | True | True | True |
| N8 | end/end | 80938203.7 | 1.0022631 | True | True | True |
| N8 | end/mid | 80938203.7 | 1.1109637 | True | True | True |
| N8 | mid/mid | 80938203.7 | 1.2661088 | True | True | True |
| N16 | end/end | 57231952.7 | 1.0586211 | True | True | True |
| N16 | end/mid | 57231952.7 | 1.1995541 | True | True | True |
| N16 | mid/mid | 57231952.7 | 1.4493714 | True | True | True |
| N4_N16 | end/end | 71539940.9 | 0.77556775 | True | True | True |
| N4_N16 | end/mid | 71539940.9 | 0.85669721 | True | True | True |
| N4_N16 | mid/mid | 71539940.9 | 0.93736686 | True | True | True |
| N16_N4 | end/end | 71539940.9 | 0.76505916 | True | True | True |
| N16_N4 | end/mid | 71539940.9 | 0.81565082 | True | True | True |
| N16_N4 | mid/mid | 71539940.9 | 0.9185117 | True | True | True |
| floor | end/end | 282167444 | 0.99642882 | True | True | True |
| floor | end/mid | 282167444 | 0.99642882 | True | True | True |
| floor | mid/mid | 282167444 | 0.99642882 | True | True | True |

Literal transport anchors: 13/13.
Literal min/floor/two-chain branches: 18/18.
Proposed scientific test passed: True.

<!-- END CANDIDATE NUMBERS -->

## Mutation results and their controls

Actual mutations compile modified target source **in memory** and leave the
product file untouched. The unmodified target must first pass literal transport
and `min`/floor/two-chain rule assertions. Failure of already-rejected scientific
adoption is never counted as a mutation kill. Mutation tables state the lane
that rejects a mutant. Exact rule assertions alone do not validate the Rouse
physical approximation. Separately, qualified oracle-backed positive controls
measure rate-budget discrimination without relying on a rejected target.
The six requested rate mutations must each reject at least one qualified row
whose unmodified target rate passed; the strict mutation flag enforces this
in addition to the literal-rule controls.

<!-- BEGIN MUTATION NUMBERS -->

Target baseline passes literal-input and kernel-rule controls before each mutation. Scientific adoption failure is never counted as a mutation kill.

| Actual target mutation | Changes observable | Rejected | Failed anchors | Failed literal branches | Failed qualified rate checks | New rate rejections from passing baseline rows |
|---|---|---|---:|---:|---:|---:|---:|
| capture_radius_x0.5 | True | True | 0 | 18 | 17 | 17 |
| capture_radius_x2 | True | True | 0 | 18 | 18 | 18 |
| candidate_diffusion_x4 | True | True | 0 | 18 | 18 | 18 |
| shared_chain_diffusivity_x4 | True | True | 4 | 18 | 18 | 18 |
| drop_one_chain_diffusion | True | True | 0 | 18 | 16 | 16 |
| length_min_to_max | True | True | 0 | 6 | 6 | 6 |
| length_min_to_i | True | True | 0 | 3 | 3 | 3 |
| length_min_to_j | True | True | 0 | 3 | 3 | 3 |
| length_min_to_mean | True | True | 0 | 6 | 3 | 3 |
| drop_sigma0_floor | True | True | 0 | 3 | 0 | 0 |
| swap_reference_prefactors_in_candidate | True | True | 0 | 10 | 4 | 4 |

| Rate control against qualified oracle | Positive controls | Rejected by rate budget |
|---|---:|---|
| radius_x0.5 | 18 | True |
| radius_x2 | 18 | True |
| diffusion_x4 | 18 | True |
| drop_one_equal_chain | 12 | True |

| Oracle prefactor swap witness | Original end > mid by 4 SE | Swap rejected |
|---|---|---|
| N4 | True | True |
| N8 | True | True |
| N16 | True | True |
| N4_N16 | True | True |
| N16_N4 | True | True |
| floor | False | False |

Actual sampler noise/coordinates x2 rejected: True.

The required prefactor mutation rescales candidate end/end by C_mid/C_end and mid/mid by C_end/C_mid, using the independent reference. A bare target label swap is a separate observational no-op; it is reported as an expected survivor, never substituted for the requested prefactor perturbation.

<!-- END MUTATION NUMBERS -->

## Executable proposed scientific test and Verifier

This is a proposed fixture-consuming test; no repository test file is changed.
The verifier extracts and executes this exact body. A failed oracle or candidate
causes a real assertion failure, preserved in JUnit and the terminal log.

```python
def test_met_kernel_reference(tmp_path):
    root = Path(os.environ["MET_KERNEL_REFERENCE_PACK"])
    sys.path.insert(0, str(root))
    import checks
    output = tmp_path / "fresh-reference"
    subprocess.run([
        sys.executable, "-B", "-u", str(root / "reference/run.py"),
        "--workers", os.environ.get("MET_KERNEL_REFERENCE_WORKERS", "8"),
        "--output", str(output),
    ], check=True)
    reference = json.loads((output / "results.json").read_text())
    import importlib.util
    spec = importlib.util.spec_from_file_location("test_rouse_inputs", root / "reference/run.py")
    oracle = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(oracle)
    cfg, digest = oracle.read_parameters()
    assert reference["parameters_sha256"] == digest
    report = checks.target_checks(checks.load_target(), cfg, reference)
    checks.assert_scientific_acceptance(report)
```

From the named worktree run:

```bash
taskset -c 6,16,17,19,21,23,24,26 /home/alon/anaconda3/envs/rmg_env/bin/python -B -u test/rmgpy/kmc/fixtures/i030_kernel_reference/verify_pack.py --workers 8 --require-all-mutations --require-acceptance --output /home/alon/runs/i030-met-kernel-reference/rework/verifications > >(tee -a /home/alon/runs/i030-met-kernel-reference/rework/verifier-final-stdout.log) 2> >(tee -a /home/alon/runs/i030-met-kernel-reference/rework/verifier-final-stderr.log >&2)
```

The default exit status checks **report reproduction**: it reruns every physical
ensemble, sample covariance, bead contacts, sampler audits and mutation controls,
checks every deterministic scientific JSON value/contact hash and exact rendered table,
validates and retains both actual generation histories,
and classifies the real adoption pytest outcome. It is not an adoption pass.
`--require-acceptance` additionally exits nonzero for a failed adoption test.
`--mutations-only --require-all-mutations` exits nonzero if any requested target
substantive mutation survives. The extra label-only no-op is not a required
perturbation; the substantive reference-prefactor swap is.
`mutation_tests.py --require-all` exposes the same strict mutation requirement.
There is no expected-failure decorator or cached simulation in full verification.

Every verification has a UUID-bearing private staging directory. Its complete
results, runtime, candidate, mutations, audits, version/hash manifest and pytest
output are published together by one same-filesystem atomic directory rename.
No shared `reproduced/` files or mutable latest pointer exists. A failed run's
pending staging directory remains distinguishable from a completed publication.
`COMPLETE.json` hashes every file in the publication. Concurrent runs cannot
mix their JSON and Markdown artifacts. Verifier outputs never write into the
fixture. `--assemble` is an explicit pack-authoring helper.

## Finding-by-finding response to round 1

1. **P1-1, shared physical inputs:** replaced runtime target mapping with cited
   literal JSON inputs and a pre-simulation hash. The oracle has no target access.
   Separate literal transport/Rg assertions kill the shared diffusivity x4 bug.
   Original benchmark constants have snapshot provenance; their experimental
   material accuracy is not independently established by this pack.
2. **P1-2, discrimination:** removed the factor-two model allowance; the numerical
   envelope sums measured statistical, step, contact, adaptive, box and plateau
   bounds. Qualified rate controls must reject 0.5/2/4 mutations under the fixed
   total cap. Actual source mutations have passing literal-rule controls; the
   requested reference-prefactor perturbation is implemented in the candidate.
   A bare label-only swap is separately reported as a no-op, consistent with the
   target docstring. The independent oracle-prefactor ordering witness and
   qualified rate-budget controls have distinct, explicitly reported roles.
3. **P2-3, numerical qualification:** enlarged to three boxes and later spatially
   mixed windows; show correction fractions, counts and measured budgets rather
   than approving overlap alone. Qualify the largest-box contact-refined oracle's
   own plateau. Refine adaptive observations independently.
   Any failed cap remains visible and blocks scientific adoption. All 18 final
   rows pass the fixed caps; the bounds remain empirical, as discussed below.
4. **P2-4, missing kernel branches:** include 4/16 and 16/4, plus the N=1 floor;
   actual `min->max`, `min->i`, `min->j`, `min->mean`, and floor-removal mutations
   must fail passing literal-rule controls. Unequal orientations are separately
   simulated because their end/mid observables need not coincide.
5. **P2-5, finite-N prefactors:** add N=4 and8 to N=16. State explicitly that this
   short-chain series does not establish universal asymptotic site amplitudes.
   No theory-derived factor-two allowance or claimed general prefactor accuracy
   remains. Entangled amplitude acceptance is outside the pack's demonstrated
   coverage; the old report-only tube calculation was removed.
6. **P2-6, sampler/contact validation:** measure sampled covariances from the actual
   trajectory transition, compare to an independently constructed bond Laplacian,
   compare contact distributions to a separate bead integrator, and check the
   absorbing-sphere continuum limit. A live sampler-noise mutation tests that
   covariance qualification actually discriminates a faulty trajectory sampler.
7. **P3-7, provenance race:** hash the exact bytes before parsing/simulation,
   fingerprint the program, record dependency versions, and preserve that snapshot
   hash. The old post-simulation parameter fingerprint is gone.
8. **P3-8, concurrent publication:** UUID-specific staging plus one atomic complete
   directory publication replaces separate writes to a shared output directory.

## Limitations, evidence, compute and remaining work

The strict full Verifier command above exited **0** on 2026-10-04. Its exact
proposed scientific pytest reported `1 passed in 10369.43s (2:52:49)`, with zero
failures, errors or skips in the retained JUnit. The final output was:

```text
REPRODUCTION PASS: all 108 reported ensembles, contact hashes, covariance samples, audits and tables reproduced exactly.
PROPOSED SCIENTIFIC TEST: PASS; all required target mutations killed: True
ATOMIC PUBLICATION: /home/alon/runs/i030-met-kernel-reference/rework/verifications/verified-20261004T023214-94dde23c1cca4886820d0af07213cfb9
```

That publication contains all **318 newly simulated physical ensembles** and
12,591,104 trajectories, pooled into 108 reported rows including the floor's
explicit site-label aliases. Every scientific JSON value and rendered table
matched the frozen fixture. The two actual generation histories were separately
validated and retained rather than required to have fictitiously identical
input/program hashes. An independent post-run check of `COMPLETE.json` verified
all 331 published file hashes, with zero mismatches. The Python/NumPy/SciPy/pytest
versions were 3.9.23 / 1.26.4 / 1.13.1 / 8.4.1 on Linux x86_64, using PCG64.
`verification.json` records the exact parsed-input, program and target hashes.

The standalone command below also exited **0**, and its complete JSON exactly
matched the frozen mutation report and the full Verifier's mutation report:

```bash
taskset -c 6,16,17,19,21,23,24,26 /home/alon/anaconda3/envs/rmg_env/bin/python -B -u test/rmgpy/kmc/fixtures/i030_kernel_reference/mutation_tests.py --require-all --output /home/alon/runs/i030-met-kernel-reference/rework/mutations-final.json > >(tee -a /home/alon/runs/i030-met-kernel-reference/rework/mutations-final-stdout.log) 2> >(tee -a /home/alon/runs/i030-met-kernel-reference/rework/mutations-final-stderr.log >&2)
```

All 11 substantive target mutations were rejected after passing baseline
controls. The six requested rate mutations produced respectively 17,18,18,18,16
and4 new qualified rate rejections; neither a baseline failure nor a literal
anchor failure was substituted for this rate discrimination. Both actual
sampler mutations were rejected. The separate bare label swap survived, as
expected for the current common target multiplier. Covariance samples (75),
stationary-variance samples (75), independent Laplacian checks (18), fixed-grid
bead contacts (3), adaptive bead contacts (3), and the Brownian-sphere continuum
limit all passed their declared checks and reproduced exactly.

All 18 final reference rows qualified and all 18 candidate comparisons fell
within their measured budgets; the 13 physical-input anchors and 18 literal
kernel branches passed. The measured total log budgets are
**0.2815176766–0.4578509059**, permitting factors **1.32513942–1.58067332**.
The final largest-box correction removes about **7.9–20.6%**, below the fixed
25% cap. These results support the explicitly proposed finite-model numerical
test. They do not establish a universal site prefactor, an experimental melt
parameterization, or a rigorous error bound for unobserved limiting regimes.
In particular, N=16 mid/mid is **44.94% above** the independent reference yet
passes its empirical factor-1.5807 band. No stronger accuracy claim is justified.

The retained runtime records give the following measured compute. CPU hours
are child plus parent CPU for the reference/authoring process, and wall hours
are measured separately for each run; they exclude preliminary discarded work
and the additional Verifier-controller/audit/standalone-mutation overhead.

| Run | Physical ensembles generated | Wall hours | CPU hours | Max child RSS MiB | Max parent RSS MiB |
|---|---:|---:|---:|---:|---:|
| `measurement-04` (three seed stages) | 288 | 3.53633 | 20.03041 | 59.82422 | 88.90234 |
| `measurement-05` (five extra floor stages) | 30 | 0.03975 | 0.15208 | 52.00391 | 88.25391 |
| Strict fresh Verifier reference | 318 | 2.88028 | 19.57490 | 60.03125 | 90.36719 |

These successful reference runs used **39.75738 measured CPU hours** in total.
Each used eight workers, one numerical-library thread per worker, and the same
eight allowed CPUs listed in the command and runtime records. Per-process RSS
peaks are measured; an instantaneous aggregate RSS was not recorded. The
discarded preliminary run and serialization-failed allocation remain documented
in their logs, and their unmetered CPU is not included in this total. No other
session's jobs were stopped or altered.

No product tuning, product-code/test edits, push, or owner adoption is performed.
Owner approval and the second adversarial review remain. There is no outstanding
failure of the stated numerical criteria. The contract's required substantive
reference-prefactor swap is testable and rejected; only the separately labelled
bare target-label swap is a no-op. No unsatisfied or demonstrably incorrect
requirement in the dispatch is claimed.
