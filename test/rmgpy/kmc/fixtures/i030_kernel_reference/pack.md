# MET-KERNEL-REFERENCE: second rework, A-03

This research fixture changes no product code or existing repository test. Its
target is the enclosing checkout's `rmgpy/kmc/met.py`. Work began at
`586db906519928e1c1c5936b0abf560aff4b696b` on `i030-met-kernel-reference`.
A-03 scratch and logs are under
`/home/alon/runs/i030-met-kernel-reference/rework2/`. The manager-owned verifier
under `/home/alon/runs/verify-i030r2/` is not inspected, changed or stopped.
No prohibited dataset, catalog or pyrolysis literature is used.

## Meaning, policy chronology and adoption status

The observable is first microscopic contact of two initially stationary,
free-draining discrete Gaussian Rouse chains in a periodic cube, for end/end,
end/mid and mid/mid sites. Initial separation is uniform conditional on lying
outside the physical sink. This is an unlike-labelled pair coefficient even
for end/end. It is a finite-model benchmark, not a validation outcome.

Reproduction, numerical qualification, strict zero-model-tolerance comparison,
and binding adoption are separate. The strict comparison retains its two
unequal-length end/end failures. The owner approved **0% model tolerance on
2026-10-04** and explicitly recorded these two known deviations of the min-chain
capture rule. The binding build test requires all sixteen other rows to pass
and both named deviations to match their reported records. It emits both records
through stdout and warnings and retains `owner_adoption.json`; an additional
failure, incomplete coverage or changed deviation record fails the build.

`acceptance_policy.json` was committed in
`d86a38844f58c614c076c3cceca11bb433cdc9e3`. The self-recorded chronology places
the declaration before corrected A-03 measurements and candidate comparisons;
Git establishes the declaration's bytes in that commit, not the actual start
times of those executions. That original proposal declared **0% additional model
tolerance**, equivalently log allowance0 and factor1. The later owner ruling
approves that same value and adds the explicit known-deviation records.
There is no independently established
closure-error bound relating the analytic coil-capture radius to microscopic
first contact in these short chains. A positive allowance would be an
engineering choice, not a physics-derived bound. Zero avoids inventing one;
the owner can explicitly adopt a different policy in a later decision.
The declaration also fixes the source-only numerical method and allocation.

| Owner-approved known deviation | Measured candidate/reference | Approval date |
|---|---:|---|
| N4_N16 end/end | 0.7684941164680252 | 2026-10-04 |
| N16_N4 end/end | 0.7682664950598593 | 2026-10-04 |

Record comparisons use relative tolerance2e-12 for floating point rounding.
Neither the numerical band nor the zero model allowance is enlarged. These rows
must still be reported as scientific failures. The production run will carry
a **0.77x sensitivity case on the unequal-length end/end termination channel**
of the min-chain capture rule, covering both named orientations. This sensitivity
scenario is an owner-directed production check, not an asymptotic prefactor claim.

`decision_order.json` records UTC timestamps, the declaration commit, exact
policy/input/program hashes, campaign start, and comparison start. The Verifier
checks their order against the original declaration's Git bytes and verifies
that its numerical policy and model allowance remain unchanged by the later
approval metadata. This is **not a blinded choice**: A-02's results, including
the 44.94% N16 mid/mid excess, were already known and are explicitly recorded
as such. These checks establish agreement with the committed policy and internal
consistency of the recorded timestamps. They do not independently establish
when measurements or comparisons actually began, or ignorance of earlier results.

## Literal physical inputs and independence

`reference/parameters.json` alone supplies the physical inputs. The oracle
never imports, loads, reads or hashes `met.py`. Target access is confined to
`checks.py` outside `reference/`. Exact parameter bytes are read, hashed and
parsed before simulation; their in-memory fingerprint is retained. Program
hashes and Python/NumPy/SciPy, PCG64, platform/machine versions are retained;
the Verifier additionally records pytest and fingerprints every active
checking program, test, policy and chronology file.

The literals are benchmark conventions inherited from the frozen target
snapshot `52aa9addb`, **not independently measured melt properties**. Numerical
independence from the current target does not establish epistemic independence
of that original material convention. Separately,13 transport/geometry/SI
anchors and18 literal min/floor/two-chain rules are compared with the real API
at relative tolerance2e-12. Shared transport mutations cannot move the oracle.

At700K, `D0=2.755363626147875e-8 m2/s`,
`sigma0=6.766081442101e-10 m`, and
`Rg^2/N=0.434e-20*104.15/6 m2`. The JSON cites the frozen H/WLF evaluation
and geometry convention. Boltzmann and Avogadro constants use the
[BIPM SI defining values](https://www.bipm.org/en/publications/si-brochure).
For each finite chain, literal bond lengths enforce

```text
D_CM(N)=D0/N
Rg^2(N)=b_N^2*(N^2-1)/(6*N)
zeta=kB*T/D0
U=3*kB*T/(2*b_N^2)*sum(|r[n+1]-r[n]|^2)
```

The free-draining Gaussian dynamics follow
[Rouse (1953)](https://doi.org/10.1063/1.1699180). Bond length changes slightly
with N to match the frozen finite-chain size convention; this is not a
fixed-bond asymptotic sequence or an atomistic chain model. Coverage is equal
N=4,8,16, both4/16 and16/4, and the N=1 point-bead floor. At N=1 the physical
Rg is0, while the distinct kernel size convention is nonzero; sigma0 exceeds
twice that kernel size. Coincident N=1 site labels explicitly alias one
physical ensemble per setting. Unequal orientations are independent because
their end/mid observables need not coincide. Only the unentangled Rouse arm
is compared; no entangled amplitude is inferred.

## Actual dynamics and exact adaptive guard

Both chains' independent, unwrapped free-end modes are evolved with exact
conditional OU transitions; relative COM Brownian diffusion is D_i+D_j.
Minimum-image separation is used only for contact. End is bead0 and mid is
bead`N//2-1`. With free-end orthonormal weights,

```text
lambda_p=4*sin^2(pi*p/(2*N))
omega_p=3*D0*lambda_p/b_N^2
Var(q_p,axis)=b_N^2/(3*lambda_p)
q_p(t+dt)=exp(-omega_p*dt)*q_p(t)
          +sqrt(Var(q_p)*(1-exp(-2*omega_p*dt)))*Z
```

The proposed observation interval obeys the local Brownian gap scale, a
maximum step and a declared minimum contact step. For the current gapg,
the exact conditional mean displacement is
`sum_p w_p*q_p*(exp(-omega_p*dt)-1)`. `guard_observation_step()` evaluates that
quantity at the actual interval. Every interval above its floor that violates
`norm(mean_shift)<=g/safety` is halved and **re-evaluated**, until the inequality
holds or the interval reaches the floor. It assumes neither linear OU decay
nor monotonicity of the norm of a sum of modes. The final time truncation is
included before guarding; an end-of-run remainder below dt_min is permitted.
The guard's inequality is not guaranteed at the declared minimum itself.

`guardTest.py` drives equilibrium states for every case and pair through the
actual helper at base/contact/adaptive settings. It asserts the exact bound
or a floor interval. A separate regression reproduces failure of the former
one-shot linear shrink. Other tests verify the64-form identity and the
whole/half-window covariance used below.

Contact is detected at the numerical surface
`sigma0+0.5825971579390107*sqrt(2*D_local*dt_min)`;
`D_local=2*D0`. Initial exclusion uses physical sigma0. The sign follows
[Gobet–Menozzi](https://arxiv.org/abs/0706.4042); their Euler theorem does not
prove a bound for exact OU transitions with adaptive observation. Minimum
contact steps and the far-step maximum/safety are separately refined. Matching
a covariance alone is not first-passage validation.

## Complete sampler qualification inside the proposed test

The directly runnable `proposed_adoption_test.py` and the identical body below call
`audit_sampler.run_audits(cfg)` **before** the full trajectory campaign. All four
audits execute: stationary variance of actual OU innovations, independent bond
Laplacian geometry/diffusion, a Brownian-sphere continuum contact check, and
adaptive contact statistics against independently evolved Euler beads.
The reference campaign additionally measures actual trajectory covariance at
the declared lags and fixed-grid bead/contact statistics. Counts, SEs, seeds,
step settings and contact fingerprints are retained, not just theoretical
mode identities.

The small-N Euler contact comparisons use independent Gaussian-bond initial
states and force/noise evolution,5 combined SE and a separately reported
Euler weak-error scale. N=1 continuum contact uses the exact infinite-sphere
hitting probability, with an explicit upper bound on the omitted cube tail.
These checks do not establish all large-N/adaptive continuum remainders.

The new mutation doubles OU innovation noise **only for modes belonging to
length16 chains**, including both unequal orientations; length4/8 modes remain
unchanged. It executes the exact adoption-test body in real pytest. A valid
kill is one stationary-variance assertion failure before any campaign output.
Its positive control is the complete unmodified audit phase. If the full
unmodified scientific comparison fails, that failure is never counted as a
mutation kill. The original doubled-noise and doubled-coordinate controls
remain separate. Actual mutant JUnit, generated test and terminal log are
retained in the Verifier publication.

## Transient bins: restored, report-only

Dedicated L/scale=6 ensembles restore exactly the original four reduced-time
bins `[0,.02]`, `[.02,.05]`, `[.05,.1]`, `[.1,.2]`. Times use each case's
`scale^2/D_unit`; unequal-case transient edges are not rescaled to a late
spatial-mixing window. The input JSON declares counts, seeds and contact steps
before measurement. Integer survivors are pooled across independent stages.

For counts n1,n2, bin widthDelta and cube volumeV,
`k_bin=V/Delta*ln(n1/n2)` and
`SE=V/Delta*sqrt((n1-n2)/(n1*n2))`. Survival has binomial SE. These are
**finite-box conditional hazards**, using the original cube-volume convention
despite conditioning initial positions outside the sink. No late Green
correction or plateau qualification is applied to them. Zero-event bins have
rate/SE0 and an undefined candidate ratio, reported honestly.

Every candidate diagnostic is marked REPORT-ONLY. Its original agreement
criterion is20% of the reference plus4SE; it does not change qualification,
the long-time comparison, any mutation kill, or the approved model tolerance.

## Numerical uncertainty and owner-approved model tolerance

The long-time observable is the largest-box contact-refined setting`box2`.
Conditional late counts give

```text
k_box=L^3/Delta*ln(n_start/n_end)
SE_box=L^3/Delta*sqrt((n_start-n_end)/(n_start*n_end))
k_inf=k_box/(1+G*k_box)
SE_inf=SE_box/(1+G*k_box)^2
G=2.837297479/(4*pi*D_relative*L)
```

The Ewald constant is independently recomputed. This is a leading far-field
Green approximation, not an exact finite-polymer reduction. Three boxes20,28,36
and windows scaled by relative COM diffusion start after at least4 spatial
mixing times. The final largest-box correction is capped at25%; the primary
oracle's own two late halves supply the plateau check.

Each individual audit keeps the A-02 qualification ceiling and upper contrast
`U(A,B)=abs(ln(kA/kB))+4*hypot(sA/kA,sB/kB)`. The limits are statistics.16,
step/contact/adaptive.12 each, box.16, plateau.20, combined numerical.55,
at least4 mixing times, and at most25% primary correction. These are accuracy
targets, not fitted model allowances. Failed qualification stays visible and
blocks the binding adoption test. The point-bead primary also checks the
analytic Smoluchowski sphere coefficient.

The **acceptance envelope no longer sums marginal sampling margins repeatedly**.
Let x be the log rates at base,step,contact,adaptive,box,box2 and the two box2
halves. Define the observed audit sum

```text
F(x)=|x_base-x_step|+|x_step-x_contact|+|x_contact-x_adaptive|
     +max(|x_contact-x_box|,|x_box-x_box2|)+|x_half1-x_half2|
    =max_j a_j*x                     (64 signed linear forms)
z_joint=normal.isf(2*normal.sf(4)/64)
B_audit=max_j(a_j*x+z_joint*sqrt(a_j*Sigma*a_j'))
B_numerical=4*SE_box2/k_box2+B_audit
```

Sigma uses the delta-method covariance of conditional log hazards. Distinct
ensembles are independent. Disjoint half-window conditional-hazard increments
have zero first-order covariance; the whole raw hazard is their
duration-weighted mean. Its nonlinear Green/log derivatives preserve covariance
between the primary whole-window estimate and both halves. Ignoring that shared
primary estimate or adding a separate4SE margin for every reuse is not the
chosen statistical accounting. Bonferroni bounds each of the64 linear forms
with total one-sided family tail allowance`2*normal.sf(4)` for one case/pair;
the primary4SE term is additional and conservative. This is a **normal/delta
approximation**, not a finite-sample theorem or simultaneous statement over
all18 rows. Individual marginal cap checks remain separately enforced.

`uncertainty_probe.json` re-examines fixed A-02 means under this predeclared
method; it is never an A-03 measurement or acceptance input. For the historical
N16 mid/mid row, the old0.4578509 sum becomes0.3168077: observed contrasts0.0747255,
joint audit sampling0.1936973, primary statistics0.0483849. This changes accounting,
**not the SEs or replica counts**. Sampling still dominates; more replicas would
shrink it but could not establish the missing continuum/volume remainder.
The already-known44.94% excess has absolute log deviation0.3711299, exceeding
that historical numerical envelope by0.0543223 with proposed model allowance0.
These counterfactual old-data values are separate from the fresh tables below.

For every fresh row, the second band is the single, explicit policy value
`B_model=0` (owner-approved on 2026-10-04). The strict scientific criterion is
`abs(ln(k_candidate/k_reference))<=B_numerical+B_model`.
Reports show both allowances, numerical consumption, residual model allowance
required, and remaining excess. Consumption means deterministic bookkeeping,
using numerical allowance first; it does not identify a physical model-error
component. No tolerance is changed after comparison.

## Fresh reference measurements

<!-- BEGIN REPRODUCED NUMBERS -->

All +/- values are one Monte Carlo standard error. Individual cap checks and the joint numerical band are empirical finite-setting diagnostics, not rigorous continuum-error bounds.

| Case | Pair | Reference / (D_unit length_unit) | Reference (m3 mol-1 s-1) | stat | step cap | contact cap | adaptive cap | box cap | plateau cap | Numerical log | Factor | Qualified |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| N4 | end/end | 51.7130491 +/- 0.4783846 | 117760285 +/- 1089371 | 0.03700 | 0.05801 | 0.04777 | 0.04401 | 0.06748 | 0.07915 | 0.23470 | 1.26453 | YES |
| N4 | end/mid | 48.767603 +/- 0.4702778 | 111052954 +/- 1070910 | 0.03857 | 0.05274 | 0.05163 | 0.04590 | 0.06324 | 0.08004 | 0.23926 | 1.27031 | YES |
| N4 | mid/mid | 45.9238113 +/- 0.4614147 | 104577108 +/- 1050728 | 0.04019 | 0.05689 | 0.06562 | 0.04502 | 0.07197 | 0.08657 | 0.26989 | 1.30982 | YES |
| N8 | end/end | 48.7477555 +/- 0.4705174 | 78494337.9 +/- 757633.8 | 0.03861 | 0.05987 | 0.04603 | 0.04661 | 0.06083 | 0.08382 | 0.22229 | 1.24893 | YES |
| N8 | end/mid | 45.1002992 +/- 0.4586319 | 72621151.2 +/- 738495.7 | 0.04068 | 0.07798 | 0.06391 | 0.04647 | 0.07507 | 0.11344 | 0.30761 | 1.36016 | YES |
| N8 | mid/mid | 39.8362115 +/- 0.4403997 | 64144841.4 +/- 709137.9 | 0.04422 | 0.08503 | 0.04885 | 0.06424 | 0.06197 | 0.13127 | 0.30948 | 1.36272 | YES |
| N16 | end/end | 48.1768072 +/- 0.4683963 | 54853800.6 +/- 533312.9 | 0.03889 | 0.07304 | 0.04759 | 0.07058 | 0.08033 | 0.09897 | 0.30881 | 1.36180 | YES |
| N16 | end/mid | 41.9693573 +/- 0.4483644 | 47786038.4 +/- 510504.8 | 0.04273 | 0.06318 | 0.05646 | 0.06959 | 0.07243 | 0.09197 | 0.29820 | 1.34743 | YES |
| N16 | mid/mid | 35.6021739 +/- 0.4232569 | 40536404.6 +/- 481917.6 | 0.04755 | 0.10688 | 0.08066 | 0.06867 | 0.09372 | 0.11177 | 0.40225 | 1.49518 | YES |
| N4_N16 | end/end | 40.879853 +/- 0.3180851 | 93091071.7 +/- 724339.2 | 0.03112 | 0.04524 | 0.03517 | 0.03419 | 0.05151 | 0.07234 | 0.19194 | 1.21159 | YES |
| N4_N16 | end/mid | 36.2421646 +/- 0.3088361 | 82530187.6 +/- 703277.6 | 0.03409 | 0.05748 | 0.03868 | 0.04464 | 0.05412 | 0.09874 | 0.23317 | 1.26260 | YES |
| N4_N16 | mid/mid | 33.7092263 +/- 0.3028701 | 76762213.4 +/- 689691.8 | 0.03594 | 0.05076 | 0.04041 | 0.04086 | 0.06674 | 0.07386 | 0.20910 | 1.23256 | YES |
| N16_N4 | end/end | 40.8919649 +/- 0.3178811 | 93118652.7 +/- 723874.6 | 0.03109 | 0.04850 | 0.03673 | 0.04169 | 0.04853 | 0.07626 | 0.20218 | 1.22406 | YES |
| N16_N4 | end/mid | 38.6342476 +/- 0.3137709 | 87977408 +/- 714515 | 0.03249 | 0.05529 | 0.03662 | 0.04590 | 0.06009 | 0.09222 | 0.22309 | 1.24994 | YES |
| N16_N4 | mid/mid | 33.9333407 +/- 0.3032155 | 77272563.9 +/- 690478.5 | 0.03574 | 0.04830 | 0.04080 | 0.04133 | 0.06161 | 0.07586 | 0.20594 | 1.22868 | YES |
| floor | end/end | 25.2228164 +/- 0.2271163 | 283178726 +/- 2549854 | 0.03602 | 0.06040 | 0.04209 | 0.04919 | 0.06587 | 0.09298 | 0.24593 | 1.27881 | YES |
| floor | end/mid | 25.2228164 +/- 0.2271163 | 283178726 +/- 2549854 | 0.03602 | 0.06040 | 0.04209 | 0.04919 | 0.06587 | 0.09298 | 0.24593 | 1.27881 | YES |
| floor | mid/mid | 25.2228164 +/- 0.2271163 | 283178726 +/- 2549854 | 0.03602 | 0.06040 | 0.04209 | 0.04919 | 0.06587 | 0.09298 | 0.24593 | 1.27881 | YES |

| Case | Pair | Run | Events | k_box +/- SE | k_infinite +/- SE | Correction removed | Late half 1 | Late half 2 | Mixing times at start |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| N4 | end/end | base | 5223 | 71.75227 +/- 0.993684 | 51.068701 +/- 0.503369 | 0.28826 | 50.5914 +/- 0.6995 | 51.54241 +/- 0.72361 | 4.73741 |
| N4 | end/end | step | 10555 | 72.75642 +/- 0.708803 | 51.575328 +/- 0.356178 | 0.29112 | 51.801 +/- 0.49427 | 51.34884 +/- 0.5131 | 4.73741 |
| N4 | end/end | contact | 10447 | 71.887273 +/- 0.703931 | 51.137052 +/- 0.356203 | 0.28865 | 51.1707 +/- 0.49459 | 51.10338 +/- 0.51276 | 4.73741 |
| N4 | end/end | box | 8563 | 63.952409 +/- 0.691246 | 50.842714 +/- 0.436894 | 0.20499 | 50.9434 +/- 0.61259 | 50.74192 +/- 0.62311 | 4.83409 |
| N4 | end/end | box2 | 8204 | 61.722388 +/- 0.681494 | 51.713049 +/- 0.478385 | 0.16217 | 51.84596 +/- 0.67332 | 51.58 +/- 0.67975 | 4.87388 |
| N4 | end/end | adaptive | 10512 | 72.364537 +/- 0.706419 | 51.378096 +/- 0.356096 | 0.29001 | 51.67554 +/- 0.49415 | 51.07924 +/- 0.51303 | 4.73741 |
| N4 | end/mid | base | 4931 | 66.330368 +/- 0.945286 | 48.260976 +/- 0.500415 | 0.27242 | 49.1923 +/- 0.69526 | 47.316 +/- 0.72041 | 4.73741 |
| N4 | end/mid | step | 9875 | 66.5099 +/- 0.669789 | 48.355947 +/- 0.35405 | 0.27295 | 48.45111 +/- 0.49227 | 48.26064 +/- 0.50904 | 4.73741 |
| N4 | end/mid | contact | 9751 | 65.60241 +/- 0.664824 | 47.874455 +/- 0.354058 | 0.27023 | 47.54842 +/- 0.49254 | 48.19885 +/- 0.50858 | 4.73741 |
| N4 | end/mid | box | 8213 | 60.59454 +/- 0.668747 | 48.697319 +/- 0.431922 | 0.19634 | 48.5876 +/- 0.60555 | 48.80692 +/- 0.61604 | 4.83409 |
| N4 | end/mid | box2 | 7717 | 57.572135 +/- 0.655415 | 48.767603 +/- 0.470278 | 0.15293 | 48.69704 +/- 0.66156 | 48.83813 +/- 0.66856 | 4.87388 |
| N4 | end/mid | adaptive | 9707 | 65.244884 +/- 0.662692 | 47.68377 +/- 0.353965 | 0.26916 | 47.56312 +/- 0.49239 | 47.80419 +/- 0.50858 | 4.73741 |
| N4 | mid/mid | base | 4701 | 61.82586 +/- 0.902302 | 45.831435 +/- 0.495836 | 0.25870 | 44.89797 +/- 0.69024 | 46.75181 +/- 0.71124 | 4.73741 |
| N4 | mid/mid | step | 9340 | 61.509465 +/- 0.636858 | 45.657338 +/- 0.350897 | 0.25772 | 45.14305 +/- 0.48853 | 46.16764 +/- 0.50356 | 4.73741 |
| N4 | mid/mid | contact | 9106 | 59.743325 +/- 0.626446 | 44.676972 +/- 0.350326 | 0.25218 | 44.72127 +/- 0.48799 | 44.63265 +/- 0.50279 | 4.73741 |
| N4 | mid/mid | box | 7712 | 56.07125 +/- 0.638594 | 45.73242 +/- 0.424808 | 0.18439 | 45.5088 +/- 0.59565 | 45.95555 +/- 0.6058 | 4.83409 |
| N4 | mid/mid | box2 | 7259 | 53.650101 +/- 0.629733 | 45.923811 +/- 0.461415 | 0.14401 | 46.06577 +/- 0.65004 | 45.78171 +/- 0.65504 | 4.87388 |
| N4 | mid/mid | adaptive | 9101 | 59.796462 +/- 0.627176 | 44.706681 +/- 0.350577 | 0.25235 | 44.96706 +/- 0.48836 | 44.44527 +/- 0.50322 | 4.73741 |
| N8 | end/end | base | 5003 | 67.793185 +/- 0.959187 | 49.030738 +/- 0.501728 | 0.27676 | 48.55736 +/- 0.69774 | 49.50064 +/- 0.72078 | 4.73741 |
| N8 | end/end | step | 9902 | 66.902131 +/- 0.672825 | 48.562947 +/- 0.354513 | 0.27412 | 48.52336 +/- 0.49292 | 48.60251 +/- 0.50963 | 4.73741 |
| N8 | end/end | contact | 9990 | 67.352651 +/- 0.674373 | 48.79989 +/- 0.354021 | 0.27546 | 48.74213 +/- 0.49219 | 48.85759 +/- 0.50896 | 4.73741 |
| N8 | end/end | box | 8307 | 61.325196 +/- 0.672974 | 49.168111 +/- 0.432601 | 0.19824 | 48.86284 +/- 0.60611 | 49.47245 +/- 0.61734 | 4.83409 |
| N8 | end/end | box2 | 7704 | 57.544476 +/- 0.655652 | 48.747755 +/- 0.470517 | 0.15287 | 48.90842 +/- 0.66263 | 48.5869 +/- 0.66819 | 4.87388 |
| N8 | end/end | adaptive | 10027 | 67.881089 +/- 0.678417 | 49.076702 +/- 0.354609 | 0.27702 | 48.78861 +/- 0.49305 | 49.3635 +/- 0.50963 | 4.73741 |
| N8 | end/mid | base | 4490 | 58.470768 +/- 0.873098 | 43.961479 +/- 0.493549 | 0.24815 | 43.68676 +/- 0.68759 | 44.23507 +/- 0.70804 | 4.73741 |
| N8 | end/mid | step | 9218 | 60.323004 +/- 0.628678 | 45.000353 +/- 0.34986 | 0.25401 | 44.7816 +/- 0.48723 | 45.21839 +/- 0.5021 | 4.73741 |
| N8 | end/mid | contact | 9015 | 58.768524 +/- 0.619315 | 44.129584 +/- 0.349206 | 0.24909 | 43.65059 +/- 0.48641 | 44.60515 +/- 0.50094 | 4.73741 |
| N8 | end/mid | box | 7428 | 53.840775 +/- 0.624796 | 44.237691 +/- 0.421794 | 0.17836 | 44.67093 +/- 0.59315 | 43.8026 +/- 0.59989 | 4.83409 |
| N8 | end/mid | box2 | 7129 | 52.529562 +/- 0.622176 | 45.100299 +/- 0.458632 | 0.14143 | 44.37492 +/- 0.64296 | 45.82185 +/- 0.65406 | 4.87388 |
| N8 | end/mid | adaptive | 9040 | 58.905349 +/- 0.6199 | 44.206689 +/- 0.349131 | 0.24953 | 44.50655 +/- 0.48648 | 43.90547 +/- 0.50102 | 4.73741 |
| N8 | mid/mid | base | 3963 | 49.785654 +/- 0.791173 | 38.864033 +/- 0.482124 | 0.21937 | 38.68285 +/- 0.67299 | 39.04474 +/- 0.69046 | 4.73741 |
| N8 | mid/mid | step | 8181 | 51.392427 +/- 0.568443 | 39.836283 +/- 0.341544 | 0.22486 | 39.96768 +/- 0.4769 | 39.70464 +/- 0.4891 | 4.73741 |
| N8 | mid/mid | contact | 8175 | 51.415444 +/- 0.568907 | 39.850111 +/- 0.341753 | 0.22494 | 39.86858 +/- 0.47708 | 39.83164 +/- 0.48947 | 4.73741 |
| N8 | mid/mid | box | 6698 | 47.552592 +/- 0.5811 | 39.902284 +/- 0.409164 | 0.16088 | 40.16669 +/- 0.57571 | 39.6372 +/- 0.58159 | 4.83409 |
| N8 | mid/mid | box2 | 6266 | 45.523065 +/- 0.575114 | 39.836211 +/- 0.4404 | 0.12492 | 40.68551 +/- 0.62416 | 38.98172 +/- 0.62141 | 4.87388 |
| N8 | mid/mid | adaptive | 8048 | 50.404335 +/- 0.562092 | 39.240019 +/- 0.340667 | 0.22150 | 39.31296 +/- 0.47575 | 39.16701 +/- 0.48776 | 4.73741 |
| N16 | end/end | base | 4920 | 65.694344 +/- 0.937255 | 47.923396 +/- 0.498766 | 0.27051 | 47.26951 +/- 0.69396 | 48.57073 +/- 0.71604 | 4.73741 |
| N16 | end/end | step | 9555 | 63.779973 +/- 0.652925 | 46.896558 +/- 0.353002 | 0.26471 | 46.68099 +/- 0.49124 | 47.11141 +/- 0.50697 | 4.73741 |
| N16 | end/end | contact | 9504 | 63.355208 +/- 0.650309 | 46.666504 +/- 0.35283 | 0.26341 | 46.32202 +/- 0.49106 | 47.00918 +/- 0.50659 | 4.73741 |
| N16 | end/end | box | 8180 | 59.920334 +/- 0.662636 | 48.260919 +/- 0.429851 | 0.19458 | 48.91503 +/- 0.60408 | 47.6025 +/- 0.61176 | 4.83409 |
| N16 | end/end | box2 | 7625 | 56.750553 +/- 0.649946 | 48.176807 +/- 0.468396 | 0.15108 | 47.66527 +/- 0.65749 | 48.68642 +/- 0.66721 | 4.87388 |
| N16 | end/end | adaptive | 9833 | 65.857506 +/- 0.664624 | 48.010166 +/- 0.353209 | 0.27100 | 48.05763 +/- 0.49121 | 47.96267 +/- 0.50771 | 4.73741 |
| N16 | end/mid | base | 4297 | 54.997911 +/- 0.839426 | 41.968961 +/- 0.488817 | 0.23690 | 42.2731 +/- 0.68197 | 41.66345 +/- 0.70063 | 4.73741 |
| N16 | end/mid | step | 8544 | 54.56544 +/- 0.590612 | 41.716653 +/- 0.345212 | 0.23547 | 41.16551 +/- 0.48112 | 42.26335 +/- 0.49492 | 4.73741 |
| N16 | end/mid | contact | 8641 | 55.274162 +/- 0.594924 | 42.129637 +/- 0.345615 | 0.23781 | 42.54584 +/- 0.4822 | 41.71085 +/- 0.49539 | 4.73741 |
| N16 | end/mid | box | 6943 | 49.680638 +/- 0.596303 | 41.389974 +/- 0.413888 | 0.16688 | 41.39359 +/- 0.58135 | 41.38636 +/- 0.58928 | 4.83409 |
| N16 | end/mid | box2 | 6608 | 48.330179 +/- 0.59457 | 41.969357 +/- 0.448364 | 0.13161 | 41.83273 +/- 0.63088 | 42.10585 +/- 0.63726 | 4.87388 |
| N16 | end/mid | adaptive | 8379 | 53.66092 +/- 0.586503 | 41.185891 +/- 0.345502 | 0.23248 | 41.25471 +/- 0.48207 | 41.117 +/- 0.4951 | 4.73741 |
| N16 | mid/mid | base | 3497 | 42.681083 +/- 0.721971 | 34.394746 +/- 0.468849 | 0.19415 | 34.91934 +/- 0.65742 | 33.86627 +/- 0.66875 | 4.73741 |
| N16 | mid/mid | step | 7336 | 44.91075 +/- 0.524526 | 35.828158 +/- 0.333822 | 0.20224 | 35.56775 +/- 0.4663 | 36.08761 +/- 0.47772 | 4.73741 |
| N16 | mid/mid | contact | 7123 | 43.40065 +/- 0.5144 | 34.86051 +/- 0.331876 | 0.19677 | 35.43434 +/- 0.46531 | 34.28202 +/- 0.47343 | 4.73741 |
| N16 | mid/mid | box | 6053 | 42.281524 +/- 0.543506 | 36.123426 +/- 0.396717 | 0.14565 | 35.93487 +/- 0.55706 | 36.31165 +/- 0.56497 | 4.83409 |
| N16 | mid/mid | box2 | 5584 | 40.076514 +/- 0.536329 | 35.602174 +/- 0.423257 | 0.11164 | 35.30524 +/- 0.59494 | 35.89848 +/- 0.60216 | 4.87388 |
| N16 | mid/mid | adaptive | 7217 | 44.225439 +/- 0.520757 | 35.390658 +/- 0.333479 | 0.19977 | 35.67059 +/- 0.46687 | 35.10961 +/- 0.47636 | 4.73741 |
| N4_N16 | end/end | base | 6507 | 63.339725 +/- 0.786555 | 40.291259 +/- 0.318272 | 0.36389 | 40.33685 +/- 0.4384 | 40.24561 +/- 0.46157 | 4.73741 |
| N4_N16 | end/end | step | 12907 | 62.703234 +/- 0.552848 | 40.032764 +/- 0.22535 | 0.36155 | 40.12805 +/- 0.31037 | 39.93722 +/- 0.32688 | 4.73741 |
| N4_N16 | end/end | contact | 12888 | 62.382141 +/- 0.550413 | 39.901639 +/- 0.22519 | 0.36037 | 39.73338 +/- 0.31082 | 40.0691 +/- 0.32578 | 4.73741 |
| N4_N16 | end/end | box | 10975 | 54.859928 +/- 0.523865 | 40.519886 +/- 0.285789 | 0.26139 | 41.00292 +/- 0.3992 | 40.03274 +/- 0.4092 | 4.83409 |
| N4_N16 | end/end | box2 | 10439 | 51.428467 +/- 0.503422 | 40.879853 +/- 0.318085 | 0.20511 | 40.67326 +/- 0.44638 | 41.0859 +/- 0.45324 | 4.87388 |
| N4_N16 | end/end | adaptive | 12893 | 62.605469 +/- 0.552283 | 39.992891 +/- 0.225373 | 0.36119 | 40.28915 +/- 0.30993 | 39.69413 +/- 0.32754 | 4.73741 |
| N4_N16 | end/mid | base | 5880 | 53.227576 +/- 0.694981 | 35.947104 +/- 0.316977 | 0.32465 | 36.26738 +/- 0.43788 | 35.62406 +/- 0.45873 | 4.73741 |
| N4_N16 | end/mid | step | 11516 | 52.136933 +/- 0.486405 | 35.446337 +/- 0.224828 | 0.32013 | 35.24885 +/- 0.31156 | 35.6428 +/- 0.32408 | 4.73741 |
| N4_N16 | end/mid | contact | 11546 | 52.356259 +/- 0.487821 | 35.547578 +/- 0.224876 | 0.32104 | 35.48198 +/- 0.3114 | 35.61306 +/- 0.32447 | 4.73741 |
| N4_N16 | end/mid | box | 9713 | 46.81501 +/- 0.475149 | 35.95613 +/- 0.280289 | 0.23195 | 36.02484 +/- 0.39235 | 35.88734 +/- 0.4004 | 4.83409 |
| N4_N16 | end/mid | box2 | 9220 | 44.297323 +/- 0.461376 | 36.242165 +/- 0.308836 | 0.18184 | 36.79349 +/- 0.43522 | 35.68709 +/- 0.43831 | 4.87388 |
| N4_N16 | end/mid | adaptive | 11724 | 53.061538 +/- 0.490641 | 35.871298 +/- 0.224232 | 0.32397 | 35.94708 +/- 0.31019 | 35.79537 +/- 0.32394 | 4.73741 |
| N4_N16 | mid/mid | base | 5486 | 48.271299 +/- 0.652368 | 33.61611 +/- 0.31638 | 0.30360 | 34.13524 +/- 0.43782 | 33.08989 +/- 0.45723 | 4.73741 |
| N4_N16 | mid/mid | step | 10936 | 47.954004 +/- 0.459009 | 33.461922 +/- 0.223498 | 0.30221 | 33.5786 +/- 0.30982 | 33.34489 +/- 0.32227 | 4.73741 |
| N4_N16 | mid/mid | contact | 10853 | 47.780417 +/- 0.45909 | 33.377308 +/- 0.224027 | 0.30144 | 33.59687 +/- 0.31046 | 33.1565 +/- 0.32318 | 4.73741 |
| N4_N16 | mid/mid | box | 8913 | 42.113879 +/- 0.446181 | 33.116815 +/- 0.275904 | 0.21364 | 33.35273 +/- 0.38686 | 32.87999 +/- 0.39353 | 4.83409 |
| N4_N16 | mid/mid | box2 | 8553 | 40.571202 +/- 0.438727 | 33.709226 +/- 0.30287 | 0.16913 | 33.67583 +/- 0.42585 | 33.74261 +/- 0.43078 | 4.87388 |
| N4_N16 | mid/mid | adaptive | 10810 | 47.587352 +/- 0.45814 | 33.282981 +/- 0.224109 | 0.30059 | 33.15057 +/- 0.31099 | 33.41494 +/- 0.32269 | 4.73741 |
| N16_N4 | end/end | base | 6536 | 63.428675 +/- 0.785914 | 40.327234 +/- 0.317688 | 0.36421 | 40.93812 +/- 0.43559 | 39.70556 +/- 0.46331 | 4.73741 |
| N16_N4 | end/end | step | 12879 | 62.466915 +/- 0.551356 | 39.936305 +/- 0.225355 | 0.36068 | 39.96792 +/- 0.31057 | 39.90466 +/- 0.32666 | 4.73741 |
| N16_N4 | end/end | contact | 12956 | 62.948291 +/- 0.553965 | 40.132512 +/- 0.225168 | 0.36245 | 40.33006 +/- 0.30984 | 39.93385 +/- 0.32698 | 4.73741 |
| N16_N4 | end/end | box | 11021 | 55.049273 +/- 0.524577 | 40.623088 +/- 0.285661 | 0.26206 | 40.72024 +/- 0.39908 | 40.52577 +/- 0.40887 | 4.83409 |
| N16_N4 | end/end | box2 | 10457 | 51.447637 +/- 0.503176 | 40.891965 +/- 0.317881 | 0.20517 | 40.60367 +/- 0.44598 | 41.17922 +/- 0.45304 | 4.87388 |
| N16_N4 | end/end | adaptive | 12796 | 61.994589 +/- 0.548945 | 39.742724 +/- 0.225599 | 0.35893 | 39.86598 +/- 0.31075 | 39.61904 +/- 0.32724 | 4.73741 |
| N16_N4 | end/mid | base | 6112 | 57.072798 +/- 0.73104 | 37.660698 +/- 0.318317 | 0.34013 | 36.85143 +/- 0.44194 | 38.45242 +/- 0.45727 | 4.73741 |
| N16_N4 | end/mid | step | 12369 | 58.306237 +/- 0.525022 | 38.193855 +/- 0.225286 | 0.34494 | 37.79589 +/- 0.3119 | 38.5875 +/- 0.32484 | 4.73741 |
| N16_N4 | end/mid | contact | 12453 | 58.604338 +/- 0.525931 | 38.321545 +/- 0.224882 | 0.34610 | 38.3771 +/- 0.31038 | 38.2659 +/- 0.32555 | 4.73741 |
| N16_N4 | end/mid | box | 10678 | 52.479794 +/- 0.508042 | 39.206535 +/- 0.283552 | 0.25292 | 39.21471 +/- 0.39638 | 39.19836 +/- 0.40558 | 4.83409 |
| N16_N4 | end/mid | box2 | 9855 | 47.924103 +/- 0.482809 | 38.634248 +/- 0.313771 | 0.19385 | 38.10617 +/- 0.4399 | 39.15888 +/- 0.44742 | 4.87388 |
| N16_N4 | end/mid | adaptive | 12245 | 57.502655 +/- 0.520381 | 37.847392 +/- 0.225433 | 0.34181 | 38.38603 +/- 0.31029 | 37.30067 +/- 0.3275 | 4.73741 |
| N16_N4 | mid/mid | base | 5389 | 47.3764 +/- 0.645987 | 33.179652 +/- 0.316842 | 0.29966 | 32.70603 +/- 0.4402 | 33.64756 +/- 0.45537 | 4.73741 |
| N16_N4 | mid/mid | step | 10823 | 47.482285 +/- 0.456852 | 33.231551 +/- 0.223776 | 0.30013 | 33.5532 +/- 0.31004 | 32.90721 +/- 0.32294 | 4.73741 |
| N16_N4 | mid/mid | contact | 10850 | 47.669927 +/- 0.45809 | 33.323353 +/- 0.223851 | 0.30096 | 33.33286 +/- 0.31047 | 33.31384 +/- 0.32257 | 4.73741 |
| N16_N4 | mid/mid | box | 9067 | 42.727373 +/- 0.448823 | 33.495003 +/- 0.275818 | 0.21608 | 33.33724 +/- 0.38622 | 33.65236 +/- 0.39383 | 4.83409 |
| N16_N4 | mid/mid | box2 | 8624 | 40.896286 +/- 0.440419 | 33.933341 +/- 0.303216 | 0.17026 | 33.85901 +/- 0.42621 | 34.00761 +/- 0.43138 | 4.87388 |
| N16_N4 | mid/mid | adaptive | 10888 | 47.901638 +/- 0.459517 | 33.436417 +/- 0.223893 | 0.30198 | 33.03326 +/- 0.31095 | 33.83541 +/- 0.32194 | 4.73741 |
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

| Case | Pair | Observed audit sum | Joint audit sampling margin | Primary statistics | Numerical band | Old sum of marginal bounds |
|---|---|---:|---:|---:|---:|---:|
| N4 | end/end | 0.04522461 | 0.1524735 | 0.03700301 | 0.2347011 | 0.3334227 |
| N4 | end/mid | 0.03589912 | 0.1647893 | 0.03857297 | 0.2392614 | 0.3321246 |
| N4 | mid/mid | 0.05571143 | 0.1739903 | 0.04018958 | 0.2698913 | 0.3662679 |
| N8 | end/end | 0.03529171 | 0.1483878 | 0.03860833 | 0.2222879 | 0.3357654 |
| N8 | end/mid | 0.09604049 | 0.170888 | 0.04067662 | 0.3076051 | 0.4175512 |
| N8 | mid/mid | 0.08492047 | 0.1803408 | 0.04422104 | 0.3094823 | 0.4355903 |
| N16 | end/end | 0.1097558 | 0.1601644 | 0.03888977 | 0.30881 | 0.4094049 |
| N16 | end/mid | 0.06275717 | 0.1927081 | 0.04273255 | 0.2981978 | 0.3963658 |
| N16 | mid/mid | 0.1355534 | 0.2191401 | 0.04755405 | 0.4022476 | 0.5092593 |
| N4_N16 | end/end | 0.03747109 | 0.1233405 | 0.03112389 | 0.1919355 | 0.2695553 |
| N4_N16 | end/mid | 0.0679057 | 0.1311796 | 0.03408583 | 0.2331711 | 0.3277466 |
| N4_N16 | mid/mid | 0.02967079 | 0.1434874 | 0.03593913 | 0.2090973 | 0.3085675 |
| N16_N4 | end/end | 0.05062734 | 0.1204542 | 0.03109472 | 0.2021762 | 0.2827984 |
| N16_N4 | end/mid | 0.07992745 | 0.1106791 | 0.0324863 | 0.2230929 | 0.3225927 |
| N16_N4 | mid/mid | 0.02508965 | 0.1451067 | 0.03574249 | 0.2059388 | 0.303645 |
| floor | end/end | 0.0535733 | 0.1563381 | 0.0360176 | 0.245929 | 0.3465407 |
| floor | end/mid | 0.0535733 | 0.1563381 | 0.0360176 | 0.245929 | 0.3465407 |
| floor | mid/mid | 0.0535733 | 0.1563381 | 0.0360176 | 0.245929 | 0.3465407 |

Dedicated transient ensembles: finite-box conditional hazards, REPORT-ONLY. Times are in each case time unit; no Green correction or acceptance is applied.

| Case | Pair | t1 | t2 | S(t2) +/- SE | k_bin reduced +/- SE | k_bin (m3 mol-1 s-1) +/- SE |
|---|---|---:|---:|---:|---:|---:|
| N4 | end/end | 0 | 0.02 | 0.94710286 +/- 0.000713886 | 586.95376 +/- 8.14059 | 1.3366035e+09 +/- 1.85376e+07 |
| N4 | end/end | 0.02 | 0.05 | 0.90698242 +/- 0.000926395 | 311.6494 +/- 4.96286 | 7.0968398e+08 +/- 1.13014e+07 |
| N4 | end/end | 0.05 | 0.1 | 0.85635376 +/- 0.00111863 | 248.13867 +/- 3.51779 | 5.6505817e+08 +/- 8.01067e+06 |
| N4 | end/end | 0.1 | 0.2 | 0.77462769 +/- 0.00133263 | 216.65028 +/- 2.41811 | 4.933532e+08 +/- 5.50648e+06 |
| N4 | end/mid | 0 | 0.02 | 0.94799805 +/- 0.000708154 | 576.75064 +/- 8.0676 | 1.3133691e+09 +/- 1.83714e+07 |
| N4 | end/mid | 0.02 | 0.05 | 0.9092509 +/- 0.000916173 | 300.46589 +/- 4.86879 | 6.8421703e+08 +/- 1.10871e+07 |
| N4 | end/mid | 0.05 | 0.1 | 0.85994466 +/- 0.00110688 | 240.85308 +/- 3.45997 | 5.4846752e+08 +/- 7.87899e+06 |
| N4 | end/mid | 0.1 | 0.2 | 0.78127035 +/- 0.00131847 | 207.24508 +/- 2.35748 | 4.7193582e+08 +/- 5.36843e+06 |
| N4 | mid/mid | 0 | 0.02 | 0.94876099 +/- 0.000703223 | 568.0624 +/- 8.00498 | 1.2935843e+09 +/- 1.82288e+07 |
| N4 | mid/mid | 0.02 | 0.05 | 0.91125488 +/- 0.000906999 | 290.40672 +/- 4.78299 | 6.6131042e+08 +/- 1.08918e+07 |
| N4 | mid/mid | 0.05 | 0.1 | 0.86432902 +/- 0.00109219 | 228.39466 +/- 3.36314 | 5.2009737e+08 +/- 7.6585e+06 |
| N4 | mid/mid | 0.1 | 0.2 | 0.79251099 +/- 0.00129335 | 187.37382 +/- 2.23071 | 4.2668523e+08 +/- 5.07974e+06 |
| N8 | end/end | 0 | 0.02 | 0.95782471 +/- 0.000641042 | 465.37735 +/- 7.2281 | 7.4935732e+08 +/- 1.16388e+07 |
| N8 | end/end | 0.02 | 0.05 | 0.92263794 +/- 0.000852107 | 269.48081 +/- 4.58224 | 4.3392188e+08 +/- 7.37839e+06 |
| N8 | end/end | 0.05 | 0.1 | 0.8764445 +/- 0.00104956 | 221.89037 +/- 3.29314 | 3.5729107e+08 +/- 5.30266e+06 |
| N8 | end/end | 0.1 | 0.2 | 0.80167643 +/- 0.00127175 | 192.60354 +/- 2.24732 | 3.1013299e+08 +/- 3.61866e+06 |
| N8 | end/mid | 0 | 0.02 | 0.95776367 +/- 0.000641485 | 466.06558 +/- 7.23356 | 7.5046551e+08 +/- 1.16476e+07 |
| N8 | end/mid | 0.02 | 0.05 | 0.92529297 +/- 0.000838561 | 248.33267 +/- 4.39566 | 3.9986884e+08 +/- 7.07795e+06 |
| N8 | end/mid | 0.05 | 0.1 | 0.88408407 +/- 0.00102102 | 196.81166 +/- 3.09248 | 3.1690897e+08 +/- 4.97956e+06 |
| N8 | end/mid | 0.1 | 0.2 | 0.81822713 +/- 0.00123003 | 167.21033 +/- 2.07867 | 2.6924448e+08 +/- 3.3471e+06 |
| N8 | mid/mid | 0 | 0.02 | 0.95932007 +/- 0.000630066 | 448.52948 +/- 7.09327 | 7.2222863e+08 +/- 1.14217e+07 |
| N8 | mid/mid | 0.02 | 0.05 | 0.92947388 +/- 0.000816597 | 227.5637 +/- 4.20137 | 3.6642634e+08 +/- 6.76511e+06 |
| N8 | mid/mid | 0.05 | 0.1 | 0.89245605 +/- 0.0009881 | 175.57097 +/- 2.91066 | 2.8270691e+08 +/- 4.68679e+06 |
| N8 | mid/mid | 0.1 | 0.2 | 0.83595785 +/- 0.00118109 | 141.26203 +/- 1.89583 | 2.2746215e+08 +/- 3.0527e+06 |
| N16 | end/end | 0 | 0.02 | 0.96385701 +/- 0.000595296 | 397.57306 +/- 6.67028 | 4.5267411e+08 +/- 7.59473e+06 |
| N16 | end/end | 0.02 | 0.05 | 0.9316508 +/- 0.000804837 | 244.69125 +/- 4.34894 | 2.7860387e+08 +/- 4.95168e+06 |
| N16 | end/end | 0.05 | 0.1 | 0.88881429 +/- 0.00100264 | 203.34131 +/- 3.13381 | 2.315231e+08 +/- 3.56814e+06 |
| N16 | end/end | 0.1 | 0.2 | 0.81838989 +/- 0.0012296 | 178.30682 +/- 2.1436 | 2.0301899e+08 +/- 2.44069e+06 |
| N16 | end/mid | 0 | 0.02 | 0.96365356 +/- 0.000596906 | 399.85296 +/- 6.68973 | 4.5526999e+08 +/- 7.61688e+06 |
| N16 | end/mid | 0.02 | 0.05 | 0.93492635 +/- 0.000786694 | 217.90152 +/- 4.10057 | 2.4810126e+08 +/- 4.66889e+06 |
| N16 | end/mid | 0.05 | 0.1 | 0.89900716 +/- 0.00096104 | 169.24359 +/- 2.84834 | 1.9269965e+08 +/- 3.2431e+06 |
| N16 | end/mid | 0.1 | 0.2 | 0.84040324 +/- 0.00116807 | 145.60382 +/- 1.9187 | 1.6578357e+08 +/- 2.18461e+06 |
| N16 | mid/mid | 0 | 0.02 | 0.96693929 +/- 0.000570256 | 363.09132 +/- 6.36934 | 4.1341343e+08 +/- 7.25209e+06 |
| N16 | mid/mid | 0.02 | 0.05 | 0.94191488 +/- 0.000746024 | 188.78981 +/- 3.80648 | 2.1495485e+08 +/- 4.33403e+06 |
| N16 | mid/mid | 0.05 | 0.1 | 0.91162109 +/- 0.000905307 | 141.22283 +/- 2.58799 | 1.607954e+08 +/- 2.94666e+06 |
| N16 | mid/mid | 0.1 | 0.2 | 0.86527507 +/- 0.00108897 | 112.70229 +/- 1.6699 | 1.283221e+08 +/- 1.90134e+06 |
| N4_N16 | end/end | 0 | 0.02 | 0.94609578 +/- 0.000720267 | 598.4438 +/- 8.22209 | 1.3627685e+09 +/- 1.87232e+07 |
| N4_N16 | end/end | 0.02 | 0.05 | 0.90726725 +/- 0.000925121 | 301.72863 +/- 4.88413 | 6.8709252e+08 +/- 1.11221e+07 |
| N4_N16 | end/end | 0.05 | 0.1 | 0.85509237 +/- 0.00112271 | 255.86309 +/- 3.57318 | 5.8264812e+08 +/- 8.1368e+06 |
| N4_N16 | end/end | 0.1 | 0.2 | 0.77234904 +/- 0.00133738 | 219.82951 +/- 2.43849 | 5.0059292e+08 +/- 5.5529e+06 |
| N4_N16 | end/mid | 0 | 0.02 | 0.94596354 +/- 0.000721099 | 599.9535 +/- 8.23274 | 1.3662063e+09 +/- 1.87475e+07 |
| N4_N16 | end/mid | 0.02 | 0.05 | 0.90635173 +/- 0.000929207 | 307.99138 +/- 4.93598 | 7.0135398e+08 +/- 1.12402e+07 |
| N4_N16 | end/mid | 0.05 | 0.1 | 0.85780843 +/- 0.0011139 | 237.80152 +/- 3.44286 | 5.4151854e+08 +/- 7.84003e+06 |
| N4_N16 | end/mid | 0.1 | 0.2 | 0.77993774 +/- 0.00132135 | 205.56007 +/- 2.35034 | 4.6809873e+08 +/- 5.35216e+06 |
| N4_N16 | mid/mid | 0 | 0.02 | 0.94837443 +/- 0.000705727 | 572.46357 +/- 8.03675 | 1.3036066e+09 +/- 1.83012e+07 |
| N4_N16 | mid/mid | 0.02 | 0.05 | 0.91053263 +/- 0.000910321 | 293.1815 +/- 4.80723 | 6.6762912e+08 +/- 1.0947e+07 |
| N4_N16 | mid/mid | 0.05 | 0.1 | 0.86362712 +/- 0.00109456 | 228.47894 +/- 3.36511 | 5.2028928e+08 +/- 7.66298e+06 |
| N4_N16 | mid/mid | 0.1 | 0.2 | 0.79378255 +/- 0.00129041 | 182.15612 +/- 2.19898 | 4.1480356e+08 +/- 5.00748e+06 |
| N16_N4 | end/end | 0 | 0.02 | 0.94847616 +/- 0.000705069 | 571.30519 +/- 8.0284 | 1.3009688e+09 +/- 1.82822e+07 |
| N16_N4 | end/end | 0.02 | 0.05 | 0.90884399 +/- 0.000918019 | 307.319 +/- 4.92394 | 6.9982284e+08 +/- 1.12127e+07 |
| N16_N4 | end/end | 0.05 | 0.1 | 0.85742188 +/- 0.00111516 | 251.61144 +/- 3.53941 | 5.7296631e+08 +/- 8.0599e+06 |
| N16_N4 | end/end | 0.1 | 0.2 | 0.77629598 +/- 0.00132912 | 214.6958 +/- 2.40512 | 4.8890249e+08 +/- 5.47691e+06 |
| N16_N4 | end/mid | 0 | 0.02 | 0.94621785 +/- 0.000719497 | 597.05042 +/- 8.21224 | 1.3595955e+09 +/- 1.87008e+07 |
| N16_N4 | end/mid | 0.02 | 0.05 | 0.90718587 +/- 0.000925485 | 303.30341 +/- 4.89681 | 6.9067858e+08 +/- 1.1151e+07 |
| N16_N4 | end/mid | 0.05 | 0.1 | 0.85835775 +/- 0.0011121 | 239.01002 +/- 3.45025 | 5.4427051e+08 +/- 7.85687e+06 |
| N16_N4 | end/mid | 0.1 | 0.2 | 0.78125 +/- 0.00131851 | 203.31166 +/- 2.33608 | 4.6297867e+08 +/- 5.31969e+06 |
| N16_N4 | mid/mid | 0 | 0.02 | 0.94856771 +/- 0.000704477 | 570.26276 +/- 8.02088 | 1.298595e+09 +/- 1.8265e+07 |
| N16_N4 | mid/mid | 0.02 | 0.05 | 0.91205851 +/- 0.000903281 | 282.59298 +/- 4.7174 | 6.4351708e+08 +/- 1.07424e+07 |
| N16_N4 | mid/mid | 0.05 | 0.1 | 0.86636353 +/- 0.00108524 | 222.04606 +/- 3.31338 | 5.0564041e+08 +/- 7.54518e+06 |
| N16_N4 | mid/mid | 0.1 | 0.2 | 0.79476929 +/- 0.00128812 | 186.30589 +/- 2.22145 | 4.2425337e+08 +/- 5.05866e+06 |
| floor | end/end | 0 | 0.02 | 0.98464966 +/- 0.000240121 | 167.06928 +/- 2.63373 | 1.8757012e+09 +/- 2.95692e+07 |
| floor | end/end | 0.02 | 0.05 | 0.97349548 +/- 0.000313731 | 82.027553 +/- 1.51696 | 9.2093038e+08 +/- 1.7031e+07 |
| floor | end/end | 0.05 | 0.1 | 0.95873642 +/- 0.000388475 | 65.996628 +/- 1.06103 | 7.4094981e+08 +/- 1.19122e+07 |
| floor | end/end | 0.1 | 0.2 | 0.93395996 +/- 0.000485063 | 56.554457 +/- 0.701762 | 6.3494175e+08 +/- 7.87874e+06 |
| floor | end/mid | 0 | 0.02 | 0.98464966 +/- 0.000240121 | 167.06928 +/- 2.63373 | 1.8757012e+09 +/- 2.95692e+07 |
| floor | end/mid | 0.02 | 0.05 | 0.97349548 +/- 0.000313731 | 82.027553 +/- 1.51696 | 9.2093038e+08 +/- 1.7031e+07 |
| floor | end/mid | 0.05 | 0.1 | 0.95873642 +/- 0.000388475 | 65.996628 +/- 1.06103 | 7.4094981e+08 +/- 1.19122e+07 |
| floor | end/mid | 0.1 | 0.2 | 0.93395996 +/- 0.000485063 | 56.554457 +/- 0.701762 | 6.3494175e+08 +/- 7.87874e+06 |
| floor | mid/mid | 0 | 0.02 | 0.98464966 +/- 0.000240121 | 167.06928 +/- 2.63373 | 1.8757012e+09 +/- 2.95692e+07 |
| floor | mid/mid | 0.02 | 0.05 | 0.97349548 +/- 0.000313731 | 82.027553 +/- 1.51696 | 9.2093038e+08 +/- 1.7031e+07 |
| floor | mid/mid | 0.05 | 0.1 | 0.95873642 +/- 0.000388475 | 65.996628 +/- 1.06103 | 7.4094981e+08 +/- 1.19122e+07 |
| floor | mid/mid | 0.1 | 0.2 | 0.93395996 +/- 0.000485063 | 56.554457 +/- 0.701762 | 6.3494175e+08 +/- 7.87874e+06 |

| Pair | Bead contacts | Mode contacts | Replicas | Combined SE of probabilities | Allowed difference | Pass |
|---|---:|---:|---:|---:|---:|---|
| end/end | 4173 | 4093 | 32768 | 0.002593692 | 0.01373256 | True |
| end/mid | 4072 | 4035 | 32768 | 0.002572203 | 0.01360662 | True |
| mid/mid | 3932 | 3878 | 32768 | 0.002531162 | 0.01337578 | True |

Sampled covariance: 75 checks; maximum absolute discrepancy / SE = 3.0215409 (limit 6).
Independent Ewald constant: 2.83729747948.

<!-- END REPRODUCED NUMBERS -->

## Two-part band and transient candidate diagnostics

<!-- BEGIN CANDIDATE NUMBERS -->

| Case | Pair | Candidate (m3 mol-1 s-1) | Candidate / reference | Absolute log deviation | Numerical band | Model band (owner-approved) | Numerical consumed | Model consumed | Residual model required | Excess | Qualified | Within zero-model band |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| N4 | end/end | 114463905 | 0.97200771 | 0.028392 | 0.234701 | 0.000000 | 0.028392 | 0.000000 | 0.000000 | 0.000000 | True | True |
| N4 | end/mid | 114463905 | 1.0307146 | 0.030252 | 0.239261 | 0.000000 | 0.030252 | 0.000000 | 0.000000 | 0.000000 | True | True |
| N4 | mid/mid | 114463905 | 1.0945407 | 0.090335 | 0.269891 | 0.000000 | 0.090335 | 0.000000 | 0.000000 | 0.000000 | True | True |
| N8 | end/end | 80938203.7 | 1.0311343 | 0.030659 | 0.222288 | 0.000000 | 0.030659 | 0.000000 | 0.000000 | 0.000000 | True | True |
| N8 | end/mid | 80938203.7 | 1.1145266 | 0.108430 | 0.307605 | 0.000000 | 0.108430 | 0.000000 | 0.000000 | 0.000000 | True | True |
| N8 | mid/mid | 80938203.7 | 1.2618038 | 0.232542 | 0.309482 | 0.000000 | 0.232542 | 0.000000 | 0.000000 | 0.000000 | True | True |
| N16 | end/end | 57231952.7 | 1.0433544 | 0.042441 | 0.308810 | 0.000000 | 0.042441 | 0.000000 | 0.000000 | 0.000000 | True | True |
| N16 | end/mid | 57231952.7 | 1.197671 | 0.180379 | 0.298198 | 0.000000 | 0.180379 | 0.000000 | 0.000000 | 0.000000 | True | True |
| N16 | mid/mid | 57231952.7 | 1.4118655 | 0.344912 | 0.402248 | 0.000000 | 0.344912 | 0.000000 | 0.000000 | 0.000000 | True | True |
| N4_N16 | end/end | 71539940.9 | 0.76849412 | 0.263322 | 0.191935 | 0.000000 | 0.191935 | 0.000000 | 0.071387 | 0.071387 | True | False |
| N4_N16 | end/mid | 71539940.9 | 0.86683362 | 0.142908 | 0.233171 | 0.000000 | 0.142908 | 0.000000 | 0.000000 | 0.000000 | True | True |
| N4_N16 | mid/mid | 71539940.9 | 0.93196819 | 0.070457 | 0.209097 | 0.000000 | 0.070457 | 0.000000 | 0.000000 | 0.000000 | True | True |
| N16_N4 | end/end | 71539940.9 | 0.7682665 | 0.263619 | 0.202176 | 0.000000 | 0.202176 | 0.000000 | 0.061442 | 0.061442 | True | False |
| N16_N4 | end/mid | 71539940.9 | 0.81316263 | 0.206824 | 0.223093 | 0.000000 | 0.206824 | 0.000000 | 0.000000 | 0.000000 | True | True |
| N16_N4 | mid/mid | 71539940.9 | 0.92581296 | 0.077083 | 0.205939 | 0.000000 | 0.077083 | 0.000000 | 0.000000 | 0.000000 | True | True |
| floor | end/end | 282167444 | 0.99642882 | 0.003578 | 0.245929 | 0.000000 | 0.003578 | 0.000000 | 0.000000 | 0.000000 | True | True |
| floor | end/mid | 282167444 | 0.99642882 | 0.003578 | 0.245929 | 0.000000 | 0.003578 | 0.000000 | 0.000000 | 0.000000 | True | True |
| floor | mid/mid | 282167444 | 0.99642882 | 0.003578 | 0.245929 | 0.000000 | 0.003578 | 0.000000 | 0.000000 | 0.000000 | True | True |

Consumption is bookkeeping: numerical allowance is used first, then model allowance. It does not estimate a physical model-error component.

| Case | Pair | t1 | t2 | Candidate / finite-box reference | 20% + 4 SE diagnostic (REPORT-ONLY) |
|---|---|---:|---:|---:|---|
| N4 | end/end | 0 | 0.02 | 0.08563789 | False (report-only) |
| N4 | end/end | 0.02 | 0.05 | 0.1612886 | False (report-only) |
| N4 | end/end | 0.05 | 0.1 | 0.2025701 | False (report-only) |
| N4 | end/end | 0.1 | 0.2 | 0.2320121 | False (report-only) |
| N4 | end/mid | 0 | 0.02 | 0.08715289 | False (report-only) |
| N4 | end/mid | 0.02 | 0.05 | 0.1672918 | False (report-only) |
| N4 | end/mid | 0.05 | 0.1 | 0.2086977 | False (report-only) |
| N4 | end/mid | 0.1 | 0.2 | 0.2425413 | False (report-only) |
| N4 | mid/mid | 0 | 0.02 | 0.08848585 | False (report-only) |
| N4 | mid/mid | 0.02 | 0.05 | 0.1730865 | False (report-only) |
| N4 | mid/mid | 0.05 | 0.1 | 0.2200817 | False (report-only) |
| N4 | mid/mid | 0.1 | 0.2 | 0.2682631 | False (report-only) |
| N8 | end/end | 0 | 0.02 | 0.1080102 | False (report-only) |
| N8 | end/end | 0.02 | 0.05 | 0.1865271 | False (report-only) |
| N8 | end/end | 0.05 | 0.1 | 0.226533 | False (report-only) |
| N8 | end/end | 0.1 | 0.2 | 0.260979 | False (report-only) |
| N8 | end/mid | 0 | 0.02 | 0.1078507 | False (report-only) |
| N8 | end/mid | 0.02 | 0.05 | 0.2024119 | False (report-only) |
| N8 | end/mid | 0.05 | 0.1 | 0.2553989 | False (report-only) |
| N8 | end/mid | 0.1 | 0.2 | 0.3006123 | False (report-only) |
| N8 | mid/mid | 0 | 0.02 | 0.1120673 | False (report-only) |
| N8 | mid/mid | 0.02 | 0.05 | 0.2208853 | False (report-only) |
| N8 | mid/mid | 0.05 | 0.1 | 0.2862972 | False (report-only) |
| N8 | mid/mid | 0.1 | 0.2 | 0.3558315 | False (report-only) |
| N16 | end/end | 0 | 0.02 | 0.1264308 | False (report-only) |
| N16 | end/end | 0.02 | 0.05 | 0.2054241 | False (report-only) |
| N16 | end/end | 0.05 | 0.1 | 0.2471976 | False (report-only) |
| N16 | end/end | 0.1 | 0.2 | 0.2819044 | False (report-only) |
| N16 | end/mid | 0 | 0.02 | 0.1257099 | False (report-only) |
| N16 | end/mid | 0.02 | 0.05 | 0.2306798 | False (report-only) |
| N16 | end/mid | 0.05 | 0.1 | 0.2970008 | False (report-only) |
| N16 | end/mid | 0.1 | 0.2 | 0.3452209 | False (report-only) |
| N16 | mid/mid | 0 | 0.02 | 0.1384376 | False (report-only) |
| N16 | mid/mid | 0.02 | 0.05 | 0.266251 | False (report-only) |
| N16 | mid/mid | 0.05 | 0.1 | 0.3559303 | False (report-only) |
| N16 | mid/mid | 0.1 | 0.2 | 0.4460023 | False (report-only) |
| N4_N16 | end/end | 0 | 0.02 | 0.05249603 | False (report-only) |
| N4_N16 | end/end | 0.02 | 0.05 | 0.1041198 | False (report-only) |
| N4_N16 | end/end | 0.05 | 0.1 | 0.1227841 | False (report-only) |
| N4_N16 | end/end | 0.1 | 0.2 | 0.1429104 | False (report-only) |
| N4_N16 | end/mid | 0 | 0.02 | 0.05236394 | False (report-only) |
| N4_N16 | end/mid | 0.02 | 0.05 | 0.1020026 | False (report-only) |
| N4_N16 | end/mid | 0.05 | 0.1 | 0.1321099 | False (report-only) |
| N4_N16 | end/mid | 0.1 | 0.2 | 0.1528309 | False (report-only) |
| N4_N16 | mid/mid | 0 | 0.02 | 0.05487847 | False (report-only) |
| N4_N16 | mid/mid | 0.02 | 0.05 | 0.1071552 | False (report-only) |
| N4_N16 | mid/mid | 0.05 | 0.1 | 0.1375003 | False (report-only) |
| N4_N16 | mid/mid | 0.1 | 0.2 | 0.172467 | False (report-only) |
| N16_N4 | end/end | 0 | 0.02 | 0.05498974 | False (report-only) |
| N16_N4 | end/end | 0.02 | 0.05 | 0.1022258 | False (report-only) |
| N16_N4 | end/end | 0.05 | 0.1 | 0.1248589 | False (report-only) |
| N16_N4 | end/end | 0.1 | 0.2 | 0.1463276 | False (report-only) |
| N16_N4 | end/mid | 0 | 0.02 | 0.05261855 | False (report-only) |
| N16_N4 | end/mid | 0.02 | 0.05 | 0.1035792 | False (report-only) |
| N16_N4 | end/mid | 0.05 | 0.1 | 0.1314419 | False (report-only) |
| N16_N4 | end/mid | 0.1 | 0.2 | 0.154521 | False (report-only) |
| N16_N4 | mid/mid | 0 | 0.02 | 0.05509026 | False (report-only) |
| N16_N4 | mid/mid | 0.02 | 0.05 | 0.1111702 | False (report-only) |
| N16_N4 | mid/mid | 0.05 | 0.1 | 0.1414838 | False (report-only) |
| N16_N4 | mid/mid | 0.1 | 0.2 | 0.1686255 | False (report-only) |
| floor | end/end | 0 | 0.02 | 0.150433 | False (report-only) |
| floor | end/end | 0.02 | 0.05 | 0.3063939 | False (report-only) |
| floor | end/end | 0.05 | 0.1 | 0.3808186 | False (report-only) |
| floor | end/end | 0.1 | 0.2 | 0.4443989 | False (report-only) |
| floor | end/mid | 0 | 0.02 | 0.150433 | False (report-only) |
| floor | end/mid | 0.02 | 0.05 | 0.3063939 | False (report-only) |
| floor | end/mid | 0.05 | 0.1 | 0.3808186 | False (report-only) |
| floor | end/mid | 0.1 | 0.2 | 0.4443989 | False (report-only) |
| floor | mid/mid | 0 | 0.02 | 0.150433 | False (report-only) |
| floor | mid/mid | 0.02 | 0.05 | 0.3063939 | False (report-only) |
| floor | mid/mid | 0.05 | 0.1 | 0.3808186 | False (report-only) |
| floor | mid/mid | 0.1 | 0.2 | 0.4443989 | False (report-only) |

Literal transport anchors: 13/13.
Literal min/floor/two-chain branches: 18/18.
Strict zero-model scientific comparison passed: False.
Model tolerance status: APPROVED by owner on 2026-10-04.

<!-- END CANDIDATE NUMBERS -->

## Every mutation result

<!-- BEGIN MUTATION NUMBERS -->

Target baseline passes literal-input and kernel-rule controls before each mutation. Scientific adoption failure is never counted as a mutation kill.
Unmodified strict scientific comparison passed: False; qualified baseline rate rows passing: 16/18.

| Actual target mutation | Changes observable | Rejected | Failed anchors | Failed literal branches | Failed qualified rate checks | New rate rejections from passing baseline rows |
|---|---|---|---:|---:|---:|---:|---:|
| capture_radius_x0.5 | True | True | 0 | 18 | 17 | 15 |
| capture_radius_x2 | True | True | 0 | 18 | 18 | 16 |
| candidate_diffusion_x4 | True | True | 0 | 18 | 18 | 16 |
| shared_chain_diffusivity_x4 | True | True | 4 | 18 | 18 | 16 |
| drop_one_chain_diffusion | True | True | 0 | 18 | 17 | 15 |
| length_min_to_max | True | True | 0 | 6 | 6 | 4 |
| length_min_to_i | True | True | 0 | 3 | 4 | 2 |
| length_min_to_j | True | True | 0 | 3 | 4 | 2 |
| length_min_to_mean | True | True | 0 | 6 | 5 | 4 |
| drop_sigma0_floor | True | True | 0 | 3 | 2 | 0 |
| swap_reference_prefactors_in_candidate | True | True | 0 | 10 | 4 | 2 |
| swap_reference_prefactors_N4_only_diagnostic | True | True | 0 | 2 | 2 | 0 |

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

| Sampler mutation | Rejected | Failure stage |
|---|---|---|
| double_actual_OU_noise | True | actual sampled stationary variance failed |
| double_sampled_internal_coordinates | True | sampled covariance failed |
| double_OU_noise_N16_only_in_exact_adoption_test | True | stationary_variance |

The N16-only noise mutant executes the exact adoption-test body and must fail its stationary-variance audit before any trajectory generation. Its positive control is the complete unmodified audit phase, not an already-failing scientific comparison.

The required prefactor mutation rescales candidate end/end by C_mid/C_end and mid/mid by C_end/C_mid, using the independent reference. A bare target label swap is a separate observational no-op; it is reported as an expected survivor, never substituted for the requested prefactor perturbation.

<!-- END MUTATION NUMBERS -->

Formula/anchor conformance and independent qualified rate discrimination are
reported separately. Global reference-prefactor reversal creates only two new
scientific failures, N8 and N16 mid/mid; its other failures already existed in
the baseline. N4-only reversal and floor removal create no new scientific
failures. Their rejection by exact formula checks establishes conformance,
not independent physical discrimination. The N4-only reversal remains an
additional diagnostic. A bare target-class label swap is an expected survivor because
the target has a common multiplier. The required reference-derived prefactor
swap is substantive and rescaling the candidate, not relabelling it.

## Directly runnable proposed test and full Verifier

The default `*Test.py` discovery collects the build entry at
`test/rmgpy/kmc/kernelReferenceBuildTest.py`. It validates its fixed sibling
`fixtures/i030_kernel_reference` path during collection and fails with
FileNotFoundError if a required fixture file is missing. The multi-hour test
is marked slow and runs with `RMG_KMC_SLOW=1`; otherwise only that expensive
execution is explicitly skipped after fixture validation. The entry forces
the fixed fixture path when delegating to `proposed_adoption_test.py`.
The owner-approved gate asserts the sixteen ordinary rows and reports exactly
the two named known deviations, including their measured ratios and approval date.

The standalone `proposed_adoption_test.py` is still runnable by explicit file
path. Its root defaults to its own directory; an explicit unprovisioned override
skips with a reason. The collected build entry cannot use that provisioning skip.
Acceptance and audit helpers now use explicit failures rather than removable
assert statements. Complete sampler qualification refuses optimized execution,
so adoption under `python -O` or `PYTHONOPTIMIZE` stops before a campaign and
cannot bypass the unchanged reference's removable assertions.
The Verifier compares this actual file's function byte-for-byte with the body
below, then executes the repository file with `MET_KERNEL_REFERENCE_PACK`
**unset**. It retains complete sampler-audit results written inside that test.

```python
def test_met_kernel_reference(tmp_path):
    import importlib.util
    import json
    import os
    from pathlib import Path
    import subprocess
    import sys
    import warnings
    sys.dont_write_bytecode = True
    root = Path(os.environ.get("MET_KERNEL_REFERENCE_PACK", str(Path(__file__).resolve().parent)))
    if not (root / "reference/run.py").is_file():
        import pytest
        pytest.skip("MET kernel reference fixture is not provisioned at " + str(root))
    sys.path.insert(0, str(root))
    import checks
    import audit_sampler
    spec = importlib.util.spec_from_file_location("test_rouse_inputs", root / "reference/run.py")
    oracle = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(oracle)
    cfg, digest = oracle.read_parameters()
    # Complete preflight, including innovation-noise stationary variance.
    audits = audit_sampler.run_audits(cfg)
    (tmp_path / "sampler_audits.json").write_text(json.dumps(audits, indent=2) + "\n")
    output = tmp_path / "fresh-reference"
    subprocess.run([
        sys.executable, "-B", "-u", str(root / "reference/run.py"),
        "--workers", os.environ.get("MET_KERNEL_REFERENCE_WORKERS", "8"),
        "--output", str(output),
    ], check=True)
    reference = json.loads((output / "results.json").read_text())
    assert reference["parameters_sha256"] == digest
    report = checks.target_checks(checks.load_target(), cfg, reference)
    (tmp_path / "candidate_results.json").write_text(json.dumps(report, indent=2) + "\n")
    adoption = checks.assert_owner_approved_adoption(report, cfg)
    (tmp_path / "owner_adoption.json").write_text(json.dumps(adoption, indent=2) + "\n")
    for deviation in adoption["known_deviations"]:
        message = "OWNER-APPROVED KNOWN DEVIATION: " + json.dumps(deviation, sort_keys=True)
        print(message, flush=True)
        warnings.warn(message, RuntimeWarning)
```

Commands use `/home/alon/anaconda3/envs/rmg_env/bin/python`; long runs tee both
streams and pin all descendants to eight CPUs. The full reproduction command is:

```bash
taskset -c 6,16,17,19,21,23,24,26 /home/alon/anaconda3/envs/rmg_env/bin/python -B -u test/rmgpy/kmc/fixtures/i030_kernel_reference/verify_pack.py --workers 8 --require-all-mutations --output /home/alon/runs/i030-met-kernel-reference/rework2/verifications > >(tee -a /home/alon/runs/i030-met-kernel-reference/rework2/verifier-stdout.log) 2> >(tee -a /home/alon/runs/i030-met-kernel-reference/rework2/verifier-stderr.log >&2)
```

The default exit checks exact numerical reproduction and the owner-approved
binding test: sixteen passing rows and two explicitly reported deviations.
`--require-acceptance` requires this binding outcome. The strict zero-model
scientific comparison still reports its two failures. `mutation_tests.py --require-all` requires substantive target controls,
their independent rate controls and every sampler mutation. Baseline science
failure is never counted as a kill. `--mutations-only --require-all-mutations`
is a separate fast Verifier mode after the report is frozen.

Every authoring and full Verifier campaign generates all371 physical ensembles
fresh, with no reuse input. The normal reference has seven settings, three
seed stages for each non-floor case and eight for the floor, totaling14,327,808
trajectories. Actual histories must use current input/reference-program bytes.
The old A-02 provenance archive and `extend_floor.py` are historical, unused in
A-03, and do not authorize replay of the former incorrect guard.

Each verification uses a UUID staging directory. Its complete JSON, Markdown,
runtime, hashes/versions, exact test, mutant evidence and JUnit are published
by one same-filesystem directory rename. `COMPLETE.json` hashes every published
file. Failed pending runs are distinguishable from complete publications.
No shared latest pointer or fixture output mutation exists.

## Responses to dispatch items1–6

1. **Adaptive guard:** repeated exact evaluation and halving, with an equilibrium
   regression across every case/class and an explicit minimum-step exception.
2. **Transient diagnostics:** restored four short-time bins, SEs and original
  20%+4SE candidate diagnostics; all report-only.
3. **Adoption sampler qualification:** all four independent audits, including
   stationary variance, execute inside the exact test. Length16-only noise is
   rejected by that actual pytest before any campaign, with a passing audit
   positive control; other sampler controls remain separate.
4. **Two-part tolerance:** joint numerical sampling envelope plus explicit
   owner-approved model allowance0; the declaration's committed bytes and self-recorded
   measurement/comparison order are checked, with every row's consumption reported.
5. **Adoptability:** default discovery collects a build entry with a fixed fixture
   path and fatal missing fixtures. Its full campaign is explicitly opt-in;
   acceptance gates remain binding under optimization, which full adoption refuses.
6. **Interpretation:** finite N<=16 only, no asymptotic/material validation claim;
   shared continuum bias, volume/time remainders, source-convention dependence,
   weak contact/coil separation and conformance-only mutation limits are explicit.

## Limits, reproduced evidence, compute and remaining work

**Limits.** Microscopic capture/coil ratios at equal N4/8/16 are approximately
1.233/.872/.616; increasing coverage to shorter chains does not create scale
separation. Finite bond mapping changes with N. Coil-scale reaction theory
([O'Shaughnessy–Vavylonis](https://arxiv.org/abs/cond-mat/9805331)) does not supply
these discrete amplitudes or a model tolerance. The empirical adjacent-setting
contrasts do not bound a common continuum bias. Taking maximum adjacent-box
contrast is not an extrapolation or a demonstrated infinite-volume remainder.
The finite-window plateau does not prove an infinite-time limit. Enlarged
boxes and the independent contact audits reduce identifiable defects but do
not remove these interpretation limits. Computational repeatability is not
physical correctness. Exact numerical equality is conditional on the recorded
software/platform, not promised across architectures or versions.

**Retained execution evidence (A-03, before round-3 gate changes).** The retained
external logs report that the full fresh Verifier command above exited0 at
`d94913cefb528544c5da6c0a6001ecb2c8c50bdb`.
The scientific tables and JSON reports are committed. Full-run JUnit, runtime
totals and publication fingerprints are external evidence at the paths below,
not independently established by those committed tables or self-recorded timestamps.
The round-3 follow-up runs fast regressions only; it does not rerun the campaign
or certify a fresh full run of the changed gate code.
It reproduced every scientific JSON value, trajectory/contact hash and rendered
table from371 freshly generated physical ensembles (14,327,808 trajectories),
pooled into126 reported ensembles including the declared floor-label aliases.
All18/18 reference rows qualified; 16/18 candidate rows fell
within the numerical plus proposed-zero-model band. The real proposed scientific
pytest outcome was **FAIL**, retained in `adoption-junit.xml` (one test;
failures=1, errors=0, skips=0).
This is neither owner approval nor an expected failure. The exact test ran from
the repository fixture with the pack-root variable unset, and complete sampler
audits written inside it matched the independently retained audit report.

The full Verifier reran all21 equilibrium-guard/uncertainty tests:21 passed,
zero failures/errors/skips. The historical A-02 uncertainty probe reproduced
exactly from its actual Git commit. Plain repository pytest collection worked
without a root variable; an explicitly unprovisioned override separately
produced one skip with the reason `MET kernel reference fixture is not provisioned`.

The atomic publication is
`/home/alon/runs/i030-met-kernel-reference/rework2/verifications/verified-20261004T131106-53f863d1e01b4372bd3aaad771765915`.
Every one of its390 `COMPLETE.json` file fingerprints was independently
checked, with zero mismatches. `verification.json` retains current input/program,
policy/chronology and unchanged target fingerprints, Python/NumPy/SciPy/pytest
versions, and the real generation history. The target SHA-256 remains
`1a8176c2e49c3fcfb5909c73bf01ad8d7230f2e5ad724e7461984a7b8388a7c1`.

The standalone `mutation_tests.py --require-all --output
/home/alon/runs/i030-met-kernel-reference/rework2/mutation-final.json` also exited0,
with both streams tee'd to the correspondingly named stdout/stderr logs. Its
complete JSON matched the full Verifier and committed mutation report.
Passing literal controls precede every target mutation; the full baseline
scientific comparison passed=False, so a pre-existing scientific
failure is never counted as a new rate rejection.

| Target source mutation | Rejected by anchors/formula | New qualified rate rejections from passing rows |
|---|---|---:|
| capture_radius_x0.5 | True | 15 |
| capture_radius_x2 | True | 16 |
| candidate_diffusion_x4 | True | 16 |
| shared_chain_diffusivity_x4 | True | 16 |
| drop_one_chain_diffusion | True | 15 |
| length_min_to_max | True | 4 |
| length_min_to_i | True | 2 |
| length_min_to_j | True | 2 |
| length_min_to_mean | True | 4 |
| drop_sigma0_floor | True | 0 |
| swap_reference_prefactors_in_candidate | True | 2 |
| swap_reference_prefactors_N4_only_diagnostic | True | 0 |

| Actual sampler mutation | Rejected | Failure stage |
|---|---|---|
| double_actual_OU_noise | True | actual sampled stationary variance failed |
| double_sampled_internal_coordinates | True | sampled covariance failed |
| double_OU_noise_N16_only_in_exact_adoption_test | True | stationary_variance |

The bare label-only diagnostic is an expected survivor=True.
The N16-only noise mutant's actual pytest/JUnit and generated exact-body test
are retained under `mutation-artifacts/`; it fails at stationary variance before
any trajectory output. Its unmodified complete audit phase passes. Model-band
consumption and every transient diagnostic remain in the tables and full JSON.

Measured numerical log bands range **0.1919354733–0.4022475628**;
the historical proposed model allowance was **0** for every row; that value
is now owner-approved with the two explicit deviation records.
Largest-box corrections remove **7.9096–20.5173%**.
Numerical uncertainty is the declared joint finite-setting normal/delta
envelope; no continuum/common-bias remainder or universal prefactor is certified.

| A-03 run | Wall hours | CPU hours (parent+children) | Child peak RSS MiB | Parent peak RSS MiB |
|---|---:|---:|---:|---:|
| Corrected author campaign | 3.74803 | 25.41922 | 87.26172 | 117.52344 |
| Fresh Verifier reference | 3.43442 | 24.22893 | 86.79297 | 117.65625 |

These reference processes used **49.64815 measured CPU hours**,
each with eight workers, one numerical-library thread per worker, and the same
eight pinned CPUs. The total excludes controller, unit/audit/standalone-mutation
overhead and prior attempts. Per-process RSS peaks are recorded; instantaneous
aggregate RSS was not measured. Python/NumPy/SciPy/pytest versions were
3.9.23/1.26.4/1.13.1/8.4.1 on Linux x86_64.

**Owner ruling and remaining work:** approval of zero model tolerance and the
two named deviations is recorded for 2026-10-04. The strict scientific table
continues to show those two failures; binding adoption requires sixteen passing
rows and the exact reported deviation records. The production run will carry
the documented 0.77x termination-channel sensitivity case. This follow-up uses
fast tests and the frozen reference; it does not claim a newly rerun full campaign.
Production execution remains with the manager. No additional owner approval
is pending for this ruling.
