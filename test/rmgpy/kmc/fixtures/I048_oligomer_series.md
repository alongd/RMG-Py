# I048 — polystyrene oligomer per-unit increment probe

This report distinguishes numerical replay of a declared approximate partition
from physical convergence and transfer to a radical polymer chain. No correction
is adopted. Only probe scripts and this report are changed.
The [reproduction commands](i048_probe/README.md) include the full
`verify_results.py --replay` source audit and numerical replay.

## Declared reductions and reference definitions

The previous probe's methyl cap convention is retained:
`CH3–[CH(Ph)–CH2]_(n−1)–CH(Ph)–CH3`, formula `C(8n+1)H(8n+4)`.
This differs from literal hydrogen-capped `H–(CH2–CHPh)n–H` by the same
terminal carbon at every length; successive increments still add C8H8.
The existing n=2 and 3 molecular evidence and the four one-ring reference
molecules are copied from the prior named scratch directory, with byte hashes
and source provenance. The n=4 and 5 inputs and searches are new.

All distinct stereo classes through n=5 are enumerated from all oriented
binary assignments and weighted as a frozen unbiased atactic chain.
Enantiomer equivalence does not imply an equilibrium mixture of sequence
classes. External rotational symmetry, local optical factors, global mirror
partners and chain configuration entropy are kept separate. Existing RMG
propagation configuration entropy is retained once.

The original reduction declaration is retained in
`method_plan_before_local_clarification.json`; `method_plan_i048.json` contains
the same reductions plus the recorded partition/reference clarifications.
After replaying the earlier short-fragment cycle, the entropy comparison was
clarified to use its absolute increment and the H reference was held constant
at 298 K. These clarifications preserve H/S consistency and precede any Tc
calculation; they do not change sampling, candidate choices, or energy levels.
The reductions declared before increment calculations
are: CREST/GFN2 `--quick`, 6 kcal/mol search window; up to 32 candidates within
12 kJ/mol for each new sequence, chosen by GFN2 energy; tight extreme xTB
optimization and Hessians, symmetry deduplication and proper/reflection orbit
completion; 1024 and 2048 correlated-importance points per basin; DFT on the
lowest three chemical minima for each larger sequence. Other retained minima
use GFN2 energies plus the functional offset of the lowest GFN2 conformer.
The declared levels are PBE-D3(BJ) and BLYP-D3(BJ)/def2-SVP//GFN2, density
fitting, grid level 3, SCF tolerance 1e-9 Eh, no optional three-body D3 term.
No reduction or ranking uses Tc or a thermochemical comparator.
The same already-declared sparse energy rule was subsequently evaluated on
the complete short-chain evidence as a diagnostic of the treatment boundary.
It does not replace any production energy or increment; all retained prior
and new source evidence is unchanged.

The exact prior fixed-angle coupled-rotor model is reused. Every acyclic
heavy-atom single bond is a torsional coordinate. Projection removes that
subspace before non-torsional harmonic vibrations are counted. The rigid
torsional Hessian supplies the Pitzer–Gwinn harmonic reference. Periodic
Voronoi basins partition the full torsional domain once, with methyl and
phenyl fundamental periods. Severe contacts below 0.45 Å have zero weight
and remain in the Monte Carlo denominator; other unresolved SCC failures
stop the producer. The stochastic CREST search is not proved converged.

## Reference and the meaning of local

Electronic energies of different atom counts cannot be subtracted from GAV
formation enthalpies directly. Following the prior probe, use the balanced
cycle `C_n + cumene + n-propylbenzene → C_(n+1) + ethylbenzene + ethane`.
Its additive source vector cancels at every length. The balanced QM-minus-GAV
H residual at 298.15 K anchors the formation-enthalpy increment to the pinned
GAV one-ring reference; the anchor is held constant and the chain thermal
increment supplies its temperature dependence. Absolute entropy increments
need no electronic-energy reference and are compared directly to oriented
GAV. This preserves the H/S derivative identity. The temperature-dependent
balanced-cycle H/S residuals are retained as separate diagnostics reproducing
the prior probe's comparison. This is not an independently determined absolute
formation enthalpy or an explicit styrene-addition calculation.

The raw chain partition includes ideal-gas translation and external rotation
at 1 bar. The primary local partition omits both factors before basin
probabilities are calculated. Thus it describes internal chain conformations,
rather than retaining an inertia bias from gas rotational weights. The gas
one-ring reference used for the H anchor remains fixed. Component subtraction at the gas weights
and removal from every reference molecule are reported as diagnostics.
GAV chain entropies are normalized to an oriented fixed sequence by undoing
finite-end symmetry and removing its constitutional optical factors; no new
R ln 2 is inserted. JSON preserves every molecular component and channel.
The identical additive increment at each length defines the group model's
per-unit baseline. The local residual compares that fixed increment with
the internal QM chain increment; the raw residual retains the measured
finite-chain external factors.
The total residual also contains errors in intrarepeat vibrations and torsions.
The balanced-cycle diagnostics retain the original GAV cycle reference; their
different quantum reference partitions describe separate comparisons. The
direct local increment relative to the fixed GAV per-unit baseline supplies
the conditional chain-transfer calculation.

Whole-chain translation and rotation change as mass and moments of inertia
change; their increments tend to zero in a long-chain limit. It is therefore
the chain external increment that can be removed to estimate that limit.
The prior cycle's negative translation/rotation terms largely belong to its
small-molecule reference. Deleting those terms changes the reference, rather
than merely removing a finite-chain artifact. This is an evidenced problem
with the motivating interpretation in the dispatch.

The conditional gas roots solve
`ΔH_p(T)+δH(T) − T[ΔS_p(T)+δS(T)] − RT ln(c0 RT/p0) = 0`,
with c0=1000 mol/m³ and p0=100000 Pa. Only the measured 298–800 K domain is
used. The pinned propagation H/S are evaluated on a one-kelvin grid and
interpolated for root finding. Each finite-fragment root is a comparison
under its transfer assumption, not a validated long-chain ceiling.

## Measured evidence

<!-- BEGIN I048:measured -->
Pinned baseline Tc: **710.020462 K** at 1 mol/L. Database commit `4a12d36fcdc193ede82c8d1ab5c1653495d445bc`; snapshot `61e2c7357bb4945dce6e6a1a39eccf216e945a01771e461391e7975165ea29eb`, 65 allowlisted files.

### Sequence and molecular coverage

Thermochemical coverage: **25/25 species**. Missing completed thermochemistry: none.
An increment is calculated only when every frozen class at both lengths and every reference is complete. Missing classes are never omitted and the remaining weights are never renormalized.
New-case source coverage: 16/16 completed searches and 16/16 selected minimum pools. The table retains available search and minimum evidence even when a case lacks completed thermochemistry.

| n | Class | Frozen atactic weight | Oriented assignments | CREST candidates | Within 12 kJ/mol | Checked candidates | Rotor wells | Minimum ESS at 298 K | Largest-population basin ESS: PBE / BLYP |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 2 | diphenylpentane_meso | 0.500000 | 2 | 25 | 14 | 14 | 10 | 92.29 | 445.17 / 588.45 |
| 2 | diphenylpentane_racemo | 0.500000 | 2 | 25 | 13 | 13 | 13 | 47.04 | 642.18 / 642.18 |
| 3 | triphenylheptane_iso | 0.250000 | 2 | 115 | 31 | 31 | 22 | 11.99 | 112.87 / 112.87 |
| 3 | triphenylheptane_hetero | 0.500000 | 4 | 142 | 53 | 53 | 26 | 12.76 | 199.43 / 199.43 |
| 3 | triphenylheptane_syndio | 0.250000 | 2 | 110 | 46 | 46 | 31 | 2.94 | 117.53 / 117.53 |
| 4 | ps4_0000 | 0.125000 | 2 | 351 | 87 | 32 | 24 | 3.07 | 51.60 / 51.60 |
| 4 | ps4_0001 | 0.250000 | 4 | 339 | 113 | 32 | 7 | 2.42 | 2.42 / 2.42 |
| 4 | ps4_0010 | 0.250000 | 4 | 353 | 122 | 32 | 18 | 7.46 | 22.87 / 22.87 |
| 4 | ps4_0011 | 0.125000 | 2 | 364 | 98 | 32 | 19 | 3.66 | 4.63 / 4.63 |
| 4 | ps4_0101 | 0.125000 | 2 | 259 | 96 | 32 | 17 | 4.76 | 40.05 / 40.05 |
| 4 | ps4_0110 | 0.125000 | 2 | 216 | 72 | 32 | 31 | 3.45 | 8.23 / 8.23 |
| 5 | ps5_00000 | 0.062500 | 2 | 325 | 73 | 32 | 36 | 1.03 | 1.03 / 1.03 |
| 5 | ps5_00001 | 0.125000 | 4 | 390 | 129 | 32 | 13 | 3.39 | 13.73 / 13.73 |
| 5 | ps5_00010 | 0.125000 | 4 | 682 | 125 | 32 | 23 | 1.05 | 1.05 / 1.05 |
| 5 | ps5_00011 | 0.125000 | 4 | 382 | 117 | 32 | 15 | 1.55 | 1.55 / 1.55 |
| 5 | ps5_00100 | 0.062500 | 2 | 554 | 129 | 32 | 29 | 2.08 | 6.61 / 6.61 |
| 5 | ps5_00101 | 0.125000 | 4 | 525 | 185 | 32 | 17 | 1.78 | 12.85 / 12.85 |
| 5 | ps5_00110 | 0.125000 | 4 | 398 | 134 | 32 | 19 | 1.16 | 1.16 / 1.16 |
| 5 | ps5_01001 | 0.125000 | 4 | 377 | 113 | 32 | 22 | 2.93 | 49.60 / 49.60 |
| 5 | ps5_01010 | 0.062500 | 2 | 824 | 266 | 32 | 19 | 3.37 | 38.64 / 38.64 |
| 5 | ps5_01110 | 0.062500 | 2 | 407 | 122 | 32 | 34 | 1.19 | 32.03 / 40.55 |

All stereochemical classes are covered through n=5. Complete lengths are averaged with their frozen assignment weights, without equilibrium diastereomer mixing entropy. The population spread below measures variation among oriented addition channels; it is not a standard error on these exact weighted means.

Monte Carlo standard errors condition on the retained wells and declared partition model. They do not quantify omitted conformers, candidate truncation, rigidity of the coupled-rotor potential, or the electronic fallback. Sequence population spread is reported separately. These quantities alone therefore do not provide a total uncertainty or prove physical convergence; low effective sample sizes also limit the linearized Monte Carlo error estimate.
The minimum effective sample size can belong to a basin of small probability. The largest-population basin column and the per-basin internal probabilities saved in JSON provide the population context; none of these diagnostics changes the declared quadrature.

### Increment and residual table at 298.15 K

H increments use the one-ring formation-enthalpy anchor at 298 K described above, followed by the chain thermal increment. S increments are absolute molecular differences; GAV S and the raw comparison receive matching oriented-end normalization. The component tables retain the original molecular raw sums. Local values come from the internal chain partition. Balanced-cycle residuals are separately retained in JSON and in the component diagnostics.

| Step | Level | GAV ΔH (kJ/mol) | GAV oriented ΔS (J/mol/K) | Calibrated raw ΔH | Raw oriented ΔS | Calibrated local ΔH | Internal ΔS | Raw δH | Raw δS | Local δH | Local δS | Local MC SE H/S | Sequence population SD H/S |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| 2→3 | pbe | 67.381 | 191.908 | 57.507 | 159.751 | 57.292 | 145.815 | -9.874 | -32.158 | -10.088 | -46.093 | 0.179 / 0.532 | 4.797 / 5.279 |
| 2→3 | blyp | 67.381 | 191.908 | 51.758 | 158.507 | 51.750 | 144.958 | -15.623 | -33.401 | -15.630 | -46.950 | 0.183 / 0.545 | 4.511 / 5.401 |
| 3→4 | pbe | 67.381 | 191.908 | 65.492 | 159.010 | 65.638 | 149.460 | -1.889 | -32.899 | -1.743 | -42.448 | 0.380 / 1.429 | 8.869 / 6.386 |
| 3→4 | blyp | 67.381 | 191.908 | 57.843 | 160.556 | 57.929 | 150.606 | -9.537 | -31.352 | -9.452 | -41.302 | 0.377 / 1.418 | 6.513 / 7.175 |
| 4→5 | pbe | 67.381 | 191.908 | 56.426 | 156.932 | 56.429 | 148.631 | -10.955 | -34.976 | -10.951 | -43.278 | 0.617 / 2.559 | 11.797 / 9.303 |
| 4→5 | blyp | 67.381 | 191.908 | 51.454 | 158.408 | 51.441 | 150.023 | -15.926 | -33.500 | -15.939 | -41.885 | 0.616 / 2.557 | 9.048 / 7.181 |

### Temperature-dependent raw and local corrections

Only new n=4 and 5 quadratures change in the 1024/2048 comparison; the reused shorter ensembles retain their original production counts.

| Step | Level | T (K) | GAV ΔH / oriented ΔS | Calibrated raw ΔH / local ΔH | Raw δH (kJ/mol) | Raw δS (J/mol/K) | Local δH | Local δS | MC SE H/S/G (kJ/mol, J/mol/K, kJ/mol) | Sequence SD H/S | Change from 1024 to 2048 points: local H/S |
| --- | --- | ---: | --- | --- | ---: | ---: | ---: | ---: | --- | --- | --- |
| 2→3 | pbe | 298.15 | 67.381 / 191.908 | 57.507 / 57.292 | -9.874 | -32.158 | -10.088 | -46.093 | 0.179 / 0.532 / 0.106 | 4.797 / 5.279 | 0.000 / 0.000 |
| 2→3 | pbe | 300.00 | 67.606 / 192.663 | 57.757 / 57.541 | -9.849 | -32.077 | -10.066 | -46.018 | 0.180 / 0.534 / 0.106 | 4.797 / 5.276 | 0.000 / 0.000 |
| 2→3 | pbe | 400.00 | 81.874 / 233.428 | 73.108 / 72.831 | -8.766 | -28.924 | -9.043 | -43.043 | 0.210 / 0.605 / 0.126 | 4.801 / 5.004 | 0.000 / 0.000 |
| 2→3 | pbe | 500.00 | 99.882 / 273.470 | 91.793 / 91.506 | -8.089 | -27.402 | -8.375 | -41.545 | 0.245 / 0.609 / 0.160 | 4.635 / 4.413 | 0.000 / 0.000 |
| 2→3 | pbe | 600.00 | 120.986 / 311.871 | 113.440 / 113.171 | -7.545 | -26.411 | -7.814 | -40.523 | 0.315 / 0.653 / 0.197 | 4.362 / 3.795 | 0.000 / 0.000 |
| 2→3 | pbe | 700.00 | 144.490 / 348.063 | 137.633 / 137.396 | -6.858 | -25.357 | -7.094 | -39.418 | 0.383 / 0.694 / 0.239 | 4.081 / 3.316 | 0.000 / 0.000 |
| 2→3 | pbe | 800.00 | 170.026 / 382.131 | 163.979 / 163.780 | -6.047 | -24.274 | -6.247 | -38.286 | 0.447 / 0.723 / 0.284 | 3.825 / 2.936 | 0.000 / 0.000 |
| 2→3 | blyp | 298.15 | 67.381 / 191.908 | 51.758 / 51.750 | -15.623 | -33.401 | -15.630 | -46.950 | 0.183 / 0.545 / 0.106 | 4.511 / 5.401 | 0.000 / 0.000 |
| 2→3 | blyp | 300.00 | 67.606 / 192.663 | 52.006 / 51.996 | -15.600 | -33.325 | -15.611 | -46.884 | 0.183 / 0.548 / 0.106 | 4.510 / 5.404 | 0.000 / 0.000 |
| 2→3 | blyp | 400.00 | 81.874 / 233.428 | 67.451 / 67.303 | -14.422 | -29.920 | -14.571 | -43.881 | 0.209 / 0.603 / 0.126 | 4.387 / 5.215 | 0.000 / 0.000 |
| 2→3 | blyp | 500.00 | 99.882 / 273.470 | 86.377 / 86.157 | -13.505 | -27.863 | -13.725 | -41.985 | 0.237 / 0.583 / 0.159 | 4.030 / 4.415 | 0.000 / 0.000 |
| 2→3 | blyp | 600.00 | 120.986 / 311.871 | 108.258 / 108.027 | -12.727 | -26.444 | -12.959 | -40.590 | 0.308 / 0.623 / 0.193 | 3.563 / 3.517 | 0.000 / 0.000 |
| 2→3 | blyp | 700.00 | 144.490 / 348.063 | 132.636 / 132.427 | -11.855 | -25.104 | -12.064 | -39.215 | 0.379 / 0.669 / 0.231 | 3.131 / 2.818 | 0.000 / 0.000 |
| 2→3 | blyp | 800.00 | 170.026 / 382.131 | 159.113 / 158.940 | -10.914 | -23.846 | -11.086 | -37.908 | 0.445 / 0.707 / 0.273 | 2.779 / 2.319 | 0.000 / 0.000 |
| 3→4 | pbe | 298.15 | 67.381 / 191.908 | 65.492 / 65.638 | -1.889 | -32.899 | -1.743 | -42.448 | 0.380 / 1.429 / 0.163 | 8.869 / 6.386 | -0.351 / -0.772 |
| 3→4 | pbe | 300.00 | 67.606 / 192.663 | 65.742 / 65.890 | -1.864 | -32.815 | -1.716 | -42.360 | 0.390 / 1.456 / 0.164 | 8.890 / 6.367 | -0.391 / -0.904 |
| 3→4 | pbe | 400.00 | 81.874 / 233.428 | 82.088 / 82.305 | 0.214 | -26.881 | 0.431 | -36.227 | 0.939 / 2.885 / 0.305 | 10.054 / 6.394 | -1.827 / -5.190 |
| 3→4 | pbe | 500.00 | 99.882 / 273.470 | 101.358 / 101.633 | 1.477 | -23.999 | 1.752 | -33.214 | 0.989 / 2.809 / 0.545 | 10.333 / 5.694 | -1.483 / -4.459 |
| 3→4 | pbe | 600.00 | 120.986 / 311.871 | 122.431 / 122.722 | 1.446 | -24.030 | 1.736 | -33.215 | 1.001 / 2.534 / 0.767 | 10.183 / 4.888 | -0.778 / -3.173 |
| 3→4 | pbe | 700.00 | 144.490 / 348.063 | 145.450 / 145.706 | 0.959 | -24.780 | 1.215 | -34.017 | 0.999 / 2.185 / 0.955 | 10.135 / 5.008 | -0.191 / -2.266 |
| 3→4 | pbe | 800.00 | 170.026 / 382.131 | 170.529 / 170.715 | 0.503 | -25.390 | 0.688 | -34.721 | 1.135 / 1.909 / 1.103 | 10.244 / 5.446 | 0.315 / -1.589 |
| 3→4 | blyp | 298.15 | 67.381 / 191.908 | 57.843 / 57.929 | -9.537 | -31.352 | -9.452 | -41.302 | 0.377 / 1.418 / 0.166 | 6.513 / 7.175 | -0.347 / -0.737 |
| 3→4 | blyp | 300.00 | 67.606 / 192.663 | 58.093 / 58.180 | -9.514 | -31.274 | -9.427 | -41.218 | 0.387 / 1.444 / 0.167 | 6.531 / 7.153 | -0.387 / -0.869 |
| 3→4 | blyp | 400.00 | 81.874 / 233.428 | 74.300 / 74.482 | -7.574 | -25.733 | -7.392 | -35.406 | 0.942 / 2.881 / 0.303 | 7.524 / 6.635 | -1.842 / -5.207 |
| 3→4 | blyp | 500.00 | 99.882 / 273.470 | 93.364 / 93.643 | -6.518 | -23.307 | -6.238 | -32.763 | 0.992 / 2.809 / 0.542 | 7.674 / 5.570 | -1.486 / -4.453 |
| 3→4 | blyp | 600.00 | 120.986 / 311.871 | 114.192 / 114.520 | -6.794 | -23.786 | -6.465 | -33.149 | 1.000 / 2.527 / 0.764 | 7.413 / 4.486 | -0.763 / -3.133 |
| 3→4 | blyp | 700.00 | 144.490 / 348.063 | 136.977 / 137.297 | -7.514 | -24.896 | -7.194 | -34.271 | 0.997 / 2.176 / 0.951 | 7.327 / 4.499 | -0.166 / -2.210 |
| 3→4 | blyp | 800.00 | 170.026 / 382.131 | 161.865 / 162.132 | -8.162 | -25.762 | -7.894 | -35.207 | 1.133 / 1.902 / 1.098 | 7.468 / 4.986 | 0.340 / -1.534 |
| 4→5 | pbe | 298.15 | 67.381 / 191.908 | 56.426 / 56.429 | -10.955 | -34.976 | -10.951 | -43.278 | 0.617 / 2.559 / 0.281 | 11.797 / 9.303 | 1.047 / 5.088 |
| 4→5 | pbe | 300.00 | 67.606 / 192.663 | 56.644 / 56.646 | -10.963 | -35.001 | -10.960 | -43.307 | 0.629 / 2.585 / 0.285 | 11.810 / 9.201 | 1.085 / 5.213 |
| 4→5 | pbe | 400.00 | 81.874 / 233.428 | 71.405 / 71.345 | -10.469 | -33.723 | -10.529 | -42.210 | 1.283 / 3.929 / 0.500 | 12.410 / 7.076 | 2.390 / 9.099 |
| 4→5 | pbe | 500.00 | 99.882 / 273.470 | 90.798 / 90.783 | -9.084 | -30.606 | -9.099 | -38.998 | 1.467 / 3.875 / 0.786 | 12.242 / 10.325 | 2.983 / 10.433 |
| 4→5 | pbe | 600.00 | 120.986 / 311.871 | 112.663 / 112.745 | -8.323 | -29.217 | -8.241 | -37.433 | 1.603 / 3.482 / 1.043 | 11.347 / 10.491 | 2.615 / 9.793 |
| 4→5 | pbe | 700.00 | 144.490 / 348.063 | 137.830 / 137.981 | -6.661 | -26.688 | -6.510 | -34.796 | 2.078 / 3.578 / 1.256 | 10.376 / 8.839 | 1.653 / 8.312 |
| 4→5 | pbe | 800.00 | 170.026 / 382.131 | 166.194 / 166.360 | -3.832 | -22.916 | -3.666 | -31.004 | 2.705 / 3.976 / 1.459 | 11.114 / 9.452 | 0.707 / 7.050 |
| 4→5 | blyp | 298.15 | 67.381 / 191.908 | 51.454 / 51.441 | -15.926 | -33.500 | -15.939 | -41.885 | 0.616 / 2.557 / 0.281 | 9.048 / 7.181 | 1.081 / 5.114 |
| 4→5 | blyp | 300.00 | 67.606 / 192.663 | 51.678 / 51.663 | -15.929 | -33.508 | -15.943 | -41.897 | 0.628 / 2.583 / 0.284 | 9.057 / 7.048 | 1.121 / 5.248 |
| 4→5 | blyp | 400.00 | 81.874 / 233.428 | 66.429 / 66.383 | -15.445 | -32.213 | -15.491 | -40.699 | 1.289 / 3.948 / 0.500 | 9.835 / 5.823 | 2.533 / 9.456 |
| 4→5 | blyp | 500.00 | 99.882 / 273.470 | 85.446 / 85.464 | -14.435 | -29.935 | -14.418 | -38.281 | 1.462 / 3.868 / 0.788 | 10.004 / 10.033 | 3.038 / 10.596 |
| 4→5 | blyp | 600.00 | 120.986 / 311.871 | 107.098 / 107.206 | -13.888 | -28.941 | -13.780 | -37.122 | 1.607 / 3.482 / 1.044 | 9.230 / 10.368 | 2.622 / 9.868 |
| 4→5 | blyp | 700.00 | 144.490 / 348.063 | 132.192 / 132.363 | -12.299 | -26.525 | -12.128 | -34.609 | 2.096 / 3.599 / 1.255 | 8.219 / 8.882 | 1.634 / 8.346 |
| 4→5 | blyp | 800.00 | 170.026 / 382.131 | 160.513 / 160.693 | -9.513 | -22.811 | -9.334 | -30.883 | 2.709 / 3.985 / 1.460 | 9.155 / 9.641 | 0.640 / 7.021 |

### Entropy component increments of the chains

These are chain increments before the balanced-reference subtraction. Local partitions remove external factors before reweighting conformer basins; their internal components need not equal the gas-weight component subtraction. No class-mixing bit is added.

| Step | Level | T (K) | Partition | Translation | External rotation | Non-torsional vibration | Coupled rotors | Basin mixing | Sum (J/mol/K) |
| --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2→3 | pbe | 298.15 | raw | 4.756 | 11.353 | 72.566 | 68.622 | 5.336 | 162.632 |
| 2→3 | pbe | 298.15 | local | 0.000 | 0.000 | 72.535 | 68.325 | 4.955 | 145.815 |
| 2→3 | pbe | 300.00 | raw | 4.756 | 11.357 | 73.172 | 68.829 | 5.353 | 163.467 |
| 2→3 | pbe | 300.00 | local | 0.000 | 0.000 | 73.140 | 68.528 | 4.976 | 146.645 |
| 2→3 | pbe | 400.00 | raw | 4.756 | 11.572 | 107.236 | 78.345 | 5.477 | 207.386 |
| 2→3 | pbe | 400.00 | local | 0.000 | 0.000 | 107.182 | 77.861 | 5.342 | 190.385 |
| 2→3 | pbe | 500.00 | raw | 4.756 | 11.714 | 142.357 | 85.178 | 4.945 | 248.950 |
| 2→3 | pbe | 500.00 | local | 0.000 | 0.000 | 142.293 | 84.601 | 5.031 | 231.925 |
| 2→3 | pbe | 600.00 | raw | 4.756 | 11.803 | 176.956 | 90.341 | 4.484 | 288.341 |
| 2→3 | pbe | 600.00 | local | 0.000 | 0.000 | 176.889 | 89.754 | 4.706 | 271.348 |
| 2→3 | pbe | 700.00 | raw | 4.756 | 11.859 | 210.239 | 94.489 | 4.245 | 325.588 |
| 2→3 | pbe | 700.00 | local | 0.000 | 0.000 | 210.172 | 93.942 | 4.531 | 308.645 |
| 2→3 | pbe | 800.00 | raw | 4.756 | 11.893 | 241.908 | 98.003 | 4.179 | 360.739 |
| 2→3 | pbe | 800.00 | local | 0.000 | 0.000 | 241.844 | 97.518 | 4.483 | 343.845 |
| 2→3 | blyp | 298.15 | raw | 4.756 | 11.594 | 72.666 | 68.101 | 4.273 | 161.389 |
| 2→3 | blyp | 298.15 | local | 0.000 | 0.000 | 72.729 | 68.264 | 3.965 | 144.958 |
| 2→3 | blyp | 300.00 | raw | 4.756 | 11.594 | 73.268 | 68.288 | 4.313 | 162.220 |
| 2→3 | blyp | 300.00 | local | 0.000 | 0.000 | 73.331 | 68.446 | 4.002 | 145.779 |
| 2→3 | blyp | 400.00 | raw | 4.756 | 11.676 | 107.200 | 77.184 | 5.575 | 206.390 |
| 2→3 | blyp | 400.00 | local | 0.000 | 0.000 | 107.236 | 77.057 | 5.254 | 189.547 |
| 2→3 | blyp | 500.00 | raw | 4.756 | 11.779 | 142.269 | 83.955 | 5.730 | 248.489 |
| 2→3 | blyp | 500.00 | local | 0.000 | 0.000 | 142.284 | 83.635 | 5.566 | 231.485 |
| 2→3 | blyp | 600.00 | raw | 4.756 | 11.859 | 176.850 | 89.303 | 5.540 | 288.308 |
| 2→3 | blyp | 600.00 | local | 0.000 | 0.000 | 176.851 | 88.902 | 5.528 | 271.281 |
| 2→3 | blyp | 700.00 | raw | 4.756 | 11.912 | 210.129 | 93.715 | 5.329 | 325.841 |
| 2→3 | blyp | 700.00 | local | 0.000 | 0.000 | 210.121 | 93.310 | 5.417 | 308.848 |
| 2→3 | blyp | 800.00 | raw | 4.756 | 11.945 | 241.800 | 97.489 | 5.178 | 361.167 |
| 2→3 | blyp | 800.00 | local | 0.000 | 0.000 | 241.786 | 97.120 | 5.317 | 344.223 |
| 3→4 | pbe | 298.15 | raw | 3.435 | 5.138 | 72.913 | 72.765 | 3.318 | 157.569 |
| 3→4 | pbe | 298.15 | local | 0.000 | 0.000 | 72.976 | 73.068 | 3.417 | 149.460 |
| 3→4 | pbe | 300.00 | raw | 3.435 | 5.134 | 73.515 | 73.040 | 3.283 | 158.407 |
| 3→4 | pbe | 300.00 | local | 0.000 | 0.000 | 73.579 | 73.345 | 3.379 | 150.303 |
| 3→4 | pbe | 400.00 | raw | 3.435 | 4.968 | 107.308 | 88.361 | 1.035 | 205.106 |
| 3→4 | pbe | 400.00 | local | 0.000 | 0.000 | 107.412 | 88.800 | 0.989 | 197.201 |
| 3→4 | pbe | 500.00 | raw | 3.435 | 4.834 | 142.169 | 99.008 | -1.415 | 248.031 |
| 3→4 | pbe | 500.00 | local | 0.000 | 0.000 | 142.304 | 99.598 | -1.645 | 240.256 |
| 3→4 | pbe | 600.00 | raw | 3.435 | 4.732 | 176.574 | 104.975 | -3.316 | 286.400 |
| 3→4 | pbe | 600.00 | local | 0.000 | 0.000 | 176.732 | 105.637 | -3.713 | 278.656 |
| 3→4 | pbe | 700.00 | raw | 3.435 | 4.663 | 209.735 | 108.647 | -4.638 | 321.842 |
| 3→4 | pbe | 700.00 | local | 0.000 | 0.000 | 209.910 | 109.294 | -5.158 | 314.046 |
| 3→4 | pbe | 800.00 | raw | 3.435 | 4.621 | 241.344 | 111.410 | -5.509 | 355.301 |
| 3→4 | pbe | 800.00 | local | 0.000 | 0.000 | 241.531 | 111.985 | -6.106 | 347.410 |
| 3→4 | blyp | 298.15 | raw | 3.435 | 5.348 | 73.137 | 74.351 | 2.845 | 159.115 |
| 3→4 | blyp | 298.15 | local | 0.000 | 0.000 | 73.172 | 74.462 | 2.972 | 150.606 |
| 3→4 | blyp | 300.00 | raw | 3.435 | 5.346 | 73.742 | 74.632 | 2.795 | 159.949 |
| 3→4 | blyp | 300.00 | local | 0.000 | 0.000 | 73.777 | 74.745 | 2.923 | 151.445 |
| 3→4 | blyp | 400.00 | raw | 3.435 | 5.224 | 107.629 | 90.134 | -0.168 | 206.254 |
| 3→4 | blyp | 400.00 | local | 0.000 | 0.000 | 107.700 | 90.378 | -0.055 | 198.023 |
| 3→4 | blyp | 500.00 | raw | 3.435 | 5.096 | 142.531 | 100.725 | -3.065 | 248.722 |
| 3→4 | blyp | 500.00 | local | 0.000 | 0.000 | 142.639 | 101.150 | -3.081 | 240.708 |
| 3→4 | blyp | 600.00 | raw | 3.435 | 4.984 | 176.945 | 106.415 | -5.136 | 286.644 |
| 3→4 | blyp | 600.00 | local | 0.000 | 0.000 | 177.084 | 106.943 | -5.305 | 278.722 |
| 3→4 | blyp | 700.00 | raw | 3.435 | 4.900 | 210.102 | 109.732 | -6.442 | 321.726 |
| 3→4 | blyp | 700.00 | local | 0.000 | 0.000 | 210.264 | 110.263 | -6.735 | 313.792 |
| 3→4 | blyp | 800.00 | raw | 3.435 | 4.843 | 241.699 | 112.151 | -7.199 | 354.928 |
| 3→4 | blyp | 800.00 | local | 0.000 | 0.000 | 241.879 | 112.618 | -7.573 | 346.924 |
| 4→5 | pbe | 298.15 | raw | 2.690 | 7.068 | 74.239 | 76.444 | -2.068 | 158.373 |
| 4→5 | pbe | 298.15 | local | 0.000 | 0.000 | 74.234 | 76.434 | -2.037 | 148.631 |
| 4→5 | pbe | 300.00 | raw | 2.690 | 7.068 | 74.845 | 76.570 | -2.070 | 159.103 |
| 4→5 | pbe | 300.00 | local | 0.000 | 0.000 | 74.839 | 76.558 | -2.040 | 149.356 |
| 4→5 | pbe | 400.00 | raw | 2.690 | 7.092 | 108.920 | 84.051 | -1.607 | 201.146 |
| 4→5 | pbe | 400.00 | local | 0.000 | 0.000 | 108.908 | 83.903 | -1.593 | 191.218 |
| 4→5 | pbe | 500.00 | raw | 2.690 | 7.117 | 144.036 | 92.048 | -1.586 | 244.305 |
| 4→5 | pbe | 500.00 | local | 0.000 | 0.000 | 144.024 | 91.978 | -1.530 | 234.472 |
| 4→5 | pbe | 600.00 | raw | 2.690 | 7.108 | 178.562 | 97.248 | -1.514 | 284.094 |
| 4→5 | pbe | 600.00 | local | 0.000 | 0.000 | 178.550 | 97.347 | -1.460 | 274.438 |
| 4→5 | pbe | 700.00 | raw | 2.690 | 7.080 | 211.758 | 102.480 | -1.192 | 322.816 |
| 4→5 | pbe | 700.00 | local | 0.000 | 0.000 | 211.743 | 102.688 | -1.164 | 313.268 |
| 4→5 | pbe | 800.00 | raw | 2.690 | 7.052 | 243.355 | 108.814 | -1.255 | 360.656 |
| 4→5 | pbe | 800.00 | local | 0.000 | 0.000 | 243.338 | 109.037 | -1.247 | 351.128 |
| 4→5 | blyp | 298.15 | raw | 2.690 | 7.095 | 74.293 | 76.736 | -0.965 | 159.849 |
| 4→5 | blyp | 298.15 | local | 0.000 | 0.000 | 74.285 | 76.697 | -0.958 | 150.023 |
| 4→5 | blyp | 300.00 | raw | 2.690 | 7.095 | 74.899 | 76.876 | -0.965 | 160.596 |
| 4→5 | blyp | 300.00 | local | 0.000 | 0.000 | 74.891 | 76.833 | -0.957 | 150.766 |
| 4→5 | blyp | 400.00 | raw | 2.690 | 7.125 | 108.982 | 84.800 | -0.941 | 202.656 |
| 4→5 | blyp | 400.00 | local | 0.000 | 0.000 | 108.970 | 84.633 | -0.874 | 192.730 |
| 4→5 | blyp | 500.00 | raw | 2.690 | 7.137 | 144.073 | 92.500 | -1.423 | 244.976 |
| 4→5 | blyp | 500.00 | local | 0.000 | 0.000 | 144.063 | 92.452 | -1.326 | 235.189 |
| 4→5 | blyp | 600.00 | raw | 2.690 | 7.117 | 178.582 | 97.487 | -1.506 | 284.370 |
| 4→5 | blyp | 600.00 | local | 0.000 | 0.000 | 178.571 | 97.611 | -1.434 | 274.748 |
| 4→5 | blyp | 700.00 | raw | 2.690 | 7.084 | 211.768 | 102.665 | -1.228 | 322.979 |
| 4→5 | blyp | 700.00 | local | 0.000 | 0.000 | 211.754 | 102.894 | -1.195 | 313.454 |
| 4→5 | blyp | 800.00 | raw | 2.690 | 7.053 | 243.359 | 108.965 | -1.306 | 360.761 |
| 4→5 | blyp | 800.00 | local | 0.000 | 0.000 | 243.343 | 109.205 | -1.299 | 351.249 |

### Enthalpy component increments of the chains

The electronic column includes the constant formation-enthalpy reference determined by the balanced cycle at 298 K. The uncalibrated electronic increment and all balanced-cycle components are also preserved in JSON. Translation and rotation each have zero H increment for one nonlinear chain growing into another. Basin mixing contributes entropy rather than a separate enthalpy; electronic ensemble averaging includes population energy shifts.

| Step | Level | T (K) | Partition | Electronic + formation anchor | Translation | External rotation | Non-torsional vibration | Coupled rotors | Basin mixing | Sum (kJ/mol) |
| --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2→3 | pbe | 298.15 | raw | -315.843 | 0.000 | 0.000 | 365.605 | 7.745 | 0.000 | 57.507 |
| 2→3 | pbe | 298.15 | local | -316.036 | 0.000 | 0.000 | 365.609 | 7.720 | 0.000 | 57.292 |
| 2→3 | pbe | 300.00 | raw | -315.827 | -0.000 | -0.000 | 365.786 | 7.799 | 0.000 | 57.757 |
| 2→3 | pbe | 300.00 | local | -316.022 | 0.000 | 0.000 | 365.789 | 7.773 | 0.000 | 57.541 |
| 2→3 | pbe | 400.00 | raw | -315.235 | 0.000 | 0.000 | 377.657 | 10.686 | 0.000 | 73.108 |
| 2→3 | pbe | 400.00 | local | -315.463 | 0.000 | 0.000 | 377.661 | 10.633 | 0.000 | 72.831 |
| 2→3 | pbe | 500.00 | raw | -315.004 | 0.000 | 0.000 | 393.396 | 13.401 | 0.000 | 91.793 |
| 2→3 | pbe | 500.00 | local | -315.213 | 0.000 | 0.000 | 393.400 | 13.319 | 0.000 | 91.506 |
| 2→3 | pbe | 600.00 | raw | -314.906 | -0.000 | -0.000 | 412.359 | 15.988 | 0.000 | 113.440 |
| 2→3 | pbe | 600.00 | local | -315.083 | 0.000 | 0.000 | 412.363 | 15.892 | 0.000 | 113.171 |
| 2→3 | pbe | 700.00 | raw | -314.819 | 0.000 | 0.000 | 433.933 | 18.519 | 0.000 | 137.633 |
| 2→3 | pbe | 700.00 | local | -314.968 | 0.000 | 0.000 | 433.937 | 18.427 | 0.000 | 137.396 |
| 2→3 | pbe | 800.00 | raw | -314.706 | 0.000 | 0.000 | 457.637 | 21.049 | 0.000 | 163.979 |
| 2→3 | pbe | 800.00 | local | -314.833 | 0.000 | 0.000 | 457.641 | 20.972 | 0.000 | 163.780 |
| 2→3 | blyp | 298.15 | raw | -321.520 | 0.000 | 0.000 | 365.621 | 7.657 | 0.000 | 51.758 |
| 2→3 | blyp | 298.15 | local | -321.507 | 0.000 | 0.000 | 365.625 | 7.632 | 0.000 | 51.750 |
| 2→3 | blyp | 300.00 | raw | -321.504 | -0.000 | -0.000 | 365.801 | 7.710 | 0.000 | 52.006 |
| 2→3 | blyp | 300.00 | local | -321.493 | 0.000 | 0.000 | 365.805 | 7.684 | 0.000 | 51.996 |
| 2→3 | blyp | 400.00 | raw | -320.745 | 0.000 | 0.000 | 377.667 | 10.530 | 0.000 | 67.451 |
| 2→3 | blyp | 400.00 | local | -320.832 | 0.000 | 0.000 | 377.672 | 10.463 | 0.000 | 67.303 |
| 2→3 | blyp | 500.00 | raw | -320.238 | 0.000 | 0.000 | 393.402 | 13.213 | 0.000 | 86.377 |
| 2→3 | blyp | 500.00 | local | -320.354 | 0.000 | 0.000 | 393.407 | 13.104 | 0.000 | 86.157 |
| 2→3 | blyp | 600.00 | raw | -319.921 | -0.000 | -0.000 | 412.362 | 15.817 | 0.000 | 108.258 |
| 2→3 | blyp | 600.00 | local | -320.023 | 0.000 | 0.000 | 412.367 | 15.683 | 0.000 | 108.027 |
| 2→3 | blyp | 700.00 | raw | -319.701 | 0.000 | 0.000 | 433.934 | 18.403 | 0.000 | 132.636 |
| 2→3 | blyp | 700.00 | local | -319.776 | 0.000 | 0.000 | 433.939 | 18.263 | 0.000 | 132.427 |
| 2→3 | blyp | 800.00 | raw | -319.526 | 0.000 | 0.000 | 457.637 | 21.002 | 0.000 | 159.113 |
| 2→3 | blyp | 800.00 | local | -319.572 | 0.000 | 0.000 | 457.641 | 20.871 | 0.000 | 158.940 |
| 3→4 | pbe | 298.15 | raw | -308.511 | 0.000 | 0.000 | 365.459 | 8.544 | 0.000 | 65.492 |
| 3→4 | pbe | 298.15 | local | -308.366 | 0.000 | 0.000 | 365.456 | 8.548 | 0.000 | 65.638 |
| 3→4 | pbe | 300.00 | raw | -308.525 | -0.000 | -0.000 | 365.639 | 8.628 | 0.000 | 65.742 |
| 3→4 | pbe | 300.00 | local | -308.378 | 0.000 | 0.000 | 365.637 | 8.631 | 0.000 | 65.890 |
| 3→4 | pbe | 400.00 | raw | -309.334 | 0.000 | 0.000 | 377.541 | 13.881 | 0.000 | 82.088 |
| 3→4 | pbe | 400.00 | local | -309.106 | 0.000 | 0.000 | 377.537 | 13.874 | 0.000 | 82.305 |
| 3→4 | pbe | 500.00 | raw | -310.076 | -0.000 | 0.000 | 393.295 | 18.140 | 0.000 | 101.358 |
| 3→4 | pbe | 500.00 | local | -309.807 | 0.000 | 0.000 | 393.290 | 18.150 | 0.000 | 101.633 |
| 3→4 | pbe | 600.00 | raw | -310.713 | -0.000 | -0.000 | 412.268 | 20.877 | 0.000 | 122.431 |
| 3→4 | pbe | 600.00 | local | -310.436 | 0.000 | 0.000 | 412.263 | 20.895 | 0.000 | 122.722 |
| 3→4 | pbe | 700.00 | raw | -311.283 | 0.000 | 0.000 | 433.851 | 22.881 | 0.000 | 145.450 |
| 3→4 | pbe | 700.00 | local | -311.017 | 0.000 | 0.000 | 433.847 | 22.876 | 0.000 | 145.706 |
| 3→4 | pbe | 800.00 | raw | -311.797 | 0.000 | 0.000 | 457.563 | 24.764 | 0.000 | 170.529 |
| 3→4 | pbe | 800.00 | local | -311.552 | 0.000 | 0.000 | 457.559 | 24.708 | 0.000 | 170.715 |
| 3→4 | blyp | 298.15 | raw | -316.174 | 0.000 | 0.000 | 365.449 | 8.569 | 0.000 | 57.843 |
| 3→4 | blyp | 298.15 | local | -316.085 | 0.000 | 0.000 | 365.448 | 8.566 | 0.000 | 57.929 |
| 3→4 | blyp | 300.00 | raw | -316.189 | -0.000 | -0.000 | 365.630 | 8.652 | 0.000 | 58.093 |
| 3→4 | blyp | 300.00 | local | -316.098 | 0.000 | 0.000 | 365.629 | 8.648 | 0.000 | 58.180 |
| 3→4 | blyp | 400.00 | raw | -317.133 | 0.000 | 0.000 | 377.533 | 13.899 | 0.000 | 74.300 |
| 3→4 | blyp | 400.00 | local | -316.935 | 0.000 | 0.000 | 377.532 | 13.886 | 0.000 | 74.482 |
| 3→4 | blyp | 500.00 | raw | -318.087 | -0.000 | 0.000 | 393.290 | 18.161 | 0.000 | 93.364 |
| 3→4 | blyp | 500.00 | local | -317.810 | 0.000 | 0.000 | 393.288 | 18.166 | 0.000 | 93.643 |
| 3→4 | blyp | 600.00 | raw | -318.922 | -0.000 | -0.000 | 412.266 | 20.848 | 0.000 | 114.192 |
| 3→4 | blyp | 600.00 | local | -318.603 | 0.000 | 0.000 | 412.263 | 20.860 | 0.000 | 114.520 |
| 3→4 | blyp | 700.00 | raw | -319.635 | 0.000 | 0.000 | 433.851 | 22.760 | 0.000 | 136.977 |
| 3→4 | blyp | 700.00 | local | -319.300 | 0.000 | 0.000 | 433.849 | 22.748 | 0.000 | 137.297 |
| 3→4 | blyp | 800.00 | raw | -320.234 | 0.000 | 0.000 | 457.564 | 24.535 | 0.000 | 161.865 |
| 3→4 | blyp | 800.00 | local | -319.899 | 0.000 | 0.000 | 457.562 | 24.470 | 0.000 | 162.132 |
| 4→5 | pbe | 298.15 | raw | -317.045 | 0.000 | 0.000 | 365.571 | 7.900 | 0.000 | 56.426 |
| 4→5 | pbe | 298.15 | local | -317.031 | 0.000 | 0.000 | 365.575 | 7.885 | 0.000 | 56.429 |
| 4→5 | pbe | 300.00 | raw | -317.042 | -0.000 | -0.000 | 365.751 | 7.934 | 0.000 | 56.644 |
| 4→5 | pbe | 300.00 | local | -317.028 | 0.000 | 0.000 | 365.756 | 7.919 | 0.000 | 56.646 |
| 4→5 | pbe | 400.00 | raw | -316.513 | 0.000 | 0.000 | 377.642 | 10.276 | 0.000 | 71.405 |
| 4→5 | pbe | 400.00 | local | -316.536 | 0.000 | 0.000 | 377.646 | 10.235 | 0.000 | 71.345 |
| 4→5 | pbe | 500.00 | raw | -315.950 | 0.000 | -0.000 | 393.398 | 13.350 | 0.000 | 90.798 |
| 4→5 | pbe | 500.00 | local | -315.970 | 0.000 | 0.000 | 393.402 | 13.351 | 0.000 | 90.783 |
| 4→5 | pbe | 600.00 | raw | -315.748 | -0.000 | -0.000 | 412.372 | 16.039 | 0.000 | 112.663 |
| 4→5 | pbe | 600.00 | local | -315.753 | 0.000 | 0.000 | 412.376 | 16.122 | 0.000 | 112.745 |
| 4→5 | pbe | 700.00 | raw | -315.727 | 0.000 | 0.000 | 433.956 | 19.600 | 0.000 | 137.830 |
| 4→5 | pbe | 700.00 | local | -315.722 | 0.000 | 0.000 | 433.960 | 19.743 | 0.000 | 137.981 |
| 4→5 | pbe | 800.00 | raw | -315.778 | 0.000 | 0.000 | 457.669 | 24.303 | 0.000 | 166.194 |
| 4→5 | pbe | 800.00 | local | -315.767 | 0.000 | 0.000 | 457.673 | 24.454 | 0.000 | 166.360 |
| 4→5 | blyp | 298.15 | raw | -322.065 | 0.000 | 0.000 | 365.566 | 7.954 | 0.000 | 51.454 |
| 4→5 | blyp | 298.15 | local | -322.063 | 0.000 | 0.000 | 365.571 | 7.934 | 0.000 | 51.441 |
| 4→5 | blyp | 300.00 | raw | -322.061 | -0.000 | -0.000 | 365.746 | 7.992 | 0.000 | 51.678 |
| 4→5 | blyp | 300.00 | local | -322.059 | 0.000 | 0.000 | 365.751 | 7.971 | 0.000 | 51.663 |
| 4→5 | blyp | 400.00 | raw | -321.649 | 0.000 | 0.000 | 377.637 | 10.440 | 0.000 | 66.429 |
| 4→5 | blyp | 400.00 | local | -321.656 | 0.000 | 0.000 | 377.641 | 10.397 | 0.000 | 66.383 |
| 4→5 | blyp | 500.00 | raw | -321.374 | 0.000 | -0.000 | 393.396 | 13.425 | 0.000 | 85.446 |
| 4→5 | blyp | 500.00 | local | -321.369 | 0.000 | 0.000 | 393.400 | 13.433 | 0.000 | 85.464 |
| 4→5 | blyp | 600.00 | raw | -321.342 | -0.000 | -0.000 | 412.372 | 16.068 | 0.000 | 107.098 |
| 4→5 | blyp | 600.00 | local | -321.329 | 0.000 | 0.000 | 412.375 | 16.159 | 0.000 | 107.206 |
| 4→5 | blyp | 700.00 | raw | -321.415 | 0.000 | 0.000 | 433.957 | 19.650 | 0.000 | 132.192 |
| 4→5 | blyp | 700.00 | local | -321.400 | 0.000 | 0.000 | 433.961 | 19.801 | 0.000 | 132.363 |
| 4→5 | blyp | 800.00 | raw | -321.527 | 0.000 | 0.000 | 457.670 | 24.370 | 0.000 | 160.513 |
| 4→5 | blyp | 800.00 | local | -321.510 | 0.000 | 0.000 | 457.674 | 24.529 | 0.000 | 160.693 |

### Balanced-reference diagnostics and finite-fragment conditional ceilings

| Step | Level | Chain translation at 298 | Chain rotation at 298 | Reference-side translation | Reference-side rotation | Balanced raw δS | Balanced local δS | Component-subtracted balanced δS | All-cycle internal δS | Raw conditional Tc (K) | Local conditional Tc (K) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| 2→3 | pbe | 4.756 | 11.353 | -18.828 | -56.650 | -30.161 | -44.096 | -43.388 | 31.419 | 623.293 ± 1.606 (MC) | 564.542 ± 1.256 (MC); channel SD 22.492 K (linearized) |
| 2→3 | blyp | 4.756 | 11.594 | -18.828 | -56.628 | -31.307 | -44.856 | -44.775 | 30.658 | 663.244 ± 1.684 (MC) | 599.615 ± 1.331 (MC); channel SD 16.420 K (linearized) |
| 3→4 | pbe | 3.435 | 5.138 | -18.828 | -56.650 | -30.902 | -40.451 | -40.915 | 35.064 | 564.858 ± 5.382 (MC) | 525.383 ± 4.337 (MC); channel SD 64.383 K (linearized) |
| 3→4 | blyp | 3.435 | 5.348 | -18.828 | -56.628 | -29.259 | -39.208 | -39.482 | 36.306 | 630.171 ± 6.484 (MC) | 584.853 ± 5.311 (MC); channel SD 44.930 K (linearized) |
| 4→5 | pbe | 2.690 | 7.068 | -18.828 | -56.650 | -32.979 | -41.281 | -41.297 | 34.234 | 615.992 ± 8.160 (MC) | 579.754 ± 6.974 (MC); channel SD 77.275 K (linearized) |
| 4→5 | blyp | 2.690 | 7.095 | -18.828 | -56.628 | -31.406 | -39.791 | -39.751 | 35.723 | 659.521 ± 9.024 (MC) | 620.165 ± 7.748 (MC); channel SD 63.528 K (linearized) |

### Functional spread and length dependence

| Step | T (K) | PBE−BLYP local δH (kJ/mol) | PBE−BLYP local δS (J/mol/K) |
| --- | ---: | ---: | ---: |
| 2→3 | 298.15 | 5.542 | 0.857 |
| 2→3 | 300.00 | 5.545 | 0.866 |
| 2→3 | 400.00 | 5.529 | 0.838 |
| 2→3 | 500.00 | 5.349 | 0.440 |
| 2→3 | 600.00 | 5.145 | 0.067 |
| 2→3 | 700.00 | 4.970 | -0.203 |
| 2→3 | 800.00 | 4.839 | -0.378 |
| 3→4 | 298.15 | 7.709 | -1.146 |
| 3→4 | 300.00 | 7.710 | -1.142 |
| 3→4 | 400.00 | 7.823 | -0.822 |
| 3→4 | 500.00 | 7.990 | -0.452 |
| 3→4 | 600.00 | 8.202 | -0.066 |
| 3→4 | 700.00 | 8.409 | 0.254 |
| 3→4 | 800.00 | 8.582 | 0.486 |
| 4→5 | 298.15 | 4.988 | -1.392 |
| 4→5 | 300.00 | 4.983 | -1.410 |
| 4→5 | 400.00 | 4.962 | -1.512 |
| 4→5 | 500.00 | 5.319 | -0.717 |
| 4→5 | 600.00 | 5.539 | -0.311 |
| 4→5 | 700.00 | 5.618 | -0.187 |
| 4→5 | 800.00 | 5.667 | -0.121 |

| Level | T (K) | Terminal change: (4→5)−(3→4), local δH (kJ/mol) | Terminal local δS change (J/mol/K) | Covariance-aware MC SE H/S | Paired sequence population SD H/S |
| --- | ---: | ---: | ---: | --- | --- |
| pbe | 298.15 | -9.209 | -0.829 | 0.865 / 3.494 | 16.573 / 12.697 |
| pbe | 300.00 | -9.244 | -0.947 | 0.886 / 3.546 | 16.603 / 12.568 |
| pbe | 400.00 | -10.960 | -5.983 | 2.051 / 6.309 | 18.577 / 10.653 |
| pbe | 500.00 | -10.851 | -5.784 | 2.232 / 6.174 | 18.847 / 12.610 |
| pbe | 600.00 | -9.977 | -4.218 | 2.326 / 5.547 | 17.770 / 11.180 |
| pbe | 700.00 | -7.725 | -0.779 | 2.659 / 5.141 | 17.032 / 9.227 |
| pbe | 800.00 | -4.355 | 3.717 | 3.294 / 5.096 | 17.822 / 10.385 |
| blyp | 298.15 | -6.488 | -0.583 | 0.858 / 3.474 | 12.889 / 11.406 |
| blyp | 300.00 | -6.516 | -0.679 | 0.879 / 3.525 | 12.915 / 11.256 |
| blyp | 400.00 | -8.099 | -5.293 | 2.059 / 6.316 | 14.769 / 9.821 |
| blyp | 500.00 | -8.179 | -5.518 | 2.235 / 6.176 | 14.990 / 12.287 |
| blyp | 600.00 | -7.315 | -3.973 | 2.329 / 5.545 | 13.795 / 11.012 |
| blyp | 700.00 | -4.934 | -0.338 | 2.673 / 5.152 | 12.926 / 9.196 |
| blyp | 800.00 | -1.440 | 4.325 | 3.296 / 5.100 | 13.847 / 10.397 |
| gfn2 | 298.15 | 0.192 | 2.826 | 0.868 / 3.501 | 1.837 / 6.450 |
| gfn2 | 300.00 | 0.176 | 2.772 | 0.891 / 3.555 | 1.838 / 6.368 |
| gfn2 | 400.00 | -1.200 | -1.192 | 2.057 / 6.325 | 4.720 / 10.946 |
| gfn2 | 500.00 | -1.532 | -1.972 | 2.222 / 6.157 | 5.936 / 13.239 |
| gfn2 | 600.00 | -0.968 | -0.975 | 2.330 / 5.558 | 5.356 / 11.505 |
| gfn2 | 700.00 | 1.136 | 2.233 | 2.680 / 5.182 | 4.844 / 9.562 |
| gfn2 | 800.00 | 4.364 | 6.540 | 3.290 / 5.118 | 7.156 / 10.589 |
| pbe_sparse | 298.15 | -2.539 | -0.075 | 0.861 / 3.483 | 16.167 / 9.741 |
| pbe_sparse | 300.00 | -2.559 | -0.143 | 0.882 / 3.535 | 16.197 / 9.631 |
| pbe_sparse | 400.00 | -3.629 | -3.299 | 2.049 / 6.303 | 18.089 / 10.661 |
| pbe_sparse | 500.00 | -3.220 | -2.416 | 2.231 / 6.172 | 18.436 / 13.468 |
| pbe_sparse | 600.00 | -2.315 | -0.787 | 2.325 / 5.547 | 17.426 / 12.130 |
| pbe_sparse | 700.00 | -0.170 | 2.488 | 2.659 / 5.143 | 16.528 / 10.048 |
| pbe_sparse | 800.00 | 3.021 | 6.746 | 3.294 / 5.100 | 17.026 / 10.722 |
| blyp_sparse | 298.15 | -2.120 | 1.366 | 0.853 / 3.461 | 12.252 / 7.666 |
| blyp_sparse | 300.00 | -2.134 | 1.317 | 0.874 / 3.511 | 12.282 / 7.532 |
| blyp_sparse | 400.00 | -3.173 | -1.702 | 2.057 / 6.310 | 14.382 / 10.054 |
| blyp_sparse | 500.00 | -3.154 | -1.689 | 2.235 / 6.177 | 14.841 / 13.314 |
| blyp_sparse | 600.00 | -2.494 | -0.512 | 2.329 / 5.548 | 13.759 / 12.028 |
| blyp_sparse | 700.00 | -0.441 | 2.621 | 2.673 / 5.157 | 12.759 / 10.061 |
| blyp_sparse | 800.00 | 2.701 | 6.812 | 3.297 / 5.107 | 13.410 / 10.847 |

### Convergence verdict and transfer

The terminal length changes are listed with covariance-aware Monte Carlo errors and paired sequence spread. Three finite increments, sparse electronic re-ranking, candidate truncation and one doubled quadrature do not establish an asymptotic plateau with a bounded total error. No validated long-chain correction or Tc is adopted; the terminal 4→5 result and its conditional roots remain finite-fragment estimates.
At 298.15 K, comparing (3→4) with (2→3): PBE: production length change 8.346 kJ/mol and 3.645 J/mol/K; uniform sparse diagnostic -4.347 kJ/mol and -0.560 J/mol/K; BLYP: production length change 6.179 kJ/mol and 5.648 J/mol/K; uniform sparse diagnostic -2.450 kJ/mol and 0.017 J/mol/K. The enthalpy trend changes sign under the already-declared energy-treatment diagnostic. The length trend therefore cannot be interpreted independently of that approximation boundary; neither diagnostic replaces the production result.
The internal chain increment is the candidate for transfer because whole-chain translation and external rotation do not accompany local growth of a macroscopic chain. Removing these factors before basin reweighting differs from subtracting gas-weighted component totals. The reported correction also contains intrarepeat QM-versus-GAV differences; it is not an isolated adjacent-phenyl pair interaction.
No n=6 calculation was allocated: required n≤5 searches, both frozen quadratures, composite energies and verification had priority within the original 48-hour cap. The owner-authorized extension explicitly prohibits n=6. This budget decision uses measured computational cost and outstanding required work, without inspecting a ceiling temperature.

### Electronic approximation sensitivity

The pure-GFN2 comparison changes electronic energies only, retaining the same geometries, Hessians and quadratures. It is a sensitivity diagnostic, not an error bar on the composite. The range of DFT−GFN2 offsets among the few calculated minima tests the constant-offset assumption locally; it cannot bound offsets for omitted conformers.

| Class | PBE / BLYP calculated offset range (kJ/mol) | Approximate internal population at 298 K: PBE / BLYP | At 800 K: PBE / BLYP |
| --- | --- | --- | --- |
| ps4_0000 | 0.945 / 0.712 | 0.858395 / 0.861938 | 0.977770 / 0.978097 |
| ps4_0001 | 1.054 / 1.160 | 0.330815 / 0.317068 | 0.270984 / 0.270166 |
| ps4_0010 | 1.511 / 1.586 | 0.695798 / 0.691775 | 0.859443 / 0.859328 |
| ps4_0011 | 1.482 / 1.200 | 0.932638 / 0.926908 | 0.988585 / 0.988141 |
| ps4_0101 | 1.792 / 1.924 | 0.852208 / 0.848670 | 0.982686 / 0.982569 |
| ps4_0110 | 14.360 / 11.016 | 0.770934 / 0.770148 | 0.799855 / 0.782054 |
| ps5_00000 | 0.209 / 0.145 | 0.977707 / 0.977833 | 0.998094 / 0.998094 |
| ps5_00001 | 1.402 / 1.364 | 0.701958 / 0.685163 | 0.790667 / 0.788842 |
| ps5_00010 | 1.134 / 1.215 | 0.355514 / 0.361055 | 0.853574 / 0.854699 |
| ps5_00011 | 0.432 / 0.682 | 0.896741 / 0.901810 | 0.739181 / 0.746013 |
| ps5_00100 | 0.312 / 0.483 | 0.851144 / 0.850797 | 0.967512 / 0.967482 |
| ps5_00101 | 1.699 / 1.765 | 0.695863 / 0.693138 | 0.990673 / 0.990625 |
| ps5_00110 | 2.527 / 2.919 | 0.835443 / 0.839860 | 0.978955 / 0.979221 |
| ps5_01001 | 16.556 / 13.380 | 0.094320 / 0.327201 | 0.805892 / 0.831499 |
| ps5_01010 | 0.259 / 0.253 | 0.932967 / 0.933318 | 0.995736 / 0.995749 |
| ps5_01110 | 16.356 / 13.271 | 0.819221 / 0.821571 | 0.986725 / 0.985312 |

| Step | T (K) | Pure-GFN2 local δH / δS (kJ/mol, J/mol/K) | PBE−GFN2 local δH / δS | BLYP−GFN2 local δH / δS |
| --- | ---: | --- | --- | --- |
| 2→3 | 298.15 | -10.289 / -41.462 | 0.200 / -4.631 | -5.342 / -5.488 |
| 2→3 | 300.00 | -10.273 / -41.410 | 0.207 / -4.608 | -5.338 / -5.474 |
| 2→3 | 400.00 | -9.640 / -39.560 | 0.597 / -3.483 | -4.932 / -4.321 |
| 2→3 | 500.00 | -9.231 / -38.646 | 0.856 / -2.899 | -4.493 / -3.340 |
| 2→3 | 600.00 | -8.826 / -37.908 | 1.011 / -2.615 | -4.133 / -2.682 |
| 2→3 | 700.00 | -8.229 / -36.994 | 1.134 / -2.425 | -3.835 / -2.221 |
| 2→3 | 800.00 | -7.514 / -36.038 | 1.267 / -2.248 | -3.572 / -1.870 |
| 3→4 | 298.15 | -11.387 / -42.867 | 9.644 / 0.418 | 1.935 / 1.564 |
| 3→4 | 300.00 | -11.376 / -42.830 | 9.659 / 0.470 | 1.949 / 1.611 |
| 3→4 | 400.00 | -9.881 / -38.599 | 10.312 / 2.372 | 2.489 / 3.194 |
| 3→4 | 500.00 | -8.836 / -36.215 | 10.587 / 3.001 | 2.597 / 3.453 |
| 3→4 | 600.00 | -8.869 / -36.254 | 10.606 / 3.039 | 2.404 / 3.105 |
| 3→4 | 700.00 | -9.318 / -36.945 | 10.533 / 2.929 | 2.124 / 2.674 |
| 3→4 | 800.00 | -9.740 / -37.510 | 10.429 / 2.789 | 1.846 / 2.303 |
| 4→5 | 298.15 | -11.195 / -40.040 | 0.244 / -3.237 | -4.745 / -1.845 |
| 4→5 | 300.00 | -11.200 / -40.058 | 0.240 / -3.249 | -4.743 / -1.839 |
| 4→5 | 400.00 | -11.081 / -39.791 | 0.552 / -2.420 | -4.410 / -0.908 |
| 4→5 | 500.00 | -10.368 / -38.188 | 1.269 / -0.811 | -4.050 / -0.093 |
| 4→5 | 600.00 | -9.837 / -37.229 | 1.596 / -0.204 | -3.943 / 0.107 |
| 4→5 | 700.00 | -8.182 / -34.712 | 1.672 / -0.083 | -3.946 / 0.103 |
| 4→5 | 800.00 | -5.376 / -30.970 | 1.710 / -0.034 | -3.958 / 0.088 |

The following diagnostic applies the same lowest-three-plus-offset rule to the complete n=2 and 3 electronic evidence, keeping the production one-ring reference fixed. The larger-chain energies already obey that rule. Production minus this uniform sparse diagnostic therefore measures sensitivity to the energy-treatment boundary; it does not change the reported production increments or select another method.

| Step | Level | T (K) | Production−uniform sparse local δH (kJ/mol) | Production−uniform sparse local δS (J/mol/K) |
| --- | --- | ---: | ---: | ---: |
| 2→3 | pbe | 298.15 | -6.023 | -3.450 |
| 2→3 | pbe | 300.00 | -6.025 | -3.460 |
| 2→3 | pbe | 400.00 | -6.096 | -3.673 |
| 2→3 | pbe | 500.00 | -6.050 | -3.575 |
| 2→3 | pbe | 600.00 | -5.916 | -3.333 |
| 2→3 | pbe | 700.00 | -5.728 | -3.044 |
| 2→3 | pbe | 800.00 | -5.508 | -2.750 |
| 2→3 | blyp | 298.15 | -4.261 | -3.682 |
| 2→3 | blyp | 300.00 | -4.261 | -3.683 |
| 2→3 | blyp | 400.00 | -4.139 | -3.345 |
| 2→3 | blyp | 500.00 | -3.927 | -2.873 |
| 2→3 | blyp | 600.00 | -3.693 | -2.447 |
| 2→3 | blyp | 700.00 | -3.458 | -2.085 |
| 2→3 | blyp | 800.00 | -3.228 | -1.777 |
| 3→4 | pbe | 298.15 | 6.670 | 0.754 |
| 3→4 | pbe | 300.00 | 6.685 | 0.803 |
| 3→4 | pbe | 400.00 | 7.331 | 2.684 |
| 3→4 | pbe | 500.00 | 7.631 | 3.368 |
| 3→4 | pbe | 600.00 | 7.662 | 3.431 |
| 3→4 | pbe | 700.00 | 7.555 | 3.267 |
| 3→4 | pbe | 800.00 | 7.376 | 3.028 |
| 3→4 | blyp | 298.15 | 4.368 | 1.949 |
| 3→4 | blyp | 300.00 | 4.382 | 1.996 |
| 3→4 | blyp | 400.00 | 4.926 | 3.592 |
| 3→4 | blyp | 500.00 | 5.025 | 3.830 |
| 3→4 | blyp | 600.00 | 4.820 | 3.461 |
| 3→4 | blyp | 700.00 | 4.493 | 2.959 |
| 3→4 | blyp | 800.00 | 4.141 | 2.488 |
| 4→5 | pbe | 298.15 | 0.000 | 0.000 |
| 4→5 | pbe | 300.00 | 0.000 | 0.000 |
| 4→5 | pbe | 400.00 | 0.000 | 0.000 |
| 4→5 | pbe | 500.00 | 0.000 | 0.000 |
| 4→5 | pbe | 600.00 | 0.000 | 0.000 |
| 4→5 | pbe | 700.00 | 0.000 | 0.000 |
| 4→5 | pbe | 800.00 | 0.000 | 0.000 |
| 4→5 | blyp | 298.15 | 0.000 | 0.000 |
| 4→5 | blyp | 300.00 | 0.000 | 0.000 |
| 4→5 | blyp | 400.00 | 0.000 | 0.000 |
| 4→5 | blyp | 500.00 | 0.000 | 0.000 |
| 4→5 | blyp | 600.00 | 0.000 | 0.000 |
| 4→5 | blyp | 700.00 | 0.000 | 0.000 |
| 4→5 | blyp | 800.00 | 0.000 | 0.000 |

### Literature cross-check

**Complete copies of the two early RIS papers could not be fetched.** The reached primary indexed excerpts support the specific parameters below; a full RIS entropy-per-dyad reconstruction was not reproduced. No literature parameter was fitted.

Yoon, Sundararajan and Flory give rounded relative statistical-weight prefactors 0.8 and 1.3. Their well-shape entropy terms, calculated here as R ln(prefactor), are -1.855321 and 2.181420 J/mol/K. These small local terms are not an absolute per-dyad entropy or a GAV residual. [Original study, p. 781](https://electronicsandbooks.com/edt/manual/Magazine/M/Macromolecules/1975%20%28Vol%208%29/No06%28691-959%29/776.pdf).

Relative RIS weights determine state-probability ratios. Multiplying every transfer-matrix weight by the same temperature-independent constant leaves those probabilities unchanged but adds R ln(constant) per step to the partition-derived entropy. An absolute intrawell reference and matching basin definitions are therefore needed to compare them with the present continuous coupled-rotor and whole QM-minus-GAV increment. A discrete-state Shannon entropy alone does not supply that reference. This is a mathematical reference limitation, independent of the unavailable full-paper copies.

Williams and Flory report a relative conformer preference of −700 cal/mol with an opposing entropy difference of −1.4 cal/mol/K: -2.928800 kJ/mol and -5.857600 J/mol/K after unit conversion. Phenyl rotational restriction has a negative relative entropy contribution. This relative-conformer quantity does not establish the sign or magnitude of an entire chain-increment correction against GAV. [Original study, discussion after eq. 28](https://www.electronicsandbooks.com/edt/manual/Magazine/J/Journal%20of%20the%20American%20Chemical%20Society%20US/1969%20%20%28vol%20091%29/12%20%20%283111-3408%29/3111-3118.pdf).

Khare and Paulaitis explicitly study coupled phenyl/backbone motions in polystyrene hexamers, supporting the use of coupled torsions. The accessible abstract supplies no absolute entropy-per-dyad comparator. The literature therefore supports the mechanism of coupling and relative entropy penalties, while agreement in sign and magnitude with the computed QM-minus-GAV chain increment is not established from the sources reached. [Primary university record and abstract](https://pure.johnshopkins.edu/en/publications/molecular-simulations-of-cooperative-ring-flip-motions-in-single-).

### Recorded computation cost

An eight-hour per-search timeout was an execution limit added by this worker, rather than the dispatched total-wall limit. Only each recorded GNU timeout wrapper was held while its unchanged CREST child continued. Original-deadline guardians preserve the actual wrapper status and scientific-child evidence separately. Guard outcomes: ps5_00110: scientific completed=True, wrapper exit=124; ps5_01001: scientific completed=True, wrapper exit=124. Successful acceptance requires time-v exit 0, normal CREST termination and source hashes; an unfinished child at the original cap is not accepted. A controlled sleep process reproduced a successful child with wrapper status 124. Preparation-supervisor reload evidence and prior logs remain in scratch. Neither a wrapper hold nor a lane recovery extends the original 48-hour deadline.

Through the recorded cutoff: 63.611 h wall since preparation; 300.416 CPU h in 462 recorded jobs; maximum single-job RSS 2.791 GiB. recorded leaf scientific jobs through this cutoff; enclosing supervisors excluded to avoid double counting CPU; build and small authoring/audit processes excluded; copied prior results incurred no new quantum cost.
The original scientific deadline is 2026-10-05 12:16:40 UTC. The cost cutoff is 15.611 h after that deadline. Report rendering and cached numerical audits after the cap do not allocate additional quantum time.
The owner subsequently authorized a hard 16-hour extension from ruling time, with absolute deadline 2026-10-06T06:17:16+00:00 (epoch 1791267436), to finish n=5 from local baseline `0cda64e8817941902deda3585ef2b8f7c51dd4f7` (22/25 cases). This authorization supersedes the original production cutoff for the resumed jobs; the original start and scientific method declaration are unchanged. The two interrupted searches restarted with the installed CREST 3.0.2 GFN2/quick/6 kcal/mol/four-thread settings. A requested native restart did not recover their checkpoint stages: actual streams show new metadynamics, so these are fresh same-settings restarts. Prior trees and launch receipts are preserved. No per-case eight-hour timer or n=6 allocation is used. Extension leaf scientific CPU recorded since its initialization: 62.596 h. Cached rendering/audits after the new cap allocate no quantum time.

7239 successful resource snapshots observed an aggregate calculation RSS peak of 8.673 GiB with threads confined to the eight declared physical cores. Monitoring exceptions and a gap are described below; these periodic observations are not a continuous peak-memory or affinity proof. Production DFT uses three lanes on 3/3/2 distinct physical cores, native OpenMP bounded by each lane, and one thread in each BLAS pool. The workspace setting is 5000 MB per job; observed RSS is measured separately.

Completed cases received minima checks and electronic calculations while remaining searches continued on the same eight pinned physical cores. The primary pipeline prepared one case at a time. A supplemental producer later used a second preparation lane, with minimum checks on core 10 and DFT on cores 6, 8 and 10. Per-case process locks serialize minimum and electronic writes. The supplemental producer finishes its active case when all searches complete; the bulk electronic stage waits on its process lock before launching the declared three lanes. The already-active ps4_0000 was excluded from supplemental scheduling because it began before the case locks were installed. This scheduling overlap did not change candidate selection, energy levels, sequence weights or sampling counts.
The final bulk electronic stage distributes independent conformer/functional single points across those same three lanes, including when only one case remains. It holds the case locks, retains the frozen lowest-three selection, and applies the existing offset and symmetry-image rules only after all direct points complete. The saved task-to-lane matrix and producing jobs are audited; a development check also compares cached energy-file hashes before and after this scheduling refactor.

The primary queue was restarted after its preparation CLI was confirmed to have no calculation children and to be waiting on the supplemental-owned ps4_0110 minimum lock. The interrupted waiting-stage exit 143, timing record and original logs were preserved; scientific producers continued. The restarted queue checks case ownership and its children return status 75 if another producer acquired the case. Such a case remains pending. `workflow_queue_restart.json` records this change; the original wall deadline is retained.

The first prepared tetramer also received its declared rotor integrations on core 6 during the searches. A per-case process lock protects these integrations and later cached bulk replay from overlapping writes. This added no physical core or scientific reduction.

A later early-rotor producer integrated prepared cases on core 12 during the searches. It finishes its active case when all searches complete, then releases a process lock required by the bulk eight-lane rotor stage. Both declared sampling counts are retained.

Later primary minimum preparation divided independent candidates between two single-core lanes; bulk preparation assigned the same eight cores among cases still missing their pools. The original supplemental process retained its single check lane. Recorded stage layouts: 8 case(s) with 2 lane(s). Each candidate batch and its core is recorded in `checks/execution.json`; the Verifier checks the complete unchanged candidate set and single-thread execution metadata. Cached candidates retain their original computation metadata. No selection, threshold, Hessian treatment or energy level changes.

One development partial-thermochemistry replay omitted explicit thread environment settings; its actual BLAS use was not measured. The monitor stopped during that replay on an affinity escape whose argv/mask were not captured. Its first restart identified an unpinned tee logger, which was corrected. Quantum-job environments inspected during this period had the caps set. The partial analysis was repeated with caps, and analysis entry points now set them before scientific imports. `development_thread_cap_exception.json` records the exceptions and monitoring gap; the final report uses capped reproductions.
<!-- END I048:measured -->
