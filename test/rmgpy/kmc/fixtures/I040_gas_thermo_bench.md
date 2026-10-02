# I040 — independent gas thermochemistry benchmark

Independent small-molecule data constrain RMG's additive propagation increment, but do not establish the real oligomer gas thermochemistry. The balanced surrogate's enthalpy discrepancy lies within the sum of its heterogeneous source uncertainty magnitudes. Entropy, the complete heat-capacity increment, and two-ring nonadditivity remain less securely constrained. The surrogate is an exact identity inside this pinned additive model; transferring its experimental error to an oligomer remains conditional. No correction is selected, and no result is ranked against a ceiling-temperature comparator.

## Reproduction and verifier

Run from `/home/alon/Code/RMG-Py-kmc-i040-gas-thermo-bench`, branch `i040-gas-thermo-bench`, base `0bbfe8ea33781b34ec6a29d12512aaa8f73f1e45`. Use the supplied environment. Build products are ignored. Both streams must be persisted.

```bash
mkdir -p /home/alon/runs/i040-gas-thermo-bench/build
PYTHONPATH=$PWD /home/alon/anaconda3/envs/rmg_env/bin/python setup.py build_ext --inplace > >(tee -a /home/alon/runs/i040-gas-thermo-bench/build/stdout.log) 2> >(tee -a /home/alon/runs/i040-gas-thermo-bench/build/stderr.log >&2)
mkdir -p /home/alon/runs/i040-gas-thermo-bench/verify
PYTHONPATH=$PWD /home/alon/anaconda3/envs/rmg_env/bin/python test/rmgpy/kmc/fixtures/i040_probe/run_probe.py --scratch /home/alon/runs/i040-gas-thermo-bench/verify > >(tee -a /home/alon/runs/i040-gas-thermo-bench/verify/stdout.log) 2> >(tee -a /home/alon/runs/i040-gas-thermo-bench/verify/stderr.log >&2)
```

There is one committed script. Every invocation materializes the pinned allowlist via I034's `git show` helper, reestimates all species, recompiles the baseline, verifies reaction atom balance, compares exact net source vectors and thermo over the Cp grid, checks roots by bisection, Brent and RMG Kc, and compares the complete generated report byte for byte. It writes `results.json` under scratch. `--write-report` is the authoring mode; omit it for verification. Literature values are transcribed inputs in the script, with original citations below; rerunning verifies their arithmetic and model comparisons, not new measurements or a live website scrape.

Database `4a12d36fcdc193ede82c8d1ab5c1653495d445bc`: 65 allowlisted files, SHA256 `61e2c7357bb4945dce6e6a1a39eccf216e945a01771e461391e7975165ea29eb`. No database mutation, rate-tree generation, prohibited dataset access, product-code changes, or QM runs.

## State, structures, and baseline

All RMG tables use 298.15 K, gas ideal standard pressure 100000 Pa. Concentration is 1000 mol/m³. H is formation enthalpy or reaction enthalpy in kJ/mol, S and Cp in J/mol/K. The Ince radical comparison alone uses its published 298.00 K anchor; its RMG comparison is evaluated at that same temperature.

Fresh compiled baseline: ΔH° = -80.039992; ΔS° = -147.413327; ΔCp°(298.15) = -0.474905. Continuous gas Tc = **710.020462 K**, compiler-grid Tc = 710.248673 K. The selected end radical is primary, not benzylic.

```text
C=Cc1ccccc1 + [CH2]C(CCc1ccccc1)c1ccccc1 -> [CH2]C(CC(CCc1ccccc1)c1ccccc1)c1ccccc1
```

Tc solves `H°(T) - T S°(T) - RT ln(c0 RT/p0) = 0`, equivalently `Kc(T)c0 = 1`. Thus `Hc = H°+RT` and `Sc = S°+R[ln(c0 RT/p0)+1]` for this association. Gas entropy errors are never inserted into a liquid-monomer standard.

Implicit sensitivities at the uncorrected root are dTc/d(δH) = -9.825656 K/(kJ/mol), dTc/d(δS) = 6.976416 K/(J/mol/K); actual shifts below use nonlinear roots, not this linear approximation.

## Main independent-reference table

δ means reference minus RMG. Quoted uncertainties are source magnitudes, not RMG confidence intervals. A blank reference means no defensible independent value was retrieved in this probe. ATcT is the explicitly versioned archive retrieved, not asserted to be today's latest network.

| Species | RMG H | Reference H ± u | δH | RMG S | Reference S ± u | δS | Sources |
| --- | --- | --- | --- | --- | --- | --- | --- |
| styrene | 147.421 | 148.59 ± 0.51 | +1.169 | 345.085 | 345.10 ± 2.10 | +0.015 | [A-sty](https://atct.anl.gov/Thermochemical%20Data/version%201.140/species/?species_number=312)/[N-sty](https://webbook.nist.gov/cgi/cbook.cgi?ID=C100425&Mask=1) |
| ethylbenzene | 29.098 | 29.98 ± 0.52 | +0.882 | 360.514 | 360.60 ± 0.50 | +0.086 | [A-ethyl](https://atct.anl.gov/Thermochemical%20Data/version%201.140/species/?species_number=395)/[N-ethyl](https://webbook.nist.gov/cgi/cbook.cgi?ID=C100414&Mask=1) |
| cumene | 2.658 | 3.90 ± 1.10 | +1.242 | 388.712 | 386.53 (u missing) | -2.182 | [N-cum](https://webbook.nist.gov/cgi/cbook.cgi?ID=C98828&Mask=1)/[N-cum](https://webbook.nist.gov/cgi/cbook.cgi?ID=C98828&Mask=1) |
| 1,3-diphenylpropane | 122.918 | 114.70 ± 4.10 | -8.218 | 518.461 | 498.90 (u missing) | -19.561 | [V24](https://www.mdpi.com/2673-9801/4/3/15)/[V24](https://www.mdpi.com/2673-9801/4/3/15) |
| 1,3-diphenylbutane | 96.479 | — | unidentified | 558.185 | — | unidentified | missing |
| 2,4-diphenylpentane | 70.039 | — | unidentified | 586.383 | — | unidentified | missing |
| 1-phenylethyl radical† | 179.540 | 174.90 (u missing) | -4.640 | 340.134 | 359.00 (u missing) | +18.866 | [I17](https://biblio.ugent.be/publication/8525511/file/8525521.pdf)/[I17](https://biblio.ugent.be/publication/8525511/file/8525521.pdf) |
| cumyl radical‡ | 135.279 | 149.49 (u missing) | +14.213 | 366.370 | 405.13 (u missing) | +38.759 | [J21](https://pmc.ncbi.nlm.nih.gov/articles/PMC8343544/#tbl3)/[J21](https://pmc.ncbi.nlm.nih.gov/articles/PMC8343544/#tbl3) |
| ethane | -85.346 | -84.05 ± 0.12 | +1.296 | 230.465 | 229.16 ± 0.10 | -1.305 | [A-index](https://atct.anl.gov/Thermochemical%20Data/version%201.140/)/[N-ethane](https://cccbdb.nist.gov/exp2x.asp?casno=74840&charge=0) |
| n-propylbenzene | 8.474 | 7.82 ± 0.84 | -0.654 | 399.939 | 397.86 (u missing) | -2.079 | [N-nprop](https://webbook.nist.gov/cgi/cbook.cgi?ID=C103651&Mask=1)/[N-nprop](https://webbook.nist.gov/cgi/cbook.cgi?ID=C103651&Mask=1) |
| P1 primary radical | 234.101 | — | unidentified | 374.804 | — | unidentified | missing |
| P2 primary radical | 301.482 | — | unidentified | 573.731 | — | unidentified | missing |
| P3 primary radical | 368.862 | — | unidentified | 771.402 | — | unidentified | missing |
| B2 benzylic radical | 246.940 | — | unidentified | 537.870 | — | unidentified | missing |
| B3 benzylic radical | 314.320 | — | unidentified | 735.542 | — | unidentified | missing |


† The displayed RMG radical row is at the reference anchor, whereas the Cp table below also exposes the standard RMG anchor. The radical enthalpy uses G4 with published bond-additivity correction; the raw G4 value is also retained below. No uncertainty is supplied for that individual calculation. Its molecular entropy, rather than the intrinsic group entropy, is used.

V24's diphenylpropane result is corrected **G3MP2**, with a published expanded uncertainty (95%); it does not meet the requested G4/CBS/W1 tier. Its gas entropy has no supplied uncertainty. This entry is supporting evidence, not a substitute for a higher-level independent benchmark. The RMG benzylic correction's CBS-QB3 fit already includes phenylethyl chemistry; reusing that training fit as independent evidence would be circular. I17 provides a separate G4 calculation.

The radical raw G4 H is 178.1 kJ/mol, versus BAC H 174.9 kJ/mol. Their difference is a method variant, not a confidence interval.

‡ Cumyl H/S are **derived conditionally**, not tabulated absolute values: `Hf(radical)=Hf(parent)+BDE−Hf(H)` and `S(radical)=S(parent)+(BDE−BDFE)/T−S(H)`. These use the independent NIST parent and atomic-H anchors with J21's gas CBS values, assuming their dissociation thermochemistry is at the nominal anchor with matching gas pressure conventions. The accessible main paper verifies that the BDE/BDFE difference largely reflects gas H-atom entropy, but its SI temperature/rotor details were not retrieved. Its molecule-specific method uncertainties are absent; the W1BD benchmark mean error is not a species uncertainty. These derived entropies do not establish a converged hindered-rotor/conformer entropy. [J21](https://pmc.ncbi.nlm.nih.gov/articles/PMC8343544/#tbl3), [N-H](https://webbook.nist.gov/cgi/cbook.cgi?ID=C12385136&Mask=1).

The source entropy standard pressures are not explicit on every accessible compilation page. Published S values are treated as nominal one-bar values; this is an unresolved convention uncertainty for legacy entries. If every entry instead used one atmosphere, each S would increase by 0.109443 J/mol/K on conversion to one bar. The surrogate has Δν = −1, so its reference ΔS would decrease by that amount; mixed-source pressure conventions require source-specific resolution. No such normalization is silently applied.

## Heat capacities

Each cell is RMG Cp at the column temperature. NIST recommended Cp arrays are comparison inputs, not independently measured values with assigned uncertainties.

| Species | 298.15 K | 300 K | 400 K | 500 K | 600 K | 700 K | 800 K | 1000 K | 1500 K |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| styrene | 122.106 | 122.800 | 160.331 | 192.213 | 218.154 | 237.526 | 256.898 | 284.177 | 324.553 |
| ethylbenzene | 128.520 | 129.286 | 170.665 | 206.522 | 236.187 | 258.613 | 281.039 | 312.963 | 359.698 |
| cumene | 150.274 | 151.168 | 199.493 | 241.333 | 276.060 | 302.106 | 328.151 | 365.180 | 420.366 |
| 1,3-diphenylpropane | 228.397 | 229.785 | 304.804 | 368.903 | 421.203 | 460.324 | 499.444 | 553.962 | 631.910 |
| 1,3-diphenylbutane | 250.151 | 251.668 | 333.632 | 403.714 | 461.077 | 503.816 | 546.556 | 606.178 | 692.578 |
| 2,4-diphenylpentane | 271.905 | 273.550 | 362.460 | 438.525 | 500.950 | 547.309 | 593.668 | 658.394 | 753.246 |
| 1-phenylethyl radical | 128.430 | 129.098 | 165.225 | 197.322 | 225.498 | 247.856 | 270.214 | 301.217 | 336.127 |
| cumyl radical | 151.485 | 152.298 | 196.230 | 234.890 | 267.441 | 292.064 | 316.687 | 351.833 | 407.020 |
| ethane | 51.542 | 51.798 | 65.605 | 78.659 | 90.291 | 99.621 | 108.951 | 123.595 | 147.109 |
| n-propylbenzene | 151.420 | 152.298 | 199.744 | 241.040 | 275.307 | 301.332 | 327.356 | 364.594 | 419.320 |
| P1 primary radical | 125.344 | 126.064 | 164.975 | 198.531 | 226.145 | 246.982 | 267.818 | 297.315 | 340.201 |
| P2 primary radical | 247.948 | 249.408 | 328.360 | 395.681 | 450.784 | 491.829 | 532.874 | 590.111 | 676.511 |
| P3 primary radical | 369.579 | 371.790 | 491.327 | 592.873 | 675.674 | 737.033 | 798.391 | 883.326 | 1009.390 |
| B2 benzylic radical | 250.061 | 251.480 | 328.192 | 394.514 | 450.388 | 493.060 | 535.731 | 594.432 | 669.007 |
| B3 benzylic radical | 371.692 | 373.862 | 491.159 | 591.705 | 675.278 | 738.263 | 801.248 | 887.646 | 1001.886 |


Reference Cp and δCp for **styrene** (N-sty; uCp missing):

| T / K | reference Cp | RMG Cp | δCp |
| --- | --- | --- | --- |
| 298.15 | 120.190 | 122.106 | -1.916 |
| 300 | 120.940 | 122.800 | -1.860 |
| 400 | 159.790 | 160.331 | -0.541 |
| 500 | 192.590 | 192.213 | +0.377 |
| 600 | 219.000 | 218.154 | +0.846 |
| 700 | 240.400 | 237.526 | +2.874 |
| 800 | 258.000 | 256.898 | +1.102 |
| 1000 | 285.200 | 284.177 | +1.023 |
| 1500 | 325.200 | 324.553 | +0.647 |


Reference Cp and δCp for **ethylbenzene** (N-ethyl; uCp missing):

| T / K | reference Cp | RMG Cp | δCp |
| --- | --- | --- | --- |
| 298.15 | 127.400 | 128.520 | -1.120 |
| 300 | 128.190 | 129.286 | -1.096 |
| 400 | 169.950 | 170.665 | -0.715 |
| 500 | 206.580 | 206.522 | +0.058 |
| 600 | 236.750 | 236.187 | +0.563 |
| 700 | 261.510 | 258.613 | +2.897 |
| 800 | 282.080 | 281.039 | +1.041 |
| 1000 | 314.040 | 312.963 | +1.077 |
| 1500 | 361.270 | 359.698 | +1.572 |


Reference Cp and δCp for **1-phenylethyl radical** (I17; uCp missing):

| T / K | reference Cp | RMG Cp | δCp |
| --- | --- | --- | --- |
| 300 | 124.600 | 129.098 | -4.498 |
| 400 | 164.600 | 165.225 | -0.625 |
| 500 | 199.000 | 197.322 | +1.678 |
| 600 | 227.100 | 225.498 | +1.602 |
| 800 | 269.100 | 270.214 | -1.114 |
| 1000 | 298.800 | 301.217 | -2.417 |
| 1500 | 343.100 | 336.127 | +6.973 |


No retrieved gas Cp reference for cumene, n-propylbenzene, the three diphenyl compounds or cumyl is represented as a zero error. These missing functions prevent a fully independent high-temperature propagation prediction.

## Independent radical bond cycles

J21's primary table was visually read from its archived image. The CBS columns below are calculated gas dissociation quantities; the Luo column is a compilation reproduced by the authors, without an uncertainty on these rows. It is retained as a conflicting comparison, not chosen by agreement with RMG. H/S derived from different parent anchors are not the internally consistent absolute CBS species values. All cycle values assume the nominal anchor; the entropy-derived column is conditional as above.

| Parent → radical + H | RMG BDE / kJ | CBS BDE / kJ | RMG BDFE / kJ | CBS BDFE / kJ | Derived CBS radical H | Derived CBS radical S | Derived Luo radical H |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ethylbenzene | 368.462881 | 366.518400 | 340.312596 | 330.954400 | 178.500400 | 365.165240 | 169.295600 |
| cumene | 350.622500 | 363.589600 | 323.076539 | 323.841600 | 149.491600 | 405.128445 | 134.010800 |


The cumyl enthalpy method spread is substantial and cannot be settled by assigning a fabricated uncertainty or by choosing whichever radical value agrees with a desired Tc. Benzyl_T still has coefficient zero in the compiled step. Its possible species-level error therefore does not shift that step when confined to this group; a correction to shared saturated groups remains a different inference.

## Model reactions and group cancellation

Reaction coefficients below are products positive, reactants negative. All are checked for atom balance. Only the two complete small-molecule reference reactions have independent numeric ΔH/ΔS. Missing product-radical or oligomer data prevents independent values for the direct addition reactions.

| Reaction | RMG ΔH | RMG ΔS | Reference ΔH ± source-magnitude bound | Reference ΔS |
| --- | --- | --- | --- | --- |
| closed ethyl | -80.039992 | -147.413327 | missing: 1,3-diphenylbutane | unidentified |
| closed cumene | -80.039992 | -147.413327 | missing: 2,4-diphenylpentane | unidentified |
| smallest benzylic | -80.039992 | -147.413327 | missing: 1-phenylethyl radical, B2 benzylic radical | unidentified |
| smallest primary | -80.039846 | -146.157638 | missing: P1 primary radical, P2 primary radical | unidentified |
| compiled P2 to P3 | -80.039992 | -147.413327 | missing: P2 primary radical, P3 primary radical | unidentified |
| benzylic B2 to B3 | -80.039992 | -147.413327 | missing: B2 benzylic radical, B3 benzylic radical | unidentified |
| balanced group surrogate | -80.039992 | -147.413327 | -82.80 ± 3.09 | -150.470000 (u incomplete) |
| two-ring nonadditivity | -0.000000 | 11.526306 | 7.15 ± 5.58 | 30.400000 (u incomplete) |


**closed ethyl:** `-1 ethylbenzene; -1 styrene; +1 1,3-diphenylbutane`

Net source vector: `{"group:Cb-(Cds-Cds)": -1, "group:Cb-Cs": 1, "group:Cds-CdsCbH": -1, "group:Cds-CdsHH": -1, "group:Cs-CbCsCsH": 1, "group:Cs-CsCsHH": 1}`.

**closed cumene:** `-1 cumene; -1 styrene; +1 2,4-diphenylpentane`

Net source vector: `{"group:Cb-(Cds-Cds)": -1, "group:Cb-Cs": 1, "group:Cds-CdsCbH": -1, "group:Cds-CdsHH": -1, "group:Cs-CbCsCsH": 1, "group:Cs-CsCsHH": 1}`.

**smallest benzylic:** `-1 1-phenylethyl radical; -1 styrene; +1 B2 benzylic radical`

Net source vector: `{"group:Cb-(Cds-Cds)": -1, "group:Cb-Cs": 1, "group:Cds-CdsCbH": -1, "group:Cds-CdsHH": -1, "group:Cs-CbCsCsH": 1, "group:Cs-CsCsHH": 1}`.

**smallest primary:** `-1 P1 primary radical; -1 styrene; +1 P2 primary radical`

Net source vector: `{"group:Cb-(Cds-Cds)": -1, "group:Cb-Cs": 1, "group:Cds-CdsCbH": -1, "group:Cds-CdsHH": -1, "group:Cs-CbCsCsH": 1, "group:Cs-CsCsHH": 1, "radical:Isobutyl": 1, "radical:RCCJ": -1}`.

**compiled P2 to P3:** `-1 P2 primary radical; -1 styrene; +1 P3 primary radical`

Net source vector: `{"group:Cb-(Cds-Cds)": -1, "group:Cb-Cs": 1, "group:Cds-CdsCbH": -1, "group:Cds-CdsHH": -1, "group:Cs-CbCsCsH": 1, "group:Cs-CsCsHH": 1}`.

**benzylic B2 to B3:** `-1 B2 benzylic radical; -1 styrene; +1 B3 benzylic radical`

Net source vector: `{"group:Cb-(Cds-Cds)": -1, "group:Cb-Cs": 1, "group:Cds-CdsCbH": -1, "group:Cds-CdsHH": -1, "group:Cs-CbCsCsH": 1, "group:Cs-CsCsHH": 1}`.

**balanced group surrogate:** `-1 ethylbenzene; -1 ethane; -1 styrene; +1 cumene; +1 n-propylbenzene`

Net source vector: `{"group:Cb-(Cds-Cds)": -1, "group:Cb-Cs": 1, "group:Cds-CdsCbH": -1, "group:Cds-CdsHH": -1, "group:Cs-CbCsCsH": 1, "group:Cs-CsCsHH": 1}`.

**two-ring nonadditivity:** `-1 1,3-diphenylpropane; -1 ethane; +1 ethylbenzene; +1 n-propylbenzene`

Net source vector: `{}`.

The balanced surrogate is **ethylbenzene + ethane + styrene → cumene + n-propylbenzene**. Its non-symmetry source vector, ΔH(T), and ΔCp(T) exactly equal the compiled increment on every reported temperature. Compiled ΔS minus surrogate ΔS is a temperature-independent +0.000000 J/mol/K, explicitly retained. Group identity does not prove equality of real-molecule nonadditivity, conformer populations, or stereochemical ensembles.

The two-ring nonadditivity reaction is **diphenylpropane + ethane → ethylbenzene + n-propylbenzene**. Its entire additive source vector vanishes, so its nonzero supporting reference enthalpy suggests physics invisible to those groups, subject to the lower-level diphenylpropane calculation. It cannot be assigned uniquely to a propagation group or used as a unique correction. Diphenylbutane/pentane stereochemistry is unspecified in these RMG graphs; RMG's local optical factors are not an independently established meso/racemic equilibrium ensemble.

The smallest primary addition is explicitly checked separately: its radical environment changes and it need not reproduce the saturated long-chain increment. Benzylic self-similar propagation cancels the same Benzyl_S correction at both ends. The compiled primary propagation cancels Isobutyl at both ends. Consequently an error confined to either radical correction gives **δΔH = 0, δΔS = 0, δTc = 0** for its self-similar step. The cumyl Benzyl_T correction has coefficient zero in the compiled step. Whole-molecule radical discrepancies do not identify errors in the individual HBI group; library, saturated skeleton, symmetry and Cp contributions must be separated first.

## Conditional effects on compiled propagation and Tc

Rows change one published surrogate input at a time, keeping the other inputs and RMG ΔCp fixed. The styrene rows are direct one-species replacements in the actual compiled reaction. The other rows transfer through the proved additive identity and therefore require absence of additional oligomer nonadditivity. Their signs follow stoichiometry. They are controls, not independent additive uncertainties; do not add the combined row again. Tc intervals are nonlinear endpoint envelopes from the published input uncertainty magnitude alone. They exclude model error and omitted Cp/nonadditivity/stereo/pressure uncertainties. An unavailable bound is printed explicitly.

| Input replaced | Axis | δΔH | δΔS | Tc / K | δTc / K | Input-bound Tc / K |
| --- | --- | --- | --- | --- | --- | --- |
| ethylbenzene | H | -0.881927 | +0.000000 | 718.694324 | +8.673862 | [713.578, 723.816] |
| ethylbenzene | S | +0.000000 | -0.086235 | 709.419404 | -0.601058 | [705.956, 712.920] |
| ethane | H | -1.295870 | +0.000000 | 722.771270 | +12.750808 | [721.589, 723.954] |
| ethane | S | +0.000000 | +1.304897 | 719.251791 | +9.231328 | [718.535, 719.970] |
| styrene | H | -1.169368 | +0.000000 | 721.524958 | +11.504496 | [716.504, 726.552] |
| styrene | S | +0.000000 | -0.015160 | 709.914719 | -0.105743 | [695.587, 724.894] |
| cumene | H | +1.241544 | +0.000000 | 697.838050 | -12.182412 | [687.072, 708.630] |
| cumene | S | +0.000000 | -2.181609 | 695.144266 | -14.876196 | u missing |
| n-propylbenzene | H | -0.654388 | +0.000000 | 716.454850 | +6.434388 | [708.197, 724.728] |
| n-propylbenzene | S | +0.000000 | -2.078567 | 695.831788 | -14.188674 | u missing |


Combining the complete surrogate references gives δΔH = **-2.760008 kJ/mol**, δΔS = **-3.056673 J/mol/K**, retaining the model symmetry offset. Enthalpy alone gives Tc 737.221104 K; entropy alone 689.363999 K; both 715.649201 K. This is a small-molecule transfer exercise with unchanged RMG ΔCp, not an independent oligomer Tc estimate.

The sum of published H uncertainty magnitudes is ±3.090 kJ/mol, giving the enthalpy-only Tc envelope [706.779, 767.867] K. ATcT H uncertainties use its network confidence convention; legacy NIST magnitudes are not all harmonized confidence intervals. This triangle-inequality envelope is conditional on interpreting each supplied magnitude as a bound, not a guaranteed statistical interval. The ATcT styrene–ethylbenzene correlation is +0.433; with both coefficients negative its covariance contribution is positive. An illustrative RSS using that pair and treating all other sources as uncorrelated is 1.640172 kJ/mol, **not** a defensible combined confidence interval. The S envelope is unavailable because cumene/n-propylbenzene lack reported uS.

For the styrene and ethylbenzene Cp discrepancies, the code integrates piecewise-linear published Cp minus RMG Cp from the anchor, separately for δH(T) and δS(T), then applies their negative surrogate coefficients. Styrene transfers directly; ethylbenzene remains conditional on the group identity. This is additional to the constant anchor errors. Separate H-only or S-only Cp rows are accounting controls; the thermodynamically consistent correction includes both. No Cp uncertainty envelope can be supplied.

| Input: Cp contribution | δΔH at baseline Tc / kJ/mol | δΔS at baseline Tc | Tc / K | δTc / K |
| --- | --- | --- | --- | --- |
| styrene:H | -0.143355 | +0.000000 | 711.467412 | +1.446950 |
| styrene:S | +0.000000 | -0.044304 | 709.719502 | -0.300960 |
| styrene:HS | -0.143355 | -0.044304 | 711.119572 | +1.099110 |
| ethylbenzene:H | -0.106675 | +0.000000 | 711.097326 | +1.076864 |
| ethylbenzene:S | +0.000000 | -0.010601 | 709.948429 | -0.072033 |
| ethylbenzene:HS | -0.106675 | -0.010601 | 710.994645 | +0.974183 |


H anchor + S anchor + Cp correction for styrene alone gives Tc 722.511182 K; for ethylbenzene alone, via the surrogate, 719.057533 K. Combining all surrogate anchor errors with only the known styrene/ethylbenzene Cp discrepancies gives 717.659206 K. All of these controls retain the other species' RMG Cp functions; a complete independent ΔCp remains unavailable.

Phenylethyl's benchmark discrepancy can test Benzyl_S only after subtracting the independently assessed saturated parent and consistent entropy conventions. Even such an identified Benzyl_S correction cancels in self-similar benzylic propagation, and Benzyl_S is absent from the compiled primary propagation. Diphenylpropane's whole-molecule residual has no unique transfer coefficient; diphenylbutane, diphenylpentane and the actual oligomer radicals have no complete independent references here. Cumyl has a conditional cycle reference with unresolved method/entropy uncertainties and no unique shared-group attribution. For those whole-molecule residuals δΔH, δΔS and δTc are **unidentified**, rather than set to zero.

## Exact species provenance

These are the actual comments, source weights, aliases, graph symmetries and optical-site counts used by RMG. Each H/S decomposition is checked against the estimated species at the I039 decomposition temperatures. Full source entry descriptions and quantity uncertainty strings are also saved in scratch results.json.

**styrene** — `C=Cc1ccccc1`; σ = 2; optical half-factor sites = 0.

```text
Thermo group additivity estimation: group(Cb-(Cds-Cds)) + group(Cb-H) + group(Cb-H) + group(Cds-CdsCbH) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cds-CdsHH) + ring(Benzene)
{"group:Cb-(Cds-Cds)": 1, "group:Cb-H": 5, "group:Cds-CdsCbH": 1, "group:Cds-CdsHH": 1, "ring:Benzene": 1}
```

**ethylbenzene** — `CCc1ccccc1`; σ = 6; optical half-factor sites = 0.

```text
Thermo group additivity estimation: group(Cs-CbCsHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene)
{"group:Cb-Cs": 1, "group:Cb-H": 5, "group:Cs-CbCsHH": 1, "group:Cs-CsHHH": 1, "ring:Benzene": 1}
```

**cumene** — `CC(C)c1ccccc1`; σ = 18; optical half-factor sites = 0.

```text
Thermo group additivity estimation: group(Cs-CbCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene)
{"group:Cb-Cs": 1, "group:Cb-H": 5, "group:Cs-CbCsCsH": 1, "group:Cs-CsHHH": 2, "ring:Benzene": 1}
```

**1,3-diphenylpropane** — `c1ccc(CCCc2ccccc2)cc1`; σ = 8; optical half-factor sites = 0.

```text
Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CbCsHH) + group(Cs-CbCsHH) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene) + ring(Benzene)
{"group:Cb-Cs": 2, "group:Cb-H": 10, "group:Cs-CbCsHH": 2, "group:Cs-CsCsHH": 1, "ring:Benzene": 2}
```

**1,3-diphenylbutane** — `CC(CCc1ccccc1)c1ccccc1`; σ = 6; optical half-factor sites = 1.

```text
Thermo group additivity estimation: group(Cs-CbCsCsH) + group(Cs-CsCsHH) + group(Cs-CbCsHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene) + ring(Benzene)
{"group:Cb-Cs": 2, "group:Cb-H": 10, "group:Cs-CbCsCsH": 1, "group:Cs-CbCsHH": 1, "group:Cs-CsCsHH": 1, "group:Cs-CsHHH": 1, "ring:Benzene": 2}
```

**2,4-diphenylpentane** — `CC(CC(C)c1ccccc1)c1ccccc1`; σ = 18; optical half-factor sites = 2.

```text
Thermo group additivity estimation: group(Cs-CbCsCsH) + group(Cs-CbCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene) + ring(Benzene)
{"group:Cb-Cs": 2, "group:Cb-H": 10, "group:Cs-CbCsCsH": 2, "group:Cs-CsCsHH": 1, "group:Cs-CsHHH": 2, "ring:Benzene": 2}
```

**1-phenylethyl radical** — `C[CH]c1ccccc1`; σ = 6; optical half-factor sites = 0.

```text
Thermo group additivity estimation: group(Cs-CbCsHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene) + radical(Benzyl_S)
{"HBI:subtract_H_atoms": 1, "group:Cb-Cs": 1, "group:Cb-H": 5, "group:Cs-CbCsHH": 1, "group:Cs-CsHHH": 1, "radical:Benzyl_S": 1, "ring:Benzene": 1}
```

**cumyl radical** — `C[C](C)c1ccccc1`; σ = 18; optical half-factor sites = 0.

```text
Thermo group additivity estimation: group(Cs-CbCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene) + radical(Benzyl_T)
{"HBI:subtract_H_atoms": 1, "group:Cb-Cs": 1, "group:Cb-H": 5, "group:Cs-CbCsCsH": 1, "group:Cs-CsHHH": 2, "radical:Benzyl_T": 1, "ring:Benzene": 1}
```

**ethane** — `CC`; σ = 18; optical half-factor sites = 0.

```text
Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH)
{"group:Cs-CsHHH": 2}
```

**n-propylbenzene** — `CCCc1ccccc1`; σ = 6; optical half-factor sites = 0.

```text
Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CbCsHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene)
{"group:Cb-Cs": 1, "group:Cb-H": 5, "group:Cs-CbCsHH": 1, "group:Cs-CsCsHH": 1, "group:Cs-CsHHH": 1, "ring:Benzene": 1}
```

**P1 primary radical** — `[CH2]Cc1ccccc1`; σ = 4; optical half-factor sites = 0.

```text
Thermo group additivity estimation: group(Cs-CbCsHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene) + radical(RCCJ)
{"HBI:subtract_H_atoms": 1, "group:Cb-Cs": 1, "group:Cb-H": 5, "group:Cs-CbCsHH": 1, "group:Cs-CsHHH": 1, "radical:RCCJ": 1, "ring:Benzene": 1}
```

**P2 primary radical** — `[CH2]C(CCc1ccccc1)c1ccccc1`; σ = 4; optical half-factor sites = 1.

```text
Thermo group additivity estimation: group(Cs-CbCsCsH) + group(Cs-CsCsHH) + group(Cs-CbCsHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene) + ring(Benzene) + radical(Isobutyl)
{"HBI:subtract_H_atoms": 1, "group:Cb-Cs": 2, "group:Cb-H": 10, "group:Cs-CbCsCsH": 1, "group:Cs-CbCsHH": 1, "group:Cs-CsCsHH": 1, "group:Cs-CsHHH": 1, "radical:Isobutyl": 1, "ring:Benzene": 2}
```

**P3 primary radical** — `[CH2]C(CC(CCc1ccccc1)c1ccccc1)c1ccccc1`; σ = 4; optical half-factor sites = 2.

```text
Thermo group additivity estimation: group(Cs-CbCsCsH) + group(Cs-CbCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CbCsHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene) + ring(Benzene) + ring(Benzene) + radical(Isobutyl)
{"HBI:subtract_H_atoms": 1, "group:Cb-Cs": 3, "group:Cb-H": 15, "group:Cs-CbCsCsH": 2, "group:Cs-CbCsHH": 1, "group:Cs-CsCsHH": 2, "group:Cs-CsHHH": 1, "radical:Isobutyl": 1, "ring:Benzene": 3}
```

**B2 benzylic radical** — `CC(C[CH]c1ccccc1)c1ccccc1`; σ = 6; optical half-factor sites = 1.

```text
Thermo group additivity estimation: group(Cs-CbCsCsH) + group(Cs-CsCsHH) + group(Cs-CbCsHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene) + ring(Benzene) + radical(Benzyl_S)
{"HBI:subtract_H_atoms": 1, "group:Cb-Cs": 2, "group:Cb-H": 10, "group:Cs-CbCsCsH": 1, "group:Cs-CbCsHH": 1, "group:Cs-CsCsHH": 1, "group:Cs-CsHHH": 1, "radical:Benzyl_S": 1, "ring:Benzene": 2}
```

**B3 benzylic radical** — `CC(CC(C[CH]c1ccccc1)c1ccccc1)c1ccccc1`; σ = 6; optical half-factor sites = 2.

```text
Thermo group additivity estimation: group(Cs-CbCsCsH) + group(Cs-CbCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CbCsHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene) + ring(Benzene) + ring(Benzene) + radical(Benzyl_S)
{"HBI:subtract_H_atoms": 1, "group:Cb-Cs": 3, "group:Cb-H": 15, "group:Cs-CbCsCsH": 2, "group:Cs-CbCsHH": 1, "group:Cs-CsCsHH": 2, "group:Cs-CsHHH": 1, "radical:Benzyl_S": 1, "ring:Benzene": 3}
```

Resolved source entries (quoted quantities are database inputs, not newly estimated error bars):

| Matched entry | Resolved data entry | H298 | S298 | Attribution |
| --- | --- | --- | --- | --- |
| group:Cb-(Cds-Cds) | Cb-(Cds-Cds) | 5.69 +\|- 0.2 kcal/mol | -7.8 +\|- 0.1 cal/(mol*K) | Cb-Cd STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714 |
| group:Cb-Cs | Cb-Cs | 5.51 +\|- 0.13 kcal/mol | -7.69 +\|- 0.1 cal/(mol*K) | Cb-Cs BENSON |
| group:Cb-H | Cb-H | 3.3 +\|- 0.11 kcal/mol | 11.53 +\|- 0.12 cal/(mol*K) | Cb-H BENSON |
| group:Cds-CdsCbH | Cds-CdsCbH | 6.78 +\|- 0.2 kcal/mol | 6.38 +\|- 0.1 cal/(mol*K) | Cd-CbH BENSON |
| group:Cds-CdsHH | Cds-CdsHH | 6.26 +\|- 0.19 kcal/mol | 27.61 +\|- 0.1 cal/(mol*K) | Cd-HH BENSON |
| group:Cs-CbCsCsH | Cs-CbCsCsH | -0.98 +\|- 0.27 kcal/mol | -12.15 +\|- 0.15 cal/(mol*K) | Cs-CbCsCsH BENSON |
| group:Cs-CbCsHH | Cs-CbCsHH | -4.86 +\|- 0.2 kcal/mol | 9.34 +\|- 0.19 cal/(mol*K) | Cs-CbCsHH BENSON |
| group:Cs-CsCsHH | Cs-CsCsHH | -4.93 +\|- 0.05 kcal/mol | 9.42 +\|- 0.13 cal/(mol*K) | Cs-CsCsHH BENSON |
| group:Cs-CsHHH | Cs-CsHHH | -10.2 +\|- 0.12 kcal/mol | 30.41 +\|- 0.08 cal/(mol*K) | Cs-CsHHH BENSON |
| radical:Benzyl_S | Benzyl_S | 88.064 +\|- 2.4 kcal/mol | -4.8554 cal/(mol*K) | Fitted From Calculations from Hexylbenzene Library, Lawrence Lai |
| radical:Benzyl_T | Benzyl_T | 83.8 kcal/mol | -5.34 cal/(mol*K) | LAY et al. |
| radical:Isobutyl | Isobutyl | 101.1 kcal/mol | 2.91 cal/(mol*K) | LAY et al. |
| radical:RCCJ | RCCJ | 101.1 +\|- 0.2 kcal/mol | 2.61 cal/(mol*K) | LAY et al. CHEN & BOZZELLI # |
| ring:Benzene | Benzene | 0 kcal/mol | 0 cal/(mol*K) | Aromatic |


Shared Benson/Stein group input uncertainties and radical fit method estimates have unknown covariance and validation scope. They cannot be treated as independent species uncertainties or combined with the external errors as if each molecule were an independent measurement. Common groups cancel before any error propagation.

## Literature gaps and costed QM plan — description only

First calculate the actual gas compounds and reaction increments, without using any ceiling-temperature target for selection or calibration. Small-species set: styrene, ethylbenzene, cumene, n-propylbenzene, ethane, 1-phenylethyl radical and cumyl radical. Large-species set: the three diphenyl compounds, P1, P2, P3 primary radicals and B2 benzylic radical. These cover the requested species, both smallest additions and the actual compiled increment. B3 can be added if benzylic length convergence is needed; it is not budgeted initially.

Use CREST/GFN2-xTB for conformer discovery (both doublet and singlet as appropriate), followed by dispersion-aware ωB97X-D/def2-TZVP optimization/frequencies and connectivity checks. Enumerate stereoisomers before conformers, explicitly including meso and racemic diphenylpentane; report each stereoisomer and a defined ensemble, rather than silently averaging across fixed tacticities. Include electronic spin degeneracy, standard pressure, rotational symmetry, and enantiomer degeneracy consistently. Deduplicate minima, retain low-energy families and expand the energy window until ensemble H/S/Cp converge. Initial conformer quotas in the budget are starting estimates, not proof of convergence.

Use canonical G4 for the small set, with raw atomization and optional BAC results reported separately. For the large set use TightPNO DLPNO-CCSD(T1) single points at def2-TZVPP/QZVPP with documented CBS extrapolation, plus canonical G4 cross-checks on representative large conformers. Check unrestricted spin contamination and open-shell diagnostics; if multi-reference character or PNO/basis sensitivity exceeds the planned uncertainty budget, escalate the affected pair rather than forcing an increment. Reaction energies should be computed directly from consistently treated species and isodesmic cycles; external ATcT/NIST anchors test method bias.

Treat methyl, phenyl and backbone torsions with relaxed periodic scans at the DFT level and replace corresponding harmonic modes by hindered-rotor partition functions. Test coupled backbone/phenyl torsions with two-dimensional grids on representative minima, and avoid double counting conformers and rotor states. Compute Boltzmann-weighted H/S/Cp over the relevant gas-temperature range. Convergence targets: sub-kJ/mol reaction-enthalpy changes and about one J/mol/K reaction-entropy changes under conformer/rotor/basis expansion; these are proposed numerical goals, not achieved uncertainties. Use the observed convergence and cross-method/cycle residuals for a covariance-aware uncertainty budget.

Budget below is **estimated serial core-hours**, not measured wall time, and includes no QM executed in this probe. A single core-hour is one CPU core for one hour; parallel speedup/memory limits need a pilot. The costed subset has finite scan/conformer quotas, so the budget is conditional on convergence.

| Planned task | Jobs | Estimated core-h/job | Estimated core-h |
| --- | --- | --- | --- |
| small-species G4, optimization/frequency included | 21 | 600 | 12600 |
| large-species DFT optimization/frequency | 84 | 80 | 6720 |
| large-species DLPNO-CCSD(T1)/CBS single-point pair | 84 | 500 | 42000 |
| small-species 24-point relaxed rotor scans | 504 | 0.5 | 252 |
| large-species 24-point relaxed rotor scans | 1008 | 2 | 2016 |
| coupled-rotor 12-by-12 grids | 576 | 2 | 1152 |
| small-species conformer search | 7 | 20 | 140 |
| large-species conformer search | 7 | 100 | 700 |
| canonical G4 cross-checks on large representative conformers | 6 | 1500 | 9000 |


Estimated initial total: **74,580 core-hours**. A factor-of-three planning range is 24,860–223,740 core-hours. At ideal use of 32 cores this corresponds to 2330.6 aggregate node-hours; it is not an elapsed-time prediction. Run a small fixed pilot to replace unit costs before authorizing this QM campaign. No QM run is authorized by this report.

## Limits of the dispatch premise

The stated gas ceiling is reproduced, but a request to propagate *each molecular benchmark error* uniquely into propagation overstates identifiability. Shared groups, cancelling radical corrections, unknown two-ring nonadditivity and unspecified stereochemistry prevent that inference. The tables show direct replacements, explicitly conditional surrogate transfers, exact zero coefficients, and unresolved effects separately. No complete independent reference exists in the retrieved evidence for the actual P2/P3 reaction. Its total uncertainty and an independently validated gas Tc remain unsettled.

A quoted liquid-monomer ceiling is not an independent validation of an ideal-gas molecular benchmark. The source-pressure ambiguity and legacy uncertainty conventions also prevent calling all supplied S/H values a common confidence interval. None of these limitations requires changing the task's implementation scope; the report and verifier are complete, while the scientific gaps require the described additional calculations or better references.

## Primary references

- **A-sty**: [ATcT v1.140 (2024), phenylethene, species 312; H298, not H0](https://atct.anl.gov/Thermochemical%20Data/version%201.140/species/?species_number=312).

- **A-ethyl**: [ATcT v1.140 (2024), ethylbenzene, species 395; H298](https://atct.anl.gov/Thermochemical%20Data/version%201.140/species/?species_number=395).

- **A-index**: [ATcT v1.140 (2024), ethane, H298](https://atct.anl.gov/Thermochemical%20Data/version%201.140/).

- **N-sty**: [NIST WebBook, styrene: Pitzer et al. 1946 S; TRC 1997 recommended gas Cp](https://webbook.nist.gov/cgi/cbook.cgi?ID=C100425&Mask=1).

- **N-ethyl**: [NIST WebBook, ethylbenzene: Miller 1978 S; TRC 1997 recommended gas Cp](https://webbook.nist.gov/cgi/cbook.cgi?ID=C100414&Mask=1).

- **N-cum**: [NIST WebBook, cumene: Prosen et al. 1945 H; Kishimoto et al. 1973 S](https://webbook.nist.gov/cgi/cbook.cgi?ID=C98828&Mask=1).

- **N-nprop**: [NIST WebBook, n-propylbenzene: Prosen et al. 1945/1946 H; Messerly et al. 1965 S](https://webbook.nist.gov/cgi/cbook.cgi?ID=C103651&Mask=1).

- **N-ethane**: [NIST CCCBDB, ethane, Gurvich compilation S298](https://cccbdb.nist.gov/exp2x.asp?casno=74840&charge=0).

- **V24**: [Verevkin et al., Oxygen 2024, 4, 266–285, DOI 10.3390/oxygen4030015; Tables 5,7: corrected G3MP2](https://www.mdpi.com/2673-9801/4/3/15).

- **I17**: [Ince et al., AIChE J. 2017, DOI 10.1002/aic.15588; SI Table S1, T6/53: G4/BAC H, molecular S, Cp](https://biblio.ugent.be/publication/8525511/file/8525521.pdf).

- **J21**: [Salamone et al., JACS 2021, 143, 11759–11776, DOI 10.1021/jacs.1c05566; Table 3, rows 11,12, gas (RO)CBS-QB3 BDE/BDFE](https://pmc.ncbi.nlm.nih.gov/articles/PMC8343544/#tbl3).

- **N-H**: [NIST WebBook, atomic hydrogen, CODATA 1984 H298 and one-bar S298](https://webbook.nist.gov/cgi/cbook.cgi?ID=C12385136&Mask=1).


V24 accessible primary-paper copy: [PDF](https://pdfs.semanticscholar.org/d7d0/6da8ff5248dc2515cbc4ccd8b8cc19c9d56b.pdf). Its gas H is corrected G3MP2 and the gas S uncertainty is not the uncertainty of vaporization entropy. Ince's group-fit residual statistics are not uncertainty bounds on a specific G4 species.

Legacy diphenylpropane combustion compilations contain conflicting derived liquid values and lack a complete vetted gas conversion/uncertainty chain in the retrieved record; they were not substituted for a gas reference. Missing references are retrieval/validation gaps in this probe, not assertions that no reference can exist. Pyrolysis results and prohibited datasets were not used.
