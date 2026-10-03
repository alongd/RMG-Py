# I043 — local adjacent-phenyl nonadditivity probe

Status: composite ensemble calculations complete; full replay verification
passed. The measured tables distinguish sampled candidates from validated
minima. Propagation corrections and ceilings are conditional transfers from
closed-shell capped molecules, with no correction adopted in product code.

The baseline, molecular sampling, and computational pilots use local resources
only. Database materialization uses the committed I034 allowlist and `git show`
at `4a12d36fcdc193ede82c8d1ab5c1653495d445bc`; the database remains read-only.
No excluded dataset, pyrolysis result, product code, database tree, or product
test is used or changed. Worktree: `RMG-Py-kmc-i043-diphenyl-nonadditivity`,
branch `i043-diphenyl-nonadditivity`, base `0bbfe8ea3`.

## Theory selected before the corrected ceiling

The initial selection was PBE-D3(BJ)/def2-SVP and BLYP-D3(BJ)/def2-SVP, with
full DFT optimization/Hessians. The ethane pilot was completed; the large pilot
was launched to measure local cost. After observing the large-job cost, the
protocol was revised to the corresponding composite levels:

- PBE-D3(BJ)/def2-SVP//GFN2-xTB.
- BLYP-D3(BJ)/def2-SVP//GFN2-xTB.

The composite protocol uses tightly optimized GFN2-xTB minima and their xTB
Hessians, with DFT electronic energies for re-ranking and separate DFT pilots.
It is not full DFT minimum thermochemistry. Dispersion uses the installed
PySCF D3(BJ) implementation, without the optional three-body term. DFT uses
density fitting, integration grid level 3, and SCF tolerance 1e-9 Eh.
No choice depends on agreement with any literature ceiling. The initial and
revised plans are persisted as `method_plan_initial.json` and `method_plan.json`
in scratch. The reproduced *uncorrected* RMG baseline predates the revision;
no corrected Tc was calculated before either selection.
The geometry/frequency revision therefore followed a reproduced baseline Tc.
The electronic functionals and basis were already fixed; a literal requirement
to freeze the entire final composite protocol before *any* Tc was not met.

The optional full-DFT trimer pilot reached an optimized geometry, but its
Hessian was stopped when observed trimer single-point costs threatened the
24-hour budget for the required ensemble calculations. The optimization record
and both streams remain in scratch, with `frequencies_complete: false`.
It is excluded from the table of frequency-validated DFT pilots. No trimer DFT
frequency validation is claimed. All retained composite minima have their
tight xTB Hessians and frequency checks.

CREST uses GFN2-xTB, `--quick`, a 6 kcal/mol candidate window, and eight threads.
The installed help explicitly describes `--quick` as a crude ensemble search:
these searches alone do not establish sampling convergence. The seed 43002
fixes RDKit starting embeddings only; CREST dynamics retain their stochastic
default. Candidates within 12 kJ/mol are the planned refinement set, selected
without a thermochemical comparator. Symmetry duplicates and unstable candidates
must be removed after tighter refinement, not treated as additional entropy.
Missing reflection partners are added only for achiral molecular
stereoisomers, by reflecting coordinates and permuting atoms back onto the
specified stereo graph. The Hessian receives the same orthogonal transformation;
the mapping also preserves signed volumes at unlabeled CH2 and CH3 atoms,
keeping the atom-label orientation of the torsional coordinate manifold.
electronic energies are inherited by exact reflection symmetry and checked with
a fresh DFT example. Chiral sequence classes do not receive the opposite
enantiomer. These are conformational wells, not new chain configurations.
The full periodic torsion domain also needs the proper atom-permutation
images removed by graph-symmetry RMS deduplication. Angular coordinates are
closed under these operations within 0.05 rad after methyl/phenyl folding;
their Hessians and electronic energies transform exactly by atom permutation.
They are integration wells, not additional chemical conformers. An explicit
orbit audit found omitted wells for diphenylpropane and the racemo dimer;
their quadratures were rebuilt before corrected Tc. The report distinguishes
chemical minima from rotor-space images. The verifier independently checks
proper angular closure within 0.05 rad and oriented reflection coverage
within 0.10 rad, consistent with the finite geometry deduplication tolerance.
These checks guard against missing symmetry-related probability mass.

An explicit CPU affinity confines calculation processes to eight physical cores;
thread variables alone did not prevent nested runtime oversubscription.
The abandoned initial attempts predated affinity enforcement, so their physical
core use cannot be certified as at most eight. Subsequent checks inspect every
worker thread's affinity, including the BLAS/OpenMP pools; current retained
computation uses the same eight physical cores across simultaneous jobs.
An audit also found that four auxiliary `tee` loggers had unrestricted affinity,
although the computational threads were confined. Those loggers were confined
to the same eight cores, and the reproduction wrapper below constrains logging
as well as computation. Earlier full-process confinement is not claimed.
A repeated trimer single point with one BLAS thread reproduced its energy,
but took longer under the current shared load. This was not a controlled
speed comparison and supplied no evidence for changing production settings.

The initial diphenylpropane and cumene searches emitted repeated OpenBLAS
warnings that the OpenMP loop could hang. Those incomplete attempts were stopped
and preserved under `failed_openblas/`. Searches were restarted with
`OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and eight CREST threads. The restart
completed the diphenylpropane search without those warnings. The failed attempts
are not molecular thermochemistry evidence.

## Rotor/conformer partition model

Every acyclic heavy-atom single bond is a torsional coordinate: backbone,
phenyl attachment, and methyl rotation. The Cartesian xTB Hessian is projected
orthogonally to translations, external rotations, and the torsional tangent
subspace. Only the remaining normal modes enter the harmonic partition.
The torsional Schur complement checks minimum stability. A Pitzer–Gwinn
quantum correction uses the generalized eigenfrequencies of the rigid
torsional curvature and projected kinetic metric, matching the actual
fixed-angle classical potential in its harmonic limit. Its oscillator
partition is not added as another independent vibration. Using relaxed
Schur frequencies with the rigid potential would mismatch that limit; this
was corrected before any corrected Tc was evaluated.

Torsional potentials are direct GFN2-xTB energies on a fixed-bond-length,
fixed-bond-angle manifold constructed from each tightly optimized conformer.
Methyl and phenyl coordinates use their fundamental periods. A periodic
Voronoi partition in the full torsional coordinate space assigns every point
to exactly one retained conformer. Thus the rotor integral and the conformer
sum partition the same coordinate space rather than counting a whole rotor
again for every well. Couplings between torsions are included in the sampled
electronic potential; they are not approximated as a product of independent
one-dimensional sine potentials.

For basin i with d torsions and the projected kinetic metric I_i,

```text
q_tor,classical,i(T) = (2π k_B T)^(d/2) sqrt(det I_i) / h^d
                     × ∫_(basin i) exp[−(V_GFN(q)−V_GFN,min,i)/(RT)] dq
q_tor,i = q_tor,classical,i × ∏_a [x_a / (2 sinh(x_a/2))]
x_a = h c ν_rigid-tor,a / (k_B T)
```

The initial integral uses importance sampling from a wrapped local Gaussian
plus a uniform tail. A failed numerical-precision preflight for the added-unit
cycle motivated a correlated proposal before corrected Tc was calculated.
The replacement draws a 50% sequential wrapped conditional Gaussian, a 40%
product mixture with broader per-coordinate tails, and a 10% uniform component.
The Gaussian covariance derives from the same rigid torsional curvature at
550 K. Each conditional wrapped Gaussian integrates to one; its strictly
triangular dependence preserves normalization by successive integration.
The product component uses local fractions 0.8 for methyl, 0.7 for phenyl,
and 0.3 for backbone torsions. Both proposals have positive density everywhere.
The potential, retained conformers, Voronoi partition, and physical partition
function are unchanged. Proposal choices and final point counts appear below.
The verifier independently reconstructs every random draw, density, and basin
mask, and numerically checks conditional normalization.

The pilot quadrature used 1024 points per basin. After inspecting
numerical precision, production starts at 4096 points per basin. The earlier
pilot predates the final symmetry partition and is retained as development
history. The final diphenylpropane partition is compared at the base and
increased point counts. Seeds are
fixed per molecule and conformer. An SCC failure on a strained sample is
retried from a fresh electronic state, with more iterations; a looser
intermediate solve, if needed, is always followed by the original accuracy.
The geometry, electronic temperature, and final Hamiltonian/accuracy are
unchanged. Failed points are not silently dropped. Restart counts are
recorded, and an unresolved failure writes its geometry and stops the run. Effective sample sizes and standard errors
are retained; a small point count is not assumed to mean a converged integral.
Standard errors describe this numerical integration only. A low effective
sample size makes asymptotic errors an imperfect convergence diagnostic.
Before calculating corrected Tc, the cycle-level numerical targets are
0.5 kJ/mol in H, 1 J/mol/K in S, and 0.5 kJ/mol in G at each reported
temperature, for both electronic levels. Sample counts are increased per
molecule when those targets are missed; proposal efficiency is also examined.
Neither is selected using Tc. An asymptotic error alone is insufficient when
a few rare importance weights dominate; low effective sample counts and
changes between retained quadratures remain visible in the report.

The implemented surface assigns infinite potential to severe overlaps with any
interatomic distance below 0.45 Å. These points have zero Boltzmann weight and
remain in the Monte Carlo denominator. This is an explicit steric regularization,
not an evaluated GFN2 energy or an SCC failure. The cutoff has been present in
the production integrator from the initial quadratures; documenting it does not
change the sampled Hamiltonian. The verifier reconstructs the geometry of every
accepted zero-weight point and confirms the cutoff accounts for it. Any other
non-finite electronic result is an error. The cutoff is a shared approximation;
representative boundary probes are not a proof about the entire excluded region.
Diphenylpropane receives an 8192-point quadrature with a retained base-count
comparison. The complete correlated preflight still missed the added-unit
targets at 700 and 800 K. Variance contributions identified heterotactic wells
22, 28, 31 and 46 for 32768 points each, before any corrected Tc was calculated.
These allocations improve numerical integration; they do not refine the
electronic Hamiltonian or select a result by agreement with a ceiling.

For either DFT electronic level, conformer probabilities are proportional to
`q_thermal,i(T) exp[−(E_DFT,i−E_DFT,min)/(RT)]`. The ensemble H is the weighted
enthalpy; S includes `−R Σ_i p_i ln p_i` once.
The Shannon term is over integration basins. Proper atom-label images carry
coordinate multiplicity canceled by the external symmetry divisor; they do
not create chemical conformers or an extra configurational entropy term.
Translation and external rotation use the ideal-gas one-bar standard.
The symmetry calculation counts graph automorphisms that preserve the
signed volume at every tetrahedral atom, including formally achiral CH2 and
CH3 atoms, then factors out internal rotor permutations. This excludes
reflections and improper ligand exchanges. A simple ratio of stereo to
constitutional automorphisms fails on the trimers: the RMG non-optical factor
already excludes some improper end exchanges. RMG's optical half factors are
audited separately and are not inserted into these fixed-sequence partitions. A specified enantiomer is calculated; its mirror is not inserted
again as a configurational entropy term. Atactic averages are conditional
averages over the fixed sequence weights above.

**Physical approximation limits:** the kinetic metric and non-torsional
vibrations are evaluated at each minimum; bond-angle relaxation along torsions
is omitted. The Pitzer–Gwinn factor is a semiclassical approximation, not an
exact coupled quantum rotor solution. The two composite levels share xTB
geometries, Hessians and rotor potentials, and share a small DFT basis. Their
spread cannot bound those common errors, basis-set error, intramolecular BSSE,
missed conformers, finite end effects, or transfer from closed-shell fragments
to radical polymer chains. A method-spread envelope is not a physical confidence
interval. These limits are relevant especially for correlated backbone and
phenyl motion at high temperature.

## Measured results

<!-- BEGIN I043:measured -->
| Molecule / stereo sequence | Atoms | CREST candidates | Within 12 kJ/mol | Stereo checked | Chemical minima incl. mirrors | Rotor-space wells | Mirror partners | Proper-permutation images | Lowest nonzero frequency (cm⁻¹) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| ethane | 8 | 2 | 2 | 2 | 1 | 1 | 0 | 0 | 301.840 |
| ethylbenzene | 18 | 3 | 3 | 3 | 1 | 1 | 0 | 0 | 41.140 |
| n_propylbenzene | 21 | 3 | 3 | 3 | 3 | 3 | 1 | 0 | 38.740 |
| diphenylpropane | 31 | 17 | 12 | 17 | 9 | 13 | 2 | 4 | 13.520 |
| diphenylpentane_meso | 37 | 25 | 14 | 25 | 10 | 10 | 2 | 0 | 10.500 |
| diphenylpentane_racemo | 37 | 25 | 13 | 25 | 8 | 13 | 0 | 5 | 7.290 |
| triphenylheptane_iso | 53 | 115 | 31 | 115 | 22 | 22 | 3 | 0 | 8.270 |
| triphenylheptane_syndio | 53 | 110 | 46 | 110 | 31 | 31 | 6 | 0 | 4.300 |
| triphenylheptane_hetero | 53 | 142 | 53 | 142 | 26 | 26 | 0 | 0 | 3.060 |
| cumene | 21 | 1 | 1 | 1 | 1 | 1 | 0 | 0 | 42.880 |

| DFT pilot | Level | Energy (Eh) | Lowest real frequency (cm⁻¹) | Largest imaginary frequency (cm⁻¹) |
| --- | --- | ---: | ---: | ---: |
| ethane_pbe | pbe-D3(BJ)/def2-SVP | -79.645081059 | 309.533 | 0.000 |

Fresh PBE single-point reflection check, diphenylpropane: energy difference -0.000047 kJ/mol (acceptance tolerance 0.02 kJ/mol).

Fresh PBE trimer single-point check with BLAS threading changed: energy difference -0.000000013 kJ/mol (acceptance tolerance 0.02 kJ/mol).

| Reproduced RMG quantity | Value |
| --- | ---: |
| Gas continuous Tc at 1 mol/L (K) | 710.020462 |
| Compiled grid Tc (K) | 710.248673 |
| Propagation ΔH° at 298.15 K (kJ/mol) | -80.039992 |
| Propagation ΔS° at 298.15 K (J/mol/K) | -147.413327 |
| diphenylpropane_isodesmic, ΔH° at 298.15 K (kJ/mol) | -0.000000 |
| diphenylpropane_isodesmic, ΔS° at 298.15 K (J/mol/K) | 11.526306 |
| meso_iso_increment, ΔH° at 298.15 K (kJ/mol) | 0.000000 |
| meso_iso_increment, ΔS° at 298.15 K (J/mol/K) | 0.000000 |
| atactic_added_unit, ΔH° at 298.15 K (kJ/mol) | -0.000000 |
| atactic_added_unit, ΔS° at 298.15 K (J/mol/K) | 0.000000 |

Steric boundary controls evaluate the uncut GFN2 Hamiltonian at retained points nearest the 0.45 Å boundary. These individual points do not bound the whole excluded region.

| Molecule / well | Side | Closest atoms | Distance (Å) | V above minimum (kJ/mol) | ln Boltzmann factor at 800 K | Status |
| --- | --- | --- | ---: | ---: | ---: | --- |
| diphenylpropane / 0 | below | H–H | 0.449485 | 1037.870 | -156.034 | finite |
| diphenylpropane / 0 | above | C–H | 0.452929 | 6110.852 | -918.707 | finite |
| triphenylheptane_iso / 16 | below | H–H | 0.446794 | 3512.042 | -528.001 | finite |
| triphenylheptane_iso / 16 | above | H–H | 0.451219 | 1202.898 | -180.844 | finite |

| Molecule | PBE electronic order, low to high | BLYP electronic order, low to high | PBE minimum E (Eh) | BLYP minimum E (Eh) |
| --- | --- | --- | ---: | ---: |
| ethane | 0 | 0 | -79.643073669 | -79.707750545 |
| ethylbenzene | 0 | 0 | -310.247029691 | -310.515669404 |
| n_propylbenzene | 0, 10000, 1 | 0, 10000, 1 | -349.478468139 | -349.778674059 |
| diphenylpropane | 3, 5, 7, 8, 9, 10009, 100001, 100003, 10, 0, 10000, 100000, 100002 | 5, 3, 7, 8, 0, 10000, 100000, 100002, 9, 10009, 100001, 100003, 10 | -580.084829786 | -580.589689980 |
| diphenylpentane_meso | 6, 3, 2, 7, 10007, 10, 13, 0, 10000, 1 | 0, 10000, 1, 6, 10, 3, 2, 13, 7, 10007 | -658.546220660 | -659.115577716 |
| diphenylpentane_racemo | 3, 0, 4, 6, 100000, 11, 100004, 9, 100002, 10, 100003, 7, 100001 | 0, 3, 4, 6, 100000, 7, 100001, 11, 100004, 9, 100002, 10, 100003 | -658.548890289 | -659.117608497 |
| triphenylheptane_iso | 11, 9, 15, 10015, 8, 19, 16, 13, 10012, 12, 17, 14, 10014, 22, 21, 23, 25, 27, 5, 0, 1, 4 | 11, 5, 0, 9, 15, 10015, 19, 8, 16, 1, 13, 4, 17, 12, 10012, 14, 10014, 22, 21, 23, 25, 27 | -967.614919661 | -968.452985526 |
| triphenylheptane_syndio | 1, 0, 5, 2, 9, 8, 10012, 12, 10011, 11, 10, 26, 27, 16, 10016, 14, 17, 18, 22, 19, 23, 10023, 25, 21, 10021, 41, 42, 38, 35, 33, 10033 | 1, 0, 5, 2, 11, 10011, 8, 9, 16, 10016, 12, 10012, 10, 14, 26, 27, 17, 18, 19, 22, 21, 10021, 23, 10023, 25, 33, 10033, 42, 38, 41, 35 | -967.620322758 | -968.458023477 |
| triphenylheptane_hetero | 20, 21, 28, 25, 22, 37, 38, 33, 41, 6, 9, 3, 0, 12, 16, 8, 31, 52, 50, 40, 46, 49, 48, 45, 30, 35 | 20, 21, 28, 25, 22, 9, 37, 33, 3, 38, 6, 0, 12, 16, 8, 41, 31, 52, 50, 46, 30, 35, 40, 49, 45, 48 | -967.620056633 | -968.457851458 |
| cumene | 0 | 0 | -349.478316956 | -349.778617106 |

### Conditional effect on compiled propagation (main table)

H and S anchors below are at 298.15 K. Tc uses the temperature-dependent correction, one-bar gas thermo and 1 mol/L monomer. Each row is a separate transfer hypothesis; rows are not added together.

| Finding | Electronic level | δΔH° (kJ/mol) | δΔS° (J/mol/K) | Corrected ΔH° (kJ/mol) | Corrected ΔS° (J/mol/K) | Gas Tc (K) | Rotor MC SE of Tc (K) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| diphenylpropane_isodesmic | PBE | -3.701 | -8.448 | -83.741 | -155.862 | 694.948 | 0.979 |
| diphenylpropane_isodesmic | BLYP | -5.360 | -10.118 | -85.400 | -157.531 | 704.227 | 1.050 |
| atactic_added_unit | PBE | -9.874 | -30.161 | -89.914 | -177.574 | 630.316 | 1.659 |
| atactic_added_unit | BLYP | -15.623 | -31.307 | -95.663 | -178.721 | 670.640 | 1.768 |

| Finding | Two-method δΔH° range (kJ/mol) | Two-method δΔS° range (J/mol/K) | Two-method Tc range (K) |
| --- | ---: | ---: | ---: |
| diphenylpropane_isodesmic | -5.360 to -3.701 | -10.118 to -8.448 | 694.948 to 704.227 |
| atactic_added_unit | -15.623 to -9.874 | -31.307 to -30.161 | 630.316 to 670.640 |

### Balanced cycle thermochemistry

| Cycle | Level | T (K) | ΔH° (kJ/mol) | ΔS° (J/mol/K) | Rotor MC SE H (kJ/mol) | Rotor MC SE S (J/mol/K) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| diphenylpropane_isodesmic | PBE | 298.15 | 3.701 | 19.975 | 0.066 | 0.179 |
| diphenylpropane_isodesmic | PBE | 300.00 | 3.692 | 19.943 | 0.067 | 0.180 |
| diphenylpropane_isodesmic | PBE | 400.00 | 3.307 | 18.821 | 0.094 | 0.245 |
| diphenylpropane_isodesmic | PBE | 500.00 | 3.038 | 18.220 | 0.121 | 0.296 |
| diphenylpropane_isodesmic | PBE | 600.00 | 2.809 | 17.801 | 0.143 | 0.328 |
| diphenylpropane_isodesmic | PBE | 700.00 | 2.600 | 17.479 | 0.162 | 0.347 |
| diphenylpropane_isodesmic | PBE | 800.00 | 2.410 | 17.224 | 0.178 | 0.358 |
| diphenylpropane_isodesmic | BLYP | 298.15 | 5.360 | 21.644 | 0.068 | 0.191 |
| diphenylpropane_isodesmic | BLYP | 300.00 | 5.344 | 21.589 | 0.069 | 0.192 |
| diphenylpropane_isodesmic | BLYP | 400.00 | 4.639 | 19.541 | 0.096 | 0.255 |
| diphenylpropane_isodesmic | BLYP | 500.00 | 4.155 | 18.454 | 0.123 | 0.306 |
| diphenylpropane_isodesmic | BLYP | 600.00 | 3.783 | 17.774 | 0.146 | 0.336 |
| diphenylpropane_isodesmic | BLYP | 700.00 | 3.485 | 17.313 | 0.165 | 0.354 |
| diphenylpropane_isodesmic | BLYP | 800.00 | 3.241 | 16.986 | 0.180 | 0.364 |
| atactic_added_unit | PBE | 298.15 | -9.874 | -33.043 | 0.178 | 0.568 |
| atactic_added_unit | PBE | 300.00 | -9.852 | -32.972 | 0.179 | 0.571 |
| atactic_added_unit | PBE | 400.00 | -8.988 | -30.446 | 0.222 | 0.660 |
| atactic_added_unit | PBE | 500.00 | -8.602 | -29.570 | 0.274 | 0.701 |
| atactic_added_unit | PBE | 600.00 | -8.384 | -29.172 | 0.349 | 0.760 |
| atactic_added_unit | PBE | 700.00 | -8.133 | -28.786 | 0.417 | 0.800 |
| atactic_added_unit | PBE | 800.00 | -7.760 | -28.290 | 0.480 | 0.824 |
| atactic_added_unit | BLYP | 298.15 | -15.623 | -34.189 | 0.182 | 0.582 |
| atactic_added_unit | BLYP | 300.00 | -15.604 | -34.124 | 0.183 | 0.584 |
| atactic_added_unit | BLYP | 400.00 | -14.659 | -31.386 | 0.219 | 0.652 |
| atactic_added_unit | BLYP | 500.00 | -14.042 | -29.998 | 0.263 | 0.664 |
| atactic_added_unit | BLYP | 600.00 | -13.598 | -29.185 | 0.338 | 0.719 |
| atactic_added_unit | BLYP | 700.00 | -13.166 | -28.521 | 0.410 | 0.766 |
| atactic_added_unit | BLYP | 800.00 | -12.667 | -27.855 | 0.475 | 0.801 |

Cycle entropy components below are in J/mol/K. Basin mixing includes coordinate images whose multiplicity is canceled by external symmetry. Translation/rotation are capped-molecule contributions and are not automatically transferable local interactions.

| Cycle | Level | T (K) | Translation | External rotation | Non-torsional vibrations | Coupled rotors | Basin mixing | Total raw ΔS° |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| diphenylpropane_isodesmic | PBE | 298.15 | 9.616 | 40.834 | -3.292 | -17.506 | -9.678 | 19.975 |
| diphenylpropane_isodesmic | PBE | 800.00 | 9.616 | 40.551 | -5.086 | -18.694 | -9.162 | 17.224 |
| diphenylpropane_isodesmic | BLYP | 298.15 | 9.616 | 41.260 | -2.462 | -16.074 | -10.696 | 21.644 |
| diphenylpropane_isodesmic | BLYP | 800.00 | 9.616 | 40.829 | -4.542 | -18.640 | -10.276 | 16.986 |
| atactic_added_unit | PBE | 298.15 | -14.072 | -45.298 | 9.867 | 19.882 | -3.422 | -33.043 |
| atactic_added_unit | PBE | 800.00 | -14.072 | -44.782 | 10.123 | 24.805 | -4.365 | -28.290 |
| atactic_added_unit | BLYP | 298.15 | -14.072 | -45.034 | 10.112 | 19.446 | -4.641 | -34.189 |
| atactic_added_unit | BLYP | 800.00 | -14.072 | -44.721 | 10.072 | 24.313 | -3.446 | -27.855 |

Diphenylpropane quadrature comparison: 4096 versus 8192 points per basin on the same final partition. Other species are unchanged; these are first-cycle changes, not new chemistry.

| Level | T (K) | ΔH° change (kJ/mol) | ΔS° change (J/mol/K) | Higher-count MC SE H (kJ/mol) | Higher-count MC SE S (J/mol/K) |
| --- | ---: | ---: | ---: | ---: | ---: |
| PBE | 298.15 | 0.023 | 0.058 | 0.066 | 0.179 |
| PBE | 300.00 | 0.023 | 0.057 | 0.067 | 0.180 |
| PBE | 400.00 | 0.023 | 0.058 | 0.094 | 0.245 |
| PBE | 500.00 | 0.023 | 0.059 | 0.121 | 0.296 |
| PBE | 600.00 | 0.017 | 0.048 | 0.143 | 0.328 |
| PBE | 700.00 | 0.002 | 0.025 | 0.162 | 0.347 |
| PBE | 800.00 | -0.020 | -0.005 | 0.178 | 0.358 |
| BLYP | 298.15 | 0.014 | -0.000 | 0.068 | 0.191 |
| BLYP | 300.00 | 0.013 | -0.002 | 0.069 | 0.192 |
| BLYP | 400.00 | -0.004 | -0.053 | 0.096 | 0.255 |
| BLYP | 500.00 | -0.005 | -0.055 | 0.123 | 0.306 |
| BLYP | 600.00 | -0.000 | -0.046 | 0.146 | 0.336 |
| BLYP | 700.00 | 0.001 | -0.044 | 0.165 | 0.354 |
| BLYP | 800.00 | -0.005 | -0.052 | 0.180 | 0.364 |

| Conditional capped-fragment channel | Level | ΔH°298 (kJ/mol) | ΔS°298 (J/mol/K) | Atactic weight |
| --- | --- | ---: | ---: | ---: |
| meso_to_iso | PBE | -2.873 | -33.926 | 0.25 |
| meso_to_iso | BLYP | -9.545 | -33.090 | 0.25 |
| meso_to_hetero | PBE | -16.088 | -45.160 | 0.25 |
| meso_to_hetero | BLYP | -22.170 | -46.948 | 0.25 |
| racemo_to_syndio | PBE | -11.104 | -27.472 | 0.25 |
| racemo_to_syndio | BLYP | -16.049 | -28.503 | 0.25 |
| racemo_to_hetero | PBE | -9.430 | -25.613 | 0.25 |
| racemo_to_hetero | BLYP | -14.728 | -28.215 | 0.25 |

The inherited corrected-G3MP2 comparison is +7.15 ± 5.58 kJ/mol and 30.4 J/mol/K at 298.15 K.

| DPP cycle level | H minus inherited G3MP2 (kJ/mol) | S minus inherited G3MP2 (J/mol/K) | S minus pinned RMG (J/mol/K) |
| --- | ---: | ---: | ---: |
| PBE | -3.449 | -10.425 | 8.448 |
| BLYP | -1.790 | -8.756 | 10.118 |

### Molecular ensemble H(T), S(T)

H* is H minus the tabulated minimum electronic energy of that molecule at that electronic level; it includes zero-point and thermal energy and is not a formation enthalpy. Full H (J/mol) = minimum E (Eh) × 2625502.602032719 + 1000 H*. Reaction cycles above use full electronic plus thermal H, not differences between H* columns.

Srep uses a specified enantiomer when the molecule is chiral. Srac includes its equal mirror population at the same total gas pressure; it equals Srep for achiral molecules. The fixed-sequence cycles use Srep; including Srac consistently leaves the weighted added-unit residual unchanged.

| Molecule | T (K) | PBE H* (kJ/mol) | PBE Srep (J/mol/K) | PBE Srac (J/mol/K) | BLYP H* (kJ/mol) | BLYP Srep (J/mol/K) | BLYP Srac (J/mol/K) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| ethane | 298.15 | 206.729 | 228.933 | 228.933 | 206.729 | 228.933 | 228.933 |
| ethane | 300.00 | 206.824 | 229.252 | 229.252 | 206.824 | 229.252 | 229.252 |
| ethane | 400.00 | 212.591 | 245.758 | 245.758 | 212.591 | 245.758 | 245.758 |
| ethane | 500.00 | 219.587 | 261.319 | 261.319 | 219.587 | 261.319 | 261.319 |
| ethane | 600.00 | 227.749 | 276.170 | 276.170 | 227.749 | 276.170 | 276.170 |
| ethane | 700.00 | 236.963 | 290.354 | 290.354 | 236.963 | 290.354 | 290.354 |
| ethane | 800.00 | 247.113 | 303.894 | 303.894 | 247.113 | 303.894 | 303.894 |
| ethylbenzene | 298.15 | 425.039 | 359.015 | 359.015 | 425.039 | 359.015 | 359.015 |
| ethylbenzene | 300.00 | 425.279 | 359.818 | 359.818 | 425.279 | 359.818 | 359.818 |
| ethylbenzene | 400.00 | 440.371 | 402.951 | 402.951 | 440.371 | 402.951 | 402.951 |
| ethylbenzene | 500.00 | 459.303 | 445.051 | 445.051 | 459.303 | 445.051 | 445.051 |
| ethylbenzene | 600.00 | 481.485 | 485.410 | 485.410 | 481.485 | 485.410 | 485.410 |
| ethylbenzene | 700.00 | 506.350 | 523.692 | 523.692 | 506.350 | 523.692 | 523.692 |
| ethylbenzene | 800.00 | 533.446 | 559.844 | 559.844 | 533.446 | 559.844 | 559.844 |
| n_propylbenzene | 298.15 | 502.927 | 396.659 | 396.659 | 503.058 | 396.562 | 396.562 |
| n_propylbenzene | 300.00 | 503.206 | 397.592 | 397.592 | 503.338 | 397.496 | 397.496 |
| n_propylbenzene | 400.00 | 520.756 | 447.752 | 447.752 | 520.902 | 447.696 | 447.696 |
| n_propylbenzene | 500.00 | 542.811 | 496.791 | 496.791 | 542.967 | 496.759 | 496.759 |
| n_propylbenzene | 600.00 | 568.671 | 543.845 | 543.845 | 568.834 | 543.825 | 543.825 |
| n_propylbenzene | 700.00 | 597.668 | 588.487 | 588.487 | 597.836 | 588.475 | 588.475 |
| n_propylbenzene | 800.00 | 629.267 | 630.647 | 630.647 | 629.439 | 630.640 | 630.640 |
| diphenylpropane | 298.15 | 723.852 | 506.765 | 506.765 | 724.140 | 504.998 | 504.998 |
| diphenylpropane | 300.00 | 724.285 | 508.214 | 508.214 | 724.581 | 506.472 | 506.472 |
| diphenylpropane | 400.00 | 751.546 | 586.123 | 586.123 | 752.174 | 585.349 | 585.349 |
| diphenylpropane | 500.00 | 785.805 | 662.303 | 662.303 | 786.659 | 662.036 | 662.036 |
| diphenylpropane | 600.00 | 825.914 | 735.284 | 735.284 | 826.918 | 735.292 | 735.292 |
| diphenylpropane | 700.00 | 870.770 | 804.346 | 804.346 | 871.869 | 804.500 | 804.500 |
| diphenylpropane | 800.00 | 919.507 | 869.372 | 869.372 | 920.662 | 869.603 | 869.603 |
| diphenylpentane_meso | 298.15 | 875.259 | 560.182 | 560.182 | 877.633 | 559.384 | 559.384 |
| diphenylpentane_meso | 300.00 | 875.771 | 561.893 | 561.893 | 878.152 | 561.117 | 561.117 |
| diphenylpentane_meso | 400.00 | 908.025 | 654.066 | 654.066 | 910.643 | 653.989 | 653.989 |
| diphenylpentane_meso | 500.00 | 948.655 | 744.408 | 744.408 | 951.377 | 744.565 | 744.565 |
| diphenylpentane_meso | 600.00 | 996.291 | 831.085 | 831.085 | 999.059 | 831.326 | 831.326 |
| diphenylpentane_meso | 700.00 | 1049.620 | 913.190 | 913.190 | 1052.401 | 913.453 | 913.453 |
| diphenylpentane_meso | 800.00 | 1107.600 | 990.549 | 990.549 | 1110.379 | 990.809 | 990.809 |
| diphenylpentane_racemo | 298.15 | 875.610 | 540.635 | 546.398 | 875.523 | 540.651 | 546.414 |
| diphenylpentane_racemo | 300.00 | 876.134 | 542.387 | 548.150 | 876.045 | 542.395 | 548.158 |
| diphenylpentane_racemo | 400.00 | 909.350 | 637.288 | 643.052 | 909.106 | 636.856 | 642.619 |
| diphenylpentane_racemo | 500.00 | 951.283 | 730.534 | 736.297 | 950.944 | 729.883 | 735.646 |
| diphenylpentane_racemo | 600.00 | 1000.161 | 819.480 | 825.243 | 999.852 | 818.881 | 824.644 |
| diphenylpentane_racemo | 700.00 | 1054.491 | 903.133 | 908.896 | 1054.289 | 902.698 | 908.461 |
| diphenylpentane_racemo | 800.00 | 1113.241 | 981.523 | 987.287 | 1113.165 | 981.257 | 987.020 |
| triphenylheptane_iso | 298.15 | 1250.953 | 721.931 | 721.931 | 1250.775 | 721.871 | 721.871 |
| triphenylheptane_iso | 300.00 | 1251.716 | 724.485 | 724.485 | 1251.538 | 724.422 | 724.422 |
| triphenylheptane_iso | 400.00 | 1299.513 | 861.127 | 861.127 | 1299.284 | 860.918 | 860.918 |
| triphenylheptane_iso | 500.00 | 1358.864 | 993.123 | 993.123 | 1358.622 | 992.882 | 992.882 |
| triphenylheptane_iso | 600.00 | 1427.983 | 1118.894 | 1118.894 | 1427.766 | 1118.699 | 1118.699 |
| triphenylheptane_iso | 700.00 | 1505.255 | 1237.863 | 1237.863 | 1505.086 | 1237.741 | 1237.741 |
| triphenylheptane_iso | 800.00 | 1589.325 | 1350.030 | 1350.030 | 1589.208 | 1349.978 | 1349.978 |
| triphenylheptane_syndio | 298.15 | 1250.250 | 708.838 | 708.838 | 1250.056 | 707.726 | 707.726 |
| triphenylheptane_syndio | 300.00 | 1251.030 | 711.448 | 711.448 | 1250.834 | 710.329 | 710.329 |
| triphenylheptane_syndio | 400.00 | 1299.932 | 851.243 | 851.243 | 1299.678 | 849.949 | 849.949 |
| triphenylheptane_syndio | 500.00 | 1360.694 | 986.378 | 986.378 | 1360.459 | 985.126 | 985.126 |
| triphenylheptane_syndio | 600.00 | 1431.298 | 1114.861 | 1114.861 | 1431.148 | 1113.760 | 1113.760 |
| triphenylheptane_syndio | 700.00 | 1509.874 | 1235.844 | 1235.844 | 1509.876 | 1234.976 | 1234.976 |
| triphenylheptane_syndio | 800.00 | 1594.928 | 1349.331 | 1349.331 | 1595.128 | 1348.727 | 1348.727 |
| triphenylheptane_hetero | 298.15 | 1251.224 | 710.697 | 716.460 | 1250.925 | 708.014 | 713.777 |
| triphenylheptane_hetero | 300.00 | 1251.987 | 713.248 | 719.011 | 1251.691 | 710.576 | 716.339 |
| triphenylheptane_hetero | 400.00 | 1299.809 | 849.941 | 855.704 | 1299.840 | 848.193 | 853.956 |
| triphenylheptane_hetero | 500.00 | 1359.687 | 983.092 | 988.855 | 1360.203 | 982.422 | 988.185 |
| triphenylheptane_hetero | 600.00 | 1429.634 | 1110.370 | 1116.133 | 1430.640 | 1110.595 | 1116.358 |
| triphenylheptane_hetero | 700.00 | 1507.753 | 1230.645 | 1236.408 | 1509.150 | 1231.475 | 1237.238 |
| triphenylheptane_hetero | 800.00 | 1592.615 | 1343.871 | 1349.634 | 1594.270 | 1345.048 | 1350.811 |
| cumene | 298.15 | 502.111 | 386.964 | 386.964 | 502.111 | 386.964 | 386.964 |
| cumene | 300.00 | 502.396 | 387.917 | 387.917 | 502.396 | 387.917 | 387.917 |
| cumene | 400.00 | 520.190 | 438.790 | 438.790 | 520.190 | 438.790 | 438.790 |
| cumene | 500.00 | 542.363 | 488.098 | 488.098 | 542.363 | 488.098 | 488.098 |
| cumene | 600.00 | 568.276 | 535.249 | 535.249 | 568.276 | 535.249 | 535.249 |
| cumene | 700.00 | 597.300 | 579.933 | 579.933 | 597.300 | 579.933 | 579.933 |
| cumene | 800.00 | 628.921 | 622.121 | 622.121 | 628.921 | 622.121 | 622.121 |

| Molecule | Points per basin | Importance proposal | External symmetry divisor | Internal rotor divisor | PBE conformer S298 (J/mol/K) | BLYP conformer S298 (J/mol/K) | Lowest basin ESS at 298 K, PBE/BLYP |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| ethane | 4096 | diagonal | 6 | 3 | -0.000 | -0.000 | 3528.4/3528.4 |
| ethylbenzene | 4096 | diagonal | 1 | 6 | -0.000 | -0.000 | 2664.6/2664.6 |
| n_propylbenzene | 4096 | diagonal | 1 | 6 | 8.757 | 8.914 | 2279.3/2279.3 |
| diphenylpropane | 8192 | diagonal | 2 | 4 | 18.436 | 19.609 | 79.4/79.4 |
| diphenylpentane_meso | 4096 | correlated | 1 | 36 | 15.383 | 17.946 | 92.3/92.3 |
| diphenylpentane_racemo | 4096 | correlated | 2 | 36 | 8.814 | 9.317 | 47.0/47.0 |
| triphenylheptane_iso | 4096–16384 | correlated | 1 | 72 | 23.030 | 24.039 | 12.0/12.0 |
| triphenylheptane_syndio | 4096 | correlated | 1 | 72 | 16.543 | 16.103 | 2.9/2.9 |
| triphenylheptane_hetero | 4096–32768 | correlated | 1 | 72 | 15.082 | 15.737 | 12.8/12.8 |
| cumene | 4096 | diagonal | 1 | 18 | -0.000 | -0.000 | 2694.0/2694.0 |

| Steric zero-weight points among selected proposal draws | Count | All draws |
| --- | ---: | ---: |
| ethane | 0 | 4096 |
| ethylbenzene | 0 | 4096 |
| n_propylbenzene | 14 | 12288 |
| diphenylpropane | 86 | 106496 |
| diphenylpentane_meso | 146 | 40960 |
| diphenylpentane_racemo | 137 | 53248 |
| triphenylheptane_iso | 542 | 102400 |
| triphenylheptane_syndio | 716 | 126976 |
| triphenylheptane_hetero | 1437 | 221184 |
| cumene | 0 | 4096 |

| Extra point allocation | Well index | Points |
| --- | ---: | ---: |
| triphenylheptane_iso | 16 | 16384 |
| triphenylheptane_hetero | 22 | 32768 |
| triphenylheptane_hetero | 28 | 32768 |
| triphenylheptane_hetero | 31 | 32768 |
| triphenylheptane_hetero | 46 | 32768 |

Same-proposal convergence: 4096 points in every well versus the selected allocations. The physical partitions and random-stream prefixes are identical.

| Cycle | Level | Maximum absolute change in H over 298–800 K (kJ/mol) | In S (J/mol/K) | In G (kJ/mol) |
| --- | --- | ---: | ---: | ---: |
| diphenylpropane_isodesmic | pbe | 0.023 | 0.059 | 0.017 |
| diphenylpropane_isodesmic | blyp | 0.014 | 0.055 | 0.037 |
| atactic_added_unit | pbe | 0.642 | 1.010 | 0.166 |
| atactic_added_unit | blyp | 0.777 | 1.228 | 0.206 |

R ln 2 = 5.763153 J/mol/K remains in the original compiler once.
The global mirror contribution excluded from both fixed-sequence averages would change the added-unit entropy by 0.000000 J/mol/K if included consistently.
The QM capped molecules' external end-symmetry increment is 2.881576 J/mol/K; the GAV end increment is 5.763153, and its local optical increment is 0.000000. Both sides are normalized before transfer. Their difference adds 2.881576 to the raw QM-minus-GAV entropy residual.
The excluded class-mixing increment is 2.881576 J/mol/K. QM end symmetry + class mixing + global mirror contributions give 5.763153 J/mol/K. GAV end + local optical increments give the same one new-centre bit. No independent R ln 2 is added to the compiler or the fixed-sequence averages.
<!-- END I043:measured -->

Absolute electronic energies in this table are computational pilot results,
not formation enthalpies or molecular free energies. CREST counts are candidate
counts before tight optimization, deduplication, and minima validation.

## Balanced comparisons and stereo accounting

I039's standard state is retained: ideal gas, 100000 Pa; monomer concentration
1000 mol/m³. The unchanged propagation root solves
`ΔH°(T) − T ΔS°(T) − RT ln(c0 RT/p0) = 0`.

The first diagnostic cycle is
`diphenylpropane + ethane → ethylbenzene + n-propylbenzene`.
The pinned additive model's enthalpy and heat-capacity contributions cancel
exactly. Its nonzero reaction entropy is symmetry bookkeeping, not an estimated
torsional interaction. The dispatch's supporting corrected-G3MP2 reference
(+7.15 ± 5.58 kJ/mol; 30.4 J/mol/K) is an inherited comparison input, not a new
QM result reproduced here. It is not used to choose a method or a conformer.

The trimer-minus-dimer electronic energy is an unbalanced molecular increment
unless the added C8H8 fragment has an explicit reference. To isolate the
additive one-ring increment without importing a styrene energy, add a supporting
cumene calculation and use
`dimer + cumene + n-propylbenzene → trimer + ethylbenzene + ethane`.
The script verifies that this cycle also has a zero net additive source vector.
Its residual is the incremental nonadditivity relative to the one-ring
surrogate, subject to cap, symmetry, and closed-shell-to-radical transfer errors.
The contract's nine listed molecules alone do not supply this balanced reference.
The growing-fragment cycle also contains interactions involving the third
phenyl. It does not isolate a uniquely transferable pair potential or establish
the long-chain limit.

For a frozen, unbiased atactic sequence, dyad weights are meso/racemo = 1/2,
1/2; triad weights are mm/rr/mr-or-rm = 1/4, 1/4, 1/2. These weights must
not be replaced by an equilibrium Boltzmann mixture of diastereomers without
changing the represented physical ensemble. Meso molecules and symmetric ends
also require distinguishing oriented stereosequences from distinct molecular
isomers. A molecule's external rotational symmetry number, the multiplicity
of its mirror partner, and the chain's sequence entropy are different factors.
RMG's existing incremental R ln 2 is retained once, as established in I041;
an additional independent R ln 2 must not be inserted into propagation.

The gas-capped fixed-class average contains an external end-exchange
entropy increment of half R ln 2. The dyad-to-triad class-mixing increment
is the other half; the weighted global mirror increment is zero. Together
these recover the one new-centre bit. The GAV capped baseline has a full
R ln 2 end-exchange increment and zero local optical increment: all five
dyad/triad stereo graphs receive two constitutional optical half factors.
Both QM and GAV must be normalized to an oriented fixed sequence before
their difference is transferred. Removing only the QM end term would be
inconsistent. The matched correction is
`δS = S_QM,raw − S_GAV,raw − S_QM,end + S_GAV,end + S_GAV,opt`,
which here is the raw residual plus half R ln 2. Equivalently it is the
QM fixed-class cycle with its class-mixing increment, minus the constitution-
lumped GAV cycle. The original compiler retains its full bit once; no
independent additional R ln 2 is inserted. This normalization was corrected
from a one-sided end subtraction before any corrected Tc was computed.
The report distinguishes the raw capped residual from this matched transfer.
Other finite-cap effects, including translation,
rotation and changes in local geometry, remain part of the transfer limit.

The named I040 and I041 report files are absent in this base worktree. Their
local committed contents were read with:

```bash
git show 5563c0ab78:test/rmgpy/kmc/fixtures/I040_gas_thermo_bench.md
git show 4f283ceb64:test/rmgpy/kmc/fixtures/I041_config_entropy.md
```

## Reproduction commands

Run from the dispatched worktree in bash. Both streams of every computational
job are persisted. Output stays under the named scratch directory.

```bash
cd /home/alon/Code/RMG-Py-kmc-i043-diphenyl-nonadditivity
I043_SCRATCH=/home/alon/runs/i043-diphenyl-nonadditivity
I043_PYTHON=/home/alon/anaconda3/envs/rmg_env/bin/python
I043_QUANTUM_PYTHON=/home/alon/anaconda3/envs/pyscf_env/bin/python
I043_PROBE=$PWD/test/rmgpy/kmc/fixtures/i043_probe
I043_CPUS=0,2,4,6,8,10,12,14
export PYTHONPATH=$PWD
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR=$I043_SCRATCH/mpl
mkdir -p "$I043_SCRATCH/logs"
i043_run() {
  local I043_LOG=$1
  shift
  taskset -c "$I043_CPUS" "$@" > >(taskset -c "$I043_CPUS" tee -a "$I043_SCRATCH/logs/$I043_LOG.stdout.log") \
    2> >(taskset -c "$I043_CPUS" tee -a "$I043_SCRATCH/logs/$I043_LOG.stderr.log" >&2)
}
i043_run build "$I043_PYTHON" setup.py build_ext --inplace
i043_run prepare "$I043_PYTHON" "$I043_PROBE/run_ensembles.py" --prepare-only
i043_run ensembles "$I043_PYTHON" "$I043_PROBE/run_ensembles.py" --hours 18
i043_run baseline "$I043_PYTHON" "$I043_PROBE/baseline.py" --scratch "$I043_SCRATCH/baseline"
# Authoring prefetched these minima alongside the earlier electronic jobs.
i043_run prefetch "$I043_PYTHON" "$I043_PROBE/run_composite.py" --minima-only --hours 4 \
  --species triphenylheptane_syndio triphenylheptane_hetero cumene
i043_run checks "$I043_PYTHON" "$I043_PROBE/check_xtb.py" --species \
  ethane ethylbenzene n_propylbenzene diphenylpropane \
  diphenylpentane_meso diphenylpentane_racemo \
  triphenylheptane_iso triphenylheptane_syndio triphenylheptane_hetero cumene
i043_run composite "$I043_PYTHON" "$I043_PROBE/run_composite.py" --hours 18
i043_run symmetry "$I043_PYTHON" "$I043_PROBE/complete_symmetry.py"
i043_run rotors "$I043_PYTHON" "$I043_PROBE/rotors.py" --samples 4096 --species \
  ethane ethylbenzene n_propylbenzene diphenylpropane \
  diphenylpentane_meso diphenylpentane_racemo \
  triphenylheptane_iso triphenylheptane_syndio triphenylheptane_hetero cumene
i043_run rotors-dpp "$I043_PYTHON" "$I043_PROBE/rotors.py" --samples 8192 --species diphenylpropane
i043_run rotors-correlated "$I043_PYTHON" "$I043_PROBE/rotors.py" --samples 4096 \
  --proposal correlated --species triphenylheptane_iso diphenylpentane_meso \
  diphenylpentane_racemo triphenylheptane_syndio triphenylheptane_hetero
i043_run rotors-iso16 "$I043_PYTHON" "$I043_PROBE/rotors.py" --samples 16384 \
  --proposal correlated --species triphenylheptane_iso --indices 16
i043_run rotors-hetero-extra "$I043_PYTHON" "$I043_PROBE/rotors.py" --samples 32768 \
  --proposal correlated --species triphenylheptane_hetero --indices 22 28 31 46
i043_run steric-dpp "$I043_PYTHON" "$I043_PROBE/rotors.py" --boundary-probe \
  --samples 8192 --species diphenylpropane --indices 0
i043_run steric-iso "$I043_PYTHON" "$I043_PROBE/rotors.py" --boundary-probe \
  --samples 16384 --proposal correlated --species triphenylheptane_iso --indices 16
I043_THERMO_ARGS=(--samples 4096 --samples-for diphenylpropane=8192
  --samples-in triphenylheptane_iso:16=16384
  --samples-in triphenylheptane_hetero:22=32768
  --samples-in triphenylheptane_hetero:28=32768
  --samples-in triphenylheptane_hetero:31=32768
  --samples-in triphenylheptane_hetero:46=32768
  --proposal-for triphenylheptane_iso=correlated
  --proposal-for diphenylpentane_meso=correlated
  --proposal-for diphenylpentane_racemo=correlated
  --proposal-for triphenylheptane_syndio=correlated
  --proposal-for triphenylheptane_hetero=correlated)
i043_run preflight "$I043_PYTHON" "$I043_PROBE/thermochemistry.py" "${I043_THERMO_ARGS[@]}" --preflight
i043_run thermo "$I043_PYTHON" "$I043_PROBE/thermochemistry.py" "${I043_THERMO_ARGS[@]}"
i043_run pilot-ethane "$I043_QUANTUM_PYTHON" "$I043_PROBE/refine.py" \
  --xyz "$I043_SCRATCH/ensembles/ethane/input.xyz" \
  --output "$I043_SCRATCH/pilot/ethane_pbe.json" --level pbe
i043_run pilot-trimer "$I043_QUANTUM_PYTHON" "$I043_PROBE/refine.py" \
  --xyz "$I043_SCRATCH/ensembles/triphenylheptane_iso/input.xyz" \
  --output "$I043_SCRATCH/pilot/triphenylheptane_iso_pbe.json" --level pbe --opt-only
i043_run reflection-check "$I043_QUANTUM_PYTHON" "$I043_PROBE/refine.py" \
  --xyz "$I043_SCRATCH/xtb_checks/diphenylpropane/10000/xtbopt.xyz" \
  --output "$I043_SCRATCH/reflection_validation/diphenylpropane_pbe.json" \
  --level pbe --single-point
i043_run thread-check "$I043_QUANTUM_PYTHON" "$I043_PROBE/refine.py" \
  --xyz "$I043_SCRATCH/xtb_checks/triphenylheptane_iso/0000/xtbopt.xyz" \
  --output "$I043_SCRATCH/thread_probe/iso_0000_pbe_blas1.json" --level pbe --single-point
# Author the measured block after all producers complete:
i043_run report "$I043_PYTHON" "$I043_PROBE/summarize.py" --scratch "$I043_SCRATCH" --write-report
# Replay every producing CLI and independently check the completed report:
i043_run verify "$I043_PYTHON" "$I043_PROBE/verify_results.py" --replay --samples 4096
```

The search and tight xTB checks reuse completed scratch jobs. Delete nothing
to obtain fresh results: use `--scratch` with a new subdirectory under the
dispatched scratch path, generate inputs there, and adjust all subsequent paths.
Searches are stochastic, so rerunning a saved-state audit verifies recorded
samples and arithmetic, not reproducibility or completeness of a fresh search.
`summarize.py --write-report` is the authoring command; without that flag it
checks the report tables. The verifier reads the recorded per-molecule quadrature counts and freshly
recalculates the pinned baseline, ethane DFT geometry/Hessian, ethane rotor
quadrature, a diphenylpropane reflection energy, and a trimer single-point
energy. Other heavy calculations
are replayed from retained inputs and audited against their source logs and
Hessians; they are not claimed as independent fresh quantum calculations.
The original trimer cost pilot was launched without `--opt-only` and interrupted
during its Hessian calculation. The command above reproduces its optimization
only; the unfinished Hessian is excluded from verification. `taskset` constrains
the full process before numerical libraries create worker threads. The listed
logical CPUs represent eight distinct physical cores on this workstation.
Computational timings are observations on a shared
workstation and are not deterministic verifier targets.

## Interpretation and remaining scientific work

Both electronic levels give a positive diphenylpropane breakup enthalpy where
the pinned additive groups give zero. Their enthalpies lie inside the inherited
G3MP2 uncertainty band, but their breakup entropies are below its quoted value.
Negating the QM-minus-GAV cycle makes both compiled propagation anchors more
negative. Their opposing effects on the ceiling partly cancel; the net
conditional gas ceiling decreases at both levels.

The frozen atactic growing-fragment cycle gives a stronger enthalpy and entropy
correction and a wider functional spread in the ceiling. It is a distinct
transfer hypothesis, not an extra term to add to the diphenylpropane-cycle
correction. Its residual contains third-ring interactions and finite-fragment
translation, rotation, geometry and end effects. The component table makes
these entropy contributions explicit; subtracting translation alone would not
establish a consistently mass-matched local interaction model.

The final numerical precision gate passes at every tabulated temperature for
both cycles and both levels. The retained convergence comparison shows the
changes from adding points as well as standard errors. This establishes
numerical reproducibility for the defined composite partition, not converged
CREST sampling or a physical confidence interval. Both functionals share the
xTB geometry, frequency and rotor approximations, the small basis and the
steric regularization.

Independent sampling-convergence checks, fuller electronic/rotor refinement,
longer fragments and radical-chain transfer remain scientific follow-ups.
Neither the inherited G3MP2 comparison nor proximity to a ceiling comparator
validates that transfer. No product correction, solvation combination or owner
choice has been made.

## Verifier output

Executed from the dispatched worktree with `PYTHONPATH=$PWD`, one BLAS thread,
and the eight-core affinity listed above:

```bash
/home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/fixtures/i043_probe/verify_results.py \
  --scratch /home/alon/runs/i043-diphenyl-nonadditivity --replay --samples 4096
```

Exit status: 0. Both streams are retained in
`/home/alon/runs/i043-diphenyl-nonadditivity/verify-full.stdout.log` and
`verify-full.stderr.log`. Actual final output:

```text
I043 same-partition quadrature comparison and deterministic energy prefixes verified
I043 all increased-well allocations and same-proposal convergence reproduced
I043 all molecular tables, balanced cycle sums, correction signs and corrected Tc residuals verified
I043 verification: retained pool hashes/stereo and composite source energies verified
I043 H/S thermodynamic derivative identities independently verified for 10 complete ensembles
I043 rotor/vibration mode counts and the Pitzer-Gwinn harmonic reference limit independently verified
I043 full Cartesian Hessian spectra agree with xTB CLI frequencies; importance draws/densities/basin masks independently reconstructed
I043 every producing CLI replayed; pinned baseline, ethane DFT optimization/Hessian, ethane rotor quadrature, a DPP reflection energy and a trimer energy freshly recalculated
I043 retained-pool report reproduction passed; fresh stochastic ensemble identity is not asserted
```

The replay also freshly reproduced the representative steric boundary controls.
Heavy ensembles, original conformer electronic energies and the other rotor
integrals were reused and audited against their retained source data. Full
fresh DFT trimer frequencies and fresh-search convergence are not claimed.
