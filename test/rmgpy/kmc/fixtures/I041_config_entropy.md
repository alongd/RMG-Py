# Configurational and optical entropy of constitution-lumped PS propagation

For the **unbiased, equally weighted stereosequence model represented by the
current constitution-only compiler**, retain the net **R ln 2 once** in
propagation. It is already in RMG's group-estimated entropy. Adding another
configurational term double counts it; deleting it changes the represented
reaction channels. Freezing the chain's existing stereocentres does not itself
justify that deletion. This conclusion follows from state and channel counting,
without a literature ceiling-temperature comparator.

The gas ceiling remains the reproduced baseline. The suppression result is an
unapplied control, or the result for a different model that permits only one
prescribed new stereochemical channel. Neither number establishes the physical
ceiling of a melt or a tacticity-dependent molecular model. No product code,
thermo, database or test is changed, and no correction is selected for the owner.

## Reproduction, inputs and scope

Worktree: `/home/alon/Code/RMG-Py-kmc-i041-config-entropy`; branch:
`i041-config-entropy`; base: `0bbfe8ea33781b34ec6a29d12512aaa8f73f1e45`.
Read-only database commit: `4a12d36fcdc193ede82c8d1ab5c1653495d445bc`.
The antecedent is [I039_tc_gap_probe.md](I039_tc_gap_probe.md). Its primary-end
species, product-size convention, monomer concentration, standard states and
benzylic-end structural control are reused. The pinned allowlist materializer
and baseline selectors come from [i034_probe/run_probe.py](i034_probe/run_probe.py).
Only git-show snapshots of the permitted database files are loaded; no rate tree
is regenerated. No excluded dataset is opened, and no literature pyrolysis
result is used. Literature retrieval is limited to statistical mechanics,
thermochemistry and polymerization thermodynamics.

<!-- BEGIN I041:constants -->
| Quantity | Value |
| --- | --- |
| R (J/mol/K) | 8.31447200 |
| Gas reference pressure (Pa) | 100000 |
| Monomer concentration (mol/m³) | 1000 |
| Pinned database files | 65 |
| Python / RDKit | 3.9.23 / 2025.03.5 |
<!-- END I041:constants -->

Run these commands from the worktree, in bash. All long-run streams are
persisted. A fresh `reproduce` directory prevents dependence on authoring output.
The scripts need no network; literature references support the interpretation
and are not fitted or used as numerical inputs.

```bash
cd /home/alon/Code/RMG-Py-kmc-i041-config-entropy
mkdir -p /home/alon/runs/i041-config-entropy/build
PYTHONPATH=$PWD /home/alon/anaconda3/envs/rmg_env/bin/python \
  setup.py build_ext --inplace \
  > >(tee -a /home/alon/runs/i041-config-entropy/build/stdout.log) \
  2> >(tee -a /home/alon/runs/i041-config-entropy/build/stderr.log >&2)

mkdir -p /home/alon/runs/i041-config-entropy/reproduce
export PYTHONPATH=$PWD
export PYTHONHASHSEED=0
export PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR=/home/alon/runs/i041-config-entropy/mpl
/home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/fixtures/i041_probe/run_probe.py \
  --scratch /home/alon/runs/i041-config-entropy/reproduce/database \
  --output /home/alon/runs/i041-config-entropy/reproduce/results.json \
  > >(tee -a /home/alon/runs/i041-config-entropy/reproduce/stdout.log) \
  2> >(tee -a /home/alon/runs/i041-config-entropy/reproduce/stderr.log >&2)
/home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/fixtures/i041_probe/render_tables.py \
  /home/alon/runs/i041-config-entropy/reproduce/results.json \
  > >(tee -a /home/alon/runs/i041-config-entropy/reproduce/render-stdout.log) \
  2> >(tee -a /home/alon/runs/i041-config-entropy/reproduce/render-stderr.log >&2)
/home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/fixtures/i041_probe/verify_results.py \
  /home/alon/runs/i041-config-entropy/reproduce/results.json \
  --snapshot /home/alon/runs/i041-config-entropy/reproduce/database \
  > >(tee -a /home/alon/runs/i041-config-entropy/reproduce/verify-stdout.log) \
  2> >(tee -a /home/alon/runs/i041-config-entropy/reproduce/verify-stderr.log >&2)
```

`render_tables.py --update` is the authoring command: with the same results
argument it fills the marked blocks. Its default invocation above only checks.
Every numerical table below comes from the scripts. JSON also preserves exact
graphs, individual atom/bond symmetry factors, numerical resonance-hybrid bond
orders, source-group weights, species H/S and canonical stereoisomer lists.

`verify_results.py` re-estimates the species from saved adjacency lists,
independently calls the uncorrected group estimator, checks individual database
source terms, enumerates stereo assignments with a separate Cartesian algorithm,
calculates Kc and roots from H−TS, and checks the compiled rate ratios. It checks
the database bytes against git show at the pin and product-source hashes, then
compares every marked report block. Its explicit microchannel model tests both
detailed balance and closure of population derivatives for nonuniform fixed
prefixes. These checks reproduce the ideal counting model; they cannot validate
uncomputed stereosequence-dependent molecular energies.

## 1. Exactly what RMG does

Code locations refer to the dispatched, unchanged base sources:

| Stage | Code path | Actual behavior |
| --- | --- | --- |
| Public thermo estimate | `rmgpy/data/thermo.py:1256`, `:1273`, `:1424`, `:1429` | Searches libraries first; the selected participants fall back to GAV. Applies −R ln σ to the returned entropy once. Library/direct-QM thermo already contains its own conventions. |
| Radical HBI | `rmgpy/data/thermo.py:2061`, `:2074`, `:2106`, `:2125`, `:2138`, `:2181` | Saturates the radical, adds the radical group and subtracts the H-atom enthalpy. A library/QM saturated symmetry is undone before radical symmetry is applied; raw GAV has no such symmetry to undo. All selected participants here use GAV. |
| Species symmetry | `rmgpy/species.py:665`, `:675`, `:677`, `:688`, `:694` | Calculates on a resonance hybrid; the docstring's older website/testing description does not describe its actual use by the public thermo estimator. |
| Molecule symmetry | `rmgpy/molecule/molecule.py:2275`, `:2285`, `:2292` | Caches the graph-algorithm result. |
| Tetrahedral chirality | `rmgpy/molecule/symmetry.py:60`, `:67`, `:89`, `:100` | Splits ligands after removing an acyclic atom, compares their constitutional graphs, and returns a factor of 0.5 for four single bonds to four different ligands. It does not enumerate R/S diastereomers or assign their energies. |
| Other factors | `rmgpy/molecule/symmetry.py:126`, `:151`, `:223`, `:546` | Multiplies acyclic atom, acyclic bond, axis and cyclic factors. A primary radical with two equivalent H ligands contributes an atom factor of 2. |
| Equilibrium constant | `rmgpy/reaction.py:814`, `:843`, `:891` | Forms exp(−ΔGp°/RT), then multiplies by (p°/RT)^Δν for Kc. |
| Reverse compilation | `rmgpy/kmc/compiler.py:1466`, `:1478`, `:1482` | Retrieves each participant's thermo and divides the family forward rate by its gas Kc. |
| Separate statmech convention | `rmgpy/statmech/conformer.pyx:132`, `:142`, `:166`, `:171` | A specified conformer's optical-isomer multiplicity multiplies its partition function and adds R ln multiplicity to S. This path is not run by the GAV baseline. |

The implementation and [RMG's symmetry/chirality documentation](https://reactionmechanismgenerator.github.io/RMG-Py/users/rmg/thermo.html#symmetry-and-chirality)
agree on local half factors. RMG's σ includes internal as well as external
symmetry; it is not simply the external rotational point-group number of one
optimized conformer. The local optical algorithm approximates a whole
constitution-lumped stereochemical ensemble. A conformer's enantiomer multiplier
in Arkane is a different object and does not enumerate that ensemble's
diastereomeric sequences.

### Exact proxy counts and contributions

The selected reaction is P(L−1)• + styrene ⇌ P(L)•. Its radical is terminal
primary carbon. H capping leaves a terminal CH2Ph unit without an asymmetric
centre; hence the proxy is not one tetrahedral centre per nominal repeat unit.
Propagation nevertheless adds one half factor. The exact participants are:

<!-- BEGIN I041:species -->
| L | Side | SMILES | Atom product | Bond product | Axis | Cyclic | RMG σ | Half factors m | σ without optical | Enumerated stereoisomers g |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3 | reactant | `[CH2]C(CCc1ccccc1)c1ccccc1` | 1.000000 | 1.000000 | 1.000000 | 4.000000 | 4.000000 | 1 | 8.000000 | 2 |
| 3 | reactant | `C=Cc1ccccc1` | 1.000000 | 1.000000 | 1.000000 | 2.000000 | 2.000000 | 0 | 2.000000 | 1 |
| 3 | product | `[CH2]C(CC(CCc1ccccc1)c1ccccc1)c1ccccc1` | 0.500000 | 1.000000 | 1.000000 | 8.000000 | 4.000000 | 2 | 16.000000 | 4 |
| 4 | reactant | `[CH2]C(CC(CCc1ccccc1)c1ccccc1)c1ccccc1` | 0.500000 | 1.000000 | 1.000000 | 8.000000 | 4.000000 | 2 | 16.000000 | 4 |
| 4 | reactant | `C=Cc1ccccc1` | 1.000000 | 1.000000 | 1.000000 | 2.000000 | 2.000000 | 0 | 2.000000 | 1 |
| 4 | product | `[CH2]C(CC(CC(CCc1ccccc1)c1ccccc1)c1ccccc1)c1ccccc1` | 0.250000 | 1.000000 | 1.000000 | 16.000000 | 4.000000 | 3 | 32.000000 | 8 |
| 5 | reactant | `[CH2]C(CC(CC(CCc1ccccc1)c1ccccc1)c1ccccc1)c1ccccc1` | 0.250000 | 1.000000 | 1.000000 | 16.000000 | 4.000000 | 3 | 32.000000 | 8 |
| 5 | reactant | `C=Cc1ccccc1` | 1.000000 | 1.000000 | 1.000000 | 2.000000 | 2.000000 | 0 | 2.000000 | 1 |
| 5 | product | `[CH2]C(CC(CC(CC(CCc1ccccc1)c1ccccc1)c1ccccc1)c1ccccc1)c1ccccc1` | 0.125000 | 1.000000 | 1.000000 | 32.000000 | 4.000000 | 4 | 64.000000 | 16 |
<!-- END I041:species -->

For the primary-end chain, the atom product is the terminal-radical factor
times the local half factors. Each phenyl ring contributes a cyclic factor;
the bond and axis factors are unity in these participants. Removing the half
factors gives σ_nonopt. Independent stereo enumeration agrees with the full
assignment count for these asymmetric-ended chains. Explicit hydrogens are used
in both enumerators: otherwise RDKit can fail to distinguish CH2• and CH3
ligands on an adjacent centre. The radical's presence and hydrogen count are
preserved; no isotope or other physical alteration is made.

For L=3, S_sym = −R ln σ_nonopt + m R ln 2:

<!-- BEGIN I041:species_entropy -->
| L=3 side | S without optical symmetry (J/mol/K) | S optical (J/mol/K) | Total symmetry S (J/mol/K) |
| --- | --- | --- | --- |
| reactant chain | -17.289458 | 5.763153 | -11.526306 |
| styrene | -5.763153 | 0.000000 | -5.763153 |
| product chain | -23.052611 | 11.526306 | -11.526306 |
<!-- END I041:species_entropy -->

The reaction's nonoptical contribution cancels exactly: its product symmetry
equals the product of the chain and monomer nonoptical symmetries. The optical
contribution is the new centre's R ln 2. The HBI radical and remote-group terms
are included in the source sum, not added again as entropy corrections.

## 2. Entropy bookkeeping from species and channels

### A constitution lump and its stereosequences

Let α label a stereochemical sequence of an oriented chain, and let qα(T)
include that sequence's conformers, internal rotations, vibrations, electronic
states and appropriate symmetry divisions. Constitution-lumping gives

```text
Q_L(T) = Σα q_L,α(T)
μ°_L,lump = −RT ln Q_L + common translational/reference terms
```

The common terms must be kept consistently in reactants and products. If all
g_L sequences have equal single-sequence free energy, Q_L = g_L q_L,single,
so μ°_L,lump = μ°_L,single − RT ln g_L. Therefore

```text
K_lump / K_single = g_L / (g_(L−1) g_M)
ΔS_stereo = R ln[g_L / (g_(L−1) g_M)]
ΔH_stereo = 0, ΔCp_stereo = 0                 [constant equal weights]
```

For the asymmetric proxies, independent binary assignments give g=2^m and
g_M=1. The new centre supplies R ln 2. This is an entropy of stereochemical
**sequences**, not an extra rotational conformer and not a further vibrational
mode.

For a chain with at least one asymmetric centre, the full assignment count can
also be written as a global mirror-pair factor times the number of relative
stereosequences. The global mirror factor is the same before and after
propagation and cancels. The additional relative-sequence choice remains.
Calling the local RMG term “optical” must not lead to removing it because a
whole-chain enantiomer factor would cancel. In this proxy it represents that
relative-sequence increment too. For symmetric ends, rings, meso structures or
correlated stereocentres, 2^m need not be the number of distinct isomers: use the
actual partition sum, not this local rule.

The more general relative factor f(T) between a full stereochemical sum and a
specified single-sequence reference gives

```text
δG = −RT ln f(T)
δS = R ln f(T) + RT d ln f(T)/dT
δH = RT² d ln f(T)/dT
```

For a reaction, f is the ratio of product and reactant relative factors. Unequal
diastereomer/conformer energies can therefore change H, S and Cp as well as the
constant entropy. The current GAV optical factor computes none of those energy
differences. [Flory's stereochemical-equilibrium paper](https://doi.org/10.1021/ja00984a007)
explicitly builds the equilibrium sum from the conformer partition functions
of the stereoisomers and includes neighbor effects; it does not support assuming
that arbitrary vinyl-polymer diastereomers are energetically identical. His
*Statistical Mechanics of Chain Molecules* (Interscience, 1969) is the relevant
RIS framework. The retrieved [Flory Nobel lecture](https://www.nobelprize.org/uploads/2018/06/flory-lecture.pdf)
also identifies bond rotations and their correlations as the basis of chain
conformational averages. These motivate a missing calculation, not a fitted
entropy correction.

### Atactic growth and fixed centres during thermal chain reactions

For random, unbiased creation of a new centre, the conditional sequence entropy
is

```text
h_stereo = −R ⟨Σs p(s | existing history) ln p(s | existing history)⟩
         = R ln 2                            [independent equal probabilities]
```

This is the incremental entropy; the entropy of the existing history is not
added repeatedly. For correlated or biased growth, the conditional entropy is
smaller. A nonequilibrium growth probability is not automatically a Boltzmann
weight or an equilibrium Kc. Sequence-dependent free energies and rates must
be specified before replacing the equilibrium sum with such probabilities.
The [Yoshida–Aoki–Sasanuma polymer-statistics study](https://doi.org/10.1021/acs.macromol.0c02063)
treats conformational and tacticity entropies separately using RIS and stochastic
sequence descriptions. Its material is not PS, and no numerical parameter from
it is transferred here.

Now hold every old centre fixed. For each prefix α, propagation has daughters
αR and αS. Let each addition channel have rate a and each daughter have its
unique reverse rate b, with equal per-channel K_single = a/b. Then

```text
Pα + M ⇌ PαR       a / b = K_single
Pα + M ⇌ PαS       a / b = K_single
total forward = 2a; reverse per chain = b
K_lump = 2a/b = 2 K_single
```

At equilibrium c(αs)=K_single c(α)c(M), so summing over either daughter and any
population of fixed prefixes gives the same lumped K. At arbitrary populations,
the summed derivative is 2a c(M)Σαc(α) − bΣα,s c(αs), so lumping closes exactly
under the specified sequence-independent rates. There is no requirement that
an old centre racemize while the chain exists. Depropagation removes its actual
last unit through one channel; subsequent propagation can choose either
daughter again. The Verifier checks this with nonuniform fixed-prefix
populations as well as equilibrium fluxes.

Thus **a fixed old stereochemical history and a single permitted future growth
channel are different constraints**. If both daughter channels remain available,
old frozen centres do not halve K_lump. If the entire new assignment is
prescribed, use K_single together with a single-channel forward rate a, not the
pooled forward 2a. Suppressing the entropy while leaving the pooled forward
unchanged instead doubles the reverse rate. The control table below explicitly
labels that fixed-forward arithmetic.

For reaction of a specific pre-existing chain during pyrolysis, retain each
unchanged stereochemical label along the trajectory. It provides no new local
degeneracy merely because it is omitted from a data structure. Reactions that
destroy a centre, generate one from a planar radical, or join chains need the
appropriate conditional channels. Treat a quenched stereo distribution as
conserved species/population labels when their rates differ; it cannot in general
be replaced by an equilibrium partition sum with a single memoryless rate.
Here the equal-rate binary model is exactly lumpable; real PS tacticity-dependent
rate equality has not been demonstrated. No pyrolysis observation is used to
choose between these descriptions.

### Does measured ΔS_p already include configuration?

An equilibrium ΔS_p for styrene propagating to a **specified atactic ensemble**
is the temperature derivative of its total propagation free energy. It already
contains the accessible conformational and sequence contributions of that
ensemble and its activity convention. One must not append R ln 2 to an already
complete atactic propagation entropy. This is a thermodynamic deduction, not
evidence that every historically tabulated entropy contains every residual term.
[Dainton and Ivin's original ceiling paper](https://doi.org/10.1038/162705a0)
identifies the propagation step at the prevailing monomer concentration as the
relevant thermodynamic reaction. [Ivin's review](https://doi.org/10.1002/(SICI)1099-0518(20000615)38:12%3C2137::AID-POLA20%3E3.0.CO;2-D)
emphasizes activity, solvent and polymer-state dependence. Neither author
defines a universal entropy independent of the polymer ensemble.

Calorimetry alone is different: S(T)=S(0)+∫Cp/T dT. A temperature-independent
sequence entropy is in the integration constant and is invisible to Cp.
Frozen random sequences can carry residual entropy even when they do not
interconvert calorimetrically. [Temperley's original residual-entropy paper](https://nvlpubs.nist.gov/nistpubs/jres/56/jresv56n2p55_A1b.pdf)
discusses random sequence and other residual contributions separately from
thermal entropy. [Warfield and Petree](https://doi.org/10.1002/pol.1961.1205516208)
derive PS/styrene thermal functions from heat capacities; their original text
expressly anticipates residual entropy in atactic PS. The retrieved
[original-paper text](https://ro.scribd.com/document/271669276/Thermodynamic-Properties-of-Polystyrene-and-Styrene)
uses functions relative to absolute zero. Its calorimetric difference should
not be silently equated with an independently established full stereochemical
entropy.

[Odian, *Principles of Polymerization*, section 3-9b-1 and Table 3-15](https://www.eng.uc.edu/~beaucag/Classes/Properties/Books/George%20Odian%20-%20Principles%20of%20Polymerization-Wiley-Interscience%20%282004%29.pdf)
defines H/S as the monomer-to-polymer-repeat differences appropriate to
propagation, and gives condensed-polymer and monomer-concentration footnotes.
It presents total polymerization quantities, not a labelled single-R/S-chain
entropy from which one should automatically restore R ln 2. However, the
retrieved table does not identify the styrene entry's detailed tacticity,
residual-entropy treatment, or measurement-to-compilation provenance. **That
historical entry's inclusion of the complete stereosequence residual term is
not established here.** A claim of certain inclusion, or of a certainly missing
R ln 2, would exceed the retrieved evidence. None of these literature H/S values
is inserted into the probe or used to select a model root.

### Chain ends, symmetry, conformers and rate degeneracy

Consistent capping and oriented, chemically unequal ends are important: they
avoid an artificial end-exchange identification and let remote end terms cancel
in an incremental propagation reaction. Do not identify a chiral molecule with
its mirror by treating a reflection as a proper rotation. If changing caps
changes equivalence of ligands, recalculate the actual stereo and nonoptical
symmetries. An end-exchange correction is not automatically an extensive term
per repeat. The standard benzylic-end control is reproduced below; its GAV
increment is identical, without implying that its kinetics is identical.

Conformational entropy is already partly represented in empirical group S/Cp
and their non-nearest-neighbor corrections. The exact net source contributions
are preserved in JSON. There is no explicit rotor or RIS calculation for these
oligomers, so no independently justified extra “chain conformational entropy”
number can be added. The missing physical correction is a calculation of
sequence-dependent conformer partition sums, with consistent group references.

Family reaction-path degeneracy counts reaction opportunities and belongs in
rates. Species symmetry and sequence counts belong in equilibrium free energies.
They must be consistent through detailed balance, not multiplied into Kc again
because a forward reaction has multiple paths. The compiled forward is a
constitution-level family estimate; interpreting it as the total rate is the
consistent use of the lumped Kc. This work does not independently validate its
stereoselectivity or absolute kinetic accuracy.

## 3. Terms that belong in the current compiled Kc

Within the stated ideal atactic constitution lump, use the estimated source
H/S/Cp, the actual nonoptical symmetry correction, **one** net R ln 2 for the
new stereochemical choice, and the gas pressure-to-concentration conversion.
No second enantiomer/tacticity term is needed. No solvation or melt activity is
computed in this work.

There is one doublet radical chain on each side of propagation, so its spin
degeneracy increment cancels. Radical-specific thermochemistry is already in
the GAV/HBI source route; a separate spin-entropy addition is not made. Rotors
and vibrations remain part of the source estimate, rather than independently
computed new correction terms.

The main table isolates every added bookkeeping contribution. “Source” includes
the actual GAV/ring/radical/HBI sum before the public symmetry correction.
Entropy columns without an explicit unit are in J/mol/K.

<!-- BEGIN I041:main -->
| T (K) | Source ΔH, 1 bar (kJ/mol) | Source ΔS (J/mol/K) | Nonoptical ΔS | Stereo ΔS | 1 M conversion ΔH (kJ/mol) | 1 M conversion ΔS | Total ΔH, 1 M (kJ/mol) | Total ΔS, 1 M (J/mol/K) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 600 | -78.963438 | -150.978316 | 0.000000 | 5.763153 | 4.988683 | 40.822038 | -73.974755 | -104.393125 |
| 700 | -78.242744 | -149.869224 | 0.000000 | 5.763153 | 5.820130 | 42.103719 | -72.422613 | -102.002352 |
| 800 | -77.427910 | -148.782564 | 0.000000 | 5.763153 | 6.651578 | 43.213962 | -70.776332 | -99.805448 |
<!-- END I041:main -->

For the association Δν=−1, with p° and c° fixed:

```text
α(T) = c°RT/p°
Kc [m³/mol] = exp(−ΔGp°/RT) RT/p°
ΔGc° = ΔGp° − RT ln α
ΔHc° = ΔHp° + RT
ΔSc° = ΔSp° + R(ln α + 1)
ceiling: Kc c(M) = 1
```

The RT enthalpy and the additional R in the entropy are derivatives of the
standard-state conversion. They leave H−TS equal to the direct Kc balance.
Writing only an effective entropy ΔSp°+R ln α while retaining ΔHp° produces
the same crossing, but is not the fixed-concentration-standard H/S pair.
Conversion is bookkeeping for the same gas state; it does not turn gas proxy
thermo into melt thermo.

The constant stereochemical term has zero enthalpy and heat capacity. Its
free-energy and rate-ratio effects are:

<!-- BEGIN I041:optical_energy -->
| T (K) | Stereo ΔH (kJ/mol) | Stereo −TΔS (kJ/mol) | Kc with / without term | Reverse with / without term, fixed total forward | Kc (m³/mol) |
| --- | --- | --- | --- | --- | --- |
| 600 | 0.000000 | -3.457892 | 2.000000 | 0.500000 | 0.00970770959 |
| 700 | 0.000000 | -4.034207 | 2.000000 | 0.500000 | 0.00119174657 |
| 800 | 0.000000 | -4.610522 | 2.000000 | 0.500000 | 0.000255817858 |
<!-- END I041:optical_energy -->

The resulting roots and unapplied controls are computed from the same molecular
thermo and compiler grid; no comparator enters the solve:

<!-- BEGIN I041:cases -->
| Convention/control | ΔΔH (kJ/mol) | ΔΔS (J/mol/K) | Kc / baseline | Reverse / baseline at fixed forward | Continuous gas Tc (K) | Grid gas Tc (K) | Continuous ΔTc (K) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| atactic constitution lump (retain once) | 0.000000 | 0.000000 | 1.000000 | 1.000000 | 710.020462 | 710.248673 | +0.000000 |
| one prescribed stereo channel / suppression control | 0.000000 | -5.763153 | 0.500000 | 2.000000 | 672.121955 | 672.223381 | -37.898507 |
| extra R ln 2 / double-counting control | 0.000000 | 5.763153 | 2.000000 | 0.500000 | 752.854028 | 752.945430 | +42.833565 |
<!-- END I041:cases -->

The first row is the accounting answer for the existing model. The second row
is not a correction licensed by frozen old centres. The third illustrates the
consequence of adding configuration again. The grid uses interpolation of
log rate ratios; its difference from the continuous root is numerical, not an
extra physical correction. The same accounting over longer proxies and the
antecedent's benzylic end is:

<!-- BEGIN I041:sizes -->
| Product L | Net half factors | Nonoptical ΔS (J/mol/K) | Stereo ΔS (J/mol/K) | Continuous gas Tc (K) | Grid gas Tc (K) |
| --- | --- | --- | --- | --- | --- |
| 3 | 1 | 0.000000 | 5.763153 | 710.020462 | 710.248673 |
| 4 | 1 | 0.000000 | 5.763153 | 710.020462 | 710.248673 |
| 5 | 1 | 0.000000 | 5.763153 | 710.020462 | 710.248673 |
| benzylic end control | 1 | 0.000000 | 5.763153 | 710.020462 | not compiled |
<!-- END I041:sizes -->

Consequences for propagation: retaining the existing term changes neither ΔH,
ΔS nor Tc; conditional suppression changes only ΔS, and hence Kc and Tc; adding
it a second time likewise changes only ΔS. The numerical changes are given in
the control table. Nonoptical symmetry contributes no propagation increment
in these graphs. A new, energy-weighted stereo/conformer calculation could
change ΔH, ΔS, Cp and Tc together; their physical changes are currently unknown
and cannot be inferred from these suppression controls.

## 4. Other families

The same accounting applies wherever a compiled reverse uses species Kc.
For independent, equally weighted binary centres, the optical component is
Δm R ln 2; unchanged centres cancel, destroying a centre reverses the sign,
and forming multiple independent centres multiplies the factors. The
temperature-dependent free-energy magnitude for the representative counts is:

<!-- BEGIN I041:general -->
| Δm | Stereo ΔS (J/mol/K) | K multiplier | −TΔS at 600 K (kJ/mol) | at 700 K | at 800 K |
| --- | --- | --- | --- | --- | --- |
| -2 | -11.526306 | 0.250000 | 6.915783 | 8.068414 | 9.221045 |
| -1 | -5.763153 | 0.500000 | 3.457892 | 4.034207 | 4.610522 |
| 0 | 0.000000 | 1.000000 | -0.000000 | -0.000000 | -0.000000 |
| 1 | 5.763153 | 2.000000 | -3.457892 | -4.034207 | -4.610522 |
| 2 | 11.526306 | 4.000000 | -6.915783 | -8.068414 | -9.221045 |
<!-- END I041:general -->

With the corresponding forward held fixed, reverse multipliers are the inverse
of the K multipliers. These are **optical components**, not full family
equilibrium constants. Recalculate nonoptical symmetry and all chemical thermo
when the reaction changes the molecule count or chain ends.

The first two controls below are the exact antecedent I034 selections. A third
abstraction is selected structurally at a backbone CH(Ph): the product radical
has three carbon ligands, one phenyl, so the tetrahedral centre is destroyed.
The homolysis graph is the antecedent's compiled recombination pair in the
dissociation direction. No control is selected by its thermo or ceiling.

<!-- BEGIN I041:families -->
| Selected forward reaction | Constitutional SMILES: reactants → products | Δm | RMG stereo ΔS (J/mol/K) | RMG stereo K multiplier | Enumerated g-product/g-reactant | Nonoptical ΔS (J/mol/K) | Full gas ΔH at 700 K (kJ/mol) | Full gas ΔS at 700 K (J/mol/K) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H_baseline | `[CH2]C(CCc1ccccc1)c1ccccc1` + `C=Cc1ccccc1` → `CC(CCc1ccccc1)c1ccccc1` + `C=Cc1[c]cccc1` | 0 | 0.000000 | 1.000000 | 1.000000 | 2.391925 | 49.676678 | -3.702089 |
| H_chain | `CC(CC(CCc1ccccc1)c1ccccc1)c1ccccc1` + `[CH2]C(CC(CCc1ccccc1)c1ccccc1)c1ccccc1` → `CC(CC(CCc1ccccc1)c1ccccc1)c1ccccc1` + `CC(CC(C[CH]c1ccccc1)c1ccccc1)c1ccccc1` | 0 | 0.000000 | 1.000000 | 1.000000 | -3.371228 | -54.546448 | -35.786774 |
| H_stereocentre | `CC(CC(CCc1ccccc1)c1ccccc1)c1ccccc1` + `[CH2]C(CC(CCc1ccccc1)c1ccccc1)c1ccccc1` → `CC(CC(CCc1ccccc1)c1ccccc1)c1ccccc1` + `CC(C[C](CCc1ccccc1)c1ccccc1)c1ccccc1` | -1 | -5.763153 | 0.500000 | 0.500000 | -3.371228 | -71.582480 | -41.864868 |
| homolysis | `CC(CC(CCc1ccccc1)c1ccccc1)c1ccccc1` → `[CH](CCc1ccccc1)c1ccccc1` + `[CH2]C(C)c1ccccc1` | -1 | -5.763153 | 0.500000 | 0.500000 | -5.763153 | 314.885753 | 154.097177 |
<!-- END I041:families -->

`H_baseline` transfers from an aromatic C–H to the primary radical and preserves
the chain's stereocentres. `H_chain` selects a terminal benzylic CH2 donor;
it also preserves the existing asymmetric centres. Their optical increments
are zero despite nonzero chemical and nonoptical entropy changes. Thus
“benzylic abstraction” alone does not mean a centre was destroyed.
`H_stereocentre` tests that different case directly. Its reverse hydrogenation
restores a centre and has the inverse optical factor under the ideal lump.

For the displayed homolysis, the summed number of stereocentres decreases;
the inverse recombination gains that optical increment. The enumerated count
ratio independently confirms the local counting in these selected graphs.
The nonoptical correction is also nonzero and must remain separate. Joining
two planar, inequivalent radical centres can instead create two independent
centres; identical fragments, end exchange, meso products, correlations or
stereospecific channels can invalidate a naive factor-per-centre rule.

Remote fixed labels cancel only when their conditional channel mapping is
preserved. Cutting/joining constitution lumps can expose correlations between
fragments that a memoryless state does not store. Sequence-dependent rates
would need explicit stereo labels or conditional rates; the simple factor
table is not proof of exact lumpability for every compiled family. The probe
checks representative graphs and the propagation channel model, not a census
or validation of every family and stereochemical state.

## Contract corrections, limits and falsification

- **The antecedent's Odian anchor claim is wrong.** I039 says the reference
  temperature is unspecified. The independently retrieved Table 3-15 heading
  says “Enthalpy and Entropy of Polymerization at 25°C.” Its phase and
  concentration footnotes still matter. The prior H/S root additionally relied
  on assigning the model's Cp continuation to the literature pair; the heading
  correction does not validate that continuation. No historical H/S root is
  recomputed or used to calibrate this report.
- **“Optical isomer entropy” needs a definition.** The observed code counts
  local constitutional half factors. In these chains the net term describes
  a relative stereosequence choice, even though the whole-chain mirror factor
  cancels. It is not a measurement of molecular tacticity free energies.
- **A unique physical answer is not determined by constitution alone.** The
  ideal equally weighted atactic accounting is reproducible, but real
  diastereomer weights, frozen populations, conformer correlations and
  stereoselective channel rates are not specified by the dispatch or estimated
  by this GAV term. A physical gas or melt Tc correction cannot be assigned a
  sign, size or confidence interval from it. There is no basis here to remove
  the term solely because an old chain's centres are fixed.

**Not retrieved:** complete publisher PDFs for the historical Dainton–Ivin
articles, Flory's monograph or stereochemical-equilibrium paper, the complete
Odian PDF, or the measurement/provenance chain behind its styrene entropy. The
Odian table heading, values and footnotes were available as indexed original
book text; Flory's original abstract was available as indexed paper text; the
RMG docs and Ivin review were readable. Temperley's PDF exceeded the retrieval
tool's size limit, but its original sequence-entropy discussion was indexed.
The original Flory lecture was available as indexed text, not as a successfully
opened complete PDF.
Warfield–Petree's publisher abstract and a transcription of the original paper
were available; that is enough to identify the calorimetric/residual distinction,
not to complete an independent audit of their measurements. The modern
polymer-statistics abstract was read; its full paper and correction were not
retrieved, and no detailed numerical equation from it is relied on.

Searched specifically for local RMG optical counting, Flory stereochemical
partition sums and RIS, Dainton–Ivin propagation thermodynamics, random-sequence
residual entropy, and Odian/Warfield–Petree entropy provenance. Primary authored
texts and official documentation support the cited claims; publisher metadata
alone is not treated as proof of an unread derivation. The equations and
microchannel model in this report are explicit derivations under stated
assumptions, not quotations from an unavailable source.

The ideal synthesis would fail physically if independent oligomer calculations
found unequal daughter free energies or if rates depended on the frozen prefix
so that summing microscopic trajectories no longer closes. The appropriate
next scientific work is a sequence-resolved conformer/free-energy and channel
calculation, with matched chain ends and reference states. Such evidence can
change both energy and entropy; increasing the length of the same additive
proxy or deleting a constant by comparison to a ceiling cannot substitute for it.

## Reproduced Verifier output

```text
I041 pinned git-show snapshot and product-source hashes verified
I041 species thermo, source terms, symmetry primitives and independent stereo enumeration verified
I041 L=3,4,5 rate/Kc ratios, concentration conversion and all gas Tc controls reproduced
I041 selected H-abstraction, homolysis and benzylic-end numbers reproduced
I041 fixed-prefix microchannel detailed balance and exact population lumping verified
I041 all 9 numeric report blocks verified
```

The final commit SHA is reported in the worker closeout to avoid embedding a
self-referential hash in this document.
