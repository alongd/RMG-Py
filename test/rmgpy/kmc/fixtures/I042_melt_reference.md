# Styrene propagation in a PS melt: reference and activity probe

Two published PC-SAFT parameter sets give reproducible, positive ceiling shifts
under the dispatch's fixed-volume convention. These are **conditional surrogate
calculations**, not validated radical transfer thermochemistry. Their discrepancy
and their predicted pressure dependence prevent selecting either as a melt
correction. FH and an offline UNIFAC calculation quantify mixing activities but
cannot, by themselves, supply the missing gas-to-condensed reference terms.

**Not retrieved:** a numerical, source-verifiable styrene/PS literature χ(T)
curve; radical-specific PC-SAFT/SAFT/COSMO parameters; numerical styrene/PS
sorption isotherms from the direct historical study; a high-temperature bulk
polymerization equilibrium measurement; or a complete styrene/PS SAFT-γ Mie
parameter set. These gaps are not replaced with fitted constants. The independent
measurements located below check phase partitioning, not the complete radical
propagation free energy. No model is ranked by agreement with a ceiling comparator.

## Reproduction and scope

Repository `/home/alon/Code/RMG-Py-kmc-i042-melt-reference`, local branch
`i042-melt-reference`, starting commit `0bbfe8ea33781b34ec6a29d12512aaa8f73f1e45`.
Database commit `4a12d36fcdc193ede82c8d1ab5c1653495d445bc` is materialized with
the existing I034 `git show` allowlist. The database remains read-only; no tree
regeneration occurs. Only this report and `i042_probe/` are deliverables. No
product-code, thermo, database or test edits, cluster calculations, push, or
owner decision occur. No excluded local dataset is opened or searched for data;
incidental off-scope web search hits are discarded and supply no evidence.

Every command below is run from the named worktree. Logs contain both streams;
the scratch results and reference downloads are outside the repository. Failed
development attempts remain in the earlier scratch logs. `final-*` logs correspond
to the final reproduction. Reproduction needs no network and uses the named Python.

```bash
mkdir -p /home/alon/runs/i042-melt-reference
PYTHONPATH=$PWD /home/alon/anaconda3/envs/rmg_env/bin/python setup.py build_ext --inplace \
  > >(tee -a /home/alon/runs/i042-melt-reference/build-stdout.log) \
  2> >(tee -a /home/alon/runs/i042-melt-reference/build-stderr.log >&2)
PYTHONPATH=$PWD /home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/fixtures/i042_probe/models.py \
  > >(tee -a /home/alon/runs/i042-melt-reference/final-models-stdout.log) \
  2> >(tee -a /home/alon/runs/i042-melt-reference/final-models-stderr.log >&2)
MPLCONFIGDIR=/home/alon/runs/i042-melt-reference/mpl PYTHONPATH=$PWD \
  /home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/fixtures/i042_probe/run_probe.py \
  --scratch /home/alon/runs/i042-melt-reference/final \
  > >(tee -a /home/alon/runs/i042-melt-reference/final-probe-stdout.log) \
  2> >(tee -a /home/alon/runs/i042-melt-reference/final-probe-stderr.log >&2)
PYTHONPATH=$PWD /home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/fixtures/i042_probe/render_tables.py \
  /home/alon/runs/i042-melt-reference/final/results.json \
  > >(tee -a /home/alon/runs/i042-melt-reference/final-render-stdout.log) \
  2> >(tee -a /home/alon/runs/i042-melt-reference/final-render-stderr.log >&2)
MPLCONFIGDIR=/home/alon/runs/i042-melt-reference/mpl PYTHONPATH=$PWD \
  /home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/fixtures/i042_probe/verify_results.py \
  /home/alon/runs/i042-melt-reference/final/results.json \
  --snapshot /home/alon/runs/i042-melt-reference/final/database \
  > >(tee -a /home/alon/runs/i042-melt-reference/final-verify-stdout.log) \
  2> >(tee -a /home/alon/runs/i042-melt-reference/final-verify-stderr.log >&2)
```

`render_tables.py --update` authors the numeric blocks; its default mode checks
them without writing. `models.py` is also an executable external-regression check.
The verifier compares each snapshot file with the actual pinned git object,
re-estimates gas thermo, independently computes concentration equilibria and
LSER algebra, checks PC-SAFT derivatives by finite differences and pressure
constraints, checks UNIFAC Gibbs-Duhem consistency, exercises actual RMG validity
guards, and checks every numeric report block. It verifies the model arithmetic
and the parameter transcription, not the underlying experiments or transferability.

## State, activity and the equilibrium actually computed

The I039 graphs and conventions are reused, with the exact I034 compiler-selected
pair, not a newly chosen reaction:

```text
[CH2]C(CCc1ccccc1)c1ccccc1 + C=Cc1ccccc1
    <=> [CH2]C(CC(CCc1ccccc1)c1ccccc1)c1ccccc1
```

For `Rn + M ⇌ Rn+1`, define the equal-concentration transfer chemical potentials
`u_i = μ_i(melt) − μ_i(ideal gas at the same concentration)` and
`ΔΔG = u_n+1 − u_n − u_M`. For the matched per-end propagation/depropagation
comparison with unit chain-population ratio:

```text
Kc,melt(T, composition) = Kc,gas(T) exp[−ΔΔG/(RT)]
ln(kprop cM / kdeprop) = ln[Kc,gas(T) cM] − ΔΔG/(RT)
Tc: the right-hand side crosses zero
```

The actual monomer concentration is `cM = NM/(NA V)` in mol/m³, as reproduced
from `rmgpy/kmc/ssa.py::record_propensity`. The compiler's default ceiling uses
a fixed monomer concentration specified independently of the SSA's evolving
populations. A closed simulation's free-monomer count changes with reaction
progress, so its instantaneous ceiling condition also changes. For unequal
chain populations, their concentration ratio belongs in the reaction quotient.

Under the dispatch density convention, `V = Neq MW/(ρ NA)` and
`cM = f_free ρ/MW`, where `Neq` counts conserved styrene equivalents, including
free monomer (use MW in kg/mol with ρ in kg/m³). Density specifies total
equivalents, not the free fraction. The
reported fixed total density includes monomer; equal monomer/repeat lattice
volumes are an explicit approximation, not a measured high-temperature EOS.
Numerical constants are RMG's own R and NA; molar mass is calculated from C8H8
using conventional C/H atomic weights. Carrier mass and its size controls are
illustrative mathematical inputs, not an accessed MWD dataset.

<!-- BEGIN I042:definitions -->
| Quantity | Value |
| --- | --- |
| Gas constant (J/mol/K) | 8.314472 |
| Avogadro constant (1/mol) | 6.02214179e+23 |
| Styrene-equivalent molar mass (g/mol) | 104.152000 |
| Fixed total mass density (kg/m³) | 1050.000 |
| Conserved repeat-equivalent concentration (mol/L) | 10.081419 |
| Volume per repeat equivalent (m³) | 1.64712796e-28 |
| Illustrative carrier DP | 960.135187 |
| Concentration reference (mol/m³) | 1000.000 |
| Short / long radical surrogate masses (g/mol) | 209.312000 / 313.464000 |
| Pinned database files | 65 |
| RMG molecular solvent entries | 203 |
| Polystyrene-labelled solvent entries | 0 |
<!-- END I042:definitions -->

The following are prescribed count/volume scenarios, not inferred equilibrium
compositions. Changing concentration gives an effective entropy shift
`R ln(cM/c0)` with zero enthalpy shift at fixed concentration. It does not edit
intrinsic gas thermo. A pure-monomer-reference activity from FH or UNIFAC is
dimensionless and cannot be substituted for a concentration in mol/m³.

<!-- BEGIN I042:concentrations -->
| Free fraction of conserved equivalents | Monomer (mol/L) | Gas Tc (K) | Gas shift (K) | Effective ΔS concentration shift (J/mol/K) | FH χ=0 activity | UNIFAC activity at 700 K |
| --- | --- | --- | --- | --- | --- | --- |
| 0.001000 | 0.010081 | 519.282773 | -190.737689 | -38.222 | 0.002713 | 0.003021 |
| 0.010000 | 0.100814 | 599.129142 | -110.891320 | -19.077 | 0.026885 | 0.029880 |
| 0.099192 | 1.000000 | 710.020462 | 0.000000 | 0.000 | 0.243942 | 0.266002 |
| 0.100000 | 1.008142 | 710.491162 | 0.470700 | 0.067 | 0.245730 | 0.267908 |
| 0.500000 | 5.040710 | 819.630984 | 109.610522 | 13.449 | 0.823931 | 0.845228 |
<!-- END I042:concentrations -->

For reference, the original gas thermo and gas-model equilibrium monomer
concentrations (`1/Kc,gas`) are:

<!-- BEGIN I042:gas -->
| T (K) | Gas ΔH, 1 M (kJ/mol) | Gas ΔS, 1 M (J/mol/K) | Gas ΔG, 1 M (kJ/mol) | Gas equilibrium monomer (mol/L) |
| --- | --- | --- | --- | --- |
| 600 | -73.975 | -104.393 | -11.339 | 0.103011 |
| 700 | -72.423 | -102.002 | -1.021 | 0.839105 |
| 800 | -70.776 | -99.805 | 9.068 | 3.909031 |
<!-- END I042:gas -->

## Main transfer and ceiling table

All roots are continuous crossings, not compiler-grid interpolations. PC-SAFT
uses the dispatch density and the prescribed monomer concentration with trace
short/long radicals in a monodisperse carrier. Linear RMG alternatives use their
unchanged database coefficients. No number here is fitted to a ceiling.

The PC-SAFT H/S columns in this main table are **effective derivatives along the
fixed-density path**, `S_path = −dΔΔG/dT` and `H_path = ΔΔG + T S_path`. The
separate isobaric derivatives below prevent mistaking these for fixed-pressure
transfer enthalpies/entropies. Linear LSER's H/S are the anchor coefficients by
construction, with zero transfer heat-capacity correction.

<!-- BEGIN I042:main -->
| Model at 1 M | ΔΔG at 600 / 700 / 800 K (kJ/mol) | ΔΔH at 700 K (kJ/mol) | ΔΔS at 700 K (J/mol/K) | Continuous Tc (K) | Shift from gas (K) | Status |
| --- | --- | --- | --- | --- | --- | --- |
| Gas baseline | 0 / 0 / 0 | 0 | 0 | 710.020462 | 0 | gas reference |
| PC-SAFT A | -2.125 / -3.153 / -4.143 | 3.897 | 10.071 | 745.614387 | 35.593925 | radical / T extrapolation |
| PC-SAFT B | -1.922 / -4.610 / -7.161 | 13.670 | 26.114 | 774.568977 | 64.548514 | radical / T extrapolation |
| RMG linear benzene | -2.432 / -3.718 / -5.004 | 5.285 | 12.862 | 753.526132 | 43.505670 | 298 K anchor extrapolation |
| RMG linear dodecane | -4.434 / -6.364 / -8.293 | 7.140 | 19.292 | 790.383770 | 80.363308 | 298 K anchor extrapolation |
| RMG linear hexadecane | -3.949 / -5.519 / -7.088 | 5.467 | 15.694 | 776.532314 | 66.511852 | 298 K anchor extrapolation |
| RMG linear toluene | -4.046 / -5.786 / -7.525 | 6.393 | 17.398 | 781.324909 | 71.304447 | 298 K anchor extrapolation |
<!-- END I042:main -->

These numerical alternatives do not validate a PS melt reference. Choosing a
molecular solvent merely changes an extrapolated LSER surrogate. PC-SAFT
calculates a polymer-specific residual chemical potential, but lacks independent
radical and high-temperature mixture validation. FH/UNIFAC total transfer and Tc
are unavailable for the distinct reason explained below; they are not zero.

## PC-SAFT: complete residual surrogate, with stated closure

Set A is transcribed from [López-Domínguez et al., Can. J. Chem. Eng. (2023),
Table 2 and §3.1](https://doi.org/10.1002/cjce.24908). Styrene parameters were
fitted to vapor pressure; its styrene/PS binary parameter was fitted to a ternary
CO2 phase boundary. The table calls m/M dimensionless and mixes K/°C in the
energy heading; here m/M has units mol/g and ε/k is interpreted in kelvin, as
required by the EOS. These are transcription/unit interpretations, not refits.

Set B is transcribed from [Marijke Aerts, TU/e thesis (2012), Chapter 2,
Table 2, printed p.15](https://doi.org/10.6100/IR723144). Its source is the VLXE
database and cited phase-equilibrium fits; it concerns pretreated polymers.
Neither source supplies parameters for the specific radical ends. The PS values
trace to [Gross & Sadowski, polymer PC-SAFT (2002)](https://doi.org/10.1021/ie010954d),
whose full text was not retrieved in this run; the actual transcriptions above
are the evidence used here.

<!-- BEGIN I042:pc_parameters -->
| Set | m(styrene) | σ(styrene), Å | ε/k(styrene), K | m(PS)/M, mol/g | σ(PS), Å | ε/k(PS), K | k(styrene,PS) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| PC-SAFT A | 2.395496 | 3.780000 | 297.660000 | 0.019000 | 4.107100 | 267.000000 | 0.005000 |
| PC-SAFT B | 3.080000 | 3.712000 | 295.980000 | 0.019000 | 4.107000 | 267.000000 | 0.000000 |
<!-- END I042:pc_parameters -->

The implementation uses the non-associating hard-sphere, chain and dispersion
terms of [Gross & Sadowski (2001)](https://doi.org/10.1021/ie0003887). All universal
coefficients and component regression values are embedded in `models.py`, from
the primary [FeOs dispersion implementation](https://github.com/feos-org/feos-pcsaft/blob/main/src/eos/dispersion.rs),
[hard-sphere implementation](https://github.com/feos-org/feos-pcsaft/blob/main/src/eos/hard_sphere.rs),
[chain implementation](https://github.com/feos-org/feos-pcsaft/blob/main/src/eos/hard_chain.rs)
and [propane parameters](https://github.com/feos-org/feos-pcsaft/blob/main/src/parameters.rs).
No external EOS package or network is needed for reproduction.

For each species, compute `u_i/(RT) = ∂[Ares/(RT V)]/∂c_i` at fixed T,V and all
other concentrations. This is a residual chemical potential; a conventional
fixed-pressure fugacity coefficient would additionally require its `ln Z`
conversion. Oligomer radicals are assigned the PS segment diameter and energy,
with segment count proportional to graph-derived mass, and the same published
monomer/PS cross parameter. No radical-end term, tacticity, or oligomer-specific
fit is invented. This closure is a calculable saturated-segment approximation.

At each state, the effective isobaric transfer derivative holds bath composition
fixed locally. It uses the EOS pressure derivative to change density, while the
ideal-gas concentration reference stays fixed. With `l = ln ρ`,
`S_isobar = −G_T + G_l P_T/P_l`. Neither derivative assumes a linear H−TS law.

<!-- BEGIN I042:pc_transfer -->
| Set | T (K) | ΔΔG (kJ/mol) | H along isochore (kJ/mol) | S along isochore (J/mol/K) | H along isobar (kJ/mol) | S along isobar (J/mol/K) | EOS P (MPa) | Monomer γ relative to ideal gas at equal c | Kc multiplier for pair |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PC-SAFT A | 600 | -2.125 | 4.172 | 10.495 | 2.238 | 7.271 | 106.023 | 3.341487 | 1.531160 |
| PC-SAFT A | 700 | -3.153 | 3.897 | 10.071 | 1.715 | 6.954 | 146.261 | 9.620477 | 1.718865 |
| PC-SAFT A | 800 | -4.143 | 3.662 | 9.756 | 1.213 | 6.695 | 184.283 | 20.055054 | 1.864300 |
| PC-SAFT B | 600 | -1.922 | 14.723 | 27.742 | 7.413 | 15.559 | 129.232 | 5.044102 | 1.470082 |
| PC-SAFT B | 700 | -4.610 | 13.670 | 26.114 | 5.540 | 14.500 | 172.458 | 19.142718 | 2.207955 |
| PC-SAFT B | 800 | -7.161 | 12.825 | 24.983 | 3.808 | 13.712 | 213.272 | 48.224491 | 2.934855 |
<!-- END I042:pc_transfer -->

The monomer γ column is `exp(u_M/RT)` relative to an ideal gas at equal
concentration. It is not a pure-liquid activity. The pair multiplier includes
both radical terms as well as the monomer term. Applying only monomer γ would
give a different, incomplete propagation correction.

If applied consistently to the gas reverse-rate reference, the effective
propagation H/S would become:

<!-- BEGIN I042:corrected -->
| Set | T (K) | Effective ΔH, isochore (kJ/mol) | Effective ΔS, isochore (J/mol/K) | Effective ΔH, isobar (kJ/mol) | Effective ΔS, isobar (J/mol/K) |
| --- | --- | --- | --- | --- | --- |
| PC-SAFT A | 600 | -69.803 | -93.898 | -71.737 | -97.122 |
| PC-SAFT A | 700 | -68.525 | -91.931 | -70.707 | -95.048 |
| PC-SAFT A | 800 | -67.114 | -90.049 | -69.564 | -93.111 |
| PC-SAFT B | 600 | -59.252 | -76.651 | -66.562 | -88.835 |
| PC-SAFT B | 700 | -58.753 | -75.888 | -66.882 | -87.502 |
| PC-SAFT B | 800 | -57.952 | -74.823 | -66.968 | -86.093 |
<!-- END I042:corrected -->

All displayed fixed-density states have positive local mechanical and
carrier/monomer composition stability derivatives in the calculated EOS. This
does not establish global phase stability, real radical properties, or experimental
validity. Their calculated pressures show that holding the supplied density at
high T describes a compressed model state. The following dense-branch pressure
control holds monomer concentration fixed and solves for density at each T.
Its tabulated H/S are local fixed-composition isobaric derivatives at each state,
not a derivative along the entire fixed-monomer-concentration sequence.

<!-- BEGIN I042:isobar -->
| Set | T (K) | EOS density at 1 bar (kg/m³) | ΔΔG (kJ/mol) | Isobar ΔΔH (kJ/mol) | Isobar ΔΔS (J/mol/K) |
| --- | --- | --- | --- | --- | --- |
| PC-SAFT A | 600 | 876.109 | -0.868 | 1.186 | 3.423 |
| PC-SAFT A | 700 | 782.612 | -1.024 | -0.466 | 0.796 |
| PC-SAFT A | 800 | 645.242 | -0.774 | -3.269 | -3.119 |
| PC-SAFT B | 600 | 858.172 | 2.154 | 8.357 | 10.339 |
| PC-SAFT B | 700 | 766.800 | 1.358 | 6.071 | 6.734 |
| PC-SAFT B | 800 | 637.992 | 1.012 | 3.104 | 2.614 |
<!-- END I042:isobar -->

The published vapor-pressure calibration and direct ternary measurement domains
are encoded in `run_probe.py::literature_domains`: the monomer fit spans
243.15–636.15 K; the styrene/PS cross parameter is calibrated at 338.15 K. Thus
the upper part of the requested window extrapolates even the monomer fit, and
the whole window extrapolates the polymer-mixture cross parameter and radical
closure. EOS evaluation has no molecular-solvent critical-temperature guard,
which is useful computationally but supplies no validation.

Prespecified density and carrier-size controls show the state dependence. They
change one named physical/model input; they are neither uncertainty bounds nor
parameter alternatives chosen by comparator distance.

<!-- BEGIN I042:controls -->
| Set | Fixed density (kg/m³) | Carrier DP | ΔΔG at 700 K (kJ/mol) | Tc (K) |
| --- | --- | --- | --- | --- |
| PC-SAFT A | 900 | 960.135 | -1.876 | 731.030573 |
| PC-SAFT A | 1050 | 100.000 | -3.137 | 745.445811 |
| PC-SAFT A | 1050 | 10000.000 | -3.154 | 745.632120 |
| PC-SAFT B | 900 | 960.135 | -0.573 | 719.630691 |
| PC-SAFT B | 1050 | 100.000 | -4.621 | 774.737624 |
| PC-SAFT B | 1050 | 10000.000 | -4.609 | 774.551254 |
<!-- END I042:controls -->

## FH, UNIFAC, SAFT-γ Mie and COSMO alternatives

For an incompressible monomer/PS lattice with `φ = cM/ceq` and a
composition-independent χ(T), FH gives
`ln aM = ln φ + (1−1/DP)(1−φ) + χ(T)(1−φ)^2` relative to pure condensed monomer.
This addresses the activity term, not the pure gas-to-melt reaction reference.
The direct styrene/PS sorption study below is a possible source of χ(T), but its
numeric data and parameter curve were not retrieved. Toluene/PS χ, styrenic
copolymer χ, and viscosity Huggins constants are not substituted for styrene/PS.

To quantify what is available, the FH table uses the **athermal χ=0 limiting
control**, not a literature fit. For the same trace short/long pair, assigning
equal-volume repeat segments gives
`ΔmixG/(RT) = ln(3/2) + ln(ceq/c0) − 1 + χ(T)(2φ−1)` after removing the ideal
concentration quotient. For `χ=Aχ+Bχ/T`, its effective mixing enthalpy is
`R Bχ(2φ−1)` and mixing entropy is
`−R[ln(3/2)+ln(ceq/c0)−1+Aχ(2φ−1)]`. Unknown coefficients remain unknown.
For long polymer chains the finite-chain `ln(3/2)` increment is replaced by
`ln[(n+1)/n]`. Composition-dependent integral χ requires differentiating its
free energy; simply inserting it into this activity law would be unsupported.

Original UNIFAC runs offline here from embedded published groups and pair
coefficients. The source for every R, Q and unlike interaction is the primary
[DDBST published original UNIFAC tables](https://www.ddbst.com/published-parameters-unifac.html).
No fitted styrene/PS parameter is added. Bath repeat counts are
`CH2 + ACCH + 5 ACH`; monomer is `CH2=CH + AC + 5 ACH`. End-saturated short and
long proxies have respectively
`CH3 + CH2 + ACCH2 + ACCH + 10 ACH` and
`CH3 + 2 CH2 + ACCH2 + 2 ACCH + 15 ACH`. Radical thermochemistry is not replaced
with saturated thermo; only this activity model has an explicit missing-radical
closure. The molecular-volume and surface combinatorial term uses coordination
number z=10. A free-volume extension is not supplied or validated.

<!-- BEGIN I042:unifac_parameters -->
| Subgroup | R | Q | Main group |
| --- | --- | --- | --- |
| 1 | 0.9011 | 0.8480 | 1 |
| 2 | 0.6744 | 0.5400 | 1 |
| 3 | 0.4469 | 0.2280 | 1 |
| 5 | 1.3454 | 1.1760 | 2 |
| 9 | 0.5313 | 0.4000 | 3 |
| 10 | 0.3652 | 0.1200 | 3 |
| 12 | 1.0396 | 0.6600 | 4 |
| 13 | 0.8121 | 0.3480 | 4 |
<!-- END I042:unifac_parameters -->

<!-- BEGIN I042:unifac_interactions -->
| Main group i | a(i,1), K | a(i,2), K | a(i,3), K | a(i,4), K |
| --- | --- | --- | --- | --- |
| 1 | 0.000 | 86.020 | 61.130 | 76.500 |
| 2 | -35.360 | 0.000 | 38.810 | 74.150 |
| 3 | -11.120 | 3.446 | 0.000 | 167.000 |
| 4 | -69.700 | -113.600 | -146.800 | 0.000 |
<!-- END I042:unifac_interactions -->

UNIFAC computes mole-fraction γ. The conversion of its excess reaction term to
the fixed concentration standard adds `RT ln(c_molecular,total/c0)` for this
association. The reported monomer activity is `xM γM`. The equivalent χ is the
FH coefficient that gives the same monomer activity at this single composition;
it is a transformation of the UNIFAC prediction, **not measured literature χ(T)**.

<!-- BEGIN I042:mixing -->
| Mixing model | T (K) | ΔmixG, 1 M (kJ/mol) | ΔmixH (kJ/mol) | ΔmixS (J/mol/K) | Pure-monomer activity | χ equivalent | Total transfer / Tc |
| --- | --- | --- | --- | --- | --- | --- | --- |
| FH χ=0 control | 600 | 8.561 | 0.000 | -14.269 | 0.243942 | 0.000000 | reference transfer missing |
| FH χ=0 control | 700 | 9.988 | 0.000 | -14.269 | 0.243942 | 0.000000 | reference transfer missing |
| FH χ=0 control | 800 | 11.415 | 0.000 | -14.269 | 0.243942 | 0.000000 | reference transfer missing |
| Original UNIFAC | 600 | 7.987 | 0.001 | -13.311 | 0.266016 | 0.106754 | reference transfer missing |
| Original UNIFAC | 700 | 9.319 | -0.006 | -13.320 | 0.266002 | 0.106690 | reference transfer missing |
| Original UNIFAC | 800 | 10.651 | -0.010 | -13.326 | 0.265963 | 0.106507 | reference transfer missing |
<!-- END I042:mixing -->

These mixing G/H/S must be added to independently supplied pure-condensed
reaction reference chemical potentials, including the radical-chain increment.
Adding them directly to gas propagation thermo would omit that reference
transfer; adding them to a PC-SAFT residual calculation would double-count its
mixing physics. Therefore no gas-to-melt ΔΔG or total Tc shift is claimed for
FH/UNIFAC. Their numerical χ sensitivities are:

<!-- BEGIN I042:fh_sensitivity -->
| T (K) | ∂ΔmixG/∂χ (kJ/mol) | ∂ΔmixH/∂Bχ (J/mol/K) |
| --- | --- | --- |
| 600 | -3.999 | -6.665 |
| 700 | -4.666 | -6.665 |
| 800 | -5.332 | -6.665 |
<!-- END I042:fh_sensitivity -->

[Jiménez-Serratos et al., Macromolecules (2017)](https://doi.org/10.1021/acs.macromol.6b02072)
provide a SAFT-γ-based coarse-grained atactic PS force field with alkane-like
backbone and toluene-like side branches, validated against melt and alkane-solution
properties. This is a relevant polymer alternative, but not a retrieved
styrene-plus-radical parameterization. The publisher record and author-manuscript
abstract were retrieved; full parameter tables were not available through the
retrieval paths attempted. No numeric SAFT-γ transfer or ceiling is fabricated.
The soft-SAFT literature also contains PS parameters, but these cannot be
transplanted into PC-SAFT: the underlying segment potentials/EOS differ.

[SCM's polymer COSMO-RS tutorial](https://www.scm.com/doc/Tutorials/COSMO-RS/COSMO-RS_polymers.html)
describes polymer central-trimer profiles, a polymer combinatorial correction,
density and average chain mass inputs. In principle this supports an offline
calculation once profiles and software are supplied. The executable probes
`amspython`, `COSMOtherm`, and `cosmors` all return unavailable on this session's
PATH; no usable radical/polymer profile inputs are supplied in the named fixtures.
This is scoped availability evidence, not a machine-wide or literature-wide
claim of absence. No profile-generation/QM job or software installation occurs.
UNIFAC is the demonstrated offline group-contribution alternative.

## RMG solvation alternatives and their limits

The pinned database coefficients below are the parameter source for each RMG
molecular-solvent calculation. Abraham G coefficients and Mintz H coefficients
are used with RMG's own anchor constants. The solvent entry metadata and all
solute descriptors are preserved in the read-only snapshot/results. The common
radical lookup matches the generic node without data and retains saturated
descriptors; the zero radical correction is an estimator fallback, not a
measured cancellation.

<!-- BEGIN I042:rmg_parameters -->
| Solvent | s_g | b_g | e_g | l_g | a_g | c_g | s_h | b_h | e_h | l_h | a_h | c_h |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| benzene | 1.074900 | 0.174920 | -0.325850 | 1.013560 | 0.566830 | 0.115970 | -13.857800 | -2.606620 | 6.442200 | -8.520470 | -10.365510 | -4.568820 |
| dodecane | 0.000020 | -0.216320 | -0.004630 | 0.982610 | 0.000000 | 0.117780 | -2.326550 | 1.767700 | 2.173210 | -9.255440 | 5.823610 | -6.244010 |
| ethylbenzene | 0.991020 | -0.245810 | -0.245760 | 1.004850 | 0.260710 | 0.209890 | unavailable | unavailable | unavailable | unavailable | unavailable | unavailable |
| hexadecane | 0.000000 | -0.026440 | 0.008640 | 0.996360 | 0.016780 | 0.006480 | -1.091090 | -0.030740 | 1.610900 | -9.318610 | -1.950260 | -4.515770 |
| toluene | 0.904490 | 0.118270 | -0.317550 | 1.032170 | 0.587440 | 0.081150 | -13.636300 | -3.514800 | 6.208150 | -8.365460 | -8.104620 | -5.656010 |
<!-- END I042:rmg_parameters -->

RMG's linear correction reproduces the toluene reference but imposes no
high-temperature fluid validity constraint. Its other available option is the
temperature-dependent K-factor model, which requires complete LSER anchors and
a CoolProp fluid with saturation properties, and refuses temperatures at or
above the solvent critical point. It is not a polymer EOS. Ethylbenzene lacks
Mintz coefficients; hexadecane has linear coefficients but no configured
CoolProp fluid name. The probe actually calls these APIs and scans each available
liquid domain for a ceiling crossing.

<!-- BEGIN I042:rmg_validity -->
| Solvent | Molecular critical T (K) | T (K) | K-factor G / H / S (kJ/mol / kJ/mol / J/mol/K) or failure | K-factor Tc (K) |
| --- | --- | --- | --- | --- |
| benzene | 562.020 | 600 | above/equal solvent critical T | no root in liquid domain |
| benzene | 562.020 | 700 | above/equal solvent critical T | no root in liquid domain |
| benzene | 562.020 | 800 | above/equal solvent critical T | no root in liquid domain |
| dodecane | 658.100 | 600 | -1.365 / -3.274 / -3.183 | no root in liquid domain |
| dodecane | 658.100 | 700 | above/equal solvent critical T | no root in liquid domain |
| dodecane | 658.100 | 800 | above/equal solvent critical T | no root in liquid domain |
| ethylbenzene | 617.120 | all | missing enthalpy coefficients | unavailable |
| hexadecane | unavailable | 600 | no CoolProp fluid name | unavailable |
| hexadecane | unavailable | 700 | no CoolProp fluid name | unavailable |
| hexadecane | unavailable | 800 | no CoolProp fluid name | unavailable |
| toluene | 591.750 | 600 | above/equal solvent critical T | no root in liquid domain |
| toluene | 591.750 | 700 | above/equal solvent critical T | no root in liquid domain |
| toluene | 591.750 | 800 | above/equal solvent critical T | no root in liquid domain |
<!-- END I042:rmg_validity -->

The one computed high-temperature K-factor row is an equal-concentration
correction differentiated along its saturation path, not a fixed-pressure melt
enthalpy. None of these K-factor models supplies a valid ceiling in its molecular
liquid domain. No polystyrene-labelled solvent entry is found in the inspected
solvent library; this does not exclude all imaginable polymer models in RMG.

## Independent equilibrium checks located

| Source | Evidence retrieved | What it checks / limitation |
| --- | --- | --- |
| [Miltz and Rosen-Doddy, styrene sorption on PS, European Polymer Journal (1986)](https://doi.org/10.1016/0014-3057(86)90200-4) | Publisher indexed abstract: low-concentration inverse-gas-chromatography isotherms are linear; sorption H/S differ between below and above glass transition. Full numeric isotherms were unavailable. | Direct monomer/PS interaction evidence. Suitable to constrain dilute monomer chemical potential and its T derivative once original numbers/standard are obtained. It does not determine radical-chain transfer. |
| [Görnert and Sadowski, ATR-FTIR phase equilibrium (2007)](https://doi.org/10.1002/masy.200751328) | Publisher abstract reports styrene/PS/CO2 coexistence at 338.15 K, at 10 and 15 MPa, with 6 and 105 kg/mol PS, compared to cloud points. No curves are digitized here. | Direct mixture equilibrium evidence and source for the fitted A binary interaction. CO2 and much lower T distinguish it from the requested neat melt; fitting data are not independent validation of A. |
| [Aerts thesis, Chapter 3](https://doi.org/10.6100/IR723144) | Reports measured residual-monomer partitioning differing from its PC-SAFT predictions for untreated latex products. | Independent warning against assuming parameter-set transferability between polymer preparation states. CO2-swollen latex is not the neat high-T melt. No measured ratio is used as a fitted parameter here. |
| I039 thermochemistry/standard-state diagnosis, `I039_tc_gap_probe.md` | Existing phase-specific polymerization calorimetry/entropy references are read for conventions. | These constrain condensed reference thermo; a compiled ceiling or a textbook equilibrium calculation is not a new equilibrium measurement. No high-T bulk styrene/PS propagation-equilibrium dataset was retrieved. |

The original direct sorption article's missing numeric data are consequential:
the present work cannot compare its predicted monomer transfer magnitude or slope
against those measurements. The scientific falsification check for either PC-SAFT
surrogate is an independently specified melt state and a measured monomer
chemical potential/temperature derivative, followed by radical-chain-increment
evidence. Agreement of one ceiling cannot validate both pieces separately.

## What this changes, and contract corrections

* Concentration changes the reaction quotient and effective entropy term. The
  conserved repeat concentration is not a free-monomer concentration or activity.
  The state must specify free-monomer count or composition and volume; those are
  not fully specified by density. The supplied compiler ceiling is a particular
  prescribed concentration, not a universal closed-melt equilibrium.
* PC-SAFT changes transfer G, its T curvature and pressure/density response. Its
  resulting effective compiled H/S and Tc are given above, conditional on the
  radical closure. The two retrieved parameter sets produce different corrections
  without any Tc fitting. The isobar control even changes the sign of one set's
  G correction, so a density/pressure mismatch cannot be silently ignored.
* FH/UNIFAC change condensed mixing chemical potentials. They leave the missing
  condensed pure-reference H/S unresolved; hence no complete compiled Tc change
  can be assigned from them alone. χ(T) without its standard, composition range
  and reference phase would not fix that omission.
* The toluene number is the existing diagnostic linear LSER correction. Inspection
  and reproduction of the gas compiler rates shows it is not, by itself, an
  implemented validated PS-melt reverse-rate reference. Molecular critical points
  bound the saturation-based K-factor alternative, not PS melt existence.
* Asking for physically validated transfer H/S and a unique physical Tc for
  **every** listed model exceeds the available inputs. Some models provide only
  mixing terms; others lack numerical parameters or radical-end coverage. Reporting
  unavailable values and runnable conditional surrogates is the reproducible
  result, not evidence that missing literature data do not exist.

No correction is chosen or applied. Further scientific work requires retrieval
of the direct sorption numbers/χ(T), independently specified melt PVT and monomer
composition, and independently justified radical transfer/reference increments.
Those are evidence needs, not requests to choose a model on the owner's behalf.

## Verifier output reproduced

```text
I042 pinned database, source hashes, gas thermo, rates and volume/concentration roots verified
I042 FeOs three component regressions and UNIFAC pure-component limits reproduced
I042 PC-SAFT transfer, isochoric/isobaric derivatives, pressure and roots independently checked
I042 FH/UNIFAC mixing quantities, LSER parameters and K-factor validity guards verified
I042 all 15 numeric report blocks verified
```

The final commit SHA is provided in the worker closeout rather than embedded in
this report, avoiding a self-referential hash.
