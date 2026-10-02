# Compiled styrene propagation: kinetics and equilibrium benchmark

The selected compiled channel uses the addition family's **generic default**,
not a styrene-specific rule or a matched training measurement. Its forward rate
is millions of times the commonly cited bulk styrene PLP-SEC benchmark in the
overlapping measured temperature range. The discrepancy remains large when
the benchmark Arrhenius expression is extrapolated to the requested higher
temperatures. Both the prefactor and activation energy contribute.

This finding changes the credibility of the propagation and unzipping
**timescales**. It does not, by itself, change propagation ΔH, ΔS, Kc or the
ceiling temperature: replacing the forward rate while retaining detailed
balance scales the reverse rate by the same factor. No correction, model
choice, new event compilation, or rate-tree regeneration is made here.

**Not retrieved:** the full concentration tables and solvent assignment of the
original styrene equilibrium measurements, a phase-matched radical
depropagation rate measurement, numerical high-temperature ESR/EPR rate points,
or the revised styrene Arrhenius parameters and confidence region from the
later IUPAC reanalysis. The numerical comparisons below distinguish direct
artifact reproduction, published benchmark transcriptions, literature-quoted
equilibrium concentrations and conditional cross-phase calculations. They are
not a validation of melt reverse kinetics.

## Reproduction and scope

Worktree `/home/alon/Code/RMG-Py-kmc-i044-kp-benchmark`, branch
`i044-kp-benchmark`, initial commit `0bbfe8ea33781b34ec6a29d12512aaa8f73f1e45`.
Read-only database commit
`4a12d36fcdc193ede82c8d1ab5c1653495d445bc`.
The supplied event artifact is read in place and its content hash verified.
Only this report and `i044_probe/` are deliverables. The database is materialized
by the existing `i034_probe/run_probe.py` allowlist and `git show` into scratch;
the probe does not invoke the I034 compilation routine. Species adjacency
lists come directly from the artifact. Thermo estimation and the benzylic
structure control follow the I039 diagnosis and conventions.

No excluded local dataset is opened. No raw molecular-weight-distribution
dataset or pyrolysis/TGA/yield result is used. Literature retrieval concerns
published propagation coefficients, uncertainty statements, polymerization
equilibrium and solution effects. PLP-SEC-derived coefficients are the
explicitly requested kinetic evidence; they are not raw SEC distributions.

Run from the worktree root, with both streams persisted for each command:

```bash
mkdir -p /home/alon/runs/i044-kp-benchmark
PYTHONPATH=$PWD /home/alon/anaconda3/envs/rmg_env/bin/python \
  setup.py build_ext --inplace \
  > >(tee -a /home/alon/runs/i044-kp-benchmark/build.stdout.log) \
  2> >(tee -a /home/alon/runs/i044-kp-benchmark/build.stderr.log >&2)
MPLCONFIGDIR=/home/alon/runs/i044-kp-benchmark/mpl PYTHONPATH=$PWD \
  /home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/fixtures/i044_probe/run_probe.py \
  --scratch /home/alon/runs/i044-kp-benchmark/probe \
  --output /home/alon/runs/i044-kp-benchmark/results.json \
  > >(tee -a /home/alon/runs/i044-kp-benchmark/probe.stdout.log) \
  2> >(tee -a /home/alon/runs/i044-kp-benchmark/probe.stderr.log >&2)
PYTHONPATH=$PWD /home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/fixtures/i044_probe/render_tables.py \
  /home/alon/runs/i044-kp-benchmark/results.json \
  > >(tee -a /home/alon/runs/i044-kp-benchmark/render.stdout.log) \
  2> >(tee -a /home/alon/runs/i044-kp-benchmark/render.stderr.log >&2)
MPLCONFIGDIR=/home/alon/runs/i044-kp-benchmark/mpl PYTHONPATH=$PWD \
  /home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/fixtures/i044_probe/verify_results.py \
  /home/alon/runs/i044-kp-benchmark/results.json \
  --scratch /home/alon/runs/i044-kp-benchmark/probe \
  > >(tee -a /home/alon/runs/i044-kp-benchmark/verify.stdout.log) \
  2> >(tee -a /home/alon/runs/i044-kp-benchmark/verify.stderr.log >&2)
```

`render_tables.py --update` was used to author the numerical blocks. Its
default invocation checks them without writing. The verifier independently
parses the pinned rule file, re-estimates thermo, calculates dimensional Kc
from species H/S and the gas reference pressure, evaluates Arrhenius ratios
and table interpolation, exercises the public propensity normalization, and
checks every report block. It checks source/snapshot hashes and both pairs
without regenerating the artifact. Reproduction requires no network; embedded
literature constants reproduce transcription and arithmetic, not historical
experiments or unavailable confidence contours.

<!-- BEGIN I044:provenance -->
| Quantity | Reproduced value |
| --- | --- |
| Artifact records | 14998 |
| Artifact SHA256 | 0883a1292e20708a17ee8dc7960cf5648b354e0c24ce3373dee83813cea85c9d |
| Pinned snapshot files | 65 |
| Pinned snapshot SHA256 | 61e2c7357bb4945dce6e6a1a39eccf216e945a01771e461391e7975165ea29eb |
| R (J/mol/K) | 8.31447200 |
| Artifact-grid crossing at 1 M (K) | 710.248673 |
| Continuous thermo crossing at 1 M (K) | 710.020462 |
<!-- END I044:provenance -->

The artifact records RMG-Py commit
`14bf074b691188cf883e8b8329de35b53953cd47`; it predates the worktree's base.
Its compiler source hash equals this checkout's compiler hash. This is an
existing-artifact audit, not a claim that the entire old event set was
recompiled at the current commit. Product sources used by the probe have
their current hashes persisted in `results.json`.

## Compiled record and RMG source

The artifact's `ps_ceiling_pairs` gives:

| Proxy | Propagation event ID | Depropagation event ID |
| --- | --- | --- |
| shorter | `evt_d65182629cef788d6cd5c3196710354d09059ea75aacf0bc5a6f8ee55f979b30` | `evt_d4a1849ae3044f5d7cca470bb5fd722649308408551ae5630d836ca51b078b38` |
| longer (`end_radical@5`) | `evt_d9cc55110401febd6b70a921bc739e8c891419d8307517e96da6ad6717f9970a` | `evt_06e0bff7b719ebc4005b626ffe11c7b6214591a32e408a46f4a521a12f4f7e39` |

For the shorter pair the exact chemical graphs, expressed as SMILES, are

```text
[CH2]C(CCc1ccccc1)c1ccccc1 + C=Cc1ccccc1
  ⇌ [CH2]C(CC(CCc1ccccc1)c1ccccc1)c1ccccc1
```

Both ends here are **terminal primary radicals**. Normal head-to-tail styrene
chain propagation has a terminal benzylic radical. The comparison is thus a
benchmark of the compiled proxy's usefulness for PS propagation, not an
experimental measurement of this exact molecular reaction.

Both pairs are enabled `R_Addition_MultipleBond` records with template
`Cds-CbH_Cds-HH;CsJ-CsHH`. The alkene node describes styrene; the radical node
describes an alkyl primary end. In the pinned
`input/kinetics/families/R_Addition_MultipleBond/rules.py`, the only entry is
`R_R;YJ`, index `3000`, rank `0`, short description `Default`. Its
ArrheniusEP alpha and temperature exponent are zero, so the E0 parameter is
the effective Ea and there is no reaction-enthalpy correction to this rate.

<!-- BEGIN I044:parameters -->
| Quantity | Recovered compiled source | Common IUPAC bulk benchmark |
| --- | --- | --- |
| A (L/mol/s) | 1e+10 | 4.27e+07 |
| Ea (kJ/mol) | 2.092 | 32.500 |
| Temperature exponent n | 0 | 0 |
| Source range (K) | 300.00–1500.00 | 261.15–363.15 |
| Source rank / matched training entry | 0 / none | PLP-SEC multi-laboratory benchmark |
| Loaded addition rules / training entries | 1 / 2962 | not applicable |
<!-- END I044:parameters -->

`family.get_kinetics(..., return_all_kinetics=False)` applied directly to the
artifact graphs/template returns source `rate rules`, entry `None` and the
comment:

```text
Estimated using template [R_R;YJ] for rate rule [Cds-CbH_Cds-HH;CsJ-CsHH]
Euclidian distance = 6.708203932499369
family: R_Addition_MultipleBond
```

The loaded training depository contains entries, but none matches this
reaction. No training entry contributes to this default's data; its source
file has no experimental reference or training-source attribution. Loading
training data alone does not build or populate missing specific rate rules.
The `None` entry/rank in artifact `rate_source` reflects ancestor estimation;
the default rule's own rank can still be recovered from the database. The
artifact's `uncertainty=UNKNOWN` is appropriate, rather than a numerical
uncertainty inferred from the benchmark discrepancy.

The recovered expression, with A and Ea in the units above, is

`k_fwd(T) = A exp[−Ea/(R T)]`.

The independent check reproduces every forward grid value in both pairs.
Their forward and reverse tables match between the shorter and longer proxy
within numerical precision. This establishes absence of chain-length
dependence in these compiled tables, not absence of physical dependence.

### Per-site degeneracy and SSA convention

Both selected forward and reverse records have `raw_path_degeneracy=1`,
`degeneracy=1`, `ssa_multiplier=1`. The forward is bimolecular, with
`distinct-pair N_A*N_B`, and the reverse is unimolecular. Here the pair-count
symbols mean participant counts; denote Avogadro's constant by NAv below to
avoid confusing it with the radical count.

`compiler.py:_build_linked_family_pair` estimates one direction with the
RMG merged-path degeneracy, then builds the reverse table by division by the
reference-thermo Kc. `_record` stores the degeneracy and the SSA multiplier.
`ssa.py:record_rate` multiplies the interpolated table only by
`ssa_multiplier`; it does not multiply again by `degeneracy`.
`event_record.py:ssa_multiplier_for` distinguishes distinct and identical
reactants; the identical-pair correction does not apply to radical + styrene.

For eligible sites in distinct components, volume V in m³ and rate in
m³/mol/s, the public `record_propensity` evaluates

```text
a_forward = k_fwd × N_radical × N_styrene / (NAv × V)
a_reverse = k_rev × N_eligible_radical_ends
```

Per radical, with concentration cM = N_styrene/(NAv V), the forward frequency
is `k_fwd cM`. Concentrations in mol/L must be converted to mol/m³ for the
compiled SI rate. No repeat-unit-count degeneracy or additional factor from
the two chain ends is inserted. The verifier exercises the public function
on supplied eligible candidates; it does not claim to run a full SSA
trajectory or independently verify production site discovery.

## IUPAC benchmark and main comparison

The commonly cited fit in [IUPAC's 2019 summary, Table 1](https://doi.org/10.1515/pac-2018-1108)
is the bulk, low-conversion styrene PLP-SEC benchmark, with the parameters
and temperature window shown above. The [original multi-laboratory report
by Buback et al. (1995)](https://pure.tue.nl/ws/portalfiles/portal/1333387/617494.pdf)
specifies ambient pressure and an upper measured range extending to 93°C;
the later tabulation rounds the window to 90°C. We use the later, commonly
cited rounded Arrhenius parameters throughout, without claiming greater
precision than their printed values.

<!-- BEGIN I044:uncertainty -->
| Source | Stated uncertainty | Meaning |
| --- | --- | --- |
| Original IUPAC EVM fit (1995) | ±0.5 K, ±10% k_p; 95% joint A/Ea region | Assumed input measurement errors; A and Ea are strongly correlated |
| Interlaboratory reanalysis (2022) | SD ln(k_p at 25°C) = 0.08; SD Ea = 1.4 kJ/mol; correlation = 0.04; 81 studies | Pooled error across monomers; not a styrene-specific extrapolation confidence band |
<!-- END I044:uncertainty -->

The original errors are assumptions used in an error-in-variables fit, not
independent parameter error bars or a guaranteed prediction interval.
The [2022 reanalysis by Beuermann et al.](https://doi.org/10.1039/D2PY00147K)
includes systematic interlaboratory variation and reports larger joint
confidence regions and slightly revised fits. Its abstract provides the
pooled errors displayed above; full revised styrene parameters were not
retrieved. The old rounded fit is used because the dispatch explicitly asks
for the commonly cited benchmark. No confidence interval is assigned to its
long extrapolation or to the compiled default.

### Measured benchmark range

<!-- BEGIN I044:measured -->
| T (K) | Recovered rule (L/mol/s) | SSA table (L/mol/s) | IUPAC Arrhenius (L/mol/s) | Rule / IUPAC | SSA / IUPAC |
| --- | --- | --- | --- | --- | --- |
| 261.15 | 3.81568e+09 | refused | 13.4892 | 2.8287e+08 | refused |
| 273.15 | 3.98065e+09 | refused | 26.0352 | 1.52895e+08 | refused |
| 298.15 | 4.30029e+09 | refused | 86.4332 | 4.97528e+07 | refused |
| 300.00 | 4.32273e+09 | 4.32273e+09 | 93.7113 | 4.61281e+07 | 4.61281e+07 |
| 323.15 | 4.59041e+09 | 4.58884e+09 | 238.324 | 1.92612e+07 | 1.92546e+07 |
| 350.00 | 4.87296e+09 | 4.87296e+09 | 602.794 | 8.08396e+06 | 8.08396e+06 |
| 363.15 | 5.00147e+09 | 4.99736e+09 | 903.236 | 5.53728e+06 | 5.53273e+06 |
| 366.15 | 5.02995e+09 | 5.02618e+09 | 986.511 | 5.09872e+06 | 5.09491e+06 |
<!-- END I044:measured -->

`Recovered rule` is the analytic expression underlying the compiled table.
Below the table's lower boundary the SSA refuses evaluation; the displayed
rule comparison there is an extrapolation of the recovered rule, also below
its stated source range. It is not a runnable compiled rate. Inside the
overlapping window, the SSA columns use its actual interpolation. The final
row shows the original report's upper endpoint using the same rounded fit.

### Higher temperatures and requested extrapolation

<!-- BEGIN I044:high_temperature -->
| T (K) | Recovered rule (L/mol/s) | SSA table (L/mol/s) | IUPAC Arrhenius (L/mol/s) | Rule / IUPAC | SSA / IUPAC |
| --- | --- | --- | --- | --- | --- |
| 393.15 | 5.27301e+09 | 5.27022e+09 | 2053.56 | 2.56775e+06 | 2.56638e+06 |
| 403.15 | 5.35739e+09 | 5.35603e+09 | 2627.91 | 2.03865e+06 | 2.03813e+06 |
| 600.00 | 6.57475e+09 | 6.57475e+09 | 63257.2 | 103937 | 103937 |
| 650.00 | 6.79029e+09 | 6.79029e+09 | 104412 | 65033.7 | 65033.7 |
| 700.00 | 6.98066e+09 | 6.98066e+09 | 160435 | 43510.9 | 43510.9 |
| 750.00 | 7.14995e+09 | 7.14995e+09 | 232795 | 30713.5 | 30713.5 |
| 800.00 | 7.30145e+09 | 7.30145e+09 | 322433 | 22644.9 | 22644.9 |
<!-- END I044:high_temperature -->

Every IUPAC value in this table is **Arrhenius extrapolation**, not a
measurement. The compiled values are within its grid and the default's
nominal source range, which does not establish physical validity of this
generic rule for styrene. No retrieved evidence validates either kinetic
expression throughout the requested high-temperature interval.

For the analytic ratio, the separable contributions are

`k_rule/k_p,IUPAC = (A_rule/A_IUPAC) exp[(Ea_IUPAC − Ea_rule)/(R T)]`.

The actual SSA ratio additionally includes the tabulation interpolation
factor. The difference in temperature dependence is therefore a kinetic
barrier mismatch, not a constant degeneracy or unit-conversion error.

<!-- BEGIN I044:contributions -->
| T (K) | A-factor ratio | Ea exponential ratio | SSA interpolation / rule |
| --- | --- | --- | --- |
| 261.15 | 234.192 | 1.20785e+06 | refused |
| 300.00 | 234.192 | 196967 | 1 |
| 363.15 | 234.192 | 23644.2 | 0.999178 |
| 600.00 | 234.192 | 443.81 | 1 |
| 700.00 | 234.192 | 185.792 | 1 |
| 800.00 | 234.192 | 96.6936 | 1 |
<!-- END I044:contributions -->

### Higher-temperature measurements, chain length and solvent effects

<!-- BEGIN I044:other_evidence -->
| Study | Retrieved numerical evidence | Retrieval / interpretation |
| --- | --- | --- |
| Yamada et al. 1992 (ESR) | 273.15–403.15 K; Ea = 39.7 kJ/mol | Publisher abstract; A and high-T rate values not retrieved |
| Zetterlund et al. 2004 (EPR) | 393.15 K; diffusion effects near 80% conversion | Publisher abstract; no digitized rate values |
| Olaj et al. 2002 (PLP chain length) | 298.15–343.15 K; observed variation 25–35%; extrapolated reduction 40–60%; half-change chain length order 100 | Abstract; short/infinite-chain limits depend on modeling function |
| DMF solution PLP (2000) | 313.15 K, 1.00 M; k_p ≈ 0.75 × bulk | Acetonitrile also reduced k_p, with a minimum at intermediate dilution |
<!-- END I044:other_evidence -->

[Yamada, Kageoka and Otsu (1992)](https://doi.org/10.1007/BF00944835)
measured absolute propagation coefficients using ESR over a wider
temperature interval than the commonly quoted PLP range. Their reported
activation energy differs from the IUPAC value; it must not be combined with
the IUPAC prefactor to fabricate a new fit. The retrieved publisher abstract
does not provide its prefactor or tabulated rate values.
[Zetterlund, Yamauchi and Yamada (2004)](https://doi.org/10.1002/macp.200300148)
studied bulk propagation to high conversion using EPR and near-IR. Its
strong diffusion effects at high conversion make it a different regime
from a low-conversion benchmark. No higher-temperature styrene
homopropagation PLP numerical series was retrieved in this search; this is
a retrieval limit, not a claim that none exists. Neither ESR nor EPR is
relabeled as PLP evidence.

[Olaj et al. (2002)](https://doi.org/10.1021/ma011215b)
reported chain-length dependence for styrene and MMA bulk PLP. Their
abstract distinguishes actual coefficient variation from extrapolated
zero/infinite-chain limits and suggests local monomer displacement by the
chain as an interpretation. Those model-dependent limits cannot be assigned
as a measured correction for the compiled short proxy.
[Beuermann's PLP-method modeling study (2002)](https://doi.org/10.1021/ma020437m)
also shows that changing pulse repetition conditions can produce apparent
coefficient variation with a fixed input propagation coefficient. A
repetition-rate trend alone therefore does not establish physical
chain-length dependence. No distribution dataset or simulation output from
that study is used here; this is its published methodological conclusion.
The [DMF/acetonitrile solution PLP study (2000)](https://doi.org/10.1016/S0014-3057(00)00021-5)
demonstrates solvent and dilution effects directly, including a reduction
in DMF and a nonmonotonic concentration dependence in acetonitrile.

These effects matter when identifying the appropriate target for short
oligomer radicals in solution or bulk, and when local monomer activity
differs from overall concentration. They are far smaller in the retrieved
studies than the compiled default discrepancy. They cannot repair the
radical-end mismatch, justify this default rate, or establish a liquid/melt
rate at the requested extrapolation temperatures. The exact equality of
the two artifact proxy tables is a property of a common default rule.

## Independent polymerization-equilibrium evidence

[Bywater and Worsfold (1962), *Anionic polymerization of styrene
(thermodynamics)*](https://doi.org/10.1002/pol.1962.1205816633), studied
butyllithium-initiated living styrene chains in benzene and cyclohexane
solution over 100–150°C. The [retrieved original abstract, mirrored in a
research-paper record](https://www.researchgate.net/publication/225109668_Thermodynamics_of_polymerization_with_special_emphasis_on_living_polymers),
reports measurable residual equilibrium monomer and different concentrations
in the two solvents, with matching free energies after heats-of-solution
corrections. This is genuine polymerization-equilibrium evidence, but its
full numerical tables and pressure/polymer-concentration conditions were
not accessible here.

[Ivin's 2000 thermodynamics paper](https://doi.org/10.1002/(SICI)1099-0518(20000615)38:12%3C2137::AID-POLA20%3E3.0.CO;2-D)
quotes the styrene equilibrium concentration at the elevated temperature
shown below and cites that original work. It also gives an approximate
room-temperature concentration. The retrieved excerpt does not assign the
elevated-temperature value to either solvent. Accordingly, the first is a
**literature-quoted experimental result with incomplete phase provenance**;
the second is an **approximate literature estimate**, not an independently
established measurement. No pressure, monomer activity coefficient, or
uncertainty is invented for either.

For gas-phase radical oligomers with unit adjacent-chain concentration
ratio, `[M]eq = 1/Kc` in mol/m³. Convert this to mol/L for comparison.
The reverse-rate table must obey `k_rev = k_fwd/Kc` at the grid nodes.

<!-- BEGIN I044:reverse -->
| T (K) | Compiled k_rev (s⁻¹) | Gas Kc (m³/mol) | [M]eq from thermo (mol/L) | [M]eq from SSA (mol/L) | IUPAC k_p / same gas Kc (s⁻¹) |
| --- | --- | --- | --- | --- | --- |
| 300.00 | 100.65 | 42948.2 | 2.32839e-08 | 2.32839e-08 | 2.18196e-06 |
| 350.00 | 9521.82 | 511.767 | 1.95401e-06 | 1.95401e-06 | 0.00117787 |
| 383.15 | 97568.3 | 51.928 | 1.92574e-05 | 1.88263e-05 | 0.0305077 |
| 600.00 | 6.77271e+08 | 0.00970771 | 0.103011 | 0.103011 | 6516.18 |
| 650.00 | 2.17608e+09 | 0.00312043 | 0.320469 | 0.320469 | 33460.7 |
| 700.00 | 5.8575e+09 | 0.00119175 | 0.839105 | 0.839105 | 134621 |
| 750.00 | 1.36909e+10 | 0.000522241 | 1.91483 | 1.91483 | 445762 |
| 800.00 | 2.85416e+10 | 0.000255818 | 3.90903 | 3.90903 | 1.2604e+06 |
<!-- END I044:reverse -->

The column replacing just k_fwd with IUPAC while retaining gas Kc is an
**unapplied diagnostic**. It quantifies the kinetics timescale sensitivity
without replacing thermo. At off-grid temperatures, forward and reverse
tables are each interpolated linearly in ln(k) versus T. Their ratio gives
the SSA's interpolated equilibrium, which need not exactly equal a fresh
continuous thermo calculation. Both are displayed rather than silently
substituting one for the other.

<!-- BEGIN I044:equilibrium -->
| T (K) | Literature [M]eq (mol/L) | Gas-model [M]eq (mol/L) | Model / literature | Conditional ΔG addition (kJ/mol) |
| --- | --- | --- | --- | --- |
| 383.15 | 0.00012 | 1.92574e-05 | 0.160479 | +5.828520 |
| 298.15 | 1e-06 | 1.91986e-08 | 0.0191986 | +9.799125 |
<!-- END I044:equilibrium -->

The gas model has a smaller equilibrium monomer concentration than the
quoted condensed-solution references. The last column is
`R T ln([M]eq,literature/[M]eq,gas)`, the positive reaction-free-energy
addition that **would** reconcile concentrations if the states and chain
activities could first be made identical. It is not an identified gas
thermochemistry error or a correction applied to this model. Differences in
solvent, polymer activity, charged versus radical ends and reference state
remain unresolved. The two references are not a justified common-state
van't Hoff series; no enthalpy/entropy fit is made from them.

A conditional conversion of these concentrations to radical reverse rates
is useful to expose the distinction between equilibrium and kinetics:
`k_rev,conditional = k_p,IUPAC [M]eq,literature`. This assumes that the
living-chain equilibrium transfers to the radical system, and that the
bulk kinetic extrapolation transfers to the solution. Neither assumption
is established. The elevated-temperature benchmark value is outside its
measured window; the room-temperature analytic default is outside the
compiled table. **The resulting rates are not measured depropagation.**

<!-- BEGIN I044:conditional_reverse -->
| T (K) | Model rule k_rev (s⁻¹) | SSA k_rev (s⁻¹) | Conditional k_p,IUPAC × literature [M]eq (s⁻¹) | Model rule / conditional |
| --- | --- | --- | --- | --- |
| 383.15 | 99862.7 | 97568.3 | 0.190105 | 525304 |
| 298.15 | 82.5596 | refused | 8.64332e-05 | 955183 |
<!-- END I044:conditional_reverse -->

No directly measured styrene radical-depropagation coefficient with matched
T, phase and concentration was retrieved from polymerization studies.
Equilibrium evidence constrains a rate **ratio** or chemical potential;
it cannot validate k_rev without phase-matched forward kinetics. Conversely,
a propagation benchmark alone cannot determine Kc.

## What the findings would change in ΔH, ΔS and Tc

| Finding / independently defined change | Propagation ΔH and ΔS | Tc | Kinetic effect |
| --- | --- | --- | --- |
| Replace the generic default by a phase- and end-matched k_p, keeping the same Kc and detailed balance | Both unchanged | Unchanged, including the SSA table crossing if both directions use the same grid scaling | Forward and reverse scale together; benchmark sensitivity is quantified above |
| Use a benzylic rather than terminal-primary end | Actual end chemistry can change both; the present I039 additive thermo control reproduces no shift, as shown below | No shift in this implemented control; physical finite-chain effects remain unresolved | Requires the corresponding kinetic channel; the primary-end rate cannot be declared its measured benchmark |
| Add independently measured chain-length/solvent kinetic dependence while keeping Kc fixed | Both unchanged | Unchanged | Can alter local propagation and reverse timescales; a concentration-dependent apparent k_p is not itself a thermo correction |
| Replace gas Kc with independent phase-matched equilibrium evidence | Constrains ΔG = ΔH − TΔS; individual changes remain unidentified by the retrieved points | Direction/size unresolved without common states and temperature dependence; a positive ΔG addition persisting near the root would lower it at fixed monomer activity | Changes reverse/forward ratio; also requires correctly matched forward kinetics |
| Correct monomer activity or polymer activity ratio | Changes the equilibrium condition; a consistent coordinate transformation alone does not change intrinsic physical thermo | Can change the physical crossing when chemical potential changes | Changes the balance of forward and reverse propensities at the specified state |

The structure-only control uses the I039 SMILES
`CC(Ph)C[CH](Ph) + styrene ⇌ CC(Ph)CC(Ph)C[CH](Ph)`, generates resonance
structures, and estimates its thermo without reaction generation or
compilation. It is selected by normal styrene chain-end chemistry, not by
distance to any ceiling comparator.

The tiny signed differences that round to zero below are floating-point
roundoff; the control does not resolve a thermodynamic shift.

<!-- BEGIN I044:structure_control -->
| I039 benzylic end control | Reproduced value |
| --- | --- |
| Continuous Tc at 1 M (K) | 710.020462 |
| ΔH change at 298.15 K (J/mol) | -0.000000 |
| ΔS change at 298.15 K (J/mol/K) | +0.000000 |
| Tc change from compiled primary end (K) | -0.000000 |
<!-- END I044:structure_control -->

Thus improved **forward kinetics alone** provides no evidence that the gas
ceiling will move downward, or that a solvation correction will then point
in the desired direction. Better species/phase thermochemistry could move
it, but that hypothesis requires independent chemical-potential data. The
probe's conclusions would fail if a phase- and end-matched measurement
supported the generic default's magnitude and temperature dependence, or
if the artifact numbers failed the pinned-source reproduction. The latter
check passes; no former evidence was found in this retrieval.

## Contract corrections and remaining scientific limits

- Calling this a compiled *styrene propagation rate* can imply a
  styrene-specific estimate that the pinned source does not supply. The
  template selects styrene chemistry, but its kinetics falls back to a
  generic, unreferenced default; a loaded training set is not a generated
  rate tree. Fixing that limitation is a separate product/database task,
  outside this probe's scope.
- A bulk benzylic-chain benchmark and a primary-end gas oligomer channel are
  different chemical and phase targets. The ratio measures a proxy/model
  discrepancy, not a statistical error bar for the exact molecular reaction.
- Part of the benchmark temperature window lies outside the compiled
  table. Reproducing analytic rule extrapolation there is possible;
  pretending the SSA can evaluate those temperatures would be wrong.
- The independent equilibrium experiment located uses living anionic
  chains in solution. The available quotation does not specify the exact
  solvent assignment. It cannot support a claimed direct numerical
  verification of radical melt k_rev, ΔH, ΔS or Tc. Its full tables and
  activities remain needed for a common-state comparison.
- The later IUPAC reanalysis makes the original confidence assumptions an
  incomplete account of interlaboratory uncertainty. Neither a narrow
  independent A/Ea error bar nor a high-temperature prediction band can be
  claimed from the retrieved material.
- The requested evidence is not fully on the named filesystem: the
  benchmark papers and polymerization equilibrium records required
  literature retrieval. Access failures are recorded above. No inaccessible
  result is presented as reproduced.

There are no user gates or pending choices. The deliverable is a completed
diagnostic record with explicit scientific retrieval limits, rather than a
selection of a correction on the owner's behalf. The remaining scientific
work is obtaining the original equilibrium tables and phase/activity
conditions, quantitative higher-temperature forward data, and matched
radical depropagation evidence. No excluded validation result was used to
choose, tune or rank any change.

## Reproduced Verifier output

```text
I044 artifact, 65 pinned database files and product-source hashes verified
I044 both pairs: default Arrhenius, dimensional Kc and reverse rates independently verified
I044 benchmark ratios, A/Ea factors, equilibrium comparisons and crossings independently verified
I044 public SSA propensity count/volume normalization verified
I044 all 11 numeric report blocks verified
```

The final commit SHA is reported in the worker's closeout, avoiding a
self-referential commit hash in this report.
