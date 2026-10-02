# Polystyrene ceiling-temperature gap: diagnosis only

The prior propagation baselines reproduce. Longer proxies give the same gas
ceiling within numerical precision. The toluene correction raises the ceiling
because its positive reaction entropy correction outweighs its positive
enthalpy correction. The arithmetic gap is reproducible; its decomposition
into independent **physical errors with confidence intervals is not
identifiable** from these group estimates and the retrieved literature.

The commonly quoted comparator is consistent with a compiled **liquid-monomer**
thermodynamic reference. It must not be treated as a measured equilibrium of
these gas-phase radical oligomers at the compiler's monomer concentration.
Converting pressure and concentration standards consistently leaves a physical
equilibrium unchanged. Changing monomer concentration or phase changes the
equilibrium itself. Both distinctions matter here.

## Reproduction and scope

Branch: `i039-tc-probe`. Database: read-only pinned commit
`4a12d36fcdc193ede82c8d1ab5c1653495d445bc`. Only this report and the scripts in
`i039_probe/` are deliverables. No product, database, test or rate-tree change is
made. No excluded local dataset is accessed, and no reaction-yield or pyrolysis
result is used. Literature retrieval is limited to polymerization thermodynamics
and its standard states. Search results outside that scope are not followed or
used. The database is materialized in scratch using the existing committed
`i034_probe/run_probe.py` allowlist and `git show` for every file.

Run from `/home/alon/Code/RMG-Py-kmc-i039-tcprobe`:

```bash
mkdir -p /tmp/i039-tc-probe-build /tmp/i039-reproduce
PYTHONPATH=$PWD /home/alon/anaconda3/envs/rmg_env/bin/python \
  setup.py build_ext --inplace \
  > >(tee -a /tmp/i039-tc-probe-build/stdout.log) \
  2> >(tee -a /tmp/i039-tc-probe-build/stderr.log >&2)
MPLCONFIGDIR=/tmp/i039-reproduce/mpl PYTHONPATH=$PWD \
  /home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/fixtures/i039_probe/run_probe.py \
  --output /tmp/i039-reproduce/results.json \
  > >(tee -a /tmp/i039-reproduce/stdout.log) \
  2> >(tee -a /tmp/i039-reproduce/stderr.log >&2)
PYTHONPATH=$PWD /home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/fixtures/i039_probe/render_tables.py \
  /tmp/i039-reproduce/results.json \
  > >(tee -a /tmp/i039-reproduce/render-stdout.log) \
  2> >(tee -a /tmp/i039-reproduce/render-stderr.log >&2)
MPLCONFIGDIR=/tmp/i039-reproduce/mpl PYTHONPATH=$PWD \
  /home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/fixtures/i039_probe/verify_results.py \
  /tmp/i039-reproduce/results.json \
  > >(tee -a /tmp/i039-reproduce/verify-stdout.log) \
  2> >(tee -a /tmp/i039-reproduce/verify-stderr.log >&2)
```

`render_tables.py --update` creates the numeric blocks during authoring; its
default mode checks them without modifying the report. The verifier re-estimates
species from adjacency lists, checks source and snapshot hashes, independently
calculates equilibrium constants, crossings, LSER corrections and the conditional
allocation, then checks every numeric block. Literature constants are embedded
with citations; reproduction needs no network. This reproduces transcription and
arithmetic, not the underlying historical measurements.

<!-- BEGIN I039:definitions -->
| Quantity | Value |
| --- | --- |
| R (J/mol/K) | 8.31447200 |
| Gas reference pressure (Pa) | 100000 |
| Monomer concentration / concentration reference (mol/m³) | 1000 |
| Conditional thermo anchor (K) | 298.15 |
| Toluene molecular critical temperature (K) | 591.750000 |
| Gas grid minus continuous crossing (K) | +0.228211 |
| Toluene grid minus continuous crossing (K) | +0.166692 |
| Model L=4 minus L=5 continuous crossing (K) | -0.000000 |
| Unapplied optical-term-suppressed control, continuous Tc (K) | 672.121955 |
| Persisted pinned database files | 65 |
<!-- END I039:definitions -->

### Reproduced baselines

The original propagation pair and its rate tables are selected by the existing
I034 baseline routine, which first checks the I032 gas numbers. The solvation
baseline uses precisely the prior linear anchor correction and compiler grid.

<!-- BEGIN I039:baselines -->
| Model at 1 mol/L | Compiler-grid Tc (K) | Continuous Tc (K) | Continuous shift (K) |
| --- | --- | --- | --- |
| gas | 710.248673 | 710.020462 | 0.000000 |
| toluene / linear H−TS | 781.491601 | 781.324909 | +71.304447 |
<!-- END I039:baselines -->

## Equilibrium, participants and conventions

The selected event is `P(L−1)• + styrene ⇌ P(L)•`. The nominal proxy size is
the **product** chain size; the bimolecular declaration uses the shorter chain.
These are H-capped, constitution-only oligomers. The radical in the selected pair
is terminal primary carbon (`[CH2]`), not a terminal benzylic styryl radical.
The exact molecular graphs, selected thermo sources and symmetry quantities are:

<!-- BEGIN I039:species -->
| L | Side | Species (SMILES) | σ including optical | Optical half factors | Thermo source |
| --- | --- | --- | --- | --- | --- |
| 3 | reactant | `[CH2]C(CCc1ccccc1)c1ccccc1` | 4.000000 | 1 | Thermo group additivity estimation: group(Cs-CbCsCsH) + group(Cs-CsCsHH) + group(Cs-CbCsHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene) + ring(Benzene) + radical(Isobutyl) |
| 3 | reactant | `C=Cc1ccccc1` | 2.000000 | 0 | Thermo group additivity estimation: group(Cb-(Cds-Cds)) + group(Cb-H) + group(Cb-H) + group(Cds-CdsCbH) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cds-CdsHH) + ring(Benzene) |
| 3 | product | `[CH2]C(CC(CCc1ccccc1)c1ccccc1)c1ccccc1` | 4.000000 | 2 | Thermo group additivity estimation: group(Cs-CbCsCsH) + group(Cs-CbCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CbCsHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene) + ring(Benzene) + ring(Benzene) + radical(Isobutyl) |
| 4 | reactant | `[CH2]C(CC(CCc1ccccc1)c1ccccc1)c1ccccc1` | 4.000000 | 2 | Thermo group additivity estimation: group(Cs-CbCsCsH) + group(Cs-CbCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CbCsHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene) + ring(Benzene) + ring(Benzene) + radical(Isobutyl) |
| 4 | reactant | `C=Cc1ccccc1` | 2.000000 | 0 | Thermo group additivity estimation: group(Cb-(Cds-Cds)) + group(Cb-H) + group(Cb-H) + group(Cds-CdsCbH) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cds-CdsHH) + ring(Benzene) |
| 4 | product | `[CH2]C(CC(CC(CCc1ccccc1)c1ccccc1)c1ccccc1)c1ccccc1` | 4.000000 | 3 | Thermo group additivity estimation: group(Cs-CbCsCsH) + group(Cs-CbCsCsH) + group(Cs-CbCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CbCsHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene) + ring(Benzene) + ring(Benzene) + ring(Benzene) + radical(Isobutyl) |
| 5 | reactant | `[CH2]C(CC(CC(CCc1ccccc1)c1ccccc1)c1ccccc1)c1ccccc1` | 4.000000 | 3 | Thermo group additivity estimation: group(Cs-CbCsCsH) + group(Cs-CbCsCsH) + group(Cs-CbCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CbCsHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene) + ring(Benzene) + ring(Benzene) + ring(Benzene) + radical(Isobutyl) |
| 5 | reactant | `C=Cc1ccccc1` | 2.000000 | 0 | Thermo group additivity estimation: group(Cb-(Cds-Cds)) + group(Cb-H) + group(Cb-H) + group(Cds-CdsCbH) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cds-CdsHH) + ring(Benzene) |
| 5 | product | `[CH2]C(CC(CC(CC(CCc1ccccc1)c1ccccc1)c1ccccc1)c1ccccc1)c1ccccc1` | 4.000000 | 4 | Thermo group additivity estimation: group(Cs-CbCsCsH) + group(Cs-CbCsCsH) + group(Cs-CbCsCsH) + group(Cs-CbCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CbCsHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene) + ring(Benzene) + ring(Benzene) + ring(Benzene) + ring(Benzene) + radical(Isobutyl) |
<!-- END I039:species -->

For this association, Δν = −1. RMG calculates

```text
Kp° = exp(−ΔGp° / RT)
Kc [m³/mol] = Kp° RT / p°
ceiling: Kc [M] = 1
α(T) = c°RT / p°
ΔGc° = ΔGp° − RT ln α
ΔHc° = ΔHp° + RT
ΔSc° = ΔSp° + R(ln α + 1)
```

The last two equations are the derivatives at fixed concentration standard.
Writing an effective balance entropy `ΔSp° + R ln α` while keeping `ΔHp°`
unchanged gives the same free energy, but that effective entropy is **not** the
thermodynamic `ΔSc°`. Omitting this distinction misassigns enthalpy and entropy.
The forward/reverse-rate ratio is checked against Kc on every compiler grid point.
Kc's SI dimensions do not make its gas thermochemistry a bulk-liquid reference.

<!-- BEGIN I039:sizes -->
| L | Chain reaction | Grid Tc (K) | Continuous Tc (K) | Same-state 1 M Tc (K) | 1-bar monomer Tc (K) | ΔH at anchor, gas (kJ/mol) | ΔS at anchor, gas (J/mol/K) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 3 | P2• + styrene ⇌ P3• | 710.248673 | 710.020462 | 710.020462 | 543.888008 | -80.039992 | -147.413327 |
| 4 | P3• + styrene ⇌ P4• | 710.248673 | 710.020462 | 710.020462 | 543.888008 | -80.039992 | -147.413327 |
| 5 | P4• + styrene ⇌ P5• | 710.248673 | 710.020462 | 710.020462 | 543.888008 | -80.039992 | -147.413327 |
<!-- END I039:sizes -->

The same-state pressure-to-concentration conversion reproduces the ceiling with
zero coordinate shift. The separate column with monomer at the gas reference
pressure changes the **physical monomer condition**. It is not a correction to
the compiler. A pure-liquid reference additionally requires a monomer transfer
chemical potential or activity and polymer-phase chemical-potential increments.
Neither a unit conversion nor an arbitrary density substitution supplies those.

### What the literature comparator means

The liquid and gas columns of the *Kirk-Othmer* compilation reproduce the quoted
liquid-monomer ceiling by the low-temperature H/S ratio. Their gas column means
gas **monomer** polymerizing to condensed polymer; it does not mean gas-phase
polymer chains. The retrieved table specifies the thermodynamic anchor temperature
and pressure. It supplies no high-temperature heat-capacity integration or
uncertainty. Thus the comparator is a compilation/extrapolation, not evidence of
a directly measured, pressure-independent melt equilibrium.
[Kirk-Othmer, “Xylylene Polymers,” Table 2](https://catalogimages.wiley.com/images/db/pdf/9781118766927.excerpt.pdf).

Odian's table explicitly uses a unit-concentration monomer entropy and gives a
different H/S pair. Its footnote does not specify a reference temperature.
Cowie's table explicitly uses pure liquid monomer and reports a different ceiling.
These are distinct reference cases, not interchangeable error bars.
[Odian, Table 3-15](https://www.eng.uc.edu/~beaucag/Classes/Properties/Books/George%20Odian%20-%20Principles%20of%20Polymerization-Wiley-Interscience%20%282004%29.pdf),
[Cowie & Arrighi, Table 3.6](https://www.eng.uc.edu/~beaucag/Classes/Properties/Books/Arrighi%2C%20Valeria_%20Cowie%2C%20J.M.G%20-%20Polymers%20_%20Chemistry%20and%20Physics%20of%20Modern%20Materials%2C%20Third%20Edition-CRC%20Press%20%282007%29.pdf).

The calorimetric and entropy evidence is kept phase-specific. Roberts et al.
measure liquid styrene to solid PS and separately account for dissolution in
styrene. Warfield & Petree estimate a condensed-state entropy using measured heat
capacities plus a low-temperature model; their result is not a gas-oligomer
entropy or a measured high-temperature melt entropy.
[Roberts, Walton & Jessup](https://nvlpubs.nist.gov/nistpubs/jres/38/jresv38n6p627_A1b.pdf),
[Warfield & Petree](https://onlinelibrary.wiley.com/doi/abs/10.1002/pol.1961.1205516208).

<!-- BEGIN I039:literature -->
| Reference | ΔHp (kJ/mol) | ΔSp (J/mol/K) | Temperature | State / uncertainty |
| --- | --- | --- | --- | --- |
| Roberts et al.: solid PS | -69.79 | not measured here | 298.15 K | pure liquid styrene → solid PS; published ±0.66 kJ/mol |
| Roberts et al.: solution | -73.38 | not measured here | 298.15 K | PS in styrene, 6.9 wt%; published ±0.69 kJ/mol |
| Warfield & Petree | not a formation enthalpy datum | -111.67096 | 298.16 K | condensed-state third-law estimate; entropy loss 26.69 cal/mol/K; uncertainty not supplied |
| Odian Table 3-15 | -73 | -104 | not specified in table footnote | liquid-monomer ΔH; 1 M monomer ΔS; uncertainties not supplied |
| Odian constant-H/S calculation | — | — | 701.923077 K | 1 M ceiling; calculation, not a measured ceiling |
| Cowie Table 3.6 | -68.5 | not specified | 583 K | pure liquid monomer ceiling; uncertainty not supplied |
| Kirk-Othmer liquid-monomer column | -69.9 | -104.6 | 298.15 K | liquid monomer / condensed polymer reference, 101.3 kPa; tabulated Tc=395 °C; no error bars |
| Kirk-Othmer gas-monomer column | -113.4 | -212.1 | 298.15 K | gas monomer / condensed polymer reference, 101.3 kPa; tabulated Tc=262 °C; no error bars |
| Dispatch's literature comparator | not supplied | not supplied | 668 K | consistent with rounded Kirk-Othmer liquid-monomer value; source identity not supplied in dispatch |
<!-- END I039:literature -->

**Common monomer convention:** the compilation's gas-monomer H/S pair is
normalized to the compiler's pressure reference and monomer concentration below.
The pressure change uses `Sbar = Ssource + R ln(pbar/psource)`; its concentration
case solves `H − TSbar − RT ln(c°RT/pbar) = 0`. This shares the monomer convention
with the model. The polymer-phase difference and frozen literature heat capacities
remain. A complete pure-liquid-to-melt common-state comparison is not established.

<!-- BEGIN I039:literature_conventions -->
| Recomputation / comparison | Temperature or difference (K) |
| --- | --- |
| Compiled literature liquid H/S ratio (frozen anchor) | 668.260038 |
| Compiled literature gas H/S ratio at source pressure | 534.653465 |
| Compiled literature gas reference at 1 bar | 534.382894 |
| Compiled literature gas reference, monomer at 1 M | 632.601028 |
| Monomer-condition shift in literature gas case: 1 bar → 1 M | +98.218134 |
| Model minus compiled literature on gas-monomer 1 M convention | +77.419434 |
| Model monomer-condition shift: 1 bar → 1 M | +166.132455 |
<!-- END I039:literature_conventions -->

The formal same-state convention contribution is exactly zero. The displayed
monomer-condition shifts show why changing a reference case can nevertheless
produce a substantial temperature difference. The fraction of the original
liquid-reference gap due **only** to concentration/phase cannot be measured from
these inputs: it needs an independently specified activity/transfer model.

## Gas enthalpy and entropy decomposition

All selected baseline species use group additivity: there are no hits in the
loaded `primaryThermoLibrary`. The scripts parse RMG's source comments, resolve
group-data aliases in the pinned snapshot and reconstruct each species at every
reported temperature. Full per-species source weights and values are in the
reproduced JSON. The net reaction contributions are:

<!-- BEGIN I039:groups -->
| Net source | Net weight | ΔH at anchor (kJ/mol) | ΔS at anchor (J/mol/K) | ΔH at 700 K (kJ/mol) | ΔS at 700 K (J/mol/K) |
| --- | --- | --- | --- | --- | --- |
| HBI:subtract_H_atoms | 0 | 0.000000 | 0.000000 | 0.000000 | 0.000000 |
| group:Cb-(Cds-Cds) | -1 | -23.809208 | 32.627657 | -31.102485 | 17.505967 |
| group:Cb-Cs | 1 | 23.055510 | -32.169357 | 29.177461 | -19.658043 |
| group:Cb-H | 0 | 0.000000 | -0.000000 | 0.000000 | -0.000000 |
| group:Cds-CdsCbH | -1 | -28.370303 | -26.703257 | -39.338568 | -49.114701 |
| group:Cds-CdsHH | -1 | -26.195026 | -115.530927 | -38.603361 | -140.801963 |
| group:Cs-CbCsCsH | 1 | -4.097279 | -50.825397 | 8.647175 | -24.892196 |
| group:Cs-CbCsHH | 0 | 0.000000 | 0.000000 | 0.000000 | 0.000000 |
| group:Cs-CsCsHH | 1 | -20.623686 | 39.424802 | -7.022965 | 67.091710 |
| group:Cs-CsHHH | 0 | 0.000000 | 0.000000 | 0.000000 | 0.000000 |
| optical:atom_half_factors | applied per species | 0.000000 | 5.763153 | 0.000000 | 5.763153 |
| radical:Isobutyl | 0 | 0.000000 | 0.000000 | 0.000000 | 0.000000 |
| ring:Benzene | 0 | 0.000000 | 0.000000 | 0.000000 | 0.000000 |
| symmetry:without_optical_factor | applied per species | 0.000000 | 0.000000 | 0.000000 | 0.000000 |
<!-- END I039:groups -->

The same primary-radical HBI group occurs on both chains and cancels, as does
the subtraction of atomic-H enthalpy. There is no identified radical-energy
offset in this selected pair. Benzene ring corrections and unchanged capping
groups also cancel. The remaining bond-environment groups determine the
reaction enthalpy, heat capacity and most of the entropy. Group uncertainties
are not supplied as a usable correlated reaction-error model; absent quantity
uncertainty is not evidence of exact thermochemistry.

RMG applies chirality through atom symmetry half factors. The script counts
those factors on the same resonance hybrid used for species symmetry, factors
them out of σ, and reports `−R ln σ_nonoptical` and the optical entropy separately.
This is RMG's constitutional optical counting, not a tacticity-resolved polymer
partition function. Libraries would retain their own entropy corrections; none
is present in these baseline participants. No molecular hindered-rotor
calculation is run by this group estimate.

The definitions table also reports a control that suppresses just the net optical
entropy term. This measures that term's temperature sensitivity. It does not
justify removing configurational entropy: the appropriate counting depends on
the specified polymer stereochemical ensemble. This control is not applied.

<!-- BEGIN I039:thermo -->
| L | T (K) | ΔH gas, 1 bar (kJ/mol) | ΔS gas, 1 bar (J/mol/K) | ΔH gas, 1 M (kJ/mol) | ΔS gas, 1 M (J/mol/K) | ΔCp gas, 1 bar (J/mol/K) |
| --- | --- | --- | --- | --- | --- | --- |
| 3 | 298.15 | -80.039992 | -147.413327 | -77.561032 | -112.405874 | -0.474905 |
| 3 | 600.00 | -78.963438 | -145.215163 | -73.974755 | -104.393125 | 6.736240 |
| 3 | 668.00 | -78.483608 | -144.458222 | -72.929541 | -102.743555 | 7.376392 |
| 3 | 700.00 | -78.242744 | -144.106072 | -72.422613 | -102.002352 | 7.677640 |
| 3 | 800.00 | -77.427910 | -143.019411 | -70.776332 | -99.805448 | 8.619040 |
| 4 | 298.15 | -80.039992 | -147.413327 | -77.561032 | -112.405874 | -0.474905 |
| 4 | 600.00 | -78.963438 | -145.215163 | -73.974755 | -104.393125 | 6.736240 |
| 4 | 668.00 | -78.483608 | -144.458222 | -72.929541 | -102.743555 | 7.376392 |
| 4 | 700.00 | -78.242744 | -144.106072 | -72.422613 | -102.002352 | 7.677640 |
| 4 | 800.00 | -77.427910 | -143.019411 | -70.776332 | -99.805448 | 8.619040 |
| 5 | 298.15 | -80.039992 | -147.413327 | -77.561032 | -112.405874 | -0.474905 |
| 5 | 600.00 | -78.963438 | -145.215163 | -73.974755 | -104.393125 | 6.736240 |
| 5 | 668.00 | -78.483608 | -144.458222 | -72.929541 | -102.743555 | 7.376392 |
| 5 | 700.00 | -78.242744 | -144.106072 | -72.422613 | -102.002352 | 7.677640 |
| 5 | 800.00 | -77.427910 | -143.019411 | -70.776332 | -99.805448 | 8.619040 |
<!-- END I039:thermo -->

Directly calling the mismatch with condensed-state measurements an enthalpy or
entropy **error** would confound phase transfer, proxy chemistry, polymer activity
and heat-capacity extrapolation. The table permits the comparison without hiding
those missing terms.

### Structural control

The terminal-benzylic control is chosen from the usual styrene connectivity:
`CC(Ph)C[CH](Ph) + styrene ⇌ CC(Ph)CC(Ph)C[CH](Ph)`. Resonance structures are
generated before estimation, following the family-generation baseline's aromatic
representation. This is a diagnostic comparison only; it is not compiled into
the event set and is not chosen by its proximity to any ceiling.

<!-- BEGIN I039:chain_end -->
| Structural control | Continuous Tc (K) | ΔH at anchor, 1 bar (kJ/mol) | ΔS at anchor, 1 bar (J/mol/K) | Shift from compiled primary-end pair (K) |
| --- | --- | --- | --- | --- |
| terminal benzylic radical | 710.020462 | -80.039992 | -147.413327 | -0.000000 |
<!-- END I039:chain_end -->

The control's exact graphs, chosen source groups and species decomposition are
persisted in JSON. Its result tests whether moving the radical end changes this
particular thermodynamic increment; it does not validate the original end's
propagation kinetics.

## Solvation enthalpy and entropy decomposition

The prior toluene model adds `ΔΔG(T) = ΔΔHanchor − T ΔΔSanchor` to the reaction.
It assumes zero solvation heat-capacity correction. The anchor LSERs are
equal-concentration gas-to-liquid transfer descriptors. Applying them after the
gas concentration-standard conversion is consistent; calling the result bulk
PS thermo is an additional physical approximation.

The pinned numerical implementation uses its own approximate constants in the
Abraham free-energy expression and separate Mintz enthalpy coefficients. The
probe reproduces those constants exactly, including the solvent intercept for
the molecule-count change. First, decomposition by LSER descriptor:

<!-- BEGIN I039:solvation -->
| Toluene LSER term | Reaction descriptor / intercept count | ΔΔH (kJ/mol) | ΔΔS (J/mol/K) | ΔΔG at anchor 298 K (kJ/mol) |
| --- | --- | --- | --- | --- |
| A | 0.00000000 | -0.000000 | 0.000000 | -0.000000 |
| B | -0.00670000 | 0.023549 | 0.063852 | 0.004521 |
| E | -0.21318000 | -1.323453 | -3.144947 | -0.386259 |
| L | -0.13081000 | 1.094286 | 1.086888 | 0.770393 |
| S | -0.06913000 | 0.942677 | 1.966126 | 0.356772 |
| intercept | -1.00000000 | 5.656010 | 17.426109 | 0.463030 |
| Sum | — | 6.393069 | 17.398028 | 1.208457 |
<!-- END I039:solvation -->

Second, decomposition by the source groups that generated those descriptors:

<!-- BEGIN I039:solvation_groups -->
| Net solute descriptor source | Weight | ΔΔH (kJ/mol) | ΔΔS (J/mol/K) |
| --- | --- | --- | --- |
| group:Cb-(Cds-Cd) | -1 | 6.794926 | 7.221209 |
| group:Cb-Cs | 1 | -7.829633 | -8.646055 |
| group:Cb-H | 0 | -0.000000 | -0.000000 |
| group:Cds-CdsCbH | -1 | 2.461608 | 1.025232 |
| group:Cds-CdsHH | -1 | 3.528467 | 4.270046 |
| group:Cs-CbCsCsH | 1 | -0.097650 | 0.194288 |
| group:Cs-CbCsHH | 0 | 0.000000 | 0.000000 |
| group:Cs-CsCsHH | 1 | -4.120658 | -4.092802 |
| group:Cs-CsHHH | 0 | 0.000000 | 0.000000 |
| ring:Benzene | 0 | 0.000000 | 0.000000 |
| LSER intercept (molecule-count change) | -1 | 5.656010 | 17.426109 |
<!-- END I039:solvation_groups -->

The individual descriptor sources and library-hit status are:

<!-- BEGIN I039:solvation_sources -->
| Species | S | B | E | L | A | Source |
| --- | --- | --- | --- | --- | --- | --- |
| `C=Cc1ccccc1` | 0.616880 | 0.185010 | 0.835790 | 4.102870 | 0.000000 | group(Cb-(Cds-Cd)) + group(Cb-H) + group(Cb-H) + group(Cds-CdsCbH) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cds-CdsHH) + ring(Benzene) |
| `[CH2]C(CC(CCc1ccccc1)c1ccccc1)c1ccccc1` | 1.637730 | 0.524030 | 1.866740 | 11.892790 | 0.000000 | group(Cs-CbCsCsH) + group(Cs-CbCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CbCsHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene) + ring(Benzene) + ring(Benzene) |
| `[CH2]C(CCc1ccccc1)c1ccccc1` | 1.089980 | 0.345720 | 1.244130 | 7.920730 | 0.000000 | group(Cs-CbCsCsH) + group(Cs-CsCsHH) + group(Cs-CbCsHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-Cs) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene) + ring(Benzene) |
<!-- END I039:solvation_sources -->

The carbon-radical lookup fails at the generic radical node and the public
estimator retains the saturated-analogue descriptors. This independently
reproduces the I034 limitation; it is not a measured zero radical-solvation
correction. The two omitted radical descriptors cancel only under that implemented
assumption. They need not cancel physically in a polymer environment.

<!-- BEGIN I039:radicals -->
| Radical | Matched node | Lookup result | Descriptors minus saturated |
| --- | --- | --- | --- |
| `[CH2]C(CC(CCc1ccccc1)c1ccccc1)c1ccccc1` | R_rad | 'Node has no parent with data in database.' | A=+0.000000, B=+0.000000, E=+0.000000, L=+0.000000, S=+0.000000, V=-0.021500 |
| `[CH2]C(CCc1ccccc1)c1ccccc1` | R_rad | 'Node has no parent with data in database.' | A=+0.000000, B=+0.000000, E=+0.000000, L=+0.000000, S=+0.000000, V=-0.021500 |
<!-- END I039:radicals -->

The net entropy shift comes almost entirely from the LSER intercept associated
with the molecule-count change; the source-group terms largely cancel. This
identifies which part of the surrogate drives the shift, without establishing
that it is an erroneous intercept or a PS melt property.

The positive solvation enthalpy makes association less exothermic, while the
positive solvation entropy makes its entropy loss smaller. At the reported
temperature window the net correction stabilizes propagation, moving the root
upward. The effective corrected concentration-standard thermodynamics are:

<!-- BEGIN I039:corrected_thermo -->
| T (K) | Toluene-corrected ΔH, 1 M (kJ/mol) | Toluene-corrected ΔS, 1 M (J/mol/K) | ΔΔGsolv (kJ/mol) |
| --- | --- | --- | --- |
| 298.15 | -71.167963 | -95.007846 | 1.205847 |
| 600.00 | -67.581686 | -86.995097 | -4.045748 |
| 668.00 | -66.536472 | -85.345527 | -5.228813 |
| 700.00 | -66.029545 | -84.604325 | -5.785550 |
| 800.00 | -64.383263 | -82.407421 | -7.525353 |
<!-- END I039:corrected_thermo -->

Unapplied component controls isolate the enthalpy, entropy, intercept and
descriptor effects. They use the same gas reaction and concentration and are
not ranked by closeness to a literature ceiling.

<!-- BEGIN I039:solvation_controls -->
| Unapplied solvation component control | Continuous Tc (K) | Shift from gas (K) |
| --- | --- | --- |
| both | 781.324909 | +71.304447 |
| descriptor terms only | 702.590676 | -7.429786 |
| enthalpy only | 647.645182 | -62.375280 |
| entropy only | 859.300652 | +149.280190 |
| intercept only | 790.504984 | +80.484522 |
<!-- END I039:solvation_controls -->

The root lies above the molecular toluene critical temperature. The linear
surrogate returns a number there because it has no molecular-fluid temperature
guard. That is not evidence that real liquid toluene or a PS melt has that
thermodynamic correction. The high-temperature K-factor model cannot supply
this root in its valid molecular-liquid domain. Source evidence is the pinned
`solvation.py`, `thermoengine.py`, and the earlier I034 baseline; no fitted
continuation is introduced.

## Proxy-size effect

The longer-proxy table above is compiled independently for each size using the
same family selection and monomer concentration. The added remote repeat units
cancel exactly in these additive reaction increments, including the net optical
factor. No size trend or failure is observed in this range.

This is a convergence statement about the **implemented estimator**, not a
bound on actual finite-chain entropy. Whole-chain translational entropy,
conformer populations, correlated torsions, tacticity and polymer activities are
not established by increasing the length of an additive proxy. The observed
zero trend cannot be used as a physical zero-uncertainty result.

## Additive gap accounting and its uncertainties

**Conditional arithmetic, not uniquely identified physical errors.** To supply
an auditable table whose rows sum to the displayed compiled gap, use the longest
gas proxy, the compiler's monomer concentration, and Odian's explicitly
unit-concentration H/S pair. Assume that pair is anchored at the script's stated
reference temperature and carries the model's concentration-standard heat-capacity
curve. The table footnote does not establish the anchor; this is an explicit
diagnostic assumption, not a measurement or parameter applied to product code.

Construct four roots by replacing neither, either, or both of the anchor H and S
with the model values while holding that heat-capacity curve fixed. Average the
two orders of H/S replacement to distribute their nonlinear interaction equally.
This allocation rule is symmetric and reproducible; thermodynamics does not
make it a unique causal explanation.

<!-- BEGIN I039:counterfactuals -->
| Anchor H | Anchor S | Continuous Tc (K) |
| --- | --- | --- |
| Odian (conditional anchor) | Odian (conditional anchor) | 725.120661 |
| Odian (conditional anchor) | model | 665.429449 |
| model | Odian (conditional anchor) | 774.437841 |
| model | model | 710.020462 |
<!-- END I039:counterfactuals -->

<!-- BEGIN I039:attribution -->
| Contribution to compiled toluene Tc minus dispatch comparator | ΔT (K) | Rounding-only sensitivity (K) | Scientific uncertainty |
| --- | --- | --- | --- |
| Convention (same physical state) | +0.000000 | ±0.000000 | exact coordinate identity |
| Enthalpy error (conditional allocation) | +46.954097 | ±5.291279 | not identifiable; no confidence interval |
| Entropy error (conditional allocation) | -62.054295 | ±4.351408 | not identifiable; no confidence interval |
| Proxy size (L=5 → L=3, observed model shift) | -0.000000 | ±0.000000 | not identifiable; no confidence interval |
| Solvation (plus grid interpolation) | +71.471139 | ±0.000000 | not identifiable; no confidence interval |
| Residual: physics not represented / unresolved reference | +57.120661 | ±9.332613 | not identifiable; no confidence interval |
| Sum | +113.491601 | correlated contributions; do not add error bars | not identifiable |
<!-- END I039:attribution -->

The convention row is the exact same-state coordinate identity. The proxy row
is the observed continuous-root difference. The solvation row includes the
displayed grid interpolation offset so that the sum targets the **compiled**
corrected crossing. The residual is labelled “physics not represented” as
requested, but also contains the unresolved literature reference and assumed
heat-capacity continuation: **its magnitude is not a demonstrated physical
missing-energy term**. Assigning it solely to melt physics would be unsupported.

Sensitivity columns perturb the printed Odian anchors over half their last
printed unit. They are rounding-only deviations, not statistical error bars.
Contributions are correlated and their rounding effects cancel in the sum.
No physical confidence interval can be given for the H/S allocations, proxy
bias, solvation transferability or residual from the available evidence. In
particular, an exactly computed coordinate identity and an observed size
difference do not provide a bound on the missing phase/activity physics.

## A-priori options, not applied or ranked

| Identified issue | Independently justifiable option | What it changes / evidence required |
| --- | --- | --- |
| Monomer condition and standard state | Specify the actual monomer activity/concentration and pressure; transform gas and condensed references consistently. | Changes the physical ceiling when monomer chemical potential changes. A coordinate conversion alone leaves the root unchanged. Needs independently specified reactor/reference thermodynamics. |
| Condensed versus gas thermochemistry | Use phase-matched polymerization calorimetry and entropy/heat-capacity measurements as constraints on a reference thermo model. | Changes H, S and their temperature dependence. Roberts' solid and solution values must remain separate; neither is a direct high-temperature radical-oligomer group value. No group is replaced by an amount chosen to reproduce a ceiling. |
| Local groups and stereochemistry | Benchmark the surviving alkene/alkyl/phenyl groups using independent molecular thermochemistry or QM; represent a specified tacticity and conformer ensemble. | Changes the local H/S increment and optical/conformer counting. A pure calorimetric polymerization enthalpy cannot uniquely repair one group without disentangling phase and other groups. |
| Proxy representation and ends | Use chemically justified radical ends and direct oligomer thermochemistry; examine explicit rotor/conformer increments over increasing chain lengths. | Can change kinetics and finite-chain thermo. The present additive estimator shows no length trend, so extending it alone supplies no demonstrated thermo correction. |
| Radical solvation | Obtain carbon-radical transfer free energies by an independently validated QM/explicit-environment protocol or non-pyrolysis equilibrium measurements. | Changes the HBI descriptors currently absent, potentially their reaction difference. Requires a genuine radical reference, not the public estimator's silent saturated fallback. |
| Melt surrogate and temperature dependence | Calculate polymer-specific chemical potentials/activities, or validate a segment-reference solvation model against independent solution/melt equilibrium thermo. | Changes transfer H/S, heat-capacity curvature and concentration/activity dependence. Molecular toluene critical behavior and an unvalidated linear extension do not determine a PS melt model. |

No option is fitted to the quoted ceiling or to excluded validation results.
No correction is applied and no option is selected on the owner's behalf.

## Contract corrections, retrieval limits and remaining scientific work

- The specified literature temperature is not a universal PS constant. A
  liquid-monomer compilation supports its approximate magnitude, but the gas
  oligomer/concentration comparison differs in phase and physical monomer
  condition. A complete common **liquid/melt** comparison still needs chemical
  potentials or activities.
- The prior “compiled toluene” number is an I034 diagnostic correction to gas
  rate ratios. `compiler.py` itself uses gas-phase equilibrium constants; its
  non-ring record status still records missing condensed-phase reference
  thermochemistry. The number does not establish an implemented, validated melt
  reverse-rate model.
- The compiled baseline selects a terminal-primary radical pair. The normal
  benzylic structural control is explicitly checked; source-group cancellation
  and the measured control result, rather than a structural suspicion alone,
  determine what can be attributed thermodynamically.
- A unique six-way causal allocation with numerical physical uncertainties is
  underdetermined. H/S allocations depend on reference phase, anchor temperature,
  heat-capacity continuation and interaction-allocation rule. The conditional
  table above is the reproducible arithmetic that can be supplied; reporting it
  as measured error contributions would be wrong.

**Not retrieved:** the original measurement/provenance chain and full standard-state
footnotes behind every compiled ceiling; matched high-temperature melt propagation
H/S/Cp, polymer activities, or radical transfer uncertainties. The historical
Dainton–Ivin publisher records were located but their full articles were not
accessible through the retrieval tools. Odian, Cowie and Kirk-Othmer values and
footnotes were available as indexed excerpts; PDF retrieval was unavailable.
Warfield's publisher abstract and Roberts' original NBS abstract supplied the
stated empirical values. This is not an exhaustive literature review or a claim
that missing measurements do not exist.

The decisive further measurement is phase-matched propagation chemical potentials
with a specified monomer activity and temperature dependence. Such evidence could
change the sign or size of the inferred thermo mismatch and invalidate the
conditional allocation without changing any reproduced model number here.

## Reproduced verifier output

```text
I039 pinned snapshot and product-source hashes verified
I039 propagation baselines, L=3,4,5 thermo and equilibrium roots independently reproduced
I039 toluene LSER and conditional additive gap accounting independently reproduced
I039 all 17 numeric report blocks verified
```

The final commit SHA is reported in the worker's closeout rather than embedded
in this report, which avoids a self-referential commit hash.
