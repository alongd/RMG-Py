# A-priori PS-melt stand-in for solvation-corrected reverse rates

## Scope and recommendation

**Recommendation to the owner, not an adopted choice:** use the pinned **toluene**
Abraham/Mintz coefficients with the **constant-solvation-enthalpy/entropy,
linear H−TS model**. Interpret this as a noncritical, aromatic segment-reference
approximation, **not literal liquid toluene at the melt temperature**. Freeze the
database revision, descriptors, reference-state convention and model before any
comparison with validation data. Use **benzene with the same linear model** as
the sensitivity arm. Neither arm is a validated PS-melt thermochemical model.

The rationale is chemical and implementation-based: toluene has an alkylated
aromatic ring, no hydrogen-bond donor or strong acceptor, and the complete
Gibbs/enthalpy coefficient sets required by the existing RMG linear model.
Ethylbenzene is closer to a capped PS repeat unit, but its pinned solvent entry
lacks every Mintz enthalpy coefficient. Its apparent chemical advantage cannot
justify inventing those coefficients or borrowing another solvent's enthalpy.
The near-critical K-factor model cannot span the requested melt interval for
any of the aromatic candidates examined. These are database/model facts, not a
ranking against a reaction, mass-loss or product dataset.

No product code, database, compiler, executor or test is changed. This report
and its two scripts are the entire deliverable. The worktree starts at
`ff8e75ce48b0ff4e5ce8a3a28b985932f784149c`; the database is read only at
`4a12d36fcdc193ede82c8d1ab5c1653495d445bc`. Nothing under `polymers/` or a
`catalog/` directory is used. No pyrolysis, TGA, MWD or product-yield dataset is
consulted. No rate tree is regenerated. No molecular-solvent critical point is
assigned to the PS melt.

## Reproduction and Verifier

Run **both committed scripts** from `/home/alon/Code/RMG-Py-kmc-i034-meltprobe`:

```bash
mkdir -p /tmp/i034-reproduce
MPLCONFIGDIR=/tmp/i034-reproduce/mpl PYTHONPATH=$PWD \
  /home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/fixtures/i034_probe/run_probe.py \
  --database /home/alon/Code/RMG-database \
  --scratch /tmp/i034-reproduce/database \
  --output /tmp/i034-reproduce/results.json \
  > >(tee -a /tmp/i034-reproduce/stdout.log) \
  2> >(tee -a /tmp/i034-reproduce/stderr.log >&2)
MPLCONFIGDIR=/tmp/i034-reproduce/mpl PYTHONPATH=$PWD \
  /home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/fixtures/i034_probe/verify_results.py \
  /tmp/i034-reproduce/results.json \
  --snapshot /tmp/i034-reproduce/database \
  > >(tee -a /tmp/i034-reproduce/verify-stdout.log) \
  2> >(tee -a /tmp/i034-reproduce/verify-stderr.log >&2)
```

The probe materializes only the named solvation/thermo/family files from the
pinned Git objects with `git show`, into scratch; it neither trusts a potentially
dirty database checkout nor writes to it. Its initial assertions reproduce the
prior gas numbers **before** performing solvation calculations. The prior script
is read from the fixed `i032-probe` commit
`413887456db0502e6a452da74d0ed4d85127bb53`; its proxy construction and selection
criteria are reused without opening the other campaign file it references.

The verifier reloads the scratch thermo and solvation database, re-estimates
the selected species, independently evaluates the LSER equations and the
piecewise K-factor model from CoolProp densities, recomputes the ceiling roots,
and compares every generated numerical report block exactly. Its default is
read-only. `--update-report-blocks` regenerates the marked tables after the same
checks; it is a report-maintenance option, not required for reproduction.
Unrounded species corrections, error messages, descriptor provenance and
ceiling residuals are in the reproduced JSON, rather than committed scratch
results. Literature measurements, bibliographic identifiers and code references
are explicitly attributed; they are not claimed as measurements by this probe.

Actual verifier output:

```text
I-034 gas baseline, descriptors, independent LSER/K-factor calculations and ceilings verified
I-034 52 corrected reaction-temperature rows and 15 report blocks reproduced
```

<!-- BEGIN I034:provenance -->
| Quantity | Reproduced value |
| --- | --- |
| Python | 3.9.23 |
| CoolProp | 6.7.0 |
| R (J/mol/K) | 8.314472 |
| Pinned snapshot files | 65 |
| Snapshot SHA256 | `61e2c7357bb4945dce6e6a1a39eccf216e945a01771e461391e7975165ea29eb` |
| Prior probe script SHA256 | `0ca66f203db5f24810d622f7ede5e8c92368134b77514f22c8cf4a55fd144e03` |
| Monomer concentration (mol/m³) | 1000.0 |
| Default ceiling grid (K) | 300–1000, step 25 |
<!-- END I034:provenance -->

## 1. Pinned database and implemented models

### Solvent inventory

For these neutral organic probes, the Abraham gas-to-solvent Gibbs LSER needs
the six solvent `s_g,b_g,e_g,l_g,a_g,c_g` coefficients and the solute
`S,B,E,L,A` descriptors. The Mintz enthalpy LSER additionally needs the six
`s_h,b_h,e_h,l_h,a_h,c_h` coefficients. `V` is available as a McGowan volume
but does not enter these gas-to-solvent LSER equations. A CoolProp solvent
identity **is not required by the linear H−TS model**; it is required by the
temperature-dependent K-factor model. Viscosity parameters, dielectric constant
and solvent `alpha,beta` are not inputs to the equilibrium correction below.

The following inventory is exhaustive for the pinned solvent library. Coefficient
availability is not a claim of applicability to PS or to high temperatures.
The legacy `c` label is reported as found, without repairing the database.

<!-- BEGIN I034:inventory -->
Pinned solvent entries: **203**. Solvent-library lines mentioning polymer/polystyrene/melt: **0**.

**Complete Abraham + Mintz (69):** `1,1-dichloroethane`, `1,2-dichlorobenzene`, `1,2-dichloroethane`, `1,4-dioxane`, `1-methylpyrrolidin-2-one`, `1-octanol`, `2-(2-hydroxyethoxy)ethanol`, `2-chlorobutane`, `2-methoxyethanol`, `2-methylbutan-2-ol`, `2-methylpropan-1-ol`, `2-methylpropan-2-ol`, `2-propan-2-yloxypropane`, `4-methyl-1,3-dioxolan-2-one`, `N,N-dibutylbutan-1-amine`, `N,N-dioctyloctan-1-amine`, `N,N-dipentylpentan-1-amine`, `acetic acid`, `acetonitrile`, `benzene`, `bis(2-ethylhexyl) hexanedioate`, `bis(2-ethylhexyl) phthalate`, `butan-2-ol`, `butanol`, `c`, `carbontet`, `chlorobenzene`, `chloroform`, `cyclohexane`, `cyclopentanone`, `decane`, `dibutylether`, `dichloromethane`, `diethyl benzene-1,2-dicarboxylate`, `dimethylformamide`, `dimethylsulfoxide`, `dodecan-1-ol`, `dodecane`, `ethane-1,2-diol`, `ethanol`, `ethylacetate`, `ethylene carbonate`, `ethylene carbonate dimethyl carbonate 50:50`, `formamide`, `heptane`, `hexadecane`, `hexan-1-ol`, `hexane`, `isooctane`, `methanol`, `morpholine-4-carbaldehyde`, `nitrobenzene`, `nonane`, `octane`, `pentan-1-ol`, `pentan-2-ol`, `pentan-3-ol`, `pentane`, `phenylmethanol`, `propan-1-ol`, `propan-2-ol`, `propan-2-one`, `propane-1,2-diol`, `squalane`, `tetradecane`, `toluene`, `triethylene glycol`, `undecane`, `water`.

**Abraham only (132):** `1,1,2,2-tetrabromoethane`, `1,1,2,2-tetrachloroethane`, `1,1,2,2-tetrachloroethene`, `1,1-dibromoethane`, `1,2,3,4-tetrahydronaphthalene`, `1,2,4-trimethylbenzene`, `1,2-xylene`, `1,3-dimethylimidazolidin-2-one`, `1,3-xylene`, `1,4-xylene`, `1-bromonaphthalene`, `1-chlorobutane`, `1-chlorohexane`, `1-chloronaphthalene`, `1-ethoxy-2-(2-ethoxyethoxy)ethane`, `1-ethylpyrrolidin-2-one`, `1-iodohexadecane`, `1-methoxy-2-(2-methoxyethoxy)ethane`, `1-methylpiperidin-2-one`, `1-nitro-2-octoxybenzene`, `1-nitroethane`, `1-nitropropane`, `1-phenylethanone`, `1-propoxypropane`, `1-tert-butoxy-2-propanol`, `1H-indene`, `2,2-dichloroacetic acid`, `2-aminoethanol`, `2-butoxyethanol`, `2-ethoxyethanol`, `2-ethylhexan-1-ol`, `2-methoxy-2-methylbutane`, `2-methoxy-2-methylpropane`, `2-methylbutan-1-ol`, `2-methylpentan-1-ol`, `2-methylpyridine`, `2-phenylacetonitrile`, `2-propan-2-yloxyethanol`, `2-propoxyethanol`, `2-sulfanylethanol`, `3-(3-hydroxypropoxy)propan-1-ol`, `3-methoxybutan-1-ol`, `3-methylbutan-1-ol`, `3-methylphenol`, `3-methylthiolane 1,1-dioxide`, `4-methylpentan-2-ol`, `4-methylpentan-2-one`, `N,N-dibutylformamide`, `N,N-diethylacetamide`, `N,N-diethylethanamine`, `N,N-dimethylacetamide`, `N-ethylacetamide`, `N-ethylformamide`, `N-methylacetamide`, `N-methylformamide`, `acetonitrile_40_water_60`, `acetonitrile_60_water_40`, `aniline`, `anisole`, `benzonitrile`, `benzyl acetate`, `bromobenzene`, `bromoethane`, `bromoform`, `butan-2-one`, `butane-1,4-diol`, `butanenitrile`, `butyl acetate`, `butylbenzene`, `carbon disulfide`, `cis-1,2-dichloroethene`, `cis-decaline`, `cumene`, `cyclohexanol`, `cyclohexanone`, `cyclohexylbenzene`, `cyclohexylcyclohexane`, `dec-1-ene`, `deca-1,9-diene`, `decalin`, `decan-1-ol`, `deuterated water`, `dibutyl benzene-1,2-dicarboxylate`, `diiodomethane`, `dinonyl benzene-1,2-dicarboxylate`, `ethoxybenzene`, `ethoxyethane`, `ethyl butanoate`, `ethylbenzene`, `fluorobenzene`, `furan-2-carbaldehyde`, `furan-2-ylmethanol`, `heptan-1-ol`, `heptan-2-one`, `hex-1-ene`, `hexadec-1-ene`, `hexafluorobenzene`, `hexamethylphosphoramide`, `hexanedinitrile`, `hexyl acetate`, `iodobenzene`, `methanol_50_water_50`, `methyl acetate`, `nitromethane`, `nonan-1-ol`, `oct-1-ene`, `octamethylcyclotetrasiloxane`, `oxolan-2-one`, `pentadecane`, `pentan-2-one`, `pentane-1,5-diol`, `pentyl acetate`, `perflexane`, `perfluorooctane`, `propane-1,3-diol`, `propanenitrile`, `propyl acetate`, `pyridine`, `quinoline`, `tert-butylbenzene`, `tetraethylene glycol`, `tetraglyme`, `thiodiglycol`, `trans-1,2-dichloroethene`, `trans-decalin`, `tributyl phosphate`, `tricaprylin`, `tridecane`, `triethyl phosphate`, `triglyme`, `triolein`, `undecan-1-ol`.

**Mintz only (2):** `diethyl carbonate`, `dimethyl carbonate`.

**Complete LSERs + mapped CoolProp identity** (the independent reference-solute call at the anchor succeeds for every row):

| Solvent | CoolProp identity | Tc (K) |
| --- | --- | --- |
| 1,2-dichloroethane | Dichloroethane | 561.60000 |
| benzene | benzene | 562.02000 |
| cyclohexane | CycloHexane | 553.60000 |
| decane | decane | 617.70000 |
| dodecane | Dodecane | 658.10000 |
| ethanol | ethanol | 514.71000 |
| heptane | Heptane | 540.13000 |
| hexane | Hexane | 507.82000 |
| methanol | Methanol | 512.50000 |
| nonane | nonane | 594.55000 |
| octane | Octane | 568.74000 |
| pentane | Pentane | 469.70000 |
| propan-2-one | Acetone | 508.10000 |
| toluene | toluene | 591.75000 |
| undecane | Undecane | 638.80000 |
| water | water | 647.09600 |
<!-- END I034:inventory -->

The solvent-library keyword scan and complete label inventory find **no PS,
polymer-melt or polymer reference entry**. Long hydrocarbon solvents such as
hexadecane and squalane do exist, but they are discrete molecular solvents,
not a polymer equation of state. RMG's ring/aromatic *solute* groups must not
be confused with an aromatic-polymer *solvent* reference.

### Validity and out-of-range behavior

| Model | Evidence and nominal domain | What code actually does |
| --- | --- | --- |
| Abraham ΔG / Mintz ΔH anchor | `rmgpy/data/solvation.py`, `calc_g`, `calc_h`, `calc_s`, `get_solvation_correction`; anchor at 298 K [1] | No temperature argument. Always returns the anchor values. Calling this at another reactor temperature does not make it a temperature-dependent prediction. Missing coefficients cause arithmetic `TypeError` in the direct correction path. |
| Linear H−TS, constant ΔHsolv/ΔSsolv (ΔCpsolv = 0) | `documentation/source/users/rmg/liquids.rst`, “linear extrapolation”; documentation says a reasonable approximation only up to approximately 400 K, with increasing deviation farther from the anchor | `rmgpy/thermo/thermoengine.py`, `process_thermo_data`, shifts Wilhoit H0 and S0 without changing Cp. No Tc check, clamp, warning or high-T validity enforcement; the solvation part silently extrapolates. All requested melt temperatures are extrapolations, including those above the stand-in's real Tc. |
| Chung/Japas/Harvey K-factor | `get_T_dep_solvation_energy_from_LSER_298`, `get_Kfactor_parameters`, `get_Kfactor`; room temperature to below solvent Tc along its saturation curve [2] | Requires complete anchor coefficients and a CoolProp identity. Explicit `InputError` for T ≥ Tc; no clamp or above-Tc extrapolation. Density/property failures, including below the fluid's lower calculable range, become `DatabaseError`. There is no explicit room-temperature lower guard: sub-anchor calls succeed when saturation properties exist, but that is extrapolation, not validation. |

The K-factor model fits both the anchor value and a finite-difference anchor
gradient, matching two branches at `0.75 Tc`; it needs saturation properties
at the anchor, the adjacent temperature, the transition and its adjacent
temperature. Hence `T < Tc` is **necessary, not sufficient**. A solvent whose
properties cannot be evaluated at the fitting temperatures also fails. The
underlying saturated vapor/liquid reference is not automatically the same
pressure/reference as a reactor PS melt.

### Candidate screen fixed without validation data

The screen covers the simplest aromatic references (benzene/toluene), a
capped-repeat-unit analogue (ethylbenzene), and two nonaromatic hydrocarbon
controls (dodecane/hexadecane). Hexadecane's critical temperature is an explicitly
embedded NIST literature datum [7], not an RMG/CoolProp result. Its pinned
`name_in_coolprop` is absent; knowing Tc does not supply the missing EOS.

<!-- BEGIN I034:candidates -->
| Stand-in | LSERs | Tc (K) | CoolProp mapping | K-factor at 600/700/800 K | Linear H−TS at 600–800 K |
| --- | --- | --- | --- | --- | --- |
| benzene | G + H | 562.02000 | benzene | InputError/InputError/InputError | explicit extrapolation at all requested T |
| dodecane | G + H | 658.10000 | Dodecane | computed/InputError/InputError | explicit extrapolation at all requested T |
| ethylbenzene | G only | 617.12000 | EthylBenzene | DatabaseError/DatabaseError/DatabaseError | unavailable |
| hexadecane | G + H | 722 ± 4 [NIST] | none | DatabaseError/DatabaseError/DatabaseError | explicit extrapolation at all requested T |
| toluene | G + H | 591.75000 | toluene | InputError/InputError/InputError | explicit extrapolation at all requested T |
<!-- END I034:candidates -->

Benzene and toluene have real critical points below even the lowest requested
temperature. Ethylbenzene's lowest requested temperature is below Tc but its
enthalpy coefficients are missing, so its K-factor call still fails there.
Dodecane's K-factor is usable at the lowest requested temperature only.
Hexadecane is below its literature Tc at the lower two requested temperatures
and above it at the highest, but has no implemented temperature-dependent path
even below Tc. **No candidate's molecular K-factor path covers the full melt
window.** Switching or clamping between models across Tc would not fix that.

### Executed boundary calls

The reference solute is the first species inserted from the selected homolysis
reaction, and is the same for every stand-in. The JSON preserves its identity
and all full error messages. Direct anchor/H−TS calls on ethylbenzene fail
because Mintz coefficients are absent; the guarded K-factor wrapper reports
`DatabaseError` before reaching the critical-temperature check.

<!-- BEGIN I034:boundaries -->
Fixed reference solute: `CC(CC(CCc1ccccc1)c1ccccc1)c1ccccc1`.

| Stand-in | Call | T (K) | Observed result (one fixed reference solute) |
| --- | --- | --- | --- |
| benzene | above_Tc | 563.02000 | InputError |
| benzene | at_Tc | 562.02000 | InputError |
| benzene | below_298 | 297.00000 | Gsolv=-76.676135 kJ/mol |
| benzene | below_triple | 277.67400 | DatabaseError |
| benzene | linear_at_800 | 800.00000 | Gsolv=-6.797573 kJ/mol |
| benzene | near_Tc | 562.01000 | Gsolv=-1.104339 kJ/mol |
| dodecane | above_Tc | 659.10000 | InputError |
| dodecane | at_Tc | 658.10000 | InputError |
| dodecane | below_298 | 297.00000 | Gsolv=-66.817700 kJ/mol |
| dodecane | below_triple | 262.60000 | DatabaseError |
| dodecane | linear_at_800 | 800.00000 | Gsolv=+15.029375 kJ/mol |
| dodecane | near_Tc | 658.09000 | Gsolv=-2.152120 kJ/mol |
| ethylbenzene | above_Tc | 618.12000 | DatabaseError |
| ethylbenzene | at_Tc | 617.12000 | DatabaseError |
| ethylbenzene | below_298 | 297.00000 | DatabaseError |
| ethylbenzene | below_triple | 177.20000 | DatabaseError |
| ethylbenzene | linear_at_800 | 800.00000 | TypeError |
| ethylbenzene | near_Tc | 617.11000 | DatabaseError |
| hexadecane | 298 | 298.00000 | DatabaseError |
| hexadecane | linear_at_800 | 800.00000 | Gsolv=+10.627920 kJ/mol |
| toluene | above_Tc | 592.75000 | InputError |
| toluene | at_Tc | 591.75000 | InputError |
| toluene | below_298 | 297.00000 | Gsolv=-76.068374 kJ/mol |
| toluene | below_triple | 177.00000 | DatabaseError |
| toluene | linear_at_800 | 800.00000 | Gsolv=-5.509485 kJ/mol |
| toluene | near_Tc | 591.74000 | Gsolv=-0.786598 kJ/mol |
<!-- END I034:boundaries -->

The near-Tc function is not a divergent ΔG pole: its high-T branch drives
K-factor toward unity as the saturated liquid density approaches the critical
density, and ΔGsolv toward zero as liquid and vapor merge. Derivatives and
critical fluid properties can be nonanalytic/ill-conditioned. Those features,
and the code's explicit cutoff, belong to the molecular stand-in, not the melt.
A critical-point “collapse of solvation” is not a physical prediction for PS.

## 2. Solutes, baseline and reaction definitions

### Descriptor estimation

`get_solute_data` first attempts a matching solute library. For radicals it
tries a hydrogen-saturated library species plus site-specific HBI radical
corrections; otherwise it uses group additivity on the saturated structure,
then removes the saturating H and applies the radical group corrections.
The fallback includes main heavy-atom groups, halogen corrections when
applicable, acyclic/cyclic long-distance interactions, and ring/polycyclic
corrections. Aromatic atom/group matching and the benzene ring correction are
present. **However, successful API calls do not imply that a carbon-radical
correction was available: the missing correction is silently skipped here.**
McGowan volume is set from atom sizes and bond counts, not fitted per reaction.
These are the algorithms in the pinned code; [1] describes the related solute
group-contribution method, not a validation of PS radicals at melt temperatures.

The baseline H-abstraction in the prior probe is **radical + styrene**, producing
an aryl radical. It is not a chain-to-chain benzylic transfer. To answer the
actual chain-transfer question without losing the prior numerical cross-check,
this probe additionally generates end-radical + pristine-chain H-abstraction,
selects a benzylic product radical, and chooses the lexicographically smallest
product-SMILES pair satisfying that criterion. It is deliberately not selected
for a favorable Kc or correction. Both H reactions are reported throughout.

The propagation reaction is reconstructed from the earliest ceiling-crossing
pair in an addition-only real-compiler artifact containing the end-radical and
end-radical + styrene proxies, with exactly their bounded participants. These
are the only declared bimolecular propagation inputs; the gas ceiling assertion
confirms that omitting irrelevant families does not change the prior crossing.
Do not replace its short end-radical proxy with an arbitrary longer chain and
claim the prior ceiling has been reproduced.

<!-- BEGIN I034:reactions -->
**H (prior styrene pair):** `[CH2]C(CCc1ccccc1)c1ccccc1 + C=Cc1ccccc1 -> CC(CCc1ccccc1)c1ccccc1 + C=Cc1[c]cccc1` (Δn = +0).

**H (chain-to-chain):** `[CH2]C(CC(CCc1ccccc1)c1ccccc1)c1ccccc1 + CC(CC(CCc1ccccc1)c1ccccc1)c1ccccc1 -> CC(CC(CCc1ccccc1)c1ccccc1)c1ccccc1 + CC(CC(C[CH]c1ccccc1)c1ccccc1)c1ccccc1` (Δn = +0).

**homolysis:** `CC(CC(CCc1ccccc1)c1ccccc1)c1ccccc1 -> [CH](CCc1ccccc1)c1ccccc1 + [CH2]C(C)c1ccccc1` (Δn = +1).

**propagation:** `[CH2]C(CCc1ccccc1)c1ccccc1 + C=Cc1ccccc1 -> [CH2]C(CC(CCc1ccccc1)c1ccccc1)c1ccccc1` (Δn = -1).
<!-- END I034:reactions -->

<!-- BEGIN I034:solutes -->
Selected unique species: **9**; failures: **0**. Default descriptor lookup and forced group additivity both succeed on every selected species.

| Species (SMILES) | S | B | E | L | A | V | Default origin |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `C=Cc1[c]cccc1` | 0.616880 | 0.185010 | 0.835790 | 4.102870 | 0.000000 | 0.933700 | group additivity |
| `C=Cc1ccccc1` | 0.616880 | 0.185010 | 0.835790 | 4.102870 | 0.000000 | 0.955200 | group additivity |
| `CC(CC(CCc1ccccc1)c1ccccc1)c1ccccc1` | 1.637730 | 0.524030 | 1.866740 | 11.892790 | 0.000000 | 2.777400 | group additivity |
| `CC(CC(C[CH]c1ccccc1)c1ccccc1)c1ccccc1` | 1.637730 | 0.524030 | 1.866740 | 11.892790 | 0.000000 | 2.755900 | group additivity |
| `CC(CCc1ccccc1)c1ccccc1` | 1.089980 | 0.345720 | 1.244130 | 7.920730 | 0.000000 | 1.887800 | group additivity |
| `[CH2]C(C)c1ccccc1` | 0.547750 | 0.178310 | 0.622610 | 4.281160 | 0.000000 | 1.117600 | group additivity |
| `[CH2]C(CC(CCc1ccccc1)c1ccccc1)c1ccccc1` | 1.637730 | 0.524030 | 1.866740 | 11.892790 | 0.000000 | 2.755900 | group additivity |
| `[CH2]C(CCc1ccccc1)c1ccccc1` | 1.089980 | 0.345720 | 1.244130 | 7.920730 | 0.000000 | 1.866300 | group additivity |
| `[CH](CCc1ccccc1)c1ccccc1` | 1.084460 | 0.334820 | 1.243040 | 7.588240 | 0.000000 | 1.725400 | group additivity |
<!-- END I034:solutes -->

### Silent carbon-radical correction gap

Both displayed H transfers have vanishing reaction solvation correction to
floating-point precision. This was investigated rather than interpreted as
physical evidence. The pinned `input/solvation/groups/radical.py` mainly corrects
the H-bonding descriptor A for oxygen/nitrogen radicals; carbon radicals descend
to `R_rad`, whose data is absent and which has no ancestor with data.
An explicit correction lookup raises `KeyError`, but
`estimate_radical_solute_data_via_hbi` catches that exception and proceeds.

<!-- BEGIN I034:radicals -->
Selected carbon-radical species lacking a tabulated radical correction: **6**.

| Radical (SMILES) | Saturated analogue | Matched node | Explicit correction lookup | max abs(ΔS,ΔB,ΔE,ΔL,ΔA) | ΔV |
| --- | --- | --- | --- | --- | --- |
| `C=Cc1[c]cccc1` | `C=Cc1ccccc1` | R_rad | KeyError: no parent with data | 0.000e+00 | -0.021500 |
| `CC(CC(C[CH]c1ccccc1)c1ccccc1)c1ccccc1` | `CC(CC(CCc1ccccc1)c1ccccc1)c1ccccc1` | R_rad | KeyError: no parent with data | 0.000e+00 | -0.021500 |
| `[CH2]C(C)c1ccccc1` | `CC(C)c1ccccc1` | R_rad | KeyError: no parent with data | 0.000e+00 | -0.021500 |
| `[CH2]C(CC(CCc1ccccc1)c1ccccc1)c1ccccc1` | `CC(CC(CCc1ccccc1)c1ccccc1)c1ccccc1` | R_rad | KeyError: no parent with data | 0.000e+00 | -0.021500 |
| `[CH2]C(CCc1ccccc1)c1ccccc1` | `CC(CCc1ccccc1)c1ccccc1` | R_rad | KeyError: no parent with data | 0.000e+00 | -0.021500 |
| `[CH](CCc1ccccc1)c1ccccc1` | `c1ccc(CCCc2ccccc2)cc1` | R_rad | KeyError: no parent with data | 0.000e+00 | -0.021500 |
<!-- END I034:radicals -->

Thus each selected carbon radical inherits the saturated analogue's S/B/E/L/A.
Its V changes because an H atom is removed, but V is not used by these
gas-to-solvent LSERs. For these H transfers the saturated-analogue descriptors
and solvent intercepts cancel on the two sides. The K-factor model likewise
receives identical anchor inputs for radical and saturated analogue. The zero
correction is consequently an **estimator limitation/cancellation**, not a
measured absence of radical-solvent effects. The diagnostic correction-lookup
failures above are distinct from the absence of fatal descriptor API failures.
Neither fitting new carbon-radical increments nor changing the estimator is
authorized by this task.

Full descriptor/library/HBI/ring comments and forced-group estimates are
preserved for every species in JSON. Success proves descriptor coverage of
these selected reaction participants only, **not** the accuracy of those
descriptors, and not coverage of every reversible reaction in a future event
inventory. The earlier bounded-set coverage claim is not substituted for an
actual coverage measurement here.

### Gas-phase numbers reproduced first

<!-- BEGIN I034:baseline -->
| T (K) | Kc,association (m³/mol) | k_homolysis (s⁻¹) | Kc,H prior | Kc,H chain | Kc,propagation (m³/mol) |
| --- | --- | --- | --- | --- | --- |
| 600 | 1.168950e+18 | 1.875979e-13 | 3.027770e-05 | 7.581626e+02 | 9.707710e-03 |
| 700 | 1.631045e+14 | 2.024922e-09 | 1.258300e-04 | 1.588368e+02 | 1.191747e-03 |
| 800 | 2.175175e+11 | 2.043769e-06 | 3.650914e-04 | 4.931260e+01 | 2.558179e-04 |
<!-- END I034:baseline -->

The homolysis is displayed as closed-shell chain → two radicals. The prior
backbone Kc number uses the **opposite association direction**; it is the
reciprocal of the homolysis Kc. Both conventions are retained explicitly.

<!-- BEGIN I034:gas_magnitudes -->
| Reaction | T (K) | ΔrGgas (kJ/mol) | Kc,gas (SI concentration convention) |
| --- | --- | --- | --- |
| H (prior styrene pair) | 600 | 51.907743 | 3.027770e-05 |
| H (chain-to-chain) | 600 | -33.079449 | 7.581626e+02 |
| homolysis | 600 | 222.498441 | 8.554688e-19 |
| propagation | 600 | 8.165660 | 9.707710e-03 |
| H (prior styrene pair) | 700 | 52.268140 | 1.258300e-04 |
| H (chain-to-chain) | 700 | -29.495706 | 1.588368e+02 |
| homolysis | 700 | 207.017729 | 6.131040e-15 |
| propagation | 700 | 22.631506 | 1.191747e-03 |
| H (prior styrene pair) | 800 | 52.649650 | 3.650914e-04 |
| H (chain-to-chain) | 800 | -25.929044 | 4.931260e+01 |
| homolysis | 800 | 191.670939 | 4.597330e-12 |
| propagation | 800 | 36.987619 | 2.558179e-04 |
<!-- END I034:gas_magnitudes -->

These ΔrGgas values come from RMG gas thermo. Dimensional Kc is not simply
`exp(−ΔrGgas/RT)` for a reaction changing molecule count. RMG
`Reaction.get_equilibrium_constant(type="Kc")` uses
`Kc = exp(−ΔrGgas/RT) [P0/(RT)]^Δn`, with `P0 = 1e5 Pa` and SI concentration
units. Homolysis Kc has mol/m³ units; association and propagation Kc have
m³/mol units; the H-transfer Kc values are dimensionless. The verifier checks
this relation independently.

## 3. Corrections and Kc shifts

At each temperature and **for the displayed reaction direction**, define

```text
ΔΔGsolv(T) = Σproducts Gsolv,i(T) − Σreactants Gsolv,i(T)
Kc,melt(T) / Kc,gas(T) = exp[−ΔΔGsolv(T)/(RT)]
```

This is the RMG equimolar gas-to-solvent (Ben-Naim) solvation convention
documented in `liquids.rst`. Keep the inherited concentration Kc convention;
do not apply an additional pressure-to-concentration standard-state conversion
to the solvation correction. Species free-energy corrections are applied
consistently to all participants, including the radical proxies and styrene.
The same species convention also guarantees reciprocal corrections for the
reverse reaction. Chain activities and finite-concentration activity
coefficients are **not** supplied by this stand-in calculation.

### Anchor reaction corrections

<!-- BEGIN I034:anchors -->
| Stand-in | Reaction | ΔΔG298 (kJ/mol) | ΔΔH298 (kJ/mol) | ΔΔS298 (J/mol/K) |
| --- | --- | --- | --- | --- |
| benzene | H (prior styrene pair) | 0.000000 | 0.000000 | 0.000000 |
| benzene | H (chain-to-chain) | 0.000000 | 0.000000 | 0.000000 |
| benzene | homolysis | -0.483730 | -4.271641 | -12.711111 |
| benzene | propagation | 1.452531 | 5.285489 | 12.862272 |
| dodecane | H (prior styrene pair) | 0.000000 | 0.000000 | 0.000000 |
| dodecane | H (chain-to-chain) | 0.000000 | 0.000000 | 0.000000 |
| dodecane | homolysis | -0.554378 | -6.036319 | -18.395777 |
| dodecane | propagation | 1.391544 | 7.140420 | 19.291532 |
| hexadecane | H (prior styrene pair) | 0.000000 | 0.000000 | 0.000000 |
| hexadecane | H (chain-to-chain) | 0.000000 | 0.000000 | 0.000000 |
| hexadecane | homolysis | 0.094409 | -4.293206 | -14.723541 |
| hexadecane | propagation | 0.790138 | 5.466959 | 15.694030 |
| toluene | H (prior styrene pair) | 0.000000 | 0.000000 | 0.000000 |
| toluene | H (chain-to-chain) | 0.000000 | 0.000000 | 0.000000 |
| toluene | homolysis | -0.291408 | -5.353525 | -16.986971 |
| toluene | propagation | 1.208457 | 6.393069 | 17.398028 |
<!-- END I034:anchors -->

For ethylbenzene only the Gibbs anchor is evaluable; no enthalpy/entropy or
high-T correction is manufactured from it. Keeping that Gibbs value constant
at all requested temperatures would introduce an additional unsupported model,
not recover the implemented Abraham/Mintz temperature dependence.

<!-- BEGIN I034:ethylbenzene_anchor -->
| Reaction | Ethylbenzene ΔΔG298 (kJ/mol) |
| --- | --- |
| H (prior styrene pair) | 0.000000 |
| H (chain-to-chain) | 0.000000 |
| homolysis | -1.049096 |
| propagation | 2.030173 |
<!-- END I034:ethylbenzene_anchor -->

### Linear H−TS: explicit high-temperature extrapolation

<!-- BEGIN I034:linear_shifts -->
| Stand-in | Reaction | T (K) | ΔΔGsolv (kJ/mol) | Kc,melt/Kc,gas | Kc,melt (SI) | abs(ΔΔGsolv)/abs(ΔrGgas) |
| --- | --- | --- | --- | --- | --- | --- |
| benzene | H (prior styrene pair) | 600 | +0.000000 | 1.000000e+00 | 3.027770e-05 | 0.000000 |
| benzene | H (chain-to-chain) | 600 | +0.000000 | 1.000000e+00 | 7.581626e+02 | 0.000000 |
| benzene | homolysis | 600 | +3.355025 | 5.104170e-01 | 4.366458e-19 | 0.015079 |
| benzene | propagation | 600 | -2.431875 | 1.628205e+00 | 1.580614e-02 | 0.297817 |
| benzene | H (prior styrene pair) | 700 | +0.000000 | 1.000000e+00 | 1.258300e-04 | 0.000000 |
| benzene | H (chain-to-chain) | 700 | +0.000000 | 1.000000e+00 | 1.588368e+02 | 0.000000 |
| benzene | homolysis | 700 | +4.626136 | 4.516485e-01 | 2.769075e-15 | 0.022347 |
| benzene | propagation | 700 | -3.718102 | 1.894272e+00 | 2.257493e-03 | 0.164289 |
| benzene | H (prior styrene pair) | 800 | +0.000000 | 1.000000e+00 | 3.650914e-04 | 0.000000 |
| benzene | H (chain-to-chain) | 800 | +0.000000 | 1.000000e+00 | 4.931260e+01 | 0.000000 |
| benzene | homolysis | 800 | +5.897248 | 4.120569e-01 | 1.894362e-12 | 0.030768 |
| benzene | propagation | 800 | -5.004329 | 2.121986e+00 | 5.428418e-04 | 0.135297 |
| dodecane | H (prior styrene pair) | 600 | +0.000000 | 1.000000e+00 | 3.027770e-05 | 0.000000 |
| dodecane | H (chain-to-chain) | 600 | +0.000000 | 1.000000e+00 | 7.581626e+02 | 0.000000 |
| dodecane | homolysis | 600 | +5.001147 | 3.669615e-01 | 3.139241e-19 | 0.022477 |
| dodecane | propagation | 600 | -4.434499 | 2.432481e+00 | 2.361382e-02 | 0.543067 |
| dodecane | H (prior styrene pair) | 700 | +0.000000 | 1.000000e+00 | 1.258300e-04 | 0.000000 |
| dodecane | H (chain-to-chain) | 700 | +0.000000 | 1.000000e+00 | 1.588368e+02 | 0.000000 |
| dodecane | homolysis | 700 | +6.840725 | 3.087091e-01 | 1.892708e-15 | 0.033044 |
| dodecane | propagation | 700 | -6.363652 | 2.984364e+00 | 3.556605e-03 | 0.281186 |
| dodecane | H (prior styrene pair) | 800 | +0.000000 | 1.000000e+00 | 3.650914e-04 | 0.000000 |
| dodecane | H (chain-to-chain) | 800 | +0.000000 | 1.000000e+00 | 4.931260e+01 | 0.000000 |
| dodecane | homolysis | 800 | +8.680302 | 2.711728e-01 | 1.246671e-12 | 0.045288 |
| dodecane | propagation | 800 | -8.292805 | 3.478992e+00 | 8.899883e-04 | 0.224205 |
| hexadecane | H (prior styrene pair) | 600 | +0.000000 | 1.000000e+00 | 3.027770e-05 | 0.000000 |
| hexadecane | H (chain-to-chain) | 600 | +0.000000 | 1.000000e+00 | 7.581626e+02 | 0.000000 |
| hexadecane | homolysis | 600 | +4.540919 | 4.024260e-01 | 3.442629e-19 | 0.020409 |
| hexadecane | propagation | 600 | -3.949459 | 2.207109e+00 | 2.142598e-02 | 0.483667 |
| hexadecane | H (prior styrene pair) | 700 | +0.000000 | 1.000000e+00 | 1.258300e-04 | 0.000000 |
| hexadecane | H (chain-to-chain) | 700 | +0.000000 | 1.000000e+00 | 1.588368e+02 | 0.000000 |
| hexadecane | homolysis | 700 | +6.013273 | 3.558716e-01 | 2.181863e-15 | 0.029047 |
| hexadecane | propagation | 700 | -5.518862 | 2.581155e+00 | 3.076082e-03 | 0.243857 |
| hexadecane | H (prior styrene pair) | 800 | +0.000000 | 1.000000e+00 | 3.650914e-04 | 0.000000 |
| hexadecane | H (chain-to-chain) | 800 | +0.000000 | 1.000000e+00 | 4.931260e+01 | 0.000000 |
| hexadecane | homolysis | 800 | +7.485627 | 3.245255e-01 | 1.491951e-12 | 0.039055 |
| hexadecane | propagation | 800 | -7.088265 | 2.902730e+00 | 7.425702e-04 | 0.191639 |
| toluene | H (prior styrene pair) | 600 | -0.000000 | 1.000000e+00 | 3.027770e-05 | 0.000000 |
| toluene | H (chain-to-chain) | 600 | +0.000000 | 1.000000e+00 | 7.581626e+02 | 0.000000 |
| toluene | homolysis | 600 | +4.838658 | 3.791108e-01 | 3.243175e-19 | 0.021747 |
| toluene | propagation | 600 | -4.045748 | 2.250123e+00 | 2.184354e-02 | 0.495459 |
| toluene | H (prior styrene pair) | 700 | -0.000000 | 1.000000e+00 | 1.258300e-04 | 0.000000 |
| toluene | H (chain-to-chain) | 700 | +0.000000 | 1.000000e+00 | 1.588368e+02 | 0.000000 |
| toluene | homolysis | 700 | +6.537355 | 3.252271e-01 | 1.993980e-15 | 0.031579 |
| toluene | propagation | 700 | -5.785550 | 2.702179e+00 | 3.220313e-03 | 0.255641 |
| toluene | H (prior styrene pair) | 800 | -0.000000 | 1.000000e+00 | 3.650914e-04 | 0.000000 |
| toluene | H (chain-to-chain) | 800 | +0.000000 | 1.000000e+00 | 4.931260e+01 | 0.000000 |
| toluene | homolysis | 800 | +8.236052 | 2.899026e-01 | 1.332778e-12 | 0.042970 |
| toluene | propagation | 800 | -7.525353 | 3.099881e+00 | 7.930049e-04 | 0.203456 |
<!-- END I034:linear_shifts -->

The “melt” label in these tables means the chosen surrogate/reference model,
not a measured PS equilibrium constant. For homolysis a negative correction
stabilizes the separated radicals relative to the parent and increases the
homolysis Kc; the association Kc shift is reciprocal. If the forward rate is
held fixed, its reverse rate changes by the inverse of the displayed ratio.
For computing homolysis as the reverse of association, the homolysis rate
instead changes by the displayed **homolysis-direction** ratio. This corrects
equilibrium thermodynamics only; it does not calculate a solvent correction
to the forward barrier or a radical cage escape probability.

### K-factor: usable points only, without clamping

<!-- BEGIN I034:kfactor_shifts -->
| Stand-in | Reaction | T (K) | ΔΔGsolv (kJ/mol) | Kc,melt/Kc,gas | Kc,melt (SI) | abs(ΔΔGsolv)/abs(ΔrGgas) |
| --- | --- | --- | --- | --- | --- | --- |
| dodecane | H (prior styrene pair) | 600 | +0.000000 | 1.000000e+00 | 3.027770e-05 | 0.000000 |
| dodecane | H (chain-to-chain) | 600 | +0.000000 | 1.000000e+00 | 7.581626e+02 | 0.000000 |
| dodecane | homolysis | 600 | +1.759378 | 7.028063e-01 | 6.012289e-19 | 0.007907 |
| dodecane | propagation | 600 | -1.364836 | 1.314671e+00 | 1.276245e-02 | 0.167143 |
<!-- END I034:kfactor_shifts -->

The candidate and boundary tables give the failed calls as well as the
computed ones. There are no omitted “melt values” at failed temperatures.
Using dodecane does not solve the target's high-temperature validity problem;
the calculation is a diagnostic of the molecular-solvent model at its one
available requested-temperature point, not an all-window sensitivity arm.

### Comparison with gas-phase magnitudes

The complete gas ΔrG and Kc tables above allow comparisons for every stand-in.
For the recommendation the absolute free-energy fraction is shown explicitly:

<!-- BEGIN I034:magnitude_comparison -->
| Toluene / linear H−TS, reaction | T (K) | abs(ΔΔGsolv) / abs(ΔrGgas) | Kc,melt / Kc,gas |
| --- | --- | --- | --- |
| H (prior styrene pair) | 600 | 0.000000 | 1.000000e+00 |
| H (chain-to-chain) | 600 | 0.000000 | 1.000000e+00 |
| homolysis | 600 | 0.021747 | 3.791108e-01 |
| propagation | 600 | 0.495459 | 2.250123e+00 |
| H (prior styrene pair) | 700 | 0.000000 | 1.000000e+00 |
| H (chain-to-chain) | 700 | 0.000000 | 1.000000e+00 |
| homolysis | 700 | 0.031579 | 3.252271e-01 |
| propagation | 700 | 0.255641 | 2.702179e+00 |
| H (prior styrene pair) | 800 | 0.000000 | 1.000000e+00 |
| H (chain-to-chain) | 800 | 0.000000 | 1.000000e+00 |
| homolysis | 800 | 0.042970 | 2.899026e-01 |
| propagation | 800 | 0.203456 | 3.099881e+00 |
<!-- END I034:magnitude_comparison -->

Homolysis remains dominated by the gas-phase bond/radical thermochemistry,
yet the positive surrogate correction suppresses its equilibrium constant and,
when derived as the reverse of association, its homolysis rate. The propagation
correction is a larger fraction of its gas free-energy difference and shifts
the crossing appreciably. The unchanged H Kc values are not evidence of melt
insensitivity: they inherit the silent carbon-radical correction gap above.

### Ceiling temperatures and interpretation

The compiler crossing is a linear interpolation of
`ln(kprop[M]/kdeprop)` on its existing grid. The continuous ceiling solves
`Kc,prop(T) [M] = 1` using the same participants and concentration; it is also
reported to avoid misrepresenting interpolation as an exact thermodynamic
root. A ceiling needs a monomer concentration/activity: it is not a solvent's
intrinsic temperature property. The gas and corrected comparisons hold that
input fixed. No monomer fugacity, reactor density or polymer activity is fitted.

<!-- BEGIN I034:ceilings -->
| Reference / model | Compiler-grid crossing (K) | Continuous Kc[M]=1 root (K) | Continuous shift from gas (K) |
| --- | --- | --- | --- |
| gas | 710.248673 | 710.020462 | 0.000000 |
| benzene / linear H−TS | 753.636245 | 753.526132 | +43.505670 |
| dodecane / linear H−TS | 790.590329 | 790.383770 | +80.363308 |
| hexadecane / linear H−TS | 776.583313 | 776.532314 | +66.511852 |
| toluene / linear H−TS | 781.491601 | 781.324909 | +71.304447 |
<!-- END I034:ceilings -->

None of these molecular K-factor models supplies a physically applicable
ceiling for the entire melt interval. Benzene/toluene fail the temperature
domain, ethylbenzene lacks enthalpy data, and the dodecane molecular reference
ends below the relevant linear-model crossing. Hexadecane has no CoolProp
mapping. No above-Tc continuation or stitched K-factor/linear ceiling is used.

## 4. Physical case and literature evidence

**Not established by the retrieved literature:** absolute solvation free energies
of these PS-chain carbon radicals in a neat PS melt over the requested
temperature interval; an uncertainty bound for transferring molecular-solvent
LSERs to that melt; or a radical-solvation dataset that validates the proposed
linear extrapolation. Retrieved transport/structure/spin-probe studies provide
chemical context, **not those missing thermodynamic values**.

| Candidate | Chemical case for a PS segment | What the evidence does and does not establish |
| --- | --- | --- |
| Toluene | Alkylated aromatic hydrocarbon; retains aromatic dispersion/polarisability without introducing strong polarity or H-bonding. Less alkyl backbone content than a capped repeat. | The toluene/PS oligomer study [4] directly compares solution and melt packing and finds only limited conformational changes upon solvation. That supports local chemical compatibility, not equality of transfer free energies. Complete pinned LSERs make it executable without new fitted parameters. |
| Benzene | Aromatic, non-H-bonding hydrocarbon; omits the saturated backbone/alkyl substitution, emphasizing the phenyl environment. | An intentional aromatic-composition sensitivity arm. Its structural similarity is a chemical inference, not proof that its bulk liquid is a PS melt. The modeled differences from toluene test one surrogate uncertainty, not the entire melt uncertainty. |
| Ethylbenzene | A hydrogen-capped PS repeat unit has ethylbenzene connectivity: one phenyl ring plus the saturated two-carbon segment. Better composition match among this set. | PS/ethylbenzene atomistic/coarse-grained simulations and NMR measurements [3] support a meaningful polymer/penetrant environment and demonstrate transport differences. This does not establish high-T carbon-radical solvation. Missing pinned Mintz coefficients rule it out for the specified models; no substituted coefficients are proposed. |
| Dodecane | Nonpolar and dispersive, but no aromatic ring or aromatic electronic response. A hydrocarbon contrast rather than a close PS reference. | Its larger Tc buys only one available requested-temperature K-factor point. The absence of phenyl content is a reason not to rank it above toluene merely because a numerical correction is smaller or a temperature limit is higher. |
| Hexadecane | Long nonpolar hydrocarbon; more segment-like size/cohesion than the small solvents, but again no aromatic content. | Complete anchor LSERs allow a linear diagnostic without CoolProp. The longer hydrocarbon does not become a polymer EOS, and its literature critical temperature still lies inside the target window [7]. |

Small-molecule solubility in PS melts is commonly represented by Henry
coefficients with temperature dependence: [5] explicitly models dissolved
gases/co-blowing agents in the melt. This shows the need for a polymer-specific
chemical-potential/solubility reference and pressure dependence; its fitted
coefficients are **not imported or used here**. Diffusion measurements in [3]
likewise do not determine equilibrium ΔG by themselves. A molecular stand-in
cannot simply substitute its saturated fluid density for a PS melt density.

For radicals, [6] studies stable nitroxide spin labels in PS and distinguishes
the local mobility/free-volume environments of chain ends and interior sites.
This is evidence that the polymer microenvironment is not a single uniform
molecular-solvent cage. It is **not** a carbon-radical solvation-free-energy
measurement; reorientational mobility is not an equilibrium stabilization
energy. No radical lifetime, cage parameter or chain-end/interior energy is
extracted from that paper for this probe. RMG radical HBI groups supply the
numerical radical descriptors here; their presence is not validation in PS.

### Recommendation, uncertainty and falsification

Fix **toluene / linear H−TS** a priori if the owner accepts a provisional
segment-reference approximation. It uses both independent anchor quantities,
preserves thermodynamic reciprocity, needs no invented enthalpy/EOS parameter,
and does not force a spurious molecular critical-point collapse into the melt.
The numerical shifts and ceilings are consequences of this chemically chosen
model; **they are not why it is selected**. Keep benzene / linear H−TS as the
single primary sensitivity arm, with otherwise identical species and standard
states. The hydrocarbon controls above are diagnostic comparisons, not added
fit knobs.

The main uncertainty is not numerical precision: it is neglect of solvation
heat-capacity curvature far from the anchor, molecular-solvent to polymer
transferability, finite proxy end/capping effects, radical HBI accuracy, and
polymer/monomer activities and melt density. In particular, the carbon-radical
correction gap is demonstrated, not just an unspecified HBI error bar.
Cancellation between reactants and
products helps extensive descriptor terms, but does not cancel solvent
intercepts when molecule count changes, nor guarantee cancellation of cavity
and connectivity entropy. The benzene/toluene difference is **not a confidence
interval**, and an absolute high-temperature melt error cannot be bounded from
the evidence retrieved here.

**Kill-test of the recommendation:** independently obtained carbon-radical and
closed-shell chemical potentials in an explicit PS melt, on these same proxy
structures and stated standard states, could reverse the reaction correction
or reveal significant temperature curvature. Such evidence would invalidate
the proposed ΔCpsolv = 0 segment reference rather than license adjusting its
coefficients against validation yields. No such expensive calculation is
performed in this diagnosis-only task. If the owner requires demonstrated
validity over the entire window rather than a transparent stand-in, **none of
the screened pinned models qualifies**; a polymer-specific a-priori model is
remaining scientific work, not a hidden claim that this probe completed it.

## 5. Contract corrections and limits

- Critical temperature alone is not the blocker for ethylbenzene: the pinned
  solvent also lacks its Mintz coefficients. Its lowest requested temperature
  is below Tc, yet the wrapper returns `DatabaseError`.
- A CoolProp identity is required for K-factor, not for every solvation model.
  The standard thermoengine already implements a documented linear H−TS
  extrapolation without an EOS or Tc guard. It can return numbers throughout
  the window, but documentation does not validate them there.
- The prior H-abstraction baseline is not a chain-to-chain transfer. Both that
  exact baseline and an independently, structurally selected chain example
  are retained, so reproducing one is not passed off as probing the other.
- The prior ceiling is a compiler-grid interpolation at the configured monomer
  concentration, not a concentration-independent or exact continuous root.
  Both definitions are reproduced separately and corrected consistently.
- The K-factor ΔG itself is not singular at Tc. Its limiting fluid reference
  collapses toward zero solvation, and the implementation rejects the endpoint
  and higher temperatures. This is a stand-in limitation, not a PS critical
  point.
- Descriptor success is bounded to the displayed pairs. Neither a full future
  reversible-inventory coverage result nor validated melt accuracy follows
  from these successful calls.
- In particular, all selected carbon-radical correction lookups fail at the
  generic `R_rad` node and the public estimator silently continues. The reported
  zero H-transfer solvation shifts are explained by saturated-analogue descriptor
  cancellation, not supported radical thermodynamics in PS.

## Sources and retrieval scope

All sources are primary papers, primary publisher abstracts, RMG source/docs
or the NIST thermophysical compilation. Retrieved on 2026-09-30. Searches were
for the RMG temperature model, PS/ethylbenzene transport, PS/toluene structure,
PS small-molecule solubility and radical spin-probe environments. Bibliographic
search results mentioning excluded validation topics were not followed or used
as evidence. The restricted thermophysical/transport literature search is not
a claim of exhaustive absence of radical-solvation data. No experimental
reaction-validation dataset or campaign result is consulted.

1. Y. Chung, F. H. Vermeire, H. Wu, P. J. Walker, M. H. Abraham, W. H. Green,
   “Group Contribution and Machine Learning Approaches to Predict Abraham
   Solute Parameters, Solvation Free Energy, and Solvation Enthalpy,”
   *J. Chem. Inf. Model.* **62** (2022), 433–446.
   DOI: `10.1021/acs.jcim.1c01103`. Publisher abstract and pinned library
   documentation retrieved; not used to assert high-temperature melt accuracy.
2. Y. Chung, R. J. Gillis, W. H. Green, “Temperature-dependent vapor–liquid
   equilibria and solvation free energy estimation from minimal data,”
   *AIChE J.* **66** (2020), e16976. DOI: `10.1002/aic.16976`.
   Publisher abstract and RMG implementation/documentation retrieved. The
   proposed model is for dilute molecular-fluid solutes below solvent Tc.
3. V. A. Harmandaris et al., “Ethylbenzene Diffusion in Polystyrene: United Atom
   Atomistic/Coarse Grained Simulations and Experiments,” *Macromolecules*
   **40** (2007), 7026–7035. DOI: `10.1021/ma070201o`.
   Publisher abstract retrieved; no experimental table or polymer distribution
   is used. Supports polymer/penetrant transport context, not radical ΔG.
4. B. Bayramoglu, R. Faller, “Structural properties of polystyrene oligomers in
   different environments: a molecular dynamics study,” *Phys. Chem. Chem.
   Phys.* **13** (2011), 18107–18114. DOI: `10.1039/C1CP21724K`.
   Publisher abstract retrieved; supports packing/structural differences.
5. R. Breuer et al., “Modeling flow and cell formation in foam sheet extrusion
   of polystyrene with CO2 and co-blowing agents. Part I: Material model,”
   *Polym. Eng. Sci.* **61** (2021), 2799–2813.
   DOI: `10.1002/pen.25800`. Publisher's small-molecule solubility discussion
   retrieved; no coefficients, fitted data, foam-output data or reaction yields
   are imported into the probe or used to rank stand-ins.
6. Y. Miwa et al., “Influence of Chain End and Molecular Weight on Molecular
   Motion of Polystyrene, Revealed by the ESR Selective Spin-Label Method,”
   *Macromolecules* **36** (2003), 3235–3239.
   DOI: `10.1021/ma030026l`. Publisher abstract only; no MWD dataset is opened,
   searched or used. Supports local spin-label environments, not melt-radical
   solvation energies.
7. NIST Chemistry WebBook, SRD 69, **hexadecane**, CAS `544-76-3`, phase-change
   table, critical-temperature average. Bibliographic source: NIST entry
   `C544763`, phase-change data. The external Tc datum is explicitly embedded
   with attribution in `run_probe.py`; the script does not fetch literature or
   depend on network access during reproduction.

Local implementation evidence: `rmgpy/data/solvation.py`,
`rmgpy/thermo/thermoengine.py`, `rmgpy/reaction.py`,
`documentation/source/users/rmg/liquids.rst`, `rmgpy/kmc/compiler.py`, and the
pinned `input/solvation/libraries/solvent.py`. Numeric tables and all code-based
claims are reproducible without network access.
