# I-032 — why the compiled PS event set has no initiation or chain H-abstraction

## Scope and reproduction

This is a diagnosis of RMG-Py `0a56f94a68a463adf438a365d9f932df83acf326`
against RMG-database `4a12d36fcdc193ede82c8d1ab5c1653495d445bc`. No
product code or database file was changed. The full probe covered 195 unique
species and found no gas-thermo or solute-descriptor estimation failure on that
bounded set.

Run both committed scripts from the repository root:

```bash
mkdir -p /tmp/i032-reproduce
MPLCONFIGDIR=/tmp/i032-reproduce/mpl PYTHONPATH=$PWD \
  /home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/fixtures/i032_probe/run_probe.py \
  --database /home/alon/Code/RMG-database \
  --output /tmp/i032-reproduce/results.json \
  > >(tee -a /tmp/i032-reproduce/stdout.log) \
  2> >(tee -a /tmp/i032-reproduce/stderr.log >&2)
PYTHONPATH=$PWD /home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/fixtures/i032_probe/verify_results.py \
  /tmp/i032-reproduce/results.json
```

The reproduced verifier output is:

```text
I-032 report claims reproduced from run_probe.py output
```

## 1. Exact compiled inventory

The artifact has the active family list `Disproportionation`,
`R_Addition_MultipleBond`, `R_Recombination`, and `intra_H_migration`. Its 314
records are:

| Family | Inventory class | Enabled | Irreversible | Refused | Total |
|---|---:|---:|---:|---:|---:|
| Disproportionation | `R1:quinoid_disproportionation` | 0 | 180 | 0 | 180 |
| Disproportionation | null | 0 | 1 | 0 | 1 |
| R_Addition_MultipleBond | null | 0 | 22 | 0 | 22 |
| R_Recombination | `R1:J_ring` | 6 | 0 | 0 | 6 |
| R_Recombination | null | 0 | 1 | 0 | 1 |
| intra_H_migration | null | 0 | 104 | 0 | 104 |
| **Total** |  | **6** | **308** | **0** | **314** |

By inventory class this is 180 quinoid-disproportionation, 6 ring-capture, and
128 null-class records. All 308 non-ring records are marked `irreversible`; none
has `reverse_of`. The six enabled ring records are three forward/reverse pairs.
There are no enabled or irreversible non-ring records that consume zero radicals
and produce radicals; the three zero-radical initiators are all ring
dissociations and require an existing `J_ring` participant. Thus an inventory
with neither radicals nor a ring junction has no firing chemical channel.

Relevant code: the family candidate appears at `compiler.py:33`; the
reversible-status overwrite starts at `compiler.py:1434`; generation and forward
record creation are at `compiler.py:1681` and `compiler.py:1683`; and only
`R1:J_ring` enters reverse construction at `compiler.py:1684`.

## 2. Why H_Abstraction is empty

`H_Abstraction` is a global candidate (`compiler.py:33`), but no declaration
returned by `ps_proxy_set` (beginning at `compiler.py:546`) schedules it. The pristine and
junction-radical proxies schedule no family; the three unimolecular radical
proxies schedule addition and intra-H migration; the radical pairs schedule
disproportionation/recombination; and radical + styrene schedules only addition.
Consequently `discover_family_reactions` intersects the global filter with each
proxy's `family_candidates` at `compiler.py:1617` and never calls RMG for
H-abstraction. The artifact therefore excludes it as “no compatible bounded L=3
proxy declaration in the PS family filter.” This is the named drop step:
**per-proxy family scheduling**, before reaction generation or record creation.
Forced generation proves the family itself is not empty: the declared inputs
produce 47 H-abstractions (14 on end-radical + end-radical, 28 on
junction-radical + end-radical, and 5 on end-radical + styrene; all other
declared inputs produce zero). None pairs a radical with a closed-shell polymer
chain: styrene is the only closed-shell co-reactant. Adding the missing pristine
chain donor produces 44, 14, and 14 reactions with the interior-, end-, and
junction-radical proxies, respectively.

## 3. Why there is no backbone homolysis

The pristine proxy schedules an empty family list. Forcing every candidate
family on that same closed-shell L=3 molecule produces 23 `R_Recombination`
dissociations and zero reactions from each of `Disproportionation`,
`H_Abstraction`, `R_Addition_MultipleBond`, and `intra_H_migration`; no other
candidate initiation family fires. The selected backbone example splits
`CC(CC(CCc1ccccc1)c1ccccc1)c1ccccc1` into
`[CH](CCc1ccccc1)c1ccccc1` and `[CH2]C(C)c1ccccc1`. RMG stores that family
reaction in the association direction with `is_forward=False`, and
`_orient_to_proxy` correctly reverses it to the pristine input direction. The
loss therefore occurs at the same scheduling step as above, not inside RMG.
For non-ring reactions the compile loop emits only the direction presented by
the scheduled proxy (`compiler.py:1681-1683`) and creates no generic partner;
the radical-pair recombination proxy consequently cannot supply homolysis as a
reverse either. The precise diagnosis is the conjunction of an empty pristine
schedule and a ring-only reverse special case—not a failure of RMG to find
backbone cleavage.

## 4. Why C3 and C9 pass

C3 (beginning at `compilerRealTest.py:618`, with its comparison at line 627)
defines `required` as the independent oracle
generated from each proxy's own scheduled family list, then checks only that its
`(site_type, family, template)` keys are a subset of the compiled keys. It uses
`ps_proxy_set(PS_PROXY_UNITS)` with `PS_PROXY_UNITS = 3`; it neither constructs
the design-ruling `(2r + 3)`-mers nor runs the full family universe. It therefore
shares both omissions it is meant to detect. The reproduced faithful oracle uses
`r = 1`, full 5-unit polymer reactants (including the end-radical reactant paired
with styrene), and all five candidate families. It generates 2,195 reactions and
finds **128 missing unique keys**. By `(site type, family)`, those keys are:
`doubly_featured` addition 1, recombination 6, migration 24;
`end_radical` recombination 7, migration 4; `interior_radical` recombination 7,
migration 13; `junction_radical` addition 1, recombination 4, migration 7;
`end_radical+end_radical` H-abstraction 5 and addition 4;
`end_radical+styrene` disproportionation 9 and H-abstraction 3;
`junction_radical+end_radical` disproportionation 5, H-abstraction 11, and
addition 8; and pristine recombination/dissociation 9. The separate radical +
pristine donor probe is not one of the eight declared site inputs, so it exposes
another catalogue-input omission but is not counted among those 128 C3 keys.
C9 (beginning at
`compilerRealTest.py:634`) checks that
the irreversible list and ceiling-pair list are nonempty, recomputes the stored
ceiling from the same two compiled rate tables, and validates ring metadata. It
does not compare a literature value and does not require a reverse for every
reversible reaction. A faithful detailed-balance C9 would flag all 308 unlinked
non-ring records, including recombination without homolysis, as well as the
missing H-transfer chemistry; it would also perform the promised independent
literature comparison.

## 5. Other anomaly: frontier refusal is overwritten

Immediately before `compiler.py:1434` a frontier channel first becomes
`refused`, but the reversible/non-ring branch beginning at that line
overwrites that status with `irreversible`. All eight present production proxies
have `frontier=False`, so the bug is latent in this artifact. A synthetic
frontier copy of the real recombination proxy reproduces status `irreversible`
with reason “condensed-phase reference thermochemistry unavailable,” and
`_is_forward_met_record` admits it because the consumer excludes only literal
`refused` records (`met.py:1250`). Such a frontier channel would enter the
normal SSA/channel total instead of the refused propensity, so its hazard would
not contribute to the design-ruling leak numerator.

## Gas-phase reverse-rate probe

RMG's reverse path obtains every species' gas-phase thermo, evaluates
`get_equilibrium_constant(T, type="Kc")` at `compiler.py:1327`, and uses
`k_reverse = k_forward / Kc`. The table uses RMG's Kc convention. The backbone
Kc is written for radical + radical association, so its reverse is unimolecular
homolysis; the H-abstraction Kc is written for the displayed forward transfer,
so both rates are bimolecular.

| T (K) | backbone Kc association | k_homolysis (s⁻¹) | A=10¹⁵ s⁻¹, Ea=320–290 kJ/mol bracket (s⁻¹) | H-abstraction Kc | H-abstraction k_reverse (m³ mol⁻¹ s⁻¹) |
|---:|---:|---:|---:|---:|---:|
| 600 | 1.16895×10¹⁸ | 1.87598×10⁻¹³ | 1.38698×10⁻¹³–5.67218×10⁻¹¹ | 3.02777×10⁻⁵ | 1.50483 |
| 700 | 1.63104×10¹⁴ | 2.02492×10⁻⁹ | 1.32365×10⁻⁹–2.29275×10⁻⁷ | 1.25830×10⁻⁴ | 1.19997 |
| 800 | 2.17518×10¹¹ | 2.04377×10⁻⁶ | 1.27806×10⁻⁶–1.16229×10⁻⁴ | 3.65091×10⁻⁴ | 1.01581 |

The homolysis result lies inside the simple 290–320 kJ/mol bond-cleavage bracket
at all three temperatures. For reference, the selected forward H-abstraction
rates are 4.55627×10⁻⁵, 1.50993×10⁻⁴, and 3.70864×10⁻⁴ m³ mol⁻¹ s⁻¹ at 600,
700, and 800 K. With the compiler's configured monomer concentration, the
gas-phase propagation/depropagation pair implies a PS ceiling temperature of
**710.2487 K**. This is an internally computed value, not the literature check
that C9 claims to perform.

## Options for reverse and initiation rates

| Option | How rates would be obtained | What still lacks thermo/reference state | Irreversible share among current logical compiled pairs | Validation-data calibration |
|---|---|---|---:|---|
| **(a) Gas-phase RMG Kc for every pair** | Estimate/library gas thermo for both sides; calculate Kc(T); emit and link `k_rev = k_fwd/Kc`, exactly as the ring path already does. Schedule pristine `R_Recombination` or create every generic reverse so backbone homolysis is present; add a radical + closed-chain H-abstraction proxy. | Numerically, none of the 195 bounded species failed RMG thermo. Scientifically, every species still lacks a condensed-PS reference state, so this is explicitly a gas-phase approximation. | **0/311** among the current 308 non-ring singleton pairs plus 3 ring pairs, if every reversible record gets a partner. Newly scheduled channels would also be paired, so they do not change the 0% share. | No campaign fit is introduced. The pinned RMG library/rule training is inherited; current C8 establishes source reproducibility, not absence of overlap with validation data. |
| **(b) Gas thermo + a-priori PS-melt solvent stand-in** | Add species solvation free energies to gas free energies before constructing Kc(T), then derive the reverse. The pinned database supports Abraham gas-to-solvent Gibbs and Mintz enthalpy LSERs at 298 K, plus a CoolProp K-factor temperature model below a molecular solvent's critical temperature. Solutes require `S,B,E,L,A,V`; a stand-in solvent requires the six `*_g` and six `*_h` coefficients plus `name_in_coolprop`. `alpha,beta` exist for an intrinsic H-abstraction correction but are not used by RMG. | All 195 bounded solutes received descriptors, including radicals. What is absent is a PS-melt solvent entry and a defensible 600–800 K polymer reference: the model is infinite-dilution solvation in a molecular solvent, and its temperature path requires a CoolProp solvent below Tc. A numerical stand-in is therefore not the same as condensed-PS thermo. | **0/311** if the owner accepts and fully specifies the surrogate; newly scheduled channels would likewise be paired. Otherwise the unsupported pairs remain unresolved. | The stand-in and descriptors must be fixed a priori; none may be selected or tuned against validation data. Existing database-training overlap remains the same unresolved provenance question as in (a). |
| **(c) Status quo; initiation elsewhere** | Leave compiled non-ring records one-way. Supply a separately sourced elementary initiation event/rate or an a-priori radical seed in the executor/input layer. A moment-model deck homolysis kernel is forbidden by the design ruling and cannot be used. H-abstraction also remains absent unless supplied separately. | All reverse condensed-phase states remain absent; an external initiation source needs its own declared phase and provenance. A finite initial seed alone is not sustained initiation and the radical population can still extinguish. | **308/311 (99.04%)** current logical pairs remain irreversible: 308 non-ring singleton pairs and 3 enabled ring pairs. | The external rate/seed must be fixed independently of validation data. No value was fitted in this probe. |

Option (a) is already numerically executable on the bounded species and gives the
rate table above, but deciding among these options is outside this report.

## Contract corrections and limits

1. “H_Abstraction ... yields zero records” is true of the final artifact but is
   misleading as a generation diagnosis: RMG produces 47 reactions when forced
   on the declared reactants. The family is never scheduled.
2. The compiler does not generally “rely on reverses.” It compiles the direction
   initiated by each proxy, including direct depropagation, and then adds an
   explicit reverse only for ring capture. Pristine homolysis would compile in
   the correct direction if `R_Recombination` were scheduled there.
3. The literal claim about a completed kMC run cannot be rerun on this branch:
   this worktree contains the compiler, state, MET, and volatility modules but no
   SSA executor. The static artifact result is conclusive—zero non-ring
   zero-radical initiation channels—but it is not an executed trajectory.
4. C3 and C9 do not implement their design-ruling definitions, as detailed
   above. In particular, C9 contains no literature comparison.
5. The requested zero-fit contamination ruling cannot be completed without
   inspecting the forbidden validation datasets. No such dataset was opened.
   The present provenance tests reproduce RMG source nodes and training
   contributors; they do not test whether those sources overlap the validation
   corpus.
