# I-103 — Measurement report: is Polymer thermo per-site or per-chain?

**Verdict, one sentence with the numbers:** the value is **per-chain (extensive)** — the
whole-molecule group-additivity thermo grows **perfectly linearly** with degree of polymerization
at a constant per-repeat-unit increment of **ΔH298 = −41.254 kJ/mol, ΔS298 = +78.827 J/mol·K,
ΔCp(1000 K) = +103.261 J/mol·K**, so the "per-site" documentation is **not implemented**.

This is a measurement only. No production file was modified.

## What varied, and what was held fixed

Chains were built with the class's **own** shared stitching recipe,
`Polymer._capped_chain_species(dp)` — `end_groups[0]--(monomer)×dp--end_groups[1]`, the identical
recipe the proxy uses (`_stitch_trimer` is this recipe frozen at dp=3). The pool was fixed for all
four cases: PE, monomer `[CH2][CH2]`, end caps `['[CH3]', '[H]']`, `Mn=5000, Mw=6000,
initial_mass=1.0`. `Mn`/`Mw`/`initial_mass` are pool-distribution moments the constructor demands;
they never enter `_capped_chain_species`. **The only quantity that varied across the four cases is
`dp`** — confirmed by the formula progression C5H12 → C7H16 → C9H20 → C13H28 (exactly +C2H4 per DP)
and by heavy-atom count +2.00 per DP.

## Backend MRO (environment measured in)

`rmgpy_env`, thermo obtained through:

```
rmgpy.data.thermo.ThermoDatabase
builtins.object
```

## Route to synchronous thermo

`ThermoDatabase.get_thermo_data(species)`, called directly on the constructed `Species`. **Not** the
`thermoengine.submit()` path (async; does not populate `spc.thermo` synchronously — the trap that
stopped the previous attempt). Evidence the data is real GA output, not a default or null: every
returned object is a `ThermoData` whose comment is `Thermo group additivity estimation:
group(Cs-CsCsHH) + ...`, with the count of `Cs-CsCsHH` backbone groups increasing by exactly one
per added repeat unit.

## Table — DP ∈ {2,3,4,6}

| DP | formula | heavy | H298 [kJ/mol] | S298 [J/mol·K] | Cp(1000 K) [J/mol·K] |
|---:|:--------|------:|--------------:|---------------:|---------------------:|
|  2 | C5H12   |     5 |      −147.235 |        348.679 |              278.487 |
|  3 | C7H16   |     7 |      −188.489 |        427.505 |              381.748 |
|  4 | C9H20   |     9 |      −229.743 |        506.332 |              485.009 |
|  6 | C13H28  |    13 |      −312.252 |        663.985 |              691.532 |

## First differences (per added repeat unit)

| DP→DP | ΔH298 [kJ/mol] | ΔS298 [J/mol·K] | ΔCp(1000 K) [J/mol·K] | Δheavy |
|:------|---------------:|----------------:|----------------------:|-------:|
| 2→3   |        −41.254 |          78.827 |               103.261 |   2.00 |
| 3→4   |        −41.254 |          78.827 |               103.261 |   2.00 |
| 4→6   |        −41.254 |          78.827 |               103.261 |   2.00 |

The increments are **identical to three decimals across all three intervals** — the relationship is
linear with no meaningful offset. This is the expected signature of additive group-additivity thermo
over a growing molecule: extensive, not intensive.

## Verdict

**Extensive / per-chain.** The returned thermo scales with molecule size; it is not flat in DP. The
per-DP increment a rescaling to a true per-site basis would have to divide out is
**ΔH298 = −41.254 kJ/mol, ΔS298 = +78.827 J/mol·K, ΔCp(1000 K) = +103.261 J/mol·K per repeat unit.**

### Consequence for the accessor (reported, not a ruling)

`Polymer`'s live proxy is a **fixed trimer** (`_stitch_trimer` = `_capped_chain_species(3)`), so
`Polymer.get_thermo_data` returns exactly the **DP=3 row** above — a per-**chain**, three-site
whole-molecule quantity — with no division by site count. Whether that is a defect or a naming slip
is explicitly **not** this ticket's call; the number is now measured. What the number says: the
value is per-chain, and the "per-site basis" asserted in the four docstrings and the CRITICAL NOTE
is not what the code computes.

## What I verified vs. inferred

- **Verified by running code:** the four thermo triples; the constant linear increments; that the
  data is genuine GA output (comments); that the construction varies only DP (formula/heavy-atom
  progression); the backend MRO; that `get_thermo_data` is synchronous and returned non-null.
- **Verified by reading source:** `self.thermo = proxy.thermo` is a bare assignment;
  `get_thermo_data` has no arithmetic; the four accessors delegate unchanged; the four docstrings and
  CRITICAL NOTE all say "per-site basis"; no division by site count exists; `_stitch_trimer` is the
  dp=3 instance of `_capped_chain_species`.
- **Inferred (not separately executed):** that the fixed-trimer proxy's own thermo equals the DP=3
  row — it follows from both using the identical stitching recipe, but I measured via
  `_capped_chain_species(3)`, not by extracting `_stitch_trimer`'s Species and re-running.
- **Could not hold fixed by construction:** nothing bearing on the result. Formula and heavy-atom
  count necessarily grow with DP — that IS the independent variable — and no other structural
  feature (branching, aromatics, end-cap identity) changed.
