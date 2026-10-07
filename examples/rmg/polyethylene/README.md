# Polyethylene method-of-moments example

This is the second linear-homopolymer example for the method-of-moments
polymer reactor. It uses the same generic `polymer()` declaration and radical
QSSA channel as the polystyrene example; no polyethylene-specific solver code
is required.

## Setup and run

Use RMG-Py `polymer-mom` and the public RMG-database `polymer` revision
`cd86d4e1c`. Build RMG-Py and place the checkout on `PYTHONPATH`. In the run
directory, configure the database through an `rmgrc` file:

```text
database.directory : /absolute/path/to/RMG-database/input
```

Then run, from a separate output directory:

```bash
python /absolute/path/to/RMG-Py/rmg.py /absolute/path/to/RMG-Py/examples/rmg/polyethylene/input.py
```

The deck declares 50 g of H-terminated PE (`Mn = 5000 g/mol`, `Mw = 6000
g/mol`) in a nominal 950 kg/m3 polymer phase, with inert N2 at 1000 K and 1
bar. N2 remains in `initialMoles` as a gas and collider, but only PE is listed
in `polymer_phase.species`; consequently the 50 g PE charge, not the N2 mass,
sets the condensed volume. Ethylene is the real gas-phase `monomer_product`
routed by the chain-end depropagation channel.

## RMG-derived channel map

The parameters below are modified-Arrhenius fits to 15 rate coefficients at
50 K intervals from 300 through 1000 K. The source is the read-only public
RMG-database polymer revision `cd86d4e1c`; the snapshot contains no `.git`
metadata with which to verify that declaration in place. Thermodynamic
reverses use thermo from the libraries named in `input.py`. No experimental product, mass-loss,
TGA, or molecular-weight-distribution data were used.

The shared `examples/rmg/polymer_rate_estimation.py` helper sends every
surrogate through `CoreEdgeReactionModel.make_new_reaction`. This is the
production model-generation seam: it generates production thermo, applies the
two-direction choice for own-reverse families, converts kinetics, and applies
`fix_barrier_height` (including the pressure-dependence positive-barrier pass).
The helper records both requested and final orientations and evaluates the
requested physical direction without adding degeneracy a second time.

| Channel | RMG family and small-molecule surrogate | A (SI) | n | Ea (J/mol) | Maximum pointwise fit error |
|---|---|---:|---:|---:|---:|
| Backbone homolysis / initiation | Thermodynamic reverse of `R_Recombination`: n-hexane central bond -> 2 n-propyl; per-bond rate multiplied by 2 PE backbone C-C bonds per repeat | 9.231130303e26 s^-1 | -2.95577727 | 373487.479 | 1.640% |
| Chain-end beta-scission / depropagation | Thermodynamic reverse of `R_Addition_MultipleBond`: 1-hexyl -> ethylene + 1-butyl; training reaction 2905 | 2.611138149e9 s^-1 | 1.11590093 | 124806.109 | 0.909% |
| Radical termination | `R_Recombination` training reaction 156 (2 1-hexyl -> n-dodecane) plus enthalpy-converted `Disproportionation` (2 1-hexyl -> n-hexane + 1-hexene) | 3.615275991e6 m3 mol^-1 s^-1 | 0.14991913 | 0.000 | 0.070% |
| H transfer | Production-processed `H_Abstraction` from n-octane interior sites plus all six `intra_H_migration` paths of 1-octyl; converted to the solver's pseudo-first-order convention | 5.573253174e-9 s^-1 | 5.74859437 | 26820.365 | 3.073% |

Production-path re-estimation changed the checked-in rows as follows (A, n,
Ea in the units above):

| Channel | Previous row | Production-path row |
|---|---|---|
| Initiation | 1.085018209e27, -2.97626501, 373605.266 | 9.231130303e26, -2.95577727, 373487.479 |
| Depropagation | 3.051203080e9, 1.09630553, 124921.405 | 2.611138149e9, 1.11590093, 124806.109 |
| Termination | 3.615275991e6, 0.14991913, 0.000 | unchanged |
| Transfer | 1.312764315e-3, 4.11498243, 33194.258 | 5.573253174e-9, 5.74859437, 26820.365 |

The initiation fit first derives a per-backbone-bond coefficient from the
n-hexane central bond. The solver multiplies initiation by `mu1 - mu0`, which
counts repeat-to-repeat links, while a capped PE chain of degree `d` has
`2d - 1` backbone C-C bonds. The deck's factor of two therefore uses
`2(mu1-mu0)` instead of `2mu1-mu0`, omitting exactly one bond per chain.
Initially the deficit is `mu0/(2mu1-mu0) = 0.2813%`; it grows if the mean
degree of polymerization falls. A styrene repeat also contributes two
backbone C-C bonds, so the same finite-chain approximation applies to PS; the
existing PS deck is intentionally unchanged in this task.

The unconstrained termination fit returned `Ea = -24.884 J/mol`. Because the
QSSA input contract requires nonnegative activation energies, its reported row
is the least-squares boundary fit with `Ea` fixed at zero; the maximum rate-point
error remains 0.070%.

For the intermolecular part of transfer, the three secondary-octyl product
isomers cover six interior CH2 sites, or three ethylene-repeat equivalents.
Their summed bimolecular coefficient is divided by three and premultiplied by
the declared initial repeat concentration, 33864.2717058 mol/m3 (950 kg/m3 /
the exact RMG repeat mass 28.0531649478 g/mol). The unimolecular
`intra_H_migration` sum is then added before fitting. This matches the
`polymer()` contract: initiation and depropagation are first-order,
termination is bimolecular, and transfer is already pseudo-first-order. Its
coefficient remains fixed at the initial repeat concentration as `mu1`
decreases.

A declaration probe gives `V_poly = 5.26315789474e-5 m3`, N2 classified as
gas, and `mu1/V_poly = 33864.2717058 mol/m3`. The generic fix is to keep N2 in
`initialMoles` but out of `polymer_phase.species`; no solver change is needed.
The PS deck was not changed here. Its analogous corrected premultiplier would
be `10081.7203363 mol/m3`.

RMG evaluates fitted Arrhenius rows with `R = 8.314472 J mol-1 K-1`, while the
QSSA law is pinned to `8.314`. The checked-in CSV therefore reports both fits;
the maximum solver-law errors are 2.365%, 1.190%, 0.070%, and 3.010% for
initiation, depropagation, termination, and transfer, respectively.

The full [rate-point table](rate_estimation/rate_points.csv) and
[reproduction script](rate_estimation/estimate_rates.py) are kept with this
example. From the RMG-Py repository root, regenerate all fit artifacts with:

```bash
PYTHONPATH=. python examples/rmg/polyethylene/rate_estimation/estimate_rates.py \
    --database ../RMG-database/input --output-dir /tmp/pe-rate-estimation
```

Omit `--database` to use `database.directory` from `rmgrc`. The checked-in CSV
can then be compared directly with `/tmp/pe-rate-estimation/rate_points.csv`.

## Reproduced run

The corrected deck completed from branch `i073-mom-pp` with the database input
snapshot declared as RMG-database `cd86d4e1c`.

- completion marker: `MODEL GENERATION COMPLETED` (exit 0)
- wall time and peak RSS: 45.59 s and 855116 kB
- final model: 8 core species / 0 core reactions; 22 edge species / 22 edge
  reactions
- condensed volume: 5.26315789473684e-5 m3; N2 is gas (`True`)
- requested termination time: 0.1 s; actual final solver time after crossing:
  0.113253992628224 s
- ethylene: core species index 7; final amount 6.98414356410886e-4 mol
- initial moments `(mu0, mu1, mu2)`: `(0.01, 1.78233008978049,
  381.204065872431)`
- final moments `(mu0, mu1, mu2)`: `(0.00999999999925203,
  1.78162868726585, 380.954790310106)`
- `mu1` loss at 0.113253992628224 s: 7.01402514641902e-4 mol of
  repeat units
- condensed PE mass: 50.0 g initially and 49.9803234395620 g finally
- ethylene core product accounts for 0.0195927331422908 g. The remaining
  condensed-plus-core-gas deficit is 8.38272957452042e-5 g; the two
  pressure-network leak terms remove 2.98815823101585e-6 mol of ethylene
  equivalent, accounting for 8.38272957447481e-5 g. The remaining integration
  error is 4.56178064073276e-16 g.

The prescribed polymer regression selection completed with 1519 passes and 1
expected failure.

## What the solver does not express

The radical-QSSA channel carries one active chain-end radical population. It
does not retain radical position, distinguish primary end radicals from
secondary mid-chain radicals, or create a new live mid-chain radical after H
transfer. Its single transfer coefficient therefore lumps intra- and
intermolecular routes as a net loss of active ends. The legacy summed
termination block also does not retain separate recombination and
disproportionation products or unsaturated-end inventories. Finally, this
channel emits only the configured monomer from beta-scission; it cannot
represent other PE beta-scission fragments, long-chain branching,
crosslinking, or explicit radical-pair/cage dynamics. Those omissions are
model limitations, not zero-rate claims.
