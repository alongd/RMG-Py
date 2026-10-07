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
bar. Ethylene is the real gas-phase `monomer_product` routed by the chain-end
depropagation channel.

## RMG-derived channel map

The parameters below are modified-Arrhenius fits to 15 rate coefficients at
50 K intervals from 300 through 1000 K. The source is the read-only public
RMG-database polymer revision `cd86d4e1c`; thermodynamic reverses use thermo
from the libraries named in `input.py`. No experimental product, mass-loss,
TGA, or molecular-weight-distribution data were used.

Each contributing reaction uses the same source priority as model generation:
matching training/depository kinetics before a rate-rule estimate. Every
ArrheniusEP or ArrheniusBM expression is converted with that reaction's 298 K
enthalpy before rates are evaluated and fitted.

| Channel | RMG family and small-molecule surrogate | A (SI) | n | Ea (J/mol) | Maximum pointwise fit error |
|---|---|---:|---:|---:|---:|
| Backbone homolysis / initiation | Thermodynamic reverse of `R_Recombination`: n-hexane central bond -> 2 n-propyl; per-bond rate multiplied by 2 PE backbone C-C bonds per repeat | 1.085018209e27 s^-1 | -2.97626501 | 373605.266 | 1.775% |
| Chain-end beta-scission / depropagation | Thermodynamic reverse of `R_Addition_MultipleBond`: 1-hexyl -> ethylene + 1-butyl; training reaction 2905 | 3.051203080e9 s^-1 | 1.09630553 | 124921.405 | 1.131% |
| Radical termination | `R_Recombination` training reaction 156 (2 1-hexyl -> n-dodecane) plus enthalpy-converted `Disproportionation` (2 1-hexyl -> n-hexane + 1-hexene) | 3.615275991e6 m3 mol^-1 s^-1 | 0.14991913 | 0.000 | 0.070% |
| H transfer | Enthalpy-converted `H_Abstraction` from n-octane interior sites plus all six `intra_H_migration` paths of 1-octyl; converted to the solver's pseudo-first-order convention | 1.312764315e-3 s^-1 | 4.11498243 | 33194.258 | 2.829% |

The initiation fit first derives a per-backbone-bond coefficient from the
n-hexane central bond. The solver multiplies initiation by `mu1 - mu0`, which
counts repeat-to-repeat links, while a capped PE chain of degree `d` has
`2d - 1` backbone C-C bonds. The deck therefore uses two bonds per ethylene
repeat (the large-chain conversion gives `2d - 2`, leaving one backbone bond
per chain unrepresented by this moment closure). A styrene repeat also
contributes two backbone C-C bonds, so the same factor applies to the PS
approximation; the existing PS deck is intentionally unchanged in this task.

The unconstrained termination fit returned `Ea = -24.884 J/mol`. Because the
QSSA input contract requires nonnegative activation energies, its reported row
is the least-squares boundary fit with `Ea` fixed at zero; the maximum rate-point
error remains 0.070%.

For the intermolecular part of transfer, the three secondary-octyl product
isomers cover six interior CH2 sites, or three ethylene-repeat equivalents.
Their summed bimolecular coefficient is divided by three and premultiplied by
the declared initial repeat concentration, 33864.2776785 mol/m3 (950 kg/m3 /
28.05316 g/mol). The unimolecular `intra_H_migration` sum is then added before
fitting. This matches the `polymer()` contract: initiation and depropagation
are first-order, termination is bimolecular, and transfer is already
pseudo-first-order.

The full rate-point table, source comments, thermochemistry provenance, and
reproduction script for this example run are in
`/home/alon/runs/i057-mom-pe/rate_estimation/`.

## Reproduced run

The corrected deck was reproduced from the `i057-mom-pe` branch (based on
RMG-Py `0693d2c3a`) with RMG-database `cd86d4e1c`. The split-stream,
tee-captured run is in
`/home/alon/runs/i057-mom-pe/rework-1/pe-run2/`.

- completion marker: `MODEL GENERATION COMPLETED` (exit 0)
- wall time and peak RSS: 36.15 s and 854580 kB
- final model: 8 core species / 0 core reactions; 22 edge species / 22 edge
  reactions
- requested termination time: 0.1 s; actual final solver time after crossing:
  0.15338938826752 s
- ethylene: core species index 7; final amount 6.95528733240413e-4 mol
- initial moments `(mu0, mu1, mu2)`: `(0.01, 1.78233008978049,
  381.204065872431)`
- final moments `(mu0, mu1, mu2)`: `(0.00999999999925401,
  1.78163052862199, 380.955444591674)`
- `mu1` loss at 0.15338938826752 s: 6.99561158501671e-4 mol of
  repeat units
- condensed PE mass: 50.0 g initially and 49.9803750954295 g finally;
  released ethylene mass is 0.0195117788381906 g, leaving a numerical
  closure residual of -1.13122290997580e-4 g (2.26e-6 of the initial PE
  charge)

The polymer test selection on this base produced 1511 passes, 1 expected
failure, and exactly the three pre-existing failures identified by the I-057
manager; there were no additional failures attributable to this example or
its generic fix.

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
