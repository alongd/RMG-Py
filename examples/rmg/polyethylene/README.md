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

| Channel | RMG family and small-molecule surrogate | A (SI) | n | Ea (J/mol) | Maximum pointwise fit error |
|---|---|---:|---:|---:|---:|
| Backbone homolysis / initiation | Thermodynamic reverse of `R_Recombination`: n-hexane central bond -> 2 n-propyl; one central bond maps to one breakable PE backbone bond | 5.425091046e26 s^-1 | -2.97626501 | 390905.266 | 1.775% |
| Chain-end beta-scission / depropagation | Thermodynamic reverse of `R_Addition_MultipleBond`: 1-hexyl -> ethylene + 1-butyl; exact rule from training reaction 2903 | 4.087800877e9 s^-1 | 1.09830553 | 126440.197 | 1.131% |
| Radical termination | Sum of `R_Recombination` (2 1-hexyl -> n-dodecane) and `Disproportionation` (2 1-hexyl -> n-hexane + 1-hexene) | 1.340078977e10 m3 mol^-1 s^-1 | -1.15292140 | 17212.164 | 0.182% |
| H transfer | `H_Abstraction` from n-octane interior sites plus all six `intra_H_migration` paths of 1-octyl; converted to the solver's pseudo-first-order convention | 2.347940229e2 s^-1 | 2.61167734 | 40786.396 | 1.026% |

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

The deck was reproduced from RMG-Py `0693d2c3a` with the generic singlet-
carbene proxy-bridge fix in this branch and RMG-database `cd86d4e1c`. The
tee-captured clean-worktree run is in `/home/alon/runs/i057-mom-pe/run5/`.

- completion marker: `MODEL GENERATION COMPLETED` (exit 0)
- wall time and peak RSS: 40.70 s and 855688 kB
- final model: 8 core species / 0 core reactions; 22 edge species / 22 edge
  reactions
- ethylene: core species index 7; final amount 1.94100127334012e-5 mol
- initial moments `(mu0, mu1, mu2)`: `(0.01, 1.78233008978049,
  381.204065872431)`
- final moments `(mu0, mu1, mu2)`: `(0.00999999999997937,
  1.7823105812549, 381.197131292584)`
- `mu1` loss: 1.95085255840777e-5 mol of repeat units
- condensed PE mass: 50.0 g initially and 49.9994527241139 g finally;
  released ethylene mass is 5.44512288848575e-4 g, leaving a numerical
  closure residual of -2.76359725035532e-6 g (5.53e-8 of the initial PE
  charge)

The polymer test selection on this base produced 1510 passes, 1 expected
failure, and the three pre-existing failures identified by the I-057 manager;
there were no failures attributable to this example or its generic fix.

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
