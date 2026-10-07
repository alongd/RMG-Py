# Polypropylene method-of-moments example

This is a generic method-of-moments example for H-terminated polypropylene
(PP). It uses the same `polymer()` declaration and radical-QSSA channel as the
polyethylene and polystyrene examples. A direct declaration probe created the
PP proxy, moment species, and propylene emission target without any solver or
input-API changes; no PP-specific code path is required.

## Setup and run

Use RMG-Py `polymer-mom` at `7628cb1fb` and the public RMG-database `polymer`
revision `cd86d4e1c`. Build RMG-Py and place the checkout on `PYTHONPATH`. In a
separate run directory, configure the database through `rmgrc`:

```text
database.directory : /absolute/path/to/RMG-database/input
```

Then run:

```bash
python /absolute/path/to/RMG-Py/rmg.py \
    /absolute/path/to/RMG-Py/examples/rmg/polypropylene/input.py
```

The deck declares 50 g of H-terminated PP (`Mn = 5000 g/mol`, `Mw = 6000
g/mol`) in a nominal 900 kg/m3 polymer phase, with inert N2 at 1000 K and 1
bar. Propylene is the real gas-phase `monomer_product` routed by chain-end
depropagation.

## RMG-derived channel map

The parameters below are modified-Arrhenius fits to 15 rate coefficients at
50 K intervals from 300 through 1000 K. They use only RMG family kinetics and
RMG thermo from the read-only database snapshot. No experimental product,
mass-loss, TGA, yield, or molecular-weight-distribution data were consulted.

Each reaction uses model generation's source priority: matching
training/depository kinetics before a rate-rule estimate. Every ArrheniusEP or
ArrheniusBM expression is converted with that reaction's 298 K enthalpy before
rates are evaluated and fitted.

| Channel | RMG family and small-molecule surrogate | A (SI) | n | Ea (J/mol) | Maximum pointwise fit error |
|---|---|---:|---:|---:|---:|
| Backbone homolysis / initiation | Thermodynamic reverse of `R_Recombination`: 2,4-dimethylpentane -> isobutyl + isopropyl; rate-rule entry 111, degeneracy 1; per-bond rate multiplied by 2 PP backbone C-C bonds per repeat | 8.720264978e26 s^-1 | -2.62997442 | 367563.716 | 2.027% |
| Chain-end beta-scission / depropagation | Thermodynamic reverse of `R_Addition_MultipleBond`: 4-methyl-2-pentyl -> propylene + isopropyl; training reaction 239, degeneracy 1 | 1.143987264e12 s^-1 | 0.43823008 | 105795.438 | 1.015% |
| Radical termination | `R_Recombination` of two 4-methyl-2-pentyl radicals (entry 160, degeneracy 0.5) plus both `Disproportionation` products (entries 226 and 20, degeneracies 2 and 3) | 1.201553179e8 m3 mol^-1 s^-1 | -0.35029901 | 0.000 | 0.511% |
| H transfer | Tertiary-site `H_Abstraction` from 2,4,6-trimethylheptane plus the tertiary `R5H_CCC` 1,5-shift of a PP-trimer-length end radical; converted to the solver's pseudo-first-order convention | 5.095722545e-3 s^-1 | 3.79856774 | 34881.095 | 7.756% |

The initiation fit first derives a single secondary--tertiary backbone-bond
coefficient. The solver multiplies initiation by `mu1 - mu0`, which counts
repeat-to-repeat links, while a long PP chain contributes two backbone C-C
bonds per propylene repeat. The deck therefore multiplies the per-bond
coefficient by two.

The unconstrained termination fit returned `Ea = -174.324 J/mol`. Because the
QSSA input contract requires nonnegative activation energies, the reported row
is the least-squares boundary fit with `Ea` fixed at zero.

For intermolecular transfer, the two tertiary-abstraction product rows carry
reaction-path degeneracies 1 and 2, spanning three tertiary C-H sites. Their
summed bimolecular coefficient is divided by three repeat equivalents and
premultiplied by the declared initial repeat concentration, 21387.9648496
mol/m3 (900 kg/m3 / 42.07974 g/mol). The unimolecular tertiary 1,5-shift is
then added before fitting. Thus transfer is already pseudo-first-order when it
enters `polymer()`.

The full [rate-point table](rate_estimation/rate_points.csv) and
[reproduction script](rate_estimation/estimate_rates.py) are kept with the
example. From the RMG-Py repository root, regenerate the fit artifacts with:

```bash
PYTHONPATH=. python examples/rmg/polypropylene/rate_estimation/estimate_rates.py \
    --database ../RMG-database/input --output-dir /tmp/pp-rate-estimation
```

Omit `--database` to use `database.directory` from `rmgrc`. Compare the
regenerated CSV directly with the checked-in file.

## Reproduced run

The deck completed from RMG-Py `7628cb1fb` with the database input snapshot
configured for RMG-database `cd86d4e1c`.

- completion marker: `MODEL GENERATION COMPLETED` (exit 0)
- plain-deck wall time and peak RSS: 36.86 s and 855276 kB
- final model: 9 core species / 1 core reaction; 30 edge species / 30 edge
  reactions
- requested termination time: 0.1 s; actual retained solver time after
  crossing: 0.102198246642922 s
- propylene: core species index 7; final amount 0.484674829734725 mol (the
  core also contains 2.41875077515731e-12 mol of its `C[C]C` isomer)
- initial moments `(mu0, mu1, mu2)`: `(0.01, 1.18822005985366,
  169.424029276636)`
- final moments `(mu0, mu1, mu2)`: `(0.00991025651065357,
  0.701759601969518, 77.9008474455832)`
- `mu1` loss: 0.486460457884141 mol of repeat units
- condensed PP mass: 50.0 g initially and 29.5298668016069 g finally
- propylene-formula core products account for 20.3949944169785 g; the
  remaining condensed-plus-core-gas deficit is 0.0751387814146582 g. The
  pressure-network leak path removes 0.00178562814699718 mol of repeat-unit
  equivalent, accounting for 0.0751387814146604 g. The remaining integration
  error is -2.15105711021124e-15 g.

The only styrene-specific text encountered in the generic input path is an
illustrative `monomer_product` example in `rmgpy/rmg/input.py:299-300`, repeated
in its validation messages at `rmgpy/rmg/input.py:421-422` and
`rmgpy/rmg/input.py:446-448`; none is a conditional branch. The generic PP
declaration probe and completed deck found no PE-, PS-, or PP-specific solver
branch, so no generic code change was needed.

The prescribed polymer regression selection completed with 1519 passes and 1
expected failure.

## What the solver does not express

The radical-QSSA channel carries one active chain-end radical population. It
does not retain radical position, distinguish primary end radicals from
secondary or tertiary mid-chain radicals, or explicitly represent PP
tacticity. Its single transfer coefficient lumps intra- and intermolecular
routes as a net loss of active ends rather than creating a new live mid-chain
radical. The summed termination block does not retain separate recombination
and disproportionation products or unsaturated-end inventories. Finally, the
channel emits only configured propylene from beta-scission; it cannot express
other PP beta-scission fragments, methyl-side-group chemistry, long-chain
branching, crosslinking, or explicit radical-pair/cage dynamics. These are
model limitations, not zero-rate claims.
