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
bar. N2 remains in `initialMoles` as a gas and collider, but only PP is listed
in `polymer_phase.species`; consequently the 50 g PP charge, not the N2 mass,
sets the condensed volume. Propylene is the real gas-phase `monomer_product`
routed by chain-end depropagation.

## RMG-derived channel map

The parameters below are modified-Arrhenius fits to 15 rate coefficients at
50 K intervals from 300 through 1000 K. They use only RMG family kinetics and
RMG thermo from the read-only database snapshot. The snapshot is declared as
database revision `cd86d4e1c`, but contains no `.git` metadata with which to
verify that declaration in place. No experimental product,
mass-loss, TGA, yield, or molecular-weight-distribution data were consulted.

The shared `examples/rmg/polymer_rate_estimation.py` helper sends every
surrogate through `CoreEdgeReactionModel.make_new_reaction`. This is the
production model-generation seam: it generates production thermo, applies the
two-direction choice for own-reverse families, converts kinetics, and applies
`fix_barrier_height` (including the pressure-dependence positive-barrier pass).
The helper records both requested and final orientations and evaluates the
requested physical direction without adding degeneracy a second time.

| Channel | RMG family and small-molecule surrogate | A (SI) | n | Ea (J/mol) | Maximum pointwise fit error |
|---|---|---:|---:|---:|---:|
| Backbone homolysis / initiation | Thermodynamic reverse of `R_Recombination`: 2,4-dimethylpentane -> isobutyl + isopropyl; rate-rule entry 111, degeneracy 1; per-bond rate multiplied by 2 PP backbone C-C bonds per repeat | 9.026483352e26 s^-1 | -2.63612559 | 367560.033 | 1.939% |
| Chain-end beta-scission / depropagation | Thermodynamic reverse of `R_Addition_MultipleBond`: 4-methyl-2-pentyl -> propylene + isopropyl; training reaction 239, degeneracy 1 | 8.798203848e11 s^-1 | 0.47152570 | 105603.715 | 0.673% |
| Radical termination | `R_Recombination` of two 4-methyl-2-pentyl radicals (entry 160, degeneracy 0.5) plus both `Disproportionation` products (entries 226 and 20, degeneracies 2 and 3) | 1.201553179e8 m3 mol^-1 s^-1 | -0.35029901 | 0.000 | 0.511% |
| H transfer | Tertiary-site `H_Abstraction` from 2,4,6-trimethylheptane plus the tertiary `R5H_CCC` 1,5-shift of a PP-trimer-length end radical; converted to the solver's pseudo-first-order convention | 5.095727471e-3 s^-1 | 3.79856761 | 34881.095 | 7.756% |

Production-path re-estimation changed the checked-in rows as follows (A, n,
Ea in the units above):

| Channel | Previous row | Production-path row |
|---|---|---|
| Initiation | 8.720264978e26, -2.62997442, 367563.716 | 9.026483352e26, -2.63612559, 367560.033 |
| Depropagation | 1.143987264e12, 0.43823008, 105795.438 | 8.798203848e11, 0.47152570, 105603.715 |
| Termination | 1.201553179e8, -0.35029901, 0.000 | unchanged |
| Transfer | 5.095722545e-3, 3.79856774, 34881.095 | 5.095727471e-3, 3.79856761, 34881.095 |

The initiation fit first derives a single secondary--tertiary backbone-bond
coefficient. The solver multiplies initiation by `mu1 - mu0`, then the deck's
factor of two gives `2(mu1-mu0)` sites. A finite capped PP chain has
`2mu1-mu0` backbone bonds, so the moment approximation omits exactly one bond
per chain. Initially the deficit is `mu0/(2mu1-mu0) = 0.4226%`; it grows if
the mean degree of polymerization falls.

The unconstrained termination fit returned `Ea = -174.324 J/mol`. Because the
QSSA input contract requires nonnegative activation energies, the reported row
is the least-squares boundary fit with `Ea` fixed at zero.

For intermolecular transfer, the two tertiary-abstraction product rows carry
reaction-path degeneracies 1 and 2, spanning three tertiary C-H sites. Their
summed bimolecular coefficient is divided by three repeat equivalents and
premultiplied by the declared initial repeat concentration, 21387.9610774
mol/m3 (900 kg/m3 / the exact RMG repeat mass 42.0797474217 g/mol). The
unimolecular tertiary 1,5-shift is then added before fitting. Thus transfer is
already pseudo-first-order when it enters `polymer()` and remains fixed at the
initial repeat concentration as `mu1` decreases.

A declaration probe gives `V_poly = 5.55555555556e-5 m3`, N2 classified as
gas, and `mu1/V_poly = 21387.9610774 mol/m3`. With N2 incorrectly included in
the phase list the volume was `8.637e-5 m3`, so the old premultiplier did not
describe the simulated phase. The generic fix is to keep N2 in
`initialMoles` but out of `polymer_phase.species`; no solver change is needed.
The PS deck was not changed here. The analogous corrected PS premultiplier
would be `10081.7203363 mol/m3`.

RMG evaluates fitted Arrhenius rows with `R = 8.314472 J mol-1 K-1`, while the
QSSA law is pinned to `8.314`. The checked-in CSV therefore reports both fits;
the maximum solver-law errors are 2.756%, 0.836%, 0.511%, and 7.829% for
initiation, depropagation, termination, and transfer, respectively.

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

The corrected deck completed from branch `i073-mom-pp` with the database input
snapshot declared as RMG-database `cd86d4e1c`.

- completion marker: `MODEL GENERATION COMPLETED` (exit 0)
- wall time and peak RSS: 49.28 s and 854880 kB
- final model: 9 core species / 1 core reaction; 30 edge species / 30 edge
  reactions
- condensed volume: 5.55555555555556e-5 m3; N2 is gas (`True`)
- requested termination time: 0.1 s; actual retained solver time after
  crossing: 0.101260273647456 s
- propylene: core species index 7; final amount 0.474941710101248 mol (the
  core also contains 2.37017802239565e-12 mol of its `C[C]C` isomer)
- initial moments `(mu0, mu1, mu2)`: `(0.01, 1.18822005985366,
  169.424029276636)`
- final moments `(mu0, mu1, mu2)`: `(0.00992450145021161,
  0.711548409580510, 79.2860171982533)`
- `mu1` loss: 0.476671650273150 mol of repeat units
- condensed PP mass: 50.0 g initially and 29.9417773534367 g finally
- propylene-formula core products account for 19.9854272011748 g; the
  remaining condensed-plus-core-gas deficit is 0.0727954453884649 g. The
  pressure-network leak path removes 0.00172994016953104 mol of repeat-unit
  equivalent, accounting for 0.0727954453884620 g. The remaining integration
  error is 2.96984659087229e-15 g.

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
