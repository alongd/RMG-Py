# I-059 compiled-event semantic equivalence

## Channel identity

A channel is a directed, rooted, atom-mapped rewrite.  The checker constructs one
colored graph containing every explicit reactant and product atom and bond,
participant-membership nodes, the heavy-atom map, and the record's rooted
`bond_ops`.  Reactant and product sides have different colors, so direction is
part of the identity.  Participant positions, event ids, site-type/family names,
templates, canonical indices, provenance, and record order are not.

Candidate identities are indexed with a 256-bit Weisfeiler-Lehman fingerprint and
are then checked by exact colored-graph isomorphism.  Exact isomorphism makes the
identity invariant to participant permutations and graph automorphisms; the
fingerprint is only an acceleration index and cannot make two non-isomorphic
rewrites equal.

## Compared channel data

All records with the same rewrite are one physical channel.  The checker compares:

- the sum of tabulated rate coefficients at every temperature;
- the sum of `k * ssa_multiplier`, the stochastic physical propensity
  coefficient, at every temperature;
- summed raw and effective degeneracy and summed symmetry-adjusted degeneracy;
- the SSA multiplier set and reactant-pair convention;
- rate order, rate units, and rate availability;
- enabled/refused/irreversible status and reciprocal `reverse_of` linkage, with
  links resolved to chemistry rather than event ids; and
- which product graphs remain in the melt and which leave as coproducts.

This aggregation makes a correctly partitioned split channel equal to its
unsplit form while detecting an over-counted or under-counted split.

## Tolerance

The command-line `--rtol` and Python `rtol` default to `1e-9`.  It applies to
rates, propensities, and degeneracy accounting with zero absolute tolerance.
Compiled values are deterministic double-precision tables; `1e-9` admits only
roundoff or harmless summation-order differences and is far smaller than any
kinetic-model uncertainty.  Temperature grids, discrete fields, graph identity,
and linkage are exact.

## Known limits

- Comparison is limited to serialized information.  It cannot detect chemistry,
  stereochemistry, phase state, or provenance that the event record omitted.
- Rates are checked only on the stored temperature grid.  The checker cannot
  assess behavior between grid points or outside the table.
- `rewrite_record/0.1` identifies coproducts by formula rather than by a product
  graph handle.  If non-isomorphic products share one formula and only some of
  them leave, fate is inherently ambiguous; the checker follows the current
  writer convention that the first product stays and later products leave.
- Heavy atoms have an explicit atom map.  Hydrogen correspondence is visible
  only through the rooted rewrite operations and the before/after graphs.
- Site-type, family, template, event-id, provenance, and record-order differences
  are intentionally invisible because they are not chemical channel identity.
