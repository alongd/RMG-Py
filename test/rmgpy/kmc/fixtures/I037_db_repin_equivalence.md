# I-037 database repin equivalence

## Inputs

The authenticated old-pin artifact is
`0883a1292e20708a17ee8dc7960cf5648b354e0c24ce3373dee83813cea85c9d`.
Its database provenance is `4a12d36fcdc193ede82c8d1ab5c1653495d445bc`.
The public-pin artifact is
`56add5736c2701ac9262a9b27953cbd1fa2898cf980912f50567dd968c4e4e7f`,
compiled against `cd86d4e1c187a132109e16cd86f624ed9fb217df`.

Both artifacts record compiler-sources SHA-256
`3d9bdaf6d05448b584976831847d00d05cb153681beff629562eb4dd2bb6db77`.
That fingerprint also matches the compiler sources in this worktree.

## Canonical comparison

Each artifact contains 14,998 records. After removing record identity fields
(`canonical_index`, `event_id`, record provenance, reverse links, and derived
junction reverse handles) and all database-SHA provenance fields, their record
multisets are equal. The SHA-256 of the sorted normalized-record multiset is
`0fdc3faef10a0dd9e58e71d4d663cf7e8bafa37830ed6b60d8f2272071d40ad9` for
both artifacts.

The only raw differences are the old versus public database provenance and
event identifiers derived from it, including references in `discovery` and
`ps_ceiling_pairs`. The top-level RMG-Py provenance also differs because the
authenticated old artifact was compiled from a different RMG-Py revision.

## Database inspection

The database diff contains two development-status reference edits and five
Evidence-of-record prose edits in training reactions. It contains no rate,
structure, family-definition, or executable-behavior change.

The retained I-032, I-034, and I-039 fixtures are historical records at the old
pin and are intentionally not rewritten.
