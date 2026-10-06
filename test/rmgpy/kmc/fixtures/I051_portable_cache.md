# I-051 portable generation cache

The portable identity is the SHA-256 of a canonical JSON object containing:

* the SHA-256 of tracked `rmgpy/` files except `rmgpy/kmc/`;
* a SHA-256 of canonical materialized `input/` content (the same for a git
  checkout and an archive copy);
* Python and RDKit versions; and
* `PYTHONHASHSEED`, retained because the serialized generated-reaction payload
  includes order-sensitive atom identifiers.

The per-request reaction key remains content-addressed by reactants, products,
families, resonance, and generator source. Cache directories are under
`portable/<sha256(identity)>/`, independent of repository history.

Commands, from a checkout, are:

```bash
python -m portable_cache identity CACHE DATABASE
python -m portable_cache migrate CACHE DATABASE
python -m portable_cache export CACHE DATABASE ARCHIVE
python -m portable_cache import DESTINATION DATABASE ARCHIVE
```

Migration copies only entries from old `<repository>-<database>/<hash-seed>`
directories whose database identity, origin commit's non-kMC `rmgpy` tree,
and hash seed all match. It reports and skips the rest, and never removes or
moves anything. Export writes one gzip tar archive
containing generated reactions, compiled artifacts, and `manifest.json`.
Import refuses the archive unless its manifest exactly matches the current
identity.

Compiled artifacts use a sibling identity that adds the tracked
`rmgpy/kmc/` tree hash to the generation identity. Compiler changes therefore
invalidate artifacts without invalidating generated-reaction entries.
The artifact key also contains the effective database SHA, family universe,
temperature grid, and PLP-SEC selection; lookup uses the manifest's exact
artifact filename rather than an arbitrary JSON file.
