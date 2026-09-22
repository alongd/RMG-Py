# Round 102 — one open for the signature and the content, one enumeration for the fields

Base `i221-saturation-no-atomtype@10ef0059a`. Probe:
`docs/i221-saturation/probes/round102_one_open_probe.py`.

| check | result | commit |
|---|---|---|
| round-102 probe | **3 of 3 findings reproduce**, 4 of 4 controls hold | `10ef0059a` |
| round-102 probe | **0 of 3**, 4 of 4 controls hold, 2 noted | this branch |
| `quarantineTest.py` | **5 failed**, 124 passed, 7 skipped | `10ef0059a` |
| `quarantineTest.py` | **136 passed** | this branch |
| `test/rmgpy/data` | 504 passed, 7 skipped | this branch |
| `test/rmgpy/rmg` | 149 passed, 2 skipped, 1 failed + 1 error (`ps2pdf`) | this branch |
| RMG-database `test/` | 320 passed, 4 errors — unchanged since round 95 | this branch |

Both diagnoses are correct and both are repaired. **One claim in the brief does not reproduce**, and
it is the one the brief was most specific about; §1 sets it out.

---

## 1. HIGH 1 — the defect is real, the exploitation is not, and `ctime` is why

**Confirmed.** `_manifest_signature` stats the path, `load_family_quarantine` then opens it
independently, and a rename between them stores B's quarantine under A's signature. Measured at
`10ef0059a` with the rename driven from inside the stat:

```
first -> reason = 'manifest B'
the cache holds B's quarantine under:
    cached identity = (A's mtime_ns, ctime_ns, size, ino)
    B's identity    = (B's ...)
a key and a value that describe two different files
```

**But the consequence the brief specifies — "restoring A produces a cache *hit* returning B" — does
not happen at this base, and the reason is round 99's own repair.** `os.rename` moves the inode's
`st_ctime_ns`; round 99 put `st_ctime_ns` into the identity; so the restored A no longer matches the
stale key and the lookup **misses** and re-reads A. Measured, not argued — `resolve_quarantine`
returns `'manifest A'`, and there is a test pinning exactly that
(`test_restoring_the_original_is_a_cache_miss`).

So: the brief's *"`st_ctime_ns` does not help: it is computed from the first open"* is the one thing
I disagree with. It is computed from the first open, and that is precisely what makes the stale key
fail to match the restored file. On an ordinary POSIX filesystem this is not exploitable by rename.

**Repaired anyway, and the reasoning is not "defence in depth".** A cache whose key describes a file
that was never parsed is wrong on its own terms, the failure mode is silent, and the residual
depends on the *same* filesystem property as round 99's declined-digest MEDIUM — coarse or synthetic
`ctime`. Where that property fails, this becomes exploitable exactly as described. The fix costs
nothing: the `fstat` was already being taken to confirm the file is regular, so the identity is now
read off **that** descriptor and `resolve_quarantine` re-keys the cache entry with it. The check and
the thing used are one object.

This is the change I named at the end of round 101 as the next thing worth doing in this area.

## 2. HIGH 2 — confirmed, including the test critique, with one thing the brief could not have known

**Confirmed in both files.** The template shape forwarded `degeneracy` and `electrons` and not
`elementary_high_p`, `allow_pdep_route` or `allow_max_rate_violation`; **all three** shapes dropped
`allow_max_rate_violation`. And the round-95 test did assign those three onto the loader's output,
so the converter was only ever asked to carry what the fixture had put there — the
fixture-assigns-what-production-never-sets shape, exactly as named.

**What the brief could not have known:** *"fix the loader to forward every field"* cannot be done by
adding constructor arguments. `TemplateReaction.__init__` has **no** `elementary_high_p`,
`allow_pdep_route` or `allow_max_rate_violation` parameter — they are base-`Reaction` attributes,
not constructor parameters of that subclass. Adding them would mean editing
`rmgpy/data/kinetics/family.py`, which is **outside this round's gates**. So the fields are carried
onto the constructed object by `_carry_entry_fields()` immediately after construction, in the
loader, with the reason in the code.

**The enumeration, and how I know it is complete.** Not a list — a partition, derived from the class:

```
_REACTION_FIELDS = every public attribute of Reaction presenting as a getset_descriptor
                   (what a Cython `cdef public` looks like, which separates state from
                    methods and from the cimported types that also appear in dir())
```

23 fields. **14 carried** from `entry.item`:

`allow_max_rate_violation, allow_pdep_route, degeneracy, duplicate, electrons,
elementary_high_p, is_forward, network_kinetics, pairs, protons, rank, reversible,
specific_collider, transition_state`

**9 excluded, each with a reason in `_NOT_CARRIED_FROM_ENTRY`:** `reactants`/`products` (the
constructor copies by slice so a model edit cannot write back into the database entry), `kinetics`
(comes from `entry.data`), `index` (a library position is not a model index), `label`, `comment`,
`SurfaceArrhenius`/`SurfaceChargeTransfer` (cimported types, not state), `k_effective_cache` (a
memo). **0 unclassified.**

A field added to `Reaction` tomorrow is carried by default; it can only be left behind by someone
writing down why. `test_every_field_of_reaction_is_carried_or_excluded_with_a_reason` pins the
partition, refuses a stale exclusion naming a field that no longer exists, and refuses a discovery
that has silently stopped finding fields. It **caught one of my own exclusions with a two-word
reason** while I was writing it.

The tests now take the flags through the loader: `_library_declaring` puts them on `entry.item` —
where `KineticsLibrary.load_entry` puts them when a library file declares them — and nothing writes
to the loader's output. `test_the_library_file_format_can_declare_all_three` drives the real
`load_entry` to show that is where they arrive from.

## 3. The two censuses

**(a) Names resolved more than once where one resolution would do — 5 sites in `quarantine.py`.**

| site | status |
|---|---|
| `_manifest_signature`'s `lstat` vs `_read_manifest`'s descent | **was the defect.** Closed on the value side: the cache is keyed on the `fstat` identity of the file parsed. The `lstat` remains for the cache *lookup*, which is irreducible — checking a cache without reading the file requires a stat |
| `_manifest_signature`'s own `lstat` vs its `os.path.isdir` | **residual.** Two resolutions in one function, of the manifest and of its directory. It only classifies; it never reads or executes |
| `_family_directory`'s `islink` + two `realpath` | **by design, documented since round 101.** Diagnostic, so the message can say *"it contains a path separator"* rather than `ELOOP`; the structural closure is the descent |
| `_database_has_any_quarantine`'s `stat` + `listdir` + `lexists` per name | **out of scope, named.** N+2 resolutions, but it decides only whether to warn |
| the descent itself | not a duplicate — one resolution per component, each pinned |

**(b) Reactions constructed from `entry.item` — 5 sites.**

| site | status |
|---|---|
| `library.py:354` (LibraryReaction, "Originally from…") | **fixed** |
| `library.py:373` (TemplateReaction, rate rule) | **fixed** — the one that dropped three |
| `library.py:389` (LibraryReaction, ordinary/pdep) | **fixed** |
| `database.py:493` `generate_reactions_from_library` | **out of gates, named.** Names 8 fields and drops the rest — including **`electrons`**, the charge-balance field round 95 was about, plus all three flags, `network_kinetics`, `pairs`, `transition_state`, `index`. This is the most consequential unfixed instance in either census |
| `database.py:607` (group/molecule entries, family side) | **out of gates, named.** Builds template-side reactions, a different shape |

Both counts are in the probe as a NOTE, so they are re-measured on every run rather than asserted
once in prose.

## 4. What I could not reach

* **No real concurrent attacker.** The rename is driven from inside the stat. The argument that the
  repair closes it is structural — one descriptor — not statistical.
* **`ctime` on anything but local ext4.** §1's whole residual is a filesystem property, and it is
  the same one that bounds the digest MEDIUM declined in round 101. Nothing here measures NFS,
  overlayfs, or a container layer, and on a filesystem with synthetic `ctime` the exploitation in
  the brief becomes live again.
* **Whether any shipped library entry sets `elementary_high_p`.** The defect is that the flag was
  dropped; I did not sweep the databases for entries that declare it, so I cannot say whether any
  real mechanism has been mis-routed by this.
* **`load_entry` cannot declare `electrons`.** The library file format has no such keyword, so a
  plasma library cannot state a signed electron count through the ordinary parser at all. In gates,
  but outside the brief — named rather than changed, because it is a format change and belongs to
  its own round.
* **`database.py`'s two sites** — out of gates, so measured and named, not repaired.

## 5. Named, not fixed

* `database.py:493` dropping `electrons` (above). Recommended as the next round's HIGH.
* The four `test_argon_metastable_thermo.py` errors, unchanged at 320 passed / 4 errors since
  round 95.
* `Plasma_Electron_Impact_Ionization/quarantine.py` still declares `requiresEngineCallSites`.
* `test_make_profile_graph` — `ps2pdf` absent, identical at base.
* Round 101 reported the data suite at 495 passed / **8** skipped and flagged the skip count as
  unexplained. It is 504 / **7** here. The 8th skip was nondeterministic and not caused by that
  round's change; recorded so the discrepancy does not go on being carried forward.
