# Round 101 — one anchor, and everything below it pinned

Base `i221-saturation-no-atomtype@d57121cb3`. Probe:
`docs/i221-saturation/probes/round101_descent_probe.py`.

| check | result | commit |
|---|---|---|
| round-101 probe | **6 of 6 findings reproduce**, 5 of 5 controls hold | `d57121cb3` |
| round-101 probe | **0 of 6**, 5 of 5 controls hold | this branch |
| `quarantineTest.py` | **10 failed**, 111 passed, 7 skipped | `d57121cb3` |
| `quarantineTest.py` | **128 passed** | this branch |
| `test/rmgpy/data` | 495 passed, 8 skipped | this branch |
| `test/rmgpy/rmg` | 149 passed, 2 skipped, 1 failed + 1 error (`ps2pdf`) | this branch |
| RMG-database `test/` | 320 passed, 4 errors — **unchanged since round 95** | this branch |

All three HIGH confirmed by reading before any probe was written, and all three reproduce. The
diagnosis is right and the framing is right: `O_NOFOLLOW` on one component is not containment of a
path.

---

## The shape of the repair

There is always **exactly one anchor that must be followed**, because a path cannot be descended
without something to descend from, and that first open resolves whatever it is given. The question
is only *which* directory gets to be it. Round 99 made it the family's parent —
`<database>/kinetics/families` — which is reachable from a shipped `family:` label. It is now
`settings['database.directory']`, which is local configuration, and every component below it
(`kinetics`, `families`, the label, the manifest) is opened with `O_NOFOLLOW` through the previous
descriptor.

That makes the rule one sentence: **nothing below the anchor may be a link.** The new helper is
`_anchor_and_components()`; `_read_manifest` descends what it returns.

**One correction to the brief.** "Either open the families root once with `O_NOFOLLOW`" is not
available: the families root is *below* the database root, so `O_NOFOLLOW` on it would refuse a
database whose own internal layout is linked, and more importantly it leaves `kinetics` unpinned —
which is HIGH 1 one level up, reproduced here as a separate finding. The second half of your
sentence, "resolve the root once and keep the descriptor", is what is implemented, with the root
being the *database* root rather than the families root.

## HIGH 1 — a linked `kinetics/families`, and why no check could catch it

```
resolve_quarantine('A_Family_Reached_Through_A_Linked_Parent') -> (KineticsQuarantine, True)
reason = 'EXECUTED_FROM_OUTSIDE_VIA_FAMILIES'
```

No race needed: the link is simply there. The part worth keeping is *why* `_family_directory`
approves. It compares `realpath(family_path)` against `realpath(families_root)` — and when
`families` is itself the link, **both resolve into the same foreign tree**, so the prefix test
passes and reports containment in a directory that is not the database. A containment check that
resolves its own root through the thing being attacked cannot fail. That is the argument for making
containment structural rather than comparative, and it is the second time this round-family has
produced it.

Reproduced one level up at `kinetics` as a separate finding, deliberately: a repair that names the
level the review happened to name is not a repair.

## HIGH 2 — refuse, and the empty path was worse than described

**Option taken: refuse when containment cannot be guaranteed.** Both sub-cases:

* **No `dir_fd`.** The fallback opened the joined pathname and followed every link in it.
  `_read_manifest` now returns a refusal with the reason *"this platform provides no directory
  descriptors, so no component of the path can be opened without following links"*. A hypothetical
  platform loses the quarantine gate and is told so, which is the failure mode this module exists to
  prefer. On every platform this ships on, `os.open in os.supports_dir_fd` is true, so the branch is
  unreachable in practice — it is a refusal, not a regression.
* **Empty `name`.** Refused, and it was worse than "reachable". `os.path.split('')` is `('', '')`
  and `os.path.join('', 'quarantine.py')` is a **relative** name, so `load_family_quarantine(label,
  '')` opened `quarantine.py` **out of the process's working directory and executed it**. Measured:

  ```
  load_family_quarantine('A_Family', '') -> KineticsQuarantine ...
  reason = 'EXECUTED_FROM_THE_WORKING_DIRECTORY'
  ```

  `os.path.relpath('')` raises `ValueError`, so the first version of the refusal *raised out of the
  gate* — caught by the test, and guarded: a gate that raises on the path it is refusing has failed
  in the direction this module exists to avoid.

## HIGH 3 — the loaded branch, which round 99 did not test

Confirmed exactly as stated, including the mechanism. With the permission error on
`kinetics/families`, `os.path.isdir` cannot look and returns False, so the family reads as "not
where this database would put it", and a **loaded** family answered from its own attribute with
`answered=True`. For a family whose attribute is `None` that is a clean bill of health produced by
a permission error.

```
resolve_quarantine('A_Loaded_Family_Behind_An_Unreadable_Parent') -> (None, True)
```

`kind` is now consulted **first**, before existence. Both loaded cases are pinned — attribute
`None` and attribute set — because refusing the second is the non-obvious half: the object's answer
is not *wrong* there, and it is still refused, since answering from it is the staleness round 95
removed. A control holds the other direction: a loaded family whose attribute is `None` still picks
up a manifest that is on disk and readable.

Round 99's test could not see any of this — it registered no loaded family, and it made the
*family* directory unreadable rather than its parent, which leaves `isdir` returning True.

## LOW — presence, not truthiness

`manifest.get(...) or ()` collapsed "declared empty" into "not declared", so
`requiresEngineCallsInSource = ()` beside a populated `requiresEngineCallSites` bypassed the
ambiguity refusal — exactly the manifest someone writes halfway through the rename. Now keyed on a
module-level `_ABSENT` sentinel. Parametrised over all three combinations of empty/populated.

One addition beyond the brief, same defect: a field **declared empty on its own** now raises too
(*"declares a call-in-source requirement that names no module, so it pins nothing while reading as a
pin"*). That is this campaign's oldest finding — a pin nothing reads is a comment — and the
truthiness collapse was hiding it.

## MEDIUM — the content digest, declined, with the exposure bounded

**Not taken.** The cache exists to avoid reading and re-executing the manifest on every lookup, and
`resolve_quarantine` is called once per authorship question per reaction admission. A digest
requires reading the file to compute it, so adding one **converts an O(1) metadata check into a
file read on every admission** and removes the cache's entire reason for existing. That is a real
cost against a residual I can bound.

**The exposure, precisely.** A stale answer requires *all* of: an equal-length edit, in place (so
the inode is unchanged), with the mtime restored, on a filesystem whose `ctime` is coarse,
synthetic, or attribute-cached. On ordinary POSIX filesystems `ctime` moves on any write and no
userspace call can set it back, so the window is closed; the residual is entirely the filesystem
property. The blast radius is one family's quarantine state until the next signature change or
process restart, and `touch` on the manifest ends it.

**What I would do instead, and did not do because it is outside this round's verifier:** have
`_read_manifest` return the identity from `os.fstat` **on the descriptor it actually read**, and
cache under that, instead of the separate `os.lstat` in `_manifest_signature`. Today the signature
and the read are two independent resolutions of the same file — the round-99 defect, still present
in miniature between those two calls. That is a strictly better use of the same effort than a
digest, and it is the thing I would put in the next round if you want this area pushed further.

## What I could not reach

* **No real concurrent attacker.** Every escape here is static — the link is simply present — which
  is stronger evidence than a race, but it means I have not shown timing behaviour under load.
* **`ctime` on anything but local ext4.** The MEDIUM's residual is exactly a filesystem property
  and nothing here measures it on NFS, overlayfs, or a container layer.
* **Whether the loaded-family branch is reachable in a real run** with an unreadable
  `kinetics/families`. The repair is correct regardless; I have no measurement of a deployment in
  that state.
* **`_database_has_any_quarantine`** (`os.listdir` + `lexists`) still addresses manifests by name.
  It only decides whether to warn, so it cannot cause an execution — but it is now the last
  name-addressed path in the module, and it would be a one-line target for the same descent.
* **The data suite's skip count moved 7 → 8.** No skip is in the gated files (all are pre-existing
  `WIP`/PR skips in `familyTest`, `kineticsTest`, `solvationTest`, `thermoTest`, `transportTest`),
  and I did not chase which of the two blank-reason `thermoTest` skips toggled.

## Named, not fixed

* The four `test_argon_metastable_thermo.py` errors — unchanged at 320 passed / 4 errors since
  round 95. Out of gates.
* `Plasma_Electron_Impact_Ionization/quarantine.py` still declares `requiresEngineCallSites`;
  honoured via the alias, warned at load. Out of gates.
* `test_make_profile_graph` — `ps2pdf` absent from this box, identical at base.
* `test_a_root_family_path_is_refused` passes on **both** engines, because no `/quarantine.py`
  exists on this box to be executed. It is a pin, not a demonstration, and is labelled as such.
