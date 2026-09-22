# Round 99 — five things the gate took from a name instead of from a fact

Base `i221-saturation-no-atomtype@18ea8278c`. Every measurement below names the commit it was
taken at. Probe: `docs/i221-saturation/probes/round99_manifest_toctou_probe.py`.

| check | result | commit |
|---|---|---|
| round-99 probe | **5 of 5 findings reproduce**, 8 of 8 controls hold | `18ea8278c` |
| round-99 probe | **0 of 5**, 8 of 8 controls hold, 4 noted | this branch |
| `quarantineTest.py` | **9 failed**, 100 passed, 7 skipped | `18ea8278c` |
| `quarantineTest.py` | **116 passed** | this branch |
| `test/rmgpy/data` | 484 passed, 7 skipped | this branch |
| `test/rmgpy/rmg` | 149 passed, 2 skipped, 1 failed + 1 error (`ps2pdf`, §7) | this branch |
| RMG-database `test/` | 320 passed, 4 errors — **unchanged from round 95** | this branch |

The round arrived in two parts. The first named a TOCTOU on the manifest path; the second, which
superseded it, named the `lstat` collapse, the cache identity, the field rename and three
test-quality gaps. Both are closed here. §1 is the one the second brief did not ask about.

---

## 1. The family directory was checked by name and opened by name (HIGH, from the first brief)

`_family_directory()` resolved the path with `os.path.realpath`, checked containment with
`startswith`, and returned `os.path.join(families_root, label)` — the **unresolved** string.
`_read_manifest()` then opened that string with `O_DIRECTORY` and no `O_NOFOLLOW`. Check and use
were two independent resolutions of one name, and the file at the end of it is `exec()`d.

Reproduced at `18ea8278c` by swapping from inside the check itself — `os.path.realpath` is allowed
to compute the true in-tree answer, and the directory becomes a link out of the database before
that answer is returned:

```
resolve_quarantine('A_Family_Swapped_Under_The_Check') -> (KineticsQuarantine, True)
reason = 'EXECUTED_FROM_OUTSIDE_THE_DATABASE_VIA_RACE'
the swap fired inside the containment check: True
```

**Round 95's docstring on `_read_manifest` claimed "the descriptor pins the directory that was
validated". That was false, and it was mine.** A descriptor pins whatever the name resolved to at
open time; that is not the same thing as what was checked. The prose was written from the intent
rather than from the code, which is how this campaign's overstatements keep getting made.

**The closure is structural, not a third check.** The label is already a single path component by
the time it arrives — `_family_directory` refuses separators — so the whole descent is two steps:
open the parent, open the label through it with `O_NOFOLLOW`, open the manifest through *that* with
`O_NOFOLLOW`. Neither component can be a link, neither is re-derived from a string after being
examined, and "the executed file is a direct child of a direct child of the families root" holds by
construction with no window in it.

**The manager proposed `fstat` the handle to confirm it is inside the resolved families root. That
does not work, and the difference matters.** An `fstat` returns `st_dev`/`st_ino` — identity, not
location. A file has no location in its inode; the only way to know a file is *inside* a directory
is to have reached it *through* that directory. Hence the descent rather than a post-hoc check.

**The deliberate narrowing:** a family directory that is itself a symlink is now refused even when
it points inside the database, because admitting it means deciding containment by resolving a name a
second time. Measured before taking it: **zero symlinks anywhere under `input/kinetics/families`**
in either database on this box (`RMG-database` 131 entries, the i221 worktree 142). The narrowing
stops at the label — a database root that is itself a symlink still resolves, pinned by
`test_a_database_root_that_is_itself_a_symlink_still_resolves`, because `database.directory` is
local configuration and not something a shipped `family:` line can steer.

**NUL and newline, confirmed not re-fixed.** Both are refused by name at `18ea8278c` already, as
controls in the probe rather than findings: `(None, False)` with `it contains the control character
'\x00'` / `'\n'`, nothing raised. Round 95 closed these and they stayed closed.

## 2. An unreadable manifest was reported as an absent one (HIGH)

Confirmed exactly as stated. `_manifest_signature` caught `(OSError, ValueError)` and returned
`'absent'`, so with the family directory still visible `resolve_quarantine` produced:

```
with the family directory unreadable : resolve_quarantine -> (None, True)
                      readable again : resolve_quarantine -> (KineticsQuarantine, True)
```

`(None, True)` is "this family carries no manifest" — a clean bill of health issued by the failure
of the check itself, demonstrated with a permission error rather than a deleted file. `kind` now has
a fourth value, `'unreadable'`, which takes the `answered=False` path and is **not cached**.

**One part of the brief does not reproduce, and this is the correction.** The review said the clean
answer is "cached for the rest of the run". Measured, it is not: restoring the permission moves the
signature from `(True, None, 'absent')` to `(True, <identity>, 'regular')`, which misses the cache
and forces a re-read. The entry is wrong only while the error lasts — during which the answer would
be recomputed identically anyway. What *is* true, and what the probe and test now pin instead, is
narrower: the answer is **entered into the cache at all**, under a key naming a directory and a
label, neither of which is what produced it. That is round 95's rule (`an answer that came from a
failure is not a property of the key`) and it is why the repair returns uncached regardless.

`ValueError` keeps its own branch and stays `'absent'`: an embedded NUL cannot name a file on any
filesystem, so there is nothing to be unsure about.

## 3. The cache identity was entirely restorable (MEDIUM)

`(mtime_ns, size, inode)` survives an equal-length in-place edit with `os.utime` putting the
modification time back. Reproduced at `18ea8278c`:

```
before the edit : reason = 'reason A -- xxxxxxxxxxxxxxxxxxxxxxxxxxxx'
after the edit  : reason = 'reason A -- xxxxxxxxxxxxxxxxxxxxxxxxxxxx'
(mtime_ns, size, inode) identical across the edit; st_ctime_ns moved
```

Identity is now `(mtime_ns, ctime_ns, size, inode)`. `ctime` is the right answer rather than a
content hash: no userspace call can set the inode change time, it moves on any write, and the
`lstat` is already being taken so it costs nothing. **Hashing was considered and rejected** — it
turns every lookup of every family into a file read, to cover a case `ctime` already covers anywhere
`ctime` is real. Where it is not real (a filesystem reporting it as a copy of `mtime`) this is no
weaker than it was; that is the residual and it is stated rather than papered over.

Both the probe and the test assert that the edit really did preserve the old three fields before
concluding anything — a harness that fails to construct its own case otherwise reports the defect
absent.

## 4. `requiresEngineCallSites` — the option taken, plainly (MEDIUM)

**The rename, at the field and at the docstring.** Round 95 renamed only the internal helper
(`_live_gate_calls` → `_gate_calls_in_source`) and argued the manifest field should keep its name.
That argument was made once and has now been overruled, so:

* the field is **`requiresEngineCallsInSource`**;
* `_check_engine_requirements`'s docstring states it is a **static syntactic presence check**, and
  enumerates that it does not show the call executes, dominates the admission it guards, receives
  the right arguments, or propagates what it raises;
* `requiresEngineCommit`'s refusal message now names the new field too.

**The old spelling is still read, with a warning, and this is a deliberate departure.** A renamed
field that silently stops being read disarms a pin in a repository nobody is watching — which is the
failure this whole field group exists to prevent, and `Plasma_Electron_Impact_Ionization/quarantine.py`
declares the old name in an RMG-database repository outside this round's gates. Declaring **both** is
refused rather than resolved by guessing. If the intent is that the old name must die, say so and I
will drop the alias — but the database manifest has to be updated in the same change, and that needs
the gates widened.

The rename is measured as **behaviour**, not as prose: at `18ea8278c` a manifest declaring
`requiresEngineCallsInSource = ("os",)` loads unrefused, because the field is not read at all — a
pin that is a comment. Here it is refused. (`test_the_new_spelling_is_honoured` passes on *both*
engines, for opposite reasons, which is exactly why its negative sibling exists.)

**The fifth counterexample is recorded, not repaired.** A local
`check_quarantine = lambda *a, **k: None` shadowing the import and then called satisfies the check.
It is in the probe as a NOTE and was already enumerated in `_gate_calls_in_source`'s docstring at
round 95. No sixth strengthening, per the ruling.

## 5. The three test-quality gaps (MEDIUM)

**(a) Red states that collapse into import errors.** Checked and reported rather than asserted. This
round's red state at `18ea8278c` is `9 failed, 100 passed, 7 skipped` — collection succeeded and
every failure is an assertion with a message. **Every red state this round is behavioural.** That is
possible because the new tests import no new engine symbol: they use `resolve_quarantine`,
`load_family_quarantine`, `settings` and the file's own helpers, all of which exist at the base. The
one module-level import this round adds to the engine is `_LEGACY_CALL_SITES_WARNED`, and
`_clear_gate_caches` reads it with `getattr(module, name, None)` precisely so an older engine gives a
missing cache rather than a setup error.

**(b) "Rates untouched" on one component of a four-component vector.** It checked
`lmbd_i_coefs.value_si[0]` only. It now checks the length is 4, the whole vector against
`[51487.7, -0.166019, -0.00176034, 4.42738e-07]`, and `beta`, `wr`, `wp` besides. A refit moves the
later coefficients hardest, so `[0]` alone passed for three of the four ways the entry could be
rewritten.

**(c) `all(...)` that could pass vacuously.** `test_the_affected_set_is_enumerated_from_the_database`
now asserts the enumeration is non-empty *before* the class check, so it cannot pass over an empty
list even if the recorded counts change. `TestTheAffectedSetCannotPassVacuously` drives
`affected_entries` to empty on both branches (no rules, no training depository) and asserts the
vacuous pass **explicitly**, so it is on the record rather than implicit; a third test shows the
enumeration is selective within each branch.

These three are gap-filling, not repairs: (b) and (c) are green on both engines, because the rates
really are untouched and the enumeration really is selective. Only their *ability to fail* changed,
and that is stated rather than dressed up as a red state.

## 6. Adjudications accepted without further argument

The `answered=False` fail-open policy is unchanged. The "0 tests ran" report is environmental. Both
were ruled on; neither was touched.

## 7. Named, not fixed

* **The four `test_argon_metastable_thermo.py` errors** in the RMG-database repo. Unchanged from
  round 95 at 320 passed / 4 errors; this round adds nothing to them. Out of gates.
* **`Plasma_Electron_Impact_Ionization/quarantine.py`** still declares `requiresEngineCallSites`.
  Honoured via the alias, warned about at load. Out of gates.
* **`test_make_profile_graph`** fails and errors in teardown because `ps2pdf` is absent from this
  box. Identical at the base.
* `test/rmgpy/data/rmgTest.py` and `test/rmgpy/rmg/rmgTest.py` share a basename, so the two trees
  cannot be collected in one run.

## 8. What I could not reach

* **No genuine two-process race.** §1's window is made deterministic from inside the check. The
  argument that the repair closes it is structural — no name is resolved twice — not statistical,
  and I did not build a real racer to confirm the base is exploitable under ordinary scheduling.
* **Whether any shipped database carries a hostile `family:` line, a symlinked family directory, or
  an unreadable manifest.** Measured only that there are no symlinks under `input/kinetics/families`
  in the two databases on this box.
* **`ctime` on the filesystems this actually ships on.** Verified on the local ext4 only. A
  filesystem that reports `ctime` as a copy of `mtime` leaves §3 no better than before, and nothing
  here detects that.
* **Whether the permission-error path is reachable in a real deployment.** The repair is correct
  regardless, but I have no measurement of a database directory that a running RMG cannot read.
* **The `_database_has_any_quarantine` scan** (`os.listdir` + `lexists`) still addresses manifests
  by name. It cannot cause an execution — it only decides whether to bother looking — so it is out
  of scope here, but it is the last name-addressed path in the module.
