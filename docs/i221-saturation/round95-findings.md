# Round 95 — the reader that ignored the carrier, and the file at the end of the path

Branch `i221-saturation-no-atomtype`, base `9b89a937c`, on plasma head `98d465d3b`.
Probe: `docs/i221-saturation/probes/round95_carrier_probe.py`. Logs: `logs/round95-*`.

Every finding in the review reproduced. **10 of 10 at `9b89a937c`, 8 of 8 controls holding; 0 of 10
after, 9 of 9 controls holding.** The review was right on every count, including the LOW.

Two things the review did not ask for and that matter more than some of what it did:

1. **Closing the loaded-family staleness makes four tests in the RMG-database repository fail** —
   and the reason they fail is that the base was handing out a clean bill of health for a family
   that is quarantined on disk. §6 has the measurement. Those tests are not repaired here: the
   database repository is outside this round's gates.
2. **My first version of the repair introduced a regression that the existing suite caught**, by
   caching an answer that came from a family *object* under a key that names only a directory and a
   label. §7.

The one place I think the review's framing is wrong is `requiresEngineCallSites`, and §5 makes that
case rather than quietly doing something else.

## 1. The HIGH — the guarantee held for one shape out of two

`authoring_families()` short-circuited for anything whose `library` was `None` and returned the
single `reaction.family` slot. `library.py` writes that slot once per `family:` line, so only the
last survived; naming a quarantined family and then an ordinary *loaded* one meant the reaction was
never converted to a `LibraryReaction` either, so the shape round 89 fixed was never entered.
Reproduced through the real loader and `add_reaction_to_core()`:

```
shape built              : TemplateReaction
longDesc named           : A_Family_Quarantined_On_Disk then An_Ordinary_Loaded_Family
.family slot (last wins) : 'An_Ordinary_Loaded_Family'
authoring_families()     : ['An_Ordinary_Loaded_Family']
entry attached           : True
ADMITTED_TO_CORE         : True
```

**Repair, in the reader.** `authoring_families` now consults provenance for *every* shape — the
kinetics comment and the entry's `long_desc`, which is where round 92 put the carrier — and *adds*
the slot to that list for a template reaction rather than replacing it. The result is a superset of
what the slot alone said, in the direction that matters: an extra label can cause a false refusal,
which is loud and is corrected by deleting a line; a missing one admits a rate the database has
refused.

**The parser was deliberately left alone, and this is a disagreement worth stating.** The review
said collecting every label in `library.py` is "necessary and not sufficient". It is not necessary,
and there is nowhere legal to put the result. A second label would have to go into the `family`
slot — whose meaning is fixed, and which `model.py` compares against the loaded families and
`electron_placement.py` reads — or into a new authored-family attribute, which the round-89 addendum
forbids by name. The entry already carries every label, because it carries the text the labels were
parsed out of. Reading it *is* the complete fix, and
`test_the_family_slot_is_not_normalised` pins the slot's behaviour so that a later round cannot
"tidy" this into the forbidden repair.

**The electron-placement hazard is honoured.** `LibraryReaction.__init__` is untouched; the slot
still holds the library label after a conversion, asserted by a probe control and by two tests.

## 2. MEDIUM 3 — the conversion, and how I know the field list is complete

Both conversion sites are now one function, `rmgpy/rmg/model.py::as_library_reaction`, and the
fields it carries are **enumerated from `LibraryReaction.__init__`'s signature at import time**, not
written out:

```python
_CONVERTED_FIELDS = tuple(
    name for name in inspect.signature(LibraryReaction.__init__).parameters
    if name not in ('self', 'library')
)
```

That is the answer to "say how you determined the list is complete". A hand-written list falling
behind the class it builds *is* this defect; a test or a fix that hard-codes its own list reproduces
it one layer up. `LibraryReaction` takes **17** parameters. The old conversion passed **8**
(`reactants`, `products`, `library`, `specific_collider`, `kinetics`, `duplicate`, `reversible`,
`entry`). The **9** it dropped were `index`, `network_kinetics`, `transition_state`, `degeneracy`,
`pairs`, `allow_pdep_route`, `elementary_high_p`, `allow_max_rate_violation`, `electrons` — and a
`TemplateReaction` holds all nine, measured, so all nine are now carried.

Three of them had teeth:

- **`electrons`.** The review is right that this is a charge-balance bug in this deck and not merely
  a semantics bug: `Ar + e- => Ar+ + 2e-` arrived in the core claiming zero net electrons, and
  `electron_placement.py` reads the field.
- **`degeneracy`,** which multiplies the rate.
- **`elementary_high_p`,** which *both call sites read off the converted object* a few statements
  later to decide whether the reaction enters a pressure-dependent network. The conversion was
  clearing a flag its own caller then tested.

`library` is replaced on purpose. `index` is excluded from the probe's comparison and the reason is
shown rather than asserted: `make_new_reaction` numbers a reaction as it admits it, so the source's
`-1` legitimately becomes `1` in the core, and a separate control drives the conversion alone to
show the index coming through.

## 3. The manifest file, which is what gets executed

Round 92 constrained the family *directory*. The file inside it is what `exec()` runs.

`_read_manifest()` now opens the directory first and reads the manifest **through that descriptor**,
with `O_NOFOLLOW`:

- the descriptor pins the directory that was validated, so a rename between validation and read
  cannot redirect it — the directory-swap TOCTOU the review named;
- `O_NOFOLLOW` refuses a symlinked manifest outright instead of resolving it;
- the result must be a regular file, because a FIFO would block a run forever on the read and a
  directory would raise from inside a gate whose purpose is to keep running and refuse.

`_manifest_signature()` now uses `os.lstat`, so a manifest that is a *link* is seen as a link and
reported as `'other'` rather than as an ordinary file — the distinction that stops "refused as a
file" being read as "carries no manifest", which would be a clean bill of health issued by the check
that refused it. A refused manifest leaves the family **unanswered**.

Control characters are refused in `_family_directory`: a NUL reached `os.stat` and came back out as
an uncaught `ValueError`, killing the run from inside the gate, and a newline selected a real
directory whose manifest was executed. Shown refusing an absolute label, `..`, a symlinked family
directory, a symlinked manifest, a NUL label, a newline label and a manifest that is a directory —
and shown still resolving an ordinary one.

On a platform without `os.supports_dir_fd` the symlinked-manifest hole still closes and the
directory-swap race does not. Stated in the docstring rather than assumed away.

## 4. Cache invalidation, both halves

**Loaded families.** `resolve_quarantine` no longer returns `family.quarantine` before looking at
disk. One path, signature-validated, for loaded and unloaded alike, using the directory the family
recorded at load time (`quarantine_path`, new in `family.py`) so that a family loaded from somewhere
other than the configured database directory is re-read against *its own* manifest. A manifest
added, edited or removed after load is now seen, in both directions. What this found in the shipped
database is §6.

**The suppression cache.** `_DISK_ANY_QUARANTINE_CACHE` is keyed on the families root's signature,
which a write *inside* an existing family directory does not move. Only the positive is cached now,
and the asymmetry is the repair: a cached positive can at worst keep a warning switched on after the
last manifest is deleted, and this suppression only ever decides whether a warning prints, never
whether a rate is refused. `_warn_unanswered` de-duplicates *before* it asks, so the re-scan the
missing negative costs happens once per distinct unanswered question, not once per edge reaction per
iteration.

## 5. `requiresEngineCallSites` — the option taken, and where I think the review is wrong

**Option taken: rename and correct the documentation. The check is not strengthened a fifth time.**

The internal checker is renamed `_live_gate_calls` → `_gate_calls_in_source`, and its docstring now
enumerates what it does not verify — execution, dominance, arguments, exception propagation,
short-circuiting, local shadowing, lambdas — in the review's own words, with the reason the list is
enumerated rather than gestured at: an incomplete list is how the overstatement kept coming back.
`False and check_quarantine(r)` is accepted, the probe measures that it is accepted, and the probe
records it as a **noted fact rather than a finding**, because this round deliberately does not
change it.

**Where I disagree: the manifest field keeps its name.** The review offered "rename it", and for the
*internal* checker that is exactly right. For the public field it would be a live regression rather
than a documentation fix: `requiresEngineCallSites` is declared by
`Cation_R_Recombination/quarantine.py` in the RMG-database repository, which this round's gates
forbid me to touch, so renaming the field here would silently disarm the pin on the branch I cannot
edit. And the field name, read literally, is the one part of this that was never overstated — it
says call sites are required, and call sites in the source are exactly what is checked. The word
that was false was "reachable branch", in the docstring, and that is gone.

The ledger prose in the database repository still describes the field in terms the check did not
earn, and now under-describes it. Named in §8, not fixed.

## 6. What closing the loaded-family staleness cost, and why it is right anyway

Running the RMG-database suite against this engine gives **320 passed, 4 errors** where the same
suite against `9b89a937c` gave **324 passed**
(`logs/round95-database-suite.log`; the 324 is `docs/argon-metastable-thermo/logs/round92-db-suite.log`
in the database repository, measured at database `f5b0e83d4`).

All four errors are one module-scoped fixture in `test/test_argon_metastable_thermo.py`, which
drives the quarantined `Plasma_Electron_Impact_Ionization` estimate through the real model-admission
path in order to measure what it delivers. The gate now refuses it at `apply_kinetics_to_reaction`.

**That refusal is correct, and the base's admission was a fail-open.** Instrumenting
`resolve_quarantine` during that fixture at `9b89a937c`:

```
loaded=True  fam.quarantine=<KineticsQuarantine ... QUARANTINED ...>  manifest_on_disk=True  answer=(<KineticsQuarantine ...>, True)
loaded=True  fam.quarantine=None                                      manifest_on_disk=True  answer=(None, True)
```

The second lookup answers **"this family carries no manifest" while the manifest is on disk in that
family's own directory** — a clean bill of health for a quarantined family, from the shipped
database, in its own test run. The family object is a second instance whose `quarantine` attribute
was never populated, and the base trusted the object. This is the review's cache MEDIUM reproducing
against real data, and it is a stronger reproduction than the synthetic one in the probe.

**Not repaired here: the database repository is outside this round's gates.** The closer is small
and belongs in that repo — those four tests need to observe the estimate *before* model admission
(`family.get_kinetics(...)`) rather than through `apply_kinetics_to_reaction`, which is the boundary
the quarantine exists to hold. Until then the database branch has four failing tests against this
engine, and I am flagging that rather than leaving it to be discovered at merge.

## 7. A regression I introduced, and the test that caught it

The first version of the loaded-family repair cached the fallback answer — the one taken from
`family.quarantine` when the family's directory is not where this database would put it — in
`_DISK_QUARANTINE_CACHE`, whose key is `(directory, label)`. That key does not name the family
object, so the next lookup of the same label *with no database loaded at all* inherited a quarantine
from an object that was no longer there, and `TestTheGateFires::test_an_unloaded_database_does_not_crash_the_gate`
went red. The answer is now returned uncached, with the reason in the code.

Recording it because it is the same defect class as the round's HIGH — a value keyed on something
that does not determine it — committed by me, in the repair for it, and caught by a test written two
rounds earlier rather than by my own probe.

## 8. Named, not fixed

- **The four database tests in §6.** RMG-database repo, out of gates.
- **`Cation_R_Recombination`'s manifest prose**, which describes `requiresEngineCallSites` in terms
  the check does not earn. Same repo, same reason.
- **`test_make_profile_graph`** fails and errors in teardown at `9b89a937c` *and* here, because
  `ps2pdf` is not installed on this box — round 92's note said graphviz `dot`, which is the wrong
  tool; the traceback is `subprocess.check_call(['ps2pdf', ...])` returning 1. Not this branch's,
  and the correction is recorded so the next round does not re-diagnose it.
- **`test/rmgpy/data/rmgTest.py` and `test/rmgpy/rmg/rmgTest.py` cannot be collected in one pytest
  run** — same basename, no `__init__.py`. Pre-existing; the two trees are measured separately
  below, as round 92 did.
- **`yaml_cantera2.py:524-526`**, I-260's untouched twin. Owner's ruling stands: name it and stop.
- **I-234's false claim at `PlasmaRadiativeRecombination/reactions.py:347`.** Another ticket's file.

## 9. What I could not reach

- **Proof that a gate call executes.** Unchanged from round 92, and now stated in the code rather
  than implied against.
- **A genuine two-process race.** Both TOCTOU reproductions delete the file from inside the stat
  that precedes the read. That exercises the real code path with real filesystem state and is not a
  real race. Round 95 moved the hook from `os.path.exists` to `os.lstat`/`os.stat` because the read
  no longer calls `exists` — a race harness that no longer sits on the code path is a test that
  passes while measuring nothing, and it very nearly became one here.
- **Whether any shipped database contains a hostile `family:` line.** Still a hardening, not a
  response to a sighting; I did not survey the databases.
- **Whether the `electrons` now carried across the conversion changes any *result*.** It changes
  what the object reports. No plasma deck was re-run this round, and the round-89 addendum's warning
  about guessing in `electron_placement` is why I did not go further than carrying the field the
  source already held.

## 10. Measurements, each named against its commit

| check | result | commit |
|---|---|---|
| round-95 probe | 10 of 10 findings, 8 of 8 controls hold | `9b89a937c` |
| round-95 probe | **0 of 10**, 9 of 9 controls hold | this branch |
| `quarantineTest.py` (whole file) | **14 failed**, 80 passed, 7 skipped | `9b89a937c` |
| `quarantineTest.py` (whole file) | **101 passed** | this branch |
| `test/rmgpy/data` | **469 passed, 7 skipped** | this branch |
| `test/rmgpy/rmg` + `electronPlacementTest` | 212 passed, 2 skipped, **1 failed** (`ps2pdf`, §8) | this branch |
| RMG-database `test/` | **320 passed, 4 errors** (§6) | this branch |
| RMG-database `test/` | 324 passed | `9b89a937c` |

The 3 tests that pass at `9b89a937c` and here are the positive controls: the slot is not normalised,
the library shape already refused both labels, and an ordinary label still resolves.
