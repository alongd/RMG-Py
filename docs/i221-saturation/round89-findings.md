# Round 89 — the key was right and three paths never delivered it

Branch `i221-saturation-no-atomtype`, base `18b94de25`, commits `e01381f06` (reproduction) and
`ed2011090` (repair). Nothing pushed; the RMG-database repo and the
`i221-argon-metastable-thermo` branch are untouched.

All seven findings were reproduced before anything was built on them, through **production
objects** — a real `KineticsLibrary`, a real `Entry`, a real `CoreEdgeReactionModel` — because the
defect this round rests on is a test that hand-built the field production omits.
`docs/i221-saturation/probes/round89_delivery_probe.py`, one file run at both ends:

| | at `18b94de25` | at `ed2011090` |
|---|---|---|
| findings reproduced | **7 of 7** | **0 of 7** |
| controls holding | 5 of 5 | 5 of 5 |

---

## HIGH 3 — the loader discarded the provenance

`KineticsLibrary.get_library_reactions()` constructed every `LibraryReaction` without
`entry=entry`, at both sites. `authoring_family()` reads the kinetics comment *and* the entry's
`long_desc`; in production the second half was unreachable. Fixed by passing the entry, with a
comment at the loop saying why it is load-bearing rather than bookkeeping.

Two things follow that are worth more than the one-line fix.

**`__reduce__` already carried it** (`library.py:122`), so provenance survives the pickle the
parallel path uses. That is now asserted rather than assumed —
`test_the_entry_survives_the_pickle_the_parallel_path_uses` — because a field that arrives and then
vanishes at a process boundary would reproduce this defect somewhere much harder to see.

**The warning was actively misleading**, exactly as reported. `_warn_unattributable()` told the
user the entry recorded no authorship and recommended adding the `family:` line that was already
there. The red probe log captures it verbatim. It is no longer emitted for that case because the
line is now read; the message itself needed no change.

## HIGH 2 — two `None`s, one value

`get_quarantine()` returned `None` for "this family carries no manifest" and for "I could not find
this family", and `add_seed_mechanism_to_core` (`model.py:1743`) makes the second case *routine*:
it converts a reaction whose family is unavailable into a library reaction rather than loading the
family, and logs that as a warning.

**Separated by `resolve_quarantine(label) -> (quarantine, answered)`.** And the more useful half of
the repair is that most of these questions turn out to be answerable: a manifest is a sidecar file
in the family's directory, so reading it needs no family object. When the family is not loaded, the
gate now reads `<database>/kinetics/families/<label>/quarantine.py` directly, cached per (database
directory, label). A routine seed conversion no longer disables the gate.

**What is left, and the argument for it.** A family this database does not contain at all is
genuinely unanswerable. That case is **admitted with its own warning, not refused**, and the choice
is deliberate:

- *Refusing* would stop ordinary runs whose seeds name foreign families — the normal case, not the
  exotic one — and would break the module's standing promise that a database with no manifest
  anywhere behaves precisely as it did before.
- *Loading the family* would mean database mutation from inside an admission gate: wrong layer,
  unbounded cost, and it can fail for reasons that have nothing to do with the question.
- *Warning* says the thing that is true: the authorship is here, the answer is not, the rate is
  being admitted, and if that family is quarantined in the database it came from, this run cannot
  know it.

The message is deliberately not the unattributable one, and a test asserts the two are
distinguishable (`"records no authoring family" not in caplog.text`).

**The bound on the warning**, which is itself a judgement: it is suppressed when *no* quarantine is
loaded anywhere, so a database carrying no manifest gains no output. A test pins that
(`test_a_database_with_no_quarantine_at_all_stays_silent`).

## HIGH 1 — the enumeration of admission paths

Now in the module docstring of `rmgpy/data/kinetics/quarantine.py`, where a reviewer can check it
against the code, and asserted by `test_the_enumeration_in_the_docstring_names_every_gate_call`,
which walks the AST of `CoreEdgeReactionModel` and fails if the gated set moves without the table
moving with it.

| path | gate |
|---|---|
| generated reaction, kinetics estimated | `apply_kinetics_to_reaction` — primary, upstream of core and edge |
| core admission | `add_reaction_to_core` |
| edge admission | `add_reaction_to_edge` |
| **pressure-dependent network** | **`add_reaction_to_unimolecular_networks` — added this round** |
| seed mechanism | covered: reaches `add_reaction_to_core` |
| reaction library | covered: reaches `add_reaction_to_edge` |
| library-to-output selection | covered upstream, deliberately not gated again — `add_reaction_library_to_output` only re-selects reactions already in `self.edge.reactions`, each of which passed the edge gate |

The pdep path was uncovered because a path reaction goes to a network *instead of* core or edge and
`generate_kinetics=False` skips estimation, so both backstops are bypassed by construction. A
network is also the last place a bad rate is still recognisable: after the k(T,P) fit its
provenance is averaged away.

## MEDIUM — the call-site check pinned an import

`requiresEngineCallSites` verified that a module *binds* the gate. Deleting every call in
`model.py` while keeping the import still passed. It now also requires a call in the module's
source, via `ast`. **Bound stated in the code**: this proves a call exists in the source, not that
it executes on a given run or covers every path. The name is kept because it now does pin what it
is named for; the residual gap is a statement, not a silence.

## MEDIUM — forgeable authorship, and the direction of its failure

Only the first `family:` line was read, so one prepended line shadowed the genuine one — a bypass
costing an attacker a single comment line. All declared labels are now consulted, so an added line
can only **widen** what is checked, never narrow it.

The remaining trust in free text is bounded and the bound is worth stating plainly: a forged or
stale line cannot invent a quarantine, because the manifest must exist and its
`appliesToKineticsClass` criterion must match. What it can do is cause a **false refusal** of an
independent rate of the quarantined class — loud, naming the file it came from, corrected by
deleting a line. The failure in the other direction is a meaningless number in a mechanism that
reports success. Biasing toward the loud error is the whole design.

---

## Addendum — the other reader of `Reaction.family`, and what I verified I did not break

**The two contracts can both be served, and are. No design question needs escalating.**

`rmgpy/electron_placement.py` keys `FAMILY_ELECTRON_PLACEMENT` on `Reaction.family` and *depends*
on the repurposing the quarantine gate had to stop trusting: `PlasmaElectronImpactIonization` is a
kinetics **library** label in that table deliberately, and the argon ionisation channel resolves
through it. The two readers want opposite things from one slot.

They coexist because `authoring_family()` never touches the slot — it reads provenance *beside* it.
That was already the design; what this round adds is proof and a guard, since my change is to the
loader that fills both.

**What I verified, and how** (`TestTheFamilySlotIsLeftAloneForItsOtherReader`, four tests, all
through the real loader and the real resolver):

1. `get_library_reactions()` on a library labelled `PlasmaElectronImpactIonization` still yields
   `reaction.family == "PlasmaElectronImpactIonization"`. Adding `entry=entry` disturbs nothing:
   `LibraryReaction.__init__` still executes `self.family = library`, and the post-construction
   `rxn.family = self.label` on the auto-generated branch is untouched.
2. **End to end through `resolve_electron_placement`**: the ionisation reaction loaded that way
   still finds its `(1, 2)` declaration and produces a balanced view — one electron on the reactant
   side, two on the product side, `electrons == 0`, canonical reaction unmutated.
3. **Both readers on one object, disagreeing correctly**: `reaction.family ==
   "PlasmaElectronImpactIonization"` (library, for placement) while `authoring_family(reaction) ==
   "Plasma_Electron_Impact_Ionization"` (family, for the gate). Neither is a normalisation of the
   other.
4. A static guard that the whole quarantine module never *assigns* to a `.family` attribute — it
   may only read the slot. That is the shape any future "cleanup" would take.

**The guard is not vacuous, and this is the part worth reading.** I applied the tempting cleanup in
a throwaway worktree — `LibraryReaction` stops writing the library label into `family` — and
measured:

```
FAILED  ...::test_placement_still_resolves_for_the_argon_ionisation_channel
        ElectronPlacementError: Reaction Ar => Arp carries no family attribution
FAILED  ...::test_the_library_label_still_lands_in_the_family_slot
FAILED  ...::test_the_two_readers_disagree_about_the_same_object_and_both_are_right
4 failed, 8 passed
```

Note **what still passed under the sabotage: five of the eight quarantine-gate tests**, including
every refusal test. The cleanup would have looked correct from inside this ticket — the gate keeps
working, because the gate deliberately stopped reading that slot — while argon ionisation placement
died. That is precisely why the hazard needed flagging from outside, and it is now the failure this
branch produces first.

Nothing was normalised, no `LibraryReaction` behaviour changed, and no separate authored-family
attribute was added. `test/rmgpy/electronPlacementTest.py` and `quarantineTest.py` run together:
**130 passed**.

---

## The audit you asked for: which other new tests build a shape production never produces

Every fixture used by the round-87 and round-89 tests, checked against the production constructor.
**One was defect-shaped, and it is the one you found.** Two more are duck-types that I am calling
acceptable, with the reason.

| fixture | production equivalent | verdict |
|---|---|---|
| `make_library_reaction(..., long_desc=...)` — assigns `reaction.entry` by hand | `get_library_reactions` — **did not** pass `entry` | **was the defect.** Kept as a unit test of the lookup, with a docstring saying so and pointing at the loader test that now anchors it |
| `make_reaction()` — `TemplateReaction(family=...)`, `.template` and `.kinetics` assigned | `__generate_reactions` sets `family` and `template`; `apply_kinetics_to_reaction` sets `kinetics` | real shape |
| `_library()` (new) — real `KineticsLibrary`, real `Entry`, `get_library_reactions()` | is the production path | real |
| manifest tests — real files, `load_family_quarantine` | is the production function | real |
| `CoreEdgeReactionModel()` + `add_reaction_to_*` | production methods, called directly | real |
| `saturatedStructureTest` — real `Molecule`, real `TransportDatabase` | production | real |
| `registered` / `_register_families` — fake `_Family` with `.quarantine` and `.label` | `KineticsFamily.__init__` sets `self.quarantine = None`; `load` sets it from the sidecar (`family.py:623,701`) | **duck-type, accepted.** The attribute and its two possible values match production exactly, and `TestTheRealQuarantinedFamily` anchors the same lookups against the real shipped family, so the fake is not the only evidence |
| `_Database` / `_Kinetics` fakes for `rmgpy.data.rmg.database` | `.kinetics.families` dict | **duck-type, accepted**, same reasoning; `TestTheDoorLeftOpen` exercises the real database |

The generalisable form, since this is the third variant of the same class on this campaign: **a
fixture that builds an object shape is a claim about production that nothing checks.** The two
duck-types above are safe only because a real-database test covers the same code path; the `entry`
fixture was not, and there was no such anchor.

---

## What I could not reach

- **`test/rmgpy/rmg/mainTest.py::TestProfiling::test_make_profile_graph` fails here**, at the fix
  and, I believe, at the base. It is environmental: **graphviz `dot` is not installed** (`ps2pdf`
  and `gs` are), so the `.ps2` intermediate is never produced and `ps2pdf` fails on a missing file.
  The code path is profiling output and is disjoint from anything this round touched. I did **not**
  get a clean before/after on this one test: the throwaway base worktree has no configured database
  directory, so its `mainTest` run gives 3 failed / 16 errors for unrelated setup reasons and is
  not a fair comparison. Stated rather than smoothed over.
- **No full RMG run.** The gate additions are tested at the unit level and through production
  constructors, not by executing the metastable deck end to end. Given the deck at
  `/home/alon/runs/ar5torr-i261-20260922-125141/input.py` is staged and this branch is the last
  blocker, that run is the obvious next check and it is not mine to start without the word.
- **The disk fallback is cached for the process lifetime.** If a run changes
  `settings['database.directory']` mid-flight, or a manifest is written while a run is in progress,
  the cached answer goes stale. I judged that out of scope — RMG does not do either — but it is a
  real bound and I am naming it rather than leaving it to be found.
- **`_counts_calls_to` is static.** A module that calls the gate inside `if False:` satisfies it.
  Closing that needs runtime instrumentation, which does not belong in a database loader.
- **Round 87's `TestProvenanceNotTheFamilySlot` still cannot be shown red at its own base**, since
  it imports a symbol that base lacks; that was stated in round 87 and has not changed.

## Measurements, each against its base

Engine `i221-saturation-no-atomtype` at `ed2011090`, plasma head `40e21b495`.

| check | result |
|---|---|
| probe at `18b94de25` | 7/7 reproduced, 5/5 controls hold |
| probe at `ed2011090` | 0/7 reproduced, 5/5 controls hold |
| new tests at `18b94de25` | **8 failed** on assertions, 47 passed, 7 skipped |
| `quarantineTest.py` at `ed2011090` | **63 passed** |
| `test/rmgpy/data` at `ed2011090` | **430 passed**, 8 skipped |
| `test/rmgpy/rmg` at `ed2011090` | 149 passed, 1 failed (graphviz missing, above), 2 skipped |

The two new tests that pass at *both* ends are the negative controls —
`test_a_family_in_the_database_without_a_manifest_is_a_clean_answer` and
`test_a_database_with_no_quarantine_at_all_stays_silent` — which must, or the repair would be
refusing things it should not.
