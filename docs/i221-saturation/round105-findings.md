# Round 105 — a carry that was a transformation, and an enumeration of the wrong class

Base `i221-saturation-no-atomtype@309acc0a4`. Probe:
`docs/i221-saturation/probes/round105_degeneracy_and_enumeration_probe.py`.

| check | result | commit |
|---|---|---|
| round-105 probe | **5 of 5 findings reproduce**, 6 of 6 controls hold | `309acc0a4` |
| round-105 probe | **0 of 5**, 6 of 6 controls hold, 2 noted | this branch |
| `quarantineTest.py` (new classes only) | **21 failed**, 12 passed | `309acc0a4` |
| `quarantineTest.py` | **163 passed** (was 136) | this branch |
| `test/rmgpy/data` | 530 passed, 8 skipped | this branch |
| `test/rmgpy/rmg` | 149 passed, 2 skipped, 1 failed + 1 error (`ps2pdf`) | this branch |
| RMG-database `test/` | 320 passed, 4 errors — unchanged since round 95 | this branch |

All three diagnoses are correct and all three are repaired. The HIGH is **worse than the
brief states**, and §1 gives the measurement. §4 is a census the brief did not ask for and
it found one more live instance, outside gates.

---

## 1. HIGH — confirmed, and the exposure is larger than the table shows

**Confirmed.** `_carry_entry_fields` assigned `degeneracy` through the property.
`Reaction.degeneracy`'s setter (`rmgpy/reaction.py:356`) is not an assignment: with
kinetics attached it multiplies the rate by a ratio and appends to the kinetics comment.
The attached object **is** `entry.data`, so both edits land on the shared database object.

The brief's table sweeps one shape. The ratio has two branches, and which one fires
depends on what the *constructor* was given:

```python
if self._degeneracy < 2:
    degeneracy_ratio = new                      # not (new / old)
else:
    degeneracy_ratio = (new * 1.0) / self._degeneracy
```

Only the **template** shape passes `degeneracy=` to its constructor. The two
`LibraryReaction` shapes pass none, so the old value is 1, the `< 2` branch fires, and the
ratio is the whole new degeneracy. Measured at `309acc0a4` through the real loader:

| degeneracy | template | "Originally from…" | ordinary / pdep |
|---|---|---|---|
| 0.5 | A ×0.5 | A ×0.5 | A ×0.5 |
| 1.0 | — | — | — |
| 1.5 | A ×1.5 | A ×1.5 | A ×1.5 |
| 2.0 | — | **A ×2** | **A ×2** |
| 3.0 | — | **A ×3** | **A ×3** |

and `entry.data.A` moves with it in every cell. So it is not only non-integer
degeneracies: on two of the three shapes **every** declared degeneracy other than 1
rescales the rate, and the brief's `3.0` — the value round 102's test used — trebles it.

**The comment is edited in all fifteen combinations**, including where the ratio is
exactly 1, because that branch of the setter does not depend on the ratio. This campaign
reads authorship out of a kinetics comment (`authoring_families`), so appending
`Multiplied by reaction path degeneracy 3.0` to a shared library entry is not cosmetic.

### Which fix, and why

**Write the storage, not the property** — `_CARRIED_THROUGH_STORAGE = {'degeneracy':
'_degeneracy'}`, consulted by the carry loop. The reasoning, against the two alternatives
the brief named:

* *Set it before the kinetics is attached.* This is the correct **rule**, and writing the
  storage is how the rule is obeyed here: the three shapes attach `kinetics=entry.data` in
  the constructor call itself, so "before" does not exist at this call site unless the
  loader starts building reactions with no rate and attaching one later — three call sites
  rearranged, and an object that is briefly invalid, to avoid one setter. Note that
  `Reaction.__init__` does exactly what this fix does with its own `degeneracy` argument
  (`self._degeneracy = degeneracy`). This is the constructor's own path, not a back door.
* *Deep-copy the kinetics before reassigning.* It protects the database object and leaves
  **the rate still wrong** — the copy is what gets rescaled and the copy is what the
  reaction uses. It would also end `reaction.kinetics is entry.data`, which the loader has
  always maintained, at a cost of one kinetics object per library reaction.

Carrying a value is not the same operation as changing it. A library file that declares
`degeneracy=3` is stating the rate it wants; re-deriving the rate from a number that was
already true applies it twice. The code says so, in
`rmgpy/data/kinetics/library.py:_CARRIED_THROUGH_STORAGE`.

**The workaround has a negative control.** `test_the_storage_table_names_a_setter_that_
really_transforms` asserts that assigning through the property *still* changes the rate.
If upstream ever makes the setter a plain assignment, that test fails and says to delete
the special case, so a back door with no reason left cannot survive unnoticed.

## 2. MEDIUM — the enumeration, and a fourth dropped field

**Confirmed, and it drops four, not three.** `as_library_reaction` enumerated
`LibraryReaction.__init__`'s parameters. Measured at base:

```
is_forward   True -> False
rank         7 -> None
comment      'kept?' -> ''
label        'Lip + CH3 <=> CH3Li' -> ''
```

`label` is the fourth. And the brief's diagnosis of the test is exact: round 102's test
iterated the same constructor signature, so the enumeration and its check shared a blind
spot precisely.

What the conversion preserves now comes from the **source's class** —
`REACTION_STATE_FIELDS`, the discovery `library.py` already had, which is why it is public
now rather than duplicated. `model.py` keeps its own exclusion dict with its own reasons,
because the policy differs per call site even where the enumeration does not: the loader
does not carry `index`/`label`/`comment` (a library position is not a model index), the
conversion carries all three (it is the same reaction, relabelled).

`family` is in neither list, deliberately. It is not `Reaction` state at all, and on a
`LibraryReaction` that slot means the LIBRARY — carrying the source's family label here
would undo the conversion. `test_the_conversion_still_replaces_the_library_slot` pins it.

**The acceptance "a test that would catch a newly-added `Reaction` field being dropped"**
is `test_every_carried_field_arrives`: it iterates the carried set derived from the class,
requires a marker value for each, and asserts each marker arrives. A field added to
`Reaction` tomorrow enters that test automatically and fails it until someone gives it a
marker and shows it survives — and the marker table is checked for equality against the
discovered set in both directions, so it cannot fall behind quietly.

### A third thing, found while writing that test

**`protons` was in the carried set and cannot be assigned at all.** It is a read-only
property derived from the charge balance, so `setattr` raised `AttributeError` and both
carries swallowed it. Round 102's findings listed it among the 14 carried fields; it never
was. It is now excluded with its reason, and
`test_every_field_called_carried_can_actually_be_assigned` iterates both policies and
fails on any carried field that cannot take an assignment — carried has to mean carried.

## 3. MEDIUM — `O_NOFOLLOW`

**Confirmed and refused.** With the constant deleted from `os`, the descent at base
followed a symlinked `kinetics/families` and executed a manifest outside the database,
returning `(<KineticsQuarantine …>, True)` — a clean answer produced with no containment at
all. It now refuses with its own message, beside the `dir_fd` refusal round 101 added.

**`O_DIRECTORY` is deliberately still optional, and here is the argument.** Its absence is
not a containment hole: a component that is not a directory fails the *next* `openat` with
`ENOTDIR`, and the manifest itself is checked with `stat.S_ISREG` on the descriptor that
was read. The two things `O_DIRECTORY` refuses early are both still refused, one step
later. `O_NOFOLLOW` has no such backstop, which is the difference.

## 4. The census the brief did not ask for

Every `.degeneracy = ` in `rmgpy/` — 9 sites. The property that decides all of them is
positional: **assigning `degeneracy` is safe exactly while `kinetics` is None**, measured
in the probe (10 → 10 before attaching, 10 → 40 after).

| site | kind | verdict |
|---|---|---|
| `species.py:1018` | `TransitionState.__init__` | not a `Reaction` |
| `common.py:397`, `:411` | application at generation time | safe — no kinetics exist yet; this is what the setter is *for* |
| `chemkin.pyx:895` | application | already undoes the rescale by hand with `change_rate(1/degen)`, and its own comment notes the setter appends. The same repair as this round's, done lossily |
| `reaction.py:2033` | carry in `Reaction.copy()` | safe — two lines before `other.kinetics` is attached |
| `family.py:210` | carry in `TemplateReaction.copy()` | safe — same ordering |
| `model.py:1152` | carry from the reverse reaction | safe — four lines before `reaction.kinetics = kinetics` |
| `isotopes.py:452` | carry | safe — before `rxn.kinetics = new_kinetics` |
| `isotopes.py:458` | carry | **UNSAFE.** Three lines later, after the kinetics were attached: the isotopologue's rate is multiplied by the reverse degeneracy. Out of gates, named not fixed |
| `library.py` | carry | the one this round repairs |

Two things worth keeping from this. `chemkin.pyx` shows the hazard was already known and
already worked around once, by hand, in a way that loses float precision and leaves the
comment appended. And `isotopes.py` does both orderings **six lines apart** — 452 before
the attach, 458 after — which is the clearest evidence available that the safety of the
other four is positional luck that nothing in the code marks.

## 5. The red state, and which part of it proves what

21 of the 33 new tests fail at `309acc0a4`
(`logs/round105-tests-red-at-309acc0a4.log`), and they are not all the same kind:

* **17 are behavioural** — 15 sweep cases (`assert 30.0 == 10.0 ± 1.0e-05`, `assert
  'Estimated fr...egeneracy 1.0' == 'Estimated from node Root'`), the four-fields
  conversion test, and the `O_NOFOLLOW` refusal. These measure the base's behaviour.
* **4 are structural** — they name `REACTION_STATE_FIELDS`, `_CARRIED_THROUGH_STORAGE` or
  `_NOT_CARRIED_IN_CONVERSION`, which the base does not have, so they fail with
  `AttributeError`. Round 101's ruling applies: a red state that is an import error proves
  nothing about behaviour. They pin the shape of the repair, not the defect. The behaviour
  behind two of them (`protons` unassignable; the conversion's dropped fields) is
  reproduced behaviourally by the probe and by the four-fields test respectively.

One new test is **green at base by construction**:
`test_the_conversion_does_not_rescale_the_rate`, because the constructor path never
rescaled — it is a guard on the repair, not a repair. It was shown able to fail: with
`_CARRIED_THROUGH_STORAGE` emptied it goes red (the conversion carries `degeneracy` too
now, and through the property it would halve the rate at `degeneracy=0.5`).

## 6. What I could not reach

* **Whether any shipped library entry declares a degeneracy other than 1.** The corruption
  is real at the loader; I did not sweep either database for entries that would trigger
  it, so I cannot say whether a shipped mechanism has been mis-rated by this. That sweep
  is cheap and belongs to whoever decides how far back the exposure goes.
* **`isotopes.py:458`** — out of gates. Named above, not measured through its own call
  path; the classification is from the ordering, which the probe measures as a property.
* **`model.py:1152`'s safety is positional too.** It is in gates and it is safe today, but
  only because the assignment precedes the attach *in that function*. If a reaction ever
  arrives there already carrying kinetics, it rescales. I did not change it: there is no
  red state for it, and manufacturing one would mean asserting a caller contract I cannot
  measure from here.
* **`database.py:493`'s `electrons` drop**, round 102's recommendation, is still out of
  gates and still open.
* **No real concurrent attacker** for the `O_NOFOLLOW` case; the constant is deleted from
  `os` rather than the platform lacking it. The refusal is structural, so the simulation
  is exact for the branch under test.

## 7. Named, not fixed

* `isotopes.py:458` — the live instance of this round's HIGH, outside gates (§4).
* `chemkin.pyx:895` — the hand-rolled undo, outside gates; correct in effect, lossy in
  float precision, and it still leaves the comment appended.
* `database.py:493` dropping `electrons`. Recommended as the next round's HIGH for the
  third round running.
* The four `test_argon_metastable_thermo.py` errors, unchanged at 320 passed / 4 errors
  since round 95.
* `Plasma_Electron_Impact_Ionization/quarantine.py` still declares
  `requiresEngineCallSites` (honoured as a warned alias; the repo is out of gates).
* `test_make_profile_graph` — `ps2pdf` absent, identical at base.
* The data suite reports **8** skips here against round 102's 7. Round 102 had already
  recorded this count as nondeterministic; the pass totals reconcile exactly
  (511 + 27 new = 538 = 530 + 8), so the eighth skip is that same flake, not this change.
* `report.md` §19.2's SHA list names engine commits only through `9b89a937c` and is now
  four rounds stale.
