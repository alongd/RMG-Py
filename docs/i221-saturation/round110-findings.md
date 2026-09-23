# Round 110 — the round trip carried the graph and dropped the payload

Base `i221-saturation-no-atomtype@2e4ff991d`. Gates: `rmgpy/data/kinetics/quarantine.py`,
`family.py`, `library.py`, `rmgpy/rmg/model.py` and the test tree, plus `rmgpy/reaction.py`
granted in the covering message and again in the addendum.

Both HIGHs are confirmed and closed, both MEDIUMs are closed, and the addendum is closed
with one correction. **Three findings the brief did not name were found while measuring the
two it did, and two of them are this branch's own regressions**: `copy()` raised on every
reaction `generate_reactions()` returns, and on every reaction over a surface molecule.

---

## 1. The findings, measured

`docs/i221-saturation/probes/round110_transport_payload_probe.py` at `2e4ff991d`
(`logs/round110-probe-red-at-2e4ff991d.log`) — **11 of 11 findings reproduce, 7 of 7
controls hold**. At the tip (`logs/round110-probe-after.log`) — **0 of 11, 7 of 7**.

| # | finding | at `2e4ff991d` |
|---|---|---|
| 1 | `copy()` erases per-atom state the old path preserved | ids `-32768…-32765` → `-1`, `props {'inRing': False}` → `{}` on all 10 atoms |
| 2 | a fragment reaction cannot be copied at all | `KeyError: 'R'` |
| 3 | **a reaction the family generated cannot be copied** | `ReactionStateNotCarried: 'labeled_atoms'` |
| 4 | **a reaction over a surface molecule cannot be copied** | `KeyError: 'Pt'` |
| 5 | **a copied `Species` loses state its own reducer omits** | `aug_inchi`, `creation_iteration`, `explicitly_allowed`, `symmetry_number` |
| 6 | a field in neither table is silently aliased | `copy().template is reaction.template` |
| 7 | a field nobody classified is aliased rather than refused | the copy shares the original's list |
| 8 | the manifest's identity is sampled before its bytes are read | content is the rewritten version, identity is the one before it |
| 9 | `Reaction.copy` drops `allow_max_rate_violation` | and `rank`, and `is_forward` |
| 10 | `Reaction.copy` severs `pairs` | a pairs member is a species the copy does not own |
| 11 | `deepcopy(reaction)` severs the labelled atoms | 2 labels detached |

Findings 3, 4 and 5 are mine, not the brief's. 3 and 4 are **regressions round 108
introduced**: before it, `copy()` did not touch `labeled_atoms` and used
`Species.copy(deep=True)`, which handles `metal`/`facet` correctly. Finding 3 is the more
serious of the two — `family.py:2651` deletes `labeled_atoms` once the labels have been
read back, so *every* reaction the generator returns was missing the field the deepened set
demanded, and `copy()` raised on the commonest reaction in the codebase. That it went
unnoticed for two rounds is the fixture's fault: round 108's tests build their reactions
with a constructor, and a constructed `TemplateReaction` always has `labeled_atoms` because
`__init__` sets it.

### HIGH 1, on production's own object

The fixture is the family's output, from molecules that had been through production's
`generate_resonance_structures()` — which is what assigns atom ids. Nothing under test is
assigned by the test, and an anti-vacuity test
(`test_the_fixture_carries_the_state_before_anything_is_copied`) asserts that, so the
acceptance cannot pass over `-1 == -1`.

`Atom.__reduce__` rebuilds through `Atom(element, radical_electrons, charge, label)` plus a
nine-name dict; `id`, `coords` and `props` are in neither. The old body called
`Species.copy(deep=True)` → `Graph.copy(deep=True)` → `Atom.copy()`, which sets all three
explicitly. The control `the atom state is production's, and Molecule.copy(deep=True) keeps
it` measures exactly that comparison, so the loss is the new transport's and not an absence
in the fixture.

---

## 2. No single mechanism gives both properties — measured, then composed

The brief asks this to be established rather than assumed. It is, and the answer is no:

| mechanism | references between fields (round 108) | per-object state (round 110) |
|---|---|---|
| `deepcopy` with a shared memo | **no** — `Molecule.__deepcopy__` discards the memo | yes |
| `pickle` | yes | **no** — the reducers above |
| `pickle` **with a dispatch table** | yes | yes |

The third is `pickle.Pickler.dispatch_table`, consulted *before* an object's own
`__reduce__` and keyed on the exact type. So the repair is still one transport, not a round
trip followed by a repair pass: `copy_reaction` calls `complete_round_trip`, which is a
`_CompletePickler` carrying complete reducers for the five classes whose own reducers lose
something.

Two reducers, because the classes fail in two different ways:

- **`_reduce_whole_object`** — `Molecule`, `Species`, `Fragment`, `CuttingLabel`. Their
  argument lists are not merely incomplete; for `Fragment` and `CuttingLabel` they name the
  *wrong class*, and `Molecule.__reduce__` passes `metal` and `facet` into `__init__`'s
  `inchi` and `smiles`, which are the fifth and sixth positional parameters rather than the
  seventh and eighth. The state is discovered (writable descriptors ∪ instance dict) and
  travels as the reduce tuple's third element with an explicit state setter, because the
  graph is cyclic and only the state element is applied after the object is memoised.
- **`_reduce_completed`** — `Atom` only, and for one measured reason. Its reducer stores
  `atomtype.label` and restores `ATOMTYPES[label]`, so a copied atom's type is the
  **interned** object. `AtomType` inherits identity equality and `is_specific_case_of` is a
  membership test over those objects: measured at `2e4ff991d`,
  `ATOMTYPES['Cs'].is_specific_case_of(<a pickled copy of it>)` is `False`. A whole-object
  reducer for `Atom` would buy the three lost fields at the price of silent group
  mismatches — a worse defect than the one being fixed. So `Atom` keeps its own reducer and
  has what it drops added back. A control pins the de-interning, so if upstream ever gives
  `AtomType` value equality this stops being a reason and says so.

### Completeness, and why it is testable rather than asserted

Which fields `Atom.__reduce__` drops is **derived**, not listed: a probe atom carrying a
distinguishable value in every plain field is round-tripped once, and whatever comes back
changed is what gets restored. Measured: `{coords, id, ignore, props, terminal}` — the
review's three plus two transient graph flags that `Atom.copy()` also drops and that cost
nothing to carry.

The hole in any such derivation is a field with no distinguishable value in the probe, so
`test_the_probe_that_derives_the_dropped_set_can_see_every_field` asserts the probe
distinguishes every writable plain field `Atom` has, excusing exactly three by name
(`edges`, `atomtype`, `mapping` — all references, audited by round 108's `is` assertions)
and failing if that excuse goes stale. An `Atom` that grows a field turns it red and names
it.

Membership of the registry is derived the same way:
`test_every_class_that_loses_state_is_registered` round-trips every class reachable in a
reaction's deepened state through its *own* reducer and requires that anything losing state
is registered **and** that anything registered still loses something. A workaround cannot
outlive its reason in either direction.

---

## 3. The census — a count and an enumeration

Twelve classes cross a serialisation, copy or reconstruction boundary in a reaction's
deepened state on this branch. **Five lose state through their own reducer; seven lose
nothing.** Measured, not read off the source (`the census` note in both probe logs):

| class | what does not survive its own reducer |
|---|---|
| `Atom` (as production holds it) | `id`, `props` |
| `Atom` (every plain field set) | `coords`, `id`, `ignore`, `props`, `terminal` |
| `Molecule` (`metal='Pt'`) | **raises `KeyError: 'Pt'`** — `metal`/`facet` misrouted into `inchi`/`smiles` |
| `Species` | `aug_inchi`, `creation_iteration`, `explicitly_allowed`, `symmetry_number` |
| `Fragment` | **raises `KeyError: 'R'`**, and would rebuild as `Molecule` |
| `CuttingLabel` | **raises `KeyError: 'R'`**, and would rebuild as `Atom` |
| `Bond` | nothing |
| `Molecule` (ordinary) | nothing |
| `TransitionState` | nothing |
| `Arrhenius` | nothing |
| `TemplateReaction` | nothing (its losses are all inside the objects it holds) |
| `LibraryReaction` | nothing |

Plain data only. An object whose `__eq__` is identity compares unequal to its own faithful
copy — every `Atom`, `Species` and `Bond` here does — so reference fields are audited by
round 108's `is` assertions instead, and mixing the two would make the census unreadable in
both directions.

**Zero reported as zero:** seven of the twelve lose nothing, and none of those seven is
registered.

Boundaries outside the deepened state, for completeness:

| site | what does not survive |
|---|---|
| `Reaction.copy` (`reaction.py:2017`) | **was** `allow_max_rate_violation`, `rank`, `is_forward`, `k_effective_cache`, and `pairs` identity — all five closed here |
| `Reaction.__reduce__` (`reaction.py:322`) | **was** `allow_max_rate_violation`, `is_forward` — closed here |
| `deepcopy(reaction)` | **was** labelled-atom identity — closed for the two classes whose `copy()` preserves the class |
| plain `pickle.dumps(reaction)` | still the five lossy reducers — see §7 |
| `DepositoryReaction.__reduce__`, `PDepReaction.__reduce__` | hand-enumerated, out of gates, named for the third round |
| `KineticsDatabase.__reduce__` region (`database.py:493`) | `electrons`, named for the sixth round |

---

## 4. The collapsed-assertion audit

The brief asks for one line per assertion in round 108's test set that compares a collapsed
form. Round 108's set is `TestTheCopyIsOneGraphNotFourLists` (14 tests); round 107's
comparator table is included because it is the same file and the same blindness.

| where | the collapsed form | what it cannot catch |
|---|---|---|
| `test_relabelling_...is_observable_in_the_structure` | `a.label == "*9"` — one string field of the atom | anything else about that atom: **this is the assertion that sat on an atom whose `id` and `props` had just been erased and passed** |
| `test_molecule_copy_still_preserves_atom_order` | `[a.element.symbol for a in other.atoms]` — symbols instead of atoms | an atom with the right element and a reset `id`, `coords` or `props`. The most on-point one: it compares exactly the form that hid HIGH 1 |
| `test_the_two_copies_are_one_implementation` | `"copy_reaction(" in body`, `"deepcopy" not in body` — source text | whether the shared helper's transport reproduces anything; it stayed green through all five losses |
| `test_every_pair_member_is_one_of_the_copys_own_species` | `is` over species — not collapsed, but shallow | *which* species a pair member is, never what that species holds: blind to `Species` losing `symmetry_number` |
| `test_every_labelled_atom_is_inside_...` | `is` over atoms — same shape | an atom that is the right object and has lost its `id` |
| `_COMPARED_BY["reactants"]`/`["products"]` (round 107) | `[s.label for s in v]` | everything about a species except its label — the whole of finding 5 |
| `_COMPARED_BY["kinetics"]`/`["network_kinetics"]` | `(v.A.value_si, v.n.value_si)` | `Ea`, `T0`, `Tmin/Tmax`, `comment`, and the class of the rate model |
| `_COMPARED_BY["transition_state"]` | `type(v).__name__` | every field of the transition state |
| `_COMPARED_BY["entry"]` | `(v.index, v.label)` | `item`, `data`, `long_desc` — including the `family:` line the quarantine gate reads |
| `_COMPARED_BY["labeled_atoms"]` | `sorted(v)` — keys instead of values | which atoms are labelled at all |
| `_COMPARED_BY["pairs"]` | `[tuple(s.label for s in pair)]` | identity, by design — round 108 added the `is` tests for exactly this |
| `_COMPARED_BY["reverse"]` | `v.index` | the entire reverse reaction |

Two of this round's findings were invisible to all of it. The pattern is one thing: **every
collapsed comparison is a projection, and a projection is blind to everything it projects
away.** Round 108 added identity assertions because value comparison could not see
references; round 110 needed state assertions because identity comparison cannot see
payload. The general fix is not a third kind of assertion — it is that a test which
compares a *projection* must say which projection and why, the way `_COMPARED_BY` now does
in its comment and `test_the_copy_keeps_both_the_references_and_the_atom_state` does by
asserting both properties on one copy.

---

## 5. MEDIUM 1 — the partition has no default left

`template` (a list the family rewrites in place at `family.py:2643`) and
`specific_collider` (a `Species`) were in neither table, so they were aliased. Both are now
in the deepened tables, and — the structural half — **there is no default**:
`copy_reaction` computes `state_fields(reaction) - deepened - _COPIED_BY_REFERENCE` and
raises `ReactionStateUnclassified` naming any field left over.

`_COPIED_BY_REFERENCE` carries seventeen entries with a reason each, including the two the
addendum asked for by name: `allow_max_rate_violation` ("dropped by `Reaction.copy` and by
`Reaction.__reduce__` and has been lost in production by both") and `is_forward` ("the one
whose default is not its common value"). `entry` and `reverse` say why they are shared
rather than deepened.

The brief is right that the old future-field test demanded only a marker, which an alias
satisfies: `test_every_reproduced_field_arrives` compares through `_COMPARED_BY`, and a
shared list compares equal to itself. It is left as it is — it audits `__reduce__`, where
aliasing is not possible — and the partition test
(`test_the_two_tables_are_a_partition_of_what_a_reaction_holds`) plus the refusal are what
close the copy side.

---

## 6. MEDIUM 2 — the identity brackets the read

`os.fstat(fd)` now happens before **and** after `handle.read()`, on the same descriptor,
and the two must agree. `st_atime` is deliberately not part of the identity: reading the
file changes it, and an identity containing it could never be compared with itself across a
read.

The test rewrites the manifest **in place** between the two calls — same inode, new bytes —
by wrapping `os.fstat` so the first call fires the rewrite after stat-ing, and requires
`(None, None)`. It asserts the rewrite actually fired, so it cannot pass by never entering
the window. A companion test pins that an unchanged manifest is still read, and a third
pins the *shape*: two `os.fstat(` calls and no `os.stat(` anywhere in `_read_manifest`, so
a later repair cannot satisfy the first test by re-stat'ing the name and reintroducing the
two-resolutions defect round 95 closed.

---

## 7. The addendum — one site was already closed, one was live

**`get_library_reactions` already carries `allow_max_rate_violation` on this branch, and
has since round 102.** Measured at both tips: an entry whose `item` carries the flag
produces a `LibraryReaction` carrying it. The review is right that the constructor calls
omit it — all three of them — and right for any branch without round 102's derived carry;
here `_carry_entry_fields` supplies it on the next line from a field set discovered from
`Reaction` rather than written down. That is what the derived carry is *for*, and it is the
cleanest evidence so far that it works: a field nobody remembered arrived anyway. Recorded
as a control (`get_library_reactions already carries allow_max_rate_violation here`) and as
a green-at-base test, not claimed as a repair.

**`Reaction.copy` was live**, and dropped `rank` and `is_forward` beside it — the same
hand-enumeration missing the same way. All three are carried now, `k_effective_cache` is
initialised (`__new__` leaves a `cdef public dict` unset and reading one raises, so the copy
could not be asked for a rate coefficient), `pairs` is remapped onto the copy's own species
instead of deep-copied separately, and `Reaction.__reduce__` carries the two fields it
stopped short of — they are the next two parameters of `__init__`, so the positional tuple
needed only extending.

`rmgpy/reaction.py` is a **compiled** module: the `.so` was backed up and rebuilt with
`python setup.py build_ext --inplace`. The first attempt failed with *"closures inside cpdef
functions not yet supported"* — the remap was written as comprehensions — and is written as
loops with a comment saying why.

### Why `__deepcopy__` is on the two subclasses and not on `Reaction`

The covering message put the deep-copy entry point in scope if the census showed it severed
the same references. It does (finding 11). But the obvious one-liner on `Reaction` would
have been a **regression**, and the measurement says so in two ways: `deepcopy` today
preserves `pairs` identity (one memo, and `Species` has no `__deepcopy__` override) where
`Reaction.copy` severed it, and `deepcopy` today preserves the subclass where
`Reaction.copy` returns a base `Reaction` — so `DepositoryReaction` and `PDepReaction`
would have lost their class. So: `pairs` is fixed in `Reaction.copy` first, and
`__deepcopy__` is defined on `TemplateReaction` and `LibraryReaction` only, whose `copy()`
preserves the class. Each memoises itself, so two references to one reaction inside one
`deepcopy` still come out as one copy; they are not memo-consistent with objects copied
outside, which is `copy()`'s own limitation and is stated in the docstring.

---

## 8. Red state

`logs/round110-tests-red-at-2e4ff991d.log` — 23 new tests, **18 fail at the base and 5
pass**:

- **14 behavioural** — the atom state and identity on one copy; the family's own output;
  `deepcopy`; the fragment reaction; the surface reaction; the `Species` fields; `template`
  and `specific_collider` aliasing; the unclassified-field refusal (which names only
  `ReactionStateNotCarried`, the base class the base already had, so its red state is the
  missing refusal rather than a missing import); the manifest rewritten mid-read; and the
  four `Reaction.copy`/`__reduce__` ones.
- **4 structural** — the registry census, the probe-completeness test, the partition test,
  and the `_COPIED_BY_REFERENCE` naming test; each fails on a name the repair adds, and
  each says so in its own docstring.
- **5 green at the base by construction** — the anti-vacuity fixture check, the interned
  atom types, `deepcopy` of a fragment working, an unchanged manifest still being read, and
  the loader already carrying the flag.

## 9. Counts, each against its commit

| suite | at `2e4ff991d` | at the tip |
|---|---|---|
| `round110_transport_payload_probe.py` | 11/11 findings, 7/7 controls | **0/11**, 7/7 |
| `test/rmgpy/data/kinetics/quarantineTest.py` | 18 failed / 200 passed | **218 passed** |
| `test/rmgpy/data/` | 563 passed, 7 skipped (round 108) | **586 passed, 7 skipped** |
| `test/rmgpy/rmg/` + `reactionTest.py` + `reactionChargeTransferTest.py` + `reactionMarcusTest.py` | 305 passed, 4 skipped, **1 failed + 1 error** | 305 passed, 4 skipped, **1 failed + 1 error** |
| RMG-database `test/` | 320 passed, 4 errors (round 108) | **320 passed, 4 errors** |

`586 = 563 + 23`, the 23 being this round's new tests; nothing else moved.

The rmg/reaction line is byte-for-byte the same at both tips, which is the point of running
it: `rmgpy/reaction.py` is a compiled module and this round changed it, so the suites that
exercise `Reaction` directly were measured at the base for comparison rather than assumed
(`logs/round110-rmg-suite-red-at-2e4ff991d.log` against `logs/round110-rmg-suite-after.log`).
The one failure and one error are both `mainTest.py::TestProfiling::test_make_profile_graph`,
which needs `ps2pdf` and does not have it; identical at the base, as it has been for eight
rounds. The database suite's four errors are round 95's intended ones.

## 10. What I could not reach

- **Whether a real run ever consumed a copy with erased atom ids.** The mechanism is
  proven and the consumers are named (`isotopes.py:394` copies whatever it is handed), but
  I did not run an isotope job and watch a resonance correspondence fail.
- **Whether `copy()` was ever called on a fragment reaction in anger.** Fragment
  chemistry is a separate workflow that this branch does not exercise; the finding is that
  a working feature was broken, not that someone hit it.
- **The residue: a plain `pickle.dumps(reaction)` still loses per-atom state.**
  `__reduce__` cannot dictate how the objects nested inside its state are pickled, so only
  a pickler carrying the dispatch table — which is what `copy()` uses — gets the complete
  reducers. Closing it needs either a change to `rmgpy/molecule/` (forbidden this round) or
  a process-wide `copyreg.pickle()` registration, which would change the behaviour of code
  that never imported this module. **This is my recommendation for the next round, and it
  belongs to the owner**: the correct fix really is in `Atom.__reduce__` and
  `Molecule.__reduce__`, and the brief asked to be told if so.
- **Whether `Molecule.__reduce__`'s misrouted `metal`/`facet` has corrupted a surface
  model.** Any code path that pickles a surface molecule hits it — it is upstream and
  predates this branch entirely — but I did not survey surface workflows.
- **Cost in a real generation loop.** Measured on the family's own reaction, 300
  iterations: `copy()` is **0.104 ms against 0.069 ms** at the base — ~50% slower, which
  is what the completeness costs — and `deepcopy(reaction)` is **0.142 ms against
  0.282 ms**, twice as fast, because it now routes through `copy()` instead of recursing
  field by field. A microbenchmark on one reaction, not a profile of a run.
- **`Species.is_solvent`, `henry_law_constant_data` and
  `liquid_volumetric_mass_transfer_coefficient_data`** are also absent from
  `Species.__reduce__`; they are carried now by construction, but I did not construct a
  case where losing them mattered.

## 11. Named, not fixed

- `rmgpy/molecule/molecule.py:140` — `Atom.__reduce__` omits `id`, `coords`, `props`,
  `terminal`, `ignore`. The root of finding 1.
- `rmgpy/molecule/molecule.py:1120` — `Molecule.__reduce__` passes `metal`/`facet` into
  `inchi`/`smiles`. The root of finding 4, and a live corruption of any pickled surface
  molecule anywhere in RMG.
- `rmgpy/molecule/molecule.py:1064` — `Molecule.__deepcopy__` discards the memo.
- `rmgpy/molecule/fragment.py` — `Fragment` and `CuttingLabel` inherit reducers that name
  the wrong class. The root of finding 2.
- `rmgpy/species.py:192` — `Species.__reduce__` omits four fields (seven counting the
  three in §10).
- `rmgpy/data/kinetics/depository.py:83`, `rmgpy/rmg/pdep.py:90` — hand-enumerated
  `__reduce__`, third round named.
- `rmgpy/tools/isotopes.py:458` — round 105's HIGH, still live.
- `rmgpy/data/kinetics/database.py:493` — drops `electrons`; sixth round named.
- `rmgpy/data/kinetics/family.py:3854`, `rmgpy/data/kinetics/database.py:755` —
  `deepcopy(reaction)`: closed for `TemplateReaction` and `LibraryReaction` by the new
  `__deepcopy__`, still severed for any other subclass.
- `report.md` §19.2 (database repo) — seven engine commits stale. Update before any merge,
  not before.

**Fixed this round and previously named:** `rmgpy/rmg/model.py:104`'s two stale exclusion
reasons ("a cimported type" → the measured wording), which were out of round 107's gates
and are in this round's.
